/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined(__LIBTORCH)

#include <ATen/Parallel.h>
#include <c10/core/DeviceGuard.h>
#include <torch/csrc/api/include/torch/cuda.h>
#include <torch/script.h>

#include "offload/offload_library.h"

#include <cassert>

#include <cfenv>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

typedef torch::Tensor torch_c_tensor_t;
typedef c10::Dict<std::string, torch::Tensor> torch_c_dict_t;
typedef torch::jit::Module torch_c_model_t;

class TorchFloatingPointMaskGuard {
public:
  TorchFloatingPointMaskGuard() : active_(std::feholdexcept(&env_) == 0) {}
  ~TorchFloatingPointMaskGuard() {
    if (active_) {
      std::feclearexcept(FE_ALL_EXCEPT);
      std::fesetenv(&env_);
    }
  }

private:
  std::fenv_t env_;
  bool active_;
};

/*******************************************************************************
 * \brief Internal helper for selecting the CUDA device when available.
 * \author Ole Schuett
 ******************************************************************************/
static bool use_cuda_if_available = true;

static bool get_positive_int_env(const char *name, int &value) {
  const char *raw = std::getenv(name);
  if (raw == nullptr || raw[0] == '\0') {
    return false;
  }
  char *end = nullptr;
  const long parsed = std::strtol(raw, &end, 10);
  if (end == raw || *end != '\0' || parsed <= 0 || parsed > INT_MAX) {
    return false;
  }
  value = static_cast<int>(parsed);
  return true;
}

static void initialize_torch_threads_from_env() {
  static bool initialized = false;
  if (initialized) {
    return;
  }
  initialized = true;

  int num_threads = 0;
  if (get_positive_int_env("CP2K_TORCH_NUM_THREADS", num_threads)) {
    at::set_num_threads(num_threads);
  }
  if (get_positive_int_env("CP2K_TORCH_NUM_INTEROP_THREADS", num_threads)) {
    at::set_num_interop_threads(num_threads);
  }
}

static torch::Device get_device() {
  initialize_torch_threads_from_env();
  if (!use_cuda_if_available || !torch::cuda::is_available()) {
    return torch::kCPU;
  }
  const auto device_count = torch::cuda::device_count();
  if (device_count <= 0) {
    return torch::kCPU;
  }
  const int chosen_device = offload_get_chosen_device();
  const int device = (chosen_device >= 0) ? chosen_device : 0;
  assert(device < device_count);
  return torch::Device(torch::kCUDA, device);
}

static torch::Device get_device_with_guard(c10::OptionalDeviceGuard &guard) {
  const auto device = get_device();
  if (device.is_cuda()) {
    guard.reset_device(device);
  }
  return device;
}

static void set_jit_fusion_strategy() {
  // JIT Fusion strategy optimization, hardcode dynamic 10, see also
  // https://github.com/mir-group/pair_nequip_allegro.git
  torch::jit::FusionStrategy strategy = {
      {torch::jit::FusionBehavior::DYNAMIC, 10}};
  torch::jit::setFusionStrategy(strategy);
}

static void copy_string_to_c_buffer(const std::string &source, char **content,
                                    int *length) {
  *length = source.length();
  *content = (char *)malloc(source.length() + 1); // +1 for null terminator
  strcpy(*content, source.c_str());
}

static bool can_load_directly_to_device(const torch::Device &device) {
  return !device.is_cuda() || device.index() == 0 ||
         torch::cuda::device_count() == 1;
}

static torch::jit::Module load_module_for_device(const char *filename,
                                                 const torch::Device &device) {
  if (can_load_directly_to_device(device)) {
    return torch::jit::load(filename, device);
  }
  auto model = torch::jit::load(filename, torch::kCPU);
  model.to(device);
  return model;
}

/*******************************************************************************
 * \brief Internal helper for creating a Torch tensor from an array.
 * \author Ole Schuett
 ******************************************************************************/
static torch_c_tensor_t *tensor_from_array(const torch::Dtype dtype,
                                           const bool req_grad, const int ndims,
                                           const int64_t sizes[],
                                           void *source) {
  initialize_torch_threads_from_env();
  const auto opts = torch::TensorOptions().dtype(dtype).requires_grad(req_grad);
  const auto sizes_ref = c10::IntArrayRef(sizes, ndims);
  return new torch_c_tensor_t(torch::from_blob(source, sizes_ref, opts));
}

static bool tensor_matches(const torch_c_tensor_t *tensor,
                           const torch::Dtype dtype,
                           const torch::Device &device, const int ndims,
                           const int64_t sizes[]) {
  if (tensor == nullptr || !tensor->defined() ||
      tensor->scalar_type() != dtype || tensor->device() != device ||
      tensor->ndimension() != ndims) {
    return false;
  }
  for (int i = 0; i < ndims; i++) {
    if (tensor->size(i) != sizes[i]) {
      return false;
    }
  }
  return tensor->is_contiguous();
}

static void reset_tensor_from_array_double(torch_c_tensor_t **tensor,
                                           const bool req_grad, const int ndims,
                                           const int64_t sizes[],
                                           double source[]) {
  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  const auto sizes_ref = c10::IntArrayRef(sizes, ndims);
  if (!tensor_matches(*tensor, torch::kFloat64, device, ndims, sizes)) {
    delete (*tensor);
    const auto opts =
        torch::TensorOptions().dtype(torch::kFloat64).device(device);
    *tensor = new torch_c_tensor_t(torch::empty(sizes_ref, opts).detach());
  }
  const auto source_tensor = torch::from_blob(
      source, sizes_ref, torch::TensorOptions().dtype(torch::kFloat64));
  {
    torch::NoGradGuard no_grad;
    (*tensor)->copy_(source_tensor);
    (*tensor)->mutable_grad() = torch::Tensor();
  }
  (*tensor)->set_requires_grad(req_grad);
}

/*******************************************************************************
 * \brief Internal helper for getting the data_ptr and sizes of a Torch tensor.
 * \author Ole Schuett
 ******************************************************************************/
static void *get_data_ptr(const torch_c_tensor_t *tensor,
                          const torch::Dtype dtype, const int ndims,
                          int64_t sizes[]) {
  assert(tensor->scalar_type() == dtype);
  assert(tensor->ndimension() == ndims);
  for (int i = 0; i < ndims; i++) {
    sizes[i] = tensor->size(i);
  }

  assert(tensor->is_contiguous());
  return tensor->data_ptr();
};

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 * \brief Creates a Torch tensor from an array of int32s.
 *        The passed array has to outlive the tensor!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_from_array_int32(torch_c_tensor_t **tensor,
                                     const bool req_grad, const int ndims,
                                     const int64_t sizes[], int32_t source[]) {
  *tensor = tensor_from_array(torch::kInt32, req_grad, ndims, sizes, source);
}

/*******************************************************************************
 * \brief Creates a Torch tensor from an array of floats.
 *        The passed array has to outlive the tensor!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_from_array_float(torch_c_tensor_t **tensor,
                                     const bool req_grad, const int ndims,
                                     const int64_t sizes[], float source[]) {
  *tensor = tensor_from_array(torch::kFloat32, req_grad, ndims, sizes, source);
}

/*******************************************************************************
 * \brief Creates a Torch tensor from an array of int64s.
 *        The passed array has to outlive the tensor!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_from_array_int64(torch_c_tensor_t **tensor,
                                     const bool req_grad, const int ndims,
                                     const int64_t sizes[], int64_t source[]) {
  *tensor = tensor_from_array(torch::kInt64, req_grad, ndims, sizes, source);
}

/*******************************************************************************
 * \brief Creates a Torch tensor from an array of doubles.
 *        The passed array has to outlive the tensor!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_from_array_double(torch_c_tensor_t **tensor,
                                      const bool req_grad, const int ndims,
                                      const int64_t sizes[], double source[]) {
  *tensor = tensor_from_array(torch::kFloat64, req_grad, ndims, sizes, source);
}

/*******************************************************************************
 * \brief Releases a string returned from the Torch C API.
 ******************************************************************************/
void torch_c_free_string(char *content) { free(content); }

/*******************************************************************************
 * \brief Reuses or creates a device tensor and copies double data into it.
 ******************************************************************************/
void torch_c_tensor_reset_from_array_double(torch_c_tensor_t **tensor,
                                            const bool req_grad,
                                            const int ndims,
                                            const int64_t sizes[],
                                            double source[]) {
  reset_tensor_from_array_double(tensor, req_grad, ndims, sizes, source);
}

/*******************************************************************************
 * \brief Creates an expanded tensor view along one singleton dimension.
 ******************************************************************************/
void torch_c_tensor_expand_dim(const torch_c_tensor_t *tensor,
                               const int64_t dim, const int64_t size,
                               torch_c_tensor_t **result) {
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  assert(*result == NULL);
  assert(dim >= 0);
  assert(dim < tensor->dim());
  std::vector<int64_t> sizes(tensor->sizes().begin(), tensor->sizes().end());
  assert(sizes[dim] == 1);
  sizes[dim] = size;
  *result = new torch_c_tensor_t(tensor->expand(sizes));
}

/*******************************************************************************
 * \brief Creates a tensor view narrowed along one dimension.
 ******************************************************************************/
void torch_c_tensor_narrow(const torch_c_tensor_t *tensor, const int64_t dim,
                           const int64_t start_index, const int64_t length,
                           torch_c_tensor_t **result) {
  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  assert(*result == NULL);
  assert(dim >= 0);
  assert(start_index >= 0);
  assert(length >= 0);
  assert(dim < tensor->ndimension());
  assert(start_index + length <= tensor->size(dim));
  *result =
      new torch_c_tensor_t(tensor->narrow(dim, start_index, length).to(device));
}

/*******************************************************************************
 * \brief Returns the data_ptr and sizes of a Torch tensor of int32s.
 *        The returned pointer is only valide during the tensor's live time!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_data_ptr_int32(const torch_c_tensor_t *tensor,
                                   const int ndims, int64_t sizes[],
                                   int32_t **data_ptr) {
  *data_ptr = (int32_t *)get_data_ptr(tensor, torch::kInt32, ndims, sizes);
}

/*******************************************************************************
 * \brief Returns the data_ptr and sizes of a Torch tensor of floats.
 *        The returned pointer is only valide during the tensor's lifetime!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_data_ptr_float(const torch_c_tensor_t *tensor,
                                   const int ndims, int64_t sizes[],
                                   float **data_ptr) {
  *data_ptr = (float *)get_data_ptr(tensor, torch::kFloat32, ndims, sizes);
}

/*******************************************************************************
 * \brief Returns the data_ptr and sizes of a Torch tensor of int64s.
 *        The returned pointer is only valide during the tensor's live time!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_data_ptr_int64(const torch_c_tensor_t *tensor,
                                   const int ndims, int64_t sizes[],
                                   int64_t **data_ptr) {
  *data_ptr = (int64_t *)get_data_ptr(tensor, torch::kInt64, ndims, sizes);
}

/*******************************************************************************
 * \brief Returns the data_ptr and sizes of a Torch tensor of doubles.
 *        The returned pointer is only valide during the tensor's live time!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_data_ptr_double(const torch_c_tensor_t *tensor,
                                    const int ndims, int64_t sizes[],
                                    double **data_ptr) {
  *data_ptr = (double *)get_data_ptr(tensor, torch::kFloat64, ndims, sizes);
}

/*******************************************************************************
 * \brief Runs autograd on a Torch tensor.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_backward(const torch_c_tensor_t *tensor,
                             const torch_c_tensor_t *outer_grad) {
  TorchFloatingPointMaskGuard fpe_guard;
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  tensor->backward(*outer_grad);
}

/*******************************************************************************
 * \brief Runs autograd on a scalar Torch tensor.
 ******************************************************************************/
void torch_c_tensor_backward_scalar(const torch_c_tensor_t *tensor) {
  TorchFloatingPointMaskGuard fpe_guard;
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  tensor->backward();
}

/*******************************************************************************
 * \brief Moves a tensor to the active device and makes it an autograd leaf.
 ******************************************************************************/
void torch_c_tensor_to_device_leaf(torch_c_tensor_t **tensor,
                                   const bool req_grad) {
  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  auto moved = (*tensor)->to(device).detach();
  moved.set_requires_grad(req_grad);
  delete (*tensor);
  *tensor = new torch_c_tensor_t(moved);
}

/*******************************************************************************
 * \brief Select whether Torch wrappers should use CUDA when available.
 ******************************************************************************/
void torch_c_use_cuda(const bool use_cuda) { use_cuda_if_available = use_cuda; }

/*******************************************************************************
 * \brief Returns the gradient of a Torch tensor which was computed by autograd.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_grad(const torch_c_tensor_t *tensor,
                         torch_c_tensor_t **grad) {
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  const torch::Tensor maybe_grad = tensor->grad();
  assert(maybe_grad.defined());
  *grad = new torch_c_tensor_t(maybe_grad.cpu().contiguous());
}

/*******************************************************************************
 * \brief Releases a Torch tensor and all its ressources.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_release(torch_c_tensor_t *tensor) { delete (tensor); }

/*******************************************************************************
 * \brief Creates an empty Torch dictionary.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_dict_create(torch_c_dict_t **dict_out) {
  assert(*dict_out == NULL);
  *dict_out = new c10::Dict<std::string, torch::Tensor>();
}

/*******************************************************************************
 * \brief Clones a Torch dictionary.
 ******************************************************************************/
void torch_c_dict_clone(const torch_c_dict_t *dict, torch_c_dict_t **dict_out) {
  assert(*dict_out == NULL);
  torch_c_dict_t *clone = new c10::Dict<std::string, torch::Tensor>();
  for (const auto &entry : *dict) {
    clone->insert(entry.key(), entry.value());
  }
  *dict_out = clone;
}

/*******************************************************************************
 * \brief Inserts a Torch tensor into a Torch dictionary.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_dict_insert(const torch_c_dict_t *dict, const char *key,
                         const torch_c_tensor_t *tensor) {
  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  dict->insert(key, tensor->to(device));
}

/*******************************************************************************
 * \brief Retrieves a Torch tensor from a Torch dictionary.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_dict_get(const torch_c_dict_t *dict, const char *key,
                      torch_c_tensor_t **tensor) {
  assert(dict->contains(key));
  *tensor = new torch_c_tensor_t(dict->at(key).cpu().contiguous());
}

/*******************************************************************************
 * \brief Releases a Torch dictionary and all its ressources.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_dict_release(torch_c_dict_t *dict) { delete (dict); }

/*******************************************************************************
 * \brief Loads a Torch model from given "*.pth" file.
 *        In Torch lingo models are called modules.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_load(torch_c_model_t **model_out, const char *filename) {
  assert(*model_out == NULL);
  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  set_jit_fusion_strategy();
  torch::jit::Module *model = new torch::jit::Module();
  *model = load_module_for_device(filename, device);
  model->eval(); // Set to evaluation mode to disable gradients, drop-out, etc.
  *model_out = model;
}

/*******************************************************************************
 * \brief Evaluates the given Torch model.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_forward(torch_c_model_t *model, const torch_c_dict_t *inputs,
                           torch_c_dict_t *outputs) {

  TorchFloatingPointMaskGuard fpe_guard;
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  auto untyped_output = model->forward({*inputs}).toGenericDict();
  outputs->clear();
  for (const auto &entry : untyped_output) {
    outputs->insert(entry.key().toStringView(), entry.value().toTensor());
  }
}

/*******************************************************************************
 * \brief Evaluates a TorchScript model method expecting argument "mol".
 ******************************************************************************/
void torch_c_model_forward_mol_tensor(torch_c_model_t *model,
                                      const char *method_name,
                                      const torch_c_dict_t *inputs,
                                      torch_c_tensor_t **output) {

  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  assert(*output == NULL);
  *output = new torch_c_tensor_t(
      model->get_method(method_name)({*inputs}).toTensor());
}

/*******************************************************************************
 * \brief Returns the weighted sum of two Torch tensors.
 ******************************************************************************/
void torch_c_tensor_weighted_sum(const torch_c_tensor_t *values,
                                 const torch_c_tensor_t *weights,
                                 torch_c_tensor_t **result) {
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  const auto weights_on_device = weights->to(values->device());
  *result = new torch_c_tensor_t((*values * weights_on_device).sum());
}

/*******************************************************************************
 * \brief Returns a scalar double value from a Torch tensor.
 ******************************************************************************/
double torch_c_tensor_item_double(const torch_c_tensor_t *tensor) {
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  return tensor->item<double>();
}

/*******************************************************************************
 * \brief Releases a Torch model and all its ressources.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_release(torch_c_model_t *model) { delete (model); }

/*******************************************************************************
 * \brief Reads metadata entry from given "*.pth" file.
 *        In Torch lingo they are called extra files.
 *        The returned char array has to be deallocated by caller!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_read_metadata(const char *filename, const char *key,
                                 char **content, int *length) {

  std::unordered_map<std::string, std::string> extra_files = {{key, ""}};
  torch::jit::load(filename, torch::kCPU, extra_files);
  const std::string &content_str = extra_files[key];
  copy_string_to_c_buffer(content_str, content, length);
}

/*******************************************************************************
 * \brief Returns true iff the Torch CUDA backend is available.
 * \author Ole Schuett
 ******************************************************************************/
bool torch_c_cuda_is_available() { return torch::cuda::is_available(); }

/*******************************************************************************
 * \brief Return the number of CUDA devices visible to Torch.
 ******************************************************************************/
int torch_c_cuda_device_count() {
  return torch::cuda::is_available() ? torch::cuda::device_count() : 0;
}

/*******************************************************************************
 * \brief Set whether to allow TF32.
 *        Needed due to changes in defaults from pytorch 1.7 to 1.11 to >=1.12
 *        See https://pytorch.org/docs/stable/notes/cuda.html
 * \author Gabriele Tocci
 ******************************************************************************/
void torch_c_allow_tf32(const bool allow_tf32) {

  at::globalContext().setAllowTF32CuBLAS(allow_tf32);
  at::globalContext().setAllowTF32CuDNN(allow_tf32);
}

/******************************************************************************
 * \brief Freeze the Torch model: generic optimization that speeds up model.
 *        See https://pytorch.org/docs/stable/generated/torch.jit.freeze.html
 * \author Gabriele Tocci
 ******************************************************************************/
void torch_c_model_freeze(torch_c_model_t *model) {

  *model = torch::jit::freeze(*model);
}

/*******************************************************************************
 * \brief Retrieves an int64 attribute. Must be called before model freeze.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_get_attr_int64(const torch_c_model_t *model, const char *key,
                                  int64_t *dest) {
  *dest = model->attr(key).toInt();
}

/*******************************************************************************
 * \brief Retrieves a double attribute. Must be called before model freeze.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_get_attr_double(const torch_c_model_t *model,
                                   const char *key, double *dest) {
  *dest = model->attr(key).toDouble();
}

/*******************************************************************************
 * \brief Retrieves a string attribute. Must be called before model freeze.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_get_attr_string(const torch_c_model_t *model,
                                   const char *key, char *dest) {
  const std::string &str = model->attr(key).toStringRef();
  assert(str.size() < 80); // default_string_length
  for (int i = 0; i < str.size(); i++) {
    dest[i] = str[i];
  }
}

/*******************************************************************************
 * \brief Retrieves a list attribute's size. Must be called before model freeze.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_get_attr_list_size(const torch_c_model_t *model,
                                      const char *key, int *size) {
  *size = model->attr(key).toList().size();
}

/*******************************************************************************
 * \brief Retrieves a single item from a string list attribute.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_get_attr_strlist(const torch_c_model_t *model,
                                    const char *key, const int index,
                                    char *dest) {
  const auto list = model->attr(key).toList();
  const std::string &str = list[index].toStringRef();
  assert(str.size() < 80); // default_string_length
  for (int i = 0; i < str.size(); i++) {
    dest[i] = str[i];
  }
}

#ifdef __cplusplus
}
#endif

#endif // defined(__LIBTORCH)

// EOF
