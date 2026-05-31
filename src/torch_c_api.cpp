/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined(__LIBTORCH)

#include <c10/core/DeviceGuard.h>
#include <torch/csrc/api/include/torch/cuda.h>
#include <torch/script.h>

#include <cctype>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

typedef torch::Tensor torch_c_tensor_t;
typedef c10::Dict<std::string, torch::Tensor> torch_c_dict_t;
typedef torch::jit::Module torch_c_model_t;

/*******************************************************************************
 * \brief Internal helper for selecting the CUDA device when available.
 * \author Ole Schuett
 ******************************************************************************/
static bool use_cuda_if_available = true;
static int selected_cuda_device = -1;
static int selected_cuda_visible_device = -1;

static torch::Device get_device() {
  if (!use_cuda_if_available || !torch::cuda::is_available()) {
    return torch::kCPU;
  }
  const auto device_count = torch::cuda::device_count();
  if (device_count <= 0) {
    return torch::kCPU;
  }
  const int device = (selected_cuda_device >= 0) ? selected_cuda_device : 0;
  return torch::Device(torch::kCUDA, device);
}

static torch::Device get_device_with_guard(c10::OptionalDeviceGuard &guard) {
  const auto device = get_device();
  if (device.is_cuda()) {
    guard.reset_device(device);
  }
  return device;
}

static int local_rank_from_env() {
  const char *rank_env_names[] = {
      "OMPI_COMM_WORLD_LOCAL_RANK", "PMI_LOCAL_RANK",  "PMIX_LOCAL_RANK",
      "MV2_COMM_WORLD_LOCAL_RANK",  "MPI_LOCALRANKID", "SLURM_LOCALID"};
  for (const auto *name : rank_env_names) {
    const char *value = std::getenv(name);
    if (value == nullptr || value[0] == '\0') {
      continue;
    }
    char *end = nullptr;
    const long rank = std::strtol(value, &end, 10);
    if (end != value && rank >= 0) {
      return static_cast<int>(rank);
    }
  }
  return 0;
}

static bool parse_nonnegative_int(const std::string &text, int &value) {
  if (text.empty()) {
    return false;
  }
  char *end = nullptr;
  const long parsed = std::strtol(text.c_str(), &end, 10);
  if (end != text.c_str() + text.size() || parsed < 0) {
    return false;
  }
  value = static_cast<int>(parsed);
  return true;
}

static std::string trim_string(const std::string &text) {
  std::size_t begin = 0;
  while (begin < text.size() &&
         std::isspace(static_cast<unsigned char>(text[begin]))) {
    ++begin;
  }
  std::size_t end = text.size();
  while (end > begin &&
         std::isspace(static_cast<unsigned char>(text[end - 1]))) {
    --end;
  }
  return text.substr(begin, end - begin);
}

static std::vector<std::string> cuda_visible_device_tokens() {
  std::vector<std::string> tokens;
  const char *visible = std::getenv("CUDA_VISIBLE_DEVICES");
  if (visible == nullptr || visible[0] == '\0') {
    return tokens;
  }
  std::string list(visible);
  std::size_t begin = 0;
  while (begin <= list.size()) {
    const std::size_t comma = list.find(',', begin);
    const std::size_t end = (comma == std::string::npos) ? list.size() : comma;
    const std::string token = trim_string(list.substr(begin, end - begin));
    if (!token.empty()) {
      tokens.push_back(token);
    }
    if (comma == std::string::npos) {
      break;
    }
    begin = comma + 1;
  }
  if (tokens.size() == 1 && tokens[0] == "-1") {
    tokens.clear();
  }
  return tokens;
}

static int cuda_visible_device_for_log(const int logical_device) {
  const std::vector<std::string> tokens = cuda_visible_device_tokens();
  if (tokens.empty()) {
    return logical_device;
  }
  if (logical_device < 0 || logical_device >= static_cast<int>(tokens.size())) {
    return -1;
  }
  int visible_device = -1;
  parse_nonnegative_int(tokens[logical_device], visible_device);
  return visible_device;
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

static torch::jit::Module load_module_for_device(
    const char *filename, const torch::Device &device,
    std::unordered_map<std::string, std::string> &extra_files) {
  if (can_load_directly_to_device(device)) {
    return torch::jit::load(filename, device, extra_files);
  }
  auto model = torch::jit::load(filename, torch::kCPU, extra_files);
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

static void grad_select_to_array_double(const torch_c_tensor_t *tensor,
                                        const torch_c_tensor_t *indices,
                                        const int ndims, const int64_t sizes[],
                                        double target[]) {
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  const torch::Tensor grad = tensor->grad();
  assert(grad.defined());
  assert(grad.scalar_type() == torch::kFloat64);
  assert(indices->scalar_type() == torch::kInt64);
  assert(indices->ndimension() == 1);

  const auto selected = grad.index_select(ndims - 1, indices->to(grad.device()))
                            .cpu()
                            .contiguous();
  assert(selected.ndimension() == ndims);
  for (int i = 0; i < ndims; i++) {
    assert(selected.size(i) == sizes[i]);
  }
  std::memcpy(target, selected.data_ptr<double>(), selected.nbytes());
}

static void grad_to_array_double(const torch_c_tensor_t *tensor,
                                 const int ndims, const int64_t sizes[],
                                 double target[]) {
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  const torch::Tensor grad = tensor->grad();
  assert(grad.defined());
  assert(grad.scalar_type() == torch::kFloat64);

  const auto selected = grad.cpu().contiguous();
  assert(selected.ndimension() == ndims);
  for (int i = 0; i < ndims; i++) {
    assert(selected.size(i) == sizes[i]);
  }
  std::memcpy(target, selected.data_ptr<double>(), selected.nbytes());
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
 * \brief Copies selected rows from a double tensor's gradient to a host array.
 ******************************************************************************/
void torch_c_tensor_grad_select_to_array_double(const torch_c_tensor_t *tensor,
                                                const torch_c_tensor_t *indices,
                                                const int ndims,
                                                const int64_t sizes[],
                                                double target[]) {
  grad_select_to_array_double(tensor, indices, ndims, sizes, target);
}

/*******************************************************************************
 * \brief Copies a double tensor's gradient to a host array.
 ******************************************************************************/
void torch_c_tensor_grad_to_array_double(const torch_c_tensor_t *tensor,
                                         const int ndims, const int64_t sizes[],
                                         double target[]) {
  grad_to_array_double(tensor, ndims, sizes, target);
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
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  tensor->backward(*outer_grad);
}

/*******************************************************************************
 * \brief Runs autograd on a scalar Torch tensor.
 ******************************************************************************/
void torch_c_tensor_backward_scalar(const torch_c_tensor_t *tensor) {
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
 * \brief Loads a Torch model and two extra-file metadata entries in one pass.
 ******************************************************************************/
void torch_c_model_load_metadata2(torch_c_model_t **model_out,
                                  const char *filename, const char *key1,
                                  const char *key2, char **content1,
                                  int *length1, char **content2, int *length2) {
  assert(*model_out == NULL);
  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  set_jit_fusion_strategy();
  std::unordered_map<std::string, std::string> extra_files = {{key1, ""},
                                                              {key2, ""}};
  torch::jit::Module *model = new torch::jit::Module();
  *model = load_module_for_device(filename, device, extra_files);
  model->eval(); // Set to evaluation mode to disable gradients, drop-out, etc.
  *model_out = model;
  copy_string_to_c_buffer(extra_files[key1], content1, length1);
  copy_string_to_c_buffer(extra_files[key2], content2, length2);
}

/*******************************************************************************
 * \brief Evaluates the given Torch model.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_forward(torch_c_model_t *model, const torch_c_dict_t *inputs,
                           torch_c_dict_t *outputs) {

  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  auto untyped_output = model->forward({*inputs}).toGenericDict();
  outputs->clear();
  for (const auto &entry : untyped_output) {
    outputs->insert(entry.key().toStringView(), entry.value().toTensor());
  }
}

/*******************************************************************************
 * \brief Evaluates a TorchScript model method expecting keyword argument "mol".
 ******************************************************************************/
void torch_c_model_forward_mol_tensor(torch_c_model_t *model,
                                      const char *method_name,
                                      const torch_c_dict_t *inputs,
                                      torch_c_tensor_t **output) {

  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  std::vector<c10::IValue> args;
  std::unordered_map<std::string, c10::IValue> kwargs;
  kwargs["mol"] = *inputs;
  *output = new torch_c_tensor_t(
      model->get_method(method_name)(args, kwargs).toTensor());
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
 * \brief Returns the number of CUDA devices visible to Torch.
 ******************************************************************************/
int torch_c_cuda_device_count() {
  return torch::cuda::is_available()
             ? static_cast<int>(torch::cuda::device_count())
             : 0;
}

/*******************************************************************************
 * \brief Select the CUDA device used by subsequent Torch wrapper calls.
 *
 * A negative device maps the MPI-local rank to a visible CUDA device. This
 *keeps multi-rank CP2K jobs from collapsing all native SKALA Torch work onto
 *device 0 when CUDA_VISIBLE_DEVICES exposes more than one GPU.
 ******************************************************************************/
int torch_c_cuda_select_device(const int requested_device) {
  const int device_count = torch_c_cuda_device_count();
  if (device_count <= 0) {
    selected_cuda_device = -1;
    selected_cuda_visible_device = -1;
    return -1;
  }
  int device = requested_device;
  if (requested_device < 0) {
    device = local_rank_from_env() % device_count;
  }
  if (device < 0 || device >= device_count) {
    return -2;
  }
  selected_cuda_device = device;
  selected_cuda_visible_device =
      cuda_visible_device_for_log(selected_cuda_device);
  return selected_cuda_device;
}

/*******************************************************************************
 * \brief Returns the physical CUDA device made visible to Torch, if known.
 ******************************************************************************/
int torch_c_cuda_visible_device() { return selected_cuda_visible_device; }

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
