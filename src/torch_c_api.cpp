/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2024 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined(__LIBTORCH)

#include <torch/csrc/api/include/torch/cuda.h>
#include <torch/script.h>

typedef c10::Dict<std::string, torch::Tensor> torch_c_dict_t;
typedef torch::jit::Module torch_c_model_t;

/*******************************************************************************
 * \brief Internal helper for selecting the CUDA device when available.
 * \author Ole Schuett
 ******************************************************************************/
static torch::Device get_device() {
  return (torch::cuda::is_available()) ? torch::kCUDA : torch::kCPU;
}

/*******************************************************************************
 * \brief Internal helper for retrieving arrays from Torch dictionary.
 * \author Ole Schuett
 ******************************************************************************/
template <typename T>
static void torch_c_dict_get(const torch_c_dict_t *dict, const char *key,
                             const int ndims, int64_t sizes[], T **dest) {

  assert(dict->contains(key));
  const torch::Tensor tensor = dict->at(key).cpu();

  assert(tensor.ndimension() == ndims);
  int64_t size_flat = 1;
  for (int i = 0; i < ndims; i++) {
    sizes[i] = tensor.size(i);
    size_flat *= sizes[i];
  }
  *dest = (T *)malloc(size_flat * sizeof(T));

  const torch::Tensor tensor_flat = tensor.flatten();
  const auto accessor = tensor_flat.accessor<T, 1>();
  for (int i = 0; i < size_flat; i++) {
    (*dest)[i] = accessor[i];
  }
};

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 * \brief Inserts array of floats into Torch dictionary.
 *        The passed array has to outlive the dictionary!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_dict_insert_float(torch_c_dict_t *dict, const char *key,
                               const int ndims, const int64_t sizes[],
                               float source[]) {
  const auto options = torch::TensorOptions().dtype(torch::kFloat32);
  const auto sizes_ref = c10::IntArrayRef(sizes, ndims);
  const torch::Tensor tensor = torch::from_blob(source, sizes_ref, options);
  dict->insert(key, tensor.to(get_device()));
}

/*******************************************************************************
 * \brief Inserts array of int64s into Torch dictionary.
 *        The passed array has to outlive the dictionary!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_dict_insert_int64(torch_c_dict_t *dict, const char *key,
                               const int ndims, const int64_t sizes[],
                               int64_t source[]) {
  const auto options = torch::TensorOptions().dtype(torch::kInt64);
  const auto sizes_ref = c10::IntArrayRef(sizes, ndims);
  const torch::Tensor tensor = torch::from_blob(source, sizes_ref, options);
  dict->insert(key, tensor.to(get_device()));
}

/*******************************************************************************
 * \brief Inserts array of doubles into Torch dictionary.
 *        The passed array has to outlive the dictionary!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_dict_insert_double(torch_c_dict_t *dict, const char *key,
                                const int ndims, const int64_t sizes[],
                                double source[]) {
  const auto options = torch::TensorOptions().dtype(torch::kFloat64);
  const auto sizes_ref = c10::IntArrayRef(sizes, ndims);
  const torch::Tensor tensor = torch::from_blob(source, sizes_ref, options);
  dict->insert(key, tensor.to(get_device()));
}

/*******************************************************************************
 * \brief Retrieves array of floats from Torch dictionary.
 *        The returned array has to be deallocated by caller!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_dict_get_float(const torch_c_dict_t *dict, const char *key,
                            const int ndims, int64_t sizes[], float **dest) {

  torch_c_dict_get<float>(dict, key, ndims, sizes, dest);
}

/*******************************************************************************
 * \brief Retrieves array of int64s from Torch dictionary.
 *        The returned array has to be deallocated by caller!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_dict_get_int64(const torch_c_dict_t *dict, const char *key,
                            const int ndims, int64_t sizes[], int64_t **dest) {

  torch_c_dict_get<int64_t>(dict, key, ndims, sizes, dest);
}

/*******************************************************************************
 * \brief Retrieves array of doubles from Torch dictionary.
 *        The returned array has to be deallocated by caller!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_dict_get_double(const torch_c_dict_t *dict, const char *key,
                             const int ndims, int64_t sizes[], double **dest) {

  torch_c_dict_get<double>(dict, key, ndims, sizes, dest);
}

/*******************************************************************************
 * \brief Creates an empty Torch dictionary.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_dict_create(torch_c_dict_t **dict_out) {
  assert(*dict_out == NULL);
  *dict_out = new c10::Dict<std::string, torch::Tensor>();
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

  torch::jit::Module *model = new torch::jit::Module();
  *model = torch::jit::load(filename, get_device());
  model->eval();

  assert(*model_out == NULL);
  *model_out = model;
}

/*******************************************************************************
 * \brief Evaluates the given Torch model.
 *        In Torch lingo this operation is called forward().
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_eval(torch_c_model_t *model, const torch_c_dict_t *inputs,
                        torch_c_dict_t *outputs) {

  auto untyped_output = model->forward({*inputs}).toGenericDict();

  outputs->clear();
  for (const auto &entry : untyped_output) {
    outputs->insert(entry.key().toStringView(), entry.value().toTensor());
  }
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
  *length = content_str.length();
  *content = (char *)malloc(content_str.length() + 1); // +1 for null terminator
  strcpy(*content, content_str.c_str());
}

/*******************************************************************************
 * \brief Returns true iff the Torch CUDA backend is available.
 * \author Ole Schuett
 ******************************************************************************/
bool torch_c_cuda_is_available() { return torch::cuda::is_available(); }

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

#ifdef __cplusplus
}
#endif

#endif // defined(__LIBTORCH)

// EOF
