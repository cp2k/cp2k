/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined(__LIBTORCH)

#include <ATen/Parallel.h>
#include <algorithm>
#include <atomic>
#if defined(__LIBTORCH_CUDA)
#include <ATen/cuda/CUDAContextLight.h>
#include <c10/cuda/CUDAAllocatorConfig.h>
#include <c10/cuda/CUDACachingAllocator.h>
#include <c10/cuda/CUDAStream.h>
#endif
#include <c10/core/DeviceGuard.h>
#include <torch/csrc/api/include/torch/cuda.h>
#include <torch/csrc/autograd/autograd.h>
#include <torch/csrc/jit/passes/freeze_module.h>
#include <torch/csrc/jit/passes/inliner.h>
#include <torch/csrc/jit/passes/subgraph_rewrite.h>
#include <torch/csrc/jit/runtime/graph_executor.h>
#include <torch/script.h>

#include "offload/offload_library.h"

#if defined(__OPENBLAS)
#include <cblas.h>
#endif

#include <cfenv>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <limits>
#include <mutex>
#include <optional>
#include <set>
#include <string>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#if defined(__OPENBLAS)
// PyTorch's oneMKL batch ABI is not compatible with OpenBLAS's same-named
// entry points. Expand grouped GEMMs into the portable CBLAS interface.
extern "C" void
cblas_sgemm_batch(const CBLAS_ORDER order, const CBLAS_TRANSPOSE *trans_a,
                  const CBLAS_TRANSPOSE *trans_b, const int *m, const int *n,
                  const int *k, const float *alpha, const float **a,
                  const int *lda, const float **b, const int *ldb,
                  const float *beta, float **c, const int *ldc,
                  const int group_count, const int *group_size) {
  int offset = 0;
  for (int group = 0; group < group_count; ++group) {
    for (int operation = 0; operation < group_size[group]; ++operation) {
      const int index = offset + operation;
      cblas_sgemm(order, trans_a[group], trans_b[group], m[group], n[group],
                  k[group], alpha[group], a[index], lda[group], b[index],
                  ldb[group], beta[group], c[index], ldc[group]);
    }
    offset += group_size[group];
  }
}

extern "C" void
cblas_dgemm_batch(const CBLAS_ORDER order, const CBLAS_TRANSPOSE *trans_a,
                  const CBLAS_TRANSPOSE *trans_b, const int *m, const int *n,
                  const int *k, const double *alpha, const double **a,
                  const int *lda, const double **b, const int *ldb,
                  const double *beta, double **c, const int *ldc,
                  const int group_count, const int *group_size) {
  int offset = 0;
  for (int group = 0; group < group_count; ++group) {
    for (int operation = 0; operation < group_size[group]; ++operation) {
      const int index = offset + operation;
      cblas_dgemm(order, trans_a[group], trans_b[group], m[group], n[group],
                  k[group], alpha[group], a[index], lda[group], b[index],
                  ldb[group], beta[group], c[index], ldc[group]);
    }
    offset += group_size[group];
  }
}
#endif

typedef torch::Tensor torch_c_tensor_t;
typedef c10::Dict<std::string, torch::Tensor> torch_c_dict_t;

// Wraps a TorchScript module together with its resolved method handles so the
// forward path does not have to lock and look up a global method cache.
struct torch_c_model_t {
  torch::jit::Module module;
  mutable std::mutex method_mutex;
  bool has_primary_method = false;
  std::string primary_method_name;
  std::optional<torch::jit::Method> primary_method;
  std::unordered_map<std::string, torch::jit::Method> extra_methods;
};

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
 * \brief Aborts the run with a message instead of letting a C++ exception
 *        unwind across the BIND(C) boundary.
 ******************************************************************************/
[[noreturn]] static void torch_c_api_abort(const std::string &message) {
  std::fprintf(stderr, "TORCH| %s\n", message.c_str());
  std::fflush(stderr);
  std::abort();
}

static void torch_c_check(const bool condition, const char *message) {
  if (!condition) {
    torch_c_api_abort(message);
  }
}

// Pinned-host budget (per rank), overridable via CP2K_TORCH_PINNED_MB.
static int64_t torch_pinned_budget_bytes() {
  static const int64_t budget = []() -> int64_t {
    int64_t mb = 256;
    const char *env = std::getenv("CP2K_TORCH_PINNED_MB");
    if (env != nullptr && env[0] != '\0') {
      char *end = nullptr;
      const long value = std::strtol(env, &end, 10);
      if (end != env && value > 0) {
        mb = static_cast<int64_t>(value);
      }
    }
    return mb * 1024 * 1024;
  }();
  return budget;
}

static void warn_pageable_pinned_once(const char *what) {
  static std::once_flag flag;
  std::call_once(flag, [what]() {
    std::fprintf(
        stderr,
        "TORCH| %s exceeds CP2K_TORCH_PINNED_MB; using pageable memory.\n",
        what);
  });
}

#define TORCH_C_TRY try {
#define TORCH_C_CATCH(NAME)                                                    \
  }                                                                            \
  catch (const std::exception &e) {                                            \
    torch_c_api_abort(std::string(NAME) + ": " + e.what());                    \
  }                                                                            \
  catch (...) {                                                                \
    torch_c_api_abort(std::string(NAME) + ": unknown exception");              \
  }

/*******************************************************************************
 * \brief Internal helper for selecting the CUDA device when available.
 * \author Ole Schuett
 ******************************************************************************/
static std::atomic<bool> use_cuda_if_available{true};

struct model_metadata_entry {
  std::string stamp;
  std::unordered_map<std::string, std::string> files;
};
static std::mutex model_metadata_mutex;
static std::unordered_map<std::string, model_metadata_entry>
    model_metadata_cache;

static std::set<int> pending_grad_readback_devices;
static std::set<const torch_c_tensor_t *> pending_grad_tensors;
static std::mutex pending_grad_mutex;

static void record_pending_grads(const torch_c_tensor_t *const tensors[],
                                 const int count, const int device) {
  std::lock_guard<std::mutex> lock(pending_grad_mutex);
  pending_grad_readback_devices.insert(device);
  for (int i = 0; i < count; ++i)
    pending_grad_tensors.insert(tensors[i]);
}

static std::set<int> take_pending_devices() {
  std::lock_guard<std::mutex> lock(pending_grad_mutex);
  std::set<int> devices = std::move(pending_grad_readback_devices);
  pending_grad_readback_devices.clear();
  pending_grad_tensors.clear();
  return devices;
}

static bool pending_grad_copy(const torch_c_tensor_t *tensor) {
  std::lock_guard<std::mutex> lock(pending_grad_mutex);
  return pending_grad_tensors.count(tensor) != 0;
}

static void enable_expandable_cuda_segments_if_unconfigured() {
#if defined(__LIBTORCH_CUDA)
  static const bool initialized = []() {
    const char *legacy_config = std::getenv("PYTORCH_CUDA_ALLOC_CONF");
    const char *config = std::getenv("PYTORCH_ALLOC_CONF");
    if ((legacy_config == nullptr || legacy_config[0] == '\0') &&
        (config == nullptr || config[0] == '\0')) {
      c10::cuda::CUDACachingAllocator::setAllocatorSettings(
          "expandable_segments:True");
    }
    return true;
  }();
  (void)initialized;
#endif
}

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
  static std::once_flag init_flag;
  std::call_once(init_flag, []() {
    int num_threads = 0;
    if (get_positive_int_env("CP2K_TORCH_NUM_THREADS", num_threads)) {
      at::set_num_threads(num_threads);
    }
    if (get_positive_int_env("CP2K_TORCH_NUM_INTEROP_THREADS", num_threads)) {
      at::set_num_interop_threads(num_threads);
    }
  });
}

extern "C" void torch_c_set_num_threads(const int num_threads) {
  TORCH_C_TRY
  int env_threads = 0;
  if (num_threads <= 0)
    return;
  if (get_positive_int_env("CP2K_TORCH_NUM_THREADS", env_threads))
    return;
  at::set_num_threads(num_threads);
  TORCH_C_CATCH("torch_c_set_num_threads")
}

extern "C" int torch_c_get_num_threads(void) {
  int num_threads = 0;
  TORCH_C_TRY
  initialize_torch_threads_from_env();
  num_threads = at::get_num_threads();
  TORCH_C_CATCH("torch_c_get_num_threads")
  return num_threads;
}

static bool cuda_device_available(int &device_count) {
  static std::once_flag probe_flag;
  static bool available = false;
  static int count = 0;
  std::call_once(probe_flag, []() {
    available = torch::cuda::is_available();
    count = available ? torch::cuda::device_count() : 0;
    available = available && count > 0;
  });
  device_count = count;
  return available;
}

static torch::Device get_device() {
  initialize_torch_threads_from_env();
  int device_count = 0;
  if (!use_cuda_if_available.load() || !cuda_device_available(device_count)) {
    return torch::kCPU;
  }
  enable_expandable_cuda_segments_if_unconfigured();
  const int chosen_device = offload_get_chosen_device();
  const int device = (chosen_device >= 0) ? chosen_device : 0;
  torch_c_check(device < device_count,
                "Torch CUDA device index is out of range");
  return torch::Device(torch::kCUDA, device);
}

static torch::Device get_device_with_guard(c10::OptionalDeviceGuard &guard) {
  const auto device = get_device();
  if (device.is_cuda()) {
    guard.reset_device(device);
  }
  return device;
}

// Wait for all deferred device-to-host readbacks recorded by the grad-pool path
// before their pinned buffers may be reused or freed.
static void synchronize_pending_grad_copies() {
#if defined(__LIBTORCH_CUDA)
  if (!torch::cuda::is_available()) {
    take_pending_devices();
    return;
  }
  c10::OptionalDeviceGuard guard;
  const auto active = get_device_with_guard(guard);
  c10::cuda::getCurrentCUDAStream().synchronize();
  for (const int device_index : take_pending_devices()) {
    if (active.is_cuda() && static_cast<int>(active.index()) == device_index) {
      continue; // Already synchronized above.
    }
    c10::DeviceGuard device_guard(torch::Device(
        torch::kCUDA, static_cast<c10::DeviceIndex>(device_index)));
    c10::cuda::getCurrentCUDAStream().synchronize();
  }
#else
  take_pending_devices();
#endif
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
  torch_c_check(content != nullptr && length != nullptr,
                "Torch string output must not be NULL");
  torch_c_check(source.size() <= static_cast<size_t>(INT_MAX),
                "Torch string exceeds the maximum supported length");
  char *buffer = static_cast<char *>(std::malloc(source.size() + 1));
  torch_c_check(buffer != nullptr, "Torch string allocation failed");
  std::memcpy(buffer, source.c_str(), source.size() + 1);
  *content = buffer;
  *length = static_cast<int>(source.size());
}

static bool can_load_directly_to_device(const torch::Device &device) {
  return !device.is_cuda() || device.index() == 0 ||
         torch::cuda::device_count() == 1;
}

static torch::jit::Module load_module_for_device(
    const char *filename, const torch::Device &device,
    std::unordered_map<std::string, std::string> *extra_files = nullptr,
    const bool load_on_cpu = false) {
  if (!load_on_cpu && can_load_directly_to_device(device)) {
    if (extra_files != nullptr) {
      return torch::jit::load(filename, device, *extra_files);
    }
    return torch::jit::load(filename, device);
  }
  auto model = (extra_files != nullptr)
                   ? torch::jit::load(filename, torch::kCPU, *extra_files)
                   : torch::jit::load(filename, torch::kCPU);
  if (!load_on_cpu) {
    model.to(device);
  }
  return model;
}

static void remap_device_constants(torch::jit::Block *block,
                                   const torch::Device &device) {
  for (torch::jit::Node *node : block->nodes()) {
    if (node->kind() == torch::jit::prim::Constant &&
        node->outputs().size() == 1 &&
        node->output()->type()->kind() == c10::TypeKind::DeviceObjType &&
        node->hasAttribute(torch::jit::attr::value) &&
        node->kindOf(torch::jit::attr::value) == torch::jit::AttributeKind::s) {
      node->s_(torch::jit::attr::value, device.str());
    }
    bool has_tensor_input = false;
    for (const torch::jit::Value *input : node->inputs()) {
      has_tensor_input |= input->type()->kind() == c10::TypeKind::TensorType;
    }
    const auto *schema = node->maybeSchema();
    if (!has_tensor_input && schema != nullptr) {
      const auto &arguments = schema->arguments();
      for (size_t i = 0; i < arguments.size() && i < node->inputs().size();
           ++i) {
        const auto input_value = torch::jit::toIValue(node->input(i));
        if (arguments[i].name() == "device" && input_value.has_value() &&
            input_value->isNone()) {
          torch::jit::WithInsertPoint insertion_guard(node);
          auto *device_value =
              node->owningGraph()->insertConstant(torch::jit::IValue(device));
          node->replaceInput(i, device_value);
        }
      }
    }
    for (torch::jit::Block *nested_block : node->blocks()) {
      remap_device_constants(nested_block, device);
    }
  }
}

static void remap_model_device_constants(torch::jit::Module &model,
                                         const torch::Device &device) {
  for (const auto &attribute : model.named_attributes(false)) {
    const auto &value = attribute.value;
    if (value.isTensor()) {
      model.setattr(attribute.name, value.toTensor().to(device));
    } else if (value.isTensorList()) {
      auto tensors = value.toTensorList();
      for (size_t i = 0; i < tensors.size(); ++i) {
        tensors.set(i, tensors.get(i).to(device));
      }
      model.setattr(attribute.name, tensors);
    }
  }
  for (const auto &method : model.get_methods()) {
    remap_device_constants(method.graph()->block(), device);
  }
  for (auto child : model.children()) {
    remap_model_device_constants(child, device);
  }
}

static void collect_block_device_constants(torch::jit::Block *block,
                                           bool &has_cuda, bool &has_cpu) {
  for (torch::jit::Node *node : block->nodes()) {
    if (node->kind() == torch::jit::prim::Constant &&
        node->outputs().size() == 1 &&
        node->output()->type()->kind() == c10::TypeKind::DeviceObjType &&
        node->hasAttribute(torch::jit::attr::value) &&
        node->kindOf(torch::jit::attr::value) == torch::jit::AttributeKind::s) {
      const std::string &value = node->s(torch::jit::attr::value);
      if (value.rfind("cuda", 0) == 0) {
        has_cuda = true;
      } else if (value.rfind("cpu", 0) == 0) {
        has_cpu = true;
      }
    }
    for (torch::jit::Block *nested_block : node->blocks()) {
      collect_block_device_constants(nested_block, has_cuda, has_cpu);
    }
  }
}

static void collect_model_device_constants(torch::jit::Module &model,
                                           bool &has_cuda, bool &has_cpu) {
  for (const auto &method : model.get_methods()) {
    collect_block_device_constants(method.graph()->block(), has_cuda, has_cpu);
  }
  for (auto child : model.children()) {
    collect_model_device_constants(child, has_cuda, has_cpu);
  }
}

static void promote_float32_constants(torch::jit::Block *block) {
  for (torch::jit::Node *node : block->nodes()) {
    if (node->kind() == torch::jit::prim::Constant &&
        node->hasAttribute(torch::jit::attr::value) &&
        node->kindOf(torch::jit::attr::value) == torch::jit::AttributeKind::t) {
      const auto tensor = node->t(torch::jit::attr::value);
      if (tensor.defined() && tensor.scalar_type() == torch::kFloat32) {
        const auto promoted = tensor.to(torch::kFloat64);
        node->t_(torch::jit::attr::value, promoted);
        node->output()->setType(c10::TensorType::create(promoted));
      }
    }
    const auto *schema = node->maybeSchema();
    if (schema != nullptr) {
      const auto &arguments = schema->arguments();
      for (size_t i = 0; i < arguments.size() && i < node->inputs().size();
           ++i) {
        const auto value = torch::jit::toIValue(node->input(i));
        if (arguments[i].name() == "dtype" && value.has_value() &&
            value->isInt() &&
            value->toInt() == static_cast<int64_t>(torch::kFloat32)) {
          // A dtype constant can also be used as a dimension or an index.
          // Replace this edge only; leave other uses and integer tensors
          // intact.
          torch::jit::WithInsertPoint insertion_guard(node);
          auto *dtype = node->owningGraph()->insertConstant(
              static_cast<int64_t>(torch::kFloat64));
          node->replaceInput(i, dtype);
        }
      }
    }
    for (torch::jit::Block *nested : node->blocks()) {
      promote_float32_constants(nested);
    }
  }
}

static void promote_model_float32(torch::jit::Module &model) {
  for (const auto &attribute : model.named_attributes(false)) {
    const auto &value = attribute.value;
    if (value.isTensor() && value.toTensor().defined() &&
        value.toTensor().scalar_type() == torch::kFloat32) {
      model.setattr(attribute.name, value.toTensor().to(torch::kFloat64));
    } else if (value.isTensorList()) {
      auto tensors = value.toTensorList();
      for (size_t i = 0; i < tensors.size(); ++i) {
        if (tensors.get(i).defined() &&
            tensors.get(i).scalar_type() == torch::kFloat32) {
          tensors.set(i, tensors.get(i).to(torch::kFloat64));
        }
      }
      model.setattr(attribute.name, tensors);
    }
  }
  for (const auto &method : model.get_methods()) {
    torch::jit::Inline(*method.graph());
    promote_float32_constants(method.graph()->block());
  }
  for (auto child : model.children()) {
    promote_model_float32(child);
  }
}

static void initialize_cuda_linalg(const torch::Device &device) {
  static std::mutex mutex;
  static std::set<int> initialized_devices;
  std::lock_guard<std::mutex> lock(mutex);
  if (!initialized_devices.insert(device.index()).second) {
    return;
  }
  const auto options =
      torch::TensorOptions().dtype(torch::kFloat32).device(device);
  static_cast<void>(at::linalg_eigh(torch::eye(1, options)));
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

struct grad_pinned_entry {
  torch::Tensor tensor;
  uint64_t last_use = 0;
  bool pinned = false;
  int64_t bytes = 0;
};

static std::mutex grad_pinned_mutex;
static std::unordered_map<const torch_c_tensor_t *, grad_pinned_entry>
    *grad_pinned_pool = nullptr;
static int64_t grad_pinned_bytes = 0;
static uint64_t grad_pinned_clock = 0;

// A returned gradient wrapper shares storage with the pool entry, so it can
// outlive that entry. Keep the bytes of such aliased pinned buffers accounted
// until the last reference to the storage disappears instead of subtracting
// them when the map entry is erased.
struct retired_pinned_entry {
  c10::weak_intrusive_ptr<c10::TensorImpl> weak;
  int64_t bytes = 0;
};
static std::vector<retired_pinned_entry> grad_pinned_retired;

// Release the bytes of retired pinned buffers whose storage has been freed.
// Must be called with grad_pinned_mutex held.
static void reclaim_retired_pinned_locked() {
  if (grad_pinned_retired.empty())
    return;
  auto out = grad_pinned_retired.begin();
  for (auto it = grad_pinned_retired.begin(); it != grad_pinned_retired.end();
       ++it) {
    if (it->weak.expired()) {
      grad_pinned_bytes -= it->bytes;
    } else {
      if (out != it) {
        *out = std::move(*it);
      }
      ++out;
    }
  }
  grad_pinned_retired.erase(out, grad_pinned_retired.end());
}

// Drop one pool entry. Must be called with grad_pinned_mutex held. Pinned bytes
// are only released when erasing drops the last reference to the storage;
// otherwise a caller-held gradient still keeps the allocation alive.
static void drop_grad_pool_entry_locked(
    std::unordered_map<const torch_c_tensor_t *, grad_pinned_entry>::iterator
        it) {
  if (it->second.pinned) {
    const auto &t = it->second.tensor;
    if (!t.defined() || t.getIntrusivePtr().use_count() == 1) {
      grad_pinned_bytes -= it->second.bytes;
    } else {
      grad_pinned_retired.push_back(retired_pinned_entry{
          c10::weak_intrusive_ptr<c10::TensorImpl>(t.getIntrusivePtr()),
          it->second.bytes});
    }
  }
  grad_pinned_pool->erase(it);
}

static void erase_grad_pool_entry(const torch_c_tensor_t *tensor) {
  if (tensor == nullptr)
    return;
  std::lock_guard<std::mutex> lock(grad_pinned_mutex);
  reclaim_retired_pinned_locked();
  if (grad_pinned_pool != nullptr) {
    const auto it = grad_pinned_pool->find(tensor);
    if (it != grad_pinned_pool->end()) {
      if (pending_grad_copy(tensor)) {
        // A deferred readback still targets this pinned buffer. Wait for it so
        // the buffer is not released mid-copy. Lock order matches
        // evict_grad_pool_locked: grad_pinned_mutex -> pending_grad_mutex.
        synchronize_pending_grad_copies();
      }
      drop_grad_pool_entry_locked(it);
    }
  }
}

// Evict idle, non-pending entries (oldest first) until `need` more bytes fit in
// the pinned budget. Must be called with grad_pinned_mutex held.
static void evict_grad_pool_locked(const int64_t need) {
  if (grad_pinned_pool == nullptr) {
    return;
  }
  reclaim_retired_pinned_locked();
  const int64_t budget = torch_pinned_budget_bytes();
  while (grad_pinned_bytes + need > budget && !grad_pinned_pool->empty()) {
    auto victim = grad_pinned_pool->end();
    uint64_t oldest = UINT64_MAX;
    for (auto it = grad_pinned_pool->begin(); it != grad_pinned_pool->end();
         ++it) {
      const auto &entry = it->second;
      if (!entry.tensor.defined()) {
        continue;
      }
      if (entry.tensor.getIntrusivePtr().use_count() != 1) {
        continue; // still referenced by a caller
      }
      if (pending_grad_copy(it->first)) {
        continue; // deferred readback in flight
      }
      if (entry.last_use < oldest) {
        oldest = entry.last_use;
        victim = it;
      }
    }
    if (victim == grad_pinned_pool->end()) {
      break; // nothing safe to evict
    }
    drop_grad_pool_entry_locked(victim);
  }
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

static bool prepare_target_tensor(torch_c_tensor_t **tensor,
                                  const bool req_grad,
                                  const torch::Device &device, const int ndims,
                                  const int64_t sizes[]) {
  const auto sizes_ref = c10::IntArrayRef(sizes, ndims);
  const bool requires_grad_changed = *tensor != nullptr &&
                                     (*tensor)->defined() &&
                                     (*tensor)->requires_grad() != req_grad;
  const bool reallocated =
      requires_grad_changed ||
      !tensor_matches(*tensor, torch::kFloat64, device, ndims, sizes);
  if (reallocated) {
    erase_grad_pool_entry(*tensor);
    delete (*tensor);
    const auto opts =
        torch::TensorOptions().dtype(torch::kFloat64).device(device);
    *tensor = new torch_c_tensor_t(torch::empty(sizes_ref, opts).detach());
  }
  return reallocated;
}

/*******************************************************************************
 * \brief Internal helper for getting the data_ptr and sizes of a Torch tensor.
 * \author Ole Schuett
 ******************************************************************************/
static void *get_data_ptr(const torch_c_tensor_t *tensor,
                          const torch::Dtype dtype, const int ndims,
                          int64_t sizes[]) {
  torch_c_check(tensor->scalar_type() == dtype,
                "Torch tensor has an unexpected dtype");
  torch_c_check(tensor->ndimension() == ndims,
                "Torch tensor has an unexpected number of dimensions");
  for (int i = 0; i < ndims; i++) {
    sizes[i] = tensor->size(i);
  }

  torch_c_check(tensor->is_contiguous(), "Torch tensor is not contiguous");
  torch_c_check(tensor->is_cpu(), "Torch tensor is not on the CPU");
  return tensor->data_ptr();
};
/*******************************************************************************
 * \brief Copies a model attribute string into a Fortran character buffer of the
 *        given length, blank-filling and rejecting strings that do not fit.
 ******************************************************************************/
static void copy_attr_string_to_fortran(const std::string &source, char *dest,
                                        const int dest_len) {
  torch_c_check(dest != nullptr && dest_len >= 0,
                "Torch attribute destination must not be NULL");
  if (static_cast<int>(source.size()) > dest_len) {
    torch_c_api_abort("Torch attribute string exceeds destination length");
  }
  std::memset(dest, ' ', dest_len);
  std::memcpy(dest, source.data(), source.size());
}

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
  TORCH_C_TRY
  *tensor = tensor_from_array(torch::kInt32, req_grad, ndims, sizes, source);
  TORCH_C_CATCH("torch_c_tensor_from_array_int32")
}

/*******************************************************************************
 * \brief Creates a Torch tensor from an array of floats.
 *        The passed array has to outlive the tensor!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_from_array_float(torch_c_tensor_t **tensor,
                                     const bool req_grad, const int ndims,
                                     const int64_t sizes[], float source[]) {
  TORCH_C_TRY
  *tensor = tensor_from_array(torch::kFloat32, req_grad, ndims, sizes, source);
  TORCH_C_CATCH("torch_c_tensor_from_array_float")
}

/*******************************************************************************
 * \brief Creates a Torch tensor from an array of int64s.
 *        The passed array has to outlive the tensor!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_from_array_int64(torch_c_tensor_t **tensor,
                                     const bool req_grad, const int ndims,
                                     const int64_t sizes[], int64_t source[]) {
  TORCH_C_TRY
  *tensor = tensor_from_array(torch::kInt64, req_grad, ndims, sizes, source);
  TORCH_C_CATCH("torch_c_tensor_from_array_int64")
}

/*******************************************************************************
 * \brief Creates a Torch tensor from an array of doubles.
 *        The passed array has to outlive the tensor!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_from_array_double(torch_c_tensor_t **tensor,
                                      const bool req_grad, const int ndims,
                                      const int64_t sizes[], double source[]) {
  TORCH_C_TRY
  *tensor = tensor_from_array(torch::kFloat64, req_grad, ndims, sizes, source);
  TORCH_C_CATCH("torch_c_tensor_from_array_double")
}

/*******************************************************************************
 * \brief Releases a string returned from the Torch C API.
 ******************************************************************************/
void torch_c_free_string(char *content) { free(content); }

/*******************************************************************************
 * \brief Reusable pinned host staging buffers for batched host-to-device
 * copies.
 ******************************************************************************/
struct staging_slot {
  torch::Tensor tensor;
  bool pinned = false;
  bool in_use = false;
};

static std::mutex pinned_pool_mutex;
static std::vector<staging_slot> *pinned_pool = nullptr;
static int64_t pinned_bytes = 0;
static thread_local int pinned_current_slot = -1;
static std::atomic<uint64_t> pinned_generation{0};
static thread_local uint64_t pinned_generation_local = 0;

static void sync_pinned_generation_locked() {
  const auto generation = pinned_generation.load(std::memory_order_relaxed);
  if (generation != pinned_generation_local) {
    pinned_current_slot = -1;
    pinned_generation_local = generation;
  }
}

static std::vector<staging_slot> &get_pinned_pool() {
  if (pinned_pool == nullptr) {
    pinned_pool = new std::vector<staging_slot>();
  }
  return *pinned_pool;
}

static void release_current_slot_locked() {
  if (pinned_current_slot >= 0 && pinned_pool != nullptr &&
      pinned_current_slot < static_cast<int>(pinned_pool->size())) {
    (*pinned_pool)[pinned_current_slot].in_use = false;
  }
  pinned_current_slot = -1;
}

static torch::Tensor make_host_buffer(const int64_t n_doubles, bool &pinned,
                                      const bool allow_pinned) {
  const auto opts = torch::TensorOptions().dtype(torch::kFloat64);
  if (allow_pinned && torch::cuda::is_available()) {
    try {
      auto buffer = torch::empty({n_doubles}, opts.pinned_memory(true));
      pinned = true;
      return buffer;
    } catch (const c10::Error &) {
    }
  }
  pinned = false;
  return torch::empty({n_doubles}, opts);
}

void torch_c_pinned_reserve(const int64_t n_doubles) {
  TORCH_C_TRY
  std::lock_guard<std::mutex> lock(pinned_pool_mutex);
  sync_pinned_generation_locked();
  auto &slots = get_pinned_pool();
  release_current_slot_locked();
  int chosen = -1;
  for (int i = 0; i < static_cast<int>(slots.size()); ++i) {
    if (!slots[i].in_use && slots[i].tensor.defined() &&
        slots[i].tensor.numel() >= n_doubles) {
      chosen = i;
      break;
    }
  }
  if (chosen < 0) {
    for (int i = 0; i < static_cast<int>(slots.size()); ++i) {
      if (!slots[i].in_use) {
        chosen = i;
        break;
      }
    }
  }
  if (chosen < 0) {
    slots.emplace_back();
    chosen = static_cast<int>(slots.size()) - 1;
  }
  if (!slots[chosen].tensor.defined() ||
      slots[chosen].tensor.numel() < n_doubles) {
    // Grow to at least the requested size (never shrink) to avoid re-pinning
    // churn when successive uploads vary in size.
    const int64_t old_numel =
        slots[chosen].tensor.defined() ? slots[chosen].tensor.numel() : 0;
    const int64_t target = std::max<int64_t>(n_doubles, old_numel);
    const int64_t old_bytes =
        (slots[chosen].tensor.defined() && slots[chosen].pinned)
            ? old_numel * static_cast<int64_t>(sizeof(double))
            : 0;
    const int64_t want_bytes = target * static_cast<int64_t>(sizeof(double));
    const bool allow_pinned =
        pinned_bytes - old_bytes + want_bytes <= torch_pinned_budget_bytes();
    bool pinned = false;
    slots[chosen].tensor = make_host_buffer(target, pinned, allow_pinned);
    slots[chosen].pinned = pinned;
    pinned_bytes = pinned_bytes - old_bytes + (pinned ? want_bytes : 0);
#if defined(__LIBTORCH_CUDA)
    if (!pinned) {
      warn_pageable_pinned_once("Torch staging upload");
    }
#endif
  }
  slots[chosen].in_use = true;
  pinned_current_slot = chosen;
  TORCH_C_CATCH("torch_c_pinned_reserve")
}

/*******************************************************************************
 * \brief Returns the current thread's staging buffer to the pool.
 ******************************************************************************/
void torch_c_pinned_finish() {
  TORCH_C_TRY
  std::lock_guard<std::mutex> lock(pinned_pool_mutex);
  sync_pinned_generation_locked();
  release_current_slot_locked();
  TORCH_C_CATCH("torch_c_pinned_finish")
}

/*******************************************************************************
 * \brief Releases all staging buffers (e.g. at environment teardown).
 ******************************************************************************/
void torch_c_pinned_release() {
  TORCH_C_TRY {
    synchronize_pending_grad_copies();
    std::lock_guard<std::mutex> lock(pinned_pool_mutex);
    sync_pinned_generation_locked();
    release_current_slot_locked();
    int in_use = 0;
    if (pinned_pool != nullptr) {
      for (const auto &slot : *pinned_pool) {
        in_use += slot.in_use ? 1 : 0;
      }
    }
    if (in_use == 0) {
      delete pinned_pool;
      pinned_pool = nullptr;
      pinned_bytes = 0;
      pinned_generation.fetch_add(1, std::memory_order_relaxed);
    }
  }
  TORCH_C_CATCH("torch_c_pinned_release")
}

/*******************************************************************************
 * \brief Uploads one double array through this thread's staging buffer.
 *        The host copy into the staging buffer is synchronous; the device copy
 *        is asynchronous (when pinned) and must be followed by
 *        torch_c_stream_synchronize() and then torch_c_pinned_finish().
 ******************************************************************************/
bool torch_c_tensor_upload_pinned(torch_c_tensor_t **tensor,
                                  const int64_t offset, const bool req_grad,
                                  const int ndims, const int64_t sizes[],
                                  double source[]) {
  TORCH_C_TRY
  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  const auto sizes_ref = c10::IntArrayRef(sizes, ndims);
  int64_t count = 1;
  for (int i = 0; i < ndims; i++) {
    torch_c_check(sizes[i] >= 0,
                  "Torch tensor upload has a negative dimension");
    if (sizes[i] == 0) {
      count = 0;
      continue;
    }
    torch_c_check(count <= std::numeric_limits<int64_t>::max() / sizes[i],
                  "Torch tensor upload element count overflows");
    count *= sizes[i];
  }
  const bool reallocated =
      prepare_target_tensor(tensor, req_grad, device, ndims, sizes);
  const auto source_tensor = torch::from_blob(
      source, sizes_ref, torch::TensorOptions().dtype(torch::kFloat64));
  if (device.is_cuda()) {
    torch::Tensor staging;
    bool staged_pinned = false;
    {
      std::lock_guard<std::mutex> lock(pinned_pool_mutex);
      sync_pinned_generation_locked();
      auto &slots = get_pinned_pool();
      torch_c_check(pinned_current_slot >= 0 &&
                        pinned_current_slot < static_cast<int>(slots.size()),
                    "Torch staging buffer was not reserved for this thread");
      auto &slot = slots[pinned_current_slot];
      staged_pinned = slot.pinned;
      if (staged_pinned) {
        torch_c_check(slot.tensor.defined(),
                      "Torch pinned staging buffer is undefined");
        const int64_t capacity = slot.tensor.numel();
        torch_c_check(offset >= 0 && offset <= capacity &&
                          count <= capacity - offset,
                      "Torch pinned staging buffer is too small");
        staging = slot.tensor;
      }
    }
    torch::NoGradGuard no_grad;
    if (staged_pinned) {
      auto staging_slice = staging.narrow(0, offset, count).view(sizes_ref);
      staging_slice.copy_(source_tensor);
      (*tensor)->copy_(staging_slice, /*non_blocking=*/true);
    } else {
      // Pageable staging would only add a redundant host-to-host copy; upload
      // straight from the source buffer (synchronous host-to-device).
      (*tensor)->copy_(source_tensor);
    }
    (*tensor)->mutable_grad() = torch::Tensor();
  } else {
    torch::NoGradGuard no_grad;
    (*tensor)->copy_(source_tensor);
    (*tensor)->mutable_grad() = torch::Tensor();
  }
  (*tensor)->set_requires_grad(req_grad);
  return reallocated;
  TORCH_C_CATCH("torch_c_tensor_upload_pinned")
}

/*******************************************************************************
 * \brief Waits for all pending device copies issued by the upload/readback
 * paths.
 ******************************************************************************/
void torch_c_stream_synchronize() {
  TORCH_C_TRY
  synchronize_pending_grad_copies();
  TORCH_C_CATCH("torch_c_stream_synchronize")
}

/*******************************************************************************
 * \brief Creates an expanded tensor view along one singleton dimension.
 ******************************************************************************/
void torch_c_tensor_expand_dim(const torch_c_tensor_t *tensor,
                               const int64_t dim, const int64_t size,
                               torch_c_tensor_t **result) {
  TORCH_C_TRY
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  torch_c_check(*result == NULL, "Torch tensor output must be NULL");
  torch_c_check(dim >= 0, "Torch tensor dimension must be non-negative");
  torch_c_check(dim < tensor->dim(), "Torch tensor dimension is out of range");
  std::vector<int64_t> sizes(tensor->sizes().begin(), tensor->sizes().end());
  torch_c_check(sizes[dim] == 1,
                "Torch tensor expansion requires a singleton dimension");
  torch_c_check(size >= 0, "Torch tensor expansion size must be non-negative");
  sizes[dim] = size;
  *result = new torch_c_tensor_t(tensor->expand(sizes));
  TORCH_C_CATCH("torch_c_tensor_expand_dim")
}

/*******************************************************************************
 * \brief Creates a tensor view narrowed along one dimension.
 *        The result shares storage with the input tensor and stays on its
 *        device; the input must already be on the active device.
 ******************************************************************************/
void torch_c_tensor_narrow(const torch_c_tensor_t *tensor, const int64_t dim,
                           const int64_t start_index, const int64_t length,
                           torch_c_tensor_t **result) {
  TORCH_C_TRY
  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  torch_c_check(*result == NULL, "Torch tensor output must be NULL");
  torch_c_check(dim >= 0, "Torch tensor dimension must be non-negative");
  torch_c_check(start_index >= 0,
                "Torch tensor start index must be non-negative");
  torch_c_check(length >= 0, "Torch tensor length must be non-negative");
  torch_c_check(dim < tensor->ndimension(),
                "Torch tensor dimension is out of range");
  torch_c_check(start_index <= tensor->size(dim) &&
                    length <= tensor->size(dim) - start_index,
                "Torch tensor narrow range is out of bounds");
  torch_c_check(
      tensor->device() == device,
      "Torch tensor narrow requires the tensor to be on the active device");
  *result = new torch_c_tensor_t(tensor->narrow(dim, start_index, length));
  TORCH_C_CATCH("torch_c_tensor_narrow")
}

/*******************************************************************************
 * \brief Deep-copies a tensor into an independent tensor on the current device.
 ******************************************************************************/
void torch_c_tensor_clone(const torch_c_tensor_t *tensor,
                          torch_c_tensor_t **result) {
  TORCH_C_TRY
  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  torch_c_check(*result == NULL, "Torch tensor output must be NULL");
  *result = new torch_c_tensor_t(tensor->detach().to(device).clone());
  TORCH_C_CATCH("torch_c_tensor_clone")
}

/*******************************************************************************
 * \brief Returns the data_ptr and sizes of a Torch tensor of int32s.
 *        The returned pointer is only valide during the tensor's live time!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_data_ptr_int32(const torch_c_tensor_t *tensor,
                                   const int ndims, int64_t sizes[],
                                   int32_t **data_ptr) {
  TORCH_C_TRY
  *data_ptr = (int32_t *)get_data_ptr(tensor, torch::kInt32, ndims, sizes);
  TORCH_C_CATCH("torch_c_tensor_data_ptr_int32")
}

/*******************************************************************************
 * \brief Returns the data_ptr and sizes of a Torch tensor of floats.
 *        The returned pointer is only valide during the tensor's lifetime!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_data_ptr_float(const torch_c_tensor_t *tensor,
                                   const int ndims, int64_t sizes[],
                                   float **data_ptr) {
  TORCH_C_TRY
  *data_ptr = (float *)get_data_ptr(tensor, torch::kFloat32, ndims, sizes);
  TORCH_C_CATCH("torch_c_tensor_data_ptr_float")
}

/*******************************************************************************
 * \brief Returns the data_ptr and sizes of a Torch tensor of int64s.
 *        The returned pointer is only valide during the tensor's live time!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_data_ptr_int64(const torch_c_tensor_t *tensor,
                                   const int ndims, int64_t sizes[],
                                   int64_t **data_ptr) {
  TORCH_C_TRY
  *data_ptr = (int64_t *)get_data_ptr(tensor, torch::kInt64, ndims, sizes);
  TORCH_C_CATCH("torch_c_tensor_data_ptr_int64")
}

/*******************************************************************************
 * \brief Returns the data_ptr and sizes of a Torch tensor of doubles.
 *        The returned pointer is only valide during the tensor's live time!
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_data_ptr_double(const torch_c_tensor_t *tensor,
                                    const int ndims, int64_t sizes[],
                                    double **data_ptr) {
  TORCH_C_TRY
  *data_ptr = (double *)get_data_ptr(tensor, torch::kFloat64, ndims, sizes);
  TORCH_C_CATCH("torch_c_tensor_data_ptr_double")
}

/*******************************************************************************
 * \brief Runs autograd on a Torch tensor.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_backward(const torch_c_tensor_t *tensor,
                             const torch_c_tensor_t *outer_grad) {
  TORCH_C_TRY
  TorchFloatingPointMaskGuard fpe_guard;
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  tensor->backward(*outer_grad);
  TORCH_C_CATCH("torch_c_tensor_backward")
}

/*******************************************************************************
 * \brief Runs autograd on a scalar Torch tensor.
 ******************************************************************************/
void torch_c_tensor_backward_scalar(const torch_c_tensor_t *tensor) {
  TORCH_C_TRY
  TorchFloatingPointMaskGuard fpe_guard;
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  tensor->backward();
  TORCH_C_CATCH("torch_c_tensor_backward_scalar")
}

/*******************************************************************************
 * \brief Moves a tensor to the active device and makes it an autograd leaf.
 ******************************************************************************/
void torch_c_tensor_to_device_leaf(torch_c_tensor_t **tensor,
                                   const bool req_grad) {
  TORCH_C_TRY
  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  auto moved = (*tensor)->to(device).detach();
  moved.set_requires_grad(req_grad);
  erase_grad_pool_entry(*tensor);
  delete (*tensor);
  *tensor = new torch_c_tensor_t(moved);
  TORCH_C_CATCH("torch_c_tensor_to_device_leaf")
}

/*******************************************************************************
 * \brief Sets a leaf tensor's requires_grad flag and clears any stale grad.
 ******************************************************************************/
void torch_c_tensor_set_requires_grad(torch_c_tensor_t *tensor,
                                      const bool requires_grad) {
  TORCH_C_TRY
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  torch_c_check(tensor != nullptr && tensor->defined(),
                "Torch tensor set_requires_grad received an undefined tensor");
  tensor->set_requires_grad(requires_grad);
  tensor->mutable_grad() = torch::Tensor();
  TORCH_C_CATCH("torch_c_tensor_set_requires_grad")
}

/*******************************************************************************
 * \brief Select whether Torch wrappers should use CUDA when available.
 ******************************************************************************/
void torch_c_use_cuda(const bool use_cuda) {
  use_cuda_if_available.store(use_cuda);
}

/*******************************************************************************
 * \brief Materializes a contiguous host-owned copy of an autograd gradient.
 ******************************************************************************/
static torch::Tensor materialize_host_gradient(const torch::Tensor &source,
                                               const torch::Tensor &like) {
  if (!source.defined()) {
    return torch::zeros(
        like.sizes(), like.options().device(torch::kCPU).requires_grad(false));
  }
  auto host =
      torch::empty(source.sizes(),
                   source.options().device(torch::kCPU).requires_grad(false));
  torch::NoGradGuard no_grad;
  host.copy_(source);
  return host;
}

static void grad_batch_impl(const torch_c_tensor_t *const *tensors,
                            const int count, torch_c_tensor_t **grads,
                            const bool synchronize) {
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  torch_c_check(count >= 0, "Torch gradient batch size must be non-negative");
  if (count == 0)
    return;
  for (int i = 0; i < count; ++i) {
    torch_c_check(tensors[i] != nullptr,
                  "Torch gradient input must not be NULL");
    torch_c_check(grads[i] == nullptr, "Torch gradient outputs must be NULL");
  }

  std::vector<torch::Tensor> sources(count);
  bool all_cuda_same_device = true;
  int source_device = -1;
  for (int i = 0; i < count; ++i) {
    sources[i] = tensors[i]->grad();
    if (!sources[i].defined()) {
      all_cuda_same_device = false;
      continue;
    }
    if (!sources[i].is_cuda()) {
      all_cuda_same_device = false;
      continue;
    }
    const int index = static_cast<int>(sources[i].get_device());
    if (source_device < 0) {
      source_device = index;
    } else if (source_device != index) {
      all_cuda_same_device = false;
    }
  }

  if (all_cuda_same_device) {
    c10::DeviceGuard source_guard(torch::Device(
        torch::kCUDA, static_cast<c10::DeviceIndex>(source_device)));
    std::vector<torch::Tensor> hosts(count);
    {
      std::lock_guard<std::mutex> lock(grad_pinned_mutex);
      if (grad_pinned_pool == nullptr) {
        grad_pinned_pool = new std::unordered_map<const torch_c_tensor_t *,
                                                  grad_pinned_entry>();
      }
      ++grad_pinned_clock;
      reclaim_retired_pinned_locked();
      for (int i = 0; i < count; ++i) {
        const auto sizes = sources[i].sizes();
        auto it = grad_pinned_pool->find(tensors[i]);
        bool reuse = false;
        if (it != grad_pinned_pool->end()) {
          const auto &entry = it->second;
          reuse = entry.tensor.defined() && entry.tensor.sizes() == sizes &&
                  entry.tensor.scalar_type() == sources[i].scalar_type() &&
                  entry.tensor.getIntrusivePtr().use_count() == 1 &&
                  !pending_grad_copy(tensors[i]);
        }
        if (reuse) {
          it->second.last_use = grad_pinned_clock;
          hosts[i] = it->second.tensor;
          continue;
        }
        if (it != grad_pinned_pool->end()) {
          if (pending_grad_copy(tensors[i])) {
            // Deferred readback into the old pinned buffer is still in flight.
            // Wait before dropping the last reference to it.
            synchronize_pending_grad_copies();
          }
          drop_grad_pool_entry_locked(it);
        }
        const int64_t need =
            sources[i].numel() * static_cast<int64_t>(sizeof(double));
        evict_grad_pool_locked(need);
        const bool allow_pinned =
            grad_pinned_bytes + need <= torch_pinned_budget_bytes() &&
            torch::cuda::is_available();
        bool pinned = false;
        torch::Tensor host;
        if (allow_pinned) {
          try {
            host = torch::empty(
                sizes,
                sources[i].options().device(torch::kCPU).pinned_memory(true));
            pinned = true;
          } catch (const c10::Error &) {
            pinned = false;
          }
        }
        if (!host.defined()) {
          host = torch::empty(sizes, sources[i].options().device(torch::kCPU));
        }
        if (!pinned) {
          warn_pageable_pinned_once("Torch gradient readback");
        }
        auto &entry = (*grad_pinned_pool)[tensors[i]];
        entry.tensor = host;
        entry.pinned = pinned;
        entry.bytes = pinned ? need : 0;
        entry.last_use = grad_pinned_clock;
        if (pinned) {
          grad_pinned_bytes += need;
        }
        hosts[i] = host;
      }
    }
    for (int i = 0; i < count; ++i) {
      hosts[i].copy_(sources[i], true);
    }
    if (synchronize) {
#if defined(__LIBTORCH_CUDA)
      c10::cuda::getCurrentCUDAStream(
          static_cast<c10::DeviceIndex>(source_device))
          .synchronize();
#else
      (void)source_device;
#endif
    } else {
      record_pending_grads(tensors, count, source_device);
    }
    for (int i = 0; i < count; ++i) {
      grads[i] = new torch_c_tensor_t(std::move(hosts[i]));
    }
  } else {
    for (int i = 0; i < count; ++i) {
      grads[i] = new torch_c_tensor_t(
          materialize_host_gradient(sources[i], *tensors[i]));
    }
  }
}

void torch_c_tensor_grad_batch(const torch_c_tensor_t *const *tensors,
                               const int count, torch_c_tensor_t **grads) {
  TORCH_C_TRY
  grad_batch_impl(tensors, count, grads, /*synchronize=*/true);
  TORCH_C_CATCH("torch_c_tensor_grad_batch")
}

void torch_c_tensor_grad_batch_deferred(const torch_c_tensor_t *const *tensors,
                                        const int count,
                                        torch_c_tensor_t **grads) {
  TORCH_C_TRY
  grad_batch_impl(tensors, count, grads, /*synchronize=*/false);
  TORCH_C_CATCH("torch_c_tensor_grad_batch_deferred")
}

/*******************************************************************************
 * \brief Releases the cached pinned gradient buffers (e.g. at environment
 * teardown).
 ******************************************************************************/
void torch_c_grad_pool_clear() {
  TORCH_C_TRY {
    synchronize_pending_grad_copies();
    std::lock_guard<std::mutex> lock(grad_pinned_mutex);
    delete grad_pinned_pool;
    grad_pinned_pool = nullptr;
    grad_pinned_retired.clear();
    grad_pinned_bytes = 0;
  }
  TORCH_C_CATCH("torch_c_grad_pool_clear")
}

/*******************************************************************************
 * \brief Releases a Torch tensor and all its ressources.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_tensor_release(torch_c_tensor_t *tensor) {
  if (tensor == nullptr)
    return;
  erase_grad_pool_entry(tensor);
  delete (tensor);
}

/*******************************************************************************
 * \brief Creates an empty Torch dictionary.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_dict_create(torch_c_dict_t **dict_out) {
  TORCH_C_TRY
  torch_c_check(*dict_out == NULL, "Torch dictionary output must be NULL");
  *dict_out = new c10::Dict<std::string, torch::Tensor>();
  TORCH_C_CATCH("torch_c_dict_create")
}

/*******************************************************************************
 * \brief Clones a Torch dictionary.
 ******************************************************************************/
void torch_c_dict_clone(const torch_c_dict_t *dict, torch_c_dict_t **dict_out) {
  TORCH_C_TRY
  torch_c_check(*dict_out == NULL, "Torch dictionary output must be NULL");
  torch_c_dict_t *clone = new c10::Dict<std::string, torch::Tensor>();
  for (const auto &entry : *dict) {
    clone->insert(entry.key(), entry.value());
  }
  *dict_out = clone;
  TORCH_C_CATCH("torch_c_dict_clone")
}

/*******************************************************************************
 * \brief Inserts a Torch tensor into a Torch dictionary.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_dict_insert(const torch_c_dict_t *dict, const char *key,
                         const torch_c_tensor_t *tensor) {
  TORCH_C_TRY
  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  if (tensor->device() == device) {
    dict->insert(key, *tensor);
  } else {
    dict->insert(key, tensor->to(device));
  }
  TORCH_C_CATCH("torch_c_dict_insert")
}

/*******************************************************************************
 * \brief Inserts or replaces a Torch tensor in a Torch dictionary.
 ******************************************************************************/
void torch_c_dict_set(const torch_c_dict_t *dict, const char *key,
                      const torch_c_tensor_t *tensor) {
  TORCH_C_TRY
  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  auto *mutable_dict = const_cast<torch_c_dict_t *>(dict);
  if (tensor->device() == device) {
    mutable_dict->insert_or_assign(key, *tensor);
  } else {
    mutable_dict->insert_or_assign(key, tensor->to(device));
  }
  TORCH_C_CATCH("torch_c_dict_set")
}

/*******************************************************************************
 * \brief Retrieves a Torch tensor from a Torch dictionary.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_dict_get(const torch_c_dict_t *dict, const char *key,
                      torch_c_tensor_t **tensor) {
  TORCH_C_TRY
  torch_c_check(dict->contains(key), "Torch dictionary key not found");
  *tensor = new torch_c_tensor_t(dict->at(key).cpu().contiguous());
  TORCH_C_CATCH("torch_c_dict_get")
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
  TORCH_C_TRY
  torch_c_check(*model_out == NULL, "Torch model output must be NULL");
  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  set_jit_fusion_strategy();
  torch_c_model_t *model = new torch_c_model_t();
  model->module = load_module_for_device(filename, device);
  model->module.eval(); // Set inference behavior for modules such as dropout.
  *model_out = model;
  TORCH_C_CATCH("torch_c_model_load")
}

/*******************************************************************************
 * \brief Loads a Torch model and reads two metadata entries.
 ******************************************************************************/
void torch_c_model_load_with_metadata(torch_c_model_t **model_out,
                                      const char *filename, const char *key1,
                                      const char *key2, char **content1,
                                      int *length1, char **content2,
                                      int *length2) {
  TORCH_C_TRY
  torch_c_check(*model_out == NULL, "Torch model output must be NULL");
  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  std::unordered_map<std::string, std::string> extra_files = {{key1, ""},
                                                              {key2, ""}};
  set_jit_fusion_strategy();
  torch_c_model_t *model = new torch_c_model_t();
  // Deduplicate Skala constants on CPU without repeated CUDA synchronization.
  model->module =
      load_module_for_device(filename, device, &extra_files, device.is_cuda());
  model->module.eval(); // Set inference behavior for modules such as dropout.
  *model_out = model;
  copy_string_to_c_buffer(extra_files[key1], content1, length1);
  copy_string_to_c_buffer(extra_files[key2], content2, length2);
  TORCH_C_CATCH("torch_c_model_load_with_metadata")
}

/*******************************************************************************
 * \brief Maps serialized TorchScript device constants to the active device.
 ******************************************************************************/
void torch_c_model_remap_device_constants(torch_c_model_t *model) {
  TORCH_C_TRY
  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  remap_model_device_constants(model->module, device);
  if (device.is_cuda()) {
    initialize_cuda_linalg(device);
  }

  bool has_cuda = false;
  bool has_cpu = false;
  collect_model_device_constants(model->module, has_cuda, has_cpu);
  if (device.is_cuda()) {
    torch_c_check(!has_cpu,
                  "SKALA TorchScript model still targets the CPU "
                  "after remapping; use the CUDA export for GPU runs");
  } else {
    torch_c_check(!has_cuda, "SKALA TorchScript model still targets CUDA after "
                             "remapping; use the CPU export or a CUDA build");
  }
  TORCH_C_CATCH("torch_c_model_remap_device_constants")
}

/*******************************************************************************
 * \brief Promote float32 model state and explicit numerical casts to float64.
 *        Trained values are retained exactly; integer/index operations are
 * kept.
 ******************************************************************************/
void torch_c_model_promote_float32(torch_c_model_t *model) {
  TORCH_C_TRY
  TorchFloatingPointMaskGuard fpe_guard;
  torch::NoGradGuard no_grad;
  promote_model_float32(model->module);
  TORCH_C_CATCH("torch_c_model_promote_float32")
}

/*******************************************************************************
 * \brief Prepare a model for higher input derivatives.
 ******************************************************************************/
void torch_c_model_prepare_higher_derivatives(torch_c_model_t *model) {
  TORCH_C_TRY
  // Preserve parameter values while avoiding fp32 cancellation and saturated
  // SiLU backward in higher input derivatives.
  torch_c_model_promote_float32(model);
  torch::jit::SubgraphRewriter rewriter;
  rewriter.RegisterRewritePattern(
      "graph(%x):\n %y = aten::silu(%x)\n return (%y)", R"IR(
graph(%x):
  %s = aten::sigmoid(%x)
  %y = aten::mul(%x, %s)
  return (%y))IR");
  std::unordered_set<torch::jit::Graph *> prepared;
  for (const auto &module : model->module.modules()) {
    for (const auto &method : module.get_methods()) {
      auto graph = method.graph();
      if (prepared.insert(graph.get()).second) {
        rewriter.runOnGraph(graph);
      }
    }
  }
  TORCH_C_CATCH("torch_c_model_prepare_higher_derivatives")
}

/*******************************************************************************
 * \brief Applies the Hessian of a scalar tensor to the given direction tensors.
 ******************************************************************************/
void torch_c_tensor_hessian_vector(const torch_c_tensor_t *tensor,
                                   const int count,
                                   const torch_c_tensor_t *const inputs[],
                                   const torch_c_tensor_t *const directions[],
                                   torch_c_tensor_t *responses[],
                                   const bool to_host, int *status) {
  *status = 0;
  for (int i = 0; i < count; i++)
    responses[i] = nullptr;
  try {
    TorchFloatingPointMaskGuard fpe_guard;
    c10::OptionalDeviceGuard guard;
    get_device_with_guard(guard);
    TORCH_CHECK(
        count > 0 && tensor->numel() == 1,
        "Hessian-vector evaluation requires a scalar and input tensors");
    std::vector<torch::Tensor> variables;
    for (int i = 0; i < count; i++) {
      TORCH_CHECK(inputs[i]->requires_grad() &&
                      inputs[i]->sizes() == directions[i]->sizes() &&
                      inputs[i]->device() == directions[i]->device(),
                  "Hessian-vector inputs and directions must have matching "
                  "shapes and devices");
      variables.push_back(*inputs[i]);
    }
    // Accumulate the Hessian-vector product one input at a time.  This keeps at
    // most one first-order (create_graph) result alive instead of holding the
    // first-order gradients for every input simultaneously, lowering peak
    // memory for the full-grid response at the cost of one backward per input.
    // A zero direction contributes nothing, so every input is included without
    // a per-input device synchronization to detect all-zero directions.
    std::vector<torch::Tensor> second(count);
    if (tensor->requires_grad()) {
      for (int i = 0; i < count; ++i) {
        std::vector<torch::Tensor> first = torch::autograd::grad(
            {*tensor}, {variables[i]}, {}, true, true, true);
        if (first.empty() || !first[0].defined() || !first[0].requires_grad())
          continue;
        auto term = (first[0] * directions[i]->detach()).sum();
        std::vector<torch::Tensor> partial =
            torch::autograd::grad({term}, variables, {}, true, false, true);
        term = torch::Tensor();
        first.clear();
        for (int j = 0; j < count; j++) {
          if (partial[j].defined())
            second[j] =
                second[j].defined() ? second[j] + partial[j] : partial[j];
        }
      }
    }
    for (int i = 0; i < count; i++) {
      auto response = second[i].defined() ? second[i].detach()
                                          : torch::zeros_like(*inputs[i]);
      responses[i] = new torch_c_tensor_t(to_host ? response.cpu().contiguous()
                                                  : response);
    }
  } catch (const std::exception &error) {
    fprintf(stderr, "Torch Hessian evaluation failed: %s\n", error.what());
    for (int i = 0; i < count; i++) {
      delete responses[i];
      responses[i] = nullptr;
    }
    *status = 1;
  }
}

/*******************************************************************************
 * \brief Disables gradients for inference-only model parameters.
 ******************************************************************************/
void torch_c_model_disable_parameter_gradients(torch_c_model_t *model) {
  TORCH_C_TRY
  torch::NoGradGuard no_grad;
  for (auto parameter : model->module.parameters()) {
    parameter.set_requires_grad(false);
  }
  TORCH_C_CATCH("torch_c_model_disable_parameter_gradients")
}

/*******************************************************************************
 * \brief Evaluates the given Torch model.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_forward(torch_c_model_t *model, const torch_c_dict_t *inputs,
                           torch_c_dict_t *outputs) {
  TORCH_C_TRY

  TorchFloatingPointMaskGuard fpe_guard;
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  auto untyped_output = model->module.forward({*inputs}).toGenericDict();
  outputs->clear();
  for (const auto &entry : untyped_output) {
    outputs->insert(entry.key().toStringView(), entry.value().toTensor());
  }
  TORCH_C_CATCH("torch_c_model_forward")
}

/*******************************************************************************
 * \brief Evaluates a TorchScript model method expecting argument "mol".
 ******************************************************************************/
// Returns the cached method handle for a frozen model. The primary method
// (resolved at freeze time) is served lock-free on the hot path; other names
// fall back to a per-model map guarded by the model's own mutex.
static torch::jit::Method get_model_method(torch_c_model_t *model,
                                           const char *method_name) {
  if (model->has_primary_method && model->primary_method_name == method_name) {
    return *model->primary_method;
  }
  std::lock_guard<std::mutex> lock(model->method_mutex);
  const std::string key(method_name);
  auto it = model->extra_methods.find(key);
  if (it == model->extra_methods.end()) {
    it = model->extra_methods.emplace(key, model->module.get_method(key)).first;
  }
  return it->second;
}

void torch_c_model_forward_mol_tensor(torch_c_model_t *model,
                                      const char *method_name,
                                      const torch_c_dict_t *inputs,
                                      torch_c_tensor_t **output) {
  TORCH_C_TRY

  c10::OptionalDeviceGuard guard;
  const auto device = get_device_with_guard(guard);
  // Avoid repeated TorchScript profiling of the frozen Skala model on CUDA.
  std::optional<torch::jit::GraphOptimizerEnabledGuard> graph_optimization;
  if (device.is_cuda()) {
    graph_optimization.emplace(false);
  }
  torch_c_check(*output == NULL, "Torch tensor output must be NULL");
  auto method = get_model_method(model, method_name);
  *output = new torch_c_tensor_t(method({*inputs}).toTensor());
  TORCH_C_CATCH("torch_c_model_forward_mol_tensor")
}

static void flatten_weighted_pair(const torch_c_tensor_t *values,
                                  const torch_c_tensor_t *weights,
                                  torch::Tensor &flat_values,
                                  torch::Tensor &flat_weights) {
  flat_values = values->flatten();
  flat_weights = weights->to(values->device()).to(values->dtype()).flatten();
  torch_c_check(flat_values.numel() == flat_weights.numel(),
                "Torch weighted sum: value/weight size mismatch");
}

/*******************************************************************************
 * \brief Returns the weighted sum of two Torch tensors.
 ******************************************************************************/
void torch_c_tensor_weighted_sum(const torch_c_tensor_t *values,
                                 const torch_c_tensor_t *weights,
                                 torch_c_tensor_t **result) {
  TORCH_C_TRY
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  torch::Tensor flat_values, flat_weights;
  flatten_weighted_pair(values, weights, flat_values, flat_weights);
  *result = new torch_c_tensor_t(torch::dot(flat_values, flat_weights));
  TORCH_C_CATCH("torch_c_tensor_weighted_sum")
}

/*******************************************************************************
 * \brief Compute per-segment weighted sums and also return the total weighted
 *        sum as a differentiable scalar tensor (the sum of the segment dots).
 *        Avoids a second full-length reduction and weight conversion.
 ******************************************************************************/
void torch_c_tensor_weighted_sum_segments_total(const torch_c_tensor_t *values,
                                                const torch_c_tensor_t *weights,
                                                const int64_t *offsets,
                                                const int64_t *lengths,
                                                const int nseg, double *result,
                                                torch_c_tensor_t **total) {
  TORCH_C_TRY
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  torch_c_check(nseg >= 0,
                "Torch weighted sum segment count must be non-negative");
  torch_c_check(values != nullptr && weights != nullptr && offsets != nullptr &&
                    lengths != nullptr && result != nullptr && total != nullptr,
                "Torch weighted sum received a NULL argument");
  torch::Tensor flat_values, flat_weights;
  flatten_weighted_pair(values, weights, flat_values, flat_weights);
  if (nseg == 0) {
    *total = new torch_c_tensor_t(flat_values.sum() * 0.0);
    return;
  }
  const int64_t numel = flat_values.numel();
  std::vector<torch::Tensor> dots;
  dots.reserve(nseg);
  for (int i = 0; i < nseg; ++i) {
    torch_c_check(offsets[i] >= 0 && lengths[i] >= 0 && offsets[i] <= numel &&
                      lengths[i] <= numel - offsets[i],
                  "Torch weighted sum range is out of bounds");
    dots.push_back(torch::dot(flat_values.narrow(0, offsets[i], lengths[i]),
                              flat_weights.narrow(0, offsets[i], lengths[i])));
  }
  auto stacked = torch::stack(dots);
  *total = new torch_c_tensor_t(stacked.sum());
  const auto host = stacked.to(torch::kFloat64).to(torch::kCPU).contiguous();
  std::memcpy(result, host.data_ptr<double>(), nseg * sizeof(double));
  TORCH_C_CATCH("torch_c_tensor_weighted_sum_segments_total")
}

/*******************************************************************************
 * \brief Adds a device scalar to an accumulator without a host synchronization.
 ******************************************************************************/
void torch_c_tensor_accumulate(torch_c_tensor_t **accumulator,
                               const torch_c_tensor_t *value) {
  TORCH_C_TRY
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  torch::NoGradGuard no_grad;
  torch_c_check(value != nullptr && value->defined(),
                "Torch accumulator received an undefined value");
  if (*accumulator == nullptr) {
    *accumulator = new torch_c_tensor_t(value->detach().clone());
  } else {
    torch_c_check((*accumulator)->defined() &&
                      (*accumulator)->scalar_type() == value->scalar_type() &&
                      (*accumulator)->device() == value->device() &&
                      (*accumulator)->sizes().equals(value->sizes()),
                  "Torch accumulator value shape/dtype/device mismatch");
    (*accumulator)->add_(value->detach());
  }
  TORCH_C_CATCH("torch_c_tensor_accumulate")
}

/*******************************************************************************
 * \brief Returns a scalar double value from a Torch tensor.
 ******************************************************************************/
double torch_c_tensor_item_double(const torch_c_tensor_t *tensor) {
  TORCH_C_TRY
  c10::OptionalDeviceGuard guard;
  get_device_with_guard(guard);
  return tensor->item<double>();
  TORCH_C_CATCH("torch_c_tensor_item_double")
}

/*******************************************************************************
 * \brief Releases a Torch model and all its ressources.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_release(torch_c_model_t *model) { delete (model); }

/*******************************************************************************
 * \brief Writes a size+mtime fingerprint of a file into a caller buffer.
 ******************************************************************************/
static std::string file_stamp(const char *filename) {
  struct stat file_stat{};
  if (stat(filename, &file_stat) != 0) {
    return std::string();
  }
#if defined(__linux__)
  return std::to_string(file_stat.st_size) + ":" +
         std::to_string(file_stat.st_mtim.tv_sec) + ":" +
         std::to_string(file_stat.st_mtim.tv_nsec);
#else
  return std::to_string(file_stat.st_size) + ":" +
         std::to_string(file_stat.st_mtime);
#endif
}

void torch_c_file_stamp(const char *filename, char *buf, int *length) {
  TORCH_C_TRY
  const std::string stamp = file_stamp(filename);
  if (length != nullptr) {
    if (buf != nullptr && *length > 0) {
      // Report the number of bytes actually written so callers can detect a
      // truncation instead of assuming the full stamp was copied.
      std::memset(buf, ' ', *length);
      const int n = std::min<int>(static_cast<int>(stamp.size()), *length);
      std::memcpy(buf, stamp.data(), n);
      *length = n;
    } else {
      // No buffer requested: return the required size.
      *length = static_cast<int>(stamp.size());
    }
  }
  TORCH_C_CATCH("torch_c_file_stamp")
}

void torch_c_model_read_metadata(const char *filename, const char *key,
                                 char **content, int *length) {
  TORCH_C_TRY

  const std::string filename_key(filename);
  const std::string stamp = file_stamp(filename);
  {
    std::lock_guard<std::mutex> lock(model_metadata_mutex);
    auto file_it = model_metadata_cache.find(filename_key);
    if (file_it != model_metadata_cache.end() &&
        file_it->second.stamp == stamp) {
      auto key_it = file_it->second.files.find(key);
      if (key_it != file_it->second.files.end()) {
        copy_string_to_c_buffer(key_it->second, content, length);
        return;
      }
    }
  }

  std::unordered_map<std::string, std::string> extra_files = {{key, ""}};
  torch::jit::load(filename, torch::kCPU, extra_files);
  const std::string content_str = extra_files[key];
  {
    std::lock_guard<std::mutex> lock(model_metadata_mutex);
    auto &entry = model_metadata_cache[filename_key];
    if (entry.stamp != stamp) {
      // The file changed; drop metadata cached for the previous revision so it
      // cannot be served under the new stamp.
      entry.files.clear();
      entry.stamp = stamp;
    }
    entry.files[key] = content_str;
  }
  copy_string_to_c_buffer(content_str, content, length);
  TORCH_C_CATCH("torch_c_model_read_metadata")
}

/*******************************************************************************
 * \brief Loads several TorchScript extra files (metadata keys) in one archive
 *        pass and stores them in the cache. Keys are newline-separated.
 ******************************************************************************/
void torch_c_model_prefetch_metadata(const char *filename, const char *keys) {
  TORCH_C_TRY
  const std::string filename_key(filename);
  const std::string stamp = file_stamp(filename);
  std::unordered_map<std::string, std::string> extra_files;
  const std::string all(keys);
  std::size_t start = 0;
  while (start <= all.size()) {
    const std::size_t end = all.find('\n', start);
    const std::string key = all.substr(
        start, end == std::string::npos ? std::string::npos : end - start);
    if (!key.empty()) {
      extra_files[key] = "";
    }
    if (end == std::string::npos) {
      break;
    }
    start = end + 1;
  }
  if (extra_files.empty()) {
    return;
  }
  torch::jit::load(filename, torch::kCPU, extra_files);
  std::lock_guard<std::mutex> lock(model_metadata_mutex);
  auto &entry = model_metadata_cache[filename_key];
  if (entry.stamp != stamp) {
    entry.files.clear();
    entry.stamp = stamp;
  }
  for (const auto &kv : extra_files) {
    entry.files[kv.first] = kv.second;
  }
  TORCH_C_CATCH("torch_c_model_prefetch_metadata")
}

/*******************************************************************************
 * \brief Clears the cached TorchScript model metadata.
 ******************************************************************************/
void torch_c_model_metadata_clear() {
  TORCH_C_TRY {
    std::lock_guard<std::mutex> lock(model_metadata_mutex);
    model_metadata_cache.clear();
  }
  TORCH_C_CATCH("torch_c_model_metadata_clear")
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
  int device_count = 0;
  cuda_device_available(device_count);
  return device_count;
}

/*******************************************************************************
 * \brief Set whether to allow TF32.
 *        Needed due to changes in defaults from pytorch 1.7 to 1.11 to >=1.12
 *        See https://pytorch.org/docs/stable/notes/cuda.html
 * \author Gabriele Tocci
 ******************************************************************************/
void torch_c_allow_tf32(const bool allow_tf32) {
  TORCH_C_TRY

  at::globalContext().setAllowTF32CuBLAS(allow_tf32);
  at::globalContext().setAllowTF32CuDNN(allow_tf32);
  TORCH_C_CATCH("torch_c_allow_tf32")
}

/******************************************************************************
 * \brief Freeze the Torch model: generic optimization that speeds up model.
 *        See https://pytorch.org/docs/stable/generated/torch.jit.freeze.html
 * \author Gabriele Tocci
 ******************************************************************************/
void torch_c_model_freeze(torch_c_model_t *model) {
  TORCH_C_TRY

  {
    std::lock_guard<std::mutex> lock(model->method_mutex);
    model->has_primary_method = false;
    model->primary_method_name.clear();
    model->primary_method.reset();
    model->extra_methods.clear();
  }
  model->module = torch::jit::freeze(model->module);
  TORCH_C_CATCH("torch_c_model_freeze")
}

/******************************************************************************
 * \brief Freeze a Torch model while preserving one exported method.
 ******************************************************************************/
void torch_c_model_freeze_preserving_method(torch_c_model_t *model,
                                            const char *method_name) {
  TORCH_C_TRY

  const std::vector<std::string> preserved_methods = {method_name};
  {
    std::lock_guard<std::mutex> lock(model->method_mutex);
    model->has_primary_method = false;
    model->primary_method_name.clear();
    model->primary_method.reset();
    model->extra_methods.clear();
  }
  torch::jit::freeze_module_inplace(&model->module, preserved_methods);
  // Resolve the preserved method now so the forward path is lock-free.
  std::lock_guard<std::mutex> lock(model->method_mutex);
  model->primary_method = model->module.get_method(method_name);
  model->primary_method_name = method_name;
  model->has_primary_method = true;
  TORCH_C_CATCH("torch_c_model_freeze_preserving_method")
}

/*******************************************************************************
 * \brief Retrieves an int64 attribute. Must be called before model freeze.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_get_attr_int64(const torch_c_model_t *model, const char *key,
                                  int64_t *dest) {
  TORCH_C_TRY
  *dest = model->module.attr(key).toInt();
  TORCH_C_CATCH("torch_c_model_get_attr_int64")
}

/*******************************************************************************
 * \brief Retrieves a double attribute. Must be called before model freeze.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_get_attr_double(const torch_c_model_t *model,
                                   const char *key, double *dest) {
  TORCH_C_TRY
  *dest = model->module.attr(key).toDouble();
  TORCH_C_CATCH("torch_c_model_get_attr_double")
}

/*******************************************************************************
 * \brief Retrieves a string attribute. Must be called before model freeze.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_get_attr_string(const torch_c_model_t *model,
                                   const char *key, char *dest,
                                   const int dest_len) {
  TORCH_C_TRY
  copy_attr_string_to_fortran(model->module.attr(key).toStringRef(), dest,
                              dest_len);
  TORCH_C_CATCH("torch_c_model_get_attr_string")
}

/*******************************************************************************
 * \brief Retrieves a list attribute's size. Must be called before model freeze.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_get_attr_list_size(const torch_c_model_t *model,
                                      const char *key, int *size) {
  TORCH_C_TRY
  const auto list_size = model->module.attr(key).toList().size();
  torch_c_check(list_size <= static_cast<size_t>(INT_MAX),
                "Torch attribute list is too long");
  *size = static_cast<int>(list_size);
  TORCH_C_CATCH("torch_c_model_get_attr_list_size")
}

/*******************************************************************************
 * \brief Retrieves a single item from a string list attribute.
 * \author Ole Schuett
 ******************************************************************************/
void torch_c_model_get_attr_strlist(const torch_c_model_t *model,
                                    const char *key, const int index,
                                    char *dest, const int dest_len) {
  TORCH_C_TRY
  const auto list = model->module.attr(key).toList();
  torch_c_check(index >= 0 && static_cast<size_t>(index) < list.size(),
                "Torch attribute list index is out of bounds");
  copy_attr_string_to_fortran(list[index].toStringRef(), dest, dest_len);
  TORCH_C_CATCH("torch_c_model_get_attr_strlist")
}

#ifdef __cplusplus
}
#endif

#endif // defined(__LIBTORCH)

// EOF
