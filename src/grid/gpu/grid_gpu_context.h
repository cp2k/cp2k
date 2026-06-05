/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/*
 * Authors :
 - Dr Mathieu Taillefumier (ETH Zurich / CSCS)
 - Advanced Micro Devices, Inc.
*/

#ifndef GRID_GPU_CONTEXT_H
#define GRID_GPU_CONTEXT_H

#ifdef __OFFLOAD_HIP
#include <hip/hip_runtime_api.h>
#else
#include <cuda_runtime.h>
#endif
#include <array>
#include <vector>

extern "C" {
#include "../common/grid_basis_set.h"
#include "../common/grid_constants.h"
}

#include "../../offload/offload_library.h"
#include "../../offload/offload_runtime.h"

namespace rocm_backend {
// a little helper class in the same spirit than std::vector. it must exist
// somewhere. Maybe possible to get the same thing with std::vector and
// specific allocator.

inline size_t round_up16(size_t n) { return (n + 15) & ~size_t(15); }

class smem_parameters;
template <typename T> class gpu_vector {
  size_t allocated_size_{0};
  size_t current_size_{0};
  bool internal_allocation_{false};
  T *device_ptr_{nullptr};
  T *host_ptr_{nullptr};

  void move_from(gpu_vector &&other) noexcept {
    allocated_size_ = other.allocated_size_;
    current_size_ = other.current_size_;
    internal_allocation_ = other.internal_allocation_;
    device_ptr_ = other.device_ptr_;
    host_ptr_ = other.host_ptr_;

    other.allocated_size_ = 0;
    other.current_size_ = 0;
    other.internal_allocation_ = false;
    other.device_ptr_ = nullptr;
    other.host_ptr_ = nullptr;
  }

public:
  gpu_vector() = default;

  gpu_vector(const gpu_vector &) = delete;
  gpu_vector &operator=(const gpu_vector &) = delete;

  gpu_vector(gpu_vector &&other) noexcept { move_from(std::move(other)); }

  gpu_vector &operator=(gpu_vector &&other) noexcept {
    if (this != &other) {
      reset();
      move_from(std::move(other));
    }
    return *this;
  }

  // size is the number of elements not the memory size
  explicit gpu_vector(const size_t size__) {
    allocated_size_ = (size__ < 16) ? 16 : round_up16(size__);
    current_size_ = size__;
    internal_allocation_ = true;

#ifndef __OFFLOAD_UNIFIED_MEMORY
    offloadMalloc((void **)&device_ptr_, sizeof(T) * allocated_size_);
#else
    hipMallocManaged((void **)&device_ptr_, sizeof(T) * allocated_size_);
#endif
    assert(device_ptr_ != nullptr);
  }

  gpu_vector(const size_t size__, void *ptr__) {
    allocated_size_ = size__;
    current_size_ = size__;
    internal_allocation_ = false;
    device_ptr_ = static_cast<T *>(ptr__);
  }
  ~gpu_vector() { reset(); }

  inline size_t size() const { return current_size_; }

  inline void copy_to_gpu(const T *data__) {
    assert(device_ptr_ != nullptr);
    assert(data__ != nullptr);
    offloadMemcpyHtoD(device_ptr_, data__, sizeof(T) * current_size_);
  }

  inline void copy_to_gpu(const T *data__, offloadStream_t &stream__) {
    assert(device_ptr_ != nullptr);
    assert(data__ != nullptr);
    offloadMemcpyAsyncHtoD(device_ptr_, data__, sizeof(T) * current_size_,
                           stream__);
  }

  inline void copy_associated_host_to_gpu(offloadStream_t &stream__) {
    assert(device_ptr_ != nullptr);
    assert(host_ptr_ != nullptr);
    // If the second assert fails it means that the object was created without
    // host buffer. It should not happen in the current scenario

    offloadMemcpyAsyncHtoD(device_ptr_, host_ptr_, sizeof(T) * current_size_,
                           stream__);
  }

  inline void copy_from_gpu(T *data__, offloadStream_t &stream__) {
    assert(device_ptr_ != nullptr);
    assert(data__ != nullptr);
    offloadMemcpyAsyncDtoH(data__, device_ptr_, sizeof(T) * current_size_,
                           stream__);
  }

  inline void copy_gpu_to_associated_host(offloadStream_t &stream__) {
    assert(device_ptr_ != nullptr);
    assert(host_ptr_ != nullptr);
    // If the second assert fails it means that the object was created without
    // host buffer. It should not happen in the current scenario

    offloadMemcpyAsyncDtoH(host_ptr_, device_ptr_, sizeof(T) * current_size_,
                           stream__);
  }

  inline void zero(offloadStream_t &stream__) {
    assert(device_ptr_ != nullptr);
    // zero device grid buffers
    offloadMemsetAsync(device_ptr_, 0, sizeof(T) * current_size_, stream__);
  }

  inline void associate(void *host_ptr__, void *device_ptr__,
                        const size_t size__) {
    assert(host_ptr__ != nullptr);
    assert(device_ptr__ != nullptr);
    reset();
    internal_allocation_ = false;
    allocated_size_ = size__;
    current_size_ = size__;
    device_ptr_ = static_cast<T *>(device_ptr__);
    host_ptr_ = static_cast<T *>(host_ptr__);
  }

  inline void zero() {
    // zero device grid buffers
    offloadMemset(device_ptr_, 0, sizeof(T) * current_size_);
  }

  inline void copy_to_gpu(const std::vector<T> &data__) {
    assert(data__.size() == current_size_);
    // if it fails it means that the vector on the gpu does not have the right
    // size. two option then
    // - resize the gpu vector
    // - or the cpu vector and gpu vector are not representing the quantity.
    assert(device_ptr_ != nullptr);
    offloadMemcpyHtoD(device_ptr_, data__.data(), sizeof(T) * data__.size());
  }

  inline void resize(const size_t new_size__) {
    if (!internal_allocation_) {
      allocated_size_ = 0;
      device_ptr_ = nullptr;
      host_ptr_ = nullptr;
    }

    assert(new_size__ != 0);
    if (allocated_size_ < new_size__) {
      if (internal_allocation_ && device_ptr_ != nullptr)
        offloadFree(device_ptr_);
      device_ptr_ = nullptr;
      allocated_size_ = (new_size__ < 16) ? 16 : round_up16(new_size__);
      offloadMalloc((void **)&device_ptr_, sizeof(T) * allocated_size_);
      internal_allocation_ = true;
      host_ptr_ = nullptr;
    }
    assert(device_ptr_ != nullptr);
    current_size_ = new_size__;
  }

  // does not invalidate the pointer. The memory is still allocated
  inline void clear() { current_size_ = 0; }

  // reset the class and free memory
  inline void reset() {
    if (internal_allocation_ && (device_ptr_ != nullptr)) {
      offloadFree(device_ptr_);
    }

    allocated_size_ = 0;
    current_size_ = 0;
    device_ptr_ = nullptr;
    host_ptr_ = nullptr;
    internal_allocation_ = false;
  }

  inline T *data() { return device_ptr_; }
  inline const T *data() const { return device_ptr_; }
};

template <typename T> class grid_info {
private:
  std::array<int, 3> full_size_;
  std::array<int, 3> local_size_;
  // origin of the local part of the grid in grid point
  std::array<int, 3> lower_corner_;
  std::array<int, 3> border_width_;
  std::array<T, 9> dh_;
  std::array<T, 9> dh_inv_;
  bool orthorhombic_{false};
  bool is_distributed_{false};
  gpu_vector<T> grid_;

public:
  grid_info(const grid_info &) = delete;
  grid_info &operator=(const grid_info &) = delete;

  grid_info(grid_info &&) noexcept = default;
  grid_info &operator=(grid_info &&) noexcept = default;

  grid_info(){};

  grid_info(const int *full_size__, const int *local_size__,
            const int *border_width__) {
    int roffset__[3] = {0, 0, 0};
    initialize(full_size__, local_size__, roffset__, border_width__);
  }

  ~grid_info() = default;

  inline void copy_to_gpu(const T *data, offloadStream_t &stream) {
    assert(data != nullptr);
    grid_.copy_to_gpu(data, stream);
  }

  inline void copy_to_gpu(offloadStream_t &stream) {
    grid_.copy_associated_host_to_gpu(stream);
  }

  inline void reset() { grid_.reset(); }

  /*
   * We do not allocate memory as the buffer is always coming from the outside
   * world. We only initialize the sizes, etc...
   */
  inline void resize(const int *full_size__, const int *local_size__,
                     const int *const roffset__,
                     const int *const border_width__) {
    initialize(full_size__, local_size__, roffset__, border_width__);
  }

  inline size_t size() const { return grid_.size(); }

  inline void zero(offloadStream_t &stream) { grid_.zero(stream); }

  inline void set_lattice_vectors(const T *dh__, const T *dh_inv__) {
    for (int i = 0; i < 9; ++i) {
      dh_[i] = dh__[i];
      dh_inv_[i] = dh_inv__[i];
    }
  }

  inline void is_distributed(const bool distributed__) {
    is_distributed_ = distributed__;
  }

  void check_orthorhombicity(const bool ortho) {
    if (ortho) {
      orthorhombic_ = true;
      return;
    }
    T norm1, norm2, norm3;
    bool orthogonal[3] = {false, false, false};
    norm1 = dh_[0] * dh_[0] + dh_[1] * dh_[1] + dh_[2] * dh_[2];
    norm2 = dh_[3] * dh_[3] + dh_[4] * dh_[4] + dh_[5] * dh_[5];
    norm3 = dh_[6] * dh_[6] + dh_[7] * dh_[7] + dh_[8] * dh_[8];

    norm1 = 1.0 / sqrt(norm1);
    norm2 = 1.0 / sqrt(norm2);
    norm3 = 1.0 / sqrt(norm3);

    /* x z */
    orthogonal[0] =
        ((fabs(dh_[0] * dh_[6] + dh_[1] * dh_[7] + dh_[2] * dh_[8]) * norm1 *
          norm3) < 1e-12);
    /* y z */
    orthogonal[1] =
        ((fabs(dh_[3] * dh_[6] + dh_[4] * dh_[7] + dh_[5] * dh_[8]) * norm2 *
          norm3) < 1e-12);
    /* x y */
    orthogonal[2] =
        ((fabs(dh_[0] * dh_[3] + dh_[1] * dh_[4] + dh_[2] * dh_[5]) * norm1 *
          norm2) < 1e-12);

    orthorhombic_ = orthogonal[0] && orthogonal[1] && orthogonal[2];
  }

  inline void copy_to_host(T *data__, offloadStream_t &stream) {
    assert(data__ != nullptr);
    grid_.copy_from_gpu(data__, stream);
  }

  inline void copy_to_host(offloadStream_t &stream) {
    grid_.copy_gpu_to_associated_host(stream);
  }

  inline void associate(void *host_ptr__, void *device_ptr__,
                        const size_t size__) {
    assert(host_ptr__ != nullptr);
    assert(device_ptr__ != nullptr);
    grid_.associate(host_ptr__, device_ptr__, size__);
  }
  inline bool is_distributed() { return is_distributed_; }

  inline int full_size(const int i) {
    assert(i < 3);
    return full_size_[i];
  }

  inline int local_size(const int i) {
    assert(i < 3);
    return local_size_[i];
  }

  inline int lower_corner(const int i) {
    assert(i < 3);
    return lower_corner_[i];
  }

  inline int border_width(const int i) {
    assert(i < 3);
    return border_width_[i];
  }

  inline T *data() { return grid_.data(); }
  inline const T *data() const { return grid_.data(); }

  inline gpu_vector<T> &grid() { return grid_; }
  inline const gpu_vector<T> &grid() const { return grid_; }

  inline T *dh() { return dh_.data(); }
  inline const T *dh() const { return dh_.data(); }

  inline T *dh_inv() { return dh_inv_.data(); }
  inline const T *dh_inv() const { return dh_inv_.data(); }

  inline bool is_orthorhombic() const { return orthorhombic_; }
  inline bool is_distributed() const { return is_distributed_; }

private:
  void initialize(const int *const full_size__, const int *const local_size__,
                  const int *const roffset__, const int *const border_width__) {
    // the calling code store things like this cube[z][y][x] (in fortran
    // cube(x,y,z)) so all sizes are [x,y,z] while we are working in C/C++ so we
    // have to permute the indices to get this right.

    full_size_[2] = full_size__[0];
    full_size_[1] = full_size__[1];
    full_size_[0] = full_size__[2];

    local_size_[2] = local_size__[0];
    local_size_[1] = local_size__[1];
    local_size_[0] = local_size__[2];

    lower_corner_[0] = roffset__[2];
    lower_corner_[1] = roffset__[1];
    lower_corner_[2] = roffset__[0];

    is_distributed_ = (full_size_[2] != local_size_[2]) ||
                      (full_size_[1] != local_size_[1]) ||
                      (full_size_[0] != local_size_[0]);

    border_width_[2] = border_width__[0];
    border_width_[1] = border_width__[1];
    border_width_[0] = border_width__[2];
  }
};

/*******************************************************************************
 * \brief Internal representation of a task.
 ******************************************************************************/
struct task_info {
  int level;
  int iatom;
  int jatom;
  int iset;
  int jset;
  int ipgf;
  int jpgf;
  int ikind, jkind;
  int border_mask;
  int block_num;
  int max_cab_size{0};
  double radius;
  double ra[3], rb[3], rp[3];
  double rab2;
  double zeta, zetb, zetp, prefactor, off_diag_twice;
  double rab[3];
  int lp_max{0};
  size_t coef_offset{0};
  int la_max, lb_max, la_min, lb_min, first_coseta, first_cosetb, ncoseta,
      ncosetb, nsgfa, nsgfb, nsgf_seta, nsgf_setb, maxcoa, maxcob;
  int sgfa, sgfb, subblock_offset;
  double3 roffset;
  int3 cube_size;
  int3 lb_cube;
  int3 cube_center;
  double discrete_radius;
  bool apply_border_mask;
  bool block_transposed;
  bool skip_task;
};

/*******************************************************************************
 * \brief Parameters of the collocate kernel.
 ******************************************************************************/

struct kernel_params {
  int first_task{0};
  // max size of cab.
  int cab_size_{0};
  int grid_full_size_[3] = {0, 0, 0};
  int grid_local_size_[3] = {0, 0, 0};
  int grid_lower_corner_[3] = {0, 0, 0};
  int grid_border_width_[3] = {0, 0, 0};
  double dh_[9];
  double dh_inv_[9];
  task_info *tasks;
  int *block_offsets{nullptr};
  char la_min_diff{0};
  char lb_min_diff{0};
  char la_max_diff{0};
  char lb_max_diff{0};
  enum grid_func func;
  double *ptr_dev[7] = {nullptr, nullptr, nullptr, nullptr,
                        nullptr, nullptr, nullptr};
  double **sphi_dev{nullptr};
  int ntasks{0};
  int *task_sorted_by_blocks_dev{nullptr};
  int *sorted_blocks_offset_dev{nullptr};
  int *num_tasks_per_block_dev{nullptr};
  int *cab_block_offset_dev{nullptr};
};

/* regroup all information about the context. */
class context_info {
private:
  int device_id_{-1};
  int lmax_{0};
  unsigned int checksum_{0};

public:
  int ntasks{0};
  int nlevels{0};
  int natoms{0};
  int nkinds{0};
  int nblocks{0};
  int cab_size_per_block_{0};
  std::vector<double *> sphi;
  std::vector<offloadStream_t> level_streams;
  offloadStream_t main_stream;
  int stats[2][20]; // [has_border_mask][lp]
  // all these tables are on the gpu. we can resize them copy to them and copy
  // from them
  gpu_vector<int> block_offsets_dev;
  gpu_vector<int> cab_block_offset_dev;
  gpu_vector<double> coef_dev_;
  gpu_vector<double> cab_dev_;
  gpu_vector<double> pab_block_;
  gpu_vector<double> hab_block_;
  gpu_vector<double> forces_;
  gpu_vector<double> virial_;
  gpu_vector<task_info> tasks_dev;
  gpu_vector<int> num_tasks_per_block_dev_;
  std::vector<grid_info<double>> grid_;
  std::vector<int> number_of_tasks_per_level_;
  std::vector<int> first_task_per_level_;
  std::vector<int> sphi_size;
  gpu_vector<double *> sphi_dev;
  gpu_vector<int> task_sorted_by_blocks_dev, sorted_blocks_offset_dev;
  bool calculate_forces{false};
  bool calculate_virial{false};
  bool compute_tau{false};
  bool apply_border_mask{false};
  // Default constructor
  context_info() = default;

  explicit context_info(int device_id__) {
    device_id_ = (device_id__ < 0) ? 0 : device_id__;
  }

  ~context_info() { clear(); }

  // Non-copyable
  context_info(const context_info &) = delete;
  context_info &operator=(const context_info &) = delete;

  // Movable
  context_info(context_info &&other) noexcept
      : device_id_{other.device_id_}, lmax_{other.lmax_},
        checksum_{other.checksum_}, ntasks{other.ntasks},
        nlevels{other.nlevels}, natoms{other.natoms}, nkinds{other.nkinds},
        nblocks{other.nblocks}, cab_size_per_block_{other.cab_size_per_block_},
        sphi{std::move(other.sphi)},
        level_streams{std::move(other.level_streams)},
        main_stream{other.main_stream},
        block_offsets_dev{std::move(other.block_offsets_dev)},
        cab_block_offset_dev{std::move(other.cab_block_offset_dev)},
        coef_dev_{std::move(other.coef_dev_)},
        cab_dev_{std::move(other.cab_dev_)},
        pab_block_{std::move(other.pab_block_)},
        hab_block_{std::move(other.hab_block_)},
        forces_{std::move(other.forces_)}, virial_{std::move(other.virial_)},
        tasks_dev{std::move(other.tasks_dev)},
        num_tasks_per_block_dev_{std::move(other.num_tasks_per_block_dev_)},
        grid_{std::move(other.grid_)},
        number_of_tasks_per_level_{std::move(other.number_of_tasks_per_level_)},
        first_task_per_level_{std::move(other.first_task_per_level_)},
        sphi_size{std::move(other.sphi_size)},
        sphi_dev{std::move(other.sphi_dev)},
        task_sorted_by_blocks_dev{std::move(other.task_sorted_by_blocks_dev)},
        sorted_blocks_offset_dev{std::move(other.sorted_blocks_offset_dev)},
        calculate_forces{other.calculate_forces},
        calculate_virial{other.calculate_virial},
        compute_tau{other.compute_tau},
        apply_border_mask{other.apply_border_mask} {

    // copy stats
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 20; ++j)
        stats[i][j] = other.stats[i][j];

    // leave other in a safe state (no ownership)
    other.device_id_ = -1;
    other.lmax_ = 0;
    other.checksum_ = 0;
    other.ntasks = 0;
    other.nlevels = 0;
    other.natoms = 0;
    other.nkinds = 0;
    other.nblocks = 0;
    other.main_stream = {};
    other.calculate_forces = false;
    other.calculate_virial = false;
    other.compute_tau = false;
    other.apply_border_mask = false;
  }

  context_info &operator=(context_info &&other) noexcept {
    if (this != &other) {
      clear(); // free current GPU resources and streams

      device_id_ = other.device_id_;
      lmax_ = other.lmax_;
      checksum_ = other.checksum_;
      ntasks = other.ntasks;
      nlevels = other.nlevels;
      natoms = other.natoms;
      nkinds = other.nkinds;
      nblocks = other.nblocks;
      sphi = std::move(other.sphi);
      level_streams = std::move(other.level_streams);
      main_stream = other.main_stream;

      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 20; ++j)
          stats[i][j] = other.stats[i][j];

      block_offsets_dev = std::move(other.block_offsets_dev);
      cab_block_offset_dev = std::move(other.cab_block_offset_dev);
      coef_dev_ = std::move(other.coef_dev_);
      cab_dev_ = std::move(other.cab_dev_);
      pab_block_ = std::move(other.pab_block_);
      hab_block_ = std::move(other.hab_block_);
      forces_ = std::move(other.forces_);
      virial_ = std::move(other.virial_);
      tasks_dev = std::move(other.tasks_dev);
      num_tasks_per_block_dev_ = std::move(other.num_tasks_per_block_dev_);
      grid_ = std::move(other.grid_);
      number_of_tasks_per_level_ = std::move(other.number_of_tasks_per_level_);
      first_task_per_level_ = std::move(other.first_task_per_level_);
      sphi_size = std::move(other.sphi_size);
      sphi_dev = std::move(other.sphi_dev);
      task_sorted_by_blocks_dev = std::move(other.task_sorted_by_blocks_dev);
      sorted_blocks_offset_dev = std::move(other.sorted_blocks_offset_dev);

      calculate_forces = other.calculate_forces;
      calculate_virial = other.calculate_virial;
      compute_tau = other.compute_tau;
      apply_border_mask = other.apply_border_mask;

      // reset other
      other.device_id_ = -1;
      other.lmax_ = 0;
      other.checksum_ = 0;
      other.ntasks = 0;
      other.nlevels = 0;
      other.natoms = 0;
      other.nkinds = 0;
      other.nblocks = 0;
      other.main_stream = {};
      other.calculate_forces = false;
      other.calculate_virial = false;
      other.compute_tau = false;
      other.apply_border_mask = false;
    }
    return *this;
  }

  void clear() {
    if (device_id_ >= 0) {
      offload_set_chosen_device(device_id_);
      offload_activate_chosen_device();
    }

    tasks_dev.reset();
    block_offsets_dev.reset();
    cab_block_offset_dev.reset();
    coef_dev_.reset();
    cab_dev_.reset();
    task_sorted_by_blocks_dev.reset();
    sorted_blocks_offset_dev.reset();
    sphi_dev.reset();
    forces_.reset();
    virial_.reset();

    for (auto &phi : sphi) {
      if (phi != nullptr)
        offloadFree(phi);
    }
    sphi.clear();

    if (main_stream) {
      offloadStreamDestroy(main_stream);
      main_stream = {};
    }

    for (auto &stream : level_streams) {
      if (stream)
        offloadStreamDestroy(stream);
    }
    level_streams.clear();

    for (auto &g : grid_) {
      g.reset();
    }
    grid_.clear();
  }

  int lmax() const { return lmax_; }

  void initialize_basis_sets(const grid_basis_set **basis_sets,
                             const int nkinds__) {
    nkinds = nkinds__;
    if (nkinds__ > (int)sphi.size()) {
      for (auto &phi : sphi)
        if (phi != nullptr) {
          offloadFree(phi);
        }

      sphi_dev.resize(nkinds__);

      sphi.resize(nkinds__, nullptr);
      sphi_size.clear();
      sphi_size.resize(nkinds__, 0);
      sphi_dev.resize(nkinds__);
    }

    // Upload basis sets to device.
    for (int i = 0; i < nkinds__; i++) {
      const auto &basis_set = basis_sets[i];
      if (sphi_size[i] < basis_set->nsgf * basis_set->maxco) {
        offloadMalloc((void **)&sphi[i],
                      basis_set->nsgf * basis_set->maxco * sizeof(double));
        sphi_size[i] = basis_set->nsgf * basis_set->maxco;
      }
      //      offloadMemset(sphi[i], 0, sizeof(double) * sphi_size[i]);
      offloadMemcpyHtoD(sphi[i], basis_set->sphi,
                        basis_set->nsgf * basis_set->maxco * sizeof(double));
    }
    sphi_dev.copy_to_gpu(sphi);
    // Find largest angular momentum.
    lmax_ = 0;
    for (int ikind = 0; ikind < nkinds; ikind++) {
      for (int iset = 0; iset < basis_sets[ikind]->nset; iset++) {
        lmax_ = std::max(lmax_, basis_sets[ikind]->lmax[iset]);
      }
    }
  }

  void create_streams() {
    // allocate main hip stream
    offloadStreamCreate(&main_stream);

    // allocate one hip stream per grid level
    if ((int)level_streams.size() < nlevels) {
      level_streams.resize(nlevels);
      for (auto &stream : level_streams) {
        offloadStreamCreate(&stream);
      }
    }
  }

  void synchronize(offloadStream_t &stream) {
    offloadStreamSynchronize(stream);
  }

  void synchornize() {
    // wait for all the streams to finish
    offloadDeviceSynchronize();
  }

  void set_device() {
    offload_set_chosen_device(device_id_);
    offload_activate_chosen_device();
  }

  void collocate_one_grid_level(const int level, const enum grid_func func,
                                int *lp_diff);
  void integrate_one_grid_level(const int level, int *lp_diff);
  void compute_hab_coefficients();
  /* basic checksum computation for simple verification that the object is sane
   */
  void compute_checksum() { checksum_ = compute_checksum_(); }
  void verify_checksum() {
    if (checksum_ != compute_checksum_()) {
      fprintf(stderr, "This object does not seem to have the right structure.\n"
                      "A casting went wrong or the object is corrupted\n");
      abort();
    }
  }
  void calculate_all_coefficients(const enum grid_func func, int *lp_diff);

private:
  kernel_params set_kernel_parameters(const int level,
                                      const smem_parameters &smem_params);
  unsigned int compute_checksum_() {
    return natoms ^ ntasks ^ nlevels ^ nkinds ^ nblocks ^ 0x4F2C5D1A;
  }
};
} // namespace rocm_backend
#endif
