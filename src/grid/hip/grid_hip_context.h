/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/*
 * Authors :
   - Dr Mathieu Taillefumier (ETH Zurich / CSCS)
   - Advanced Micro Devices, Inc.
*/

#ifndef GRID_HIP_CONTEXT_H
#define GRID_HIP_CONTEXT_H

#include <hip/hip_runtime_api.h>
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

class smem_parameters;
template <typename T> class gpu_vector {
  size_t allocated_size_{0};
  size_t current_size_{0};
  T *ptr = nullptr;

public:
  gpu_vector() {
    ptr = nullptr;
    allocated_size_ = 0;
    current_size_ = 0;
  }
  // size is the number of elements not the memory size
  gpu_vector(const size_t size_) {
    if (size_ < 16) {
      allocated_size_ = 16;
    } else {
      allocated_size_ = (size_ / 16 + 1) * 16;
    }
    current_size_ = size_;

    offloadMalloc((void **)&ptr, sizeof(T) * allocated_size_);
  }

  ~gpu_vector() { reset(); }

  inline size_t size() { return current_size_; }

  inline void copy_to_gpu(const T *data__) {
    offloadMemcpyHtoD(ptr, data__, sizeof(T) * current_size_);
  }

  inline void copy_to_gpu(const T *data__, offloadStream_t &stream) {
    offloadMemcpyAsyncHtoD(ptr, data__, sizeof(T) * current_size_, stream);
  }

  inline void copy_from_gpu(T *data__, offloadStream_t &stream) {
    offloadMemcpyAsyncDtoH(data__, ptr, sizeof(T) * current_size_, stream);
  }

  inline void zero(offloadStream_t &stream) {
    // zero device grid buffers
    offloadMemsetAsync(ptr, 0, sizeof(T) * current_size_, stream);
  }

  inline void zero() {
    // zero device grid buffers
    offloadMemset(ptr, 0, sizeof(T) * current_size_);
  }

  inline void copy_to_gpu(const std::vector<T> &data__) {
    assert(data__.size() == current_size_);
    // if it fails it means that the vector on the gpu does not have the right
    // size. two option then
    // - resize the gpu vector
    // - or the cpu vector and gpu vector are not representing the quantity.

    offloadMemcpyHtoD(ptr, data__.data(), sizeof(T) * data__.size());
  }

  inline void resize(const size_t new_size_) {
    if (allocated_size_ < new_size_) {
      if (ptr != nullptr)
        offloadFree(ptr);
      allocated_size_ = (new_size_ / 16 + (new_size_ % 16 != 0)) * 16;
      offloadMalloc((void **)&ptr, sizeof(T) * allocated_size_);
    }
    current_size_ = new_size_;
  }

  // does not invalidate the pointer. The memory is still allocated
  inline void clear() { current_size_ = 0; }

  // reset the class and free memory
  inline void reset() {
    if (ptr != nullptr)
      offloadFree(ptr);

    allocated_size_ = 0;
    current_size_ = 0;
    ptr = nullptr;
  }

  inline T *data() { return ptr; }
};

template <typename T> class grid_info {
  int full_size_[3] = {0, 0, 0};
  int local_size_[3] = {0, 0, 0};
  // origin of the local part of the grid in grid point
  int lower_corner_[3] = {0, 0, 0};
  int border_width_[3] = {0, 0, 0};
  double dh_[9];
  double dh_inv_[9];
  bool orthorhombic_{false};
  bool is_distributed_{false};
  gpu_vector<T> grid_;

public:
  grid_info(){};

  grid_info(const int *full_size__, const int *local_size__,
            const int *border_width__) {
    initialize(full_size__, local_size__, border_width__);
  }

  ~grid_info() { grid_.reset(); };

  inline T *data() { return grid_.data(); }

  inline void copy_to_gpu(const T *data, offloadStream_t &stream) {
    grid_.copy_to_gpu(data, stream);
  }

  inline void reset() { grid_.reset(); }

  inline void resize(const int *full_size__, const int *local_size__,
                     const int *const roffset__,
                     const int *const border_width__) {
    initialize(full_size__, local_size__, roffset__, border_width__);
  }

  inline size_t size() const { return grid_.size(); }

  inline void zero(offloadStream_t &stream) { grid_.zero(stream); }
  inline gpu_vector<T> &grid() { return grid_; }
  inline void set_lattice_vectors(const double *dh__, const double *dh_inv__) {
    memcpy(dh_, dh__, sizeof(double) * 9);
    memcpy(dh_inv_, dh_inv__, sizeof(double) * 9);
  }

  inline T *dh() { return dh_; }

  inline T *dh_inv() { return dh_inv_; }

  inline bool is_orthorhombic() { return orthorhombic_; }

  inline void is_distributed(const bool distributed__) {
    is_distributed_ = distributed__;
  }

  void check_orthorhombicity(const bool ortho) {
    if (ortho) {
      orthorhombic_ = true;
      return;
    }
    double norm1, norm2, norm3;
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

  inline void copy_to_host(double *data__, offloadStream_t &stream) {
    grid_.copy_from_gpu(data__, stream);
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

    grid_.resize(local_size_[0] * local_size_[1] * local_size_[2]);
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
  int smem_alpha_offset{0};
  int smem_cab_offset{0};
  int first_task{0};
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
  double *ptr_dev[6] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  double **sphi_dev{nullptr};
  int ntasks{0};
  int *task_sorted_by_blocks_dev{nullptr};
  int *sorted_blocks_offset_dev{nullptr};
  int *num_tasks_per_block_dev{nullptr};
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
  std::vector<double *> sphi;
  std::vector<offloadStream_t> level_streams;
  offloadStream_t main_stream;
  int stats[2][20]; // [has_border_mask][lp]
  // all these tables are on the gpu. we can resize them copy to them and copy
  // from them
  gpu_vector<int> block_offsets_dev;
  gpu_vector<double> coef_dev_;
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

  context_info() {}
  context_info(const int device_id__) {
    if (device_id__ < 0)
      device_id_ = 0;
    else
      device_id_ = device_id__;
  }
  ~context_info() { clear(); }

  void clear() {
    hipSetDevice(device_id_);
    tasks_dev.reset();
    block_offsets_dev.reset();
    coef_dev_.reset();
    task_sorted_by_blocks_dev.reset();
    sorted_blocks_offset_dev.reset();
    sphi_dev.reset();
    forces_.reset();
    virial_.reset();
    for (auto &phi : sphi)
      if (phi != nullptr)
        offloadFree(phi);
    sphi.clear();

    offloadStreamDestroy(main_stream);

    for (int i = 0; i < nlevels; i++) {
      offloadStreamDestroy(level_streams[i]);
    }
    level_streams.clear();

    for (auto &grid : grid_) {
      grid.reset();
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
      offloadMemset(sphi[i], 0, sizeof(double) * sphi_size[i]);
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

  void set_device() { hipSetDevice(device_id_); }

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

private:
  kernel_params set_kernel_parameters(const int level,
                                      const smem_parameters &smem_params);
  unsigned int compute_checksum_() {
    return natoms ^ ntasks ^ nlevels ^ nkinds ^ nblocks ^ 0x4F2C5D1A;
  }
};
} // namespace rocm_backend
#endif
