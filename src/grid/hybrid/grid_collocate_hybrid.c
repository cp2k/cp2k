/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifdef __GRID_CUDA

#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/grid_common.h"
#include "../common/grid_library.h"
#include "../cpu/coefficients.h"
#include "../cpu/collocation_integration.h"
#include "../cpu/cpu_private_header.h"
#include "../cpu/grid_context_cpu.h"
#include "../cpu/grid_prepare_pab_dgemm.h"
#include "../cpu/tensor_local.h"
#include "../cpu/utils.h"

typedef struct grid_context_ grid_hybrid_task_list;

void initialize_grid_parameters_on_gpu_step1(void *const ctx, const int level);

void compute_collocation_gpu(pgf_list_gpu *handler);

void rotate_to_cartesian_harmonics(const grid_basis_set *ibasis,
                                   const grid_basis_set *jbasis,
                                   const int iatom, const int jatom,
                                   const int iset, const int jset,
                                   double *const block, tensor *work,
                                   tensor *pab);

static pgf_list_gpu *create_worker_list(const int number_of_workers,
                                        const int batch_size, int device_id,
                                        const bool apply_cutoff);

void grid_collocate(collocation_integration *const handler,
                    const bool use_ortho, const double zetp, const double rp[3],
                    const double radius);
static void destroy_worker_list(pgf_list_gpu *const list);

static void reset_list_gpu(pgf_list_gpu *const list) {
  cudaSetDevice(list->device_id);
  list->list_length = 0;
  list->coef_dynamic_alloc_size_gpu_ = 0;
}

static void my_worker_is_running(pgf_list_gpu *const my_worker) {
  if (!my_worker)
    return;

  if (my_worker->running) {
    cudaSetDevice(my_worker->device_id);
    cudaStreamSynchronize(my_worker->stream);
    reset_list_gpu(my_worker);
    my_worker->running = false;
  }
}

static void add_orbital_to_list(pgf_list_gpu *const list, const int cmax,
                                const _task *task, const tensor *const coef) {
  assert(list->batch_size > list->list_length);
  list->coef_offset_cpu_[0] = 0;
  const int lp = coef->size[2] - 1;
  const int coef_alloc_size_ = ((lp + 1) * (lp + 2) * (lp + 3)) / 6;
  if (list->list_length > 0) {
    list->coef_offset_cpu_[list->list_length] =
        list->coef_offset_cpu_[list->list_length - 1] +
        list->coef_previous_alloc_size_;
    list->coef_dynamic_alloc_size_gpu_ += coef_alloc_size_;
  } else {
    list->coef_dynamic_alloc_size_gpu_ = coef_alloc_size_;
  }
  list->coef_previous_alloc_size_ = coef_alloc_size_;

  if (list->coef_dynamic_alloc_size_gpu_ > list->coef_alloc_size_gpu_) {
    list->durty = true;
    list->coef_cpu_ =
        realloc(list->coef_cpu_,
                sizeof(double) * (list->coef_dynamic_alloc_size_gpu_ + 4096));
    list->coef_alloc_size_gpu_ = list->coef_dynamic_alloc_size_gpu_ + 4096;
  }

  if (list->lmax < lp)
    list->lmax = lp;
  transform_xyz_to_triangular(
      coef, list->coef_cpu_ + list->coef_offset_cpu_[list->list_length]);

  list->cmax = imax(list->cmax, cmax);
  memcpy(&list->task_list_cpu_[list->list_length], task, sizeof(_task));

  list->list_length++;
}

static pgf_list_gpu *create_worker_list(const int number_of_workers,
                                        const int batch_size,
                                        const int device_id,
                                        const bool apply_cutoff) {
  pgf_list_gpu *list =
      (pgf_list_gpu *)malloc(sizeof(pgf_list_gpu) * number_of_workers);
  memset(list, 0, sizeof(pgf_list_gpu) * number_of_workers);
  for (int i = 0; i < number_of_workers; i++) {
    list[i].device_id = device_id;
    list[i].lmax = -1;
    list[i].batch_size = batch_size;
    list[i].list_length = 0;
    list[i].coef_dynamic_alloc_size_gpu_ = 0;
    list[i].task_list_cpu_ = (_task *)malloc(sizeof(_task) * list->batch_size);
    list[i].coef_offset_cpu_ = (int *)malloc(sizeof(int) * list->batch_size);
    list[i].coef_cpu_ =
        (double *)malloc(sizeof(double) * list->batch_size * 8 * 8 * 8);
    list[i].coef_alloc_size_gpu_ = list->batch_size * 8 * 8 * 8;
    cudaSetDevice(list[i].device_id);

    cudaMalloc((void **)&list[i].coef_offset_gpu_,
               sizeof(int) * list->batch_size);
    cudaMalloc((void **)&list[i].coef_gpu_,
               sizeof(double) * list->coef_alloc_size_gpu_);
    cudaMalloc((void **)&list[i].task_list_gpu_,
               sizeof(_task) * list->batch_size);

    cudaStreamCreate(&list[i].stream);
    cublasCreate(&list[i].blas_handle);
    cublasSetStream(list[i].blas_handle, list[i].stream);
    cudaEventCreate(&list[i].event);
    list[i].next = list + i + 1;
    list[i].running = false;
    list[i].apply_cutoff = apply_cutoff;
  }

  list[number_of_workers - 1].next = NULL;

  return list;
}

static void destroy_worker_list(pgf_list_gpu *const lst) {
  cudaSetDevice(lst->device_id);
  cudaFree(lst->coef_offset_gpu_);
  cudaFree(lst->coef_gpu_);
  cudaFree(lst->data_gpu_);
  cudaFree(lst->task_list_gpu_);
  cudaStreamDestroy(lst->stream);
  cublasDestroy(lst->blas_handle);
  cudaEventDestroy(lst->event);
  free(lst->coef_offset_cpu_);
  free(lst->coef_cpu_);
  free(lst->task_list_cpu_);
}

static void
initialize_grid_parameters_on_gpu(collocation_integration *const handler) {
  for (int worker = 0; worker < handler->worker_list_size; worker++) {
    assert(handler->worker_list[worker].device_id >= 0);
    cudaSetDevice(handler->worker_list[worker].device_id);

    handler->worker_list[worker].grid_size.x = handler->grid.size[2];
    handler->worker_list[worker].grid_size.y = handler->grid.size[1];
    handler->worker_list[worker].grid_size.z = handler->grid.size[0];

    handler->worker_list[worker].grid_full_size.x = handler->grid.full_size[2];
    handler->worker_list[worker].grid_full_size.y = handler->grid.full_size[1];
    handler->worker_list[worker].grid_full_size.z = handler->grid.full_size[0];

    handler->worker_list[worker].window_size.x = handler->grid.window_size[2];
    handler->worker_list[worker].window_size.y = handler->grid.window_size[1];
    handler->worker_list[worker].window_size.z = handler->grid.window_size[0];

    handler->worker_list[worker].window_shift.x = handler->grid.window_shift[2];
    handler->worker_list[worker].window_shift.y = handler->grid.window_shift[1];
    handler->worker_list[worker].window_shift.z = handler->grid.window_shift[0];

    handler->worker_list[worker].grid_lower_corner_position.x =
        handler->grid.lower_corner[2];
    handler->worker_list[worker].grid_lower_corner_position.y =
        handler->grid.lower_corner[1];
    handler->worker_list[worker].grid_lower_corner_position.z =
        handler->grid.lower_corner[0];

    if (handler->worker_list[worker].data_gpu_ == NULL) {
      cudaMalloc((void **)&handler->worker_list[worker].data_gpu_,
                 sizeof(double) * handler->grid.alloc_size_);
      handler->worker_list[worker].data_gpu_old_size_ =
          handler->grid.alloc_size_;
    } else {
      if (handler->worker_list[worker].data_gpu_old_size_ <
          handler->grid.alloc_size_) {
        cudaFree(handler->worker_list[worker].data_gpu_);
        cudaMalloc((void **)&handler->worker_list[worker].data_gpu_,
                   sizeof(double) * handler->grid.alloc_size_);
        handler->worker_list[worker].data_gpu_old_size_ =
            handler->grid.alloc_size_;
      }
      handler->worker_list[worker].zeroing_grid = true;
    }

    cudaMemsetAsync(handler->worker_list[worker].data_gpu_, 0,
                    sizeof(double) * handler->grid.alloc_size_,
                    handler->worker_list[worker].stream);
    reset_list_gpu(handler->worker_list + worker);
  }
}

void release_gpu_resources(collocation_integration *handler) {
  for (int i = 0; i < handler->worker_list_size; i++) {
    destroy_worker_list(handler->worker_list + i);
  }
  free(handler->worker_list);
}

static void initialize_worker_list_on_gpu(collocation_integration *handle,
                                          const int device_id,
                                          const int number_of_gaussian,
                                          const int number_of_workers) {
  assert(handle != NULL);
  handle->worker_list_size = number_of_workers;
  handle->number_of_gaussian = number_of_gaussian;

  handle->lmax = -1;

  /* we can inclrease this afterwards */
  /* right now only one list */

  handle->worker_list = create_worker_list(
      number_of_workers, number_of_gaussian, device_id, handle->apply_cutoff);
}

//******************************************************************************
// \brief Collocate a range of tasks which are destined for the same grid level.
// \author Ole Schuett
//******************************************************************************
static void collocate_one_grid_level_gpu(grid_context *const ctx,
                                         const int *border_width,
                                         const enum grid_func func,
                                         const int level,
                                         const grid_buffer *pab_blocks) {
  // how many threads have participated in the omp region
  tensor *const grid = &ctx->grid[level];
  const int max_threads = omp_get_max_threads();

  // Using default(shared) because with GCC 9 the behavior around const changed:
  // https://www.gnu.org/software/gcc/gcc-9/porting_to.html
#pragma omp parallel default(shared)
  {
    const int thread_id = omp_get_thread_num();

    struct collocation_integration_ *handler = ctx->handler[thread_id];

    tensor_copy(&handler->grid, grid);

    handler->func = func;
    grid_prepare_get_ldiffs_dgemm(handler->func, handler->lmin_diff,
                                  handler->lmax_diff);
    initialize_basis_vectors(handler, (const double(*)[3])grid->dh,
                             (const double(*)[3])grid->dh_inv);

    handler->apply_cutoff = ctx->apply_cutoff;

    for (int d = 0; d < 3; d++)
      handler->orthogonal[d] = handler->grid.orthogonal[d];

    tensor_copy(&handler->grid, grid);
    initialize_grid_parameters_on_gpu(handler);

    if (max_threads == 1) {
      // only one thread in the omp region.
      handler->grid.data = grid->data;
    } else {
      // one handler per omp thread
      handler->grid.data =
          ((double *)ctx->scratch) + thread_id * grid->alloc_size_;
      memset(handler->grid.data, 0, sizeof(double) * handler->grid.alloc_size_);
    }

    /* I have two workers per thread. one doing computation while the other get
     * initialized. Then I do a ping pong between them */
    pgf_list_gpu *my_worker_1 = &handler->worker_list[0];
    pgf_list_gpu *my_worker_2 = &handler->worker_list[1];
    pgf_list_gpu *current_worker = my_worker_1;
    /* the one doing computations on gpu */
    pgf_list_gpu *running_worker = NULL;

    reset_list_gpu(my_worker_1);
    reset_list_gpu(my_worker_2);
    my_worker_1->running = false;
    my_worker_2->running = false;
    memcpy(&my_worker_1->border_width, border_width, sizeof(int3));
    memcpy(&my_worker_2->border_width, border_width, sizeof(int3));
    tensor work, pab, pab_prep;

    // Allocate pab matrix for re-use across tasks.
    initialize_tensor_2(&pab, ctx->maxco, ctx->maxco);
    alloc_tensor(&pab);

    initialize_tensor_2(&work, ctx->maxco, ctx->maxco);
    alloc_tensor(&work);

    initialize_tensor_2(&pab_prep, ctx->maxco, ctx->maxco);
    alloc_tensor(&pab_prep);

    const _task *prev_task = NULL;
#pragma omp for schedule(static)
    for (int itask = 0; itask < ctx->tasks_per_level[level]; itask++) {
      // Define some convenient aliases.
      _task *const task = &ctx->tasks[level][itask];
      if (task->level != level) {
        printf("level %d, %d\n", task->level, level);
        abort();
      }

      compute_coefficients(ctx, handler, prev_task, task, pab_blocks, &pab,
                           &work, &pab_prep);

      {
        int cubecenter[3];
        int cube_size[3];
        int lb_cube[3], ub_cube[3];
        double roffset[3];
        double disr_radius;

        /* seting up the cube parameters */
        int cmax = compute_cube_properties(
            handler->orthogonal[0] && handler->orthogonal[1] &&
                handler->orthogonal[2],
            task->radius, handler->dh, handler->dh_inv, task->rp, &disr_radius,
            roffset, cubecenter, lb_cube, ub_cube, cube_size);
        task->l1_plus_l2_ = handler->coef.size[2] - 1;
        add_orbital_to_list(current_worker, cmax, task, &handler->coef);
      }
      /* The list is full so we can start computation on GPU */
      if (current_worker->batch_size == current_worker->list_length) {
        current_worker->running = true;
        compute_collocation_gpu(current_worker);

        /* We need to be sure that the other queue is empty or wait for it to
         * finish */
        my_worker_is_running(running_worker);

        /* now swap the two queues */
        if (!running_worker) {
          running_worker = current_worker;
          current_worker = my_worker_2;
        } else {
          pgf_list_gpu *tmp = running_worker;
          running_worker = current_worker;
          current_worker = tmp;
        }
      }
    }

    /* We may have a partial batch to do so we must complete them before doing
     * anything else */
    if (current_worker->list_length) {
      /* ensure that the potential work running on the stream is done
       * before running the next batch */
      current_worker->running = true;
      compute_collocation_gpu(current_worker);
    }

    my_worker_is_running(current_worker);
    my_worker_is_running(running_worker);

    // Merge worker grids to worker 0

    double alpha = 1.0;
    cublasDaxpy(handler->worker_list[0].blas_handle, handler->grid.alloc_size_,
                &alpha, handler->worker_list[1].data_gpu_, 1,
                handler->worker_list[0].data_gpu_, 1);

    free(pab.data);
    free(pab_prep.data);
    free(work.data);

    handler->grid.data = NULL;
  }

  if (ctx->number_of_devices == 1) {
    double alpha = 1.0;
    cudaStreamSynchronize(ctx->handler[0]->worker_list[0].stream);
    for (int thread = 1; thread < max_threads; thread++) {
      cudaStreamSynchronize(ctx->handler[thread]->worker_list[0].stream);
      cublasDaxpy(ctx->handler[0]->worker_list[0].blas_handle,
                  ctx->handler[0]->grid.alloc_size_, &alpha,
                  ctx->handler[thread]->worker_list[0].data_gpu_, 1,
                  ctx->handler[0]->worker_list[0].data_gpu_, 1);
    }
    cudaStreamSynchronize(ctx->handler[0]->worker_list->stream);
    cudaMemcpy(ctx->grid[level].data, ctx->handler[0]->worker_list->data_gpu_,
               sizeof(double) * ctx->handler[0]->grid.alloc_size_,
               cudaMemcpyDeviceToHost);
  } else {
    int *dev_thread = (int *)malloc(sizeof(int) * ctx->number_of_devices);

    for (int dev = 0; dev < ctx->number_of_devices; dev++)
      dev_thread[dev] = -1;

    for (int thread = 1; thread < max_threads; thread++) {
      if (dev_thread[ctx->handler[thread]->worker_list[0].device_id] == -1)
        dev_thread[ctx->handler[thread]->worker_list[0].device_id] = thread;
    }

    for (int thread = 1; thread < max_threads; thread++) {
      double alpha = 1.0;
      /* two different handlers use the same gpu */
      cudaSetDevice(ctx->handler[thread]->worker_list[0].device_id);
      const int thread_reduce =
          dev_thread[ctx->handler[thread]->worker_list[0].device_id];
      cublasDaxpy(ctx->handler[thread_reduce]->worker_list[0].blas_handle,
                  ctx->handler[thread_reduce]->grid.alloc_size_, &alpha,
                  ctx->handler[thread]->worker_list[0].data_gpu_, 1,
                  ctx->handler[thread_reduce]->worker_list[0].data_gpu_, 1);
    }

    /* final reduction over different gpus. Need explicit copy to cpu */

    for (int dev = 0; dev < ctx->number_of_devices; dev++) {
      cudaSetDevice(dev);
      const int thread_ = dev_thread[dev];
      cudaMemcpy(ctx->scratch, ctx->handler[thread_]->worker_list[0].data_gpu_,
                 sizeof(double) * ctx->grid[level].alloc_size_,
                 cudaMemcpyDeviceToHost);
      cblas_daxpy(ctx->grid[level].alloc_size_, 1.0, ctx->scratch, 1,
                  ctx->grid[level].data, 1);
    }
  }
}

void grid_hybrid_collocate_task_list(
    void *const ptr, const bool orthorhombic, const enum grid_func func,
    const int nlevels, const int npts_global[nlevels][3],
    const int npts_local[nlevels][3], const int shift_local[nlevels][3],
    const int border_width[nlevels][3], const double dh[nlevels][3][3],
    const double dh_inv[nlevels][3][3], const grid_buffer *pab_blocks,
    double *grid[nlevels]) {
  assert(ptr != NULL);
  grid_context *const ctx = (grid_context *const)ptr;
  const int max_threads = omp_get_max_threads();

  /* need to set the grid parameters first i need the total size for allocating
   * scratch memory */
  if (orthorhombic) {
    ctx->orthorhombic = true;
  } else {
    ctx->orthorhombic = false;
  }

  for (int level = 0; level < ctx->nlevels; level++) {
    set_grid_parameters(&ctx->grid[level], orthorhombic, npts_global[level],
                        npts_local[level], shift_local[level],
                        border_width[level], dh[level], dh_inv[level],
                        grid[level]);
  }

  if (ctx->scratch == NULL) {
    int max_size = ctx->grid[0].alloc_size_;

    for (int x = 1; x < nlevels; x++) {
      max_size = imax(ctx->grid[x].alloc_size_, max_size);
    }

    max_size = ((max_size / 4096) + 1) * 4096;

    ctx->scratch = memalign(4096, sizeof(double) * max_size * max_threads);
  }

  /*
    Each handler is assigned to a given device (possibly different) and a
    given omp thread so the code can already run on single node multiple
    devices if we treat the grid levels sequentially.

    it would be possible to add extra parallelization over the grid level but
    it will require some changes in the cuda code since constant memory on the
    device is used to store the grid size. Removing this is not complicated
    though, one extra parameter to the function call
  */

  int length_queue = 0;

  for (int level = 0; level < ctx->nlevels; level++) {
    length_queue = imax(length_queue, ctx->tasks_per_level[level]);
  }

  length_queue /= (2 * max_threads);

  length_queue++;

  length_queue = imin(length_queue, ctx->queue_length);

  for (int thread = 0; thread < max_threads; thread++) {
    initialize_worker_list_on_gpu(
        ctx->handler[thread], ctx->device_id[thread % ctx->number_of_devices],
        length_queue, /* basically the number of gaussian we can treat */
        2);           /* number of workers */
  }

  /* We can not dispatch the different grids over different GPU with this
   * implementation. */
  for (int level = 0; level < ctx->nlevels; level++) {
    bool orthogonal[3];
    verify_orthogonality(dh[level], orthogonal);
    if (orthogonal[0] && orthogonal[1] && orthogonal[2] && (!ctx->orthorhombic))
      ctx->orthorhombic = true;
    else
      ctx->orthorhombic = false;
    initialize_grid_parameters_on_gpu_step1(ctx, level);
    collocate_one_grid_level_gpu(ctx, border_width[level], func, level,
                                 pab_blocks);
  }

  /* data is *not* allocated when I create the ctx so need to set it up to null.
   */
#pragma omp parallel for
  for (int thread = 0; thread < max_threads; thread++)
    release_gpu_resources(ctx->handler[thread]);

  free(ctx->scratch);
  ctx->scratch = NULL;
}

/*******************************************************************************
 * \brief Allocates a task list for the hybrid backend.
 *        See grid_task_list.h for details.
 ******************************************************************************/
void grid_hybrid_create_task_list(
    const int ntasks, const int nlevels, const int natoms, const int nkinds,
    const int nblocks, const int block_offsets[nblocks],
    const double atom_positions[natoms][3], const int atom_kinds[natoms],
    const grid_basis_set *basis_sets[nkinds], const int level_list[ntasks],
    const int iatom_list[ntasks], const int jatom_list[ntasks],
    const int iset_list[ntasks], const int jset_list[ntasks],
    const int ipgf_list[ntasks], const int jpgf_list[ntasks],
    const int border_mask_list[ntasks], const int block_num_list[ntasks],
    const double radius_list[ntasks], const double rab_list[ntasks][3],
    grid_hybrid_task_list **task_list) {

  if (*task_list == NULL) {
    *task_list = create_grid_context_cpu(
        ntasks, nlevels, natoms, nkinds, nblocks, block_offsets, atom_positions,
        atom_kinds, basis_sets, level_list, iatom_list, jatom_list, iset_list,
        jset_list, ipgf_list, jpgf_list, border_mask_list, block_num_list,
        radius_list, rab_list);
  } else {
    update_grid_context_cpu(
        ntasks, nlevels, natoms, nkinds, nblocks, block_offsets, atom_positions,
        atom_kinds, basis_sets, level_list, iatom_list, jatom_list, iset_list,
        jset_list, ipgf_list, jpgf_list, border_mask_list, block_num_list,
        radius_list, rab_list, *task_list);
  }

  /* does not allocate anything on the GPU. */
  /* allocation only occurs when the collocate (or integrate) function is
   * called. Resources are released before exiting the function */

  const grid_library_config config = grid_library_get_config();

  /* I do *not* store the address of config.device_id */
  initialize_grid_context_on_gpu(*task_list, 1 /* number of devices */,
                                 &config.device_id);

  update_queue_length(*task_list, config.queue_length);

  if (config.apply_cutoff) {
    apply_cutoff(*task_list);
  }
}

/*******************************************************************************
 * \brief Deallocates given task list, basis_sets have to be freed separately.
 ******************************************************************************/
void grid_hybrid_free_task_list(grid_hybrid_task_list *ptr) {
  destroy_grid_context_cpu(ptr);
}

#endif
