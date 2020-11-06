/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifdef __COLLOCATE_GPU

#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/grid_common.h"
#include "../cpu/coefficients.h"
#include "../cpu/collocation_integration.h"
#include "../cpu/grid_prepare_pab_dgemm.h"
#include "../cpu/private_header.h"
#include "../cpu/tensor_local.h"
#include "../cpu/utils.h"
#include "../ref/grid_ref_task_list.h"
void initialize_grid_parameters_on_gpu_step1(void *const ctx, const int level);
double compute_coefficients(grid_context *const ctx,
                            struct collocation_integration_ *handler,
                            const _task *task, tensor *const pab,
                            tensor *const work, tensor *const subblock,
                            tensor *const pab_prep, int *const prev_block_num,
                            int *const prev_iset, int *const prev_jset,
                            double *rp);

void compute_collocation_gpu(pgf_list_gpu *handler);
void rotate_to_cartesian_harmonics(const grid_basis_set *ibasis,
                                   const grid_basis_set *jbasis,
                                   const int iatom, const int jatom,
                                   const int iset, const int jset,
                                   double *const block, tensor *work,
                                   tensor *pab);

void initialize_grid_parameters_on_gpu(collocation_integration *handler);

pgf_list_gpu *create_worker_list(const int number_of_workers,
                                 const int batch_size, int device_id,
                                 const bool apply_cutoff);

void grid_collocate(collocation_integration *const handler,
                    const bool use_ortho, const double zetp, const double rp[3],
                    const double radius);
void destroy_worker_list(pgf_list_gpu *const list);

void reset_list_gpu(pgf_list_gpu *const list) {
  cudaSetDevice(list->device_id);
  list->list_length = 0;
  list->coef_dynamic_alloc_size_gpu_ = 0;
}

void my_worker_is_running(pgf_list_gpu *const my_worker) {
  if (!my_worker)
    return;

  if (my_worker->running) {
    cudaSetDevice(my_worker->device_id);
    cudaStreamSynchronize(my_worker->stream);
    reset_list_gpu(my_worker);
    my_worker->running = false;
  }
}

void add_orbital_to_list(pgf_list_gpu *const list, const int lp,
                         const double rp[3], const double radius,
                         const double zetp, const tensor *const coef) {
  assert(list->batch_size > list->list_length);

  list->lmax_cpu_[list->list_length] = lp;
  list->rp_cpu_[list->list_length].x = rp[0];
  list->rp_cpu_[list->list_length].y = rp[1];
  list->rp_cpu_[list->list_length].z = rp[2];

  list->radius_cpu_[list->list_length] = radius;

  list->zeta_cpu_[list->list_length] = zetp;
  list->coef_offset_cpu_[0] = 0;
  if (list->list_length > 0) {
    list->coef_offset_cpu_[list->list_length] =
        list->coef_offset_cpu_[list->list_length - 1] +
        list->coef_previous_alloc_size_;
    list->coef_dynamic_alloc_size_gpu_ += coef->alloc_size_;
  } else {
    list->coef_dynamic_alloc_size_gpu_ = coef->alloc_size_;
  }
  list->coef_previous_alloc_size_ = coef->alloc_size_;

  if (list->coef_dynamic_alloc_size_gpu_ > list->coef_alloc_size_gpu_) {
    list->durty = true;
    list->coef_cpu_ =
        realloc(list->coef_cpu_,
                sizeof(double) * (list->coef_dynamic_alloc_size_gpu_ + 4096));
    list->coef_alloc_size_gpu_ = list->coef_dynamic_alloc_size_gpu_ + 4096;
  }

  if (list->lmax < lp)
    list->lmax = lp;
  memcpy(list->coef_cpu_ + list->coef_offset_cpu_[list->list_length],
         coef->data, sizeof(double) * coef->alloc_size_);
  list->list_length++;
}

pgf_list_gpu *create_worker_list(const int number_of_workers,
                                 const int batch_size, const int device_id,
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
    list[i].lmax_cpu_ = (int *)malloc(sizeof(int) * list->batch_size);
    list[i].coef_offset_cpu_ = (int *)malloc(sizeof(int) * list->batch_size);
    list[i].rp_cpu_ = (double3 *)malloc(sizeof(double3) * list->batch_size);
    list[i].radius_cpu_ = (double *)malloc(sizeof(double) * list->batch_size);
    list[i].zeta_cpu_ = (double *)malloc(sizeof(double) * list->batch_size);
    list[i].coef_cpu_ =
        (double *)malloc(sizeof(double) * list->batch_size * 8 * 8 * 8);
    list[i].coef_alloc_size_gpu_ = list->batch_size * 8 * 8 * 8;
    cudaSetDevice(list[i].device_id);
    cudaMalloc((void **)&list[i].radius_gpu_,
               sizeof(double) * list->batch_size);
    cudaMalloc((void **)&list[i].coef_offset_gpu_,
               sizeof(int) * list->batch_size);
    cudaMalloc((void **)&list[i].lmax_gpu_, sizeof(int) * list->batch_size);
    cudaMalloc((void **)&list[i].rp_gpu_, sizeof(double3) * list->batch_size);
    cudaMalloc((void **)&list[i].zeta_gpu_, sizeof(double) * list->batch_size);
    cudaMalloc((void **)&list[i].coef_gpu_,
               sizeof(double) * list->coef_alloc_size_gpu_);
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

inline void destroy_worker_list(pgf_list_gpu *const lst) {
  cudaSetDevice(lst->device_id);
  cudaFree(lst->coef_offset_gpu_);
  cudaFree(lst->radius_gpu_);
  cudaFree(lst->lmax_gpu_);
  cudaFree(lst->rp_gpu_);
  cudaFree(lst->zeta_gpu_);
  cudaFree(lst->coef_gpu_);

  cudaFree(lst->data_gpu_);

  cudaStreamDestroy(lst->stream);
  cublasDestroy(lst->blas_handle);
  cudaEventDestroy(lst->event);
  free(lst->lmax_cpu_);
  free(lst->coef_offset_cpu_);
  free(lst->rp_cpu_);
  free(lst->radius_cpu_);
  free(lst->zeta_cpu_);
  free(lst->coef_cpu_);
}

inline void release_gpu_resources(collocation_integration *handler) {
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
                                         const int func, const int level) {
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
    initialize_basis_vectors(handler, grid->dh, grid->dh_inv);

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

    tensor work, pab, pab_prep, subblock;

    // Allocate pab matrix for re-use across tasks.
    initialize_tensor_2(&pab, ctx->maxco, ctx->maxco);
    alloc_tensor(&pab);

    initialize_tensor_2(&work, ctx->maxco, ctx->maxco);
    alloc_tensor(&work);

    initialize_tensor_2(&pab_prep, ctx->maxco, ctx->maxco);
    alloc_tensor(&pab_prep);

    initialize_tensor_2(&subblock, ctx->maxco, ctx->maxco);
    alloc_tensor(&subblock);

    // Initialize variables to detect when a new subblock has to be fetched.
    int prev_block_num = -1, prev_iset = -1, prev_jset = -1;

#pragma omp for schedule(static)
    for (int itask = 0; itask < ctx->tasks_per_level[level]; itask++) {
      // Define some convenient aliases.
      const _task *task = &ctx->tasks[level][itask];
      if (task->level != level) {
        printf("level %d, %d\n", task->level, level);
        abort();
      }
      double rp[3];
      double zetp = compute_coefficients(ctx, handler, task, &pab, &work,
                                         &subblock, &pab_prep, &prev_block_num,
                                         &prev_iset, &prev_jset, rp);

      add_orbital_to_list(current_worker, handler->coef.size[2] - 1, rp,
                          task->radius, zetp, &handler->coef);

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

    free(subblock.data);
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

void grid_collocate_task_list_gpu(
    void *const ptr, const bool orthorhombic, const int func, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    double *grid[nlevels]) {
  assert(ptr != NULL);
  grid_context *const ctx = (grid_context *const)ptr;
  const int max_threads = omp_get_max_threads();

  /* need to set the grid parameters first i need the total size for allocating
   * scratch memory */
  for (int level = 0; level < ctx->nlevels; level++) {
    set_grid_parameters(ctx, level, npts_global[level], npts_local[level],
                        shift_local[level], border_width[level], dh[level],
                        dh_inv[level], grid[level]);
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

  for (int thread = 0; thread < max_threads; thread++) {
    initialize_worker_list_on_gpu(
        ctx->handler[thread], ctx->device_id[thread % ctx->number_of_devices],
        length_queue, /* basically the number of gaussian we can treat */
        2);           /* number of workers */
  }

  /* We can not dispatch the different grids over different GPU with this
   * implementation. */
  for (int level = 0; level < ctx->nlevels; level++) {
    initialize_grid_parameters_on_gpu_step1(ctx, level);
    collocate_one_grid_level_gpu(ctx, func, level);
  }

  /* data is *not* allocated when I create the ctx so need to set it up to null.
   */
#pragma omp parallel for
  for (int thread = 0; thread < max_threads; thread++)
    release_gpu_resources(ctx->handler[thread]);

  free(ctx->scratch);
  ctx->scratch = NULL;
}
#endif
