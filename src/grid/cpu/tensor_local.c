/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include "tensor_local.h"
#include "../common/grid_common.h"
#include "utils.h"

size_t realloc_tensor(tensor *t) {
  assert(t);

  if (t->alloc_size_ == 0) {
    /* there is a mistake somewhere. We can not have t->old_alloc_size_ != 0 and
     * no allocation */
    abort();
  }

  if ((t->old_alloc_size_ >= t->alloc_size_) && (t->data != NULL))
    return t->alloc_size_;

  if ((t->old_alloc_size_ < t->alloc_size_) && (t->data != NULL)) {
    free(t->data);
  }

  t->data = NULL;

  if (t->data == NULL) {
    t->data = memalign(4096, sizeof(double) * t->alloc_size_);
    if (!t->data)
      abort();
    t->old_alloc_size_ = t->alloc_size_;
  }

  return t->alloc_size_;
}

void alloc_tensor(tensor *t) {
  if (t == NULL) {
    abort();
  }

  t->data = memalign(4096, sizeof(double) * t->alloc_size_);
  if (!t->data)
    abort();
  t->old_alloc_size_ = t->alloc_size_;
}

void tensor_copy(tensor *const b, const tensor *const a) {
  memcpy(b, a, sizeof(tensor));
}

void compute_block_dimensions(const int *const grid_size, int *const blockDim) {
  int block_size_test[8] = {2, 3, 4, 5, 6, 7, 8, 10};
  bool block_divided[8];

  blockDim[0] = -1;
  blockDim[1] = -1;
  blockDim[2] = -1;

  for (int d = 0; d < 3; d++) {
    for (int s = 0; s < 8; s++)
      block_divided[s] = (grid_size[d] % block_size_test[s] == 0);

    if (block_divided[7]) {
      blockDim[d] = 10;
      continue;
    }

    if (block_divided[6]) {
      blockDim[d] = 8;
      continue;
    }

    if (block_divided[5]) {
      blockDim[d] = 7;
      continue;
    }

    if (block_divided[4]) {
      blockDim[d] = 6;
      continue;
    }

    if (block_divided[3]) {
      blockDim[d] = 5;
      continue;
    }

    if (block_divided[2]) {
      blockDim[d] = 4;
      continue;
    }

    if (block_divided[1]) {
      blockDim[d] = 3;
      continue;
    }

    if (block_divided[0]) {
      blockDim[d] = 2;
      continue;
    }
  }
}

void decompose_grid_to_blocked_grid(const tensor *gr,
                                    struct tensor_ *block_grid) {
  if ((gr == NULL) || (block_grid == NULL)) {
    abort();
  }

  tensor tmp;
  initialize_tensor_3(&tmp, block_grid->blockDim[0], block_grid->blockDim[1],
                      block_grid->blockDim[2]);

  int lower_corner[3], upper_corner[3];

  for (int z = 0; z < block_grid->size[0]; z++) {
    lower_corner[0] = z * block_grid->blockDim[0];
    upper_corner[0] = lower_corner[0] + imin(gr->size[0] - lower_corner[0],
                                             block_grid->blockDim[0]);
    for (int y = 0; y < block_grid->size[1]; y++) {
      lower_corner[1] = y * block_grid->blockDim[1];
      upper_corner[1] = lower_corner[1] + imin(gr->size[1] - lower_corner[1],
                                               block_grid->blockDim[1]);
      for (int x = 0; x < block_grid->size[2]; x++) {
        tmp.data = &idx4(block_grid[0], z, y, x, 0);
        lower_corner[2] = x * block_grid->blockDim[2];
        upper_corner[2] = lower_corner[2] + imin(gr->size[2] - lower_corner[2],
                                                 block_grid->blockDim[2]);

        extract_sub_grid(lower_corner, upper_corner, NULL,
                         gr, // original grid
                         &tmp);
      }
    }
  }
}

/* recompose the natural grid from the block decomposed grid. Result is copied
 * to the grid gr */

void recompose_grid_from_blocked_grid(const struct tensor_ *block_grid,
                                      tensor *gr) {
  tensor tmp;
  initialize_tensor_3(&tmp, block_grid->blockDim[0], block_grid->blockDim[1],
                      block_grid->blockDim[2]);
  int lower_corner[3], upper_corner[3];
  for (int z = 0; z < block_grid->size[0]; z++) {
    lower_corner[0] = z * block_grid->blockDim[0];
    upper_corner[0] = lower_corner[0] + imin(gr->size[0] - lower_corner[0],
                                             block_grid->blockDim[0]);
    for (int y = 0; y < block_grid->size[1]; y++) {
      lower_corner[1] = y * block_grid->blockDim[1];
      upper_corner[1] = lower_corner[1] + imin(gr->size[1] - lower_corner[1],
                                               block_grid->blockDim[1]);
      for (int x = 0; x < block_grid->size[2]; x++) {
        tmp.data = &idx4(block_grid[0], z, y, x, 0);
        lower_corner[2] = x * block_grid->blockDim[2];
        upper_corner[2] = lower_corner[2] + imin(gr->size[2] - lower_corner[2],
                                                 block_grid->blockDim[2]);

        const int sizex = upper_corner[2] - lower_corner[2];
        const int sizey = upper_corner[1] - lower_corner[1];
        const int sizez = upper_corner[0] - lower_corner[0];

        for (int z = 0; z < sizez; z++) {
          for (int y = 0; y < sizey; y++) {
            double *__restrict__ dst =
                &idx3(gr[0], lower_corner[0] + z, lower_corner[1] + y,
                      lower_corner[2]);
            double *__restrict__ src = &idx3(tmp, z, y, 0);
            for (int x = 0; x < sizex; x++) {
              dst[x] = src[x];
            }
          }
        }
      }
    }
  }
}

/* recompose the natural grid from the block decomposed grid and add the result
 * to the grid gr */
void add_blocked_tensor_to_tensor(const struct tensor_ *block_grid,
                                  tensor *gr) {
  tensor tmp;
  initialize_tensor_3(&tmp, block_grid->blockDim[0], block_grid->blockDim[1],
                      block_grid->blockDim[2]);
  int lower_corner[3], upper_corner[3];

  for (int z = 0; z < block_grid->size[0]; z++) {
    lower_corner[0] = z * block_grid->blockDim[0];
    upper_corner[0] = lower_corner[0] + imin(gr->size[0] - lower_corner[0],
                                             block_grid->blockDim[0]);
    for (int y = 0; y < block_grid->size[1]; y++) {
      lower_corner[1] = y * block_grid->blockDim[1];
      upper_corner[1] = lower_corner[1] + imin(gr->size[1] - lower_corner[1],
                                               block_grid->blockDim[1]);
      for (int x = 0; x < block_grid->size[2]; x++) {
        tmp.data = &idx4(block_grid[0], z, y, x, 0);
        lower_corner[2] = x * block_grid->blockDim[2];
        upper_corner[2] = lower_corner[2] + imin(gr->size[2] - lower_corner[2],
                                                 block_grid->blockDim[2]);

        const int sizex = upper_corner[2] - lower_corner[2];
        const int sizey = upper_corner[1] - lower_corner[1];
        const int sizez = upper_corner[0] - lower_corner[0];

        for (int z1 = 0; z1 < sizez; z1++) {
          for (int y1 = 0; y1 < sizey; y1++) {
            double *__restrict__ dst =
                &idx3(gr[0], lower_corner[0] + z1, lower_corner[1] + y1,
                      lower_corner[2]);
            double *__restrict__ src = &idx3(tmp, z1, y1, 0);
            for (int x1 = 0; x1 < sizex; x1++) {
              dst[x1] += src[x1];
            }
          }
        }
      }
    }
  }
}

/* recompose the natural grid from the block decomposed grid and add the result
 * to the grid gr. The blocked tensor coordinates are in the yxz format */
void add_transpose_blocked_tensor_to_tensor(const struct tensor_ *block_grid,
                                            tensor *gr) {
  tensor tmp;
  initialize_tensor_3(&tmp, block_grid->blockDim[0], block_grid->blockDim[1],
                      block_grid->blockDim[2]);
  int lower_corner[3], upper_corner[3];

  for (int y = 0; y < block_grid->size[1]; y++) {
    lower_corner[1] = y * block_grid->blockDim[1];
    upper_corner[1] = lower_corner[1] + imin(gr->size[1] - lower_corner[1],
                                             block_grid->blockDim[1]);
    for (int x = 0; x < block_grid->size[2]; x++) {
      lower_corner[2] = x * block_grid->blockDim[2];
      upper_corner[2] = lower_corner[2] + imin(gr->size[2] - lower_corner[2],
                                               block_grid->blockDim[2]);
      for (int z = 0; z < block_grid->size[0]; z++) {
        lower_corner[0] = z * block_grid->blockDim[0];
        upper_corner[0] = lower_corner[0] + imin(gr->size[0] - lower_corner[0],
                                                 block_grid->blockDim[0]);
        tmp.data = &idx4(block_grid[0], y, x, z, 0);

        const int sizex = upper_corner[2] - lower_corner[2];
        const int sizey = upper_corner[1] - lower_corner[1];
        const int sizez = upper_corner[0] - lower_corner[0];

        for (int z1 = 0; z1 < sizez; z1++) {
          for (int y1 = 0; y1 < sizey; y1++) {
            double *__restrict__ dst =
                &idx3(gr[0], lower_corner[0] + z1, lower_corner[1] + y1,
                      lower_corner[2]);
            double *__restrict__ src = &idx3(tmp, z1, y1, 0);
            for (int x1 = 0; x1 < sizex; x1++) {
              dst[x1] += src[x1];
            }
          }
        }
      }
    }
  }
}

void compute_block_boundaries(const int *blockDim, const int *lb_grid,
                              const int *grid_size,
                              const int *blocked_grid_size, const int *period,
                              const int *cube_center, const int *cube_size,
                              const int *lower_boundaries_cube,
                              int *lower_block_corner, int *upper_block_corner,
                              int *pol_offsets) {
  int position[3];
  return_cube_position(grid_size, lb_grid, cube_center, lower_boundaries_cube,
                       period, position);
  pol_offsets[0] = 0;
  pol_offsets[1] = 0;
  pol_offsets[2] = 0;

  for (int axis = 0; axis < 3; axis++) {
    int tmp = position[axis];
    int blockidx = tmp / blockDim[axis];
    lower_block_corner[axis] = blockidx;
    pol_offsets[axis] = tmp - blockidx * blockDim[axis];
    tmp = position[axis] + cube_size[axis];
    if ((grid_size[axis] != period[axis]) && (tmp > grid_size[axis])) {
      upper_block_corner[axis] = blocked_grid_size[axis];
    } else {
      upper_block_corner[axis] =
          tmp / blockDim[axis] + ((tmp - blockidx * blockDim[axis]) != 0);
    }
  }

  return;
}

/* initialize a tensor structure for a tensor of dimension dim <= 4 */

void initialize_tensor_blocked(struct tensor_ *a, const int dim,
                               const int *const sizes,
                               const int *const blockDim) {
  assert(a != NULL);

  a->block = (void *)malloc(sizeof(struct tensor_));
  memset(a->block, 0, sizeof(struct tensor_));

  switch (dim) {
  case 4:
    initialize_tensor_4((struct tensor_ *)a->block, blockDim[0], blockDim[1],
                        blockDim[2], blockDim[3]);
    break;
  case 3:
    initialize_tensor_3((struct tensor_ *)a->block, blockDim[0], blockDim[1],
                        blockDim[2]);
    break;
  case 2:
    initialize_tensor_2((struct tensor_ *)a->block, blockDim[0], blockDim[1]);
    break;
  default:
    printf("We should not be here");
    assert(0);
    break;
  }

  assert(a->block->alloc_size_ != 0);
  a->dim_ = dim + 1;

  for (int d = 0; d < dim; d++)
    a->blockDim[d] = blockDim[d];

  for (int d = 0; d < a->dim_ - 1; d++) {
    a->size[d] = sizes[d] / a->blockDim[d] + (sizes[d] % a->blockDim[d] != 0);
    a->unblocked_size[d] = sizes[d];
  }

  a->size[dim] = a->block->alloc_size_;
  // we need proper alignment here. But can be done later
  /* a->ld_ = (sizes[a->dim_ - 1] / 32 + 1) * 32; */
  a->ld_ = a->block->alloc_size_;
  switch (a->dim_) {
  case 5: {
    a->offsets[0] = a->ld_ * a->size[1] * a->size[2] * a->size[3];
    a->offsets[1] = a->ld_ * a->size[1] * a->size[2];
    a->offsets[2] = a->ld_ * a->size[2];
    a->offsets[3] = a->ld_;
    break;
  }
  case 4: {
    a->offsets[0] = a->ld_ * a->size[1] * a->size[2];
    a->offsets[1] = a->ld_ * a->size[2];
    a->offsets[2] = a->ld_;
    break;
  }
  case 3: {
    a->offsets[0] = a->ld_ * a->size[1];
    a->offsets[1] = a->ld_;
  } break;
  case 2: { // matrix case
    a->offsets[0] = a->ld_;
  } break;
  case 1:
    break;
  }

  a->alloc_size_ = a->offsets[0] * a->size[0];
  a->blocked_decomposition = true;
  assert(a->alloc_size_ != 0);
  return;
}
