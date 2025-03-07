/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#include "grid_dgemm_tensor_local.h"
#include "../common/grid_common.h"
#include "grid_dgemm_utils.h"

size_t realloc_tensor(tensor *t) {
  assert(t != NULL);

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
    t->data = malloc(sizeof(double) * t->alloc_size_);
    assert(t->data != NULL);
    t->old_alloc_size_ = t->alloc_size_;
  }

  return t->alloc_size_;
}

void alloc_tensor(tensor *t) {
  assert(t != NULL);

  t->data = malloc(sizeof(double) * t->alloc_size_);
  assert(t->data != NULL);
  t->old_alloc_size_ = t->alloc_size_;
}

void tensor_copy(tensor *const b, const tensor *const a) {
  memcpy(b, a, sizeof(tensor));
}
