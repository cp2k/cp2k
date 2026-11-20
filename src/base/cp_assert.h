/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef CP_ASSERT_H
#define CP_ASSERT_H

#include <stdio.h>
#include <stdlib.h>

static inline void cp_assert_fail(const char *condition, const char *file,
                                  int line, const char *function) {
  fprintf(stderr, "%s:%d: %s: assertion failed: %s\n", file, line, function,
          condition);
  abort();
}

/* Keep CP2K's runtime assertions active in all build configurations. */
#define cp_assert(condition)                                                   \
  do {                                                                         \
    if (!(condition)) {                                                        \
      cp_assert_fail(#condition, __FILE__, __LINE__, __func__);                \
    }                                                                          \
  } while (0)

#endif
