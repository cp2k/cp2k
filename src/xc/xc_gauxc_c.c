/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <stddef.h>
#include <stdio.h>

void print_c_string(const char *c_message) {
  if (c_message != NULL) {
    printf("%s", c_message);
  }
}
