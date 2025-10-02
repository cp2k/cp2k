/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <stdlib.h>

// Since __DATA_DIR can be arbitrarily long we must not let the preprocessor
// expand that macro in Fortran code as it could exceed the line length limit.
static const char *data_dir = __DATA_DIR;

/*******************************************************************************
 * \brief Returns path of data directory if set, otherwise an empty string.
 * \author Ole Schuett
 ******************************************************************************/
const char *get_data_dir() {
  const char *overwrite = getenv("CP2K_DATA_DIR");
  if (overwrite != NULL) {
    return overwrite;
  }
  return data_dir;
}

// EOF
