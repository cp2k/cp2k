/* This file is part of tblite.
 * SPDX-Identifier: LGPL-3.0-or-later
 *
 * tblite is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tblite is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tblite.  If not, see <https://www.gnu.org/licenses/>.
**/

/**@dir include/tblite
 * @brief
 * Collection of public API bindings for the tblite library.
 */

/**@file tblite.h
 * @brief
 * Provides convenience functionality for working with the C API bindings for tblite
 * by including all relevant headers and defining general type generic macros for
 * working with the object handles.
 */

#pragma once

#include "tblite/error.h"
#include "tblite/container.h"
#include "tblite/context.h"
#include "tblite/double_dictionary.h"
#include "tblite/structure.h"
#include "tblite/calculator.h"
#include "tblite/solvation.h"
#include "tblite/result.h"
#include "tblite/table.h"
#include "tblite/param.h"
#include "tblite/version.h"

/// Generic macro to free an object handle.
///
/// @param ptr: Opaque object handle
#define tblite_delete(ptr) _Generic((ptr), \
                       tblite_error: tblite_delete_error, \
                   tblite_container: tblite_delete_container, \
                     tblite_context: tblite_delete_context, \
                   tblite_structure: tblite_delete_structure, \
           tblite_double_dictionary: tblite_delete_double_dictionary,\
                  tblite_calculator: tblite_delete_calculator, \
                      tblite_result: tblite_delete_result, \
                       tblite_table: tblite_delete_table, \
                       tblite_param: tblite_delete_param \
                                   )(&ptr)

/// Generic macro to check for error status
///
/// @param ptr: Opaque object handle
#define tblite_check(ptr) _Generic((ptr), \
                      tblite_error: tblite_check_error, \
                    tblite_context: tblite_check_context \
                                  )(ptr)
