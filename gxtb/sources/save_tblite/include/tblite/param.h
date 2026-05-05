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

/**@file tblite/param.h
 * @brief
 * Provides a representation of a parametrization of an xTB Hamiltonian which can be used
 * to instantiate a #tblite_calculator object.
 *
 * The parametrization data itself can be represented as a #tblite_table data structure,
 * which provides the possibility to customize the parametrization programmatically.
 */

#pragma once

#include "tblite/macros.h"
#include "tblite/error.h"
#include "tblite/table.h"

/// Parametrization records
typedef struct _tblite_param* tblite_param;

/// Create new parametrization records object
///
/// @return Parametrization records
TBLITE_API_ENTRY tblite_param TBLITE_API_CALL
tblite_new_param(void);

/// Delete a parametrization records object
///
/// @param param: Parametrization records
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_param(tblite_param* param);

/// Load parametrization records from data table
///
/// @param error: Handle for error messages
/// @param param: Parametrization records
/// @param table: Table data structure
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_load_param(tblite_error error,
                  tblite_param param,
                  tblite_table table);

/// Dump parametrization records to data table
///
/// @param error: Handle for error messages
/// @param param: Parametrization records
/// @param table: Table data structure
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_dump_param(tblite_error error,
                  tblite_param param,
                  tblite_table table);

/// Export GFN2-xTB parametrization records
///
/// @param error: Handle for error messages
/// @param param: Parametrization records
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_export_gfn2_param(tblite_error error,
                         tblite_param param);

/// Export GFN1-xTB parametrization records
///
/// @param error: Handle for error messages
/// @param param: Parametrization records
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_export_gfn1_param(tblite_error error,
                         tblite_param param);

/// Export IPEA1-xTB parametrization records
///
/// @param error: Handle for error messages
/// @param param: Parametrization records
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_export_ipea1_param(tblite_error error,
                          tblite_param param);
