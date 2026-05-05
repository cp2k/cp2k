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

/**@file tblite/error.h
 * @brief
 * Provides a light-weight error handler for communicating with the library.
 *
 * The library provides two different kinds handlers, a light error handle type
 * #tblite_error is defined here.
 *
 * The error handle is used in the context of simple tasks and requires
 * only small overhead to construct and use.
 * It is mainly used in the context of retrieving data or building structures.
 */

#pragma once

#include "tblite/macros.h"

/// Error instance
typedef struct _tblite_error* tblite_error;

/// Create new error handle object
TBLITE_API_ENTRY tblite_error TBLITE_API_CALL
tblite_new_error(void);

/// Delete an error handle object
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_error(tblite_error* /* error */);

/// Check error handle status
TBLITE_API_ENTRY int TBLITE_API_CALL
tblite_check_error(tblite_error /* error */);

/// Clear error handle status
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_clear_error(tblite_error /* error */);

/// Get error message from error handle
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_error(tblite_error /* error */,
                 char* /* buffer */,
                 const int* /* buffersize */);

/// Set error message to error handle
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_set_error(tblite_error /* error */,
                 char* /* buffer */,
                 const int* /* buffersize */);
