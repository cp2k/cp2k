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

/**@file tblite/context.h
 * @brief
 * Provides a persistent configuration object to modify the behaviour of a calculation.
 * Acts as an error handler.
 *
 * The environment context #tblite_context can be considered a persistent setup for all
 * calculations performed with the library, it is usually used together with calculator
 * objects #tblite_calculator.
 * While the error handle can only contain a single error, multiple errors can be
 * accumulated in a context object, which allows storing more complex error information
 * like they can occur in an actual calculation.
 */

#pragma once

#include "tblite/macros.h"
#include "tblite/error.h"

/// Context manager for the library usage
typedef struct _tblite_context* tblite_context;

/// Define callback function for use in custom logger
typedef void (*tblite_logger_callback)(tblite_error, char*, int, void*);

#ifdef TBLITE_CFFI
extern "Python" void TBLITE_API_CALL
logger_callback(tblite_error error, char* msg, int len, void* user_data);
#endif

/// Create new calculation environment object
///
/// @return New context handle
TBLITE_API_ENTRY tblite_context TBLITE_API_CALL
tblite_new_context(void);

/// Delete a calculation environment object
///
/// @param ctx: Context handle
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_context(tblite_context* ctx);

/// Check calculation environment status
///
/// @param ctx: Context handle
/// @return Current status, 0: okay, 1: error
TBLITE_API_ENTRY int TBLITE_API_CALL
tblite_check_context(tblite_context ctx);

/// Get error message from calculation environment
///
/// @param ctx: Context handle
/// @param buffer: Character buffer for writing error message to
/// @param buffersize: Maximum length of buffer (optional)
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_context_error(tblite_context ctx,
                         char* buffer,
                         const int* buffersize);

/// Set custom logger function
///
/// @param ctx: Context handle
/// @param callback: Procedure pointer implementing logger
/// @param userdata: Passthrough data pointer
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_set_context_logger(tblite_context ctx,
                          tblite_logger_callback callback,
                          void* userdata);

/// Enable colorful output
///
/// @param ctx: Context handle
/// @param color: Set color support, 0: disabled, 1: enabled
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_set_context_color(tblite_context ctx,
                         int color);

/// Set verbosity level of printout
///
/// @param ctx: Context handle
/// @param verbosity: Printout verbosity
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_set_context_verbosity(tblite_context ctx,
                             int verbosity);
