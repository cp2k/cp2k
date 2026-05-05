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

/**@file tblite/table.h
 * @brief
 * Provides a representation of a generic table data structure.
 *
 * Used to mirror the data available in the #tblite_param object. It aims to provide a
 * programmatic accessible representation of the parametrization records.
 */

#pragma once

#include "tblite/macros.h"
#include "tblite/error.h"

/// Handle for holding a table data structure
///
/// The table can either own its data or reference another table. Table references
/// can be created from existing table data structures or by adding new table entries
/// into an existing table, which returns a reference to the newly created table.
typedef struct _tblite_table* tblite_table;

/// Create new data table object
///
/// @param table: Table object to reference in new table (optional)
/// @return: New table data structure
TBLITE_API_ENTRY tblite_table TBLITE_API_CALL
tblite_new_table(tblite_table* table);

/// Delete a data table object
///
/// @param table: Table object to be deleted
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_table(tblite_table* table);

/// Set floating point number to data table
///
/// @param error: Error handle
/// @param table: Table data structure
/// @param key: Key to set value at
/// @param value: Double value or array to add to table
/// @param n: Number of entries to add, 0 for adding scalars
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_set_double(tblite_error error,
                        tblite_table table,
                        char key[],
                        double* value,
                        int n);

/// Set integer number to data table (use int64_t rather than long)
///
/// @param error: Error handle
/// @param table: Table data structure
/// @param key: Key to set value at
/// @param value: Integer value or array to add to table
/// @param n: Number of entries to add, 0 for adding scalars
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_set_int64_t(tblite_error error,
                         tblite_table table,
                         char key[],
                         int64_t* value,
                         int n);

/// Set boolean value to data table
///
/// @param error: Error handle
/// @param table: Table data structure
/// @param key: Key to set value at
/// @param value: Boolean value or array to add to table
/// @param n: Number of entries to add, 0 for adding scalars
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_set_bool(tblite_error error,
                      tblite_table table,
                      char key[],
                      bool* value,
                      int n);

/// Set character string to data table
///
/// @param error: Error handle
/// @param table: Table data structure
/// @param key: Key to set value at
/// @param value: Character value or array to add to table
/// @param n: Number of entries to add, 0 for adding scalars
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_table_set_char(tblite_error error,
                      tblite_table table,
                      char key[],
                      char (* value)[],
                      int n);

/// Create new subtable in existing data table
///
/// @param error: Error handle
/// @param table: Table data structure
/// @param key: Key to add new subtable at
/// @return New table data structure
TBLITE_API_ENTRY tblite_table TBLITE_API_CALL
tblite_table_add_table(tblite_error error,
                       tblite_table table,
                       char key[]);

/*
 * Type generic macros
 */

/// Generic setter based on the type of the value
///
/// @param error: Error handle
/// @param table: Table data structure
/// @param key: Key to set value at
/// @param value: Value to set at key
#define tblite_table_set_value(error, table, key, value, ...) \
    _Generic((value), \
     double*: tblite_table_set_double, \
    int64_t*: tblite_table_set_int64_t, \
       bool*: tblite_table_set_bool, \
   char(*)[]: tblite_table_set_char \
            )((error), (table), (key), (value), __VA_ARGS__)
