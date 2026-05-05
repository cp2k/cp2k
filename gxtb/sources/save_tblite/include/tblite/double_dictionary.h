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

/**@file tblite/double_dictionary.h
 * @brief
 * Provides access to the double dictionary class in tblite, which gathers post processing values.
 *
 * Provides access to a double dictionary class, the entries of the dictionary can be retrieved by index.
 * The number of indices can be requested, making asimple itareation over thee possible
 */
#pragma once
#include "tblite/macros.h"

// Double Dictionary Container
typedef struct _tblite_double_dictionary* tblite_double_dictionary;

/// Get number of entries
///
/// @param error: Error handler
/// @param dict: Double dictionary instance
TBLITE_API_ENTRY int TBLITE_API_CALL
tblite_get_n_entries_dict(tblite_error error,
                          tblite_double_dictionary dict);

/// Get the array associated with an entry by index, together with its dimensions
///
/// @param error: Error handler
/// @param dict: Double dictionary instance
/// @param index: Index of the entry for which to retrieve the label
/// @param array: Array associated to the entry addressed by the index
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_array_entry_index(tblite_error error,
                            tblite_double_dictionary dict,
                            const int* index,
                            double* array);

/// Get the array associated with an entry by index, together with its dimensions
///
/// @param error: Error handler
/// @param dict: Double dictionary instance
/// @param index: Index of the entry for which to retrieve the label
/// @param dim1: 1st dimension of the associated tensor to the index
/// @param dim2: 2nd dimension of the associated tensor to the index
/// @param dim3: 3rd dimension of the associated tensor to the index
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_array_size_index(tblite_error error,
                            tblite_double_dictionary dict,
                            const int* index,
                            int* dim1,
                            int* dim2,
                            int* dim3);

/// Get the array associated with an entry by index, together with its dimensions
///
/// @param error: Error handler
/// @param dict: Double dictionary instance
/// @param label: label of the entry for which to retrieve the label
/// @param array: Array associated to the entry addressed by the label
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_array_entry_label(tblite_error error,
                            tblite_double_dictionary dict,
                            char* label,
                            double* array);

/// Get the array associated with an entry by label, together with its dimensions
///
/// @param error: Error handler
/// @param dict: Double dictionary instance
/// @param label: label of the entry for which to retrieve the label
/// @param dim1: 1st dimension of the associated tensor to the label
/// @param dim2: 2nd dimension of the associated tensor to the label
/// @param dim3: 3rd dimension of the associated tensor to the label
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_array_size_label(tblite_error error,
                            tblite_double_dictionary dict,
                            char* label,
                            int* dim1,
                            int* dim2,
                            int* dim3);

/// Get label of an entry by index
///
/// @param error: Error handler
/// @param dict: Double dictionary instance
/// @param index: Index of the entry for which to retrieve the label
/// @param label: Label which is retrieved
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_label_entry_index(tblite_error error,
                             tblite_double_dictionary dict,
                             const int* index,
                             char* label,
                             const int* buffersize);

/// Delete dictionary
///
/// @param dict: Double dictionary instance
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_double_dictionary(tblite_double_dictionary* dict);