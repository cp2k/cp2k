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

/**@file tblite/version.h
 * @brief
 * Provides access to the version, compatibility and features exported by this API.
 */

#pragma once

#include "tblite/macros.h"

/// Retrieve version of library used
///
/// @return Compact version number in the format 10000 * major + 100 * minor + patch
TBLITE_API_ENTRY int TBLITE_API_CALL
tblite_get_version(void);
