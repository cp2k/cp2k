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

/**@file tblite/structure.h
 * @brief
 * The structure data #tblite_structure is used to represent the system of interest
 * in the library.
 *
 * It contains immutable system specific information like the
 * number of atoms, the unique atom groups and the boundary conditions as well as
 * mutable geometry data like cartesian coordinates and lattice parameters.
 *
 * To change immutable parameters of the #tblite_structure data the object should
 * be reinstantiated as well as all dependent objects, like #tblite_calculator or
 * #tblite_result instances.
 */

#pragma once

#include "tblite/macros.h"
#include "tblite/error.h"

/*
 * Molecular structure data class
**/

/// Molecular structure data class.
///
/// The structure data is used to represent the system of interest in the library.
/// It contains immutable system specific information like the number of atoms,
/// the unique atom groups and the boundary conditions as well as mutable geometry
/// data like cartesian coordinates and lattice parameters.
typedef struct _tblite_structure* tblite_structure;

/// Create new molecular structure data
///
/// @param error: Handle for error messages
/// @param natoms: Number of atoms
/// @param numbers: Atomic numbers for each atom, shape [natoms]
/// @param positions: Cartesian coordinates in Bohr for each atom, shape [natoms][3]
/// @param charge: Total charge of the system, (optional)
/// @param uhf: Number of unpaired electrons, (optional)
/// @param lattice: Lattice parameters in Bohr, shape [3][3], (optional)
/// @param periodic: Periodic dimensions, shape [3], (optional)
/// @return New molecular structure data
TBLITE_API_ENTRY tblite_structure TBLITE_API_CALL
tblite_new_structure(tblite_error error,
                     const int natoms,
                     const int* numbers,
                     const double* positions,
                     const double* charge,
                     const int* uhf,
                     const double* lattice,
                     const bool* periodic);

/// Delete molecular structure data
///
/// @param mol: Molecular structure data
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_structure(tblite_structure* mol);

/// Update coordinates and lattice parameters (quantities in Bohr)
///
/// @param error: Handle for error messages
/// @param mol: Molecular structure data
/// @param positions: Cartesian coordinates in Bohr for each atom, shape [natoms][3]
/// @param lattice: Lattice parameters in Bohr, shape [3][3], (optional)
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_update_structure_geometry(tblite_error error,
                                 tblite_structure mol,
                                 const double* positions,
                                 const double* lattice);

/// Update total charge in structure object
///
/// @param error: Handle for error messages
/// @param mol: Molecular structure data
/// @param charge: Total charge of the system
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_update_structure_charge(tblite_error error,
                               tblite_structure mol,
                               const double* charge);

/// Update number of unpaired electrons in structure object
///
/// @param error: Handle for error messages
/// @param mol: Molecular structure data
/// @param uhf: Number of unpaired electrons
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_update_structure_uhf(tblite_error error,
                            tblite_structure mol,
                            const int* uhf);
