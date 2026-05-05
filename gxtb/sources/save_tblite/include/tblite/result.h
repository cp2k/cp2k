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

/**@file tblite/result.h
 * @brief
 * Provides a storage container #tblite_result for all persistent quantities produced
 * in a calculation by #tblite_get_singlepoint. The data stored in #tblite_result can
 * be used as restart data for subsequent calculations.
 *
 * The individual entries of the #tblite_result container can be queried and retrieved.
 * For some quantities, like the count of the shells and orbitals duplicated queries
 * are implemented to allow the usage of the #tblite_result container without requiring
 * access to the #tblite_calculator which produced it.
 */

#pragma once

#include "tblite/macros.h"
#include "tblite/error.h"
#include "tblite/double_dictionary.h"

/// Container to for storing and handling calculation results
typedef struct _tblite_result* tblite_result;

/// Create new result container
///
/// @return New result container
TBLITE_API_ENTRY tblite_result TBLITE_API_CALL
tblite_new_result(void);

/// Create new result container from existing container
///
/// @param res: Existing result container
/// @return Copy of the results in a new container
TBLITE_API_ENTRY tblite_result TBLITE_API_CALL
tblite_copy_result(tblite_result res);

/// Delete a calculation environment object
///
/// @param res: Result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_result(tblite_result* res);

/// Retrieve number of atoms from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param nat: Number of atoms
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_number_of_atoms(tblite_error error,
                                  tblite_result res,
                                  int* nat);

/// Retrieve number of shells from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param nsh: Number of shells
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_number_of_shells(tblite_error error,
                                   tblite_result res,
                                   int* nsh);

/// Retrieve number of orbitals from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param nao: Number of orbitals
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_number_of_orbitals(tblite_error error,
                                     tblite_result res,
                                     int* nao);

/// Retrieve number of spins from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param nspin: Number of spins
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_number_of_spins(tblite_error error,
                                  tblite_result res,
                                  int* nspin);

/// Retrieve energy from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param energy: Total energy
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_energy(tblite_error error,
                         tblite_result res,
                         double* energy);

/// Retrieve atom-resolved energies from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param energies: Atom-resolved energies, shape [nat]
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_energies(tblite_error error,
                           tblite_result res,
                           double* energies);

/// Retrieve gradient from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param gradient: Cartesian gradient, shape [nat][3]
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_gradient(tblite_error error,
                           tblite_result res,
                           double* gradient);

/// Retrieve virial from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param sigma: Strain derivatives, shape [3][3]
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_virial(tblite_error error,
                         tblite_result res,
                         double* sigma);

/// Retrieve atomic charges from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param charges: Atomic partial charges, shape [nat]
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_charges(tblite_error error,
                          tblite_result res,
                          double* charges);

/// Retrieve bond orders from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param mbo: Bond orders, shape [nat][nat]
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_bond_orders(tblite_error error,
                              tblite_result res,
                              double* mbo);

/// Retrieve dipole moment from result container (order x, y, z)
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param dipole: Total dipole moment, shape [3]
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_dipole(tblite_error error,
                         tblite_result res,
                         double dipole[3]);

/// Retrieve traceless quadrupole moment from result container (packed xx, xy, yy, xz, yz, zz)
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param quadrupole: Total traceless quadrupole moment, shape [6]
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_quadrupole(tblite_error error,
                             tblite_result res,
                             double quadrupole[6]);

/// Retrieve orbital energies from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param emo: Eigenvalues for each orbital
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_orbital_energies(tblite_error error,
                                   tblite_result res,
                                   double* emo);

/// Retrieve orbital occupations from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param occ: Occupation numbers for each orbital
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_orbital_occupations(tblite_error error,
                                      tblite_result res,
                                      double* occ);

/// Retrieve orbital coefficients from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param cmat: Orbital coefficient matrix, shape [nspin][nao][nao]
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_orbital_coefficients(tblite_error error,
                                       tblite_result res,
                                       double* cmat);

/// Retrieve density matrix from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param pmat: Density matrix, shape [nspin][nao][nao]
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_density_matrix(tblite_error error,
                                 tblite_result res,
                                 double* pmat);

/// Retrieve overlap matrix from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param smat: Overlap matrix, shape [nao][nao]
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_overlap_matrix(tblite_error error,
                                 tblite_result res,
                                 double* smat);

/// Retrieve Hamiltonian matrix from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param hmat: Hamiltonian matrix, shape [nao][nao]
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_result_hamiltonian_matrix(tblite_error error,
                                     tblite_result res,
                                     double* hmat);

/// Retrieve Hamiltonian matrix from result container
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param dict: Pointer to dictionary object
TBLITE_API_ENTRY tblite_double_dictionary TBLITE_API_CALL
tblite_get_post_processing_dict(tblite_error error,
                                tblite_result res);

/// Save wavefunction to file
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param filename: Name of the file to save the wavefunction to
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_save_result_wavefunction(tblite_error error,
                                tblite_result res,
                                const char* filename);

/// Load wavefunction from file
///
/// @param error: Handle for error messages
/// @param res: Result container
/// @param filename: Name of the file to load the wavefunction from
TBLITE_API_CALL void TBLITE_API_CALL
tblite_load_result_wavefunction(tblite_error error,
                                tblite_result res,
                                const char* filename);