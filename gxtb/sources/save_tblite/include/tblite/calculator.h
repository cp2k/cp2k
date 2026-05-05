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

/**@file tblite/calculator.h
 * @brief
 * Provides a single point calculator for performing extended tight-binding computations.
 *
 * Provides a parametrized type #tblite_calculator storing the method parameters required
 * for the evaluation of single point calculations. The calculator is parametrized for
 * a method as well as for the molecular structure data #tblite_structure it was created
 * with. Changes in the #tblite_structure required to reinstantiate the structure data
 * also require to reinstantiate the #tblite_calculator object.
 */

#pragma once

#include "tblite/macros.h"
#include "tblite/context.h"
#include "tblite/structure.h"
#include "tblite/result.h"
#include "tblite/param.h"

/// Possible initial guess of the wavefunction.
typedef enum {
   /// Use superposition of atomic densities as guess.
   TBLITE_GUESS_SAD = 0,
   /// Use partial charges obtained by electronegativity equilibration as guess.
   TBLITE_GUESS_EEQ = 1,
} tblite_guess;

/// Single point calculator
typedef struct _tblite_calculator* tblite_calculator;

/// Construct calculator with GFN2-xTB parametrisation loaded
///
/// @param ctx: Context handle
/// @param mol: Molecular structure data
/// @return New calculator instance
TBLITE_API_ENTRY tblite_calculator TBLITE_API_CALL
tblite_new_gfn2_calculator(tblite_context ctx,
                           tblite_structure mol);

/// Construct calculator with GFN1-xTB parametrisation loaded
///
/// @param ctx: Context handle
/// @param mol: Molecular structure data
/// @return New calculator instance
TBLITE_API_ENTRY tblite_calculator TBLITE_API_CALL
tblite_new_gfn1_calculator(tblite_context ctx,
                           tblite_structure mol);

/// Construct calculator with IPEA1-xTB parametrisation loaded
///
/// @param ctx: Context handle
/// @param mol: Molecular structure data
/// @return New calculator instance
TBLITE_API_ENTRY tblite_calculator TBLITE_API_CALL
tblite_new_ipea1_calculator(tblite_context ctx,
                            tblite_structure mol);

/// Construct calculator from parametrization records
///
/// @param ctx: Context handle
/// @param mol: Molecular structure data
/// @return New calculator instance
/// @param param: Parametrization records
TBLITE_API_ENTRY tblite_calculator TBLITE_API_CALL
tblite_new_xtb_calculator(tblite_context ctx,
                          tblite_structure mol,
                          tblite_param param);

/// Delete calculator
///
/// @param calc: Calculator instance
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_calculator(tblite_calculator* calc);

/// Set calculation accuracy for the calculator object
///
/// @param ctx: Context handle
/// @param calc: Calculator instance
/// @param acc: Accuracy value for numerical thresholds
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_set_calculator_accuracy(tblite_context ctx,
                               tblite_calculator calc,
                               double acc);

/// Set maximum number of allowed iterations in calculator object
///
/// @param ctx: Context handle
/// @param calc: Calculator instance
/// @param max_iter: Maximum number of allowed iterations
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_set_calculator_max_iter(tblite_context ctx,
                               tblite_calculator calc,
                               int max_iter);

/// Set damping parameter for mixer in calculator object
///
/// @param ctx: Context handle
/// @param calc: Calculator instance
/// @param damping: Value for mixer damping
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_set_calculator_mixer_damping(tblite_context ctx,
                                    tblite_calculator calc,
                                    double damping);

/// Set initial guess for creating new wavefunction objects
///
/// @param ctx: Context handle
/// @param calc: Calculator instance
/// @param guess: Guess for initializing the wavefunction
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_set_calculator_guess(tblite_context ctx,
                            tblite_calculator calc,
                            tblite_guess guess);

/// Set electronic temperature for the calculator object (in Hartree)
///
/// @param ctx: Context handle
/// @param calc: Calculator instance
/// @param etemp: Electronic temperature in Hartree
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_set_calculator_temperature(tblite_context ctx,
                                  tblite_calculator calc,
                                  double etemp);

/// Set the flag in the calculator to retain the integral matrices
///
/// @param ctx: Context handle
/// @param calc: Calculator instance
/// @param save_integrals: Flag to enable storing of integrals
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_set_calculator_save_integrals(tblite_context ctx,
                                     tblite_calculator calc,
                                     int save_integrals);

/// Query calculator for the number of shells
///
/// @param ctx: Context handle
/// @param calc: Calculator instance
/// @param nsh: Number of shells
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_calculator_shell_count(tblite_context ctx,
                                  tblite_calculator calc,
                                  int* nsh);

/// Query calculator for index mapping from shells to atomic centers
///
/// @param ctx: Context handle
/// @param calc: Calculator instance
/// @param sh2at: Index mapping from shells to atomic centers
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_calculator_shell_map(tblite_context ctx,
                                tblite_calculator calc,
                                int* sh2at);

/// Query calculator for angular momenta of shells
///
/// @param ctx: Context handle
/// @param calc: Calculator instance
/// @param am: Angular momenta of each shell
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_calculator_angular_momenta(tblite_context ctx,
                                      tblite_calculator calc,
                                      int* am);

/// Query calculator for the number of orbitals
///
/// @param ctx: Context handle
/// @param calc: Calculator instance
/// @param nao: Number of atomic orbitals
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_calculator_orbital_count(tblite_context ctx,
                                    tblite_calculator calc,
                                    int* nao);

/// Query calculator for index mapping from atomic orbitals to shells
///
/// @param ctx: Context handle
/// @param calc: Calculator instance
/// @param ao2sh: Index mapping from atomic orbitals to shells
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_calculator_orbital_map(tblite_context ctx,
                                  tblite_calculator calc,
                                  int* ao2sh);

/// Perform single point calculation
///
/// @param ctx: Context handle
/// @param mol: Molecular structure data
/// @param calc: Calculator instance
/// @param res: Result container
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_singlepoint(tblite_context ctx,
                       tblite_structure mol,
                       tblite_calculator calc,
                       tblite_result res);

/// Push Back new conatiner to post processing construct
///
/// @param post_proc: Post Processing instance 
/// @param charptr: String of the post processing desired
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_push_back_post_processing_str(tblite_context ctx,
                                     tblite_calculator calc, 
                                     char* charptr);

/// Push Back new conatiner to post processing construct
///
/// @param post_proc: Post Processing instance 
/// @param param: Param instance containing post processing information
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_push_back_post_processing_param(tblite_context ctx,
                                       tblite_calculator calc,
                                       tblite_param param);