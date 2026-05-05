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

/**@file tblite/solvation.h
 * @brief
 * Provides solvation specific interaction containers which can be added
 * to a #tblite_calculator.
 */

#pragma once

#include "tblite/macros.h"
#include "tblite/structure.h"
#include "tblite/context.h"
#include "tblite/container.h"

/// Reference state enumerator
enum tblite_ref_solvation_state{
    tblite_state_gsolv = 1,
    tblite_state_bar1mol = 2,
    tblite_state_reference = 3,
};

enum tblite_born_kernel{
    tblite_born_still = 1,
    tblite_born_p16 = 2,
};

enum tblite_solvation_param{
    tblite_solvation_gbe = 10,
    tblite_solvation_alpb_gfn1 = 11,
    tblite_solvation_alpb_gfn2 = 12,
    tblite_solvation_gb = 20,
    tblite_solvation_gbsa_gfn1 = 21,
    tblite_solvation_gbsa_gfn2 = 22,
};

/// Create new CPCM implicit solvation container using internal parameters
///
/// @param error: Error handle
/// @param mol: Molecular structure data
/// @param eps: epsilon value for solvent
/// @return New interaction container
TBLITE_API_ENTRY tblite_container TBLITE_API_CALL
tblite_new_cpcm_solvation_epsilon(tblite_error error,
                                  tblite_structure mol,
                                  double eps);

/// Create new ALPB implicit solvation container using internal parameters
///
/// @param error: Error handle
/// @param mol: Molecular structure data
/// @param eps: epsilon value of solvent
/// @param version: type of solvation model to use (ALPB=1, GBSA=0)
/// @param born: type of Born kernel to use (Still=1, P16=2)
/// @return New interaction container
TBLITE_API_ENTRY tblite_container TBLITE_API_CALL
tblite_new_gb_solvation_epsilon(tblite_error error,
                                tblite_structure mol,
                                double eps,
                                int version,
                                int born);

/// Create new ALPB implicit solvation container using internal parameters
///
/// @param error: Error handle
/// @param mol: Molecular structure data
/// @param solvent: Solvent name to be described
/// @param version: Parametrization (ALPB(GFN1)=11, ALPB(GFN2)=12, GBSA(GFN1)=21, GBSA(GFN2)=22)
/// @param refstate: Reference state of solution (gsolv=1, bar1mol=2, reference=3)
/// @return New interaction container
TBLITE_API_ENTRY tblite_container TBLITE_API_CALL
tblite_new_alpb_solvation_solvent(tblite_error error,
                                  tblite_structure mol,
                                  char* solvent,
                                  int version,
                                  int refstate);