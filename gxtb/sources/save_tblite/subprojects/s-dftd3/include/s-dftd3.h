/* This file is part of s-dftd3.
 * SPDX-Identifier: LGPL-3.0-or-later
 *
 * s-dftd3 is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * s-dftd3 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.
**/
#pragma once

#ifdef __cplusplus
#define SDFTD3_API_ENTRY extern "C"
#else
#define SDFTD3_API_ENTRY extern
#ifndef SDFTD3_CFFI
#include <stdbool.h>
#endif
#endif
#define SDFTD3_API_CALL
#define SDFTD3_API_SUFFIX__V_0_2
#define SDFTD3_API_SUFFIX__V_0_3
#define SDFTD3_API_SUFFIX__V_0_4
#define SDFTD3_API_SUFFIX__V_0_5
#define SDFTD3_API_SUFFIX__V_1_3

/// Error handle class
typedef struct _dftd3_error* dftd3_error;

/// Molecular structure data class
typedef struct _dftd3_structure* dftd3_structure;

/// Dispersion model class
typedef struct _dftd3_model* dftd3_model;

/// Counter-poisecorrection parameters class
typedef struct _dftd3_gcp* dftd3_gcp;

/// Damping parameter class
typedef struct _dftd3_param* dftd3_param;

/*
 * Type generic macro for convenience
**/

#ifndef __cplusplus
#define dftd3_delete(ptr) _Generic((ptr), \
                       dftd3_error: dftd3_delete_error, \
                   dftd3_structure: dftd3_delete_structure, \
                       dftd3_model: dftd3_delete_model, \
                       dftd3_param: dftd3_delete_param, \
                         dftd3_gcp: dftd3_delete_gcp \
                                  )(&ptr)
#endif

/*
 * Global API queries
**/

/// Obtain library version as major * 10000 + minor + 100 + patch
SDFTD3_API_ENTRY int SDFTD3_API_CALL
dftd3_get_version() SDFTD3_API_SUFFIX__V_0_2;

/*
 * Error handle class
**/

/// Create new error handle object
SDFTD3_API_ENTRY dftd3_error SDFTD3_API_CALL
dftd3_new_error() SDFTD3_API_SUFFIX__V_0_2;

/// Check error handle status
SDFTD3_API_ENTRY int SDFTD3_API_CALL
dftd3_check_error(dftd3_error /* error */) SDFTD3_API_SUFFIX__V_0_2;

/// Get error message from error handle
SDFTD3_API_ENTRY void SDFTD3_API_CALL
dftd3_get_error(dftd3_error /* error */,
                char* /* buffer */,
                const int* /* buffersize */) SDFTD3_API_SUFFIX__V_0_2;

/// Delete error handle object
SDFTD3_API_ENTRY void SDFTD3_API_CALL
dftd3_delete_error(dftd3_error* /* error */) SDFTD3_API_SUFFIX__V_0_2;

/*
 * Molecular structure data class
**/

/// Create new molecular structure data (quantities in Bohr)
SDFTD3_API_ENTRY dftd3_structure SDFTD3_API_CALL
dftd3_new_structure(dftd3_error /* error */,
                    const int /* natoms */,
                    const int* /* numbers [natoms] */,
                    const double* /* positions [natoms][3] */,
                    const double* /* lattice [3][3] */,
                    const bool* /* periodic [3] */) SDFTD3_API_SUFFIX__V_0_2;

/// Delete molecular structure data
SDFTD3_API_ENTRY void SDFTD3_API_CALL
dftd3_delete_structure(dftd3_structure* /* mol */) SDFTD3_API_SUFFIX__V_0_2;

/// Update coordinates and lattice parameters (quantities in Bohr)
SDFTD3_API_ENTRY void SDFTD3_API_CALL
dftd3_update_structure(dftd3_error /* error */,
                       dftd3_structure /* mol */,
                       const double* /* positions [natoms][3] */,
                       const double* /* lattice [3][3] */) SDFTD3_API_SUFFIX__V_0_2;

/*
 * Dispersion model class
**/

/// Create new D3 dispersion model
SDFTD3_API_ENTRY dftd3_model SDFTD3_API_CALL
dftd3_new_d3_model(dftd3_error /* error */,
                   dftd3_structure /* mol */) SDFTD3_API_SUFFIX__V_0_2;

/// Set realspace cutoffs (quantities in Bohr)
SDFTD3_API_ENTRY void SDFTD3_API_CALL
dftd3_set_model_realspace_cutoff(dftd3_error /* error */,
                                 dftd3_model /* model */,
                                 double /* disp2 */,
                                 double /* disp3 */,
                                 double /* cn */) SDFTD3_API_SUFFIX__V_0_5;

/// Delete dispersion model
SDFTD3_API_ENTRY void SDFTD3_API_CALL
dftd3_delete_model(dftd3_model* /* disp */) SDFTD3_API_SUFFIX__V_0_2;

/*
 * Damping parameter class
**/

/// Create new zero damping parameters
SDFTD3_API_ENTRY dftd3_param SDFTD3_API_CALL
dftd3_new_zero_damping(dftd3_error /* error */,
                       double /* s6 */,
                       double /* s8 */,
                       double /* s9 */,
                       double /* rs6 */,
                       double /* rs8 */,
                       double /* alp */) SDFTD3_API_SUFFIX__V_0_4;

/// Load zero damping parameters from internal storage
SDFTD3_API_ENTRY dftd3_param SDFTD3_API_CALL
dftd3_load_zero_damping(dftd3_error /* error */,
                        char* /* method */,
                        bool /* atm */) SDFTD3_API_SUFFIX__V_0_4;

/// Create new rational damping parameters
SDFTD3_API_ENTRY dftd3_param SDFTD3_API_CALL
dftd3_new_rational_damping(dftd3_error /* error */,
                           double /* s6 */,
                           double /* s8 */,
                           double /* s9 */,
                           double /* a1 */,
                           double /* a2 */,
                           double /* alp */) SDFTD3_API_SUFFIX__V_0_4;

/// Load rational damping parameters from internal storage
SDFTD3_API_ENTRY dftd3_param SDFTD3_API_CALL
dftd3_load_rational_damping(dftd3_error /* error */,
                            char* /* method */,
                            bool /* atm */) SDFTD3_API_SUFFIX__V_0_4;

/// Create new modified zero damping parameters
SDFTD3_API_ENTRY dftd3_param SDFTD3_API_CALL
dftd3_new_mzero_damping(dftd3_error /* error */,
                        double /* s6 */,
                        double /* s8 */,
                        double /* s9 */,
                        double /* rs6 */,
                        double /* rs8 */,
                        double /* alp */,
                        double /* bet */) SDFTD3_API_SUFFIX__V_0_4;

/// Load modified zero damping parameters from internal storage
SDFTD3_API_ENTRY dftd3_param SDFTD3_API_CALL
dftd3_load_mzero_damping(dftd3_error /* error */,
                         char* /* method */,
                         bool /* atm */) SDFTD3_API_SUFFIX__V_0_4;

/// Create new modified rational damping parameters
SDFTD3_API_ENTRY dftd3_param SDFTD3_API_CALL
dftd3_new_mrational_damping(dftd3_error /* error */,
                            double /* s6 */,
                            double /* s8 */,
                            double /* s9 */,
                            double /* a1 */,
                            double /* a2 */,
                            double /* alp */) SDFTD3_API_SUFFIX__V_0_4;

/// Load modified rational damping parameters from internal storage
SDFTD3_API_ENTRY dftd3_param SDFTD3_API_CALL
dftd3_load_mrational_damping(dftd3_error /* error */,
                             char* /* method */,
                             bool /* atm */) SDFTD3_API_SUFFIX__V_0_4;

/// Create new optimized power damping parameters
SDFTD3_API_ENTRY dftd3_param SDFTD3_API_CALL
dftd3_new_optimizedpower_damping(dftd3_error /* error */,
                                 double /* s6 */,
                                 double /* s8 */,
                                 double /* s9 */,
                                 double /* a1 */,
                                 double /* a2 */,
                                 double /* alp */,
                                 double /* bet */) SDFTD3_API_SUFFIX__V_0_5;

/// Load optimized power damping parameters from internal storage
SDFTD3_API_ENTRY dftd3_param SDFTD3_API_CALL
dftd3_load_optimizedpower_damping(dftd3_error /* error */,
                                  char* /* method */,
                                  bool /* atm */) SDFTD3_API_SUFFIX__V_0_5;

/// Create new CSO damping parameters
SDFTD3_API_ENTRY dftd3_param SDFTD3_API_CALL
dftd3_new_cso_damping(dftd3_error /* error */,
                      double /* s6 */,
                      double /* s9 */,
                      double /* a1 */,
                      double /* a2 */,
                      double /* a3 */,
                      double /* a4 */,
                      double /* alp */) SDFTD3_API_SUFFIX__V_1_3;

/// Load CSO damping parameters from internal storage
SDFTD3_API_ENTRY dftd3_param SDFTD3_API_CALL
dftd3_load_cso_damping(dftd3_error /* error */,
                       char* /* method */,
                       bool /* atm */) SDFTD3_API_SUFFIX__V_1_3;

/// Delete damping parameters
SDFTD3_API_ENTRY void SDFTD3_API_CALL
dftd3_delete_param(dftd3_param* /* param */) SDFTD3_API_SUFFIX__V_0_2;

/*
 * Counter-poise correction parameters
**/

/// Load geometric counter-poise parameters from internal storage
SDFTD3_API_ENTRY dftd3_gcp SDFTD3_API_CALL
dftd3_load_gcp_param(dftd3_error /* error */,
                     dftd3_structure /* mol */,
                     char* /* method */,
                     char* /* basis */) SDFTD3_API_SUFFIX__V_1_3;

/// Set realspace cutoffs (quantities in Bohr)
SDFTD3_API_ENTRY void SDFTD3_API_CALL
dftd3_set_gcp_realspace_cutoff(dftd3_error /* error */,
                               dftd3_gcp /* gcp */,
                               double /* bas */,
                               double /* srb */) SDFTD3_API_SUFFIX__V_1_3;

/// Delete counter-poise parameters
SDFTD3_API_ENTRY void SDFTD3_API_CALL
dftd3_delete_gcp(dftd3_gcp* /* gcp */) SDFTD3_API_SUFFIX__V_1_3;

/*
 * Perform dispersion calculations
**/

/// Evaluate the dispersion energy and its derivatives
SDFTD3_API_ENTRY void SDFTD3_API_CALL
dftd3_get_dispersion(dftd3_error /* error */,
                     dftd3_structure /* mol */,
                     dftd3_model /* disp */,
                     dftd3_param /* param */,
                     double* /* energy */,
                     double* /* gradient[n][3] */,
                     double* /* sigma[3][3] */) SDFTD3_API_SUFFIX__V_0_2;

/// Evaluate the pairwise representation of the dispersion energy
SDFTD3_API_ENTRY void SDFTD3_API_CALL
dftd3_get_pairwise_dispersion(dftd3_error /* error */,
                              dftd3_structure /* mol */,
                              dftd3_model /* disp */,
                              dftd3_param /* param */,
                              double* /* pair_energy2[n][n] */,
                              double* /* pair_energy3[n][n] */) SDFTD3_API_SUFFIX__V_0_5;

/*
 * Perform geometric counterpoise calculations
**/

/// Evaluate the dispersion energy and its derivatives
SDFTD3_API_ENTRY void SDFTD3_API_CALL
dftd3_get_counterpoise(dftd3_error /* error */,
                       dftd3_structure /* mol */,
                       dftd3_gcp /* gcp */,
                       double* /* energy */,
                       double* /* gradient[n][3] */,
                       double* /* sigma[3][3] */) SDFTD3_API_SUFFIX__V_1_3;