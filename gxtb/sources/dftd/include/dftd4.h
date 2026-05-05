/* This file is part of dftd4.
 * SPDX-Identifier: LGPL-3.0-or-later
 *
 * dftd4 is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * dftd4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with dftd4.  If not, see <https://www.gnu.org/licenses/>.
**/
#pragma once

#ifdef __cplusplus
#define DFTD4_API_ENTRY extern "C"
#else
#define DFTD4_API_ENTRY extern
#ifndef DFTD4_CFFI
#include <stdbool.h>
#endif
#endif
#define DFTD4_API_CALL
#define DFTD4_API_SUFFIX__V_3_0
#define DFTD4_API_SUFFIX__V_3_1
#define DFTD4_API_SUFFIX__V_3_2
#define DFTD4_API_SUFFIX__V_3_3
#define DFTD4_API_SUFFIX__V_3_4
#define DFTD4_API_SUFFIX__V_3_5
#define DFTD4_API_SUFFIX__V_4_0
#define DFTD4_API_SUFFIX__V_5_0

/// Error handle class
typedef struct _dftd4_error* dftd4_error;

/// Molecular structure data class
typedef struct _dftd4_structure* dftd4_structure;

/// Dispersion model class
typedef struct _dftd4_model* dftd4_model;

/// Dispersion model enumerator
enum dftd4_dispersion_model{
    dftd4_model_d4 = 1,
    dftd4_model_d4s = 2,
};

/// Damping function class
typedef struct _dftd4_damping* dftd4_damping;

/// Two-body damping function enumerator
enum dftd4_damping_twobody{
    dftd4_damping_twobody_rational = 1,
    dftd4_damping_twobody_screened = 2,
    dftd4_damping_twobody_zero = 3,
    dftd4_damping_twobody_mzero = 4,
    dftd4_damping_twobody_optpower = 5,
    dftd4_damping_twobody_cso = 6,
    dftd4_damping_twobody_koide = 7,
};

/// Three-body damping function enumerator
enum dftd4_damping_threebody{
    dftd4_damping_threebody_none = -1,
    dftd4_damping_threebody_rational = 1,
    dftd4_damping_threebody_screened = 2,
    dftd4_damping_threebody_zero = 3,
    dftd4_damping_threebody_zero_avg = 4,
};

/// Damping parameter class
typedef struct _dftd4_param* dftd4_param;

/*
 * Type generic macro for convenience
**/

#define dftd4_delete(ptr) _Generic((ptr), \
                       dftd4_error: dftd4_delete_error, \
                   dftd4_structure: dftd4_delete_structure, \
                       dftd4_model: dftd4_delete_model, \
                       dftd4_param: dftd4_delete_param, \
                     dftd4_damping: dftd4_delete_damping \
                                  )(&ptr)

/*
 * Global API queries
**/

/// Obtain library version as major * 10000 + minor + 100 + patch
DFTD4_API_ENTRY int DFTD4_API_CALL
dftd4_get_version(void) DFTD4_API_SUFFIX__V_3_0;

/*
 * Error handle class
**/

/// Create new error handle object
DFTD4_API_ENTRY dftd4_error DFTD4_API_CALL
dftd4_new_error(void) DFTD4_API_SUFFIX__V_3_0;

/// Check error handle status
DFTD4_API_ENTRY int DFTD4_API_CALL
dftd4_check_error(dftd4_error /* error */) DFTD4_API_SUFFIX__V_3_0;

/// Get error message from error handle
DFTD4_API_ENTRY void DFTD4_API_CALL
dftd4_get_error(dftd4_error /* error */,
                char* /* buffer */,
                const int* /* buffersize */) DFTD4_API_SUFFIX__V_3_0;

/// Delete error handle object
DFTD4_API_ENTRY void DFTD4_API_CALL
dftd4_delete_error(dftd4_error* /* error */) DFTD4_API_SUFFIX__V_3_0;

/*
 * Molecular structure data class
**/

/// Create new molecular structure data (quantities in Bohr)
DFTD4_API_ENTRY dftd4_structure DFTD4_API_CALL
dftd4_new_structure(dftd4_error /* error */,
                    const int /* natoms */,
                    const int* /* numbers [natoms] */,
                    const double* /* positions [natoms][3] */,
                    const double* /* charge */,
                    const double* /* lattice [3][3] */,
                    const bool* /* periodic [3] */) DFTD4_API_SUFFIX__V_3_0;

/// Delete molecular structure data
DFTD4_API_ENTRY void DFTD4_API_CALL
dftd4_delete_structure(dftd4_structure* /* mol */) DFTD4_API_SUFFIX__V_3_0;

/// Update coordinates and lattice parameters (quantities in Bohr)
DFTD4_API_ENTRY void DFTD4_API_CALL
dftd4_update_structure(dftd4_error /* error */,
                       dftd4_structure /* mol */,
                       const double* /* positions [natoms][3] */,
                       const double* /* lattice [3][3] */) DFTD4_API_SUFFIX__V_3_0;

/*
 * Dispersion model class
**/

/// Create new D4 dispersion model
DFTD4_API_ENTRY dftd4_model DFTD4_API_CALL
dftd4_new_d4_model(dftd4_error /* error */,
                   dftd4_structure /* mol */) DFTD4_API_SUFFIX__V_3_0;

/// Create new D4 dispersion model
DFTD4_API_ENTRY dftd4_model DFTD4_API_CALL
dftd4_new_d4s_model(dftd4_error /* error */,
                    dftd4_structure /* mol */) DFTD4_API_SUFFIX__V_4_0;

/// Create new D4 dispersion model
DFTD4_API_ENTRY dftd4_model DFTD4_API_CALL
dftd4_custom_d4_model(dftd4_error /* error */,
                      dftd4_structure /* mol */,
                      double /* ga */,
                      double /* gc */,
                      double /* wf */) DFTD4_API_SUFFIX__V_3_1;

/// Create new D4 dispersion model
DFTD4_API_ENTRY dftd4_model DFTD4_API_CALL
dftd4_custom_d4s_model(dftd4_error /* error */,
                       dftd4_structure /* mol */,
                       double /* ga */,
                       double /* gc */) DFTD4_API_SUFFIX__V_4_0;


/// Delete dispersion model
DFTD4_API_ENTRY void DFTD4_API_CALL
dftd4_delete_model(dftd4_model* /* disp */) DFTD4_API_SUFFIX__V_3_0;

/*
 * Damping function class
**/

/// Create a new damping function with specified two-body and three-body damping
DFTD4_API_ENTRY dftd4_damping DFTD4_API_CALL
dftd4_new_damping(dftd4_error /* error */,
                  int /* damping_2b_id */,
                  int /* damping_3b_id */) DFTD4_API_SUFFIX__V_5_0;

/// Create a new default damping function for a dispersion model
DFTD4_API_ENTRY dftd4_damping DFTD4_API_CALL
dftd4_new_default_damping(dftd4_error /* error */,
                          dftd4_model /* disp */) DFTD4_API_SUFFIX__V_5_0;

/// Check the availability of the damping parameters for the use damping function
DFTD4_API_ENTRY void DFTD4_API_CALL
dftd4_check_params(dftd4_error /* error */,
                   dftd4_damping /* damp */,
                   dftd4_param /* param */) DFTD4_API_SUFFIX__V_5_0;

/// Delete damping parameters
DFTD4_API_ENTRY void DFTD4_API_CALL
dftd4_delete_damping(dftd4_damping* /* damp */) DFTD4_API_SUFFIX__V_5_0;

/*
 * Damping parameter class
**/

/// Create new damping parameters
DFTD4_API_ENTRY dftd4_param DFTD4_API_CALL
dftd4_new_param(double /* s6 */,
                double /* s8 */,
                double /* s9 */,
                double /* a1 */,
                double /* a2 */,
                double /* a3 */,
                double /* a4 */,
                double /* rs6 */,
                double /* rs8 */,
                double /* rs9 */,
                double /* alp */,
                double /* bet */) DFTD4_API_SUFFIX__V_5_0;

/// Load damping parameters from internal storage
DFTD4_API_ENTRY dftd4_param DFTD4_API_CALL
dftd4_load_param(dftd4_error /* error */,
                 char* /* method */,
                 int /* model_id */,
                 int /* damping_2b_id */,
                 int /* damping_3b_id */) DFTD4_API_SUFFIX__V_5_0;

/// Load default damping parameters for the dispersion model from internal storage
DFTD4_API_ENTRY dftd4_param DFTD4_API_CALL
dftd4_load_default_param(dftd4_error /* error */,
                         char* /* method */,
                         dftd4_model /* disp */) DFTD4_API_SUFFIX__V_5_0;

/// Delete damping parameters
DFTD4_API_ENTRY void DFTD4_API_CALL
dftd4_delete_param(dftd4_param* /* param */) DFTD4_API_SUFFIX__V_3_0;

/*
 * Perform dispersion calculations
**/

/// Evaluate properties related to the dispersion model
DFTD4_API_ENTRY void DFTD4_API_CALL
dftd4_get_properties(dftd4_error /* error */,
                     dftd4_structure /* mol */,
                     dftd4_model /* disp */,
                     double* /* cn[n] */,
                     double* /* charges[n] */,
                     double* /* c6[n*n] */,
                     double* /* alpha[n] */,
                     double* /* alphaqq[n] */) DFTD4_API_SUFFIX__V_3_1;

/// Evaluate the dispersion energy and its derivative
DFTD4_API_ENTRY void DFTD4_API_CALL
dftd4_get_dispersion(dftd4_error /* error */,
                     dftd4_structure /* mol */,
                     dftd4_model /* disp */,
                     dftd4_damping /* damp */,
                     dftd4_param /* param */,
                     double* /* energy */,
                     double* /* gradient[n][3] */,
                     double* /* sigma[3][3] */) DFTD4_API_SUFFIX__V_5_0;

/// Evaluate the dispersion hessian numerically
DFTD4_API_ENTRY void DFTD4_API_CALL
dftd4_get_numerical_hessian(dftd4_error /* error */,
                            dftd4_structure /* mol */,
                            dftd4_model /* disp */,
                            dftd4_damping /* damp */,
                            dftd4_param /* param */,
                            double* /* hess[n][3][n][3] */) DFTD4_API_SUFFIX__V_5_0;

/// Evaluate the pairwise representation of the dispersion energy
DFTD4_API_ENTRY void DFTD4_API_CALL
dftd4_get_pairwise_dispersion(dftd4_error /* error */,
                              dftd4_structure /* mol */,
                              dftd4_model /* disp */,
                              dftd4_damping /* damp */,
                              dftd4_param /* param */,
                              double* /* pair_energy2[n][n] */,
                              double* /* pair_energy3[n][n] */) DFTD4_API_SUFFIX__V_5_0;
