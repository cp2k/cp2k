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

#include <stdio.h>
#include <stdlib.h>

#include "dftd4.h"

static inline void
show_error(dftd4_error error)
{
    char message[512];
    dftd4_get_error(error, message, NULL);
    printf("[Message] %s\n", message);
}

int test_uninitialized_error(void)
{
    printf("Start test: uninitialized error\n");
    dftd4_error error = NULL;
    return dftd4_check_error(error) ? 0 : 1;
}

int test_uninitialized_structure(void)
{
    printf("Start test: uninitialized structure\n");
    dftd4_error error = NULL;
    dftd4_structure mol = NULL;

    error = dftd4_new_error();

    double xyz[6] = { 0.0 };
    dftd4_update_structure(error, mol, xyz, NULL);
    if (!dftd4_check_error(error))
        goto unexpected;

    show_error(error);

    dftd4_delete(error);
    return 0;

unexpected:
    printf("[Fatal] Unexpected pass for unititalized-structure test\n");
    dftd4_delete(error);
    return 1;
}

int test_example(void)
{
    printf("Start test: example\n");
    int const natoms = 7;
    int const nat_sq = natoms * natoms;
    int const nat3 = natoms * 3;
    int const nat3_sq = nat3 * nat3;
    int const attyp[7] = { 6, 6, 6, 1, 1, 1, 1 };
    double const coord[21] = {
        +0.00000000000000, +0.00000000000000, -1.79755622305860,
        +0.00000000000000, +0.00000000000000, +0.95338756106749,
        +0.00000000000000, +0.00000000000000, +3.22281255790261,
        -0.96412815539807, -1.66991895015711, -2.53624948351102,
        -0.96412815539807, +1.66991895015711, -2.53624948351102,
        +1.92825631079613, +0.00000000000000, -2.53624948351102,
        +0.00000000000000, +0.00000000000000, +5.23010455462158 };
    double energy;
    double sigma[9];
    double* pair_disp2;
    double* pair_disp3;
    double* gradient;
    double* hessian;
    double* c6;

    pair_disp2 = (double*)malloc(nat_sq * sizeof(double));
    pair_disp3 = (double*)malloc(nat_sq * sizeof(double));
    gradient = (double*)malloc(nat3 * sizeof(double));
    hessian = (double*)malloc(nat3_sq * sizeof(double));
    c6 = (double*)malloc(nat_sq * sizeof(double));

    dftd4_error error = NULL;

    if (dftd4_get_version() <= 0) {
        goto err;
    }

    dftd4_structure mol;
    dftd4_model disp;
    dftd4_damping damp;
    dftd4_param param;

    error = dftd4_new_error();
    if (!error) {
        goto err;
    }

    mol = dftd4_new_structure(error, natoms, attyp, coord, NULL, NULL, NULL);
    if (dftd4_check_error(error)) {
        goto err;
    };
    if (!mol) {
        goto err;
    }

    disp = dftd4_new_d4_model(error, mol);
    if (dftd4_check_error(error)) {
        goto err;
    }
    if (!disp) {
        goto err;
    }

    // C6 coefficients
    dftd4_get_properties(error, mol, disp, NULL, NULL, c6, NULL, NULL);
    if (dftd4_check_error(error)) {
        goto err;
    }

    // PBE-D4
    param = dftd4_new_param(1.0, 0.95948085, 0.0, 0.38574991, 4.80688534, 
        0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 0.0);
    if (dftd4_check_error(error)) {
        goto err;
    }
    if (!param) {
        goto err;
    }

    // Setup default damping functions for DFT-D4
    damp = dftd4_new_default_damping(error, disp);
    if (dftd4_check_error(error)) {
       goto err;
    }
    dftd4_get_dispersion(error, mol, disp, damp, param, &energy, NULL, NULL);
    if (dftd4_check_error(error)) {
        goto err;
    }
    dftd4_get_dispersion(error, mol, disp, damp, param, &energy, gradient, sigma);
    if (dftd4_check_error(error)) {
        goto err;
    }

    dftd4_get_numerical_hessian(error, mol, disp, damp, param, hessian);
    if (dftd4_check_error(error)) {
        goto err;
    }
    dftd4_get_pairwise_dispersion(error, mol, disp, damp, param, pair_disp2, pair_disp3);
    if (dftd4_check_error(error)) {
        goto err;
    }
    dftd4_delete(param);

    // DSD-BLYP-D4-ATM
    param = dftd4_load_default_param(error, "dsdblyp", disp);
    if (dftd4_check_error(error)) {
        goto err;
    }
    if (!param) {
        goto err;
    }

    dftd4_get_dispersion(error, mol, disp, damp, param, &energy, NULL, NULL);
    if (dftd4_check_error(error)) {
        goto err;
    }
    dftd4_get_dispersion(error, mol, disp, damp, param, &energy, gradient, sigma);
    if (dftd4_check_error(error)) {
        goto err;
    }
    dftd4_get_numerical_hessian(error, mol, disp, damp, param, hessian);
    if (dftd4_check_error(error)) {
        goto err;
    }
    dftd4_delete(param);
    dftd4_delete(disp);

    // D4S - r2SCAN-3c
    disp = dftd4_custom_d4s_model(error, mol, 2.0, 1.0);
    if (dftd4_check_error(error)) {
        goto err;
    }
    if (!disp) {
        goto err;
    }

    param = dftd4_load_default_param(error, "r2scan_3c", disp);
    if (dftd4_check_error(error)) {
        goto err;
    }
    if (!param) {
        goto err;
    }

    dftd4_get_dispersion(error, mol, disp, damp, param, &energy, NULL, NULL);
    if (dftd4_check_error(error)) {
        goto err;
    }

    dftd4_get_dispersion(error, mol, disp, damp, param, &energy, gradient, sigma);
    if (dftd4_check_error(error)) {
        goto err;
    }
    dftd4_get_numerical_hessian(error, mol, disp, damp, param, hessian);
    if (dftd4_check_error(error)) {
        goto err;
    }
    dftd4_delete(param);
    dftd4_delete(damp);

    // modified TPSSh-D4(sc)-ATM(sc)
    param = dftd4_new_param(1.0, 1.76596355, 0.5, 0.42822303, 4.54257102, 0.6, 0.6, 
      0.0, 0.0, 0.0, 16.0, 0.0);
    if (dftd4_check_error(error)) {
        goto err;
    }
    if (!param) {
        goto err;
    }

    // Setup screened damping functions
    damp = dftd4_new_damping(error, dftd4_damping_twobody_screened, dftd4_damping_threebody_screened);
    if (dftd4_check_error(error)) {
        goto err;
    }
    dftd4_check_params(error, damp, param);
    if (dftd4_check_error(error)) {
        goto err;
    }
    dftd4_get_dispersion(error, mol, disp, damp, param, &energy, NULL, NULL);
    if (dftd4_check_error(error)) {
        goto err;
    }    
    dftd4_get_dispersion(error, mol, disp, damp, param, &energy, gradient, sigma);
    if (dftd4_check_error(error)) {
        goto err;
    }
    dftd4_delete(param);
    dftd4_delete(damp);
    dftd4_delete(disp);

    // Attempt to create custom d4 model
    disp = dftd4_custom_d4_model(error, mol, 2.0, 1.0, 4.0);
    if (dftd4_check_error(error)) {
        goto err;
    }
    if (!disp) {
        goto err;
    }

    dftd4_delete(disp);

    // Attempt to create default d4s model
    disp = dftd4_new_d4s_model(error, mol);
    if (dftd4_check_error(error)) {
        goto err;
    }
    if (!disp) {
        goto err;
    }

    dftd4_delete(disp);
    dftd4_delete(mol);
    dftd4_delete(error);

    if (param) {
        goto err;
    }
    if (disp) {
        goto err;
    }
    if (mol) {
        goto err;
    }
    if (error) {
        goto err;
    }

    free(pair_disp2);
    free(pair_disp3);
    free(gradient);
    free(hessian);
    free(c6);

    return 0;

err:
    if (dftd4_check_error(error)) {
        char message[512];
        dftd4_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    dftd4_delete(param);
    dftd4_delete(disp);
    dftd4_delete(mol);
    dftd4_delete(error);

    free(pair_disp2);
    free(pair_disp3);
    free(gradient);
    free(hessian);
    free(c6);

    return 1;
}

int test_atm_toggle(void)
{
    printf("Start test: atm toggle check\n");
    
    // Setup molecule (same as test_example)
    int const natoms = 7;
    int const attyp[7] = { 6, 6, 6, 1, 1, 1, 1 };
    double const coord[21] = {
        +0.00000000000000, +0.00000000000000, -1.79755622305860,
        +0.00000000000000, +0.00000000000000, +0.95338756106749,
        +0.00000000000000, +0.00000000000000, +3.22281255790261,
        -0.96412815539807, -1.66991895015711, -2.53624948351102,
        -0.96412815539807, +1.66991895015711, -2.53624948351102,
        +1.92825631079613, +0.00000000000000, -2.53624948351102,
        +0.00000000000000, +0.00000000000000, +5.23010455462158 };

    double energy_atm_on = 0.0;
    double energy_atm_off = 0.0;
    
    dftd4_error error = dftd4_new_error();

    dftd4_structure mol = dftd4_new_structure(error, natoms, attyp, coord, NULL, NULL, NULL);
    if (!mol || dftd4_check_error(error)) goto err;

    dftd4_model disp = dftd4_new_d4_model(error, mol);
    if (!disp || dftd4_check_error(error)) goto err;

    // Calculate with ATM term
    dftd4_damping damp = dftd4_new_default_damping(error, disp);
    damp = dftd4_new_damping(error, dftd4_damping_twobody_rational, 
      dftd4_damping_threebody_zero_avg);
    dftd4_param param_on = dftd4_load_param(error, "pbe", dftd4_model_d4, 
      dftd4_damping_twobody_rational, dftd4_damping_threebody_zero_avg);
    if (!param_on || dftd4_check_error(error)) goto err;
    
    dftd4_get_dispersion(error, mol, disp, damp, param_on, &energy_atm_on, NULL, NULL);
    if (dftd4_check_error(error)) goto err;
    dftd4_delete(damp);
    dftd4_delete(param_on);

    // Calculate without ATM
    damp = dftd4_new_damping(error, dftd4_damping_twobody_rational, 
      dftd4_damping_threebody_none);
    dftd4_param param_off = dftd4_load_param(error, "pbe", dftd4_model_d4, 
      dftd4_damping_twobody_rational, dftd4_damping_threebody_none);
    if (!param_off || dftd4_check_error(error)) goto err;

    dftd4_get_dispersion(error, mol, disp, damp, param_off, &energy_atm_off, NULL, NULL);
    if (dftd4_check_error(error)) goto err;
    dftd4_delete(param_off);

    // If the bug exists (https://github.com/dftd4/dftd4/issues/333), 'false' 
    // is treated as 'true', and energies will be identical. The difference 
    // will be small though because the ATM term itself is small.
    printf("Energy (ATM=on):  %.8f\n", energy_atm_on);
    printf("Energy (ATM=off): %.8f\n", energy_atm_off);

    double diff = energy_atm_on - energy_atm_off;
    if (diff < 0) diff = -diff;

    if (diff < 1e-10) {
        printf("[Fatal] ATM toggle failed. Energies are identical (Diff: %e)\n", diff);
        goto err;
    }

    // Cleanup
    dftd4_delete(disp);
    dftd4_delete(mol);
    dftd4_delete(error);
    return 0;

err:
    if (dftd4_check_error(error)) {
        show_error(error);
    }
    dftd4_delete(disp);
    dftd4_delete(mol);
    dftd4_delete(error);
    return 1;
}

int main(void)
{
    int stat = 0;
    stat += test_uninitialized_error();
    stat += test_uninitialized_structure();
    stat += test_example();
    stat += test_atm_toggle();

    return stat == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
