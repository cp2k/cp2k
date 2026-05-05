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

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tblite.h"

#define check(x, ...)           \
    _Generic((x),               \
             int                \
             : check_int,       \
             int*               \
             : check_int_array, \
             double             \
             : check_double,    \
             double*            \
             : check_double_array)(x, __VA_ARGS__)

static inline bool
check_int(int actual, int expected, const char* msg)
{
    if (expected == actual) {
        return true;
    }
    fprintf(stderr, "[Fatal] %s: expected %d, got %d\n", msg, expected, actual);
    return false;
}

static inline bool
check_double(double actual, double expected, double tol, const char* msg)
{
    if (fabs(expected - actual) < tol) {
        return true;
    }
    fprintf(stderr, "[Fatal] %s: expected %3.7f, got %3.7f\n", msg, expected, actual);
    return false;
}

static inline bool
check_double_array(double* actual, double* expected, double tol, int ndim, const char* msg)
{
    for (int i = 0; i != ndim; i++) {
        if (!check_double(actual[i], expected[i], tol, msg)) {
            return false;
        }
    }
    return true;
}

static inline bool
check_int_array(const int* actual, const int* expected, int ndim, const char* msg)
{
    for (int i = 0; i != ndim; i++) {
        if (!check(actual[i], expected[i], msg)) {
            return false;
        }
    }
    return true;
}

static inline double
norm2(int n, double* vec)
{
    double norm = 0.0;
    int i;
    for (i = 0; i != n; i++)
        norm += vec[i] * vec[i];
    return sqrt(norm);
}

static inline double
sum(int n, double* vec)
{
    double val = 0.0;
    int i;
    for (i = 0; i != n; i++)
        val += vec[i];
    return val;
}

#define show(x) \
    _Generic((x), \
        tblite_error: show_error, \
        tblite_context: show_context_error \
    )((x))

static inline void
show_error(tblite_error error)
{
    char message[512];
    tblite_get_error(error, message, NULL);
    printf("[Message] %s\n", message);
    tblite_clear_error(error);
}

static inline void
show_context_error(tblite_context ctx)
{
    char message[512];
    tblite_get_context_error(ctx, message, NULL);
    printf("[Message] %s\n", message);
}

static inline tblite_structure
get_structure_1(tblite_error error)
{
    int natoms = 22;
    int num[22] = { 7, 1, 6, 1, 6, 6, 1, 1, 1, 8, 6, 6, 1, 1, 1, 8, 7, 1, 6, 1, 1, 1 };
    double xyz[66] = {
        2.65893135608838,
        -2.39249423371715,
        -3.66065400053935,
        3.49612941769371,
        -0.88484673975624,
        -2.85194146578362,
        -0.06354076626069,
        -2.63180732150005,
        -3.28819116275323,
        -1.07444177498884,
        -1.92306930149582,
        -4.93716401361053,
        -0.83329925447427,
        -5.37320588052218,
        -2.81379718546920,
        -0.90691285352090,
        -1.04371377845950,
        -1.04918016247507,
        -2.86418317801214,
        -5.46484901686185,
        -2.49961410229771,
        -0.34235262692151,
        -6.52310417728877,
        -4.43935278498325,
        0.13208660968384,
        -6.10946566962768,
        -1.15032982743173,
        -2.96060093623907,
        0.01043357425890,
        -0.99937552379387,
        3.76519127865000,
        -3.27106236675729,
        -5.83678272799149,
        6.47957316843231,
        -2.46911747464509,
        -6.21176914665408,
        7.32688324906998,
        -1.67889171278096,
        -4.51496113512671,
        6.54881843238363,
        -1.06760660462911,
        -7.71597456720663,
        7.56369260941896,
        -4.10015651865148,
        -6.82588105651977,
        2.64916867837331,
        -4.60764575400925,
        -7.35167957128511,
        0.77231592220237,
        -0.92788783332000,
        0.90692539619101,
        2.18437036702702,
        -2.20200039553542,
        0.92105755612696,
        0.01367202674183,
        0.22095199845428,
        3.27728206652909,
        1.67849497305706,
        0.53855308534857,
        4.43416031916610,
        -0.89254709011762,
        2.01704896333243,
        2.87780123699499,
        -1.32658751691561,
        -0.95404596601807,
        4.30967630773603,
    };

    return tblite_new_structure(
        error,
        natoms,
        num,
        xyz,
        NULL,
        NULL,
        NULL,
        NULL);
}

static inline tblite_structure
get_structure_2(tblite_error error)
{
    int natoms = 17;
    int num[17] = { 14, 14, 14, 14, 14, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    double xyz[51] = {
        +1.19317075263240,
        +2.10734582303615,
        -3.94069302326331,
        +1.32974747794452,
        +5.74383319654784,
        -1.40646299357997,
        +1.36246249005040,
        -1.51437759946738,
        -1.38835921669240,
        -1.86383433000121,
        -1.43677027343543,
        +1.65343257587405,
        -1.84846922430006,
        -5.11625835608814,
        +4.12639172175940,
        -1.21073388457482,
        +2.08204436847949,
        -5.39986048218966,
        +3.32212225093368,
        +2.13580609828116,
        -5.77918927003773,
        +3.84368159195704,
        -1.60306480080261,
        -0.06684474846163,
        +1.13319183423377,
        -3.84842795885741,
        -2.94128196218086,
        -4.34788687075438,
        -1.16479029631036,
        +0.36221328927203,
        -1.50552338174053,
        +0.81829285234045,
        +3.29535409206112,
        -3.87198157702809,
        -5.06545254280118,
        +6.07282119819235,
        -2.23981149763562,
        -7.36061632570679,
        +2.48584729172346,
        +0.63031417166410,
        -5.37888394172819,
        +5.41803105663996,
        +1.37439281298519,
        +8.08088049648159,
        -2.96190852429966,
        -0.92666191369478,
        +5.81989514634268,
        +0.26140479183806,
        +3.62581929732841,
        +5.70054411368805,
        +0.20910420334480,
    };

    return tblite_new_structure(
        error,
        natoms,
        num,
        xyz,
        NULL,
        NULL,
        NULL,
        NULL);
}

static inline tblite_structure
get_structure_3(tblite_error error, bool charge)
{
    int natoms = 5;
    int num[5] = { 6, 1, 1, 1, 1 };
    double charge_cation = 1.0;
    int uhf_cation = 1;
    double xyz_cation[15] = {
        0.00084274718833,
        -0.00015740451438,
        0.00006401702382,
        -1.55577637922762,
        1.24722361298503,
        -0.70873815063999,
        1.24747864259928,
        1.55718783736300,
        0.70764030938274,
        1.55611738014180,
        -1.24919736758084,
        -0.70740441319739,
        -1.24866239070178,
        -1.55505667825279,
        0.70843823743083,
    };
    double charge_neutral = 0.0;
    int uhf_neutral = 0;
    double xyz_neutral[3 * 5] = {
        0.00000000000000,
        0.00000000000000,
        0.00000000000000,
        1.18771160655551,
        -1.18771160655551,
        1.18771160655551,
        -1.18771160655551,
        1.18771160655551,
        1.18771160655551,
        -1.18771160655551,
        -1.18771160655551,
        -1.18771160655551,
        1.18771160655551,
        1.18771160655551,
        -1.18771160655551,
    };

    return tblite_new_structure(
        error,
        natoms,
        num,
        charge ? xyz_cation : xyz_neutral,
        charge ? &charge_cation : &charge_neutral,
        charge ? &uhf_cation : &uhf_neutral,
        NULL,
        NULL);
}

static inline tblite_structure
get_structure_4(tblite_error error)
{
    int natoms = 12;
    int num[12] = { 6, 6, 6, 6, 8, 8, 8, 8, 8, 8, 8, 8 };
    double xyz[36] = {
        0.00000000000000,
        0.00000000000000,
        0.00000000000000,
        5.47113387782808,
        0.00000000000000,
        5.47113387782808,
        5.47113387782808,
        5.47113387782808,
        0.00000000000000,
        0.00000000000000,
        5.47113387782808,
        5.47113387782808,
        1.27461999191594,
        1.27461999191594,
        1.27461999191594,
        6.74556489717295,
        1.27461999191594,
        4.19651388591214,
        6.74556489717295,
        4.19651388591214,
        9.66764776374022,
        9.66764776374022,
        6.74556489717295,
        4.19651388591214,
        9.66764776374022,
        9.66764776374022,
        9.66764776374022,
        4.19651388591214,
        9.66764776374022,
        6.74556489717295,
        4.19651388591214,
        6.74556489717295,
        1.27461999191594,
        1.27461999191594,
        4.19651388591214,
        6.74556489717295,
    };
    double lattice[9] = {
        10.94216438765978,
        0.00000000000000,
        0.00000000000000,
        0.00000000000000,
        10.94216438765978,
        0.00000000000000,
        0.00000000000000,
        0.00000000000000,
        10.94216438765978,
    };
    bool periodic[3] = { true };

    return tblite_new_structure(
        error,
        natoms,
        num,
        xyz,
        NULL,
        NULL,
        lattice,
        periodic);
}

static inline tblite_structure
get_structure_5(tblite_error error)
{
    int natoms = 21;
    int num[21] = { 24, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 6, 6, 6, 1, 6, 1, 6, 1, 1, 1 };
    double xyz[3 * 21] = {
        +0.00000000000000, +0.00000000000000, -0.06044684528305,
        +0.00000000000000, +3.19613712523833, +2.30877824528580,
        +2.18828801115897, +3.32943780995850, +0.70249948585735,
        +1.33235791539260, +3.55640652898451, -1.83908673090077,
        -1.33235791539260, +3.55640652898451, -1.83908673090077,
        -2.18828801115897, +3.32943780995850, +0.70249948585735,
        +0.00000000000000, +3.10509505378016, +4.34935395653655,
        +4.13810718850644, +3.28428734944129, +1.31235006648465,
        +2.52190264478215, +3.60569548880831, -3.50208900904436,
        -2.52190264478215, +3.60569548880831, -3.50208900904436,
        -4.13810718850644, +3.28428734944129, +1.31235006648465,
        +2.18828801115897, -3.32943780995850, +0.70249948585735,
        +0.00000000000000, -3.19613712523833, +2.30877824528580,
        +1.33235791539260, -3.55640652898451, -1.83908673090077,
        +4.13810718850644, -3.28428734944129, +1.31235006648465,
        -2.18828801115897, -3.32943780995850, +0.70249948585735,
        +0.00000000000000, -3.10509505378016, +4.34935395653655,
        -1.33235791539260, -3.55640652898451, -1.83908673090077,
        +2.52190264478215, -3.60569548880831, -3.50208900904436,
        -4.13810718850644, -3.28428734944129, +1.31235006648465,
        -2.52190264478215, -3.60569548880831, -3.50208900904436,
    };

    return tblite_new_structure(
        error,
        natoms,
        num,
        xyz,
        NULL,
        NULL,
        NULL,
        NULL);
}

int test_version(void)
{
    printf("Start test: version\n");
    return tblite_get_version() > 0 ? 0 : 1;
}

int test_uninitialized_error(void)
{
    printf("Start test: uninitialized error\n");
    tblite_error error = NULL;
    return tblite_check(error) ? 0 : 1;
}

int test_uninitialized_context(void)
{
    printf("Start test: uninitialized context\n");
    tblite_context ctx = NULL;
    return tblite_check(ctx) ? 0 : 1;
}

int test_uninitialized_structure(void)
{
    printf("Start test: uninitialized structure\n");
    tblite_error error = NULL;
    tblite_structure mol = NULL;

    error = tblite_new_error();

    double xyz[6] = { 0.0 };
    tblite_update_structure_geometry(error, mol, xyz, NULL);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    double charge = 0.0;
    tblite_update_structure_charge(error, mol, &charge);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    int uhf = 1;
    tblite_update_structure_uhf(error, mol, &uhf);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_delete(error);
    return 0;

unexpected:
    printf("[Fatal] Unexpected pass for unititalized-structure test\n");
    tblite_delete(error);
    return 1;
}

int test_uninitialized_result(void)
{
    printf("Start test: uninitialized result\n");
    tblite_error error = NULL;
    tblite_result res = NULL;

    error = tblite_new_error();

    int natoms;
    tblite_get_result_number_of_atoms(error, res, &natoms);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    int nshells;
    tblite_get_result_number_of_shells(error, res, &nshells);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    int norbs;
    tblite_get_result_number_of_orbitals(error, res, &norbs);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    int nspin;
    tblite_get_result_number_of_spins(error, res, &nspin);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    double energy;
    tblite_get_result_energy(error, res, &energy);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    double gradient[12];
    tblite_get_result_gradient(error, res, gradient);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    double sigma[9];
    tblite_get_result_virial(error, res, sigma);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    double charges[7];
    tblite_get_result_charges(error, res, charges);
    if (!tblite_check(error))
        goto unexpected;

    double mbo[20];
    tblite_get_result_bond_orders(error, res, mbo);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    show(error);

    double dipole[3];
    tblite_get_result_dipole(error, res, dipole);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    double quadrupole[6];
    tblite_get_result_quadrupole(error, res, quadrupole);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    double emo[7];
    tblite_get_result_orbital_energies(error, res, emo);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    double occ[7];
    tblite_get_result_orbital_occupations(error, res, occ);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    double mat[49];
    tblite_get_result_orbital_coefficients(error, res, mat);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_get_result_density_matrix(error, res, mat);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_get_result_overlap_matrix(error, res, mat);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_get_result_hamiltonian_matrix(error, res, mat);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_load_result_wavefunction(error, res, ".uninit.npz");
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    res = tblite_new_result();
    tblite_load_result_wavefunction(error, res, ".uninit.npz");
    if (!tblite_check(error))
        goto unexpected;

    tblite_delete(res);
    show(error);

    tblite_save_result_wavefunction(error, res, ".uninit.npz");
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    res = tblite_new_result();
    tblite_save_result_wavefunction(error, res, ".uninit.npz");
    if (!tblite_check(error))
        goto unexpected;

    tblite_delete(res);
    show(error);

    tblite_delete(error);
    tblite_delete(res);
    return 0;

unexpected:
    printf("[Fatal] Unexpected pass for unititalized-result test\n");
    tblite_delete(error);
    return 1;
}

int test_uninitialized_calculator(void)
{
    printf("Start test: uninitialized calculator\n");
    tblite_context ctx = NULL;
    tblite_calculator calc = NULL;
    tblite_container cont = NULL;

    ctx = tblite_new_context();

    tblite_set_calculator_accuracy(ctx, calc, 1.0);
    if (!tblite_check(ctx))
        goto unexpected;

    show(ctx);

    tblite_set_calculator_temperature(ctx, calc, 0.0);
    if (!tblite_check(ctx))
        goto unexpected;

    show(ctx);

    tblite_set_calculator_max_iter(ctx, calc, 12);
    if (!tblite_check(ctx))
        goto unexpected;

    show(ctx);

    tblite_set_calculator_mixer_damping(ctx, calc, 0.2);
    if (!tblite_check(ctx))
        goto unexpected;

    show(ctx);

    tblite_set_calculator_guess(ctx, calc, TBLITE_GUESS_SAD);
    if (!tblite_check(ctx))
        goto unexpected;

    show(ctx);

    tblite_set_calculator_save_integrals(ctx, calc, 1);
    if (!tblite_check(ctx))
        goto unexpected;

    show(ctx);

    tblite_calculator_push_back(ctx, calc, &cont);
    if (!tblite_check(ctx))
        goto unexpected;

    show(ctx);

    int nsh;
    tblite_get_calculator_shell_count(ctx, calc, &nsh);
    if (!tblite_check(ctx))
        goto unexpected;

    show(ctx);

    int sh2at[5];
    tblite_get_calculator_shell_map(ctx, calc, sh2at);
    if (!tblite_check(ctx))
        goto unexpected;

    show(ctx);

    int am[5];
    tblite_get_calculator_angular_momenta(ctx, calc, am);
    if (!tblite_check(ctx))
        goto unexpected;

    show(ctx);

    int nao;
    tblite_get_calculator_orbital_count(ctx, calc, &nao);
    if (!tblite_check(ctx))
        goto unexpected;

    show(ctx);

    int ao2sh[9];
    tblite_get_calculator_orbital_map(ctx, calc, ao2sh);
    if (!tblite_check(ctx))
        goto unexpected;

    show(ctx);

    tblite_delete(ctx);
    return 0;

unexpected:
    printf("[Fatal] Unexpected pass for unititalized-calculator test\n");
    tblite_delete(ctx);
    return 1;
}

int test_uninitialized_table(void)
{
    printf("Start test: uninitialized table\n");
    tblite_error error = NULL;
    tblite_param param = NULL;
    tblite_table table = NULL;
    char key[] = "some-key";

    error = tblite_new_error();
    param = tblite_new_param();

    tblite_load_param(error, param, table);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_dump_param(error, param, table);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    double dval = 0.0;
    tblite_table_set_value(error, table, key, &dval, 0);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    int64_t ival = 0;
    tblite_table_set_value(error, table, key, &ival, 0);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    bool lval = 0;
    tblite_table_set_value(error, table, key, &lval, 0);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    char cval[] = "some-val";
    tblite_table_set_value(error, table, key, &cval, 0);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_table tval = tblite_table_add_table(error, table, key);
    tblite_delete(tval);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_delete(error);
    tblite_delete(param);
    return 0;

unexpected:
    printf("[Fatal] Unexpected pass for unititalized-table test\n");
    tblite_delete(error);
    tblite_delete(param);
    return 1;
}

int test_uninitialized_param(void)
{
    printf("Start test: uninitialized param\n");
    tblite_error error = NULL;
    tblite_param param = NULL;
    tblite_table table = NULL;

    error = tblite_new_error();
    table = tblite_new_table(NULL);

    tblite_load_param(error, param, table);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_dump_param(error, param, table);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_export_gfn2_param(error, param);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_export_gfn1_param(error, param);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_export_ipea1_param(error, param);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_delete(error);
    tblite_delete(table);
    return 0;

unexpected:
    printf("[Fatal] Unexpected pass for unititalized-param test\n");
    tblite_delete(error);
    tblite_delete(table);
    return 1;
}

int test_error_setter(void)
{
    printf("Start test: error set\n");
    tblite_error error = tblite_new_error();
    tblite_set_error(error, "some error message", NULL);
    return tblite_check(error) ? 0 : 1;
}

int test_empty_result(void)
{
    printf("Start test: empty result\n");
    tblite_error error = NULL;
    tblite_result res = NULL;

    error = tblite_new_error();
    res = tblite_new_result();

    int nspin;
    tblite_get_result_number_of_spins(error, res, &nspin);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    double energy;
    tblite_get_result_energy(error, res, &energy);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    double gradient[12];
    tblite_get_result_gradient(error, res, gradient);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    double sigma[12];
    tblite_get_result_virial(error, res, sigma);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    double energies[4];
    tblite_get_result_energies(error, res, energies);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    double mat[10];
    tblite_get_result_density_matrix(error, res, mat);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_get_result_overlap_matrix(error, res, mat);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_get_result_hamiltonian_matrix(error, res, mat);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_delete(error);
    return 0;

unexpected:
    printf("[Fatal] Unexpected pass for empty-result test\n");
    tblite_delete(error);
    return 1;
}

int test_invalid_structure(void)
{
    printf("Start test: invalid structure\n");
    tblite_error error = NULL;
    tblite_structure mol = NULL;

    int natoms = 2;
    int num[2] = { 1, 1 };
    double xyz[6] = { 0.0 };

    error = tblite_new_error();

    mol = tblite_new_structure(error, natoms, num, xyz, NULL, NULL, NULL, NULL);
    if (!tblite_check(error))
        goto unexpected;

    show(error);

    tblite_delete(error);
    tblite_delete(mol);
    return 0;

unexpected:
    printf("[Fatal] Unexpected pass for invalid-structure test\n");
    tblite_delete(error);
    tblite_delete(mol);
    return 1;
}

int test_table_builder(void)
{
    printf("Start test: table-builder\n");
    tblite_error error = NULL;
    tblite_table table = NULL, child1 = NULL, child2 = NULL, child3 = NULL;
    double dval;
    double darr[2];
    char cval[4];
    char carr[3][2];
    int64_t iarr[2];

    error = tblite_new_error();

    table = tblite_new_table(NULL);
    child1 = tblite_table_add_table(error, table, "hamiltonian");
    if (tblite_check(error))
        goto err;
    if (!child1)
        goto err;
    child2 = tblite_table_add_table(error, child1, "xtb");
    if (tblite_check(error))
        goto err;
    if (!child2)
        goto err;

    dval = 0.5;
    tblite_table_set_value(error, child2, "wexp", &dval, 0);
    if (tblite_check(error))
        goto err;

    dval = 2.0e-2;
    tblite_table_set_value(error, child2, "enscale", &dval, 0);
    if (tblite_check(error))
        goto err;

    strcpy(cval, "gfn");
    tblite_table_set_value(error, child2, "cn", &cval, 0);
    if (tblite_check(error))
        goto err;

    child3 = tblite_table_add_table(error, child2, "shell");
    if (tblite_check(error))
        goto err;
    if (!child3)
        goto err;

    dval = 1.85;
    tblite_table_set_value(error, child3, "ss", &dval, 0);
    if (tblite_check(error))
        goto err;

    dval = 2.23;
    tblite_table_set_value(error, child3, "pp", &dval, 0);
    if (tblite_check(error))
        goto err;

    tblite_delete(child1);
    tblite_delete(child2);
    tblite_delete(child3);

    child1 = tblite_table_add_table(error, table, "element");
    if (tblite_check(error))
        goto err;
    if (!child1)
        goto err;
    child2 = tblite_table_add_table(error, child1, "C");
    if (tblite_check(error))
        goto err;
    if (!child2)
        goto err;

    strcpy(carr[0], "2s");
    strcpy(carr[1], "2p");
    tblite_table_set_value(error, child2, "shells", carr, 2);
    if (tblite_check(error))
        goto err;

    darr[0] = -13.970922;
    darr[1] = -10.063292;
    tblite_table_set_value(error, child2, "levels", darr, 2);
    if (tblite_check(error))
        goto err;

    darr[0] = 2.096432;
    darr[1] = 1.8;
    tblite_table_set_value(error, child2, "slater", darr, 2);
    if (tblite_check(error))
        goto err;

    iarr[0] = 4;
    iarr[1] = 4;
    tblite_table_set_value(error, child2, "ngauss", iarr, 2);
    if (tblite_check(error))
        goto err;

    tblite_delete(error);
    tblite_delete(table);
    tblite_delete(child1);
    tblite_delete(child2);
    tblite_delete(child3);

    return 0;

err:
    tblite_delete(error);
    tblite_delete(table);
    tblite_delete(child1);
    tblite_delete(child2);
    tblite_delete(child3);
    return 1;
}

int test_param_load(void)
{
    printf("Start test: param-load\n");
    tblite_error error = NULL;
    tblite_param param = NULL;
    tblite_table table = NULL;

    error = tblite_new_error();
    param = tblite_new_param();

    tblite_export_gfn1_param(error, param);
    if (tblite_check(error))
        goto err;

    tblite_export_gfn2_param(error, param);
    if (tblite_check(error))
        goto err;

    table = tblite_new_table(NULL);
    tblite_dump_param(error, param, table);
    if (tblite_check(error))
        goto err;

    tblite_delete(param);
    param = tblite_new_param();

    tblite_load_param(error, param, table);
    if (tblite_check(error))
        goto err;

    tblite_export_ipea1_param(error, param);
    if (tblite_check(error))
        goto err;

    tblite_delete(error);
    tblite_delete(table);
    tblite_delete(param);

    return 0;

err:
    tblite_delete(error);
    tblite_delete(table);
    tblite_delete(param);
    return 1;
}

int test_calc_restart(void)
{
    printf("Start test: calculator-restart\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_result res1 = NULL, res2 = NULL;

    const double thr = 5.0e-7;
    double energy;
    double energies[22];

    error = tblite_new_error();
    ctx = tblite_new_context();
    res1 = tblite_new_result();

    mol = get_structure_1(error);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_gfn1_calculator(ctx, mol);
    if (!calc)
        goto err;

    tblite_set_calculator_accuracy(ctx, calc, 2.0);
    tblite_set_calculator_max_iter(ctx, calc, 50);
    tblite_set_calculator_mixer_damping(ctx, calc, 0.2);
    tblite_set_calculator_temperature(ctx, calc, 0.0);

    tblite_get_singlepoint(ctx, mol, calc, res1);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res1, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -34.98079463818, thr, "energy error"))
        goto err;

    // reset calculator
    tblite_delete(calc);
    calc = tblite_new_gfn2_calculator(ctx, mol);

    tblite_get_singlepoint(ctx, mol, calc, res1);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res1, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -32.96247211794, thr, "energy error"))
        goto err;

    tblite_get_result_energies(error, res1, energies);
    if (tblite_check(error))
        goto err;

    if (!check(energy, sum(22, energies), thr, "energy error"))
        goto err;

    res2 = tblite_copy_result(res1);
    tblite_delete(res1);
    tblite_set_calculator_max_iter(ctx, calc, 3);

    tblite_get_singlepoint(ctx, mol, calc, res2);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res2, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -32.96247199299, thr, "energy error"))
        goto err;

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(res2);
    return 0;

err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(res1);
    tblite_delete(res2);
    return 1;
}

void example_callback(tblite_error error, char* msg, int len, void* udata)
{
    printf("[callback] %.*s\n", len, msg);
}

int test_callback(void)
{
    printf("Start test: callback\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_result res = NULL;

    const double thr = 5.0e-7;
    int natoms, nshells, norb;
    double energy;

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();

    tblite_logger_callback callback = example_callback;
    tblite_set_context_logger(ctx, callback, NULL);
    tblite_set_context_color(ctx, 1);

    mol = get_structure_1(error);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_gfn2_calculator(ctx, mol);

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_number_of_atoms(error, res, &natoms);
    if (tblite_check(error))
        goto err;

    if (!check(natoms, 22, "dimension error"))
        goto err;

    tblite_get_result_number_of_shells(error, res, &nshells);
    if (tblite_check(error))
        goto err;

    if (!check(nshells, 32, "dimension error"))
        goto err;

    tblite_get_result_number_of_orbitals(error, res, &norb);
    if (tblite_check(error))
        goto err;

    if (!check(norb, 52, "dimension error"))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -32.96247211794, thr, "energy error"))
        goto err;

    tblite_set_context_logger(ctx, NULL, NULL);
    tblite_set_context_color(ctx, 0);
    tblite_set_calculator_max_iter(ctx, calc, 3);

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -32.96247195792, thr, "energy error"))
        goto err;

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(res);
    return 0;

err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(res);
    return 1;
}

int test_gfn2_si5h12(void)
{
    printf("Start test: GFN2-xTB\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_result res = NULL;
    tblite_param param = NULL;

    const double thr = 5.0e-7;
    double energy;
    double gradient[51] = { 0.0 };
    int nshells, nao;
    int sh2at[27], am[27], ao2sh[57];
    const int ref_sh2at[27] = {
        0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16
    };
    const int ref_am[27] = {
        0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    const int ref_ao2sh[57] = {
        0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 4, 4, 4, 5, 5, 5, 5, 5, 6, 7, 7, 7, 8, 8, 8, 8, 8, 9,
        10, 10, 10, 11, 11, 11, 11, 11, 12, 13, 13, 13, 14, 14, 14, 14, 14, 15, 16, 17, 18, 19,
        20, 21, 22, 23, 24, 25, 26
    };

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();
    param = tblite_new_param();

    mol = get_structure_2(error);
    if (tblite_check(error))
        goto err;

    tblite_export_gfn2_param(error, param);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_xtb_calculator(ctx, mol, param);
    if (!calc)
        goto err;

    tblite_get_calculator_shell_count(ctx, calc, &nshells);
    if (tblite_check(ctx))
        goto err;

    if (!check(nshells, 27, "dimension error"))
        goto err;

    tblite_get_calculator_shell_map(ctx, calc, sh2at);
    if (tblite_check(ctx))
        goto err;

    if (!check(sh2at, ref_sh2at, nshells, "shell mapping mismatch"))
        goto err;

    tblite_get_calculator_angular_momenta(ctx, calc, am);
    if (tblite_check(ctx))
        goto err;

    if (!check(am, ref_am, nshells, "angular momenta mismatch"))
        goto err;

    tblite_get_calculator_orbital_count(ctx, calc, &nao);
    if (tblite_check(ctx))
        goto err;

    if (!check(nao, 57, "dimension error"))
        goto err;

    tblite_get_calculator_orbital_map(ctx, calc, ao2sh);
    if (tblite_check(ctx))
        goto err;

    if (!check(ao2sh, ref_ao2sh, nao, "orbital mapping mismatch"))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -14.7042098026524, thr, "energy error"))
        goto err;

    tblite_get_result_gradient(error, res, gradient);

    if (!check(norm2(3 * 17, gradient), 0.0283783086422, thr, "gradient error"))
        goto err;

    double xyz_2[3 * 17] = {
        +1.79735083062742,
        -4.42665143396774,
        +0.00000000000000,
        +5.47735160077691,
        -1.95257593092707,
        +0.00000000000000,
        -1.87902462147155,
        -1.89828871441757,
        +0.00000000000000,
        -0.96120535979079,
        +2.44146318756260,
        +0.00000000000000,
        -4.66879760562471,
        +4.87274102808145,
        +0.00000000000000,
        +1.79476797598653,
        -6.08411129731253,
        -2.27178821177621,
        +1.79476797598653,
        -6.08411129731253,
        +2.27178821177621,
        -3.43213700587973,
        -2.48796105064581,
        +2.26922051661074,
        -3.43213700587973,
        -2.48796105064581,
        -2.26922051661074,
        +0.58302096074950,
        +3.05710624847502,
        -2.26863281575249,
        +0.58302096074950,
        +3.05710624847502,
        +2.26863281575249,
        -4.10211039714465,
        +7.62280080181539,
        +0.00000000000000,
        -6.19547097682738,
        +4.27265764500722,
        -2.27794345710076,
        -6.19547097682738,
        +4.27265764500722,
        +2.27794345710076,
        +7.78350996936708,
        -3.55410309152337,
        +0.00000000000000,
        +5.52628183760127,
        -0.31038446883577,
        -2.27576215506289,
        +5.52628183760127,
        -0.31038446883577,
        +2.27576215506289,
    };
    tblite_update_structure_geometry(error, mol, xyz_2, NULL);
    if (tblite_check(error))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -14.7039683607488, thr, "energy error"))
        goto err;

    tblite_get_result_gradient(error, res, gradient);
    if (tblite_check(error))
        goto err;

    if (!check(norm2(3 * 17, gradient), 0.0279191562471, thr, "gradient error"))
        goto err;

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(res);
    tblite_delete(param);
    return 0;

err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(res);
    tblite_delete(param);
    return 1;
}

int test_ipea1_ch4(void)
{
    printf("Start test: IPEA1-xTB\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_result res = NULL;

    const double thr = 5.0e-7;
    int norbs;
    double charge_neutral = 0.0;
    int uhf_neutral = 0;
    double xyz_neutral[15] = {
        0.00000000000000,
        0.00000000000000,
        0.00000000000000,
        1.18771160655551,
        -1.18771160655551,
        1.18771160655551,
        -1.18771160655551,
        1.18771160655551,
        1.18771160655551,
        -1.18771160655551,
        -1.18771160655551,
        -1.18771160655551,
        1.18771160655551,
        1.18771160655551,
        -1.18771160655551,
    };
    double overlap[13 * 13];
    double energy;

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();

    tblite_set_context_verbosity(ctx, 0);
    if (tblite_check(ctx))
        goto err;

    mol = get_structure_3(error, true);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_ipea1_calculator(ctx, mol);
    if (!calc)
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_set_context_verbosity(ctx, 2);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -4.027631971332, thr, "energy error"))
        goto err;

    tblite_update_structure_geometry(error, mol, xyz_neutral, NULL);
    if (tblite_check(error))
        goto err;

    tblite_update_structure_charge(error, mol, &charge_neutral);
    if (tblite_check(error))
        goto err;

    tblite_update_structure_uhf(error, mol, &uhf_neutral);
    if (tblite_check(error))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -4.670465980661, thr, "energy error"))
        goto err;

    tblite_get_result_overlap_matrix(error, res, overlap);
    if (!tblite_check(error))
        goto err;

    show(error);

    tblite_get_result_hamiltonian_matrix(error, res, overlap);
    if (!tblite_check(error))
        goto err;

    show(error);

    tblite_set_calculator_save_integrals(ctx, calc, 1);
    if (tblite_check(ctx))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_number_of_orbitals(error, res, &norbs);
    if (tblite_check(error))
        goto err;

    if (!check(norbs, 13, "dimension error"))
        goto err;

    tblite_get_result_overlap_matrix(error, res, overlap);
    if (tblite_check(error))
        goto err;

    if (!check(overlap[0], 1.0, thr, "overlap error"))
        goto err;

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(res);
    return 0;

err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(res);
    return 1;
}

int test_gfn1_co2(void)
{
    printf("Start test: GFN1-xTB\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_result res = NULL;

    const double thr = 5.0e-7;
    double energy;

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();

    mol = get_structure_4(error);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_gfn1_calculator(ctx, mol);
    if (!calc)
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -46.203659007308, thr, "energy error"))
        goto err;

    tblite_save_result_wavefunction(error, res, ".test-co2-gfn1.npz");
    if (tblite_check(error))
        goto err;

    tblite_delete(res);
    res = tblite_new_result();

    tblite_load_result_wavefunction(error, res, ".test-co2-gfn1.npz");
    if (tblite_check(error))
        goto err;

    tblite_set_calculator_max_iter(ctx, calc, 3);
    if (tblite_check(ctx))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -46.203659007308, thr, "energy error"))
        goto err;

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(res);
    return 0;

err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(res);
    return 1;
}

int test_gfn2_convergence(void)
{
    printf("Start test: GFN2-xTB (convergence)\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_result res = NULL;
    tblite_param param = NULL;

    const double thr = 5.0e-7;
    int natoms = 2;
    int num[2] = { 3, 8 };
    double xyz[3 * 2] = {
        +0.00000000000000,
        +0.00000000000000,
        +1.50105302628963,
        +0.00000000000000,
        +0.00000000000000,
        -1.50105302628963,
    };
    double energy;

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();
    param = tblite_new_param();

    mol = tblite_new_structure(error, natoms, num, xyz, NULL, NULL, NULL, NULL);
    if (tblite_check(error))
        goto err;

    tblite_export_gfn2_param(error, param);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_xtb_calculator(ctx, mol, param);
    if (!calc)
        goto err;

    tblite_set_calculator_guess(ctx, calc, TBLITE_GUESS_EEQ);
    if (tblite_check(ctx))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -4.228326553369, thr, "energy error"))
        goto err;

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(res);
    tblite_delete(param);
    return 0;

err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(res);
    tblite_delete(param);
    return 1;
}

int test_spgfn1()
{
    printf("Start test: spGFN1-xTB\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_container cont = NULL;
    tblite_result res = NULL;

    const double thr = 5.0e-7;

    double energy;

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();

    mol = get_structure_5(error);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_gfn1_calculator(ctx, mol);
    if (!calc)
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -28.349613833732931, thr, "energy error"))
        goto err;

    cont = tblite_new_spin_polarization(ctx, mol, calc, 1.0);
    if (tblite_check(ctx))
        goto err;

    tblite_calculator_push_back(ctx, calc, &cont);
    if (tblite_check(ctx))
        goto err;

    if (!!cont)
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    int nspin;
    tblite_get_result_number_of_spins(error, res, &nspin);
    if (tblite_check(error))
        goto err;

    if (!check(nspin, 2, "nspin error"))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -28.349613833732931, thr, "energy error"))
        goto err;

    int uhf = 2;
    tblite_update_structure_uhf(error, mol, &uhf);
    if (tblite_check(error))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;
                       
    if (!check(energy, -28.376485853038645, thr, "energy error"))
        goto err;

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 0;

err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 1;
}

int get_array_dimension(int ndim1, int ndim2, int ndim3){
    int array_size = ndim1;

    if (ndim2 == 0){
        return array_size;
    }

    array_size = ndim1 * ndim2;

    if (ndim3 == 0){
        return array_size;
    }

    array_size = ndim1 * ndim2 * ndim3;
    return array_size;

}

int test_dict_api()
{
    printf("Start test: dictionary entries\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_container cont = NULL;
    tblite_result res = NULL;
    tblite_double_dictionary dict = NULL;
    char label[40];

    const double thr = 5.0e-7;

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();

    mol = get_structure_5(error);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_gfn1_calculator(ctx, mol);
    if (!calc)
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    dict = tblite_get_post_processing_dict(error, res);

    int n_dict_entries = 0;
    n_dict_entries = tblite_get_n_entries_dict(error, dict);
    if (!check_int(3, n_dict_entries, "Check number of entries in dict")) {
        goto err;
    }

    int ndim1 = 0;
    int ndim2 = 0;
    int ndim3 = 0;

    tblite_get_array_size_index(error, dict, &n_dict_entries, &ndim1, &ndim2, &ndim3);

    if (!check_int(ndim1, 6, "Check dimension of quadrupole tensor")) {
        goto err;
    }
    if (!check_int(ndim2, 0, "Check dimension of quadrupole tensor")) {
        goto err;
    }
    if (!check_int(ndim3, 0, "Check dimension of quadrupole tensor")) {
        goto err;
    }
    tblite_get_label_entry_index(error, dict, &n_dict_entries, label, NULL);
    if (!(strcmp(label, "molecular-quadrupole") == 0)){
        goto err;
    }

    int ndim1_label = 0;
    int ndim2_label = 0;
    int ndim3_label = 0;

    tblite_get_array_size_label(error, dict, label, &ndim1_label, &ndim2_label, &ndim3_label);

    if (!check_int(ndim1, ndim1_label, "Comparing dimension of quadrupole tensor")) {
        goto err;
    }
    if (!check_int(ndim2, ndim2_label, "Comparing dimension of quadrupole tensor")) {
        goto err;
    }
    if (!check_int(ndim3, ndim3_label, "Comparing dimension of quadrupole tensor")) {
        goto err;
    }

    n_dict_entries = 2;

    tblite_get_array_size_index(error, dict, &n_dict_entries, &ndim1, &ndim2, &ndim3);

    if (!check_int(ndim1, 3, "Check dimension of dipole tensor")) {
        goto err;
    }
    if (!check_int(ndim2, 0, "Check dimension of dipole tensor")) {
        goto err;
    }
    if (!check_int(ndim3, 0, "Check dimension of dipole tensor")) {
        goto err;
    }
    tblite_get_label_entry_index(error, dict, &n_dict_entries, label, NULL);
    if (!(strcmp(label, "molecular-dipole") == 0)){
        goto err;
    }

    tblite_get_array_size_label(error, dict, label, &ndim1_label, &ndim2_label, &ndim3_label);

    if (!check_int(ndim1, ndim1_label, "Check dimension of dipole tensor")) {
        goto err;
    }
    if (!check_int(ndim2, ndim2_label, "Check dimension of dipole tensor")) {
        goto err;
    }
    if (!check_int(ndim3, ndim3_label, "Check dimension of dipole tensor")) {
        goto err;
    }

    n_dict_entries = 1;

    tblite_get_array_size_index(error, dict, &n_dict_entries, &ndim1, &ndim2, &ndim3);

    if (!check_int(ndim1, 21, "Check dimension of wbo tensor")) {
        goto err;
    }
    if (!check_int(ndim2, 21, "Check dimension of wbo tensor")) {
        goto err;
    }
    if (!check_int(ndim3, 1, "Check dimension of wbo tensor")) {
        goto err;
    }
    tblite_get_label_entry_index(error, dict, &n_dict_entries, label, NULL);
    if (!(strcmp(label, "bond-orders") == 0)){
        goto err;
    }

    tblite_get_array_size_label(error, dict, label, &ndim1_label, &ndim2_label, &ndim3_label);

    if (!check_int(ndim1, ndim1_label, "Check dimension of wbo tensor")) {
        goto err;
    }
    if (!check_int(ndim2, ndim2_label, "Check dimension of wbo tensor")) {
        goto err;
    }
    if (!check_int(ndim3, ndim3_label, "Check dimension of wbo tensor")) {
        goto err;
    }

    cont = tblite_new_spin_polarization(ctx, mol, calc, 1.0);
    if (tblite_check(ctx))
        goto err;

    tblite_calculator_push_back(ctx, calc, &cont);
    if (tblite_check(ctx))
        goto err;

    if (!!cont)
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;
    dict= NULL;
    dict = tblite_get_post_processing_dict(error, res);
    if (tblite_check(error))
        goto err;

    ndim1 = 0;
    ndim2 = 0;
    ndim3 = 0;
    int index = 1;
    ndim1_label = 0;
    ndim2_label = 0;
    ndim3_label = 0;

    tblite_get_array_size_index(error, dict, &index, &ndim1, &ndim2, &ndim3);
    if (tblite_check(error))
        goto err;
    tblite_get_label_entry_index(error, dict, &index, label, NULL);
    if (tblite_check(error))
        goto err;

    if (!(strcmp(label, "bond-orders") == 0)){
        goto err;
    }

    tblite_get_array_size_label(error, dict, label, &ndim1_label, &ndim2_label, &ndim3_label);

    if (!check_int(ndim1, ndim1_label, "Check dimension of wbo tensor")) {
        goto err;
    }
    if (!check_int(ndim2, ndim2_label, "Check dimension of wbo tensor")) {
        goto err;
    }
    if (!check_int(ndim3, ndim3_label, "Check dimension of wbo tensor")) {
        goto err;
    }
    int dim_tot = 0;
    dim_tot = get_array_dimension(ndim1, ndim2, ndim3);
    double* wbo_index = (double*) malloc((dim_tot)*sizeof(double));
    double* wbo_label = (double*) malloc((dim_tot)*sizeof(double));

    tblite_get_array_entry_index(error, dict, &index, wbo_index);
    tblite_get_array_entry_label(error, dict, "bond-orders", wbo_label);
    if (!check(wbo_index, wbo_label, thr, dim_tot, "Comparing the wbo arrays retrieved by index and label")){
       goto err;
    }
    free(wbo_index);
    free(wbo_label);
    index = 2;

    tblite_get_array_size_index(error, dict, &index, &ndim1, &ndim2, &ndim3);
    if (tblite_check(error))
        goto err;
    tblite_get_label_entry_index(error, dict, &index, label, NULL);
    if (tblite_check(error))
        goto err;

    if (!(strcmp(label, "molecular-dipole") == 0)){
        goto err;
    }

    tblite_get_array_size_label(error, dict, label, &ndim1_label, &ndim2_label, &ndim3_label);

    if (!check_int(ndim1, ndim1_label, "Check dimension of dipole tensor")) {
        goto err;
    }
    if (!check_int(ndim2, ndim2_label, "Check dimension of dipole tensor")) {
        goto err;
    }
    if (!check_int(ndim3, ndim3_label, "Check dimension of dipole tensor")) {
        goto err;
    }
    dim_tot = get_array_dimension(ndim1, ndim2, ndim3);
    wbo_index = (double*) malloc((dim_tot)*sizeof(double));
    wbo_label = (double*) malloc((dim_tot)*sizeof(double));

    tblite_get_array_entry_index(error, dict, &index, wbo_index);
    tblite_get_array_entry_label(error, dict, "molecular-dipole", wbo_label);
    if (!check(wbo_index, wbo_label, thr, dim_tot, "Comparing the dipole arrays retrieved by index and label")){
       goto err;
    }
    free(wbo_index);
    free(wbo_label);
    int uhf = 2;
    tblite_update_structure_uhf(error, mol, &uhf);
    if (tblite_check(error))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    dict = tblite_get_post_processing_dict(error, res);
    if (tblite_check(error))
        goto err;
    index = 1;
      tblite_get_array_size_index(error, dict, &index, &ndim1, &ndim2, &ndim3);
    if (tblite_check(error))
        goto err;
    tblite_get_label_entry_index(error, dict, &index, label, NULL);
    if (tblite_check(error))
        goto err;

    if (!(strcmp(label, "bond-orders") == 0)){
        goto err;
    }

    tblite_get_array_size_label(error, dict, label, &ndim1_label, &ndim2_label, &ndim3_label);

    if (!check_int(ndim1, ndim1_label, "Check dimension of wbo tensor")) {
        goto err;
    }
    if (!check_int(ndim2, ndim2_label, "Check dimension of wbo tensor")) {
        goto err;
    }
    if (!check_int(ndim3, ndim3_label, "Check dimension of wbo tensor")) {
        goto err;
    }
    dim_tot = 0;
    dim_tot = get_array_dimension(ndim1, ndim2, ndim3);
    wbo_index = (double*) malloc((dim_tot)*sizeof(double));
    wbo_label = (double*) malloc((dim_tot)*sizeof(double));

    tblite_get_array_entry_index(error, dict, &index, wbo_index);
    tblite_get_array_entry_label(error, dict, "bond-orders", wbo_label);
    if (!check(wbo_index, wbo_label, thr, dim_tot, "Comparing the wbo arrays retrieved by index and label")){
       goto err;
    }
    free(wbo_index);
    free(wbo_label);

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    tblite_delete(dict);
    return 0;

    err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    tblite_delete(dict);
    return 1;
}

int test_uninitialized_dict()
{
    printf("Start test: uninitialized dictionary\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_container cont = NULL;
    tblite_result res = NULL;
    tblite_double_dictionary dict = NULL;

    error = tblite_new_error();
    ctx = tblite_new_context();

    mol = get_structure_5(error);
    if (tblite_check(error))
        goto unexpected;

    calc = tblite_new_gfn1_calculator(ctx, mol);
    if (!calc)
        goto unexpected;

    dict = tblite_get_post_processing_dict(error, res);
    if (!tblite_check(error))
        goto unexpected;

    res = tblite_new_result();
    dict = tblite_get_post_processing_dict(error, res);
    if (!tblite_check(error))
        goto unexpected;

    tblite_get_singlepoint(ctx, mol, calc, res);

    error = NULL;
    dict = tblite_get_post_processing_dict(error, res);
    if (!tblite_check(error))
        goto unexpected;
    error = tblite_new_error();

    int n_entries = 0;
    n_entries = tblite_get_n_entries_dict(error, dict);
    if (!check(0, n_entries, "Number of entries of empty dictionary unequal to 0"))
        goto unexpected;
    int dim1 = 0;
    int dim2 = 0;
    int dim3 = 0;
    const int index = 1;

    tblite_get_array_size_index(error, dict, &index, &dim1, &dim2, &dim3);
    if (!check_int(dim1, 0, "Check dimension of empty dict")) {
        goto unexpected;
    }
    if (!check_int(dim2, 0, "Check dimension of empty dict")) {
        goto unexpected;
    }
    if (!check_int(dim3, 0, "Check dimension of empty dict")) {
        goto unexpected;
    }

    double array[10];

    tblite_get_array_entry_index(error, dict, &index, array);

    error = tblite_new_error();
    dict = tblite_get_post_processing_dict(error, res);

    error = NULL;
    n_entries = 0;
    n_entries = tblite_get_n_entries_dict(error, dict);
    if (!check(0, n_entries, "Number of entries of empty dictionary unequal to 0"))
        goto unexpected;
    dim1 = 0;
    dim2 = 0;
    dim3 = 0;

    tblite_get_array_size_index(error, dict, &index, &dim1, &dim2, &dim3);
    if (!tblite_check(error)){
        goto unexpected;
    }
    if (!check_int(dim1, 0, "Check dimension of empty dict")) {
        goto unexpected;
    }
    if (!check_int(dim2, 0, "Check dimension of empty dict")) {
        goto unexpected;
    }
    if (!check_int(dim3, 0, "Check dimension of empty dict")) {
        goto unexpected;
    }

    tblite_get_array_entry_index(error, dict, &index, array);
    if (!tblite_check(error)){
        goto unexpected;
    }
    tblite_get_array_size_label(error, dict, "wbo", &dim1, &dim2, &dim3);
    if (!tblite_check(error)){
        goto unexpected;
    }
    if (!check_int(dim1, 0, "Check dimension of empty dict")) {
        goto unexpected;
    }
    if (!check_int(dim2, 0, "Check dimension of empty dict")) {
        goto unexpected;
    }
    if (!check_int(dim3, 0, "Check dimension of empty dict")) {
        goto unexpected;
    }

    tblite_get_array_entry_label(error, dict, "wbo", array);
    if (!tblite_check(error)){
        goto unexpected;
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);

    return 0;

    unexpected:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 1;
}

int test_h2plus_wbo()
{
    printf("Test H2+ WBO and molecular multipoles\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_result res = NULL;
    tblite_double_dictionary dict = NULL;

    const double thr = 5.0e-4;

    int natoms = 2;
    int num[2] = { 1, 1 };
    double xyz[2*3] = {
        +0.00000000, +0.000000000, +0.472429040,
        +0.00000000, +0.000000000, -0.472429040,
    };

    error = tblite_new_error();
    ctx = tblite_new_context();
    double charge = +1.0;
    int uhf = 1;

    mol = tblite_new_structure(error, natoms, num, xyz, &charge, &uhf, NULL, NULL);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_gfn2_calculator(ctx, mol);
    if (tblite_check(ctx))
        goto err;

    res = tblite_new_result();

    tblite_push_back_post_processing_str(ctx, calc, "molmom");
    if (tblite_check(ctx))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    dict = tblite_get_post_processing_dict(error, res);
    if (tblite_check(error))
        goto err;

    int n_entries = 0;
    n_entries = tblite_get_n_entries_dict(error, dict);
    if (!check_int(n_entries, 2,"Check number of entries in dict")) {
        goto err;
    }

    double dipm[3] = {0.000000, 0.000000, 0.000000};
    int index = 1;
    int ndim1 = 0;
    int ndim2 = 0;
    int ndim3 = 0;

    tblite_get_array_size_index(error, dict, &index, &ndim1, &ndim2, &ndim3);
    if (!check_int(ndim1, 3, "Check dimension of dipole tensor")) {
        goto err;
    }
    if (!check_int(ndim2, 0, "Check dimension of dipole tensor")) {
        goto err;
    }
    if (!check_int(ndim3, 0, "Check dimension of dipole tensor")) {
        goto err;
    }
    double* dipm_dict = (double*) malloc((get_array_dimension(ndim1, ndim2, ndim3))*sizeof(double));

    tblite_get_array_entry_index(error, dict, &index, dipm_dict);
    double norm2_comp ;
    double norm2_exp ;
    norm2_comp = norm2(ndim1, dipm_dict);
    norm2_exp = norm2(ndim1, dipm);
    if (!check_double(norm2_comp, norm2_exp, thr, "Norm2 of mol. dipole not mathcing with expected")) {
        free(dipm_dict);
        goto err;
    }
    free(dipm_dict);

    double qp[6] = {-0.1519, 0.0000, -0.1519, 0.0000, 0.0000, 0.3038};
    index = 2;

    tblite_get_array_size_index(error, dict, &index, &ndim1, &ndim2, &ndim3);
    if (!check_int(ndim1, 6, "Check dimension of quadrupole tensor")) {
        goto err;
    }
    if (!check_int(ndim2, 0, "Check dimension of quadrupole tensor")) {
        goto err;
    }
    if (!check_int(ndim3, 0, "Check dimension of quadrupole tensor")) {
        goto err;
    }
    double* qp_dict = (double*) malloc((get_array_dimension(ndim1, ndim2, ndim3))*sizeof(double));

    tblite_get_array_entry_index(error, dict, &index, qp_dict);

    norm2_comp = norm2(ndim1, qp_dict);
    norm2_exp = norm2(ndim1, qp);
    if (!check_double(norm2_comp, norm2_exp, thr, "Norm2 of mol. dipole not mathcing with expected")) {
        free(qp_dict);
        goto err;
    }
    free(qp_dict);

    calc = NULL;

    calc = tblite_new_gfn2_calculator(ctx, mol);
    if (tblite_check(ctx))
        goto err;

    res = tblite_new_result();

    tblite_push_back_post_processing_str(ctx, calc, "bond-orders");
    if (tblite_check(ctx))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    dict = tblite_get_post_processing_dict(error, res);
    if (tblite_check(error))
        goto err;

    n_entries = 0;
    n_entries = tblite_get_n_entries_dict(error, dict);
    if (!check_int(n_entries, 1, "Check number of entries in dict")) {
        goto err;
    }

    double wbo[4] = {0.000000000000000E+000, 0.500000000000000, 0.500000000000000,
  0.000000000000000E+000};
    index = 1;
    tblite_get_array_size_index(error, dict, &index, &ndim1, &ndim2, &ndim3);

    double* wbo_comp = (double*) malloc((get_array_dimension(ndim1, ndim2, ndim3))*sizeof(double));
    tblite_get_array_entry_label(error, dict, "bond-orders", wbo_comp);
    if (!check(wbo_comp, wbo, thr, natoms*natoms, "WBO for H2+ are not as expected should be 0.5!")){
        goto err;
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(res);
    tblite_delete(dict);
    return 0;
err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(res);
    tblite_delete(dict);
    return 1;
}

int my_strcmp(const char *a, const char *b)
{
  for(; *a && *b && *a == *b; ++a, ++b)
    ;
  if(*a < *b)
    return -1;
  return *a > *b;
}

int test_post_processing_api()
{
    printf("Start test: Post Processing\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_container cont = NULL;
    tblite_result res = NULL;
    tblite_double_dictionary dict = NULL;
    tblite_table table = NULL;
    tblite_param param = NULL;
    tblite_table table_extra = NULL;
    tblite_table bottom = NULL;
    const double thr = 5.0e-7;

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();

    mol = get_structure_5(error);
    if (tblite_check(error))
        goto err;

    tblite_push_back_post_processing_str(ctx, calc, "bond-orders");
    if (!tblite_check(ctx))
        goto err;

    show(ctx);

    calc = tblite_new_gfn1_calculator(ctx, mol);
    if (!calc)
        goto err;

    tblite_push_back_post_processing_str(ctx, calc, "bond-orders");

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    dict = tblite_get_post_processing_dict(error, res);
    if (!dict)
        goto err;

    int n_dict_entries = 0;
    n_dict_entries = tblite_get_n_entries_dict(error, dict);
    if (!check_int(n_dict_entries, 1, "Check number of entries in dict, double addition of wbo")) {
        goto err;
    }
    double dipole[3];
    double qp[6];
    tblite_get_result_dipole(error, res, dipole);
    if (!tblite_check(error))
        goto err;
    show(error);
    tblite_get_result_quadrupole(error, res, qp);
    if (!tblite_check(error))
        goto err;
    show(error);
    double dipm_dict[3];
    double qp_dict[6];
    int ndim1;
    int ndim2;
    int ndim3;
    tblite_get_array_size_label(error, dict, "molecular-dipole", &ndim1, &ndim2, &ndim3);
    if (!tblite_check(error))
        goto err;
    show(error);
    tblite_get_array_size_label(error, dict, "molecular-quadrupole", &ndim1, &ndim2, &ndim3);
    if (!tblite_check(error))
        goto err;
    show(error);

    tblite_get_array_entry_label(error, dict, "molecular-dipole", dipm_dict);
    if (!tblite_check(error))
        goto err;
    show(error);
    tblite_get_array_entry_label(error, dict, "molecular-quadrupole", qp_dict);
    if (!tblite_check(error))
        goto err;
    show(error);

    dict = NULL;

    tblite_push_back_post_processing_str(ctx, calc, "molmom");
    // Test that adding another wbo container will lead to nothing, i.e. there will be 3 entries in the dictionary
    tblite_push_back_post_processing_str(ctx, calc, "bond-orders");

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    dict = tblite_get_post_processing_dict(error, res);
    if (tblite_check(error))
        goto err;
    if (!dict)
        goto err;

    n_dict_entries = 0;
    n_dict_entries = tblite_get_n_entries_dict(error, dict);
    if (tblite_check(error))
        goto err;
    if (!check_int(n_dict_entries, 3, "Check number of entries in dict")) {
        goto err;
    }

    ndim1 = 0;
    ndim2 = 0;
    ndim3 = 0;

    tblite_get_array_size_index(error, dict, &n_dict_entries, &ndim1, &ndim2, &ndim3);
    if (tblite_check(error))
        goto err;
    if (!check_int(ndim1, 6, "Check dimension of quadrupole tensor")) {
        goto err;
    }
    if (!check_int(ndim2, 0, "Check dimension of quadrupole tensor")) {
        goto err;
    }
    if (!check_int(ndim3, 0, "Check dimension of quadrupole tensor")) {
        goto err;
    }

    n_dict_entries = 2;

    tblite_get_array_size_index(error, dict, &n_dict_entries, &ndim1, &ndim2, &ndim3);
    if (tblite_check(error))
        goto err;
    if (!check_int(ndim1, 3, "Check dimension of dipole tensor")) {
        goto err;
    }
    if (!check_int(ndim2, 0, "Check dimension of dipole tensor")) {
        goto err;
    }
    if (!check_int(ndim3, 0, "Check dimension of dipole tensor")) {
        goto err;
    }

    n_dict_entries = 1;

    tblite_get_array_size_index(error, dict, &n_dict_entries, &ndim1, &ndim2, &ndim3);
    if (tblite_check(error))
        goto err;

    if (!check_int(ndim1, 21, "Check dimension of wbo tensor")) {
        goto err;
    }
    if (!check_int(ndim2, 21, "Check dimension of wbo tensor")) {
        goto err;
    }
    if (!check_int(ndim3, 1, "Check dimension of wbo tensor")) {
        goto err;
    }

    n_dict_entries = 4;

    tblite_get_array_size_index(error, dict, &n_dict_entries, &ndim1, &ndim2, &ndim3);
    if (!tblite_check(error))
        goto err;
    show(error);

    if (!check_int(ndim1, 0, "Check dimension of wbo tensor")) {
        goto err;
    }
    if (!check_int(ndim2, 0, "Check dimension of wbo tensor")) {
        goto err;
    }
    if (!check_int(ndim3, 0, "Check dimension of wbo tensor")) {
        goto err;
    }

    tblite_get_array_entry_label(error, dict, "molecular-dipole", dipm_dict);
    tblite_get_array_entry_label(error, dict, "molecular-quadrupole", qp_dict);
    tblite_get_result_dipole(error, res, dipole);
    if (tblite_check(error))
        goto err;
    tblite_get_result_quadrupole(error, res, qp);
    if (tblite_check(error))
        goto err;

    if (!check(dipm_dict, dipole, thr, 3, "Dipole from dict lookup and API are unequal!")){
        goto err;
    }
    if (!check(qp_dict, qp, thr, 6, "Dipole from dict lookup and API are unequal!")){
        goto err;
    }

    bottom = tblite_new_table(NULL);
    table_extra = tblite_new_table(NULL);
    table = tblite_new_table(NULL);
    param = tblite_new_param();
    tblite_export_gfn1_param(error, param);
    tblite_dump_param(error, param, table);

    table_extra = tblite_table_add_table(error, table, "post-processing");
    if (tblite_check(error))
        goto err;
    bottom = tblite_table_add_table(error, table_extra, "molecular-multipole");
    if (tblite_check(error))
        goto err;
    bool bo = true;
    tblite_table_set_bool(error, bottom, "dipole", &bo, 0);
    if (tblite_check(error))
        goto err;
    tblite_delete(param);
    tblite_push_back_post_processing_param(ctx, calc, param);
    if (!tblite_check(ctx))
        goto err;
    show(ctx);

    param = tblite_new_param();
    tblite_load_param(error, param, table);
    table = tblite_new_table(NULL);
    tblite_dump_param(error, param, table);
     if (tblite_check(error))
        goto err;
    tblite_delete(calc);
    calc = NULL;

    tblite_push_back_post_processing_param(ctx, calc, param);
    if (!tblite_check(ctx))
        goto err;
    show(ctx);

    calc = tblite_new_gfn1_calculator(ctx, mol);
    if (!calc)
        goto err;

    tblite_push_back_post_processing_param(ctx, calc, param);
    if (tblite_check(ctx))
        goto err;

    tblite_delete(res);
    res = tblite_new_result();
    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;
    dict = NULL;
    dict = tblite_get_post_processing_dict(error, res);
    if (tblite_check(error))
        goto err;
    if (!dict)
        goto err;

    n_dict_entries = 0;
    n_dict_entries = tblite_get_n_entries_dict(error, dict);
    if (tblite_check(error))
        goto err;
    if (!check_int(n_dict_entries, 1, "Check number of entries in dict, using param for push_back")) {
        goto err;
    }
    calc = tblite_new_gfn1_calculator(ctx, mol);
    cont = tblite_new_spin_polarization(ctx, mol, calc, 1.0);
    if (tblite_check(ctx))
        goto err;

    tblite_calculator_push_back(ctx, calc, &cont);
    if (tblite_check(ctx))
        goto err;

    if (!!cont)
        goto err;
    tblite_set_context_verbosity(ctx, 2);
    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;
    dict = NULL;
    dict = tblite_get_post_processing_dict(error, res);
    if (tblite_check(error))
        goto err;
    if (!dict)
        goto err;

    int index = 1;
    tblite_get_array_size_index(error, dict, &index, &ndim1, &ndim2, &ndim3);

    double* wbo_index = (double*) malloc((get_array_dimension(ndim1, ndim2, ndim3))*sizeof(double));
    tblite_get_array_entry_index(error, dict, &index, wbo_index);
    tblite_get_array_size_label(error, dict, "bond-orders", &ndim1, &ndim2, &ndim3);
    double* wbo_label = (double*) malloc((get_array_dimension(ndim1, ndim2, ndim3))*sizeof(double));
    tblite_get_array_entry_label(error, dict, "bond-orders", wbo_label);
    if (!check(wbo_label, wbo_index, thr, get_array_dimension(ndim1, ndim2, ndim3), "WBO for spin-polarized hamiltonian are not equal!")){
        goto err;
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    tblite_delete(dict);
    tblite_delete(param);
    tblite_delete(table);
    tblite_delete(bottom);
    tblite_delete(table_extra);
    return 0;

    err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    tblite_delete(dict);
    tblite_delete(param);
    tblite_delete(table);
    tblite_delete(bottom);
    tblite_delete(table_extra);
    return 1;
}

int test_xtbml_api()
{
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_container cont = NULL;
    tblite_result res = NULL;
    tblite_double_dictionary dict = NULL;
    tblite_table table = NULL;
    tblite_table table_extra = NULL;
    tblite_table bottom = NULL;
    tblite_param param = NULL;

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();

    mol = get_structure_5(error);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_gfn1_calculator(ctx, mol);
    if (!calc)
        goto err;

    tblite_push_back_post_processing_str(ctx, calc, "xtbml");

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    dict = tblite_get_post_processing_dict(error, res);
    if (!dict)
        goto err;

    int n_dict_entries = 0;
    n_dict_entries = tblite_get_n_entries_dict(error, dict);
    if (!check_int(n_dict_entries, 38, "Check number of entries in dict, double addition of wbo")) {
        goto err;
    }

    double CN_array[21];
    tblite_get_array_entry_label(error, dict, "CN_A", CN_array);
    if (tblite_check(error))
        goto err;


    bottom = tblite_new_table(NULL);
    table_extra = tblite_new_table(NULL);
    table = tblite_new_table(NULL);
    param = tblite_new_param();
    
    tblite_export_gfn1_param(error, param);
    if (tblite_check(error))
        goto err;
    tblite_dump_param(error, param, table);

    table_extra = tblite_table_add_table(error, table, "post-processing"); 
    if (tblite_check(error))
        goto err;
    bottom = tblite_table_add_table(error, table_extra, "xtbml");
    if (tblite_check(error))
        goto err;
    bool bo = true;
    tblite_table_set_bool(error, bottom, "geometry", &bo, 0);
    tblite_table_set_bool(error, bottom, "density", &bo, 0);
    tblite_table_set_bool(error, bottom, "orbital", &bo, 0);
    tblite_table_set_bool(error, bottom, "energy", &bo, 0);
    tblite_table_set_bool(error, bottom, "convolution", &bo, 0);
    tblite_table_set_bool(error, bottom, "tensorial-output", &bo, 0);

    tblite_delete(param);
    param = tblite_new_param();
    tblite_load_param(error, param, table);
    tblite_delete(calc);
    tblite_delete(dict);

    calc = tblite_new_gfn1_calculator(ctx, mol);
    if (!calc)
        goto err;

    tblite_push_back_post_processing_param(ctx, calc, param);
    if (tblite_check(ctx))
        goto err;

    tblite_delete(res);
    res = tblite_new_result();
    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    dict = tblite_get_post_processing_dict(error, res);
    if (!dict)
        goto err;

    n_dict_entries = 0;
    n_dict_entries = tblite_get_n_entries_dict(error, dict);
    if (!check_int(n_dict_entries, 101, "Check number of entries in dict, double addition of wbo")) {
        goto err;
    }
    
    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    tblite_delete(dict);
    tblite_delete(param);
    tblite_delete(table);
    tblite_delete(bottom);
    tblite_delete(table_extra);
    return 0;
    err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    tblite_delete(dict);
    tblite_delete(param);
    tblite_delete(table);
    return 1;
}

int test_uninitialized_solvation()
{
    printf("Start test: Uninitialized Solvation\n");
    tblite_error error = NULL;
    tblite_structure mol = NULL;
    tblite_container cont = NULL;

    error = tblite_new_error();
    // check that we get an immediate return if error is not associated
    cont = tblite_new_gb_solvation_epsilon(error, mol, 0.0, 1, 1);
    if (cont != NULL) goto unexpected;
    error = tblite_new_error();

    //check if it fails when mol is not associated
    cont = tblite_new_cpcm_solvation_epsilon(error, mol, 0.0);
    if (!tblite_check(error))
        goto unexpected;
    show(error);

    tblite_delete(error);
    tblite_delete(mol);
    tblite_delete(cont);
    return 0;

unexpected:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(mol);
    tblite_delete(cont);
    return 1;
}

int test_solvation_cpcm_eps()
{
    printf("Start test: CPCM Solvation\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_container cont = NULL;
    tblite_result res = NULL;
    const double thr = 1.0e-6;

    double energy;

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();
    mol = get_structure_5(error);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_gfn2_calculator(ctx, mol);
    if (!calc)
        goto err;

    cont = tblite_new_cpcm_solvation_epsilon(error, mol, 7.0);
    if (tblite_check(error))
        goto err;

    tblite_calculator_push_back(ctx, calc, &cont);
    if (tblite_check(ctx))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -28.43287176929, thr, "energy error"))
        goto err;

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 0;
err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 1;
}

int test_solvation_alpb_eps()
{
    printf("Start test: ALPB Solvation\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_container cont = NULL;
    tblite_result res = NULL;
    const double thr = 1.0e-6;

    double energy;

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();
    mol = get_structure_5(error);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_gfn2_calculator(ctx, mol);
    if (!calc)
        goto err;

    cont = tblite_new_gb_solvation_epsilon(error, mol, 7.0, tblite_solvation_gbe, tblite_born_p16);
    if (tblite_check(error))
        goto err;

    tblite_calculator_push_back(ctx, calc, &cont);
    if (tblite_check(ctx))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -28.43674134364, thr, "energy error"))
        goto err;

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 0;
err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 1;
}

int test_solvation_alpb_gfn2()
{
    printf("Start test: ALPB(GFN2) Solvation\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_container cont = NULL;
    tblite_result res = NULL;
    const double thr = 1.0e-6;

    double energy;

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();
    mol = get_structure_5(error);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_gfn2_calculator(ctx, mol);
    if (!calc)
        goto err;

    cont = tblite_new_alpb_solvation_solvent(error, mol, "ethanol", tblite_solvation_alpb_gfn2, tblite_state_bar1mol);
    if (tblite_check(error))
        goto err;

    tblite_calculator_push_back(ctx, calc, &cont);
    if (tblite_check(ctx))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -28.4455122, thr, "energy error"))
        goto err;

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 0;
err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 1;
}

int test_solvation_gbsa_gfn2()
{
    printf("Start test: GBSA(GFN2) Solvation\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_container cont = NULL;
    tblite_result res = NULL;
    const double thr = 1.0e-6;

    double energy;

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();
    mol = get_structure_5(error);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_gfn2_calculator(ctx, mol);
    if (!calc)
        goto err;

    cont = tblite_new_alpb_solvation_solvent(error, mol, "water", tblite_solvation_gbsa_gfn2, tblite_state_gsolv);
    if (tblite_check(error))
        goto err;

    tblite_calculator_push_back(ctx, calc, &cont);
    if (tblite_check(ctx))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -28.439916755536, thr, "energy error"))
        goto err;

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 0;
err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 1;
}

int test_solvation_gb_eps()
{
    printf("Start test: GB Solvation\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_container cont = NULL;
    tblite_result res = NULL;
    const double thr = 1.0e-6;

    double energy;

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();
    mol = get_structure_5(error);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_gfn2_calculator(ctx, mol);
    if (!calc)
        goto err;

    cont = tblite_new_gb_solvation_epsilon(error, mol, 7.0, tblite_solvation_gb, tblite_born_still);
    if (tblite_check(error))
        goto err;

    tblite_calculator_push_back(ctx, calc, &cont);
    if (tblite_check(ctx))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -28.4367683, thr, "energy error"))
        goto err;

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 0;
err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 1;
}

int test_solvation_alpb_gfn1()
{
    printf("Start test: ALPB(GFN1) Solvation\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_container cont = NULL;
    tblite_result res = NULL;
    const double thr = 1.0e-6;

    double energy;

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();
    mol = get_structure_5(error);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_gfn1_calculator(ctx, mol);
    if (!calc)
        goto err;

    cont = tblite_new_alpb_solvation_solvent(error, mol, "ethanol", tblite_solvation_alpb_gfn1, tblite_state_bar1mol);
    if (tblite_check(error))
        goto err;

    tblite_calculator_push_back(ctx, calc, &cont);
    if (tblite_check(ctx))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -28.3572423, thr, "energy error"))
        goto err;

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 0;
err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 1;
}

int test_solvation_gbsa_gfn1()
{
    printf("Start test: GBSA(GFN1) Solvation\n");
    tblite_error error = NULL;
    tblite_context ctx = NULL;
    tblite_structure mol = NULL;
    tblite_calculator calc = NULL;
    tblite_container cont = NULL;
    tblite_result res = NULL;
    const double thr = 1.0e-6;

    double energy;

    error = tblite_new_error();
    ctx = tblite_new_context();
    res = tblite_new_result();
    mol = get_structure_5(error);
    if (tblite_check(error))
        goto err;

    calc = tblite_new_gfn1_calculator(ctx, mol);
    if (!calc)
        goto err;

    cont = tblite_new_alpb_solvation_solvent(error, mol, "thf", tblite_solvation_gbsa_gfn1, tblite_state_reference);
    if (tblite_check(error))
        goto err;

    tblite_calculator_push_back(ctx, calc, &cont);
    if (tblite_check(ctx))
        goto err;

    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check(ctx))
        goto err;

    tblite_get_result_energy(error, res, &energy);
    if (tblite_check(error))
        goto err;

    if (!check(energy, -28.3593652, thr, "energy error"))
        goto err;

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 0;
err:
    if (tblite_check(error)) {
        char message[512];
        tblite_get_error(error, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    if (tblite_check(ctx)) {
        char message[512];
        tblite_get_context_error(ctx, message, NULL);
        printf("[Fatal] %s\n", message);
    }

    tblite_delete(error);
    tblite_delete(ctx);
    tblite_delete(mol);
    tblite_delete(calc);
    tblite_delete(cont);
    tblite_delete(res);
    return 1;
}

int main(void)
{
    int stat = 0;
    stat += test_version();
    stat += test_uninitialized_error();
    stat += test_uninitialized_context();
    stat += test_uninitialized_structure();
    stat += test_uninitialized_result();
    stat += test_uninitialized_calculator();
    stat += test_uninitialized_table();
    stat += test_uninitialized_param();
    stat += test_uninitialized_dict();
    stat += test_uninitialized_solvation();
    stat += test_empty_result();
    stat += test_invalid_structure();
    stat += test_table_builder();
    stat += test_param_load();
    stat += test_gfn2_si5h12();
    stat += test_ipea1_ch4();
    stat += test_gfn1_co2();
    stat += test_gfn2_convergence();
    stat += test_spgfn1();
    stat += test_calc_restart();
    stat += test_callback();
    stat += test_dict_api();
    stat += test_post_processing_api();
    stat += test_h2plus_wbo();
    stat += test_solvation_gb_eps();
    stat += test_solvation_alpb_eps();
    stat += test_solvation_cpcm_eps();
    stat += test_solvation_gbsa_gfn2();
    stat += test_solvation_gbsa_gfn1();
    stat += test_solvation_alpb_gfn2();
    stat += test_solvation_alpb_gfn1();
    stat += test_xtbml_api();
    printf("Test finished with %d errors.\n", stat);
    return stat > 0 ? 1 : 0;
}
