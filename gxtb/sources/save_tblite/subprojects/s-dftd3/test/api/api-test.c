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
#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "dftd3.h"

static inline void
show_error(dftd3_error error)
{
   char message[512];
   dftd3_get_error(error, message, NULL);
   printf("[Message] %s\n", message);
}

static inline dftd3_structure
get_test_structure(dftd3_error error)
{
   int const natoms = 7;
   int const attyp[7] = {6,6,6,1,1,1,1};
   double const coord[21] =
      {0.00000000000000, 0.00000000000000,-1.79755622305860,
       0.00000000000000, 0.00000000000000, 0.95338756106749,
       0.00000000000000, 0.00000000000000, 3.22281255790261,
      -0.96412815539807,-1.66991895015711,-2.53624948351102,
      -0.96412815539807, 1.66991895015711,-2.53624948351102,
       1.92825631079613, 0.00000000000000,-2.53624948351102,
       0.00000000000000, 0.00000000000000, 5.23010455462158};
   return dftd3_new_structure(error, natoms, attyp, coord, NULL, NULL);
}

int
test_version (void)
{
   printf("Start test: version\n");
   return dftd3_get_version() > 0 ? 0 : 1;
}

int
test_uninitialized_error (void)
{
   printf("Start test: uninitialized error\n");
   dftd3_error error = NULL;
   return dftd3_check_error(error) ? 0 : 1;
}

int
test_uninitialized_structure (void)
{
   printf("Start test: uninitialized structure\n");
   dftd3_error error = NULL;
   dftd3_structure mol = NULL;

   error = dftd3_new_error();

   double xyz[6] = {0.0};
   dftd3_update_structure(error, mol, xyz, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_delete(error);
   return 0;

unexpected:
   printf("[Fatal] Unexpected pass for uninitialized-structure test\n");
   dftd3_delete(error);
   return 1;
}

int
test_uninitialized_param (void)
{
   printf("Start test: uninitialized parameters\n");
   dftd3_error error = NULL;
   dftd3_structure mol = NULL;
   dftd3_model disp = NULL;
   dftd3_param param = NULL;
   double energy;

   error = dftd3_new_error();

   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_get_pairwise_dispersion(error, mol, disp, param, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   mol = get_test_structure(error);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_get_pairwise_dispersion(error, mol, disp, param, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   disp = dftd3_new_d3_model(error, mol);

   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_get_pairwise_dispersion(error, mol, disp, param, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_delete(error);
   dftd3_delete(disp);
   dftd3_delete(mol);
   return 0;

unexpected:
   printf("[Fatal] Unexpected pass for uninitialized-parameters test\n");
   dftd3_delete(error);
   dftd3_delete(disp);
   dftd3_delete(mol);
   return 1;
}

int
test_uninitialized_model (void)
{
   printf("Start test: uninitialized model\n");
   dftd3_error error = NULL;
   dftd3_model disp = NULL;

   error = dftd3_new_error();

   dftd3_set_model_realspace_cutoff(error, disp, 0.0, 0.0, 0.0);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_delete(error);
   dftd3_delete(disp);
   return 0;

unexpected:
   printf("[Fatal] Unexpected pass for uninitialized-model test\n");
   dftd3_delete(error);
   dftd3_delete(disp);
   return 1;
}

int
test_uninitialized_gcp (void)
{
   printf("Start test: uninitialized counter-poise\n");
   dftd3_error error = NULL;
   dftd3_structure mol = NULL;
   dftd3_gcp gcp = NULL;
   double energy;

   error = dftd3_new_error();

   dftd3_set_gcp_realspace_cutoff(error, gcp, 0.0, 0.0);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_get_counterpoise(error, mol, gcp, &energy, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   mol = get_test_structure(error);

   dftd3_get_counterpoise(error, mol, gcp, &energy, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_delete(error);
   dftd3_delete(mol);
   return 0;

unexpected:
   printf("[Fatal] Unexpected pass for uninitialized-counter-poise test\n");
   dftd3_delete(error);
   dftd3_delete(mol);
   return 1;
}

int
test_invalid_structure (void)
{
   printf("Start test: invalid structure\n");
   dftd3_error error = NULL;
   dftd3_structure mol = NULL;

   int natoms = 2;
   int num[2] = {1, 1};
   double xyz[6] = {0.0};

   error = dftd3_new_error();

   mol = dftd3_new_structure(error, natoms, num, xyz, NULL, NULL);
   if (!dftd3_check_error(error)) goto unexpected;

   show_error(error);

   dftd3_delete(error);
   dftd3_delete(mol);
   return 0;

unexpected:
   printf("[Fatal] Unexpected pass for invalid-structure test\n");
   dftd3_delete(error);
   dftd3_delete(mol);
   return 1;
}

int
test_d3 (void) {
   double energy;
   double pair_disp2[49];
   double pair_disp3[49];
   double gradient[21];
   double sigma[9];

   dftd3_error error;
   dftd3_structure mol;
   dftd3_model disp;
   dftd3_param param;

   error = dftd3_new_error();
   assert(!!error);

   mol = get_test_structure(error);
   if (dftd3_check_error(error)) {return 1;};
   assert(!!mol);

   disp = dftd3_new_d3_model(error, mol);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!disp);

   // PBE-D3(BJ)
   param = dftd3_new_rational_damping(error, 1.0, 0.7875, 0.0, 0.4289, 4.4407, 14.0);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!param);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_pairwise_dispersion(error, mol, disp, param, pair_disp2, pair_disp3);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(param);

   // RPBE-D3(0)
   param = dftd3_load_zero_damping(error, "rpbe", false);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!param);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(param);

   dftd3_set_model_realspace_cutoff(error, disp, 50.0, 30.0, 25.0);
   if (dftd3_check_error(error)) {return 1;}

   // DSD-BLYP-D3(BJ)-ATM
   param = dftd3_load_rational_damping(error, "dsdblyp", true);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!param);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(param);

   // BLYP-D3(0)-ATM
   param = dftd3_new_zero_damping(error, 1.0, 1.682, 1.0, 1.094, 1.0, 14.0);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!param);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(param);

   // B3LYP-D3(CSO)
   param = dftd3_new_cso_damping(error, 1.0, 0.0, 0.86, 2.5, 0.0, 6.25, 14.0);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!param);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_pairwise_dispersion(error, mol, disp, param, pair_disp2, pair_disp3);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(param);

   // PBE-D3(CSO)-ATM
   param = dftd3_load_cso_damping(error, "pbe", true);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!param);
   dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(param);

   dftd3_delete(disp);
   dftd3_delete(mol);
   dftd3_delete(error);

   assert(!param);
   assert(!disp);
   assert(!mol);
   assert(!error);

   return 0;
}

int
test_gcp (void) {
   double energy;
   double pair_disp2[49];
   double pair_disp3[49];
   double gradient[21];
   double sigma[9];

   dftd3_error error;
   dftd3_structure mol;
   dftd3_gcp gcp;

   error = dftd3_new_error();
   assert(!!error);

   mol = get_test_structure(error);
   if (dftd3_check_error(error)) {return 1;};
   assert(!!mol);

   gcp = dftd3_load_gcp_param(error, mol, "pbeh3c", NULL);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!gcp);

   dftd3_get_counterpoise(error, mol, gcp, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_counterpoise(error, mol, gcp, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(gcp);

   gcp = dftd3_load_gcp_param(error, mol, "b973c", NULL);
   if (dftd3_check_error(error)) {return 1;}
   assert(!!gcp);

   dftd3_set_gcp_realspace_cutoff(error, gcp, 50.0, 30.0);
   if (dftd3_check_error(error)) {return 1;}

   dftd3_get_counterpoise(error, mol, gcp, &energy, NULL, NULL);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_get_counterpoise(error, mol, gcp, &energy, gradient, sigma);
   if (dftd3_check_error(error)) {return 1;}
   dftd3_delete(gcp);

   dftd3_delete(mol);
   dftd3_delete(error);

   assert(!gcp);
   assert(!mol);
   assert(!error);

   return 0;
}

int
main (void)
{
   int stat = 0;
   stat += test_version();
   stat += test_uninitialized_error();
   stat += test_uninitialized_structure();
   stat += test_uninitialized_model();
   stat += test_uninitialized_param();
   stat += test_uninitialized_gcp();
   stat += test_invalid_structure();
   stat += test_d3();
   stat += test_gcp();
   return stat;
}
