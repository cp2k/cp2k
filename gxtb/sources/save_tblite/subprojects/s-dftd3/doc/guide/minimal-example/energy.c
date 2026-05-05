#include <stdio.h>
#include <stdbool.h>

#include "dftd3.h"

int main(void)
{
  dftd3_error error = dftd3_new_error();
  dftd3_structure mol = NULL;
  dftd3_model d3 = NULL;
  dftd3_param param = NULL;

  int nat = 5;
  int num[5] = {6, 1, 1, 1, 1};
  double xyz[15] = {  // coordinates in Bohr
     0.0000000, -0.0000000,  0.0000000,
    -1.1922080,  1.1922080,  1.1922080,
     1.1922080, -1.1922080,  1.1922080,
    -1.1922080, -1.1922080, -1.1922080,
     1.1922080,  1.1922080, -1.1922080};
  mol = dftd3_new_structure(error, nat, num, xyz, NULL, NULL);
  if (dftd3_check_error(error)) goto handle_error;

  char method[5] = "PBE0";
  param = dftd3_load_rational_damping(error, method, false);
  if (dftd3_check_error(error)) goto handle_error;

  d3 = dftd3_new_d3_model(error, mol);
  if (dftd3_check_error(error)) goto handle_error;

  double energy;
  dftd3_get_dispersion(error, mol, d3, param, &energy, NULL, NULL);
  if (dftd3_check_error(error)) goto handle_error;

  printf("Dispersion energy for %s-D3(BJ) is %13.10lf Hartree\n", method, energy);

  dftd3_delete(error);
  dftd3_delete(mol);
  dftd3_delete(d3);
  dftd3_delete(param);
  return 0;

handle_error:
  char msg[512];
  dftd3_get_error(error, msg, NULL);
  printf("Error: %s\n", msg);

  dftd3_delete(error);
  dftd3_delete(mol);
  dftd3_delete(d3);
  dftd3_delete(param);
  return 1;
}
