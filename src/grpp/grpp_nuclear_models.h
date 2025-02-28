/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: MIT                                              */
/*----------------------------------------------------------------------------*/

/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#ifndef TEST_LIBGRPP_F90_X_NUCLEAR_MODELS_H
#define TEST_LIBGRPP_F90_X_NUCLEAR_MODELS_H

#define FERMI_UNITS_TO_ATOMIC (1.0 / 52917.7210903)

double libgrpp_estimate_nuclear_rms_radius_johnson_1985(int A);

double libgrpp_estimate_nuclear_rms_radius_golovko_2008(int A);

int libgrpp_estimate_fermi_model_parameters(double R_rms, double *c, double *a);

double libgrpp_charge_density_ball(double r, double Z, double R_rms);

double libgrpp_charge_density_gaussian(double r, double Z, double R_rms);

double libgrpp_charge_density_fermi(double r, double Z, double c, double a);

double libgrpp_charge_density_fermi_bubble(double r, double Z, double c,
                                           double a, double k);

double libgrpp_rms_radius_fermi(int Z, double c, double a);

double libgrpp_rms_radius_fermi_bubble(int Z, double c, double a, double k);

double libgrpp_coulomb_potential_point(double r, double Z);

double libgrpp_coulomb_potential_ball(double r, double Z, double R_rms);

double libgrpp_coulomb_potential_gaussian(double r, double Z, double R_rms);

double libgrpp_coulomb_potential_fermi(double r, double Z, double c, double a);

double libgrpp_coulomb_potential_fermi_bubble(double r, double Z, double c,
                                              double a, double k);

#endif // TEST_LIBGRPP_F90_X_NUCLEAR_MODELS_H
