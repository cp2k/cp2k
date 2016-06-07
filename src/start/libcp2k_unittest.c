/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2016  CP2K developers group                         *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libcp2k.h"

/**
 * \brief Unit test of the C-interface provided via libcp2k.h
 * \author Ole Schuett
 */

int main(int argc, char** argv){

    printf("Unit test starts ...\n");

    // test cp2k_get_version()
    printf("Testing cp_c_get_version(): ");
    char version_str[100];
    cp2k_get_version(version_str, 100);
    printf("%s.\n", version_str);


    // create simple input file
    const char* inp_fn = "H2.inp";
    FILE* f = fopen(inp_fn, "w");
    fprintf(f,"&FORCE_EVAL\n");
    fprintf(f,"  METHOD Quickstep\n");
    fprintf(f,"  &DFT\n");
    fprintf(f,"    BASIS_SET_FILE_NAME BASIS_SET\n");
    fprintf(f,"    POTENTIAL_FILE_NAME POTENTIAL\n");
    fprintf(f,"    LSD\n");
    fprintf(f,"    &MGRID\n");
    fprintf(f,"      CUTOFF 140\n");
    fprintf(f,"    &END MGRID\n");
    fprintf(f,"    &QS\n");
    fprintf(f,"      EPS_DEFAULT 1.0E-8\n");
    fprintf(f,"    &END QS\n");
    fprintf(f,"    &SCF\n");
    fprintf(f,"      EPS_DIIS 0.1\n");
    fprintf(f,"      EPS_SCF 1.0E-4\n");
    fprintf(f,"      MAX_DIIS 4\n");
    fprintf(f,"      MAX_SCF 3\n");
    fprintf(f,"      SCF_GUESS atomic\n");
    fprintf(f,"    &END SCF\n");
    fprintf(f,"    &XC\n");
    fprintf(f,"      &XC_FUNCTIONAL Pade\n");
    fprintf(f,"      &END XC_FUNCTIONAL\n");
    fprintf(f,"    &END XC\n");
    fprintf(f,"  &END DFT\n");
    fprintf(f,"  &SUBSYS\n");
    fprintf(f,"    &CELL\n");
    fprintf(f,"      ABC 8.0 4.0 4.0\n");
    fprintf(f,"    &END CELL\n");
    fprintf(f,"    &COORD\n");
    fprintf(f,"    H     0.000000  0.000000  0.000000\n");
    fprintf(f,"    H     1.000000  0.000000  0.000000\n");
    fprintf(f,"    &END COORD\n");
    fprintf(f,"    &KIND H\n");
    fprintf(f,"      BASIS_SET DZV-GTH-PADE\n");
    fprintf(f,"      POTENTIAL GTH-PADE-q1\n");
    fprintf(f,"    &END KIND\n");
    fprintf(f,"  &END SUBSYS\n");
    fprintf(f,"&END FORCE_EVAL\n");
    fprintf(f,"&GLOBAL\n");
    fprintf(f,"  PRINT_LEVEL SILENT\n");
    fprintf(f,"  PROJECT libcp2k_unittest_H2\n");
    fprintf(f,"&END GLOBAL\n");
    fclose(f);

    // use input file to create a force environment
    force_env_t force_env;
    cp2k_init();
    cp2k_create_force_env(&force_env, inp_fn, "__STD_OUT__");
    cp2k_calc_energy_force(force_env);

    // check energy
    double energy;
    cp2k_get_potential_energy(force_env, &energy);
    printf("\n ENERGY: %.12f\n",energy);
    if(fabs(-1.118912797546392 - energy) / fabs(energy) > 1e-13){
        printf("Wrong energy\n");
        return(-1);
    }

    // clean up
    cp2k_finalize();
    remove(inp_fn);

    printf("Unit test finished, found no errors\n");
    return(0);
}

//EOF
