/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2016  CP2K developers group                         *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "cp2k_c_interface.h"

int main(int argc, char** argv){

    printf("Unit test starts...\n");

    printf("Testing cp_c_get_version(): ");
    char version_str[100];
    cp_c_get_version(version_str, 100);
    printf("%s.\n", version_str);

    printf("Unit test finished, found no errors\n");
    return(0);
}

//EOF
