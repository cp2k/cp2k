/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#ifndef MA_TOPOLOGY_CRAY_H
#define MA_TOPOLOGY_CRAY_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#if defined __SEASTAR || __GEMINI

int extract_topology();
int remove_topology();
#endif
#endif
