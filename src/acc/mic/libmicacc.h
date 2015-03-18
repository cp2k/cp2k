/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015 the CP2K developers group                      *
 *****************************************************************************/
/* Hans Pabst (Intel Corp.)
******************************************************************************/
#ifndef LIBMICACC_H
#define LIBMICACC_H

#if defined(__ACC) && defined(__ACC_MIC)

#include "../include/acc.h"
#include <libxstream.h>

#endif // defined(__ACC) && defined(__ACC_MIC)
#endif // LIBMICACC_H
