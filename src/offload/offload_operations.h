/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2022 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

#ifndef OFFLOAD_OPERATIONS_H
#define OFFLOAD_OPERATIONS_H

#if defined __OFFLOAD_HIP
#include "offload_hip_internal.h"
#endif

#if defined __OFFLOAD_CUDA
#include "offload_cuda_internal.h"
#endif

#endif
