/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2019  CP2K developers group                         *
 *****************************************************************************/

// Declare channels for kernel to kernel communication

/******************************************************************************
 *  Authors: Arjun Ramaswami
 *****************************************************************************/
#ifndef FFT_CHANNELS_H
#define FFT_CHANNELS_H

#pragma OPENCL EXTENSION cl_intel_channels : enable

channel cmplex chaninfft[8] __attribute__((depth(8)));
channel cmplex chanoutfft[8] __attribute__((depth(8)));

channel cmplex chaninfft2[8] __attribute__((depth(8)));
channel cmplex chanoutfft2[8] __attribute__((depth(8)));

channel cmplex chaninfetch[8] __attribute__((depth(8)));
#endif
