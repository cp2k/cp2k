/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#ifndef PW_CUDA_H
#define PW_CUDA_H
/******************************************************************************
 *  Authors: Benjamin G Levine, Andreas Gloess
 *
 *  2012/05/18                 Refacturing - original files:
 *                              - cuda_tools/cuda_pw_cu.cu
 *
 *****************************************************************************/
#if defined ( __PW_CUDA )

/* Double precision complex procedures */
extern "C" void pw_cuda_cfffg_z_ (const double          *din,
                                        cuDoubleComplex *zout,
                                  const int             *ghatmap,
                                  const int             *npts,
                                  const int              ngpts,
                                  const double           scale);


extern "C" void pw_cuda_sfffc_z_ (const cuDoubleComplex *zin,
                                        double          *dout,
                                  const int             *ghatmap,
                                  const int             *npts,
                                  const int              ngpts,
                                  const int              nmaps,
                                  const double           scale);


extern "C" void pw_cuda_cff_z_   (const double          *din,
                                        cuDoubleComplex *zout,
                                  const int             *npts);


extern "C" void pw_cuda_ffc_z_   (const cuDoubleComplex *zin,
                                        double          *dout,
                                  const int             *npts);


extern "C" void pw_cuda_cf_z_    (const double          *din,
                                        cuDoubleComplex *zout,
                                  const int             *npts);


extern "C" void pw_cuda_fc_z_    (const cuDoubleComplex *zin,
                                        double          *dout,
                                  const int             *npts);


extern "C" void pw_cuda_f_z_     (const cuDoubleComplex *zin,
                                        cuDoubleComplex *zout,
                                  const int              dir,
                                  const int              n,
                                  const int              m);


extern "C" void pw_cuda_fg_z_    (const cuDoubleComplex *zin,
                                        cuDoubleComplex *zout,
                                  const int             *ghatmap,
                                  const int             *npts,
                                  const int              mmax,
                                  const int              ngpts,
                                  const double           scale);


extern "C" void pw_cuda_sf_z_    (const cuDoubleComplex *zin,
                                        cuDoubleComplex *zout,
                                  const int             *ghatmap,
                                  const int             *npts,
                                  const int              mmax,
                                  const int              ngpts,
                                  const int              nmaps,
                                  const double           scale);
#endif
#endif
