!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
MODULE machine_internal

#if defined ( __AIX )
#include "machine_aix.f90"

#elif defined ( __IRIX )
#include "machine_irix.f90"

#elif defined ( __G95 )
#include "machine_g95.f90"

#elif defined ( __GFORTRAN )
#include "machine_gfortran.f90"

#elif defined ( __CCE )
#include "machine_cce.f90"

#elif defined ( __INTEL )
#include "machine_intel.f90"

! must be before __PGI because pgf90 defines __PGI
#elif defined ( __XT3 )
#include "machine_xt3.f90"

#elif defined ( __XT5 )
#include "machine_xt5.f90"

#elif defined ( __PGI )
#include "machine_pgi.f90"

#elif defined ( __T3E )
#include "machine_t3e.f90"

#elif defined ( __DEC )
#include "machine_dec.f90"

#elif defined ( __ABSOFT )
#include "machine_absoft.f90"

#elif defined ( __SUN )
#include "machine_sun.f90"

#elif defined ( __NAG )
#include "machine_nag.f90"

#elif defined ( __SX )
#include "machine_sx.f90"

#else
#error "Missing machine-type (e.g. -D__GFORTRAN)."
#endif

END MODULE machine_internal
