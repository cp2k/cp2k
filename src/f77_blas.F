!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief This module is the real poisoning for BLAS. To avoid I/O issues at
!>        compiling time it is better to avoid the inclusion of the f77_blas_extra
!>        and f77_blas_generic into this module and include them directly in the
!>        needed modules. This reduces a lot the compile time on several machines..
!>        See f77_blas_netlib, f77_blas_poison and f77_blas_generic for further
!>        comments.
!> \author Joost VandeVondele
!> \par History
!>      - Teodoro Laino [tlaino] University of Zurich - 10.2008
! *****************************************************************************
MODULE f77_blas
! In the needed module, please includes f77_blas_generic and f77_blas_extra
! USE f77_blas_extra
! USE f77_blas_generic
  USE f77_blas_poison

END MODULE f77_blas
