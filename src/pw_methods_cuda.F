!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \note
!> This module contains routines necessary to operate on plan waves on NVIDIA
!> gpus using CUDA.
!> \par History
!>      BGL (06-Mar-2008) : Created
!> \author Benjamin G. Levine
! *****************************************************************************
MODULE pw_methods_cuda
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             sp
  USE pw_types,                        ONLY: RECIPROCALSPACE,&
                                             pw_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: cuda_pw_fft_wrap_c1dr3d, cuda_pw_fft_wrap_r3dc1d

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_methods_cuda'
  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.FALSE.

CONTAINS

! *****************************************************************************
!> \brief perform an fft followed by a gather on the gpu
!> \author Benjamin G Levine
! *****************************************************************************
  SUBROUTINE cuda_pw_fft_wrap_r3dc1d(pw1, pw2, dir, n, scale)
    TYPE(pw_type)                            :: pw1, pw2
    INTEGER                                  :: dir
    INTEGER, DIMENSION(:), POINTER           :: n
    REAL(kind=dp)                            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'cuda_pw_fft_wrap_r3dc1d', &
      routineP = moduleN//':'//routineN

    COMPLEX(kind=sp), DIMENSION(:), POINTER  :: csp
    COMPLEX(kind=sp), DIMENSION(:, :, :), &
      POINTER                                :: csp3
    INTEGER                                  :: ierr, l1, l2, l3, ngpts
    INTEGER, DIMENSION(:), POINTER           :: mapl, mapm, mapn
    INTEGER, DIMENSION(:, :), POINTER        :: ghat
    LOGICAL                                  :: failure
    REAL(kind=sp)                            :: fscale

    failure = .FALSE.
#if defined ( __CUDA )
    ngpts = SIZE ( pw2%pw_grid%gsq  )
    ALLOCATE ( csp3( n(1), n(2), n(3) ), STAT = ierr )
    ALLOCATE ( csp( ngpts ), STAT = ierr )

    l1 = LBOUND(pw1%cr3d,1)
    l2 = LBOUND(pw1%cr3d,2)
    l3 = LBOUND(pw1%cr3d,3)

    csp3(1:,1:,1:)=pw1%cr3d(l1:,l2:,l3:)
    fscale=scale

    mapl => pw2%pw_grid%mapl%pos
    mapm => pw2%pw_grid%mapm%pos
    mapn => pw2%pw_grid%mapn%pos
    ghat => pw2%pw_grid%g_hat

    !fsign=-dir

    CALL pw_fft_wrap_fg_cu(dir, csp3, csp, fscale, n, &
      mapl, mapm, mapn, ngpts, ghat, l1, l2, l3)
    
    pw2%cc(:)=csp(:)

    DEALLOCATE( csp3, STAT = ierr )
    DEALLOCATE( csp, STAT = ierr )
    pw2 % in_space = RECIPROCALSPACE
#endif

  END SUBROUTINE cuda_pw_fft_wrap_r3dc1d


! *****************************************************************************
!> \brief perform an fft followed by a gather on the gpu
!> \author Benjamin G Levine
! *****************************************************************************
  SUBROUTINE cuda_pw_fft_wrap_c1dr3d(pw1, pw2, dir, n, scale)
    TYPE(pw_type)                            :: pw1, pw2
    INTEGER                                  :: dir
    INTEGER, DIMENSION(:), POINTER           :: n
    REAL(kind=dp)                            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'cuda_pw_fft_wrap_c1dr3d', &
      routineP = moduleN//':'//routineN

    COMPLEX(kind=sp), DIMENSION(:), POINTER  :: csp
    COMPLEX(kind=sp), DIMENSION(:, :, :), &
      POINTER                                :: csp3
    INTEGER                                  :: ierr, l1, l2, l3, ngpts
    INTEGER, DIMENSION(:), POINTER           :: mapl, mapm, mapn
    INTEGER, DIMENSION(:, :), POINTER        :: ghat
    LOGICAL                                  :: failure
    REAL(kind=sp)                            :: fscale

    failure = .FALSE.
#if defined ( __CUDA )
    ngpts = SIZE ( pw1%pw_grid%gsq  )
    ALLOCATE ( csp3( n(1), n(2), n(3) ), STAT = ierr )
    ALLOCATE ( csp( ngpts ), STAT = ierr )

    l1 = LBOUND(pw2%cr3d,1)
    l2 = LBOUND(pw2%cr3d,2)
    l3 = LBOUND(pw2%cr3d,3)

    csp(:)=pw1%cc(:)

    fscale=scale

    mapl => pw1%pw_grid%mapl%pos
    mapm => pw1%pw_grid%mapm%pos
    mapn => pw1%pw_grid%mapn%pos
    ghat => pw1%pw_grid%g_hat

    CALL pw_fft_wrap_sf_cu(dir, csp, csp3, fscale, n, &
      mapl, mapm, mapn, ngpts, ghat, l1, l2, l3)
    
    pw2%cr3d(l1:,l2:,l3:)=csp3(1:,1:,1:)

    DEALLOCATE( csp3, STAT = ierr )
    DEALLOCATE( csp, STAT = ierr )
#endif

  END SUBROUTINE cuda_pw_fft_wrap_c1dr3d

END MODULE pw_methods_cuda

