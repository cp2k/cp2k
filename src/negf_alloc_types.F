!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Allocatable vectors for NEGF based quantum transport calculations.
! **************************************************************************************************

MODULE negf_alloc_types
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_alloc_types'

   PUBLIC :: negf_allocatable_ivector, negf_allocatable_rvector

! **************************************************************************************************
!> \brief Allocatable 1-D integer vector
! **************************************************************************************************
   TYPE negf_allocatable_ivector
      !> allocatable 1-D integer vector
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: vector
   END TYPE negf_allocatable_ivector

! **************************************************************************************************
!> \brief Allocatable 1-D real vector
! **************************************************************************************************
   TYPE negf_allocatable_rvector
      !> allocatable 1-D real vector
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: vector
   END TYPE negf_allocatable_rvector

END MODULE negf_alloc_types
