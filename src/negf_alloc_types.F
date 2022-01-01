!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
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

   PUBLIC :: negf_allocatable_ivector, negf_allocatable_rvector, &
             negf_allocatable_imatrix, negf_allocatable_rmatrix

! **************************************************************************************************
!> \brief Allocatable 1-D integer vector
! **************************************************************************************************
   TYPE negf_allocatable_ivector
      !> allocatable 1-D real vector
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: vector
   END TYPE negf_allocatable_ivector

! **************************************************************************************************
!> \brief Allocatable 1-D real vector
! **************************************************************************************************
   TYPE negf_allocatable_rvector
      !> allocatable 1-D real vector
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: vector
   END TYPE negf_allocatable_rvector

! **************************************************************************************************
!> \brief Allocatable 2-D integer matrix
! **************************************************************************************************
   TYPE negf_allocatable_imatrix
      !> allocatable 2-D integer matrix
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: matrix
   END TYPE negf_allocatable_imatrix

! **************************************************************************************************
!> \brief Allocatable 2-D real matrix
! **************************************************************************************************
   TYPE negf_allocatable_rmatrix
      !> allocatable 2-D real matrix
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: matrix
   END TYPE negf_allocatable_rmatrix

END MODULE negf_alloc_types
