!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief DBT tall-and-skinny base types.
!>        Mostly wrappers around existing DBM routines.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_tas_types
   USE dbm_api,                         ONLY: dbm_distribution_obj,&
                                              dbm_iterator,&
                                              dbm_type
   USE dbt_tas_global,                  ONLY: dbt_tas_distribution,&
                                              dbt_tas_rowcol_data
   USE kinds,                           ONLY: dp,&
                                              int_8
#include "../../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_tas_types'

   PUBLIC :: &
      dbt_tas_distribution_type, &
      dbt_tas_iterator, &
      dbt_tas_split_info, &
      dbt_tas_type, &
      dbt_tas_mm_storage

   ! info on MPI Cartesian grid that is split on MPI subgroups.
   ! info on distribution of matrix rows / columns to different subgroups.
   TYPE dbt_tas_split_info
      INTEGER :: mp_comm ! global communicator
      INTEGER, DIMENSION(2) :: pdims ! dimensions of process grid
      INTEGER :: igroup ! which subgroup do I belong to
      INTEGER :: ngroup ! how many groups in total
      INTEGER :: split_rowcol ! split row or column?
      INTEGER :: pgrid_split_size ! how many process rows/cols in subgroups
      INTEGER :: group_size ! group size (how many cores) of subgroups
      INTEGER :: mp_comm_group ! sub communicator
      INTEGER, ALLOCATABLE :: ngroup_opt ! optimal number of groups (split factor)
      LOGICAL, DIMENSION(2) :: strict_split = [.FALSE., .FALSE.]
      ! if .true., split factor should not be modified (2 parameters for current and general settings)
      INTEGER, POINTER :: refcount => NULL() ! lightweight reference counting for communicators
   END TYPE

   TYPE dbt_tas_distribution_type
      TYPE(dbt_tas_split_info) :: info
      TYPE(dbm_distribution_obj) :: dbm_dist
      CLASS(dbt_tas_distribution), ALLOCATABLE :: row_dist
      CLASS(dbt_tas_distribution), ALLOCATABLE :: col_dist
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:) :: local_rowcols
   END TYPE

   ! storage for batched matrix multiplication
   TYPE dbt_tas_mm_storage
      TYPE(dbt_tas_type), POINTER :: store_batched => NULL() ! intermediate replicated matrix
      TYPE(dbt_tas_type), POINTER :: store_batched_repl => NULL() ! intermediate replicated matrix
      LOGICAL :: batched_out ! whether replicated matrix has been changed in mm and should be copied to actual matrix
      LOGICAL :: batched_trans
      REAL(dp) :: batched_beta
   END TYPE

   ! type for tall-and-skinny matrices
   TYPE dbt_tas_type
      TYPE(dbt_tas_distribution_type)  :: dist
      CLASS(dbt_tas_rowcol_data), ALLOCATABLE :: row_blk_size
      CLASS(dbt_tas_rowcol_data), ALLOCATABLE :: col_blk_size

      TYPE(dbm_type) :: matrix ! matrix on subgroup
      INTEGER(KIND=int_8) :: nblkrows ! total number of rows
      INTEGER(KIND=int_8) :: nblkcols ! total number of columns
      INTEGER(KIND=int_8) :: nblkrowscols_split ! nblkrows or nblkcols depending on which is splitted
      INTEGER(KIND=int_8) :: nfullrows ! total number of full (not blocked) rows
      INTEGER(KIND=int_8) :: nfullcols ! total number of full (not blocked) columns
      LOGICAL :: valid = .FALSE. ! has been created?

      ! storage and flags for batched matrix multiplication
      INTEGER :: do_batched = 0 ! state flag for batched multiplication
      TYPE(dbt_tas_mm_storage), ALLOCATABLE :: mm_storage ! storage for batched processing of matrix matrix multiplication.
      LOGICAL :: has_opt_pgrid = .FALSE. ! whether pgrid was automatically optimized
   END TYPE

   TYPE dbt_tas_iterator
      TYPE(dbt_tas_split_info) :: info
      TYPE(dbt_tas_distribution_type) :: dist
      TYPE(dbm_iterator) :: iter
   END TYPE dbt_tas_iterator

END MODULE
