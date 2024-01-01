!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
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
   USE message_passing,                 ONLY: mp_cart_type
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

! **************************************************************************************************
!> \brief info on MPI Cartesian grid that is split on MPI subgroups.
!>        info on distribution of matrix rows / columns to different subgroups.
!> \var mp_comm           global communicator
!> \var pdims             dimensions of process grid
!> \var igroup            which subgroup do I belong to
!> \var ngroup            how many groups in total
!> \var split_rowcol      split row or column?
!> \var pgrid_split_size  how many process rows/cols in subgroups
!> \var group_size        group size (how many cores) of subgroups
!> \var mp_comm_group     sub communicator
!> \var ngroup_opt        optimal number of groups (split factor)
!> \var strict_split      if .true., split factor should not be modified
!>                        (2 parameters for current and general settings)
!> \var refcount          lightweight reference counting for communicators
! **************************************************************************************************
   TYPE dbt_tas_split_info
      TYPE(mp_cart_type)    :: mp_comm
      INTEGER, DIMENSION(2) :: pdims = [-1, -1]
      INTEGER               :: igroup = -1
      INTEGER               :: ngroup = -1
      INTEGER               :: split_rowcol = -1
      INTEGER               :: pgrid_split_size = -1
      INTEGER               :: group_size = -1
      TYPE(mp_cart_type)    :: mp_comm_group
      INTEGER, ALLOCATABLE  :: ngroup_opt
      LOGICAL, DIMENSION(2) :: strict_split = [.FALSE., .FALSE.]
      INTEGER, POINTER      :: refcount => NULL()
   END TYPE

   TYPE dbt_tas_distribution_type
      TYPE(dbt_tas_split_info)                       :: info
      TYPE(dbm_distribution_obj)                     :: dbm_dist
      CLASS(dbt_tas_distribution), ALLOCATABLE       :: row_dist
      CLASS(dbt_tas_distribution), ALLOCATABLE       :: col_dist
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:) :: local_rowcols
   END TYPE

! **************************************************************************************************
!> \brief storage for batched matrix multiplication
!> \var store_batched       intermediate replicated matrix
!> \var store_batched_repl  intermediate replicated matrix
!> \var batched_out         whether replicated matrix has been changed in mm...
!>                          and should be copied to actual matrix
! **************************************************************************************************
   TYPE dbt_tas_mm_storage
      TYPE(dbt_tas_type), POINTER :: store_batched => NULL()
      TYPE(dbt_tas_type), POINTER :: store_batched_repl => NULL()
      LOGICAL                     :: batched_out = .FALSE.
      LOGICAL                     :: batched_trans = .FALSE.
      REAL(dp)                    :: batched_beta = 1.0_dp
   END TYPE

! **************************************************************************************************
!> \brief type for tall-and-skinny matrices
!> \var matrix              matrix on subgroup
!> \var nblkrows            total number of rows
!> \var nblkcols            total number of columns
!> \var nblkrowscols_split  nblkrows or nblkcols depending on which is splitted
!> \var nfullrows           total number of full (not blocked) rows
!> \var nfullcols           total number of full (not blocked) columns
!> \var valid               has been created?
!> \var do_batched          state flag for batched multiplication
!> \var mm_storage          storage for batched processing of matrix matrix multiplication.
!> \var has_opt_pgrid       whether pgrid was automatically optimized
! **************************************************************************************************
   TYPE dbt_tas_type
      TYPE(dbt_tas_distribution_type)         :: dist
      CLASS(dbt_tas_rowcol_data), ALLOCATABLE :: row_blk_size
      CLASS(dbt_tas_rowcol_data), ALLOCATABLE :: col_blk_size
      TYPE(dbm_type)                          :: matrix
      INTEGER(KIND=int_8)                     :: nblkrows = -1
      INTEGER(KIND=int_8)                     :: nblkcols = -1
      INTEGER(KIND=int_8)                     :: nblkrowscols_split = -1
      INTEGER(KIND=int_8)                     :: nfullrows = -1
      INTEGER(KIND=int_8)                     :: nfullcols = -1
      LOGICAL                                 :: valid = .FALSE.
      INTEGER                                 :: do_batched = 0
      TYPE(dbt_tas_mm_storage), ALLOCATABLE   :: mm_storage
      LOGICAL                                 :: has_opt_pgrid = .FALSE.
   END TYPE

   TYPE dbt_tas_iterator
      TYPE(dbt_tas_distribution_type), POINTER :: dist => NULL()
      TYPE(dbm_iterator)                       :: iter
   END TYPE dbt_tas_iterator

END MODULE
