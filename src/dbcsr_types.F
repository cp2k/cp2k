!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR data types
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.95
!>
!> <b>Modification history:</b>
!> - Created 2009-05-12
! *****************************************************************************
MODULE dbcsr_types
  USE array_types,                     ONLY: array_i1d_obj
  USE btree_I8_k_CP2D_v,               ONLY: btree_2d_c => btree
  USE btree_I8_k_DP2D_v,               ONLY: btree_2d_d => btree
  USE btree_I8_k_SP2D_v,               ONLY: btree_2d_s => btree
  USE btree_I8_k_ZP2D_v,               ONLY: btree_2d_z => btree
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             real_4,&
                                             real_8
  USE message_passing,                 ONLY: mp_window_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_types'

  PUBLIC :: dbcsr_obj,&
            dbcsr_type,&
            dbcsr_data_obj,&
            dbcsr_data_area_type,&
            dbcsr_work_type,&
            dbcsr_array_type,&
            dbcsr_mp_obj, dbcsr_mp_type,&
            dbcsr_distribution_obj, dbcsr_distribution_type,&
            dbcsr_imagedistribution_type,&
            dbcsr_iterator,&
            dbcsr_1d_set_obj,&
            dbcsr_2d_set_obj,&
            dbcsr_mutable_obj, dbcsr_mutable_type
  PUBLIC :: dbcsr_meta_size
  PUBLIC :: dbcsr_slot_size,&
            dbcsr_slot_row_p,&
            dbcsr_slot_col_i,&
            dbcsr_slot_blk_p,&
            dbcsr_slot_nblks,&
            dbcsr_slot_nze,&
            dbcsr_num_slots
  PUBLIC :: dbcsr_type_real_4, dbcsr_type_real_8,&
            dbcsr_type_complex_4, dbcsr_type_complex_8,&
            dbcsr_type_real_default, dbcsr_type_complex_default
  PUBLIC :: dbcsr_type_normal, dbcsr_type_symmetric,&
            dbcsr_type_antisymmetric, dbcsr_type_transposed


  INTEGER, PARAMETER :: dbcsr_meta_size = 13 ! The number of meta fields


! *****************************************************************************
!> \brief A processor (process) grid distribution
!> \var myproc                my processor (process) number
!> \var numproc               number of processors (processes)
!> \var myprow, mypcol        my process grid row and column
!> \var mp_group              message-passing group ID
!> \var pgrid                 processor grid
!> \var refcount              reference counter
!> \var numthreads, mythread  number of threads
! *****************************************************************************
  TYPE dbcsr_mp_type
     INTEGER                            :: mynode, numnodes, myprow, mypcol
     INTEGER                            :: mp_group
     INTEGER, DIMENSION(:,:), POINTER   :: pgrid
     INTEGER                            :: numthreads
     INTEGER                            :: refcount
  END TYPE dbcsr_mp_type

! *****************************************************************************
!> \brief A processor (process) grid distribution
!> \var myproc                my processor (process) number
!> \var numproc               number of processors (processes)
!> \var myprow, mypcol        my process grid row and column
!> \var pgrid                 processor grid
!> \var refcount              reference counter
! *****************************************************************************
  TYPE dbcsr_mp_obj
     TYPE(dbcsr_mp_type), POINTER :: mp
  END TYPE dbcsr_mp_obj

! *****************************************************************************
!> \brief Matrix distribution on the processor grid
!> \var row_dist, col_dist    standard row and column distributions of matrix
!>                            elements' rows and columns into processor grid
!>                            rows and columns
!> \var local_rows            list of rows local to the processor grid row
!> \var local_cols            list of columns local to the processor grid
!>                            column
!> \var thread_dist           thread distribution (of the rows)
!> \var mp_env                mulitprocessor environment on which the
!>                            distribution is based
!> \var refcount              reference counter
!> \var dist2d                associated distribution_2d type
! *****************************************************************************
  TYPE dbcsr_distribution_type
     TYPE(array_i1d_obj)               :: row_dist, col_dist
     TYPE(array_i1d_obj)               :: local_rows, local_cols
     TYPE(array_i1d_obj)               :: thread_dist
     TYPE(dbcsr_mp_obj)                :: mp_env
     INTEGER                           :: refcount
  END TYPE dbcsr_distribution_type
  TYPE dbcsr_distribution_obj
     TYPE(dbcsr_distribution_type), POINTER  :: d
  END TYPE dbcsr_distribution_obj

! *****************************************************************************
!> \brief Image distributions are used to map incompatible processor row and
!>        column distributions.
!>
!> Used to ease storage or transfer between two different-sizes
!> sets. For example, if there are 4 real processor rows that are
!> mapped to 8 "virtual" processor rows, then there are two images for
!> every real processor row.
!> \var main                            the main distribution
!> \var row_image, col_image            distribution of matrix elements' rows
!>                                      and columns into image  rows and
!>                                      columns
!> \var row_decimation, col_decimation  Number of imaged rows/columns mapped
!>                                      to a real row/column
!> \var row_multiplicity, col_multiplicity   Number of real rows/columns mapped
!>                                      to a virtual row/column
! *****************************************************************************
  TYPE dbcsr_imagedistribution_type
     TYPE(dbcsr_distribution_obj)       :: main
     TYPE(array_i1d_obj)                :: row_image, col_image
     INTEGER                            :: row_decimation, col_decimation
     INTEGER                            :: row_multiplicity, col_multiplicity
  END TYPE dbcsr_imagedistribution_type


! *****************************************************************************
!> \brief Stores actual data
!> \var r_sp        stores real values in single precision
!> \var r_dp        stores real values in double precision
!> \var c_sp        stores complex values in single precision
!> \var c_dp        stores complex values in double precision
!> \var refcount    reference counter for current structure
!> \var special     these data are specially-allocated
!> \var data_type   which of the data types is actually used
! *****************************************************************************
  TYPE dbcsr_data_area_type
     REAL(KIND=real_4), DIMENSION(:), POINTER    :: r_sp
     REAL(KIND=real_8), DIMENSION(:), POINTER    :: r_dp
     COMPLEX(KIND=real_4), DIMENSION(:), POINTER :: c_sp
     COMPLEX(KIND=real_8), DIMENSION(:), POINTER :: c_dp
     INTEGER                                  :: refcount
     LOGICAL                                  :: special
     CHARACTER                                :: data_type
  END TYPE dbcsr_data_area_type

  TYPE dbcsr_data_obj
     TYPE(dbcsr_data_area_type), POINTER           :: d
  END TYPE dbcsr_data_obj

  CHARACTER, PARAMETER          :: dbcsr_type_real_4 = 'r'
  CHARACTER, PARAMETER          :: dbcsr_type_real_8 = 'R'
  CHARACTER, PARAMETER          :: dbcsr_type_complex_4 = 'c'
  CHARACTER, PARAMETER          :: dbcsr_type_complex_8 = 'C'

#ifndef __SGL
  CHARACTER, PARAMETER          :: dbcsr_type_real_default = 'R'
  CHARACTER, PARAMETER          :: dbcsr_type_complex_default = 'C'
#else
  CHARACTER, PARAMETER          :: dbcsr_type_real_default = 'r'
  CHARACTER, PARAMETER          :: dbcsr_type_complex_default = 'c'
#endif

! *****************************************************************************
!> \brief Data area with random acces reads, insertions, and deletions.
!> \var tree_s, tree_d, tree_c, tree_z  Data types for the various types
!> \var refcount                        Reference counter
!> \var data_type                       The data type that is stored
! *****************************************************************************
  TYPE dbcsr_mutable_type
     TYPE(btree_2d_s)          :: btree_s
     TYPE(btree_2d_d)          :: btree_d
     TYPE(btree_2d_c)          :: btree_c
     TYPE(btree_2d_z)          :: btree_z
     INTEGER                   :: refcount
     CHARACTER                 :: data_type
  END TYPE dbcsr_mutable_type

! *****************************************************************************
!> \brief Object for the mutable data type
! *****************************************************************************
  TYPE dbcsr_mutable_obj
     TYPE(dbcsr_mutable_type), POINTER :: m
  END TYPE dbcsr_mutable_obj

! *****************************************************************************
!> \brief Used for assembling a real matrix.
!> \var m           the final bcsr matrix
!> \var data_area   holds actual values.
!> \var row_i       the row index of all of the blocks.
!> \var col_i       the column index of all of the blocks.
!> \var blk_p       the pointer into the data array of this block.
!> \var lastblk     index of the last block entered into the row_i, col_i, and
!>                  blk_p data structures
!> \var datasize    the actual size of data present in the data element
!> \var tree        tree used to index data blocks (alternative to the row_i,
!>                  col_i, and blk_p indices when index is scattered).
! *****************************************************************************
  TYPE dbcsr_work_type
     TYPE(dbcsr_data_obj)                    :: data_area
     INTEGER, DIMENSION(:), POINTER          :: row_i,&
          &                                     col_i,&
          &                                     blk_p
     INTEGER        :: lastblk, datasize
     !TYPE(btree_2d_d)        :: tree
     TYPE(dbcsr_mutable_obj)                 :: mutable
  END TYPE dbcsr_work_type

! *****************************************************************************
!> \brief The BCSR sparse matrix type.
!> \par arrays data and index hold the bulk of the data.
!> \var initialized whether the matrix is initialized (i.e., it is not random
!>                  garbage from memory)
!> \var valid       whether the matrix is valid (consistent)
!> \var name        name of the matrix
!> \var index(:)    agglomeration of the indices and offsets
!>                            of pointers into this array.
!> \par the pointers row_p, col_i, blk_p point into the index array.
!> \var row_p       points into the col_i and blk_p arrays, each element
!>                  (1:nblkrows_total+1) points to the previous row's last
!>                  element. So each rows has elements
!>                  row_p(row)+1:row_p(row+1).
!> \var col_i       the global blocked column number of this block.
!> \var blk_p       the pointer into the data array of this block.
!> \var row_blk_size          sizes (rows in a block) of blocked rows
!> \var col_blk_size          sizes (columns in a block) of blocked columns
!> \var max_rbs, max_cbs      maximal row and column sizes
!> \var dist2d      distribution_2d used by this matrix
!> \var dist        distribution used by this matrix
!> \var nblks       number of blocks locally present
!> \var nze         number of non-zero elements locally present
!> \var nblkrows_total, nblk_cols_total size of entire matrix in blocked rows
!>                  and blocked columns
!> \var nblkrows_local, nblk_cols_local size of local part of matrix in blocked
!>                  rows and blocked columns
!> \var nfullrows_total, nfull_cols_total size of entire matrix in full rows
!>                  and full columns
!> \var nfullrows_local, nfull_cols_local size of local part of matrix in full
!>                  rows and full columns
!> \var matrix_type           'N' for normal, 'T' for transposed, 'S' for
!>                            symmetric, and 'A' for antisymmetric
!> \var data_type             'r'/'R' for single/double precision real or
!>                            'c'/'C' for single/double precision complex data
!> \var special_memory        uses special (MPI-allocated) memory
!> \var refcount              reference count
!> \var work_mutable          uses the mutable data for working and not the
!>                            append-only data
! *****************************************************************************
  TYPE dbcsr_type
     LOGICAL                            :: valid, initialized
     CHARACTER(LEN=default_string_length) :: name
     !REAL(kind=dp), DIMENSION(:), POINTER :: DATA
     TYPE(dbcsr_data_obj)               :: data_area
     INTEGER, DIMENSION(:), POINTER     :: index
     INTEGER, DIMENSION(:), POINTER     :: row_p,&
          &                                col_i,&
          &                                blk_p
     TYPE(array_i1d_obj)                   row_blk_size,&
          &                                col_blk_size
     INTEGER                            :: max_rbs, max_cbs
     TYPE(distribution_2d_type), POINTER :: dist2d
     INTEGER                 :: nblks, nze
     INTEGER                 :: nblkrows_total, nblkcols_total
     INTEGER                 :: nfullrows_total, nfullcols_total
     INTEGER                 :: nblkrows_local, nblkcols_local
     INTEGER                 :: nfullrows_local, nfullcols_local
     CHARACTER               :: matrix_type, data_type
     LOGICAL                 :: special_memory
     !TYPE(dbcsr_work_type)   :: w
     TYPE(dbcsr_work_type), DIMENSION(:), POINTER :: wms
     TYPE(dbcsr_distribution_obj)       :: dist
     INTEGER                 :: refcount
     LOGICAL                 :: work_mutable
  END TYPE dbcsr_type

  CHARACTER, PARAMETER        :: dbcsr_type_normal = 'N'
  CHARACTER, PARAMETER        :: dbcsr_type_symmetric = 'S'
  CHARACTER, PARAMETER        :: dbcsr_type_antisymmetric = 'A'
  CHARACTER, PARAMETER        :: dbcsr_type_transposed = 'T'

! *****************************************************************************
!> These specify which array index in the index array is the start of the
!> specified variable. For example, row_p => index(cp_bcsr_slot_row_p).
! *****************************************************************************
  INTEGER, PARAMETER          :: dbcsr_slot_size  = 1 ! Size of the assigned
                                                      ! values in the index
                                                      ! array.
  INTEGER, PARAMETER          :: dbcsr_slot_row_p = dbcsr_meta_size+2
  INTEGER, PARAMETER          :: dbcsr_slot_col_i = dbcsr_meta_size+4
  INTEGER, PARAMETER          :: dbcsr_slot_blk_p = dbcsr_meta_size+6
  INTEGER, PARAMETER          :: dbcsr_slot_nblks = 2
  INTEGER, PARAMETER          :: dbcsr_slot_nze   = 3
  INTEGER, PARAMETER          :: dbcsr_num_slots  = dbcsr_meta_size+7


! *****************************************************************************
!> \brief The BCSR sparse matrix object.
!> \var m  The DBCSR array.
! *****************************************************************************
  TYPE dbcsr_obj
     TYPE(dbcsr_type), POINTER :: m => NULL()
  END TYPE dbcsr_obj

! *****************************************************************************
!> \brief Pointer to a work matrix.
!> \var w           the work matrix
! *****************************************************************************
  TYPE dbcsr_work_type_p
     TYPE(dbcsr_work_type), POINTER     :: w
  END TYPE dbcsr_work_type_p

! *****************************************************************************
!> \brief A 1-D set of DBCSR matrices
!> \var s           the set of matrices
!> \var refcount    reference counter
! *****************************************************************************
  TYPE dbcsr_1d_set_obj
     TYPE(dbcsr_obj), DIMENSION(:), POINTER :: s
     INTEGER, POINTER                       :: refcount
  END TYPE dbcsr_1d_set_obj

! *****************************************************************************
!> \brief A 2-D set of DBCSR matrices
!> \var s           the set of matrices
!> \var refcount    reference counter
! *****************************************************************************
  TYPE dbcsr_2d_set_obj
     TYPE(dbcsr_obj), DIMENSION(:,:), POINTER :: s
     INTEGER, POINTER                         :: refcount
  END TYPE dbcsr_2d_set_obj

! *****************************************************************************
!> \brief An array of DBCSR matrices
!> \var mats        the matrices
!> \var data_windows, index_windows     RMA windows
!> \var image_dist  image distribution
! *****************************************************************************
  TYPE dbcsr_array_type
     TYPE(dbcsr_obj), DIMENSION(:), POINTER :: mats
     TYPE(mp_window_type), DIMENSION(:), POINTER :: data_windows, index_windows
     TYPE(dbcsr_imagedistribution_type) :: image_dist
  END TYPE dbcsr_array_type

! *****************************************************************************
!> \brief An iterator over a DBCSR matrix.
!> \var matrix           the matrix
!> \var buffer_r_sp      Real single precision buffer
!> \var buffer_r_dp      Real double precision buffer
!> \var buffer_c_sp      Complex single precision buffer
!> \var buffer_c_dp      Complex double precision buffer
!> \var row              Current row
!> \var pos              Current position
!> \var row_size         Size of current row
!> \var rbs, cbs         Pointer to row and column size array
!> \var contiguous_pointers   Whether pointers to data should be contiguous
!>                            in memory.
! *****************************************************************************
  TYPE dbcsr_iterator
     TYPE(dbcsr_obj)                               :: matrix
     REAL(KIND=real_4), DIMENSION(:,:), POINTER    :: buffer_r_sp
     REAL(KIND=real_8), DIMENSION(:,:), POINTER    :: buffer_r_dp
     COMPLEX(KIND=real_4), DIMENSION(:,:), POINTER :: buffer_c_sp
     COMPLEX(KIND=real_8), DIMENSION(:,:), POINTER :: buffer_c_dp
     INTEGER            :: row, row_size
     INTEGER            :: pos
     INTEGER, DIMENSION(:), POINTER :: rbs, cbs
     LOGICAL            :: contiguous_pointers
  END TYPE dbcsr_iterator

END MODULE dbcsr_types
