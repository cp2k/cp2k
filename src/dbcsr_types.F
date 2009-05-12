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
!> - Created
! *****************************************************************************
MODULE dbcsr_types
  USE array_types,                     ONLY: array_i1d_obj
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_window_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_types'

  PUBLIC :: dbcsr_type,&
            dbcsr_work_type,&
            dbcsr_array_type,&
            dbcsr_mp_obj, dbcsr_mp_type,&
            dbcsr_distribution_obj, dbcsr_distribution_type,&
            dbcsr_rptdistribution_type
  PUBLIC :: dbcsr_meta_size
  PUBLIC :: dbcsr_slot_size,&
            dbcsr_slot_row_p,&
            dbcsr_slot_col_i,&
            dbcsr_slot_blk_p,&
            dbcsr_slot_nblks,&
            dbcsr_slot_nze,&
            dbcsr_num_slots

  INTEGER, PARAMETER :: dbcsr_meta_size = 13 ! The number of meta fields


! *****************************************************************************
!> \brief A processor (process) grid distribution
!> \var myproc                my processor (process) number
!> \var numproc               number of processors (processes)
!> \var myprow, mypcol        my process grid row and column
!> \var mp_group              message-passing group ID
!> \var pgrid                 processor grid
!> \var refcount              reference counter
! *****************************************************************************
  TYPE dbcsr_mp_type
     INTEGER                            :: myproc, numproc, myprow, mypcol
     INTEGER                            :: mp_group
     INTEGER, DIMENSION(:,:), POINTER   :: pgrid
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
!> \var mp_env                mulitprocessor environment on which the
!>                            distribution is based
!> \var refcount              reference counter
! *****************************************************************************
  TYPE dbcsr_distribution_type
     TYPE(array_i1d_obj)               :: row_dist, col_dist
     TYPE(dbcsr_mp_obj)                :: mp_env
     INTEGER                           :: refcount
  END TYPE dbcsr_distribution_type
  TYPE dbcsr_distribution_obj
     TYPE(dbcsr_distribution_type), POINTER  :: d
  END TYPE dbcsr_distribution_obj

! *****************************************************************************
!> \brief A subdivided processor grid distribution: the repetitions of a
!>        distribution
!> \var main                            the main distribution
!> \var nrow_repeats, ncol_repeats      number of row and column repetitions
!> \var row_repetition, col_repetition  distribution of matrix elements' rows
!>                                      and columns into subrows and subcolumns
! *****************************************************************************
  TYPE dbcsr_rptdistribution_type
     TYPE(dbcsr_distribution_obj)       :: main
     INTEGER                            :: nrow_repeats, ncol_repeats
     TYPE(array_i1d_obj)                :: row_repetition, col_repetition
  END TYPE dbcsr_rptdistribution_type

! *****************************************************************************
!> \brief The BCSR sparse matrix type.
!> \par arrays data and index hold the bulk of the data.
!> \var valid       whether the matrix is valid (consistent)
!> \var name        name of the matrix
!> \var data(:)     actual matrix elements.
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
!> \var special_memory        uses special (MPI-allocated) memory
! *****************************************************************************
  TYPE dbcsr_type
     LOGICAL        :: valid
     CHARACTER(LEN=40)        :: name
     REAL(kind=dp), DIMENSION(:), POINTER    :: DATA
     INTEGER, DIMENSION(:), POINTER     :: index
     INTEGER, DIMENSION(:), POINTER     :: row_p,&
          &                                col_i,&
          &                                blk_p,&
          &                                row_blk_size,&
          &                                col_blk_size
     TYPE(distribution_2d_type), POINTER :: dist2d
     TYPE(dbcsr_distribution_obj)        :: dist
     INTEGER                 :: nblks, nze
     INTEGER                 :: nblkrows_total, nblkcols_total
     INTEGER                 :: nfullrows_total, nfullcols_total
     INTEGER                 :: nblkrows_local, nblkcols_local
     INTEGER                 :: nfullrows_local, nfullcols_local
     CHARACTER               :: matrix_type
     LOGICAL                 :: special_memory
  END TYPE dbcsr_type

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
!> \brief Used for assembling a real matrix.
!> \var m           the final bcsr matrix
!> \var data        holds actual values.
!> \var shared_data data is shared with the final matrix
!> \var row_i       the row index of all of the blocks.
!> \var col_i       the column index of all of the blocks.
!> \var blk_p       the pointer into the data array of this block.
!> \var lastblk     index of the last block entered into the row_i, col_i, and
!>                  blk_p data structures
!> \var datasize    the actual size of data present in the data element
!> \var src         source processor
! *****************************************************************************
  TYPE dbcsr_work_type
     REAL(KIND=dp), DIMENSION(:), POINTER    :: DATA
     LOGICAL                            :: shared_data
     INTEGER, DIMENSION(:), POINTER     :: row_i,&
          &                                col_i,&
          &                                blk_p
     INTEGER        :: lastblk, datasize
  END TYPE dbcsr_work_type


! *****************************************************************************
!> \brief Pointer to a DBCSR matrix.
!> \var m           the matrix
! *****************************************************************************
  TYPE dbcsr_type_p
     TYPE(dbcsr_type), POINTER          :: m
  END TYPE dbcsr_type_p

! *****************************************************************************
!> \brief Pointer to a work matrix.
!> \var w           the work matrix
! *****************************************************************************
  TYPE dbcsr_work_type_p
     TYPE(dbcsr_work_type), POINTER     :: w
  END TYPE dbcsr_work_type_p

! *****************************************************************************
!> \brief An array of DBCSR matrices
!> \var mats        the matrices
!> \var data_windows, index_windows     RMA windows
!> \var rdist       distribution repetition
! *****************************************************************************
  TYPE dbcsr_array_type
     TYPE(dbcsr_type), DIMENSION(:), POINTER :: mats
     TYPE(mp_window_type), DIMENSION(:), POINTER :: data_windows, index_windows
     TYPE(dbcsr_rptdistribution_type) :: rdist
  END TYPE dbcsr_array_type



! *****************************************************************************
!> \brief Data used for fetching from remote processes
!> \var remote                metadata about the matrices present on other
!>                            processors. Useful for their indices
!> \var window_index          MPI window of indices for this matrix
!> \var window_data           MPI window of data for this matrix
!> \var fetch_slots           metadata for buffers, index slots as indexed by
!>                            the variables fetch_slot_*
!> \var index                 index for data in buffer
!> \var buffer                storage for actual data that is fetched
!> \var special_mem           buffer and index use specially allocated memory,
!>                            i.e., MPI allocation
!> \var transposed            This fetch is associated only with transposed
!>                            blocks
!> \var extents               allocation tracking for buffer_xfer and
!>                            buffer_use
!> \var unicast_groups        MPI groups containing a single processor, one for
!>                            every processor in the system
!> \var wait        window must be waited upon
!> \var waitfor     the processor for which we must wait (i.e., unlock)
! *****************************************************************************
!  TYPE fetch_env
!     TYPE(cp_para_env_type), POINTER    :: para_env
!     TYPE(cp_bcsr_work_type_p), DIMENSION(:), POINTER  :: remote
!     INTEGER        :: window_meta,&
!          &            window_index,&
!          &            window_data
!     !
!     INTEGER, DIMENSION(:,:), POINTER   :: fetch_slots, extents
!     INTEGER, DIMENSION(:), POINTER     :: index, fs
!     REAL(KIND=dp), DIMENSION(:), POINTER    :: buffer
!     LOGICAL                            :: special_memory, transposed
!     !
!     TYPE(cp_bcsr_type), POINTER        :: m
!     INTEGER, DIMENSION(:,:), POINTER   :: blacs2mpi
!     INTEGER, DIMENSION(:), POINTER     :: unicast_groups
!     LOGICAL                            :: specialmem_buffer, specialmem_index
!     LOGICAL, POINTER                   :: wait
!     LOGICAL, DIMENSION(:), POINTER     :: waitfor
!     !INTEGER, DIMENSION(:)              :: data_windows
!  END TYPE fetch_env
!! -----------------------------------------------------------------------------
!  
!  INTEGER, PARAMETER          :: fetch_slot_index = 1 ! My own index, 0 if
!                                                      ! empty
!  INTEGER, PARAMETER          :: fetch_slot_fetch_info = 2 ! Type of fetch,
!                                                           ! see fetch_info_*
!  INTEGER, PARAMETER          :: fetch_slot_from = 3 ! Source processor
!  INTEGER, PARAMETER          :: fetch_slot_row = 4 ! Row/starting row of data
!  INTEGER, PARAMETER          :: fetch_slot_col = 5 ! Column/starting block of
!                                                    ! data
!  INTEGER, PARAMETER          :: fetch_slot_tr = 6  ! Data is transposed.
!  ! Number of rows, cols, both blocked or regular, depending on
!  ! multiplication type.
!  INTEGER, PARAMETER          :: fetch_slot_full_rows = 7 ! "y" size
!  INTEGER, PARAMETER          :: fetch_slot_full_cols = 8 ! "x" size
!  INTEGER, PARAMETER          :: fetch_slot_blocked_rows = 9
!  INTEGER, PARAMETER          :: fetch_slot_blocked_cols = 10
!  INTEGER, PARAMETER          :: fetch_slot_icols = 11 ! Number of indexed
!                                                       ! (blocked) columns
!                                                       ! (entries in col_i
!                                                       ! and row_p)
!  INTEGER, PARAMETER          :: fetch_slot_extent = 12 ! Corresponding (data)
!                                                        ! extent in buffer
!  INTEGER, PARAMETER          :: fetch_slot_data_begin = 13
!  INTEGER, PARAMETER          :: fetch_slot_data_end = 14
!  INTEGER, PARAMETER          :: fetch_slot_hint = 15 ! Hint for next transfer
!  INTEGER, PARAMETER          :: fetch_slot_offset = 16 ! Offset at source
!                                                        ! processor
!  ! Single or multiple blocks
!  INTEGER, PARAMETER          :: fetch_slots_num = 16
!  !
!  ! Standard fetch (xfer) or already fetched/data present in local buffer (use)
!  INTEGER, PARAMETER          :: fetch_info_normal = 0
!  ! Data is local (xfer)/data is present in local buffer (use)
!  INTEGER, PARAMETER          :: fetch_info_nofetch = 1
!  ! Empty due to sparsity, no fetch, no use
!  INTEGER, PARAMETER          :: fetch_info_empty = 2
!  ! Data is mine, no fetch
!  INTEGER, PARAMETER          :: fetch_info_local = 3
!
!  ! Selection of matrix A or matrix B to simplify some "if" statements.
!  LOGICAL, PARAMETER          :: fetch_matrix_a = .TRUE.
!  LOGICAL, PARAMETER          :: fetch_matrix_b = .FALSE.
!  

END MODULE dbcsr_types
