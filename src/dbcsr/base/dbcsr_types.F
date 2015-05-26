!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
  USE btree_i8_k_cp2d_v,               ONLY: btree_2d_c => btree
  USE btree_i8_k_dp2d_v,               ONLY: btree_2d_d => btree
  USE btree_i8_k_sp2d_v,               ONLY: btree_2d_s => btree
  USE btree_i8_k_zp2d_v,               ONLY: btree_2d_z => btree
  USE dbcsr_data_types,                ONLY: &
       dbcsr_data_area_type, dbcsr_data_obj, dbcsr_datatype_sizeof, &
       dbcsr_memtype_default, dbcsr_memtype_type, dbcsr_scalar_type, &
       dbcsr_type_complex_4, dbcsr_type_complex_4_2d, dbcsr_type_complex_8, &
       dbcsr_type_complex_8_2d, dbcsr_type_complex_default, dbcsr_type_int_4, &
       dbcsr_type_real_4, dbcsr_type_real_4_2d, dbcsr_type_real_8, &
       dbcsr_type_real_8_2d, dbcsr_type_real_default
  USE kinds,                           ONLY: default_string_length,&
                                             int_8

  !$ USE OMP_LIB

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_types'

  PUBLIC :: dbcsr_obj,&
            dbcsr_type,&
            dbcsr_scalar_type,&
            dbcsr_data_obj,&
            dbcsr_data_area_type,&
            dbcsr_work_type,&
            dbcsr_2d_array_type,&
            dbcsr_array_type,&
            dbcsr_mp_obj,&
            dbcsr_distribution_obj,&
            dbcsr_imagedistribution_type,&
            dbcsr_imagedistribution_obj,&
            dbcsr_iterator,&
            dbcsr_mutable_obj,&
            dbcsr_obj_type_p
  PUBLIC :: dbcsr_meta_size
  PUBLIC :: dbcsr_slot_size,&
            dbcsr_slot_row_p,&
            dbcsr_slot_col_i,&
            dbcsr_slot_blk_p,&
            dbcsr_slot_thr_c,&
            dbcsr_slot_coo_l,&
            dbcsr_slot_nblks,&
            dbcsr_slot_nze,&
            dbcsr_slot_dense,&
            dbcsr_slot_nblkrows_total,&
            dbcsr_slot_nblkcols_total,&
            dbcsr_slot_nfullrows_total,&
            dbcsr_slot_nfullcols_total,&
            dbcsr_slot_nblkrows_local,&
            dbcsr_slot_nblkcols_local,&
            dbcsr_slot_nfullrows_local,&
            dbcsr_slot_nfullcols_local,&
            dbcsr_slot_type,&
            dbcsr_slot_home_prow,&
            dbcsr_slot_home_pcol,&
            dbcsr_slot_home_rowi,&
            dbcsr_slot_home_coli,&
            dbcsr_slot_home_vprow,&
            dbcsr_slot_home_vpcol,&
            dbcsr_num_slots

  PUBLIC :: dbcsr_mpi_size_limits

  PUBLIC :: dbcsr_type_real_4, dbcsr_type_real_8,&
            dbcsr_type_complex_4, dbcsr_type_complex_8,&
            dbcsr_type_real_default, dbcsr_type_complex_default,&
            dbcsr_type_real_4_2d, dbcsr_type_real_8_2d,&
            dbcsr_type_complex_4_2d, dbcsr_type_complex_8_2d,&
            dbcsr_type_int_4
  PUBLIC :: dbcsr_datatype_sizeof
  PUBLIC :: dbcsr_memtype_type,&
            dbcsr_memtype_default

  PUBLIC :: dbcsr_type_invalid, dbcsr_type_no_symmetry, dbcsr_type_symmetric,&
            dbcsr_type_antisymmetric, dbcsr_type_hermitian, dbcsr_type_antihermitian
  PUBLIC :: dbcsr_no_transpose, dbcsr_transpose, dbcsr_conjugate_transpose
  PUBLIC :: dbcsr_repl_none, dbcsr_repl_row, dbcsr_repl_col, dbcsr_repl_full

  PUBLIC :: dbcsr_filter_frobenius
  PUBLIC :: dbcsr_norm_frobenius, dbcsr_norm_maxabsnorm,&
            dbcsr_norm_gershgorin, dbcsr_norm_column

  PUBLIC :: dbcsr_func_inverse, dbcsr_func_tanh, dbcsr_func_dtanh,&
            dbcsr_func_ddtanh, dbcsr_func_artanh, dbcsr_func_inverse_special,&
            dbcsr_func_spread_from_zero,&
            dbcsr_func_sin,&
            dbcsr_func_dsin,&
            dbcsr_func_ddsin,&
            dbcsr_func_asin,&
            dbcsr_func_truncate,&
            dbcsr_func_cos

  PUBLIC :: dbcsr_mpi_statistics_type

! *****************************************************************************
!> \brief A processor (process) grid distribution
!> \var myproc                my processor (process) number
!> \var numproc               number of processors (processes)
!> \var myprow, mypcol        my process grid row and column
!> \var mp_group              message-passing group ID
!> \var pgrid                 processor grid
!> \var refcount              reference counter
!> \var subgroups_defined     whether the subgroups are defined
!> \var prow_group            per-process-row communicator
!> \var pcol_group            pre-process-column communicator
! *****************************************************************************
  TYPE dbcsr_mp_type
     INTEGER                            :: mynode, numnodes, myprow, mypcol
     INTEGER                            :: mp_group
     INTEGER, DIMENSION(:,:), POINTER   :: pgrid => Null()
     INTEGER                            :: refcount
     LOGICAL                            :: subgroups_defined
     INTEGER                            :: prow_group, pcol_group
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
     TYPE(dbcsr_mp_type), POINTER :: mp => Null()
  END TYPE dbcsr_mp_obj

! *****************************************************************************
!> \brief Matrix distribution on the processor grid
!> \var row_dist, col_dist    standard row and column distributions of matrix
!>                            elements' rows and columns into processor grid
!>                            rows and columns
!> \var local_rows            list of rows local to the processor grid row
!> \var local_cols            list of columns local to the processor grid
!>                            column
!> \var other_l_rows          local rows for each process row
!> \var other_l_cols          local columns for each process column
!> \var has_other_l_rows      other_rows is defined
!> \var has_other_l_cols      other_cols is defined
!> \var global_row_map        mapping from rows to sequence in local rows
!>                            (global to local mapping)
!> \var global_col_map        mapping from rows to sequence in local columns
!>                            (global to local mapping)
!> \var has_global_row_map    whether other_row_map is defined
!> \var has_global_col_map    whether other_col_map is defined
!> \var row_map               distribution map for rows
!> \var col_map               distribution map for columns
!> \var thread_dist           thread distribution (of the rows)
!> \var numthreads            number of threads in the environment
!> \var mp_env                mulitprocessor environment on which the
!>                            distribution is based
!> \var refcount              reference counter
! *****************************************************************************
  TYPE dbcsr_distribution_type
     TYPE(array_i1d_obj)               :: row_dist_block, col_dist_block
     TYPE(array_i1d_obj)               :: local_rows, local_cols
     TYPE(array_i1d_obj)               :: row_dist_cluster, col_dist_cluster
     TYPE(array_i1d_obj), DIMENSION(:), POINTER :: local_cluster_rows => Null(),&
                                                   local_cluster_cols => Null()
     TYPE(array_i1d_obj), DIMENSION(:), POINTER :: other_l_rows => Null(),&
                                                   other_l_cols => Null()
     LOGICAL                           :: has_other_l_rows, has_other_l_cols
     TYPE(array_i1d_obj)               :: global_row_map, global_col_map
     LOGICAL                           :: has_global_row_map, has_global_col_map
     TYPE(array_i1d_obj)               :: row_map, col_map
     LOGICAL                           :: has_thread_dist
     TYPE(array_i1d_obj)               :: thread_dist
     INTEGER                           :: num_threads
     TYPE(dbcsr_mp_obj)                :: mp_env
     INTEGER                           :: refcount
  END TYPE dbcsr_distribution_type
  TYPE dbcsr_distribution_obj
     TYPE(dbcsr_distribution_type), POINTER  :: d => Null()
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
!> \var refcount                        count of references
! *****************************************************************************
  TYPE dbcsr_imagedistribution_type
     TYPE(dbcsr_distribution_obj)       :: main
     TYPE(array_i1d_obj)                :: row_image, col_image
     INTEGER                            :: row_decimation, col_decimation
     INTEGER                            :: row_multiplicity, col_multiplicity
     !
     TYPE(array_i1d_obj)                :: vrow_dist, vcol_dist
     !
     TYPE(array_i1d_obj), DIMENSION(:), POINTER :: other_vl_rows => Null(),&
                                                   other_vl_cols => Null()
     TYPE(array_i1d_obj)               :: global_vrow_map, global_vcol_map
     LOGICAL                           :: has_other_vl_rows, has_other_vl_cols,&
                                          has_global_vrow_map, has_global_vcol_map
     INTEGER :: id
     INTEGER :: refcount
  END TYPE dbcsr_imagedistribution_type

  TYPE dbcsr_imagedistribution_obj
     TYPE(dbcsr_imagedistribution_type), POINTER :: i => Null()
  END TYPE dbcsr_imagedistribution_obj

  ! Different method for dbcsr_filter
  INTEGER, PARAMETER            :: dbcsr_filter_frobenius = 1

  ! Different norm for dbcsr_norm
  INTEGER, PARAMETER            :: dbcsr_norm_frobenius = 1
  INTEGER, PARAMETER            :: dbcsr_norm_maxabsnorm = 2
  INTEGER, PARAMETER            :: dbcsr_norm_gershgorin = 3
  INTEGER, PARAMETER            :: dbcsr_norm_column = 4

! *****************************************************************************
!> \brief Buffer for blocks
!> \var refcount        Reference counter
!> \var dirty           Whether any buffers are dirty
!> \var buffers         Buffers
!> \var rcb             Row and column and offset and dimensions of data in the
!>                      buffer
!> \var main            Main memory
!> \var backing         Backing memory (in lieu of main memory)
!> \var data_type       Data type used for the buffers
! *****************************************************************************
  TYPE dbcsr_block_buffer_type
     INTEGER                                     :: refcount
     LOGICAL, DIMENSION(:), POINTER              :: dirty => Null()
     TYPE(dbcsr_data_obj), DIMENSION(:), POINTER :: buffers  => Null()
     INTEGER, DIMENSION(:,:), POINTER            :: rcb  => Null()
     TYPE(dbcsr_data_obj)                        :: main
     TYPE(dbcsr_data_obj), DIMENSION(:), POINTER :: backing  => Null()
     INTEGER                                     :: data_type
  END TYPE dbcsr_block_buffer_type

! *****************************************************************************
!> \brief Object for the buffer of blocks
!> \var t               Block buffer
! *****************************************************************************
  TYPE dbcsr_block_buffer_obj
     TYPE(dbcsr_block_buffer_type), POINTER :: b => Null()
  END TYPE dbcsr_block_buffer_obj


! *****************************************************************************
!> \brief An iterator over a DBCSR matrix.
!> \note This is briefly changed to allow being included in the dbcsr_type type
!>       What is written here is what the structure should be and not what it
!>       is.
!> \var matrix           the matrix
!> \var buffer_2d        Buffers for repointing 2d pointers (1 per thread)
!> \var row              Current row (per thread)
!> \var pos              Current position (per thread)
!> \var row_size         Size of current row
!> \var rbs, cbs         Pointer to row and column size array
!> \var roff, coff       Pointer to row and column offset array
!> \var local_indexing   The matrix has local indexing
!> \var contiguous_pointers   Whether pointers to data should be contiguous
!>                            in memory.
!> \var shared           Iterators share matrix
!> \var dynamic          Ignores the thread distribution (FCFS by block)
!> \var dynamic_byrows   Ignores the thread distribution (FCFS by row)
!> \var common_pos       Position when in mixed mode (row or block
!>                       depending in dynamic_byrows
!> \var local_rows       Mapping of local rows to global rows (if local
!>                       indexing is enabled)
!> \var global_rows      Mapping of global rows to local rows (if local
!>                       indexing is enabled)
! *****************************************************************************
  TYPE dbcsr_iterator
     TYPE(dbcsr_obj), POINTER                       :: matrix => Null()
     TYPE(dbcsr_block_buffer_obj)                   :: buffer_2d
     INTEGER                                        :: pos, row, row_size, &
          row_offset
     INTEGER, DIMENSION(:), POINTER :: rbs  => Null(),&
                                       cbs  => Null(),&
                                       roff => Null(),&
                                       coff => Null()
     LOGICAL            :: local_indexing
     LOGICAL            :: contiguous_pointers, transpose, read_only
     LOGICAL            :: shared, dynamic, dynamic_byrows
     INTEGER, POINTER   :: common_pos => Null()
     ! Copies from the matrix.
     INTEGER            :: nblks, nblkrows_total, max_r_size, max_c_size
     INTEGER, DIMENSION(:), POINTER :: row_p       => Null(),&
                                       col_i       => Null(),&
                                       blk_p       => Null(),&
                                       tdist       => Null(),&
                                       local_rows  => Null(),&
                                       global_rows => Null()
     TYPE(dbcsr_data_obj) :: data_area
  END TYPE dbcsr_iterator


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
     INTEGER                   :: data_type
  END TYPE dbcsr_mutable_type

! *****************************************************************************
!> \brief Object for the mutable data type
! *****************************************************************************
  TYPE dbcsr_mutable_obj
     TYPE(dbcsr_mutable_type), POINTER :: m => Null()
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
     INTEGER, DIMENSION(:), POINTER          :: row_i => Null(),&
                                                col_i => Null(),&
                                                blk_p => Null()
     INTEGER        :: lastblk, datasize
     !TYPE(btree_2d_d)        :: tree
     TYPE(dbcsr_mutable_obj)                 :: mutable
  END TYPE dbcsr_work_type

! *****************************************************************************
!> \brief The BCSR sparse matrix type.
!> \par arrays data and index hold the bulk of the data.
!> \var serial_number   a unique number of each created matrix
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
!> \var coo_l       coordinate list (used for direct indexing)
!> \var col_i       the global blocked column number of this block.
!> \var blk_p       the pointer into the data array of this block.
!> \var thr_c       elements/thread for list index
!> \var row_blk_size          sizes (rows in a block) of blocked rows
!> \var col_blk_size          sizes (columns in a block) of blocked columns
!> \var row_blk_offset          row offset (size = nrow+1)
!> \var col_blk_offset          col offset (size = ncol+1)
!> \var max_rbs, max_cbs      maximal row and column sizes
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
!> \var symmetry              matrix has symmetry
!> \var negate_real           symmetry is realized by negating the real part
!> \var negate_imaginary      symmetry is realized by negating complex part
!>                            (i.e., antisymmetric)
!> \var data_type             'r'/'R' for single/double precision real or
!>                            'c'/'C' for single/double precision complex data
!> \var replication_type      multi-process replication used in the matrix
!> \var transpose             matrix is logically transposed
!> \var bscs                  BCS Column instead of BCS Row
!> \var local_indexing        Local indexing of rows instead of global indexing.
!> \var local_rows            Map of local to global rows when local indexing
!>                            is enabled
!> \var local_rows            Map of global to local rows when local indexing
!>                            is enabled
!> \var index_memory_type     memory type for the index
!> \var data_memory_type      memory type for data
!> \var buffers               Block buffers
!> \var refcount              reference count
!> \var work_mutable          uses the mutable data for working and not the
!>                            append-only data
! *****************************************************************************
  TYPE dbcsr_type
     INTEGER                            :: serial_number
     LOGICAL                            :: valid
     LOGICAL                            :: initialized = .FALSE.
     CHARACTER(LEN=default_string_length) :: name
     TYPE(dbcsr_data_obj)               :: data_area
     INTEGER, DIMENSION(:), POINTER     :: index => Null()
     INTEGER, DIMENSION(:), POINTER     :: row_p => Null(),&
                                           col_i => Null(),&
                                           blk_p => Null(),&
                                           thr_c => Null(),&
                                           coo_l => Null()
     TYPE(array_i1d_obj)                   row_blk_size,&
                                           col_blk_size,&
                                           row_blk_offset,&
                                           col_blk_offset,&
                                           local_rows, global_rows,&
                                           local_cols, global_cols
     LOGICAL                               has_local_rows, has_global_rows,&
                                           has_local_cols, has_global_cols
     INTEGER                            :: max_rbs, max_cbs

     INTEGER                             :: sparsity_id
     INTEGER                             :: id_nr ! use in sm_pool

     INTEGER                 :: nblks, nze
     INTEGER                 :: nblkrows_total, nblkcols_total
     INTEGER                 :: nfullrows_total, nfullcols_total
     INTEGER                 :: nblkrows_local, nblkcols_local
     INTEGER                 :: nfullrows_local, nfullcols_local
     INTEGER                 :: data_type
     CHARACTER               :: replication_type
     LOGICAL                 :: symmetry, negate_real, negate_imaginary
     LOGICAL                 :: bcsc, local_indexing, list_indexing
     TYPE(dbcsr_memtype_type):: data_memory_type, index_memory_type
     TYPE(dbcsr_block_buffer_obj) :: buffers
     TYPE(dbcsr_work_type), DIMENSION(:), POINTER :: wms => Null()
     TYPE(dbcsr_distribution_obj)       :: dist
     INTEGER                 :: refcount
     LOGICAL                 :: work_mutable
  END TYPE dbcsr_type


  CHARACTER, PARAMETER        :: dbcsr_type_invalid = '0'
  CHARACTER, PARAMETER        :: dbcsr_type_no_symmetry = 'N'
  CHARACTER, PARAMETER        :: dbcsr_type_symmetric = 'S'
  CHARACTER, PARAMETER        :: dbcsr_type_antisymmetric = 'A'
  CHARACTER, PARAMETER        :: dbcsr_type_hermitian = 'H'
  CHARACTER, PARAMETER        :: dbcsr_type_antihermitian = 'K'

  !
  ! multiply transpositions
  CHARACTER, PARAMETER        :: dbcsr_no_transpose = 'N'
  CHARACTER, PARAMETER        :: dbcsr_transpose = 'T'
  CHARACTER, PARAMETER        :: dbcsr_conjugate_transpose = 'C'

  CHARACTER, PARAMETER        :: dbcsr_repl_none = 'N'
  CHARACTER, PARAMETER        :: dbcsr_repl_row = 'R'
  CHARACTER, PARAMETER        :: dbcsr_repl_col = 'C'
  CHARACTER, PARAMETER        :: dbcsr_repl_full = 'A'

  ! 
  ! Function types 
  INTEGER, PARAMETER          :: dbcsr_func_inverse = 0 
  INTEGER, PARAMETER          :: dbcsr_func_tanh = 1 
  INTEGER, PARAMETER          :: dbcsr_func_dtanh = 2 
  INTEGER, PARAMETER          :: dbcsr_func_ddtanh = 3 
  INTEGER, PARAMETER          :: dbcsr_func_artanh = 4 
  INTEGER, PARAMETER          :: dbcsr_func_inverse_special = 5 
  INTEGER, PARAMETER          :: dbcsr_func_spread_from_zero = 6 
  INTEGER, PARAMETER          :: dbcsr_func_sin = 7  
  INTEGER, PARAMETER          :: dbcsr_func_dsin = 8 
  INTEGER, PARAMETER          :: dbcsr_func_ddsin = 9 
  INTEGER, PARAMETER          :: dbcsr_func_asin = 10 
  INTEGER, PARAMETER          :: dbcsr_func_cos = 11 
  INTEGER, PARAMETER          :: dbcsr_func_truncate = 12 


! *****************************************************************************
!> These specify which array index in the index array is the start of the
!> specified variable. For example, row_p => index(cp_bcsr_slot_row_p).
! *****************************************************************************
  INTEGER, PARAMETER          :: dbcsr_slot_size  = 1 ! Size of the assigned
                                                      ! values in the index
                                                      ! array.
  INTEGER, PARAMETER          :: dbcsr_slot_nblks = 2
  INTEGER, PARAMETER          :: dbcsr_slot_nze   = 3
  INTEGER, PARAMETER          :: dbcsr_slot_dense = 4
  INTEGER, PARAMETER          :: dbcsr_slot_nblkrows_total   = 5
  INTEGER, PARAMETER          :: dbcsr_slot_nblkcols_total   = 6
  INTEGER, PARAMETER          :: dbcsr_slot_nfullrows_total  = 7
  INTEGER, PARAMETER          :: dbcsr_slot_nfullcols_total  = 8
  INTEGER, PARAMETER          :: dbcsr_slot_nblkrows_local   = 9
  INTEGER, PARAMETER          :: dbcsr_slot_nblkcols_local   = 10
  INTEGER, PARAMETER          :: dbcsr_slot_nfullrows_local  = 11
  INTEGER, PARAMETER          :: dbcsr_slot_nfullcols_local  = 12
  INTEGER, PARAMETER          :: dbcsr_slot_type = 13
  INTEGER, PARAMETER          :: dbcsr_slot_home_prow = 14
  INTEGER, PARAMETER          :: dbcsr_slot_home_pcol = 15
  INTEGER, PARAMETER          :: dbcsr_slot_home_rowi = 16
  INTEGER, PARAMETER          :: dbcsr_slot_home_coli = 17
  INTEGER, PARAMETER          :: dbcsr_slot_home_vprow = 18
  INTEGER, PARAMETER          :: dbcsr_slot_home_vpcol = 19
  ! The number of meta fields.  Its value should be the index of the
  ! last slot listed above.
  INTEGER, PARAMETER          :: dbcsr_meta_size = 19
  INTEGER, PARAMETER          :: dbcsr_slot_row_p = dbcsr_meta_size+2
  INTEGER, PARAMETER          :: dbcsr_slot_col_i = dbcsr_meta_size+4
  INTEGER, PARAMETER          :: dbcsr_slot_blk_p = dbcsr_meta_size+6
  INTEGER, PARAMETER          :: dbcsr_slot_thr_c = dbcsr_meta_size+8
  INTEGER, PARAMETER          :: dbcsr_slot_coo_l = dbcsr_meta_size+10
  INTEGER, PARAMETER          :: dbcsr_num_slots  = dbcsr_meta_size+11 ! previous + 1
  
  ! MPI message size limits (in bytes): 128, 8192, 32KB, 128KB, 4MB, 16MB 
  INTEGER(KIND=int_8), DIMENSION(6), PARAMETER :: dbcsr_mpi_size_limits = & 
       (/ 2**7, 2**13, 2**15, 2**17, 2**22, 2**24 /)

! *****************************************************************************
!> \brief The BCSR sparse matrix object.
!> \var m  The DBCSR array.
! *****************************************************************************
  TYPE dbcsr_obj
     TYPE(dbcsr_type) :: m
  END TYPE dbcsr_obj

! *****************************************************************************
!> \brief Pointer to a work matrix.
!> \var w           the work matrix
! *****************************************************************************
  TYPE dbcsr_work_type_p
     TYPE(dbcsr_work_type), POINTER     :: w => Null()
  END TYPE dbcsr_work_type_p

! *****************************************************************************
!> \brief Pointer to a object.
!> \var matrix         the dbcsr_object
! *****************************************************************************
  TYPE dbcsr_obj_type_p
     TYPE(dbcsr_obj), POINTER           :: matrix => Null()
  END TYPE dbcsr_obj_type_p

! *****************************************************************************
!> \brief A 1-D array of DBCSR matrices
!> \var s           the array of matrices
!> \var refcount    reference counter
! *****************************************************************************
  TYPE dbcsr_1d_array_obj
     TYPE(dbcsr_obj), DIMENSION(:), POINTER :: s => Null()
     INTEGER, POINTER                       :: refcount => Null()
  END TYPE dbcsr_1d_array_obj

! *****************************************************************************
!> \brief A 2-D array of DBCSR matrices
!> \var s           the array of matrices
!> \var refcount    reference counter
! *****************************************************************************
  TYPE dbcsr_2d_array_obj
     TYPE(dbcsr_obj), DIMENSION(:,:), POINTER :: s => Null()
     INTEGER, POINTER                         :: refcount => Null()
  END TYPE dbcsr_2d_array_obj

! *****************************************************************************
!> \brief An array of DBCSR matrices
!> \var mats        the matrices
!> \var data_windows, index_windows     RMA windows
!> \var image_dist  image distribution
! *****************************************************************************
  TYPE dbcsr_array_type
     TYPE(dbcsr_obj), DIMENSION(:), POINTER :: mats => Null()
     TYPE(dbcsr_imagedistribution_type) :: image_dist
  END TYPE dbcsr_array_type

! *****************************************************************************
!> \brief A 2-d array of DBCSR matrices
!> \var mats        the matrices
!> \var image_dist  image distribution
! *****************************************************************************
  TYPE dbcsr_2d_array_type
     TYPE(dbcsr_obj), DIMENSION(:,:), POINTER :: mats => Null()
     TYPE(dbcsr_imagedistribution_obj) :: image_dist
  END TYPE dbcsr_2d_array_type

! *****************************************************************************
!> \brief DBCSR MPI statistics
!> \var ...
!> \var ...
! *****************************************************************************
  TYPE dbcsr_mpi_statistics_type
     INTEGER                               :: col_nimages, row_nimages
     INTEGER                               :: nexchanged
     INTEGER                               :: nfiltered
     INTEGER                               :: nprocs
     ! rank 1: 1=right, 2=left
     ! rank 2: 1=total, 2=min, 3=max
     INTEGER(KIND=int_8), DIMENSION(2,3) :: data_size
     INTEGER(KIND=int_8), DIMENSION(2,3) :: index_size
     ! message size breakdown
     INTEGER(KIND=int_8), DIMENSION(2,SIZE(dbcsr_mpi_size_limits)+1) :: data_size_breakdown
  END TYPE dbcsr_mpi_statistics_type
  
END MODULE dbcsr_types
