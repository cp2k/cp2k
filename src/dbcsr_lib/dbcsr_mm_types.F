!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Datatypes used by the dbcsr matrix-matrix multiplication machinery.
!>	
!> \author  Urban Borstnik
!>
!> <b>Modification history:</b>
!>  - 2013-01    reorganized code (Ole Schuett)
! *****************************************************************************
MODULE dbcsr_mm_types

  USE dbcsr_kinds,                     ONLY: dp,&
                                             int_1,&
                                             int_4,&
                                             int_8,&
                                             real_8,&
                                             sp
  USE dbcsr_pq_types,                  ONLY: dbcsr_pq_type,&
                                             dbcsr_ps_set_type
  USE dbcsr_types,                     ONLY: dbcsr_work_type

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_types'

  PUBLIC :: carrier_type, hash_table_type, ele_type

! *****************************************************************************
!> \brief Types needed for the hashtable.
! *****************************************************************************
  TYPE ele_type
     INTEGER :: c=0
     INTEGER :: p=0
  END TYPE ele_type
  
  TYPE hash_table_type
     TYPE(ele_type), DIMENSION(:), POINTER :: table
     INTEGER :: nele=0
     INTEGER :: nmax=0
     INTEGER :: prime=0
  END TYPE hash_table_type

  
! *****************************************************************************
!> \brief Used to carry data among the various calls.  Each thread has
!>        its own private copy.
!> \var id                Unique ID of each carrier_type instance
!> \var c_has_symmetry    The product matrix has symmetry
!> \var keep_sparsity     Sparsity of C matrix should be kept
!> \var use_eps           Use on-the-fly filtering
!> \var param_sets        Set of parameter stacks
!> \var local_indexing    The A and B matrix indices and C matrix work indices
!>                        are not global but local to process rows and columns
!> \var m_sizes           Block sizes of A and C matrix rows, indexed locally
!> \var n_sizes           Block sizes of B and C matrix columns, indexed locally
!> \var k_sizes           Block sizes of A matrix columns and B matrix rows,
!>                        indexed locally
!> \var c_local_rows      C and A matrix local rows.  Map from
!>                        local row (index) to global row (value).
!> \var c_local_cols      C and B matrix local columns.  Map from local column
!>                        (index) to global column (value).
!> \var k_locals          A matrix local columns and B matrix local rows.  Map
!>                        from local row/column (index) to global row/column
!>                        (value).
!> \var c_global_rows     C and A matrix global rows.  Map from global rows
!>                        (index) to local rows (value).
!> \var c_global_cols     C and B matrix global columns.  Map from global
!>                        columns (index) to local columns (value).
!> \var max_m             Maximum size of C or A matrix row block size.
!> \var max_n             Maximum size of C or B matrix column block size.
!> \var max_k             Maximum size of A matrix colum and B matrix row block
!>                        size.
!> \var m_size_maps       Map from an A matrix row block size to an index of
!>                        the most common A matrix row block sizes.
!> \var n_size_maps       Map from a B matrix column block size to an index of
!>                        the most common B matrix column block sizes.
!> \var k_size_maps       Map from an A matrix row or B matrix column block
!>                        size to an index of the most row/column block sizes.
!> \var m_size_maps_size  Size of the m_size_maps array (i.e., the maximum block
!>                        size + 1).
!> \var n_size_maps_size  Size of the n_size_maps array (i.e., the maximum block
!>                        size + 1).
!> \var k_size_maps_size  Size of the k_size_maps array (i.e., the maximum block
!>                        size + 1).
!> \var nm_stacks         The number of most common m blocks sizes.
!> \var nn_stacks         The number of most common n blocks sizes.
!> \var nk_stacks         The number of most common k blocks sizes.
!> \var stack_map         Map from most common (n,k,m) block sizes to a stack
!>                        number within a stack group.
!> \var default_stack     Stack to use when stack members do not have common
!>                        m, n, and k sizes.
!> \var row_max_epss      Maximum eps to be used for one row.
!> \var a_norms           Norms of A matrix blocks.
!> \var b_norms           Norms of B matrix blocks.
!> \var my_wm             Work matrix associated with this thread
!> \var lastblk           Number of elements in the work matrix
!>
!> \var datasize          Data size of the work matrix (view of
!>                        same work matrix variable).
!> \var original_lastblk  Number of work matrix blocks before addition
!> \var flop              flop count
!> \var t_index           Time for indexing
!> \var t_gemm            Calculation time
! *****************************************************************************
  TYPE carrier_type
     INTEGER :: id
     LOGICAL :: c_has_symmetry, keep_sparsity, use_eps
     LOGICAL :: local_indexing
     TYPE(hash_table_type), DIMENSION(:), POINTER  :: c_hashes
     TYPE(dbcsr_pq_type), POINTER                  :: queue
     TYPE(dbcsr_ps_set_type)                       :: param_sets
     INTEGER, DIMENSION(:), POINTER :: m_sizes, n_sizes, k_sizes
     INTEGER, DIMENSION(:), POINTER :: c_local_rows, c_local_cols, k_locals,&
                                       c_global_rows, c_global_cols
     INTEGER                        :: max_m, max_n, max_k
     INTEGER(KIND=int_4), DIMENSION(:), POINTER :: m_size_maps,&
                                                   n_size_maps,&
                                                   k_size_maps
     INTEGER                        :: m_size_maps_size,&
                                       n_size_maps_size,&
                                       k_size_maps_size
     INTEGER                        :: nm_stacks, nn_stacks, nk_stacks
     INTEGER(KIND=int_1), DIMENSION(:,:,:), POINTER :: stack_map
     INTEGER                                        :: default_stack
     REAL(KIND=sp), DIMENSION(:), POINTER :: row_max_epss, a_norms, b_norms
     REAL(KIND=real_8)     :: eps
     TYPE(dbcsr_work_type) :: my_wm
     INTEGER               :: lastblk, datasize
     INTEGER               :: original_lastblk
     INTEGER(kind=int_8)   :: flop
     REAL(KIND=dp)         :: t_index, t_gemm
     INTEGER, DIMENSION(:,:), POINTER :: right_data_sr, right_data_rr,&
                                         left_data_sr, left_data_rr,&
                                         right_index_sr, right_index_rr,&
                                         left_index_sr, left_index_rr
  END TYPE carrier_type

  
  
END MODULE dbcsr_mm_types
