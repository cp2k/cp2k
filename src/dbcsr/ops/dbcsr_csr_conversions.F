!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief DBCSR to CSR matrix format conversion
!> \par History
!>       2014.10 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************

MODULE dbcsr_csr_conversions

  USE array_types,                     ONLY: array_data,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release
  USE dbcsr_data_methods,              ONLY: dbcsr_data_clear_pointer,&
                                             dbcsr_data_init,&
                                             dbcsr_data_new,&
                                             dbcsr_data_release
  USE dbcsr_data_types,                ONLY: dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8,&
                                             dbcsr_type_real_default
  USE dbcsr_dist_methods,              ONLY: dbcsr_distribution_col_dist,&
                                             dbcsr_distribution_init,&
                                             dbcsr_distribution_mp,&
                                             dbcsr_distribution_new,&
                                             dbcsr_distribution_release
  USE dbcsr_error_handling,            ONLY: &
       dbcsr_assert, dbcsr_caller_error, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_fatal_level, dbcsr_internal_error, &
       dbcsr_warning_level, dbcsr_wrong_args_error
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_sizes, dbcsr_distribution, dbcsr_get_data_type, &
       dbcsr_get_num_blocks, dbcsr_get_nze, dbcsr_has_symmetry, dbcsr_init, &
       dbcsr_is_initialized, dbcsr_name, dbcsr_nblkcols_total, &
       dbcsr_nblkrows_total, dbcsr_nfullrows_local, dbcsr_release, &
       dbcsr_row_block_sizes, dbcsr_valid_index
  USE dbcsr_mp_methods,                ONLY: dbcsr_mp_group,&
                                             dbcsr_mp_mynode,&
                                             dbcsr_mp_new,&
                                             dbcsr_mp_numnodes,&
                                             dbcsr_mp_release
  USE dbcsr_operations,                ONLY: dbcsr_add,&
                                             dbcsr_copy,&
                                             dbcsr_get_info,&
                                             dbcsr_norm
  USE dbcsr_transformations,           ONLY: dbcsr_complete_redistribute,&
                                             dbcsr_desymmetrize_deep
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_distribution_obj,&
                                             dbcsr_iterator,&
                                             dbcsr_mp_obj,&
                                             dbcsr_norm_maxabsnorm,&
                                             dbcsr_obj
  USE dbcsr_work_operations,           ONLY: dbcsr_create
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8,&
                                             real_4,&
                                             real_8,&
                                             sp
  USE message_passing,                 ONLY: mp_environ,&
                                             mp_gather,&
                                             mp_recv,&
                                             mp_send,&
                                             mp_sum

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_csr_conversions'

  LOGICAL, PARAMETER, PRIVATE          :: careful_mod = .FALSE.

  INTEGER, PARAMETER, PUBLIC           :: csr_dbcsr_blkrow_dist=1, csr_eqrow_ceil_dist=2, csr_eqrow_floor_dist=3

! *****************************************************************************
!> \brief Mapping data relating local CSR indices to local indices of a block-row 
!>          distributed (BRD) DBCSR matrix, and containing the block structure 
!>          of the original DBCSR matrix from which the CSR matrix was created.
!> \param csr_to_brd_ind    csr_to_brd_ind(csr_ind) gives the location of a 
!>                          matrix element with CSR index csr_ind (location in 
!>                          nzval_local) inside the data_area of the corresponding 
!>                          BRD matrix. If an element of the DBCSR matrix is 
!>                          treated as 0 in the CSR format, the index of this value
!>                          is not in csr_to_brd_ind.
!> \param brd_to_csr_ind    same as csr_to_brd_ind but inverse mapping. If a 
!>                          given DBCSR index dbcsr_ind points to a zero element, 
!>                          then brd_to_csr_ind(dbcsr_ind) is -1.
!> \param brd_mat DBCSR     BRD matrix acting as an intermediate step in any 
!>                          conversion from and to DBCSR format.
!> \param dbcsr_*           data from original DBCSR matrix (not block-row 
!>                          distributed), representing the original block 
!>                          structure.
!> \param has_dbcsr_block_data    whether dbcsr_* fields are defined
!> \par History
!>      10.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  TYPE csr_mapping_data
    PRIVATE
    INTEGER, DIMENSION(:), POINTER           :: csr_to_brd_ind => NULL(),&
                                                brd_to_csr_ind => NULL()
    TYPE(dbcsr_obj)                          :: brd_mat
    
    LOGICAL                                  :: has_dbcsr_block_data = .FALSE.
    INTEGER                                  :: dbcsr_nblkcols_total,&
                                                dbcsr_nblkrows_total,&
                                                dbcsr_nblks_local
    INTEGER, DIMENSION(:), POINTER           :: dbcsr_row_p, dbcsr_col_i, &
                                                dbcsr_row_blk_size, dbcsr_col_blk_size
  ENDTYPE

! *****************************************************************************
!> \brief Data type of CSR matrices
!> \param r_sp            real, single precision data array
!> \param r_dp            real, double precision data array
!> \param c_sp            complex, single precision data array
!> \param c_sp            complex, double precision data array
!> \param data_type       data type of CSR matrix
!> \par History
!>      11.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  TYPE csr_data_area_type
     REAL(KIND=real_4), DIMENSION(:), POINTER      :: r_sp  => Null()
     REAL(KIND=real_8), DIMENSION(:), POINTER      :: r_dp  => Null()
     COMPLEX(KIND=real_4), DIMENSION(:), POINTER   :: c_sp  => Null()
     COMPLEX(KIND=real_8), DIMENSION(:), POINTER   :: c_dp  => Null()
     INTEGER                                       :: data_type  = -1
  ENDTYPE

! *****************************************************************************
!> \brief Type for CSR matrices
!> \param nrows_total      total number of rows
!> \param ncols_total      total number of columns
!> \param nze_total        total number of nonzero elements
!> \param nze_local        local number of nonzero elements
!> \param nrows_local      local number of rows
!> \param nzval_local      values of local non-zero elements, row-wise ordering.
!> \param rowptr_local     indices of elements inside nzval_local starting a row
!> \param colind_local     column indices of elements inside nzval_local
!> \param dbcsr_mapping    mapping data relating indices of nzval_local to indices 
!>                         of a block-row distributed DBCSR matrix
!> \param has_mapping      whether dbcsr_mapping is defined
!> \param valid            whether essential data (excluding dbcsr_mapping) is 
!>                         completely allocated
!> \param has_indices      whether rowptr_local and colind_local are defined
!> \param mp_group         message-passing group ID
!> \par History
!>      10.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  TYPE csr_type
    INTEGER                                  :: nrows_total, ncols_total,&
                                                nze_local, nrows_local, mp_group
    INTEGER(KIND=int_8)                      :: nze_total
    INTEGER, DIMENSION(:), POINTER           :: rowptr_local => NULL(),&
                                                colind_local => NULL(),&
                                                nzerow_local => NULL()
    TYPE(csr_data_area_type)                 :: nzval_local
    TYPE(csr_mapping_data)                   :: dbcsr_mapping
    LOGICAL                                  :: has_mapping = .FALSE.
    LOGICAL                                  :: valid = .FALSE.
    LOGICAL                                  :: has_indices = .FALSE.
  ENDTYPE csr_type

  PUBLIC :: csr_type,convert_csr_to_dbcsr,&
            csr_create_from_dbcsr,&
            csr_destroy, csr_conversion_test,&
            convert_dbcsr_to_csr, csr_write,&
            csr_create, csr_create_new, csr_create_template,&
            csr_copy

  INTERFACE csr_create
     MODULE PROCEDURE csr_create_new, csr_create_template
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief Create a new CSR matrix and allocate all internal data (excluding dbcsr_mapping)
!> \param[out] csr_mat       CSR matrix to return
!> \param[in] nrows_total    total number of rows
!> \param[in] ncols_total    total number of columns
!> \param[in] nze_total      total number of non-zero elements
!> \param[in] nze_local      local number of non-zero elements
!> \param[in] nrows_local    local number of rows
!> \param mp_group ...
!> \param[in] data_type      data type of the CSR matrix (default real double prec.)
!> \param error ...
!> \par History
!>      10.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE csr_create_new(csr_mat, nrows_total, ncols_total, nze_total,&
                            nze_local, nrows_local, mp_group, data_type, error) 
    TYPE(csr_type), INTENT(OUT)              :: csr_mat
    INTEGER, INTENT(IN)                      :: nrows_total, ncols_total
    INTEGER(KIND=int_8)                      :: nze_total
    INTEGER, INTENT(IN)                      :: nze_local, nrows_local, &
                                                mp_group
    INTEGER, INTENT(IN), OPTIONAL            :: data_type
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csr_create_new', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert(nrows_total, "GE", nrows_local,& 
      dbcsr_fatal_level, dbcsr_wrong_args_error,routineN,&
      "local number of rows must not exceed total number of rows",&
      __LINE__, error=error)

    CALL dbcsr_assert(nze_total .GE. nze_local,&
      dbcsr_fatal_level, dbcsr_wrong_args_error,routineN,&
      "local number of non-zero elements must not exceed total number of non-zero elements",&
      __LINE__, error=error)

    CALL dbcsr_assert(INT(nrows_total,kind=int_8)*INT(ncols_total,kind=int_8) .GE. nze_total,&
      dbcsr_fatal_level, dbcsr_wrong_args_error,routineN, &
      "Total number of non-zero elements must not exceed total matrix size",&
      __LINE__, error=error)

    CALL dbcsr_assert(nrows_local*ncols_total, "GE", nze_local,&
      dbcsr_fatal_level, dbcsr_wrong_args_error,routineN, &
      "Local number of non-zero elements must not exceed local matrix size",&
      __LINE__, error=error)

    csr_mat%ncols_total = ncols_total
    csr_mat%nrows_total = nrows_total
    csr_mat%nze_total = nze_total
    csr_mat%nze_local = nze_local
    ALLOCATE( csr_mat%colind_local(nze_local) )
    csr_mat%nrows_local = nrows_local
    ALLOCATE( csr_mat%rowptr_local(nrows_local + 1) )
    ALLOCATE( csr_mat%nzerow_local(nrows_local) )

    IF (PRESENT(data_type)) THEN
      csr_mat%nzval_local%data_type = data_type
    ELSE
      csr_mat%nzval_local%data_type = dbcsr_type_real_default
    ENDIF

    SELECT CASE (csr_mat%nzval_local%data_type)
      CASE (dbcsr_type_real_4)
        ALLOCATE(csr_mat%nzval_local%r_sp(nze_local))
      CASE (dbcsr_type_real_8)
        ALLOCATE(csr_mat%nzval_local%r_dp(nze_local))
      CASE (dbcsr_type_complex_4)
        ALLOCATE(csr_mat%nzval_local%c_sp(nze_local))
      CASE (dbcsr_type_complex_8)
        ALLOCATE(csr_mat%nzval_local%c_dp(nze_local))
      CASE DEFAULT
        CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
                          routineN, "Invalid matrix type",__LINE__,error)
    END SELECT

    csr_mat%mp_group = mp_group

    csr_mat%valid = .TRUE.
    csr_mat%has_mapping = .FALSE.
    csr_mat%has_indices = .FALSE.

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE csr_create_new

! *****************************************************************************
!> \brief Create a new CSR matrix and allocate all internal data using 
!>        an existing CSR matrix. Copies the indices but no actual matrix data.
!> \param[out] matrix_b Target CSR matrix
!> \param[int] matrix_a Source CSR matrix
!> \param error ...
!> \par History
!>      11.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE csr_create_template(matrix_b, matrix_a, error)
    TYPE(csr_type), INTENT(OUT)              :: matrix_b
    TYPE(csr_type), INTENT(IN)               :: matrix_a
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csr_create_template', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    TYPE(csr_mapping_data)                   :: map

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert(matrix_a%valid,dbcsr_fatal_level,dbcsr_wrong_args_error,routineN,&
      "Source CSR matrix must be created.",__LINE__, error)

    CALL csr_create_new(matrix_b, matrix_a%nrows_total, matrix_a%ncols_total,&
      matrix_a%nze_total, matrix_a%nze_local, matrix_a%nrows_local,&
      matrix_a%mp_group, matrix_a%nzval_local%data_type, error)

    matrix_b%mp_group = matrix_a%mp_group
    matrix_b%has_mapping = matrix_a%has_mapping
    matrix_b%has_indices = matrix_a%has_indices

    IF(matrix_a%has_indices) THEN
      matrix_b%rowptr_local = matrix_a%rowptr_local
      matrix_b%nzerow_local = matrix_a%nzerow_local
      matrix_b%colind_local = matrix_a%colind_local
    ENDIF
    
    IF (matrix_a%has_mapping) THEN
      map = matrix_a%dbcsr_mapping
      ALLOCATE(matrix_b%dbcsr_mapping%csr_to_brd_ind(SIZE(map%csr_to_brd_ind)))
      ALLOCATE(matrix_b%dbcsr_mapping%brd_to_csr_ind(SIZE(map%brd_to_csr_ind)))     
      matrix_b%dbcsr_mapping%csr_to_brd_ind = map%csr_to_brd_ind
      matrix_b%dbcsr_mapping%brd_to_csr_ind = map%brd_to_csr_ind
      matrix_b%dbcsr_mapping%has_dbcsr_block_data = map%has_dbcsr_block_data
      IF (matrix_b%dbcsr_mapping%has_dbcsr_block_data) THEN
        matrix_b%dbcsr_mapping%dbcsr_nblkcols_total = map%dbcsr_nblkcols_total
        matrix_b%dbcsr_mapping%dbcsr_nblkrows_total = map%dbcsr_nblkrows_total
        ALLOCATE(matrix_b%dbcsr_mapping%dbcsr_row_blk_size(map%dbcsr_nblkrows_total))
        ALLOCATE(matrix_b%dbcsr_mapping%dbcsr_col_blk_size(map%dbcsr_nblkcols_total))
        ALLOCATE(matrix_b%dbcsr_mapping%dbcsr_row_p(map%dbcsr_nblkrows_total+1))
        ALLOCATE(matrix_b%dbcsr_mapping%dbcsr_col_i(map%dbcsr_nblks_local))
        matrix_b%dbcsr_mapping%dbcsr_nblks_local = map%dbcsr_nblks_local
        matrix_b%dbcsr_mapping%dbcsr_row_p = map%dbcsr_row_p
        matrix_b%dbcsr_mapping%dbcsr_col_i = map%dbcsr_col_i
        matrix_b%dbcsr_mapping%dbcsr_row_blk_size = map%dbcsr_row_blk_size
        matrix_b%dbcsr_mapping%dbcsr_col_blk_size = map%dbcsr_col_blk_size
      ENDIF

      CALL dbcsr_init(matrix_b%dbcsr_mapping%brd_mat)
      CALL dbcsr_copy(matrix_b%dbcsr_mapping%brd_mat, map%brd_mat, error=error)

    matrix_b%valid = .TRUE.
      
    ENDIF

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE csr_create_template

! *****************************************************************************
!> \brief Copy a CSR matrix.
!> \param[out] matrix_b Target CSR matrix
!> \param[in] matrix_a Source CSR matrix
!> \param error ...
!> \par History
!>      10.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE csr_copy(matrix_b, matrix_a, error)
    TYPE(csr_type), INTENT(OUT)              :: matrix_b
    TYPE(csr_type), INTENT(IN)               :: matrix_a
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csr_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert(matrix_a%valid,dbcsr_fatal_level,dbcsr_wrong_args_error,routineN,&
      "Source CSR matrix must be created.",__LINE__, error)
    CALL dbcsr_assert(matrix_a%has_indices,dbcsr_fatal_level,dbcsr_wrong_args_error,routineN,&
      "Source CSR matrix must have indices.",__LINE__, error)

    IF (.NOT. matrix_b%valid) THEN
      CALL csr_create(matrix_b, matrix_a, error)
    ENDIF

    CALL dbcsr_assert(matrix_b%nzval_local%data_type, "EQ",&
      matrix_a%nzval_local%data_type,dbcsr_fatal_level,dbcsr_wrong_args_error,&
      routineN,"Matrices don't have the same datatype.",__LINE__, error)

    CALL dbcsr_assert(matrix_b%nze_local, "EQ",&
      matrix_a%nze_local,dbcsr_fatal_level,dbcsr_wrong_args_error,routineN,&
      "Matrices don't have the same number of local non-zero elements.",__LINE__, error)

    CALL dbcsr_assert(matrix_b%nrows_local, "EQ", matrix_a%nrows_local,&
      dbcsr_fatal_level,dbcsr_wrong_args_error,routineN,&
      "Matrices don't have the same number of local rows.",__LINE__, error)

    IF (.NOT. matrix_b%has_indices) THEN
      matrix_b%rowptr_local = matrix_a%rowptr_local
      matrix_b%nzerow_local = matrix_a%nzerow_local
      matrix_b%colind_local = matrix_a%colind_local
    ENDIF

    SELECT CASE (matrix_a%nzval_local%data_type)
       CASE (dbcsr_type_real_4)
         matrix_b%nzval_local%r_sp = matrix_a%nzval_local%r_sp
       CASE (dbcsr_type_real_8)
         matrix_b%nzval_local%r_dp = matrix_a%nzval_local%r_dp
       CASE (dbcsr_type_complex_4)
         matrix_b%nzval_local%c_sp = matrix_a%nzval_local%c_sp
       CASE (dbcsr_type_complex_8)
         matrix_b%nzval_local%c_dp = matrix_a%nzval_local%c_dp
     END SELECT

  CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE csr_copy

! *****************************************************************************
!> \brief create a vector containing the number of non-zero elements in each
!> row of a CSR matrix
!> \param[in] csr_mat   CSR matrix
!> \param[out] nzerow  number of non-zero elements in each row
!> \param error ...
!> \par History
!>      11.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE csr_create_nzerow(csr_mat, nzerow, error)
    TYPE(csr_type), INTENT(IN)               :: csr_mat
    INTEGER, DIMENSION(:), INTENT(INOUT), &
      POINTER                                :: nzerow
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csr_create_nzerow', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, k

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert(csr_mat%valid,dbcsr_fatal_level,dbcsr_wrong_args_error,routineN,&
      "CSR matrix must be created.",__LINE__, error)

    DO k=1,csr_mat%nrows_local
      nzerow(k) = csr_mat%rowptr_local(k+1) - csr_mat%rowptr_local(k) 
    ENDDO

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE csr_create_nzerow
   
! *****************************************************************************
!> \brief destroy a CSR matrix
!> \param[in,out] csr_mat ...
!> \param error ...
!> \par History
!>      10.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE csr_destroy(csr_mat, error)
    TYPE(csr_type), INTENT(INOUT)            :: csr_mat
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csr_destroy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    TYPE(csr_mapping_data)                   :: map

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert(csr_mat%valid,dbcsr_fatal_level,dbcsr_wrong_args_error,routineN,& 
      "CSR matrix must be created before destroying it.",__LINE__, error)

    IF (ASSOCIATED(csr_mat%rowptr_local)) DEALLOCATE(csr_mat%rowptr_local)
    IF (ASSOCIATED(csr_mat%nzerow_local)) DEALLOCATE(csr_mat%nzerow_local)
    IF (ASSOCIATED(csr_mat%colind_local)) DEALLOCATE(csr_mat%colind_local)

    IF (csr_mat%has_mapping) THEN
      map = csr_mat%dbcsr_mapping
      IF (ASSOCIATED(map%csr_to_brd_ind))&
        DEALLOCATE(map%csr_to_brd_ind)
      IF (ASSOCIATED(map%brd_to_csr_ind))&
        DEALLOCATE(map%brd_to_csr_ind)
      IF (ASSOCIATED(map%dbcsr_row_blk_size))& 
        DEALLOCATE(map%dbcsr_row_blk_size)
      IF (ASSOCIATED(map%dbcsr_col_blk_size))& 
        DEALLOCATE(map%dbcsr_col_blk_size)
      IF (ASSOCIATED(map%dbcsr_row_p))&
        DEALLOCATE(map%dbcsr_row_p)
      IF (ASSOCIATED(map%dbcsr_col_i))&
        DEALLOCATE(map%dbcsr_col_i)
  
      IF (dbcsr_is_initialized(map%brd_mat))&
        CALL dbcsr_release(map%brd_mat)
    ENDIF

    IF (ASSOCIATED(csr_mat%nzval_local%r_dp)) &
      DEALLOCATE(csr_mat%nzval_local%r_dp)
    IF (ASSOCIATED(csr_mat%nzval_local%r_sp)) &
      DEALLOCATE(csr_mat%nzval_local%r_sp)
    IF (ASSOCIATED(csr_mat%nzval_local%c_sp)) &
      DEALLOCATE(csr_mat%nzval_local%c_sp)
    IF (ASSOCIATED(csr_mat%nzval_local%c_dp)) &
      DEALLOCATE(csr_mat%nzval_local%c_dp)

    csr_mat%has_mapping = .FALSE.
    csr_mat%valid = .FALSE.
    csr_mat%dbcsr_mapping%has_dbcsr_block_data = .FALSE.
    csr_mat%has_indices = .FALSE.
    csr_mat%nzval_local%data_type = -1

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE csr_destroy
 
! *****************************************************************************
!> \brief Allocate the internals of a CSR matrix using data from a block-row 
!>        distributed DBCSR matrix
!> \param[in] brd_mat block-row-distributed DBCSR matrix
!> \param[out] csr_mat CSR matrix
!> \param[in] dbcsr_zero_ind  indices of elements in the data_area of the DBCSR matrix 
!>                            that shall be treated as zeroes in the CSR matrix (optional).
!> \param error ...
!> \par History
!>      10.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE csr_create_from_brd( brd_mat, csr_mat, dbcsr_zero_ind, error) 
    TYPE(dbcsr_obj), INTENT(IN)              :: brd_mat
    TYPE(csr_type), INTENT(OUT)              :: csr_mat
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: dbcsr_zero_ind
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csr_create_from_brd', &
      routineP = moduleN//':'//routineN

    INTEGER :: data_type, error_handler, mp_group, nfullcols_total, &
      nfullrows, nfullrows_total, nze_local
    INTEGER(KIND=int_8)                      :: nze_total
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: zero_ind
    INTEGER, DIMENSION(:), POINTER           :: cdist, csr_index, dbcsr_index
    TYPE(array_i1d_obj)                      :: col_dist
    TYPE(dbcsr_distribution_obj)             :: dist_current

    CALL dbcsr_error_set(routineN, error_handler, error)
    NULLIFY(dbcsr_index, csr_index, cdist)

    IF (PRESENT(dbcsr_zero_ind)) THEN
      ALLOCATE(zero_ind(SIZE(dbcsr_zero_ind)))
      zero_ind(:) = dbcsr_zero_ind(:)
    ELSE
      ALLOCATE(zero_ind(0))
    ENDIF
      
    dist_current = dbcsr_distribution(brd_mat)

    mp_group = dbcsr_mp_group(dbcsr_distribution_mp(dist_current))
    col_dist = dbcsr_distribution_col_dist(dist_current)
    cdist => array_data(col_dist)

    CALL dbcsr_assert(ALL(cdist .EQ. 0), dbcsr_fatal_level,dbcsr_wrong_args_error,&
      routineN, "DBCSR matrix not block-row distributed.",__LINE__, error)
 
    ! Calculate mapping between BRD and CSR indices
    CALL csr_get_dbcsr_mapping (brd_mat, dbcsr_index, csr_index, nze_local, &
                                zero_ind, error=error)

    CALL dbcsr_get_info(brd_mat, nfullrows_total = nfullrows_total,&
                        nfullcols_total = nfullcols_total)

    ! Sum up local number of non-zero elements to get total number
    nze_total = nze_local
    CALL mp_sum(nze_total,mp_group) 

    nfullrows = dbcsr_nfullrows_local(brd_mat)
    data_type = dbcsr_get_data_type(brd_mat)

    ! Allocate CSR matrix
    CALL csr_create_new(csr_mat, nfullrows_total, nfullcols_total, nze_total,&
                        nze_local, nfullrows, mp_group, data_type, error)

    csr_mat%dbcsr_mapping%brd_to_csr_ind => csr_index
    csr_mat%dbcsr_mapping%csr_to_brd_ind => dbcsr_index

    csr_mat%has_mapping = .TRUE.
    csr_mat%dbcsr_mapping%brd_mat = brd_mat

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE csr_create_from_brd

! *****************************************************************************
!> \brief create the mapping information between a block-row distributed DBCSR 
!>        matrix and the corresponding CSR matrix
!> \param[in] brd_mat          the block-row distributed DBCSR matrix
!> \param[out] dbcsr_index     csr to dbcsr index mapping
!> \param[out] csr_index       dbcsr to csr index mapping
!> \param[out] csr_nze_local   number of local non-zero elements
!> \param[in] dbcsr_zero_ind   see documentation of csr_create_from_brd
!> \param error ...
!> \par History
!>      10.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE csr_get_dbcsr_mapping (brd_mat, dbcsr_index, csr_index, csr_nze_local, &
                                    dbcsr_zero_ind, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: brd_mat
    INTEGER, DIMENSION(:), INTENT(OUT), &
      POINTER                                :: dbcsr_index, csr_index
    INTEGER, INTENT(OUT)                     :: csr_nze_local
    INTEGER, DIMENSION(:), INTENT(IN)        :: dbcsr_zero_ind
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csr_get_dbcsr_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, blkcol, blkrow, col_blk_size, csr_ind, dbcsr_ind, el_sum, &
      error_handler, fullcol_sum_blkrow, l, m, n, nblkrows_total, nze, &
      prev_blk, prev_blkrow, prev_row_blk_size, row_blk_offset, row_blk_size
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nfullcol_blkrow
    INTEGER, DIMENSION(:), POINTER           :: dbcsr_index_nozeroes
    LOGICAL                                  :: tr
    TYPE(dbcsr_iterator)                     :: iter

    CALL dbcsr_error_set(routineN, error_handler, error)

    NULLIFY(dbcsr_index, csr_index) 

    CALL dbcsr_get_info(brd_mat, nblkrows_total = nblkrows_total)
    nze = dbcsr_get_nze(brd_mat)

    ALLOCATE(nfullcol_blkrow(nblkrows_total))
    ALLOCATE(dbcsr_index(nze)) 
    ALLOCATE(csr_index(nze))
 
    CALL dbcsr_iterator_start(iter, brd_mat, read_only = .TRUE.)
    nfullcol_blkrow = 0 ! number of non-zero full columns in each block row
    prev_blk = 0
 
    DO WHILE (dbcsr_iterator_blocks_left(iter))
      CALL dbcsr_iterator_next_block(iter, blkrow, blkcol, blk, transposed = tr,&
        col_size = col_blk_size)

      CALL dbcsr_assert(blk, "EQ", prev_blk + 1,&
        dbcsr_fatal_level, dbcsr_internal_error,routineN,&
        "iterator is required to traverse the blocks in a row-wise fashion",&
        __LINE__, error=error)

      prev_blk = blk

      nfullcol_blkrow(blkrow) = nfullcol_blkrow(blkrow) + col_blk_size
      CALL dbcsr_assert(.NOT. tr,dbcsr_fatal_level, dbcsr_wrong_args_error,routineN,&
        "DBCSR block data must not be transposed",__LINE__, error=error)
    ENDDO
    CALL dbcsr_iterator_stop(iter)

    el_sum = 0 ! number of elements above current block row

    prev_blkrow = 0 ! store number and size of previous block row
    prev_row_blk_size = 0

    CALL dbcsr_iterator_start(iter, brd_mat, read_only = .TRUE.)
    DO WHILE (dbcsr_iterator_blocks_left(iter)) 

      CALL dbcsr_iterator_next_block(iter, blkrow, blkcol, blk, transposed = tr, &
        row_size = row_blk_size, col_size = col_blk_size, row_offset = row_blk_offset)

      IF(blkrow .GT. prev_blkrow) THEN ! new block row
        IF (prev_blkrow .GT. 0) THEN
          el_sum = el_sum + nfullcol_blkrow(prev_blkrow)*prev_row_blk_size
        ENDIF
         
        ! number of non-zero full columns on the left of current block:
        fullcol_sum_blkrow = 0         
        
        dbcsr_ind = el_sum
      ENDIF
      DO n = 1, col_blk_size !nr of columns
        DO m = 1, row_blk_size !nr of rows
          dbcsr_ind = dbcsr_ind + 1
          csr_ind = (m-1)*nfullcol_blkrow(blkrow) + fullcol_sum_blkrow + n + el_sum
          dbcsr_index(csr_ind) = dbcsr_ind
          csr_index(dbcsr_ind) = csr_ind
        ENDDO
      ENDDO
      fullcol_sum_blkrow = fullcol_sum_blkrow + col_blk_size
      prev_blkrow = blkrow 
      prev_row_blk_size = row_blk_size
    ENDDO
    CALL dbcsr_iterator_stop(iter)
 
    IF (SIZE(dbcsr_zero_ind) .GT. 0) THEN
      ALLOCATE(dbcsr_index_nozeroes(nze-SIZE(dbcsr_zero_ind)))
      m=0 ! csr index if zeroes are excluded from CSR data
      DO l=1,nze ! csr index if zeroes are included in CSR data
        IF (ANY(dbcsr_zero_ind .EQ. dbcsr_index(l))) THEN
          csr_index(dbcsr_index(l)) = -1
        ELSE
          m = m+1
          dbcsr_index_nozeroes(m) = dbcsr_index(l)
          csr_index(dbcsr_index(l)) = m
        ENDIF
      ENDDO
      DEALLOCATE(dbcsr_index)
      dbcsr_index => dbcsr_index_nozeroes
    ENDIF

    IF (SIZE(dbcsr_zero_ind) .GT. 0) THEN
      csr_nze_local = m
    ELSE
      csr_nze_local = nze
    ENDIF

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE csr_get_dbcsr_mapping


! *****************************************************************************
!> \brief Copies data from a CSR matrix to a block-row distributed DBCSR matrix. 
!>        The DBCSR matrix must have a block structure consistent with the CSR matrix.
!> \param[in,out] brd_mat    block-row distributed DBCSR matrix
!> \param[in] csr_mat        CSR matrix
!> \param error ...
!> \par History
!>      10.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE convert_csr_to_brd (brd_mat, csr_mat, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: brd_mat
    TYPE(csr_type), INTENT(IN)               :: csr_mat
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'convert_csr_to_brd', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_type, error_handler, &
                                                ind, k, nze

    CALL dbcsr_error_set(routineN, error_handler, error)

    data_type = dbcsr_get_data_type(brd_mat) 
    nze = dbcsr_get_nze(brd_mat)
    CALL dbcsr_data_release(brd_mat%m%data_area)
    CALL dbcsr_data_new (brd_mat%m%data_area, data_type, nze)

    SELECT CASE (data_type)
      CASE (dbcsr_type_real_4)
        brd_mat%m%data_area%d%r_sp = 0.0_sp
      CASE (dbcsr_type_real_8)
        brd_mat%m%data_area%d%r_dp = 0.0_dp
      CASE (dbcsr_type_complex_4)
        brd_mat%m%data_area%d%c_sp = 0.0_sp
      CASE (dbcsr_type_complex_8)
        brd_mat%m%data_area%d%c_dp = 0.0_dp
      END SELECT

    DO k=1,csr_mat%nze_local
      ind = csr_mat%dbcsr_mapping%csr_to_brd_ind(k)
      SELECT CASE (data_type)
        CASE (dbcsr_type_real_4)
          brd_mat%m%data_area%d%r_sp(ind) = csr_mat%nzval_local%r_sp(k)
        CASE (dbcsr_type_real_8)
          brd_mat%m%data_area%d%r_dp(ind) = csr_mat%nzval_local%r_dp(k)
        CASE (dbcsr_type_complex_4)
          brd_mat%m%data_area%d%c_sp(ind) = csr_mat%nzval_local%c_sp(k)
        CASE (dbcsr_type_complex_8)
          brd_mat%m%data_area%d%c_dp(ind) = csr_mat%nzval_local%c_dp(k)
      END SELECT
    ENDDO

   CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE convert_csr_to_brd

! *****************************************************************************
!> \brief Convert a block-row distributed DBCSR matrix to a CSR matrix
!>        The DBCSR matrix must have a block structure consistent with the CSR matrix.
!> \param[in] brd_mat        block-row distributed DBCSR matrix 
!> \param[in,out] csr_mat    CSR matrix 
!> \param error ...
!> \par History
!>      10.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE convert_brd_to_csr (brd_mat, csr_mat, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: brd_mat
    TYPE(csr_type), INTENT(INOUT)            :: csr_mat
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'convert_brd_to_csr', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, blkcol, blkrow, col_blk_offset, col_blk_size, csr_ind, &
      data_type, dbcsr_ind, el_sum, error_handler, ind_blk_data, k, &
      local_row_ind, m, n, nblkrows_total, node_row_offset, prev_blkrow, &
      prev_row_blk_size, row_blk_offset, row_blk_size
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nfullcol_blkrow
    INTEGER, DIMENSION(:), POINTER           :: colind, csr_index, &
                                                dbcsr_index, nzerow, rowptr
    LOGICAL                                  :: new_ind, tr
    TYPE(dbcsr_data_obj)                     :: block
    TYPE(dbcsr_iterator)                     :: iter

    CALL dbcsr_error_set(routineN, error_handler, error)
    NULLIFY(rowptr, colind, dbcsr_index, csr_index)

    dbcsr_index => csr_mat%dbcsr_mapping%csr_to_brd_ind
    csr_index => csr_mat%dbcsr_mapping%brd_to_csr_ind 

    
    ! CSR indices are not recalculated if indices are already defined
    new_ind = .NOT. (csr_mat%has_indices)

    CALL  dbcsr_assert(csr_mat%has_mapping,dbcsr_fatal_level, dbcsr_wrong_args_error,routineN,&
        "DBCSR mapping of CSR matrix must be defined",__LINE__, error=error)

    ! Calculate mapping between CSR matrix and DBCSR matrix if not yet defined
    !IF (.NOT. csr_mat%has_mapping ) THEN 
    !  CALL csr_get_dbcsr_mapping (brd_mat, dbcsr_index, csr_index, nze, error=error)
    !ENDIF
   
    CALL dbcsr_get_info(brd_mat, nblkrows_total = nblkrows_total)
    ALLOCATE( nfullcol_blkrow(nblkrows_total))
   
    ! iteration over blocks without touching data, 
    ! in order to get number of non-zero full columns in each block row
    CALL dbcsr_iterator_start(iter, brd_mat, read_only = .TRUE.)
    blkrow = 0
    nfullcol_blkrow = 0 ! number of non-zero full columns in each block row
    data_type = dbcsr_get_data_type(brd_mat)

    DO WHILE (dbcsr_iterator_blocks_left(iter))
      CALL dbcsr_iterator_next_block(iter, blkrow, blkcol, blk, col_size = col_blk_size,&
        row_offset= row_blk_offset)
      nfullcol_blkrow(blkrow) = nfullcol_blkrow(blkrow) + col_blk_size
      IF(blk .EQ. 1) THEN
        node_row_offset=row_blk_offset
      ENDIF 
    ENDDO

    CALL dbcsr_iterator_stop(iter)

    ! Copy data from BRD matrix to CSR matrix and calculate CSR indices
    prev_blkrow = 0
    prev_row_blk_size = 0
    el_sum = 0 ! number of elements above current block row
    colind => csr_mat%colind_local
    rowptr => csr_mat%rowptr_local
    nzerow => csr_mat%nzerow_local
    CALL dbcsr_data_init(block)
    CALL dbcsr_data_new(block, data_type)

    CALL dbcsr_iterator_start(iter, brd_mat, read_only = .TRUE.)

    IF (new_ind) rowptr(:) = 0 ! initialize to 0 in order to check which rows are 0 at a later time
    DO WHILE (dbcsr_iterator_blocks_left(iter))  
      CALL dbcsr_iterator_next_block(iter, blkrow, blkcol, block, tr,&
        col_size = col_blk_size, row_size = row_blk_size, row_offset = row_blk_offset,&
        col_offset = col_blk_offset)

      CALL dbcsr_assert(.NOT. tr,dbcsr_fatal_level, dbcsr_wrong_args_error,routineN,&
        "DBCSR block data must not be transposed",__LINE__, error=error)

      IF(blkrow > prev_blkrow) THEN ! new block row
        local_row_ind = row_blk_offset-node_row_offset ! actually: local row index - 1
        IF (prev_blkrow .GT. 0) THEN
          el_sum = el_sum + nfullcol_blkrow(prev_blkrow)*prev_row_blk_size 
        ENDIF
        dbcsr_ind = el_sum
      ENDIF
      DO n = 1, col_blk_size !nr of columns
        DO m = 1, row_blk_size !nr of rows
          dbcsr_ind = dbcsr_ind + 1
          csr_ind = csr_index(dbcsr_ind) ! get CSR index for current element
          IF (csr_ind .GT. 0) THEN ! is non-zero element if csr_ind > 0
            IF(new_ind) THEN
              ! Calculate CSR column index
              colind(csr_ind) = col_blk_offset + n - 1
              ! Calculate CSR row pointer
              ! (not accounting for zero elements inside non-zero blocks)
              IF (rowptr(local_row_ind + m) .LE. 0) rowptr(local_row_ind + m) =&
                rowptr(local_row_ind + m) + el_sum + 1 + nfullcol_blkrow(blkrow)*(m-1)
            ENDIF
              ind_blk_data = (m+row_blk_size*(n-1)) ! index of data inside DBCSR blocks
            SELECT CASE (csr_mat%nzval_local%data_type)
              CASE (dbcsr_type_real_4)
                csr_mat%nzval_local%r_sp(csr_ind) = block%d%r_sp(ind_blk_data)
              CASE (dbcsr_type_real_8)
                csr_mat%nzval_local%r_dp(csr_ind) = block%d%r_dp(ind_blk_data)
              CASE (dbcsr_type_complex_4)
                csr_mat%nzval_local%c_sp(csr_ind) = block%d%c_sp(ind_blk_data)
              CASE (dbcsr_type_complex_8)
                csr_mat%nzval_local%c_dp(csr_ind) = block%d%c_dp(ind_blk_data)
            END SELECT
          ELSE ! is zero element if ind = -1
            ! CSR row pointer has to be corrected if element is zero 
            ! (subtract 1 from all subsequent row pointers)
            IF(new_ind) rowptr(local_row_ind + m+1 :) = rowptr(local_row_ind + m+1 :) - 1
          ENDIF
        ENDDO
      ENDDO
      prev_blkrow = blkrow
      prev_row_blk_size = row_blk_size
    ENDDO
    
    IF (new_ind) THEN
      ! Repeat previous row pointer for row pointers to zero rows
      IF (csr_mat%nrows_local .GT. 0) rowptr(1) = 1
      DO k=1, csr_mat%nrows_local
        IF (rowptr(k) .LE. 0) rowptr(k) = rowptr(k-1)
      ENDDO

      rowptr(csr_mat%nrows_local+1) = csr_mat%nze_local + 1
    ENDIF

    CALL csr_create_nzerow(csr_mat, nzerow, error)

    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_data_clear_pointer (block)
    CALL dbcsr_data_release (block)

    IF (new_ind) csr_mat%has_indices = .TRUE.

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE convert_brd_to_csr

! *****************************************************************************
!> \brief create CSR matrix including dbcsr_mapping from arbitrary DBCSR matrix
!>        in order to prepare conversion.
!> \param[in] dbcsr_mat ...
!> \param[out] csr_mat ...
!> \param[in] dist_format     how to distribute CSR rows over processes:
!>                            csr_dbcsr_blkrow_dist: the number of rows per process 
!>                            is adapted to the row block sizes in the DBCSR format 
!>                            such that blocks are not splitted over different processes.
!>                            csr_eqrow_ceil_dist: each process holds ceiling(N/P) CSR rows.
!>                            csr_eqrow_floor_dist: each process holds floor(N/P) CSR rows.
!> \param[in] numnodes        number of nodes to use for distributing CSR matrix
!>                            (optional, default is number of nodes used for 
!>                            DBCSR matrix)
!> \param error ...
!> \par History
!>      10.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE csr_create_from_dbcsr(dbcsr_mat, csr_mat, dist_format, numnodes, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: dbcsr_mat
    TYPE(csr_type), INTENT(OUT)              :: csr_mat
    INTEGER, INTENT(IN)                      :: dist_format
    INTEGER, INTENT(IN), OPTIONAL            :: numnodes
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csr_create_from_dbcsr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dbcsr_numnodes, &
                                                error_handler, &
                                                nblkcols_total, &
                                                nblkrows_total, nblks_local, &
                                                num_p
    LOGICAL                                  :: equal_dist, floor_dist
    TYPE(dbcsr_obj)                          :: brd_mat, dbcsr_mat_nosym

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert(dbcsr_valid_index(dbcsr_mat), dbcsr_fatal_level,&
    dbcsr_wrong_args_error, routineN, "Invalid DBCSR matrix", __LINE__, error=error)

    SELECT CASE (dist_format)
    CASE (csr_dbcsr_blkrow_dist)
      equal_dist = .FALSE.
      floor_dist = .FALSE.
    CASE (csr_eqrow_ceil_dist)
      equal_dist = .TRUE.
      floor_dist = .FALSE.
    CASE (csr_eqrow_floor_dist)
      equal_dist = .TRUE.
      floor_dist = .TRUE.
    END SELECT

    ! Conversion does not support matrices in symmetric format, therefore desymmetrize
    IF(dbcsr_has_symmetry(dbcsr_mat)) THEN
      CALL dbcsr_init(dbcsr_mat_nosym)
      CALL dbcsr_desymmetrize_deep(dbcsr_mat,dbcsr_mat_nosym,untransposed_data = .TRUE., error=error)
    ELSE
      dbcsr_mat_nosym = dbcsr_mat
    ENDIF

    dbcsr_numnodes = dbcsr_mp_numnodes(dbcsr_distribution_mp(dbcsr_distribution(dbcsr_mat)))
    IF (PRESENT(numnodes)) THEN
      CALL dbcsr_assert(numnodes .LE. dbcsr_numnodes, dbcsr_fatal_level, &
        dbcsr_wrong_args_error, routineN, &
        "Number of nodes used for CSR matrix must not exceed total number of nodes", &
        __LINE__, error=error)

      num_p = numnodes
    ELSE
      num_p = dbcsr_numnodes
    ENDIF

    CALL dbcsr_create_brd (dbcsr_mat_nosym, brd_mat, equal_dist, floor_dist,&
                           num_p, error=error)

    ! Create CSR matrix from BRD matrix
    CALL csr_create_from_brd(brd_mat, csr_mat, error=error)

    ! Store DBCSR block data inside CSR matrix 
    ! (otherwise, this data is lost when converting from DBCSR to CSR)
    nblks_local = dbcsr_get_num_blocks(dbcsr_mat_nosym)
    nblkrows_total = dbcsr_nblkrows_total(dbcsr_mat_nosym)
    nblkcols_total = dbcsr_nblkcols_total(dbcsr_mat_nosym)

    csr_mat%dbcsr_mapping%dbcsr_nblkcols_total = nblkcols_total
    csr_mat%dbcsr_mapping%dbcsr_nblkrows_total = nblkrows_total
    csr_mat%dbcsr_mapping%dbcsr_nblks_local = nblks_local
    ALLOCATE(csr_mat%dbcsr_mapping%dbcsr_row_p(nblkrows_total+1))
    csr_mat%dbcsr_mapping%dbcsr_row_p = dbcsr_mat_nosym%m%row_p
    ALLOCATE(csr_mat%dbcsr_mapping%dbcsr_col_i(nblks_local))
    csr_mat%dbcsr_mapping%dbcsr_col_i = dbcsr_mat_nosym%m%col_i

    ALLOCATE(csr_mat%dbcsr_mapping%dbcsr_row_blk_size(nblkrows_total))
    ALLOCATE(csr_mat%dbcsr_mapping%dbcsr_col_blk_size(nblkcols_total))

    csr_mat%dbcsr_mapping%dbcsr_row_blk_size = array_data(dbcsr_row_block_sizes(dbcsr_mat_nosym))
    csr_mat%dbcsr_mapping%dbcsr_col_blk_size = array_data(dbcsr_col_block_sizes(dbcsr_mat_nosym))

    csr_mat%dbcsr_mapping%has_dbcsr_block_data = .TRUE.

    IF(dbcsr_has_symmetry(dbcsr_mat)) CALL dbcsr_release(dbcsr_mat_nosym)
 
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE csr_create_from_dbcsr

! *****************************************************************************
!> \brief Check if two DBCSR matrices have the same block structure (including
!>        sparsity pattern)
!> \param[in] matrix_a  ...
!> \param[in] matrix_b ...
!> \retval is_equal whether matrix_a and matrix_b have the same block structure
!> \par History
!>      11.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  FUNCTION dbcsr_has_same_block_structure(matrix_a, matrix_b) RESULT (is_equal)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    LOGICAL                                  :: is_equal

    is_equal = .TRUE.
  
    IF (dbcsr_nblkcols_total(matrix_a) .NE. dbcsr_nblkcols_total(matrix_b)) is_equal = .FALSE.
    IF (dbcsr_nblkrows_total(matrix_a) .NE. dbcsr_nblkrows_total(matrix_b)) is_equal = .FALSE.
    IF ((matrix_a%m%nblks) .NE. (matrix_b%m%nblks)) is_equal = .FALSE.
    IF (ANY(matrix_a%m%row_p .NE. matrix_b%m%row_p)) is_equal = .FALSE.
    IF (ANY(matrix_a%m%col_i .NE. matrix_b%m%col_i)) is_equal = .FALSE.
    IF (ANY(array_data(dbcsr_row_block_sizes(matrix_a)) .NE.&
      array_data(dbcsr_row_block_sizes(matrix_b)))) is_equal = .FALSE.
    IF (ANY(array_data(dbcsr_row_block_sizes(matrix_a)) .NE.&
      array_data(dbcsr_row_block_sizes(matrix_b)))) is_equal = .FALSE.

  END FUNCTION dbcsr_has_same_block_structure

! *****************************************************************************
!> \brief Helper function to assert that a given CSR matrix and a given DBCSR 
!>        matrix are consistent before doing the conversion
!> \param[in] csr_mat ...
!> \param[in] dbcsr_mat ...
!> \param error ...
!> \par History
!>      11.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE csr_assert_consistency_with_dbcsr(csr_mat, dbcsr_mat, error)
    TYPE(csr_type), INTENT(IN)               :: csr_mat
    TYPE(dbcsr_obj), INTENT(IN)              :: dbcsr_mat
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'csr_assert_consistency_with_dbcsr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    TYPE(csr_mapping_data)                   :: map

    CALL dbcsr_error_set(routineN, error_handler, error)
    map = csr_mat%dbcsr_mapping
    IF (map%has_dbcsr_block_data) THEN
      CALL dbcsr_assert(map%dbcsr_nblkcols_total .EQ.&
        dbcsr_nblkcols_total(dbcsr_mat), dbcsr_fatal_level,&
        dbcsr_wrong_args_error,routineN,&
        "field nblkcols_total of DBCSR matrix not consistent with CSR matrix",&
        __LINE__, error=error)
      CALL dbcsr_assert(map%dbcsr_nblkrows_total .EQ.&
        dbcsr_nblkrows_total(dbcsr_mat),dbcsr_fatal_level,&
        dbcsr_wrong_args_error,routineN,&
        "field nblkrows_total of DBCSR matrix not consistent with CSR matrix",&
        __LINE__, error=error)
      CALL dbcsr_assert(map%dbcsr_nblks_local .EQ. dbcsr_mat%m%nblks,&
        dbcsr_fatal_level, dbcsr_wrong_args_error,routineN,&
        "field nblks of DBCSR matrix not consistent with CSR matrix",&
        __LINE__, error=error)
      CALL dbcsr_assert(ALL(map%dbcsr_row_p .EQ. dbcsr_mat%m%row_p),&
        dbcsr_fatal_level, dbcsr_wrong_args_error,routineN,&
        "field row_p of DBCSR matrix not consistent with CSR matrix",&
        __LINE__, error=error)
      CALL dbcsr_assert(ALL(map%dbcsr_col_i .EQ. dbcsr_mat%m%col_i),&
        dbcsr_fatal_level, dbcsr_wrong_args_error,routineN,&
        "field dbcsr_col_i of DBCSR matrix not consistent with CSR matrix",&
        __LINE__, error=error)
      CALL dbcsr_assert((ALL(map%dbcsr_row_blk_size .EQ.&
        array_data(dbcsr_row_block_sizes(dbcsr_mat)))),&
        dbcsr_fatal_level, dbcsr_wrong_args_error,routineN,&
        "field row_blk_size of DBCSR matrix not consistent with CSR matrix",&
        __LINE__, error=error)
      CALL dbcsr_assert((ALL(map%dbcsr_col_blk_size .EQ. &
        array_data(dbcsr_col_block_sizes(dbcsr_mat)))),&
        dbcsr_fatal_level, dbcsr_wrong_args_error,routineN,&
        "field col_blk_size of DBCSR matrix not consistent with CSR matrix",&
        __LINE__, error=error)
    ELSE
      CALL dbcsr_assert(.FALSE., dbcsr_warning_level, dbcsr_caller_error, routineN,&
        "Can not assert consistency of the matrices as no block data stored in CSR matrix.",&
        __LINE__, error=error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE csr_assert_consistency_with_dbcsr

! *****************************************************************************
!> \brief Convert a DBCSR matrix to a CSR matrix.
!> \param[in] dbcsr_mat DBCSR matrix to convert
!> \param[in,out] csr_mat correctly allocated CSR matrix
!> \param error ...
!> \par History
!>      10.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE convert_dbcsr_to_csr(dbcsr_mat, csr_mat, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: dbcsr_mat
    TYPE(csr_type), INTENT(INOUT)            :: csr_mat
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'convert_dbcsr_to_csr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    TYPE(dbcsr_obj)                          :: dbcsr_mat_nosym

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert(dbcsr_valid_index(dbcsr_mat), dbcsr_fatal_level,&
    dbcsr_wrong_args_error, routineN, "Invalid DBCSR matrix", __LINE__, error=error)
    CALL dbcsr_assert(dbcsr_get_data_type(dbcsr_mat),"EQ",&
      csr_mat%nzval_local%data_type, dbcsr_fatal_level, dbcsr_wrong_args_error,&
      routineN, "DBCSR and CSR matrix must have same type", __LINE__, error=error)

    CALL dbcsr_assert(csr_mat%has_mapping,dbcsr_fatal_level,&
        dbcsr_wrong_args_error,routineN,&
        "CSR_mat must contain mapping to DBCSR format",&
        __LINE__, error=error)

    IF(dbcsr_has_symmetry(dbcsr_mat)) THEN
      CALL dbcsr_init(dbcsr_mat_nosym)
      CALL dbcsr_desymmetrize_deep(dbcsr_mat,dbcsr_mat_nosym,untransposed_data = .TRUE., error=error)
    ELSE
      dbcsr_mat_nosym = dbcsr_mat
    ENDIF

    CALL csr_assert_consistency_with_dbcsr(csr_mat, dbcsr_mat_nosym, error)

    ! 1) DBCSR -> BRD
    CALL dbcsr_complete_redistribute(dbcsr_mat_nosym,csr_mat%dbcsr_mapping%brd_mat,error=error)
    ! 2) BRD -> CSR
    CALL convert_brd_to_csr(csr_mat%dbcsr_mapping%brd_mat, csr_mat, error=error)

    IF(dbcsr_has_symmetry(dbcsr_mat)) CALL dbcsr_release(dbcsr_mat_nosym)

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE convert_dbcsr_to_csr

! *****************************************************************************
!> \brief convert a CSR matrix to a DBCSR matrix
!> \param[in,out] dbcsr_mat correctly allocated DBCSR matrix
!> \param[in,out] csr_mat CSR matrix to convert
!> \param error ...
!> \par History
!>      10.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE convert_csr_to_dbcsr(dbcsr_mat, csr_mat, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: dbcsr_mat
    TYPE(csr_type), INTENT(INOUT)            :: csr_mat
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'convert_csr_to_dbcsr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    TYPE(dbcsr_obj)                          :: dbcsr_mat_nosym

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert(dbcsr_valid_index(dbcsr_mat), dbcsr_fatal_level,&
    dbcsr_wrong_args_error, routineN, "Invalid DBCSR matrix", __LINE__, error=error)

    CALL dbcsr_assert(dbcsr_get_data_type(dbcsr_mat),"EQ",&
      csr_mat%nzval_local%data_type, dbcsr_fatal_level,&
      dbcsr_wrong_args_error, routineN, "DBCSR and CSR matrix must have same type",&
      __LINE__, error=error)

    CALL dbcsr_assert(csr_mat%has_mapping,dbcsr_fatal_level,&
        dbcsr_wrong_args_error,routineN,&
        "CSR_mat must contain mapping to DBCSR format",&
        __LINE__, error=error)

    ! Desymmetrize to assert that DBCSR matrix has sparsity pattern consistent with CSR matrix
    IF(dbcsr_has_symmetry(dbcsr_mat)) THEN
      CALL dbcsr_init(dbcsr_mat_nosym)
      CALL dbcsr_desymmetrize_deep(dbcsr_mat,dbcsr_mat_nosym,untransposed_data = .TRUE., error=error)
    ELSE
      dbcsr_mat_nosym = dbcsr_mat
    ENDIF

    CALL csr_assert_consistency_with_dbcsr(csr_mat, dbcsr_mat_nosym, error)
    
    IF(dbcsr_has_symmetry(dbcsr_mat)) CALL dbcsr_release(dbcsr_mat_nosym)

    ! 1) CSR -> BRD
    CALL convert_csr_to_brd(csr_mat%dbcsr_mapping%brd_mat, csr_mat, error)

    ! 2) BRD -> DBCSR
    CALL dbcsr_complete_redistribute(csr_mat%dbcsr_mapping%brd_mat, dbcsr_mat, error=error )

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE convert_csr_to_dbcsr

! *****************************************************************************
!> \brief Test the conversion by converting to CSR format and converting back. 
!>        The maximum norm of the differences between the original and the 
!>        back-converted matrix is calculated.
!> \param[in] dbcsr_mat ...
!> \param[out] csr_mat ...
!> \param[out] norm ...
!> \param error ...
!> \par History
!>      10.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE csr_conversion_test( dbcsr_mat, csr_mat, norm, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: dbcsr_mat
    TYPE(csr_type), INTENT(OUT)              :: csr_mat
    REAL(KIND=real_8), INTENT(OUT)           :: norm
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csr_conversion_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    TYPE(dbcsr_obj)                          :: dbcsr_mat_conv

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL csr_create_from_dbcsr(dbcsr_mat, csr_mat, csr_eqrow_ceil_dist, error=error)

    CALL convert_dbcsr_to_csr(dbcsr_mat, csr_mat, error=error)

    CALL dbcsr_init(dbcsr_mat_conv)
    CALL dbcsr_copy(dbcsr_mat_conv, dbcsr_mat, error=error)

    CALL convert_csr_to_dbcsr(dbcsr_mat_conv, csr_mat, error)

    CALL dbcsr_add(dbcsr_mat_conv, dbcsr_mat, 1.0_dp, -1.0_dp, error)

    CALL dbcsr_norm(dbcsr_mat_conv, dbcsr_norm_maxabsnorm, norm_scalar=norm, error=error)
    CALL dbcsr_release(dbcsr_mat_conv)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE csr_conversion_test
 
! *****************************************************************************
!> \brief Write a CSR matrix to file
!> \param[in] csr_mat ...
!> \param[in] unit_nr           unit number to which output is written
!> \param[in] upper_triangle    If true (default: false), write only upper 
!>                              triangular part of matrix
!> \param error ...
!> \par History
!>      11.2014 created [Patrick Seewald]
!> \author Patrick Seewald
! *****************************************************************************
  SUBROUTINE csr_write(csr_mat, unit_nr, upper_triangle, error)
    TYPE(csr_type), INTENT(IN)               :: csr_mat
    INTEGER, INTENT(IN)                      :: unit_nr
    LOGICAL, INTENT(IN), OPTIONAL            :: upper_triangle
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csr_write', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: data_format
    COMPLEX(KIND=real_4), ALLOCATABLE, &
      DIMENSION(:)                           :: nzval_to_master_c_sp
    COMPLEX(KIND=real_8), ALLOCATABLE, &
      DIMENSION(:)                           :: nzval_to_master_c_dp
    INTEGER                                  :: error_handler, i, k, l, m, &
                                                mynode, numnodes, rowind, &
                                                tag1, tag2, tag3
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: colind_to_master, &
                                                nzerow_to_master, &
                                                sizes_numrowlocal, &
                                                sizes_nzelocal
    LOGICAL                                  :: ut
    REAL(KIND=real_4), ALLOCATABLE, &
      DIMENSION(:)                           :: nzval_to_master_r_sp
    REAL(KIND=real_8), ALLOCATABLE, &
      DIMENSION(:)                           :: nzval_to_master_r_dp

    CALL dbcsr_error_set(routineN, error_handler, error)

    IF (PRESENT(upper_triangle)) THEN
      ut = upper_triangle
    ELSE
      ut = .FALSE.
    ENDIF

    CALL dbcsr_assert(csr_mat%valid,dbcsr_fatal_level,&
        dbcsr_wrong_args_error,routineN,&
        "can not write invalid CSR matrix",&
        __LINE__, error=error)

    tag1 = 0
    tag2 = 1
    tag3 = 2

    CALL mp_environ(numnodes, mynode, csr_mat%mp_group)

    ! gather sizes (number of local non-zero elements and number of local rows)
    IF (mynode .EQ. 0) THEN
      ALLOCATE(sizes_nzelocal(numnodes))
      ALLOCATE(sizes_numrowlocal(numnodes))
    ENDIF

    CALL mp_gather(csr_mat%nze_local,sizes_nzelocal,0,csr_mat%mp_group)
    CALL mp_gather(csr_mat%nrows_local, sizes_numrowlocal,0,csr_mat%mp_group)

    ! for each node, send matrix data to node 0 (master) and write data
    DO i=0, numnodes-1
      IF (mynode .EQ.0) THEN ! allocations for receiving data from node i
        ALLOCATE(colind_to_master(sizes_nzelocal(i+1)))
        ALLOCATE(nzerow_to_master(sizes_numrowlocal(i+1)))

        SELECT CASE (csr_mat%nzval_local%data_type)
          CASE (dbcsr_type_real_4)
            data_format = "(2(I8),E23.6E2)"
            ALLOCATE(nzval_to_master_r_sp(sizes_nzelocal(i+1)))
          CASE (dbcsr_type_real_8)
            data_format = "(2(I8),E23.14E3)"
            ALLOCATE(nzval_to_master_r_dp(sizes_nzelocal(i+1)))
          CASE (dbcsr_type_complex_4)
            data_format = "(2(I8),2(E23.6E2))"
            ALLOCATE(nzval_to_master_c_sp(sizes_nzelocal(i+1)))
          CASE (dbcsr_type_complex_8)
            data_format = "(2(I8),2(E23.14E3))"
            ALLOCATE(nzval_to_master_c_dp(sizes_nzelocal(i+1)))
        END SELECT
      ENDIF

      IF (mynode .EQ. 0) THEN ! receive at node 0
        IF (i .EQ. 0) THEN  ! data from node 0, no need for mpi routines
          colind_to_master(:) = csr_mat%colind_local(:)
          nzerow_to_master(:) = csr_mat%nzerow_local(:)
          SELECT CASE (csr_mat%nzval_local%data_type)
            CASE (dbcsr_type_real_4)
              nzval_to_master_r_sp(:) = csr_mat%nzval_local%r_sp(:)
            CASE (dbcsr_type_real_8)
              nzval_to_master_r_dp(:) = csr_mat%nzval_local%r_dp(:)
            CASE (dbcsr_type_complex_4)
              nzval_to_master_c_sp(:) = csr_mat%nzval_local%c_sp(:)
            CASE (dbcsr_type_complex_8)
              nzval_to_master_c_dp(:) = csr_mat%nzval_local%c_dp(:)
          END SELECT
        ELSE ! receive data from nodes with rank > 0
          CALL mp_recv(colind_to_master,i,tag1,csr_mat%mp_group)
          CALL mp_recv(nzerow_to_master,i,tag2,csr_mat%mp_group)
          SELECT CASE (csr_mat%nzval_local%data_type)
            CASE (dbcsr_type_real_4)
              CALL mp_recv(nzval_to_master_r_sp, i, tag3, csr_mat%mp_group)
            CASE (dbcsr_type_real_8)
              CALL mp_recv(nzval_to_master_r_dp, i, tag3, csr_mat%mp_group)
            CASE (dbcsr_type_complex_4)
              CALL mp_recv(nzval_to_master_c_sp, i, tag3, csr_mat%mp_group)
            CASE (dbcsr_type_complex_8)
              CALL mp_recv(nzval_to_master_c_dp, i, tag3, csr_mat%mp_group)
          END SELECT
        ENDIF
      ENDIF

      IF ((mynode .EQ. i) .AND. (i .NE. 0)) THEN ! send from nodes with rank > 0
        CALL mp_send(csr_mat%colind_local,0,tag1,csr_mat%mp_group)
        CALL mp_send(csr_mat%nzerow_local,0,tag2,csr_mat%mp_group)
        SELECT CASE (csr_mat%nzval_local%data_type)
          CASE (dbcsr_type_real_4)
            CALL mp_send(csr_mat%nzval_local%r_sp, 0, tag3, csr_mat%mp_group)
          CASE (dbcsr_type_real_8)
            CALL mp_send(csr_mat%nzval_local%r_dp, 0, tag3, csr_mat%mp_group)
          CASE (dbcsr_type_complex_4)
            CALL mp_send(csr_mat%nzval_local%c_sp, 0, tag3, csr_mat%mp_group)
          CASE (dbcsr_type_complex_8)
            CALL mp_send(csr_mat%nzval_local%c_dp, 0, tag3, csr_mat%mp_group)
        END SELECT
      ENDIF

      IF (mynode .EQ. 0) THEN ! write data received at node 0
      !WRITE(unit_nr,"(A27)") "#row ind, col ind, value"
        m = 0
        DO k=1, sizes_numrowlocal(i+1)
          rowind = k + SUM(sizes_numrowlocal(1:i)) ! row index: local to global
          DO l=1, nzerow_to_master(k)
            m=m+1
            IF ((.NOT. ut) .OR. (k .LE. colind_to_master(m))) THEN
              SELECT CASE (csr_mat%nzval_local%data_type)
                CASE (dbcsr_type_real_4)
                  WRITE(unit_nr,data_format) rowind, colind_to_master(m), &
                                             nzval_to_master_r_sp(m)
                CASE (dbcsr_type_real_8)
                  WRITE(unit_nr,data_format) rowind, colind_to_master(m), &
                                             nzval_to_master_r_dp(m)
                CASE (dbcsr_type_complex_4)
                  WRITE(unit_nr,data_format) rowind, colind_to_master(m), &
                                             nzval_to_master_c_sp(m)
                CASE (dbcsr_type_complex_8)
                  WRITE(unit_nr,data_format) rowind, colind_to_master(m), &
                                             nzval_to_master_c_dp(m)
              END SELECT
            ENDIF
          ENDDO
        ENDDO

        DEALLOCATE(colind_to_master)
        DEALLOCATE(nzerow_to_master)
  
        SELECT CASE (csr_mat%nzval_local%data_type)
          CASE (dbcsr_type_real_4)
            DEALLOCATE(nzval_to_master_r_sp)
          CASE (dbcsr_type_real_8)
            DEALLOCATE(nzval_to_master_r_dp)
          CASE (dbcsr_type_complex_4)
            DEALLOCATE(nzval_to_master_c_sp)
          CASE (dbcsr_type_complex_8)
            DEALLOCATE(nzval_to_master_c_dp)
        END SELECT
      ENDIF
    ENDDO

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE csr_write
 
! *****************************************************************************
!> \brief Converts a DBCSR matrix to a block row distributed matrix.
!> \param[in] dbcsr_mat        DBCSR matrix to be converted
!> \param[out] brd_mat         converted matrix
!> \param[in] equal_dist       see documentation of csr_create_from_dbcsr
!> \param[in] floor_dist       see documentation of csr_create_from_dbcsr
!> \param[in] numnodes         number of nodes to use for block row distribution
!> \param[in,out] error        error
!>
!> \par History
!>      04.2013 created [Hossein Bani-Hashemian]
!>      11.2014 modified [Patrick Seewald]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dbcsr_create_brd (dbcsr_mat, brd_mat, equal_dist, floor_dist, numnodes, error)

    TYPE(dbcsr_obj), INTENT(IN)              :: dbcsr_mat
    TYPE(dbcsr_obj), INTENT(OUT)             :: brd_mat
    LOGICAL, INTENT(IN)                      :: equal_dist, floor_dist
    INTEGER, INTENT(IN)                      :: numnodes
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_create_brd', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: matrix_type
    CHARACTER(LEN=default_string_length)     :: matrix_name
    INTEGER :: chunk_size, data_type, error_handler, i, k, l, m, mp_group, &
      mynode, nblkcols_total, nblkrows_total, nfullrows_local, &
      nfullrows_total, node_size, numnodes_total, remainder, row_index, &
      split_row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rdist_tmp, &
                                                row_blk_size_new_tmp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: pgrid
    INTEGER, DIMENSION(:), POINTER           :: cdist, rdist, &
                                                row_blk_size_data, &
                                                row_blk_size_new_data
    TYPE(array_i1d_obj)                      :: brd_col_dist, brd_row_dist, &
                                                col_blk_size, row_blk_size, &
                                                row_blk_size_new
    TYPE(dbcsr_distribution_obj)             :: dist_current, dist_new
    TYPE(dbcsr_mp_obj)                       :: mp_obj_current, mp_obj_new

    CALL dbcsr_error_set(routineN, error_handler, error)
 
    NULLIFY(row_blk_size_data, rdist, row_blk_size_new_data)
    CALL dbcsr_get_info(dbcsr_mat,&
                        nblkrows_total=nblkrows_total,&
                        nblkcols_total=nblkcols_total,&
                        nfullrows_total=nfullrows_total,&
                        row_blk_size=row_blk_size,&
                        col_blk_size=col_blk_size,& 
                        matrix_type=matrix_type,&
                        data_type = data_type)

    matrix_name = dbcsr_name(dbcsr_mat)

    ALLOCATE (cdist(nblkcols_total))
    cdist = 0

    CALL array_nullify (brd_row_dist)
    CALL array_nullify (brd_col_dist)

    dist_current = dbcsr_distribution(dbcsr_mat)
    mp_obj_current = dbcsr_distribution_mp(dist_current)
    mp_group = dbcsr_mp_group(mp_obj_current)
    mynode = dbcsr_mp_mynode(mp_obj_current)
    numnodes_total = dbcsr_mp_numnodes(mp_obj_current)

    ALLOCATE (pgrid(numnodes_total,1))

    IF (equal_dist) THEN ! Equally distribute rows over processors -> cut blocks

      row_blk_size_data => array_data(row_blk_size)

      ! Calculate the number of rows a processor can hold
      IF (floor_dist) THEN
        nfullrows_local = FLOOR(REAL(nfullrows_total)/numnodes)
      ELSE 
        nfullrows_local = CEILING(REAL(nfullrows_total)/numnodes)
      ENDIF

      ! allocate maximum amount of memory possibly needed
      ALLOCATE(rdist_tmp(nblkrows_total + numnodes - 1)) ! row distribution
      ALLOCATE(row_blk_size_new_tmp(nblkrows_total + numnodes - 1)) ! new sizes of block rows

      k=0 ! counter for block rows
      m=0 ! node counter
      node_size = nfullrows_local ! space available on current node in number of rows
      IF (node_size .GT. 0) THEN
        DO l=1,nblkrows_total
          split_row = row_blk_size_data(l) ! size of current block row (number of rows)
          DO WHILE(split_row .GE. node_size) ! cut block row and send it to two nodes
            k = k + 1
            m = m + 1
            row_blk_size_new_tmp(k) = node_size ! size of first part of block row
            rdist_tmp(k) = m-1 ! send first part to node m
            split_row = split_row - node_size ! size of remaining part of block rows
            node_size = nfullrows_local ! space available on next node
            IF(floor_dist .AND. (m .EQ. numnodes - 1)) THEN ! send all remaining rows to last node
              node_size = nfullrows_total-(numnodes - 1) * node_size
            ENDIF
          ENDDO
          IF (split_row .GT. 0) THEN ! enough space left on next node for remaining rows
            k = k + 1
            row_blk_size_new_tmp(k) = split_row ! size of remaining part of block row
            rdist_tmp(k) = m ! send to next node
            node_size = node_size - split_row ! remaining space on next node
          ENDIF
        ENDDO
      ELSE ! send everything to last node if node_size = 0 
        rdist_tmp(1:nblkrows_total) = numnodes-1
        row_blk_size_new_tmp(1:nblkrows_total) = row_blk_size_data ! row blocks unchanged
        k = nblkrows_total
      ENDIF

      ! Copy data to correctly allocated variables
      ALLOCATE(row_blk_size_new_data(k))
      row_blk_size_new_data = row_blk_size_new_tmp(1:k)
      CALL array_new(row_blk_size_new,row_blk_size_new_data, gift = .TRUE.)
      ALLOCATE(rdist(k))
      rdist = rdist_tmp(1:k)

    ELSE ! Leave block rows intact (do not cut)
      ALLOCATE (rdist(nblkrows_total)) 
      rdist = 0
      IF (numnodes .GT. nblkrows_total) THEN
        rdist = (/(i, i = 0, nblkrows_total-1)/)
      ELSE
        chunk_size = FLOOR(REAL(nblkrows_total)/numnodes)
        remainder = MOD(nblkrows_total,numnodes)
        row_index = 0
        DO i = 1, chunk_size*numnodes, chunk_size
          rdist(i:i+chunk_size-1) = row_index
          row_index = row_index + 1
        END DO
        IF (remainder .GT. 0) rdist(chunk_size*numnodes+1:nblkrows_total) = row_index-1
      END IF
      row_blk_size_new = row_blk_size
    ENDIF

    CALL array_new (brd_row_dist, rdist, gift = .TRUE.)
    CALL array_new (brd_col_dist, cdist, gift = .TRUE.)

    pgrid(:,:) = RESHAPE((/(i, i = 0, numnodes_total-1)/),(/numnodes_total,1/)) 
    CALL dbcsr_mp_new(mp_obj_new, pgrid, mp_group, mynode, numnodes=numnodes_total)
    CALL dbcsr_distribution_init(dist_new)
    CALL dbcsr_distribution_new(dist_new,mp_obj_new,brd_row_dist,brd_col_dist)

    CALL dbcsr_init(brd_mat)
    CALL dbcsr_create (brd_mat, TRIM(matrix_name)//" row-block distributed",& 
      dist_new, matrix_type, row_blk_size_new, col_blk_size, error=error) 
    CALL dbcsr_complete_redistribute(dbcsr_mat,brd_mat,error=error)
 
    DEALLOCATE (pgrid)

    CALL array_release(brd_row_dist)
    CALL array_release(brd_col_dist)
    
    IF(equal_dist) CALL array_release(row_blk_size_new)   
 
    CALL dbcsr_distribution_release(dist_new)
    CALL dbcsr_mp_release(mp_obj_new)
    
    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_create_brd

END MODULE dbcsr_csr_conversions
