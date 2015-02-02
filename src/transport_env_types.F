!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief CP2K transport environment and related C-interoperable types 
!> \par History
!>       05.2013 created C-interoperable matrices [Hossein Bani-Hashemian]
!>       07.2013 created transport_env [Hossein Bani-Hashemian]
!>       11.2014 revised into CSR matrices [Hossein Bani-Hashemian]
!>       12.2014 merged csr_interop and transport [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE transport_env_types

  USE cp_dbcsr_interface,              ONLY: cp_csr_destroy,&
                                             cp_dbcsr_release,&
                                             cp_dbcsr_type,&
                                             csr_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop

  USE, INTRINSIC :: ISO_C_BINDING,     ONLY: C_DOUBLE, C_INT, &
                                             C_NULL_PTR, C_NULL_FUNPTR, &
                                             C_PTR, C_FUNPTR, &
                                             C_F_POINTER
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'transport_env_types'

  PUBLIC :: transport_env_type, cp2k_transport_parameters
  PUBLIC :: cp2k_csr_interop_type

! DO NOT change the ORDERING or the NAMES in the following data type
  TYPE, BIND(C) :: cp2k_transport_parameters
     INTEGER(C_INT)  :: n_occ
     INTEGER(C_INT)  :: n_atoms
     REAL(C_DOUBLE)  :: evoltfactor
     INTEGER(C_INT)  :: method
     INTEGER(C_INT)  :: bandwidth
     INTEGER(C_INT)  :: n_cells
     INTEGER(C_INT)  :: n_abscissae
     INTEGER(C_INT)  :: n_kpoint
     INTEGER(C_INT)  :: num_interval
     INTEGER(C_INT)  :: num_contacts
     INTEGER(C_INT)  :: ndof
     INTEGER(C_INT)  :: tasks_per_point
     INTEGER(C_INT)  :: cores_per_node
     REAL(C_DOUBLE)  :: colzero_threshold
     REAL(C_DOUBLE)  :: eps_limit
     REAL(C_DOUBLE)  :: eps_decay
     REAL(C_DOUBLE)  :: eps_singularity_curvatures
     REAL(C_DOUBLE)  :: eps_mu
     REAL(C_DOUBLE)  :: eps_eigval_degen
     REAL(C_DOUBLE)  :: energy_interval
     REAL(C_DOUBLE)  :: min_interval
     REAL(C_DOUBLE)  :: temperature
  END TYPE cp2k_transport_parameters

  TYPE transport_env_type
     TYPE(C_FUNPTR)                   :: ext_c_method_ptr
     TYPE(cp2k_transport_parameters)  :: params 
     TYPE(cp_dbcsr_type)              :: template_matrix_sym
     TYPE(cp_dbcsr_type)              :: template_matrix_nosym
     TYPE(csr_type)                   :: s_matrix 
     TYPE(csr_type)                   :: ks_matrix 
     TYPE(csr_type)                   :: p_matrix
  END TYPE transport_env_type

! DO NOT change the ORDERING or the NAMES in the following data type
  TYPE, BIND(C) :: cp2k_csr_interop_type
    INTEGER(C_INT) :: nrows_total
    INTEGER(C_INT) :: ncols_total
    INTEGER(C_INT) :: nze_total
    INTEGER(C_INT) :: nze_local
    INTEGER(C_INT) :: nrows_local
    INTEGER(C_INT) :: data_type
    INTEGER(C_INT) :: first_row
    TYPE(C_PTR)    :: rowptr_local
    TYPE(C_PTR)    :: colind_local
    TYPE(C_PTR)    :: nzerow_local
    TYPE(C_PTR)    :: nzvals_local
  ENDTYPE cp2k_csr_interop_type

  PUBLIC :: csr_interop_nullify, &
            csr_interop_matrix_get_info
  PUBLIC :: transport_env_release

CONTAINS

! *****************************************************************************
!> \brief releases the transport_env
!> \param[inout] transport_env the transport_env to be released
!> \param[inout] error         CP2K error
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE transport_env_release(transport_env,error)
    TYPE(transport_env_type), POINTER        :: transport_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'transport_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.

    CPPostcondition(ASSOCIATED(transport_env),cp_failure_level,routineP,error,failure)

    CALL cp_csr_destroy(transport_env%s_matrix, error)
    CALL cp_csr_destroy(transport_env%ks_matrix, error)
    CALL cp_csr_destroy(transport_env%p_matrix, error)
    CALL cp_dbcsr_release(transport_env%template_matrix_sym, error)
    CALL cp_dbcsr_release(transport_env%template_matrix_nosym, error)

    transport_env%ext_c_method_ptr = C_NULL_FUNPTR

    DEALLOCATE(transport_env, stat=istat)
    CPPostcondition(istat==0, cp_failure_level, routineP, error, failure)

    CALL timestop(handle)

  END SUBROUTINE transport_env_release

! *****************************************************************************
!> \brief nullifies (and zeroizes) a C-interoperable CSR matrix
!> \param[inout] csr_interop_mat the matrix to be nullified
!> \param[inout] error           CP2K error
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE csr_interop_nullify(csr_interop_mat, error)

    TYPE(cp2k_csr_interop_type), &
      INTENT(INOUT)                          :: csr_interop_mat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csr_interop_nullify', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    csr_interop_mat%nrows_total = 0
    csr_interop_mat%ncols_total = 0
    csr_interop_mat%nze_total   = 0
    csr_interop_mat%nze_local   = 0
    csr_interop_mat%nrows_local = 0
    csr_interop_mat%data_type   = 0
    csr_interop_mat%first_row   = 0
    csr_interop_mat%rowptr_local = C_NULL_PTR
    csr_interop_mat%colind_local = C_NULL_PTR
    csr_interop_mat%nzerow_local = C_NULL_PTR
    csr_interop_mat%nzvals_local = C_NULL_PTR

    CALL timestop(handle)

  END SUBROUTINE csr_interop_nullify

! *****************************************************************************
!> \brief gets the fields of a C-interoperable CSR matrix
!> \param[in] csr_interop_mat C-interoperable CSR matrix
!> \param[out] nrows_total     total number of rows
!> \param[out] ncols_total     total number of columns
!> \param[out] nze_local       number of local nonzero elements
!> \param[out] nze_total       total number of nonzero elements
!> \param[out] nrows_local     number of local rows
!> \param[out] data_type       data type
!> \param[out] first_row       index of the first row (C indexing)
!> \param[out] rowptr_local    row pointer (local - Fortran indexing)
!> \param[out] colind_local    column index (local - Fortran indexing)
!> \param[out] nzerow_local    number of nunzeros per row (index-i, local - Fortran indexing)
!> \param[out] nzvals_local    nonzero elements (local)
!> \param[inout] error           CP2K error
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE csr_interop_matrix_get_info(csr_interop_mat, &
                        nrows_total, ncols_total, nze_local, nze_total, nrows_local, data_type, &
                        first_row, rowptr_local, colind_local, nzerow_local, nzvals_local, error)

    TYPE(cp2k_csr_interop_type), INTENT(IN)  :: csr_interop_mat
    INTEGER, INTENT(OUT), OPTIONAL           :: nrows_total, ncols_total, &
                                                nze_local, nze_total, &
                                                nrows_local, data_type, &
                                                first_row
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL, POINTER                      :: rowptr_local, colind_local, &
                                                nzerow_local
    REAL(dp), DIMENSION(:), INTENT(OUT), &
      OPTIONAL, POINTER                      :: nzvals_local
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csr_interop_matrix_get_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    IF (PRESENT(nrows_total)) nrows_total = csr_interop_mat%nrows_total
    IF (PRESENT(ncols_total)) ncols_total = csr_interop_mat%ncols_total
    IF (PRESENT(nze_local)) nze_local = csr_interop_mat%nze_local
    IF (PRESENT(nze_total)) nze_total = csr_interop_mat%nze_total
    IF (PRESENT(nrows_local)) nrows_local = csr_interop_mat%nrows_local
    IF (PRESENT(data_type)) data_type = csr_interop_mat%data_type
    IF (PRESENT(first_row)) first_row = csr_interop_mat%first_row

    IF (PRESENT(rowptr_local)) CALL C_F_POINTER(csr_interop_mat%rowptr_local, rowptr_local, [nrows_local+1])
    IF (PRESENT(colind_local)) CALL C_F_POINTER(csr_interop_mat%colind_local, colind_local, [nze_local])
    IF (PRESENT(nzerow_local)) CALL C_F_POINTER(csr_interop_mat%nzerow_local, nzerow_local, [nrows_local])
    IF (PRESENT(nzvals_local)) CALL C_F_POINTER(csr_interop_mat%nzvals_local, nzvals_local, [nze_local])

    CALL timestop(handle)

  END SUBROUTINE csr_interop_matrix_get_info

END MODULE transport_env_types


