!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines for DFT+NEGF calculations (coupling with the quantum transport code OMEN)
!> \par History
!>       12.2012 created external_scf_method [Hossein Bani-Hashemian]
!>       05.2013 created rotines to work with C-interoperable matrices [Hossein Bani-Hashemian]
!>       07.2013 created transport_env routines [Hossein Bani-Hashemian]
!>       11.2014 switch to CSR matrices [Hossein Bani-Hashemian]
!>       12.2014 merged [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE transport

  USE cp_dbcsr_interface,              ONLY: &
       cp_convert_csr_to_dbcsr, cp_convert_dbcsr_to_csr, cp_csr_create, &
       cp_csr_create_from_dbcsr, cp_dbcsr_copy, cp_dbcsr_copy_into_existing, &
       cp_dbcsr_desymmetrize, cp_dbcsr_has_symmetry, cp_dbcsr_init, &
       cp_dbcsr_type, csr_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_allgather
  USE physcon,                         ONLY: evolt
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE transport_env_types,             ONLY: cp2k_csr_interop_type,&
                                             cp2k_transport_parameters,&
                                             csr_interop_matrix_get_info,&
                                             csr_interop_nullify,&
                                             transport_env_type

  USE, INTRINSIC :: ISO_C_BINDING,     ONLY: C_ASSOCIATED, C_LOC, C_F_PROCPOINTER, &
                                             C_PTR, C_NULL_PTR, C_DOUBLE, C_INT
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'transport'

  PUBLIC :: transport_env_create, transport_initialize, external_scf_method 

!> interface between C/C++ and FORTRAN
  INTERFACE c_func_interface
! *****************************************************************************
!> \brief C routine that takes the S and H matrices as input and outputs a P matrix
!> \param cp2k_transport_params transport parameters read form a CP2K input file
!> \param s_mat  C-interoperable overlap matrix 
!> \param ks_mat C-interoperable Kohn-Sham matrix
!> \param p_mat  C-interoperable density matrix
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
    SUBROUTINE c_scf_routine(cp2k_transport_params, s_mat, ks_mat, p_mat) BIND(C)
      IMPORT :: C_INT, C_PTR, cp2k_csr_interop_type, cp2k_transport_parameters
      IMPLICIT NONE
      TYPE(cp2k_transport_parameters), VALUE, INTENT(IN) :: cp2k_transport_params
      TYPE(cp2k_csr_interop_type), VALUE, INTENT(IN)     :: s_mat
      TYPE(cp2k_csr_interop_type), VALUE, INTENT(IN)     :: ks_mat
      TYPE(cp2k_csr_interop_type), INTENT(INOUT)         :: p_mat
    END SUBROUTINE c_scf_routine
  END INTERFACE c_func_interface

CONTAINS

! *****************************************************************************
!> \brief creates the transport environment
!> \param[inout] qs_env the qs_env containing the transport_env
!> \param[inout] error  CP2K error
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE transport_env_create(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'transport_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: input
    TYPE(transport_env_type), POINTER        :: transport_env

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env, transport_env=transport_env, input=input, error=error)

    CPPostcondition(.NOT.ASSOCIATED(transport_env),cp_failure_level,routineP,error,failure)

    ALLOCATE(transport_env,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL transport_init_read_input(input,transport_env,error)
    CALL cp_dbcsr_init(transport_env%template_matrix_sym, error)
    CALL cp_dbcsr_init(transport_env%template_matrix_nosym, error)
    CALL set_qs_env(qs_env,transport_env=transport_env,error=error)

    CALL timestop(handle)

  END SUBROUTINE transport_env_create

! *****************************************************************************
!> \brief intitializes all fields of transport_env using the parameters read from 
!>        the corresponding input section 
!> \param[inout] input         the input file 
!> \param[inout] transport_env the transport_env to be initialized
!> \param[inout] error         CP2K error
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE transport_init_read_input(input,transport_env,error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(transport_env_type), INTENT(INOUT)  :: transport_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'transport_init_read_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(section_vals_type), POINTER         :: transport_section

    CALL timeset(routineN,handle)

    transport_section => section_vals_get_subs_vals(input,"DFT%TRANSPORT",error=error)
    
!> read from input
    CALL section_vals_val_get(transport_section,"DENSITY_MATRIX_CONSTRUCTION_METHOD",i_val=transport_env%params%method,error=error)
    CALL section_vals_val_get(transport_section,"BANDWIDTH",i_val=transport_env%params%bandwidth,error=error)
    CALL section_vals_val_get(transport_section,"N_CELLS",i_val=transport_env%params%n_cells,error=error)
    CALL section_vals_val_get(transport_section,"N_ABSCISSAE",i_val=transport_env%params%n_abscissae,error=error)
    CALL section_vals_val_get(transport_section,"N_KPOINTS",i_val=transport_env%params%n_kpoint,error=error)
    CALL section_vals_val_get(transport_section,"NUM_INTERVAL",i_val=transport_env%params%num_interval,error=error)
    CALL section_vals_val_get(transport_section,"NUM_CONTACTS",i_val=transport_env%params%num_contacts,error=error)
    CALL section_vals_val_get(transport_section,"N_DOF",i_val=transport_env%params%ndof,error=error)
    CALL section_vals_val_get(transport_section,"TASKS_PER_POINT",i_val=transport_env%params%tasks_per_point,error=error)
    CALL section_vals_val_get(transport_section,"CORES_PER_NODE",i_val=transport_env%params%cores_per_node,error=error)
    CALL section_vals_val_get(transport_section,"COLZERO_THRESHOLD",r_val=transport_env%params%colzero_threshold,error=error)
    CALL section_vals_val_get(transport_section,"EPS_LIMIT",r_val=transport_env%params%eps_limit,error=error)
    CALL section_vals_val_get(transport_section,"EPS_DECAY",r_val=transport_env%params%eps_decay,error=error)
    CALL section_vals_val_get(transport_section,"EPS_SINGULARITY_CURVATURES",&
                                                r_val=transport_env%params%eps_singularity_curvatures,error=error)
    CALL section_vals_val_get(transport_section,"EPS_MU",r_val=transport_env%params%eps_mu,error=error)
    CALL section_vals_val_get(transport_section,"EPS_EIGVAL_DEGEN",r_val=transport_env%params%eps_eigval_degen,error=error)
    CALL section_vals_val_get(transport_section,"ENERGY_INTERVAL",r_val=transport_env%params%energy_interval,error=error)
    CALL section_vals_val_get(transport_section,"MIN_INTERVAL",r_val=transport_env%params%min_interval,error=error)
    CALL section_vals_val_get(transport_section,"TEMPERATURE",r_val=transport_env%params%temperature,error=error)

    CALL timestop(handle)

  END SUBROUTINE transport_init_read_input

! *****************************************************************************
!> \brief initializes the transport environment
!> \param[inout] transport_env the transport env to be initialized 
!> \param[in]    template_matrix   template matrix to keep the sparsity of matrices fixed
!> \param[inout] error             CP2K error
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE transport_initialize(transport_env, template_matrix, error)
    TYPE(transport_env_type), INTENT(INOUT)  :: transport_env
    TYPE(cp_dbcsr_type), INTENT(IN)          :: template_matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'transport_initialize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)  

    IF (cp_dbcsr_has_symmetry(template_matrix)) THEN
       CALL cp_dbcsr_copy(transport_env%template_matrix_sym, template_matrix, error=error)
       CALL cp_dbcsr_desymmetrize(transport_env%template_matrix_sym, transport_env%template_matrix_nosym, error)
    ELSE 
       CALL cp_dbcsr_copy(transport_env%template_matrix_nosym, template_matrix, error=error)
       CALL cp_dbcsr_copy(transport_env%template_matrix_sym, template_matrix, error=error)
    END IF

    CALL cp_csr_create_from_dbcsr(transport_env%template_matrix_nosym, transport_env%s_matrix, 1, error=error)
    CALL cp_convert_dbcsr_to_csr(transport_env%template_matrix_nosym, transport_env%s_matrix, error=error)

    CALL cp_csr_create(transport_env%ks_matrix, transport_env%s_matrix, error)
    CALL cp_csr_create(transport_env%p_matrix , transport_env%s_matrix, error)

    CALL timestop(handle)
 
  END SUBROUTINE transport_initialize

! *****************************************************************************
!> \brief SCF calcualtion with an externally evaluated density matrix
!> \param[inout] transport_env  transport environment 
!> \param[in]    matrix_s       DBCSR overlap matrix
!> \param[in]    matrix_ks      DBCSR Kohn-Sham matrix
!> \param[inout] matrix_p       DBCSR density matrix
!> \param[in]    nelectron_spin number of electrons 
!> \param[in]    natoms         number of atoms 
!> \param[inout] error          CP2K error
!> \par History
!>       12.2012 created [Hossein Bani-Hashemian]
!>       12.2014 revised [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
   SUBROUTINE external_scf_method(transport_env, matrix_s, matrix_ks, matrix_p, nelectron_spin, natoms, error)

    TYPE(transport_env_type), INTENT(INOUT)  :: transport_env
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_s, matrix_ks
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_p
    INTEGER, INTENT(IN)                      :: nelectron_spin, natoms
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'external_scf_method', &
      routineP = moduleN//':'//routineN

    TYPE(cp2k_csr_interop_type)              :: ks_mat, p_mat, s_mat

    PROCEDURE(c_scf_routine), POINTER        :: c_method
    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    CALL C_F_PROCPOINTER(transport_env%ext_c_method_ptr, c_method)
    CALL cp_assert(C_ASSOCIATED(transport_env%ext_c_method_ptr),cp_fatal_level,cp_assertion_failed,routineP,&
                  "MISSING C/C++ ROUTINE: The TRANSPORT section is meant to be used together with an external "//&
                  "program, e.g. the quantum transport code OMEN, that provides CP2K with a density matrix.",&
                  error, failure)

    transport_env%params%n_occ = nelectron_spin
    transport_env%params%n_atoms = natoms
    transport_env%params%evoltfactor = evolt

    CALL csr_interop_nullify(s_mat, error)
    CALL csr_interop_nullify(ks_mat, error)
    CALL csr_interop_nullify(p_mat, error)

    CALL cp_dbcsr_copy_into_existing(transport_env%template_matrix_sym, matrix_s, error) 
    CALL convert_dbcsr_to_csr_interop(transport_env%template_matrix_sym, transport_env%s_matrix, s_mat, error=error)

    CALL cp_dbcsr_copy_into_existing(transport_env%template_matrix_sym, matrix_ks, error)
    CALL convert_dbcsr_to_csr_interop(transport_env%template_matrix_sym, transport_env%ks_matrix, ks_mat, error=error)

    CALL cp_dbcsr_copy_into_existing(transport_env%template_matrix_sym, matrix_s, error)
    CALL convert_dbcsr_to_csr_interop(transport_env%template_matrix_sym, transport_env%p_matrix, p_mat, error=error)

    CALL c_method(transport_env%params, s_mat, ks_mat, p_mat)

    CALL convert_csr_interop_to_dbcsr(p_mat, transport_env%p_matrix, transport_env%template_matrix_nosym, error)
    CALL cp_dbcsr_copy(matrix_p, transport_env%template_matrix_nosym, error=error)

    CALL timestop(handle)

  END SUBROUTINE external_scf_method

! *****************************************************************************
!> \brief converts a DBCSR matrix to a C-interoperable CSR matrix
!> \param[in]    dbcsr_mat  DBCSR matrix to be converted
!> \param[inout] csr_mat    auxiliary CSR matrix
!> \param[inout] csr_interop_mat C-interoperable CSR matrix
!> \param[inout] error      CP2K error
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE convert_dbcsr_to_csr_interop(dbcsr_mat, csr_mat, csr_interop_mat, error)

    TYPE(cp_dbcsr_type), INTENT(IN)          :: dbcsr_mat
    TYPE(csr_type), INTENT(INOUT)            :: csr_mat
    TYPE(cp2k_csr_interop_type), INTENT(INOUT)      :: csr_interop_mat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'convert_dbcsr_to_csr_interop', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, mepos, num_pe, mp_group
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nrows_local_all, first_row_all
    INTEGER(C_INT), DIMENSION(:), POINTER    :: colind_local, rowptr_local, nzerow_local
    REAL(C_DOUBLE), DIMENSION(:), POINTER    :: nzvals_local
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    mp_group = logger%para_env%group
    mepos = logger%para_env%mepos
    num_pe = logger%para_env%num_pe

!> dbcsr to csr
    CALL cp_convert_dbcsr_to_csr(dbcsr_mat, csr_mat, error=error)

!> csr to csr_interop
    rowptr_local => csr_mat%rowptr_local
    colind_local => csr_mat%colind_local
    nzerow_local => csr_mat%nzerow_local
    nzvals_local => csr_mat%nzval_local%r_dp !> support real double percision for now

    IF (SIZE(rowptr_local) .EQ. 0) THEN
       csr_interop_mat%rowptr_local = C_NULL_PTR
    ELSE
       csr_interop_mat%rowptr_local = C_LOC(rowptr_local(1))
    END IF

    IF (SIZE(colind_local) .EQ. 0) THEN
       csr_interop_mat%colind_local = C_NULL_PTR
    ELSE
       csr_interop_mat%colind_local = C_LOC(colind_local(1))
    END IF

    IF (SIZE(nzerow_local) .EQ. 0) THEN
       csr_interop_mat%nzerow_local = C_NULL_PTR
    ELSE
       csr_interop_mat%nzerow_local = C_LOC(nzerow_local(1))
    END IF

    IF (SIZE(nzvals_local) .EQ. 0) THEN
       csr_interop_mat%nzvals_local = C_NULL_PTR
    ELSE
       csr_interop_mat%nzvals_local = C_LOC(nzvals_local(1))
    END IF

    ALLOCATE(nrows_local_all(0:num_pe-1), first_row_all(0:num_pe-1))
    CALL mp_allgather(csr_mat%nrows_local, nrows_local_all, mp_group)
    CALL cumsum_i(nrows_local_all, first_row_all)

    IF (mepos .EQ. 0) THEN
       csr_interop_mat%first_row = 0
    ELSE
       csr_interop_mat%first_row = first_row_all(mepos - 1)
    END IF
    csr_interop_mat%nrows_total = csr_mat%nrows_total
    csr_interop_mat%ncols_total = csr_mat%ncols_total
    csr_interop_mat%nze_local = csr_mat%nze_local
    IF (csr_mat%nze_total>HUGE(csr_interop_mat%nze_total)) THEN
        ! overflow in nze
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,CPSourceFileRef,&
             error, failure)
    ENDIF
    csr_interop_mat%nze_total = INT(csr_mat%nze_total,KIND=KIND(csr_interop_mat%nze_total))
    csr_interop_mat%nrows_local = csr_mat%nrows_local
    csr_interop_mat%data_type = csr_mat%nzval_local%data_type

    CALL timestop(handle)

    CONTAINS
! *****************************************************************************
!> \brief cumulative sum of a 1d array of integers
!> \param[in]  arr    input array
!> \param[out] cumsum cumulative sum of the input array
! *****************************************************************************
    SUBROUTINE cumsum_i(arr,cumsum)
    INTEGER, DIMENSION(:), INTENT(IN)        :: arr
    INTEGER, DIMENSION(SIZE(arr)), &
      INTENT(OUT)                            :: cumsum

    INTEGER                                  :: i

       cumsum(1)=arr(1)
       DO i=2,SIZE(arr)
          cumsum(i)=cumsum(i-1)+arr(i)
       END DO
    END SUBROUTINE cumsum_i

  END SUBROUTINE convert_dbcsr_to_csr_interop

! *****************************************************************************
!> \brief converts a C-interoperable CSR matrix to a DBCSR matrix
!> \param[in] csr_interop_mat C-interoperable CSR matrix
!> \param[inout] csr_mat         auxiliary CSR matrix
!> \param[inout] dbcsr_mat       DBCSR matrix
!> \param[inout] error           CP2K error
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE convert_csr_interop_to_dbcsr(csr_interop_mat, csr_mat, dbcsr_mat, error)

    TYPE(cp2k_csr_interop_type), INTENT(IN)  :: csr_interop_mat
    TYPE(csr_type), INTENT(INOUT)            :: csr_mat
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: dbcsr_mat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'convert_csr_interop_to_dbcsr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_type, handle, &
                                                ncols_total, nrows_local, &
                                                nrows_total, nze_local, &
                                                nze_total
    INTEGER, DIMENSION(:), POINTER           :: colind_local, nzerow_local, &
                                                rowptr_local
    REAL(dp), DIMENSION(:), POINTER          :: nzvals_local

    CALL timeset(routineN,handle)

!> csr_interop to csr
    CALL csr_interop_matrix_get_info(csr_interop_mat, &
                   nrows_total = nrows_total, ncols_total = ncols_total, nze_local = nze_local, &
                   nze_total = nze_total, nrows_local = nrows_local, data_type = data_type,&
                   rowptr_local = rowptr_local, colind_local = colind_local, &
                   nzerow_local = nzerow_local, nzvals_local = nzvals_local, &
                   error=error)

    csr_mat%nrows_total = nrows_total
    csr_mat%ncols_total = ncols_total
    csr_mat%nze_local = nze_local
    csr_mat%nze_total = nze_total
    csr_mat%nrows_local = nrows_local
    csr_mat%nzval_local%data_type = data_type

    csr_mat%rowptr_local = rowptr_local
    csr_mat%colind_local = colind_local
    csr_mat%nzerow_local = nzerow_local
    csr_mat%nzval_local%r_dp = nzvals_local

!> csr to dbcsr
    CALL cp_convert_csr_to_dbcsr(dbcsr_mat, csr_mat, error=error)

    CALL timestop(handle)

  END SUBROUTINE convert_csr_interop_to_dbcsr

END MODULE transport


