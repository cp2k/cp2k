!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utility functions for the perturbation calculations.
!> \note
!>      - routines are programmed with spins in mind
!>        but are as of now not tested with them
!> \par History
!>      22-08-2002, TCH, started development
! *****************************************************************************
MODULE qs_p_env_methods
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_allocate_matrix_set,&
                                             cp_dbcsr_copy,&
                                             cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_scale,&
                                             cp_dbcsr_set,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                             cp_fm_symm,&
                                             cp_fm_triangular_multiply
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             cp_fm_pool_type,&
                                             fm_pool_create_fm,&
                                             fm_pool_get_el_struct,&
                                             fm_pool_give_back_fm,&
                                             fm_pools_create_fm_vect
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_retain,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_copy,&
                                             cp_fm_vect_dealloc
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE hartree_local_methods,           ONLY: init_coulomb_local
  USE hartree_local_types,             ONLY: hartree_local_create
  USE input_constants,                 ONLY: ot_precond_none
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE preconditioner_types,            ONLY: init_preconditioner
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kpp1_env_methods,             ONLY: kpp1_calc_k_p_p1,&
                                             kpp1_calc_k_p_p1_fdiff,&
                                             kpp1_create,&
                                             kpp1_did_change
  USE qs_kpp1_env_types,               ONLY: qs_kpp1_env_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_linres_types,                 ONLY: linres_control_type
  USE qs_local_rho_types,              ONLY: local_rho_set_create
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_rho0_methods,                 ONLY: init_rho0
  USE qs_rho_atom_methods,             ONLY: allocate_rho_atom_internals
  USE qs_rho_methods,                  ONLY: qs_rho_rebuild,&
                                             qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_type
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_p_env_methods'
  INTEGER, PRIVATE, SAVE :: last_p_env_id=0
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.

  PRIVATE
  PUBLIC :: p_env_create, p_env_psi0_changed, p_op_l2_fawzi
  PUBLIC :: p_op_l1, p_op_l2, p_preortho, p_postortho,&
            p_env_did_change

CONTAINS

! *****************************************************************************
!> \brief allocates and initializes the perturbation environment (no setup)
!> \param p_env the environment to initialize
!> \param qs_env the qs_environment for the system
!> \param kpp1_env the environment that builds the second order
!>        perturbation kernel
!> \param p1_option ...
!> \param psi0d ...
!> \param orthogonal_orbitals if the orbitals are orthogonal
!> \param linres_control ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE p_env_create(p_env, qs_env, kpp1_env, p1_option, &
       psi0d, orthogonal_orbitals,linres_control, error)

    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_kpp1_env_type), OPTIONAL, &
      POINTER                                :: kpp1_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: p1_option
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: psi0d
    LOGICAL, INTENT(in), OPTIONAL            :: orthogonal_orbitals
    TYPE(linres_control_type), OPTIONAL, &
      POINTER                                :: linres_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n_ao, n_mo, n_spins, &
                                                natom, spin, stat
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools, mo_mo_fm_pools
    TYPE(cp_fm_type), POINTER                :: qs_env_c

! code

    CALL timeset(routineN,handle)

    n_spins = qs_env%dft_control%nspins

    failure=.FALSE.
    NULLIFY( ao_mo_fm_pools, mo_mo_fm_pools, matrix_s)

    CALL get_qs_env(qs_env, matrix_s=matrix_s, error=error)

    ALLOCATE(p_env,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       NULLIFY(p_env%kpp1, &
            p_env%p1, &
            p_env%m_epsilon, &
            p_env%psi0d, &
            p_env%S_psi0, &
            p_env%kpp1_env,&
            p_env%rho1,&
            p_env%rho1_xc,&
            p_env%Smo_inv,&
            p_env%local_rho_set, &
            p_env%hartree_local,&
            p_env%PS_psi0,&
            p_env%preconditioner,&
            p_env%ev_h0)

       p_env%ref_count=1
       last_p_env_id=last_p_env_id+1
       p_env%id_nr=last_p_env_id
       p_env%iter=0

       p_env%new_preconditioner = .TRUE.
       p_env%only_energy = .FALSE.
       p_env%os_valid = .FALSE.
       p_env%ls_count = 0
       p_env%delta=0.0_dp
       p_env%gnorm=0.0_dp
       p_env%gnorm_old=0.0_dp
       p_env%etotal=0.0_dp
       p_env%gradient=0.0_dp

       IF (PRESENT(kpp1_env)) THEN
          p_env%kpp1_env => kpp1_env
       ELSE
          CALL kpp1_create(p_env%kpp1_env, qs_env=qs_env, error=error)
       END IF

       IF (PRESENT(p1_option)) THEN
          p_env%p1 => p1_option
       ELSE

          CALL cp_dbcsr_allocate_matrix_set(p_env%p1,n_spins,error)
          DO spin=1, n_spins
             ALLOCATE(p_env%p1(spin)%matrix)
             CALL cp_dbcsr_init(p_env%p1(spin)%matrix,error=error)
             CALL cp_dbcsr_copy(p_env%p1(spin)%matrix,matrix_s(1)%matrix,&
                  name="p_env"//TRIM(ADJUSTL(cp_to_string(p_env%id_nr)))//&
                  "%p1-"//TRIM(ADJUSTL(cp_to_string(spin))),error=error)
             CALL cp_dbcsr_set(p_env%p1(spin)%matrix,0.0_dp,error=error)
          END DO
       END IF

       CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools,&
            mo_mo_fm_pools=mo_mo_fm_pools,&
            error=error)

       p_env%n_mo=0
       p_env%n_ao=0
       DO spin=1, n_spins
          IF (PRESENT(psi0d)) THEN
             CALL cp_fm_get_info(psi0d(spin)%matrix, &
                  ncol_global=n_mo, nrow_global=n_ao,error=error)
          ELSE
             CALL get_mo_set(qs_env%mos(spin)%mo_set, mo_coeff=qs_env_c)
             CALL cp_fm_get_info(qs_env_c, &
                  ncol_global=n_mo, nrow_global=n_ao,error=error)
          END IF
          p_env%n_mo(spin) = n_mo
          p_env%n_ao(spin) = n_ao
       END DO

       p_env%orthogonal_orbitals=.FALSE.
       IF (PRESENT(orthogonal_orbitals)) &
            p_env%orthogonal_orbitals=orthogonal_orbitals

       CALL fm_pools_create_fm_vect(ao_mo_fm_pools,elements=p_env%S_psi0,&
            name="p_env"//TRIM(ADJUSTL(cp_to_string(p_env%id_nr)))//"%S_psi0",&
            error=error)

       ! alloc m_epsilon
       CALL fm_pools_create_fm_vect(mo_mo_fm_pools,elements=p_env%m_epsilon,&
            name="p_env"//TRIM(ADJUSTL(cp_to_string(p_env%id_nr)))//&
            "%m_epsilon",&
            error=error)

       ! alloc Smo_inv
       IF (.NOT. p_env%orthogonal_orbitals) THEN
          CALL fm_pools_create_fm_vect(mo_mo_fm_pools,elements=p_env%Smo_inv,&
               name="p_env"//TRIM(ADJUSTL(cp_to_string(p_env%id_nr)))//&
               "%Smo_inv",&
               error=error)
       END IF

       IF (PRESENT(psi0d)) THEN
          IF (ASSOCIATED(psi0d)) THEN
             CALL cp_fm_vect_copy(psi0d,p_env%psi0d,error=error)
          END IF
       ELSE IF (.NOT.p_env%orthogonal_orbitals) THEN
          CALL fm_pools_create_fm_vect(ao_mo_fm_pools,&
               elements=p_env%psi0d,&
               name="p_env"//TRIM(ADJUSTL(cp_to_string(p_env%id_nr)))//&
               "%psi0d", &
               error=error)
       END IF

       !----------------------!
       ! GAPW initializations !
       !----------------------!
       IF (qs_env%dft_control%qs_control%gapw) THEN
          CALL local_rho_set_create(p_env%local_rho_set,error=error)
          CALL allocate_rho_atom_internals(qs_env, p_env%local_rho_set%rho_atom_set,error=error)
          CALL init_rho0(qs_env, qs_env%dft_control%qs_control%gapw_control, &
                         .TRUE., p_env%local_rho_set, error)
          CALL hartree_local_create(p_env%hartree_local,error=error)
          CALL get_qs_env(qs_env=qs_env,natom=natom,error=error)
          CALL init_coulomb_local(p_env%hartree_local,natom,error=error)
       END IF

       !------------------------!
       ! LINRES initializations !
       !------------------------!
       IF(PRESENT(linres_control)) THEN

         IF(linres_control%preconditioner_type/= ot_precond_none) THEN
         ! Initialize the preconditioner matrix
           IF (.NOT.ASSOCIATED(p_env%preconditioner)) THEN

             ALLOCATE(p_env%preconditioner(n_spins), STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO spin =1,n_spins
               CALL init_preconditioner(p_env%preconditioner(spin),&
                    para_env=qs_env%para_env,blacs_env=qs_env%blacs_env,error=error)
             END DO
             p_env%os_valid = .FALSE.

             CALL fm_pools_create_fm_vect(ao_mo_fm_pools,elements=p_env%PS_psi0,&
                  name="p_env"//TRIM(ADJUSTL(cp_to_string(p_env%id_nr)))//"%PS_psi0",&
                  error=error)

             ALLOCATE(p_env%ev_h0(n_spins),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!             CALL fm_pools_create_fm_vect(ao_mo_fm_pools,elements=p_env%ev_h0,&
!                  name="p_env"//TRIM(ADJUSTL(cp_to_string(p_env%id_nr)))//"%ev_h0",&
!                  error=error)
           END IF
         END IF

       END IF

    END IF

    CALL timestop(handle)

  END SUBROUTINE p_env_create

! *****************************************************************************
!> \brief checks that the intenal storage is allocated, and allocs it if needed
!> \param p_env the environment to check
!> \param qs_env the qs environment this p_env lives in
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      private routine
! *****************************************************************************
  SUBROUTINE p_env_check_i_alloc(p_env, qs_env, error)
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_env_check_i_alloc', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=25)                        :: name
    INTEGER                                  :: handle, ispin, nspins
    LOGICAL                                  :: failure, gapw_xc
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(dft_control,matrix_s)

    CPPrecondition(ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(p_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env, dft_control=dft_control,error=error)
       gapw_xc = dft_control%qs_control%gapw_xc
       IF (.NOT.ASSOCIATED(p_env%kpp1)) THEN
          CALL get_qs_env(qs_env, matrix_s=matrix_s, error=error)
          nspins=dft_control%nspins

          CALL cp_dbcsr_allocate_matrix_set(p_env%kpp1,nspins,error)
             name="p_env"//cp_to_string(p_env%id_nr)//"%kpp1-"
             CALL compress(name,full=.TRUE.)
             DO ispin=1,nspins
                ALLOCATE(p_env%kpp1(ispin)%matrix)
                CALL cp_dbcsr_init(p_env%kpp1(ispin)%matrix,error=error)
                CALL cp_dbcsr_copy(p_env%kpp1(ispin)%matrix,matrix_s(1)%matrix,&
                     name=TRIM(name)//ADJUSTL(cp_to_string(ispin)),error=error)
             END DO

          CALL qs_rho_rebuild(p_env%rho1, qs_env=qs_env, error=error)
          IF(gapw_xc) THEN
             CALL qs_rho_rebuild(p_env%rho1_xc,qs_env=qs_env,error=error)
          END IF

       END IF

       IF (.NOT.ASSOCIATED(p_env%rho1)) THEN
          CALL qs_rho_rebuild(p_env%rho1, qs_env=qs_env, error=error)
          IF(gapw_xc) THEN
             CALL qs_rho_rebuild(p_env%rho1_xc,qs_env=qs_env,error=error)
          END IF
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE p_env_check_i_alloc

! *****************************************************************************
!> \brief function to advise of changes either in the grids or in the the
!>      matrix structure (to update p_env's internal cache)
!> \param p_env the p_env to inform of the changes
!> \param s_struct_changed true if the structure of the s matrix has changed
!> \param grid_changed true if the grids have changed and have to be rebuilt
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE p_env_did_change(p_env, s_struct_changed, grid_changed,error)
    TYPE(qs_p_env_type), POINTER             :: p_env
    LOGICAL, INTENT(in), OPTIONAL            :: s_struct_changed, grid_changed
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_env_did_change', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(p_env%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.not.failure) THEN
     IF (PRESENT(s_struct_changed)) THEN
        IF (s_struct_changed) THEN
           IF (ASSOCIATED(p_env%kpp1)) THEN
              CALL cp_dbcsr_deallocate_matrix_set(p_env%kpp1,error=error)
           END IF
        END IF
     END IF

     IF (PRESENT(grid_changed)) THEN
        IF (grid_changed) THEN
           IF(ASSOCIATED(p_env%rho1_xc)) THEN
             CALL qs_rho_release(p_env%rho1_xc,error=error)
           END IF
           CALL qs_rho_release(p_env%rho1,error=error)
        END IF
     END IF

     CPPrecondition(ASSOCIATED(p_env%kpp1_env),cp_failure_level,routineP,error,failure)
  END IF
  IF (.not.failure) THEN
     CALL kpp1_did_change(p_env%kpp1_env,s_struct_changed=s_struct_changed,&
          grid_changed=grid_changed,error=error)
  END IF
END SUBROUTINE p_env_did_change

! *****************************************************************************
!> \brief To be called after the value of psi0 has changed.
!>      Recalculates the quantities S_psi0 and m_epsilon.
!> \param p_env the perturbation environment to set
!> \param qs_env ...
!> \param psi0 the value of psi0, if not given defaults to the qs_env mos
!> \param Hrho_psi0d is given, then the partial result Hrho_psi0d is stored in
!>        that vector
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE p_env_psi0_changed(p_env, qs_env, psi0, Hrho_psi0d, error)

    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: psi0
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(inout), OPTIONAL                :: Hrho_psi0d
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_env_psi0_changed', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, lfomo, n_spins, nmo, &
                                                output_unit, spin, stat
    LOGICAL                                  :: failure, was_present
    REAL(KIND=dp)                            :: maxocc
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: my_psi0
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: input, lr_section

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(ao_mo_fm_pools,mos,my_psi0,matrix_s,mos, para_env,ks_env, rho, &
         logger, input,lr_section, energy, matrix_ks)
    logger => cp_error_get_logger(error)
    n_spins = qs_env%dft_control%nspins

    CPPrecondition(ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(p_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       p_env%iter=p_env%iter+1

       CALL get_qs_env(qs_env,&
                       ks_env=ks_env,&
                       mos=mos,&
                       matrix_s=matrix_s,&
                       matrix_ks=matrix_ks,&
                       para_env=para_env,&
                       rho=rho,&
                       input=input,&
                       energy=energy,&
                       error=error)

       CALL mpools_get(qs_env%mpools,&
            ao_mo_fm_pools=ao_mo_fm_pools, error=error)
       ! def my_psi0
       IF (PRESENT(psi0)) THEN
          CALL cp_fm_vect_copy(psi0,my_psi0,error=error)
       ELSE
          ALLOCATE(my_psi0(n_spins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (.NOT.failure) THEN
             DO spin=1,n_spins
                NULLIFY(my_psi0(spin)%matrix)
                CALL get_mo_set(mos(spin)%mo_set,&
                     mo_coeff=my_psi0(spin)%matrix)
                CALL cp_fm_retain(my_psi0(spin)%matrix,error=error)
             END DO
          END IF
       END IF

    END IF

    IF (.NOT.failure) THEN
       lr_section => section_vals_get_subs_vals(input,"PROPERTIES%LINRES",error=error)
       ! def psi0d
       IF (p_env%orthogonal_orbitals) THEN
          IF (ASSOCIATED(p_env%psi0d)) THEN
             CALL cp_fm_vect_dealloc(p_env%psi0d,error=error)
          END IF
          p_env%psi0d => my_psi0
       ELSE

          DO spin=1,n_spins
             ! m_epsilon=choleski_decomposition(my_psi0^T S my_psi0)^-1
             ! could be optimized by combining next two calls
             CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix,&
                  my_psi0(spin)%matrix,&
                  p_env%S_psi0(spin)%matrix,&
                  ncol=p_env%n_mo(spin), alpha=1.0_dp,error=error)
             CALL cp_gemm(transa='T',transb='N',n=p_env%n_mo(spin),&
                  m=p_env%n_mo(spin),k=p_env%n_ao(spin),alpha=1.0_dp,&
                  matrix_a=my_psi0(spin)%matrix, &
                  matrix_b=p_env%S_psi0(spin)%matrix,&
                  beta=0.0_dp,matrix_c=p_env%m_epsilon(spin)%matrix,&
                  error=error)
             CALL cp_fm_cholesky_decompose(p_env%m_epsilon(spin)%matrix,&
                  n=p_env%n_mo(spin),error=error)

             ! Smo_inv= (my_psi0^T S my_psi0)^-1
             CALL cp_fm_set_all(p_env%Smo_inv(spin)%matrix,0.0_dp,1.0_dp,error=error)
             ! faster using cp_fm_cholesky_invert ?
             CALL cp_fm_triangular_multiply(&
                  triangular_matrix=p_env%m_epsilon(spin)%matrix,&
                  matrix_b=p_env%Smo_inv(spin)%matrix,side='R',&
                  invert_tr=.TRUE., n_rows=p_env%n_mo(spin),&
                  n_cols=p_env%n_mo(spin),error=error)
             CALL cp_fm_triangular_multiply(&
                  triangular_matrix=p_env%m_epsilon(spin)%matrix,&
                  matrix_b=p_env%Smo_inv(spin)%matrix,side='R',&
                  transpose_tr=.TRUE.,&
                  invert_tr=.TRUE., n_rows=p_env%n_mo(spin),&
                  n_cols=p_env%n_mo(spin),error=error)

             ! psi0d=my_psi0 (my_psi0^T S my_psi0)^-1
             ! faster using cp_fm_cholesky_invert ?
             CALL cp_fm_to_fm(my_psi0(spin)%matrix,&
                  p_env%psi0d(spin)%matrix, error=error)
             CALL cp_fm_triangular_multiply(&
                  triangular_matrix=p_env%m_epsilon(spin)%matrix,&
                  matrix_b=p_env%psi0d(spin)%matrix,side='R',&
                  invert_tr=.TRUE., n_rows=p_env%n_ao(spin),&
                  n_cols=p_env%n_mo(spin),error=error)
             CALL cp_fm_triangular_multiply(&
                  triangular_matrix=p_env%m_epsilon(spin)%matrix,&
                  matrix_b=p_env%psi0d(spin)%matrix,side='R',&
                  transpose_tr=.TRUE.,&
                  invert_tr=.TRUE., n_rows=p_env%n_ao(spin),&
                  n_cols=p_env%n_mo(spin),error=error)

             ! updates P
             CALL get_mo_set(mos(spin)%mo_set,lfomo=lfomo,&
                  nmo=nmo,maxocc=maxocc)
             IF (lfomo>nmo) THEN
                CALL cp_dbcsr_set(qs_env%rho%rho_ao(spin)%matrix,0.0_dp,error=error)
                CALL cp_dbcsr_plus_fm_fm_t(qs_env%rho%rho_ao(spin)%matrix,&
                     matrix_v=my_psi0(spin)%matrix,&
                     matrix_g=p_env%psi0d(spin)%matrix,&
                     ncol=p_env%n_mo(spin),error=error)
                CALL cp_dbcsr_scale(qs_env%rho%rho_ao(spin)%matrix,alpha_scalar=maxocc,&
                     error=error)
             ELSE
                CALL cp_unimplemented_error(fromWhere=routineP,&
                     message="symmetrized onesided smearing to do",&
                     error=error)
             END IF
          END DO

          ! updates rho
          CALL qs_rho_update_rho(rho_struct=qs_env%rho, &
               qs_env=qs_env, error=error)

          ! tells ks_env that p changed
          CALL qs_ks_did_change(ks_env=ks_env,&
               rho_changed=.TRUE., error=error)

       END IF

       ! updates K (if necessary)
       CALL qs_ks_update_qs_env(ks_env=ks_env,qs_env=qs_env,&
            error=error)
       output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
            extension=".linresLog",error=error)
       IF (output_unit>0) THEN
          CALL section_vals_get(lr_section,explicit=was_present,error=error)
          IF(was_present) THEN
            WRITE(UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))") &
            "Total energy ground state:                     ", energy%total
          END IF
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
            "PRINT%PROGRAM_RUN_INFO",error=error)
       !-----------------------------------------------------------------------|
       ! calculates                                                            |
       ! m_epsilon = - psi0d^T times K times psi0d                             |
       !           = - [K times psi0d]^T times psi0d (because K is symmetric)  |
       !-----------------------------------------------------------------------|
       DO spin=1, n_spins
          ! S_psi0 = k times psi0d
          CALL cp_dbcsr_sm_fm_multiply(matrix_ks(spin)%matrix, &
               p_env%psi0d(spin)%matrix, &
               p_env%S_psi0(spin)%matrix, p_env%n_mo(spin),error=error)
          IF (PRESENT(Hrho_psi0d)) THEN
             CALL cp_fm_scale_and_add(alpha=0.0_dp,matrix_a=Hrho_psi0d(spin)%matrix,&
                  beta=1.0_dp,matrix_b=p_env%S_psi0(spin)%matrix,&
                  error=error)
          END IF
          ! m_epsilon = -1 times S_psi0^T times psi0d
          CALL cp_gemm('T', 'N', &
               p_env%n_mo(spin), p_env%n_mo(spin), p_env%n_ao(spin), &
               -1.0_dp, p_env%S_psi0(spin)%matrix, p_env%psi0d(spin)%matrix, &
               0.0_dp, p_env%m_epsilon(spin)%matrix,error=error)
!   DO i =1,size(p_env%m_epsilon(spin)%matrix%local_data,1)
!    write(*,'(I4,4f12.6)') i,(p_env%m_epsilon(spin)%matrix%local_data(i,j), j=1,4)
!   end do
! stop 'o'
       END DO

       !----------------------------------|
       ! calculates S_psi0 = S * my_psi0  |
       !----------------------------------|
       ! calculating this reduces the mat mult without storing a full aoxao
       ! matrix (for P). If nspin>1 you might consider calculating it on the
       ! fly to spare some memory
       CALL get_qs_env(qs_env, matrix_s=matrix_s,error=error)
       DO spin=1,n_spins
          CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, &
               my_psi0(spin)%matrix, &
               p_env%S_psi0(spin)%matrix, &
               p_env%n_mo(spin),error=error)
       END DO

       ! releases my_psi0
       IF (p_env%orthogonal_orbitals) THEN
          NULLIFY(my_psi0)
       ELSE
          CALL cp_fm_vect_dealloc(my_psi0,error=error)
       END IF

       ! tells kpp1_env about the change of psi0
       CALL kpp1_did_change(p_env%kpp1_env,psi0_changed=.TRUE.,error=error)

    END IF

    CALL timestop(handle)

  END SUBROUTINE p_env_psi0_changed

! *****************************************************************************
!> \brief Evaluates Fv (S_mo)^-1 - Sv(epsilon) and stores it in res
!> \param p_env perturbation calculation environment
!> \param qs_env the qs_env that is perturbed by this p_env
!> \param v the matrix to operate on
!> \param res the result
!> \param error error handling object (optional)
!> \par History
!>      10.2002, TCH, extracted single spin calculation
!> \author Thomas Chassaing
! *****************************************************************************
  SUBROUTINE p_op_l1(p_env, qs_env, v, res, error)

    ! argument
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(in)                             :: v
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(inout)                          :: res
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_op_l1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n_spins, spin
    LOGICAL                                  :: failure

! code

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(p_env%ref_count>0,cp_failure_level,routineP,error,failure)

    n_spins = qs_env%dft_control%nspins
    DO spin=1, n_spins
       CALL p_op_l1_spin(p_env, qs_env, spin, v(spin)%matrix,&
            res(spin)%matrix,error)
    END DO

  END SUBROUTINE p_op_l1

! *****************************************************************************
!> \brief Evaluates Fv (S_mo)^-1 - Sv(epsilon) and stores it in res
!>      for a given spin
!> \param p_env perturbation calculation environment
!> \param qs_env the qs_env that is perturbed by this p_env
!> \param spin the spin to calculate (1 or 2 normally)
!> \param v the matrix to operate on
!> \param res the result
!> \param error error handling object (optional)
!> \par History
!>      10.2002, TCH, created
!> \author Thomas Chassaing
!> \note
!>      Same as p_op_l1 but takes a spin as additional argument.
! *****************************************************************************
  SUBROUTINE p_op_l1_spin(p_env, qs_env, spin, v, res, error)

    ! argument
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: spin
    TYPE(cp_fm_type), POINTER                :: v, res
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_op_l1_spin', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ncol
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: k_p
    TYPE(cp_fm_pool_type), POINTER           :: maxao_maxmo_fm_pool
    TYPE(cp_fm_type), POINTER                :: tmp
    TYPE(cp_para_env_type), POINTER          :: para_env

! code

    CALL timeset(routineN,handle)

    NULLIFY(tmp, matrix_ks, matrix_s)
    failure = .FALSE.

    CPPrecondition(ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(p_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(0<spin,cp_failure_level,routineP,error,failure)
    CPPrecondition(spin<=qs_env%dft_control%nspins,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env, &
                    para_env=para_env,&
                    matrix_s=matrix_s,&
                    matrix_ks=matrix_ks,&
                    error=error)
    CALL mpools_get(qs_env%mpools, maxao_maxmo_fm_pool=maxao_maxmo_fm_pool,&
         error=error)
    CALL fm_pool_create_fm(maxao_maxmo_fm_pool, tmp, error=error)

    k_p => matrix_ks(spin)%matrix
    CALL cp_fm_get_info(v, ncol_global=ncol,error=error)

    IF (p_env%orthogonal_orbitals) THEN
       CALL cp_dbcsr_sm_fm_multiply(k_p, v, res, ncol ,error=error)
    ELSE
       CALL cp_dbcsr_sm_fm_multiply(k_p, v, tmp, ncol ,error=error)
       CALL cp_fm_symm('R', 'U', p_env%n_ao(spin), p_env%n_mo(spin), 1.0_dp, &
         p_env%Smo_inv(spin)%matrix, tmp, 0.0_dp, res,error=error)
    END IF

    CALL cp_fm_symm('R', 'U', p_env%n_ao(spin), p_env%n_mo(spin), 1.0_dp, &
         p_env%m_epsilon(spin)%matrix, v, 0.0_dp, tmp,error=error)
    CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, tmp, &
         res, p_env%n_mo(spin), alpha=1.0_dp,beta=1.0_dp,error=error)
    CALL fm_pool_give_back_fm(maxao_maxmo_fm_pool,tmp,error=error); NULLIFY(tmp)

    CALL timestop(handle)

  END SUBROUTINE p_op_l1_spin

! *****************************************************************************
!> \brief evaluates res = alpha kpp1(v)*psi0 + beta res
!>      with kpp1 evaluated with p=qs_env%rho%rho_ao, p1=p1
!> \param p_env the perturbation environment
!> \param qs_env the qs_env that is perturbed by this p_env
!> \param p1 direction in which evaluate the second derivative
!> \param res place where to store the result
!> \param alpha scale factor of the result (defaults to 1.0)
!> \param beta scale factor of the old values (defaults to 0.0)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.09.2002 adapted for new qs_p_env_type (TC)
!>      03.2003 extended for p1 not taken from v (TC)
!> \author fawzi
!> \note
!>      qs_env%rho must be up to date
!>      it would be better to pass rho1, not p1
! *****************************************************************************
  SUBROUTINE p_op_l2(p_env, qs_env, p1, res, alpha, beta, error)

    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p1
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(INOUT)                          :: res
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: alpha, beta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_op_l2', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: fdiff = .FALSE.

    INTEGER                                  :: handle, ispin, n_spins
    INTEGER, SAVE                            :: iter = 0
    LOGICAL                                  :: failure, gapw, gapw_xc
    REAL(KIND=dp)                            :: my_alpha, my_beta
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)

    failure=.FALSE.

    gapw=qs_env%dft_control%qs_control%gapw
    gapw_xc=qs_env%dft_control%qs_control%gapw_xc
    my_alpha=1.0_dp
    IF (PRESENT(alpha)) my_alpha=alpha
    my_beta=0.0_dp
    IF (PRESENT(beta)) my_beta=beta

    iter=iter+1
    NULLIFY(dft_control, rho)

    CPPrecondition(ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(p_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       CALL p_env_check_i_alloc(p_env,qs_env=qs_env,error=error)

       CALL get_qs_env(qs_env, dft_control=dft_control, rho=rho, error=error)
       n_spins = dft_control%nspins

       DO ispin=1, SIZE(p1)
          ! hack to avoid crashes in ep regs
          IF(.not.ASSOCIATED(p_env%rho1%rho_ao(ispin)%matrix,p1(ispin)%matrix)) THEN
             CALL cp_dbcsr_copy(p_env%rho1%rho_ao(ispin)%matrix,p1(ispin)%matrix,error=error)
          ENDIF
       ENDDO
       CALL qs_rho_update_rho(rho_struct=p_env%rho1, &
                              qs_env=qs_env, error=error)

       IF (fdiff) THEN
          CALL kpp1_calc_k_p_p1_fdiff(kpp1_env=p_env%kpp1_env,qs_env=qs_env,&
                                      k_p_p1=p_env%kpp1, rho=rho ,rho1=p_env%rho1,&
                                      error=error)
       ELSE
          CALL kpp1_calc_k_p_p1(kpp1_env=p_env%kpp1_env, p_env=p_env, qs_env=qs_env,&
                                k_p_p1=p_env%kpp1, rho=rho ,rho1=p_env%rho1, rho1_xc=p_env%rho1,&
                                error=error)
       END IF

       DO ispin=1,n_spins
          CALL cp_dbcsr_sm_fm_multiply(p_env%kpp1(ispin)%matrix,&
                                 p_env%psi0d(ispin)%matrix, res(ispin)%matrix,&
                                 ncol=p_env%n_mo(ispin),&
                                 alpha=my_alpha,beta=my_beta,error=error)
       END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE p_op_l2

! *****************************************************************************
!> \brief does a preorthogonalization of the given matrix:
!>      v = (I-PS)v
!> \param p_env the perturbation environment
!> \param qs_env the qs_env that is perturbed by this p_env
!> \param v matrix to orthogonalize
!> \param n_cols the number of columns of C to multiply (defaults to size(v,2))
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.09.2002 adapted for new qs_p_env_type (TC)
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE p_preortho(p_env, qs_env, v, n_cols, error)

    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(inout)                          :: v
    INTEGER, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: n_cols
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_preortho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cols, handle, max_cols, &
                                                maxnmo, n_spins, nmo2, spin, &
                                                v_cols, v_rows
    LOGICAL                                  :: failure
    TYPE(cp_fm_pool_type), POINTER           :: maxmo_maxmo_fm_pool
    TYPE(cp_fm_struct_type), POINTER         :: maxmo_maxmo_fmstruct, &
                                                tmp_fmstruct
    TYPE(cp_fm_type), POINTER                :: tmp_matrix
    TYPE(dft_control_type), POINTER          :: dft_control

! code

    CALL timeset(routineN,handle)

    failure=.FALSE.

    NULLIFY(tmp_matrix,maxmo_maxmo_fm_pool, maxmo_maxmo_fmstruct,tmp_fmstruct,&
         dft_control)

    CPPrecondition(ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(p_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env, dft_control=dft_control,&
            error=error)
       CALL mpools_get(qs_env%mpools,maxmo_maxmo_fm_pool=maxmo_maxmo_fm_pool,&
            error=error)
       n_spins = qs_env%dft_control%nspins
       maxmo_maxmo_fmstruct => fm_pool_get_el_struct(maxmo_maxmo_fm_pool,error=error)
       CALL cp_fm_struct_get(maxmo_maxmo_fmstruct, nrow_global=nmo2,ncol_global=maxnmo,error=error)
       CPPrecondition(SIZE(v)>=n_spins,cp_failure_level,routineP,error,failure)
    END IF
    IF (.not.failure) THEN
       ! alloc tmp storage
       IF (PRESENT(n_cols)) THEN
          max_cols=MAXVAL(n_cols(1:n_spins))
       ELSE
          max_cols=0
          DO spin=1, n_spins
             CALL cp_fm_get_info(v(spin)%matrix, ncol_global=v_cols,error=error)
             max_cols = MAX(max_cols, v_cols)
          END DO
       END IF
       IF (max_cols <= nmo2) THEN
          CALL fm_pool_create_fm(maxmo_maxmo_fm_pool,tmp_matrix,error=error)
       ELSE
          CALL cp_fm_struct_create(tmp_fmstruct,nrow_global=max_cols,&
               ncol_global=maxnmo,template_fmstruct=maxmo_maxmo_fmstruct,&
               error=error)
          CALL cp_fm_create(tmp_matrix,matrix_struct=tmp_fmstruct,&
               error=error)
          CALL cp_fm_struct_release(tmp_fmstruct,error=error)
       END IF

       DO spin=1, n_spins

          CALL cp_fm_get_info(v(spin)%matrix, &
               nrow_global=v_rows, ncol_global=v_cols,error=error)
          CPPrecondition(v_rows>=p_env%n_ao(spin),cp_failure_level,routineP,error,failure)
          cols = v_cols
          IF (PRESENT(n_cols)) THEN
             CPPrecondition(n_cols(spin)<=cols,cp_failure_level,routineP,error,failure)
             cols=n_cols(spin)
          END IF
          CPPrecondition(cols<=max_cols,cp_failure_level,routineP,error,failure)

          IF (.NOT. failure) THEN
             ! tmp_matrix = v^T (S psi0)
             CALL cp_gemm(transa='T',transb='N',m=cols,n=p_env%n_mo(spin),&
                  k=p_env%n_ao(spin),alpha=1.0_dp,matrix_a=v(spin)%matrix,&
                  matrix_b=p_env%S_psi0(spin)%matrix,beta=0.0_dp,&
                  matrix_c=tmp_matrix,error=error)
             ! v = v- psi0d tmp_matrix^T = v - psi0d psi0^T S v
             CALL cp_gemm(transa='N',transb='T',m=p_env%n_ao(spin),n=cols,&
                  k=p_env%n_mo(spin),alpha=-1.0_dp,&
                  matrix_a=p_env%psi0d(spin)%matrix,matrix_b=tmp_matrix,&
                  beta=1.0_dp,matrix_c=v(spin)%matrix,error=error)
          END IF

       END DO

       IF (max_cols <= nmo2) THEN
          CALL fm_pool_give_back_fm(maxmo_maxmo_fm_pool,tmp_matrix,&
               error=error)
       ELSE
          CALL cp_fm_release(tmp_matrix,error=error)
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE p_preortho

! *****************************************************************************
!> \brief does a postorthogonalization on the given matrix vector:
!>      v = (I-SP) v
!> \param p_env the perturbation environment
!> \param qs_env the qs_env that is perturbed by this p_env
!> \param v matrix to orthogonalize
!> \param n_cols the number of columns of C to multiply (defaults to size(v,2))
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE p_postortho(p_env, qs_env, v, n_cols, error)

    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(inout)                          :: v
    INTEGER, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: n_cols
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_postortho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cols, handle, max_cols, &
                                                maxnmo, n_spins, nmo2, spin, &
                                                v_cols, v_rows
    LOGICAL                                  :: failure
    TYPE(cp_fm_pool_type), POINTER           :: maxmo_maxmo_fm_pool
    TYPE(cp_fm_struct_type), POINTER         :: maxmo_maxmo_fmstruct, &
                                                tmp_fmstruct
    TYPE(cp_fm_type), POINTER                :: tmp_matrix
    TYPE(dft_control_type), POINTER          :: dft_control

! code

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(tmp_matrix,maxmo_maxmo_fm_pool, maxmo_maxmo_fmstruct,tmp_fmstruct,&
         dft_control)

    CPPrecondition(ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(p_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env, dft_control=dft_control,error=error)
       CALL mpools_get(qs_env%mpools, maxmo_maxmo_fm_pool=maxmo_maxmo_fm_pool,&
            error=error)
       n_spins = dft_control%nspins
       maxmo_maxmo_fmstruct => fm_pool_get_el_struct(maxmo_maxmo_fm_pool,error=error)
       CALL cp_fm_struct_get(maxmo_maxmo_fmstruct, nrow_global=nmo2,ncol_global=maxnmo,error=error)
       CPPrecondition(SIZE(v)>=n_spins,cp_failure_level,routineP,error,failure)
    END IF
    IF (.not.failure) THEN
       ! alloc tmp storage
       IF (PRESENT(n_cols)) THEN
          max_cols=MAXVAL(n_cols(1:n_spins))
       ELSE
          max_cols=0
          DO spin=1, n_spins
             CALL cp_fm_get_info(v(spin)%matrix, ncol_global=v_cols,error=error)
             max_cols = MAX(max_cols, v_cols)
          END DO
       END IF
       IF (max_cols <= nmo2) THEN
          CALL fm_pool_create_fm(maxmo_maxmo_fm_pool,tmp_matrix,error=error)
       ELSE
          CALL cp_fm_struct_create(tmp_fmstruct,nrow_global=max_cols,&
               ncol_global=maxnmo,template_fmstruct=maxmo_maxmo_fmstruct,&
               error=error)
          CALL cp_fm_create(tmp_matrix,matrix_struct=tmp_fmstruct,&
               error=error)
          CALL cp_fm_struct_release(tmp_fmstruct,error=error)
       END IF

       DO spin=1, n_spins

          CALL cp_fm_get_info(v(spin)%matrix, &
               nrow_global=v_rows, ncol_global=v_cols,error=error)
          CPPrecondition(v_rows>=p_env%n_ao(spin),cp_failure_level,routineP,error,failure)
          cols = v_cols
          IF (PRESENT(n_cols)) THEN
             CPPrecondition(n_cols(spin)<=cols,cp_failure_level,routineP,error,failure)
             cols=n_cols(spin)
          END IF
          CPPrecondition(cols<=max_cols,cp_failure_level,routineP,error,failure)

          IF (.NOT. failure) THEN
             ! tmp_matrix = v^T psi0d
             CALL cp_gemm(transa='T',transb='N',m=cols,n=p_env%n_mo(spin),&
                  k=p_env%n_ao(spin),alpha=1.0_dp,matrix_a=v(spin)%matrix,&
                  matrix_b=p_env%psi0d(spin)%matrix,beta=0.0_dp,&
                  matrix_c=tmp_matrix,error=error)
             ! v = v- (S psi0) tmp_matrix^T = v - S psi0 psi0d^T v
             CALL cp_gemm(transa='N',transb='T',m=p_env%n_ao(spin),n=cols,&
                  k=p_env%n_mo(spin),alpha=-1.0_dp,&
                  matrix_a=p_env%S_psi0(spin)%matrix,matrix_b=tmp_matrix,&
                  beta=1.0_dp,matrix_c=v(spin)%matrix,error=error)
          END IF

       END DO

       IF (max_cols <= nmo2) THEN
          CALL fm_pool_give_back_fm(maxmo_maxmo_fm_pool,tmp_matrix,&
               error=error)
       ELSE
          CALL cp_fm_release(tmp_matrix,error=error)
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE p_postortho

! *****************************************************************************
!> \brief evaluates res = alpha kpp1(v)*psi0 + beta res
!> \param p_env the perturbation environment
!> \param qs_env the qs_env that is perturbed by this p_env
!> \param v values of the coeffs of the correcting orbitals
!> \param res place where to store the result
!> \param alpha scale factor of the result (defaults to 1.0)
!> \param beta scale factor of the old values (defaults to 0.0)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.09.2002 adapted for new qs_p_env_type (TC)
!> \author Fawzi Mohamed
!> \note
!>      to be removed
! *****************************************************************************
  SUBROUTINE p_op_l2_fawzi(p_env, qs_env, v, res, alpha, beta, error)

    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(in)                             :: v
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(inout)                          :: res
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: alpha, beta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_op_l2_fawzi', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, lfomo, &
                                                n_spins, nmo
    INTEGER, SAVE                            :: iter = 0
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: maxocc, my_alpha, my_beta
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, p1_new
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)

    failure=.FALSE.

    my_alpha=1.0_dp
    IF (PRESENT(alpha)) my_alpha=alpha
    my_beta=0.0_dp
    IF (PRESENT(beta)) my_beta=beta

    iter=iter+1
    NULLIFY(p1_new, matrix_s, dft_control, para_env, rho, mos)

    CPPrecondition(ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(p_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       CALL p_env_check_i_alloc(p_env,qs_env=qs_env,error=error)

       CALL get_qs_env(qs_env,dft_control=dft_control,&
            para_env=para_env, rho=rho, matrix_s=matrix_s,mos=mos,&
            error=error)
       n_spins=dft_control%nspins

       p1_new => p_env%rho1%rho_ao

    END IF
    IF (.not.failure) THEN
       DO ispin=1,n_spins
          CALL get_mo_set(mos(ispin)%mo_set,lfomo=lfomo,nmo=nmo,&
               maxocc=maxocc)
          IF (lfomo>nmo) THEN
             CALL cp_dbcsr_set(p1_new(ispin)%matrix,0.0_dp,error=error)
             ! optimize?
             CALL cp_dbcsr_plus_fm_fm_t(p1_new(ispin)%matrix,&
                  matrix_v=v(ispin)%matrix,&
                  matrix_g=p_env%psi0d(ispin)%matrix,&
                  ncol=p_env%n_mo(ispin),error=error)
             CALL cp_dbcsr_plus_fm_fm_t(p1_new(ispin)%matrix,&
                  matrix_v=p_env%psi0d(ispin)%matrix,&
                  matrix_g=v(ispin)%matrix,&
                  ncol=p_env%n_mo(ispin),error=error)
             CALL cp_dbcsr_scale(p1_new(ispin)%matrix,alpha_scalar=maxocc,error=error)
          ELSE
             CALL cp_unimplemented_error(fromWhere=routineP,&
                  message="symmetrized onesided smearing to do",&
                  error=error)
          END IF
          IF (debug_this_module) THEN
!FM             CALL cp_sm_output(logger=logger, &
!FM                  outputName="p_env"//TRIM(ADJUSTL(cp_to_string(p_env%id_nr)))//&
!FM                  "%p1_new-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
!FM                  fromWhere=routineP, iter=iter, matrix=p1_new(ispin)%matrix,&
!FM                  para_env=qs_env%para_env, error=error)
          END IF
       END DO
       CALL p_op_l2(p_env, qs_env, p1_new, res, alpha=my_alpha, beta=my_beta, error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE p_op_l2_fawzi

END MODULE qs_p_env_methods
