!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for a linear scaling quickstep SCF run based on the density
!>        matrix
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE dm_ls_scf
  USE bibliography,                    ONLY: Kolafa2004,&
                                             Lin2009,&
                                             Lin2013,&
                                             Niklasson2003,&
                                             Niklasson2014,&
                                             Shao2003,&
                                             VandeVondele2012,&
                                             cite_reference
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_add_on_diag, cp_dbcsr_arnoldi_extremal, &
       cp_dbcsr_binary_read, cp_dbcsr_binary_write, cp_dbcsr_checksum, &
       cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_distribution, &
       cp_dbcsr_filter, cp_dbcsr_frobenius_norm, cp_dbcsr_get_occupation, &
       cp_dbcsr_init, cp_dbcsr_multiply, cp_dbcsr_p_type, cp_dbcsr_release, &
       cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_type, dbcsr_type_no_symmetry
  USE cp_external_control,             ONLY: external_control
  USE cp_para_env,                     ONLY: cp_para_env_retain
  USE dm_ls_chebyshev,                 ONLY: compute_chebyshev
  USE dm_ls_scf_curvy,                 ONLY: deallocate_curvy_data,&
                                             dm_ls_curvy_optimization
  USE dm_ls_scf_methods,               ONLY: apply_matrix_preconditioner,&
                                             compute_homo_lumo,&
                                             compute_matrix_preconditioner,&
                                             density_matrix_sign,&
                                             density_matrix_sign_fixed_mu,&
                                             density_matrix_tc2,&
                                             density_matrix_trs4
  USE dm_ls_scf_qs,                    ONLY: ls_scf_dm_to_ks,&
                                             ls_scf_init_qs,&
                                             ls_scf_qs_atomic_guess,&
                                             matrix_ls_create,&
                                             matrix_ls_to_qs,&
                                             matrix_qs_to_ls,&
                                             rho_mixing_ls_init
  USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
  USE input_constants,                 ONLY: &
       ls_cluster_atomic, ls_cluster_molecular, ls_s_inversion_hotelling, &
       ls_s_inversion_none, ls_s_inversion_sign_sqrt, &
       ls_s_preconditioner_atomic, ls_s_preconditioner_molecular, &
       ls_s_preconditioner_none, ls_scf_ns, ls_scf_pexsi, ls_scf_tc2, &
       ls_scf_trs4
  USE input_enumeration_types,         ONLY: enum_i2c,&
                                             enumeration_type
  USE input_keyword_types,             ONLY: keyword_get,&
                                             keyword_type
  USE input_section_types,             ONLY: section_get_keyword,&
                                             section_release,&
                                             section_type,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_retain,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE iterate_matrix,                  ONLY: invert_Hotelling,&
                                             matrix_sqrt_newton_schulz,&
                                             purify_mcweeny
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE mathlib,                         ONLY: binomial
  USE molecule_types_new,              ONLY: molecule_of_atom,&
                                             molecule_type
  USE particle_types,                  ONLY: particle_type
  USE pexsi_methods,                   ONLY: density_matrix_pexsi,&
                                             pexsi_finalize_scf,&
                                             pexsi_init_read_input,&
                                             pexsi_init_scf,&
                                             pexsi_set_convergence_tolerance,&
                                             pexsi_to_qs
  USE pexsi_types,                     ONLY: lib_pexsi_init
  USE qs_density_mixing_types,         ONLY: create_mixing_section,&
                                             mixing_storage_create
  USE qs_diis,                         ONLY: qs_diis_b_clear_sparse,&
                                             qs_diis_b_create_sparse,&
                                             qs_diis_b_step_4lscf
  USE qs_diis_types,                   ONLY: qs_diis_b_release_sparse,&
                                             qs_diis_buffer_type_sparse
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_scf_post_gpw,                 ONLY: qs_scf_post_moments,&
                                             write_mo_free_results
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE transport,                       ONLY: external_scf_method,&
                                             transport_initialize
  USE transport_env_types,             ONLY: transport_env_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dm_ls_scf'

  PUBLIC :: ls_scf,ls_scf_create

CONTAINS

! *****************************************************************************
!> \brief perform an linear scaling scf procedure: entry point
!>
!> \param qs_env ...
!> \param error ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(ls_scf_env_type), POINTER           :: ls_scf_env

    CALL timeset(routineN,handle)

    ! get scf env
    CALL get_qs_env(qs_env,ls_scf_env=ls_scf_env,error=error)

    ! initialize the scf phase
    CALL ls_scf_init_scf(qs_env,ls_scf_env,error)

    ! perform the actual scf
    CALL ls_scf_main(qs_env,ls_scf_env,error)

    ! do post scf processing
    CALL ls_scf_post(qs_env,ls_scf_env,error)

    CALL timestop(handle)

  END SUBROUTINE ls_scf

! *****************************************************************************
!> \brief Creation and basic initialization of the LS type.
!> \param qs_env ...
!> \param error ...
!> \par History
!>       2012.11 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_create(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat, unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(ls_scf_env_type), POINTER           :: ls_scf_env
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: input, mixing_section

    CALL timeset(routineN,handle)

    failure=.FALSE.

    CALL cite_reference(VandeVondele2012)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ALLOCATE(ls_scf_env,stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)


    ! get basic quantities from the qs_env
    CALL get_qs_env(qs_env,nelectron_total=ls_scf_env%nelectron_total, &
                           matrix_s=matrix_s,&
                           dft_control=dft_control,&
                           particle_set=particle_set,&
                           molecule_set=molecule_set,&
                           input=input,&
                           has_unit_metric=ls_scf_env%has_unit_metric,&
                           para_env=ls_scf_env%para_env,&
                           do_transport=ls_scf_env%do_transport,&
                           nelectron_spin=ls_scf_env%nelectron_spin,error=error)

    ! copy some basic stuff
    ls_scf_env%nspins=dft_control%nspins
    ls_scf_env%natoms=SIZE(particle_set,1)
    CALL cp_para_env_retain(ls_scf_env%para_env,error)

    ! initialize block to group to defined molecules
    ALLOCATE(ls_scf_env%ls_mstruct%atom_to_molecule(ls_scf_env%natoms),stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL molecule_of_atom(molecule_set,atom_to_mol=ls_scf_env%ls_mstruct%atom_to_molecule,error=error)

    ! parse the ls_scf section and set derived quantities
    CALL ls_scf_init_read_write_input(input,ls_scf_env,unit_nr,error)

    ! set up the buffer for the history of matrices
    ls_scf_env%scf_history%nstore=ls_scf_env%extrapolation_order
    ls_scf_env%scf_history%istore=0
    ALLOCATE(ls_scf_env%scf_history%matrix(ls_scf_env%nspins,ls_scf_env%scf_history%nstore))

    NULLIFY(ls_scf_env%mixing_store)
    mixing_section => section_vals_get_subs_vals(input,"DFT%LS_SCF%RHO_MIXING",error=error)
    CALL mixing_storage_create(ls_scf_env%mixing_store, mixing_section, &
                               ls_scf_env%density_mixing_method, &
                               dft_control%qs_control%cutoff, error)

    ! initialize PEXSI
    IF (ls_scf_env%purification_method .EQ. ls_scf_pexsi) THEN
      CALL cp_assert(dft_control%qs_control%eps_filter_matrix .EQ. 0.0_dp, &
                     cp_failure_level, cp_assertion_failed, routineP, &
                     "EPS_FILTER_MATRIX must be set to 0 for PEXSI.", &
                     only_ionode=.TRUE., error = error)
       CALL lib_pexsi_init(ls_scf_env%pexsi, ls_scf_env%para_env%group, ls_scf_env%nspins, error)
    ENDIF
 
    ! put the ls_scf_env in qs_env
    CALL set_qs_env(qs_env,ls_scf_env=ls_scf_env,error=error)

    CALL timestop(handle)

  END SUBROUTINE ls_scf_create

! *****************************************************************************
!> \brief initialization needed for scf
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param error ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_init_scf(qs_env,ls_scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_init_scf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin, unit_nr
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! get basic quantities from the qs_env
    CALL get_qs_env(qs_env,nelectron_total=ls_scf_env%nelectron_total, &
                           matrix_s=matrix_s,&
                           ks_env=ks_env,&
                           dft_control=dft_control,&
                           molecule_set=molecule_set,&
                           input=input,&
                           has_unit_metric=ls_scf_env%has_unit_metric,&
                           para_env=ls_scf_env%para_env,&
                           nelectron_spin=ls_scf_env%nelectron_spin,error=error)

    ! some basic initialization of the QS side of things
    CALL ls_scf_init_qs(qs_env,error)

    ! create the matrix template for use in the ls procedures
    CALL matrix_ls_create(matrix_ls=ls_scf_env%matrix_s,matrix_qs=matrix_s(1)%matrix,&
                          ls_mstruct=ls_scf_env%ls_mstruct,error=error)

    nspin=ls_scf_env%nspins
    ALLOCATE(ls_scf_env%matrix_p(nspin))
    DO ispin=1,nspin
       CALL cp_dbcsr_init(ls_scf_env%matrix_p(ispin),error=error)
       CALL cp_dbcsr_create(ls_scf_env%matrix_p(ispin),template=ls_scf_env%matrix_s,&
                            matrix_type=dbcsr_type_no_symmetry, error=error)
    ENDDO

    ALLOCATE(ls_scf_env%matrix_ks(nspin))
    DO ispin=1,nspin
       CALL cp_dbcsr_init(ls_scf_env%matrix_ks(ispin),error=error)
       CALL cp_dbcsr_create(ls_scf_env%matrix_ks(ispin),template=ls_scf_env%matrix_s,&
                            matrix_type=dbcsr_type_no_symmetry, error=error)
    ENDDO

    ! set up matrix S, and needed functions of S
    CALL ls_scf_init_matrix_s(matrix_s(1)%matrix,ls_scf_env,error)

    ! get the initial guess for the SCF
    CALL ls_scf_initial_guess(qs_env,ls_scf_env,error)

    IF (ls_scf_env%do_rho_mixing) THEN
      CALL rho_mixing_ls_init(qs_env, ls_scf_env, error)
    ENDIF

    IF (ls_scf_env%purification_method .EQ. ls_scf_pexsi) THEN
       CALL pexsi_init_scf(ks_env, ls_scf_env%pexsi, matrix_s(1)%matrix, error)
    ENDIF

    IF (qs_env%do_transport) THEN
       CALL transport_initialize(ks_env,qs_env%transport_env, matrix_s(1)%matrix, error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE ls_scf_init_scf

! *****************************************************************************
!> \brief deal with the scf initial guess
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param error ...
!> \par History
!>       2012.11 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_initial_guess(qs_env,ls_scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_initial_guess', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: aspc_guess = 2, &
                                                atomic_guess = 1, &
                                                restart_guess = 3

    CHARACTER(LEN=default_path_length)       :: file_name, project_name
    INTEGER                                  :: handle, iaspc, &
                                                initial_guess_type, ispin, &
                                                istore, naspc, unit_nr
    REAL(KIND=dp)                            :: alpha, cs_pos
    TYPE(cp_dbcsr_type)                      :: matrix_tmp1
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF


    IF (unit_nr>0) WRITE(unit_nr,'()')
    ! if there is no history go for the atomic guess, otherwise extrapolate the dm history
    IF (ls_scf_env%scf_history%istore==0) THEN
       IF (ls_scf_env%restart_read) THEN
         initial_guess_type=restart_guess
       ELSE
         initial_guess_type=atomic_guess
       ENDIF
    ELSE
       initial_guess_type=aspc_guess
    ENDIF

    ! how to get the initial guess
    SELECT CASE(initial_guess_type)
    CASE(atomic_guess)
       CALL ls_scf_qs_atomic_guess(qs_env,ls_scf_env%energy_init,error)
    CASE(restart_guess)
       project_name = logger%iter_info%project_name
       DO ispin=1,SIZE(ls_scf_env%matrix_p)
          WRITE(file_name,'(A,I0,A)') TRIM(project_name)//"_LS_DM_SPIN_",ispin,"_RESTART.dm"
          CALL cp_dbcsr_binary_read(file_name, distribution=cp_dbcsr_distribution(ls_scf_env%matrix_p(1)), &
                                    matrix_new=ls_scf_env%matrix_p(ispin), error=error)
          cs_pos = cp_dbcsr_checksum (ls_scf_env%matrix_p(ispin), pos=.TRUE., error=error)
          IF (unit_nr>0) THEN
             WRITE(unit_nr,'(T2,A,E20.8)') "Read restart DM "//TRIM(file_name)//" with checksum: ",cs_pos
          ENDIF
       ENDDO

       ! directly go to computing the corresponding energy and ks matrix
       CALL ls_scf_dm_to_ks(qs_env,ls_scf_env,ls_scf_env%energy_init,iscf=0,error=error)
    CASE(aspc_guess)
       CALL cite_reference(Kolafa2004)
       naspc=MIN(ls_scf_env%scf_history%istore,ls_scf_env%scf_history%nstore)
       DO ispin=1,SIZE(ls_scf_env%matrix_p)
          ! actual extrapolation
          CALL cp_dbcsr_set(ls_scf_env%matrix_p(ispin),0.0_dp,error=error)
          DO iaspc=1,naspc
             alpha=(-1.0_dp)**(iaspc + 1)*REAL(iaspc,KIND=dp)*&
                  binomial(2*naspc,naspc - iaspc)/binomial(2*naspc - 2,naspc -1)
             istore=MOD(ls_scf_env%scf_history%istore-iaspc,ls_scf_env%scf_history%nstore)+1
             CALL cp_dbcsr_add(ls_scf_env%matrix_p(ispin), ls_scf_env%scf_history%matrix(ispin,istore), 1.0_dp, alpha, error=error)
          ENDDO
       ENDDO
    END SELECT

    ! which cases need getting purified and non-orthogonal ?
    SELECT CASE(initial_guess_type)
    CASE(atomic_guess, restart_guess)
      ! do nothing
    CASE(aspc_guess)
       ! purification can't be done on the pexsi matrix, which is not necessarily idempotent,
       ! and not stored in an ortho basis form
       IF (.NOT.(ls_scf_env%purification_method .EQ. ls_scf_pexsi)) THEN
         DO ispin=1,SIZE(ls_scf_env%matrix_p)
            ! linear combination of P's is not idempotent. A bit of McWeeny is needed to ensure it is again
            IF(SIZE(ls_scf_env%matrix_p)==1)CALL cp_dbcsr_scale(ls_scf_env%matrix_p(ispin),0.5_dp,error=error)
            CALL purify_mcweeny(ls_scf_env%matrix_p(ispin:ispin),ls_scf_env%eps_filter,3,error)
            IF(SIZE(ls_scf_env%matrix_p)==1)CALL cp_dbcsr_scale(ls_scf_env%matrix_p(ispin),2.0_dp,error=error)

            IF (ls_scf_env%use_s_sqrt) THEN
               ! need to get P in the non-orthogonal basis if it was stored differently
               CALL cp_dbcsr_init(matrix_tmp1,error=error)
               CALL cp_dbcsr_create(matrix_tmp1,template=ls_scf_env%matrix_s,&
                                    matrix_type=dbcsr_type_no_symmetry,error=error)
               CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt_inv, ls_scf_env%matrix_p(ispin),&
                                            0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)
               CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s_sqrt_inv, &
                                            0.0_dp, ls_scf_env%matrix_p(ispin)  , &
                                            filter_eps=ls_scf_env%eps_filter,error=error)
               CALL cp_dbcsr_release(matrix_tmp1,error=error)

               IF (ls_scf_env%has_s_preconditioner) THEN
                   CALL apply_matrix_preconditioner(ls_scf_env%matrix_p(ispin),"forward", &
                                  ls_scf_env%matrix_bs_sqrt,ls_scf_env%matrix_bs_sqrt_inv,error)
               ENDIF
            ENDIF
         ENDDO
       ENDIF

       ! compute corresponding energy and ks matrix
       CALL ls_scf_dm_to_ks(qs_env,ls_scf_env,ls_scf_env%energy_init,iscf=0,error=error)
    END SELECT

    IF (unit_nr>0) THEN
       WRITE(unit_nr,'(T2,A,F20.9)') "Energy with the initial guess:",ls_scf_env%energy_init
       WRITE(unit_nr,'()')
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE ls_scf_initial_guess

! *****************************************************************************
!> \brief store a history of matrices for later use in ls_scf_initial_guess
!> \param ls_scf_env ...
!> \param error ...
!> \par History
!>       2012.11 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_store_result(ls_scf_env,error)
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_store_result', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name, project_name
    INTEGER                                  :: handle, ispin, istore, &
                                                nelectron_spin_real, unit_nr
    REAL(KIND=dp)                            :: cs_pos, homo_spin, lumo_spin, &
                                                mu_spin
    TYPE(cp_dbcsr_type)                      :: matrix_ks_deviation_tmp, &
                                                matrix_p_tmp, matrix_tmp1
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    IF (ls_scf_env%restart_write) THEN
       DO ispin=1,SIZE(ls_scf_env%matrix_p)
          project_name = logger%iter_info%project_name
          WRITE(file_name,'(A,I0,A)') TRIM(project_name)//"_LS_DM_SPIN_",ispin,"_RESTART.dm"
          cs_pos = cp_dbcsr_checksum (ls_scf_env%matrix_p(ispin), pos=.TRUE., error=error)
          IF (unit_nr>0) THEN
             WRITE(unit_nr,'(T2,A,E20.8)') "Writing restart DM "//TRIM(file_name)//" with checksum: ",cs_pos
          ENDIF

          ! perform trs4 once to get a full dm required by other purification methods
          IF (ls_scf_env%do_transport .OR. ls_scf_env%purification_method .EQ. ls_scf_pexsi) THEN
             IF (unit_nr>0) THEN
                WRITE(unit_nr,'(T2,A)') "DM has the sparsity of S! "
                WRITE(unit_nr,'(T2,A)') "TRS4 will be performed once to prepare "//&
                                        "a full DM to store in the restart file."
             ENDIF

             CALL cp_dbcsr_init(matrix_p_tmp,error=error)
             CALL cp_dbcsr_create(matrix_p_tmp,template=ls_scf_env%matrix_ks(ispin),error=error)
             CALL cp_dbcsr_set(matrix_p_tmp,0.0_dp,error=error)

             CALL cp_dbcsr_init(matrix_ks_deviation_tmp,error=error)
             CALL cp_dbcsr_create(matrix_ks_deviation_tmp,template=ls_scf_env%matrix_ks(ispin),error=error)
             CALL cp_dbcsr_set(matrix_ks_deviation_tmp,0.0_dp,error=error)

             nelectron_spin_real=ls_scf_env%nelectron_spin(ispin)
             IF (ls_scf_env%nspins==1) nelectron_spin_real=nelectron_spin_real/2
             CALL density_matrix_trs4(matrix_p_tmp, ls_scf_env%matrix_ks(ispin),  ls_scf_env%matrix_s_sqrt_inv,&
                                      nelectron_spin_real, 1.0E-06_dp, homo_spin, lumo_spin, mu_spin, &
                                      dynamic_threshold=.FALSE., matrix_ks_deviation=matrix_ks_deviation_tmp, &
                                      eps_lanczos=1.0E-03_dp, max_iter_lanczos=128, error=error)
             IF (ls_scf_env%nspins==1) CALL cp_dbcsr_scale(matrix_p_tmp, 2.0_dp, error=error)

             CALL cp_dbcsr_binary_write(matrix_p_tmp,file_name,error)

             CALL cp_dbcsr_release(matrix_p_tmp,error=error)
             CALL cp_dbcsr_release(matrix_ks_deviation_tmp,error=error)
          ELSE
             CALL cp_dbcsr_binary_write(ls_scf_env%matrix_p(ispin),file_name,error)
          END IF 
       ENDDO
    END IF

    IF (ls_scf_env%scf_history%nstore>0) THEN
       ls_scf_env%scf_history%istore=ls_scf_env%scf_history%istore+1
       DO ispin=1,SIZE(ls_scf_env%matrix_p)
          istore=MOD(ls_scf_env%scf_history%istore-1,ls_scf_env%scf_history%nstore)+1
          IF (ls_scf_env%scf_history%istore<=ls_scf_env%scf_history%nstore) &
             CALL cp_dbcsr_init(ls_scf_env%scf_history%matrix(ispin,istore),error=error)
          CALL cp_dbcsr_copy(ls_scf_env%scf_history%matrix(ispin,istore), ls_scf_env%matrix_p(ispin), error=error) 

          ! if we have the sqrt around, we use it to go to the orthogonal basis
          IF (ls_scf_env%use_s_sqrt) THEN
             ! usualy sqrt(S) * P * sqrt(S) should be available, or could be stored at least,
             ! so that the next multiplications could be saved.
             CALL cp_dbcsr_init(matrix_tmp1,error=error)
             CALL cp_dbcsr_create(matrix_tmp1,template=ls_scf_env%matrix_s,&
                                  matrix_type=dbcsr_type_no_symmetry,error=error)

             IF (ls_scf_env%has_s_preconditioner) THEN
                 CALL apply_matrix_preconditioner(ls_scf_env%scf_history%matrix(ispin,istore),"backward", &
                                ls_scf_env%matrix_bs_sqrt,ls_scf_env%matrix_bs_sqrt_inv,error)
             ENDIF
             CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt, ls_scf_env%scf_history%matrix(ispin,istore),&
                                          0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)
             CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s_sqrt, &
                                          0.0_dp, ls_scf_env%scf_history%matrix(ispin,istore) , &
                                          filter_eps=ls_scf_env%eps_filter,error=error)
             CALL cp_dbcsr_release(matrix_tmp1,error=error)
          ENDIF

       ENDDO
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE ls_scf_store_result

! *****************************************************************************
!> \brief initialize S matrix related properties (sqrt, inverse...)
!>        Might be factored-out since this seems common code with the other SCF.
!> \param matrix_s ...
!> \param ls_scf_env ...
!> \param error ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_init_matrix_S(matrix_s,ls_scf_env,error)
    TYPE(cp_dbcsr_type)                      :: matrix_s
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_init_matrix_S', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, unit_nr
    REAL(KIND=dp)                            :: frob_matrix, frob_matrix_base
    TYPE(cp_dbcsr_type)                      :: matrix_tmp1, matrix_tmp2
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! make our own copy of S
    IF (ls_scf_env%has_unit_metric) THEN
       CALL cp_dbcsr_set(ls_scf_env%matrix_s,0.0_dp,error=error)
       CALL cp_dbcsr_add_on_diag(ls_scf_env%matrix_s,1.0_dp,error=error)
    ELSE
       CALL matrix_qs_to_ls(ls_scf_env%matrix_s,matrix_s,ls_scf_env%ls_mstruct,error=error)
    ENDIF

    CALL cp_dbcsr_filter(ls_scf_env%matrix_s,ls_scf_env%eps_filter,error=error)

    ! needs a preconditioner for S
    IF (ls_scf_env%has_s_preconditioner) THEN
       CALL cp_dbcsr_init(ls_scf_env%matrix_bs_sqrt,error=error)
       CALL cp_dbcsr_create(ls_scf_env%matrix_bs_sqrt,template=ls_scf_env%matrix_s,&
                            matrix_type=dbcsr_type_no_symmetry, error=error)
       CALL cp_dbcsr_init(ls_scf_env%matrix_bs_sqrt_inv,error=error)
       CALL cp_dbcsr_create(ls_scf_env%matrix_bs_sqrt_inv,template=ls_scf_env%matrix_s,&
                            matrix_type=dbcsr_type_no_symmetry, error=error)
       CALL compute_matrix_preconditioner(ls_scf_env%matrix_s,&
                            ls_scf_env%s_preconditioner_type, ls_scf_env%ls_mstruct, &
                            ls_scf_env%matrix_bs_sqrt,ls_scf_env%matrix_bs_sqrt_inv,&
                            ls_scf_env%eps_filter,ls_scf_env%sign_sqrt_order,&
                            ls_scf_env%eps_lanczos,ls_scf_env%max_iter_lanczos,error)
    ENDIF

    ! precondition S
    IF (ls_scf_env%has_s_preconditioner) THEN
       CALL apply_matrix_preconditioner(ls_scf_env%matrix_s,"forward", &
                            ls_scf_env%matrix_bs_sqrt,ls_scf_env%matrix_bs_sqrt_inv,error)
    ENDIF

    ! compute sqrt(S) and inv(sqrt(S))
    IF (ls_scf_env%use_s_sqrt) THEN
        CALL cp_dbcsr_init(ls_scf_env%matrix_s_sqrt,error=error)
        CALL cp_dbcsr_init(ls_scf_env%matrix_s_sqrt_inv,error=error)
        CALL cp_dbcsr_create(ls_scf_env%matrix_s_sqrt,template=ls_scf_env%matrix_s,&
                             matrix_type=dbcsr_type_no_symmetry,error=error)
        CALL cp_dbcsr_create(ls_scf_env%matrix_s_sqrt_inv,template=ls_scf_env%matrix_s,&
                             matrix_type=dbcsr_type_no_symmetry,error=error)


        CALL matrix_sqrt_Newton_Schulz(ls_scf_env%matrix_s_sqrt,ls_scf_env%matrix_s_sqrt_inv,&
                                       ls_scf_env%matrix_s,ls_scf_env%eps_filter,&
                                       ls_scf_env%sign_sqrt_order, &
                                       ls_scf_env%eps_lanczos, ls_scf_env%max_iter_lanczos, error)

        IF (.TRUE.) THEN
           CALL cp_dbcsr_init(matrix_tmp1,error=error)
           CALL cp_dbcsr_create(matrix_tmp1,template=ls_scf_env%matrix_s,&
                                matrix_type=dbcsr_type_no_symmetry,error=error)
           CALL cp_dbcsr_init(matrix_tmp2,error=error)
           CALL cp_dbcsr_create(matrix_tmp2,template=ls_scf_env%matrix_s,&
                                matrix_type=dbcsr_type_no_symmetry,error=error)

           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt_inv, ls_scf_env%matrix_s,&
                                        0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)

           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s_sqrt_inv, &
                                        0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter,error=error)

           frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp2)
           CALL cp_dbcsr_add_on_diag(matrix_tmp2,-1.0_dp,error=error)
           frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp2)
           IF (unit_nr>0) THEN
              WRITE(unit_nr,*) "Error for (inv(sqrt(S))*S*inv(sqrt(S))-I)",frob_matrix/frob_matrix_base
           ENDIF

           CALL cp_dbcsr_release(matrix_tmp1,error=error)
           CALL cp_dbcsr_release(matrix_tmp2,error=error)
        ENDIF
    ENDIF

    ! compute the inverse of S
    IF (ls_scf_env%needs_s_inv) THEN
        CALL cp_dbcsr_init(ls_scf_env%matrix_s_inv,error=error)
        CALL cp_dbcsr_create(ls_scf_env%matrix_s_inv,template=ls_scf_env%matrix_s,&
                             matrix_type=dbcsr_type_no_symmetry,error=error)
        IF (.NOT.ls_scf_env%use_s_sqrt) THEN
           CALL invert_Hotelling(ls_scf_env%matrix_s_inv,ls_scf_env%matrix_s,ls_scf_env%eps_filter,error=error)
        ELSE
           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt_inv, ls_scf_env%matrix_s_sqrt_inv, &
                                        0.0_dp, ls_scf_env%matrix_s_inv, filter_eps=ls_scf_env%eps_filter,error=error)
        ENDIF
        IF (.TRUE.) THEN
           CALL cp_dbcsr_init(matrix_tmp1,error=error)
           CALL cp_dbcsr_create(matrix_tmp1,template=ls_scf_env%matrix_s,&
                                matrix_type=dbcsr_type_no_symmetry,error=error)
           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_inv, ls_scf_env%matrix_s,&
                                        0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)
           frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp1)
           CALL cp_dbcsr_add_on_diag(matrix_tmp1,-1.0_dp,error=error)
           frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp1)
           IF (unit_nr>0) THEN
              WRITE(unit_nr,*) "Error for (inv(S)*S-I)",frob_matrix/frob_matrix_base
           ENDIF
           CALL cp_dbcsr_release(matrix_tmp1,error=error)
        ENDIF
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE ls_scf_init_matrix_s

! *****************************************************************************
!> \brief parse the input section, no need to pass it around
!> \param input ...
!> \param ls_scf_env ...
!> \param unit_nr ...
!> \param error ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_init_read_write_input(input,ls_scf_env,unit_nr,error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(ls_scf_env_type), INTENT(INOUT)     :: ls_scf_env
    INTEGER, INTENT(IN)                      :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_init_read_write_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: mu
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER         :: chebyshev_section, &
                                                curvy_section, &
                                                ls_scf_section, &
                                                mixing_section, pexsi_section

    CALL timeset(routineN,handle)
    CALL cite_reference(VandeVondele2012)
    failure=.FALSE.
    ls_scf_section => section_vals_get_subs_vals(input,"DFT%LS_SCF",error=error)
    curvy_section => section_vals_get_subs_vals(ls_scf_section,"CURVY_STEPS",error=error)

    ! should come from input
    CALL section_vals_val_get(ls_scf_section,"LS_DIIS",l_val=ls_scf_env%ls_diis,error=error)
    CALL section_vals_val_get(ls_scf_section,"INI_DIIS",i_val=ls_scf_env%iter_ini_diis,error=error)
    CALL section_vals_val_get(ls_scf_section,"MAX_DIIS",i_val=ls_scf_env%max_diis,error=error)
    CALL section_vals_val_get(ls_scf_section,"NMIXING",i_val=ls_scf_env%nmixing,error=error)
    CALL section_vals_val_get(ls_scf_section,"EPS_DIIS",r_val=ls_scf_env%eps_diis,error=error)
    CALL section_vals_val_get(ls_scf_section,"EPS_SCF",r_val=ls_scf_env%eps_scf,error=error)
    CALL section_vals_val_get(ls_scf_section,"EPS_FILTER",r_val=ls_scf_env%eps_filter,error=error)
    CALL section_vals_val_get(ls_scf_section,"MU",r_val=mu,error=error)
    CALL section_vals_val_get(ls_scf_section,"FIXED_MU",l_val=ls_scf_env%fixed_mu,error=error)
    ls_scf_env%mu_spin=mu
    CALL section_vals_val_get(ls_scf_section,"MIXING_FRACTION",r_val=ls_scf_env%mixing_fraction,error=error)
    CALL section_vals_val_get(ls_scf_section,"MAX_SCF",i_val=ls_scf_env%max_scf,error=error)
    CALL section_vals_val_get(ls_scf_section,"S_PRECONDITIONER",i_val=ls_scf_env%s_preconditioner_type,error=error)
    CALL section_vals_val_get(ls_scf_section,"MATRIX_CLUSTER_TYPE",i_val=ls_scf_env%ls_mstruct%cluster_type,error=error)
    CALL section_vals_val_get(ls_scf_section,"SINGLE_PRECISION_MATRICES",l_val=ls_scf_env%ls_mstruct%single_precision,error=error)
    CALL section_vals_val_get(ls_scf_section,"S_INVERSION",i_val=ls_scf_env%s_inversion_type,error=error)
    CALL section_vals_val_get(ls_scf_section,"REPORT_ALL_SPARSITIES",l_val=ls_scf_env%report_all_sparsities,error=error)
    CALL section_vals_val_get(ls_scf_section,"PERFORM_MU_SCAN",l_val=ls_scf_env%perform_mu_scan,error=error)
    CALL section_vals_val_get(ls_scf_section,"PURIFICATION_METHOD",i_val=ls_scf_env%purification_method,error=error)
    CALL section_vals_val_get(ls_scf_section,"DYNAMIC_THRESHOLD",l_val=ls_scf_env%dynamic_threshold,error=error)
    CALL section_vals_val_get(ls_scf_section,"NON_MONOTONIC",l_val=ls_scf_env%non_monotonic,error=error)
    CALL section_vals_val_get(ls_scf_section,"SIGN_SQRT_ORDER",i_val=ls_scf_env%sign_sqrt_order,error=error)
    CALL section_vals_val_get(ls_scf_section,"EXTRAPOLATION_ORDER",i_val=ls_scf_env%extrapolation_order,error=error)
    CALL section_vals_val_get(ls_scf_section,"RESTART_READ",l_val=ls_scf_env%restart_read,error=error)
    CALL section_vals_val_get(ls_scf_section,"RESTART_WRITE",l_val=ls_scf_env%restart_write,error=error)
    CALL section_vals_val_get(ls_scf_section,"EPS_LANCZOS",r_val=ls_scf_env%eps_lanczos,error=error)
    CALL section_vals_val_get(ls_scf_section,"MAX_ITER_LANCZOS",i_val=ls_scf_env%max_iter_lanczos,error=error)

    CALL section_vals_get(curvy_section, explicit=ls_scf_env%curvy_steps, error=error)
    CALL section_vals_val_get(curvy_section,"LINE_SEARCH",i_val=ls_scf_env%curvy_data%line_search_type,error=error)
    CALL section_vals_val_get(curvy_section,"N_BCH_HISTORY",i_val=ls_scf_env%curvy_data%n_bch_hist,error=error)
    CALL section_vals_val_get(curvy_section,"MIN_HESSIAN_SHIFT",r_val=ls_scf_env%curvy_data%min_shift,error=error)
    CALL section_vals_val_get(curvy_section,"FILTER_FACTOR",r_val=ls_scf_env%curvy_data%filter_factor,error=error)
    CALL section_vals_val_get(curvy_section,"FILTER_FACTOR_SCALE",r_val=ls_scf_env%curvy_data%scale_filter,error=error)
    CALL section_vals_val_get(curvy_section,"MIN_FILTER",r_val=ls_scf_env%curvy_data%min_filter,error=error)

    ls_scf_env%extrapolation_order=MAX(0,ls_scf_env%extrapolation_order)

    chebyshev_section => section_vals_get_subs_vals(input,"DFT%LS_SCF%CHEBYSHEV",error=error)
    CALL section_vals_get(chebyshev_section,explicit=ls_scf_env%chebyshev%compute_chebyshev,error=error)
    IF (ls_scf_env%chebyshev%compute_chebyshev) THEN
      CALL section_vals_val_get(chebyshev_section,"N_CHEBYSHEV",i_val=ls_scf_env%chebyshev%n_chebyshev,error=error)
      CALL section_vals_val_get(chebyshev_section,"DOS%N_GRIDPOINTS",i_val=ls_scf_env%chebyshev%n_gridpoint_dos,error=error)

      ls_scf_env%chebyshev%print_key_dos  => &
          section_vals_get_subs_vals(chebyshev_section,"DOS",error=error)
      CALL section_vals_retain(ls_scf_env%chebyshev%print_key_dos,error=error)

      ls_scf_env%chebyshev%print_key_cube => &
          section_vals_get_subs_vals(chebyshev_section,"PRINT_SPECIFIC_E_DENSITY_CUBE",error=error)
      CALL section_vals_retain(ls_scf_env%chebyshev%print_key_cube,error=error)
    ENDIF

    mixing_section => section_vals_get_subs_vals(input,"DFT%LS_SCF%RHO_MIXING",error=error)
    CALL section_vals_get(mixing_section, explicit=ls_scf_env%do_rho_mixing, error=error)

    CALL section_vals_val_get(mixing_section,"METHOD",i_val=ls_scf_env%density_mixing_method,error=error)
    CALL cp_assert(.NOT.(ls_scf_env%ls_diis.AND.ls_scf_env%do_rho_mixing),cp_failure_level,cp_assertion_failed,&
                   routineP,"LS_DIIS and RHO_MIXING are not compatible."//&
                   CPSourceFileRef,only_ionode=.TRUE.)

    pexsi_section => section_vals_get_subs_vals(input,"DFT%LS_SCF%PEXSI",error=error)
    CALL section_vals_get(pexsi_section, error=error)

    IF (ls_scf_env%purification_method .EQ. ls_scf_pexsi) THEN
      CALL pexsi_init_read_input(pexsi_section,ls_scf_env%pexsi, error)
      IF (.NOT. ls_scf_env%restart_write) THEN
        ! Turn off S inversion (not used for PEXSI).
        ! Methods such as purification must thus be avoided... which is OK, as the density matrix computed in pexsi is
        ! a finite temperature one, and thus not idempotent
        ls_scf_env%s_inversion_type = ls_s_inversion_none
        ! PEXSI needs the sparsity pattern of S to be fixed by the upper bound ~ Int[|phi_a||phi_b|], 
        ! they can not be filtered based on magnitude, as small elements in S (e.g. symmetry) do not necessarily
        ! correspond to small elements in the density matrix, with non-zero contributions to the total density.
        ! the easiest way to make sure S is untouched is to set eps_filter to zero (which should be inconsequential,
        ! as a run based on pexsi should execute exactly zero multiplications)
        ls_scf_env%eps_filter = 0.0_dp
      END IF
    ENDIF

    ! Turn off S inversion and set eps_filter to zero for transport
    IF (ls_scf_env%do_transport .AND. .NOT. ls_scf_env%restart_write) THEN
      ls_scf_env%s_inversion_type = ls_s_inversion_none 
      ls_scf_env%eps_filter = 0.0_dp
    END IF

    SELECT CASE(ls_scf_env%s_inversion_type)
    CASE(ls_s_inversion_sign_sqrt)
        ls_scf_env%needs_s_inv=.TRUE.
        ls_scf_env%use_s_sqrt=.TRUE.
    CASE(ls_s_inversion_hotelling)
        ls_scf_env%needs_s_inv=.TRUE.
        ls_scf_env%use_s_sqrt=.FALSE.
    CASE(ls_s_inversion_none)
        ls_scf_env%needs_s_inv=.FALSE.
        ls_scf_env%use_s_sqrt=.FALSE.
    CASE DEFAULT
        CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    SELECT CASE(ls_scf_env%s_preconditioner_type)
    CASE(ls_s_preconditioner_none)
      ls_scf_env%has_s_preconditioner=.FALSE.
    CASE DEFAULT
      ls_scf_env%has_s_preconditioner=.TRUE.
    END SELECT

    ! verify some requirements for the curvy steps
    CALL cp_assert(.NOT.(ls_scf_env%curvy_steps .AND. ls_scf_env%purification_method .EQ. ls_scf_pexsi), cp_failure_level, &
                       cp_assertion_failed, routineP, "CURVY_STEPS can not be used together with PEXSI. ", &
                       only_ionode=.TRUE., error=error)
    CALL cp_assert(.NOT.(ls_scf_env%curvy_steps .AND. ls_scf_env%do_transport), cp_failure_level, &
                       cp_assertion_failed, routineP, "CURVY_STEPS can not be used together with TRANSPORT. ", &
                       only_ionode=.TRUE., error=error)
    CALL cp_assert(.NOT.(ls_scf_env%curvy_steps.AND.ls_scf_env%has_s_preconditioner),cp_failure_level,cp_assertion_failed,&
                      routineP,"S Preconditioning not implemented in combination with CURVY_STEPS. "//&
                       CPSourceFileRef,&
                       only_ionode=.TRUE.)
    CALL cp_assert(.NOT.(ls_scf_env%curvy_steps.AND..NOT.ls_scf_env%use_s_sqrt),cp_failure_level,cp_assertion_failed,&
                      routineP,"CURVY_STEPS requires the use of the sqrt inversion. "//&
                       CPSourceFileRef,&
                       only_ionode=.TRUE.)

    ! an undocumented feature ... allows for just doing the initial guess, no expensive stuff
    IF (ls_scf_env%max_scf<0) THEN
        ls_scf_env%needs_s_inv=.FALSE.
        ls_scf_env%use_s_sqrt=.FALSE.
      ls_scf_env%has_s_preconditioner=.FALSE.
    ENDIF

    IF (unit_nr>0) THEN
       WRITE(unit_nr,'()')
       WRITE(unit_nr,'(T2,A,A,A)') REPEAT("-",30)," Linear scaling SCF ",REPEAT("-",29)
       WRITE(unit_nr,'(T2,A,T38,E20.3)') "eps_scf:",ls_scf_env%eps_scf
       WRITE(unit_nr,'(T2,A,T38,E20.3)') "eps_filter:",ls_scf_env%eps_filter
       IF (ls_scf_env%do_rho_mixing) THEN
         IF (ls_scf_env%density_mixing_method>0) THEN
           NULLIFY(section)
           CALL create_mixing_section(section,ls_scf=.TRUE.,error=error)
           keyword => section_get_keyword(section,"METHOD",error=error)
           CALL keyword_get(keyword,enum=enum,error=error)
           WRITE (unit_nr,"(T2,A,T38,A20)")&
                  "Density mixing in g-space:",ADJUSTR(TRIM(enum_i2c(enum,&
                                  ls_scf_env%density_mixing_method,error=error)))
           CALL section_release(section,error=error)
         END IF
       ELSE
         WRITE(unit_nr,'(T2,A,T38,E20.3)') "mixing_fraction:",ls_scf_env%mixing_fraction
       ENDIF
       WRITE(unit_nr,'(T2,A,T38,I20)') "max_scf:",ls_scf_env%max_scf
       IF (ls_scf_env%ls_diis) THEN
       WRITE(unit_nr,'(T2,A,T38,I20)')   "DIIS: max_diis:",ls_scf_env%max_diis
       WRITE(unit_nr,'(T2,A,T38,E20.3)') "DIIS: eps_diis:",ls_scf_env%eps_diis
       WRITE(unit_nr,'(T2,A,T38,I20)')   "DIIS: ini_diis:",ls_scf_env%iter_ini_diis
       WRITE(unit_nr,'(T2,A,T38,I20)')   "DIIS: nmixing:" ,ls_scf_env%nmixing
       ENDIF
       WRITE(unit_nr,'(T2,A,T38,L20)') "fixed chemical potential (mu)",ls_scf_env%fixed_mu
       WRITE(unit_nr,'(T2,A,T38,L20)') "has unit metric:",ls_scf_env%has_unit_metric
       WRITE(unit_nr,'(T2,A,T38,L20)') "Computing inv(S):",ls_scf_env%needs_s_inv
       WRITE(unit_nr,'(T2,A,T38,L20)') "Computing sqrt(S):",ls_scf_env%use_s_sqrt
       WRITE(unit_nr,'(T2,A,T38,L20)') "Computing s preconditioner ",ls_scf_env%has_s_preconditioner
       WRITE(unit_nr,'(T2,A,T38,I20)') "sign sqrt order:",ls_scf_env%sign_sqrt_order
       WRITE(unit_nr,'(T2,A,T38,I20)') "Extrapolation order:",ls_scf_env%extrapolation_order
       WRITE(unit_nr,'(T2,A,T38,L20)') "Use single precision matrices",ls_scf_env%ls_mstruct%single_precision

       SELECT CASE(ls_scf_env%s_preconditioner_type)
       CASE(ls_s_preconditioner_none)
           WRITE(unit_nr,'(T2,A,T38,A20)') "S preconditioner type ","NONE"
       CASE(ls_s_preconditioner_atomic)
           WRITE(unit_nr,'(T2,A,T38,A20)') "S preconditioner type ","ATOMIC"
       CASE(ls_s_preconditioner_molecular)
           WRITE(unit_nr,'(T2,A,T38,A20)') "S preconditioner type ","MOLECULAR"
       END SELECT

       IF (ls_scf_env%curvy_steps) THEN
          WRITE(unit_nr,'(T2,A,T38,A30)') "Using curvy steps to optimize the density matrix"
          CALL cite_reference(Shao2003)
       ENDIF

       SELECT CASE(ls_scf_env%purification_method)
       CASE(ls_scf_ns)
           WRITE(unit_nr,'(T2,A,T38,A30)') "Purification method","Newton-Schulz sign iter"
       CASE(ls_scf_tc2)
           CALL cite_reference(Niklasson2014)
           WRITE(unit_nr,'(T2,A,T38,A30)') "Purification method","Trace conserving 2nd  order"
       CASE(ls_scf_trs4)
           CALL cite_reference(Niklasson2003)
           WRITE(unit_nr,'(T2,A,T38,A30)') "Purification method","Trace resetting 4th order"
       CASE(ls_scf_pexsi)
           CALL cite_reference(Lin2009)
           CALL cite_reference(Lin2013)
           WRITE(unit_nr,'(T2,A,T38,A20)') "Purification method","PEXSI"
       CASE DEFAULT
           CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT

       SELECT CASE(ls_scf_env%ls_mstruct%cluster_type)
       CASE(ls_cluster_atomic)
           WRITE(unit_nr,'(T2,A,T38,A20)') "Cluster type","ATOMIC"
       CASE(ls_cluster_molecular)
           WRITE(unit_nr,'(T2,A,T38,A20)') "Cluster type","MOLECULAR"
       END SELECT

       IF (ls_scf_env%chebyshev%compute_chebyshev) THEN
          WRITE(unit_nr,'(T2,A,T38,A20)') "Computing Chebyshev","TRUE"
          WRITE(unit_nr,'(T2,A,T38,I20)') "N_CHEBYSHEV:",ls_scf_env%chebyshev%n_chebyshev
          WRITE(unit_nr,'(T2,A,T38,I20)') "N_GRIDPOINT_DOS:",ls_scf_env%chebyshev%n_gridpoint_dos
       ELSE
          WRITE(unit_nr,'(T2,A,T38,A20)') "Computing Chebyshev","FALSE"
       ENDIF

       WRITE(unit_nr,'(T2,A)') REPEAT("-",79)
       WRITE(unit_nr,'()')
       CALL m_flush(unit_nr)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE ls_scf_init_read_write_input

! *****************************************************************************
!> \brief Main SCF routine. Can we keep it clean ?
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param error ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_main(qs_env,ls_scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_main', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iscf, ispin, &
                                                nelectron_spin_real, nmixing, &
                                                nspin, unit_nr
    LOGICAL                                  :: check_convergence, diis_step, &
                                                do_transport, should_stop
    REAL(KIND=dp)                            :: energy_diff, energy_new, &
                                                energy_old, eps_diis, t1, t2
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_dbcsr_type), ALLOCATABLE, &
      DIMENSION(:)                           :: matrix_ks_deviation, &
                                                matrix_mixing_old
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_diis_buffer_type_sparse), &
      POINTER                                :: diis_buffer
    TYPE(transport_env_type), POINTER        :: transport_env

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    nspin=ls_scf_env%nspins

    ! old quantities, useful for mixing
    ALLOCATE(matrix_mixing_old(nspin), matrix_ks_deviation(nspin))
    DO ispin=1,nspin
       CALL cp_dbcsr_init(matrix_mixing_old(ispin),error=error)
       CALL cp_dbcsr_create(matrix_mixing_old(ispin),template=ls_scf_env%matrix_ks(ispin), error=error)

       CALL cp_dbcsr_init(matrix_ks_deviation(ispin),error=error)
       CALL cp_dbcsr_create(matrix_ks_deviation(ispin),template=ls_scf_env%matrix_ks(ispin), error=error)
       CALL cp_dbcsr_set(matrix_ks_deviation(ispin),0.0_dp,error=error)
    ENDDO
    ls_scf_env%homo_spin(:) = 0.0_dp
    ls_scf_env%lumo_spin(:) = 0.0_dp

    energy_old=0.0_dp
    IF(ls_scf_env%scf_history%istore>0)energy_old=ls_scf_env%energy_init
    check_convergence=.TRUE.
    iscf=0
    IF (ls_scf_env%ls_diis) THEN
     diis_step = .FALSE.
     eps_diis  = ls_scf_env%eps_diis
     nmixing   = ls_scf_env%nmixing
     NULLIFY(diis_buffer)
     IF (.NOT.ASSOCIATED(diis_buffer)) THEN
        CALL qs_diis_b_create_sparse(diis_buffer, &
                                     nbuffer=ls_scf_env%max_diis,  &
                                     error=error)
     END IF
     CALL qs_diis_b_clear_sparse(diis_buffer,error=error)
     CALL get_qs_env(qs_env,matrix_s=matrix_s,error=error)
    END IF

    CALL get_qs_env(qs_env, transport_env=transport_env, do_transport=do_transport, error=error)

    ! the real SCF loop
    DO

      ! check on max SCF or timing/exit
      CALL external_control(should_stop,"SCF",start_time=qs_env%start_time,target_time=qs_env%target_time,error=error)
      IF (should_stop .OR. iscf>=ls_scf_env%max_scf) THEN
         IF (unit_nr>0) WRITE(unit_nr,'(T2,A)') "SCF not converged! "
         EXIT
      ENDIF

      t1 = m_walltime()
      iscf=iscf+1

      ! first get a copy of the current KS matrix
      CALL get_qs_env(qs_env, matrix_ks=matrix_ks, error=error)
      DO ispin=1,nspin
         CALL matrix_qs_to_ls(ls_scf_env%matrix_ks(ispin),matrix_ks(ispin)%matrix,ls_scf_env%ls_mstruct,error=error)
         IF (ls_scf_env%has_s_preconditioner) THEN
             CALL apply_matrix_preconditioner(ls_scf_env%matrix_ks(ispin),"forward", &
                            ls_scf_env%matrix_bs_sqrt,ls_scf_env%matrix_bs_sqrt_inv,error)
         ENDIF
         CALL cp_dbcsr_filter(ls_scf_env%matrix_ks(ispin),ls_scf_env%eps_filter,error=error)
      ENDDO
      ! run curvy steps if required. Needs an idempotent DM (either perification or restart)
      IF((iscf>1.OR.ls_scf_env%scf_history%istore>0).AND.ls_scf_env%curvy_steps)THEN
         CALL dm_ls_curvy_optimization(ls_scf_env,energy_old,check_convergence,error)
      ELSE
         ! turn the KS matrix in a density matrix
         DO ispin=1,nspin
            IF (ls_scf_env%do_rho_mixing) THEN
              CALL cp_dbcsr_copy(matrix_mixing_old(ispin), ls_scf_env%matrix_ks(ispin), error=error)
            ELSE
             IF (iscf==1) THEN
                ! initialize the mixing matrix with the current state if needed
                CALL cp_dbcsr_copy(matrix_mixing_old(ispin), ls_scf_env%matrix_ks(ispin), error=error)
             ELSE
              IF (ls_scf_env%ls_diis) THEN ! ------- IF-DIIS+MIX--- START
               IF (diis_step.and.(iscf-1).ge.ls_scf_env%iter_ini_diis) THEN
                IF (unit_nr>0) THEN
                 WRITE(unit_nr,'(A61)') &
                        '*************************************************************'
                 WRITE(unit_nr,'(A50,2(I3,A1),L1,A1)') &
                        " Using DIIS mixed KS:  (iscf,INI_DIIS,DIIS_STEP)=(" , &
                        iscf,",",ls_scf_env%iter_ini_diis,",",diis_step,")"
                 WRITE(unit_nr,'(A52)') &
                        " KS_nw= DIIS-Linear-Combination-Previous KS matrices"
                 WRITE(unit_nr,'(61A)') &
                        "*************************************************************"
                ENDIF
                CALL cp_dbcsr_copy(matrix_mixing_old(ispin),   & ! out
                                   ls_scf_env%matrix_ks(ispin),& ! in
                                   error=error)
               ELSE
                IF (unit_nr>0) THEN
                 WRITE(unit_nr,'(A57)') &
                       "*********************************************************"
                 WRITE(unit_nr,'(A23,F5.3,A25,I3)') &
                       " Using MIXING_FRACTION=",ls_scf_env%mixing_fraction, &
                       " to mix KS matrix:  iscf=",iscf 
                 WRITE(unit_nr,'(A7,F5.3,A6,F5.3,A7)') &
                       " KS_nw=",ls_scf_env%mixing_fraction,"*KS + ", &
                       1.0_dp-ls_scf_env%mixing_fraction,"*KS_old"
                 WRITE(unit_nr,'(A57)') &
                       "*********************************************************"
                ENDIF
                ! perform the mixing of ks matrices
                CALL cp_dbcsr_add(matrix_mixing_old(ispin)   ,       &
                                  ls_scf_env%matrix_ks(ispin),       &
                                  1.0_dp-ls_scf_env%mixing_fraction, &
                                  ls_scf_env%mixing_fraction,        &
                                  error=error)
               ENDIF
              ELSE ! otherwise
                IF (unit_nr>0) THEN
                 WRITE(unit_nr,'(A57)') &
                       "*********************************************************"
                 WRITE(unit_nr,'(A23,F5.3,A25,I3)') &
                       " Using MIXING_FRACTION=", ls_scf_env%mixing_fraction, &
                       " to mix KS matrix:  iscf=",iscf
                 WRITE(unit_nr,'(A7,F5.3,A6,F5.3,A7)') &
                       " KS_nw=",ls_scf_env%mixing_fraction,"*KS + ", &
                       1.0_dp-ls_scf_env%mixing_fraction,"*KS_old"
                 WRITE(unit_nr,'(A57)') &
                       "*********************************************************"
                ENDIF
                ! perform the mixing of ks matrices
                CALL cp_dbcsr_add(matrix_mixing_old(ispin)   ,       &
                                  ls_scf_env%matrix_ks(ispin),       &
                                  1.0_dp-ls_scf_env%mixing_fraction, &
                                  ls_scf_env%mixing_fraction,        &
                                  error=error)
              ENDIF ! ------- IF-DIIS+MIX--- END
             ENDIF
            ENDIF

            ! compute the density matrix that matches it
            ! we need the proper number of states
            nelectron_spin_real=ls_scf_env%nelectron_spin(ispin)
            IF (ls_scf_env%nspins==1) nelectron_spin_real=nelectron_spin_real/2

            IF (do_transport) THEN
               CALL cp_assert(.NOT.(ls_scf_env%has_s_preconditioner), cp_failure_level, cp_assertion_failed,&
                              routineP, "NOT YET IMPLEMENTED with S preconditioner. ", only_ionode=.TRUE., error=error)
               CALL cp_assert(ls_scf_env%ls_mstruct%cluster_type .EQ. ls_cluster_atomic, cp_failure_level,&
                              cp_assertion_failed, routineP, "NOT YET IMPLEMENTED with molecular clustering. ",& 
                              only_ionode=.TRUE., error=error)

               ! get the current Kohn-Sham matrix (ks) and return matrix_p evaluated using an external C routine
               CALL external_scf_method(transport_env, ls_scf_env%matrix_s, matrix_mixing_old(ispin), &
                                        ls_scf_env%matrix_p(ispin), nelectron_spin_real, ls_scf_env%natoms, error)

            ELSE 
               SELECT CASE(ls_scf_env%purification_method)
               CASE(ls_scf_ns)
                 CALL density_matrix_sign(ls_scf_env%matrix_p(ispin),ls_scf_env%mu_spin(ispin), ls_scf_env%fixed_mu, &
                                          matrix_mixing_old(ispin),ls_scf_env%matrix_s, ls_scf_env%matrix_s_inv, &
                                          nelectron_spin_real,ls_scf_env%eps_filter,error)
               CASE(ls_scf_tc2)
                 CALL density_matrix_tc2(ls_scf_env%matrix_p(ispin), matrix_mixing_old(ispin),  ls_scf_env%matrix_s_sqrt_inv,&
                                          nelectron_spin_real, ls_scf_env%eps_filter, ls_scf_env%homo_spin(ispin),&
                                          ls_scf_env%lumo_spin(ispin), non_monotonic=ls_scf_env%non_monotonic, &
                                          eps_lanczos=ls_scf_env%eps_lanczos, max_iter_lanczos=ls_scf_env%max_iter_lanczos,&
                                          error=error)
               CASE(ls_scf_trs4)
                 CALL density_matrix_trs4(ls_scf_env%matrix_p(ispin), matrix_mixing_old(ispin),  ls_scf_env%matrix_s_sqrt_inv,&
                                          nelectron_spin_real, ls_scf_env%eps_filter, ls_scf_env%homo_spin(ispin),&
                                          ls_scf_env%lumo_spin(ispin), ls_scf_env%mu_spin(ispin), &
                                          dynamic_threshold=ls_scf_env%dynamic_threshold,&
                                          matrix_ks_deviation=matrix_ks_deviation(ispin), &
                                          eps_lanczos=ls_scf_env%eps_lanczos, max_iter_lanczos=ls_scf_env%max_iter_lanczos,&
                                          error=error)
               CASE(ls_scf_pexsi)
                 CALL cp_assert(.NOT.(ls_scf_env%has_s_preconditioner),cp_failure_level,cp_assertion_failed, routineP, &
                                "S preconditioning not implemented in combination with the PEXSI library. ", &
                                only_ionode=.TRUE., error=error)
                 CALL cp_assert(ls_scf_env%ls_mstruct%cluster_type .EQ. ls_cluster_atomic, cp_failure_level, cp_assertion_failed, & 
                                routineP,"Molecular clustering not implemented in combination with the PEXSI library. ", &
                                only_ionode=.TRUE., error=error)
                 CALL cp_assert(.NOT. ls_scf_env%ls_mstruct%single_precision,cp_failure_level,cp_assertion_failed, routineP, &
                                "PEXSI library requires double precision datatype. ", only_ionode=.TRUE.,error=error)
                 CALL density_matrix_pexsi(ls_scf_env%pexsi, ls_scf_env%matrix_p(ispin), ls_scf_env%pexsi%matrix_w(ispin),&
                                           ls_scf_env%pexsi%kTS(ispin), matrix_mixing_old(ispin), ls_scf_env%matrix_s,& 
                                           nelectron_spin_real,ls_scf_env%mu_spin(ispin),iscf,ispin,error)
               END SELECT
            END IF
 
            IF (ls_scf_env%has_s_preconditioner) THEN
                CALL apply_matrix_preconditioner(ls_scf_env%matrix_p(ispin),"forward", &
                               ls_scf_env%matrix_bs_sqrt,ls_scf_env%matrix_bs_sqrt_inv,error)
            ENDIF
            CALL cp_dbcsr_filter(ls_scf_env%matrix_p(ispin),ls_scf_env%eps_filter,error=error)
   
            IF (ls_scf_env%nspins==1) CALL cp_dbcsr_scale(ls_scf_env%matrix_p(ispin),2.0_dp,error=error)

         ENDDO
      END IF

      ! compute the corresponding new energy KS matrix and new energy
      CALL ls_scf_dm_to_ks(qs_env,ls_scf_env,energy_new,iscf,error)

      IF (ls_scf_env%purification_method .EQ. ls_scf_pexsi) THEN 
        CALL pexsi_to_qs(ls_scf_env, qs_env, kTS = ls_scf_env%pexsi%kTS, error = error)
      ENDIF

      ! report current SCF loop
      energy_diff=energy_new-energy_old
      energy_old=energy_new

      t2 = m_walltime()
      IF (unit_nr>0) THEN
         WRITE(unit_nr,*)
         WRITE(unit_nr,'(T2,A,I6,F20.9,F20.9,F12.6)') "SCF",iscf, energy_new,energy_diff, t2-t1
         WRITE(unit_nr,*)
         CALL m_flush(unit_nr)
      ENDIF

      ! exit criterium on the energy only for the time being
      IF (check_convergence.AND.ABS(energy_diff)<ls_scf_env%eps_scf*ls_scf_env%nelectron_total) EXIT

      IF (ls_scf_env%ls_diis) THEN
! diis_buffer, buffer with 1) Kohn-Sham history matrix, 
!                          2) KS error history matrix (f=KPS-SPK),
!                          3) B matrix (for finding DIIS weighting coefficients)
       CALL qs_diis_b_step_4lscf(diis_buffer,qs_env,ls_scf_env,unit_nr, & 
                                 iscf,diis_step,eps_diis,nmixing,matrix_s(1)%matrix, & 
                                 ls_scf_env%eps_filter,error)               
      ENDIF

      IF (ls_scf_env%purification_method .EQ. ls_scf_pexsi) THEN
        CALL pexsi_set_convergence_tolerance(ls_scf_env%pexsi, energy_diff, &
                                             ls_scf_env%eps_scf*ls_scf_env%nelectron_total, &
                                             ! initialize in second scf step of first SCF cycle:
                                             (iscf.EQ.2).AND.(ls_scf_env%scf_history%istore.EQ.0), &
                                             check_convergence, error)
      ENDIF

    ENDDO

    ! free storage
    IF (ls_scf_env%ls_diis) THEN
       CALL qs_diis_b_release_sparse(diis_buffer,error)
    ENDIF
    DO ispin=1,nspin
       CALL cp_dbcsr_release(matrix_mixing_old(ispin),error=error)
       CALL cp_dbcsr_release(matrix_ks_deviation(ispin),error=error)
    ENDDO
    DEALLOCATE(matrix_mixing_old, matrix_ks_deviation)

    CALL timestop(handle)

  END SUBROUTINE ls_scf_main


! *****************************************************************************
!> \brief after SCF we have a density matrix, and the self consistent KS matrix
!>        analyze its properties.
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param error ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_post(qs_env,ls_scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_post', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, unit_nr
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_w
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! store the matrix for a next scf run
    CALL ls_scf_store_result(ls_scf_env,error)

    ! write homo and lumo energy (if not already part of the output)
    IF (ls_scf_env%curvy_steps) THEN
       CALL post_scf_homo_lumo(ls_scf_env,error)
    ENDIF

    ! compute the matrix_w if associated
    CALL get_qs_env(qs_env,matrix_w=matrix_w,error=error)
    IF (ASSOCIATED(matrix_w)) THEN
      IF (ls_scf_env%purification_method .EQ. ls_scf_pexsi) THEN
        CALL pexsi_to_qs(ls_scf_env, qs_env, matrix_w = ls_scf_env%pexsi%matrix_w, error=error)
      ELSE
        CALL calculate_w_matrix(matrix_w,ls_scf_env,error=error) 
      ENDIF
    ENDIF

    ! compute properties

    IF (ls_scf_env%perform_mu_scan) CALL post_scf_mu_scan(ls_scf_env,error)

    IF (ls_scf_env%report_all_sparsities) CALL post_scf_sparsities(ls_scf_env,error)

    CALL write_mo_free_results(qs_env,error)

    IF (ls_scf_env%chebyshev%compute_chebyshev) CALL compute_chebyshev(qs_env,ls_scf_env,error)

    IF (.TRUE.) CALL post_scf_experiment(ls_scf_env,error)

    CALL qs_scf_post_moments(qs_env%input, logger, qs_env, unit_nr, error=error)

    ! clean up used data

    CALL cp_dbcsr_release(ls_scf_env%matrix_s,error=error)
    CALL deallocate_curvy_data(ls_scf_env%curvy_data,error)

    IF (ls_scf_env%has_s_preconditioner) THEN
       CALL cp_dbcsr_release(ls_scf_env%matrix_bs_sqrt,error=error)
       CALL cp_dbcsr_release(ls_scf_env%matrix_bs_sqrt_inv,error=error)
    ENDIF

    IF (ls_scf_env%needs_s_inv) THEN
        CALL cp_dbcsr_release(ls_scf_env%matrix_s_inv,error=error)
    ENDIF

    IF (ls_scf_env%use_s_sqrt) THEN
        CALL cp_dbcsr_release(ls_scf_env%matrix_s_sqrt,error=error)
        CALL cp_dbcsr_release(ls_scf_env%matrix_s_sqrt_inv,error=error)
    ENDIF

    DO ispin=1,SIZE(ls_scf_env%matrix_p)
       CALL cp_dbcsr_release(ls_scf_env%matrix_p(ispin),error=error)
    ENDDO
    DEALLOCATE(ls_scf_env%matrix_p)

    DO ispin=1,SIZE(ls_scf_env%matrix_ks)
       CALL cp_dbcsr_release(ls_scf_env%matrix_ks(ispin),error=error)
    ENDDO
    DEALLOCATE(ls_scf_env%matrix_ks)

    IF (ls_scf_env%purification_method .EQ. ls_scf_pexsi) &
      CALL pexsi_finalize_scf(ls_scf_env%pexsi, ls_scf_env%mu_spin, error)

    CALL timestop(handle)

  END SUBROUTINE ls_scf_post

! *****************************************************************************
!> \brief Compute the HOMO LUMO energies post SCF
!> \param ls_scf_env ...
!> \param error ...
!> \par History
!>       2013.06 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE post_scf_homo_lumo(ls_scf_env,error)
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'post_scf_homo_lumo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin, unit_nr
    LOGICAL                                  :: converged, failure
    REAL(KIND=dp)                            :: eps_max, eps_min, homo, lumo
    TYPE(cp_dbcsr_type)                      :: matrix_k, matrix_p
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    IF (unit_nr>0) WRITE(unit_nr,'(T2,A)') ""

    ! TODO: remove these limitations
    CPPrecondition(.NOT.ls_scf_env%has_s_preconditioner,cp_failure_level,routineP,error,failure)
    CPPrecondition(ls_scf_env%use_s_sqrt,cp_failure_level,routineP,error,failure)

    nspin=ls_scf_env%nspins

    CALL cp_dbcsr_init(matrix_p,error=error)
    CALL cp_dbcsr_create(matrix_p,template=ls_scf_env%matrix_p(1), matrix_type=dbcsr_type_no_symmetry, error=error)

    CALL cp_dbcsr_init(matrix_k,error=error)
    CALL cp_dbcsr_create(matrix_k,template=ls_scf_env%matrix_p(1), matrix_type=dbcsr_type_no_symmetry, error=error)

    DO ispin=1,nspin
       ! ortho basis ks
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt_inv, ls_scf_env%matrix_ks(ispin),&
                              0.0_dp, matrix_k, filter_eps=ls_scf_env%eps_filter,error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_k, ls_scf_env%matrix_s_sqrt_inv, &
                              0.0_dp, matrix_k, filter_eps=ls_scf_env%eps_filter,error=error)

       ! extremal eigenvalues ks
       CALL cp_dbcsr_arnoldi_extremal(matrix_k, eps_max, eps_min, max_iter=ls_scf_env%max_iter_lanczos, &
                               threshold=ls_scf_env%eps_lanczos, converged=converged, error=error)

       ! ortho basis p
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt, ls_scf_env%matrix_p(ispin),&
                              0.0_dp, matrix_p, filter_eps=ls_scf_env%eps_filter,error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_p, ls_scf_env%matrix_s_sqrt, &
                              0.0_dp, matrix_p, filter_eps=ls_scf_env%eps_filter,error=error)
       IF(nspin==1)CALL cp_dbcsr_scale(matrix_p,0.5_dp,error=error)

       ! go compute homo lumo
       CALL compute_homo_lumo(matrix_k,matrix_p,eps_min,eps_max,ls_scf_env%eps_filter, &
                             ls_scf_env%max_iter_lanczos,ls_scf_env%eps_lanczos,homo,lumo,unit_nr,error)

    ENDDO

    CALL cp_dbcsr_release(matrix_p,error=error)
    CALL cp_dbcsr_release(matrix_k,error=error)

    CALL timestop(handle)

  END SUBROUTINE post_scf_homo_lumo

! *****************************************************************************
!> \brief Compute the density matrix for various values of the chemical potential
!> \param ls_scf_env ...
!> \param error ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE post_scf_mu_scan(ls_scf_env,error)
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'post_scf_mu_scan', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, imu, ispin, &
                                                nelectron_spin_real, nmu, &
                                                nspin, unit_nr
    REAL(KIND=dp)                            :: mu, t1, t2, trace
    TYPE(cp_dbcsr_type)                      :: matrix_p
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    nspin=ls_scf_env%nspins

    CALL cp_dbcsr_init(matrix_p,error=error)
    CALL cp_dbcsr_create(matrix_p,template=ls_scf_env%matrix_p(1),error=error)

    nmu=10
    DO imu=0,nmu

       t1 = m_walltime()

       mu = -0.4_dp+imu*(0.1_dp+0.4_dp)/nmu

       IF (unit_nr>0) WRITE(unit_nr,*) "------- starting with mu ",mu

       DO ispin=1,nspin
          ! we need the proper number of states
          nelectron_spin_real=ls_scf_env%nelectron_spin(ispin)
          IF (ls_scf_env%nspins==1) nelectron_spin_real=nelectron_spin_real/2

          CALL density_matrix_sign_fixed_mu(matrix_p,trace,mu, &
                                            ls_scf_env%matrix_ks(ispin),ls_scf_env%matrix_s,&
                                            ls_scf_env%matrix_s_inv,ls_scf_env%eps_filter,error)
       ENDDO

       t2 = m_walltime()

       IF (unit_nr>0) WRITE(unit_nr,*) " obtained " ,mu,trace,t2-t1

    ENDDO

    CALL cp_dbcsr_release(matrix_p,error=error)

    CALL timestop(handle)

  END SUBROUTINE post_scf_mu_scan

! *****************************************************************************
!> \brief Report on the sparsities of various interesting matrices.
!>
!> \param ls_scf_env ...
!> \param error ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE post_scf_sparsities(ls_scf_env,error)
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'post_scf_sparsities', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: title
    INTEGER                                  :: handle, ispin, nspin, unit_nr
    TYPE(cp_dbcsr_type)                      :: matrix_tmp1, matrix_tmp2
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    nspin=ls_scf_env%nspins

    IF (unit_nr>0) THEN
       WRITE(unit_nr,'()')
       WRITE(unit_nr,'(T2,A,E17.3)') "Sparsity reports for eps_filter: ", ls_scf_env%eps_filter
       WRITE(unit_nr,'()')
    ENDIF

    CALL report_matrix_sparsity(ls_scf_env%matrix_s,unit_nr,"overlap matrix (S)", &
                                ls_scf_env%eps_filter,error)

    DO ispin=1,nspin
       WRITE(title,'(A,I3)') "Kohn-Sham matrix (H) for spin ",ispin
       CALL report_matrix_sparsity(ls_scf_env%matrix_ks(ispin),unit_nr,title, &
                                   ls_scf_env%eps_filter,error)
    ENDDO

    CALL cp_dbcsr_init(matrix_tmp1,error=error)
    CALL cp_dbcsr_create(matrix_tmp1,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_tmp2,error=error)
    CALL cp_dbcsr_create(matrix_tmp2,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)

    DO ispin=1,nspin
       WRITE(title,'(A,I3)') "Density matrix (P) for spin ",ispin
       CALL report_matrix_sparsity(ls_scf_env%matrix_p(ispin),unit_nr,title, &
                                   ls_scf_env%eps_filter,error)

       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s, ls_scf_env%matrix_p(ispin), &
                              0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)

       WRITE(title,'(A,I3,A)') "S * P(",ispin,")"
       CALL report_matrix_sparsity(matrix_tmp1,unit_nr,title, ls_scf_env%eps_filter,error)

       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s, &
                              0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter, error=error)
       WRITE(title,'(A,I3,A)') "S * P(",ispin,") * S"
       CALL report_matrix_sparsity(matrix_tmp2,unit_nr,title, ls_scf_env%eps_filter,error)
    ENDDO

    IF (ls_scf_env%needs_s_inv) THEN
       CALL report_matrix_sparsity(ls_scf_env%matrix_s_inv,unit_nr,"inv(S)", &
                                   ls_scf_env%eps_filter,error)
       DO ispin=1,nspin
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_inv, ls_scf_env%matrix_ks(ispin), &
                                 0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)

          WRITE(title,'(A,I3,A)') "inv(S) * H(",ispin,")"
          CALL report_matrix_sparsity(matrix_tmp1,unit_nr,title, ls_scf_env%eps_filter,error)
       ENDDO
    ENDIF

    IF (ls_scf_env%use_s_sqrt) THEN

       CALL report_matrix_sparsity(ls_scf_env%matrix_s_sqrt,unit_nr,"sqrt(S)", &
                                   ls_scf_env%eps_filter,error)
       CALL report_matrix_sparsity(ls_scf_env%matrix_s_sqrt_inv,unit_nr,"inv(sqrt(S))", &
                                   ls_scf_env%eps_filter,error)

       DO ispin=1,nspin
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt_inv, ls_scf_env%matrix_ks(ispin), &
                                 0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s_sqrt_inv, &
                                 0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter,error=error)
          WRITE(title,'(A,I3,A)') "inv(sqrt(S)) * H(",ispin,") * inv(sqrt(S))"
          CALL report_matrix_sparsity(matrix_tmp2,unit_nr,title, ls_scf_env%eps_filter,error)
       ENDDO

       DO ispin=1,nspin
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt, ls_scf_env%matrix_p(ispin), &
                                 0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_s_sqrt, &
                                 0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter,error=error)
          WRITE(title,'(A,I3,A)') "sqrt(S) * P(",ispin,") * sqrt(S)"
          CALL report_matrix_sparsity(matrix_tmp2,unit_nr,title, ls_scf_env%eps_filter,error)
       ENDDO

    ENDIF

    CALL cp_dbcsr_release(matrix_tmp1,error=error)
    CALL cp_dbcsr_release(matrix_tmp2,error=error)

    CALL timestop(handle)

  END SUBROUTINE post_scf_sparsities

! *****************************************************************************
!> \brief Helper routine to report on the sparsity of a single matrix,
!>        for several filtering values
!> \param matrix ...
!> \param unit_nr ...
!> \param title ...
!> \param eps ...
!> \param error ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE report_matrix_sparsity(matrix,unit_nr,title,eps,error)
    TYPE(cp_dbcsr_type)                      :: matrix
    INTEGER                                  :: unit_nr
    CHARACTER(LEN=*)                         :: title
    REAL(KIND=dp)                            :: eps
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'report_matrix_sparsity', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: eps_local, occ
    TYPE(cp_dbcsr_type)                      :: matrix_tmp

    CALL timeset(routineN,handle)
    CALL cp_dbcsr_init(matrix_tmp,error=error)
    CALL cp_dbcsr_create(matrix_tmp,template=matrix,name=TRIM(title),error=error)
    CALL cp_dbcsr_copy(matrix_tmp,matrix,name=TRIM(title),error=error)

    IF (unit_nr>0) THEN
        WRITE(unit_nr,'(T2,A)') "Sparsity for : "//TRIM(title)
    ENDIF

    eps_local=eps
    DO
      IF (eps_local>1.1_dp) EXIT
      CALL cp_dbcsr_filter(matrix_tmp,eps_local,error=error)
      occ=cp_dbcsr_get_occupation(matrix_tmp)
      IF (unit_nr>0) WRITE(unit_nr,'(T2,F16.12,A3,F16.12)') eps_local," : ",occ
      eps_local=eps_local*10
    ENDDO

    CALL cp_dbcsr_release(matrix_tmp,error=error)

    CALL timestop(handle)

  END SUBROUTINE report_matrix_sparsity

! *****************************************************************************
!> \brief Compute matrix_w as needed for the forces
!> \param matrix_w ...
!> \param ls_scf_env ...
!> \param error ...
!> \par History
!>       2010.11 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE calculate_w_matrix(matrix_w,ls_scf_env,error)
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_w
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_w_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    REAL(KIND=dp)                            :: scaling
    TYPE(cp_dbcsr_type)                      :: matrix_tmp1, matrix_tmp2, &
                                                matrix_tmp3

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_init(matrix_tmp1,error=error)
    CALL cp_dbcsr_create(matrix_tmp1,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_tmp2,error=error)
    CALL cp_dbcsr_create(matrix_tmp2,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_tmp3,error=error)
    CALL cp_dbcsr_create(matrix_tmp3,template=ls_scf_env%matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)

    IF (ls_scf_env%nspins==1) THEN
      scaling=0.5_dp
    ELSE
      scaling=1.0_dp
    ENDIF


    DO ispin=1,ls_scf_env%nspins

      CALL cp_dbcsr_copy(matrix_tmp3,ls_scf_env%matrix_ks(ispin),error=error)
      IF (ls_scf_env%has_s_preconditioner) THEN
         CALL apply_matrix_preconditioner(matrix_tmp3,"backward", &
                       ls_scf_env%matrix_bs_sqrt,ls_scf_env%matrix_bs_sqrt_inv,error)
      ENDIF
      CALL cp_dbcsr_filter(matrix_tmp3,ls_scf_env%eps_filter,error=error)

      CALL cp_dbcsr_multiply("N", "N", scaling, ls_scf_env%matrix_p(ispin), matrix_tmp3, &
                              0.0_dp, matrix_tmp1, filter_eps=ls_scf_env%eps_filter,error=error)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ls_scf_env%matrix_p(ispin),&
                              0.0_dp, matrix_tmp2, filter_eps=ls_scf_env%eps_filter,error=error)
      CALL matrix_ls_to_qs(matrix_w(ispin)%matrix, matrix_tmp2, ls_scf_env%ls_mstruct, error=error)
    ENDDO

    CALL cp_dbcsr_release(matrix_tmp1,error=error)
    CALL cp_dbcsr_release(matrix_tmp2,error=error)
    CALL cp_dbcsr_release(matrix_tmp3,error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_w_matrix

! *****************************************************************************
!> \brief a place for quick experiments
!> \param ls_scf_env ...
!> \param error ...
!> \par History
!>       2010.11 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE post_scf_experiment(ls_scf_env,error)
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'post_scf_experiment', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, unit_nr
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    CALL timestop(handle)
  END SUBROUTINE post_scf_experiment

END MODULE dm_ls_scf
