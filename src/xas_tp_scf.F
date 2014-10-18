!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief xas_scf for the tp method
!>       It is repeaated for every atom that have to be excited
!> \par History
!>      created 05.2005
!> \author MI (05.2005)
! *****************************************************************************
MODULE xas_tp_scf
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_p_type
  USE cp_external_control,             ONLY: external_control
  USE cp_fm_types,                     ONLY: cp_fm_get_submatrix,&
                                             cp_fm_init_random,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: ot_precond_full_kinetic,&
                                             ot_precond_solver_default,&
                                             xas_dscf,&
                                             xas_scf_general
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE message_passing,                 ONLY: mp_sync
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE preconditioner,                  ONLY: make_preconditioner
  USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                             init_preconditioner,&
                                             preconditioner_type
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_density_mixing_types,         ONLY: broyden_mixing_new_nr,&
                                             broyden_mixing_nr,&
                                             direct_mixing_nr,&
                                             gspace_mixing_nr,&
                                             multisecant_mixing_nr,&
                                             no_mixing_nr,&
                                             pulay_mixing_nr
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_gspace_mixing,                ONLY: gspace_mixing,&
                                             self_consistency_check
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_loc_methods,                  ONLY: qs_loc_driver
  USE qs_loc_types,                    ONLY: get_qs_loc_env,&
                                             localized_wfn_control_type,&
                                             qs_loc_env_new_type
  USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             set_mo_occupation,&
                                             set_mo_set
  USE qs_ot_eigensolver,               ONLY: ot_eigensolver
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf,                          ONLY: scf_env_cleanup,&
                                             scf_env_do_scf
  USE qs_scf_diagonalization,          ONLY: do_general_diag
  USE qs_scf_initialization,           ONLY: qs_scf_env_initialize
  USE qs_scf_methods,                  ONLY: scf_env_density_mixing
  USE qs_scf_output,                   ONLY: qs_scf_print_summary
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE scf_control_types,               ONLY: scf_control_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_control,                     ONLY: xas_control_type
  USE xas_env_types,                   ONLY: get_xas_env,&
                                             set_xas_env,&
                                             xas_environment_type
  USE xas_restart,                     ONLY: xas_write_restart
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tp_scf'

  ! *** Public subroutines ***

  PUBLIC :: xas_do_tp_scf, xes_scf_once

CONTAINS

! *****************************************************************************
!> \brief perform an scf loop to calculate the xas spectrum
!>      given by the excitation of a inner state of a selected atom
!>      by using the transition potential method
!> \param dft_control ...
!> \param xas_env the environment for XAS  calculations
!> \param iatom ...
!> \param scf_env the scf_env where to perform the scf procedure
!> \param qs_env the qs_env, the scf_env and xas_env live in
!> \param xas_section ...
!> \param scf_section ...
!> \param converged ...
!> \param should_stop ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE xas_do_tp_scf (dft_control,xas_env,iatom,scf_env,qs_env,&
       xas_section,scf_section,converged,should_stop,error)

    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(xas_environment_type), POINTER      :: xas_env
    INTEGER, INTENT(IN)                      :: iatom
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: xas_section, scf_section
    LOGICAL, INTENT(OUT)                     :: converged, should_stop
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xas_do_tp_scf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, handle2, ispin, &
                                                iter_count, output_unit
    LOGICAL                                  :: diis_step, energy_only, &
                                                exit_loop, failure, gapw
    REAL(KIND=dp)                            :: t1, t2
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(xas_control_type), POINTER          :: xas_control

    CALL timeset(routineN,handle)
    NULLIFY(xas_control,matrix_s,matrix_ks,para_env)
    NULLIFY(rho,energy,scf_control,logger, ks_env,mos,atomic_kind_set)
    NULLIFY(qs_charges)

    logger => cp_error_get_logger(error)
    t1 = m_walltime()
    failure=.FALSE.
    converged = .TRUE.

    CPPrecondition(ASSOCIATED(xas_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(xas_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         scf_control=scf_control,&
         matrix_s=matrix_s,energy=energy,&
         qs_charges=qs_charges,&
         ks_env=ks_env,para_env=para_env,&
         error=error)

    energy_only = .FALSE.
    output_unit=cp_print_key_unit_nr(logger,xas_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".xasLog",error=error)
    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") "XAS_TP_SCF WAVEFUNCTION OPTIMIZATION"
    END IF

    !   GAPW method must be used
    gapw = dft_control%qs_control%gapw
    CPPrecondition(gapw,cp_failure_level,routineP,error,failure)
    xas_control => dft_control%xas_control

    CALL cp_add_iter_level(logger%iter_info,"XAS_SCF",error=error)

    CALL get_qs_env(qs_env,matrix_ks=matrix_ks,rho=rho,mos=mos,error=error)

    iter_count = 0
    diis_step = .FALSE.

    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,&
            FMT="(/,T3,A,T12,A,T31,A,T40,A,T60,A,T75,A/,T3,A)")&
            "Step","Update method","Time","Convergence","Total energy","Change",&
            REPEAT("-",77)
    END IF

    !   *** SCF loop ***

    energy%tot_old = 0.0_dp
    scf_loop: DO
       CALL timeset(routineN//"_inner_loop",handle2)

       exit_loop = .FALSE.
       IF (output_unit > 0) CALL m_flush(output_unit)

       iter_count = iter_count + 1
       CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=iter_count,error=error)

       ! ** here qs_env%rho%rho_r and qs_env%rho%rho_g should be up to date

       CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., just_energy=energy_only, error=error)

       scf_env%mixing_store%alpha =  xas_env%mixing_store%alpha
       scf_env%mixing_store%iter_method = xas_env%mixing_store%iter_method

       SELECT CASE (xas_control%scf_method)
       CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown scf method method for core level spectroscopy"//&
               cp_to_string(xas_control%scf_method),error,failure)

       CASE(xas_scf_general) ! diagonalisation (default)
          scf_env%iter_count=iter_count
          CALL do_general_diag(scf_env,mos,matrix_ks,&
                               matrix_s,scf_control,scf_section, &
                               diis_step,xas_env,error)

       END SELECT

       SELECT CASE(xas_env%mixing_method)
       CASE(direct_mixing_nr)
          CALL scf_env_density_mixing(scf_env%p_mix_new,&
               xas_env%mixing_store, rho%rho_ao, para_env, scf_env%iter_delta, scf_env%iter_count, &
               diis=diis_step, error=error)
       CASE(gspace_mixing_nr,pulay_mixing_nr,broyden_mixing_nr,&
                 broyden_mixing_new_nr,multisecant_mixing_nr)
          ! Compute the difference p_out-p_in
          CALL self_consistency_check(rho%rho_ao,scf_env%p_delta,para_env,scf_env%p_mix_new,&
                delta= scf_env%iter_delta, error=error)
       CASE(no_mixing_nr)
       CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown scf mixing method: "//&
               cp_to_string(xas_env%mixing_method),error,failure)
       END SELECT

       t2 = m_walltime()

       IF (output_unit>0.and.scf_env%print_iter_line) THEN
          WRITE (UNIT=output_unit,&
               FMT="(T2,I5,1X,A,T20,E8.2,1X,F6.1,1X,F14.8,1X,F20.10,1X,ES9.2)")&
               iter_count,TRIM(scf_env%iter_method),&
               scf_env%iter_param,t2 - t1,scf_env%iter_delta,energy%total,&
               energy%total-energy%tot_old

       END IF
       energy%tot_old = energy%total

       ! ** convergence check
       CALL external_control(should_stop,"XASSCF",target_time=qs_env%target_time,&
           start_time=qs_env%start_time,error=error)
       IF (scf_env%iter_delta < xas_control%eps_scf) THEN
          IF (output_unit>0) THEN
             WRITE(UNIT=output_unit,FMT="(/,T3,A,I5,A/)")&
                  "*** SCF run converged in ",iter_count," steps ***"
          END IF
          exit_loop = .TRUE.
       ELSE IF (should_stop.OR. iter_count == xas_control%max_scf) THEN
          IF (output_unit>0) THEN
             WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
                  "*** SCF run NOT converged ***"
          END IF
          converged = .FALSE.
          exit_loop = .TRUE.
       END IF
       !   *** Exit if we have finished with the SCF inner loop ***
       IF (exit_loop) THEN
          ! now, print out energies and charges corresponding to the obtained wfn
          ! (this actually is not 100% consistent at this point)!
          CALL qs_scf_print_summary(output_unit,qs_env,error=error)
          CALL cp_iterate(logger%iter_info,last=.TRUE.,iter_nr=iter_count,error=error)
       END IF

       ! ** Write restart file **
       CALL xas_write_restart(xas_env, xas_section, qs_env, xas_control%xas_method,&
            iatom,error=error)

       IF (exit_loop) THEN
          CALL timestop(handle2)
          EXIT scf_loop
       END IF

       IF (.NOT.BTEST(cp_print_key_should_output(logger%iter_info,&
         xas_section,"PRINT%ITERATION_INFO/TIME_CUMUL",error=error),cp_p_file)) t1 = m_walltime()

       !   *** mixing methods have the new density matrix in p_mix_new
       IF (xas_env%mixing_method > 0) THEN
          DO ispin=1,dft_control%nspins
             CALL cp_dbcsr_copy(rho%rho_ao(ispin)%matrix,scf_env%p_mix_new(ispin)%matrix,&
                  error=error)
          END DO
       ENDIF

       ! ** update qs_env%rho
       CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
       IF(xas_env%mixing_method>=gspace_mixing_nr) THEN
          CALL gspace_mixing(qs_env, scf_env, xas_env%mixing_store, rho, para_env, error=error)
       END IF

       CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)
       CALL timestop(handle2)

    END DO scf_loop

    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Ionization potential of the excited atom:      ",xas_env%IP_energy
       CALL m_flush(output_unit)
    END IF

    CALL mp_sync(para_env%group)
    CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)

    CALL cls_prepare_states(xas_control,xas_env,qs_env,iatom,xas_section,output_unit,error=error)

    CALL mp_sync(para_env%group)

    CALL cp_print_key_finished_output(output_unit,logger,xas_section,&
         "PRINT%PROGRAM_RUN_INFO", error=error)
    CALL cp_rm_iter_level(logger%iter_info,"XAS_SCF",error=error)

    CALL timestop(handle)

  END SUBROUTINE xas_do_tp_scf
! *****************************************************************************
!> \brief  Post processing of the optimized wfn in XAS scheme, as preparation for
!>         the calculation of the spectrum
!> \param xas_control ...
!> \param xas_env the environment for XAS  calculations
!> \param qs_env the qs_env, the scf_env and xas_env live in
!> \param iatom index of the excited atom
!> \param xas_section ...
!> \param output_unit ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE cls_prepare_states(xas_control,xas_env,qs_env,iatom,xas_section,output_unit,error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(xas_environment_type), POINTER      :: xas_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iatom
    TYPE(section_vals_type), POINTER         :: xas_section
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cls_prepare_states', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ikind, isgf, ispin, istat, istate, j, my_kind, &
      my_state, nao, natom, nexc_search, nmo, nvirtual2, uno_iter, xas_estate
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf
    INTEGER, DIMENSION(:), POINTER           :: mykind_of_kind
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:, :), POINTER       :: centers_wfn
    REAL(KIND=dp)                            :: component, dist, max_overlap, &
                                                ra(3), rac(3), rc(3), &
                                                sto_state_overlap, uno_eps
    REAL(KIND=dp), DIMENSION(:), POINTER     :: all_evals, eigenvalues, &
                                                uno_evals
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer, vecbuffer2
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: stogto_overlap
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: kinetic, matrix_ks, matrix_s
    TYPE(cp_fm_type), POINTER                :: all_vectors, excvec_coeff, &
                                                excvec_overlap, mo_coeff, &
                                                uno_orbs
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(preconditioner_type), POINTER       :: local_preconditioner
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(section_vals_type), POINTER         :: loc_section, print_loc_section

    CALL timeset(routineN,handle)

    NULLIFY(atomic_kind,dft_control,matrix_s,matrix_ks, qs_kind_set, kinetic)
    NULLIFY(cell,particle_set,local_preconditioner,vecbuffer,vecbuffer2)
    NULLIFY(dft_control,loc_section,mos,mo_coeff,eigenvalues)
    NULLIFY(centers_wfn,mykind_of_kind,qs_loc_env,localized_wfn_control,stogto_overlap)
    NULLIFY(all_evals,all_vectors,excvec_coeff,excvec_overlap,uno_evals,para_env,blacs_env)

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(xas_env),cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env,&
                    cell=cell,&
                    dft_control=dft_control,&
                    matrix_ks=matrix_ks,&
                    matrix_s=matrix_s,&
                    kinetic=kinetic,&
                    mos=mos,&
                    particle_set=particle_set,&
                    qs_kind_set=qs_kind_set,&
                    para_env=para_env,&
                    blacs_env=blacs_env,&
                    error=error)

    ! Some elements from the xas_env
    CALL get_xas_env(xas_env=xas_env,&
         all_vectors=all_vectors,all_evals=all_evals,&
         excvec_coeff=excvec_coeff,&
         nvirtual2=nvirtual2,xas_estate=xas_estate,&
         excvec_overlap=excvec_overlap,nexc_search=nexc_search,error=error)
    CPPrecondition(ASSOCIATED(excvec_overlap),cp_failure_level,routineP,error,failure)

    CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff,nao=nao,&
         eigenvalues=eigenvalues)

    ALLOCATE(vecbuffer(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    vecbuffer = 0.0_dp
    ALLOCATE(vecbuffer2(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    vecbuffer2 = 0.0_dp
    natom=SIZE(particle_set,1)
    ALLOCATE (first_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, error=error)
    ALLOCATE(centers_wfn(3,nexc_search),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    centers_wfn=0.0_dp

    ! Possible only for emission only
    IF(.NOT.xas_control%scf_method==xas_scf_general) THEN
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") " Eigenstates are derived "//&
          "from the MOs optimized by OT. Follows localization of the core states"//&
          " to identify the excited orital. "
       END IF

       CALL get_xas_env(xas_env=xas_env, &
            mykind_of_kind=mykind_of_kind, qs_loc_env=qs_loc_env,&
            stogto_overlap=stogto_overlap,error=error)
       CALL get_qs_loc_env(qs_loc_env=qs_loc_env,&
            localized_wfn_control=localized_wfn_control,error=error)
       loc_section => section_vals_get_subs_vals(xas_section,"LOCALIZE",error=error)
       print_loc_section => section_vals_get_subs_vals(loc_section,"PRINT",error=error)
       CALL qs_loc_driver(qs_env,qs_loc_env,loc_section,print_loc_section,myspin=1,error=error)
       ra(1:3) = particle_set(iatom)%r(1:3)

       NULLIFY(atomic_kind)
       atomic_kind =>  particle_set(iatom)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            kind_number=ikind)
       my_kind = mykind_of_kind(ikind)


       CALL cp_fm_get_submatrix(mo_coeff,vecbuffer2,1,my_state,&
            nao,1,transpose=.TRUE.,error=error)

       ! Rotate the wfn to get the eigenstate of the KS hamiltonian
       ! Only ispin=1 should be needed
       DO ispin=1,dft_control%nspins
          CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff,nmo=nmo,&
               eigenvalues=eigenvalues)
          CALL calculate_subspace_eigenvalues(mo_coeff,&
               matrix_ks(ispin)%matrix,eigenvalues, &
               do_rotation=.TRUE.,error=error)
       END DO  ! ispin

       !Search for the core state to be excited
       max_overlap = 0.0_dp
       DO istate = 1,nexc_search
          centers_wfn(1,istate) =  localized_wfn_control%centers_set(1)%array(1,istate)
          centers_wfn(2,istate) =  localized_wfn_control%centers_set(1)%array(2,istate)
          centers_wfn(3,istate) =  localized_wfn_control%centers_set(1)%array(3,istate)

          rc(1:3) = centers_wfn(1:3,istate)
          rac = pbc(ra,rc,cell)
          dist = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)

          IF(dist < 1.0_dp) THEN
             CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,istate,&
                  nao,1,transpose=.TRUE.,error=error)
             sto_state_overlap=0.0_dp
             DO i = 1,SIZE(stogto_overlap(my_kind)%array,1)
                component = 0.0_dp
                DO j = 1,SIZE(stogto_overlap(my_kind)%array,2)
                   isgf = first_sgf(iatom) + j - 1
                   component = component + stogto_overlap(my_kind)%array(i,j)*vecbuffer(1,isgf)
                END DO  ! j size
                sto_state_overlap = sto_state_overlap + &
                     component * component
             END DO  ! i size
             IF(sto_state_overlap .GT. max_overlap) THEN
                max_overlap = sto_state_overlap
                my_state = istate
             END IF
          END IF
          xas_estate = my_state
       END DO  !  istate

       CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff)
       CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,xas_estate,&
            nao,1,transpose=.TRUE.,error=error)
       CALL cp_fm_set_submatrix(excvec_coeff,vecbuffer2,1,1,&
            nao,1,transpose=.TRUE.,error=error)
       !
    END IF

    CALL mp_sync(para_env%group)
    !Calculate the virtual states from the KS matrix matrix_ks(1)
    IF(nvirtual2 .GT. 0) THEN
       NULLIFY(mo_coeff)
       CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff,nmo=nmo)
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,/,T2,A,I5,A,I6,A)") " Calculation of ", nvirtual2,&
            " additional virtual states of the subspace complementary to the "//&
            " lowest ", nmo, " states"
       END IF

       NULLIFY(uno_orbs,uno_evals,local_preconditioner)
       ALLOCATE(local_preconditioner,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL init_preconditioner(local_preconditioner,para_env=para_env,&
                                blacs_env=blacs_env,error=error)

       CALL make_preconditioner(local_preconditioner, &
            precon_type=ot_precond_full_kinetic, &
            solver_type=ot_precond_solver_default, &
            matrix_h=matrix_ks(1)%matrix, &
            matrix_s=matrix_s(1)%matrix,&
            matrix_t=kinetic(1)%matrix, &
            convert_precond_to_dbcsr=.TRUE.,&
            mo_set=mos(1)%mo_set,energy_gap=0.2_dp,error=error)

       CALL get_xas_env(xas_env=xas_env,unoccupied_orbs=uno_orbs,&
            unoccupied_evals=uno_evals,unoccupied_eps=uno_eps,unoccupied_max_iter=uno_iter,error=error)
       CALL cp_fm_init_random(uno_orbs,nvirtual2,error=error)

       CALL ot_eigensolver(matrix_h=matrix_ks(1)%matrix,matrix_s=matrix_s(1)%matrix, &
            matrix_c_fm=uno_orbs,matrix_orthogonal_space_fm=mo_coeff,&
            preconditioner=local_preconditioner,eps_gradient=uno_eps,&
            iter_max=uno_iter,size_ortho_space=nmo,error=error)

       CALL calculate_subspace_eigenvalues(uno_orbs,matrix_ks(1)%matrix,&
            uno_evals,do_rotation=.TRUE.,error=error)
       CALL destroy_preconditioner(local_preconditioner,error=error)

       DEALLOCATE(local_preconditioner,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL mp_sync(para_env%group)
    ! Prapare arrays for the calculation of the spectrum
    IF(.NOT. xas_control%xas_method == xas_dscf) THEN
       ! Copy the final vectors in the array
       NULLIFY(all_vectors,all_evals)
       CALL get_xas_env(xas_env=xas_env,all_vectors=all_vectors,&
            all_evals=all_evals,error=error)
       CALL get_mo_set(mos(1)%mo_set, eigenvalues=eigenvalues,mo_coeff=mo_coeff,&
            nmo=nmo)

       CALL cp_fm_to_fm(mo_coeff,all_vectors,ncol=nmo,&
               source_start=1,target_start=1)
       DO istate = 1,nmo
             all_evals(istate) = eigenvalues(istate)
       ENDDO
       IF(nvirtual2 .GT. 0) THEN
          CALL cp_fm_to_fm(uno_orbs,all_vectors,ncol=nvirtual2,&
               source_start=1,target_start=1+nmo)
         DO istate = 1,nvirtual2
           all_evals(istate+nmo) = uno_evals(istate)
         END DO
       END IF
    END IF

    DEALLOCATE(vecbuffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(vecbuffer2,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(centers_wfn,first_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE cls_prepare_states

! *****************************************************************************
!> \brief  SCF for emission spectra calculations: vacancy in valence
!> \param qs_env the qs_env, the scf_env and xas_env live in
!> \param xas_env the environment for XAS  calculations
!> \param scf_section ...
!> \param converged ...
!> \param should_stop ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE xes_scf_once(qs_env,xas_env,scf_section,converged,should_stop,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(xas_environment_type), POINTER      :: xas_env
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(OUT)                     :: converged, should_stop
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xes_scf_once', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, istate, &
                                                nelectron, nmo, nvirtual, &
                                                nvirtual2, output_unit
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: all_evals, eigenvalues
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_fm_type), POINTER                :: all_vectors, mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control

    failure = .FALSE.
    NULLIFY(dft_control,scf_control,matrix_ks,mos,para_env)
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    CALL timeset(routineN,handle)

    CPPrecondition(ASSOCIATED(xas_env),cp_failure_level,routineP,error,failure)
    IF(.NOT. failure) THEN
       CALL get_qs_env(qs_env=qs_env,dft_control=dft_control, scf_control=scf_control,&
            matrix_ks=matrix_ks,mos=mos,para_env=para_env,error=error)


       IF(xas_env%homo_occ == 0) THEN
         CALL get_mo_set(mos(1)%mo_set,nelectron=nelectron)
         nelectron = nelectron - 1
         CALL set_mo_set(mos(1)%mo_set,nelectron=nelectron,error=error)
         CALL set_mo_occupation(mo_set=mos(1)%mo_set,smear=xas_env%smear,&
                             error=error)


         NULLIFY(scf_env)
         CALL get_qs_env(qs_env,scf_env=scf_env,error=error)
         CALL qs_scf_env_initialize(scf_env=scf_env,qs_env=qs_env,&
                                  do_xas_tp=.TRUE., error=error)

         CALL scf_env_do_scf(scf_env=scf_env, qs_env=qs_env,  &
              converged=converged, should_stop=should_stop, error=error)
         CALL scf_env_cleanup(scf_env,qs_env=qs_env,error=error)

       END IF

       !   The eigenstate of the KS Hamiltonian are nedeed
       NULLIFY(mo_coeff,eigenvalues)
       IF(scf_control%use_ot) THEN
          IF (output_unit>0) THEN
             WRITE (UNIT=output_unit,FMT="(/,T10,A,/)")&
                  "Get eigenstates and eigenvalues from ground state MOs"
          END IF
          DO ispin = 1,SIZE(mos)
             CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff,&
                  eigenvalues=eigenvalues)
             CALL calculate_subspace_eigenvalues(mo_coeff,&
                  matrix_ks(ispin)%matrix,eigenvalues, &
                  do_rotation=.TRUE.,error=error)
          END DO
       END IF
       NULLIFY(all_vectors,all_evals)
       CALL get_xas_env(xas_env=xas_env,all_vectors=all_vectors,&
            all_evals=all_evals,nvirtual2=nvirtual2,error=error)
       CALL get_mo_set(mos(1)%mo_set, eigenvalues=eigenvalues,mo_coeff=mo_coeff, nmo=nmo)

       CALL cp_fm_to_fm(mo_coeff,all_vectors,ncol=nmo,&
            source_start=1,target_start=1)
       DO istate = 1,nmo
          all_evals(istate) = eigenvalues(istate)
       ENDDO

       IF(nvirtual2/=0) THEN
          IF (output_unit>0) THEN
             WRITE (UNIT=output_unit,FMT="(/,T10,A,/)")&
                  "WARNING: for this XES calculation additional unoccupied MOs are not needed"
          END IF
          nvirtual2=0
          nvirtual = nmo
          CALL set_xas_env(xas_env=xas_env,nvirtual=nvirtual,nvirtual2=nvirtual2,error=error)
       END IF

       IF(xas_env%homo_occ == 0) THEN
          nelectron = nelectron + 1
          CALL set_mo_set(mos(1)%mo_set,nelectron=nelectron,error=error)
       END IF
    END IF  ! failure
    CALL timestop(handle)

  END SUBROUTINE xes_scf_once

END MODULE xas_tp_scf
