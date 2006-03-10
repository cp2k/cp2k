!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2006  CP2K developers group
!-----------------------------------------------------------------------------!

!!****** cp2k/qs_scf [1.0] *
!!
!!   NAME
!!     qs_scf
!!
!!   FUNCTION
!!     Routines for the Quickstep SCF run.
!!
!!   AUTHOR
!!     Matthias Krack (30.04.2001)
!!
!!   MODIFICATION HISTORY
!!     - Joost VandeVondele (02.2002)
!!          added code for: incremental (pab and gvg) update
!!                           initialisation (init_cube, l_info)
!!     - Joost VandeVondele (02.2002)
!!          called the poisson code of the classical part
!!          this takes into account the spherical cutoff and allows for
!!          isolated systems
!!     - Joost VandeVondele (02.2002)
!!          added multiple grid feature
!!          changed to spherical cutoff consistently (?)
!!          therefore removed the gradient correct functionals
!!     - updated with the new QS data structures (10.04.02,MK)
!!     - copy_matrix replaced by transfer_matrix (11.04.02,MK)
!!     - nrebuild_rho and nrebuild_gvg unified (12.04.02,MK)
!!     - set_mo_occupation for smearing of the MO occupation numbers
!!       (17.04.02,MK)
!!     - MO level shifting added (22.04.02,MK)
!!     - Usage of TYPE mo_set_p_type
!!     - Joost VandeVondele (05.2002)
!!           added cholesky based diagonalisation
!!     - 05.2002 added pao method [fawzi]
!!     - parallel FFT (JGH 22.05.2002)
!!     - 06.2002 moved KS matrix construction to qs_build_KS_matrix.F [fawzi]
!!     - started to include more LSD (01.2003,Joost VandeVondele)
!!     - 02.2003 scf_env [fawzi]
!!     - got rid of nrebuild (01.2004, Joost VandeVondele)
!!     - 10.2004 removed pao [fawzi]
!!     - 03.2006 large cleaning action [Joost VandeVondele]
!!
!!   SOURCE
!******************************************************************************

MODULE qs_scf
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
  USE cp_fm_diag,                      ONLY: cp_fm_power
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_get_el_struct
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm
  USE cp_sm_pool_types,                ONLY: cp_sm_pool_p_type,&
                                             sm_pools_create_matrix_vect,&
                                             sm_pools_flush_cache,&
                                             sm_pools_give_back_matrix_vect
  USE global_types,                    ONLY: global_environment_type
  USE harris_env_types,                ONLY: harris_env_type
  USE harris_functional,               ONLY: harris_eigenvalue_calculation,&
                                             harris_eigenvalue_trace_KS_Pmix,&
                                             harris_energy_correction,&
                                             harris_postprocessing
  USE input_constants,                 ONLY: core_guess,&
                                             densities_guess,&
                                             ot_precond_full_all,&
                                             ot_precond_none,&
                                             qs_debug_run
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: evolt,&
                                             kcalmol
  USE preconditioner,                  ONLY: destroy_preconditioner,&
                                             init_preconditioner,&
                                             make_preconditioner
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_diis,                         ONLY: qs_diis_b_clear,&
                                             qs_diis_b_create,&
                                             qs_diis_b_step
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_reorthogonalize_vectors,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_initial_guess,                ONLY: calculate_first_density_matrix
  USE qs_ks_methods,                   ONLY: qs_ks_create,&
                                             qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release
  USE qs_loc_methods,                  ONLY: loc_initialize_rho
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                             calculate_magnitude,&
                                             calculate_orthonormality,&
                                             calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             init_mo_set,&
                                             mo_set_p_type,&
                                             set_mo_set,&
                                             write_mo_set
  USE qs_ot,                           ONLY: qs_ot_new_preconditioner
  USE qs_ot_scf,                       ONLY: ot_scf_destroy,&
                                             ot_scf_init,&
                                             ot_scf_mini,&
                                             ot_scf_read_input
  USE qs_overlap,                      ONLY: write_fm_with_basis_info,&
                                             write_sparse_matrix
  USE qs_rho_atom_methods,             ONLY: calculate_rho_atom_coeff
  USE qs_rho_methods,                  ONLY: duplicate_rho_type,&
                                             qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_methods,                  ONLY: eigensolver,&
                                             scf_env_density_mixing,&
                                             simple_eigensolver
  USE qs_scf_post,                     ONLY: scf_post_calculation
  USE qs_scf_types,                    ONLY: ao_mixing_nr,&
                                             general_diag_method_nr,&
                                             no_mixing_nr,&
                                             ot_method_nr,&
                                             qs_scf_env_type,&
                                             scf_env_create,&
                                             scf_env_release,&
                                             special_diag_method_nr
  USE qs_wf_history_methods,           ONLY: wfi_extrapolate,&
                                             wfi_get_method_label,&
                                             wfi_update
  USE qs_wf_history_types,             ONLY: wfi_frozen_method_nr,&
                                             wfi_use_guess_method_nr
  USE scf_control_types,               ONLY: scf_control_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type,&
                                             real_matrix_type,&
                                             transfer_matrix
  USE termination,                     ONLY: external_control
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf'

! *** Public subroutines ***

  PUBLIC :: scf, scf_env_cleanup, &
            init_scf_run, init_scf_loop

! *****************************************************************************

CONTAINS

!!****f* qs_scf/scf [1.0] *
!!
!!   NAME
!!     scf
!!
!!   SYNOPSIS
!!     Subroutine scf(qs_env, globenv, error)
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Type(global_environment_type), pointer:: globenv
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine scf
!!
!!   FUNCTION
!!     perform an scf procedure in the given qs_env
!!
!!   NOTES
!!     globenv should be removed
!!
!!   ARGUMENTS
!!     - qs_env: the qs_environment where to perform the scf procedure
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 introduced scf_env, moved real work to scf_env_do_scf [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf(qs_env,globenv,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: converged, failure, &
                                                should_stop
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: dft_section, input, &
                                                scf_section

    NULLIFY(scf_env)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL get_qs_env(qs_env,scf_env=scf_env,error=error,input=input, &
                       dft_control=dft_control, scf_control=scf_control)

       dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
       scf_section => section_vals_get_subs_vals(dft_section,"SCF",error=error)

       IF (.NOT.ASSOCIATED(scf_env)) THEN ! i.e. for MD this is associated on the second step (it so seems)
          CALL scf_env_create(scf_env, error=error)
          CALL set_qs_env(qs_env,scf_env=scf_env,error=error)
          CALL scf_env_release(scf_env)

          ! set some of the methods that might be used in this SCF.
          ! this might not yet be the ideal place to set this kind of stuff (who knows?)
          CALL get_qs_env(qs_env,scf_env=scf_env,scf_control=scf_control, &
                          dft_control=dft_control, error=error)

          scf_env%method=general_diag_method_nr 
          IF (dft_control%qs_control%semi_empirical) THEN
            IF (dft_control%qs_control%se_control%orthogonal_basis) &
              scf_env%method=special_diag_method_nr
          END IF
          IF (scf_control%use_ot) scf_env%method=ot_method_nr
  
          SELECT CASE(scf_env%method)
          CASE(ot_method_nr)
             scf_env%mixing_method=no_mixing_nr
          CASE DEFAULT
             scf_env%mixing_method=ao_mixing_nr
          END SELECT
          scf_env%p_mix=scf_control%p_mix

       END IF

       CALL init_scf_run(scf_env=scf_env, qs_env=qs_env, globenv=globenv, &
                                      scf_section=scf_section, error=error)

       CALL scf_env_do_scf(scf_env=scf_env, qs_env=qs_env, globenv=globenv, &
                           converged=converged, should_stop=should_stop, error=error)

       !   *** add the converged wavefunction to the wavefunction history
       IF (ASSOCIATED(qs_env%wf_history)) THEN
           CALL wfi_update(qs_env%wf_history,qs_env=qs_env,dt=1.0_dp, error=error)
       END IF
       
       ! *** compute properties that depend on the converged wavefunction
       IF ( .NOT. ( should_stop .OR. dft_control%qs_control%semi_empirical ) ) THEN
              CALL scf_post_calculation(dft_section, scf_env, qs_env, globenv, error)
       END IF
       
       ! *** cleanup
       CALL scf_env_cleanup(scf_env,qs_env=qs_env,error=error)

    END IF

  END SUBROUTINE scf
!***************************************************************************

!!****f* qs_scf/scf_env_do_scf [1.0] *
!!
!!   NAME
!!     scf_env_do_scf
!!
!!   SYNOPSIS
!!     Subroutine scf_env_do_scf(scf_env, qs_env, globenv, error)
!!       Type(qs_scf_env_type), Pointer:: scf_env
!!       Type(global_environment_type), pointer:: globenv
!!       Type(qs_environment_type), Pointer:: qs_env
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine scf_env_do_scf
!!
!!   FUNCTION
!!     perform an scf loop
!!
!!   NOTES
!!     globenv should be removed
!!
!!   ARGUMENTS
!!     - scf_env: the scf_env where to perform the scf procedure
!!     - qs_env: the qs_env, the scf_env lives in
!!     - converged: will be true / false if converged is reached
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     long history, see cvs and qs_scf module history
!!     02.2003 introduced scf_env [fawzi]
!!     09.2005 Frozen density approximation [TdK]
!!
!!*** **********************************************************************
  SUBROUTINE scf_env_do_scf(scf_env,qs_env,globenv,converged,should_stop,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(OUT)                     :: converged, should_stop
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scf_env_do_scf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, handle2, i, ispin, &
                                                output_unit
    LOGICAL                                  :: diis_step, energy_only, &
                                                failure, gapw, gapw_xc, &
                                                harris_flag, has_unit_metric, &
                                                use_jacobi
    REAL(KIND=dp)                            :: t1, t2
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho, rho_xc
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: dft_section, input, &
                                                scf_section

    CALL timeset(routineN,"I"," ",handle)

    failure=.FALSE.
    converged=.TRUE.

    NULLIFY(dft_control,matrix_h,matrix_s,matrix_ks,rho,rho_xc,energy, &
            scf_control,logger,qs_charges,ks_env,mos,atomic_kind_set, &
            particle_set,harris_env,dft_section,input,&
            scf_section)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)

    logger => cp_error_get_logger(error)
    t1 = m_walltime()

    CALL get_qs_env(qs_env=qs_env,input=input, dft_control=dft_control, &
                    scf_control=scf_control,error=error)

    dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
    scf_section => section_vals_get_subs_vals(dft_section,"SCF",error=error)

    output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".scfLog",error=error)

    IF (output_unit>0) WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") &
            "SCF WAVEFUNCTION OPTIMIZATION"

    CALL init_scf_loop(scf_env=scf_env, qs_env=qs_env, globenv=globenv, &
                        scf_section=scf_section, error=error)

    ! short cut flags setting the different methods
    gapw = dft_control%qs_control%gapw
    gapw_xc = dft_control%qs_control%gapw_xc
    harris_flag = qs_env%use_harris
    has_unit_metric=.FALSE.
    IF (dft_control%qs_control%semi_empirical) &
       has_unit_metric=dft_control%qs_control%se_control%orthogonal_basis

    CALL get_qs_env(qs_env=qs_env,&
                    matrix_h=matrix_h,&
                    matrix_s=matrix_s,energy=energy,&
                    particle_set=particle_set,&
                    qs_charges=qs_charges,&
                    ks_env=ks_env, &
                    harris_env=harris_env,&
                    atomic_kind_set=atomic_kind_set,&
                    matrix_ks=matrix_ks,rho=rho,rho_xc=rho_xc,mos=mos, &
                    error=error)

    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,&
            FMT="(/,T3,A,T9,A,T36,A,T49,A,T68,A,/,T3,A)")&
            "Step","Update method","Time","Convergence","Total energy",&
            REPEAT("-",77)
    END IF

!   Some flags needed to be set at the beginning of the loop
    diis_step = .FALSE.
    use_jacobi = .FALSE.
    energy_only = .FALSE.

!   *** SCF loop, optimisation of the wfn coefficients
!   *** qs_env%rho%rho_r and qs_env%rho%rho_g should be up to date here
    scf_env%iter_count = 0

    scf_loop: DO

       CALL timeset("scf_iter","I"," ",handle2)
       scf_env%iter_count = scf_env%iter_count + 1

       IF (output_unit > 0) CALL m_flush(output_unit)

       CALL qs_ks_update_qs_env(ks_env,qs_env=qs_env,&
            calculate_forces=.FALSE.,just_energy=energy_only,error=error)

       ! print 'heavy weight' or relatively expensive quantities
       CALL qs_scf_loop_print(qs_env,scf_env,globenv,error)

       scf_env%p_mix = scf_control%p_mix
       scf_env%iter_param = 0.0_dp
       IF (scf_env%mixing_method==ao_mixing_nr) scf_env%iter_param = scf_env%p_mix

       ! this takes known energy and derivatives and produces new wfns and or density matrix
       SELECT CASE (scf_env%method)
       CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                 routineP,"unknown scf method: "//&
                 cp_to_string(scf_env%method),error,failure)
       CASE(general_diag_method_nr) 
          CALL qs_scf_loop_do_general_diag(scf_env,mos,rho,matrix_ks,matrix_s,scf_control,scf_section, &
                                       diis_step,use_jacobi,error)
       CASE(special_diag_method_nr)
          CALL qs_scf_loop_do_special_diag(scf_env,mos,rho,matrix_ks,scf_control,scf_section, &
                                       diis_step,error)
       CASE(ot_method_nr) 
          CALL qs_scf_loop_do_ot(scf_env,mos,rho,qs_env%mo_derivs,energy%total, & 
                                 matrix_s, energy_only, has_unit_metric,error)
       END SELECT

       ! another heavy weight print object, print controlled by dft_section
       DO ispin=1,SIZE(mos)
          CALL write_mo_set(mos(ispin)%mo_set,atomic_kind_set,&
                            particle_set,4,6,dft_section,globenv)
       ENDDO

       ! ** calculation of the harris energy correction *** !
       IF (harris_flag) THEN
         CALL harris_energy_correction(qs_env, harris_env, &
                                       globenv=globenv, fast=.TRUE., error=error)
         IF (scf_env%method .NE. ot_method_nr) THEN
           CALL harris_eigenvalue_trace_KS_Pmix(scf_env, qs_env, harris_env, error=error)
         ELSE
           CALL harris_eigenvalue_calculation(qs_env=qs_env, harris_env=harris_env, error=error)
         END IF
         IF (globenv%run_type_id == qs_debug_run) THEN
           energy%total = harris_env%harris_energy%Eharris
         END IF
       END IF

       SELECT CASE(scf_env%mixing_method)
       CASE(ao_mixing_nr)
          CALL scf_env_density_mixing(scf_env%p_mix_new,scf_env%mixing_method,&
                         scf_env%p_mix,scf_env%iter_delta, qs_env=qs_env,error=error)
       CASE(no_mixing_nr)
       CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                 routineP,"unknown scf mixing method: "//&
                 cp_to_string(scf_env%mixing_method),error,failure)
       END SELECT

       t2 = m_walltime()

       IF ((output_unit>0).and.scf_env%print_iter_line) THEN
          WRITE (UNIT=output_unit,&
               FMT="(T2,I5,2X,A,T22,E10.2,T32,F8.2,T40,2F20.10)")&
               scf_env%iter_count,TRIM(scf_env%iter_method),&
               scf_env%iter_param,t2 - t1,scf_env%iter_delta,energy%total
       END IF

       IF (harris_flag) THEN
          CALL harris_postprocessing(harris_env, error=error)
       END IF

       ! ** can we exit the SCF loop  ?
       CALL external_control(should_stop,"SCF",globenv,error)
       IF (scf_env%iter_delta < scf_control%eps_scf) THEN
          IF (output_unit>0) THEN
             WRITE(UNIT=output_unit,FMT="(/,T3,A,I5,A/)")&
                  "*** SCF run converged in ",scf_env%iter_count," steps ***"
          END IF
          CALL timestop(0.0_dp,handle2)
          EXIT scf_loop
       ELSE IF (should_stop.OR.  scf_env%iter_count == scf_control%max_scf) THEN
          IF (output_unit>0) THEN
             WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
                  "*** SCF run NOT converged ***"
          END IF
          converged=.FALSE.
          CALL timestop(0.0_dp,handle2)
          EXIT scf_loop
       END IF

       ! ** Write restart file **
       IF (scf_control%write_restart_each>0.AND.&
            MODULO(scf_env%iter_count,scf_control%write_restart_each)==0) THEN
          CALL write_mo_set(mos,atomic_kind_set,particle_set,globenv,&
               id_nr=1)
       END IF

       IF (.NOT.BTEST(cp_print_key_should_output(logger%iter_info,&
            scf_section,"PRINT%ITERATION_INFO/TIME_CUMUL",error=error),cp_p_file)) t1 = m_walltime()  

       ! ** mixing methods have the new density matrix in p_mix_new **
       IF (scf_env%mixing_method.EQ.ao_mixing_nr) THEN
           DO ispin=1,dft_control%nspins
              CALL transfer_matrix(scf_env%p_mix_new(ispin)%matrix,&
                                       rho%rho_ao(ispin)%matrix)
           END DO
       ENDIF

       ! ** update qs_env%rho
       CALL qs_rho_update_rho(rho, qs_env=qs_env, gapw=gapw, error=error)
       IF(gapw ) THEN
         CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
       ENDIF
       IF(gapw_xc ) THEN
         CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
         CALL qs_rho_update_rho(rho_xc, qs_env=qs_env, gapw=gapw_xc, error=error)
       ENDIF

       CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)

       CALL timestop(0.0_dp,handle2)

    END DO scf_loop

!   *** now, print out energies and charges corresponding to the obtained wfn
!   (this actually seems not 100% consistent at this point)
    CALL qs_scf_print_summary(output_unit,rho,qs_charges,energy,scf_env, &
                              dft_control,qs_env%qmmm,gapw,gapw_xc)

    IF (harris_flag) THEN
      energy%total = harris_env%harris_energy%Eharris
    END IF

!   *** Write restart file ***
    IF (scf_control%write_restart_each>=0) THEN
       CALL write_mo_set(mos,atomic_kind_set,particle_set,globenv,&
            id_nr=1)
    END IF

!   *** mixing methods need to undo mixing of the density matrix (restore original density) ***
    IF (scf_env%mixing_method.EQ.ao_mixing_nr) THEN
        CALL scf_env_density_mixing(scf_env%p_mix_new,scf_env%mixing_method,&
                                    scf_env%p_mix,scf_env%iter_delta,&
                                    qs_env=qs_env,invert=.TRUE.,error=error)
        DO ispin=1,dft_control%nspins
           CALL transfer_matrix(scf_env%p_mix_new(ispin)%matrix,&
                              rho%rho_ao(ispin)%matrix)
        END DO
    ENDIF

!   *** update rspace rho since the mo changed
!   *** this might not always be needed (i.e. no post calculation / no forces )
!   *** but guarantees that rho and wfn are consistent at this point
    CALL qs_rho_update_rho(rho, qs_env=qs_env, gapw = gapw, error=error)
    IF(gapw) THEN
       CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
    END IF
    IF(gapw_xc ) THEN
       CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
       CALL qs_rho_update_rho(rho_xc, qs_env=qs_env, gapw=gapw_xc, error=error)
    ENDIF

    CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)


    CALL timestop(0.0_dp,handle)

  END SUBROUTINE scf_env_do_scf

!!****f* qs_scf/qs_scf_loop_do_general_diag *
!!
!!   NAME
!!      qs_scf_loop_do_general_diag
!!
!!   FUNCTION
!!      the inner loop of scf, specific to diagonalization with S matrix
!!      basically, in goes the ks matrix out goes a new p matrix
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE qs_scf_loop_do_general_diag(scf_env,mos,rho,matrix_ks,matrix_s,scf_control,scf_section, &
                                       diis_step,use_jacobi,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step, use_jacobi
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ispin
    LOGICAL                                  :: do_level_shift
    REAL(kind=dp)                            :: diis_error

   DO ispin=1,SIZE(matrix_ks)
    CALL copy_sm_to_fm(matrix_ks(ispin)%matrix,scf_env%scf_work1(ispin)%matrix)
   ENDDO

   IF (scf_env%iter_count > 1) THEN
      CALL qs_diis_b_step(scf_env%scf_diis_buffer,mos,scf_env%scf_work1,&
           scf_env%scf_work2,scf_env%iter_delta,diis_error,diis_step,&
           scf_control%eps_diis,s_matrix=matrix_s,scf_section=scf_section,&
           error=error)
   END IF

   do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
      ((scf_control%density_guess == core_guess).OR.(scf_env%iter_count > 1)))

   IF (diis_step) THEN
      scf_env%p_mix = 1.0_dp
      scf_env%iter_param = diis_error
      IF (use_jacobi) THEN
         scf_env%iter_method = "DIIS/Jacobi"
      ELSE
         scf_env%iter_method = "DIIS/Diag."
      END IF
   ELSE
      IF (use_jacobi) THEN
         scf_env%iter_method = "Mixing/Jacobi"
      ELSE
         scf_env%iter_method = "Mixing/Diag."
      END IF
   END IF

   IF ((scf_env%iter_count > 1).AND.(scf_env%iter_delta < scf_control%eps_jacobi)) THEN
      use_jacobi = .TRUE.
   ELSE
      use_jacobi = .FALSE.
   END IF

   scf_env%iter_delta=0.0_dp

   DO ispin=1,SIZE(matrix_ks)

      CALL eigensolver(scf_env%scf_work1(ispin)%matrix, mos(ispin)%mo_set,scf_env%ortho,&
           scf_env%scf_work2,&
           do_level_shift,scf_control%level_shift,&
           use_cholesky=scf_control%use_cholesky,&
           work_syevx=scf_control%work_syevx,&
           use_jacobi=use_jacobi,&
           jacobi_threshold=scf_control%jacobi_threshold,&
           smear=scf_control%smear, &
           error=error)

      CALL calculate_density_matrix(mos(ispin)%mo_set,&
           scf_env%p_mix_new(ispin)%matrix,error=error)

   ENDDO

END SUBROUTINE qs_scf_loop_do_general_diag

!!****f* qs_scf/qs_scf_loop_do_special_diag *
!!
!!   NAME
!!      qs_scf_loop_do_special_diag
!!
!!   FUNCTION
!!      the inner loop of scf, specific to diagonalization without S matrix
!!      basically, in goes the ks matrix out goes a new p matrix
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE qs_scf_loop_do_special_diag(scf_env,mos,rho,matrix_ks,scf_control,scf_section, &
                                       diis_step,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ispin
    LOGICAL                                  :: do_level_shift
    REAL(kind=dp)                            :: diis_error

   DO ispin=1,SIZE(matrix_ks)
     CALL copy_sm_to_fm(matrix_ks(ispin)%matrix,scf_env%scf_work1(ispin)%matrix)
   ENDDO
   IF (scf_env%iter_count > 1) THEN
      CALL qs_diis_b_step(scf_env%scf_diis_buffer,mos,scf_env%scf_work1,&
           scf_env%scf_work2,scf_env%iter_delta,diis_error,diis_step,&
           scf_control%eps_diis,scf_section=scf_section,error=error)
   END IF
   do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
      ((scf_control%density_guess == core_guess).OR.(scf_env%iter_count > 1)))
   IF (diis_step) THEN
      scf_env%p_mix = 1.0_dp
      scf_env%iter_param = diis_error
      scf_env%iter_method = "DIIS/Diag."
   ELSE
      scf_env%iter_method = "Mixing/Diag."
   END IF
   scf_env%iter_delta=0.0_dp
   DO ispin=1,SIZE(matrix_ks)
      CALL simple_eigensolver(scf_env%scf_work1(ispin)%matrix,&
           mos(ispin)%mo_set,scf_env%scf_work2,&
           do_level_shift,scf_control%level_shift,&
           work_syevx=scf_control%work_syevx,&
           smear=scf_control%smear,error=error)
      CALL calculate_density_matrix(mos(ispin)%mo_set,&
           scf_env%p_mix_new(ispin)%matrix,error=error)
   ENDDO

END SUBROUTINE qs_scf_loop_do_special_diag

!!****f* qs_scf/qs_scf_loop_do_ot *
!!
!!   NAME
!!      qs_scf_loop_do_ot
!!
!!   FUNCTION
!!      the inner loop of scf, specific to using to the orbital transformation method
!!      basically, in goes the ks matrix out goes a new p matrix
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE qs_scf_loop_do_ot(scf_env,mos,rho,mo_derivs,total_energy,matrix_s,energy_only,has_unit_metric,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs
    REAL(KIND=dp), INTENT(IN)                :: total_energy
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    LOGICAL, INTENT(OUT)                     :: energy_only
    LOGICAL, INTENT(IN)                      :: has_unit_metric
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ispin
    TYPE(real_matrix_type), POINTER          :: orthogonality_metric

! in case of LSD the first spin qs_ot_env will drive the minimization
! in the case of a restricted calculation, it will make sure the spin orbitals are equal

   IF (has_unit_metric) THEN
      NULLIFY(orthogonality_metric)
   ELSE
      orthogonality_metric=>matrix_s(1)%matrix
   ENDIF

   CALL ot_scf_mini(mos,mo_derivs,orthogonality_metric, &
                    total_energy, energy_only,scf_env%iter_delta, &
                    scf_env%qs_ot_env)

   DO ispin=1,SIZE(mos)
      CALL calculate_density_matrix(mos(ispin)%mo_set,&
            rho%rho_ao(ispin)%matrix,error=error)
   ENDDO

   scf_env%iter_method=scf_env%qs_ot_env(1)%OT_METHOD_FULL
   scf_env%iter_param=scf_env%qs_ot_env(1)%ds_min

END SUBROUTINE qs_scf_loop_do_ot

!!****f* qs_scf/qs_scf_print_summary *
!!
!!   NAME
!!     qs_scf_print_summary
!!   FUNCTION
!!     writes rather detailed summary of densities and energies 
!!     after the SCF
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE qs_scf_print_summary(output_unit,rho,qs_charges,energy,scf_env,dft_control,qmmm,gapw,gapw_xc)

    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(dft_control_type), POINTER          :: dft_control
    LOGICAL, INTENT(IN)                      :: qmmm, gapw, gapw_xc

    INTEGER                                  :: ispin
    REAL(kind=dp)                            :: tot1_h, tot1_s

    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
            "Total electronic density (r-space): ",&
            SUM(rho%tot_rho_r),&
            SUM(rho%tot_rho_r)+ REAL(scf_env%nelectron,dp),&
            "Total core charge density (r-space):",&
            qs_charges%total_rho_core_rspace,&
            qs_charges%total_rho_core_rspace - REAL(scf_env%nelectron+dft_control%charge,dp)
       IF(gapw) THEN
          tot1_h =  qs_charges%total_rho1_hard(1)
          tot1_s =  qs_charges%total_rho1_soft(1)
          DO ispin=2,dft_control%nspins
            tot1_h = tot1_h + qs_charges%total_rho1_hard(ispin)
            tot1_s = tot1_s + qs_charges%total_rho1_soft(ispin)
          END DO
          WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
                "Hard and soft densities (Lebedev):",&
                tot1_h, tot1_s
          WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
               "Total Rho_soft + Rho1_hard - Rho1_soft (r-space): ",&
               SUM(rho%tot_rho_r)+ tot1_h - tot1_s ,&
               "Total charge density (r-space):      ",&
               SUM(rho%tot_rho_r)+ tot1_h - tot1_s &
               + qs_charges%total_rho_core_rspace,&
               "Total Rho_soft + Rho0_soft (g-space):",&
               qs_charges%total_rho_gspace
       ELSE
         WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Total charge density (r-space):     ",&
            SUM(rho%tot_rho_r)+&
            qs_charges%total_rho_core_rspace,&
            "Total charge density (g-space):     ",qs_charges%total_rho_gspace
       END IF
       IF (dft_control%qs_control%semi_empirical) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Core-core repulsion energy [eV]:               ",energy%core_overlap*evolt,&
            "Core Hamiltonian energy [eV]:                  ",energy%core*evolt,&
            "Two-electron integral energy [eV]:             ",energy%hartree*evolt,&
            "Electronic energy [eV]:                        ",&
                                              (energy%core+0.5_dp*energy%hartree)*evolt
       ELSE
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Overlap energy of the core charge distribution:",energy%core_overlap,&
            "Self energy of the core charge distribution:   ",energy%core_self,&
            "Core Hamiltonian energy:                       ",energy%core,&
            "Hartree energy:                                ",energy%hartree,&
            "Exchange-correlation energy:                   ",energy%exc
         IF(gapw) THEN
           WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
              "GAPW| Exc from hard and soft atomic rho1:      ",energy%exc1,&
              "GAPW| local Eh = 1 center integrals:           ",energy%hartree_1c,&
              "GAPW| local Eh = 2 center integrals:           ",energy%hartree_2c,&
              "GAPW| local Eh = 3 centers integrals:          ",energy%hartree_3c
         END IF
         IF(gapw_xc) THEN
           WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
              "GAPW_XC| Exc from hard and soft atomic rho1:      ",energy%exc1
         END IF
       END IF
       IF  (qmmm) THEN
          WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
               "QM/MM Electrostatic energy:                    ",energy%qmmm_el
       END IF
       IF (dft_control%qs_control%mulliken_restraint) THEN
          WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                "Mulliken restraint energy: ",energy%mulliken
       END IF
       IF (dft_control%qs_control%semi_empirical) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Total energy [eV]:                             ",energy%total*evolt
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Atomic reference energy [eV]:                  ",energy%core_self*evolt,&
            "Heat of formation [kcal/mol]:                  ",&
            (energy%total+energy%core_self)*kcalmol
       ELSE
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Total energy:                                  ",energy%total
       END IF
       CALL m_flush(output_unit)
    END IF

END SUBROUTINE qs_scf_print_summary

!!****f* qs_scf/qs_scf_loop_print *
!!
!!   NAME
!!      qs_scf_loop_print
!!
!!   FUNCTION
!!      collects the 'heavy duty' printing tasks out of the SCF loop
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE qs_scf_loop_print(qs_env,scf_env,globenv,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_loop_print', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, iw
    REAL(KIND=dp)                            :: mo_mag_max, mo_mag_min, &
                                                orthonormality
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_p, matrix_s
    TYPE(section_vals_type), POINTER         :: dft_section, input, &
                                                scf_section

    logger => cp_error_get_logger(error)
    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env,input=input, dft_control=dft_control, &
                    error=error)

    dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
    scf_section => section_vals_get_subs_vals(dft_section,"SCF",error=error)


    CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, &
                    matrix_ks=matrix_ks,rho=rho,mos=mos, &
                    error=error)

    matrix_p => rho%rho_ao

    DO ispin=1,dft_control%nspins

      IF (BTEST(cp_print_key_should_output(logger%iter_info,&
           dft_section,"PRINT%AO_MATRICES/DENSITY",error=error),cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%AO_MATRICES/DENSITY",&
              extension=".Log",error=error)
         CALL write_sparse_matrix(matrix_p(ispin)%matrix,4,6,qs_env,globenv,output_unit=iw)
         CALL cp_print_key_finished_output(iw,logger,dft_section,&
              "PRINT%AO_MATRICES/DENSITY", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info,&
           dft_section,"PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",error=error),cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",&
              extension=".Log",error=error)
         IF (dft_control%qs_control%semi_empirical) THEN
            CALL write_sparse_matrix(matrix_ks(ispin)%matrix,4,6,qs_env,globenv,&
                 scale=evolt,output_unit=iw)
         ELSE
            CALL write_sparse_matrix(matrix_ks(ispin)%matrix,4,6,qs_env,globenv,output_unit=iw)
         END IF
         CALL cp_print_key_finished_output(iw,logger,dft_section,&
              "PRINT%AO_MATRICES/KOHN_SHAM_MATRIX", error=error)
      END IF

    ENDDO

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         scf_section,"PRINT%MO_ORTHONORMALITY",error=error),cp_p_file)) THEN
       IF(scf_env%method==special_diag_method_nr) THEN
         CALL calculate_orthonormality(orthonormality,mos)
       ELSE
         CALL calculate_orthonormality(orthonormality,mos,matrix_s(1)%matrix)
       END IF
       iw=cp_print_key_unit_nr(logger,scf_section,"PRINT%MO_ORTHONORMALITY",&
            extension=".scfLog",error=error)
       IF (iw>0) THEN
          WRITE(iw,'(T8,A,T60,E20.4)')  &
                 " Maximum deviation from MO S-orthonormality",orthonormality
       ENDIF
       CALL cp_print_key_finished_output(iw,logger,scf_section,&
            "PRINT%MO_ORTHONORMALITY", error=error)
    ENDIF

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         scf_section,"PRINT%MO_MAGNITUDE",error=error),cp_p_file)) THEN
       CALL calculate_magnitude(mos,mo_mag_min,mo_mag_max)
       iw=cp_print_key_unit_nr(logger,scf_section,"PRINT%MO_MAGNITUDE",&
            extension=".scfLog",error=error)
       IF (iw>0) THEN
          WRITE(iw,'(T8,A,T40,2E20.4)')  &
                 " Minimum/Maximum MO magnitude ",mo_mag_min,mo_mag_max
       ENDIF
       CALL cp_print_key_finished_output(iw,logger,scf_section,&
            "PRINT%MO_MAGNITUDE", error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_loop_print

!!****f* qs_scf/init_scf_loop *
!!
!!   NAME
!!       init_scf_loop
!!
!!   FUNCTION
!!       inits those objects needed if you want to restart the scf with, say
!!       only a new initial guess, or different density functional or ...
!!       this will happen just before the scf loop starts
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE init_scf_loop(scf_env,qs_env,globenv,scf_section,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_scf_loop', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, &
                                                number_of_OT_envs, stat
    LOGICAL                                  :: do_rotation, failure, &
                                                has_unit_metric
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(real_matrix_type), POINTER          :: orthogonality_metric
    TYPE(scf_control_type), POINTER          :: scf_control

    CALL timeset("init_scf_loop","I"," ",handle)

    NULLIFY(scf_control,matrix_s,matrix_ks,dft_control,mos,&
         mo_coeff)

    failure=.FALSE.
    
    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,&
         scf_control=scf_control,&
         dft_control=dft_control,&
         mos=mos,matrix_ks=matrix_ks,&
         matrix_s=matrix_s, error=error)

    SELECT CASE (scf_env%method)
    CASE DEFAULT
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown scf method method:"//&
               cp_to_string(scf_env%method),error,failure)
    CASE (general_diag_method_nr,special_diag_method_nr)
       IF (.NOT.ASSOCIATED(scf_env%scf_diis_buffer)) THEN       
          CALL qs_diis_b_create(scf_env%scf_diis_buffer,&
               nbuffer=scf_control%max_diis,&
               error=error)
       END IF
       CALL qs_diis_b_clear(scf_env%scf_diis_buffer,error=error)
    CASE (ot_method_nr)
        has_unit_metric=.FALSE.
        IF (dft_control%qs_control%semi_empirical) THEN
            IF (dft_control%qs_control%se_control%orthogonal_basis) THEN
                has_unit_metric=.TRUE.
            ENDIF
        END IF
       IF (.NOT.ASSOCIATED(scf_env%qs_ot_env)) THEN

          ! restricted calculations require just on set of OT orbitals
          number_of_OT_envs=dft_control%nspins
          IF (dft_control%restricted) number_of_OT_envs=1

          ALLOCATE(scf_env%qs_ot_env(number_of_OT_envs),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          ! XXX Joost XXX should disentangle reading input from this part
          CALL ot_scf_read_input(scf_env%qs_ot_env,globenv)

          ! keep a note that we are restricted
          IF (dft_control%restricted) THEN
              scf_env%qs_ot_env(1)%restricted=.TRUE.
              ! requires rotation
              CPPrecondition(scf_env%qs_ot_env(1)%settings%do_rotation,cp_failure_level,routineP,error,failure) 
          ELSE 
              scf_env%qs_ot_env(:)%restricted=.FALSE.
          ENDIF

          ! might need the KS matrix to init properly
          CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,&
               error=error, calculate_forces=.FALSE.,just_energy=.FALSE.)

          IF (.NOT.ASSOCIATED(scf_env%ot_preconditioner)) THEN
               IF (scf_env%qs_ot_env(1)%settings%preconditioner_type == ot_precond_full_all) THEN
                  ALLOCATE(scf_env%ot_preconditioner(dft_control%nspins), stat=stat)
               ELSE
                  ALLOCATE(scf_env%ot_preconditioner(1), stat=stat)
               ENDIF
               CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
               DO i=1,SIZE(scf_env%ot_preconditioner)
                  ALLOCATE(scf_env%ot_preconditioner(i)%preconditioner)
                  CALL init_preconditioner(scf_env%ot_preconditioner(i)%preconditioner,&
                                           para_env=qs_env%para_env,&
                                           blacs_env=qs_env%blacs_env)
               ENDDO
          END IF

          IF (scf_env%qs_ot_env(1)%settings%preconditioner_type /= ot_precond_none) THEN 
             !
             ! preconditioning still needs to be done correctly with has_unit_metric
             ! notice that a big part of the preconditioning (S^-1) is fine anyhow
             !
             CPPrecondition(.NOT. has_unit_metric,cp_failure_level,routineP,error,failure) 
             ! constuct the preconditioner based on the first spin set. 
             ! this might or might not be reasonable. Presumably it is.

             IF (scf_env%qs_ot_env(1)%settings%preconditioner_type == ot_precond_full_all) THEN
                ! this will rotate the MOs to be eigen states, which is not compatible with rotation
                do_rotation=scf_env%qs_ot_env(1)%settings%do_rotation

               DO ispin=1,dft_control%nspins
                  CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
                  CALL calculate_subspace_eigenvalues(mo_coeff,matrix_ks(ispin)%matrix,&
                            para_env=mo_coeff%matrix_struct%para_env, &
                            do_rotation = .TRUE., co_rotate=qs_env%mo_derivs(ispin)%matrix)
                  CALL make_preconditioner(scf_env%ot_preconditioner(ispin)%preconditioner, &
                         scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                         matrix_ks(1)%matrix,matrix_s(1)%matrix,&
                         qs_env%kinetic(1)%matrix, &
                         mo_coeff,scf_env%qs_ot_env(1)%settings%energy_gap) 
               ENDDO

             ELSE

               CALL get_mo_set(mo_set=mos(1)%mo_set, mo_coeff=mo_coeff)
               CALL make_preconditioner(scf_env%ot_preconditioner(1)%preconditioner, &
                       scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                       matrix_ks(1)%matrix,matrix_s(1)%matrix,&
                       qs_env%kinetic(1)%matrix, &
                       mo_coeff,scf_env%qs_ot_env(1)%settings%energy_gap) 

             ENDIF

          ENDIF

          IF (has_unit_metric) THEN
             NULLIFY(orthogonality_metric)
          ELSE
             orthogonality_metric=>matrix_s(1)%matrix
          ENDIF

          CALL ot_scf_init(mos,orthogonality_metric,scf_env%qs_ot_env)
          
          IF (scf_env%qs_ot_env(1)%settings%preconditioner_type /= ot_precond_none) THEN
             IF (scf_env%qs_ot_env(1)%settings%preconditioner_type == ot_precond_full_all) THEN
                DO ispin=1,SIZE(scf_env%qs_ot_env)
                    CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin),&
                                                  scf_env%ot_preconditioner(ispin)%preconditioner)
                ENDDO
          ELSE
                DO ispin=1,SIZE(scf_env%qs_ot_env)
                    CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin),&
                                                 scf_env%ot_preconditioner(1)%preconditioner)
                ENDDO
             ENDIF
          ENDIF
       ENDIF
    END SELECT

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE init_scf_loop

!!****f* qs_scf/init_scf_run *
!!
!!   NAME
!!      init_scf_run
!!
!!   FUNCTION
!!      performs those initialisations that need to be done only once
!!      (e.g. that only depend on the atomic positions)
!!      this will be called in scf
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE init_scf_run(scf_env, qs_env, globenv, scf_section,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_scf_run', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homo, i, ispin, iw, &
                                                nao, ndep, nelectron_spin, &
                                                output_unit
    LOGICAL                                  :: failure, &
                                                gth_potential_present, &
                                                id_equal, &
                                                my_transition_potential
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control

    CALL timeset("init_scf_run","I"," ",handle)

    NULLIFY(scf_control, atomic_kind_set, matrix_h, matrix_s, matrix_ks, &
         dft_control, mos, rho)
    failure=.FALSE.
    my_transition_potential = .FALSE.
    logger => cp_error_get_logger(error)
    
    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL scf_env_check_i_alloc(scf_env=scf_env, qs_env=qs_env,&
         error=error)

    CALL get_qs_env(qs_env=qs_env,&
         scf_control=scf_control,&
         dft_control=dft_control,&
         atomic_kind_set=atomic_kind_set,&
         mos=mos,matrix_ks=matrix_ks, rho=rho, &
         matrix_h=matrix_h,matrix_s=matrix_s, error=error)

    ! some special checks to eliminate the first problems with restricted
    ! should move earlier I think
    IF (dft_control%restricted) THEN
       IF (scf_env%method .NE. ot_method_nr) THEN
          IF (output_unit>0) WRITE(output_unit,*) "OT only for restricted"
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
       CPPrecondition(dft_control%nspins.EQ.2,cp_failure_level,routineP,error,failure)
    ENDIF

    ! do some assertions here on these matrices having the same structure
    id_equal=ALL((/((matrix_s(i)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id),i=1,SIZE(matrix_s))/))
    CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    id_equal=ALL((/((rho%rho_ao(i)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id),i=1,SIZE(rho%rho_ao))/))
    CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    id_equal=ALL((/((matrix_h(i)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id),i=1,SIZE(matrix_h))/))
    CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    id_equal=ALL((/((matrix_ks(i)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id),i=1,SIZE(matrix_ks))/))
    CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         gth_potential_present=gth_potential_present)

! updates the total number of electrons
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         nelectron=scf_env%nelectron)
    scf_env%nelectron = scf_env%nelectron - dft_control%charge
! print occupation numbers
    output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".scfLog",error=error)
    IF (output_unit>0) THEN
       DO ispin=1,dft_control%nspins
          CALL get_mo_set(mo_set=mos(ispin)%mo_set,&
               homo=homo, &
               nelectron=nelectron_spin, &
               nao=nao)
          IF (dft_control%nspins > 1) THEN
            WRITE (UNIT=output_unit,FMT="(/,T3,A,I2)") "Spin",ispin
          END IF
          WRITE (UNIT=output_unit,FMT="(/,(T3,A,T31,I10))")&
            "Number of electrons:",nelectron_spin,&
            "Number of occupied orbitals:",homo
       ENDDO
       WRITE (UNIT=output_unit,FMT="((T3,A,T31,I10))")&
         "Number of orbital functions:",nao
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
         "PRINT%PROGRAM_RUN_INFO", error=error)

! calc ortho matrix
    ndep = 0
    IF (scf_env%method /= ot_method_nr .AND. &
        scf_env%method /= special_diag_method_nr ) THEN
       CALL copy_sm_to_fm(matrix_s(1)%matrix,scf_env%ortho)
       IF (scf_control%use_cholesky ) THEN
          CALL cp_fm_cholesky_decompose(scf_env%ortho)
       ELSE 
          CALL cp_fm_power(scf_env%ortho,scf_env%scf_work2,-0.5_dp,scf_control%eps_eigval,ndep)
       ENDIF
       IF (BTEST(cp_print_key_should_output(logger%iter_info,&
            qs_env%input,"DFT%PRINT%AO_MATRICES/ORTHO",error=error),cp_p_file)) THEN
          iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/ORTHO",&
               extension=".Log",error=error)
          CALL write_fm_with_basis_info(scf_env%ortho,4,6,qs_env,globenv,output_unit=iw)
          CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
               "DFT%PRINT%AO_MATRICES/ORTHO", error=error)
       END IF
    ENDIF
    output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
               extension=".scfLog",error=error)
    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(T3,A,I10)")&
                  "Number of independent orbital functions:",nao - ndep
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
          "PRINT%PROGRAM_RUN_INFO", error=error)

    my_transition_potential = (dft_control%do_xas_calculation .AND. dft_control%xas_estate>0)
    ! initializes rho and the mos
    IF( my_transition_potential) THEN
      ! initialize the density with the localized mos
        CALL loc_initialize_rho(qs_env,error=error)  
    ELSE
      CALL scf_env_initial_rho_setup(scf_env,qs_env=qs_env,&
           globenv=globenv, scf_section=scf_section, error=error)
    END IF

    ! Frozen density approximation
    IF (ASSOCIATED(qs_env%wf_history)) THEN
      IF (qs_env%wf_history%interpolation_method_nr==wfi_frozen_method_nr) THEN
        IF (.NOT. ASSOCIATED(qs_env%wf_history%past_states(1)%snapshot)) THEN
          CALL wfi_update(qs_env%wf_history, qs_env=qs_env, dt=1.0_dp, &
                          error=error)
          CALL duplicate_rho_type(rho_input=rho, &
                                  rho_output=qs_env%wf_history%past_states(1)%snapshot%rho_frozen, &
                                  qs_env=qs_env, error=error)
        END IF
      END IF
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE init_scf_run

!!****f* qs_scf/scf_env_check_i_alloc [1.0] *
!!
!!   NAME
!!     scf_env_check_i_alloc
!!
!!   FUNCTION
!!     checks the allocation status of the needed matrixes, and if necessary
!!     allocate them
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - scf_env: the scf_env to be checked
!!     - qs_env: the qs_env, the scf_env lives in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf_env_check_i_alloc(scf_env,qs_env,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scf_env_check_i_alloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homo, ispin, nao, &
                                                nrow_block, stat
    LOGICAL                                  :: failure, &
                                                gth_potential_present, &
                                                my_transition_potential, &
                                                uniform_occupation
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_fm_struct_type), POINTER         :: ao_ao_fmstruct, ao_mo_fmstruct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: S_sm_pools
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I"," ",handle)

    NULLIFY(matrix_ks, ao_mo_fm_pools, matrix_s, ao_mo_fmstruct, ao_ao_fmstruct,&
            dft_control, mos, occupation_numbers, S_sm_pools, ks_env)
    NULLIFY(atomic_kind_set, mo_derivs, mo_coeff)

    failure=.FALSE.

    my_transition_potential = .FALSE.
    uniform_occupation = .TRUE.

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         mos=mos,&
         matrix_ks=matrix_ks,&
         ks_env=ks_env,&
         atomic_kind_set=atomic_kind_set,&
         matrix_s=matrix_s)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         gth_potential_present=gth_potential_present)
    CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools,&
         S_sm_pools=S_sm_pools,&
         error=error)

    ! very first tests for xas

    my_transition_potential = (dft_control%do_xas_calculation .AND. dft_control%xas_estate>0)
!   *** finish initialization of the MOs ***
    CPPrecondition(ASSOCIATED(mos),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       DO ispin=1,SIZE(mos)
          CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,homo=homo)
          IF (.NOT.ASSOCIATED(mo_coeff)) THEN
             CALL init_mo_set(mos(ispin)%mo_set,&
                  ao_mo_fm_pools(ispin)%pool,&
                  name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                  "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                  error=error)
          END IF
          IF(my_transition_potential .AND. ispin==1) THEN
            CALL get_mo_set(mos(ispin)%mo_set,&
                 occupation_numbers=occupation_numbers)
            occupation_numbers(dft_control%xas_estate) = &
                 dft_control%xas_control%occ_estate
            occupation_numbers(homo) = dft_control%xas_control%occ_homo
            uniform_occupation = .FALSE.

            CALL set_mo_set(mos(ispin)%mo_set,&
                 uniform_occupation=uniform_occupation)
          END IF
       END DO
    END IF

!   *** get the mo_derivs OK if needed ***
    IF (qs_env%requires_mo_derivs) THEN
     CALL get_qs_env(qs_env,mo_derivs=mo_derivs)
     IF (.NOT.ASSOCIATED(mo_derivs)) THEN
       IF (dft_control%restricted) THEN ! right now, there might be more mos than needed derivs
          ALLOCATE(mo_derivs(1))
          CALL get_mo_set(mos(1)%mo_set,mo_coeff=mo_coeff)
          CALL cp_fm_create(mo_derivs(1)%matrix,mo_coeff%matrix_struct)
       ELSE
          ALLOCATE(mo_derivs(dft_control%nspins))
          DO ispin=1,dft_control%nspins
             CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
             CALL cp_fm_create(mo_derivs(ispin)%matrix,mo_coeff%matrix_struct)
          ENDDO
       ENDIF
       CALL set_qs_env(qs_env,mo_derivs=mo_derivs)
     ENDIF
    ELSE
       ! nothing should be done
    ENDIF


!   *** Allocate the distributed SCF matrices ***

    IF ((.NOT.ASSOCIATED(scf_env%scf_work1)).OR.&
        (.NOT.ASSOCIATED(scf_env%scf_work2)).OR.&
        (.NOT.ASSOCIATED(scf_env%ortho))) THEN

       ao_mo_fmstruct => fm_pool_get_el_struct(ao_mo_fm_pools(1)%pool,&
            error=error)
       CALL cp_fm_struct_get(ao_mo_fmstruct, nrow_block=nrow_block,&
            error=error)
       CALL get_mo_set(mos(1)%mo_set,nao=nao)
       CALL cp_fm_struct_create(fmstruct=ao_ao_fmstruct,&
            nrow_block=nrow_block,ncol_block=nrow_block,&
            nrow_global=nao, ncol_global=nao,&
            template_fmstruct=ao_mo_fmstruct, error=error)


       IF (.NOT.ASSOCIATED(scf_env%scf_work1)) THEN
          IF (scf_env%method/=ot_method_nr) THEN
             ALLOCATE(scf_env%scf_work1(dft_control%nspins), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(scf_env%scf_work1)) THEN
             DO ispin=1,SIZE(scf_env%scf_work1)
                NULLIFY(scf_env%scf_work1(ispin)%matrix)
                CALL cp_fm_create(scf_env%scf_work1(ispin)%matrix,&
                     matrix_struct=ao_ao_fmstruct,&
                     name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
                     "WORK_MATRIX-1-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                     error=error)
             ENDDO
          ENDIF
       END IF
       IF ((.NOT.ASSOCIATED(scf_env%ortho)).AND.&
           ((scf_env%method /= ot_method_nr).AND.&
            (scf_env%method /= special_diag_method_nr))) THEN
          CALL cp_fm_create(scf_env%ortho,&
               matrix_struct=ao_ao_fmstruct,&
               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
               "ORTHO_MATRIX",&
               error=error)
       END IF
       IF ((.NOT.ASSOCIATED(scf_env%scf_work2)).AND.&
           (scf_env%method /= ot_method_nr)) THEN
          CALL cp_fm_create(scf_env%scf_work2,&
               matrix_struct=ao_ao_fmstruct,&
               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
               "WORK_MATRIX-2",&
               error=error)
       END IF
       CALL cp_fm_struct_release(ao_ao_fmstruct,error=error)

    END IF

!   *** Allocate matrix_ks and put it in the QS environment ***

    IF (.not.ASSOCIATED(matrix_ks)) THEN
       CALL sm_pools_create_matrix_vect(S_sm_pools,matrix_ks,&
            name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
            "KOHN-SHAM_MATRIX",&
            error=error)
       CALL set_qs_env(qs_env=qs_env,&
            matrix_ks=matrix_ks)
    END IF

!   *** allocate p_mix_new ***
    IF (scf_env%mixing_method.EQ.ao_mixing_nr) THEN
       IF (.NOT.ASSOCIATED(scf_env%p_mix_new)) THEN
          CALL sm_pools_create_matrix_vect(S_sm_pools,scf_env%p_mix_new,&
               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
               "DENSITY",&
               error=error)
       END IF
    END IF

!   *** allocate the ks env **
    IF (.NOT.ASSOCIATED(ks_env)) THEN
       CALL qs_ks_create(ks_env,qs_env=qs_env)
       CALL set_qs_env(qs_env, ks_env=ks_env)
       CALL qs_ks_release(ks_env)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE scf_env_check_i_alloc

!!****f* qs_scf/scf_env_initial_rho_setup [1.0] *
!!
!!   NAME
!!     scf_env_initial_rho_setup
!!
!!   FUNCTION
!!     Initializes rho and the mos, so that an scf cycle can start
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - scf_env: the scf env in which to do the scf
!!     - qs_env: the qs env the scf_env lives in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf_env_initial_rho_setup(scf_env, qs_env, globenv, scf_section, error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_initial_rho_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: extrapolation_method_nr, &
                                                handle, ispin, nmo, &
                                                output_unit
    LOGICAL                                  :: failure, gapw, gapw_xc, &
                                                orthogonal_wf
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,"I"," ",handle)
    failure=.FALSE.
    NULLIFY(mo_coeff)
    gapw = qs_env%dft_control%qs_control%gapw
    gapw_xc = qs_env%dft_control%qs_control%gapw_xc
    logger => cp_error_get_logger(error)
    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
       extrapolation_method_nr=wfi_use_guess_method_nr
       IF (ASSOCIATED(qs_env%wf_history)) THEN
          CALL wfi_extrapolate(globenv,qs_env%wf_history, &
               qs_env=qs_env, dt=1.0_dp, &
               extrapolation_method_nr=extrapolation_method_nr,&
               orthogonal_wf=orthogonal_wf, gapw=gapw, gapw_xc=gapw_xc, error=error)
          ! wfi_use_guess_method_nr the wavefunctions are not yet initialized
          IF ((.NOT.orthogonal_wf).AND.&
              (scf_env%method == ot_method_nr).AND.&
              (.NOT.(extrapolation_method_nr == wfi_use_guess_method_nr))) THEN
             DO ispin=1,SIZE(qs_env%mos)
                CALL get_mo_set(qs_env%mos(ispin)%mo_set, &
                     mo_coeff=mo_coeff, nmo=nmo)
                CALL qs_env_reorthogonalize_vectors(qs_env, &
                     v_matrix=mo_coeff, ispin=ispin, n_col=nmo,&
                     error=error)
             END DO
          END IF
       END IF
       output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
            extension=".scfLog",error=error)
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,T3,A)")&
               "Extrapolation method: "//&
               TRIM(wfi_get_method_label(extrapolation_method_nr,error=error))
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
            "PRINT%PROGRAM_RUN_INFO", error=error)
       IF (extrapolation_method_nr==wfi_use_guess_method_nr) THEN
          CALL calculate_first_density_matrix(scf_env=scf_env,qs_env=qs_env,&
               globenv=globenv, error=error)
          IF (.NOT.(qs_env%scf_control%density_guess==densities_guess)) THEN
            CALL qs_rho_update_rho(qs_env%rho,qs_env=qs_env, gapw=gapw, error=error)
            IF(gapw) THEN
               CALL calculate_rho_atom_coeff(qs_env,qs_env%rho%rho_ao)
            ENDIF
            IF(gapw_xc) THEN
               CALL qs_rho_update_rho(qs_env%rho_xc,qs_env=qs_env, gapw=gapw_xc, error=error)
               CALL calculate_rho_atom_coeff(qs_env,qs_env%rho%rho_ao)
            END IF
            CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.,&
                 error=error)
          END IF
       END IF
    END IF

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE scf_env_initial_rho_setup
!***************************************************************************

!!****f* qs_scf/scf_env_cleanup [1.0] *
!!
!!   NAME
!!     scf_env_cleanup
!!
!!   FUNCTION
!!     perform cleanup operations (like releasing temporary storage)
!!     at the end of the scf
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf_env_cleanup(scf_env,qs_env,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, stat
    LOGICAL                                  :: failure
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: S_sm_pools

    CALL timeset(routineN,"I"," ",handle)

    failure=.FALSE.
    NULLIFY(S_sm_pools)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL mpools_get(qs_env%mpools,S_sm_pools=S_sm_pools,error=error)

!   *** Release SCF work storage ***

       IF (ASSOCIATED(scf_env%scf_work1)) THEN
          DO ispin=1,SIZE(scf_env%scf_work1)
             CALL cp_fm_release(scf_env%scf_work1(ispin)%matrix)
          ENDDO
          DEALLOCATE(scf_env%scf_work1)
       ENDIF
       IF (ASSOCIATED(scf_env%scf_work2)) CALL cp_fm_release(scf_env%scf_work2)
       IF (ASSOCIATED(scf_env%ortho)) CALL cp_fm_release(scf_env%ortho)

       IF (ASSOCIATED(scf_env%p_mix_new)) THEN
          CALL sm_pools_give_back_matrix_vect(S_sm_pools,scf_env%p_mix_new,&
               error=error)
          CALL sm_pools_flush_cache(S_sm_pools,error=error)
       END IF

! *** method dependent cleanup
       SELECT CASE(scf_env%method)
       CASE(ot_method_nr)
          DO ispin=1,SIZE(scf_env%ot_preconditioner)
             CALL destroy_preconditioner(scf_env%ot_preconditioner(ispin)%preconditioner) 
             DEALLOCATE(scf_env%ot_preconditioner(ispin)%preconditioner)
          ENDDO
          DEALLOCATE(scf_env%ot_preconditioner,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DO ispin=1,SIZE(scf_env%qs_ot_env)
             CALL ot_scf_destroy(scf_env%qs_ot_env(ispin))
          ENDDO
          DEALLOCATE(scf_env%qs_ot_env,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       CASE(general_diag_method_nr)
          !
       CASE(special_diag_method_nr)
          !
       CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown scf method method:"//&
               cp_to_string(scf_env%method),error,failure)
       END SELECT

    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE scf_env_cleanup
!***************************************************************************

END MODULE qs_scf
