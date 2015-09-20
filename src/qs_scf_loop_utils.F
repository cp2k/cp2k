!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief Utility routines for qs_scf
! *****************************************************************************
MODULE qs_scf_loop_utils
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_name,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_external_control,             ONLY: external_control
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE kpoint_types,                    ONLY: kpoint_type
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
  USE qs_fb_env_methods,               ONLY: fb_env_do_diag
  USE qs_gspace_mixing,                ONLY: gspace_mixing,&
                                             self_consistency_check
  USE qs_ks_types,                     ONLY: qs_ks_did_change,&
                                             qs_ks_env_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_occupation,                ONLY: set_mo_occupation
  USE qs_mo_types,                     ONLY: mo_set_p_type
  USE qs_ot_scf,                       ONLY: ot_scf_destroy,&
                                             ot_scf_mini
  USE qs_outer_scf,                    ONLY: outer_loop_gradient
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE qs_scf_diagonalization,          ONLY: do_block_davidson_diag,&
                                             do_block_krylov_diag,&
                                             do_general_diag,&
                                             do_general_diag_kp,&
                                             do_ot_diag,&
                                             do_roks_diag,&
                                             do_scf_diag_subspace,&
                                             do_special_diag
  USE qs_scf_methods,                  ONLY: scf_env_density_mixing
  USE qs_scf_output,                   ONLY: qs_scf_print_summary
  USE qs_scf_types,                    ONLY: block_davidson_diag_method_nr,&
                                             block_krylov_diag_method_nr,&
                                             filter_matrix_diag_method_nr,&
                                             general_diag_method_nr,&
                                             ot_diag_method_nr,&
                                             ot_method_nr,&
                                             qs_scf_env_type,&
                                             special_diag_method_nr
  USE scf_control_types,               ONLY: scf_control_type,&
                                             smear_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_loop_utils'

  PUBLIC :: qs_scf_set_loop_flags,&
            qs_scf_new_mos, qs_scf_new_mos_kp,&
            qs_scf_density_mixing, qs_scf_check_inner_exit, &
            qs_scf_check_outer_exit, qs_scf_inner_finalize, qs_scf_rho_update

CONTAINS


! *****************************************************************************
!> \brief computes properties for a given hamiltonian using the current wfn
!> \param scf_env ...
!> \param diis_step ...
!> \param energy_only ...
!> \param just_energy ...
!> \param exit_inner_loop ...
! *****************************************************************************
  SUBROUTINE qs_scf_set_loop_flags(scf_env,diis_step,&
                      energy_only,just_energy,exit_inner_loop)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    LOGICAL                                  :: diis_step, energy_only, &
                                                just_energy, exit_inner_loop

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_set_loop_flags', &
      routineP = moduleN//':'//routineN

! Some flags needed to be set at the beginning of the loop

    diis_step = .FALSE.
    energy_only = .FALSE.
    just_energy = .FALSE.

    ! SCF loop, optimisation of the wfn coefficients
    ! qs_env%rho%rho_r and qs_env%rho%rho_g should be up to date here

    scf_env%iter_count = 0
    exit_inner_loop    = .FALSE.

   END SUBROUTINE qs_scf_set_loop_flags


! *****************************************************************************
!> \brief takes known energy and derivatives and produces new wfns
!>        and or density matrix
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param scf_section ...
!> \param diis_step ...
!> \param energy_only ...
! *****************************************************************************
  SUBROUTINE qs_scf_new_mos(qs_env,scf_env,scf_control,scf_section,diis_step,&
                            energy_only)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL                                  :: diis_step, energy_only

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_new_mos', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin
    LOGICAL                                  :: failure, has_unit_metric, &
                                                skip_diag_sub
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho

    failure=.FALSE.
    NULLIFY(energy, ks_env, matrix_ks, matrix_s, rho, mos, dft_control)

    CALL get_qs_env(qs_env=qs_env,&
                    matrix_s=matrix_s,energy=energy,&
                    ks_env=ks_env, &
                    matrix_ks=matrix_ks,rho=rho,mos=mos, &
                    dft_control=dft_control, &
                    has_unit_metric=has_unit_metric)
    scf_env%iter_param = 0.0_dp


    SELECT CASE (scf_env%method)
    CASE DEFAULT
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
            routineP,"unknown scf method: "//&
            cp_to_string(scf_env%method),failure)

    ! *************************************************************************
    ! Filter matrix diagonalisation: ugly implementation at this point of time
    ! *************************************************************************
    CASE(filter_matrix_diag_method_nr)

       CALL fb_env_do_diag(scf_env%filter_matrix_env, qs_env, &
                            matrix_ks, matrix_s, scf_section, diis_step)

    ! Diagonlization in non orthonormal case
    CASE(general_diag_method_nr)
       IF (dft_control%roks) THEN
          CALL do_roks_diag(scf_env,mos,matrix_ks,matrix_s,&
                                        scf_control,scf_section,diis_step,&
                                        has_unit_metric)
       ELSE
          CALL do_general_diag(scf_env,mos,matrix_ks,&
                               matrix_s,scf_control,scf_section, &
                               diis_step)
          IF(scf_control%do_diag_sub) THEN
            skip_diag_sub = (scf_env%subspace_env%eps_diag_sub > 0.0_dp) .AND. &
                   (scf_env%iter_count==1  .OR. scf_env%iter_delta>scf_env%subspace_env%eps_diag_sub )
            IF( .NOT. skip_diag_sub) THEN
               CALL do_scf_diag_subspace(qs_env,scf_env,scf_env%subspace_env,mos,rho,&
                    ks_env,scf_section,scf_control)
            END IF
          END IF
       END IF
    ! Diagonlization in orthonormal case
    CASE(special_diag_method_nr)
       IF (dft_control%roks) THEN
          CALL do_roks_diag(scf_env,mos,matrix_ks,matrix_s,&
                                        scf_control,scf_section,diis_step,&
                                        has_unit_metric)
       ELSE
          CALL do_special_diag(scf_env,mos,matrix_ks,&
                                           scf_control,scf_section, &
                                           diis_step)
       END IF
    ! OT diagonlization
    CASE(ot_diag_method_nr)
       CALL do_ot_diag(scf_env,mos,matrix_ks,matrix_s,&
                            scf_control,scf_section,diis_step)
    ! Block Krylov diagonlization
    CASE(block_krylov_diag_method_nr)
       IF((scf_env%krylov_space%eps_std_diag > 0.0_dp) .AND.&
          (scf_env%iter_count==1 .OR. scf_env%iter_delta>scf_env%krylov_space%eps_std_diag)) THEN
          IF(scf_env%krylov_space%always_check_conv) THEN
             CALL do_block_krylov_diag(scf_env,mos,matrix_ks,&
                  scf_control, scf_section, check_moconv_only=.TRUE.)
          END IF
          CALL do_general_diag(scf_env,mos,matrix_ks,&
               matrix_s,scf_control,scf_section, diis_step)
       ELSE
         CALL do_block_krylov_diag(scf_env,mos,matrix_ks, &
                                   scf_control, scf_section)
       END IF
       IF(scf_control%do_diag_sub) THEN
          skip_diag_sub = (scf_env%subspace_env%eps_diag_sub > 0.0_dp) .AND. &
                   (scf_env%iter_count==1  .OR. scf_env%iter_delta>scf_env%subspace_env%eps_diag_sub )
          IF( .NOT. skip_diag_sub) THEN
               CALL do_scf_diag_subspace(qs_env,scf_env,scf_env%subspace_env,mos,rho,&
                    ks_env,scf_section, scf_control)
          END IF
       END IF
    ! Block Davidson diagonlization
    CASE(block_davidson_diag_method_nr)
       CALL do_block_davidson_diag(qs_env,scf_env,mos,matrix_ks,matrix_s,scf_control,&
           scf_section,.FALSE.)
    ! OT without diagonlization. Needs special treatment for SCP runs
    CASE(ot_method_nr)
       CALL qs_scf_loop_do_ot(qs_env,scf_env,scf_control%smear,mos,rho,&
            qs_env%mo_derivs,energy%total, &
            matrix_s, energy_only=energy_only,has_unit_metric=has_unit_metric)
    END SELECT

    energy%kTS = 0.0_dp
    energy%efermi = 0.0_dp
    CALL get_qs_env(qs_env,mos=mos)
    DO ispin=1,SIZE(mos)
       energy%kTS = energy%kTS + mos(ispin)%mo_set%kTS
       energy%efermi = energy%efermi + mos(ispin)%mo_set%mu
    ENDDO
    energy%efermi = energy%efermi /REAL(SIZE(mos),KIND=dp)


  END SUBROUTINE qs_scf_new_mos


! *****************************************************************************
!> \brief Updates MOs and density matrix using diagonalization
!>        Kpoint code
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param diis_step ...
! *****************************************************************************
  SUBROUTINE qs_scf_new_mos_kp(qs_env,scf_env,scf_control,diis_step)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control
    LOGICAL                                  :: diis_step

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_new_mos_kp', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, has_unit_metric
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(kpoint_type), POINTER               :: kpoints

    failure=.FALSE.

    NULLIFY(dft_control, kpoints, matrix_ks, matrix_s)

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    kpoints=kpoints)
    scf_env%iter_param = 0.0_dp

    CALL cp_assert(.NOT.dft_control%roks,cp_failure_level,cp_assertion_failed,&
         routineP,"KP code: ROKS method not available: ",failure)

    SELECT CASE (scf_env%method)
       CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"KP code: Unknown scf method: "//&
               cp_to_string(scf_env%method),failure)
       CASE(general_diag_method_nr)
          ! Diagonlization in non orthonormal case
          CALL get_qs_env(qs_env,matrix_ks_kp=matrix_ks,matrix_s_kp=matrix_s)
          CALL do_general_diag_kp(scf_env,matrix_ks,matrix_s,kpoints,&
               scf_control,diis_step)
       CASE(special_diag_method_nr)
          CALL get_qs_env(qs_env=qs_env,has_unit_metric=has_unit_metric)
          CPASSERT(has_unit_metric)
          ! Diagonlization in orthonormal case
!         CALL do_special_diag(scf_env,mos,matrix_ks,&
!                              scf_control,scf_section,diis_step)
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"KP code: Scf method not available: "//&
               cp_to_string(scf_env%method),failure)
       CASE(ot_diag_method_nr,&
            block_krylov_diag_method_nr,&
            block_davidson_diag_method_nr,&
            ot_method_nr)
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"KP code: Scf method not available: "//&
               cp_to_string(scf_env%method),failure)
    END SELECT

  END SUBROUTINE qs_scf_new_mos_kp


! *****************************************************************************
!> \brief the inner loop of scf, specific to using to the orbital transformation method
!>       basically, in goes the ks matrix out goes a new p matrix
!> \param qs_env ...
!> \param scf_env ...
!> \param smear ...
!> \param mos ...
!> \param rho ...
!> \param mo_derivs ...
!> \param total_energy ...
!> \param matrix_s ...
!> \param energy_only ...
!> \param has_unit_metric ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
!>      2013    moved from qs_scf [Florian Schiffmann]
! *****************************************************************************
  SUBROUTINE qs_scf_loop_do_ot(qs_env,scf_env,smear,mos,rho,mo_derivs,total_energy,&
                               matrix_s,energy_only,has_unit_metric)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(smear_type), POINTER                :: smear
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs
    REAL(KIND=dp), INTENT(IN)                :: total_energy
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    LOGICAL, INTENT(INOUT)                   :: energy_only
    LOGICAL, INTENT(IN)                      :: has_unit_metric

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_loop_do_ot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao
    TYPE(cp_dbcsr_type), POINTER             :: orthogonality_metric

    CALL timeset(routineN,handle)
    NULLIFY(rho_ao)
    failure = .FALSE.

    CALL qs_rho_get(rho, rho_ao=rho_ao)

    IF (has_unit_metric) THEN
      NULLIFY(orthogonality_metric)
    ELSE
       orthogonality_metric=>matrix_s(1)%matrix
    END IF

    ! in case of LSD the first spin qs_ot_env will drive the minimization
    ! in the case of a restricted calculation, it will make sure the spin orbitals are equal

    CALL ot_scf_mini(mos,mo_derivs,smear,orthogonality_metric, &
                    total_energy,energy_only,scf_env%iter_delta, &
                    scf_env%qs_ot_env,qs_env%input)


    DO ispin=1,SIZE(mos)
       CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,&
                              smear=smear)
    ENDDO

    DO ispin=1,SIZE(mos)
      CALL calculate_density_matrix(mos(ispin)%mo_set,&
                                    rho_ao(ispin)%matrix,&
                                    use_dbcsr=.TRUE.)
    END DO

    scf_env%iter_method=scf_env%qs_ot_env(1)%OT_METHOD_FULL
    scf_env%iter_param=scf_env%qs_ot_env(1)%ds_min
    qs_env%broyden_adaptive_sigma=scf_env%qs_ot_env(1)%broyden_adaptive_sigma

    CALL timestop(handle)

  END SUBROUTINE qs_scf_loop_do_ot


! *****************************************************************************
!> \brief Performs the requested density mixing if any needed
!> \param scf_env   Holds SCF environment information
!> \param rho       All data for the electron density
!> \param para_env  Parallel environment
!> \param diis_step Did we do a DIIS step?
!> \param do_kpoints Is this a kpoint run?
! *****************************************************************************
 SUBROUTINE qs_scf_density_mixing(scf_env,rho,para_env,diis_step,do_kpoints)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL                                  :: diis_step
    LOGICAL, INTENT(IN)                      :: do_kpoints

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_density_mixing', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: rho_ao_kp

    failure=.FALSE.
    NULLIFY(rho_ao_kp)

    CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)

    SELECT CASE(scf_env%mixing_method)
    CASE(direct_mixing_nr)
       CALL scf_env_density_mixing(scf_env%p_mix_new,&
            scf_env%mixing_store, rho_ao_kp, para_env, scf_env%iter_delta, scf_env%iter_count, &
            diis=diis_step)
    CASE(gspace_mixing_nr,pulay_mixing_nr,broyden_mixing_nr,&
         broyden_mixing_new_nr,multisecant_mixing_nr)
       ! Compute the difference p_out-p_in
       CALL self_consistency_check(rho_ao_kp,scf_env%p_delta,para_env,scf_env%p_mix_new,&
            delta=scf_env%iter_delta)
    CASE(no_mixing_nr)
    CASE DEFAULT
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
            routineP,"unknown scf mixing method: "//&
            cp_to_string(scf_env%mixing_method),failure)
    END SELECT

  END SUBROUTINE qs_scf_density_mixing


! *****************************************************************************
!> \brief checks whether exit conditions for outer loop are satisfied
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param should_stop ...
!> \param outer_loop_converged ...
!> \param exit_outer_loop ...
! *****************************************************************************
  SUBROUTINE qs_scf_check_outer_exit(qs_env,scf_env,scf_control,should_stop,&
                                outer_loop_converged,exit_outer_loop)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control
    LOGICAL                                  :: should_stop, &
                                                outer_loop_converged, &
                                                exit_outer_loop

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_check_outer_exit', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: outer_loop_eps

    outer_loop_converged=.TRUE.
    IF (scf_control%outer_scf%have_scf) THEN
       ! We have an outer SCF loop...
       scf_env%outer_scf%iter_count=scf_env%outer_scf%iter_count+1
       outer_loop_converged=.FALSE.

       CALL outer_loop_gradient(qs_env,scf_env)
       outer_loop_eps=SQRT(SUM(scf_env%outer_scf%gradient(:,scf_env%outer_scf%iter_count)**2))/ &
            SIZE(scf_env%outer_scf%gradient,1)
       IF (outer_loop_eps<scf_control%outer_scf%eps_scf) outer_loop_converged=.TRUE.
    END IF

    exit_outer_loop=should_stop .OR. outer_loop_converged.OR.&
             scf_env%outer_scf%iter_count>scf_control%outer_scf%max_scf

  END SUBROUTINE qs_scf_check_outer_exit


! *****************************************************************************
!> \brief checks whether exit conditions for inner loop are satisfied
!> \param qs_env ...
!> \param scf_env ...
!> \param scf_control ...
!> \param should_stop ...
!> \param exit_inner_loop ...
!> \param inner_loop_converged ...
!> \param output_unit ...
! *****************************************************************************
  SUBROUTINE qs_scf_check_inner_exit(qs_env,scf_env,scf_control,should_stop,&
                               exit_inner_loop,inner_loop_converged,output_unit)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control
    LOGICAL                                  :: should_stop, exit_inner_loop, &
                                                inner_loop_converged
    INTEGER                                  :: output_unit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_check_inner_exit', &
      routineP = moduleN//':'//routineN

     inner_loop_converged=.FALSE.
     exit_inner_loop = .FALSE.

     CALL external_control(should_stop,"SCF",target_time=qs_env%target_time, &
          start_time=qs_env%start_time)
     IF (scf_env%iter_delta < scf_control%eps_scf) THEN
        IF (output_unit>0) THEN
           WRITE(UNIT=output_unit,FMT="(/,T3,A,I5,A/)")&
                "*** SCF run converged in ",scf_env%iter_count," steps ***"
        END IF
        inner_loop_converged=.TRUE.
        exit_inner_loop = .TRUE.
     ELSE IF (should_stop .OR.  scf_env%iter_count >= scf_control%max_scf) THEN
        IF (output_unit>0) THEN
           WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
                "*** SCF run NOT converged ***"
        END IF
        inner_loop_converged=.FALSE.
        exit_inner_loop = .TRUE.
     END IF

  END SUBROUTINE qs_scf_check_inner_exit


! *****************************************************************************
!> \brief undoing density mixing. Important upon convergence
!> \param scf_env ...
!> \param rho ...
!> \param dft_control ...
!> \param para_env ...
!> \param diis_step ...
! *****************************************************************************
  SUBROUTINE qs_scf_undo_mixing(scf_env,rho,dft_control,para_env,diis_step)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL                                  :: diis_step

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_undo_mixing', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ic, ispin, nc
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: rho_ao_kp

     failure = .FALSE.
     NULLIFY(rho_ao_kp)

     IF (scf_env%mixing_method>0) THEN
        CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
        nc = SIZE(scf_env%p_mix_new,2)
        SELECT CASE(scf_env%mixing_method)
        CASE(direct_mixing_nr)
          CALL scf_env_density_mixing(scf_env%p_mix_new,scf_env%mixing_store,&
                                         rho_ao_kp,para_env,scf_env%iter_delta,&
                                         scf_env%iter_count,diis=diis_step,&
                                         invert=.TRUE.)
          DO ic=1,nc
             DO ispin=1,dft_control%nspins
                CALL cp_dbcsr_copy(rho_ao_kp(ispin,ic)%matrix,scf_env%p_mix_new(ispin,ic)%matrix,&
                                   name=TRIM(cp_dbcsr_name(rho_ao_kp(ispin,ic)%matrix)))
             END DO
          END DO
        CASE(gspace_mixing_nr,pulay_mixing_nr,broyden_mixing_nr,&
             broyden_mixing_new_nr,multisecant_mixing_nr)
          DO ic=1,nc
             DO ispin=1,dft_control%nspins
                CALL cp_dbcsr_copy(rho_ao_kp(ispin,ic)%matrix,scf_env%p_mix_new(ispin,ic)%matrix,&
                                name=TRIM(cp_dbcsr_name(rho_ao_kp(ispin,ic)%matrix)))
             END DO
          END DO
        END SELECT
     ENDIF
  END SUBROUTINE qs_scf_undo_mixing


! *****************************************************************************
!> \brief Performes the updates rho (takes care of mixing as well)
!> \param rho ...
!> \param qs_env ...
!> \param scf_env ...
!> \param ks_env ...
!> \param mix_rho ...
! *****************************************************************************
  SUBROUTINE qs_scf_rho_update(rho,qs_env,scf_env,ks_env,mix_rho)
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    LOGICAL, INTENT(IN)                      :: mix_rho

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_rho_update', &
      routineP = moduleN//':'//routineN

    TYPE(cp_para_env_type), POINTER          :: para_env

    NULLIFY(para_env)
    CALL get_qs_env(qs_env, para_env=para_env)
    ! ** update qs_env%rho
    CALL qs_rho_update_rho(rho, qs_env=qs_env)
    ! ** Density mixing through density matrix or on the reciprocal space grid (exclusive)
    IF(mix_rho) THEN
      CALL gspace_mixing(qs_env, scf_env%mixing_method, scf_env%mixing_store, rho, &
                         para_env, scf_env%iter_count)

    END IF
    CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.)

  END SUBROUTINE qs_scf_rho_update


! *****************************************************************************
!> \brief Performes the necessary steps before leaving innner scf loop
!> \param scf_env ...
!> \param qs_env ...
!> \param diis_step ...
!> \param output_unit ...
! *****************************************************************************
  SUBROUTINE qs_scf_inner_finalize(scf_env,qs_env,diis_step,output_unit)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL                                  :: diis_step
    INTEGER, INTENT(IN)                      :: output_unit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_inner_finalize', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: do_kpoints, failure
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho

    failure=.FALSE.
    NULLIFY(energy, rho, dft_control, ks_env)

    CALL get_qs_env(qs_env=qs_env, energy=energy,ks_env=ks_env,&
                    rho=rho, dft_control=dft_control,para_env=para_env,&
                    do_kpoints=do_kpoints)

    CALL cleanup_scf_loop(scf_env)

    ! now, print out energies and charges corresponding to the obtained wfn
    ! (this actually is not 100% consistent at this point)!
    CALL qs_scf_print_summary(output_unit,qs_env)

    CALL qs_scf_undo_mixing (scf_env,rho,dft_control,para_env,diis_step)

    !   *** update rspace rho since the mo changed
    !   *** this might not always be needed (i.e. no post calculation / no forces )
    !   *** but guarantees that rho and wfn are consistent at this point
    CALL qs_scf_rho_update(rho,qs_env,scf_env,ks_env,mix_rho=.FALSE.)

  END SUBROUTINE qs_scf_inner_finalize


! *****************************************************************************
!> \brief perform cleanup operations at the end of an scf loop
!> \param scf_env ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE cleanup_scf_loop(scf_env)
    TYPE(qs_scf_env_type), POINTER           :: scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'cleanup_scf_loop', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure=.FALSE.

    CPASSERT(ASSOCIATED(scf_env))
    CPASSERT(scf_env%ref_count>0)
! *** method dependent cleanup
       SELECT CASE(scf_env%method)
       CASE(ot_method_nr)
          DO ispin=1,SIZE(scf_env%qs_ot_env)
             CALL ot_scf_destroy(scf_env%qs_ot_env(ispin))
          ENDDO
          DEALLOCATE(scf_env%qs_ot_env,stat=stat)
          CPASSERT(stat==0)
       CASE(ot_diag_method_nr)
   !
       CASE(general_diag_method_nr)
   !
       CASE(special_diag_method_nr)
   !
       CASE(block_krylov_diag_method_nr,block_davidson_diag_method_nr)
   !
       CASE(filter_matrix_diag_method_nr)
   !
       CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown scf method method:"//&
               cp_to_string(scf_env%method),failure)
       END SELECT


    CALL timestop(handle)

  END SUBROUTINE cleanup_scf_loop

END MODULE qs_scf_loop_utils
