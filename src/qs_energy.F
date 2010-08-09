!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Perform a QUICKSTEP wavefunction optimization (single point)
!> \par History
!>      none
!> \author MK (29.10.2002)
! *****************************************************************************
MODULE qs_energy
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_control_utils,                ONLY: read_becke_section,&
                                             read_ddapc_section
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_obj_p_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_to_fm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE efield_utils,                    ONLY: calculate_ecore_efield
  USE et_coupling,                     ONLY: calc_et_coupling
  USE f77_blas
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_diag,&
                                             ehrenfest,&
                                             use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_get
  USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                             pw_pool_type
  USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_dftb_dispersion,              ONLY: calculate_dftb_dispersion
  USE qs_dftb_matrices,                ONLY: build_dftb_matrices
  USE qs_dispersion_pairpot,           ONLY: calculate_dispersion_pairpot
  USE qs_environment_methods,          ONLY: qs_env_update_s_mstruct
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_external_potential,           ONLY: external_c_potential,&
                                             external_e_potential
  USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
  USE qs_ks_scp_methods,               ONLY: qs_ks_scp_update
  USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
  USE qs_overlap,                      ONLY: build_overlap_matrix
  USE qs_scf,                          ONLY: scf
  USE rt_propagation_methods,          ONLY: propagation_step,&
                                             put_data_to_history,&
                                             s_matrices_create,&
                                             update_core_and_matrices
  USE rt_propagation_output,           ONLY: rt_prop_output
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE rt_propagation_utils,            ONLY: calc_S_derivs
  USE rt_propagator_init,              ONLY: init_emd_propagators
  USE scp_dispersion,                  ONLY: scp_nddo_dispersion,&
                                             scp_qs_dispersion
  USE se_core_core,                    ONLY: se_core_core_interaction
  USE se_core_matrix,                  ONLY: build_se_core_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_methods,                     ONLY: xas
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy'

  PUBLIC :: qs_energies

CONTAINS

! *****************************************************************************
!> \brief   Driver routine for QUICKSTEP single point wavefunction optimization.
!> \author  MK
!> \date    29.10.2002
!> \par History
!>          - consistent_energies option added (25.08.2005, TdK)
!>          - introduced driver for energy in order to properly decide between 
!>            SCF or RTP (fschiff 02.09)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE qs_energies (qs_env, globenv, consistent_energies, calc_forces, &
                          error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(IN), OPTIONAL            :: consistent_energies, &
                                                calc_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: run_rtp

    CALL get_qs_env(qs_env=qs_env,run_rtp=run_rtp,error=error)
    IF(run_rtp)THEN
       CALL qs_energies_rtp(qs_env,globenv,error=error)
    ELSE
       CALL qs_energies_scf(qs_env, globenv, consistent_energies, calc_forces, &
                          error)
    END IF

  END SUBROUTINE qs_energies
    

! *****************************************************************************
!> \brief   QUICKSTEP single point wavefunction optimization.
!> \author  MK
!> \date    29.10.2002
!> \par History
!>          - consistent_energies option added (25.08.2005, TdK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE qs_energies_scf (qs_env, globenv, consistent_energies, calc_forces, &
                          error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(IN), OPTIONAL            :: consistent_energies, &
                                                calc_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies_scf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nder
    LOGICAL                                  :: do_et, my_calc_forces
    TYPE(cp_dbcsr_obj_p_type), &
      DIMENSION(:), POINTER                  :: matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(section_vals_type), POINTER         :: input, rest_b_section

    CALL timeset(routineN,handle)

    my_calc_forces = .FALSE.
    IF( PRESENT( calc_forces ) ) my_calc_forces = calc_forces
    para_env=>qs_env%para_env

    CALL get_qs_env(qs_env=qs_env,input=input, dft_control=dft_control, &
                    distribution_2d=distribution_2d, error=error)
    
    CALL build_qs_neighbor_lists(qs_env,para_env,force_env_section=input,error=error)
    qs_env%dft_control%qs_control%becke_control%need_pot=.TRUE.

    ! *** Calculate the overlap and the core Hamiltonian integral matrix ***
    IF ( dft_control%qs_control%semi_empirical ) THEN
       CALL build_se_core_matrix(qs_env=qs_env, para_env=para_env,&
                                 calculate_forces=.FALSE.,error=error)
       CALL qs_env_update_s_mstruct(qs_env,error=error)
       CALL se_core_core_interaction(qs_env, para_env, calculate_forces=.FALSE., error=error)
       IF ( dft_control%qs_control%se_control%scp ) THEN
          CALL scp_nddo_dispersion (qs_env,calc_forces,error)
       END IF
    ELSEIF ( dft_control%qs_control%dftb ) THEN
       CALL build_dftb_matrices(qs_env=qs_env, para_env=para_env,&
                              calculate_forces=.FALSE.,error=error)
       CALL calculate_dftb_dispersion(qs_env=qs_env, para_env=para_env,&
                              calculate_forces=.FALSE.,error=error)
       CALL qs_env_update_s_mstruct(qs_env,error=error)
    ELSE
       CALL build_core_hamiltonian_matrix(qs_env=qs_env, globenv=globenv,&
                                          calculate_forces=.FALSE.,error=error)
       ! *** In case of denisty_fitting, calculate the corresponding overlap_matrices
       IF( dft_control%do_admm) THEN
         NULLIFY(matrix_s_aux_fit,matrix_s_aux_fit_vs_orb)
         IF( my_calc_forces ) THEN
           nder = 1
         ELSE
           nder = 0
         END IF
         CALL get_qs_env(qs_env=qs_env, matrix_s_aux_fit=matrix_s_aux_fit, error=error)
         CALL build_overlap_matrix(qs_env,para_env,nderivative=nder,matrix_s=matrix_s_aux_fit,&
                                   matrix_name="AUX_FIT_OVERLAP",&
                                   basis_set_id_a=use_aux_fit_basis_set,&
                                   basis_set_id_b=use_aux_fit_basis_set, &
                                   neighbor_list_sab=qs_env%sab_aux_fit,&
                                   error=error)
         CALL set_qs_env(qs_env=qs_env,matrix_s_aux_fit=matrix_s_aux_fit,error=error)
         CALL get_qs_env(qs_env=qs_env, matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb, error=error)
         CALL build_overlap_matrix(qs_env,para_env,nderivative=nder,matrix_s=matrix_s_aux_fit_vs_orb,&
                                   matrix_name="MIXED_OVERLAP",&
                                   basis_set_id_a=use_aux_fit_basis_set,&
                                   basis_set_id_b=use_orb_basis_set, &
                                   neighbor_list_sab=qs_env%sab_aux_fit_vs_orb,&
                                   error=error)
         CALL set_qs_env(qs_env=qs_env,matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb,error=error)
       END IF

       CALL qs_env_update_s_mstruct(qs_env,error=error)
       CALL calculate_ecore_self(qs_env,error=error)
       CALL calculate_ecore_efield(qs_env,calculate_forces=.FALSE.,error=error)
       CALL calculate_ecore_overlap(qs_env, para_env, calculate_forces=.FALSE.,error=error)
       CALL external_c_potential(qs_env,calculate_forces=.FALSE.,error=error)
       CALL external_e_potential(qs_env,error=error)

       IF ( dft_control%scp ) THEN
          IF(dft_control%scp_control%dispersion ) THEN
             CALL scp_qs_dispersion ( qs_env, calc_forces, error )
          END IF
       ELSE
          ! Add possible pair potential dispersion energy - Evaluate first so we can print
          ! energy info at the end of the SCF
          CALL calculate_dispersion_pairpot(qs_env,calc_forces,error)
       END IF

    END IF

    ! *** Perform a SCF run ***
    CALL scf(qs_env=qs_env,  error=error)

    IF (PRESENT(consistent_energies)) THEN
      IF (consistent_energies) THEN
        IF ( dft_control%scp ) THEN
           CALL qs_ks_scp_update ( qs_env, just_energy=.TRUE., error=error )
        END IF
        CALL qs_ks_update_qs_env(ks_env=qs_env%ks_env, &
                                 qs_env=qs_env, &
                                 calculate_forces=.FALSE., &
                                 just_energy=.TRUE.,error=error)
      END IF
    END IF

    ! **********  Calculate the electron transfer coupling elements********
    do_et=.FALSE.
    do_et=qs_env%dft_control%qs_control%et_coupling_calc
    IF(do_et)THEN
       qs_env%et_coupling%energy=qs_env%energy%total
       qs_env%et_coupling%keep_matrix=.TRUE.
       qs_env%et_coupling%first_run=.TRUE.
       CALL qs_ks_update_qs_env(ks_env=qs_env%ks_env, &
                                qs_env=qs_env, &
                                calculate_forces=.FALSE., &
                                just_energy=.TRUE.,error=error)
       qs_env%et_coupling%first_run=.FALSE.
       IF(qs_env%dft_control%qs_control%ddapc_restraint)THEN
          rest_b_section =>  section_vals_get_subs_vals(input,"PROPERTIES%ET_COUPLING%DDAPC_RESTRAINT_B",&
                                                        error=error)
          CALL read_ddapc_section(qs_control=dft_control%qs_control,&
                                  ddapc_restraint_section=rest_b_section,error=error)
       END IF
       IF(qs_env%dft_control%qs_control%becke_restraint)THEN
          rest_b_section => section_vals_get_subs_vals(input,"PROPERTIES%ET_COUPLING%BECKE_RESTRAINT_B",&
                                                       error=error)
          CALL read_becke_section(qs_control=dft_control%qs_control,&
                                     becke_section=rest_b_section,error=error)
       END IF
       CALL scf(qs_env=qs_env, error=error)
       qs_env%et_coupling%keep_matrix=.TRUE.

       CALL qs_ks_update_qs_env(ks_env=qs_env%ks_env, &
                                qs_env=qs_env, &
                                calculate_forces=.FALSE., &
                                just_energy=.TRUE.,error=error) 
       CALL calc_et_coupling(qs_env,error)
       IF(qs_env%dft_control%qs_control%becke_restraint)THEN
          CALL pw_env_get(qs_env%pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
          CALL pw_pool_give_back_pw(auxbas_pw_pool,&
               qs_env%dft_control%qs_control%becke_control%becke_pot%pw,error=error)
          qs_env%dft_control%qs_control%becke_control%need_pot=.TRUE.
       END IF
    END IF
    
    !Properties
    IF(dft_control%do_xas_calculation) THEN
      CALL xas(qs_env, dft_control, error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_energies_scf

! *****************************************************************************
!> \brief   energy calculation for real time propagation (iteratet through integrator)
!> \author  Florian Schiffmann
!> \date    02.2009
!> \version 1.0
! *****************************************************************************
  SUBROUTINE qs_energies_rtp(qs_env,globenv,error)
    
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies_rtp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: aspc_order, handle, i, ihist
    LOGICAL                                  :: diagonalize
    REAL(Kind=dp)                            :: delta_iter, eps_ener
    TYPE(cp_dbcsr_obj_p_type), &
      DIMENSION(:), POINTER                  :: s_mat
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new, exp_H_old, &
                                                mos_new, mos_old
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL timeset(routineN,handle)

    eps_ener=qs_env%dft_control%rtp_control%eps_ener
    aspc_order=qs_env%dft_control%rtp_control%aspc_order
    IF(qs_env%dft_control%rtp_control%initial_step)THEN
       CALL init_emd_propagators(qs_env,error)
    ELSE
       IF(qs_env%rtp%iter==1)THEN
          CALL update_core_and_matrices(qs_env,globenv,error)
          CALL calculate_ecore_efield(qs_env,calculate_forces=.FALSE.,error=error)
          CALL external_c_potential(qs_env,calculate_forces=.FALSE.,error=error)
          CALL external_e_potential(qs_env,error=error)
       END IF
       CALL get_qs_env(qs_env=qs_env,rtp=rtp,matrix_s=s_mat,error=error)

       diagonalize=(qs_env%dft_control%rtp_control%mat_exp==do_diag)
       CALL s_matrices_create (s_mat,rtp,diagonalize,error)
       CALL calc_S_derivs(qs_env,error)

       CALL get_rtp(rtp=rtp,&
                 exp_H_old=exp_H_old,&
                 exp_H_new=exp_H_new,&
                 mos_old=mos_old,&
                 mos_new=mos_new,&
                 error=error)
       
       eps_ener=qs_env%dft_control%rtp_control%eps_ener
       
       CALL propagation_step(qs_env,delta_iter, error)
       rtp%energy_new=qs_env%energy%total
       rtp%converged=(delta_iter.LT.eps_ener)

       IF(rtp%converged)THEN
          ihist=MOD(rtp%istep,aspc_order)+1
          DO i=1,SIZE(mos_new)
             CALL cp_fm_to_fm(mos_new(i)%matrix,mos_old(i)%matrix,error)
             CALL cp_fm_to_fm(exp_H_new(i)%matrix,exp_H_old(i)%matrix,error)
          END DO
          IF(rtp%history%mos_or_H==2)&
               CALL  put_data_to_history(qs_env,exp_H=exp_H_new,ihist=ihist,error=error)
          IF(rtp%history%mos_or_H==1)&
               CALL  put_data_to_history(qs_env,mos=mos_new,s_mat=s_mat,ihist=ihist,error=error)
       END IF
       CALL rt_prop_output(qs_env,ehrenfest,delta_iter=delta_iter,error=error)
       rtp%energy_old=rtp%energy_new       
    END IF
    CALL timestop(handle)
  END SUBROUTINE qs_energies_rtp


END MODULE qs_energy
