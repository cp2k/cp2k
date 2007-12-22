!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_energy [1.0] *
!!
!!   NAME
!!     qs_energy
!!
!!   FUNCTION
!!     Perform a QUICKSTEP wavefunction optimization (single point)
!!
!!   AUTHOR
!!     MK (29.10.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_energy
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_control_utils,                ONLY: read_becke_section,&
                                             read_ddapc_section
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE et_coupling,                     ONLY: calc_et_coupling
  USE f77_blas
  USE global_types,                    ONLY: global_environment_type
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
  USE qs_environment_methods,          ONLY: qs_env_update_s_mstruct
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
  USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
  USE qs_scf,                          ONLY: scf
  USE qs_semi_empirical_hamiltonian,   ONLY: build_se_core_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_methods,                     ONLY: xas
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy'

  PUBLIC :: qs_energies

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE qs_energies (qs_env, globenv, consistent_energies, calc_forces, &
                          error)
!   Purpose: QUICKSTEP single point wavefunction optimization.

!   History: - Creation (29.10.2002,MK)
!            - consistent_energies option added (25.08.2005, TdK)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(IN), OPTIONAL            :: consistent_energies, &
                                                calc_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_energies', &
      routine = "SUBROUTINE "//routineN//" (MODULE "//moduleN//")"

    INTEGER                                  :: handle
    LOGICAL                                  :: do_et, right_iteration_level
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(section_vals_type), POINTER         :: input, rest_b_section

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    para_env=>qs_env%para_env

    CALL get_qs_env(qs_env=qs_env,input=input, dft_control=dft_control, &
                    error=error)

    CALL build_qs_neighbor_lists(qs_env,para_env,force_env_section=input,error=error)
    qs_env%dft_control%qs_control%becke_control%need_pot=.TRUE.

    ! *** Calculate the overlap and the core Hamiltonian integral matrix ***
    IF ( dft_control%qs_control%semi_empirical ) THEN
       CALL build_se_core_matrix(qs_env=qs_env, para_env=para_env,&
                                 calculate_forces=.FALSE.,error=error)
       CALL qs_env_update_s_mstruct(qs_env,error=error)
    ELSEIF ( dft_control%qs_control%dftb ) THEN
       CALL build_dftb_matrices(qs_env=qs_env, para_env=para_env,&
                              calculate_forces=.FALSE.,error=error)
       CALL calculate_dftb_dispersion(qs_env=qs_env, para_env=para_env,&
                              calculate_forces=.FALSE.,error=error)
       CALL qs_env_update_s_mstruct(qs_env,error=error)
    ELSE
       CALL build_core_hamiltonian_matrix(qs_env=qs_env, globenv=globenv,&
                                          calculate_forces=.FALSE.,error=error)
       CALL qs_env_update_s_mstruct(qs_env,error=error)
       CALL calculate_ecore_self(qs_env,error=error)
       CALL calculate_ecore_overlap(qs_env, para_env, &
                                      calculate_forces=.FALSE.,error=error)
    END IF

    ! *** Perform a SCF run ***
    CALL scf(qs_env=qs_env, globenv=globenv, error=error)

    IF (PRESENT(consistent_energies)) THEN
      IF (consistent_energies) THEN
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
       IF(qs_env%dft_control%qs_control%becke_constraint)THEN
          rest_b_section => section_vals_get_subs_vals(input,"PROPERTIES%ET_COUPLING%BECKE_CONSTRAINT_B",&
                                                       error=error)
          CALL read_becke_section(qs_control=dft_control%qs_control,&
                                     becke_section=rest_b_section,error=error)
       END IF
       CALL scf(qs_env=qs_env, globenv=globenv, error=error)
       qs_env%et_coupling%keep_matrix=.TRUE.

       CALL qs_ks_update_qs_env(ks_env=qs_env%ks_env, &
                                qs_env=qs_env, &
                                calculate_forces=.FALSE., &
                                just_energy=.TRUE.,error=error) 
       CALL calc_et_coupling(qs_env,error)
       IF(qs_env%dft_control%qs_control%becke_constraint)THEN
          CALL pw_env_get(qs_env%pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
          CALL pw_pool_give_back_pw(auxbas_pw_pool,&
               qs_env%dft_control%qs_control%becke_control%becke_pot%pw,error=error)
          qs_env%dft_control%qs_control%becke_control%need_pot=.TRUE.
       END IF
    END IF
    



    !Properties
    right_iteration_level = .TRUE.
    !right_iteration_level = (MODULO(iter_count,dft_control%xas_control%each_step)==0)
    IF(dft_control%do_xas_calculation .AND. right_iteration_level ) THEN
      CALL xas(qs_env, dft_control%xas_control, globenv,error=error)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_energies

! *****************************************************************************

END MODULE qs_energy
