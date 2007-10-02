!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/init_extended_system_types [1.0] *
!!
!!   NAME
!!     init_extended_system_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     cjm, FEB 20 2001:  added subroutine initialize_extended_parameters
!!     cjm, MAY 03 2001:  reorganized and added separtate routines for
!!                        nhc_part, nhc_baro, nhc_ao, npt
!!   SOURCE
!******************************************************************************

MODULE init_extended_system_types

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE global_types,                    ONLY: global_environment_type
  USE init_extended_system_variables,  ONLY: init_barostat_variables,&
                                             init_nhc_variables
  USE input_constants,                 ONLY: nph_uniaxial_damped_ensemble,&
                                             nph_uniaxial_ensemble,&
                                             npt_f_ensemble,&
                                             npt_i_ensemble
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_remove_values,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: global_constraint_type,&
                                             molecule_type
  USE nhc_mapping,                     ONLY: nhc_to_barostat_mapping,&
                                             nhc_to_coef_mapping,&
                                             nhc_to_particle_mapping,&
                                             nhc_to_shell_mapping
  USE termination,                     ONLY: stop_program
  USE thermostat_types,                ONLY: thermostats_info_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PUBLIC :: initialize_nhc_part, initialize_nhc_baro, initialize_npt, &
            initialize_nhc_coef, initialize_nhc_shell
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'init_extended_system_types'

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** init_extended_system_types/initialize_npt [1.0] *
!!
!!   NAME
!!     initialize_npt_type
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE initialize_npt ( simpar, para_env, globenv, npt_info, box, work_section, error)

    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(npt_info_type), DIMENSION(:, :), &
      POINTER                                :: npt_info
    TYPE(cell_type), INTENT(inout)           :: box
    TYPE(section_vals_type), POINTER         :: work_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_npt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ind, j, stat
    LOGICAL                                  :: explicit, failure, restart
    REAL(KIND=dp)                            :: temp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: buffer
    TYPE(section_vals_type), POINTER         :: work_section2

    NULLIFY(work_section2)
    explicit = .FALSE.
    failure=.FALSE.
    restart = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(npt_info),cp_fatal_level,routineP,error,failure)
    ! first allocating the npt_info_type if requested
    SELECT CASE ( simpar % ensemble )
    CASE ( npt_i_ensemble )
       ALLOCATE ( npt_info ( 1, 1 ), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       npt_info  ( :, : ) % eps = LOG ( box % deth ) / 3.0_dp
       temp = simpar % temp_ext

    CASE ( npt_f_ensemble )
       ALLOCATE ( npt_info ( 3, 3 ), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       temp = simpar % temp_ext

    CASE ( nph_uniaxial_ensemble )
       ALLOCATE ( npt_info ( 1, 1 ), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       temp = simpar % temp_ext

    CASE (  nph_uniaxial_damped_ensemble )
       ALLOCATE ( npt_info ( 1, 1 ), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       temp = simpar % temp_ext

    CASE  DEFAULT
       ! Do nothing..
       NULLIFY(npt_info)
    END SELECT

    IF (ASSOCIATED(npt_info)) THEN
       IF (ASSOCIATED(work_section)) THEN
          work_section2 => section_vals_get_subs_vals(work_section,"VELOCITY",error=error)
          CALL section_vals_get(work_section2, explicit=explicit, error=error)
          restart=explicit
          work_section2 => section_vals_get_subs_vals(work_section,"MASS",error=error)
          CALL section_vals_get(work_section2, explicit=explicit, error=error)
          CALL cp_assert(restart.EQV.explicit,cp_failure_level,cp_assertion_failed,routineP,&
               "You need to define both VELOCITY and MASS section (or none) in the BAROSTAT section",&
               error=error, failure=failure)
          restart=explicit.AND.restart
       END IF
       

       IF ( restart ) THEN
          CALL section_vals_val_get(work_section,"VELOCITY%_DEFAULT_KEYWORD_",r_vals=buffer,error=error)
          ind = 0
          DO i = 1, SIZE(npt_info,1)
             DO j = 1, SIZE(npt_info,2)
                ind = ind + 1
                npt_info ( i, j ) % v = buffer(ind)
             ENDDO
          ENDDO
          CALL section_vals_val_get(work_section,"MASS%_DEFAULT_KEYWORD_",r_vals=buffer,error=error)
          ind = 0
          DO i = 1, SIZE(npt_info,1)
             DO j = 1, SIZE(npt_info,2)
                ind = ind + 1
                npt_info ( i, j ) % mass = buffer(ind)
             ENDDO
          ENDDO
       ELSE
          CALL init_barostat_variables ( npt_info, simpar % tau_cell, temp,  &
               box % deth, simpar % nfree, simpar % ensemble, simpar % cmass,&
               para_env, globenv ,error=error)
       END IF

    END IF
  END SUBROUTINE initialize_npt

!******************************************************************************
!!****** initialize_extended_types/initialize_nhc_baro [1.0] *
!!
!!   NAME
!!     initialize_nhc_baro
!!
!!   FUNCTION
!!     fire up the thermostats, if NPT
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE initialize_nhc_baro ( simpar, para_env, globenv, nhc, nose_section, save_mem, error )

    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    TYPE(section_vals_type), POINTER         :: nose_section
    LOGICAL, INTENT(IN)                      :: save_mem
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_nhc_baro', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure, restart
    REAL(KIND=dp)                            :: temp

    failure = .FALSE.
    restart = .FALSE.
    CALL nhc_to_barostat_mapping ( simpar, nhc, para_env)
    ! Set up the Yoshida weights
    IF ( nhc % nyosh > 0 ) THEN
       ALLOCATE ( nhc % dt_yosh ( 1 : nhc%nyosh ), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       CALL set_yoshida_coef ( nhc, simpar % dt )
    ENDIF

     CALL restart_nose(nhc,nose_section,save_mem,restart,error)

     IF(.NOT. restart) THEN
       ! Initializing thermostat forces and velocities for the Nose-Hoover
       ! Chain variables
       SELECT CASE (simpar % ensemble )
       CASE DEFAULT
          temp = simpar % temp_ext
       END SELECT

       IF ( nhc % nhc_len /= 0 ) &
            CALL init_nhc_variables ( nhc, simpar % tau_nhc,  &
            temp, para_env, globenv ,error=error)
     END IF

  END SUBROUTINE initialize_nhc_baro

!******************************************************************************
!!****** initialize_extended_types/initialize_nhc_part [1.0] *
!!
!!   NAME
!!     initialize_nhc_part
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE initialize_nhc_part ( thermostats_info, box, simpar, local_molecules,&
       molecule, molecule_kind_set, para_env, globenv, nhc, print_section, nose_section,&
       gci, save_mem, error)

    TYPE(thermostats_info_type), POINTER     :: thermostats_info
    TYPE(cell_type), POINTER                 :: box
    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    TYPE(section_vals_type), POINTER         :: print_section, nose_section
    TYPE(global_constraint_type), POINTER    :: gci
    LOGICAL, INTENT(IN)                      :: save_mem
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_nhc_part', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure, restart

    failure = .FALSE.
    restart = .FALSE.
    ! fire up the thermostats, if not NVE

    CALL nhc_to_particle_mapping ( thermostats_info, simpar, local_molecules,&
         molecule, molecule_kind_set, nhc, para_env )

    ! Set up the Yoshida weights
    IF ( nhc % nyosh > 0 ) THEN
       ALLOCATE ( nhc % dt_yosh ( 1 : nhc%nyosh ), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       CALL set_yoshida_coef ( nhc, simpar % dt )
    ENDIF

    CALL restart_nose(nhc,nose_section,save_mem,restart,error)

    IF(.NOT. restart) THEN
       ! Initializing thermostat forces and velocities for the Nose-Hoover
       ! Chain variables
       IF ( nhc % nhc_len /= 0 ) &
            CALL init_nhc_variables ( nhc, simpar % tau_nhc,  &
            simpar % temp_ext, para_env, globenv ,error=error)

    ENDIF

  END SUBROUTINE initialize_nhc_part

!******************************************************************************
!!****** initialize_extended_types/initialize_nhc_coef [1.0] *
!!
!!   NAME
!!     initialize_nhc_coef
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!   04-02-04 gt: rewritten using new dynamical_coef_type
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE initialize_nhc_coef ( dyn_coeff_set, para_env, globenv, nhc, &
             coeff_section, dt, save_mem, nose_section, error )

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    TYPE(section_vals_type), POINTER         :: coeff_section
    REAL(KIND=dp), INTENT(IN)                :: dt
    LOGICAL, INTENT(IN)                      :: save_mem
    TYPE(section_vals_type), POINTER         :: nose_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_nhc_coef', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure, restart
    REAL(KIND=dp)                            :: kg_coeff_target_temp, tau_nhc

    failure = .FALSE.
    restart = .FALSE.
    CALL section_vals_val_get(nose_section,"timecon",r_val=tau_nhc,error=error)
    CALL section_vals_val_get(coeff_section,"temperature",r_val=kg_coeff_target_temp,error=error)

    ! fire up the thermostats
    CALL nhc_to_coef_mapping ( dyn_coeff_set, nhc, coeff_section, para_env ,error=error)

    IF ( nhc % nhc_len /= 0 ) THEN
       CALL restart_nose(nhc, nose_section, save_mem, restart, error)
       
       IF(.NOT. restart) THEN
          ! Initializing thermostat forces and velocities for the Nose-Hoover
          ! Chain variables
          CALL init_nhc_variables (nhc, tau_nhc, kg_coeff_target_temp, para_env, globenv,error=error)
          ! Set up the Yoshida weights
          IF ( nhc % nyosh > 0 ) THEN
             ALLOCATE ( nhc % dt_yosh ( 1 : nhc%nyosh ), STAT = stat )
             CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
             CALL set_yoshida_coef ( nhc, dt )
          ENDIF
       END IF
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "nve for coefs implemented but not tested")
    END IF

  END SUBROUTINE initialize_nhc_coef

!******************************************************************************
!!****** initialize_extended_types/initialize_nhc_shell [1.0] *
!!
!!   NAME
!!     initialize_nhc_shell
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE initialize_nhc_shell(simpar, para_env, globenv, nhc, atomic_kind_set, &
             local_particles, nose_section, save_mem, error)
 
    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(section_vals_type), POINTER         :: nose_section
    LOGICAL, INTENT(IN)                      :: save_mem
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_nhc_shell', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure, restart

    failure = .FALSE.
    CALL nhc_to_shell_mapping(simpar, nhc, atomic_kind_set, local_particles, para_env)

    IF(simpar%shell_nose) THEN
       restart = .FALSE.
       ! Set up the Yoshida weights
       IF ( nhc % nyosh > 0 ) THEN
          ALLOCATE ( nhc % dt_yosh ( 1 : nhc%nyosh ), STAT = stat )
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
          CALL set_yoshida_coef ( nhc, simpar % dt )
       ENDIF
       
       CALL restart_nose(nhc, nose_section, save_mem, restart, error)
       
       IF(.NOT. restart) THEN
          ! Initialize thermostat forces and velocities
          ! Chain variables
          IF ( nhc % nhc_len /= 0 ) THEN
             CALL init_nhc_variables (nhc, simpar%shell_tau_nhc, simpar%temp_sh_ext, &
                  para_env, globenv, error=error)
          END IF
       END IF
    END IF
  END SUBROUTINE  initialize_nhc_shell
  
!!*****
!******************************************************************************
!!****** extended_system_initialize/set_yoshida_coef [1.0] *
!!
!!   NAME
!!     set_yoshida_coef
!!
!!   FUNCTION
!!     This lists the coefficients for the Yoshida method (higher
!!     order integrator used in NVT)
!!
!!   AUTHOR
!!
!!   CREATION DATE
!!     14-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE set_yoshida_coef ( nhc, dt )

    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    REAL(KIND=dp), INTENT(IN)                :: dt

    CHARACTER(len=*), PARAMETER :: routineN = 'set_yoshida_coef', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp), DIMENSION(nhc%nyosh)      :: yosh_wt

    SELECT CASE (nhc % nyosh)
    CASE DEFAULT
       CALL stop_program ( 'yoshida_coef', 'Value not available' )
    CASE (1)
       yosh_wt(1) = 1.0_dp
    CASE (3)
       yosh_wt(1) = 1.0_dp/(2.0_dp-(2.0_dp)**(1.0_dp/3.0_dp))
       yosh_wt(2) = 1.0_dp - 2.0_dp*yosh_wt(1)
       yosh_wt(3) = yosh_wt(1)
    CASE (5)
       yosh_wt(1) = 1.0_dp/(4.0_dp-(4.0_dp)**(1.0_dp/3.0_dp))
       yosh_wt(2) = yosh_wt(1)
       yosh_wt(4) = yosh_wt(1)
       yosh_wt(5) = yosh_wt(1)
       yosh_wt(3) = 1.0_dp - 4.0_dp*yosh_wt(1)
    CASE (7)
       yosh_wt(1) = .78451361047756_dp
       yosh_wt(2) = .235573213359357_dp
       yosh_wt(3) = -1.17767998417887_dp
       yosh_wt(4) = 1.0_dp - 2.0_dp*(yosh_wt(1)+yosh_wt(2)+yosh_wt(3))
       yosh_wt(5) = yosh_wt(3)
       yosh_wt(6) = yosh_wt(2)
       yosh_wt(7) = yosh_wt(1)
    CASE (9)
       yosh_wt(1) = 0.192_dp
       yosh_wt(2) = 0.554910818409783619692725006662999_dp
       yosh_wt(3) = 0.124659619941888644216504240951585_dp
       yosh_wt(4) = -0.843182063596933505315033808282941_dp
       yosh_wt(5) = 1.0_dp - 2.0_dp*(yosh_wt(1)+yosh_wt(2)+&
            yosh_wt(3)+yosh_wt(4))
       yosh_wt(6) = yosh_wt(4)
       yosh_wt(7) = yosh_wt(3)
       yosh_wt(8) = yosh_wt(2)
       yosh_wt(9) = yosh_wt(1)
    CASE (15)
       yosh_wt(1) = 0.102799849391985_dp
       yosh_wt(2) = -0.196061023297549e1_dp
       yosh_wt(3) = 0.193813913762276e1_dp
       yosh_wt(4) = -0.158240635368243_dp
       yosh_wt(5) = -0.144485223686048e1_dp
       yosh_wt(6) = 0.253693336566229_dp
       yosh_wt(7) = 0.914844246229740_dp
       yosh_wt(8) = 1.0_dp - 2.0_dp*(yosh_wt(1)+yosh_wt(2)+&
            yosh_wt(3)+yosh_wt(4)+yosh_wt(5)+yosh_wt(6)+yosh_wt(7))
       yosh_wt(9) = yosh_wt(7)
       yosh_wt(10) = yosh_wt(6)
       yosh_wt(11) = yosh_wt(5)
       yosh_wt(12) = yosh_wt(4)
       yosh_wt(13) = yosh_wt(3)
       yosh_wt(14) = yosh_wt(2)
       yosh_wt(15) = yosh_wt(1)
    END SELECT
    nhc % dt_yosh = dt * yosh_wt / REAL ( nhc % nc,KIND=dp)

  END SUBROUTINE set_yoshida_coef

!******************************************************************************
!!****** initialize_extended_types/restart_nose [1.0] *
!!
!!   NAME
!!     restart_nose 
!!
!!   FUNCTION
!!     read coordinate, velocities, forces and masses of the
!!     thermostat from restart file
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!    24-07-07 created 
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE restart_nose(nhc,nose_section,save_mem,restart,error)

    TYPE(lnhc_parameters_type), &
      INTENT(INOUT)                          :: nhc
    TYPE(section_vals_type), POINTER         :: nose_section
    LOGICAL, INTENT(IN)                      :: save_mem
    LOGICAL, INTENT(OUT)                     :: restart
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'restart_nose', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ind, j
    LOGICAL                                  :: explicit
    REAL(KIND=dp), DIMENSION(:), POINTER     :: buffer
    TYPE(section_vals_type), POINTER         :: work_section

    NULLIFY(buffer, work_section)
    explicit = .FALSE.
    restart = .FALSE.
    
    IF (ASSOCIATED(nose_section)) THEN
       work_section => section_vals_get_subs_vals(nose_section,"VELOCITY",error=error)
       CALL section_vals_get(work_section, explicit=explicit, error=error)
       restart=explicit
       work_section => section_vals_get_subs_vals(nose_section,"COORD",error=error)
       CALL section_vals_get(work_section, explicit=explicit, error=error)
       CALL cp_assert(restart.or..not.explicit,cp_failure_level,cp_assertion_failed,routineP,&
            "You need to define both VELOCITY and COORD and MASS and FORCE section (or none) in the NOSE section",&
            error=error)
       restart=explicit.and.restart
       work_section => section_vals_get_subs_vals(nose_section,"MASS",error=error)
       CALL section_vals_get(work_section, explicit=explicit, error=error)
       CALL cp_assert(restart.or..not.explicit,cp_failure_level,cp_assertion_failed,routineP,&
            "You need to define both VELOCITY and COORD and MASS and FORCE section (or none) in the NOSE section",&
            error=error)
       restart=explicit.and.restart
       work_section => section_vals_get_subs_vals(nose_section,"FORCE",error=error)
       CALL section_vals_get(work_section, explicit=explicit, error=error)
       CALL cp_assert(restart.or..not.explicit,cp_failure_level,cp_assertion_failed,routineP,&
            "You need to define both VELOCITY and COORD and MASS and FORCE section (or none) in the NOSE section",&
            error=error)
       restart=explicit.and.restart
    END IF

    IF (restart) THEN
       CALL section_vals_val_get(nose_section,"COORD%_DEFAULT_KEYWORD_",r_vals=buffer,error=error)
       DO i = 1, SIZE ( nhc % nvt, 2)
          ind = nhc % INDEX(i)
          ind = (ind-1) * nhc % nhc_len
          DO j = 1, SIZE(nhc % nvt, 1)
             ind = ind + 1
             nhc % nvt(j,i) % eta = buffer ( ind )
          END DO
       END DO
       CALL section_vals_val_get(nose_section,"VELOCITY%_DEFAULT_KEYWORD_",r_vals=buffer,error=error)
       DO i = 1, SIZE ( nhc % nvt, 2)
          ind = nhc % INDEX(i)
          ind = (ind-1) * nhc % nhc_len
          DO j = 1, SIZE(nhc % nvt, 1)
             ind = ind + 1
             nhc % nvt(j,i) % v = buffer ( ind )
          END DO
       END DO
       CALL section_vals_val_get(nose_section,"MASS%_DEFAULT_KEYWORD_",r_vals=buffer,error=error)
       DO i = 1, SIZE ( nhc % nvt, 2)
          ind = nhc % INDEX(i)
          ind = (ind-1) * nhc % nhc_len
          DO j = 1, SIZE(nhc % nvt, 1)
             ind = ind + 1
             nhc % nvt(j,i) % mass = buffer ( ind )
          END DO
       END DO
       CALL section_vals_val_get(nose_section,"FORCE%_DEFAULT_KEYWORD_",r_vals=buffer,error=error)
       DO i = 1, SIZE ( nhc % nvt, 2)
          ind = nhc % INDEX(i)
          ind = (ind-1) * nhc % nhc_len
          DO j = 1, SIZE(nhc % nvt, 1)
             ind = ind + 1
             nhc % nvt(j,i) % f = buffer ( ind )
          END DO
       END DO

       IF(save_mem) THEN
         NULLIFY(work_section)
         work_section => section_vals_get_subs_vals(nose_section,"COORD",error=error)
         CALL section_vals_remove_values(work_section, error)
         NULLIFY(work_section)
         work_section => section_vals_get_subs_vals(nose_section,"VELOCITY",error=error)
         CALL section_vals_remove_values(work_section, error)
         NULLIFY(work_section)
         work_section => section_vals_get_subs_vals(nose_section,"FORCE",error=error)
         CALL section_vals_remove_values(work_section, error)
         NULLIFY(work_section)
         work_section => section_vals_get_subs_vals(nose_section,"MASS",error=error)
         CALL section_vals_remove_values(work_section, error)
       END IF
     ENDIF
  END SUBROUTINE restart_nose

!!*****
!******************************************************************************

END MODULE init_extended_system_types

!******************************************************************************
