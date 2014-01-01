!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Setup of regions with different temperature
!> \author MI
! *****************************************************************************
MODULE thermal_region_utils

  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_constants,                 ONLY: npt_f_ensemble,&
                                             npt_i_ensemble,&
                                             nvt_ensemble
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE particle_list_types,             ONLY: particle_list_type
  USE simpar_types,                    ONLY: simpar_type
  USE thermal_region_types,            ONLY: allocate_thermal_regions,&
                                             release_thermal_regions,&
                                             thermal_region_type,&
                                             thermal_regions_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: create_thermal_regions, print_thermal_regions

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermal_region_utils'

CONTAINS

! *****************************************************************************
!> \brief create thermal_regions
!> \author
! *****************************************************************************
  SUBROUTINE create_thermal_regions(thermal_regions, md_section, simpar, force_env, error)
    TYPE(thermal_regions_type), POINTER      :: thermal_regions
    TYPE(section_vals_type), POINTER         :: md_section
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_thermal_regions', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, il, ipart, ireg, nlist, &
                                                nregions
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    LOGICAL                                  :: apply_thermostat, explicit, &
                                                failure
    REAL(KIND=dp)                            :: temp, temp_tol
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: region_sections, &
                                                thermal_region_section
    TYPE(thermal_region_type), POINTER       :: t_region

    failure =.FALSE.

    NULLIFY(region_sections,t_region,thermal_region_section, particles,subsys,tmplist)
    CALL allocate_thermal_regions(thermal_regions,error)
    thermal_region_section =>  section_vals_get_subs_vals(md_section,"THERMAL_REGION",error=error)
    CALL section_vals_get(thermal_region_section, explicit=explicit, error=error)
    IF (explicit)THEN
       apply_thermostat =(simpar%ensemble == nvt_ensemble)  .OR.&
                         (simpar%ensemble == npt_f_ensemble).OR.&
                         (simpar%ensemble == npt_i_ensemble)
       IF(apply_thermostat) THEN
         CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
             "With the chosen ensemble the temperature is "//&
             "controlled by thermostats. The definition of different thermal "//&
             "regions might result inconsistent with the presence of thermostats."//&
CPSourceFileRef,&
             only_ionode=.TRUE.)
       END IF
       IF(simpar%temp_tol>0.0_dp) THEN
         CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
             "Control of the global temperature by rescaling of the velocity"//&
             "is not consistent with the presence of different thermal regions."//&
             "The temperature of different regions is rescaled separatedly."//&
CPSourceFileRef,&
             only_ionode=.TRUE.)
       END IF
       CALL section_vals_val_get(thermal_region_section,"FORCE_RESCALING",l_val=thermal_regions%force_rescaling , error=error)
       region_sections =>  section_vals_get_subs_vals(thermal_region_section,"DEFINE_REGION",error=error)
       CALL section_vals_get(region_sections, n_repetition=nregions, error=error)
       IF(nregions>0)THEN
         thermal_regions%nregions = nregions
         thermal_regions%section => thermal_region_section
         ALLOCATE(thermal_regions%thermal_region(nregions))
         CALL force_env_get(force_env, subsys=subsys,error=error)
         CALL cp_subsys_get(subsys, particles=particles,error=error)
         DO  ireg = 1,nregions
           NULLIFY(t_region)
           t_region => thermal_regions%thermal_region(ireg)
           CALL section_vals_val_get(region_sections,"LIST",i_rep_section=ireg,n_rep_val=nlist,error=error)
           NULLIFY(t_region%part_index)
           t_region%npart = 0
           DO il = 1,nlist

             CALL section_vals_val_get(region_sections,"LIST", i_rep_section=ireg,i_rep_val=il,i_vals=tmplist, error=error)
             CALL reallocate(t_region%part_index,1,t_region%npart+SIZE(tmplist))
             DO i =1,SIZE(tmplist)
               ipart = tmplist(i)
               CPPostcondition(((ipart>0).AND.(ipart<=particles%n_els)),cp_failure_level,routineP,error,failure)
               t_region%npart = t_region%npart + 1
               t_region%part_index(t_region%npart) = ipart
               particles%els(ipart)%t_region_index = ireg
             END DO
           END DO
           CALL section_vals_val_get(region_sections,"TEMPERATURE",&
               i_rep_section=ireg,r_val=temp , error=error)
           t_region%temp_expected = temp
           CALL section_vals_val_get(region_sections,"TEMP_TOL", &
               i_rep_section=ireg,r_val=temp_tol , error=error)
           t_region%temp_tol = temp_tol
         END DO
         simpar%do_thermal_region = .TRUE.
       ELSE
         CALL release_thermal_regions(thermal_regions,error)
         simpar%do_thermal_region = .FALSE.
       END IF
    ELSE
      CALL release_thermal_regions(thermal_regions,error)
      simpar%do_thermal_region = .FALSE.
    END IF

  END SUBROUTINE create_thermal_regions

! *****************************************************************************
!> \brief create thermal_regions
!> \author
! *****************************************************************************
  SUBROUTINE print_thermal_regions(thermal_regions, itimes, time, pos, act , error)
    TYPE(thermal_regions_type), POINTER      :: thermal_regions
    INTEGER, INTENT(IN)                      :: itimes
    REAL(KIND=dp), INTENT(IN)                :: time
    CHARACTER(LEN=default_string_length)     :: pos, act
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_thermal_regions', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: fmd
    INTEGER                                  :: ireg, istat, nregions, unit
    LOGICAL                                  :: failure, new_file
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: temp
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    IF(ASSOCIATED(thermal_regions)) THEN
      print_key => section_vals_get_subs_vals(thermal_regions%section,"PRINT%TEMPERATURE",error=error)
      IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
          unit = cp_print_key_unit_nr(logger,thermal_regions%section,"PRINT%TEMPERATURE",&
                 extension=".tregion",file_position=pos, &
                 file_action=act, is_new_file=new_file, error=error)
          IF(unit > 0) THEN
            IF (new_file) THEN
               WRITE(unit,'(A)')"# Temperature per Region"
               WRITE(unit,'("#",3X,A,2X,A,13X,A)')"Step Nr.","Time[fs]","Temp.[K] ...."
            END IF
            nregions = thermal_regions%nregions
            ALLOCATE (temp(0:nregions), STAT=istat)
            CPPrecondition(istat==0,cp_fatal_level,routineP,error,failure)
            temp = 0.0_dp
            temp(0) =  thermal_regions%temp_reg0
            DO ireg = 1,nregions
              temp(ireg) = thermal_regions%thermal_region(ireg)%temperature
            END DO
            fmd="(I10,F20.3,"//TRIM(ADJUSTL(cp_to_string(nregions+1)))//"F20.6)"
            fmd=TRIM(fmd)
            WRITE(UNIT=unit,FMT=fmd) itimes,time, temp(0:nregions)
            DEALLOCATE(temp, STAT=istat)
            CPPrecondition(istat==0,cp_fatal_level,routineP,error,failure)
          END IF
          CALL cp_print_key_finished_output(unit,logger,thermal_regions%section,"PRINT%TEMPERATURE", error=error)
      END IF
    END IF
  END SUBROUTINE print_thermal_regions

END MODULE thermal_region_utils
