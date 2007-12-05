!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****** cp2k/force_fields [1.0] *
!!
!!   NAME
!!     force_fields
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     Subroutine input_torsions changed (DG) 05-Dec-2000
!!     Output formats changed (DG) 05-Dec-2000
!!     JGH (26-01-2002) : force field parameters stored in tables, not in
!!       matrices. Input changed to have parameters labeled by the position
!!       and not atom pairs (triples etc)
!!     Teo (11.2005) : Moved all information on force field  pair_potential to
!!                     a much lighter memory structure
!!
!!   SOURCE
!******************************************************************************

MODULE force_fields

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE ewald_environment_types,         ONLY: ewald_environment_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type
  USE force_field_types,               ONLY: deallocate_ff_type,&
                                             force_field_type,&
                                             init_ff_type
  USE force_fields_input,              ONLY: read_force_field_charmm,&
                                             read_force_field_gromos,&
                                             read_force_field_section
  USE force_fields_util,               ONLY: clean_intra_force_kind,&
                                             force_field_pack,&
                                             force_field_qeff_output
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_fields'

  PRIVATE
  PUBLIC :: force_field_control

!!*****
!******************************************************************************

CONTAINS
  
!******************************************************************************
!!****** force_fields/force_field_control [1.0] *
!!
!!   NAME
!!     force_field_control
!!
!!   FUNCTION
!!     1. If reading in from external file, make sure its there first
!!     2. Read in the force_field from the corresponding locations
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE force_field_control (atomic_kind_set,particle_set,&
       molecule_kind_set,molecule_set,&
       ewald_env,fist_nonbond_env,&
       root_section,para_env,globenv, qmmm, qmmm_env, subsys_section, mm_section,&
       shell_particle_set, core_particle_set, cell, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(IN), OPTIONAL            :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(section_vals_type), POINTER         :: subsys_section, mm_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: shell_particle_set, &
                                                core_particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw
    LOGICAL                                  :: failure, found
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(force_field_type)                   :: ff_type

    CALL timeset(routineN,'I','',handle)
    failure = .FALSE.
    logger => cp_error_get_logger(error)


    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO",&
         extension=".mmLog",error=error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Initialize the ff_type structure type
    !-----------------------------------------------------------------------------
    CALL init_ff_type(ff_type,error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Read in the force field section in the input file if any
    !-----------------------------------------------------------------------------
    CALL read_force_field_section(ff_type,para_env,globenv,mm_section,error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. If reading in from external file, make sure its there first
    !-----------------------------------------------------------------------------
    IF(ff_type%ff_type/="OFF") THEN
       IF(para_env%ionode) THEN
          INQUIRE(FILE=ff_type%ff_file_name,EXIST=found)
          IF(found) THEN
          ELSE
             CALL stop_program ("force_field_control","force field file missing")
          END IF
       END IF
    END IF
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Read in the force field from the corresponding locations
    !-----------------------------------------------------------------------------
    IF(ff_type%ff_type/="OFF") THEN
       SELECT CASE (ff_type%ff_type)
       CASE ("CHM")
          CALL read_force_field_charmm(ff_type,para_env,mm_section,error)
       CASE ("G87","G96")
          CALL read_force_field_gromos(ff_type,para_env,mm_section,error)
       CASE DEFAULT
          CALL stop_program ("force_field_control","force field type not implemented")
       END SELECT
    END IF
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Pack all force field info into different structures
    !-----------------------------------------------------------------------------
    CALL force_field_pack (particle_set,atomic_kind_set,&
         molecule_kind_set,molecule_set,&
         ewald_env,fist_nonbond_env,&
         ff_type,root_section, qmmm, qmmm_env, mm_section, subsys_section, &
         shell_particle_set=shell_particle_set, core_particle_set=core_particle_set,&
         cell=cell, error=error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 6. Output total system charge assigned to qeff
    !-----------------------------------------------------------------------------
    CALL force_field_qeff_output (particle_set,atomic_kind_set,&
         molecule_kind_set,molecule_set,mm_section, ff_type, error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 7. Clean up "UNSET" bond,bend,UB,TORSION,IMPR,ONFO kinds
    !-----------------------------------------------------------------------------
    CALL clean_intra_force_kind (molecule_kind_set,mm_section,error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 8. Cleanup the ff_type structure type
    !-----------------------------------------------------------------------------
    CALL deallocate_ff_type(ff_type, error=error)

    CALL cp_print_key_finished_output(iw,logger,mm_section,&
         "PRINT%FF_INFO",error=error)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE force_field_control

END MODULE force_fields
!******************************************************************************
