
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!***** cp2k/reftraj_util [1.0] *
!!
!!   NAME
!!     reftraj_util
!!
!!   FUNCTION
!!     Initialize the analysis of trajectories to be done 
!!     by activating the REFTRAJ ensemble
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     Created 10-07 [MI]
!!
!!   SOURCE
!******************************************************************************


MODULE reftraj_util

  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE md,                              ONLY: simulation_parameters_type
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE reftraj_types,                   ONLY: reftraj_type
  USE termination,                     ONLY: stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE


      PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'reftraj_util'

      
      PUBLIC ::   initialize_reftraj
      

!!***
! *****************************************************************************

CONTAINS


!!****f* reftraj_types/initialize_reftraj [1.0] *
!!
!!   NAME
!!    initialize_reftraj 
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     10.2007 created 
!!
!!*** **********************************************************************

  SUBROUTINE initialize_reftraj(reftraj,reftraj_section,md_env,error)

    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(section_vals_type), POINTER         :: reftraj_section
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_reftraj', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=10)                        :: AA
    CHARACTER(LEN=80)                        :: title
    CHARACTER(LEN=default_string_length)     :: name
    INTEGER                                  :: cellunit, iatom, iter, natom, &
                                                Natom_read, nskip, trjunit
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: hx, hy, hz, x, y, z
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

    failure = .FALSE.

    NULLIFY (subsys, force_env, particles, particle_set, simpar)
    CALL get_md_env(md_env=md_env, force_env=force_env, para_env=para_env,&
         simpar=simpar, error=error)
    CALL force_env_get(force_env=force_env,subsys=subsys,error=error)
    CALL cp_subsys_get(subsys=subsys, particles=particles,error=error)
    natom = particles%n_els
    particle_set => particles%els

    trjunit = reftraj%info%trjfile_unit

    IF (para_env%mepos==para_env%source) THEN
    REWIND(trjunit)
    !Check the first snapshot in the trajectory file
    READ(trjunit,*,ERR=999,END=998) natom_read    
    CPPostcondition(natom_read==natom,cp_failure_level,routineP,error,failure)
    READ(trjunit,'(A)',ERR=999,END=998) title
    DO iatom = 1,natom
      READ(trjunit,*,ERR=999,END=998) AA, x, y, z
      name=TRIM(particle_set(iatom)%atomic_kind%element_symbol)
      CPPostcondition((TRIM(AA)==name),cp_failure_level,routineP,error,failure)
    END DO
    REWIND(trjunit)

    
    IF(reftraj%info%first_snapshot>0) THEN
      nskip = reftraj%info%first_snapshot
      DO iter = 1,nskip
        READ(trjunit,*,ERR=999,END=998) natom_read
        READ(trjunit,'(A)',ERR=999,END=998) title
        DO iatom = 1,natom
          READ(trjunit,*,ERR=999,END=998)  AA, x, y, z 
        END DO
      END DO      
    END IF
    IF(reftraj%info%variable_volume) THEN
      cellunit = reftraj%info%cellfile_unit
      IF(reftraj%info%first_snapshot>1) THEN
        DO iter = 1,nskip
          READ(cellunit,'(A)',ERR=997,END=996) title
          READ(cellunit,*,ERR=997,END=996) AA, hx,hy,hz
        END DO
      END IF
    END IF
    END IF

    IF(reftraj%info%last_snapshot>0) THEN
      simpar%nsteps = (reftraj%info%last_snapshot - reftraj%info%first_snapshot)
    END IF

    IF(reftraj%info%eval_ef) THEN
      ! set up and printout
    END IF
    
    IF(reftraj%info%msd) THEN
      ! set up and printout
      CALL initialize_msd_reftraj(reftraj,reftraj_section,error=error)
    END IF
    
    RETURN
996 CONTINUE ! end of file
    CALL stop_program(routineN,moduleN,__LINE__,&
         "End of reference cell file reached")
997 CONTINUE ! error
    CALL stop_program(routineN,moduleN,__LINE__,&
         "Error reading reference cell file")
998 CONTINUE ! end of file
    CALL stop_program(routineN,moduleN,__LINE__,&
         "End of reference trajectory file reached")
999 CONTINUE ! error
    CALL stop_program(routineN,moduleN,__LINE__,&
         "Error reading reference trajectory file")

  END SUBROUTINE initialize_reftraj

!!****f* reftraj_types/initialize_msd_reftraj [1.0] *
!!
!!   NAME
!!    initialize_msd_reftraj 
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     10.2007 created 
!!
!!*** **********************************************************************

  SUBROUTINE initialize_msd_reftraj(reftraj,reftraj_section,error)
    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(section_vals_type), POINTER         :: reftraj_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_msd_reftraj', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.


  
  END SUBROUTINE initialize_msd_reftraj

END MODULE reftraj_util

