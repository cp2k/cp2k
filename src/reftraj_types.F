!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!***** cp2k/reftraj_types [1.0] *
!!
!!   NAME
!!     reftraj_types
!!
!!   FUNCTION
!!    initialization of the reftraj structure used to analyse 
!!    previously generated trajectories 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     Created 10-07 [MI]
!!
!!   SOURCE
!******************************************************************************


MODULE reftraj_types

  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE


    PRIVATE
    PUBLIC :: reftraj_type, reftraj_info_type, reftraj_msd_type,&
              create_reftraj,&
              release_reftraj,&
              retain_reftraj


!!****f* reftraj_types/reftraj_type [1.0] *
!!
!!   NAME
!!     reftraj_type, reftraj_msd_type, reftraj_info_type
!!
!!   FUNCTION
!!      parameters related to the analysis of previously generated trajecorties
!! 
!!   AUTHOR
!!     MI
!!
!!*** **********************************************************************
  TYPE reftraj_info_type
     INTEGER :: trjfile_unit
     INTEGER :: cellfile_unit
     INTEGER :: first_snapshot
     INTEGER :: last_snapshot
     LOGICAL :: eval_ef
     LOGICAL :: variable_volume
     LOGICAL :: msd
  END TYPE reftraj_info_type

  TYPE reftraj_msd_type
     INTEGER :: ref0_unit
     LOGICAL :: msd_kind, msd_molecule, msd_region
     REAL(KIND=dp) :: ref0_com(3)
     REAL(KIND=dp), POINTER, DIMENSION(:,:) :: ref0_pos
     REAL(KIND=dp), POINTER, DIMENSION(:,:) :: ref0_com_molecule 
     REAL(KIND=dp), POINTER, DIMENSION(:,:) :: val_msd_kind
     REAL(KIND=dp), POINTER, DIMENSION(:,:) :: val_msd_molecule
     REAL(KIND=dp), POINTER, DIMENSION(:,:) :: val_msd_region
  END TYPE reftraj_msd_type

 
  TYPE reftraj_type
     INTEGER :: ref_count
     INTEGER :: itimes
     INTEGER :: itimes0
     INTEGER :: isnap
     INTEGER :: natom
     LOGICAL :: init
     REAL(KIND=dp) :: epot, epot0, time, time0
     TYPE(reftraj_info_type), POINTER   :: info
     TYPE(reftraj_msd_type), POINTER    :: msd
  END TYPE reftraj_type
  
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'reftraj_types'

!!***
! *****************************************************************************

CONTAINS


  SUBROUTINE create_reftraj(reftraj,reftraj_section,error)
  
    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(section_vals_type), POINTER         :: reftraj_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_reftraj', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(.NOT. ASSOCIATED(reftraj),cp_failure_level,routineP,error,failure)
    ALLOCATE(reftraj, stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    reftraj%ref_count = 1

    NULLIFY(reftraj%info)
    NULLIFY(reftraj%msd)
    

    ALLOCATE(reftraj%info, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    
    ! Some Initialization
    reftraj%info%trjfile_unit = -1
    CALL section_vals_val_get(reftraj_section,"TRJFILENAME",c_val=filename,error=error)
    CALL open_file(TRIM(filename),unit_number=reftraj%info%trjfile_unit)
    !MI: It would be good to check the topology to be sure that it is the same of
    ! the input

    
    CALL section_vals_val_get(reftraj_section,"VARIABLE_VOLUME",l_val=reftraj%info%variable_volume,error=error)

    reftraj%info%cellfile_unit = -1
    IF(reftraj%info%variable_volume) THEN
       CALL section_vals_val_get(reftraj_section,"CELLFILENAME",c_val=filename,error=error)
       CALL open_file(TRIM(filename),unit_number=reftraj%info%cellfile_unit)
    END IF

    CALL section_vals_val_get(reftraj_section,"FIRST_SNAPSHOT",i_val=reftraj%info%first_snapshot,error=error)
    CALL section_vals_val_get(reftraj_section,"LAST_SNAPSHOT",i_val=reftraj%info%last_snapshot,error=error)

    CALL section_vals_val_get(reftraj_section,"EVAL_ENERGY_FORCES",l_val=reftraj%info%eval_ef,error=error)
  
    CALL section_vals_val_get(reftraj_section,"MSD%_SECTION_PARAMETERS_",&
            l_val=reftraj%info%msd,error=error)

  END SUBROUTINE create_reftraj

!!****f* reftraj_types/retain_reftraj [1.0] *
!!
!!   NAME
!!    retain_reftraj 
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
  SUBROUTINE retain_reftraj(reftraj, error)

    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'retain_reftraj', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    
    IF (ASSOCIATED(reftraj)) THEN
       CPPrecondition(reftraj%ref_count>0,cp_failure_level,routineP,error,failure)
       reftraj%ref_count=reftraj%ref_count+1
    END IF

  END SUBROUTINE retain_reftraj


!!****f* reftraj_types/release_reftraj [1.0] *
!!
!!   NAME
!!    release_reftraj 
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
  SUBROUTINE release_reftraj(reftraj, error)

    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_reftraj', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF(ASSOCIATED(reftraj)) THEN
       CPPrecondition(reftraj%ref_count>0,cp_failure_level,routineP,error,failure)
       reftraj%ref_count=reftraj%ref_count-1
       IF(reftraj%ref_count<1) THEN
          IF ( reftraj%info%trjfile_unit /= -1 ) THEN 
              CALL close_file(unit_number=reftraj%info%trjfile_unit)
          END IF
          IF ( reftraj%info%cellfile_unit /= -1 ) THEN 
              CALL close_file(unit_number=reftraj%info%cellfile_unit)
          END IF
          IF(ASSOCIATED(reftraj%info)) THEN
             DEALLOCATE (reftraj%info, STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF(ASSOCIATED(reftraj%msd)) THEN
             DEALLOCATE(reftraj%msd%ref0_pos,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
             IF(reftraj%msd%msd_kind) THEN
                 DEALLOCATE(reftraj%msd%val_msd_kind,STAT=istat)
                 CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
             END IF
             IF(reftraj%msd%msd_molecule) THEN
                 DEALLOCATE(reftraj%msd%val_msd_molecule,STAT=istat)
                 CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                 DEALLOCATE(reftraj%msd%ref0_com_molecule,STAT=istat)
                 CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
             END IF

             DEALLOCATE (reftraj%msd, STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
          
          DEALLOCATE( reftraj , stat=istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF

    END IF
  END SUBROUTINE release_reftraj

END MODULE reftraj_types
