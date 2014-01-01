!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief initialization of the reftraj structure used to analyse
!>     previously generated trajectories
!> \par History
!>      Created 10-07 [MI]
!> \author MI
! *****************************************************************************
MODULE reftraj_types

  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp
#include "cp_common_uses.h"

  IMPLICIT NONE

    PRIVATE
    PUBLIC :: reftraj_type, reftraj_info_type, reftraj_msd_type,&
              create_reftraj, release_reftraj, retain_reftraj

! *****************************************************************************
!> \brief parameters related to the analysis of previously generated trajecorties
!> \author MI
! *****************************************************************************
  TYPE reftraj_info_type
     INTEGER                                  :: first_snapshot
     INTEGER                                  :: last_snapshot
     INTEGER                                  :: stride
     LOGICAL                                  :: eval_ef
     LOGICAL                                  :: variable_volume
     LOGICAL                                  :: msd
     TYPE(cp_parser_type), POINTER            :: traj_parser
     TYPE(cp_parser_type), POINTER            :: cell_parser
  END TYPE reftraj_info_type

! *****************************************************************************
  TYPE reftraj_msd_type
     LOGICAL                                  :: disp_atom, msd_kind, msd_molecule, msd_region
     INTEGER                                  :: num_disp_atom, ref0_unit
     INTEGER, POINTER, DIMENSION(:)           :: disp_atom_index
     REAL(KIND=dp)                            :: disp_atom_tol, drcom(3), ref0_com(3), total_mass
     REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: disp_atom_dr
     REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: ref0_pos
     REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: ref0_com_molecule
     REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: val_msd_kind
     REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: val_msd_molecule
     REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: val_msd_region
  END TYPE reftraj_msd_type

! *****************************************************************************
  TYPE reftraj_type
     INTEGER                                  :: ref_count
     INTEGER                                  :: itimes
     INTEGER                                  :: itimes0
     INTEGER                                  :: isnap
     INTEGER                                  :: natom
     LOGICAL                                  :: init
     REAL(KIND=dp)                            :: epot, epot0, time, time0
     TYPE(reftraj_info_type), POINTER         :: info
     TYPE(reftraj_msd_type), POINTER          :: msd
  END TYPE reftraj_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'reftraj_types'

CONTAINS

! *****************************************************************************
  SUBROUTINE create_reftraj(reftraj,reftraj_section,para_env,error)

    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(section_vals_type), POINTER         :: reftraj_section
    TYPE(cp_para_env_type), POINTER          :: para_env
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
    NULLIFY(reftraj%info%traj_parser)
    NULLIFY(reftraj%info%cell_parser)

    ! Initialize parser for trajectory
    CALL section_vals_val_get(reftraj_section,"TRAJ_FILE_NAME",c_val=filename,error=error)
    CALL parser_create(reftraj%info%traj_parser,filename,para_env=para_env,error=error)

    CALL section_vals_val_get(reftraj_section,"VARIABLE_VOLUME",l_val=reftraj%info%variable_volume,&
         error=error)
    IF(reftraj%info%variable_volume) THEN
       ! In case requested initialize parser for cell
       CALL section_vals_val_get(reftraj_section,"CELL_FILE_NAME",c_val=filename,error=error)
       CALL parser_create(reftraj%info%cell_parser,filename,para_env=para_env,error=error)
    END IF

    CALL section_vals_val_get(reftraj_section,"FIRST_SNAPSHOT",i_val=reftraj%info%first_snapshot,error=error)
    CALL section_vals_val_get(reftraj_section,"LAST_SNAPSHOT",i_val=reftraj%info%last_snapshot,error=error)
    CALL section_vals_val_get(reftraj_section,"STRIDE",i_val=reftraj%info%stride,error=error)
    CALL section_vals_val_get(reftraj_section,"EVAL_ENERGY_FORCES",l_val=reftraj%info%eval_ef,error=error)

    CALL section_vals_val_get(reftraj_section,"MSD%_SECTION_PARAMETERS_",&
            l_val=reftraj%info%msd,error=error)

  END SUBROUTINE create_reftraj

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created
!> \author MI
! *****************************************************************************
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

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created
!> \author MI
! *****************************************************************************
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
          CALL parser_release(reftraj%info%traj_parser,error=error)
          CALL parser_release(reftraj%info%cell_parser,error=error)
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
             IF(reftraj%msd%disp_atom) THEN
                 DEALLOCATE(reftraj%msd%disp_atom_index,STAT=istat)
                 CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                 DEALLOCATE(reftraj%msd%disp_atom_dr,STAT=istat)
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
