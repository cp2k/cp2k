!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Barostat structure: module containing barostat available for MD
!> \author teo [tlaino] - University of Zurich - 09.2007
! *****************************************************************************
MODULE barostat_types
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE extended_system_init,            ONLY: initialize_npt
  USE extended_system_types,           ONLY: npt_info_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: npe_f_ensemble,&
                                             nph_uniaxial_damped_ensemble,&
                                             nph_uniaxial_ensemble,&
                                             npt_f_ensemble,&
                                             npt_i_ensemble
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: barostat_type,&
            create_barostat_type,&
            release_barostat_type,&
            retain_barostat_type

! *****************************************************************************
  TYPE barostat_type
     INTEGER :: id_nr, ref_count, type_of_barostat
     REAL(KIND=dp) :: temp_ext
     TYPE ( npt_info_type ), POINTER :: npt (:,:)
  END TYPE barostat_type

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'barostat_types'
  INTEGER, PRIVATE, SAVE               :: last_barostat_id_nr=0

CONTAINS

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2007 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE create_barostat_type( barostat, motion_section, force_env,&
       simpar, para_env, globenv, error )
    TYPE(barostat_type), POINTER             :: barostat
    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_barostat_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: check, explicit, failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: barostat_section

    failure =.FALSE.
    check = .NOT.ASSOCIATED(barostat)
    CPPrecondition(check,cp_failure_level,routineP,error,failure)
    barostat_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT",error=error)
    CALL section_vals_get(barostat_section, explicit=explicit, error=error)
    IF (simpar%ensemble == npt_i_ensemble .OR. &
        simpar%ensemble == npt_f_ensemble .OR. &
        simpar%ensemble == npe_f_ensemble .OR. &
        simpar%ensemble == nph_uniaxial_ensemble .OR. &
        simpar%ensemble == nph_uniaxial_damped_ensemble) THEN
       ALLOCATE(barostat, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       last_barostat_id_nr = last_barostat_id_nr + 1
       barostat%id_nr      = last_barostat_id_nr
       barostat%ref_count  =  1
       NULLIFY(barostat%npt)
       CALL force_env_get( force_env, cell=cell, error=error )

       CALL section_vals_val_get(barostat_section,"TEMPERATURE",r_val=barostat%temp_ext,error=error)
       IF(barostat%temp_ext < 0.0_dp) barostat%temp_ext = simpar%temp_ext
       simpar%temp_baro_ext =  barostat%temp_ext 
       ! Initialize or possibly restart Barostat
       CALL initialize_npt (simpar, para_env, globenv, barostat%npt,&
            cell, work_section=barostat_section, error=error)

       ! If none of the possible barostat has been allocated let's deallocate
       ! the full structure
       IF(.NOT.ASSOCIATED(barostat%npt)) THEN
          CALL  release_barostat_type(barostat, error)
       END IF
    END IF
    
  END SUBROUTINE create_barostat_type

! *****************************************************************************
!> \brief retains the given barostat
!> \par History
!>      09.2007 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE retain_barostat_type(barostat,error)
    TYPE(barostat_type), POINTER             :: barostat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'retain_barostat_type', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    
    IF (.NOT. failure) THEN
       IF (ASSOCIATED(barostat)) THEN
          CPPrecondition(barostat%ref_count>0,cp_failure_level,routineP,error,failure)
          barostat%ref_count=barostat%ref_count+1
       END IF
    END IF
  END SUBROUTINE retain_barostat_type

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2007 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE release_barostat_type(barostat, error)
    TYPE(barostat_type), POINTER             :: barostat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_barostat_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: check, failure

    failure =.FALSE.
    IF (ASSOCIATED(barostat)) THEN
       check = barostat%ref_count>0
       CPPrecondition(check,cp_failure_level,routineP,error,failure)
       barostat%ref_count=barostat%ref_count-1
       IF (barostat%ref_count<1) THEN
          IF ( ASSOCIATED ( barostat%npt ) ) THEN
             DEALLOCATE (barostat%npt , stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          DEALLOCATE(barostat, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF

  END SUBROUTINE release_barostat_type

END MODULE barostat_types
