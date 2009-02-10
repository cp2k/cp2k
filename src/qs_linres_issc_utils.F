!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Chemical shift calculation by dfpt
!>      Initialization of the issc_env, creation of the special neighbor lists
!>      Perturbation Hamiltonians by application of the p and rxp oprtators to  psi0
!>      Write output
!>      Deallocate everything
!> \note
!>      The psi0 should be localized
!>      the Sebastiani method works within the assumption that the orbitals are
!>      completely contained in the simulation box
! *****************************************************************************
MODULE qs_linres_issc_utils

  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_types,                 ONLY: issc_env_create,&
                                             issc_env_type,&
                                             linres_control_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: issc_env_cleanup, issc_env_init

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_issc_utils'

CONTAINS

! *****************************************************************************
!> \brief Initialize the issc environment
! *****************************************************************************
  SUBROUTINE issc_env_init(issc_env,qs_env,error)
    !
    TYPE(issc_env_type)                      :: issc_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'issc_env_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, output_unit
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(section_vals_type), POINTER         :: issc_section, lr_section

!

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(linres_control)
    NULLIFY(logger,issc_section)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         &                             extension=".linresLog",error=error)
    
    IF(issc_env%ref_count /= 0) THEN
       CALL issc_env_cleanup(issc_env,error=error)
    ENDIF

    IF(output_unit>0) THEN
       WRITE(output_unit,"(/,T20,A,/)") "*** Start indirect spin-spin coupling Calculation ***"
       WRITE(output_unit,"(T10,A,/)") "Inizialization of the ISSC environment"
    ENDIF

    CALL issc_env_create(issc_env,error=error)
    !
    issc_section => section_vals_get_subs_vals(qs_env%input, &
         &          "PROPERTIES%LINRES%SPINSPIN",error=error)
    !CALL section_vals_val_get(nmr_section,"INTERPOLATE_SHIFT",l_val=nmr_env%interpolate_shift,error=error)
    !CALL section_vals_val_get(nmr_section,"SHIFT_GAPW_RADIUS",r_val=nmr_env%shift_gapw_radius,error=error)

    CALL get_qs_env(qs_env=qs_env,&
         &          linres_control=linres_control,&
         &          error=error)
    !
    ! check that the psi0 are localized and you have all the centers
    CPPrecondition(linres_control%localized_psi0,cp_warning_level,routineP,error,failure)
    IF(failure.AND.(output_unit>0)) THEN
       WRITE(output_unit,'(A)') &
            ' To get indirect spin-spin coupling parameters within PBC you need to localize zero order orbitals '
    ENDIF
    !
    ! Conversion factors
    IF(output_unit>0) THEN
       !WRITE(output_unit,"(T2,A,T60,ES15.6)")&
       !     & "  NMR|  shift gapw radius (a.u.) ",nmr_env%shift_gapw_radius
    ENDIF

    !
    ! Initialize the issc tensor
    !ALLOCATE(issc_env%issc(3,3,natom),&
    !     &   issc_env%issc_loc(3,3,natom),STAT=istat)
    !CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    !issc_env%issc = 0.0_dp
    !issc_env%issc_loc = 0.0_dp

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         &                            "PRINT%PROGRAM_RUN_INFO",error=error)

    CALL timestop(handle)

  END SUBROUTINE issc_env_init

! *****************************************************************************
!> \brief Deallocate the issc environment
!> \par History
! *****************************************************************************
  SUBROUTINE issc_env_cleanup(issc_env,error)

    TYPE(issc_env_type)                      :: issc_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'issc_env_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    issc_env%ref_count = issc_env%ref_count - 1
    IF(issc_env%ref_count == 0 ) THEN
       IF(ASSOCIATED(issc_env%issc)) THEN
          DEALLOCATE(issc_env%issc,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       IF(ASSOCIATED(issc_env%issc_loc)) THEN
          DEALLOCATE(issc_env%issc_loc,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
    ENDIF ! ref count

  END SUBROUTINE issc_env_cleanup

END MODULE qs_linres_issc_utils
