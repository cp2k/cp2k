
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****m* cp2k/qs_linres_methods *
!!
!!   NAME 
!!     qs_linres_methods
!!
!!   FUNCTION
!!     localize wavefunctions
!!     linear response scf
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     created 07-2005 [MI]
!!
!!   SOURCE
!!****

MODULE qs_linres_methods

  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type 
  USE kinds,                           ONLY: dp
  USE input_constants,                 ONLY: do_loc_none,state_loc_all,&
                                             op_loc_berry
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_loc_control,                  ONLY: localized_wfn_control_type,&
                                             set_loc_centers,&
                                             set_loc_wfn_lists
  USE qs_loc_methods,                  ONLY: qs_loc_driver
  USE qs_loc_types,                    ONLY: qs_loc_env_release,&
                                             qs_loc_env_retain,&
                                             qs_loc_env_create,&
                                             qs_loc_env_destroy,&
                                             qs_loc_env_new_type
  USE qs_loc_utils,                    ONLY: qs_loc_env_init
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE


  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: linres_localize

  CHARACTER(LEN=*), PARAMETER :: moduleN='qs_linres_methods'

!!***
! *****************************************************************************

CONTAINS
                                      
! *****************************************************************************
!!****f* cp2k/qs_linres_methods/linres_localize
!!
!!   NAME
!!     linres_localize
!!
!!   FUNCTION
!!     Find the centers and spreads of the wfn, 
!!     if required apply a localization algorithm
!!
!!   NOTES
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     07.2005 created [MI]
!!
!****************************************************************************

  SUBROUTINE linres_localize(localized_wfn_control,qs_env,nspins,globenv,centers_only,error)

    TYPE(localized_wfn_control_type), POINTER   :: localized_wfn_control
    TYPE(qs_environment_type), POINTER          :: qs_env
    INTEGER, INTENT(IN)                         :: nspins
    LOGICAL, INTENT(IN), OPTIONAL               :: centers_only
    TYPE(global_environment_type),  &
      POINTER                                   :: globenv
    TYPE(cp_error_type), & 
         INTENT(INOUT), OPTIONAL                :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'boh', &
                                   routineP = moduleN//'/'//routineN

    INTEGER :: ispin, istat, n_mo(2), nelectron, nmoloc(2)
    LOGICAL :: failure, ionode, my_centers_only
    REAL(dp) :: maxocc
    TYPE(cp_logger_type), POINTER               :: logger
    TYPE(mo_set_p_type), DIMENSION(:), POINTER  :: mos
    TYPE(qs_loc_env_new_type), POINTER          :: qs_loc_env     

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%source==logger%para_env%mepos
    my_centers_only = .FALSE.
    IF(PRESENT(centers_only)) my_centers_only = centers_only
    nmoloc(1:2) = 0
    NULLIFY(mos,qs_loc_env)

    IF(my_centers_only) THEN
       localized_wfn_control%set_of_states = state_loc_all
       localized_wfn_control%localization_method = do_loc_none
       localized_wfn_control%operator_type = op_loc_berry
    END IF

    CALL get_qs_env(qs_env=qs_env, mos=mos)
    DO ispin = 1,nspins
       CALL get_mo_set(mos(ispin)%mo_set,nmo=n_mo(ispin),nelectron=nelectron,maxocc=maxocc)
       IF(localized_wfn_control%set_of_states == state_loc_all) THEN
         nmoloc(ispin) = NINT(nelectron/maxocc)
       ELSE
         nmoloc(ispin) = MIN(localized_wfn_control%nloc_states(1) ,n_mo(ispin))
       END IF
    END DO  ! ispin
    CALL set_loc_wfn_lists(localized_wfn_control,&
         nmoloc,n_mo,nspins,error=error)
    CALL set_loc_centers(localized_wfn_control,nmoloc,&
         nspins,globenv,error=error) 
    CALL qs_loc_env_create(qs_loc_env,error=error)
    CALL qs_loc_env_init(qs_loc_env,localized_wfn_control,qs_env,error=error)
    ! The localized orbitals are copied in the mos, 
    ! The orbital centers are stored in linres_control%localized_wfn_control
    CALL qs_loc_driver(qs_env,qs_loc_env,ionode,error=error)
    CALL qs_loc_env_destroy(qs_loc_env, error=error)
    DEALLOCATE(qs_loc_env,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE linres_localize
  
! *****************************************************************************

END MODULE qs_linres_methods
