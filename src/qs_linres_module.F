!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****m* cp2k/qs_linres_module *
!!
!!   NAME
!!     qs_linres_module
!!
!!   FUNCTION
!!     Contains the setup for  the calculation of properties by linear response
!!     by the application of second order density functional perturbation theory.
!!     The knowledge of the ground state energy, density and wavefunctions is assumed.
!!     Uses the self consistent approach.
!!     Properties that can be calculated : none
!!
!!   AUTHOR
!!      MI 
!!
!!   MODIFICATION HISTORY
!!      created 06-2005 [MI]
!!
!!   SOURCE
!!****

MODULE qs_linres_module

  USE cp_control_types,                ONLY: dft_control_type
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
  USE input_constants,                 ONLY: lr_none, lr_chemshift,&
                                             state_loc_all
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_get_subs_vals2,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_linres_nmr_utils,             ONLY: nmr_do_nmr,&
                                             nmr_env_cleanup,&
                                             nmr_env_init
  USE qs_linres_methods,               ONLY: linres_localize
  USE qs_linres_nmr_op,                ONLY: nmr_operators
  USE qs_linres_types,                 ONLY: linres_control_create,&
                                             linres_control_release,&
                                             linres_control_retain,&
                                             linres_control_type,&
                                             nmr_env_type
  USE qs_loc_control,                  ONLY: read_loc_section,&
                                             set_loc_centers,&
                                             set_loc_wfn_lists
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_p_env_methods,                ONLY: p_env_create,&
                                             p_env_psi0_changed
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: linres_calculation

  CHARACTER(LEN=*), PARAMETER :: moduleN='qs_linres_module'

CONTAINS

! *****************************************************************************
!!****f* cp2k/qs_linres_module/linres_calculation
!!
!!   NAME
!!     linres_calculation
!!
!!   FUNCTION
!!     Driver for the linear response calculatios
!!     Initialize the perturbation environment
!!     Define which properties is to be calculated
!!     Start up the optimization of the response density and wfn
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     06.2005 created [MI]
!!
!****************************************************************************

  SUBROUTINE linres_calculation(qs_env, globenv, error)

    TYPE(qs_environment_type), POINTER          :: qs_env
    TYPE(global_environment_type),  &
      POINTER                                   :: globenv
    TYPE(cp_error_type), &
         INTENT(INOUT), OPTIONAL                :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_calculation', &
                                   routineP = moduleN//'/'//routineN
 
    INTEGER                                     :: handle,output_unit
    LOGICAL                                     :: ionode, failure, &
                                                   lr_calculation, was_present
    TYPE(cp_logger_type), POINTER               :: logger    
    TYPE(dft_control_type), POINTER             :: dft_control
    TYPE(linres_control_type), POINTER          :: linres_control
    TYPE(nmr_env_type)                          :: nmr_env
    TYPE(qs_p_env_type), POINTER                :: p_env
    TYPE(section_vals_type), POINTER            :: lr_section, prop_section

    CALL timeset(routineN,"I"," ",handle)    
    failure = .FALSE. 
    lr_calculation = .FALSE. 
    was_present = .FALSE.

    NULLIFY(dft_control,p_env,linres_control,logger,lr_section,prop_section)

    logger => cp_error_get_logger(error)
    ionode = (logger%para_env%mepos==logger%para_env%source)
    IF(ionode) output_unit = cp_logger_get_default_unit_nr(logger)

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    ! The transition potential method to calculate XAS needs LSD
    lr_section =>  section_vals_get_subs_vals(globenv%input_file, &
                                                 "FORCE_EVAL%PROPERTIES%LINRES")
    call section_vals_get(lr_section,explicit=lr_calculation,error=error)

    IF (lr_calculation) THEN


       CALL linres_init(lr_section,p_env,qs_env,ionode,globenv,error=error)

       CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, &
            linres_control=linres_control)

       !The type of perturbation has not been defined yet
       linres_control%property = lr_none

       prop_section => section_vals_get_subs_vals(lr_section,"NMR",error=error)
       CALL section_vals_get(prop_section,explicit=was_present,error=error)
       IF(was_present)THEN
         linres_control%property = lr_chemshift
         linres_control%do_kernel = .FALSE.

         IF(.NOT.linres_control%localized_psi0) THEN
            CALL stop_program(routineP,&
             "Are you sure that you want to calculate the chemical shift without localized psi0??")

           CALL linres_localize(linres_control%localized_wfn_control,qs_env,&
                 dft_control%nspins,globenv,centers_only=.TRUE.,error=error)
         END IF

         !Initialize the nmr environment
         CALL nmr_env_init(nmr_env,qs_env,globenv,error=error)
         CALL nmr_operators(nmr_env,qs_env,error=error)
         CALL nmr_do_nmr(nmr_env,p_env,qs_env,globenv,error=error)
         CALL nmr_env_cleanup(nmr_env,error=error)

       END IF

      ! Other possible LR calculations can be introduced here
      
       CALL p_env_release(p_env)


    ELSE
       IF(ionode)  WRITE (output_unit, "(2X,A)") "",&
       "-----------------------------------------------------------------------------",&
       "-              No LR calculation has been specified in the input            -",&
       "                        cp2k is going to stop, bye bye                       ",&
       "-----------------------------------------------------------------------------",&
       ""    
    END IF 

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE linres_calculation

! *****************************************************************************
!!****f* cp2k/qs_linres_module/linres_init
!!
!!   NAME
!!     linres_init
!!
!!   FUNCTION
!!     Initialize some general settings like the p_env
!!     Localize the psi0 if required
!!
!!   NOTES
!!     - The localization should probably be always for all the occupied states
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     06.2005 created [MI]
!!
!****************************************************************************

  SUBROUTINE linres_init(lr_section,p_env,qs_env,ionode,globenv,error)

    TYPE(section_vals_type), POINTER            :: lr_section
    TYPE(qs_environment_type), POINTER          :: qs_env
    TYPE(qs_p_env_type), POINTER                :: p_env
    LOGICAL, INTENT(IN)                         :: ionode
    TYPE(global_environment_type),  &
      POINTER                                   :: globenv
    TYPE(cp_error_type), & 
         INTENT(INOUT), OPTIONAL                :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_init', &
                                   routineP = moduleN//'/'//routineN
 
    INTEGER :: istat, ispin, n_mo(2),nelectron,nmoloc(2)
    LOGICAL :: failure
    REAL(dp) :: maxocc
    TYPE(dft_control_type), POINTER             :: dft_control
    TYPE(linres_control_type), POINTER          :: linres_control
    TYPE(section_vals_type), POINTER            :: loc_section

    failure = .FALSE.

    NULLIFY(dft_control,linres_control,loc_section)

    CPPrecondition(.NOT.ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN


       CALL linres_control_create(linres_control,error=error)
       CALL set_qs_env(qs_env=qs_env, linres_control=linres_control)
       CALL linres_control_release(linres_control,error=error)
       CALL get_qs_env(qs_env=qs_env, linres_control=linres_control,&
                       dft_control=dft_control)

       ! Localized Psi0 are required when the position operator has to be defined (nmr)
       CALL section_vals_val_get(lr_section,"LOCALIZE_PSI0",&
           l_val=linres_control%localized_psi0,error=error)
       IF(linres_control%localized_psi0) THEN
         loc_section =>section_vals_get_subs_vals(lr_section,"LOCALIZE",error=error)
         CALL read_loc_section(linres_control%localized_wfn_control,loc_section,error=error)
         CALL linres_localize(linres_control%localized_wfn_control,qs_env,&
               dft_control%nspins,globenv,error=error)
       END IF

       CALL section_vals_val_get(lr_section,"OPT_METHOD",i_val=linres_control%opt_method,error=error)
       CALL section_vals_val_get(lr_section,"MAX_SCF",i_val=linres_control%max_scf,error=error)
       CALL section_vals_val_get(lr_section,"EPS_SCF",r_val=linres_control%eps_scf,error=error)
       CALL section_vals_val_get(lr_section,"LINESEARCH",i_val=linres_control%ls_method)
       CALL section_vals_val_get(lr_section,"STEPSIZE",r_val=linres_control%ds_min)

       CALL section_vals_val_get(lr_section,"PRECONDITIONER",i_val=linres_control%preconditioner_type)
       CALL section_vals_val_get(lr_section,"ENERGY_GAP",r_val=linres_control%energy_gap)
       
      !------------------!
      ! create the p_env !
      !------------------!
       CALL p_env_create(p_env, qs_env, orthogonal_orbitals=.TRUE.,linres_control=linres_control,error=error)

       ! update the m_epsilon matrix
       CALL p_env_psi0_changed(p_env,qs_env,error=error)
    END IF

  END SUBROUTINE linres_init


! *****************************************************************************
END MODULE qs_linres_module


