!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Contains the setup for  the calculation of properties by linear response
!>      by the application of second order density functional perturbation theory.
!>      The knowledge of the ground state energy, density and wavefunctions is assumed.
!>      Uses the self consistent approach.
!>      Properties that can be calculated : none
!> \par History
!>       created 06-2005 [MI]
!> \author MI
! *****************************************************************************
MODULE qs_linres_module
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             use_qmmm,&
                                             use_qs_force
  USE input_constants,                 ONLY: lr_current,&
                                             lr_none,&
                                             ot_precond_full_all,&
                                             ot_precond_full_kinetic,&
                                             ot_precond_full_single,&
                                             ot_precond_full_single_inverse,&
                                             ot_precond_none,&
                                             ot_precond_s_inverse
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE qmmm_types,                      ONLY: primary_subsys,&
                                             qs_subsys
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_linres_current,               ONLY: current_build_chi,&
                                             current_build_current
  USE qs_linres_current_utils,         ONLY: current_env_cleanup,&
                                             current_env_init,&
                                             current_response
  USE qs_linres_epr_nablavks,          ONLY: epr_nablavks
  USE qs_linres_epr_ownutils,          ONLY: epr_g_print,&
                                             epr_g_so,&
                                             epr_g_soo,&
                                             epr_g_zke,&
                                             epr_ind_magnetic_field
  USE qs_linres_epr_utils,             ONLY: epr_env_cleanup,&
                                             epr_env_init
  USE qs_linres_issc_utils,            ONLY: issc_env_cleanup,&
                                             issc_env_init,&
                                             issc_issc,&
                                             issc_print,&
                                             issc_response
  USE qs_linres_methods,               ONLY: linres_localize
  USE qs_linres_nmr_shift,             ONLY: nmr_shift,&
                                             nmr_shift_print
  USE qs_linres_nmr_utils,             ONLY: nmr_env_cleanup,&
                                             nmr_env_init
  USE qs_linres_op,                    ONLY: current_operators,&
                                             issc_operators,&
                                             polar_operators
  USE qs_linres_polar_utils,           ONLY: polar_env_cleanup,&
                                             polar_env_init,&
                                             polar_polar,&
                                             polar_print,&
                                             polar_response
  USE qs_linres_types,                 ONLY: current_env_type,&
                                             epr_env_type,&
                                             issc_env_type,&
                                             linres_control_create,&
                                             linres_control_release,&
                                             linres_control_type,&
                                             nmr_env_type,&
                                             polar_env_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: mo_set_p_type
  USE qs_p_env_methods,                ONLY: p_env_create,&
                                             p_env_psi0_changed
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: linres_calculation, linres_calculation_low

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_module'

CONTAINS

! *****************************************************************************
!> \brief Driver for the linear response calculatios
!> \param force_env ...
!> \param error ...
!> \par History
!>      06.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE linres_calculation(force_env, error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_calculation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(qs_environment_type), POINTER       :: qs_env

    CALL timeset(routineN,handle)

    NULLIFY(qs_env)

    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)

    SELECT CASE(force_env%in_use)
    CASE(use_qs_force)
       CALL force_env_get(force_env,&
                          qs_env=qs_env,error=error)
    CASE(use_qmmm)
       CALL force_env_get(force_env%sub_force_env(primary_subsys)%force_env%sub_force_env(qs_subsys)%force_env,&
                          qs_env=qs_env,error=error)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Doesnt recognize this force_env.")
    END SELECT

    qs_env%linres_run = .TRUE.

    CALL linres_calculation_low(qs_env,  error=error)

    CALL timestop(handle)

  END SUBROUTINE linres_calculation

! *****************************************************************************
!> \brief Linear response can be called as run type or as post scf calculation
!>      Initialize the perturbation environment
!>      Define which properties is to be calculated
!>      Start up the optimization of the response density and wfn
!> \param qs_env ...
!> \param error ...
!> \par History
!>      06.2005 created [MI]
!>      02.2013 added polarizability section [SL]
!> \author MI
! *****************************************************************************
  SUBROUTINE linres_calculation_low(qs_env,  error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_calculation_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, output_unit
    LOGICAL                                  :: epr_present, failure, &
                                                issc_present, lr_calculation, &
                                                nmr_present, polar_present
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(section_vals_type), POINTER         :: lr_section, prop_section

    CALL timeset(routineN,handle)
    failure = .FALSE.
    lr_calculation = .FALSE.
    nmr_present = .FALSE.
    epr_present = .FALSE.
    issc_present = .FALSE.
    polar_present = .FALSE.

    NULLIFY(dft_control,p_env,linres_control,logger,prop_section,lr_section)
    logger => cp_error_get_logger(error)

    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)
    CALL section_vals_get(lr_section,explicit=lr_calculation,error=error)

    IF(lr_calculation) THEN
       CALL linres_init(lr_section,p_env,qs_env,error=error)
       output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
            extension=".linresLog",error=error)
       CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, &
            linres_control=linres_control,error=error)

       !The type of perturbation has not been defined yet
       linres_control%property = lr_none
       !
       ! We do NMR or EPR, then compute the current response
       prop_section => section_vals_get_subs_vals(lr_section,"NMR",error=error)
       CALL section_vals_get(prop_section,explicit=nmr_present,error=error)
       prop_section => section_vals_get_subs_vals(lr_section,"EPR",error=error)
       CALL section_vals_get(prop_section,explicit=epr_present,error=error)

       IF(nmr_present.OR.epr_present) THEN
          CALL  nmr_epr_linres(linres_control,qs_env,p_env,dft_control, &
                nmr_present,epr_present,output_unit, error=error)
       ENDIF
       !
       ! We do the indirect spin-spin coupling calculation
       prop_section => section_vals_get_subs_vals(lr_section,"SPINSPIN",error=error)
       CALL section_vals_get(prop_section,explicit=issc_present,error=error)

       IF(issc_present) THEN
         CALL issc_linres(linres_control,qs_env, p_env, dft_control, error=error)
       ENDIF

       !
       ! We do the polarizability calculation
       prop_section => section_vals_get_subs_vals(lr_section,"POLAR",error=error)
       CALL section_vals_get(prop_section,explicit=polar_present,error=error)

       IF(polar_present) THEN
         CALL polar_linres(qs_env,p_env,error=error)
       END IF

       ! Other possible LR calculations can be introduced here

       CALL p_env_release(p_env,error=error)

       IF(output_unit>0) THEN
         WRITE (UNIT=output_unit,FMT="(/,T3,A,/,T25,A,/,T3,A,/)")&
            REPEAT("=",77),&
            "ENDED LINRES CALCULATION",&
            REPEAT("=",77)
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
            "PRINT%PROGRAM_RUN_INFO",error=error)
!!    ELSE
!       output_unit = cp_logger_get_default_io_unit(logger)
!       IF(output_unit>0) THEN
!          WRITE (output_unit, "(2X,A)") "",&
!               "-----------------------------------------------------------------------------",&
!               "-              No LR calculation has been specified in the input            -",&
!               "                        cp2k is going to stop, bye bye                       ",&
!               "-----------------------------------------------------------------------------",&
!               ""
!       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE linres_calculation_low

! *****************************************************************************
!> \brief Initialize some general settings like the p_env
!>      Localize the psi0 if required
!> \param lr_section ...
!> \param p_env ...
!> \param qs_env ...
!> \param error ...
!> \par History
!>      06.2005 created [MI]
!> \author MI
!> \note
!>      - The localization should probably be always for all the occupied states
! *****************************************************************************
  SUBROUTINE linres_init(lr_section,p_env,qs_env,error)

    TYPE(section_vals_type), POINTER         :: lr_section
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, output_unit
    LOGICAL                                  :: do_it, failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, rho_ao
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: loc_section

    failure = .FALSE.

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)
    NULLIFY(dft_control, linres_control, loc_section, rho, mos, matrix_ks, rho_ao)

    CPPrecondition(.NOT.ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL linres_control_create(linres_control,error=error)
       CALL set_qs_env(qs_env=qs_env, linres_control=linres_control,error=error)
       CALL linres_control_release(linres_control,error=error)
       CALL get_qs_env(qs_env=qs_env, linres_control=linres_control,&
                       dft_control=dft_control,matrix_ks=matrix_ks,mos=mos,rho=rho,error=error)
       CALL qs_rho_get(rho, rho_ao=rho_ao, error=error)

       ! Localized Psi0 are required when the position operator has to be defined (nmr)
       loc_section =>section_vals_get_subs_vals(lr_section,"LOCALIZE",error=error)
       CALL section_vals_val_get(loc_section,"_SECTION_PARAMETERS_",&
               l_val=linres_control%localized_psi0,error=error)
       IF(linres_control%localized_psi0) THEN
         IF(output_unit>0) THEN
           WRITE (UNIT=output_unit,FMT="(/,T3,A,A)")&
                  "Localization of the ground state orbitals",&
                  " before starting the linear response calculation"
         END IF

         CALL linres_localize(qs_env, linres_control,dft_control%nspins,error=error)

         DO ispin=1,dft_control%nspins
            CALL calculate_density_matrix(mos(ispin)%mo_set,rho_ao(ispin)%matrix,error=error)
         ENDDO
! ** update qs_env%rho
         CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
       END IF

       CALL section_vals_val_get(lr_section,"RESTART",l_val=linres_control%linres_restart,error=error)
       CALL section_vals_val_get(lr_section,"MAX_ITER",i_val=linres_control%max_iter,error=error)
       CALL section_vals_val_get(lr_section,"EPS",r_val=linres_control%eps,error=error)
       CALL section_vals_val_get(lr_section,"RESTART_EVERY",i_val=linres_control%restart_every,error=error)
       CALL section_vals_val_get(lr_section,"PRECONDITIONER",i_val=linres_control%preconditioner_type,error=error)
       CALL section_vals_val_get(lr_section,"ENERGY_GAP",r_val=linres_control%energy_gap,error=error)

       IF(output_unit>0) THEN
         WRITE (UNIT=output_unit,FMT="(/,T3,A,/,T25,A,/,T3,A,/)")&
            REPEAT("=",77),&
            "START LINRES CALCULATION",&
            REPEAT("=",77)

         WRITE (UNIT=output_unit,FMT="(/,T10,A,/)")&
            "Properties to be Calulated:"
         CALL section_vals_val_get(lr_section,"NMR%_SECTION_PARAMETERS_",&
            l_val=do_it,error=error)
         IF(do_it) WRITE (UNIT=output_unit,FMT="(T45,A)") &
             "NMR Chemical Shift"

         IF(linres_control%localized_psi0) WRITE (UNIT=output_unit,FMT="(T2,A,T65,A)")&
            "LINRES|"," LOCALIZED PSI0"

         WRITE(UNIT=output_unit,FMT="(T2,A,T60,A)")&
              "LINRES| Optimization algorithm"," Conjugate Gradients"

         SELECT CASE(linres_control%preconditioner_type)
         CASE(ot_precond_none)
            WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
                 "LINRES| Preconditioner","                NONE"
         CASE(ot_precond_full_single)
            WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
                 "LINRES| Preconditioner","         FULL_SINGLE"
         CASE(ot_precond_full_kinetic)
            WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
                 "LINRES| Preconditioner","        FULL_KINETIC"
         CASE(ot_precond_s_inverse)
            WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
                 "LINRES| Preconditioner","      FULL_S_INVERSE"
         CASE(ot_precond_full_single_inverse)
            WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
                 "LINRES| Preconditioner"," FULL_SINGLE_INVERSE"
         CASE(ot_precond_full_all)
            WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
                 "LINRES| Preconditioner","            FULL_ALL"
         CASE DEFAULT
            CALL stop_program(routineN,moduleN,__LINE__,"Preconditioner NYI")
         END SELECT

         WRITE (UNIT=output_unit,FMT="(T2,A,T72,ES8.1)")&
              "LINRES| EPS",linres_control%eps
         WRITE (UNIT=output_unit,FMT="(T2,A,T72,I8)")&
              "LINRES| MAX_ITER",linres_control%max_iter
       END IF

      !------------------!
      ! create the p_env !
      !------------------!
       CALL p_env_create(p_env, qs_env, orthogonal_orbitals=.TRUE.,linres_control=linres_control,error=error)

       ! update the m_epsilon matrix
       CALL p_env_psi0_changed(p_env,qs_env,error=error)

       ! calculate eigenvectros and eigenvalues of K
       p_env%os_valid = .FALSE.
       p_env%new_preconditioner = .TRUE.
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

  END SUBROUTINE linres_init

! *****************************************************************************
!> \brief ...
!> \param linres_control ...
!> \param qs_env ...
!> \param p_env ...
!> \param dft_control ...
!> \param nmr_present ...
!> \param epr_present ...
!> \param output_unit ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE nmr_epr_linres(linres_control,qs_env,p_env,dft_control,nmr_present,epr_present,output_unit, error)

    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(dft_control_type), POINTER          :: dft_control
    LOGICAL                                  :: nmr_present, epr_present
    INTEGER                                  :: output_unit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nmr_epr_linres', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iB
    LOGICAL                                  :: do_qmmm
    TYPE(current_env_type)                   :: current_env
    TYPE(epr_env_type)                       :: epr_env
    TYPE(nmr_env_type)                       :: nmr_env

    linres_control%property = lr_current
    IF(.NOT.linres_control%localized_psi0) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            "Are you sure that you want to calculate the chemical "//&
            "shift without localized psi0?")
       CALL linres_localize(qs_env, linres_control,&
            dft_control%nspins,centers_only=.TRUE.,error=error)
    ENDIF
    IF(dft_control%nspins/=2.AND.epr_present) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            "LSD is needed to perform a g tensor calculation!")
    ENDIF
    !
    !Initialize the current environment
    do_qmmm = .FALSE.
    current_env%ref_count=0
    IF(qs_env%qmmm) do_qmmm = .TRUE.
    current_env%do_qmmm = do_qmmm
    !current_env%prop='nmr'
    CALL current_env_init(current_env,qs_env,error=error)
    CALL current_operators(current_env,qs_env,error=error)
    CALL current_response(current_env,p_env,qs_env,error)
    !
    IF(current_env%all_pert_op_done) THEN
       !Initialize the nmr environment
       IF(nmr_present) THEN
          nmr_env%ref_count=0
          CALL nmr_env_init(nmr_env,qs_env,error=error)
       ENDIF
       !
       !Initialize the epr environment
       IF(epr_present) THEN
          epr_env%ref_count=0
          CALL epr_env_init(epr_env,qs_env,error=error)
          CALL epr_g_zke(epr_env,qs_env,error=error)
          CALL epr_nablavks(epr_env,qs_env,error=error)
       ENDIF
       !
       ! Build the rs_gauge if needed
       !CALL current_set_gauge(current_env,qs_env,error=error)
       !
       ! Loop over field direction
       DO iB = 1,3
          !
          ! Build current response and succeptibility
          CALL current_build_current(current_env,qs_env,iB,error=error)
          CALL current_build_chi(current_env,qs_env,iB,error=error)
          !
          ! Compute NMR shift
          IF(nmr_present) THEN
             CALL nmr_shift(nmr_env,current_env,qs_env,iB,error=error)
          ENDIF
          !
          ! Compute EPR
          IF(epr_present) THEN
             CALL epr_ind_magnetic_field(epr_env,current_env,qs_env,iB,error=error)
             CALL epr_g_so(epr_env,current_env,qs_env,iB,error=error)
             CALL epr_g_soo(epr_env,current_env,qs_env,iB,error=error)
          ENDIF
       ENDDO
       !
       ! Finalized the nmr environment
       IF(nmr_present) THEN
          CALL nmr_shift_print(nmr_env,current_env,qs_env,error=error)
          CALL nmr_env_cleanup(nmr_env,error=error)
       ENDIF
       !
       ! Finalized the epr environment
       IF(epr_present) THEN
          CALL epr_g_print(epr_env,qs_env,error=error)
          CALL epr_env_cleanup(epr_env,error=error)
       ENDIF
       !
    ELSE
      IF(output_unit>0) THEN
         WRITE(output_unit,"(T10,A,/T20,A,/)")&
           "CURRENT: Not all responses to perturbation operators could be calculated.",&
           " Hence: NO nmr and NO epr possible."
      END IF
    END IF
    ! Finalized the current environment
    CALL current_env_cleanup(current_env,qs_env,error=error)


  END SUBROUTINE  nmr_epr_linres

! *****************************************************************************
!> \brief ...
!> \param linres_control ...
!> \param qs_env ...
!> \param p_env ...
!> \param dft_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE issc_linres(linres_control,qs_env, p_env, dft_control, error)

    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'issc_linres', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom
    LOGICAL                                  :: do_qmmm
    TYPE(current_env_type)                   :: current_env
    TYPE(issc_env_type)                      :: issc_env

    linres_control%property = lr_current
    IF(.NOT.linres_control%localized_psi0) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            "Are you sure that you want to calculate the chemical "//&
            "shift without localized psi0?")
       CALL linres_localize(qs_env,linres_control,&
            dft_control%nspins,centers_only=.TRUE.,error=error)
    ENDIF
    !
    !Initialize the current environment
    do_qmmm = .FALSE.
    current_env%ref_count=0
    IF(qs_env%qmmm) do_qmmm = .TRUE.
    current_env%do_qmmm = do_qmmm
    !current_env%prop='issc'
    !CALL current_env_init(current_env,qs_env,error=error)
    !CALL current_response(current_env,p_env,qs_env,error)
    !
    !Initialize the issc environment
    issc_env%ref_count=0
    CALL issc_env_init(issc_env,qs_env,error=error)
    !
    ! Loop over atoms
    DO iatom = 1,issc_env%issc_natms
       CALL issc_operators(issc_env,qs_env,iatom,error)
       CALL issc_response(issc_env,p_env,qs_env,error)
       CALL issc_issc(issc_env,qs_env,iatom,error=error)
    ENDDO
    !
    ! Finalized the issc environment
    CALL issc_print(issc_env,qs_env,error)
    CALL issc_env_cleanup(issc_env,error)

  END SUBROUTINE issc_linres


! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE polar_linres(qs_env,p_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'polar_linres', &
      routineP = moduleN//':'//routineN

    TYPE(polar_env_type)                     :: polar_env

    polar_env%ref_count=0
    CALL polar_env_init(polar_env,qs_env,error=error)
    CALL polar_operators(polar_env,qs_env,error=error)
    CALL polar_response(polar_env,p_env,qs_env,error=error)
    CALL polar_polar(polar_env,qs_env,error=error)
    !
    CALL polar_print(polar_env,qs_env,error)
    CALL polar_env_cleanup(polar_env,error)

  END SUBROUTINE polar_linres


END MODULE qs_linres_module

