!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004  CP2K developers group                          !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****s* cp2k/cp_control_utils [1.0] *
!!
!!   NAME
!!     cp_control_utils
!!
!!   FUNCTION
!!     Utilities to set up the control types
!!
!!   SOURCE
!******************************************************************************

MODULE cp_control_utils
  USE cp_control_types,                ONLY: dft_control_create,&
                                             dft_control_type,&
                                             qs_control_type,&
                                             tddfpt_control_type,&
                                             wave_control_type,&
                                             xc_deriv_pw,&
                                             xc_deriv_spline2,&
                                             xc_deriv_spline2_smooth,&
                                             xc_deriv_spline3,&
                                             xc_deriv_spline3_smooth,&
                                             xc_new_f_routine,&
                                             xc_old_f_routine,&
                                             xc_rho_nn10,&
                                             xc_rho_no_smooth,&
                                             xc_rho_spline2_smooth,&
                                             xc_rho_spline3_smooth,&
                                             xc_test_f_routine,&
                                             xc_test_lsd_f_routine
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_check,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_reset,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE qs_parser,                       ONLY: finish_parser,&
                                             read_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE qs_wf_history_types,             ONLY: wfi_linear_p_method_nr,&
                                             wfi_linear_ps_method_nr,&
                                             wfi_linear_wf_method_nr,&
                                             wfi_ps_method_nr,&
                                             wfi_use_guess_method_nr,&
                                             wfi_use_prev_p_method_nr,&
                                             wfi_use_prev_rho_r_method_nr,&
                                             wfi_use_prev_wf_method_nr
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE xc_derivatives,                  ONLY: xc_functional_get_info
  USE xc_functionals,                  ONLY: xc_get_functional

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: moduleN = "cp_control_utils"

! *** Public subroutines ***

  PUBLIC :: read_dft_control,&
            read_qs_control,&
            read_tddfpt_control,&
            read_wave_control, &
            write_dft_control,&
            write_qs_control,&
            write_wave_control

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE read_dft_control(dft_control,globenv, error)

!   Purpose: Read the DFT control parameters.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************
!!>----------------------------------------------------------------------------
!!  SECTION: &dft ... &end                                                    !
!!                                                                            !
!!  functional       [ ... ]                                                  !
!!  x-functional     [ ... ]                                                  !
!!  c-functional     [ ... ]                                                  !
!!  ke-functional    [ ... ]                                                  !
!!  density_cutoff   den_cut                                                  !
!!  gradient_cutoff  grad_cut                                                 !
!!  spin             nspins                                                   !
!!  lsd                                                                       !
!!  grid             plane_wave     [cutoff ecut, mesh nx ny nz]              !
!!  charge           charge                                                   !
!!  xc_smooth_rho    [ NONE, NN10, SPLINE2, SPLINE3 ](in order of smoothness) !
!!  xc_deriv         [ PW, SPLINE3, SPLINE2, SPLINE2_SMOOTH, SPLINE3_SMOOTH ] !
!!                   (in order of growing smoothness)
!!  functional_routine [ old, new, test ]
!!                                                                            !
!!<----------------------------------------------------------------------------

    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE read_dft_control (MODULE dft_types)"

    CHARACTER(LEN=3)                         :: test_result
    CHARACTER(LEN=60)                        :: keyword, str1, str2, string
    INTEGER                                  :: i
    LOGICAL                                  :: section_found, xc_fun_set

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(dft_control)) THEN
       CALL dft_control_create(dft_control, error=error)
    END IF

!   *** Load the default values ***

    dft_control % nspins = 1
    dft_control % xc_rho_smooth_id = xc_rho_no_smooth
    dft_control % xc_deriv_method_id = xc_deriv_pw
    dft_control % derive_function_rho = .FALSE.
    dft_control % integration_grid % pw_default = .TRUE.
    dft_control % integration_grid % pw_cutoff  = -1._dp
    dft_control % integration_grid % pw_mesh    = -1_dp
    dft_control % do_tddfpt_calculation = .FALSE.
    dft_control % charge = 0
    dft_control % multiplicity = 0
    dft_control % use_kinetic_energy_density = .FALSE.

    DO i=1,SIZE(dft_control % xc_control % functionals)
       dft_control % xc_control % functionals ( i ) = ""
    END DO
    dft_control % xc_control % density_cut = 1.0E-10_dp
    dft_control % xc_control % density_smooth_cut_range = 10.0_dp
    dft_control % xc_control % gradient_cut = 1.0E-10_dp
    dft_control % xc_control % gradient_functionals = .FALSE.
    dft_control % xc_control % crossterms = .FALSE.
    dft_control % xc_control % functional_routine = xc_old_f_routine

    SELECT CASE (globenv%run_type)
    CASE ("WAVEFUNCTION OPTIMIZATION","MONTE CARLO")
      dft_control%forces = .FALSE.
    CASE ("MOLECULAR DYNAMICS","GEOMETRY OPTIMIZATION")
      dft_control%forces = .TRUE.
    CASE DEFAULT
      CALL stop_program(routine,"unknown run type "//TRIM(globenv%run_type))
    END SELECT

!   *** Read the input section ***

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="DFT",&
                      section_found=section_found)

    IF (section_found) THEN
      DO
        test_result = test_object()
        IF (test_result /= "EOL") THEN
          CALL read_object(keyword)
          CALL stop_parser(routine=routine,message="INVALID_KEYWORD")
        END IF
        test_result = test_object(newline=.TRUE.)
        IF (test_result == "EOS") THEN
          EXIT
        ELSE IF (test_result == "EOF") THEN
          CALL stop_parser(routine=routine,message="EOF")
        ELSE
          CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
          SELECT CASE (TRIM(keyword))
          CASE ("EXCHANGE-CORRELATION-FUNCTIONAL","XC-FUNCTIONAL","XC-FUN",&
                "XCFUN","FUNCTIONAL","FUN")
            CALL read_object(string,lower_to_upper=.TRUE.)
            CALL xc_get_functional ( 0, string, dft_control % xc_control % functionals, &
                                        dft_control % xc_control % gradient_functionals, &
                                        dft_control % xc_control % crossterms )
          CASE ("EXCHANGE-FUNCTIONAL","X-FUNCTIONAL","X-FUN","XFUN")
            CALL read_object(string,lower_to_upper=.TRUE.)
            CALL xc_get_functional ( 1, string, dft_control % xc_control % functionals, &
                                        dft_control % xc_control % gradient_functionals, &
                                        dft_control % xc_control % crossterms )
          CASE ("CORRELATION-FUNCTIONAL","C-FUNCTIONAL","C-FUN","CFUN")
            CALL read_object(string,lower_to_upper=.TRUE.)
            CALL xc_get_functional ( 2, string, dft_control % xc_control % functionals, &
                                        dft_control % xc_control % gradient_functionals, &
                                        dft_control % xc_control % crossterms )
          CASE ("KINETIC-ENERGY-FUNCTIONAL","KE-FUNCTIONAL","KE-FUN","KEFUN")
            CALL read_object(string,lower_to_upper=.TRUE.)
            CALL xc_get_functional ( 3, string, dft_control % xc_control % functionals, &
                                        dft_control % xc_control % gradient_functionals, &
                                        dft_control % xc_control % crossterms )
          CASE ("DENSITY_CUTOFF")
            CALL read_object(dft_control % xc_control % density_cut)
          CASE ("DENSITY_SMOOTH_CUTOFF_RANGE")
            CALL read_object(dft_control % xc_control % density_smooth_cut_range)
          CASE ("GRADIENT_CUTOFF")
            CALL read_object(dft_control % xc_control % gradient_cut)
          CASE ("GRID")
            CALL read_object(str1,lower_to_upper=.TRUE.)
            SELECT CASE (TRIM(str1))
            CASE ( "PLANE_WAVE", "PW" )
              dft_control % integration_grid % pw_default = .FALSE.
              CALL read_object(str2,lower_to_upper=.TRUE.)
              SELECT CASE (TRIM(str2))
              CASE ( "MESH" )
                CALL read_object(dft_control % integration_grid % pw_mesh(1))
                CALL read_object(dft_control % integration_grid % pw_mesh(2))
                CALL read_object(dft_control % integration_grid % pw_mesh(3))
              CASE ( "CUTOFF" )
                CALL read_object(dft_control % integration_grid % pw_cutoff)
              END SELECT
            END SELECT
          CASE ("MULTIPLICITY","MULTIP")
             CALL read_object(dft_control%multiplicity)
          CASE ("NSPINS","SPINS","SPIN")
            CALL read_object(dft_control%nspins)
          CASE ("SPIN_POLARISATION","SPIN_POLARIZATION","LSD")
            dft_control%nspins = 2
          CASE ("CHARGE")
             CALL read_object(dft_control%charge)
          CASE ("KINETIC_ENERGY_DENSITY","TAU")
             dft_control % use_kinetic_energy_density = .TRUE.
          CASE ("EXCITATIONS")
             CALL read_object(string, lower_to_upper=.TRUE.)
             SELECT CASE (TRIM(string))
             CASE ("TDLR","TDDFPT")
                dft_control % do_tddfpt_calculation = .TRUE.
             END SELECT
          CASE ("XC_SMOOTH_RHO")
             CALL read_object(str1,lower_to_upper=.TRUE.)
             SELECT CASE(str1)
             CASE ("NONE")
                dft_control % xc_rho_smooth_id= xc_rho_no_smooth
             CASE ("SPLINE2")
                dft_control % xc_rho_smooth_id= xc_rho_spline2_smooth
             CASE ("SPLINE3")
                dft_control % xc_rho_smooth_id= xc_rho_spline3_smooth
             CASE ("NN10")
                dft_control % xc_rho_smooth_id= xc_rho_nn10
             CASE default
                CALL stop_program(routine,"unknown XC_SMOOTH_RHO option: "//&
                     TRIM(str1))
             END SELECT
          CASE ("XC_DERIV")
             CALL read_object(str1,lower_to_upper=.TRUE.)
             SELECT CASE(str1)
             CASE ("PW")
                dft_control % xc_deriv_method_id = xc_deriv_pw
             CASE ("SPLINE2")
                dft_control % xc_deriv_method_id = xc_deriv_spline2
             CASE ("SPLINE3")
                dft_control % xc_deriv_method_id = xc_deriv_spline3
             CASE ("SPLINE2_SMOOTH")
                dft_control % xc_deriv_method_id = xc_deriv_spline2_smooth
             CASE ("SPLINE3_SMOOTH")
                dft_control % xc_deriv_method_id = xc_deriv_spline3_smooth
             CASE default
                CALL stop_program(routine,"unknown XC_DERIV option: "//&
                     TRIM(str1))
             END SELECT
          CASE ("DERIVE_FUNCTION_RHO")
             CALL read_object(dft_control%derive_function_rho)
          CASE ("FUNCTIONAL_ROUTINE")
             CALL read_object(str1,lower_to_upper=.TRUE.)
             SELECT CASE(str1)
             CASE ("NEW")
                dft_control % xc_control % functional_routine = &
                     xc_new_f_routine
             CASE ("OLD")
                dft_control % xc_control % functional_routine = &
                     xc_old_f_routine
             CASE ("TEST")
                dft_control % xc_control % functional_routine = &
                     xc_test_f_routine
             CASE ("TEST_LSD")
                dft_control % xc_control % functional_routine = &
                     xc_test_lsd_f_routine
             CASE DEFAULT
                CALL stop_parser(routine=routine,message="INVALID_KEYWORD for FUNCTIONAL_ROUTINE "//str1)
             END SELECT
          CASE DEFAULT
            CALL stop_parser(routine=routine,message="INVALID_KEYWORD "//keyword)
          END SELECT
        END IF
      END DO
    END IF

    xc_fun_set=.FALSE.
    DO i=1,SIZE( dft_control % xc_control % functionals )
       IF (dft_control % xc_control % functionals(i)=="") THEN
          dft_control % xc_control % functionals(i)="NONE"
       ELSE
          xc_fun_set=.TRUE.
       END IF
    END DO
    IF (.NOT.xc_fun_set) THEN
       dft_control % xc_control % functionals ( 1 ) = "SLATER"
       dft_control % xc_control % functionals ( 2 ) = "PZ"
    END IF

    CALL finish_parser()

  END SUBROUTINE read_dft_control

! *****************************************************************************
!!>----------------------------------------------------------------------------
!!  SECTION: &qs ... &end                                                     !
!!                                                                            !
!!  CUTOFF                 : ecut                                             !
!!  MULTI-GRID             : ngrid e1 e2 e3 ...                               !
!!  COMMENSURATE           : grids are commensurate in realspace              !
!!  REALSPACE_MGRIDS       : use realspace techniques for multigriding        !
!!  MAP_PAA                : map diagonal blocks on the highest grid          !
!!  MAP_CONSISTENT         : fully consistent realspace integration           !
!!  EPS_DEFAULT            : eps                                              !
!!  EPS_CORE_CHARGE        : eps                                              !
!!  EPS_GVG_RSPACE         : eps                                              !
!!  EPS_PGF_AUX            : eps                                              !
!!  EPS_PGF_ORB            : eps                                              !
!!  EPS_PPL                : eps                                              !
!!  EPS_PPNL               : eps                                              !
!!  EPS_RHO                : eps                                              !
!!  EPS_RHO_GSPACE         : eps                                              !
!!  EPS_RHO_RSPACE         : eps                                              !
!!  PROGRESSION_FACTOR     : eprog                                            !
!!  RELATIVE_CUTOFF        : relcut                                           !
!!  NGRID_LEVEL            : ngrid                                            !
!!  METHOD                 : [GPW, GAPW]                                      !
!!  PAO                                                                       !
!!  wf_interpolation       : [use_guess,use_prev_p,use_prev_rho_r,            !
!!                            linear_wf,linear_p,linear_ps, ps, use_prev_wf]  !
!!                                                                            !
!!<----------------------------------------------------------------------------

  SUBROUTINE read_qs_control(qs_control,globenv)

!   Purpose: Read the QUICKSTEP/KG control parameters.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_control_type), INTENT(OUT)       :: qs_control
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE read_qs_control (MODULE dft_types)"

    CHARACTER(LEN=3)                         :: test_result
    CHARACTER(LEN=60)                        :: keyword, string, string2
    INTEGER                                  :: igrid_level, istat, &
                                                ngrid_level
    LOGICAL                                  :: cutoff_set, multgrid_set, &
                                                ngrid_set, progression_set, &
                                                section_found
    REAL(KIND=dp)                            :: value

!   ---------------------------------------------------------------------------
!   *** Load the default values ***
! is now set as would result from EPS_DEFAULT=1.0E-10

    qs_control%eps_core_charge = 1.0E-12_dp
    qs_control%eps_pgf_aux = 1.0E-5_dp
    qs_control%eps_pgf_orb = 1.0E-5_dp
    qs_control%eps_ppl = 1.0E-12_dp
    qs_control%eps_ppnl = 1.0E-12_dp
    qs_control%eps_rho_gspace = 1.0E-10_dp
    qs_control%eps_rho_rspace = 1.0E-10_dp
    qs_control%eps_gvg_rspace = 1.0E-5_dp
    qs_control%progression_factor = 3.0_dp
    qs_control%relative_cutoff = 15.0_dp ! to 15 a.u. 30 Ry
    qs_control%gapw = .FALSE.
    qs_control%gpw = .TRUE.
    qs_control%pao = .FALSE.
    qs_control%pw_grid_opt%spherical = .FALSE.
    qs_control%pw_grid_opt%fullspace = .TRUE.
    qs_control%rs_distributed = .FALSE.
    qs_control%wf_interpolation_method_nr = wfi_linear_p_method_nr
    qs_control%commensurate_mgrids = .FALSE.
    qs_control%realspace_mgrids = .FALSE.
    qs_control%map_paa = .FALSE.
    qs_control%map_consistent = .FALSE.
    qs_control%polarization = .FALSE.
!
    cutoff_set = .FALSE.
    multgrid_set = .FALSE.
    ngrid_set = .FALSE.
    progression_set = .FALSE.

!   *** Read the input section ***

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="QS",&
                      section_found=section_found)

    IF (section_found) THEN
       qs_control%method = "GPW"
    ELSE
       CALL finish_parser()
       CALL start_parser(file_name=globenv%input_file_name,&
                        globenv=globenv,&
                        start_section_label="KG",&
                        section_found=section_found)
       qs_control%method = "KG_NOPOL"
    END IF

    IF (section_found) THEN
      DO
        test_result = test_object()
        IF (test_result /= "EOL") THEN
          CALL read_object(keyword)
          CALL stop_parser(routine=routine,message="INVALID_KEYWORD")
        END IF
        test_result = test_object(newline=.TRUE.)
        IF (test_result == "EOS") THEN
          EXIT
        ELSE IF (test_result == "EOF") THEN
          CALL stop_parser(routine=routine,message="EOF")
        ELSE
          CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
          SELECT CASE (TRIM(keyword))
          CASE ("DENSITY_CUTOFF","E_CUTOFF","ECUTOFF","CUTOFF","ECUT")
            CALL read_object(value)
            qs_control%cutoff = MAX(0.01_dp,0.5_dp*value) 
                                        ! ok input is in Ry, value is in a.u.
            cutoff_set = .TRUE.
          CASE ("REALSPACE_MGRIDS","REALSPACE_MGRID")
            qs_control%realspace_mgrids = .TRUE.
          CASE ("COMMENSURATE")
            qs_control%commensurate_mgrids = .TRUE.
          CASE ("MAP_PAA")
            qs_control%map_paa = .TRUE.
          CASE ("MAP_CONSISTENT")
            qs_control%map_consistent = .TRUE.
          CASE ("MULTI_GRID","MULTGRID","MGRID")
            CALL read_object(igrid_level)
            IF ( ngrid_set ) THEN
              IF (igrid_level /= ngrid_level) THEN
                CALL stop_parser(routine,"Inconsistent values for number of "//&
                                         "multi grids")
              END IF
            END IF
            ngrid_level = igrid_level
            ngrid_set = .TRUE.
            ALLOCATE (qs_control%e_cutoff(ngrid_level),STAT=istat)
            IF (istat /= 0) CALL stop_memory(routine,"qs_control%e_cutoff",&
                                             ngrid_level*dp_size)
            DO igrid_level=1,ngrid_level
              CALL read_object(value)
              qs_control%e_cutoff(igrid_level) = MAX(0.01_dp,0.5_dp*value) 
                                                ! ok input is in Ry, value is in a.u.
            END DO
            multgrid_set = .TRUE.
          CASE ("EPS_DEFAULT")
            CALL read_object(value)
            ! now tries to give parameters that such that energies 
            ! are approx accurate with EPS_DEFAULT 
            ! gradients with sqrt(EPS_DEFAULT)
            ! setting EPS_DEFAULT = ESP_SCF**2 should allow for convergence 
            ! the /100 is used to set some of the parameters to a stricter value
            ! because they seem computationally not very expensive
            qs_control%eps_core_charge = value/100.0_dp
            qs_control%eps_pgf_aux = SQRT(value) 
            qs_control%eps_pgf_orb = SQRT(value)
            qs_control%eps_ppl = value/100.0_dp
            qs_control%eps_ppnl = value/100.0_dp
            qs_control%eps_rho_gspace = value
            qs_control%eps_rho_rspace = value
            qs_control%eps_gvg_rspace = SQRT(value)
          CASE ("EPS_CORE_CHARGE")
            CALL read_object(qs_control%eps_core_charge)
          CASE ("EPS_GVG_RSPACE","EPS_GVG")
            CALL read_object(qs_control%eps_gvg_rspace)
          CASE ("EPS_PGF_AUX")
            CALL read_object(qs_control%eps_pgf_aux)
          CASE ("EPS_PGF_ORB")
            CALL read_object(qs_control%eps_pgf_orb)
          CASE ("EPS_PPL")
            CALL read_object(qs_control%eps_ppl)
          CASE ("EPS_PPNL")
            CALL read_object(qs_control%eps_ppnl)
          CASE ("EPS_RHO")
            CALL read_object(value)
            qs_control%eps_rho_gspace = value
            qs_control%eps_rho_rspace = value
          CASE ("EPS_RHO_GSPACE")
            CALL read_object(qs_control%eps_rho_gspace)
          CASE ("EPS_RHO_RSPACE")
            CALL read_object(qs_control%eps_rho_rspace)
          CASE ("PROGRESSION_FACTOR","PROFAC")
            CALL read_object(value)
            progression_set=.TRUE.
            qs_control%progression_factor = MAX(0.01_dp,value)
          CASE ("RELATIVE_CUTOFF","REL_CUTOFF")
            CALL read_object(value)
            qs_control%relative_cutoff = MAX(5.0_dp,0.5*value) 
                                        ! now input is in Ry, value is in a.u.
          CASE ("NGRID_LEVEL","NGRID")
            CALL read_object(ngrid_level)
            ngrid_set = .TRUE.
          CASE ("METHOD")
            CALL read_object(string,lower_to_upper=.TRUE.)
            SELECT CASE (string)
            CASE ("GAPW")
              qs_control%method = "GAPW"
              qs_control%gapw = .TRUE.
            CASE ("GPW")
              qs_control%method = "GPW"
              qs_control%gpw = .TRUE.
            CASE ("KG_POL")
              qs_control%method = "KG_POL"
              qs_control%polarization = .TRUE.
            CASE ("KG_NOPOL")
              qs_control%method = "KG_NOPOL"
              qs_control%polarization = .FALSE.
            CASE ("KG_GPW")
              qs_control%method = "KG_GPW"
              qs_control%polarization = .FALSE.
            CASE DEFAULT
              CALL stop_parser(routine,"Invalid method <"//TRIM(keyword)//&
                                       "> specified for METHOD")
            END SELECT
          CASE ("PAO")
            qs_control%pao = .TRUE.
          CASE ("PW_GRID")
            CALL read_object(string,lower_to_upper=.TRUE.)
            SELECT CASE (string)
            CASE ("SPHERICAL")
              qs_control%pw_grid_opt%spherical = .TRUE.
              qs_control%pw_grid_opt%fullspace = .FALSE.
            CASE ("NON-SPHERICAL")
              qs_control%pw_grid_opt%spherical = .FALSE.
              CALL read_object(string2,lower_to_upper=.TRUE.)
              SELECT CASE (string2)
              CASE ("FULLSPACE")
                qs_control%pw_grid_opt%fullspace = .TRUE.
              CASE ("HALFSPACE")
                qs_control%pw_grid_opt%fullspace = .FALSE.
              CASE DEFAULT
                CALL stop_parser(routine,"Invalid method <"//TRIM(keyword)//&
                                         "> specified for PW_GRID NON-SPHERICAL")
              END SELECT
            CASE DEFAULT
              CALL stop_parser(routine,"Invalid method <"//TRIM(keyword)//&
                                       "> specified for PW_GRID")
            END SELECT
          CASE ("RS_GRID")
            CALL read_object(string,lower_to_upper=.TRUE.)
            SELECT CASE (string)
            CASE ("DISTRIBUTED")
                 qs_control%rs_distributed = .TRUE.
            CASE ("REPLICATED")
                 qs_control%rs_distributed = .FALSE.
            CASE DEFAULT
                CALL stop_parser(routine,"Invalid method <"//TRIM(string)//&
                                         "> specified for RS_GRID")
            END SELECT
         CASE ( 'WF_INTERPOLATION','INTERPOLATION','EXTRAPOLATION' )
            CALL read_object(string,lower_to_upper=.TRUE.)
            SELECT CASE ( string )
            CASE DEFAULT
               CALL stop_parser ( routine, &
                    'unknown wf interpolatiom method '//string )
            CASE ( 'USE_GUESS' )
               qs_control % wf_interpolation_method_nr = &
                    wfi_use_guess_method_nr
            CASE ( 'USE_PREV_P')
               qs_control % wf_interpolation_method_nr = &
                    wfi_use_prev_p_method_nr
            CASE ( 'USE_PREV_WF')
               qs_control % wf_interpolation_method_nr = &
                    wfi_use_prev_wf_method_nr
            CASE ( 'USE_PREV_RHO_R')
               qs_control % wf_interpolation_method_nr = &
                    wfi_use_prev_rho_r_method_nr
            CASE ( 'LINEAR_WF')
               qs_control % wf_interpolation_method_nr = &
                    wfi_linear_wf_method_nr
            CASE ( 'LINEAR_P')
               qs_control % wf_interpolation_method_nr = &
                    wfi_linear_p_method_nr
            CASE ( 'LINEAR_PS')
               qs_control % wf_interpolation_method_nr = &
                    wfi_linear_ps_method_nr
            CASE ( 'PS')
               qs_control % wf_interpolation_method_nr = &
                    wfi_ps_method_nr
               CALL read_object(qs_control % wf_extrapolation_order)
            END SELECT
          CASE DEFAULT
            CALL stop_parser(routine=routine,message="INVALID_KEYWORD "&
                 //keyword)
          END SELECT
        END IF
      END DO

    END IF

    CALL finish_parser()

    IF (qs_control%realspace_mgrids .AND. .NOT. qs_control%commensurate_mgrids) THEN
       CALL stop_program(routine,"Realspace multigrids require the COMMENSURATE keyword")
    ENDIF

    IF (qs_control%commensurate_mgrids) THEN
      IF (progression_set) CALL stop_program(routine,&
           "Do not specify a progression factor for commensurate grids (NYI)")
      qs_control%progression_factor=4.0_dp ! this will result in approx. 
                                           ! correct cutoffs for the multigrids
    ENDIF

    ! set defaults for cutoff and multigrids if not yet done
    IF ( .NOT. cutoff_set ) THEN
      IF ( .NOT. multgrid_set ) THEN
        qs_control%cutoff = MAX(0.01_dp,0.5_dp*280._dp) ! give a 280 Ry default cutoff. 
                                                        ! should be very fine
      ELSE
        qs_control%cutoff = qs_control%e_cutoff(1)
      END IF
    END IF
    IF ( .NOT. ngrid_set ) ngrid_level = 4
    IF ( .NOT. multgrid_set ) THEN
      ALLOCATE (qs_control%e_cutoff(ngrid_level),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"qs_control%e_cutoff",&
                                       ngrid_level*dp_size)
      qs_control%e_cutoff(1) = qs_control%cutoff
      DO igrid_level=2,ngrid_level
        qs_control%e_cutoff(igrid_level) = qs_control%e_cutoff(igrid_level-1) &
            / qs_control%progression_factor
      END DO
    ELSE
      IF (qs_control%commensurate_mgrids) THEN
          CALL stop_program(routine,&
               "Do not specify cutoffs for the commensurate grids (NYI)")
      ENDIF
      ! set cutoff to smallest value in multgrid available with >= cutoff
      DO igrid_level=ngrid_level,1,-1
        IF ( qs_control%cutoff <= qs_control%e_cutoff(igrid_level) ) THEN
          qs_control%cutoff = qs_control%e_cutoff(igrid_level)
          EXIT
        END IF
        ! set largest grid value to cutoff
        IF ( igrid_level == 1) THEN
          qs_control%cutoff = qs_control%e_cutoff(1)
        END IF
      END DO
    END IF
    ! check that multigrids are ordered
    DO igrid_level=2,ngrid_level
      IF ( qs_control%e_cutoff(igrid_level) > &
           qs_control%e_cutoff(igrid_level-1) ) THEN
         CALL stop_program (routine,"Multi grids not ordered")
      END IF
    END DO

  END SUBROUTINE read_qs_control

! *****************************************************************************

!!>----------------------------------------------------------------------------
!!  SECTION: &tddfpt ... &end                                                 !
!!                                                                            !
!!  nev                    : (integer) number of electrons                    !
!!  nreortho               : (integer) number of reorthogonalization steps    !
!!  kernel                 : (y/n) compute the kernel? (for debugging)        !
!!  res_e_type             : (singlets/triplets) for restricted calculation   !
!!  functional             : combined functional to use in kernel             !
!!  exchange-functional    : exchange-functional to use in kernel             !
!!  correlation-functional : correlation-functional to use in kernel          !
!!  kinetic-energy-functional : kinetic-energy-functional to use in kernel    !
!!                                                                            !
!!<----------------------------------------------------------------------------

  SUBROUTINE read_tddfpt_control(t_control, glob_env)


    TYPE(tddfpt_control_type)                :: t_control
    TYPE(global_environment_type), &
      INTENT(in)                             :: glob_env

    CHARACTER(len=*), PARAMETER :: routineN = "read_tddfpt_control", &
      routineP = moduleN//"/"//routineN

    CHARACTER(LEN=60)                        :: field, keyword
    LOGICAL                                  :: section_found

! -------------------------------------------------------------------
!----------------!
! do the parsing !
!----------------!

    section_found = .FALSE.
        CALL start_parser(file_name=glob_env%input_file_name,&
                      globenv=glob_env,&
                      start_section_label="TDDFPT",&
                      section_found=section_found)

    DO WHILE (section_found.AND.(test_object(newline=.TRUE.) /= "EOS"))

       CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)

       SELECT CASE (TRIM(keyword))

       CASE ("NEV", "N_EV", "EV")
          CALL read_object(t_control%n_ev)
          IF (t_control%n_ev < 1) &
               CALL stop_parser(routineP, "number of eignvalues must be >= 1")

       CASE ("NREORTHO", "REORTHO", "N_REORTHO", &
             "REORTHOGONALIZATIONS")
          CALL read_object(t_control%n_reortho)
          IF (t_control%n_reortho < 1) &
               CALL stop_parser(routineP, "no. of reortho. steps must  be >= 1")

       CASE ("KERNEL")
          CALL read_object(field, lower_to_upper=.TRUE.)
          SELECT CASE (TRIM(field))
          CASE ("Y", "YES") 
             t_control%do_kernel = .TRUE.
          CASE ("N", "NO")
             t_control%do_kernel = .FALSE.
          CASE DEFAULT
             t_control%do_kernel = .TRUE.
          END SELECT

       CASE ("RES_E_TYPE", &
             "RESTRICTED_EXCITATIONS_TYPE")
          CALL read_object(field, lower_to_upper=.TRUE.)
          SELECT CASE (TRIM(field))
          CASE ("SINGLET", "SINGLETS", "S")
             t_control%restricted_excitations_type = "SINGLET"
          CASE ("TRIPLET", "TRIPLETS", "T")
             t_control%restricted_excitations_type = "TRIPLET"
          CASE DEFAULT
             t_control%restricted_excitations_type = "SINGLET"
          END SELECT

       CASE ("EXCHANGE-CORRELATION-FUNCTIONAL","XC-FUNCTIONAL","XC-FUN",&
             "XCFUN","FUNCTIONAL","FUN")
         CALL read_object(field, lower_to_upper=.TRUE.)
         CALL xc_get_functional (0, field, t_control%xc_control%functionals, &
                                 t_control%xc_control%gradient_functionals, &
                                 t_control%xc_control%crossterms )
       CASE ("EXCHANGE-FUNCTIONAL","X-FUNCTIONAL","X-FUN","XFUN")
         CALL read_object(field, lower_to_upper=.TRUE.)
         CALL xc_get_functional (1, field, t_control%xc_control%functionals, &
                                 t_control%xc_control%gradient_functionals, &
                                 t_control%xc_control%crossterms )
       CASE ("CORRELATION-FUNCTIONAL","C-FUNCTIONAL","C-FUN","CFUN")
         CALL read_object(field, lower_to_upper=.TRUE.)
         CALL xc_get_functional (2, field, t_control%xc_control%functionals, &
                                 t_control%xc_control%gradient_functionals, &
                                 t_control%xc_control%crossterms )
       CASE ("KINETIC-ENERGY-FUNCTIONAL","KE-FUNCTIONAL","KE-FUN","KEFUN")
         CALL read_object(field, lower_to_upper=.TRUE.)
         CALL xc_get_functional (3, field, t_control%xc_control%functionals, &
                                 t_control%xc_control%gradient_functionals, &
                                 t_control%xc_control%crossterms )
       CASE DEFAULT
          CALL stop_parser(routineP, "unrecognized option <"//TRIM(keyword)//">")

      END SELECT

   END DO
   
   CALL finish_parser()

  END SUBROUTINE read_tddfpt_control

!   ***************************************************************************

  SUBROUTINE read_wave_control(wave_control,globenv)

!   Purpose: Read the WAVE control parameters.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************
!!>----------------------------------------------------------------------------
!!  SECTION: &wave ... &end                                                   !
!!                                                                            !
!!  ecut             e_cutoff_wf                                              !
!!  dual             e_cutoff_dual                                            !
!!                                                                            !
!!<----------------------------------------------------------------------------

    TYPE(wave_control_type), INTENT(OUT)     :: wave_control
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE read_wave_control (MODULE dft_types)"

    CHARACTER(LEN=60)                        :: keyword
    LOGICAL                                  :: section_found

!   ---------------------------------------------------------------------------
!   *** Load the default values ***

    wave_control%e_cutoff_wf = -1.0_dp
    wave_control%e_cutoff_dual = 4.0_dp

!   *** Read the input section ***

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="WAVE",&
                      section_found=section_found)

    DO WHILE (section_found.AND.(test_object(newline=.TRUE.) /= "EOS"))
      CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
      SELECT CASE (TRIM(keyword))
      CASE ("CUTOFF_DUAL","DUAL")
        CALL read_object(wave_control%e_cutoff_dual)
      CASE ("CUTOFF_WF","CUTOFF")
        CALL read_object(wave_control%e_cutoff_wf)
      END SELECT
    END DO

    CALL finish_parser()

  END SUBROUTINE read_wave_control

! *****************************************************************************

  SUBROUTINE write_dft_control(dft_control,globenv, error)

    !   Purpose: Write the DFT control parameters to the output unit.

    !   ***************************************************************************

    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=160)                       :: reference
    CHARACTER(len=20)                        :: tmpStr
    INTEGER                                  :: i, output_unit
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: sub_error

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN
    failure=.FALSE.

    IF (globenv%print%dft_control_parameters) THEN

       output_unit = globenv%scr

       IF (dft_control%forces) THEN
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T78,A)")&
               "DFT| Force calculation","YES"
       ELSE
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T78,A)")&
               "DFT| Force calculation"," NO"
       END IF

       IF ( dft_control % nspins == 1 ) THEN
          !do nothing
       ELSE IF ( dft_control % nspins == 2 ) THEN
          WRITE ( output_unit, '( A, T60, A )' ) " DFT| ","spin polarized system"
       ELSE
          WRITE ( output_unit, '( A, T60, A, T77, I4 )' ) " DFT| ","sets of spin states:",&
               dft_control % nspins
       END IF

       IF ( dft_control % charge /= 0 ) THEN
          WRITE ( output_unit, '( A, T60, A, T77, I4 )' ) " DFT| ","charge:",&
               dft_control % charge
       END IF

       WRITE (UNIT=output_unit,FMT="(T2,A,T67,E14.5)")&
            "DFT| Cutoffs: density ",dft_control%xc_control%density_cut,&
            "DFT|          gradient",dft_control%xc_control%gradient_cut,&
            "DFT|          cutoff_smoothing_range",dft_control%xc_control%density_smooth_cut_range

       tmpStr=xc_get_rho_smooth_label(dft_control%xc_rho_smooth_id)
       WRITE ( output_unit, '( A, T61, A )' ) &
            " DFT| XC density smoothing ",ADJUSTR(tmpStr)
       tmpStr=xc_get_deriv_method_label(dft_control%xc_deriv_method_id)
       WRITE ( output_unit, '( A, T61, A )' ) &
            " DFT| XC derivatives ",ADJUSTR(tmpStr)
       IF ( dft_control % derive_function_rho )  &
            WRITE ( output_unit, '( A, T79, A )' ) &
            " DFT| Density function derivative ","ON"

       IF ( dft_control % integration_grid % pw_default ) THEN
          WRITE ( output_unit, '( A, T74, A )' ) " DFT| XC mesh (PW density)", "default"
       ELSE
          IF ( dft_control % integration_grid % pw_mesh ( 1 ) > 0 ) &
               WRITE ( output_unit, '( A, T51, 3I10 )' ) &
               " DFT| XC mesh (PW density)", dft_control % integration_grid % pw_mesh
          IF ( dft_control % integration_grid % pw_cutoff > 0 ) &
               WRITE ( output_unit, '( A, T67, F10.1, "a.u." )' ) &
               " DFT| XC mesh (PW density)", dft_control % integration_grid % pw_cutoff
       END IF

       !TC      CALL xc_set(dft_control%xc_control%density_cut, dft_control%xc_control%gradient_cut)
       !      CALL xc_print_info ( output_unit, dft_control % xc_control % functionals, &
       !                          dft_control % nspins /= 1 )
       IF (.NOT.failure) THEN
          CALL cp_error_init(sub_error,stop_level=cp_fatal_level, &
               print_level=cp_fatal_level,template_error=error)
          DO i=1,SIZE(dft_control % xc_control % functionals)
             IF ( dft_control % xc_control % functionals(i) /= "NONE") THEN
                CALL xc_functional_get_info(&
                     functional=dft_control % xc_control % functionals(i),&
                     lsd=(dft_control % nspins /= 1),&
                     reference=reference, error=sub_error)
                CALL cp_error_check(error,failure)

                WRITE (output_unit, fmt="(' FUNCTIONAL| ',a,':')")&
                     TRIM(dft_control % xc_control % functionals(i))
                IF (failure) THEN
                   WRITE (output_unit, fmt="(' FUNCTIONAL| ***UNKNOWN***')")
                   failure=.FALSE.
                   CALL cp_error_reset(sub_error)
                ELSE
                   WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)")&
                        TRIM(reference)
                END IF
             END IF
          END DO
          CALL cp_error_dealloc_ref(sub_error)
       END IF
    END IF
  END SUBROUTINE write_dft_control

! *****************************************************************************

  SUBROUTINE write_qs_control(qs_control,globenv)

!   Purpose: Write the DFT control parameters to the output unit.

!   ***************************************************************************

    TYPE(qs_control_type), INTENT(IN)        :: qs_control
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: igrid_level, ngrid_level, &
                                                output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%dft_control_parameters) THEN
      ngrid_level = SIZE(qs_control%e_cutoff)
      output_unit = globenv%scr
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T71,A)")&
        "QS| Method:",ADJUSTR(qs_control%method)
      IF ( qs_control% pw_grid_opt % spherical ) THEN
         WRITE (UNIT=output_unit,FMT="(T2,A,T61,A)")&
            "QS| Density plane wave grid type"," SPHERICAL HALFSPACE"
      ELSE IF ( qs_control% pw_grid_opt % fullspace ) THEN
         WRITE (UNIT=output_unit,FMT="(T2,A,T57,A)")&
            "QS| Density plane wave grid type"," NON-SPHERICAL FULLSPACE"
      ELSE
         WRITE (UNIT=output_unit,FMT="(T2,A,T57,A)")&
            "QS| Density plane wave grid type"," NON-SPHERICAL HALFSPACE"
      END IF
      WRITE (UNIT=output_unit,FMT="(T2,A,T71,I10)")&
        "QS| Number of grid levels:",SIZE(qs_control%e_cutoff)
      IF (ngrid_level == 1) THEN
        WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
          "QS| Density cutoff [a.u.]:",qs_control%e_cutoff(1)
      ELSE
        WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
          "QS| Density cutoff [a.u.]:",qs_control%cutoff
        WRITE (UNIT=output_unit,FMT="(T2,A)") "QS| Using commensurate multigrids"
        WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
             "QS| Multi grid cutoff [a.u.]: 1) grid level",qs_control%e_cutoff(1)
        WRITE (UNIT=output_unit,FMT="(T2,A,I3,A,T71,F10.1)")&
                ("QS|                         ",igrid_level,") grid level",&
              qs_control%e_cutoff(igrid_level),&
              igrid_level=2,SIZE(qs_control%e_cutoff))
      END IF
      IF (qs_control%pao) THEN
         WRITE (UNIT=output_unit,FMT="(T2,A)") "QS| PAO active"
      END IF
      WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
        "QS| Grid level progression factor:",qs_control%progression_factor
      WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
        "QS| Relative density cutoff [a.u.]:",qs_control%relative_cutoff
      IF (qs_control%map_consistent) THEN
         WRITE (UNIT=output_unit,FMT="(T2,A)")&
            "QS| Consistent realspace mapping and integration "
      ENDIF
      WRITE (UNIT=output_unit,FMT="(T2,A,T73,ES8.1)")&
        "QS| Interaction thresholds: eps_pgf_orb:",&
        qs_control%eps_pgf_orb,&
        "QS|                         eps_core_charge:",&
        qs_control%eps_core_charge,&
        "QS|                         eps_rho_gspace:",&
        qs_control%eps_rho_gspace,&
        "QS|                         eps_rho_rspace:",&
        qs_control%eps_rho_rspace,&
        "QS|                         eps_gvg_rspace:",&
        qs_control%eps_gvg_rspace,&
        "QS|                         eps_ppl:",&
        qs_control%eps_ppl,&
        "QS|                         eps_ppnl:",&
        qs_control%eps_ppnl
    END IF

  END SUBROUTINE write_qs_control

! *****************************************************************************

  SUBROUTINE write_wave_control(wave_control,globenv)

!   Purpose: Write the DFT control parameters to the output unit.

!   ***************************************************************************

    TYPE(wave_control_type), INTENT(IN)      :: wave_control
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%dft_control_parameters) THEN
      output_unit = globenv%scr
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T71,F10.3)")&
        "WAVE| Plane wave cutoff [a.u.]:",wave_control%e_cutoff_wf
      WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.2)")&
        "WAVE| Density cutoff (DUAL):",wave_control%e_cutoff_dual
     END IF

  END SUBROUTINE write_wave_control

! *****************************************************************************

!!****f* cp_control_types/xc_get_rho_smooth_label [1.0] *
!!
!!   NAME
!!     xc_get_rho_smooth_label
!!
!!   FUNCTION
!!     returns a string that describes the smoothing of rho
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - xc_rho_smooth_id: the id that represent the smoothing
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
FUNCTION xc_get_rho_smooth_label(xc_rho_smooth_id,error) RESULT(res)
    INTEGER, INTENT(in)                      :: xc_rho_smooth_id
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    CHARACTER(len=10)                        :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_get_rho_smooth_label', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  SELECT CASE(xc_rho_smooth_id)
  CASE (xc_rho_no_smooth)
     res="NONE"
  CASE (xc_rho_spline2_smooth)
     res="SPLINE2"
  CASE (xc_rho_spline3_smooth)
     res="SPLINE3"
  CASE (xc_rho_nn10)
     res="NN10"
  CASE default
     WRITE (res,"('UKN',i6)") xc_rho_smooth_id
  END SELECT
END FUNCTION xc_get_rho_smooth_label
!***************************************************************************

!!****f* cp_control_types/xc_get_deriv_method_label [1.0] *
!!
!!   NAME
!!     xc_get_deriv_method_label
!!
!!   FUNCTION
!!     returns a string that describes the derivative used in the xc 
!!     calculation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - xc_deriv_method_id: the id that represent the derivative method
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2002 created [fawzi]
!!
!!*** **********************************************************************
FUNCTION xc_get_deriv_method_label(xc_deriv_method_id,error) RESULT(res)
    INTEGER, INTENT(in)                      :: xc_deriv_method_id
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    CHARACTER(len=20)                        :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_get_deriv_method_label', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  SELECT CASE(xc_deriv_method_id)
  CASE (xc_deriv_pw)
     res="PW"
  CASE (xc_deriv_spline2)
     res="SPLINE2"
  CASE (xc_deriv_spline3)
     res="SPLINE3"
  CASE (xc_deriv_spline2_smooth)
     res="SPLINE2_SMOOTH"
  CASE (xc_deriv_spline3_smooth)
     res="SPLINE3_SMOOTH"
  CASE default
     WRITE (res,"('UKN',i6)") xc_deriv_method_id
  END SELECT
END FUNCTION xc_get_deriv_method_label
!***************************************************************************

END MODULE cp_control_utils
