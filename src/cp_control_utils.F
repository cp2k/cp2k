!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
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

  USE band,                            ONLY: band_structure_dealloc_ref,&
                                             band_structure_preinit,&
                                             band_structure_type
  USE cp_control_types,                ONLY: dft_control_type, xc_grid_type,&
                                             kg_control_type, qs_control_type,&
                                             wave_control_type, tddfpt_control_type, &
                                             pw_grid_option, xc_control_type, &
                                             xc_rho_no_smooth, xc_rho_spline2_smooth,&
                                             xc_rho_spline3_smooth, xc_rho_nn10, &
                                             xc_deriv_pw, xc_deriv_spline2,&
                                             xc_deriv_spline3, xc_deriv_spline2_smooth,&
                                             xc_deriv_spline3_smooth, &
                                             dft_control_create
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error,&
                                             cp_unimplemented_error_nr
  USE cp_fm_types,                     ONLY: cp_fm_p_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE pw_grid_types,                   ONLY: pw_grid_type
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
  USE xc_functionals,                  ONLY: xc_get_functional,&
                                             xc_print_info,&
                                             xc_set
  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: moduleN = "cp_control_utils"

! *** Public subroutines ***

  PUBLIC :: read_dft_control,&
            read_kg_control,&
            read_qs_control,&
            read_tddfpt_control,&
            read_wave_control

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
!!                                                                            !
!!<----------------------------------------------------------------------------

    TYPE(dft_control_type), POINTER           :: dft_control
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(cp_error_type), INTENT(inout), OPTIONAL  :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_dft_control (MODULE dft_types)"

!   *** Local variables ***

    CHARACTER(LEN=60) :: keyword, string, str1, str2, str3
    CHARACTER(LEN=3)  :: test_result
    INTEGER           :: istat, i ,tmp
    LOGICAL           :: section_found
    REAL (KIND = dp)       :: ecut, mesh(3)

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

    dft_control % xc_control % functionals ( 1 ) = "SLATER"
    dft_control % xc_control % functionals ( 2 ) = "PZ"
    dft_control % xc_control % functionals ( 3 ) = "NONE"
    dft_control % xc_control % density_cut = 1.0E-10_dp
    dft_control % xc_control % gradient_cut = 1.0E-8_dp
    dft_control % xc_control % gradient_functionals = .FALSE.
    dft_control % xc_control % crossterms = .FALSE.

    SELECT CASE (globenv%run_type)
    CASE ("WAVEFUNCTION OPTIMIZATION")
      dft_control%forces = .FALSE.
    CASE DEFAULT
      dft_control%forces = .TRUE.
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
          CASE DEFAULT
            CALL stop_parser(routine=routine,message="INVALID_KEYWORD")
          END SELECT
        END IF
      END DO
    END IF

    CALL finish_parser()

  END SUBROUTINE read_dft_control

! *****************************************************************************
!!>----------------------------------------------------------------------------
!!  SECTION: &kg ... &end                                                     !
!!                                                                            !
!!  CUTOFF                 : ecut                                             !
!!  MULTI-GRID             : ngrid e1 e2 e3 ...                               !
!!  EPS_DEFAULT            : eps                                              !
!!  EPS_CORE_CHARGE        : eps                                              !
!!  EPS_GVG_RSPACE         : eps                                              !
!!  EPS_PGF_AUX            : eps                                              !
!!  EPS_PGF_ORB            : eps                                              !
!!  EPS_RHO                : eps                                              !
!!  EPS_RHO_GSPACE         : eps                                              !
!!  EPS_RHO_RSPACE         : eps                                              !
!!  PROGRESSION_FACTOR     : eprog                                            !
!!  RELATIVE_CUTOFF        : relcut                                           !
!!  NGRID_LEVEL            : ngrid                                            !
!!  METHOD                 : [KG_POL,KG_NOPOL]                                !
!!                                                                            !
!!<----------------------------------------------------------------------------

  SUBROUTINE read_kg_control(kg_control,globenv)

!   Purpose: Read the KIM-GORDON control parameters.

!   History: - Creation (24.9.2002)

!   ***************************************************************************

    TYPE(kg_control_type), INTENT(OUT)        :: kg_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_kg_control (MODULE dft_types)"

!   *** Local variables ***

    CHARACTER(LEN=60) :: keyword,string,string2
    CHARACTER(LEN=3)  :: test_result
    REAL(KIND = dp)          :: value
    INTEGER           :: igrid_level,istat,ngrid_level
    LOGICAL           :: section_found
    LOGICAL           :: cutoff_set, multgrid_set, ngrid_set

!   ---------------------------------------------------------------------------

!   *** Load the default values ***

    kg_control%method = "KG_NOPOL"
    kg_control%eps_core_charge = 1.0E-12_dp
    kg_control%eps_pgf_aux = 1.0E-6_dp
    kg_control%eps_pgf_orb = 1.0E-6_dp
    kg_control%eps_ppl = 1.0E-6_dp
    kg_control%eps_rho_gspace = 1.0E-8_dp
    kg_control%eps_rho_rspace = 1.0E-8_dp
    kg_control%eps_gvg_rspace = 1.0E-6_dp
    kg_control%progression_factor = 2.0_dp
    kg_control%relative_cutoff = 25.0_dp
    kg_control%polarization = .FALSE.
    kg_control%pw_grid_opt% spherical = .FALSE.
    kg_control%pw_grid_opt%fullspace = .TRUE. 
    kg_control%rs_distributed = .FALSE.
!
    cutoff_set = .FALSE.
    multgrid_set = .FALSE.
    ngrid_set = .FALSE.

!   *** Read the input section ***

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="KG",&
                      section_found=section_found)

    DO WHILE (section_found.AND.(test_object(newline=.TRUE.) /= "EOS"))
      CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
      SELECT CASE (TRIM(keyword))
      CASE ("DENSITY_CUTOFF","E_CUTOFF","ECUTOFF","CUTOFF","ECUT")
        CALL read_object(value)
        kg_control%cutoff = MAX(0.01_dp,0.5_dp*value)
        cutoff_set = .TRUE.
      CASE ("MULTI_GRID","MULTGRID","MGRID")
        CALL read_object(igrid_level)
        IF ( ngrid_set ) THEN
           IF ( igrid_level /= ngrid_level ) CALL stop_parser(routine,&
              "Inconsistent values for number of multi grids")
        END IF
        ngrid_level = igrid_level
        ngrid_set = .TRUE.
        ALLOCATE (kg_control%e_cutoff(ngrid_level),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"kg_control%e_cutoff",&
                                         ngrid_level*dp_size)
        DO igrid_level=1,ngrid_level
          CALL read_object(value)
          kg_control%e_cutoff(igrid_level) = MAX(0.01_dp,0.5_dp*value)
        END DO
        multgrid_set = .TRUE.
      CASE ("EPS_DEFAULT")
        CALL read_object(value)
        kg_control%eps_core_charge = value/100
        kg_control%eps_pgf_aux = value
        kg_control%eps_pgf_orb = value
        kg_control%eps_ppl = value/100
        kg_control%eps_rho_gspace = value
        kg_control%eps_rho_rspace = value
        kg_control%eps_gvg_rspace = value
      CASE ("EPS_CORE_CHARGE")
        CALL read_object(kg_control%eps_core_charge)
      CASE ("EPS_GVG_RSPACE","EPS_GVG")
        CALL read_object(kg_control%eps_gvg_rspace)
      CASE ("EPS_PGF_AUX")
        CALL read_object(kg_control%eps_pgf_aux)
      CASE ("EPS_PGF_ORB")
        CALL read_object(kg_control%eps_pgf_orb)
      CASE ("EPS_PPL")
        CALL read_object(kg_control%eps_ppl)
      CASE ("EPS_RHO")
        CALL read_object(value)
        kg_control%eps_rho_gspace = value
        kg_control%eps_rho_rspace = value
      CASE ("EPS_RHO_GSPACE")
        CALL read_object(kg_control%eps_rho_gspace)
      CASE ("EPS_RHO_RSPACE")
        CALL read_object(kg_control%eps_rho_rspace)
      CASE ("PROGRESSION_FACTOR","PROFAC")
        CALL read_object(value)
        kg_control%progression_factor = MAX(0.01_dp,value)
      CASE ("RELATIVE_CUTOFF","REL_CUTOFF")
        CALL read_object(value)
        kg_control%relative_cutoff = MAX(10.0_dp,value)
      CASE ("NGRID_LEVEL","NGRID")
        CALL read_object(ngrid_level)
        ngrid_set = .TRUE.
      CASE ("METHOD")
        CALL read_object(string,lower_to_upper=.TRUE.)
        SELECT CASE (string)
        CASE ("KG_POL")
          kg_control%method = "KG_POL"
          kg_control%polarization = .TRUE.
        CASE ("KG_NOPOL")
          kg_control%method = "KG_NOPOL"
          kg_control%polarization = .FALSE.
        CASE DEFAULT
          CALL stop_parser(routine,&
                           "Invalid method <"//TRIM(keyword)//&
                           "> specified for KG")
        END SELECT
      CASE ("PW_GRID")
        CALL read_object(string,lower_to_upper=.TRUE.)
        SELECT CASE (string)
        CASE ("SPHERICAL")
          kg_control%pw_grid_opt%spherical = .TRUE.
          kg_control%pw_grid_opt%fullspace = .FALSE.
        CASE ("NON-SPHERICAL")
          kg_control%pw_grid_opt%spherical = .FALSE.
          CALL read_object(string2,lower_to_upper=.TRUE.)
          SELECT CASE (string2)
          CASE ("FULLSPACE")
            kg_control%pw_grid_opt%fullspace = .TRUE.
          CASE ("HALFSPACE")
            kg_control%pw_grid_opt%fullspace = .FALSE.
          CASE DEFAULT
            CALL stop_parser(routine,&
                           "Invalid method <"//TRIM(keyword)//&
                           "> specified for KG")
          END SELECT

          CASE ("RS_GRID")
            CALL read_object(string,lower_to_upper=.TRUE.)
            SELECT CASE (string)
            CASE ("DISTRIBUTED")
                 kg_control%rs_distributed = .TRUE.
            CASE ("REPLICATED")
                 kg_control%rs_distributed = .FALSE.
            CASE DEFAULT
                CALL stop_parser(routine,"Invalid method <"//TRIM(string)//&
                                         "> specified for RS_GRID")
            END SELECT
        CASE DEFAULT
          CALL stop_parser(routine,&
                           "Invalid method <"//TRIM(keyword)//&
                           "> specified for KG")
        END SELECT
      END SELECT
    END DO

    CALL finish_parser()

!   set defaults for cutoff and multigrids if not yet done
    IF ( .NOT. cutoff_set ) THEN
      IF ( .NOT. multgrid_set ) THEN
        kg_control%cutoff = MAX(0.01_dp,0.5_dp*20._dp)
      ELSE
        kg_control%cutoff = kg_control%e_cutoff(1)
      END IF
    END IF
    IF ( .NOT. ngrid_set ) ngrid_level = 1
    IF ( .NOT. multgrid_set ) THEN
      ALLOCATE (kg_control%e_cutoff(ngrid_level),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"kg_control%e_cutoff",&
                                       ngrid_level*dp_size)
      kg_control%e_cutoff(1) = kg_control%cutoff
      DO igrid_level=2,ngrid_level
        kg_control%e_cutoff(igrid_level) = kg_control%e_cutoff(igrid_level-1) &
            / kg_control%progression_factor
      END DO
    ELSE
      ! set cutoff to smallest value in multgrid available with >= cutoff
      DO igrid_level=ngrid_level,1,-1
        IF ( kg_control%cutoff <= kg_control%e_cutoff(igrid_level) ) THEN
          kg_control%cutoff = kg_control%e_cutoff(igrid_level)
          EXIT
        END IF
        ! set largest grid value to cutoff
        IF ( igrid_level == 1) THEN
          kg_control%cutoff = kg_control%e_cutoff(1)
        END IF
      END DO
    END IF
    ! check that multigrids are ordered
    DO igrid_level=2,ngrid_level
      IF ( kg_control%e_cutoff(igrid_level) > kg_control%e_cutoff(igrid_level-1) ) THEN
         CALL stop_program (routine,"Multi grids not ordered")
      END IF
    END DO

  END SUBROUTINE read_kg_control

! *****************************************************************************
!!>----------------------------------------------------------------------------
!!  SECTION: &qs ... &end                                                     !
!!                                                                            !
!!  CUTOFF                 : ecut                                             !
!!  MULTI-GRID             : ngrid e1 e2 e3 ...                               !
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

!   Purpose: Read the QUICKSTEP control parameters.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_control_type), INTENT(OUT)        :: qs_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_qs_control (MODULE dft_types)"

!   *** Local variables ***

    CHARACTER(LEN=60) :: keyword,string,string2
    CHARACTER(LEN=3)  :: test_result
    REAL(KIND = dp)          :: value
    INTEGER           :: igrid_level,istat,ngrid_level
    LOGICAL           :: section_found
    LOGICAL           :: cutoff_set, multgrid_set, ngrid_set

!   ---------------------------------------------------------------------------

!   *** Load the default values ***

    qs_control%method = "GPW"
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
!
    cutoff_set = .FALSE.
    multgrid_set = .FALSE.
    ngrid_set = .FALSE.

!   *** Read the input section ***

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="QS",&
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
          CASE ("DENSITY_CUTOFF","E_CUTOFF","ECUTOFF","CUTOFF","ECUT")
            CALL read_object(value)
            qs_control%cutoff = MAX(0.01_dp,0.5_dp*value) ! ok input is in Ry, value is in a.u.
            cutoff_set = .TRUE.
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
              qs_control%e_cutoff(igrid_level) = MAX(0.01_dp,0.5_dp*value) ! ok input is in Ry, value is in a.u.
            END DO
            multgrid_set = .TRUE.
          CASE ("EPS_DEFAULT")
            CALL read_object(value)
            ! now tries to give parameters that such that energies are approx accurate with EPS_DEFAULT 
            ! gradients with sqrt(EPS_DEFAULT)
            ! setting EPS_DEFAULT = ESP_SCF**2 should allow for convergence 
            ! the /100 is used to set some of the parameters to a stricter value
            ! because they seem computationally not very expensive
            qs_control%eps_core_charge = value/100.0_dp
            qs_control%eps_pgf_aux = SQRT(value) ! unknown to me
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
            qs_control%progression_factor = MAX(0.01_dp,value)
          CASE ("RELATIVE_CUTOFF","REL_CUTOFF")
            CALL read_object(value)
            qs_control%relative_cutoff = MAX(5.0_dp,0.5*value) ! now input is in Ry, value is in a.u.
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

!   set defaults for cutoff and multigrids if not yet done
    IF ( .NOT. cutoff_set ) THEN
      IF ( .NOT. multgrid_set ) THEN
        qs_control%cutoff = MAX(0.01_dp,0.5_dp*280._dp) ! give a 280 Ry default cutoff. should be very fine
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
      IF ( qs_control%e_cutoff(igrid_level) > qs_control%e_cutoff(igrid_level-1) ) THEN
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

    IMPLICIT NONE

    ! arguments
    TYPE(tddfpt_control_type)                 :: t_control
    TYPE(global_environment_type), INTENT(in) :: glob_env

    ! locals
    LOGICAL                                   :: section_found
    CHARACTER(LEN=60)                         :: keyword, field
    INTEGER                                   :: error, length

    CHARACTER(len=*), PARAMETER               :: routineN = "read_tddfpt_control", &
                                                 routineP = moduleN//"/"//routineN

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

    TYPE(wave_control_type), INTENT(OUT)      :: wave_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_wave_control (MODULE dft_types)"

!   *** Local variables ***

    CHARACTER(LEN=60) :: keyword,string
    LOGICAL           :: section_found

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

END MODULE cp_control_utils
