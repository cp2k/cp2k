!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****s* cp2k/dft_types [1.0] *
!!
!!   NAME
!!     dft_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH,MK
!!
!!   MODIFICATION HISTORY
!!     - qs_control_type extended (18.01.2002,MK)
!!     - added pao (18.04.02,fawzi)
!!     - new xc code (31.05.02,jgh)
!!     - qs section revised (29.08.02,jgh)
!!     - kg section added (24.09.02,gt)
!!     - pw_grid options added (22.02.03,jgh)
!!
!!   SOURCE
!******************************************************************************

MODULE dft_types
  USE band,                            ONLY: band_structure_dealloc_ref,&
                                             band_structure_preinit,&
                                             band_structure_type
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

  TYPE xc_grid_type
     LOGICAL                   :: pw_default
!FM     TYPE ( pw_grid_type )     :: pw_grid
     INTEGER                   :: pw_mesh ( 3 )
     REAL (KIND = dp)               :: pw_cutoff
  END TYPE xc_grid_type

! *** Control parameters for pw grids ***

  TYPE pw_grid_option
     LOGICAL :: spherical
     LOGICAL :: fullspace
  END TYPE pw_grid_option

! *** Control parameters for a KIM-GORDON calculation ***

  TYPE kg_control_type
    CHARACTER(LEN=10)               :: method
    REAL(KIND = dp)                 :: eps_core_charge,&
                                       eps_pgf_aux,&
                                       eps_pgf_orb,&
                                       eps_ppl,&
                                       eps_rho_gspace,&
                                       eps_rho_rspace,&
                                       eps_gvg_rspace,&
                                       progression_factor,&
                                       relative_cutoff
    LOGICAL                         :: polarization
    REAL(KIND = dp)                        :: cutoff
    REAL(KIND = dp), DIMENSION(:), POINTER :: e_cutoff
    TYPE ( pw_grid_option )         :: pw_grid_opt
    LOGICAL                         :: rs_distributed
  END TYPE kg_control_type

! *** Control parameters for a QUICKSTEP calculation ***

! eps_pgf_orb: Cutoff value for the interaction of the primitive Gaussian-type
!              functions (primitive basis functions).

  TYPE qs_control_type
!MK PRIVATE would be nice
    CHARACTER(LEN=10)               :: method
    REAL(KIND = dp)                 :: eps_core_charge,&
                                       eps_pgf_aux,&
                                       eps_pgf_orb,&
                                       eps_ppl,&
                                       eps_ppnl,&
                                       eps_rho_gspace,&
                                       eps_rho_rspace,&
                                       eps_gvg_rspace,&
                                       progression_factor,&
                                       relative_cutoff
    LOGICAL                         :: gapw,gpw,pao
    REAL(KIND = dp)                        :: cutoff
    REAL(KIND = dp), DIMENSION(:), POINTER :: e_cutoff
    TYPE ( pw_grid_option )         :: pw_grid_opt
    LOGICAL                         :: rs_distributed
    INTEGER                         :: wf_interpolation_method_nr
    INTEGER                         :: wf_extrapolation_order
  END TYPE qs_control_type

! *** Control parameters for a WAVE calculation ***

  TYPE wave_control_type
     REAL(KIND = dp) :: e_cutoff_wf
     REAL(KIND = dp) :: e_cutoff_dual
  END TYPE wave_control_type

! *** Control parameters for a dft xc functionals ***

  TYPE xc_control_type
     CHARACTER(LEN=40), DIMENSION(3) :: functionals
     LOGICAL, DIMENSION(3)           :: gradient_functionals, crossterms
     REAL(KIND = dp)                        :: density_cut,&
                                        gradient_cut
  END TYPE xc_control_type

  ! *** Control parameters for a TIME-DEPENDENT PERTURBATION calculation ***

!!****s* dft_types/tddfpt_control_type
!!
!! NAME
!!   tddfpt_control_type
!!
!! FUNCTION
!!   A type that holds controling information for a time-dependent
!!   perturbation calculation.
!!
!! ATTRIBUTES
!!   - n_ev       : number of eigenvalues to calculate
!!   - n_reortho  : how many time to reorthogonalize (in the lanczos algorithm)
!!   - do_kernel  : wether to evaluate the kernel (this is a debugging option)
!!   - restricted_excitations_type : { SINGLET | TRIPLET } which excitations
!!                  to calculate
!!   - functional : the functional to use in the perturbation kernel
!!   - gradient_functionals:
!!   - crossterms :
!!   - lumos_eigenvalues : holds the eigenvalues of the lumos (if calculated in QS)
!!   - lumos      : holds the first few lumos (if calculated in QS)
!!
!! NOTES  
!!   The lumos are helpfull in choosing a initial vector for the TDDFPT 
!!   calculation, since they can be used to construct the solutions of the
!!   TDDFPT operator without the perturbation kernel.
!!  
!!***
  TYPE tddfpt_control_type
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: lumos
     TYPE(xc_control_type)                     :: xc_control
     INTEGER                                   :: n_ev
     INTEGER                                   :: n_reortho
     LOGICAL                                   :: do_kernel
     CHARACTER(LEN=7)                          :: restricted_excitations_type
     REAL(KIND = dp), DIMENSION(:,:), POINTER    :: lumos_eigenvalues
  END TYPE tddfpt_control_type

  ! *** Control parameters for a DFT calculation ***

  TYPE dft_control_type
     TYPE(band_structure_type)       :: band
     TYPE(kg_control_type)           :: kg_control
     TYPE(qs_control_type)           :: qs_control
     TYPE(tddfpt_control_type)       :: tddfpt_control
     TYPE(wave_control_type)         :: wave_control
     TYPE(xc_grid_type)              :: integration_grid
     TYPE(xc_control_type)           :: xc_control
     INTEGER                         :: nspins,&
                                        charge,&
                                        multiplicity,&
                                        xc_rho_smooth_id,&
                                        xc_deriv_method_id,&
                                        ref_count,&
                                        id_nr
     LOGICAL                         :: forces,&
                                        do_tddfpt_calculation, &
                                        derive_function_rho, &
                                        use_kinetic_energy_density
  END TYPE dft_control_type

  CHARACTER(LEN=*), PARAMETER :: moduleN = "dft_types"

  INTEGER, PARAMETER, PUBLIC :: xc_rho_no_smooth=0, xc_rho_spline2_smooth=1,&
       xc_rho_spline3_smooth=2, xc_rho_nn10=3
  INTEGER, PARAMETER, PUBLIC :: xc_deriv_pw=0, xc_deriv_spline2=1,&
       xc_deriv_spline3=2, xc_deriv_spline2_smooth=3,&
       xc_deriv_spline3_smooth=4
  INTEGER, SAVE :: last_dft_control_id=0

! *** Public data types ***

  PUBLIC :: dft_control_type, xc_grid_type,&
            kg_control_type, qs_control_type,&
            wave_control_type, tddfpt_control_type, &
            pw_grid_option, xc_control_type

! *** Public subroutines ***

  PUBLIC :: read_dft_control,&
            read_kg_control,&
            read_qs_control,&
            read_tddfpt_control,&
            read_wave_control,&
            write_dft_control,&
            write_kg_control,&
            write_qs_control,&
            write_wave_control,&
            dft_control_retain,&
            dft_control_release,&
            dft_control_create

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

  SUBROUTINE write_dft_control(dft_control,globenv)

!   Purpose: Write the DFT control parameters to the output unit.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************

    TYPE(dft_control_type), POINTER           :: dft_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: output_unit
    CHARACTER(len=20) :: tmpStr

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

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
        "DFT|          gradient",dft_control%xc_control%gradient_cut

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

      CALL xc_set(dft_control%xc_control%density_cut, dft_control%xc_control%gradient_cut)
      CALL xc_print_info ( output_unit, dft_control % xc_control % functionals, &
                           dft_control % nspins /= 1 )

    END IF

  END SUBROUTINE write_dft_control

! *****************************************************************************

  SUBROUTINE write_kg_control(kg_control,globenv)

!   Purpose: Write the DFT control parameters to the output unit.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************

    TYPE(kg_control_type), INTENT(IN)         :: kg_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: igrid_level,ngrid_level,output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%dft_control_parameters) THEN
      ngrid_level = SIZE(kg_control%e_cutoff)
      output_unit = globenv%scr
      WRITE (UNIT=output_unit,FMT="(/,T2,A,T71,A)")&
        "KG| Method:",ADJUSTR(kg_control%method)
      WRITE (UNIT=output_unit,FMT="(T2,A,T71,I10)")&
        "KG| Number of grid levels:",SIZE(kg_control%e_cutoff)
      IF ( kg_control% pw_grid_opt % spherical ) THEN
         WRITE (UNIT=output_unit,FMT="(T2,A,T61,A)")&
            "QS| Density plane wave grid type"," SPHERICAL HALFSPACE"
      ELSE IF ( kg_control% pw_grid_opt % fullspace ) THEN
         WRITE (UNIT=output_unit,FMT="(T2,A,T57,A)")&
            "QS| Density plane wave grid type"," NON-SPHERICAL FULLSPACE"
      ELSE
         WRITE (UNIT=output_unit,FMT="(T2,A,T57,A)")&
            "QS| Density plane wave grid type"," NON-SPHERICAL HALFSPACE"
      END IF
      IF (ngrid_level == 1) THEN
        WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
          "KG| Density cutoff [a.u.]:",kg_control%e_cutoff(1)
      ELSE
        WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
          "KG| Density cutoff [a.u.]:",kg_control%cutoff
        WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
          "KG| Multi grid cutoff [a.u.]: 1) grid level",kg_control%e_cutoff(1)
        WRITE (UNIT=output_unit,FMT="(T2,A,I3,A,T71,F10.1)")&
          ("KG|                         ",igrid_level,") grid level",&
           kg_control%e_cutoff(igrid_level),&
           igrid_level=2,SIZE(kg_control%e_cutoff))
      END IF
      WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
        "KG| Grid level progression factor:",kg_control%progression_factor
      WRITE (UNIT=output_unit,FMT="(T2,A,T71,F10.1)")&
        "KG| Relative density cutoff [a.u.]:",kg_control%relative_cutoff
      WRITE (UNIT=output_unit,FMT="(T2,A,T73,ES8.1)")&
        "KG| Interaction thresholds: eps_pgf_orb:",&
        kg_control%eps_pgf_orb,&
        "KG|                         eps_core_charge:",&
        kg_control%eps_core_charge,&
        "KG|                         eps_rho_gspace:",&
        kg_control%eps_rho_gspace,&
        "KG|                         eps_rho_rspace:",&
        kg_control%eps_rho_rspace,&
        "KG|                         eps_gvg_rspace:",&
        kg_control%eps_gvg_rspace  
    END IF

  END SUBROUTINE write_kg_control

! *****************************************************************************

  SUBROUTINE write_qs_control(qs_control,globenv)

!   Purpose: Write the DFT control parameters to the output unit.

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_control_type), INTENT(IN)         :: qs_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: igrid_level,ngrid_level,output_unit

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

!   History: - Creation (21.01.2002,MK)

!   ***************************************************************************

    TYPE(wave_control_type), INTENT(IN)       :: wave_control
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: output_unit

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

!!****f* dft_types/xc_get_rho_smooth_label [1.0] *
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
  INTEGER, INTENT(in) :: xc_rho_smooth_id
  CHARACTER(len=10) :: res
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  CHARACTER(len=*), PARAMETER :: routineN='xc_get_rho_smooth_label',&
        routineP=moduleN//':'//routineN
  LOGICAL :: failure

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

!!****f* dft_types/xc_get_deriv_method_label [1.0] *
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
  INTEGER, INTENT(in) :: xc_deriv_method_id
  CHARACTER(len=20) :: res
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  CHARACTER(len=*), PARAMETER :: routineN='xc_get_deriv_method_label',&
        routineP=moduleN//':'//routineN
  LOGICAL :: failure
  
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

!!****f* dft_types/dft_control_create [1.0] *
!!
!!   NAME
!!     dft_control_create
!!
!!   SYNOPSIS
!!     Subroutine dft_control_create(dft_control, error)
!!       Type(dft_control_type), Pointer:: dft_control
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine dft_control_create
!!
!!   FUNCTION
!!     allocates and perform a very basic initialization
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - dft_control: the object to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE dft_control_create(dft_control, error)
  TYPE(dft_control_type), POINTER :: dft_control
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='dft_control_create',&
       routineP=moduleN//':'//routineN
  INTEGER :: stat

  failure=.FALSE.

  CPPrecondition(.NOT.ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)  
  IF (.NOT. failure) THEN
     ALLOCATE (dft_control,STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.not.failure) THEN
     dft_control%ref_count=1
     last_dft_control_id=last_dft_control_id+1
     dft_control%id_nr=last_dft_control_id
     CALL band_structure_preinit(dft_control%band,error=error)
     CALL kg_control_init(dft_control%kg_control, error=error)
     CALL qs_control_init(dft_control%qs_control, error=error)
     CALL tddfpt_control_init(dft_control%tddfpt_control, error=error)
     CALL wave_control_init(dft_control%wave_control, error=error)
     CALL xc_grid_init(dft_control%integration_grid, error=error)
     CALL xc_control_init(dft_control%xc_control, error=error)
  END IF
END SUBROUTINE dft_control_create
!***************************************************************************

!!****f* dft_types/dft_control_retain [1.0] *
!!
!!   NAME
!!     dft_control_retain
!!
!!   SYNOPSIS
!!     Subroutine dft_control_retain(dft_control, error)
!!       Type(dft_control_type), Pointer:: dft_control
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine dft_control_retain
!!
!!   FUNCTION
!!     retains the given dft_control
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - dft_control: the dft_control to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE dft_control_retain(dft_control,error)
  TYPE(dft_control_type), POINTER :: dft_control
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='dft_control_retain',&
        routineP=moduleN//':'//routineN

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(dft_control%ref_count>0,cp_failure_level,routineP,error)
     dft_control%ref_count=dft_control%ref_count+1
  END IF
END SUBROUTINE dft_control_retain
!***************************************************************************

!!****f* dft_types/dft_control_release [1.0] *
!!
!!   NAME
!!     dft_control_release
!!
!!   SYNOPSIS
!!     Subroutine dft_control_release(dft_control, error)
!!       Type(dft_control_type), Pointer:: dft_control
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine dft_control_release
!!
!!   FUNCTION
!!     releases the given dft_control
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - dft_control: the dft_control to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE dft_control_release(dft_control,error)
  TYPE(dft_control_type), POINTER :: dft_control
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='dft_control_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat

  failure=.FALSE.
  IF (ASSOCIATED(dft_control)) THEN
     CPPreconditionNoFail(dft_control%ref_count>0,cp_failure_level,routineP,error)
     dft_control%ref_count=dft_control%ref_count-1
     IF (dft_control%ref_count==0) THEN
        CALL band_structure_dealloc_ref(dft_control%band,error=error)
        CALL kg_control_dealloc_ref(dft_control%kg_control, error=error)
        CALL qs_control_dealloc_ref(dft_control%qs_control, error=error)
        CALL tddfpt_control_dealloc_ref(dft_control%tddfpt_control, error=error)
        CALL wave_control_dealloc_ref(dft_control%wave_control, error=error)
        CALL xc_grid_dealloc_ref(dft_control%integration_grid, error=error)
        CALL xc_control_dealloc_ref(dft_control%xc_control, error=error)
        DEALLOCATE(dft_control, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(dft_control)
END SUBROUTINE dft_control_release
!***************************************************************************

SUBROUTINE kg_control_init(kg_control, error)
  TYPE(kg_control_type), INTENT(out) :: kg_control
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  NULLIFY(kg_control%e_cutoff)
END SUBROUTINE kg_control_init

SUBROUTINE kg_control_dealloc_ref(kg_control, error)
  TYPE(kg_control_type), INTENT(inout) :: kg_control
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  CHARACTER(len=*), PARAMETER :: routineN = "kg_control_dealloc_ref", &
       routineP = moduleN//"/"//routineN
  INTEGER :: stat

  IF (ASSOCIATED(kg_control%e_cutoff)) THEN
     DEALLOCATE(kg_control%e_cutoff,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
END SUBROUTINE kg_control_dealloc_ref

SUBROUTINE qs_control_init(qs_control, error)
  TYPE(qs_control_type), INTENT(out) :: qs_control
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  NULLIFY(qs_control%e_cutoff)
END SUBROUTINE qs_control_init

SUBROUTINE qs_control_dealloc_ref(qs_control, error)
  TYPE(qs_control_type), INTENT(inout) :: qs_control
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  INTEGER :: stat
    CHARACTER(len=*), PARAMETER :: routineN = "qs_control_dealloc_ref", &
         routineP = moduleN//"/"//routineN

  IF (ASSOCIATED(qs_control%e_cutoff)) THEN
     DEALLOCATE(qs_control%e_cutoff,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
END SUBROUTINE qs_control_dealloc_ref

SUBROUTINE tddfpt_control_init(tddfpt_control, error)
  TYPE(tddfpt_control_type), INTENT(out) :: tddfpt_control
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  NULLIFY(tddfpt_control%lumos,tddfpt_control%lumos_eigenvalues)
END SUBROUTINE tddfpt_control_init

SUBROUTINE tddfpt_control_dealloc_ref(tddfpt_control, error)
  TYPE(tddfpt_control_type), INTENT(inout) :: tddfpt_control
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  INTEGER :: stat
  CHARACTER(len=*), PARAMETER :: routineN = "tddfpt_control_dealloc_ref", &
       routineP = moduleN//"/"//routineN
  
  IF (ASSOCIATED(tddfpt_control%lumos)) THEN
     DEALLOCATE(tddfpt_control%lumos,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
  IF (ASSOCIATED(tddfpt_control%lumos_eigenvalues)) THEN
     DEALLOCATE(tddfpt_control%lumos_eigenvalues,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
END SUBROUTINE tddfpt_control_dealloc_ref

SUBROUTINE wave_control_init(wave_control, error)
  TYPE(wave_control_type), INTENT(out) :: wave_control
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  

END SUBROUTINE wave_control_init

SUBROUTINE wave_control_dealloc_ref(wave_control, error)
  TYPE(wave_control_type), INTENT(inout) :: wave_control
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  INTEGER :: stat

END SUBROUTINE wave_control_dealloc_ref

SUBROUTINE xc_grid_init(xc_grid, error)
  TYPE(xc_grid_type), INTENT(out) :: xc_grid
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  

END SUBROUTINE xc_grid_init

SUBROUTINE xc_grid_dealloc_ref(xc_grid, error)
  TYPE(xc_grid_type), INTENT(inout) :: xc_grid
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  INTEGER :: stat

END SUBROUTINE xc_grid_dealloc_ref

SUBROUTINE xc_control_init(xc_control, error)
  TYPE(xc_control_type), INTENT(out) :: xc_control
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  

END SUBROUTINE xc_control_init

SUBROUTINE xc_control_dealloc_ref(xc_control, error)
  TYPE(xc_control_type), INTENT(inout) :: xc_control
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  INTEGER :: stat

END SUBROUTINE xc_control_dealloc_ref

END MODULE dft_types
