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
                                             gapw_control_create,&
                                             qs_control_type,&
                                             tddfpt_control_type,&
                                             wave_control_type,&
                                             xc_debug_new_routine,&
                                             xc_deriv_fd,&
                                             xc_deriv_nn10_smooth,&
                                             xc_deriv_nn50_smooth,&
                                             xc_deriv_pw,&
                                             xc_deriv_spline2,&
                                             xc_deriv_spline2_smooth,&
                                             xc_deriv_spline3,&
                                             xc_deriv_spline3_smooth,&
                                             xc_new_f_routine,&
                                             xc_old_f_routine,&
                                             xc_rho_nn10,&
                                             xc_rho_nn50,&
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
  USE cp_keywords,                     ONLY: keyword_type, &
                                             list_keywords, &
                                             find_keyword
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
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_setall,&
                                             xc_rho_cflags_type

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

! *** set up the &QS &END keyword section

PRIVATE :: init_qs_keywords
INTEGER, PARAMETER, PRIVATE :: qs_key_CUTOFF=1
INTEGER, PARAMETER, PRIVATE :: qs_key_REALSPACE_MGRIDS=2
INTEGER, PARAMETER, PRIVATE :: qs_key_COMMENSURATE=3
INTEGER, PARAMETER, PRIVATE :: qs_key_MAP_PAA=4
INTEGER, PARAMETER, PRIVATE :: qs_key_MAP_CONSISTENT=5
INTEGER, PARAMETER, PRIVATE :: qs_key_MULTIGRID=6
INTEGER, PARAMETER, PRIVATE :: qs_key_EPS_DEFAULT=7
INTEGER, PARAMETER, PRIVATE :: qs_key_EPS_CORE_CHARGE=8
INTEGER, PARAMETER, PRIVATE :: qs_key_EPS_GVG_RSPACE=9
INTEGER, PARAMETER, PRIVATE :: qs_key_EPS_PGF_AUX=10
INTEGER, PARAMETER, PRIVATE :: qs_key_EPS_PGF_ORB=11
INTEGER, PARAMETER, PRIVATE :: qs_key_EPS_PPL=12
INTEGER, PARAMETER, PRIVATE :: qs_key_EPS_PPNL=13
INTEGER, PARAMETER, PRIVATE :: qs_key_EPS_RHO=14
INTEGER, PARAMETER, PRIVATE :: qs_key_EPS_RHO_GSPACE=15
INTEGER, PARAMETER, PRIVATE :: qs_key_EPS_RHO_RSPACE=16
INTEGER, PARAMETER, PRIVATE :: qs_key_PROGRESSION_FACTOR=17
INTEGER, PARAMETER, PRIVATE :: qs_key_REL_CUTOFF=18
INTEGER, PARAMETER, PRIVATE :: qs_key_NGRID=19
INTEGER, PARAMETER, PRIVATE :: qs_key_EPSFIT=20
INTEGER, PARAMETER, PRIVATE :: qs_key_EPSISO=21
INTEGER, PARAMETER, PRIVATE :: qs_key_EPSRHO0=22
INTEGER, PARAMETER, PRIVATE :: qs_key_QUADRATURE=23
INTEGER, PARAMETER, PRIVATE :: qs_key_LMAXN1=24
INTEGER, PARAMETER, PRIVATE :: qs_key_LMAXN0=25
INTEGER, PARAMETER, PRIVATE :: qs_key_LADDN0=26
INTEGER, PARAMETER, PRIVATE :: qs_key_ALPHA0_H=27
INTEGER, PARAMETER, PRIVATE :: qs_key_ALPHA0_S=28
INTEGER, PARAMETER, PRIVATE :: qs_key_PAO=29
INTEGER, PARAMETER, PRIVATE :: qs_key_PW_GRID=30
INTEGER, PARAMETER, PRIVATE :: qs_key_RS_GRID=31
INTEGER, PARAMETER, PRIVATE :: qs_key_EXTRAPOLATION=32
INTEGER, PARAMETER, PRIVATE :: qs_key_METHOD=33
INTEGER, PARAMETER, PRIVATE :: qs_key_MULLIKEN_RESTRAINT=34
INTEGER, PARAMETER, PRIVATE :: N_qs_keywords = 34
TYPE(keyword_type), PRIVATE, DIMENSION(1:N_qs_keywords) :: qs_keywords


!!***
! *****************************************************************************

CONTAINS

  ! init the qs_keywords
  SUBROUTINE init_qs_keywords()
     qs_keywords(qs_key_CUTOFF)=keyword_type("DENSITY_CUTOFF", &
                                       "DENSITY_CUTOFF real",&
                                       "Sets the density cutoff in Ry",&
                                       "CUTOFF","ECUT","E_CUTOFF")
     qs_keywords(qs_key_REALSPACE_MGRIDS)=keyword_type("REALSPACE_MGRIDS", &
                                       "REALSPACE_MGRIDS",&
                                       "Use realspace techniques for multigrids",&
                                       "REALSPACE_MGRID","","")
     qs_keywords(qs_key_COMMENSURATE)=keyword_type("COMMENSURATE", &
                                       "COMMENSURATE",&
                                       "Have realspace multigrids grids that are commensurate in real space",&
                                       "","","")
     qs_keywords(qs_key_MAP_PAA)=keyword_type("MAP_PAA", &
                                       "MAP_PAA",&
                                       "Put diagonal density matrix elements on the highest cutoff multigrid",&
                                       "","","")
     qs_keywords(qs_key_MAP_CONSISTENT)=keyword_type("MAP_CONSISTENT", &
                                       "MAP_CONSISTENT",&
                                       "Guarantee an exact derivative (Hks) of the energy",&
                                       "","","")
     qs_keywords(qs_key_MULTIGRID)=keyword_type("MULTIGRID", &
                                       "MULTIGRID integer real ... real",&
                                       "Set the cutoffs of the multigrids",&
                                       "MULTI_GRID","MGRID","")
     qs_keywords(qs_key_EPS_DEFAULT)=keyword_type("EPS_DEFAULT", &
                                       "EPS_DEFAULT real",&
                                       "Try setting all EPS_xxx to values leading to an energy correct up to EPS_DEFAULT",&
                                       "","","")
     qs_keywords(qs_key_EPS_CORE_CHARGE)=keyword_type("EPS_CORE_CHARGE", &
                                       "EPS_CORE_CHARGE real",&
                                       "Precision for mapping the core charges",&
                                       "","","")
     qs_keywords(qs_key_EPS_GVG_RSPACE)=keyword_type("EPS_GVG_RSPACE", &
                                       "EPS_GVG_RSPACE real",&
                                       "Sets precision of the realspace KS matrix element integration",&
                                       "EPS_GVG","","")
     qs_keywords(qs_key_EPS_PGF_AUX)=keyword_type("EPS_PGF_AUX", &
                                       "EPS_PGF_AUX real",&
                                       "XXXXXXXXXXX set precision of the ...",&
                                       "","","")
     qs_keywords(qs_key_EPS_PGF_ORB)=keyword_type("EPS_PGF_ORB", &
                                       "EPS_PGF_ORB real",&
                                       "Sets precision of the overlap matrix elements",&
                                       "","","")
     qs_keywords(qs_key_EPS_PPL)=keyword_type("EPS_PPL", &
                                       "EPS_PPL real",&
                                       "Sets precision of the local part of the pseudo potential",&
                                       "","","")
     qs_keywords(qs_key_EPS_PPNL)=keyword_type("EPS_PPNL", &
                                       "EPS_PPNL real",&
                                       "Sets precision of the non-local part of the pseudo potential",&
                                       "","","")
     qs_keywords(qs_key_EPS_RHO)=keyword_type("EPS_RHO", &
                                       "EPS_RHO real",&
                                       "Sets precision of the density mapping on the grids",&
                                       "","","")
     qs_keywords(qs_key_EPS_RHO_GSPACE)=keyword_type("EPS_RHO_GSPACE", &
                                       "EPS_RHO_GSPACE real",&
                                       "Sets precision of the density mapping in gspace. Unused ?",&
                                       "","","")
     qs_keywords(qs_key_EPS_RHO_RSPACE)=keyword_type("EPS_RHO_RSPACE", &
                                       "EPS_RHO_RSPACE real",&
                                       "Sets precision of the density mapping in rspace",&
                                       "","","")
     qs_keywords(qs_key_PROGRESSION_FACTOR)=keyword_type("PROGRESSION_FACTOR", &
                                       "PROGRESSION_FACTOR real",&
                                       "Subsequent multigrids have a cutoff that is profac smaller then the next bigger grid",&
                                       "RROFAC","","")
     qs_keywords(qs_key_REL_CUTOFF)=keyword_type("REL_CUTOFF", &
                                       "REL_CUTOFF real",&
                                       "Determines the grid at which a Gaussian is mapped,"// &
                                                " giving the cutoff used for a gaussian with alpha=1",&
                                       "RELATIVE_CUTOFF","","")
     qs_keywords(qs_key_NGRID)=keyword_type("NGRID", &
                                       "NGRID integer",&
                                       "The number of multigrids that should be used",&
                                       "NGRID_LEVEL","","")
     qs_keywords(qs_key_EPSFIT)=keyword_type("EPSFIT", &
                                       "EPSFIT real",&
                                       "GAPW : precision to give the extention of a hard gaussian ",&
                                       "EPS_FIT","","")
     qs_keywords(qs_key_EPSISO)=keyword_type("EPSISO", &
                                       "EPSISO real",&
                                       "GAPW : precision to determine an isolated projector",&
                                       "EPS_ISO","","")
     qs_keywords(qs_key_EPSRHO0)=keyword_type("EPSRHO0", &
                                       "EPSRHO0 real",&
                                       "GAPW : precision to determine the range of V(rho0-rho0soft)",&
                                       "EPSVRHO0","EPS_VRHO0","")
     qs_keywords(qs_key_QUADRATURE)=keyword_type("QUADRATURE", &
                                       "QUADRATURE {GC_SIMPLE|GC_TRANSFORMED|GC_LOG}",&
                                       "GAPW: algorithm to construct the atomic radial grids",&
                                       "","","")
     qs_keywords(qs_key_LMAXN1)=keyword_type("LMAXN1", &
                                       "LMAXN1 integer",&
                                       "GAPW : max L number for espansion of the atomic densities in spherical gaussians",&
                                       "LMAXRHO1","","")
     qs_keywords(qs_key_LMAXN0)=keyword_type("LMAXN0", &
                                       "LMAXN0 integer",&
                                       "GAPW : max L number for the expansion compensation densities in spherical gaussians",&
                                       "LMAXRHO0","","")
     qs_keywords(qs_key_LADDN0)=keyword_type("LADDN0", &
                                       "LADDN0 integer",&
                                       "GAPW : XXXXXXXXXX",&
                                       "","","")
     qs_keywords(qs_key_ALPHA0_H)=keyword_type("ALPHA0_HARD", &
                                       "ALPHA0_HARD real ",&
                                       "GAPW : Exponent for hard compensation charge",&
                                       "ALPHA0_H","","")
     qs_keywords(qs_key_ALPHA0_S)=keyword_type("ALPHA0_SOFT", &
                                       "ALPHA0_SOFT real ",&
                                       "GAPW : Exponent for soft compensation charge",&
                                       "ALPHA0_S","","")
     qs_keywords(qs_key_PAO)=keyword_type("PAO", &
                                       "PAO",&
                                       "Use the PAO method",&
                                       "","","")
     qs_keywords(qs_key_PW_GRID)=keyword_type("PW_GRID", &
                                       "PW_GRID {SPHERICAL|NON-SPHERICAL|FULLSPACE|HALFSPACE}",&
                                       "What kind of PW_GRID should be employed",&
                                       "","","")
     qs_keywords(qs_key_RS_GRID)=keyword_type("RS_GRID", &
                                       "RS_GRID {DISTRIBUTED|REPLICATED}",&
                                       "Parallellisation strategy of the realspace grids",&
                                       "","","")
     qs_keywords(qs_key_EXTRAPOLATION)=keyword_type("EXTRAPOLATION", &
                                       "EXTRAPOLATION {PS N|USE_GUESS|USE_PREV_P|USE_PREV_WF|"//&
                                                 "USE_PREV_RHO_R|LINEAR_WF|LINEAR_P|LINEAR_PS}",&
                                       "Extrapolation strategy for e.g. MD",&
                                       "INTERPOLATION","WF_INTERPOLATION","")
     qs_keywords(qs_key_METHOD)=keyword_type("METHOD", &
                                       "METHOD {GAPW|GPW|KG_POL|KG_NOPOL|KG_GPW}",&
                                       "What method should be employed",&
                                       "","","")
     qs_keywords(qs_key_MULLIKEN_RESTRAINT)=keyword_type("MULLIKEN_RESTRAINT", &
                                       "MULLIKEN_RESTRAINT",&
                                       "Employ a restraint on the mulliken charges",&
                                       "","","")
  END SUBROUTINE

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
!!  tau_cutoff       tau_cut                                                  !
!!  spin             nspins                                                   !
!!  lsd                                                                       !
!!  grid             plane_wave     [cutoff ecut, mesh nx ny nz]              !
!!  charge           charge                                                   !
!!  xc_smooth_rho    [ NONE, NN10, SPLINE2, SPLINE3 ](in order of smoothness) !
!!  xc_deriv         [ PW, SPLINE3, SPLINE2, NN10_SMOOTH, SPLINE2_SMOOTH (NN6)!
!!                     SPLINE3_SMOOTH (NN4) ] (in order of growing smoothness)!
!!  functional_routine [ old, new, test ]                                     !
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
    LOGICAL                                  :: failure, section_found, &
                                                xc_fun_set
    TYPE(cp_error_type)                      :: sub_error
    TYPE(xc_rho_cflags_type)                 :: needs

!   ---------------------------------------------------------------------------

    failure=.FALSE.
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
    dft_control % xc_control % density_smooth_cut_range = 0.0_dp
    dft_control % xc_control % gradient_cut = 1.0E-10_dp
    dft_control % xc_control % tau_cut = 1.0e-10_dp
    dft_control % xc_control % gradient_functionals = .FALSE.
    dft_control % xc_control % crossterms = .FALSE.
    dft_control % xc_control % functional_routine = xc_new_f_routine

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
          CASE ("TAU_CUTOFF")
            CALL read_object(dft_control % xc_control % tau_cut)
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
             CASE ("SPLINE2","NN6")
                dft_control % xc_rho_smooth_id= xc_rho_spline2_smooth
             CASE ("SPLINE3","NN4")
                dft_control % xc_rho_smooth_id= xc_rho_spline3_smooth
             CASE ("NN10")
                dft_control % xc_rho_smooth_id= xc_rho_nn10
             CASE ("NN50")
                dft_control % xc_rho_smooth_id= xc_rho_nn50
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
             CASE ("FD")
                dft_control % xc_deriv_method_id = xc_deriv_fd
             CASE ("NN10_SMOOTH")
                dft_control % xc_deriv_method_id = xc_deriv_nn10_smooth
             CASE ("NN50_SMOOTH")
                dft_control % xc_deriv_method_id = xc_deriv_nn50_smooth
             CASE ("SPLINE2_SMOOTH","NN6_SMOOTH")
                dft_control % xc_deriv_method_id = xc_deriv_spline2_smooth
             CASE ("SPLINE3_SMOOTH","NN4_SMOOTH")
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
             CASE ("DEBUG","DEBUG_NEW")
                dft_control % xc_control % functional_routine = &
                     xc_debug_new_routine                
             CASE DEFAULT
                CALL stop_parser(routine=routine,&
                  message="INVALID_KEYWORD for FUNCTIONAL_ROUTINE "//str1)
             END SELECT
          CASE DEFAULT
            CALL stop_parser(routine=routine,&
                 message="INVALID_KEYWORD "//keyword)
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
       dft_control % xc_control % functionals ( 1 ) = "XALPHA"
       dft_control % xc_control % functionals ( 2 ) = "PZ81_DMC"
    END IF    
    
    IF (.NOT.failure) THEN
       ! assess the needs of the requested functionals
       CALL xc_rho_cflags_setall(needs,.FALSE.,error=error)
       CALL cp_error_init(sub_error,stop_level=cp_fatal_level, &
            print_level=cp_fatal_level,template_error=error)
       DO i=1,SIZE(dft_control % xc_control % functionals)
          IF ( dft_control % xc_control % functionals(i) /= "NONE") THEN
             CALL xc_functional_get_info(&
                  functional=dft_control % xc_control % functionals(i),&
                  lsd=(dft_control % nspins /= 1),&
                  needs=needs, error=sub_error)
             CALL cp_error_check(sub_error,failure)
             
             IF (failure) THEN
                failure=.FALSE.
                CALL cp_error_reset(sub_error)
             END IF
          END IF
       END DO
       CALL cp_error_dealloc_ref(sub_error)
       
       ! sets the flag that the kinetic energy needs to be calculated
       ! if some functional needs it
       dft_control % use_kinetic_energy_density=&
            dft_control % use_kinetic_energy_density .OR.&
            needs%tau .OR. needs%tau_spin
       ! fills in the old structure (I would like to remove it)
       dft_control % xc_control % gradient_functionals=&
            dft_control % xc_control % gradient_functionals.OR.&
            needs%norm_drho .OR. needs%norm_drho_spin .OR.&
            needs%drho .OR. needs%drho_spin
       ! presence of crossterms cannot be determined before evaluating
       ! the functional
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

    TYPE(qs_control_type), INTENT(INOUT)     :: qs_control
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE read_qs_control (MODULE dft_types)"

    CHARACTER(LEN=3)                         :: test_result
    CHARACTER(LEN=60)                        :: keyword, string, string2
    INTEGER                                  :: igrid_level, istat, &
                                                ivalue, ngrid_level
    LOGICAL                                  :: cutoff_set, multgrid_set, &
                                                ngrid_set, progression_set, &
                                                section_found
    REAL(KIND=dp)                            :: value

!   ---------------------------------------------------------------------------
!   *** Load the default values ***
! is now set as would result from EPS_DEFAULT=1.0E-10

    qs_control%method = "GPW"
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
    qs_control%semi_empirical = .FALSE.
    qs_control%pw_grid_opt%spherical = .FALSE.
    qs_control%pw_grid_opt%fullspace = .TRUE.
    qs_control%rs_distributed = .FALSE.
    qs_control%wf_interpolation_method_nr = wfi_linear_p_method_nr
    qs_control%commensurate_mgrids = .FALSE.
    qs_control%realspace_mgrids = .FALSE.
    qs_control%map_paa = .FALSE.
    qs_control%map_consistent = .FALSE.
    qs_control%polarization = .FALSE.
    qs_control%mulliken_restraint = .FALSE.
    qs_control%wf_extrapolation_order = 0

    CALL gapw_control_create(qs_control%gapw_control)

    qs_control%gapw_control%eps_fit =  1.0E-4_dp
    qs_control%gapw_control%eps_iso =  1.0E-12_dp
    qs_control%gapw_control%eps_Vrho0  = 1.0E-6_dp
    qs_control%gapw_control%ladd_rho0 = 99
    qs_control%gapw_control%lmax_rho0 = 2
    qs_control%gapw_control%lmax_sphere = 6 
    qs_control%gapw_control%quadrature = 3
    qs_control%gapw_control%quadrature_name = "LogTrans-Gauss-Chebyshev"
    qs_control%gapw_control%lrho1_eq_lrho0 = .FALSE.
    qs_control%gapw_control%alpha0_hard_from_input = .FALSE.
    qs_control%gapw_control%alpha0_soft_from_input = .FALSE.
    qs_control%gapw_control%alpha0_hard = 0.0_dp
    qs_control%gapw_control%alpha0_soft = 0.0_dp

    qs_control%mulliken_restraint_control%strength=0.1_dp
    qs_control%mulliken_restraint_control%target=1.0_dp

!
    cutoff_set = .FALSE.
    multgrid_set = .FALSE.
    ngrid_set = .FALSE.
    progression_set = .FALSE.

!   *** Read the input section ***
    CALL init_qs_keywords()
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
          SELECT CASE (find_keyword(qs_keywords,keyword))
          CASE (qs_key_CUTOFF)
            CALL read_object(value)
            qs_control%cutoff = MAX(0.01_dp,0.5_dp*value) 
                                        ! ok input is in Ry, value is in a.u.
            cutoff_set = .TRUE.
          CASE (qs_key_REALSPACE_MGRIDS)
            qs_control%realspace_mgrids = .TRUE.
          CASE (qs_key_COMMENSURATE)
            qs_control%commensurate_mgrids = .TRUE.
          CASE (qs_key_MAP_PAA)
            qs_control%map_paa = .TRUE.
          CASE (qs_key_MAP_CONSISTENT)
            qs_control%map_consistent = .TRUE.
          CASE (qs_key_MULTIGRID)
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
          CASE (qs_key_EPS_DEFAULT)
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
          CASE (qs_key_EPS_CORE_CHARGE)
            CALL read_object(qs_control%eps_core_charge)
          CASE (qs_key_EPS_GVG_RSPACE)
            CALL read_object(qs_control%eps_gvg_rspace)
          CASE (qs_key_EPS_PGF_AUX)
            CALL read_object(qs_control%eps_pgf_aux)
          CASE (qs_key_EPS_PGF_ORB)
            CALL read_object(qs_control%eps_pgf_orb)
          CASE (qs_key_EPS_PPL)
            CALL read_object(qs_control%eps_ppl)
          CASE (qs_key_EPS_PPNL)
            CALL read_object(qs_control%eps_ppnl)
          CASE (qs_key_EPS_RHO)
            CALL read_object(value)
            qs_control%eps_rho_gspace = value
            qs_control%eps_rho_rspace = value
          CASE (qs_key_EPS_RHO_GSPACE)
            CALL read_object(qs_control%eps_rho_gspace)
          CASE (qs_key_EPS_RHO_RSPACE)
            CALL read_object(qs_control%eps_rho_rspace)
          CASE (qs_key_PROGRESSION_FACTOR)
            CALL read_object(value)
            progression_set=.TRUE.
            qs_control%progression_factor = MAX(0.01_dp,value)
          CASE (qs_key_REL_CUTOFF)
            CALL read_object(value)
            qs_control%relative_cutoff = MAX(5.0_dp,0.5*value) 
                                        ! now input is in Ry, value is in a.u.
          CASE (qs_key_NGRID)
            CALL read_object(ngrid_level)
            ngrid_set = .TRUE.
          CASE (qs_key_METHOD)
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
            CASE ("EHT")
              qs_control%method = "EHT"
              qs_control%semi_empirical = .TRUE.
              qs_control%se_control%orthogonal_basis = .FALSE.
            CASE ("MNDO")
              qs_control%method = "MNDO"
              qs_control%semi_empirical = .TRUE.
              qs_control%se_control%orthogonal_basis = .TRUE.
            CASE ("AM1")
              qs_control%method = "AM1"
              qs_control%semi_empirical = .TRUE.
              qs_control%se_control%orthogonal_basis = .TRUE.
            CASE ("PM3")
              qs_control%method = "PM3"
              qs_control%semi_empirical = .TRUE.
              qs_control%se_control%orthogonal_basis = .TRUE.
            CASE DEFAULT
              CALL stop_parser(routine,"Invalid method <"//TRIM(keyword)//&
                                       "> specified for METHOD")
            END SELECT
          CASE (qs_key_EPSFIT)
            CALL read_object(value)
            qs_control%gapw_control%eps_fit = value
          CASE (qs_key_EPSISO)
            CALL read_object(value)
            qs_control%gapw_control%eps_iso = value
          CASE (qs_key_EPSRHO0)
            CALL read_object(value)
            qs_control%gapw_control%eps_Vrho0 = value
          CASE (qs_key_QUADRATURE)
            CALL read_object(string,lower_to_upper=.TRUE.)
            SELECT CASE (string)
            CASE ("GC_SIMPLE")
              qs_control%gapw_control%quadrature=1
              qs_control%gapw_control%quadrature_name="Gauss-Chebyshev"
            CASE ("GC_TRANSFORMED")
              qs_control%gapw_control%quadrature=2
              qs_control%gapw_control%quadrature_name="Trans_Gauss-Chebyshev"
            CASE ("GC_LOG")
              qs_control%gapw_control%quadrature=3
              qs_control%gapw_control%quadrature_name="LogTrans-Gauss-Chebyshev"
            CASE DEFAULT
              CALL stop_parser(routine,"quadrature definition")
            END SELECT
          CASE (qs_key_LMAXN1)
            CALL read_object(ivalue)
            qs_control%gapw_control%lmax_sphere = ivalue
          CASE (qs_key_LMAXN0)
            CALL read_object(ivalue)
            qs_control%gapw_control%lmax_rho0 = ivalue
          CASE (qs_key_LADDN0)
            CALL read_object(ivalue)
            qs_control%gapw_control%ladd_rho0 = ivalue
          CASE (qs_key_ALPHA0_H)
            CALL read_object(value)
            qs_control%gapw_control%alpha0_hard = value
            qs_control%gapw_control%alpha0_hard_from_input = .TRUE.
          CASE (qs_key_ALPHA0_S)
            CALL read_object(value)
            qs_control%gapw_control%alpha0_soft = value
            qs_control%gapw_control%alpha0_soft_from_input = .TRUE.
          CASE (qs_key_PAO)
            qs_control%pao = .TRUE.
          CASE (qs_key_MULLIKEN_RESTRAINT)
            qs_control%mulliken_restraint = .TRUE.
            SELECT CASE (test_object())
            CASE ("STR")
                CALL read_object(string,lower_to_upper=.TRUE.)
                SELECT CASE (string)
                CASE ("STRENGTH")
                   CALL read_object(qs_control%mulliken_restraint_control%strength, &
                                    lower_to_upper=.TRUE.)
                CASE ("TARGET")
                   CALL read_object(qs_control%mulliken_restraint_control%target, &
                                    lower_to_upper=.TRUE.)
                CASE DEFAULT
                  CALL stop_parser(routine,"unexpected string (" &
                                   //TRIM(string)//") following MULLIKEN_RESTRAINT")
                END SELECT
            CASE ("EOL")
                ! do nothing
            CASE DEFAULT
              CALL stop_parser(routine,"unexpected data following MULLIKEN_RESTRAINT")
            END SELECT
          CASE (qs_key_PW_GRID)
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
          CASE (qs_key_RS_GRID)
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
         CASE (qs_key_EXTRAPOLATION)
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
            CALL list_keywords(qs_keywords,globenv%scr)
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
        IF ( qs_control%semi_empirical ) THEN
          qs_control%cutoff = 5._dp ! to be adjusted later
        ELSE
          qs_control%cutoff = MAX(0.01_dp,0.5_dp*280._dp) ! give a 280 Ry default cutoff. 
                                                          ! should be very fine
        END IF
      ELSE
        qs_control%cutoff = qs_control%e_cutoff(1)
      END IF
    END IF
    IF ( .NOT. ngrid_set ) THEN
      IF ( qs_control%semi_empirical ) THEN
        ngrid_level = 1
      ELSE
        ngrid_level = 4
      END IF
    END IF
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
    ! do not print dft section for semi empirical methods
    IF(qs_control%semi_empirical) THEN
      globenv%print%dft_control_parameters=.FALSE.
    END IF

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

       CASE ("INVERT_S", "INV_S", "S-1")
          CALL read_object(field, lower_to_upper=.TRUE.)
          SELECT CASE (TRIM(field))
          CASE ("Y", "YES") 
             t_control%invert_S = .TRUE.
          CASE ("N", "NO")
             t_control%invert_S = .FALSE.
          CASE DEFAULT
             t_control%invert_S = .FALSE.
          END SELECT

       CASE ("RES_ETYPE", &
             "RESTRICTED_EXCITATIONS_TYPE")
          CALL read_object(field, lower_to_upper=.TRUE.)
          SELECT CASE (TRIM(field))
          CASE ("SINGLET", "SINGLETS", "S")
             t_control%res_etype = "SINGLET"
          CASE ("TRIPLET", "TRIPLETS", "T")
             t_control%res_etype = "TRIPLET"
          CASE DEFAULT
             t_control%res_etype = "SINGLET"
          END SELECT

       CASE ("DIAG_METHOD")
          CALL read_object(field, lower_to_upper=.TRUE.)
          SELECT CASE (TRIM(field))
          CASE ("DAVIDSON")
             t_control%diag_method = "DAVIDSON"
          CASE ("LANCZOS")
             t_control%diag_method = "LANCZOS"
          CASE DEFAULT
             t_control%diag_method = "DAVIDSON"
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
            "DFT|          tau     ",dft_control%xc_control%tau_cut,&
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
          WRITE (output_unit, fmt="(' FUNCTIONAL| ROUTINE=',a)")&
               TRIM(xc_get_routine_label(dft_control % xc_control&
               % functional_routine))
          CALL cp_error_init(sub_error,stop_level=cp_fatal_level, &
               print_level=cp_fatal_level,template_error=error)
          DO i=1,SIZE(dft_control % xc_control % functionals)
             IF ( dft_control % xc_control % functionals(i) /= "NONE") THEN
                CALL xc_functional_get_info(&
                     functional=dft_control % xc_control % functionals(i),&
                     lsd=(dft_control % nspins /= 1),&
                     reference=reference, error=sub_error)
                CALL cp_error_check(sub_error,failure)

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
        IF (qs_control%commensurate_mgrids) &
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
      IF (qs_control%gapw) THEN
         WRITE (UNIT=output_unit,FMT="(T2,A,T73,ES8.1)")&
           "QS| GAPW|                   eps_fit:",&
           qs_control%gapw_control%eps_fit,&
           "QS| GAPW|                   eps_iso:",&
           qs_control%gapw_control%eps_iso
         WRITE (UNIT=output_unit,FMT="(T2,A,T55,A30)")&
           "QS| GAPW|   atom-r-grid: quadrature:",&
           qs_control%gapw_control%quadrature_name
         WRITE (UNIT=output_unit,FMT="(T2,A,T71,I10)")&
           "QS| GAPW|      atom-s-grid:  max l :",&
           qs_control%gapw_control%lmax_sphere ,& 
           "QS| GAPW|      max_l_rho0 :",&
            qs_control%gapw_control%lmax_rho0
      END IF
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
  CASE (xc_rho_nn50)
     res="NN50"
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
  CASE (xc_deriv_nn10_smooth)
     res="NN10_SMOOTH"
  CASE (xc_deriv_nn50_smooth)
     res="NN50_SMOOTH"
  CASE (xc_deriv_fd)
     res="FD"
  CASE default
     WRITE (res,"('UKN',i6)") xc_deriv_method_id
  END SELECT
END FUNCTION xc_get_deriv_method_label
!***************************************************************************

!!****f* cp_control_types/xc_get_routine_label [1.0] *
!!
!!   NAME
!!     xc_get_routine_label
!!
!!   FUNCTION
!!     returns a string that describes the functional routine used in the xc 
!!     calculation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - xc_functional_routine_id: the id that represent the functional 
!!       routine
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
FUNCTION xc_get_routine_label(xc_functional_routine_id,error) RESULT(res)
    INTEGER, INTENT(in)                      :: xc_functional_routine_id
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    CHARACTER(len=20)                        :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_get_routine_label', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  SELECT CASE(xc_functional_routine_id)
  CASE(xc_old_f_routine)
     res="OLD"
  CASE(xc_new_f_routine)
     res="NEW"
  CASE(xc_test_f_routine)
     res="TEST (new vs old)"
  CASE(xc_test_lsd_f_routine)
     res="TEST_LSD"
  CASE(xc_debug_new_routine)
     res="DEBUG_NEW"
  CASE default
     WRITE (res,"('UKN',i6)") xc_functional_routine_id
  END SELECT
END FUNCTION xc_get_routine_label
!***************************************************************************

END MODULE cp_control_utils
