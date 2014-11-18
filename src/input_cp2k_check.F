!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief checks the input and perform some automatic "magic" on it
!> \par History
!>      01.2006 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE input_cp2k_check
  USE cp2k_info,                       ONLY: id_cp2k_version
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_set_create,&
                                             cp_unit_set_release,&
                                             cp_unit_set_type
  USE input_constants,                 ONLY: &
       do_region_global, do_thermo_al, do_thermo_csvr, do_thermo_gle, &
       do_thermo_nose, do_thermo_same_as_part, id_release_version, &
       npt_f_ensemble, npt_i_ensemble, vdw_nl_LMKLL, xc_funct_b3lyp, &
       xc_funct_beefvdw, xc_funct_blyp, xc_funct_bp, xc_funct_hcth120, &
       xc_funct_no_shortcut, xc_funct_olyp, xc_funct_pade, xc_funct_pbe, &
       xc_funct_pbe0, xc_funct_tpss, xc_funct_xwpbe, xc_none, &
       xc_vdw_fun_nonloc
  USE input_keyword_types,             ONLY: keyword_type
  USE input_parsing,                   ONLY: section_vals_parse
  USE input_section_types,             ONLY: &
       section_type, section_vals_check_release, section_vals_create, &
       section_vals_get, section_vals_get_subs_vals, &
       section_vals_get_subs_vals3, section_vals_release, &
       section_vals_remove_values, section_vals_set_subs_vals, &
       section_vals_type, section_vals_val_get, section_vals_val_set, &
       section_vals_val_unset
  USE input_val_types,                 ONLY: logical_t
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_input_constants,              ONLY: do_vwn5
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_check'

  PUBLIC :: check_cp2k_input, xc_functionals_expand, remove_restart_info

CONTAINS

! *****************************************************************************
!> \brief performs further checks on an input that parsed successfully
!> \param input_declaration ...
!> \param input_file the parsed input
!> \param para_env ...
!> \param output_unit ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
!> \note
!>      at the moment does nothing
! *****************************************************************************
  SUBROUTINE check_cp2k_input(input_declaration,input_file,para_env,output_unit,error)
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(section_vals_type), POINTER         :: input_file
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN), OPTIONAL            :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'check_cp2k_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iforce_eval, &
                                                nforce_eval
    LOGICAL                                  :: enable_unsupported_features, &
                                                explicit, explicit_mix, &
                                                failure
    TYPE(section_vals_type), POINTER         :: section, section1, section2, &
                                                section3, sections

    CALL timeset(routineN,handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(input_file),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(input_file%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       ! ext_restart
       IF(PRESENT(output_unit)) &
          CALL handle_ext_restart(input_declaration, input_file,para_env,output_unit,error)

       ! checks on force_eval section
       sections => section_vals_get_subs_vals(input_file,"FORCE_EVAL",error=error)
       CALL section_vals_get(sections, n_repetition=nforce_eval, error=error)

       ! multiple force_eval only if present RESPA or MIXED calculation is performed
       section2 => section_vals_get_subs_vals(input_file,"MOTION%MD%RESPA",error=error)
       CALL section_vals_get(section2,explicit=explicit,error=error)
       DO iforce_eval=1,nforce_eval
          section3 =>section_vals_get_subs_vals(sections,"MIXED",&
               i_rep_section=iforce_eval,error=error)
          CALL section_vals_get(section3,explicit=explicit_mix,error=error)
          IF(explicit_mix)EXIT
       END DO

       IF ((explicit.AND.(nforce_eval==1)).OR.(.NOT.explicit.AND.(nforce_eval>1)))THEN
          IF((explicit_mix.AND.(nforce_eval==1)).OR.(.NOT.explicit_mix.AND.(nforce_eval>1)))THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP,"Error multiple force_env without RESPA or MIXED, or RESPA with one single "//&
                  " or MIXED with only two force_env section.",error,failure)
          END IF
       END IF
       DO iforce_eval = 1,nforce_eval
          section => section_vals_get_subs_vals3(sections,"DFT",i_rep_section=iforce_eval,error=error)
          ! xc: expand and fix default for tddfpt
          section1 => section_vals_get_subs_vals(section,"XC",error=error)
          section2 => section_vals_get_subs_vals(section,"XC%XC_FUNCTIONAL",error=error)
          CALL xc_functionals_expand(section2,section1,error=error)
          section1 => section_vals_get_subs_vals(section,"TDDFPT%XC",error=error)
          section2 => section_vals_get_subs_vals(section,"TDDFPT%XC%XC_FUNCTIONAL",error=error)
          CALL section_vals_get(section2,explicit=explicit,error=error)
          IF (explicit) THEN
             CALL xc_functionals_expand(section2,section1,error=error)
          ELSE
             section2 => section_vals_get_subs_vals(section,"XC%XC_FUNCTIONAL",error=error)
             CALL section_vals_set_subs_vals(section,"TDDFPT%XC%XC_FUNCTIONAL",section2,error=error)
          END IF
       END DO

       ! Check for unsupported features for release version
       IF (id_cp2k_version==id_release_version) THEN
          CALL section_vals_val_get(input_file,"GLOBAL%ENABLE_UNSUPPORTED_FEATURES",&
               l_val=enable_unsupported_features, error=error)
          CALL section_vals_check_release(input_file,enable_unsupported_features,error)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE check_cp2k_input

! *****************************************************************************
!> \brief expand a shortcutted functional section
!> \param functionals the functional section to expand
!> \param xc_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE xc_functionals_expand(functionals,xc_section,error)
    TYPE(section_vals_type), POINTER         :: functionals, xc_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_functionals_expand', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: shortcut
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (.NOT. failure) THEN
       CALL section_vals_val_get(functionals,"_SECTION_PARAMETERS_",&
            i_val=shortcut,error=error)
       SELECT CASE(shortcut)
       CASE(xc_funct_no_shortcut, xc_none)
          ! nothing to expand
       CASE(xc_funct_pbe0)
          CALL section_vals_val_set(functionals,"PBE%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"PBE%SCALE_X",&
               r_val=0.75_dp,error=error)
          CALL section_vals_val_set(functionals,"PBE%SCALE_C",&
               r_val=1.0_dp,error=error)
          ! Hartree Fock Exact Exchange
          CALL section_vals_val_set(xc_section,"HF%FRACTION",&
               r_val=0.25_dp,error=error)
          CALL section_vals_val_set(functionals,"_SECTION_PARAMETERS_",&
               i_val=xc_funct_no_shortcut,error=error)
       CASE(xc_funct_beefvdw)       
          CALL section_vals_val_set(functionals,"PBE%_SECTION_PARAMETERS_",& !40% PBEc
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"PBE%SCALE_C",&
               r_val=0.3998335231_dp,error=error)
          CALL section_vals_val_set(functionals,"PBE%SCALE_X",&  !no PBEx
               r_val=0.0000000000_dp,error=error)

               !PW92 correlation functional from libxc is required.
               !The cp2k-native PW92 gives disagreeing results (in the 0.01E_H
               !decimal) and yields inconsistent forces in a DEBUG run.
               !(rk, 6.3.2014)
          CALL section_vals_val_set(functionals,"LIBXC%_SECTION_PARAMETERS_",& !60%LDA
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"LIBXC%FUNCTIONAL",&
               c_val="XC_LDA_C_PW",error=error)
          CALL section_vals_val_set(functionals,"LIBXC%SCALE",&
               r_val=0.6001664769_dp,error=error)

          CALL section_vals_val_set(functionals,"BEEF%_SECTION_PARAMETERS_",& !BEEF exchange
               l_val=.TRUE.,error=error)

               !NONLOCAL, LMKLL.
          CALL section_vals_val_set(xc_section,"VDW_POTENTIAL%DISPERSION_FUNCTIONAL",&
               i_val=xc_vdw_fun_nonloc,error=error)
          CALL section_vals_val_set(xc_section,"VDW_POTENTIAL%NON_LOCAL%TYPE",&
               i_val=vdw_nl_LMKLL,error=error)
          CALL section_vals_val_set(functionals,"_SECTION_PARAMETERS_",&
               i_val=xc_funct_no_shortcut,error=error)
       CASE(xc_funct_b3lyp)
          CALL section_vals_val_set(functionals,"BECKE88%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"BECKE88%SCALE_X",&
               r_val=0.72_dp,error=error)
          CALL section_vals_val_set(functionals,"LYP%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"LYP%SCALE_C",&
               r_val=0.81_dp,error=error)
          CALL section_vals_val_set(functionals,"VWN%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"VWN%FUNCTIONAL_TYPE",&
               i_val=do_vwn5,error=error)
          CALL section_vals_val_set(functionals,"VWN%SCALE_C",&
               r_val=0.19_dp,error=error)
          CALL section_vals_val_set(functionals,"XALPHA%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"XALPHA%SCALE_X",&
               r_val=0.08_dp,error=error)
          ! Hartree Fock Exact Exchange
          CALL section_vals_val_set(xc_section,"HF%FRACTION",&
               r_val=0.20_dp,error=error)
          CALL section_vals_val_set(functionals,"_SECTION_PARAMETERS_",&
               i_val=xc_funct_no_shortcut,error=error)
       CASE(xc_funct_blyp)
          CALL section_vals_val_set(functionals,"BECKE88%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"LYP%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"_SECTION_PARAMETERS_",&
               i_val=xc_funct_no_shortcut,error=error)
       CASE(xc_funct_bp)
          CALL section_vals_val_set(functionals,"BECKE88%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"P86C%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"_SECTION_PARAMETERS_",&
               i_val=xc_funct_no_shortcut,error=error)
       CASE(xc_funct_pade)
          CALL section_vals_val_set(functionals,"PADE%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"_SECTION_PARAMETERS_",&
               i_val=xc_funct_no_shortcut,error=error)
       CASE(xc_funct_pbe)
          CALL section_vals_val_set(functionals,"PBE%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"_SECTION_PARAMETERS_",&
               i_val=xc_funct_no_shortcut,error=error)
       CASE(xc_funct_xwpbe)
          CALL section_vals_val_set(functionals,"XWPBE%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"_SECTION_PARAMETERS_",&
               i_val=xc_funct_no_shortcut,error=error)
       CASE(xc_funct_tpss)
          CALL section_vals_val_set(functionals,"TPSS%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"_SECTION_PARAMETERS_",&
               i_val=xc_funct_no_shortcut,error=error)
       CASE(xc_funct_olyp)
          CALL section_vals_val_set(functionals,"OPTX%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"LYP%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"_SECTION_PARAMETERS_",&
               i_val=xc_funct_no_shortcut,error=error)
       CASE(xc_funct_hcth120)
          CALL section_vals_val_set(functionals,"HCTH%_SECTION_PARAMETERS_",&
               l_val=.TRUE.,error=error)
          CALL section_vals_val_set(functionals,"HCTH%PARAMETER_SET",&
               i_val=120,error=error)
          CALL section_vals_val_set(functionals,"_SECTION_PARAMETERS_",&
               i_val=xc_funct_no_shortcut,error=error)
       CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown shortcut "//TRIM(ADJUSTL(cp_to_string(shortcut))),&
               error,failure)
       END SELECT
    END IF
  END SUBROUTINE xc_functionals_expand

! *****************************************************************************
!> \brief Replaces the requested sections in the input with those found
!>      in the external restart (EXT_RESTART%RESTART_FILE_NAME).
!> \param input_declaration ...
!> \param input_file the input file to initialize
!> \param para_env ...
!> \param output_unit ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE handle_ext_restart(input_declaration, input_file,para_env,output_unit,error)
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(section_vals_type), POINTER         :: input_file
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'handle_ext_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(default_path_length)           :: binary_restart_file, &
                                                r_file_path
    CHARACTER(default_string_length)         :: path
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: restarted_infos
    INTEGER                                  :: ensemble, handle, i_rep_val, &
                                                iforce_eval, myi, n_rep_val, &
                                                nforce_eval1, nforce_eval2, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: ivec, iwalkers_status, iwork, &
                                                rwalkers_status
    LOGICAL                                  :: bsse_check, check, explicit1, &
                                                explicit2, failure, flag, &
                                                flag2, qmmm_check, &
                                                subsys_check
    REAL(KIND=dp)                            :: myt
    REAL(KIND=dp), DIMENSION(:), POINTER     :: vec, work
    TYPE(cp_error_type)                      :: sub_error
    TYPE(cp_parser_type), POINTER            :: cpparser
    TYPE(cp_unit_set_type), POINTER          :: default_units
    TYPE(section_vals_type), POINTER         :: r_section, rep_sections, &
                                                restart_file, section, &
                                                section1, section2, &
                                                sections1, sections2

    NULLIFY(restarted_infos,iwalkers_status,rwalkers_status,vec,ivec,work,iwork)
    CALL timeset(routineN,handle)
    failure=.FALSE.
    IF (.NOT. failure) THEN
       ! Handle restart file
       r_section => section_vals_get_subs_vals(input_file, "EXT_RESTART",error=error)
       CALL section_vals_val_get(r_section,"RESTART_FILE_NAME", c_val=r_file_path, error=error)
       CALL section_vals_val_get(r_section,"BINARY_RESTART_FILE_NAME",c_val=binary_restart_file,&
                                 error=error)

       IF (r_file_path/=" ") THEN
          CALL cp_error_init(sub_error,template_error=error,stop_level=cp_fatal_level)
          ! parse the input
          NULLIFY(default_units,restart_file)
          CALL section_vals_create(restart_file,input_declaration, error=sub_error)
          NULLIFY(cpparser)
          CALL parser_create(cpparser,file_name=r_file_path,para_env=para_env,error=sub_error)
          CALL cp_unit_set_create(default_units, "OUTPUT",error=sub_error)
          CALL section_vals_parse(restart_file,cpparser,root_section=.FALSE.,&
               default_units=default_units,error=sub_error)
          CALL cp_unit_set_release(default_units,error=sub_error)
          CALL parser_release(cpparser,error=sub_error)
          CALL cp_error_dealloc_ref(sub_error)

          ! Restart and input files same number of force_env sections
          sections1  => section_vals_get_subs_vals(restart_file,"FORCE_EVAL",error=error)
          CALL section_vals_get(sections1, n_repetition=nforce_eval1, error=error)
          sections2  => section_vals_get_subs_vals(input_file,"FORCE_EVAL",error=error)
          CALL section_vals_get(sections2, n_repetition=nforce_eval2, error=error)
          IF (nforce_eval1/=nforce_eval2) THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP,"Restart and input file MUST have the number of force_env sections",&
                  error,failure)
          END IF
          ! Handle default restarts
          CALL handle_defaults_restart(r_section, error=error)

          ! Real restart of force_evals
          DO iforce_eval = 1, nforce_eval1
             section1 => section_vals_get_subs_vals3(sections1,"SUBSYS",&
                  i_rep_section=iforce_eval,error=error)
             section2 => section_vals_get_subs_vals3(sections2,"SUBSYS",&
                  i_rep_section=iforce_eval,error=error)
             ! Some care needs to be handled when treating multiple force_eval
             ! Both subsys need to be consistently associated or not
             ! Mixed stuff will be rejected for safety reason..
             subsys_check = (ASSOCIATED(section1).EQV.ASSOCIATED(section2))
             IF (subsys_check) THEN
                IF (ASSOCIATED(section1)) THEN
                   CALL section_vals_val_get(r_section,"RESTART_CELL",l_val=flag,error=error)
                   IF (flag) THEN
                      section  => section_vals_get_subs_vals(section1,"CELL",error=error)
                      CALL section_vals_set_subs_vals(section2,"CELL",section,error=error)
                      CALL set_restart_info("CELL",restarted_infos,error)
                   END IF

                   CALL section_vals_val_get(r_section,"RESTART_POS",l_val=flag,error=error)
                   IF (flag) THEN
                      section => section_vals_get_subs_vals(section1,"COORD",error=error)
                      CALL section_vals_set_subs_vals(section2,"COORD",section,error=error)
                      CALL set_restart_info("COORDINATES",restarted_infos,error)
                      ! Copy over also the information on the multiple_unit_cell
                      CALL section_vals_val_get(section1,"TOPOLOGY%MULTIPLE_UNIT_CELL",i_vals=ivec,&
                                                error=error)
                      ALLOCATE(iwork(3),stat=stat)
                      CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                      iwork=ivec
                      CALL section_vals_val_set(section2,"TOPOLOGY%MULTIPLE_UNIT_CELL",i_vals_ptr=iwork,&
                                                error=error)
                   END IF

                   CALL section_vals_val_get(r_section,"RESTART_RANDOMG",l_val=flag,error=error)
                   IF (flag) THEN
                      section => section_vals_get_subs_vals(section1,"RNG_INIT",error=error)
                      CALL section_vals_set_subs_vals(section2,"RNG_INIT",section,error=error)
                      CALL set_restart_info("RANDOM NUMBER GENERATOR",restarted_infos,error)
                   END IF

                   CALL section_vals_val_get(r_section,"RESTART_VEL",l_val=flag,error=error)
                   IF (flag) THEN
                      section => section_vals_get_subs_vals(section1,"VELOCITY",error=error)
                      CALL section_vals_set_subs_vals(section2,"VELOCITY",section,error=error)
                      CALL set_restart_info("VELOCITIES",restarted_infos,error)
                   END IF

                   ! Core-Shell information "restarted" only when strictly necessary
                   CALL section_vals_val_get(r_section,"RESTART_SHELL_POS",l_val=flag,error=error)
                   IF (flag) THEN
                      section => section_vals_get_subs_vals(section1,"SHELL_COORD",error=error)
                      CALL section_vals_set_subs_vals(section2,"SHELL_COORD",section,error=error)
                      IF  (check_restart(section1, section2, "SHELL_COORD",error)) &
                           CALL set_restart_info("SHELL COORDINATES",restarted_infos,error)
                   END IF
                   CALL section_vals_val_get(r_section,"RESTART_CORE_POS",l_val=flag,error=error)
                   IF (flag) THEN
                      section => section_vals_get_subs_vals(section1,"CORE_COORD",error=error)
                      CALL section_vals_set_subs_vals(section2,"CORE_COORD",section,error=error)
                      IF  (check_restart(section1, section2, "CORE_COORD",error)) &
                           CALL set_restart_info("CORE COORDINATES",restarted_infos,error)
                   END IF
                   CALL section_vals_val_get(r_section,"RESTART_SHELL_VELOCITY",l_val=flag,error=error)
                   IF (flag) THEN
                      section => section_vals_get_subs_vals(section1,"SHELL_VELOCITY",error=error)
                      CALL section_vals_set_subs_vals(section2,"SHELL_VELOCITY",section,error=error)
                      IF  (check_restart(section1, section2, "SHELL_VELOCITY",error)) &
                           CALL set_restart_info("SHELL VELOCITIES",restarted_infos,error)
                   END IF
                   CALL section_vals_val_get(r_section,"RESTART_CORE_VELOCITY",l_val=flag,error=error)
                   IF (flag) THEN
                      section => section_vals_get_subs_vals(section1,"CORE_VELOCITY",error=error)
                      CALL section_vals_set_subs_vals(section2,"CORE_VELOCITY",section,error=error)
                      IF  (check_restart(section1, section2, "CORE_VELOCITY",error)) &
                           CALL set_restart_info("CORE VELOCITIES",restarted_infos,error)
                   END IF
                END IF
             ELSE
                CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                     routineP,"Error while reading the restart file. Two force_eval have incompatible"//&
                     " subsys.One of them has an allocated subsys while the other has not! Check your"//&
                     " input file or whether the restart file is compatible with the input!",&
                     error,failure)
             END IF
             ! QMMM restarts
             CALL section_vals_val_get(r_section,"RESTART_QMMM",l_val=flag,error=error)
             section1 => section_vals_get_subs_vals3(sections1,"QMMM",i_rep_section=iforce_eval,error=error)
             section2 => section_vals_get_subs_vals3(sections2,"QMMM",i_rep_section=iforce_eval,error=error)
             CALL section_vals_get(section1,explicit=explicit1,error=error)
             CALL section_vals_get(section2,explicit=explicit2,error=error)
             qmmm_check = (explicit1.AND.explicit2)
             IF (flag.AND.qmmm_check) THEN
                CALL set_restart_info("QMMM TRANSLATION VECTOR",restarted_infos,error)
                CALL section_vals_val_get(section1,"INITIAL_TRANSLATION_VECTOR",r_vals=vec,error=error)
                ALLOCATE(work(3),stat=stat)
                CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                work = vec
                CALL section_vals_val_set(section2,"INITIAL_TRANSLATION_VECTOR",r_vals_ptr=work,error=error)
             END IF
             ! BSSE restarts
             CALL section_vals_val_get(r_section,"RESTART_BSSE",l_val=flag,error=error)
             section1 => section_vals_get_subs_vals3(sections1,"BSSE",i_rep_section=iforce_eval,error=error)
             section2 => section_vals_get_subs_vals3(sections2,"BSSE",i_rep_section=iforce_eval,error=error)
             CALL section_vals_get(section1,explicit=explicit1,error=error)
             CALL section_vals_get(section2,explicit=explicit2,error=error)
             bsse_check = (explicit1.AND.explicit2)
             IF (flag.AND.bsse_check) THEN
                section  => section_vals_get_subs_vals(section1,"FRAGMENT_ENERGIES",error=error)
                CALL section_vals_set_subs_vals(section2,"FRAGMENT_ENERGIES",section,error=error)
                CALL set_restart_info("BSSE FRAGMENT ENERGIES",restarted_infos,error)
             END IF
          END DO

          CALL section_vals_val_get(r_section,"RESTART_COUNTERS",l_val=flag,error=error)
          IF (flag) THEN
             IF (check_restart(input_file, restart_file, "MOTION%MD",error)) THEN
                CALL section_vals_val_get(restart_file,"MOTION%MD%STEP_START_VAL",i_val=myi,error=error)
                CALL section_vals_val_set(input_file,"MOTION%MD%STEP_START_VAL",i_val=myi,error=error)
                CALL section_vals_val_get(restart_file,"MOTION%MD%TIME_START_VAL",r_val=myt,error=error)
                CALL section_vals_val_set(input_file,"MOTION%MD%TIME_START_VAL",r_val=myt,error=error)
                CALL section_vals_val_get(restart_file,"MOTION%MD%ECONS_START_VAL",r_val=myt,error=error)
                CALL section_vals_val_set(input_file,"MOTION%MD%ECONS_START_VAL",r_val=myt,error=error)
                CALL set_restart_info("MD COUNTERS",restarted_infos,error)
             END IF
             !
             IF (check_restart(input_file, restart_file, "MOTION%GEO_OPT",error)) THEN
                ! GEO_OPT
                CALL section_vals_val_get(restart_file,"MOTION%GEO_OPT%STEP_START_VAL",i_val=myi,error=error)
                CALL section_vals_val_set(input_file,"MOTION%GEO_OPT%STEP_START_VAL",i_val=myi,error=error)
                CALL set_restart_info("GEO_OPT COUNTERS",restarted_infos,error)
                ! ROT_OPT
                IF (check_restart(input_file, restart_file, "MOTION%GEO_OPT%TRANSITION_STATE%DIMER%ROT_OPT",error)) THEN
                   CALL section_vals_val_get(restart_file,"MOTION%GEO_OPT%TRANSITION_STATE%DIMER%ROT_OPT%STEP_START_VAL",&
                        i_val=myi,error=error)
                   CALL section_vals_val_set(input_file,"MOTION%GEO_OPT%TRANSITION_STATE%DIMER%ROT_OPT%STEP_START_VAL",&
                        i_val=myi,error=error)
                   CALL set_restart_info("ROT_OPT COUNTERS",restarted_infos,error)
                END IF
             END IF
             !
             IF (check_restart(input_file, restart_file, "MOTION%GEO_OPT",error)) THEN
                ! CELL_OPT
                CALL section_vals_val_get(restart_file,"MOTION%CELL_OPT%STEP_START_VAL",i_val=myi,error=error)
                CALL section_vals_val_set(input_file,"MOTION%CELL_OPT%STEP_START_VAL",i_val=myi,error=error)
                CALL set_restart_info("CELL_OPT COUNTERS",restarted_infos,error)
             END IF
             !
             IF (check_restart(input_file, restart_file, "OPTIMIZE_INPUT",error)) THEN
                CALL section_vals_val_get(restart_file,"OPTIMIZE_INPUT%ITER_START_VAL",i_val=myi,error=error)
                CALL section_vals_val_set(input_file,"OPTIMIZE_INPUT%ITER_START_VAL",i_val=myi,error=error)
                CALL set_restart_info("OPTIMIZE_INPUT ITERATION NUMBER",restarted_infos,error)
             END IF
             !
             IF (check_restart(input_file, restart_file, "MOTION%PINT",error)) THEN
                ! PINT
                CALL section_vals_val_get(restart_file,"MOTION%PINT%ITERATION",i_val=myi,error=error)
                CALL section_vals_val_set(input_file,"MOTION%PINT%ITERATION",i_val=myi,error=error)
                CALL set_restart_info("PINT ITERATION NUMBER",restarted_infos,error)
             END IF
             !
             CALL section_vals_val_get(r_section,"RESTART_METADYNAMICS",l_val=flag2,error=error)
             IF (flag2.AND.check_restart(input_file, restart_file, "MOTION%FREE_ENERGY%METADYN",error)) THEN
                CALL section_vals_val_get(restart_file,&
                     "MOTION%FREE_ENERGY%METADYN%STEP_START_VAL",i_val=myi,error=error)
                CALL section_vals_val_set(input_file,&
                     "MOTION%FREE_ENERGY%METADYN%STEP_START_VAL",i_val=myi,error=error)
                CALL section_vals_val_get(restart_file,&
                     "MOTION%FREE_ENERGY%METADYN%NHILLS_START_VAL",i_val=myi,error=error)
                CALL section_vals_val_set(input_file,&
                     "MOTION%FREE_ENERGY%METADYN%NHILLS_START_VAL",i_val=myi,error=error)
                !RG Adaptive hills
                CALL section_vals_val_get(restart_file,&
                     "MOTION%FREE_ENERGY%METADYN%OLD_HILL_NUMBER",i_val=myi,error=error)
                CALL section_vals_val_set(input_file,&
                     "MOTION%FREE_ENERGY%METADYN%OLD_HILL_NUMBER",i_val=myi,error=error)
                CALL section_vals_val_get(restart_file,&
                     "MOTION%FREE_ENERGY%METADYN%OLD_HILL_STEP",i_val=myi,error=error)
                CALL section_vals_val_set(input_file,&
                     "MOTION%FREE_ENERGY%METADYN%OLD_HILL_STEP",i_val=myi,error=error)
                !RG Adaptive hills
                CALL set_restart_info("METADYNAMIC COUNTERS",restarted_infos,error)
             END IF
          END IF

          CALL section_vals_val_get(r_section,"RESTART_AVERAGES",l_val=flag,error=error)
          IF (flag) THEN
             IF (check_restart(input_file, restart_file, "MOTION%MD",error)) THEN
                rep_sections => section_vals_get_subs_vals(restart_file,"MOTION%MD%AVERAGES%RESTART_AVERAGES",error=error)
                CALL section_vals_set_subs_vals(input_file,"MOTION%MD%AVERAGES%RESTART_AVERAGES",rep_sections,error=error)
                CALL set_restart_info("MD AVERAGES",restarted_infos,error)
             END IF
          END IF

          CALL section_vals_val_get(r_section,"RESTART_BAND",l_val=flag,error=error)
          IF (flag.AND.check_restart(input_file, restart_file, "MOTION%BAND",error)) THEN
             rep_sections => section_vals_get_subs_vals(restart_file,"MOTION%BAND%REPLICA",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%BAND%REPLICA",rep_sections,error=error)
             CALL set_restart_info("BAND CALCULATION",restarted_infos,error)
          END IF

          CALL section_vals_val_get(r_section,"RESTART_OPTIMIZE_INPUT_VARIABLES",l_val=flag,error=error)
          IF (flag.AND.check_restart(input_file, restart_file, "OPTIMIZE_INPUT%VARIABLE",error)) THEN
             rep_sections => section_vals_get_subs_vals(restart_file,"OPTIMIZE_INPUT%VARIABLE",error=error)
             CALL section_vals_set_subs_vals(input_file,"OPTIMIZE_INPUT%VARIABLE",rep_sections,error=error)
             CALL set_restart_info("OPTIMIZE_INPUT: VARIABLES",restarted_infos,error)
          END IF

          CALL section_vals_val_get(r_section,"RESTART_BAROSTAT",l_val=flag,error=error)
          IF (flag.AND.check_restart(input_file, restart_file, "MOTION%MD%BAROSTAT",error)) THEN
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%MD%BAROSTAT%MASS",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%MD%BAROSTAT%MASS",&
                  section,error=error)
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%MD%BAROSTAT%VELOCITY",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%MD%BAROSTAT%VELOCITY",&
                  section,error=error)
             CALL set_restart_info("BAROSTAT",restarted_infos,error)
          END IF

          flag = check_restart(input_file, restart_file, "MOTION%MD", error)
          IF(flag) THEN
             CALL section_vals_val_get(input_file,"MOTION%MD%ENSEMBLE",i_val=ensemble,error=error)
             IF(ensemble == npt_i_ensemble .OR. ensemble==npt_f_ensemble) THEN
                CALL section_vals_val_get(r_section,"RESTART_BAROSTAT_THERMOSTAT",l_val=flag,error=error)
                check = check_restart(input_file, restart_file, "MOTION%MD%BAROSTAT", error)
                CALL restart_thermostat(flag, input_file, restart_file, "MOTION%MD%BAROSTAT%THERMOSTAT", &
                     check=check, error=error)
                IF (flag.AND.check) CALL set_restart_info("THERMOSTAT OF BAROSTAT",restarted_infos,error)
             END IF
          END IF

          check = check_restart(input_file, restart_file, "MOTION%MD%SHELL", error)
          IF(check) THEN
             CALL section_vals_val_get(r_section,"RESTART_SHELL_THERMOSTAT",l_val=flag,error=error)
             CALL restart_thermostat(flag, input_file, restart_file, "MOTION%MD%SHELL%THERMOSTAT", error=error)
             CALL set_restart_info("SHELL THERMOSTAT",restarted_infos,error)
          END IF

          CALL section_vals_val_get(r_section,"RESTART_THERMOSTAT",l_val=flag,error=error)
          CALL restart_thermostat(flag,input_file, restart_file, "MOTION%MD%THERMOSTAT", error=error)
          IF (flag) CALL set_restart_info("PARTICLE THERMOSTAT",restarted_infos,error)

          CALL section_vals_val_get(r_section,"RESTART_CONSTRAINT",l_val=flag,error=error)
          IF (flag.AND.check_restart(input_file, restart_file, "MOTION%CONSTRAINT",error)) THEN
             section => section_vals_get_subs_vals(restart_file,"MOTION%CONSTRAINT",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%CONSTRAINT",section,error=error)
             CALL set_restart_info("CONSTRAINTS/RESTRAINTS",restarted_infos,error)
          END IF

          CALL section_vals_val_get(r_section,"RESTART_METADYNAMICS",l_val=flag,error=error)
          IF (flag.AND.check_restart(input_file, restart_file, "MOTION%FREE_ENERGY%METADYN",error)) THEN
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_POS",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_POS",&
                  section,error=error)
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_SCALE",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_SCALE",&
                  section,error=error)
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_HEIGHT",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_HEIGHT",&
                  section,error=error)
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_INVDT",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_INVDT",&
                  section,error=error)
             ! Extended Lagrangian
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%FREE_ENERGY%METADYN%EXT_LAGRANGE_SS0",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%FREE_ENERGY%METADYN%EXT_LAGRANGE_SS0",&
                  section,error=error)
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%FREE_ENERGY%METADYN%EXT_LAGRANGE_VVP",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%FREE_ENERGY%METADYN%EXT_LAGRANGE_VVP",&
                  section,error=error)
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%FREE_ENERGY%METADYN%EXT_LAGRANGE_SS",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%FREE_ENERGY%METADYN%EXT_LAGRANGE_SS",&
                  section,error=error)
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%FREE_ENERGY%METADYN%EXT_LAGRANGE_FS",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%FREE_ENERGY%METADYN%EXT_LAGRANGE_FS",&
                  section,error=error)
             CALL set_restart_info("METADYNAMICS",restarted_infos,error)
          END IF

          CALL section_vals_val_get(r_section,"RESTART_WALKERS",l_val=flag,error=error)
          IF (flag.AND.check_restart(input_file, restart_file, "MOTION%FREE_ENERGY%METADYN%MULTIPLE_WALKERS",error)) THEN
             CALL section_vals_val_get(restart_file,"MOTION%FREE_ENERGY%METADYN%MULTIPLE_WALKERS%WALKERS_STATUS",&
                  i_vals=rwalkers_status, error=error)
             ALLOCATE(iwalkers_status(SIZE(rwalkers_status)),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             iwalkers_status = rwalkers_status
             CALL section_vals_val_set(input_file,"MOTION%FREE_ENERGY%METADYN%MULTIPLE_WALKERS%WALKERS_STATUS",&
                  i_vals_ptr=iwalkers_status, error=error)
             CALL set_restart_info("WALKERS INFO",restarted_infos,error)
          END IF

          CALL section_vals_val_get(r_section,"RESTART_DIMER",l_val=flag,error=error)
          IF (flag.AND.check_restart(input_file, restart_file, "MOTION%GEO_OPT%TRANSITION_STATE%DIMER",error)) THEN
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%GEO_OPT%TRANSITION_STATE%DIMER%DIMER_VECTOR",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%GEO_OPT%TRANSITION_STATE%DIMER%DIMER_VECTOR",&
                  section,error=error)
             CALL set_restart_info("DIMER TRANSITION STATE SEARCH",restarted_infos,error)
          END IF

          CALL section_vals_val_get(r_section,"CUSTOM_PATH",n_rep_val=n_rep_val,error=error)
          DO i_rep_val=1,n_rep_val
             CALL section_vals_val_get(r_section,"CUSTOM_PATH",i_rep_val=i_rep_val,c_val=path,error=error)
             IF (path/=" ") THEN
                section => section_vals_get_subs_vals(restart_file,path,error=error)
                CALL section_vals_set_subs_vals(input_file,path,section,error=error)
                CALL set_restart_info("USER RESTART: "//TRIM(path),restarted_infos,error)
             END IF
          END DO

          CALL section_vals_val_get(r_section,"RESTART_RTP",l_val=flag,error=error)
!          IF(flag.AND.check_restart(input_file, restart_file, "FORCE_EVAL%DFT%REAL_TIME_PROPAGATION",error)) THEN
          IF(flag)THEN
               section => section_vals_get_subs_vals(restart_file,&
                  "FORCE_EVAL%DFT%REAL_TIME_PROPAGATION",error=error)
               CALL section_vals_val_get(section,"INITIAL_WFN",i_val=myi,error=error)
               CALL section_vals_val_set(input_file,"FORCE_EVAL%DFT%REAL_TIME_PROPAGATION%INITIAL_WFN",&
                  i_val=myi,error=error)
               CALL set_restart_info("REAL TIME PROPAGATION",restarted_infos,error)
          END IF

          ! PIMD
          CALL section_vals_val_get(r_section,"RESTART_PINT_POS",l_val=flag,error=error)
          IF(flag) THEN
               section => section_vals_get_subs_vals(restart_file,"MOTION%PINT%BEADS%COORD",error=error)
               CALL section_vals_set_subs_vals(input_file,"MOTION%PINT%BEADS%COORD",section,error=error)
               CALL set_restart_info("PINT BEAD POSITIONS",restarted_infos,error)
          END IF
          CALL section_vals_val_get(r_section,"RESTART_PINT_VEL",l_val=flag,error=error)
          IF(flag) THEN
               section => section_vals_get_subs_vals(restart_file,"MOTION%PINT%BEADS%VELOCITY",error=error)
               CALL section_vals_set_subs_vals(input_file,"MOTION%PINT%BEADS%VELOCITY",section,error=error)
               CALL set_restart_info("PINT BEAD VELOCITIES",restarted_infos,error)
          END IF
          CALL section_vals_val_get(r_section,"RESTART_PINT_NOSE",l_val=flag,error=error)
          IF(flag) THEN
               section => section_vals_get_subs_vals(restart_file,"MOTION%PINT%NOSE%COORD",error=error)
               CALL section_vals_set_subs_vals(input_file,"MOTION%PINT%NOSE%COORD",section,error=error)
               section => section_vals_get_subs_vals(restart_file,"MOTION%PINT%NOSE%VELOCITY",error=error)
               CALL section_vals_set_subs_vals(input_file,"MOTION%PINT%NOSE%VELOCITY",section,error=error)
               CALL set_restart_info("PINT NOSE THERMOSTAT",restarted_infos,error)
          END IF
          CALL section_vals_val_get(r_section,"RESTART_PINT_GLE",l_val=flag,error=error)
          IF(flag) THEN
               section => section_vals_get_subs_vals(restart_file,"MOTION%PINT%GLE",error=error)
               CALL section_vals_set_subs_vals(input_file,"MOTION%PINT%GLE",section,error=error)
               CALL set_restart_info("PINT GLE THERMOSTAT",restarted_infos,error)
          END IF


          ! PIMC
          !
          CALL section_vals_val_get(r_section,"RESTART_HELIUM_POS",l_val=flag,error=error)
          IF(flag) THEN
               CALL section_vals_val_get(input_file,"MOTION%PINT%HELIUM%NUM_ENV",&
                    explicit=explicit1, error=error)
               IF ( .NOT. explicit1 ) THEN
                    CALL section_vals_val_get(restart_file,"MOTION%PINT%HELIUM%NUM_ENV",i_val=myi,error=error)
                    CALL section_vals_val_set(input_file,"MOTION%PINT%HELIUM%NUM_ENV",i_val=myi,error=error)
               END IF
               section => section_vals_get_subs_vals(restart_file,"MOTION%PINT%HELIUM%COORD",error=error)
               CALL section_vals_set_subs_vals(input_file,"MOTION%PINT%HELIUM%COORD",section,error=error)
               CALL set_restart_info("HELIUM BEAD POSITIONS",restarted_infos,error)
          END IF
          !
          CALL section_vals_val_get(r_section,"RESTART_HELIUM_PERMUTATION",l_val=flag,error=error)
          IF(flag) THEN
               CALL section_vals_val_get(input_file,"MOTION%PINT%HELIUM%NUM_ENV",&
                    explicit=explicit1, error=error)
               IF ( .NOT. explicit1 ) THEN
                    CALL section_vals_val_get(restart_file,"MOTION%PINT%HELIUM%NUM_ENV",i_val=myi,error=error)
                    CALL section_vals_val_set(input_file,"MOTION%PINT%HELIUM%NUM_ENV",i_val=myi,error=error)
               END IF
               section => section_vals_get_subs_vals(restart_file,"MOTION%PINT%HELIUM%PERM",error=error)
               CALL section_vals_set_subs_vals(input_file,"MOTION%PINT%HELIUM%PERM",section,error=error)
               CALL set_restart_info("HELIUM PERMUTATION STATE",restarted_infos,error)
          END IF
          !
          CALL section_vals_val_get(r_section,"RESTART_HELIUM_FORCE",l_val=flag,error=error)
          IF(flag) THEN
               CALL section_vals_val_get(input_file,"MOTION%PINT%HELIUM%NUM_ENV",&
                    explicit=explicit1, error=error)
               IF ( .NOT. explicit1 ) THEN
                    CALL section_vals_val_get(restart_file,"MOTION%PINT%HELIUM%NUM_ENV",i_val=myi,error=error)
                    CALL section_vals_val_set(input_file,"MOTION%PINT%HELIUM%NUM_ENV",i_val=myi,error=error)
               END IF
               section => section_vals_get_subs_vals(restart_file,"MOTION%PINT%HELIUM%FORCE",error=error)
               CALL section_vals_set_subs_vals(input_file,"MOTION%PINT%HELIUM%FORCE",section,error=error)
               CALL set_restart_info("HELIUM FORCES ON SOLUTE",restarted_infos,error)
          END IF
          !
          CALL section_vals_val_get(r_section,"RESTART_HELIUM_RNG",l_val=flag,error=error)
          IF(flag) THEN
               CALL section_vals_val_get(input_file,"MOTION%PINT%HELIUM%NUM_ENV",&
                    explicit=explicit1, error=error)
               IF ( .NOT. explicit1 ) THEN
                    CALL section_vals_val_get(restart_file,"MOTION%PINT%HELIUM%NUM_ENV",i_val=myi,error=error)
                    CALL section_vals_val_set(input_file,"MOTION%PINT%HELIUM%NUM_ENV",i_val=myi,error=error)
               END IF
               section => section_vals_get_subs_vals(restart_file,"MOTION%PINT%HELIUM%RNG_STATE",error=error)
               CALL section_vals_set_subs_vals(input_file,"MOTION%PINT%HELIUM%RNG_STATE",section,error=error)
               CALL set_restart_info("HELIUM RNG STATE",restarted_infos,error)
          END IF
          !
          !
          CALL section_vals_val_get(r_section,"RESTART_HELIUM_DENSITIES",l_val=flag,error=error)
          IF(flag) THEN
               CALL section_vals_val_get(input_file,"MOTION%PINT%HELIUM%NUM_ENV",&
                    explicit=explicit1, error=error)
               IF ( .NOT. explicit1 ) THEN
                    CALL section_vals_val_get(restart_file,"MOTION%PINT%HELIUM%NUM_ENV",i_val=myi,error=error)
                    CALL section_vals_val_set(input_file,"MOTION%PINT%HELIUM%NUM_ENV",i_val=myi,error=error)
               END IF
               section => section_vals_get_subs_vals(restart_file,"MOTION%PINT%HELIUM%RHO",error=error)
               CALL section_vals_set_subs_vals(input_file,"MOTION%PINT%HELIUM%RHO",section,error=error)
               CALL set_restart_info("HELIUM DENSITIES",restarted_infos,error)
          END IF
          !
          CALL section_vals_val_set(r_section,"RESTART_FILE_NAME", c_val=" ",error=error)
          CALL section_vals_release(restart_file,error=error)
          CALL release_restart_info(restarted_infos, r_file_path, binary_restart_file,&
                                    output_unit, error)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE handle_ext_restart

! *****************************************************************************
!> \brief store information on the restarted quantities
!> \param label ...
!> \param restarted_infos ...
!> \param error ...
!> \author Teodoro Laino [tlaino] 09.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE set_restart_info(label, restarted_infos, error)

    CHARACTER(LEN=*), INTENT(IN)             :: label
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: restarted_infos
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_restart_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: isize
    LOGICAL                                  :: failure

    failure = .FALSE.
    isize = 0
    IF (ASSOCIATED(restarted_infos)) isize = SIZE(restarted_infos)
    isize = isize + 1
    CALL reallocate(restarted_infos,1,isize)
    restarted_infos(isize) = TRIM(label)

  END SUBROUTINE set_restart_info

! *****************************************************************************
!> \brief dumps on output the information on the information effectively restarted
!> \param restarted_infos ...
!> \param r_file_path ...
!> \param binary_restart_file ...
!> \param output_unit ...
!> \param error ...
!> \author Teodoro Laino [tlaino] 09.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE release_restart_info(restarted_infos, r_file_path,&
                                  binary_restart_file, output_unit, error)
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: restarted_infos
    CHARACTER(LEN=*), INTENT(IN)             :: r_file_path, &
                                                binary_restart_file
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_restart_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (output_unit>0.AND.ASSOCIATED(restarted_infos)) THEN
       WRITE(output_unit,'(1X,79("*"))')
       WRITE(output_unit,'(1X,"*",T30,A,T80,"*")')" RESTART INFORMATION "
       WRITE(output_unit,'(1X,79("*"))')
       WRITE(output_unit,'(1X,"*",T80,"*")')
       i = 1
       WRITE(output_unit,'(1X,"*",A,T26,A,T80,"*")')"    RESTART FILE NAME:  ",&
            r_file_path(53*(i-1)+1:53*i)
       DO i = 2, CEILING(REAL(LEN_TRIM(r_file_path),KIND=dp)/53.0_dp)
          WRITE(output_unit,'(T1,1X,"*",T26,A,T80,"*")')r_file_path(53*(i-1)+1:53*i)
       END DO
       IF (LEN_TRIM(binary_restart_file) > 0) THEN
          i = 1
          WRITE (output_unit,'(1X,"*",A,T26,A,T80,"*")')"  BINARY RESTART FILE:  ",&
            binary_restart_file(53*(i-1)+1:53*i)
          DO i = 2, CEILING(REAL(LEN_TRIM(binary_restart_file),KIND=dp)/53.0_dp)
             WRITE(output_unit,'(T1,1X,"*",T26,A,T80,"*")') binary_restart_file(53*(i-1)+1:53*i)
          END DO
       END IF
       WRITE(output_unit,'(1X,"*",T80,"*")')
       WRITE(output_unit,'(1X,"*", A,T80,"*")')" RESTARTED QUANTITIES:  "
       DO j = 1, SIZE(restarted_infos)
          DO i = 1, CEILING(REAL(LEN_TRIM(restarted_infos(j)),KIND=dp)/53.0_dp)
             WRITE(output_unit,'(T1,1X,"*",T26,A,T80,"*")')restarted_infos(j)(53*(i-1)+1:53*i)
          END DO
       END DO
       WRITE(output_unit,'(1X,79("*"),/)')
    END IF
    IF (ASSOCIATED(restarted_infos)) THEN
       DEALLOCATE(restarted_infos,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE release_restart_info


! *****************************************************************************
!> \brief Possibly restart thermostats information
!> \param flag ...
!> \param input_file the input file to initialize
!> \param restart_file ...
!> \param path ...
!> \param check ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] 10.2007- University of Zurich
! *****************************************************************************
  SUBROUTINE restart_thermostat(flag, input_file, restart_file, path, check, error)
    LOGICAL, INTENT(IN)                      :: flag
    TYPE(section_vals_type), POINTER         :: input_file, restart_file
    CHARACTER(LEN=*), INTENT(IN)             :: path
    LOGICAL, INTENT(IN), OPTIONAL            :: check
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'restart_thermostat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: input_region, input_type, &
                                                restart_region, restart_type
    LOGICAL                                  :: check_loc, failure, &
                                                skip_other_checks
    TYPE(section_vals_type), POINTER         :: section

    failure   = .FALSE.
    check_loc = check_restart(input_file, restart_file, TRIM(path), error)
    skip_other_checks = PRESENT(check)
    IF (skip_other_checks) check_loc = check
    IF (flag.AND.check_loc) THEN
       ! Let's check if the thermostat type is different otherwise it does not make any
       ! sense to do any kind of restart
       CALL section_vals_val_get(input_file,TRIM(path)//"%TYPE",i_val=input_type,error=error)
       CALL section_vals_val_get(restart_file,TRIM(path)//"%TYPE",i_val=restart_type,error=error)

       IF (input_type==do_thermo_same_as_part) THEN
          CALL section_vals_val_get(input_file,"MOTION%MD%THERMOSTAT%TYPE",i_val=input_type,error=error)
       END IF

       IF (skip_other_checks) THEN
          input_region   = do_region_global
          restart_region = do_region_global
       ELSE
          ! Also the regions must be the same..
          CALL section_vals_val_get(input_file,TRIM(path)//"%REGION",i_val=input_region,error=error)
          CALL section_vals_val_get(restart_file,TRIM(path)//"%REGION",i_val=restart_region,error=error)
       END IF

       IF ((input_type==restart_type).AND.(input_region==restart_region)) THEN
          SELECT CASE(input_type)
          CASE(do_thermo_nose)
             section => section_vals_get_subs_vals(restart_file,TRIM(path)//"%NOSE%COORD",error=error)
             CALL section_vals_set_subs_vals(input_file,TRIM(path)//"%NOSE%COORD",section,error=error)

             section => section_vals_get_subs_vals(restart_file,TRIM(path)//"%NOSE%VELOCITY",error=error)
             CALL section_vals_set_subs_vals(input_file,TRIM(path)//"%NOSE%VELOCITY",section,error=error)

             section => section_vals_get_subs_vals(restart_file,TRIM(path)//"%NOSE%MASS",error=error)
             CALL section_vals_set_subs_vals(input_file,TRIM(path)//"%NOSE%MASS",section,error=error)

             section => section_vals_get_subs_vals(restart_file,TRIM(path)//"%NOSE%FORCE",error=error)
             CALL section_vals_set_subs_vals(input_file,TRIM(path)//"%NOSE%FORCE",section,error=error)
          CASE(do_thermo_csvr)
             section => section_vals_get_subs_vals(restart_file,TRIM(path)//"%CSVR%THERMOSTAT_ENERGY",error=error)
             CALL section_vals_set_subs_vals(input_file,TRIM(path)//"%CSVR%THERMOSTAT_ENERGY",section,error=error)
             section => section_vals_get_subs_vals(restart_file,TRIM(path)//"%CSVR%RNG_INIT",error=error)
             CALL section_vals_set_subs_vals(input_file,TRIM(path)//"%CSVR%RNG_INIT",section,error=error)
          CASE(do_thermo_gle)
             section => section_vals_get_subs_vals(restart_file,TRIM(path)//"%GLE%THERMOSTAT_ENERGY",error=error)
             CALL section_vals_set_subs_vals(input_file,TRIM(path)//"%GLE%THERMOSTAT_ENERGY",section,error=error)
             section => section_vals_get_subs_vals(restart_file,TRIM(path)//"%GLE%RNG_INIT",error=error)
             CALL section_vals_set_subs_vals(input_file,TRIM(path)//"%GLE%RNG_INIT",section,error=error)
             section => section_vals_get_subs_vals(restart_file,TRIM(path)//"%GLE%S",error=error)
             CALL section_vals_set_subs_vals(input_file,TRIM(path)//"%GLE%S",section,error=error)
          CASE(do_thermo_al)
             section => section_vals_get_subs_vals(restart_file,TRIM(path)//"%AD_LANGEVIN%CHI",error=error)
             CALL section_vals_set_subs_vals(input_file,TRIM(path)//"%AD_LANGEVIN%CHI",section,error=error)
             section => section_vals_get_subs_vals(restart_file,TRIM(path)//"%AD_LANGEVIN%MASS",error=error)
             CALL section_vals_set_subs_vals(input_file,TRIM(path)//"%AD_LANGEVIN%MASS",section,error=error)
          END SELECT
       ELSE
          CALL cp_assert((input_type==restart_type),cp_warning_level,cp_assertion_failed,routineP,&
                        "Requested to restart thermostat: "//TRIM(path)//". The thermostat "//&
                        "specified in the input file and the information present in the restart "//&
                        "file do not match the same type of thermostat! Restarting is not possible! "//&
                        "Thermostat will not be restarted! "//&
CPSourceFileRef,&
                        only_ionode=.TRUE.)
          CALL cp_assert((input_region==restart_region),cp_warning_level,cp_assertion_failed,routineP,&
                        "Requested to restart thermostat: "//TRIM(path)//". The thermostat "//&
                        "specified in the input file and the information present in the restart "//&
                        "file do not match the same type of REGION! Restarting is not possible! "//&
                        "Thermostat will not be restarted! "//&
CPSourceFileRef,&
                        only_ionode=.TRUE.)
       END IF
    END IF
  END SUBROUTINE restart_thermostat

! *****************************************************************************
!> \brief Checks if there are the proper conditions to do a restart
!> \param input_file the input file to initialize
!> \param restart_file ...
!> \param tag_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval do_restart ...
!> \author teo
! *****************************************************************************
  FUNCTION check_restart(input_file, restart_file, tag_section, error) RESULT(do_restart)
    TYPE(section_vals_type), POINTER         :: input_file, restart_file
    CHARACTER(LEN=*), INTENT(IN)             :: tag_section
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: do_restart

    CHARACTER(len=*), PARAMETER :: routineN = 'check_restart', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: explicit1, explicit2, failure
    TYPE(section_vals_type), POINTER         :: work_section

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(work_section)
    IF (.NOT. failure) THEN
       work_section => section_vals_get_subs_vals(input_file,TRIM(tag_section),error=error)
       CALL section_vals_get(work_section,explicit=explicit1,error=error)
       work_section => section_vals_get_subs_vals(restart_file,TRIM(tag_section),error=error)
       CALL section_vals_get(work_section,explicit=explicit2,error=error)

       do_restart = explicit1.AND.explicit2
    END IF
    CALL timestop(handle)
  END FUNCTION check_restart

! *****************************************************************************
!> \brief Removes section used to restart a calculation from an
!>      input file in memory
!> \param input_file the input file to initialize
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE remove_restart_info(input_file,error)
    TYPE(section_vals_type), POINTER         :: input_file
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'remove_restart_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iforce_eval, &
                                                nforce_eval1
    LOGICAL                                  :: explicit1, failure
    TYPE(section_vals_type), POINTER         :: md_section, motion_section, &
                                                section1, section_to_delete, &
                                                sections1, work_section

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(work_section)
    IF (.NOT. failure) THEN
       section_to_delete => section_vals_get_subs_vals(input_file,"EXT_RESTART",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       sections1  => section_vals_get_subs_vals(input_file,"FORCE_EVAL",error=error)
       CALL section_vals_get(sections1, n_repetition=nforce_eval1, error=error)

       DO iforce_eval = 1, nforce_eval1
          section1 => section_vals_get_subs_vals3(sections1,"SUBSYS",i_rep_section=iforce_eval,error=error)
          section_to_delete => section_vals_get_subs_vals(section1,"COORD",error=error)
          CALL section_vals_remove_values(section_to_delete,error)
          section_to_delete => section_vals_get_subs_vals(section1,"VELOCITY",error=error)
          CALL section_vals_remove_values(section_to_delete,error)
       END DO

       motion_section => section_vals_get_subs_vals(input_file,"MOTION",error=error)
       md_section => section_vals_get_subs_vals(motion_section,"MD",error=error)
       CALL section_vals_get(md_section,explicit=explicit1,error=error)
       IF (explicit1) THEN
          CALL section_vals_val_unset(md_section,"STEP_START_VAL",error=error)
          CALL section_vals_val_unset(md_section,"TIME_START_VAL",error=error)
          CALL section_vals_val_unset(md_section,"ECONS_START_VAL",error=error)
       END IF
       work_section => section_vals_get_subs_vals(motion_section,"FREE_ENERGY%METADYN",error=error)
       CALL section_vals_get(work_section,explicit=explicit1,error=error)
       IF (explicit1) THEN
          CALL section_vals_val_unset(motion_section,"FREE_ENERGY%METADYN%STEP_START_VAL",error=error)
          CALL section_vals_val_unset(motion_section,"FREE_ENERGY%METADYN%NHILLS_START_VAL",error=error)
       END IF
       section_to_delete => section_vals_get_subs_vals(motion_section,"BAND%REPLICA",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(md_section,"AVERAGES%RESTART_AVERAGES",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(md_section,"THERMOSTAT%NOSE%COORD",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(md_section,"THERMOSTAT%NOSE%VELOCITY",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(md_section,"THERMOSTAT%NOSE%MASS",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(md_section,"THERMOSTAT%NOSE%FORCE",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(md_section,"BAROSTAT%MASS",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(md_section,"BAROSTAT%VELOCITY",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(md_section,"BAROSTAT%THERMOSTAT%NOSE%COORD",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(md_section,"BAROSTAT%THERMOSTAT%NOSE%VELOCITY",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(md_section,"BAROSTAT%THERMOSTAT%NOSE%MASS",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(md_section,"BAROSTAT%THERMOSTAT%NOSE%FORCE",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(md_section,"SHELL%THERMOSTAT%NOSE%COORD",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(md_section,"SHELL%THERMOSTAT%NOSE%VELOCITY",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(md_section,"SHELL%THERMOSTAT%NOSE%MASS",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(md_section,"SHELL%THERMOSTAT%NOSE%FORCE",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       ! Constrained/Restrained section
       section_to_delete => section_vals_get_subs_vals(motion_section,"CONSTRAINT%FIX_ATOM_RESTART",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(motion_section,"CONSTRAINT%COLVAR_RESTART",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       ! Free energies restarts
       section_to_delete => section_vals_get_subs_vals(motion_section,"FREE_ENERGY%METADYN%SPAWNED_HILLS_POS",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(motion_section,"FREE_ENERGY%METADYN%SPAWNED_HILLS_SCALE",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(motion_section,"FREE_ENERGY%METADYN%SPAWNED_HILLS_HEIGHT",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(motion_section,"FREE_ENERGY%METADYN%SPAWNED_HILLS_INVDT",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(motion_section,"FREE_ENERGY%METADYN%EXT_LAGRANGE_SS0",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(motion_section,"FREE_ENERGY%METADYN%EXT_LAGRANGE_VVP",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(motion_section,"FREE_ENERGY%METADYN%EXT_LAGRANGE_SS",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(motion_section,"FREE_ENERGY%METADYN%EXT_LAGRANGE_FS",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
    END IF
    CALL timestop(handle)
  END SUBROUTINE remove_restart_info

! *****************************************************************************
!> \brief This subroutine controls the defaults for the restartable quantities..
!> \param r_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo - University of Zurich - 09.2007 [tlaino]
! *****************************************************************************
  SUBROUTINE handle_defaults_restart(r_section,error)
    TYPE(section_vals_type), POINTER         :: r_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'handle_defaults_restart', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ik, nval
    LOGICAL                                  :: failure, restart_default
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section

    CALL timeset(routineN,handle)
    failure=.FALSE.
    IF (.NOT. failure) THEN
       NULLIFY(keyword, section)
       CALL section_vals_get(r_section, section=section, error=error)
       CALL section_vals_val_get(r_section,"RESTART_DEFAULT",l_val=restart_default,error=error)
       DO ik=-1,section%n_keywords
          keyword => section%keywords(ik)%keyword
          IF (ASSOCIATED(keyword)) THEN
             IF (keyword%type_of_var==logical_t.AND.keyword%names(1)(1:8)=="RESTART_") THEN
                IF (TRIM(keyword%names(1))=="RESTART_DEFAULT") CYCLE
                CALL section_vals_val_get(r_section,keyword%names(1),n_rep_val=nval,error=error)
                IF (nval==0) THEN
                   ! User didn't specify any value, use the value of the RESTART_DEFAULT keyword..
                   CALL section_vals_val_set(r_section,keyword%names(1),l_val=restart_default,error=error)
                END IF
             END IF
          END IF
       END DO
    END IF
   CALL timestop(handle)

  END SUBROUTINE handle_defaults_restart

END MODULE input_cp2k_check
