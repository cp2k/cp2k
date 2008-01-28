!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_cp2k_check *
!!
!!   NAME
!!     input_cp2k_check
!!
!!   FUNCTION
!!     checks the input and perform some automatic "magic" on it
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     01.2006 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE input_cp2k_check
  USE cp2k_info,                       ONLY: id_cp2k_version
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser,                       ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_set_create,&
                                             cp_unit_set_release,&
                                             cp_unit_set_type
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_region_global, do_thermo_csvr, do_thermo_nose, &
       do_thermo_same_as_part, do_vwn5, id_release_version, xc_funct_b3lyp, &
       xc_funct_blyp, xc_funct_bp, xc_funct_hcth120, xc_funct_no_shortcut, &
       xc_funct_olyp, xc_funct_pade, xc_funct_pbe, xc_funct_pbe0, &
       xc_funct_tpss, xc_funct_xwpbe, xc_none
  USE input_cp2k,                      ONLY: create_cp2k_root_section
  USE input_keyword_types,             ONLY: keyword_type
  USE input_parsing,                   ONLY: section_vals_parse
  USE input_section_types,             ONLY: &
       section_release, section_type, section_vals_check_release, &
       section_vals_create, section_vals_get, section_vals_get_subs_vals, &
       section_vals_get_subs_vals3, section_vals_release, &
       section_vals_remove_values, section_vals_set_subs_vals, &
       section_vals_type, section_vals_val_get, section_vals_val_set, &
       section_vals_val_unset
  USE input_val_types,                 ONLY: logical_t
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_check'

  PUBLIC :: check_cp2k_input, xc_functionals_expand, remove_restart_info

!!***
!****************************************************************************
CONTAINS

!!****f* input_cp2k_check/check_cp2k_input *
!!
!!   NAME
!!     check_cp2k_input
!!
!!   FUNCTION
!!     performs further checks on an input that parsed successfully
!!
!!   NOTES
!!     at the moment does nothing
!!
!!   ARGUMENTS
!!     - input_file: the parsed input
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE check_cp2k_input(input_file,para_env,error)
    TYPE(section_vals_type), POINTER         :: input_file
    TYPE(cp_para_env_type), POINTER          :: para_env
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
       CALL handle_ext_restart(input_file,para_env,error)

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
!***************************************************************************

!!****f* input_cp2k_check/xc_functionals_expand *
!!
!!   NAME
!!     xc_functionals_expand
!!
!!   FUNCTION
!!     expand a shortcutted functional section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - functionals: the functional section to expand
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
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
!***************************************************************************

!!****f* input_cp2k_check/handle_ext_restart *
!!
!!   NAME
!!     handle_ext_restart
!!
!!   FUNCTION
!!     Replaces the requested sections in the input with those found
!!     in the external restart (EXT_RESTART%RESTART_FILE_NAME).
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - input_file: the input file to initialize
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE handle_ext_restart(input_file,para_env,error)
    TYPE(section_vals_type), POINTER         :: input_file
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'handle_ext_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(default_path_length)           :: r_file_path
    CHARACTER(default_string_length)         :: path
    INTEGER                                  :: handle, i_rep_val, &
                                                iforce_eval, myi, n_rep_val, &
                                                nforce_eval1, nforce_eval2
    LOGICAL                                  :: bsse_check, check, explicit1, &
                                                explicit2, failure, flag, &
                                                flag2, qmmm_check, &
                                                subsys_check
    REAL(KIND=dp)                            :: myt
    REAL(KIND=dp), DIMENSION(:), POINTER     :: vec, work
    TYPE(cp_error_type)                      :: sub_error
    TYPE(cp_parser_type), POINTER            :: cpparser
    TYPE(cp_unit_set_type), POINTER          :: default_units
    TYPE(section_type), POINTER              :: input_structure
    TYPE(section_vals_type), POINTER         :: r_section, rep_sections, &
                                                restart_file, section, &
                                                section1, section2, &
                                                sections1, sections2

    CALL timeset(routineN,handle)
   
    failure=.FALSE. 
    IF (.NOT. failure) THEN
       r_section => section_vals_get_subs_vals(input_file, "EXT_RESTART",error=error)
       CALL section_vals_val_get(r_section,"RESTART_FILE_NAME", c_val=r_file_path, error=error)
       
       IF (r_file_path/=" ") THEN
          CALL cp_error_init(sub_error,template_error=error,stop_level=cp_fatal_level)
          ! parse the input
          NULLIFY(input_structure, default_units)
          CALL create_cp2k_root_section(input_structure,error=sub_error)
          NULLIFY(restart_file)
          CALL section_vals_create(restart_file,input_structure, error=sub_error)
          CALL section_release(input_structure,error=sub_error)
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
                   END IF
                   
                   CALL section_vals_val_get(r_section,"RESTART_POS",l_val=flag,error=error)
                   IF (flag) THEN
                      section => section_vals_get_subs_vals(section1,"COORD",error=error)
                      CALL section_vals_set_subs_vals(section2,"COORD",section,error=error)
                   END IF
                   
                   CALL section_vals_val_get(r_section,"RESTART_RANDOMG",l_val=flag,error=error)
                   IF (flag) THEN
                      section => section_vals_get_subs_vals(section1,"RNG_INIT",error=error)
                      CALL section_vals_set_subs_vals(section2,"RNG_INIT",section,error=error)
                   END IF
                   
                   CALL section_vals_val_get(r_section,"RESTART_VEL",l_val=flag,error=error)
                   IF (flag) THEN
                      section => section_vals_get_subs_vals(section1,"VELOCITY",error=error)
                      CALL section_vals_set_subs_vals(section2,"VELOCITY",section,error=error)
                   END IF
                   
                   CALL section_vals_val_get(r_section,"RESTART_SHELL_POS",l_val=flag,error=error)
                   IF (flag) THEN
                      section => section_vals_get_subs_vals(section1,"SHELL_COORD",error=error)
                      CALL section_vals_set_subs_vals(section2,"SHELL_COORD",section,error=error)
                   END IF
  
                   CALL section_vals_val_get(r_section,"RESTART_CORE_POS",l_val=flag,error=error)
                   IF (flag) THEN
                      section => section_vals_get_subs_vals(section1,"CORE_COORD",error=error)
                      CALL section_vals_set_subs_vals(section2,"CORE_COORD",section,error=error)
                   END IF
                   
                   CALL section_vals_val_get(r_section,"RESTART_SHELL_VELOCITY",l_val=flag,error=error)
                   IF (flag) THEN
                      section => section_vals_get_subs_vals(section1,"SHELL_VELOCITY",error=error)
                      CALL section_vals_set_subs_vals(section2,"SHELL_VELOCITY",section,error=error)
                   END IF
  
                   CALL section_vals_val_get(r_section,"RESTART_CORE_VELOCITY",l_val=flag,error=error)
                   IF (flag) THEN
                      section => section_vals_get_subs_vals(section1,"CORE_VELOCITY",error=error)
                      CALL section_vals_set_subs_vals(section2,"CORE_VELOCITY",section,error=error)
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
                CALL section_vals_val_get(section1,"INITIAL_TRANSLATION_VECTOR",r_vals=vec,error=error)
                ALLOCATE(work(3))
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
             END IF
          END DO

          CALL section_vals_val_get(r_section,"RESTART_COUNTERS",l_val=flag,error=error)
          IF (flag) THEN
             IF (check_restart(input_file, restart_file, "MOTION%MD",error)) THEN
                CALL section_vals_val_get(restart_file,"MOTION%MD%STEP_START_VAL",i_val=myi,error=error)
                CALL section_vals_val_get(restart_file,"MOTION%MD%TIME_START_VAL",r_val=myt,error=error)
                CALL section_vals_val_set(input_file,"MOTION%MD%STEP_START_VAL",i_val=myi,error=error)
                CALL section_vals_val_set(input_file,"MOTION%MD%TIME_START_VAL",r_val=myt,error=error)
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
             END IF
          END IF

          CALL section_vals_val_get(r_section,"RESTART_BAND",l_val=flag,error=error)
          IF (flag.AND.check_restart(input_file, restart_file, "MOTION%BAND",error)) THEN
             rep_sections => section_vals_get_subs_vals(restart_file,"MOTION%BAND%REPLICA",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%BAND%REPLICA",rep_sections,error=error)
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
          END IF

          CALL section_vals_val_get(r_section,"RESTART_BAROSTAT_THERMOSTAT",l_val=flag,error=error)
          check = check_restart(input_file, restart_file, "MOTION%MD%BAROSTAT", error)
          CALL restart_thermostat(flag, input_file, restart_file, "MOTION%MD%BAROSTAT%THERMOSTAT", check, error)

          CALL section_vals_val_get(r_section,"RESTART_KG_COEFF_POS",l_val=flag,error=error)
          IF (flag.AND.check_restart(input_file, restart_file, "MOTION%MD%KG_COEFF",error)) THEN
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%MD%KG_COEFF%RESTART_INDEX",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%MD%KG_COEFF%RESTART_INDEX",&
                  section,error=error)
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%MD%KG_COEFF%COORD",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%MD%KG_COEFF%COORD",&
                  section,error=error)
          END IF

          CALL section_vals_val_get(r_section,"RESTART_KG_COEFF_VEL",l_val=flag,error=error)
          IF (flag.AND.check_restart(input_file, restart_file, "MOTION%MD%KG_COEFF",error)) THEN
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%MD%KG_COEFF%RESTART_INDEX",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%MD%KG_COEFF%RESTART_INDEX",&
                  section,error=error)
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%MD%KG_COEFF%VELOCITY",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%MD%KG_COEFF%VELOCITY",&
                  section,error=error)
          END IF

          CALL section_vals_val_get(r_section,"RESTART_KG_THERMOSTAT",l_val=flag,error=error)
          CALL restart_thermostat(flag, input_file, restart_file, "MOTION%MD%KG_COEFF%THERMOSTAT", error=error)

          check = check_restart(input_file, restart_file, "MOTION%MD%SHELL", error)
          IF(check) THEN
            CALL section_vals_val_get(r_section,"RESTART_SHELL_THERMOSTAT",l_val=flag,error=error)
            CALL restart_thermostat(flag, input_file, restart_file, "MOTION%MD%SHELL%THERMOSTAT", error=error)
          END IF
          
          CALL section_vals_val_get(r_section,"RESTART_THERMOSTAT",l_val=flag,error=error)
          CALL restart_thermostat(flag,input_file, restart_file, "MOTION%MD%THERMOSTAT", error=error)

          CALL section_vals_val_get(r_section,"RESTART_CONSTRAINT",l_val=flag,error=error)
          IF (flag.AND.check_restart(input_file, restart_file, "MOTION%CONSTRAINT",error)) THEN
             section => section_vals_get_subs_vals(restart_file,"MOTION%CONSTRAINT",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%CONSTRAINT",section,error=error)
          END IF

          CALL section_vals_val_get(r_section,"RESTART_METADYNAMICS",l_val=flag,error=error)
          IF (flag.AND.check_restart(input_file, restart_file, "MOTION%FREE_ENERGY%METADYN",error)) THEN
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_POS",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_POS",&
                  section,error=error)
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_AMPL",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_AMPL",&
                  section,error=error)
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_HEIGHT",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_HEIGHT",&
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
          END IF

          CALL section_vals_val_get(r_section,"RESTART_DIMER",l_val=flag,error=error)
          IF (flag.AND.check_restart(input_file, restart_file, "MOTION%GEO_OPT%TRANSITION_STATE%DIMER",error)) THEN
             section => section_vals_get_subs_vals(restart_file,&
                  "MOTION%GEO_OPT%TRANSITION_STATE%DIMER%DIMER_VECTOR",error=error)
             CALL section_vals_set_subs_vals(input_file,"MOTION%GEO_OPT%TRANSITION_STATE%DIMER%DIMER_VECTOR",&
                  section,error=error)
          END IF

          CALL section_vals_val_get(r_section,"CUSTOM_PATH",n_rep_val=n_rep_val,error=error)
          DO i_rep_val=1,n_rep_val
             CALL section_vals_val_get(r_section,"CUSTOM_PATH",i_rep_val=i_rep_val,c_val=path,error=error)
             IF (path/=" ") THEN
                section => section_vals_get_subs_vals(restart_file,path,error=error)
                CALL section_vals_set_subs_vals(input_file,path,section,error=error)
             END IF
          END DO
          CALL section_vals_val_set(r_section,"RESTART_FILE_NAME", c_val=" ",error=error)
          CALL section_vals_release(restart_file,error=error)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE handle_ext_restart

!!****f* input_cp2k_check/restart_thermostat *
!!
!!   NAME
!!     restart_thermostat
!!
!!   FUNCTION
!!     Possibly restart thermostats information
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - input_file: the input file to initialize
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] 10.2007- University of Zurich
!!
!!*** **********************************************************************
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

!!****f* input_cp2k_check/check_restart *
!!
!!   NAME
!!     check_restart
!!
!!   FUNCTION
!!     Checks if there are the proper conditions to do a restart
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - input_file: the input file to initialize
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
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

!!****f* input_cp2k_check/remove_restart_info *
!!
!!   NAME
!!     remove_restart_info
!!
!!   FUNCTION
!!     Removes section used to restart a calculation from an
!!     input file in memory
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - input_file: the input file to initialize
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE remove_restart_info(input_file,error)
    TYPE(section_vals_type), POINTER         :: input_file
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'remove_restart_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iforce_eval, &
                                                nforce_eval1
    LOGICAL                                  :: explicit1, failure
    TYPE(section_vals_type), POINTER         :: section1, section_to_delete, &
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
          section1 => section_vals_get_subs_vals3(sections1,"SUBSYS",&
               i_rep_section=iforce_eval,error=error)
          section_to_delete => section_vals_get_subs_vals(section1,"COORD",error=error)
          CALL section_vals_remove_values(section_to_delete,error)
          section_to_delete => section_vals_get_subs_vals(section1,"VELOCITY",error=error)
          CALL section_vals_remove_values(section_to_delete,error)
       END DO
       work_section => section_vals_get_subs_vals(input_file,"MOTION%MD",error=error)
       CALL section_vals_get(work_section,explicit=explicit1,error=error)
       IF (explicit1) THEN
          CALL section_vals_val_unset(input_file,"MOTION%MD%STEP_START_VAL",error=error)
          CALL section_vals_val_unset(input_file,"MOTION%MD%TIME_START_VAL",error=error)
       END IF
       work_section => section_vals_get_subs_vals(input_file,"MOTION%FREE_ENERGY%METADYN",error=error)
       CALL section_vals_get(work_section,explicit=explicit1,error=error)
       IF (explicit1) THEN
          CALL section_vals_val_unset(input_file,"MOTION%FREE_ENERGY%METADYN%STEP_START_VAL",error=error)
          CALL section_vals_val_unset(input_file,"MOTION%FREE_ENERGY%METADYN%NHILLS_START_VAL",error=error)
       END IF
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%BAND%REPLICA",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%BAND%REPLICA",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%BAROSTAT%MASS",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%BAROSTAT%VELOCITY",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%BAROSTAT%THERMOSTAT%NOSE%COORD",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%BAROSTAT%THERMOSTAT%NOSE%VELOCITY",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%BAROSTAT%THERMOSTAT%NOSE%MASS",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%BAROSTAT%THERMOSTAT%NOSE%FORCE",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%KG_COEFF%RESTART_INDEX",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%KG_COEFF%COORD",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%KG_COEFF%RESTART_INDEX",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%KG_COEFF%VELOCITY",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%KG_COEFF%THERMOSTAT%NOSE%COORD",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%KG_COEFF%THERMOSTAT%NOSE%VELOCITY",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%KG_COEFF%THERMOSTAT%NOSE%MASS",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%KG_COEFF%THERMOSTAT%NOSE%FORCE",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%THERMOSTAT%SHELL%NOSE%COORD",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%THERMOSTAT%SHELL%NOSE%VELOCITY",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%THERMOSTAT%SHELL%NOSE%MASS",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%THERMOSTAT%SHELL%NOSE%FORCE",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%THERMOSTAT%NOSE%COORD",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%THERMOSTAT%NOSE%VELOCITY",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%THERMOSTAT%NOSE%MASS",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%MD%THERMOSTAT%NOSE%FORCE",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       ! Constrained/Restrained section not reset.. to be done externally..
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_POS",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_AMPL",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%FREE_ENERGY%METADYN%SPAWNED_HILLS_HEIGHT",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%FREE_ENERGY%METADYN%EXT_LAGRANGE_SS0",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
       section_to_delete => section_vals_get_subs_vals(input_file,&
            "MOTION%FREE_ENERGY%METADYN%EXT_LAGRANGE_VVP",error=error)
       CALL section_vals_remove_values(section_to_delete,error)
    END IF
    CALL timestop(handle)
  END SUBROUTINE remove_restart_info

!!****f* input_cp2k_check/handle_defaults_restart *
!!
!!   NAME
!!     handle_defaults_restart
!!
!!   FUNCTION
!!     This subroutine controls the defaults for the restartable quantities..
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - input_file: the input file to initialize
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo - University of Zurich - 09.2007 [tlaino]
!!
!!*** **********************************************************************
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
