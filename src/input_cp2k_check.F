!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
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
  USE kinds, ONLY: dp, default_string_length
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
  USE input_section_types, ONLY: section_release,&
       section_type,section_vals_create,&
       section_vals_type, section_vals_get_subs_vals,&
       section_vals_val_get, section_vals_get, section_vals_set_subs_vals,&
       section_vals_val_set, section_vals_release
  USE input_cp2k,                      ONLY: create_cp2k_root_section
  USE input_constants, ONLY: xc_funct_no_shortcut, xc_funct_blyp,&
       xc_funct_pade, xc_funct_pbe, xc_funct_tpss, xc_funct_olyp, xc_funct_hcth120
  USE cp_parser,                       ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE input_parsing,                   ONLY: section_vals_parse
  USE cp_units, ONLY: cp_unit_set_type, cp_unit_set_create, cp_unit_set_release
#include "cp_common_uses.h"

  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='input_cp2k_check'

  public :: check_cp2k_input, xc_functionals_expand

!!***
!****************************************************************************
contains


!!****f* input_cp2k/check_cp2k_input *
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
    TYPE(section_vals_type), pointer :: input_file
    TYPE(cp_para_env_type), pointer :: para_env
    type(cp_error_type), optional, intent(inout) :: error

    LOGICAL :: failure, explicit
    character(len=*), parameter :: routineN='check_cp2k_input',&
         routineP=moduleN//':'//routineN
    TYPE(section_vals_type), pointer :: section

    failure=.false.

    CPPrecondition(ASSOCIATED(input_file),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(input_file%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       ! ext_restart
       CALL handle_ext_restart(input_file,para_env,error)

       ! xc: expand and fix default for tddfpt
       section => section_vals_get_subs_vals(input_file,&
            "FORCE_EVAL%DFT%XC%XC_FUNCTIONAL",error=error)
       CALL xc_functionals_expand(section,error=error)
       section => section_vals_get_subs_vals(input_file,&
            "FORCE_EVAL%DFT%KG%XC_FUNCTIONAL",error=error)
       CALL xc_functionals_expand(section,error=error)
       section => section_vals_get_subs_vals(input_file,&
            "FORCE_EVAL%DFT%TDDFPT%XC%XC_FUNCTIONAL",error=error)
       CALL section_vals_get(section,explicit=explicit,error=error)
       IF (explicit) THEN
          CALL xc_functionals_expand(section,error=error)
       ELSE
          section => section_vals_get_subs_vals(input_file,&
               "FORCE_EVAL%DFT%XC%XC_FUNCTIONAL",error=error)
          CALL section_vals_set_subs_vals(input_file,&
               "FORCE_EVAL%DFT%TDDFPT%XC%XC_FUNCTIONAL",&
               section,error=error)
       END IF

    END IF
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
SUBROUTINE xc_functionals_expand(functionals,error)
  TYPE(section_vals_type), pointer :: functionals
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='xc_functionals_expand',&
        routineP=moduleN//':'//routineN
  integer :: shortcut
  failure=.false.
  
  IF (.NOT. failure) THEN
     CALL section_vals_val_get(functionals,"_SECTION_PARAMETERS_",&
          i_val=shortcut,error=error)
     SELECT CASE(shortcut)
     CASE(xc_funct_no_shortcut)
        ! nothing to expand
     CASE(xc_funct_blyp)
        CALL section_vals_val_set(functionals,"BECKE88%_SECTION_PARAMETERS_",&
             l_val=.TRUE.,error=error)
        CALL section_vals_val_set(functionals,"LYP%_SECTION_PARAMETERS_",&
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
!!     in the external restart (EXT_RESTART%EXTERNAL_FILE).
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
  TYPE(section_vals_type), pointer :: input_file
  TYPE(cp_para_env_type), pointer :: para_env
  type(cp_error_type), optional, intent(inout) :: error

  LOGICAL :: failure, flag, explicit1, explicit2
  character(len=*), parameter :: routineN='handle_ext_restart',&
       routineP=moduleN//':'//routineN
  TYPE(section_vals_type), POINTER :: restart_input, r_section, section, work_section
  CHARACTER(default_string_length) :: r_file_path, path
  TYPE(section_type), POINTER              :: input_structure
  TYPE(cp_unit_set_type), POINTER :: default_units
  TYPE(cp_parser_type), POINTER            :: cpparser
  INTEGER :: n_rep_val, i_rep_val
  REAL(KIND=dp) :: myt
  INTEGER       :: myi

  failure=.false.
  NULLIFY(work_section)
  IF (.NOT. failure) THEN
     r_section => section_vals_get_subs_vals(input_file,&
          "EXT_RESTART",error=error)
     CALL section_vals_val_get(r_section,"EXTERNAL_FILE", &
          c_val=r_file_path, error=error)

     IF (r_file_path/=" ") THEN
        ! parse the input
        NULLIFY(input_structure, default_units)
        CALL create_cp2k_root_section(input_structure,error=error)
        NULLIFY(restart_input)
        CALL section_vals_create(restart_input,input_structure, error=error)
        CALL section_release(input_structure,error=error)
        NULLIFY(cpparser)
        CALL parser_create(cpparser,file_name=r_file_path,para_env=para_env,&
             error=error)
        CALL cp_unit_set_create(default_units, "OUTPUT",error=error)
        CALL section_vals_parse(restart_input,cpparser,root_section=.FALSE.,&
             default_units=default_units,error=error)
        CALL cp_unit_set_release(default_units,error=error)
        CALL parser_release(cpparser,error=error)

        CALL section_vals_val_get(r_section,"RESTART_CELL",l_val=flag,&
             error=error)
        IF (flag) THEN
           section => section_vals_get_subs_vals(restart_input,&
                "FORCE_EVAL%SUBSYS%CELL",error=error)
           CALL section_vals_set_subs_vals(input_file,"FORCE_EVAL%SUBSYS%CELL",&
                section,error=error)
        END IF
        CALL section_vals_val_get(r_section,"RESTART_POS",l_val=flag,&
             error=error)
        IF (flag) THEN
           section => section_vals_get_subs_vals(restart_input,&
                "FORCE_EVAL%SUBSYS%COORD",error=error)
           CALL section_vals_set_subs_vals(input_file,"FORCE_EVAL%SUBSYS%COORD",&
                section,error=error)
        END IF
        CALL section_vals_val_get(r_section,"RESTART_RANDOMG",l_val=flag,&
             error=error)
        IF (flag) THEN
           section => section_vals_get_subs_vals(restart_input,&
                "FORCE_EVAL%SUBSYS%RNG_INIT",error=error)
           CALL section_vals_set_subs_vals(input_file,"FORCE_EVAL%SUBSYS%RNG_INIT",&
                section,error=error)
        END IF
        CALL section_vals_val_get(r_section,"RESTART_VEL",l_val=flag,&
             error=error)
        IF (flag) THEN
           section => section_vals_get_subs_vals(restart_input,&
                "FORCE_EVAL%SUBSYS%VELOCITY",error=error)
           CALL section_vals_set_subs_vals(input_file,"FORCE_EVAL%SUBSYS%VELOCITY",&
                section,error=error)
        END IF
        CALL section_vals_val_get(r_section,"RESTART_COUNTERS",l_val=flag,&
             error=error)
        IF (flag) THEN
           work_section => section_vals_get_subs_vals(input_file,"MOTION%MD",error=error)
           CALL section_vals_get(work_section,explicit=explicit1,error=error)
           work_section => section_vals_get_subs_vals(restart_input,"MOTION%MD",error=error)
           CALL section_vals_get(work_section,explicit=explicit2,error=error)
           IF (explicit1.AND.explicit2) THEN
              CALL section_vals_val_get(restart_input,"MOTION%MD%STEP_START_VAL",i_val=myi,error=error)
              CALL section_vals_val_get(restart_input,"MOTION%MD%TIME_START_VAL",r_val=myt,error=error)
              CALL section_vals_val_set(input_file,"MOTION%MD%STEP_START_VAL",i_val=myi,error=error)
              CALL section_vals_val_set(input_file,"MOTION%MD%TIME_START_VAL",r_val=myt,error=error)
           END IF
           !
           work_section => section_vals_get_subs_vals(input_file,"MOTION%MD%HMC",error=error)
           CALL section_vals_get(work_section,explicit=explicit1,error=error)
           work_section => section_vals_get_subs_vals(restart_input,"MOTION%MD%HMC",error=error)
           CALL section_vals_get(work_section,explicit=explicit2,error=error)
           IF (explicit1.AND.explicit2) THEN
              CALL section_vals_val_get(restart_input,"MOTION%MD%HMC%MD_APPROACH_START_VAL",i_val=myi,error=error)
              CALL section_vals_val_set(input_file,"MOTION%MD%HMC%MD_APPROACH_START_VAL",i_val=myi,error=error)
              CALL section_vals_val_get(restart_input,"MOTION%MD%HMC%ACCEPTED_MD_START_VAL",i_val=myi,error=error)
              CALL section_vals_val_set(input_file,"MOTION%MD%HMC%ACCEPTED_MD_START_VAL",i_val=myi,error=error)
           END IF
        END IF
        CALL section_vals_val_get(r_section,"RESTART_BAROSTAT",l_val=flag,&
             error=error)
        IF (flag) THEN
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%BAROSTAT%MASS",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%BAROSTAT%MASS",&
                section,error=error)
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%BAROSTAT%VELOCITY",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%BAROSTAT%VELOCITY",&
                section,error=error)
        END IF
        CALL section_vals_val_get(r_section,"RESTART_BAROSTAT_NOSE",l_val=flag,&
             error=error)
        IF (flag) THEN
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%BAROSTAT%NOSE%COORD",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%BAROSTAT%NOSE%COORD",&
                section,error=error)
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%BAROSTAT%NOSE%VELOCITY",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%BAROSTAT%NOSE%VELOCITY",&
                section,error=error)
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%BAROSTAT%NOSE%MASS",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%BAROSTAT%NOSE%MASS",&
                section,error=error)
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%BAROSTAT%NOSE%FORCE",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%BAROSTAT%NOSE%FORCE",&
                section,error=error)
        END IF
        CALL section_vals_val_get(r_section,"RESTART_KG_COEFF_POS",l_val=flag,&
             error=error)
        IF (flag) THEN
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%KG_COEFF%RESTART_INDEX",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%KG_COEFF%RESTART_INDEX",&
                section,error=error)
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%KG_COEFF%COORD",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%KG_COEFF%COORD",&
                section,error=error)
        END IF
        CALL section_vals_val_get(r_section,"RESTART_KG_COEFF_VEL",l_val=flag,&
             error=error)
        IF (flag) THEN
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%KG_COEFF%RESTART_INDEX",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%KG_COEFF%RESTART_INDEX",&
                section,error=error)
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%KG_COEFF%VELOCITY",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%KG_COEFF%VELOCITY",&
                section,error=error)
        END IF
        CALL section_vals_val_get(r_section,"RESTART_KG_NOSE",l_val=flag,&
             error=error)
        IF (flag) THEN
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%KG_COEFF%NOSE%COORD",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%KG_COEFF%NOSE%COORD",&
                section,error=error)
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%KG_COEFF%NOSE%VELOCITY",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%KG_COEFF%NOSE%VELOCITY",&
                section,error=error)
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%KG_COEFF%NOSE%MASS",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%KG_COEFF%NOSE%MASS",&
                section,error=error)
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%KG_COEFF%NOSE%FORCE",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%KG_COEFF%NOSE%FORCE",&
                section,error=error)
        END IF
        CALL section_vals_val_get(r_section,"RESTART_NOSE",l_val=flag,&
             error=error)
        IF (flag) THEN
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%NOSE%COORD",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%NOSE%COORD",&
                section,error=error)
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%NOSE%VELOCITY",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%NOSE%VELOCITY",&
                section,error=error)
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%NOSE%MASS",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%NOSE%MASS",&
                section,error=error)
           section => section_vals_get_subs_vals(restart_input,&
                "MOTION%MD%NOSE%FORCE",error=error)
           CALL section_vals_set_subs_vals(input_file,"MOTION%MD%NOSE%FORCE",&
                section,error=error)
        END IF
        CALL section_vals_val_get(r_section,"CUSTOM_PATH",n_rep_val=n_rep_val,&
             error=error)
        DO i_rep_val=1,n_rep_val
           CALL section_vals_val_get(r_section,"CUSTOM_PATH",&
                i_rep_val=i_rep_val,c_val=path,&
                error=error)
           IF (path/=" ") THEN
              section => section_vals_get_subs_vals(restart_input,&
                   path,error=error)
              CALL section_vals_set_subs_vals(input_file,path,&
                   section,error=error)
           END IF
        END DO
        CALL section_vals_val_set(r_section,"EXTERNAL_FILE",&
             c_val=" ",error=error)
        CALL section_vals_release(restart_input,error=error)
     END IF
  END IF
END SUBROUTINE handle_ext_restart
!***************************************************************************

END MODULE input_cp2k_check
