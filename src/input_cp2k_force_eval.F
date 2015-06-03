!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief builds the input structure for the FORCE_EVAL section of cp2k
!> \par History
!>      06.2004 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE input_cp2k_force_eval
  USE cp_output_handling,              ONLY: add_last_numeric,&
                                             cp_print_key_section_create,&
                                             high_print_level,&
                                             low_print_level,&
                                             medium_print_level
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE input_constants,                 ONLY: &
       do_eip, do_fist, do_mixed, do_qmmm, do_qs, do_stress_analytical, &
       do_stress_diagonal_anal, do_stress_diagonal_numer, do_stress_none, &
       do_stress_numerical, numerical
  USE input_cp2k_dft,                  ONLY: create_bsse_section,&
                                             create_dft_section
  USE input_cp2k_eip,                  ONLY: create_eip_section
  USE input_cp2k_mixed,                ONLY: create_mix_section
  USE input_cp2k_mm,                   ONLY: create_mm_section
  USE input_cp2k_properties_dft,       ONLY: create_properties_section
  USE input_cp2k_qmmm,                 ONLY: create_qmmm_section
  USE input_cp2k_subsys,               ONLY: create_subsys_section
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: char_t,&
                                             integer_t,&
                                             lchar_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_force_eval'

  PUBLIC :: create_force_eval_section

CONTAINS

! *****************************************************************************
!> \brief creates the force_eval section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE create_force_eval_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_force_eval_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="force_eval",&
         description="parameters needed to calculate energy and forces and"//&
         " describe the system you want to analyze.",&
         n_keywords=1, n_subsections=9, repeats=.TRUE., required=.TRUE.,&
         error=error)

    NULLIFY(subsection)
    NULLIFY(keyword)
    CALL keyword_create(keyword, name="METHOD",&
         description="Which method should be used to compute forces",&
         usage="METHOD <STRING>",&
         enum_c_vals=s2a("QS",&
                         "FIST",&
                         "QMMM",&
                         "EIP",&
                         "QUICKSTEP",&
                         "MIXED"),&
         enum_desc=s2a("Electronic structure methods (DFT, ...)",&
                       "Molecular Mechanics",&
                       "Hybrid quantum classical",&
                       "Empirical Interatomic Potential",&
                       "Electronic structure methods (DFT, ...)",&
                       "Use a combination of two of the above"),&
         enum_i_vals=(/do_qs, do_fist, do_qmmm, do_eip, do_qs, do_mixed/),&
         default_i_val=do_qs, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="STRESS_TENSOR",&
         description="Controls the calculation of the stress tensor. The combinations defined below"//&
         " are not implemented for all methods.",&
         usage="stress_tensor (NONE|ANALYTICAL|NUMERICAL|DIAGONAL_ANA|DIAGONAL_NUM)", &
         default_i_val=do_stress_none,&
         enum_c_vals=s2a("NONE","ANALYTICAL","NUMERICAL","DIAGONAL_ANALYTICAL","DIAGONAL_NUMERICAL"),&
         enum_i_vals=(/do_stress_none,do_stress_analytical,do_stress_numerical,&
                       do_stress_diagonal_anal,do_stress_diagonal_numer/),&
         enum_desc=s2a("Do not compute stress tensor",&
                       "Compute the stress tensor analytically (if available).",&
                       "Compute the stress tensor numerically.",&
                       "Compute the diagonal part only of the stress tensor analytically (if available).",&
                       "Compute the diagonal part only of the stress tensor numerically"),&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL create_ext_pot_section(subsection, error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_rescale_force_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_mix_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_dft_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_mm_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_qmmm_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_eip_section(subsection, error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection, error=error)

    CALL create_bsse_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_subsys_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_properties_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_f_env_print_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

  END SUBROUTINE create_force_eval_section

! *****************************************************************************
!> \brief Creates the section for applying an external potential
!> \param section ...
!> \param error ...
!> \date 03.2008
!> \author teo
! *****************************************************************************
  SUBROUTINE create_ext_pot_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ext_pot_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="EXTERNAL_POTENTIAL",&
            description="Section controlling the presence of an external potential dependent "//&
            " on the atomic positions (X,Y,Z)",&
            n_keywords=7, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS_LIST",&
            description="Specifies the atoms on which the external potential will act",&
            usage="ATOMS_LIST {INT} {INT} ..",required=.FALSE.,repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FUNCTION",&
            description="Specifies the functional form in mathematical notation. Variables must be the atomic "//&
            "coordinates (X,Y,Z).",usage="FUNCTION  X^2+Y^2+Z^2+LOG(ABS(X+Y))", required=.TRUE.,&
            type_of_var=lchar_t, n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARAMETERS",&
            description="Defines the parameters of the functional form",&
            usage="PARAMETERS a b D", required=.TRUE., type_of_var=char_t,&
            n_var=-1, repeats=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VALUES",&
            description="Defines the values of  parameter of the functional form",&
            usage="VALUES ", required=.TRUE., type_of_var=real_t,&
            n_var=-1, repeats=.TRUE., unit_str="internal_cp2k", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="UNITS",&
            description="Optionally, allows to define valid CP2K unit strings for each parameter value. "//&
                        "It is assumed that the corresponding parameter value is specified in this unit.",&
            usage="UNITS angstrom eV*angstrom^-1 angstrom^1 K", required=.FALSE., type_of_var=char_t,&
            n_var=-1, repeats=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DX",&
            description="Parameter used for computing the derivative with the Ridders' method.",&
            usage="DX <REAL>", default_r_val=0.1_dp, unit_str="bohr", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ERROR_LIMIT",&
            description="Checks that the error in computing the derivative is not larger than "//&
            "the value set; in case error is larger a warning message is printed.",&
            usage="ERROR_LIMIT <REAL>", default_r_val=1.0E-12_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


    END IF
  END SUBROUTINE create_ext_pot_section

! *****************************************************************************
!> \brief Creates the section controlling the rescaling of forces
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_rescale_force_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_rescale_force_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="RESCALE_FORCES",&
            description="Section controlling the rescaling of forces. Useful when"//&
            " starting from quite bad geometries with unphysically large forces.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="MAX_FORCE",&
            description="Specify the Maximum Values of the force. If the force"//&
            " of one atom exceed this value it's rescaled to the MAX_FORCE"//&
            " value.",&
            default_r_val=cp_unit_to_cp2k(value=50.0_dp,&
                                          unit_str="kcalmol*angstrom^-1",&
                                          error=error),&
             unit_str="hartree*bohr^-1",&
             error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_rescale_force_section

! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE create_f_env_print_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_f_env_print_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    NULLIFY (keyword)
    NULLIFY (print_key)

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       CALL section_create(section,name="PRINT",&
            description="Properties that you want to output and that are common to all methods",&
            n_keywords=0, n_subsections=5, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of basic information generated by force_eval", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"FORCES",&
            description="Controls the printing of the forces after each force evaluation",&
            print_level=high_print_level, filename="__STD_OUT__",error=error)
       CALL keyword_create(keyword=keyword,&
                           name="NDIGITS",&
                           description="Specifies the number of digits used "//&
                                       "for the printing of the forces",&
                           usage="NDIGITS 6",&
                           default_i_val=8,&
                           required=.FALSE.,&
                           repeats=.FALSE.,&
                           error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"GRID_INFORMATION",&
            description="Controls the printing of information regarding the PW and RS grid structures.",&
            print_level=medium_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"TOTAL_NUMBERS",&
            description="Controls the printing of the total number of atoms, kinds,...",&
            print_level=medium_print_level, filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DISTRIBUTION",&
            description="Controls the printing of the distribution of molecules, atoms, ...",&
            print_level=medium_print_level, filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DISTRIBUTION2D",&
            description="Controls the printing of the distribution of matrix blocks,...",&
            print_level=high_print_level, filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DISTRIBUTION1D",&
            description="Each node prints out its distribution info ...",&
            print_level=high_print_level, filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"STRESS_TENSOR",&
            description="Controls the printing of the stress tensor",&
            print_level=high_print_level, filename="__STD_OUT__",error=error)
       CALL keyword_create(keyword=keyword,&
                           name="NDIGITS",&
                           description="Specifies the number of digits used "//&
                                       "for the printing of the stress tensor",&
                           usage="NDIGITS 6",&
                           default_i_val=8,&
                           required=.FALSE.,&
                           repeats=.FALSE.,&
                           error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF

  END SUBROUTINE create_f_env_print_section

END MODULE input_cp2k_force_eval
