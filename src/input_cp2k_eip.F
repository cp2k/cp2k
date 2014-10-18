!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Creates the EIP section of the input
!> \par History
!>      03.2006 created
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
MODULE input_cp2k_eip
  USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                             high_print_level,&
                                             medium_print_level
  USE f77_blas
  USE input_constants,                 ONLY: use_bazant_eip,&
                                             use_lenosky_eip
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: enum_t
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_eip'

  PUBLIC :: create_eip_section

CONTAINS

! *****************************************************************************
!> \brief Create the input section for EIP
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2006 created
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE create_eip_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_eip_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

!   ------------------------------------------------------------------------

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="EIP", &
            description="This section contains all information to run an "//&
            "Empirical Interatomic Potential (EIP) calculation.", &
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE., &
            error=error)

       NULLIFY(subsection, keyword)

       CALL keyword_create(keyword, name="EIP_MODEL", &
            description="Selects the empirical interaction potential model", &
            usage="EIP_MODEL BAZANT", required=.TRUE., type_of_var=enum_t, &
            n_var=1, repeats=.FALSE., variants=(/"EIP-MODEL"/), &
            enum_c_vals=s2a("BAZANT", "EDIP", "LENOSKY"), &
            enum_i_vals=(/use_bazant_eip, use_bazant_eip, use_lenosky_eip/), &
            default_i_val=use_lenosky_eip, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL create_eip_print_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error)
       CALL section_release(subsection, error=error)

    END IF

  END SUBROUTINE create_eip_section

! *****************************************************************************
!> \brief Creates the print section for the eip subsection
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2006 created
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE create_eip_print_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_eip_print_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: print_key

!   ------------------------------------------------------------------------

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section, name="PRINT", &
            description="Section of possible print options in EIP code.", &
            n_keywords=0, n_subsections=6, repeats=.FALSE., required=.TRUE., &
            error=error)

       NULLIFY(print_key)

       CALL cp_print_key_section_create(print_key, "ENERGIES", &
            description="Controls the printing of the EIP energies.", &
            print_level=medium_print_level, filename="__STD_OUT__", &
            error=error)
       CALL section_add_subsection(section, print_key, error=error)
       CALL section_release(print_key, error=error)

       CALL cp_print_key_section_create(print_key, "ENERGIES_VAR", &
            description="Controls the printing of the variance of the EIP energies.", &
            print_level=high_print_level, filename="__STD_OUT__", &
            error=error)
       CALL section_add_subsection(section, print_key, error=error)
       CALL section_release(print_key, error=error)

       CALL cp_print_key_section_create(print_key, "FORCES", &
            description="Controls the printing of the EIP forces.", &
            print_level=medium_print_level, filename="__STD_OUT__", &
            error=error)
       CALL section_add_subsection(section, print_key, error=error)
       CALL section_release(print_key, error=error)

       CALL cp_print_key_section_create(print_key, "COORD_AVG", &
            description="Controls the printing of the average coordination number.", &
            print_level=high_print_level, filename="__STD_OUT__", &
            error=error)
       CALL section_add_subsection(section, print_key, error=error)
       CALL section_release(print_key, error=error)

       CALL cp_print_key_section_create(print_key, "COORD_VAR", &
            description="Controls the printing of the variance of the coordination number.", &
            print_level=high_print_level, filename="__STD_OUT__", &
            error=error)
       CALL section_add_subsection(section, print_key, error=error)
       CALL section_release(print_key, error=error)

       CALL cp_print_key_section_create(print_key, "COUNT", &
            description="Controls the printing of the number of function calls.", &
            print_level=high_print_level, filename="__STD_OUT__", &
            error=error)
       CALL section_add_subsection(section, print_key, error=error)
       CALL section_release(print_key, error=error)

    END IF
  END SUBROUTINE create_eip_print_section

END MODULE input_cp2k_eip
