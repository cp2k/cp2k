!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief function that builds the resp section of the input
!> \par History
!>      02.2007 created
!> \author Joost VandeVondele
! *****************************************************************************
MODULE input_cp2k_resp
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE input_constants
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: integer_t,&
                                             real_t
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_resp'

  PUBLIC :: create_resp_section

CONTAINS

! *****************************************************************************
!> \brief Creates the RESP section
!> \param section the section to create
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE create_resp_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_resp_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="RESP",&
            description="Requests a RESP fit of charges. The Hartree potential should match "//&
                        "the one of an isolated system in those points of space that are used for fitting "//&
                        "(i.e. isolated poisson solver, and big cells). All restraints are harmonic!",&
            n_keywords=2, n_subsections=2, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="STRIDE",&
            description="specifies the stride going over the potential grid.",&
            usage="STRIDE 2",&
            default_i_val=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

       CALL keyword_create(keyword, name="INTEGER_TOTAL_CHARGE",&
            description="specifies the stride going over the potential grid.",&
            usage="INTEGER_TOTAL_CHARGE TRUE",&
            default_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

       CALL keyword_create(keyword, name="RESTRAIN_HEAVIES_TO_ZERO",&
            description="Restrain non-hydrogen atoms to zero.",&
            usage="RESTRAIN_HEAVIES_TO_ZERO FALSE",&
            default_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

       CALL keyword_create(keyword, name="RESTRAIN_HEAVIES_STRENGTH",&
            description="The strenght of the restraint on the heavy atoms",&
            usage="RESTRAIN_HEAVIES_STRENGTH 0.0001 ",&
            default_r_val=1.0E-6_dp ,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMIN",&
            description="specifies the minimum distance a fit point is away from an atom. ",&
            usage="RMIN 2.1",&
            default_r_val=cp_unit_to_cp2k(value=2.1_dp, unit_str="angstrom", error=error),&
            unit_str='angstrom', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

       CALL keyword_create(keyword, name="RMAX",&
            description="specifies the minimum distance a fit point is away from an atom. ",&
            usage="RMAX 2.5",&
            default_r_val=cp_unit_to_cp2k(value=2.5_dp, unit_str="angstrom", error=error),&
            unit_str='angstrom', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

       CALL create_constraint_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_restraint_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of basic information during the run", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section, print_key, error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_resp_section

! *****************************************************************************
!> \brief specifies constraints to be satisfied in a resp fit
!> \param section the section to create
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE create_constraint_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_constraint_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="CONSTRAINT",&
            description="specifies a linear constraint on the fitted charges."//&
                        "This can be used to give equal values to equivalent atoms." //&
                        "sum over atom_list c_i * q_i = t",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="TARGET",&
            description="the target value for the constraint",&
            usage="TARGET 0.0",&
            n_var=1,default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="ATOM_LIST",&
            description="Defines the list of atoms involved in this constraint",&
            usage="ATOM_LIST 3 4",&
            type_of_var=integer_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="ATOM_COEF",&
            description="Defines the the coefficient of the atom in this this linear constraint",&
            usage="ATOM_COEF 1.0 -1.0",&
            type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_constraint_section

! *****************************************************************************
!> \brief specifies restraints to be added to a resp fit
!> \param section the section to create
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE create_restraint_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_restraint_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="RESTRAINT",&
            description="specifies a restraint on the fitted charges."//&
                        "This can be used to restrain values to zero." //&
                        "s*(sum over atom_list c_i * q_i - t)**2",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="TARGET",&
            description="the target value for the restraint",&
            usage="TARGET 0.0",&
            n_var=1,default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STRENGTH",&
            description="the target value for the constraint",&
            usage="STRENGTH 0.001",&
            n_var=1,default_r_val=0.001_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="ATOM_LIST",&
            description="Defines the list of atoms involved in this restraint",&
            usage="ATOM_LIST 3 4",&
            type_of_var=integer_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="ATOM_COEF",&
            description="Defines the the coefficient of the atom in this this linear restraint",&
            usage="ATOM_COEF 1.0 -1.0",&
            type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_restraint_section

END MODULE input_cp2k_resp
