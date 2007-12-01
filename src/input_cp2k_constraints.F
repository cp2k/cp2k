!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_cp2k_constraints *
!!
!!   NAME
!!     input_cp2k_constraints
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     teo & fawzi
!!
!!   MODIFICATION HISTORY
!!     10.2005 split input_cp2k into smaller modules [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE input_cp2k_constraints
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE input_constants
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
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_constraints'

PUBLIC :: create_constraint_section

!***
!****************************************************************************
CONTAINS

!!****f* input_cp2k_constraints/create_constraint_section [1.0] *
!!
!!   NAME
!!     create_constraint_section
!!
!!   FUNCTION
!!     Create the constraint section. This section is useful to impose
!!     constraints
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_constraint_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_constraint_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="constraint",&
            description="Section specifying information regarding how to impose constraints"// &
            " on the system.",&
            n_keywords=0, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection, keyword, print_key)
       CALL keyword_create(keyword, name="SHAKE",&
            description="Set the tolerance for the shake constraint algorithm.",&
            usage="SHAKE <REAL>",&
            default_r_val=1.0E-6_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CONSTRAIN_INIT",&
            description="Apply constraints to the initial position and velocities."//&
            " Default is to apply constraints only after the first MD steps.",&
            usage="CONSTRAIN_INIT <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_hbonds_section(subsection, error=error)
       CALL restraint_info_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_g3x3_section(subsection, error=error)
       CALL restraint_info_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_g4x6_section(subsection, error=error)
       CALL restraint_info_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_collective_section(subsection, error=error)
       CALL restraint_info_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_fixed_atom_section(subsection, error=error)
       CALL restraint_info_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_f_a_rest_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_clv_rest_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL cp_print_key_section_create(print_key,"constraint_info",&
            description="Prints information about iterative constraints solutions",&
            print_level=high_print_level, filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"lagrange_multipliers",&
            description="Prints out the lagrange multipliers of the specified constraints during an MD.",&
            print_level=high_print_level, filename="",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF

  END SUBROUTINE create_constraint_section
!***************************************************************************

!!****f* input_cp2k_constraints/create_clv_rest_section [1.0] *
!!
!!   NAME
!!     create_clv_rest_section
!!
!!   FUNCTION
!!     Create the restart section for colvar restraints
!!     This section will be only used for restraint restarts. 
!!     Constraints are handled automatically
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino 08.2006
!!
!!*** **********************************************************************
  SUBROUTINE create_clv_rest_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_clv_rest_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY( keyword)
       CALL section_create(section,name="COLVAR_RESTART",&
            description="Specify restart position only for COLVAR restraints.",&
            n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="The restarting values for COLVAR restraints."//&
            " The order is an internal order. So if you decide to modify these values by hand"//&
            " first think what you're doing!",repeats=.TRUE.,&
            usage="<REAL>", type_of_var=real_t, n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
    END IF

  END SUBROUTINE create_clv_rest_section

!!****f* input_cp2k_constraints/create_f_a_rest_section [1.0] *
!!
!!   NAME
!!     create_f_a_rest_section
!!
!!   FUNCTION
!!     Create the restart section for fixed atoms restraints
!!     This section will be only used for restraint restarts. 
!!     Constraints are handled automatically
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino 08.2006
!!
!!*** **********************************************************************
  SUBROUTINE create_f_a_rest_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_f_a_rest_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY( keyword)
       CALL section_create(section,name="FIX_ATOM_RESTART",&
            description="Specify restart position only for FIXED_ATOMS restraints.",&
            n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="The restarting position of fixed atoms for restraints."//&
            " The order is an internal order. So if you decide to modify these values by hand"//&
            " first think what you're doing!",repeats=.TRUE.,&
            usage="<REAL> <REAL> <REAL>", type_of_var=real_t, n_var=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_f_a_rest_section

!!****f* input_cp2k_constraints/restraint_info_section [1.0] *
!!
!!   NAME
!!     restraint_info_section
!!
!!   FUNCTION
!!     Create the restraint info section in the constraint section
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino 08.2006
!!
!!*** **********************************************************************
  SUBROUTINE restraint_info_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'restraint_info_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(subsection, keyword)
       CALL section_create(subsection,name="RESTRAINT",&
            description="Activate and specify information on restraint instead of constraint",&
            n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="K",&
            description="Specifies the force constant for the harmonic restraint",&
            usage="K {real}",&
            type_of_var=real_t, default_r_val=0.0_dp, unit_str="internal_cp2k",error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)       
    END IF

  END SUBROUTINE restraint_info_section

!!****f* input_cp2k_constraints/create_collective_section [1.0] *
!!
!!   NAME
!!     create_collecive_section
!!
!!   FUNCTION
!!     Create the constraint section for collective constraints
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Joost VandeVondele [01.2006]
!!
!!*** **********************************************************************
  SUBROUTINE create_collective_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_collective_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="COLLECTIVE",&
            description="Used to constraint collective (general) degrees of freedom, "//&
            "writing langrangian multipliers to file.",&
            n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="COLVAR",&
            description="Specifies the index (in input file order) of the type of colvar to constrain.",&
            usage="COLVAR {int}",&
            type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOLECULE",&
            description="Specifies the index of the molecule kind (in input file order)"//&
            "on which the constraint will be applied."//&
            " MOLECULE and MOLNAME keyword exclude themself mutually.",&
            usage="MOLECULE {integer}", required=.TRUE., n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOLNAME",&
            variants=(/"SEGNAME"/),&
            description="Specifies the name of the molecule on which the constraint will be applied.",&
            usage="MOLNAME {character}", required=.TRUE., n_var=1, type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INTERMOLECULAR",&
            description="Specify if the constraint/restraint is intermolecular.",&
            usage="INTERMOLECULAR <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       

       CALL keyword_create(keyword, name="TARGET",&
            description="Specifies the target value of the constrained collective"//&
            " variable (units depend on the colvar).",&
            usage="TARGET {real}",&
            type_of_var=real_t, unit_str="internal_cp2k",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXCLUDE_QM",&
            description="Does not apply the constraint to the QM region within a QM/MM calculation",&
            usage="EXCLUDE_QM <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXCLUDE_MM",&
            description="Does not apply the constraint to the MM region within a QM/MM calculation",&
            usage="EXCLUDE_MM <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_collective_section

!***************************************************************************
!!****f* input_cp2k_constraints/create_fixed_atom_section [1.0] *
!!
!!   NAME
!!     create_fixed_atom_section
!!
!!   FUNCTION
!!     Create the constraint section that fixes atoms
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_fixed_atom_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_fixed_atom_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="fixed_atoms",&
            description="This section is used to constraint the overall atomic position (X,Y,Z).",&
            n_keywords=3, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       ! Integer
       CALL keyword_create(keyword, name="LIST",&
            description="Specifies a list of atoms to freeze.",&
            usage="LIST {integer} {integer} .. {integer}", required=.TRUE.,repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RANGE",&
            description="Specifies the starting and the ending index of atoms range you wish to freeze.",&
            usage="RANGE {integer} {integer}", required=.TRUE.,repeats=.TRUE.,&
            n_var=2, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOLNAME",&
            variants=(/"SEGNAME"/),&
            description="Specifies the name of the molecule to fix",&
            usage="MOLNAME WAT MEOH", required=.TRUE., repeats=.TRUE.,&
            n_var=-1,type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MM_SUBSYS",&
            variants=(/"PROTEIN"/),&
            description="In a QM/MM run all  MM atoms are fixed according to the argument.",&
            usage="MM_SUBSYS (NONE|ATOMIC|MOLECULAR)",&
            enum_c_vals=s2a( "NONE","ATOMIC","MOLECULAR"),&
            enum_i_vals=(/do_constr_none,do_constr_atomic,do_constr_molec/),&
            enum_desc=s2a("fix nothing",&
            "only the MM atoms itself",&
            "the full molecule/residue that contains a MM atom (i.e. some QM atoms might be fixed as well)"),&
            default_i_val=do_constr_none,repeats=.FALSE.,required=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="QM_SUBSYS",&
            description="In a QM/MM run all QM atoms are fixed according to the argument.",&
            usage="QM_SUBSYS (NONE|ATOMIC|MOLECULAR)",&
            enum_c_vals=s2a( "NONE","ATOMIC","MOLECULAR"),&
            enum_desc=s2a("fix nothing",&
            "only the QM atoms itself",&
            "the full molecule/residue that contains a QM atom (i.e. some MM atoms might be fixed as well)"),&
            enum_i_vals=(/do_constr_none,do_constr_atomic,do_constr_molec/),&
            default_i_val=do_constr_none,repeats=.FALSE.,required=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXCLUDE_QM",&
            description="Does not apply the constraint to the QM region within a QM/MM calculation."//&
            " This keyword is active only together with MOLNAME",&
            usage="EXCLUDE_QM <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXCLUDE_MM",&
            description="Does not apply the constraint to the MM region within a QM/MM calculation."//&
            " This keyword is active only together with MOLNAME",&
            usage="EXCLUDE_MM <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_fixed_atom_section

!***************************************************************************
!!****f* input_cp2k_constraints/create_g3x3_section [1.0] *
!!
!!   NAME
!!     create_g3x3_section
!!
!!   FUNCTION
!!     Create the constraint section specialized on g3x3 constraints
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_g3x3_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_g3x3_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="g3x3",&
            description="This section is used to set 3x3 (3 atoms and 3 distances) constraints.",&
            n_keywords=3, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       ! Integer
       CALL keyword_create(keyword, name="MOLECULE",&
            variants=(/"MOL"/),&
            description="Specifies the molecule kind number on which constraint will be applied."//&
            " MOLECULE and MOLNAME keyword exclude themself mutually.",&
            usage="MOL {integer}", required=.TRUE., n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOLNAME",&
            variants=(/"SEGNAME"/),&
            description="Specifies the name of the molecule on which the constraint will be applied.",&
            usage="MOLNAME {character}", required=.TRUE., n_var=1, type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INTERMOLECULAR",&
            description="Specify if the constraint/restraint is intermolecular.",&
            usage="INTERMOLECULAR <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Atoms' index on which apply the constraint", usage="ATOMS 1 3 6",&
            n_var=-1,type_of_var=integer_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Real
       CALL keyword_create(keyword, name="DISTANCES",&
            description="The constrained distances' values.",&
            usage="DISTANCES {real} {real} {real}", type_of_var=real_t,&
            unit_str="internal_cp2k",n_var=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)              

       ! Logical
       CALL keyword_create(keyword, name="EXCLUDE_QM",&
            description="Does not apply the constraint to the QM region within a QM/MM calculation",&
            usage="EXCLUDE_QM <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXCLUDE_MM",&
            description="Does not apply the constraint to the MM region within a QM/MM calculation",&
            usage="EXCLUDE_MM <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF

  END SUBROUTINE create_g3x3_section

!***************************************************************************
!!****f* input_cp2k_constraints/create_hbonds_section [1.0] *
!!
!!   NAME
!!     create_hbonds_section
!!
!!   FUNCTION
!!     Create the constraint section specialized on H BONDS constraints
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_hbonds_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_hbonds_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="HBONDS",&
            description="This section is used to set bonds constraints involving Hydrogen atoms",&
            n_keywords=3, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       ! Character
       CALL keyword_create(keyword, name="ATOM_TYPE",&
            description="Defines the atoms' type forming a bond with an hydrogen. If not specified "//&
            " the default bond value of the first molecule is used as constraint target",&
            usage="ATOMS <CHARACTER>",&
            n_var=-1,type_of_var=char_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOLECULE",&
            description="Specifies the indexes of the molecule kind (in input file order)"//&
            "on which the constraint will be applied."//&
            " MOLECULE and MOLNAME keyword exclude themself mutually.",&
            usage="MOLECULE {integer} .. {integer} ", required=.TRUE., n_var=-1,&
            type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOLNAME",&
            variants=(/"SEGNAME"/),&
            description="Specifies the names of the molecule on which the constraint will be applied.",&
            usage="MOLNAME {character} .. {character} ", required=.TRUE., n_var=-1,&
            type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXCLUDE_QM",&
            description="Does not shake HBONDS in the QM region within a QM/MM calculation",&
            usage="EXCLUDE_QM <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXCLUDE_MM",&
            description="Does not shake HBONDS in the MM region within a QM/MM calculation",&
            usage="EXCLUDE_MM <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Real
       CALL keyword_create(keyword, name="TARGETS",&
            description="The constrained distances' values  for the types defines in ATOM_TYPE.",&
            usage="TARGETS {real} {real} {real}", type_of_var=real_t, n_var=-1,&
            unit_str="internal_cp2k", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)              

    END IF

  END SUBROUTINE create_hbonds_section


!***************************************************************************
!!****f* input_cp2k_constraints/create_g4x6_section [1.0] *
!!
!!   NAME
!!     create_g4x6_section
!!
!!   FUNCTION
!!     Create the constraint section specialized on g4x6 constraints
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_g4x6_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_g4x6_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="g4x6",&
            description="This section is used to set 4x6 (4 atoms and 6 distances) constraints.",&
            n_keywords=3, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       ! Integer
       CALL keyword_create(keyword, name="MOLECULE",&
            variants=(/"MOL"/),&
            description="Specifies the molecule number on which constraint will be applied."//&
            " MOLECULE and MOLNAME keyword exclude themself mutually.",&
            usage="MOL {integer}", required=.TRUE., n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOLNAME",&
            variants=(/"SEGNAME"/),&
            description="Specifies the name of the molecule on which the constraint will be applied.",&
            usage="MOLNAME {character}", required=.TRUE., n_var=1, type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INTERMOLECULAR",&
            description="Specify if the constraint/restraint is intermolecular.",&
            usage="INTERMOLECULAR <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Atoms' index on which apply the constraint", usage="ATOMS 1 3 6 4",&
            n_var=4,type_of_var=integer_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Real
       CALL keyword_create(keyword, name="DISTANCES",&
            description="The constrained distances' values.",&
            usage="DISTANCES {real} {real} {real} {real} {real} {real}",&
            type_of_var=real_t, n_var=6,unit_str="internal_cp2k",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)              

       ! Logical
       CALL keyword_create(keyword, name="EXCLUDE_QM",&
            description="Does not apply the constraint to the QM region within a QM/MM calculation",&
            usage="EXCLUDE_QM <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXCLUDE_MM",&
            description="Does not apply the constraint to the MM region within a QM/MM calculation",&
            usage="EXCLUDE_MM <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF

  END SUBROUTINE create_g4x6_section

END MODULE input_cp2k_constraints
