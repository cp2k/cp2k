!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_cp2k_mm *
!!
!!   NAME
!!     input_cp2k_mm
!!
!!   FUNCTION
!!     creates the mm section of the input
!!
!!   NOTES
!!     moved out of input_cp2k
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2004 created
!!
!!   SOURCE
!****************************************************************************
MODULE input_cp2k_mm
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE input_constants
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_poisson,                   ONLY: create_poisson_section
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
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_mm'

  PUBLIC :: create_mm_section
  PUBLIC :: create_NONBONDED14_section, create_LJ_section, create_Williams_section, create_Goodwin_section
!***
!****************************************************************************
CONTAINS

!***************************************************************************

!!****f* input_cp2k_mm/create_mm_section [1.0] *
!!
!!   NAME
!!     create_mm_section
!!
!!   FUNCTION
!!     Create the input section for FIST.. Come on.. Let's get woohooo
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
  SUBROUTINE create_mm_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mm_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="mm",&
            description="This section contains all information to run a MM calculation.",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)

       CALL create_forcefield_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
       
       CALL create_poisson_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_print_mm_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_mm_section

!***************************************************************************
!!****f* input_cp2k_mm/create_print_mm_section [1.0] *
!!
!!   NAME
!!     create_print_mm_section
!!
!!   FUNCTION
!!     Create the print mm section
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
  SUBROUTINE create_print_mm_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_print_mm_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options in MM code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(print_key)

       CALL cp_print_key_section_create(print_key,"DERIVATIVES",&
            description="Controls the printing of derivatives.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DIPOLE",&
            description="Controls the printing of dipole information. Computed with"//&
            " the Berry phase formalism. Results in atomic Units.", &
            print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"NEIGHBOR_LISTS",&
            description="Activates the printing of the neighbor lists.", &
            print_level=high_print_level,each=(/1/),filename="",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_BANNER",&
            description="Controls the printing of the banner of the MM program",&
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of information regarding the run.",&
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
       
    CALL cp_print_key_section_create(print_key, "ff_info", description=&
         "controls the printing of information in the forcefield settings", &
         print_level=high_print_level,filename="__STD_OUT__",&
         error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_print_mm_section

!***************************************************************************
!!****f* input_cp2k_mm/create_forcefield_section [1.0] *
!!
!!   NAME
!!     create_forcefield_section
!!
!!   FUNCTION
!!     Create the forcefield section. This section is useful to set up the
!!     proper force_field for FIST calculations
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
  SUBROUTINE create_forcefield_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_forcefield_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="FORCEFIELD",&
            description="Section specifying information regarding how to set up properly"// &
            " a force_field for the classical calculations.",&
            n_keywords=2, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection,keyword)

       CALL keyword_create(keyword, name="PARMTYPE",&
            description="Specifies the type of parameters file to use.",&
            usage="PARMTYPE {OFF|CHM}",default_c_val="OFF",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARMFILE",&
            description="Specifies the filename that contains the parameters of the FF.",&
            usage="PARMFILE {FILENAME}",default_lc_val=" ",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT_NB",&
            description="Cutoff radius for nonbonded interactions",&
            usage="RCUT_NB {real}", type_of_var=real_t, default_r_val=-1.0_dp,&
            unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)                     

       CALL keyword_create(keyword, name="VDW_SCALE14",&
            description="Scaling factor for the VDW 1-4 ",&
            usage="VDW_SCALE14 1.0", default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EI_SCALE14",&
            description="Scaling factor for the electrostatics 1-4 ",&
            usage="EI_SCALE14 1.0", default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       !
       ! subsections
       !
       CALL create_NONBONDED_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_NONBONDED14_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_CHARGE_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_BOND_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_BEND_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_TORSION_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_forcefield_section

!***************************************************************************
!!****f* input_cp2k_mm/create_TORSION_section [1.0] *
!!
!!   NAME
!!     create_torsion_section
!!
!!   FUNCTION
!!     This section specifies the torsion of the MM atoms
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
  SUBROUTINE create_TORSION_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_TORSION_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="TORSION",&
            description="specifies the torsion potential of the MM system.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the tors.",&
            usage="ATOMS {KIND1} {KIND2} {KIND3} {KIND4}", required=.TRUE., type_of_var=char_t,&
            n_var=4, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KIND",&
            description="kind of tors",&
            usage="KIND {KIND1}", required=.FALSE., type_of_var=lchar_t,&
            default_lc_val="CHARMM",&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K",&
            description="Defines the force constant of the potential",&
            usage="K {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="hartree",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PHI0",&
            description="Defines the phase of the potential",&
            usage="PHI0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="rad", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="M",&
            description="Defines the multiplicity of the potential",&
            usage="M {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_TORSION_section

!***************************************************************************
!!****f* input_cp2k_mm/create_BEND_section [1.0] *
!!
!!   NAME
!!     create_bend_section
!!
!!   FUNCTION
!!     This section specifies the bend of the MM atoms
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
  SUBROUTINE create_BEND_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_BEND_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BEND",&
            description="specifies the bend potential of the MM system.",&
            n_keywords=1, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword,subsection)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the bend.",&
            usage="ATOMS {KIND1} {KIND2} {KIND3}", required=.TRUE., type_of_var=char_t,&
            n_var=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KIND",&
            description="kind of bend",&
            usage="KIND {KIND1}", required=.FALSE., type_of_var=lchar_t,&
            default_lc_val="CHARMM",&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K",&
            description="Defines the force constant of the potential",&
            usage="K {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="hartree*rad^-2",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="THETA0",&
            description="Defines the equilibrium angle.",&
            usage="THETA0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str='rad',error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(subsection,name="UB",&
            description="specifies the UB potential of the MM system.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="K",&
            description="force constant",&
            usage="K {real}", required=.TRUE., type_of_var=real_t,&
            default_r_val=0.0_dp,&
            n_var=-1, unit_str="hartree*bohr^-2",error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="R0",&
            description="Defines the equilibrium distance.",&
            usage="R0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="bohr", error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_BEND_section

!***************************************************************************
!!****f* input_cp2k_mm/create_BOND_section [1.0] *
!!
!!   NAME
!!     create_bond_section
!!
!!   FUNCTION
!!     This section specifies the bond of the MM atoms
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
  SUBROUTINE create_BOND_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_BOND_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BOND",&
            description="specifies the bond potential the MM system",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the bond.",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KIND",&
            description="kind of bond",&
            usage="KIND {KIND1}", required=.FALSE., type_of_var=lchar_t,&
            default_lc_val="CHARMM",&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="K",&
            description="Defines the force constant of the potential. For QUARTIC potentials"//&
            " three numbers are expected.",&
            usage="K {real}", required=.TRUE., type_of_var=real_t,&
            n_var=-1, unit_str="hartree*bohr^-2",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="R0",&
            description="Defines the equilibrium distance.",&
            usage="R0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="bohr", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_BOND_section

!***************************************************************************
!!****f* input_cp2k_mm/create_charge_section [1.0] *
!!
!!   NAME
!!     create_charge_section
!!
!!   FUNCTION
!!     This section specifies the charge of the MM atoms
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
  SUBROUTINE create_charge_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_charge_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="charge",&
            description="This section specifies the charge of the MM atoms",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOM",&
            description="Defines the atomic kind of the charge.",&
            usage="ATOM {KIND1}", required=.TRUE., type_of_var=char_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE",&
            description="Defines the charge of the MM atom in electron charge unit.",&
            usage="CHARGE {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_charge_section

!***************************************************************************
!!****f* input_cp2k_mm/create_NONBONDED14_section [1.0] *
!!
!!   NAME
!!     create_NONBONDED14_section
!!
!!   FUNCTION
!!     This section specifies the input parameters for 1-4 NON-BONDED 
!!     Interactions
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
  SUBROUTINE create_NONBONDED14_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_NONBONDED14_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="nonbonded14",&
            description="This section specifies the input parameters for 1-4 NON-BONDED interactions.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)
       
       CALL create_LJ_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Williams_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Goodwin_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_NONBONDED14_section

!***************************************************************************
!!****f* input_cp2k_mm/create_NONBONDED_section [1.0] *
!!
!!   NAME
!!     create_NONBONDED_section
!!
!!   FUNCTION
!!     This section specifies the input parameters for 1-4 NON-BONDED 
!!     Interactions
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
  SUBROUTINE create_NONBONDED_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_NONBONDED_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="nonbonded",&
            description="This section specifies the input parameters for 1-4 NON-BONDED interactions.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)
       
       CALL create_LJ_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_Williams_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_EAM_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)       

       CALL create_Goodwin_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_IPBV_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_BMHFT_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_NONBONDED_section

!***************************************************************************
!!****f* input_cp2k_mm/create_EAM_section [1.0] *
!!
!!   NAME
!!     create_EAM_section
!!
!!   FUNCTION
!!     This section specifies the input parameters for EAM  potential type
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
  SUBROUTINE create_EAM_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_EAM_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="EAM",&
            description="This section specifies the input parameters for EAM potential type.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARMFILE",&
            description="Specifies the filename that contains the parameters of the EAM FF.",&
            usage="PARMFILE {FILENAME}",default_lc_val=" ",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_EAM_section

!***************************************************************************
!!****f* input_cp2k_mm/create_LJ_section [1.0] *
!!
!!   NAME
!!     create_LJ_section
!!
!!   FUNCTION
!!     This section specifies the input parameters for Lennard-Jones potential type
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
  SUBROUTINE create_LJ_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_LJ_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="lennard-jones",&
            description="This section specifies the input parameters for LENNARD-JONES potential type.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPSILON",&
            description="Defines the EPSILON parameter of the LJ potential",&
            usage="EPSILON {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SIGMA",&
            description="Defines the SIGMA parameter of the LJ potential",&
            usage="SIGMA {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the LJ potential",&
            usage="RCUT {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_LJ_section


!***************************************************************************
!!****f* input_cp2k_mm/create_Williams_section [1.0] *
!!
!!   NAME
!!     create_Williams_section
!!
!!   FUNCTION
!!     This section specifies the input parameters for Williams potential type
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
  SUBROUTINE create_Williams_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_Williams_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="williams",&
            description="This section specifies the input parameters for WILLIAMS potential type.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="A",&
            description="Defines the A parameter of the Williams potential",&
            usage="A {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B",&
            description="Defines the B parameter of the Williams potential",&
            usage="B {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom^-1",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="C",&
            description="Defines the C parameter of the Williams potential",&
            usage="C {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e*angstrom^6",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the Williams potential",&
            usage="RCUT {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_Williams_section

!***************************************************************************
!!****f* input_cp2k_mm/create_Goodwin_section [1.0] *
!!
!!   NAME
!!     create_Goodwin_section
!!
!!   FUNCTION
!!     This section specifies the input parameters for Goodwin potential type
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
  SUBROUTINE create_Goodwin_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_Goodwin_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="goodwin",&
            description="This section specifies the input parameters for GOODWIN potential type.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VR0",&
            description="Defines the VR0 parameter of the Goodwin potential",&
            usage="VR0 {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="K_e",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="D",&
            description="Defines the D parameter of the Goodwin potential",&
            usage="D {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DC",&
            description="Defines the DC parameter of the Goodwin potential",&
            usage="DC {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="M",&
            description="Defines the M parameter of the Goodwin potential",&
            usage="M {real}", required=.TRUE., type_of_var=integer_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MC",&
            description="Defines the MC parameter of the Goodwin potential",&
            usage="MC {real}", required=.TRUE., type_of_var=integer_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the Goodwin potential",&
            usage="RCUT {real}", required=.TRUE., type_of_var=real_t,&
            n_var=1, unit_str="angstrom", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_Goodwin_section

!***************************************************************************
!!****f* input_cp2k_mm/create_ipbv_section [1.0] *
!!
!!   NAME
!!     create_ipbv_section
!!
!!   FUNCTION
!!     This section specifies the input parameters for IPBV potential type
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
  SUBROUTINE create_ipbv_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ipbv_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="ipbv",&
            description="This section specifies the input parameters for IPBV potential type.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the IPBV nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the IPBV potential",&
            usage="RCUT {real}", required=.TRUE., default_r_val=7.8_dp,&
            unit_str="angstrom", n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_ipbv_section

!***************************************************************************
!!****f* input_cp2k_mm/create_BMHFT_section [1.0] *
!!
!!   NAME
!!     create_BMHFT_section
!!
!!   FUNCTION
!!     This section specifies the input parameters for BMHFT potential type
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
  SUBROUTINE create_BMHFT_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_BMHFT_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BMHFT",&
            description="This section specifies the input parameters for BMHFT potential type.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kind involved in the BMHFT nonbond potential",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="RCUT",&
            description="Defines the cutoff parameter of the BMHFT potential",&
            usage="RCUT {real}", default_r_val=7.8_dp,&
            unit_str="angstrom", n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_BMHFT_section

!***************************************************************************

END MODULE input_cp2k_mm
