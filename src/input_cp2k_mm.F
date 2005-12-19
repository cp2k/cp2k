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
  use kinds, only: dp
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_vals_type,&
       section_vals_get_subs_vals, section_vals_get
  USE input_constants
  USE input_keyword_types, ONLY: keyword_type, keyword_create,&
       keyword_release
  USE input_val_types, ONLY: val_type, no_t, logical_t, integer_t, real_t, char_t,&
       enum_t
  USE cp_output_handling, only: cp_print_key_section_create

  USE string_utilities, ONLY: s2a
  USE input_cp2k_dft, only: create_mt_section
#include "cp_common_uses.h"

  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='m'

  public :: create_mm_section
  PUBLIC :: create_NONBONDED14_section, create_LJ_section, create_Williams_section, create_Goodwin_section
!***
!****************************************************************************
contains

!!****f* input_cp2k_mm/create_poisson_mm_section [1.0] *
!!
!!   NAME
!!     create_poisson_mm_section
!!
!!   FUNCTION
!!     Creates the Poisson section
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
  SUBROUTINE create_poisson_mm_section(section,error)
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_poisson_mm_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword
    TYPE(section_type), pointer :: subsection

    failure=.false.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="poisson_mm",&
            description="Sets up the poisson resolutor.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, subsection)
       CALL keyword_create(keyword, name="POISSON_SOLVER",&
            variants=(/"POISSON", "PSOLVER"/),& 
            description="Specify which kind of solver to use to solve the Poisson equation.",&
            usage="POISSON_SOLVER char",&
            enum_c_vals=s2a( "PERIODIC", "ANALYTIC", "MT"),&
            enum_i_vals=(/ use_periodic, use_analytic, use_mt/),&
            default_i_val=use_periodic, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

       CALL keyword_create(keyword, name="PERIODIC",&
            description="Specify the directions on wich apply PBC ",&
            usage="PERIODIC (x|y|z|xy|xz|yz|xyz|none)",&
            enum_c_vals=s2a( "x","y","z","xy","xz","yz","xyz","none"),&
            enum_i_vals=(/ use_perd_x,  use_perd_y,   use_perd_z,&
                           use_perd_xy, use_perd_xz, use_perd_yz,&
                           use_perd_xyz, use_perd_none /),&
            default_i_val=use_perd_xyz, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)     

       CALL create_mt_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_ewald_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_poisson_mm_section
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
    TYPE(section_type), pointer :: section
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='create_mm_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), pointer :: subsection

    failure=.false.

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
       
       CALL create_poisson_mm_section(subsection,error)
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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_print_mm_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: print_key
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options in MM code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(print_key)

       CALL cp_print_key_section_create(print_key,"ATOMIC_COORDINATES",&
            description="Controls the printing of atomic coordinates during MM.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"BASIC_DATA_TYPES",&
            description="Controls the printing of the basic data types.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"BLACS_INFORMATION",&
            description="Controls the printing of information about the specified BLACS context.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CELL_PARAMETERS",&
            description="Controls the printing of cell information", &
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CPUTIME",&
            description="Controls the printing of cputime information along the MM execution.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DERIVATIVES",&
            description="Controls the printing of derivatives.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DIPOLE",&
            description="Controls the printing of dipole information. Computed with"//&
            " the Berry phase formalism. OutPut in atomic Units.", &
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DISTRIBUTION",&
            description="Controls the printing of the distribution of the simple pair neighbour list.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"INTERATOMIC_DISTANCES",&
            description="Controls the printing of interatomic distances.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"KIND_INFORMATION",&
            description="Controls the printing of information on atomic kinds.", &
            print_level=low_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MEMORY",&
            description="Controls the printing of information on memory usage", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"NEIGHBOUR_LISTS",&
            description="Controls the printing of the set of neighbor lists.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PHYSICAL_CONSTANTS",&
            description="Controls the printing of physical constants used in the code.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"POTENTIAL",&
            description="Controls the printing of the potential generated by the MM "//&
            "charge distribution on a cube file.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PRECISION",&
            description="Controls the precision of the results printed.",&
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
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

       CALL cp_print_key_section_create(print_key,"PW_GRID_INFORMATION",&
            description="Controls the printing of information regarding the pw_grid structure.",&
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"TOTAL_NUMBERS",&
            description="Controls the printing of the total number of kinds, atoms, etc.",&
            print_level=low_print_level,each=(/1/),filename="__STD_OUT__",&
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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_forcefield_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    TYPE(section_type), POINTER :: subsection
    
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
            usage="PARMFILE {FILENAME}",default_c_val=" ",error=error)
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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_TORSION_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="TORSION",&
            description="This section specifies the torsion potential of the MM system."//&
            "JORGENSEN type potential supported.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the bond.",&
            usage="ATOMS {KIND1} {KIND2} {KIND3} {KIND4}", required=.TRUE., type_of_var=char_t,&
            n_var=4, error=error)
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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_BEND_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BEND",&
            description="This section specifies the bend potential of the MM system."//&
            "HARMONIC type potential supported.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the bond.",&
            usage="ATOMS {KIND1} {KIND2} {KIND3}", required=.TRUE., type_of_var=char_t,&
            n_var=3, error=error)
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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_BOND_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BOND",&
            description="This section specifies the potential of the bond of the MM system",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the bond.",&
            usage="ATOMS {KIND1} {KIND2}", required=.TRUE., type_of_var=char_t,&
            n_var=2, error=error)
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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_charge_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="charge",&
            description="This section specifies the charge of the MM atoms",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_NONBONDED14_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    
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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_NONBONDED_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    
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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_EAM_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
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
            usage="PARMFILE {FILENAME}",default_c_val=" ",error=error)
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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_LJ_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_Williams_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="williams",&
            description="This section specifies the input parameters for WILLIAMS potential type.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_Goodwin_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="goodwin",&
            description="This section specifies the input parameters for GOODWIN potential type.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_ipbv_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="ipbv",&
            description="This section specifies the input parameters for IPBV potential type.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_BMHFT_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), POINTER :: keyword
    
    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BMHFT",&
            description="This section specifies the input parameters for BMHFT potential type.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
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

!!****f* input_cp2k_mm/create_ewald_section *
!!
!!   NAME
!!     create_ewald_section
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: will contain the ewald section 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************

SUBROUTINE create_ewald_section(section,error)
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='create_ewald_section',&
        routineP=moduleN//':'//routineN
  TYPE(keyword_type), pointer :: keyword

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  if (.not. failure) then
    CALL section_create(section,name="ewald",&
         description="ewald parameters",&
         n_keywords=7, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)

    nullify(keyword)
    CALL keyword_create(keyword, name="EWALD_TYPE",&
         description="The type of ewald you want to perform."//&
         " EWALD is the standard non-fft based ewald."//&
         " SPME is the smooth particle mesh using beta-Euler splines."//&
         " PME is the particle mesh using fft interpolation.",&
         usage="EWALD_TYPE {EWALD|SPME|PME}",&
         default_i_val=do_ewald_ewald,&
         enum_c_vals=ewald_name_types,&
         enum_i_vals=ewald_name_index,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="alpha",&
         description="alpha parameter associated with Ewald (EWALD|PME|SPME)",&
         usage="alpha .44", default_r_val=.44_dp, unit_str='angstrom^-1',error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="gmax",&
         description="maximum number of grid points (SPME and EWALD)",&
         usage="gmax [integer]", default_i_val=25,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ns_max",&
         description="number of grid points on small mesh (PME only)",&
         usage="ns_max 12", default_i_val=12,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="o_spline",&
         description="order of the beta-Euler spline (SPME only)",&
         usage="o_spline 6", default_i_val=6,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="epsilon",&
         description="tolerance of gaussians for fft interpolation (PME only)",&
         usage="epsilon 1e-6", default_r_val=1.e-6_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

  end if
END SUBROUTINE create_ewald_section

!***************************************************************************

END MODULE input_cp2k_mm
