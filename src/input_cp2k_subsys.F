!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_cp2k_subsys *
!!
!!   NAME
!!     input_cp2k_subsys
!!
!!   FUNCTION
!!     builds the subsystem section of the input
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     teo & fawzi
!!
!!   MODIFICATION HISTORY
!!     10.2005 split input_cp2k [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE input_cp2k_subsys
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE input_constants
  USE input_cp2k_colvar,               ONLY: create_colvar_section
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
  USE physcon,                         ONLY: bohr
  USE string_utilities,                ONLY: s2a
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_subsys'

PUBLIC :: create_subsys_section
!***
!****************************************************************************
CONTAINS

!!****f* input_cp2k_subsys/create_cell_section *
!!
!!   NAME
!!     create_cell_section
!!
!!   FUNCTION
!!     creates the cell section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_cell_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_cell_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"cell",&
            description="parameters needed to set up the cell",&
            n_keywords=6, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, subsection)
! Real
       CALL keyword_create(keyword, name="A",&
            description="Read the hmat(1,1); hmat(2,1); hmat(3,1)",&
            usage="A real real real",&
            n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B",&
            description="Read the hmat(1,2); hmat(2,2); hmat(3,2)",&
            usage="B real real real",&
            n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="C",&
            description="Read the hmat(1,3); hmat(2,3); hmat(3,3)",&
            usage="C real real real",&
            n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ABC",&
            description="Read the hmat(1,1); hmat(2,2); hmat(3,3)",&
            usage="ABC real real real",&
            n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


! Character
       CALL keyword_create(keyword, name="UNITS",&
            variants=(/"UNIT"/),& 
            description="Specify the global units used in the input.",&
            usage="UNITS ANGSTROM",&
            enum_c_vals=s2a( "ANGSTROM","BOHR","SCALED_ANGSTROM","SCALED_BOHR"),&
            enum_i_vals=(/ use_unit_angstrom, use_unit_bohr,&
            use_unit_scaled_angstrom, use_unit_scaled_bohr/),&
            default_i_val=use_unit_bohr, error=error)
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


!Integers

       CALL keyword_create(keyword, name="SUBCELLS",&
            description="Read the subcells numbers.",&
            usage="SUBCELLS integer",&
            n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

!Subsections

       CALL create_cell_ref_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_cell_section

!***************************************************************************
!!****f* input_cp2k/create_rng_section [1.0] *
!!
!!   NAME
!!     create_rng_section
!!
!!   FUNCTION
!!     Creates the random number restart section
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
  SUBROUTINE create_rng_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_rng_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="RNG_INIT",&
            description="Information to initialize the parallel random number generator streams",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify an initial RNG stream",repeats=.TRUE.,&
            type_of_var=lchar_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_rng_section
!***************************************************************************

!!****f* input_cp2k_subsys/create_cell_ref_section *
!!
!!   NAME
!!     create_cell_ref_section
!!
!!   FUNCTION
!!     creates the cell section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_cell_ref_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_cell_ref_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"cell_ref",&
            description="parameters needed to set up the reference cell",&
            n_keywords=6, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
! Real
       CALL keyword_create(keyword, name="A",&
            description="Read the hmat(1,1); hmat(2,1); hmat(3,1)",&
            usage="A real real real",&
            n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B",&
            description="Read the hmat(1,2); hmat(2,2); hmat(3,2)",&
            usage="B real real real",&
            n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="C",&
            description="Read the hmat(1,3); hmat(2,3); hmat(3,3)",&
            usage="C real real real",&
            n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ABC",&
            description="Read the hmat(1,1); hmat(2,2); hmat(3,3)",&
            usage="ABC real real real",&
            n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


! Character
       CALL keyword_create(keyword, name="UNITS",&
            variants=(/"UNIT"/),& 
            description="Specify the global units used in the input.",&
            usage="UNITS ANGSTROM",&
            enum_c_vals=s2a( "ANGSTROM","BOHR","SCALED_ANGSTROM","SCALED_BOHR"),&
            enum_i_vals=(/ use_unit_angstrom, use_unit_bohr,&
            use_unit_scaled_angstrom, use_unit_scaled_bohr/),&
            default_i_val=use_unit_bohr, error=error)
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

!Integers

       CALL keyword_create(keyword, name="SUBCELLS",&
            description="Read the subcells numbers.",&
            usage="SUBCELLS integer",&
            n_var=1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_cell_ref_section

!!****f* input_cp2k_subsys/create_subsys_section [1.0] *
!!
!!   NAME
!!     create_subsys_section
!!
!!   FUNCTION
!!     creates the structure of a subsys, i.e. a full set of 
!!     atoms+mol+bounds+cell
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
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE create_subsys_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_subsys_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="subsys",&
            description="a subsystem: coordinates, topology, molecules and cell",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection,print_key,keyword)

       CALL create_rng_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_cell_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_coord_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_velocity_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)       

       CALL create_kind_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_constraint_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_topology_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="print",&
            description="Controls printings related to the subsys",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL cp_print_key_section_create(print_key,"atomic_coordinates",&
            description="controls the output of the atomic coordinates "//&
            "(check IF trajectory is not what you want)",&
            print_level=medium_print_level, filename="__STD_OUT__",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
       
       CALL cp_print_key_section_create(print_key,"interatomic_distances",&
            description="controls the output of the interatomic distances",&
            print_level=debug_print_level, filename="__STD_OUT__",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
       
       CALL cp_print_key_section_create(print_key, "topology_info", description=&
            "controls the printing of information in the topology settings", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL keyword_create(keyword,"pdb_info",&
            description="if the printkey is active prints information regarding PDB",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"xyz_info",&
            description="if the printkey is active prints information regarding XYZ",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"psf_info",&
            description="if the printkey is active prints information regarding PSF",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"util_info",&
            description="if the printkey is active prints information regarding topology utilities",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)  

       CALL cp_print_key_section_create(print_key,"cell",&
            description="controls the output of the cell parameters",&
            print_level=medium_print_level, filename="__STD_OUT__",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"subcell",&
            description="controls the output of the subcell parameters"//&
            " during the construction of the neighbouring lists.",&
            print_level=high_print_level, filename="__STD_OUT__",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
     
       CALL cp_print_key_section_create(print_key,"kinds",&
            description="controls the output of information on the kinds",&
            print_level=medium_print_level, filename="__STD_OUT__",error=error)
       CALL keyword_create(keyword, name="potential",&
            description="If the printkey is activated controls the printing of the"//&
            " elp_potential, gth_potential or all electron potential information",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL keyword_create(keyword, name="basis_set",&
            description="If the printkey is activated controls the printing of basis set information",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
       
       CALL cp_print_key_section_create(print_key,"molecules",&
            description="controls the output of information on the molecules",&
            print_level=medium_print_level, filename="__STD_OUT__",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"radii",&
            description="controls the output of radii information",&
            print_level=high_print_level, filename="__STD_OUT__",error=error)
       CALL keyword_create(keyword, name="core_charges_radii",&
            description="If the printkey is activated controls the printing of the radii of the core charges",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL keyword_create(keyword, name="pgf_radii",&
            description="If the printkey is activated controls the printing of the core gaussian radii",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error) 
       
       CALL keyword_create(keyword, name="set_radii",&
            description="If the printkey is activated controls the printing of the set_radii",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       
       CALL keyword_create(keyword, name="kind_radii",&
            description="If the printkey is activated controls the printing of the kind_radii",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       
       CALL keyword_create(keyword, name="core_charge_radii",&
            description="If the printkey is activated controls the printing of the core_charge_radii",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="gth_ppl_radii",&
            description="If the printkey is activated controls the printing of the "//&
            "gth pseudo potential local radii",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL keyword_create(keyword, name="gth_ppnl_radii",&
            description="If the printkey is activated controls the printing of the "//&
            "gth pseudo potential non local radii",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL keyword_create(keyword, name="gapw_prj_radii",&
            description="If the printkey is activated controls the printing of the gapw projector radii",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)       
       
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)       
       
    END IF
  END SUBROUTINE create_subsys_section

!***************************************************************************
!!****f* input_cp2k_subsys/create_velocity_section [1.0] *
!!
!!   NAME
!!     create_velocity_section
!!
!!   FUNCTION
!!     Creates the velocity section
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
  SUBROUTINE create_velocity_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_velocity_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="velocity",&
            description="The velocities for simple systems xyz format by default",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="The atoms velocities in xyz format. The same order as"//&
            " the atom's order is assumed.",repeats=.TRUE.,required=.TRUE.,&
            usage="<REAL> <REAL> <REAL>", type_of_var=real_t, n_var=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_velocity_section

!***************************************************************************
!!****f* input_cp2k_subsys/create_coord_section [1.0] *
!!
!!   NAME
!!     create_coord_section
!!
!!   FUNCTION
!!     Creates the &COORD section
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
  SUBROUTINE create_coord_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_coord_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="coord",&
            description="The coordinates for simple systems (like the QM ones)"//&
            " xyz format by default. More complex systems should be given with"//&
            " an external pdb file.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="The atoms coordinates in xyz format.",repeats=.TRUE.,&
            usage="H 0.7 0.3 0.3", type_of_var=lchar_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_coord_section


!***************************************************************************
!!****f* input_cp2k_subsys/create_kind_section [1.0] *
!!
!!   NAME
!!     create_kind_section
!!
!!   FUNCTION
!!     Creates the QM/MM section
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
  SUBROUTINE create_kind_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_kind_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="kind",&
            description="The description of the kind of the atoms (mostly for QM)",&
            n_keywords=5, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="The name of the kind described in this section.",&
            usage="H", default_c_val="DEFAULT", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AUX_BASIS_SET",&
            variants=s2a("AUXILIARY_BASIS_SET", "AUX_BASIS"),&
            description="The auxliliary basis set (GTO type)",&
            usage="AUX_BASIS_SET DZVP", default_c_val=" ", &
            n_var=1, required=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AUX_BASIS_NORMALIZATION",&
            variants=s2a( "AUXILIARY_BASIS_NORMALIZATION", "AUX_BASIS_NORM"),&
            description="The normalization of the auxliliary basis set",&
            usage="AUX_BASIS_NORMALIZATION NO", &
            enum_c_vals=s2a( "NONE","NO","WFN","WAVEFUNCTION","FUNCTION",&
            "DENS","DENSITY","SQUARE","UNDEFINED"),&
            enum_i_vals=(/0,0,1,1,1,2,2,2,-1/),&
            default_i_val=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BASIS_SET",&
            variants=s2a("ORBITAL_BASIS_SET","ORB_BASIS"),&
            description="The basis set (GTO type)",&
            usage="BASIS_SET DZVP", default_c_val=" ", &
            n_var=1, required=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BASIS_NORMALIZATION",&
            variants=s2a("ORBITAL_BASIS_NORMALIZATION","ORB_BASIS_NORM"),&
            description="The normalization of the auxliliary basis set",&
            usage="AUX_BASIS_NORMALIZATION NO", &
            enum_c_vals=s2a( "NONE", "NO", "WFN", "WAVEFUNCTION",&
            "FUNCTION", "DENS", "DENSITY", "SQUARE", "UNDEFINED"),&
            enum_i_vals=(/0,0,1,1,1,2,2,2,-1/),&
            default_i_val=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORE_CORRECTION",&
            description="Corrects the effective nuclear charge",&
            usage="CORE_CORRECTION 1.0", n_var=1, required=.FALSE., &
            default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ELEMENT",&
            variants=(/ "ELEMENT_SYMBOL" /),&
            description="The element of the actual kind "//&
            "(if not given it is inferred from the kind name)",&
            usage="ELEMENT O", type_of_var=char_t,n_var=1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MASS",&
            variants=s2a("ATOMIC_MASS", "ATOMIC_WEIGHT", "WEIGHT"),&
            description="The mass of the atom "//&
            "(if negative or non present it is inferred from the element symbol)",&
            usage="MASS 2.0", type_of_var=real_t,n_var=1,&
            required=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="POTENTIAL",&
            variants=(/ "POT" /),&
            description="The element of the actual kind "//&
            "(if not given it is inferred from the kind name)",&
            usage="POTENTIAL ghost", type_of_var=char_t, n_var=1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="HARD_EXP_RADIUS",&
            description="The region where the hard density is supposed to be confined"//&
                        "(GAPW)(in Bohr, default is 1.2 for H and 1.512 otherwise)",&
            usage="HARD_EXP_RADIUS 0.9", type_of_var=real_t,n_var=1,&
            required=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_RAD_LOCAL",&
            description="Max radius for the basis functions used to"//&
            " generate the local projectors in GAPW [Bohr]",&
            usage="MAX_RAD_LOCAL 15.0", default_r_val=13.0_dp*bohr,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RHO0_EXP_RADIUS",&
            description="the radius which defines the atomic region where "//&
            "the hard compensation density is confined."//&
            "should be less than HARD_EXP_RADIUS (GAPW)(Bohr, default equals HARD_EXP_RADIUS)",&
            usage="RHO_EXP_RADIUS 0.9", type_of_var=real_t,n_var=1,&
            required=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LEBEDEV_GRID",&
            description="The number of points for the angular part of "//&
            "the local grid (GAPW)",&
            usage="LEBEDEV_GRID 40", default_i_val=50,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RADIAL_GRID",&
            description="The number of points for the radial part of "//&
            "the local grid (GAPW)",&
            usage="RADIAL_GRID 70", default_i_val=50,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Logicals
       CALL keyword_create(keyword, name="GPW_TYPE",&
          description="Force one type to be treated by the GPW scheme,"//&
          " whatever are its primitives, even if the GAPW method is used",&
          usage="GPW_TYPE",&
          default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_kind_section

!***************************************************************************

!!****f* input_cp2k_subsys/create_constraint_section [1.0] *
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_constraint_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="constraint",&
            description="Section specifying information regarding how to impose constraints"// &
            " on the system.",&
            n_keywords=0, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection, keyword)
       !
       ! subsections
       !
       CALL keyword_create(keyword, name="SHAKE",&
            description="Set the tolerance for the shake constraint algorithm.",&
            usage="SHAKE real",&
            default_r_val=1.0E-6_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_hbonds_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_g3x3_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_g4x6_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_fixed_atom_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_collective_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_constraint_section
!***************************************************************************

!!****f* input_cp2k_subsys/create_collective_section [1.0] *
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_collective_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="COLLECTIVE",&
            description="This section is used to constraint collective (general) degrees of freedom"//&
            "using general routines, and writing langrangian multipliers to file ",&
            n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="COLVAR",&
             description="Specifies which colvar to use ",&
             usage="COLVAR {int}",&
             type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOLECULE",&
            variants=(/"MOL"/),&
            description="Specifies the molecule number on which constraint will be applied",&
            usage="MOL {integer}", required=.TRUE., n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TARGET",&
            description="Specifies the target value of the constrained collective variable (units depend on the colvar) ",&
            usage="TARGET {real}",&
            type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_collective_section

!***************************************************************************
!!****f* input_cp2k_subsys/create_fixed_atom_section [1.0] *
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
            usage="LIST {integer} {integer} .. {integer}", required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RANGE",&
            description="Specifies the starting and the ending index of atoms range you wish to freeze.",&
            usage="RANGE {integer} {integer}", required=.TRUE.,&
            n_var=2, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOLNAME",&
            description="Specifies the name of the molecule to fix",&
            usage="MOLNAME WAT MEOH", required=.TRUE., repeats=.TRUE.,&
            n_var=-1,type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MM_SUBSYS",&
            variants=(/"PROTEIN"/),&
            description="In a QM/MM run all  MM atoms are fixed. Only the QM fragment is allowed to move."//&
            " Default .FALSE.",&
            usage="MM_SUBSYS logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="QM_SUBSYS",&
            description="In a QM/MM run all QM atoms are fixed. Only the MM fragment is allowed to move."//&
            " Default .FALSE.",&
            usage="QM_SUBSYS logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_fixed_atom_section

!***************************************************************************
!!****f* input_cp2k_subsys/create_g3x3_section [1.0] *
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
            description="Specifies the molecule number on which constraint will be applied",&
            usage="MOL {integer}", required=.TRUE., n_var=1, type_of_var=integer_t, error=error)
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
            usage="DISTANCES {real} {real} {real}", type_of_var=real_t, n_var=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)              

    END IF

  END SUBROUTINE create_g3x3_section

!***************************************************************************
!!****f* input_cp2k_subsys/create_hbonds_section [1.0] *
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_hbonds_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="HBONDS",&
            description="This section is used to set bonds constraints involving Hydrogen atoms",&
            n_keywords=3, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
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

       ! Real
       CALL keyword_create(keyword, name="TARGETS",&
            description="The constrained distances' values  for the types defines in ATOM_TYPE.",&
            usage="TARGETS {real} {real} {real}", type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)              

    END IF

  END SUBROUTINE create_hbonds_section


!***************************************************************************
!!****f* input_cp2k_subsys/create_g4x6_section [1.0] *
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
            description="Specifies the molecule number on which constraint will be applied",&
            usage="MOL {integer}", required=.TRUE., n_var=1, type_of_var=integer_t, error=error)
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
            type_of_var=real_t, n_var=6,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)              

    END IF

  END SUBROUTINE create_g4x6_section

!***************************************************************************
!!****f* input_cp2k_subsys/create_topology_section [1.0] *
!!
!!   NAME
!!     create_topology_section
!!
!!   FUNCTION
!!     Create the topology section for FIST.. and the base is running running...
!!     Contains all information regarding topology to be read in input file..
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
  SUBROUTINE create_topology_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_topology_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="topology",&
            description="Section specifying information regarding how to handle the topology"// &
            " for classical runs.",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       ! Logical
       CALL keyword_create(keyword, name="CHARGE_OCCUP",&
            variants=(/ "CHARGE_O" /),&
            description="Read MM charges from the OCCUP field of PDB file.",&
            usage="CHARGE_OCCUP logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE_BETA",&
            variants=(/ "CHARGE_B" /),&
            description="Read MM charges from the BETA field of PDB file.",&
            usage="CHARGE_BETA logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARA_RES",&
            description="For a protein, each residue is now considered a molecule",&
            usage="PARA_RES logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DUMP_TOP",&
            description="This keyword enables the dumping of the topology as a PSF file.",&
            usage="DUMP_TOP logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="REORDER",&
            description="Reorder a list of atomic coordinates into order so it can be packed correctly.",&
            usage="REORDER logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOL_CHECK",&
            description="Check molecules have the same number of atom and names.",&
            usage="MOL_CHECK logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Character
       CALL keyword_create(keyword, name="COORD_FILE",&
            variants=s2a("COORD_FI", "COORD_FILE_NAME"),&
            description="Specifies the filename that contains coordinates.",&
            usage="COORD_FILE FILENAME",required=.TRUE.,type_of_var=lchar_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CONN_FILE",&
            variants=(/"CONN_FIL"/),&
            description="Specifies the filename that contains the molecular connectivity.",&
            usage="CONN_FILE FILENAME",required=.TRUE.,type_of_var=lchar_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="COORDINATE",&
            variants=s2a( "COORDINA", "COORD"),&
            description="Set up the way in which coordinates will be read.",&
            usage="COORDINATE (PDB||XYZ||OFF)", &
            enum_c_vals=s2a( "PDB","XYZ","OFF"),&
            enum_i_vals=(/do_coord_pdb, do_coord_xyz, do_coord_off/),&
            default_i_val=do_coord_off,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
 
       CALL keyword_create(keyword, name="CONNECTIVITY",&
            variants=(/ "CONNECTI"/),&
            description="Set up the way in which connectivities will be read. "// &
            "Default generates topology",&
            usage="CONNECTIVITY (PSF||GENERATE||OFF)", &
            enum_c_vals=s2a( "PSF","GENERATE","OFF"),&
            enum_i_vals=(/do_conn_psf, do_conn_generate,&
                          do_conn_off/),&
            default_i_val=do_conn_generate,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXCLUDE_VDW",&
            description="Specifies which kind of Van der Waals interaction to skip.",&
            usage="EXCLUDE_VDW (1-2||1-3||1-4)", &
            enum_c_vals=s2a( "1-2", "1-3", "1-4"),&
            enum_i_vals=(/do_skip_vdw_12, do_skip_vdw_13, do_skip_vdw_14/),&
            default_i_val=do_skip_vdw_14,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXCLUDE_EI",&
            description="Specifies which kind of Electrostatic interaction to skip.",&
            usage="EXCLUDE_EI (1-2||1-3||1-4)", &
            enum_c_vals=s2a( "1-2", "1-3", "1-4"),&
            enum_i_vals=(/do_skip_ei_12, do_skip_ei_13, do_skip_ei_14/),&
            default_i_val=do_skip_ei_14,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
 
       CALL keyword_create(keyword, name="BONDPARM",&
            description="Used in conjunction with BONDPARM_FACTOR to "//&
                        "help determine wheather there is bonding "//&
                        "between two atoms based on a distance criteria. "//&
                        "Can use covalent radii information or VDW radii information",&
            usage="BONDPARM (COVALENT||VDW)", &
            enum_c_vals=s2a( "COVALENT", "VDW"),&
            enum_i_vals=(/do_bondparm_covalent, do_bondparm_vdw/),&
            default_i_val=do_bondparm_covalent,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NHCOPT",&
            description="Determines the number of DOF each NHC is attached to.",&
            usage="NHCOPT (GLOBAL||MOLECULE||MASSIVE||NONE)", &
            enum_c_vals=s2a( "GLOBAL", "MOLECULE", "MASSIVE", "NONE"),&
            enum_i_vals=(/do_nhcopt_global, do_nhcopt_molecule,&
                          do_nhcopt_massive, do_nhcopt_none/),&
            default_i_val=do_nhcopt_global,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       ! Integer
       CALL keyword_create(keyword, name="EXCLUDE_CHECK",&
            variants=(/"EXCLUDE_CH"/),&
            description="Used in conjunction with PARA_RES. Determines "//&
                        "the number of molecules ahead and behind that "//&
                        "needs to be checked in order to generate the "//&
                        "proper exclusion list. 3 is safe for all "//&
                        "molecules situations but 1 should work for "//&
                        "all proteins.",&
            usage="exclude_check {integer}", default_i_val=3,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Real
       CALL keyword_create(keyword, name="BONDPARM_FACTOR",&
            description="Used in conjunction with BONDPARM to help "//&
                        "determine wheather there is bonding between "//&
                        "two atoms based on a distance criteria.",&
            usage="bondparm_factor {real}", default_r_val=1.1_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_topology_section

END MODULE input_cp2k_subsys
