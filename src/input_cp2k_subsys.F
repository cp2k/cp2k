!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
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
  USE cp_units,                        ONLY: cp_unit_to_cp2k
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

PUBLIC :: create_subsys_section,&
          create_cell_section,&
          create_structure_data
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
  SUBROUTINE create_cell_section(section,periodic,error)
    TYPE(section_type), POINTER              :: section
    INTEGER, INTENT(IN), OPTIONAL            :: periodic
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_cell_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: my_periodic
    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    my_periodic = use_perd_xyz
    IF (PRESENT(periodic)) my_periodic = periodic
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"CELL",&
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

       CALL keyword_create(keyword, name="SUBCELLS",&
            description="Read the subcells numbers.",&
            usage="SUBCELLS <real>",&
            n_var=1,default_r_val=2.0_dp,error=error)
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
            description="Specify the directions on wich apply PBC. Important notice, "//&
                        " this only applies to the generation of the pair lists."//&
                        " See the POISSON section to specify the periodicity used for the electrostatics."//&
                        " Typically the settings should be the same.",&
            usage="PERIODIC (x|y|z|xy|xz|yz|xyz|none)",&
            enum_c_vals=s2a( "x","y","z","xy","xz","yz","xyz","none"),&
            enum_i_vals=(/ use_perd_x,  use_perd_y,   use_perd_z,&
                           use_perd_xy, use_perd_xz, use_perd_yz,&
                           use_perd_xyz, use_perd_none /),&
            default_i_val=my_periodic, error=error)
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
    TYPE(cp_error_type), INTENT(inout)       :: error

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
    TYPE(cp_error_type), INTENT(inout)       :: error

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

       CALL keyword_create(keyword, name="SUBCELLS",&
            description="Read the subcells numbers.",&
            usage="SUBCELLS <real>",&
            n_var=1,default_r_val=2.0_dp,error=error)
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
    TYPE(cp_error_type), INTENT(inout)       :: error

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
            n_keywords=0, n_subsections=9, repeats=.FALSE., required=.TRUE.,&
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

       CALL create_topology_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_shell_coord_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_shell_vel_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_core_coord_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_core_vel_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="print",&
            description="Controls printings related to the subsys",&
            n_keywords=0, n_subsections=9, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL cp_print_key_section_create(print_key,"atomic_coordinates",&
            description="controls the output of the atomic coordinates when setting up the"//&
            "force environment. For printing coordinates during MD or GEO refer to the keyword"//&
            " trajectory.",&
            print_level=medium_print_level, filename="__STD_OUT__",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL create_structure_data(print_key, error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)           

       CALL cp_print_key_section_create(print_key,"interatomic_distances",&
            description="controls the output of the interatomic distances when setting up the"//&
            "force environment",&
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
       CALL keyword_create(keyword,"g96_info",&
            description="if the printkey is active prints information regarding G96",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"gtop_info",&
            description="if the printkey is active prints information regarding GROMOS Topology",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"util_info",&
            description="if the printkey is active prints information regarding topology utilities",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
       CALL keyword_create(keyword,"generate_info",&
            description="if the printkey is active prints information regarding topology generation",&
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

!!****f* input_cp2k_subsys/create_structure_data [1.0] *
!!
!!   NAME
!!     create_structure_data
!!
!!   FUNCTION
!!     creates structure data section for output.. both subsys  (for initialization)
!!     and motion section..
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
!!     
!!
!!*** **********************************************************************
  SUBROUTINE create_structure_data(print_key,error)
    TYPE(section_type), POINTER              :: print_key
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_structure_data', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(print_key),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(keyword)
       CALL cp_print_key_section_create(print_key,&
                                        name="STRUCTURE_DATA",&
                                        description="Request the printing of special "//&
                                                    "structure data during a structure "//&
                                                    "optimization (in MOTION%PRINT) or when "//&
                                                    "setting up a subsys (in SUBSYS%PRINT).",&
                                        print_level=high_print_level,&
                                        filename="__STD_OUT__",&
                                        error=error)
       CALL keyword_create(keyword,&
                           name="POSITION",&
                           variants=(/"POS"/),&
                           description="Print the position vector of "//&
                                       "atom a specified by its index",&
                           usage="POSITION {integer}",&
                           n_var=1,&
                           repeats=.TRUE.,&
                           required=.FALSE.,&
                           type_of_var=integer_t,&
                           error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,&
                           name="DISTANCE",&
                           variants=(/"DIS"/),&
                           description="Print the distance between the "//&
                                       "atoms a and b specified by their "//&
                                       "indices",&
                           usage="DISTANCE {integer} {integer}",&
                           n_var=2,&
                           repeats=.TRUE.,&
                           required=.FALSE.,&
                           type_of_var=integer_t,&
                           error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,&
                           name="ANGLE",&
                           variants=(/"ANG"/),&
                           description="Print the angle formed by the "//&
                                       "atoms specified by their indices",&
                           usage="ANGLE {integer}  {integer} {integer}",&
                           n_var=3,&
                           repeats=.TRUE.,&
                           required=.FALSE.,&
                           type_of_var=integer_t,&
                           error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,&
                           name="DIHEDRAL_ANGLE",&
                           variants=s2a("DIHEDRAL","DIH"),&
                           description="rint the dihedral angle between the "//&
                                       "planes defined by the atoms (a,b,c) and "//&
                                       "the atoms (b,c,d) specified by their indices",&
                           usage="DIHEDRAL_ANGLE {integer}  {integer} {integer} {integer}",&
                           n_var=4,&
                           repeats=.TRUE.,&
                           required=.FALSE.,&
                           type_of_var=integer_t,&
                           error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_structure_data

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
    TYPE(cp_error_type), INTENT(inout)       :: error

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
!!****f* input_cp2k_subsys/create_shell_vel_section [1.0] *
!!
!!   NAME
!!     create_shell_vel_section
!!
!!   FUNCTION
!!     Creates the shell velocity section
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
  SUBROUTINE create_shell_vel_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_shell_vel_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="shell_velocity",&
            description="The velocities of shells for shell-model potentials, "//&
            "in xyz format  ",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="The shell velocities in xyz format. The same order as"//&
            " the shell coordinate is assumed.",repeats=.TRUE.,required=.TRUE.,&
            usage="<REAL> <REAL> <REAL> ", type_of_var=real_t, n_var=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_shell_vel_section

!***************************************************************************
!!****f* input_cp2k_subsys/create_core_vel_section [1.0] *
!!
!!   NAME
!!     create_core_vel_section
!!
!!   FUNCTION
!!     Creates the shell velocity section
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
  SUBROUTINE create_core_vel_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_core_vel_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="core_velocity",&
            description="The velocities of cores for shell-model potentials, "//&
            "in xyz format  ",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="The core velocities in xyz format. The same order as"//&
            " the shell coordinate is assumed.",repeats=.TRUE.,required=.TRUE.,&
            usage="<REAL> <REAL> <REAL> ", type_of_var=real_t, n_var=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_core_vel_section

!***************************************************************************
!!****f* input_cp2k_subsys/create_potential_section [1.0] *
!!
!!   NAME
!!     create_potential_section
!!
!!   FUNCTION
!!     Creates the &POTENTIAL section
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
  SUBROUTINE create_potential_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_potential_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    IF (.NOT. failure) THEN
       CALL section_create(section,name="potential",&
            description="Section used to specify Potentials.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="CP2K Pseudo Potential Standard Format (GTH, ALL or ELP)",&
            repeats=.TRUE.,type_of_var=lchar_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_potential_section

!***************************************************************************
!!****f* input_cp2k_subsys/create_basis_section [1.0] *
!!
!!   NAME
!!     create_basis_section
!!
!!   FUNCTION
!!     Creates the &BASIS section
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
  SUBROUTINE create_basis_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_basis_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    IF (.NOT. failure) THEN
       CALL section_create(section,name="basis",&
            description="Section used to specify a general basis set for QM calculations.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="CP2K Basis Set Standard Format",repeats=.TRUE.,&
            type_of_var=lchar_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_basis_section

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
    TYPE(cp_error_type), INTENT(inout)       :: error

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
!!****f* input_cp2k_subsys/create_shell_coord_section [1.0] *
!!
!!   NAME
!!     create_shell_coord_section
!!
!!   FUNCTION
!!     Creates the &SHELL_COORD section
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
  SUBROUTINE create_shell_coord_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_shell_coord_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="shell_coord",&
            description="The shell coordinates for the shell-model potentials"//&
            " xyz format with an additional column for the index of the corresponding particle",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="The shell coordinates in xyz format.",repeats=.TRUE.,&
            usage="at_name rx ry rz at_index", type_of_var=lchar_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_shell_coord_section

!***************************************************************************
!!****f* input_cp2k_subsys/create_core_coord_section [1.0] *
!!
!!   NAME
!!     create_core_coord_section
!!
!!   FUNCTION
!!     Creates the &core_COORD section
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
  SUBROUTINE create_core_coord_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_core_coord_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="core_coord",&
            description="The core coordinates for the shell-model potentials"//&
            " xyz format with an additional column for the index of the corresponding particle",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="The core coordinates in xyz format.",repeats=.TRUE.,&
            usage="at_name rx ry rz at_index", type_of_var=lchar_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_core_coord_section


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
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_kind_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="kind",&
            description="The description of the kind of the atoms (mostly for QM)",&
            n_keywords=17, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
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

       CALL keyword_create(keyword, name="ELEC_CONF",&
            description="Specifies the electronic configration used in construction the "// &
                        "atomic initial guess (see the pseudo potential file for the default values.",&
            usage="ELEC_COND n_elec(s)  n_elec(p)  n_elec(d)  ... ", required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
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
            description="The name of the pseudopotential for the defined kind.",&
            usage="POTENTIAL <PSEUDO-POTENTIAL-NAME>", default_c_val=" ", n_var=1,&
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

       CALL keyword_create(keyword,&
                           name="GHOST",&
                           description="This keyword makes all atoms of this kind "//&
                                        "ghost atoms, i.e. without pseudo or nuclear charge."//&
                                       "Useful to just have the basis set at that position (BSSE calculations).",&
                           usage="GHOST",&
                           default_l_val=.FALSE.,&
                           lone_keyword_l_val=.TRUE.,&
                           error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsection)
       CALL create_basis_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_potential_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
    END IF

  END SUBROUTINE create_kind_section

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
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_topology_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection, &
                                                subsubsection, &
                                                subsubsubsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="TOPOLOGY",&
            description="Section specifying information regarding how to handle the topology"// &
            " for classical runs.",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, print_key)
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
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOL_CHECK",&
            description="Check molecules have the same number of atom and names.",&
            usage="MOL_CHECK logical",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,"DUMP_PDB",&
            description="controls the dumping of the PDB at the starting geometry",&
            print_level=debug_print_level, filename="topology_pdb.pdb",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DUMP_PSF",&
            description="controls the dumping of the PSF connectivity",&
            print_level=debug_print_level, filename="topology_psf.psf",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL keyword_create(keyword, name="USE_G96_VELOCITY",&
            description="Use the velocities in the G96 coordinate files as the starting velocity",&
            usage="USE_G96_VELOCITY logical",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Character
       CALL keyword_create(keyword, name="COORD_FILE_NAME",&
            variants=s2a("COORD_FILE"),&
            description="Specifies the filename that contains coordinates.",&
            usage="COORD_FILE FILENAME",required=.TRUE.,type_of_var=lchar_t,&
            default_lc_val="UNSET",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COORDINATE",&
            description="Set up the way in which coordinates will be read.",&
            usage="COORDINATE (PDB||XYZ||G96||OFF)", &
            enum_c_vals=s2a( "PDB","XYZ","G96","OFF"),&
            enum_i_vals=(/do_coord_pdb, do_coord_xyz, do_coord_g96, do_coord_off/),&
            enum_desc=s2a("Coordinates provided through a PDB file format",&
                          "Coordinates provided through an XYZ file format",&
                          "Coordinates provided through a GROMOS96 file format",&
                          "Coordinates read in the &COORD section of the input file"),&
            default_i_val=do_coord_off,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL connectivity_framework(section, do_conn_generate, error)

       CALL keyword_create(keyword, name="EXCLUDE_VDW",&
            description="Specifies which kind of Van der Waals interaction to skip.",&
            usage="EXCLUDE_VDW (1-2||1-3||1-4)", &
            enum_c_vals=s2a( "1-2", "1-3", "1-4"),&
            enum_i_vals=(/do_skip_12, do_skip_13, do_skip_14/),&
            default_i_val=do_skip_13,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXCLUDE_EI",&
            description="Specifies which kind of Electrostatic interaction to skip.",&
            usage="EXCLUDE_EI (1-2||1-3||1-4)", &
            enum_c_vals=s2a( "1-2", "1-3", "1-4"),&
            enum_i_vals=(/do_skip_12, do_skip_13, do_skip_14/),&
            default_i_val=do_skip_13,error=error)
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

 
       NULLIFY(subsection, subsubsection)
       CALL section_create(subsection,name="GENERATE",&
            description="Setup of keywords controlling the generation of the connectivity",&
            n_keywords=2, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="REORDER",&
            description="Reorder a list of atomic coordinates into order so it can be packed correctly.",&
            usage="REORDER <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CREATE_MOLECULES",&
            description="Create molecules names and definition. Can be used to override the "//&
            " molecules specifications of a possible input connectivity or to create molecules"//&
            " specifications for file types as XYZ, missing of molecules definitions.",&
            usage="CREATE_MOLECULES <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
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
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BONDPARM_FACTOR",&
            description="Used in conjunction with BONDPARM to help "//&
                        "determine wheather there is bonding between "//&
                        "two atoms based on a distance criteria.",&
            usage="bondparm_factor {real}", default_r_val=1.1_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
             keyword=keyword,&
             name="BONDLENGTH_MAX",&
             description="Maximum distance to generate neighbor lists to build connectivity",&
             usage="BONDLENGTH_MAX <real>",&
             default_r_val=cp_unit_to_cp2k(value=3.0_dp,&
                                           unit_str="angstrom",&
                                           error=error),&
             unit_str="angstrom",&
             error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
             keyword=keyword,&
             name="BONDLENGTH_MIN",&
             description="Minimum distance to generate neighbor lists to build connectivity",&
             usage="BONDLENGTH_MIN <real>",&
             default_r_val=cp_unit_to_cp2k(value=0.01_dp,&
                                           unit_str="angstrom",&
                                           error=error),&
             unit_str="angstrom",&
             error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsubsection)
       ! BONDS
       CALL section_create(subsubsection,name="BOND",&
            description="Section used to add/remove  bonds in the connectivity."//&
            " Useful for systems with a complex connectivity, difficult to find out automatically.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the bond",&
            usage="&OT (ADD|REMOVE)",&
            enum_c_vals=s2a("ADD","REMOVE"),&
            enum_i_vals=(/do_add,do_remove/),&
            default_i_val=do_add,&
            error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies two atomic index united by a covalent bond",&
            usage="ATOMS {integer} {integer}", type_of_var=integer_t, n_var=2,&
            required=.TRUE., repeats=.TRUE.,error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)      
       
       CALL section_add_subsection(subsection, subsubsection, error=error)
       CALL section_release(subsubsection,error=error)
       NULLIFY(subsubsection)
       ! ANGLES
       CALL section_create(subsubsection,name="ANGLE",&
            description="Section used to add/remove angles in the connectivity."//&
            " Useful for systems with a complex connectivity, difficult to find out automatically.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the bond",&
            usage="&OT (ADD|REMOVE)",&
            enum_c_vals=s2a("ADD","REMOVE"),&
            enum_i_vals=(/do_add,do_remove/),&
            default_i_val=do_add,&
            error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies two atomic index united by a covalent bond",&
            usage="ATOMS {integer} {integer} {integer} ", type_of_var=integer_t, n_var=3,&
            required=.TRUE., repeats=.TRUE.,error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)      
       
       CALL section_add_subsection(subsection, subsubsection, error=error)
       CALL section_release(subsubsection,error=error)
       NULLIFY(subsubsection)
       ! TORSIONS
       CALL section_create(subsubsection,name="TORSION",&
            description="Section used to add/remove torsion in the connectivity."//&
            " Useful for systems with a complex connectivity, difficult to find out automatically.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the bond",&
            usage="&OT (ADD|REMOVE)",&
            enum_c_vals=s2a("ADD","REMOVE"),&
            enum_i_vals=(/do_add,do_remove/),&
            default_i_val=do_add,&
            error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies two atomic index united by a covalent bond",&
            usage="ATOMS {integer} {integer} {integer} {integer} ", type_of_var=integer_t, n_var=4,&
            required=.TRUE., repeats=.TRUE.,error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)      
       
       CALL section_add_subsection(subsection, subsubsection, error=error)
       CALL section_release(subsubsection,error=error)
       NULLIFY(subsubsection)
       ! IMPROPERS
       CALL section_create(subsubsection,name="IMPROPER",&
            description="Section used to add/remove improper in the connectivity."//&
            " Useful for systems with a complex connectivity, difficult to find out automatically.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the bond",&
            usage="&OT (ADD|REMOVE)",&
            enum_c_vals=s2a("ADD","REMOVE"),&
            enum_i_vals=(/do_add,do_remove/),&
            default_i_val=do_add,&
            error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
       
       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies two atomic index united by a covalent bond",&
            usage="ATOMS {integer} {integer} {integer} {integer} ", type_of_var=integer_t, n_var=4,&
            required=.TRUE., repeats=.TRUE.,error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)      
       
       CALL section_add_subsection(subsection, subsubsection, error=error)
       CALL section_release(subsubsection,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)


       NULLIFY(subsection, subsubsection, subsubsubsection)
       CALL section_create(subsection,name="MOL_SET",&
            description="Specify the connectivity of a full system specifyin the connectivity"//&
            " of fragments of the system",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL section_create(subsubsection,name="MOLECULE",&
            description="Specify information about the connectivity of single molecules",&
            n_keywords=2, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="NMOL",&
            description="number of molecules ",&
            usage="NMOL {integer}", default_i_val=1,error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL connectivity_framework(subsubsection, do_conn_psf, error)

       CALL section_add_subsection(subsection, subsubsection, error=error)
       CALL section_release(subsubsection,error=error)

       CALL section_create(subsubsection,name="MERGE_MOLECULES",&
            description="Enables the creation of connecting bridges (bonds, angles, torsions, impropers)"//&
            " between the two or more molecules defined with independent connectivity.",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL section_create(subsubsubsection,name="bonds",&
            description="Defines new bonds",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)      
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Two integer indexes per line defining the new bond."//&
            " Indexes must be relative to the full system and not to the single molecules",&
            repeats=.TRUE.,required=.TRUE.,&
            usage="<INTEGER> <INTEGER>", type_of_var=integer_t, n_var=2, error=error)
       CALL section_add_keyword(subsubsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
       CALL section_add_subsection(subsubsection, subsubsubsection, error=error)
       CALL section_release(subsubsubsection,error=error)

       CALL section_create(subsubsubsection,name="angles",&
            description="Defines new angles",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)      
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Three integer indexes per line defining the new angle"//&
            " Indexes must be relative to the full system and not to the single molecules",repeats=.TRUE.,&
            usage="<INTEGER> <INTEGER> <INTEGER>", type_of_var=integer_t, n_var=3, error=error)
       CALL section_add_keyword(subsubsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
       CALL section_add_subsection(subsubsection, subsubsubsection, error=error)
       CALL section_release(subsubsubsection,error=error)

       CALL section_create(subsubsubsection,name="torsions",&
            description="Defines new torsions",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)      
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Four integer indexes per line defining the new torsion"//&
            " Indexes must be relative to the full system and not to the single molecules",repeats=.TRUE.,&
            usage="<INTEGER> <INTEGER> <INTEGER> <INTEGER>", type_of_var=integer_t, n_var=4, error=error)
       CALL section_add_keyword(subsubsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
       CALL section_add_subsection(subsubsection, subsubsubsection, error=error)
       CALL section_release(subsubsubsection,error=error)

       CALL section_create(subsubsubsection,name="impropers",&
            description="Defines new impropers",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)      
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Four integer indexes per line defining the new improper"//&
            " Indexes must be relative to the full system and not to the single molecules",repeats=.TRUE.,&
            usage="<INTEGER> <INTEGER> <INTEGER> <INTEGER>", type_of_var=integer_t, n_var=4, error=error)
       CALL section_add_keyword(subsubsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
       CALL section_add_subsection(subsubsection, subsubsubsection, error=error)
       CALL section_release(subsubsubsection,error=error)

       CALL section_add_subsection(subsection, subsubsection, error=error)
       CALL section_release(subsubsection,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_topology_section

!***************************************************************************
!!****f* input_cp2k_subsys/connectivity_framework [1.0] *
!!
!!   NAME
!!     connectivity_framework
!!
!!   FUNCTION
!!     Specify keywords used to define connectivity 
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
  SUBROUTINE connectivity_framework(section,default,error)
    TYPE(section_type), POINTER              :: section
    INTEGER, INTENT(IN)                      :: default
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'connectivity_framework', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="CONN_FILE_NAME",&
            variants=(/"CONN_FILE"/),&
            description="Specifies the filename that contains the molecular connectivity.",&
            usage="CONN_FILE_NAME FILENAME",required=.TRUE.,type_of_var=lchar_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="CONNECTIVITY",&
            description="Ways to determine and generate a molecules. "// &
            "Default is to use GENERATE",&
            usage="CONNECTIVITY (PSF||UPSF||MOL_SET||GENERATE||OFF||G87||G96)", &
            enum_c_vals=s2a("PSF","UPSF","MOL_SET","GENERATE","OFF","G87","G96"),&
            enum_i_vals=(/do_conn_psf,&
                          do_conn_psf_u,&
                          do_conn_mol_set,&
                          do_conn_generate,&
                          do_conn_off,&
                          do_conn_g87,do_conn_g96/),&
            enum_desc=s2a("Use  a PSF file to determine the connectivity."//&
                                  " (support standard CHARMM/XPLOR and EXT CHARMM)",&
                          "Read a PSF file in an unformatted way (useful for not so standard PSF).",&
                          "Use multiple PSF (for now...) files to generate the whole sytem.",&
                          "Use a simple distance criteria. (Look at keyword BONDPARM)",&
                          "Do not generate molecules. (Use for QS or ill defined systems)",&
                          "Use GROMOS G87 topology file.",&
                          "Use GROMOS G96 topology file."),&
            default_i_val=default,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE connectivity_framework

END MODULE input_cp2k_subsys
