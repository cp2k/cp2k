!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief builds the subsystem section of the input
!> \par History
!>      10.2005 split input_cp2k [fawzi]
!> \author teo & fawzi
! *****************************************************************************
MODULE input_cp2k_subsys
  USE bibliography,                    ONLY: Goedecker1996,&
                                             Hartwigsen1998,&
                                             Krack2005,&
                                             VandeVondele2005a,&
                                             VandeVondele2007
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE input_constants
  USE input_cp2k_colvar,               ONLY: create_colvar_section
  USE input_cp2k_mm,                   ONLY: create_neighbor_lists_section
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
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE physcon,                         ONLY: bohr
  USE string_utilities,                ONLY: s2a
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_subsys'

PUBLIC :: create_subsys_section,&
          create_cell_section,&
          create_structure_data_section,&
          create_rng_section
!***
CONTAINS

! *****************************************************************************
!> \brief creates the cell section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  RECURSIVE SUBROUTINE create_cell_section(section,periodic,label,error)
    TYPE(section_type), POINTER              :: section
    INTEGER, INTENT(IN), OPTIONAL            :: periodic
    CHARACTER(LEN=*), OPTIONAL               :: label
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_cell_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: my_label
    INTEGER                                  :: my_periodic
    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    my_periodic = use_perd_xyz
    my_label    = "CELL"
    IF (PRESENT(periodic)) my_periodic = periodic
    IF (PRESENT(label))    my_label = label
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,TRIM(my_label),&
            description="Input parameters needed to set up the "//TRIM(my_label)//".",&
            n_keywords=6, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, subsection)
       CALL keyword_create(keyword, name="A",&
            description="Specify the Cartesian components for the cell vector A. "//&
                        "This defines the first column of the h matrix.",&
            usage="A  10.000  0.000  0.000",unit_str="angstrom",&
            n_var=3,type_of_var=real_t,repeats=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="B",&
            description="Specify the Cartesian components for the cell vector B. "//&
                        "This defines the second column of the h matrix.",&
            usage="B   0.000 10.000  0.000", unit_str="angstrom",&
            n_var=3,type_of_var=real_t,repeats=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="C",&
            description="Specify the Cartesian components for the cell vector C. "//&
                        "This defines the third column of the h matrix.",&
            usage="C   0.000  0.000 10.000", unit_str="angstrom",&
            n_var=3,type_of_var=real_t,repeats=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ABC",&
            description="Specify the lengths of the cell vectors A, B, and C, which"//&
            " defines the diagonal elements of h matrix for an orthorhombic cell."//&
            " For non-orthorhombic cells it is possible either to specify the angles "//&
            "ALPHA, BETA, GAMMA via ALPHA_BETA_GAMMA keyword or alternatively use the keywords "//&
            "A, B, and C. The convention is that A lies along the X-axis, B is in the XY plane.",&
            usage="ABC 10.000 10.000 10.000", unit_str="angstrom",&
            n_var=3,type_of_var=real_t,repeats=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALPHA_BETA_GAMMA",&
            variants=(/"ANGLES"/),&
            description="Specify the angles between the vectors A, B and C when using the ABC keyword. "//&
            "The convention is that A lies along the X-axis, B is in the XY plane. "//&
            "ALPHA is the angle between B and C, BETA is the angle between A and C and "//&
            "GAMMA the angle between A and B.",&
            usage="ALPHA_BETA_GAMMA [deg] 90.0 90.0 120.0", unit_str="deg",&
            n_var=3,default_r_vals=(/cp_unit_to_cp2k(value=90.0_dp,unit_str="deg",error=error),&
                                     cp_unit_to_cp2k(value=90.0_dp,unit_str="deg",error=error),&
                                     cp_unit_to_cp2k(value=90.0_dp,unit_str="deg",error=error)/),&
                                     repeats=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CELL_FILE_NAME",&
            description="Possibility to read the cell from an external file ",&
            repeats=.FALSE., usage="CELL_FILE_NAME <CHARACTER>",required=.TRUE.,&
            type_of_var=lchar_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CELL_FILE_FORMAT",&
            description="Specify the format of the cell file (if used)",&
            usage="CELL_FILE_FORMAT (CP2K|XSC)", &
            enum_c_vals=s2a("CP2K","XSC"),&
            enum_i_vals=(/do_cell_cp2k,do_cell_xsc/),&
            enum_desc=s2a("Cell info in the CP2K native format.",&
                          "Cell info in the XSC format (NAMD)" ),&
            default_i_val=do_cell_cp2k,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PERIODIC",&
            description="Specify the directions for which periodic boundary conditions (PBC) will be applied. "//&
                        "Important notice: This applies to the generation of the pair lists as well as to the "//&
                        "application of the PBCs to positions. "//&
                        "See the POISSON section to specify the periodicity used for the electrostatics. "//&
                        "Typically the settings should be the same.",&
            usage="PERIODIC (x|y|z|xy|xz|yz|xyz|none)",&
            enum_c_vals=s2a( "x","y","z","xy","xz","yz","xyz","none"),&
            enum_i_vals=(/ use_perd_x,  use_perd_y,   use_perd_z,&
                           use_perd_xy, use_perd_xz, use_perd_yz,&
                           use_perd_xyz, use_perd_none /),&
            default_i_val=my_periodic, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MULTIPLE_UNIT_CELL",&
            description="Specifies the numbers of repetition in space (X, Y, Z) of the defined cell, "//&
            "assuming it as a unit cell. This keyword affects only the CELL specification. The same keyword "//&
            "in SUBSYS%TOPOLOGY%MULTIPLE_UNIT_CELL should be modified in order to affect the coordinates "//&
            "specification.", usage="MULTIPLE_UNIT_CELL 1 1 1", &
            n_var=3,default_i_vals=(/1,1,1/),repeats=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       IF (TRIM(my_label)=="CELL") THEN
          CALL create_cell_section(subsection,periodic,"CELL_REF",error=error)
          CALL section_add_subsection(section,subsection,error=error)
          CALL section_release(subsection,error=error)
       END IF
    END IF
  END SUBROUTINE create_cell_section

! *****************************************************************************
!> \brief Creates the random number restart section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
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

! *****************************************************************************
!> \brief creates the structure of a subsys, i.e. a full set of
!>      atoms+mol+bounds+cell
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE create_subsys_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_subsys_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="subsys",&
            description="a subsystem: coordinates, topology, molecules and cell",&
            n_keywords=0, n_subsections=9, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)

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

       CALL create_colvar_section(section=subsection,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_multipole_section(subsection, error=error)
       CALL section_add_subsection(section,subsection,error=error)
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

       CALL create_subsys_print_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
    END IF

  END SUBROUTINE create_subsys_section

! *****************************************************************************
!> \brief Creates the subsys print section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_subsys_print_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_subsys_print_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    NULLIFY(print_key, keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Controls printings related to the subsys",&
            n_keywords=0, n_subsections=9, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL cp_print_key_section_create(print_key,"atomic_coordinates",&
            description="controls the output of the atomic coordinates when setting up the"//&
            "force environment. For printing coordinates during MD or GEO refer to the keyword"//&
            " trajectory.",unit_str="angstrom",&
            print_level=medium_print_level, filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL create_structure_data_section(print_key, error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"interatomic_distances",&
            description="controls the output of the interatomic distances when setting up the"//&
            "force environment",unit_str="angstrom",&
            print_level=debug_print_level, filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key, "topology_info", description=&
            "controls the printing of information in the topology settings", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL keyword_create(keyword,"xtl_info",&
            description="Prints information when parsing XTL files.",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"cif_info",&
            description="Prints information when parsing CIF files.",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"pdb_info",&
            description="Prints information when parsing PDB files.",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"xyz_info",&
            description="Prints information when parsing XYZ files.",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"psf_info",&
            description="Prints information when parsing PSF files.",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"amber_info",&
            description="Prints information when parsing ABER topology files.",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"g96_info",&
            description="Prints information when parsing G96 files.",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"crd_info",&
            description="Prints information when parsing CRD files.",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"gtop_info",&
            description="Prints information when parsing GROMOS topology files.",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"util_info",&
            description="Prints information regarding topology utilities",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"generate_info",&
            description="Prints information regarding topology generation",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"cell",&
            description="controls the output of the cell parameters",&
            print_level=medium_print_level, filename="__STD_OUT__",&
            unit_str="angstrom",error=error)
       CALL section_add_subsection(section,print_key,error=error)
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
       CALL keyword_create(keyword, name="se_parameters",&
            description="If the printkey is activated controls the printing of the semi-empirical parameters.",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SYMMETRY",&
            description="controls the output of symmetry information",&
            print_level=debug_print_level+1, filename="__STD_OUT__",error=error)
       CALL keyword_create(keyword, name="MOLECULE",&
            description="Assume the system is an isolated molecule",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="EPS_GEO",&
            description="Accuracy required for symmetry detection",&
            default_r_val=1.e-4_dp, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="STANDARD_ORIENTATION",&
            description="Print molecular coordinates in standard orientation",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="INERTIA",&
            description="Print molecular inertia tensor",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="SYMMETRY_ELEMENTS",&
            description="Print symmetry elements",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="ALL",&
            description="Print all symmetry information",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="ROTATION_MATRICES",&
            description="All the rotation matrices of the point group",&
            default_l_val=.FALSE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="CHECK_SYMMETRY",&
            description="Check if calculated symmetry has expected value."//&
            " Use either Schoenfliess or Hermann-Maugin symbols",&
            default_c_val="NONE", error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"molecules",&
            description="controls the output of information on the molecules",&
            print_level=medium_print_level, filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"radii",&
            description="controls the output of radii information",unit_str="angstrom",&
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

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_subsys_print_section

! *****************************************************************************
!> \brief Creates the multipole section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_multipole_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_multipole_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="multipoles",&
            description="Specifies the dipoles and quadrupoles for particles.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, subsection)
       CALL section_create(subsection,name="dipoles",&
            description="Specifies the dipoles of the particles.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="The dipole components for each atom.",repeats=.TRUE.,&
            usage="dx dy dz", type_of_var=real_t, n_var=3, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="quadrupoles",&
            description="Specifies the quadrupoles of the particles.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="The quadrupole components for each atom.",repeats=.TRUE.,&
            usage="qxx qxy qxz qyy qyz qzz", type_of_var=real_t, n_var=6, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_multipole_section

! *****************************************************************************
!> \brief creates structure data section for output.. both subsys  (for initialization)
!>      and motion section..
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE create_structure_data_section(print_key,error)
    TYPE(section_type), POINTER              :: print_key
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_structure_data_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(print_key),cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN

       NULLIFY(keyword)

       CALL cp_print_key_section_create(print_key,name="STRUCTURE_DATA",&
            description="Request the printing of special structure data during a structure "//&
            "optimization (in MOTION%PRINT) or when setting up a subsys (in SUBSYS%PRINT).",&
            print_level=high_print_level,filename="__STD_OUT__",unit_str="angstrom",error=error)

       CALL keyword_create(keyword, name="POSITION", variants=(/"POS"/),&
            description="Print the position vector of atom a specified by its index",&
            usage="POSITION {integer}",n_var=1,repeats=.TRUE.,required=.FALSE.,&
            type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="DISTANCE",variants=(/"DIS"/),&
            description="Print the distance between the atoms a and b specified by their indices",&
            usage="DISTANCE {integer} {integer}",n_var=2,repeats=.TRUE.,required=.FALSE.,&
            type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="ANGLE",variants=(/"ANG"/),&
            description="Print the angle formed by the atoms specified by their indices",&
            usage="ANGLE {integer} {integer} {integer}",n_var=3, repeats=.TRUE.,&
            required=.FALSE., type_of_var=integer_t, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DIHEDRAL_ANGLE", variants=s2a("DIHEDRAL","DIH"),&
            description="Print the dihedral angle between the planes defined by the atoms (a,b,c) and "//&
            "the atoms (b,c,d) specified by their indices",&
            usage="DIHEDRAL_ANGLE {integer}  {integer} {integer} {integer}",n_var=4,&
            repeats=.TRUE.,required=.FALSE.,type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF

  END SUBROUTINE create_structure_data_section

! *****************************************************************************
!> \brief Creates the velocity section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
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

! *****************************************************************************
!> \brief Creates the shell velocity section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
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

! *****************************************************************************
!> \brief Creates the shell velocity section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
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

! *****************************************************************************
!> \brief Creates the &POTENTIAL section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
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

! *****************************************************************************
!> \brief Creates the &BASIS section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
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

! *****************************************************************************
!> \brief Creates the &GEMINAL section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_geminal_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_geminal_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    IF (.NOT. failure) THEN
       CALL section_create(section,name="geminal",&
            description="Section used to specify a geminal basis set for QM calculations.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="CP2K Basis Set Standard Format",repeats=.TRUE.,&
            type_of_var=lchar_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_geminal_section

! *****************************************************************************
!> \brief Creates the &COORD section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
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
       CALL keyword_create(keyword, name="UNIT",&
             description='Specify the unit of measurement for the coordinates in input'//&
             "All available CP2K units can be used.",&
             usage="UNIT angstrom",default_c_val="angstrom",supported_feature=.TRUE.,&
             error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="SCALED",&
             description='Specify if the coordinateds in input are scaled.',&
             usage="SCALED F",default_l_val=.FALSE.,supported_feature=.TRUE.,&
             lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="ATOM_KIND  X  Y  Z  <MOLNAME>. MOLNAME is optional. "//&
            "If not provided the molecule name is internally created. All other "//&
            "fields after MOLNAME are simply ignored.",repeats=.TRUE.,&
            usage="H 0.7 0.3 0.3", type_of_var=lchar_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_coord_section

! *****************************************************************************
!> \brief Creates the &SHELL_COORD section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
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
       CALL keyword_create(keyword, name="UNIT",&
             description='Specify the unit of measurement for the coordinates in input'//&
             "All available CP2K units can be used.",&
             usage="UNIT angstrom",default_c_val="angstrom",supported_feature=.TRUE.,&
             error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="SCALED",&
             description='Specify if the coordinateds in input are scaled.',&
             usage="SCALED F",default_l_val=.FALSE.,supported_feature=.TRUE.,&
             lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="The shell coordinates in xyz format.",repeats=.TRUE.,&
            usage="at_name rx ry rz at_index", type_of_var=lchar_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_shell_coord_section

! *****************************************************************************
!> \brief Creates the &core_COORD section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
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
       CALL keyword_create(keyword, name="UNIT",&
             description='Specify the unit of measurement for the coordinates in input'//&
             "All available CP2K units can be used.",&
             usage="UNIT angstrom",default_c_val="angstrom",supported_feature=.TRUE.,&
             error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="SCALED",&
             description='Specify if the coordinateds in input are scaled.',&
             usage="SCALED F",default_l_val=.FALSE.,supported_feature=.TRUE.,&
             lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="The core coordinates in xyz format.",repeats=.TRUE.,&
            usage="at_name rx ry rz at_index", type_of_var=lchar_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_core_coord_section

! *****************************************************************************
!> \brief Creates the QM/MM section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
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
            n_keywords=18, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
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

       CALL keyword_create(keyword, name="AUX_FIT_BASIS_SET",&
            variants=s2a("AUXILIARY_FIT_BASIS_SET", "AUX_FIT_BASIS"),&
            description="The auxliliary basis set (GTO type) for wfn fitting",&
            usage="AUX_FIT_BASIS_SET DZVP", default_c_val=" ", &
            n_var=1, required=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="AUX_BASIS_FIT_NORMALIZATION",&
            variants=s2a( "AUXILIARY_BASIS_FIT_NORMALIZATION", "AUX_BASIS_FIT_NORM"),&
            description="The normalization of the auxliliary basis set for wfn fitting",&
            usage="AUX_BASIS_FIT_NORMALIZATION NO", &
            enum_c_vals=s2a( "NONE","NO","WFN","WAVEFUNCTION","FUNCTION",&
            "DENS","DENSITY","SQUARE","UNDEFINED"),&
            enum_i_vals=(/0,0,1,1,1,2,2,2,-1/),&
            default_i_val=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BASIS_SET",&
            variants=s2a("ORBITAL_BASIS_SET","ORB_BASIS"),&
            description="The primary Gaussian basis set (NONE implies no basis used, meaningful with GHOST)",&
            usage="BASIS_SET DZVP", default_c_val=" ", &
            citations=(/VandeVondele2005a,VandeVondele2007/),&
            n_var=1, required=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GEMINAL_BASIS_SET",&
            description="The Geminal Gaussian basis set to be used in HF exchange fitting",&
            usage="GEMINAL_BASIS_SET NAME", default_c_val=" ", &
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
            usage="POTENTIAL <PSEUDO-POTENTIAL-NAME>", default_c_val="GTH", n_var=1,&
            citations=(/Goedecker1996, Hartwigsen1998, Krack2005/),&
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

       CALL keyword_create(keyword, name="ALPHA_SCP",&
            description="The polarizability for scalar-isotropic polarization "//&
            "using SCP with FIST as the driver",&
            usage="ALPHA_SCP 0.2", default_r_val=1.0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="I_SCP",&
            description="The dispersion parameter for scalar-isotropic polarization "//&
            "using SCP with FIST as the driver",&
            usage="I_SCP 0.2", default_r_val=1.0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Logicals
       CALL keyword_create(keyword, name="SE_P_ORBITALS_ON_H",&
          description="Forces the usage of p-orbitals on H for SEMI-EMPIRICAL calculations. "//&
          " This keyword applies only when the KIND is specifying an Hydrogen element. In all "//&
          " other cases is simply ignored. ",&
          usage="SE_P_ORBITALS_ON_H",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GPW_TYPE",&
          description="Force one type to be treated by the GPW scheme,"//&
          " whatever are its primitives, even if the GAPW method is used",&
          usage="GPW_TYPE",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,&
                           name="GHOST",&
                           description="This keyword makes all atoms of this kind "//&
                                        "ghost atoms, i.e. without pseudo or nuclear charge."//&
                                       "Useful to just have the basis set at that position (BSSE calculations),"//&
                                       "or to have a non-interacting particle with BASIS_SET NONE",&
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

       CALL create_geminal_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_potential_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_dft_plus_u_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_kind_section

! *****************************************************************************
!> \brief Create the topology section for FIST.. and the base is running running...
!>      Contains all information regarding topology to be read in input file..
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_topology_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_topology_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

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

       CALL keyword_create(keyword, name="CHARGE_EXTENDED",&
            description="Read MM charges from the very last field of PDB file (starting from column 81)."//&
            " No limitations of number of digits.",&
            usage="CHARGE_EXTENDED logical",&
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
            usage="COORD_FILE_NAME <FILENAME>",required=.TRUE.,type_of_var=lchar_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COORD_FILE_FORMAT",&
            variants=s2a("COORDINATE"),&
            description="Set up the way in which coordinates will be read.",&
            usage="COORD_FILE_FORMAT (OFF|PDB|XYZ|G96|CRD|CIF|XTL)", &
            enum_c_vals=s2a( "OFF","PDB","XYZ","G96","CRD","CIF","XTL"),&
            enum_i_vals=(/do_coord_off, do_coord_pdb, do_coord_xyz, do_coord_g96, do_coord_crd,&
            do_coord_cif,do_coord_xtl/),&
            enum_desc=s2a(&
            "Coordinates read in the &COORD section of the input file",&
            "Coordinates provided through a PDB file format",&
            "Coordinates provided through an XYZ file format",&
            "Coordinates provided through a GROMOS96 file format",&
            "Coordinates provided through an AMBER file format",&
            "Coordinates provided through a CIF (Crystallographic Information File) file format",&
            "Coordinates provided through a XTL (MSI native) file format"),&
            default_i_val=do_coord_off,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL connectivity_framework(section, do_conn_generate, error)

       CALL keyword_create(keyword, name="EXCLUDE_VDW",&
            description="Specifies which kind of Van der Waals interaction to skip.",&
            usage="EXCLUDE_VDW (1-3||1-4)", &
            enum_c_vals=s2a( "1-3", "1-4"),&
            enum_i_vals=(/do_skip_13, do_skip_14/),&
            default_i_val=do_skip_13,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXCLUDE_EI",&
            description="Specifies which kind of Electrostatic interaction to skip.",&
            usage="EXCLUDE_EI (1-3||1-4)", &
            enum_c_vals=s2a( "1-3", "1-4"),&
            enum_i_vals=(/do_skip_13, do_skip_14/),&
            default_i_val=do_skip_13,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MULTIPLE_UNIT_CELL",&
            description="Specifies the numbers of repetition in space (X, Y, Z) of the defined cell, "//&
            "assuming it as a unit cell. This keyword affects only the coordinates specification. The same keyword "//&
            "in SUBSYS%CELL%MULTIPLE_UNIT_CELL should be modified in order to affect the cell "//&
            "specification.", usage="MULTIPLE_UNIT_CELL 1 1 1", &
            n_var=3,default_i_vals=(/1,1,1/),repeats=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,"DUMP_PDB",&
            description="controls the dumping of the PDB at the starting geometry",&
            print_level=debug_print_level, filename="dump",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DUMP_PSF",&
            description="controls the dumping of the PSF connectivity",&
            print_level=debug_print_level, filename="dump",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       NULLIFY(subsection)
       CALL create_center_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection, error=error)

       CALL create_generate_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_molset_section(subsection, error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
    END IF

  END SUBROUTINE create_topology_section

! *****************************************************************************
!> \brief Specify keywords used to center molecule in the box
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - University of Zurich - 06.2009
! *****************************************************************************
  SUBROUTINE create_center_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_center_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(keyword)
       CALL section_create(section,"CENTER_COORDINATES",&
            description="Allows centering the coordinates of the system in the box. "//&
            "The centering point can be defined by the user.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="Controls the activation of the centering method",&
            usage="&CENTER_COORDINATES T",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CENTER_POINT",&
            description="Specify the point used for centering the coordinates. Default is to "//&
            "center the system in cell/2. ", type_of_var=real_t, n_var=3,&
            repeats=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_center_section

! *****************************************************************************
!> \brief Specify keywords used to setup several molecules with few connectivity files
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - University of Zurich - 08.2008
! *****************************************************************************
  SUBROUTINE create_molset_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_molset_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection, subsubsection

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(keyword, subsection, subsubsection)
       CALL section_create(section,name="MOL_SET",&
            description="Specify the connectivity of a full system specifyin the connectivity"//&
            " of fragments of the system",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       ! MOLECULES
       CALL section_create(subsection,name="MOLECULE",&
            description="Specify information about the connectivity of single molecules",&
            n_keywords=2, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="NMOL",&
            description="number of molecules ",&
            usage="NMOL {integer}", default_i_val=1,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL connectivity_framework(subsection, do_conn_psf, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! MERGE MOLECULES
       CALL section_create(subsection,name="MERGE_MOLECULES",&
            description="Enables the creation of connecting bridges (bonds, angles, torsions, impropers)"//&
            " between the two or more molecules defined with independent connectivity.",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL section_create(subsubsection,name="bonds",&
            description="Defines new bonds",n_keywords=2, n_subsections=0, repeats=.FALSE.,&
            required=.TRUE.,error=error)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Two integer indexes per line defining the new bond."//&
            " Indexes must be relative to the full system and not to the single molecules",&
            repeats=.TRUE.,required=.TRUE.,&
            usage="<INTEGER> <INTEGER>", type_of_var=integer_t, n_var=2, error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(subsection, subsubsection, error=error)
       CALL section_release(subsubsection,error=error)

       CALL section_create(subsubsection,name="angles",&
            description="Defines new angles",n_keywords=2, n_subsections=0,&
            repeats=.FALSE., required=.TRUE.,error=error)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Three integer indexes per line defining the new angle"//&
            " Indexes must be relative to the full system and not to the single molecules",repeats=.TRUE.,&
            usage="<INTEGER> <INTEGER> <INTEGER>", type_of_var=integer_t, n_var=3, error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(subsection, subsubsection, error=error)
       CALL section_release(subsubsection,error=error)

       CALL section_create(subsubsection,name="torsions",&
            description="Defines new torsions",n_keywords=2, n_subsections=0,&
            repeats=.FALSE., required=.TRUE.,error=error)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Four integer indexes per line defining the new torsion"//&
            " Indexes must be relative to the full system and not to the single molecules",repeats=.TRUE.,&
            usage="<INTEGER> <INTEGER> <INTEGER> <INTEGER>", type_of_var=integer_t, n_var=4, error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(subsection, subsubsection, error=error)
       CALL section_release(subsubsection,error=error)

       CALL section_create(subsubsection,name="impropers",&
            description="Defines new impropers",n_keywords=2, n_subsections=0,&
            repeats=.FALSE., required=.TRUE.,error=error)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Four integer indexes per line defining the new improper"//&
            " Indexes must be relative to the full system and not to the single molecules",repeats=.TRUE.,&
            usage="<INTEGER> <INTEGER> <INTEGER> <INTEGER>", type_of_var=integer_t, n_var=4, error=error)
       CALL section_add_keyword(subsubsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(subsection, subsubsection, error=error)
       CALL section_release(subsubsection,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_molset_section

! *****************************************************************************
!> \brief Specify keywords used to generate connectivity
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - University of Zurich - 08.2008
! *****************************************************************************
  SUBROUTINE create_generate_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_generate_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(keyword,subsection)
       CALL section_create(section,name="GENERATE",&
            description="Setup of keywords controlling the generation of the connectivity",&
            n_keywords=2, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="REORDER",&
            description="Reorder a list of atomic coordinates into order so it can be packed correctly.",&
            usage="REORDER <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CREATE_MOLECULES",&
            description="Create molecules names and definition. Can be used to override the "//&
            " molecules specifications of a possible input connectivity or to create molecules"//&
            " specifications for file types as XYZ, missing of molecules definitions.",&
            usage="CREATE_MOLECULES <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
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

       CALL keyword_create(keyword, name="BONDPARM_FACTOR",&
            description="Used in conjunction with BONDPARM to help "//&
                        "determine wheather there is bonding between "//&
                        "two atoms based on a distance criteria.",&
            usage="bondparm_factor {real}", default_r_val=1.1_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,name="BONDLENGTH_MAX",&
             description="Maximum distance to generate neighbor lists to build connectivity",&
             usage="BONDLENGTH_MAX <real>",default_r_val=cp_unit_to_cp2k(value=3.0_dp,&
             unit_str="angstrom",error=error),unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,name="BONDLENGTH_MIN",&
             description="Minimum distance to generate neighbor lists to build connectivity",&
             usage="BONDLENGTH_MIN <real>",default_r_val=cp_unit_to_cp2k(value=0.01_dp,&
             unit_str="angstrom",error=error),unit_str="angstrom",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! BONDS
       CALL section_create(subsection,name="BOND",&
            description="Section used to add/remove  bonds in the connectivity."//&
            " Useful for systems with a complex connectivity, difficult to find out automatically.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the bond",&
            usage="&BOND (ADD|REMOVE)",&
            enum_c_vals=s2a("ADD","REMOVE"),&
            enum_i_vals=(/do_add,do_remove/),&
            default_i_val=do_add,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies two atomic index united by a covalent bond",&
            usage="ATOMS {integer} {integer}", type_of_var=integer_t, n_var=2,&
            required=.TRUE., repeats=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! ANGLES
       CALL section_create(subsection,name="ANGLE",&
            description="Section used to add/remove angles in the connectivity."//&
            " Useful for systems with a complex connectivity, difficult to find out automatically.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the bond",&
            usage="&ANGLE (ADD|REMOVE)",&
            enum_c_vals=s2a("ADD","REMOVE"),&
            enum_i_vals=(/do_add,do_remove/),&
            default_i_val=do_add,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies two atomic index united by a covalent bond",&
            usage="ATOMS {integer} {integer} {integer} ", type_of_var=integer_t, n_var=3,&
            required=.TRUE., repeats=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! TORSIONS
       CALL section_create(subsection,name="TORSION",&
            description="Section used to add/remove torsion in the connectivity."//&
            " Useful for systems with a complex connectivity, difficult to find out automatically.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the bond",&
            usage="&TORSION (ADD|REMOVE)",&
            enum_c_vals=s2a("ADD","REMOVE"),&
            enum_i_vals=(/do_add,do_remove/),&
            default_i_val=do_add,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies two atomic index united by a covalent bond",&
            usage="ATOMS {integer} {integer} {integer} {integer} ", type_of_var=integer_t, n_var=4,&
            required=.TRUE., repeats=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! IMPROPERS
       CALL section_create(subsection,name="IMPROPER",&
            description="Section used to add/remove improper in the connectivity."//&
            " Useful for systems with a complex connectivity, difficult to find out automatically.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the bond",&
            usage="&IMPROPER (ADD|REMOVE)",&
            enum_c_vals=s2a("ADD","REMOVE"),&
            enum_i_vals=(/do_add,do_remove/),&
            default_i_val=do_add,&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies two atomic index united by a covalent bond",&
            usage="ATOMS {integer} {integer} {integer} {integer} ", type_of_var=integer_t, n_var=4,&
            required=.TRUE., repeats=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! ISOLATED ATOMS
       CALL section_create(subsection,name="ISOLATED_ATOMS",&
            description=" This section specifies the  atoms that one considers isolated. Useful when present "//&
            " ions in solution.",n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword, name="LIST",&
            description="Specifies a list of atomic indexes of the isolated ion",&
            usage="LIST {integer}", type_of_var=integer_t, n_var=-1,&
            required=.TRUE., repeats=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! Neighbor lists keys and printing handling the construction of NL for the connectivity
       CALL create_neighbor_lists_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_gen_print_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_generate_section

! *****************************************************************************
!> \brief Create the print gen section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_gen_print_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_gen_print_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options in GENERATE code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"NEIGHBOR_LISTS",&
            description="Activates the printing of the neighbor lists used"//&
            " for generating the connectivity.", print_level=high_print_level,&
            filename="", unit_str="angstrom", error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SUBCELL",&
            description="Activates the printing of the subcells used for the"//&
            "generation of neighbor lists for connectivity.", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_gen_print_section

! *****************************************************************************
!> \brief Specify keywords used to define connectivity
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
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
            usage="CONN_FILE_NAME <FILENAME>",required=.TRUE.,type_of_var=lchar_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CONN_FILE_FORMAT",&
            variants=(/"CONNECTIVITY"/),&
            description="Ways to determine and generate a molecules. "// &
            "Default is to use GENERATE",&
            usage="CONN_FILE_FORMAT (PSF|UPSF|MOL_SET|GENERATE|OFF|G87|G96|AMBER)", &
            enum_c_vals=s2a("PSF","UPSF","MOL_SET","GENERATE","OFF","G87","G96","AMBER"),&
            enum_i_vals=(/do_conn_psf,&
                          do_conn_psf_u,&
                          do_conn_mol_set,&
                          do_conn_generate,&
                          do_conn_off,&
                          do_conn_g87,&
                          do_conn_g96,&
                          do_conn_amb7/),&
            enum_desc=s2a("Use  a PSF file to determine the connectivity."//&
                                  " (support standard CHARMM/XPLOR and EXT CHARMM)",&
                          "Read a PSF file in an unformatted way (useful for not so standard PSF).",&
                          "Use multiple PSF (for now...) files to generate the whole sytem.",&
                          "Use a simple distance criteria. (Look at keyword BONDPARM)",&
                          "Do not generate molecules. (Use for QS or ill defined systems)",&
                          "Use GROMOS G87 topology file.",&
                          "Use GROMOS G96 topology file.",&
                          "Use AMBER topology file for reading connectivity (compatible starting from AMBER V.7)"),&
            default_i_val=default,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE connectivity_framework

! *****************************************************************************
!> \brief      Create CP2K input section for the DFT+U method parameters
!> \author     Matthias Krack (MK)
!> \date       01.11.2007
!> \version    1.0
! *****************************************************************************
  SUBROUTINE create_dft_plus_u_section(section,error)

    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_dft_plus_u_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure = .FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

      CALL section_create(section,&
                          name="DFT_PLUS_U",&
                          description="Define the parameters for a DFT+U run",&
                          n_keywords=3,&
                          n_subsections=0,&
                          repeats=.FALSE.,&
                          required=.FALSE.,&
                          error=error)

      NULLIFY (keyword)

      CALL keyword_create(keyword,&
                          name="_SECTION_PARAMETERS_",&
                          description="controls the activation of the DFT+U section",&
                          usage="&DFT_PLUS_U ON",&
                          default_l_val=.FALSE.,&
                          lone_keyword_l_val=.TRUE.,&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="L",&
                          description="Angular momentum quantum number of the"//&
                                      "orbitals to which the correction is applied",&
                          required=.TRUE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=integer_t,&
                          default_i_val=-1,&
                          usage="L 2",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="U_MINUS_J",&
                          description="Effective parameter U(eff) = U - J",&
                          required=.TRUE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=0.0_dp,&
                          unit_str="au_e",&
                          usage="U_MINUS_J [eV] 1.4",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_dft_plus_u_section

END MODULE input_cp2k_subsys
