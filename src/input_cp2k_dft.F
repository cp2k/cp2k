!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief function that build the dft section of the input
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE input_cp2k_dft
  USE bibliography,                    ONLY: &
       Avezac2005, Blochl1995, Brelaz1979, Dewar1977, Dewar1985, Dudarev1997, &
       Dudarev1998, Elstner1998, Guidon2010, Hu2007, Hunt2003, Iannuzzi2005, &
       Iannuzzi2006, Iannuzzi2007, Kolafa2004, Krack2000, Krack2002, &
       Kunert2003, Lippert1997, Lippert1999, Perdew1981, Porezag1995, &
       Repasky2002, Rocha2006, Schenter2008, Seifert1996, Stewart1982, &
       Stewart1989, Stewart2007, Thiel1992, Tozer1996, VandeVondele2003, &
       VandeVondele2005a, VandeVondele2005b, VandeVondele2006, Weber2008, &
       Zhao1994, Zhechkov2005
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE input_constants
  USE input_cp2k_almo,                 ONLY: create_almo_scf_section
  USE input_cp2k_atom,                 ONLY: create_powell_section
  USE input_cp2k_distribution,         ONLY: create_distribution_section
  USE input_cp2k_ls,                   ONLY: create_ls_scf_section
  USE input_cp2k_mm,                   ONLY: create_dipoles_section,&
                                             create_neighbor_lists_section
  USE input_cp2k_motion,               ONLY: add_format_keyword
  USE input_cp2k_poisson,              ONLY: create_poisson_section
  USE input_cp2k_rsgrid,               ONLY: create_rsgrid_section
  USE input_cp2k_xc,                   ONLY: create_xc_section
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
                                             logical_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_dft'

  PUBLIC :: create_dft_section,create_relativistic_section
  PUBLIC :: create_bsse_section, create_qs_section
  PUBLIC :: create_scf_section, create_ep_section
  PUBLIC :: create_interp_section, create_localize_section, print_wanniers

CONTAINS

! *****************************************************************************
!> \brief creates an ep section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE create_ep_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ep_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, print_section

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="ep",&
         description="parameter needed by an ep calculation",&
         repeats=.FALSE., required=.TRUE.,&
         error=error)

    NULLIFY(keyword)
    CALL keyword_create(keyword, name="EPS_LIN_SOLV",&
         description="Requested convergence of the linear solver (for psi1)",&
         usage="EPS_LIN_SOLV",&
         default_r_val=1.e-5_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="ROTATE",&
         description="If rotations from a unique set of coefficients should"//&
         " be used or if single molecule optimizations generate it",&
         usage="ROTATE",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="START_COEFFS",&
         description="Starting coefficients for roatation based ep",&
         usage="START_COEFFS 1.0 0.0 ...",&
          type_of_var=real_t, n_var=-1,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="E0_ONLY",&
         description="If only e0 should be calculated",&
         usage="E0_ONLY",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="COMP_INPUT",&
         description="Path to the input to be used for the component of the main system",&
         usage="COMP_INPUT comp.inp",&
         default_lc_val="comp.inp",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="AT_PER_MOL",&
         description="Number of atoms in each molecule (at the moment only uniform system cam be handled)",&
         usage="at_per_mol 3",&
         default_i_val=3,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    NULLIFY(print_section)
    CALL section_create(print_section,name="print",&
         description="Section of possible print options for an EP calculation.",&
         repeats=.FALSE., required=.TRUE.,error=error)

    NULLIFY(print_key)
    CALL cp_print_key_section_create(print_key,"ENERGY",&
         description="Controls the printing of the various energies",&
         print_level=medium_print_level,filename="",&
         error=error)
    CALL section_add_subsection(print_section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"RUN_INFO",&
         description="Prints various informations on the progress of the ep calculation",&
         print_level=medium_print_level,filename="__STD_OUT__",&
         error=error)
    CALL keyword_create(keyword, name="LIN_SOLV",&
         description="print the linear solver progress",&
         usage="lin_solv",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL section_add_subsection(print_section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"EP_MATRIXES",&
         description="Prints various matrixes of an ep calculation",&
         print_level=debug_print_level,filename="__STD_OUT__",&
         error=error)

    CALL keyword_create(keyword, name="psi0",&
         description="print the value of the psi0 matrix",&
         usage="psi0",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="psi0_blocks",&
         description="print the value of the blocs of the psi0 matrix",&
         usage="psi0_blocks",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="psi1",&
         description="print the value of the psi1 matrix",&
         usage="psi1",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="mo_s_inv",&
         description="print the value of the inverse of the mo overlap matrix",&
         usage="mo_s_inv",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL section_add_subsection(print_section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL cp_print_key_section_create(print_key,"EP_RHO_CUBE",&
         description="Prints out a cube file with the calculated rho",&
         print_level=debug_print_level,filename="EP_RHO",&
         error=error)

    CALL keyword_create(keyword, name="rho0",&
         description="print the cube of rho0,the psi0 derived density matrix",&
         usage="rho0 OFF",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="rho1",&
         description="print the value of the blocs of rho1, the psi1 (and psi0) derived density matrix",&
         usage="rho1 F",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="stride",&
         description="The stride (X,Y,Z) used to write the cube file "//&
         "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
         " 1 number valid for all components.",&
         usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="mo_s_inv",&
         description="print the value of the inverse of the mo overlap matrix",&
         usage="mo_s_inv",&
         default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL section_add_subsection(print_section,print_key,error=error)
    CALL section_release(print_key,error=error)

    CALL section_add_subsection(section,print_section,error=error)
    CALL section_release(print_section,error=error)

  END SUBROUTINE create_ep_section

! *****************************************************************************
!> \brief creates the dft section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE create_dft_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_dft_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="dft",&
         description="parameter needed by dft programs",&
         n_keywords=3, n_subsections=4, repeats=.FALSE., required=.TRUE.,&
         error=error)

    NULLIFY(keyword)
    CALL keyword_create(keyword, name="BASIS_SET_FILE_NAME",&
         description="Name of the basis set file, may include a path",&
         usage="BASIS_SET_FILE_NAME <FILENAME>",&
         type_of_var=lchar_t,repeats=.TRUE.,&
         default_lc_val="BASIS_SET",n_var=1,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="POTENTIAL_FILE_NAME",&
         description="Name of the pseudo potential file, may include a path",&
         usage="POTENTIAL_FILE_NAME <FILENAME>",&
         default_lc_val="POTENTIAL",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="GEMINAL_FILE_NAME",&
         description="Name of the geminal basis set file, may include a path",&
         usage="GEMINAL_FILE_NAME <FILENAME>",&
         default_lc_val="BASIS_GEMINAL",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="WFN_RESTART_FILE_NAME",&
         variants=(/"RESTART_FILE_NAME"/),&
         description="Name of the wavefunction restart file, may include a path."//&
         " If no file is specified, the default is to open the file as generated by the wfn restart print key.",&
         usage="WFN_RESTART_FILE_NAME <FILENAME>",&
         type_of_var=lchar_t,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword=keyword,&
         name="UKS",&
         variants=s2a("UNRESTRICTED_KOHN_SHAM",&
         "LSD",&
         "SPIN_POLARIZED"),&
         description="Requests a spin-polarized calculation using alpha "//&
         "and beta orbitals, i.e. no spin restriction is applied",&
         usage="LSD",&
         default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword=keyword,&
         name="ROKS",&
         variants=(/"RESTRICTED_OPEN_KOHN_SHAM"/),&
         description="Requests a restricted open Kohn-Sham calculation",&
         usage="ROKS",&
         default_l_val=.FALSE.,&
         lone_keyword_l_val=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword=keyword,&
         name="MULTIPLICITY",&
         variants=(/"MULTIP"/),&
         description="Two times the total spin plus one. "//&
         "Specify 3 for a triplet, 4 for a quartet, "//&
         "and so on. Default is 1 (singlet) for an "//&
         "even number and 2 (doublet) for an odd number "//&
         "of electrons.",&
         usage="MULTIPLICITY 3",&
         default_i_val=0,& ! this default value is just a flag to get the above
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="CHARGE",&
         description="The total charge of the system",&
         usage="CHARGE -1",&
         default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="EXCITATIONS",&
         description="If excitations should be calculated",&
         usage="EXCITATIONS",&
         enum_c_vals=s2a("NONE","TDLR","TDDFPT"),&
         enum_i_vals=(/ no_excitations, tddfpt_excitations, &
         tddfpt_excitations/),&
         default_i_val=no_excitations, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword,&
                        name="PLUS_U_METHOD",&
                        description="Method employed for the calculation of the DFT+U contribution",&
                        required=.FALSE.,&
                        repeats=.FALSE.,&
                        enum_c_vals=s2a("LOWDIN","MULLIKEN","MULLIKEN_CHARGES"),&
                        enum_i_vals=(/plus_u_lowdin,plus_u_mulliken,plus_u_mulliken_charges/),&
                        enum_desc=s2a("Method based on Lowdin population analysis "//&
                                      "(computationally expensive, since the diagonalization of the "//&
                                      "overlap matrix is required, but possibly more robust than Mulliken)",&
                                      "Method based on Mulliken population analysis using the net AO and "//&
                                      "overlap populations (computationally cheap method)",&
                                      "Method based on Mulliken gross orbital populations (GOP)"),&
                        n_var=1,&
                        default_i_val=plus_u_mulliken,&
                        usage="METHOD Lowdin",&
                        error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword,&
                        name="RELAX_MULTIPLICITY",&
                        variants=(/"RELAX_MULTIP"/),&
                        description="Do not enforce the occupation of alpha and beta MOs due to the initially "//&
                                    "defined multiplicity, but rather follow an Aufbau principle. "//&
                                    "A threshold value greater than zero activates this option. "//&
                                    "Larger threshold values increase the probability for a spin flip. "//&
                                    "This option is only valid for unrestricted (i.e. spin polarised) "//&
                                    "Kohn-Sham (UKS) calculations.",&
                        usage="RELAX_MULTIPLICITY 0.00001",&
                        required=.FALSE.,&
                        repeats=.FALSE.,&
                        default_r_val=0.0_dp,&
                        error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="SUBCELLS",&
         description="Read the grid size for subcell generation in the construction of "//&
         "neighbor lists.", usage="SUBCELLS 1.5",&
         n_var=1,default_r_val=2.0_dp, error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    NULLIFY(subsection)
    CALL create_scf_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_ls_scf_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_almo_scf_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_kg_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_admm_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_qs_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_tddfpt_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_mgrid_section(subsection,error=error)
    CALL section_add_subsection(section, subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_xc_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_relativistic_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_sic_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_low_spin_roks_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_efield_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_per_efield_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_ext_pot_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)
! ZMP sections to include the external density or v_xc potential
    CALL create_ext_den_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_ext_vxc_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_poisson_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_implicit_solv_section(subsection, error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)

    CALL create_density_fitting_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_xas_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_localize_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_rtp_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_print_dft_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)
  END SUBROUTINE create_dft_section

! *****************************************************************************
!> \brief Implicit Solvation Model
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author tlaino
! *****************************************************************************
  SUBROUTINE create_implicit_solv_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_implicit_solv_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    NULLIFY(keyword, subsection, print_key)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="SCRF",&
            description="Adds an implicit solvation model to the DFT calculation."//&
            " Know also as Self Consistent Reaction Field.",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword,name="EPS_OUT",&
            description="Value of the dielectric constant outside the sphere",&
            usage="EPS_OUT <REAL>",&
            default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="LMAX",&
            description="Maximum value of L used in the multipole expansion",&
            usage="LMAX <INTEGER>",&
            default_i_val=3,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_sphere_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing basic info about the method", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_implicit_solv_section

! *****************************************************************************
!> \brief Create Sphere cavity
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author tlaino
! *****************************************************************************
  SUBROUTINE create_sphere_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_sphere_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    NULLIFY(keyword, subsection)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="SPHERE",&
            description="Treats the implicit solvent environment like a sphere",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword,name="RADIUS",&
            description="Value of the spherical cavity in the dielectric medium",&
            usage="RADIUS <REAL>",&
            unit_str="angstrom",&
            type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_center_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_sphere_section

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author tlaino
! *****************************************************************************
  SUBROUTINE create_center_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_center_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="CENTER",&
            description="Defines the center of the sphere.",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword,name="XYZ",&
            description="Coordinates of the center of the sphere",&
            usage="XYZ <REAL> <REAL> <REAL>",&
            unit_str="angstrom",&
            type_of_var=real_t, n_var=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="ATOM_LIST",&
            description="Defines a list of atoms to define the center of the sphere",&
            usage="ATOM_LIST <INTEGER> .. <INTEGER>",&
            type_of_var=integer_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="WEIGHT",&
            description="Defines the weight used to define the center of the sphere"//&
            " (if ATOM_LIST is provided)",&
            usage="WEIGHT (UNIT|MASS)",&
            enum_c_vals=(/"UNIT","MASS"/),&
            enum_i_vals=(/weight_type_unit,weight_type_mass/),&
            default_i_val=weight_type_unit,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="FIXED",&
            description="Specify if the center of the sphere should be fixed or"//&
            " allowed to move",&
            usage="FIXED <LOGICAL>",&
            default_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_center_section

! *****************************************************************************
!> \brief parameters fo the localization of wavefunctions
!> \par History
!>      03.2005 created [MI]
! *****************************************************************************
  SUBROUTINE create_localize_section(section, error)

    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_localize_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, print_section, &
                                                subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       NULLIFY(keyword, print_key)
       CALL section_create(section,name="LOCALIZE",&
            description="Use one of the available methods to define the localization "//&
            " and possibly to optimize it to a minimum or a maximum.",&
            n_keywords=8, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the MOS localization procedure",&
            usage="&LOCALIZE T",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_ITER",&
            description="Maximum number of iterations used for localization methods",&
            usage="MAX_ITER 2000", default_i_val=10000, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_CRAZY_ANGLE",&
            description="Largest allowed angle for the crazy rotations algorithm (smaller is slower but more stable).",&
            usage="MAX_CRAZY_ANGLE 0.1", unit_str="rad", default_r_val=0.2_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CRAZY_SCALE",&
            description="scale angles",&
            usage="CRAZY_SCALE 0.9", default_r_val=1.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CRAZY_USE_DIAG",&
            description="Use diagonalization (slow) or pade based calculation of matrix exponentials.",&
            usage="CRAZY_USE_DIAG ", default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="USE_HISTORY",&
            description="Generate an improved initial guess based on a history of results, which is useful during MD."//&
                        "Will only work if the number of states to be localized remains constant.",&
            usage="USE_HISTORY ", default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_OCCUPATION",&
            description="Tolerance in the occupation number to select only fully occupied orbitals for the rotation",&
            usage="EPS_OCCUPATION 1.E-5", default_r_val=1.0E-8_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="OUT_ITER_EACH",&
            description="Every how many iterations of the localization algorithm"//&
            "(Jacobi) the tolerance value is printed out",&
            usage="OUT_ITER_EACH 100", default_i_val=100, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_LOCALIZATION",&
            description="Tolerance used in the convergence criterium of the localization methods.",&
            usage="EPS_LOCALIZATION 1.0E-2", default_r_val=1.0E-4_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="MIN_OR_MAX",&
            description="Requires the maximization of the spread of the wfn",&
            usage="MIN_OR_MAX (SPREADMIN|SPREADMAX)",&
            enum_c_vals=(/"SPREADMIN","SPREADMAX"/),&
            enum_i_vals=(/do_loc_min, do_loc_max/),&
            default_i_val=do_loc_min,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="METHOD",&
            description="Method of optimization if any",&
            usage="METHOD (JACOBI|CRAZY|DIRECT|L1SD|NONE)",&
            enum_c_vals=s2a("NONE","JACOBI","CRAZY","L1SD","DIRECT"),&
            enum_i_vals=(/do_loc_none, &
                          do_loc_jacobi, &
                          do_loc_crazy, &
                          do_loc_l1_norm_sd, &
                          do_loc_direct/),&
            enum_desc=s2a("No localization is applied",&
            "Using 2 x 2 rotations of the orbitals, slow but robust",&
            "A new fast method is applied, might be slightly less robust than jacobi, but usually much faster",&
            "Steepest descent minimization of an approximate l1 norm",&
            "Using a direct minimisation approach"),&
            default_i_val=do_loc_jacobi,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="JACOBI_FALLBACK",&
            description="Use Jacobi method in case no convergence was achieved"//&
            " by using the crazy rotations method.",&
            usage="JACOBI_FALLBACK", default_l_val=.TRUE., &
            lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART",&
            description="Restart the localization from a set of orbitals"//&
            " read from a localization restart file.",&
            usage="RESTART", default_l_val=.FALSE., &
            lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LOCHOMO_RESTART_FILE_NAME",&
            description="File name where to read the MOS from"//&
            "which to restart the localization procedure for occupied states",&
            usage="LOCHOMO_RESTART_FILE_NAME <FILENAME>",&
            type_of_var=lchar_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LOCLUMO_RESTART_FILE_NAME",&
            description="File name where to read the MOS from"//&
            "which to restart the localization procedure for unoccupied states",&
            usage="LOCLUMO_RESTART_FILE_NAME <FILENAME>",&
            type_of_var=lchar_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="OPERATOR",&
            description="Type of opertator which defines the spread functional",&
            usage="OPERATOR (BERRY|BOYS|PIPEK)",&
            enum_c_vals=s2a("BERRY","BOYS","PIPEK"),&
            enum_i_vals=(/op_loc_berry, op_loc_boys, op_loc_pipek/),&
            default_i_val=op_loc_berry,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="LIST",&
            description="Indexes of the occupied wfn to be localized"//&
            "This keyword can be repeated several times"//&
            "(useful if you have to specify many indexes).",&
            usage="LIST 1 2",&
            n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="LIST_UNOCCUPIED",&
            description="Indexes of the unoccupied states to be localized, "//&
            "up to now only valid in combination with GPW. "//&
            "This keyword has to be present if unoccupied states should be localized. "//&
            "This keyword can be repeated several times"//&
            "(useful if you have to specify many indexes).",&
            usage="LIST 1 2",&
            n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="STATES",&
            description="Which states to localize, LUMO up to now only available in GPW",&
            usage="STATES (HOMO|LUMO|ALL)",&
            enum_c_vals=s2a("OCCUPIED","UNOCCUPIED","ALL"),&
            enum_i_vals=(/do_loc_homo, do_loc_lumo,do_loc_both/),&
            default_i_val=do_loc_homo,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="ENERGY_RANGE",&
            description="Select the orbitals to be localized within the given energy range."//&
            "This type of selection cannot be added on top of the selection through a LIST. It reads to reals that are"//&
            " lower and higher boundaries of the energy range.",&
            usage=" ENERGY_RANGE lower_bound {real}, higher_bound {real}", &
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=2,default_r_vals=(/0._dp,0._dp/),unit_str='eV',&
            type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(print_section)
       CALL section_create(print_section,name="PRINT",&
            description="Collects all printing options related to the Wannier centers and "//&
            "properties computed with Wannier centers.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)
       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing basic info about the method", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(print_section,print_key,error=error)
       CALL section_release(print_key,error=error)
       ! Add printing of wannier infos
       CALL print_wanniers(print_section, error)
       NULLIFY(subsection)
       ! Total Dipoles with wannier
       CALL create_dipoles_section(subsection,"TOTAL_DIPOLE",debug_print_level+1,error)
       CALL section_add_subsection(print_section, subsection, error=error)
       CALL section_release(subsection,error=error)
       ! Molecular Dipoles with wannier
       CALL create_dipoles_section(subsection,"MOLECULAR_DIPOLES",debug_print_level+1,error)
       CALL section_add_subsection(print_section, subsection, error=error)
       CALL section_release(subsection,error=error)
       ! Molecular States with wannier
       CALL create_molecular_states_section(subsection,error)
       CALL section_add_subsection(print_section, subsection, error=error)
       CALL section_release(subsection,error=error)
       ! Wannier States with wannier
       CALL create_wannier_states_section(subsection,error)
       CALL section_add_subsection(print_section, subsection, error=error)
       CALL section_release(subsection,error=error)
       CALL section_add_subsection(section,print_section,error=error)
       CALL section_release(print_section,error=error)


    END IF
  END SUBROUTINE create_localize_section

! *****************************************************************************
!> \brief Controls the printing of the basic info coming from the LOCALIZE
!>        section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE print_wanniers(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_wanniers', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    NULLIFY(print_key, keyword)
    CALL cp_print_key_section_create(print_key,"WANNIER_CUBES",&
         description="Controls the printing of the wannier functions ", &
         print_level=high_print_level,add_last=add_last_numeric,filename="",&
         error=error)
    CALL keyword_create(keyword, name="stride",&
         description="The stride (X,Y,Z) used to write the cube file "//&
         "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
         " 1 number valid for all components.",&
         usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword,name="CUBES_LU_BOUNDS",&
         variants=(/"CUBES_LU"/),&
         description="The lower and upper index of the states to be printed as cube",&
         usage="CUBES_LU_BOUNDS integer integer",&
         n_var=2,default_i_vals=(/0,-2/), type_of_var=integer_t,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword,name="CUBES_LIST",&
         description="Indexes of the states to be printed as cube files"//&
         "This keyword can be repeated several times"//&
         "(useful if you have to specify many indexes).",&
         usage="CUBES_LIST 1 2",&
         n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="APPEND",&
            description="append the cube files when they already exist",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    NULLIFY(print_key)
    CALL cp_print_key_section_create(print_key,"WANNIER_CENTERS",&
         description="Controls the printing of the wannier functions", &
         print_level=high_print_level,add_last=add_last_numeric,filename="",&
         unit_str="angstrom",error=error)

    CALL keyword_create(keyword, name="IONS+CENTERS",&
         description="prints out the wannier centers together with the particles",&
         usage="IONS+CENTERS", default_l_val=.FALSE., &
         lone_keyword_l_val=.TRUE., error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL add_format_keyword(keyword, print_key, pos=.TRUE.,& 
         description="Specifies the format of the output file when IONS+CENTERS is enabled.",&
         error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    NULLIFY(print_key)
    CALL cp_print_key_section_create(print_key,"WANNIER_SPREADS",&
         description="Controls the printing of the wannier functions", &
         print_level=high_print_level,add_last=add_last_numeric,filename="",&
         error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

    NULLIFY(print_key)
    CALL cp_print_key_section_create(print_key,"LOC_RESTART",&
         description="Controls the printing of restart file for localized MOS", &
         print_level=high_print_level,add_last=add_last_numeric,filename="",&
         error=error)
    CALL section_add_subsection(section,print_key,error=error)
    CALL section_release(print_key,error=error)

  END SUBROUTINE print_wanniers

! *****************************************************************************
!> \brief Create the print dft section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_print_dft_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_print_dft_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options in DFT code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(print_key, keyword, subsection)

       CALL cp_print_key_section_create(print_key,"PROGRAM_BANNER",&
            description="Controls the printing of the banner of the MM program",&
            print_level=silent_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"KINETIC_ENERGY",&
            description="Controls the printing of the kinetic energy",&
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DERIVATIVES",&
            description="Print all derivatives after the DFT calculation", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL keyword_create(keyword=keyword,&
            name="ndigits",&
            description="Specify the number of digits used to print derivatives",&
            default_i_val=6,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key_section=print_key,&
            name="neighbor_lists",&
            description="Controls the printing of the neighbor lists",&
            print_level=debug_print_level, filename="", unit_str="angstrom",&
            error=error)
       CALL keyword_create(keyword=keyword,&
            name="sab_orb",&
            description="Activates the printing of the orbital "//&
            "orbital neighbor lists, "//&
            "i.e. the overlap neighbor lists",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sab_aux_fit",&
            description="Activates the printing of the orbital "//&
            "orbital neighbor lists wavefunction fitting basis, "//&
            "i.e. the overlap neighbor lists",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sab_aux_fit_vs_orb",&
            description="Activates the printing of the orbital "//&
            "orbital mixed neighbor lists of wavefunction fitting basis, "//&
            "and the orbital basis, i.e. the overlap neighbor lists",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sab_scp",&
            description="Activates the printing of the vdW SCP "//&
            "neighbor lists ",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sab_vdw",&
            description="Activates the printing of the vdW "//&
            "neighbor lists (from DFT, DFTB, SE), "//&
            "i.e. the dispersion neighbor lists",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sab_cn",&
            description="Activates the printing of the "//&
            "neighbor lists used for coordination numbers in vdW DFT-D3",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sac_ae",&
            description="Activates the printing of the orbital "//&
            "nuclear attraction neighbor lists (erfc potential)",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sac_ppl",&
            description="Activates the printing of the orbital "//&
            "GTH-PPL neighbor lists (local part of the "//&
            "Goedecker-Teter-Hutter pseudo potentials)",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sap_ppnl",&
            description="Activates the printing of the orbital "//&
            "GTH-PPNL neighbor lists (non-local part of the"//&
            "Goedecker-Teter-Hutter pseudo potentials)",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sap_oce",&
            description="Activates the printing of the orbital "//&
            "PAW-projector neighbor lists (only GAPW)",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sab_se",&
            description="Activates the printing of the two-center "//&
            "neighbor lists for Coulomb type interactions in NDDO ",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sab_lrc",&
            description="Activates the printing of the long-range SE correction "//&
            "neighbor lists (only when doing long-range SE with integral scheme KDSO and KDSO-d)",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sab_tbe",&
            description="Activates the printing of the DFTB Ewald "//&
            "neighbor lists ",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sab_core",&
            description="Activates the printing of core interaction "//&
            "neighbor lists ",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="soo_list",&
            description="Activates the printing of RI orbital-orbital "//&
            "neighbor lists ",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sip_list",&
            description="Activates the printing of RI basis-projector interaction "//&
            "neighbor lists ",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SUBCELL",&
            description="Activates the printing of the subcells used for the"//&
            "generation of neighbor lists.", unit_str="angstrom",&
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"AO_MATRICES",&
            description="Controls the printing of the ao (i.e. contracted gaussian) matrices (debug).", &
            print_level=debug_print_level,filename="__STD_OUT__",&
            error=error)
       CALL keyword_create(keyword, name="CORE_HAMILTONIAN",&
            description="If the printkey is activated controls the printing of the hamiltonian matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="DENSITY",&
            description="If the printkey is activated controls the printing of the density (P) matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINETIC_ENERGY",&
            description="If the printkey is activated controls the printing of the kinetic energy matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="KOHN_SHAM_MATRIX",&
            description="If the printkey is activated controls the printing of the kohn-sham matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="MATRIX_VXC",&
            description="If the printkey is activated compute and print the matrix of the exchange and correlation potential."//&
                        "Only the GGA part for GPW is printed",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="ORTHO",&
            description="If the printkey is activated controls the printing of the orthogonalization matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="OVERLAP",&
            description="If the printkey is activated controls the printing of the overlap matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FERMI_CONTACT",&
            description="If the printkey is activated controls the printing of the Fermi contact matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="PSO",&
            description="If the printkey is activated controls the printing of the paramagnetic spin-orbit matrices",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="EFG",&
            description="If the printkey is activated controls the printing of the electric field gradient matrices",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="POTENTIAL_ENERGY",&
            description="If the printkey is activated controls the printing of the potential energy matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="OCE_HARD",&
            description="If the printkey is activated controls the printing of the OCE HARD matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="OCE_SOFT",&
            description="If the printkey is activated controls the printing of the OCE SOFT matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="W_MATRIX",&
            description="If the printkey is activated controls the printing of the w matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="W_MATRIX_AUX_FIT",&
            description="If the printkey is activated controls the printing of the w matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="DERIVATIVES",&
            description="If the printkey is activated controls the printing "//&
            "of derivatives (for the matrixes that support this)",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MO",&
            description="Controls the printing of the molecular orbitals."//&
            "Note that this is only functional with diagonalization based methods, in particular not with OT (see MO_CUBES)", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL keyword_create(keyword, name="Cartesian",&
            description="If the printkey is activated controls the printing of the mo in the cartesian basis",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="EIGENVALUES",variants=s2a("EIGVALS"),&
            description="If the printkey is activated controls the printing of the eigenvalues of the mos",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="EIGENVECTORS",variants=s2a("EIGVECS"),&
            description="If the printkey is activated controls the printing of the eigenvectors of the mos",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="OCCUPATION_NUMBERS",variants=s2a("OCCNUMS"),&
            description="If the printkey is activated controls the printing of the occupation numbers of the mos",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword=keyword, name="NDIGITS",&
            description="Specify the number of digits used to print the MO eigenvalues and occupation numbers",&
            default_i_val=6, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,&
                           name="MO_INDEX_RANGE",&
                           variants=s2a("MO_RANGE","RANGE"),&
                           description="Allows to print only a subset of the MO eigenvectors or eigenvalues. "//&
                                       "The indices of the first and the last MO have to be specified",&
                           required=.FALSE.,&
                           repeats=.FALSE.,&
                           n_var=2,&
                           type_of_var=integer_t,&
                           default_i_vals=(/0,0/),&
                           usage="MO_INDEX_RANGE 10 15",&
                           error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL create_mo_cubes_section(print_key,error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL create_stm_section(print_key,error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)


       CALL create_wfn_mix_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="GAPW",&
            description="Controls the printing of some gapw related information (debug).",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)
       CALL cp_print_key_section_create(print_key,"projectors",&
            description="If the printkey is activated controls if information on"//&
            " the projectors is printed.",&
            print_level=debug_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
       CALL cp_print_key_section_create(print_key,"rho0_information",&
            description="If the printkey is activated controls if information on rho0 is printed.",&
            print_level=debug_print_level,filename="__STD_OUT__",unit_str="angstrom",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL cp_print_key_section_create(print_key,"dft_control_parameters",&
            description="Controls the printing of dft control parameters.", &
            print_level=medium_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"OPTICAL_CONDUCTIVITY",&
            description="Controls the printing of the optical conductivity.", &
            print_level=high_print_level,filename="optical_conductivity",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,&
            name="E_DENSITY_CUBE",&
            description="Controls the printing of cube files with "//&
                        "the electronic density and, for LSD "//&
                        "calculations, the spin density",&
            print_level=high_print_level,filename="",error=error)
       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components.",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,&
            name="TOTAL_DENSITY",&
            description="Print the total electronic density in the case "//&
                        "of a GAPW run. This keyword has only an effect, "//&
                        "if PAW atoms are present. The default is to print "//&
                        "only the soft part of the electronic density and "//&
                        "to ignore the hard part. NOTE: The total density in "//&
                        "real space might exhibit unphysical features like "//&
                        "spikes due to the finite and thus truncated g vector "//&
                        "expansion",&
            usage="TOTAL_DENSITY {logical}",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="APPEND",&
            description="append the cube files when they already exist",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="XRD_INTERFACE",&
            description="It activates the print out of exponents and coefficients for the "//&
           " Gaussian expansion of the core densities, based on atom calculations for each kind."//&
           " The resulting core dansities are needed to compute the form factors."//&
           " If GAPW the local densities are also given in terms of a Gaussian expansion,"//& 
           " by fitting the difference between local-fhard and local-soft density for each atom."//&
           " In this case the keyword TOTAL_DENSITY is set to FALSE",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NGAUSS",&
            description="Number of Gaussian functions used in the expansion of atomic (core) density",&
            usage="NGAUSS 10",n_var=1,default_i_val=12, type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"tot_density_cube",&
            description="Controls printing of cube files with "//&
            "the total density (electrons+atomic core). Note that "//&
            "the value of the total density is positive where the "//&
            "electron density dominates and negative where the core is.",&
            print_level=high_print_level,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components.",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="APPEND",&
            description="append the cube files when they already exist",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"v_hartree_cube",&
            description="Controls the printing of a cube file with eletrostatic "//&
            " potential generated by the total density (electrons+ions). It is "//&
            " valid only for QS with GPW formalism .", &
            print_level=high_print_level,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components.",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="APPEND",&
            description="append the cube files when they already exist",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)


! ZMP adding the print section for  the v_xc cube 
       CALL cp_print_key_section_create(print_key,"v_xc_cube",&
            description="Controls the printing of a cube file with xc "//&
            " potential generated by the ZMP method (for the moment). It is "//&
            " valid only for QS with GPW formalism .", &
            print_level=high_print_level,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components.",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="APPEND",&
            description="append the cube files when they already exist",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"efield_cube",&
            description="Controls the printing of cube files with electric "//&
            " field generated by the total density (electrons+ions). It is "//&
            " valid only for QS with GPW formalism .", &
            print_level=high_print_level,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components.",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="APPEND",&
            description="append the cube files when they already exist",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"ELF_CUBE",&
            description="Controls printing of cube files with "//&
            "the electron localization function (ELF). Note that "//&
            "the value of ELF is defined between 0 and 1:  Pauli kinetic energy density normalized "//&
            " by the kinetic energy density of a uniform el. gas of same density.",&
            print_level=high_print_level,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components.",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="APPEND",&
            description="append the cube files when they already exist",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="density_cutoff",&
            description=" ",&
            usage="density_cutoff 0.0001",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            type_of_var=real_t,&
            default_r_val=1.0e-10_dp, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)


       CALL cp_print_key_section_create(print_key,"PDOS",&
            description="Print out the DOS projected per kind and per angular momentum  ",&
            print_level=high_print_level,common_iter_levels=1,filename="",&
            error=error)
       CALL keyword_create(keyword, name="COMPONENTS",&
            description="Print out pdos distinguishing all angular momentum components.",&
            usage="COMPONENTS", default_l_val=.FALSE., &
            lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="APPEND",&
            description="Append the pdos obtained at different iterations to the pdos  output file."//&
            "By defaut the file is overwritten",&
            usage="APPEND", default_l_val=.FALSE., &
            lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="NLUMO",&
            description="Number of virtual orbitals to be added to the MO set (-1=all)",&
            usage="NLUMO integer",default_i_val=0, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="OUT_EACH_MO",&
            description="Output on the status of the calculation every OUT_EACH_MO states. If -1 no output",&
            usage="OUT_EACH_MO integer",default_i_val=-1, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
!
       CALL section_create(subsection,name="LDOS",&
             description="Controls the printing of local PDOS, projected on subsets"//&
             " of atoms given through lists",&
             n_keywords=4, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
             error=error)
       CALL keyword_create(keyword, name="COMPONENTS",&
            description="Print out pdos distinguishing all angular momentum components.",&
            usage="COMPONENTS", default_l_val=.FALSE., &
            lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LIST",&
            description="Specifies a list of indexes of atoms where to project the DOS  ",&
            usage="LIST {integer}  {integer} ..  {integer} ",type_of_var=integer_t,&
            n_var=-1, required=.TRUE., repeats=.TRUE., error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(print_key,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="R_LDOS",&
             description="Controls the printing of local PDOS, projected on 3D volume in real space,"//&
             " the volume is defined in terms of position with respect to atoms in the lists",&
             n_keywords=4, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
             error=error)

       CALL keyword_create(keyword, name="LIST",&
            description="Specifies a list of indexes of atoms used to define the real space volume  ",&
            usage="LIST {integer}  {integer} ..  {integer} ",type_of_var=integer_t,&
            n_var=-1, required=.TRUE., repeats=.TRUE., error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="XRANGE",&
            description="range of positions in Cartesian direction x: all grid points within "//&
            " this range from at least one atom of the list are considered",&
            usage="XRANGE -10.0 10.0",unit_str="angstrom",n_var=2,type_of_var=real_t,  error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="YRANGE",&
            description="range of positions in Cartesian direction y: all grid points within "//&
            " this range from at least one atom of the list are considered",&
            usage="YRANGE -10.0 10.0",unit_str="angstrom",n_var=2,type_of_var=real_t,  error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="ZRANGE",&
            description="range of positions in Cartesian direction z: all grid points within "//&
            " this range from at least one atom of the list are considered",&
            usage="ZRANGE -10.0 10.0",unit_str="angstrom",n_var=2,type_of_var=real_t,  error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ERANGE",&
            description="only project states with the eigenvalues in the given interval. "//&
                        "Default is all states.",&
                        usage="ERANGE -1.0 1.0",unit_str="hartree",n_var=2,type_of_var=real_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL section_add_subsection(print_key,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

!
       !Printing of Moments
       CALL create_dipoles_section(print_key,"MOMENTS",high_print_level,error)
       CALL keyword_create(keyword=keyword,&
            name="MAX_MOMENT",&
            description="Maximum moment to be calculated",&
            usage="MAX_MOMENT {integer}",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            type_of_var=integer_t,&
            default_i_val=1,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword=keyword,&
            name="MAGNETIC",&
            description="Calculate also magnetic moments",&
            usage="MAGNETIC yes",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section, print_key, error=error)
       CALL section_release(print_key,error=error)

       ! Mulliken population analysis
       CALL cp_print_key_section_create(print_key,"MULLIKEN",&
            description="Controls the printing of the Mulliken (spin) population analysis", &
            print_level=medium_print_level,filename="__STD_OUT__",&
            common_iter_levels=1, error=error)
       CALL keyword_create(&
            keyword=keyword,&
            name="PRINT_GOP",&
            description="Print the gross orbital populations (GOP) in addition to the gross atomic populations (GAP) "//&
                        "and net charges",&
            usage="PRINT_GOP yes",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(&
            keyword=keyword,&
            name="PRINT_ALL",&
            description="Print all information including the full net AO and overlap population matrix",&
            usage="PRINT_ALL yes",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       ! Lowdin population analysis (fairly expensive to compute, so only at high)
       CALL cp_print_key_section_create(print_key,"LOWDIN",&
            description="Controls the printing of the Lowdin (spin) population analysis", &
            print_level=high_print_level,filename="__STD_OUT__",&
            common_iter_levels=1, error=error)
       CALL keyword_create(&
            keyword=keyword,&
            name="PRINT_GOP",&
            description="Print the orbital populations in addition to the atomic populations and net charges",&
            usage="PRINT_GOP yes",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(&
            keyword=keyword,&
            name="PRINT_ALL",&
            description="Print all information including the full symmetrically orthogonalised density matrix",&
            usage="PRINT_ALL yes",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(&
            print_key_section=print_key,&
            name="XRAY_DIFFRACTION_SPECTRUM",&
            description="Calculate and print the coherent X-ray "//&
            "diffraction spectrum",&
            print_level=debug_print_level,&
            filename="",&
            citations=(/Krack2000,Krack2002/),&
            error=error)
       CALL keyword_create(&
            keyword=keyword,&
            name="Q_MAX",&
            variants=(/"Q_MAXIMUM"/),&
            description="Maximum Q value calculated for the spectrum",&
            usage="Q_MAX {real}",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=20.0_dp,&
            unit_str="angstrom^-1",&
            error=error),&
            unit_str="angstrom^-1",&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key_section=print_key,&
            name="ELECTRIC_FIELD_GRADIENT",&
            description="Calculate and print the electric field gradients"//&
            "at atomic positions",&
            print_level=debug_print_level,&
            filename="__STD_OUT__",error=error)

       CALL keyword_create(keyword=keyword,&
            name="INTERPOLATION",&
            description="Use interpolation method from real space grid",&
            usage="INTERPOLATION {logical}",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="GSPACE_SMOOTHING",&
            description="Use a G-space smoothing function",&
            usage="GSPACE_SMOOTHING cutoff {real}, width {real}", &
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=2,default_r_vals=(/-1._dp,-1._dp/),&
            type_of_var=real_t, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="DEBUG",&
            description="Print additional debug output",&
            usage="DEBUG {logical}",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_gspace_interp_section(subsection,error=error)
       CALL section_add_subsection(print_key, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key_section=print_key,&
            name="BASIS_MOLOPT_QUANTITIES",&
            description="Print the two quantities needed in the basis molopt generation:"//&
            " total energy and condition number of the overlap matrix (S matrix)",&
            print_level=high_print_level,&
            filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key_section=print_key,&
            name="HYPERFINE_COUPLING_TENSOR",&
            description="Calculate and print the EPR hyperfine coupling tensor"//&
            " at atomic positions",&
            print_level=debug_print_level,&
            filename="__STD_OUT__",error=error)

       CALL keyword_create(keyword=keyword,&
            name="INTERACTION_RADIUS",&
            description="Radius of interaction for EPR hyperfine tensor calculation",&
            usage="INTERACTION_RADIUS radius {real}",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,default_r_val=10._dp,&
            type_of_var=real_t, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key_section=print_key,&
            name="OPTIMIZE_GEMINALS",&
            description="Optimize the parameters of the geminal basis set",&
            print_level=debug_print_level,&
            filename="__STD_OUT__",error=error)

       CALL create_powell_section(subsection,error)
       CALL section_add_subsection(print_key, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(&
            print_key_section=print_key,&
            name="PLUS_U",&
            description="Controls the printing for the DFT+U methods",&
            print_level=high_print_level,&
            filename="__STD_OUT__",&
            each_iter_names=s2a("QS_SCF"),&
            each_iter_values=(/0/),&
            citations=(/Dudarev1997,Dudarev1998/),&
            error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_print_dft_section

! *****************************************************************************
!> \brief creates the input section for dealing with homo lumos, including dumping cubes
! *****************************************************************************
  SUBROUTINE create_mo_cubes_section(print_key,error)
    TYPE(section_type), POINTER              :: print_key
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mo_cubes_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword

    NULLIFY(keyword)

    CALL cp_print_key_section_create(print_key,"MO_CUBES",&
         description="Controls the printing of cubes of the molecular orbitals.", &
         print_level=high_print_level,filename="",&
         error=error)
    CALL keyword_create(keyword, name="stride",&
         description="The stride (X,Y,Z) used to write the cube file "//&
         "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
         " 1 number valid for all components.",&
         usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="write_cube",&
         description="If the MO cube file should be written. If false, the eigenvalues are still computed."//&
                     " Can also be useful in combination with STM calculations",&
         default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="nlumo",&
         description="If the printkey is activated controls the number of lumos"//&
                     "that are printed and dumped as a cube (-1=all)",&
         default_i_val=0, error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="nhomo",&
         description="If the printkey is activated controls the number of homos that dumped as a cube (-1=all),"//&
                     " eigenvalues are always all dumped",&
         default_i_val=1, error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

!    CALL keyword_create(keyword, name="STM",&
!         description="do stm",&
!         default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
!    CALL section_add_keyword(print_key,keyword,error=error)
!    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="APPEND",&
         description="append the cube files when they already exist",&
         default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

!    CALL keyword_create(keyword, name="STM_BIAS",&
!         description="Bias energy for scanning tunneling microscopy (STM) image generation."//&
!                     "Orbital densities are summed according to the bias energy."//&
!                     "For negative values, states in the range ]EF+bias,EF] are summed,"//&
!                     "While positive values sum states in the range [EF,EF+bias[."//&
!                     "If postive biases are used, sufficiently many unoccupied stated"//&
!                     " (see ADDED_MOS and NLUMO ) should be computed.",&
!         n_var=-1,type_of_var=real_t, default_r_vals=(/0.0_dp/), unit_str='eV', error=error)
!    CALL section_add_keyword(print_key,keyword,error=error)
!    CALL keyword_release(keyword,error=error)

!    CALL keyword_create(keyword, name="TH_TORB",&
!         description="Tip orbital symmetry in Tersoff-Hamann approximation to compute STM images",&
!         required=.FALSE.,&
!         repeats=.TRUE.,&
!         default_i_val=orb_s,&
!         usage="TH_TORB s dz2",&
!         enum_c_vals=s2a("S","PX","PY","PZ","DXY","DYZ","DZX","DX2","DY2","DZ2"),&
!         enum_i_vals=(/orb_s,orb_px,orb_py,orb_pz,orb_dxy,orb_dyz,orb_dzx,orb_dx2,orb_dy2,orb_dz2/),&
!         enum_desc=s2a("s orbital","px orbital","py orbital","pz orbital",&
!         "dxy orbital","dyz orbital","dzx orbital","x^2 orbital","y^2 orbital","z^2 orbital"),&
!         error=error)
!    CALL section_add_keyword(print_key,keyword,error=error)
!    CALL keyword_release(keyword,error=error)

  END SUBROUTINE create_mo_cubes_section

  SUBROUTINE create_stm_section(print_key,error)
    TYPE(section_type), POINTER              :: print_key
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_stm_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword

    NULLIFY(keyword)

    CALL cp_print_key_section_create(print_key,"STM",&
         description="Controls the printing of cubes for the generation of STM images.", &
         print_level=high_print_level,filename="",&
         error=error)
    CALL keyword_create(keyword, name="stride",&
         description="The stride (X,Y,Z) used to write the cube file "//&
         "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
         " 1 number valid for all components.",&
         usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="nlumo",&
         description="If the printkey is activated controls the number of additional lumos"//&
                     " that are computed to be able to reproduce STM images obtained"//&
                     "  from positive bias (imaging unoccupied states)",&
         default_i_val=0, error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)


    CALL keyword_create(keyword, name="BIAS",&
         description="Bias energy for scanning tunneling microscopy (STM) image generation."//&
                     "Orbital densities are summed according to the bias energy."//&
                     "For negative values, states in the range ]EF+bias,EF] are summed,"//&
                     "While positive values sum states in the range [EF,EF+bias[."//&
                     "If postive biases are used, sufficiently many unoccupied stated"//&
                     " (see ADDED_MOS and NLUMO ) should be computed.",&
         n_var=-1,type_of_var=real_t, default_r_vals=(/0.0_dp/), unit_str='eV', error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="TH_TORB",&
         description="Tip orbital symmetry in Tersoff-Hamann approximation to compute STM images",&
         required=.FALSE.,&
         repeats=.TRUE.,&
         default_i_val=orb_s,&
         usage="TH_TORB s dz2",&
         enum_c_vals=s2a("S","PX","PY","PZ","DXY","DYZ","DZX","DX2","DY2","DZ2"),&
         enum_i_vals=(/orb_s,orb_px,orb_py,orb_pz,orb_dxy,orb_dyz,orb_dzx,orb_dx2,orb_dy2,orb_dz2/),&
         enum_desc=s2a("s orbital","px orbital","py orbital","pz orbital",&
         "dxy orbital","dyz orbital","dzx orbital","x^2 orbital","y^2 orbital","z^2 orbital"),&
         error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="REF_ENERGY",&
         description="By default the reference energy is the Fermi energy. In order to compare"//&
                     " with STS experiments, where specific energy ranges are addressed, here"//&
                     "  one can set a different reference energy."//&
                     " The energy range is anyway controlled by the BIAS",&
         type_of_var=real_t,default_r_val=0.0_dp, unit_str='eV',error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="APPEND",&
         description="append the cube files when they already exist",&
         default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
    CALL section_add_keyword(print_key,keyword,error=error)
    CALL keyword_release(keyword,error=error)



  END SUBROUTINE create_stm_section

  SUBROUTINE create_wfn_mix_section(section, error)

    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_wfn_mix_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    NULLIFY(subsection)
    NULLIFY(keyword)

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="WFN_MIX",&
            description="A section that allows manipulation of the MO coeffs,"//&
                        " e.g. for changing a ground state into an excited state."//&
                        "Starting from a copy of the original MOs, changes can be made"//&
                        "by adding linear combinations of HOMO/LUMO of the original MOs to the result MOs",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE., error=error)

       CALL keyword_create(keyword, name="OVERWRITE_MOS",&
            description="If the keyword is active molecular orbitals in memory will be replaced by the mixed wfn."//&
                        " In combination with RTP or EMD no restart will be required to use the mixed wfn.",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(subsection,name="UPDATE",&
            description="update a result MO with with a linear combination of of original MOs."//&
                        " This section can be repeated to build arbitrary linear combinations using repeatedly y=a*y+b*x.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE., error=error)

       CALL keyword_create(keyword, name="RESULT_MO_INDEX",&
            description="Index of the MO (y) to be modified. Counting down in energy with HOMO=1",&
            usage="RESULT_MO_INDEX 1", type_of_var=integer_t,default_i_val=0, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESULT_MARKED_STATE",&
            description="Specifies the MO according to "//&
                        "the marks set in MOLECULAR_STATES. The value corresponds to the repetition "//&
                        " of MARK_STATES in MOLECULAR_STATES",&
            usage="ORIG_MARKED_STATE 1", type_of_var=integer_t,default_i_val=0, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESULT_SPIN_INDEX",&
            description="Spin of the MO (y) to be modified.",&
            enum_c_vals=s2a("Alpha","Beta"),&
            enum_i_vals=(/ 1, 2/),&  ! direct index in array
            default_i_val=1,&
            enum_desc=s2a("Majority spin","Minority spin"), error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESULT_SCALE",&
            description="Scaling factor of the result variable (a).",&
            usage="RESULT_SCALE 0.0", type_of_var=real_t, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ORIG_MO_INDEX",&
            description="Index of the original MO (x). "//&
                        "Counting down in energy with HOMO=1 or up from LUMO=1, depending on ORIG_IS_VIRTUAL.",&
            usage="ORIG_MO_INDEX 1", type_of_var=integer_t,default_i_val=0, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ORIG_MARKED_STATE",&
            description="Specifies the MO according to "//&
                        "the marks set in MOLECULAR_STATES. The value corresponds to the repetition "//&
                        " of MARK_STATES in MOLECULAR_STATES",&
            usage="ORIG_MARKED_STATE 1", type_of_var=integer_t,default_i_val=0, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ORIG_SPIN_INDEX",&
            description="Spin of the MO (x) to be modified.",&
            enum_c_vals=s2a("Alpha","Beta"),&
            enum_i_vals=(/ 1, 2/),&  ! direct index in array
            default_i_val=1,&
            enum_desc=s2a("Majority spin","Minority spin"), error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ORIG_SCALE",&
            description="Scaling factor of the result variable (b).",&
            usage="ORIG_SCALE 0.0", type_of_var=real_t, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ORIG_IS_VIRTUAL",&
            description="The original MO (x) is a LUMO.",&
            usage="ORIG_IS_VIRTUAL", type_of_var=logical_t, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_wfn_mix_section

! *****************************************************************************
!> \brief creates the input section for the molecular states
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_molecular_states_section(print_key,error)
    TYPE(section_type), POINTER              :: print_key
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_molecular_states_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key2

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(print_key),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(print_key2, keyword)
       CALL cp_print_key_section_create(print_key,"MOLECULAR_STATES",&
            description="Controls printing of molecular states ",&
            print_level=high_print_level,filename=" ",citations=(/Hunt2003/), error=error)

       CALL keyword_create(keyword, name="CUBE_EVAL_RANGE",&
            description="only write cubes if the eigenvalues of the corresponding molecular states lie in the given interval. "//&
                        "Default is all states.",&
                        usage="CUBE_EVAL_RANGE -1.0 1.0",unit_str="hartree",n_var=2,type_of_var=real_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="MARK_STATES",&
            description="Can be used to mark given molecular states."//&
            " Sets a mark to both, occupied and unoccupied states. "//&
            "Occupied states are counted beginning with HOMO=1, "//&
            "unoccupied states are counted beginning with LUMO=1, "//&
            "This is only meaningful in combination with WFN_MIX. "//&
            "First integer specifies the molecule, second integer specifies the state.",&
            usage="MARK_STATES integer integer",&
            n_var=2,default_i_vals=(/-1,-1/), type_of_var=integer_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key2,"cubes",&
            description="Controls the printing of cube files", &
            print_level=high_print_level,filename="",error=error)
       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components.",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key2,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(print_key,print_key2,error=error)
       CALL section_release(print_key2,error=error)
    END IF
  END SUBROUTINE create_molecular_states_section


  SUBROUTINE create_wannier_states_section(print_key,error)
    TYPE(section_type), POINTER              :: print_key
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_wannier_states_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key2

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(print_key),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(print_key2, keyword)
       CALL cp_print_key_section_create(print_key,"WANNIER_STATES",&
            description="Controls printing of molecular states ",&
            print_level=high_print_level,filename=" ", error=error)

       CALL keyword_create(keyword, name="CUBE_EVAL_RANGE",&
            description="only write cubes if the eigenvalues of the corresponding molecular states lie in the given interval. "//&
                        "Default is all states.",&
                        usage="CUBE_EVAL_RANGE -1.0 1.0",unit_str="hartree",n_var=2,type_of_var=real_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="MARK_STATES",&
            description="Can be used to mark given molecular states."//&
            " Sets a mark to both, occupied and unoccupied states. "//&
            "Occupied states are counted beginning with HOMO=1, "//&
            "unoccupied states are counted beginning with LUMO=1, "//&
            "This is only meaningful in combination with WFN_MIX. "//&
            "First integer specifies the molecule, second integer specifies the state.",&
            usage="MARK_STATES integer integer",&
            n_var=2,default_i_vals=(/-1,-1/), type_of_var=integer_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key2,"cubes",&
            description="Controls the printing of cube files", &
            print_level=high_print_level,filename="",error=error)
       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components.",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key2,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(print_key,print_key2,error=error)
       CALL section_release(print_key2,error=error)
    END IF
  END SUBROUTINE create_wannier_states_section

! *****************************************************************************
!> \brief creates the input section for the qs part
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_qs_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qs_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"qs",&
            description="parameters needed to set up the Quickstep framework",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, subsection)

       ! Reals
       CALL keyword_create(keyword, name="EPS_DEFAULT",&
            description="Try setting all EPS_xxx to values leading to an energy correct up to EPS_DEFAULT",&
            usage="EPS_DEFAULT real", default_r_val=1.0E-10_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_CORE_CHARGE",&
            description="Precision for mapping the core charges.Overrides EPS_DEFAULT/100.0 value",&
            usage="EPS_CORE_CHARGE real", type_of_var=real_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_GVG_RSPACE",&
            variants=(/"EPS_GVG"/),&
            description="Sets precision of the realspace KS matrix element integration. Overrides SQRT(EPS_DEFAULT) value",&
            usage="EPS_GVG_RSPACE real",type_of_var=real_t ,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_PGF_ORB",&
            description="Sets precision of the overlap matrix elements. Overrides SQRT(EPS_DEFAULT) value",&
            usage="EPS_PGF_ORB real",type_of_var=real_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_KG_ORB",&
            description="Sets precision used in coloring the subsets for the Kim-Gordon method. Overrides SQRT(EPS_DEFAULT) value",&
            usage="EPS_KG_ORB 1.0E-8",&
            type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_3C_REDUCE",&
            description="GAPW: Factor to reduce the precision in the construction"//&
            " of the 3 center lists for the calculation of the OCE coefficients.",&
            usage="EPS_3C_REDUCE real", default_r_val=1.0_dp,&
            type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_PPL",&
            description="Adjusts the precision for the local part of the pseudo potential. ",&
            usage="EPS_PPL real", type_of_var=real_t, default_r_val=1.0E-2_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_PPNL",&
            description="Sets precision of the non-local part of the pseudo potential. Overrides sqrt(EPS_DEFAULT) value",&
            usage="EPS_PPNL real", type_of_var=real_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_CPC",&
            description="Sets precision of the GAPW projection. Overrides EPS_DEFAULT value",&
            usage="EPS_CPC real", type_of_var=real_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_RHO",&
            description="Sets precision of the density mapping on the grids.Overrides EPS_DEFAULT value",&
            usage="EPS_RHO real",type_of_var=real_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_RHO_RSPACE",&
            description="Sets precision of the density mapping in rspace.Overrides EPS_DEFAULT value."//&
            ".Overrides EPS_RHO value",&
            usage="EPS_RHO_RSPACE real",type_of_var=real_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_RHO_GSPACE",&
            description="Sets precision of the density mapping in gspace.Overrides EPS_DEFAULT value."//&
            ".Overrides EPS_RHO value",&
            usage="EPS_RHO_GSPACE real",type_of_var=real_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_FILTER_MATRIX",&
            description="Sets the threshold for filtering matrix elements.",&
            usage="EPS_FILTER_MATRIX 1.0E-6", type_of_var=real_t,default_r_val=0.0E0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPSFIT",&
            variants=(/"EPS_FIT"/),&
            description="GAPW and LRIPAW: precision to give the extention of a hard gaussian ",&
            usage="EPSFIT real", default_r_val=1.0E-4_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPSISO",&
            variants=(/"EPS_ISO"/),&
            description="GAPW : precision to determine an isolated projector",&
            usage="EPSISO real", default_r_val=1.0E-12_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPSSVD",&
            variants=(/"EPS_SVD"/),&
            description="GAPW : tolerance used in the singular value decomposition of the projector matrix",&
            usage="EPS_SVD real", default_r_val=1.0E-8_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPSRHO0",&
            variants=s2a("EPSVRHO0","EPS_VRHO0"),&
            description="GAPW : precision to determine the range of V(rho0-rho0soft)",&
            usage="EPSRHO0 real", default_r_val=1.0E-6_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALPHA0_HARD",&
            variants=s2a("ALPHA0_H","ALPHA0"),&
            description="GAPW and LRIPAW : Exponent for hard compensation charge",&
            usage="ALPHA0_HARD real", default_r_val=0.0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FORCE_PAW",&
            description="Use the GAPW scheme also for atoms with soft basis sets, i.e. "//&
            " the local densities are computed even if hard and soft should be equal. "//&
            "If this keyword is not set to true, those atoms with soft basis sets are treated by a GPW scheme, i.e. "//&
            "the corresponding density contribution goes on the global grid and is expanded in PW. "//&
            " This option nullifies the effect of the GPW_TYPE in the atomic KIND",&
            usage="FORCE_PAW",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_RAD_LOCAL",&
            description="GAPW : maximum radius of gaussian functions"//&
            " included in the generation of projectors",&
            usage="MAX_RAD_LOCAL real", default_r_val=25.0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Logicals
       CALL keyword_create(keyword, name="SCP",&
            description="Introduce additional self-consistent polarization through"//&
            " additional response basis functions (read in through AUX_BASIS.",&
            usage="SCP",lone_keyword_l_val=.TRUE.,&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LS_SCF",&
            description="Perform a linear scaling SCF",&
            usage="LS_SCF",lone_keyword_l_val=.TRUE.,&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALMO_SCF",&
            description="Perform ALMO SCF",&
            usage="ALMO_SCF",lone_keyword_l_val=.TRUE.,&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KG_METHOD",&
            description="Use a Kim-Gordon-like scheme.",&
            usage="KG_METHOD",lone_keyword_l_val=.TRUE.,&
            default_l_val=.FALSE.,citations=(/Iannuzzi2006, Brelaz1979/),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAP_CONSISTENT",&
            description="Compute the exact derivative (Hks) of the energy with respect to the density matrix. "//&
                        "This is slightly more expensive than using an approximate computation, "//&
                        "but consistent mapping can improve the stability of the SCF procedure, "//&
                        "especially for a tight EPS_SCF and a less tight EPS_DEFAULT.",&
            usage="MAP_CONSISTENT FALSE",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHECK_BCSR_CODE",&
            description="Check the BCSR code on actual data, once per QS run.",&
            usage="CHECK_BCSR_CODE TRUE",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BCSR_CODE",&
            description="Selects BCSR pathway.",&
            usage="BCSR_CODE 0",&
            default_i_val=0,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Integers
       CALL keyword_create(keyword, name="LMAXN1",&
            variants=(/"LMAXRHO1"/),&
            description="GAPW : max L number for espansion of the atomic densities in spherical gaussians",&
            usage="LMAXN1 integer",&
            default_i_val=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LMAXN0",&
            variants=(/"LMAXRHO0"/),&
            description="GAPW : max L number for the expansion compensation densities in spherical gaussians",&
            usage="LMAXN0 integer",&
            default_i_val=2,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LADDN0",&
            description="GAPW : integer added to the max L of the basis set, used to determine the "//&
            "maximum value of L for the compensation charge density.",&
            usage="LADDN0 integer",&
            default_i_val=99,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Characters
       CALL keyword_create(keyword, name="QUADRATURE",&
            description="GAPW and LRIPAW: algorithm to construct the atomic radial grids",&
            usage="QUADRATURE GC_SIMPLE",&
            enum_c_vals=s2a("GC_SIMPLE","GC_TRANSFORMED","GC_LOG"),&
            enum_i_vals=(/ do_gapw_gcs,do_gapw_gct,do_gapw_log/),&
            enum_desc=s2a("Gauss-Chebyshev",&
                          "Trans_Gauss-Chebyshev",&
                          "LogTrans-Gauss-Chebyshev"),&
            default_i_val=do_gapw_log, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PW_GRID",&
            description="What kind of PW_GRID should be employed",&
            usage="PW_GRID NS-FULLSPACE",&
            enum_c_vals=s2a("SPHERICAL","NS-FULLSPACE","NS-HALFSPACE"),&
            enum_desc=s2a("- not tested"," tested"," - not tested"),&
            enum_i_vals=(/ do_pwgrid_spherical, do_pwgrid_ns_fullspace,do_pwgrid_ns_halfspace/),&
            default_i_val=do_pwgrid_ns_fullspace, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PW_GRID_LAYOUT",&
            description="Force a particular real-space layout for the plane waves grids. "//&
                        "Numbers &le; 0 mean that this dimension is free, incorrect layouts will be ignored. "//&
                        "The default (/-1,-1/) causes CP2K to select a good value, "//&
                        "i.e. plane distributed for large grids, more general distribution for small grids.",&
            usage="PW_GRID_LAYOUT 4 16",&
            repeats=.FALSE.,required=.FALSE.,n_var=2,&
            default_i_vals=(/-1,-1/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PW_GRID_BLOCKED",&
            description="Can be used to set the distribution in g-space for the pw grids and their FFT.",&
            usage="PW_GRID_BLOCKED FREE",&
            enum_c_vals=s2a("FREE","TRUE","FALSE"),&
            enum_desc=s2a("CP2K will select an appropriate value","blocked","not blocked"),&
            enum_i_vals=(/do_pw_grid_blocked_free,do_pw_grid_blocked_true,do_pw_grid_blocked_false/),&
            default_i_val=do_pw_grid_blocked_free, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXTRAPOLATION",&
            variants=s2a("INTERPOLATION","WF_INTERPOLATION"),&
            description="Extrapolation strategy for the wavefunction during e.g. MD."//&
                        "PS and ASPC are recommended, see also EXTRAPOLATION_ORDER.",&
            citations=(/Kolafa2004,VandeVondele2005a/),&
            usage="EXTRAPOLATION PS",&
            enum_c_vals=s2a("USE_GUESS","USE_PREV_P","USE_PREV_RHO_R","LINEAR_WF",&
            "LINEAR_P","LINEAR_PS","USE_PREV_WF","PS","FROZEN","ASPC"),&
            enum_desc=s2a("Use the method specified with SCF_GUESS, i.e. no extrapolation",&
            "Use the previous density matrix",&
            "Use the previous density in real space",&
            "Linear extrapolation of the wavefunction",&
            "Linear extrapolation of the density matrix",&
            "Linear extrapolation of the density matrix times the overlap matrix",&
            "Use the previous wavefunction",&
            "Higher order extrapolation of the density matrix times the overlap matrix",&
            "Frozen ...",&
            "Always stable predictor corrector, similar to PS, but going for MD stability instead of intial guess accuracy."),&
            enum_i_vals=(/&
            wfi_use_guess_method_nr,&
            wfi_use_prev_p_method_nr,&
            wfi_use_prev_rho_r_method_nr,&
            wfi_linear_wf_method_nr,&
            wfi_linear_p_method_nr,&
            wfi_linear_ps_method_nr,&
            wfi_use_prev_wf_method_nr,&
            wfi_ps_method_nr,&
            wfi_frozen_method_nr,&
            wfi_aspc_nr/),&
            default_i_val=wfi_aspc_nr, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXTRAPOLATION_ORDER",&
            description="Order for the PS or ASPC extrapolation (typically 2-4). "//&
                        "Higher order might bring more accuracy, but comes, "//&
                        "for large systems, also at some cost. "//&
                        "In some cases, a high order extrapolation is not stable,"//&
                        " and the order needs to be reduced.",&
            usage="EXTRAPOLATION_ORDER {integer}",default_i_val=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="METHOD",&
            description="Specifies the electronic structure method that should be employed",&
            usage="METHOD GAPW",&
            enum_c_vals=s2a("GAPW","GAPW_XC","GPW","LRIPAW",&
            "MNDO","MNDOD","AM1","PM3","PM6","PDG","RM1","PNNL","DFTB","OFGPW","SCPTB"),&
            enum_desc=s2a("Gaussian and augmented plane waves method",&
                          "Gaussian and augmented plane waves method only for XC",&
                          "Gaussian and plane waves method",&
                          "Local resolution of identity projector augmented wave method",&
                          "MNDO semiempirical","MNDO-d semiempirical","AM1 semiempirical",&
                          "PM3 semiempirical","PM6 semiempirical","PDG semiempirical",&
                          "RM1 semiempirical",&
                          "PNNL semiempirical",&
                          "DFTB Density Functional based Tight-Binding",&
                          "OFGPW Orbital-free GPW method",&
                          "SCPTB Self-Consistent-Polarization Tight-Binding"),&
            enum_i_vals=(/ do_method_gapw, do_method_gapw_xc, do_method_gpw, do_method_lripaw,&
            do_method_mndo, do_method_mndod, do_method_am1, do_method_pm3, &
            do_method_pm6, do_method_pdg, do_method_rm1, do_method_pnnl, do_method_dftb, &
            do_method_ofgpw, do_method_scptb/),&
            citations=(/Lippert1997,Lippert1999,Krack2000,VandeVondele2005a,&
            VandeVondele2006,Dewar1977,Dewar1985,Rocha2006,Stewart1989,Thiel1992,&
            Repasky2002,Stewart2007,Schenter2008/),&
            default_i_val=do_method_gpw, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORE_PPL",&
            description="Specifies the method used to calculate the local pseudopotential contribution.",&
            usage="CORE_PPL ANALYTIC",&
            enum_c_vals=s2a("ANALYTIC","GRID"),&
            enum_desc=s2a("Analytic integration of integrals",&
                          "Numerical integration on real space grid. Lumped together with core charge"),&
            enum_i_vals=(/ do_ppl_analytic, do_ppl_grid/), &
            default_i_val=do_ppl_analytic, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_distribution_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_dftb_control_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_scptb_control_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_se_control_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_mulliken_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_ddapc_restraint_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_becke_restraint_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_s2_restraint_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_harris_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
    END IF

  END SUBROUTINE create_qs_section

! *****************************************************************************
  SUBROUTINE create_dftb_control_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_dftb_control_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"DFTB",&
            description="Parameters needed to set up the DFTB methods",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            citations=(/Porezag1995, Seifert1996, Elstner1998, Zhechkov2005/),&
            error=error)

       NULLIFY(subsection)
       CALL create_dftb_parameter_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="self_consistent",&
            description="Use self-consistent method",&
            citations=(/Elstner1998/),&
            usage="SELF_CONSISTENT",default_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="orthogonal_basis",&
            description="Assume orthogonal basis set",&
            usage="ORTHOGONAL_BASIS",default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="do_ewald",&
            description="Use Ewald type method instead of direct sum for Coulomb interaction",&
            usage="DO_EWALD",default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="dispersion",&
            description="Use dispersion correction",&
            citations=(/Zhechkov2005/),lone_keyword_l_val=.TRUE.,&
            usage="DISPERSION",default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="HB_SR_GAMMA",&
            description="Uses a modified version for the GAMMA within the SCC-DFTB scheme, "//&
            "specifically tuned for hydrogen bonds.",&
            citations=(/Hu2007/),lone_keyword_l_val=.TRUE.,&
            usage="HB_SR_GAMMA",default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_disp",&
            description="Define accuracy of dispersion interaction",&
            usage="EPS_DISP",default_r_val=0.0001_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_dftb_control_section

! *****************************************************************************
  SUBROUTINE create_scptb_control_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_scptb_control_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"SCPTB",&
            description="Parameters needed to set up the SCPTB methods",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword,name="PARAMETER_FILE_NAME",&
            description="Specify file that contains the atomic parameters",&
            usage="PARAMETER_FILE_NAME filename",&
            n_var=1,type_of_var=char_t,default_c_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="DISPERSION_PARAMETER_FILE",&
            description="Specify file that contains the atomic dispersion parameters",&
            usage="DISPERSION_PARAMETER_FILE filename",&
            n_var=1,type_of_var=char_t,default_c_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DISPERSION",&
            description="Use dispersion correction",&
            lone_keyword_l_val=.TRUE.,&
            usage="DISPERSION",default_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DISPERSION_RADIUS",&
            description="Define radius of dispersion interaction",&
            usage="DISPERSION_RADIUS",default_r_val=15._dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COORDINATION_CUTOFF",&
            description="Define cutoff for coordination number calculation",&
            usage="COORDINATION_CUTOFF",default_r_val=1.e-6_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="D3_SCALING",&
            description="Scaling parameters (s6,sr6,s8) for the D3 dispersion method,",&
            usage="D3_SCALING 1.0 1.0 1.0", n_var=3, default_r_vals=(/1.0_dp,1.0_dp,1.0_dp/),required=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STO_NG",&
            description="Provides the order of the Slater orbital expansion of Gaussian-Type Orbitals.",&
            usage="STO_NG",default_i_val=6, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PAIR_CUTOFF",&
            description="Define cutoff for pair potential calculation",&
            usage="PAIR_CUTOFF",default_r_val=1.e-8_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="do_ewald",&
            description="Use Ewald type method instead of direct sum for Coulomb interaction",&
            usage="DO_EWALD",default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="do_scc",&
            description="Use self consistent charge method. Can be used together with DO_SCP to get TB method",&
            usage="DO_SCC",default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="do_scp",&
            description="Use SCP method. Can be used to switch off SCP to get a SCC-DFTB method",&
            usage="DO_SCP",default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_scptb_control_section

! *****************************************************************************
  SUBROUTINE create_admm_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_admm_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"AUXILIARY_DENSITY_MATRIX_METHOD",&
            description="Parameters needed to set wavefunction fitting",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            citations=(/Guidon2010/),&
            error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="METHOD",&
            description="Method that shall be used for wavefunction fitting. Use BASIS_PROJECTION for MD.",&
            usage="METHOD BASIS_PROJECTION",&
            enum_c_vals=s2a("BASIS_PROJECTION","BLOCK_DENSITY_MATRIX"),&
            enum_i_vals=(/do_admm_basis_set_projection, do_admm_block_density_matrix/),&
            enum_desc=s2a("Construct auxiliary density matrix from auxiliary basis", &
                          "Construct auxiliary density from a blocked Fock matrix"), &
            citations=(/Guidon2010/),&
            default_i_val=do_admm_basis_set_projection, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="ADMM_PURIFICATION_METHOD",&
            description="Method that shall be used for wavefunction fitting. Use MO_DIAG for MD.",&
            usage="ADMM_PURIFICATION_METHOD NONE",&
            enum_c_vals=s2a("NONE","CAUCHY","CAUCHY_SUBSPACE","MO_DIAG", "MO_NO_DIAG"),&
            enum_i_vals=(/do_admm_purify_none, do_admm_purify_cauchy, do_admm_purify_cauchy_subspace,&
                          do_admm_purify_mo_diag, do_admm_purify_mo_no_diag/),&
            enum_desc=s2a("Do not apply any purification", &
                          "Apply McWeeny purification via general Cauchy representation", &
                          "Apply McWeeny purification via Cauchy representation in occupied subspace", &
                          "Calculate MO derivatives via Cauchy representation by diagonalization", &
                          "Calculate MO derivatives via Cauchy representation by inversion"), &
            citations=(/Guidon2010/),&
            default_i_val=do_admm_purify_mo_diag, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsection)
       CALL create_admm_block_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_admm_section

! *****************************************************************************
  SUBROUTINE create_admm_block_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_admm_block_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"BLOCK_DENSITY_MATRIX_METHOD",&
            description="Parameters needed to set wavefunction fitting",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="BASIS_PROJECTION",&
            description="Defines usage of auxiliary basis set",&
            usage="BASIS_PROJECTION ON",&
            enum_c_vals=s2a("ON","OFF"),&
            enum_i_vals=(/do_admm_block_aux_basis_on, do_admm_block_aux_basis_off/),&
            enum_desc=s2a("Use auxiliary basis", &
                          "Use primary basis"), &
            default_i_val=do_admm_block_aux_basis_on, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="PURIFICATION",&
            description="Method that shall be used for Purification",&
            usage="PURIFICATION ON",&
            enum_c_vals=s2a("OFF","FULL","BLOCKED"),&
            enum_i_vals=(/do_admm_block_purify_off, do_admm_block_purify_full, do_admm_block_purify_blocked/),&
            enum_desc=s2a("Do not apply any purification", &
                          "Apply full McWeeny purification via general Cauchy representation",&
                          "Apply blocked McWeeny purification via general Cauchy representation"), &
            default_i_val=do_admm_block_purify_full, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(keyword, subsection)
       ! FRAGMENT SECTION
       CALL section_create(subsection,name="BLOCK",&
            description="Specify the atom number belonging to this fragment.",&
            n_keywords=2, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="LIST",&
            description="Specifies a list of atoms.",&
            usage="LIST {integer} {integer} .. {integer}", required=.FALSE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_admm_block_section

! *****************************************************************************
  SUBROUTINE create_dftb_parameter_section(section, error)

    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_dftb_parameter_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,name="PARAMETER",&
            description="Information on where to find DFTB parameters",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword,name="SK_FILE",&
            description="Define parameter file for atom pair",&
            usage="SK_FILE a1 a2 filename",&
            n_var=3,type_of_var=char_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="PARAM_FILE_PATH",&
            description="Specify the directory with the DFTB parameter files. "//&
            "Used in combination with the filenames specified in the file "//&
            "given in PARAM_FILE_NAME.", usage="PARAM_FILE_PATH pathname",&
            n_var=1,type_of_var=char_t,default_c_val="./",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="PARAM_FILE_NAME",&
            description="Specify file that contains the names of "//&
            "Slater-Koster tables: A plain text file, each line has the "//&
            'format "ATOM1 ATOM2 filename.spl".',&
            usage="PARAM_FILE_NAME filename",&
            n_var=1,type_of_var=char_t,default_c_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="UFF_FORCE_FIELD",&
            description="Name of file with UFF parameters that will be used "//&
            "for the dispersion correction. Needs to be specified when "//&
            "DISPERSION==.TRUE., otherwise cp2k crashes with a Segmentation "//&
            "Fault.", usage="UFF_FORCE_FIELD filename",required=.FALSE.,&
            n_var=1,type_of_var=char_t,default_c_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_dftb_parameter_section

! *****************************************************************************
  SUBROUTINE create_se_control_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_se_control_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"SE",&
            description="Parameters needed to set up the Semi-empirical methods",&
            n_keywords=8, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="ORTHOGONAL_BASIS",&
            description="Assume orthogonal basis set. This flag is overwritten by "//&
                        "methods with fixed orthogonal/non-orthogonal basis set.",&
            usage="ORTHOGONAL_BASIS",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STO_NG",&
            description="Provides the order of the Slater orbital expansion of Gaussian-Type Orbitals.",&
            usage="STO_NG",default_i_val=6, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ANALYTICAL_GRADIENTS",&
            description="Nuclear Gradients are computed analytically or numerically",&
            usage="ANALYTICAL_GRADIENTS",default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELTA",&
            description="Step size in finite difference force calculation",&
            usage="DELTA {real} ",default_r_val=1.e-6_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INTEGRAL_SCREENING",&
            description="Specifies the functional form for the ",&
            usage="INTEGRAL_SCREENING (KDSO|KDSO-D|SLATER)",&
            enum_c_vals=s2a("KDSO","KDSO-D","SLATER"),&
            enum_i_vals=(/ do_se_IS_kdso, do_se_IS_kdso_d, do_se_IS_slater/),&
            enum_desc=s2a("Uses the standard NDDO Klopman-Dewar-Sabelli-Ohno equation "//&
            "for the screening of the Coulomb interactions.",&
            "Uses a modified Klopman-Dewar-Sabelli-Ohno equation, dumping the screening "//&
            "parameter for the Coulomb interactions.",&
            "Uses an exponential Slater-type function for modelling the Coulomb interactions."),&
            default_i_val=do_se_IS_kdso, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PERIODIC",&
            description="Specifies the type of treatment for the electrostatic long-range part "//&
            "in semi-empirical calculations.",&
            usage="PERIODIC (NONE|EWALD|EWALD_R3|EWALD_GKS)",&
            enum_c_vals=s2a("NONE","EWALD","EWALD_R3","EWALD_GKS"),&
            enum_i_vals=(/ do_se_lr_none, do_se_lr_ewald, do_se_lr_ewald_r3, do_se_lr_ewald_gks/),&
            enum_desc=s2a("The long-range part is not explicitly treaten. The interaction "//&
            "depends uniquely on the Cutoffs used for the calculation.",&
            "Enables the usage of Multipoles Ewald summation schemes. The short-range part "//&
            "is tapered according the value of RC_COULOMB.",&
            "Enables the usage of Multipoles Ewald summation schemes together with a long-range "//&
            "treatment for the 1/R^3 term, deriving from the short-range component. This option "//&
            "is active only for K-DSO integral screening type.",&
            "Use Ewald directly in Coulomb integral evaluation, works only with Slater screening"),&
            default_i_val=do_se_lr_none, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SCP",&
            description="Perform a SCP-NDDO calculation",&
            usage="SCP",default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FORCE_KDSO-D_EXCHANGE",&
            description="This keywords forces the usage of the KDSO-D integral screening "//&
            "for the Exchange integrals (default is to apply the screening only to the "//&
            "Coulomb integrals.",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsection)
       CALL create_coulomb_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_exchange_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_screening_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_lr_corr_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_neighbor_lists_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_se_memory_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_se_print_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_se_ga_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
    END IF

  END SUBROUTINE create_se_control_section

! *****************************************************************************
!> \brief Create the COULOMB se section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino]
!> \date  03.2009
! *****************************************************************************
  SUBROUTINE create_lr_corr_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_lr_corr_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="LR_CORRECTION",&
            description="Setup parameters for the evaluation of the long-range correction term in SE "//&
            "calculations.", n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="CUTOFF",&
            description="Atomic Cutoff Radius Cutoff for the evaluation of the long-ranbe correction integrals. ",&
            usage="CUTOFF {real} ",unit_str="angstrom",&
            default_r_val=cp_unit_to_cp2k(value=6.0_dp,unit_str="angstrom",error=error),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RC_TAPER",&
            description="Atomic Cutoff Radius Cutoff for Tapering the long-range correction integrals. "//&
            "If not specified it assumes the same value specified for the CUTOFF.",&
            usage="RC_TAPER {real} ",unit_str="angstrom",type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RC_RANGE",&
            description="Range of cutoff switch function (tapering): 0.5*(1-TANH((r-r0)/RC_RANGE)), "//&
            "where r0=2.0*RC_TAPER-20.0*RC_RANGE.",&
            usage="RC_RANGE {real} ",unit_str="angstrom",default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_lr_corr_section

! *****************************************************************************
!> \brief Create the COULOMB se section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino]
!> \date  03.2009
! *****************************************************************************
  SUBROUTINE create_coulomb_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_coulomb_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="COULOMB",&
            description="Setup parameters for the evaluation of the COULOMB term in SE "//&
            "calculations.", n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="CUTOFF",&
            description="Atomic Cutoff Radius Cutoff for the evaluation of the  Coulomb integrals. "//&
            "For non-periodic calculation the default value is exactly the full cell dimension, in order "//&
            "to evaluate all pair interactions. Instead, for periodic calculations the default numerical value is used." ,&
            usage="CUTOFF {real} ",unit_str="angstrom",&
            default_r_val=cp_unit_to_cp2k(value=12.0_dp,unit_str="angstrom",error=error),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RC_TAPER",&
            description="Atomic Cutoff Radius Cutoff for Tapering Coulomb integrals. "//&
            "If not specified it assumes the same value specified for the CUTOFF.",&
            usage="RC_TAPER {real} ",unit_str="angstrom",type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RC_RANGE",&
            description="Range of cutoff switch function (tapering): 0.5*(1-TANH((r-r0)/RC_RANGE)), "//&
            "where r0=2.0*RC_TAPER-20.0*RC_RANGE.",&
            usage="RC_RANGE {real} ",unit_str="angstrom",default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_coulomb_section

! *****************************************************************************
!> \brief Create the EXCHANGE se section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino]
!> \date  03.2009
! *****************************************************************************
  SUBROUTINE create_exchange_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_exchange_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="EXCHANGE",&
            description="Setup parameters for the evaluation of the EXCHANGE and "//&
            " core Hamiltonian terms in SE calculations.", n_keywords=0, n_subsections=1,&
            repeats=.FALSE., required=.FALSE.,error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="CUTOFF",&
            description="Atomic Cutoff Radius Cutoff for the evaluation of the Exchange integrals. "//&
            "For non-periodic calculation the default value is exactly the full cell dimension, in order "//&
            "to evaluate all pair interactions. Instead, for periodic calculations the default is the "//&
            "minimum value between 1/4 of the cell dimension and the value specified in input (either"//&
            " explicitly defined or the default numerical value).",&
            usage="CUTOFF {real} ",unit_str="angstrom",&
            default_r_val=cp_unit_to_cp2k(value=12.0_dp,unit_str="angstrom",error=error),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RC_TAPER",&
            description="Atomic Cutoff Radius Cutoff for Tapering Exchange integrals. "//&
            "If not specified it assumes the same value specified for the CUTOFF.",&
            usage="RC_TAPER {real} ",unit_str="angstrom",type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RC_RANGE",&
            description="Range of cutoff switch function (tapering): 0.5*(1-TANH((r-r0)/RC_RANGE)), "//&
            "where r0=2.0*RC_TAPER-20.0*RC_RANGE.",&
            usage="RC_RANGE {real} ",unit_str="angstrom",default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_exchange_section

! *****************************************************************************
!> \brief Create the SCREENING se section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino]
!> \date  03.2009
! *****************************************************************************
  SUBROUTINE create_screening_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_screening_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="SCREENING",&
            description="Setup parameters for the tapering of the Coulomb/Exchange Screening in "//&
            "KDSO-D integral scheme,", n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="RC_TAPER",&
            description="Atomic Cutoff Radius Cutoff for Tapering the screening term. ",&
            usage="RC_TAPER {real} ",unit_str="angstrom",&
            default_r_val=cp_unit_to_cp2k(value=12.0_dp,unit_str="angstrom",error=error),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RC_RANGE",&
            description="Range of cutoff switch function (tapering): 0.5*(1-TANH((r-r0)/RC_RANGE)), "//&
            "where r0=2*RC_TAPER-20*RC_RANGE.",&
            usage="RC_RANGE {real} ",unit_str="angstrom",default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_screening_section

! *****************************************************************************
!> \brief Create the print se section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_se_print_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_se_print_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options in SE code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"NEIGHBOR_LISTS",&
            description="Activates the printing of the neighbor lists used"//&
            " for the periodic SE calculations.", &
            print_level=high_print_level,filename="",unit_str="angstrom",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SUBCELL",&
            description="Activates the printing of the subcells used for the"//&
            "generation of neighbor lists for periodic SE.", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"EWALD_INFO",&
            description="Activates the printing of the information for "//&
            "Ewald multipole summation in periodic SE.", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_se_print_section

! *****************************************************************************
!> \brief creates the input section for use with the GA part of the code
!> \param section the section to create
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! *****************************************************************************
  SUBROUTINE create_se_ga_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_se_ga_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"GA",&
            description="Sets up memory parameters for the storage of the integrals",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(&
            keyword=keyword,&
            name="NCELLS",&
            description="Defines the number of linked cells for the neighbor list. "//&
                        "Default value is number of processors",&
            usage="NCELLS 10",&
            default_i_val=0,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_se_ga_section

! *****************************************************************************
!> \brief creates the input section for the se-memory part integral storage
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! *****************************************************************************
  SUBROUTINE create_se_memory_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_se_memory_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"MEMORY",&
            description="Sets up memory parameters for the storage of the integrals",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(&
            keyword=keyword,&
            name="EPS_STORAGE",&
            description="Storage threshold for compression is EPS_STORAGE",&
            usage="EPS_STORAGE 1.0E-10",&
            default_r_val=1.0E-10_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="MAX_MEMORY",&
            description="Defines the maximum amount of memory [MB] used to store precomputed "//&
                         "(possibly compressed) two-electron two-center integrals",&
            usage="MAX_MEMORY 256",&
            default_i_val=50,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COMPRESS",&
            description="Enables the compression of the integrals in memory.",&
            usage="COMPRESS <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF

  END SUBROUTINE create_se_memory_section


! *****************************************************************************
  SUBROUTINE create_mulliken_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mulliken_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"MULLIKEN_RESTRAINT",&
            description="Use mulliken charges in a restraint (check code for details)",&
            n_keywords=7, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="STRENGTH",&
            description="force constant of the restraint",&
            usage="STRENGTH {real} ",default_r_val=0.1_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TARGET",&
            description="target value of the restraint",&
            usage="TARGET {real} ",default_r_val=1._dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the list of atoms that is summed in the restraint",&
            usage="ATOMS {integer} {integer} .. {integer}", required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF

  END SUBROUTINE create_mulliken_section

! *****************************************************************************
  SUBROUTINE create_density_fitting_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_density_fitting_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    NULLIFY(keyword, print_key)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"DENSITY_FITTING",&
            description="Setup parameters for density fitting (Bloechl charges or density derived "//&
            " atomic point charges (DDAPC) charges)",&
            n_keywords=7, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            citations=(/Blochl1995/),&
            error=error)

       CALL keyword_create(keyword, name="NUM_GAUSS",&
            description="Specifies the numbers of gaussian used to fit the QM density for each atomic site.",&
            usage="NUM_GAUSS {integer}", required=.FALSE.,&
            n_var=1, type_of_var=integer_t, default_i_val=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PFACTOR",&
            description="Specifies the progression factor for the gaussian exponent for each atomic site.",&
            usage="PFACTOR {real}", required=.FALSE.,&
            n_var=1, type_of_var=real_t, default_r_val=1.5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MIN_RADIUS",&
            description="Specifies the smallest radius of the gaussian used in the fit. All other radius are"//&
            " obtained with the progression factor.",&
            usage="MIN_RADIUS {real}", required=.FALSE.,&
            unit_str="angstrom",n_var=1, type_of_var=real_t, default_r_val=0.5_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RADII",&
            description="Specifies all the radius of the gaussian used in the fit for each atomic site. The use"//&
            " of this keyword disables all other keywords of this section.",&
            usage="RADII {real} {real} .. {real}", required=.FALSE.,&
            unit_str="angstrom",n_var=-1, type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

        CALL keyword_create(keyword, name="GCUT",&
            description="Cutoff for charge fit in G-space.",&
            usage="GCUT {real}", required=.FALSE.,&
            n_var=1, type_of_var=real_t, default_r_val=SQRT(6.0_dp),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing of basic information during the run", &
            print_level=high_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)

       CALL keyword_create(keyword, name="CONDITION_NUMBER",&
            description="Prints information regarding the condition numbers of the A matrix (to be inverted)",&
            usage="ANALYTICAL_GTERM <LOGICAL>",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF

  END SUBROUTINE create_density_fitting_section

! *****************************************************************************
  SUBROUTINE create_ddapc_restraint_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_ddapc_restraint_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    NULLIFY(keyword, print_key)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"DDAPC_RESTRAINT",&
            description="Use DDAPC charges in a restraint (check code for details)",&
            n_keywords=7, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="TYPE_OF_DENSITY",&
            description="Specifies the type of density used for the fitting",&
            usage="TYPE_OF_DENSITY (FULL|SPIN)",&
            enum_c_vals=s2a("FULL","SPIN"),&
            enum_i_vals=(/ do_full_density, do_spin_density/),&
            enum_desc=s2a("Full density","Spin density"),&
            default_i_val=do_full_density, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STRENGTH",&
            description="force constant of the restraint",&
            usage="STRENGTH {real} ",default_r_val=0.1_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TARGET",&
            description="target value of the restraint",&
            usage="TARGET {real} ",default_r_val=1._dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the list of atoms that is summed in the restraint",&
            usage="ATOMS {integer} {integer} .. {integer}", required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="COEFF",&
            description="Defines the the coefficient of the atom in the atom list (default is one) ",&
            usage="COEFF 1.0 -1.0",&
            type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FUNCTIONAL_FORM",&
            description="Specifies the functional form of the term added",&
            usage="FUNCTIONAL_FORM RESTRAINT",&
            enum_c_vals=s2a("RESTRAINT","CONSTRAINT"),&
            enum_i_vals=(/ do_ddapc_restraint, do_ddapc_constraint/),&
            enum_desc=s2a("Harmonic potential: s*(q-t)**2","Constraint form: s*(q-t)"),&
            default_i_val=do_ddapc_restraint, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing basic info about the method", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF

  END SUBROUTINE create_ddapc_restraint_section

! *****************************************************************************
  SUBROUTINE create_becke_restraint_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_becke_restraint_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    NULLIFY(keyword, print_key)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"BECKE_RESTRAINT",&
            description="Use Becke weight population in a restraint/constraint ",&
            n_keywords=7, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="STRENGTH",&
            description="force constant of the restraint",&
            usage="STRENGTH {real} ",default_r_val=0.1_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TARGET",&
            description="target value of the restraint",&
            usage="TARGET {real} ",default_r_val=1._dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the list of atoms that is summed in the restraint",&
            usage="ATOMS {integer} {integer} .. {integer}", required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="COEFF",&
            description="Defines the the coefficient of the atom in the atom list (default is one)",&
            usage="COEFF 1.0 -1.0",&
            type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FUNCTIONAL_FORM",&
            description="Specifies the functional form of the term added",&
            usage="FUNCTIONAL_FORM RESTRAINT",&
            enum_c_vals=s2a("RESTRAINT","CONSTRAINT"),&
            enum_i_vals=(/ do_ddapc_restraint, do_ddapc_constraint/),&
            enum_desc=s2a("Harmonic potential: s*(q-t)**2","Constraint form: s*(q-t)"),&
            default_i_val=do_ddapc_restraint, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TYPE_OF_DENSITY",&
            description="Specifies the type of density used for the fitting",&
            usage="TYPE_OF_DENSITY (FULL|SPIN)",&
            enum_c_vals=s2a("FULL","SPIN"),&
            enum_i_vals=(/ do_full_density, do_spin_density/),&
            enum_desc=s2a("Full density","Spin density"),&
            default_i_val=do_full_density, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing basic info about the method", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF

  END SUBROUTINE create_becke_restraint_section

! *****************************************************************************
  SUBROUTINE create_s2_restraint_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_s2_restraint_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       CALL section_create(section,"S2_RESTRAINT",&
            description="Use S2 in a re/constraint (OT only)",&
            n_keywords=7, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="STRENGTH",&
            description="force constant of the restraint",&
            usage="STRENGTH {real} ",default_r_val=0.1_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TARGET",&
            description="target value of the restraint",&
            usage="TARGET {real} ",default_r_val=1._dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FUNCTIONAL_FORM",&
            description="Specifies the functional form of the term added",&
            usage="FUNCTIONAL_FORM RESTRAINT",&
            enum_c_vals=s2a("RESTRAINT","CONSTRAINT"),&
            enum_i_vals=(/ do_s2_restraint, do_s2_constraint/),&
            enum_desc=s2a("Harmonic potential: s*(q-t)**2","Constraint form: s*(q-t)"),&
            default_i_val=do_s2_restraint, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_s2_restraint_section

! *****************************************************************************
!> \brief creates the input section for the tddfpt part
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_tddfpt_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_tddfpt_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"tddfpt",&
            description="parameters needed to set up the Time Dependent Density Functional PT",&
            n_keywords=5, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            citations=(/Iannuzzi2005/),error=error)

       NULLIFY(subsection,keyword)

       ! Integer
       CALL keyword_create(keyword, name="MAX_KV",&
            variants=s2a("MAX_VECTORS"),&
            description=" maximal number of Krylov space vectors",&
            usage="MAX_KV someInteger>0",&
            n_var=1,type_of_var=integer_t,&
            default_i_val=60, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTARTS",&
            variants=s2a("N_RESTARTS"),&
            description=" maximal number subspace search restarts",&
            usage="RESTARTS someInteger>0",&
            n_var=1,type_of_var=integer_t,&
            default_i_val=5, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NEV",&
            variants=s2a("N_EV", "EV"),&
            description=" number of excitations to calculate",&
            usage="NEV someInteger>0",&
            n_var=1,type_of_var=integer_t,&
            default_i_val=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NLUMO",&
            description=" number of additional unoccupied orbitals ",&
            usage="NLUMO 10",&
            n_var=1,type_of_var=integer_t,&
            default_i_val=5, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NREORTHO",&
            variants=s2a("N_REORTHO","REORTHO","REORTHOGONALITAZIONS"),&
            description=" number of reorthogonalization steps",&
            usage="NREORTHO someInteger>0",&
            n_var=1,type_of_var=integer_t,&
            default_i_val=2, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Logical
       CALL keyword_create(keyword, name="KERNEL",&
            variants=s2a("DO_KERNEL"),&
            description="compute the kernel (debug purpose only)",&
            usage="KERNEL logical_value",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LSD_SINGLETS",&
            description="compute singlets using lsd vxc kernel",&
            usage="LSD_SINGLETS logical_value",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INVERT_S",&
            variants=s2a("INVERT_OVERLAP"),&
            description="use the inverse of the overlap matrix",&
            usage="INVERT_S logical_value",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PRECONDITIONER",&
            variants=s2a("PRECOND"),&
            description="use the preconditioner (only for Davidson)",&
            usage="PRECONDITIONER logical_value",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Character
       CALL keyword_create(keyword, name="RES_ETYPE",&
            variants=s2a("RESTRICTED_EXCITATIONS_TYPE", "RES_E_TYPE"),&
            description="(singlets/triplets) for restricted calculation",&
            usage="RES_ETYPE T",&
            enum_c_vals=s2a("S","SINGLET","SINGLETS","T","TRIPLET","TRIPLETS"),&
            enum_i_vals=(/ tddfpt_singlet, tddfpt_singlet, tddfpt_singlet,&
            tddfpt_triplet, tddfpt_triplet, tddfpt_triplet/),&
            default_i_val=tddfpt_singlet, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DIAG_METHOD",&
            variants=s2a("DIAGONALIZATION_METHOD", "METHOD"),&
            description="Diagonalization method used in tddfpt",&
            usage="DIAG_METHOD DAVIDSON",&
            enum_c_vals=s2a("DAVIDSON","LANCZOS"),&
            enum_i_vals=(/ tddfpt_davidson, tddfpt_lanczos/),&
            default_i_val=tddfpt_davidson, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="OE_CORR",&
            variants=s2a("ORBITAL_EIGENVALUES_CORRECTION"),&
            description="Which type of orbital eigenvalue correction to use\n"//&
            "(to yield better HOMO-LUMO energies)",&
            usage="OE_CORR SAOP",&
            enum_c_vals=s2a("NONE", "LB", "LB_ALPHA", "LB94", "GLLB", "GLB", "SAOP","SIC"),&
            enum_i_vals=(/ oe_none, oe_lb, oe_lb, oe_lb, oe_gllb, oe_gllb, oe_saop, oe_sic /),&
            default_i_val=oe_none, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Real
       CALL keyword_create(keyword, name="CONVERGENCE",&
            variants=s2a("CONV"),&
            description="The convergence of the eigenvalues",&
            usage="CONVERGENCE 1.0E-6 ",&
            n_var=1,type_of_var=real_t,&
            default_r_val=1.0e-5_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_xc_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_sic_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
    END IF

  END SUBROUTINE create_tddfpt_section

! *****************************************************************************
!> \brief creates the input section for the relativistic part
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author jens
! *****************************************************************************
  SUBROUTINE create_relativistic_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_relativistic_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"relativistic",&
            description="parameters needed and setup for relativistic calculations",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="method",&
            description="type of relativistic correction used",&
            usage="method DKH", default_i_val=rel_none,&
            enum_c_vals=s2a("NONE","DKH","ZORA"),&
            enum_i_vals=(/ rel_none, rel_dkh, rel_zora /),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DKH_order",&
            description="The order of the DKH transformation ",&
            usage="DKH_order 2", default_i_val=2,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ZORA_type",&
            description="Type of ZORA method to be used",&
            usage="ZORA_type scMP", default_i_val=rel_zora_full,&
            enum_c_vals=s2a("FULL","MP","scMP"), &
            enum_desc=s2a("Full ZORA method (not implemented)",&
                          "ZORA with atomic model potential",&
                          "Scaled ZORA with atomic model potential"),&
            enum_i_vals=(/ rel_zora_full, rel_zora_mp, rel_sczora_mp/),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="transformation",&
            description="type of DKH transformation, full: use full matrix transformation,"//&
             " molecule: use transformation blocked by molecule, atom: use atomic blocks ",&
            usage="transformation full", default_i_val=rel_trans_full,&
            enum_c_vals=s2a("FULL","MOLECULE","ATOM"), &
            enum_i_vals=(/ rel_trans_full, rel_trans_molecule, rel_trans_atom/),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="z_cutoff",&
            description="The minimal atomic number considered for atom transformation",&
            usage="z_cutoff 50", default_i_val=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="potential",&
            description="External potential used in DKH transformation, full 1/r or erfc(r)/r",&
            usage="POTENTIAL {FULL,ERFC}", default_i_val=rel_pot_erfc,&
            enum_c_vals=s2a("FULL","ERFC"),&
            enum_i_vals=(/ rel_pot_full, rel_pot_erfc /),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_relativistic_section

! *****************************************************************************
!> \brief creates the structure of the section with the dft scf parameters
!> \param section will contain the scf section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE create_scf_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_scf_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
    NULLIFY(print_key)

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"scf",&
            description="parameters needed perform an scf run",&
            n_keywords=24, n_subsections=3, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY (subsection)

       CALL create_ot_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_diagonalization_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_outer_scf_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_smear_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_mixing_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       NULLIFY (keyword)

       CALL keyword_create(keyword, name="MAX_ITER_LUMO",&
            variants=(/"MAX_ITER_LUMOS"/),&
            description="The maximum number of iteration for the lumo computation",&
            usage="MAX_ITER_LUMO 100", default_i_val=299,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_LUMO",&
            variants=(/"EPS_LUMOS"/),&
            description="target accuracy of the computation of the lumo energy",&
            usage="EPS_LUMO 1.e-6", default_r_val=1.0e-5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_SCF",&
            description="Maximum number of SCF iteration to be performed for one optimization",&
            usage="MAX_SCF 200", default_i_val=50,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_SCF_HISTORY", variants=(/"MAX_SCF_HIST"/), &
            description="Maximum number of SCF iterations after the history pipeline is filled", &
            usage="MAX_SCF_HISTORY 1", required=.FALSE., default_i_val=0,lone_keyword_i_val=1, &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="MAX_DIIS",&
            variants=(/"MAX_DIIS_BUFFER_SIZE"/),&
            description="Maximum number of DIIS vectors to be used",&
            usage="MAX_DIIS 3", default_i_val=4,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LEVEL_SHIFT",&
            variants=(/"LSHIFT"/),&
            description="Use level shifting to improve convergence",&
            usage="LEVEL_SHIFT 0.1", default_r_val=0._dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_SCF",&
            description="target accuracy for the scf convergence",&
            usage="EPS_SCF 1.e-6", default_r_val=1.e-5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_SCF_HISTORY", variants=(/"EPS_SCF_HIST"/), &
            description="target accuracy for the scf convergence after the history pipeline is filled",&
            usage="EPS_SCF_HISTORY 1.e-5", default_r_val=0.0_dp,lone_keyword_r_val=1.0e-5_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHOLESKY",&
            description="If the cholesky method should be used for computing "//&
            "the inverse of S, and in this case calling which Lapack routines",&
            usage="CHOLESKY REDUCE", default_i_val=cholesky_reduce,&
            enum_c_vals=s2a("OFF","REDUCE","RESTORE","INVERSE"),&
            enum_desc=s2a("The cholesky algorithm is not used","Reduce is called",&
            "Reduce is replaced by two restore",&
            "Restore uses operator multiply by inverse of the triangular matrix"),&
            enum_i_vals=(/cholesky_off,cholesky_reduce,cholesky_restore,cholesky_inverse/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_EIGVAL",&
            description="Throw away linear combinations of basis functions with a small eigenvalue in S",&
            usage="EPS_EIGVAL 1.0", default_r_val=1.0e-5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_DIIS",&
            description="Threshold on the convergence to start using DIAG/DIIS",&
            usage="EPS_DIIS 5.0e-2", default_r_val=0.1_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SCF_GUESS",&
            description="Change the initial guess for the wavefunction.",&
            usage="SCF_GUESS RESTART", default_i_val=atomic_guess,&
            enum_c_vals=s2a("ATOMIC","RESTART","RANDOM","CORE","DENSITIES",&
            "HISTORY_RESTART","MOPAC","SPARSE","ALMO_SCF","NONE"),&
            enum_desc=s2a("Generate an atomic density using the atomic code", &
            "Use the RESTART file as an initial guess (and ATOMIC if not present).", &
            "Use random wavefunction coefficients.", &
            "Diagonalize the core hamiltonian for an initial guess.", &
            "Use the aux_basis_set for collocation.", &
            "Extrapolated from previous RESTART files.", &
            "Use same guess as MOPAC for semi-empirical methods or a simple diagonal density matrix for other methods", &
            "Generate a sparse wavefunction using the atomic code (for OT based methods)", &
            "Accurate guess for molecular systems based on ALMO methods", &
            "Skip initial guess (only for NON-SCC DFTB)."), &
            enum_i_vals=(/atomic_guess,restart_guess,random_guess,core_guess,&
            densities_guess,history_guess,mopac_guess,sparse_guess,almo_scf_guess,no_guess/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NROW_BLOCK",&
            description="sets the number of rows in a scalapack block",&
            usage="NROW_BLOCK 31", default_i_val=32,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NCOL_BLOCK",&
            description="Sets the number of columns in a scalapack block",&
            usage="NCOL_BLOCK 31", default_i_val=32,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ADDED_MOS",&
            description="Number of additional MOS added for each spin",&
            usage="ADDED_MOS", default_i_val=0,n_var=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="ROKS_SCHEME",&
            description="Selects the ROKS scheme when ROKS is applied.",&
            usage="ROKS_SCHEME HIGH-SPIN",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            enum_c_vals=s2a("GENERAL","HIGH-SPIN"),&
            enum_i_vals=(/general_roks,high_spin_roks/),&
            default_i_val=high_spin_roks,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="ROKS_F",&
            variants=(/"F_ROKS"/),&
            description="Allows to define the parameter f for the "//&
            "general ROKS scheme.",&
            usage="ROKS_PARAMETER 1/2",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            type_of_var=real_t,&
            default_r_val=0.5_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="ROKS_PARAMETERS",&
            variants=(/"ROKS_PARAMETER"/),&
            description="Allows to define all parameters for the high-spin "//&
            "ROKS scheme explicitly. "//&
            "The full set of 6 parameters has to be specified "//&
            "in the order acc, bcc, aoo, boo, avv, bvv",&
            usage="ROKS_PARAMETERS 1/2 1/2 1/2 1/2 1/2 1/2",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=6,&
            type_of_var=real_t,&
            default_r_vals=(/-0.5_dp,1.5_dp,0.5_dp,0.5_dp,1.5_dp,-0.5_dp/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(subsection,"PRINT","printing of information during the scf",&
            error=error, repeats=.FALSE., required=.FALSE.)

        CALL cp_print_key_section_create(print_key,"RESTART",&
            description="Controls the dumping of the MO restart file during scf."//&
                        "By default keeps a short history of three restarts."//&
                        "See also RESTART_HISTORY", &
            print_level=low_print_level, common_iter_levels=3,&
            each_iter_names=s2a("QS_SCF"),each_iter_values=(/20/), &
            add_last=add_last_numeric,filename="RESTART",error=error)
       CALL keyword_create(keyword, name="BACKUP_COPIES",&
             description="Specifies the maximum index of backup copies.",&
             usage="BACKUP_COPIES {int}",&
             default_i_val=3, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"RESTART_HISTORY",&
            description="Dumps unique MO restart files during the run keeping all of them.",&
            print_level=low_print_level, common_iter_levels=0,&
            each_iter_names=s2a("__ROOT__","MD","GEO_OPT","ROT_OPT","NEB","METADYNAMICS","QS_SCF"),&
            each_iter_values=(/500,500,500,500,500,500,500/), &
            filename="RESTART",error=error)
       CALL keyword_create(keyword, name="BACKUP_COPIES",&
             description="Specifies the maximum index of backup copies.",&
             usage="BACKUP_COPIES {int}",&
             default_i_val=3, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"iteration_info",&
            description="Controls the printing of basic iteration information during the scf.", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL keyword_create(keyword, name="time_cumul",&
            description="If the printkey is activated switches the printing of timings"//&
            " to cumulative (over the scf).",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing of basic information during the SCF", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MO_ORTHONORMALITY",&
            description="Controls the printing relative to the orthonormality of MOs (CT S C).", &
            print_level=high_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MO_MAGNITUDE",&
            description="Prints the min/max eigenvalues of the overlap of the MOs without S (CT C).", &
            print_level=high_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"detailed_energy",&
            description="Controls the printing of detailed energy information.", &
            print_level=high_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"diis_info",&
            description="Controls the printing of diis information.", &
            print_level=high_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"total_densities",&
            description="Controls the printing of total densities.", &
            print_level=medium_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"Lanczos",&
            description="Controls the printing of information on Lanczos refinement iterations.", &
            print_level=medium_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"DIAG_SUB_SCF",&
            description="Controls the printing of information on subspace diagonalization internal loop. ", &
            print_level=medium_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"Davidson",&
            description="Controls the printing of information on Davidson iterations.", &
            print_level=medium_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_scf_section

! *****************************************************************************
!> \brief creates the KG section
!> \author Martin Haeufel [2012.07]
! *****************************************************************************
  SUBROUTINE create_kg_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_kg_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"KG_METHOD",&
            description="Specifies the parameters for a Kim-Gordon-like partitioning"//&
                        " into molecular subunits",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            citations=(/Iannuzzi2006, Brelaz1979/),error=error)
      
       NULLIFY (keyword, subsection, print_key)

       CALL keyword_create(keyword, name="COLORING_METHOD",&
            description="Which algorithm to use for coloring.",&
            usage="COLORING_METHOD GREEDY",&
            default_i_val=kg_color_dsatur,&
            enum_c_vals=s2a( "DSATUR", "GREEDY"),&
            enum_desc=s2a("Maximum degree of saturation, relatively accurate",&
                          "Greedy, fast coloring, less accurate"),&
            enum_i_vals=(/kg_color_dsatur, kg_color_greedy /),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(subsection, name="PRINT",&
            description="Print section",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL cp_print_key_section_create(print_key,"NEIGHBOR_LISTS",&
            description="Controls the printing of the neighbor lists.", &
            print_level=low_print_level, filename="__STD_OUT__", unit_str="angstrom",&
            error=error)

       CALL keyword_create(keyword=keyword,&
            name="SAB_ORB_FULL",&
            description="Activates the printing of the full orbital "//&
            "orbital neighbor lists.", &
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="SAB_ORB_MOLECULAR",&
            description="Activates the printing of the orbital "//&
            "orbital neighbor lists for molecular subsets.",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection, print_key, error=error)
       CALL section_release(print_key, error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection, error=error)

    END IF

  END SUBROUTINE create_kg_section

! *****************************************************************************
!> \brief creates the structure of the section with scf parameters
!>      controlling an other loop
!> \param section will contain the scf section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Joost VandeVondele [2006.03]
! *****************************************************************************
  SUBROUTINE create_outer_scf_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_outer_scf_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"OUTER_SCF",&
            description="parameters controlling the outer SCF loop",&
            n_keywords=9, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the outer SCF loop",&
            usage="&OUTER_SCF ON",default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TYPE",&
            description="Specifies which kind of outer SCF should be employed",&
            usage="TYPE DDAPC_CONSTRAINT ",&
            default_i_val=outer_scf_none,&
            enum_c_vals=s2a( "DDAPC_CONSTRAINT", "S2_CONSTRAINT", "BECKE_CONSTRAINT",&
             "SCP", "NONE" ),&
            enum_desc=s2a("Enforce a constraint on the DDAPC, requires the corresponding section", &
            "Enforce a constraint on the S2, requires the corresponding section", &
            "Enforce a constraint on the Becke weight population,requires the corresponding section", &
            "Enforce outer loop optimization for SCP, requires the corresponding metod (-DFT or -NDDO)", &
            "Do nothing in the outer loop, useful for resetting the inner loop,"),&
            enum_i_vals=(/outer_scf_ddapc_constraint,outer_scf_s2_constraint,&
            outer_scf_becke_constraint,outer_scf_scp,outer_scf_none/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="OPTIMIZER",&
            description="Method used to bring the outer loop to a stationary point",&
            usage="OPTIMIZER SD",&
            default_i_val=outer_scf_optimizer_none,&
            enum_c_vals=s2a("SD","DIIS","NONE","BISECT"),&
            enum_desc=s2a("Takes steps in the direction of the gradient, multiplied by step_size", &
            "Uses a Direct Inversion in the Iterative Subspace method", &
            "Do nothing, useful only with the none type",&
            "Bisection on the gradient, useful for difficult one dimensional cases"),&
            enum_i_vals=(/outer_scf_optimizer_sd,outer_scf_optimizer_diis,outer_scf_optimizer_none,outer_scf_optimizer_bisect/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BISECT_TRUST_COUNT",&
            description="Maximum number of times the same point will be used in bisection,"//&
                        " a small number guards against the effect of wrongly converged states.", &
            usage="BISECT_TRUST_COUNT 5", default_i_val=10,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_SCF",&
            description="The target gradient of the outer scf variables. "//&
            "Notice that the EPS_SCF of the inner loop also determines "//&
            "the value that can be reached in the outer loop, "//&
            "typically EPS_SCF of the outer loop must be smaller "//&
            "than EPS_SCF of the inner loop.", &
            usage="EPS_SCF 1.0E-6 ", default_r_val=1.0E-5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DIIS_BUFFER_LENGTH",&
            description="Maximum number of DIIS vectors used ", &
            usage="DIIS_BUFFER_LENGTH 5", default_i_val=3,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXTRAPOLATION_ORDER",&
            description="Number of past states used in the extrapolation of the variables during e.g. MD", &
            usage="EXTRAPOLATION_ORDER 5", default_i_val=3,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_SCF",&
            description="The maximum number of outer loops ", &
            usage="MAX_SCF 20", default_i_val=50,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STEP_SIZE",&
            description="The initial step_size used in the optimizer (currently steepest descent)."//&
            "Note that in cases where a sadle point is sought for (DDAPC_CONSTRAINT),"//&
            " this can be negative", &
            usage="STEP_SIZE -1.0", default_r_val=0.5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_outer_scf_section


! *****************************************************************************
!> \brief Create the BSSE section for counterpoise correction
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_bsse_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_bsse_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BSSE",&
            description="This section is used to set up the BSSE calculation. "//&
            "It also requires that for each atomic kind X a kind X_ghost is present, "//&
            "with the GHOST keyword specified, in addition to the other required fields.",&
            n_keywords=3, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, subsection)
       ! FRAGMENT SECTION
       CALL section_create(subsection,name="FRAGMENT",&
            description="Specify the atom number belonging to this fragment.",&
            n_keywords=2, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="LIST",&
            description="Specifies a list of atoms.",&
            usage="LIST {integer} {integer} .. {integer}", required=.FALSE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! CONFIGURATION SECTION
       CALL section_create(subsection,name="CONFIGURATION",&
            description="Specify additional parameters for the combinatorial configurations.",&
            n_keywords=2, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="GLB_CONF",&
            description="Specifies the global configuration using 1 or 0.",&
            usage="GLB_CONF {integer} {integer} .. {integer}", required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SUB_CONF",&
            description="Specifies the subconfiguration using 1 or 0 belonging to the global configuration.",&
            usage="SUB_CONF {integer} {integer} .. {integer}", required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="MULTIPLICITY",&
            variants=(/"MULTIP"/),&
            description="Specify for each fragment the multiplicity. Two times the total spin plus one. "//&
            "Specify 3 for a triplet, 4 for a quartet,and so on. Default is 1 (singlet) for an "//&
            "even number and 2 (doublet) for an odd number of electrons.",&
            usage="MULTIPLICITY 3",&
            default_i_val=0,& ! this default value is just a flag to get the above
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE",&
            description="The total charge for each fragment.",&
            usage="CHARGE -1",&
            default_i_val=0,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="FRAGMENT_ENERGIES",&
            description="This section contains the energies of the fragments already"//&
            " computed. It is useful as a summary and specifically for restarting BSSE runs.",&
            n_keywords=2, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="The energy computed for each fragment",repeats=.TRUE.,&
            usage="{REAL}", type_of_var=real_t, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_print_bsse_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_bsse_section

! *****************************************************************************
!> \brief Create the print bsse section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_print_bsse_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_print_bsse_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options in BSSE code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of information regarding the run.",&
            print_level=low_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"RESTART",&
            description="Controls the dumping of the restart file during BSSE runs."//&
            "By default the restart is updated after each configuration calculation is "//&
            " completed.", &
            print_level=silent_print_level, common_iter_levels=0,  &
            add_last=add_last_numeric,filename="",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_print_bsse_section

! *****************************************************************************
!> \brief creates the interpolation section for the periodic QM/MM
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author tlaino
! *****************************************************************************
  SUBROUTINE create_gspace_interp_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_gspace_interp_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="interpolator",&
            description="controls the interpolation for the G-space term",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, print_key)

       CALL keyword_create(keyword, name="aint_precond",&
            description="the approximate inverse to use to get the starting point"//&
            " for the linear solver of the spline3 methods",&
            usage="kind spline3",&
            default_i_val=precond_spl3_aint,&
            enum_c_vals=s2a( "copy","spl3_nopbc_aint1","spl3_nopbc_precond1",&
            "spl3_nopbc_aint2","spl3_nopbc_precond2","spl3_nopbc_precond3"),&
            enum_i_vals=(/no_precond,precond_spl3_aint, precond_spl3_1,&
            precond_spl3_aint2, precond_spl3_2, precond_spl3_3/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="precond",&
            description="The preconditioner used"//&
            " for the linear solver of the spline3 methods",&
            usage="kind spline3",&
            default_i_val=precond_spl3_3,&
            enum_c_vals=s2a("copy","spl3_nopbc_aint1","spl3_nopbc_precond1",&
            "spl3_nopbc_aint2","spl3_nopbc_precond2","spl3_nopbc_precond3"),&
            enum_i_vals=(/no_precond,precond_spl3_aint, precond_spl3_1,&
            precond_spl3_aint2, precond_spl3_2, precond_spl3_3/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_x",&
            description="accuracy on the solution for spline3 the interpolators",&
            usage="eps_x 1.e-15", default_r_val=1.e-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_r",&
            description="accuracy on the residual for spline3 the interpolators",&
            usage="eps_r 1.e-15", default_r_val=1.e-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="max_iter",&
            variants=(/'maxiter'/),&
            description="the maximum number of iterations",&
            usage="max_iter 200", default_i_val=100, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"conv_info",&
            description="if convergence information about the linear solver"//&
            " of the spline methods should be printed", &
            print_level=medium_print_level,each_iter_names=s2a("SPLINE_FIND_COEFFS"),&
            each_iter_values=(/10/),filename="__STD_OUT__",&
            add_last=add_last_numeric,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_gspace_interp_section

! *****************************************************************************
!> \brief This section specifies the flags for the harris functional
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE create_harris_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_harris_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="HARRIS",&
            description="This section specifies the flags for the calculation of the harris functional",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ACTIVATE", &
            description="Activates the harris functional.", &
            usage="ACTIVATE TRUE", required=.FALSE., &
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="HARRIS_ENERGY",&
            description="Overrides the QS Energy.", &
            usage="HARRIS_ENERGY TRUE", required=.FALSE., &
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

    END IF
  END SUBROUTINE create_harris_section

! *****************************************************************************
!> \brief creates the multigrid
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE create_mgrid_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_mgrid_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="mgrid",&
            description="multigrid information",&
            n_keywords=5, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="NGRIDS",&
            description="The number of multigrids to use",&
            usage="ngrids 1", default_i_val=4, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword, name="cutoff",&
            description="The cutoff of the finest grid level. Default value for "//&
            "SE or DFTB calculation is 1.0 [Ry].",&
            usage="cutoff 300",default_r_val=cp_unit_to_cp2k(value=280.0_dp,&
            unit_str="Ry",error=error), n_var=1, unit_str="Ry", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="progression_factor",&
            description="Factor used to find the cutoff of the multigrids that"//&
            " where not given explicitly",&
            usage="progression_factor <integer>", default_r_val=3._dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="commensurate",&
            description="If the grids should be commensurate. If true overrides "//&
            "the progression factor and the cutoffs of the sub grids",&
            usage="commensurate", default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="realspace",&
            description="If both rho and rho_gspace are needed ",&
            usage="realspace", default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="REL_CUTOFF",&
            variants=(/"RELATIVE_CUTOFF"/),&
            description="Determines the grid at which a Gaussian is mapped,"//&
            " giving the cutoff used for a gaussian with alpha=1."//&
            " A value 50+-10Ry might be required for highly accurate results, "//&
            " Or for simulations with a variable cell."//&
            " Versions prior to 2.3 used a default of 30Ry.",&
            usage="RELATIVE_CUTOFF real", default_r_val=20.0_dp,&
            unit_str="Ry",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MULTIGRID_SET",&
            description="Activate a manual setting of the multigrids",&
            usage="MULTIGRID_SET", default_l_val=.FALSE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SKIP_LOAD_BALANCE_DISTRIBUTED",&
            description="Skip load balancing on distributed multigrids, which might be memory intensive."//&
                        "If not explicitly specified, runs using more than 1024 MPI tasks will default to .TRUE.",&
            usage="SKIP_LOAD_BALANCE_DISTRIBUTED", default_l_val=.FALSE., lone_keyword_l_val=.TRUE., &
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="MULTIGRID_CUTOFF",&
            variants=(/"CUTOFF_LIST"/),&
            description="List of cutoff values to set up multigrids manually",&
            usage="MULTIGRID_CUTOFF 200.0 100.0 ", n_var=-1, type_of_var=real_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsection)
       CALL create_rsgrid_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_interp_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
    END IF
  END SUBROUTINE create_mgrid_section

! *****************************************************************************
!> \brief creates the interpolation section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author tlaino
! *****************************************************************************
  SUBROUTINE create_interp_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_interp_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="interpolator",&
            description="kind of interpolation used between the multigrids",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, print_key)

       CALL keyword_create(keyword, name="kind",&
            description="the interpolator to use",&
            usage="kind spline3",&
            default_i_val=pw_interp,&
            enum_c_vals=s2a("pw","spline3_nopbc","spline3"),&
            enum_i_vals=(/pw_interp,&
            spline3_nopbc_interp,spline3_pbc_interp/),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="safe_computation",&
            description="if a non unrolled calculation is to be performed in parallel",&
            usage="safe_computation OFF",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="aint_precond",&
            description="the approximate inverse to use to get the starting point"//&
            " for the linear solver of the spline3 methods",&
            usage="aint_precond copy",&
            default_i_val=precond_spl3_aint,&
            enum_c_vals=s2a( "copy","spl3_nopbc_aint1","spl3_nopbc_aint2",&
            "spl3_nopbc_precond1","spl3_nopbc_precond2","spl3_nopbc_precond3"),&
            enum_i_vals=(/no_precond,precond_spl3_aint, precond_spl3_aint2,&
            precond_spl3_1,precond_spl3_2,precond_spl3_3/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="precond",&
            description="The preconditioner used"//&
            " for the linear solver of the spline3 methods",&
            usage="PRECOND copy",&
            default_i_val=precond_spl3_3,&
            enum_c_vals=s2a( "copy","spl3_nopbc_aint1","spl3_nopbc_aint2",&
            "spl3_nopbc_precond1","spl3_nopbc_precond2","spl3_nopbc_precond3"),&
            enum_i_vals=(/no_precond,precond_spl3_aint, precond_spl3_aint2,&
            precond_spl3_1,precond_spl3_2,precond_spl3_3/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_x",&
            description="accuracy on the solution for spline3 the interpolators",&
            usage="eps_x 1.e-15", default_r_val=1.e-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_r",&
            description="accuracy on the residual for spline3 the interpolators",&
            usage="eps_r 1.e-15", default_r_val=1.e-10_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="max_iter",&
            variants=(/'maxiter'/),&
            description="the maximum number of iterations",&
            usage="max_iter 200", default_i_val=100, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"conv_info",&
            description="if convergence information about the linear solver"//&
            " of the spline methods should be printed", &
            print_level=medium_print_level,each_iter_names=s2a("SPLINE_FIND_COEFFS"),&
            each_iter_values=(/10/),filename="__STD_OUT__",&
            add_last=add_last_numeric,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_interp_section

! *****************************************************************************
!> \brief creates the sic (self interaction correction) section
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE create_sic_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_sic_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"sic",&
            description="parameters for the self interaction correction",&
            n_keywords=6, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            citations=(/VandeVondele2005b,Perdew1981,Avezac2005/),&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="SIC_SCALING_A",&
            description="Scaling of the coulomb term in sic [experimental]",&
            usage="SIC_SCALING_A 0.5",&
            citations=(/VandeVondele2005b/),&
            default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SIC_SCALING_B",&
            description="Scaling of the xc term in sic [experimental]",&
            usage="SIC_SCALING_B 0.5",&
            citations=(/VandeVondele2005b/),&
            default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SIC_METHOD",&
            description="Method used to remove the self interaction",&
            usage="SIC_METHOD MAURI_US",&
            default_i_val=sic_none,&
            enum_c_vals=s2a( "NONE", "MAURI_US", "MAURI_SPZ", "AD", "EXPLICIT_ORBITALS"),&
            enum_i_vals=(/sic_none,sic_mauri_us,sic_mauri_spz,sic_ad,sic_eo/),&
            enum_desc=s2a("Do not apply a sic correction",&
                          "Employ a (scaled) correction proposed by Mauri and co-workers"//&
                                 " on the spin density / doublet unpaired orbital",&
                          "Employ a (scaled) Perdew-Zunger expression"//&
                                 " on the spin density / doublet unpaired orbital",&
                          "The average density correction",&
                          "(scaled) Perdew-Zunger correction explicitly on a set of orbitals."),&
            citations=(/VandeVondele2005b,Perdew1981,Avezac2005/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ORBITAL_SET",&
            description="Type of orbitals treated with the SIC",&
            usage="ORBITAL_SET ALL",&
            default_i_val=sic_list_unpaired,&
            enum_c_vals=s2a("UNPAIRED","ALL"),&
            enum_desc=s2a("correction for the unpaired orbitals only, requires a restricted open shell calculation",&
                          "correction for all orbitals, requires a LSD or ROKS calculation"),&
            enum_i_vals=(/sic_list_unpaired,sic_list_all/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_sic_section

! *****************************************************************************
!> \brief creates the low spin roks section
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE create_low_spin_roks_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_low_spin_roks_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"LOW_SPIN_ROKS",&
            description="Specify the details of the low spin ROKS method."//&
                        "In particular, one can specify various terms added to the energy of the high spin roks configuration"//&
                        " with a energy scaling factor, and a prescription of the spin state.",&
            n_keywords=6, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword,name="ENERGY_SCALING",&
            description="The scaling factors for each term added to the total energy."//&
                        "This list should contain one number for each term added to the total energy.",&
            usage="ENERGY_SCALING 1.0 -1.0 ",&
            n_var=-1,type_of_var=real_t,repeats=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,name="SPIN_CONFIGURATION",&
            description="for each singly occupied orbital, specify if this should be an alpha (=1) or a beta (=2) orbital"//&
                        "This keyword should be repeated, each repetition corresponding to an additional term." ,&
            usage="SPIN_CONFIGURATION 1 2",&
            n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_low_spin_roks_section

! *****************************************************************************
!> \brief Creates the section for applying an electrostatic external potential
!> \author teo
!> \date 12.2009
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
            description="Section controlling the presence of an electrostatic "//&
            "external potential dependent on the atomic positions (X,Y,Z)",&
            n_keywords=7, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="FUNCTION",&
            description="Specifies the functional form in mathematical notation. Variables must be the atomic "//&
            "coordinates (X,Y,Z) of the grid.",usage="FUNCTION  X^2+Y^2+Z^2+LOG(ABS(X+Y))", required=.TRUE.,&
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
            description="Defines the values of parameter of the functional form",&
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

       CALL keyword_create(keyword, name="STATIC",&
            description="Specifies the external potential as STATIC or time dependent. At the moment "//&
                        "only static potentials are implemented.",&
            usage="STATIC T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DX",&
            description="Parameter used for computing the derivative with the Ridders method.",&
            usage="DX <REAL>", default_r_val=0.1_dp, unit_str="bohr", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ERROR_LIMIT",&
            description="Checks that the error in computing the derivative is not larger than "//&
            "the value set. In case prints a warning message.",&
            usage="ERROR_LIMIT <REAL>", default_r_val=1.0E-12_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_ext_pot_section

! *****************************************************************************
!> \brief ZMP Creates the section for reading user supplied external density
!> \author D. Varsano [daniele.varsano@nano.cnr.it] 
!> \date 03.2011
! *****************************************************************************
  SUBROUTINE create_ext_den_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ext_den_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="EXTERNAL_DENSITY",&
            description="Section for the use of the ZMP technique on external densities.",&
            n_keywords=4, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            citations=(/Zhao1994, Tozer1996/),&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="FILE_DENSITY",&
            description="Specifies the filename containing the target density in *.cube format."//&
                        "In the MGRID section it must be imposed NGRID 1, as it works with only"//&
                        "one grid. The number of points in each direction, and the spacing must"//&
                        "be previously defined choosing the plane waves cut-off in section MGRID"//& 
                        "keyword CUTOFF, and the cube dimention in section SUBSYS / CELL / keyword ABC",&
            usage="DENSITY_FILE_NAME <FILENAME>",&
            type_of_var=char_t,default_c_val="RHO_O.dat", n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LAMBDA",&
            description="Lagrange multiplier defined in the constraint ZMP method. When starting, use"//&
                        "small values when starting from scratch (around 5,10). Then gradually increase"//&
                        "the values depending, restarting from the previous calculation with the smaller"//&
                        "value. To choose the progressive values of LAMBDA look at the convergence of the"//&
                        " eigenvalues.",&
            usage="DX <REAL>", default_r_val=10.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ZMP_CONSTRAINT",&
            description="Specify which kind of constraint to solve the ZMP equation. The COULOMB default"//&
                        "option is more stable.",&
            usage="ZMP_CONSTRAINT <CHAR>",&
            enum_c_vals=s2a("COULOMB", "DIFF", "NONE"),&
            enum_i_vals=(/use_coulomb, use_diff, use_no/),&
            enum_desc=s2a("Coulomb constraint, integral of [rho_0(r)-rho(r)]/|r-r'|",&
                          "Simple constraint, [rho_0(r)-rho(r)]",&
                          "No constrain imposed"),&
            default_i_val=use_coulomb, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FERMI_AMALDI",&
            description="Add the Fermi-Amaldi contribution to the Hartree potential."//&
                        "It leads to a more stable convergence.",&
            usage="FERMI_AMALDI <LOGICAL>",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_ext_den_section

! *****************************************************************************
!> \brief ZMP Creates the section for creating the external v_xc 
!> \author D. Varsano [daniele.varsano@nano.cnr.it] 
!> \date 03.2011
! *****************************************************************************
  SUBROUTINE create_ext_vxc_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ext_vxc_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="EXTERNAL_VXC",&
            description="SCF convergence with external v_xc calculated through previous ZMP"//&
                        "calculation",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="FILE_VXC",&
            description="The *.cube filename containing the v_xc potential. This works only"//&
                        "with NGRID 1 imposed in the MGRID section. The number of points in each"//&
                        "direction, and the spacing must equal to those previously used in the ZMP"//&
                        "calculation and defined through the plane wave cut-off and the cube dimension"//&
                        "respectively set in section MGRID / keyword CUTOFF, and in section SUBSYS /"//&
                        "CELL / keyword ABC",&
            usage="FILE_VXC <FILENAME>",&
            type_of_var=char_t,default_c_val="VXC_O.dat", n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_ext_vxc_section


! *****************************************************************************
!> \brief creates the section for static periodic fields
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Florian Schiffmann
! *****************************************************************************
  SUBROUTINE create_per_efield_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_per_efield_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"PERIODIC_EFIELD",&
            description="parameters for finite periodic electric field computed using "//&
            "the Berry phase approach. IMPORTANT: Can only be used in combination "//&
            " with OT. Can not be used in combination with RTP or EMD.",&
            n_keywords=6, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="INTENSITY",&
            description="Intensity of the electric field in a.u",&
            usage="INTENSITY  0.001",&
            default_r_val=0._dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="POLARISATION",&
            description="Polarisation vector of electric field",&
            usage="POLARISIATION  0.0 0.0 1.0",&
            repeats=.FALSE.,required=.FALSE.,n_var=3,&
            type_of_var=real_t,default_r_vals=(/0.0_dp,0.0_dp,1.0_dp/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_per_efield_section
! *****************************************************************************
!> \brief creates the section for time dependent nonperiodic fields
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Florian Schiffmann
! *****************************************************************************
  SUBROUTINE create_efield_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_efield_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"EFIELD",&
            description="parameters for finite, time dependent, nonperiodic electric fields",&
            n_keywords=6, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword,subsection)

       CALL keyword_create(keyword, name="INTENSITY",&
            description="Intensity of the electric field in W*cm-2 which corresponds "//&
            "to a maximal amplitude in a.u. of sqrt(I/(3.50944*10^16))",&
            usage="INTENSITY  0.001",&
            default_r_val=0._dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="POLARISATION",&
            description="Polarisation vector of electric field",&
            usage="POLARISATION  0.0 0.0 1.0",&
            repeats=.FALSE.,required=.FALSE.,n_var=3,&
            type_of_var=real_t,default_r_vals=(/0.0_dp,0.0_dp,1.0_dp/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="wavelength",&
            description="Wavelength of efield field",&
            usage="Wavelength  1.E0",&
            default_r_val=0._dp,unit_str="nm",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PHASE",&
            description="phase offset of the cosine given in multiples of pi",&
            usage="Phase  1.E0",&
            default_r_val=0._dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ENVELOP",&
            description="Shape of the efield pulse",&
            usage="ENVELOP CONSTANT",&
            default_i_val=constant_env,&
            enum_c_vals=s2a( "CONSTANT", "GAUSSIAN","RAMP" ),&
            enum_desc=s2a("No envelop function is applied to the strength",&
                          "A Gaussian function is used as envelop ",&
                          "Linear tune in/out of the field"),&
            enum_i_vals=(/constant_env,gaussian_env,ramp_env/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_constant_env_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_ramp_env_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_gaussian_env_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_efield_section

! *****************************************************************************
!> \brief makes the orbital transformation section
!> \par History
!>      11.2004 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE create_ot_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ot_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"OT",&
            description="Sets the various options for the orbital transformation (OT) method. "//&
            "Default settings already provide an efficient, yet robust method. "//&
            "Most systems benefit from using the FULL_ALL preconditioner "//&
            "combined with a small value (0.001) of ENERGY_GAP."//&
            "Well-behaved systems might benefit from using a DIIS minimizer.",&
            n_keywords=27, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            citations=(/VandeVondele2003,Weber2008/), error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the ot method",&
            usage="&OT T",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALGORITHM",&
            description="Algorithm to be used for OT",&
            usage="ALGORITHM STRICT",&
            default_i_val=ot_algo_taylor_or_diag,&
            enum_c_vals=s2a( "STRICT", "IRAC" ),&
            enum_desc=s2a("Strict orthogonality: Taylor or diagonalization based algorithm.",&
                          "Orbital Transformation based Iterative Refinement "//&
                          "of the Approximative Congruence transformation (OT/IR)."),&
            enum_i_vals=(/ot_algo_taylor_or_diag,ot_algo_irac/),&
            citations=(/VandeVondele2003,VandeVondele2005a,Weber2008/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="IRAC_DEGREE",&
            description="The refinement polynomial degree (2, 3 or 4).",&
            usage="IRAC_DEGREE 4",&
            default_i_val=4,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_IRAC",&
            description="Maximum allowed refinement iteration.",&
            usage="MAX_IRAC 5",&
            default_i_val=50,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MIXED_PRECISION",&
            description="Uses a mixed precision algorithm." //&
            "With a well behaved basis set (i.e. condition number less than 1/eps_sp)"//&
            "it provides double precision accuracy results and up to a 2 fold speedup for building and "//&
            "applying the preconditioner.",&
            usage="MIXED_PRECISION T",&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ORTHO_IRAC",&
            description="The orthogonality method.",&
            usage="ORTHO_IRAC POLY",&
            default_i_val=ot_chol_irac,&
            enum_c_vals=s2a( "CHOL", "POLY", "LWDN"),&
            enum_desc=s2a("Cholesky.","Polynomial.","Loewdin."),&
            enum_i_vals=(/ot_chol_irac,ot_poly_irac,ot_lwdn_irac/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_IRAC_FILTER_MATRIX",&
            description="Sets the threshold for filtering the matrices.",&
            usage="EPS_IRAC_FILTER_MATRIX 1.0E-5",&
            default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_IRAC",&
            description="Targeted accuracy during the refinement iteration.",&
            usage="EPS_IRAC 1.0E-5",&
            default_r_val=1.0E-10_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_IRAC_QUICK_EXIT",&
            description="Only one extra refinement iteration is "//&
            "done when the norm is below this value.",&
            usage="EPS_IRAC_QUICK_EXIT 1.0E-2",&
            default_r_val=1.0E-5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_IRAC_SWITCH",&
            description="The algorithm switches to the polynomial "//&
            "refinement when the norm is below this value.",&
            usage="EPS_IRAC_SWITCH 1.0E-3",&
            default_r_val=1.0E-2_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ON_THE_FLY_LOC",&
            description="On the fly localization of the molecular orbitals. "//&
            "Can only be used with OT/IRAC.",&
            usage="ON_THE_FLY_LOC T",&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MINIMIZER",&
            description="Minimizer to be used with the OT method",&
            usage="MINIMIZER DIIS",&
            default_i_val=ot_mini_cg,&
            enum_c_vals=s2a( "SD", "CG", "DIIS","BROYDEN"),&
            enum_desc=s2a("Steepest descent: not recommended","Conjugate Gradients: most reliable, use for difficult systems."//&
            " The total energy should decrease at every OT CG step if the line search is appropriate.", &
            "Direct inversion in the iterative subspace: less reliable than CG, but sometimes about 50% faster",&
            "Broyden mixing approximating the inverse Hessian"),&
            enum_i_vals=(/ot_mini_sd,ot_mini_cg,ot_mini_diis,ot_mini_broyden/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SAFE_DIIS",&
            variants=(/"SAFER_DIIS"/),&
            description="Reject DIIS steps if they point away from the"//&
            " minimum, do SD in that case.",&
            usage="SAFE_DIIS ON", default_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="N_HISTORY_VEC",&
            variants=s2a("NDIIS","N_DIIS","N_BROYDEN"),&
            description="Number of history vectors to be used with DIIS or BROYDEN",&
            usage="N_DIIS 4",&
            default_i_val=7,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BROYDEN_BETA",&
            description="Underrelaxation for the broyden mixer",&
            usage="BROYDEN_BETA 0.9",&
            default_r_val=0.9_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BROYDEN_GAMMA",&
            description="Backtracking parameter",&
            usage="BROYDEN_GAMMA 0.5",&
            default_r_val=0.5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BROYDEN_SIGMA",&
            description="Curvature of energy functional.",&
            usage="BROYDEN_SIGMA 0.25",&
            default_r_val=0.25_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BROYDEN_ETA",&
            description="Dampening of estimated energy curvature.",&
            usage="BROYDEN_ETA 0.7",&
            default_r_val=0.7_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BROYDEN_OMEGA",&
            description="Growth limit of curvature.",&
            usage="BROYDEN_OMEGA 1.1",&
            default_r_val=1.1_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BROYDEN_SIGMA_DECREASE",&
            description="Reduction of curvature on bad approximation.",&
            usage="BROYDEN_SIGMA_DECREASE 0.7",&
            default_r_val=0.7_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BROYDEN_SIGMA_MIN",&
            description="Minimum adaptive curvature.",&
            usage="BROYDEN_SIGMA_MIN 0.05",&
            default_r_val=0.05_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BROYDEN_FORGET_HISTORY",&
            description="Forget history on bad approximation", &
            usage="BROYDEN_FORGET_HISTORY OFF", default_l_val=.FALSE., &
            lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BROYDEN_ADAPTIVE_SIGMA",&
            description="Enable adaptive curvature estimation", &
            usage="BROYDEN_ADAPTIVE_SIGMA ON", default_l_val=.TRUE., &
            lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BROYDEN_ENABLE_FLIP",&
            description="Ensure positive definite update", &
            usage="BROYDEN_ENABLE_FLIP ON", default_l_val=.TRUE., &
            lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LINESEARCH",&
            variants=(/"LINE_SEARCH"/),&
            description="1D line search algorithm to be used with the OT minimizer,"//&
            " in increasing order of robustness and cost. MINIMIZER CG combined with"//&
            " LINESEARCH GOLD should always find an electronic minimum. "//&
            " Whereas the 2PNT minimizer is almost always OK, 3PNT might be needed for systems"//&
            " in which successive OT CG steps do not decrease the total energy.",&
            usage="LINESEARCH GOLD",&
            default_i_val=ls_2pnt,&
            enum_c_vals=s2a( "NONE", "2PNT", "3PNT","GOLD"),&
            enum_desc=s2a("take fixed lenght steps","extrapolate based on 2 points", &
            "... or on 3 points","perform 1D golden section search of the minimum (very expensive)"),&
            enum_i_vals=(/ls_none,ls_2pnt,ls_3pnt,ls_gold/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STEPSIZE",&
            description="Initial stepsize used for the line search, sometimes this parameter can be reduced to stablize DIIS"//&
            " or to improve the CG behavior in the first few steps",&
            usage="STEPSIZE 0.4",&
            default_r_val=0.15_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GOLD_TARGET",&
            description="Target relative uncertainty in the location of the minimum for LINESEARCH GOLD",&
            usage="GOLD_TARGET 0.1",&
            default_r_val=0.01_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PRECONDITIONER",&
            description="Type of preconditioner to be used with all minimization schemes. "//&
            "They differ in effectiveness, cost of construction, cost of application. "//&
            "Properly preconditioned minimization can be orders of magnitude faster than doing nothing.",&
            usage="PRECONDITIONER FULL_ALL",&
            default_i_val=ot_precond_full_kinetic,&
            enum_c_vals=s2a("FULL_ALL","FULL_SINGLE_INVERSE","FULL_SINGLE","FULL_KINETIC","FULL_S_INVERSE",&
            "SPARSE_DIAG","SPARSE_KINETIC",&
            "SPARSE_STABLE_APPROXIMATE_INVERSE",&
            "NONE"),&
            enum_desc=s2a("Most effective state selective preconditioner based on diagonalization, "//&
            "requires the ENERGY_GAP parameter to be an underestimate of the HOMO-LUMO gap. "//&
            "This preconditioner is recommended for almost all systems, except very large systems where "//&
            "make_preconditioner would dominate the total computational cost.",&
            "Based on H-eS cholesky inversion, similar to FULL_SINGLE in preconditioning efficiency "//&
            "but cheaper to construct, "//&
            "might be somewhat less robust. Recommended for large systems.",&
            "Based on H-eS diagonalisation, not as good as FULL_ALL, but somewhat cheaper to apply. ",&
            "Cholesky inversion of S and T, fast construction, robust, and relatively good, "//&
            "use for very large systems.",&
            "Cholesky inversion of S, not as good as FULL_KINETIC, yet equally expensive.",&
            "Only based on atomic blocks, very cheap in construction but quite poor.",&
            "Toward -sparse- linear solver for T-eS.",&
            "Stabilized Approximate INVerse preconditioning using DBCSR sparse matrices (not working yet).",&
            "skip preconditioning"),&
            enum_i_vals=(/ot_precond_full_all,ot_precond_full_single_inverse,ot_precond_full_single, &
            ot_precond_full_kinetic,ot_precond_s_inverse,&
            ot_precond_sparse_diag,ot_precond_sparse_kinetic,&
            ot_precond_sparse_kinetic_sainv,&
            ot_precond_none/),&
            citations=(/VandeVondele2003/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PRECOND_SOLVER",&
            description="How the preconditioner is applied to the residual.",&
            usage="PRECOND_SOLVER DIRECT",&
            default_i_val=ot_precond_solver_default,&
            enum_c_vals=s2a( "DEFAULT", "DIRECT", "INVERSE_CHOLESKY","SAINV"),&
            enum_desc=s2a("the default","Cholesky decomposition followed by triangular solve "//&
            "(works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)",&
            "Cholesky decomposition followed by explicit inversion "//&
            "(works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)",&
            "SAINV inversion (works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)"),&
            enum_i_vals=(/ot_precond_solver_default,&
            ot_precond_solver_direct,&
            ot_precond_solver_inv_chol,&
            ot_precond_solver_sainv/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ENERGY_GAP",&
            description="Should be an estimate for the energy gap [a.u.] (HOMO-LUMO) and is used in preconditioning, "//&
            "especially effective with the FULL_ALL preconditioner, in which case it should be an underestimate "//&
            "of the gap (0.001 doing normally fine). For the other preconditioners, making this value larger (0.2)"//&
            " will tame the preconditioner in case of poor initial guesses.",&
            usage="ENERGY_GAP 0.001",&
            default_r_val=0.2_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_TAYLOR",&
            variants=(/"EPSTAYLOR"/),&
            description="Target accuracy of the taylor expansion for the matrix functions, should normally be kept as is.",&
            usage="EPS_TAYLOR 1.0E-15",&
            default_r_val=1.0E-16_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_TAYLOR",&
            description="Maximum order of the Taylor expansion before diagonalisation is prefered, for large parallel runs"//&
            " a slightly higher order could sometimes result in a small speedup.",&
            usage="MAX_TAYLOR 5",&
            default_i_val=4,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ROTATION",&
            description="Introduce additional variables so that rotations of the occupied"//&
            " subspace are allowed as well, only needed for cases where the energy is not invariant under "//&
            " a rotation of the occupied subspace such as non-singlet restricted calculations "//&
            " or fractional occupations.",&
            usage="ROTATION",lone_keyword_l_val=.TRUE.,&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SCP_NDDO",&
            description="Introduce additional self-consistent polarization through"//&
            " response basis set = orbital basis set for NDDO.)",&
            usage="SCP_NDDO",lone_keyword_l_val=.TRUE.,&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ENERGIES",&
            description="Optimize orbital energies for use in Fermi-Dirac smearing "//&
                        "(requires ROTATION and FD smearing to be active).",&
            usage="ENERGIES",lone_keyword_l_val=.TRUE.,&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="OCCUPATION_PRECONDITIONER",&
            description="Preconditioner with the occupation numbers (FD smearing)",&
            usage="OCCUPATION_PRECONDITIONER",lone_keyword_l_val=.TRUE.,&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NONDIAG_ENERGY",&
            description="Add a non-diagonal energy penalty (FD smearing)",&
            usage="NONDIAG_ENERGY",lone_keyword_l_val=.TRUE.,&
            default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NONDIAG_ENERGY_STRENGTH",&
            description="The prefactor for the non-diagonal energy penalty (FD smearing)",&
            usage="NONDIAG_ENERGY_STRENGTH", default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


    END IF
  END SUBROUTINE create_ot_section

! *****************************************************************************
!> \brief creates the diagonalization section
!> \par History
!>      10.2008 created [JGH]
! *****************************************************************************
  SUBROUTINE create_diagonalization_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_diagonalization_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"DIAGONALIZATION",&
            description="Set up type and parameters for Kohn-Sham matrix diagonalization.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the diagonalization method",&
            usage="&DIAGONALIZATION T",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALGORITHM",&
            description="Algorithm to be used for diagonalization",&
            usage="ALGORITHM STANDARD",&
            default_i_val=diag_standard,&
            enum_c_vals=s2a( "STANDARD","OT", "LANCZOS","DAVIDSON" ),&
            enum_desc=s2a("Standard diagonalization: LAPACK methods or Jacobi.",&
                          "Iterative diagonalization using OT method",&
                          "Block Krylov-space approach to self-consistent diagonalisation",&
                          "Preconditioned blocked Davidson"),&
            enum_i_vals=(/diag_standard, diag_ot, diag_block_krylov, diag_block_davidson/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="JACOBI_THRESHOLD",&
            description="Controls the accuracy of the pseudo-diagonalization method using Jacobi rotations",&
            usage="JACOBI_THRESHOLD 1.0E-6",&
            default_r_val=1.0E-7_dp,&
            citations=(/Stewart1982/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_JACOBI",&
            description="Below this threshold value for the SCF convergence the pseudo-diagonalization "//&
                        "method using Jacobi rotations is activated. This method is much faster than a "//&
                        "real diagonalization and it is even speeding up while achieving full convergence."//&
                        "However, it needs a pre-converged wavefunction obtained by at least one real "//&
                        "diagonalization which is further optimized while keeping the original eigenvalue "//&
                        "spectrum. The MO eigenvalues are NOT updated. The method might be useful to speed "//&
                        "up calculations for large systems e.g. using a semi-empirical method.",&
            usage="EPS_JACOBI 1.0E-5",&
            default_r_val=0.0_dp,&
            citations=(/Stewart1982/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_ADAPT",&
            description="Required accuracy in iterative diagonalization as compared to current SCF convergence",&
            usage="EPS_ADAPT 0.01",&
            default_r_val=0._dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_ITER",&
            description="Maximum number of iterations in iterative diagonalization",&
            usage="MAX_ITER 20",&
            default_i_val=2,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_ITER",&
            description="Required accuracy in iterative diagonalization",&
            usage="EPS_ITER 1.e-8",&
            default_r_val=1.e-8_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       NULLIFY(subsection)
       CALL create_ot_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)


       NULLIFY(subsection)
       CALL create_krylov_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_diag_subspace_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_davidson_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_diagonalization_section


  SUBROUTINE create_davidson_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_davidson_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"DAVIDSON",&
            description=" ",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="PRECONDITIONER",&
            description="Type of preconditioner to be used with all minimization schemes. ",&
            usage="PRECONDITIONER FULL_ALL",&
            default_i_val=ot_precond_full_all,&
            enum_c_vals=s2a("FULL_ALL","FULL_SINGLE_INVERSE","NONE"),&
            enum_desc=s2a("Most effective state selective preconditioner based on diagonalization ",&
            "Based on H-eS cholesky inversion, similar to FULL_SINGLE in preconditioning efficiency "//&
            "but cheaper to construct, might be somewhat less robust. Recommended for large systems.",&
            "skip preconditioning"),&
            enum_i_vals=(/ot_precond_full_all,ot_precond_full_single_inverse,ot_precond_none/),&
            citations=(/VandeVondele2003/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PRECOND_SOLVER",&
            description="How the preconditioner is applied to the residual.",&
            usage="PRECOND_SOLVER DIRECT",&
            default_i_val=ot_precond_solver_default,&
            enum_c_vals=s2a( "DEFAULT", "DIRECT", "INVERSE_CHOLESKY","SAINV"),&
            enum_desc=s2a("the default","Cholesky decomposition followed by triangular solve "//&
            "(works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)",&
            "Cholesky decomposition followed by explicit inversion "//&
            "(works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)",&
            "SAINV inversion (works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)"),&
            enum_i_vals=(/ot_precond_solver_default,&
            ot_precond_solver_direct,&
            ot_precond_solver_inv_chol,&
            ot_precond_solver_sainv/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ENERGY_GAP",&
            description="Should be an estimate for the energy gap [a.u.] (HOMO-LUMO) and is used in preconditioning, "//&
            "especially effective with the FULL_ALL preconditioner, in which case it should be an underestimate "//&
            "of the gap (0.001 doing normally fine). For the other preconditioners, making this value larger (0.2)"//&
            " will tame the preconditioner in case of poor initial guesses.",&
            usage="ENERGY_GAP 0.001",&
            default_r_val=0.2_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NEW_PREC_EACH",&
            description="Number of SCF iterations after which a new Preconditioner is computed",&
            usage="NEW_PREC_EACH 10", default_i_val=20,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FIRST_PREC",&
            description="First SCF iteration at which a Preconditioner is employed",&
            usage="FIRST_PREC 1", default_i_val=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CONV_MOS_PERCENT",&
            description="Minimal percent of MOS that have to converge within the Davidson loop"//&
            " before the SCF iteration is completed and a new Hamiltonian is computed",&
            usage="CONV_MOS_PERCENT 0.8", default_r_val=0.5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SPARSE_MOS",&
            description="Use MOS as sparse matrix and avoid as much as possible multiplications with full matrices",&
            usage="SPARSE_MOS",default_l_val=.TRUE.,&
            lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_davidson_section

  SUBROUTINE create_krylov_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_krylov_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"KRYLOV",&
            description=" ",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="NKRYLOV",&
            description="Dimension of the Krylov space used for the Lanczos refinement",&
            usage="NKRYLOV 20",&
            default_i_val=4,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NBLOCK",&
            description="Size of the block of vectors refined simultaneously by the Lanczos procedure",&
            usage="NBLOCK 1",&
            default_i_val=32,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_KRYLOV",&
            description="Convergence criterion for the MOs",&
            usage="EPS_KRYLOV 0.00001",&
            default_r_val=0.0000001_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_STD_DIAG",&
            description="Level of convergence to be reached before starting the Lanczos procedure."//&
            " Above this threshold a standard diagonalization method is used. "//&
            " If negative Lanczos is started at the first iteration",&
            usage="EPS_STD_DIAG 0.001",&
            default_r_val=-1.0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHECK_MOS_CONV",&
            description="This requires to check the convergence of MOS also when standard "//&
            "diagonalization steps are performed, if the block krylov approach is active.",&
            usage="CHECK_MOS_CONV T",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


    END IF
  END SUBROUTINE create_krylov_section


  SUBROUTINE create_diag_subspace_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_diag_subspace_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"DIAG_SUB_SCF",&
            description="Activation of self-consistenf subspace refinement by diagonalization "//&
            "of H by adjusting the occupation but keeping the MOS unchanged.",&
            n_keywords=2, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of inner SCF loop to refine occupations in MOS subspace",&
            usage="&DIAG_SUB_SCF T",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_ITER",&
            description="Maximum number of iterations for the SCF inner loop",&
            usage="MAX_ITER 20",&
            default_i_val=2,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_ENE",&
            description="Required energy accuracy for convergence of subspace diagonalization",&
            usage="EPS_ENE 1.e-8",&
            default_r_val=1.e-4_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_ADAPT_SCF",&
            description="Required density matrix accuracy as compared to current SCF convergence",&
            usage="EPS_ADAPT_SCF 1.e-1",&
            default_r_val=1._dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_SKIP_SUB_DIAG",&
            description="Level of convergence to be reached before starting the internal loop of subspace rotations."//&
            " Above this threshold only the outer diagonalization method is used. "//&
            " If negative the subspace rotation is started at the first iteration",&
            usage="EPS_SKIP_SUB_DIAG 0.001",&
            default_r_val=-1.0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_mixing_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
    END IF
  END SUBROUTINE create_diag_subspace_section

! *****************************************************************************
!> \brief makes the input section for core-level spectroscopy simulations
!> \par History
!>      03.2005 created [MI]
! *****************************************************************************
  SUBROUTINE create_xas_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_xas_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"xas",&
            description="Sets the method of choice to calculate core-level excitation spectra. "//&
            "The occupied states from  which we calculate the "//&
            "excitation should be specified. "//&
            "Localization of the orbitals may be useful.",&
            n_keywords=10, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            citations=(/Iannuzzi2007/),&
            error=error)

       NULLIFY(keyword,subsection,print_key)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of core-level spectroscopy simulations",&
            usage="&XAS T",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="METHOD",&
            variants=(/"XAS_METHOD"/),&
            description="Method to be used to calculate core-level excitation spectra",&
            usage="METHOD TP_HH",&
            default_i_val=xas_none,&
            enum_c_vals=s2a( "NONE", "TP_HH", "TP_FH", "TP_VAL", "TP_XHH","TP_XFH", "DSCF"),&
            enum_desc=s2a(&
            "No core electron spectroscopy", "Transition potential half-hole",&
            "Transition potential full-hole", "Hole in homo for X-ray emission only " , &
            "Transition potential excited half-hole", &
            "Transition potential excited full-hole " , &
            "DSCF calculations to compute the first (core)excited state"),&
            enum_i_vals=(/xas_none,xas_tp_hh,xas_tp_fh,xes_tp_val,xas_tp_xhh,&
            xas_tp_xfh,xas_dscf/), error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="XES_CORE",&
            description="occupation of the core state in XES calculation by TP_VAL."//&
            "The homo is emptied by the same amount",&
            usage="XES_CORE 0.5",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DIPOLE_FORM",&
            variants=(/"DIP_FORM"/),&
            description="Type of integral to get the oscillator strengths "//&
            "in the diipole approximation",&
            usage="DIPOLE_FORM string",&
            default_i_val=xas_dip_vel,&
            enum_c_vals=s2a( "LENGTH","VELOCITY" ),&
            enum_desc=s2a("Length form &lang; i | e r | j &rang;",&
                          "Velocity form &lang; i | d/dr | j &rang;"),&
            enum_i_vals=(/xas_dip_len2,xas_dip_vel/),error=error)
       CALL section_add_keyword(section,keyword,error=error)

       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="SCF_OPTIMIZER",&
            description="Optimization algorithm: diagonalization or OT",&
            usage="SCF_OPTIMIZER GENERAL",&
            default_i_val=xas_scf_general,&
            enum_c_vals=s2a( "DEFAULT", "GENERAL"),&
            enum_desc=s2a("same as in std SCF", "diagonalization"),&
            enum_i_vals=(/xas_scf_default,xas_scf_general/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAXSTEP",&
            description="Max # of steps in the cls-scf for one excitation",&
            usage="MAXSTEP 150",&
            default_i_val=150,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CONVERGENCE",&
            variants=(/"CONV"/),&
            description="Convergence criterion for the xas-scf",&
            usage="CONVERGENCE 0.00005",&
            default_r_val=0.5E-6_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_DIIS",&
            description="treshold on the convergence to start"//&
            " using DIAG/DIIS for the cls-scf"//&
            "if default, the scf_control value is used",&
            usage="EPS_DIIS 0.5",&
            default_r_val=-1._dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STATE_TYPE",&
            variants=(/"TYPE"/),&
            description="Type of the orbitas that are excited for the xas spectra calculation",&
            usage="STATE_TYPE 1S",&
            default_i_val=xas_1s_type,&
            enum_c_vals=s2a("1S","2S","2P"),&
            enum_desc=s2a("1s orbitals","2s orbitals","2p orbitals"),&
            enum_i_vals=(/xas_1s_type,xas_2s_type,xas_2p_type/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STATE_SEARCH",&
            description="# of states where to look for the one to be excited",&
            usage="STATE_SEARCH 1",&
            default_i_val=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="ATOMS_LIST",&
            variants=(/"AT_LIST"/),&
            description="Indexes of the atoms to be excited"//&
            "This keyword can be repeated several times"//&
            "(useful if you have to specify many indexes).",&
            usage="ATOMS_LIST {integer}  {integer} ..  {integer} ",&
            n_var=-1,type_of_var=integer_t,required=.TRUE.,repeats=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ADDED_MOS",&
            description="Number of additional MOS added spin up only",&
            usage="ADDED_MOS {integer}", default_i_val=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_ITER_ADDED",&
            description="maximum number of iteration in calculation of added orbitals",&
            usage="MAX_ITER_ADDED 100", default_i_val=2999,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_ADDED",&
            description="target accuracy incalculation of the added orbitals",&
            usage="EPS_ADDED 1.e-6", default_r_val=1.0e-5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NGAUSS",&
            description="Number of gto's for the expansion of the sto"//&
            "of the type given by STATE_TYPE",&
            usage="NGAUSS {integer}", default_i_val=3,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART",&
            description="Restart the excited state if the restart file exists",&
            usage="RESTART",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WFN_RESTART_FILE_NAME",&
            variants=(/"RESTART_FILE_NAME"/),&
            description="Root of the file names where to read the MOS from"//&
            "which to restart the calculation of the core level excited states",&
            usage="WFN_RESTART_FILE_NAME <FILENAME>",&
            type_of_var=lchar_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_localize_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,"PRINT",&
            "printing of information during the core-level spectroscopy simulation",&
            error=error,repeats=.FALSE.,required=.FALSE.)

       ! Add printing of wannier infos
       CALL print_wanniers(subsection, error)

       CALL cp_print_key_section_create(print_key,"iteration_info",&
            description="Controls the printing of basic iteration information during the xas scf.", &
            print_level=low_print_level,filename="__STD_OUT__",&
            error=error)
       CALL keyword_create(keyword, name="time_cumul",&
            description="If the printkey is activated switches the printing of timings"//&
            " to cumulative (over the scf).",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing of basic iteration information in CLS", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"XES_SPECTRUM",&
            description="Controls the dumping of the CLS output files containing the emission spectra",&
            print_level=low_print_level,common_iter_levels=3,filename="",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"XAS_SPECTRUM",&
            description="Controls the dumping of the CLS output files containing the absorption spectra",&
            print_level=low_print_level,common_iter_levels=3,filename="",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"RESTART",&
            description="Controls the dumping of MO restart file during the scf"//&
            "of a Core-Level-Spectroscopy calculation. For each new excited atom,"//&
            "one different restart file is dumped. These restart files should be"//&
            "employed only to restart the same type of CLS calculation, "//&
            "i.e. with the same core potential.", &
            print_level=low_print_level,common_iter_levels=3,each_iter_names=s2a("XAS_SCF"),&
            add_last=add_last_numeric,each_iter_values=(/3/),filename="",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CLS_FUNCTION_CUBES",&
            description="Controls the printing of the relaxed orbitals ", &
            print_level=high_print_level,common_iter_levels=3,add_last=add_last_numeric,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components.",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="CUBES_LU_BOUNDS",&
            variants=(/"CUBES_LU"/),&
            description="The lower and upper index of the states to be printed as cube",&
            usage="CUBES_LU_BOUNDS integer integer",&
            n_var=2,default_i_vals=(/0,-2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="CUBES_LIST",&
            description="Indexes of the states to be printed as cube files"//&
            "This keyword can be repeated several times"//&
            "(useful if you have to specify many indexes).",&
            usage="CUBES_LIST 1 2",&
            n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="APPEND",&
            description="append the cube files when they already exist",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_mixing_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_smear_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_xas_section

! *****************************************************************************
!> \brief      Create CP2K input section for the smearing of occupation numbers
!> \author     Matthias Krack (MK)
!> \date       27.08.2008
!> \version    1.0
! *****************************************************************************
  SUBROUTINE create_smear_section(section,error)

    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_smear_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure = .FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

      CALL section_create(section,&
                          name="SMEAR",&
                          description="Define the smearing of the MO occupation numbers",&
                          n_keywords=6,&
                          n_subsections=0,&
                          repeats=.FALSE.,&
                          required=.FALSE.,&
                          error=error)

      NULLIFY (keyword)

      CALL keyword_create(keyword,&
                          name="_SECTION_PARAMETERS_",&
                          description="Controls the activation of smearing",&
                          usage="&SMEAR ON",&
                          default_l_val=.FALSE.,&
                          lone_keyword_l_val=.TRUE.,&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="METHOD",&
                          description="Smearing method to be applied",&
                          usage="METHOD Fermi_Dirac",&
                          default_i_val=smear_energy_window,&
                          enum_c_vals=s2a("FERMI_DIRAC","ENERGY_WINDOW","LIST"),&
                          enum_i_vals=(/smear_fermi_dirac,smear_energy_window,smear_list/),&
                          enum_desc=s2a("Fermi-Dirac distribution defined by the keyword ELECTRONIC_TEMPERATURE",&
                                        "Energy window defined by the keyword WINDOW_SIZE",&
                                        "Use a fixed list of occupations"),&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="LIST",&
                          description="A list of fractional occupations to use. Must match the number of states "//&
                                      "and sum up to the correct number of electrons",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=-1,&
                          type_of_var=real_t,&
                          usage="LIST 2.0 0.6666 0.6666 0.66666 0.0 0.0",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="ELECTRONIC_TEMPERATURE",&
                          variants=s2a("ELEC_TEMP","TELEC"),&
                          description="Electronic temperature in the case of Fermi-Dirac smearing",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=cp_unit_to_cp2k(value=300.0_dp, unit_str="K", error=error),&
                          unit_str="K",&
                          usage="ELECTRONIC_TEMPERATURE [K] 300",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="EPS_FERMI_DIRAC",&
                          description="Accuracy checks on occupation numbers use this as a tolerance",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=1.0E-10_dp,&
                          usage="EPS_FERMI_DIRAC 1.0E-6",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="WINDOW_SIZE",&
                          description="Size of the energy window centred at the Fermi level",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=0.0_dp,&
                          unit_str="au_e",&
                          usage="WINDOW_SIZE [eV] 0.3",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="FIXED_MAGNETIC_MOMENT",&
                          description="Imposed difference between the numbers of electrons of spin up "//&
                                      "and spin down: m = n(up) - n(down). A negative value (default) allows "//&
                                      "for a change of the magnetic moment. -1 specifically keeps an integer "//&
                                      "number of spin up and spin down electrons.",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=-100.0_dp,&
                          usage="FIXED_MAGNETIC_MOMENT 1.5",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_smear_section

! *****************************************************************************
!> \brief      Create CP2K input section for the mixing of the density matrix to
!>             be used only with diagonalization methods, i.e. not with OT
!> \author     MI
!> \date       20.02.2009
!> \version    1.0
! *****************************************************************************
  SUBROUTINE create_mixing_section(section,error)

    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_mixing_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure = .FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

      CALL section_create(section,&
                          name="MIXING",&
                          description="Define type and parameters for mixing"//&
                          "procedures to be applied to the density matrix. Normally, "//&
                          "only one type of mixing method should be accepted. The mixing "//&
                          "procedures activated by this section are only active for diagonalization "//&
                          "methods, i.e. not with minimization methods based on OT.",&
                          n_keywords=16,&
                          n_subsections=0,&
                          repeats=.FALSE.,&
                          required=.FALSE.,&
                          error=error)

      NULLIFY (keyword)

      CALL keyword_create(keyword,&
                          name="_SECTION_PARAMETERS_",&
                          description="Controls the activation of the mixing procedure",&
                          usage="&MIXING ON",&
                          default_l_val=.TRUE.,&
                          lone_keyword_l_val=.TRUE.,&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="METHOD",&
                          description="Mixing method to be applied",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          usage="METHOD KERKER_MIXING",&
                          default_i_val=direct_p_mix,&
                          enum_c_vals=s2a("NONE",&
                                          "DIRECT_P_MIXING",&
                                          "KERKER_MIXING",&
                                          "PULAY_MIXING",&
                                          "BROYDEN_MIXING",&
                                          "BROYDEN_MIXING_NEW",&
                                          "MULTISECANT_MIXING"),&
                          enum_i_vals=(/no_mix,direct_p_mix,kerker_mix,pulay_mix,broy_mix,&
                          broy_mix_new,multisec_mix/),&
                          enum_desc=s2a("No mixing is applied",&
                             "Direct mixing of new and old density matrices",&
                             "Mixing of the potential in reciprocal space using the Kerker damping",&
                             "Pulay mixing","Broyden mixing","Broyden mixing second version",&
                             "Multisecant scheme for mixing" ),&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="ALPHA",&
                          description="Fraction of new density to be included",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=0.4_dp,&
                          usage="ALPHA 0.2",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="BETA",&
                          description="Denominator parameter in Kerker damping "//&
                             "introduced to suppress charge sloshing: rho_mix(g) ="//&
                             "rho_in(g) + alpha*g^2/(g^2 + beta^2)*(rho_out(g)-"//&
                             "rho_in(g))",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=0.5_dp,&
                          unit_str="bohr^-1",&
                          usage="BETA 1.5",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="PULAY_ALPHA",&
                          description="Fraction of new density to be added to the Pulay expansion",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=0.0_dp,&
                          usage="PULAY_ALPHA 0.2",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="PULAY_BETA",&
                          description="Fraction of residual contribution to be added to Pulay expansion",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=1.0_dp,&
                          usage="PULAY_BETA 0.2",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="NMIXING",&
            description="Minimal number of density mixing (should be greater than 0),"//&
            "before starting DIIS",&
            usage="NMIXING 1", default_i_val=2,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="NBUFFER",&
                          variants=s2a("NPULAY","NBROYDEN","NMULTISECANT"),&
                          description="Number of previous steps stored for the actual mixing scheme",&
                          usage="NBUFFER 2", default_i_val=4,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="BROY_W0",&
                          description=" w0 parameter used in Broyden mixing",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=0.01_dp,&
                          usage="BROY_W0 0.03",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="BROY_WREF",&
                          description="",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=100.0_dp,&
                          usage="BROY_WREF 0.2",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="BROY_WMAX",&
                          description="",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=30.0_dp,&
                          usage="BROY_WMAX 10.0",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="REGULARIZATION",&
                          description="Regularization parameter to stabilize "//&
                          "the inversion of the residual matrix {Yn^t Yn} in the "//&
                          "multisecant mixing scheme (noise)",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=0.00001_dp,&
                          usage="REGULARIZATION 0.000001",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="MAX_STEP",&
                          description="Upper bound for the magnitude of the "//&
                          "unpredicted step size in the update by the "//&
                          "multisecant mixing scheme",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=0.1_dp,&
                          usage="MAX_STEP .2",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,&
                          name="R_FACTOR",&
                          description="Control factor for the magnitude of the "//&
                          "unpredicted step size in the update by the "//&
                          "multisecant mixing scheme",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=0.05_dp,&
                          usage="R_FACTOR .12",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="NSKIP",&
                          variants=(/"NSKIP_MIXING"/),&
            description="Number of initial iteration for which the mixing is skipped",&
            usage="NSKIP 10", default_i_val=0,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="N_SIMPLE_MIX",&
                          variants=(/"NSIMPLEMIX"/),&
            description="Number of kerker damping iterations before starting other mixing procedures",&
            usage="NSIMPLEMIX", default_i_val=0,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="KERKER_MIN",&
                          description="Minimal Kerker damping factor: MAX(g^2/(g^2 + beta^2),KERKER_MIN)",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=0.1_dp,&
                          usage="KERKER_MIN 0.05",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="MAX_GVEC_EXP",&
                          description="Restricts the G-space mixing to lower part of G-vector spectrum,"//&
                          " up to a G0, by assigning the exponent of the Gaussian that can be "//&
                          "represented by vectors smaller than G0 within a certain accuracy. ",&
                          required=.FALSE.,&
                          repeats=.FALSE.,&
                          n_var=1,&
                          type_of_var=real_t,&
                          default_r_val=-1._dp,&
                          usage="MAX_GVEC_EXP 3.",&
                          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_mixing_section

! *****************************************************************************
  SUBROUTINE create_rtp_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_rtp_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, print_section

    failure=.FALSE.
    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"REAL_TIME_PROPAGATION",&
            description="Parameters needed to set up the real time propagation"//&
            " for the electron dynamics",&
            n_keywords=4, n_subsections=4, repeats=.FALSE., required=.FALSE.,&
            citations=(/Kunert2003/),&
            error=error)

       CALL keyword_create(keyword, name="MAX_ITER",&
            description="Maximal number of iterations for the self consistent propagator loop.",&
            usage="MAX_ITER 10",&
            default_i_val=10,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_ITER",&
            description="Convergence criterium for the self consistent propagator loop.",&
            usage="EPS_ITER 1.0E-5",&
            default_r_val=1.0E-7_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ASPC_ORDER",&
            description="Speciefies how many steps will be used for extrapolation. "//&
            "One will be always used which is means X(t+dt)=X(t)",&
            usage="ASPC_ORDER 3",&
            default_i_val=3,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXTRAPOLATION",&
            description="Controls which quantity is extrapolated",&
            usage="EXTRAPOLATION MOS",default_i_val=extr_mos, &
            enum_c_vals=s2a("MOS","S_KS"),&
            enum_i_vals=(/extr_mos,extr_s_ks/),&
            enum_desc=s2a("Direct extrapolation on the MO matrices",&
            "Extrapolation on the Kohn-Sham matrix, resp. the exponential in case of ETRS"),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAT_EXP",&
            description="Which method should be used to calculate the exponential"//&
            " in the propagator. For Ehrenfest MD only the Taylor method works, "//&
            "for real time propagation diagonalization works as well.",&
            usage="MAT_EXP TAYLOR",default_i_val=do_arnoldi, &
            enum_c_vals=s2a("TAYLOR","DIAG","PADE","ARNOLDI"),&
            enum_i_vals=(/do_taylor,do_diag,do_pade,do_arnoldi/),&
            enum_desc=s2a("exponential is evaluated using scaling and squaring in combination"//&
            " with a taylor expansion of the exponential.",&
            "uses the eigenvalues and eigenvectors to calculate the exponential.",&
            "uses scaling and squaring together with the pade approximation",&
            "uses arnoldi subspace algorithm to compute exp(H)*MO directly, can't be used in "//&
            "combination with Crank Nicholson"),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SC_CHECK_START",&
            description="Speciefies how many iteration steps will be done without "//&
            "a check for self consistency. Can save some time in big calculations.",&
            usage="SC_CHECK_START 3",&
            default_i_val=0,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXP_ACCURACY",&
            description="Accuracy for the taylor and pade approximation. "//&
            "This is only an upper bound bound since the norm used for the guess "//&
            "is an upper bound for the needed one.",&
            usage="EXP_ACCURACY 1.0E-6",&
            default_r_val=1.0E-9_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PROPAGATOR",&
            description="Which propagator should be used for the orbitals",&
            usage="PROPAGATOR ETRS",default_i_val=do_etrs, &
            enum_c_vals=s2a("ETRS","CN","EM"),&
            enum_i_vals=(/do_etrs,do_cn,do_em/),&
            enum_desc=s2a("enforced time reversible symmetry",&
            "Crank Nicholson propagator",&
            "Exponential midpoint propagator"),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INITIAL_WFN",&
            description="Controls the initial WFN used for propagation.",&
            usage="INITIAL_WFN SCF_WFN",default_i_val=use_scf_wfn, &
            enum_c_vals=s2a("SCF_WFN","RESTART_WFN","RT_RESTART"),&
            enum_i_vals=(/use_scf_wfn,use_restart_wfn,use_rt_restart/),&
            enum_desc=s2a("An SCF run is performed to get the initial state.",&
            "A wavefunction from a previous scf is propageted. Especially useful,"//&
            " if electronic constraints or restraints are used in the previous calculation, "//&
            "since these do not work in the rtp scheme.",&
            "use the wavefunction of a real time propagation/ehrenfest run"),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="APPLY_DELTA_PULSE",&
            description="Applies a delta kick to the initial wfn (only RTP for now - the EMD "//&
            " case is not yet implemented).",&
            usage="APPLY_DELTA_PULSE",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PERIODIC",&
            description="Apply a delta-kick that is compatible with periodic boundary conditions"//&
            " for any value of DELTA_PULSE_SCALE. Uses perturbation theory for the preparation of"//&
            " the initial wfn. Note that the pulse is only applied when INITIAL_WFN is set to SCF_WFN,"//&
            " and not for restarts (RT_RESTART).",&
            usage="PERIODIC",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELTA_PULSE_DIRECTION",&
            description="Direction of the applied electric field. The k vector is given as"//&
            " 2*Pi*[i,j,k]*inv(h_mat), which for PERIODIC .FALSE. yields exp(ikr) periodic with"//&
            " the unit cell, only if DELTA_PULSE_SCALE is set to unity. For an orthorhombic cell"//&
            " [1,0,0] yields [2*Pi/L_x,0,0]. For small cells, this results in a very large kick.",&
            usage="DELTA_PULSE_DIRECTION 1 1 1",n_var=3,default_i_vals=(/1,0,0/),&
            type_of_var=integer_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELTA_PULSE_SCALE",&
            description="Scale the k vector, which for PERIODIC .FALSE. results in exp(ikr) no"//&
            " longer being periodic with the unit cell. The norm of k is the strength of the"//&
            " applied electric field in atomic units.",&
            usage="DELTA_PULSE_SCALE 0.01 ",n_var=1,default_r_val=0.001_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="HFX_BALANCE_IN_CORE",&
            description="If HFX is used, this keyword forces a redistribution/recalculation"//&
            " of the integrals, balanced with respect to the in core steps.",&
            usage="HFX_BALANCE_IN_CORE",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       NULLIFY(print_section)
       CALL section_create(print_section,name="PRINT",&
            description="Section of possible print options for an RTP runs",&
            repeats=.FALSE., required=.TRUE.,error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing within real time propagation and Eherenfest dynamics",&
            print_level=low_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(print_section,print_key,error=error)
       CALL section_release(print_key,error=error)


       CALL cp_print_key_section_create(print_key,"RESTART",&
            description="Controls the dumping of the MO restart file during rtp."//&
                        "By default keeps a short history of three restarts."//&
                        "See also RESTART_HISTORY", &
            print_level=low_print_level, common_iter_levels=3,&
            each_iter_names=s2a("MD"),each_iter_values=(/20/), &
            add_last=add_last_numeric,filename="RESTART",error=error)
       CALL keyword_create(keyword, name="BACKUP_COPIES",&
            description="Specifies the maximum index of backup copies.",&
            usage="BACKUP_COPIES {int}",&
            default_i_val=3, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(print_section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"RESTART_HISTORY",&
            description="Dumps unique MO restart files during the run keeping all of them.",&
            print_level=low_print_level, common_iter_levels=0,&
            each_iter_names=s2a("MD"),&
            each_iter_values=(/500/), &
            filename="RESTART",error=error)
       CALL keyword_create(keyword, name="BACKUP_COPIES",&
            description="Specifies the maximum index of backup copies.",&
            usage="BACKUP_COPIES {int}",&
            default_i_val=3, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(print_section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section,print_section,error=error)
       CALL section_release(print_section,error=error)


    END IF
  END SUBROUTINE create_rtp_section

! *****************************************************************************
  SUBROUTINE create_constant_env_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_constant_env_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"CONSTANT_ENV",&
            description="parameters for a constant envelop",&
            n_keywords=6, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="START_STEP",&
            description="First step the field is applied ",&
            usage="START_STEP 0",&
            default_i_val=0,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="END_STEP",&
            description="Last step the field is applied",&
            usage="END_STEP 2",&
            default_i_val=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF

  END SUBROUTINE create_constant_env_section

  SUBROUTINE create_gaussian_env_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_gaussian_env_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"GAUSSIAN_ENV",&
            description="parameters for a gaussian envelop",&
            n_keywords=6, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="T0",&
            description="Center of the gaussian envelop (maximum of the gaussian)",&
            usage="T0 2.0E0",&
            default_r_val=0.0E0_dp,&
            unit_str="fs",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SIGMA",&
            description="Width of the gaussian ",&
            usage="SIGMA 2.0E0",&
            default_r_val=-1.0E0_dp,&
            unit_str="fs",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF

  END SUBROUTINE create_gaussian_env_section
! *****************************************************************************
  SUBROUTINE create_ramp_env_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ramp_env_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"RAMP_ENV",&
            description="Parameters for an trapeziodal envelop ",&
            n_keywords=6, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="START_STEP_IN",&
            description="Step when the electric field starts to be applied ",&
            usage="START_STEP_IN 0",&
            default_i_val=0,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="END_STEP_IN",&
            description="Step when the field reaches the full strength",&
            usage="END_STEP_IN 2",&
            default_i_val=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="START_STEP_OUT",&
            description="Step when the field starts to vanish ",&
            usage="START_STEP 0",&
            default_i_val=0,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="END_STEP_OUT",&
            description="Step when the field disappears",&
            usage="END_TIME 2",&
            default_i_val=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF

  END SUBROUTINE create_ramp_env_section


END MODULE input_cp2k_dft
