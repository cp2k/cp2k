!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_cp2k_dft *
!!
!!   NAME
!!     input_cp2k_dft
!!
!!   FUNCTION
!!     function that build the dft section of the input
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     10.2005 moved out of input_cp2k [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE input_cp2k_dft
  USE bibliography,                    ONLY: &
       Avezac2005, Becke1988, Elstner1998, Goedecker1996, Heyd2004, &
       Iannuzzi2005, Iannuzzi2006, Iannuzzi2007, Krack2002, Lee1988, &
       Lippert1997, Lippert1999, Ortiz1994, Perdew1981, Perdew1996, &
       Porezag1995, Seifert1996, VandeVondele2003, VandeVondele2005a, &
       VandeVondele2005b, Vosko1980, Zhang1998, Zhechkov2005
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE input_constants
  USE input_cp2k_distribution,         ONLY: create_distribution_section
  USE input_cp2k_hfx,                  ONLY: create_hf_section
  USE input_cp2k_poisson,              ONLY: create_poisson_section
  USE input_cp2k_resp,                 ONLY: create_resp_section
  USE input_cp2k_rsgrid,               ONLY: create_rsgrid_section
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
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_dft'

  PUBLIC :: create_dft_section
  PUBLIC :: create_bsse_section, create_qs_section
  PUBLIC :: create_properties_section, create_scf_section, create_ep_section
  !***
  !****************************************************************************
CONTAINS

!!****f* input_cp2k_dft/create_ep_section *
!!
!!   NAME
!!     create_ep_section
!!
!!   FUNCTION
!!     creates an ep section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
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
         description="How many points are printed out (the bigger the smaller the cube file",&
         usage="stride 3",&
         default_i_val=2,error=error)
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
  !**************************************************************************

!!****f* input_cp2k_dft/create_dft_section *
!!
!!   NAME
!!     create_dft_section
!!
!!   FUNCTION
!!     creates the dft section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
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
         default_lc_val="BASIS_SET",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="POTENTIAL_FILE_NAME",&
         description="Name of the pseudo potential file, may include a path",&
         usage="POTENTIAL_FILE_NAME <FILENAME>",&
         default_lc_val="POTENTIAL",error=error)
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

    NULLIFY(subsection)
    CALL create_scf_section(subsection,error)
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

    CALL create_sic_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_field_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_poisson_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_implicit_solv_section(subsection, error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)       

    CALL create_density_fitting_section(subsection, error)
    CALL section_add_subsection(section,subsection,error=error)
    CALL section_release(subsection,error=error)       

    CALL create_xas_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_kg_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

    CALL create_print_dft_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
    CALL section_release(subsection,error=error)

  END SUBROUTINE create_dft_section
  !***************************************************************************

!!****f* input_cp2k_dft/create_implicit_solv_section *
!!
!!   NAME
!!     create_implicit_solv_section
!!
!!   FUNCTION
!!     Implicit Solvation Model
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     tlaino
!!
!!*** **********************************************************************
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

!!****f* input_cp2k_dft/create_sphere_section *
!!
!!   NAME
!!     create_sphere_section
!!
!!   FUNCTION
!!     Create Sphere cavity
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     tlaino
!!
!!*** **********************************************************************
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

!!****f* input_cp2k_dft/create_center_section *
!!
!!   NAME
!!     create_center_section
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     tlaino
!!
!!*** **********************************************************************
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

!****f* input_cp2k_dft/create_localize_section *
!!
!!   NAME
!!     create_localize_section
!!   FUNCTION
!!     parameters fo the localization of wavefunctions
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2005 created [MI]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE create_localize_section(section, error)

    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_localize_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       NULLIFY(keyword, print_key, subsection)
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
            enum_c_vals=s2a("NONE","JACOBI","CRAZY","DIRECT","L1SD"),&
            enum_i_vals=(/do_loc_none, do_loc_jacobi, do_loc_crazy, do_loc_direct,do_loc_l1_norm_sd/),&
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

       CALL keyword_create(keyword,name="OPERATOR",&
            description="Type of opertator which defines the spread functional",&
            usage="OPERATOR (BERRY|BOYS|PIPEK)",&
            enum_c_vals=s2a("BERRY","BOYS","PIPEK"),&
            enum_i_vals=(/op_loc_berry, op_loc_boys, op_loc_pipek/),&
            default_i_val=op_loc_berry,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       
       CALL keyword_create(keyword, name="RANGE",&
            description="The lower and upper index of the orbitals to be localized",&
            usage="RANGE {integer} {integer}",&
            n_var=2,default_i_vals=(/0,0/),type_of_var=integer_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword,name="LIST",&
            description="Indexes of the wfn to be localized"//&
            "This keyword can be repeated several times"//&
            "(useful if you have to specify many indexes).",&
            usage="LIST 1 2",&
            n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(subsection,"PRINT",&
            "printing of information about the localization of the orbitals",&
            error=error,repeats=.FALSE.,required=.FALSE.)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing basic info about the method", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"WANNIER_CUBES",&
            description="Controls the printing of the wannier functions ", &
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="the stride used to write the cube file "//&
            "(the bigger, the smaller the cube gets)",&
            default_i_val=2, error=error)
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
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"WANNIER_CENTERS",&
            description="Controls the printing of the wannier functions", &
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)

       CALL keyword_create(keyword, name="IONS+CENTERS",&
            description="prints out the wannier centers together with the particles",&
            usage="IONS+CENTERS", default_l_val=.FALSE., &
            lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="format",&
            description="Specifies the format of the output file when IONS+CENTERS is enabled.",&
            usage="FORMAT (ATOMIC|XMOL|DCD)",&
            default_i_val=dump_xmol,&
            enum_c_vals=output_name_types,&
            enum_i_vals=output_name_index,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"WANNIER_SPREADS",&
            description="Controls the printing of the wannier functions", &
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_localize_section

  !***************************************************************************
  !!****f* input_cp2k_dft/create_print_dft_section [1.0] *
  !!
  !!   NAME
  !!     create_print_dft_section
  !!
  !!   FUNCTION
  !!     Create the print dft section
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
  SUBROUTINE create_print_dft_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_print_dft_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, print_key2, &
                                                subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options in DFT code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(print_key, keyword, print_key2, subsection)

       CALL cp_print_key_section_create(print_key,"PROGRAM_BANNER",&
            description="Controls the printing of the banner of the MM program",&
            print_level=silent_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)       

       CALL cp_print_key_section_create(print_key,"derivatives",&
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
            print_level=debug_print_level,&
            filename="",&
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
            name="sac_ppnl",&
            description="Activates the printing of the orbital(a) "//&
            "GTH-PPNL neighbor lists (non-local part of the"//&
            "Goedecker-Teter-Hutter pseudo potentials)",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sbc_ppnl",&
            description="Activates the printing of the orbital(b) "//&
            "GTH-PPNL neighbor lists (non-local part of the"//&
            "Goedecker-Teter-Hutter pseudo potentials)",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sac_oce",&
            description="Activates the printing of the orbital "//&
            "PAW-projector neighbor lists (only GAPW)",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sbc_oce",&
            description="Activates the printing of the orbital "//&
            "PAW-projector neighbor lists (only GAPW)",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

       CALL keyword_create(keyword=keyword,&
            name="sab_2c",&
            description="Activates the printing of the two-center "//&
            "orbital V(Hartree_local) neighbor lists "//&
            "(only GAPW)",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

       CALL keyword_create(keyword=keyword,&
            name="sac_3c",&
            description="Activates the printing of the three-center "//&
            "orbital(a) V(Hartree_local) neighbor lists "//&
            "(only GAPW)",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="sbc_3c",&
            description="Activates the printing of the three-center "//&
            "orbital(b) V(Hartree_local) neighbor lists "//&
            "(only GAPW)",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SUBCELL",&
            description="Activates the printing of the subcells used for the"//&
            "generation of neighbor lists.", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"ao_matrices",&
            description="Controls the printing of the ao (i.e. contracted gaussian) matrices (debug).", &
            print_level=debug_print_level,filename="__STD_OUT__",&
            error=error)
       CALL keyword_create(keyword, name="core_hamiltonian",&
            description="If the printkey is activated controls the printing of the hamiltonian matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="density",&
            description="If the printkey is activated controls the printing of the density (P) matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

       CALL keyword_create(keyword, name="kinetic_energy",&
            description="If the printkey is activated controls the printing of the kinetic energy matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL keyword_create(keyword, name="kohn_sham_matrix",&
            description="If the printkey is activated controls the printing of the kohn-sham matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL keyword_create(keyword, name="ortho",&
            description="If the printkey is activated controls the printing of the orthogonalization matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL keyword_create(keyword, name="overlap",&
            description="If the printkey is activated controls the printing of the overlap matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL keyword_create(keyword, name="potential_energy",&
            description="If the printkey is activated controls the printing of the potential energy matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="oce_hard",&
            description="If the printkey is activated controls the printing of the OCE HARD matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="oce_soft",&
            description="If the printkey is activated controls the printing of the OCE SOFT matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="w_matrix",&
            description="If the printkey is activated controls the printing of the w matrix",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="derivatives",&
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
       CALL keyword_create(keyword, name="cartesian",&
            description="If the printkey is activated controls the printing of the mo in the cartesian basis",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="eigenvalues",&
            description="If the printkey is activated controls the printing of the eigenvalues of the mos",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL keyword_create(keyword, name="eigenvectors",&
            description="If the printkey is activated controls the printing of the eigenvectors of the mos",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL keyword_create(keyword, name="occupation_numbers",&
            description="If the printkey is activated controls the printing of the occupation numbers of the mos",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MO_CUBES",&
            description="Controls the printing of cubes of the molecular orbitals (squared).", &
            print_level=high_print_level,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="the stride used to write the cube file "//&
            "(the bigger, the smaller the cube gets)",&
            default_i_val=2, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL keyword_create(keyword, name="write_cube",&
            description="If the cube file should be written. If false, the eigenvalues are still computed.",&
            default_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL keyword_create(keyword, name="nlumo",&
            description="If the printkey is activated controls the number of lumos"//&
                        "that are printed and dumped as a cube (-1=all)",&
            default_i_val=5, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL keyword_create(keyword, name="nhomo",&
            description="If the printkey is activated controls the number of homos that dumped as a cube (-1=all),"//&
                        " eigenvalues are always all dumped",&
            default_i_val=1, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)


       CALL cp_print_key_section_create(print_key,"gapw",&
            description="Controls the printing of some gapw related information (debug).", &
            print_level=debug_print_level,filename="__STD_OUT__",&
            error=error)
       CALL keyword_create(keyword, name="projectors",&
            description="If the printkey is activated controls if information on"//&
            " the projectors is printed.",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="rho0_information",&
            description="If the printkey is activated controls if information on rho0 is printed.",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

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

       CALL cp_print_key_section_create(print_key,"e_density_cube",&
            description="Controls the printing of cube files with "//&
            "the total electronic density and, for LSD "//&
            "calculations, the spin density", &
            print_level=high_print_level,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="the stride used to write the cube file "//&
            "(the bigger, the smaller the cube gets)",&
            default_i_val=2, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"tot_density_cube",&
            description="Controls the printing of cube files with "//&
            "the total density (electronic+core)", &
            print_level=high_print_level,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="the stride used to write the cube file "//&
            "(the bigger, the smaller the cube gets)",&
            default_i_val=2, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"v_hartree_cube",&
            description="Controls the printing of a cube file with eletrostatic potential (electrons+core).", &
            print_level=high_print_level,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="the stride used to write the cube file "//&
            "(the bigger, the smaller the cube gets)",&
            default_i_val=2, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"LOCALIZATION",&
            Description="Controls the printing of the Wannier centers.", &
            print_level=high_print_level,filename="__STD_OUT__", error=error)

       CALL create_localize_section(subsection,error)
       CALL section_add_subsection(print_key, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL keyword_create(keyword, name="TOTAL_DIPOLE",&
            description="Print out the total dipole calcuated from the Wannier centers. "//&
            "This only makes sense for isolated systems.",&
            usage="TOTAL_DIPOLE", default_l_val=.FALSE., &
            lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOLECULAR_DIPOLES",&
            description="Print out molecular dipoles calculated from Wannier centers",&
            usage="MOLECULAR_DIPOLES", default_l_val=.TRUE., &
            lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       !Printing of Moments
       CALL cp_print_key_section_create(print_key,"MOMENTS",&
            description="Controls printing of electric moments", &
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)

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
            name="PERIODIC",&
            description="Use Berry phase formula (PERIODIC=T) or simple operator (PERIODIC=F)",&
            usage="PERIODIC {logical}",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE.,& 
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword, name="REFERENCE",&
            variants=s2a("REF"),&
            description="Define the reference point for electrical moment calculations.",&
            usage="REFERENCE COM",&
            enum_c_vals=s2a("COM","COAC","USER_DEFINED","ZERO","DEFAULT"),&
            enum_desc=s2a("Use Center of Mass",&
                          "Use Center of Atomic Charges",&
                          "Use User Defined Point (Keyword:REF_POINT)",&
                          "Use Origin of Coordinate System",&
                          "Use Default (depends on method: periodic or localized)"),&
            enum_i_vals=(/ use_mom_ref_com,&
                           use_mom_ref_coac,&
                           use_mom_ref_user,&
                           use_mom_ref_zero,&
                           use_mom_ref_nondef /),&
            default_i_val=use_mom_ref_nondef, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword, name="REFERENCE_POINT",&
            variants=s2a("REF_POINT"),&
            description="Fixed reference point for electrical moment calculations.",&
            usage="REFERENCE_POINT x y z", &
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=3,default_r_vals=(/0._dp,0._dp,0._dp/),&
            type_of_var=real_t, &
            unit_str='bohr',&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
       !End of Printing of Moments

       CALL cp_print_key_section_create(print_key,"MOLECULAR_STATES",&
            description="Controls printing of molecular states (requires LOCALIZATION keyword)", &
            print_level=high_print_level,filename=" ",&
            error=error)

       CALL cp_print_key_section_create(print_key2,"cubes",&
            description="Controls the printing of cube files in localization", &
            print_level=high_print_level,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="the stride used to write the cube file "//&
            "(the bigger, the smaller the cube gets)",&
            default_i_val=2, error=error)
       CALL section_add_keyword(print_key2,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL section_add_subsection(print_key,print_key2,error=error)
       CALL section_release(print_key2,error=error)       

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MULLIKEN",&
            description="Controls the printing of the Mulliken (spin) population analysis.", &
            print_level=medium_print_level,filename="__STD_OUT__",&
            common_iter_levels=1, error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       ! fairly expensive to compute, so only at high
       CALL cp_print_key_section_create(print_key,"LOWDIN",&
            description="Controls the printing of the Lowdin charges.", &
            print_level=high_print_level,filename="__STD_OUT__",&
            common_iter_levels=1, error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(&
            print_key_section=print_key,&
            name="XRAY_DIFFRACTION_SPECTRUM",&
            description="Calculate and print the coherent X-ray "//&
            "diffraction spectrum",&
            print_level=debug_print_level,&
            filename="",&
            citations=(/Krack2002/),&
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

    END IF

  END SUBROUTINE create_print_dft_section

  !!****f* input_cp2k_dft/create_qs_section [1.0] *
  !!
  !!   NAME
  !!     create_qs_section
  !!
  !!   FUNCTION
  !!     creates the input section for the qs part
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
            n_keywords=34, n_subsections=2, repeats=.FALSE., required=.FALSE.,&
            citations=(/Lippert1997,Lippert1999,Krack2002,VandeVondele2005a/),&
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

       CALL keyword_create(keyword, name="EPS_3C_REDUCE",&
            description="GAPW: Factor to reduce the precision in the construction"//&
            " of the 3 center lists for the calculation of the OCE coefficients.",&
            usage="EPS_3C_REDUCE real", default_r_val=1.0_dp,&
            type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_PPL",&
            description="Sets precision of the local part of the pseudo potential.Overrides EPS_DEFAULT/100.0 value",&
            usage="EPS_PPL real", type_of_var=real_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_PPNL",&
            description="Sets precision of the non-local part of the pseudo potential.Overrides EPS_DEFAULT/100.0 value",&
            usage="EPS_PPNL real", type_of_var=real_t,&
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

       CALL keyword_create(keyword, name="EPS_KG_GPW",&
            description="Sets precision of KG_GPW single molecule density integration."//&
            " Overrides SQRT(EPS_DEFAULT)*1000.0 value",&
            usage="EPS_KG_GPW real", type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPSFIT",&
            variants=(/"EPS_FIT"/),&
            description="GAPW : precision to give the extention of a hard gaussian ",&
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
            variants=(/"ALPHA0_H"/),&
            description="GAPW : Exponent for hard compensation charge",&
            usage="ALPHA0_HARD real", default_r_val=0.0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALPHA0_SOFT",&
            variants=(/"ALPHA0_S"/),&
            description="GAPW : Exponent for soft compensation charge",&
            usage="ALPHA0_SOFT real", default_r_val=0.0_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FULL_GAPW",&
            description="Perform a full GAPW calculations, with all the terms"//&
            " if alpha0_h/=alpha0_s (much more expensive). Otherwise the "//&
            "compensation densities are set equal and the 2/3 center terms cancel out",&
            usage="FULL_GAPW",&
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
       CALL keyword_create(keyword, name="MAP_PAA",&
            description="Put diagonal density matrix elements on the highest cutoff multigrid",&
            usage="MAP_PAA",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
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

       NULLIFY(subsection)

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
            description="GAPW : XXXXXXXXXX",&
            usage="LADDN0 integer",&
            default_i_val=99,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Characters
       CALL keyword_create(keyword, name="QUADRATURE",&
            description="GAPW: algorithm to construct the atomic radial grids",&
            usage="QUADRATURE GC_SIMPLE",&
            enum_c_vals=s2a("GC_SIMPLE","GC_TRANSFORMED","GC_LOG"),&
            enum_i_vals=(/ do_gapw_gcs,do_gapw_gct,do_gapw_log/),&
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

       CALL keyword_create(keyword, name="EXTRAPOLATION",&
            variants=s2a("INTERPOLATION","WF_INTERPOLATION"),&
            description="Extrapolation strategy for the wavefunction during e.g. MD."//&
                        "PS and ASPC are recommended, see also EXTRAPOLATION_ORDER.",&
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
            enum_i_vals=(/ do_extrpl_use_guess,& 
            do_extrpl_use_prev_p,&
            do_extrpl_use_prev_rho_r,&
            do_extrpl_linear_wf,&
            do_extrpl_linear_p,&
            do_extrpl_linear_ps,&
            do_extrpl_use_prev_wf,&
            do_extrpl_ps,&
            do_extrpl_frozen,&
            do_extrpl_aspc /),&
            default_i_val=do_extrpl_ps, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXTRAPOLATION_ORDER",&
            description="Order for the PS extrapolation (typically 2-4). "//&
                        "higher order might bring more accuracy, but comes, "//&
                        "for large systems, also at some cost. "//&
                        "In some cases, a high order extrapolation is not stable,"//&
                        " and the order needs to be reduced.",&
            usage="EXTRAPOLATION_ORDER {integer}",default_i_val=3, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="METHOD",&
            description="Specifies the electronic structure method that should be employed",&
            usage="METHOD GAPW",&
            enum_c_vals=s2a("GAPW","GAPW_XC","GPW","KG_POL","KG_NOPOL","KG_GPW",&
            "MNDO","AM1","PM3","PDG","RM1","DFTB"),&
            enum_desc=s2a("Gaussian and augmented plane waves method",&
                          "Gaussian and augmented plane waves method only for XC",&
                          "Gaussian and plane waves method",&
                          "Polarisable Kim-Gordon method",&
                          "Non-polarisable Kim-Gordon method", " DFT Embedding method ",&
                          "MNDO semiempirical","AM1 semiempirical",&
                          "PM3 semiempirical","PDG semiempirical","RM1 semiempirical",&
                          "DFTB: Density Functional based Tight-Binding"),&
            enum_i_vals=(/ do_method_gapw, do_method_gapw_xc, do_method_gpw, &
            do_method_kg_pol, do_method_kg_nopol, do_method_kg_gpw, &
            do_method_mndo, do_method_am1, do_method_pm3, do_method_pdg, &
            do_method_rm1,do_method_dftb/),&
            citations=(/Lippert1997,Lippert1999,VandeVondele2005a,Iannuzzi2006/),&
            default_i_val=do_method_gpw, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsection)
       CALL create_distribution_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_dftb_control_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_se_control_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_mulliken_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_ddapc_restraint_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_becke_constraint_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_s2_restraint_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_harris_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error) 

    END IF

  END SUBROUTINE create_qs_section

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
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            citations=(/Porezag1995, Seifert1996, Elstner1998, Zhechkov2005/),&
            error=error)

       NULLIFY(subsection)
       CALL create_dftb_parameter_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="self_consistent",&
            description="Use self_consistent method",&
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
            usage="DO_EWALD",default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="dispersion",&
            description="Use dispersion correction",&
            citations=(/Zhechkov2005/),&
            usage="DISPERSION",default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="eps_disp",&
            description="Define accuracy of dispersion interaction",&
            usage="EPS_DISP",default_r_val=0.0001_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_dftb_control_section

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
            description="Specify the directory with the DFTB parameter files",&
            usage="PARAM_FILE_PATH pathname",&
            n_var=1,type_of_var=char_t,default_c_val="./",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="PARAM_FILE_NAME",&
            description="Specify file that contains the SK_FILE names",&
            usage="PARAM_FILE_NAME filename",&
            n_var=1,type_of_var=char_t,default_c_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="UFF_FORCE_FIELD",&
            description="Name of file with UFF parameters",&
            usage="UFF_FORCE_FIELD filename",required=.FALSE.,&
            n_var=1,type_of_var=char_t,default_c_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_dftb_parameter_section

  SUBROUTINE create_se_control_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_se_control_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"SE",&
            description="Parameters needed to set up the Semi-empirical methods",&
            n_keywords=7, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="orthogonal_basis",&
            description="Assume orthogonal basis set. This flag is overwritten by "//&
                        "methods with fixed orthogonal/non-orthogonal basis set.",&
            usage="ORTHOGONAL_BASIS",default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="analytical_gradients",&
            description="Nuclear Gradients are computed analytically or numerically",&
            usage="ANALYTICAL_GRADIENTS",default_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DELTA",&
            description="Step size in finite difference force calculation",&
            usage="DELTA {real} ",default_r_val=1.e-6_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RC_INTERACTION",&
            description="Cutoff radius for localization",&
            usage="RC_INTERACTION {real} ",default_r_val=24._dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RC_COULOMB",&
            description="Cutoff radius for Coulomb interactions",&
            usage="RC_COULOMB {real} ",default_r_val=24._dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RC_RANGE",&
            description="Range of cutoff switch function (tapering)",&
            usage="RC_RANGE {real} ",default_r_val=8._dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="do_ewald",&
            description="Use an Ewald type method for the calculation of the "//&
                        "Coulomb energy",&
            usage="DO_EWALD",default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_se_control_section

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
            description="Setup parameters for density fitting (Bloechl charges)",&
            n_keywords=7, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
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

  SUBROUTINE create_becke_constraint_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_becke_constraint_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    NULLIFY(keyword, print_key)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"BECKE_CONSTRAINT",&
            description="Use Becke weight population in a constraint ",&
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
            description="Defines the the coefficient of the atom in the atom list (default is one), currently DDAPC only ",&
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

  END SUBROUTINE create_becke_constraint_section

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

  !!****f* xc_derivatives/create_xc_fun_section [1.0] *
  !!   NAME
  !!     create_xc_fun_section
  !!
  !!   FUNCTION
  !!     creates the structure of the section needed to select the xc functional
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - section: the section that will be created
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE create_xc_fun_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_xc_fun_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="xc_functional",&
            description="The xc functional to use",&
            n_keywords=0, n_subsections=4, repeats=.FALSE., required=.FALSE.,&
            citations=(/Ortiz1994,Becke1988,Perdew1996,Zhang1998,Lee1988, &
                        Heyd2004,Vosko1980, Goedecker1996,Perdew1981/),&
            error=error)

       NULLIFY(subsection,keyword)
       CALL keyword_create(keyword,name="_SECTION_PARAMETERS_",&
            description="shorcut for the most commont functionals combinations",&
            usage="&xc_functional BLYP",&
            enum_c_vals=s2a("B3LYP","PBE0","BLYP","BP","PADE","PBE","TPSS","HCTH120","OLYP","NO_SHORTCUT","NONE"),&
            enum_i_vals=(/xc_funct_b3lyp,xc_funct_pbe0,xc_funct_blyp,xc_funct_bp,xc_funct_pade,xc_funct_pbe,&
            xc_funct_tpss,xc_funct_hcth120,xc_funct_olyp,xc_funct_no_shortcut,xc_none/),&
            default_i_val=xc_funct_no_shortcut,lone_keyword_i_val=xc_funct_no_shortcut,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(subsection,name="BECKE88",&
            description="Uses the Becke 88 exchange functional",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            citations=(/Becke1988/),&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"scale_x",&
            description="scales the exchange part of the functional",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="LYP_ADIABATIC",&
            description="Uses the LYP correlation functional in an adiabatic fashion",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            citations=(/Lee1988/),&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"LAMBDA",&
            description="Defines the parameter of the adiabatic curve.",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="BECKE88_LR",&
            description="Uses the Becke 88 longrange exchange functional",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            citations=(/Becke1988/),&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"scale_x",&
            description="scales the exchange part of the functional",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"OMEGA",&
            description="Potential parameter in erf(omega*r)/r",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="LYP",&
            description="Uses the LYP functional",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            citations=(/Lee1988/),&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"scale_c",&
            description="scales the correlation part of the functional",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="PADE",&
            description="Uses the PADE functional",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            citations=(/Goedecker1996/),&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="HCTH",&
            description="Uses the HCTH class of functionals",&
            n_keywords=0, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword,name="PARAMETER_SET",&
            description="Which version of the parameters should be used",&
            usage="PARAMETER_SET 407",&
            enum_c_vals=(/"93 ","120","147","407"/),&
            enum_i_vals=(/93,120,147,407/),&
            default_i_val=120,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="OPTX",&
            description="Uses the OPTX functional",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"scale_x",&
            description="scales the exchange part of the functional",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="CS1",&
            description="Uses the CS1 functional",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="XGGA",&
            description="Uses one of the XGGA functionals (optimized versions of "//&
            "some of these functionals might be available outside this section).",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,name="FUNCTIONAL",&
            description="Which one of the XGGA functionals should be used",&
            usage="FUNCTIONAL PW86X",&
            enum_c_vals=(/&
            "BECKE88X",&
            "PW86X   ",&
            "PW91X   ",&
            "PBEX    ",&
            "REV_PBEX",&
            "OPTX    "/),&
            enum_i_vals=(/xgga_b88x,xgga_pw86,xgga_pw91,xgga_pbex,xgga_revpbe,xgga_opt/),&
            default_i_val=xgga_b88x,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="KE_GGA",&
            description="Uses one of the KE_GGA functionals (optimized versions of "//&
            "some of these functionals might be available outside this section).",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,name="FUNCTIONAL",&
            description="Which one of the KE_GGA functionals should be used",&
            usage="FUNCTIONAL LLP",&
            enum_c_vals=(/"OL1 ","OL2 ","LLP ","PW86","PW91","LC  ","T92 ","PBE "/),&
            enum_i_vals=(/ke_ol1,ke_ol2,ke_llp,ke_pw86,ke_pw91,ke_lc,ke_t92,ke_pbe/),&
            default_i_val=ke_ol1,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="P86C",&
            description="Uses the P86C functional",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"scale_c",&
            description="scales the correlation part of the functional",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="PW92",&
            description="Uses the PerdewWang correlation functional.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"SCALE",&
            description="Scaling of the energy functional",&
            default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,name="PARAMETRIZATION",&
            description="Which one of parametrizations should be used",&
            usage="PARAMETRIZATION DMC",&
            enum_c_vals=(/&
            "ORIGINAL",&
            "DMC     ",&
            "VMC     "/),&
            enum_i_vals=(/c_pw92,c_pw92dmc,c_pw92vmc/),&
            default_i_val=c_pw92,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="PZ81",&
            description="Uses the PZ functional.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            citations=(/Perdew1981,Ortiz1994/), error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,name="PARAMETRIZATION",&
            description="Which one of parametrizations should be used",&
            usage="PARAMETRIZATION DMC",&
            enum_c_vals=(/&
            "ORIGINAL",&
            "DMC     ",&
            "VMC     "/),&
            enum_i_vals=(/c_pz,c_pzdmc,c_pzvmc/),&
            default_i_val=pz_orig,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"scale_c",&
            description="scales the correlation part of the functional",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="TFW",&
            description="Uses the TFW functional",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="TF",&
            description="Uses the TF functional",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="VWN",&
            description="Uses the VWN functional",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            citations=(/Vosko1980/),&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"scale_c",&
            description="scales the correlation part of the functional",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FUNCTIONAL_TYPE",&
         description="Which version of the VWN functional should be used",&
         usage="FUNCTIONAL_TYPE VWN5",&
         enum_c_vals=s2a("VWN5","VWN3"),&
         enum_i_vals=(/do_vwn5, do_vwn3 /),&
         enum_desc=s2a("This is the recommended (correct) version of the VWN functional", &
                       "This version is the default in Gaussian, but not recommended."//&
                       "Notice that it is also employed in Gaussian's default version of B3LYP"),&
         default_i_val=do_vwn5, error=error)

       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="XALPHA",&
            description="Uses the XALPHA (SLATER) functional.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,name="XA",&
            description="Value of the xa parameter (this does not change the exponent, "//&
            "just the mixing)",&
            usage="XA 0.7", default_r_val=2._dp/3._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"scale_x",&
            description="scales the exchange part of the functional",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="TPSS",&
            description="Uses the TPSS functional",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"scale_x",&
            description="scales the exchange part of the functional",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"scale_c",&
            description="scales the correlation part of the functional",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,name="PBE",&
            description="Uses the PBE functional",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            citations=(/Perdew1996,Zhang1998/),&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"parametrization",&
            description="switches between the original PBE and the revised PBE "//&
            "proposed ba Yingkai Zhang and Weitao Yang",&
            enum_i_vals=(/xc_pbe_orig,xc_rev_pbe/),enum_c_vals=(/"ORIG  ","revPBE"/),&
            default_i_val=xc_pbe_orig,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"scale_x",&
            description="scales the exchange part of the functional",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"scale_c",&
            description="scales the correlation part of the functional",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
      
       CALL section_create(subsection,name="XWPBE",&
            description="Uses the short range PBE functional",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            citations=(/Heyd2004/),&
            error=error)
       CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
            description="activates the functional",&
            lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"scale_x",&
            description="scales the exchange part of the functional",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"scale_x0",&
            description="scales the exchange part of the original hole PBE-functional",&
            default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,"omega",&
            description="screening parameter",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_xc_fun_section

  !!****f* xc_derivatives/create_xc_potential_section [1.0] *
  !!   NAME
  !!     create_xc_potential_section
  !!
  !!   FUNCTION
  !!     creates the structure of the section needed to select an xc potential
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - section: the section that will be created
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     thomas chassaing
  !!
  !!*** **********************************************************************
  SUBROUTINE create_xc_potential_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_xc_potential_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="xc_potential",&
            description="The xc potential to use (CAREFUL: xc potential here refers "//&
            "to potentials that are not derived from an xc functional, but rather are "//&
            "modelled directly. Therefore there is no consistent xc energy available. "//&
            "To still get an energy expression, see ENERGY below",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(subsection,keyword)
       CALL section_create(subsection,name="SAOP",&
            description="Uses the SAOP potential",&
            n_keywords=3, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
       CALL keyword_create(keyword,name="ALPHA",&
            description="Value of the alpha parameter (default = 1.19).",&
            usage="ALPHA 1.19", default_r_val=1.19_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,name="BETA",&
            description="Value of the beta parameter (default = 0.01).",&
            usage="BETA 0.01", default_r_val=0.01_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword,name="K_RHO",&
            description="Value of the K_rho parameter (default = 0.42).",&
            usage="ALPHA 0.42", default_r_val=0.42_dp,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL keyword_create(keyword, name="ENERGY",&
            description="How to determine the total energy.",&
            usage="ENERGY [NONE,XC_FUNCTIONAL,SUM_EIGENVALUES",&
            enum_c_vals=s2a("NONE","XC_FUNCTIONAL","FUNCTIONAL","SUM_EIGENVALUES",&
            "SOE"),&
            enum_i_vals=(/ xc_pot_energy_none, &
            xc_pot_energy_xc_functional,&
            xc_pot_energy_xc_functional,&
            xc_pot_energy_sum_eigenvalues,&
            xc_pot_energy_sum_eigenvalues /),&
            default_i_val=xc_pot_energy_none, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_xc_potential_section

  !***************************************************************************
  !!****f* input_cp2k_dft/create_tddfpt_section [1.0] *
  !!
  !!   NAME
  !!     create_tddfpt_section
  !!
  !!   FUNCTION
  !!     creates the input section for the tddfpt part
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
            n_keywords=5, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            citations=(/Iannuzzi2005/),&
            error=error)

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

  !!****f* input_cp2k_dft/create_xc_section [1.0] *
  !!
  !!   NAME
  !!     create_xc_section
  !!
  !!   FUNCTION
  !!     creates the input section for the xc part
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
  SUBROUTINE create_xc_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_xc_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"xc",&
            description="parameters needed calculate the xc potential",&
            n_keywords=5, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection,keyword)

       CALL keyword_create(keyword, name="density_cutoff",&
            description="The cutoff on the density used by the xc calculation",&
            usage="density_cutoff 1.e-11", default_r_val=1.0e-10_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="gradient_cutoff",&
            description="The cutoff on the gradient of the density used by the "//&
            "xc calculation",&
            usage="gradient_cutoff 1.e-11", default_r_val=1.0e-10_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DENSITY_SMOOTH_CUTOFF_RANGE",&
            description="Parameter for the smoothing procedure in"//&
            "xc calculation",&
            usage="gradient_cutoff {real}", default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="tau_cutoff",&
            description="The cutoff on tau used by the xc calculation",&
            usage="tau_cutoff 1.e-11", default_r_val=1.0e-10_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FUNCTIONAL_ROUTINE",&
            description="Select the code for xc calculation",&
            usage="FUNCTIONAL_ROUTINE NEW", default_i_val=xc_new_f_routine,&
            enum_c_vals=s2a("NEW","OLD","TEST_LSD","DEBUG"),&
            enum_i_vals=(/ xc_new_f_routine, xc_old_f_routine, &
            xc_test_lsd_f_routine, xc_debug_new_routine/),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(subsection, name="xc_grid",&!FM to do
            description="The xc parameters used when calculating the xc on the grid",&
            n_keywords=5, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)
       CALL keyword_create(keyword, name="xc_smooth_rho",&
            description="The density smoothing used for the xc calculation",&
            usage="xc_smooth_rho nn10", default_i_val=xc_rho_no_smooth,&
            enum_c_vals=s2a("NONE","NN50","NN10","SPLINE2","NN6","SPLINE3","NN4"),&
            enum_i_vals=(/ xc_rho_no_smooth, xc_rho_nn50, xc_rho_nn10, &
            xc_rho_spline2_smooth,xc_rho_spline2_smooth,xc_rho_spline3_smooth,&
            xc_rho_spline3_smooth/),&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="xc_deriv",&
            description="The method used to compute the derivatives",&
            usage="xc_deriv NN10_SMOOTH", default_i_val=xc_deriv_pw,&
            enum_c_vals=s2a("PW","SPLINE3","SPLINE2","NN50_SMOOTH","NN10_SMOOTH",&
            "SPLINE2_SMOOTH", "NN6_SMOOTH", "SPLINE3_SMOOTH", "NN4_SMOOTH"),&
            enum_i_vals=(/xc_deriv_pw, xc_deriv_spline3, xc_deriv_spline2,&
            xc_deriv_nn50_smooth, xc_deriv_nn10_smooth, xc_deriv_spline2_smooth,&
            xc_deriv_spline2_smooth, xc_deriv_spline3_smooth, xc_deriv_spline3_smooth/),&
            error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="use_finer_grid",&
            description="Uses a finer grid only to calculate the xc",&
            usage="use_finer_grid", default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_xc_fun_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_hf_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL create_xc_potential_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_xc_section
  !***************************************************************************

  !!****f* input_cp2k_dft/create_scf_section *
  !!
  !!   NAME
  !!     create_scf_section
  !!
  !!   FUNCTION
  !!     creates the structure of the section with the dft scf parameters
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - section: will contain the scf section
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
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
            n_keywords=25, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)
       CALL create_ot_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_outer_scf_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="MAX_ITER_LUMO",&
            variants=(/"MAX_ITER_LUMOS"/),&
            description="The maximum number of iteration for the lumo computation",&
            usage="MAX_ITER_LUMO 100", default_i_val=2999,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_LUMO",&
            variants=(/"EPS_LUMOS"/),&
            description="target accuracy of the computation of the lumo energy",&
            usage="EPS_LUMO 1.e-6", default_r_val=1.0e-5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WORK_SYEVX",&
            description="limit the memory usage of SYEVX to a fraction of the maximum size",&
            usage="work_syevx 0.8", default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SMEAR",&
            description="Use a finite temperature like smearing of the occupation numbers",&
            usage="SMEAR 0.1", default_r_val=0._dp,error=error)
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
            "the inverse of S",&
            usage="CHOLESKY OFF", default_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="JACOBI_THRESHOLD",&
            description="XXXXXXXX switch to jacobi if ...",&
            usage="JACOBI_THRESHOLD 1.0e-6", default_r_val=1.0E-7_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_JACOBI",&
            description="XXXXXXXX target jacobi precision",&
            usage="EPS_JACOBI 1.e-5", default_r_val=0._dp,error=error)
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

       CALL keyword_create(keyword, name="MIXING",&
            variants=(/"DENSITY_MIXING"/),&
            description="Fraction of new density to be mixed in",&
            usage="MIXING 0.1", default_r_val=0.4_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(&
            keyword=keyword,&
            name="GSPACE_MIXING",&
            variants=(/"GSPACE_DENSITY_MIXING"/),&
            description="Controls the mixing of the old and the new density "//&
            "in the reciprocal space combined with the Kerker "//&
            "damping in order to suppress charge sloshing:"//&
            "rho_mix(g) = rho_old(g) + A*g^2/(g^2 + B^2)*(rho_new(g) - rho_old(g))",&
            usage="GSPACE_MIXING {real A} {real B}",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=2,&
            type_of_var=real_t,&
            default_r_vals=(/1.0_dp,0.0_dp/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SCF_GUESS",&
            description="Change the initial guess for the wavefunction.",&
            usage="SCF_GUESS RESTART", default_i_val=atomic_guess,&
            enum_c_vals=s2a("ATOMIC","RESTART","RANDOM","CORE","DENSITIES","HISTORY_RESTART","NONE"),&
            enum_desc=s2a("Generate a trial atomic density based on the SZV part of the basis", &
            "Use the RESTART file as an initial guess (and ATOMIC if not present).", &
            "Use random wavefunction coefficients.", &
            "Diagonalize the core hamiltonian for an initial guess.", &
            "Use the aux_basis_set for collocation.", &
            "Extrapolated from previous RESTART files.", &
            "Skip initial guess (only for NON-SCC DFTB)."), &
            enum_i_vals=(/atomic_guess,restart_guess,random_guess,core_guess,&
            densities_guess,history_guess,no_guess/),&
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
            print_level=silent_print_level, common_iter_levels=3,&
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
            print_level=silent_print_level, common_iter_levels=0,&
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

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_scf_section
  !***************************************************************************

  !!****f* input_cp2k_dft/create_outer_scf_section *
  !!
  !!   NAME
  !!     create_outer_scf_section
  !!
  !!   FUNCTION
  !!     creates the structure of the section with scf parameters
  !!     controlling an other loop
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - section: will contain the scf section
  !!     - error: variable to control error logging, stopping,...
  !!       see module cp_error_handling
  !!
  !!   AUTHOR
  !!     Joost VandeVondele [2006.03]
  !!
  !!*** **********************************************************************
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
             "NONE" ),&
            enum_desc=s2a("Enforce a constraint on the DDAPC, requires the corresponding section", &
            "Enforce a constraint on the S2, requires the corresponding section", &
            "Enforce a constraint on the Becke weight population,requires the corresponding section", &
            "Do nothing in the outer loop, useful for resetting the inner loop,"),&
            enum_i_vals=(/outer_scf_ddapc_constraint,outer_scf_s2_constraint,&
            outer_scf_becke, outer_scf_none/),&
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

  !!****f* input_cp2k_dft/create_properties_section [1.0] *
  !!
  !!   NAME
  !!     create_properties_section
  !!
  !!   FUNCTION
  !!     Create the PROPERTIES section 
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
  SUBROUTINE create_properties_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_properties_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="PROPERTIES",&
            description="This section is used to set up the PROPERTIES calculation.",&
            n_keywords=3, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY (subsection)

       CALL create_linres_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
       
       CALL create_et_coupling_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
       
       CALL create_resp_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL cp_print_key_section_create(subsection,name="FIT_CHARGE",&
            description="This section is used to print the density derived atomic point charges."//&
            "The fit of the charges is controlled through the DENSITY_FITTING section",&
            print_level=high_print_level,filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_properties_section

  !***************************************************************************
  !!****f* input_cp2k_dft/create_bsse_section [1.0] *
  !!
  !!   NAME
  !!     create_bsse_section
  !!
  !!   FUNCTION
  !!     Create the BSSE section for counterpoise correction
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

       CALL keyword_create(keyword, name="RANGE",&
            description="Specifies the starting and the ending index of atoms range.",&
            usage="RANGE {integer} {integer}", required=.FALSE.,&
            n_var=2, type_of_var=integer_t, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! CONFIGURATION SECTION
       CALL section_create(subsection,name="CONFIGURATION",&
            description="Specify additional parameters for the combinatorial configurations.",&
            n_keywords=2, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
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


       CALL create_print_bsse_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_bsse_section

  !***************************************************************************
  !!****f* input_cp2k_dft/create_print_bsse_section [1.0] *
  !!
  !!   NAME
  !!     create_print_bsse_section
  !!
  !!   FUNCTION
  !!     Create the print bsse section
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
            print_level=silent_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_print_bsse_section

  !!****f* input_cp2k_qmmm/create_gspace_interp_section *
  !!
  !!   NAME
  !!     create_gspace_interp_section
  !!
  !!   FUNCTION
  !!     creates the interpolation section for the periodic QM/MM 
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     tlaino
  !!
  !!*** **********************************************************************
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
  !***************************************************************************

  !!****f* input_cp2k_dft/create_harris_section [1.0] *
  !!
  !!   NAME
  !!     create_harris_section
  !!
  !!   FUNCTION
  !!     This section specifies the flags for the harris functional
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
  !!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
  !!
  !!*** **********************************************************************
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

  !***************************************************************************

  !!****f* input_cp2k_dft/create_mgrid_section *
  !!
  !!   NAME
  !!     create_mgrid_section
  !!
  !!   FUNCTION
  !!     creates the multigrid
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
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

       CALL keyword_create(&
            keyword=keyword,&
            name="cutoff",&
            description="The cutoff (from the highest grid to the lowest)",&
            usage="cutoff 300",&
            default_r_val=cp_unit_to_cp2k(value=280.0_dp,&
            unit_str="Ry",&
            error=error),&
            n_var=-1,&
            unit_str="Ry",&
            error=error)
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
            " giving the cutoff used for a gaussian with alpha=1",&
            usage="RELATIVE_CUTOFF real", default_r_val=15.0_dp,&
            unit_str="Ry",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MULTIGRID_SET",&
            description="Activate a manual setting of the multigrids",&
            usage="MULTIGRID_SET", default_l_val=.FALSE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="MULTIGRID_CUTOFF",&
            variants=(/"CUTOFF_LIST"/),&
            description="List of cutoff values to set up multigrids manually",&
            usage="MULTIGRID_CUTOFF 200.0 100.0 ",&
            n_var=-1,type_of_var=real_t,repeats=.TRUE.,error=error)
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
  !***************************************************************************


  !!****f* input_cp2k_dft/create_interp_section *
  !!
  !!   NAME
  !!     create_interp_section
  !!
  !!   FUNCTION
  !!     creates the interpolation section
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     tlaino
  !!
  !!*** **********************************************************************
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

  !!****f* input_cp2k_dft/create_sic_section *
  !!
  !!   NAME
  !!     create_sic_section
  !!
  !!   FUNCTION
  !!     creates the sic (self interaction correction) section
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
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

  !!****f* input_cp2k_dft/create_field_section *
  !!
  !!   NAME
  !!     create_field_section
  !!
  !!   FUNCTION
  !!     creates the section for finite fields
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     JGH
  !!
  !!*** **********************************************************************
  SUBROUTINE create_field_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_field_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"field",&
            description="parameters for finite fields",&
            n_keywords=6, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ELECTRICAL_FIELD_STRENGTH",&
            variants=(/"EFIELD_STRENGTH"/),&
            description="Strength of electrical field",&
            usage="ELECTRICAL_FIELD_STRENGTH  1.0",&
            default_r_val=0._dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ELECTRICAL_FIELD_POLARISATION",&
            variants=(/"EFIELD_POLARISATION"/),&
            description="Polarisation vector of electrical field",&
            usage="ELECTRICAL_FIELD_POLARISIATION  0.0 0.0 1.0",&
            repeats=.FALSE.,required=.FALSE.,n_var=3,&
            type_of_var=real_t,default_r_vals=(/0.0_dp,0.0_dp,1.0_dp/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BERRY_PHASE",&
            description="Use the Berry phase definition for the polarisation",&
            usage="BERRY_PHASE ON/OFF", default_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


    END IF
  END SUBROUTINE create_field_section

  !!****f* input_cp2k_dft/create_ot_section *
  !!
  !!   NAME
  !!     create_ot_section
  !!   FUNCTION
  !!     makes the orbital transformation section
  !!   NOTES
  !!
  !!   INPUTS
  !!
  !!   MODIFICATION HISTORY
  !!     11.2004 created [Joost VandeVondele]
  !!
  !!   SOURCE
  !!*** **********************************************************************
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
            n_keywords=6, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            citations=(/VandeVondele2003/), error=error)

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
                          "of the Approximative Congruence transformation (OT/IRAC)."),&
            enum_i_vals=(/ot_algo_taylor_or_diag,ot_algo_irac/),&
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
            enum_c_vals=s2a( "SD", "CG", "DIIS"),&
            enum_desc=s2a("Steepest descent: not recommended","Conjugate Gradients: most reliable, use for difficult systems."//&
            " The total energy should decrease at every OT CG step if the line search is appropriate.", &
            "Direct inversion in the iterative subspace: less reliable than CG, but sometimes about 50% faster"),&
            enum_i_vals=(/ot_mini_sd,ot_mini_cg,ot_mini_diis/),&
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

       CALL keyword_create(keyword, name="N_DIIS",&
            variants=(/"NDIIS"/),&
            description="Number of history vectors to be used with DIIS",&
            usage="N_DIIS 4",&
            default_i_val=7,error=error)
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
            enum_c_vals=s2a("FULL_ALL","FULL_SINGLE_INVERSE","FULL_SINGLE","FULL_KINETIC","FULL_S_INVERSE","SPARSE_DIAG","NONE"),&
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
            "skip preconditioning"),&
            enum_i_vals=(/ot_precond_full_all,ot_precond_full_single_inverse,ot_precond_full_single, &
            ot_precond_full_kinetic,ot_precond_s_inverse,&
            ot_precond_sparse_diag,ot_precond_none/),&
            citations=(/VandeVondele2003/),&
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

    END IF
  END SUBROUTINE create_ot_section

  !!****f* input_cp2k_dft/create_xas_section *
  !!
  !!   NAME
  !!     create_xas_section
  !!   FUNCTION
  !!     makes the input section for core-level spectroscopy simulations
  !!
  !!   NOTES
  !!
  !!   INPUTS
  !!
  !!   MODIFICATION HISTORY
  !!     03.2005 created [MI]
  !!
  !!   SOURCE
  !!*** **********************************************************************
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
            n_keywords=10, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
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
            enum_desc=s2a("No XAS is calculated", "Transition potential halfhole",&
            "Transition potential fullhole", "Hole in homo for XES " , &
            "Transition potential excited halfhole", &
            "Transition potential excited fullhole " , &
            "DSCF calculations for the first excited state"),&
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
            enum_desc=s2a("Length form <i|e r |j>","Velocity form <i|d/dr|j>"),&
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

       CALL keyword_create(keyword, name="MIXING",&
            description="fraction of new density to be mixed in"//&
            " the xas-scf,"//&
            "if default, the scf_control value is used",&
            usage="MIXING 0.5",&
            default_r_val=-1._dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STATE_TYPE",&
            variants=(/"TYPE"/),&
            description="Type of the orbitas that are excited for the xas spectra calculation",&
            usage="STATE_TYPE 1S",&
            default_i_val=xas_1s_type,&
            enum_c_vals=s2a( "1S", "2S"),&
            enum_desc=s2a("1s orbitals", "2s orbitals"),&
            enum_i_vals=(/xas_1s_type,xas_2s_type/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STATE_SEARCH",&
            description="# of states where to look for the one to be excited",&
            usage="STATE_SEARCH 1",&
            default_i_val=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SMEAR",&
            description="smearing coefficient",&
            usage="SMEAR 0.5",&
            default_r_val=0._dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="ATOMS_LU_BOUNDS",&
            variants=(/"AT_LU"/),&
            description="The lower and upper index of the atoms to be excited",&
            usage="ATOMS_LU_BOUNDS integer integer",&
            n_var=2,default_i_vals=(/0,0/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="ATOMS_LIST",&
            variants=(/"AT_LIST"/),&
            description="Indexes of the atoms to be excited"//&
            "This keyword can be repeated several times"//&
            "(useful if you have to specify many indexes).",&
            usage="ATOMS_LIST 1 2",&
            n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
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

       CALL keyword_create(keyword, name="EMISSION",&
            description="Calculate also the emission spectrum",&
            usage="EMISSION",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
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

       NULLIFY(subsection)
       CALL create_localize_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL section_create(subsection,"PRINT",&
            "printing of information during the core-level spectroscopy simulation",&
            error=error,repeats=.FALSE.,required=.FALSE.)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing of basic iteration information in CLS", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SPECTRUM",&
            description="Controls the dumping of the CLS output files containing the spectra",&
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
            description="the stride used to write the cube file "//&
            "(the bigger, the smaller the cube gets)",&
            default_i_val=2, error=error)
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
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_xas_section

  !!*** **********************************************************************

  !!****f* input_cp2k_dft/create_linres_section [1.0] *
  !!
  !!   NAME
  !!     create_linres_section
  !!
  !!   FUNCTION
  !!     creates the input structure used to activate
  !!     a linear response calculation 
  !!     Available properties : none
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
  !!     MI
  !!
  !!*** **********************************************************************
  SUBROUTINE create_linres_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_linres_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
    NULLIFY(keyword,subsection,print_key)

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="linres",&
            description="The linear response is used to calculate one of the "//&
            " following properties: chemicalshift, ...(to be implemented) ",&
            n_keywords=5, n_subsections=2, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="EPS_SCF",&
            description="target accuracy for the scf convergence",&
            usage="EPS_SCF 1.e-6", default_r_val=1.e-5_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_SCF",&
            description="Maximum number of SCF iteration to be performed for one optimization",&
            usage="MAX_SCF 200", default_i_val=50,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="OPT_METHOD",&
            variants=(/"MINIMIZER"/),&
            description="Define the algorithms used to optimize the Psi1 in the scf iteration",&
            usage="OPT_METHOD CG",&
            default_i_val=lr_cg,&
            enum_c_vals=s2a( "CG", "DIIS", "SD"),&
            enum_desc=s2a("Conjugate Gradients", "DIIS (WARNING this option doesnt work yet)",&
            "Steepesst Descent"),&
            enum_i_vals=(/lr_cg,lr_diis,lr_sd/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="N_DIIS",&
            variants=(/"NDIIS"/),&
            description="Number of history vectors to be used with DIIS "//&
            "(WARNING this option doesnt work yet)",&
            usage="N_DIIS 4",&
            default_i_val=7,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LINESEARCH",&
            variants=(/"LINE_SEARCH"/),&
            description="1D line search algorithm,"//&
            " in increasing order of robustness and cost.",&
            usage="LINESEARCH GOLD",&
            default_i_val=ls_none,&
            enum_c_vals=s2a( "NONE", "2PNT", "3PNT","GOLD"),&
            enum_desc=s2a("take fixed lenght steps","extrapolate based on 2 points", &
            "... or on 3 points","perform 1D golden section search of the minimum "//&
            "(WARNING this option doesnt work yet)"),&
            enum_i_vals=(/ls_none,ls_2pnt,ls_3pnt,ls_gold/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STEPSIZE",&
            description="Initial stepsize used for the line search",&
            usage="STEPSIZE 0.4",&
            default_r_val=0.15_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PRECONDITIONER",&
            description="Type of preconditioner to be used with all minimization schemes."//&
            "They differ in effectiveness, cost of construction, cost of application.",&
            usage="PRECONDITIONER FULL_S_INVERSE",&
            default_i_val=lr_precond_none,&
            enum_c_vals=s2a("FULL_ALL","FULL_S_INVERSE","FULL_H0","NONE"),&
            enum_desc=s2a("State selective, similar to full_single",&
            "Cholesky inversion of S", "full H0","skip preconditioning"),&
            enum_i_vals=(/lr_precond_full_all,lr_precond_s_inverse,&
            lr_precond_h0_inverse,lr_precond_none/), error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ENERGY_GAP",&
            description="Energy gap estimate [a.u.] for preconditioning",&
            usage="ENERGY_GAP 0.1",&
            default_r_val=0.2_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHECK_ORTHO",&
            description="Require the check of the orthogonality condition",&
            usage="CHECK_ORTHO",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsection)
       CALL create_localize_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)


       NULLIFY(subsection)
       CALL create_nmr_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_epr_section(subsection,error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection,"PRINT","printing of information during the XAS calculation",&
            error=error,repeats=.FALSE.,required=.FALSE.)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing of basic iteration information during the XAS", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)       

    END IF
  END SUBROUTINE create_linres_section

  SUBROUTINE create_nmr_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_nmr_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
    NULLIFY(keyword,print_key, subsection)

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="nmr",&
            description="The chemical shift is calculated by DFPT. This feature is not operative yet, DO NOT USE IT ",&
            n_keywords=5, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="controls the activation of the nmr calculation",&
            usage="&NMR T",&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="GAUGE",&
            description="The gauge used to compute the chemical shifts within GAPW.",&
            usage="GAUGE R",&
            default_i_val=nmr_gauge_r_and_step_function,&
            enum_c_vals=s2a("R", "R_AND_BADER", "R_AND_STEP_FUNCTION"),&
            enum_desc=s2a("Position gauge (doesnt work well).",&
                          "Position and Bader gauges for the soft and the local parts, respectively.",&
                          "Position and step function for the soft and the local parts, respectively."),&
            enum_i_vals=(/nmr_gauge_r,nmr_gauge_r_and_Bader,nmr_gauge_r_and_step_function/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FULL",&
            description="Calculate the orbital dependent correction "//&
            " to get the full current density", usage="FULL",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NICS",&
            description="Calculate the chemical shift in a set of points  "//&
            " given from an external file", usage="NICS",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NICS_FILE_NAME",&
            description="Name of the file with the NICS points coordinates",&
            usage="NICS_FILE_NAME nics_file",&
            default_lc_val="nics_file",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="RESTART_NMR",&
            description="Restart the NMR calculation from a previous run",&
            usage="RESTART_NMR",default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsection)
       CALL section_create(subsection,name="PRINT",&
            description="print results of nmr calculation",&
            error=error,repeats=.FALSE.,required=.FALSE.)

       CALL cp_print_key_section_create(print_key,"CURRENT_CUBES",&
            description="Controls the printing of the current density ", &
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="the stride used to write the cube file "//&
            "(the bigger, the smaller the cube gets)",&
            default_i_val=2, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"RESPONSE_FUNCTION_CUBES",&
            description="Controls the printing of the response functions ", &
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="the stride used to write the cube file "//&
            "(the bigger, the smaller the cube gets)",&
            default_i_val=2, error=error)
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
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"CHI_TENSOR",&
            description="Controls the printing of susceptibility",&
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SHIELDING_TENSOR",&
            description="Controls the printing of the chemical shift",&
            print_level=low_print_level,add_last=add_last_numeric,filename="",&
            error=error)

       CALL keyword_create(keyword,name="ATOMS_LU_BOUNDS",&
            variants=(/"ATOMS_LU"/),&
            description="The lower and upper atomic index for which the tensor is printed",&
            usage="ATOMS_LU_BOUNDS integer integer",&
            n_var=2,default_i_vals=(/0,-2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS_LIST",&
            description="list of atoms for which the shift is printed into a file ",&
            usage="LIST_ATOMS 1 2",n_var=-1,&
            type_of_var=integer_t,repeats=.TRUE., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_interp_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_nmr_section

  SUBROUTINE create_epr_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_epr_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection, &
                                                subsubsection

    failure=.FALSE.
    NULLIFY(keyword,print_key, subsection, subsubsection)

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="epr",&
            description="The g tensor is calculated by DFPT ",&
            n_keywords=5, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
           description="controls the activation of the epr calculation",&
           usage="&EPR T",&
           default_l_val=.FALSE.,&
           lone_keyword_l_val=.TRUE.,&
           error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="GAUGE",&
            description="The gauge used to compute the g tensor within GAPW.",&
            usage="GAUGE R",&
            default_i_val=epr_gauge_r_and_step_function,&
            enum_c_vals=s2a("R", "R_AND_BADER", "R_AND_STEP_FUNCTION"),&
            enum_desc=s2a("Position gauge (doesnt work well).",&
                          "Position and Bader gauges for the soft and the local parts, respectively.",&
                          "Position and step function for the soft and the local parts, respectively."),&
            enum_i_vals=(/epr_gauge_r,epr_gauge_r_and_Bader,epr_gauge_r_and_step_function/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FULL",&
            description="Calculate the orbital dependent correction "//&
            " to get the full current density", usage="FULL",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_EPR",&
            description="Restart the EPR calculation from a previous run",&
            usage="RESTART_EPR",default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsection)
       CALL section_create(subsection,name="PRINT",&
            description="print results of epr calculation",&
            error=error,repeats=.FALSE.,required=.FALSE.)
 
       CALL cp_print_key_section_create(print_key,"CURRENT_CUBES",&
            description="Controls the printing of the current density ", &
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="the stride used to write the cube file "//&
            "(the bigger, the smaller the cube gets)",&
            default_i_val=2, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
 
       CALL cp_print_key_section_create(print_key,"NABLAVKS_CUBES",&
            description="Controls the printing of the components of nabla v_ks ", &
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="the stride used to write the cube file "//&
            "(the bigger, the smaller the cube gets)",&
            default_i_val=2, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"G_TENSOR",&
            description="Controls the printing of the g tensor",&
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL create_xc_section(subsubsection,error)
       CALL section_add_subsection(print_key, subsubsection, error=error)
       CALL section_release(subsubsection,error=error)

       CALL keyword_create(keyword, name="GAPW_MAX_ALPHA",&
            description="Maximum alpha of GTH potentials allowed on the soft grids ",&
            usage="GAPW_MAX_ALPHA real", default_r_val=5.0_dp,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SOO_RHO_HARD",&
            description="Whether or not to include the atomic parts of the density "//&
            "in the SOO part of the g tensor", usage="SOO_RHO_HARD", &
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
 
       CALL cp_print_key_section_create(print_key,"RESPONSE_FUNCTION_CUBES",&
            description="Controls the printing of the response functions ", &
            print_level=high_print_level,add_last=add_last_numeric,filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="the stride used to write the cube file "//&
            "(the bigger, the smaller the cube gets)",&
            default_i_val=2, error=error)
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
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
 
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
 
       NULLIFY(subsection)
       CALL create_interp_section(subsection,error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
 
     END IF

  END SUBROUTINE create_epr_section


  !!****f* input_cp2k_dft/create_kg_section *
  !!
  !!   NAME
  !!     create_kg_section
  !!
  !!   FUNCTION
  !!     creates a section with the parameters specific to the kg method
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - section:the section to create
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
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
       CALL section_create(section,name="kg",&
            description="Parameters specific to KG methods",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY ( keyword, subsection, print_key )

       CALL section_create(subsection,name="MOL_CHARGE",&
            description="This section is used to set a fixed charge for a molecule kind.",&
            n_keywords=3, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       ! Integer
       CALL keyword_create(keyword, name="MOLECULE",&
            variants=(/"MOL"/),&
            description="Specifies the number of the molecule kind on which the charge is applied."//&
            " MOLECULE and MOLNAME keyword exclude themself mutually.",&
            usage="MOL {integer}", required=.TRUE., n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MOLNAME",&
            variants=(/"SEGNAME"/),&
            description="Specifies the name of the molecule kind on which the charge is applied.",&
            usage="MOLNAME {character}", required=.TRUE.,&
            n_var=1, type_of_var=char_t, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CHARGE",&
            description="Value of the charge for this molecule kind",&
            usage="CHARGE {integer}",default_i_val=0,&
            required=.TRUE.,n_var=1,type_of_var=integer_t,error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL section_create(subsection,name="PRINT",&
            description="print results of kg_gpw calculation",&
            error=error,repeats=.FALSE.,required=.FALSE.)

       CALL cp_print_key_section_create(print_key,"kg_info",&
            description="",&
            print_level=medium_print_level, filename="__STD_OUT__",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MOMENTS",&
            description="",print_level=high_print_level, filename="",error=error)

       CALL keyword_create(keyword, name="NMOMENTS",&
              description="Moment order",&
              usage="NMOMENTS {integer}",default_i_val=2,&
              required=.FALSE.,n_var=1,type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
    END IF

  END SUBROUTINE create_kg_section


  !!****f* input_cp2k_dft/create_et_coupling_section
  !!
  !!   NAME
  !!     create_et_coupling_section
  !!
  !!   FUNCTION
  !!     creates the section for electron transfer coupling
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fschiff
  !!
  !!*** **********************************************************************
  SUBROUTINE create_et_coupling_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_et_coupling_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"ET_COUPLING",&
            description="specifies the two constraints/restraints for extracting ET coupling elements",&
            n_keywords=1, n_subsections=4, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(subsection)
       CALL create_restraint_A(subsection,"DDAPC_RESTRAINT_A",error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
       
       NULLIFY(subsection)
       CALL create_restraint_A(subsection,"DDAPC_RESTRAINT_B",error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
        
       NULLIFY(subsection)
       CALL create_restraint_A(subsection,"BECKE_CONSTRAINT_A",error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(subsection)
       CALL create_restraint_A(subsection,"BECKE_CONSTRAINT_B",error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL keyword_create(keyword, name="TYPE_OF_CONSTRAINT",&
            description="Specifies the type of constraint",&
            usage="TYPE_OF_CONSTRAINT DDAPC",&
            enum_c_vals=s2a("NONE","DDAPC","BECKE"),&
            enum_i_vals=(/ do_no_et, do_et_ddapc, do_et_becke/),&
            enum_desc=s2a("NONE","ddapc_restraint","Sperical potential"),&
            default_i_val=do_no_et,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing basic info about the method", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_et_coupling_section 

  SUBROUTINE create_restraint_A(section,section_name,error)
    TYPE(section_type), POINTER              :: section
    CHARACTER(len=*), INTENT(in)             :: section_name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_restraint_A', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    NULLIFY(keyword, print_key)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,TRIM(ADJUSTL(section_name)),&
            description="Use DDAPC charges in a restraint (check code for details)",&
            n_keywords=7, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
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
            description="Defines the the coefficient of the atom in the atom list (default is one), currently DDAPC only ",&
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

  END SUBROUTINE create_restraint_A

  !***************************************************************************

END MODULE input_cp2k_dft

!***************************************************************************
