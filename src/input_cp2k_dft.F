!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
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
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE cp_units,                        ONLY: cp_unit_to_cp2k
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
  USE input_val_types,                 ONLY: integer_t,&
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
  PUBLIC :: create_bsse_section
  PUBLIC :: create_properties_section, create_scf_section
!***
!****************************************************************************
CONTAINS

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
         usage="BASIS_SET_FILE_NAME ../BASIS_SET",&
         default_lc_val="BASIS_SET",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="POTENTIAL_FILE_NAME",&
         description="Name of the pseudo potential files, may include a path",&
         usage="POTENTIAL_FILE_NAME ../POTENTIAL",&
         default_lc_val="POTENTIAL",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword, name="RESTART_FILE_NAME",&
         description="Name of the restart file, may include a path",&
         usage="RESTART_FILE_NAME ../RESTART",&
         default_lc_val="RESTART",error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL keyword_create(keyword=keyword,&
                        name="UKS",&
                        variants=s2a("UNRESTRICTED_KOHN_SHAM",&
                                     "LSD",&
                                     "SPIN_POLARIZED"),&
                        description="Requests a spin-polarized calculation using alpha "//&
                                    "and beta orbitals, i.e. no spin restriction is applied",&
                        usage="lsd",&
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
    CALL keyword_create(keyword, name="RESTRICTED",&
         description="If the calculation should be (spin) restricted",&
         usage="RESTRICTED",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
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

    CALL create_poisson_section(subsection,error)
    CALL section_add_subsection(section, subsection, error=error)
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_localize_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

      NULLIFY(keyword)
      CALL section_create(section,name="LOCALIZE",&
           description="Use one of the available methods to define the localization "//&
           " and possibly to optimize it to a minimum or a maximum.",&
           n_keywords=8, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
           error=error)

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
           usage="METHOD (JACOBI|CRAZY|DIRECT|NONE)",&
           enum_c_vals=s2a("NONE","JACOBI","CRAZY","DIRECT"),&
           enum_i_vals=(/do_loc_jacobi, do_loc_crazy, do_loc_direct, do_loc_none/),&
           default_i_val=do_loc_jacobi,error=error)
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
 
      CALL keyword_create(keyword,name="SET_OF_STATES",&
           description="Selected states to be localized",&
           usage="SET_OF_STATES (ALL|UPTON|FROMITOJ|LIST|NONE)",&
           enum_c_vals=s2a("ALL","UPTON","FROMITOJ","LIST","NONE"),&
           enum_i_vals=(/state_loc_all, state_loc_upton, state_loc_fromitoj,& 
           state_loc_list,state_loc_none/),&
           default_i_val=state_loc_all,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)


      CALL keyword_create(keyword,name="LU_BOUNDS",&
           description="The lower and upper index of the set of wfn to be localized",&
           usage="LU_BOUNDS integer integer",&
           n_var=2,default_i_vals=(/0,0/),type_of_var=integer_t,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,name="UP_BOUND",&
           description="Highest index of  wfn to be localized",&
           usage="UP_BOUND  integer",&
           default_i_val=0,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)


      CALL keyword_create(keyword,name="LISTSTATES",&
           description="Indexes of the wfn to be localized"//&
           "This keyword can be repeated several times"//&
           "(useful if you have to specify many indexes).",&
           usage="LISTSTATES 1 2",&
           n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)
 
      CALL keyword_create(keyword,name="EACH_STEP",&
           description="Evry how many steps the localization calculation"//&
           " is performed  along an iterative procedure (e.g. MD)",&
           usage="EACH_STEP integer",&
           default_i_val=1,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="FNAME",&
           description="Root of the file names where to store the spectra",&
           usage="FNAME loc_orb",&
           default_lc_val="loc_orb",error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_print_dft_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, print_key2

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options in DFT code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(print_key, keyword, print_key2)

       CALL cp_print_key_section_create(print_key,"PROGRAM_BANNER",&
            description="Controls the printing of the banner of the MM program",&
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)       

       CALL cp_print_key_section_create(print_key,"derivatives",&
            description="Print all derivatives after the DFT calculation", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)       

       CALL cp_print_key_section_create(print_key_section=print_key,&
                                        name="neighbor_lists",&
                                        description="Controls the printing of the neighbor lists",&
                                        print_level=debug_print_level,&
                                        each=(/1/),&
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

       CALL cp_print_key_section_create(print_key,"ao_matrices",&
            description="Controls the printing of the ao (i.e. contracted gaussian) matrices (debug).", &
            print_level=debug_print_level,each=(/1/),filename="__STD_OUT__",&
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
            description="Controls the printing of the molecular orbitals .", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
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
            print_level=high_print_level,each=(/1/),filename="",&
            error=error)
       CALL keyword_create(keyword, name="stride",&
            description="the stride used to write the cube file "//&
            "(the bigger, the smaller the cube gets)",&
            default_i_val=2, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL keyword_create(keyword, name="nlumo",&
            description="If the printkey is activated controls the number of lumos that are printed (-1=all)",&
            default_i_val=5, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)    
       CALL keyword_create(keyword, name="nhomo",&
            description="If the printkey is activated controls the number of homos that are calculated and printed",&
            default_i_val=1, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)


       CALL cp_print_key_section_create(print_key,"gapw",&
            description="Controls the printing of some gapw related information (debug).", &
            print_level=debug_print_level,each=(/1/),filename="__STD_OUT__",&
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
            print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"OPTICAL_CONDUCTIVITY",&
            description="Controls the printing of the optical conductivity.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"e_density_cube",&
            description="Controls the printing of a cube file with electronic density.", &
            print_level=high_print_level,each=(/1/),filename="",&
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
            print_level=high_print_level,each=(/1/),filename="",&
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
            print_level=high_print_level,each=(/1/),common_iter_levels=1,&
            filename="__STD_OUT__", error=error)
       CALL keyword_create(keyword, name="LOC_TYPE",&
            description="Specify the method to use for the localization.",&
            usage="LOC_TYPE (JACOBI|CRAZY|DIRECT|NONE)",&
            enum_c_vals=s2a("JACOBI","CRAZY","DIRECT","NONE"),&
            enum_i_vals=(/ do_loc_jacobi, do_loc_crazy, do_loc_direct, do_loc_none/),&
            default_i_val=do_loc_jacobi, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="TOTAL_DIPOLE",&
            description="prints out the total dipole of the wannier centers.",&
            usage="TOTAL_DIPOLE", default_l_val=.false., &
            lone_keyword_l_val=.true., error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key2,"cubes",&
            description="Controls the printing of cube files in localization", &
            print_level=high_print_level,each=(/1/),filename="",&
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

       CALL cp_print_key_section_create(print_key,"MOLECULAR_STATES",&
            description="Controls the printing of the molecular states", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            error=error)

       CALL cp_print_key_section_create(print_key2,"cubes",&
            description="Controls the printing of cube files in localization", &
            print_level=high_print_level,each=(/1/),filename="",&
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

       CALL cp_print_key_section_create(print_key,"DIPOLE",&
            description="Calculates the dipole, for the moment it only supports"//&
            "the stupid non periodic dipole, should be extended soon", &
            print_level=high_print_level,each=(/1/),filename=" ",&
            common_iter_levels=1,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MULLIKEN",&
            description="Controls the printing of the Mulliken Population analysis.", &
            print_level=medium_print_level,each=(/1/),filename="__STD_OUT__",&
            common_iter_levels=1, error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       ! fairly expensive to compute, so only at high
       CALL cp_print_key_section_create(print_key,"LOWDIN",&
            description="Controls the printing of the Lowdin charges.", &
            print_level=high_print_level,each=(/1/),filename="__STD_OUT__",&
            common_iter_levels=1, error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(&
             print_key_section=print_key,&
             name="XRAY_DIFFRACTION_SPECTRUM",&
             description="Calculate and print the coherent X-ray "//&
                         "diffraction spectrum",&
             print_level=debug_print_level,&
             each=(/1/),&
             filename="",&
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
          description="Guarantee an exact derivative (Hks) of the energy",&
          usage="MAP_CONSISTENT",&
          default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="LOCALIZE_WFN",&
          description="Check and possibly optimize the localization of the wavefunctions",&
          usage="LOCALIZE_WFN",&
          default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error) 

     NULLIFY(subsection)
     CALL create_localize_section(subsection,error)
     CALL section_add_subsection(section, subsection, error=error)
     CALL section_release(subsection,error=error)


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

     CALL keyword_create(keyword, name="RS_GRID",&
          description="Parallellisation strategy of the realspace grids",&
          usage="RS_GRID DISTRIBUTED",&
          enum_c_vals=s2a("DISTRIBUTED","REPLICATED"),&
          enum_desc=s2a("should be faster in principle, is often not, rarely tested", "default"),&
          enum_i_vals=(/ do_rsgrid_distributed, do_rsgrid_replicated/),&
          default_i_val=do_rsgrid_replicated, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EXTRAPOLATION",&
          variants=s2a("INTERPOLATION","WF_INTERPOLATION"),&
          description="Extrapolation strategy for e.g. MD",&
          usage="EXTRAPOLATION PS",&     ! XXXXXXXXXXX how to add the order ? i.e. PS 3
          enum_c_vals=s2a("USE_GUESS","USE_PREV_P","USE_PREV_RHO_R","LINEAR_WF",&
                          "LINEAR_P","LINEAR_PS","USE_PREV_WF","PS","FROZEN","ASPC"),&
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
          default_i_val=do_extrpl_linear_p, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EXTRAPOLATION_ORDER",&
          description="Order for the PS extrapolation",&
          usage="EXTRAPOLATION_ORDER {integer}",default_i_val=0, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)



     CALL keyword_create(keyword, name="METHOD",&
          description="What method should be employed",&
          usage="METHOD GAPW",&
          enum_c_vals=s2a("GAPW","GAPW_XC","GPW","KG_POL","KG_NOPOL","KG_GPW","EHT",&
                          "MNDO","AM1","PM3","PDG"),&
          enum_desc=s2a("Gaussian and augmented plane waves method",&
                        "Gaussian and augmented plane waves method only for XC",&
                        "Gaussian and plane waves method",&
                        "Polarisable Kim-Gordon method",&
                        "Non-polarisable Kim-Gordon method", " DFT Embedding method ",&
                        "EHT semiempirical","MNDO semiempirical","AM1 semiempirical",&
                        "PM3 semiempirical","PDG semiempirical"),&
          enum_i_vals=(/ do_method_gapw, do_method_gapw_xc, do_method_gpw, do_method_kg_pol,&
                         do_method_kg_nopol, do_method_kg_gpw, do_method_eht,&
                         do_method_mndo, do_method_am1, do_method_pm3, do_method_pdg/),&
          default_i_val=do_method_gpw, error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

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
     CALL create_s2_restraint_section(subsection,error)
     CALL section_add_subsection(section, subsection, error=error)
     CALL section_release(subsection,error=error)

     NULLIFY(subsection)
     CALL create_harris_section(subsection,error)
     CALL section_add_subsection(section, subsection, error=error)
     CALL section_release(subsection,error=error) 

  END IF

END SUBROUTINE create_qs_section

SUBROUTINE create_se_control_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
           description="Assume orthogonal basis set",&
           usage="ORTHOGONAL_BASIS ",default_l_val=.FALSE., error=error)
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

    END IF

END SUBROUTINE create_se_control_section

SUBROUTINE create_mulliken_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
           usage="LIST {integer} {integer} .. {integer}", required=.TRUE.,&
           n_var=-1, type_of_var=integer_t, error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

    END IF

END SUBROUTINE create_mulliken_section

SUBROUTINE create_ddapc_restraint_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

      CALL keyword_create(keyword, name="FUNCTIONAL_FORM",&
           description="Specifies the functional form of the term added",&
           usage="FUNCTIONAL_FORM RESTRAINT",&
           enum_c_vals=s2a("RESTRAINT","CONSTRAINT"),&
           enum_i_vals=(/ do_ddapc_restraint, do_ddapc_constraint/),&
           enum_desc=s2a("Harmonic potential: s*(q-t)**2","Constraint form: s*(q-t)"),&
           default_i_val=do_ddapc_restraint, error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

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
            print_level=high_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)       

    END IF

END SUBROUTINE create_ddapc_restraint_section

SUBROUTINE create_s2_restraint_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
         error=error)
     
     NULLIFY(subsection,keyword)
     CALL keyword_create(keyword,name="_SECTION_PARAMETERS_",&
          description="shorcut for the most commont functionals combinations",&
          usage="&xc_functional BLYP",&
          enum_c_vals=s2a("BLYP","PADE","PBE","TPSS","HCTH120","OLYP","NO_SHORTCUT"),&
          enum_i_vals=(/xc_funct_blyp,xc_funct_pade,xc_funct_pbe,&
          xc_funct_tpss,xc_funct_hcth120,xc_funct_olyp,xc_funct_no_shortcut/),&
          default_i_val=xc_funct_no_shortcut,lone_keyword_i_val=xc_funct_no_shortcut,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL section_create(subsection,name="BECKE88",&
          description="Uses the Becke 88 exchange functional",&
          n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
     CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
          description="activates the functional",&
          lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
     CALL section_add_keyword(subsection,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)

     CALL section_create(subsection,name="LYP",&
          description="Uses the LYP functional",&
          n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
          error=error)
     CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
          description="activates the functional",&
          lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
     CALL section_add_keyword(subsection,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)

     CALL section_create(subsection,name="PADE",&
          description="Uses the PADE functional",&
          n_keywords=0, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
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
          error=error)
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
          error=error)
     CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
          description="activates the functional",&
          lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
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
          error=error)
     CALL keyword_create(keyword,"_SECTION_PARAMETERS_",&
          description="activates the functional",&
          lone_keyword_l_val=.TRUE.,default_l_val=.FALSE.,error=error)
     CALL section_add_keyword(subsection,keyword,error=error)
     CALL keyword_release(keyword,error=error)
     CALL keyword_create(keyword,"parametrization",&
          description="switches between the original pbe and the revised pbe "//&
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
            "SPLINE2_SMOOTH", "NN6_SMOOTH", "SPLINE3_SMOOTH", "NN4_SMOOTH","FD"),&
            enum_i_vals=(/xc_deriv_pw, xc_deriv_spline3, xc_deriv_spline2,&
            xc_deriv_nn50_smooth, xc_deriv_nn10_smooth, xc_deriv_spline2_smooth,&
            xc_deriv_spline2_smooth, xc_deriv_spline3_smooth, xc_deriv_spline3_smooth,xc_deriv_fd/),&
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
            n_keywords=24, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(subsection)
       CALL create_ot_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_outer_scf_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="WRITE_RESTART_EACH",&
            description="How often the wavefunction restart file should be "//&
            "written during an scf (0=only at the end, -1=never).",usage="WRITE_RESTART_EACH 5",&
            default_i_val=10,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

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

       CALL keyword_create(keyword, name="MAX_SCF_ASPC", &
            description="Maximum number of SCF iterations to be performed for ASPC", &
            usage="MAX_SCF_ASPC 1", required=.FALSE., default_i_val=1,error=error)
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

       CALL keyword_create(keyword, name="SCF_GUESS",&
            description="Change the initial guess for the wavefunction.",&
            usage="SCF_GUESS RESTART", default_i_val=atomic_guess,&
            enum_c_vals=s2a("ATOMIC","RESTART","RANDOM","CORE","DENSITIES"),&
            enum_i_vals=(/atomic_guess,restart_guess,random_guess,core_guess,&
            densities_guess/),&
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
                           name="ROKS_PARAMETERS",&
                           variants=(/"ROKS_PARAMETER"/),&
                           description="Allows to define all ROKS parameters explicitly. "//&
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
            error=error,repeats=.FALSE.,required=.FALSE.)

       CALL cp_print_key_section_create(print_key,"iteration_info",&
            description="Controls the printing of basic iteration information during the scf.", &
            print_level=low_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
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
            print_level=low_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error) 

       CALL cp_print_key_section_create(print_key,"MO_ORTHONORMALITY",&
            description="Controls the printing relative to the orthonormality of MOs (CT S C).", &
            print_level=high_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)      

       CALL cp_print_key_section_create(print_key,"MO_MAGNITUDE",&
            description="Prints the min/max eigenvalues of the overlap of the MOs without S (CT C).", &
            print_level=high_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)      

       CALL cp_print_key_section_create(print_key,"detailed_energy",&
            description="Controls the printing of detailed energy information.", &
            print_level=medium_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"diis_info",&
            description="Controls the printing of diis information.", &
            print_level=high_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"total_densities",&
            description="Controls the printing of total densities.", &
            print_level=medium_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_outer_scf_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"OUTER_SCF",&
            description="parameters controlling the outer scf loop ",&
            n_keywords=8, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="TYPE",&
            description="Specifies which kind of outer scf should be employed",&
            usage="TYPE DDAPC_CONSTRAINT ",&
            default_i_val=outer_scf_none,&
            enum_c_vals=s2a( "DDAPC_CONSTRAINT", "S2_CONSTRAINT", "NONE" ),&
            enum_desc=s2a("Enforce a constraint on the DDAPC, requires the corresponding section", &
                          "Enforce a constraint on the S2, requires the corresponding section", &
                          "Do nothing in the outer loop, useful for resetting the inner loop,"),&
            enum_i_vals=(/outer_scf_ddapc_constraint,outer_scf_s2_constraint,outer_scf_none/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="OPTIMIZER",&
            description="Method used to bring the outer loop to a stationary point",&
            usage="OPTIMIZER SD",&
            default_i_val=outer_scf_optimizer_none,&
            enum_c_vals=s2a("SD","DIIS","NONE"),&
            enum_desc=s2a("Takes steps in the direction of the gradient, multiplied by step_size", &
                          "Uses a Direct Inversion in the Iterative Subspace method", &
                          "Do nothing, useful only with the none type"),&
            enum_i_vals=(/outer_scf_optimizer_sd,outer_scf_optimizer_diis,outer_scf_optimizer_none/),&
            error=error)
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

       CALL create_fitcharge_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
      
    END IF

  END SUBROUTINE create_properties_section
!***************************************************************************

!!****f* input_cp2k_dft/create_fitcharge_section [1.0] *
!!
!!   NAME
!!     create_fitcharge_section
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
  SUBROUTINE create_fitcharge_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_fitcharge_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL cp_print_key_section_create(section,name="FIT_CHARGE",&
            description="This section is used to set up calculation of density derived atomic point charges.",&
            print_level=high_print_level,filename="__STD_OUT__",error=error)

       NULLIFY(keyword, print_key )

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
            n_var=1, type_of_var=real_t, default_r_val=0.5_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       
       CALL keyword_create(keyword, name="RADII",&
            description="Specifies all the radius of the gaussian used in the fit for each atomic site. The use"//&
            " of this keyword disables all other keywords of this section.",&
            usage="RADII {real} {real} .. {real}", required=.FALSE.,&
            n_var=-1, type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GCUT",&
            description="Cutoff for fit in G-space.",&
            usage="GCUT {real}", required=.FALSE.,&
            n_var=1, type_of_var=real_t, default_r_val=SQRT(6.0_dp),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing of basic information during the run", &
            print_level=low_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF

  END SUBROUTINE create_fitcharge_section

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_bsse_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="BSSE",&
            description="This section is used to set up the BSSE calculation.",&
            n_keywords=3, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, subsection)
       CALL section_create(subsection,name="FRAGMENT",&
            description="Specify the atom number belonging to this fragment.",&
            n_keywords=2, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       ! Integer
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
            print_level=silent_print_level,each=(/1/),filename="__STD_OUT__",&
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
            print_level=medium_print_level,each=(/10/),filename="__STD_OUT__",&
            add_last=(/1/),error=error)
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_harris_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="HARRIS",&
            description="This section specifies the flags for the calculation of the harris functional",&
            n_keywords=3, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       
       CALL keyword_create(keyword, name="ACTIVATE",&
            description="Activates the harris functional.",&
            usage="ACTIVATE TRUE", required=.FALSE.,&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., &
            error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="FILL_WF_HISTORY", &
            description="Converges to self-consistency until the " //&
            "WF-history is filled.", usage="FILL_WF_HISTORY TRUE", &
            required=.FALSE., default_l_val=.FALSE., &
            lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="HARRIS_MAX_SCF", &
            description="The maximum number of SCF cycles after the wf_history is filled.", &
            usage="HARRIS_MAX_SCF 1", required=.FALSE., &
            default_i_val=1, lone_keyword_i_val=1, &
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
       CALL keyword_create(keyword, name="ngrids",&
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
            enum_c_vals=s2a("pw","linear","spline3_nopbc","spline3"),&
            enum_i_vals=(/pw_interp,linear_interp,&
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
            print_level=medium_print_level,each=(/10/),filename="__STD_OUT__",&
            add_last=(/1/),error=error)
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
         error=error)

     NULLIFY(keyword)

     CALL keyword_create(keyword, name="SIC_SCALING_A",&
          description="Scaling of the coulomb term in sic [experimental]",&
          usage="SIC_SCALING_A 0.5",&
          default_r_val=1.0_dp,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="SIC_SCALING_B",&
          description="Scaling of the xc term in sic [experimental]",&
          usage="SIC_SCALING_B 0.5",&
          default_r_val=1.0_dp,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="SIC_METHOD",&
          description="Method used to remove the self interaction",&
          usage="SIC_METHOD MAURI_US",&
          default_i_val=sic_none,&
          enum_c_vals=s2a( "NONE", "MAURI_US", "MAURI_SPZ", "AD"),&
          enum_i_vals=(/sic_none,sic_mauri_us,sic_mauri_spz,sic_ad/),&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

  END IF
END SUBROUTINE create_sic_section

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_ot_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

  failure=.FALSE.
  
  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL section_create(section,"ot",&
         description="Sets the various options for the orbital transformation (OT) method. "//&
                     "Default settings already provide an efficient, yet robust method. "//&
                     "Well-behaved systems might benefit from using a DIIS minimizer.",&
         n_keywords=6, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
         error=error)

     NULLIFY(keyword)

     CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
          description="controls the activation of the ot method",&
          usage="&OT T",&
          default_l_val=.FALSE.,&
          lone_keyword_l_val=.TRUE.,&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="MINIMIZER",&
          description="Minimizer to be used with the OT method",&
          usage="MINIMIZER DIIS",&
          default_i_val=ot_mini_cg,&
          enum_c_vals=s2a( "SD", "CG", "DIIS"),&
          enum_desc=s2a("Steepest descent","Conjugate Gradients", &
                        "Direct inversion in the iterative subspace"),&
          enum_i_vals=(/ot_mini_sd,ot_mini_cg,ot_mini_diis/),&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="SAFE_DIIS",&
            variants=(/"SAFER_DIIS"/),&
            description="Reject DIIS steps if they point away from the"//&
            " minimum, do SD in that case ",&
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
                      " LINESEARCH GOLD should always find an electronic minimum.",&
          usage="LINESEARCH GOLD",&
          default_i_val=ot_ls_2pnt,&
          enum_c_vals=s2a( "NONE", "2PNT", "3PNT","GOLD"),&
          enum_desc=s2a("take fixed lenght steps","extrapolate based on 2 points", &
                        "... or on 3 points","perform 1D golden section search of the minimum"),&
          enum_i_vals=(/ot_ls_none,ot_ls_2pnt,ot_ls_3pnt,ot_ls_gold/),&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="STEPSIZE",&
          description="Initial stepsize used for the line search",&
          usage="STEPSIZE 0.4",&
          default_r_val=0.15_dp,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="GOLD_TARGET",&
          description="Target relative uncertainty in the location of the minimum",&
          usage="GOLD_TARGET 0.1",&
          default_r_val=0.01_dp,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="PRECONDITIONER",&
          description="Type of preconditioner to be used with all minimization schemes."//&
                      "They differ in effectiveness, cost of construction, cost of application.",&
          usage="PRECONDITIONER FULL_SINGLE",&
          default_i_val=ot_precond_full_kinetic,&
          enum_c_vals=s2a("FULL_ALL","FULL_SINGLE","FULL_SINGLE_INVERSE","FULL_KINETIC","FULL_S_INVERSE","SPARSE_DIAG","NONE"),&
          enum_desc=s2a("State selective, similar to full_single","Based on diagonalisation", "Based on H-eS cholesky inversion", &
                        "Cholesky inversion of S and T", &
                        "Cholesky inversion of S","Only based on atomic blocks", &
                        "skip preconditioning"),&
          enum_i_vals=(/ot_precond_full_all,ot_precond_full_single,ot_precond_full_single_inverse, &
                        ot_precond_full_kinetic,ot_precond_s_inverse,&
                        ot_precond_sparse_diag,ot_precond_none/),&
          error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="ENERGY_GAP",&
          description="Energy gap estimate [a.u.] for preconditioning",&
          usage="EPS_TAYLOR 0.1",&
          default_r_val=0.2_dp,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="EPS_TAYLOR",&
          variants=(/"EPSTAYLOR"/),&
          description="Target accuracy of the taylor expansion for the matrix functions",&
          usage="EPS_TAYLOR 1.0E-15",&
          default_r_val=1.0E-16_dp,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="MAX_TAYLOR",&
          description="Maximum order of the Taylor expansion before diagonalisation is prefered",&
          usage="MAX_TAYLOR 5",&
          default_i_val=4,error=error)
     CALL section_add_keyword(section,keyword,error=error)
     CALL keyword_release(keyword,error=error)

     CALL keyword_create(keyword, name="ROTATION",&
          description="Introduce additional variables so that rotations of the occupied"//&
                      " subspace are allowed as well",&
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
!!     makes the x-ray absorption section
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_xas_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
  
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CALL section_create(section,"xas",&
         description="Sets the method of choice to calculate the XAS. "//&
                     "The occupied states from  which we calculate the "//&
                     "excitation should be specified. "//&
                     "Localization of the orbitals may be useful.",&
         n_keywords=10, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
         error=error)

      NULLIFY(keyword,subsection,print_key)

      CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
          description="controls the activation of the xes/xas calculation",&
          usage="&XAS T",&
          default_l_val=.FALSE.,&
          lone_keyword_l_val=.TRUE.,&
          error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)


      CALL keyword_create(keyword, name="METHOD",&
           variants=(/"XAS_METHOD"/),&
           description="Method to be used to calculate XAS",&
           usage="METHOD TP_HH",&
           default_i_val=xas_none,&
           enum_c_vals=s2a( "NONE", "TP_HH", "TP_FH", "TP_VAL", "TP_XHH","TP_XFH", "TDDFT","DSCF"),&
           enum_desc=s2a("No XAS is calculated", "Transition potential halfhole",&
                         "Transition potential fullhole", "Hole in homo for XES " , &
                         "Transition potential excited halfhole", &
                         "Transition potential excited fullhole " , &
                         "Excitations by TDDFT","DSCF calculations for the first excited state"),&
           enum_i_vals=(/xas_none,xas_tp_hh,xas_tp_fh,xes_tp_val,xas_tp_xhh,&
           xas_tp_xfh,xas_tddft,xas_dscf/), error=error)
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
           default_i_val=xas_scf_default,&
           enum_c_vals=s2a( "DEFAULT", "GENERAL", "OT"),&
           enum_desc=s2a("same as in std SCF", "diagonalization",&
           "adapted orbital transformation"),&
           enum_i_vals=(/xas_scf_default,xas_scf_general,xas_scf_ot/),&
           error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)
  
      CALL keyword_create(keyword, name="MAXSTEP",&
           description="Max # of steps in the xas-scf for one excitation",&
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
           " using DIAG/DIIS for the xas-scf"//&
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
 
      CALL keyword_create(keyword, name="EACH_STEP",&
           description="How often  the xas spectra are computed",&
           usage="EACH_STEP 1",&
           default_i_val=1,error=error)
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

      CALL keyword_create(keyword, name="FNAME",&
           description="Root of the file names where to store the spectra",&
           usage="FNAME xas_1s_o",&
           default_lc_val="xas_out",error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="SMEAR",&
           description="smearing coefficient",&
           usage="SMEAR 0.5",&
           default_r_val=0._dp,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="FERMI",&
           description="Fermi distributionfor the occupation numbers",&
           usage="FERMI",&
           default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
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

      CALL keyword_create(keyword, name="LOCALIZE_WFN",&
           description="Check and possibly optimize the localization of the wavefunctions",&
           usage="LOCALIZE_WFN",&
           default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error) 

      CALL keyword_create(keyword, name="CUBES",&
           description="Print some of the wavefunctions as cube files",&
           usage="CUBES",&
           default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error) 

      CALL keyword_create(keyword,name="CUBES_LU_BOUNDS",&
           variants=(/"CUBES_LU"/),&
           description="The lower and upper index of the states to be printed as cube files",&
           usage="CUBES_LU_BOUNDS integer integer",&
           n_var=2,default_i_vals=(/0,-2/), type_of_var=integer_t,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,name="CUBES_LIST",&
           description="Indexes of the states to be printed as cube files"//&
           "This keyword can be repeated several times"//&
           "(useful if you have to specify many indexes).",&
           usage="CUBES_LIST 1 2",&
           n_var=-1,type_of_var=integer_t,repeats=.TRUE.,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword,name="CUBES_STRIDE",&
           description="the stride used to write the cube file ",&
           usage="CUBES_STRIDE 2",&
           default_i_val=2,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL create_localize_section(subsection,error)
      CALL section_add_subsection(section, subsection, error=error)
      CALL section_release(subsection,error=error)

      CALL section_create(subsection,"PRINT","printing of information during the XAS calculation",&
           error=error,repeats=.FALSE.,required=.FALSE.)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing of basic iteration information during the XAS", &
            print_level=low_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
            error=error)
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
           enum_desc=s2a("Conjugate Gradients", "DIIS", "Steepesst Descent"),&
           enum_i_vals=(/lr_cg,lr_diis,lr_sd/),&
           error=error)
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
           description="1D line search algorithm,"//&
                      " in increasing order of robustness and cost.",&
           usage="LINESEARCH GOLD",&
           default_i_val=ot_ls_none,&
           enum_c_vals=s2a( "NONE", "2PNT", "3PNT","GOLD"),&
           enum_desc=s2a("take fixed lenght steps","extrapolate based on 2 points", &
                "... or on 3 points","perform 1D golden section search of the minimum"),&
           enum_i_vals=(/lr_ls_none,lr_ls_2pnt,lr_ls_3pnt,lr_ls_gold/),&
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

      CALL keyword_create(keyword, name="LOCALIZE_PSI0",&
           description="Maximize the localization of the ground state wavefunctions",&
           usage="LOCALIZE_PSI0",&
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

       CALL section_create(subsection,"PRINT","printing of information during the XAS calculation",&
            error=error,repeats=.FALSE.,required=.FALSE.)

       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing of basic iteration information during the XAS", &
            print_level=low_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)       

    END IF
  END SUBROUTINE create_linres_section

  SUBROUTINE create_nmr_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
            description="The chemical shift is calculated by DFPT ",&
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

       CALL keyword_create(keyword, name="CURRENT_DENSITY",&
            description="Calculate the current density and the "//&
            " induced magnetic field on a mesh in real space",&
            usage="CURRENT_DENSITY",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
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

       CALL section_create(subsection,name="PRINT",&
            description="print results of nmr calculation",&
            error=error,repeats=.FALSE.,required=.FALSE.)

         CALL cp_print_key_section_create(print_key,"CURRENT_CUBES",&
              description="Controls the printing of the current density ", &
              print_level=high_print_level,each=(/1/),add_last=(/1/),filename="",&
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
              print_level=high_print_level,each=(/1/),add_last=(/1/),filename="",&
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



    END IF

  END SUBROUTINE create_nmr_section

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_kg_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: subsection

  failure=.FALSE.

  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
       CALL section_create(section,name="kg",&
            description="The parameters specific to kg",&
            error=error)

     NULLIFY(subsection)
     CALL create_xc_fun_section(subsection,error=error)
     subsection%description="kinetic functional for the coupling of molecules"
     CALL section_add_subsection(section,subsection,error=error)
     CALL section_release(subsection,error=error)
  END IF
END SUBROUTINE create_kg_section
!***************************************************************************

END MODULE input_cp2k_dft
