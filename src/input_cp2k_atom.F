!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief builds the input structure for the ATOM module
!> \author jgh
! *****************************************************************************
MODULE input_cp2k_atom
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE f77_blas
  USE input_constants
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
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_atom'

  PUBLIC :: create_atom_section, create_powell_section

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief Creates the input section for the atom code
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author jgh
! *****************************************************************************
  SUBROUTINE create_atom_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_atom_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="ATOM",&
            description="Section handling input for atomic calculations.",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            supported_feature=.FALSE.,error=error)
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMIC_NUMBER",&
            description="Specify the atomic number",&
            default_i_val=1, supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ELEMENT",&
            description="Specify the element to be calculated",&
            usage="ELEMENT char",n_var=1,type_of_var=char_t,&
            default_c_val="H", supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RUN_TYPE",&
            description="Type of run that you want to perform "//&
                "[ENERGY,BASIS_OPTIMIZATION,PSEUDOPOTENTIAL_OPTIMIZATION,,...] ",&
            usage="RUN_TYPE ENERGY",&
            default_i_val=atom_energy_run,&
            enum_c_vals=s2a("NONE","ENERGY","BASIS_OPTIMIZATION","PSEUDOPOTENTIAL_OPTIMIZATION"),&
            enum_i_vals=(/atom_no_run, atom_energy_run, atom_basis_run, atom_pseudo_run/),&
            supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COULOMB_INTEGRALS",&
            description="Method to calculate Coulomb integrals",&
            usage="COULOMB_INTEGRALS analytic",&
            default_i_val=do_numeric,&
            enum_c_vals=(/"ANALYTIC                 ",&
                          "SEMI_ANALYTIC            ",&
                          "NUMERIC                  "/),&
            enum_i_vals= (/ do_analytic, do_semi_analytic, do_numeric /),&
            supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXCHANGE_INTEGRALS",&
            description="Method to calculate Exchange integrals",&
            usage="EXCHANGE_INTEGRALS analytic",&
            default_i_val=do_numeric,&
            enum_c_vals=(/"ANALYTIC                 ",&
                          "SEMI_ANALYTIC            ",&
                          "NUMERIC                  "/),&
            enum_i_vals= (/ do_analytic, do_semi_analytic, do_numeric /),&
            supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CORE",&
            description="Specifies the core electrons for a pseudopotential",&
            usage="CORE 1s2 ... or CORE [Ne] or CORE none for 0 electron cores", required=.TRUE., repeats=.FALSE.,&
            n_var=-1,type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ELECTRON_CONFIGURATION",&
            description="Specifies the electron configuration. "//&
            "Optional the multiplicity (m) and a core state [XX] can be declared",&
            usage="ELECTRON_CONFIGURATION (1) [Ne] 3s2 ... ", required=.TRUE., repeats=.TRUE.,&
            n_var=-1,type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_ANGULAR_MOMENTUM",&
            description="Specifies the largest angular momentum calculated [0-3]",&
            usage="MAX_ANGULAR_MOMENTUM 3", required=.FALSE., repeats=.FALSE.,&
            default_i_val=3, supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CALCULATE_STATES",&
            description="Specifies the number of states calculated per l value",&
            usage="CALCULATE_STATES  5 5 5 3 ", required=.FALSE., repeats=.FALSE.,&
            default_i_val=0, n_var=-1,type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_atom_print_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_atom_aebasis_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_atom_ppbasis_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_atom_method_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_optimization_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_potential_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_powell_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_atom_section

! *****************************************************************************
!> \brief Create the print atom section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author jgh
! *****************************************************************************
  SUBROUTINE create_atom_print_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_atom_print_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="print",&
            description="Section of possible print options specific of the ATOM code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(print_key, keyword)

       CALL cp_print_key_section_create(print_key,"PROGRAM_BANNER",&
            description="Controls the printing of the banner of the ATOM program",&
            print_level=silent_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"METHOD_INFO",&
            description="Controls the printing of method information",&
            print_level=medium_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"BASIS_SET",&
            description="Controls the printing of the basis sets",&
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"POTENTIAL",&
            description="Controls the printing of the potentials",&
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"FIT_DENSITY",&
            description="Fit the total electronic density to a linear combination of Gaussian functions",&
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL keyword_create(keyword, name="NUM_GTO",&
            description="Number of Gaussian type functions for density fit",&
            usage="NUM_GTO integer ",type_of_var=integer_t,&
            default_i_val=40,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"FIT_KGPOT",&
            description="Fit an approximation to the non-additive kinetic energy potential used in KG",&
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL keyword_create(keyword, name="NUM_GAUSSIAN",&
            description="Number of Gaussian terms for the fit",&
            usage="NUM_GAUSSIAN integer ",type_of_var=integer_t,&
            default_i_val=1,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="NUM_POLYNOM",&
            description="Number of terms in the polynomial expansion",&
            usage="NUM_POLYNOM integer ",type_of_var=integer_t,&
            default_i_val=4,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"RESPONSE_BASIS",&
            description="Calculate a response basis set contraction scheme",&
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL keyword_create(keyword, name="DELTA_CHARGE",&
            description="Variation of charge used in finite difference calculation",&
            usage="DELTA_CHARGE real ",type_of_var=real_t,&
            default_r_val=0.05_dp,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="DERIVATIVES",&
            description="Number of wavefunction derivatives to calculate",&
            usage="DERIVATIVES integer ",type_of_var=integer_t,&
            default_i_val=2,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"SCF_INFO",&
            description="Controls the printing of SCF information",&
            print_level=medium_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"ORBITALS",&
            description="Controls the printing of the optimized orbitals information",&
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"FIT_PSEUDO",&
            description="Controls the printing of FIT PSEUDO task",&
            print_level=medium_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"FIT_BASIS",&
            description="Controls the printing of FIT BASIS task",&
            print_level=medium_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"UPF_FILE",&
            description="Write GTH pseudopotential in UPF format",&
            print_level=high_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"ENERGIES_MINUS_KINETIC",&
            description="Print out the total energy and orbital energies without "  //&
            "the kinetic energy component. Useful for atomic calculations used" //&
            "during SCPTB parametrization",&
            print_level=debug_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
    END IF
  END SUBROUTINE create_atom_print_section

! *****************************************************************************
!> \brief Create the all-electron basis section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author jgh
! *****************************************************************************
  SUBROUTINE create_atom_aebasis_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_atom_aebasis_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="AE_BASIS",&
            description="Section of basis set information for all-electron calculations.",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL atom_basis_section(section,error)

    END IF
  END SUBROUTINE create_atom_aebasis_section

! *****************************************************************************
!> \brief Create the pseudopotential basis section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author jgh
! *****************************************************************************
  SUBROUTINE create_atom_ppbasis_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_atom_ppbasis_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="PP_BASIS",&
            description="Section of basis set information for pseudopotential calculations.",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL atom_basis_section(section,error)

    END IF
  END SUBROUTINE create_atom_ppbasis_section

! *****************************************************************************
!> \brief Keywords in the atom basis section
!> \param section the section to fill
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author jgh
! *****************************************************************************
  SUBROUTINE atom_basis_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_basis_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="BASIS_TYPE",&
            description="Basis set type",&
            usage="BASIS_TYPE GAUSSIAN",&
            default_i_val=gaussian,&
            enum_c_vals=(/"GAUSSIAN                 ",&
                          "GEOMETRICAL_GTO          ",&
                          "CONTRACTED_GTO           ",&
                          "SLATER                   ",&
                          "NUMERICAL                "/),&
            enum_i_vals= (/ gaussian, geometrical_gto, contracted_gto, slater, numerical /),&
            supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NUM_GTO",&
            description="Number of Gaussian type functions for s, p, d, ...",&
            usage="NUM_GTO 5 5 5 ",n_var=-1,type_of_var=integer_t,&
            default_i_val=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NUM_SLATER",&
            description="Number of Slater type functions for s, p, d, ...",&
            usage="NUM_SLATER 5 5 5 ",n_var=-1,type_of_var=integer_t,&
            default_i_val=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="START_INDEX",&
            description="Starting index for Geometrical Basis sets",&
            usage="START_INDEX 0 2 5 4 ",n_var=-1,type_of_var=integer_t,&
            default_i_val=0,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="S_EXPONENTS",&
            description="Exponents for s functions",&
            usage="S_EXPONENTS 1.0 2.0 ... ",n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="P_EXPONENTS",&
            description="Exponents for p functions",&
            usage="P_EXPONENTS 1.0 2.0 ... ",n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="D_EXPONENTS",&
            description="Exponents for d functions",&
            usage="D_EXPONENTS 1.0 2.0 ... ",n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="F_EXPONENTS",&
            description="Exponents for f functions",&
            usage="F_EXPONENTS 1.0 2.0 ... ",n_var=-1,type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="S_QUANTUM_NUMBERS",&
            description="Main quantum numbers for s functions",&
            usage="S_QUANTUM_NUMBERS 1 2 ... ",n_var=-1,type_of_var=integer_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="P_QUANTUM_NUMBERS",&
            description="Main quantum numbers for p functions",&
            usage="P_QUANTUM_NUMBERS 2 3 ... ",n_var=-1,type_of_var=integer_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="D_QUANTUM_NUMBERS",&
            description="Main quantum numbers for d functions",&
            usage="D_QUANTUM_NUMBERS 3 4 ... ",n_var=-1,type_of_var=integer_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="F_QUANTUM_NUMBERS",&
            description="Main quantum numbers for f functions",&
            usage="F_QUANTUM_NUMBERS 4 5 ... ",n_var=-1,type_of_var=integer_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GEOMETRICAL_FACTOR",&
            description="Geometrical basis: factor C in a*C^k",&
            usage="GEOMETRICAL_FACTOR real",&
            default_r_val=2.6_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GEO_START_VALUE",&
            description="Geometrical basis: starting value a in a*C^k",&
            usage="GEO_START_VALUE real",&
            default_r_val=0.016_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BASIS_SET_FILE_NAME",&
            description="Name of the basis set file, may include a path",&
            usage="BASIS_SET_FILE_NAME <FILENAME>",&
            default_lc_val="BASIS_SET",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BASIS_SET",&
            variants=s2a("ORBITAL_BASIS_SET","ORB_BASIS"),&
            description="The contracted Gaussian basis set",&
            usage="BASIS_SET DZVP", default_c_val=" ", &
            n_var=1, required=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="QUADRATURE",&
            description="Algorithm to construct the atomic radial grids",&
            usage="QUADRATURE GC_SIMPLE",&
            enum_c_vals=s2a("GC_SIMPLE","GC_TRANSFORMED","GC_LOG"),&
            enum_i_vals=(/ do_gapw_gcs,do_gapw_gct,do_gapw_log/),&
            default_i_val=do_gapw_log, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GRID_POINTS",&
            description="Number of radial grid points",&
            usage="GRID_POINTS integer",&
            default_i_val=400,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_EIGENVALUE",&
            description="Cutoff of overlap matrix eigenvalues included into basis",&
            usage="EPS_EIGENVALUE real",&
            default_r_val=1.e-12_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsection)
       CALL create_basis_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE atom_basis_section

! *****************************************************************************
!> \brief Create the method section for Atom calculations
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author jgh
! *****************************************************************************
  SUBROUTINE create_atom_method_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_atom_method_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    NULLIFY(subsection,keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="METHOD",&
            description="Section of information on method to use.",&
            n_keywords=0, n_subsections=2, repeats=.TRUE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword, name="METHOD_TYPE",&
            description="Type of electronic structure method to be used",&
            usage="METHOD_TYPE KOHN-SHAM",&
            default_i_val=do_rks_atom,&
            enum_c_vals=(/"KOHN-SHAM                ",&
                          "RKS                      ",&
                          "UKS                      ",&
                          "HARTREE-FOCK             ",&
                          "RHF                      ",&
                          "UHF                      ",&
                          "ROHF                     "/),&
            enum_i_vals= (/ do_rks_atom, do_rks_atom, do_uks_atom, do_rhf_atom, &
                            do_rhf_atom, do_uhf_atom, do_rohf_atom /),&
            supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RELATIVISTIC",&
            description="Type of scalar relativistic method to be used",&
            usage="RELATIVISTIC DKH(2)",&
            default_i_val=do_nonrel_atom,&
            enum_c_vals=(/"OFF                      ",&
                          "ZORA(MP)                 ",&
                          "scZORA(MP)               ",&
                          "DKH(0)                   ",&
                          "DKH(1)                   ",&
                          "DKH(2)                   ",&
                          "DKH(3)                   ",&
                          "DKH(4)                   ",&
                          "DKH(5)                   "/),&
            enum_i_vals= (/ do_nonrel_atom, do_zoramp_atom, do_sczoramp_atom, do_dkh0_atom,&
            do_dkh1_atom, do_dkh2_atom, do_dkh3_atom, do_dkh4_atom, do_dkh5_atom /),&
            supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_xc_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

! ZMP creating zubsection for the zmp calculations
       CALL create_zmp_section(subsection,error) 
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
                          
       CALL create_external_vxc(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_atom_method_section

! *****************************************************************************
!> \brief Create the ZMP subsection for Atom calculations
!>  
!> \author D. Varsano [daniele.varsano@nano.cnr.it] 
! *****************************************************************************
  SUBROUTINE create_zmp_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_zmp_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    NULLIFY(subsection,keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="ZMP",&
             description="Section used to specify ZMP Potentials.",&
             n_keywords=3, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
             error=error)

       CALL keyword_create(keyword, name="FILE_DENSITY",&
            description="Specifies the filename containing the target density ",&
            usage="FILE_DENSITY <FILENAME>",&
            type_of_var=char_t,default_c_val="RHO_O.dat", n_var=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GRID_TOL",&
            description="Tolerance in the equivalence of read-grid in ZMP method",&
            usage="GRID_TOL <REAL>", default_r_val=1.E-12_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LAMBDA",&
            description="Parameter used for the constraint in ZMP method",&
            usage="LAMBDA <REAL>", default_r_val=10.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DM",&
            description="read external density from density matrix",&
            usage="DM <LOGICAL>", type_of_var=logical_t,default_l_val=.FALSE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_zmp_restart_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_zmp_section

! *****************************************************************************
!> \brief Create the ZMP restart subsection for Atom calculations
!>  
!> \author D. Varsano [daniele.varsano@nano.cnr.it] 
! *****************************************************************************
  SUBROUTINE create_zmp_restart_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_zmp_restart_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="RESTART",&
             description="Section used to specify the restart option in the ZMP"//&
                         "procedure, and the file that must be read.",&
             n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
             error=error)

       CALL keyword_create(keyword, name="FILE_RESTART",&
            description="Specifies the filename containing the restart file density ",&
            usage="FILE_RESTART <FILENAME>",&
            type_of_var=char_t,default_c_val="RESTART.wfn", n_var=-1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF

  END SUBROUTINE create_zmp_restart_section

! *****************************************************************************
!> \brief Subroutine to create the external v_xc potential 
!>  
!> \author D. Varsano [daniele.varsano@nano.cnr.it] 
! *****************************************************************************
  SUBROUTINE create_external_vxc(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_external_vxc', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CALL section_create(section,name="EXTERNAL_VXC",&
            description="Section used to specify exernal VXC Potentials.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

      CALL keyword_create(keyword, name="FILE_VXC",&
            description="Specifies the filename containing the external vxc ",&
            usage="FILE_VXC <FILENAME>",&
            type_of_var=char_t,default_c_val="VXC.dat", n_var=-1, error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="GRID_TOL",&
            description="Tolerance in the equivalence of read-grid in ZMP method",&
            usage="GRID_TOL <REAL>", default_r_val=1.E-12_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_external_vxc

! *****************************************************************************
!> \brief Create the optimization section for Atom calculations
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author jgh
! *****************************************************************************
  SUBROUTINE create_optimization_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_optimization_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="OPTIMIZATION",&
            description="Section of information on optimization thresholds and algorithms.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="MAX_ITER",&
            description="Maximum number of iterations for optimization",&
            usage="MAX_ITER  50", default_i_val=200,&
            supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_SCF",&
            description="Convergence criterion for SCF",&
            usage="EPS_SCF 1.e-10", default_r_val=1.e-6_dp,&
            supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DAMPING",&
            description="Damping parameter for extrapolation method",&
            usage="DAMPING 0.4", default_r_val=0.4_dp,&
            supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_DIIS",&
            description="Starting DIIS method at convergence to EPS_DIIS",&
            usage="EPS_DIIS  0.01", default_r_val=10000._dp,&
            supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="N_DIIS",&
            description="Maximum number of DIIS vectors",&
            usage="N_DIIS  6", default_i_val=5,&
            supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_optimization_section

! *****************************************************************************
!> \brief Create the potential section for Atom calculations
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author jgh
! *****************************************************************************
  SUBROUTINE create_potential_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_potential_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="POTENTIAL",&
            description="Section of information on potential.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.FALSE.,&
            error=error)

       CALL keyword_create(keyword, name="CONFINEMENT",&
            description="Definition of confinement potential",&
            usage="CONFINEMENT prefactor range exponent", &
            default_r_vals=(/ 0._dp, 4._dp, 2._dp /),&
            repeats=.FALSE., n_var=-1,&
            supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PSEUDO_TYPE",&
            description="Pseudopotential type",&
            usage="PSEUDO_TYPE GTH",&
            default_i_val=no_pseudo,&
            enum_c_vals=(/"NONE                     ",&
                          "GTH                      "/),&
            enum_i_vals= (/ no_pseudo, gth_pseudo /),&
            supported_feature=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="POTENTIAL_FILE_NAME",&
            description="Name of the pseudo potential file, may include a path",&
            usage="POTENTIAL_FILE_NAME <FILENAME>",&
            default_lc_val="POTENTIAL",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="POTENTIAL_NAME",&
            variants=(/ "POT_NAME" /),&
            description="The name of the pseudopotential for the defined kind.",&
            usage="POTENTIAL_NAME <PSEUDO-POTENTIAL-NAME>", default_c_val=" ", n_var=1,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(subsection)
       CALL create_gthpotential_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_potential_section

! *****************************************************************************
!> \brief Creates the &GTH_POTENTIAL section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_gthpotential_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_gthpotential_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    IF (.NOT. failure) THEN
       CALL section_create(section,name="GTH_POTENTIAL",&
            description="Section used to specify Potentials.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword)
       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="CP2K Pseudo Potential Standard Format (GTH, ALL or KG)",&
            repeats=.TRUE.,type_of_var=lchar_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_gthpotential_section

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
!> \brief Creates the &POWELL section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_powell_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_powell_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    IF (.NOT. failure) THEN
       CALL section_create(section,name="powell",&
            description="Section defines basic parameters for Powell optimization",&
            n_keywords=4, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="ACCURACY",&
            description="Final accuracy requested in optimization (RHOEND)",&
            usage="ACCURACY 0.00001",&
            default_r_val=1.e-6_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STEP_SIZE",&
            description="Initial step size for search algorithm (RHOBEG)",&
            usage="STEP_SIZE 0.005",&
            default_r_val=0.005_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MAX_FUN",&
            description="Maximum number of function evaluations",&
            usage="MAX_FUN 1000",&
            default_i_val=5000,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WEIGHT_POT_VIRTUAL",&
            description="Weight for virtual states in pseudopotential optimization",&
            usage="WEIGHT_POT_VIRTUAL 1.0",&
            default_r_val=0.01_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WEIGHT_POT_SEMICORE",&
            description="Weight for semi core states in pseudopotential optimization",&
            usage="WEIGHT_POT_SEMICORE 1.0",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WEIGHT_POT_VALENCE",&
            description="Weight for valence states in pseudopotential optimization",&
            usage="WEIGHT_POT_VALENCE 1.0",&
            default_r_val=100.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WEIGHT_POT_NODE",&
            description="Weight for node mismatch in pseudopotential optimization",&
            usage="WEIGHT_POT_NODE 1.0",&
            default_r_val=1000.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WEIGHT_ELECTRON_CONFIGURATION",&
            description="Weight for different electronic states in optimization",&
            usage="WEIGHT_ELECTRON_CONFIGURATION 1.0 0.1 ...",&
            n_var=-1,type_of_var=real_t, default_r_val=1.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WEIGHT_METHOD",&
            description="Weight for different methods in optimization",&
            usage="WEIGHT_METHOD 1.0 0.1 ...",&
            n_var=-1,type_of_var=real_t, default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TARGET_POT_VIRTUAL",&
            description="Target accuracy for virtual state eigenvalues in pseudopotential optimization",&
            usage="TARGET_POT_VIRTUAL 0.0001",&
            default_r_val=0.001_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TARGET_POT_VALENCE",&
            description="Target accuracy for valence state eigenvalues in pseudopotential optimization",&
            usage="TARGET_POT_VALENCE 0.0001",&
            default_r_val=0.00001_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="TARGET_POT_SEMICORE",&
            description="Target accuracy for semicore state eigenvalues in pseudopotential optimization",&
            usage="TARGET_POT_SEMICORE 0.01",&
            default_r_val=0.001_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WEIGHT_PSIR0",&
            description="Weight for the wavefunctions at r=0 (only occupied states)",&
            usage="WEIGHT_PSIR0 0.01",&
            default_r_val=0._dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCOV_MULTIPLICATION",&
            description="Multiply Rcov integration limit for charge conservation",&
            usage="RCOV_MULTIPLICATION  1.10",&
            default_r_val=1._dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_powell_section

! *****************************************************************************

END MODULE input_cp2k_atom
