!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief builds the input structure for the MIXED environment
!> \par History
!>      10.2008 created [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! **************************************************************************************************
MODULE input_cp2k_mixed
   USE bibliography,                    ONLY: Holmberg2017,&
                                              Holmberg2018,&
                                              Mavros2015,&
                                              Migliore2009
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              low_print_level,&
                                              medium_print_level
   USE input_constants,                 ONLY: mix_cdft,&
                                              mix_coupled,&
                                              mix_generic,&
                                              mix_linear_combination,&
                                              mix_minimum,&
                                              mix_restrained
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
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_mixed'

   PUBLIC :: create_mix_section

CONTAINS

! **************************************************************************************************
!> \brief Create the input section for MIXED.
!> \param section the section to create
!> \author fschiff
! **************************************************************************************************
   SUBROUTINE create_mix_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: sub2section, sub3section, subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="MIXED", &
                          description="This section contains all information to run with a hamiltonian "// &
                          "defined by a mixing of force_evals", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      NULLIFY (keyword, subsection)

      CALL keyword_create( &
         keyword, __LOCATION__, name="MIXING_TYPE", &
         description="The type of mixing to be employed", &
         usage="MIXING_TYPE LINEAR_COMBINATION", &
         default_i_val=mix_linear_combination, &
         enum_c_vals=s2a("LINEAR_COMBINATION", &
                         "MINIMUM", &
                         "COUPLED", &
                         "RESTRAINT", &
                         "GENMIX", &
                         "MIXED_CDFT"), &
         enum_desc=s2a("Linear combination of force envs (support only 2 force_evals)", &
                       "Use the force env with the minimum energy (support only 2 force_evals)", &
                       "Consider the force envs as a two state system with a given"// &
                       " coupling matrix element (support only 2 force_evals)", &
                       "Use the difference between the energy of the force envs as a"// &
                       " restraint on the first (support only 2 force_evals)", &
                       "Defines a user-driven generica coupling (support for an unlimited number of force_eval)", &
                       "Consider each force env as a CDFT state (supports an unlimited number of force_eval "// &
                       "for calculation of CDFT properties, but only two states can be mixed for forces)."), &
         enum_i_vals=(/mix_linear_combination, mix_minimum, mix_coupled, mix_restrained, mix_generic, &
                       mix_cdft/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GROUP_PARTITION", &
                          description="gives the exact number of processors for each group."// &
                          " If not specified processors allocated will be equally distributed for"// &
                          " the specified subforce_eval, trying to build a number of groups equal to the"// &
                          " number of subforce_eval specified.", &
                          usage="group_partition  2 2 4 2 4 ", type_of_var=integer_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NGROUPS", variants=(/"NGROUP"/), &
                          description="Gives the wanted number of groups. If not specified the number"// &
                          " of groups is set to the number of subforce_eval defined.", &
                          usage="ngroups 4", type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Double force_eval
      CALL section_create(subsection, __LOCATION__, name="LINEAR", &
                          description="Linear combination between two force_eval:  F= lambda F1 + (1-lambda) F2", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="LAMBDA", &
                          description="Specify the mixing parameter lambda in the formula.", &
                          usage="lambda <REAL>", type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
      ! Mixed CDFT section
      CALL section_create(subsection, __LOCATION__, name="MIXED_CDFT", &
                          description="Calculate properties involving multiple constrained states. "// &
                          "Each repetition of the FORCE_EVAL section defines a new CDFT state that is "// &
                          "included in the simulation. The DFT&QS&CDFT section must be active in each "// &
                          "FORCE_EVAL and it must be consistently defined. When the keyword "// &
                          "MIXED&NGROUPS is set to a value 2 or larger, the CDFT states are solved in "// &
                          "parallel, whereas when it is set to 1, the states are solved in serial. "// &
                          "During MD, the system can be translated using only two of the CDFT states, "// &
                          "which are selected with the keyword FORCE_STATES. The forces are determined "// &
                          "by the linear combination F= lambda F1 + (1-lambda) F2.", &
                          n_keywords=11, n_subsections=2, repeats=.FALSE., citations=(/Holmberg2017, Holmberg2018/))

      CALL keyword_create(keyword, __LOCATION__, name="LAMBDA", &
                          description="Specify the mixing parameter lambda in the formula.", &
                          usage="lambda <REAL>", type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FORCE_STATES", &
                          description="Defines the CDFT states used to translate the system. ", &
                          usage="FORCE_STATES 1 1", n_var=2, &
                          default_i_vals=(/1, 2/), type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COUPLING", &
                          description="Parameter determining how often the CDFT electronic coupling element "// &
                          "is calculated. Use a negative number to disable and 0 means every step. By default, "// &
                          "the coupling is calculated by rotating the CDFT states to eigenstates of the weight "// &
                          "function matrix when a single constraint is active and the constraint definitions are "// &
                          "identical in both CDFT states. Otherwise uses Lowdin orthogonalization. For more than "// &
                          "two CDFT states, the couplings are not computed pairwise and the values might "// &
                          "deviate from values computed separately for each unique CDFT state pair.", &
                          usage="COUPLING <INT>", &
                          default_i_val=-1, &
                          type_of_var=integer_t, n_var=1)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARALLEL_BUILD", &
                          description="Build CDFT weight function and gradients in parallel on all "// &
                          "N MPI processors before starting the CDFT SCF calculations of the 2 "// &
                          "involved CDFT states in parallel on N/2 processors. Supports only Becke "// &
                          "constraints that are identical in both states. Limited to 1 "// &
                          "charge constraint per state (different target values). "// &
                          "The keyword MIXED&NGROUPS must be set to 2.", &
                          usage="PARALLEL_BUILD TRUE", type_of_var=logical_t, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DLB", &
                          description="Controls the activation of dynamic load balancing during a mixed CDFT calculation."// &
                          " Requires Gaussian cavity confinement. Works only in conjunction with keyword PARALLEL_BUILD.", &
                          usage="DLB", type_of_var=logical_t, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="METRIC", variants=(/"COUPLING_METRIC"/), &
                          description="Compute reliability metric for the CDFT electronic coupling element by "// &
                          "diagonalizing the difference density matrix.", &
                          usage="METRIC", type_of_var=logical_t, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE., &
                          citations=(/Mavros2015/))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WFN_OVERLAP", &
                          description="Compute the CDFT electronic coupling element using the wavefunction overlap "// &
                          "method in addition to the standard method defined by the keyword COUPLING. "// &
                          "In this method, the unconstrained KS ground state wavefunction (WFN_RESTART_FILE_NAME) "// &
                          "is represented as a linear combination of the CDFT states. For more than two CDFT states, "// &
                          "the coupling is computed pairwise for every state pair (contrary to other coupling methods).", &
                          usage="WFN_OVERLAP", type_of_var=logical_t, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE., &
                          citations=(/Migliore2009/))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LOWDIN", &
                          description="Compute the CDFT electronic coupling element using Lowdin orthogonalization. "// &
                          "This is the default behavior with multiple constraints and nonidentical constraints. "// &
                          "By activating this keyword, this method is also used to compute the coupling "// &
                          "when a single constraint is active in addition to the standard method.", &
                          usage="LOWDIN", type_of_var=logical_t, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CI", variants=(/"CONFIGURATION_INTERACTION"/), &
                          description="Perform a CDFT configuration interaction calculation (CDFT-CI). "// &
                          "The CI vector is expanded in the basis of the CDFT states. Diagonalizes the "// &
                          "nonorthogonal diabatic CDFT Hamiltonian. The energies and expansion coefficients "// &
                          "of the CDFT-CI states are outputted. Keyword COUPLING must be active "// &
                          "to use this feature.", &
                          usage="LOWDIN", type_of_var=logical_t, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NONORTHOGONAL_COUPLING", &
                          variants=(/"NONORTHO_COUPLING"/), &
                          description="Print out the nonorthogonal diabatic CDFT coupling between states, "// &
                          "as it appears in the mixed CDFT Hamiltonian before orthogonalization (coupling "// &
                          "calculations) and CDFT-CI. Useful for (re)constructing the Hamiltonian for additional "// &
                          "analysis. This is the CDFT interaction energy between states.", &
                          usage="NONORTHOGONAL_COUPLING", type_of_var=logical_t, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SCALE_WITH_OCCUPATION_NUMBERS", &
                          description="Scale molecular orbitals with occupation numbers before calculating "// &
                          "the electronic coupling. Affects only simulations which employ MO smearing. "// &
                          "Disabling this keyword in conjunction with a properly selected EPS_OCCUPIED "// &
                          "threshold might be useful in systems with a large number of fractionally "// &
                          "occupied orbitals.", &
                          usage="SCALE_WITH_OCCUPATION_NUMBERS FALSE", type_of_var=logical_t, &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WFN_RESTART_FILE_NAME", &
                          description="Name of the wavefunction restart file that defines the unconstrained"// &
                          " KS ground state, which is used to compute the electronic coupling with"// &
                          " the wavefunction overlap method. May include a path.", &
                          usage="WFN_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_SVD", &
                          description="Determines the matrix inversion solver needed to compute the coupling."// &
                          " Default value implies LU decomposition, while values between 0.0 and 1.0"// &
                          " imply SVD decomposition. For SVD, the value acts as a threshold"// &
                          " for screening singular values so that only values above it are included"// &
                          " in the matrix pseudoinverse.", &
                          usage="EPS_SVD <REAL>", type_of_var=real_t, &
                          default_r_val=0.0_dp, repeats=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_OCCUPIED", &
                          description="Threshold for determining which molecular orbitals are considered occupied"// &
                          " when fractional and/or empty orbitals are employed. Can and usually should be less than"// &
                          " the threshold EPS_FERMI_DIRAC defined in section SCF&SMEAR. Note that the number occupied"// &
                          " MOs should be constant in each CDFT state, since the CDFT coupling is only defined between"// &
                          " states in the same spin state. Fractionally occupied MOs might exhibit linear dependencies"// &
                          " and a singular value decomposition (EPS_SVD) can be used for removing these.", &
                          usage="EPS_OCCUPIED <REAL>", type_of_var=real_t, &
                          default_r_val=1.0E-6_dp, repeats=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LOAD_SCALE", &
                          description="Control parameter for dynamic load balancing during a mixed CDFT calculation."// &
                          " See code for details. Works only in conjunction with keyword PARALLEL_BUILD.", &
                          usage="LOAD_SCALE <REAL>", type_of_var=real_t, &
                          default_r_val=2.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MORE_WORK", &
                          description="Control parameter for dynamic load balancing during a mixed CDFT calculation."// &
                          " See code for details. Works only in conjunction with keyword PARALLEL_BUILD.", &
                          usage="MORE_WORK <INT>", type_of_var=integer_t, &
                          default_i_val=0, repeats=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VERY_OVERLOADED", &
                          description="Control parameter for dynamic load balancing during a mixed CDFT calculation."// &
                          " See code for details. Works only in conjunction with keyword PARALLEL_BUILD.", &
                          usage="VERY_OVERLOADED <REAL>", type_of_var=real_t, &
                          default_r_val=0.0_dp, repeats=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BLOCK_DIAGONALIZE", &
                          description="Block diagonalize the CDFT Hamiltonian. Control settings should be given in "// &
                          "section &BLOCK_DIAGONALIZE. All requested electronic couplings are printed out after "// &
                          "block diagonalization. When CDFT-CI and block diagonalization are both requested, "// &
                          "the CI calculation is performed using the block diagonalized Hamiltonian.", &
                          usage="BLOCK_DIAGONALIZE", type_of_var=logical_t, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      NULLIFY (sub2section)
      CALL create_mixed_cdft_block_section(sub2section)
      CALL section_add_subsection(subsection, sub2section)
      CALL section_release(sub2section)

      CALL create_print_mixed_cdft_section(sub2section)
      CALL section_add_subsection(subsection, sub2section)
      CALL section_release(sub2section)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
      !
      CALL section_create(subsection, __LOCATION__, name="COUPLING", &
                          description="Coupling between two force_eval: E=(E1+E2 - sqrt((E1-E2)**2+4*H12**2))/2", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="COUPLING_PARAMETER", &
                          description="Coupling parameter H12 used in the coupling", &
                          usage="COUPLING_PARAMETER <REAL>", type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="RESTRAINT", &
                          description="Restraint between two force_eval: E = E1 + k*(E1-E2-t)**2", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="RESTRAINT_TARGET", &
                          description="Target value of the restraint (t) ", &
                          usage="RESTRAINT_TARGET <REAL>", type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTRAINT_STRENGTH", &
                          description="Strength of the restraint (k) in "// &
                          "k*(E1-E2-t)**2", &
                          usage="RESTRAINT_STRENGTH <REAL>", type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! Multiple force_eval
      CALL section_create(subsection, __LOCATION__, name="GENERIC", &
                          description="User driven coupling between two or more force_eval.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="MIXING_FUNCTION", &
                          description="Specifies the mixing functional form in mathematical notation.", &
                          usage="MIXING_FUNCTION (E1+E2-LOG(E1/E2))", type_of_var=lchar_t, &
                          n_var=1)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VARIABLES", &
                          description="Defines the variables of the functional form. To allow an efficient"// &
                          " mapping the order of the energy variables will be considered identical to the"// &
                          " order of the force_eval in the force_eval_order list.", &
                          usage="VARIABLES x", type_of_var=char_t, &
                          n_var=-1)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARAMETERS", &
                          description="Defines the parameters of the functional form", &
                          usage="PARAMETERS a b D", type_of_var=char_t, &
                          n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VALUES", &
                          description="Defines the values of parameter of the functional form", &
                          usage="VALUES ", type_of_var=real_t, &
                          n_var=-1, repeats=.TRUE., unit_str="internal_cp2k")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="UNITS", &
                          description="Optionally, allows to define valid CP2K unit strings for each parameter value. "// &
                          "It is assumed that the corresponding parameter value is specified in this unit.", &
                          usage="UNITS angstrom eV*angstrom^-1 angstrom^1 K", type_of_var=char_t, &
                          n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DX", &
                          description="Parameter used for computing the derivative with the Ridders' method.", &
                          usage="DX <REAL>", default_r_val=0.1_dp, unit_str="bohr")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ERROR_LIMIT", &
                          description="Checks that the error in computing the derivative is not larger than "// &
                          "the value set; in case error is larger a warning message is printed.", &
                          usage="ERROR_LIMIT <REAL>", default_r_val=1.0E-12_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! Mapping of atoms
      NULLIFY (sub2section, sub3section)
      CALL section_create(subsection, __LOCATION__, name="MAPPING", &
                          description="Defines the mapping of atoms for the different force_eval with the mixed force_eval."// &
                          " The default is to have a mapping 1-1 between atom index (i.e. all force_eval share the same"// &
                          " geometrical structure). The mapping is based on defining fragments and the mapping the "// &
                          " fragments between the several force_eval and the mixed force_eval", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      ! Mixed force_eval
      CALL section_create(sub2section, __LOCATION__, name="FORCE_EVAL_MIXED", &
                          description="Defines the fragments for the mixed force_eval (reference)", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      CALL section_create(sub3section, __LOCATION__, name="FRAGMENT", &
                          description="Fragment definition", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Defines the index of the fragment defined", &
                          usage="<INTEGER>", type_of_var=integer_t, n_var=1)
      CALL section_add_keyword(sub3section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_DEFAULT_KEYWORD_", &
                          description="Starting and ending atomic index defining one fragment must be provided", &
                          usage="{Integer} {Integer}", type_of_var=integer_t, n_var=2, repeats=.TRUE.)
      CALL section_add_keyword(sub3section, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(sub2section, sub3section)
      CALL section_release(sub3section)
      CALL section_add_subsection(subsection, sub2section)
      CALL section_release(sub2section)

      ! All other force_eval
      CALL section_create(sub2section, __LOCATION__, name="FORCE_EVAL", &
                          description="Defines the fragments and the mapping for each force_eval (an integer index (ID) "// &
                          "needs to be provided as parameter)", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create( &
         keyword, __LOCATION__, name="DEFINE_FRAGMENTS", &
         description="Specify the fragments definition of the force_eval through the fragments of the"// &
         " force_eval_mixed. This avoids the pedantic definition of the fragments for the force_eval,"// &
         " assuming the order of the fragments for the specified force_eval is the same as the sequence "// &
         " of integers provided. Easier to USE should be preferred to the specification of the single fragments.", &
         usage="DEFINE_FRAGMENTS <INTEGER> .. <INTEGER>", type_of_var=integer_t, n_var=-1)
      CALL section_add_keyword(sub2section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Defines the index of the force_eval for which fragments and mappings are provided", &
                          usage="<INTEGER>", type_of_var=integer_t, n_var=1)
      CALL section_add_keyword(sub2section, keyword)
      CALL keyword_release(keyword)

      CALL section_create(sub3section, __LOCATION__, name="FRAGMENT", &
                          description="Fragment definition", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Defines the index of the fragment defined", &
                          usage="<INTEGER>", type_of_var=integer_t, n_var=1)
      CALL section_add_keyword(sub3section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_DEFAULT_KEYWORD_", &
                          description="Starting and ending atomic index defining one fragment must be provided", &
                          usage="{Integer} {Integer}", type_of_var=integer_t, n_var=2, repeats=.FALSE.)
      CALL section_add_keyword(sub3section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAP", &
                          description="Provides the index of the fragment of the MIXED force_eval mapped on the"// &
                          " locally defined fragment.", &
                          usage="MAP <INTEGER>", type_of_var=integer_t, n_var=1, repeats=.FALSE.)
      CALL section_add_keyword(sub3section, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(sub2section, sub3section)
      CALL section_release(sub3section)
      CALL section_add_subsection(subsection, sub2section)
      CALL section_release(sub2section)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_print_mix_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
   END SUBROUTINE create_mix_section

! **************************************************************************************************
!> \brief Create the print section for mixed
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_print_mix_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="print", &
                          description="Section of possible print options in MIXED env.", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

      NULLIFY (print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "PROGRAM_RUN_INFO", &
                                       description="Controls the printing of information during the evaluation of "// &
                                       "the mixed environment. ", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "DIPOLE", &
                                       description="Controls the printing of dipole information. "// &
                                       "Requires the DIPOLE calculation be active for all subforce_eval.", &
                                       print_level=medium_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)
   END SUBROUTINE create_print_mix_section

! **************************************************************************************************
!> \brief Create the print section specific to mixed CDFT (forked from print_mix_section)
!> \param section the section to create
!> \author Nico Holmberg [06.2017]
! **************************************************************************************************
   SUBROUTINE create_print_mixed_cdft_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="print", &
                          description="Section of possible print options for the mixed CDFT environment.", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

      NULLIFY (print_key, keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "PROGRAM_RUN_INFO", &
                                       description="Controls the printing of information during the evaluation of "// &
                                       "the mixed CDFT environment. ", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")

      CALL keyword_create(keyword, __LOCATION__, name="MO_OVERLAP_MATRIX", &
                          description="Controls the printing of the MO overlap matrices between CDFT states. "// &
                          "The matrices are printed out in plain text.", &
                          usage="MO_OVERLAP_MATRIX TRUE", type_of_var=logical_t, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MO_OVERLAP_EIGENVALUES", &
                          description="Controls the printing of the eigenvalues/singular values of the CDFT MO overlap "// &
                          "matrices. The product of the eigenvalues/singular values is the CDFT MO overlap. "// &
                          "Useful mainly for checking which singular values will get screened for a particular EPS_SVD.", &
                          usage="MO_OVERLAP_EIGENVALUES TRUE", type_of_var=logical_t, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_print_mixed_cdft_section
! **************************************************************************************************
!> \brief Creates the control section used to setup block diagonalization of the mixed
!>        CDFT Hamiltonian matrix
!> \param section the section to create
!> \author Nico Holmberg [11.2017]
! **************************************************************************************************
   SUBROUTINE create_mixed_cdft_block_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="BLOCK_DIAGONALIZE", &
                          description="Control section to setup block diagonalization of the mixed CDFT Hamiltonian. "// &
                          "Constructs a new Hamiltonian by diagonalizing the initial matrix within each block and "// &
                          "by rotating the off-diagonal blocks (which represent the interactions between different "// &
                          "blocks) by the eigenvectors of diagonal blocks.", &
                          n_keywords=2, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="BLOCK", &
                          description="Defines which CDFT states are included in a block. Each repetition of this keyword "// &
                          "defines a new block. The Hamiltonian matrix elements of the requested states are collected "// &
                          "into a new matrix and subsequently diagonalized. The eigenvectors of this matrix are used to "// &
                          "rotate the matrix blocks describing the interactions between blocks.", &
                          usage="BLOCK 1 2", repeats=.TRUE., &
                          type_of_var=integer_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="IGNORE_EXCITED", &
                          description="Ignore excited states related to each block when constructing the new mixed "// &
                          "CDFT Hamiltonian. This reduces the dimensionality of the Hamiltonian.", &
                          usage="IGNORE_EXCITED FALSE", type_of_var=logical_t, &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RECURSIVE_DIAGONALIZATION", &
                          description="Perform block diagonalization recursively until only two blocks remain. "// &
                          "For example, if the elements of a 8x8 matrix are first collected into 4 blocks "// &
                          "(using keyword BLOCK), this keyword will transform the matrix to a 2x2 matrix "// &
                          "(8x8 -> 4x4 -> 2x2). In this example, the blocks of the 2x2 matrix would be "// &
                          "assembled from the first and last 2 blocks of the 4x4 matrix.", &
                          usage="RECURSIVE_DIAGONALIZATION TRUE", type_of_var=logical_t, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_mixed_cdft_block_section

END MODULE input_cp2k_mixed
