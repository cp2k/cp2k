!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the semi empirical section of the input
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!>      07.2024 moved out of input_cp2k [JGH]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_se
   USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                              high_print_level
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE input_constants,                 ONLY: &
        do_se_is_kdso, do_se_is_kdso_d, do_se_is_slater, do_se_lr_ewald, do_se_lr_ewald_gks, &
        do_se_lr_ewald_r3, do_se_lr_none, gaussian, numerical, slater
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
                                              real_t
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_se'

   PUBLIC :: create_se_control_section

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_se_control_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="SE", &
                          description="Parameters needed to set up the Semi-empirical methods", &
                          n_keywords=8, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="ORTHOGONAL_BASIS", &
                          description="Assume orthogonal basis set. This flag is overwritten by "// &
                          "methods with fixed orthogonal/non-orthogonal basis set.", &
                          usage="ORTHOGONAL_BASIS", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STO_NG", &
                          description="Provides the order of the Slater orbital expansion of Gaussian-Type Orbitals.", &
                          usage="STO_NG", default_i_val=6)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ANALYTICAL_GRADIENTS", &
                          description="Nuclear Gradients are computed analytically or numerically", &
                          usage="ANALYTICAL_GRADIENTS", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DELTA", &
                          description="Step size in finite difference force calculation", &
                          usage="DELTA {real} ", default_r_val=1.e-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INTEGRAL_SCREENING", &
                          description="Specifies the functional form for the ", &
                          usage="INTEGRAL_SCREENING (KDSO|KDSO-D|SLATER)", &
                          enum_c_vals=s2a("KDSO", "KDSO-D", "SLATER"), &
                          enum_i_vals=(/do_se_IS_kdso, do_se_IS_kdso_d, do_se_IS_slater/), &
                          enum_desc=s2a("Uses the standard NDDO Klopman-Dewar-Sabelli-Ohno equation "// &
                                        "for the screening of the Coulomb interactions.", &
                                        "Uses a modified Klopman-Dewar-Sabelli-Ohno equation, dumping the screening "// &
                                        "parameter for the Coulomb interactions.", &
                                        "Uses an exponential Slater-type function for modelling the Coulomb interactions."), &
                          default_i_val=do_se_IS_kdso)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PERIODIC", &
                          description="Specifies the type of treatment for the electrostatic long-range part "// &
                          "in semi-empirical calculations.", &
                          usage="PERIODIC (NONE|EWALD|EWALD_R3|EWALD_GKS)", &
                          enum_c_vals=s2a("NONE", "EWALD", "EWALD_R3", "EWALD_GKS"), &
                          enum_i_vals=(/do_se_lr_none, do_se_lr_ewald, do_se_lr_ewald_r3, do_se_lr_ewald_gks/), &
                          enum_desc=s2a("The long-range part is not explicitly treaten. The interaction "// &
                                        "depends uniquely on the Cutoffs used for the calculation.", &
                                        "Enables the usage of Multipoles Ewald summation schemes. The short-range part "// &
                                        "is tapered according the value of RC_COULOMB.", &
                                        "Enables the usage of Multipoles Ewald summation schemes together with a long-range "// &
                                        "treatment for the 1/R^3 term, deriving from the short-range component. This option "// &
                                        "is active only for K-DSO integral screening type.", &
                                        "Use Ewald directly in Coulomb integral evaluation, works only with Slater screening"), &
                          default_i_val=do_se_lr_none)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FORCE_KDSO-D_EXCHANGE", &
                          description="This keywords forces the usage of the KDSO-D integral screening "// &
                          "for the Exchange integrals (default is to apply the screening only to the "// &
                          "Coulomb integrals.", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISPERSION", &
                          description="Use dispersion correction", &
                          lone_keyword_l_val=.TRUE., &
                          usage="DISPERSION", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISPERSION_PARAMETER_FILE", &
                          description="Specify file that contains the atomic dispersion parameters", &
                          usage="DISPERSION_PARAMETER_FILE filename", &
                          n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISPERSION_RADIUS", &
                          description="Define radius of dispersion interaction", &
                          usage="DISPERSION_RADIUS", default_r_val=15._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COORDINATION_CUTOFF", &
                          description="Define cutoff for coordination number calculation", &
                          usage="COORDINATION_CUTOFF", default_r_val=1.e-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="D3_SCALING", &
                          description="Scaling parameters (s6,sr6,s8) for the D3 dispersion method,", &
                          usage="D3_SCALING 1.0 1.0 1.0", n_var=3, default_r_vals=(/0.0_dp, 0.0_dp, 0.0_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL create_coulomb_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_exchange_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_screening_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_lr_corr_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_neighbor_lists_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_se_memory_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_se_print_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_se_ga_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_se_control_section

! **************************************************************************************************
!> \brief Create the COULOMB se section
!> \param section the section to create
!> \date  03.2009
!> \author Teodoro Laino [tlaino]
! **************************************************************************************************
   SUBROUTINE create_lr_corr_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="LR_CORRECTION", &
                          description="Setup parameters for the evaluation of the long-range correction term in SE "// &
                          "calculations.", n_keywords=0, n_subsections=1, repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF", &
                          description="Atomic Cutoff Radius Cutoff for the evaluation of the long-ranbe correction integrals. ", &
                          usage="CUTOFF {real} ", unit_str="angstrom", &
                          default_r_val=cp_unit_to_cp2k(value=6.0_dp, unit_str="angstrom"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_TAPER", &
                          description="Atomic Cutoff Radius Cutoff for Tapering the long-range correction integrals. "// &
                          "If not specified it assumes the same value specified for the CUTOFF.", &
                          usage="RC_TAPER {real} ", unit_str="angstrom", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_RANGE", &
                          description="Range of cutoff switch function (tapering): 0.5*(1-TANH((r-r0)/RC_RANGE)), "// &
                          "where r0=2.0*RC_TAPER-20.0*RC_RANGE.", &
                          usage="RC_RANGE {real} ", unit_str="angstrom", default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_lr_corr_section

! **************************************************************************************************
!> \brief Create the COULOMB se section
!> \param section the section to create
!> \date  03.2009
!> \author Teodoro Laino [tlaino]
! **************************************************************************************************
   SUBROUTINE create_coulomb_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="COULOMB", &
                          description="Setup parameters for the evaluation of the COULOMB term in SE "// &
                          "calculations.", n_keywords=0, n_subsections=1, repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="CUTOFF", &
         description="Atomic Cutoff Radius Cutoff for the evaluation of the  Coulomb integrals. "// &
         "For non-periodic calculation the default value is exactly the full cell dimension, in order "// &
         "to evaluate all pair interactions. Instead, for periodic calculations the default numerical value is used.", &
         usage="CUTOFF {real} ", unit_str="angstrom", &
         default_r_val=cp_unit_to_cp2k(value=12.0_dp, unit_str="angstrom"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_TAPER", &
                          description="Atomic Cutoff Radius Cutoff for Tapering Coulomb integrals. "// &
                          "If not specified it assumes the same value specified for the CUTOFF.", &
                          usage="RC_TAPER {real} ", unit_str="angstrom", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_RANGE", &
                          description="Range of cutoff switch function (tapering): 0.5*(1-TANH((r-r0)/RC_RANGE)), "// &
                          "where r0=2.0*RC_TAPER-20.0*RC_RANGE.", &
                          usage="RC_RANGE {real} ", unit_str="angstrom", default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_coulomb_section

! **************************************************************************************************
!> \brief Create the EXCHANGE se section
!> \param section the section to create
!> \date  03.2009
!> \author Teodoro Laino [tlaino]
! **************************************************************************************************
   SUBROUTINE create_exchange_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="EXCHANGE", &
                          description="Setup parameters for the evaluation of the EXCHANGE and "// &
                          "core Hamiltonian terms in SE calculations.", n_keywords=0, n_subsections=1, &
                          repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF", &
                          description="Atomic Cutoff Radius Cutoff for the evaluation of the Exchange integrals. "// &
                          "For non-periodic calculation the default value is exactly the full cell dimension, in order "// &
                          "to evaluate all pair interactions. Instead, for periodic calculations the default is the "// &
                          "minimum value between 1/4 of the cell dimension and the value specified in input (either"// &
                          " explicitly defined or the default numerical value).", &
                          usage="CUTOFF {real} ", unit_str="angstrom", &
                          default_r_val=cp_unit_to_cp2k(value=12.0_dp, unit_str="angstrom"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_TAPER", &
                          description="Atomic Cutoff Radius Cutoff for Tapering Exchange integrals. "// &
                          "If not specified it assumes the same value specified for the CUTOFF.", &
                          usage="RC_TAPER {real} ", unit_str="angstrom", type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_RANGE", &
                          description="Range of cutoff switch function (tapering): 0.5*(1-TANH((r-r0)/RC_RANGE)), "// &
                          "where r0=2.0*RC_TAPER-20.0*RC_RANGE.", &
                          usage="RC_RANGE {real} ", unit_str="angstrom", default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_exchange_section

! **************************************************************************************************
!> \brief Create the SCREENING se section
!> \param section the section to create
!> \date  03.2009
!> \author Teodoro Laino [tlaino]
! **************************************************************************************************
   SUBROUTINE create_screening_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="SCREENING", &
                          description="Setup parameters for the tapering of the Coulomb/Exchange Screening in "// &
                          "KDSO-D integral scheme,", n_keywords=0, n_subsections=1, repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_TAPER", &
                          description="Atomic Cutoff Radius Cutoff for Tapering the screening term. ", &
                          usage="RC_TAPER {real} ", unit_str="angstrom", &
                          default_r_val=cp_unit_to_cp2k(value=12.0_dp, unit_str="angstrom"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RC_RANGE", &
                          description="Range of cutoff switch function (tapering): 0.5*(1-TANH((r-r0)/RC_RANGE)), "// &
                          "where r0=2*RC_TAPER-20*RC_RANGE.", &
                          usage="RC_RANGE {real} ", unit_str="angstrom", default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_screening_section

! **************************************************************************************************
!> \brief Create the print se section
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_se_print_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="print", &
                          description="Section of possible print options in SE code.", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "NEIGHBOR_LISTS", &
                                       description="Activates the printing of the neighbor lists used"// &
                                       " for the periodic SE calculations.", &
                                       print_level=high_print_level, filename="", unit_str="angstrom")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "SUBCELL", &
                                       description="Activates the printing of the subcells used for the "// &
                                       "generation of neighbor lists for periodic SE.", &
                                       print_level=high_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "EWALD_INFO", &
                                       description="Activates the printing of the information for "// &
                                       "Ewald multipole summation in periodic SE.", &
                                       print_level=high_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_se_print_section

! **************************************************************************************************
!> \brief creates the input section for use with the GA part of the code
!> \param section the section to create
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! **************************************************************************************************
   SUBROUTINE create_se_ga_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="GA", &
                          description="Sets up memory parameters for the storage of the integrals", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="NCELLS", &
         description="Defines the number of linked cells for the neighbor list. "// &
         "Default value is number of processors", &
         usage="NCELLS 10", &
         default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
   END SUBROUTINE create_se_ga_section

! **************************************************************************************************
!> \brief creates the input section for the se-memory part integral storage
!> \param section the section to create
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! **************************************************************************************************
   SUBROUTINE create_se_memory_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="MEMORY", &
                          description="Sets up memory parameters for the storage of the integrals", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      NULLIFY (keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="EPS_STORAGE", &
         description="Storage threshold for compression is EPS_STORAGE", &
         usage="EPS_STORAGE 1.0E-10", &
         default_r_val=1.0E-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="MAX_MEMORY", &
         description="Defines the maximum amount of memory [MB] used to store precomputed "// &
         "(possibly compressed) two-electron two-center integrals", &
         usage="MAX_MEMORY 256", &
         default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COMPRESS", &
                          description="Enables the compression of the integrals in memory.", &
                          usage="COMPRESS <LOGICAL>", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_se_memory_section

END MODULE input_cp2k_se
