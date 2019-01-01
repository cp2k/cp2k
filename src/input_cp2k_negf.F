!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Input section for NEGF based quantum transport calculations.
! **************************************************************************************************

MODULE input_cp2k_negf
   USE bibliography,                    ONLY: Bailey2006,&
                                              Papior2017
   USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                              high_print_level
   USE input_constants,                 ONLY: negfint_method_cc,&
                                              negfint_method_simpson
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
                                              real_t
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: kelvin
   USE qs_density_mixing_types,         ONLY: create_mixing_section
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_negf'

   PUBLIC :: create_negf_section

CONTAINS

! **************************************************************************************************
!> \brief Create NEGF input section.
!> \param section input section
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE create_negf_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_negf_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, "NEGF", &
                          description="Parameters which control quantum transport calculation"// &
                          " based on Non-Equilibrium Green's Function method.", &
                          citations=(/Bailey2006, Papior2017/), &
                          n_keywords=18, n_subsections=6, repeats=.FALSE.)

      NULLIFY (keyword, print_key, subsection)

      CALL create_contact_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_atomlist_section(subsection, "SCATTERING_REGION", "Defines atoms which form the scattering region.")
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! mixing section
      CALL create_mixing_section(subsection, ls_scf=.FALSE.)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, name="DISABLE_CACHE", &
                          description="Do not keep contact self-energy matrices for future reuse", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! convergence thresholds
      CALL keyword_create(keyword, name="EPS_DENSITY", &
                          description="Target accuracy for electronic density.", &
                          n_var=1, type_of_var=real_t, default_r_val=1.0e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="EPS_GREEN", &
                          description="Target accuracy for surface Green's functions.", &
                          n_var=1, type_of_var=real_t, default_r_val=1.0e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="EPS_SCF", &
                          description="Target accuracy for SCF convergence.", &
                          n_var=1, type_of_var=real_t, default_r_val=1.0e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="EPS_GEO", &
                          description="Accuracy in mapping atoms between different force environments.", &
                          n_var=1, type_of_var=real_t, unit_str="angstrom", &
                          default_r_val=1.0e-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="ENERGY_LBOUND", &
                          description="Lower bound energy of the conductance band.", &
                          n_var=1, type_of_var=real_t, unit_str="hartree", &
                          default_r_val=-5.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="ETA", &
                          description="Infinitesimal offset from the real axis.", &
                          n_var=1, type_of_var=real_t, unit_str="hartree", &
                          default_r_val=1.0e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="HOMO_LUMO_GAP", &
                          description="The gap between the HOMO and some fictitious LUMO. This option is used as"// &
                          " an initial offset to determine the actual Fermi level of bulk contacts."// &
                          " It does not need to be exact HOMO-LUMO gap, just some value to start with.", &
                          n_var=1, type_of_var=real_t, unit_str="hartree", &
                          default_r_val=0.2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="DELTA_NPOLES", &
                          description="Number of poles of Fermi function to consider.", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="GAMMA_KT", &
                          description="Offset from the axis (in terms of k*T)"// &
                          " where poles of the Fermi function reside.", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=20)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="INTEGRATION_METHOD", &
                          description="Method to integrate Green's functions along a closed-circuit contour.", &
                          default_i_val=negfint_method_cc, &
                          enum_c_vals=s2a("CLENSHAW-CURTIS", "SIMPSON"), &
                          enum_desc=s2a( &
                          "Adaptive Clenshaw-Curtis quadrature method", &
                          "Adaptive Simpson method"), &
                          enum_i_vals=(/negfint_method_cc, negfint_method_simpson/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="INTEGRATION_MIN_POINTS", &
                          description="Initial (minimal) number of grid point for adaptive numerical integration.", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=16)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="INTEGRATION_MAX_POINTS", &
                          description="Maximal number of grid point for adaptive numerical integration.", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=512)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="MAX_SCF", &
                          description="Maximum number of SCF iterations to be performed.", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=30)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="NPROC_POINT", &
                          description="Number of MPI processes to be used per energy point. "// &
                          " Default is to use all processors (0).", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="V_SHIFT", &
                          description="Initial value of the Hartree potential shift", &
                          n_var=1, type_of_var=real_t, unit_str="hartree", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="V_SHIFT_OFFSET", &
                          description="Initial offset to determine the optimal shift in Hartree potential.", &
                          n_var=1, type_of_var=real_t, default_r_val=0.10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="V_SHIFT_MAX_ITERS", &
                          description="Maximal number of iteration to determine the optimal shift in Hartree potential.", &
                          n_var=1, type_of_var=integer_t, default_i_val=30)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! PRINT subsection
      CALL section_create(subsection, "PRINT", "Print properties for the scattering region.", &
                          repeats=.FALSE.)

      CALL create_print_dos_section(print_key, "DOS", "the Density of States (DOS)")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL create_print_dos_section(print_key, "TRANSMISSION", "the Transmission Coefficient")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_negf_section

! **************************************************************************************************
!> \brief Create NEGF%CONTACT input section.
!> \param section input section
!> \par History
!>    * 09.2017 split from create_negf_section() [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE create_contact_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_contact_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, "CONTACT", &
                          description="Section defining the contact region of NEGF setup.", &
                          n_keywords=5, n_subsections=3, repeats=.TRUE.)

      NULLIFY (keyword, print_key, subsection)

      CALL create_atomlist_section(subsection, "BULK_REGION", &
                                   "the bulk contact adjacent to the screening region.")
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_atomlist_section(subsection, "SCREENING_REGION", &
                                   "the given contact adjacent to the scattering region.")
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, name="FORCE_EVAL_SECTION", &
                          description=" Index of the FORCE_EVAL section which will be used for bulk calculation.", &
                          n_var=1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="ELECTRIC_POTENTIAL", &
                          description="External electrostatic potential applied to the given contact.", &
                          n_var=1, type_of_var=real_t, unit_str="hartree", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="FERMI_LEVEL", &
                          description="Contact Fermi level at the given temperature."// &
                          " If this keyword is not given explicitly, the Fermi level"// &
                          " will be automatically computed prior the actual NEGF calculation.", &
                          n_var=1, type_of_var=real_t, unit_str="hartree", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="REFINE_FERMI_LEVEL", &
                          description="Compute the Fermi level using the value from the FERMI_LEVEL keyword"// &
                          " as a starting point. By default the Fermi level is computed only"// &
                          " when the keyword FERMI_LEVEL is not given explicitly.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="TEMPERATURE", &
                          description="Electronic temperature.", &
                          n_var=1, type_of_var=real_t, unit_str="K", &
                          default_r_val=300.0_dp/kelvin)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! PRINT subsection
      CALL section_create(subsection, "PRINT", "Print properties for the given contact.", &
                          repeats=.FALSE.)

      CALL create_print_dos_section(print_key, "DOS", "the Density of States (DOS)")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_contact_section

! **************************************************************************************************
!> \brief Create an atomic list section.
!> \param section              NEGF section
!> \param name                 name of the new section
!> \param description          section description
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE create_atomlist_section(section, name, description)
      TYPE(section_type), POINTER                        :: section
      CHARACTER(len=*), INTENT(in)                       :: name, description

      CHARACTER(len=*), PARAMETER :: routineN = 'create_atomlist_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, TRIM(ADJUSTL(name)), &
                          description="Atoms belonging to "//TRIM(ADJUSTL(description)), &
                          n_keywords=2, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, name="LIST", &
                          description="Specifies a list of atoms.", &
                          usage="LIST {integer} {integer} .. {integer}", repeats=.TRUE., &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="MOLNAME", &
                          description="Specifies a list of named molecular fragments.", &
                          usage="MOLNAME WAT MEOH", repeats=.TRUE., &
                          n_var=-1, type_of_var=char_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
   END SUBROUTINE create_atomlist_section

! **************************************************************************************************
!> \brief Create the DOS print section.
!> \param section              section to create
!> \param name                 name of the new section
!> \param description          section description
!> \par History
!>    * 11.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE create_print_dos_section(section, name, description)
      TYPE(section_type), POINTER                        :: section
      CHARACTER(len=*), INTENT(in)                       :: name, description

      CHARACTER(len=*), PARAMETER :: routineN = 'create_print_dos_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CALL cp_print_key_section_create(section, TRIM(ADJUSTL(name)), &
                                       description="Controls the printing of "//TRIM(ADJUSTL(description))//".", &
                                       print_level=high_print_level, filename="__STD_OUT__")
      NULLIFY (keyword)

      CALL keyword_create(keyword, name="FROM_ENERGY", &
                          description="Energy point to start with.", &
                          n_var=1, type_of_var=real_t, unit_str="hartree", &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="TILL_ENERGY", &
                          description="Energy point to end with.", &
                          n_var=1, type_of_var=real_t, unit_str="hartree", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="N_GRIDPOINTS", &
                          description="Number of points to compute.", &
                          n_var=1, type_of_var=integer_t, default_i_val=201)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
   END SUBROUTINE create_print_dos_section
END MODULE input_cp2k_negf
