!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the active space section of the input
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!>      07.2024 moved out of input_cp2k_dft [JGH]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_as
   USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                              debug_print_level,&
                                              high_print_level,&
                                              low_print_level,&
                                              medium_print_level,&
                                              silent_print_level
   USE input_constants,                 ONLY: &
        casci_canonical, dmft_model, eri_method_full_gpw, eri_method_gpw_ht, eri_operator_coulomb, &
        eri_operator_erf, eri_operator_erfc, eri_operator_gaussian, eri_operator_trunc, &
        eri_operator_yukawa, gaussian, hf_model, manual_selection, mao_projection, no_solver, &
        qiskit_solver, rsdft_model, wannier_projection
   USE input_cp2k_loc,                  ONLY: create_localize_section
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_as'

   PUBLIC :: create_active_space_section

CONTAINS

! **************************************************************************************************
!> \brief      Create CP2K input section for the calculation of an active space Hamiltonian
!> \param section ...
!> \par History:
!>      - Creation 06.04.2016
!> \author     JHU
! **************************************************************************************************
   SUBROUTINE create_active_space_section(section)

      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, name="ACTIVE_SPACE", &
                          description="Define parameters and method to calculate an electronic active space", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword, subsection, print_key)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="_SECTION_PARAMETERS_", &
                          description="Controls the activation of the ACTIVE_SPACE section", &
                          usage="&ACTIVE_SPACE ON", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MODEL", &
                          description="Model used to calculate active space interaction Hamiltonian.", &
                          usage="MODEL HARTREE-FOCK", &
                          default_i_val=hf_model, &
                          enum_c_vals=s2a("HARTREE-FOCK", "RSDFT", "DMFT"), &
                          enum_i_vals=(/hf_model, rsdft_model, dmft_model/), &
                          enum_desc=s2a("Hartree-Fock model for interaction Hamiltonian", &
                                        "Range-separated DFT model for interaction Hamiltonian", &
                                        "DMFT model Hamiltonian"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ACTIVE_ELECTRONS", &
                          description="The number of active electrons in the CAS space", &
                          usage="ACTIVE_ELECTRONS 4", n_var=1, default_i_val=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ACTIVE_ORBITALS", &
                          description="The number of active orbitals defining the CAS space.", &
                          usage="ACTIVE_ORBITALS 2", n_var=1, default_i_val=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ACTIVE_ORBITAL_INDICES", &
                          description="The indices of the active orbitals. Requires ORBITAL_SELECTION MANUAL!", &
                          usage="ACTIVE_ORBITAL_INDICES 2 3 {...}", n_var=-1, default_i_vals=(/-1/), &
                          type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ISOLATED_SYSTEM", &
                          description="System is treated without any periodic boundary conditions.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "FCIDUMP", &
                                       description="Controls the writing of a file in FCIDUMP format.", &
                                       print_level=high_print_level, filename="")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL keyword_create(keyword, __LOCATION__, name="ORBITAL_SELECTION", &
                          description="Method used to select active space orbitals.", &
                          usage="ORBITAL_SELECTION CANONICAL", &
                          default_i_val=casci_canonical, &
                          enum_c_vals=s2a("CANONICAL", "WANNIER_PROJECTION", "MAO", "MANUAL"), &
                          enum_i_vals=(/casci_canonical, wannier_projection, mao_projection, manual_selection/), &
                          enum_desc=s2a("Select orbitals using energy ordering of canoncial orbitals", &
                                        "Select orbitals from projected Wannier functions", &
                                        "Select orbitals from modified atomic orbitals", &
                                        "Select orbitals manually via ACTIVE_ORBITAL_INDICES"))

      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SUBSPACE_ATOM", &
                          description="Number of atom that defines the subspace to be projected on.", &
                          usage="SUBSPACE_ATOM x", default_i_val=-1, &
                          type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SUBSPACE_SHELL", &
                          description="Shell definition for subsapce.", &
                          usage="SUBSPACE_SHELL 3d4s", default_c_val="X", &
                          type_of_var=char_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SCF_EMBEDDING", &
                          description="Whether to turn on the self-consistent embedding scheme", &
                          default_l_val=.FALSE., lone_keyword_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="QCSCHEMA", &
                          description="Name of the QCSchema file, may include a path", &
                          usage="QCSCHEMA <FILENAME>", &
                          type_of_var=lchar_t, repeats=.FALSE., &
                          default_lc_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="AS_SOLVER", &
                          description="The external active space solver for the embedding approach", &
                          usage="AS_SOLVER QISKIT", &
                          default_i_val=no_solver, &
                          enum_c_vals=s2a("NONE", "QISKIT"), &
                          enum_i_vals=(/no_solver, qiskit_solver/), &
                          enum_desc=s2a("NO solver, used to produce FCIDUMP/QCSchema files", &
                                        "QISKIT active space solver"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ITER", &
                          description="Energy convergence threshold of the DFT embedding scheme.", &
                          usage="EPS_ITER 1.0E-6 ", type_of_var=real_t, &
                          default_r_val=1.0E-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Max number of iterations for the DFT embedding scheme.", &
                          usage="MAX_ITER 50", type_of_var=integer_t, &
                          default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_print_orb_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_eri_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_eri_gpw(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_localize_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_socket_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_active_space_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_socket_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="SOCKET", &
                          description="Parameters to set up the socket communicating to the external active space solver.", &
                          n_keywords=3, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="INET", &
                          description="Use an INET socket rather than a UNIX socket.", &
                          usage="INET <LOGICAL>", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PORT", &
                          description="Port number for the socket client.", &
                          usage="port <INTEGER>", &
                          default_i_val=12345)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HOST", &
                          description="Host name for the socket client.", &
                          usage="host <HOSTNAME>", &
                          default_c_val="embedding_socket")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_socket_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_print_orb_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PRINT_ORBITAL_CUBES", &
                          description="Controls printing of active orbital cube files.", &
                          n_keywords=5, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="FILENAME", &
                          description="Body of Filename for the cube files.", &
                          usage="FILENAME {name}", default_c_val="ActiveOrbital", &
                          type_of_var=char_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ALIST", &
                          description="List of alpha orbitals to be printed. -1 defaults to all values", &
                          usage="ALIST {1 2 3 ...}", n_var=-1, default_i_vals=(/-1/), &
                          lone_keyword_i_val=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BLIST", &
                          description="List of beta orbitals to be printed. -1 defaults to all values", &
                          usage="ALIST {1 2 3 ...}", n_var=-1, default_i_vals=(/-1/), &
                          lone_keyword_i_val=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file"// &
                          " (larger values result in smaller cube files)."// &
                          " You can provide 3 numbers (for X,Y,Z) or 1 number valid for all components", &
                          usage="STRIDE {2 2 2}", n_var=-1, default_i_vals=(/2, 2, 2/), &
                          type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STOP_AFTER_CUBES", &
                          description="Whether to stop the computation after printing the cubes.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_print_orb_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_eri_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="ERI", &
                          description="Parameters for the electron repulsion integrals.", &
                          n_keywords=5, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="METHOD", &
                          description="Method used in ERI calculation.", &
                          usage="METHOD FULL_GPW", &
                          enum_c_vals=s2a("FULL_GPW", "GPW_HALF_TRANSFORM"), &
                          enum_i_vals=(/eri_method_full_gpw, eri_method_gpw_ht/), &
                          enum_desc=s2a("Use the GPW approach with MOs", &
                                        "Use the GPW approach for half-transformed MO ERIs"), &
                          default_i_val=eri_method_full_gpw)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OPERATOR", &
                          description="Operator used in ERI calculation.", &
                          usage="OPERATOR <1/R>", &
                          enum_c_vals=s2a("<1/R>", "<EXP(-A*R)/R>", "<ERF(A*R)/R>", &
                                          "<ERFC(A*R)/R>", "<EXP(-A*R2)/R>", "<H(A-R)/R>"), &
                          enum_i_vals=(/eri_operator_coulomb, eri_operator_yukawa, &
                                        eri_operator_erf, eri_operator_erfc, eri_operator_gaussian, eri_operator_trunc/), &
                          enum_desc=s2a("Coulomb operator", &
                                        "Yukawa potential operator", &
                                        "Error function potential operator", &
                                        "Complementary error function potential operator", &
                                        "Gaussian potential operator", &
                                        "Truncated Coulomb potential"), &
                          default_i_val=eri_operator_coulomb)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OPERATOR_PARAMETER", &
                          description="Range parameter for ERI operator.", &
                          usage="OPERATOR_PARAMETER 4.0", type_of_var=real_t, &
                          default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PERIODICITY", &
                          description="Periodicity used for operators in ERI calclulation.", &
                          usage="PERIODICITY {1 1 1}", n_var=-1, default_i_vals=(/1, 1, 1/), &
                          type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_RADIUS", &
                          description="Radius of operator interactions in non-periodic cases.", &
                          usage="CUTOFF_RADIUS 20", type_of_var=real_t, &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EPS_INTEGRAL", &
         description="Accuracy of ERIs that will be stored.", &
         usage="EPS_INTEGRAL 1.0E-10 ", type_of_var=real_t, &
         default_r_val=1.0E-12_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_eri_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_eri_gpw(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="ERI_GPW", &
                          description="Parameters for the GPW approach to electron repulsion integrals.", &
                          n_keywords=5, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="EPS_GRID", &
                          description="Determines a threshold for the GPW based integration", &
                          usage="EPS_GRID 1.0E-9 ", type_of_var=real_t, &
                          default_r_val=1.0E-8_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          description="Determines a threshold for the sparse matrix multiplications if METHOD "// &
                          "GPW_HALF_TRANSFORM is used", &
                          usage="EPS_FILTER 1.0E-9 ", type_of_var=real_t, &
                          default_r_val=1.0E-9_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF", &
                          description="The cutoff of the finest grid level in the GPW integration.", &
                          usage="CUTOFF 300", type_of_var=real_t, &
                          default_r_val=300.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REL_CUTOFF", &
                          variants=(/"RELATIVE_CUTOFF"/), &
                          description="Determines the grid at which a Gaussian is mapped.", &
                          usage="REL_CUTOFF 50", type_of_var=real_t, &
                          default_r_val=50.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STORE_WFN", &
                          variants=(/"STORE_WAVEFUNCTION"/), &
                          description="Store wavefunction in real space representation for integration.", &
                          usage="STORE_WFN  T", type_of_var=logical_t, &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GROUP_SIZE", &
                          description="Sets the size of a subgroup for ERI calculation, "// &
                          "each of which with a full set of work grids, arrays or orbitals "// &
                          "depending on the method of grids (work grids, arrays, orbitals). "// &
                          "Small numbers reduce communication but increase the memory demands. "// &
                          "A negative number indicates all processes (default).", &
                          usage="GROUP_SIZE 2", type_of_var=integer_t, &
                          default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PRINT_LEVEL", &
                          variants=(/"IOLEVEL"/), &
                          description="How much output is written by the individual groups.", &
                          usage="PRINT_LEVEL HIGH", &
                          default_i_val=silent_print_level, enum_c_vals= &
                          s2a("SILENT", "LOW", "MEDIUM", "HIGH", "DEBUG"), &
                          enum_desc=s2a("Almost no output", &
                                        "Little output", "Quite some output", "Lots of output", &
                                        "Everything is written out, useful for debugging purposes only"), &
                          enum_i_vals=(/silent_print_level, low_print_level, medium_print_level, &
                                        high_print_level, debug_print_level/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_eri_gpw

END MODULE input_cp2k_as
