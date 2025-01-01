!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the XAS section of the input
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!>      07.2024 moved out of input_cp2k_dft [JGH]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_xas
   USE bibliography,                    ONLY: Iannuzzi2007,&
                                              Shigeta2001
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              debug_print_level,&
                                              high_print_level,&
                                              low_print_level
   USE input_constants,                 ONLY: &
        do_potential_coulomb, do_potential_id, do_potential_short, do_potential_truncated, &
        gaussian, ot_mini_cg, ot_mini_diis, tddfpt_singlet, tddfpt_spin_cons, tddfpt_spin_flip, &
        tddfpt_triplet, xas_1s_type, xas_2p_type, xas_2s_type, xas_3d_type, xas_3p_type, &
        xas_3s_type, xas_4d_type, xas_4f_type, xas_4p_type, xas_4s_type, xas_dip_len, xas_dip_vel, &
        xas_dscf, xas_none, xas_not_excited, xas_tdp_by_index, xas_tdp_by_kind, xas_tp_fh, &
        xas_tp_flex, xas_tp_hh, xas_tp_xfh, xas_tp_xhh, xes_tp_val
   USE input_cp2k_loc,                  ONLY: create_localize_section,&
                                              print_wanniers
   USE input_cp2k_print_dft,            ONLY: create_pdos_section
   USE input_cp2k_scf,                  ONLY: create_scf_section
   USE input_cp2k_xc,                   ONLY: create_xc_fun_section
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
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_xas'

   PUBLIC :: create_xas_section, create_xas_tdp_section

CONTAINS

! **************************************************************************************************
!> \brief makes the input section for core-level spectroscopy simulations
!> \param section ...
!> \par History
!>      03.2005 created [MI]
! **************************************************************************************************
   SUBROUTINE create_xas_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="xas", &
                          description="Sets the method of choice to calculate core-level excitation spectra. "// &
                          "The occupied states from  which we calculate the "// &
                          "excitation should be specified. "// &
                          "Localization of the orbitals may be useful.", &
                          n_keywords=10, n_subsections=1, repeats=.FALSE., &
                          citations=(/Iannuzzi2007/))

      NULLIFY (keyword, subsection, print_key)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of core-level spectroscopy simulations", &
                          usage="&XAS T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="METHOD", &
                          variants=(/"XAS_METHOD"/), &
                          description="Method to be used to calculate core-level excitation spectra", &
                          usage="METHOD TP_HH", &
                          default_i_val=xas_none, &
                          enum_c_vals=s2a("NONE", "TP_HH", "TP_FH", "TP_VAL", "TP_XHH", "TP_XFH", "DSCF", "TP_FLEX"), &
                          enum_desc=s2a( &
                          "No core electron spectroscopy", "Transition potential half-hole", &
                          "Transition potential full-hole", "Hole in homo for X-ray emission only ", &
                          "Transition potential excited half-hole", &
                          "Transition potential excited full-hole ", &
                          "DSCF calculations to compute the first (core)excited state", &
                          "Transition potential with generalized core occupation and total number of electrons"), &
                          enum_i_vals=(/xas_none, xas_tp_hh, xas_tp_fh, xes_tp_val, xas_tp_xhh, &
                                        xas_tp_xfh, xas_dscf, xas_tp_flex/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XAS_CORE", &
                          description="Occupation of the core state in XAS calculation by TP_FLEX.", &
                          usage="XAS_CORE 0.5", &
                          default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XAS_TOT_EL", &
                          description="Total number of electrons for spin channel alpha, in XAS calculation by TP_FLEX. "// &
                          "If it is a negative value, the number of electrons is set to GS number of electrons "// &
                          "minus the amount subtracted from the core state", &
                          usage="XAS_TOT_EL 10", &
                          default_r_val=-1._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XES_CORE", &
                          description="Occupation of the core state in XES calculation by TP_VAL. "// &
                          "The HOMO is emptied by the same amount.", &
                          usage="XES_CORE 0.5", &
                          default_r_val=1._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XES_EMPTY_HOMO", &
                          description="Set the occupation of the HOMO in XES calculation by TP_VAL. "// &
                          "The HOMO can be emptied or not, if the core is still full.", &
                          usage="XES_EMPTY_HOMO", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DIPOLE_FORM", &
                          variants=(/"DIP_FORM"/), &
                          description="Type of integral to get the oscillator strengths "// &
                          "in the diipole approximation", &
                          usage="DIPOLE_FORM string", &
                          default_i_val=xas_dip_vel, &
                          enum_c_vals=s2a("LENGTH", "VELOCITY"), &
                          enum_desc=s2a("Length form &lang; i | e r | j &rang;", &
                                        "Velocity form &lang; i | d/dr | j &rang;"), &
                          enum_i_vals=(/xas_dip_len, xas_dip_vel/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

!   replace the specialized keyword with standard scf section
!    scf_env is added to xas_env

      NULLIFY (subsection)
      CALL create_scf_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, __LOCATION__, name="STATE_TYPE", &
                          variants=(/"TYPE"/), &
                          description="Type of the orbitals that are excited for the xas spectra calculation", &
                          usage="STATE_TYPE 1S", &
                          default_i_val=xas_1s_type, &
                          enum_c_vals=s2a("1S", "2S", "2P", "3S", "3P", "3D", "4S", "4P", "4D", "4F"), &
                          enum_desc=s2a("1s orbitals", "2s orbitals", "2p orbitals", "3s orbitals", "3p orbitals", &
                                        "3d orbitals", "4s orbitals", "4p orbitals", "4d orbitals", "4f orbitals"), &
                          enum_i_vals=(/xas_1s_type, xas_2s_type, xas_2p_type, xas_3s_type, xas_3p_type, xas_3d_type, &
                                        xas_4s_type, xas_4p_type, xas_4d_type, xas_4f_type/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STATE_SEARCH", &
                          description="# of states where to look for the one to be excited", &
                          usage="STATE_SEARCH 1", &
                          default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SPIN_CHANNEL", &
                          description="# Spin channel of the excited orbital", &
                          usage="SPIN_CHANNEL 1", &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMS_LIST", &
                          variants=(/"AT_LIST"/), &
                          description="Indexes of the atoms to be excited. "// &
                          "This keyword can be repeated several times "// &
                          "(useful if you have to specify many indexes).", &
                          usage="ATOMS_LIST {integer}  {integer} ..  {integer} ", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OVERLAP_THRESHOLD", &
                          description="Threshold for including more than one initial core excited state "// &
                          "per atom. The threshold is taken relative to the maximum overlap.", &
                          usage="OVERLAP_THRESHOLD 8.e-1", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORBITAL_LIST", &
                          variants=(/"ORBITAL_LIST"/), &
                          description="Indices of the localized orbitals to be excited. "// &
                          "This keyword can be repeated several times "// &
                          "(useful if you have to specify many indexes).", &
                          usage="ORBITAL_LIST {integer}  {integer} ..  {integer} ", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ADDED_MOS", &
                          description="Number of additional MOS added spin up only", &
                          usage="ADDED_MOS {integer}", default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER_ADDED", &
                          description="maximum number of iteration in calculation of added orbitals", &
                          usage="MAX_ITER_ADDED 100", default_i_val=2999)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ADDED", &
                          description="target accuracy incalculation of the added orbitals", &
                          usage="EPS_ADDED 1.e-6", default_r_val=1.0e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NGAUSS", &
                          description="Number of gto's for the expansion of the STO "// &
                          "of the type given by STATE_TYPE", &
                          usage="NGAUSS {integer}", default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTART", &
                          description="Restart the excited state if the restart file exists", &
                          usage="RESTART", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WFN_RESTART_FILE_NAME", &
                          variants=(/"RESTART_FILE_NAME"/), &
                          description="Root of the file names where to read the MOS from "// &
                          "which to restart the calculation of the core level excited states", &
                          usage="WFN_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_localize_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="PRINT", &
                          description="printing of information during the core-level spectroscopy simulation", &
                          repeats=.FALSE.)

      ! Add printing of wannier infos
      CALL print_wanniers(subsection)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "iteration_info", &
                                       description="Controls the printing of basic iteration information during the XAS SCF.", &
                                       print_level=low_print_level, filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, name="time_cumul", &
                          description="If the printkey is activated switches the printing of timings"// &
                          " to cumulative (over the SCF).", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "program_run_info", &
                                       description="Controls the printing of basic iteration information in CLS", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "XES_SPECTRUM", &
                                       description="Controls the dumping of the CLS output files containing the emission spectra", &
                                       print_level=low_print_level, common_iter_levels=3, filename="")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create( &
         print_key, __LOCATION__, "XAS_SPECTRUM", &
         description="Controls the dumping of the CLS output files containing the absorption spectra", &
         print_level=low_print_level, common_iter_levels=3, filename="")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL create_pdos_section(print_key)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESTART", &
                                       description="Controls the dumping of MO restart file during the SCF. "// &
                                       "of a Core-Level-Spectroscopy calculation. For each new excited atom, "// &
                                       "one different restart file is dumped. These restart files should be "// &
                                       "employed only to restart the same type of CLS calculation, "// &
                                       "i.e. with the same core potential.", &
                                       print_level=low_print_level, common_iter_levels=3, each_iter_names=s2a("XAS_SCF"), &
                                       add_last=add_last_numeric, each_iter_values=(/3/), filename="")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "FULL_RESTART", &
                                       description="Controls the dumping of a standard MO restart file "// &
                                       "where coefficients and occupation numbers are those of the TP scheme, "// &
                                       "i.e. with emptied core state.", &
                                       print_level=high_print_level, common_iter_levels=3, each_iter_names=s2a("XAS_SCF"), &
                                       add_last=add_last_numeric, each_iter_values=(/3/), filename="")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "CLS_FUNCTION_CUBES", &
                                       description="Controls the printing of the relaxed orbitals ", &
                                       print_level=high_print_level, common_iter_levels=3, add_last=add_last_numeric, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUBES_LU_BOUNDS", &
                          variants=(/"CUBES_LU"/), &
                          description="The lower and upper index of the states to be printed as cube", &
                          usage="CUBES_LU_BOUNDS integer integer", &
                          n_var=2, default_i_vals=(/0, -2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUBES_LIST", &
                          description="Indexes of the states to be printed as cube files "// &
                          "This keyword can be repeated several times "// &
                          "(useful if you have to specify many indexes).", &
                          usage="CUBES_LIST 1 2", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_xas_section

! **************************************************************************************************
!> \brief makes the input section for core-level spectroscopy simulations using
!>        linear response TDDFT
!> \param section ...
!> \par History
!>      11.2017 created [AB]
! **************************************************************************************************
   SUBROUTINE create_xas_tdp_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection, subsubsection, &
                                                            subsubsubsection

      NULLIFY (keyword, print_key, subsection, subsubsection, subsubsubsection)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="XAS_TDP", &
                          description="XAS simulations using linear-response TDDFT. Excitation from "// &
                          "specified core orbitals are considered one at a time. In case of high "// &
                          "symmetry structures, donor core orbitals should be localized.", &
                          n_keywords=19, n_subsections=4, repeats=.FALSE.)

      NULLIFY (keyword, subsection, print_key)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of XAS simulations with linear "// &
                          "response TDDFT.", &
                          usage="&TDP_XAS {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHECK_ONLY", &
                          description="This keyword defines whether the full calculation should "// &
                          "be done or not. If set to .TRUE., only the determination "// &
                          "of donor MOs is conducted. This run option allows for "// &
                          "cheap verification of the input parameters", &
                          usage="CHECK_ONLY {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE., &
                          repeats=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTART_FROM_FILE", &
                          variants=s2a("RESTART_FILENAME", "RST_FILENAME", "RESTART_FILE", "RST_FILE"), &
                          description="By providing a RESTART file containing the linear-response "// &
                          "orbitals and excitations energies from a previous calculation, "// &
                          "all computations are skipped except for the corresponding "// &
                          "PDOS and/or CUBE file printing as defined in the PRINT "// &
                          "subsection. Basis sets and geometry need to be consistent.", &
                          usage="RESTART_FROM_FILE <FILENAME>", &
                          type_of_var=char_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EXCITATIONS", &
                          variants=(/"EXCITATION"/), &
                          description="Specify the type of excitation to consider. In case of a "// &
                          "resctricted closed-shell ground state calculation, "// &
                          "RCS_SINGLET or/and RCS_TRIPLET can be chosen. In case of a "// &
                          "open-shell ground state calculation (either UKS or ROKS), "// &
                          "standard spin conserving excitation (OS_SPIN_CONS) or/and "// &
                          "spin-flip excitation (OS_SPIN_FLIP) can be chosen.", &
                          usage="EXCITATIONS {string}", &
                          repeats=.TRUE., &
                          default_i_val=tddfpt_singlet, &
                          enum_c_vals=s2a("RCS_SINGLET", "RCS_TRIPLET", "OS_SPIN_CONS", "OS_SPIN_FLIP"), &
                          enum_desc=s2a("Singlet excitation on top of restricted closed-shell ground state", &
                                        "Triplet excitation on top of restricted closed-shell ground state", &
                                        "Spin-conserving excitations on top of open-shell ground state", &
                                        "Spin-flip excitation on top of open-shell ground state"), &
                          enum_i_vals=(/tddfpt_singlet, tddfpt_triplet, tddfpt_spin_cons, tddfpt_spin_flip/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_PGF_XAS", &
                          variants=s2a("EPS_PGF", "EPS_PGF_XAS_TDP"), &
                          description="The threshold used to determine the spacial extent of all "// &
                          "primitive Gaussian functions used for the construction "// &
                          "of neighbor lists in the XAS_TDP method. "// &
                          "By default, takes the value of QS%EPS_PGF_ORB. Useful if "// &
                          "the former value is tiny due to possible ground state HFX "// &
                          "contributions.", &
                          usage="EPS_PGS_XAS {real}", &
                          type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          variants=s2a("EPS_FILTER_MATRIX"), &
                          description="The threshold used for sparse matrix operations", &
                          usage="EPS_FILTER {real}", &
                          type_of_var=real_t, &
                          default_r_val=1.0E-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DIPOLE_FORM", &
                          variants=(/"DIP_FORM"/), &
                          description="Type of integral to get the oscillator strengths "// &
                          "in the dipole approximation", &
                          usage="DIPOLE_FORM {string}", &
                          default_i_val=xas_dip_vel, &
                          enum_c_vals=s2a("LENGTH", "VELOCITY"), &
                          enum_desc=s2a("Length form &lang; 0 | e r | n &rang;", &
                                        "Velocity form &lang; 0 | d/dr | n &rang;"), &
                          enum_i_vals=(/xas_dip_len, xas_dip_vel/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="QUADRUPOLE", &
                          variants=s2a("DO_QUADRUPOLE", "DO_QUAD", "QUAD"), &
                          description="Compute the electric quadrupole contribution to the "// &
                          "oscillator strenghts (in the length representation with "// &
                          "the origin set on the relevant excited atom)", &
                          usage="QUADRUPOLE {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XYZ_DIPOLE", &
                          variants=s2a("DIPOLE_XYZ"), &
                          description="Whether the detailed contributions of the dipole oscillator "// &
                          "strengths along the X,Y,Z directions should be printed.", &
                          usage="XYZ_DIPOLE {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

!  the GW2X correction subsection
      CALL section_create(subsection, __LOCATION__, name="GW2X", &
                          description="Specifications for the GW2X calculation of core "// &
                          "ionization potentials. Note that donor states need to be actively "// &
                          "localized using the LOCALIZE keyword in DONOR_STATES. N_SEARCH "// &
                          "should be kept to the minimum, such that only core states are localized.", &
                          citations=(/Shigeta2001/), &
                          n_keywords=8, &
                          n_subsections=0, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Enables the GW2X correction of the core ionization potentials", &
                          usage="&GW2X {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XPS_ONLY", &
                          description="If set to .TRUE., only run GW2X calculations for XPS "// &
                          "spectroscopy and ignore all XAS calculations. It is still "// &
                          "required to define the DONOR_STATES and KERNEL%EXACT_EXCHANGE "// &
                          "subsections.", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BATCH_SIZE", &
                          description="MOs batch size for batched tensor contraction. Larger "// &
                          "size is faster, but uses more memory. Default should be safe "// &
                          "in most cases.", &
                          default_i_val=64)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_GW2X", &
                          description="Convergence threshold for GW2X iterations (in eV)", &
                          default_r_val=1.E-2_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="C_OS", &
                          description="Opposite-spin scling factor. SCS => 6/5, SOS => 1.3", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="C_SS", &
                          description="Same-spin scling factor. SCS => 1/3, SOS => 0.0", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_GW2X_ITER", &
                          description="Maximum number of iterations for GW2X", &
                          default_i_val=10)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PSEUDO_CANONICAL", &
                          variants=(/"PSEUDO_CANO"/), &
                          description="Whether the pseudo-canonical version of GW2X should be used "// &
                          "(versus only using the diagonal of the generalized Fock matrix)", &
                          default_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

!  The donor state subsection

      CALL section_create(subsection, __LOCATION__, name="DONOR_STATES", &
                          description="Specifications for the donor states from which core "// &
                          "electrons are excited", &
                          n_keywords=6, &
                          n_subsections=0, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="DEFINE_EXCITED", &
                          description="Whether the atoms to be excited should be defined by "// &
                          "a list of atom indices or by a list of atom kinds.", &
                          usage="DEFINE_EXCITED {string}", &
                          default_i_val=xas_tdp_by_index, &
                          enum_c_vals=s2a("BY_INDEX", "BY_KIND"), &
                          enum_i_vals=(/xas_tdp_by_index, xas_tdp_by_kind/), &
                          enum_desc=s2a("Excited atoms are defined by a list of indices", &
                                        "Excited atoms are defined by a list of atomic kinds"))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOM_LIST", &
                          variants=(/"AT_LIST"/), &
                          description="Indices of the atoms to be excited. "// &
                          "Keyword only taken into account if DEFINE_EXCITED = BY_INDEX", &
                          usage="ATOM_LIST {integer}  {integer} ..  {integer} ", &
                          n_var=-1, type_of_var=integer_t, repeats=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="KIND_LIST", &
                          description="Kind of atoms to be excited. "// &
                          "All atoms of the specified kinds are considered. "// &
                          "Keyword only taken into account if DEFINE_EXCITED = BY_KIND", &
                          usage="KIND_LIST {string}  {string} ..  {string} ", &
                          n_var=-1, type_of_var=char_t, repeats=.FALSE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STATE_TYPES", &
                          variants=(/"TYPES"/), &
                          description="Types of orbitals that are excited, for each atom/kind, "// &
                          "in order to do LR-TDDFT driven xas spectra calculation. "// &
                          "This keyword MUST have the same number of entries as the relevant "// &
                          "KIND_LIST or ATOM_LIST. The order of the specified state types must "// &
                          "correspond to the order of the relevant kinds/indices. "// &
                          "This keyword can be repeated, useful when multiple orbital types "// &
                          "should be excited for specific kinds/atoms.", &
                          n_var=-1, default_i_val=xas_not_excited, repeats=.TRUE., &
                          usage="STATE_TYPES {string}  {string} .. {string}", &
                          enum_c_vals=s2a("1S", "2S", "2P", "NE"), &
                          enum_desc=s2a("1s orbital", "2s orbital", "2p orbitals", "not excited"), &
                          enum_i_vals=(/xas_1s_type, xas_2s_type, xas_2p_type, xas_not_excited/))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="N_SEARCH", &
                          description="Number of MOs (per spin) to search to find specified donor core "// &
                          "orbitals, starting from the lowest in energy and upward. By default, "// &
                          "all HOMOs are searched. If the LOCALIZE keyword is used, "// &
                          "then all searched states are first localized.", &
                          usage="N_SEARCH {integer}", &
                          default_i_val=-1)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LOCALIZE", &
                          variants=s2a("LOC", "DO_LOC"), &
                          description="Whether the N_SEARCH potential donor states should be "// &
                          "actively localized. Necessary in case of excited atoms "// &
                          "equivalent under symmetry or GW2X correction.", &
                          usage="LOCALIZE {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
! End of the donor states subsection

! The OT solver subsection
      CALL section_create(subsection, __LOCATION__, name="OT_SOLVER", &
                          description="Specifications for the iterative OT solver. Note: only "// &
                          "available within the Tamm-Dancoff approximation. Recommanded if excitations "// &
                          "from multiple donor states take place.", &
                          n_keywords=4, &
                          n_subsections=0, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Enables the usage of the OT iterator solver", &
                          usage="&OT_SOLVER {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximum number of iterations allowed for the OT solver", &
                          usage="MAX_ITER {integer}", &
                          default_i_val=50)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ITER", &
                          description="Convergence threshold for the OT solver", &
                          usage="EPS_ITER {double}", &
                          default_r_val=1.0E-4_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MINIMIZER", &
                          description="Minimizer to be used with the OT solver", &
                          usage="MINIMIZER DIIS", &
                          default_i_val=ot_mini_diis, &
                          enum_c_vals=s2a("CG", "DIIS"), &
                          enum_desc=s2a("Conjugated gradient: safer", &
                                        "Direct inversion of the iterative subspace: faster"), &
                          enum_i_vals=(/ot_mini_cg, ot_mini_diis/))
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
! End of the OT solver subsection

      CALL keyword_create(keyword, __LOCATION__, name="SPIN_ORBIT_COUPLING", &
                          variants=(/"SOC"/), &
                          description="Whether spin-orbit coupling should be added. "// &
                          "Note: only applies for spin-restricted calculations with "// &
                          "singlet and triplet excitations OR spin-unrestricted "// &
                          "calculations with both spin-conserving and spin-flip.", &
                          usage="SOC {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TAMM_DANCOFF", &
                          variants=(/"TDA"/), &
                          description="Whether the Tamm-Dancoff approximation should be used.", &
                          usage="TAMM_DANCOFF {logical}", &
                          default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GRID", &
                          variants=(/"ATOMIC_GRID"/), &
                          description="Specification of the atomic angular and radial grids for "// &
                          "a given atomic kind. This keyword can/should be repeated "// &
                          "for each excited kind. The default grid dimensions are "// &
                          "those set for the GAPW ground state calculation. These "// &
                          "grids are used for the xc-kernel integration. "// &
                          "Usage: GRID < KIND > < LEBEDEV_GRID > < RADIAL_GRID >", &
                          usage="GRID {string} {integer} {integer}", &
                          n_var=3, type_of_var=char_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="N_EXCITED", &
                          variants=(/"N_ROOTS"/), &
                          description="The number of excited states to compute per donor "// &
                          "molecular orbital. (e.g. if 2p excitations, "// &
                          "3*N_EXCITED excited states are considered). "// &
                          "If N_EXCITED is set to -1, all excitations are considered", &
                          usage="N_EXCITED {integer}", &
                          default_i_val=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ENERGY_RANGE", &
                          variants=s2a("E_RANGE"), &
                          description="The energy range in eV for which excitations are considered. "// &
                          "Only excitated states within the range of: first excitation "// &
                          "energy + ENERGY_RANGE are kept. If ENERGY_RANGE "// &
                          "and N_EXCITED are specified, the former has priority. "// &
                          "Negative values are ignored and N_EXCITED takes over.", &
                          usage="ENERGY_RANGE {real}", &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

!  The KERNEL subsection
      CALL section_create(subsection, __LOCATION__, name="KERNEL", &
                          description="Defines how the kernel is built in terms of functionals.", &
                          n_keywords=1, &
                          n_subsections=1, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="RI_REGION", &
                          variants=(/"RI_RADIUS"/), &
                          description="The region defined by a sphere of the given radius around "// &
                          "the excited atom defining which RI_XAS basis elements are "// &
                          "considered for the RI projection of the density. Each basis "// &
                          "element which center is in this region is taken. The density "// &
                          "for a given excited atom is expressed as : "// &
                          "sum_ijkl P_ij (ijk) S_kl^-1 xi_l, where P_ij is the density "// &
                          "matrix, i,j span the orbital basis and k,l the RI_XAS basis "// &
                          "in the region. The larger the radius, the more basis "// &
                          "functions to expand the density. However, it is assumed "// &
                          "that it is a small number and the code does not scale well "// &
                          "as the number of RI basis elements gets too large. "// &
                          "Expressed in Angstrom. If the radius is set to 0.0, only "// &
                          "the RI basis elements centered on the excited atom are used.", &
                          usage="RI_REGION {real}", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

!  The XC_FUNCTIONAL subsubsection
      CALL create_xc_fun_section(subsubsection)
      CALL section_add_subsection(subsection, subsubsection)
      CALL section_release(subsubsection)

! The EXACT_EXCHANGE subsubsection
      CALL section_create(subsubsection, __LOCATION__, name="EXACT_EXCHANGE", &
                          description="Whether exact-exchange should be added to the kernel and "// &
                          "if so, with which fraction and operator.", &
                          n_keywords=7, &
                          n_subsections=1, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Enables the addition of exact exchange to the kernel.", &
                          usage="&EXACT_EXCHANGE {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="POTENTIAL_TYPE", &
                          variants=s2a("OP", "OPERATOR", "POTENTIAL"), &
                          description="The type of operator used for exact exchange. The standard "// &
                          "Coulomb operator cannot be used in periodic systems.", &
                          usage="OPERATOR {string}", &
                          repeats=.FALSE., &
                          default_i_val=do_potential_coulomb, &
                          enum_c_vals=s2a("COULOMB", "TRUNCATED", "SHORTRANGE"), &
                          enum_desc=s2a("Standard Coulomb operator: 1/r", &
                                        "Truncated Coulomb operator: 1/r if r < R_c, 0 otherwise ", &
                                        "Short range: erfc(omega*r)/r"), &
                          enum_i_vals=(/do_potential_coulomb, do_potential_truncated, &
                                        do_potential_short/))
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_RADIUS", &
                          variants=s2a("R_C", "RC", "RANGE"), &
                          description="The cutoff radius (in Angstrom) for the truncated Coulomb operator.", &
                          usage="CUTOFF_RADIUS {double}", &
                          default_r_val=0.0_dp, &
                          repeats=.FALSE.)
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="T_C_G_DATA", &
                          description="Location of the file t_c_g.dat that contains the data for the "// &
                          "evaluation of the truncated gamma function ", &
                          usage="T_C_G_DATA {string}", &
                          default_c_val="t_c_g.dat")
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OMEGA", &
                          description="The range parameter for the short range operator (in 1/a0).", &
                          usage="OMEGA {double}", &
                          default_r_val=0.0_dp, &
                          repeats=.FALSE.)
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_RANGE", &
                          description="The threshold to determine the effective range of the short range "// &
                          "operator: erfc(omega*eff_range)/eff_range = EPS_RANGE", &
                          usage="EPS_RANGE = {double}", &
                          default_r_val=1.0E-6_dp, &
                          repeats=.FALSE.)
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_SCREENING", &
                          variants=s2a("EPS_SCREEN"), &
                          description="A threshold to determine which primitive 3-center integrals "// &
                          "are kept for contraction, as the latter operation can be "// &
                          "expensive (especially for large basis sets ). "// &
                          "If |(ab|c)| < EPS_SCREENNING, it is discarded.", &
                          default_r_val=1.0E-8_dp, &
                          repeats=.FALSE.)
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SCALE", &
                          variants=s2a("FRACTION"), &
                          description="Scaling of the exact exchange contribution.", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(subsubsection, keyword)
      CALL keyword_release(keyword)

      !The RI metric subsection
      CALL section_create(subsubsubsection, __LOCATION__, name="RI_METRIC", &
                          description="This subsection allows for the definition of an exchange "// &
                          "RI metric that is different from the main exchange potential. "// &
                          "By default (i.e. if this subsection is ignored), the "// &
                          "exchange kernel is computed in the V approximation: "// &
                          "(ab|ij) = (ab|P) V^-1 (Q|ij), where V = (P|Q). With a RI "// &
                          "metric, we have a 2 step RI involving the metric potential "// &
                          "for the 3-center integrals: "// &
                          "(ab|ij) = (ab!P) (P!Q)^-1 (Q|R) (R!S)^-1 (S!ij), where | "// &
                          "stands for the exchange potential and ! for the metric "// &
                          "potential. This allows for drastic screening of the "// &
                          "3-center integrals by selecting shorter range metric.", &
                          n_keywords=5, &
                          n_subsections=0, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Enables the use of a RI metric.", &
                          usage="&EXACT_EXCHANGE {logical}", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsubsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="POTENTIAL_TYPE", &
                          variants=s2a("OP", "OPERATOR", "POTENTIAL"), &
                          description="The type of operator used for the metric.", &
                          usage="OPERATOR {string}", &
                          repeats=.FALSE., &
                          default_i_val=do_potential_id, &
                          enum_c_vals=s2a("OVERLAP", "TRUNCATED", "SHORTRANGE"), &
                          enum_desc=s2a("Overlap operator (=identity)", &
                                        "Truncated Coulomb operator: 1/r if r < R_c, 0 otherwise ", &
                                        "Short range: erfc(omega*r)/r"), &
                          enum_i_vals=(/do_potential_id, do_potential_truncated, &
                                        do_potential_short/))
      CALL section_add_keyword(subsubsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_RADIUS", &
                          variants=s2a("R_C", "RC", "RANGE"), &
                          description="The cutoff radius (in Angstrom) for the truncated Coulomb operator.", &
                          usage="CUTOFF_RADIUS {double}", &
                          default_r_val=0.0_dp, &
                          repeats=.FALSE.)
      CALL section_add_keyword(subsubsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="T_C_G_DATA", &
                          description="Location of the file t_c_g.dat that contains the data for the "// &
                          "evaluation of the truncated gamma function ", &
                          usage="T_C_G_DATA {string}", &
                          default_c_val="t_c_g.dat")
      CALL section_add_keyword(subsubsubsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OMEGA", &
                          description="The range parameter for the short range operator (in 1/a0).", &
                          usage="OMEGA {double}", &
                          default_r_val=0.0_dp, &
                          repeats=.FALSE.)
      CALL section_add_keyword(subsubsubsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(subsubsection, subsubsubsection)
      CALL section_release(subsubsubsection)

      CALL section_add_subsection(subsection, subsubsection)
      CALL section_release(subsubsection)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
! End of Kernel subsection

      CALL section_create(subsection, __LOCATION__, "PRINT", "Controls the printing of information during "// &
                          "XAS TDP calculations", repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="SPECTRUM", &
                                       description="Controles the dumping of the XAS TDP spectrum in output files", &
                                       print_level=low_print_level, filename="", common_iter_levels=3)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, name="RESTART_WFN", &
                                       description="Controles the dumping of a MO restart file for a given "// &
                                       "excited state index. Only for K-edge RKS calculations. "// &
                                       "Can be repeated to get multiple *.wfn files at once.", &
                                       print_level=debug_print_level, filename="", common_iter_levels=1)
      CALL keyword_create(keyword, __LOCATION__, name="EXCITED_STATE_INDEX", variants=(/"INDEX"/), &
                          description="The index of the excited state that should be dumped", &
                          usage="INDEX {int}", default_i_val=1, repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL create_pdos_section(print_key)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "CUBES", &
                                       description="Controls the printing of the linear-response orbitals "// &
                                       "as *.cube files.", &
                                       print_level=high_print_level, common_iter_levels=1, &
                                       add_last=add_last_numeric, filename="")

      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUBES_LU_BOUNDS", &
                          variants=(/"CUBES_LU"/), &
                          description="The lower and upper index of the excited states to be printed as cube", &
                          usage="CUBES_LU_BOUNDS integer integer", &
                          n_var=2, default_i_vals=(/1, 0/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUBES_LIST", &
                          description="Indexes of the excited states to be printed as cube files "// &
                          "This keyword can be repeated several times "// &
                          "(useful if you have to specify many indexes).", &
                          usage="CUBES_LIST 1 2", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESTART", &
                                       description="Controls the dumping of LR-orbitals coefficients "// &
                                       "and corresponding excitation energies such that "// &
                                       "the program can be restarted for PDOS or CUBE "// &
                                       "printing without the heavy computing.", &
                                       print_level=high_print_level, filename="", common_iter_levels=3)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_xas_tdp_section

END MODULE input_cp2k_xas
