!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the dft section of the input
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_dft
   USE basis_set_types,                 ONLY: basis_sort_default,&
                                              basis_sort_zet
   USE bibliography,                    ONLY: &
        Andermatt2016, Andreussi2012, Avezac2005, Bengtsson1999, Blochl1995, Brelaz1979, &
        Fattebert2002, Guidon2010, Iannuzzi2006, Kunert2003, Marek2025, Merlot2014, Perdew1981, &
        VandeVondele2005b, Yin2017
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              high_print_level,&
                                              low_print_level,&
                                              medium_print_level,&
                                              silent_print_level
   USE cp_spline_utils,                 ONLY: pw_interp,&
                                              spline3_nopbc_interp,&
                                              spline3_pbc_interp
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE input_constants,                 ONLY: &
        admm1_type, admm2_type, admmp_type, admmq_type, admms_type, do_admm_aux_exch_func_bee, &
        do_admm_aux_exch_func_bee_libxc, do_admm_aux_exch_func_default, &
        do_admm_aux_exch_func_default_libxc, do_admm_aux_exch_func_none, &
        do_admm_aux_exch_func_opt, do_admm_aux_exch_func_opt_libxc, do_admm_aux_exch_func_pbex, &
        do_admm_aux_exch_func_pbex_libxc, do_admm_aux_exch_func_sx_libxc, &
        do_admm_basis_projection, do_admm_blocked_projection, do_admm_blocking_purify_full, &
        do_admm_charge_constrained_projection, do_admm_exch_scaling_merlot, &
        do_admm_exch_scaling_none, do_admm_purify_cauchy, do_admm_purify_cauchy_subspace, &
        do_admm_purify_mcweeny, do_admm_purify_mo_diag, do_admm_purify_mo_no_diag, &
        do_admm_purify_none, do_admm_purify_none_dm, do_arnoldi, do_bch, do_cn, do_em, do_etrs, &
        do_exact, do_pade, do_taylor, ehrenfest, gaussian, kg_color_dsatur, kg_color_greedy, &
        kg_tnadd_atomic, kg_tnadd_embed, kg_tnadd_embed_ri, kg_tnadd_none, no_admm_type, &
        numerical, plus_u_lowdin, plus_u_mulliken, plus_u_mulliken_charges, real_time_propagation, &
        rel_dkh, rel_none, rel_pot_erfc, rel_pot_full, rel_sczora_mp, rel_trans_atom, &
        rel_trans_full, rel_trans_molecule, rel_zora, rel_zora_full, rel_zora_mp, &
        rtp_bse_ham_g0w0, rtp_bse_ham_ks, rtp_method_bse, rtp_method_tddft, sccs_andreussi, &
        sccs_derivative_cd3, sccs_derivative_cd5, sccs_derivative_cd7, sccs_derivative_fft, &
        sccs_fattebert_gygi, sic_ad, sic_eo, sic_list_all, sic_list_unpaired, sic_mauri_spz, &
        sic_mauri_us, sic_none, slater, use_mom_ref_coac, use_mom_ref_com, use_mom_ref_user, &
        use_mom_ref_zero, use_restart_wfn, use_rt_restart, use_scf_wfn, weight_type_mass, &
        weight_type_unit
   USE input_cp2k_almo,                 ONLY: create_almo_scf_section
   USE input_cp2k_as,                   ONLY: create_active_space_section
   USE input_cp2k_ec,                   ONLY: create_ec_section
   USE input_cp2k_exstate,              ONLY: create_exstate_section
   USE input_cp2k_external,             ONLY: create_ext_den_section,&
                                              create_ext_pot_section,&
                                              create_ext_vxc_section
   USE input_cp2k_field,                ONLY: create_efield_section,&
                                              create_per_efield_section
   USE input_cp2k_harris,               ONLY: create_harris_section
   USE input_cp2k_kpoints,              ONLY: create_kpoints_section,&
                                              create_kpoint_set_section
   USE input_cp2k_loc,                  ONLY: create_localize_section
   USE input_cp2k_ls,                   ONLY: create_ls_scf_section
   USE input_cp2k_poisson,              ONLY: create_poisson_section
   USE input_cp2k_print_dft,            ONLY: create_print_dft_section
   USE input_cp2k_projection_rtp,       ONLY: create_projection_rtp_section
   USE input_cp2k_qs,                   ONLY: create_lrigpw_section,&
                                              create_qs_section
   USE input_cp2k_rsgrid,               ONLY: create_rsgrid_section
   USE input_cp2k_scf,                  ONLY: create_scf_section
   USE input_cp2k_smeagol,              ONLY: create_dft_smeagol_section
   USE input_cp2k_transport,            ONLY: create_transport_section
   USE input_cp2k_xas,                  ONLY: create_xas_section,&
                                              create_xas_tdp_section
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
   USE physcon,                         ONLY: evolt,&
                                              femtoseconds
   USE pw_spline_utils,                 ONLY: no_precond,&
                                              precond_spl3_1,&
                                              precond_spl3_2,&
                                              precond_spl3_3,&
                                              precond_spl3_aint,&
                                              precond_spl3_aint2
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_dft'

   PUBLIC :: create_dft_section
   PUBLIC :: create_bsse_section
   PUBLIC :: create_interp_section
   PUBLIC :: create_mgrid_section

CONTAINS

! **************************************************************************************************
!> \brief creates the dft section
!> \param section the section to be created
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE create_dft_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="DFT", &
                          description="Parameter needed by LCAO DFT programs", &
                          n_keywords=3, n_subsections=4, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="BASIS_SET_FILE_NAME", &
                          description="Name of the basis set file, may include a path", &
                          usage="BASIS_SET_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t, repeats=.TRUE., &
                          default_lc_val="BASIS_SET", n_var=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="POTENTIAL_FILE_NAME", &
                          description="Name of the pseudo potential file, may include a path", &
                          usage="POTENTIAL_FILE_NAME <FILENAME>", &
                          default_lc_val="POTENTIAL")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WFN_RESTART_FILE_NAME", &
                          variants=["RESTART_FILE_NAME"], &
                          description="Name of the wavefunction restart file, may include a path."// &
                          " If no file is specified, the default is to open the file as generated by the wfn restart print key.", &
                          usage="WFN_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="UKS", &
                          variants=s2a("UNRESTRICTED_KOHN_SHAM", &
                                       "LSD", &
                                       "SPIN_POLARIZED"), &
                          description="Requests a spin-polarized calculation using alpha "// &
                          "and beta orbitals, i.e. no spin restriction is applied", &
                          usage="LSD", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="ROKS", &
                          variants=["RESTRICTED_OPEN_KOHN_SHAM"], &
                          description="Requests a restricted open Kohn-Sham calculation", &
                          usage="ROKS", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="MULTIPLICITY", &
                          variants=["MULTIP"], &
                          description="Two times the total spin plus one. "// &
                          "Specify 3 for a triplet, 4 for a quartet, "// &
                          "and so on. Default is 1 (singlet) for an "// &
                          "even number and 2 (doublet) for an odd number "// &
                          "of electrons.", &
                          usage="MULTIPLICITY 3", &
                          default_i_val=0) ! this default value is just a flag to get the above
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="CHARGE", &
                          description="The total charge of the system", &
                          usage="CHARGE -1", &
                          default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="PLUS_U_METHOD", &
                          description="Method employed for the calculation of the DFT+U contribution", &
                          repeats=.FALSE., &
                          enum_c_vals=s2a("LOWDIN", "MULLIKEN", "MULLIKEN_CHARGES"), &
                          enum_i_vals=[plus_u_lowdin, plus_u_mulliken, plus_u_mulliken_charges], &
                          enum_desc=s2a("Method based on Lowdin population analysis "// &
                                        "(computationally expensive, since the diagonalization of the "// &
                                        "overlap matrix is required, but possibly more robust than Mulliken)", &
                                        "Method based on Mulliken population analysis using the net AO and "// &
                                        "overlap populations (computationally cheap method)", &
                                        "Method based on Mulliken gross orbital populations (GOP)"), &
                          n_var=1, &
                          default_i_val=plus_u_mulliken, &
                          usage="PLUS_U_METHOD Lowdin")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="RELAX_MULTIPLICITY", &
                          variants=["RELAX_MULTIP"], &
                          description="Tolerance in Hartrees. Do not enforce the occupation "// &
                          "of alpha and beta MOs due to the initially "// &
                          "defined multiplicity, but rather follow the Aufbau principle. "// &
                          "A value greater than zero activates this option. "// &
                          "If alpha/beta MOs differ in energy less than this tolerance, "// &
                          "then alpha-MO occupation is preferred even if it is higher "// &
                          "in energy (within the tolerance). "// &
                          "Such spin-symmetry broken (spin-polarized) occupation is used "// &
                          "as SCF input, which (is assumed to) bias the SCF "// &
                          "towards a spin-polarized solution. "// &
                          "Thus, larger tolerance increases chances of ending up "// &
                          "with spin-polarization. "// &
                          "This option is only valid for unrestricted (i.e. spin polarised) "// &
                          "Kohn-Sham (UKS) calculations. It also needs non-zero "// &
                          "[ADDED_MOS](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.ADDED_MOS) to actually affect the calculations, "// &
                          "which is why it is not expected to work with [OT](#CP2K_INPUT.FORCE_EVAL.DFT.SCF.OT) "// &
                          "and may raise errors when used with OT. "// &
                          "For more details see [this discussion](https://github.com/cp2k/cp2k/issues/4389).", &
                          usage="RELAX_MULTIPLICITY 0.00001", &
                          repeats=.FALSE., &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SUBCELLS", &
                          description="Read the grid size for subcell generation in the construction of "// &
                          "neighbor lists.", usage="SUBCELLS 1.5", &
                          n_var=1, default_r_val=2.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="AUTO_BASIS", &
                          description="Specify size of automatically generated auxiliary (RI) basis sets: "// &
                          "Options={small,medium,large,huge}", &
                          usage="AUTO_BASIS {basis_type} {basis_size}", &
                          type_of_var=char_t, repeats=.TRUE., n_var=-1, default_c_vals=["X", "X"])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SURFACE_DIPOLE_CORRECTION", &
                          variants=s2a("SURFACE_DIPOLE", &
                                       "SURF_DIP"), &
                          description="For slab calculations with asymmetric geometries, activate the correction of "// &
                          "the electrostatic potential with "// &
                          "by compensating for the surface dipole. Implemented only for slabs with normal "// &
                          "parallel to one Cartesian axis. The normal direction is given by the keyword SURF_DIP_DIR", &
                          usage="SURF_DIP", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE., &
                          citations=[Bengtsson1999])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SURF_DIP_DIR", &
                          description="Cartesian axis parallel to surface normal.", &
                          enum_c_vals=s2a("X", "Y", "Z"), &
                          enum_i_vals=[1, 2, 3], &
                          enum_desc=s2a("Along x", "Along y", "Along z"), &
                          n_var=1, &
                          default_i_val=3, &
                          usage="SURF_DIP_DIR Z")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SURF_DIP_POS", &
                          description="This keyword assigns an user defined position in Angstroms "// &
                          "in the direction normal to the surface (given by SURF_DIP_DIR). "// &
                          "The default value is -1.0_dp which appplies the correction at a position "// &
                          "that has minimum electron density on the grid.", &
                          usage="SURF_DIP_POS -1.0_dp", &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SURF_DIP_SWITCH", &
                          description="WARNING: Experimental feature under development that will help the "// &
                          "user to switch parameters to facilitate  SCF convergence. In its current form the "// &
                          "surface dipole correction is switched off if the calculation does not converge in "// &
                          "(0.5*MAX_SCF + 1) outer_scf steps. "// &
                          "The default value is .FALSE.", &
                          usage="SURF_DIP_SWITCH .TRUE.", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="CORE_CORR_DIP", &
                          description="If the total CORE_CORRECTION is non-zero and surface dipole "// &
                          "correction is switched on, presence of this keyword will adjust electron "// &
                          "density via MO occupation to reflect the total CORE_CORRECTION. "// &
                          "The default value is .FALSE.", &
                          usage="CORE_CORR_DIP .TRUE.", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SORT_BASIS", &
                          description="Sort basis sets according to a certain criterion. ", &
                          enum_c_vals=s2a("DEFAULT", "EXP"), &
                          enum_i_vals=[basis_sort_default, basis_sort_zet], &
                          enum_desc=s2a("don't sort", "sort w.r.t. exponent"), &
                          default_i_val=basis_sort_default, &
                          usage="SORT_BASIS EXP")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL create_scf_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_ls_scf_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_almo_scf_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_kg_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_harris_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_ec_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_exstate_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_admm_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_qs_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_mgrid_section(subsection, create_subsections=.TRUE.)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_xc_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_relativistic_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_sic_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_low_spin_roks_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_efield_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_per_efield_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_ext_pot_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_transport_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! ZMP sections to include the external density or v_xc potential
      CALL create_ext_den_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_ext_vxc_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_poisson_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_kpoints_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_kpoint_set_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_implicit_solv_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_density_fitting_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_xas_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_xas_tdp_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_localize_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_rtp_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_print_dft_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_sccs_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_active_space_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_dft_smeagol_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_hairy_probes_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_dft_section

! **************************************************************************************************
!> \brief Hairy Probe DFT Model
!> \param section ...
!> \author Margherita Buraschi
! **************************************************************************************************

   SUBROUTINE create_hairy_probes_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, &
                          name="HAIRY_PROBES", &
                          description="Sets up a Hairy Probe calculation. ", &
                          n_keywords=0, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="_SECTION_PARAMETERS_", &
                          description="Controls the activation of hairy probe", &
                          usage="&HAIRY_PROBES ON", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOM_IDS", &
                          description="Indexes of the atoms to which the probes are attached.", &
                          usage="ATOM_IDS <INTEGER> .. <INTEGER>", &
                          type_of_var=integer_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="T", &
                          description="Electronic temperature [K]", &
                          usage="T <REAL>", &
                          default_r_val=cp_unit_to_cp2k(value=300.0_dp, unit_str="K"), &
                          unit_str="K")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MU", &
                          description="Chemical potential of the electrons in the probes [eV] ", &
                          usage="MU <REAL>", &
                          default_r_val=cp_unit_to_cp2k(value=0.0_dp, unit_str="eV"), &
                          unit_str="eV")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ALPHA", &
                          description="Parameter for solution probes ", &
                          usage="ALPHA <REAL>", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_HP", &
                          description=" Tolerance for accuracy checks on occupation numbers "// &
                          "calculated using hair-probes. ", &
                          usage="EPS_HP <REAL>", &
                          default_r_val=1.0E-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
   END SUBROUTINE create_hairy_probes_section
!####################################################################################

! **************************************************************************************************
!> \brief Implicit Solvation Model
!> \param section ...
!> \author tlaino
! **************************************************************************************************
   SUBROUTINE create_implicit_solv_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      NULLIFY (keyword, subsection, print_key)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="SCRF", &
                          description="Adds an implicit solvation model to the DFT calculation."// &
                          " Know also as Self Consistent Reaction Field.", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_OUT", &
                          description="Value of the dielectric constant outside the sphere", &
                          usage="EPS_OUT <REAL>", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LMAX", &
                          description="Maximum value of L used in the multipole expansion", &
                          usage="LMAX <INTEGER>", &
                          default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_sphere_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "program_run_info", &
                                       description="Controls the printing basic info about the method", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_implicit_solv_section

! **************************************************************************************************
!> \brief Create Sphere cavity
!> \param section ...
!> \author tlaino
! **************************************************************************************************
   SUBROUTINE create_sphere_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      NULLIFY (keyword, subsection)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="SPHERE", &
                          description="Treats the implicit solvent environment like a sphere", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="RADIUS", &
                          description="Value of the spherical cavity in the dielectric medium", &
                          usage="RADIUS <REAL>", &
                          unit_str="angstrom", &
                          type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_center_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_sphere_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
!> \author tlaino
! **************************************************************************************************
   SUBROUTINE create_center_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="CENTER", &
                          description="Defines the center of the sphere.", &
                          n_keywords=0, n_subsections=0, repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="XYZ", &
                          description="Coordinates of the center of the sphere", &
                          usage="XYZ <REAL> <REAL> <REAL>", &
                          unit_str="angstrom", &
                          type_of_var=real_t, n_var=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOM_LIST", &
                          description="Defines a list of atoms to define the center of the sphere", &
                          usage="ATOM_LIST <INTEGER> .. <INTEGER>", &
                          type_of_var=integer_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WEIGHT_TYPE", &
                          description="Defines the weight used to define the center of the sphere"// &
                          " (if ATOM_LIST is provided)", &
                          usage="WEIGHT_TYPE (UNIT|MASS)", &
                          enum_c_vals=["UNIT", "MASS"], &
                          enum_i_vals=[weight_type_unit, weight_type_mass], &
                          default_i_val=weight_type_unit)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FIXED", &
                          description="Specify if the center of the sphere should be fixed or"// &
                          " allowed to move", &
                          usage="FIXED <LOGICAL>", &
                          default_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_center_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_admm_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="AUXILIARY_DENSITY_MATRIX_METHOD", &
                          description="Parameters needed for the ADMM method.", &
                          n_keywords=1, n_subsections=1, repeats=.FALSE., &
                          citations=[Guidon2010])

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="ADMM_TYPE", &
         description="Type of ADMM (sort name) as refered in literature. "// &
         "This sets values for METHOD, ADMM_PURIFICATION_METHOD, and EXCH_SCALING_MODEL", &
         enum_c_vals=s2a("NONE", "ADMM1", "ADMM2", "ADMMS", "ADMMP", "ADMMQ"), &
         enum_desc=s2a("No short name is used, use specific definitions (default)", &
                       "ADMM1 method from Guidon2010", &
                       "ADMM2 method from Guidon2010", &
                       "ADMMS method from Merlot2014", &
                       "ADMMP method from Merlot2014", &
                       "ADMMQ method from Merlot2014"), &
         enum_i_vals=[no_admm_type, admm1_type, admm2_type, admms_type, admmp_type, admmq_type], &
         default_i_val=no_admm_type, &
         citations=[Guidon2010, Merlot2014])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="ADMM_PURIFICATION_METHOD", &
         description="Method that shall be used for wavefunction fitting. Use MO_DIAG for MD.", &
         enum_c_vals=s2a("NONE", "CAUCHY", "CAUCHY_SUBSPACE", "MO_DIAG", "MO_NO_DIAG", "MCWEENY", "NONE_DM"), &
         enum_i_vals=[do_admm_purify_none, do_admm_purify_cauchy, do_admm_purify_cauchy_subspace, &
                      do_admm_purify_mo_diag, do_admm_purify_mo_no_diag, &
                      do_admm_purify_mcweeny, do_admm_purify_none_dm], &
         enum_desc=s2a("Do not apply any purification", &
                       "Perform purification via general Cauchy representation", &
                       "Perform purification via Cauchy representation in occupied subspace", &
                       "Calculate MO derivatives via Cauchy representation by diagonalization", &
                       "Calculate MO derivatives via Cauchy representation by inversion", &
                       "Perform original McWeeny purification via matrix multiplications", &
                       "Do not apply any purification, works directly with density matrix"), &
         default_i_val=do_admm_purify_mo_diag)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="METHOD", &
         description="Method that shall be used for wavefunction fitting. Use BASIS_PROJECTION for MD.", &
         enum_c_vals=s2a("BASIS_PROJECTION", "BLOCKED_PROJECTION_PURIFY_FULL", "BLOCKED_PROJECTION", &
                         "CHARGE_CONSTRAINED_PROJECTION"), &
         enum_i_vals=[do_admm_basis_projection, do_admm_blocking_purify_full, do_admm_blocked_projection, &
                      do_admm_charge_constrained_projection], &
         enum_desc=s2a("Construct auxiliary density matrix from auxiliary basis.", &
                       "Construct auxiliary density from a blocked Fock matrix,"// &
                       " but use the original matrix for purification.", &
                       "Construct auxiliary density from a blocked Fock matrix.", &
                       "Construct auxiliary density from auxiliary basis enforcing charge constrain."), &
         default_i_val=do_admm_basis_projection)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="EXCH_SCALING_MODEL", &
         description="Scaling of the exchange correction calculated by the auxiliary density matrix.", &
         enum_c_vals=s2a("NONE", "MERLOT"), &
         enum_i_vals=[do_admm_exch_scaling_none, do_admm_exch_scaling_merlot], &
         enum_desc=s2a("No scaling is enabled, refers to methods ADMM1, ADMM2 or ADMMQ.", &
                       "Exchange scaling according to Merlot (2014)"), &
         default_i_val=do_admm_exch_scaling_none)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="EXCH_CORRECTION_FUNC", &
         description="Exchange functional which is used for the ADMM correction. "// &
         "LibXC implementations require linking with LibXC", &
         enum_c_vals=s2a("DEFAULT", "PBEX", "NONE", "OPTX", "BECKE88X", &
                         "PBEX_LIBXC", "BECKE88X_LIBXC", "OPTX_LIBXC", "DEFAULT_LIBXC", "LDA_X_LIBXC"), &
         enum_i_vals=[do_admm_aux_exch_func_default, do_admm_aux_exch_func_pbex, &
                      do_admm_aux_exch_func_none, do_admm_aux_exch_func_opt, do_admm_aux_exch_func_bee, &
                      do_admm_aux_exch_func_pbex_libxc, do_admm_aux_exch_func_bee_libxc, &
                      do_admm_aux_exch_func_opt_libxc, do_admm_aux_exch_func_default_libxc, &
                      do_admm_aux_exch_func_sx_libxc], &
         enum_desc=s2a("Use PBE-based corrections according to the chosen interaction operator.", &
                       "Use PBEX functional for exchange correction.", &
                       "No correction: X(D)-x(d)-> 0.", &
                       "Use OPTX functional for exchange correction.", &
                       "Use Becke88X functional for exchange correction.", &
                       "Use PBEX functional (LibXC implementation) for exchange correction.", &
                       "Use Becke88X functional (LibXC implementation) for exchange correction.", &
                       "Use OPTX functional (LibXC implementation) for exchange correction.", &
                       "Use PBE-based corrections (LibXC where possible) to the chosen interaction operator.", &
                       "Use Slater X functional (LibXC where possible) for exchange correction."), &
         default_i_val=do_admm_aux_exch_func_default)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="optx_a1", &
                          description="OPTX a1 coefficient", &
                          default_r_val=1.05151_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="optx_a2", &
                          description="OPTX a2 coefficient", &
                          default_r_val=1.43169_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="optx_gamma", &
                          description="OPTX gamma coefficient", &
                          default_r_val=0.006_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BLOCK_LIST", &
                          description="Specifies a list of atoms.", &
                          usage="BLOCK_LIST {integer} {integer} .. {integer}", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          description="Define accuracy of DBCSR operations", &
                          usage="EPS_FILTER", default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_admm_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_density_fitting_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      NULLIFY (keyword, print_key)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="DENSITY_FITTING", &
                          description="Setup parameters for density fitting (Bloechl charges or density derived "// &
                          "atomic point charges (DDAPC) charges)", &
                          n_keywords=7, n_subsections=0, repeats=.FALSE., &
                          citations=[Blochl1995])

      CALL keyword_create(keyword, __LOCATION__, name="NUM_GAUSS", &
                          description="Specifies the numbers of gaussian used to fit the QM density for each atomic site.", &
                          usage="NUM_GAUSS {integer}", &
                          n_var=1, type_of_var=integer_t, default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PFACTOR", &
                          description="Specifies the progression factor for the gaussian exponent for each atomic site.", &
                          usage="PFACTOR {real}", &
                          n_var=1, type_of_var=real_t, default_r_val=1.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MIN_RADIUS", &
                          description="Specifies the smallest radius of the gaussian used in the fit. All other radius are"// &
                          " obtained with the progression factor.", &
                          usage="MIN_RADIUS {real}", &
                          unit_str="angstrom", n_var=1, type_of_var=real_t, default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RADII", &
                          description="Specifies all the radius of the gaussian used in the fit for each atomic site. The use"// &
                          " of this keyword disables all other keywords of this section.", &
                          usage="RADII {real} {real} .. {real}", &
                          unit_str="angstrom", n_var=-1, type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GCUT", &
                          description="Cutoff for charge fit in G-space.", &
                          usage="GCUT {real}", &
                          n_var=1, type_of_var=real_t, default_r_val=SQRT(6.0_dp))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "program_run_info", &
                                       description="Controls the printing of basic information during the run", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")

      CALL keyword_create(keyword, __LOCATION__, name="CONDITION_NUMBER", &
                          description="Prints information regarding the condition numbers of the A matrix (to be inverted)", &
                          usage="CONDITION_NUMBER <LOGICAL>", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_density_fitting_section

! **************************************************************************************************
!> \brief creates the input section for the relativistic part
!> \param section the section to create
!> \author jens
! **************************************************************************************************
   SUBROUTINE create_relativistic_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="relativistic", &
                          description="parameters needed and setup for relativistic calculations", &
                          n_keywords=5, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="method", &
                          description="type of relativistic correction used", &
                          usage="method (NONE|DKH|ZORA)", default_i_val=rel_none, &
                          enum_c_vals=s2a("NONE", "DKH", "ZORA"), &
                          enum_i_vals=[rel_none, rel_dkh, rel_zora], &
                          enum_desc=s2a("Use no relativistic correction", &
                                        "Use Douglas-Kroll-Hess method", &
                                        "Use ZORA method"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DKH_order", &
                          description="The order of the DKH transformation ", &
                          usage="DKH_order 2", default_i_val=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ZORA_type", &
                          description="Type of ZORA method to be used", &
                          usage="ZORA_type scMP", default_i_val=rel_zora_full, &
                          enum_c_vals=s2a("FULL", "MP", "scMP"), &
                          enum_desc=s2a("Full ZORA method (not implemented)", &
                                        "ZORA with atomic model potential", &
                                        "Scaled ZORA with atomic model potential"), &
                          enum_i_vals=[rel_zora_full, rel_zora_mp, rel_sczora_mp])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="transformation", &
                          description="Type of DKH transformation", &
                          usage="transformation (FULL|MOLECULE|ATOM)", default_i_val=rel_trans_atom, &
                          enum_c_vals=s2a("FULL", "MOLECULE", "ATOM"), &
                          enum_i_vals=[rel_trans_full, rel_trans_molecule, rel_trans_atom], &
                          enum_desc=s2a("Use full matrix transformation", &
                                        "Use transformation blocked by molecule", &
                                        "Use atomic blocks"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="z_cutoff", &
                          description="The minimal atomic number considered for atom transformation", &
                          usage="z_cutoff 50", default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="potential", &
                          description="External potential used in DKH transformation, full 1/r or erfc(r)/r", &
                          usage="POTENTIAL {FULL,ERFC}", default_i_val=rel_pot_erfc, &
                          enum_c_vals=s2a("FULL", "ERFC"), &
                          enum_i_vals=[rel_pot_full, rel_pot_erfc])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_relativistic_section

! **************************************************************************************************
!> \brief creates the KG section
!> \param section ...
!> \author Martin Haeufel [2012.07]
! **************************************************************************************************
   SUBROUTINE create_kg_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="KG_METHOD", &
                          description="Specifies the parameters for a Kim-Gordon-like partitioning"// &
                          " into molecular subunits", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE., &
                          citations=[Iannuzzi2006, Brelaz1979, Andermatt2016])

      NULLIFY (keyword, subsection, print_key)

      ! add a XC section
      CALL create_xc_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! add LRI section
      CALL create_lrigpw_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, __LOCATION__, name="COLORING_METHOD", &
                          description="Which algorithm to use for coloring.", &
                          usage="COLORING_METHOD GREEDY", &
                          default_i_val=kg_color_dsatur, &
                          enum_c_vals=s2a("DSATUR", "GREEDY"), &
                          enum_desc=s2a("Maximum degree of saturation, relatively accurate", &
                                        "Greedy, fast coloring, less accurate"), &
                          enum_i_vals=[kg_color_dsatur, kg_color_greedy])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TNADD_METHOD", &
                          description="Algorithm to use for the calculation of the nonadditive kinetic energy.", &
                          usage="TNADD_METHOD ATOMIC", &
                          default_i_val=kg_tnadd_embed, &
                          enum_c_vals=s2a("EMBEDDING", "RI_EMBEDDING", "ATOMIC", "NONE"), &
                          enum_desc=s2a("Use full embedding potential (see Iannuzzi et al)", &
                                        "Use full embedding potential with RI density fitting", &
                                        "Use sum of atomic model potentials", &
                                        "Do not use kinetic energy embedding"), &
                          enum_i_vals=[kg_tnadd_embed, kg_tnadd_embed_ri, kg_tnadd_atomic, kg_tnadd_none])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INTEGRATION_GRID", &
                          description="Grid [small,medium,large,huge]to be used for the TNADD integration.", &
                          usage="INTEGRATION_GRID MEDIUM", &
                          default_c_val="MEDIUM")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL section_create(subsection, __LOCATION__, name="PRINT", &
                          description="Print section", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "NEIGHBOR_LISTS", &
                                       description="Controls the printing of the neighbor lists.", &
                                       print_level=low_print_level, filename="__STD_OUT__", unit_str="angstrom")

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SAB_ORB_FULL", &
                          description="Activates the printing of the full orbital "// &
                          "orbital neighbor lists.", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SAB_ORB_MOLECULAR", &
                          description="Activates the printing of the orbital "// &
                          "orbital neighbor lists for molecular subsets.", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SAC_KIN", &
                          description="Activates the printing of the orbital "// &
                          "atomic potential neighbor list.", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_kg_section

! **************************************************************************************************
!> \brief Create the BSSE section for counterpoise correction
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_bsse_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="BSSE", &
                          description="This section is used to set up the BSSE calculation. "// &
                          "It also requires that for each atomic kind X a kind X_ghost is present, "// &
                          "with the GHOST keyword specified, in addition to the other required fields.", &
                          n_keywords=3, n_subsections=1, repeats=.FALSE.)

      NULLIFY (keyword, subsection)
      ! FRAGMENT SECTION
      CALL section_create(subsection, __LOCATION__, name="FRAGMENT", &
                          description="Specify the atom number belonging to this fragment.", &
                          n_keywords=2, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="LIST", &
                          description="Specifies a list of atoms.", &
                          usage="LIST {integer} {integer} .. {integer}", &
                          repeats=.TRUE., n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! CONFIGURATION SECTION
      CALL section_create(subsection, __LOCATION__, name="CONFIGURATION", &
                          description="Specify additional parameters for the combinatorial configurations. "// &
                          "Use this section to manually specify charge and multiplicity of the fragments "// &
                          "and their combinations.", &
                          n_keywords=2, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="GLB_CONF", &
                          description="Specifies the global configuration using 1 or 0 for each fragment. "// &
                          "1 specifies the respective fragment as used, 0 as unused.", &
                          usage="GLB_CONF {integer} {integer} .. {integer}", &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SUB_CONF", &
                          description="Specifies the subconfiguration using 1 or 0 belonging to the global configuration. "// &
                          "1 specifies the respective fragment as real, 0 as ghost.", &
                          usage="SUB_CONF {integer} {integer} .. {integer}", &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="MULTIPLICITY", &
                          variants=["MULTIP"], &
                          description="Specify for each fragment the multiplicity. Two times the total spin plus one. "// &
                          "Specify 3 for a triplet, 4 for a quartet,and so on. Default is 1 (singlet) for an "// &
                          "even number and 2 (doublet) for an odd number of electrons.", &
                          usage="MULTIPLICITY 3", &
                          default_i_val=0) ! this default value is just a flag to get the above
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHARGE", &
                          description="The total charge for each fragment.", &
                          usage="CHARGE -1", &
                          default_i_val=0)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, name="FRAGMENT_ENERGIES", &
                          description="This section contains the energies of the fragments already"// &
                          " computed. It is useful as a summary and specifically for restarting BSSE runs.", &
                          n_keywords=2, n_subsections=0, repeats=.TRUE.)
      CALL keyword_create(keyword, __LOCATION__, name="_DEFAULT_KEYWORD_", &
                          description="The energy computed for each fragment", repeats=.TRUE., &
                          usage="{REAL}", type_of_var=real_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_print_bsse_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_bsse_section

! **************************************************************************************************
!> \brief Create the print bsse section
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_print_bsse_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="print", &
                          description="Section of possible print options in BSSE code.", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "PROGRAM_RUN_INFO", &
                                       description="Controls the printing of information regarding the run.", &
                                       print_level=low_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESTART", &
                                       description="Controls the dumping of the restart file during BSSE runs. "// &
                                       "By default the restart is updated after each configuration calculation is "// &
                                       "completed.", &
                                       print_level=silent_print_level, common_iter_levels=0, &
                                       add_last=add_last_numeric, filename="")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_print_bsse_section

! **************************************************************************************************
!> \brief input section for optional parameters for RIGPW
!> \param section the section to create
!> \author JGH [06.2017]
! **************************************************************************************************
   SUBROUTINE create_rigpw_section(section)
      TYPE(section_type), POINTER                        :: section

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="RIGPW", &
                          description="This section specifies optional parameters for RIGPW.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

!     CALL keyword_create(keyword, __LOCATION__, name="RI_OVERLAP_MATRIX", &
!                         description="Specifies whether to calculate the inverse or the "// &
!                         "pseudoinverse of the overlap matrix of the auxiliary "// &
!                         "basis set. Calculating the pseudoinverse is necessary "// &
!                         "for very large auxiliary basis sets, but more expensive. "// &
!                         "Using the pseudoinverse, consistent forces are not "// &
!                         "guaranteed yet.", &
!                         usage="RI_OVERLAP_MATRIX INVERSE", &
!                         enum_c_vals=s2a("INVERSE", "PSEUDO_INVERSE_SVD", "PSEUDO_INVERSE_DIAG", &
!                                         "AUTOSELECT"), &
!                         enum_desc=s2a("Calculate inverse of the overlap matrix.", &
!                                       "Calculate the pseuodinverse of the overlap matrix "// &
!                                       "using singular value decomposition.", &
!                                       "Calculate the pseudoinverse of the overlap matrix "// &
!                                       "by prior diagonalization.", &
!                                       "Choose automatically for each pair whether to "// &
!                                       "calculate the inverse or pseudoinverse based on the "// &
!                                       "condition number of the overlap matrix for each pair. "// &
!                                       "Calculating the pseudoinverse is much more expensive."), &
!                         enum_i_vals=(/do_lri_inv, do_lri_pseudoinv_svd, &
!                                       do_lri_pseudoinv_diag, do_lri_inv_auto/), &
!                         default_i_val=do_lri_inv)
!     CALL section_add_keyword(section, keyword)
!     CALL keyword_release(keyword)

   END SUBROUTINE create_rigpw_section

! **************************************************************************************************
!> \brief creates the multigrid
!> \param section             input section to create
!> \param create_subsections  indicates whether or not subsections INTERPOLATOR and RS_GRID
!>                            should be created
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE create_mgrid_section(section, create_subsections)
      TYPE(section_type), POINTER                        :: section
      LOGICAL, INTENT(in)                                :: create_subsections

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="mgrid", &
                          description="multigrid information", &
                          n_keywords=5, n_subsections=1, repeats=.FALSE.)
      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="NGRIDS", &
                          description="The number of multigrids to use", &
                          usage="ngrids 1", default_i_val=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="cutoff", &
                          description="The cutoff of the finest grid level. Default value for "// &
                          "SE or DFTB calculation is 1.0 [Ry].", &
                          usage="cutoff 300", default_r_val=cp_unit_to_cp2k(value=280.0_dp, &
                                                                            unit_str="Ry"), n_var=1, unit_str="Ry")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="progression_factor", &
                          description="Factor used to find the cutoff of the multigrids that"// &
                          " where not given explicitly", &
                          usage="progression_factor <integer>", default_r_val=3._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="commensurate", &
                          description="If the grids should be commensurate. If true overrides "// &
                          "the progression factor and the cutoffs of the sub grids", &
                          usage="commensurate", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="realspace", &
                          description="If both rho and rho_gspace are needed ", &
                          usage="realspace", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REL_CUTOFF", &
                          variants=["RELATIVE_CUTOFF"], &
                          description="Determines the grid at which a Gaussian is mapped,"// &
                          " giving the cutoff used for a gaussian with alpha=1."// &
                          " A value 50+-10Ry might be required for highly accurate results,"// &
                          " Or for simulations with a variable cell."// &
                          " Versions prior to 2.3 used a default of 30Ry.", &
                          usage="RELATIVE_CUTOFF real", default_r_val=20.0_dp, &
                          unit_str="Ry")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIGRID_SET", &
                          description="Activate a manual setting of the multigrids", &
                          usage="MULTIGRID_SET", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SKIP_LOAD_BALANCE_DISTRIBUTED", &
                          description="Skips load balancing on distributed multigrids.  "// &
                          "Memory usage is O(p) so may be used "// &
                          "for all but the very largest runs.", &
                          usage="SKIP_LOAD_BALANCE_DISTRIBUTED", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
!    CALL keyword_create(keyword, __LOCATION__, name="SKIP_LOAD_BALANCE_DISTRIBUTED",&
!         description="Skip load balancing on distributed multigrids, which might be memory intensive."//&
!                     "If not explicitly specified, runs using more than 1024 MPI tasks will default to .TRUE.",&
!         usage="SKIP_LOAD_BALANCE_DISTRIBUTED", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)

      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIGRID_CUTOFF", &
                          variants=["CUTOFF_LIST"], &
                          description="List of cutoff values to set up multigrids manually", &
                          usage="MULTIGRID_CUTOFF 200.0 100.0 ", &
                          n_var=-1, &
                          type_of_var=real_t, &
                          unit_str="Ry")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      IF (create_subsections) THEN
         NULLIFY (subsection)
         CALL create_rsgrid_section(subsection)
         CALL section_add_subsection(section, subsection)
         CALL section_release(subsection)

         NULLIFY (subsection)
         CALL create_interp_section(subsection)
         CALL section_add_subsection(section, subsection)
         CALL section_release(subsection)
      END IF
   END SUBROUTINE create_mgrid_section

! **************************************************************************************************
!> \brief creates the interpolation section
!> \param section ...
!> \author tlaino
! **************************************************************************************************
   SUBROUTINE create_interp_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="interpolator", &
                          description="kind of interpolation used between the multigrids", &
                          n_keywords=5, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword, print_key)

      CALL keyword_create(keyword, __LOCATION__, name="kind", &
                          description="the interpolator to use", &
                          usage="kind spline3", &
                          default_i_val=pw_interp, &
                          enum_c_vals=s2a("pw", "spline3_nopbc", "spline3"), &
                          enum_i_vals=[pw_interp, &
                                       spline3_nopbc_interp, spline3_pbc_interp])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="safe_computation", &
                          description="if a non unrolled calculation is to be performed in parallel", &
                          usage="safe_computation OFF", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="aint_precond", &
                          description="the approximate inverse to use to get the starting point"// &
                          " for the linear solver of the spline3 methods", &
                          usage="aint_precond copy", &
                          default_i_val=precond_spl3_aint, &
                          enum_c_vals=s2a("copy", "spl3_nopbc_aint1", "spl3_nopbc_aint2", &
                                          "spl3_nopbc_precond1", "spl3_nopbc_precond2", "spl3_nopbc_precond3"), &
                          enum_i_vals=[no_precond, precond_spl3_aint, precond_spl3_aint2, &
                                       precond_spl3_1, precond_spl3_2, precond_spl3_3])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="precond", &
                          description="The preconditioner used"// &
                          " for the linear solver of the spline3 methods", &
                          usage="PRECOND copy", &
                          default_i_val=precond_spl3_3, &
                          enum_c_vals=s2a("copy", "spl3_nopbc_aint1", "spl3_nopbc_aint2", &
                                          "spl3_nopbc_precond1", "spl3_nopbc_precond2", "spl3_nopbc_precond3"), &
                          enum_i_vals=[no_precond, precond_spl3_aint, precond_spl3_aint2, &
                                       precond_spl3_1, precond_spl3_2, precond_spl3_3])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="eps_x", &
                          description="accuracy on the solution for spline3 the interpolators", &
                          usage="eps_x 1.e-15", default_r_val=1.e-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="eps_r", &
                          description="accuracy on the residual for spline3 the interpolators", &
                          usage="eps_r 1.e-15", default_r_val=1.e-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="max_iter", &
                          variants=['maxiter'], &
                          description="the maximum number of iterations", &
                          usage="max_iter 200", default_i_val=100)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "conv_info", &
                                       description="if convergence information about the linear solver"// &
                                       " of the spline methods should be printed", &
                                       print_level=medium_print_level, each_iter_names=s2a("SPLINE_FIND_COEFFS"), &
                                       each_iter_values=[10], filename="__STD_OUT__", &
                                       add_last=add_last_numeric)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_interp_section

! **************************************************************************************************
!> \brief creates the sic (self interaction correction) section
!> \param section ...
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE create_sic_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="sic", &
                          description="parameters for the self interaction correction", &
                          n_keywords=6, n_subsections=0, repeats=.FALSE., &
                          citations=[VandeVondele2005b, Perdew1981, Avezac2005])

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SIC_SCALING_A", &
                          description="Scaling of the coulomb term in sic [experimental]", &
                          usage="SIC_SCALING_A 0.5", &
                          citations=[VandeVondele2005b], &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SIC_SCALING_B", &
                          description="Scaling of the xc term in sic [experimental]", &
                          usage="SIC_SCALING_B 0.5", &
                          citations=[VandeVondele2005b], &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SIC_METHOD", &
                          description="Method used to remove the self interaction", &
                          usage="SIC_METHOD MAURI_US", &
                          default_i_val=sic_none, &
                          enum_c_vals=s2a("NONE", "MAURI_US", "MAURI_SPZ", "AD", "EXPLICIT_ORBITALS"), &
                          enum_i_vals=[sic_none, sic_mauri_us, sic_mauri_spz, sic_ad, sic_eo], &
                          enum_desc=s2a("Do not apply a sic correction", &
                                        "Employ a (scaled) correction proposed by Mauri and co-workers"// &
                                        " on the spin density / doublet unpaired orbital", &
                                        "Employ a (scaled) Perdew-Zunger expression"// &
                                        " on the spin density / doublet unpaired orbital", &
                                        "The average density correction", &
                                        "(scaled) Perdew-Zunger correction explicitly on a set of orbitals."), &
                          citations=[VandeVondele2005b, Perdew1981, Avezac2005])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORBITAL_SET", &
                          description="Type of orbitals treated with the SIC", &
                          usage="ORBITAL_SET ALL", &
                          default_i_val=sic_list_unpaired, &
                          enum_c_vals=s2a("UNPAIRED", "ALL"), &
                          enum_desc=s2a("correction for the unpaired orbitals only, requires a restricted open shell calculation", &
                                        "correction for all orbitals, requires a LSD or ROKS calculation"), &
                          enum_i_vals=[sic_list_unpaired, sic_list_all])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_sic_section

! **************************************************************************************************
!> \brief creates the low spin roks section
!> \param section ...
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE create_low_spin_roks_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="LOW_SPIN_ROKS", &
                          description="Specify the details of the low spin ROKS method. "// &
                          "In particular, one can specify various terms added to the energy of the high spin roks configuration"// &
                          " with a energy scaling factor, and a prescription of the spin state.", &
                          n_keywords=6, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ENERGY_SCALING", &
                          description="The scaling factors for each term added to the total energy. "// &
                          "This list should contain one number for each term added to the total energy.", &
                          usage="ENERGY_SCALING 1.0 -1.0 ", &
                          n_var=-1, type_of_var=real_t, repeats=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create( &
         keyword, __LOCATION__, name="SPIN_CONFIGURATION", &
         description="For each singly occupied orbital, specify if this should be an alpha (=1) or a beta (=2) orbital. "// &
         "This keyword should be repeated, each repetition corresponding to an additional term.", &
         usage="SPIN_CONFIGURATION 1 2", &
         n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_low_spin_roks_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_rtp_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, print_section, subsection

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="REAL_TIME_PROPAGATION", &
                          description="Parameters needed to set up the real time propagation"// &
                          " for the electron dynamics. This currently works only in the NVE ensemble.", &
                          n_keywords=4, n_subsections=4, repeats=.FALSE., &
                          citations=[Kunert2003, Andermatt2016])

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximal number of iterations for the self consistent propagator loop.", &
                          usage="MAX_ITER 10", &
                          default_i_val=10)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ITER", &
                          description="Convergence criterion for the self consistent propagator loop.", &
                          usage="EPS_ITER 1.0E-5", &
                          default_r_val=1.0E-7_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ASPC_ORDER", &
                          description="Speciefies how many steps will be used for extrapolation. "// &
                          "One will be always used which is means X(t+dt)=X(t)", &
                          usage="ASPC_ORDER 3", &
                          default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAT_EXP", &
                          description="Which method should be used to calculate the exponential"// &
                          " in the propagator. It is recommended to use BCH when employing density_propagation "// &
                          "and ARNOLDI otherwise.", &
                          usage="MAT_EXP TAYLOR", default_i_val=do_arnoldi, &
                          enum_c_vals=s2a("TAYLOR", "PADE", "ARNOLDI", "BCH", "EXACT"), &
                          enum_i_vals=[do_taylor, do_pade, do_arnoldi, do_bch, do_exact], &
                          enum_desc=s2a("exponential is evaluated using scaling and squaring in combination"// &
                                        " with a taylor expansion of the exponential.", &
                                        "uses scaling and squaring together with the pade approximation", &
                                        "uses arnoldi subspace algorithm to compute exp(H)*MO directly, can't be used in "// &
                                        "combination with Crank Nicholson or density propagation", &
                                        "Uses a Baker-Campbell-Hausdorff expansion to propagate the density matrix,"// &
                                        " only works for density propagation", &
                                        "Uses diagonalisation of the exponent matrices to determine the "// &
                                        "matrix exponential exactly. Only implemented for GWBSE."))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DENSITY_PROPAGATION", &
                          description="The density matrix is propagated instead of the molecular orbitals. "// &
                          "This can allow a linear scaling simulation. The density matrix is filtered with "// &
                          "the threshold based on the EPS_FILTER keyword from the LS_SCF section", &
                          usage="DENSITY_PROPAGATION .TRUE.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SC_CHECK_START", &
                          description="Speciefies how many iteration steps will be done without "// &
                          "a check for self consistency. Can save some time in big calculations.", &
                          usage="SC_CHECK_START 3", &
                          default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EXP_ACCURACY", &
                          description="Accuracy for the taylor and pade approximation. "// &
                          "This is only an upper bound bound since the norm used for the guess "// &
                          "is an upper bound for the needed one.", &
                          usage="EXP_ACCURACY 1.0E-6", &
                          default_r_val=1.0E-9_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PROPAGATOR", &
                          description="Which propagator should be used for the orbitals", &
                          usage="PROPAGATOR ETRS", default_i_val=do_etrs, &
                          enum_c_vals=s2a("ETRS", "CN", "EM"), &
                          enum_i_vals=[do_etrs, do_cn, do_em], &
                          enum_desc=s2a("enforced time reversible symmetry", &
                                        "Crank Nicholson propagator", &
                                        "Exponential midpoint propagator"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INITIAL_WFN", &
                          description="Controls the initial WFN used for propagation. "// &
                          "Note that some energy contributions may not be "// &
                          "initialized in the restart cases, for instance "// &
                          "electronic entropy energy in the case of smearing.", &
                          usage="INITIAL_WFN SCF_WFN", default_i_val=use_scf_wfn, &
                          enum_c_vals=s2a("SCF_WFN", "RESTART_WFN", "RT_RESTART"), &
                          enum_i_vals=[use_scf_wfn, use_restart_wfn, use_rt_restart], &
                          enum_desc=s2a("An SCF run is performed to get the initial state.", &
                                        "A wavefunction from a previous SCF is propagated. Especially useful,"// &
                                        " if electronic constraints or restraints are used in the previous calculation, "// &
                                        "since these do not work in the rtp scheme.", &
                                        "use the wavefunction of a real time propagation/ehrenfest run"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="APPLY_WFN_MIX_INIT_RESTART", &
                          description="If set to True and in the case of INITIAL_WFN=RESTART_WFN, call the "// &
                          "DFT%PRINT%WFN_MIX section to mix the read initial wfn. The starting wave-function of the "// &
                          "RTP will be the mixed one. Setting this to True without a defined WFN_MIX section will "// &
                          "not do anything as defining a WFN_MIX section without this keyword for RTP run with "// &
                          "INITIAL_WFN=RESTART_WFN. Note that if INITIAL_WFN=SCF_WFN, this keyword is not needed to "// &
                          "apply the mixing defined in the WFN_MIX section. Default is False.", &
                          usage="APPLY_WFN_MIX_INIT_RESTART", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="APPLY_DELTA_PULSE", &
                          description="Applies a delta kick to the initial wfn (only RTP for now - the EMD"// &
                          " case is not yet implemented). Only work for INITIAL_WFN=SCF_WFN", &
                          usage="APPLY_DELTA_PULSE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="APPLY_DELTA_PULSE_MAG", &
                          description="Applies a magnetic delta kick to the initial wfn (only RTP for now - the EMD"// &
                          " case is not yet implemented). Only work for INITIAL_WFN=SCF_WFN", &
                          usage="APPLY_DELTA_PULSE_MAG", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VELOCITY_GAUGE", &
                          description="Perform propagation in the velocity gauge using the explicit vector potential"// &
                          " only a constant vector potential as of now (corresonding to a delta-pulse)."// &
                          " uses DELTA_PULSE_SCALE and DELTA_PULSE_DIRECTION to define the vector potential", &
                          usage="VELOCITY_GAUGE T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAUGE_ORIG", &
                          description="Define gauge origin for magnetic perturbation", &
                          usage="GAUGE_ORIG COM", &
                          enum_c_vals=s2a("COM", "COAC", "USER_DEFINED", "ZERO"), &
                          enum_desc=s2a("Use Center of Mass", &
                                        "Use Center of Atomic Charges", &
                                        "Use User Defined Point (Keyword:REF_POINT)", &
                                        "Use Origin of Coordinate System"), &
                          enum_i_vals=[use_mom_ref_com, &
                                       use_mom_ref_coac, &
                                       use_mom_ref_user, &
                                       use_mom_ref_zero], &
                          default_i_val=use_mom_ref_com)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAUGE_ORIG_MANUAL", &
                          description="Manually defined gauge origin for magnetic perturbation [in Bohr!]", &
                          usage="GAUGE_ORIG_MANUAL x y z", &
                          repeats=.FALSE., &
                          n_var=3, default_r_vals=[0._dp, 0._dp, 0._dp], &
                          type_of_var=real_t, &
                          unit_str='bohr')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VG_COM_NL", &
                          description="apply gauge transformed non-local potential term"// &
                          " only affects VELOCITY_GAUGE=.TRUE.", &
                          usage="VG_COM_NL T", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COM_NL", &
                          description="Include non-local commutator for periodic delta pulse."// &
                          " only affects PERIODIC=.TRUE.", &
                          usage="COM_NL", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LEN_REP", &
                          description="Use length representation delta pulse (in conjunction with PERIODIC T)."// &
                          " This corresponds to a 1st order perturbation in the length gauge."// &
                          " Note that this is NOT compatible with a periodic calculation!"// &
                          " Uses the reference point defined in DFT%PRINT%MOMENTS ", &
                          usage="LEN_REP T", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PERIODIC", &
                          description="Apply a delta-kick that is compatible with periodic boundary conditions"// &
                          " for any value of DELTA_PULSE_SCALE. Uses perturbation theory for the preparation of"// &
                          " the initial wfn with the velocity operator as perturbation."// &
                          " If LEN_REP is .FALSE. this corresponds to a first order velocity gauge."// &
                          " Note that the pulse is only applied when INITIAL_WFN is set to SCF_WFN,"// &
                          " and not for restarts (RT_RESTART).", &
                          usage="PERIODIC", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DELTA_PULSE_DIRECTION", &
                          description="Direction of the applied electric field. The k vector is given as"// &
                          " 2*Pi*[i,j,k]*inv(h_mat), which for PERIODIC .FALSE. yields exp(ikr) periodic with"// &
                          " the unit cell, only if DELTA_PULSE_SCALE is set to unity. For an orthorhombic cell"// &
                          " [1,0,0] yields [2*Pi/L_x,0,0]. For small cells, this results in a very large kick.", &
                          usage="DELTA_PULSE_DIRECTION 1 1 1", n_var=3, default_i_vals=[1, 0, 0], &
                          type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DELTA_PULSE_SCALE", &
                          description="Scale the k vector, which for PERIODIC .FALSE. results in exp(ikr) no"// &
                          " longer being periodic with the unit cell. The norm of k is the strength of the"// &
                          " applied electric field in atomic units.", &
                          usage="DELTA_PULSE_SCALE 0.01 ", n_var=1, default_r_val=0.001_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HFX_BALANCE_IN_CORE", &
                          description="If HFX is used, this keyword forces a redistribution/recalculation"// &
                          " of the integrals, balanced with respect to the in core steps.", &
                          usage="HFX_BALANCE_IN_CORE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MCWEENY_MAX_ITER", &
                          description="Determines the maximum amount of McWeeny steps used after each converged"// &
                          " step in density propagation", &
                          usage="MCWEENY_MAX_ITER 2", default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="ACCURACY_REFINEMENT", &
         description="If using density propagation some parts should be calculated with a higher accuracy than the rest"// &
         " to reduce numerical noise. This factor determines by how much the filtering threshold is"// &
         " reduced for these calculations.", &
         usage="ACCURACY_REFINEMENT", default_i_val=100)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MCWEENY_EPS", &
                          description="Threshold after which McWeeny is terminated", &
                          usage="MCWEENY_EPS 0.00001", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (print_section)
      CALL section_create(print_section, __LOCATION__, name="PRINT", &
                          description="Section of possible print options for an RTP runs", &
                          repeats=.FALSE.)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "PROGRAM_RUN_INFO", &
                                       description="Controls the printing within real time propagation and Eherenfest dynamics", &
                                       print_level=low_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESTART", &
                                       description="Controls the dumping of the MO restart file during rtp. "// &
                                       "By default keeps a short history of three restarts. "// &
                                       "See also RESTART_HISTORY. In density propagation this controls the printing of "// &
                                       "density matrix.", &
                                       print_level=low_print_level, common_iter_levels=3, &
                                       each_iter_names=s2a("MD"), each_iter_values=[20], &
                                       add_last=add_last_numeric, filename="RESTART")
      CALL keyword_create(keyword, __LOCATION__, name="BACKUP_COPIES", &
                          description="Specifies the maximum number of backup copies.", &
                          usage="BACKUP_COPIES {int}", &
                          default_i_val=1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESTART_HISTORY", &
                                       description="Dumps unique MO restart files during the run keeping all of them. "// &
                                       "In density propagation it dumps the density matrix instead", &
                                       print_level=low_print_level, common_iter_levels=0, &
                                       each_iter_names=s2a("MD"), &
                                       each_iter_values=[500], &
                                       filename="RESTART")
      CALL keyword_create(keyword, __LOCATION__, name="BACKUP_COPIES", &
                          description="Specifies the maximum number of backup copies.", &
                          usage="BACKUP_COPIES {int}", &
                          default_i_val=1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "FIELD", &
                                       description="Print the time-dependent field applied during an EMD simulation in "// &
                                       "atomic unit.", &
                                       print_level=high_print_level, common_iter_levels=-1, &
                                       each_iter_names=s2a("MD"), &
                                       each_iter_values=[1], &
                                       filename="applied_field")
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      CALL create_projection_rtp_section(print_key)
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "CURRENT_INT", &
                                       description="Print the integral of the current density (only if the"// &
                                       " imaginary part of the density is NOT zero.", &
                                       print_level=high_print_level, common_iter_levels=1, &
                                       each_iter_names=s2a("MD"), &
                                       each_iter_values=[1], &
                                       filename="rtp_j_int")
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "CURRENT", &
                                       description="Print the current during an EMD simulation to cube files.", &
                                       print_level=high_print_level, common_iter_levels=0, &
                                       each_iter_names=s2a("MD"), &
                                       each_iter_values=[20], &
                                       filename="current")
      CALL keyword_create(keyword, __LOCATION__, name="BACKUP_COPIES", &
                          description="Specifies the maximum number of backup copies.", &
                          usage="BACKUP_COPIES {int}", &
                          default_i_val=1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=[2, 2, 2], type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      ! Marek : Add print option for ASCII density files - DEVELPMENT ONLY?
      CALL cp_print_key_section_create(print_key, __LOCATION__, "DENSITY_MATRIX", &
                                       description="Prints the density matrix at iterations in clear text to a file", &
                                       print_level=high_print_level, common_iter_levels=0, &
                                       each_iter_names=s2a("MD"), &
                                       each_iter_values=[1], &
                                       filename="rho")
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)
      ! Marek : Moments ASCII print
      CALL cp_print_key_section_create(print_key, __LOCATION__, "MOMENTS", &
                                       description="Prints the time-dependent electronic moments at "// &
                                       "iterations in clear text to a file.", &
                                       print_level=high_print_level, common_iter_levels=1, &
                                       each_iter_names=s2a("MD"), &
                                       each_iter_values=[1], &
                                       filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE", &
                          variants=s2a("REF"), &
                          description="Define the reference point for the calculation of the electrostatic moment.", &
                          usage="REFERENCE COM", &
                          enum_c_vals=s2a("COM", "COAC", "USER_DEFINED", "ZERO"), &
                          enum_desc=s2a("Use Center of Mass", &
                                        "Use Center of Atomic Charges", &
                                        "Use User Defined Point (Keyword:REFERENCE_POINT)", &
                                        "Use Origin of Coordinate System"), &
                          enum_i_vals=[use_mom_ref_com, &
                                       use_mom_ref_coac, &
                                       use_mom_ref_user, &
                                       use_mom_ref_zero], &
                          default_i_val=use_mom_ref_coac)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE_POINT", &
                          variants=s2a("REF_POINT"), &
                          description="Fixed reference point for the calculations of the electrostatic moment.", &
                          usage="REFERENCE_POINT x y z", &
                          repeats=.FALSE., &
                          n_var=3, default_r_vals=[0._dp, 0._dp, 0._dp], &
                          type_of_var=real_t, &
                          unit_str='angstrom')
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)
      ! Marek : Fourier transform of MOMENTS ASCII print
      CALL cp_print_key_section_create(print_key, __LOCATION__, "MOMENTS_FT", &
                                       description="Prints the calculated Fourier transform of "// &
                                       "time-dependent moments. For calculations with real time pulse (not delta kick) "// &
                                       "can be supplied with starting time.", &
                                       print_level=medium_print_level, common_iter_levels=0, &
                                       each_iter_names=s2a("MD"), &
                                       each_iter_values=[1], &
                                       filename="MOMENTS_FT")
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)
      ! Marek : Chosen element of (Fourier transformed) polarizability tensor (energy dependent) - text format
      CALL cp_print_key_section_create(print_key, __LOCATION__, "POLARIZABILITY", &
                                       description="Prints the chosen element of the energy dependent polarizability tensor "// &
                                       "to a specified file. The tensor is calculated as ratio of "// &
                                       "Fourier transform of the dipole "// &
                                       "moment trace and Fourier transform of the applied field "// &
                                       "(for delta kick, constant real field is applied.", &
                                       print_level=medium_print_level, common_iter_levels=0, &
                                       each_iter_names=s2a("MD"), &
                                       each_iter_values=[1], &
                                       filename="POLARIZABILITY")
      CALL keyword_create(keyword, __LOCATION__, "ELEMENT", &
                          description="Specifies the element of polarizability which is to be printed out "// &
                          "(indexing starts at 1). If not explicitly provided, RTBSE code tries to guess "// &
                          "the optimal values - for applied electric field (both delta pulse and RT field) "// &
                          "with only a single non-zero cartesian component, prints the 3 trivially available elements.", &
                          type_of_var=integer_t, default_i_vals=[1, 1], n_var=2, usage="ELEMENT 1 1", repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "E_CONSTITUENTS", &
                                       description="Print the energy constituents (relevant to RTP) which make up "// &
                                       "the Total Energy", &
                                       print_level=high_print_level, common_iter_levels=1, &
                                       each_iter_names=s2a("MD"), &
                                       each_iter_values=[1], &
                                       filename="rtp")
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)

      CALL section_add_subsection(section, print_section)
      CALL section_release(print_section)

      ! RTBSE subsection
      NULLIFY (subsection)
      CALL create_rtbse_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
      ! FT subsection
      CALL create_ft_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_rtp_section
! **************************************************************************************************
!> \brief Creates the subsection for specialized options of RTBSE code
!> \param section The created RTBSE section
!> \author Stepan Marek
! **************************************************************************************************
   SUBROUTINE create_rtbse_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, name="RTBSE", &
                          description="Controls options for the real-time Bethe-Salpeter (RTBSE) propagation. "// &
                          "Note that running RTBSE requires previous low-scaling "// &
                          "[GW](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW) calculation. Also note that "// &
                          "designating this section as RTBSE run but choosing run type ENERGY leads to potential "// &
                          "deallocation errors.", &
                          repeats=.FALSE., citations=[Marek2025])

      ! Marek : Controlling flow to RTBSE
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Which method is used for the time propagation of electronic structure. "// &
                          "By default, use the TDDFT method. Can also choose RT-BSE method, which propagates the lesser "// &
                          "Green's function instead of density matrix/molecular orbitals.", &
                          usage="&RTBSE TDDFT", &
                          default_i_val=rtp_method_tddft, &
                          lone_keyword_i_val=rtp_method_bse, &
                          enum_c_vals=s2a("TDDFT", "RTBSE"), &
                          enum_i_vals=[rtp_method_tddft, rtp_method_bse], &
                          enum_desc=s2a("Use TDDFT for density matrix/MO propagation.", &
                                        "Use RT-BSE for Green's function propagation"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Marek : Development option - run GWBSE starting from the KS Hamiltonian
      CALL keyword_create(keyword, __LOCATION__, name="RTBSE_HAMILTONIAN", &
                          description="Which Hamiltonian to use as the single-particle Hamiltonian"// &
                          " in the Green's propagator.", &
                          usage="RTBSE_HAMILTONIAN G0W0", &
                          default_i_val=rtp_bse_ham_g0w0, &
                          enum_c_vals=s2a("KS", "G0W0"), &
                          enum_i_vals=[rtp_bse_ham_ks, rtp_bse_ham_g0w0], &
                          enum_desc=s2a("Use Kohn-Sham Hamiltonian for Green's propagation.", &
                                        "Use G0W0 Hamiltonian for Green's function propagation"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_rtbse_section
! **************************************************************************************************
!> \brief Creates the subsection for Fourier transform options applicable to RTP output
!> \param ft_section The created FT section
!> \date 11.2025
!> \author Stepan Marek
! **************************************************************************************************
   SUBROUTINE create_ft_section(ft_section)
      TYPE(section_type), POINTER                        :: ft_section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(ft_section))

      ! Create the section itself
      CALL section_create(ft_section, __LOCATION__, &
                          name="FT", &
                          description="Define parameters for Fourier transforms used in RTP outputs.", &
                          repeats=.FALSE.)

      ! Start time keyword
      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, "START_TIME", &
                          description="The starting time from which damping is applied and from which on the trace is "// &
                          "considered for the Fourier transform (Fourier transform is used for the calculation of "// &
                          "MOMENTS_FT and POLARIZABILITY). Useful for real-time pulse - "// &
                          "one can specify the center of the pulse as the starting point.", &
                          type_of_var=real_t, &
                          unit_str="fs", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(ft_section, keyword)
      CALL keyword_release(keyword)

      ! Damping keyword
      CALL keyword_create(keyword, __LOCATION__, "DAMPING", &
                          description="Numerical Fourier transform (required for calculation of "// &
                          "MOMENTS_FT and POLARIZABILITY) can oscillate "// &
                          "when the final time trace values are far away from zero. "// &
                          "This keyword controls the exponential damping added to the Fourier transform "// &
                          "(Fourier transform is used for calculation of MOMENTS_FT and POLARIZABILITY). "// &
                          "For negative values (the default), calculates the damping at the run time so that the last point "// &
                          "in the time trace is reduced by factor e^(-4). When set manually, determines the time in which "// &
                          "the moments trace is reduced by factor of e^(-1), except when set to zero, in which case "// &
                          "the damping is not applied.", &
                          type_of_var=real_t, &
                          unit_str="fs", &
                          default_r_val=-1.0_dp/femtoseconds)
      CALL section_add_keyword(ft_section, keyword)
      CALL keyword_release(keyword)

      ! Create the Padé subsection
      NULLIFY (subsection)
      CALL section_create(subsection, __LOCATION__, name="PADE", &
                          description="Defines the parameters for the Padé interpolation of the "// &
                          "Fourier transforms used in the output of RTP. Only available with the GreenX library linked to CP2K.", &
                          repeats=.FALSE.)

      ! Explicit presence of the section turns on the Padé interpolation
      CALL keyword_create(keyword, __LOCATION__, "_SECTION_PARAMETERS_", &
                          description="Turns on the Padé interpolation", &
                          type_of_var=logical_t, &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      ! Minimum interpolated energy
      CALL keyword_create(keyword, __LOCATION__, "E_MIN", &
                          description="The minimum energy of the Padé interpolation output.", &
                          type_of_var=real_t, &
                          unit_str="eV", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      ! Maximum interpolated energy
      CALL keyword_create(keyword, __LOCATION__, "E_MAX", &
                          description="The maximum energy of the Padé interpolation output.", &
                          type_of_var=real_t, &
                          unit_str="eV", &
                          default_r_val=100.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      ! Energy resolution
      CALL keyword_create(keyword, __LOCATION__, "E_STEP", &
                          description="The energy resolution of the Padé interpolation output.", &
                          type_of_var=real_t, &
                          unit_str="eV", &
                          default_r_val=0.02_dp/evolt)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      ! Minimum fitting energy
      CALL keyword_create(keyword, __LOCATION__, "FIT_E_MIN", &
                          description="The lower boundary in energy for the points "// &
                          "used in the fitting of Padé parameters. If negative, uses "// &
                          "value of E_MIN (default).", &
                          type_of_var=real_t, &
                          unit_str="eV", &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      ! Maximum fitting energy
      CALL keyword_create(keyword, __LOCATION__, "FIT_E_MAX", &
                          description="The upper boundary in energy for the points "// &
                          "used in the fitting of Padé parameters. If negative, uses "// &
                          "the value of E_MAX (default).", &
                          type_of_var=real_t, &
                          unit_str="eV", &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      ! Add the Padé subsection
      CALL section_add_subsection(ft_section, subsection)
      CALL section_release(subsection)
   END SUBROUTINE create_ft_section

! **************************************************************************************************
!> \brief      Create CP2K input section for the SCCS model
!> \param section ...
!> \par History:
!>      - Creation (10.10.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE create_sccs_section(section)

      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, &
                          name="SCCS", &
                          description="Define the parameters for self-consistent continuum solvation (SCCS) model", &
                          citations=[Fattebert2002, Andreussi2012, Yin2017], &
                          n_keywords=8, &
                          n_subsections=2, &
                          repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="_SECTION_PARAMETERS_", &
                          description="Controls the activation of the SCCS section", &
                          usage="&SCCS ON", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="ALPHA", &
                          description="Solvent specific tunable parameter for the calculation of "// &
                          "the repulsion term $G^\text{rep} = \alpha S$ "// &
                          "where $S$ is the (quantum) surface of the cavity", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.0_dp, &
                          unit_str="mN/m")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="BETA", &
                          description="Solvent specific tunable parameter for the calculation of "// &
                          "the dispersion term $G^\text{dis} = \beta V$ "// &
                          "where $V$ is the (quantum) volume of the cavity", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.0_dp, &
                          unit_str="GPa")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="DELTA_RHO", &
                          description="Numerical increment for the calculation of the (quantum) "// &
                          "surface of the solute cavity", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=2.0E-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="DERIVATIVE_METHOD", &
                          description="Method for the calculation of the numerical derivatives on the real-space grids", &
                          usage="DERIVATIVE_METHOD cd5", &
                          repeats=.FALSE., &
                          n_var=1, &
                          default_i_val=sccs_derivative_fft, &
                          enum_c_vals=s2a("FFT", "CD3", "CD5", "CD7"), &
                          enum_i_vals=[sccs_derivative_fft, &
                                       sccs_derivative_cd3, &
                                       sccs_derivative_cd5, &
                                       sccs_derivative_cd7], &
                          enum_desc=s2a("Fast Fourier transformation", &
                                        "3-point stencil central differences", &
                                        "5-point stencil central differences", &
                                        "7-point stencil central differences"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="RELATIVE_PERMITTIVITY", &
                          variants=s2a("DIELECTRIC_CONSTANT", "EPSILON_RELATIVE", "EPSILON_SOLVENT"), &
                          description="Relative permittivity (dielectric constant) of the solvent (medium)", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=80.0_dp, &
                          usage="RELATIVE_PERMITTIVITY 78.36")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="EPS_SCCS", &
                          variants=s2a("EPS_ITER", "TAU_POL"), &
                          description="Tolerance for the convergence of the polarisation density, "// &
                          "i.e. requested accuracy for the SCCS iteration cycle", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=1.0E-6_dp, &
                          usage="EPS_ITER 1.0E-7")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="EPS_SCF", &
                          description="The SCCS iteration cycle is activated only if the SCF iteration cycle "// &
                          "is converged to this threshold value", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.5_dp, &
                          usage="EPS_SCF 1.0E-2")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="GAMMA", &
                          variants=s2a("SURFACE_TENSION"), &
                          description="Surface tension of the solvent used for the calculation of "// &
                          "the cavitation term $G^\text{cav} = \gamma S$ "// &
                          "where $S$ is the (quantum) surface of the cavity", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.0_dp, &
                          unit_str="mN/m")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="MAX_ITER", &
                          description="Maximum number of SCCS iteration steps performed to converge "// &
                          "within the given tolerance", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=integer_t, &
                          default_i_val=100, &
                          usage="MAX_ITER 50")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="METHOD", &
                          description="Method used for the smoothing of the dielectric function", &
                          usage="METHOD Fattebert-Gygi", &
                          default_i_val=sccs_andreussi, &
                          enum_c_vals=s2a("ANDREUSSI", "FATTEBERT-GYGI"), &
                          enum_i_vals=[sccs_andreussi, sccs_fattebert_gygi], &
                          enum_desc=s2a("Smoothing function proposed by Andreussi et al.", &
                                        "Smoothing function proposed by Fattebert and Gygi"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="MIXING", &
                          variants=["ETA"], &
                          description="Mixing parameter (Hartree damping) employed during the iteration procedure", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.6_dp, &
                          usage="MIXING 0.2")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)

      CALL section_create(subsection, __LOCATION__, &
                          name="ANDREUSSI", &
                          description="Define the parameters of the dielectric smoothing function proposed by "// &
                          "Andreussi et al.", &
                          citations=[Andreussi2012], &
                          n_keywords=2, &
                          n_subsections=0, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="RHO_MAX", &
                          description="Maximum density value used for the smoothing of the dielectric function", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.0035_dp, &
                          usage="RHO_MAX 0.01")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="RHO_MIN", &
                          description="Minimum density value used for the smoothing of the dielectric function", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.0001_dp, &
                          usage="RHO_MIN 0.0003")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL section_create(subsection, __LOCATION__, &
                          name="FATTEBERT-GYGI", &
                          description="Define the parameters of the dielectric smoothing function proposed by "// &
                          "Fattebert and Gygi", &
                          citations=[Fattebert2002], &
                          n_keywords=2, &
                          n_subsections=0, &
                          repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="BETA", &
                          description="Parameter &beta; changes the width of the interface solute-solvent", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=1.7_dp, &
                          usage="BETA 1.3")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="RHO_ZERO", &
                          variants=["RHO0"], &
                          description="Parameter $\rho_0$ defines the critical density in the middle "// &
                          "of the interface solute-solvent", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.0006_dp, &
                          usage="RHO_ZERO 0.0004")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_sccs_section

END MODULE input_cp2k_dft
