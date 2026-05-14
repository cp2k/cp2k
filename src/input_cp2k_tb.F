!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the dft section of the input
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_tb
   USE bibliography,                    ONLY: Elstner1998,&
                                              Grimme2017,&
                                              Hu2007,&
                                              Katbashev2025,&
                                              Porezag1995,&
                                              Seifert1996,&
                                              Zhechkov2005
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE eeq_input,                       ONLY: create_eeq_control_section
   USE input_constants,                 ONLY: &
        dispersion_d2, dispersion_d3, dispersion_d3bj, dispersion_uff, gfn1xtb, gfn2xtb, &
        gfn_tblite, ipea1xtb, slater, tblite_cli_born_kernel_auto, tblite_cli_born_kernel_p16, &
        tblite_cli_born_kernel_still, tblite_cli_solution_state_bar1mol, &
        tblite_cli_solution_state_gsolv, tblite_cli_solution_state_reference, &
        tblite_cli_solvation_alpb, tblite_cli_solvation_cpcm, tblite_cli_solvation_gb, &
        tblite_cli_solvation_gbe, tblite_cli_solvation_gbsa, tblite_guess_ceh, tblite_guess_eeq, &
        tblite_guess_sad, tblite_mixer_damping_default, tblite_mixer_iterations_default, &
        tblite_mixer_max_weight_default, tblite_mixer_memory_inherit, &
        tblite_mixer_min_weight_default, tblite_mixer_omega0_default, &
        tblite_mixer_weight_factor_default, tblite_scc_mixer_auto, tblite_scc_mixer_cp2k, &
        tblite_scc_mixer_none, tblite_scc_mixer_tblite, tblite_solver_gvd, tblite_solver_gvr
   USE input_cp2k_mm,                   ONLY: create_GENPOT_section
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_tb'

   PUBLIC :: create_dftb_control_section, create_xtb_control_section

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_dftb_control_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="DFTB", &
                          description="Parameters needed to set up the DFTB methods", &
                          n_keywords=1, n_subsections=1, repeats=.FALSE., &
                          citations=[Porezag1995, Seifert1996, Elstner1998, Zhechkov2005])

      NULLIFY (subsection)
      CALL create_dftb_parameter_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_tblite_mixer_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="self_consistent", &
                          description="Use self-consistent method", &
                          citations=[Elstner1998], &
                          usage="SELF_CONSISTENT", default_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="orthogonal_basis", &
                          description="Assume orthogonal basis set", &
                          usage="ORTHOGONAL_BASIS", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="do_ewald", &
                          description="Deprecated debug override for the DFTB Ewald treatment. "// &
                          "By default the setting is inferred from SUBSYS/CELL/PERIODIC. "// &
                          "If set to TRUE, the DFT/POISSON section must use a periodic Poisson solver.", &
                          usage="DO_EWALD", default_l_val=.FALSE., lone_keyword_l_val=.TRUE., &
                          deprecation_notice="Use SUBSYS/CELL/PERIODIC to select periodicity. "// &
                          "DO_EWALD is retained only as a debug override for the DFTB Ewald path.")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="dispersion", &
                          description="Use dispersion correction", &
                          citations=[Zhechkov2005], lone_keyword_l_val=.TRUE., &
                          usage="DISPERSION", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DIAGONAL_DFTB3", &
                          description="Use a diagonal version of the 3rd order energy correction (DFTB3) ", &
                          lone_keyword_l_val=.TRUE., &
                          usage="DIAGONAL_DFTB3", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HB_SR_GAMMA", &
                          description="Uses a modified version for the GAMMA within the SCC-DFTB scheme, "// &
                          "specifically tuned for hydrogen bonds.", &
                          citations=[Hu2007], lone_keyword_l_val=.TRUE., &
                          usage="HB_SR_GAMMA", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="eps_disp", &
                          description="Define accuracy of dispersion interaction", &
                          usage="EPS_DISP", default_r_val=0.0001_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL add_tb_scc_mixer_keywords(section, "CP2K-internal DFTB", &
                                     "AUTO keeps the existing CP2K density/SCF mixing path for DFTB.")

   END SUBROUTINE create_dftb_control_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_xtb_control_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="xTB", &
                          description="Parameters needed to set up the xTB methods", &
                          n_keywords=1, n_subsections=1, repeats=.FALSE., &
                          citations=[GRIMME2017])

      NULLIFY (subsection)
      CALL create_xtb_parameter_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_xtb_nonbonded_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_eeq_control_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_tblite_mixer_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_xtb_tblite_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="GFN_TYPE", &
                          description="Selects the xTB backend. GFN_TYPE 1 uses the CP2K-internal GFN1-xTB "// &
                          "implementation. GFN_TYPE 0 uses the CP2K-internal GFN0-xTB implementation, which "// &
                          "has no SCC variables to mix, so XTB/SCC_MIXER AUTO is treated as NONE. "// &
                          "GFN_TYPE TBLITE uses the tblite backend and requires an XTB/TBLITE "// &
                          "section, where the actual tblite method is selected.", &
                          usage="GFN_TYPE (0|1|TBLITE)", &
                          enum_c_vals=s2a("0", "1", "TBLITE"), &
                          enum_i_vals=[0, gfn1xtb, gfn_tblite], &
                          enum_desc=s2a("Use the CP2K-internal GFN0-xTB implementation.", &
                                        "Use the CP2K-internal GFN1-xTB implementation.", &
                                        "Use the CP2K/tblite backend; requires XTB/TBLITE."), &
                          default_i_val=gfn1xtb)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL add_tb_scc_mixer_keywords(section, "xTB", &
                                     "AUTO is treated as NONE for CP2K-internal GFN0-xTB, keeps the CP2K "// &
                                     "charge mixer for CP2K-internal GFN1-xTB, and uses "// &
                                     "tblite's native SCC mixer for XTB/GFN_TYPE TBLITE.")

      CALL keyword_create(keyword, __LOCATION__, name="DO_EWALD", &
                          description="Deprecated debug override for the xTB Ewald treatment. "// &
                          "By default the setting is inferred from SUBSYS/CELL/PERIODIC. "// &
                          "If set to TRUE for CP2K-internal xTB, the DFT/POISSON section must use a periodic "// &
                          "Poisson solver. For CP2K/tblite, the periodic long-range treatment is performed "// &
                          "inside tblite.", &
                          usage="DO_EWALD", default_l_val=.FALSE., lone_keyword_l_val=.TRUE., &
                          deprecation_notice="Use SUBSYS/CELL/PERIODIC to select periodicity. "// &
                          "DO_EWALD is retained only as a debug override for the CP2K-internal xTB Ewald path.")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STO_NG", &
                          description="Provides the order of the Slater orbital expansion in GTOs.", &
                          usage="STO_NG 3", default_i_val=6)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HYDROGEN_STO_NG", &
                          description="Number of GTOs for Hydrogen basis expansion.", &
                          usage="HYDROGEN_STO_NG 3", default_i_val=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="USE_HALOGEN_CORRECTION", &
                          description="Use XB interaction term", &
                          usage="USE_HALOGEN_CORRECTION T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DO_NONBONDED", &
                          description="Controls the computation of real-space "// &
                          "(short-range) nonbonded interactions as correction to xTB.", &
                          usage="DO_NONBONDED T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VDW_POTENTIAL", &
                          description="vdW potential to be used: NONE, DFTD3, DFTD4. "// &
                          "Defaults: DFTD3(gfn1), DFTD4(gfn0, gfn2).", &
                          usage="VDW_POTENTIAL type", default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COULOMB_INTERACTION", &
                          description="Use Coulomb interaction terms (electrostatics + TB3); for debug only", &
                          usage="COULOMB_INTERACTION T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COULOMB_LR", &
                          description="Use Coulomb LR (1/r) interaction terms; for debug only", &
                          usage="COULOMB_LR T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TB3_INTERACTION", &
                          description="Use TB3 interaction terms; for debug only", &
                          usage="TB3_INTERACTION T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHECK_ATOMIC_CHARGES", &
                          description="Stop calculation if atomic charges are outside chemical range.", &
                          usage="CHECK_ATOMIC_CHARGES T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VARIATIONAL_DIPOLE", &
                          description="gfn0-xTB use dipole definition from energy derivative.", &
                          usage="VARIATIONAL_DIPOLE T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_PAIRPOTENTIAL", &
                          description="Accuracy for the repulsive pair potential.", &
                          usage="EPS_PAIRPOTENTIAL 1.0E-8", default_r_val=1.0e-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EN_SHIFT_TYPE", &
                          description="Shift function for electronegativity in EEQ method. "// &
                          "[Select/Molecule/Crystal] Default Select from periodicity.", &
                          usage="EN_SHIFT_TYPE [Select/Molecule/Crystal]", &
                          n_var=1, type_of_var=char_t, default_c_val="Molecule")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_xtb_control_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_dftb_parameter_section(section)

      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, name="PARAMETER", &
                          description="Information on where to find DFTB parameters", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SK_FILE", &
                          description="Define parameter file for atom pair", &
                          usage="SK_FILE a1 a2 filename", &
                          n_var=3, type_of_var=char_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARAM_FILE_PATH", &
                          description="Specify the directory with the DFTB parameter files. "// &
                          "Used in combination with the filenames specified in the file "// &
                          "given in PARAM_FILE_NAME.", usage="PARAM_FILE_PATH pathname", &
                          n_var=1, type_of_var=char_t, default_c_val="./")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARAM_FILE_NAME", &
                          description="Specify file that contains the names of "// &
                          "Slater-Koster tables: A plain text file, each line has the "// &
                          'format "ATOM1 ATOM2 filename.spl".', &
                          usage="PARAM_FILE_NAME filename", &
                          n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISPERSION_TYPE", &
                          description="Use dispersion correction of the specified type."// &
                          " Dispersion correction has to be switched on in the DFTB section.", &
                          usage="DISPERSION_TYPE (UFF|D3|D3(BJ)|D2)", &
                          enum_c_vals=s2a("UFF", "D3", "D3(BJ)", "D2"), &
                          enum_i_vals=[dispersion_uff, dispersion_d3, dispersion_d3bj, dispersion_d2], &
                          enum_desc=s2a("Uses the UFF force field for a pair potential dispersion correction.", &
                                        "Uses the Grimme D3 method (simplified) for a pair potential dispersion correction.", &
                                        "Uses the Grimme D3 method (simplified) with Becke-Johnson attenuation.", &
                                        "Uses the Grimme D2 method for pair potential dispersion correction."), &
                          default_i_val=dispersion_uff)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="UFF_FORCE_FIELD", &
                          description="Name of file with UFF parameters that will be used "// &
                          "for the dispersion correction. Needs to be specified when "// &
                          "DISPERSION==.TRUE., otherwise cp2k crashes with a Segmentation "// &
                          "Fault.", usage="UFF_FORCE_FIELD filename", &
                          n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISPERSION_PARAMETER_FILE", &
                          description="Specify file that contains the atomic dispersion "// &
                          "parameters for the D3 method", &
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
                          usage="D3_SCALING 1.0 1.0 1.0", n_var=3, default_r_vals=[0.0_dp, 0.0_dp, 0.0_dp])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="D3BJ_SCALING", &
                          description="Scaling parameters (s6,a1,s8,a2) for the D3(BJ) dispersion method,", &
                          usage="D3BJ_SCALING 1.0 1.0 1.0 1.0", n_var=4, &
                          default_r_vals=[0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="D2_SCALING", &
                          description="Scaling parameter for the D2 dispersion method,", &
                          usage="D2_SCALING 1.0", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="D2_EXP_PRE", &
                          description="Exp prefactor for damping for the D2 dispersion method,", &
                          usage="D2_EXP_PRE 2.0", default_r_val=2.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HB_SR_PARAM", &
                          description="Uses a modified version for the GAMMA within the SCC-DFTB scheme, "// &
                          "specifically tuned for hydrogen bonds. Specify the exponent used in the exponential.", &
                          usage="HB_SR_PARAM {real}", default_r_val=4.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_dftb_parameter_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_xtb_parameter_section(section)

      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, name="PARAMETER", &
                          description="Information on and where to find xTB parameters", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="PARAM_FILE_PATH", &
                          description="Specify the directory with the xTB parameter file. ", &
                          usage="PARAM_FILE_PATH pathname", &
                          n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARAM_FILE_NAME", &
                          description="Specify file that contains all xTB default parameters. ", &
                          usage="PARAM_FILE_NAME filename", &
                          n_var=1, type_of_var=char_t, default_c_val="xTB_parameters")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DISPERSION_PARAMETER_FILE", &
                          description="Specify file that contains the atomic dispersion "// &
                          "parameters for the D3 method", &
                          usage="DISPERSION_PARAMETER_FILE filename", &
                          n_var=1, type_of_var=char_t, default_c_val="dftd3.dat")
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

      CALL keyword_create(keyword, __LOCATION__, name="D3BJ_SCALING", &
                          description="Scaling parameters (s6,s8) for the D3 dispersion method.", &
                          usage="D3BJ_SCALING 1.0 2.4", n_var=2, default_r_vals=[1.0_dp, 2.4_dp])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="D3BJ_PARAM", &
                          description="Becke-Johnson parameters (a1, a2 for the D3 dispersion method.", &
                          usage="D3BJ_PARAM 0.63 5.0", n_var=2, default_r_vals=[0.63_dp, 5.0_dp])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HUCKEL_CONSTANTS", &
                          description="Huckel parameters (s, p, d, sp, 2sH).", &
                          usage="HUCKEL_CONSTANTS 1.85 2.25 2.00 2.08 2.85", n_var=5, &
                          default_r_vals=[1.85_dp, 2.25_dp, 2.00_dp, 2.08_dp, 2.85_dp])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COULOMB_CONSTANTS", &
                          description="Scaling parameters for Coulomb interactions (electrons, nuclei).", &
                          usage="COULOMB_CONSTANTS 2.00 1.50", n_var=2, &
                          default_r_vals=[2.00_dp, 1.50_dp])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CN_CONSTANTS", &
                          description="Scaling parameters for Coordination number correction term.", &
                          usage="CN_CONSTANTS 0.006 -0.003 -0.005", n_var=3, &
                          default_r_vals=[0.006_dp, -0.003_dp, -0.005_dp])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EN_CONSTANTS", &
                          description="Scaling parameters for electronegativity correction term.", &
                          usage="EN_CONSTANTS -0.007 0.000 0.000", n_var=3, &
                          default_r_vals=[-0.007_dp, 0.000_dp, 0.000_dp])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BEN_CONSTANT", &
                          description="Scaling parameter for electronegativity correction term.", &
                          usage="BEN_CONSTANT 4.0", n_var=1, &
                          default_r_val=4.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ENSCALE", &
                          description="Scaling parameter repulsive energy (dEN in exponential).", &
                          usage="ENSCALE 0.01", n_var=1, &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="HALOGEN_BINDING", &
                          description="Scaling parameters for electronegativity correction term.", &
                          usage="HALOGEN_BINDING 1.30 0.44", n_var=2, default_r_vals=[1.30_dp, 0.44_dp])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="KAB_PARAM", &
                          description="Specifies the specific Kab value for types A and B.", &
                          usage="KAB_PARAM kind1 kind2 value ", repeats=.TRUE., &
                          n_var=-1, type_of_var=char_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="XB_RADIUS", &
                          description="Specifies the radius [Bohr] of the XB pair interaction in xTB.", &
                          usage="XB_RADIUS 20.0 ", repeats=.FALSE., &
                          n_var=1, default_r_val=20.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COULOMB_SR_CUT", &
                          description="Maximum range of short range part of Coulomb interaction.", &
                          usage="COULOMB_SR_CUT 20.0 ", repeats=.FALSE., &
                          n_var=1, default_r_val=20.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COULOMB_SR_EPS", &
                          description="Cutoff for short range part of Coulomb interaction.", &
                          usage="COULOMB_SR_EPS 1.E-3 ", repeats=.FALSE., &
                          n_var=1, default_r_val=1.0E-03_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SRB_PARAMETER", &
                          description="SRB parameters (ksrb, esrb, gscal, c1, c2, shift).", &
                          usage="SRB_PARAMETER -0.0129 3.48 0.51 -1.71 2.11 0.0537", n_var=6, &
                          default_r_vals=[-0.0129_dp, 3.4847_dp, 0.5097_dp, &
                                          -1.70549806_dp, 2.10878369_dp, 0.0537_dp])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_xtb_parameter_section
! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_xtb_nonbonded_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="NONBONDED", &
                          description="This section specifies the input parameters for NON-BONDED interactions.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      NULLIFY (subsection)

      CALL create_GENPOT_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="DX", &
                          description="Parameter used for computing the derivative with the Ridders' method.", &
                          usage="DX <REAL>", default_r_val=0.1_dp, unit_str="bohr")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ERROR_LIMIT", &
                          description="Checks that the error in computing the derivative is not larger than "// &
                          "the value set; in case error is larger a warning message is printed.", &
                          usage="ERROR_LIMIT <REAL>", default_r_val=1.0E-12_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_xtb_nonbonded_section
! **************************************************************************************************
!> \brief Creates the &TBLITE section
!> \param section the section to create
!> \author JVP
! **************************************************************************************************
   SUBROUTINE create_xtb_tblite_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="TBLITE", &
                          description="Section used to specify options for an xTB computation using tblite. "// &
                          "This section is required when XTB/GFN_TYPE is set to TBLITE. "// &
                          "Periodicity follows SUBSYS/CELL/PERIODIC; the deprecated XTB/DO_EWALD keyword "// &
                          "is read consistently with CP2K-internal xTB for diagnostics/debugging. "// &
                          "Initial densities and restarts are controlled by the regular CP2K SCF settings, "// &
                          "for example DFT/SCF/SCF_GUESS, not by a separate tblite guess in this section. "// &
                          "Spin-polarized CP2K/tblite calculations are selected with the regular DFT/UKS keyword "// &
                          "or its LSD/SPIN_POLARIZED aliases; DFT/MULTIPLICITY sets the number of unpaired "// &
                          "electrons. UKS uses tblite's spin-polarized Hamiltonian; AUTO selects the native tblite "// &
                          "SCC mixer for XTB/GFN_TYPE TBLITE, while SCC_MIXER CP2K mixes both total-charge and "// &
                          "magnetization SCC variables through CP2K's TB charge-mixing infrastructure. Native "// &
                          "tblite SCC convergence follows ACCURACY; SCC_MIXER CP2K follows the regular CP2K "// &
                          "SCF convergence settings such as DFT/SCF/EPS_SCF. With QS/LS_SCF, CP2K disables the "// &
                          "separate SCC-variable mixer and lets LS_SCF optimize the density matrix directly. "// &
                          "The same direct SCC-variable update is available for Gamma-point OT calculations, "// &
                          "where smearing and ADDED_MOS are not supported.", &
                          n_keywords=1, n_subsections=1, repeats=.FALSE., citations=[Katbashev2025])

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Marks the tblite options section as present. The actual switch to the tblite "// &
                          "backend is XTB/GFN_TYPE TBLITE; this section then collects the numerical model, SCC "// &
                          "mixing, and optional native-reference settings used by that backend.", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="METHOD", &
                          description="Selects the xTB Hamiltonian model used by tblite. Without PARAM this chooses "// &
                          "one of tblite's built-in parametrizations, including the basis, reference occupations, "// &
                          "repulsion, dispersion, electrostatic, and coordination-number dependent terms. When PARAM "// &
                          "is given, METHOD remains the declared model family that the custom parametrization is "// &
                          "expected to implement and that CP2K uses for method-dependent setup and diagnostics.", &
                          usage="METHOD (GFN1|GFN2|IPEA1)", &
                          enum_c_vals=s2a("GFN1", "GFN2", "IPEA1"), &
                          enum_i_vals=[gfn1xtb, gfn2xtb, ipea1xtb], &
                          enum_desc=s2a("Uses the GFN1-XTB method by Grimme.", &
                                        "Uses the GFN2-XTB method by Grimme.", &
                                        "Uses the IPEA1 method by Grimme."), &
                          default_i_val=gfn2xtb)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARAM", &
                          description="Custom tblite parametrization file. The file is read by tblite and used to "// &
                          "construct the xTB calculator instead of the built-in parameter table for METHOD. It may "// &
                          "therefore change model constants such as Hamiltonian, repulsion, dispersion, electrostatic, "// &
                          "and coordination-number response parameters. It does not change the CP2K charge, spin, "// &
                          "cell, k-point, or SCF setup. The parametrization must be compatible with the selected "// &
                          "METHOD; CP2K does not reinterpret or validate the physical content of the file.", &
                          usage="PARAM filename", n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ACCURACY", &
                          description="Numerical accuracy factor for the CP2K/tblite backend. This scales tblite's "// &
                          "internal accuracy-dependent thresholds. With the native tblite SCC mixer, this is the "// &
                          "user-visible convergence control and corresponds to native tblite's --acc setting. With "// &
                          "SCC_MIXER CP2K, SCF convergence instead follows DFT/SCF/EPS_SCF; ACCURACY still controls "// &
                          "the tblite numerical setup, such as the distance cutoff used to build the xTB basis and "// &
                          "neighbor data. "// &
                          "The default matches native tblite. Smaller values make the calculation tighter and usually "// &
                          "more expensive; larger values loosen the numerical setup and can affect energies, forces, "// &
                          "and stress. Native REFERENCE_CLI diagnostics inherit this value so the production and "// &
                          "external reference calculations use the same accuracy setting.", &
                          usage="ACCURACY 1.0", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE_CLI", &
                          description="Enable an optional diagnostic comparison against an external native tblite "// &
                          "command line calculation. CP2K writes the current geometry, charge, spin, method, optional "// &
                          "PARAM file, electronic temperature, and requested reference settings to temporary files, "// &
                          "runs the external program, and compares the resulting energy, gradient, and optionally "// &
                          "virial with the CP2K/tblite values. This check is meant for debugging backend consistency, "// &
                          "not for production sampling. The native tblite CLI exposes only an electronic temperature "// &
                          "for the occupation broadening; if CP2K uses another smearing method, the reference run "// &
                          "prints a warning and uses Fermi-Dirac smearing with the CP2K electronic temperature. "// &
                          "If true, the XTB/TBLITE/REFERENCE_CLI section must be present. "// &
                          "The check is also enabled whenever that section is present. It is skipped for CP2K "// &
                          "multi-k-point calculations, because the native tblite CLI does not reproduce CP2K KPOINTS "// &
                          "sampling.", &
                          usage="REFERENCE_CLI", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL create_xtb_reference_cli_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_xtb_tblite_section

! **************************************************************************************************
!> \brief Creates the common &TBLITE_MIXER section for TB SCC-variable mixing.
!> \param section section to create
! **************************************************************************************************
   SUBROUTINE create_tblite_mixer_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="TBLITE_MIXER", &
                          description="Controls tblite's modified Broyden SCC-variable mixer when SCC_MIXER is set "// &
                          "to TBLITE, either explicitly or through AUTO where AUTO selects tblite. The mixer "// &
                          "variables are the self-consistent "// &
                          "xTB/DFTB populations and, for methods that require them, atomic multipoles. ITERATIONS "// &
                          "limits the number of tblite SCC-mixer updates, while MEMORY controls how many previous "// &
                          "Broyden update vectors are retained. These are independent quantities when MEMORY is set "// &
                          "explicitly; otherwise MEMORY follows ITERATIONS, matching the tblite library default. "// &
                          "DAMPING controls the size of the retained population/multipole update. OMEGA0, "// &
                          "MIN_WEIGHT, MAX_WEIGHT, and WEIGHT_FACTOR expose tblite's modified-Broyden weighting "// &
                          "constants. ITERATIONS is the tblite-side SCC-mixer update limit and defaults to the "// &
                          "native tblite value of 250. When SCC_MIXER is explicitly set to TBLITE, CP2K also "// &
                          "sets the regular DFT/SCF/MAX_SCF limit to ITERATIONS so that the tblite SCC cycle is "// &
                          "not truncated by CP2K's default SCF limit. With SCC_MIXER AUTO or CP2K, DFT/SCF/MAX_SCF "// &
                          "keeps its normal CP2K meaning. This section configures the tblite SCC mixer, not the "// &
                          "regular CP2K charge mixer. With SCC_MIXER CP2K, the CP2K-side SCC-variable mixer is "// &
                          "controlled by DFT/SCF/MIXING. With QS/LS_SCF, this section is ignored and CP2K resets "// &
                          "SCC_MIXER to NONE internally, because LS_SCF performs the density-matrix optimization "// &
                          "instead of a separate SCC-variable mixing step. Native tblite "// &
                          "REFERENCE_CLI diagnostics inherit this section. ITERATIONS and SOLVER are representable "// &
                          "by tblite run directly; CP2K reports inherited MEMORY, DAMPING, or Broyden-weight settings "// &
                          "that native tblite CLI cannot reproduce.", &
                          n_keywords=1, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="ITERATIONS", &
                          description="Maximum number of tblite SCC-mixer updates. This is the tblite-side SCC "// &
                          "iteration limit, matching native tblite's --iterations setting for reference CLI "// &
                          "diagnostics. For explicit SCC_MIXER TBLITE runs, CP2K synchronizes DFT/SCF/MAX_SCF to "// &
                          "this value. For SCC_MIXER AUTO or CP2K, DFT/SCF/MAX_SCF remains the normal CP2K "// &
                          "SCF-loop iteration limit. ITERATIONS is independent of MEMORY, which controls the "// &
                          "Broyden history length.", &
                          usage="ITERATIONS 250", default_i_val=tblite_mixer_iterations_default)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MEMORY", &
                          description="Number of previous modified-Broyden update vectors kept by the tblite SCC "// &
                          "mixer. This is a history length, not an iteration limit. Larger values retain more "// &
                          "history and use more memory; smaller values reduce the history and can be more robust for "// &
                          "difficult SCC cycles. A value of zero means library default: inherit ITERATIONS, as native "// &
                          "tblite does when it constructs the Broyden mixer from calc%max_iter.", &
                          usage="MEMORY 250", default_i_val=tblite_mixer_memory_inherit)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SOLVER", &
                          description="Electronic eigensolver selected for native tblite SCC calculations where "// &
                          "tblite owns the electronic solve, and forwarded to native REFERENCE_CLI diagnostics as "// &
                          "tblite run --solver. GVD is tblite's default LAPACK divide-and-conquer generalized "// &
                          "symmetric eigensolver. GVR selects LAPACK's relatively robust generalized symmetric "// &
                          "eigensolver. Production CP2K/tblite calculations still use CP2K's electronic "// &
                          "diagonalization machinery.", &
                          usage="SOLVER (GVD|GVR)", &
                          enum_c_vals=s2a("GVD", "GVR"), &
                          enum_i_vals=[tblite_solver_gvd, tblite_solver_gvr], &
                          enum_desc=s2a("Use LAPACK's divide-and-conquer generalized symmetric eigensolver.", &
                                        "Use LAPACK's relatively robust generalized symmetric eigensolver."), &
                          default_i_val=tblite_solver_gvd)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DAMPING", &
                          description="Damping parameter for the tblite modified Broyden SCC mixer. It controls the "// &
                          "size of the population/multipole update retained by the native tblite mixer between SCC "// &
                          "iterations. Smaller values damp oscillatory SCC cycles more strongly but may slow "// &
                          "convergence; larger values can accelerate well-behaved cases but may destabilize difficult "// &
                          "ones. For CP2K-internal GFN1-xTB UKS k-point calculations with explicit SCC_MIXER TBLITE, "// &
                          "CP2K uses a conservative damping of 0.25 unless DAMPING is set explicitly.", &
                          usage="DAMPING 0.4", default_r_val=tblite_mixer_damping_default)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OMEGA0", &
                          description="Regularization weight added to the diagonal of tblite's modified-Broyden "// &
                          "linear system. The default matches the hard-coded value in native tblite and damps nearly "// &
                          "singular history updates without changing the first simple-damping step.", &
                          usage="OMEGA0 0.01", default_r_val=tblite_mixer_omega0_default)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MIN_WEIGHT", &
                          description="Lower bound for the dynamic Broyden history weight derived from the current "// &
                          "SCC residual. The default matches native tblite and prevents tiny weights from removing "// &
                          "history vectors from the update.", &
                          usage="MIN_WEIGHT 1.0", default_r_val=tblite_mixer_min_weight_default)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_WEIGHT", &
                          description="Upper bound for the dynamic Broyden history weight derived from the current "// &
                          "SCC residual. The default matches native tblite and caps the influence of very small "// &
                          "residual norms.", &
                          usage="MAX_WEIGHT 100000.0", default_r_val=tblite_mixer_max_weight_default)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WEIGHT_FACTOR", &
                          description="Residual-to-weight scaling factor used by tblite's modified Broyden mixer. "// &
                          "For a residual norm above WEIGHT_FACTOR/MAX_WEIGHT, the dynamic weight is approximately "// &
                          "WEIGHT_FACTOR divided by that residual norm; otherwise MAX_WEIGHT is used before applying "// &
                          "MIN_WEIGHT. The default matches native tblite.", &
                          usage="WEIGHT_FACTOR 0.01", default_r_val=tblite_mixer_weight_factor_default)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_tblite_mixer_section

! **************************************************************************************************
!> \brief Adds common SCC-mixer keywords for TB methods.
!> \param section section receiving the keywords
!> \param TARGET method/path description for the manual text
!> \param auto_desc AUTO behavior for this section
! **************************************************************************************************
   SUBROUTINE add_tb_scc_mixer_keywords(section, TARGET, auto_desc)
      TYPE(section_type), POINTER                        :: section
      CHARACTER(LEN=*), INTENT(IN)                       :: target, auto_desc

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="SCC_MIXER", &
                          description="Selects the SCC variable mixer for "//TARGET//" calculations. "// &
                          auto_desc//" The SCC variables are the model populations and, where required by the "// &
                          "selected method, atomic multipoles that feed the self-consistent xTB potential. TBLITE "// &
                          "uses tblite's modified Broyden mixer for those variables inside the CP2K SCF cycle. CP2K "// &
                          "maps the variables onto CP2K's charge-mixing infrastructure. When TBLITE is active, the "// &
                          "native SCC convergence follows XTB/TBLITE/ACCURACY and the tblite-side SCC update limit "// &
                          "follows TBLITE_MIXER/ITERATIONS. If SCC_MIXER is explicitly set to TBLITE, CP2K also "// &
                          "sets DFT/SCF/MAX_SCF to TBLITE_MIXER/ITERATIONS. With SCC_MIXER AUTO or CP2K, MAX_SCF "// &
                          "keeps its normal CP2K meaning. With QS/LS_SCF, explicit SCC_MIXER settings are reset to "// &
                          "NONE internally and CP2K emits a warning; LS_SCF controls the density-matrix optimization. "// &
                          "When CP2K is active, convergence follows the normal CP2K SCF "// &
                          "controls, in particular DFT/SCF/EPS_SCF and DFT/SCF/MIXING. For GFN2/tblite, the CP2K "// &
                          "mixer acts on the full SCC-variable set: shell charges plus atomic dipole and quadrupole "// &
                          "variables. CP2K-internal GFN0-xTB has no SCC variables to mix, so AUTO and NONE are "// &
                          "equivalent there and explicit CP2K or TBLITE mixer choices are reset to NONE. For "// &
                          "CP2K-internal GFN1 and DFTB/DFTB3, the SCC variables are only charge or shell-charge like. The "// &
                          "initial SCC variables are seeded from CP2K's current density, so DFT/SCF/SCF_GUESS and "// &
                          "DFT/WFN_RESTART_FILE_NAME define the starting point. NONE performs direct updates without "// &
                          "mixing and is intended for debugging, since it can be unstable.", &
                          usage="SCC_MIXER (AUTO|TBLITE|CP2K|NONE)", &
                          enum_c_vals=s2a("AUTO", "TBLITE", "CP2K", "NONE"), &
                          enum_i_vals=[tblite_scc_mixer_auto, tblite_scc_mixer_tblite, &
                                       tblite_scc_mixer_cp2k, tblite_scc_mixer_none], &
                          enum_desc=s2a("Default method-dependent choice.", &
                                        "Use the tblite modified Broyden SCC mixer.", &
                                        "Use CP2K charge mixing for SCC variables.", &
                                        "Do not mix SCC variables; for debugging."), &
                          default_i_val=tblite_scc_mixer_auto)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE add_tb_scc_mixer_keywords

! **************************************************************************************************
!> \brief Creates the &REFERENCE_CLI section used to compare against native tblite.
!> \param section the section to create
! **************************************************************************************************
   SUBROUTINE create_xtb_reference_cli_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="REFERENCE_CLI", &
                          description="Run the native tblite command line interface as a reference "// &
                          "for CP2K/tblite energies, nuclear gradients, and optionally virials. "// &
                          "The section is a lone section: if it is present, the reference check is enabled. By default "// &
                          "the native CLI calculation inherits the CP2K/tblite method, PARAM file, ACCURACY, charge, "// &
                          "spin state, spin-polarized Hamiltonian setting, electronic temperature, "// &
                          "TBLITE_MIXER/ITERATIONS, and TBLITE_MIXER/SOLVER from the surrounding input. CP2K forwards "// &
                          "DFT/UKS, including its LSD and SPIN_POLARIZED aliases, as native tblite "// &
                          "--spin-polarized; DFT/MULTIPLICITY is forwarded separately as --spin. "// &
                          "REFERENCE_CLI/EFIELD can add a homogeneous electric field to the "// &
                          "external native tblite command only. REFERENCE_CLI/IMPLICIT_SOLVATION can add native "// &
                          "tblite implicit-solvent settings to the external reference command only. "// &
                          "Native tblite run does not expose independent options for mixer MEMORY, DAMPING, OMEGA0, "// &
                          "MIN_WEIGHT, MAX_WEIGHT, or WEIGHT_FACTOR; CP2K warns if these inherited production settings "// &
                          "cannot be reproduced and aborts when STOP_ON_ERROR is enabled. "// &
                          "Native CLI verbosity follows GLOBAL/PRINT_LEVEL: SILENT is forwarded as --silent, and "// &
                          "HIGH or DEBUG are forwarded as --verbose. "// &
                          "It remains a diagnostic one-shot calculation from temporary input files and is independent "// &
                          "of CP2K's SCF history, diagonalization, and k-point machinery. If RESTART "// &
                          "is set, the named native tblite restart file is passed to the external command. Non-Fermi "// &
                          "CP2K smearing methods cannot be represented by native tblite and are replaced by a "// &
                          "Fermi-Dirac electronic temperature in the reference command with an explicit warning.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Activates the native tblite CLI reference check when the section is present. "// &
                          "The check compares CP2K/tblite against a freshly executed external tblite calculation.", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PROGRAM_NAME", &
                          description="Executable name or path for the native tblite command line program used for "// &
                          "the reference calculation. This allows testing the CP2K/tblite implementation against a "// &
                          "specific tblite build without changing the CP2K executable.", &
                          usage="PROGRAM_NAME tblite", n_var=1, type_of_var=char_t, default_c_val="tblite")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GUESS", &
                          description="Initial population guess used only for the native tblite reference command. "// &
                          "SAD uses the superposition of atomic densities, matching the default of tblite run. EEQ "// &
                          "uses electronegativity equilibration charges, while CEH uses the charge-equilibration "// &
                          "Hamiltonian guess. This keyword does not affect the production CP2K/tblite calculation; "// &
                          "there, initial densities and restarts remain controlled by CP2K's SCF machinery.", &
                          usage="GUESS (SAD|EEQ|CEH)", &
                          enum_c_vals=s2a("SAD", "EEQ", "CEH"), &
                          enum_i_vals=[tblite_guess_sad, tblite_guess_eeq, tblite_guess_ceh], &
                          enum_desc=s2a("Use the superposition of atomic densities guess.", &
                                        "Use the electronegativity equilibration guess.", &
                                        "Use the charge-equilibration Hamiltonian guess."), &
                          default_i_val=tblite_guess_sad)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WORK_DIRECTORY", &
                          description="Directory used for the temporary native tblite CLI input and output files. "// &
                          "The CP2K/tblite result is not taken from this directory; the files are only used to run and "// &
                          "inspect the independent reference calculation.", &
                          usage="WORK_DIRECTORY .", n_var=1, type_of_var=char_t, default_c_val=".")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PREFIX", &
                          description="Prefix for the temporary native tblite CLI input and output files. Changing "// &
                          "the prefix is useful when several reference checks are run in the same directory or when "// &
                          "KEEP_FILES is used to inspect the generated geometry, gradient, JSON, or log files.", &
                          usage="PREFIX tblite-reference", n_var=1, type_of_var=char_t, &
                          default_c_val="tblite-reference")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INPUT_FORMAT", &
                          description="Input format hint passed to native tblite run as --input. The default is GEN, "// &
                          "matching the temporary DFTB+ GEN geometry written by CP2K for the reference calculation. "// &
                          "Only change this together with a compatible externally inspected or modified input file "// &
                          "workflow; the normal CP2K-generated reference geometry is GEN.", &
                          usage="INPUT_FORMAT gen", n_var=1, type_of_var=char_t, default_c_val="gen")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTART", &
                          description="Native tblite restart file passed to the external reference command. When this "// &
                          "keyword is absent, CP2K passes --no-restart so each diagnostic reference is a clean one-shot "// &
                          "run from the generated geometry and selected CLI options. When a filename is given, CP2K "// &
                          "passes --restart with that file. This affects only the external native tblite run; it does "// &
                          "not read or write CP2K wavefunction or density restart data and is therefore independent of "// &
                          "DFT/SCF/SCF_GUESS RESTART.", &
                          usage="RESTART filename", n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="POST_PROCESSING", &
                          description="Additional native tblite post-SCF analysis passed as tblite run "// &
                          "--post-processing. The value can be a built-in post-processing module name or a TOML file. "// &
                          "Native tblite already enables Mayer-Wiberg bond-order analysis internally by default. "// &
                          "The additional built-in module names currently exposed by tblite include molmom for "// &
                          "molecular moments and xtbml for xTB-ML descriptors. If the value names a TOML file, tblite "// &
                          "reads its [post-processing] configuration and can select the post-processing modules more "// &
                          "finely. This keyword affects only the external REFERENCE_CLI diagnostic run; CP2K does not "// &
                          "parse the generated post-processing data for the energy/force/virial comparison.", &
                          usage="POST_PROCESSING molmom", n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="POST_PROCESSING_OUTPUT", &
                          description="Filename passed to native tblite as --post-processing-output. Native tblite "// &
                          "documents this as the NPZ output file for post-processing data, with tblite-data.npz as its "// &
                          "CLI default. In the locally checked tblite 0.5.0 source this option appears in the help text "// &
                          "but is not wired in the command-line parser, so CP2K emits it only when explicitly requested "// &
                          "and prints a warning before running the reference command. Use this keyword only with tblite "// &
                          "CLI builds that actually accept the option. Unless KEEP_FILES is enabled, CP2K removes this "// &
                          "explicit output file after the reference check.", &
                          usage="POST_PROCESSING_OUTPUT tblite-data.npz", n_var=1, type_of_var=char_t, &
                          default_c_val="", lone_keyword_c_val="tblite-data.npz")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GRAD", &
                          description="Filename passed to native tblite as --grad for the tagged ASCII gradient output. "// &
                          "The native file contains the total energy, energy contributions, analytical molecular "// &
                          "gradient, and virial. CP2K reads this file for the REFERENCE_CLI energy, force, and virial "// &
                          "checks. If the keyword is omitted, CP2K writes a temporary PREFIX.grad file in "// &
                          "WORK_DIRECTORY. If the keyword is present without a filename, CP2K uses tblite's CLI default "// &
                          "name tblite.txt. Unless KEEP_FILES is enabled, CP2K removes the selected file after the "// &
                          "reference check.", &
                          usage="GRAD tblite.txt", n_var=1, type_of_var=char_t, default_c_val="", &
                          lone_keyword_c_val="tblite.txt")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JSON", &
                          description="Filename passed to native tblite as --json for JSON result output. Without "// &
                          "--grad the native JSON usually contains the energy and energy contributions; because CP2K's "// &
                          "REFERENCE_CLI path also requests --grad, native tblite can additionally include derivative "// &
                          "data such as gradient and virial when supported by the CLI version. CP2K currently reads the "// &
                          "tagged GRAD file for comparisons and keeps the JSON file as an inspectable native artifact. "// &
                          "If the keyword is omitted, CP2K writes a temporary PREFIX.json file in WORK_DIRECTORY. If the "// &
                          "keyword is present without a filename, CP2K uses tblite's CLI default name tblite.json. "// &
                          "Unless KEEP_FILES is enabled, CP2K removes the selected file after the reference check.", &
                          usage="JSON tblite.json", n_var=1, type_of_var=char_t, default_c_val="", &
                          lone_keyword_c_val="tblite.json")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ELECTRONIC_TEMPERATURE_GUESS", &
                          description="Electronic temperature used only for the native tblite reference command's "// &
                          "initial CEH population guess. CP2K forwards positive values as tblite run "// &
                          "--etemp-guess. The value does not affect CP2K/tblite production calculations and is "// &
                          "independent of DFT/SCF/SMEAR/ELECTRONIC_TEMPERATURE, which controls the final electronic "// &
                          "occupation temperature passed as --etemp. A value of zero leaves the native tblite default "// &
                          "for the guess temperature untouched. Positive values are valid only together with GUESS CEH; "// &
                          "CP2K aborts for SAD or EEQ because native tblite uses this option for the CEH guess.", &
                          usage="ELECTRONIC_TEMPERATURE_GUESS [K] 4000", &
                          default_r_val=cp_unit_to_cp2k(value=0.0_dp, unit_str="K"), unit_str="K")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EFIELD", &
                          description="Homogeneous electric field passed only to the native tblite reference command "// &
                          "as tblite run --efield. The three components are given in V/Angstrom, matching the native "// &
                          "tblite CLI convention. This keyword does not affect the production CP2K/tblite "// &
                          "calculation because that path currently does not construct tblite's external electric "// &
                          "field container. If enabled together with energy, force, or virial checks, the comparison "// &
                          "therefore intentionally diagnoses the native field response rather than matching a CP2K "// &
                          "production field calculation.", &
                          usage="EFIELD 0.0 0.0 0.0", repeats=.FALSE., n_var=3, type_of_var=real_t, &
                          default_r_vals=[0.0_dp, 0.0_dp, 0.0_dp])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="KEEP_FILES", &
                          description="Keep the temporary native tblite CLI input and output files after the "// &
                          "comparison. This is mainly a debugging aid for inspecting the generated GEN geometry, the "// &
                          "native tblite gradient/virial file, the JSON output, and the CLI log.", &
                          usage="KEEP_FILES T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ERROR_LIMIT", &
                          description="Tolerance for reporting reference CLI deviations. The value is applied to the "// &
                          "absolute CP2K-minus-native differences printed by the enabled energy, force, and virial "// &
                          "checks. Exceeding this limit produces a diagnostic warning, or an abort when STOP_ON_ERROR "// &
                          "is enabled.", &
                          usage="ERROR_LIMIT 1.0E-8", default_r_val=1.0E-8_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STOP_ON_ERROR", &
                          description="Abort the CP2K run if any enabled reference CLI comparison exceeds "// &
                          "ERROR_LIMIT. When disabled, CP2K prints the deviation but continues, which is useful for "// &
                          "exploratory diagnostics and regression-test probes.", &
                          usage="STOP_ON_ERROR T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHECK_ENERGY", &
                          description="Compare the CP2K/tblite total xTB energy with the energy obtained from the "// &
                          "native tblite reference run. This checks the Hamiltonian setup, SCC solution, static "// &
                          "repulsion/dispersion terms, charge and spin state, electronic temperature, and optional "// &
                          "PARAM file consistency for the one-shot reference geometry.", &
                          usage="CHECK_ENERGY T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHECK_FORCES", &
                          description="Compare CP2K/tblite nuclear gradients with the gradients written by the "// &
                          "native tblite reference run. This is a sensitive check of derivative terms, including "// &
                          "repulsion, dispersion, electrostatic response, and coordination-number response. The "// &
                          "comparison uses the same geometry and model setup as the energy reference.", &
                          usage="CHECK_FORCES T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHECK_VIRIAL", &
                          description="Compare CP2K/tblite virial with the virial written by the native tblite "// &
                          "reference run. This is a diagnostic for cell derivatives and stress-related terms. It is "// &
                          "disabled by default because virial conventions, finite periodicity, and the absence of "// &
                          "native CLI k-point sampling make this comparison more restrictive than the energy and "// &
                          "force checks.", &
                          usage="CHECK_VIRIAL T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL create_xtb_reference_cli_guess_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_xtb_reference_cli_param_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_xtb_reference_cli_fit_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_xtb_reference_cli_tagdiff_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_xtb_reference_cli_solvation_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_xtb_reference_cli_section

! **************************************************************************************************
!> \brief Creates the &REFERENCE_CLI%GUESS_CLI section for native tblite guess.
!> \param section section to create
! **************************************************************************************************
   SUBROUTINE create_xtb_reference_cli_guess_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="GUESS_CLI", &
                          description="Run the native tblite guess subcommand as an auxiliary REFERENCE_CLI "// &
                          "diagnostic. The command executes tblite guess on the generated GEN geometry by default and "// &
                          "prints the initializing charges and dipole moments. CP2K forwards the surrounding charge "// &
                          "and multiplicity as native --charge and --spin. GLOBAL/PRINT_LEVEL controls native "// &
                          "--silent or --verbose exactly as for the main REFERENCE_CLI run command. This section is "// &
                          "diagnostic only; CP2K does not use the generated guess data to seed the production "// &
                          "calculation.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Activates the native tblite guess auxiliary command when the section is "// &
                          "present.", lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="METHOD", &
                          description="Initial-population method passed to native tblite guess as --method. SAD uses "// &
                          "superposed atomic densities, EEQ uses electronegativity equilibration, and CEH uses the "// &
                          "charge-extended Huckel guess. The default matches native tblite guess.", &
                          usage="METHOD (SAD|EEQ|CEH)", &
                          enum_c_vals=s2a("SAD", "EEQ", "CEH"), &
                          enum_i_vals=[tblite_guess_sad, tblite_guess_eeq, tblite_guess_ceh], &
                          enum_desc=s2a("Use the superposition of atomic densities guess.", &
                                        "Use the electronegativity equilibration guess.", &
                                        "Use the charge-extended Huckel guess."), &
                          default_i_val=tblite_guess_ceh)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ELECTRONIC_TEMPERATURE_GUESS", &
                          description="Electronic temperature for native tblite guess --etemp-guess. This option is "// &
                          "used by the CEH guess in native tblite. A value of zero leaves tblite's own default "// &
                          "untouched.", usage="ELECTRONIC_TEMPERATURE_GUESS [K] 4000", &
                          default_r_val=cp_unit_to_cp2k(value=0.0_dp, unit_str="K"), unit_str="K")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SOLVER", &
                          description="Electronic solver passed to native tblite guess as --solver. GVD is the "// &
                          "default generalized eigenvalue decomposition solver; GVR selects the alternative solver "// &
                          "exposed by tblite.", usage="SOLVER (GVD|GVR)", &
                          enum_c_vals=s2a("GVD", "GVR"), enum_i_vals=[tblite_solver_gvd, tblite_solver_gvr], &
                          enum_desc=s2a("Use the GVD electronic solver.", "Use the GVR electronic solver."), &
                          default_i_val=tblite_solver_gvd)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EFIELD", &
                          description="Homogeneous electric field passed to native tblite guess as --efield. The "// &
                          "three components are given in V/Angstrom, matching native tblite.", &
                          usage="EFIELD 0.0 0.0 0.0", repeats=.FALSE., n_var=3, type_of_var=real_t, &
                          default_r_vals=[0.0_dp, 0.0_dp, 0.0_dp])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GRAD", &
                          description="Request native tblite guess --grad, evaluating the analytic gradient of the "// &
                          "guessed charges.", usage="GRAD T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JSON", &
                          description="Filename passed to native tblite guess as --json. If the keyword is present "// &
                          "without a filename, CP2K uses tblite's CLI default name tblite.json. If omitted, CP2K does "// &
                          "not request JSON output for the auxiliary guess command.", &
                          usage="JSON tblite.json", n_var=1, type_of_var=char_t, default_c_val="", &
                          lone_keyword_c_val="tblite.json")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INPUT_FORMAT", &
                          description="Input format hint passed to native tblite guess as --input. The default is GEN, "// &
                          "matching the temporary DFTB+ GEN geometry written by CP2K for the enclosing reference "// &
                          "calculation.", &
                          usage="INPUT_FORMAT gen", n_var=1, type_of_var=char_t, default_c_val="gen")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INPUT_FILE", &
                          description="Positional input structure for native tblite guess. The file format must match "// &
                          "GUESS_CLI/INPUT_FORMAT. If omitted, CP2K uses the temporary GEN geometry written for the "// &
                          "enclosing REFERENCE_CLI calculation.", &
                          usage="INPUT_FILE structure.gen", n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_xtb_reference_cli_guess_section

! **************************************************************************************************
!> \brief Creates the &REFERENCE_CLI%PARAM_CLI section for native tblite param.
!> \param section section to create
! **************************************************************************************************
   SUBROUTINE create_xtb_reference_cli_param_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PARAM_CLI", &
                          description="Run the native tblite param subcommand as an auxiliary REFERENCE_CLI "// &
                          "diagnostic. The command inspects, imports, exports, or converts tblite parameter data. "// &
                          "If no INPUT_FILE is supplied, CP2K selects the surrounding XTB/TBLITE/METHOD as native "// &
                          "--method. If METHOD is specified here, it overrides only this auxiliary command.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Activates the native tblite param auxiliary command when the section is "// &
                          "present.", lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="METHOD", &
                          description="Base parametrization passed to native tblite param as --method. If omitted, "// &
                          "CP2K uses the surrounding XTB/TBLITE/METHOD.", &
                          usage="METHOD (GFN1|GFN2|IPEA1)", &
                          enum_c_vals=s2a("GFN1", "GFN2", "IPEA1"), &
                          enum_i_vals=[gfn1xtb, gfn2xtb, ipea1xtb], &
                          enum_desc=s2a("Use the GFN1-xTB base parametrization.", &
                                        "Use the GFN2-xTB base parametrization.", &
                                        "Use the IPEA1 base parametrization."), &
                          default_i_val=gfn2xtb)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OUTPUT", &
                          description="Output parameter file passed to native tblite param as --output. If omitted, "// &
                          "native tblite writes to its default stream or file for the selected operation.", &
                          usage="OUTPUT parameters.toml", n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INPUT_FILE", &
                          description="Optional positional input parameter file for native tblite param. If omitted, "// &
                          "the command operates on the selected built-in METHOD.", &
                          usage="INPUT_FILE parameters.toml", n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_xtb_reference_cli_param_section

! **************************************************************************************************
!> \brief Creates the &REFERENCE_CLI%FIT_CLI section for native tblite fit.
!> \param section section to create
! **************************************************************************************************
   SUBROUTINE create_xtb_reference_cli_fit_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="FIT_CLI", &
                          description="Run the native tblite fit subcommand as an auxiliary REFERENCE_CLI diagnostic. "// &
                          "The command optimizes or dry-runs tblite parameters from a starting PARAM_FILE and a fit "// &
                          "settings INPUT_FILE. Because native fitting can be expensive and writes its own outputs, "// &
                          "the section is inactive unless explicitly present.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Activates the native tblite fit auxiliary command when the section is present.", &
                          lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARAM_FILE", &
                          description="Starting parameter file passed as the first positional argument to native "// &
                          "tblite fit.", usage="PARAM_FILE start.toml", n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INPUT_FILE", &
                          description="Fit-settings input file passed as the second positional argument to native "// &
                          "tblite fit.", usage="INPUT_FILE fit.toml", n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DRY_RUN", &
                          description="Pass native tblite fit --dry-run. This expands and checks the fit setup "// &
                          "without starting the optimization.", usage="DRY_RUN T", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COPY", &
                          description="Filename passed to native tblite fit as --copy. Native tblite writes the full "// &
                          "representation of the fit input to this file, with defaults filled in and masks expanded.", &
                          usage="COPY expanded-fit.toml", n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_xtb_reference_cli_fit_section

! **************************************************************************************************
!> \brief Creates the &REFERENCE_CLI%TAGDIFF_CLI section for native tblite tagdiff.
!> \param section section to create
! **************************************************************************************************
   SUBROUTINE create_xtb_reference_cli_tagdiff_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="TAGDIFF_CLI", &
                          description="Run the native tblite tagdiff auxiliary command. It compares an ACTUAL tagged "// &
                          "data file against a REFERENCE tagged data file and reports differences using the keys from "// &
                          "the reference file.", n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Activates the native tblite tagdiff auxiliary command when the section is "// &
                          "present.", lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ACTUAL", &
                          description="Actual tagged data file passed as the first positional argument to native "// &
                          "tblite tagdiff.", usage="ACTUAL actual.tag", n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE", &
                          description="Reference tagged data file passed as the second positional argument to native "// &
                          "tblite tagdiff. Native tblite uses only keys present in this reference file.", &
                          usage="REFERENCE reference.tag", n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FIT", &
                          description="Pass native tblite tagdiff --fit, producing output suitable for tblite fit.", &
                          usage="FIT T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_xtb_reference_cli_tagdiff_section

! **************************************************************************************************
!> \brief Creates the &REFERENCE_CLI%IMPLICIT_SOLVATION section for native tblite.
!> \param section section to create
! **************************************************************************************************
   SUBROUTINE create_xtb_reference_cli_solvation_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="IMPLICIT_SOLVATION", &
                          description="Native tblite implicit-solvation settings for the external REFERENCE_CLI "// &
                          "diagnostic command. This section does not affect the production CP2K/tblite calculation. "// &
                          "The section is active when present; no additional enable keyword is needed. "// &
                          "CP2K aborts for fully three-dimensional periodic cells and prints a warning for one- or "// &
                          "two-dimensional periodic cells because native tblite's implicit-solvent models are mainly "// &
                          "molecular continuum-solvation diagnostics. "// &
                          "MODEL selects the native tblite run solvation switch: ALPB maps to --alpb, GBSA to --gbsa, "// &
                          "GBE to --gbe, GB to --gb, and CPCM to --cpcm. SOLVENT is the argument passed to that switch. "// &
                          "For ALPB and GBSA it must be a named solvent with tblite parameters; for GBE, GB, and CPCM "// &
                          "it can be either a named solvent or a dielectric constant written as a token such as 78.36. "// &
                          "BORN_KERNEL maps to --born-kernel and applies only to ALPB, GBSA, GBE, and GB; AUTO leaves "// &
                          "the native tblite model-dependent default untouched. SOLUTION_STATE maps to --solv-state "// &
                          "and is only supported for the parametrized named solvent ALPB/GBSA models; GSOLV leaves "// &
                          "the native tblite default untouched.", &
                          n_keywords=4, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="MODEL", &
                          description="Implicit-solvation model passed to native tblite run. ALPB selects the "// &
                          "parametrized analytical linearized Poisson-Boltzmann model with CDS and solution-state "// &
                          "shift support. GBSA selects the parametrized generalized Born surface-area model with CDS "// &
                          "and solution-state shift support. GBE selects finite-epsilon generalized Born/ALPB "// &
                          "electrostatics without the named-solvent CDS and shift terms. GB selects generalized Born "// &
                          "electrostatics without the named-solvent CDS and shift terms. CPCM selects tblite's "// &
                          "conductor-like polarizable continuum model.", &
                          usage="MODEL (ALPB|GBSA|GBE|GB|CPCM)", &
                          enum_c_vals=s2a("ALPB", "GBSA", "GBE", "GB", "CPCM"), &
                          enum_i_vals=[tblite_cli_solvation_alpb, tblite_cli_solvation_gbsa, &
                                       tblite_cli_solvation_gbe, tblite_cli_solvation_gb, &
                                       tblite_cli_solvation_cpcm], &
                          enum_desc=s2a("Use native tblite --alpb.", &
                                        "Use native tblite --gbsa.", &
                                        "Use native tblite --gbe.", &
                                        "Use native tblite --gb.", &
                                        "Use native tblite --cpcm."), &
                          default_i_val=tblite_cli_solvation_alpb)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SOLVENT", &
                          description="Solvent argument passed to the selected native tblite solvation switch. "// &
                          "For MODEL ALPB and MODEL GBSA this must be a named solvent supported by tblite's "// &
                          "parametrized solvation data, for example water, acetone, methanol, thf, or toluene. "// &
                          "For MODEL GBE, MODEL GB, and MODEL CPCM this can also be a dielectric constant, for "// &
                          "example 78.36. Native tblite normalizes names by ignoring case, spaces, commas, and "// &
                          "hyphens before lookup.", &
                          usage="SOLVENT water", n_var=1, type_of_var=char_t, default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BORN_KERNEL", &
                          description="Born interaction kernel passed to native tblite as --born-kernel. AUTO leaves "// &
                          "the native default untouched: tblite uses p16 for ALPB/GBE and still for GB/GBSA. P16 "// &
                          "selects the Lange p16 kernel, while STILL selects the classical Still kernel. This keyword "// &
                          "is not valid with MODEL CPCM because CPCM does not use the Born-kernel machinery.", &
                          usage="BORN_KERNEL (AUTO|P16|STILL)", &
                          enum_c_vals=s2a("AUTO", "P16", "STILL"), &
                          enum_i_vals=[tblite_cli_born_kernel_auto, tblite_cli_born_kernel_p16, &
                                       tblite_cli_born_kernel_still], &
                          enum_desc=s2a("Use tblite's model-dependent default kernel.", &
                                        "Use the Lange p16 Born kernel.", &
                                        "Use the classical Still Born kernel."), &
                          default_i_val=tblite_cli_born_kernel_auto)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SOLUTION_STATE", &
                          description="Solution-state convention passed to native tblite as --solv-state. GSOLV is "// &
                          "the default solvation free-energy convention. BAR1MOL applies the 1 bar ideal-gas to "// &
                          "1 mol/L liquid-solution shift. REFERENCE requests tblite's reference-state convention. "// &
                          "Native tblite supports non-default solution-state shifts only for named parametrized "// &
                          "ALPB/GBSA solvents, not for GBE, GB, or CPCM.", &
                          usage="SOLUTION_STATE (GSOLV|BAR1MOL|REFERENCE)", &
                          enum_c_vals=s2a("GSOLV", "BAR1MOL", "REFERENCE"), &
                          enum_i_vals=[tblite_cli_solution_state_gsolv, tblite_cli_solution_state_bar1mol, &
                                       tblite_cli_solution_state_reference], &
                          enum_desc=s2a("Use tblite's default solvation free-energy convention.", &
                                        "Apply the 1 bar ideal-gas to 1 mol/L solution shift.", &
                                        "Use tblite's reference-state convention."), &
                          default_i_val=tblite_cli_solution_state_gsolv)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_xtb_reference_cli_solvation_section

END MODULE input_cp2k_tb
