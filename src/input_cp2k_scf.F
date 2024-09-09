!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the scf section of the input
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!>      07.2024 moved out of input_cp2k_dft [JGH]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_scf
   USE bibliography,                    ONLY: Becke1988b,&
                                              Holmberg2017,&
                                              Holmberg2018,&
                                              Schiffmann2015,&
                                              Stewart1982,&
                                              VandeVondele2003,&
                                              VandeVondele2005a,&
                                              Weber2008
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              high_print_level,&
                                              low_print_level
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE input_constants,                 ONLY: &
        atomic_guess, becke_cutoff_element, becke_cutoff_global, broyden_type_1, &
        broyden_type_1_explicit, broyden_type_1_explicit_ls, broyden_type_1_ls, broyden_type_2, &
        broyden_type_2_explicit, broyden_type_2_explicit_ls, broyden_type_2_ls, &
        cdft_alpha_constraint, cdft_beta_constraint, cdft_charge_constraint, &
        cdft_magnetization_constraint, cholesky_dbcsr, cholesky_inverse, cholesky_off, &
        cholesky_reduce, cholesky_restore, core_guess, diag_block_davidson, diag_block_krylov, &
        diag_filter_matrix, diag_ot, diag_standard, gaussian, general_roks, high_spin_roks, &
        history_guess, jacobian_fd1, jacobian_fd1_backward, jacobian_fd1_central, jacobian_fd2, &
        jacobian_fd2_backward, ls_2pnt, ls_3pnt, ls_gold, ls_none, mopac_guess, no_guess, &
        numerical, ot_algo_irac, ot_algo_taylor_or_diag, ot_chol_irac, ot_lwdn_irac, &
        ot_mini_broyden, ot_mini_cg, ot_mini_diis, ot_mini_sd, ot_poly_irac, ot_precond_full_all, &
        ot_precond_full_kinetic, ot_precond_full_single, ot_precond_full_single_inverse, &
        ot_precond_none, ot_precond_s_inverse, ot_precond_solver_default, &
        ot_precond_solver_direct, ot_precond_solver_inv_chol, ot_precond_solver_update, &
        outer_scf_basis_center_opt, outer_scf_becke_constraint, outer_scf_cdft_constraint, &
        outer_scf_ddapc_constraint, outer_scf_hirshfeld_constraint, outer_scf_none, &
        outer_scf_optimizer_bisect, outer_scf_optimizer_broyden, outer_scf_optimizer_diis, &
        outer_scf_optimizer_newton, outer_scf_optimizer_newton_ls, outer_scf_optimizer_none, &
        outer_scf_optimizer_sd, outer_scf_optimizer_secant, outer_scf_s2_constraint, &
        radius_covalent, radius_default, radius_single, radius_user, radius_vdw, random_guess, &
        restart_guess, shape_function_density, shape_function_gaussian, smear_energy_window, &
        smear_fermi_dirac, smear_list, sparse_guess
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: integer_t,&
                                              real_t
   USE kinds,                           ONLY: dp
   USE qs_density_mixing_types,         ONLY: create_mixing_section
   USE qs_fb_input,                     ONLY: create_filtermatrix_section
   USE qs_mom_types,                    ONLY: create_mom_section
   USE string_utilities,                ONLY: newline,&
                                              s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_scf'

   PUBLIC :: create_scf_section, create_cdft_control_section

CONTAINS

! **************************************************************************************************
!> \brief creates the structure of the section with the DFT SCF parameters
!> \param section will contain the SCF section
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE create_scf_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      NULLIFY (print_key)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="scf", &
                          description="Parameters needed to perform an SCF run.", &
                          n_keywords=18, n_subsections=7, repeats=.FALSE.)

      NULLIFY (subsection)

      CALL create_ot_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_diagonalization_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_outer_scf_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_smear_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_mixing_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_mom_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER_LUMO", &
                          variants=(/"MAX_ITER_LUMOS"/), &
                          description="Maximum number of iterations for the calculation of the LUMO energies "// &
                          "with the OT eigensolver.", &
                          usage="MAX_ITER_LUMO 100", default_i_val=299)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_LUMO", &
                          variants=(/"EPS_LUMOS"/), &
                          description="Target accuracy for the calculation of the LUMO energies with the OT eigensolver.", &
                          usage="EPS_LUMO 1.0E-6", default_r_val=1.0E-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_SCF", &
                          description="Maximum number of SCF iteration to be performed for one optimization", &
                          usage="MAX_SCF 200", default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_SCF_HISTORY", variants=(/"MAX_SCF_HIST"/), &
                          description="Maximum number of SCF iterations after the history pipeline is filled", &
                          usage="MAX_SCF_HISTORY 1", default_i_val=0, lone_keyword_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_DIIS", &
                          variants=(/"MAX_DIIS_BUFFER_SIZE"/), &
                          description="Maximum number of DIIS vectors to be used", &
                          usage="MAX_DIIS 3", default_i_val=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LEVEL_SHIFT", &
                          variants=(/"LSHIFT"/), &
                          description="Use level shifting to improve convergence", &
                          unit_str="au_e", &
                          usage="LEVEL_SHIFT 0.1", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_SCF", &
                          description="Target accuracy for the SCF convergence.", &
                          usage="EPS_SCF 1.e-6", default_r_val=1.e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_SCF_HISTORY", variants=(/"EPS_SCF_HIST"/), &
                          description="Target accuracy for the SCF convergence after the history pipeline is filled.", &
                          usage="EPS_SCF_HISTORY 1.e-5", default_r_val=0.0_dp, lone_keyword_r_val=1.0e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHOLESKY", &
                          description="If the cholesky method should be used for computing "// &
                          "the inverse of S, and in this case calling which Lapack routines", &
                          usage="CHOLESKY REDUCE", default_i_val=cholesky_restore, &
                          enum_c_vals=s2a("OFF", "REDUCE", "RESTORE", "INVERSE", "INVERSE_DBCSR"), &
                          enum_desc=s2a("The cholesky algorithm is not used", "Reduce is called", &
                                        "Reduce is replaced by two restore", &
                                        "Restore uses operator multiply by inverse of the triangular matrix", &
                                        "Like inverse, but matrix stored as dbcsr, sparce matrix algebra used when possible"), &
                          enum_i_vals=(/cholesky_off, cholesky_reduce, cholesky_restore, cholesky_inverse, cholesky_dbcsr/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_EIGVAL", &
                          description="Throw away linear combinations of basis functions with a small eigenvalue in S", &
                          usage="EPS_EIGVAL 1.0", default_r_val=1.0e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_DIIS", &
                          description="Threshold on the convergence to start using DIAG/DIIS", &
                          usage="EPS_DIIS 5.0e-2", default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="SCF_GUESS", &
         description="Change the initial guess for the wavefunction.", &
         usage="SCF_GUESS RESTART", default_i_val=atomic_guess, &
         enum_c_vals=s2a("ATOMIC", "RESTART", "RANDOM", "CORE", &
                         "HISTORY_RESTART", "MOPAC", "SPARSE", "NONE"), &
         enum_desc=s2a("Generate an atomic density using the atomic code", &
                       "Use the RESTART file as an initial guess (and ATOMIC if not present).", &
                       "Use random wavefunction coefficients.", &
                       "Diagonalize the core hamiltonian for an initial guess.", &
                       "Extrapolated from previous RESTART files.", &
                       "Use same guess as MOPAC for semi-empirical methods or a simple diagonal density matrix for other methods", &
                       "Generate a sparse wavefunction using the atomic code (for OT based methods)", &
                       "Skip initial guess (only for non-self consistent methods)."), &
         enum_i_vals=(/atomic_guess, restart_guess, random_guess, core_guess, &
                       history_guess, mopac_guess, sparse_guess, no_guess/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NROW_BLOCK", &
                          description="sets the number of rows in a scalapack block", &
                          usage="NROW_BLOCK 31", default_i_val=32)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NCOL_BLOCK", &
                          description="Sets the number of columns in a scalapack block", &
                          usage="NCOL_BLOCK 31", default_i_val=32)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ADDED_MOS", &
                          description="Number of additional MOS added for each spin. Use -1 to add all available. "// &
                          "alpha/beta spin can be specified independently (if spin-polarized calculation requested).", &
                          usage="ADDED_MOS", default_i_val=0, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="ROKS_SCHEME", &
                          description="Selects the ROKS scheme when ROKS is applied.", &
                          usage="ROKS_SCHEME HIGH-SPIN", &
                          repeats=.FALSE., &
                          n_var=1, &
                          enum_c_vals=s2a("GENERAL", "HIGH-SPIN"), &
                          enum_i_vals=(/general_roks, high_spin_roks/), &
                          default_i_val=high_spin_roks)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="ROKS_F", &
                          variants=(/"F_ROKS"/), &
                          description="Allows to define the parameter f for the "// &
                          "general ROKS scheme.", &
                          usage="ROKS_PARAMETER 1/2", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="ROKS_PARAMETERS", &
                          variants=(/"ROKS_PARAMETER"/), &
                          description="Allows to define all parameters for the high-spin "// &
                          "ROKS scheme explicitly. "// &
                          "The full set of 6 parameters has to be specified "// &
                          "in the order acc, bcc, aoo, boo, avv, bvv", &
                          usage="ROKS_PARAMETERS 1/2 1/2 1/2 1/2 1/2 1/2", &
                          repeats=.FALSE., &
                          n_var=6, &
                          type_of_var=real_t, &
                          default_r_vals=(/-0.5_dp, 1.5_dp, 0.5_dp, 0.5_dp, 1.5_dp, -0.5_dp/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="IGNORE_CONVERGENCE_FAILURE", &
                          description="If true, only a warning is issued if an SCF "// &
                          "iteration has not converged. By default, a run is aborted "// &
                          "if the required convergence criteria have not been achieved.", &
                          usage="IGNORE_CONVERGENCE_FAILURE logical_value", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL section_create(subsection, __LOCATION__, name="PRINT", &
                          description="Printing of information during the SCF.", repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESTART", &
                                       description="Controls the dumping of the MO restart file during SCF. "// &
                                       "By default keeps a short history of three restarts. "// &
                                       "See also RESTART_HISTORY", &
                                       print_level=low_print_level, common_iter_levels=3, &
                                       each_iter_names=s2a("QS_SCF"), each_iter_values=(/20/), &
                                       add_last=add_last_numeric, filename="RESTART")
      CALL keyword_create(keyword, __LOCATION__, name="BACKUP_COPIES", &
                          description="Specifies the maximum number of backup copies.", &
                          usage="BACKUP_COPIES {int}", &
                          default_i_val=1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create( &
         print_key, __LOCATION__, "RESTART_HISTORY", &
         description="Dumps unique MO restart files during the run keeping all of them.", &
         print_level=low_print_level, common_iter_levels=0, &
         each_iter_names=s2a("__ROOT__", "MD", "GEO_OPT", "ROT_OPT", "NEB", "METADYNAMICS", "QS_SCF"), &
         each_iter_values=(/500, 500, 500, 500, 500, 500, 500/), &
         filename="RESTART")
      CALL keyword_create(keyword, __LOCATION__, name="BACKUP_COPIES", &
                          description="Specifies the maximum number of backup copies.", &
                          usage="BACKUP_COPIES {int}", &
                          default_i_val=1)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "iteration_info", &
                                       description="Controls the printing of basic iteration information during the SCF.", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, name="time_cumul", &
                          description="If the printkey is activated switches the printing of timings"// &
                          " to cumulative (over the SCF).", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "program_run_info", &
                                       description="Controls the printing of basic information during the SCF.", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "MO_ORTHONORMALITY", &
                                       description="Controls the printing relative to the orthonormality of MOs (CT S C).", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "MO_MAGNITUDE", &
                                       description="Prints the min/max eigenvalues of the overlap of the MOs without S (CT C).", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "detailed_energy", &
                                       description="Controls the printing of detailed energy information.", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "diis_info", &
                                       description="Controls the printing of DIIS information.", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "total_densities", &
                                       description="Controls the printing of total densities.", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "Lanczos", &
                                       description="Controls the printing of information on Lanczos refinement iterations.", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create( &
         print_key, __LOCATION__, "DIAG_SUB_SCF", &
         description="Controls the printing of information on subspace diagonalization internal loop. ", &
         print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "Davidson", &
                                       description="Controls the printing of information on Davidson iterations.", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "FILTER_MATRIX", &
                                       description="Controls the printing of information on Filter Matrix method.", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(subsection, print_key)
      CALL section_release(print_key)

      CALL keyword_create(keyword, __LOCATION__, name="DM_RESTART_WRITE", &
                          description="Write the density matrix into a binary file at the end of the SCF.", &
                          usage="DM_RESTART_WRITE", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_scf_section

! **************************************************************************************************
!> \brief creates the structure of the section with SCF parameters
!>      controlling an other loop
!> \param section will contain the SCF section
!> \author Joost VandeVondele [2006.03]
! **************************************************************************************************
   SUBROUTINE create_outer_scf_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="OUTER_SCF", &
                          description="parameters controlling the outer SCF loop", &
                          n_keywords=13, n_subsections=1, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of the outer SCF loop", &
                          usage="&OUTER_SCF ON", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! add CDFT_OPT section
      NULLIFY (subsection)
      CALL create_cdft_opt_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, __LOCATION__, name="TYPE", &
                          description="Specifies which kind of outer SCF should be employed", &
                          usage="TYPE DDAPC_CONSTRAINT ", &
                          default_i_val=outer_scf_none, &
                          enum_c_vals=s2a("DDAPC_CONSTRAINT", "S2_CONSTRAINT", &
                                          "BASIS_CENTER_OPT", "CDFT_CONSTRAINT", "NONE"), &
                          enum_desc=s2a("Enforce a constraint on the DDAPC, requires the corresponding section", &
                                        "Enforce a constraint on the S2, requires the corresponding section", &
                                        "Optimize positions of basis functions, if atom types FLOATING_BASIS_CENTER "// &
                                        "are defined", &
                                        "Enforce a constraint on a generic CDFT weight population. "// &
                                        "Requires the corresponding section QS&CDFT"// &
                                        " which determines the type of weight used.", &
                                        "Do nothing in the outer loop, useful for resetting the inner loop,"), &
                          enum_i_vals=(/outer_scf_ddapc_constraint, outer_scf_s2_constraint, &
                                        outer_scf_basis_center_opt, outer_scf_cdft_constraint, outer_scf_none/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OPTIMIZER", &
                          description="Method used to bring the outer loop to a stationary point", &
                          usage="OPTIMIZER SD", &
                          default_i_val=outer_scf_optimizer_none, &
                          enum_c_vals=s2a("SD", "DIIS", "NONE", "BISECT", "BROYDEN", "NEWTON", "SECANT", "NEWTON_LS"), &
                          enum_desc=s2a("Takes steps in the direction of the gradient, multiplied by step_size", &
                                        "Uses a Direct Inversion in the Iterative Subspace method", &
                                        "Do nothing, useful only with the none type", &
                                        "Bisection of the gradient, useful for difficult one dimensional cases", &
                                        "Broyden's method. Variant defined in BROYDEN_TYPE.", &
                                        "Newton's method. Only compatible with CDFT constraints.", &
                                        "Secant method. Only for one dimensional cases. See Broyden for "// &
                                        "multidimensional cases.", &
                                        "Newton's method with backtracking line search to find the optimal step size. "// &
                                        "Only compatible with CDFT constraints. Starts from the regular Newton solution "// &
                                        "and successively reduces the step size until the L2 norm of the CDFT gradient "// &
                                        "decreases or MAX_LS steps is reached. Potentially very expensive because "// &
                                        "each iteration performs a full SCF calculation."), &
                          enum_i_vals=(/outer_scf_optimizer_sd, outer_scf_optimizer_diis, outer_scf_optimizer_none, &
                                        outer_scf_optimizer_bisect, outer_scf_optimizer_broyden, &
                                        outer_scf_optimizer_newton, outer_scf_optimizer_secant, &
                                        outer_scf_optimizer_newton_ls/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BISECT_TRUST_COUNT", &
                          description="Maximum number of times the same point will be used in bisection,"// &
                          " a small number guards against the effect of wrongly converged states.", &
                          usage="BISECT_TRUST_COUNT 5", default_i_val=10)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_SCF", &
                          description="The target gradient of the outer SCF variables. "// &
                          "Notice that the EPS_SCF of the inner loop also determines "// &
                          "the value that can be reached in the outer loop, "// &
                          "typically EPS_SCF of the outer loop must be smaller "// &
                          "than or equal to EPS_SCF of the inner loop.", &
                          usage="EPS_SCF 1.0E-6 ", default_r_val=1.0E-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DIIS_BUFFER_LENGTH", &
                          description="Maximum number of DIIS vectors used ", &
                          usage="DIIS_BUFFER_LENGTH 5", default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EXTRAPOLATION_ORDER", &
                          description="Number of past states used in the extrapolation of the variables during e.g. MD", &
                          usage="EXTRAPOLATION_ORDER 5", default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_SCF", &
                          description="The maximum number of outer loops ", &
                          usage="MAX_SCF 20", default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STEP_SIZE", &
                          description="The initial step_size used in the optimizer (currently steepest descent). "// &
                          "Note that in cases where a sadle point is sought for (constrained DFT),"// &
                          " this can be negative. For Newton and Broyden optimizers, use a value less/higher than "// &
                          "the default 1.0 (in absolute value, the sign is not significant) to active an under/overrelaxed "// &
                          "optimizer.", &
                          usage="STEP_SIZE -1.0", default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_outer_scf_section

! **************************************************************************************************
!> \brief makes the orbital transformation section
!> \param section ...
!> \par History
!>      11.2004 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE create_ot_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="OT", &
                          description="Sets the various options for the orbital transformation (OT) method. "// &
                          "Default settings already provide an efficient, yet robust method. "// &
                          "Most systems benefit from using the FULL_ALL preconditioner "// &
                          "combined with a small value (0.001) of ENERGY_GAP. "// &
                          "Well-behaved systems might benefit from using a DIIS minimizer. "//newline//newline// &
                          "**Advantages:** "// &
                          "It's fast, because no expensive diagonalisation is performed. "// &
                          "If preconditioned correctly, method guaranteed to find minimum. "//newline//newline// &
                          "**Disadvantages:** "// &
                          "Sensitive to preconditioning. A good preconditioner can be expensive. "// &
                          "No smearing, or advanced SCF mixing possible: POOR convergence for metallic systems.", &
                          n_keywords=27, n_subsections=0, repeats=.FALSE., &
                          citations=(/VandeVondele2003, Weber2008/))

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of the ot method", &
                          usage="&OT T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ALGORITHM", &
                          description="Algorithm to be used for OT", &
                          usage="ALGORITHM STRICT", &
                          default_i_val=ot_algo_taylor_or_diag, &
                          enum_c_vals=s2a("STRICT", "IRAC"), &
                          enum_desc=s2a("Strict orthogonality: Taylor or diagonalization based algorithm.", &
                                        "Orbital Transformation based Iterative Refinement "// &
                                        "of the Approximative Congruence transformation (OT/IR)."), &
                          enum_i_vals=(/ot_algo_taylor_or_diag, ot_algo_irac/), &
                          citations=(/VandeVondele2003, VandeVondele2005a, Weber2008/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="IRAC_DEGREE", &
                          description="The refinement polynomial degree (2, 3 or 4).", &
                          usage="IRAC_DEGREE 4", &
                          default_i_val=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_IRAC", &
                          description="Maximum allowed refinement iteration.", &
                          usage="MAX_IRAC 5", &
                          default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ORTHO_IRAC", &
                          description="The orthogonality method.", &
                          usage="ORTHO_IRAC POLY", &
                          default_i_val=ot_chol_irac, &
                          enum_c_vals=s2a("CHOL", "POLY", "LWDN"), &
                          enum_desc=s2a("Cholesky.", "Polynomial.", "Loewdin."), &
                          enum_i_vals=(/ot_chol_irac, ot_poly_irac, ot_lwdn_irac/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_IRAC_FILTER_MATRIX", &
                          description="Sets the threshold for filtering the matrices.", &
                          usage="EPS_IRAC_FILTER_MATRIX 1.0E-5", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_IRAC", &
                          description="Targeted accuracy during the refinement iteration.", &
                          usage="EPS_IRAC 1.0E-5", &
                          default_r_val=1.0E-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_IRAC_QUICK_EXIT", &
                          description="Only one extra refinement iteration is "// &
                          "done when the norm is below this value.", &
                          usage="EPS_IRAC_QUICK_EXIT 1.0E-2", &
                          default_r_val=1.0E-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_IRAC_SWITCH", &
                          description="The algorithm switches to the polynomial "// &
                          "refinement when the norm is below this value.", &
                          usage="EPS_IRAC_SWITCH 1.0E-3", &
                          default_r_val=1.0E-2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ON_THE_FLY_LOC", &
                          description="On the fly localization of the molecular orbitals. "// &
                          "Can only be used with OT/IRAC.", &
                          usage="ON_THE_FLY_LOC T", &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="MINIMIZER", &
         description="Minimizer to be used with the OT method", &
         usage="MINIMIZER DIIS", &
         default_i_val=ot_mini_cg, &
         enum_c_vals=s2a("SD", "CG", "DIIS", "BROYDEN"), &
         enum_desc=s2a("Steepest descent: not recommended", "Conjugate Gradients: most reliable, use for difficult systems."// &
                       " The total energy should decrease at every OT CG step if the line search is appropriate.", &
                       "Direct inversion in the iterative subspace: less reliable than CG, but sometimes about 50% faster", &
                       "Broyden mixing approximating the inverse Hessian"), &
         enum_i_vals=(/ot_mini_sd, ot_mini_cg, ot_mini_diis, ot_mini_broyden/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SAFE_DIIS", &
                          variants=(/"SAFER_DIIS"/), &
                          description="Reject DIIS steps if they point away from the"// &
                          " minimum, do SD in that case.", &
                          usage="SAFE_DIIS ON", default_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="N_HISTORY_VEC", &
                          variants=s2a("NDIIS", "N_DIIS", "N_BROYDEN"), &
                          description="Number of history vectors to be used with DIIS or BROYDEN", &
                          usage="N_DIIS 4", &
                          default_i_val=7)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_BETA", &
                          description="Underrelaxation for the broyden mixer", &
                          usage="BROYDEN_BETA 0.9", &
                          default_r_val=0.9_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_GAMMA", &
                          description="Backtracking parameter", &
                          usage="BROYDEN_GAMMA 0.5", &
                          default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_SIGMA", &
                          description="Curvature of energy functional.", &
                          usage="BROYDEN_SIGMA 0.25", &
                          default_r_val=0.25_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_ETA", &
                          description="Dampening of estimated energy curvature.", &
                          usage="BROYDEN_ETA 0.7", &
                          default_r_val=0.7_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_OMEGA", &
                          description="Growth limit of curvature.", &
                          usage="BROYDEN_OMEGA 1.1", &
                          default_r_val=1.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_SIGMA_DECREASE", &
                          description="Reduction of curvature on bad approximation.", &
                          usage="BROYDEN_SIGMA_DECREASE 0.7", &
                          default_r_val=0.7_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_SIGMA_MIN", &
                          description="Minimum adaptive curvature.", &
                          usage="BROYDEN_SIGMA_MIN 0.05", &
                          default_r_val=0.05_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_FORGET_HISTORY", &
                          description="Forget history on bad approximation", &
                          usage="BROYDEN_FORGET_HISTORY OFF", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_ADAPTIVE_SIGMA", &
                          description="Enable adaptive curvature estimation", &
                          usage="BROYDEN_ADAPTIVE_SIGMA ON", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_ENABLE_FLIP", &
                          description="Ensure positive definite update", &
                          usage="BROYDEN_ENABLE_FLIP ON", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LINESEARCH", &
                          variants=(/"LINE_SEARCH"/), &
                          description="1D line search algorithm to be used with the OT minimizer,"// &
                          " in increasing order of robustness and cost. MINIMIZER CG combined with"// &
                          " LINESEARCH GOLD should always find an electronic minimum."// &
                          " Whereas the 2PNT minimizer is almost always OK, 3PNT might be needed for systems"// &
                          " in which successive OT CG steps do not decrease the total energy.", &
                          usage="LINESEARCH GOLD", &
                          default_i_val=ls_2pnt, &
                          enum_c_vals=s2a("NONE", "2PNT", "3PNT", "GOLD"), &
                          enum_desc=s2a("take fixed length steps", "extrapolate based on 2 points", &
                                        "... or on 3 points", "perform 1D golden section search of the minimum (very expensive)"), &
                          enum_i_vals=(/ls_none, ls_2pnt, ls_3pnt, ls_gold/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="STEPSIZE", &
         description="Initial stepsize used for the line search, sometimes this parameter can be reduced to stabilize DIIS"// &
         " or to improve the CG behavior in the first few steps."// &
         " The optimal value depends on the quality of the preconditioner."// &
         " A negative values leaves the choice to CP2K depending on the preconditioner.", &
         usage="STEPSIZE 0.4", &
         default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GOLD_TARGET", &
                          description="Target relative uncertainty in the location of the minimum for LINESEARCH GOLD", &
                          usage="GOLD_TARGET 0.1", &
                          default_r_val=0.01_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="PRECONDITIONER", &
         description="Type of preconditioner to be used with all minimization schemes. "// &
         "They differ in effectiveness, cost of construction, cost of application. "// &
         "Properly preconditioned minimization can be orders of magnitude faster than doing nothing.", &
         usage="PRECONDITIONER FULL_ALL", &
         default_i_val=ot_precond_full_kinetic, &
         enum_c_vals=s2a("FULL_ALL", "FULL_SINGLE_INVERSE", "FULL_SINGLE", "FULL_KINETIC", "FULL_S_INVERSE", &
                         "NONE"), &
         enum_desc=s2a("Most effective state selective preconditioner based on diagonalization, "// &
                       "requires the ENERGY_GAP parameter to be an underestimate of the HOMO-LUMO gap. "// &
                       "This preconditioner is recommended for almost all systems, except very large systems where "// &
                       "make_preconditioner would dominate the total computational cost.", &
                       "Based on H-eS cholesky inversion, similar to FULL_SINGLE in preconditioning efficiency "// &
                       "but cheaper to construct, "// &
                       "might be somewhat less robust. Recommended for large systems.", &
                       "Based on H-eS diagonalisation, not as good as FULL_ALL, but somewhat cheaper to apply. ", &
                       "Cholesky inversion of S and T, fast construction, robust, and relatively good, "// &
                       "use for very large systems.", &
                       "Cholesky inversion of S, not as good as FULL_KINETIC, yet equally expensive.", &
                       "skip preconditioning"), &
         enum_i_vals=(/ot_precond_full_all, ot_precond_full_single_inverse, ot_precond_full_single, &
                       ot_precond_full_kinetic, ot_precond_s_inverse, ot_precond_none/), &
         citations=(/VandeVondele2003, Weber2008, Schiffmann2015/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHOLESKY", &
                          description="If FULL_ALL the cholesky decomposition of the S matrix is used. "// &
                          "Options on the algorithm to be used.", &
                          usage="CHOLESKY REDUCE", default_i_val=cholesky_reduce, &
                          enum_c_vals=s2a("OFF", "REDUCE", "RESTORE", "INVERSE", "INVERSE_DBCSR"), &
                          enum_desc=s2a("The cholesky algorithm is not used", "Reduce is called", &
                                        "Reduce is replaced by two restore", &
                                        "Restore uses operator multiply by inverse of the triangular matrix", &
                                        "Like inverse, but matrix stored as dbcsr, sparce matrix algebra used when possible"), &
                          enum_i_vals=(/cholesky_off, cholesky_reduce, cholesky_restore, cholesky_inverse, cholesky_dbcsr/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="PRECOND_SOLVER", &
         description="How the preconditioner is applied to the residual.", &
         usage="PRECOND_SOLVER DIRECT", &
         default_i_val=ot_precond_solver_default, &
         enum_c_vals=s2a("DEFAULT", "DIRECT", "INVERSE_CHOLESKY", "INVERSE_UPDATE"), &
         enum_desc=s2a("the default", "Cholesky decomposition followed by triangular solve "// &
                       "(works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)", &
                       "Cholesky decomposition followed by explicit inversion "// &
                       "(works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)", &
                       "Performs a Hotelling update of the inverse if a previous preconditioner is present. "// &
                       "Mainly useful for GPU accelerated systems (works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)"), &
         enum_i_vals=(/ot_precond_solver_default, &
                       ot_precond_solver_direct, &
                       ot_precond_solver_inv_chol, &
                       ot_precond_solver_update/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="ENERGY_GAP", &
         description="Should be an estimate for the energy gap [a.u.] (HOMO-LUMO) and is used in preconditioning, "// &
         "especially effective with the FULL_ALL preconditioner, in which case it should be an underestimate "// &
         "of the gap (can be a small number, e.g. 0.002)."// &
         " FULL_SINGLE_INVERSE takes it as lower bound (values below 0.05 can cause stability issues)."// &
         " In general, higher values will tame the preconditioner in case of poor initial guesses."// &
         " A negative value will leave the choice to CP2K depending on type of preconditioner.", &
         usage="ENERGY_GAP 0.001", &
         default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EPS_TAYLOR", &
         variants=(/"EPSTAYLOR"/), &
         description="Target accuracy of the taylor expansion for the matrix functions, should normally be kept as is.", &
         usage="EPS_TAYLOR 1.0E-15", &
         default_r_val=1.0E-16_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="MAX_TAYLOR", &
         description="Maximum order of the Taylor expansion before diagonalisation is preferred, for large parallel runs"// &
         " a slightly higher order could sometimes result in a small speedup.", &
         usage="MAX_TAYLOR 5", &
         default_i_val=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ROTATION", &
                          description="Introduce additional variables so that rotations of the occupied"// &
                          " subspace are allowed as well, only needed for cases where the energy is not invariant under"// &
                          " a rotation of the occupied subspace such as non-singlet restricted calculations"// &
                          " or fractional occupations.", &
                          usage="ROTATION", lone_keyword_l_val=.TRUE., &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ENERGIES", &
                          description="Optimize orbital energies for use in Fermi-Dirac smearing "// &
                          "(requires ROTATION and FD smearing to be active).", &
                          usage="ENERGIES", lone_keyword_l_val=.TRUE., &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OCCUPATION_PRECONDITIONER", &
                          description="Preconditioner with the occupation numbers (FD smearing)", &
                          usage="OCCUPATION_PRECONDITIONER", lone_keyword_l_val=.TRUE., &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NONDIAG_ENERGY", &
                          description="Add a non-diagonal energy penalty (FD smearing)", &
                          usage="NONDIAG_ENERGY", lone_keyword_l_val=.TRUE., &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NONDIAG_ENERGY_STRENGTH", &
                          description="The prefactor for the non-diagonal energy penalty (FD smearing)", &
                          usage="NONDIAG_ENERGY_STRENGTH", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_ot_section

! **************************************************************************************************
!> \brief creates the diagonalization section
!> \param section ...
!> \par History
!>      10.2008 created [JGH]
! **************************************************************************************************
   SUBROUTINE create_diagonalization_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="DIAGONALIZATION", &
                          description="Set up type and parameters for Kohn-Sham matrix diagonalization.", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of the diagonalization method", &
                          usage="&DIAGONALIZATION T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ALGORITHM", &
                          description="Algorithm to be used for diagonalization", &
                          usage="ALGORITHM STANDARD", &
                          default_i_val=diag_standard, &
                          enum_c_vals=s2a("STANDARD", "OT", "LANCZOS", "DAVIDSON", "FILTER_MATRIX"), &
                          enum_desc=s2a("Standard diagonalization: LAPACK methods or Jacobi.", &
                                        "Iterative diagonalization using OT method", &
                                        "Block Krylov-space approach to self-consistent diagonalisation", &
                                        "Preconditioned blocked Davidson", &
                                        "Filter matrix diagonalization"), &
                          enum_i_vals=(/diag_standard, diag_ot, diag_block_krylov, diag_block_davidson, &
                                        diag_filter_matrix/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBI_THRESHOLD", &
                          description="Controls the accuracy of the pseudo-diagonalization method using Jacobi rotations", &
                          usage="JACOBI_THRESHOLD 1.0E-6", &
                          default_r_val=1.0E-7_dp, &
                          citations=(/Stewart1982/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_JACOBI", &
                          description="Below this threshold value for the SCF convergence the pseudo-diagonalization "// &
                          "method using Jacobi rotations is activated. This method is much faster than a "// &
                          "real diagonalization and it is even speeding up while achieving full convergence. "// &
                          "However, it needs a pre-converged wavefunction obtained by at least one real "// &
                          "diagonalization which is further optimized while keeping the original eigenvalue "// &
                          "spectrum. The MO eigenvalues are NOT updated. The method might be useful to speed "// &
                          "up calculations for large systems e.g. using a semi-empirical method.", &
                          usage="EPS_JACOBI 1.0E-5", &
                          default_r_val=0.0_dp, &
                          citations=(/Stewart1982/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ADAPT", &
                          description="Required accuracy in iterative diagonalization as compared to current SCF convergence", &
                          usage="EPS_ADAPT 0.01", &
                          default_r_val=0._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximum number of iterations in iterative diagonalization", &
                          usage="MAX_ITER 20", &
                          default_i_val=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ITER", &
                          description="Required accuracy in iterative diagonalization", &
                          usage="EPS_ITER 1.e-8", &
                          default_r_val=1.e-8_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL create_ot_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_krylov_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_diag_subspace_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_davidson_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_filtermatrix_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_diagonalization_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_davidson_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="DAVIDSON", &
                          description=" ", &
                          n_keywords=2, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="PRECONDITIONER", &
         description="Type of preconditioner to be used with all minimization schemes. ", &
         usage="PRECONDITIONER FULL_ALL", &
         default_i_val=ot_precond_full_all, &
         enum_c_vals=s2a("FULL_ALL", "FULL_SINGLE_INVERSE", "NONE"), &
         enum_desc=s2a("Most effective state selective preconditioner based on diagonalization ", &
                       "Based on H-eS cholesky inversion, similar to FULL_SINGLE in preconditioning efficiency "// &
                       "but cheaper to construct, might be somewhat less robust. Recommended for large systems.", &
                       "skip preconditioning"), &
         enum_i_vals=(/ot_precond_full_all, ot_precond_full_single_inverse, ot_precond_none/), &
         citations=(/VandeVondele2003/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PRECOND_SOLVER", &
                          description="How the preconditioner is applied to the residual.", &
                          usage="PRECOND_SOLVER DIRECT", &
                          default_i_val=ot_precond_solver_default, &
                          enum_c_vals=s2a("DEFAULT", "DIRECT", "INVERSE_CHOLESKY"), &
                          enum_desc=s2a("the default", "Cholesky decomposition followed by triangular solve "// &
                                        "(works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)", &
                                        "Cholesky decomposition followed by explicit inversion "// &
                                        "(works for FULL_KINETIC/SINGLE_INVERSE/S_INVERSE)"), &
                          enum_i_vals=(/ot_precond_solver_default, &
                                        ot_precond_solver_direct, &
                                        ot_precond_solver_inv_chol/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="ENERGY_GAP", &
         description="Should be an estimate for the energy gap [a.u.] (HOMO-LUMO) and is used in preconditioning, "// &
         "especially effective with the FULL_ALL preconditioner, in which case it should be an underestimate "// &
         "of the gap (0.001 doing normally fine). For the other preconditioners, making this value larger (0.2)"// &
         " will tame the preconditioner in case of poor initial guesses.", &
         usage="ENERGY_GAP 0.001", &
         default_r_val=0.2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NEW_PREC_EACH", &
                          description="Number of SCF iterations after which a new Preconditioner is computed", &
                          usage="NEW_PREC_EACH 10", default_i_val=20)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FIRST_PREC", &
                          description="First SCF iteration at which a Preconditioner is employed", &
                          usage="FIRST_PREC 1", default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CONV_MOS_PERCENT", &
                          description="Minimal percent of MOS that have to converge within the Davidson loop"// &
                          " before the SCF iteration is completed and a new Hamiltonian is computed", &
                          usage="CONV_MOS_PERCENT 0.8", default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SPARSE_MOS", &
                          description="Use MOS as sparse matrix and avoid as much as possible multiplications with full matrices", &
                          usage="SPARSE_MOS", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_davidson_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_krylov_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="KRYLOV", &
                          description=" ", &
                          n_keywords=2, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NKRYLOV", &
                          description="Dimension of the Krylov space used for the Lanczos refinement", &
                          usage="NKRYLOV 20", &
                          default_i_val=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NBLOCK", &
                          description="Size of the block of vectors refined simultaneously by the Lanczos procedure", &
                          usage="NBLOCK 1", &
                          default_i_val=32)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_KRYLOV", &
                          description="Convergence criterion for the MOs", &
                          usage="EPS_KRYLOV 0.00001", &
                          default_r_val=0.0000001_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_STD_DIAG", &
                          description="Level of convergence to be reached before starting the Lanczos procedure."// &
                          " Above this threshold a standard diagonalization method is used."// &
                          " If negative Lanczos is started at the first iteration", &
                          usage="EPS_STD_DIAG 0.001", &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHECK_MOS_CONV", &
                          description="This requires to check the convergence of MOS also when standard "// &
                          "diagonalization steps are performed, if the block krylov approach is active.", &
                          usage="CHECK_MOS_CONV T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_krylov_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_diag_subspace_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="DIAG_SUB_SCF", &
                          description="Activation of self-consistenf subspace refinement by diagonalization "// &
                          "of H by adjusting the occupation but keeping the MOS unchanged.", &
                          n_keywords=2, n_subsections=1, repeats=.FALSE.)

      NULLIFY (keyword, subsection)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of inner SCF loop to refine occupations in MOS subspace", &
                          usage="&DIAG_SUB_SCF T", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximum number of iterations for the SCF inner loop", &
                          usage="MAX_ITER 20", &
                          default_i_val=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ENE", &
                          description="Required energy accuracy for convergence of subspace diagonalization", &
                          usage="EPS_ENE 1.e-8", &
                          default_r_val=1.e-4_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ADAPT_SCF", &
                          description="Required density matrix accuracy as compared to current SCF convergence", &
                          usage="EPS_ADAPT_SCF 1.e-1", &
                          default_r_val=1._dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EPS_SKIP_SUB_DIAG", &
         description="Level of convergence to be reached before starting the internal loop of subspace rotations."// &
         " Above this threshold only the outer diagonalization method is used."// &
         " If negative the subspace rotation is started at the first iteration", &
         usage="EPS_SKIP_SUB_DIAG 0.001", &
         default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_mixing_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
   END SUBROUTINE create_diag_subspace_section

! **************************************************************************************************
!> \brief      Create CP2K input section for the smearing of occupation numbers
!> \param section ...
!> \date       27.08.2008
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE create_smear_section(section)

      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, &
                          name="SMEAR", &
                          description="Define the smearing of the MO occupation numbers", &
                          n_keywords=6, &
                          n_subsections=0, &
                          repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="_SECTION_PARAMETERS_", &
                          description="Controls the activation of smearing", &
                          usage="&SMEAR ON", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="METHOD", &
                          description="Smearing method to be applied", &
                          usage="METHOD Fermi_Dirac", &
                          default_i_val=smear_energy_window, &
                          enum_c_vals=s2a("FERMI_DIRAC", "ENERGY_WINDOW", "LIST"), &
                          enum_i_vals=(/smear_fermi_dirac, smear_energy_window, smear_list/), &
                          enum_desc=s2a("Fermi-Dirac distribution defined by the keyword ELECTRONIC_TEMPERATURE", &
                                        "Energy window defined by the keyword WINDOW_SIZE", &
                                        "Use a fixed list of occupations"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="LIST", &
                          description="A list of fractional occupations to use. Must match the number of states "// &
                          "and sum up to the correct number of electrons", &
                          repeats=.FALSE., &
                          n_var=-1, &
                          type_of_var=real_t, &
                          usage="LIST 2.0 0.6666 0.6666 0.66666 0.0 0.0")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="ELECTRONIC_TEMPERATURE", &
                          variants=s2a("ELEC_TEMP", "TELEC"), &
                          description="Electronic temperature in the case of Fermi-Dirac smearing", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=cp_unit_to_cp2k(value=300.0_dp, unit_str="K"), &
                          unit_str="K", &
                          usage="ELECTRONIC_TEMPERATURE [K] 300")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="EPS_FERMI_DIRAC", &
                          description="Accuracy checks on occupation numbers use this as a tolerance", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=1.0E-10_dp, &
                          usage="EPS_FERMI_DIRAC 1.0E-6")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="WINDOW_SIZE", &
                          description="Size of the energy window centred at the Fermi level", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=0.0_dp, &
                          unit_str="au_e", &
                          usage="WINDOW_SIZE [eV] 0.3")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FIXED_MAGNETIC_MOMENT", &
                          description="Imposed difference between the numbers of electrons of spin up "// &
                          "and spin down: m = n(up) - n(down). A negative value (default) allows "// &
                          "for a change of the magnetic moment. -1 specifically keeps an integer "// &
                          "number of spin up and spin down electrons.", &
                          repeats=.FALSE., &
                          n_var=1, &
                          type_of_var=real_t, &
                          default_r_val=-100.0_dp, &
                          usage="FIXED_MAGNETIC_MOMENT 1.5")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_smear_section

! **************************************************************************************************
!> \brief Creates the input section for defining CDFT constraints.
!> \param section the section to create
! **************************************************************************************************
   SUBROUTINE create_cdft_control_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: group_section, print_key, subsection

      NULLIFY (keyword, subsection, group_section, print_key)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="CDFT", &
                          description="Parameters needed to set up a constrained DFT calculation."// &
                          " Each repetition of the ATOM_GROUP section defines a new constraint."// &
                          " The constraint(s) is (are) converged in a separate external SCF loop with settings"// &
                          " read from the OUTER_SCF section. Supported constraints: Becke and Gaussian"// &
                          " Hirshfeld (partial).", n_keywords=8, n_subsections=2, &
                          repeats=.FALSE., citations=(/Holmberg2017, Holmberg2018/))

      NULLIFY (subsection, keyword)
      CALL create_outer_scf_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_becke_constraint_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_hirshfeld_constraint_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, __LOCATION__, name="TYPE_OF_CONSTRAINT", &
                          description="Specifies the type of constraint used.", &
                          usage="TYPE_OF_CONSTRAINT (NONE|HIRSHFELD|BECKE)", &
                          enum_c_vals=s2a("NONE", "HIRSHFELD", "BECKE"), &
                          enum_i_vals=(/outer_scf_none, outer_scf_hirshfeld_constraint, &
                                        outer_scf_becke_constraint/), &
                          enum_desc=s2a("No constraint (disables section).", &
                                        "Gaussian Hirshfeld constraint. Partial implementation: no forces. "// &
                                        "Requires corresponding section. Not as extensively tested.", &
                                        "Becke constraint. Requires corresponding section."), &
                          citations=(/Becke1988b/), &
                          default_i_val=outer_scf_none)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STRENGTH", &
                          description="Constraint force constants (Lagrange multipliers). "// &
                          "Give one value per constraint group.", &
                          type_of_var=real_t, n_var=-1, &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TARGET", &
                          description="Constraint target values. Give one value per constraint group. "// &
                          "The target value is the desired number of valence electrons, spin moment, or the number of "// &
                          "alpha or beta electrons on the atoms that define the constraint, suitably multiplied by "// &
                          "atomic coefficients in case a relative constraint between two sets of atoms is employed. "// &
                          "Note that core charges are not subtracted from the target value.", &
                          usage="TARGET {real}", repeats=.FALSE., &
                          type_of_var=real_t, n_var=-1, &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMIC_CHARGES", &
                          description="Calculate atomic CDFT charges with selected weight function"// &
                          " (Z = Z_core - Z_CDFT). With fragment based constraints, charges are"// &
                          " relative to the fragment reference state i.e. Z = Z_CDFT -"// &
                          " Z_frag_reference. Note: if the number of atoms is greater than the"// &
                          " default pw_pool max cache, calculation of atomic CDFT charges"// &
                          " will prompt harmless warnings during deallocation of atomic grids.", &
                          usage="ATOMIC_CHARGES", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FRAGMENT_A_FILE_NAME", variants=(/"FRAGMENT_A_FILE"/), &
                          description="Name of the reference total electron density cube file for fragment A."// &
                          " May include a path. The reference electron density needs to be outputted"// &
                          " on the same grid as the full system (same cutoff and cell, output stride 1).", &
                          usage="FRAGMENT_A_FILE_NAME <FILENAME>", &
                          default_lc_val="fragment_a.cube")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FRAGMENT_B_FILE_NAME", variants=(/"FRAGMENT_B_FILE"/), &
                          description="Name of the reference total electron density cube file for fragment B."// &
                          " May include a path. The reference electron density needs to be outputted"// &
                          " on the same grid as the full system (same cutoff and cell, output stride 1).", &
                          usage="FRAGMENT_B_FILE_NAME <FILENAME>", &
                          default_lc_val="fragment_b.cube")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FRAGMENT_A_SPIN_FILE", &
                          variants=(/"FRAGMENT_A_SPIN_FILE_NAME"/), &
                          description="Name of the reference spin density cube file for fragment A."// &
                          " May include a path. The reference spin density needs to be outputted"// &
                          " on the same grid as the full system (same cutoff and cell, output stride 1).", &
                          usage="FRAGMENT_A_FILE_NAME <FILENAME>", &
                          default_lc_val="fragment_a_spin.cube")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FRAGMENT_B_SPIN_FILE", &
                          variants=(/"FRAGMENT_B_SPIN_FILE_NAME"/), &
                          description="Name of the reference spin density cube file for fragment B."// &
                          " May include a path. The reference spin density needs to be outputted"// &
                          " on the same grid as the full system (same cutoff and cell, output stride 1).", &
                          usage="FRAGMENT_B_FILE_NAME <FILENAME>", &
                          default_lc_val="fragment_b_spin.cube")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FLIP_FRAGMENT_A", &
                          description="Logical which determines if the reference spin difference density "// &
                          "(rho_alpha-rho_beta) for fragment A should be flipped. With default (off) "// &
                          "value, the fragment is constrained to have more alpha than beta electrons "// &
                          "if the isolated fragment has unpaired electrons. Useful in conjunction with "// &
                          "FLIP_FRAGMENT_B.", &
                          usage="FLIP_FRAGMENT_A", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FLIP_FRAGMENT_B", &
                          description="Logical which determines if the reference spin difference density "// &
                          "(rho_alpha-rho_beta) for fragment B should be flipped. With default (off) "// &
                          "value, the fragment is constrained to have more alpha than beta electrons "// &
                          "if the isolated fragment has unpaired electrons. Useful in conjunction with "// &
                          "FLIP_FRAGMENT_A.", &
                          usage="FLIP_FRAGMENT_B", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "PROGRAM_RUN_INFO", &
                                       description="Controls the printing of basic info about the method.", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")

      CALL section_create(subsection, __LOCATION__, name="WEIGHT_FUNCTION", &
                          description="Controls the printing of cube files with "// &
                          "the CDFT weight function(s). Intended for single-point testing. "// &
                          "In multistep simulations, generated cube files are overwritten each step.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(print_key, subsection)
      CALL section_release(subsection)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL section_create(group_section, __LOCATION__, name="ATOM_GROUP", &
                          description="Define a group of atoms for use in a CDFT constraint. Each repetition of "// &
                          "this section creates a new constraint.", &
                          n_keywords=4, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMS", &
                          description="Specifies the list of atoms that are included in the constraint group.", &
                          usage="ATOMS {integer} {integer} .. {integer}", &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(group_section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COEFF", &
                          description="Defines coefficients for the atoms in the list of atoms. Accepts values +/-1.0.", &
                          usage="COEFF 1.0 -1.0", repeats=.TRUE., &
                          type_of_var=real_t, n_var=-1)
      CALL section_add_keyword(group_section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CONSTRAINT_TYPE ", &
                          description="Determines what type of constraint to apply. ", &
                          usage="CONSTRAINT_TYPE (CHARGE|MAGNETIZATION|ALPHA|BETA)", &
                          enum_c_vals=s2a("CHARGE", "MAGNETIZATION", "ALPHA", "BETA"), &
                          enum_i_vals=(/cdft_charge_constraint, cdft_magnetization_constraint, &
                                        cdft_alpha_constraint, cdft_beta_constraint/), &
                          enum_desc=s2a("Total charge density constraint (rho_alpha + rho_beta).", &
                                        "Magnetization density constraint (rho_alpha - rho_beta).", &
                                        "Alpha spin density constraint.", &
                                        "Beta spin density constraint."), &
                          default_i_val=cdft_charge_constraint)
      CALL section_add_keyword(group_section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FRAGMENT_CONSTRAINT", &
                          description="Use a fragment based constraint. "// &
                          "Takes as input the electron densities of two isolated fragments in the "// &
                          "same geometry that they have in the full system. "// &
                          "The isolated fragment densities are read from cube files defined in FRAGMENT_{A,B}_FILE. "// &
                          "For magnetization density constraints, additional files containing the spin difference "// &
                          "densities must be defined with the keywords FRAGMENT_{A,B}_SPIN_FILE. "// &
                          "With this keyword active, the target value of the constraint is calculated from the "// &
                          "the superposition of the isolated fragment densities. Supports only static calculations.", &
                          usage="FRAGMENT_CONSTRAINT", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(group_section, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, group_section)
      CALL section_release(group_section)

      CALL section_create(group_section, __LOCATION__, name="DUMMY_ATOMS", &
                          description="Define an extra group of atoms for which only atomic CDFT charges "// &
                          "should be computed. The section cannot contain any constraint "// &
                          "atoms that were included in section ATOM_GROUP.", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMS", &
                          description="Specifies the list of atoms that are included in the DUMMY_ATOMS group.", &
                          usage="ATOMS {integer} {integer} .. {integer}", &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(group_section, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, group_section)
      CALL section_release(group_section)

      CALL keyword_create(keyword, __LOCATION__, name="REUSE_PRECOND", &
                          description="Reuse a previously built OT preconditioner between subsequent CDFT SCF iterations "// &
                          "if the inner OT SCF loop converged in PRECOND_FREQ steps or less. Intended mainly for MD "// &
                          "simulations with the FULL_ALL preconditioner to speed up the final iterations of the CDFT SCF loop.", &
                          usage="REUSE_PRECOND yes", repeats=.FALSE., n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PRECOND_FREQ", &
                          description="See REUSE_PRECOND.", &
                          usage="PRECOND_FREQ {int}", default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_REUSE", &
                          description="Determines how many times a previously built preconditioner can be reused.", &
                          usage="MAX_REUSE {int}", default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PURGE_HISTORY", &
                          description="Purge wavefunction and constraint history to improve SCF convergence during MD."// &
                          " Counts how often the convergence of the first CDFT SCF iteration takes 2 or more outer SCF"// &
                          " iterations and purges the history if the counter exceeds PURGE_FREQ, and PURGE_OFFSET"// &
                          " MD steps have passed since the last purge."// &
                          " The counter is zeroed after each purge.", &
                          usage="PURGE_HISTORY yes", repeats=.FALSE., n_var=1, &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PURGE_FREQ", &
                          description="See PURGE_HISTORY.", &
                          usage="PURGE_FREQ {int} ", default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PURGE_OFFSET", &
                          description="See PURGE_HISTORY.", &
                          usage="PURGE_OFFSET {int} ", default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COUNTER", &
                          description="A counter to track the total number of energy evaluations. Needed by"// &
                          " some optimizers to print information. Useful mainly for restarts.", &
                          usage="COUNTER {int} ", default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="IN_MEMORY", &
                          description="Precompute gradients due to constraint during"// &
                          " initial formation of constraint and store them in memory. Does"// &
                          " nothing if forces are not calculated.", &
                          usage="IN_MEMORY", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_cdft_control_section

! **************************************************************************************************
!> \brief Creates the input section for defining Gaussian Hirshfeld CDFT constraints.
!> \param section the section to create
! **************************************************************************************************
   SUBROUTINE create_hirshfeld_constraint_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="HIRSHFELD_CONSTRAINT", &
                          description="Parameters for CDFT with a Gaussian Hirshfeld constraint.", &
                          n_keywords=11, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="SHAPE_FUNCTION", &
                          description="Type of shape function used for Hirshfeld partitioning.", &
                          usage="SHAPE_FUNCTION {Gaussian,Density}", repeats=.FALSE., n_var=1, &
                          default_i_val=shape_function_gaussian, &
                          enum_c_vals=s2a("GAUSSIAN", "DENSITY"), &
                          enum_desc=s2a("One Gaussian per atom with radius determined by the keyword GAUSSIAN_SHAPE.", &
                                        "Atomic density expanded in terms of multiple Gaussians."), &
                          enum_i_vals=(/shape_function_gaussian, shape_function_density/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAUSSIAN_SHAPE", &
                          description="Specifies the type of Gaussian used for SHAPE_FUNCTION GAUSSIAN.", &
                          usage="GAUSSIAN_SHAPE (SINGLE|VDW|COVALENT|USER)", &
                          enum_c_vals=s2a("DEFAULT", "SINGLE", "VDW", "COVALENT", "USER"), &
                          enum_i_vals=(/radius_default, radius_single, radius_vdw, radius_covalent, radius_user/), &
                          enum_desc=s2a("Use covalent radii (in angstrom) to construct Gaussians, but fixed"// &
                                        " 1.0_dp radius for elements with a radius larger than this value.", &
                                        "Single Gaussian for all atom types with radius given by GAUSSIAN_RADIUS.", &
                                        "Use van der Waals radii to construct Gaussians.", &
                                        "Use covalent radii to construct Gaussians.", &
                                        "Use user defined radii (keyword ATOMIC_RADII) to construct Gaussians."), &
                          default_i_val=radius_default)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GAUSSIAN_RADIUS", &
                          description="Radius parameter controlling the creation of Gaussians.", &
                          usage="GAUSSIAN_RADIUS <REAL>", &
                          unit_str="angstrom", &
                          default_r_val=cp_unit_to_cp2k(3.0_dp, "angstrom"), &
                          type_of_var=real_t, n_var=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMIC_RADII", &
                          description="Defines custom radii to setup the spherical Gaussians. "// &
                          "Give one value per element in the same order as they "// &
                          "appear in the input coordinates.", &
                          usage="ATOMIC_RADII {real} {real} {real}", repeats=.FALSE., &
                          unit_str="angstrom", &
                          type_of_var=real_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="USE_BOHR", &
                          description="Convert the Gaussian radius from angstrom to bohr. This results in a larger "// &
                          "Gaussian than without unit conversion.", &
                          usage="CAVITY_USE_BOHR TRUE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PRINT_DENSITY", &
                          description="Logical to control printing of Hirshfeld densities to .cube file.", &
                          usage="PRINT_DENSITY TRUE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMS_MEMORY", &
                          description="Number of atomic gradients to store in memory.", &
                          usage="ATOMS_MEMORY", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=80)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="USE_ATOMIC_CUTOFF", &
                          description="Logical to control use of ATOMIC_CUTOFF.", &
                          usage="USE_ATOMIC_CUTOFF TRUE", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_CUTOFF", &
                          description="Numerical cutoff for calculation of weight function.", &
                          usage="EPS_CUTOFF {real} ", default_r_val=1.0e-12_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMIC_CUTOFF", &
                          description="Numerical cutoff for calculation of Hirshfeld densities.", &
                          usage="ATOMIC_CUTOFF {real} ", default_r_val=1.0e-12_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_hirshfeld_constraint_section

! **************************************************************************************************
!> \brief Create input section to define CDFT constraint settings specific to Becke weight function.
!> \param section the section to create
! **************************************************************************************************
   SUBROUTINE create_becke_constraint_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="BECKE_CONSTRAINT", &
                          description="Define settings influencing the construction of the Becke weight function.", &
                          n_keywords=13, repeats=.FALSE., citations=(/Becke1988b/))

      CALL keyword_create(keyword, __LOCATION__, name="ADJUST_SIZE", &
                          description="Adjust Becke cell boundaries with atomic"// &
                          " radii to generate a heteronuclear cutoff profile. These"// &
                          " radii are defined with the keyword ATOMIC_RADII.", &
                          usage="ADJUST_SIZE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMIC_RADII", &
                          description="Defines atomic radii to generate a heteronuclear cutoff profile."// &
                          " Give one value per element in the same order as they"// &
                          " appear in the input coordinates.", &
                          usage="ATOMIC_RADII {real} {real} {real}", repeats=.FALSE., &
                          unit_str="angstrom", &
                          type_of_var=real_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SHOULD_SKIP", &
                          description="If grid point is farther than GLOBAL_CUTOFF from all constraint atoms, "// &
                          "move directly to next grid point, thus saving computational resources.", &
                          usage="SHOULD_SKIP", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CAVITY_CONFINE", &
                          description="Activates Gaussian cavity confinement. The constraint is evaluated only inside "// &
                          "the cavity. The cavity is formed by summing spherical Gaussians centered on the constraint atoms.", &
                          usage="CAVITY_CONFINE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CAVITY_SHAPE", &
                          description="Specifies the type of Gaussian cavity used.", &
                          usage="CAVITY_SHAPE (SINGLE|VDW|COVALENT|USER)", &
                          enum_c_vals=s2a("DEFAULT", "SINGLE", "VDW", "COVALENT", "USER"), &
                          enum_i_vals=(/radius_default, radius_single, radius_vdw, radius_covalent, radius_user/), &
                          enum_desc=s2a("Use covalent radii (in angstrom) to construct Gaussians, but fixed"// &
                                        " 1.0_dp radius for elements with a radius larger than this value.", &
                                        "Single Gaussian for all atom types with radius given by CAVITY_RADIUS.", &
                                        "Use van der Waals radii to construct Gaussians.", &
                                        "Use covalent radii to construct Gaussians.", &
                                        "Use user defined radii (keyword ATOMIC_RADII) to construct Gaussians."), &
                          default_i_val=radius_default)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CAVITY_USE_BOHR", &
                          description="Convert the cavity radius from angstrom to bohr. This results in a larger"// &
                          " confinement cavity than without unit conversion.", &
                          usage="CAVITY_USE_BOHR TRUE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CAVITY_PRINT", &
                          description="Print cavity in Gaussian cube file format. Currently, printing options"// &
                          " are hardcoded.", &
                          usage="CAVITY_PRINT", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CAVITY_RADIUS", &
                          description="Radius parameter controlling the creation of Gaussian cavity confinement.", &
                          usage="CAVITY_RADIUS <REAL>", &
                          unit_str="angstrom", &
                          default_r_val=cp_unit_to_cp2k(3.0_dp, "angstrom"), &
                          type_of_var=real_t, n_var=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_CAVITY", &
                          description="Density threshold for cavity creation. Grid points where the Gaussian"// &
                          " density falls below the threshold are ignored.", &
                          usage="EPS_CAVITY {real} ", default_r_val=1.0e-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CUTOFF_TYPE", &
                          description="Specifies the type of cutoff used when building the Becke weight function.", &
                          usage="CUTOFF_TYPE (GLOBAL|ELEMENT)", &
                          enum_c_vals=s2a("GLOBAL", "ELEMENT"), &
                          enum_i_vals=(/becke_cutoff_global, becke_cutoff_element/), &
                          enum_desc=s2a("Use a single value for all elements. Read from GLOBAL_CUTOFF.", &
                                        "Use a different value for all elements. Values read from ELEMENT_CUTOFF."), &
                          default_i_val=becke_cutoff_global)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GLOBAL_CUTOFF", &
                          description="Parameter used to select which atoms contribute to the"// &
                          " weight function at each real space grid point.", &
                          usage="GLOBAL_CUTOFF <REAL>", &
                          unit_str="angstrom", &
                          default_r_val=cp_unit_to_cp2k(3.1750632515_dp, "angstrom"), &
                          type_of_var=real_t, n_var=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ELEMENT_CUTOFF", &
                          description="Defines element specific cutoffs to decide which atoms contribute to the"// &
                          " weight function at each real space grid point. Give one value per element in the same"// &
                          " order as they appear in the coordinates.", &
                          usage="ELEMENT_CUTOFF {real} {real} {real}", repeats=.FALSE., &
                          unit_str="angstrom", &
                          type_of_var=real_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="IN_MEMORY", &
                          description="Precompute gradients due to Becke constraint during"// &
                          " initial formation of constraint and store them in memory. Useful"// &
                          " in combination with confinement, memory intensive otherwise. Does"// &
                          " nothing if forces are not calculated.", &
                          usage="IN_MEMORY", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_becke_constraint_section

! **************************************************************************************************
!> \brief creates the input section for parameters related to CDFT specific optimizers
!> \param section the section to be created
!> \par History
!>      03.2018 separated from create_outer_scf_section [Nico Holmberg]
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE create_cdft_opt_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="CDFT_OPT", &
                          description="Parameters controlling optimization methods that are compatible "// &
                          "only with CDFT based constraints (i.e. CDFT SCF is active). Specifically, "// &
                          "the control parameters for the Broyden and Newton optimizers are defined in this "// &
                          "section.", &
                          n_keywords=10, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BROYDEN_TYPE", &
                          description="Specifies the Broyden optimizer variant to use.", &
                          usage="BROYDEN_TYPE BT1", &
                          default_i_val=broyden_type_1, &
                          enum_c_vals=s2a("BT1", "BT1_EXPLICIT", "BT2", "BT2_EXPLICIT", &
                                          "BT1_LS", "BT1_EXPLICIT_LS", "BT2_LS", "BT2_EXPLICIT_LS"), &
                          enum_desc=s2a("Broyden's first method, also known as the good method. The initial Jacobian"// &
                                        " is built from MD history if available. Otherwise switches to SD for one"// &
                                        " SCF iteration until a Jacobian can be built from the SCF history.", &
                                        "Same as BT1, but computes the explicit Jacobian with finite differences. "// &
                                        "Requires a CDFT SCF procedure to be active.", &
                                        "Same as BT1, but uses Broyden's second method, also known as the bad method.", &
                                        "Same as BT1_EXPLICIT, but using Broyden's second method.", &
                                        "Same as BT1, but uses backtracking line search for optimizing the step size "// &
                                        "(see optimizer NEWTON_LS).", &
                                        "Same as BT1_EXPLICIT, but uses backtracking line search for optimizing the step size.", &
                                        "Same as BT2, but uses backtracking line search for optimizing the step size.", &
                                        "Same as BT2_EXPLICIT, but uses backtracking line search for optimizing the step size."), &
                          enum_i_vals=(/broyden_type_1, broyden_type_1_explicit, broyden_type_2, &
                                        broyden_type_2_explicit, broyden_type_1_ls, broyden_type_1_explicit_ls, &
                                        broyden_type_2_ls, broyden_type_2_explicit_ls/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBIAN_TYPE", &
                          description="Finite difference method used to calculate the inverse Jacobian "// &
                          "needed by some optimizers. Compatible only with CDFT constraints.", &
                          usage="JACOBIAN_TYPE FD1", &
                          default_i_val=jacobian_fd1, &
                          enum_c_vals=s2a("FD1", "FD1_BACKWARD", "FD2", "FD2_BACKWARD", "FD1_CENTRAL"), &
                          enum_desc=s2a("First order forward difference (one extra energy evaluation per constraint).", &
                                        "First order backward difference (one extra energy evaluation per constraint).", &
                                        "Second order forward difference (two extra energy evaluations per constraint).", &
                                        "Second order backward difference (two extra energy evaluations per constraint).", &
                                        "First order central difference (two extra energy evaluations per constraint)."), &
                          enum_i_vals=(/jacobian_fd1, jacobian_fd1_backward, jacobian_fd2, &
                                        jacobian_fd2_backward, jacobian_fd1_central/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBIAN_STEP", &
                          description="Step size to use in the calculation of the inverse Jacobian with finite differences. "// &
                          "Expects one value for all constraints, or one value per constraint.", &
                          usage="JACOBIAN_STEP 5.0E-3 ", n_var=-1, default_r_val=5.0E-3_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBIAN_FREQ", &
                          description="Defines parameters that control how often the explicit Jacobian is built,"// &
                          " which is needed by some optimizers. Expects two values. The first value"// &
                          " determines how many consecutive CDFT SCF iterations should skip a rebuild,"// &
                          " whereas the latter how many MD steps. The values can be zero (meaning never"// &
                          " rebuild) or positive. Both values cannot be zero.", &
                          usage="JACOBIAN_FREQ 1 1", n_var=2, &
                          default_i_vals=(/1, 1/), type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBIAN_RESTART", &
                          description="Restart the inverse Jacobian using the vector defined with keyword JACOBIAN_VECTOR.", &
                          usage="JACOBIAN_RESTART TRUE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBIAN_VECTOR", &
                          description="Defines the inverse Jacobian matrix. Useful for restarting calculations. "// &
                          "Expects n^2 values where n is the total number of constraints. "// &
                          "The matrix should be given in row major order.", &
                          usage="JACOBIAN_VECTOR 1.0 0.0", n_var=-1, type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_LS", &
                          description="The maximum number of backtracking line search steps to perform.", &
                          usage="MAX_LS 5", default_i_val=5)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FACTOR_LS", &
                          description="Control parameter for backtracking line search. The step size is reduced by "// &
                          "this factor on every line search iteration. Value must be between 0 and 1 (exclusive).", &
                          usage="FACTOR_LS 0.5", default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CONTINUE_LS", &
                          description="Continue backtracking line search until MAX_LS steps are reached or the "// &
                          "norm of the CDFT gradient no longer decreases. Default (false) behavior exits the "// &
                          "line search procedure on the first step that the gradient decreases.", &
                          usage="CONTINUE_LS TRUE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_cdft_opt_section

END MODULE input_cp2k_scf
