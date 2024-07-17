!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the dft section of the input
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_opt
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              high_print_level,&
                                              silent_print_level
   USE input_constants,                 ONLY: &
        do_lri_opt_all, do_lri_opt_coeff, do_lri_opt_exps, embed_diff, embed_fa, &
        embed_grid_angstrom, embed_grid_bohr, embed_level_shift, embed_none, embed_quasi_newton, &
        embed_resp, embed_steep_desc, gaussian
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: integer_t,&
                                              lchar_t,&
                                              logical_t
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_opt'

   PUBLIC :: create_optimize_lri_basis_section
   PUBLIC :: create_optimize_embed
   PUBLIC :: create_optimize_dmfet

CONTAINS

! **************************************************************************************************
!> \brief input section for optimization of the auxililary basis for LRIGPW
!> \param section the section to create
!> \author Dorothea Golze [05.2014]
! **************************************************************************************************
   SUBROUTINE create_optimize_lri_basis_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="OPTIMIZE_LRI_BASIS", &
                          description="This section specifies the parameters for optimizing "// &
                          "the lri auxiliary basis sets for LRIGPW. The Powell optimizer is used.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword, subsection)

      CALL keyword_create(keyword, __LOCATION__, name="ACCURACY", &
                          description="Target accuracy for the objective function (RHOEND)", &
                          usage="ACCURACY 5.0E-4", default_r_val=1.0E-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_FUN", &
                          description="Maximum number of function evaluations", &
                          usage="MAX_FUN 200", default_i_val=4000)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STEP_SIZE", &
                          description="Initial step size for search algorithm (RHOBEG)", &
                          usage="STEP_SIZE 1.0E-1", default_r_val=5.0E-2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CONDITION_WEIGHT", &
                          description="This keyword allows to give different weight "// &
                          "factors to the condition number (LOG(cond) is used).", &
                          usage="CONDITION_WEIGHT 1.0E-4", default_r_val=1.0E-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="USE_CONDITION_NUMBER", &
                          description="Determines whether condition number should be part "// &
                          "of optimization or not", &
                          usage="USE_CONDITION_NUMBER", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GEOMETRIC_SEQUENCE", &
                          description="Exponents are assumed to be a geometric sequence. "// &
                          "Only the minimal and maximal exponents of one set are optimized and "// &
                          "the other exponents are obtained by geometric progression.", &
                          usage="GEOMETRIC_SEQUENCE", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DEGREES_OF_FREEDOM", &
                          description="Specifies the degrees of freedom in the basis "// &
                          "optimization.", &
                          usage="DEGREES_OF_FREEDOM ALL", &
                          enum_c_vals=s2a("ALL", "COEFFICIENTS", "EXPONENTS"), &
                          enum_desc=s2a("Set all parameters in the basis to be variable.", &
                                        "Set all coefficients in the basis set to be variable.", &
                                        "Set all exponents in the basis to be variable."), &
                          enum_i_vals=(/do_lri_opt_all, do_lri_opt_coeff, do_lri_opt_exps/), &
                          default_i_val=do_lri_opt_exps)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_constrain_exponents_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_optimize_lri_basis_section

! **************************************************************************************************
!> \brief Input for DFT embedding
!> \param section ...
!> \author Vladimir Rybkin [08.2017]
! **************************************************************************************************
   SUBROUTINE create_optimize_embed(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="OPT_EMBED", &
                          description="This section specifies optional parameters for DFT embedding potential optimization.", &
                          n_keywords=19, n_subsections=4, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REG_LAMBDA", &
                          description="Parameter for Yang's regularization "// &
                          "involving kinetic matrix.", &
                          usage="REG_LAMBDA 0.0001", default_r_val=0.0001_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="N_ITER", &
                          description="Maximum number of iterations "// &
                          "in the optimization procedure.", &
                          usage="N_ITER 75", default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TRUST_RAD", &
                          description="Maximum number of iterations "// &
                          "in the optimization procedure.", &
                          usage="TRUST_RAD 0.5", default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DENS_CONV_MAX", &
                          description="Convergence criterion for "// &
                          "the maximum electron density difference.", &
                          usage="DENS_CONV_MAX 0.01", default_r_val=0.01_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DENS_CONV_INT", &
                          description="Convergence criterion for "// &
                          "the integrated electron density difference.", &
                          usage="DENS_CONV_INT 0.1", default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SPIN_DENS_CONV_MAX", &
                          description="Convergence criterion for "// &
                          "the maximum electron density difference.", &
                          usage="DENS_CONV_MAX 0.01", default_r_val=0.01_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SPIN_DENS_CONV_INT", &
                          description="Convergence criterion for "// &
                          "the integrated electron density difference.", &
                          usage="DENS_CONV_INT 0.1", default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OPTIMIZER", &
                          description="Optimize embedding potential.", &
                          usage="OPTIMIZER  LEVEL_SHIFT", &
                          default_i_val=embed_steep_desc, &
                          enum_c_vals=s2a("STEEPEST_DESCENT", "QUASI_NEWTON", "LEVEL_SHIFT"), &
                          enum_desc=s2a("Steepest descent.", "Quasi-Newton.", "Level shift."), &
                          enum_i_vals=(/embed_steep_desc, embed_quasi_newton, embed_level_shift/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GRID_OPT", &
                          description="Optimize embedding potential on the grid. ", &
                          usage="GRID_OPT .TRUE.", &
                          default_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LEEUWEN-BAERENDS", &
                          description="Van Leeuwen-Baerends iterative update. Alternative to Wu-Yang "// &
                          "optimizer. Use only with ADD_CONTST_POT.", &
                          usage="LEEUWEN-BAERENDS .TRUE.", &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FAB", &
                          description="Finzel-Ayers-Bultinck iterative update. Generally, not reliable. ", &
                          usage="FAB .TRUE.", &
                          default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VW_CUTOFF", &
                          description="Cutoff for von Weizsacker potential in "// &
                          "the FAB optimization procedure.", &
                          usage="VW_CUTOFF 0.01", default_r_val=0.01_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VW_SMOOTH_CUT_RANGE", &
                          description="Smooth cutoff range for von Weizsacker potential in "// &
                          "the FAB optimization procedure.", &
                          usage="VW_CUTOFF 1.0", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="POT_GUESS", &
                          description="Specifies the guess of the embedding  "// &
                          "potential. For optimization in finite basis (not grid optimization) "// &
                          "in is a constant part to be added to the one in finite basis. ", &
                          usage="DEGREES_OF_FREEDOM ALL", &
                          enum_c_vals=s2a("NONE", "DIFF", "Fermi_Amaldi", "RESP"), &
                          enum_desc=s2a("Initial guess is zero grid.", &
                                        "Initial density difference. A euristic but working approach.", &
                                        "Fermi-Amaldi potential. More rigorous than DIFF, although less efficient.", &
                                        "Coulomb interaction between the subsystem using RESP charges)"// &
                                        " on the total system."), &
                          enum_i_vals=(/embed_none, embed_diff, embed_fa, embed_resp/), &
                          default_i_val=embed_none)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHARGE_DISTR_WIDTH", &
                          description="Width of the Gaussian representing "// &
                          "point charges. To be used with ADD_COULOMB_POT.", &
                          usage="CHARGE_DISTR_WIDTH 3.000", default_r_val=1.12490_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="READ_EMBED_POT", &
                          description="Read the embedding potential  "// &
                          "restart vector as a guess.", &
                          usage="READ_EMBED_POT .FALSE.", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="READ_EMBED_POT_CUBE", &
                          description="Read the embedding potential  "// &
                          "(restart) from the cube file.", &
                          usage="READ_EMBED_POT_CUBE .FALSE.", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EMBED_RESTART_FILE_NAME", &
                          description="Root of the file name where to read the embedding "// &
                          "potential guess.", &
                          usage="EMBED_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EMBED_CUBE_FILE_NAME", &
                          description="Root of the file name where to read the embedding "// &
                          "potential (guess) as a cube.", &
                          usage="EMBED_CUBE_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EMBED_SPIN_CUBE_FILE_NAME", &
                          description="Root of the file name where to read the spin part "// &
                          "of the embedding potential (guess) as a cube.", &
                          usage="EMBED_SPIN_CUBE_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_print_embed_diff(section)

      CALL create_print_embed_pot_cube(section)

      CALL create_print_embed_restart_vec(section)

      CALL create_print_simple_grid(section)

   END SUBROUTINE create_optimize_embed

! **************************************************************************************************
!> \brief Input for density matrix functional embedding, DMFET
!> \param section ...
!> \author Vladimir Rybkin [08.2018]
! **************************************************************************************************
   SUBROUTINE create_optimize_dmfet(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="OPT_DMFET", &
                          description="This section specifies optional parameters for DMFET matrix potential optimization.", &
                          n_keywords=8, n_subsections=4, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="N_ITER", &
                          description="Maximum number of iterations "// &
                          "in the optimization procedure.", &
                          usage="N_ITER 75", default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TRUST_RAD", &
                          description="Step length "// &
                          "in the optimization procedure.", &
                          usage="TRUST_RAD 0.5", default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DM_CONV_MAX", &
                          description="Convergence criterion for "// &
                          "the maximum element of density matrix difference.", &
                          usage="DM_CONV_MAX 0.01", default_r_val=0.01_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DM_CONV_INT", &
                          description="Convergence criterion for "// &
                          "the total density matrix difference.", &
                          usage="DM_CONV_INT 0.1", default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BETA_DM_CONV_MAX", &
                          description="Convergence criterion for "// &
                          "the maximum element of the beta-spin density "// &
                          "matrix difference.", &
                          usage="DM_CONV_MAX 0.01", default_r_val=0.01_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BETA_DM_CONV_INT", &
                          description="Convergence criterion for "// &
                          "the total beta-spin density matrix difference.", &
                          usage="DM_CONV_INT 0.1", default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="READ_DMFET_POT", &
                          description="Read the matrix embedding potential "// &
                          "(restart) from the cube file.", &
                          usage="READ_DMFET_POT .FALSE.", default_l_val=.FALSE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DMFET_RESTART_FILE_NAME", &
                          description="Root of the file name where to read the matrix "// &
                          "potential guess.", &
                          usage="DMFET_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_optimize_dmfet

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_print_embed_diff(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      NULLIFY (print_key, keyword)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "EMBED_DENS_DIFF", &
                                       description="Controls the printing of cube files with "// &
                                       "embedding densisty differences", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_print_embed_diff

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_print_embed_pot_cube(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      NULLIFY (print_key, keyword)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "EMBED_POT_CUBE", &
                                       description="Controls the printing of cube files with "// &
                                       "with embedding potential", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 1 1 1", n_var=-1, default_i_vals=(/1, 1, 1/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_print_embed_pot_cube

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_print_embed_restart_vec(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(section_type), POINTER                        :: print_key

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "EMBED_POT_VECTOR", &
                                       description="Controls the printing of cube files with "// &
                                       "with embedding potential", &
                                       print_level=silent_print_level, add_last=add_last_numeric, filename="")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_print_embed_restart_vec

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_print_simple_grid(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      NULLIFY (print_key, keyword)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "WRITE_SIMPLE_GRID", &
                                       description="Controls the printing of simple grid "// &
                                       "files with embedding potential: X Y Z value", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")

      CALL keyword_create(keyword, __LOCATION__, name="STRIDE", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 1 1 1", n_var=-1, default_i_vals=(/1, 1, 1/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="UNITS", &
                          description="Units of the volumetric file: Angstrom or Bohr.", &
                          usage="UNITS BOHR", &
                          default_i_val=embed_grid_bohr, &
                          enum_c_vals=s2a("BOHR", "ANGSTROM"), &
                          enum_desc=s2a("Atomic units: Bohr", "Metric units: Angstrom."), &
                          enum_i_vals=(/embed_grid_bohr, embed_grid_angstrom/))
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FOLD_COORD", &
                          description="Activates printing folded coordinates corresponding "// &
                          "to the simple grid. Used as input for external programs.", &
                          usage="FOLD_COORD .TRUE.", n_var=1, type_of_var=logical_t, &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)

      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_print_simple_grid

! **************************************************************************************************
!> \brief input section for constraints for auxiliary basis set optimization
!> \param section the section to create
!> \author Dorothea Golze [11.2014]
! **************************************************************************************************
   SUBROUTINE create_constrain_exponents_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CALL section_create(section, __LOCATION__, name="CONSTRAIN_EXPONENTS", &
                          description="specifies constraints for the exponents of the "// &
                          "lri auxiliary basis sets in the optimization.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SCALE", &
                          description="Defines the upper and lower boundaries as "// &
                          "(1+scale)*exp and (1-scale)*exp. Fermi-like constraint "// &
                          "function", &
                          usage="SCALE 0.3", default_r_val=0.3_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FERMI_EXP", &
                          description="Exponent in the fermi-like constraint function. ", &
                          usage="FERMI_EXP 2.63", default_r_val=2.63391_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_constrain_exponents_section

END MODULE input_cp2k_opt
