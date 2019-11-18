!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
MODULE input_cp2k_loc
   USE bibliography,                    ONLY: Hunt2003
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              debug_print_level,&
                                              high_print_level,&
                                              low_print_level
   USE input_constants,                 ONLY: &
        do_loc_both, do_loc_crazy, do_loc_direct, do_loc_homo, do_loc_jacobi, do_loc_l1_norm_sd, &
        do_loc_lumo, do_loc_max, do_loc_min, do_loc_none, do_loc_scdm, op_loc_berry, op_loc_boys, &
        op_loc_pipek
   USE input_cp2k_mm,                   ONLY: create_dipoles_section
   USE input_cp2k_motion_print,         ONLY: add_format_keyword
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
                                              real_t
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_loc'

   PUBLIC :: create_localize_section, print_wanniers

CONTAINS

! **************************************************************************************************
!> \brief parameters fo the localization of wavefunctions
!> \param section ...
!> \par History
!>      03.2005 created [MI]
! **************************************************************************************************

   SUBROUTINE create_localize_section(section)

      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_localize_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, print_section, subsection

      CPASSERT(.NOT. ASSOCIATED(section))

      NULLIFY (keyword, print_key)
      CALL section_create(section, __LOCATION__, name="LOCALIZE", &
                          description="Use one of the available methods to define the localization "// &
                          " and possibly to optimize it to a minimum or a maximum.", &
                          n_keywords=8, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="controls the activation of the MOS localization procedure", &
                          usage="&LOCALIZE T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximum number of iterations used for localization methods", &
                          usage="MAX_ITER 2000", default_i_val=10000)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="MAX_CRAZY_ANGLE", &
         description="Largest allowed angle for the crazy rotations algorithm (smaller is slower but more stable).", &
         usage="MAX_CRAZY_ANGLE 0.1", unit_str="rad", default_r_val=0.2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CRAZY_SCALE", &
                          description="scale angles", &
                          usage="CRAZY_SCALE 0.9", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CRAZY_USE_DIAG", &
                          description="Use diagonalization (slow) or pade based calculation of matrix exponentials.", &
                          usage="CRAZY_USE_DIAG ", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="USE_HISTORY", &
         description="Generate an improved initial guess based on a history of results, which is useful during MD."// &
         "Will only work if the number of states to be localized remains constant.", &
         usage="USE_HISTORY ", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="EPS_OCCUPATION", &
         description="Tolerance in the occupation number to select only fully occupied orbitals for the rotation", &
         usage="EPS_OCCUPATION 1.E-5", default_r_val=1.0E-8_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OUT_ITER_EACH", &
                          description="Every how many iterations of the localization algorithm"// &
                          "(Jacobi) the tolerance value is printed out", &
                          usage="OUT_ITER_EACH 100", default_i_val=100)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_LOCALIZATION", &
                          description="Tolerance used in the convergence criterion of the localization methods.", &
                          usage="EPS_LOCALIZATION 1.0E-2", default_r_val=1.0E-4_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MIN_OR_MAX", &
                          description="Requires the maximization of the spread of the wfn", &
                          usage="MIN_OR_MAX (SPREADMIN|SPREADMAX)", &
                          enum_c_vals=(/"SPREADMIN", "SPREADMAX"/), &
                          enum_i_vals=(/do_loc_min, do_loc_max/), &
                          default_i_val=do_loc_min)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, name="METHOD", &
         description="Method of optimization if any", &
         usage="METHOD (JACOBI|CRAZY|DIRECT|L1SD|SCDM|NONE)", &
         enum_c_vals=s2a("NONE", "JACOBI", "CRAZY", "L1SD", "DIRECT", "SCDM"), &
         enum_i_vals=(/do_loc_none, &
                       do_loc_jacobi, &
                       do_loc_crazy, &
                       do_loc_l1_norm_sd, &
                       do_loc_direct, do_loc_scdm/), &
         enum_desc=s2a("No localization is applied", &
                       "Using 2 x 2 rotations of the orbitals, slow but robust", &
                       "A new fast method is applied, might be slightly less robust than jacobi, but usually much faster", &
                       "Steepest descent minimization of an approximate l1 norm", &
                       "Using a direct minimisation approacha", "Use QR factorization"), &
         default_i_val=do_loc_jacobi)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBI_FALLBACK", &
                          description="Use Jacobi method in case no convergence was achieved"// &
                          " by using the crazy rotations method.", &
                          usage="JACOBI_FALLBACK", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="JACOBI_REFINEMENT", &
                          description="Use Jacobi method to refine the localisation obtained by SCDM", &
                          usage="JACOBI_REFINEMENT", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTART", &
                          description="Restart the localization from a set of orbitals"// &
                          " read from a localization restart file.", &
                          usage="RESTART", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LOCHOMO_RESTART_FILE_NAME", &
                          description="File name where to read the MOS from"// &
                          "which to restart the localization procedure for occupied states", &
                          usage="LOCHOMO_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LOCLUMO_RESTART_FILE_NAME", &
                          description="File name where to read the MOS from"// &
                          "which to restart the localization procedure for unoccupied states", &
                          usage="LOCLUMO_RESTART_FILE_NAME <FILENAME>", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="OPERATOR", &
                          description="Type of opertator which defines the spread functional", &
                          usage="OPERATOR (BERRY|BOYS|PIPEK)", &
                          enum_c_vals=s2a("BERRY", "BOYS", "PIPEK"), &
                          enum_i_vals=(/op_loc_berry, op_loc_boys, op_loc_pipek/), &
                          default_i_val=op_loc_berry)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LIST", &
                          description="Indexes of the occupied wfn to be localized"// &
                          "This keyword can be repeated several times"// &
                          "(useful if you have to specify many indexes).", &
                          usage="LIST 1 2", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LIST_UNOCCUPIED", &
                          description="Indexes of the unoccupied states to be localized, "// &
                          "up to now only valid in combination with GPW. "// &
                          "This keyword has to be present if unoccupied states should be localized. "// &
                          "This keyword can be repeated several times"// &
                          "(useful if you have to specify many indexes).", &
                          usage="LIST 1 2", &
                          n_var=-1, type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STATES", &
                          description="Which states to localize, LUMO up to now only available in GPW", &
                          usage="STATES (HOMO|LUMO|ALL)", &
                          enum_c_vals=s2a("OCCUPIED", "UNOCCUPIED", "ALL"), &
                          enum_i_vals=(/do_loc_homo, do_loc_lumo, do_loc_both/), &
                          default_i_val=do_loc_homo)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create( &
         keyword, __LOCATION__, &
         name="ENERGY_RANGE", &
         description="Select the orbitals to be localized within the given energy range."// &
         "This type of selection cannot be added on top of the selection through a LIST. It reads to reals that are"// &
         " lower and higher boundaries of the energy range.", &
         usage=" ENERGY_RANGE lower_bound {real}, higher_bound {real}", &
         repeats=.FALSE., &
         n_var=2, default_r_vals=(/0._dp, 0._dp/), unit_str='eV', &
         type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (print_section)
      CALL section_create(print_section, __LOCATION__, name="PRINT", &
                          description="Collects all printing options related to the Wannier centers and "// &
                          "properties computed with Wannier centers.", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)
      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "program_run_info", &
                                       description="Controls the printing basic info about the method", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(print_section, print_key)
      CALL section_release(print_key)
      ! Add printing of wannier infos
      CALL print_wanniers(print_section)
      NULLIFY (subsection)
      ! Total Dipoles with wannier
      CALL create_dipoles_section(subsection, "TOTAL_DIPOLE", debug_print_level + 1)
      CALL section_add_subsection(print_section, subsection)
      CALL section_release(subsection)
      ! Molecular Dipoles with wannier
      CALL create_dipoles_section(subsection, "MOLECULAR_DIPOLES", debug_print_level + 1)
      CALL section_add_subsection(print_section, subsection)
      CALL section_release(subsection)
      ! Molecular Mulipole Moments with wannier
      CALL cp_print_key_section_create(subsection, __LOCATION__, name="MOLECULAR_MOMENTS", &
                                       description="Section controlling the calculation of molecular multipole moments.", &
                                       print_level=debug_print_level + 1, filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, name="ORDER", &
                          description="Maximum order of mulitpoles to be calculated.", &
                          usage=" ORDER {integer}", default_i_val=2, type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      !
      CALL section_add_subsection(print_section, subsection)
      CALL section_release(subsection)
      ! Molecular States with wannier
      CALL create_molecular_states_section(subsection)
      CALL section_add_subsection(print_section, subsection)
      CALL section_release(subsection)
      ! Wannier States with wannier
      CALL create_wannier_states_section(subsection)
      CALL section_add_subsection(print_section, subsection)
      CALL section_release(subsection)
      CALL section_add_subsection(section, print_section)
      CALL section_release(print_section)

   END SUBROUTINE create_localize_section

! **************************************************************************************************
!> \brief Controls the printing of the basic info coming from the LOCALIZE
!>        section
!> \param section ...
!> \author teo
! **************************************************************************************************
   SUBROUTINE print_wanniers(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'print_wanniers', routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(ASSOCIATED(section))
      NULLIFY (print_key, keyword)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "WANNIER_CUBES", &
                                       description="Controls the printing of the wannier functions ", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
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
                          description="Indexes of the states to be printed as cube files"// &
                          "This keyword can be repeated several times"// &
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

      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "WANNIER_CENTERS", &
                                       description="Controls the printing of the wannier functions", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="", &
                                       unit_str="angstrom")

      CALL keyword_create(keyword, __LOCATION__, name="IONS+CENTERS", &
                          description="prints out the wannier centers together with the particles", &
                          usage="IONS+CENTERS", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL add_format_keyword(keyword, print_key, pos=.TRUE., &
                              description="Specifies the format of the output file when IONS+CENTERS is enabled.")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "WANNIER_SPREADS", &
                                       description="Controls the printing of the wannier functions", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      NULLIFY (print_key)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "LOC_RESTART", &
                                       description="Controls the printing of restart file for localized MOS", &
                                       print_level=high_print_level, add_last=add_last_numeric, filename="")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE print_wanniers

! **************************************************************************************************
!> \brief creates the input section for the molecular states
!> \param print_key ...
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_molecular_states_section(print_key)
      TYPE(section_type), POINTER                        :: print_key

      CHARACTER(len=*), PARAMETER :: routineN = 'create_molecular_states_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key2

      CPASSERT(.NOT. ASSOCIATED(print_key))
      NULLIFY (print_key2, keyword)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "MOLECULAR_STATES", &
                                       description="Controls printing of molecular states ", &
                                       print_level=high_print_level, filename=" ", citations=(/Hunt2003/))

      CALL keyword_create( &
         keyword, __LOCATION__, name="CUBE_EVAL_RANGE", &
         description="only write cubes if the eigenvalues of the corresponding molecular states lie in the given interval. "// &
         "Default is all states.", &
         usage="CUBE_EVAL_RANGE -1.0 1.0", unit_str="hartree", n_var=2, type_of_var=real_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MARK_STATES", &
                          description="Can be used to mark given molecular states."// &
                          " Sets a mark to both, occupied and unoccupied states. "// &
                          "Occupied states are counted beginning with HOMO=1, "// &
                          "unoccupied states are counted beginning with LUMO=1, "// &
                          "This is only meaningful in combination with WFN_MIX. "// &
                          "First integer specifies the molecule, second integer specifies the state.", &
                          usage="MARK_STATES integer integer", &
                          n_var=2, default_i_vals=(/-1, -1/), type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL cp_print_key_section_create(print_key2, __LOCATION__, "cubes", &
                                       description="Controls the printing of cube files", &
                                       print_level=high_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key2, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_key, print_key2)
      CALL section_release(print_key2)
   END SUBROUTINE create_molecular_states_section

! **************************************************************************************************
!> \brief ...
!> \param print_key ...
! **************************************************************************************************
   SUBROUTINE create_wannier_states_section(print_key)
      TYPE(section_type), POINTER                        :: print_key

      CHARACTER(len=*), PARAMETER :: routineN = 'create_wannier_states_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key2

      CPASSERT(.NOT. ASSOCIATED(print_key))
      NULLIFY (print_key2, keyword)
      CALL cp_print_key_section_create(print_key, __LOCATION__, "WANNIER_STATES", &
                                       description="Controls printing of molecular states ", &
                                       print_level=high_print_level, filename=" ")

      CALL keyword_create( &
         keyword, __LOCATION__, name="CUBE_EVAL_RANGE", &
         description="only write cubes if the eigenvalues of the corresponding molecular states lie in the given interval. "// &
         "Default is all states.", &
         usage="CUBE_EVAL_RANGE -1.0 1.0", unit_str="hartree", n_var=2, type_of_var=real_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MARK_STATES", &
                          description="Can be used to mark given molecular states."// &
                          " Sets a mark to both, occupied and unoccupied states. "// &
                          "Occupied states are counted beginning with HOMO=1, "// &
                          "unoccupied states are counted beginning with LUMO=1, "// &
                          "This is only meaningful in combination with WFN_MIX. "// &
                          "First integer specifies the molecule, second integer specifies the state.", &
                          usage="MARK_STATES integer integer", &
                          n_var=2, default_i_vals=(/-1, -1/), type_of_var=integer_t, repeats=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL cp_print_key_section_create(print_key2, __LOCATION__, "cubes", &
                                       description="Controls the printing of cube files", &
                                       print_level=high_print_level, filename="")
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key2, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(print_key, print_key2)
      CALL section_release(print_key2)
   END SUBROUTINE create_wannier_states_section

END MODULE
