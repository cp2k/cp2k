!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility subroutines for mixed CDFT calculations
!> \par   History
!>                 separated from mixed_cdft_methods [01.2017]
!> \author Nico Holmberg [01.2017]
! **************************************************************************************************
MODULE mixed_cdft_utils
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_array_utils,                  ONLY: cp_1d_i_p_type,&
                                              cp_1d_r_p_type,&
                                              cp_2d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_retain,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr_bc
   USE cp_files,                        ONLY: open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_to_fm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_create,&
                                              cp_logger_set,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_init,          ONLY: init_input_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cube_utils,                      ONLY: init_cube_info,&
                                              return_cube_max_iradius
   USE d3_poly,                         ONLY: init_d3_poly_module
   USE dbcsr_api,                       ONLY: dbcsr_desymmetrize,&
                                              dbcsr_get_info,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_release_p,&
                                              dbcsr_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type,&
                                              multiple_fe_list
   USE gaussian_gridlevels,             ONLY: init_gaussian_gridlevel
   USE global_types,                    ONLY: global_environment_type
   USE hirshfeld_types,                 ONLY: create_hirshfeld_type,&
                                              release_hirshfeld_type,&
                                              set_hirshfeld_info
   USE input_constants,                 ONLY: becke_cutoff_element,&
                                              mixed_cdft_parallel,&
                                              mixed_cdft_parallel_nobuild,&
                                              mixed_cdft_serial,&
                                              outer_scf_becke_constraint,&
                                              outer_scf_hirshfeld_constraint,&
                                              shape_function_gaussian
   USE input_section_types,             ONLY: section_vals_duplicate,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_release,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE message_passing,                 ONLY: mp_irecv,&
                                              mp_isend,&
                                              mp_sum,&
                                              mp_wait,&
                                              mp_waitall
   USE mixed_cdft_types,                ONLY: mixed_cdft_result_type_release,&
                                              mixed_cdft_result_type_set,&
                                              mixed_cdft_settings_type,&
                                              mixed_cdft_type,&
                                              mixed_cdft_work_type_init
   USE mixed_environment_types,         ONLY: get_mixed_env,&
                                              mixed_environment_type
   USE pw_env_methods,                  ONLY: pw_env_create
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: HALFSPACE,&
                                              pw_grid_type
   USE pw_grids,                        ONLY: do_pw_grid_blocked_false,&
                                              pw_grid_create,&
                                              pw_grid_release,&
                                              pw_grid_setup
   USE pw_pool_types,                   ONLY: pw_pool_create,&
                                              pw_pool_p_type,&
                                              pw_pool_type
   USE qs_cdft_types,                   ONLY: cdft_control_create,&
                                              cdft_control_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: create_qs_kind_set,&
                                              qs_kind_type
   USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                              realspace_grid_input_type,&
                                              realspace_grid_p_type,&
                                              rs_grid_create,&
                                              rs_grid_create_descriptor,&
                                              rs_grid_print
#include "base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

! *** Public subroutines ***

   PUBLIC :: mixed_cdft_parse_settings, mixed_cdft_transfer_settings, &
             mixed_cdft_init_structures, mixed_cdft_redistribute_arrays, &
             mixed_cdft_print_couplings, map_permutation_to_states, hfun_zero, &
             mixed_cdft_release_work, mixed_cdft_read_block_diag, &
             mixed_cdft_get_blocks, mixed_cdft_diagonalize_blocks, &
             mixed_cdft_assemble_block_diag

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mixed_cdft_utils'

CONTAINS

! **************************************************************************************************
!> \brief Parse settings for mixed cdft calculation and check their consistency
!> \param force_env the force_env that holds the CDFT mixed_env
!> \param mixed_env the mixed_env that holds the CDFT states
!> \param mixed_cdft control section for mixed CDFT
!> \param settings container for settings related to the mixed CDFT calculation
!> \param natom the total number of atoms
!> \par History
!>       01.2017  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_parse_settings(force_env, mixed_env, mixed_cdft, &
                                        settings, natom)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(mixed_environment_type), POINTER              :: mixed_env
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(mixed_cdft_settings_type)                     :: settings
      INTEGER                                            :: natom

      INTEGER                                            :: i, iatom, iforce_eval, igroup, &
                                                            nforce_eval, nkinds
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: constraint_type
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: array_sizes
      LOGICAL                                            :: is_match
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(cp_1d_i_p_type), ALLOCATABLE, DIMENSION(:, :) :: atoms
      TYPE(cp_1d_r_p_type), ALLOCATABLE, DIMENSION(:, :) :: coeff
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(force_env_type), POINTER                      :: force_env_qs
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_environment_type), POINTER                 :: qs_env

      NULLIFY (dft_control, qs_env, pw_env, auxbas_pw_pool, force_env_qs, &
               cdft_control)
      ! Allocate storage for temporaries used for checking settings consistency
      settings%max_nkinds = 30
      nforce_eval = SIZE(force_env%sub_force_env)
      ALLOCATE (settings%grid_span(nforce_eval))
      ALLOCATE (settings%npts(3, nforce_eval))
      ALLOCATE (settings%cutoff(nforce_eval))
      ALLOCATE (settings%rel_cutoff(nforce_eval))
      ALLOCATE (settings%spherical(nforce_eval))
      ALLOCATE (settings%rs_dims(2, nforce_eval))
      ALLOCATE (settings%odd(nforce_eval))
      ALLOCATE (settings%atoms(natom, nforce_eval))
      IF (mixed_cdft%run_type == mixed_cdft_parallel) THEN
         ALLOCATE (settings%coeffs(natom, nforce_eval))
         settings%coeffs = 0.0_dp
      END IF
      ! Some of the checked settings are only defined for certain types of constraints
      ! We nonetheless use arrays that are large enough to contain settings for all constraints
      ! This is not completely optimal...
      ALLOCATE (settings%si(6, nforce_eval))
      ALLOCATE (settings%sb(8, nforce_eval))
      ALLOCATE (settings%sr(5, nforce_eval))
      ALLOCATE (settings%cutoffs(settings%max_nkinds, nforce_eval))
      ALLOCATE (settings%radii(settings%max_nkinds, nforce_eval))
      settings%grid_span = 0
      settings%npts = 0
      settings%cutoff = 0.0_dp
      settings%rel_cutoff = 0.0_dp
      settings%spherical = 0
      settings%is_spherical = .FALSE.
      settings%rs_dims = 0
      settings%odd = 0
      settings%is_odd = .FALSE.
      settings%atoms = 0
      settings%si = 0
      settings%sr = 0.0_dp
      settings%sb = .FALSE.
      settings%cutoffs = 0.0_dp
      settings%radii = 0.0_dp
      ! Get information from the sub_force_envs
      DO iforce_eval = 1, nforce_eval
         IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
         force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
         IF (mixed_env%do_mixed_qmmm_cdft) THEN
            qs_env => force_env_qs%qmmm_env%qs_env
         ELSE
            CALL force_env_get(force_env_qs, qs_env=qs_env)
         END IF
         CALL get_qs_env(qs_env, pw_env=pw_env, dft_control=dft_control)
         IF (.NOT. dft_control%qs_control%cdft) &
            CALL cp_abort(__LOCATION__, &
                          "A mixed CDFT simulation with multiple force_evals was requested, "// &
                          "but CDFT constraints were not active in the QS section of all force_evals!")
         cdft_control => dft_control%qs_control%cdft_control
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         settings%bo = auxbas_pw_pool%pw_grid%bounds_local
         ! Only the rank 0 process collects info about pw_grid and CDFT
         IF (force_env_qs%para_env%ionode) THEN
            ! Grid settings
            settings%grid_span(iforce_eval) = auxbas_pw_pool%pw_grid%grid_span
            settings%npts(:, iforce_eval) = auxbas_pw_pool%pw_grid%npts
            settings%cutoff(iforce_eval) = auxbas_pw_pool%pw_grid%cutoff
            settings%rel_cutoff(iforce_eval) = dft_control%qs_control%relative_cutoff
            IF (auxbas_pw_pool%pw_grid%spherical) settings%spherical(iforce_eval) = 1
            settings%rs_dims(:, iforce_eval) = auxbas_pw_pool%pw_grid%para%rs_dims
            IF (auxbas_pw_pool%pw_grid%grid_span == HALFSPACE) settings%odd(iforce_eval) = 1
            ! Becke constraint atoms/coeffs
            IF (cdft_control%natoms .GT. SIZE(settings%atoms, 1)) &
               CALL cp_abort(__LOCATION__, &
                             "More CDFT constraint atoms than defined in mixed section. "// &
                             "Use default values for MIXED\MAPPING.")
            settings%atoms(1:cdft_control%natoms, iforce_eval) = cdft_control%atoms
            IF (mixed_cdft%run_type == mixed_cdft_parallel) &
               settings%coeffs(1:cdft_control%natoms, iforce_eval) = cdft_control%group(1)%coeff
            ! Integer type settings
            IF (cdft_control%type == outer_scf_becke_constraint) THEN
               settings%si(1, iforce_eval) = cdft_control%becke_control%cutoff_type
               settings%si(2, iforce_eval) = cdft_control%becke_control%cavity_shape
            END IF
            settings%si(3, iforce_eval) = dft_control%multiplicity
            settings%si(4, iforce_eval) = SIZE(cdft_control%group)
            settings%si(5, iforce_eval) = cdft_control%type
            IF (cdft_control%type == outer_scf_hirshfeld_constraint) THEN
               settings%si(6, iforce_eval) = cdft_control%hirshfeld_control%shape_function
               settings%si(6, iforce_eval) = cdft_control%hirshfeld_control%gaussian_shape
            END IF
            ! Logicals
            IF (cdft_control%type == outer_scf_becke_constraint) THEN
               settings%sb(1, iforce_eval) = cdft_control%becke_control%cavity_confine
               settings%sb(2, iforce_eval) = cdft_control%becke_control%should_skip
               settings%sb(3, iforce_eval) = cdft_control%becke_control%print_cavity
               settings%sb(4, iforce_eval) = cdft_control%becke_control%in_memory
               settings%sb(5, iforce_eval) = cdft_control%becke_control%adjust
               settings%sb(8, iforce_eval) = cdft_control%becke_control%use_bohr
            END IF
            IF (cdft_control%type == outer_scf_hirshfeld_constraint) THEN
               settings%sb(8, iforce_eval) = cdft_control%hirshfeld_control%use_bohr
            END IF
            settings%sb(6, iforce_eval) = cdft_control%atomic_charges
            settings%sb(7, iforce_eval) = qs_env%has_unit_metric
            ! Reals
            IF (cdft_control%type == outer_scf_becke_constraint) THEN
               settings%sr(1, iforce_eval) = cdft_control%becke_control%rcavity
               settings%sr(2, iforce_eval) = cdft_control%becke_control%rglobal
               settings%sr(3, iforce_eval) = cdft_control%becke_control%eps_cavity
            END IF
            IF (cdft_control%type == outer_scf_hirshfeld_constraint) THEN
               settings%sr(2, iforce_eval) = cdft_control%hirshfeld_control%radius
            END IF
            settings%sr(4, iforce_eval) = dft_control%qs_control%eps_rho_rspace
            settings%sr(5, iforce_eval) = pw_env%cube_info(pw_env%auxbas_grid)%max_rad_ga
            IF (cdft_control%type == outer_scf_becke_constraint) THEN
               IF (cdft_control%becke_control%cutoff_type == becke_cutoff_element) THEN
                  nkinds = SIZE(cdft_control%becke_control%cutoffs_tmp)
                  IF (nkinds .GT. settings%max_nkinds) &
                     CALL cp_abort(__LOCATION__, &
                                   "More than "//TRIM(cp_to_string(settings%max_nkinds))// &
                                   " unique elements were defined in BECKE_CONSTRAINT\ELEMENT_CUTOFF. Are you sure"// &
                                   " your input is correct? If yes, please increase max_nkinds and recompile.")
                  settings%cutoffs(1:nkinds, iforce_eval) = cdft_control%becke_control%cutoffs_tmp(:)
               END IF
               IF (cdft_control%becke_control%adjust) THEN
                  CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
                  IF (.NOT. SIZE(atomic_kind_set) == SIZE(cdft_control%becke_control%radii_tmp)) &
                     CALL cp_abort(__LOCATION__, &
                                   "Length of keyword BECKE_CONSTRAINT\ATOMIC_RADII does not "// &
                                   "match number of atomic kinds in the input coordinate file.")
                  nkinds = SIZE(cdft_control%becke_control%radii_tmp)
                  IF (nkinds .GT. settings%max_nkinds) &
                     CALL cp_abort(__LOCATION__, &
                                   "More than "//TRIM(cp_to_string(settings%max_nkinds))// &
                                   " unique elements were defined in BECKE_CONSTRAINT\ATOMIC_RADII. Are you sure"// &
                                   " your input is correct? If yes, please increase max_nkinds and recompile.")
                  settings%radii(1:nkinds, iforce_eval) = cdft_control%becke_control%radii_tmp(:)
               END IF
            END IF
            IF (cdft_control%type == outer_scf_hirshfeld_constraint) THEN
               IF (ASSOCIATED(cdft_control%hirshfeld_control%radii)) THEN
                  CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
                  IF (.NOT. SIZE(atomic_kind_set) == SIZE(cdft_control%hirshfeld_control%radii)) &
                     CALL cp_abort(__LOCATION__, &
                                   "Length of keyword HIRSHFELD_CONSTRAINT&RADII does not "// &
                                   "match number of atomic kinds in the input coordinate file.")
                  nkinds = SIZE(cdft_control%hirshfeld_control%radii)
                  IF (nkinds .GT. settings%max_nkinds) &
                     CALL cp_abort(__LOCATION__, &
                                   "More than "//TRIM(cp_to_string(settings%max_nkinds))// &
                                   " unique elements were defined in HIRSHFELD_CONSTRAINT&RADII. Are you sure"// &
                                   " your input is correct? If yes, please increase max_nkinds and recompile.")
                  settings%radii(1:nkinds, iforce_eval) = cdft_control%hirshfeld_control%radii(:)
               END IF
            END IF
         END IF
      END DO
      ! Make sure the grids are consistent
      CALL mp_sum(settings%grid_span, force_env%para_env%group)
      CALL mp_sum(settings%npts, force_env%para_env%group)
      CALL mp_sum(settings%cutoff, force_env%para_env%group)
      CALL mp_sum(settings%rel_cutoff, force_env%para_env%group)
      CALL mp_sum(settings%spherical, force_env%para_env%group)
      CALL mp_sum(settings%rs_dims, force_env%para_env%group)
      CALL mp_sum(settings%odd, force_env%para_env%group)
      is_match = .TRUE.
      DO iforce_eval = 2, nforce_eval
         is_match = is_match .AND. (settings%grid_span(1) == settings%grid_span(iforce_eval))
         is_match = is_match .AND. (settings%npts(1, 1) == settings%npts(1, iforce_eval))
         is_match = is_match .AND. (settings%cutoff(1) == settings%cutoff(iforce_eval))
         is_match = is_match .AND. (settings%rel_cutoff(1) == settings%rel_cutoff(iforce_eval))
         is_match = is_match .AND. (settings%spherical(1) == settings%spherical(iforce_eval))
         is_match = is_match .AND. (settings%rs_dims(1, 1) == settings%rs_dims(1, iforce_eval))
         is_match = is_match .AND. (settings%rs_dims(2, 1) == settings%rs_dims(2, iforce_eval))
         is_match = is_match .AND. (settings%odd(1) == settings%odd(iforce_eval))
      END DO
      IF (.NOT. is_match) &
         CALL cp_abort(__LOCATION__, &
                       "Mismatch detected in the &MGRID settings of the CDFT force_evals.")
      IF (settings%spherical(1) == 1) settings%is_spherical = .TRUE.
      IF (settings%odd(1) == 1) settings%is_odd = .TRUE.
      ! Make sure CDFT settings are consistent
      CALL mp_sum(settings%atoms, force_env%para_env%group)
      IF (mixed_cdft%run_type == mixed_cdft_parallel) &
         CALL mp_sum(settings%coeffs, force_env%para_env%group)
      settings%ncdft = 0
      DO i = 1, SIZE(settings%atoms, 1)
         DO iforce_eval = 2, nforce_eval
            IF (mixed_cdft%run_type == mixed_cdft_parallel) THEN
               IF (settings%atoms(i, 1) /= settings%atoms(i, iforce_eval)) is_match = .FALSE.
               IF (settings%coeffs(i, 1) /= settings%coeffs(i, iforce_eval)) is_match = .FALSE.
            END IF
         END DO
         IF (settings%atoms(i, 1) /= 0) settings%ncdft = settings%ncdft + 1
      END DO
      IF (.NOT. is_match .AND. mixed_cdft%run_type == mixed_cdft_parallel) &
         CALL cp_abort(__LOCATION__, &
                       "Mismatch detected in the &CDFT section of the CDFT force_evals. "// &
                       "Parallel mode mixed CDFT requires identical constraint definitions in both CDFT states. "// &
                       "Switch to serial mode or disable keyword PARALLEL_BUILD if you "// &
                       "want to use nonidentical constraint definitions.")
      CALL mp_sum(settings%si, force_env%para_env%group)
      CALL mp_sum(settings%sr, force_env%para_env%group)
      DO i = 1, SIZE(settings%sb, 1)
         CALL mp_sum(settings%sb(i, 1), force_env%para_env%group)
         DO iforce_eval = 2, nforce_eval
            CALL mp_sum(settings%sb(i, iforce_eval), force_env%para_env%group)
            IF (settings%sb(i, 1) .NEQV. settings%sb(i, iforce_eval)) is_match = .FALSE.
         END DO
      END DO
      DO i = 1, SIZE(settings%si, 1)
         DO iforce_eval = 2, nforce_eval
            IF (settings%si(i, 1) /= settings%si(i, iforce_eval)) is_match = .FALSE.
         END DO
      END DO
      DO i = 1, SIZE(settings%sr, 1)
         DO iforce_eval = 2, nforce_eval
            IF (settings%sr(i, 1) /= settings%sr(i, iforce_eval)) is_match = .FALSE.
         END DO
      END DO
      IF (.NOT. is_match) &
         CALL cp_abort(__LOCATION__, &
                       "Mismatch detected in the &CDFT settings of the CDFT force_evals.")
      ! Some CDFT features are currently disabled for mixed calculations: check that these features were not requested
      IF (mixed_cdft%dlb .AND. .NOT. settings%sb(1, 1)) &
         CALL cp_abort(__LOCATION__, &
                       "Parallel mode mixed CDFT load balancing requires Gaussian cavity confinement.")
      ! Check for identical constraints in case of run type serial/parallel_nobuild
      IF (mixed_cdft%run_type /= mixed_cdft_parallel) THEN
         ! Get array sizes
         ALLOCATE (array_sizes(nforce_eval, settings%si(4, 1), 2))
         array_sizes(:, :, :) = 0
         DO iforce_eval = 1, nforce_eval
            IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
            force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
            IF (mixed_env%do_mixed_qmmm_cdft) THEN
               qs_env => force_env_qs%qmmm_env%qs_env
            ELSE
               CALL force_env_get(force_env_qs, qs_env=qs_env)
            END IF
            CALL get_qs_env(qs_env, dft_control=dft_control)
            cdft_control => dft_control%qs_control%cdft_control
            IF (force_env_qs%para_env%ionode) THEN
               DO igroup = 1, SIZE(cdft_control%group)
                  array_sizes(iforce_eval, igroup, 1) = SIZE(cdft_control%group(igroup)%atoms)
                  array_sizes(iforce_eval, igroup, 2) = SIZE(cdft_control%group(igroup)%coeff)
               END DO
            END IF
         END DO
         ! Sum up array sizes and check consistency
         CALL mp_sum(array_sizes, force_env%para_env%group)
         IF (ANY(array_sizes(:, :, 1) /= array_sizes(1, 1, 1)) .OR. &
             ANY(array_sizes(:, :, 2) /= array_sizes(1, 1, 2))) &
            mixed_cdft%identical_constraints = .FALSE.
         ! Check constraint definitions
         IF (mixed_cdft%identical_constraints) THEN
            ! Prepare temporary storage
            ALLOCATE (atoms(nforce_eval, settings%si(4, 1)))
            ALLOCATE (coeff(nforce_eval, settings%si(4, 1)))
            ALLOCATE (constraint_type(nforce_eval, settings%si(4, 1)))
            constraint_type(:, :) = 0
            DO iforce_eval = 1, nforce_eval
               DO i = 1, settings%si(4, 1)
                  NULLIFY (atoms(iforce_eval, i)%array)
                  NULLIFY (coeff(iforce_eval, i)%array)
                  ALLOCATE (atoms(iforce_eval, i)%array(array_sizes(iforce_eval, i, 1)))
                  ALLOCATE (coeff(iforce_eval, i)%array(array_sizes(iforce_eval, i, 1)))
                  atoms(iforce_eval, i)%array(:) = 0
                  coeff(iforce_eval, i)%array(:) = 0
               END DO
               ! Get constraint definitions
               IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
               force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
               IF (mixed_env%do_mixed_qmmm_cdft) THEN
                  qs_env => force_env_qs%qmmm_env%qs_env
               ELSE
                  CALL force_env_get(force_env_qs, qs_env=qs_env)
               END IF
               CALL get_qs_env(qs_env, dft_control=dft_control)
               cdft_control => dft_control%qs_control%cdft_control
               IF (force_env_qs%para_env%ionode) THEN
                  DO i = 1, settings%si(4, 1)
                     atoms(iforce_eval, i)%array(:) = cdft_control%group(i)%atoms
                     coeff(iforce_eval, i)%array(:) = cdft_control%group(i)%coeff
                     constraint_type(iforce_eval, i) = cdft_control%group(i)%constraint_type
                  END DO
               END IF
            END DO
            ! Sum up constraint definitions and check consistency
            DO i = 1, settings%si(4, 1)
               DO iforce_eval = 1, nforce_eval
                  CALL mp_sum(atoms(iforce_eval, i)%array, force_env%para_env%group)
                  CALL mp_sum(coeff(iforce_eval, i)%array, force_env%para_env%group)
                  CALL mp_sum(constraint_type(iforce_eval, i), force_env%para_env%group)
               END DO
               DO iforce_eval = 2, nforce_eval
                  DO iatom = 1, SIZE(atoms(1, i)%array)
                     IF (atoms(1, i)%array(iatom) /= atoms(iforce_eval, i)%array(iatom)) &
                        mixed_cdft%identical_constraints = .FALSE.
                     IF (coeff(1, i)%array(iatom) /= coeff(iforce_eval, i)%array(iatom)) &
                        mixed_cdft%identical_constraints = .FALSE.
                     IF (.NOT. mixed_cdft%identical_constraints) EXIT
                  END DO
                  IF (constraint_type(1, i) /= constraint_type(iforce_eval, i)) &
                     mixed_cdft%identical_constraints = .FALSE.
                  IF (.NOT. mixed_cdft%identical_constraints) EXIT
               END DO
               IF (.NOT. mixed_cdft%identical_constraints) EXIT
            END DO
            ! Deallocate temporary storage
            DO iforce_eval = 1, nforce_eval
               DO i = 1, settings%si(4, 1)
                  DEALLOCATE (atoms(iforce_eval, i)%array)
                  DEALLOCATE (coeff(iforce_eval, i)%array)
               END DO
            END DO
            DEALLOCATE (atoms)
            DEALLOCATE (coeff)
            DEALLOCATE (constraint_type)
         END IF
         DEALLOCATE (array_sizes)
      END IF
      ! Deallocate some arrays that are no longer needed
      IF (mixed_cdft%identical_constraints .AND. mixed_cdft%run_type /= mixed_cdft_parallel_nobuild) THEN
         DO iforce_eval = 1, nforce_eval
            IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
            force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
            IF (mixed_env%do_mixed_qmmm_cdft) THEN
               qs_env => force_env_qs%qmmm_env%qs_env
            ELSE
               CALL force_env_get(force_env_qs, qs_env=qs_env)
            END IF
            CALL get_qs_env(qs_env, dft_control=dft_control)
            cdft_control => dft_control%qs_control%cdft_control
            IF (mixed_cdft%run_type == mixed_cdft_parallel) THEN
               IF (.NOT. dft_control%qs_control%gapw) THEN
                  DO i = 1, SIZE(cdft_control%group)
                     DEALLOCATE (cdft_control%group(i)%coeff)
                     DEALLOCATE (cdft_control%group(i)%atoms)
                  END DO
                  IF (.NOT. cdft_control%atomic_charges) DEALLOCATE (cdft_control%atoms)
               END IF
            ELSE IF (mixed_cdft%run_type == mixed_cdft_serial) THEN
               IF (iforce_eval == 1) CYCLE
               DO igroup = 1, SIZE(cdft_control%group)
                  IF (.NOT. dft_control%qs_control%gapw) THEN
                     DEALLOCATE (cdft_control%group(igroup)%coeff)
                     DEALLOCATE (cdft_control%group(igroup)%atoms)
                  END IF
               END DO
               IF (cdft_control%type == outer_scf_becke_constraint) THEN
                  IF (.NOT. cdft_control%atomic_charges) DEALLOCATE (cdft_control%atoms)
                  IF (cdft_control%becke_control%cavity_confine) &
                     CALL release_hirshfeld_type(cdft_control%becke_control%cavity_env)
                  IF (cdft_control%becke_control%cutoff_type == becke_cutoff_element) &
                     DEALLOCATE (cdft_control%becke_control%cutoffs_tmp)
                  IF (cdft_control%becke_control%adjust) &
                     DEALLOCATE (cdft_control%becke_control%radii_tmp)
               END IF
            END IF
         END DO
      END IF

   END SUBROUTINE mixed_cdft_parse_settings

! **************************************************************************************************
!> \brief Transfer settings to mixed_cdft
!> \param force_env the force_env that holds the CDFT states
!> \param mixed_cdft the control section for mixed CDFT calculations
!> \param settings container for settings related to the mixed CDFT calculation
!> \par History
!>       01.2017  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_transfer_settings(force_env, mixed_cdft, settings)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(mixed_cdft_settings_type)                     :: settings

      INTEGER                                            :: i, nkinds
      LOGICAL                                            :: is_match
      TYPE(cdft_control_type), POINTER                   :: cdft_control

      NULLIFY (cdft_control)
      is_match = .TRUE.
      ! Transfer global settings
      mixed_cdft%multiplicity = settings%si(3, 1)
      mixed_cdft%has_unit_metric = settings%sb(7, 1)
      mixed_cdft%eps_rho_rspace = settings%sr(4, 1)
      mixed_cdft%nconstraint = settings%si(4, 1)
      settings%radius = settings%sr(5, 1)
      ! Transfer settings only needed if the constraint should be built in parallel
      IF (mixed_cdft%run_type == mixed_cdft_parallel) THEN
         IF (settings%sb(6, 1)) &
            CALL cp_abort(__LOCATION__, &
                          "Calculation of atomic Becke charges not supported with parallel mode mixed CDFT")
         IF (mixed_cdft%nconstraint /= 1) &
            CALL cp_abort(__LOCATION__, &
                          "Parallel mode mixed CDFT does not yet support multiple constraints.")

         IF (settings%si(5, 1) /= outer_scf_becke_constraint) &
            CALL cp_abort(__LOCATION__, &
                          "Parallel mode mixed CDFT does not support Hirshfeld constraints.")

         CALL cdft_control_create(mixed_cdft%cdft_control)
         cdft_control => mixed_cdft%cdft_control
         ALLOCATE (cdft_control%atoms(settings%ncdft))
         cdft_control%atoms = settings%atoms(1:settings%ncdft, 1)
         ALLOCATE (cdft_control%group(1))
         ALLOCATE (cdft_control%group(1)%atoms(settings%ncdft))
         ALLOCATE (cdft_control%group(1)%coeff(settings%ncdft))
         NULLIFY (cdft_control%group(1)%weight%pw)
         NULLIFY (cdft_control%group(1)%gradients)
         NULLIFY (cdft_control%group(1)%integrated)
         cdft_control%group(1)%atoms = cdft_control%atoms
         cdft_control%group(1)%coeff = settings%coeffs(1:settings%ncdft, 1)
         cdft_control%natoms = settings%ncdft
         cdft_control%atomic_charges = settings%sb(6, 1)
         cdft_control%becke_control%cutoff_type = settings%si(1, 1)
         cdft_control%becke_control%cavity_confine = settings%sb(1, 1)
         cdft_control%becke_control%should_skip = settings%sb(2, 1)
         cdft_control%becke_control%print_cavity = settings%sb(3, 1)
         cdft_control%becke_control%in_memory = settings%sb(4, 1)
         cdft_control%becke_control%adjust = settings%sb(5, 1)
         cdft_control%becke_control%cavity_shape = settings%si(2, 1)
         cdft_control%becke_control%use_bohr = settings%sb(8, 1)
         cdft_control%becke_control%rcavity = settings%sr(1, 1)
         cdft_control%becke_control%rglobal = settings%sr(2, 1)
         cdft_control%becke_control%eps_cavity = settings%sr(3, 1)
         nkinds = 0
         IF (cdft_control%becke_control%cutoff_type == becke_cutoff_element) THEN
            CALL mp_sum(settings%cutoffs, force_env%para_env%group)
            DO i = 1, SIZE(settings%cutoffs, 1)
               IF (settings%cutoffs(i, 1) /= settings%cutoffs(i, 2)) is_match = .FALSE.
               IF (settings%cutoffs(i, 1) /= 0.0_dp) nkinds = nkinds + 1
            END DO
            IF (.NOT. is_match) &
               CALL cp_abort(__LOCATION__, &
                             "Mismatch detected in the &BECKE_CONSTRAINT"// &
                             "&ELEMENT_CUTOFF settings of the two force_evals.")
            ALLOCATE (cdft_control%becke_control%cutoffs_tmp(nkinds))
            cdft_control%becke_control%cutoffs_tmp = settings%cutoffs(1:nkinds, 1)
         END IF
         nkinds = 0
         IF (cdft_control%becke_control%adjust) THEN
            CALL mp_sum(settings%radii, force_env%para_env%group)
            DO i = 1, SIZE(settings%radii, 1)
               IF (settings%radii(i, 1) /= settings%radii(i, 2)) is_match = .FALSE.
               IF (settings%radii(i, 1) /= 0.0_dp) nkinds = nkinds + 1
            END DO
            IF (.NOT. is_match) &
               CALL cp_abort(__LOCATION__, &
                             "Mismatch detected in the &BECKE_CONSTRAINT"// &
                             "&ATOMIC_RADII settings of the two force_evals.")
            ALLOCATE (cdft_control%becke_control%radii(nkinds))
            cdft_control%becke_control%radii = settings%radii(1:nkinds, 1)
         END IF
      END IF

   END SUBROUTINE mixed_cdft_transfer_settings

! **************************************************************************************************
!> \brief Initialize all the structures needed for a mixed CDFT calculation
!> \param force_env the force_env that holds the CDFT mixed_env
!> \param force_env_qs the force_env that holds the qs_env, which is CDFT state specific
!> \param mixed_env the mixed_env that holds the CDFT states
!> \param mixed_cdft the control section for mixed CDFT calculations
!> \param settings container for settings related to the mixed CDFT calculation
!> \par History
!>       01.2017 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_init_structures(force_env, force_env_qs, mixed_env, mixed_cdft, settings)
      TYPE(force_env_type), POINTER                      :: force_env, force_env_qs
      TYPE(mixed_environment_type), POINTER              :: mixed_env
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(mixed_cdft_settings_type)                     :: settings

      CHARACTER(len=default_path_length)                 :: c_val, input_file_path, output_file_path
      INTEGER                                            :: i, imap, iounit, j, lp, n_force_eval, &
                                                            ncpu, nforce_eval, ntargets, offset, &
                                                            req(3), unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: bounds
      INTEGER, DIMENSION(2, 3)                           :: bo, bo_mixed
      INTEGER, DIMENSION(3)                              :: higher_grid_layout
      INTEGER, DIMENSION(:), POINTER                     :: i_force_eval, mixed_rs_dims, recvbuffer, &
                                                            recvbuffer2, sendbuffer
      LOGICAL                                            :: is_match
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell_mix
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys_mix
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
         POINTER                                         :: rs_descs
      TYPE(realspace_grid_input_type)                    :: input_settings
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_grids
      TYPE(section_vals_type), POINTER :: force_env_section, force_env_sections, kind_section, &
         print_section, root_section, rs_grid_section, subsys_section

      NULLIFY (cell_mix, subsys_mix, force_env_section, subsys_section, &
               print_section, root_section, kind_section, force_env_sections, rs_grid_section, &
               auxbas_pw_pool, pw_env, pw_pools, pw_grid, sendbuffer, qs_env, &
               mixed_rs_dims, i_force_eval, recvbuffer, recvbuffer2, &
               blacs_env, globenv, atomic_kind_set, qs_kind_set, rs_descs, rs_grids)

      logger => cp_get_default_logger()
      CALL force_env_get(force_env=force_env, force_env_section=force_env_section)
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
      iounit = cp_print_key_unit_nr(logger, print_section, '', extension='.mixedLog')
      is_match = .TRUE.
      nforce_eval = SIZE(force_env%sub_force_env)
      ncpu = force_env%para_env%num_pe
      ! Get infos about the mixed subsys
      IF (.NOT. mixed_env%do_mixed_qmmm_cdft) THEN
         CALL force_env_get(force_env=force_env, &
                            subsys=subsys_mix)
      ELSE
         CALL get_qs_env(force_env_qs%qmmm_env%qs_env, &
                         cp_subsys=subsys_mix)
      END IF
      ! Init structures only needed when the CDFT states are treated in parallel
      IF (mixed_cdft%run_type == mixed_cdft_parallel) THEN
         ! Start building the mixed auxbas_pw_pool
         CALL pw_grid_create(pw_grid, force_env%para_env%group)
         CALL pw_env_create(mixed_cdft%pw_env)
         ! Decide what kind of layout to use and setup the grid
         ! Processor mappings currently supported:
         !  (2np,1)  --> (np,1)
         !  (nx,2ny) --> (nx,ny)
         !  (nx,ny)  --> (nx*ny/2,1) (required when xc_smooth is in use and with intermediate proc counts)
         !
         ! For cases 2 and 3, dlb redistributes YZ slices from overloaded processors to underloaded processors
         ! For case 1, XZ slices are redistributed
         ! TODO: Unify mappings. Now we essentially have separate code for cases 1-2 and 3.
         !       This leads to very messy code especially with dlb turned on...
         !       In terms of memory usage, it would be beneficial to replace case 1 with 3
         !       and implement a similar arbitrary mapping to replace case 2

         mixed_cdft%is_pencil = .FALSE. ! Flag to control the first two mappings
         mixed_cdft%is_special = .FALSE. ! Flag to control the last mapping
         ! With xc smoothing, the grid is always (ncpu/2,1) distributed
         ! and correct behavior cannot be guaranteed for ncpu/2 > nx, so we abort...
         IF (ncpu/2 .GT. settings%npts(1, 1)) &
            CPABORT("ncpu/2 => nx: decrease ncpu or disable xc_smoothing")
         !
         ALLOCATE (mixed_rs_dims(2))
         IF (settings%rs_dims(2, 1) /= 1) mixed_cdft%is_pencil = .TRUE.
         IF (.NOT. mixed_cdft%is_pencil .AND. ncpu .GT. settings%npts(1, 1)) mixed_cdft%is_special = .TRUE.
         IF (mixed_cdft%is_special) THEN
            mixed_rs_dims = (/-1, -1/)
         ELSE IF (mixed_cdft%is_pencil) THEN
            mixed_rs_dims = (/settings%rs_dims(1, 1), 2*settings%rs_dims(2, 1)/)
         ELSE
            mixed_rs_dims = (/2*settings%rs_dims(1, 1), 1/)
         END IF
         IF (.NOT. mixed_env%do_mixed_qmmm_cdft) THEN
            CALL force_env_get(force_env=force_env, &
                               cell=cell_mix)
         ELSE
            CALL get_qs_env(force_env_qs%qmmm_env%qs_env, &
                            cell=cell_mix)
         END IF
         CALL pw_grid_setup(cell_mix%hmat, pw_grid, grid_span=settings%grid_span(1), &
                            npts=settings%npts(:, 1), cutoff=settings%cutoff(1), &
                            spherical=settings%is_spherical, odd=settings%is_odd, &
                            fft_usage=.TRUE., ncommensurate=0, icommensurate=1, &
                            blocked=do_pw_grid_blocked_false, rs_dims=mixed_rs_dims, &
                            iounit=iounit)
         ! Check if the layout was successfully created
         IF (mixed_cdft%is_special) THEN
            IF (.NOT. pw_grid%para%rs_dims(2) /= 1) is_match = .FALSE.
         ELSE IF (mixed_cdft%is_pencil) THEN
            IF (.NOT. pw_grid%para%rs_dims(1) == mixed_rs_dims(1)) is_match = .FALSE.
         ELSE
            IF (.NOT. pw_grid%para%rs_dims(2) == 1) is_match = .FALSE.
         END IF
         IF (.NOT. is_match) &
            CALL cp_abort(__LOCATION__, &
                          "Unable to create a suitable grid distribution "// &
                          "for mixed CDFT calculations. Try decreasing the total number "// &
                          " of processors or disabling xc_smoothing.")
         DEALLOCATE (mixed_rs_dims)
         ! Create the pool
         bo_mixed = pw_grid%bounds_local
         ALLOCATE (pw_pools(1))
         NULLIFY (pw_pools(1)%pool)
         CALL pw_pool_create(pw_pools(1)%pool, pw_grid=pw_grid)
         ! Initialize Gaussian cavity confinement
         IF (mixed_cdft%cdft_control%becke_control%cavity_confine) THEN
            CALL create_hirshfeld_type(mixed_cdft%cdft_control%becke_control%cavity_env)
            CALL set_hirshfeld_info(mixed_cdft%cdft_control%becke_control%cavity_env, &
                                    shape_function_type=shape_function_gaussian, iterative=.FALSE., &
                                    radius_type=mixed_cdft%cdft_control%becke_control%cavity_shape, &
                                    use_bohr=mixed_cdft%cdft_control%becke_control%use_bohr)
         END IF
         ! Gaussian confinement/wavefunction overlap method needs qs_kind_set
         ! Gaussian cavity confinement also needs the auxbas_rs_grid
         IF (mixed_cdft%cdft_control%becke_control%cavity_confine .OR. &
             mixed_cdft%wfn_overlap_method) THEN
            print_section => section_vals_get_subs_vals(force_env_section, &
                                                        "PRINT%GRID_INFORMATION")
            ALLOCATE (mixed_cdft%pw_env%gridlevel_info)
            CALL init_gaussian_gridlevel(mixed_cdft%pw_env%gridlevel_info, &
                                         ngrid_levels=1, cutoff=settings%cutoff, &
                                         rel_cutoff=settings%rel_cutoff(1), &
                                         print_section=print_section)
            ALLOCATE (rs_descs(1))
            ALLOCATE (rs_grids(1))
            ALLOCATE (mixed_cdft%pw_env%cube_info(1))
            higher_grid_layout = (/-1, -1, -1/)
            CALL init_d3_poly_module()
            CALL init_cube_info(mixed_cdft%pw_env%cube_info(1), &
                                pw_grid%dr(:), pw_grid%dh(:, :), &
                                pw_grid%dh_inv(:, :), &
                                pw_grid%orthorhombic, settings%radius)
            NULLIFY (root_section, force_env_section, force_env_sections, rs_grid_section)
            CALL force_env_get(force_env, root_section=root_section)
            force_env_sections => section_vals_get_subs_vals(root_section, "FORCE_EVAL")
            CALL multiple_fe_list(force_env_sections, root_section, i_force_eval, n_force_eval)
            CALL section_vals_duplicate(force_env_sections, force_env_section, &
                                        i_force_eval(2), i_force_eval(2))
            rs_grid_section => section_vals_get_subs_vals(force_env_section, "DFT%MGRID%RS_GRID")
            CALL init_input_type(input_settings, &
                                 nsmax=2*MAX(1, return_cube_max_iradius(mixed_cdft%pw_env%cube_info(1))) + 1, &
                                 rs_grid_section=rs_grid_section, ilevel=1, &
                                 higher_grid_layout=higher_grid_layout)
            NULLIFY (rs_descs(1)%rs_desc)
            CALL rs_grid_create_descriptor(rs_descs(1)%rs_desc, pw_grid, input_settings)
            IF (rs_descs(1)%rs_desc%distributed) higher_grid_layout = rs_descs(1)%rs_desc%group_dim
            NULLIFY (rs_grids(1)%rs_grid)
            CALL rs_grid_create(rs_grids(1)%rs_grid, rs_descs(1)%rs_desc)
            CALL rs_grid_print(rs_grids(1)%rs_grid, iounit)
            mixed_cdft%pw_env%rs_descs => rs_descs
            mixed_cdft%pw_env%rs_grids => rs_grids
            ! qs_kind_set
            subsys_section => section_vals_get_subs_vals(force_env_sections, "SUBSYS", &
                                                         i_rep_section=i_force_eval(1))
            kind_section => section_vals_get_subs_vals(subsys_section, "KIND")
            NULLIFY (qs_kind_set)
            CALL cp_subsys_get(subsys_mix, atomic_kind_set=atomic_kind_set)
            CALL create_qs_kind_set(qs_kind_set, atomic_kind_set, kind_section, &
                                    force_env%para_env, force_env_section)
            mixed_cdft%qs_kind_set => qs_kind_set
            DEALLOCATE (i_force_eval)
            CALL section_vals_release(force_env_section)
         END IF
         CALL force_env_get(force_env=force_env, &
                            force_env_section=force_env_section)
         CALL pw_grid_release(pw_grid)
         mixed_cdft%pw_env%auxbas_grid = 1
         NULLIFY (mixed_cdft%pw_env%pw_pools)
         mixed_cdft%pw_env%pw_pools => pw_pools
         bo = settings%bo
         ! Determine which processors need to exchange data when redistributing the weight/gradient
         IF (.NOT. mixed_cdft%is_special) THEN
            ALLOCATE (mixed_cdft%dest_list(2))
            ALLOCATE (mixed_cdft%source_list(2))
            imap = force_env%para_env%mepos/2
            mixed_cdft%dest_list = (/imap, imap + force_env%para_env%num_pe/2/)
            imap = MOD(force_env%para_env%mepos, force_env%para_env%num_pe/2) + &
                   MODULO(force_env%para_env%mepos, force_env%para_env%num_pe/2)
            mixed_cdft%source_list = (/imap, imap + 1/)
            ! Determine bounds of the data that is replicated
            ALLOCATE (mixed_cdft%recv_bo(4))
            ALLOCATE (sendbuffer(2), recvbuffer(2), recvbuffer2(2))
            IF (mixed_cdft%is_pencil) THEN
               sendbuffer = (/bo_mixed(1, 2), bo_mixed(2, 2)/)
            ELSE
               sendbuffer = (/bo_mixed(1, 1), bo_mixed(2, 1)/)
            END IF
            ! Communicate bounds in steps
            CALL mp_isend(msgin=sendbuffer, dest=mixed_cdft%dest_list(1), &
                          request=req(1), comm=force_env%para_env%group)
            CALL mp_irecv(msgout=recvbuffer, source=mixed_cdft%source_list(1), &
                          request=req(2), comm=force_env%para_env%group)
            CALL mp_irecv(msgout=recvbuffer2, source=mixed_cdft%source_list(2), &
                          request=req(3), comm=force_env%para_env%group)
            CALL mp_wait(req(1))
            CALL mp_isend(msgin=sendbuffer, dest=mixed_cdft%dest_list(2), &
                          request=req(1), comm=force_env%para_env%group)
            CALL mp_waitall(req)
            mixed_cdft%recv_bo(1:2) = recvbuffer
            mixed_cdft%recv_bo(3:4) = recvbuffer2
            DEALLOCATE (sendbuffer, recvbuffer, recvbuffer2)
         ELSE
            IF (mixed_env%do_mixed_qmmm_cdft) THEN
               qs_env => force_env_qs%qmmm_env%qs_env
            ELSE
               CALL force_env_get(force_env_qs, qs_env=qs_env)
            END IF
            CALL get_qs_env(qs_env, pw_env=pw_env)
            CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
            ! work out the pw grid points each proc holds in the two (identical) parallel proc groups
            ! note we only care about the x dir since we assume the y dir is not subdivided
            ALLOCATE (bounds(0:auxbas_pw_pool%pw_grid%para%group_size - 1, 1:2))
            DO i = 0, auxbas_pw_pool%pw_grid%para%group_size - 1
               bounds(i, 1:2) = auxbas_pw_pool%pw_grid%para%bo(1:2, 1, i, 1)
               bounds(i, 1:2) = bounds(i, 1:2) - auxbas_pw_pool%pw_grid%npts(1)/2 - 1
            END DO
            ! work out which procs to send my grid points
            ! first get the number of target procs per group
            ntargets = 0
            offset = -1
            DO i = 0, auxbas_pw_pool%pw_grid%para%group_size - 1
               IF ((bounds(i, 1) .GE. bo_mixed(1, 1) .AND. bounds(i, 1) .LE. bo_mixed(2, 1)) .OR. &
                   (bounds(i, 2) .GE. bo_mixed(1, 1) .AND. bounds(i, 2) .LE. bo_mixed(2, 1))) THEN
                  ntargets = ntargets + 1
                  IF (offset == -1) offset = i
               ELSE IF (bounds(i, 2) .GT. bo_mixed(2, 1)) THEN
                  EXIT
               ELSE
                  CYCLE
               END IF
            END DO
            ALLOCATE (mixed_cdft%dest_list(ntargets))
            ALLOCATE (mixed_cdft%dest_list_bo(2, ntargets))
            ! now determine the actual grid points to send
            j = 1
            DO i = offset, offset + ntargets - 1
               mixed_cdft%dest_list(j) = i
               mixed_cdft%dest_list_bo(:, j) = (/bo_mixed(1, 1) + (bounds(i, 1) - bo_mixed(1, 1)), &
                                                 bo_mixed(2, 1) + (bounds(i, 2) - bo_mixed(2, 1))/)
               j = j + 1
            END DO
            ALLOCATE (mixed_cdft%dest_list_save(ntargets), mixed_cdft%dest_bo_save(2, ntargets))
            ! We need to store backups of these arrays since they might get reallocated during dlb
            mixed_cdft%dest_list_save = mixed_cdft%dest_list
            mixed_cdft%dest_bo_save = mixed_cdft%dest_list_bo
            ! finally determine which procs will send me grid points
            ! now we need info about y dir also
            DEALLOCATE (bounds)
            ALLOCATE (bounds(0:pw_pools(1)%pool%pw_grid%para%group_size - 1, 1:4))
            DO i = 0, pw_pools(1)%pool%pw_grid%para%group_size - 1
               bounds(i, 1:2) = pw_pools(1)%pool%pw_grid%para%bo(1:2, 1, i, 1)
               bounds(i, 3:4) = pw_pools(1)%pool%pw_grid%para%bo(1:2, 2, i, 1)
               bounds(i, 1:2) = bounds(i, 1:2) - pw_pools(1)%pool%pw_grid%npts(1)/2 - 1
               bounds(i, 3:4) = bounds(i, 3:4) - pw_pools(1)%pool%pw_grid%npts(2)/2 - 1
            END DO
            ntargets = 0
            offset = -1
            DO i = 0, pw_pools(1)%pool%pw_grid%para%group_size - 1
               IF ((bo(1, 1) .GE. bounds(i, 1) .AND. bo(1, 1) .LE. bounds(i, 2)) .OR. &
                   (bo(2, 1) .GE. bounds(i, 1) .AND. bo(2, 1) .LE. bounds(i, 2))) THEN
                  ntargets = ntargets + 1
                  IF (offset == -1) offset = i
               ELSE IF (bo(2, 1) .LT. bounds(i, 1)) THEN
                  EXIT
               ELSE
                  CYCLE
               END IF
            END DO
            ALLOCATE (mixed_cdft%source_list(ntargets))
            ALLOCATE (mixed_cdft%source_list_bo(4, ntargets))
            j = 1
            DO i = offset, offset + ntargets - 1
               mixed_cdft%source_list(j) = i
               IF (bo(1, 1) .GE. bounds(i, 1) .AND. bo(2, 1) .LE. bounds(i, 2)) THEN
                  mixed_cdft%source_list_bo(:, j) = (/bo(1, 1), bo(2, 1), &
                                                      bounds(i, 3), bounds(i, 4)/)
               ELSE IF (bo(1, 1) .GE. bounds(i, 1) .AND. bo(1, 1) .LE. bounds(i, 2)) THEN
                  mixed_cdft%source_list_bo(:, j) = (/bo(1, 1), bounds(i, 2), &
                                                      bounds(i, 3), bounds(i, 4)/)
               ELSE
                  mixed_cdft%source_list_bo(:, j) = (/bounds(i, 1), bo(2, 1), &
                                                      bounds(i, 3), bounds(i, 4)/)
               END IF
               j = j + 1
            END DO
            ALLOCATE (mixed_cdft%source_list_save(ntargets), mixed_cdft%source_bo_save(4, ntargets))
            ! We need to store backups of these arrays since they might get reallocated during dlb
            mixed_cdft%source_list_save = mixed_cdft%source_list
            mixed_cdft%source_bo_save = mixed_cdft%source_list_bo
            DEALLOCATE (bounds)
         END IF
      ELSE
         ! Create loggers to redirect the output of all CDFT states to different files
         ! even when the states are treated in serial (the initial print of QS data [basis set etc] for
         ! all states unfortunately goes to the first log file)
         CALL force_env_get(force_env, root_section=root_section)
         ALLOCATE (mixed_cdft%sub_logger(nforce_eval - 1))
         DO i = 1, nforce_eval - 1
            IF (force_env%para_env%ionode) THEN
               CALL section_vals_val_get(root_section, "GLOBAL%PROJECT_NAME", &
                                         c_val=input_file_path)
               lp = LEN_TRIM(input_file_path)
               input_file_path(lp + 1:LEN(input_file_path)) = "-r-"//ADJUSTL(cp_to_string(i + 1))
               lp = LEN_TRIM(input_file_path)
               output_file_path = input_file_path(1:lp)//".out"
               CALL open_file(file_name=output_file_path, file_status="UNKNOWN", &
                              file_action="WRITE", file_position="APPEND", &
                              unit_number=unit_nr)
            ELSE
               unit_nr = -1
            END IF
            CALL cp_logger_create(mixed_cdft%sub_logger(i)%p, &
                                  para_env=force_env%para_env, &
                                  default_global_unit_nr=unit_nr, &
                                  close_global_unit_on_dealloc=.FALSE.)
            ! Try to use better names for the local log if it is not too late
            CALL section_vals_val_get(root_section, "GLOBAL%OUTPUT_FILE_NAME", &
                                      c_val=c_val)
            IF (c_val /= "") THEN
               CALL cp_logger_set(mixed_cdft%sub_logger(i)%p, &
                                  local_filename=TRIM(c_val)//"_localLog")
            END IF
            CALL section_vals_val_get(root_section, "GLOBAL%PROJECT", c_val=c_val)
            IF (c_val /= "") THEN
               CALL cp_logger_set(mixed_cdft%sub_logger(i)%p, &
                                  local_filename=TRIM(c_val)//"_localLog")
            END IF
            mixed_cdft%sub_logger(i)%p%iter_info%project_name = c_val
            CALL section_vals_val_get(root_section, "GLOBAL%PRINT_LEVEL", &
                                      i_val=mixed_cdft%sub_logger(i)%p%iter_info%print_level)
         END DO
         IF (mixed_cdft%wfn_overlap_method) THEN
            ! qs_kind_set
            NULLIFY (root_section, force_env_section, force_env_sections, rs_grid_section)
            CALL force_env_get(force_env, root_section=root_section)
            force_env_sections => section_vals_get_subs_vals(root_section, "FORCE_EVAL")
            CALL multiple_fe_list(force_env_sections, root_section, i_force_eval, n_force_eval)
            CALL section_vals_duplicate(force_env_sections, force_env_section, &
                                        i_force_eval(2), i_force_eval(2))
            subsys_section => section_vals_get_subs_vals(force_env_sections, "SUBSYS", &
                                                         i_rep_section=i_force_eval(1))
            kind_section => section_vals_get_subs_vals(subsys_section, "KIND")
            NULLIFY (qs_kind_set)
            CALL cp_subsys_get(subsys_mix, atomic_kind_set=atomic_kind_set)
            CALL create_qs_kind_set(qs_kind_set, atomic_kind_set, kind_section, &
                                    force_env%para_env, force_env_section)
            mixed_cdft%qs_kind_set => qs_kind_set
            DEALLOCATE (i_force_eval)
            CALL section_vals_release(force_env_section)
            mixed_cdft%qs_kind_set => qs_kind_set
         END IF
         CALL force_env_get(force_env=force_env, &
                            force_env_section=force_env_section)
      END IF
      ! Deallocate settings temporaries
      DEALLOCATE (settings%grid_span)
      DEALLOCATE (settings%npts)
      DEALLOCATE (settings%spherical)
      DEALLOCATE (settings%rs_dims)
      DEALLOCATE (settings%odd)
      DEALLOCATE (settings%atoms)
      IF (mixed_cdft%run_type == mixed_cdft_parallel) &
         DEALLOCATE (settings%coeffs)
      DEALLOCATE (settings%cutoffs)
      DEALLOCATE (settings%radii)
      DEALLOCATE (settings%si)
      DEALLOCATE (settings%sr)
      DEALLOCATE (settings%sb)
      DEALLOCATE (settings%cutoff)
      DEALLOCATE (settings%rel_cutoff)
      ! Setup mixed blacs_env for redistributing arrays during ET coupling calculation
      IF (mixed_env%do_mixed_et) THEN
         NULLIFY (root_section)
         CALL force_env_get(force_env, globenv=globenv, root_section=root_section)
         CALL cp_blacs_env_create(blacs_env, force_env%para_env, globenv%blacs_grid_layout, &
                                  globenv%blacs_repeatable)
         mixed_cdft%blacs_env => blacs_env
         CALL cp_blacs_env_retain(mixed_cdft%blacs_env)
         CALL cp_blacs_env_release(blacs_env)
      END IF
      CALL cp_print_key_finished_output(iounit, logger, force_env_section, &
                                        "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")

   END SUBROUTINE mixed_cdft_init_structures

! **************************************************************************************************
!> \brief Redistribute arrays needed for an ET coupling calculation from individual CDFT states to
!>        the mixed CDFT env, that is, move the arrays to the correct blacs context. For parallel
!>        simulations, the array processor distributions also change from N to 2N processors.
!> \param force_env the force_env that holds the CDFT states
!> \par History
!>       01.2017  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_redistribute_arrays(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      INTEGER                                            :: iforce_eval, ispin, ivar, ncol_overlap, &
                                                            ncol_wmat, nforce_eval, nrow_overlap, &
                                                            nrow_wmat, nspins, nvar
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ncol_mo, nrow_mo
      LOGICAL                                            :: uniform_occupation
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: has_occupation_numbers
      TYPE(cp_1d_r_p_type), ALLOCATABLE, DIMENSION(:, :) :: occno_tmp
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: matrix_p_tmp, mixed_matrix_p_tmp, &
                                                            mixed_mo_coeff, mixed_wmat_tmp, &
                                                            mo_coeff_tmp, wmat_tmp
      TYPE(cp_fm_p_type), POINTER                        :: matrix_s_tmp, mixed_matrix_s_tmp
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_mo, fm_struct_overlap, &
                                                            fm_struct_tmp, fm_struct_wmat
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: density_matrix, w_matrix
      TYPE(dbcsr_type)                                   :: desymm_tmp
      TYPE(dbcsr_type), POINTER                          :: mixed_matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(force_env_type), POINTER                      :: force_env_qs
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(mixed_environment_type), POINTER              :: mixed_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      NULLIFY (mixed_env, mixed_cdft, qs_env, dft_control, fm_struct_mo, &
               fm_struct_wmat, fm_struct_overlap, fm_struct_tmp, &
               mixed_matrix_s_tmp, matrix_s_tmp, mixed_wmat_tmp, &
               wmat_tmp, mixed_mo_coeff, matrix_p_tmp, mixed_matrix_p_tmp, &
               mixed_matrix_s, density_matrix, blacs_env, w_matrix, force_env_qs)
      CPASSERT(ASSOCIATED(force_env))
      mixed_env => force_env%mixed_env
      nforce_eval = SIZE(force_env%sub_force_env)
      CALL get_mixed_env(mixed_env, cdft_control=mixed_cdft)
      CPASSERT(ASSOCIATED(mixed_cdft))
      CALL mixed_cdft_work_type_init(mixed_cdft%matrix)
      ! Get nspins and query for non-uniform occupation numbers
      ALLOCATE (has_occupation_numbers(nforce_eval))
      has_occupation_numbers = .FALSE.
      DO iforce_eval = 1, nforce_eval
         IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
         force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
         IF (force_env%mixed_env%do_mixed_qmmm_cdft) THEN
            qs_env => force_env_qs%qmmm_env%qs_env
         ELSE
            CALL force_env_get(force_env_qs, qs_env=qs_env)
         END IF
         CALL get_qs_env(qs_env, dft_control=dft_control)
         CPASSERT(ASSOCIATED(dft_control))
         nspins = dft_control%nspins
         IF (force_env_qs%para_env%ionode) &
            has_occupation_numbers(iforce_eval) = ALLOCATED(dft_control%qs_control%cdft_control%occupations)
      END DO
      CALL mp_sum(has_occupation_numbers(1), force_env%para_env%group)
      DO iforce_eval = 2, nforce_eval
         CALL mp_sum(has_occupation_numbers(iforce_eval), force_env%para_env%group)
         IF (has_occupation_numbers(1) .NEQV. has_occupation_numbers(iforce_eval)) &
            CALL cp_abort(__LOCATION__, &
                          "Mixing of uniform and non-uniform occupations is not allowed.")
      END DO
      uniform_occupation = .NOT. has_occupation_numbers(1)
      DEALLOCATE (has_occupation_numbers)
      ! Get number of weight functions per state as well as the type of each constraint
      nvar = SIZE(dft_control%qs_control%cdft_control%target)
      IF (.NOT. ALLOCATED(mixed_cdft%constraint_type)) THEN
         ALLOCATE (mixed_cdft%constraint_type(nvar, nforce_eval))
         mixed_cdft%constraint_type(:, :) = 0
         IF (mixed_cdft%identical_constraints) THEN
            DO ivar = 1, nvar
               mixed_cdft%constraint_type(ivar, :) = &
                  dft_control%qs_control%cdft_control%group(ivar)%constraint_type
            END DO
         ELSE
            ! Possibly couple spin and charge constraints
            DO iforce_eval = 1, nforce_eval
               IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
               IF (force_env%mixed_env%do_mixed_qmmm_cdft) THEN
                  qs_env => force_env%sub_force_env(iforce_eval)%force_env%qmmm_env%qs_env
               ELSE
                  CALL force_env_get(force_env%sub_force_env(iforce_eval)%force_env, qs_env=qs_env)
               END IF
               CALL get_qs_env(qs_env, dft_control=dft_control)
               IF (force_env%sub_force_env(iforce_eval)%force_env%para_env%ionode) THEN
                  DO ivar = 1, nvar
                     mixed_cdft%constraint_type(ivar, iforce_eval) = &
                        dft_control%qs_control%cdft_control%group(ivar)%constraint_type
                  END DO
               END IF
            END DO
            CALL mp_sum(mixed_cdft%constraint_type, force_env%para_env%group)
         END IF
      END IF
      ! Transfer data from sub_force_envs to temporaries
      ALLOCATE (mixed_cdft%matrix%mixed_mo_coeff(nforce_eval, nspins))
      mixed_mo_coeff => mixed_cdft%matrix%mixed_mo_coeff
      ALLOCATE (mixed_cdft%matrix%w_matrix(nforce_eval, nvar))
      w_matrix => mixed_cdft%matrix%w_matrix
      CALL dbcsr_init_p(mixed_cdft%matrix%mixed_matrix_s)
      mixed_matrix_s => mixed_cdft%matrix%mixed_matrix_s
      IF (mixed_cdft%calculate_metric) THEN
         ALLOCATE (mixed_cdft%matrix%density_matrix(nforce_eval, nspins))
         density_matrix => mixed_cdft%matrix%density_matrix
      END IF
      ALLOCATE (mo_coeff_tmp(nforce_eval, nspins), wmat_tmp(nforce_eval, nvar))
      ALLOCATE (matrix_s_tmp, nrow_mo(nspins), ncol_mo(nspins))
      NULLIFY (matrix_s_tmp%matrix)
      IF (mixed_cdft%calculate_metric) ALLOCATE (matrix_p_tmp(nforce_eval, nspins))
      IF (.NOT. uniform_occupation) THEN
         ALLOCATE (mixed_cdft%occupations(nforce_eval, nspins))
         ALLOCATE (occno_tmp(nforce_eval, nspins))
      END IF
      DO iforce_eval = 1, nforce_eval
         ! Temporary arrays need to be nulled on every process
         DO ivar = 1, nvar
            NULLIFY (wmat_tmp(iforce_eval, ivar)%matrix)
         END DO
         DO ispin = 1, nspins
            NULLIFY (mo_coeff_tmp(iforce_eval, ispin)%matrix)
            ! Valgrind 3.12/gfortran 4.8.4 oddly complains here (unconditional jump)
            ! if mixed_cdft%calculate_metric = .FALSE. and the need to null the array
            ! is queried with IF (mixed_cdft%calculate_metric) &
            IF (ASSOCIATED(matrix_p_tmp)) &
               NULLIFY (matrix_p_tmp(iforce_eval, ispin)%matrix)
            IF (.NOT. uniform_occupation) &
               NULLIFY (occno_tmp(iforce_eval, ispin)%array)
         END DO
         IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
         ! From this point onward, we access data local to the sub_force_envs
         ! Get qs_env
         force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
         IF (force_env%mixed_env%do_mixed_qmmm_cdft) THEN
            qs_env => force_env_qs%qmmm_env%qs_env
         ELSE
            CALL force_env_get(force_env_qs, qs_env=qs_env)
         END IF
         CALL get_qs_env(qs_env, dft_control=dft_control, blacs_env=blacs_env)
         ! Store dimensions of the transferred arrays
         CALL dbcsr_get_info(dft_control%qs_control%cdft_control%matrix_s%matrix, &
                             nfullrows_total=nrow_overlap, nfullcols_total=ncol_overlap)
         CALL dbcsr_get_info(dft_control%qs_control%cdft_control%wmat(1)%matrix, &
                             nfullrows_total=nrow_wmat, nfullcols_total=ncol_wmat)
         ! MO Coefficients
         DO ispin = 1, nspins
            CALL cp_fm_get_info(dft_control%qs_control%cdft_control%mo_coeff(ispin)%matrix, &
                                ncol_global=ncol_mo(ispin), nrow_global=nrow_mo(ispin))
            CALL cp_fm_create(matrix=mo_coeff_tmp(iforce_eval, ispin)%matrix, &
                              matrix_struct=dft_control%qs_control%cdft_control%mo_coeff(ispin)%matrix%matrix_struct, &
                              name="MO_COEFF_"//TRIM(ADJUSTL(cp_to_string(iforce_eval)))//"_" &
                              //TRIM(ADJUSTL(cp_to_string(ispin)))//"_MATRIX")
            CALL cp_fm_to_fm(dft_control%qs_control%cdft_control%mo_coeff(ispin)%matrix, &
                             mo_coeff_tmp(iforce_eval, ispin)%matrix)
            CALL cp_fm_release(dft_control%qs_control%cdft_control%mo_coeff(ispin)%matrix)
         END DO
         DEALLOCATE (dft_control%qs_control%cdft_control%mo_coeff)
         ! Matrix representation(s) of the weight function(s) (dbcsr -> fm)
         CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nrow_wmat, ncol_global=ncol_wmat, context=blacs_env, &
                                  para_env=force_env%sub_force_env(iforce_eval)%force_env%para_env, &
                                  square_blocks=.TRUE.)
         DO ivar = 1, nvar
            CALL cp_fm_create(wmat_tmp(iforce_eval, ivar)%matrix, fm_struct_tmp, name="w_matrix")
            CALL dbcsr_desymmetrize(dft_control%qs_control%cdft_control%wmat(ivar)%matrix, desymm_tmp)
            CALL copy_dbcsr_to_fm(desymm_tmp, wmat_tmp(iforce_eval, ivar)%matrix)
            CALL dbcsr_release(desymm_tmp)
            CALL dbcsr_release_p(dft_control%qs_control%cdft_control%wmat(ivar)%matrix)
         END DO
         DEALLOCATE (dft_control%qs_control%cdft_control%wmat)
         CALL cp_fm_struct_release(fm_struct_tmp)
         ! Overlap matrix is the same for all sub_force_envs, so we just copy the first one (dbcsr -> fm)
         IF (iforce_eval == 1) THEN
            CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nrow_overlap, &
                                     ncol_global=ncol_overlap, context=blacs_env, &
                                     para_env=force_env%sub_force_env(iforce_eval)%force_env%para_env)
            CALL cp_fm_create(matrix_s_tmp%matrix, fm_struct_tmp, name="s_matrix")
            CALL cp_fm_struct_release(fm_struct_tmp)
            CALL dbcsr_desymmetrize(dft_control%qs_control%cdft_control%matrix_s%matrix, desymm_tmp)
            CALL copy_dbcsr_to_fm(desymm_tmp, matrix_s_tmp%matrix)
            CALL dbcsr_release(desymm_tmp)
         END IF
         CALL dbcsr_release_p(dft_control%qs_control%cdft_control%matrix_s%matrix)
         ! Density_matrix (dbcsr -> fm)
         IF (mixed_cdft%calculate_metric) THEN
            DO ispin = 1, nspins
               ! Size AOxAO
               CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol_overlap, &
                                        ncol_global=ncol_overlap, context=blacs_env, &
                                        para_env=force_env%sub_force_env(iforce_eval)%force_env%para_env)
               CALL cp_fm_create(matrix_p_tmp(iforce_eval, ispin)%matrix, fm_struct_tmp, name="dm_matrix")
               CALL cp_fm_struct_release(fm_struct_tmp)
               CALL dbcsr_desymmetrize(dft_control%qs_control%cdft_control%matrix_p(ispin)%matrix, desymm_tmp)
               CALL copy_dbcsr_to_fm(desymm_tmp, matrix_p_tmp(iforce_eval, ispin)%matrix)
               CALL dbcsr_release(desymm_tmp)
               CALL dbcsr_release_p(dft_control%qs_control%cdft_control%matrix_p(ispin)%matrix)
            END DO
            DEALLOCATE (dft_control%qs_control%cdft_control%matrix_p)
         END IF
         ! Occupation numbers
         IF (.NOT. uniform_occupation) THEN
            DO ispin = 1, nspins
               IF (ncol_mo(ispin) /= SIZE(dft_control%qs_control%cdft_control%occupations(ispin)%array)) &
                  CPABORT("Array dimensions dont match.")
               IF (force_env_qs%para_env%ionode) THEN
                  ALLOCATE (occno_tmp(iforce_eval, ispin)%array(ncol_mo(ispin)))
                  occno_tmp(iforce_eval, ispin)%array = dft_control%qs_control%cdft_control%occupations(ispin)%array
               END IF
               DEALLOCATE (dft_control%qs_control%cdft_control%occupations(ispin)%array)
            END DO
            DEALLOCATE (dft_control%qs_control%cdft_control%occupations)
         END IF
      END DO
      ! Create needed fm structs
      CALL cp_fm_struct_create(fm_struct_wmat, nrow_global=nrow_wmat, ncol_global=ncol_wmat, &
                               context=mixed_cdft%blacs_env, para_env=force_env%para_env)
      CALL cp_fm_struct_create(fm_struct_overlap, nrow_global=nrow_overlap, ncol_global=ncol_overlap, &
                               context=mixed_cdft%blacs_env, para_env=force_env%para_env)
      ! Redistribute arrays with copy_general (this is not optimal for dbcsr matrices but...)
      ! We use this method for the serial case (mixed_cdft%run_type == mixed_cdft_serial) as well to move the arrays to the
      ! correct blacs_env, which is impossible using a simple copy of the arrays
      ALLOCATE (mixed_wmat_tmp(nforce_eval, nvar), mixed_matrix_s_tmp)
      IF (mixed_cdft%calculate_metric) &
         ALLOCATE (mixed_matrix_p_tmp(nforce_eval, nspins))
      DO iforce_eval = 1, nforce_eval
         ! MO coefficients
         DO ispin = 1, nspins
            NULLIFY (fm_struct_mo)
            CALL cp_fm_struct_create(fm_struct_mo, nrow_global=nrow_mo(ispin), ncol_global=ncol_mo(ispin), &
                                     context=mixed_cdft%blacs_env, para_env=force_env%para_env)
            NULLIFY (mixed_mo_coeff(iforce_eval, ispin)%matrix)
            CALL cp_fm_create(matrix=mixed_mo_coeff(iforce_eval, ispin)%matrix, &
                              matrix_struct=fm_struct_mo, &
                              name="MO_COEFF_"//TRIM(ADJUSTL(cp_to_string(iforce_eval)))//"_" &
                              //TRIM(ADJUSTL(cp_to_string(ispin)))//"_MATRIX")
            CALL cp_fm_copy_general(mo_coeff_tmp(iforce_eval, ispin)%matrix, &
                                    mixed_mo_coeff(iforce_eval, ispin)%matrix, &
                                    mixed_cdft%blacs_env%para_env)
            CALL cp_fm_release(mo_coeff_tmp(iforce_eval, ispin)%matrix)
            CALL cp_fm_struct_release(fm_struct_mo)
         END DO
         ! Weight
         DO ivar = 1, nvar
            NULLIFY (mixed_wmat_tmp(iforce_eval, ivar)%matrix)
            NULLIFY (w_matrix(iforce_eval, ivar)%matrix)
            CALL dbcsr_init_p(w_matrix(iforce_eval, ivar)%matrix)
            CALL cp_fm_create(matrix=mixed_wmat_tmp(iforce_eval, ivar)%matrix, &
                              matrix_struct=fm_struct_wmat, &
                              name="WEIGHT_"//TRIM(ADJUSTL(cp_to_string(iforce_eval)))//"_MATRIX")
            CALL cp_fm_copy_general(wmat_tmp(iforce_eval, ivar)%matrix, &
                                    mixed_wmat_tmp(iforce_eval, ivar)%matrix, &
                                    mixed_cdft%blacs_env%para_env)
            CALL cp_fm_release(wmat_tmp(iforce_eval, ivar)%matrix)
            ! (fm -> dbcsr)
            CALL copy_fm_to_dbcsr_bc(mixed_wmat_tmp(iforce_eval, ivar)%matrix, &
                                     w_matrix(iforce_eval, ivar)%matrix)
            CALL cp_fm_release(mixed_wmat_tmp(iforce_eval, ivar)%matrix)
         END DO
         ! Density matrix (fm -> dbcsr)
         IF (mixed_cdft%calculate_metric) THEN
            DO ispin = 1, nspins
               NULLIFY (mixed_matrix_p_tmp(iforce_eval, ispin)%matrix)
               NULLIFY (density_matrix(iforce_eval, ispin)%matrix)
               CALL dbcsr_init_p(density_matrix(iforce_eval, ispin)%matrix)
               CALL cp_fm_create(matrix=mixed_matrix_p_tmp(iforce_eval, ispin)%matrix, &
                                 matrix_struct=fm_struct_overlap, &
                                 name="DENSITY_"//TRIM(ADJUSTL(cp_to_string(iforce_eval)))//"_"// &
                                 TRIM(ADJUSTL(cp_to_string(ispin)))//"_MATRIX")
               CALL cp_fm_copy_general(matrix_p_tmp(iforce_eval, ispin)%matrix, &
                                       mixed_matrix_p_tmp(iforce_eval, ispin)%matrix, &
                                       mixed_cdft%blacs_env%para_env)
               CALL cp_fm_release(matrix_p_tmp(iforce_eval, ispin)%matrix)
               CALL copy_fm_to_dbcsr_bc(mixed_matrix_p_tmp(iforce_eval, ispin)%matrix, &
                                        density_matrix(iforce_eval, ispin)%matrix)
               CALL cp_fm_release(mixed_matrix_p_tmp(iforce_eval, ispin)%matrix)
            END DO
         END IF
      END DO
      CALL cp_fm_struct_release(fm_struct_wmat)
      DEALLOCATE (mo_coeff_tmp, wmat_tmp, mixed_wmat_tmp)
      IF (mixed_cdft%calculate_metric) THEN
         DEALLOCATE (matrix_p_tmp)
         DEALLOCATE (mixed_matrix_p_tmp)
      END IF
      ! Overlap (fm -> dbcsr)
      NULLIFY (mixed_matrix_s_tmp%matrix)
      CALL cp_fm_create(matrix=mixed_matrix_s_tmp%matrix, &
                        matrix_struct=fm_struct_overlap, &
                        name="OVERLAP_MATRIX")
      CALL cp_fm_struct_release(fm_struct_overlap)
      CALL cp_fm_copy_general(matrix_s_tmp%matrix, &
                              mixed_matrix_s_tmp%matrix, &
                              mixed_cdft%blacs_env%para_env)
      CALL cp_fm_release(matrix_s_tmp%matrix)
      CALL copy_fm_to_dbcsr_bc(mixed_matrix_s_tmp%matrix, mixed_matrix_s)
      CALL cp_fm_release(mixed_matrix_s_tmp%matrix)
      ! Occupation numbers
      IF (.NOT. uniform_occupation) THEN
         DO iforce_eval = 1, nforce_eval
            DO ispin = 1, nspins
               ALLOCATE (mixed_cdft%occupations(iforce_eval, ispin)%array(ncol_mo(ispin)))
               mixed_cdft%occupations(iforce_eval, ispin)%array = 0.0_dp
               IF (ASSOCIATED(occno_tmp(iforce_eval, ispin)%array)) THEN
                  mixed_cdft%occupations(iforce_eval, ispin)%array = occno_tmp(iforce_eval, ispin)%array
                  DEALLOCATE (occno_tmp(iforce_eval, ispin)%array)
               END IF
               CALL mp_sum(mixed_cdft%occupations(iforce_eval, ispin)%array, force_env%para_env%group)
            END DO
         END DO
         DEALLOCATE (occno_tmp)
      END IF
      DEALLOCATE (matrix_s_tmp, mixed_matrix_s_tmp)
      DEALLOCATE (ncol_mo, nrow_mo)

   END SUBROUTINE mixed_cdft_redistribute_arrays
! **************************************************************************************************
!> \brief Routine to print out the electronic coupling(s) between CDFT states.
!> \param force_env the force_env that holds the CDFT states
!> \par History
!>       11.17  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_print_couplings(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      INTEGER                                            :: iounit, ipermutation, istate, ivar, &
                                                            jstate, nforce_eval, npermutations, &
                                                            nvar
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(section_vals_type), POINTER                   :: force_env_section, print_section

      NULLIFY (print_section, mixed_cdft)

      logger => cp_get_default_logger()
      CPASSERT(ASSOCIATED(force_env))
      CALL force_env_get(force_env=force_env, &
                         force_env_section=force_env_section)
      CALL get_mixed_env(force_env%mixed_env, cdft_control=mixed_cdft)
      CPASSERT(ASSOCIATED(mixed_cdft))
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
      iounit = cp_print_key_unit_nr(logger, print_section, '', extension='.mixedLog')
      !
      CPASSERT(ALLOCATED(mixed_cdft%results%strength))
      CPASSERT(ALLOCATED(mixed_cdft%results%W_diagonal))
      CPASSERT(ALLOCATED(mixed_cdft%results%S))
      CPASSERT(ALLOCATED(mixed_cdft%results%energy))
      nforce_eval = SIZE(force_env%sub_force_env)
      nvar = SIZE(mixed_cdft%results%strength, 1)
      npermutations = nforce_eval*(nforce_eval - 1)/2 ! Size of upper triangular part
      IF (iounit > 0) THEN
         WRITE (iounit, '(/,T3,A,T66)') &
            '------------------------- CDFT coupling information --------------------------'
         WRITE (iounit, '(T3,A,T66,(3X,F12.2))') &
            'Information at step (fs):', mixed_cdft%sim_step*mixed_cdft%sim_dt
         DO ipermutation = 1, npermutations
            CALL map_permutation_to_states(nforce_eval, ipermutation, istate, jstate)
            WRITE (iounit, '(/,T3,A)') REPEAT('#', 44)
            WRITE (iounit, '(T3,A,I3,A,I3,A)') '###### CDFT states I =', istate, ' and J = ', jstate, ' ######'
            WRITE (iounit, '(T3,A)') REPEAT('#', 44)
            DO ivar = 1, nvar
               IF (ivar > 1) &
                  WRITE (iounit, '(A)') ''
               WRITE (iounit, '(T3,A,T60,(3X,I18))') 'Atomic group:', ivar
               WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
                  'Strength of constraint I:', mixed_cdft%results%strength(ivar, istate)
               WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
                  'Strength of constraint J:', mixed_cdft%results%strength(ivar, jstate)
               WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
                  'Final value of constraint I:', mixed_cdft%results%W_diagonal(ivar, istate)
               WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
                  'Final value of constraint J:', mixed_cdft%results%W_diagonal(ivar, jstate)
            END DO
            WRITE (iounit, '(/,T3,A,T60,(3X,F18.12))') &
               'Overlap between states I and J:', mixed_cdft%results%S(istate, jstate)
            WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
               'Charge transfer energy (J-I) (Hartree):', (mixed_cdft%results%energy(jstate) - mixed_cdft%results%energy(istate))
            WRITE (iounit, *)
            IF (ALLOCATED(mixed_cdft%results%rotation)) THEN
               IF (ABS(mixed_cdft%results%rotation(ipermutation))*1.0E3_dp .GE. 0.1_dp) THEN
                  WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
                     'Diabatic electronic coupling (rotation, mHartree):', &
                     ABS(mixed_cdft%results%rotation(ipermutation)*1.0E3_dp)
               ELSE
                  WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
                     'Diabatic electronic coupling (rotation, microHartree):', &
                     ABS(mixed_cdft%results%rotation(ipermutation)*1.0E6_dp)
               END IF
            END IF
            IF (ALLOCATED(mixed_cdft%results%lowdin)) THEN
               IF (ABS(mixed_cdft%results%lowdin(ipermutation))*1.0E3_dp .GE. 0.1_dp) THEN
                  WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
                     'Diabatic electronic coupling (Lowdin, mHartree):', &
                     ABS(mixed_cdft%results%lowdin(ipermutation)*1.0E3_dp)
               ELSE
                  WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
                     'Diabatic electronic coupling (Lowdin, microHartree):', &
                     ABS(mixed_cdft%results%lowdin(ipermutation)*1.0E6_dp)
               END IF
            END IF
            IF (ALLOCATED(mixed_cdft%results%wfn)) THEN
               IF (mixed_cdft%results%wfn(ipermutation)*1.0E3_dp .GE. 0.1_dp) THEN
                  WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
                     'Diabatic electronic coupling (wfn overlap, mHartree):', &
                     ABS(mixed_cdft%results%wfn(ipermutation)*1.0E3_dp)
               ELSE
                  WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
                     'Diabatic electronic coupling (wfn overlap, microHartree):', &
                     ABS(mixed_cdft%results%wfn(ipermutation)*1.0E6_dp)
               END IF
            END IF
            IF (ALLOCATED(mixed_cdft%results%nonortho)) THEN
               WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
                  'Diabatic electronic coupling (nonorthogonal, Hartree):', mixed_cdft%results%nonortho(ipermutation)
            END IF
            IF (ALLOCATED(mixed_cdft%results%metric)) THEN
               WRITE (iounit, *)
               IF (SIZE(mixed_cdft%results%metric, 2) == 1) THEN
                  WRITE (iounit, '(T3,A,T66,(3X,F12.6))') &
                     'Coupling reliability metric (0 is ideal):', mixed_cdft%results%metric(ipermutation, 1)
               ELSE
                  WRITE (iounit, '(T3,A,T54,(3X,2F12.6))') &
                     'Coupling reliability metric (0 is ideal):', &
                     mixed_cdft%results%metric(ipermutation, 1), mixed_cdft%results%metric(ipermutation, 2)
               END IF
            END IF
         END DO
         WRITE (iounit, '(T3,A)') &
            '------------------------------------------------------------------------------'
      END IF
      CALL cp_print_key_finished_output(iounit, logger, force_env_section, &
                                        "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")

   END SUBROUTINE mixed_cdft_print_couplings

! **************************************************************************************************
!> \brief Release storage reserved for mixed CDFT matrices
!> \param force_env the force_env that holds the CDFT states
!> \par History
!>       11.17  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_release_work(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft

      NULLIFY (mixed_cdft)

      CPASSERT(ASSOCIATED(force_env))
      CALL get_mixed_env(force_env%mixed_env, cdft_control=mixed_cdft)
      CPASSERT(ASSOCIATED(mixed_cdft))
      CALL mixed_cdft_result_type_release(mixed_cdft%results)

   END SUBROUTINE mixed_cdft_release_work

! **************************************************************************************************
!> \brief Given the size of a symmetric matrix and a permutation index, returns indices (i, j) of the
!>        off-diagonal element that corresponds to the permutation index. Assumes that the permutation
!>        index was computed by going through the upper triangular part of the input matrix row-by-row.
!> \param n the size of the symmetric matrix
!> \param ipermutation the permutation index
!> \param i the row index corresponding to ipermutation
!> \param j the column index corresponding to ipermutation
! **************************************************************************************************
   SUBROUTINE map_permutation_to_states(n, ipermutation, i, j)
      INTEGER, INTENT(IN)                                :: n, ipermutation
      INTEGER, INTENT(OUT)                               :: i, j

      INTEGER                                            :: kcol, kpermutation, krow, npermutations

      npermutations = n*(n - 1)/2 ! Size of upper triangular part
      IF (ipermutation > npermutations) &
         CPABORT("Permutation index out of bounds")
      kpermutation = 0
      DO krow = 1, n
         DO kcol = krow + 1, n
            kpermutation = kpermutation + 1
            IF (kpermutation == ipermutation) THEN
               i = krow
               j = kcol
               RETURN
            END IF
         END DO
      END DO

   END SUBROUTINE map_permutation_to_states

! **************************************************************************************************
!> \brief Determine confinement bounds along confinement dir (hardcoded to be z)
!>        and determine the number of nonzero entries
!>        Optionally zero entries below a given threshold
!> \param fun input 3D potential (real space)
!> \param th threshold for screening values
!> \param just_zero determines if fun should only be zeroed without returning bounds/work
!> \param bounds the confinement bounds: fun is nonzero only between these values along 3rd dimension
!> \param work an estimate of the total number of grid points where fun is nonzero
! **************************************************************************************************
   SUBROUTINE hfun_zero(fun, th, just_zero, bounds, work)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: fun
      REAL(KIND=dp), INTENT(IN)                          :: th
      LOGICAL                                            :: just_zero
      INTEGER, OPTIONAL                                  :: bounds(2), work

      INTEGER                                            :: i1, i2, i3, lb, n1, n2, n3, nzeroed, &
                                                            nzeroed_total, ub
      LOGICAL                                            :: lb_final, ub_final

      n1 = SIZE(fun, 1)
      n2 = SIZE(fun, 2)
      n3 = SIZE(fun, 3)
      nzeroed_total = 0
      IF (.NOT. just_zero) THEN
         CPASSERT(PRESENT(bounds))
         CPASSERT(PRESENT(work))
         lb = 1
         lb_final = .FALSE.
         ub_final = .FALSE.
      END IF
      DO i3 = 1, n3
         IF (.NOT. just_zero) nzeroed = 0
         DO i2 = 1, n2
            DO i1 = 1, n1
               IF (fun(i1, i2, i3) < th) THEN
                  IF (.NOT. just_zero) THEN
                     nzeroed = nzeroed + 1
                     nzeroed_total = nzeroed_total + 1
                  ELSE
                     fun(i1, i2, i3) = 0.0_dp
                  END IF
               END IF
            END DO
         END DO
         IF (.NOT. just_zero) THEN
            IF (nzeroed == (n2*n1)) THEN
               IF (.NOT. lb_final) THEN
                  lb = i3
               ELSE IF (.NOT. ub_final) THEN
                  ub = i3
                  ub_final = .TRUE.
               END IF
            ELSE
               IF (.NOT. lb_final) lb_final = .TRUE.
               IF (ub_final) ub_final = .FALSE. ! Safeguard against "holes"
            END IF
         END IF
      END DO
      IF (.NOT. just_zero) THEN
         IF (.NOT. ub_final) ub = n3
         bounds(1) = lb
         bounds(2) = ub
         bounds = bounds - (n3/2) - 1
         work = n3*n2*n1 - nzeroed_total
      END IF

   END SUBROUTINE hfun_zero

! **************************************************************************************************
!> \brief Read input section related to block diagonalization of the mixed CDFT Hamiltonian matrix.
!> \param force_env        the force_env that holds the CDFT states
!> \param blocks           list of CDFT states defining the matrix blocks
!> \param ignore_excited   flag that determines if excited states resulting from the block
!>                         diagonalization process should be ignored
!> \param nrecursion       integer that determines how many steps of recursive block diagonalization
!>                         is performed (1 if disabled)
!> \par History
!>       01.18  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_read_block_diag(force_env, blocks, ignore_excited, nrecursion)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(cp_1d_i_p_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: blocks
      LOGICAL, INTENT(OUT)                               :: ignore_excited
      INTEGER, INTENT(OUT)                               :: nrecursion

      INTEGER                                            :: i, j, k, l, nblk, nforce_eval
      INTEGER, DIMENSION(:), POINTER                     :: tmplist
      LOGICAL                                            :: do_recursive, explicit, has_duplicates
      TYPE(section_vals_type), POINTER                   :: block_section, force_env_section

      EXTERNAL                                           :: dsygv

      NULLIFY (force_env_section, block_section)
      CPASSERT(ASSOCIATED(force_env))
      nforce_eval = SIZE(force_env%sub_force_env)

      CALL force_env_get(force_env=force_env, &
                         force_env_section=force_env_section)
      block_section => section_vals_get_subs_vals(force_env_section, "MIXED%MIXED_CDFT%BLOCK_DIAGONALIZE")

      CALL section_vals_get(block_section, explicit=explicit)
      IF (.NOT. explicit) &
         CALL cp_abort(__LOCATION__, &
                       "Block diagonalization of CDFT Hamiltonian was requested, but the "// &
                       "corresponding input section is missing!")

      CALL section_vals_val_get(block_section, "BLOCK", n_rep_val=nblk)
      ALLOCATE (blocks(nblk))
      DO i = 1, nblk
         NULLIFY (blocks(i)%array)
         CALL section_vals_val_get(block_section, "BLOCK", i_rep_val=i, i_vals=tmplist)
         IF (SIZE(tmplist) < 1) &
            CPABORT("Each BLOCK must contain at least 1 state.")
         ALLOCATE (blocks(i)%array(SIZE(tmplist)))
         blocks(i)%array(:) = tmplist(:)
      END DO
      CALL section_vals_val_get(block_section, "IGNORE_EXCITED", l_val=ignore_excited)
      CALL section_vals_val_get(block_section, "RECURSIVE_DIAGONALIZATION", l_val=do_recursive)
      ! Check that the requested states exist
      DO i = 1, nblk
         DO j = 1, SIZE(blocks(i)%array)
            IF (blocks(i)%array(j) < 1 .OR. blocks(i)%array(j) > nforce_eval) &
               CPABORT("Requested state does not exist.")
         END DO
      END DO
      ! Check for duplicates
      has_duplicates = .FALSE.
      DO i = 1, nblk
         ! Within same block
         DO j = 1, SIZE(blocks(i)%array)
            DO k = j + 1, SIZE(blocks(i)%array)
               IF (blocks(i)%array(j) == blocks(i)%array(k)) has_duplicates = .TRUE.
            END DO
         END DO
         ! Within different blocks
         DO j = i + 1, nblk
            DO k = 1, SIZE(blocks(i)%array)
               DO l = 1, SIZE(blocks(j)%array)
                  IF (blocks(i)%array(k) == blocks(j)%array(l)) has_duplicates = .TRUE.
               END DO
            END DO
         END DO
      END DO
      IF (has_duplicates) CPABORT("Duplicate states are not allowed.")
      nrecursion = 1
      IF (do_recursive) THEN
         IF (MODULO(nblk, 2) /= 0) THEN
            CALL cp_warn(__LOCATION__, &
                         "Number of blocks not divisible with 2. Recursive diagonalization not possible. "// &
                         "Calculation proceeds without.")
            nrecursion = 1
         ELSE
            nrecursion = nblk/2
         END IF
         IF (nrecursion /= 1 .AND. .NOT. ignore_excited) &
            CALL cp_abort(__LOCATION__, &
                          "Keyword IGNORE_EXCITED must be active for recursive diagonalization.")
      END IF

   END SUBROUTINE mixed_cdft_read_block_diag

! **************************************************************************************************
!> \brief Assembles the matrix blocks from the mixed CDFT Hamiltonian.
!> \param mixed_cdft the env that holds the CDFT states
!> \param blocks  list of CDFT states defining the matrix blocks
!> \param H_block list of Hamiltonian matrix blocks
!> \param S_block list of overlap matrix blocks
!> \par History
!>       01.18  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_get_blocks(mixed_cdft, blocks, H_block, S_block)
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(cp_1d_i_p_type), ALLOCATABLE, DIMENSION(:)    :: blocks
      TYPE(cp_2d_r_p_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: H_block, S_block

      INTEGER                                            :: i, icol, irow, j, k, nblk

      EXTERNAL                                           :: dsygv

      CPASSERT(ASSOCIATED(mixed_cdft))

      nblk = SIZE(blocks)
      ALLOCATE (H_block(nblk), S_block(nblk))
      DO i = 1, nblk
         NULLIFY (H_block(i)%array)
         NULLIFY (S_block(i)%array)
         ALLOCATE (H_block(i)%array(SIZE(blocks(i)%array), SIZE(blocks(i)%array)))
         ALLOCATE (S_block(i)%array(SIZE(blocks(i)%array), SIZE(blocks(i)%array)))
         icol = 0
         DO j = 1, SIZE(blocks(i)%array)
            irow = 0
            icol = icol + 1
            DO k = 1, SIZE(blocks(i)%array)
               irow = irow + 1
               H_block(i)%array(irow, icol) = mixed_cdft%results%H(blocks(i)%array(k), blocks(i)%array(j))
               S_block(i)%array(irow, icol) = mixed_cdft%results%S(blocks(i)%array(k), blocks(i)%array(j))
            END DO
         END DO
         ! Check that none of the interaction energies is repulsive
         IF (ANY(H_block(i)%array .GE. 0.0_dp)) &
            CALL cp_abort(__LOCATION__, &
                          "At least one of the interaction energies within block "//TRIM(ADJUSTL(cp_to_string(i)))// &
                          " is repulsive.")
      END DO

   END SUBROUTINE mixed_cdft_get_blocks

! **************************************************************************************************
!> \brief Diagonalizes each of the matrix blocks.
!> \param blocks  list of CDFT states defining the matrix blocks
!> \param H_block list of Hamiltonian matrix blocks
!> \param S_block list of overlap matrix blocks
!> \param eigenvalues list of eigenvalues for each block
!> \par History
!>       01.18  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_diagonalize_blocks(blocks, H_block, S_block, eigenvalues)
      TYPE(cp_1d_i_p_type), ALLOCATABLE, DIMENSION(:)    :: blocks
      TYPE(cp_2d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: H_block, S_block
      TYPE(cp_1d_r_p_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: eigenvalues

      INTEGER                                            :: i, info, nblk, work_array_size
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: H_mat_copy, S_mat_copy

      EXTERNAL                                           :: dsygv

      nblk = SIZE(blocks)
      ALLOCATE (eigenvalues(nblk))
      DO i = 1, nblk
         NULLIFY (eigenvalues(i)%array)
         ALLOCATE (eigenvalues(i)%array(SIZE(blocks(i)%array)))
         eigenvalues(i)%array = 0.0_dp
         ! Workspace query
         ALLOCATE (work(1))
         info = 0
         ALLOCATE (H_mat_copy(SIZE(blocks(i)%array), SIZE(blocks(i)%array)))
         ALLOCATE (S_mat_copy(SIZE(blocks(i)%array), SIZE(blocks(i)%array)))
         H_mat_copy(:, :) = H_block(i)%array(:, :) ! Need explicit copies because dsygv destroys original values
         S_mat_copy(:, :) = S_block(i)%array(:, :)
         CALL dsygv(1, 'V', 'U', SIZE(blocks(i)%array), H_mat_copy, SIZE(blocks(i)%array), &
                    S_mat_copy, SIZE(blocks(i)%array), eigenvalues(i)%array, work, -1, info)
         work_array_size = NINT(work(1))
         DEALLOCATE (H_mat_copy, S_mat_copy)
         ! Allocate work array
         DEALLOCATE (work)
         ALLOCATE (work(work_array_size))
         work = 0.0_dp
         ! Solve Hc = eSc
         info = 0
         CALL dsygv(1, 'V', 'U', SIZE(blocks(i)%array), H_block(i)%array, SIZE(blocks(i)%array), &
                    S_block(i)%array, SIZE(blocks(i)%array), eigenvalues(i)%array, work, work_array_size, info)
         IF (info /= 0) THEN
            IF (info > SIZE(blocks(i)%array)) THEN
               CPABORT("Matrix S is not positive definite")
            ELSE
               CPABORT("Diagonalization of H matrix failed.")
            END IF
         END IF
         DEALLOCATE (work)
      END DO

   END SUBROUTINE mixed_cdft_diagonalize_blocks

! **************************************************************************************************
!> \brief Assembles the new block diagonalized mixed CDFT Hamiltonian and overlap matrices.
!> \param mixed_cdft the env that holds the CDFT states
!> \param blocks  list of CDFT states defining the matrix blocks
!> \param H_block list of Hamiltonian matrix blocks
!> \param eigenvalues list of eigenvalues for each block
!> \param n size of the new Hamiltonian and overlap matrices
!> \param iounit the output unit
!> \par History
!>       01.18  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_assemble_block_diag(mixed_cdft, blocks, H_block, eigenvalues, &
                                             n, iounit)
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(cp_1d_i_p_type), ALLOCATABLE, DIMENSION(:)    :: blocks
      TYPE(cp_2d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: H_block
      TYPE(cp_1d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: eigenvalues
      INTEGER                                            :: n, iounit

      CHARACTER(LEN=20)                                  :: ilabel, jlabel
      CHARACTER(LEN=3)                                   :: tmp
      INTEGER                                            :: i, icol, ipermutation, irow, j, k, l, &
                                                            nblk, npermutations
      LOGICAL                                            :: ignore_excited
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: H_mat, H_offdiag, S_mat, S_offdiag

      EXTERNAL                                           :: dsygv

      ALLOCATE (H_mat(n, n), S_mat(n, n))
      nblk = SIZE(blocks)
      ignore_excited = (nblk == n)
      ! The diagonal contains the eigenvalues of each block
      IF (iounit > 0) WRITE (iounit, '(/,T3,A)') "Eigenvalues of the block diagonalized states"
      H_mat(:, :) = 0.0_dp
      S_mat(:, :) = 0.0_dp
      k = 1
      DO i = 1, nblk
         IF (iounit > 0) WRITE (iounit, '(T6,A,I3)') "Block", i
         DO j = 1, SIZE(eigenvalues(i)%array)
            H_mat(k, k) = eigenvalues(i)%array(j)
            S_mat(k, k) = 1.0_dp
            k = k + 1
            IF (iounit > 0) THEN
               IF (j == 1) THEN
                  WRITE (iounit, '(T9,A,T58,(3X,F20.14))') 'Ground state energy:', eigenvalues(i)%array(j)
               ELSE
                  WRITE (iounit, '(T9,A,I2,A,T58,(3X,F20.14))') &
                     'Excited state (', j - 1, ' ) energy:', eigenvalues(i)%array(j)
               END IF
            END IF
            IF (ignore_excited .AND. j == 1) EXIT
         END DO
      END DO
      ! Transform the off-diagonal blocks using the eigenvectors of each block
      npermutations = nblk*(nblk - 1)/2
      IF (iounit > 0) WRITE (iounit, '(/,T3,A)') "Interactions between block diagonalized states"
      DO ipermutation = 1, npermutations
         CALL map_permutation_to_states(nblk, ipermutation, i, j)
         ! Get the untransformed off-diagonal block
         ALLOCATE (H_offdiag(SIZE(blocks(i)%array), SIZE(blocks(j)%array)))
         ALLOCATE (S_offdiag(SIZE(blocks(i)%array), SIZE(blocks(j)%array)))
         icol = 0
         DO k = 1, SIZE(blocks(j)%array)
            irow = 0
            icol = icol + 1
            DO l = 1, SIZE(blocks(i)%array)
               irow = irow + 1
               H_offdiag(irow, icol) = mixed_cdft%results%H(blocks(i)%array(l), blocks(j)%array(k))
               S_offdiag(irow, icol) = mixed_cdft%results%S(blocks(i)%array(l), blocks(j)%array(k))
            END DO
         END DO
         ! Check that none of the interaction energies is repulsive
         IF (ANY(H_offdiag .GE. 0.0_dp)) &
            CALL cp_abort(__LOCATION__, &
                          "At least one of the interaction energies between blocks "//TRIM(ADJUSTL(cp_to_string(i)))// &
                          " and "//TRIM(ADJUSTL(cp_to_string(j)))//" is repulsive.")
         ! Now transform: C_i^T * H * C_j
         H_offdiag(:, :) = MATMUL(H_offdiag, H_block(j)%array)
         H_offdiag(:, :) = MATMUL(TRANSPOSE(H_block(i)%array), H_offdiag)
         S_offdiag(:, :) = MATMUL(S_offdiag, H_block(j)%array)
         S_offdiag(:, :) = MATMUL(TRANSPOSE(H_block(i)%array), S_offdiag)
         ! Make sure the transformation preserves the sign of elements in the S and H matrices
         ! The S/H matrices contain only positive/negative values so that any sign flipping occurs in the
         ! same elements in both matrices
         ! Check for sign flipping using the S matrix
         IF (ANY(S_offdiag .LT. 0.0_dp)) THEN
            DO l = 1, SIZE(S_offdiag, 2)
               DO k = 1, SIZE(S_offdiag, 1)
                  IF (S_offdiag(k, l) .LT. 0.0_dp) THEN
                     S_offdiag(k, l) = -1.0_dp*S_offdiag(k, l)
                     H_offdiag(k, l) = -1.0_dp*H_offdiag(k, l)
                  END IF
               END DO
            END DO
         END IF
         IF (ignore_excited) THEN
            H_mat(i, j) = H_offdiag(1, 1)
            H_mat(j, i) = H_mat(i, j)
            S_mat(i, j) = S_offdiag(1, 1)
            S_mat(j, i) = S_mat(i, j)
         ELSE
            irow = 1
            icol = 1
            DO k = 1, i - 1
               irow = irow + SIZE(blocks(k)%array)
            END DO
            DO k = 1, j - 1
               icol = icol + SIZE(blocks(k)%array)
            END DO
            H_mat(irow:irow + SIZE(H_offdiag, 1) - 1, icol:icol + SIZE(H_offdiag, 2) - 1) = H_offdiag(:, :)
            H_mat(icol:icol + SIZE(H_offdiag, 2) - 1, irow:irow + SIZE(H_offdiag, 1) - 1) = TRANSPOSE(H_offdiag)
            S_mat(irow:irow + SIZE(H_offdiag, 1) - 1, icol:icol + SIZE(H_offdiag, 2) - 1) = S_offdiag(:, :)
            S_mat(icol:icol + SIZE(H_offdiag, 2) - 1, irow:irow + SIZE(H_offdiag, 1) - 1) = TRANSPOSE(S_offdiag)
         END IF
         IF (iounit > 0) THEN
            WRITE (iounit, '(/,T3,A)') REPEAT('#', 39)
            WRITE (iounit, '(T3,A,I3,A,I3,A)') '###### Blocks I =', i, ' and J = ', j, ' ######'
            WRITE (iounit, '(T3,A)') REPEAT('#', 39)
            WRITE (iounit, '(T3,A)') 'Interaction energies'
            DO irow = 1, SIZE(H_offdiag, 1)
               ilabel = "(ground state)"
               IF (irow > 1) THEN
                  IF (ignore_excited) EXIT
                  WRITE (tmp, '(I3)') irow - 1
                  ilabel = "(excited state "//TRIM(ADJUSTL(tmp))//")"
               END IF
               DO icol = 1, SIZE(H_offdiag, 2)
                  jlabel = "(ground state)"
                  IF (icol > 1) THEN
                     IF (ignore_excited) EXIT
                     WRITE (tmp, '(I3)') icol - 1
                     jlabel = "(excited state "//TRIM(ADJUSTL(tmp))//")"
                  END IF
                  WRITE (iounit, '(T6,A,T58,(3X,F20.14))') TRIM(ilabel)//'-'//TRIM(jlabel)//':', H_offdiag(irow, icol)
               END DO
            END DO
            WRITE (iounit, '(T3,A)') 'Overlaps'
            DO irow = 1, SIZE(H_offdiag, 1)
               ilabel = "(ground state)"
               IF (irow > 1) THEN
                  IF (ignore_excited) EXIT
                  ilabel = "(excited state)"
                  WRITE (tmp, '(I3)') irow - 1
                  ilabel = "(excited state "//TRIM(ADJUSTL(tmp))//")"
               END IF
               DO icol = 1, SIZE(H_offdiag, 2)
                  jlabel = "(ground state)"
                  IF (icol > 1) THEN
                     IF (ignore_excited) EXIT
                     WRITE (tmp, '(I3)') icol - 1
                     jlabel = "(excited state "//TRIM(ADJUSTL(tmp))//")"
                  END IF
                  WRITE (iounit, '(T6,A,T58,(3X,F20.14))') TRIM(ilabel)//'-'//TRIM(jlabel)//':', S_offdiag(irow, icol)
               END DO
            END DO
         END IF
         DEALLOCATE (H_offdiag, S_offdiag)
      END DO
      CALL mixed_cdft_result_type_set(mixed_cdft%results, H=H_mat, S=S_mat)
      ! Deallocate work
      DEALLOCATE (H_mat, S_mat)

   END SUBROUTINE mixed_cdft_assemble_block_diag

END MODULE mixed_cdft_utils
