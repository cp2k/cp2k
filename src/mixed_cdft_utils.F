!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
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
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_retain,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: becke_control_create,&
                                              dft_control_type
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
   USE input_constants,                 ONLY: becke_cavity_conf,&
                                              becke_cutoff_element,&
                                              becke_dynamic_conf,&
                                              becke_mixed_conf,&
                                              becke_none_conf,&
                                              becke_static_conf,&
                                              cdft_density_constraint,&
                                              shape_function_gaussian
   USE input_section_types,             ONLY: section_vals_duplicate,&
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
   USE mixed_cdft_types,                ONLY: mixed_cdft_settings_type,&
                                              mixed_cdft_type
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
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

! *** Public subroutines ***

   PUBLIC :: mixed_cdft_parse_settings, mixed_cdft_transfer_settings, &
             mixed_cdft_init_structures, mixed_cdft_redistribute_arrays, &
             hfun_zero

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

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_parse_settings', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iforce_eval, nforce_eval, nkinds
      LOGICAL                                            :: is_match
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(force_env_type), POINTER                      :: force_env_qs
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_environment_type), POINTER                 :: qs_env

      NULLIFY (dft_control, qs_env, pw_env, auxbas_pw_pool, force_env_qs)
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
      ALLOCATE (settings%coeffs(natom, nforce_eval))
      ALLOCATE (settings%si(17, nforce_eval))
      ALLOCATE (settings%sr(6, nforce_eval))
      ALLOCATE (settings%confine_bounds(2, nforce_eval))
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
      settings%coeffs = 0.0_dp
      settings%si = 0
      settings%sr = 0.0_dp
      settings%confine_bounds = 0.0_dp
      settings%cutoffs = 0.0_dp
      settings%radii = 0.0_dp
      ! Get information from the sub_force_envs
      ! TODO: transfer logicals as logicals instead of ints
      DO iforce_eval = 1, nforce_eval
         IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
         IF (mixed_env%do_mixed_qmmm_cdft) THEN
            force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
            qs_env => force_env_qs%qmmm_env%qs_env
         ELSE
            force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
            CALL force_env_get(force_env_qs, qs_env=qs_env)
         END IF
         CALL get_qs_env(qs_env, pw_env=pw_env, dft_control=dft_control)
         IF (.NOT. dft_control%qs_control%becke_restraint) &
            CALL cp_abort(__LOCATION__, "Mixed CDFT currently only supports the Becke constraint."// &
                          "Becke constraint must be active in QS section of both force_evals!")
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         settings%bo = auxbas_pw_pool%pw_grid%bounds_local
         ! Only the rank 0 process collects info about pw_grid and becke
         IF (force_env_qs%para_env%mepos == &
             force_env_qs%para_env%source) THEN
            settings%grid_span(iforce_eval) = auxbas_pw_pool%pw_grid%grid_span
            settings%npts(:, iforce_eval) = auxbas_pw_pool%pw_grid%npts
            settings%cutoff(iforce_eval) = auxbas_pw_pool%pw_grid%cutoff
            settings%rel_cutoff(iforce_eval) = dft_control%qs_control%relative_cutoff
            IF (auxbas_pw_pool%pw_grid%spherical) settings%spherical(iforce_eval) = 1
            settings%rs_dims(:, iforce_eval) = auxbas_pw_pool%pw_grid%para%rs_dims
            IF (auxbas_pw_pool%pw_grid%grid_span == HALFSPACE) settings%odd(iforce_eval) = 1
            IF (dft_control%qs_control%becke_control%natoms .GT. SIZE(settings%atoms, 1)) &
               CALL cp_abort(__LOCATION__, "More Becke constraint atoms than defined in mixed section."// &
                             " Use default values for MIXED\MAPPING.")
            settings%atoms(1:dft_control%qs_control%becke_control%natoms, iforce_eval) = &
               dft_control%qs_control%becke_control%atoms
            settings%coeffs(1:dft_control%qs_control%becke_control%natoms, iforce_eval) = &
               dft_control%qs_control%becke_control%coeff
            ! Only GAPW and serial calculation needs to keep a copy of these arrays for constraint integration
            IF (.NOT. dft_control%qs_control%gapw) THEN
               IF (mixed_cdft%is_parallel .OR. &
                   (.NOT. mixed_cdft%is_parallel .AND. iforce_eval /= 1)) THEN
                  DEALLOCATE (dft_control%qs_control%becke_control%atoms)
                  DEALLOCATE (dft_control%qs_control%becke_control%coeff)
               END IF
            END IF
            settings%si(1, iforce_eval) = dft_control%qs_control%becke_control%confine_method
            settings%si(2, iforce_eval) = dft_control%qs_control%becke_control%cutoff_type
            IF (dft_control%qs_control%becke_control%confine) settings%si(3, iforce_eval) = 1
            IF (dft_control%qs_control%becke_control%dynamic_confine) settings%si(4, iforce_eval) = 1
            IF (dft_control%qs_control%becke_control%cavity_confine) THEN
               settings%si(5, iforce_eval) = 1
               IF (.NOT. mixed_cdft%is_parallel .AND. iforce_eval /= 1) &
                  CALL release_hirshfeld_type(dft_control%qs_control%becke_control%cavity_env)
            END IF
            IF (dft_control%qs_control%becke_control%should_skip) settings%si(6, iforce_eval) = 1
            IF (dft_control%qs_control%becke_control%print_cavity) settings%si(7, iforce_eval) = 1
            IF (dft_control%qs_control%becke_control%in_memory) settings%si(8, iforce_eval) = 1
            settings%si(9, iforce_eval) = dft_control%qs_control%becke_control%confine_dir
            IF (dft_control%qs_control%becke_control%adjust) settings%si(10, iforce_eval) = 1
            IF (dft_control%qs_control%becke_control%atomic_charges) settings%si(11, iforce_eval) = 1
            settings%si(12, iforce_eval) = dft_control%qs_control%becke_control%cavity_shape
            settings%si(13, iforce_eval) = dft_control%multiplicity
            IF (qs_env%has_unit_metric) THEN
               settings%si(14, iforce_eval) = 1
            ELSE
               settings%si(14, iforce_eval) = 0
            END IF
            settings%si(15, iforce_eval) = dft_control%qs_control%becke_control%constraint_type
            settings%si(16, iforce_eval) = dft_control%qs_control%becke_control%combined_type
            IF (dft_control%qs_control%becke_control%use_bohr) THEN
               settings%si(17, iforce_eval) = 1
            ELSE
               settings%si(17, iforce_eval) = 0
            END IF
            settings%sr(1, iforce_eval) = dft_control%qs_control%becke_control%dynamic_radius
            settings%sr(2, iforce_eval) = dft_control%qs_control%becke_control%rcavity
            settings%sr(3, iforce_eval) = dft_control%qs_control%becke_control%rglobal
            settings%sr(4, iforce_eval) = dft_control%qs_control%becke_control%eps_cavity
            settings%sr(5, iforce_eval) = dft_control%qs_control%eps_rho_rspace
            settings%sr(6, iforce_eval) = pw_env%cube_info(pw_env%auxbas_grid)%max_rad_ga
            settings%confine_bounds(:, iforce_eval) = dft_control%qs_control%becke_control%confine_bounds
            IF (dft_control%qs_control%becke_control%cutoff_type == becke_cutoff_element) THEN
               nkinds = SIZE(dft_control%qs_control%becke_control%cutoffs_tmp)
               IF (nkinds .GT. settings%max_nkinds) &
                  CALL cp_abort(__LOCATION__, "More than "//TRIM(cp_to_string(settings%max_nkinds))// &
                                " unique elements were defined in BECKE_RESTRAINT\ELEMENT_CUTOFF. Are you sure"// &
                                " your input is correct? If yes, please increase max_nkinds and recompile.")
               settings%cutoffs(1:nkinds, iforce_eval) = dft_control%qs_control%becke_control%cutoffs_tmp(:)
               IF (.NOT. mixed_cdft%is_parallel .AND. iforce_eval /= 1) &
                  DEALLOCATE (dft_control%qs_control%becke_control%cutoffs_tmp)
            ELSE

            END IF
            IF (dft_control%qs_control%becke_control%adjust) THEN
               CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
               IF (.NOT. SIZE(atomic_kind_set) == SIZE(dft_control%qs_control%becke_control%radii_tmp)) &
                  CALL cp_abort(__LOCATION__, "Length of keyword BECKE_RESTRAINT\ATOMIC_RADII does not "// &
                                "match number of atomic kinds in the input coordinate file.")
               nkinds = SIZE(dft_control%qs_control%becke_control%radii_tmp)
               IF (nkinds .GT. settings%max_nkinds) &
                  CALL cp_abort(__LOCATION__, "More than "//TRIM(cp_to_string(settings%max_nkinds))// &
                                " unique elements were defined in BECKE_RESTRAINT\ATOMIC_RADII. Are you sure"// &
                                " your input is correct? If yes, please increase max_nkinds and recompile.")
               settings%radii(1:nkinds, iforce_eval) = dft_control%qs_control%becke_control%radii_tmp(:)
               IF (.NOT. mixed_cdft%is_parallel .AND. iforce_eval /= 1) &
                  DEALLOCATE (dft_control%qs_control%becke_control%radii_tmp)
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
      is_match = is_match .AND. (settings%grid_span(1) == settings%grid_span(2))
      is_match = is_match .AND. (settings%npts(1, 1) == settings%npts(1, 2))
      is_match = is_match .AND. (settings%cutoff(1) == settings%cutoff(2))
      is_match = is_match .AND. (settings%rel_cutoff(1) == settings%rel_cutoff(2))
      is_match = is_match .AND. (settings%spherical(1) == settings%spherical(2))
      is_match = is_match .AND. (settings%rs_dims(1, 1) == settings%rs_dims(1, 2))
      is_match = is_match .AND. (settings%rs_dims(2, 1) == settings%rs_dims(2, 2))
      is_match = is_match .AND. (settings%odd(1) == settings%odd(2))
      IF (.NOT. is_match) CPABORT("Mismatch detected in the &MGRID settings of the two force_evals.")
      IF (settings%spherical(1) == 1) settings%is_spherical = .TRUE.
      IF (settings%odd(1) == 1) settings%is_odd = .TRUE.
      ! Make sure Becke settings are consistent
      CALL mp_sum(settings%atoms, force_env%para_env%group)
      CALL mp_sum(settings%coeffs, force_env%para_env%group)
      settings%nbecke = 0
      DO i = 1, SIZE(settings%atoms, 1)
         IF (settings%atoms(i, 1) /= settings%atoms(i, 2)) is_match = .FALSE.
         IF (settings%coeffs(i, 1) /= settings%coeffs(i, 2)) is_match = .FALSE.
         IF (settings%atoms(i, 1) /= 0) settings%nbecke = settings%nbecke+1
      END DO
      IF (.NOT. is_match) &
         CPABORT("Mismatch detected in the &BECKE_RESTRAINT settings of the two force_evals.")
      CALL mp_sum(settings%si, force_env%para_env%group)
      CALL mp_sum(settings%sr, force_env%para_env%group)
      CALL mp_sum(settings%confine_bounds, force_env%para_env%group)
      DO i = 1, SIZE(settings%si, 1)
         IF (settings%si(i, 1) /= settings%si(i, 2)) is_match = .FALSE.
      END DO
      DO i = 1, SIZE(settings%sr, 1)
         IF (settings%sr(i, 1) /= settings%sr(i, 2)) is_match = .FALSE.
      END DO
      IF (.NOT. is_match) &
         CPABORT("Mismatch detected in the &BECKE_RESTRAINT settings of the two force_evals.")
      is_match = is_match .AND. (settings%confine_bounds(1, 1) == settings%confine_bounds(1, 2))
      is_match = is_match .AND. (settings%confine_bounds(2, 1) == settings%confine_bounds(2, 2))
      IF (.NOT. is_match) &
         CALL cp_abort(__LOCATION__, "Mismatch detected in the &BECKE_RESTRAINT&CONFINE_BOUNDS"// &
                       " settings of the two force_evals.")
      ! Some CDFT features are currently disabled for mixed calculations: check that these features were not requested
      SELECT CASE (settings%si (1, 1))
      CASE (becke_cavity_conf, becke_none_conf)
         ! Do nothing
      CASE (becke_static_conf, becke_dynamic_conf, becke_mixed_conf)
         IF (.NOT. mixed_cdft%is_parallel) &
            CALL cp_abort(__LOCATION__, "Serial mode mixed CDFT calculation supports only confinement types"// &
                          "NONE and CAVITY.")
         IF (settings%si(9, 1) /= 3) &
            CPABORT("Only z-direction allowed for static/dynamic confinement with mixed CDFT")
      CASE DEFAULT
         CPABORT("Unknown confinement method.")
      END SELECT
      IF (settings%si(11, 1) == 1) &
         CPABORT("Calculation of atomic Becke charges not supported with mixed CDFT")
      IF (mixed_cdft%dlb .AND. .NOT. settings%si(5, 1) == 1) &
         CPABORT("Mixed CDFT load balancing requires Gaussian cavity confinement.")

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

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_transfer_settings', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, nkinds
      LOGICAL                                            :: is_match

      is_match = .TRUE.
      ! Transfer global settings
      mixed_cdft%multiplicity = settings%si(13, 1)
      mixed_cdft%has_unit_metric = .FALSE.
      IF (settings%si(14, 1) == 1) mixed_cdft%has_unit_metric = .TRUE.
      mixed_cdft%constraint_type = settings%si(15, 1)
      mixed_cdft%combined_type = settings%si(16, 1)
      mixed_cdft%eps_rho_rspace = settings%sr(5, 1)
      settings%radius = settings%sr(6, 1)
      IF (mixed_cdft%constraint_type /= cdft_density_constraint) &
         CALL cp_abort(__LOCATION__, "Unsupported constraint type. Mixed CDFT only accepts total density "// &
                       " density constraints at this moment.")
      ! Transfer settings only needed if the constraint should be built in parallel
      IF (mixed_cdft%is_parallel) THEN
         CALL becke_control_create(mixed_cdft%becke_control)
         ALLOCATE (mixed_cdft%becke_control%atoms(settings%nbecke))
         mixed_cdft%becke_control%atoms = settings%atoms(1:settings%nbecke, 1)
         ALLOCATE (mixed_cdft%becke_control%coeff(settings%nbecke))
         mixed_cdft%becke_control%coeff = settings%coeffs(1:settings%nbecke, 1)
         mixed_cdft%becke_control%natoms = settings%nbecke
         mixed_cdft%becke_control%confine_method = settings%si(1, 1)
         mixed_cdft%becke_control%cutoff_type = settings%si(2, 1)
         IF (settings%si(3, 1) == 1) mixed_cdft%becke_control%confine = .TRUE.
         IF (settings%si(4, 1) == 1) mixed_cdft%becke_control%dynamic_confine = .TRUE.
         IF (settings%si(5, 1) == 1) mixed_cdft%becke_control%cavity_confine = .TRUE.
         IF (settings%si(6, 1) == 1) mixed_cdft%becke_control%should_skip = .TRUE.
         IF (settings%si(7, 1) == 1) mixed_cdft%becke_control%print_cavity = .TRUE.
         IF (settings%si(8, 1) == 1) mixed_cdft%becke_control%in_memory = .TRUE.
         mixed_cdft%becke_control%confine_dir = settings%si(9, 1)
         IF (settings%si(10, 1) == 1) mixed_cdft%becke_control%adjust = .TRUE.
         IF (settings%si(11, 1) == 1) mixed_cdft%becke_control%atomic_charges = .TRUE.
         mixed_cdft%becke_control%cavity_shape = settings%si(12, 1)
         IF (settings%si(17, 1) == 1) mixed_cdft%becke_control%use_bohr = .TRUE.
         mixed_cdft%becke_control%dynamic_radius = settings%sr(1, 1)
         mixed_cdft%becke_control%rcavity = settings%sr(2, 1)
         mixed_cdft%becke_control%rglobal = settings%sr(3, 1)
         mixed_cdft%becke_control%eps_cavity = settings%sr(4, 1)
         mixed_cdft%becke_control%confine_bounds = settings%confine_bounds(:, 1)
         nkinds = 0
         IF (mixed_cdft%becke_control%cutoff_type == becke_cutoff_element) THEN
            CALL mp_sum(settings%cutoffs, force_env%para_env%group)
            DO i = 1, SIZE(settings%cutoffs, 1)
               IF (settings%cutoffs(i, 1) /= settings%cutoffs(i, 2)) is_match = .FALSE.
               IF (settings%cutoffs(i, 1) /= 0.0_dp) nkinds = nkinds+1
            END DO
            IF (.NOT. is_match) &
               CALL cp_abort(__LOCATION__, "Mismatch detected in the &BECKE_RESTRAINT"// &
                             "&ELEMENT_CUTOFF settings of the two force_evals.")
            ALLOCATE (mixed_cdft%becke_control%cutoffs_tmp(nkinds))
            mixed_cdft%becke_control%cutoffs_tmp = settings%cutoffs(1:nkinds, 1)
         END IF
         nkinds = 0
         IF (mixed_cdft%becke_control%adjust) THEN
            CALL mp_sum(settings%radii, force_env%para_env%group)
            DO i = 1, SIZE(settings%radii, 1)
               IF (settings%radii(i, 1) /= settings%radii(i, 2)) is_match = .FALSE.
               IF (settings%radii(i, 1) /= 0.0_dp) nkinds = nkinds+1
            END DO
            IF (.NOT. is_match) &
               CALL cp_abort(__LOCATION__, "Mismatch detected in the &BECKE_RESTRAINT"// &
                             "&ATOMIC_RADII settings of the two force_evals.")
            ALLOCATE (mixed_cdft%becke_control%radii(nkinds))
            mixed_cdft%becke_control%radii = settings%radii(1:nkinds, 1)
         END IF
      ELSE
         SELECT CASE (settings%si (1, 1))
         CASE (becke_cavity_conf, becke_none_conf)
            ! Do nothing
         CASE DEFAULT
            CALL cp_abort(__LOCATION__, "Serial mode mixed CDFT calculation supports only confinement types"// &
                          "NONE and CAVITY.")
         END SELECT
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

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_init_structures', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_path_length)                 :: c_val, input_file_path, output_file_path
      INTEGER                                            :: i, imap, iounit, j, lp, n_force_eval, &
                                                            ncpu, ntargets, offset, req(3), unit_nr
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
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%PRINT%PROGRAM_RUN_INFO")
      iounit = cp_print_key_unit_nr(logger, print_section, '', extension='.mixedLog')
      is_match = .TRUE.
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
      IF (mixed_cdft%is_parallel) THEN
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
         ! Check if the layout was succesfully created
         IF (mixed_cdft%is_special) THEN
            IF (.NOT. pw_grid%para%rs_dims(2) /= 1) is_match = .FALSE.
         ELSE IF (mixed_cdft%is_pencil) THEN
            IF (.NOT. pw_grid%para%rs_dims(1) == mixed_rs_dims(1)) is_match = .FALSE.
         ELSE
            IF (.NOT. pw_grid%para%rs_dims(2) == 1) is_match = .FALSE.
         END IF
         IF (.NOT. is_match) &
            CALL cp_abort(__LOCATION__, "Unable to create a suitable grid distribution "// &
                          "for mixed CDFT calculations. Try decreasing the total number "// &
                          " of processors or disabling xc_smoothing.")
         DEALLOCATE (mixed_rs_dims)
         ! Create the pool
         bo_mixed = pw_grid%bounds_local
         ALLOCATE (pw_pools(1))
         NULLIFY (pw_pools(1)%pool)
         CALL pw_pool_create(pw_pools(1)%pool, pw_grid=pw_grid)
         ! Initialize Gaussian cavity confinement
         IF (mixed_cdft%becke_control%cavity_confine) THEN
            CALL create_hirshfeld_type(mixed_cdft%becke_control%cavity_env)
            CALL set_hirshfeld_info(mixed_cdft%becke_control%cavity_env, &
                                    shape_function_type=shape_function_gaussian, iterative=.FALSE., &
                                    radius_type=mixed_cdft%becke_control%cavity_shape, &
                                    use_bohr=mixed_cdft%becke_control%use_bohr)
         END IF
         ! Gaussian confinement/wavefunction overlap method needs qs_kind_set
         ! Gaussian cavity confinement also needs the auxbas_rs_grid
         IF (mixed_cdft%becke_control%cavity_confine .OR. &
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
                                 nsmax=2*MAX(1, return_cube_max_iradius(mixed_cdft%pw_env%cube_info(1)))+1, &
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
            mixed_cdft%dest_list = (/imap, imap+force_env%para_env%num_pe/2/)
            imap = MOD(force_env%para_env%mepos, force_env%para_env%num_pe/2)+ &
                   MODULO(force_env%para_env%mepos, force_env%para_env%num_pe/2)
            mixed_cdft%source_list = (/imap, imap+1/)
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
            ALLOCATE (bounds(0:auxbas_pw_pool%pw_grid%para%group_size-1, 1:2))
            DO i = 0, auxbas_pw_pool%pw_grid%para%group_size-1
               bounds(i, 1:2) = auxbas_pw_pool%pw_grid%para%bo(1:2, 1, i, 1)
               bounds(i, 1:2) = bounds(i, 1:2)-auxbas_pw_pool%pw_grid%npts(1)/2-1
            END DO
            ! work out which procs to send my grid points
            ! first get the number of target procs per group
            ntargets = 0
            offset = -1
            DO i = 0, auxbas_pw_pool%pw_grid%para%group_size-1
               IF ((bounds(i, 1) .GE. bo_mixed(1, 1) .AND. bounds(i, 1) .LE. bo_mixed(2, 1)) .OR. &
                   (bounds(i, 2) .GE. bo_mixed(1, 1) .AND. bounds(i, 2) .LE. bo_mixed(2, 1))) THEN
                  ntargets = ntargets+1
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
            DO i = offset, offset+ntargets-1
               mixed_cdft%dest_list(j) = i
               mixed_cdft%dest_list_bo(:, j) = (/bo_mixed(1, 1)+(bounds(i, 1)-bo_mixed(1, 1)), &
                                                 bo_mixed(2, 1)+(bounds(i, 2)-bo_mixed(2, 1))/)
               j = j+1
            END DO
            ALLOCATE (mixed_cdft%dest_list_save(ntargets), mixed_cdft%dest_bo_save(2, ntargets))
            ! We need to store backups of these arrays since they might get reallocated during dlb
            mixed_cdft%dest_list_save = mixed_cdft%dest_list
            mixed_cdft%dest_bo_save = mixed_cdft%dest_list_bo
            ! finally determine which procs will send me grid points
            ! now we need info about y dir also
            DEALLOCATE (bounds)
            ALLOCATE (bounds(0:pw_pools(1)%pool%pw_grid%para%group_size-1, 1:4))
            DO i = 0, pw_pools(1)%pool%pw_grid%para%group_size-1
               bounds(i, 1:2) = pw_pools(1)%pool%pw_grid%para%bo(1:2, 1, i, 1)
               bounds(i, 3:4) = pw_pools(1)%pool%pw_grid%para%bo(1:2, 2, i, 1)
               bounds(i, 1:2) = bounds(i, 1:2)-pw_pools(1)%pool%pw_grid%npts(1)/2-1
               bounds(i, 3:4) = bounds(i, 3:4)-pw_pools(1)%pool%pw_grid%npts(2)/2-1
            END DO
            ntargets = 0
            offset = -1
            DO i = 0, pw_pools(1)%pool%pw_grid%para%group_size-1
               IF ((bo(1, 1) .GE. bounds(i, 1) .AND. bo(1, 1) .LE. bounds(i, 2)) .OR. &
                   (bo(2, 1) .GE. bounds(i, 1) .AND. bo(2, 1) .LE. bounds(i, 2))) THEN
                  ntargets = ntargets+1
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
            DO i = offset, offset+ntargets-1
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
               j = j+1
            END DO
            ALLOCATE (mixed_cdft%source_list_save(ntargets), mixed_cdft%source_bo_save(4, ntargets))
            ! We need to store backups of these arrays since they might get reallocated during dlb
            mixed_cdft%source_list_save = mixed_cdft%source_list
            mixed_cdft%source_bo_save = mixed_cdft%source_list_bo
            DEALLOCATE (bounds)
         END IF
      ELSE
         ! Create a logger to redirect output of the second CDFT state to a different file than the first state
         ! even when the states are treated in serial (the initial print of QS data [basis set etc] for the
         ! second state unfortunately goes to the first log file)
         CALL force_env_get(force_env, root_section=root_section)
         IF (force_env%para_env%mepos == force_env%para_env%source) THEN
            CALL section_vals_val_get(root_section, "GLOBAL%PROJECT_NAME", &
                                      c_val=input_file_path)
            lp = LEN_TRIM(input_file_path)
            input_file_path(lp+1:LEN(input_file_path)) = "-r-2"
            lp = LEN_TRIM(input_file_path)
            output_file_path = input_file_path(1:lp)//".out"
            CALL open_file(file_name=output_file_path, file_status="UNKNOWN", &
                           file_action="WRITE", file_position="APPEND", &
                           unit_number=unit_nr)
         ELSE
            unit_nr = -1
         END IF
         CALL cp_logger_create(mixed_cdft%sub_logger, &
                               para_env=force_env%para_env, &
                               default_global_unit_nr=unit_nr, &
                               close_global_unit_on_dealloc=.FALSE.)
         ! Try to use better names for the local log if it is not too late
         CALL section_vals_val_get(root_section, "GLOBAL%OUTPUT_FILE_NAME", &
                                   c_val=c_val)
         IF (c_val /= "") THEN
            CALL cp_logger_set(mixed_cdft%sub_logger, &
                               local_filename=TRIM(c_val)//"_localLog")
         END IF
         CALL section_vals_val_get(root_section, "GLOBAL%PROJECT", c_val=c_val)
         IF (c_val /= "") THEN
            CALL cp_logger_set(mixed_cdft%sub_logger, &
                               local_filename=TRIM(c_val)//"_localLog")
         END IF
         mixed_cdft%sub_logger%iter_info%project_name = c_val
         CALL section_vals_val_get(root_section, "GLOBAL%PRINT_LEVEL", &
                                   i_val=mixed_cdft%sub_logger%iter_info%print_level)
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
      DEALLOCATE (settings%coeffs)
      DEALLOCATE (settings%cutoffs)
      DEALLOCATE (settings%radii)
      DEALLOCATE (settings%si)
      DEALLOCATE (settings%sr)
      DEALLOCATE (settings%confine_bounds)
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
                                        "MIXED%PRINT%PROGRAM_RUN_INFO")

   END SUBROUTINE mixed_cdft_init_structures

! **************************************************************************************************
!> \brief Redistribute arrays needed for an ET coupling calculation from individual CDFT states to
!>        the mixed CDFT env, that is, move the arrays to the correct blacs context. For parallel
!>        simulations, the array processor distributions also change from N to 2N processors.
!> \param force_env the force_env that holds the CDFT states
!> \param mixed_wmat_a matrix representation of the weight function of state A
!> \param mixed_wmat_b matrix representation of the weight function of state B
!> \param mixed_matrix_s atomic orbital overlap matrix
!> \param mixed_mo_coeff molecular orbitals coefficients for both states
!> \param density_matrix density matrices for both states
!> \par History
!>       01.2017  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_redistribute_arrays(force_env, mixed_wmat_a, mixed_wmat_b, &
                                             mixed_matrix_s, mixed_mo_coeff, density_matrix)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(dbcsr_type)                                   :: mixed_wmat_a, mixed_wmat_b
      TYPE(dbcsr_type), POINTER                          :: mixed_matrix_s
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: mixed_mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: density_matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_redistribute_arrays', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iforce_eval, ispin, ncol_overlap, &
                                                            ncol_wmat, nforce_eval, nrow_overlap, &
                                                            nrow_wmat, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ncol_mo, nrow_mo
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mixed_wmat_tmp, wmat_tmp
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: matrix_p_tmp, mixed_matrix_p_tmp, &
                                                            mo_coeff_tmp
      TYPE(cp_fm_p_type), POINTER                        :: matrix_s_tmp, mixed_matrix_s_tmp
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_mo, fm_struct_overlap, &
                                                            fm_struct_tmp, fm_struct_wmat
      TYPE(dbcsr_type)                                   :: desymm_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(mixed_environment_type), POINTER              :: mixed_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      NULLIFY (mixed_env, mixed_cdft, qs_env, dft_control, fm_struct_mo, &
               fm_struct_wmat, fm_struct_overlap, fm_struct_tmp, &
               mixed_matrix_s_tmp, matrix_s_tmp, mixed_wmat_tmp, &
               wmat_tmp, mixed_mo_coeff, matrix_p_tmp, mixed_matrix_p_tmp, &
               mixed_matrix_s, density_matrix, blacs_env)

      CPASSERT(ASSOCIATED(force_env))
      mixed_env => force_env%mixed_env
      nforce_eval = SIZE(force_env%sub_force_env)
      CALL get_mixed_env(mixed_env, cdft_control=mixed_cdft)
      CPASSERT(ASSOCIATED(mixed_cdft))
      ! Get nspins
      DO iforce_eval = 1, nforce_eval
         IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
         IF (force_env%mixed_env%do_mixed_qmmm_cdft) THEN
            qs_env => force_env%sub_force_env(iforce_eval)%force_env%qmmm_env%qs_env
         ELSE
            CALL force_env_get(force_env%sub_force_env(iforce_eval)%force_env, qs_env=qs_env)
         END IF
         CALL get_qs_env(qs_env, dft_control=dft_control)
         nspins = dft_control%nspins
      END DO
      ALLOCATE (mixed_mo_coeff(2, nspins))
      CALL dbcsr_init_p(mixed_matrix_s)
      IF (mixed_cdft%calculate_metric) &
         ALLOCATE (density_matrix(2, nspins))
      ! Transfer data from sub_force_envs to temporaries
      ALLOCATE (mo_coeff_tmp(2, nspins), wmat_tmp(2), matrix_s_tmp, nrow_mo(nspins), ncol_mo(nspins))
      IF (mixed_cdft%calculate_metric) ALLOCATE (matrix_p_tmp(2, nspins))
      NULLIFY (matrix_s_tmp%matrix)
      DO iforce_eval = 1, nforce_eval
         ! Temporary arrays need to be nulled on every process
         NULLIFY (wmat_tmp(iforce_eval)%matrix)
         DO ispin = 1, nspins
            NULLIFY (mo_coeff_tmp(iforce_eval, ispin)%matrix)
            IF (mixed_cdft%calculate_metric) &
               NULLIFY (matrix_p_tmp(iforce_eval, ispin)%matrix)
         END DO
         IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
         ! From this point onward, we access data local to the sub_force_envs
         ! Get qs_env
         IF (force_env%mixed_env%do_mixed_qmmm_cdft) THEN
            qs_env => force_env%sub_force_env(iforce_eval)%force_env%qmmm_env%qs_env
         ELSE
            CALL force_env_get(force_env%sub_force_env(iforce_eval)%force_env, qs_env=qs_env)
         END IF
         CALL get_qs_env(qs_env, dft_control=dft_control, blacs_env=blacs_env)
         ! Store dimensions of the transferred arrays
         CALL dbcsr_get_info(dft_control%qs_control%cdft_control%matrix_s%matrix, &
                             nfullrows_total=nrow_overlap, nfullcols_total=ncol_overlap)
         CALL dbcsr_get_info(dft_control%qs_control%cdft_control%wmat%matrix, &
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
         ! Matrix representation of weight function (dbcsr -> fm)
         CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nrow_wmat, ncol_global=ncol_wmat, context=blacs_env, &
                                  para_env=force_env%sub_force_env(iforce_eval)%force_env%para_env, &
                                  square_blocks=.TRUE.)
         CALL cp_fm_create(wmat_tmp(iforce_eval)%matrix, fm_struct_tmp, name="fm_matrix")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL dbcsr_desymmetrize(dft_control%qs_control%cdft_control%wmat%matrix, desymm_tmp)
         CALL copy_dbcsr_to_fm(desymm_tmp, wmat_tmp(iforce_eval)%matrix)
         CALL dbcsr_release(desymm_tmp)
         CALL dbcsr_release_p(dft_control%qs_control%cdft_control%wmat%matrix)
         ! Overlap matrix is the same for both sub_force_envs, so we just copy the first one (dbcsr -> fm)
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
      END DO
      ! Create needed fm structs
      CALL cp_fm_struct_create(fm_struct_wmat, nrow_global=nrow_wmat, ncol_global=ncol_wmat, &
                               context=mixed_cdft%blacs_env, para_env=force_env%para_env)
      CALL cp_fm_struct_create(fm_struct_overlap, nrow_global=nrow_overlap, ncol_global=ncol_overlap, &
                               context=mixed_cdft%blacs_env, para_env=force_env%para_env)
      ! Redistribute arrays with copy_general (this is not optimal for dbcsr matrices but...)
      ! We use this method for the serial case (mixed_cdft%is_parallel = .FALSE.) as well to move the arrays to the
      ! correct blacs_env, which is impossible using a simple copy of the arrays
      ALLOCATE (mixed_wmat_tmp(2), mixed_matrix_s_tmp)
      IF (mixed_cdft%calculate_metric) &
         ALLOCATE (mixed_matrix_p_tmp(2, nspins))
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
         NULLIFY (mixed_wmat_tmp(iforce_eval)%matrix)
         CALL cp_fm_create(matrix=mixed_wmat_tmp(iforce_eval)%matrix, &
                           matrix_struct=fm_struct_wmat, &
                           name="WEIGHT_"//TRIM(ADJUSTL(cp_to_string(iforce_eval)))//"_MATRIX")
         CALL cp_fm_copy_general(wmat_tmp(iforce_eval)%matrix, &
                                 mixed_wmat_tmp(iforce_eval)%matrix, &
                                 mixed_cdft%blacs_env%para_env)
         CALL cp_fm_release(wmat_tmp(iforce_eval)%matrix)
         ! (fm -> dbcsr)
         IF (iforce_eval == 1) THEN
            CALL copy_fm_to_dbcsr_bc(mixed_wmat_tmp(iforce_eval)%matrix, mixed_wmat_a)
         ELSE
            CALL copy_fm_to_dbcsr_bc(mixed_wmat_tmp(iforce_eval)%matrix, mixed_wmat_b)
         END IF
         CALL cp_fm_release(mixed_wmat_tmp(iforce_eval)%matrix)
         ! Density matrix (fm -> dbcsr)
         IF (mixed_cdft%calculate_metric) THEN
            DO ispin = 1, nspins
               NULLIFY (mixed_matrix_p_tmp(iforce_eval, ispin)%matrix)
               NULLIFY (density_matrix(iforce_eval, ispin)%matrix)
               CALL dbcsr_init_p(density_matrix(iforce_eval, ispin)%matrix)
               CALL cp_fm_create(matrix=mixed_matrix_p_tmp(iforce_eval, ispin)%matrix, &
                                 matrix_struct=fm_struct_overlap, &
                                 name="DENSITY_"//TRIM(ADJUSTL(cp_to_string(iforce_eval)))//"_" &
                                 //TRIM(ADJUSTL(cp_to_string(ispin)))//"_MATRIX")
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
      DEALLOCATE (matrix_s_tmp, mixed_matrix_s_tmp)
      DEALLOCATE (ncol_mo, nrow_mo)

   END SUBROUTINE mixed_cdft_redistribute_arrays

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

      CHARACTER(len=*), PARAMETER :: routineN = 'hfun_zero', routineP = moduleN//':'//routineN

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
                     nzeroed = nzeroed+1
                     nzeroed_total = nzeroed_total+1
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
         bounds = bounds-(n3/2)-1
         work = n3*n2*n1-nzeroed_total
      END IF

   END SUBROUTINE hfun_zero

END MODULE mixed_cdft_utils
