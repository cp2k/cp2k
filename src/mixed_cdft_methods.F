!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods for mixed CDFT calculations
!> \par   History
!>                 Separated CDFT routines from mixed_environment_utils
!> \author Nico Holmberg [01.2017]
! **************************************************************************************************
MODULE mixed_cdft_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_retain,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: becke_control_create,&
                                              dft_control_type
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_syevd
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr_bc,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_invert,&
                                              cp_fm_transpose
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE cp_realspace_grid_init,          ONLY: init_input_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE cube_utils,                      ONLY: init_cube_info,&
                                              return_cube_max_iradius
   USE d3_poly,                         ONLY: init_d3_poly_module
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_create, dbcsr_desymmetrize, dbcsr_get_info, dbcsr_init_p, dbcsr_p_type, &
        dbcsr_release, dbcsr_release_p, dbcsr_scale, dbcsr_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type,&
                                              multiple_fe_list,&
                                              use_qmmm,&
                                              use_qmmmx,&
                                              use_qs_force
   USE gaussian_gridlevels,             ONLY: init_gaussian_gridlevel
   USE global_types,                    ONLY: global_environment_type
   USE hirshfeld_methods,               ONLY: create_shape_function
   USE hirshfeld_types,                 ONLY: create_hirshfeld_type,&
                                              hirshfeld_type,&
                                              release_hirshfeld_type,&
                                              set_hirshfeld_info
   USE input_constants,                 ONLY: becke_cutoff_element,&
                                              becke_cutoff_global,&
                                              cdft_density_constraint,&
                                              do_spin_density,&
                                              shape_function_gaussian
   USE input_section_types,             ONLY: section_vals_duplicate,&
                                              section_vals_get_subs_vals,&
                                              section_vals_release,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp,&
                                              int_8
   USE machine,                         ONLY: m_walltime
   USE mathlib,                         ONLY: diamat_all
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_irecv,&
                                              mp_isend,&
                                              mp_sum,&
                                              mp_test,&
                                              mp_testall,&
                                              mp_wait,&
                                              mp_waitall
   USE mixed_cdft_types,                ONLY: mixed_cdft_type,&
                                              mixed_cdft_type_create
   USE mixed_environment_types,         ONLY: get_mixed_env,&
                                              mixed_environment_type,&
                                              set_mixed_env
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_methods,                  ONLY: pw_env_create
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: HALFSPACE,&
                                              pw_grid_type
   USE pw_grids,                        ONLY: do_pw_grid_blocked_false,&
                                              pw_grid_create,&
                                              pw_grid_release,&
                                              pw_grid_setup
   USE pw_methods,                      ONLY: pw_scale
   USE pw_pool_types,                   ONLY: pw_pool_create,&
                                              pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE
   USE qs_collocate_density,            ONLY: collocate_pgf_product_rspace
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: create_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_io,                        ONLY: read_mo_set,&
                                              wfn_restart_file_name
   USE qs_mo_methods,                   ONLY: make_basis_simple,&
                                              make_basis_sm
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              mo_set_p_type,&
                                              set_mo_set
   USE qs_modify_pab_block,             ONLY: FUNC_AB
   USE realspace_grid_types,            ONLY: &
        realspace_grid_desc_p_type, realspace_grid_input_type, realspace_grid_p_type, &
        realspace_grid_type, rs2pw, rs_grid_create, rs_grid_create_descriptor, rs_grid_print, &
        rs_grid_release, rs_grid_retain, rs_grid_zero, rs_pw_transfer
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mixed_cdft_methods'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

   PUBLIC :: mixed_cdft_init, &
             mixed_cdft_build_weight, &
             mixed_cdft_calculate_coupling

CONTAINS

! **************************************************************************************************
!> \brief Initialize a mixed CDFT calculation (currently only relevant with Becke constraint)
!> \param force_env the force_env that holds the CDFT states
!> \param calculate_forces determines if forces should be calculted
!> \par History
!>       01.2016  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_init(force_env, calculate_forces)
      TYPE(force_env_type), POINTER                      :: force_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_init', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, i, iforce_eval, imap, iounit, j, max_nkinds, n_force_eval, nbecke, ncpu, &
         nforce_eval, nkinds, ntargets, offset, req(3)
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: bounds
      INTEGER, DIMENSION(2, 3)                           :: bo, bo_mixed
      INTEGER, DIMENSION(3)                              :: higher_grid_layout
      INTEGER, DIMENSION(:), POINTER                     :: grid_span, i_force_eval, mixed_rs_dims, &
                                                            odd, recvbuffer, recvbuffer2, &
                                                            sendbuffer, spherical
      INTEGER, DIMENSION(:, :), POINTER                  :: atoms, npts, rs_dims, settings_int
      LOGICAL                                            :: is_match, is_odd, is_qmmm, is_spherical
      REAL(KIND=dp)                                      :: radius
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cutoff, rel_cutoff
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: coeffs, confine_bounds, cutoffs, radii, &
                                                            settings_real
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell_mix
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys_mix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(force_env_type), POINTER                      :: force_env_qs
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(mixed_environment_type), POINTER              :: mixed_env
      TYPE(particle_list_type), POINTER                  :: particles_mix
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
         md_section, mixed_section, print_section, root_section, rs_grid_section, subsys_section

      NULLIFY (cell_mix, subsys_mix, force_env_qs, force_env_section, subsys_section, &
               print_section, root_section, kind_section, force_env_sections, rs_grid_section, &
               mixed_section, md_section, mixed_env, mixed_cdft, pw_env, auxbas_pw_pool, &
               pw_pools, qs_env, pw_grid, dft_control, sendbuffer, grid_span, odd, spherical, &
               mixed_rs_dims, i_force_eval, recvbuffer, recvbuffer2, npts, rs_dims, atoms, &
               settings_int, cutoff, rel_cutoff, coeffs, settings_real, confine_bounds, &
               cutoffs, blacs_env, globenv, atomic_kind_set, qs_kind_set, rs_descs, rs_grids)

      is_qmmm = .FALSE.
      logger => cp_get_default_logger()
      CPASSERT(ASSOCIATED(force_env))
      CALL timeset(routineN, handle)
      ! First determine if the calculation is pure DFT or QMMM and find the qs force_env
      nforce_eval = SIZE(force_env%sub_force_env)
      DO iforce_eval = 1, nforce_eval
         IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
         SELECT CASE (force_env%sub_force_env (iforce_eval)%force_env%in_use)
         CASE (use_qs_force)
            force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
         CASE (use_qmmm)
            is_qmmm = .TRUE.
            ! This is really the container for QMMM
            force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
         CASE (use_qmmmx)
            CALL cp_abort(__LOCATION__, &
                          "No force mixing allowed for mixed CDFT QM/MM")
         CASE DEFAULT
            CPASSERT(.FALSE.)
         END SELECT
      END DO
      CPASSERT(ASSOCIATED(force_env_qs))
      CALL force_env_get(force_env=force_env, force_env_section=force_env_section)
      ! Get infos about the mixed subsys
      IF (.NOT. is_qmmm) THEN
         CALL force_env_get(force_env=force_env, &
                            subsys=subsys_mix, &
                            cell=cell_mix)
         CALL cp_subsys_get(subsys=subsys_mix, &
                            particles=particles_mix)
      ELSE
         CALL get_qs_env(force_env_qs%qmmm_env%qs_env, &
                         cp_subsys=subsys_mix, &
                         cell=cell_mix)
         CALL cp_subsys_get(subsys=subsys_mix, &
                            particles=particles_mix)
      END IF
      ncpu = force_env%para_env%num_pe
      mixed_env => force_env%mixed_env
      mixed_section => section_vals_get_subs_vals(force_env_section, "MIXED")
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%PRINT%PROGRAM_RUN_INFO")
      iounit = cp_print_key_unit_nr(logger, print_section, '', extension='.mixedLog')
      ! Init mixed_cdft_type
      ALLOCATE (mixed_cdft)
      CALL mixed_cdft_type_create(mixed_cdft)
      ! Store QMMM flag
      mixed_env%do_mixed_qmmm_cdft = is_qmmm
      ! Setup dynamic load balancing
      CALL section_vals_val_get(mixed_section, "LINEAR%MIXED_CDFT_DLB", l_val=mixed_cdft%dlb)
      mixed_cdft%dlb = mixed_cdft%dlb .AND. calculate_forces ! disable if forces are not needed
      IF (mixed_cdft%dlb) THEN
         ALLOCATE (mixed_cdft%dlb_control)
         NULLIFY (mixed_cdft%dlb_control%weight, mixed_cdft%dlb_control%gradients, &
                  mixed_cdft%dlb_control%cavity, mixed_cdft%dlb_control%target_list, &
                  mixed_cdft%dlb_control%bo, mixed_cdft%dlb_control%expected_work, &
                  mixed_cdft%dlb_control%prediction_error, mixed_cdft%dlb_control%sendbuff, &
                  mixed_cdft%dlb_control%recvbuff, mixed_cdft%dlb_control%recv_work_repl, &
                  mixed_cdft%dlb_control%recv_info)
         CALL section_vals_val_get(mixed_section, "LINEAR%LOAD_SCALE", r_val=mixed_cdft%dlb_control%load_scale)
         CALL section_vals_val_get(mixed_section, "LINEAR%VERY_OVERLOADED", &
                                   r_val=mixed_cdft%dlb_control%very_overloaded)
         CALL section_vals_val_get(mixed_section, "LINEAR%MORE_WORK", i_val=mixed_cdft%dlb_control%more_work)
      END IF
      ! Metric/Wavefunction overlap method
      mixed_cdft%calculate_metric = .FALSE.
      mixed_cdft%wfn_overlap_method = .FALSE.
      IF (mixed_env%do_mixed_et) THEN
         CALL section_vals_val_get(mixed_section, "LINEAR%MIXED_CDFT_METRIC", l_val=mixed_cdft%calculate_metric)
         CALL section_vals_val_get(mixed_section, "LINEAR%MIXED_CDFT_WFN_OVERLAP", l_val=mixed_cdft%wfn_overlap_method)
      END IF
      ! Inversion method
      CALL section_vals_val_get(mixed_section, "LINEAR%EPS_SVD", r_val=mixed_cdft%eps_svd)
      IF (mixed_cdft%eps_svd .LT. 0.0_dp .OR. mixed_cdft%eps_svd .GT. 1.0_dp) &
         CALL cp_abort(__LOCATION__, &
                       "Illegal value for EPS_SVD. Value must be between 0.0 and 1.0.")
      ! Start building the mixed auxbas_pw_pool
      CALL pw_grid_create(pw_grid, force_env%para_env%group)
      CALL pw_env_create(mixed_cdft%pw_env)
      ! Allocate storage for temporaries
      ALLOCATE (grid_span(nforce_eval))
      grid_span = 0
      ALLOCATE (npts(3, nforce_eval))
      npts = 0
      ALLOCATE (cutoff(nforce_eval))
      cutoff = 0.0_dp
      ALLOCATE (rel_cutoff(nforce_eval))
      rel_cutoff = 0.0_dp
      ALLOCATE (spherical(nforce_eval))
      spherical = 0
      is_spherical = .FALSE.
      ALLOCATE (rs_dims(2, nforce_eval))
      rs_dims = 0
      ALLOCATE (odd(nforce_eval))
      odd = 0
      is_odd = .FALSE.
      ALLOCATE (atoms(SIZE(particles_mix%els), nforce_eval))
      atoms = 0
      ALLOCATE (coeffs(SIZE(particles_mix%els), nforce_eval))
      coeffs = 0.0_dp
      ALLOCATE (settings_int(17, nforce_eval))
      settings_int = 0
      ALLOCATE (settings_real(6, nforce_eval))
      settings_real = 0.0_dp
      ALLOCATE (confine_bounds(2, nforce_eval))
      confine_bounds = 0.0_dp
      max_nkinds = 30
      ALLOCATE (cutoffs(max_nkinds, nforce_eval))
      cutoffs = 0.0_dp
      ALLOCATE (radii(max_nkinds, nforce_eval))
      radii = 0.0_dp
      ! Get information from the sub_force_envs
      ! TODO: transfer logicals as logicals instead of ints
      DO iforce_eval = 1, nforce_eval
         IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
         IF (is_qmmm) THEN
            qs_env => force_env_qs%qmmm_env%qs_env
         ELSE
            CALL force_env_get(force_env_qs, qs_env=qs_env)
         END IF
         CALL get_qs_env(qs_env, pw_env=pw_env, dft_control=dft_control)
         IF (.NOT. dft_control%qs_control%becke_restraint) &
            CALL cp_abort(__LOCATION__, &
                          "Mixed CDFT currently only supports the Becke constraint."// &
                          "Becke constraint must be active in QS section of both force_evals!")
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         bo = auxbas_pw_pool%pw_grid%bounds_local
         ! Only the rank 0 process collects info about pw_grid and becke
         IF (force_env_qs%para_env%mepos == &
             force_env_qs%para_env%source) THEN
            grid_span(iforce_eval) = auxbas_pw_pool%pw_grid%grid_span
            npts(:, iforce_eval) = auxbas_pw_pool%pw_grid%npts
            cutoff(iforce_eval) = auxbas_pw_pool%pw_grid%cutoff
            rel_cutoff(iforce_eval) = dft_control%qs_control%relative_cutoff
            IF (auxbas_pw_pool%pw_grid%spherical) spherical(iforce_eval) = 1
            rs_dims(:, iforce_eval) = auxbas_pw_pool%pw_grid%para%rs_dims
            IF (auxbas_pw_pool%pw_grid%grid_span == HALFSPACE) odd(iforce_eval) = 1
            IF (dft_control%qs_control%becke_control%natoms .GT. SIZE(atoms, 1)) &
               CALL cp_abort(__LOCATION__, &
                             "More Becke constraint atoms than defined in mixed section."// &
                             " Use default values for MIXED\MAPPING.")
            atoms(1:dft_control%qs_control%becke_control%natoms, iforce_eval) = dft_control%qs_control%becke_control%atoms
            coeffs(1:dft_control%qs_control%becke_control%natoms, iforce_eval) = dft_control%qs_control%becke_control%coeff
            ! Only GAPW needs to keep a copy of these arrays for constraint integration
            IF (.NOT. dft_control%qs_control%gapw) THEN
               DEALLOCATE (dft_control%qs_control%becke_control%atoms)
               DEALLOCATE (dft_control%qs_control%becke_control%coeff)
            END IF
            settings_int(1, iforce_eval) = dft_control%qs_control%becke_control%confine_method
            settings_int(2, iforce_eval) = dft_control%qs_control%becke_control%cutoff_type
            IF (dft_control%qs_control%becke_control%confine) settings_int(3, iforce_eval) = 1
            IF (dft_control%qs_control%becke_control%dynamic_confine) settings_int(4, iforce_eval) = 1
            IF (dft_control%qs_control%becke_control%cavity_confine) THEN
               settings_int(5, iforce_eval) = 1
               CALL release_hirshfeld_type(dft_control%qs_control%becke_control%cavity_env)
            END IF
            IF (dft_control%qs_control%becke_control%should_skip) settings_int(6, iforce_eval) = 1
            IF (dft_control%qs_control%becke_control%print_cavity) settings_int(7, iforce_eval) = 1
            IF (dft_control%qs_control%becke_control%in_memory) settings_int(8, iforce_eval) = 1
            settings_int(9, iforce_eval) = dft_control%qs_control%becke_control%confine_dir
            IF (dft_control%qs_control%becke_control%adjust) settings_int(10, iforce_eval) = 1
            IF (dft_control%qs_control%becke_control%atomic_charges) settings_int(11, iforce_eval) = 1
            settings_int(12, iforce_eval) = dft_control%qs_control%becke_control%cavity_shape
            settings_int(13, iforce_eval) = dft_control%multiplicity
            IF (qs_env%has_unit_metric) THEN
               settings_int(14, iforce_eval) = 1
            ELSE
               settings_int(14, iforce_eval) = 0
            END IF
            settings_int(15, iforce_eval) = dft_control%qs_control%becke_control%constraint_type
            settings_int(16, iforce_eval) = dft_control%qs_control%becke_control%combined_type
            IF (dft_control%qs_control%becke_control%use_bohr) THEN
               settings_int(17, iforce_eval) = 1
            ELSE
               settings_int(17, iforce_eval) = 0
            END IF
            settings_real(1, iforce_eval) = dft_control%qs_control%becke_control%dynamic_radius
            settings_real(2, iforce_eval) = dft_control%qs_control%becke_control%rcavity
            settings_real(3, iforce_eval) = dft_control%qs_control%becke_control%rglobal
            settings_real(4, iforce_eval) = dft_control%qs_control%becke_control%eps_cavity
            settings_real(5, iforce_eval) = dft_control%qs_control%eps_rho_rspace
            settings_real(6, iforce_eval) = pw_env%cube_info(pw_env%auxbas_grid)%max_rad_ga
            confine_bounds(:, iforce_eval) = dft_control%qs_control%becke_control%confine_bounds
            IF (dft_control%qs_control%becke_control%cutoff_type == becke_cutoff_element) THEN
               nkinds = SIZE(dft_control%qs_control%becke_control%cutoffs_tmp)
               IF (nkinds .GT. max_nkinds) &
                  CALL cp_abort(__LOCATION__, &
                                "More than "//TRIM(cp_to_string(max_nkinds))//" unique elements"// &
                                "were defined in BECKE_RESTRAINT\ELEMENT_CUTOFF. Are you sure"// &
                                " your input is correct? If yes, please increase max_nkinds and recompile.")
               cutoffs(1:nkinds, iforce_eval) = dft_control%qs_control%becke_control%cutoffs_tmp(:)
               DEALLOCATE (dft_control%qs_control%becke_control%cutoffs_tmp)
            END IF
            IF (dft_control%qs_control%becke_control%adjust) THEN
               CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
               IF (.NOT. SIZE(atomic_kind_set) == SIZE(dft_control%qs_control%becke_control%radii_tmp)) &
                  CALL cp_abort(__LOCATION__, &
                                "Length of keyword BECKE_RESTRAINT\ATOMIC_RADII does not "// &
                                "match number of atomic kinds in the input coordinate file.")
               nkinds = SIZE(dft_control%qs_control%becke_control%radii_tmp)
               IF (nkinds .GT. max_nkinds) &
                  CALL cp_abort(__LOCATION__, &
                                "More than "//TRIM(cp_to_string(max_nkinds))//" unique elements"// &
                                "were defined in BECKE_RESTRAINT\ATOMIC_RADII. Are you sure"// &
                                " your input is correct? If yes, please increase max_nkinds and recompile.")
               radii(1:nkinds, iforce_eval) = dft_control%qs_control%becke_control%radii_tmp(:)
               DEALLOCATE (dft_control%qs_control%becke_control%radii_tmp)
            END IF
         END IF
      END DO
      ! Make sure the grids are consistent
      CALL mp_sum(grid_span, force_env%para_env%group)
      CALL mp_sum(npts, force_env%para_env%group)
      CALL mp_sum(cutoff, force_env%para_env%group)
      CALL mp_sum(rel_cutoff, force_env%para_env%group)
      CALL mp_sum(spherical, force_env%para_env%group)
      CALL mp_sum(rs_dims, force_env%para_env%group)
      CALL mp_sum(odd, force_env%para_env%group)
      is_match = .TRUE.
      is_match = is_match .AND. (grid_span(1) == grid_span(2))
      is_match = is_match .AND. (npts(1, 1) == npts(1, 2))
      is_match = is_match .AND. (cutoff(1) == cutoff(2))
      is_match = is_match .AND. (rel_cutoff(1) == rel_cutoff(2))
      is_match = is_match .AND. (spherical(1) == spherical(2))
      is_match = is_match .AND. (rs_dims(1, 1) == rs_dims(1, 2))
      is_match = is_match .AND. (rs_dims(2, 1) == rs_dims(2, 2))
      is_match = is_match .AND. (odd(1) == odd(2))
      IF (.NOT. is_match) &
         CALL cp_abort(__LOCATION__, &
                       "Mismatch detected in the &MGRID settings of the two force_evals.")
      IF (spherical(1) == 1) is_spherical = .TRUE.
      IF (odd(1) == 1) is_odd = .TRUE.
      ! Make sure Becke settings are consistent and transfer them
      CALL mp_sum(atoms, force_env%para_env%group)
      CALL mp_sum(coeffs, force_env%para_env%group)
      nbecke = 0
      nkinds = 0
      DO i = 1, SIZE(atoms, 1)
         IF (atoms(i, 1) /= atoms(i, 2)) is_match = .FALSE.
         IF (coeffs(i, 1) /= coeffs(i, 2)) is_match = .FALSE.
         IF (atoms(i, 1) /= 0) nbecke = nbecke+1
      END DO
      IF (.NOT. is_match) &
         CALL cp_abort(__LOCATION__, &
                       "Mismatch detected in the &BECKE_RESTRAINT settings of the two force_evals.")
      CALL becke_control_create(mixed_cdft%becke_control)
      ALLOCATE (mixed_cdft%becke_control%atoms(nbecke))
      mixed_cdft%becke_control%atoms = atoms(1:nbecke, 1)
      ALLOCATE (mixed_cdft%becke_control%coeff(nbecke))
      mixed_cdft%becke_control%coeff = coeffs(1:nbecke, 1)
      mixed_cdft%becke_control%natoms = nbecke
      CALL mp_sum(settings_int, force_env%para_env%group)
      CALL mp_sum(settings_real, force_env%para_env%group)
      CALL mp_sum(confine_bounds, force_env%para_env%group)
      DO i = 1, SIZE(settings_int, 1)
         IF (settings_int(i, 1) /= settings_int(i, 2)) is_match = .FALSE.
         IF (i <= 5) THEN
            IF (settings_real(i, 1) /= settings_real(i, 2)) is_match = .FALSE.
         END IF
      END DO
      IF (.NOT. is_match) &
         CALL cp_abort(__LOCATION__, &
                       "Mismatch detected in the &BECKE_RESTRAINT settings of the two force_evals.")
      is_match = is_match .AND. (confine_bounds(1, 1) == confine_bounds(1, 2))
      is_match = is_match .AND. (confine_bounds(2, 1) == confine_bounds(2, 2))
      IF (.NOT. is_match) &
         CALL cp_abort(__LOCATION__, &
                       "Mismatch detected in the &BECKE_RESTRAINT&CONFINE_BOUNDS"// &
                       " settings of the two force_evals.")
      mixed_cdft%becke_control%confine_method = settings_int(1, 1)
      mixed_cdft%becke_control%cutoff_type = settings_int(2, 1)
      nkinds = 0
      IF (mixed_cdft%becke_control%cutoff_type == becke_cutoff_element) THEN
         CALL mp_sum(cutoffs, force_env%para_env%group)
         DO i = 1, SIZE(cutoffs, 1)
            IF (cutoffs(i, 1) /= cutoffs(i, 2)) is_match = .FALSE.
            IF (cutoffs(i, 1) /= 0.0_dp) nkinds = nkinds+1
         END DO
         IF (.NOT. is_match) &
            CALL cp_abort(__LOCATION__, &
                          "Mismatch detected in the &BECKE_RESTRAINT&ELEMENT_CUTOFF"// &
                          " settings of the two force_evals.")
         ALLOCATE (mixed_cdft%becke_control%cutoffs_tmp(nkinds))
         mixed_cdft%becke_control%cutoffs_tmp = cutoffs(1:nkinds, 1)
      END IF
      IF (settings_int(3, 1) == 1) mixed_cdft%becke_control%confine = .TRUE.
      IF (settings_int(4, 1) == 1) mixed_cdft%becke_control%dynamic_confine = .TRUE.
      IF (settings_int(5, 1) == 1) mixed_cdft%becke_control%cavity_confine = .TRUE.
      IF (settings_int(6, 1) == 1) mixed_cdft%becke_control%should_skip = .TRUE.
      IF (settings_int(7, 1) == 1) mixed_cdft%becke_control%print_cavity = .TRUE.
      IF (settings_int(8, 1) == 1) mixed_cdft%becke_control%in_memory = .TRUE.
      mixed_cdft%becke_control%confine_dir = settings_int(9, 1)
      IF (settings_int(10, 1) == 1) mixed_cdft%becke_control%adjust = .TRUE.
      IF (mixed_cdft%becke_control%adjust) THEN
         CALL mp_sum(radii, force_env%para_env%group)
         DO i = 1, SIZE(radii, 1)
            IF (radii(i, 1) /= radii(i, 2)) is_match = .FALSE.
            IF (radii(i, 1) /= 0.0_dp) nkinds = nkinds+1
         END DO
         IF (.NOT. is_match) &
            CALL cp_abort(__LOCATION__, &
                          "Mismatch detected in the &BECKE_RESTRAINT&ATOMIC_RADII"// &
                          " settings of the two force_evals.")
         ALLOCATE (mixed_cdft%becke_control%radii(nkinds))
         mixed_cdft%becke_control%radii = radii(1:nkinds, 1)
      END IF
      IF (settings_int(11, 1) == 1) mixed_cdft%becke_control%atomic_charges = .TRUE.
      mixed_cdft%becke_control%cavity_shape = settings_int(12, 1)
      mixed_cdft%multiplicity = settings_int(13, 1)
      mixed_cdft%has_unit_metric = .FALSE.
      IF (settings_int(14, 1) == 1) mixed_cdft%has_unit_metric = .FALSE.
      mixed_cdft%constraint_type = settings_int(15, 1)
      IF (mixed_cdft%constraint_type /= cdft_density_constraint) &
         CALL cp_abort(__LOCATION__, &
                       "Unsupported constraint type. Mixed CDFT only accepts total density "// &
                       " density constraints at this moment.")
      mixed_cdft%combined_type = settings_int(16, 1)
      IF (settings_int(17, 1) == 1) mixed_cdft%becke_control%use_bohr = .TRUE.
      mixed_cdft%becke_control%dynamic_radius = settings_real(1, 1)
      mixed_cdft%becke_control%rcavity = settings_real(2, 1)
      mixed_cdft%becke_control%rglobal = settings_real(3, 1)
      mixed_cdft%becke_control%eps_cavity = settings_real(4, 1)
      mixed_cdft%eps_rho_rspace = settings_real(5, 1)
      mixed_cdft%becke_control%confine_bounds = confine_bounds(:, 1)
      radius = settings_real(6, 1)
      ! Initialize Gaussian cavity confinement
      IF (mixed_cdft%becke_control%cavity_confine) THEN
         CALL create_hirshfeld_type(mixed_cdft%becke_control%cavity_env)
         CALL set_hirshfeld_info(mixed_cdft%becke_control%cavity_env, &
                                 shape_function_type=shape_function_gaussian, iterative=.FALSE., &
                                 radius_type=mixed_cdft%becke_control%cavity_shape, &
                                 use_bohr=mixed_cdft%becke_control%use_bohr)
      END IF
      mixed_cdft%first_iteration = .TRUE.
      ! MD related settings
      CALL force_env_get(force_env, root_section=root_section)
      md_section => section_vals_get_subs_vals(root_section, "MOTION%MD")
      CALL section_vals_val_get(md_section, "TIMESTEP", r_val=mixed_cdft%sim_dt)
      CALL section_vals_val_get(md_section, "STEP_START_VAL", i_val=mixed_cdft%sim_step)
      mixed_cdft%sim_step = mixed_cdft%sim_step-1 ! to get the first step correct
      mixed_cdft%sim_dt = cp_unit_from_cp2k(mixed_cdft%sim_dt, "fs")
      ! Some CDFT features are currently disabled for mixed calculations: check that these features were not requested
      IF (.NOT. mixed_cdft%becke_control%confine_dir == 3) &
         CALL cp_abort(__LOCATION__, &
                       "Only z-direction allowed for static/dynamic confinement with mixed CDFT")
      IF (mixed_cdft%becke_control%atomic_charges) &
         CALL cp_abort(__LOCATION__, &
                       "Calculation of atomic Becke charges not supported with mixed CDFT")
      IF (mixed_cdft%dlb .AND. .NOT. mixed_cdft%becke_control%cavity_confine) &
         CALL cp_abort(__LOCATION__, &
                       "Mixed CDFT load balancing requires Gaussian cavity confinement.")
      ! Decide what kind of layout to use and setup the grid
      ! Processor mappings currently supported:
      !  (2np,1)  --> (np,1)
      !  (nx,2ny) --> (nx,ny)
      !  (nx,ny)  --> (nx*ny/2,1) (required when xc_smooth is in use and with intermediate proc counts)
      !
      ! For cases 2 and 3, dlb redistributes YZ slices from overloaded processors to underloaded processors
      ! For case 1, XZ slices are redistributed
      ! These mappings could be unified. Now we essentially have separate code for cases 1-2 and 3.
      ! This leads to very messy code especially with dlb turned on...
      ! In terms of memory usage, it would be beneficial to replace case 1 with 3
      ! and implement a similar arbitrary mapping to replace case 2

      mixed_cdft%is_pencil = .FALSE. ! Flag to control the first two mappings
      mixed_cdft%is_special = .FALSE. ! Flag to control the last mapping
      ! With xc smoothing, the grid is always (ncpu/2,1) distributed
      ! and correct behavior cannot be guaranteed for ncpu/2 > nx, so we abort...
      IF (ncpu/2 .GT. npts(1, 1)) &
         CALL cp_abort(__LOCATION__, &
                       "ncpu/2 => nx: decrease ncpu or disable xc_smoothing")
      !
      ALLOCATE (mixed_rs_dims(2))
      IF (rs_dims(2, 1) /= 1) mixed_cdft%is_pencil = .TRUE.
      IF (.NOT. mixed_cdft%is_pencil .AND. ncpu .GT. npts(1, 1)) mixed_cdft%is_special = .TRUE.
      IF (mixed_cdft%is_special) THEN
         mixed_rs_dims = (/-1, -1/)
      ELSE IF (mixed_cdft%is_pencil) THEN
         mixed_rs_dims = (/rs_dims(1, 1), 2*rs_dims(2, 1)/)
      ELSE
         mixed_rs_dims = (/2*rs_dims(1, 1), 1/)
      END IF
      CALL pw_grid_setup(cell_mix%hmat, pw_grid, grid_span=grid_span(1), npts=npts(:, 1), &
                         cutoff=cutoff(1), spherical=is_spherical, odd=is_odd, &
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
         CALL cp_abort(__LOCATION__, &
                       "Unable to create a suitable grid distribution "// &
                       "for mixed CDFT calculations. Try decreasing the total number "// &
                       " of processors or disabling xc_smoothing.")
      ! Deallocate temporaries
      DEALLOCATE (mixed_rs_dims)
      DEALLOCATE (grid_span)
      DEALLOCATE (npts)
      DEALLOCATE (spherical)
      DEALLOCATE (rs_dims)
      DEALLOCATE (odd)
      DEALLOCATE (atoms)
      DEALLOCATE (coeffs)
      DEALLOCATE (cutoffs)
      DEALLOCATE (radii)
      DEALLOCATE (settings_int)
      DEALLOCATE (settings_real)
      DEALLOCATE (confine_bounds)
      ! Create the pool
      bo_mixed = pw_grid%bounds_local
      ALLOCATE (pw_pools(1))
      NULLIFY (pw_pools(1)%pool)
      CALL pw_pool_create(pw_pools(1)%pool, pw_grid=pw_grid)
      ! Gaussian confinement/wavefunction overlap method needs qs_kind_set
      ! Gaussian cavity confinement also needs the auxbas_rs_grid
      IF (mixed_cdft%becke_control%cavity_confine .OR. &
          mixed_cdft%wfn_overlap_method) THEN
         print_section => section_vals_get_subs_vals(force_env_section, &
                                                     "PRINT%GRID_INFORMATION")
         ALLOCATE (mixed_cdft%pw_env%gridlevel_info)
         CALL init_gaussian_gridlevel(mixed_cdft%pw_env%gridlevel_info, &
                                      ngrid_levels=1, cutoff=cutoff, &
                                      rel_cutoff=rel_cutoff(1), &
                                      print_section=print_section)
         ALLOCATE (rs_descs(1))
         ALLOCATE (rs_grids(1))
         ALLOCATE (mixed_cdft%pw_env%cube_info(1))
         higher_grid_layout = (/-1, -1, -1/)
         CALL init_d3_poly_module()
         CALL init_cube_info(mixed_cdft%pw_env%cube_info(1), &
                             pw_grid%dr(:), pw_grid%dh(:, :), &
                             pw_grid%dh_inv(:, :), &
                             pw_grid%orthorhombic, radius)
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
      DEALLOCATE (cutoff, rel_cutoff)
      CALL pw_grid_release(pw_grid)
      mixed_cdft%pw_env%auxbas_grid = 1
      NULLIFY (mixed_cdft%pw_env%pw_pools)
      mixed_cdft%pw_env%pw_pools => pw_pools
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
      ! Write information about the mixed CDFT calculation
      IF (iounit > 0) THEN
         WRITE (iounit, *) ""
         WRITE (iounit, FMT="(T2,A,T71)") "MIXED_CDFT| Activating mixed CDFT calculation"
         WRITE (iounit, FMT="(T2,A,T71)") "MIXED_CDFT| Becke constraint is first built using all available processors"
         WRITE (iounit, FMT="(T2,A,T71)") "MIXED_CDFT| and then copied to both states with their own processor groups"
         WRITE (iounit, FMT="(T2,A,T71,L10)") &
            "MIXED_CDFT| Calculating electronic coupling between states:         ", mixed_env%do_mixed_et
         WRITE (iounit, FMT="(T2,A,T71,L10)") &
            "MIXED_CDFT| Dynamic load balancing enabled:                         ", mixed_cdft%dlb
         IF (mixed_cdft%dlb) THEN
            WRITE (iounit, FMT="(T2,A,T71)") "MIXED_CDFT| Dynamic load balancing parameters:                         "
            WRITE (iounit, FMT="(T2,A,T71,F10.2)") &
               "MIXED_CDFT| load_scale:", mixed_cdft%dlb_control%load_scale
            WRITE (iounit, FMT="(T2,A,T71,F10.2)") &
               "MIXED_CDFT| very_overloaded:", mixed_cdft%dlb_control%very_overloaded
            WRITE (iounit, FMT="(T2,A,T71,I10)") &
               "MIXED_CDFT| more_work:", mixed_cdft%dlb_control%more_work
         END IF
         IF (mixed_env%do_mixed_et) THEN
            IF (mixed_cdft%eps_svd == 0.0_dp) THEN
               WRITE (iounit, FMT="(T2,A,T71)") "MIXED_CDFT| Matrix inversions calculated with LU decomposition."
            ELSE
               WRITE (iounit, FMT="(T2,A,T71)") "MIXED_CDFT| Matrix inversions calculated with SVD decomposition."
               WRITE (iounit, FMT="(T2,A,T71,ES10.2)") "MIXED_CDFT| EPS_SVD:", mixed_cdft%eps_svd
            END IF
         END IF
      END IF
      CALL set_mixed_env(mixed_env, cdft_control=mixed_cdft)
      CALL cp_print_key_finished_output(iounit, logger, force_env_section, &
                                        "MIXED%PRINT%PROGRAM_RUN_INFO")
      CALL timestop(handle)

   END SUBROUTINE mixed_cdft_init

! **************************************************************************************************
!> \brief Build CDFT weight and gradient in parallel on all available processors
!> \param force_env the force_env that holds the CDFT states
!> \param calculate_forces determines if forces should be calculted
!> \par History
!>       01.2016  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_build_weight(force_env, calculate_forces)
      TYPE(force_env_type), POINTER                      :: force_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_build_weight', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, handle2, i, iforce_eval, ind, INDEX(6), iounit, j, lb_min, &
         my_special_work, natom, nforce_eval, recv_offset, req(3), ub_max
      INTEGER, DIMENSION(2, 3)                           :: bo
      INTEGER, DIMENSION(:), POINTER                     :: lb, req_total, sendbuffer_i, ub
      REAL(KIND=dp)                                      :: t1, t2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys_mix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(force_env_type), POINTER                      :: force_env_qs
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(mixed_environment_type), POINTER              :: mixed_env
      TYPE(particle_list_type), POINTER                  :: particles_mix
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool, mixed_auxbas_pw_pool
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: force_env_section, print_section

      TYPE buffers
         INTEGER                               :: imap(6)
         INTEGER, DIMENSION(:), &
            POINTER                             :: iv => null()
         LOGICAL, DIMENSION(:, :, :), &
            POINTER                             :: b3 => null()
         REAL(KIND=dp), POINTER, &
            DIMENSION(:, :, :)                    :: r3 => null()
         REAL(KIND=dp), POINTER, &
            DIMENSION(:, :, :, :)                  :: r4 => null()
      END TYPE buffers
      TYPE(buffers), DIMENSION(:), POINTER     :: recvbuffer

      NULLIFY (subsys_mix, force_env_qs, particles_mix, force_env_section, print_section, &
               mixed_env, mixed_cdft, pw_env, auxbas_pw_pool, mixed_auxbas_pw_pool, &
               qs_env, dft_control, sendbuffer_i, lb, ub, req_total, recvbuffer)

      logger => cp_get_default_logger()
      CPASSERT(ASSOCIATED(force_env))
      nforce_eval = SIZE(force_env%sub_force_env)
      CALL timeset(routineN, handle)
      t1 = m_walltime()
      ! Get infos about the mixed subsys
      CALL force_env_get(force_env=force_env, &
                         subsys=subsys_mix, &
                         force_env_section=force_env_section)
      CALL cp_subsys_get(subsys=subsys_mix, &
                         particles=particles_mix)
      DO iforce_eval = 1, nforce_eval
         IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
         SELECT CASE (force_env%sub_force_env (iforce_eval)%force_env%in_use)
         CASE (use_qs_force)
            force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
         CASE (use_qmmm)
            force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
         CASE DEFAULT
            CPASSERT(.FALSE.)
         END SELECT
      END DO
      IF (.NOT. force_env%mixed_env%do_mixed_qmmm_cdft) THEN
         CALL force_env_get(force_env=force_env_qs, &
                            qs_env=qs_env, &
                            subsys=subsys_mix)
         CALL cp_subsys_get(subsys=subsys_mix, &
                            particles=particles_mix)
      ELSE
         qs_env => force_env_qs%qmmm_env%qs_env
         CALL get_qs_env(qs_env, cp_subsys=subsys_mix)
         CALL cp_subsys_get(subsys=subsys_mix, &
                            particles=particles_mix)
      END IF
      mixed_env => force_env%mixed_env
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%PRINT%PROGRAM_RUN_INFO")
      iounit = cp_print_key_unit_nr(logger, print_section, '', extension='.mixedLog')
      CALL get_mixed_env(mixed_env, cdft_control=mixed_cdft)
      CPASSERT(ASSOCIATED(mixed_cdft))
      ! Calculate the Becke weight function and gradient on all np processors
      CALL pw_env_get(pw_env=mixed_cdft%pw_env, auxbas_pw_pool=mixed_auxbas_pw_pool)
      natom = SIZE(particles_mix%els)
      CALL mixed_becke_constraint(force_env, calculate_forces)
      ! Start replicating the working arrays on both np/2 processor groups
      mixed_cdft%sim_step = mixed_cdft%sim_step+1
      CALL get_qs_env(qs_env, pw_env=pw_env, dft_control=dft_control)
      CPASSERT(dft_control%qs_control%becke_restraint)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      bo = auxbas_pw_pool%pw_grid%bounds_local
      ! First determine the size of the arrays along the confinement dir
      IF (mixed_cdft%is_special) THEN
         my_special_work = 2
      ELSE
         my_special_work = 1
      END IF
      ALLOCATE (recvbuffer(SIZE(mixed_cdft%source_list)))
      ALLOCATE (req_total(my_special_work*SIZE(mixed_cdft%dest_list)+SIZE(mixed_cdft%source_list)))
      ALLOCATE (lb(SIZE(mixed_cdft%source_list)), ub(SIZE(mixed_cdft%source_list)))
      IF (mixed_cdft%becke_control%cavity_confine .AND. &
          .NOT. mixed_cdft%becke_control%dynamic_confine) THEN
         ! Only Gaussian confinement => the bounds depend on the processor and need to be communicated
         ALLOCATE (sendbuffer_i(2))
         IF (mixed_cdft%becke_control%confine_bounds_int(2) .LT. &
             mixed_auxbas_pw_pool%pw_grid%bounds_local(2, 3)) THEN
            sendbuffer_i(2) = mixed_cdft%becke_control%confine_bounds_int(2)-1
         ELSE
            sendbuffer_i(2) = mixed_cdft%becke_control%confine_bounds_int(2)
         END IF
         IF (mixed_cdft%becke_control%confine_bounds_int(1) .GT. &
             mixed_auxbas_pw_pool%pw_grid%bounds_local(1, 3)) THEN
            sendbuffer_i(1) = mixed_cdft%becke_control%confine_bounds_int(1)+1
         ELSE
            sendbuffer_i(1) = mixed_cdft%becke_control%confine_bounds_int(1)
         END IF
         DO i = 1, SIZE(mixed_cdft%source_list)
            ALLOCATE (recvbuffer(i)%iv(2))
            CALL mp_irecv(msgout=recvbuffer(i)%iv, source=mixed_cdft%source_list(i), &
                          request=req_total(i), &
                          comm=force_env%para_env%group)
         END DO
         DO i = 1, my_special_work
            DO j = 1, SIZE(mixed_cdft%dest_list)
               ind = j+(i-1)*SIZE(mixed_cdft%dest_list)+SIZE(mixed_cdft%source_list)
               CALL mp_isend(msgin=sendbuffer_i, &
                             dest=mixed_cdft%dest_list(j)+(i-1)*force_env%para_env%num_pe/2, &
                             request=req_total(ind), &
                             comm=force_env%para_env%group)
            END DO
         END DO
         CALL mp_waitall(req_total)
         ! Find the largest/smallest value of ub/lb
         DEALLOCATE (sendbuffer_i)
         lb_min = HUGE(0)
         ub_max = -HUGE(0)
         DO i = 1, SIZE(mixed_cdft%source_list)
            lb(i) = recvbuffer(i)%iv(1)
            ub(i) = recvbuffer(i)%iv(2)
            IF (lb(i) .LT. lb_min) lb_min = lb(i)
            IF (ub(i) .GT. ub_max) ub_max = ub(i)
            DEALLOCATE (recvbuffer(i)%iv)
         END DO
         ! Take into account the grids already communicated during dlb
         IF (mixed_cdft%dlb) THEN
            IF (ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
               DO j = 1, SIZE(mixed_cdft%dlb_control%recv_work_repl)
                  IF (mixed_cdft%dlb_control%recv_work_repl(j)) THEN
                     DO i = 1, SIZE(mixed_cdft%dlb_control%recvbuff(j)%buffs)
                        IF (LBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight, 3) &
                            .LT. lb_min) lb_min = LBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight, 3)
                        IF (UBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight, 3) &
                            .GT. ub_max) ub_max = UBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight, 3)
                     END DO
                  END IF
               END DO
            END IF
         END IF
      ELSE IF (mixed_cdft%becke_control%dynamic_confine .OR. &
               mixed_cdft%becke_control%confine) THEN
         ! No Gaussian confinement or mixed confinement => the bounds are the same on all processors
         ub_max = mixed_cdft%becke_control%confine_bounds_int(2)-1
         lb_min = mixed_cdft%becke_control%confine_bounds_int(1)+1
         lb = lb_min
         ub = ub_max
      ELSE
         ! No confinement
         ub_max = bo(2, 3)
         lb_min = bo(1, 3)
         lb = lb_min
         ub = ub_max
      END IF
      ! Determine the sender specific indices of grid slices that are to be received
      CALL timeset(routineN//"_comm", handle2)
      DO j = 1, SIZE(recvbuffer)
         ind = j+(j/2)
         IF (mixed_cdft%is_special) THEN
            recvbuffer(j)%imap = (/mixed_cdft%source_list_bo(1, j), mixed_cdft%source_list_bo(2, j), &
                                   mixed_cdft%source_list_bo(3, j), mixed_cdft%source_list_bo(4, j), &
                                   lb(j), ub(j)/)
         ELSE IF (mixed_cdft%is_pencil) THEN
            recvbuffer(j)%imap = (/bo(1, 1), bo(2, 1), mixed_cdft%recv_bo(ind), mixed_cdft%recv_bo(ind+1), lb(j), ub(j)/)
         ELSE
            recvbuffer(j)%imap = (/mixed_cdft%recv_bo(ind), mixed_cdft%recv_bo(ind+1), bo(1, 2), bo(2, 2), lb(j), ub(j)/)
         END IF
      END DO
      IF (mixed_cdft%dlb .AND. .NOT. mixed_cdft%is_special) THEN
         IF (mixed_cdft%dlb_control%recv_work_repl(1) .OR. mixed_cdft%dlb_control%recv_work_repl(2)) THEN
            DO j = 1, 2
               recv_offset = 0
               IF (mixed_cdft%dlb_control%recv_work_repl(j)) &
                  recv_offset = SUM(mixed_cdft%dlb_control%recv_info(j)%target_list(2, :))
               IF (mixed_cdft%is_pencil) THEN
                  recvbuffer(j)%imap(1) = recvbuffer(j)%imap(1)+recv_offset
               ELSE
                  recvbuffer(j)%imap(3) = recvbuffer(j)%imap(3)+recv_offset
               END IF
            END DO
         END IF
      END IF
      ! Transfer the arrays one-by-one and deallocate shared storage
      ! Start with  theweight function
      DO j = 1, SIZE(mixed_cdft%source_list)
         ALLOCATE (recvbuffer(j)%r3(recvbuffer(j)%imap(1):recvbuffer(j)%imap(2), &
                                    recvbuffer(j)%imap(3):recvbuffer(j)%imap(4), &
                                    recvbuffer(j)%imap(5):recvbuffer(j)%imap(6)))

         CALL mp_irecv(msgout=recvbuffer(j)%r3, source=mixed_cdft%source_list(j), &
                       request=req_total(j), comm=force_env%para_env%group)
      END DO
      DO i = 1, my_special_work
         DO j = 1, SIZE(mixed_cdft%dest_list)
            ind = j+(i-1)*SIZE(mixed_cdft%dest_list)+SIZE(mixed_cdft%source_list)
            IF (mixed_cdft%is_special) THEN
               CALL mp_isend(msgin=mixed_cdft%sendbuff(j)%weight, &
                             dest=mixed_cdft%dest_list(j)+(i-1)*force_env%para_env%num_pe/2, &
                             request=req_total(ind), &
                             comm=force_env%para_env%group)
            ELSE
               CALL mp_isend(msgin=mixed_cdft%weight, dest=mixed_cdft%dest_list(j), &
                             request=req_total(ind), comm=force_env%para_env%group)
            END IF
         END DO
      END DO
      CALL mp_waitall(req_total)
      IF (mixed_cdft%is_special) THEN
         DO j = 1, SIZE(mixed_cdft%dest_list)
            DEALLOCATE (mixed_cdft%sendbuff(j)%weight)
         END DO
      ELSE
         DEALLOCATE (mixed_cdft%weight)
      END IF
      ! In principle, we could reduce the memory footprint of becke_pot by only transfering the nonzero portion
      ! of the potential, but this would require a custom integrate_v_rspace
      CALL pw_pool_create_pw(auxbas_pw_pool, dft_control%qs_control%becke_control%becke_pot%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      dft_control%qs_control%becke_control%becke_pot%pw%cr3d = 0.0_dp
      ! Assemble the recved slices
      DO j = 1, SIZE(mixed_cdft%source_list)
         dft_control%qs_control%becke_control%becke_pot%pw%cr3d(recvbuffer(j)%imap(1):recvbuffer(j)%imap(2), &
                                                                recvbuffer(j)%imap(3):recvbuffer(j)%imap(4), &
                                                                recvbuffer(j)%imap(5):recvbuffer(j)%imap(6)) = &
            recvbuffer(j)%r3
      END DO
      ! Do the same for slices sent during dlb
      IF (mixed_cdft%dlb) THEN
         IF (ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
            DO j = 1, SIZE(mixed_cdft%dlb_control%recv_work_repl)
               IF (mixed_cdft%dlb_control%recv_work_repl(j)) THEN
                  DO i = 1, SIZE(mixed_cdft%dlb_control%recvbuff(j)%buffs)
                     index = (/LBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight, 1), &
                               UBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight, 1), &
                               LBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight, 2), &
                               UBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight, 2), &
                               LBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight, 3), &
                               UBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight, 3)/)
                     dft_control%qs_control%becke_control%becke_pot%pw%cr3d(INDEX(1):INDEX(2), &
                                                                            INDEX(3):INDEX(4), &
                                                                            INDEX(5):INDEX(6)) = &
                        mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight
                     DEALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight)
                  END DO
               END IF
            END DO
         END IF
      END IF
      ! Gaussian confinement cavity
      IF (mixed_cdft%becke_control%cavity_confine) THEN
         DO j = 1, SIZE(mixed_cdft%source_list)
            CALL mp_irecv(msgout=recvbuffer(j)%r3, source=mixed_cdft%source_list(j), &
                          request=req_total(j), comm=force_env%para_env%group)
         END DO
         DO i = 1, my_special_work
            DO j = 1, SIZE(mixed_cdft%dest_list)
               ind = j+(i-1)*SIZE(mixed_cdft%dest_list)+SIZE(mixed_cdft%source_list)
               IF (mixed_cdft%is_special) THEN
                  CALL mp_isend(msgin=mixed_cdft%sendbuff(j)%cavity, &
                                dest=mixed_cdft%dest_list(j)+(i-1)*force_env%para_env%num_pe/2, &
                                request=req_total(ind), &
                                comm=force_env%para_env%group)
               ELSE
                  CALL mp_isend(msgin=mixed_cdft%cavity, dest=mixed_cdft%dest_list(j), &
                                request=req_total(ind), comm=force_env%para_env%group)
               END IF
            END DO
         END DO
         CALL mp_waitall(req_total)
         IF (mixed_cdft%is_special) THEN
            DO j = 1, SIZE(mixed_cdft%dest_list)
               DEALLOCATE (mixed_cdft%sendbuff(j)%cavity)
            END DO
         ELSE
            DEALLOCATE (mixed_cdft%cavity)
         END IF
         ! We only need the nonzero part of the confinement cavity
         ALLOCATE (dft_control%qs_control%becke_control%cavity_mat(bo(1, 1):bo(2, 1), &
                                                                   bo(1, 2):bo(2, 2), &
                                                                   lb_min:ub_max))
         dft_control%qs_control%becke_control%cavity_mat = 0.0_dp

         DO j = 1, SIZE(mixed_cdft%source_list)
            dft_control%qs_control%becke_control%cavity_mat(recvbuffer(j)%imap(1):recvbuffer(j)%imap(2), &
                                                            recvbuffer(j)%imap(3):recvbuffer(j)%imap(4), &
                                                            recvbuffer(j)%imap(5):recvbuffer(j)%imap(6)) = recvbuffer(j)%r3
         END DO
         IF (mixed_cdft%dlb) THEN
            IF (ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
               DO j = 1, SIZE(mixed_cdft%dlb_control%recv_work_repl)
                  IF (mixed_cdft%dlb_control%recv_work_repl(j)) THEN
                     DO i = 1, SIZE(mixed_cdft%dlb_control%recvbuff(j)%buffs)
                        index = (/LBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%cavity, 1), &
                                  UBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%cavity, 1), &
                                  LBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%cavity, 2), &
                                  UBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%cavity, 2), &
                                  LBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%cavity, 3), &
                                  UBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%cavity, 3)/)
                        dft_control%qs_control%becke_control%cavity_mat(INDEX(1):INDEX(2), &
                                                                        INDEX(3):INDEX(4), &
                                                                        INDEX(5):INDEX(6)) = &
                           mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%cavity
                        DEALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%cavity)
                     END DO
                  END IF
               END DO
            END IF
         END IF
      END IF
      DO j = 1, SIZE(mixed_cdft%source_list)
         DEALLOCATE (recvbuffer(j)%r3)
      END DO
      IF (calculate_forces) THEN
         ! List of skipped grid points for other confiment methods
         IF (mixed_cdft%becke_control%should_skip .AND. .NOT. &
             mixed_cdft%becke_control%cavity_confine) THEN
            CALL mp_isend(msgin=mixed_cdft%becke_control%skip_list, dest=mixed_cdft%dest_list(1), &
                          request=req(1), comm=force_env%para_env%group)
            DO j = 1, 2
               ALLOCATE (recvbuffer(j)%b3(recvbuffer(j)%imap(1):recvbuffer(j)%imap(2), &
                                          recvbuffer(j)%imap(3):recvbuffer(j)%imap(4), &
                                          recvbuffer(j)%imap(5):recvbuffer(j)%imap(6)))
               CALL mp_irecv(msgout=recvbuffer(j)%b3, source=mixed_cdft%source_list(j), &
                             request=req(j+1), comm=force_env%para_env%group)
            END DO
            CALL mp_wait(req(1))
            CALL mp_isend(msgin=mixed_cdft%becke_control%skip_list, dest=mixed_cdft%dest_list(2), &
                          request=req(1), comm=force_env%para_env%group)
            CALL mp_waitall(req)
            DEALLOCATE (mixed_cdft%becke_control%skip_list)
            ALLOCATE (dft_control%qs_control%becke_control%skip_list(bo(1, 1):bo(2, 1), &
                                                                     bo(1, 2):bo(2, 2), &
                                                                     lb_min:ub_max))
            DO j = 1, 2
               dft_control%qs_control%becke_control%skip_list(recvbuffer(j)%imap(1):recvbuffer(j)%imap(2), &
                                                              recvbuffer(j)%imap(3):recvbuffer(j)%imap(4), &
                                                              recvbuffer(j)%imap(5):recvbuffer(j)%imap(6)) = &
                  recvbuffer(j)%b3
               DEALLOCATE (recvbuffer(j)%b3)
            END DO
         END IF
         ! Gradients of the weight function
         DO j = 1, SIZE(mixed_cdft%source_list)
            ALLOCATE (recvbuffer(j)%r4(3*natom, recvbuffer(j)%imap(1):recvbuffer(j)%imap(2), &
                                       recvbuffer(j)%imap(3):recvbuffer(j)%imap(4), &
                                       recvbuffer(j)%imap(5):recvbuffer(j)%imap(6)))
            CALL mp_irecv(msgout=recvbuffer(j)%r4, source=mixed_cdft%source_list(j), &
                          request=req_total(j), comm=force_env%para_env%group)
         END DO
         DO i = 1, my_special_work
            DO j = 1, SIZE(mixed_cdft%dest_list)
               ind = j+(i-1)*SIZE(mixed_cdft%dest_list)+SIZE(mixed_cdft%source_list)
               IF (mixed_cdft%is_special) THEN
                  CALL mp_isend(msgin=mixed_cdft%sendbuff(j)%gradients, &
                                dest=mixed_cdft%dest_list(j)+(i-1)*force_env%para_env%num_pe/2, &
                                request=req_total(ind), &
                                comm=force_env%para_env%group)
               ELSE
                  CALL mp_isend(msgin=mixed_cdft%becke_control%gradients, dest=mixed_cdft%dest_list(j), &
                                request=req_total(ind), comm=force_env%para_env%group)
               END IF
            END DO
         END DO
         CALL mp_waitall(req_total)
         IF (mixed_cdft%is_special) THEN
            DO j = 1, SIZE(mixed_cdft%dest_list)
               DEALLOCATE (mixed_cdft%sendbuff(j)%gradients)
            END DO
            DEALLOCATE (mixed_cdft%sendbuff)
         ELSE
            DEALLOCATE (mixed_cdft%becke_control%gradients)
         END IF
         ALLOCATE (dft_control%qs_control%becke_control%gradients(3*natom, bo(1, 1):bo(2, 1), &
                                                                  bo(1, 2):bo(2, 2), lb_min:ub_max))
         DO j = 1, SIZE(mixed_cdft%source_list)
            dft_control%qs_control%becke_control%gradients(:, recvbuffer(j)%imap(1):recvbuffer(j)%imap(2), &
                                                           recvbuffer(j)%imap(3):recvbuffer(j)%imap(4), &
                                                           recvbuffer(j)%imap(5):recvbuffer(j)%imap(6)) = &
               recvbuffer(j)%r4
            DEALLOCATE (recvbuffer(j)%r4)
         END DO
         IF (mixed_cdft%dlb) THEN
            IF (ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
               DO j = 1, SIZE(mixed_cdft%dlb_control%recv_work_repl)
                  IF (mixed_cdft%dlb_control%recv_work_repl(j)) THEN
                     DO i = 1, SIZE(mixed_cdft%dlb_control%recvbuff(j)%buffs)
                        index = (/LBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%gradients, 2), &
                                  UBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%gradients, 2), &
                                  LBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%gradients, 3), &
                                  UBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%gradients, 3), &
                                  LBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%gradients, 4), &
                                  UBOUND(mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%gradients, 4)/)
                        dft_control%qs_control%becke_control%gradients(:, INDEX(1):INDEX(2), &
                                                                       INDEX(3):INDEX(4), &
                                                                       INDEX(5):INDEX(6)) = &
                           mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%gradients
                        DEALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%gradients)
                     END DO
                  END IF
               END DO
            END IF
         END IF
      END IF
      ! Clean up remaining temporaries
      IF (mixed_cdft%dlb) THEN
         IF (ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
            DO j = 1, SIZE(mixed_cdft%dlb_control%recv_work_repl)
               IF (mixed_cdft%dlb_control%recv_work_repl(j)) THEN
                  IF (ASSOCIATED(mixed_cdft%dlb_control%recv_info(j)%target_list)) &
                     DEALLOCATE (mixed_cdft%dlb_control%recv_info(j)%target_list)
                  DEALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs)
               END IF
            END DO
            DEALLOCATE (mixed_cdft%dlb_control%recv_info, mixed_cdft%dlb_control%recvbuff)
         END IF
         IF (ASSOCIATED(mixed_cdft%dlb_control%target_list)) &
            DEALLOCATE (mixed_cdft%dlb_control%target_list)
         DEALLOCATE (mixed_cdft%dlb_control%recv_work_repl)
      END IF
      DEALLOCATE (recvbuffer)
      DEALLOCATE (req_total)
      DEALLOCATE (lb)
      DEALLOCATE (ub)
      CALL timestop(handle2)
      ! Set some flags so the weight is not rebuilt during SCF
      dft_control%qs_control%becke_control%external_control = .TRUE.
      dft_control%qs_control%becke_control%need_pot = .FALSE.
      ! Store the bound indices for force calculation
      IF (calculate_forces) THEN
         ! Note: these are the final bounds i.e. the first nonzero elements of the arrays...
         ! On the contrary, mixed_cdft%....%bounds include (possibly) the zero elements
         ! This is messy and needs to be fixed
         dft_control%qs_control%becke_control%confine_bounds_int(2) = ub_max
         dft_control%qs_control%becke_control%confine_bounds_int(1) = lb_min
      END IF
      CALL pw_scale(dft_control%qs_control%becke_control%becke_pot%pw, &
                    dft_control%qs_control%becke_control%becke_pot%pw%pw_grid%dvol)
      ! Set flags for ET coupling calculation
      IF (mixed_env%do_mixed_et) THEN
         IF (MODULO(mixed_cdft%sim_step, mixed_env%et_freq) == 0) THEN
            dft_control%qs_control%cdft_control%do_et = .TRUE.
            dft_control%qs_control%cdft_control%calculate_metric = mixed_cdft%calculate_metric
         END IF
      END IF
      t2 = m_walltime()
      IF (iounit > 0) THEN
         WRITE (iounit, '(A)') ' '
         WRITE (iounit, '(T2,A,F6.1,A)') 'MIXED_CDFT| Becke constraint built in ', t2-t1, ' seconds'
         WRITE (iounit, '(A)') ' '
      END IF
      CALL cp_print_key_finished_output(iounit, logger, force_env_section, &
                                        "MIXED%PRINT%PROGRAM_RUN_INFO")
      CALL timestop(handle)

   END SUBROUTINE mixed_cdft_build_weight

! **************************************************************************************************
!> \brief ET coupling adapted to mixed calculations
!> \param force_env the force_env that holds the CDFT states
!> \par History
!>       02.15  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_calculate_coupling(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_calculate_coupling', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: file_name
      INTEGER :: handle, iforce_eval, iounit, ispin, j, k, my_id, nao, ncol_local, ncol_overlap, &
         ncol_wmat, nforce_eval, nmo, nrow_local, nrow_overlap, nrow_wmat, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ncol_mo, nrow_mo
      LOGICAL                                            :: exist, natom_mismatch
      REAL(KIND=dp)                                      :: coupling_wfn, energy_a, energy_b, &
                                                            maxocc, Sda, strength_a, strength_b, &
                                                            Waa, Wbb, Wda
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: a, b, evals, metric, S_det
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: overlaps
      REAL(KIND=dp), DIMENSION(2)                        :: c, eigenv
      REAL(KIND=dp), DIMENSION(2, 2)                     :: S_mat, tmp_mat, U, W_mat
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mixed_wmat_tmp, rest_MO, wmat_tmp
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: matrix_p_tmp, mixed_matrix_p_tmp, &
                                                            mixed_mo_coeff, mo_coeff_tmp
      TYPE(cp_fm_p_type), POINTER                        :: matrix_s_tmp, mixed_matrix_s_tmp
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_mo, fm_struct_overlap, &
                                                            fm_struct_tmp, fm_struct_wmat, &
                                                            mo_mo_fmstruct
      TYPE(cp_fm_type), POINTER                          :: inverse_mat, mo_overlap, mo_tmp, SMO, &
                                                            Tinverse, tmp2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys_mix
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: density_matrix
      TYPE(dbcsr_type)                                   :: desymm_tmp, e_vectors, mixed_wmat_a, &
                                                            mixed_wmat_b
      TYPE(dbcsr_type), POINTER                          :: mixed_matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(mixed_environment_type), POINTER              :: mixed_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: force_env_section, linear_section, &
                                                            print_section

      NULLIFY (subsys_mix, force_env_section, print_section, &
               linear_section, mixed_env, mixed_cdft, qs_env, dft_control, &
               fm_struct_mo, fm_struct_wmat, fm_struct_overlap, fm_struct_tmp, &
               mo_mo_fmstruct, inverse_mat, SMO, Tinverse, tmp2, mo_overlap, &
               mo_tmp, mixed_matrix_s_tmp, matrix_s_tmp, mixed_wmat_tmp, &
               wmat_tmp, rest_MO, mixed_mo_coeff, mo_coeff_tmp, matrix_p_tmp, &
               mixed_matrix_p_tmp, mixed_matrix_s, density_matrix, blacs_env, &
               energy, mo_set, particle_set, qs_kind_set, atomic_kind_set)

      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(force_env))
      CALL timeset(routineN, handle)
      CALL force_env_get(force_env=force_env, &
                         force_env_section=force_env_section)
      mixed_env => force_env%mixed_env
      nforce_eval = SIZE(force_env%sub_force_env)
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%PRINT%PROGRAM_RUN_INFO")
      iounit = cp_print_key_unit_nr(logger, print_section, '', extension='.mixedLog')
      CALL get_mixed_env(mixed_env, cdft_control=mixed_cdft)
      CPASSERT(ASSOCIATED(mixed_cdft))
      ! Get qs_env and nspins
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
         CALL get_qs_env(qs_env, blacs_env=blacs_env)
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
      ALLOCATE (mixed_mo_coeff(2, nspins))
      ALLOCATE (mixed_wmat_tmp(2), mixed_matrix_s_tmp)
      IF (mixed_cdft%calculate_metric) THEN
         ALLOCATE (density_matrix(2, nspins), mixed_matrix_p_tmp(2, nspins))
      END IF
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
               ALLOCATE (density_matrix(iforce_eval, ispin)%matrix)
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
      CALL dbcsr_init_p(mixed_matrix_s)
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
      ! Begin coupling calculation
      ! Create the MO_MO fm struct (mo_mo_fm_pools%struct)
      ! The number of MOs/AOs is equal to the number of columns/rows of mo_coeff(:,:)%matrix

      ! TODO: Handle magnetization/combined constraints
      my_id = -1
      !my_id=dft_control%qs_control%becke_control%density_type
      ALLOCATE (rest_MO(2), a(nspins), b(nspins), S_det(nspins))
      DO ispin = 1, nspins
         NULLIFY (fm_struct_mo, mo_mo_fmstruct)
         nao = nrow_mo(ispin)
         nmo = ncol_mo(ispin)
         CALL cp_fm_struct_create(fm_struct_mo, nrow_global=nrow_mo(ispin), ncol_global=ncol_mo(ispin), &
                                  context=mixed_cdft%blacs_env, para_env=force_env%para_env)
         CALL cp_fm_struct_create(mo_mo_fmstruct, nrow_global=nmo, ncol_global=nmo, &
                                  context=mixed_cdft%blacs_env, para_env=force_env%para_env)
         CALL cp_fm_create(matrix=tmp2, matrix_struct=fm_struct_mo, &
                           name="ET_TMP"//TRIM(ADJUSTL(cp_to_string(2)))//"MATRIX")
         CALL cp_fm_struct_release(fm_struct_mo)
         CALL cp_fm_create(matrix=inverse_mat, matrix_struct=mo_mo_fmstruct, &
                           name="INVERSE"//TRIM(ADJUSTL(cp_to_string(2)))//"MATRIX")
         CALL cp_fm_create(matrix=Tinverse, matrix_struct=mo_mo_fmstruct, &
                           name="T_INVERSE"//TRIM(ADJUSTL(cp_to_string(2)))//"MATRIX")
         CALL cp_fm_create(matrix=SMO, matrix_struct=mo_mo_fmstruct, &
                           name="ET_SMO"//TRIM(ADJUSTL(cp_to_string(1)))//"MATRIX")
         DO j = 1, 2
            NULLIFY (rest_MO(j)%matrix)
            CALL cp_fm_create(matrix=rest_MO(j)%matrix, matrix_struct=mo_mo_fmstruct, &
                              name="ET_rest_MO"//TRIM(ADJUSTL(cp_to_string(j)))//"MATRIX")
         END DO
         CALL cp_fm_struct_release(mo_mo_fmstruct)
         ! calculate the MO overlap  (C_d)^T S C_a
         CALL cp_dbcsr_sm_fm_multiply(mixed_matrix_s, mixed_mo_coeff(1, ispin)%matrix, &
                                      tmp2, nmo, 1.0_dp, 0.0_dp)
         CALL cp_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, &
                      mixed_mo_coeff(2, ispin)%matrix, &
                      tmp2, 0.0_dp, SMO)
         ! calculate the MO-representation of the restraint matrix A (C_d)^T W_a C_a
         CALL cp_dbcsr_sm_fm_multiply(mixed_wmat_a, &
                                      mixed_mo_coeff(1, ispin)%matrix, &
                                      tmp2, nmo, 1.0_dp, 0.0_dp)
         CALL cp_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, &
                      mixed_mo_coeff(2, ispin)%matrix, &
                      tmp2, 0.0_dp, rest_MO(1)%matrix)
         ! calculate the MO-representation of the restraint matrix D (C_a)^T W_d C_d
         CALL cp_dbcsr_sm_fm_multiply(mixed_wmat_b, &
                                      mixed_mo_coeff(2, ispin)%matrix, &
                                      tmp2, nmo, 1.0_dp, 0.0_dp)
         CALL cp_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, &
                      mixed_mo_coeff(1, ispin)%matrix, &
                      tmp2, 0.0_dp, rest_MO(2)%matrix)
         ! Invert and calculate determinant
         CALL cp_fm_invert(SMO, inverse_mat, S_det(ispin), eps_svd=mixed_cdft%eps_svd)
         CALL cp_fm_get_info(inverse_mat, nrow_local=nrow_local, ncol_local=ncol_local)
         !
         b(ispin) = 0.0_dp
         DO j = 1, ncol_local
            DO k = 1, nrow_local
               b(ispin) = b(ispin)+rest_MO(2)%matrix%local_data(k, j)*inverse_mat%local_data(k, j)
            END DO
         END DO
         !
         CALL cp_fm_transpose(inverse_mat, Tinverse)
         a(ispin) = 0.0_dp
         DO j = 1, ncol_local
            DO k = 1, nrow_local
               a(ispin) = a(ispin)+rest_MO(1)%matrix%local_data(k, j)*Tinverse%local_data(k, j)
            END DO
         END DO
         IF ((my_id == do_spin_density) .AND. ispin == 2) THEN
            a(ispin) = -a(ispin)
            b(ispin) = -b(ispin)
         END IF
         ! Compute density matrix difference P = P_D - P_A
         IF (mixed_cdft%calculate_metric) THEN
            CALL dbcsr_add(density_matrix(1, ispin)%matrix, density_matrix(2, ispin)%matrix, -1.0_dp, 1.0_dp)
            CALL dbcsr_release_p(density_matrix(2, ispin)%matrix)
         END IF
         !
         CALL mp_sum(a(ispin), force_env%para_env%group)
         CALL mp_sum(b(ispin), force_env%para_env%group)
         ! Release work
         CALL cp_fm_release(tmp2)
         CALL cp_fm_release(rest_MO(1)%matrix)
         CALL cp_fm_release(rest_MO(2)%matrix)
         CALL cp_fm_release(SMO)
         CALL cp_fm_release(Tinverse)
         CALL cp_fm_release(inverse_mat)
      END DO
      ! solve eigenstates for the projector matrix
      IF (nspins == 2) THEN
         Sda = ABS(S_det(1)*S_det(2)) ! The sign might be wrong on some processors (used to be?)
         Wda = ((a(1)+a(2))+(b(1)+b(2)))*0.5_dp*Sda
      ELSE
         Sda = S_det(1)**2
         Wda = (a(1)+b(1))*Sda
      END IF
      ! Transfer info about the constraint calculations
      Waa = 0.0_dp
      Wbb = 0.0_dp
      strength_a = mixed_env%strength(1)
      strength_b = mixed_env%strength(2)
      energy_a = 0.0_dp
      energy_b = 0.0_dp
      DO iforce_eval = 1, nforce_eval
         IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
         CALL get_qs_env(qs_env, energy=energy)
         IF (force_env%sub_force_env(iforce_eval)%force_env%para_env%mepos == &
             force_env%sub_force_env(iforce_eval)%force_env%para_env%source) THEN
            IF (iforce_eval == 1) THEN
               Waa = dft_control%qs_control%cdft_control%value(1)
               energy_a = energy%total
            ELSE
               Wbb = dft_control%qs_control%cdft_control%value(1)
               energy_b = energy%total
            END IF
         END IF
      END DO
      CALL mp_sum(Waa, force_env%para_env%group)
      CALL mp_sum(Wbb, force_env%para_env%group)
      CALL mp_sum(energy_a, force_env%para_env%group)
      CALL mp_sum(energy_b, force_env%para_env%group)
      ! Compute metric to assess reliability of coupling
      IF (mixed_cdft%calculate_metric) THEN
         ALLOCATE (metric(nspins))
         metric = 0.0_dp
         CALL dbcsr_create(e_vectors, template=density_matrix(1, 1)%matrix)
         ! Take into account doubly occupied orbitals without LSD
         IF (nspins == 1) CALL dbcsr_scale(density_matrix(1, 1)%matrix, alpha_scalar=0.5_dp)
         DO ispin = 1, nspins
            ALLOCATE (evals(ncol_mo(ispin)))
            CALL cp_dbcsr_syevd(density_matrix(1, ispin)%matrix, e_vectors, evals, &
                                para_env=force_env%para_env, blacs_env=mixed_cdft%blacs_env)
            CALL dbcsr_release_p(density_matrix(1, ispin)%matrix)
            DO j = 1, ncol_mo(ispin)
               metric(ispin) = metric(ispin)+(evals(j)**2-evals(j)**4)
            END DO
            DEALLOCATE (evals)
         END DO
         CALL dbcsr_release(e_vectors)
         DEALLOCATE (density_matrix)
         metric = metric/4.0_dp
      END IF
      ! Compute coupling also with the wavefunction overlap method, see Migliore2009
      ! Requires the unconstrained KS ground state wavefunction as input
      IF (mixed_cdft%wfn_overlap_method) THEN
         ! Create mo_set for input wfn
         ALLOCATE (mo_set(nspins))
         IF (nspins == 2) THEN
            maxocc = 1.0_dp
         ELSE
            maxocc = 2.0_dp
         END IF
         DO ispin = 1, nspins
            nao = nrow_mo(ispin)
            nmo = ncol_mo(ispin)
            NULLIFY (mo_set(ispin)%mo_set)
            ! Only OT with fully occupied orbitals is implicitly supported (TODO: add check in _init)
            CALL allocate_mo_set(mo_set(ispin)%mo_set, nao=nao, nmo=nmo, nelectron=INT(maxocc*nmo), &
                                 n_el_f=REAL(maxocc*nmo, dp), maxocc=maxocc, &
                                 flexible_electron_count=0.0_dp)
            CALL set_mo_set(mo_set(ispin)%mo_set, uniform_occupation=.TRUE., homo=nmo)
            CALL cp_fm_create(matrix=mo_set(ispin)%mo_set%mo_coeff, &
                              matrix_struct=mixed_mo_coeff(1, ispin)%matrix%matrix_struct, &
                              name="GS_MO_COEFF"//TRIM(ADJUSTL(cp_to_string(ispin)))//"MATRIX")
            ALLOCATE (mo_set(ispin)%mo_set%eigenvalues(nmo))
            ALLOCATE (mo_set(ispin)%mo_set%occupation_numbers(nmo))
         END DO
         ! Read wfn file (note we assume that the basis set is the same)
         IF (force_env%mixed_env%do_mixed_qmmm_cdft) &
            ! This really shouldnt be a problem?
            CALL cp_abort(__LOCATION__, &
                          "QMMM + wavefunction overlap method not supported.")
         CALL force_env_get(force_env=force_env, subsys=subsys_mix)
         linear_section => section_vals_get_subs_vals(force_env_section, "MIXED%LINEAR")
         CALL cp_subsys_get(subsys_mix, atomic_kind_set=atomic_kind_set, particle_set=particle_set)
         CPASSERT(ASSOCIATED(mixed_cdft%qs_kind_set))
         IF (force_env%para_env%ionode) &
            CALL wfn_restart_file_name(file_name, exist, linear_section, logger)
         CALL mp_bcast(exist, force_env%para_env%source, force_env%para_env%group)
         CALL mp_bcast(file_name, force_env%para_env%source, force_env%para_env%group)
         IF (.NOT. exist) &
            CALL cp_abort(__LOCATION__, &
                          "User requested to restart the wavefunction from the file named: "// &
                          TRIM(file_name)//". This file does not exist. Please check the existence of"// &
                          " the file or change properly the value of the keyword WFN_RESTART_FILE_NAME in"// &
                          " section FORCE_EVAL\MIXED\LINEAR.")
         CALL read_mo_set(mo_array=mo_set, atomic_kind_set=atomic_kind_set, &
                          qs_kind_set=mixed_cdft%qs_kind_set, particle_set=particle_set, &
                          para_env=force_env%para_env, id_nr=0, multiplicity=mixed_cdft%multiplicity, &
                          dft_section=linear_section, natom_mismatch=natom_mismatch, &
                          cdft=.TRUE.)
         IF (natom_mismatch) &
            CALL cp_abort(__LOCATION__, &
                          "Restart wfn file has a wrong number of atoms")
         ! Orthonormalize wfn
         DO ispin = 1, nspins
            IF (mixed_cdft%has_unit_metric) THEN
               CALL make_basis_simple(mo_set(ispin)%mo_set%mo_coeff, ncol_mo(ispin))
            ELSE
               CALL make_basis_sm(mo_set(ispin)%mo_set%mo_coeff, ncol_mo(ispin), mixed_matrix_s)
            END IF
         END DO
         ! Calculate MO overlaps between reference state (R) and CDFT states A/D
         ALLOCATE (overlaps(2, nspins))
         overlaps = 0.0_dp
         DO ispin = 1, nspins
            nao = nrow_mo(ispin)
            nmo = ncol_mo(ispin)
            CALL cp_fm_struct_create(mo_mo_fmstruct, nrow_global=nmo, ncol_global=nmo, &
                                     context=mixed_cdft%blacs_env, para_env=force_env%para_env)
            CALL cp_fm_create(matrix=mo_overlap, matrix_struct=mo_mo_fmstruct, &
                              name="MO_OVERLAP_MATRIX")
            CALL cp_fm_create(matrix=inverse_mat, matrix_struct=mo_mo_fmstruct, &
                              name="INVERSE_MO_OVERLAP_MATRIX")
            CALL cp_fm_struct_release(mo_mo_fmstruct)
            CALL cp_fm_create(matrix=mo_tmp, &
                              matrix_struct=mixed_mo_coeff(1, ispin)%matrix%matrix_struct, &
                              name="OVERLAP_MO_COEFF")
            ! S*C_r
            CALL cp_dbcsr_sm_fm_multiply(mixed_matrix_s, mo_set(ispin)%mo_set%mo_coeff, &
                                         mo_tmp, nmo, 1.0_dp, 0.0_dp)
            DO j = 1, 2
               ! C_j^T * (S*C_r)
               CALL cp_fm_set_all(mo_overlap, alpha=0.0_dp)
               CALL cp_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, &
                            mixed_mo_coeff(j, ispin)%matrix, &
                            mo_tmp, 0.0_dp, mo_overlap)
               CALL cp_fm_invert(mo_overlap, inverse_mat, overlaps(j, ispin), eps_svd=mixed_cdft%eps_svd)
            END DO
            CALL cp_fm_release(mo_overlap)
            CALL cp_fm_release(inverse_mat)
            CALL cp_fm_release(mo_tmp)
            CALL deallocate_mo_set(mo_set(ispin)%mo_set)
         END DO
         DEALLOCATE (mo_set)
         IF (nspins == 2) THEN
            overlaps(1, 1) = ABS(overlaps(1, 1)*overlaps(1, 2)) ! A in eq. 12c
            overlaps(2, 1) = ABS(overlaps(2, 1)*overlaps(2, 2)) ! B in eq. 12c
         ELSE
            overlaps(1, 1) = overlaps(1, 1)**2
            overlaps(2, 1) = overlaps(2, 1)**2
         END IF
         ! Calculate coupling using eq. 12c
         coupling_wfn = ABS((overlaps(1, 1)*overlaps(2, 1)/(overlaps(1, 1)**2-overlaps(2, 1)**2))* &
                            (energy_b-energy_a)/(1.0_dp-Sda**2)* &
                            (1.0_dp-(overlaps(1, 1)**2+overlaps(2, 1)**2)/(2.0_dp*overlaps(1, 1)*overlaps(2, 1))*Sda))
         DEALLOCATE (overlaps)
      END IF
      ! Release remaining work
      DEALLOCATE (nrow_mo, ncol_mo)
      DO iforce_eval = 1, nforce_eval
         DO ispin = 1, nspins
            CALL cp_fm_release(mixed_mo_coeff(iforce_eval, ispin)%matrix)
         END DO
      END DO
      CALL dbcsr_release(mixed_wmat_a)
      CALL dbcsr_release(mixed_wmat_b)
      CALL dbcsr_release_p(mixed_matrix_s)
      DEALLOCATE (mixed_mo_coeff, rest_MO)
      ! construct S and W
      S_mat(1, 1) = 1.0_dp
      S_mat(2, 2) = 1.0_dp
      S_mat(2, 1) = Sda
      S_mat(1, 2) = Sda
      W_mat(1, 1) = Wbb
      W_mat(2, 2) = Waa
      W_mat(2, 1) = Wda
      W_mat(1, 2) = Wda
      ! Invert S via eigendecomposition and compute S^-(1/2)
      CALL diamat_all(S_mat, eigenv, .TRUE.)
      IF (eigenv(1) .LT. 1.0e-14_dp) eigenv(1) = 1.0e-14_dp ! Safeguard against division with 0 and negative numbers
      IF (eigenv(2) .LT. 1.0e-14_dp) eigenv(2) = 1.0e-14_dp
      U = 0.0_dp
      U(1, 1) = 1.0_dp/SQRT(eigenv(1))
      U(2, 2) = 1.0_dp/SQRT(eigenv(2))
      tmp_mat = MATMUL(U, TRANSPOSE(S_mat))
      U = MATMUL(S_mat, tmp_mat) ! S^(-1/2)
      ! Orthogonalize states i.e. solve generalized eigenvalue eq WV = SVL
      ! Convert to standard eigenvalue problem via symmetric orthogonalisation
      tmp_mat = MATMUL(W_mat, U) !
      W_mat = MATMUL(U, tmp_mat) ! W' = S^(-1/2) * W * S^(-1/2)
      CALL diamat_all(W_mat, eigenv, .TRUE.) ! Solve W'V' = AV'
      tmp_mat = MATMUL(U, W_mat) ! Reverse transformation V = S^(-1/2) V
      ! Construct final, orthogonal diabatic Hamiltonian matrix H
      W_mat(1, 1) = energy_b
      W_mat(2, 2) = energy_a
      c(1) = (energy_b+strength_b*Wbb)*Sda-strength_b*Wda ! H_AB = F_B*S_AB + V_B * W_AB
      c(2) = (energy_a+strength_a*Waa)*Sda-strength_a*Wda ! H_BA = F_A*S_BA - V_A * W_BA
      W_mat(1, 2) = (c(1)+c(2))*0.5_dp ! H''(1,2) = 0.5*(H_AB+H_AB) = H''(2,1)
      W_mat(2, 1) = W_mat(1, 2)
      S_mat = MATMUL(W_mat, (tmp_mat)) ! H'' * V
      W_mat = MATMUL(TRANSPOSE(tmp_mat), S_mat) ! H = V^T * H'' * V
      IF (iounit > 0) THEN
         WRITE (iounit, '(/,T3,A,T66)') &
            '------------------------- CDFT coupling information -------------------------'
         WRITE (iounit, '(T3,A,T66,(3X,F12.2))') &
            'Information at step (fs)          :', mixed_cdft%sim_step*mixed_cdft%sim_dt
         WRITE (iounit, '(T3,A,T60,(3X,F18.12))') 'Strength of constraint A          :', strength_a
         WRITE (iounit, '(T3,A,T60,(3X,F18.12))') 'Strength of constraint B          :', strength_b
         WRITE (iounit, '(T3,A,T60,(3X,F18.12))') 'Final value of constraint A       :', Waa
         WRITE (iounit, '(T3,A,T60,(3X,F18.12))') 'Final value of constraint B       :', Wbb
         WRITE (iounit, '(T3,A,T60,(3X,F18.12))') 'Overlap between states A and B    :', Sda
         WRITE (iounit, '(T3,A,T60,(3X,F18.12))') 'Charge transfer energy (B-A) (Ha) :', (energy_b-energy_a)
         WRITE (iounit, *)
         IF (ABS(W_mat(1, 2))*1.0E3_dp .GE. 0.1_dp) THEN
            WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
               'Diabatic electronic coupling (mHartree):', ABS(W_mat(1, 2)*1.0E3_dp)
         ELSE
            WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
               'Diabatic electronic coupling (microHartree):', ABS(W_mat(1, 2)*1.0E6_dp)
         END IF
         IF (mixed_cdft%calculate_metric) THEN
            WRITE (iounit, *)
            IF (nspins == 1) THEN
               WRITE (iounit, '(T3,A,T66,(3X,F12.6))') 'Coupling reliability metric (0 is ideal):', metric
            ELSE
               WRITE (iounit, '(T3,A,T54,(3X,2F12.6))') &
                  'Coupling reliability metric (0 is ideal):', metric(1), metric(2)
            END IF
         END IF
         IF (mixed_cdft%wfn_overlap_method) THEN
            WRITE (iounit, *)
            IF (coupling_wfn*1.0E3_dp .GE. 0.1_dp) THEN
               WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
                  'Diabatic coupling (wavefunction overlap method, mHartree):', coupling_wfn*1.0E3_dp
            ELSE
               WRITE (iounit, '(T3,A,T60,(3X,F18.12))') &
                  'Diabatic coupling (wavefunction overlap method, microHartree):', coupling_wfn*1.0E6_dp
            END IF
         END IF
         WRITE (iounit, '(T3,A)') &
            '------------------------------------------------------------------------------'
      END IF
      DEALLOCATE (a, b, S_det)
      IF (mixed_cdft%calculate_metric) DEALLOCATE (metric)

      CALL cp_print_key_finished_output(iounit, logger, force_env_section, &
                                        "MIXED%PRINT%PROGRAM_RUN_INFO")
      CALL timestop(handle)

   END SUBROUTINE mixed_cdft_calculate_coupling
! **************************************************************************************************
!> \brief Becke constraint adapted to mixed calculations, details in et_coupling.F
!> \param force_env the force_env that holds the CDFT states
!> \param calculate_forces determines if forces should be calculted
!> \par History
!>       02.2016  created [Nico Holmberg]
!>       03.2016  added dynamic load balancing (dlb)
!>                changed pw_p_type data types to rank-3 reals to accommodate dlb
!>                and to reduce overall memory footprint
!>                split to subroutines [Nico Holmberg]
!>       04.2016  introduced mixed grid mapping [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_becke_constraint(force_env, calculate_forces)
      TYPE(force_env_type), POINTER                      :: force_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_becke_constraint', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: catom
      LOGICAL                                            :: in_memory, store_vectors
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: is_constraint
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: coefficients
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: position_vecs, R12
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: pair_dist_vecs
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(mixed_environment_type), POINTER              :: mixed_env

      NULLIFY (mixed_env, mixed_cdft)
      store_vectors = .TRUE.
      logger => cp_get_default_logger()
      CALL timeset(routineN, handle)
      mixed_env => force_env%mixed_env
      CALL get_mixed_env(mixed_env, cdft_control=mixed_cdft)
      CALL mixed_becke_constraint_init(force_env, mixed_cdft, calculate_forces, &
                                       is_constraint, in_memory, store_vectors, &
                                       R12, position_vecs, pair_dist_vecs, &
                                       coefficients, catom)
      CALL mixed_becke_constraint_low(force_env, mixed_cdft, in_memory, &
                                      is_constraint, store_vectors, R12, &
                                      position_vecs, pair_dist_vecs, &
                                      coefficients, catom)
      CALL timestop(handle)

   END SUBROUTINE mixed_becke_constraint
! **************************************************************************************************
!> \brief Initialize the mixed Becke constraint calculation
!> \param force_env the force_env that holds the CDFT states
!> \param mixed_cdft container for structures needed to build the mixed CDFT constraint
!> \param calculate_forces determines if forces should be calculted
!> \param is_constraint a list used to determine which atoms in the system define the constraint
!> \param in_memory decides whether to build the weight function gradients in parallel before solving
!>                  the CDFT states or later during the SCF procedure of the individual states
!> \param store_vectors should temporary arrays be stored in memory to accelerate the calculation
!> \param R12 temporary array holding the pairwise atomic distances
!> \param position_vecs temporary array holding the pbc corrected atomic position vectors
!> \param pair_dist_vecs temporary array holding the pairwise displament vectors
!> \param coefficients array that determines how atoms should be summed to form the constraint
!> \param catom temporary array to map the global index of constraint atoms to their position
!>              in a list that holds only constraint atoms
!> \par History
!>       03.2016  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_becke_constraint_init(force_env, mixed_cdft, calculate_forces, &
                                          is_constraint, in_memory, store_vectors, &
                                          R12, position_vecs, pair_dist_vecs, coefficients, &
                                          catom)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      LOGICAL, INTENT(IN)                                :: calculate_forces
      LOGICAL, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: is_constraint
      LOGICAL, INTENT(OUT)                               :: in_memory
      LOGICAL, INTENT(IN)                                :: store_vectors
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(out)                                     :: R12, position_vecs
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(out)                                     :: pair_dist_vecs
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: coefficients
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(out)    :: catom

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_becke_constraint_init', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=2)                                   :: element_symbol
      INTEGER :: atom_a, bounds(2), dir, handle, i, iatom, iex, iforce_eval, ikind, iounit, &
         ithread, j, jatom, katom, lb_index, my_work, my_work_size, natom, nforce_eval, nkind, &
         np(3), npme, nthread, numexp, offset_dlb, tmp_index(2), ub_index, unit_nr
      INTEGER, DIMENSION(2, 3)                           :: bo, bo_conf
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores, stride
      LOGICAL                                            :: build
      REAL(kind=dp)                                      :: alpha, chi, coef, dynamic_radius, ircov, &
                                                            jrcov, lb, ra(3), tmp_const, ub, uij
      REAL(kind=dp), DIMENSION(3)                        :: cell_v, dist_vec, dr, r, r1, shift
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys_mix
      TYPE(force_env_type), POINTER                      :: force_env_qs
      TYPE(hirshfeld_type), POINTER                      :: cavity_env
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_type), POINTER                 :: rs_cavity
      TYPE(section_vals_type), POINTER                   :: force_env_section, print_section

      NULLIFY (pab, cell, force_env_qs, particle_set, force_env_section, print_section, &
               qs_kind_set, particles, subsys_mix, rs_cavity, cavity_env, auxbas_pw_pool, &
               atomic_kind_set)
      logger => cp_get_default_logger()
      nforce_eval = SIZE(force_env%sub_force_env)
      CALL timeset(routineN, handle)
      CALL force_env_get(force_env=force_env, force_env_section=force_env_section)
      IF (.NOT. force_env%mixed_env%do_mixed_qmmm_cdft) THEN
         CALL force_env_get(force_env=force_env, &
                            subsys=subsys_mix, &
                            cell=cell)
         CALL cp_subsys_get(subsys=subsys_mix, &
                            particles=particles, &
                            particle_set=particle_set)
      ELSE
         DO iforce_eval = 1, nforce_eval
            IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
            force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
         END DO
         CALL get_qs_env(force_env_qs%qmmm_env%qs_env, &
                         cp_subsys=subsys_mix, &
                         cell=cell)
         CALL cp_subsys_get(subsys=subsys_mix, &
                            particles=particles, &
                            particle_set=particle_set)
      END IF
      natom = SIZE(particles%els)
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%PRINT%PROGRAM_RUN_INFO")
      IF (.NOT. ASSOCIATED(mixed_cdft%becke_control%cutoffs)) THEN
         CALL cp_subsys_get(subsys_mix, atomic_kind_set=atomic_kind_set)
         ALLOCATE (mixed_cdft%becke_control%cutoffs(natom))
         SELECT CASE (mixed_cdft%becke_control%cutoff_type)
         CASE (becke_cutoff_global)
            mixed_cdft%becke_control%cutoffs(:) = mixed_cdft%becke_control%rglobal
         CASE (becke_cutoff_element)
            IF (.NOT. SIZE(atomic_kind_set) == SIZE(mixed_cdft%becke_control%cutoffs_tmp)) &
               CALL cp_abort(__LOCATION__, &
                             "Size of keyword BECKE_RESTRAINT\ELEMENT_CUTOFFS does"// &
                             "not match number of atomic kinds in the input coordinate file.")
            DO ikind = 1, SIZE(atomic_kind_set)
               CALL get_atomic_kind(atomic_kind_set(ikind), natom=katom, atom_list=atom_list)
               DO iatom = 1, katom
                  atom_a = atom_list(iatom)
                  mixed_cdft%becke_control%cutoffs(atom_a) = mixed_cdft%becke_control%cutoffs_tmp(ikind)
               END DO
            END DO
            DEALLOCATE (mixed_cdft%becke_control%cutoffs_tmp)
         END SELECT
      END IF
      build = .FALSE.
      IF (mixed_cdft%becke_control%adjust .AND. &
          .NOT. ASSOCIATED(mixed_cdft%becke_control%aij)) THEN
         ALLOCATE (mixed_cdft%becke_control%aij(natom, natom))
         build = .TRUE.
      END IF
      ALLOCATE (catom(mixed_cdft%becke_control%natoms))
      IF (mixed_cdft%becke_control%dynamic_confine .OR. &
          mixed_cdft%becke_control%save_pot .OR. &
          mixed_cdft%becke_control%cavity_confine .OR. &
          mixed_cdft%becke_control%should_skip .OR. &
          mixed_cdft%first_iteration) THEN
         ALLOCATE (is_constraint(natom))
         is_constraint = .FALSE.
      END IF
      in_memory = calculate_forces .AND. mixed_cdft%becke_control%in_memory
      IF (in_memory .NEQV. calculate_forces) &
         CALL cp_abort(__LOCATION__, &
                       "The flag BECKE_RESTRAINT\IN_MEMORY must be activated "// &
                       "for the calculation of mixed CDFT forces")
      IF (in_memory .OR. mixed_cdft%first_iteration) ALLOCATE (coefficients(natom))
      DO i = 1, mixed_cdft%becke_control%natoms
         catom(i) = mixed_cdft%becke_control%atoms(i)
         IF (mixed_cdft%becke_control%dynamic_confine .OR. &
             mixed_cdft%becke_control%save_pot .OR. &
             mixed_cdft%becke_control%cavity_confine .OR. &
             mixed_cdft%becke_control%should_skip .OR. &
             mixed_cdft%first_iteration) &
            is_constraint(catom(i)) = .TRUE.
         IF (in_memory .OR. mixed_cdft%first_iteration) &
            coefficients(catom(i)) = mixed_cdft%becke_control%coeff(i)
      ENDDO
      CALL pw_env_get(pw_env=mixed_cdft%pw_env, auxbas_pw_pool=auxbas_pw_pool)
      bo = auxbas_pw_pool%pw_grid%bounds_local
      np = auxbas_pw_pool%pw_grid%npts
      dr = auxbas_pw_pool%pw_grid%dr
      shift = -REAL(MODULO(np, 2), dp)*dr/2.0_dp
      IF (store_vectors) THEN
         IF (in_memory) ALLOCATE (pair_dist_vecs(3, natom, natom))
         ALLOCATE (position_vecs(3, natom))
      END IF
      DO i = 1, 3
         cell_v(i) = cell%hmat(i, i)
      END DO
      IF (mixed_cdft%becke_control%confine) THEN
         dir = mixed_cdft%becke_control%confine_dir
         lb = mixed_cdft%becke_control%confine_bounds(1)
         ub = mixed_cdft%becke_control%confine_bounds(2)
         lb = lb-cell%hmat(dir, dir)/2._dp
         ub = ub-cell%hmat(dir, dir)/2._dp
      ELSE IF (mixed_cdft%becke_control%dynamic_confine) THEN
         dir = mixed_cdft%becke_control%confine_dir
         lb = HUGE(0.0_dp)
         ub = -HUGE(0.0_dp)
         dynamic_radius = mixed_cdft%becke_control%dynamic_radius
      END IF
      ALLOCATE (R12(natom, natom))
      DO iatom = 1, natom-1
         DO jatom = iatom+1, natom
            r = particle_set(iatom)%r
            r1 = particle_set(jatom)%r
            DO i = 1, 3
               r(i) = MODULO(r(i), cell%hmat(i, i))-cell%hmat(i, i)/2._dp
               r1(i) = MODULO(r1(i), cell%hmat(i, i))-cell%hmat(i, i)/2._dp
            END DO
            dist_vec = (r-r1)-ANINT((r-r1)/cell_v)*cell_v
            IF (store_vectors) THEN
               position_vecs(:, iatom) = r(:)
               IF (iatom == 1 .AND. jatom == natom) position_vecs(:, jatom) = r1(:)
               IF (in_memory) THEN
                  pair_dist_vecs(:, iatom, jatom) = dist_vec(:)
                  pair_dist_vecs(:, jatom, iatom) = -dist_vec(:)
               END IF
            END IF
            R12(iatom, jatom) = SQRT(DOT_PRODUCT(dist_vec, dist_vec))
            R12(jatom, iatom) = R12(iatom, jatom)
            IF (mixed_cdft%becke_control%dynamic_confine .AND. &
                is_constraint(iatom) .AND. is_constraint(jatom)) THEN
               IF ((r(dir)+dynamic_radius) .GT. ub) ub = r(dir)+dynamic_radius
               IF ((r1(dir)+dynamic_radius) .GT. ub) ub = r1(dir)+dynamic_radius
               IF ((r(dir)-dynamic_radius) .LT. lb) lb = r(dir)-dynamic_radius
               IF ((r1(dir)-dynamic_radius) .LT. lb) lb = r1(dir)-dynamic_radius
            END IF
            IF (build) THEN
               CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                    kind_number=ikind)
               ircov = mixed_cdft%becke_control%radii(ikind)
               CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, &
                                    kind_number=ikind)
               jrcov = mixed_cdft%becke_control%radii(ikind)
               IF (ircov .NE. jrcov) THEN
                  chi = ircov/jrcov
                  uij = (chi-1.0_dp)/(chi+1.0_dp)
                  mixed_cdft%becke_control%aij(iatom, jatom) = uij/(uij**2-1.0_dp)
                  IF (mixed_cdft%becke_control%aij(iatom, jatom) &
                      .GT. 0.5_dp) THEN
                     mixed_cdft%becke_control%aij(iatom, jatom) = 0.5_dp
                  ELSE IF (mixed_cdft%becke_control%aij(iatom, jatom) &
                           .LT. -0.5_dp) THEN
                     mixed_cdft%becke_control%aij(iatom, jatom) = -0.5_dp
                  END IF
               ELSE
                  mixed_cdft%becke_control%aij(iatom, jatom) = 0.0_dp
               END IF
               mixed_cdft%becke_control%aij(jatom, iatom) = &
                  -mixed_cdft%becke_control%aij(iatom, jatom)
            END IF
         END DO
      END DO
      ! Dump some additional information about the calculation
      IF (mixed_cdft%first_iteration) THEN
         print_section => section_vals_get_subs_vals(force_env_section, "MIXED%PRINT%PROGRAM_RUN_INFO")
         iounit = cp_print_key_unit_nr(logger, print_section, '', extension='.mixedLog')
         IF (iounit > 0) THEN
            WRITE (iounit, '(/,T3,A,T66)') &
               '----------------------------- Becke atomic parameters -----------------------------'
            IF (mixed_cdft%becke_control%adjust) THEN
               WRITE (iounit, '(/,T3,A,A)') &
                  'Atom  Element   Coefficient', '     Cutoff (angstrom)    CDFT Radius (angstrom)'
               DO iatom = 1, natom
                  CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                       element_symbol=element_symbol, &
                                       kind_number=ikind)
                  ircov = cp_unit_from_cp2k(mixed_cdft%becke_control%radii(ikind), "angstrom")
                  IF (is_constraint(iatom)) THEN
                     coef = coefficients(iatom)
                  ELSE
                     coef = 0.0_dp
                  END IF
                  WRITE (iounit, "(i7,T15,A2,T23,F8.3,T39,F8.3,T61,F8.3,T75,F8.3)") &
                     iatom, element_symbol, coef, &
                     cp_unit_from_cp2k(mixed_cdft%becke_control%cutoffs(iatom), "angstrom"), &
                     ircov
               END DO
            ELSE
               WRITE (iounit, '(/,T3,A,A)') &
                  'Atom  Element   Coefficient', '     Cutoff (angstrom)'
               DO iatom = 1, natom
                  CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                       element_symbol=element_symbol)
                  IF (is_constraint(iatom)) THEN
                     coef = coefficients(iatom)
                  ELSE
                     coef = 0.0_dp
                  END IF
                  WRITE (iounit, "(i7,T15,A2,T23,F8.3,T39,F8.3,T61,F8.3,T75)") &
                     iatom, element_symbol, coef, &
                     cp_unit_from_cp2k(mixed_cdft%becke_control%cutoffs(iatom), "angstrom")
               END DO
            END IF
            WRITE (iounit, '(/,T3,A)') &
               '-----------------------------------------------------------------------------------'
         END IF
         CALL cp_print_key_finished_output(iounit, logger, force_env_section, &
                                           "MIXED%PRINT%PROGRAM_RUN_INFO")
         mixed_cdft%first_iteration = .FALSE.
      END IF

      IF (mixed_cdft%becke_control%dynamic_confine .OR. &
          mixed_cdft%becke_control%confine) THEN
         tmp_const = (ub-shift(dir))/dr(dir)
         IF (tmp_const .GE. 0.0_dp) THEN
            ub_index = CEILING(tmp_const)
         ELSE
            ub_index = FLOOR(tmp_const)
         END IF
         tmp_const = (lb-shift(dir))/dr(dir)
         IF (tmp_const .GE. 0.0_dp) THEN
            lb_index = CEILING(tmp_const)
         ELSE
            lb_index = FLOOR(tmp_const)
         END IF
         mixed_cdft%becke_control%confine_bounds_int(1) = lb_index
         mixed_cdft%becke_control%confine_bounds_int(2) = ub_index
      END IF
      IF (mixed_cdft%becke_control%cavity_confine) THEN
         CPASSERT(ASSOCIATED(mixed_cdft%qs_kind_set))
         cavity_env => mixed_cdft%becke_control%cavity_env
         qs_kind_set => mixed_cdft%qs_kind_set
         CALL cp_subsys_get(subsys_mix, atomic_kind_set=atomic_kind_set)
         nkind = SIZE(qs_kind_set)
         IF (.NOT. ASSOCIATED(cavity_env%kind_shape_fn)) &
            CALL create_shape_function(cavity_env, qs_kind_set, atomic_kind_set, &
                                       radius=mixed_cdft%becke_control%rcavity, &
                                       radii_list=mixed_cdft%becke_control%radii)
         NULLIFY (rs_cavity)
         CALL pw_env_get(pw_env=mixed_cdft%pw_env, auxbas_rs_grid=rs_cavity, &
                         auxbas_pw_pool=auxbas_pw_pool)
         ! be careful in parallel nsmax is choosen with multigrid in mind!
         CALL rs_grid_retain(rs_cavity)
         CALL rs_grid_zero(rs_cavity)
         ALLOCATE (pab(1, 1))
         nthread = 1
         ithread = 0
         DO ikind = 1, SIZE(atomic_kind_set)
            numexp = cavity_env%kind_shape_fn(ikind)%numexp
            IF (numexp <= 0) CYCLE
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=katom, atom_list=atom_list)
            ALLOCATE (cores(katom))
            DO iex = 1, numexp
               alpha = cavity_env%kind_shape_fn(ikind)%zet(iex)
               coef = cavity_env%kind_shape_fn(ikind)%coef(iex)
               npme = 0
               cores = 0
               DO iatom = 1, katom
                  IF (rs_cavity%desc%parallel .AND. .NOT. rs_cavity%desc%distributed) THEN
                     ! replicated realspace grid, split the atoms up between procs
                     IF (MODULO(iatom, rs_cavity%desc%group_size) == rs_cavity%desc%my_pos) THEN
                        npme = npme+1
                        cores(npme) = iatom
                     ENDIF
                  ELSE
                     npme = npme+1
                     cores(npme) = iatom
                  ENDIF
               END DO
               DO j = 1, npme
                  iatom = cores(j)
                  atom_a = atom_list(iatom)
                  pab(1, 1) = coef
                  IF (store_vectors) THEN
                     ra(:) = position_vecs(:, atom_a)+cell_v(:)/2._dp
                  ELSE
                     ra(:) = pbc(particle_set(atom_a)%r, cell)
                  END IF
                  IF (is_constraint(atom_a)) &
                     CALL collocate_pgf_product_rspace(0, alpha, 0, 0, 0.0_dp, 0, ra, &
                                                       (/0.0_dp, 0.0_dp, 0.0_dp/), 0.0_dp, 1.0_dp, pab, 0, 0, &
                                                       rs_cavity, cell, mixed_cdft%pw_env%cube_info(1), &
                                                       mixed_cdft%eps_rho_rspace, ga_gb_function=FUNC_AB, &
                                                       ithread=ithread, use_subpatch=.TRUE., &
                                                       subpatch_pattern=0_int_8, lmax_global=0)
               END DO
            END DO
            DEALLOCATE (cores)
         END DO
         DEALLOCATE (pab)
         CALL pw_pool_create_pw(auxbas_pw_pool, mixed_cdft%becke_control%cavity%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL rs_pw_transfer(rs_cavity, mixed_cdft%becke_control%cavity%pw, rs2pw)
         CALL rs_grid_release(rs_cavity)
         IF (.NOT. mixed_cdft%becke_control%dynamic_confine) THEN
            CALL hfun_zero(mixed_cdft%becke_control%cavity%pw%cr3d, &
                           mixed_cdft%becke_control%eps_cavity, &
                           just_zero=.FALSE., bounds=bounds, work=my_work)
            IF (mixed_cdft%is_pencil .OR. mixed_cdft%is_special) THEN
               my_work_size = (bounds(2)-bounds(1)+1)*(bo(2, 2)-bo(1, 2)+1)
            ELSE
               my_work_size = (bounds(2)-bounds(1)+1)*(bo(2, 1)-bo(1, 1)+1)
            END IF
            mixed_cdft%becke_control%confine_bounds_int(1) = bounds(1)
            mixed_cdft%becke_control%confine_bounds_int(2) = bounds(2)
         END IF
         IF (mixed_cdft%becke_control%print_cavity) THEN
            CALL hfun_zero(mixed_cdft%becke_control%cavity%pw%cr3d, &
                           mixed_cdft%becke_control%eps_cavity, just_zero=.TRUE.)
            NULLIFY (stride)
            ALLOCATE (stride(3))
            stride = (/2, 2, 2/)
            unit_nr = cp_print_key_unit_nr(logger, print_section, "", &
                                           middle_name="BECKE_CAVITY", &
                                           extension=".cube", file_position="REWIND", &
                                           log_filename=.FALSE.)
            IF (force_env%para_env%mepos == force_env%para_env%source .AND. unit_nr .LT. 1) &
               CALL cp_abort(__LOCATION__, &
                             "Please turn on PROGRAM_RUN_INFO to print cavity")
            CALL cp_pw_to_cube(mixed_cdft%becke_control%cavity%pw, &
                               unit_nr, "CAVITY", particles=particles, &
                               stride=stride)
            CALL cp_print_key_finished_output(unit_nr, logger, print_section, '')
            DEALLOCATE (stride)
         END IF
      END IF
      bo_conf = bo
      IF (mixed_cdft%becke_control%dynamic_confine .OR. &
          mixed_cdft%becke_control%confine) THEN
         IF (ub_index .GT. bo(2, dir)) THEN
            tmp_index(2) = bo(2, dir)-1
         ELSE
            tmp_index(2) = ub_index-1
         END IF
         IF (lb_index .LT. bo(1, dir)) THEN
            tmp_index(1) = bo(1, dir)+1
         ELSE
            tmp_index(1) = lb_index+1
         END IF
         bo_conf(2, dir) = tmp_index(2)
         bo_conf(1, dir) = tmp_index(1)
      ELSE IF (mixed_cdft%becke_control%cavity_confine) THEN
         IF (bounds(2) .LT. bo(2, 3)) bo_conf(2, 3) = bounds(2)-1
         IF (bounds(1) .GT. bo(1, 3)) bo_conf(1, 3) = bounds(1)+1
      END IF

      ! Load balance
      IF (mixed_cdft%dlb) &
         CALL mixed_becke_constraint_dlb(force_env, mixed_cdft, my_work, &
                                         my_work_size, natom, bo, bo_conf)
      ! The bounds have been finalized => time to allocate storage for working matrices
      offset_dlb = 0
      IF (mixed_cdft%dlb) THEN
         IF (mixed_cdft%dlb_control%send_work .AND. .NOT. mixed_cdft%is_special) &
            offset_dlb = SUM(mixed_cdft%dlb_control%target_list(2, :))
      END IF
      IF (mixed_cdft%becke_control%cavity_confine) THEN
         ! Get rid of the zero part of the confinement cavity (cr3d -> real(:,:,:))
         IF (mixed_cdft%is_special) THEN
            ALLOCATE (mixed_cdft%sendbuff(SIZE(mixed_cdft%dest_list)))
            DO i = 1, SIZE(mixed_cdft%dest_list)
               ALLOCATE (mixed_cdft%sendbuff(i)%cavity(mixed_cdft%dest_list_bo(1, i):mixed_cdft%dest_list_bo(2, i), &
                                                       bo(1, 2):bo(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
               mixed_cdft%sendbuff(i)%cavity = mixed_cdft%becke_control%cavity%pw%cr3d(mixed_cdft%dest_list_bo(1, i): &
                                                                                       mixed_cdft%dest_list_bo(2, i), &
                                                                                       bo(1, 2):bo(2, 2), &
                                                                                       bo_conf(1, 3):bo_conf(2, 3))
            END DO
         ELSE IF (mixed_cdft%is_pencil) THEN
            ALLOCATE (mixed_cdft%cavity(bo(1, 1)+offset_dlb:bo(2, 1), bo(1, 2):bo(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
            mixed_cdft%cavity = mixed_cdft%becke_control%cavity%pw%cr3d(bo(1, 1)+offset_dlb:bo(2, 1), &
                                                                        bo(1, 2):bo(2, 2), &
                                                                        bo_conf(1, 3):bo_conf(2, 3))
         ELSE
            ALLOCATE (mixed_cdft%cavity(bo(1, 1):bo(2, 1), bo(1, 2)+offset_dlb:bo(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
            mixed_cdft%cavity = mixed_cdft%becke_control%cavity%pw%cr3d(bo(1, 1):bo(2, 1), &
                                                                        bo(1, 2)+offset_dlb:bo(2, 2), &
                                                                        bo_conf(1, 3):bo_conf(2, 3))
         END IF
         CALL pw_pool_give_back_pw(auxbas_pw_pool, mixed_cdft%becke_control%cavity%pw)
      END IF
      IF (mixed_cdft%is_special) THEN
         DO i = 1, SIZE(mixed_cdft%dest_list)
            ALLOCATE (mixed_cdft%sendbuff(i)%weight(mixed_cdft%dest_list_bo(1, i):mixed_cdft%dest_list_bo(2, i), &
                                                    bo(1, 2):bo(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
            mixed_cdft%sendbuff(i)%weight = 0.0_dp
         END DO
      ELSE IF (mixed_cdft%is_pencil) THEN
         ALLOCATE (mixed_cdft%weight(bo(1, 1)+offset_dlb:bo(2, 1), bo(1, 2):bo(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
         mixed_cdft%weight = 0.0_dp
      ELSE
         ALLOCATE (mixed_cdft%weight(bo(1, 1):bo(2, 1), bo(1, 2)+offset_dlb:bo(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
         mixed_cdft%weight = 0.0_dp
      END IF
      IF (in_memory) THEN
         IF (mixed_cdft%is_special) THEN
            DO i = 1, SIZE(mixed_cdft%dest_list)
               ALLOCATE (mixed_cdft%sendbuff(i)%gradients(3*natom, mixed_cdft%dest_list_bo(1, i): &
                                                          mixed_cdft%dest_list_bo(2, i), &
                                                          bo(1, 2):bo(2, 2), &
                                                          bo_conf(1, 3):bo_conf(2, 3)))
               mixed_cdft%sendbuff(i)%gradients = 0.0_dp
            END DO
         ELSE IF (mixed_cdft%is_pencil) THEN
            ALLOCATE (mixed_cdft%becke_control%gradients(3*natom, bo(1, 1)+offset_dlb:bo(2, 1), &
                                                         bo(1, 2):bo(2, 2), &
                                                         bo_conf(1, 3):bo_conf(2, 3)))
            mixed_cdft%becke_control%gradients = 0.0_dp
         ELSE
            ALLOCATE (mixed_cdft%becke_control%gradients(3*natom, bo(1, 1):bo(2, 1), &
                                                         bo(1, 2)+offset_dlb:bo(2, 2), &
                                                         bo_conf(1, 3):bo_conf(2, 3)))
            mixed_cdft%becke_control%gradients = 0.0_dp
         END IF
         IF (mixed_cdft%becke_control%should_skip .AND. .NOT. &
             mixed_cdft%becke_control%cavity_confine) THEN
            ALLOCATE (mixed_cdft%becke_control%skip_list(bo_conf(1, 1):bo_conf(2, 1), &
                                                         bo_conf(1, 2):bo_conf(2, 2), &
                                                         bo_conf(1, 3):bo_conf(2, 3)))
            mixed_cdft%becke_control%skip_list = .FALSE.
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE mixed_becke_constraint_init

! **************************************************************************************************
!> \brief Setup load balancing for mixed Becke calculation
!> \param force_env the force_env that holds the CDFT states
!> \param mixed_cdft container for structures needed to build the mixed CDFT constraint
!> \param my_work an estimate of the work per processor
!> \param my_work_size size of the smallest array slice per processor. overloaded processors will
!>                     redistribute works as integer multiples of this value.
!> \param natom the total number of atoms
!> \param bo bounds of the realspace grid that holds the electron density
!> \param bo_conf same as bo, but bounds along z-direction have been compacted with confinement
!> \par History
!>       03.2016  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_becke_constraint_dlb(force_env, mixed_cdft, my_work, &
                                         my_work_size, natom, bo, bo_conf)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      INTEGER, INTENT(IN)                                :: my_work, my_work_size, natom
      INTEGER, DIMENSION(2, 3)                           :: bo, bo_conf

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_becke_constraint_dlb', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: should_deallocate = 7000, &
                                                            uninitialized = -7000

      INTEGER :: actually_sent, exhausted_work, handle, i, ind, iounit, ispecial, j, max_targets, &
         more_work, my_pos, my_req(2), my_special_work, my_target, no_overloaded, no_underloaded, &
         nsend, nsend_limit, nsend_max, offset, offset_proc, offset_special, req(4), send_total, &
         tags(2)
      INTEGER, DIMENSION(:), POINTER :: buffsize, cumulative_work, expected_work, load_imbalance, &
         nrecv, nsend_proc, req_recv, req_total, sendbuffer, should_warn, tmp, work_index, &
         work_size
      INTEGER, DIMENSION(:, :), POINTER                  :: targets, tmp_bo
      LOGICAL                                            :: consistent
      LOGICAL, DIMENSION(:), POINTER                     :: mask_recv, mask_send, touched
      REAL(kind=dp)                                      :: average_work, load_scale, &
                                                            very_overloaded, work_factor
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: cavity
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: force_env_section, print_section

      TYPE buffers
         LOGICAL, POINTER, DIMENSION(:)        :: bv
         INTEGER, POINTER, DIMENSION(:)        :: iv
      END TYPE buffers
      TYPE(buffers), POINTER, DIMENSION(:)     :: recvbuffer, sbuff
      CHARACTER(len=2)                         :: dummy

      logger => cp_get_default_logger()
      CALL timeset(routineN, handle)
      mixed_cdft%dlb_control%recv_work = .FALSE.
      mixed_cdft%dlb_control%send_work = .FALSE.
      NULLIFY (expected_work, work_index, load_imbalance, work_size, &
               cumulative_work, sendbuffer, buffsize, req_recv, req_total, &
               tmp, nrecv, nsend_proc, targets, tmp_bo, touched, &
               mask_recv, mask_send, cavity, recvbuffer, sbuff, force_env_section, &
               print_section)
      CALL force_env_get(force_env=force_env, force_env_section=force_env_section)
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%PRINT%PROGRAM_RUN_INFO")
      iounit = cp_print_key_unit_nr(logger, print_section, '', extension='.mixedLog')
      ! These numerical values control data redistribution and are system sensitive
      ! Currently they are not refined during run time which may cause crashes
      ! However, using too many processors or a confinement cavity that is too large relative to the
      ! total system volume are more likely culprits.
      load_scale = mixed_cdft%dlb_control%load_scale
      very_overloaded = mixed_cdft%dlb_control%very_overloaded
      more_work = mixed_cdft%dlb_control%more_work
      max_targets = 40
      work_factor = 0.8_dp
      ! Reset targets/sources
      IF (mixed_cdft%is_special) THEN
         DEALLOCATE (mixed_cdft%dest_list, mixed_cdft%dest_list_bo, &
                     mixed_cdft%source_list, mixed_cdft%source_list_bo)
         ALLOCATE (mixed_cdft%dest_list(SIZE(mixed_cdft%dest_list_save)), &
                   mixed_cdft%dest_list_bo(SIZE(mixed_cdft%dest_bo_save, 1), SIZE(mixed_cdft%dest_bo_save, 2)), &
                   mixed_cdft%source_list(SIZE(mixed_cdft%source_list_save)), &
                   mixed_cdft%source_list_bo(SIZE(mixed_cdft%source_bo_save, 1), SIZE(mixed_cdft%source_bo_save, 2)))
         mixed_cdft%dest_list = mixed_cdft%dest_list_save
         mixed_cdft%source_list = mixed_cdft%source_list_save
         mixed_cdft%dest_list_bo = mixed_cdft%dest_bo_save
         mixed_cdft%source_list_bo = mixed_cdft%source_bo_save
      END IF
      ALLOCATE (mixed_cdft%dlb_control%expected_work(force_env%para_env%num_pe), &
                expected_work(force_env%para_env%num_pe), &
                work_size(force_env%para_env%num_pe))
      IF (debug_this_module) THEN
         ALLOCATE (should_warn(force_env%para_env%num_pe))
         should_warn = 0
      END IF
      expected_work = 0
      expected_work(force_env%para_env%mepos+1) = my_work
      work_size = 0
      work_size(force_env%para_env%mepos+1) = my_work_size
      IF (ASSOCIATED(mixed_cdft%dlb_control%prediction_error)) THEN
         IF (mixed_cdft%is_pencil .OR. mixed_cdft%is_special) THEN
            work_size(force_env%para_env%mepos+1) = work_size(force_env%para_env%mepos+1)- &
                                                    NINT(REAL(mixed_cdft%dlb_control% &
                                                              prediction_error(force_env%para_env%mepos+1), dp)/ &
                                                         REAL(bo(2, 1)-bo(1, 1)+1, dp))
         ELSE
            work_size(force_env%para_env%mepos+1) = work_size(force_env%para_env%mepos+1)- &
                                                    NINT(REAL(mixed_cdft%dlb_control% &
                                                              prediction_error(force_env%para_env%mepos+1), dp)/ &
                                                         REAL(bo(2, 2)-bo(1, 2)+1, dp))
         END IF
      END IF
      CALL mp_sum(expected_work, force_env%para_env%group)
      CALL mp_sum(work_size, force_env%para_env%group)
      ! We store the unsorted expected work to refine the estimate on subsequent calls to this routine
      mixed_cdft%dlb_control%expected_work = expected_work
      ! Take into account the prediction error of the last step
      IF (ASSOCIATED(mixed_cdft%dlb_control%prediction_error)) &
         expected_work = expected_work-mixed_cdft%dlb_control%prediction_error
      !
      average_work = REAL(SUM(expected_work), dp)/REAL(force_env%para_env%num_pe, dp)
      ALLOCATE (work_index(force_env%para_env%num_pe), &
                load_imbalance(force_env%para_env%num_pe), &
                targets(2, force_env%para_env%num_pe))
      load_imbalance = expected_work-NINT(average_work)
      no_overloaded = 0
      no_underloaded = 0
      targets = 0
      ! Convert the load imbalance to a multiple of the actual work size
      DO i = 1, force_env%para_env%num_pe
         IF (load_imbalance(i) .GT. 0) THEN
            no_overloaded = no_overloaded+1
            ! Allow heavily overloaded processors to dump more data since most likely they have a lot of 'real' work
            IF (expected_work(i) .GT. NINT(very_overloaded*average_work)) THEN
               load_imbalance(i) = (CEILING(REAL(load_imbalance(i), dp)/REAL(work_size(i), dp))+more_work)*work_size(i)
            ELSE
               load_imbalance(i) = CEILING(REAL(load_imbalance(i), dp)/REAL(work_size(i), dp))*work_size(i)
            END IF
         ELSE
            ! Allow the underloaded processors to take load_scale amount of additional work
            ! otherwise we may be unable to exhaust all overloaded processors
            load_imbalance(i) = NINT(load_imbalance(i)*load_scale)
            no_underloaded = no_underloaded+1
         END IF
      END DO
      CALL sort(expected_work, force_env%para_env%num_pe, indices=work_index)
      ! Redistribute work in order from the most overloaded processors to the most underloaded processors
      ! Each underloaded processor is limited to one overloaded processor
      IF (load_imbalance(force_env%para_env%mepos+1) > 0) THEN
         offset = 0
         mixed_cdft%dlb_control%send_work = .TRUE.
         ! Build up the total amount of work that needs redistribution
         ALLOCATE (cumulative_work(force_env%para_env%num_pe))
         cumulative_work = 0
         DO i = force_env%para_env%num_pe, force_env%para_env%num_pe-no_overloaded+1, -1
            IF (work_index(i) == force_env%para_env%mepos+1) THEN
               EXIT
            ELSE
               offset = offset+load_imbalance(work_index(i))
               IF (i == force_env%para_env%num_pe) THEN
                  cumulative_work(i) = load_imbalance(work_index(i))
               ELSE
                  cumulative_work(i) = cumulative_work(i+1)+load_imbalance(work_index(i))
               END IF
            END IF
         END DO
         my_pos = i
         j = force_env%para_env%num_pe
         nsend_max = load_imbalance(work_index(j))/work_size(work_index(j))
         exhausted_work = 0
         ! Determine send offset by going through all processors that are more overloaded than my_pos
         DO i = 1, no_underloaded
            IF (my_pos == force_env%para_env%num_pe) EXIT
            nsend = -load_imbalance(work_index(i))/work_size(work_index(j))
            IF (nsend .LT. 1) nsend = 1
            nsend_max = nsend_max-nsend
            IF (nsend_max .LT. 0) nsend = nsend+nsend_max
            exhausted_work = exhausted_work+nsend*work_size(work_index(j))
            offset = offset-nsend*work_size(work_index(j))
            IF (offset .LT. 0) EXIT
            IF (exhausted_work .EQ. cumulative_work(j)) THEN
               j = j-1
               nsend_max = load_imbalance(work_index(j))/work_size(work_index(j))
            END IF
         END DO
         ! Underloaded processors were fully exhausted: rewind index
         ! Load balancing will fail if this happens on multiple processors
         IF (i .GT. no_underloaded) THEN
            i = no_underloaded
         END IF
         my_target = i
         DEALLOCATE (cumulative_work)
         ! Determine how much and who to send slices of my grid points
         nsend_max = load_imbalance(force_env%para_env%mepos+1)/work_size(force_env%para_env%mepos+1)
         ! This the actual number of available array slices
         IF (mixed_cdft%is_pencil .OR. mixed_cdft%is_special) THEN
            nsend_limit = bo(2, 1)-bo(1, 1)+1
         ELSE
            nsend_limit = bo(2, 2)-bo(1, 2)+1
         END IF
         IF (.NOT. mixed_cdft%is_special) THEN
            ALLOCATE (mixed_cdft%dlb_control%target_list(3, max_targets))
         ELSE
            ALLOCATE (mixed_cdft%dlb_control%target_list(3+2*SIZE(mixed_cdft%dest_list), max_targets))
            ALLOCATE (touched(SIZE(mixed_cdft%dest_list)))
            touched = .FALSE.
         END IF
         mixed_cdft%dlb_control%target_list = uninitialized
         i = 1
         ispecial = 1
         offset_special = 0
         targets(1, my_pos) = my_target
         send_total = 0
         ! Main loop. Note, we actually allow my_pos to offload more slices than nsend_max
         DO
            nsend = -load_imbalance(work_index(my_target))/work_size(force_env%para_env%mepos+1)
            IF (nsend .LT. 1) nsend = 1 ! send at least one block
            ! Prevent over redistribution: leave at least (1-work_factor)*nsend_limit slices to my_pos
            IF (nsend .GT. NINT(work_factor*nsend_limit-send_total)) THEN
               nsend = NINT(work_factor*nsend_limit-send_total)
               IF (debug_this_module) &
                  should_warn(force_env%para_env%mepos+1) = 1
            END IF
            mixed_cdft%dlb_control%target_list(1, i) = work_index(my_target)-1 ! This is the actual processor rank
            IF (mixed_cdft%is_special) THEN
               mixed_cdft%dlb_control%target_list(2, i) = 0
               actually_sent = nsend
               DO j = ispecial, SIZE(mixed_cdft%dest_list)
                  mixed_cdft%dlb_control%target_list(2, i) = mixed_cdft%dlb_control%target_list(2, i)+1
                  touched(j) = .TRUE.
                  IF (nsend .LT. mixed_cdft%dest_list_bo(2, j)-mixed_cdft%dest_list_bo(1, j)+1) THEN
                     mixed_cdft%dlb_control%target_list(3+2*j-1, i) = mixed_cdft%dest_list_bo(1, j)
                     mixed_cdft%dlb_control%target_list(3+2*j, i) = mixed_cdft%dest_list_bo(1, j)+nsend-1
                     mixed_cdft%dest_list_bo(1, j) = mixed_cdft%dest_list_bo(1, j)+nsend
                     nsend = 0
                     EXIT
                  ELSE
                     mixed_cdft%dlb_control%target_list(3+2*j-1, i) = mixed_cdft%dest_list_bo(1, j)
                     mixed_cdft%dlb_control%target_list(3+2*j, i) = mixed_cdft%dest_list_bo(2, j)
                     nsend = nsend-(mixed_cdft%dest_list_bo(2, j)-mixed_cdft%dest_list_bo(1, j)+1)
                     mixed_cdft%dest_list_bo(1:2, j) = should_deallocate
                  END IF
                  IF (nsend .LE. 0) EXIT
               END DO
               IF (mixed_cdft%dest_list_bo(1, ispecial) .EQ. should_deallocate) ispecial = j+1
               actually_sent = actually_sent-nsend
               nsend_max = nsend_max-actually_sent
               send_total = send_total+actually_sent
            ELSE
               mixed_cdft%dlb_control%target_list(2, i) = nsend
               nsend_max = nsend_max-nsend
               send_total = send_total+nsend
            END IF
            IF (nsend_max .LT. 0) nsend_max = 0
            IF (nsend_max .EQ. 0) EXIT
            IF (my_target /= no_underloaded) THEN
               my_target = my_target+1
            ELSE
               ! If multiple processors execute this block load balancing will fail
               mixed_cdft%dlb_control%target_list(2, i) = mixed_cdft%dlb_control%target_list(2, i)+nsend_max
               nsend_max = 0
               EXIT
            END IF
            i = i+1
            IF (i .GT. max_targets) &
               CALL cp_abort(__LOCATION__, &
                             "Load balancing error: increase max_targets")
         END DO
         IF (.NOT. mixed_cdft%is_special) THEN
            CALL reallocate(mixed_cdft%dlb_control%target_list, 1, 3, 1, i)
         ELSE
            CALL reallocate(mixed_cdft%dlb_control%target_list, 1, 3+2*SIZE(mixed_cdft%dest_list), 1, i)
         END IF
         targets(2, my_pos) = my_target
         ! Equalize the load on the target processors
         IF (.NOT. mixed_cdft%is_special) THEN
            IF (send_total .GT. NINT(work_factor*nsend_limit)) send_total = NINT(work_factor*nsend_limit)-1
            nsend = NINT(REAL(send_total, dp)/REAL(SIZE(mixed_cdft%dlb_control%target_list, 2), dp))
            mixed_cdft%dlb_control%target_list(2, :) = nsend
         END IF
      ELSE
         DO i = 1, no_underloaded
            IF (work_index(i) == force_env%para_env%mepos+1) EXIT
         END DO
         my_pos = i
      END IF
      CALL mp_sum(targets, force_env%para_env%group)
      IF (debug_this_module) THEN
         CALL mp_sum(should_warn, force_env%para_env%group)
         IF (ANY(should_warn == 1)) &
            CALL cp_warn(__LOCATION__, &
                         "MIXED_CDFT DLB: Attempted to redistribute more array"// &
                         " slices than actually available. Leaving a fraction of the total "// &
                         " slices on the overloaded processor. Perhaps you have set LOAD_SCALE too high?")
         DEALLOCATE (should_warn)
      END IF
      ! check that there is one-to-one mapping between over- and underloaded processors
      IF (force_env%para_env%mepos == force_env%para_env%source) THEN
         consistent = .TRUE.
         DO i = force_env%para_env%num_pe-1, force_env%para_env%num_pe-no_overloaded+1, -1
            IF (targets(1, i) .GT. no_underloaded) consistent = .FALSE.
            IF (targets(1, i) .GT. targets(2, i+1)) THEN
               CYCLE
            ELSE
               consistent = .FALSE.
            END IF
         END DO
         IF (.NOT. consistent) THEN
            IF (debug_this_module .AND. iounit > 0) THEN
               DO i = force_env%para_env%num_pe-1, force_env%para_env%num_pe-no_overloaded+1, -1
                  WRITE (iounit, '(A,I8,I8,I8,I8,I8)') &
                     'load balancing info', load_imbalance(i), work_index(i), &
                     work_size(i), targets(1, i), targets(2, i)
               END DO
            END IF
            CALL cp_abort(__LOCATION__, &
                          "Load balancing error: too much data to redistribute."// &
                          " Increase LOAD_SCALE or change the number of processors."// &
                          " If the confinement cavity occupies a large volume relative"// &
                          " to the total system volume, it might be worth disabling DLB.")
         END IF
      END IF
      ! Tell the target processors which grid points they should compute
      IF (my_pos .LE. no_underloaded) THEN
         DO i = force_env%para_env%num_pe, force_env%para_env%num_pe-no_overloaded+1, -1
            IF (targets(1, i) .LE. my_pos .AND. targets(2, i) .GE. my_pos) THEN
               mixed_cdft%dlb_control%recv_work = .TRUE.
               mixed_cdft%dlb_control%my_source = work_index(i)-1
               EXIT
            END IF
         END DO
         IF (mixed_cdft%dlb_control%recv_work) THEN
            IF (.NOT. mixed_cdft%is_special) THEN
               ALLOCATE (mixed_cdft%dlb_control%bo(12))
               CALL mp_irecv(msgout=mixed_cdft%dlb_control%bo, source=mixed_cdft%dlb_control%my_source, &
                             request=req(1), comm=force_env%para_env%group)
               CALL mp_wait(req(1))
               mixed_cdft%dlb_control%my_dest_repl = (/mixed_cdft%dlb_control%bo(11), mixed_cdft%dlb_control%bo(12)/)
               mixed_cdft%dlb_control%dest_tags_repl = (/mixed_cdft%dlb_control%bo(9), mixed_cdft%dlb_control%bo(10)/)
               ALLOCATE (mixed_cdft%dlb_control%cavity(mixed_cdft%dlb_control%bo(1):mixed_cdft%dlb_control%bo(2), &
                                                       mixed_cdft%dlb_control%bo(3):mixed_cdft%dlb_control%bo(4), &
                                                       mixed_cdft%dlb_control%bo(7):mixed_cdft%dlb_control%bo(8)))
               ALLOCATE (mixed_cdft%dlb_control%weight(mixed_cdft%dlb_control%bo(1):mixed_cdft%dlb_control%bo(2), &
                                                       mixed_cdft%dlb_control%bo(3):mixed_cdft%dlb_control%bo(4), &
                                                       mixed_cdft%dlb_control%bo(7):mixed_cdft%dlb_control%bo(8)))
               ALLOCATE (mixed_cdft%dlb_control%gradients(3*natom, &
                                                          mixed_cdft%dlb_control%bo(1):mixed_cdft%dlb_control%bo(2), &
                                                          mixed_cdft%dlb_control%bo(3):mixed_cdft%dlb_control%bo(4), &
                                                          mixed_cdft%dlb_control%bo(7):mixed_cdft%dlb_control%bo(8)))
               mixed_cdft%dlb_control%gradients = 0.0_dp
               mixed_cdft%dlb_control%weight = 0.0_dp
               CALL mp_irecv(msgout=mixed_cdft%dlb_control%cavity, source=mixed_cdft%dlb_control%my_source, &
                             request=req(1), comm=force_env%para_env%group)
               CALL mp_wait(req(1))
               DEALLOCATE (mixed_cdft%dlb_control%bo)
            ELSE
               ALLOCATE (buffsize(1))
               CALL mp_irecv(msgout=buffsize, source=mixed_cdft%dlb_control%my_source, &
                             request=req(1), comm=force_env%para_env%group)
               CALL mp_wait(req(1))
               ALLOCATE (mixed_cdft%dlb_control%bo(12*buffsize(1)))
               CALL mp_irecv(msgout=mixed_cdft%dlb_control%bo, source=mixed_cdft%dlb_control%my_source, &
                             request=req(1), comm=force_env%para_env%group)
               ALLOCATE (mixed_cdft%dlb_control%sendbuff(buffsize(1)))
               ALLOCATE (req_recv(buffsize(1)))
               DEALLOCATE (buffsize)
               CALL mp_wait(req(1))
               DO j = 1, SIZE(mixed_cdft%dlb_control%sendbuff)
                  ALLOCATE (mixed_cdft%dlb_control%sendbuff(j)%cavity(mixed_cdft%dlb_control%bo(12*(j-1)+1): &
                                                                      mixed_cdft%dlb_control%bo(12*(j-1)+2), &
                                                                      mixed_cdft%dlb_control%bo(12*(j-1)+3): &
                                                                      mixed_cdft%dlb_control%bo(12*(j-1)+4), &
                                                                      mixed_cdft%dlb_control%bo(12*(j-1)+7): &
                                                                      mixed_cdft%dlb_control%bo(12*(j-1)+8)))
                  CALL mp_irecv(msgout=mixed_cdft%dlb_control%sendbuff(j)%cavity, &
                                source=mixed_cdft%dlb_control%my_source, &
                                request=req_recv(j), comm=force_env%para_env%group)
                  ALLOCATE (mixed_cdft%dlb_control%sendbuff(j)%weight(mixed_cdft%dlb_control%bo(12*(j-1)+1): &
                                                                      mixed_cdft%dlb_control%bo(12*(j-1)+2), &
                                                                      mixed_cdft%dlb_control%bo(12*(j-1)+3): &
                                                                      mixed_cdft%dlb_control%bo(12*(j-1)+4), &
                                                                      mixed_cdft%dlb_control%bo(12*(j-1)+7): &
                                                                      mixed_cdft%dlb_control%bo(12*(j-1)+8)))
                  ALLOCATE (mixed_cdft%dlb_control%sendbuff(j)%gradients(3*natom, &
                                                                         mixed_cdft%dlb_control%bo(12*(j-1)+1): &
                                                                         mixed_cdft%dlb_control%bo(12*(j-1)+2), &
                                                                         mixed_cdft%dlb_control%bo(12*(j-1)+3): &
                                                                         mixed_cdft%dlb_control%bo(12*(j-1)+4), &
                                                                         mixed_cdft%dlb_control%bo(12*(j-1)+7): &
                                                                         mixed_cdft%dlb_control%bo(12*(j-1)+8)))
                  mixed_cdft%dlb_control%sendbuff(j)%weight = 0.0_dp
                  mixed_cdft%dlb_control%sendbuff(j)%gradients = 0.0_dp
                  mixed_cdft%dlb_control%sendbuff(j)%tag = (/mixed_cdft%dlb_control%bo(12*(j-1)+9), &
                                                             mixed_cdft%dlb_control%bo(12*(j-1)+10)/)
                  mixed_cdft%dlb_control%sendbuff(j)%rank = (/mixed_cdft%dlb_control%bo(12*(j-1)+11), &
                                                              mixed_cdft%dlb_control%bo(12*(j-1)+12)/)
               END DO
               CALL mp_waitall(req_recv)
               DEALLOCATE (req_recv)
            END IF
         END IF
      ELSE
         IF (.NOT. mixed_cdft%is_special) THEN
            offset = 0
            ALLOCATE (sendbuffer(12))
            send_total = 0
            DO i = 1, SIZE(mixed_cdft%dlb_control%target_list, 2)
               tags = (/(i-1)*3+1+force_env%para_env%mepos*6*max_targets, &
                        (i-1)*3+1+3*max_targets+force_env%para_env%mepos*6*max_targets/) ! Unique communicator tags
               mixed_cdft%dlb_control%target_list(3, i) = tags(1)
               IF (mixed_cdft%is_pencil) THEN
                  sendbuffer = (/bo_conf(1, 1)+offset, &
                                 bo_conf(1, 1)+offset+(mixed_cdft%dlb_control%target_list(2, i)-1), &
                                 bo_conf(1, 2), bo_conf(2, 2), bo(1, 3), bo(2, 3), bo_conf(1, 3), bo_conf(2, 3), &
                                 tags(1), tags(2), mixed_cdft%dest_list(1), mixed_cdft%dest_list(2)/)
               ELSE
                  sendbuffer = (/bo_conf(1, 1), bo_conf(2, 1), &
                                 bo_conf(1, 2)+offset, &
                                 bo_conf(1, 2)+offset+(mixed_cdft%dlb_control%target_list(2, i)-1), &
                                 bo(1, 3), bo(2, 3), bo_conf(1, 3), bo_conf(2, 3), tags(1), tags(2), &
                                 mixed_cdft%dest_list(1), mixed_cdft%dest_list(2)/)
               END IF
               send_total = send_total+mixed_cdft%dlb_control%target_list(2, i)-1
               CALL mp_isend(msgin=sendbuffer, dest=mixed_cdft%dlb_control%target_list(1, i), &
                             request=req(1), comm=force_env%para_env%group)
               CALL mp_wait(req(1))
               IF (mixed_cdft%is_pencil) THEN
                  ALLOCATE (cavity(bo_conf(1, 1)+offset: &
                                   bo_conf(1, 1)+offset+(mixed_cdft%dlb_control%target_list(2, i)-1), &
                                   bo_conf(1, 2):bo_conf(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
                  cavity = mixed_cdft%becke_control%cavity%pw%cr3d(bo_conf(1, 1)+offset: &
                                                                   bo_conf(1, 1)+offset+ &
                                                                   (mixed_cdft%dlb_control%target_list(2, i)-1), &
                                                                   bo_conf(1, 2):bo_conf(2, 2), &
                                                                   bo_conf(1, 3):bo_conf(2, 3))
               ELSE
                  ALLOCATE (cavity(bo_conf(1, 1):bo_conf(2, 1), &
                                   bo_conf(1, 2)+offset: &
                                   bo_conf(1, 2)+offset+(mixed_cdft%dlb_control%target_list(2, i)-1), &
                                   bo_conf(1, 3):bo_conf(2, 3)))
                  cavity = mixed_cdft%becke_control%cavity%pw%cr3d(bo_conf(1, 1):bo_conf(2, 1), &
                                                                   bo_conf(1, 2)+offset: &
                                                                   bo_conf(1, 2)+offset+ &
                                                                   (mixed_cdft%dlb_control%target_list(2, i)-1), &
                                                                   bo_conf(1, 3):bo_conf(2, 3))
               END IF
               CALL mp_isend(msgin=cavity, &
                             dest=mixed_cdft%dlb_control%target_list(1, i), &
                             request=req(1), comm=force_env%para_env%group)
               CALL mp_wait(req(1))
               offset = offset+mixed_cdft%dlb_control%target_list(2, i)
               DEALLOCATE (cavity)
            END DO
            IF (mixed_cdft%is_pencil) THEN
               mixed_cdft%dlb_control%distributed(1) = bo_conf(1, 1)
               mixed_cdft%dlb_control%distributed(2) = bo_conf(1, 1)+offset-1
            ELSE
               mixed_cdft%dlb_control%distributed(1) = bo_conf(1, 2)
               mixed_cdft%dlb_control%distributed(2) = bo_conf(1, 2)+offset-1
            END IF
            DEALLOCATE (sendbuffer)
         ELSE
            ALLOCATE (buffsize(1))
            DO i = 1, SIZE(mixed_cdft%dlb_control%target_list, 2)
               buffsize = mixed_cdft%dlb_control%target_list(2, i)
               ! Unique communicator tags (dont actually need these, should be removed)
               tags = (/(i-1)*3+1+force_env%para_env%mepos*6*max_targets, &
                        (i-1)*3+1+3*max_targets+force_env%para_env%mepos*6*max_targets/)
               DO j = 4, SIZE(mixed_cdft%dlb_control%target_list, 1)
                  IF (mixed_cdft%dlb_control%target_list(j, i) .GT. uninitialized) EXIT
               END DO
               offset_special = j
               offset_proc = j-4-(j-4)/2
               CALL mp_isend(msgin=buffsize, &
                             dest=mixed_cdft%dlb_control%target_list(1, i), &
                             request=req(1), comm=force_env%para_env%group)
               CALL mp_wait(req(1))
               ALLOCATE (sendbuffer(12*buffsize(1)))
               DO j = 1, buffsize(1)
                  sendbuffer(12*(j-1)+1:12*(j-1)+12) = (/mixed_cdft%dlb_control%target_list(offset_special+2*(j-1), i), &
                                                         mixed_cdft%dlb_control%target_list(offset_special+2*j-1, i), &
                                                         bo_conf(1, 2), bo_conf(2, 2), bo(1, 3), bo(2, 3), &
                                                         bo_conf(1, 3), bo_conf(2, 3), tags(1), tags(2), &
                                                         mixed_cdft%dest_list(j+offset_proc), &
                                                         mixed_cdft%dest_list(j+offset_proc)+force_env%para_env%num_pe/2/)
               END DO
               CALL mp_isend(msgin=sendbuffer, &
                             dest=mixed_cdft%dlb_control%target_list(1, i), &
                             request=req(1), comm=force_env%para_env%group)
               CALL mp_wait(req(1))
               DEALLOCATE (sendbuffer)
               DO j = 1, buffsize(1)
                  ALLOCATE (cavity(mixed_cdft%dlb_control%target_list(offset_special+2*(j-1), i): &
                                   mixed_cdft%dlb_control%target_list(offset_special+2*j-1, i), &
                                   bo_conf(1, 2):bo_conf(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
                  cavity = mixed_cdft%becke_control%cavity%pw%cr3d(LBOUND(cavity, 1):UBOUND(cavity, 1), &
                                                                   bo_conf(1, 2):bo_conf(2, 2), &
                                                                   bo_conf(1, 3):bo_conf(2, 3))
                  CALL mp_isend(msgin=cavity, &
                                dest=mixed_cdft%dlb_control%target_list(1, i), &
                                request=req(1), comm=force_env%para_env%group)
                  CALL mp_wait(req(1))
                  DEALLOCATE (cavity)
               END DO
            END DO
            DEALLOCATE (buffsize)
         END IF
      END IF
      DEALLOCATE (expected_work, work_size, load_imbalance, work_index, targets)
      ! Once calculated, data defined on the distributed grid points is sent directly to the processors that own the
      ! grid points after the constraint is copied onto the two processor groups, instead of sending the data back to
      ! the original owner
      IF (mixed_cdft%is_special) THEN
         my_special_work = 2
         ALLOCATE (mask_send(SIZE(mixed_cdft%dest_list)), mask_recv(SIZE(mixed_cdft%source_list)))
         ALLOCATE (nsend_proc(SIZE(mixed_cdft%dest_list)), nrecv(SIZE(mixed_cdft%source_list)))
         nrecv = 0
         nsend_proc = 0
         mask_recv = .FALSE.
         mask_send = .FALSE.
      ELSE
         my_special_work = 1
      END IF
      ALLOCATE (recvbuffer(SIZE(mixed_cdft%source_list)), sbuff(SIZE(mixed_cdft%dest_list)))
      ALLOCATE (req_total(my_special_work*SIZE(mixed_cdft%source_list)+(my_special_work**2)*SIZE(mixed_cdft%dest_list)))
      ALLOCATE (mixed_cdft%dlb_control%recv_work_repl(SIZE(mixed_cdft%source_list)))
      DO i = 1, SIZE(mixed_cdft%source_list)
         NULLIFY (recvbuffer(i)%bv, recvbuffer(i)%iv)
         ALLOCATE (recvbuffer(i)%bv(1), recvbuffer(i)%iv(3))
         CALL mp_irecv(msgout=recvbuffer(i)%bv, &
                       source=mixed_cdft%source_list(i), &
                       request=req_total(i), tag=1, comm=force_env%para_env%group)
         IF (mixed_cdft%is_special) &
            CALL mp_irecv(msgout=recvbuffer(i)%iv, &
                          source=mixed_cdft%source_list(i), &
                          request=req_total(i+SIZE(mixed_cdft%source_list)), &
                          tag=2, comm=force_env%para_env%group)
      END DO
      DO i = 1, my_special_work
         DO j = 1, SIZE(mixed_cdft%dest_list)
            IF (i == 1) THEN
               NULLIFY (sbuff(j)%iv, sbuff(j)%bv)
               ALLOCATE (sbuff(j)%bv(1))
               sbuff(j)%bv = mixed_cdft%dlb_control%send_work
               IF (mixed_cdft%is_special) THEN
                  ALLOCATE (sbuff(j)%iv(3))
                  sbuff(j)%iv(1:2) = mixed_cdft%dest_list_bo(1:2, j)
                  sbuff(j)%iv(3) = 0
                  IF (sbuff(j)%iv(1) .EQ. should_deallocate) mask_send(j) = .TRUE.
                  IF (mixed_cdft%dlb_control%send_work) THEN
                     sbuff(j)%bv = touched(j)
                     IF (touched(j)) THEN
                        nsend = 0
                        DO ispecial = 1, SIZE(mixed_cdft%dlb_control%target_list, 2)
                           IF (mixed_cdft%dlb_control%target_list(4+2*(j-1), ispecial) .NE. uninitialized) &
                              nsend = nsend+1
                        END DO
                        sbuff(j)%iv(3) = nsend
                        nsend_proc(j) = nsend
                     END IF
                  END IF
               END IF
            END IF
            ind = j+(i-1)*SIZE(mixed_cdft%dest_list)+my_special_work*SIZE(mixed_cdft%source_list)
            CALL mp_isend(msgin=sbuff(j)%bv, &
                          dest=mixed_cdft%dest_list(j)+(i-1)*force_env%para_env%num_pe/2, &
                          request=req_total(ind), tag=1, &
                          comm=force_env%para_env%group)
            IF (mixed_cdft%is_special) &
               CALL mp_isend(msgin=sbuff(j)%iv, &
                             dest=mixed_cdft%dest_list(j)+(i-1)*force_env%para_env%num_pe/2, &
                             request=req_total(ind+2*SIZE(mixed_cdft%dest_list)), tag=2, &
                             comm=force_env%para_env%group)
         END DO
      END DO
      CALL mp_waitall(req_total)
      DEALLOCATE (req_total)
      DO i = 1, SIZE(mixed_cdft%source_list)
         mixed_cdft%dlb_control%recv_work_repl(i) = recvbuffer(i)%bv(1)
         IF (mixed_cdft%is_special .AND. mixed_cdft%dlb_control%recv_work_repl(i)) THEN
            mixed_cdft%source_list_bo(1:2, i) = recvbuffer(i)%iv(1:2)
            nrecv(i) = recvbuffer(i)%iv(3)
            IF (recvbuffer(i)%iv(1) .EQ. should_deallocate) mask_recv(i) = .TRUE.
         END IF
         DEALLOCATE (recvbuffer(i)%bv)
         IF (ASSOCIATED(recvbuffer(i)%iv)) DEALLOCATE (recvbuffer(i)%iv)
      END DO
      DO j = 1, SIZE(mixed_cdft%dest_list)
         DEALLOCATE (sbuff(j)%bv)
         IF (ASSOCIATED(sbuff(j)%iv)) DEALLOCATE (sbuff(j)%iv)
      END DO
      DEALLOCATE (recvbuffer)
      ! For some reason if debug_this_module is true and is_special is false, the deallocate statement
      ! on line 3433 gets executed no matter what (gfortran 5.3.0 bug?). Printing out the variable seems to fix it...
      IF (debug_this_module) THEN
         WRITE (dummy, *) mixed_cdft%is_special
      END IF
      IF (.NOT. mixed_cdft%is_special) THEN
         IF (mixed_cdft%dlb_control%send_work) THEN
            ALLOCATE (sendbuffer(6))
            IF (mixed_cdft%is_pencil) THEN
               sendbuffer = (/SIZE(mixed_cdft%dlb_control%target_list, 2), bo_conf(1, 3), bo_conf(2, 3), &
                              bo_conf(1, 1), bo_conf(1, 2), bo_conf(2, 2)/)
            ELSE
               sendbuffer = (/SIZE(mixed_cdft%dlb_control%target_list, 2), bo_conf(1, 3), bo_conf(2, 3), &
                              bo_conf(1, 2), bo_conf(1, 1), bo_conf(2, 1)/)
            END IF
         END IF

         IF (mixed_cdft%dlb_control%recv_work_repl(1) .OR. mixed_cdft%dlb_control%recv_work_repl(2)) THEN
            ALLOCATE (mixed_cdft%dlb_control%recv_info(2))
            NULLIFY (mixed_cdft%dlb_control%recv_info(1)%target_list, mixed_cdft%dlb_control%recv_info(2)%target_list)
            ALLOCATE (mixed_cdft%dlb_control%recvbuff(2))
            NULLIFY (mixed_cdft%dlb_control%recvbuff(1)%buffs, mixed_cdft%dlb_control%recvbuff(2)%buffs)
         END IF
         DO i = 1, 2
            IF (mixed_cdft%dlb_control%recv_work_repl(i) .AND. &
                mixed_cdft%dlb_control%send_work) THEN
               ALLOCATE (mixed_cdft%dlb_control%recv_info(i)%matrix_info(6))
               IF (i == 2) mixed_cdft%dlb_control%target_list(3, :) = mixed_cdft%dlb_control%target_list(3, :)+ &
                                                                      3*max_targets
               CALL mp_irecv(msgout=mixed_cdft%dlb_control%recv_info(i)%matrix_info, &
                             source=mixed_cdft%source_list(i), &
                             request=my_req(1), comm=force_env%para_env%group)
               CALL mp_isend(msgin=sendbuffer, &
                             dest=mixed_cdft%dest_list(i), &
                             request=my_req(2), comm=force_env%para_env%group)
               CALL mp_waitall(my_req)
               ALLOCATE (mixed_cdft%dlb_control%recv_info(i)% &
                         target_list(3, mixed_cdft%dlb_control%recv_info(i)%matrix_info(1)))
               CALL mp_irecv(mixed_cdft%dlb_control%recv_info(i)%target_list, &
                             source=mixed_cdft%source_list(i), &
                             request=my_req(1), comm=force_env%para_env%group)
               CALL mp_isend(msgin=mixed_cdft%dlb_control%target_list, &
                             dest=mixed_cdft%dest_list(i), &
                             request=my_req(2), comm=force_env%para_env%group)
               CALL mp_waitall(my_req)
            ELSE IF (mixed_cdft%dlb_control%recv_work_repl(i)) THEN
               ! Determine which processors will send me data
               ALLOCATE (mixed_cdft%dlb_control%recv_info(i)%matrix_info(6))
               CALL mp_irecv(msgout=mixed_cdft%dlb_control%recv_info(i)%matrix_info, &
                             source=mixed_cdft%source_list(i), &
                             request=req(1), comm=force_env%para_env%group)
               CALL mp_wait(req(1))
               ALLOCATE (mixed_cdft%dlb_control%recv_info(i)% &
                         target_list(3, mixed_cdft%dlb_control%recv_info(i)%matrix_info(1)))
               CALL mp_irecv(mixed_cdft%dlb_control%recv_info(i)%target_list, &
                             source=mixed_cdft%source_list(i), &
                             request=req(1), comm=force_env%para_env%group)
               CALL mp_wait(req(1))
            ELSE IF (mixed_cdft%dlb_control%send_work) THEN
               ! Tell my destination processors which processors will send them data
               IF (i == 2) &
                  mixed_cdft%dlb_control%target_list(3, :) = mixed_cdft%dlb_control%target_list(3, :)+3*max_targets
               CALL mp_isend(msgin=sendbuffer, &
                             dest=mixed_cdft%dest_list(i), &
                             request=req(1), comm=force_env%para_env%group)
               CALL mp_wait(req(1))
               CALL mp_isend(msgin=mixed_cdft%dlb_control%target_list, &
                             dest=mixed_cdft%dest_list(i), &
                             request=req(1), comm=force_env%para_env%group)
               CALL mp_wait(req(1))
            END IF
         END DO
         IF (ASSOCIATED(sendbuffer)) DEALLOCATE (sendbuffer)
      ELSE
         IF (mixed_cdft%dlb_control%send_work) THEN
            ALLOCATE (req_total(COUNT(mixed_cdft%dlb_control%recv_work_repl)+2*COUNT(touched)))
         ELSE IF (ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
            ALLOCATE (req_total(COUNT(mixed_cdft%dlb_control%recv_work_repl)))
         END IF
         IF (mixed_cdft%dlb_control%send_work) THEN
            ind = COUNT(mixed_cdft%dlb_control%recv_work_repl)
            DO j = 1, SIZE(mixed_cdft%dest_list)
               IF (touched(j)) THEN
                  ALLOCATE (sbuff(j)%iv(4+3*nsend_proc(j)))
                  sbuff(j)%iv(1:4) = (/bo_conf(1, 2), bo_conf(2, 2), bo_conf(1, 3), bo_conf(2, 3)/)
                  offset = 5
                  DO i = 1, SIZE(mixed_cdft%dlb_control%target_list, 2)
                     IF (mixed_cdft%dlb_control%target_list(4+2*(j-1), i) .NE. uninitialized) THEN
                        sbuff(j)%iv(offset:offset+2) = (/mixed_cdft%dlb_control%target_list(1, i), &
                                                         mixed_cdft%dlb_control%target_list(4+2*(j-1), i), &
                                                         mixed_cdft%dlb_control%target_list(4+2*j-1, i)/)
                        offset = offset+3
                     END IF
                  END DO
                  DO ispecial = 1, my_special_work
                     CALL mp_isend(msgin=sbuff(j)%iv, &
                                   dest=mixed_cdft%dest_list(j)+(ispecial-1)*force_env%para_env%num_pe/2, &
                                   request=req_total(ind+ispecial), comm=force_env%para_env%group)
                  END DO
                  ind = ind+my_special_work
               END IF
            END DO
         END IF
         IF (ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
            ALLOCATE (mixed_cdft%dlb_control%recv_info(SIZE(mixed_cdft%source_list)))
            ALLOCATE (mixed_cdft%dlb_control%recvbuff(SIZE(mixed_cdft%source_list)))
            ind = 1
            DO j = 1, SIZE(mixed_cdft%source_list)
               NULLIFY (mixed_cdft%dlb_control%recv_info(j)%target_list, &
                        mixed_cdft%dlb_control%recvbuff(j)%buffs)
               IF (mixed_cdft%dlb_control%recv_work_repl(j)) THEN
                  ALLOCATE (mixed_cdft%dlb_control%recv_info(j)%matrix_info(4+3*nrecv(j)))
                  CALL mp_irecv(mixed_cdft%dlb_control%recv_info(j)%matrix_info, &
                                source=mixed_cdft%source_list(j), &
                                request=req_total(ind), comm=force_env%para_env%group)
                  ind = ind+1
               END IF
            END DO
         END IF
         IF (ASSOCIATED(req_total)) THEN
            CALL mp_waitall(req_total)
            DEALLOCATE (req_total)
         END IF
         IF (ANY(mask_send)) THEN
            ALLOCATE (tmp(SIZE(mixed_cdft%dest_list)-COUNT(mask_send)), &
                      tmp_bo(2, SIZE(mixed_cdft%dest_list)-COUNT(mask_send)))
            i = 1
            DO j = 1, SIZE(mixed_cdft%dest_list)
               IF (.NOT. mask_send(j)) THEN
                  tmp(i) = mixed_cdft%dest_list(j)
                  tmp_bo(1:2, i) = mixed_cdft%dest_list_bo(1:2, j)
                  i = i+1
               END IF
            END DO
            DEALLOCATE (mixed_cdft%dest_list, mixed_cdft%dest_list_bo)
            ALLOCATE (mixed_cdft%dest_list(SIZE(tmp)), mixed_cdft%dest_list_bo(2, SIZE(tmp)))
            mixed_cdft%dest_list = tmp
            mixed_cdft%dest_list_bo = tmp_bo
            DEALLOCATE (tmp, tmp_bo)
         END IF
         IF (ANY(mask_recv)) THEN
            ALLOCATE (tmp(SIZE(mixed_cdft%source_list)-COUNT(mask_recv)), &
                      tmp_bo(4, SIZE(mixed_cdft%source_list)-COUNT(mask_recv)))
            i = 1
            DO j = 1, SIZE(mixed_cdft%source_list)
               IF (.NOT. mask_recv(j)) THEN
                  tmp(i) = mixed_cdft%source_list(j)
                  tmp_bo(1:4, i) = mixed_cdft%source_list_bo(1:4, j)
                  i = i+1
               END IF
            END DO
            DEALLOCATE (mixed_cdft%source_list, mixed_cdft%source_list_bo)
            ALLOCATE (mixed_cdft%source_list(SIZE(tmp)), mixed_cdft%source_list_bo(4, SIZE(tmp)))
            mixed_cdft%source_list = tmp
            mixed_cdft%source_list_bo = tmp_bo
            DEALLOCATE (tmp, tmp_bo)
         END IF
         DEALLOCATE (mask_recv, mask_send)
         DEALLOCATE (nsend_proc, nrecv)
         IF (mixed_cdft%dlb_control%send_work) THEN
            DO j = 1, SIZE(mixed_cdft%dest_list)
               IF (touched(j)) DEALLOCATE (sbuff(j)%iv)
            END DO
            IF (ASSOCIATED(touched)) DEALLOCATE (touched)
         END IF
      END IF
      DEALLOCATE (sbuff)
      CALL cp_print_key_finished_output(iounit, logger, force_env_section, &
                                        "MIXED%PRINT%PROGRAM_RUN_INFO")
      CALL timestop(handle)

   END SUBROUTINE mixed_becke_constraint_dlb

! **************************************************************************************************
!> \brief Low level routine to build mixed Becke constraint and gradients
!> \param force_env the force_env that holds the CDFT states
!> \param mixed_cdft container for structures needed to build the mixed CDFT constraint
!> \param in_memory decides whether to build the weight function gradients in parallel before solving
!>                  the CDFT states or later during the SCF procedure of the individual states
!> \param is_constraint a list used to determine which atoms in the system define the constraint
!> \param store_vectors should temporary arrays be stored in memory to accelerate the calculation
!> \param R12 temporary array holding the pairwise atomic distances
!> \param position_vecs temporary array holding the pbc corrected atomic position vectors
!> \param pair_dist_vecs temporary array holding the pairwise displament vectors
!> \param coefficients array that determines how atoms should be summed to form the constraint
!> \param catom temporary array to map the global index of constraint atoms to their position
!>              in a list that holds only constraint atoms
!> \par History
!>       03.2016  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_becke_constraint_low(force_env, mixed_cdft, in_memory, &
                                         is_constraint, store_vectors, R12, position_vecs, &
                                         pair_dist_vecs, coefficients, catom)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      LOGICAL, INTENT(IN)                                :: in_memory
      LOGICAL, ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  :: is_constraint
      LOGICAL, INTENT(IN)                                :: store_vectors
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: R12, position_vecs
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: pair_dist_vecs
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: coefficients
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  :: catom

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_becke_constraint_low', &
         routineP = moduleN//':'//routineN

      INTEGER :: dir, handle, i, iatom, icomm, iforce_eval, ind(3), index, iounit, ip, ispecial, &
         iwork, j, jatom, jcomm, k, lb_index, my_special_work, my_work, natom, nbuffs, &
         nforce_eval, np(3), nsent_total, nskipped, nwork, offset, offset_repl, ub_index
      INTEGER, DIMENSION(:), POINTER                     :: req_recv, work, work_dlb
      INTEGER, DIMENSION(:, :), POINTER                  :: nsent, req_send
      LOGICAL                                            :: completed_recv, should_communicate
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: skip_me
      LOGICAL, ALLOCATABLE, DIMENSION(:, :)              :: completed
      REAL(kind=dp)                                      :: dist1, dist2, dmyexp, my1, my1_homo, &
                                                            myexp, sum_cell_f_all, &
                                                            sum_cell_f_constr, th, tmp_const
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: cell_functions, distances, ds_dR_i, &
                                                            ds_dR_j
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: d_sum_const_dR, d_sum_Pm_dR, &
                                                            distance_vecs, dP_i_dRi
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dP_i_dRj
      REAL(kind=dp), DIMENSION(3)                        :: cell_v, dist_vec, dmy_dR_i, dmy_dR_j, &
                                                            dr, dr1_r2, dr_i_dR, dr_ij_dR, &
                                                            dr_j_dR, grid_p, r, r1, shift
      REAL(kind=dp), DIMENSION(:), POINTER               :: cutoffs
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: cavity, weight
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: gradients
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys_mix
      TYPE(force_env_type), POINTER                      :: force_env_qs
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(section_vals_type), POINTER                   :: force_env_section, print_section

      logger => cp_get_default_logger()
      NULLIFY (work, req_recv, req_send, work_dlb, nsent, cutoffs, cavity, &
               weight, gradients, cell, subsys_mix, force_env_qs, &
               particle_set, particles, auxbas_pw_pool, force_env_section, &
               print_section)
      CALL timeset(routineN, handle)
      nforce_eval = SIZE(force_env%sub_force_env)
      CALL force_env_get(force_env=force_env, force_env_section=force_env_section)
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%PRINT%PROGRAM_RUN_INFO")
      iounit = cp_print_key_unit_nr(logger, print_section, '', extension='.mixedLog')
      IF (.NOT. force_env%mixed_env%do_mixed_qmmm_cdft) THEN
         CALL force_env_get(force_env=force_env, &
                            subsys=subsys_mix, &
                            cell=cell)
         CALL cp_subsys_get(subsys=subsys_mix, &
                            particles=particles, &
                            particle_set=particle_set)
      ELSE
         DO iforce_eval = 1, nforce_eval
            IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
            force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
         END DO
         CALL get_qs_env(force_env_qs%qmmm_env%qs_env, &
                         cp_subsys=subsys_mix, &
                         cell=cell)
         CALL cp_subsys_get(subsys=subsys_mix, &
                            particles=particles, &
                            particle_set=particle_set)
      END IF
      natom = SIZE(particles%els)
      CALL pw_env_get(pw_env=mixed_cdft%pw_env, auxbas_pw_pool=auxbas_pw_pool)
      np = auxbas_pw_pool%pw_grid%npts
      dr = auxbas_pw_pool%pw_grid%dr
      shift = -REAL(MODULO(np, 2), dp)*dr/2.0_dp
      dir = mixed_cdft%becke_control%confine_dir
      lb_index = mixed_cdft%becke_control%confine_bounds_int(1)
      ub_index = mixed_cdft%becke_control%confine_bounds_int(2)
      ALLOCATE (cell_functions(natom), skip_me(natom))
      IF (store_vectors) THEN
         ALLOCATE (distances(natom))
         ALLOCATE (distance_vecs(3, natom))
      END IF
      IF (in_memory) THEN
         ALLOCATE (ds_dR_j(3))
         ALLOCATE (ds_dR_i(3))
         ALLOCATE (d_sum_Pm_dR(3, natom))
         ALLOCATE (d_sum_const_dR(3, natom))
         ALLOCATE (dP_i_dRj(3, natom, natom))
         ALLOCATE (dP_i_dRi(3, natom))
         th = 1.0e-8_dp
      END IF
      IF (mixed_cdft%dlb) THEN
         ALLOCATE (work(force_env%para_env%num_pe), work_dlb(force_env%para_env%num_pe))
         work = 0
         work_dlb = 0
      END IF
      my_work = 1
      my_special_work = 1
      ! Load balancing: allocate storage for receiving buffers and post recv requests
      IF (mixed_cdft%dlb) THEN
         IF (mixed_cdft%dlb_control%recv_work) THEN
            my_work = 2
            IF (.NOT. mixed_cdft%is_special) THEN
               ALLOCATE (req_send(2, 3))
            ELSE
               ALLOCATE (req_send(2, 3*SIZE(mixed_cdft%dlb_control%sendbuff)))
            END IF
         END IF
         IF (ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
            IF (.NOT. mixed_cdft%is_special) THEN
               offset_repl = 0
               IF (mixed_cdft%dlb_control%recv_work_repl(1) .AND. mixed_cdft%dlb_control%recv_work_repl(2)) THEN
                  ALLOCATE (req_recv(3*(SIZE(mixed_cdft%dlb_control%recv_info(1)%target_list, 2)+ &
                                        SIZE(mixed_cdft%dlb_control%recv_info(2)%target_list, 2))))
                  offset_repl = 3*SIZE(mixed_cdft%dlb_control%recv_info(1)%target_list, 2)
               ELSE IF (mixed_cdft%dlb_control%recv_work_repl(1)) THEN
                  ALLOCATE (req_recv(3*(SIZE(mixed_cdft%dlb_control%recv_info(1)%target_list, 2))))
               ELSE
                  ALLOCATE (req_recv(3*(SIZE(mixed_cdft%dlb_control%recv_info(2)%target_list, 2))))
               END IF
            ELSE
               nbuffs = 0
               offset_repl = 1
               DO j = 1, SIZE(mixed_cdft%dlb_control%recv_work_repl)
                  IF (mixed_cdft%dlb_control%recv_work_repl(j)) THEN
                     nbuffs = nbuffs+(SIZE(mixed_cdft%dlb_control%recv_info(j)%matrix_info)-4)/3
                  END IF
               END DO
               ALLOCATE (req_recv(3*nbuffs))
            END IF
            DO j = 1, SIZE(mixed_cdft%dlb_control%recv_work_repl)
               IF (mixed_cdft%dlb_control%recv_work_repl(j)) THEN
                  IF (.NOT. mixed_cdft%is_special) THEN
                     offset = 0
                     index = j+(j/2)
                     ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(SIZE(mixed_cdft%dlb_control%recv_info(j)%target_list, 2)))
                     DO i = 1, SIZE(mixed_cdft%dlb_control%recv_info(j)%target_list, 2)
                        IF (mixed_cdft%is_pencil) THEN
                           ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                     weight(mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)+offset: &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)+offset+ &
                                            (mixed_cdft%dlb_control%recv_info(j)%target_list(2, i)-1), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(5): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(6), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(2): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(3)))
                           ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                     cavity(mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)+offset: &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)+offset+ &
                                            (mixed_cdft%dlb_control%recv_info(j)%target_list(2, i)-1), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(5): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(6), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(2): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(3)))
                           ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                     gradients(3*natom, &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)+offset: &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)+offset+ &
                                               (mixed_cdft%dlb_control%recv_info(j)%target_list(2, i)-1), &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(5): &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(6), &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(2): &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(3)))
                        ELSE
                           ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                     weight(mixed_cdft%dlb_control%recv_info(j)%matrix_info(5): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(6), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)+offset: &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)+offset+ &
                                            (mixed_cdft%dlb_control%recv_info(j)%target_list(2, i)-1), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(2): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(3)))
                           ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                     cavity(mixed_cdft%dlb_control%recv_info(j)%matrix_info(5): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(6), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)+offset: &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)+offset+ &
                                            (mixed_cdft%dlb_control%recv_info(j)%target_list(2, i)-1), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(2): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(3)))
                           ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                     gradients(3*natom, &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(5): &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(6), &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)+offset: &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)+offset+ &
                                               (mixed_cdft%dlb_control%recv_info(j)%target_list(2, i)-1), &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(2): &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(3)))
                        END IF
                        CALL mp_irecv(msgout=mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%cavity, &
                                      source=mixed_cdft%dlb_control%recv_info(j)%target_list(1, i), &
                                      request=req_recv(3*(i-1)+(j-1)*offset_repl+1), &
                                      comm=force_env%para_env%group, &
                                      tag=mixed_cdft%dlb_control%recv_info(j)%target_list(3, i))
                        CALL mp_irecv(msgout=mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight, &
                                      source=mixed_cdft%dlb_control%recv_info(j)%target_list(1, i), &
                                      request=req_recv(3*(i-1)+(j-1)*offset_repl+2), &
                                      comm=force_env%para_env%group, &
                                      tag=mixed_cdft%dlb_control%recv_info(j)%target_list(3, i)+1)
                        CALL mp_irecv(msgout=mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%gradients, &
                                      source=mixed_cdft%dlb_control%recv_info(j)%target_list(1, i), &
                                      request=req_recv(3*(i-1)+(j-1)*offset_repl+3), &
                                      comm=force_env%para_env%group, &
                                      tag=mixed_cdft%dlb_control%recv_info(j)%target_list(3, i)+2)
                        offset = offset+mixed_cdft%dlb_control%recv_info(j)%target_list(2, i)
                     END DO
                     DEALLOCATE (mixed_cdft%dlb_control%recv_info(j)%matrix_info)
                  ELSE
                     ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)% &
                               buffs((SIZE(mixed_cdft%dlb_control%recv_info(j)%matrix_info)-4)/3))
                     index = 6
                     DO i = 1, SIZE(mixed_cdft%dlb_control%recvbuff(j)%buffs)
                        ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                  weight(mixed_cdft%dlb_control%recv_info(j)%matrix_info(index): &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(index+1), &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(1): &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(2), &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(3): &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)))
                        ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                  cavity(mixed_cdft%dlb_control%recv_info(j)%matrix_info(index): &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(index+1), &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(1): &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(2), &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(3): &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)))
                        ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                  gradients(3*natom, mixed_cdft%dlb_control%recv_info(j)%matrix_info(index): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(index+1), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(1): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(2), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(3): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)))
                        CALL mp_irecv(msgout=mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%cavity, &
                                      source=mixed_cdft%dlb_control%recv_info(j)%matrix_info(index-1), &
                                      request=req_recv(offset_repl), &
                                      comm=force_env%para_env%group, tag=1)
                        CALL mp_irecv(msgout=mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight, &
                                      source=mixed_cdft%dlb_control%recv_info(j)%matrix_info(index-1), &
                                      request=req_recv(offset_repl+1), &
                                      comm=force_env%para_env%group, tag=2)
                        CALL mp_irecv(msgout=mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%gradients, &
                                      source=mixed_cdft%dlb_control%recv_info(j)%matrix_info(index-1), &
                                      request=req_recv(offset_repl+2), &
                                      comm=force_env%para_env%group, tag=3)
                        index = index+3
                        offset_repl = offset_repl+3
                     END DO
                     DEALLOCATE (mixed_cdft%dlb_control%recv_info(j)%matrix_info)
                  END IF
               END IF
            END DO
         END IF
      END IF
      cutoffs => mixed_cdft%becke_control%cutoffs
      should_communicate = .FALSE.
      DO i = 1, 3
         cell_v(i) = cell%hmat(i, i)
      END DO
      DO iwork = my_work, 1, -1
         IF (iwork == 2) THEN
            IF (.NOT. mixed_cdft%is_special) THEN
               cavity => mixed_cdft%dlb_control%cavity
               weight => mixed_cdft%dlb_control%weight
               gradients => mixed_cdft%dlb_control%gradients
               ALLOCATE (completed(2, 3), nsent(2, 3))
            ELSE
               my_special_work = SIZE(mixed_cdft%dlb_control%sendbuff)
               ALLOCATE (completed(2, 3*my_special_work), nsent(2, 3*my_special_work))
            END IF
            completed = .FALSE.
            nsent = 0
         ELSE
            IF (.NOT. mixed_cdft%is_special) THEN
               weight => mixed_cdft%weight
               cavity => mixed_cdft%cavity
               gradients => mixed_cdft%becke_control%gradients
            ELSE
               my_special_work = SIZE(mixed_cdft%dest_list)
            END IF
         END IF
         DO ispecial = 1, my_special_work
            nwork = 0
            IF (mixed_cdft%is_special) THEN
               IF (iwork == 1) THEN
                  weight => mixed_cdft%sendbuff(ispecial)%weight
                  cavity => mixed_cdft%sendbuff(ispecial)%cavity
                  gradients => mixed_cdft%sendbuff(ispecial)%gradients
               ELSE
                  weight => mixed_cdft%dlb_control%sendbuff(ispecial)%weight
                  cavity => mixed_cdft%dlb_control%sendbuff(ispecial)%cavity
                  gradients => mixed_cdft%dlb_control%sendbuff(ispecial)%gradients
               END IF
            END IF
            DO k = LBOUND(weight, 1), UBOUND(weight, 1)
               IF (mixed_cdft%dlb .AND. mixed_cdft%is_pencil .AND. .NOT. mixed_cdft%is_special) THEN
                  IF (mixed_cdft%dlb_control%send_work) THEN
                     IF (k .GE. mixed_cdft%dlb_control%distributed(1) .AND. &
                         k .LE. mixed_cdft%dlb_control%distributed(2)) THEN
                        CYCLE
                     END IF
                  END IF
               END IF
               DO j = LBOUND(weight, 2), UBOUND(weight, 2)
                  IF (mixed_cdft%dlb .AND. .NOT. mixed_cdft%is_pencil .AND. .NOT. mixed_cdft%is_special) THEN
                     IF (mixed_cdft%dlb_control%send_work) THEN
                        IF (j .GE. mixed_cdft%dlb_control%distributed(1) .AND. &
                            j .LE. mixed_cdft%dlb_control%distributed(2)) THEN
                           CYCLE
                        END IF
                     END IF
                  END IF
                  ! Check if any of the buffers have become available for deallocation
                  IF (should_communicate) THEN
                     DO icomm = 1, SIZE(nsent, 2)
                        DO jcomm = 1, SIZE(nsent, 1)
                           IF (nsent(jcomm, icomm) == 1) CYCLE
                           CALL mp_test(req_send(jcomm, icomm), completed(jcomm, icomm))
                           IF (completed(jcomm, icomm)) THEN
                              nsent(jcomm, icomm) = nsent(jcomm, icomm)+1
                              nsent_total = nsent_total+1
                              IF (nsent_total == SIZE(nsent, 1)*SIZE(nsent, 2)) should_communicate = .FALSE.
                           END IF
                           IF (ALL(completed(:, icomm))) THEN
                              IF (MODULO(icomm, 3) == 1) THEN
                                 IF (.NOT. mixed_cdft%is_special) THEN
                                    DEALLOCATE (mixed_cdft%dlb_control%cavity)
                                 ELSE
                                    DEALLOCATE (mixed_cdft%dlb_control%sendbuff((icomm-1)/3+1)%cavity)
                                 END IF
                              ELSE IF (MODULO(icomm, 3) == 2) THEN
                                 IF (.NOT. mixed_cdft%is_special) THEN
                                    DEALLOCATE (mixed_cdft%dlb_control%weight)
                                 ELSE
                                    DEALLOCATE (mixed_cdft%dlb_control%sendbuff((icomm-1)/3+1)%weight)
                                 END IF
                              ELSE
                                 IF (.NOT. mixed_cdft%is_special) THEN
                                    DEALLOCATE (mixed_cdft%dlb_control%gradients)
                                 ELSE
                                    DEALLOCATE (mixed_cdft%dlb_control%sendbuff((icomm-1)/3+1)%gradients)
                                 END IF
                              END IF
                           END IF
                        END DO
                     END DO
                  END IF
                  ! Poll to prevent starvation
                  IF (ASSOCIATED(req_recv)) &
                     completed_recv = mp_testall(req_recv)
                  !
                  DO i = LBOUND(weight, 3), UBOUND(weight, 3)
                     IF (mixed_cdft%becke_control%cavity_confine) THEN
                        IF (cavity(k, j, i) < mixed_cdft%becke_control%eps_cavity) CYCLE
                     END IF
                     IF (mixed_cdft%becke_control%confine .OR. &
                         mixed_cdft%becke_control%dynamic_confine) THEN
                        ind = (/k, j, i/)
                        IF (ind(dir) .LE. lb_index .OR. ind(dir) .GE. ub_index) CYCLE
                     END IF
                     grid_p(1) = k*dr(1)+shift(1)
                     grid_p(2) = j*dr(2)+shift(2)
                     grid_p(3) = i*dr(3)+shift(3)
                     nskipped = 0
                     cell_functions = 1.0_dp
                     skip_me = .FALSE.
                     IF (store_vectors) distances = 0.0_dp
                     IF (in_memory) THEN
                        d_sum_Pm_dR = 0.0_dp
                        d_sum_const_dR = 0.0_dp
                        dP_i_dRi = 0.0_dp
                     END IF
                     DO iatom = 1, natom
                        IF (skip_me(iatom)) THEN
                           cell_functions(iatom) = 0.0_dp
                           IF (mixed_cdft%becke_control%should_skip) THEN
                              IF (is_constraint(iatom)) nskipped = nskipped+1
                              IF (nskipped == mixed_cdft%becke_control%natoms) THEN
                                 IF (in_memory) THEN
                                    IF (mixed_cdft%becke_control%cavity_confine) THEN
                                       cavity(k, j, i) = 0.0_dp
                                    ELSE
                                       mixed_cdft%becke_control%skip_list(k, j, i) = .TRUE.
                                    END IF
                                 END IF
                                 EXIT
                              END IF
                           END IF
                           CYCLE
                        END IF
                        IF (store_vectors) THEN
                           IF (distances(iatom) .EQ. 0.0_dp) THEN
                              r = position_vecs(:, iatom)
                              dist_vec = (r-grid_p)-ANINT((r-grid_p)/cell_v)*cell_v
                              dist1 = SQRT(DOT_PRODUCT(dist_vec, dist_vec))
                              distance_vecs(:, iatom) = dist_vec
                              distances(iatom) = dist1
                           ELSE
                              dist_vec = distance_vecs(:, iatom)
                              dist1 = distances(iatom)
                           END IF
                        ELSE
                           r = particle_set(iatom)%r
                           DO ip = 1, 3
                              r(ip) = MODULO(r(ip), cell%hmat(ip, ip))-cell%hmat(ip, ip)/2._dp
                           END DO
                           dist_vec = (r-grid_p)-ANINT((r-grid_p)/cell_v)*cell_v
                           dist1 = SQRT(DOT_PRODUCT(dist_vec, dist_vec))
                        END IF
                        IF (dist1 .LE. cutoffs(iatom)) THEN
                           IF (in_memory) THEN
                              IF (dist1 .LE. th) dist1 = th
                              dr_i_dR(:) = dist_vec(:)/dist1
                           END IF
                           DO jatom = 1, natom
                              IF (jatom .NE. iatom) THEN
                                 IF (jatom < iatom) THEN
                                    IF (.NOT. skip_me(jatom)) CYCLE
                                 END IF
                                 IF (store_vectors) THEN
                                    IF (distances(jatom) .EQ. 0.0_dp) THEN
                                       r1 = position_vecs(:, jatom)
                                       dist_vec = (r1-grid_p)-ANINT((r1-grid_p)/cell_v)*cell_v
                                       dist2 = SQRT(DOT_PRODUCT(dist_vec, dist_vec))
                                       distance_vecs(:, jatom) = dist_vec
                                       distances(jatom) = dist2
                                    ELSE
                                       dist_vec = distance_vecs(:, jatom)
                                       dist2 = distances(jatom)
                                    END IF
                                 ELSE
                                    r1 = particle_set(jatom)%r
                                    DO ip = 1, 3
                                       r1(ip) = MODULO(r1(ip), cell%hmat(ip, ip))-cell%hmat(ip, ip)/2._dp
                                    END DO
                                    dist_vec = (r1-grid_p)-ANINT((r1-grid_p)/cell_v)*cell_v
                                    dist2 = SQRT(DOT_PRODUCT(dist_vec, dist_vec))
                                 END IF
                                 IF (in_memory) THEN
                                    IF (store_vectors) THEN
                                       dr1_r2 = pair_dist_vecs(:, iatom, jatom)
                                    ELSE
                                       dr1_r2 = (r-r1)-ANINT((r-r1)/cell_v)*cell_v
                                    END IF
                                    IF (dist2 .LE. th) dist2 = th
                                    tmp_const = (R12(iatom, jatom)**3)
                                    dr_ij_dR(:) = dr1_r2(:)/tmp_const
                                    !derivativ w.r.t. Rj
                                    dr_j_dR = dist_vec(:)/dist2
                                    dmy_dR_j(:) = -(dr_j_dR(:)/R12(iatom, jatom)-(dist1-dist2)*dr_ij_dR(:))
                                    !derivativ w.r.t. Ri
                                    dmy_dR_i(:) = dr_i_dR(:)/R12(iatom, jatom)-(dist1-dist2)*dr_ij_dR(:)
                                 END IF
                                 my1 = (dist1-dist2)/R12(iatom, jatom)
                                 IF (mixed_cdft%becke_control%adjust) THEN
                                    my1_homo = my1
                                    my1 = my1+ &
                                          mixed_cdft%becke_control%aij(iatom, jatom)*(1.0_dp-my1**2)
                                 END IF
                                 myexp = 1.5_dp*my1-0.5_dp*my1**3
                                 IF (in_memory) THEN
                                    dmyexp = 1.5_dp-1.5_dp*my1**2
                                    tmp_const = (1.5_dp**2)*dmyexp*(1-myexp**2)* &
                                                (1.0_dp-((1.5_dp*myexp-0.5_dp*(myexp**3))**2))

                                    ds_dR_i(:) = -0.5_dp*tmp_const*dmy_dR_i(:)
                                    ds_dR_j(:) = -0.5_dp*tmp_const*dmy_dR_j(:)
                                    IF (mixed_cdft%becke_control%adjust) THEN
                                       tmp_const = 1.0_dp-2.0_dp*my1_homo*mixed_cdft%becke_control%aij(iatom, jatom)
                                       ds_dR_i(:) = ds_dR_i(:)*tmp_const
                                       ds_dR_j(:) = ds_dR_j(:)*tmp_const
                                    END IF
                                 END IF
                                 myexp = 1.5_dp*myexp-0.5_dp*myexp**3
                                 myexp = 1.5_dp*myexp-0.5_dp*myexp**3
                                 tmp_const = 0.5_dp*(1.0_dp-myexp)
                                 cell_functions(iatom) = cell_functions(iatom)*tmp_const
                                 IF (in_memory) THEN
                                    IF (ABS(tmp_const) .LE. th) tmp_const = tmp_const+th
                                    dP_i_dRi(:, iatom) = dP_i_dRi(:, iatom)+ds_dR_i(:)/tmp_const
                                    dP_i_dRj(:, iatom, jatom) = ds_dR_j(:)/tmp_const
                                 END IF

                                 IF (dist2 .LE. cutoffs(jatom)) THEN
                                    tmp_const = 0.5_dp*(1.0_dp+myexp)
                                    cell_functions(jatom) = cell_functions(jatom)*tmp_const
                                    IF (in_memory) THEN
                                       IF (ABS(tmp_const) .LE. th) tmp_const = tmp_const+th
                                       dP_i_dRj(:, jatom, iatom) = -ds_dR_i(:)/tmp_const
                                       dP_i_dRi(:, jatom) = dP_i_dRi(:, jatom)-ds_dR_j(:)/tmp_const
                                    END IF
                                 ELSE
                                    skip_me(jatom) = .TRUE.
                                 END IF
                              END IF
                           END DO
                           IF (in_memory) THEN
                              dP_i_dRi(:, iatom) = cell_functions(iatom)*dP_i_dRi(:, iatom)
                              d_sum_Pm_dR(:, iatom) = d_sum_Pm_dR(:, iatom)+dP_i_dRi(:, iatom)
                              IF (is_constraint(iatom)) &
                                 d_sum_const_dR(:, iatom) = d_sum_const_dR(:, iatom)+dP_i_dRi(:, iatom)* &
                                                            coefficients(iatom)
                              DO jatom = 1, natom
                                 IF (jatom .NE. iatom) THEN
                                    IF (jatom < iatom) THEN
                                       IF (.NOT. skip_me(jatom)) THEN
                                          dP_i_dRj(:, iatom, jatom) = cell_functions(iatom)*dP_i_dRj(:, iatom, jatom)
                                          d_sum_Pm_dR(:, jatom) = d_sum_Pm_dR(:, jatom)+dP_i_dRj(:, iatom, jatom)
                                          IF (is_constraint(iatom)) &
                                             d_sum_const_dR(:, jatom) = d_sum_const_dR(:, jatom)+ &
                                                                        dP_i_dRj(:, iatom, jatom)* &
                                                                        coefficients(iatom)
                                          CYCLE
                                       END IF
                                    END IF
                                    dP_i_dRj(:, iatom, jatom) = cell_functions(iatom)*dP_i_dRj(:, iatom, jatom)
                                    d_sum_Pm_dR(:, jatom) = d_sum_Pm_dR(:, jatom)+dP_i_dRj(:, iatom, jatom)
                                    IF (is_constraint(iatom)) &
                                       d_sum_const_dR(:, jatom) = d_sum_const_dR(:, jatom)+dP_i_dRj(:, iatom, jatom)* &
                                                                  coefficients(iatom)
                                 END IF
                              END DO
                           END IF
                        ELSE
                           cell_functions(iatom) = 0.0_dp
                           skip_me(iatom) = .TRUE.
                           IF (mixed_cdft%becke_control%should_skip) THEN
                              IF (is_constraint(iatom)) nskipped = nskipped+1
                              IF (nskipped == mixed_cdft%becke_control%natoms) THEN
                                 IF (in_memory) THEN
                                    IF (mixed_cdft%becke_control%cavity_confine) THEN
                                       cavity(k, j, i) = 0.0_dp
                                    ELSE
                                       mixed_cdft%becke_control%skip_list(k, j, i) = .TRUE.
                                    END IF
                                 END IF
                                 EXIT
                              END IF
                           END IF
                        END IF
                     END DO
                     IF (nskipped == mixed_cdft%becke_control%natoms) CYCLE
                     sum_cell_f_constr = 0.0_dp
                     DO ip = 1, mixed_cdft%becke_control%natoms
                        sum_cell_f_constr = sum_cell_f_constr+cell_functions(catom(ip))* &
                                            mixed_cdft%becke_control%coeff(ip)
                     END DO
                     sum_cell_f_all = 0.0_dp
                     nwork = nwork+1
                     DO ip = 1, natom
                        sum_cell_f_all = sum_cell_f_all+cell_functions(ip)
                     END DO
                     IF (in_memory) THEN
                        DO iatom = 1, natom
                           IF (ABS(sum_cell_f_all) .GT. 0.0_dp) THEN
                              gradients(3*(iatom-1)+1:3*(iatom-1)+3, k, j, i) = &
                                 d_sum_const_dR(:, iatom)/sum_cell_f_all-sum_cell_f_constr* &
                                 d_sum_Pm_dR(:, iatom)/(sum_cell_f_all**2)
                           END IF
                        END DO
                     END IF
                     IF (ABS(sum_cell_f_all) .GT. 0.000001) &
                        weight(k, j, i) = sum_cell_f_constr/sum_cell_f_all
                  END DO ! i
               END DO ! j
            END DO ! k
            ! Load balancing: post send requests
            IF (iwork == 2) THEN
               IF (.NOT. mixed_cdft%is_special) THEN
                  DO i = 1, SIZE(req_send, 1)
                     CALL mp_isend(msgin=mixed_cdft%dlb_control%cavity, &
                                   dest=mixed_cdft%dlb_control%my_dest_repl(i), &
                                   request=req_send(i, 1), comm=force_env%para_env%group, &
                                   tag=mixed_cdft%dlb_control%dest_tags_repl(i))
                     CALL mp_isend(msgin=mixed_cdft%dlb_control%weight, &
                                   dest=mixed_cdft%dlb_control%my_dest_repl(i), &
                                   request=req_send(i, 2), comm=force_env%para_env%group, &
                                   tag=mixed_cdft%dlb_control%dest_tags_repl(i)+1)
                     CALL mp_isend(msgin=mixed_cdft%dlb_control%gradients, &
                                   dest=mixed_cdft%dlb_control%my_dest_repl(i), &
                                   request=req_send(i, 3), comm=force_env%para_env%group, &
                                   tag=mixed_cdft%dlb_control%dest_tags_repl(i)+2)
                  END DO
                  should_communicate = .TRUE.
                  nsent_total = 0
               ELSE
                  DO i = 1, SIZE(req_send, 1)
                     CALL mp_isend(msgin=mixed_cdft%dlb_control%sendbuff(ispecial)%cavity, &
                                   dest=mixed_cdft%dlb_control%sendbuff(ispecial)%rank(i), &
                                   request=req_send(i, 3*(ispecial-1)+1), &
                                   comm=force_env%para_env%group, tag=1)
                     CALL mp_isend(msgin=mixed_cdft%dlb_control%sendbuff(ispecial)%weight, &
                                   dest=mixed_cdft%dlb_control%sendbuff(ispecial)%rank(i), &
                                   request=req_send(i, 3*(ispecial-1)+2), &
                                   comm=force_env%para_env%group, tag=2)
                     CALL mp_isend(msgin=mixed_cdft%dlb_control%sendbuff(ispecial)%gradients, &
                                   dest=mixed_cdft%dlb_control%sendbuff(ispecial)%rank(i), &
                                   request=req_send(i, 3*(ispecial-1)+3), &
                                   comm=force_env%para_env%group, tag=3)
                  END DO
                  IF (ispecial .EQ. my_special_work) THEN
                     should_communicate = .TRUE.
                     nsent_total = 0
                  END IF
               END IF
               work(mixed_cdft%dlb_control%my_source+1) = work(mixed_cdft%dlb_control%my_source+1)+nwork
               work_dlb(force_env%para_env%mepos+1) = work_dlb(force_env%para_env%mepos+1)+nwork
            ELSE
               IF (mixed_cdft%dlb) work(force_env%para_env%mepos+1) = work(force_env%para_env%mepos+1)+nwork
               IF (mixed_cdft%dlb) work_dlb(force_env%para_env%mepos+1) = work_dlb(force_env%para_env%mepos+1)+nwork
            END IF
         END DO ! ispecial
      END DO ! iwork
      ! Load balancing: wait for communication and deallocate sending buffers
      IF (mixed_cdft%dlb) THEN
         IF (mixed_cdft%dlb_control%recv_work .AND. &
             ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
            IF (should_communicate) THEN
               CALL mp_waitall(req_send)
            END IF
            CALL mp_waitall(req_recv)
            IF (ASSOCIATED(mixed_cdft%dlb_control%cavity)) &
               DEALLOCATE (mixed_cdft%dlb_control%cavity)
            IF (ASSOCIATED(mixed_cdft%dlb_control%weight)) &
               DEALLOCATE (mixed_cdft%dlb_control%weight)
            IF (ASSOCIATED(mixed_cdft%dlb_control%gradients)) &
               DEALLOCATE (mixed_cdft%dlb_control%gradients)
            IF (mixed_cdft%is_special) THEN
               DO j = 1, SIZE(mixed_cdft%dlb_control%sendbuff)
                  IF (ASSOCIATED(mixed_cdft%dlb_control%sendbuff(j)%cavity)) &
                     DEALLOCATE (mixed_cdft%dlb_control%sendbuff(j)%cavity)
                  IF (ASSOCIATED(mixed_cdft%dlb_control%sendbuff(j)%weight)) &
                     DEALLOCATE (mixed_cdft%dlb_control%sendbuff(j)%weight)
                  IF (ASSOCIATED(mixed_cdft%dlb_control%sendbuff(j)%gradients)) &
                     DEALLOCATE (mixed_cdft%dlb_control%sendbuff(j)%gradients)
               END DO
               DEALLOCATE (mixed_cdft%dlb_control%sendbuff)
            END IF
            DEALLOCATE (req_send, req_recv)
         ELSE IF (mixed_cdft%dlb_control%recv_work) THEN
            IF (should_communicate) THEN
               CALL mp_waitall(req_send)
            END IF
            IF (ASSOCIATED(mixed_cdft%dlb_control%cavity)) &
               DEALLOCATE (mixed_cdft%dlb_control%cavity)
            IF (ASSOCIATED(mixed_cdft%dlb_control%weight)) &
               DEALLOCATE (mixed_cdft%dlb_control%weight)
            IF (ASSOCIATED(mixed_cdft%dlb_control%gradients)) &
               DEALLOCATE (mixed_cdft%dlb_control%gradients)
            IF (mixed_cdft%is_special) THEN
               DO j = 1, SIZE(mixed_cdft%dlb_control%sendbuff)
                  IF (ASSOCIATED(mixed_cdft%dlb_control%sendbuff(j)%cavity)) &
                     DEALLOCATE (mixed_cdft%dlb_control%sendbuff(j)%cavity)
                  IF (ASSOCIATED(mixed_cdft%dlb_control%sendbuff(j)%weight)) &
                     DEALLOCATE (mixed_cdft%dlb_control%sendbuff(j)%weight)
                  IF (ASSOCIATED(mixed_cdft%dlb_control%sendbuff(j)%gradients)) &
                     DEALLOCATE (mixed_cdft%dlb_control%sendbuff(j)%gradients)
               END DO
               DEALLOCATE (mixed_cdft%dlb_control%sendbuff)
            END IF
            DEALLOCATE (req_send)
         ELSE IF (ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
            CALL mp_waitall(req_recv)
            DEALLOCATE (req_recv)
         END IF
      END IF
      IF (mixed_cdft%dlb) THEN
         CALL mp_sum(work, force_env%para_env%group)
         CALL mp_sum(work_dlb, force_env%para_env%group)
         IF (.NOT. ASSOCIATED(mixed_cdft%dlb_control%prediction_error)) &
            ALLOCATE (mixed_cdft%dlb_control%prediction_error(force_env%para_env%num_pe))
         mixed_cdft%dlb_control%prediction_error = mixed_cdft%dlb_control%expected_work-work
         IF (debug_this_module .AND. iounit > 0) THEN
            DO i = 1, SIZE(work, 1)
               WRITE (iounit, '(A,I10,I10,I10)') &
                  'Work', work(i), work_dlb(i), mixed_cdft%dlb_control%expected_work(i)
            END DO
         END IF
         DEALLOCATE (work, work_dlb, mixed_cdft%dlb_control%expected_work)
      END IF
      NULLIFY (gradients, weight, cavity)
      IF (ALLOCATED(coefficients)) &
         DEALLOCATE (coefficients)
      IF (in_memory) THEN
         DEALLOCATE (ds_dR_j)
         DEALLOCATE (ds_dR_i)
         DEALLOCATE (d_sum_Pm_dR)
         DEALLOCATE (d_sum_const_dR)
         DEALLOCATE (dP_i_dRj)
         DEALLOCATE (dP_i_dRi)
         NULLIFY (gradients)
         IF (store_vectors) THEN
            DEALLOCATE (pair_dist_vecs)
         END IF
      END IF
      NULLIFY (cutoffs)
      IF (ALLOCATED(is_constraint)) &
         DEALLOCATE (is_constraint)
      DEALLOCATE (catom)
      DEALLOCATE (R12)
      DEALLOCATE (cell_functions)
      DEALLOCATE (skip_me)
      IF (ALLOCATED(completed)) &
         DEALLOCATE (completed)
      IF (ASSOCIATED(nsent)) &
         DEALLOCATE (nsent)
      IF (store_vectors) THEN
         DEALLOCATE (distances)
         DEALLOCATE (distance_vecs)
         DEALLOCATE (position_vecs)
      END IF
      IF (ASSOCIATED(req_send)) &
         DEALLOCATE (req_send)
      IF (ASSOCIATED(req_recv)) &
         DEALLOCATE (req_recv)
      CALL cp_print_key_finished_output(iounit, logger, force_env_section, &
                                        "MIXED%PRINT%PROGRAM_RUN_INFO")
      CALL timestop(handle)

   END SUBROUTINE mixed_becke_constraint_low

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

END MODULE mixed_cdft_methods
