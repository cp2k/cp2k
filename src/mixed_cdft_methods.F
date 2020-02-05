!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
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
   USE cp_array_utils,                  ONLY: cp_1d_i_p_type,&
                                              cp_1d_r_p_type,&
                                              cp_2d_r_p_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_syevd
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_invert,&
                                              cp_fm_transpose
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type,&
                                              cp_fm_write_formatted
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_init_p, dbcsr_p_type, dbcsr_release, &
        dbcsr_release_p, dbcsr_scale, dbcsr_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type,&
                                              use_qmmm,&
                                              use_qmmmx,&
                                              use_qs_force
   USE grid_collocate,                  ONLY: collocate_pgf_product
   USE grid_modify_pab_block,           ONLY: FUNC_AB
   USE hirshfeld_methods,               ONLY: create_shape_function
   USE hirshfeld_types,                 ONLY: hirshfeld_type
   USE input_constants,                 ONLY: &
        becke_cutoff_element, becke_cutoff_global, cdft_alpha_constraint, cdft_beta_constraint, &
        cdft_charge_constraint, cdft_magnetization_constraint, mix_cdft, mixed_cdft_parallel, &
        mixed_cdft_parallel_nobuild, mixed_cdft_serial, outer_scf_becke_constraint
   USE input_section_types,             ONLY: section_get_lval,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
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
   USE mixed_cdft_types,                ONLY: mixed_cdft_result_type_set,&
                                              mixed_cdft_settings_type,&
                                              mixed_cdft_type,&
                                              mixed_cdft_type_create,&
                                              mixed_cdft_work_type_release
   USE mixed_cdft_utils,                ONLY: &
        hfun_zero, map_permutation_to_states, mixed_cdft_assemble_block_diag, &
        mixed_cdft_diagonalize_blocks, mixed_cdft_get_blocks, mixed_cdft_init_structures, &
        mixed_cdft_parse_settings, mixed_cdft_print_couplings, mixed_cdft_read_block_diag, &
        mixed_cdft_redistribute_arrays, mixed_cdft_release_work, mixed_cdft_transfer_settings
   USE mixed_environment_types,         ONLY: get_mixed_env,&
                                              mixed_environment_type,&
                                              set_mixed_env
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_scale
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE
   USE qs_cdft_types,                   ONLY: cdft_control_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_io,                        ONLY: read_mo_set,&
                                              wfn_restart_file_name
   USE qs_mo_methods,                   ONLY: make_basis_simple,&
                                              make_basis_sm
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              mo_set_p_type,&
                                              set_mo_set
   USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                              rs2pw,&
                                              rs_grid_release,&
                                              rs_grid_retain,&
                                              rs_grid_zero,&
                                              rs_pw_transfer
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
!> \brief Initialize a mixed CDFT calculation
!> \param force_env the force_env that holds the CDFT states
!> \param calculate_forces determines if forces should be calculated
!> \par History
!>       01.2016  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_init(force_env, calculate_forces)
      TYPE(force_env_type), POINTER                      :: force_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_init', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: et_freq, handle, iforce_eval, iounit, &
                                                            mixing_type, nforce_eval
      LOGICAL                                            :: explicit, is_parallel, is_qmmm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys_mix
      TYPE(force_env_type), POINTER                      :: force_env_qs
      TYPE(mixed_cdft_settings_type)                     :: settings
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(mixed_environment_type), POINTER              :: mixed_env
      TYPE(particle_list_type), POINTER                  :: particles_mix
      TYPE(section_vals_type), POINTER                   :: force_env_section, mapping_section, &
                                                            md_section, mixed_section, &
                                                            print_section, root_section

      NULLIFY (subsys_mix, force_env_qs, force_env_section, print_section, &
               root_section, mixed_section, md_section, mixed_env, mixed_cdft, &
               mapping_section)

      NULLIFY (settings%grid_span, settings%npts, settings%cutoff, settings%rel_cutoff, &
               settings%spherical, settings%rs_dims, settings%odd, settings%atoms, &
               settings%coeffs, settings%si, settings%sr, &
               settings%cutoffs, settings%radii)

      is_qmmm = .FALSE.
      logger => cp_get_default_logger()
      CPASSERT(ASSOCIATED(force_env))
      nforce_eval = SIZE(force_env%sub_force_env)
      CALL timeset(routineN, handle)
      CALL force_env_get(force_env=force_env, force_env_section=force_env_section)
      mixed_env => force_env%mixed_env
      mixed_section => section_vals_get_subs_vals(force_env_section, "MIXED")
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
      iounit = cp_print_key_unit_nr(logger, print_section, '', extension='.mixedLog')
      ! Check if a mixed CDFT calculation is requested
      CALL section_vals_val_get(mixed_section, "MIXING_TYPE", i_val=mixing_type)
      IF (mixing_type == mix_cdft .AND. .NOT. ASSOCIATED(mixed_env%cdft_control)) THEN
         mixed_env%do_mixed_cdft = .TRUE.
         IF (mixed_env%do_mixed_cdft) THEN
            ! Sanity check
            IF (nforce_eval .LT. 2) &
               CALL cp_abort(__LOCATION__, &
                             "Mixed CDFT calculation requires at least 2 force_evals.")
            mapping_section => section_vals_get_subs_vals(mixed_section, "MAPPING")
            CALL section_vals_get(mapping_section, explicit=explicit)
            ! The sub_force_envs must share the same geometrical structure
            IF (explicit) &
               CPABORT("Please disable section &MAPPING for mixed CDFT calculations")
            CALL section_vals_val_get(mixed_section, "MIXED_CDFT%COUPLING", i_val=et_freq)
            IF (et_freq .LT. 0) THEN
               mixed_env%do_mixed_et = .FALSE.
            ELSE
               mixed_env%do_mixed_et = .TRUE.
               IF (et_freq == 0) THEN
                  mixed_env%et_freq = 1
               ELSE
                  mixed_env%et_freq = et_freq
               END IF
            END IF
            ! Start initializing the mixed_cdft type
            ! First determine if the calculation is pure DFT or QMMM and find the qs force_env
            DO iforce_eval = 1, nforce_eval
               IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
               SELECT CASE (force_env%sub_force_env(iforce_eval)%force_env%in_use)
               CASE (use_qs_force)
                  force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
               CASE (use_qmmm)
                  is_qmmm = .TRUE.
                  ! This is really the container for QMMM
                  force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
               CASE (use_qmmmx)
                  CPABORT("No force mixing allowed for mixed CDFT QM/MM")
               CASE DEFAULT
                  CPASSERT(.FALSE.)
               END SELECT
               CPASSERT(ASSOCIATED(force_env_qs))
            END DO
            ! Get infos about the mixed subsys
            IF (.NOT. is_qmmm) THEN
               CALL force_env_get(force_env=force_env, &
                                  subsys=subsys_mix)
               CALL cp_subsys_get(subsys=subsys_mix, &
                                  particles=particles_mix)
            ELSE
               CALL get_qs_env(force_env_qs%qmmm_env%qs_env, &
                               cp_subsys=subsys_mix)
               CALL cp_subsys_get(subsys=subsys_mix, &
                                  particles=particles_mix)
            END IF
            ! Init mixed_cdft_type
            ALLOCATE (mixed_cdft)
            CALL mixed_cdft_type_create(mixed_cdft)
            mixed_cdft%first_iteration = .TRUE.
            ! Determine what run type to use
            IF (mixed_env%ngroups == 1) THEN
               ! States treated in serial, possibly copying CDFT weight function and gradients from state to state
               mixed_cdft%run_type = mixed_cdft_serial
            ELSE IF (mixed_env%ngroups == 2) THEN
               CALL section_vals_val_get(mixed_section, "MIXED_CDFT%PARALLEL_BUILD", l_val=is_parallel)
               IF (is_parallel) THEN
                  ! Treat states in parallel, build weight function and gradients in parallel before SCF process
                  mixed_cdft%run_type = mixed_cdft_parallel
                  IF (.NOT. nforce_eval == 2) &
                     CALL cp_abort(__LOCATION__, &
                                   "Parallel mode mixed CDFT calculation supports only 2 force_evals.")
               ELSE
                  ! Treat states in parallel, but each states builds its own weight function and gradients
                  mixed_cdft%run_type = mixed_cdft_parallel_nobuild
               END IF
            ELSE
               mixed_cdft%run_type = mixed_cdft_parallel_nobuild
            END IF
            ! Store QMMM flag
            mixed_env%do_mixed_qmmm_cdft = is_qmmm
            ! Setup dynamic load balancing
            CALL section_vals_val_get(mixed_section, "MIXED_CDFT%DLB", l_val=mixed_cdft%dlb)
            mixed_cdft%dlb = mixed_cdft%dlb .AND. calculate_forces ! disable if forces are not needed
            mixed_cdft%dlb = mixed_cdft%dlb .AND. (mixed_cdft%run_type == mixed_cdft_parallel) ! disable if not parallel
            IF (mixed_cdft%dlb) THEN
               ALLOCATE (mixed_cdft%dlb_control)
               NULLIFY (mixed_cdft%dlb_control%weight, mixed_cdft%dlb_control%gradients, &
                        mixed_cdft%dlb_control%cavity, mixed_cdft%dlb_control%target_list, &
                        mixed_cdft%dlb_control%bo, mixed_cdft%dlb_control%expected_work, &
                        mixed_cdft%dlb_control%prediction_error, mixed_cdft%dlb_control%sendbuff, &
                        mixed_cdft%dlb_control%recvbuff, mixed_cdft%dlb_control%recv_work_repl, &
                        mixed_cdft%dlb_control%recv_info)
               CALL section_vals_val_get(mixed_section, "MIXED_CDFT%LOAD_SCALE", &
                                         r_val=mixed_cdft%dlb_control%load_scale)
               CALL section_vals_val_get(mixed_section, "MIXED_CDFT%VERY_OVERLOADED", &
                                         r_val=mixed_cdft%dlb_control%very_overloaded)
               CALL section_vals_val_get(mixed_section, "MIXED_CDFT%MORE_WORK", &
                                         i_val=mixed_cdft%dlb_control%more_work)
            END IF
            ! Metric/Wavefunction overlap method/Lowdin orthogonalization/CDFT-CI
            mixed_cdft%calculate_metric = .FALSE.
            mixed_cdft%wfn_overlap_method = .FALSE.
            mixed_cdft%use_lowdin = .FALSE.
            mixed_cdft%do_ci = .FALSE.
            mixed_cdft%nonortho_coupling = .FALSE.
            mixed_cdft%identical_constraints = .TRUE.
            mixed_cdft%block_diagonalize = .FALSE.
            IF (mixed_env%do_mixed_et) THEN
               CALL section_vals_val_get(mixed_section, "MIXED_CDFT%METRIC", &
                                         l_val=mixed_cdft%calculate_metric)
               CALL section_vals_val_get(mixed_section, "MIXED_CDFT%WFN_OVERLAP", &
                                         l_val=mixed_cdft%wfn_overlap_method)
               CALL section_vals_val_get(mixed_section, "MIXED_CDFT%LOWDIN", &
                                         l_val=mixed_cdft%use_lowdin)
               CALL section_vals_val_get(mixed_section, "MIXED_CDFT%CI", &
                                         l_val=mixed_cdft%do_ci)
               CALL section_vals_val_get(mixed_section, "MIXED_CDFT%NONORTHOGONAL_COUPLING", &
                                         l_val=mixed_cdft%nonortho_coupling)
               CALL section_vals_val_get(mixed_section, "MIXED_CDFT%BLOCK_DIAGONALIZE", &
                                         l_val=mixed_cdft%block_diagonalize)
            END IF
            ! Inversion method
            CALL section_vals_val_get(mixed_section, "MIXED_CDFT%EPS_SVD", r_val=mixed_cdft%eps_svd)
            IF (mixed_cdft%eps_svd .LT. 0.0_dp .OR. mixed_cdft%eps_svd .GT. 1.0_dp) &
               CPABORT("Illegal value for EPS_SVD. Value must be between 0.0 and 1.0.")
            ! MD related settings
            CALL force_env_get(force_env, root_section=root_section)
            md_section => section_vals_get_subs_vals(root_section, "MOTION%MD")
            CALL section_vals_val_get(md_section, "TIMESTEP", r_val=mixed_cdft%sim_dt)
            CALL section_vals_val_get(md_section, "STEP_START_VAL", i_val=mixed_cdft%sim_step)
            mixed_cdft%sim_step = mixed_cdft%sim_step - 1 ! to get the first step correct
            mixed_cdft%sim_dt = cp_unit_from_cp2k(mixed_cdft%sim_dt, "fs")
            ! Parse constraint settings from the individual force_evals and check consistency
            CALL mixed_cdft_parse_settings(force_env, mixed_env, mixed_cdft, &
                                           settings, natom=SIZE(particles_mix%els))
            ! Transfer settings to mixed_cdft
            CALL mixed_cdft_transfer_settings(force_env, mixed_cdft, settings)
            ! Initilize necessary structures
            CALL mixed_cdft_init_structures(force_env, force_env_qs, mixed_env, mixed_cdft, settings)
            ! Write information about the mixed CDFT calculation
            IF (iounit > 0) THEN
               WRITE (iounit, *) ""
               WRITE (iounit, FMT="(T2,A,T71)") &
                  "MIXED_CDFT| Activating mixed CDFT calculation"
               WRITE (iounit, FMT="(T2,A,T71,I10)") &
                  "MIXED_CDFT| Number of CDFT states:                                  ", nforce_eval
               SELECT CASE (mixed_cdft%run_type)
               CASE (mixed_cdft_parallel)
                  WRITE (iounit, FMT="(T2,A,T71)") &
                     "MIXED_CDFT| CDFT states calculation mode: parallel with build"
                  WRITE (iounit, FMT="(T2,A,T71)") &
                     "MIXED_CDFT| Becke constraint is first built using all available processors"
                  WRITE (iounit, FMT="(T2,A,T71)") &
                     "            and then copied to both states with their own processor groups"
               CASE (mixed_cdft_serial)
                  WRITE (iounit, FMT="(T2,A,T71)") &
                     "MIXED_CDFT| CDFT states calculation mode: serial"
                  IF (mixed_cdft%identical_constraints) THEN
                     WRITE (iounit, FMT="(T2,A,T71)") &
                        "MIXED_CDFT| The constraints are built before the SCF procedure of the first"
                     WRITE (iounit, FMT="(T2,A,T71)") &
                        "            CDFT state and subsequently copied to the other states"
                  ELSE
                     WRITE (iounit, FMT="(T2,A,T71)") &
                        "MIXED_CDFT| The constraints are separately built for all CDFT states"
                  END IF
               CASE (mixed_cdft_parallel_nobuild)
                  WRITE (iounit, FMT="(T2,A,T71)") &
                     "MIXED_CDFT| CDFT states calculation mode: parallel without build"
                  WRITE (iounit, FMT="(T2,A,T71)") &
                     "MIXED_CDFT| The constraints are separately built for all CDFT states"
               CASE DEFAULT
                  CPABORT("Unknown mixed CDFT run type.")
               END SELECT
               WRITE (iounit, FMT="(T2,A,T71,L10)") &
                  "MIXED_CDFT| Calculating electronic coupling between states:         ", mixed_env%do_mixed_et
               WRITE (iounit, FMT="(T2,A,T71,L10)") &
                  "MIXED_CDFT| Calculating electronic coupling reliability metric:     ", mixed_cdft%calculate_metric
               WRITE (iounit, FMT="(T2,A,T71,L10)") &
                  "MIXED_CDFT| Configuration interaction (CDFT-CI) was requested:      ", mixed_cdft%do_ci
               WRITE (iounit, FMT="(T2,A,T71,L10)") &
                  "MIXED_CDFT| Block diagonalizing the mixed CDFT Hamiltonian:         ", mixed_cdft%block_diagonalize
               IF (mixed_cdft%run_type == mixed_cdft_parallel) THEN
                  WRITE (iounit, FMT="(T2,A,T71,L10)") &
                     "MIXED_CDFT| Dynamic load balancing enabled:                         ", mixed_cdft%dlb
                  IF (mixed_cdft%dlb) THEN
                     WRITE (iounit, FMT="(T2,A,T71)") "MIXED_CDFT| Dynamic load balancing parameters:"
                     WRITE (iounit, FMT="(T2,A,T71,F10.2)") &
                        "MIXED_CDFT| load_scale:", mixed_cdft%dlb_control%load_scale
                     WRITE (iounit, FMT="(T2,A,T71,F10.2)") &
                        "MIXED_CDFT| very_overloaded:", mixed_cdft%dlb_control%very_overloaded
                     WRITE (iounit, FMT="(T2,A,T71,I10)") &
                        "MIXED_CDFT| more_work:", mixed_cdft%dlb_control%more_work
                  END IF
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
         END IF
      END IF
      CALL cp_print_key_finished_output(iounit, logger, force_env_section, &
                                        "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
      CALL timestop(handle)

   END SUBROUTINE mixed_cdft_init

! **************************************************************************************************
!> \brief Driver routine to handle the build of CDFT weight/gradient in parallel and serial modes
!> \param force_env the force_env that holds the CDFT states
!> \param calculate_forces if forces should be calculated
!> \param iforce_eval index of the currently active CDFT state (serial mode only)
!> \par History
!>       01.2017  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_build_weight(force_env, calculate_forces, iforce_eval)
      TYPE(force_env_type), POINTER                      :: force_env
      LOGICAL, INTENT(IN)                                :: calculate_forces
      INTEGER, INTENT(IN), OPTIONAL                      :: iforce_eval

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_build_weight', &
         routineP = moduleN//':'//routineN

      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft

      NULLIFY (mixed_cdft)
      CPASSERT(ASSOCIATED(force_env))
      CALL get_mixed_env(force_env%mixed_env, cdft_control=mixed_cdft)
      CPASSERT(ASSOCIATED(mixed_cdft))
      IF (.NOT. PRESENT(iforce_eval)) THEN
         SELECT CASE (mixed_cdft%run_type)
         CASE (mixed_cdft_parallel)
            CALL mixed_cdft_build_weight_parallel(force_env, calculate_forces)
         CASE (mixed_cdft_parallel_nobuild)
            CALL mixed_cdft_set_flags(force_env)
         CASE DEFAULT
            ! Do nothing
         END SELECT
      ELSE
         SELECT CASE (mixed_cdft%run_type)
         CASE (mixed_cdft_serial)
            CALL mixed_cdft_transfer_weight(force_env, calculate_forces, iforce_eval)
         CASE DEFAULT
            ! Do nothing
         END SELECT
      END IF

   END SUBROUTINE mixed_cdft_build_weight

! **************************************************************************************************
!> \brief Build CDFT weight and gradient on 2N processors and copy it to two N processor subgroups
!> \param force_env the force_env that holds the CDFT states
!> \param calculate_forces if forces should be calculted
!> \par History
!>       01.2016  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_build_weight_parallel(force_env, calculate_forces)
      TYPE(force_env_type), POINTER                      :: force_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_build_weight_parallel', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, handle2, i, iforce_eval, ind, INDEX(6), iounit, j, lb_min, &
         my_special_work, natom, nforce_eval, recv_offset, ub_max
      INTEGER, DIMENSION(2, 3)                           :: bo
      INTEGER, DIMENSION(:), POINTER                     :: lb, req_total, sendbuffer_i, ub
      REAL(KIND=dp)                                      :: t1, t2
      TYPE(cdft_control_type), POINTER                   :: cdft_control, cdft_control_target
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
         INTEGER                                         :: imap(6)
         INTEGER, DIMENSION(:), &
            POINTER                                      :: iv => null()
         REAL(KIND=dp), POINTER, &
            DIMENSION(:, :, :)                           :: r3 => null()
         REAL(KIND=dp), POINTER, &
            DIMENSION(:, :, :, :)                        :: r4 => null()
      END TYPE buffers
      TYPE(buffers), DIMENSION(:), POINTER               :: recvbuffer

      NULLIFY (subsys_mix, force_env_qs, particles_mix, force_env_section, print_section, &
               mixed_env, mixed_cdft, pw_env, auxbas_pw_pool, mixed_auxbas_pw_pool, &
               qs_env, dft_control, sendbuffer_i, lb, ub, req_total, recvbuffer, &
               cdft_control, cdft_control_target)

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
         SELECT CASE (force_env%sub_force_env(iforce_eval)%force_env%in_use)
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
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
      iounit = cp_print_key_unit_nr(logger, print_section, '', extension='.mixedLog')
      CALL get_mixed_env(mixed_env, cdft_control=mixed_cdft)
      CPASSERT(ASSOCIATED(mixed_cdft))
      cdft_control => mixed_cdft%cdft_control
      CPASSERT(ASSOCIATED(cdft_control))
      ! Calculate the Becke weight function and gradient on all np processors
      CALL pw_env_get(pw_env=mixed_cdft%pw_env, auxbas_pw_pool=mixed_auxbas_pw_pool)
      natom = SIZE(particles_mix%els)
      CALL mixed_becke_constraint(force_env, calculate_forces)
      ! Start replicating the working arrays on both np/2 processor groups
      mixed_cdft%sim_step = mixed_cdft%sim_step + 1
      CALL get_qs_env(qs_env, pw_env=pw_env, dft_control=dft_control)
      cdft_control_target => dft_control%qs_control%cdft_control
      CPASSERT(dft_control%qs_control%cdft)
      CPASSERT(ASSOCIATED(cdft_control_target))
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      bo = auxbas_pw_pool%pw_grid%bounds_local
      ! First determine the size of the arrays along the confinement dir
      IF (mixed_cdft%is_special) THEN
         my_special_work = 2
      ELSE
         my_special_work = 1
      END IF
      ALLOCATE (recvbuffer(SIZE(mixed_cdft%source_list)))
      ALLOCATE (req_total(my_special_work*SIZE(mixed_cdft%dest_list) + SIZE(mixed_cdft%source_list)))
      ALLOCATE (lb(SIZE(mixed_cdft%source_list)), ub(SIZE(mixed_cdft%source_list)))
      IF (cdft_control%becke_control%cavity_confine) THEN
         ! Gaussian confinement => the bounds depend on the processor and need to be communicated
         ALLOCATE (sendbuffer_i(2))
         sendbuffer_i = cdft_control%becke_control%confine_bounds
         DO i = 1, SIZE(mixed_cdft%source_list)
            ALLOCATE (recvbuffer(i)%iv(2))
            CALL mp_irecv(msgout=recvbuffer(i)%iv, source=mixed_cdft%source_list(i), &
                          request=req_total(i), &
                          comm=force_env%para_env%group)
         END DO
         DO i = 1, my_special_work
            DO j = 1, SIZE(mixed_cdft%dest_list)
               ind = j + (i - 1)*SIZE(mixed_cdft%dest_list) + SIZE(mixed_cdft%source_list)
               CALL mp_isend(msgin=sendbuffer_i, &
                             dest=mixed_cdft%dest_list(j) + (i - 1)*force_env%para_env%num_pe/2, &
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
         ind = j + (j/2)
         IF (mixed_cdft%is_special) THEN
            recvbuffer(j)%imap = (/mixed_cdft%source_list_bo(1, j), mixed_cdft%source_list_bo(2, j), &
                                   mixed_cdft%source_list_bo(3, j), mixed_cdft%source_list_bo(4, j), &
                                   lb(j), ub(j)/)
         ELSE IF (mixed_cdft%is_pencil) THEN
            recvbuffer(j)%imap = (/bo(1, 1), bo(2, 1), mixed_cdft%recv_bo(ind), mixed_cdft%recv_bo(ind + 1), lb(j), ub(j)/)
         ELSE
            recvbuffer(j)%imap = (/mixed_cdft%recv_bo(ind), mixed_cdft%recv_bo(ind + 1), bo(1, 2), bo(2, 2), lb(j), ub(j)/)
         END IF
      END DO
      IF (mixed_cdft%dlb .AND. .NOT. mixed_cdft%is_special) THEN
         IF (mixed_cdft%dlb_control%recv_work_repl(1) .OR. mixed_cdft%dlb_control%recv_work_repl(2)) THEN
            DO j = 1, 2
               recv_offset = 0
               IF (mixed_cdft%dlb_control%recv_work_repl(j)) &
                  recv_offset = SUM(mixed_cdft%dlb_control%recv_info(j)%target_list(2, :))
               IF (mixed_cdft%is_pencil) THEN
                  recvbuffer(j)%imap(1) = recvbuffer(j)%imap(1) + recv_offset
               ELSE
                  recvbuffer(j)%imap(3) = recvbuffer(j)%imap(3) + recv_offset
               END IF
            END DO
         END IF
      END IF
      ! Transfer the arrays one-by-one and deallocate shared storage
      ! Start with the weight function
      DO j = 1, SIZE(mixed_cdft%source_list)
         ALLOCATE (recvbuffer(j)%r3(recvbuffer(j)%imap(1):recvbuffer(j)%imap(2), &
                                    recvbuffer(j)%imap(3):recvbuffer(j)%imap(4), &
                                    recvbuffer(j)%imap(5):recvbuffer(j)%imap(6)))

         CALL mp_irecv(msgout=recvbuffer(j)%r3, source=mixed_cdft%source_list(j), &
                       request=req_total(j), comm=force_env%para_env%group)
      END DO
      DO i = 1, my_special_work
         DO j = 1, SIZE(mixed_cdft%dest_list)
            ind = j + (i - 1)*SIZE(mixed_cdft%dest_list) + SIZE(mixed_cdft%source_list)
            IF (mixed_cdft%is_special) THEN
               CALL mp_isend(msgin=mixed_cdft%sendbuff(j)%weight, &
                             dest=mixed_cdft%dest_list(j) + (i - 1)*force_env%para_env%num_pe/2, &
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
      CALL pw_pool_create_pw(auxbas_pw_pool, cdft_control_target%group(1)%weight%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      cdft_control_target%group(1)%weight%pw%cr3d = 0.0_dp
      ! Assemble the recved slices
      DO j = 1, SIZE(mixed_cdft%source_list)
         cdft_control_target%group(1)%weight%pw%cr3d(recvbuffer(j)%imap(1):recvbuffer(j)%imap(2), &
                                                     recvbuffer(j)%imap(3):recvbuffer(j)%imap(4), &
                                                     recvbuffer(j)%imap(5):recvbuffer(j)%imap(6)) = recvbuffer(j)%r3
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
                     cdft_control_target%group(1)%weight%pw%cr3d(INDEX(1):INDEX(2), &
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
      IF (cdft_control%becke_control%cavity_confine) THEN
         DO j = 1, SIZE(mixed_cdft%source_list)
            CALL mp_irecv(msgout=recvbuffer(j)%r3, source=mixed_cdft%source_list(j), &
                          request=req_total(j), comm=force_env%para_env%group)
         END DO
         DO i = 1, my_special_work
            DO j = 1, SIZE(mixed_cdft%dest_list)
               ind = j + (i - 1)*SIZE(mixed_cdft%dest_list) + SIZE(mixed_cdft%source_list)
               IF (mixed_cdft%is_special) THEN
                  CALL mp_isend(msgin=mixed_cdft%sendbuff(j)%cavity, &
                                dest=mixed_cdft%dest_list(j) + (i - 1)*force_env%para_env%num_pe/2, &
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
         ALLOCATE (cdft_control_target%becke_control%cavity_mat(bo(1, 1):bo(2, 1), &
                                                                bo(1, 2):bo(2, 2), &
                                                                lb_min:ub_max))
         cdft_control_target%becke_control%cavity_mat = 0.0_dp

         DO j = 1, SIZE(mixed_cdft%source_list)
            cdft_control_target%becke_control%cavity_mat(recvbuffer(j)%imap(1):recvbuffer(j)%imap(2), &
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
                        cdft_control_target%becke_control%cavity_mat(INDEX(1):INDEX(2), &
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
               ind = j + (i - 1)*SIZE(mixed_cdft%dest_list) + SIZE(mixed_cdft%source_list)
               IF (mixed_cdft%is_special) THEN
                  CALL mp_isend(msgin=mixed_cdft%sendbuff(j)%gradients, &
                                dest=mixed_cdft%dest_list(j) + (i - 1)*force_env%para_env%num_pe/2, &
                                request=req_total(ind), &
                                comm=force_env%para_env%group)
               ELSE
                  CALL mp_isend(msgin=cdft_control%group(1)%gradients, dest=mixed_cdft%dest_list(j), &
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
            DEALLOCATE (cdft_control%group(1)%gradients)
         END IF
         ALLOCATE (cdft_control_target%group(1)%gradients(3*natom, bo(1, 1):bo(2, 1), &
                                                          bo(1, 2):bo(2, 2), lb_min:ub_max))
         DO j = 1, SIZE(mixed_cdft%source_list)
            cdft_control_target%group(1)%gradients(:, recvbuffer(j)%imap(1):recvbuffer(j)%imap(2), &
                                                   recvbuffer(j)%imap(3):recvbuffer(j)%imap(4), &
                                                   recvbuffer(j)%imap(5):recvbuffer(j)%imap(6)) = recvbuffer(j)%r4
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
                        cdft_control_target%group(1)%gradients(:, INDEX(1):INDEX(2), &
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
      cdft_control_target%external_control = .TRUE.
      cdft_control_target%need_pot = .FALSE.
      cdft_control_target%transfer_pot = .FALSE.
      ! Store the bound indices for force calculation
      IF (calculate_forces) THEN
         cdft_control_target%becke_control%confine_bounds(2) = ub_max
         cdft_control_target%becke_control%confine_bounds(1) = lb_min
      END IF
      CALL pw_scale(cdft_control_target%group(1)%weight%pw, &
                    cdft_control_target%group(1)%weight%pw%pw_grid%dvol)
      ! Set flags for ET coupling calculation
      IF (mixed_env%do_mixed_et) THEN
         IF (MODULO(mixed_cdft%sim_step, mixed_env%et_freq) == 0) THEN
            dft_control%qs_control%cdft_control%do_et = .TRUE.
            dft_control%qs_control%cdft_control%calculate_metric = mixed_cdft%calculate_metric
         ELSE
            dft_control%qs_control%cdft_control%do_et = .FALSE.
            dft_control%qs_control%cdft_control%calculate_metric = .FALSE.
         END IF
      END IF
      t2 = m_walltime()
      IF (iounit > 0) THEN
         WRITE (iounit, '(A)') ' '
         WRITE (iounit, '(T2,A,F6.1,A)') 'MIXED_CDFT| Becke constraint built in ', t2 - t1, ' seconds'
         WRITE (iounit, '(A)') ' '
      END IF
      CALL cp_print_key_finished_output(iounit, logger, force_env_section, &
                                        "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
      CALL timestop(handle)

   END SUBROUTINE mixed_cdft_build_weight_parallel

! **************************************************************************************************
!> \brief Transfer CDFT weight/gradient between force_evals
!> \param force_env the force_env that holds the CDFT sub_force_envs
!> \param calculate_forces if forces should be computed
!> \param iforce_eval index of the currently active CDFT state
!> \par History
!>       01.2017  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_transfer_weight(force_env, calculate_forces, iforce_eval)
      TYPE(force_env_type), POINTER                      :: force_env
      LOGICAL, INTENT(IN)                                :: calculate_forces
      INTEGER, INTENT(IN)                                :: iforce_eval

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_transfer_weight', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: bounds_of(8), handle, iatom, igroup, &
                                                            jforce_eval, nforce_eval
      LOGICAL, SAVE                                      :: first_call = .TRUE.
      TYPE(cdft_control_type), POINTER                   :: cdft_control_source, cdft_control_target
      TYPE(dft_control_type), POINTER                    :: dft_control_source, dft_control_target
      TYPE(force_env_type), POINTER                      :: force_env_qs_source, force_env_qs_target
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(mixed_environment_type), POINTER              :: mixed_env
      TYPE(pw_env_type), POINTER                         :: pw_env_source, pw_env_target
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool_source, &
                                                            auxbas_pw_pool_target
      TYPE(qs_environment_type), POINTER                 :: qs_env_source, qs_env_target

      NULLIFY (mixed_cdft, dft_control_source, dft_control_target, force_env_qs_source, &
               force_env_qs_target, pw_env_source, pw_env_target, auxbas_pw_pool_source, &
               auxbas_pw_pool_target, qs_env_source, qs_env_target, mixed_env, &
               cdft_control_source, cdft_control_target)
      mixed_env => force_env%mixed_env
      CALL get_mixed_env(mixed_env, cdft_control=mixed_cdft)
      CALL timeset(routineN, handle)
      IF (iforce_eval == 1) THEN
         jforce_eval = 1
      ELSE
         jforce_eval = iforce_eval - 1
      END IF
      nforce_eval = SIZE(force_env%sub_force_env)
      SELECT CASE (force_env%sub_force_env(jforce_eval)%force_env%in_use)
      CASE (use_qs_force, use_qmmm)
         force_env_qs_source => force_env%sub_force_env(jforce_eval)%force_env
         force_env_qs_target => force_env%sub_force_env(iforce_eval)%force_env
      CASE DEFAULT
         CPASSERT(.FALSE.)
      END SELECT
      IF (.NOT. force_env%mixed_env%do_mixed_qmmm_cdft) THEN
         CALL force_env_get(force_env=force_env_qs_source, &
                            qs_env=qs_env_source)
         CALL force_env_get(force_env=force_env_qs_target, &
                            qs_env=qs_env_target)
      ELSE
         qs_env_source => force_env_qs_source%qmmm_env%qs_env
         qs_env_target => force_env_qs_target%qmmm_env%qs_env
      END IF
      IF (iforce_eval == 1) THEN
         ! The first force_eval builds the weight function and gradients in qs_cdft_methods.F
         ! Set some flags so the constraint is saved if the constraint definitions are identical in all CDFT states
         CALL get_qs_env(qs_env_source, dft_control=dft_control_source)
         cdft_control_source => dft_control_source%qs_control%cdft_control
         cdft_control_source%external_control = .FALSE.
         cdft_control_source%need_pot = .TRUE.
         IF (mixed_cdft%identical_constraints) THEN
            cdft_control_source%transfer_pot = .TRUE.
         ELSE
            cdft_control_source%transfer_pot = .FALSE.
         END IF
         mixed_cdft%sim_step = mixed_cdft%sim_step + 1
      ELSE
         ! Transfer the constraint from the ith force_eval to the i+1th
         CALL get_qs_env(qs_env_source, dft_control=dft_control_source, &
                         pw_env=pw_env_source)
         CALL pw_env_get(pw_env_source, auxbas_pw_pool=auxbas_pw_pool_source)
         cdft_control_source => dft_control_source%qs_control%cdft_control
         CALL get_qs_env(qs_env_target, dft_control=dft_control_target, &
                         pw_env=pw_env_target)
         CALL pw_env_get(pw_env_target, auxbas_pw_pool=auxbas_pw_pool_target)
         cdft_control_target => dft_control_target%qs_control%cdft_control
         ! The constraint can be transfered only when the constraint defitions are identical in all CDFT states
         IF (mixed_cdft%identical_constraints) THEN
            ! Weight function
            DO igroup = 1, SIZE(cdft_control_target%group)
               CALL pw_pool_create_pw(auxbas_pw_pool_target, &
                                      cdft_control_target%group(igroup)%weight%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               ! We have ensured that the grids are consistent => no danger in using explicit copy
               cdft_control_target%group(igroup)%weight%pw%cr3d = cdft_control_source%group(igroup)%weight%pw%cr3d
               CALL pw_pool_give_back_pw(auxbas_pw_pool_source, cdft_control_source%group(igroup)%weight%pw)
            END DO
            ! Cavity
            IF (cdft_control_source%type == outer_scf_becke_constraint) THEN
               IF (cdft_control_source%becke_control%cavity_confine) THEN
                  CALL pw_pool_create_pw(auxbas_pw_pool_target, cdft_control_target%becke_control%cavity%pw, &
                                         use_data=REALDATA3D, in_space=REALSPACE)
                  cdft_control_target%becke_control%cavity%pw%cr3d(:, :, :) = &
                     cdft_control_source%becke_control%cavity%pw%cr3d
                  CALL pw_pool_give_back_pw(auxbas_pw_pool_source, cdft_control_source%becke_control%cavity%pw)
               END IF
            END IF
            ! Gradients
            IF (calculate_forces) THEN
               DO igroup = 1, SIZE(cdft_control_source%group)
                  bounds_of = (/LBOUND(cdft_control_source%group(igroup)%gradients, 1), &
                                UBOUND(cdft_control_source%group(igroup)%gradients, 1), &
                                LBOUND(cdft_control_source%group(igroup)%gradients, 2), &
                                UBOUND(cdft_control_source%group(igroup)%gradients, 2), &
                                LBOUND(cdft_control_source%group(igroup)%gradients, 3), &
                                UBOUND(cdft_control_source%group(igroup)%gradients, 3), &
                                LBOUND(cdft_control_source%group(igroup)%gradients, 4), &
                                UBOUND(cdft_control_source%group(igroup)%gradients, 4)/)
                  ALLOCATE (cdft_control_target%group(igroup)% &
                            gradients(bounds_of(1):bounds_of(2), bounds_of(3):bounds_of(4), &
                                      bounds_of(5):bounds_of(6), bounds_of(7):bounds_of(8)))
                  cdft_control_target%group(igroup)%gradients = cdft_control_source%group(igroup)%gradients
                  DEALLOCATE (cdft_control_source%group(igroup)%gradients)
               END DO
            END IF
            ! Atomic weight functions needed for CDFT charges
            IF (cdft_control_source%atomic_charges) THEN
               IF (.NOT. ASSOCIATED(cdft_control_target%charge)) &
                  ALLOCATE (cdft_control_target%charge(cdft_control_target%natoms))
               DO iatom = 1, cdft_control_target%natoms
                  CALL pw_pool_create_pw(auxbas_pw_pool_target, &
                                         cdft_control_target%charge(iatom)%pw, &
                                         use_data=REALDATA3D, in_space=REALSPACE)
                  cdft_control_target%charge(iatom)%pw%cr3d = cdft_control_source%charge(iatom)%pw%cr3d
                  CALL pw_pool_give_back_pw(auxbas_pw_pool_source, cdft_control_source%charge(iatom)%pw)
               END DO
            END IF
            ! Set some flags so the weight is not rebuilt during SCF
            cdft_control_target%external_control = .FALSE.
            cdft_control_target%need_pot = .FALSE.
            ! For states i+1 < nforce_eval, prevent deallocation of constraint
            IF (iforce_eval == nforce_eval) THEN
               cdft_control_target%transfer_pot = .FALSE.
            ELSE
               cdft_control_target%transfer_pot = .TRUE.
            END IF
            cdft_control_target%first_iteration = .FALSE.
         ELSE
            ! Force rebuild of constraint and dont save it
            cdft_control_target%external_control = .FALSE.
            cdft_control_target%need_pot = .TRUE.
            cdft_control_target%transfer_pot = .FALSE.
            IF (first_call) THEN
               cdft_control_target%first_iteration = .TRUE.
            ELSE
               cdft_control_target%first_iteration = .FALSE.
            END IF
         END IF
      END IF
      ! Set flags for ET coupling calculation
      IF (mixed_env%do_mixed_et) THEN
         IF (MODULO(mixed_cdft%sim_step, mixed_env%et_freq) == 0) THEN
            IF (iforce_eval == 1) THEN
               cdft_control_source%do_et = .TRUE.
               cdft_control_source%calculate_metric = mixed_cdft%calculate_metric
            ELSE
               cdft_control_target%do_et = .TRUE.
               cdft_control_target%calculate_metric = mixed_cdft%calculate_metric
            END IF
         ELSE
            IF (iforce_eval == 1) THEN
               cdft_control_source%do_et = .FALSE.
               cdft_control_source%calculate_metric = .FALSE.
            ELSE
               cdft_control_target%do_et = .FALSE.
               cdft_control_target%calculate_metric = .FALSE.
            END IF
         END IF
      END IF
      IF (iforce_eval == nforce_eval .AND. first_call) first_call = .FALSE.
      CALL timestop(handle)

   END SUBROUTINE mixed_cdft_transfer_weight

! **************************************************************************************************
!> \brief In case CDFT states are treated in paralle, sets flags so that each CDFT state
!>        builds their own weight functions and gradients
!> \param force_env the force_env that holds the CDFT sub_force_envs
!> \par History
!>       09.2018  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_set_flags(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_set_flags', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iforce_eval, nforce_eval
      LOGICAL, SAVE                                      :: first_call = .TRUE.
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(force_env_type), POINTER                      :: force_env_qs
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(mixed_environment_type), POINTER              :: mixed_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      NULLIFY (mixed_cdft, dft_control, force_env_qs, qs_env, mixed_env, cdft_control)
      mixed_env => force_env%mixed_env
      CALL get_mixed_env(mixed_env, cdft_control=mixed_cdft)
      CALL timeset(routineN, handle)
      nforce_eval = SIZE(force_env%sub_force_env)
      DO iforce_eval = 1, nforce_eval
         IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
         SELECT CASE (force_env%sub_force_env(iforce_eval)%force_env%in_use)
         CASE (use_qs_force, use_qmmm)
            force_env_qs => force_env%sub_force_env(iforce_eval)%force_env
         CASE DEFAULT
            CPASSERT(.FALSE.)
         END SELECT
         IF (.NOT. force_env%mixed_env%do_mixed_qmmm_cdft) THEN
            CALL force_env_get(force_env=force_env_qs, qs_env=qs_env)
         ELSE
            qs_env => force_env_qs%qmmm_env%qs_env
         END IF
         ! All force_evals build the weight function and gradients in qs_cdft_methods.F
         ! Update flags to match run type
         CALL get_qs_env(qs_env, dft_control=dft_control)
         cdft_control => dft_control%qs_control%cdft_control
         cdft_control%external_control = .FALSE.
         cdft_control%need_pot = .TRUE.
         cdft_control%transfer_pot = .FALSE.
         IF (first_call) THEN
            cdft_control%first_iteration = .TRUE.
         ELSE
            cdft_control%first_iteration = .FALSE.
         END IF
         ! Set flags for ET coupling calculation
         IF (mixed_env%do_mixed_et) THEN
            IF (MODULO(mixed_cdft%sim_step, mixed_env%et_freq) == 0) THEN
               cdft_control%do_et = .TRUE.
               cdft_control%calculate_metric = mixed_cdft%calculate_metric
            ELSE
               cdft_control%do_et = .FALSE.
               cdft_control%calculate_metric = .FALSE.
            END IF
         END IF
      END DO
      mixed_cdft%sim_step = mixed_cdft%sim_step + 1
      IF (first_call) first_call = .FALSE.
      CALL timestop(handle)

   END SUBROUTINE mixed_cdft_set_flags

! **************************************************************************************************
!> \brief Driver routine to calculate the electronic coupling(s) between CDFT states.
!> \param force_env the force_env that holds the CDFT states
!> \par History
!>       02.15  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_calculate_coupling(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_calculate_coupling', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CPASSERT(ASSOCIATED(force_env))
      CALL timeset(routineN, handle)
      ! Move needed arrays from individual CDFT states to the mixed CDFT env
      CALL mixed_cdft_redistribute_arrays(force_env)
      ! Calculate the mixed CDFT Hamiltonian and overlap matrices.
      ! All work matrices defined in the wavefunction basis get deallocated on exit.
      ! Any analyses which depend on these work matrices are performed within.
      CALL mixed_cdft_interaction_matrices(force_env)
      ! Calculate eletronic couplings between states (Lowdin/rotation)
      CALL mixed_cdft_calculate_coupling_low(force_env)
      ! Print out couplings
      CALL mixed_cdft_print_couplings(force_env)
      ! Block diagonalize the mixed CDFT Hamiltonian matrix
      CALL mixed_cdft_block_diag(force_env)
      ! CDFT Configuration Interaction
      CALL mixed_cdft_configuration_interaction(force_env)
      ! Clean up
      CALL mixed_cdft_release_work(force_env)
      CALL timestop(handle)

   END SUBROUTINE mixed_cdft_calculate_coupling

! **************************************************************************************************
!> \brief Routine to calculate the mixed CDFT Hamiltonian and overlap matrices.
!> \param force_env the force_env that holds the CDFT states
!> \par History
!>       11.17  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_interaction_matrices(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_interaction_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER :: check_ao(2), check_mo(2), handle, iforce_eval, ipermutation, ispin, istate, ivar, &
         j, jstate, k, moeigvalunit, mounit, nao, ncol_local, nforce_eval, nmo, npermutations, &
         nrow_local, nspins, nvar
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ncol_mo, nrow_mo
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: homo
      LOGICAL                                            :: nelectron_mismatch, print_mo, &
                                                            print_mo_eigval, should_scale, &
                                                            uniform_occupation
      REAL(KIND=dp)                                      :: c(2), eps_occupied, nelectron_tot, &
                                                            sum_a(2), sum_b(2)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coupling_nonortho, eigenv, energy, Sda
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: H_mat, S_det, S_mat, strength, tmp_mat, &
                                                            W_diagonal, Wad, Wda
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: a, b
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigval
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mo_overlap
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: mixed_mo_coeff
      TYPE(cp_fm_p_type), DIMENSION(:, :, :), POINTER    :: w_matrix_mo
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_mo, mo_mo_fmstruct
      TYPE(cp_fm_type), POINTER                          :: inverse_mat, Tinverse, tmp2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: density_matrix, density_matrix_diff, &
                                                            w_matrix
      TYPE(dbcsr_type), POINTER                          :: mixed_matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(mixed_environment_type), POINTER              :: mixed_env
      TYPE(qs_energy_type), POINTER                      :: energy_qs
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: force_env_section, mixed_cdft_section, &
                                                            print_section

      NULLIFY (force_env_section, print_section, mixed_cdft_section, &
               mixed_env, mixed_cdft, qs_env, dft_control, fm_struct_mo, &
               density_matrix_diff, mo_mo_fmstruct, inverse_mat, mo_overlap, &
               Tinverse, tmp2, w_matrix_mo, mixed_mo_coeff, mixed_matrix_s, &
               density_matrix, energy_qs, w_matrix, mo_eigval)
      logger => cp_get_default_logger()
      CPASSERT(ASSOCIATED(force_env))
      CALL timeset(routineN, handle)
      CALL force_env_get(force_env=force_env, &
                         force_env_section=force_env_section)
      mixed_env => force_env%mixed_env
      nforce_eval = SIZE(force_env%sub_force_env)
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
      IF (section_get_lval(print_section, "MO_OVERLAP_MATRIX")) THEN
         print_mo = .TRUE.
         mounit = cp_print_key_unit_nr(logger, print_section, extension='.moOverlap', on_file=.TRUE.)
      ELSE
         print_mo = .FALSE.
      END IF
      IF (section_get_lval(print_section, "MO_OVERLAP_EIGENVALUES")) THEN
         print_mo_eigval = .TRUE.
         moeigvalunit = cp_print_key_unit_nr(logger, print_section, extension='.moOverlapEigval', on_file=.TRUE.)
      ELSE
         print_mo_eigval = .FALSE.
      END IF
      CALL get_mixed_env(mixed_env, cdft_control=mixed_cdft)
      ! Get redistributed work matrices
      CPASSERT(ASSOCIATED(mixed_cdft))
      CPASSERT(ASSOCIATED(mixed_cdft%matrix%mixed_mo_coeff))
      CPASSERT(ASSOCIATED(mixed_cdft%matrix%w_matrix))
      CPASSERT(ASSOCIATED(mixed_cdft%matrix%mixed_matrix_s))
      mixed_mo_coeff => mixed_cdft%matrix%mixed_mo_coeff
      w_matrix => mixed_cdft%matrix%w_matrix
      mixed_matrix_s => mixed_cdft%matrix%mixed_matrix_s
      IF (mixed_cdft%calculate_metric) THEN
         CPASSERT(ASSOCIATED(mixed_cdft%matrix%density_matrix))
         density_matrix => mixed_cdft%matrix%density_matrix
      END IF
      ! Get number of weight functions per state
      nvar = SIZE(w_matrix, 2)
      nspins = SIZE(mixed_mo_coeff, 2)
      ! Check that the number of MOs/AOs is equal in every CDFT state
      ALLOCATE (nrow_mo(nspins), ncol_mo(nspins))
      DO ispin = 1, nspins
         CALL cp_fm_get_info(mixed_mo_coeff(1, ispin)%matrix, ncol_global=check_mo(1), nrow_global=check_ao(1))
         DO iforce_eval = 2, nforce_eval
            CALL cp_fm_get_info(mixed_mo_coeff(iforce_eval, ispin)%matrix, ncol_global=check_mo(2), nrow_global=check_ao(2))
            IF (check_ao(1) /= check_ao(2)) &
               CALL cp_abort(__LOCATION__, &
                             "The number of atomic orbitals must be the same in every CDFT state.")
            IF (check_mo(1) /= check_mo(2)) &
               CALL cp_abort(__LOCATION__, &
                             "The number of molecular orbitals must be the same in every CDFT state.")
         END DO
      END DO
      ! Allocate work
      npermutations = nforce_eval*(nforce_eval - 1)/2 ! Size of upper triangular part
      ALLOCATE (w_matrix_mo(nforce_eval, nforce_eval, nvar))
      ALLOCATE (mo_overlap(npermutations), S_det(npermutations, nspins))
      ALLOCATE (a(nspins, nvar, npermutations), b(nspins, nvar, npermutations))
      a = 0.0_dp
      b = 0.0_dp
      DO istate = 1, npermutations
         NULLIFY (mo_overlap(istate)%matrix)
      END DO
      IF (mixed_cdft%calculate_metric) THEN
         ALLOCATE (density_matrix_diff(npermutations, nspins))
         DO ispin = 1, nspins
            DO ipermutation = 1, npermutations
               NULLIFY (density_matrix_diff(ipermutation, ispin)%matrix)
               CALL dbcsr_init_p(density_matrix_diff(ipermutation, ispin)%matrix)
               CALL map_permutation_to_states(nforce_eval, ipermutation, istate, jstate)
               CALL dbcsr_copy(density_matrix_diff(ipermutation, ispin)%matrix, &
                               density_matrix(istate, ispin)%matrix, name="DENSITY_MATRIX")
            END DO
         END DO
      END IF
      ! Check for uniform occupations
      uniform_occupation = .NOT. ALLOCATED(mixed_cdft%occupations)
      should_scale = .FALSE.
      IF (.NOT. uniform_occupation) THEN
         ALLOCATE (homo(nforce_eval, nspins))
         mixed_cdft_section => section_vals_get_subs_vals(force_env_section, "MIXED%MIXED_CDFT")
         CALL section_vals_val_get(mixed_cdft_section, "EPS_OCCUPIED", r_val=eps_occupied)
         IF (eps_occupied .GT. 1.0_dp .OR. eps_occupied .LT. 0.0_dp) &
            CALL cp_abort(__LOCATION__, &
                          "Keyword EPS_OCCUPIED only accepts values between 0.0 and 1.0")
         IF (mixed_cdft%eps_svd == 0.0_dp) &
            CALL cp_warn(__LOCATION__, &
                         "The usage of SVD based matrix inversions with fractionally occupied "// &
                         "orbitals is strongly recommended to screen nearly orthogonal states.")
         CALL section_vals_val_get(mixed_cdft_section, "SCALE_WITH_OCCUPATION_NUMBERS", l_val=should_scale)
      END IF
      ! Start the actual calculation
      DO ispin = 1, nspins
         ! Create the MOxMO fm struct (mo_mo_fm_pools%struct)
         ! The number of MOs/AOs is equal to the number of columns/rows of mo_coeff(:,:)%matrix
         NULLIFY (fm_struct_mo, mo_mo_fmstruct)
         CALL cp_fm_get_info(mixed_mo_coeff(1, ispin)%matrix, ncol_global=ncol_mo(ispin), nrow_global=nrow_mo(ispin))
         nao = nrow_mo(ispin)
         IF (uniform_occupation) THEN
            nmo = ncol_mo(ispin)
         ELSE
            nmo = ncol_mo(ispin)
            ! Find indices of highest (fractionally) occupied molecular orbital
            homo(:, ispin) = nmo
            DO istate = 1, nforce_eval
               DO j = nmo, 1, -1
                  IF (mixed_cdft%occupations(istate, ispin)%array(j) .GE. eps_occupied) THEN
                     homo(istate, ispin) = j
                     EXIT
                  END IF
               END DO
            END DO
            ! Make matrices square by using the largest homo and emit warning if a state has fewer occupied MOs
            ! Although it would be possible to handle the nonsquare situation as well,
            ! all CDFT states should be in the same spin state for meaningful results
            nmo = MAXVAL(homo(:, ispin))
            ! Also check that the number of electrons is conserved (using a fixed sensible threshold)
            nelectron_mismatch = .FALSE.
            nelectron_tot = SUM(mixed_cdft%occupations(1, ispin)%array(1:nmo))
            DO istate = 2, nforce_eval
               IF (ABS(SUM(mixed_cdft%occupations(istate, ispin)%array(1:nmo)) - nelectron_tot) .GT. 1.0E-4_dp) &
                  nelectron_mismatch = .TRUE.
            END DO
            IF (ANY(homo(:, ispin) /= nmo)) THEN
               IF (ispin == 1) THEN
                  CALL cp_warn(__LOCATION__, &
                               "The number of occupied alpha MOs is not constant across all CDFT states. "// &
                               "Calculation proceeds but the results will likely be meaningless.")
               ELSE
                  CALL cp_warn(__LOCATION__, &
                               "The number of occupied beta MOs is not constant across all CDFT states. "// &
                               "Calculation proceeds but the results will likely be meaningless.")
               END IF
            ELSE IF (nelectron_mismatch) THEN
               IF (ispin == 1) THEN
                  CALL cp_warn(__LOCATION__, &
                               "The number of alpha electrons is not constant across all CDFT states. "// &
                               "Calculation proceeds but the results will likely be meaningless.")
               ELSE
                  CALL cp_warn(__LOCATION__, &
                               "The number of beta electrons is not constant across all CDFT states. "// &
                               "Calculation proceeds but the results will likely be meaningless.")
               END IF
            END IF
         END IF
         CALL cp_fm_struct_create(fm_struct_mo, nrow_global=nao, ncol_global=nmo, &
                                  context=mixed_cdft%blacs_env, para_env=force_env%para_env)
         CALL cp_fm_struct_create(mo_mo_fmstruct, nrow_global=nmo, ncol_global=nmo, &
                                  context=mixed_cdft%blacs_env, para_env=force_env%para_env)
         ! Allocate work
         CALL cp_fm_create(matrix=tmp2, matrix_struct=fm_struct_mo, &
                           name="ET_TMP_"//TRIM(ADJUSTL(cp_to_string(ispin)))//"_MATRIX")
         CALL cp_fm_struct_release(fm_struct_mo)
         CALL cp_fm_create(matrix=inverse_mat, matrix_struct=mo_mo_fmstruct, &
                           name="INVERSE_"//TRIM(ADJUSTL(cp_to_string(ispin)))//"_MATRIX")
         CALL cp_fm_create(matrix=Tinverse, matrix_struct=mo_mo_fmstruct, &
                           name="T_INVERSE_"//TRIM(ADJUSTL(cp_to_string(ispin)))//"_MATRIX")
         DO istate = 1, npermutations
            CALL cp_fm_create(matrix=mo_overlap(istate)%matrix, matrix_struct=mo_mo_fmstruct, &
                              name="MO_OVERLAP_"//TRIM(ADJUSTL(cp_to_string(istate)))//"_"// &
                              TRIM(ADJUSTL(cp_to_string(ispin)))//"_MATRIX")
         END DO
         DO ivar = 1, nvar
            DO istate = 1, nforce_eval
               DO jstate = 1, nforce_eval
                  NULLIFY (w_matrix_mo(istate, jstate, ivar)%matrix)
                  IF (istate == jstate) CYCLE
                  CALL cp_fm_create(matrix=w_matrix_mo(istate, jstate, ivar)%matrix, matrix_struct=mo_mo_fmstruct, &
                                    name="W_"//TRIM(ADJUSTL(cp_to_string(istate)))//"_"// &
                                    TRIM(ADJUSTL(cp_to_string(jstate)))//"_"// &
                                    TRIM(ADJUSTL(cp_to_string(ivar)))//"_MATRIX")
               END DO
            END DO
         END DO
         CALL cp_fm_struct_release(mo_mo_fmstruct)
         ! Remove empty MOs and (possibly) scale rest with occupation numbers
         IF (.NOT. uniform_occupation) THEN
            DO iforce_eval = 1, nforce_eval
               CALL cp_fm_to_fm(mixed_mo_coeff(iforce_eval, ispin)%matrix, tmp2, nmo, 1, 1)
               CALL cp_fm_release(mixed_mo_coeff(iforce_eval, ispin)%matrix)
               CALL cp_fm_create(mixed_mo_coeff(iforce_eval, ispin)%matrix, &
                                 matrix_struct=tmp2%matrix_struct, &
                                 name="MO_COEFF_"//TRIM(ADJUSTL(cp_to_string(iforce_eval)))//"_" &
                                 //TRIM(ADJUSTL(cp_to_string(ispin)))//"_MATRIX")
               CALL cp_fm_to_fm(tmp2, mixed_mo_coeff(iforce_eval, ispin)%matrix)
               IF (should_scale) &
                  CALL cp_fm_column_scale(mixed_mo_coeff(iforce_eval, ispin)%matrix, &
                                          mixed_cdft%occupations(iforce_eval, ispin)%array(1:nmo))
               DEALLOCATE (mixed_cdft%occupations(iforce_eval, ispin)%array)
            END DO
         END IF
         ! calculate the MO overlaps (C_j)^T S C_i
         ipermutation = 0
         DO istate = 1, nforce_eval
            DO jstate = istate + 1, nforce_eval
               ipermutation = ipermutation + 1
               CALL cp_dbcsr_sm_fm_multiply(mixed_matrix_s, mixed_mo_coeff(istate, ispin)%matrix, &
                                            tmp2, nmo, 1.0_dp, 0.0_dp)
               CALL cp_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, &
                            mixed_mo_coeff(jstate, ispin)%matrix, &
                            tmp2, 0.0_dp, mo_overlap(ipermutation)%matrix)
               IF (print_mo) &
                  CALL cp_fm_write_formatted(mo_overlap(ipermutation)%matrix, mounit, &
                                             "# MO overlap matrix (step "//TRIM(ADJUSTL(cp_to_string(mixed_cdft%sim_step)))// &
                                             "): CDFT states "//TRIM(ADJUSTL(cp_to_string(istate)))//" and "// &
                                             TRIM(ADJUSTL(cp_to_string(jstate)))//" (spin "// &
                                             TRIM(ADJUSTL(cp_to_string(ispin)))//")")
            END DO
         END DO
         ! calculate the MO-representations of the restraint matrices of all CDFT states
         DO ivar = 1, nvar
            DO jstate = 1, nforce_eval
               DO istate = 1, nforce_eval
                  IF (istate == jstate) CYCLE
                  ! State i: (C_j)^T W_i C_i
                  CALL cp_dbcsr_sm_fm_multiply(w_matrix(istate, ivar)%matrix, &
                                               mixed_mo_coeff(istate, ispin)%matrix, &
                                               tmp2, nmo, 1.0_dp, 0.0_dp)
                  CALL cp_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, &
                               mixed_mo_coeff(jstate, ispin)%matrix, &
                               tmp2, 0.0_dp, w_matrix_mo(istate, jstate, ivar)%matrix)
               END DO
            END DO
         END DO
         DO ipermutation = 1, npermutations
            ! Invert and calculate determinant of MO overlaps
            CALL map_permutation_to_states(nforce_eval, ipermutation, istate, jstate)
            IF (print_mo_eigval) THEN
               NULLIFY (mo_eigval)
               CALL cp_fm_invert(mo_overlap(ipermutation)%matrix, inverse_mat, &
                                 S_det(ipermutation, ispin), eps_svd=mixed_cdft%eps_svd, &
                                 eigval=mo_eigval)
               IF (moeigvalunit > 0) THEN
                  IF (mixed_cdft%eps_svd == 0.0_dp) THEN
                     WRITE (moeigvalunit, '(A,I2,A,I2,A,I1,A)') &
                        "# MO Overlap matrix eigenvalues for CDFT states ", istate, " and ", jstate, &
                        " (spin ", ispin, ")"
                  ELSE
                     WRITE (moeigvalunit, '(A,I2,A,I2,A,I1,A)') &
                        "# MO Overlap matrix singular values for CDFT states ", istate, " and ", jstate, &
                        " (spin ", ispin, ")"
                  END IF
                  WRITE (moeigvalunit, '(A1, A9, A12)') "#", "Index", ADJUSTL("Value")
                  DO j = 1, SIZE(mo_eigval)
                     WRITE (moeigvalunit, '(I10, F12.8)') j, mo_eigval(j)
                  END DO
               END IF
               DEALLOCATE (mo_eigval)
            ELSE
               CALL cp_fm_invert(mo_overlap(ipermutation)%matrix, inverse_mat, &
                                 S_det(ipermutation, ispin), eps_svd=mixed_cdft%eps_svd)
            END IF
            CALL cp_fm_get_info(inverse_mat, nrow_local=nrow_local, ncol_local=ncol_local)
            ! Calculate <Psi_i | w_j(r) | Psi_j> for ivar:th constraint
            DO j = 1, ncol_local
               DO k = 1, nrow_local
                  DO ivar = 1, nvar
                     b(ispin, ivar, ipermutation) = b(ispin, ivar, ipermutation) + &
                                                    w_matrix_mo(jstate, istate, ivar)%matrix%local_data(k, j)* &
                                                    inverse_mat%local_data(k, j)
                  END DO
               END DO
            END DO
            ! Calculate <Psi_j | w_i(r) | Psi_i> for ivar:th constraint
            CALL cp_fm_transpose(inverse_mat, Tinverse)
            DO j = 1, ncol_local
               DO k = 1, nrow_local
                  DO ivar = 1, nvar
                     a(ispin, ivar, ipermutation) = a(ispin, ivar, ipermutation) + &
                                                    w_matrix_mo(istate, jstate, ivar)%matrix%local_data(k, j)* &
                                                    Tinverse%local_data(k, j)
                  END DO
               END DO
            END DO
            ! Handle different constraint types
            DO ivar = 1, nvar
               SELECT CASE (mixed_cdft%constraint_type(ivar, istate))
               CASE (cdft_charge_constraint)
                  ! No action needed
               CASE (cdft_magnetization_constraint)
                  IF (ispin == 2) a(ispin, ivar, ipermutation) = -a(ispin, ivar, ipermutation)
               CASE (cdft_alpha_constraint)
                  ! Constraint applied to alpha electrons only, set integrals involving beta to zero
                  IF (ispin == 2) a(ispin, ivar, ipermutation) = 0.0_dp
               CASE (cdft_beta_constraint)
                  ! Constraint applied to beta electrons only, set integrals involving alpha to zero
                  IF (ispin == 1) a(ispin, ivar, ipermutation) = 0.0_dp
               CASE DEFAULT
                  CPABORT("Unknown constraint type.")
               END SELECT
               SELECT CASE (mixed_cdft%constraint_type(ivar, jstate))
               CASE (cdft_charge_constraint)
                  ! No action needed
               CASE (cdft_magnetization_constraint)
                  IF (ispin == 2) b(ispin, ivar, ipermutation) = -b(ispin, ivar, ipermutation)
               CASE (cdft_alpha_constraint)
                  ! Constraint applied to alpha electrons only, set integrals involving beta to zero
                  IF (ispin == 2) b(ispin, ivar, ipermutation) = 0.0_dp
               CASE (cdft_beta_constraint)
                  ! Constraint applied to beta electrons only, set integrals involving alpha to zero
                  IF (ispin == 1) b(ispin, ivar, ipermutation) = 0.0_dp
               CASE DEFAULT
                  CPABORT("Unknown constraint type.")
               END SELECT
            END DO
            ! Compute density matrix difference P = P_j - P_i
            IF (mixed_cdft%calculate_metric) &
               CALL dbcsr_add(density_matrix_diff(ipermutation, ispin)%matrix, &
                              density_matrix(jstate, ispin)%matrix, -1.0_dp, 1.0_dp)
            !
            CALL mp_sum(a(ispin, :, ipermutation), force_env%para_env%group)
            CALL mp_sum(b(ispin, :, ipermutation), force_env%para_env%group)
         END DO
         ! Release work
         CALL cp_fm_release(tmp2)
         DO ivar = 1, nvar
            DO jstate = 1, nforce_eval
               DO istate = 1, nforce_eval
                  IF (istate == jstate) CYCLE
                  CALL cp_fm_release(w_matrix_mo(istate, jstate, ivar)%matrix)
               END DO
            END DO
         END DO
         DO ipermutation = 1, npermutations
            CALL cp_fm_release(mo_overlap(ipermutation)%matrix)
         END DO
         CALL cp_fm_release(Tinverse)
         CALL cp_fm_release(inverse_mat)
      END DO
      DEALLOCATE (mo_overlap)
      DEALLOCATE (w_matrix_mo)
      IF (.NOT. uniform_occupation) THEN
         DEALLOCATE (homo)
         DEALLOCATE (mixed_cdft%occupations)
      END IF
      IF (print_mo) &
         CALL cp_print_key_finished_output(mounit, logger, force_env_section, &
                                           "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO", on_file=.TRUE.)
      IF (print_mo_eigval) &
         CALL cp_print_key_finished_output(moeigvalunit, logger, force_env_section, &
                                           "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO", on_file=.TRUE.)
      ! solve eigenstates for the projector matrix
      ALLOCATE (Wda(nvar, npermutations))
      ALLOCATE (Sda(npermutations))
      IF (.NOT. mixed_cdft%identical_constraints) ALLOCATE (Wad(nvar, npermutations))
      DO ipermutation = 1, npermutations
         IF (nspins == 2) THEN
            Sda(ipermutation) = ABS(S_det(ipermutation, 1)*S_det(ipermutation, 2))
         ELSE
            Sda(ipermutation) = S_det(ipermutation, 1)**2
         END IF
         ! Finalize <Psi_j | w_i(r) | Psi_i> by multiplication with Sda
         DO ivar = 1, nvar
            IF (mixed_cdft%identical_constraints) THEN
               Wda(ivar, ipermutation) = (SUM(a(:, ivar, ipermutation)) + SUM(b(:, ivar, ipermutation)))* &
                                         Sda(ipermutation)/2.0_dp
            ELSE
               Wda(ivar, ipermutation) = SUM(a(:, ivar, ipermutation))*Sda(ipermutation)
               Wad(ivar, ipermutation) = SUM(b(:, ivar, ipermutation))*Sda(ipermutation)
            END IF
         END DO
      END DO
      DEALLOCATE (a, b, S_det)
      ! Transfer info about the constraint calculations
      ALLOCATE (W_diagonal(nvar, nforce_eval), strength(nvar, nforce_eval), energy(nforce_eval))
      W_diagonal = 0.0_dp
      DO iforce_eval = 1, nforce_eval
         strength(:, iforce_eval) = mixed_env%strength(iforce_eval, :)
      END DO
      energy = 0.0_dp
      DO iforce_eval = 1, nforce_eval
         IF (.NOT. ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
         IF (force_env%mixed_env%do_mixed_qmmm_cdft) THEN
            qs_env => force_env%sub_force_env(iforce_eval)%force_env%qmmm_env%qs_env
         ELSE
            CALL force_env_get(force_env%sub_force_env(iforce_eval)%force_env, qs_env=qs_env)
         END IF
         CALL get_qs_env(qs_env, energy=energy_qs, dft_control=dft_control)
         IF (force_env%sub_force_env(iforce_eval)%force_env%para_env%mepos == &
             force_env%sub_force_env(iforce_eval)%force_env%para_env%source) THEN
            W_diagonal(:, iforce_eval) = dft_control%qs_control%cdft_control%value(:)
            energy(iforce_eval) = energy_qs%total
         END IF
      END DO
      CALL mp_sum(W_diagonal, force_env%para_env%group)
      CALL mp_sum(energy, force_env%para_env%group)
      CALL mixed_cdft_result_type_set(mixed_cdft%results, Wda=Wda, W_diagonal=W_diagonal, &
                                      energy=energy, strength=strength)
      IF (.NOT. mixed_cdft%identical_constraints) CALL mixed_cdft_result_type_set(mixed_cdft%results, Wad=Wad)
      ! Construct S
      ALLOCATE (S_mat(nforce_eval, nforce_eval))
      DO istate = 1, nforce_eval
         S_mat(istate, istate) = 1.0_dp
      END DO
      DO ipermutation = 1, npermutations
         CALL map_permutation_to_states(nforce_eval, ipermutation, istate, jstate)
         S_mat(istate, jstate) = Sda(ipermutation)
         S_mat(jstate, istate) = Sda(ipermutation)
      END DO
      CALL mixed_cdft_result_type_set(mixed_cdft%results, S=S_mat)
      ! Invert S via eigendecomposition and compute S^-(1/2)
      ALLOCATE (eigenv(nforce_eval), tmp_mat(nforce_eval, nforce_eval))
      CALL diamat_all(S_mat, eigenv, .TRUE.)
      tmp_mat = 0.0_dp
      DO istate = 1, nforce_eval
         IF (eigenv(istate) .LT. 1.0e-14_dp) THEN
            ! Safeguard against division with 0 and negative numbers
            eigenv(istate) = 1.0e-14_dp
            CALL cp_warn(__LOCATION__, &
                         "The overlap matrix is numerically nearly singular. "// &
                         "Calculation proceeds but the results might be meaningless.")
         END IF
         tmp_mat(istate, istate) = 1.0_dp/SQRT(eigenv(istate))
      END DO
      tmp_mat(:, :) = MATMUL(tmp_mat, TRANSPOSE(S_mat))
      S_mat(:, :) = MATMUL(S_mat, tmp_mat) ! S^(-1/2)
      CALL mixed_cdft_result_type_set(mixed_cdft%results, S_minushalf=S_mat)
      DEALLOCATE (eigenv, tmp_mat, S_mat)
      ! Construct nonorthogonal diabatic Hamiltonian matrix H''
      ALLOCATE (H_mat(nforce_eval, nforce_eval))
      IF (mixed_cdft%nonortho_coupling) ALLOCATE (coupling_nonortho(npermutations))
      DO istate = 1, nforce_eval
         H_mat(istate, istate) = energy(istate)
      END DO
      DO ipermutation = 1, npermutations
         CALL map_permutation_to_states(nforce_eval, ipermutation, istate, jstate)
         sum_a = 0.0_dp
         sum_b = 0.0_dp
         DO ivar = 1, nvar
            ! V_J * <Psi_J | w_J(r) | Psi_J>
            sum_b(1) = sum_b(1) + strength(ivar, jstate)*W_diagonal(ivar, jstate)
            ! V_I * <Psi_I | w_I(r) | Psi_I>
            sum_a(1) = sum_a(1) + strength(ivar, istate)*W_diagonal(ivar, istate)
            IF (mixed_cdft%identical_constraints) THEN
               ! V_J * W_IJ
               sum_b(2) = sum_b(2) + strength(ivar, jstate)*Wda(ivar, ipermutation)
               ! V_I * W_JI
               sum_a(2) = sum_a(2) + strength(ivar, istate)*Wda(ivar, ipermutation)
            ELSE
               ! V_J * W_IJ
               sum_b(2) = sum_b(2) + strength(ivar, jstate)*Wad(ivar, ipermutation)
               ! V_I * W_JI
               sum_a(2) = sum_a(2) + strength(ivar, istate)*Wda(ivar, ipermutation)
            END IF
         END DO
         ! Denote F_X = <Psi_X | H_X + V_X*w_X(r) | Psi_X> = E_X + V_X*<Psi_X | w_X(r) | Psi_X>
         ! H_IJ = F_J*S_IJ - V_J * W_IJ
         c(1) = (energy(jstate) + sum_b(1))*Sda(ipermutation) - sum_b(2)
         ! H_JI = F_I*S_JI - V_I * W_JI
         c(2) = (energy(istate) + sum_a(1))*Sda(ipermutation) - sum_a(2)
         ! H''(I,J) = 0.5*(H_IJ+H_JI) = H''(J,I)
         H_mat(istate, jstate) = (c(1) + c(2))*0.5_dp
         H_mat(jstate, istate) = H_mat(istate, jstate)
         IF (mixed_cdft%nonortho_coupling) coupling_nonortho(ipermutation) = H_mat(istate, jstate)
      END DO
      CALL mixed_cdft_result_type_set(mixed_cdft%results, H=H_mat)
      DEALLOCATE (H_mat, W_diagonal, Wda, strength, energy, Sda)
      IF (ALLOCATED(Wad)) DEALLOCATE (Wad)
      IF (mixed_cdft%nonortho_coupling) THEN
         CALL mixed_cdft_result_type_set(mixed_cdft%results, nonortho=coupling_nonortho)
         DEALLOCATE (coupling_nonortho)
      END IF
      ! Compute metric to assess reliability of coupling
      IF (mixed_cdft%calculate_metric) CALL mixed_cdft_calculate_metric(force_env, mixed_cdft, density_matrix_diff, ncol_mo)
      ! Compute coupling also with the wavefunction overlap method, see Migliore2009
      ! Requires the unconstrained KS ground state wavefunction as input
      IF (mixed_cdft%wfn_overlap_method) THEN
         IF (.NOT. uniform_occupation) &
            CALL cp_abort(__LOCATION__, &
                          "Wavefunction overlap method supports only uniformly occupied MOs.")
         CALL mixed_cdft_wfn_overlap_method(force_env, mixed_cdft, ncol_mo, nrow_mo)
      END IF
      ! Release remaining work
      DEALLOCATE (nrow_mo, ncol_mo)
      CALL mixed_cdft_work_type_release(mixed_cdft%matrix)
      CALL timestop(handle)

   END SUBROUTINE mixed_cdft_interaction_matrices

! **************************************************************************************************
!> \brief Routine to calculate the CDFT electronic couplings.
!> \param force_env the force_env that holds the CDFT states
!> \par History
!>       11.17  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_calculate_coupling_low(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_calculate_coupling_low', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ipermutation, istate, jstate, &
                                                            nforce_eval, npermutations, nvar
      LOGICAL                                            :: use_lowdin, use_rotation
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coupling_lowdin, coupling_rotation, &
                                                            eigenv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: tmp_mat, W_mat
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft

      NULLIFY (mixed_cdft)
      CPASSERT(ASSOCIATED(force_env))
      CALL get_mixed_env(force_env%mixed_env, cdft_control=mixed_cdft)
      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(mixed_cdft))
      CPASSERT(ALLOCATED(mixed_cdft%results%W_diagonal))
      CPASSERT(ALLOCATED(mixed_cdft%results%Wda))
      CPASSERT(ALLOCATED(mixed_cdft%results%S_minushalf))
      CPASSERT(ALLOCATED(mixed_cdft%results%H))
      ! Decide which methods to use for computing the coupling
      ! Default behavior is to use rotation when a single constraint is active, otherwise uses Lowdin orthogonalization
      ! The latter can also be explicitly requested when a single constraint is active
      ! Possibly computes the coupling additionally with the wavefunction overlap method
      nforce_eval = SIZE(mixed_cdft%results%H, 1)
      nvar = SIZE(mixed_cdft%results%Wda, 1)
      npermutations = nforce_eval*(nforce_eval - 1)/2
      ALLOCATE (tmp_mat(nforce_eval, nforce_eval))
      IF (nvar == 1 .AND. mixed_cdft%identical_constraints) THEN
         use_rotation = .TRUE.
         use_lowdin = mixed_cdft%use_lowdin
      ELSE
         use_rotation = .FALSE.
         use_lowdin = .TRUE.
      END IF
      ! Calculate coupling by rotating the CDFT states to eigenstates of the weight matrix W (single constraint only)
      IF (use_rotation) THEN
         ! Construct W
         ALLOCATE (W_mat(nforce_eval, nforce_eval), coupling_rotation(npermutations))
         ALLOCATE (eigenv(nforce_eval))
         ! W_mat(i, i) = N_i where N_i is the value of the constraint in state i
         DO istate = 1, nforce_eval
            W_mat(istate, istate) = SUM(mixed_cdft%results%W_diagonal(:, istate))
         END DO
         ! W_mat(i, j) = <Psi_i|w(r)|Psi_j>
         DO ipermutation = 1, npermutations
            CALL map_permutation_to_states(nforce_eval, ipermutation, istate, jstate)
            W_mat(istate, jstate) = SUM(mixed_cdft%results%Wda(:, ipermutation))
            W_mat(jstate, istate) = W_mat(istate, jstate)
         END DO
         ! Solve generalized eigenvalue equation WV = SVL
         ! Convert to standard eigenvalue problem via symmetric orthogonalisation
         tmp_mat(:, :) = MATMUL(W_mat, mixed_cdft%results%S_minushalf) ! W * S^(-1/2)
         W_mat(:, :) = MATMUL(mixed_cdft%results%S_minushalf, tmp_mat) ! W' = S^(-1/2) * W * S^(-1/2)
         CALL diamat_all(W_mat, eigenv, .TRUE.) ! Solve W'V' = AV'
         tmp_mat(:, :) = MATMUL(mixed_cdft%results%S_minushalf, W_mat) ! Reverse transformation V = S^(-1/2) V'
         ! Construct final, orthogonal diabatic Hamiltonian matrix H
         W_mat(:, :) = MATMUL(mixed_cdft%results%H, tmp_mat) ! H'' * V
         W_mat(:, :) = MATMUL(TRANSPOSE(tmp_mat), W_mat) ! H = V^T * H'' * V
         DO ipermutation = 1, npermutations
            CALL map_permutation_to_states(nforce_eval, ipermutation, istate, jstate)
            coupling_rotation(ipermutation) = W_mat(istate, jstate)
         END DO
         CALL mixed_cdft_result_type_set(mixed_cdft%results, rotation=coupling_rotation)
         DEALLOCATE (W_mat, coupling_rotation, eigenv)
      END IF
      ! Calculate coupling by Lowdin orthogonalization
      IF (use_lowdin) THEN
         ALLOCATE (coupling_lowdin(npermutations))
         tmp_mat(:, :) = MATMUL(mixed_cdft%results%H, mixed_cdft%results%S_minushalf) ! H'' * S^(-1/2)
         ! Final orthogonal diabatic Hamiltonian matrix H
         tmp_mat(:, :) = MATMUL(mixed_cdft%results%S_minushalf, tmp_mat) ! H = S^(-1/2) * H'' * S^(-1/2)
         DO ipermutation = 1, npermutations
            CALL map_permutation_to_states(nforce_eval, ipermutation, istate, jstate)
            coupling_lowdin(ipermutation) = tmp_mat(istate, jstate)
         END DO
         CALL mixed_cdft_result_type_set(mixed_cdft%results, lowdin=coupling_lowdin)
         DEALLOCATE (coupling_lowdin)
      END IF
      DEALLOCATE (tmp_mat)
      CALL timestop(handle)

   END SUBROUTINE mixed_cdft_calculate_coupling_low

! **************************************************************************************************
!> \brief Performs a configuration interaction calculation in the basis spanned by the CDFT states.
!> \param force_env the force_env that holds the CDFT states
!> \par History
!>       11.17  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_configuration_interaction(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_configuration_interaction', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, info, iounit, istate, ivar, &
                                                            nforce_eval, work_array_size
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenv, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: H_mat, H_mat_copy, S_mat, S_mat_copy
      REAL(KIND=dp), EXTERNAL                            :: dnrm2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(section_vals_type), POINTER                   :: force_env_section, print_section

      EXTERNAL                                           :: dsygv

      NULLIFY (logger, force_env_section, print_section, mixed_cdft)

      CPASSERT(ASSOCIATED(force_env))
      CALL get_mixed_env(force_env%mixed_env, cdft_control=mixed_cdft)
      CPASSERT(ASSOCIATED(mixed_cdft))

      IF (.NOT. mixed_cdft%do_ci) RETURN

      logger => cp_get_default_logger()
      CALL timeset(routineN, handle)
      CALL force_env_get(force_env=force_env, &
                         force_env_section=force_env_section)
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
      iounit = cp_print_key_unit_nr(logger, print_section, '', extension='.mixedLog')

      CPASSERT(ALLOCATED(mixed_cdft%results%S))
      CPASSERT(ALLOCATED(mixed_cdft%results%H))
      nforce_eval = SIZE(mixed_cdft%results%S, 1)
      ALLOCATE (S_mat(nforce_eval, nforce_eval), H_mat(nforce_eval, nforce_eval))
      ALLOCATE (eigenv(nforce_eval))
      S_mat(:, :) = mixed_cdft%results%S(:, :)
      H_mat(:, :) = mixed_cdft%results%H(:, :)
      ! Workspace query
      ALLOCATE (work(1))
      info = 0
      ALLOCATE (H_mat_copy(nforce_eval, nforce_eval), S_mat_copy(nforce_eval, nforce_eval))
      H_mat_copy(:, :) = H_mat(:, :) ! Need explicit copies because dsygv destroys original values
      S_mat_copy(:, :) = S_mat(:, :)
      CALL dsygv(1, 'V', 'U', nforce_eval, H_mat_copy, nforce_eval, S_mat_copy, nforce_eval, eigenv, work, -1, info)
      work_array_size = NINT(work(1))
      DEALLOCATE (H_mat_copy, S_mat_copy)
      ! Allocate work array
      DEALLOCATE (work)
      ALLOCATE (work(work_array_size))
      work = 0.0_dp
      ! Solve Hc = eSc
      info = 0
      CALL dsygv(1, 'V', 'U', nforce_eval, H_mat, nforce_eval, S_mat, nforce_eval, eigenv, work, work_array_size, info)
      IF (info /= 0) THEN
         IF (info > nforce_eval) THEN
            CPABORT("Matrix S is not positive definite")
         ELSE
            CPABORT("Diagonalization of H matrix failed.")
         END IF
      END IF
      ! dsygv returns eigenvectors (stored in columns of H_mat) that are normalized to H^T * S * H = I
      ! Renormalize eigenvectors to H^T * H = I
      DO ivar = 1, nforce_eval
         H_mat(:, ivar) = H_mat(:, ivar)/dnrm2(nforce_eval, H_mat(:, ivar), 1)
      END DO
      DEALLOCATE (work)
      IF (iounit > 0) THEN
         WRITE (iounit, '(/,T3,A)') '------------------ CDFT Configuration Interaction (CDFT-CI) ------------------'
         DO ivar = 1, nforce_eval
            IF (ivar == 1) THEN
               WRITE (iounit, '(T3,A,T58,(3X,F20.14))') 'Ground state energy:', eigenv(ivar)
            ELSE
               WRITE (iounit, '(/,T3,A,I2,A,T58,(3X,F20.14))') 'Excited state (', ivar - 1, ' ) energy:', eigenv(ivar)
            END IF
            DO istate = 1, nforce_eval, 2
               IF (istate == 1) THEN
                  WRITE (iounit, '(T3,A,T54,(3X,2F12.6))') &
                     'Expansion coefficients:', H_mat(istate, ivar), H_mat(istate + 1, ivar)
               ELSE IF (istate .LT. nforce_eval) THEN
                  WRITE (iounit, '(T54,(3X,2F12.6))') H_mat(istate, ivar), H_mat(istate + 1, ivar)
               ELSE
                  WRITE (iounit, '(T54,(3X,F12.6))') H_mat(istate, ivar)
               END IF
            END DO
         END DO
         WRITE (iounit, '(T3,A)') &
            '------------------------------------------------------------------------------'
      END IF
      DEALLOCATE (S_mat, H_mat, eigenv)
      CALL cp_print_key_finished_output(iounit, logger, force_env_section, &
                                        "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
      CALL timestop(handle)

   END SUBROUTINE mixed_cdft_configuration_interaction
! **************************************************************************************************
!> \brief Block diagonalizes the mixed CDFT Hamiltonian matrix.
!> \param force_env the force_env that holds the CDFT states
!> \par History
!>       11.17  created [Nico Holmberg]
!>       01.18  added recursive diagonalization
!>              split to subroutines [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_block_diag(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_block_diag', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iounit, irecursion, j, n, &
                                                            nblk, nforce_eval, nrecursion
      LOGICAL                                            :: ignore_excited
      TYPE(cp_1d_i_p_type), ALLOCATABLE, DIMENSION(:)    :: blocks
      TYPE(cp_1d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: eigenvalues
      TYPE(cp_2d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: H_block, S_block
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(section_vals_type), POINTER                   :: force_env_section, print_section

      EXTERNAL                                           :: dsygv

      NULLIFY (logger, force_env_section, print_section, mixed_cdft)

      CPASSERT(ASSOCIATED(force_env))
      CALL get_mixed_env(force_env%mixed_env, cdft_control=mixed_cdft)
      CPASSERT(ASSOCIATED(mixed_cdft))

      IF (.NOT. mixed_cdft%block_diagonalize) RETURN

      logger => cp_get_default_logger()
      CALL timeset(routineN, handle)

      CPASSERT(ALLOCATED(mixed_cdft%results%S))
      CPASSERT(ALLOCATED(mixed_cdft%results%H))
      nforce_eval = SIZE(mixed_cdft%results%S, 1)

      CALL force_env_get(force_env=force_env, &
                         force_env_section=force_env_section)
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
      iounit = cp_print_key_unit_nr(logger, print_section, '', extension='.mixedLog')
      ! Read block definitions from input
      CALL mixed_cdft_read_block_diag(force_env, blocks, ignore_excited, nrecursion)
      nblk = SIZE(blocks)
      ! Start block diagonalization
      DO irecursion = 1, nrecursion
         ! Print block definitions
         IF (iounit > 0 .AND. irecursion == 1) THEN
            WRITE (iounit, '(/,T3,A)') '-------------------------- CDFT BLOCK DIAGONALIZATION ------------------------'
            WRITE (iounit, '(T3,A)') 'Block diagonalizing the mixed CDFT Hamiltonian'
            WRITE (iounit, '(T3,A,I3)') 'Number of blocks:', nblk
            WRITE (iounit, '(T3,A,L3)') 'Ignoring excited states within blocks:', ignore_excited
            WRITE (iounit, '(/,T3,A)') 'List of CDFT states for each block'
            DO i = 1, nblk
               WRITE (iounit, '(T6,A,I3,A,6I3)') 'Block', i, ':', (blocks(i)%array(j), j=1, SIZE(blocks(i)%array))
            END DO
         END IF
         ! Recursive diagonalization: update counters and references
         IF (irecursion > 1) THEN
            nblk = nblk/2
            ALLOCATE (blocks(nblk))
            j = 1
            DO i = 1, nblk
               NULLIFY (blocks(i)%array)
               ALLOCATE (blocks(i)%array(2))
               blocks(i)%array = (/j, j + 1/)
               j = j + 2
            END DO
            ! Print info
            IF (iounit > 0) THEN
               WRITE (iounit, '(/, T3,A)') 'Recursive block diagonalization of the mixed CDFT Hamiltonian'
               WRITE (iounit, '(T6,A)') 'Block diagonalization is continued until only two matrix blocks remain.'
               WRITE (iounit, '(T6,A)') 'The new blocks are formed by collecting pairs of blocks from the previous'
               WRITE (iounit, '(T6,A)') 'block diagonalized matrix in ascending order.'
               WRITE (iounit, '(/,T3,A,I3,A,I3)') 'Recursion step:', irecursion - 1, ' of ', nrecursion - 1
               WRITE (iounit, '(/,T3,A)') 'List of old block indices for each new block'
               DO i = 1, nblk
                  WRITE (iounit, '(T6,A,I3,A,6I3)') 'Block', i, ':', (blocks(i)%array(j), j=1, SIZE(blocks(i)%array))
               END DO
            END IF
         END IF
         ! Get the Hamiltonian and overlap matrices of each block
         CALL mixed_cdft_get_blocks(mixed_cdft, blocks, H_block, S_block)
         ! Diagonalize blocks
         CALL mixed_cdft_diagonalize_blocks(blocks, H_block, S_block, eigenvalues)
         ! Assemble the block diagonalized matrices
         IF (ignore_excited) THEN
            n = nblk
         ELSE
            n = nforce_eval
         END IF
         CALL mixed_cdft_assemble_block_diag(mixed_cdft, blocks, H_block, eigenvalues, n, iounit)
         ! Deallocate work
         DO i = 1, nblk
            DEALLOCATE (H_block(i)%array)
            DEALLOCATE (S_block(i)%array)
            DEALLOCATE (eigenvalues(i)%array)
            DEALLOCATE (blocks(i)%array)
         END DO
         DEALLOCATE (H_block, S_block, eigenvalues, blocks)
      END DO ! recursion
      IF (iounit > 0) &
         WRITE (iounit, '(T3,A)') &
         '------------------------------------------------------------------------------'
      CALL cp_print_key_finished_output(iounit, logger, force_env_section, &
                                        "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
      CALL timestop(handle)

   END SUBROUTINE mixed_cdft_block_diag
! **************************************************************************************************
!> \brief Routine to calculate the CDFT electronic coupling reliability metric
!> \param force_env the force_env that holds the CDFT states
!> \param mixed_cdft the mixed_cdft env
!> \param density_matrix_diff array holding difference density matrices (P_j - P_i) for every CDFT
!>                            state permutation
!> \param ncol_mo the number of MOs per spin
!> \par History
!>       11.17  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_calculate_metric(force_env, mixed_cdft, density_matrix_diff, ncol_mo)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: density_matrix_diff
      INTEGER, DIMENSION(:)                              :: ncol_mo

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_calculate_metric', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ipermutation, ispin, j, &
                                                            nforce_eval, npermutations, nspins
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: metric
      TYPE(dbcsr_type)                                   :: e_vectors

      CALL timeset(routineN, handle)
      nforce_eval = SIZE(mixed_cdft%results%H, 1)
      npermutations = nforce_eval*(nforce_eval - 1)/2
      nspins = SIZE(density_matrix_diff, 2)
      ALLOCATE (metric(npermutations, nspins))
      metric = 0.0_dp
      CALL dbcsr_create(e_vectors, template=density_matrix_diff(1, 1)%matrix)
      DO ispin = 1, nspins
         ALLOCATE (evals(ncol_mo(ispin)))
         DO ipermutation = 1, npermutations
            ! Take into account doubly occupied orbitals without LSD
            IF (nspins == 1) &
               CALL dbcsr_scale(density_matrix_diff(ipermutation, 1)%matrix, alpha_scalar=0.5_dp)
            ! Diagonalize difference density matrix
            CALL cp_dbcsr_syevd(density_matrix_diff(ipermutation, ispin)%matrix, e_vectors, evals, &
                                para_env=force_env%para_env, blacs_env=mixed_cdft%blacs_env)
            CALL dbcsr_release_p(density_matrix_diff(ipermutation, ispin)%matrix)
            DO j = 1, ncol_mo(ispin)
               metric(ipermutation, ispin) = metric(ipermutation, ispin) + (evals(j)**2 - evals(j)**4)
            END DO
         END DO
         DEALLOCATE (evals)
      END DO
      CALL dbcsr_release(e_vectors)
      DEALLOCATE (density_matrix_diff)
      metric(:, :) = metric(:, :)/4.0_dp
      CALL mixed_cdft_result_type_set(mixed_cdft%results, metric=metric)
      DEALLOCATE (metric)
      CALL timestop(handle)

   END SUBROUTINE mixed_cdft_calculate_metric

! **************************************************************************************************
!> \brief Routine to calculate the electronic coupling according to the wavefunction overlap method
!> \param force_env the force_env that holds the CDFT states
!> \param mixed_cdft the mixed_cdft env
!> \param ncol_mo the number of MOs per spin
!> \param nrow_mo the number of AOs per spin
!> \par History
!>       11.17  created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_wfn_overlap_method(force_env, mixed_cdft, ncol_mo, nrow_mo)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(mixed_cdft_type), POINTER                     :: mixed_cdft
      INTEGER, DIMENSION(:)                              :: ncol_mo, nrow_mo

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_wfn_overlap_method', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: file_name
      INTEGER                                            :: handle, ipermutation, ispin, istate, &
                                                            jstate, nao, nforce_eval, nmo, &
                                                            npermutations, nspins
      LOGICAL                                            :: exist, natom_mismatch
      REAL(KIND=dp)                                      :: energy_diff, maxocc, Sda
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coupling_wfn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: overlaps
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: mixed_mo_coeff
      TYPE(cp_fm_struct_type), POINTER                   :: mo_mo_fmstruct
      TYPE(cp_fm_type), POINTER                          :: inverse_mat, mo_overlap_wfn, mo_tmp
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys_mix
      TYPE(dbcsr_type), POINTER                          :: mixed_matrix_s
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: force_env_section, mixed_cdft_section

      NULLIFY (mixed_cdft_section, subsys_mix, mo_set, particle_set, qs_kind_set, atomic_kind_set, &
               mixed_mo_coeff, mixed_matrix_s, force_env_section)
      logger => cp_get_default_logger()

      CALL timeset(routineN, handle)
      nforce_eval = SIZE(mixed_cdft%results%H, 1)
      npermutations = nforce_eval*(nforce_eval - 1)/2
      nspins = SIZE(nrow_mo)
      mixed_mo_coeff => mixed_cdft%matrix%mixed_mo_coeff
      mixed_matrix_s => mixed_cdft%matrix%mixed_matrix_s
      CALL force_env_get(force_env=force_env, &
                         force_env_section=force_env_section)
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
         ! Only OT with fully occupied orbitals is implicitly supported
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
      mixed_cdft_section => section_vals_get_subs_vals(force_env_section, "MIXED%MIXED_CDFT")
      CALL cp_subsys_get(subsys_mix, atomic_kind_set=atomic_kind_set, particle_set=particle_set)
      CPASSERT(ASSOCIATED(mixed_cdft%qs_kind_set))
      IF (force_env%para_env%ionode) &
         CALL wfn_restart_file_name(file_name, exist, mixed_cdft_section, logger)
      CALL mp_bcast(exist, force_env%para_env%source, force_env%para_env%group)
      CALL mp_bcast(file_name, force_env%para_env%source, force_env%para_env%group)
      IF (.NOT. exist) &
         CALL cp_abort(__LOCATION__, &
                       "User requested to restart the wavefunction from the file named: "// &
                       TRIM(file_name)//". This file does not exist. Please check the existence of"// &
                       " the file or change properly the value of the keyword WFN_RESTART_FILE_NAME in"// &
                       " section FORCE_EVAL\MIXED\MIXED_CDFT.")
      CALL read_mo_set(mo_array=mo_set, atomic_kind_set=atomic_kind_set, &
                       qs_kind_set=mixed_cdft%qs_kind_set, particle_set=particle_set, &
                       para_env=force_env%para_env, id_nr=0, multiplicity=mixed_cdft%multiplicity, &
                       dft_section=mixed_cdft_section, natom_mismatch=natom_mismatch, &
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
      ! Calculate MO overlaps between reference state (R) and CDFT state pairs I/J
      ALLOCATE (coupling_wfn(npermutations))
      ALLOCATE (overlaps(2, npermutations, nspins))
      overlaps = 0.0_dp
      DO ispin = 1, nspins
         ! Allocate work
         nao = nrow_mo(ispin)
         nmo = ncol_mo(ispin)
         CALL cp_fm_struct_create(mo_mo_fmstruct, nrow_global=nmo, ncol_global=nmo, &
                                  context=mixed_cdft%blacs_env, para_env=force_env%para_env)
         CALL cp_fm_create(matrix=mo_overlap_wfn, matrix_struct=mo_mo_fmstruct, &
                           name="MO_OVERLAP_MATRIX_WFN")
         CALL cp_fm_create(matrix=inverse_mat, matrix_struct=mo_mo_fmstruct, &
                           name="INVERSE_MO_OVERLAP_MATRIX_WFN")
         CALL cp_fm_struct_release(mo_mo_fmstruct)
         CALL cp_fm_create(matrix=mo_tmp, &
                           matrix_struct=mixed_mo_coeff(1, ispin)%matrix%matrix_struct, &
                           name="OVERLAP_MO_COEFF_WFN")
         DO ipermutation = 1, npermutations
            CALL map_permutation_to_states(nforce_eval, ipermutation, istate, jstate)
            ! S*C_r
            CALL cp_dbcsr_sm_fm_multiply(mixed_matrix_s, mo_set(ispin)%mo_set%mo_coeff, &
                                         mo_tmp, nmo, 1.0_dp, 0.0_dp)
            ! C_i^T * (S*C_r)
            CALL cp_fm_set_all(mo_overlap_wfn, alpha=0.0_dp)
            CALL cp_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, &
                         mixed_mo_coeff(istate, ispin)%matrix, &
                         mo_tmp, 0.0_dp, mo_overlap_wfn)
            CALL cp_fm_invert(mo_overlap_wfn, inverse_mat, overlaps(1, ipermutation, ispin), eps_svd=mixed_cdft%eps_svd)
            ! C_j^T * (S*C_r)
            CALL cp_fm_set_all(mo_overlap_wfn, alpha=0.0_dp)
            CALL cp_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, &
                         mixed_mo_coeff(jstate, ispin)%matrix, &
                         mo_tmp, 0.0_dp, mo_overlap_wfn)
            CALL cp_fm_invert(mo_overlap_wfn, inverse_mat, overlaps(2, ipermutation, ispin), eps_svd=mixed_cdft%eps_svd)
         END DO
         CALL cp_fm_release(mo_overlap_wfn)
         CALL cp_fm_release(inverse_mat)
         CALL cp_fm_release(mo_tmp)
         CALL deallocate_mo_set(mo_set(ispin)%mo_set)
      END DO
      DEALLOCATE (mo_set)
      DO ipermutation = 1, npermutations
         CALL map_permutation_to_states(nforce_eval, ipermutation, istate, jstate)
         IF (nspins == 2) THEN
            overlaps(1, ipermutation, 1) = ABS(overlaps(1, ipermutation, 1)*overlaps(1, ipermutation, 2)) ! A in eq. 12c
            overlaps(2, ipermutation, 1) = ABS(overlaps(2, ipermutation, 1)*overlaps(2, ipermutation, 2)) ! B in eq. 12c
         ELSE
            overlaps(1, ipermutation, 1) = overlaps(1, ipermutation, 1)**2
            overlaps(2, ipermutation, 1) = overlaps(2, ipermutation, 1)**2
         END IF
         ! Calculate coupling using eq. 12c
         ! The coupling is singular if A = B (i.e. states I/J are identical or charge in ground state is fully delocalized)
         IF (ABS(overlaps(1, ipermutation, 1) - overlaps(2, ipermutation, 1)) .LE. 1.0e-14_dp) THEN
            CALL cp_warn(__LOCATION__, &
                         "Coupling between states is singular and set to zero. "// &
                         "Potential causes: coupling is computed between identical CDFT states or the spin/charge "// &
                         "density is fully delocalized in the unconstrained ground state.")
            coupling_wfn(ipermutation) = 0.0_dp
         ELSE
            energy_diff = mixed_cdft%results%energy(jstate) - mixed_cdft%results%energy(istate)
            Sda = mixed_cdft%results%S(istate, jstate)
            coupling_wfn(ipermutation) = ABS((overlaps(1, ipermutation, 1)*overlaps(2, ipermutation, 1)/ &
                                              (overlaps(1, ipermutation, 1)**2 - overlaps(2, ipermutation, 1)**2))* &
                                             (energy_diff)/(1.0_dp - Sda**2)* &
                                             (1.0_dp - (overlaps(1, ipermutation, 1)**2 + overlaps(2, ipermutation, 1)**2)/ &
                                              (2.0_dp*overlaps(1, ipermutation, 1)*overlaps(2, ipermutation, 1))* &
                                              Sda))
         END IF
      END DO
      DEALLOCATE (overlaps)
      CALL mixed_cdft_result_type_set(mixed_cdft%results, wfn=coupling_wfn)
      DEALLOCATE (coupling_wfn)
      CALL timestop(handle)

   END SUBROUTINE mixed_cdft_wfn_overlap_method

! **************************************************************************************************
!> \brief Becke constraint adapted to mixed calculations, details in qs_cdft_methods.F
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
      INTEGER :: atom_a, bounds(2), handle, i, iatom, iex, iforce_eval, ikind, iounit, ithread, j, &
         jatom, katom, my_work, my_work_size, natom, nforce_eval, nkind, np(3), npme, nthread, &
         numexp, offset_dlb, unit_nr
      INTEGER, DIMENSION(2, 3)                           :: bo, bo_conf
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores, stride
      LOGICAL                                            :: build, mpi_io
      REAL(kind=dp)                                      :: alpha, chi, coef, ircov, jrcov, ra(3), &
                                                            uij
      REAL(kind=dp), DIMENSION(3)                        :: cell_v, dist_vec, dr, r, r1, shift
      REAL(KIND=dp), DIMENSION(:), POINTER               :: radii_list
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cdft_control_type), POINTER                   :: cdft_control
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
               atomic_kind_set, radii_list, cdft_control)
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
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
      cdft_control => mixed_cdft%cdft_control
      CPASSERT(ASSOCIATED(cdft_control))
      IF (.NOT. ASSOCIATED(cdft_control%becke_control%cutoffs)) THEN
         CALL cp_subsys_get(subsys_mix, atomic_kind_set=atomic_kind_set)
         ALLOCATE (cdft_control%becke_control%cutoffs(natom))
         SELECT CASE (cdft_control%becke_control%cutoff_type)
         CASE (becke_cutoff_global)
            cdft_control%becke_control%cutoffs(:) = cdft_control%becke_control%rglobal
         CASE (becke_cutoff_element)
            IF (.NOT. SIZE(atomic_kind_set) == SIZE(cdft_control%becke_control%cutoffs_tmp)) &
               CALL cp_abort(__LOCATION__, &
                             "Size of keyword BECKE_CONSTRAINT\ELEMENT_CUTOFFS does"// &
                             "not match number of atomic kinds in the input coordinate file.")
            DO ikind = 1, SIZE(atomic_kind_set)
               CALL get_atomic_kind(atomic_kind_set(ikind), natom=katom, atom_list=atom_list)
               DO iatom = 1, katom
                  atom_a = atom_list(iatom)
                  cdft_control%becke_control%cutoffs(atom_a) = cdft_control%becke_control%cutoffs_tmp(ikind)
               END DO
            END DO
            DEALLOCATE (cdft_control%becke_control%cutoffs_tmp)
         END SELECT
      END IF
      build = .FALSE.
      IF (cdft_control%becke_control%adjust .AND. &
          .NOT. ASSOCIATED(cdft_control%becke_control%aij)) THEN
         ALLOCATE (cdft_control%becke_control%aij(natom, natom))
         build = .TRUE.
      END IF
      ALLOCATE (catom(cdft_control%natoms))
      IF (cdft_control%save_pot .OR. &
          cdft_control%becke_control%cavity_confine .OR. &
          cdft_control%becke_control%should_skip .OR. &
          mixed_cdft%first_iteration) THEN
         ALLOCATE (is_constraint(natom))
         is_constraint = .FALSE.
      END IF
      in_memory = calculate_forces .AND. cdft_control%becke_control%in_memory
      IF (in_memory .NEQV. calculate_forces) &
         CALL cp_abort(__LOCATION__, &
                       "The flag BECKE_CONSTRAINT\IN_MEMORY must be activated "// &
                       "for the calculation of mixed CDFT forces")
      IF (in_memory .OR. mixed_cdft%first_iteration) ALLOCATE (coefficients(natom))
      DO i = 1, cdft_control%natoms
         catom(i) = cdft_control%atoms(i)
         IF (cdft_control%save_pot .OR. &
             cdft_control%becke_control%cavity_confine .OR. &
             cdft_control%becke_control%should_skip .OR. &
             mixed_cdft%first_iteration) &
            is_constraint(catom(i)) = .TRUE.
         IF (in_memory .OR. mixed_cdft%first_iteration) &
            coefficients(catom(i)) = cdft_control%group(1)%coeff(i)
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
      ALLOCATE (R12(natom, natom))
      DO iatom = 1, natom - 1
         DO jatom = iatom + 1, natom
            r = particle_set(iatom)%r
            r1 = particle_set(jatom)%r
            DO i = 1, 3
               r(i) = MODULO(r(i), cell%hmat(i, i)) - cell%hmat(i, i)/2._dp
               r1(i) = MODULO(r1(i), cell%hmat(i, i)) - cell%hmat(i, i)/2._dp
            END DO
            dist_vec = (r - r1) - ANINT((r - r1)/cell_v)*cell_v
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
            IF (build) THEN
               CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                    kind_number=ikind)
               ircov = cdft_control%becke_control%radii(ikind)
               CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, &
                                    kind_number=ikind)
               jrcov = cdft_control%becke_control%radii(ikind)
               IF (ircov .NE. jrcov) THEN
                  chi = ircov/jrcov
                  uij = (chi - 1.0_dp)/(chi + 1.0_dp)
                  cdft_control%becke_control%aij(iatom, jatom) = uij/(uij**2 - 1.0_dp)
                  IF (cdft_control%becke_control%aij(iatom, jatom) &
                      .GT. 0.5_dp) THEN
                     cdft_control%becke_control%aij(iatom, jatom) = 0.5_dp
                  ELSE IF (cdft_control%becke_control%aij(iatom, jatom) &
                           .LT. -0.5_dp) THEN
                     cdft_control%becke_control%aij(iatom, jatom) = -0.5_dp
                  END IF
               ELSE
                  cdft_control%becke_control%aij(iatom, jatom) = 0.0_dp
               END IF
               cdft_control%becke_control%aij(jatom, iatom) = &
                  -cdft_control%becke_control%aij(iatom, jatom)
            END IF
         END DO
      END DO
      ! Dump some additional information about the calculation
      IF (mixed_cdft%first_iteration) THEN
         print_section => section_vals_get_subs_vals(force_env_section, "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
         iounit = cp_print_key_unit_nr(logger, print_section, '', extension='.mixedLog')
         IF (iounit > 0) THEN
            WRITE (iounit, '(/,T3,A,T66)') &
               '-------------------------- Becke atomic parameters ---------------------------'
            IF (cdft_control%becke_control%adjust) THEN
               WRITE (iounit, '(T3,A,A)') &
                  'Atom   Element  Coefficient', '     Cutoff (angstrom)       CDFT Radius (angstrom)'
               DO iatom = 1, natom
                  CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                       element_symbol=element_symbol, &
                                       kind_number=ikind)
                  ircov = cp_unit_from_cp2k(cdft_control%becke_control%radii(ikind), "angstrom")
                  IF (is_constraint(iatom)) THEN
                     coef = coefficients(iatom)
                  ELSE
                     coef = 0.0_dp
                  END IF
                  WRITE (iounit, "(i6,T14,A2,T22,F8.3,T44,F8.3,T73,F8.3)") &
                     iatom, ADJUSTR(element_symbol), coef, &
                     cp_unit_from_cp2k(cdft_control%becke_control%cutoffs(iatom), "angstrom"), &
                     ircov
               END DO
            ELSE
               WRITE (iounit, '(T3,A,A)') &
                  'Atom   Element  Coefficient', '     Cutoff (angstrom)'
               DO iatom = 1, natom
                  CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                       element_symbol=element_symbol)
                  IF (is_constraint(iatom)) THEN
                     coef = coefficients(iatom)
                  ELSE
                     coef = 0.0_dp
                  END IF
                  WRITE (iounit, "(i6,T14,A2,T22,F8.3,T44,F8.3)") &
                     iatom, ADJUSTR(element_symbol), coef, &
                     cp_unit_from_cp2k(cdft_control%becke_control%cutoffs(iatom), "angstrom")
               END DO
            END IF
            WRITE (iounit, '(T3,A)') &
               '------------------------------------------------------------------------------'
         END IF
         CALL cp_print_key_finished_output(iounit, logger, force_env_section, &
                                           "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
         mixed_cdft%first_iteration = .FALSE.
      END IF

      IF (cdft_control%becke_control%cavity_confine) THEN
         CPASSERT(ASSOCIATED(mixed_cdft%qs_kind_set))
         cavity_env => cdft_control%becke_control%cavity_env
         qs_kind_set => mixed_cdft%qs_kind_set
         CALL cp_subsys_get(subsys_mix, atomic_kind_set=atomic_kind_set)
         nkind = SIZE(qs_kind_set)
         IF (.NOT. ASSOCIATED(cavity_env%kind_shape_fn)) THEN
            IF (ASSOCIATED(cdft_control%becke_control%radii)) THEN
               ALLOCATE (radii_list(SIZE(cdft_control%becke_control%radii)))
               DO ikind = 1, SIZE(cdft_control%becke_control%radii)
                  IF (cavity_env%use_bohr) THEN
                     radii_list(ikind) = cdft_control%becke_control%radii(ikind)
                  ELSE
                     radii_list(ikind) = cp_unit_from_cp2k(cdft_control%becke_control%radii(ikind), "angstrom")
                  END IF
               END DO
            END IF
            CALL create_shape_function(cavity_env, qs_kind_set, atomic_kind_set, &
                                       radius=cdft_control%becke_control%rcavity, &
                                       radii_list=radii_list)
            IF (ASSOCIATED(radii_list)) &
               DEALLOCATE (radii_list)
         END IF
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
                        npme = npme + 1
                        cores(npme) = iatom
                     ENDIF
                  ELSE
                     npme = npme + 1
                     cores(npme) = iatom
                  ENDIF
               END DO
               DO j = 1, npme
                  iatom = cores(j)
                  atom_a = atom_list(iatom)
                  pab(1, 1) = coef
                  IF (store_vectors) THEN
                     ra(:) = position_vecs(:, atom_a) + cell_v(:)/2._dp
                  ELSE
                     ra(:) = pbc(particle_set(atom_a)%r, cell)
                  END IF
                  IF (is_constraint(atom_a)) &
                     CALL collocate_pgf_product(0, alpha, 0, 0, 0.0_dp, 0, ra, &
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
         CALL pw_pool_create_pw(auxbas_pw_pool, cdft_control%becke_control%cavity%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL rs_pw_transfer(rs_cavity, cdft_control%becke_control%cavity%pw, rs2pw)
         CALL rs_grid_release(rs_cavity)
         CALL hfun_zero(cdft_control%becke_control%cavity%pw%cr3d, &
                        cdft_control%becke_control%eps_cavity, &
                        just_zero=.FALSE., bounds=bounds, work=my_work)
         IF (bounds(2) .LT. bo(2, 3)) THEN
            bounds(2) = bounds(2) - 1
         ELSE
            bounds(2) = bo(2, 3)
         END IF
         IF (bounds(1) .GT. bo(1, 3)) THEN
            ! In the special case bounds(1) == bounds(2) == bo(2, 3), after this check
            ! bounds(1) > bounds(2) and the subsequent gradient allocation (:, :, :, bounds(1):bounds(2))
            ! will correctly allocate a 0-sized array
            bounds(1) = bounds(1) + 1
         ELSE
            bounds(1) = bo(1, 3)
         END IF
         IF (bounds(1) > bounds(2)) THEN
            my_work_size = 0
         ELSE
            my_work_size = (bounds(2) - bounds(1) + 1)
            IF (mixed_cdft%is_pencil .OR. mixed_cdft%is_special) THEN
               my_work_size = my_work_size*(bo(2, 2) - bo(1, 2) + 1)
            ELSE
               my_work_size = my_work_size*(bo(2, 1) - bo(1, 1) + 1)
            END IF
         END IF
         cdft_control%becke_control%confine_bounds = bounds
         IF (cdft_control%becke_control%print_cavity) THEN
            CALL hfun_zero(cdft_control%becke_control%cavity%pw%cr3d, &
                           cdft_control%becke_control%eps_cavity, just_zero=.TRUE.)
            NULLIFY (stride)
            ALLOCATE (stride(3))
            stride = (/2, 2, 2/)
            mpi_io = .TRUE.
            unit_nr = cp_print_key_unit_nr(logger, print_section, "", &
                                           middle_name="BECKE_CAVITY", &
                                           extension=".cube", file_position="REWIND", &
                                           log_filename=.FALSE., mpi_io=mpi_io)
            IF (force_env%para_env%ionode .AND. unit_nr .LT. 1) &
               CALL cp_abort(__LOCATION__, &
                             "Please turn on PROGRAM_RUN_INFO to print cavity")
            CALL cp_pw_to_cube(cdft_control%becke_control%cavity%pw, &
                               unit_nr, "CAVITY", particles=particles, &
                               stride=stride, mpi_io=mpi_io)
            CALL cp_print_key_finished_output(unit_nr, logger, print_section, '', mpi_io=mpi_io)
            DEALLOCATE (stride)
         END IF
      END IF
      bo_conf = bo
      IF (cdft_control%becke_control%cavity_confine) &
         bo_conf(:, 3) = cdft_control%becke_control%confine_bounds
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
      IF (cdft_control%becke_control%cavity_confine) THEN
         ! Get rid of the zero part of the confinement cavity (cr3d -> real(:,:,:))
         IF (mixed_cdft%is_special) THEN
            ALLOCATE (mixed_cdft%sendbuff(SIZE(mixed_cdft%dest_list)))
            DO i = 1, SIZE(mixed_cdft%dest_list)
               ALLOCATE (mixed_cdft%sendbuff(i)%cavity(mixed_cdft%dest_list_bo(1, i):mixed_cdft%dest_list_bo(2, i), &
                                                       bo(1, 2):bo(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
               mixed_cdft%sendbuff(i)%cavity = cdft_control%becke_control%cavity%pw%cr3d(mixed_cdft%dest_list_bo(1, i): &
                                                                                         mixed_cdft%dest_list_bo(2, i), &
                                                                                         bo(1, 2):bo(2, 2), &
                                                                                         bo_conf(1, 3):bo_conf(2, 3))
            END DO
         ELSE IF (mixed_cdft%is_pencil) THEN
            ALLOCATE (mixed_cdft%cavity(bo(1, 1) + offset_dlb:bo(2, 1), bo(1, 2):bo(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
            mixed_cdft%cavity = cdft_control%becke_control%cavity%pw%cr3d(bo(1, 1) + offset_dlb:bo(2, 1), &
                                                                          bo(1, 2):bo(2, 2), &
                                                                          bo_conf(1, 3):bo_conf(2, 3))
         ELSE
            ALLOCATE (mixed_cdft%cavity(bo(1, 1):bo(2, 1), bo(1, 2) + offset_dlb:bo(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
            mixed_cdft%cavity = cdft_control%becke_control%cavity%pw%cr3d(bo(1, 1):bo(2, 1), &
                                                                          bo(1, 2) + offset_dlb:bo(2, 2), &
                                                                          bo_conf(1, 3):bo_conf(2, 3))
         END IF
         CALL pw_pool_give_back_pw(auxbas_pw_pool, cdft_control%becke_control%cavity%pw)
      END IF
      IF (mixed_cdft%is_special) THEN
         DO i = 1, SIZE(mixed_cdft%dest_list)
            ALLOCATE (mixed_cdft%sendbuff(i)%weight(mixed_cdft%dest_list_bo(1, i):mixed_cdft%dest_list_bo(2, i), &
                                                    bo(1, 2):bo(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
            mixed_cdft%sendbuff(i)%weight = 0.0_dp
         END DO
      ELSE IF (mixed_cdft%is_pencil) THEN
         ALLOCATE (mixed_cdft%weight(bo(1, 1) + offset_dlb:bo(2, 1), bo(1, 2):bo(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
         mixed_cdft%weight = 0.0_dp
      ELSE
         ALLOCATE (mixed_cdft%weight(bo(1, 1):bo(2, 1), bo(1, 2) + offset_dlb:bo(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
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
            ALLOCATE (cdft_control%group(1)%gradients(3*natom, bo(1, 1) + offset_dlb:bo(2, 1), &
                                                      bo(1, 2):bo(2, 2), &
                                                      bo_conf(1, 3):bo_conf(2, 3)))
            cdft_control%group(1)%gradients = 0.0_dp
         ELSE
            ALLOCATE (cdft_control%group(1)%gradients(3*natom, bo(1, 1):bo(2, 1), &
                                                      bo(1, 2) + offset_dlb:bo(2, 2), &
                                                      bo_conf(1, 3):bo_conf(2, 3)))
            cdft_control%group(1)%gradients = 0.0_dp
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
         more_work, my_pos, my_special_work, my_target, no_overloaded, no_underloaded, nsend, &
         nsend_limit, nsend_max, offset, offset_proc, offset_special, req(4), send_total, tags(2)
      INTEGER, DIMENSION(:), POINTER :: buffsize, cumulative_work, expected_work, load_imbalance, &
         nrecv, nsend_proc, req_recv, req_total, sendbuffer, should_warn, tmp, work_index, &
         work_size
      INTEGER, DIMENSION(:, :), POINTER                  :: targets, tmp_bo
      LOGICAL                                            :: consistent
      LOGICAL, DIMENSION(:), POINTER                     :: mask_recv, mask_send, touched
      REAL(kind=dp)                                      :: average_work, load_scale, &
                                                            very_overloaded, work_factor
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: cavity
      TYPE(cdft_control_type), POINTER                   :: cdft_control
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
               print_section, cdft_control)
      CALL force_env_get(force_env=force_env, force_env_section=force_env_section)
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
      iounit = cp_print_key_unit_nr(logger, print_section, '', extension='.mixedLog')
      cdft_control => mixed_cdft%cdft_control
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
      expected_work(force_env%para_env%mepos + 1) = my_work
      work_size = 0
      work_size(force_env%para_env%mepos + 1) = my_work_size
      IF (ASSOCIATED(mixed_cdft%dlb_control%prediction_error)) THEN
         IF (mixed_cdft%is_pencil .OR. mixed_cdft%is_special) THEN
            work_size(force_env%para_env%mepos + 1) = work_size(force_env%para_env%mepos + 1) - &
                                                      NINT(REAL(mixed_cdft%dlb_control% &
                                                                prediction_error(force_env%para_env%mepos + 1), dp)/ &
                                                           REAL(bo(2, 1) - bo(1, 1) + 1, dp))
         ELSE
            work_size(force_env%para_env%mepos + 1) = work_size(force_env%para_env%mepos + 1) - &
                                                      NINT(REAL(mixed_cdft%dlb_control% &
                                                                prediction_error(force_env%para_env%mepos + 1), dp)/ &
                                                           REAL(bo(2, 2) - bo(1, 2) + 1, dp))
         END IF
      END IF
      CALL mp_sum(expected_work, force_env%para_env%group)
      CALL mp_sum(work_size, force_env%para_env%group)
      ! We store the unsorted expected work to refine the estimate on subsequent calls to this routine
      mixed_cdft%dlb_control%expected_work = expected_work
      ! Take into account the prediction error of the last step
      IF (ASSOCIATED(mixed_cdft%dlb_control%prediction_error)) &
         expected_work = expected_work - mixed_cdft%dlb_control%prediction_error
      !
      average_work = REAL(SUM(expected_work), dp)/REAL(force_env%para_env%num_pe, dp)
      ALLOCATE (work_index(force_env%para_env%num_pe), &
                load_imbalance(force_env%para_env%num_pe), &
                targets(2, force_env%para_env%num_pe))
      load_imbalance = expected_work - NINT(average_work)
      no_overloaded = 0
      no_underloaded = 0
      targets = 0
      ! Convert the load imbalance to a multiple of the actual work size
      DO i = 1, force_env%para_env%num_pe
         IF (load_imbalance(i) .GT. 0) THEN
            no_overloaded = no_overloaded + 1
            ! Allow heavily overloaded processors to dump more data since most likely they have a lot of 'real' work
            IF (expected_work(i) .GT. NINT(very_overloaded*average_work)) THEN
               load_imbalance(i) = (CEILING(REAL(load_imbalance(i), dp)/REAL(work_size(i), dp)) + more_work)*work_size(i)
            ELSE
               load_imbalance(i) = CEILING(REAL(load_imbalance(i), dp)/REAL(work_size(i), dp))*work_size(i)
            END IF
         ELSE
            ! Allow the underloaded processors to take load_scale amount of additional work
            ! otherwise we may be unable to exhaust all overloaded processors
            load_imbalance(i) = NINT(load_imbalance(i)*load_scale)
            no_underloaded = no_underloaded + 1
         END IF
      END DO
      CALL sort(expected_work, force_env%para_env%num_pe, indices=work_index)
      ! Redistribute work in order from the most overloaded processors to the most underloaded processors
      ! Each underloaded processor is limited to one overloaded processor
      IF (load_imbalance(force_env%para_env%mepos + 1) > 0) THEN
         offset = 0
         mixed_cdft%dlb_control%send_work = .TRUE.
         ! Build up the total amount of work that needs redistribution
         ALLOCATE (cumulative_work(force_env%para_env%num_pe))
         cumulative_work = 0
         DO i = force_env%para_env%num_pe, force_env%para_env%num_pe - no_overloaded + 1, -1
            IF (work_index(i) == force_env%para_env%mepos + 1) THEN
               EXIT
            ELSE
               offset = offset + load_imbalance(work_index(i))
               IF (i == force_env%para_env%num_pe) THEN
                  cumulative_work(i) = load_imbalance(work_index(i))
               ELSE
                  cumulative_work(i) = cumulative_work(i + 1) + load_imbalance(work_index(i))
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
            nsend_max = nsend_max - nsend
            IF (nsend_max .LT. 0) nsend = nsend + nsend_max
            exhausted_work = exhausted_work + nsend*work_size(work_index(j))
            offset = offset - nsend*work_size(work_index(j))
            IF (offset .LT. 0) EXIT
            IF (exhausted_work .EQ. cumulative_work(j)) THEN
               j = j - 1
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
         nsend_max = load_imbalance(force_env%para_env%mepos + 1)/work_size(force_env%para_env%mepos + 1)
         ! This the actual number of available array slices
         IF (mixed_cdft%is_pencil .OR. mixed_cdft%is_special) THEN
            nsend_limit = bo(2, 1) - bo(1, 1) + 1
         ELSE
            nsend_limit = bo(2, 2) - bo(1, 2) + 1
         END IF
         IF (.NOT. mixed_cdft%is_special) THEN
            ALLOCATE (mixed_cdft%dlb_control%target_list(3, max_targets))
         ELSE
            ALLOCATE (mixed_cdft%dlb_control%target_list(3 + 2*SIZE(mixed_cdft%dest_list), max_targets))
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
            nsend = -load_imbalance(work_index(my_target))/work_size(force_env%para_env%mepos + 1)
            IF (nsend .LT. 1) nsend = 1 ! send at least one block
            ! Prevent over redistribution: leave at least (1-work_factor)*nsend_limit slices to my_pos
            IF (nsend .GT. NINT(work_factor*nsend_limit - send_total)) THEN
               nsend = NINT(work_factor*nsend_limit - send_total)
               IF (debug_this_module) &
                  should_warn(force_env%para_env%mepos + 1) = 1
            END IF
            mixed_cdft%dlb_control%target_list(1, i) = work_index(my_target) - 1 ! This is the actual processor rank
            IF (mixed_cdft%is_special) THEN
               mixed_cdft%dlb_control%target_list(2, i) = 0
               actually_sent = nsend
               DO j = ispecial, SIZE(mixed_cdft%dest_list)
                  mixed_cdft%dlb_control%target_list(2, i) = mixed_cdft%dlb_control%target_list(2, i) + 1
                  touched(j) = .TRUE.
                  IF (nsend .LT. mixed_cdft%dest_list_bo(2, j) - mixed_cdft%dest_list_bo(1, j) + 1) THEN
                     mixed_cdft%dlb_control%target_list(3 + 2*j - 1, i) = mixed_cdft%dest_list_bo(1, j)
                     mixed_cdft%dlb_control%target_list(3 + 2*j, i) = mixed_cdft%dest_list_bo(1, j) + nsend - 1
                     mixed_cdft%dest_list_bo(1, j) = mixed_cdft%dest_list_bo(1, j) + nsend
                     nsend = 0
                     EXIT
                  ELSE
                     mixed_cdft%dlb_control%target_list(3 + 2*j - 1, i) = mixed_cdft%dest_list_bo(1, j)
                     mixed_cdft%dlb_control%target_list(3 + 2*j, i) = mixed_cdft%dest_list_bo(2, j)
                     nsend = nsend - (mixed_cdft%dest_list_bo(2, j) - mixed_cdft%dest_list_bo(1, j) + 1)
                     mixed_cdft%dest_list_bo(1:2, j) = should_deallocate
                  END IF
                  IF (nsend .LE. 0) EXIT
               END DO
               IF (mixed_cdft%dest_list_bo(1, ispecial) .EQ. should_deallocate) ispecial = j + 1
               actually_sent = actually_sent - nsend
               nsend_max = nsend_max - actually_sent
               send_total = send_total + actually_sent
            ELSE
               mixed_cdft%dlb_control%target_list(2, i) = nsend
               nsend_max = nsend_max - nsend
               send_total = send_total + nsend
            END IF
            IF (nsend_max .LT. 0) nsend_max = 0
            IF (nsend_max .EQ. 0) EXIT
            IF (my_target /= no_underloaded) THEN
               my_target = my_target + 1
            ELSE
               ! If multiple processors execute this block load balancing will fail
               mixed_cdft%dlb_control%target_list(2, i) = mixed_cdft%dlb_control%target_list(2, i) + nsend_max
               nsend_max = 0
               EXIT
            END IF
            i = i + 1
            IF (i .GT. max_targets) &
               CALL cp_abort(__LOCATION__, &
                             "Load balancing error: increase max_targets")
         END DO
         IF (.NOT. mixed_cdft%is_special) THEN
            CALL reallocate(mixed_cdft%dlb_control%target_list, 1, 3, 1, i)
         ELSE
            CALL reallocate(mixed_cdft%dlb_control%target_list, 1, 3 + 2*SIZE(mixed_cdft%dest_list), 1, i)
         END IF
         targets(2, my_pos) = my_target
         ! Equalize the load on the target processors
         IF (.NOT. mixed_cdft%is_special) THEN
            IF (send_total .GT. NINT(work_factor*nsend_limit)) send_total = NINT(work_factor*nsend_limit) - 1
            nsend = NINT(REAL(send_total, dp)/REAL(SIZE(mixed_cdft%dlb_control%target_list, 2), dp))
            mixed_cdft%dlb_control%target_list(2, :) = nsend
         END IF
      ELSE
         DO i = 1, no_underloaded
            IF (work_index(i) == force_env%para_env%mepos + 1) EXIT
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
      IF (force_env%para_env%ionode) THEN
         consistent = .TRUE.
         DO i = force_env%para_env%num_pe - 1, force_env%para_env%num_pe - no_overloaded + 1, -1
            IF (targets(1, i) .GT. no_underloaded) consistent = .FALSE.
            IF (targets(1, i) .GT. targets(2, i + 1)) THEN
               CYCLE
            ELSE
               consistent = .FALSE.
            END IF
         END DO
         IF (.NOT. consistent) THEN
            IF (debug_this_module .AND. iounit > 0) THEN
               DO i = force_env%para_env%num_pe - 1, force_env%para_env%num_pe - no_overloaded + 1, -1
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
         DO i = force_env%para_env%num_pe, force_env%para_env%num_pe - no_overloaded + 1, -1
            IF (targets(1, i) .LE. my_pos .AND. targets(2, i) .GE. my_pos) THEN
               mixed_cdft%dlb_control%recv_work = .TRUE.
               mixed_cdft%dlb_control%my_source = work_index(i) - 1
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
                  ALLOCATE (mixed_cdft%dlb_control%sendbuff(j)%cavity(mixed_cdft%dlb_control%bo(12*(j - 1) + 1): &
                                                                      mixed_cdft%dlb_control%bo(12*(j - 1) + 2), &
                                                                      mixed_cdft%dlb_control%bo(12*(j - 1) + 3): &
                                                                      mixed_cdft%dlb_control%bo(12*(j - 1) + 4), &
                                                                      mixed_cdft%dlb_control%bo(12*(j - 1) + 7): &
                                                                      mixed_cdft%dlb_control%bo(12*(j - 1) + 8)))
                  CALL mp_irecv(msgout=mixed_cdft%dlb_control%sendbuff(j)%cavity, &
                                source=mixed_cdft%dlb_control%my_source, &
                                request=req_recv(j), comm=force_env%para_env%group)
                  ALLOCATE (mixed_cdft%dlb_control%sendbuff(j)%weight(mixed_cdft%dlb_control%bo(12*(j - 1) + 1): &
                                                                      mixed_cdft%dlb_control%bo(12*(j - 1) + 2), &
                                                                      mixed_cdft%dlb_control%bo(12*(j - 1) + 3): &
                                                                      mixed_cdft%dlb_control%bo(12*(j - 1) + 4), &
                                                                      mixed_cdft%dlb_control%bo(12*(j - 1) + 7): &
                                                                      mixed_cdft%dlb_control%bo(12*(j - 1) + 8)))
                  ALLOCATE (mixed_cdft%dlb_control%sendbuff(j)%gradients(3*natom, &
                                                                         mixed_cdft%dlb_control%bo(12*(j - 1) + 1): &
                                                                         mixed_cdft%dlb_control%bo(12*(j - 1) + 2), &
                                                                         mixed_cdft%dlb_control%bo(12*(j - 1) + 3): &
                                                                         mixed_cdft%dlb_control%bo(12*(j - 1) + 4), &
                                                                         mixed_cdft%dlb_control%bo(12*(j - 1) + 7): &
                                                                         mixed_cdft%dlb_control%bo(12*(j - 1) + 8)))
                  mixed_cdft%dlb_control%sendbuff(j)%weight = 0.0_dp
                  mixed_cdft%dlb_control%sendbuff(j)%gradients = 0.0_dp
                  mixed_cdft%dlb_control%sendbuff(j)%tag = (/mixed_cdft%dlb_control%bo(12*(j - 1) + 9), &
                                                             mixed_cdft%dlb_control%bo(12*(j - 1) + 10)/)
                  mixed_cdft%dlb_control%sendbuff(j)%rank = (/mixed_cdft%dlb_control%bo(12*(j - 1) + 11), &
                                                              mixed_cdft%dlb_control%bo(12*(j - 1) + 12)/)
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
               tags = (/(i - 1)*3 + 1 + force_env%para_env%mepos*6*max_targets, &
                        (i - 1)*3 + 1 + 3*max_targets + force_env%para_env%mepos*6*max_targets/) ! Unique communicator tags
               mixed_cdft%dlb_control%target_list(3, i) = tags(1)
               IF (mixed_cdft%is_pencil) THEN
                  sendbuffer = (/bo_conf(1, 1) + offset, &
                                 bo_conf(1, 1) + offset + (mixed_cdft%dlb_control%target_list(2, i) - 1), &
                                 bo_conf(1, 2), bo_conf(2, 2), bo(1, 3), bo(2, 3), bo_conf(1, 3), bo_conf(2, 3), &
                                 tags(1), tags(2), mixed_cdft%dest_list(1), mixed_cdft%dest_list(2)/)
               ELSE
                  sendbuffer = (/bo_conf(1, 1), bo_conf(2, 1), &
                                 bo_conf(1, 2) + offset, &
                                 bo_conf(1, 2) + offset + (mixed_cdft%dlb_control%target_list(2, i) - 1), &
                                 bo(1, 3), bo(2, 3), bo_conf(1, 3), bo_conf(2, 3), tags(1), tags(2), &
                                 mixed_cdft%dest_list(1), mixed_cdft%dest_list(2)/)
               END IF
               send_total = send_total + mixed_cdft%dlb_control%target_list(2, i) - 1
               CALL mp_isend(msgin=sendbuffer, dest=mixed_cdft%dlb_control%target_list(1, i), &
                             request=req(1), comm=force_env%para_env%group)
               CALL mp_wait(req(1))
               IF (mixed_cdft%is_pencil) THEN
                  ALLOCATE (cavity(bo_conf(1, 1) + offset: &
                                   bo_conf(1, 1) + offset + (mixed_cdft%dlb_control%target_list(2, i) - 1), &
                                   bo_conf(1, 2):bo_conf(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
                  cavity = cdft_control%becke_control%cavity%pw%cr3d(bo_conf(1, 1) + offset: &
                                                                     bo_conf(1, 1) + offset + &
                                                                     (mixed_cdft%dlb_control%target_list(2, i) - 1), &
                                                                     bo_conf(1, 2):bo_conf(2, 2), &
                                                                     bo_conf(1, 3):bo_conf(2, 3))
               ELSE
                  ALLOCATE (cavity(bo_conf(1, 1):bo_conf(2, 1), &
                                   bo_conf(1, 2) + offset: &
                                   bo_conf(1, 2) + offset + (mixed_cdft%dlb_control%target_list(2, i) - 1), &
                                   bo_conf(1, 3):bo_conf(2, 3)))
                  cavity = cdft_control%becke_control%cavity%pw%cr3d(bo_conf(1, 1):bo_conf(2, 1), &
                                                                     bo_conf(1, 2) + offset: &
                                                                     bo_conf(1, 2) + offset + &
                                                                     (mixed_cdft%dlb_control%target_list(2, i) - 1), &
                                                                     bo_conf(1, 3):bo_conf(2, 3))
               END IF
               CALL mp_isend(msgin=cavity, &
                             dest=mixed_cdft%dlb_control%target_list(1, i), &
                             request=req(1), comm=force_env%para_env%group)
               CALL mp_wait(req(1))
               offset = offset + mixed_cdft%dlb_control%target_list(2, i)
               DEALLOCATE (cavity)
            END DO
            IF (mixed_cdft%is_pencil) THEN
               mixed_cdft%dlb_control%distributed(1) = bo_conf(1, 1)
               mixed_cdft%dlb_control%distributed(2) = bo_conf(1, 1) + offset - 1
            ELSE
               mixed_cdft%dlb_control%distributed(1) = bo_conf(1, 2)
               mixed_cdft%dlb_control%distributed(2) = bo_conf(1, 2) + offset - 1
            END IF
            DEALLOCATE (sendbuffer)
         ELSE
            ALLOCATE (buffsize(1))
            DO i = 1, SIZE(mixed_cdft%dlb_control%target_list, 2)
               buffsize = mixed_cdft%dlb_control%target_list(2, i)
               ! Unique communicator tags (dont actually need these, should be removed)
               tags = (/(i - 1)*3 + 1 + force_env%para_env%mepos*6*max_targets, &
                        (i - 1)*3 + 1 + 3*max_targets + force_env%para_env%mepos*6*max_targets/)
               DO j = 4, SIZE(mixed_cdft%dlb_control%target_list, 1)
                  IF (mixed_cdft%dlb_control%target_list(j, i) .GT. uninitialized) EXIT
               END DO
               offset_special = j
               offset_proc = j - 4 - (j - 4)/2
               CALL mp_isend(msgin=buffsize, &
                             dest=mixed_cdft%dlb_control%target_list(1, i), &
                             request=req(1), comm=force_env%para_env%group)
               CALL mp_wait(req(1))
               ALLOCATE (sendbuffer(12*buffsize(1)))
               DO j = 1, buffsize(1)
                 sendbuffer(12*(j - 1) + 1:12*(j - 1) + 12) = (/mixed_cdft%dlb_control%target_list(offset_special + 2*(j - 1), i), &
                                                                 mixed_cdft%dlb_control%target_list(offset_special + 2*j - 1, i), &
                                                                 bo_conf(1, 2), bo_conf(2, 2), bo(1, 3), bo(2, 3), &
                                                                 bo_conf(1, 3), bo_conf(2, 3), tags(1), tags(2), &
                                                                 mixed_cdft%dest_list(j + offset_proc), &
                                                               mixed_cdft%dest_list(j + offset_proc) + force_env%para_env%num_pe/2/)
               END DO
               CALL mp_isend(msgin=sendbuffer, &
                             dest=mixed_cdft%dlb_control%target_list(1, i), &
                             request=req(1), comm=force_env%para_env%group)
               CALL mp_wait(req(1))
               DEALLOCATE (sendbuffer)
               DO j = 1, buffsize(1)
                  ALLOCATE (cavity(mixed_cdft%dlb_control%target_list(offset_special + 2*(j - 1), i): &
                                   mixed_cdft%dlb_control%target_list(offset_special + 2*j - 1, i), &
                                   bo_conf(1, 2):bo_conf(2, 2), bo_conf(1, 3):bo_conf(2, 3)))
                  cavity = cdft_control%becke_control%cavity%pw%cr3d(LBOUND(cavity, 1):UBOUND(cavity, 1), &
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
      ALLOCATE (req_total(my_special_work*SIZE(mixed_cdft%source_list) + (my_special_work**2)*SIZE(mixed_cdft%dest_list)))
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
                          request=req_total(i + SIZE(mixed_cdft%source_list)), &
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
                           IF (mixed_cdft%dlb_control%target_list(4 + 2*(j - 1), ispecial) .NE. uninitialized) &
                              nsend = nsend + 1
                        END DO
                        sbuff(j)%iv(3) = nsend
                        nsend_proc(j) = nsend
                     END IF
                  END IF
               END IF
            END IF
            ind = j + (i - 1)*SIZE(mixed_cdft%dest_list) + my_special_work*SIZE(mixed_cdft%source_list)
            CALL mp_isend(msgin=sbuff(j)%bv, &
                          dest=mixed_cdft%dest_list(j) + (i - 1)*force_env%para_env%num_pe/2, &
                          request=req_total(ind), tag=1, &
                          comm=force_env%para_env%group)
            IF (mixed_cdft%is_special) &
               CALL mp_isend(msgin=sbuff(j)%iv, &
                             dest=mixed_cdft%dest_list(j) + (i - 1)*force_env%para_env%num_pe/2, &
                             request=req_total(ind + 2*SIZE(mixed_cdft%dest_list)), tag=2, &
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
            ALLOCATE (req_total(COUNT(mixed_cdft%dlb_control%recv_work_repl) + 2))
            ALLOCATE (sendbuffer(6))
            IF (mixed_cdft%is_pencil) THEN
               sendbuffer = (/SIZE(mixed_cdft%dlb_control%target_list, 2), bo_conf(1, 3), bo_conf(2, 3), &
                              bo_conf(1, 1), bo_conf(1, 2), bo_conf(2, 2)/)
            ELSE
               sendbuffer = (/SIZE(mixed_cdft%dlb_control%target_list, 2), bo_conf(1, 3), bo_conf(2, 3), &
                              bo_conf(1, 2), bo_conf(1, 1), bo_conf(2, 1)/)
            END IF
         ELSE IF (ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
            ALLOCATE (req_total(COUNT(mixed_cdft%dlb_control%recv_work_repl)))
         END IF
         IF (ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
            ALLOCATE (mixed_cdft%dlb_control%recv_info(2))
            NULLIFY (mixed_cdft%dlb_control%recv_info(1)%target_list, mixed_cdft%dlb_control%recv_info(2)%target_list)
            ALLOCATE (mixed_cdft%dlb_control%recvbuff(2))
            NULLIFY (mixed_cdft%dlb_control%recvbuff(1)%buffs, mixed_cdft%dlb_control%recvbuff(2)%buffs)
         END IF
         ! First communicate which grid points were distributed
         IF (mixed_cdft%dlb_control%send_work) THEN
            ind = COUNT(mixed_cdft%dlb_control%recv_work_repl) + 1
            DO i = 1, 2
               CALL mp_isend(msgin=sendbuffer, &
                             dest=mixed_cdft%dest_list(i), &
                             request=req_total(ind), comm=force_env%para_env%group)
               ind = ind + 1
            END DO
         END IF
         IF (ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
            ind = 1
            DO i = 1, 2
               IF (mixed_cdft%dlb_control%recv_work_repl(i)) THEN
                  ALLOCATE (mixed_cdft%dlb_control%recv_info(i)%matrix_info(6))
                  CALL mp_irecv(msgout=mixed_cdft%dlb_control%recv_info(i)%matrix_info, &
                                source=mixed_cdft%source_list(i), &
                                request=req_total(ind), comm=force_env%para_env%group)
                  ind = ind + 1
               END IF
            END DO
         END IF
         IF (ASSOCIATED(req_total)) THEN
            CALL mp_waitall(req_total)
         END IF
         ! Now communicate which processor handles which grid points
         IF (mixed_cdft%dlb_control%send_work) THEN
            ind = COUNT(mixed_cdft%dlb_control%recv_work_repl) + 1
            DO i = 1, 2
               IF (i == 2) &
                  mixed_cdft%dlb_control%target_list(3, :) = mixed_cdft%dlb_control%target_list(3, :) + 3*max_targets
               CALL mp_isend(msgin=mixed_cdft%dlb_control%target_list, &
                             dest=mixed_cdft%dest_list(i), &
                             request=req_total(ind), comm=force_env%para_env%group)
               ind = ind + 1
            END DO
         END IF
         IF (ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
            ind = 1
            DO i = 1, 2
               IF (mixed_cdft%dlb_control%recv_work_repl(i)) THEN
                  ALLOCATE (mixed_cdft%dlb_control%recv_info(i)% &
                            target_list(3, mixed_cdft%dlb_control%recv_info(i)%matrix_info(1)))
                  CALL mp_irecv(msgout=mixed_cdft%dlb_control%recv_info(i)%target_list, &
                                source=mixed_cdft%source_list(i), &
                                request=req_total(ind), comm=force_env%para_env%group)
                  ind = ind + 1
               END IF
            END DO
         END IF
         IF (ASSOCIATED(req_total)) THEN
            CALL mp_waitall(req_total)
            DEALLOCATE (req_total)
         END IF
         IF (ASSOCIATED(sendbuffer)) DEALLOCATE (sendbuffer)
      ELSE
         IF (mixed_cdft%dlb_control%send_work) THEN
            ALLOCATE (req_total(COUNT(mixed_cdft%dlb_control%recv_work_repl) + 2*COUNT(touched)))
         ELSE IF (ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
            ALLOCATE (req_total(COUNT(mixed_cdft%dlb_control%recv_work_repl)))
         END IF
         IF (mixed_cdft%dlb_control%send_work) THEN
            ind = COUNT(mixed_cdft%dlb_control%recv_work_repl)
            DO j = 1, SIZE(mixed_cdft%dest_list)
               IF (touched(j)) THEN
                  ALLOCATE (sbuff(j)%iv(4 + 3*nsend_proc(j)))
                  sbuff(j)%iv(1:4) = (/bo_conf(1, 2), bo_conf(2, 2), bo_conf(1, 3), bo_conf(2, 3)/)
                  offset = 5
                  DO i = 1, SIZE(mixed_cdft%dlb_control%target_list, 2)
                     IF (mixed_cdft%dlb_control%target_list(4 + 2*(j - 1), i) .NE. uninitialized) THEN
                        sbuff(j)%iv(offset:offset + 2) = (/mixed_cdft%dlb_control%target_list(1, i), &
                                                           mixed_cdft%dlb_control%target_list(4 + 2*(j - 1), i), &
                                                           mixed_cdft%dlb_control%target_list(4 + 2*j - 1, i)/)
                        offset = offset + 3
                     END IF
                  END DO
                  DO ispecial = 1, my_special_work
                     CALL mp_isend(msgin=sbuff(j)%iv, &
                                   dest=mixed_cdft%dest_list(j) + (ispecial - 1)*force_env%para_env%num_pe/2, &
                                   request=req_total(ind + ispecial), comm=force_env%para_env%group)
                  END DO
                  ind = ind + my_special_work
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
                  ALLOCATE (mixed_cdft%dlb_control%recv_info(j)%matrix_info(4 + 3*nrecv(j)))
                  CALL mp_irecv(mixed_cdft%dlb_control%recv_info(j)%matrix_info, &
                                source=mixed_cdft%source_list(j), &
                                request=req_total(ind), comm=force_env%para_env%group)
                  ind = ind + 1
               END IF
            END DO
         END IF
         IF (ASSOCIATED(req_total)) THEN
            CALL mp_waitall(req_total)
            DEALLOCATE (req_total)
         END IF
         IF (ANY(mask_send)) THEN
            ALLOCATE (tmp(SIZE(mixed_cdft%dest_list) - COUNT(mask_send)), &
                      tmp_bo(2, SIZE(mixed_cdft%dest_list) - COUNT(mask_send)))
            i = 1
            DO j = 1, SIZE(mixed_cdft%dest_list)
               IF (.NOT. mask_send(j)) THEN
                  tmp(i) = mixed_cdft%dest_list(j)
                  tmp_bo(1:2, i) = mixed_cdft%dest_list_bo(1:2, j)
                  i = i + 1
               END IF
            END DO
            DEALLOCATE (mixed_cdft%dest_list, mixed_cdft%dest_list_bo)
            ALLOCATE (mixed_cdft%dest_list(SIZE(tmp)), mixed_cdft%dest_list_bo(2, SIZE(tmp)))
            mixed_cdft%dest_list = tmp
            mixed_cdft%dest_list_bo = tmp_bo
            DEALLOCATE (tmp, tmp_bo)
         END IF
         IF (ANY(mask_recv)) THEN
            ALLOCATE (tmp(SIZE(mixed_cdft%source_list) - COUNT(mask_recv)), &
                      tmp_bo(4, SIZE(mixed_cdft%source_list) - COUNT(mask_recv)))
            i = 1
            DO j = 1, SIZE(mixed_cdft%source_list)
               IF (.NOT. mask_recv(j)) THEN
                  tmp(i) = mixed_cdft%source_list(j)
                  tmp_bo(1:4, i) = mixed_cdft%source_list_bo(1:4, j)
                  i = i + 1
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
                                        "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
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

      INTEGER :: handle, i, iatom, icomm, iforce_eval, index, iounit, ip, ispecial, iwork, j, &
         jatom, jcomm, k, my_special_work, my_work, natom, nbuffs, nforce_eval, np(3), &
         nsent_total, nskipped, nwork, offset, offset_repl
      INTEGER, DIMENSION(:), POINTER                     :: req_recv, req_total, work, work_dlb
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
      TYPE(cdft_control_type), POINTER                   :: cdft_control
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
               print_section, cdft_control)
      CALL timeset(routineN, handle)
      nforce_eval = SIZE(force_env%sub_force_env)
      CALL force_env_get(force_env=force_env, force_env_section=force_env_section)
      print_section => section_vals_get_subs_vals(force_env_section, "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
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
      cdft_control => mixed_cdft%cdft_control
      CALL pw_env_get(pw_env=mixed_cdft%pw_env, auxbas_pw_pool=auxbas_pw_pool)
      np = auxbas_pw_pool%pw_grid%npts
      dr = auxbas_pw_pool%pw_grid%dr
      shift = -REAL(MODULO(np, 2), dp)*dr/2.0_dp
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
                  ALLOCATE (req_recv(3*(SIZE(mixed_cdft%dlb_control%recv_info(1)%target_list, 2) + &
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
                     nbuffs = nbuffs + (SIZE(mixed_cdft%dlb_control%recv_info(j)%matrix_info) - 4)/3
                  END IF
               END DO
               ALLOCATE (req_recv(3*nbuffs))
            END IF
            DO j = 1, SIZE(mixed_cdft%dlb_control%recv_work_repl)
               IF (mixed_cdft%dlb_control%recv_work_repl(j)) THEN
                  IF (.NOT. mixed_cdft%is_special) THEN
                     offset = 0
                     index = j + (j/2)
                     ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(SIZE(mixed_cdft%dlb_control%recv_info(j)%target_list, 2)))
                     DO i = 1, SIZE(mixed_cdft%dlb_control%recv_info(j)%target_list, 2)
                        IF (mixed_cdft%is_pencil) THEN
                           ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                     weight(mixed_cdft%dlb_control%recv_info(j)%matrix_info(4) + offset: &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(4) + offset + &
                                            (mixed_cdft%dlb_control%recv_info(j)%target_list(2, i) - 1), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(5): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(6), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(2): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(3)))
                           ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                     cavity(mixed_cdft%dlb_control%recv_info(j)%matrix_info(4) + offset: &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(4) + offset + &
                                            (mixed_cdft%dlb_control%recv_info(j)%target_list(2, i) - 1), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(5): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(6), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(2): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(3)))
                           ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                     gradients(3*natom, &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(4) + offset: &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(4) + offset + &
                                               (mixed_cdft%dlb_control%recv_info(j)%target_list(2, i) - 1), &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(5): &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(6), &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(2): &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(3)))
                        ELSE
                           ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                     weight(mixed_cdft%dlb_control%recv_info(j)%matrix_info(5): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(6), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(4) + offset: &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(4) + offset + &
                                            (mixed_cdft%dlb_control%recv_info(j)%target_list(2, i) - 1), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(2): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(3)))
                           ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                     cavity(mixed_cdft%dlb_control%recv_info(j)%matrix_info(5): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(6), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(4) + offset: &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(4) + offset + &
                                            (mixed_cdft%dlb_control%recv_info(j)%target_list(2, i) - 1), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(2): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(3)))
                           ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                     gradients(3*natom, &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(5): &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(6), &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(4) + offset: &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(4) + offset + &
                                               (mixed_cdft%dlb_control%recv_info(j)%target_list(2, i) - 1), &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(2): &
                                               mixed_cdft%dlb_control%recv_info(j)%matrix_info(3)))
                        END IF

                        CALL mp_irecv(msgout=mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%cavity, &
                                      source=mixed_cdft%dlb_control%recv_info(j)%target_list(1, i), &
                                      request=req_recv(3*(i - 1) + (j - 1)*offset_repl + 1), &
                                      comm=force_env%para_env%group, &
                                      tag=mixed_cdft%dlb_control%recv_info(j)%target_list(3, i))
                        CALL mp_irecv(msgout=mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight, &
                                      source=mixed_cdft%dlb_control%recv_info(j)%target_list(1, i), &
                                      request=req_recv(3*(i - 1) + (j - 1)*offset_repl + 2), &
                                      comm=force_env%para_env%group, &
                                      tag=mixed_cdft%dlb_control%recv_info(j)%target_list(3, i) + 1)
                        CALL mp_irecv(msgout=mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%gradients, &
                                      source=mixed_cdft%dlb_control%recv_info(j)%target_list(1, i), &
                                      request=req_recv(3*(i - 1) + (j - 1)*offset_repl + 3), &
                                      comm=force_env%para_env%group, &
                                      tag=mixed_cdft%dlb_control%recv_info(j)%target_list(3, i) + 2)
                        offset = offset + mixed_cdft%dlb_control%recv_info(j)%target_list(2, i)
                     END DO
                     DEALLOCATE (mixed_cdft%dlb_control%recv_info(j)%matrix_info)
                  ELSE
                     ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)% &
                               buffs((SIZE(mixed_cdft%dlb_control%recv_info(j)%matrix_info) - 4)/3))
                     index = 6
                     DO i = 1, SIZE(mixed_cdft%dlb_control%recvbuff(j)%buffs)
                        ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                  weight(mixed_cdft%dlb_control%recv_info(j)%matrix_info(index): &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(index + 1), &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(1): &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(2), &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(3): &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)))
                        ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                  cavity(mixed_cdft%dlb_control%recv_info(j)%matrix_info(index): &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(index + 1), &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(1): &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(2), &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(3): &
                                         mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)))
                        ALLOCATE (mixed_cdft%dlb_control%recvbuff(j)%buffs(i)% &
                                  gradients(3*natom, mixed_cdft%dlb_control%recv_info(j)%matrix_info(index): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(index + 1), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(1): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(2), &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(3): &
                                            mixed_cdft%dlb_control%recv_info(j)%matrix_info(4)))
                        CALL mp_irecv(msgout=mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%cavity, &
                                      source=mixed_cdft%dlb_control%recv_info(j)%matrix_info(index - 1), &
                                      request=req_recv(offset_repl), &
                                      comm=force_env%para_env%group, tag=1)
                        CALL mp_irecv(msgout=mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%weight, &
                                      source=mixed_cdft%dlb_control%recv_info(j)%matrix_info(index - 1), &
                                      request=req_recv(offset_repl + 1), &
                                      comm=force_env%para_env%group, tag=2)
                        CALL mp_irecv(msgout=mixed_cdft%dlb_control%recvbuff(j)%buffs(i)%gradients, &
                                      source=mixed_cdft%dlb_control%recv_info(j)%matrix_info(index - 1), &
                                      request=req_recv(offset_repl + 2), &
                                      comm=force_env%para_env%group, tag=3)
                        index = index + 3
                        offset_repl = offset_repl + 3
                     END DO
                     DEALLOCATE (mixed_cdft%dlb_control%recv_info(j)%matrix_info)
                  END IF
               END IF
            END DO
         END IF
      END IF
      cutoffs => cdft_control%becke_control%cutoffs
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
               gradients => cdft_control%group(1)%gradients
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
                              nsent(jcomm, icomm) = nsent(jcomm, icomm) + 1
                              nsent_total = nsent_total + 1
                              IF (nsent_total == SIZE(nsent, 1)*SIZE(nsent, 2)) should_communicate = .FALSE.
                           END IF
                           IF (ALL(completed(:, icomm))) THEN
                              IF (MODULO(icomm, 3) == 1) THEN
                                 IF (.NOT. mixed_cdft%is_special) THEN
                                    DEALLOCATE (mixed_cdft%dlb_control%cavity)
                                 ELSE
                                    DEALLOCATE (mixed_cdft%dlb_control%sendbuff((icomm - 1)/3 + 1)%cavity)
                                 END IF
                              ELSE IF (MODULO(icomm, 3) == 2) THEN
                                 IF (.NOT. mixed_cdft%is_special) THEN
                                    DEALLOCATE (mixed_cdft%dlb_control%weight)
                                 ELSE
                                    DEALLOCATE (mixed_cdft%dlb_control%sendbuff((icomm - 1)/3 + 1)%weight)
                                 END IF
                              ELSE
                                 IF (.NOT. mixed_cdft%is_special) THEN
                                    DEALLOCATE (mixed_cdft%dlb_control%gradients)
                                 ELSE
                                    DEALLOCATE (mixed_cdft%dlb_control%sendbuff((icomm - 1)/3 + 1)%gradients)
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
                     IF (cdft_control%becke_control%cavity_confine) THEN
                        IF (cavity(k, j, i) < cdft_control%becke_control%eps_cavity) CYCLE
                     END IF
                     grid_p(1) = k*dr(1) + shift(1)
                     grid_p(2) = j*dr(2) + shift(2)
                     grid_p(3) = i*dr(3) + shift(3)
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
                           IF (cdft_control%becke_control%should_skip) THEN
                              IF (is_constraint(iatom)) nskipped = nskipped + 1
                              IF (nskipped == cdft_control%natoms) THEN
                                 IF (in_memory) THEN
                                    IF (cdft_control%becke_control%cavity_confine) THEN
                                       cavity(k, j, i) = 0.0_dp
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
                              dist_vec = (r - grid_p) - ANINT((r - grid_p)/cell_v)*cell_v
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
                              r(ip) = MODULO(r(ip), cell%hmat(ip, ip)) - cell%hmat(ip, ip)/2._dp
                           END DO
                           dist_vec = (r - grid_p) - ANINT((r - grid_p)/cell_v)*cell_v
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
                                       dist_vec = (r1 - grid_p) - ANINT((r1 - grid_p)/cell_v)*cell_v
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
                                       r1(ip) = MODULO(r1(ip), cell%hmat(ip, ip)) - cell%hmat(ip, ip)/2._dp
                                    END DO
                                    dist_vec = (r1 - grid_p) - ANINT((r1 - grid_p)/cell_v)*cell_v
                                    dist2 = SQRT(DOT_PRODUCT(dist_vec, dist_vec))
                                 END IF
                                 IF (in_memory) THEN
                                    IF (store_vectors) THEN
                                       dr1_r2 = pair_dist_vecs(:, iatom, jatom)
                                    ELSE
                                       dr1_r2 = (r - r1) - ANINT((r - r1)/cell_v)*cell_v
                                    END IF
                                    IF (dist2 .LE. th) dist2 = th
                                    tmp_const = (R12(iatom, jatom)**3)
                                    dr_ij_dR(:) = dr1_r2(:)/tmp_const
                                    !derivativ w.r.t. Rj
                                    dr_j_dR = dist_vec(:)/dist2
                                    dmy_dR_j(:) = -(dr_j_dR(:)/R12(iatom, jatom) - (dist1 - dist2)*dr_ij_dR(:))
                                    !derivativ w.r.t. Ri
                                    dmy_dR_i(:) = dr_i_dR(:)/R12(iatom, jatom) - (dist1 - dist2)*dr_ij_dR(:)
                                 END IF
                                 my1 = (dist1 - dist2)/R12(iatom, jatom)
                                 IF (cdft_control%becke_control%adjust) THEN
                                    my1_homo = my1
                                    my1 = my1 + &
                                          cdft_control%becke_control%aij(iatom, jatom)*(1.0_dp - my1**2)
                                 END IF
                                 myexp = 1.5_dp*my1 - 0.5_dp*my1**3
                                 IF (in_memory) THEN
                                    dmyexp = 1.5_dp - 1.5_dp*my1**2
                                    tmp_const = (1.5_dp**2)*dmyexp*(1 - myexp**2)* &
                                                (1.0_dp - ((1.5_dp*myexp - 0.5_dp*(myexp**3))**2))

                                    ds_dR_i(:) = -0.5_dp*tmp_const*dmy_dR_i(:)
                                    ds_dR_j(:) = -0.5_dp*tmp_const*dmy_dR_j(:)
                                    IF (cdft_control%becke_control%adjust) THEN
                                       tmp_const = 1.0_dp - 2.0_dp*my1_homo*cdft_control%becke_control%aij(iatom, jatom)
                                       ds_dR_i(:) = ds_dR_i(:)*tmp_const
                                       ds_dR_j(:) = ds_dR_j(:)*tmp_const
                                    END IF
                                 END IF
                                 myexp = 1.5_dp*myexp - 0.5_dp*myexp**3
                                 myexp = 1.5_dp*myexp - 0.5_dp*myexp**3
                                 tmp_const = 0.5_dp*(1.0_dp - myexp)
                                 cell_functions(iatom) = cell_functions(iatom)*tmp_const
                                 IF (in_memory) THEN
                                    IF (ABS(tmp_const) .LE. th) tmp_const = tmp_const + th
                                    dP_i_dRi(:, iatom) = dP_i_dRi(:, iatom) + ds_dR_i(:)/tmp_const
                                    dP_i_dRj(:, iatom, jatom) = ds_dR_j(:)/tmp_const
                                 END IF

                                 IF (dist2 .LE. cutoffs(jatom)) THEN
                                    tmp_const = 0.5_dp*(1.0_dp + myexp)
                                    cell_functions(jatom) = cell_functions(jatom)*tmp_const
                                    IF (in_memory) THEN
                                       IF (ABS(tmp_const) .LE. th) tmp_const = tmp_const + th
                                       dP_i_dRj(:, jatom, iatom) = -ds_dR_i(:)/tmp_const
                                       dP_i_dRi(:, jatom) = dP_i_dRi(:, jatom) - ds_dR_j(:)/tmp_const
                                    END IF
                                 ELSE
                                    skip_me(jatom) = .TRUE.
                                 END IF
                              END IF
                           END DO
                           IF (in_memory) THEN
                              dP_i_dRi(:, iatom) = cell_functions(iatom)*dP_i_dRi(:, iatom)
                              d_sum_Pm_dR(:, iatom) = d_sum_Pm_dR(:, iatom) + dP_i_dRi(:, iatom)
                              IF (is_constraint(iatom)) &
                                 d_sum_const_dR(:, iatom) = d_sum_const_dR(:, iatom) + dP_i_dRi(:, iatom)* &
                                                            coefficients(iatom)
                              DO jatom = 1, natom
                                 IF (jatom .NE. iatom) THEN
                                    IF (jatom < iatom) THEN
                                       IF (.NOT. skip_me(jatom)) THEN
                                          dP_i_dRj(:, iatom, jatom) = cell_functions(iatom)*dP_i_dRj(:, iatom, jatom)
                                          d_sum_Pm_dR(:, jatom) = d_sum_Pm_dR(:, jatom) + dP_i_dRj(:, iatom, jatom)
                                          IF (is_constraint(iatom)) &
                                             d_sum_const_dR(:, jatom) = d_sum_const_dR(:, jatom) + &
                                                                        dP_i_dRj(:, iatom, jatom)* &
                                                                        coefficients(iatom)
                                          CYCLE
                                       END IF
                                    END IF
                                    dP_i_dRj(:, iatom, jatom) = cell_functions(iatom)*dP_i_dRj(:, iatom, jatom)
                                    d_sum_Pm_dR(:, jatom) = d_sum_Pm_dR(:, jatom) + dP_i_dRj(:, iatom, jatom)
                                    IF (is_constraint(iatom)) &
                                       d_sum_const_dR(:, jatom) = d_sum_const_dR(:, jatom) + dP_i_dRj(:, iatom, jatom)* &
                                                                  coefficients(iatom)
                                 END IF
                              END DO
                           END IF
                        ELSE
                           cell_functions(iatom) = 0.0_dp
                           skip_me(iatom) = .TRUE.
                           IF (cdft_control%becke_control%should_skip) THEN
                              IF (is_constraint(iatom)) nskipped = nskipped + 1
                              IF (nskipped == cdft_control%natoms) THEN
                                 IF (in_memory) THEN
                                    IF (cdft_control%becke_control%cavity_confine) THEN
                                       cavity(k, j, i) = 0.0_dp
                                    END IF
                                 END IF
                                 EXIT
                              END IF
                           END IF
                        END IF
                     END DO
                     IF (nskipped == cdft_control%natoms) CYCLE
                     sum_cell_f_constr = 0.0_dp
                     DO ip = 1, cdft_control%natoms
                        sum_cell_f_constr = sum_cell_f_constr + cell_functions(catom(ip))* &
                                            cdft_control%group(1)%coeff(ip)
                     END DO
                     sum_cell_f_all = 0.0_dp
                     nwork = nwork + 1
                     DO ip = 1, natom
                        sum_cell_f_all = sum_cell_f_all + cell_functions(ip)
                     END DO
                     IF (in_memory) THEN
                        DO iatom = 1, natom
                           IF (ABS(sum_cell_f_all) .GT. 0.0_dp) THEN
                              gradients(3*(iatom - 1) + 1:3*(iatom - 1) + 3, k, j, i) = &
                                 d_sum_const_dR(:, iatom)/sum_cell_f_all - sum_cell_f_constr* &
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
                                   tag=mixed_cdft%dlb_control%dest_tags_repl(i) + 1)
                     CALL mp_isend(msgin=mixed_cdft%dlb_control%gradients, &
                                   dest=mixed_cdft%dlb_control%my_dest_repl(i), &
                                   request=req_send(i, 3), comm=force_env%para_env%group, &
                                   tag=mixed_cdft%dlb_control%dest_tags_repl(i) + 2)
                  END DO
                  should_communicate = .TRUE.
                  nsent_total = 0
               ELSE
                  DO i = 1, SIZE(req_send, 1)
                     CALL mp_isend(msgin=mixed_cdft%dlb_control%sendbuff(ispecial)%cavity, &
                                   dest=mixed_cdft%dlb_control%sendbuff(ispecial)%rank(i), &
                                   request=req_send(i, 3*(ispecial - 1) + 1), &
                                   comm=force_env%para_env%group, tag=1)
                     CALL mp_isend(msgin=mixed_cdft%dlb_control%sendbuff(ispecial)%weight, &
                                   dest=mixed_cdft%dlb_control%sendbuff(ispecial)%rank(i), &
                                   request=req_send(i, 3*(ispecial - 1) + 2), &
                                   comm=force_env%para_env%group, tag=2)
                     CALL mp_isend(msgin=mixed_cdft%dlb_control%sendbuff(ispecial)%gradients, &
                                   dest=mixed_cdft%dlb_control%sendbuff(ispecial)%rank(i), &
                                   request=req_send(i, 3*(ispecial - 1) + 3), &
                                   comm=force_env%para_env%group, tag=3)
                  END DO
                  IF (ispecial .EQ. my_special_work) THEN
                     should_communicate = .TRUE.
                     nsent_total = 0
                  END IF
               END IF
               work(mixed_cdft%dlb_control%my_source + 1) = work(mixed_cdft%dlb_control%my_source + 1) + nwork
               work_dlb(force_env%para_env%mepos + 1) = work_dlb(force_env%para_env%mepos + 1) + nwork
            ELSE
               IF (mixed_cdft%dlb) work(force_env%para_env%mepos + 1) = work(force_env%para_env%mepos + 1) + nwork
               IF (mixed_cdft%dlb) work_dlb(force_env%para_env%mepos + 1) = work_dlb(force_env%para_env%mepos + 1) + nwork
            END IF
         END DO ! ispecial
      END DO ! iwork
      ! Load balancing: wait for communication and deallocate sending buffers
      IF (mixed_cdft%dlb) THEN
         IF (mixed_cdft%dlb_control%recv_work .AND. &
             ANY(mixed_cdft%dlb_control%recv_work_repl)) THEN
            ALLOCATE (req_total(SIZE(req_recv) + SIZE(req_send, 1)*SIZE(req_send, 2)))
            index = SIZE(req_recv)
            req_total(1:index) = req_recv
            DO i = 1, SIZE(req_send, 2)
               DO j = 1, SIZE(req_send, 1)
                  index = index + 1
                  req_total(index) = req_send(j, i)
               END DO
            END DO
            CALL mp_waitall(req_total)
            DEALLOCATE (req_total)
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
         mixed_cdft%dlb_control%prediction_error = mixed_cdft%dlb_control%expected_work - work
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
                                        "MIXED%MIXED_CDFT%PRINT%PROGRAM_RUN_INFO")
      CALL timestop(handle)

   END SUBROUTINE mixed_becke_constraint_low

END MODULE mixed_cdft_methods
