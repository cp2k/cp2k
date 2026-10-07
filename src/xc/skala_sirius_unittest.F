!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM skala_sirius_unittest
   USE cell_methods, ONLY: cell_create, init_cell
   USE cell_types, ONLY: cell_release, &
                         cell_type
   USE ieee_arithmetic, ONLY: ieee_is_finite
   USE input_cp2k_pwdft, ONLY: create_pwdft_section
   USE input_cp2k_xc, ONLY: create_xc_section
   USE input_section_types, ONLY: section_release, &
                                  section_type, &
                                  section_vals_create, &
                                  section_vals_get_subs_vals, &
                                  section_vals_release, &
                                  section_vals_type, &
                                  section_vals_val_set
   USE iso_c_binding, ONLY: C_LOC
   USE kinds, ONLY: default_path_length, &
                    default_string_length, &
                    dp
   USE libcp2k, ONLY: cp2k_finalize, &
                      cp2k_init
   USE message_passing, ONLY: mp_comm_type, &
                              mp_comm_world
   USE particle_types, ONLY: particle_type
   USE sirius, ONLY: &
      SIRIUS_NUMBER_TYPE, &
      sirius_add_atom, sirius_add_atom_type, sirius_add_atom_type_radial_function, &
      sirius_add_atom_type_aw_descriptor, sirius_add_atom_type_lo_descriptor, &
      sirius_context_handler, sirius_create_context, sirius_create_ground_state, &
      sirius_create_kset, sirius_find_ground_state, sirius_free_handler, &
      sirius_generate_density, sirius_generate_effective_potential, sirius_get_energy, sirius_ground_state_handler, &
      sirius_get_forces, sirius_get_stress_tensor, sirius_get_total_magnetization, &
      sirius_get_max_num_gkvec, sirius_set_atom_position, sirius_update_ground_state, &
      sirius_import_parameters, sirius_initialize_context, sirius_kpoint_set_handler, &
      sirius_option_set, &
      sirius_set_atom_type_dion, sirius_set_atom_type_paw, sirius_set_atom_type_radial_grid, &
      sirius_set_atom_type_radial_grid_inf, &
      sirius_set_lattice_vectors, sirius_set_parameters
#if defined(__SIRIUS_XC_DERIVATIVE)
   USE sirius, ONLY: sirius_check_external_xc_derivative
#endif
#if defined(__SIRIUS_XC_POINT_GEOMETRY)
   USE sirius, ONLY: sirius_get_xc_geometry_derivatives
#endif
#if defined(__SIRIUS_XC_POINT_RESPONSE)
   USE sirius, ONLY: sirius_check_external_xc_response
#endif
   USE sirius_interface, ONLY: cp_sirius_init
   USE skala_sirius, ONLY: skala_sirius_register_atom_grid, &
                           skala_sirius_register_points, &
                           skala_sirius_release, &
                           skala_sirius_state_type, skala_sirius_update

   IMPLICIT NONE

   INTEGER, PARAMETER :: nrad = 400, nmt = 801, nfree = 3000
   TYPE(mp_comm_type) :: comm
   TYPE(cell_type), POINTER :: cell => NULL()
   TYPE(particle_type), POINTER :: atoms(:) => NULL()
   TYPE(section_type), POINTER :: definition => NULL()
   TYPE(section_type), POINTER :: pwdft_definition => NULL()
   TYPE(section_vals_type), POINTER :: xc => NULL()
   TYPE(section_vals_type), POINTER :: pwdft => NULL(), grid => NULL()
   TYPE(skala_sirius_state_type), POINTER :: state => NULL()
   TYPE(sirius_context_handler) :: ctx
   TYPE(sirius_kpoint_set_handler) :: ks
   TYPE(sirius_ground_state_handler) :: gs
   CHARACTER(len=default_path_length) :: model, paw_setup, lapw_setup
   CHARACTER(len=32) :: side_arg, step_arg, grid_arg, model_mode, device_arg, model_device_arg, response_arg, bands_arg
   CHARACTER(len=32) :: relativity_arg, core_relativity_arg, radial_arg, charge_arg
   CHARACTER(len=3) :: processing_unit = "CPU"
   INTEGER :: i, j, l, mode, mag, nr, nb, x, y, z, ierr, side, first_mode, last_mode, geometry, natoms
   INTEGER :: updated_basis, fresh_basis, max_scf_iterations, num_bands, radial_points
   REAL(dp) :: lattice(3, 3), r(nrad), v(nrad), rho(nrad), orbital(nrad), dion(1, 1), &
               energy, exc, repeat_exc, kpoint(3, 1), reference(3, 2), fd_step, response_step, &
               rmt(nmt), rfree(nfree), rhofree(nfree), soc_dion(2, 2)
   REAL(dp) :: force_steps(5) = [0.0_dp, 2.0e-4_dp, -2.0e-4_dp, 1.0e-4_dp, -1.0e-4_dp]
   REAL(dp), PARAMETER :: direction(3) = [0.7_dp, -0.4_dp, 0.2_dp]
   REAL(dp) :: force_energy(5), force_derivative, force_fd, force_errors(2), entropy, descf, density_tol, energy_tol
   REAL(dp) :: stress(3, 3), strain(3, 3), magnetization(3)
   REAL(dp), TARGET :: extra_charge
   REAL(dp) :: updated_energy, updated_exc
   REAL(dp), ALLOCATABLE :: forces(:, :)
#if defined(__SIRIUS_XC_DERIVATIVE)
   REAL(dp) :: derivatives(3, 2), errors(2, 2)
#endif
#if defined(__SIRIUS_XC_POINT_RESPONSE)
   REAL(dp) :: response_errors(2)
#endif
#if defined(__SIRIUS_XC_POINT_GEOMETRY)
   REAL(dp) :: cell_derivative(3, 3), repeat_cell(3, 3)
   REAL(dp), ALLOCATABLE :: atom_derivative(:, :), repeat_atom(:, :)
#endif
   REAL(dp), ALLOCATABLE :: points(:, :), weights(:), descriptor_weights(:), centres(:, :)
   INTEGER, ALLOCATABLE :: sizes(:)
   LOGICAL :: converged, atom_grid, force_test, stress_test, derivative_test, update_test, spinor_test, orbital_only
   LOGICAL :: spinor_stress, soc_test, paw_soc_test, lapw_soc_test

   CALL get_command_argument(1, model)
   CALL get_command_argument(2, paw_setup)
   CALL get_command_argument(3, side_arg)
   CALL get_command_argument(4, step_arg)
   CALL get_command_argument(5, grid_arg)
   CALL get_command_argument(6, model_mode)
   CALL get_command_argument(7, device_arg)
   CALL get_command_argument(8, model_device_arg)
   CALL get_command_argument(9, response_arg)
   CALL get_command_argument(10, bands_arg)
   CALL get_command_argument(11, relativity_arg)
   CALL get_command_argument(12, lapw_setup)
   CALL get_command_argument(13, core_relativity_arg)
   CALL get_command_argument(14, radial_arg)
   CALL get_command_argument(15, charge_arg)
   extra_charge = 0.0_dp
   IF (LEN_TRIM(charge_arg) > 0) THEN
      READ (charge_arg, *, IOSTAT=ierr) extra_charge
      IF (ierr /= 0) ERROR STOP "Invalid extra charge"
      IF (.NOT. ieee_is_finite(extra_charge)) ERROR STOP "Extra charge must be finite"
   END IF
   IF (LEN_TRIM(relativity_arg) == 0) relativity_arg = "none"
   IF (LEN_TRIM(core_relativity_arg) == 0) core_relativity_arg = "none"
   SELECT CASE (TRIM(relativity_arg))
   CASE ("none", "zora", "iora", "koelling_harmon", "dirac")
   CASE DEFAULT
      ERROR STOP "Unknown LAPW valence relativity"
   END SELECT
   SELECT CASE (TRIM(core_relativity_arg))
   CASE ("none", "zora", "dirac")
   CASE DEFAULT
      ERROR STOP "Core response requires NR, ZORA or Dirac"
   END SELECT
   IF (LEN_TRIM(lapw_setup) == 0 .AND. core_relativity_arg /= "none") THEN
      ERROR STOP "Core relativity requires a supplied LAPW atom file"
   END IF
   num_bands = 4
   IF (LEN_TRIM(bands_arg) > 0) THEN
      READ (bands_arg, *, IOSTAT=ierr) num_bands
      IF (ierr /= 0) ERROR STOP "Invalid band count"
      IF (num_bands < 1) ERROR STOP "Band count must be positive"
   END IF
   IF (TRIM(device_arg) == "gpu") THEN
      processing_unit = "GPU"
   ELSE IF (LEN_TRIM(device_arg) > 0 .AND. TRIM(device_arg) /= "cpu") THEN
      ERROR STOP "Processing unit must be cpu or gpu"
   END IF
   IF (LEN_TRIM(model_device_arg) > 0 .AND. &
       TRIM(model_device_arg) /= "cpu" .AND. TRIM(model_device_arg) /= "gpu") THEN
      ERROR STOP "Neural model device must be cpu or gpu"
   END IF
   IF (LEN_TRIM(model_mode) > 0 .AND. TRIM(model_mode) /= "skala-float32") THEN
      ERROR STOP "Unknown model mode (expected skala-float32 or no mode for the analytic fixture)"
   END IF
   force_test = TRIM(grid_arg) == "atom-grid-forces" .OR. TRIM(grid_arg) == "spinor-atom-grid-forces" .OR. &
                TRIM(grid_arg) == "spinor-lapw-soc-forces"
   stress_test = TRIM(grid_arg) == "atom-grid-stress" .OR. TRIM(grid_arg) == "spinor-atom-grid-stress" .OR. &
                 TRIM(grid_arg) == "spinor-lapw-soc-stress"
   derivative_test = force_test .OR. stress_test
   update_test = TRIM(grid_arg) == "atom-grid-update"
   paw_soc_test = TRIM(grid_arg) == "spinor-paw-soc-orbitals"
   lapw_soc_test = TRIM(grid_arg) == "spinor-lapw-soc-orbitals" .OR. &
                   TRIM(grid_arg) == "spinor-lapw-soc-forces" .OR. TRIM(grid_arg) == "spinor-lapw-soc-stress"
   soc_test = TRIM(grid_arg) == "spinor-soc-orbitals" .OR. paw_soc_test
   orbital_only = TRIM(grid_arg) == "atom-grid-orbitals" .OR. &
                  TRIM(grid_arg) == "spinor-orbitals" .OR. TRIM(grid_arg) == "spinor-paw-orbitals" .OR. &
                  TRIM(grid_arg) == "spinor-atom-grid-orbitals" .OR. soc_test .OR. (lapw_soc_test .AND. .NOT. derivative_test)
   spinor_stress = TRIM(grid_arg) == "spinor-stress"
   spinor_test = TRIM(grid_arg) == "spinor" .OR. TRIM(grid_arg) == "spinor-atom-grid" .OR. &
                 TRIM(grid_arg) == "spinor-atom-grid-forces" .OR. TRIM(grid_arg) == "spinor-atom-grid-stress" .OR. &
                 (orbital_only .AND. TRIM(grid_arg) /= "atom-grid-orbitals") .OR. spinor_stress .OR. lapw_soc_test
   atom_grid = TRIM(grid_arg) == "atom-grid" .OR. TRIM(grid_arg) == "atom-grid-orbitals" .OR. &
               derivative_test .OR. update_test .OR. &
               TRIM(grid_arg) == "spinor-atom-grid" .OR. TRIM(grid_arg) == "spinor-atom-grid-orbitals" .OR. &
               spinor_stress .OR. paw_soc_test .OR. lapw_soc_test
   IF (LEN_TRIM(grid_arg) > 0 .AND. .NOT. (atom_grid .OR. spinor_test)) ERROR STOP "Unknown atom-grid test mode"
   radial_points = MERGE(30, 100, derivative_test)
   IF (spinor_stress .OR. orbital_only) radial_points = 12
   IF (LEN_TRIM(radial_arg) > 0) THEN
      IF (.NOT. atom_grid) ERROR STOP "Radial point count requires an atom grid"
      READ (radial_arg, *, IOSTAT=ierr) radial_points
      IF (ierr /= 0 .OR. radial_points < 2) ERROR STOP "Radial point count must be at least two"
   END IF
#if !defined(__SIRIUS_XC_POINT_UPDATE)
   IF (update_test .OR. derivative_test .OR. spinor_stress) ERROR STOP "SIRIUS lacks the point-quadrature update API"
#endif
#if !defined(__SIRIUS_XC_POINT_GEOMETRY)
   IF (derivative_test .OR. spinor_stress) ERROR STOP "SIRIUS lacks the point geometry callbacks"
#endif
   fd_step = 0.0_dp
   IF (LEN_TRIM(step_arg) > 0) THEN
      READ (step_arg, *, IOSTAT=ierr) fd_step
      IF (ierr /= 0) ERROR STOP "Invalid finite-difference step"
      IF (.NOT. ieee_is_finite(fd_step) .OR. fd_step <= 0.0_dp) ERROR STOP "Finite-difference step must be positive"
#if !defined(__SIRIUS_XC_DERIVATIVE)
      IF (.NOT. derivative_test) ERROR STOP "SIRIUS lacks the external XC orbital derivative diagnostic"
#endif
   END IF
   response_step = MERGE(0.0_dp, fd_step, derivative_test)
   IF (LEN_TRIM(response_arg) > 0) THEN
      READ (response_arg, *, IOSTAT=ierr) response_step
      IF (ierr /= 0) ERROR STOP "Invalid field-response step"
      IF (.NOT. ieee_is_finite(response_step) .OR. response_step <= 0.0_dp) THEN
         ERROR STOP "Field-response step must be positive"
      END IF
   END IF
   side = 6
   IF (LEN_TRIM(side_arg) > 0) THEN
      READ (side_arg, *, IOSTAT=ierr) side
      IF (ierr /= 0) ERROR STOP "Invalid quadrature dimension"
   END IF
   IF (side < 4 .OR. side > 64) ERROR STOP "Quadrature dimension must be between 4 and 64"
   IF (orbital_only .AND. fd_step <= 0.0_dp) ERROR STOP "The orbital-only mode requires a finite-difference step"
   IF (spinor_stress .AND. fd_step <= 0.0_dp) ERROR STOP "The spinor stress mode requires a finite-difference step"
   first_mode = 0
   last_mode = 1
   IF (LEN_TRIM(paw_setup) > 0) first_mode = 1
   IF (TRIM(paw_setup) == "LAPW") THEN
      first_mode = 2
      last_mode = 2
   END IF
   IF (first_mode /= 2 .AND. relativity_arg /= "none") ERROR STOP "Valence relativity requires the LAPW fixture"
   IF (relativity_arg == "dirac" .AND. lapw_soc_test) ERROR STOP "Dirac includes SOC without a separate correction"
   IF (first_mode /= 2 .AND. LEN_TRIM(lapw_setup) > 0) ERROR STOP "Atom file requires the LAPW fixture"
   IF (spinor_test) THEN
      IF (lapw_soc_test .AND. first_mode /= 2) ERROR STOP "LAPW SOC fixture requires LAPW as the second argument"
      IF (TRIM(grid_arg) == "spinor-paw-orbitals" .OR. paw_soc_test) first_mode = 1
      IF (soc_test .AND. first_mode /= 0 .AND. .NOT. paw_soc_test) THEN
         ERROR STOP "SOC fixture requires norm-conserving pseudopotentials or the synthetic PAW setup"
      END IF
      IF (paw_soc_test .AND. LEN_TRIM(paw_setup) > 0) ERROR STOP "PAW SOC fixture requires the synthetic setup"
      IF (spinor_stress .AND. first_mode /= 0) ERROR STOP "Isolated XC stress fixture requires norm-conserving PW orbitals"
      last_mode = first_mode
   END IF
   reference = 0.0_dp
   IF (LEN_TRIM(model) == 0 .OR. LEN_TRIM(model) > default_string_length) THEN
      ERROR STOP "Pass the test model path, or SKALA with SKALA_MODEL set"
   END IF
   CALL cp2k_init()
   CALL cp_sirius_init()
   comm = mp_comm_world
   CALL comm%init()
   IF (comm%mepos == 0) THEN
      IF (atom_grid) THEN
         WRITE (*, '(2A)') "Periodic atom quadrature; PAW setup: ", TRIM(paw_setup)
      ELSE
         WRITE (*, '(A,I0,2A)') "Points per block: ", side**3, "; PAW setup: ", TRIM(paw_setup)
      END IF
   END IF
   density_tol = MERGE(1.0e-10_dp, 2.0e-6_dp, derivative_test .OR. update_test)
   energy_tol = MERGE(1.0e-11_dp, 2.0e-7_dp, derivative_test .OR. update_test)
   max_scf_iterations = 150
   IF (model_mode == "skala-float32") THEN
      ! Actual-model preset. SIRIUS point callbacks promote the checkpoint to Float64.
      ! Keep the analytic double-precision checks and all comparison bounds unchanged.
      density_tol = 1.0e-8_dp
      energy_tol = 1.0e-9_dp
      max_scf_iterations = 300
      ! Use the actual-model stencil without relaxing the derivative error bound.
      IF (derivative_test) force_steps = 10.0_dp*force_steps
   END IF
   IF (derivative_test .AND. LEN_TRIM(step_arg) > 0) THEN
      ! Keep both geometry stencils inside the same LAPW point domains.
      force_steps = fd_step*[0.0_dp, 1.0_dp, -1.0_dp, 0.5_dp, -0.5_dp]
      fd_step = 0.0_dp
   END IF
   IF (orbital_only .OR. spinor_stress) max_scf_iterations = 1
   ! LAPW already has occupied auxiliary-Hamiltonian orbitals for the operator check.
   IF (orbital_only .AND. first_mode == 2) max_scf_iterations = 0
   IF (comm%mepos == 0) THEN
      WRITE (*, '(A,2ES12.3)') "SCF density and energy thresholds: ", density_tol, energy_tol
      WRITE (*, '(A,I0)') "Maximum SCF iterations: ", max_scf_iterations
   END IF
   CALL create_xc_section(definition)
   CALL section_vals_create(xc, definition)
   CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%MODEL", c_val=TRIM(model))
   CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID", l_val=.TRUE.)
   CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID_USE_CUDA", l_val=TRIM(model_device_arg) == "gpu")
   IF (atom_grid) THEN
      CALL create_pwdft_section(pwdft_definition)
      CALL section_vals_create(pwdft, pwdft_definition)
      CALL section_vals_val_set(pwdft, "SKALA_GRID%N_RADIAL", i_val=radial_points)
      CALL section_vals_val_set(pwdft, "SKALA_GRID%N_ANGULAR", i_val=MERGE(50, 110, derivative_test))
      IF (spinor_stress .OR. orbital_only) THEN
         ! These modes test the fixed-orbital chain rule, not grid convergence.
         CALL section_vals_val_set(pwdft, "SKALA_GRID%N_ANGULAR", i_val=26)
      END IF
      CALL section_vals_val_set(pwdft, "SKALA_GRID%EPS_RESPONSE", &
                                r_val=MERGE(1.0e-6_dp, 1.0e-8_dp, model_mode == "skala-float32"))
      IF (spinor_test) THEN
         CALL section_vals_val_set(pwdft, "SKALA_GRID%N_SPIN_ANGULAR", &
                                   i_val=MERGE(194, 26, extra_charge /= 0.0_dp .AND. model_mode == "skala-float32"))
      END IF
      grid => section_vals_get_subs_vals(pwdft, "SKALA_GRID")
   END IF
   lattice = 0.0_dp
   lattice(1, 1) = 8.0_dp
   lattice(2, 2) = 9.0_dp
   lattice(3, 3) = 10.0_dp
   CALL cell_create(cell, hmat=lattice, periodic=[1, 1, 1])
   natoms = MERGE(2, 1, derivative_test .OR. (update_test .AND. first_mode /= 2))
   strain = RESHAPE([0.3_dp, 0.1_dp, -0.07_dp, 0.1_dp, -0.2_dp, 0.06_dp, -0.07_dp, 0.06_dp, 0.15_dp], [3, 3])
   ALLOCATE (atoms(natoms), forces(3, natoms))
#if defined(__SIRIUS_XC_POINT_GEOMETRY)
   ALLOCATE (atom_derivative(3, natoms), repeat_atom(3, natoms))
#endif
   atoms(1)%r = [1.6_dp, 2.7_dp, 4.0_dp]
   IF (natoms == 2) atoms(2)%r = [3.84_dp, 4.59_dp, 5.9_dp]
   IF (derivative_test .AND. first_mode == 2) atoms(2)%r = [5.44_dp, 5.29_dp, 6.6_dp]
   nr = 0
   nb = 0
   IF (comm%mepos == 0) THEN
      nr = side**3
      nb = 1
   END IF
   ALLOCATE (points(3, nr), weights(nr), descriptor_weights(nr), centres(3, nb), sizes(nb))
   IF (nb > 0) THEN
      i = 0
      DO z = 0, side - 1
         DO y = 0, side - 1
            DO x = 0, side - 1
               i = i + 1
               points(:, i) = [8*(0.037_dp + REAL(x, dp)/side), 9*(0.029_dp + REAL(y, dp)/side), &
                               10*(0.041_dp + REAL(z, dp)/side)]
            END DO
         END DO
      END DO
      sizes = nr
      centres(:, 1) = atoms(1)%r
   END IF
   weights = cell%deth/REAL(side**3, dp)
   descriptor_weights = weights
   IF (nr > 0) weights(1) = 0.0_dp
   kpoint(:, 1) = [0.07_dp, 0.03_dp, 0.02_dp]

   DO mode = first_mode, last_mode
      DO mag = MERGE(3, MERGE(MIN(mode, 1), 0, derivative_test .OR. update_test), spinor_test), &
         MERGE(3, MERGE(MIN(mode, 1), 1, derivative_test .OR. update_test), spinor_test)
         DO geometry = 1, MERGE(5, MERGE(2, 1, update_test), derivative_test)
            atoms(1)%r = [1.6_dp, 2.7_dp, 4.0_dp] + force_steps(geometry)*direction
            IF (update_test) CALL set_update_geometry(geometry == 2)
            IF (stress_test) THEN
               lattice = 0.0_dp
               lattice(1, 1) = 8.0_dp
               lattice(2, 2) = 9.0_dp
               lattice(3, 3) = 10.0_dp
               CALL init_cell(cell, hmat=lattice)
               atoms(1)%r = MATMUL(lattice, [0.2_dp, 0.3_dp, 0.4_dp])
               atoms(2)%r = MATMUL(lattice, [0.48_dp, 0.51_dp, 0.59_dp])
               IF (mode == 2) atoms(2)%r = MATMUL(lattice, [0.68_dp, 5.29_dp/9.0_dp, 0.66_dp])
            END IF
            IF (.NOT. derivative_test .OR. geometry == 1) THEN
               CALL sirius_create_context(comm%get_handle(), ctx)
               IF (mode == 2) THEN
                  CALL sirius_import_parameters(ctx, &
                                                '{"parameters":{"electronic_structure_method":"full_potential_lapwlo",'// &
                                                '"use_symmetry":false,"pw_cutoff":5.0,"gk_cutoff":1.5,"num_bands":4,'// &
                                                '"xc_functionals":[],"smearing_width":0.01,"auto_rmt":0,'// &
                                                '"core_relativity":"none","valence_relativity":"none",'// &
                                                '"lmax_apw":3,"lmax_rho":6,"lmax_pot":6},'// &
                                                '"control":{"mpi_grid_dims":[1,1],"processing_unit":"'// &
                                                processing_unit//'","verbosity":1},'// &
                                                '"mixer":{"type":"linear","beta":0.2,"use_hartree":false},'// &
                                                '"iterative_solver":{"type":"exact","min_tolerance":1e-11}}')
                  IF (derivative_test) THEN
                     IF (.NOT. spinor_test) CALL sirius_import_parameters(ctx, '{"control":{"use_second_variation":false}}')
                     CALL sirius_set_parameters(ctx, pw_cutoff=9.0_dp)
                  END IF
                  IF (spinor_test) CALL sirius_import_parameters(ctx, '{"control":{"use_second_variation":true}}')
                  IF (relativity_arg == "dirac") CALL sirius_import_parameters(ctx, '{"control":{"use_second_variation":false}}')
                  CALL sirius_set_parameters(ctx, valence_rel=TRIM(relativity_arg))
                  CALL sirius_set_parameters(ctx, core_rel=TRIM(core_relativity_arg))
               ELSE
                  CALL sirius_import_parameters(ctx, &
                                                '{"parameters":{"electronic_structure_method":"pseudopotential",'// &
                                                '"use_symmetry":false,"pw_cutoff":8.0,"gk_cutoff":2.6,"num_bands":4,'// &
                                                '"xc_functionals":[],"smearing_width":0.01},'// &
                                                '"control":{"mpi_grid_dims":[1,1],"processing_unit":"'// &
                                                processing_unit//'","verbosity":1},'// &
                                                '"mixer":{"type":"linear","beta":0.2,"use_hartree":false},'// &
                                                '"iterative_solver":{"type":"exact","min_tolerance":1e-11,"num_steps":80}}')
               END IF
               CALL sirius_set_parameters(ctx, num_mag_dims=mag, num_bands=num_bands)
               CALL sirius_option_set(ctx, "parameters", "extra_charge", SIRIUS_NUMBER_TYPE, C_LOC(extra_charge), error_code=ierr)
               IF (ierr /= 0) ERROR STOP "SIRIUS could not set the extra charge"
               IF (spinor_test .AND. extra_charge /= 0.0_dp .AND. model_mode == "skala-float32") THEN
                  CALL sirius_import_parameters(ctx, '{"mixer":{"type":"broyden2","beta":0.05,"max_history":8}}')
               END IF
               IF (orbital_only .AND. mag == 1 .AND. relativity_arg /= "dirac") THEN
                  CALL sirius_import_parameters(ctx, '{"parameters":{"fixed_mag":0.5}}')
               END IF
               IF (lapw_soc_test) CALL sirius_import_parameters(ctx, '{"parameters":{"so_correction":true}}')
               IF (spinor_test .AND. mode /= 2) CALL sirius_import_parameters(ctx, '{"iterative_solver":{"type":"davidson"}}')
               CALL sirius_set_lattice_vectors(ctx, lattice(:, 1), lattice(:, 2), lattice(:, 3))
               IF (mode == 2 .AND. LEN_TRIM(lapw_setup) > 0) THEN
                  CALL sirius_add_atom_type(ctx, "He", fname=TRIM(lapw_setup))
               ELSE IF (mode == 2) THEN
                  CALL sirius_add_atom_type(ctx, "He", zn=MERGE(1, 2, spinor_test))
                  DO i = 1, nmt
                     rmt(i) = 1.0e-6_dp + (1.8_dp - 1.0e-6_dp)*(REAL(i - 1, dp)/REAL(nmt - 1, dp))**2
                  END DO
                  DO i = 1, nfree
                     rfree(i) = 1.0e-7_dp + (20.0_dp - 1.0e-7_dp)*(REAL(i - 1, dp)/REAL(nfree - 1, dp))**3
                  END DO
                  rhofree = 16.0_dp/ACOS(-1.0_dp)*EXP(-4*rfree)
                  IF (spinor_test) rhofree = EXP(-2*rfree)/ACOS(-1.0_dp)
                  CALL sirius_set_atom_type_radial_grid(ctx, "He", nmt, rmt)
                  CALL sirius_set_atom_type_radial_grid_inf(ctx, "He", nfree, rfree)
                  CALL sirius_add_atom_type_radial_function(ctx, "He", "ae_rho", rhofree, nfree)
                  DO l = 0, 3
                     CALL sirius_add_atom_type_aw_descriptor(ctx, "He", l + 1, l, -0.3_dp, 0, .FALSE.)
                     CALL sirius_add_atom_type_aw_descriptor(ctx, "He", l + 1, l, -0.3_dp, 1, .FALSE.)
                  END DO
                  CALL sirius_add_atom_type_lo_descriptor(ctx, "He", 1, 1, 0, -0.7_dp, 0, .FALSE.)
                  CALL sirius_add_atom_type_lo_descriptor(ctx, "He", 1, 1, 0, -0.7_dp, 1, .FALSE.)
               ELSE IF (mode == 1 .AND. LEN_TRIM(paw_setup) > 0) THEN
                  CALL sirius_add_atom_type(ctx, "He", fname=TRIM(paw_setup))
               ELSE
                  ! An odd electron count exercises spinor occupations without a closed-shell restriction.
                  CALL sirius_add_atom_type(ctx, "He", zn=MERGE(1, 2, spinor_test), spin_orbit=soc_test)
                  DO i = 1, nrad
                     r(i) = 20.0_dp*(REAL(i - 1, dp)/REAL(nrad - 1, dp))**2
                     IF (mode == 1) r(i) = 0.001_dp + 2.999_dp*(REAL(i - 1, dp)/REAL(nrad - 1, dp))**2
                  END DO
                  v = -REAL(MERGE(1, 2, spinor_test), dp)/SQRT(1.0_dp + r**2)
                  rho = 2.0_dp*REAL(MERGE(1, 2, spinor_test), dp)*EXP(-r**2)*r
                  IF (spinor_test .AND. mode == 1) rho = rho/(EXP(-r(1)**2) - EXP(-r(nrad)**2))
                  orbital = EXP(-r)
                  CALL sirius_set_atom_type_radial_grid(ctx, "He", nrad, r)
                  CALL sirius_add_atom_type_radial_function(ctx, "He", "vloc", v, nrad)
                  CALL sirius_add_atom_type_radial_function(ctx, "He", "ps_rho_total", rho, nrad)
                  CALL sirius_add_atom_type_radial_function(ctx, "He", "ps_atomic_wf", orbital, nrad, n=1, l=0)
                  dion = 0.0_dp
                  IF (mode == 1) THEN
                     ! Deliberately synthetic setup: tests the coupled implementation, not PAW accuracy.
                     IF (soc_test) THEN
                        DO i = 1, 2
                           l = MERGE(-1, 1, i == 1)
                           CALL sirius_add_atom_type_radial_function(ctx, "He", "beta", 0.2_dp*r**2*EXP(-2*r**2), nrad, l=l)
                           CALL sirius_add_atom_type_radial_function(ctx, "He", "ae_paw_wf", r**2*EXP(-r**2), nrad)
                           CALL sirius_add_atom_type_radial_function(ctx, "He", "ps_paw_wf", r**2*EXP(-0.8_dp*r**2), nrad)
                        END DO
                        DO j = 1, 2
                           DO i = 1, j
                              DO l = 0, 2, 2
                                 CALL sirius_add_atom_type_radial_function(ctx, "He", "q_aug", &
                                                                           0.05_dp*r**2*EXP(-2*r**2)* &
                                                                           (1.0_dp + 0.4_dp*(i - 1) + 0.2_dp*(j - 1))/ &
                                                                           REAL(l + 1, dp), nrad, l=l, idxrf1=i, idxrf2=j)
                              END DO
                           END DO
                        END DO
                     ELSE
                        CALL sirius_add_atom_type_radial_function(ctx, "He", "beta", 0.2_dp*r*EXP(-2*r**2), nrad, l=0)
                        CALL sirius_add_atom_type_radial_function(ctx, "He", "ae_paw_wf", r*EXP(-r**2), nrad)
                        CALL sirius_add_atom_type_radial_function(ctx, "He", "ps_paw_wf", r*EXP(-0.8_dp*r**2), nrad)
                        CALL sirius_add_atom_type_radial_function(ctx, "He", "q_aug", 0.05_dp*r**2*EXP(-2*r**2), nrad, &
                                                                  l=0, idxrf1=1, idxrf2=1)
                     END IF
                     CALL sirius_add_atom_type_radial_function(ctx, "He", "ae_paw_core", 0.02_dp*EXP(-3*r**2), nrad)
                     CALL sirius_add_atom_type_radial_function(ctx, "He", "ae_paw_core_tau", 0.03_dp*EXP(-3*r**2), nrad)
                     CALL sirius_add_atom_type_radial_function(ctx, "He", "ps_rho_core", 0.01_dp*EXP(-2*r**2), nrad)
                     IF (soc_test) THEN
                        soc_dion = 0.0_dp
                        soc_dion(1, 1) = 0.7_dp
                        soc_dion(2, 2) = -0.4_dp
                        CALL sirius_set_atom_type_dion(ctx, "He", 2, soc_dion)
                        CALL sirius_set_atom_type_paw(ctx, "He", 0.0_dp, [1.0_dp, 1.0_dp], 2)
                     ELSE
                        CALL sirius_set_atom_type_dion(ctx, "He", 1, dion)
                        CALL sirius_set_atom_type_paw(ctx, "He", 0.0_dp, [2.0_dp], 1)
                     END IF
                  ELSE IF (soc_test) THEN
                     CALL sirius_add_atom_type_radial_function(ctx, "He", "beta", r**2*EXP(-r**2), nrad, l=-1)
                     CALL sirius_add_atom_type_radial_function(ctx, "He", "beta", r**2*EXP(-r**2), nrad, l=1)
                     soc_dion = 0.0_dp
                     soc_dion(1, 1) = 0.7_dp
                     soc_dion(2, 2) = -0.4_dp
                     CALL sirius_set_atom_type_dion(ctx, "He", 2, soc_dion)
                  ELSE
                     CALL sirius_set_atom_type_dion(ctx, "He", 0, dion(:0, :0))
                  END IF
               END IF
               DO i = 1, natoms
                  IF (spinor_test) THEN
                     CALL sirius_add_atom(ctx, "He", MATMUL(cell%h_inv, atoms(i)%r), vector_field=[0.3_dp, -0.2_dp, 0.1_dp])
                  ELSE IF (orbital_only .AND. mag == 1) THEN
                     CALL sirius_add_atom(ctx, "He", MATMUL(cell%h_inv, atoms(i)%r), vector_field=[0.0_dp, 0.0_dp, 0.3_dp])
                  ELSE
                     CALL sirius_add_atom(ctx, "He", MATMUL(cell%h_inv, atoms(i)%r))
                  END IF
               END DO
               IF (atom_grid) THEN
                  CALL skala_sirius_register_atom_grid(ctx, xc, atoms, cell, comm, grid, state)
               ELSE
                 CALL skala_sirius_register_points(ctx, xc, atoms, cell, points, weights, descriptor_weights, sizes, centres, state)
               END IF
               CALL sirius_initialize_context(ctx)
               CALL sirius_create_kset(ctx, 1, kpoint, [1.0_dp], .TRUE., ks)
               CALL sirius_create_ground_state(ks, gs)
            END IF
            IF (derivative_test .AND. geometry > 1) THEN
               ! Preserve both density and orbital index sets throughout the strain stencil.
               IF (stress_test) THEN
                  lattice = lattice + force_steps(geometry)*MATMUL(strain, lattice)
                  CALL init_cell(cell, hmat=lattice)
                  atoms(1)%r = MATMUL(lattice, [0.2_dp, 0.3_dp, 0.4_dp])
                  atoms(2)%r = MATMUL(lattice, [0.48_dp, 0.51_dp, 0.59_dp])
                  IF (mode == 2) atoms(2)%r = MATMUL(lattice, [0.68_dp, 5.29_dp/9.0_dp, 0.66_dp])
               END IF
               CALL sirius_set_lattice_vectors(ctx, lattice(:, 1), lattice(:, 2), lattice(:, 3))
               DO i = 1, natoms
                  CALL sirius_set_atom_position(ctx, i, MATMUL(cell%h_inv, atoms(i)%r))
               END DO
               CALL skala_sirius_update(ctx, state, atoms, cell, comm)
               CALL sirius_update_ground_state(gs)
            END IF
            converged = .FALSE.
            CALL sirius_find_ground_state(gs, density_tol=density_tol, energy_tol=energy_tol, &
                                          iter_solver_tol=MERGE(1.0e-11_dp, 1.0e-9_dp, derivative_test .OR. update_test), &
                                          initial_guess=(.NOT. derivative_test .OR. geometry == 1), &
                                          max_niter=max_scf_iterations, save_state=.FALSE., converged=converged, error_code=ierr)
            IF (ierr /= 0) ERROR STOP "SIRIUS Skala SCF failed"
            IF (orbital_only .AND. mode == 2) THEN
               ! Keep density and tau tied to the current test orbitals.
               CALL sirius_generate_density(gs, add_core=.TRUE., transform_to_rg=.TRUE., efermi=.FALSE., error_code=ierr)
               IF (ierr /= 0) ERROR STOP "SIRIUS spinor density generation failed"
               CALL sirius_generate_effective_potential(gs, error_code=ierr)
               IF (ierr /= 0) ERROR STOP "SIRIUS spinor potential generation failed"
            END IF
            CALL sirius_get_energy(gs, "total", energy)
            CALL sirius_get_energy(gs, "descf", descf)
            energy = energy + descf
            CALL sirius_get_energy(gs, "exc", exc)
            IF (comm%mepos == 0) WRITE (*, '(A,2I3,L3,2ES22.12)') "Point SCF (method 0=PW 1=PAW 2=LAPW,mag,converged,E,Exc): ", &
               mode, mag, converged, energy, exc
            IF ((.NOT. converged .AND. .NOT. (orbital_only .OR. spinor_stress)) .OR. &
                .NOT. ieee_is_finite(energy) .OR. .NOT. ieee_is_finite(exc)) THEN
               ERROR STOP "SIRIUS Skala point SCF did not converge with finite energies"
            END IF
            IF (ABS(exc) < 1.0e-6_dp) ERROR STOP "The model must produce a nontrivial XC energy"
            IF (spinor_test) THEN
               CALL sirius_get_total_magnetization(gs, magnetization(1), error_code=ierr)
               IF (ierr /= 0 .OR. .NOT. ALL(ieee_is_finite(magnetization))) ERROR STOP "Invalid spinor magnetization"
               IF (extra_charge /= 0.0_dp .AND. SQRT(SUM(magnetization**2)) < 1.0e-4_dp) THEN
                  ERROR STOP "Vacuous open-shell spinor derivative check"
               END IF
               IF (comm%mepos == 0) WRITE (*, '(A,3ES19.10)') "Spinor magnetization (x,y,z): ", magnetization
            END IF
            IF (orbital_only .AND. mag == 1 .AND. relativity_arg /= "dirac") THEN
               CALL sirius_get_total_magnetization(gs, magnetization(1), error_code=ierr)
               IF (ierr /= 0 .OR. .NOT. ALL(ieee_is_finite(magnetization))) ERROR STOP "Invalid collinear magnetization"
               ! The collinear API returns its single z component in the first slot.
               IF (ABS(magnetization(1)) < 1.0e-4_dp) ERROR STOP "Vacuous collinear orbital derivative check"
               IF (comm%mepos == 0) WRITE (*, '(A,ES19.10)') "Collinear magnetization: ", magnetization(1)
            END IF
            IF (geometry == 1 .AND. .NOT. spinor_test) reference(mode + 1, mag + 1) = energy
            IF (derivative_test) THEN
               CALL sirius_get_energy(gs, "demet", entropy)
               ! At fractional occupations forces differentiate the free energy E-TS.
               force_energy(geometry) = energy + entropy
               IF (geometry == 1) THEN
                  IF (stress_test) THEN
                     CALL sirius_get_stress_tensor(gs, "total", stress, error_code=ierr)
                     IF (ierr /= 0 .OR. .NOT. ALL(ieee_is_finite(stress))) ERROR STOP "SIRIUS Skala stress failed"
                     force_derivative = cell%deth*SUM(stress*strain)
                  ELSE
                     CALL sirius_get_forces(gs, "total", forces, error_code=ierr)
                     IF (ierr /= 0 .OR. .NOT. ALL(ieee_is_finite(forces))) ERROR STOP "SIRIUS Skala total forces failed"
                     force_derivative = -DOT_PRODUCT(forces(:, 1), direction)
                     IF (MAXVAL(ABS(SUM(forces, DIM=2))) > 2.0e-6_dp) ERROR STOP "Total forces violate translation invariance"
                  END IF
               END IF
            END IF
            CALL sirius_generate_effective_potential(gs)
            CALL sirius_get_energy(gs, "exc", repeat_exc)
            IF (ABS(exc - repeat_exc) > 1.0e-7_dp) ERROR STOP "Unchanged density changed the point XC energy"
            exc = repeat_exc
#if defined(__SIRIUS_XC_POINT_GEOMETRY)
            IF (atom_grid .AND. .NOT. derivative_test .AND. .NOT. update_test .AND. &
                .NOT. spinor_stress .AND. .NOT. orbital_only) THEN
               CALL sirius_get_xc_geometry_derivatives(gs, natoms, atom_derivative, cell_derivative, error_code=ierr)
               IF (ierr /= 0) ERROR STOP "SIRIUS Skala geometry callback failed"
               IF (.NOT. ALL(ieee_is_finite(atom_derivative)) .OR. .NOT. ALL(ieee_is_finite(cell_derivative))) THEN
                  ERROR STOP "Nonfinite explicit XC geometry derivatives"
               END IF
               IF (MAXVAL(ABS(cell_derivative)) < 1.0e-8_dp) ERROR STOP "Vacuous XC cell derivative check"
               IF (comm%mepos == 0) WRITE (*, '(A,*(ES19.10))') "Explicit XC atom/cell derivatives: ", &
                  atom_derivative, cell_derivative
               CALL sirius_get_xc_geometry_derivatives(gs, natoms, repeat_atom, repeat_cell, error_code=ierr)
               IF (ierr /= 0) ERROR STOP "Repeated geometry callback failed"
               IF (comm%mepos == 0) WRITE (*, '(A,2ES19.10)') "XC geometry warmup differences: ", &
                  MAXVAL(ABS(repeat_atom - atom_derivative)), MAXVAL(ABS(repeat_cell - cell_derivative))
               ! TorchScript can specialize its graph on the first geometry-gradient evaluation.
               atom_derivative = repeat_atom
               cell_derivative = repeat_cell
               CALL sirius_get_xc_geometry_derivatives(gs, natoms, repeat_atom, repeat_cell, error_code=ierr)
               IF (ierr /= 0) ERROR STOP "Third geometry callback failed"
               IF (comm%mepos == 0) WRITE (*, '(A,2ES19.10)') "Repeated XC atom/cell derivative differences: ", &
                  MAXVAL(ABS(repeat_atom - atom_derivative)), MAXVAL(ABS(repeat_cell - cell_derivative))
               IF (MAXVAL(ABS(repeat_atom - atom_derivative)) > 1.0e-12_dp .OR. &
                   MAXVAL(ABS(repeat_cell - cell_derivative)) > 1.0e-12_dp) ERROR STOP "Geometry callback is not repeatable"
               CALL sirius_get_energy(gs, "exc", repeat_exc)
               IF (ABS(exc - repeat_exc) > 1.0e-7_dp) ERROR STOP "Geometry callback changed the ground-state XC energy"
               IF (comm%mepos == 0) atoms(1)%r(1) = 1.601_dp
               CALL sirius_get_xc_geometry_derivatives(gs, natoms, repeat_atom, repeat_cell, error_code=ierr)
               IF (ierr == 0) ERROR STOP "Geometry callback accepted a stale CP2K atom grid"
               CALL sirius_generate_effective_potential(gs, error_code=ierr)
               IF (ierr == 0) ERROR STOP "Energy callback accepted a stale CP2K atom grid"
               atoms(1)%r(1) = 1.6_dp
               CALL sirius_generate_effective_potential(gs, error_code=ierr)
               IF (ierr /= 0) ERROR STOP "Restoring the atom grid did not restore energy evaluation"
            END IF
#endif
            IF (update_test) THEN
               IF (geometry == 1) THEN
                  CALL sirius_get_max_num_gkvec(ks, updated_basis)
                  CALL set_update_geometry(.TRUE.)
                  CALL sirius_set_lattice_vectors(ctx, lattice(:, 1), lattice(:, 2), lattice(:, 3))
                  DO i = 1, natoms
                     CALL sirius_set_atom_position(ctx, i, MATMUL(cell%h_inv, atoms(i)%r))
                  END DO
                  CALL skala_sirius_update(ctx, state, atoms, cell, comm)
                  CALL sirius_update_ground_state(gs)
                  CALL skala_sirius_update(ctx, state, atoms, cell, comm)
                  converged = .FALSE.
                  CALL sirius_find_ground_state(gs, density_tol=density_tol, energy_tol=energy_tol, &
                                                iter_solver_tol=1.0e-11_dp, initial_guess=.FALSE., &
                                             max_niter=max_scf_iterations, save_state=.FALSE., converged=converged, error_code=ierr)
                  IF (ierr /= 0 .OR. .NOT. converged) ERROR STOP "Updated geometry SCF did not converge"
                  CALL sirius_get_energy(gs, "total", updated_energy)
                  CALL sirius_get_energy(gs, "descf", descf)
                  updated_energy = updated_energy + descf
                  CALL sirius_get_energy(gs, "exc", updated_exc)
                  IF (.NOT. ieee_is_finite(updated_energy) .OR. .NOT. ieee_is_finite(updated_exc)) THEN
                     ERROR STOP "Nonfinite energy after geometry update"
                  END IF
                  IF (ABS(updated_energy - energy) < 1.0e-7_dp) ERROR STOP "Vacuous geometry-update check"
#if defined(__SIRIUS_XC_POINT_GEOMETRY)
                  CALL sirius_get_xc_geometry_derivatives(gs, natoms, atom_derivative, cell_derivative, error_code=ierr)
                  IF (ierr /= 0) ERROR STOP "Updated geometry callback failed"
#endif
               ELSE
                  CALL sirius_get_max_num_gkvec(ks, fresh_basis)
                  IF (updated_basis /= fresh_basis) ERROR STOP "Geometry-update fixture crosses the orbital cutoff"
                  IF (comm%mepos == 0) WRITE (*, '(A,2I3,2ES19.10)') "Updated/fresh energy differences (method,mag,E,Exc): ", &
                     mode, mag, ABS(updated_energy - energy), ABS(updated_exc - exc)
                  IF (ABS(updated_energy - energy) > 2.0e-7_dp .OR. ABS(updated_exc - exc) > 2.0e-7_dp) THEN
                     ERROR STOP "Updated and fresh geometry energies differ"
                  END IF
#if defined(__SIRIUS_XC_POINT_GEOMETRY)
                  CALL sirius_get_xc_geometry_derivatives(gs, natoms, repeat_atom, repeat_cell, error_code=ierr)
                  IF (ierr /= 0) ERROR STOP "Fresh geometry callback failed"
                  IF (.NOT. ALL(ieee_is_finite(atom_derivative)) .OR. .NOT. ALL(ieee_is_finite(cell_derivative)) .OR. &
                      .NOT. ALL(ieee_is_finite(repeat_atom)) .OR. .NOT. ALL(ieee_is_finite(repeat_cell))) THEN
                     ERROR STOP "Nonfinite geometry derivatives in update comparison"
                  END IF
                  IF (comm%mepos == 0) WRITE (*, '(A,2ES19.10)') "Updated/fresh XC geometry derivative differences: ", &
                     MAXVAL(ABS(repeat_atom - atom_derivative)), MAXVAL(ABS(repeat_cell - cell_derivative))
                  IF (MAXVAL(ABS(repeat_atom - atom_derivative)) > 2.0e-6_dp .OR. &
                      MAXVAL(ABS(repeat_cell - cell_derivative)) > 2.0e-6_dp) ERROR STOP "Updated geometry derivatives differ"
#endif
               END IF
            END IF
#if defined(__SIRIUS_XC_DERIVATIVE)
            IF (spinor_stress) CALL check_spinor_stress()
            IF (fd_step > 0.0_dp .AND. .NOT. spinor_stress) THEN
               CALL sirius_check_external_xc_derivative(gs, fd_step, derivatives, error_code=ierr)
               IF (ierr /= 0) ERROR STOP "SIRIUS Skala orbital derivative check failed"
               DO i = 1, 2
                  errors(:, i) = ABS(derivatives(2:3, i) - derivatives(1, i))
                  IF (comm%mepos == 0) WRITE (*, '(A,I0,6ES19.10)') "XC orbital derivative (direction,h,H,FDh,FDh/2,errors): ", &
                     i, fd_step, derivatives(:, i), errors(:, i)
               END DO
               IF (.NOT. ALL(ieee_is_finite(derivatives))) ERROR STOP "Nonfinite XC orbital derivatives"
               IF (MAXVAL(errors) > 2.0e-5_dp) ERROR STOP "XC orbital derivative mismatch"
               IF (MAXVAL(ABS(derivatives(1, :))) < 1.0e-5_dp) ERROR STOP "Vacuous orbital derivative check"
               CALL sirius_get_energy(gs, "exc", repeat_exc)
               IF (exc /= repeat_exc) ERROR STOP "Orbital diagnostic changed the ground-state XC energy"
               CALL sirius_generate_effective_potential(gs)
               CALL sirius_get_energy(gs, "exc", repeat_exc)
               IF (ABS(exc - repeat_exc) > 1.0e-7_dp) ERROR STOP "Orbital diagnostic did not restore the callback"
            END IF
#endif
#if defined(__SIRIUS_XC_POINT_RESPONSE)
            IF (mode == 2 .AND. response_step > 0.0_dp) THEN
               CALL sirius_check_external_xc_response(gs, response_step, response_errors, error_code=ierr)
               IF (ierr /= 0) ERROR STOP "SIRIUS Skala field response check failed"
               IF (comm%mepos == 0) WRITE (*, '(A,2ES19.10)') "XC field/point response relative errors: ", response_errors
               IF (.NOT. ALL(ieee_is_finite(response_errors))) ERROR STOP "Nonfinite XC field response errors"
               IF (MAXVAL(response_errors) > MERGE(1.0e-4_dp, 1.0e-8_dp, model_mode == "skala-float32")) THEN
                  ERROR STOP "XC field response mismatch"
               END IF
               CALL sirius_get_energy(gs, "exc", repeat_exc)
               IF (exc /= repeat_exc) ERROR STOP "Field response check changed the ground-state XC energy"
            END IF
#endif
            IF (.NOT. derivative_test .OR. geometry == 5) THEN
               CALL sirius_free_handler(gs)
               CALL sirius_free_handler(ks)
               CALL sirius_free_handler(ctx)
               CALL skala_sirius_release(state)
               IF (ASSOCIATED(state)) ERROR STOP "Callback state was not released"
            END IF
         END DO
         IF (derivative_test) THEN
            DO i = 2, 4, 2
               force_fd = (force_energy(i) - force_energy(i + 1))/(2*force_steps(i))
               force_errors(i/2) = ABS(force_fd - force_derivative)
               IF (comm%mepos == 0) WRITE (*, '(3A,2I3,4ES19.10)') "Total ", MERGE("strain", "force ", stress_test), &
                  " derivative (method,mag,h,analytic,FD,error): ", &
                  mode, mag, force_steps(i), force_derivative, force_fd, ABS(force_fd - force_derivative)
            END DO
            IF (.NOT. ALL(ieee_is_finite(force_errors)) .OR. MAXVAL(force_errors) > 2.0e-5_dp) THEN
               ERROR STOP "Geometry derivative disagrees with the finite-difference free energy"
            END IF
         END IF
      END DO
   END DO
   IF (.NOT. derivative_test .AND. .NOT. update_test .AND. .NOT. orbital_only .AND. &
       ANY(ABS(reference(:, 1) - reference(:, 2)) > 2.0e-5_dp)) THEN
      ERROR STOP "Scalar and unpolarized collinear SCF energies differ"
   END IF
   DEALLOCATE (points, weights, descriptor_weights, centres, sizes, atoms, forces)
#if defined(__SIRIUS_XC_POINT_GEOMETRY)
   DEALLOCATE (atom_derivative, repeat_atom)
#endif
   CALL cell_release(cell)
   CALL section_vals_release(xc)
   CALL section_release(definition)
   IF (atom_grid) THEN
      NULLIFY (grid)
      CALL section_vals_release(pwdft)
      CALL section_release(pwdft_definition)
   END IF
   CALL cp2k_finalize()
CONTAINS
! **************************************************************************************************
!> \brief Check fixed-orbital XC strain derivatives against finite differences.
! **************************************************************************************************
   SUBROUTINE check_spinor_stress()
      INTEGER                                            :: j
      REAL(dp)                                           :: analytic, base(3, 3), energies(7), &
                                                            finite_difference(3), fourth_order, &
                                                            restored, steps(7), strain_error(2)

      base = lattice
      steps = fd_step*[0.0_dp, 1.0_dp, -1.0_dp, 0.5_dp, -0.5_dp, 0.25_dp, -0.25_dp]
      CALL frozen_spinor_energy(base, 0.0_dp, energies(1))
      CALL sirius_get_stress_tensor(gs, "xc", stress, error_code=ierr)
      IF (ierr /= 0 .OR. .NOT. ALL(ieee_is_finite(stress))) ERROR STOP "Spinor XC stress failed"
      analytic = cell%deth*SUM(stress*strain)
      IF (ABS(analytic) < 1.0e-6_dp) ERROR STOP "Vacuous spinor XC stress check"
      DO j = 2, 7
         CALL frozen_spinor_energy(base, steps(j), energies(j))
      END DO
      DO j = 1, 3
         finite_difference(j) = (energies(2*j) - energies(2*j + 1))/(2*steps(2*j))
         IF (comm%mepos == 0) WRITE (*, '(A,4ES19.10)') "Frozen spinor XC strain (h,analytic,FD,error): ", &
            steps(2*j), analytic, finite_difference(j), ABS(finite_difference(j) - analytic)
      END DO
      DO j = 1, 2
         fourth_order = (4*finite_difference(j + 1) - finite_difference(j))/3
         IF (comm%mepos == 0) WRITE (*, '(A,4ES19.10)') "Frozen spinor XC four-point strain (h,analytic,FD,error): ", &
            steps(2*j + 2), analytic, fourth_order, ABS(fourth_order - analytic)
         strain_error(j) = ABS(fourth_order - analytic)
      END DO
      IF (.NOT. ALL(ieee_is_finite(strain_error)) .OR. MAXVAL(strain_error) > 2.0e-5_dp) THEN
         ERROR STOP "Spinor XC stress disagrees with the fixed-orbital energy derivative"
      END IF
      CALL frozen_spinor_energy(base, 0.0_dp, restored)
      IF (ABS(restored - energies(1)) > 1.0e-7_dp) ERROR STOP "Spinor stress check did not restore the XC energy"
   END SUBROUTINE check_spinor_stress

! **************************************************************************************************
!> \brief Evaluate XC under strain without updating orbitals or occupations.
!> \param base Reference cell vectors.
!> \param step Dimensionless strain amplitude.
!> \param value XC energy at the strained geometry.
! **************************************************************************************************
   SUBROUTINE frozen_spinor_energy(base, step, value)
      REAL(dp), INTENT(IN)                               :: base(3, 3), step
      REAL(dp), INTENT(OUT)                              :: value

      lattice = base + step*MATMUL(strain, base)
      CALL init_cell(cell, hmat=lattice)
      atoms(1)%r = MATMUL(lattice, [0.2_dp, 0.3_dp, 0.4_dp])
      CALL sirius_set_lattice_vectors(ctx, lattice(:, 1), lattice(:, 2), lattice(:, 3))
      CALL sirius_set_atom_position(ctx, 1, [0.2_dp, 0.3_dp, 0.4_dp])
      CALL skala_sirius_update(ctx, state, atoms, cell, comm)
      CALL sirius_update_ground_state(gs)
      ! Do not solve the Hamiltonian or refit occupations along the strain stencil.
      CALL sirius_generate_density(gs, add_core=.TRUE., transform_to_rg=.TRUE., efermi=.FALSE., error_code=ierr)
      IF (ierr /= 0) ERROR STOP "Fixed-orbital spinor density failed"
      CALL sirius_generate_effective_potential(gs, error_code=ierr)
      IF (ierr /= 0) ERROR STOP "Fixed-orbital spinor potential failed"
      CALL sirius_get_energy(gs, "exc", value)
      IF (.NOT. ieee_is_finite(value)) ERROR STOP "Nonfinite fixed-orbital spinor XC energy"
   END SUBROUTINE frozen_spinor_energy

! **************************************************************************************************
!> \brief Set the reference geometry or a small coupled cell and atom displacement.
!> \param moved Select the displaced geometry instead of the reference geometry.
! **************************************************************************************************
   SUBROUTINE set_update_geometry(moved)
      LOGICAL, INTENT(IN)                                :: moved

      lattice = 0.0_dp
      lattice(1, 1) = 8.0_dp
      lattice(2, 2) = 9.0_dp
      lattice(3, 3) = 10.0_dp
      atoms(1)%r = [1.6_dp, 2.7_dp, 4.0_dp]
      IF (moved) THEN
         ! Cell strain must not change the cutoff-selected basis in this comparison.
         lattice(1, 1) = 8.001_dp
         lattice(1, 2) = 0.0003_dp
         atoms(1)%r = atoms(1)%r + 0.01_dp*direction
      END IF
      CALL init_cell(cell, lattice)
   END SUBROUTINE set_update_geometry
END PROGRAM skala_sirius_unittest
