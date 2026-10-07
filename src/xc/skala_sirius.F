!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Connect SIRIUS primitive XC fields to the CP2K Skala evaluator.
! **************************************************************************************************
MODULE skala_sirius
   USE ISO_C_BINDING, ONLY: C_ASSOCIATED, C_DOUBLE, C_F_POINTER, C_FUNLOC, C_INT, C_LOC, C_PTR
   USE cell_types, ONLY: cell_type
   USE ieee_arithmetic, ONLY: ieee_is_finite
   USE input_section_types, ONLY: section_vals_type, section_vals_val_get
   USE kinds, ONLY: dp
   USE lebedev, ONLY: get_number_of_lebedev_grid, lebedev_grid
   USE message_passing, ONLY: mp_comm_type
   USE particle_types, ONLY: particle_type
   USE skala_atom_grid, ONLY: skala_atom_grid_build, skala_atom_grid_pullback, skala_atom_grid_direction
   USE skala_gpw_features, ONLY: skala_grid_layout_type
   USE skala_gpw_functional, ONLY: skala_grid_eval, skala_point_eval
#if defined(__SIRIUS)
   USE sirius, ONLY: sirius_context_handler
#if defined(__SIRIUS_XC_CALLBACK) && defined(__LIBTORCH)
   USE sirius, ONLY: sirius_set_xc_callback
#endif
#if defined(__SIRIUS_XC_POINT_CALLBACK) && defined(__LIBTORCH)
   USE sirius, ONLY: sirius_get_parameters, sirius_set_xc_point_callback
#endif
#if defined(__SIRIUS_XC_SPINOR) && defined(__LIBTORCH)
   USE sirius, ONLY: sirius_set_xc_spinor_mode
#endif
#if defined(__SIRIUS_XC_POINT_GEOMETRY) && defined(__LIBTORCH)
   USE sirius, ONLY: sirius_set_xc_point_geometry_callbacks
#endif
#if defined(__SIRIUS_XC_POINT_RESPONSE) && defined(__LIBTORCH)
   USE sirius, ONLY: sirius_set_xc_point_response_callback
#endif
#if defined(__SIRIUS_XC_POINT_HESSIAN) && defined(__LIBTORCH)
   USE sirius, ONLY: sirius_set_xc_point_hessian_callback
#endif
#if defined(__SIRIUS_XC_RESPONSE_TOLERANCE) && defined(__LIBTORCH)
   USE sirius, ONLY: sirius_set_xc_response_tolerance
#endif
#if defined(__SIRIUS_XC_QUADRATURE_DIRECTION) && defined(__LIBTORCH)
   USE sirius, ONLY: sirius_set_xc_quadrature_direction_callback
#endif
#if defined(__SIRIUS_XC_POINT_UPDATE) && defined(__LIBTORCH)
   USE sirius, ONLY: sirius_update_xc_point_set
#endif
#endif
#include "../base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_sirius'
   PUBLIC :: skala_sirius_state_type, skala_sirius_release
#if defined(__SIRIUS)
   PUBLIC :: skala_sirius_register, skala_sirius_register_points, skala_sirius_register_atom_grid, skala_sirius_update
#endif

   TYPE, BIND(C) :: callback_token_type
      INTEGER(C_INT) :: id = 0
   END TYPE callback_token_type

   TYPE skala_sirius_state_type
      PRIVATE
      TYPE(callback_token_type) :: token
      TYPE(cell_type), POINTER :: cell => NULL()
      TYPE(particle_type), POINTER :: particles(:) => NULL()
      TYPE(section_vals_type), POINTER :: xc_section => NULL()
      INTEGER :: nradial = 0, nangular = 0, shell = 0
      INTEGER :: nspin_angular = 194
      REAL(dp) :: scale = 0.0_dp
      REAL(dp) :: grid_lattice(3, 3) = 0.0_dp
      REAL(dp), ALLOCATABLE :: grid_centres(:, :)
      TYPE(skala_sirius_state_type), POINTER :: next => NULL()
   END TYPE skala_sirius_state_type

   ! Tokens cross the C ABI; the Fortran pointers remain in this per-context registry.
   TYPE(skala_sirius_state_type), POINTER, SAVE :: states => NULL()
   INTEGER(C_INT), SAVE :: next_id = 0

CONTAINS

#if defined(__SIRIUS)
! **************************************************************************************************
!> \brief Construct a periodic atom grid and register the joint field callback.
!> \param context SIRIUS context.
!> \param xc_section XC input owned by the PWDFT environment.
!> \param particles Current particle set.
!> \param cell Current cell.
!> \param group Context communicator.
!> \param grid_section PW_DFT/SKALA_GRID input.
!> \param state Per-context callback state.
! **************************************************************************************************
   SUBROUTINE skala_sirius_register_atom_grid(context, xc_section, particles, cell, group, grid_section, state)
      TYPE(sirius_context_handler), INTENT(IN)           :: context
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_comm_type), INTENT(IN)                     :: group
      TYPE(section_vals_type), POINTER                   :: grid_section
      TYPE(skala_sirius_state_type), POINTER             :: state

      INTEGER                                            :: a, error, na, nr, nspin, shell, spin_grid
      INTEGER, ALLOCATABLE                               :: sizes(:)
      REAL(dp)                                           :: response_tolerance, scale
      REAL(dp), ALLOCATABLE                              :: block_centres(:, :), centres(:, :), &
                                                            descriptor(:), energy(:), points(:, :)

      IF (ANY(cell%perd /= 1)) THEN
         CPABORT("SIRIUS Skala atom quadrature requires a three-dimensional periodic cell")
      END IF
      CALL section_vals_val_get(grid_section, "N_RADIAL", i_val=nr)
      CALL section_vals_val_get(grid_section, "N_ANGULAR", i_val=na)
      CALL section_vals_val_get(grid_section, "N_SPIN_ANGULAR", i_val=nspin)
      CALL section_vals_val_get(grid_section, "IMAGE_SHELL", i_val=shell)
      CALL section_vals_val_get(grid_section, "RADIAL_SCALE", r_val=scale)
      CALL section_vals_val_get(grid_section, "EPS_RESPONSE", r_val=response_tolerance)
      IF (.NOT. ieee_is_finite(response_tolerance) .OR. response_tolerance <= 0.0_dp .OR. &
          response_tolerance >= 1.0_dp) THEN
         CPABORT("SIRIUS Skala EPS_RESPONSE must be between zero and one")
      END IF
      error = 0
      IF (nspin < 6 .OR. nspin > 974) THEN
         error = 1
      ELSE
         spin_grid = get_number_of_lebedev_grid(n=nspin)
         IF (lebedev_grid(spin_grid)%n /= nspin .OR. ANY(lebedev_grid(spin_grid)%w < 0.0_dp)) error = 1
      END IF
      CALL group%max(error)
      IF (error /= 0) THEN
         CPABORT("N_SPIN_ANGULAR requires an available positive Lebedev rule up to 974")
      END IF
      ALLOCATE (centres(3, SIZE(particles)))
      DO a = 1, SIZE(particles)
         centres(:, a) = particles(a)%r
      END DO
      CALL skala_atom_grid_build(group, cell%hmat, centres, nr, na, shell, scale, &
                                 points, energy, descriptor, sizes, block_centres)
      CALL skala_sirius_register_points(context, xc_section, particles, cell, points, energy, descriptor, &
                                        sizes, block_centres, state)
#if defined(__SIRIUS_XC_RESPONSE_TOLERANCE) && defined(__LIBTORCH)
      CALL sirius_set_xc_response_tolerance(context, response_tolerance, error)
      IF (error /= 0) CPABORT("Could not set the SIRIUS Skala response tolerance")
#endif
      state%nradial = nr
      state%nangular = na
      state%nspin_angular = nspin
      state%shell = shell
      state%scale = scale
      state%grid_lattice = cell%hmat
      CALL MOVE_ALLOC(centres, state%grid_centres)
#if defined(__SIRIUS_XC_POINT_GEOMETRY) && defined(__LIBTORCH)
      CALL sirius_set_xc_point_geometry_callbacks(context, C_FUNLOC(evaluate_skala_point_geometry), &
                                                  C_FUNLOC(pullback_skala_quadrature), C_LOC(state%token), error)
      IF (error /= 0) CPABORT("Could not register the SIRIUS Skala geometry callbacks")
#else
      MARK_USED(error)
#endif
#if defined(__SIRIUS_XC_QUADRATURE_DIRECTION) && defined(__LIBTORCH)
      CALL sirius_set_xc_quadrature_direction_callback(context, C_FUNLOC(differentiate_skala_quadrature), &
                                                       C_LOC(state%token), error)
      IF (error /= 0) CPABORT("Could not register the SIRIUS Skala quadrature direction callback")
#endif
   END SUBROUTINE skala_sirius_register_atom_grid

! **************************************************************************************************
!> \brief Register Skala before initializing the SIRIUS context.
!> \param context SIRIUS context.
!> \param xc_section XC input, owned by the PWDFT environment.
!> \param particles Particle set, owned by the PWDFT subsystem.
!> \param cell Cell, owned by the PWDFT subsystem.
!> \param state Callback state, released after destroying the SIRIUS context.
! **************************************************************************************************
   SUBROUTINE skala_sirius_register(context, xc_section, particles, cell, state)
      TYPE(sirius_context_handler), INTENT(IN) :: context
      TYPE(section_vals_type), POINTER :: xc_section
      TYPE(particle_type), POINTER :: particles(:)
      TYPE(cell_type), POINTER :: cell
      TYPE(skala_sirius_state_type), POINTER :: state
#if defined(__SIRIUS_XC_CALLBACK) && defined(__LIBTORCH)
      INTEGER :: error

      CALL create_state(xc_section, particles, cell, state)
      CALL sirius_set_xc_callback(context, C_FUNLOC(evaluate_skala_fields), C_LOC(state%token), error)
      IF (error /= 0) CPABORT("Could not register the SIRIUS Skala callback")
#else
      MARK_USED(context)
      MARK_USED(xc_section)
      MARK_USED(particles)
      MARK_USED(cell)
      MARK_USED(state)
      CPABORT("Skala with SIRIUS requires LibTorch and SIRIUS with the primitive-field XC callback API")
#endif
   END SUBROUTINE skala_sirius_register

! **************************************************************************************************
!> \brief Register a fixed atom-block quadrature for joint reconstructed Skala fields.
!> \param context SIRIUS context, before initialization.
!> \param xc_section XC input, owned by the PWDFT environment.
!> \param particles Particle set, owned by the PWDFT subsystem.
!> \param cell Cell, owned by the PWDFT subsystem.
!> \param points Cartesian points in bohr, owned once across the context communicator.
!> \param energy_weights Energy-partition quadrature weights.
!> \param descriptor_weights Nonlocal descriptor-domain quadrature weights.
!> \param atom_sizes Number of points in each complete local atom block.
!> \param atom_coords Cartesian centres of these blocks in bohr.
!> \param state Callback state, released after destroying the SIRIUS context.
! **************************************************************************************************
   SUBROUTINE skala_sirius_register_points(context, xc_section, particles, cell, points, energy_weights, &
                                           descriptor_weights, atom_sizes, atom_coords, state)
      TYPE(sirius_context_handler), INTENT(IN) :: context
      TYPE(section_vals_type), POINTER :: xc_section
      TYPE(particle_type), POINTER :: particles(:)
      TYPE(cell_type), POINTER :: cell
      REAL(dp), INTENT(IN) :: points(:, :), energy_weights(:), descriptor_weights(:), atom_coords(:, :)
      INTEGER, INTENT(IN) :: atom_sizes(:)
      TYPE(skala_sirius_state_type), POINTER :: state
#if defined(__SIRIUS_XC_POINT_CALLBACK) && defined(__LIBTORCH)
      INTEGER :: nr, nb, error, mag

      nr = SIZE(energy_weights)
      nb = SIZE(atom_sizes)
      CPASSERT(ALL(SHAPE(points) == [3, nr]))
      CPASSERT(ALL(SHAPE(atom_coords) == [3, nb]))
      CPASSERT(SIZE(descriptor_weights) == nr)
      CALL create_state(xc_section, particles, cell, state)
      CALL sirius_set_xc_point_callback(context, C_FUNLOC(evaluate_skala_points), C_LOC(state%token), nr, nb, &
                                        points, energy_weights, descriptor_weights, atom_sizes, atom_coords, error)
      IF (error /= 0) CPABORT("Could not register the SIRIUS Skala point callback")
      CALL sirius_get_parameters(context, num_mag_dims=mag, error_code=error)
      IF (error /= 0) CPABORT("Could not query the SIRIUS spin convention")
      IF (mag == 3) THEN
#if defined(__SIRIUS_XC_SPINOR)
         CALL sirius_set_xc_spinor_mode(context, .TRUE., error)
         IF (error /= 0) CPABORT("Could not register the SIRIUS Skala spinor convention")
#else
         CPABORT("Noncollinear Skala requires the SIRIUS four-component point XC API")
#endif
      END IF
#if defined(__SIRIUS_XC_POINT_HESSIAN)
      CALL sirius_set_xc_point_hessian_callback(context, C_FUNLOC(evaluate_skala_point_hessian), &
                                                C_LOC(state%token), error)
      IF (error /= 0) CPABORT("Could not register the SIRIUS Skala joint Hessian callback")
#elif defined(__SIRIUS_XC_POINT_RESPONSE)
      CALL sirius_set_xc_point_response_callback(context, C_FUNLOC(evaluate_skala_point_response), &
                                                 C_LOC(state%token), error)
      IF (error /= 0) CPABORT("Could not register the SIRIUS Skala response callback")
#endif
#else
      MARK_USED(context)
      MARK_USED(xc_section)
      MARK_USED(particles)
      MARK_USED(cell)
      MARK_USED(points)
      MARK_USED(energy_weights)
      MARK_USED(descriptor_weights)
      MARK_USED(atom_sizes)
      MARK_USED(atom_coords)
      MARK_USED(state)
      CPABORT("Reconstructed Skala fields require LibTorch and the SIRIUS point XC callback API")
#endif
   END SUBROUTINE skala_sirius_register_points
#endif

! **************************************************************************************************
!> \brief Allocate the caller-owned callback state and its C token.
!> \param xc_section XC input owned by the PWDFT environment.
!> \param particles Current particle set.
!> \param cell Current cell.
!> \param state Per-context callback state.
! **************************************************************************************************
   SUBROUTINE create_state(xc_section, particles, cell, state)
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(skala_sirius_state_type), POINTER             :: state

      CPASSERT(.NOT. ASSOCIATED(state))
      CPASSERT(next_id < HUGE(next_id))
      ALLOCATE (state)
      next_id = next_id + 1
      state%token%id = next_id
      state%xc_section => xc_section
      CPASSERT(ASSOCIATED(particles))
      CPASSERT(ASSOCIATED(cell))
      state%particles => particles
      state%cell => cell
      state%next => states
      states => state
   END SUBROUTINE create_state

#if defined(__SIRIUS)
! **************************************************************************************************
!> \brief Refresh the atom grid after setting SIRIUS geometry and before updating its ground state.
!> \param context Initialized SIRIUS context with the new atoms and lattice.
!> \param state Callback state.
!> \param particles Current particle set.
!> \param cell Current cell.
!> \param group Context communicator.
! **************************************************************************************************
   SUBROUTINE skala_sirius_update(context, state, particles, cell, group)
      TYPE(sirius_context_handler), INTENT(IN)           :: context
      TYPE(skala_sirius_state_type), POINTER             :: state
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_comm_type), INTENT(IN)                     :: group

      INTEGER                                            :: changed
#if defined(__SIRIUS_XC_POINT_UPDATE) && defined(__LIBTORCH)
      INTEGER                                            :: a, error
      INTEGER, ALLOCATABLE                               :: sizes(:)
      REAL(dp), ALLOCATABLE                              :: block_centres(:, :), centres(:, :), &
                                                            descriptor(:), energy(:), points(:, :)
#endif

      IF (.NOT. ASSOCIATED(state)) RETURN
      CPASSERT(ASSOCIATED(particles))
      CPASSERT(ASSOCIATED(cell))
      state%particles => particles
      state%cell => cell
      changed = MERGE(0, 1, grid_geometry_current(state))
      CALL group%max(changed)
      IF (changed == 0) RETURN
#if defined(__SIRIUS_XC_POINT_UPDATE) && defined(__LIBTORCH)
      CPASSERT(ALL(cell%perd == 1))
      CPASSERT(SIZE(particles) == SIZE(state%grid_centres, 2))
      ALLOCATE (centres(3, SIZE(particles)))
      DO a = 1, SIZE(particles)
         centres(:, a) = particles(a)%r
      END DO
      CALL skala_atom_grid_build(group, cell%hmat, centres, state%nradial, state%nangular, state%shell, state%scale, &
                                 points, energy, descriptor, sizes, block_centres)
      CALL sirius_update_xc_point_set(context, SIZE(energy), SIZE(sizes), points, energy, descriptor, &
                                      sizes, block_centres, error_code=error)
      IF (error /= 0) CPABORT("Could not update the SIRIUS Skala atom quadrature")
      state%grid_lattice = cell%hmat
      CALL MOVE_ALLOC(centres, state%grid_centres)
#else
      MARK_USED(context)
      CPABORT("Moving Skala atom grids require the SIRIUS point-quadrature update API")
#endif
   END SUBROUTINE skala_sirius_update
#endif

! **************************************************************************************************
!> \brief Release callback state after the associated SIRIUS context has been destroyed.
!> \param state Callback state.
! **************************************************************************************************
   SUBROUTINE skala_sirius_release(state)
      TYPE(skala_sirius_state_type), POINTER             :: state

      TYPE(skala_sirius_state_type), POINTER             :: current, previous

      IF (.NOT. ASSOCIATED(state)) RETURN
      current => states
      NULLIFY (previous)
      DO WHILE (ASSOCIATED(current))
         IF (ASSOCIATED(current, state)) EXIT
         previous => current
         current => current%next
      END DO
      CPASSERT(ASSOCIATED(current))
      IF (ASSOCIATED(previous)) THEN
         previous%next => current%next
      ELSE
         states => current%next
      END IF
      DEALLOCATE (state)
      NULLIFY (state)
   END SUBROUTINE skala_sirius_release

#if defined(__SIRIUS_XC_CALLBACK) && defined(__LIBTORCH)
! **************************************************************************************************
!> \brief SIRIUS C callback. One FFT communicator partitions a complete regular grid.
!> \param user_data Opaque callback token.
!> \param nr Number of local quadrature points.
!> \param ns Number of primitive-field spin components.
!> \param dims Global regular-grid dimensions.
!> \param z_offset Start of the local regular-grid slab.
!> \param fcomm_context Fortran context communicator handle.
!> \param fcomm_fft Fortran FFT communicator handle.
!> \param lattice Cell vectors as columns, in bohr.
!> \param weight Regular-grid volume element.
!> \param rho Density fields.
!> \param gradient Cartesian density gradients.
!> \param tau Kinetic-energy density fields.
!> \param energy Total XC energy over the context communicator.
!> \param vrho Density potential.
!> \param vgradient Density-gradient potential.
!> \param vtau Kinetic-density potential.
!> \param status Zero on success, nonzero on callback failure.
! **************************************************************************************************
   SUBROUTINE evaluate_skala_fields(user_data, nr, ns, dims, z_offset, fcomm_context, fcomm_fft, &
                                    lattice, weight, rho, gradient, tau, energy, vrho, vgradient, vtau, status) BIND(C)
      TYPE(C_PTR), VALUE                                 :: user_data
      INTEGER(C_INT), VALUE                              :: nr, ns
      INTEGER(C_INT), INTENT(IN)                         :: dims(3)
      INTEGER(C_INT), VALUE                              :: z_offset, fcomm_context, fcomm_fft
      REAL(C_DOUBLE), INTENT(IN)                         :: lattice(3, 3)
      REAL(C_DOUBLE), VALUE                              :: weight
      REAL(C_DOUBLE), INTENT(IN)                         :: rho(nr, ns), gradient(nr, 3, ns), &
                                                            tau(nr, ns)
      REAL(C_DOUBLE), INTENT(OUT)                        :: energy, vrho(nr, ns), &
                                                            vgradient(nr, 3, ns), vtau(nr, ns)
      INTEGER(C_INT), INTENT(OUT)                        :: status

      INTEGER                                            :: i, invalid, nz
      TYPE(mp_comm_type), TARGET                         :: context_comm, fft_comm
      TYPE(skala_grid_layout_type)                       :: grid
      TYPE(skala_sirius_state_type), POINTER             :: state

      CALL fft_comm%set_handle(INT(fcomm_fft))
      CALL context_comm%set_handle(INT(fcomm_context))
      energy = 0.0_dp
      vrho = 0.0_dp
      vgradient = 0.0_dp
      vtau = 0.0_dp
      status = 1
      state => callback_state(user_data)
      invalid = 0
      IF (.NOT. ASSOCIATED(state)) invalid = 1
      IF (nr < 0 .OR. ALL(ns /= [1, 2]) .OR. ANY(dims <= 0) .OR. weight <= 0.0_dp) invalid = 1
      CALL context_comm%max(invalid)
      IF (invalid /= 0) RETURN
      nz = nr/(dims(1)*dims(2))
      IF (nr /= dims(1)*dims(2)*nz .OR. z_offset < 0 .OR. z_offset + nz > dims(3)) invalid = 1
      IF (MAXVAL(ABS(lattice - state%cell%hmat)) > 1.0e-10_dp) invalid = 1
      CALL context_comm%max(invalid)
      IF (invalid /= 0) RETURN

      grid%group => fft_comm
      grid%npts = dims
      grid%bounds(1, :) = 0
      grid%bounds(2, :) = dims - 1
      grid%bounds_local = grid%bounds
      grid%bounds_local(:, 3) = [z_offset, z_offset + nz - 1]
      grid%ngpts_local = nr
      grid%dvol = weight
      DO i = 1, 3
         grid%dh(:, i) = lattice(:, i)/REAL(dims(i), dp)
      END DO
      CALL skala_grid_eval(grid, rho, gradient, tau, state%xc_section, state%particles, state%cell, &
                           energy, vrho, vgradient, vtau)
      status = 0
   END SUBROUTINE evaluate_skala_fields
#endif

! **************************************************************************************************
!> \brief Check that an atom quadrature still belongs to the current host geometry.
!> \param state Per-context callback state.
!> \return Whether the saved grid matches the current atoms and cell.
! **************************************************************************************************
   LOGICAL FUNCTION grid_geometry_current(state) RESULT(current)
      TYPE(skala_sirius_state_type), INTENT(IN)          :: state

      INTEGER                                            :: a

      current = .TRUE.
      IF (state%nradial == 0) RETURN
      current = .FALSE.
      IF (.NOT. ASSOCIATED(state%particles) .OR. .NOT. ASSOCIATED(state%cell)) RETURN
      IF (.NOT. ALLOCATED(state%grid_centres)) RETURN
      IF (SIZE(state%particles) /= SIZE(state%grid_centres, 2)) RETURN
      IF (MAXVAL(ABS(state%cell%hmat - state%grid_lattice)) > 1.0e-12_dp) RETURN
      DO a = 1, SIZE(state%particles)
         IF (MAXVAL(ABS(state%particles(a)%r - state%grid_centres(:, a))) > 1.0e-12_dp) RETURN
      END DO
      current = .TRUE.
   END FUNCTION grid_geometry_current

#if defined(__SIRIUS_XC_POINT_CALLBACK) && defined(__LIBTORCH)
! **************************************************************************************************
!> \brief Return raw discrete covectors for joint fields, without dividing by quadrature weights.
!> \param user_data Opaque callback token.
!> \param nr Number of local quadrature points.
!> \param ns Number of primitive-field spin components.
!> \param nb Number of complete local atom blocks.
!> \param fcomm Fortran context communicator handle.
!> \param coordinates Cartesian quadrature points, in bohr.
!> \param energy_weights Energy-partition quadrature weights.
!> \param descriptor_weights Nonlocal descriptor-domain weights.
!> \param atom_sizes Point count of each complete local atom block.
!> \param atom_coords Cartesian atom-block centres, in bohr.
!> \param rho Density fields.
!> \param gradient Cartesian density gradients.
!> \param tau Kinetic-energy density fields.
!> \param energy Total XC energy over the context communicator.
!> \param drho Discrete density covectors or their directional response.
!> \param dgradient Discrete gradient covectors or their directional response.
!> \param dtau Discrete kinetic-density covectors or their directional response.
!> \param status Zero on success, nonzero on callback failure.
! **************************************************************************************************
   SUBROUTINE evaluate_skala_points(user_data, nr, ns, nb, fcomm, coordinates, energy_weights, &
                                    descriptor_weights, atom_sizes, atom_coords, rho, gradient, tau, &
                                    energy, drho, dgradient, dtau, status) BIND(C)
      TYPE(C_PTR), VALUE                                 :: user_data
      INTEGER(C_INT), VALUE                              :: nr, ns, nb, fcomm
      REAL(C_DOUBLE), INTENT(IN)                         :: coordinates(3, nr), energy_weights(nr), &
                                                            descriptor_weights(nr)
      INTEGER(C_INT), INTENT(IN)                         :: atom_sizes(nb)
      REAL(C_DOUBLE), INTENT(IN)                         :: atom_coords(3, nb), rho(nr, ns), &
                                                            gradient(nr, 3, ns), tau(nr, ns)
      REAL(C_DOUBLE), INTENT(OUT)                        :: energy, drho(nr, ns), &
                                                            dgradient(nr, 3, ns), dtau(nr, ns)
      INTEGER(C_INT), INTENT(OUT)                        :: status

      INTEGER                                            :: invalid
      TYPE(mp_comm_type)                                 :: comm
      TYPE(skala_sirius_state_type), POINTER             :: state

      CALL comm%set_handle(INT(fcomm))
      energy = 0.0_dp
      drho = 0.0_dp
      dgradient = 0.0_dp
      dtau = 0.0_dp
      status = 1
      state => callback_state(user_data)
      invalid = 0
      IF (.NOT. ASSOCIATED(state)) invalid = 1
      IF (nr < 0 .OR. nb < 0 .OR. ALL(ns /= [1, 2, 4])) invalid = 1
      IF (ASSOCIATED(state)) THEN
         IF (.NOT. grid_geometry_current(state)) invalid = 1
      END IF
      CALL comm%max(invalid)
      IF (invalid /= 0) RETURN
      CALL skala_point_eval(comm, coordinates, energy_weights, descriptor_weights, atom_sizes, atom_coords, &
                            rho, gradient, tau, state%xc_section, energy, drho, dgradient, dtau, &
                            spin_angular_points=state%nspin_angular, model_fp64=.TRUE.)
      status = 0
   END SUBROUTINE evaluate_skala_points
#endif

#if defined(__LIBTORCH) && \
   (defined(__SIRIUS_XC_CALLBACK) ||defined(__SIRIUS_XC_POINT_CALLBACK) ||\
   defined(__SIRIUS_XC_POINT_GEOMETRY) ||defined(__SIRIUS_XC_POINT_RESPONSE) ||\
   defined(__SIRIUS_XC_POINT_HESSIAN) ||defined(__SIRIUS_XC_QUADRATURE_DIRECTION))
! **************************************************************************************************
!> \brief Return the registry entry without exporting Fortran pointers across the C ABI.
!> \param user_data Opaque callback token.
!> \return Matching callback state, or a null pointer for an unknown token.
! **************************************************************************************************
   FUNCTION callback_state(user_data) RESULT(state)
      TYPE(C_PTR), VALUE                                 :: user_data
      TYPE(skala_sirius_state_type), POINTER             :: state

      TYPE(callback_token_type), POINTER                 :: token

      NULLIFY (state)
      IF (.NOT. C_ASSOCIATED(user_data)) RETURN
      CALL C_F_POINTER(user_data, token)
      state => states
      DO WHILE (ASSOCIATED(state))
         IF (state%token%id == token%id) EXIT
         state => state%next
      END DO
   END FUNCTION callback_state
#endif

#if defined(__SIRIUS_XC_POINT_RESPONSE) && defined(__LIBTORCH)
! **************************************************************************************************
!> \brief Apply Skala's XC Hessian to reconstructed field directions at fixed quadrature.
!> \param user_data Opaque callback token.
!> \param nr Number of local quadrature points.
!> \param ns Number of primitive-field spin components.
!> \param nb Number of complete local atom blocks.
!> \param fcomm Fortran context communicator handle.
!> \param coordinates Cartesian quadrature points, in bohr.
!> \param energy_weights Energy-partition quadrature weights.
!> \param descriptor_weights Nonlocal descriptor-domain weights.
!> \param atom_sizes Point count of each complete local atom block.
!> \param atom_coords Cartesian atom-block centres, in bohr.
!> \param rho Density fields.
!> \param gradient Cartesian density gradients.
!> \param tau Kinetic-energy density fields.
!> \param delta_rho Density direction.
!> \param delta_gradient Density-gradient direction.
!> \param delta_tau Kinetic-density direction.
!> \param drho Discrete density covectors or their directional response.
!> \param dgradient Discrete gradient covectors or their directional response.
!> \param dtau Discrete kinetic-density covectors or their directional response.
!> \param status Zero on success, nonzero on callback failure.
! **************************************************************************************************
   SUBROUTINE evaluate_skala_point_response(user_data, nr, ns, nb, fcomm, coordinates, energy_weights, &
                                            descriptor_weights, atom_sizes, atom_coords, rho, gradient, tau, &
                                            delta_rho, delta_gradient, delta_tau, drho, dgradient, dtau, status) BIND(C)
      TYPE(C_PTR), VALUE                                 :: user_data
      INTEGER(C_INT), VALUE                              :: nr, ns, nb, fcomm
      REAL(C_DOUBLE), INTENT(IN)                         :: coordinates(3, nr), energy_weights(nr), &
                                                            descriptor_weights(nr)
      INTEGER(C_INT), INTENT(IN)                         :: atom_sizes(nb)
      REAL(C_DOUBLE), INTENT(IN) :: atom_coords(3, nb), rho(nr, ns), gradient(nr, 3, ns), &
         tau(nr, ns), delta_rho(nr, ns), delta_gradient(nr, 3, ns), delta_tau(nr, ns)
      REAL(C_DOUBLE), INTENT(OUT)                        :: drho(nr, ns), dgradient(nr, 3, ns), &
                                                            dtau(nr, ns)
      INTEGER(C_INT), INTENT(OUT)                        :: status

      INTEGER                                            :: invalid
      REAL(dp)                                           :: energy
      TYPE(mp_comm_type)                                 :: comm
      TYPE(skala_sirius_state_type), POINTER             :: state

      CALL comm%set_handle(INT(fcomm))
      state => callback_state(user_data)
      drho = 0.0_dp
      dgradient = 0.0_dp
      dtau = 0.0_dp
      status = 1
      invalid = 0
      IF (nr < 0 .OR. nb < 0 .OR. ALL(ns /= [1, 2, 4]) .OR. .NOT. ASSOCIATED(state)) invalid = 1
      IF (ASSOCIATED(state)) THEN
         IF (.NOT. grid_geometry_current(state)) invalid = 1
      END IF
      CALL comm%max(invalid)
      IF (invalid /= 0) RETURN
      CALL skala_point_eval(comm, coordinates, energy_weights, descriptor_weights, atom_sizes, atom_coords, &
                            rho, gradient, tau, state%xc_section, energy, drho, dgradient, dtau, &
                            response_density=delta_rho, response_gradient=delta_gradient, response_kinetic=delta_tau, &
                            spin_angular_points=state%nspin_angular, model_fp64=.TRUE.)
      status = 0
   END SUBROUTINE evaluate_skala_point_response
#endif

#if defined(__SIRIUS_XC_POINT_HESSIAN) && defined(__LIBTORCH)
! **************************************************************************************************
!> \brief Joint Hessian-vector product of all independent Skala point inputs.
!> \param user_data Opaque callback token.
!> \param nr Number of local quadrature points.
!> \param ns Number of primitive-field spin components.
!> \param nb Number of complete local atom blocks.
!> \param fcomm Fortran context communicator handle.
!> \param coordinates Cartesian quadrature points, in bohr.
!> \param energy_weights Energy-partition quadrature weights.
!> \param descriptor_weights Nonlocal descriptor-domain weights.
!> \param atom_sizes Point count of each complete local atom block.
!> \param atom_coords Cartesian atom-block centres, in bohr.
!> \param rho Density fields.
!> \param gradient Cartesian density gradients.
!> \param tau Kinetic-energy density fields.
!> \param delta_rho Density direction.
!> \param delta_gradient Density-gradient direction.
!> \param delta_tau Kinetic-density direction.
!> \param delta_points Quadrature-point direction.
!> \param delta_energy Energy-weight direction.
!> \param delta_descriptor Descriptor-weight direction.
!> \param delta_centres Atom-block-centre direction.
!> \param drho Discrete density covectors or their directional response.
!> \param dgradient Discrete gradient covectors or their directional response.
!> \param dtau Discrete kinetic-density covectors or their directional response.
!> \param dpoints Quadrature-point covectors or their directional response.
!> \param denergy Energy-weight covectors or their directional response.
!> \param ddescriptor Descriptor-weight covectors or their directional response.
!> \param dcentres Atom-block-centre covectors or their directional response.
!> \param status Zero on success, nonzero on callback failure.
! **************************************************************************************************
   SUBROUTINE evaluate_skala_point_hessian(user_data, nr, ns, nb, fcomm, coordinates, energy_weights, &
                                           descriptor_weights, atom_sizes, atom_coords, rho, gradient, tau, &
                                           delta_rho, delta_gradient, delta_tau, delta_points, delta_energy, &
                                           delta_descriptor, delta_centres, drho, dgradient, dtau, &
                                           dpoints, denergy, ddescriptor, dcentres, status) BIND(C)
      TYPE(C_PTR), VALUE                                 :: user_data
      INTEGER(C_INT), VALUE                              :: nr, ns, nb, fcomm
      REAL(C_DOUBLE), INTENT(IN)                         :: coordinates(3, nr), energy_weights(nr), &
                                                            descriptor_weights(nr)
      INTEGER(C_INT), INTENT(IN)                         :: atom_sizes(nb)
      REAL(C_DOUBLE), INTENT(IN) :: atom_coords(3, nb), rho(nr, ns), gradient(nr, 3, ns), &
         tau(nr, ns), delta_rho(nr, ns), delta_gradient(nr, 3, ns), delta_tau(nr, ns), &
         delta_points(3, nr), delta_energy(nr), delta_descriptor(nr), delta_centres(3, nb)
      REAL(C_DOUBLE), INTENT(OUT) :: drho(nr, ns), dgradient(nr, 3, ns), dtau(nr, ns), &
         dpoints(3, nr), denergy(nr), ddescriptor(nr), dcentres(3, nb)
      INTEGER(C_INT), INTENT(OUT)                        :: status

      INTEGER                                            :: invalid
      REAL(dp)                                           :: energy
      TYPE(mp_comm_type)                                 :: comm
      TYPE(skala_sirius_state_type), POINTER             :: state

      CALL comm%set_handle(INT(fcomm))
      state => callback_state(user_data)
      status = 1
      drho = 0.0_dp
      dgradient = 0.0_dp
      dtau = 0.0_dp
      dpoints = 0.0_dp
      denergy = 0.0_dp
      ddescriptor = 0.0_dp
      dcentres = 0.0_dp
      invalid = 0
      IF (.NOT. ASSOCIATED(state) .OR. nr < 0 .OR. nb < 0 .OR. ALL(ns /= [1, 2, 4])) invalid = 1
      IF (ASSOCIATED(state)) THEN
         IF (.NOT. grid_geometry_current(state)) invalid = 1
      END IF
      CALL comm%max(invalid)
      IF (invalid /= 0) RETURN
      CALL skala_point_eval(comm, coordinates, energy_weights, descriptor_weights, atom_sizes, atom_coords, &
                            rho, gradient, tau, state%xc_section, energy, drho, dgradient, dtau, &
                            dpoints, dcentres, denergy, ddescriptor, delta_rho, delta_gradient, delta_tau, &
                            delta_points, delta_centres, delta_energy, delta_descriptor, &
                            spin_angular_points=state%nspin_angular, model_fp64=.TRUE.)
      status = 0
   END SUBROUTINE evaluate_skala_point_hessian
#endif

#if defined(__SIRIUS_XC_POINT_GEOMETRY) && defined(__LIBTORCH)

! **************************************************************************************************
!> \brief Differentiate the model's independent fields, coordinates and quadrature weights.
!> \param user_data Opaque callback token.
!> \param nr Number of local quadrature points.
!> \param ns Number of primitive-field spin components.
!> \param nb Number of complete local atom blocks.
!> \param fcomm Fortran context communicator handle.
!> \param coordinates Cartesian quadrature points, in bohr.
!> \param energy_weights Energy-partition quadrature weights.
!> \param descriptor_weights Nonlocal descriptor-domain weights.
!> \param atom_sizes Point count of each complete local atom block.
!> \param atom_coords Cartesian atom-block centres, in bohr.
!> \param rho Density fields.
!> \param gradient Cartesian density gradients.
!> \param tau Kinetic-energy density fields.
!> \param energy Total XC energy over the context communicator.
!> \param drho Discrete density covectors or their directional response.
!> \param dgradient Discrete gradient covectors or their directional response.
!> \param dtau Discrete kinetic-density covectors or their directional response.
!> \param dpoints Quadrature-point covectors or their directional response.
!> \param denergy Energy-weight covectors or their directional response.
!> \param ddescriptor Descriptor-weight covectors or their directional response.
!> \param dcentres Atom-block-centre covectors or their directional response.
!> \param status Zero on success, nonzero on callback failure.
! **************************************************************************************************
   SUBROUTINE evaluate_skala_point_geometry(user_data, nr, ns, nb, fcomm, coordinates, energy_weights, &
                                            descriptor_weights, atom_sizes, atom_coords, rho, gradient, tau, &
                                            energy, drho, dgradient, dtau, dpoints, denergy, ddescriptor, &
                                            dcentres, status) BIND(C)
      TYPE(C_PTR), VALUE                                 :: user_data
      INTEGER(C_INT), VALUE                              :: nr, ns, nb, fcomm
      REAL(C_DOUBLE), INTENT(IN)                         :: coordinates(3, nr), energy_weights(nr), &
                                                            descriptor_weights(nr)
      INTEGER(C_INT), INTENT(IN)                         :: atom_sizes(nb)
      REAL(C_DOUBLE), INTENT(IN)                         :: atom_coords(3, nb), rho(nr, ns), &
                                                            gradient(nr, 3, ns), tau(nr, ns)
      REAL(C_DOUBLE), INTENT(OUT) :: energy, drho(nr, ns), dgradient(nr, 3, ns), dtau(nr, ns), &
         dpoints(3, nr), denergy(nr), ddescriptor(nr), dcentres(3, nb)
      INTEGER(C_INT), INTENT(OUT)                        :: status

      INTEGER                                            :: invalid
      TYPE(mp_comm_type)                                 :: comm
      TYPE(skala_sirius_state_type), POINTER             :: state

      CALL comm%set_handle(INT(fcomm))
      status = 1
      state => callback_state(user_data)
      invalid = 0
      IF (.NOT. ASSOCIATED(state)) invalid = 1
      IF (nr < 0 .OR. nb < 0 .OR. ALL(ns /= [1, 2, 4])) invalid = 1
      CALL comm%max(invalid)
      IF (invalid /= 0) RETURN
      CALL skala_point_eval(comm, coordinates, energy_weights, descriptor_weights, atom_sizes, atom_coords, &
                            rho, gradient, tau, state%xc_section, energy, drho, dgradient, dtau, &
                            dpoints, dcentres, denergy, ddescriptor, &
                            spin_angular_points=state%nspin_angular, model_fp64=.TRUE.)
      status = 0
   END SUBROUTINE evaluate_skala_point_geometry

! **************************************************************************************************
!> \brief Pull back combined model and field-sampling covectors through the periodic atom grid.
!> \param user_data Opaque callback token.
!> \param nr Number of local quadrature points.
!> \param nb Number of complete local atom blocks.
!> \param natom Number of atoms in the cell.
!> \param fcomm Fortran context communicator handle.
!> \param atom_sizes Point count of each complete local atom block.
!> \param dpoints Quadrature-point covectors or their directional response.
!> \param denergy Energy-weight covectors or their directional response.
!> \param ddescriptor Descriptor-weight covectors or their directional response.
!> \param dcentres Atom-block-centre covectors or their directional response.
!> \param datoms Cartesian atom covectors.
!> \param dlattice Cell-vector covectors.
!> \param status Zero on success, nonzero on callback failure.
! **************************************************************************************************
   SUBROUTINE pullback_skala_quadrature(user_data, nr, nb, natom, fcomm, atom_sizes, dpoints, &
                                        denergy, ddescriptor, dcentres, datoms, dlattice, status) BIND(C)
      TYPE(C_PTR), VALUE                                 :: user_data
      INTEGER(C_INT), VALUE                              :: nr, nb, natom, fcomm
      INTEGER(C_INT), INTENT(IN)                         :: atom_sizes(nb)
      REAL(C_DOUBLE), INTENT(IN)                         :: dpoints(3, nr), denergy(nr), &
                                                            ddescriptor(nr), dcentres(3, nb)
      REAL(C_DOUBLE), INTENT(OUT)                        :: datoms(3, natom), dlattice(3, 3)
      INTEGER(C_INT), INTENT(OUT)                        :: status

      INTEGER                                            :: a, invalid
      REAL(dp)                                           :: centres(3, natom)
      TYPE(mp_comm_type)                                 :: comm
      TYPE(skala_sirius_state_type), POINTER             :: state

      CALL comm%set_handle(INT(fcomm))
      status = 1
      datoms = 0.0_dp
      dlattice = 0.0_dp
      state => callback_state(user_data)
      invalid = 0
      IF (.NOT. ASSOCIATED(state)) invalid = 1
      IF (nr < 0 .OR. nb < 0 .OR. natom < 1) invalid = 1
      CALL comm%max(invalid)
      IF (invalid /= 0) RETURN
      IF (SIZE(state%particles) /= natom .OR. state%nradial < 2) invalid = 1
      CALL comm%max(invalid)
      IF (invalid /= 0) RETURN
      DO a = 1, natom
         centres(:, a) = state%particles(a)%r
      END DO
      IF (.NOT. grid_geometry_current(state)) invalid = 1
      CALL comm%max(invalid)
      IF (invalid /= 0) RETURN
      CALL skala_atom_grid_pullback(comm, state%cell%hmat, centres, state%nradial, state%nangular, state%shell, &
                                    state%scale, atom_sizes, dpoints, denergy, ddescriptor, dcentres, datoms, dlattice)
      status = 0
   END SUBROUTINE pullback_skala_quadrature
#endif

#if defined(__SIRIUS_XC_QUADRATURE_DIRECTION) && defined(__LIBTORCH)
! **************************************************************************************************
!> \brief Forward the host atom/cell direction through the complete local atom quadrature.
!> \param user_data Opaque callback token.
!> \param nr Number of local quadrature points.
!> \param nb Number of complete local atom blocks.
!> \param natom Number of atoms in the cell.
!> \param fcomm Fortran context communicator handle.
!> \param atom_sizes Point count of each complete local atom block.
!> \param datoms Cartesian atom displacement directions, in bohr.
!> \param dlattice Cell-vector directions, in bohr.
!> \param dpoints Quadrature-point covectors or their directional response.
!> \param denergy Energy-weight covectors or their directional response.
!> \param ddescriptor Descriptor-weight covectors or their directional response.
!> \param dcentres Atom-block-centre covectors or their directional response.
!> \param status Zero on success, nonzero on callback failure.
! **************************************************************************************************
   SUBROUTINE differentiate_skala_quadrature(user_data, nr, nb, natom, fcomm, atom_sizes, datoms, dlattice, &
                                             dpoints, denergy, ddescriptor, dcentres, status) BIND(C)
      TYPE(C_PTR), VALUE                                 :: user_data
      INTEGER(C_INT), VALUE                              :: nr, nb, natom, fcomm
      INTEGER(C_INT), INTENT(IN)                         :: atom_sizes(nb)
      REAL(C_DOUBLE), INTENT(IN)                         :: datoms(3, natom), dlattice(3, 3)
      REAL(C_DOUBLE), INTENT(OUT)                        :: dpoints(3, nr), denergy(nr), &
                                                            ddescriptor(nr), dcentres(3, nb)
      INTEGER(C_INT), INTENT(OUT)                        :: status

      INTEGER                                            :: a, invalid
      REAL(dp)                                           :: centres(3, natom)
      TYPE(mp_comm_type)                                 :: comm
      TYPE(skala_sirius_state_type), POINTER             :: state

      CALL comm%set_handle(INT(fcomm))
      status = 1
      dpoints = 0.0_dp
      denergy = 0.0_dp
      ddescriptor = 0.0_dp
      dcentres = 0.0_dp
      state => callback_state(user_data)
      invalid = 0
      IF (.NOT. ASSOCIATED(state) .OR. nr < 0 .OR. nb < 0 .OR. natom < 1) invalid = 1
      CALL comm%max(invalid)
      IF (invalid /= 0) RETURN
      IF (SIZE(state%particles) /= natom .OR. state%nradial < 2) invalid = 1
      IF (.NOT. grid_geometry_current(state)) invalid = 1
      CALL comm%max(invalid)
      IF (invalid /= 0) RETURN
      DO a = 1, natom
         centres(:, a) = state%particles(a)%r
      END DO
      CALL skala_atom_grid_direction(comm, state%cell%hmat, centres, state%nradial, state%nangular, state%shell, &
                                     state%scale, atom_sizes, datoms, dlattice, dpoints, denergy, ddescriptor, dcentres)
      status = 0
   END SUBROUTINE differentiate_skala_quadrature
#endif
END MODULE skala_sirius
