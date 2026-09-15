!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_twopt3d_analysis
   USE trajana_alignment,              ONLY: center_positions,&
                                              fit_rotation,&
                                              identity3
   USE trajana_cell_source,            ONLY: cell_source_type
   USE trajana_command_line,            ONLY: fail,&
                                              get_integer_option,&
                                              get_option,&
                                              get_real_option,&
                                              has_flag
   USE trajana_cube_io,                 ONLY: write_cube
   USE trajana_fft,                     ONLY: fft_any_in_place
   USE trajana_frame_controls,          ONLY: frame_selected
   USE trajana_geometry,                ONLY: minimum_image,&
                                              wrap_cartesian
   USE trajana_groups,                  ONLY: group_type,&
                                              mass_table_type,&
                                              read_groups,&
                                              read_masses,&
                                              validate_groups
   USE trajana_kinds,                   ONLY: dp
   USE trajana_molecular_motion,        ONLY: build_atom_masses,&
                                              check_homogeneous_groups,&
                                              molecule_motion
   USE trajana_selections,              ONLY: build_selection
   USE trajana_text_utils,              ONLY: lower_case
   USE trajana_trajectory_io,           ONLY: close_output,&
                                              open_output,&
                                              xyz_reader_type
   USE trajana_trajectory_types,        ONLY: frame_type
   USE trajana_twopt_thermodynamics,    ONLY: build_twopt_weights,&
                                              integrate_twopt_channel,&
                                              partition_twopt,&
                                              rotational_temperatures,&
                                              twopt_channel_type,&
                                              twopt_weights_type

   IMPLICIT NONE
   PRIVATE

   REAL(dp), PARAMETER :: bohr_per_atomic_time_to_angstrom_per_ps = 21876.91263641118_dp
   REAL(dp), PARAMETER :: light_speed_cm_s = 2.99792458E10_dp

   PUBLIC :: print_twopt3d_help, run_twopt3d

CONTAINS

   SUBROUTINE run_twopt3d()
      CHARACTER(LEN=512)                                 :: message
      CHARACTER(LEN=:), ALLOCATABLE                      :: align_spec, cell_path, cell_text, entropy_convention, &
                                                            bulk_entropy_text, grid_text, group_path, mass_path, origin_text, &
                                                            output_prefix, periodic_text, position_path, &
                                                            reference_path, spectrum_path, summary_path, &
                                                            vacf_path, velocity_path, velocity_unit, window
      INTEGER                                            :: capacity, correlation_frames, first, frames, group, &
                                                            ierr, last, minimum_origins, molecules, origins, &
                                                            output_unit, reference_atoms, rotational_dof, &
                                                            spectrum_unit, stride, vacf_unit
      INTEGER                                            :: dimensions(3)
      INTEGER, ALLOCATABLE                               :: align_index(:), counts(:), offsets(:), &
                                                            origin_voxel(:), sample_order(:)
      LOGICAL                                            :: align_found, bulk_entropy_found, current_linear, eof, origin_found, &
                                                            position_eof, reference_eof, reference_found, &
                                                            remove_drift, spectrum_requested, vacf_requested
      LOGICAL, ALLOCATABLE                               :: align_mask(:), linear(:)
      REAL(dp)                                           :: align_center(3), bulk_entropy, corrected_velocity_scale, dt_fs, &
                                                            effective_dt, excess_entropy_integral, molecular_mass, &
                                                            minus_t_delta_s_integral, rotational_symmetry, spacing, &
                                                            temperature, velocity_scale, velocity_scale_extra
      REAL(dp)                                           :: grid_origin(3), inertia_average(3), rotation_matrix(3, 3)
      REAL(dp), ALLOCATABLE                              :: align_masses(:), current_alignment(:, :), &
                                                            inertia_sum(:, :), masses(:), reference_alignment(:, :), &
                                                            rotation(:, :, :), translation(:, :, :), centers(:, :, :), &
                                                            previous_axes(:, :, :)
      REAL(dp)                                           :: current_axes(3, 3)
      TYPE(cell_source_type)                             :: cells
      TYPE(frame_type)                                   :: position_frame, reference_frame, velocity_frame
      TYPE(group_type), ALLOCATABLE                      :: groups(:)
      TYPE(mass_table_type)                              :: mass_table
      TYPE(xyz_reader_type)                              :: positions, reference, velocities

      CALL get_option("--velocity", velocity_path, eof)
      CALL get_option("--position", position_path, eof)
      CALL get_option("--groups", group_path, eof)
      CALL get_option("--mass-file", mass_path, eof)
      CALL get_option("--grid", grid_text, eof)
      CALL get_option("--origin", origin_text, origin_found)
      CALL get_option("--align-select", align_spec, align_found)
      CALL get_option("--reference", reference_path, reference_found)
      CALL get_option("--output-prefix", output_prefix, eof, "twopt3d")
      CALL get_option("--summary", summary_path, eof, "")
      CALL get_option("--spectrum", spectrum_path, spectrum_requested)
      CALL get_option("--vacf", vacf_path, vacf_requested)
      CALL get_option("--velocity-unit", velocity_unit, eof, "cp2k")
      CALL get_option("--window", window, eof, "none")
      CALL get_option("--entropy-convention", entropy_convention, eof, "lin2003")
      CALL get_option("--bulk-entropy", bulk_entropy_text, bulk_entropy_found)
      CALL get_option("--cell", cell_text, eof, "")
      CALL get_option("--cell-file", cell_path, eof, "")
      CALL get_option("--periodic", periodic_text, eof, "XYZ")
      CALL get_real_option("--dt-fs", dt_fs, -1.0_dp)
      CALL get_real_option("--temperature", temperature, -1.0_dp)
      CALL get_real_option("--spacing", spacing, -1.0_dp)
      CALL get_real_option("--velocity-scale", velocity_scale_extra, 1.0_dp)
      CALL get_real_option("--rotational-symmetry", rotational_symmetry, 1.0_dp)
      CALL get_integer_option("--correlation-frames", correlation_frames, 500)
      CALL get_integer_option("--minimum-origins", minimum_origins, 20)
      CALL get_integer_option("--first", first, 1)
      CALL get_integer_option("--last", last, HUGE(1))
      CALL get_integer_option("--stride", stride, 1)
      remove_drift = .NOT. has_flag("--keep-system-drift")

      IF (LEN_TRIM(velocity_path) == 0) CALL fail("twopt3d.x requires --velocity")
      IF (LEN_TRIM(position_path) == 0) CALL fail("twopt3d.x requires --position")
      IF (LEN_TRIM(group_path) == 0) CALL fail("twopt3d.x requires --groups")
      IF (LEN_TRIM(mass_path) == 0) CALL fail("twopt3d.x requires --mass-file")
      IF (LEN_TRIM(grid_text) == 0) CALL fail("twopt3d.x requires --grid NX NY NZ")
      READ (grid_text, *, IOSTAT=ierr) dimensions
      IF (ierr /= 0 .OR. ANY(dimensions < 1)) CALL fail("--grid expects three positive integers")
      IF (origin_found) THEN
         READ (origin_text, *, IOSTAT=ierr) grid_origin
         IF (ierr /= 0) CALL fail("--origin expects three coordinates in angstrom")
      END IF
      IF (dt_fs <= 0.0_dp) CALL fail("twopt3d.x requires a positive --dt-fs")
      IF (temperature <= 0.0_dp) CALL fail("twopt3d.x requires a positive --temperature")
      IF (spacing <= 0.0_dp) CALL fail("twopt3d.x requires a positive --spacing")
      IF (velocity_scale_extra <= 0.0_dp) CALL fail("--velocity-scale must be positive")
      IF (rotational_symmetry <= 0.0_dp) CALL fail("--rotational-symmetry must be positive")
      IF (correlation_frames < 2) CALL fail("--correlation-frames must be at least two")
      IF (minimum_origins < 1) CALL fail("--minimum-origins must be positive")
      IF (first < 1 .OR. last < first .OR. stride < 1) CALL fail("Invalid frame range")
      IF (reference_found .AND. .NOT. align_found) CALL fail("--reference requires --align-select")
      bulk_entropy = 0.0_dp
      IF (bulk_entropy_found) THEN
         READ (bulk_entropy_text, *, IOSTAT=ierr) bulk_entropy
         IF (ierr /= 0) CALL fail("--bulk-entropy expects a number in J/(mol K)")
      END IF
      IF (.NOT. origin_found) grid_origin = -0.5_dp*REAL(dimensions, dp)*spacing
      IF (LEN_TRIM(summary_path) == 0) summary_path = TRIM(output_prefix)//"-summary.dat"

      velocity_unit = lower_case(velocity_unit)
      SELECT CASE (velocity_unit)
      CASE ("cp2k", "bohr/au_time", "bohr*au_t^-1")
         velocity_scale = bohr_per_atomic_time_to_angstrom_per_ps
      CASE ("angstrom/ps", "a/ps")
         velocity_scale = 1.0_dp
      CASE ("angstrom/fs", "a/fs")
         velocity_scale = 1000.0_dp
      CASE DEFAULT
         CALL fail("--velocity-unit expects cp2k, angstrom/ps, or angstrom/fs")
      END SELECT
      corrected_velocity_scale = velocity_scale*velocity_scale_extra
      window = lower_case(window)
      IF (window /= "none" .AND. window /= "hann") CALL fail("--window expects none or hann")
      entropy_convention = lower_case(entropy_convention)
      IF (entropy_convention /= "lin2003" .AND. entropy_convention /= "rigorous") &
         CALL fail("--entropy-convention expects lin2003 or rigorous")

      CALL read_masses(mass_path, mass_table, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL read_groups(group_path, groups, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL cells%configure(cell_text, cell_path, periodic_text, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL positions%open_file(position_path, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL velocities%open_file(velocity_path, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      IF (reference_found) THEN
         CALL reference%open_file(reference_path, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         CALL reference%read_frame(reference_frame, reference_eof, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         IF (reference_eof) CALL fail("The alignment reference is empty")
         CALL reference%close_file()
      END IF

      frames = 0
      capacity = 0
      DO
         CALL positions%read_frame(position_frame, eof, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         CALL velocities%read_frame(velocity_frame, position_eof, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         IF (eof .NEQV. position_eof) CALL fail("Position and velocity trajectories have different lengths")
         IF (eof) EXIT
         CALL cells%attach(position_frame, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         IF (position_frame%number > last) EXIT
         IF (.NOT. frame_selected(position_frame%number, first, last, stride)) CYCLE

         CALL validate_synchronized_frames(position_frame, velocity_frame)
         IF (frames == 0) THEN
            CALL validate_groups(groups, position_frame%natoms, 2, ierr, message)
            IF (ierr /= 0) CALL fail(message)
            CALL build_atom_masses(position_frame, mass_table, masses, ierr, message)
            IF (ierr /= 0) CALL fail(message)
            CALL check_homogeneous_groups(groups, position_frame%label, masses, molecular_mass, ierr, message)
            IF (ierr /= 0) CALL fail(message)
            molecules = SIZE(groups)
            ALLOCATE (linear(molecules), inertia_sum(3, molecules))
            linear = .FALSE.
            inertia_sum = 0.0_dp
            ALLOCATE (previous_axes(3, 3, molecules))
            previous_axes = 0.0_dp
            IF (align_found) THEN
               CALL build_selection(align_spec, position_frame%label, align_mask, ierr, message)
               IF (ierr /= 0) CALL fail(message)
               IF (COUNT(align_mask) < 3) &
                  CALL fail("--align-select requires at least three non-collinear atoms")
               ALLOCATE (align_index(COUNT(align_mask)), align_masses(COUNT(align_mask)))
               align_index = PACK([(group, group=1, position_frame%natoms)], align_mask)
               align_masses = masses(align_index)
               ALLOCATE (reference_alignment(3, SIZE(align_index)), current_alignment(3, SIZE(align_index)))
               IF (reference_found) THEN
                  CALL select_reference(reference_frame, position_frame, align_index, reference_alignment)
               ELSE
                  reference_alignment = position_frame%value(:, align_index)
               END IF
               CALL center_positions(reference_alignment, align_masses)
               IF (.NOT. alignment_is_noncollinear(reference_alignment)) &
                  CALL fail("--align-select atoms must not all be collinear")
            ELSE
               ALLOCATE (align_index(0), align_masses(0), reference_alignment(0, 0), current_alignment(0, 0))
            END IF
            capacity = 128
            ALLOCATE (centers(3, molecules, capacity), translation(3, molecules, capacity), &
                      rotation(3, molecules, capacity))
         ELSE IF (position_frame%natoms /= SIZE(masses)) THEN
            CALL fail("The atom count changes along the trajectory")
         END IF

         IF (frames == capacity) THEN
            CALL grow_motion(centers, 2*capacity)
            CALL grow_motion(translation, 2*capacity)
            CALL grow_motion(rotation, 2*capacity)
            capacity = 2*capacity
         END IF
         frames = frames + 1
         velocity_frame%value = velocity_frame%value*corrected_velocity_scale
         IF (remove_drift) CALL remove_system_drift(velocity_frame%value, masses)

         rotation_matrix = identity3()
         align_center = 0.0_dp
         IF (align_found) THEN
            current_alignment = position_frame%value(:, align_index)
            CALL center_positions(current_alignment, align_masses, align_center)
            CALL fit_rotation(current_alignment, reference_alignment, align_masses, rotation_matrix, ierr)
            IF (ierr /= 0) CALL fail("Rotational alignment did not converge")
         END IF
         DO group = 1, molecules
            CALL molecule_motion(position_frame, velocity_frame%value, masses, groups(group), &
                                 centers(:, group, frames), translation(:, group, frames), &
                                 rotation(:, group, frames), inertia_average, current_axes, current_linear, &
                                 ierr, message)
            IF (ierr /= 0) CALL fail(message)
            IF (.NOT. current_linear) &
               CALL track_principal_axes(current_axes, previous_axes(:, :, group), rotation(:, group, frames), &
                                         frames == 1)
            IF (frames == 1) THEN
               linear(group) = current_linear
            ELSE IF (linear(group) .NEQV. current_linear) THEN
               CALL fail("A solvent molecule changes between linear and nonlinear geometry")
            END IF
            inertia_sum(:, group) = inertia_sum(:, group) + inertia_average
            CALL transform_molecule_frame(position_frame, align_found, align_center, rotation_matrix, &
                                          centers(:, group, frames), translation(:, group, frames), &
                                          rotation(:, group, frames), current_linear, ierr, message)
            IF (ierr /= 0) CALL fail(message)
         END DO
      END DO
      CALL positions%close_file()
      CALL velocities%close_file()
      CALL cells%close()
      IF (frames < correlation_frames) CALL fail("3D-2PT needs at least --correlation-frames selected frames")
      CALL shrink_motion(centers, frames)
      CALL shrink_motion(translation, frames)
      CALL shrink_motion(rotation, frames)

      IF (ANY(linear .NEQV. linear(1))) CALL fail("Analyze linear and nonlinear solvent species separately")
      IF (linear(1)) THEN
         rotational_dof = 2
      ELSE
         rotational_dof = 3
      END IF
      inertia_average = SUM(inertia_sum, DIM=2)/REAL(frames*molecules, dp)
      effective_dt = dt_fs*REAL(stride, dp)
      origins = frames - correlation_frames + 1
      CALL assign_grid_origins(centers, origins, grid_origin, spacing, dimensions, origin_voxel, counts)
      CALL order_grid_samples(origin_voxel, counts, offsets, sample_order)

      IF (spectrum_requested) THEN
         CALL open_output(spectrum_path, spectrum_unit, ierr, message)
         IF (ierr /= 0) CALL fail(message)
      ELSE
         spectrum_unit = -1
      END IF
      IF (vacf_requested) THEN
         CALL open_output(vacf_path, vacf_unit, ierr, message)
         IF (ierr /= 0) CALL fail(message)
      ELSE
         vacf_unit = -1
      END IF
      CALL analyze_grid(translation, rotation, counts, offsets, sample_order, origins, &
                        correlation_frames, effective_dt, temperature, molecular_mass, rotational_dof, &
                        inertia_average, rotational_symmetry, entropy_convention, window, minimum_origins, &
                        grid_origin, spacing, dimensions, output_prefix, spectrum_unit, vacf_unit, &
                        bulk_entropy_found, bulk_entropy, excess_entropy_integral, minus_t_delta_s_integral, &
                        ierr, message)
      IF (ierr /= 0) CALL fail(message)
      IF (spectrum_requested) CALL close_output(spectrum_unit)
      IF (vacf_requested) CALL close_output(vacf_unit)

      CALL open_output(summary_path, output_unit, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      reference_atoms = SIZE(align_index)
      CALL write_summary(output_unit, frames, origins, molecules, correlation_frames, effective_dt, temperature, &
                         dimensions, grid_origin, spacing, minimum_origins, molecular_mass, rotational_dof, &
                         inertia_average, rotational_symmetry, reference_atoms, remove_drift, entropy_convention, &
                         window, COUNT(origin_voxel == 0), output_prefix, bulk_entropy_found, bulk_entropy, &
                         excess_entropy_integral, minus_t_delta_s_integral)
      CALL close_output(output_unit)
   END SUBROUTINE run_twopt3d

   SUBROUTINE validate_synchronized_frames(position, velocity)
      TYPE(frame_type), INTENT(IN)                        :: position, velocity

      IF (position%natoms /= velocity%natoms) &
         CALL fail("Position and velocity trajectories have different atom counts")
      IF (ANY(position%label /= velocity%label)) &
         CALL fail("Position and velocity trajectories have different atom labels or ordering")
   END SUBROUTINE validate_synchronized_frames

   SUBROUTINE select_reference(reference, trajectory, selected_index, positions)
      TYPE(frame_type), INTENT(IN)                        :: reference, trajectory
      INTEGER, INTENT(IN)                                :: selected_index(:)
      REAL(dp), INTENT(OUT)                              :: positions(:, :)

      IF (reference%natoms == trajectory%natoms) THEN
         IF (ANY(reference%label(selected_index) /= trajectory%label(selected_index))) &
            CALL fail("The reference atom labels do not match the alignment selection")
         positions = reference%value(:, selected_index)
      ELSE IF (reference%natoms == SIZE(selected_index)) THEN
         IF (ANY(reference%label /= trajectory%label(selected_index))) &
            CALL fail("The reference atom labels do not match the alignment selection")
         positions = reference%value
      ELSE
         CALL fail("The reference must contain either all trajectory atoms or exactly the alignment selection")
      END IF
   END SUBROUTINE select_reference

   LOGICAL FUNCTION alignment_is_noncollinear(positions)
      REAL(dp), INTENT(IN)                               :: positions(:, :)

      INTEGER                                            :: first, second
      REAL(dp)                                           :: cross_value(3), scale

      alignment_is_noncollinear = .FALSE.
      scale = MAXVAL(SQRT(SUM(positions**2, DIM=1)))
      IF (scale <= TINY(1.0_dp)) RETURN
      DO first = 1, SIZE(positions, 2) - 1
         DO second = first + 1, SIZE(positions, 2)
            cross_value = [positions(2, first)*positions(3, second) - &
                           positions(3, first)*positions(2, second), &
                           positions(3, first)*positions(1, second) - &
                           positions(1, first)*positions(3, second), &
                           positions(1, first)*positions(2, second) - &
                           positions(2, first)*positions(1, second)]
            IF (NORM2(cross_value) > 1.0E-10_dp*scale**2) THEN
               alignment_is_noncollinear = .TRUE.
               RETURN
            END IF
         END DO
      END DO
   END FUNCTION alignment_is_noncollinear

   SUBROUTINE remove_system_drift(velocity, masses)
      REAL(dp), INTENT(INOUT)                            :: velocity(:, :)
      REAL(dp), INTENT(IN)                               :: masses(:)

      INTEGER                                            :: atom
      REAL(dp)                                           :: drift(3)

      drift = 0.0_dp
      DO atom = 1, SIZE(masses)
         drift = drift + masses(atom)*velocity(:, atom)
      END DO
      drift = drift/SUM(masses)
      DO atom = 1, SIZE(masses)
         velocity(:, atom) = velocity(:, atom) - drift
      END DO
   END SUBROUTINE remove_system_drift

   SUBROUTINE transform_molecule_frame(frame, aligned, align_center, rotation_matrix, center, translation, &
                                       rotation, linear, ierr, message)
      TYPE(frame_type), INTENT(IN)                        :: frame
      LOGICAL, INTENT(IN)                                :: aligned
      REAL(dp), INTENT(IN)                               :: align_center(3), rotation_matrix(3, 3)
      REAL(dp), INTENT(INOUT)                            :: center(3), translation(3), rotation(3)
      LOGICAL, INTENT(IN)                                :: linear
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      LOGICAL                                             :: ok
      REAL(dp)                                           :: displacement(3), transformed(3)

      ierr = 0
      message = ""
      IF (aligned) THEN
         displacement = center - align_center
         IF (frame%cell%valid) THEN
            CALL minimum_image(frame%cell, displacement, transformed, ok)
            IF (.NOT. ok) THEN
               ierr = 1
               message = "Singular cell while centering a solvent molecule"
               RETURN
            END IF
            displacement = transformed
         END IF
         center = MATMUL(rotation_matrix, displacement)
         translation = MATMUL(rotation_matrix, translation)
         IF (linear) rotation = MATMUL(rotation_matrix, rotation)
      ELSE IF (frame%cell%valid) THEN
         CALL wrap_cartesian(frame%cell, center, transformed, ok)
         IF (.NOT. ok) THEN
            ierr = 1
            message = "Singular cell while wrapping a solvent molecule"
            RETURN
         END IF
         center = transformed
      END IF
   END SUBROUTINE transform_molecule_frame

   SUBROUTINE track_principal_axes(current, previous, components, first_frame)
      REAL(dp), INTENT(INOUT)                            :: current(3, 3), previous(3, 3), components(3)
      LOGICAL, INTENT(IN)                                :: first_frame

      INTEGER                                            :: axis

      IF (.NOT. first_frame) THEN
         DO axis = 1, 3
            IF (DOT_PRODUCT(current(:, axis), previous(:, axis)) < 0.0_dp) THEN
               current(:, axis) = -current(:, axis)
               components(axis) = -components(axis)
            END IF
         END DO
      END IF
      previous = current
   END SUBROUTINE track_principal_axes

   SUBROUTINE grow_motion(values, new_capacity)
      REAL(dp), ALLOCATABLE, INTENT(INOUT)               :: values(:, :, :)
      INTEGER, INTENT(IN)                                :: new_capacity

      REAL(dp), ALLOCATABLE                              :: work(:, :, :)

      ALLOCATE (work(SIZE(values, 1), SIZE(values, 2), new_capacity))
      work(:, :, :SIZE(values, 3)) = values
      CALL MOVE_ALLOC(work, values)
   END SUBROUTINE grow_motion

   SUBROUTINE shrink_motion(values, frames)
      REAL(dp), ALLOCATABLE, INTENT(INOUT)               :: values(:, :, :)
      INTEGER, INTENT(IN)                                :: frames

      REAL(dp), ALLOCATABLE                              :: work(:, :, :)

      ALLOCATE (work(SIZE(values, 1), SIZE(values, 2), frames))
      work = values(:, :, :frames)
      CALL MOVE_ALLOC(work, values)
   END SUBROUTINE shrink_motion

   SUBROUTINE assign_grid_origins(centers, origins, origin, spacing, dimensions, voxel_of_sample, counts)
      REAL(dp), INTENT(IN)                               :: centers(:, :, :), origin(3), spacing
      INTEGER, INTENT(IN)                                :: origins, dimensions(3)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: voxel_of_sample(:), counts(:)

      INTEGER                                            :: grid_position(3), molecule, sample, time, voxel

      ALLOCATE (voxel_of_sample(origins*SIZE(centers, 2)), counts(PRODUCT(dimensions)))
      voxel_of_sample = 0
      counts = 0
      DO molecule = 1, SIZE(centers, 2)
         DO time = 1, origins
            sample = (molecule - 1)*origins + time
            grid_position = FLOOR((centers(:, molecule, time) - origin)/spacing) + 1
            IF (ANY(grid_position < 1) .OR. ANY(grid_position > dimensions)) CYCLE
            voxel = (grid_position(1) - 1)*dimensions(2)*dimensions(3) + &
                    (grid_position(2) - 1)*dimensions(3) + grid_position(3)
            voxel_of_sample(sample) = voxel
            counts(voxel) = counts(voxel) + 1
         END DO
      END DO
   END SUBROUTINE assign_grid_origins

   SUBROUTINE order_grid_samples(voxel_of_sample, counts, offsets, sample_order)
      INTEGER, INTENT(IN)                                :: voxel_of_sample(:), counts(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: offsets(:), sample_order(:)

      INTEGER                                            :: sample, voxel
      INTEGER, ALLOCATABLE                               :: cursor(:)

      ALLOCATE (offsets(0:SIZE(counts)), cursor(SIZE(counts)), sample_order(SUM(counts)))
      offsets(0) = 0
      DO voxel = 1, SIZE(counts)
         offsets(voxel) = offsets(voxel - 1) + counts(voxel)
      END DO
      cursor = offsets(:SIZE(counts) - 1)
      DO sample = 1, SIZE(voxel_of_sample)
         voxel = voxel_of_sample(sample)
         IF (voxel == 0) CYCLE
         cursor(voxel) = cursor(voxel) + 1
         sample_order(cursor(voxel)) = sample
      END DO
   END SUBROUTINE order_grid_samples

   SUBROUTINE analyze_grid(translation, rotation, counts, offsets, sample_order, origins, &
                           correlation_frames, dt_fs, temperature, molecular_mass, rotational_dof, &
                           inertia, rotational_symmetry, entropy_convention, window, minimum_origins, &
                           origin, spacing, dimensions, prefix, spectrum_unit, vacf_unit, &
                           bulk_entropy_found, bulk_entropy, excess_entropy_integral, &
                           minus_t_delta_s_integral, ierr, message)
      REAL(dp), INTENT(IN)                               :: translation(:, :, :), rotation(:, :, :), dt_fs, &
                                                            temperature, molecular_mass, inertia(3), &
                                                            rotational_symmetry, origin(3), spacing
      INTEGER, INTENT(IN)                                :: counts(:), offsets(0:), &
                                                            sample_order(:), origins, correlation_frames, &
                                                            rotational_dof, minimum_origins, dimensions(3), &
                                                            spectrum_unit, vacf_unit
      CHARACTER(LEN=*), INTENT(IN)                       :: entropy_convention, window, prefix
      LOGICAL, INTENT(IN)                                :: bulk_entropy_found
      REAL(dp), INTENT(IN)                               :: bulk_entropy
      REAL(dp), INTENT(OUT)                              :: excess_entropy_integral, minus_t_delta_s_integral
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      CHARACTER(LEN=1024)                                :: path
      INTEGER                                            :: entry, frequency, lag, molecule, sample, time, voxel
      REAL(dp)                                           :: density_value, df_cm, local_volume, negative_rotation, &
                                                            negative_translation, rotational_temperature(3), &
                                                            voxel_center(3), voxel_volume
      REAL(dp), ALLOCATABLE                              :: density(:), frequency_cm(:), origin_count(:), &
                                                            entropy_excess(:), minus_t_delta_s(:), &
                                                            minus_t_delta_s_density(:), &
                                                            rotation_corr(:), rotation_dos(:), rotation_gas(:), &
                                                            rotation_solid(:), rotation_fluidicity(:), &
                                                            rotation_entropy(:), rotation_negative(:), total_entropy(:), &
                                                            translation_corr(:), translation_dos(:), &
                                                            translation_gas(:), translation_solid(:), &
                                                            translation_fluidicity(:), translation_entropy(:), &
                                                            translation_negative(:)
      TYPE(twopt_channel_type)                           :: rotation_result, translation_result
      TYPE(twopt_weights_type)                           :: weights

      ierr = 0
      message = ""
      excess_entropy_integral = 0.0_dp
      minus_t_delta_s_integral = 0.0_dp
      voxel_volume = spacing**3
      ALLOCATE (density(SIZE(counts)), origin_count(SIZE(counts)), &
                translation_entropy(SIZE(counts)), rotation_entropy(SIZE(counts)), total_entropy(SIZE(counts)), &
                entropy_excess(SIZE(counts)), minus_t_delta_s(SIZE(counts)), &
                minus_t_delta_s_density(SIZE(counts)), &
                translation_fluidicity(SIZE(counts)), rotation_fluidicity(SIZE(counts)), &
                translation_negative(SIZE(counts)), rotation_negative(SIZE(counts)))
      origin_count = REAL(counts, dp)
      density = origin_count/(REAL(origins, dp)*voxel_volume)
      translation_entropy = 0.0_dp
      rotation_entropy = 0.0_dp
      total_entropy = 0.0_dp
      entropy_excess = 0.0_dp
      minus_t_delta_s = 0.0_dp
      minus_t_delta_s_density = 0.0_dp
      translation_fluidicity = 0.0_dp
      rotation_fluidicity = 0.0_dp
      translation_negative = 0.0_dp
      rotation_negative = 0.0_dp
      CALL rotational_temperatures(inertia, rotational_dof == 2, rotational_temperature)

      IF (spectrum_unit /= -1) THEN
         WRITE (spectrum_unit, "(A)") "# ix iy iz x[A] y[A] z[A] origins density[molecule/A^3] "// &
                                      "wavenumber[cm^-1] trans_DoS trans_gas trans_solid rot_DoS rot_gas rot_solid"
      END IF
      IF (vacf_unit /= -1) THEN
         WRITE (vacf_unit, "(A)") "# ix iy iz x[A] y[A] z[A] origins lag_time[fs] trans_VACF rot_VACF"
      END IF

      DO voxel = 1, SIZE(counts)
         IF (counts(voxel) < minimum_origins) CYCLE
         ALLOCATE (translation_corr(0:correlation_frames - 1), rotation_corr(0:correlation_frames - 1))
         translation_corr = 0.0_dp
         rotation_corr = 0.0_dp
         DO entry = offsets(voxel - 1) + 1, offsets(voxel)
            sample = sample_order(entry)
            molecule = (sample - 1)/origins + 1
            time = MOD(sample - 1, origins) + 1
            DO lag = 0, correlation_frames - 1
               translation_corr(lag) = translation_corr(lag) + &
                  DOT_PRODUCT(translation(:, molecule, time), translation(:, molecule, time + lag))
               rotation_corr(lag) = rotation_corr(lag) + &
                  DOT_PRODUCT(rotation(:, molecule, time), rotation(:, molecule, time + lag))
            END DO
         END DO
         translation_corr = translation_corr/REAL(counts(voxel), dp)
         rotation_corr = rotation_corr/REAL(counts(voxel), dp)
         CALL correlation_to_dos(translation_corr, 3, dt_fs, window, translation_dos, df_cm, &
                                 negative_translation, ierr, message)
         IF (ierr /= 0) RETURN
         CALL correlation_to_dos(rotation_corr, rotational_dof, dt_fs, window, rotation_dos, df_cm, &
                                 negative_rotation, ierr, message)
         IF (ierr /= 0) RETURN
         translation_negative(voxel) = negative_translation
         rotation_negative(voxel) = negative_rotation
         ALLOCATE (frequency_cm(0:UBOUND(translation_dos, 1)))
         frequency_cm = [(REAL(frequency, dp)*df_cm, frequency=0, UBOUND(translation_dos, 1))]

         density_value = density(voxel)
         local_volume = 1.0_dp/density_value
         CALL partition_twopt(translation_dos, df_cm, temperature, molecular_mass, 1, local_volume, &
                              translation_gas, translation_solid, translation_result)
         CALL partition_twopt(rotation_dos, df_cm, temperature, molecular_mass, 1, local_volume, &
                              rotation_gas, rotation_solid, rotation_result)
         CALL build_twopt_weights(translation_result%packing_fraction, molecular_mass, &
            translation_result%gas_dof/3.0_dp, temperature, local_volume, rotational_temperature, &
            rotational_symmetry, entropy_convention, weights)
         CALL integrate_twopt_channel(frequency_cm, df_cm, temperature, translation_gas, translation_solid, &
                                      weights%entropy_translation, translation_result)
         CALL integrate_twopt_channel(frequency_cm, df_cm, temperature, rotation_gas, rotation_solid, &
                                      weights%entropy_rotation, rotation_result)
         translation_entropy(voxel) = translation_result%entropy_quantum_j_mol_k
         rotation_entropy(voxel) = rotation_result%entropy_quantum_j_mol_k
         total_entropy(voxel) = translation_entropy(voxel) + rotation_entropy(voxel)
         translation_fluidicity(voxel) = translation_result%fluidicity
         rotation_fluidicity(voxel) = rotation_result%fluidicity
         CALL voxel_coordinates(voxel, dimensions, origin, spacing, voxel_center)

         IF (vacf_unit /= -1) THEN
            DO lag = 0, correlation_frames - 1
               WRITE (vacf_unit, "(3I8,3ES24.16,I12,3ES24.16)") &
                  voxel_indices(voxel, dimensions), voxel_center, counts(voxel), &
                  REAL(lag, dp)*dt_fs, translation_corr(lag)/translation_corr(0), &
                  rotation_corr(lag)/rotation_corr(0)
            END DO
         END IF
         IF (spectrum_unit /= -1) THEN
            DO frequency = 0, UBOUND(translation_dos, 1)
               WRITE (spectrum_unit, "(3I8,3ES24.16,I12,8ES24.16)") &
                  voxel_indices(voxel, dimensions), voxel_center, counts(voxel), density_value, &
                  frequency_cm(frequency), translation_dos(frequency), translation_gas(frequency), &
                  translation_solid(frequency), rotation_dos(frequency), rotation_gas(frequency), &
                  rotation_solid(frequency)
            END DO
         END IF
         DEALLOCATE (translation_corr, rotation_corr, translation_dos, rotation_dos, frequency_cm, &
                     translation_gas, translation_solid, rotation_gas, rotation_solid)
      END DO

      path = TRIM(prefix)//"-density.cube"
      CALL write_cube(path, "3D-2PT molecular number density [molecule/angstrom^3]", origin, spacing, &
                      dimensions, density, ierr, message)
      IF (ierr /= 0) RETURN
      path = TRIM(prefix)//"-origins.cube"
      CALL write_cube(path, "3D-2PT time-origin sample count", origin, spacing, dimensions, origin_count, &
                      ierr, message)
      IF (ierr /= 0) RETURN
      path = TRIM(prefix)//"-entropy-translation.cube"
      CALL write_cube(path, "3D-2PT translational entropy per solvent molecule [J/(mol K)]", origin, spacing, &
                      dimensions, translation_entropy, ierr, message)
      IF (ierr /= 0) RETURN
      path = TRIM(prefix)//"-entropy-rotation.cube"
      CALL write_cube(path, "3D-2PT rotational entropy per solvent molecule [J/(mol K)]", origin, spacing, &
                      dimensions, rotation_entropy, ierr, message)
      IF (ierr /= 0) RETURN
      path = TRIM(prefix)//"-entropy-total.cube"
      CALL write_cube(path, "3D-2PT intermolecular entropy per solvent molecule [J/(mol K)]", origin, spacing, &
                      dimensions, total_entropy, ierr, message)
      IF (ierr /= 0) RETURN
      path = TRIM(prefix)//"-fluidicity-translation.cube"
      CALL write_cube(path, "3D-2PT translational fluidicity", origin, spacing, dimensions, &
                      translation_fluidicity, ierr, message)
      IF (ierr /= 0) RETURN
      path = TRIM(prefix)//"-fluidicity-rotation.cube"
      CALL write_cube(path, "3D-2PT rotational fluidicity", origin, spacing, dimensions, &
                      rotation_fluidicity, ierr, message)
      IF (ierr /= 0) RETURN
      path = TRIM(prefix)//"-negative-vdos-translation.cube"
      CALL write_cube(path, "Fraction of negative translational VDoS removed before 2PT", origin, spacing, &
                      dimensions, translation_negative, ierr, message)
      IF (ierr /= 0) RETURN
      path = TRIM(prefix)//"-negative-vdos-rotation.cube"
      CALL write_cube(path, "Fraction of negative rotational VDoS removed before 2PT", origin, spacing, &
                      dimensions, rotation_negative, ierr, message)
      IF (ierr /= 0) RETURN

      IF (bulk_entropy_found) THEN
         WHERE (counts >= minimum_origins)
            entropy_excess = total_entropy - bulk_entropy
            minus_t_delta_s = -temperature*entropy_excess/1000.0_dp
            minus_t_delta_s_density = density*minus_t_delta_s
         END WHERE
         excess_entropy_integral = SUM(density*entropy_excess, MASK=counts >= minimum_origins)*voxel_volume
         minus_t_delta_s_integral = SUM(minus_t_delta_s_density)*voxel_volume

         path = TRIM(prefix)//"-entropy-excess.cube"
         CALL write_cube(path, "3D-2PT excess intermolecular entropy per solvent [J/(mol K)]", origin, spacing, &
                         dimensions, entropy_excess, ierr, message)
         IF (ierr /= 0) RETURN
         path = TRIM(prefix)//"-minus-t-delta-s.cube"
         CALL write_cube(path, "3D-2PT -T Delta S per solvent molecule [kJ/mol]", origin, spacing, &
                         dimensions, minus_t_delta_s, ierr, message)
         IF (ierr /= 0) RETURN
         path = TRIM(prefix)//"-minus-t-delta-s-density.cube"
         CALL write_cube(path, "3D-2PT density-weighted -T Delta S [kJ/(mol angstrom^3)]", origin, spacing, &
                         dimensions, minus_t_delta_s_density, ierr, message)
      END IF
   END SUBROUTINE analyze_grid

   SUBROUTINE correlation_to_dos(correlation, target_dof, dt_fs, window, dos, frequency_step, &
                                 negative_fraction, ierr, message)
      REAL(dp), INTENT(IN)                               :: correlation(0:), dt_fs
      INTEGER, INTENT(IN)                                :: target_dof
      CHARACTER(LEN=*), INTENT(IN)                       :: window
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: dos(:)
      REAL(dp), INTENT(OUT)                              :: frequency_step, negative_fraction
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      COMPLEX(dp), ALLOCATABLE                           :: work(:)
      INTEGER                                            :: frequency, lag, n, nfft
      REAL(dp)                                           :: area, denominator, taper

      ierr = 0
      message = ""
      n = SIZE(correlation)
      nfft = 2*n - 1
      ALLOCATE (work(nfft), dos(0:n - 1))
      work = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO lag = 0, n - 1
         taper = 1.0_dp
         IF (window == "hann" .AND. n > 1) &
            taper = 0.5_dp*(1.0_dp + COS(ACOS(-1.0_dp)*REAL(lag, dp)/REAL(n - 1, dp)))
         work(lag + 1) = CMPLX(correlation(lag)*taper, 0.0_dp, KIND=dp)
         IF (lag > 0) work(nfft - lag + 1) = work(lag + 1)
      END DO
      CALL fft_any_in_place(work, inverse=.FALSE.)
      DO frequency = 0, n - 1
         dos(frequency) = REAL(work(frequency + 1), dp)
      END DO
      denominator = SUM(ABS(dos))
      IF (denominator > TINY(1.0_dp)) THEN
         negative_fraction = SUM(ABS(MIN(dos, 0.0_dp)))/denominator
      ELSE
         negative_fraction = 0.0_dp
      END IF
      dos = MAX(dos, 0.0_dp)
      frequency_step = 1.0E15_dp/(REAL(nfft, dp)*dt_fs*light_speed_cm_s)
      area = frequency_step*(0.5_dp*(dos(0) + dos(UBOUND(dos, 1))) + &
                             SUM(dos(1:UBOUND(dos, 1) - 1)))
      IF (area <= TINY(1.0_dp)) THEN
         ierr = 1
         message = "A populated voxel has zero spectral power"
         RETURN
      END IF
      dos = dos*REAL(target_dof, dp)/area
   END SUBROUTINE correlation_to_dos

   PURE SUBROUTINE voxel_coordinates(voxel, dimensions, origin, spacing, center)
      INTEGER, INTENT(IN)                                :: voxel, dimensions(3)
      REAL(dp), INTENT(IN)                               :: origin(3), spacing
      REAL(dp), INTENT(OUT)                              :: center(3)

      INTEGER                                            :: indices(3)

      indices = voxel_indices(voxel, dimensions)
      center = origin + (REAL(indices, dp) - 0.5_dp)*spacing
   END SUBROUTINE voxel_coordinates

   PURE FUNCTION voxel_indices(voxel, dimensions) RESULT(indices)
      INTEGER, INTENT(IN)                                :: voxel, dimensions(3)
      INTEGER                                            :: indices(3), remainder

      indices(1) = (voxel - 1)/(dimensions(2)*dimensions(3)) + 1
      remainder = MOD(voxel - 1, dimensions(2)*dimensions(3))
      indices(2) = remainder/dimensions(3) + 1
      indices(3) = MOD(remainder, dimensions(3)) + 1
   END FUNCTION voxel_indices

   SUBROUTINE write_summary(unit, frames, origins, molecules, correlation_frames, dt_fs, temperature, &
                            dimensions, origin, spacing, minimum_origins, molecular_mass, rotational_dof, &
                            inertia, rotational_symmetry, alignment_atoms, remove_drift, entropy_convention, &
                            window, outside, prefix, bulk_entropy_found, bulk_entropy, &
                            excess_entropy_integral, minus_t_delta_s_integral)
      INTEGER, INTENT(IN)                                :: unit, frames, origins, molecules, correlation_frames, &
                                                            dimensions(3), minimum_origins, rotational_dof, &
                                                            alignment_atoms, outside
      REAL(dp), INTENT(IN)                               :: dt_fs, temperature, origin(3), spacing, molecular_mass, &
                                                            inertia(3), rotational_symmetry, bulk_entropy, &
                                                            excess_entropy_integral, minus_t_delta_s_integral
      LOGICAL, INTENT(IN)                                :: remove_drift, bulk_entropy_found
      CHARACTER(LEN=*), INTENT(IN)                       :: entropy_convention, window, prefix

      WRITE (unit, "(A)") "# Spatially resolved two-phase thermodynamics"
      WRITE (unit, "(A,I0)") "frames: ", frames
      WRITE (unit, "(A,I0)") "time_origins_per_molecule: ", origins
      WRITE (unit, "(A,I0)") "molecules: ", molecules
      WRITE (unit, "(A,I0)") "correlation_frames: ", correlation_frames
      WRITE (unit, "(A,ES24.16)") "effective_dt_fs: ", dt_fs
      WRITE (unit, "(A,ES24.16)") "temperature_K: ", temperature
      WRITE (unit, "(A,3(I0,1X))") "grid: ", dimensions
      WRITE (unit, "(A,3ES24.16)") "origin_angstrom: ", origin
      WRITE (unit, "(A,ES24.16)") "spacing_angstrom: ", spacing
      WRITE (unit, "(A,I0)") "minimum_origins: ", minimum_origins
      WRITE (unit, "(A,I0)") "outside_grid_samples: ", outside
      WRITE (unit, "(A,ES24.16)") "molecular_mass_g_mol: ", molecular_mass
      WRITE (unit, "(A,I0)") "rotational_dof_per_molecule: ", rotational_dof
      WRITE (unit, "(A,3ES24.16)") "mean_principal_inertia_amu_angstrom2: ", inertia
      WRITE (unit, "(A,ES24.16)") "rotational_symmetry: ", rotational_symmetry
      WRITE (unit, "(A,I0)") "alignment_atoms: ", alignment_atoms
      WRITE (unit, "(A,L1)") "system_drift_removed: ", remove_drift
      WRITE (unit, "(A,A)") "entropy_convention: ", TRIM(entropy_convention)
      WRITE (unit, "(A,A)") "vacf_window: ", TRIM(window)
      WRITE (unit, "(A,A)") "output_prefix: ", TRIM(prefix)
      IF (bulk_entropy_found) THEN
         WRITE (unit, "(A,ES24.16)") "bulk_intermolecular_entropy_J_mol_K: ", bulk_entropy
         WRITE (unit, "(A,ES24.16)") "integrated_excess_entropy_J_mol_K: ", excess_entropy_integral
         WRITE (unit, "(A,ES24.16)") "integrated_minus_T_delta_S_kJ_mol: ", minus_t_delta_s_integral
      ELSE
         WRITE (unit, "(A)") "bulk_intermolecular_entropy: not supplied"
      END IF
      WRITE (unit, "(A)") "note: entropy maps contain zero where origin_count < minimum_origins"
      WRITE (unit, "(A)") "note: negative finite-sampling VDoS bins are projected to zero before normalization"
      WRITE (unit, "(A)") "note: pair-energy enthalpy/free-energy maps require an external energy decomposition"
   END SUBROUTINE write_summary

   SUBROUTINE print_twopt3d_help()
      WRITE (*, "(A)") "Usage: twopt3d.x --velocity VEL.xyz --position POS.xyz --groups solvent.groups [options]"
      WRITE (*, "(A)") ""
      WRITE (*, "(A)") "Required:"
      WRITE (*, "(A)") "  --velocity FILE          synchronized CP2K XYZ velocity trajectory"
      WRITE (*, "(A)") "  --position FILE          synchronized CP2K XYZ position trajectory"
      WRITE (*, "(A)") "  --groups FILE            one homogeneous solvent molecule per line"
      WRITE (*, "(A)") "  --mass-file FILE         LABEL MASS[g/mol] table for all trajectory atoms"
      WRITE (*, "(A)") "  --dt-fs VALUE            time between input frames"
      WRITE (*, "(A)") "  --temperature VALUE      thermodynamic temperature [K]"
      WRITE (*, "(A)") "  --grid NX,NY,NZ          orthogonal grid dimensions"
      WRITE (*, "(A)") "  --spacing VALUE          isotropic voxel spacing [angstrom]"
      WRITE (*, "(A)") ""
      WRITE (*, "(A)") "Spatial frame and sampling:"
      WRITE (*, "(A)") "  --origin X,Y,Z           grid corner [angstrom] (default: grid centered at zero)"
      WRITE (*, "(A)") "  --align-select SPEC      solute atoms defining the moving reference frame"
      WRITE (*, "(A)") "  --reference FILE         optional alignment reference (default first selected frame)"
      WRITE (*, "(A)") "  --cell/--cell-file       cell for molecule reconstruction and minimum images"
      WRITE (*, "(A)") "  --correlation-frames N   local VACF points (default 500)"
      WRITE (*, "(A)") "  --minimum-origins N      samples required for an entropy voxel (default 20)"
      WRITE (*, "(A)") "  --first/--last/--stride  frame selection"
      WRITE (*, "(A)") ""
      WRITE (*, "(A)") "Thermodynamics and output:"
      WRITE (*, "(A)") "  --rotational-symmetry N  solvent rotational symmetry number (default 1)"
      WRITE (*, "(A)") "  --window none|hann       VACF lag window (default none)"
      WRITE (*, "(A)") "  --entropy-convention C   lin2003 (default) or rigorous"
      WRITE (*, "(A)") "  --bulk-entropy VALUE     bulk translation+rotation entropy [J/(mol K)]"
      WRITE (*, "(A)") "  --output-prefix PREFIX   CUBE output prefix (default twopt3d)"
      WRITE (*, "(A)") "  --summary FILE           run metadata and limitations"
      WRITE (*, "(A)") "  --spectrum FILE          optional voxel-resolved VDoS diagnostics"
      WRITE (*, "(A)") "  --vacf FILE              optional voxel-resolved VACF diagnostics"
      WRITE (*, "(A)") "  --velocity-unit UNIT     cp2k (default), angstrom/ps, or angstrom/fs"
   END SUBROUTINE print_twopt3d_help

END MODULE trajana_twopt3d_analysis
