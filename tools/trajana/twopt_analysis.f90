!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_twopt_analysis
   USE trajana_cell_source,             ONLY: cell_source_type
   USE trajana_command_line,            ONLY: fail,&
                                              get_integer_option,&
                                              get_option,&
                                              get_real_option,&
                                              has_flag
   USE trajana_frame_controls,          ONLY: frame_selected
   USE trajana_geometry,                ONLY: cell_volume,&
                                              minimum_image
   USE trajana_groups,                  ONLY: group_type,&
                                              mass_table_type,&
                                              read_groups,&
                                              read_masses,&
                                              validate_groups
   USE trajana_kinds,                   ONLY: dp
   USE trajana_text_utils,              ONLY: lower_case
   USE trajana_time_series,             ONLY: real_autocorrelation_sum,&
                                              real_power_sum
   USE trajana_trajectory_io,           ONLY: close_output,&
                                              open_output,&
                                              xyz_reader_type
   USE trajana_trajectory_types,        ONLY: frame_type
   USE trajana_twopt_thermodynamics,    ONLY: build_twopt_weights,&
                                              gas_constant,&
                                              integrate_twopt_channel,&
                                              partition_twopt,&
                                              rotational_temperatures,&
                                              twopt_channel_type,&
                                              twopt_weights_type

   IMPLICIT NONE
   PRIVATE

   REAL(dp), PARAMETER :: bohr_per_atomic_time_to_angstrom_per_ps = 21876.91263641118_dp

   PUBLIC :: run_twopt, print_twopt_help

CONTAINS

   SUBROUTINE run_twopt()
      CHARACTER(LEN=512)                                 :: message
      CHARACTER(LEN=:), ALLOCATABLE                      :: cell_path, cell_text, entropy_convention, &
                                                            group_path, mass_path, output_path, &
                                                            periodic_text, position_path, spectrum_path, &
                                                            vacf_path, velocity_path, velocity_unit, window
      INTEGER                                            :: capacity, constraints, first, frames, ierr, last, &
                                                            molecules, nfft, output_unit, rotational_dof, &
                                                            spectrum_unit, stride, vibrational_dof, volume_frames
      LOGICAL                                            :: eof, groups_found, molecular, &
                                                            position_eof, position_found, remove_drift, &
                                                            spectrum_requested, vacf_requested
      LOGICAL, ALLOCATABLE                               :: linear(:)
      REAL(dp)                                           :: classical_cv, dt_fs, effective_dt, energy, &
                                                            energy_reference, explicit_volume, mass_total, &
                                                            molecular_mass, rotational_symmetry, &
                                                            selected_volume, temperature, velocity_scale, &
                                                            velocity_scale_extra
      REAL(dp)                                           :: kinetic_temperature(3), rotational_temperature(3)
      REAL(dp), ALLOCATABLE                              :: frequency(:), inertia_sum(:, :), masses(:), &
                                                            position_values(:, :), rotational(:, :), &
                                                            rotational_dos(:), rotational_gas(:), &
                                                            rotational_solid(:), translation(:, :), &
                                                            translation_dos(:), translation_gas(:), &
                                                            translation_solid(:), vibration(:, :), &
                                                            vibration_dos(:), vibration_gas(:), &
                                                            vibration_solid(:)
      TYPE(cell_source_type)                             :: cells
      TYPE(frame_type)                                   :: position_frame, velocity_frame
      TYPE(group_type), ALLOCATABLE                      :: groups(:)
      TYPE(mass_table_type)                              :: mass_table
      TYPE(twopt_channel_type)                           :: rotation_result, total_result, &
                                                            translation_result, vibration_result
      TYPE(twopt_weights_type)                           :: weights
      TYPE(xyz_reader_type)                              :: positions, velocities

      CALL get_option("--velocity", velocity_path, eof)
      CALL get_option("--position", position_path, position_found)
      CALL get_option("--mass-file", mass_path, eof)
      CALL get_option("--groups", group_path, groups_found)
      CALL get_option("--output", output_path, eof, "-")
      CALL get_option("--spectrum", spectrum_path, spectrum_requested)
      CALL get_option("--vacf", vacf_path, vacf_requested)
      CALL get_option("--velocity-unit", velocity_unit, eof, "cp2k")
      CALL get_option("--window", window, eof, "none")
      CALL get_option("--entropy-convention", entropy_convention, eof, "lin2003")
      CALL get_option("--cell", cell_text, eof, "")
      CALL get_option("--cell-file", cell_path, eof, "")
      CALL get_option("--periodic", periodic_text, eof, "XYZ")
      CALL get_real_option("--dt-fs", dt_fs, -1.0_dp)
      CALL get_real_option("--temperature", temperature, -1.0_dp)
      CALL get_real_option("--volume", explicit_volume, -1.0_dp)
      CALL get_real_option("--velocity-scale", velocity_scale_extra, 1.0_dp)
      CALL get_real_option("--rotational-symmetry", rotational_symmetry, 1.0_dp)
      CALL get_real_option("--energy-kj-mol", energy, HUGE(1.0_dp))
      CALL get_real_option("--classical-cv-j-mol-k", classical_cv, HUGE(1.0_dp))
      CALL get_integer_option("--constraints", constraints, 0)
      CALL get_integer_option("--first", first, 1)
      CALL get_integer_option("--last", last, HUGE(1))
      CALL get_integer_option("--stride", stride, 1)
      remove_drift = .NOT. has_flag("--keep-system-drift")

      IF (LEN_TRIM(velocity_path) == 0) CALL fail("twopt.x requires --velocity")
      IF (LEN_TRIM(mass_path) == 0) CALL fail("twopt.x requires --mass-file")
      IF (dt_fs <= 0.0_dp) CALL fail("twopt.x requires a positive --dt-fs")
      IF (temperature <= 0.0_dp) CALL fail("twopt.x requires a positive --temperature")
      IF (explicit_volume >= 0.0_dp .AND. explicit_volume <= TINY(1.0_dp)) &
         CALL fail("--volume must be positive when provided")
      IF (velocity_scale_extra <= 0.0_dp) CALL fail("--velocity-scale must be positive")
      IF (rotational_symmetry <= 0.0_dp) CALL fail("--rotational-symmetry must be positive")
      IF (constraints < 0) CALL fail("--constraints cannot be negative")
      IF (first < 1 .OR. last < first .OR. stride < 1) CALL fail("Invalid frame range")

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
      velocity_scale = velocity_scale*velocity_scale_extra
      window = lower_case(window)
      IF (window /= "none" .AND. window /= "hann") CALL fail("--window expects none or hann")
      entropy_convention = lower_case(entropy_convention)
      IF (entropy_convention /= "lin2003" .AND. entropy_convention /= "rigorous") &
         CALL fail("--entropy-convention expects lin2003 or rigorous")

      CALL read_masses(mass_path, mass_table, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      IF (groups_found) THEN
         CALL read_groups(group_path, groups, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         molecular = groups_need_positions(groups)
         IF (molecular .AND. .NOT. position_found) CALL fail("Molecular 2PT requires --position")
      ELSE
         molecular = .FALSE.
         ALLOCATE (groups(0))
      END IF

      CALL cells%configure(cell_text, cell_path, periodic_text, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      IF (position_found) THEN
         CALL positions%open_file(position_path, ierr, message)
         IF (ierr /= 0) CALL fail(message)
      END IF
      CALL velocities%open_file(velocity_path, ierr, message)
      IF (ierr /= 0) CALL fail(message)

      frames = 0
      capacity = 0
      selected_volume = 0.0_dp
      volume_frames = 0
      mass_total = 0.0_dp
      molecules = 0
      rotational_dof = 0
      vibrational_dof = 0
      DO
         CALL velocities%read_frame(velocity_frame, eof, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         IF (position_found) THEN
            CALL positions%read_frame(position_frame, position_eof, ierr, message)
            IF (ierr /= 0) CALL fail(message)
            IF (eof .NEQV. position_eof) CALL fail("Position and velocity trajectories have different lengths")
            IF (.NOT. eof) THEN
               CALL cells%attach(position_frame, ierr, message)
               IF (ierr /= 0) CALL fail(message)
            END IF
         ELSE IF (.NOT. eof) THEN
            CALL cells%attach(velocity_frame, ierr, message)
            IF (ierr /= 0) CALL fail(message)
         END IF
         IF (eof) EXIT
         IF (velocity_frame%number > last) EXIT
         IF (.NOT. frame_selected(velocity_frame%number, first, last, stride)) CYCLE

         IF (frames == 0) THEN
            IF (.NOT. groups_found) CALL make_atomic_groups(velocity_frame%natoms, groups)
            CALL validate_groups(groups, velocity_frame%natoms, 1, ierr, message)
            IF (ierr /= 0) CALL fail(message)
            CALL validate_group_partition(groups, velocity_frame%natoms, ierr, message)
            IF (ierr /= 0) CALL fail(message)
            IF (position_found .AND. position_frame%natoms /= velocity_frame%natoms) &
               CALL fail("Position and velocity trajectories have different atom counts")
            CALL build_masses(velocity_frame, mass_table, masses, mass_total, ierr, message)
            IF (ierr /= 0) CALL fail(message)
            molecules = SIZE(groups)
            ALLOCATE (linear(molecules), inertia_sum(3, molecules))
            linear = .FALSE.
            inertia_sum = 0.0_dp
            capacity = 128
            ALLOCATE (translation(3*molecules, capacity))
            IF (molecular) THEN
               ALLOCATE (rotational(3*molecules, capacity), vibration(3*velocity_frame%natoms, capacity))
            ELSE
               ALLOCATE (rotational(0, capacity), vibration(0, capacity))
            END IF
         ELSE
            IF (velocity_frame%natoms /= SIZE(masses)) CALL fail("The atom count changes along the trajectory")
            IF (position_found .AND. position_frame%natoms /= SIZE(masses)) &
               CALL fail("The atom count changes along the position trajectory")
         END IF

         IF (frames == capacity) THEN
            CALL grow_series(translation, 2*capacity)
            CALL grow_series(rotational, 2*capacity)
            CALL grow_series(vibration, 2*capacity)
            capacity = 2*capacity
         END IF
         frames = frames + 1
         IF (position_found) THEN
            position_values = position_frame%value
            IF (position_frame%cell%valid) THEN
               selected_volume = selected_volume + cell_volume(position_frame%cell)
               volume_frames = volume_frames + 1
            END IF
         ELSE
            ALLOCATE (position_values(0, 0))
            IF (velocity_frame%cell%valid) THEN
               selected_volume = selected_volume + cell_volume(velocity_frame%cell)
               volume_frames = volume_frames + 1
            END IF
         END IF
         IF (position_found) THEN
            CALL decompose_frame(velocity_frame%value*velocity_scale, position_values, position_frame, &
                                 masses, groups, molecular, remove_drift, frames, translation(:, frames), &
                                 rotational(:, frames), vibration(:, frames), linear, inertia_sum, ierr, message)
         ELSE
            CALL decompose_frame(velocity_frame%value*velocity_scale, position_values, velocity_frame, &
                                 masses, groups, molecular, remove_drift, frames, translation(:, frames), &
                                 rotational(:, frames), vibration(:, frames), linear, inertia_sum, ierr, message)
         END IF
         IF (ierr /= 0) CALL fail(message)
         DEALLOCATE (position_values)
      END DO
      CALL velocities%close_file()
      IF (position_found) THEN
         CALL positions%close_file()
      END IF
      CALL cells%close()
      IF (frames < 2) CALL fail("2PT requires at least two selected frames")

      CALL shrink_series(translation, frames)
      CALL shrink_series(rotational, frames)
      CALL shrink_series(vibration, frames)
      CALL molecular_dof(groups, linear, constraints, rotational_dof, vibrational_dof, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      molecular_mass = mass_total/REAL(molecules, dp)
      IF (explicit_volume > 0.0_dp) THEN
         selected_volume = explicit_volume
      ELSE IF (volume_frames == frames) THEN
         selected_volume = selected_volume/REAL(frames, dp)
      ELSE
         CALL fail("2PT requires --volume or cell data in the position trajectory")
      END IF
      inertia_sum = inertia_sum/REAL(frames, dp)
      CALL mean_rotational_temperatures(inertia_sum, linear, rotational_temperature)

      effective_dt = dt_fs*REAL(stride, dp)
      CALL make_dos(translation, 3*molecules, effective_dt, window, translation_dos, nfft)
      CALL make_optional_dos(rotational, rotational_dof, effective_dt, window, nfft, rotational_dos)
      CALL make_optional_dos(vibration, vibrational_dof, effective_dt, window, nfft, vibration_dos)
      ALLOCATE (frequency(0:nfft/2))
      frequency = [(REAL(ierr, dp)*1.0E15_dp/(REAL(nfft, dp)*effective_dt*2.99792458E10_dp), &
                    ierr=0, nfft/2)]

      CALL partition_twopt(translation_dos, frequency(1), temperature, molecular_mass, molecules, &
                           selected_volume, translation_gas, translation_solid, translation_result)
      IF (rotational_dof > 0) THEN
         CALL partition_twopt(rotational_dos, frequency(1), temperature, molecular_mass, molecules, &
                              selected_volume, rotational_gas, rotational_solid, rotation_result)
      ELSE
         ALLOCATE (rotational_gas(0:nfft/2), rotational_solid(0:nfft/2))
         rotational_gas = 0.0_dp
         rotational_solid = 0.0_dp
      END IF
      ALLOCATE (vibration_gas(0:nfft/2), vibration_solid(0:nfft/2))
      vibration_gas = 0.0_dp
      vibration_solid = vibration_dos
      vibration_result%dof = REAL(vibrational_dof, dp)
      vibration_result%s0_cm = vibration_dos(0)

      CALL build_twopt_weights(translation_result%packing_fraction, molecular_mass, &
         translation_result%gas_dof/3.0_dp, temperature, selected_volume, rotational_temperature, &
         rotational_symmetry, entropy_convention, weights)
      CALL integrate_twopt_channel(frequency, frequency(1), temperature, translation_gas, &
                                   translation_solid, weights%entropy_translation, translation_result)
      IF (rotational_dof > 0) &
         CALL integrate_twopt_channel(frequency, frequency(1), temperature, rotational_gas, &
                                      rotational_solid, weights%entropy_rotation, rotation_result)
      IF (vibrational_dof > 0) &
         CALL integrate_twopt_channel(frequency, frequency(1), temperature, vibration_gas, &
                                      vibration_solid, 0.0_dp, vibration_result)
      translation_result%diffusion_cm2_s = translation_result%s0_cm*gas_constant*temperature/ &
                                            (12.0_dp*2.99792458E10_dp*mass_total)*1.0E5_dp
      CALL sum_channels(translation_result, rotation_result, vibration_result, total_result)

      energy_reference = 0.0_dp
      IF (energy < 0.5_dp*HUGE(1.0_dp)) THEN
         energy_reference = energy - total_result%energy_classical_kj_mol
         total_result%energy_quantum_kj_mol = total_result%energy_quantum_kj_mol + energy_reference
         total_result%energy_classical_kj_mol = energy
         total_result%free_energy_quantum_kj_mol = total_result%free_energy_quantum_kj_mol + energy_reference
         total_result%free_energy_classical_kj_mol = total_result%free_energy_classical_kj_mol + energy_reference
      END IF
      IF (classical_cv < 0.5_dp*HUGE(1.0_dp)) THEN
         total_result%cv_quantum_j_mol_k = classical_cv + &
            total_result%cv_quantum_j_mol_k - total_result%cv_classical_j_mol_k
         total_result%cv_classical_j_mol_k = classical_cv
      END IF

      kinetic_temperature(1) = series_temperature(translation, 3*molecules)
      kinetic_temperature(2) = series_temperature(rotational, rotational_dof)
      kinetic_temperature(3) = series_temperature(vibration, vibrational_dof)
      CALL open_output(output_path, output_unit, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL write_thermodynamics(output_unit, frames, effective_dt, temperature, selected_volume, &
         molecular_mass, molecules, constraints, remove_drift, entropy_convention, kinetic_temperature, &
         rotational_temperature, energy_reference, energy < 0.5_dp*HUGE(1.0_dp), &
         translation_result, rotation_result, vibration_result, total_result)
      CALL close_output(output_unit)

      IF (spectrum_requested) THEN
         CALL open_output(spectrum_path, spectrum_unit, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         CALL write_spectrum(spectrum_unit, frequency, translation_dos, translation_gas, translation_solid, &
                             rotational_dos, rotational_gas, rotational_solid, vibration_dos)
         CALL close_output(spectrum_unit)
      END IF
      IF (vacf_requested) CALL write_vacf(vacf_path, effective_dt, translation, rotational, vibration)
   END SUBROUTINE run_twopt

   SUBROUTINE print_twopt_help()
      WRITE (*, "(A)") "Usage: twopt.x --velocity PROJECT-vel-1.xyz --mass-file masses.dat [options]"
      WRITE (*, "(A)") ""
      WRITE (*, "(A)") "Required:"
      WRITE (*, "(A)") "  --velocity FILE          CP2K XYZ velocity trajectory"
      WRITE (*, "(A)") "  --mass-file FILE         LABEL MASS[g/mol] table"
      WRITE (*, "(A)") "  --dt-fs VALUE            time between input frames"
      WRITE (*, "(A)") "  --temperature VALUE      thermodynamic temperature [K]"
      WRITE (*, "(A)") "  --volume VALUE           average volume [angstrom^3], unless cell data are available"
      WRITE (*, "(A)") ""
      WRITE (*, "(A)") "Molecular decomposition:"
      WRITE (*, "(A)") "  --position FILE          synchronized CP2K XYZ positions"
      WRITE (*, "(A)") "  --groups FILE            one molecule per line: LABEL atom..."
      WRITE (*, "(A)") "  --cell/--cell-file       optional cell for reconstructing split molecules"
      WRITE (*, "(A)") "  --rotational-symmetry N  molecular rotational symmetry number (default 1)"
      WRITE (*, "(A)") "  --constraints N          constrained internal degrees of freedom"
      WRITE (*, "(A)") ""
      WRITE (*, "(A)") "Outputs and conventions:"
      WRITE (*, "(A)") "  --output FILE            thermodynamic summary (default stdout)"
      WRITE (*, "(A)") "  --spectrum FILE          total/gas/solid density of states"
      WRITE (*, "(A)") "  --vacf FILE              normalized mass-weighted VACF diagnostics"
      WRITE (*, "(A)") "  --velocity-unit UNIT     cp2k (default), angstrom/ps, or angstrom/fs"
      WRITE (*, "(A)") "  --entropy-convention C   lin2003 (default) or rigorous"
      WRITE (*, "(A)") "  --energy-kj-mol VALUE    optional mean MD energy reference"
      WRITE (*, "(A)") "  --classical-cv-j-mol-k V optional fluctuation-derived classical Cv"
   END SUBROUTINE print_twopt_help

   LOGICAL FUNCTION groups_need_positions(groups)
      TYPE(group_type), INTENT(IN)                       :: groups(:)

      INTEGER                                            :: group

      groups_need_positions = .FALSE.
      DO group = 1, SIZE(groups)
         IF (SIZE(groups(group)%atom) > 1) groups_need_positions = .TRUE.
      END DO
   END FUNCTION groups_need_positions

   SUBROUTINE make_atomic_groups(natoms, groups)
      INTEGER, INTENT(IN)                                :: natoms
      TYPE(group_type), ALLOCATABLE, INTENT(INOUT)       :: groups(:)

      INTEGER                                            :: atom

      IF (ALLOCATED(groups)) DEALLOCATE (groups)
      ALLOCATE (groups(natoms))
      DO atom = 1, natoms
         WRITE (groups(atom)%label, "(A,I0)") "A", atom
         ALLOCATE (groups(atom)%atom(1))
         groups(atom)%atom(1) = atom
      END DO
   END SUBROUTINE make_atomic_groups

   SUBROUTINE validate_group_partition(groups, natoms, ierr, message)
      TYPE(group_type), INTENT(IN)                       :: groups(:)
      INTEGER, INTENT(IN)                                :: natoms
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      INTEGER                                            :: group, item
      INTEGER, ALLOCATABLE                               :: count(:)

      ALLOCATE (count(natoms))
      count = 0
      DO group = 1, SIZE(groups)
         DO item = 1, SIZE(groups(group)%atom)
            count(groups(group)%atom(item)) = count(groups(group)%atom(item)) + 1
         END DO
      END DO
      ierr = 0
      message = ""
      IF (ANY(count /= 1)) THEN
         ierr = 1
         message = "The group file must contain every atom exactly once"
      END IF
   END SUBROUTINE validate_group_partition

   SUBROUTINE build_masses(frame, table, masses, total, ierr, message)
      TYPE(frame_type), INTENT(IN)                       :: frame
      TYPE(mass_table_type), INTENT(IN)                  :: table
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: masses(:)
      REAL(dp), INTENT(OUT)                              :: total
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      INTEGER                                            :: atom
      LOGICAL                                            :: found

      ALLOCATE (masses(frame%natoms))
      total = 0.0_dp
      ierr = 0
      message = ""
      DO atom = 1, frame%natoms
         CALL table%lookup(frame%label(atom), masses(atom), found)
         IF (.NOT. found) THEN
            ierr = 1
            message = "No mass for trajectory label "//TRIM(frame%label(atom))
            RETURN
         END IF
         total = total + masses(atom)
      END DO
   END SUBROUTINE build_masses

   SUBROUTINE grow_series(series, new_capacity)
      REAL(dp), ALLOCATABLE, INTENT(INOUT)               :: series(:, :)
      INTEGER, INTENT(IN)                                :: new_capacity

      REAL(dp), ALLOCATABLE                              :: grown(:, :)

      ALLOCATE (grown(SIZE(series, 1), new_capacity))
      IF (SIZE(series, 2) > 0) grown(:, :SIZE(series, 2)) = series
      CALL MOVE_ALLOC(grown, series)
   END SUBROUTINE grow_series

   SUBROUTINE shrink_series(series, frames)
      REAL(dp), ALLOCATABLE, INTENT(INOUT)               :: series(:, :)
      INTEGER, INTENT(IN)                                :: frames

      REAL(dp), ALLOCATABLE                              :: shrunk(:, :)

      IF (SIZE(series, 2) == frames) RETURN
      ALLOCATE (shrunk(SIZE(series, 1), frames))
      shrunk = series(:, :frames)
      CALL MOVE_ALLOC(shrunk, series)
   END SUBROUTINE shrink_series

   SUBROUTINE decompose_frame(velocity, position, position_frame, masses, groups, molecular, remove_drift, &
                              frame_index, translation, rotation, vibration, linear, inertia_sum, ierr, message)
      REAL(dp), INTENT(IN)                               :: velocity(:, :), position(:, :), masses(:)
      TYPE(frame_type), INTENT(IN)                       :: position_frame
      TYPE(group_type), INTENT(IN)                       :: groups(:)
      LOGICAL, INTENT(IN)                                :: molecular, remove_drift
      INTEGER, INTENT(IN)                                :: frame_index
      REAL(dp), INTENT(OUT)                              :: translation(:), rotation(:), vibration(:)
      LOGICAL, INTENT(INOUT)                             :: linear(:)
      REAL(dp), INTENT(INOUT)                            :: inertia_sum(:, :)
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      INTEGER                                            :: atom, axis, component, group, item
      LOGICAL                                            :: current_linear, ok
      REAL(dp)                                           :: angular_momentum(3), center_velocity(3), corrected(3, SIZE(masses)), &
                                                            displacement(3), drift(3), eigenvalues(3), eigenvectors(3, 3), &
                                                            group_mass, inertia(3, 3), omega(3), relative(3, SIZE(masses)), &
                                                            relative_velocity(3), rotational_velocity(3)
      REAL(dp)                                           :: weighted_omega(3)

      ierr = 0
      message = ""
      translation = 0.0_dp
      IF (molecular) THEN
         rotation = 0.0_dp
         vibration = 0.0_dp
      END IF
      drift = 0.0_dp
      IF (remove_drift) THEN
         DO atom = 1, SIZE(masses)
            drift = drift + masses(atom)*velocity(:, atom)
         END DO
         drift = drift/SUM(masses)
      END IF
      DO atom = 1, SIZE(masses)
         corrected(:, atom) = velocity(:, atom) - drift
      END DO

      DO group = 1, SIZE(groups)
         group_mass = SUM(masses(groups(group)%atom))
         center_velocity = 0.0_dp
         DO item = 1, SIZE(groups(group)%atom)
            atom = groups(group)%atom(item)
            center_velocity = center_velocity + masses(atom)*corrected(:, atom)
         END DO
         center_velocity = center_velocity/group_mass
         translation(3*group - 2:3*group) = SQRT(group_mass)*center_velocity
         IF (.NOT. molecular .OR. SIZE(groups(group)%atom) == 1) CYCLE

         relative = 0.0_dp
         atom = groups(group)%atom(1)
         DO item = 2, SIZE(groups(group)%atom)
            displacement = position(:, groups(group)%atom(item)) - position(:, atom)
            IF (position_frame%cell%valid) THEN
               CALL minimum_image(position_frame%cell, displacement, relative(:, groups(group)%atom(item)), ok)
               IF (.NOT. ok) THEN
                  ierr = 1
                  message = "Singular cell while reconstructing a molecule"
                  RETURN
               END IF
            ELSE
               relative(:, groups(group)%atom(item)) = displacement
            END IF
         END DO
         displacement = 0.0_dp
         DO item = 1, SIZE(groups(group)%atom)
            atom = groups(group)%atom(item)
            displacement = displacement + masses(atom)*relative(:, atom)
         END DO
         displacement = displacement/group_mass
         DO item = 1, SIZE(groups(group)%atom)
            atom = groups(group)%atom(item)
            relative(:, atom) = relative(:, atom) - displacement
         END DO

         inertia = 0.0_dp
         angular_momentum = 0.0_dp
         DO item = 1, SIZE(groups(group)%atom)
            atom = groups(group)%atom(item)
            relative_velocity = corrected(:, atom) - center_velocity
            inertia = inertia + masses(atom)*(DOT_PRODUCT(relative(:, atom), relative(:, atom))*identity3() - &
                                              outer_product(relative(:, atom), relative(:, atom)))
            angular_momentum = angular_momentum + masses(atom)*cross(relative(:, atom), relative_velocity)
         END DO
         CALL diagonalize_symmetric3(inertia, eigenvalues, eigenvectors)
         current_linear = SIZE(groups(group)%atom) == 2 .OR. &
                          eigenvalues(1) <= 1.0E-8_dp*MAX(eigenvalues(3), TINY(1.0_dp))
         IF (frame_index == 1) THEN
            linear(group) = current_linear
         ELSE IF (linear(group) .NEQV. current_linear) THEN
            ierr = 1
            message = "A molecule changes between linear and nonlinear geometry"
            RETURN
         END IF
         inertia_sum(:, group) = inertia_sum(:, group) + eigenvalues
         omega = 0.0_dp
         weighted_omega = 0.0_dp
         DO axis = 1, 3
            IF (eigenvalues(axis) > 1.0E-10_dp*MAX(eigenvalues(3), 1.0_dp)) THEN
               omega = omega + eigenvectors(:, axis)* &
                       DOT_PRODUCT(eigenvectors(:, axis), angular_momentum)/eigenvalues(axis)
               weighted_omega = weighted_omega + eigenvectors(:, axis)* &
                                DOT_PRODUCT(eigenvectors(:, axis), angular_momentum)/SQRT(eigenvalues(axis))
            END IF
         END DO
         rotation(3*group - 2:3*group) = weighted_omega
         DO item = 1, SIZE(groups(group)%atom)
            atom = groups(group)%atom(item)
            rotational_velocity = cross(omega, relative(:, atom))
            relative_velocity = corrected(:, atom) - center_velocity - rotational_velocity
            DO component = 1, 3
               vibration(3*atom - 3 + component) = SQRT(masses(atom))*relative_velocity(component)
            END DO
         END DO
      END DO
   END SUBROUTINE decompose_frame

   SUBROUTINE molecular_dof(groups, linear, constraints, rotational_dof, vibrational_dof, ierr, message)
      TYPE(group_type), INTENT(IN)                       :: groups(:)
      LOGICAL, INTENT(IN)                                :: linear(:)
      INTEGER, INTENT(IN)                                :: constraints
      INTEGER, INTENT(OUT)                               :: rotational_dof, vibrational_dof, ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      INTEGER                                            :: group, natoms, rotation

      rotational_dof = 0
      vibrational_dof = 0
      DO group = 1, SIZE(groups)
         natoms = SIZE(groups(group)%atom)
         IF (natoms == 1) THEN
            rotation = 0
         ELSE IF (linear(group)) THEN
            rotation = 2
         ELSE
            rotation = 3
         END IF
         rotational_dof = rotational_dof + rotation
         vibrational_dof = vibrational_dof + 3*natoms - 3 - rotation
      END DO
      vibrational_dof = vibrational_dof - constraints
      ierr = 0
      message = ""
      IF (vibrational_dof < 0) THEN
         ierr = 1
         message = "--constraints exceeds the available internal degrees of freedom"
      END IF
   END SUBROUTINE molecular_dof

   SUBROUTINE mean_rotational_temperatures(inertia_sum, linear, temperatures)
      REAL(dp), INTENT(IN)                               :: inertia_sum(:, :)
      LOGICAL, INTENT(IN)                                :: linear(:)
      REAL(dp), INTENT(OUT)                              :: temperatures(3)

      INTEGER                                            :: group, rotational_molecules
      LOGICAL                                            :: all_linear
      REAL(dp)                                           :: mean_inertia(3)

      mean_inertia = 0.0_dp
      rotational_molecules = 0
      all_linear = .TRUE.
      DO group = 1, SIZE(linear)
         IF (inertia_sum(3, group) <= 0.0_dp) CYCLE
         rotational_molecules = rotational_molecules + 1
         mean_inertia = mean_inertia + inertia_sum(:, group)
         all_linear = all_linear .AND. linear(group)
      END DO
      IF (rotational_molecules == 0) THEN
         temperatures = 0.0_dp
      ELSE
         mean_inertia = mean_inertia/REAL(rotational_molecules, dp)
         CALL rotational_temperatures(mean_inertia, all_linear, temperatures)
      END IF
   END SUBROUTINE mean_rotational_temperatures

   SUBROUTINE make_dos(series, target_dof, dt_fs, window, dos, nfft)
      REAL(dp), INTENT(IN)                               :: series(:, :), dt_fs
      INTEGER, INTENT(IN)                                :: target_dof
      CHARACTER(LEN=*), INTENT(IN)                       :: window
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: dos(:)
      INTEGER, INTENT(OUT)                               :: nfft

      INTEGER                                            :: frequency
      REAL(dp)                                           :: area, frequency_step, window_norm
      REAL(dp), ALLOCATABLE                              :: power(:)

      CALL real_power_sum(series, window, power, nfft, window_norm)
      DO frequency = 1, nfft/2 - 1
         power(frequency) = 2.0_dp*power(frequency)
      END DO
      frequency_step = 1.0E15_dp/(REAL(nfft, dp)*dt_fs*2.99792458E10_dp)
      area = frequency_step*(0.5_dp*(power(0) + power(nfft/2)) + SUM(power(1:nfft/2 - 1)))
      IF (area <= TINY(1.0_dp)) CALL fail("A 2PT velocity channel has zero spectral power")
      ALLOCATE (dos(0:nfft/2))
      dos = power*REAL(target_dof, dp)/area
   END SUBROUTINE make_dos

   SUBROUTINE make_optional_dos(series, target_dof, dt_fs, window, expected_nfft, dos)
      REAL(dp), INTENT(IN)                               :: series(:, :), dt_fs
      INTEGER, INTENT(IN)                                :: target_dof, expected_nfft
      CHARACTER(LEN=*), INTENT(IN)                       :: window
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: dos(:)

      INTEGER                                            :: nfft

      IF (target_dof == 0) THEN
         ALLOCATE (dos(0:expected_nfft/2))
         dos = 0.0_dp
      ELSE
         CALL make_dos(series, target_dof, dt_fs, window, dos, nfft)
         IF (nfft /= expected_nfft) CALL fail("Internal FFT size mismatch")
      END IF
   END SUBROUTINE make_optional_dos

   REAL(dp) FUNCTION series_temperature(series, dof)
      REAL(dp), INTENT(IN)                               :: series(:, :)
      INTEGER, INTENT(IN)                                :: dof

      IF (dof <= 0 .OR. SIZE(series) == 0) THEN
         series_temperature = 0.0_dp
      ELSE
         series_temperature = SUM(series**2)/REAL(SIZE(series, 2), dp)/(0.1_dp*gas_constant*REAL(dof, dp))
      END IF
   END FUNCTION series_temperature

   SUBROUTINE sum_channels(translation, rotation, vibration, total)
      TYPE(twopt_channel_type), INTENT(IN)                :: translation, rotation, vibration
      TYPE(twopt_channel_type), INTENT(OUT)               :: total

      total%dof = translation%dof + rotation%dof + vibration%dof
      total%gas_dof = translation%gas_dof + rotation%gas_dof
      total%zpe_kj_mol = translation%zpe_kj_mol + rotation%zpe_kj_mol + vibration%zpe_kj_mol
      total%energy_quantum_kj_mol = translation%energy_quantum_kj_mol + &
                                    rotation%energy_quantum_kj_mol + vibration%energy_quantum_kj_mol
      total%energy_classical_kj_mol = translation%energy_classical_kj_mol + &
                                      rotation%energy_classical_kj_mol + vibration%energy_classical_kj_mol
      total%entropy_quantum_j_mol_k = translation%entropy_quantum_j_mol_k + &
                                      rotation%entropy_quantum_j_mol_k + vibration%entropy_quantum_j_mol_k
      total%entropy_classical_j_mol_k = translation%entropy_classical_j_mol_k + &
                                        rotation%entropy_classical_j_mol_k + vibration%entropy_classical_j_mol_k
      total%free_energy_quantum_kj_mol = translation%free_energy_quantum_kj_mol + &
                                         rotation%free_energy_quantum_kj_mol + vibration%free_energy_quantum_kj_mol
      total%free_energy_classical_kj_mol = translation%free_energy_classical_kj_mol + &
                                           rotation%free_energy_classical_kj_mol + vibration%free_energy_classical_kj_mol
      total%cv_quantum_j_mol_k = translation%cv_quantum_j_mol_k + &
                                 rotation%cv_quantum_j_mol_k + vibration%cv_quantum_j_mol_k
      total%cv_classical_j_mol_k = translation%cv_classical_j_mol_k + &
                                   rotation%cv_classical_j_mol_k + vibration%cv_classical_j_mol_k
   END SUBROUTINE sum_channels

   SUBROUTINE write_thermodynamics(unit, frames, dt_fs, temperature, volume, molecular_mass, molecules, &
                                   constraints, remove_drift, convention, kinetic_temperature, &
                                   rotational_temperature, energy_reference, has_energy, translation, &
                                   rotation, vibration, total)
      INTEGER, INTENT(IN)                                :: unit, frames, molecules, constraints
      REAL(dp), INTENT(IN)                               :: dt_fs, temperature, volume, molecular_mass, &
                                                            kinetic_temperature(3), rotational_temperature(3), &
                                                            energy_reference
      LOGICAL, INTENT(IN)                                :: remove_drift, has_energy
      CHARACTER(LEN=*), INTENT(IN)                       :: convention
      TYPE(twopt_channel_type), INTENT(IN)                :: translation, rotation, vibration, total

      WRITE (unit, "(A)") "# Conventional two-phase thermodynamics (2PT)"
      WRITE (unit, "(A,I0)") "# frames: ", frames
      WRITE (unit, "(A,ES24.16)") "# frame_dt_fs: ", dt_fs
      WRITE (unit, "(A,ES24.16)") "# temperature_K: ", temperature
      WRITE (unit, "(A,ES24.16)") "# volume_angstrom3: ", volume
      WRITE (unit, "(A,I0)") "# molecules: ", molecules
      WRITE (unit, "(A,ES24.16)") "# mean_molecular_mass_g_mol: ", molecular_mass
      WRITE (unit, "(A,I0)") "# constraints: ", constraints
      WRITE (unit, "(A,L1)") "# removed_system_com_drift: ", remove_drift
      WRITE (unit, "(A,A)") "# entropy_convention: ", TRIM(convention)
      WRITE (unit, "(A,3(1X,ES24.16))") "# kinetic_temperature_K (trans rot vib):", kinetic_temperature
      WRITE (unit, "(A,3(1X,ES24.16))") "# rotational_temperature_K:", rotational_temperature
      IF (has_energy) THEN
         WRITE (unit, "(A,ES24.16)") "# potential_reference_shift_kJ_mol: ", energy_reference
      ELSE
         WRITE (unit, "(A)") "# no MD energy reference: E and A are spectral contributions only"
      END IF
      WRITE (unit, "(A)") "# channel dof gas_dof s0[cm] K fluidicity packing D[cm2/s] "// &
         "ZPE[kJ/mol] Eq[kJ/mol] Ec[kJ/mol] Sq[J/mol/K] Sc[J/mol/K] Aq[kJ/mol] Ac[kJ/mol] "// &
         "Cvq[J/mol/K] Cvc[J/mol/K]"
      CALL write_channel(unit, "translation", translation)
      CALL write_channel(unit, "rotation", rotation)
      CALL write_channel(unit, "vibration", vibration)
      CALL write_channel(unit, "total", total)
   END SUBROUTINE write_thermodynamics

   SUBROUTINE write_channel(unit, label, result)
      INTEGER, INTENT(IN)                                :: unit
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      TYPE(twopt_channel_type), INTENT(IN)                :: result

      WRITE (unit, "(A,17(1X,ES24.16))") TRIM(label), result%dof, result%gas_dof, result%s0_cm, &
         result%k_parameter, result%fluidicity, result%packing_fraction, result%diffusion_cm2_s, &
         result%zpe_kj_mol, result%energy_quantum_kj_mol, result%energy_classical_kj_mol, &
         result%entropy_quantum_j_mol_k, result%entropy_classical_j_mol_k, &
         result%free_energy_quantum_kj_mol, result%free_energy_classical_kj_mol, &
         result%cv_quantum_j_mol_k, result%cv_classical_j_mol_k
   END SUBROUTINE write_channel

   SUBROUTINE write_spectrum(unit, frequency, translation, translation_gas, translation_solid, &
                             rotation, rotation_gas, rotation_solid, vibration)
      INTEGER, INTENT(IN)                                :: unit
      REAL(dp), INTENT(IN)                               :: frequency(0:), translation(0:), &
                                                            translation_gas(0:), translation_solid(0:), &
                                                            rotation(0:), rotation_gas(0:), &
                                                            rotation_solid(0:), vibration(0:)

      INTEGER                                            :: index

      WRITE (unit, "(A)") "# wavenumber[cm^-1] total translation translation_gas translation_solid "// &
         "rotation rotation_gas rotation_solid vibration; all DoS columns in cm"
      DO index = 0, UBOUND(frequency, 1)
         WRITE (unit, "(9ES24.16)") frequency(index), translation(index) + rotation(index) + vibration(index), &
            translation(index), translation_gas(index), translation_solid(index), rotation(index), &
            rotation_gas(index), rotation_solid(index), vibration(index)
      END DO
   END SUBROUTINE write_spectrum

   SUBROUTINE write_vacf(path, dt_fs, translation, rotation, vibration)
      CHARACTER(LEN=*), INTENT(IN)                       :: path
      REAL(dp), INTENT(IN)                               :: dt_fs, translation(:, :), rotation(:, :), vibration(:, :)

      CHARACTER(LEN=512)                                 :: message
      INTEGER                                            :: ierr, lag, unit
      REAL(dp), ALLOCATABLE                              :: correlation_rotation(:), correlation_translation(:), &
                                                            correlation_vibration(:)

      CALL normalized_correlation(translation, correlation_translation)
      CALL normalized_correlation(rotation, correlation_rotation)
      CALL normalized_correlation(vibration, correlation_vibration)
      CALL open_output(path, unit, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      WRITE (unit, "(A)") "# lag[fs] normalized_translation normalized_rotation normalized_vibration"
      DO lag = 0, SIZE(translation, 2) - 1
         WRITE (unit, "(4ES24.16)") REAL(lag, dp)*dt_fs, correlation_translation(lag), &
            correlation_rotation(lag), correlation_vibration(lag)
      END DO
      CALL close_output(unit)
   END SUBROUTINE write_vacf

   SUBROUTINE normalized_correlation(series, correlation)
      REAL(dp), INTENT(IN)                               :: series(:, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: correlation(:)

      INTEGER                                            :: lag

      IF (SIZE(series, 1) == 0) THEN
         ALLOCATE (correlation(0:SIZE(series, 2) - 1))
         correlation = 0.0_dp
         RETURN
      END IF
      CALL real_autocorrelation_sum(series, SIZE(series, 2) - 1, correlation)
      DO lag = 0, UBOUND(correlation, 1)
         correlation(lag) = correlation(lag)/REAL(SIZE(series, 2) - lag, dp)
      END DO
      IF (correlation(0) > TINY(1.0_dp)) correlation = correlation/correlation(0)
   END SUBROUTINE normalized_correlation

   PURE FUNCTION identity3() RESULT(identity)
      REAL(dp)                                           :: identity(3, 3)

      identity = 0.0_dp
      identity(1, 1) = 1.0_dp
      identity(2, 2) = 1.0_dp
      identity(3, 3) = 1.0_dp
   END FUNCTION identity3

   PURE FUNCTION outer_product(a, b) RESULT(product)
      REAL(dp), INTENT(IN)                               :: a(3), b(3)
      REAL(dp)                                           :: product(3, 3)

      INTEGER                                            :: column

      DO column = 1, 3
         product(:, column) = a*b(column)
      END DO
   END FUNCTION outer_product

   PURE FUNCTION cross(a, b) RESULT(product)
      REAL(dp), INTENT(IN)                               :: a(3), b(3)
      REAL(dp)                                           :: product(3)

      product(1) = a(2)*b(3) - a(3)*b(2)
      product(2) = a(3)*b(1) - a(1)*b(3)
      product(3) = a(1)*b(2) - a(2)*b(1)
   END FUNCTION cross

   SUBROUTINE diagonalize_symmetric3(matrix, eigenvalues, eigenvectors)
      REAL(dp), INTENT(IN)                               :: matrix(3, 3)
      REAL(dp), INTENT(OUT)                              :: eigenvalues(3), eigenvectors(3, 3)

      INTEGER                                            :: iteration, p, q
      REAL(dp)                                           :: c, maximum, rotated(3, 3), s, tau, t, work(3, 3)

      work = matrix
      eigenvectors = identity3()
      DO iteration = 1, 50
         CALL largest_offdiagonal(work, p, q, maximum)
         IF (maximum <= 1.0E-14_dp*MAX(MAXVAL(ABS(work)), 1.0_dp)) EXIT
         tau = (work(q, q) - work(p, p))/(2.0_dp*work(p, q))
         IF (tau >= 0.0_dp) THEN
            t = 1.0_dp/(tau + SQRT(1.0_dp + tau*tau))
         ELSE
            t = -1.0_dp/(-tau + SQRT(1.0_dp + tau*tau))
         END IF
         c = 1.0_dp/SQRT(1.0_dp + t*t)
         s = t*c
         rotated = identity3()
         rotated(p, p) = c
         rotated(q, q) = c
         rotated(p, q) = s
         rotated(q, p) = -s
         work = MATMUL(TRANSPOSE(rotated), MATMUL(work, rotated))
         eigenvectors = MATMUL(eigenvectors, rotated)
      END DO
      eigenvalues = [work(1, 1), work(2, 2), work(3, 3)]
      CALL sort_eigensystem(eigenvalues, eigenvectors)
      eigenvalues = MAX(eigenvalues, 0.0_dp)
   END SUBROUTINE diagonalize_symmetric3

   SUBROUTINE largest_offdiagonal(matrix, p, q, maximum)
      REAL(dp), INTENT(IN)                               :: matrix(3, 3)
      INTEGER, INTENT(OUT)                               :: p, q
      REAL(dp), INTENT(OUT)                              :: maximum

      INTEGER                                            :: i, j

      p = 1
      q = 2
      maximum = ABS(matrix(p, q))
      DO i = 1, 2
         DO j = i + 1, 3
            IF (ABS(matrix(i, j)) > maximum) THEN
               maximum = ABS(matrix(i, j))
               p = i
               q = j
            END IF
         END DO
      END DO
   END SUBROUTINE largest_offdiagonal

   SUBROUTINE sort_eigensystem(values, vectors)
      REAL(dp), INTENT(INOUT)                            :: values(3), vectors(3, 3)

      INTEGER                                            :: i, j
      REAL(dp)                                           :: temporary, vector(3)

      DO i = 1, 2
         DO j = i + 1, 3
            IF (values(j) < values(i)) THEN
               temporary = values(i)
               values(i) = values(j)
               values(j) = temporary
               vector = vectors(:, i)
               vectors(:, i) = vectors(:, j)
               vectors(:, j) = vector
            END IF
         END DO
      END DO
   END SUBROUTINE sort_eigensystem

END MODULE trajana_twopt_analysis
