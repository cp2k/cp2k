!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_hbond_analysis
   USE trajana_cell_source,             ONLY: cell_source_type
   USE trajana_command_line,            ONLY: fail,&
                                              get_integer_option,&
                                              get_option,&
                                              get_real_option
   USE trajana_frame_controls,          ONLY: frame_selected
   USE trajana_geometry,                ONLY: minimum_image
   USE trajana_groups,                  ONLY: group_type,&
                                              read_groups,&
                                              validate_groups
   USE trajana_kinds,                   ONLY: dp
   USE trajana_trajectory_io,           ONLY: close_output,&
                                              open_output,&
                                              xyz_reader_type
   USE trajana_trajectory_types,        ONLY: frame_type

   IMPLICIT NONE
   PRIVATE

   TYPE :: population_type
      INTEGER, ALLOCATABLE :: acceptors(:)
      INTEGER, ALLOCATABLE :: donors(:)
      INTEGER, ALLOCATABLE :: count(:)
   END TYPE population_type

   PUBLIC :: run_hbond

CONTAINS

   SUBROUTINE run_hbond()
      CHARACTER(LEN=512)                                 :: message
      CHARACTER(LEN=:), ALLOCATABLE                      :: cell_path, cell_text, group_path, &
                                                            input_path, output_path, &
                                                            periodic_text, summary_path
      INTEGER :: acceptor, bonds, donor, expected_atoms, first, frames, free_hydrogens, group, &
         hydrogen, hydrogen_index, ierr, last, output_unit, stride, summary_unit, total_hydrogens
      INTEGER, ALLOCATABLE                               :: accepted(:), donated(:), &
                                                            hydrogen_bonds(:), offsets(:)
      LOGICAL                                            :: eof, found, initialized, ok, &
                                                            summary_requested
      REAL(dp) :: angle, angle_coefficient, angle_max, cosine, delta, limit, m2_bonds, mean_bonds, &
         oh(3), oh_length, oh_max, oo(3), oo_length, roo_zero, std_bonds
      TYPE(cell_source_type)                             :: cells
      TYPE(frame_type)                                   :: frame
      TYPE(group_type), ALLOCATABLE                      :: groups(:)
      TYPE(population_type)                              :: population
      TYPE(xyz_reader_type)                              :: trajectory

      CALL get_option("--input", input_path, found, "-")
      CALL get_option("--output", output_path, found, "-")
      CALL get_option("--summary", summary_path, summary_requested)
      CALL get_option("--groups", group_path, found)
      IF (.NOT. found) CALL fail("hbond requires --groups FILE")
      CALL get_option("--cell", cell_text, found, "")
      CALL get_option("--cell-file", cell_path, found, "")
      CALL get_option("--periodic", periodic_text, found, "XYZ")
      CALL get_real_option("--roo-zero", roo_zero, 3.3_dp)
      CALL get_real_option("--angle-coefficient", angle_coefficient, 0.00044_dp)
      CALL get_real_option("--angle-max", angle_max, 45.0_dp)
      CALL get_real_option("--oh-max", oh_max, 1.25_dp)
      CALL get_integer_option("--first", first, 1)
      CALL get_integer_option("--last", last, HUGE(1))
      CALL get_integer_option("--stride", stride, 1)
      IF (roo_zero <= 0.0_dp .OR. angle_coefficient < 0.0_dp .OR. angle_max <= 0.0_dp .OR. oh_max <= 0.0_dp) &
         CALL fail("Invalid hydrogen-bond criterion")
      IF (first < 1 .OR. last < first .OR. stride < 1) CALL fail("Invalid frame range")

      CALL read_groups(group_path, groups, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL cells%configure(cell_text, cell_path, periodic_text, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL trajectory%open_file(input_path, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL open_output(output_path, output_unit, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      WRITE (output_unit, "(A)") "# frame   hydrogen_bonds   bonds_per_group   free_donor_H_fraction"

      frames = 0
      expected_atoms = -1
      mean_bonds = 0.0_dp
      m2_bonds = 0.0_dp
      initialized = .FALSE.

      DO
         CALL trajectory%read_frame(frame, eof, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         IF (eof) EXIT
         CALL cells%attach(frame, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         IF (frame%number > last) EXIT
         IF (.NOT. frame_selected(frame%number, first, last, stride)) CYCLE
         IF (.NOT. frame%cell%valid) CALL fail("hbond requires --cell, --cell-file, or an extended-XYZ Lattice")
         IF (expected_atoms < 0) expected_atoms = frame%natoms
         IF (frame%natoms /= expected_atoms) CALL fail("The atom count changes along the trajectory")

         IF (.NOT. initialized) THEN
            CALL validate_groups(groups, frame%natoms, 2, ierr, message)
            IF (ierr /= 0) CALL fail(message)
            total_hydrogens = SUM([(SIZE(groups(group)%atom) - 1, group=1, SIZE(groups))])
            ALLOCATE (accepted(SIZE(groups)), donated(SIZE(groups)), hydrogen_bonds(total_hydrogens), offsets(SIZE(groups)))
            offsets(1) = 0
            DO group = 2, SIZE(groups)
               offsets(group) = offsets(group - 1) + SIZE(groups(group - 1)%atom) - 1
            END DO
            ALLOCATE (population%acceptors(0), population%donors(0), population%count(0))
            initialized = .TRUE.
         END IF

         accepted = 0
         donated = 0
         hydrogen_bonds = 0
         bonds = 0
         DO donor = 1, SIZE(groups)
            DO hydrogen = 2, SIZE(groups(donor)%atom)
               CALL minimum_image(frame%cell, &
                                  frame%value(:, groups(donor)%atom(hydrogen)) - frame%value(:, groups(donor)%atom(1)), oh, ok)
               IF (.NOT. ok) CALL fail("Singular cell in trajectory")
               oh_length = NORM2(oh)
               IF (oh_length <= TINY(1.0_dp) .OR. oh_length > oh_max) CYCLE
               hydrogen_index = offsets(donor) + hydrogen - 1
               DO acceptor = 1, SIZE(groups)
                  IF (acceptor == donor) CYCLE
                  CALL minimum_image(frame%cell, &
                                     frame%value(:, groups(acceptor)%atom(1)) - frame%value(:, groups(donor)%atom(1)), oo, ok)
                  IF (.NOT. ok) CALL fail("Singular cell in trajectory")
                  oo_length = NORM2(oo)
                  IF (oo_length <= TINY(1.0_dp)) CYCLE
                  cosine = MAX(-1.0_dp, MIN(1.0_dp, DOT_PRODUCT(oh, oo)/(oh_length*oo_length)))
                  angle = ACOS(cosine)*180.0_dp/ACOS(-1.0_dp)
                  limit = roo_zero - angle_coefficient*angle**2
                  IF (angle <= angle_max .AND. oo_length < limit) THEN
                     bonds = bonds + 1
                     donated(donor) = donated(donor) + 1
                     accepted(acceptor) = accepted(acceptor) + 1
                     hydrogen_bonds(hydrogen_index) = hydrogen_bonds(hydrogen_index) + 1
                  END IF
               END DO
            END DO
         END DO

         free_hydrogens = COUNT(hydrogen_bonds == 0)
         DO group = 1, SIZE(groups)
            CALL add_population(population, accepted(group), donated(group))
         END DO
         frames = frames + 1
         delta = REAL(bonds, dp) - mean_bonds
         mean_bonds = mean_bonds + delta/REAL(frames, dp)
         m2_bonds = m2_bonds + delta*(REAL(bonds, dp) - mean_bonds)
         WRITE (output_unit, "(2I12,2ES24.16)") frame%number, bonds, REAL(bonds, dp)/REAL(SIZE(groups), dp), &
            REAL(free_hydrogens, dp)/REAL(total_hydrogens, dp)
      END DO

      IF (.NOT. initialized .OR. frames == 0) CALL fail("No frames selected for hydrogen-bond analysis")
      IF (frames > 1) THEN
         std_bonds = SQRT(m2_bonds/REAL(frames - 1, dp))
      ELSE
         std_bonds = 0.0_dp
      END IF
      WRITE (output_unit, "(A,ES24.16)") "# mean hydrogen bonds per frame: ", mean_bonds
      WRITE (output_unit, "(A,ES24.16)") "# sample standard deviation: ", std_bonds
      WRITE (output_unit, "(A,ES24.16)") "# mean hydrogen bonds per group: ", mean_bonds/REAL(SIZE(groups), dp)
      CALL close_output(output_unit)

      IF (summary_requested) THEN
         CALL open_output(summary_path, summary_unit, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         WRITE (summary_unit, "(A)") "# accepted   donated   fraction_of_group_frames"
         DO group = 1, SIZE(population%count)
            WRITE (summary_unit, "(2I12,ES24.16)") population%acceptors(group), population%donors(group), &
               REAL(population%count(group), dp)/(REAL(frames, dp)*REAL(SIZE(groups), dp))
         END DO
         CALL close_output(summary_unit)
      END IF

      CALL trajectory%close_file()
      CALL cells%close()
   END SUBROUTINE run_hbond

   SUBROUTINE add_population(population, acceptors, donors)
      TYPE(population_type), INTENT(INOUT)               :: population
      INTEGER, INTENT(IN)                                :: acceptors, donors

      INTEGER                                            :: count, item
      INTEGER, ALLOCATABLE                               :: grown(:)

      DO item = 1, SIZE(population%count)
         IF (population%acceptors(item) == acceptors .AND. population%donors(item) == donors) THEN
            population%count(item) = population%count(item) + 1
            RETURN
         END IF
      END DO
      count = SIZE(population%count)
      ALLOCATE (grown(count + 1))
      IF (count > 0) grown(:count) = population%acceptors
      grown(count + 1) = acceptors
      CALL MOVE_ALLOC(grown, population%acceptors)
      ALLOCATE (grown(count + 1))
      IF (count > 0) grown(:count) = population%donors
      grown(count + 1) = donors
      CALL MOVE_ALLOC(grown, population%donors)
      ALLOCATE (grown(count + 1))
      IF (count > 0) grown(:count) = population%count
      grown(count + 1) = 1
      CALL MOVE_ALLOC(grown, population%count)
   END SUBROUTINE add_population

END MODULE trajana_hbond_analysis
