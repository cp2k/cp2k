!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_rdf_analysis
   USE trajana_cell_source,             ONLY: cell_source_type
   USE trajana_command_line,            ONLY: fail,&
                                              get_integer_option,&
                                              get_option,&
                                              get_real_option
   USE trajana_frame_controls,          ONLY: frame_selected
   USE trajana_geometry,                ONLY: cell_volume,&
                                              minimum_image,&
                                              shortest_cell_height
   USE trajana_kinds,                   ONLY: dp
   USE trajana_selections,              ONLY: build_selection
   USE trajana_trajectory_io,           ONLY: close_output,&
                                              open_output,&
                                              xyz_reader_type
   USE trajana_trajectory_types,        ONLY: frame_type

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: run_rdf

CONTAINS

   SUBROUTINE run_rdf()
      CHARACTER(LEN=512)                                 :: message
      CHARACTER(LEN=:), ALLOCATABLE                      :: cell_path, cell_text, input_path, &
                                                            output_path, periodic_text, &
                                                            selection_a, selection_b
      INTEGER                                            :: atom_a, atom_b, bin, bins, &
                                                            expected_atoms, first, frames, ierr, &
                                                            last, output_unit, overlap, stride
      LOGICAL                                            :: eof, found, initialized, ok
      LOGICAL, ALLOCATABLE                               :: selected_a(:), selected_b(:)
      REAL(dp)                                           :: bin_width, cumulative, displacement(3), &
                                                            distance, NEAREST(3), pair_count, pi, &
                                                            reference_count, rmax, shell, volume
      REAL(dp), ALLOCATABLE                              :: histogram(:), normalization(:)
      TYPE(cell_source_type)                             :: cells
      TYPE(frame_type)                                   :: frame
      TYPE(xyz_reader_type)                              :: trajectory

      CALL get_option("--input", input_path, found, "-")
      CALL get_option("--output", output_path, found, "-")
      CALL get_option("--cell", cell_text, found, "")
      CALL get_option("--cell-file", cell_path, found, "")
      CALL get_option("--periodic", periodic_text, found, "XYZ")
      CALL get_option("--select-a", selection_a, found, "all")
      CALL get_option("--select-b", selection_b, found, "all")
      CALL get_integer_option("--bins", bins, 200)
      CALL get_real_option("--rmax", rmax, -1.0_dp)
      CALL get_integer_option("--first", first, 1)
      CALL get_integer_option("--last", last, HUGE(1))
      CALL get_integer_option("--stride", stride, 1)
      IF (bins < 1 .OR. first < 1 .OR. last < first .OR. stride < 1) CALL fail("Invalid RDF options")

      CALL cells%configure(cell_text, cell_path, periodic_text, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL trajectory%open_file(input_path, ierr, message)
      IF (ierr /= 0) CALL fail(message)

      initialized = .FALSE.
      frames = 0
      reference_count = 0.0_dp
      expected_atoms = -1
      pi = ACOS(-1.0_dp)

      DO
         CALL trajectory%read_frame(frame, eof, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         IF (eof) EXIT
         CALL cells%attach(frame, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         IF (frame%number > last) EXIT
         IF (.NOT. frame_selected(frame%number, first, last, stride)) CYCLE
         IF (.NOT. frame%cell%valid) CALL fail("rdf requires --cell, --cell-file, or an extended-XYZ Lattice")
         IF (.NOT. ALL(frame%cell%periodic)) CALL fail("rdf normalization requires periodicity in XYZ")
         IF (expected_atoms < 0) expected_atoms = frame%natoms
         IF (frame%natoms /= expected_atoms) CALL fail("The atom count changes along the trajectory")

         CALL build_selection(selection_a, frame%label, selected_a, ierr, message)
         IF (ierr /= 0) CALL fail("Invalid --select-a: "//TRIM(message))
         CALL build_selection(selection_b, frame%label, selected_b, ierr, message)
         IF (ierr /= 0) CALL fail("Invalid --select-b: "//TRIM(message))

         IF (.NOT. initialized) THEN
            IF (rmax <= 0.0_dp) rmax = 0.5_dp*shortest_cell_height(frame%cell)
            IF (rmax <= 0.0_dp) CALL fail("Cannot determine a positive RDF cutoff")
            bin_width = rmax/REAL(bins, dp)
            ALLOCATE (histogram(bins), normalization(bins))
            histogram = 0.0_dp
            normalization = 0.0_dp
            initialized = .TRUE.
         END IF
         IF (2.0_dp*rmax > shortest_cell_height(frame%cell)*(1.0_dp + 100.0_dp*EPSILON(1.0_dp))) &
            CALL fail("--rmax exceeds half the shortest periodic cell height")

         volume = cell_volume(frame%cell)
         IF (volume <= TINY(1.0_dp)) CALL fail("Cell volume is zero")
         overlap = COUNT(selected_a .AND. selected_b)
         pair_count = REAL(COUNT(selected_a), dp)*REAL(COUNT(selected_b), dp) - REAL(overlap, dp)
         IF (pair_count <= 0.0_dp) CALL fail("Selections contain no distinct atom pairs")

         DO atom_a = 1, frame%natoms
            IF (.NOT. selected_a(atom_a)) CYCLE
            DO atom_b = 1, frame%natoms
               IF (.NOT. selected_b(atom_b) .OR. atom_a == atom_b) CYCLE
               displacement = frame%value(:, atom_b) - frame%value(:, atom_a)
               CALL minimum_image(frame%cell, displacement, nearest, ok)
               IF (.NOT. ok) CALL fail("Singular cell in trajectory")
               distance = NORM2(nearest)
               IF (distance >= rmax) CYCLE
               bin = INT(distance/bin_width) + 1
               histogram(bin) = histogram(bin) + 1.0_dp
            END DO
         END DO

         DO bin = 1, bins
            shell = 4.0_dp*pi*((REAL(bin, dp)*bin_width)**3 - (REAL(bin - 1, dp)*bin_width)**3)/3.0_dp
            normalization(bin) = normalization(bin) + pair_count*shell/volume
         END DO
         reference_count = reference_count + REAL(COUNT(selected_a), dp)
         frames = frames + 1
      END DO

      IF (.NOT. initialized .OR. frames == 0) CALL fail("No frames selected for RDF analysis")
      CALL open_output(output_path, output_unit, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      WRITE (output_unit, "(A)") "# r [angstrom]   g_ab(r)   coordination_number_of_b_around_a"
      WRITE (output_unit, "(A,A)") "# selection_a: ", TRIM(selection_a)
      WRITE (output_unit, "(A,A)") "# selection_b: ", TRIM(selection_b)
      cumulative = 0.0_dp
      DO bin = 1, bins
         cumulative = cumulative + histogram(bin)
         IF (normalization(bin) > 0.0_dp) THEN
            WRITE (output_unit, "(3ES24.16)") (REAL(bin, dp) - 0.5_dp)*bin_width, &
               histogram(bin)/normalization(bin), cumulative/reference_count
         ELSE
            WRITE (output_unit, "(3ES24.16)") (REAL(bin, dp) - 0.5_dp)*bin_width, 0.0_dp, &
               cumulative/reference_count
         END IF
      END DO

      CALL close_output(output_unit)
      CALL trajectory%close_file()
      CALL cells%close()
   END SUBROUTINE run_rdf

END MODULE trajana_rdf_analysis
