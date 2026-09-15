!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_molecular_motion
   USE trajana_geometry,                ONLY: minimum_image
   USE trajana_groups,                  ONLY: group_type,&
                                              mass_table_type
   USE trajana_kinds,                   ONLY: dp
   USE trajana_linear_algebra,          ONLY: diagonalize_symmetric
   USE trajana_trajectory_types,        ONLY: frame_type

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: build_atom_masses, check_homogeneous_groups, molecule_motion

CONTAINS

   SUBROUTINE build_atom_masses(frame, table, masses, ierr, message)
      TYPE(frame_type), INTENT(IN)                        :: frame
      TYPE(mass_table_type), INTENT(IN)                   :: table
      REAL(dp), ALLOCATABLE, INTENT(OUT)                  :: masses(:)
      INTEGER, INTENT(OUT)                                :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                       :: message

      INTEGER                                             :: atom
      LOGICAL                                             :: found

      ierr = 0
      message = ""
      ALLOCATE (masses(frame%natoms))
      DO atom = 1, frame%natoms
         CALL table%lookup(frame%label(atom), masses(atom), found)
         IF (.NOT. found) THEN
            ierr = 1
            message = "No mass found for trajectory label "//TRIM(frame%label(atom))
            RETURN
         END IF
      END DO
   END SUBROUTINE build_atom_masses

   SUBROUTINE check_homogeneous_groups(groups, labels, masses, molecular_mass, ierr, message)
      TYPE(group_type), INTENT(IN)                        :: groups(:)
      CHARACTER(LEN=*), INTENT(IN)                       :: labels(:)
      REAL(dp), INTENT(IN)                               :: masses(:)
      REAL(dp), INTENT(OUT)                              :: molecular_mass
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      INTEGER                                             :: atom, group, item
      LOGICAL, ALLOCATABLE                               :: used(:)

      ierr = 0
      message = ""
      molecular_mass = 0.0_dp
      ALLOCATE (used(SIZE(masses)))
      used = .FALSE.
      IF (SIZE(groups) == 0) THEN
         ierr = 1
         message = "At least one solvent molecule is required"
         RETURN
      END IF
      IF (SIZE(groups(1)%atom) < 2) THEN
         ierr = 1
         message = "3D-2PT requires molecular groups containing at least two atoms"
         RETURN
      END IF
      DO group = 1, SIZE(groups)
         IF (SIZE(groups(group)%atom) /= SIZE(groups(1)%atom)) THEN
            ierr = 1
            message = "All 3D-2PT solvent groups must have the same number of atoms"
            RETURN
         END IF
         DO item = 1, SIZE(groups(group)%atom)
            atom = groups(group)%atom(item)
            IF (atom < 1 .OR. atom > SIZE(masses)) THEN
               ierr = 1
               message = "A solvent group contains an invalid atom index"
               RETURN
            END IF
            IF (used(atom)) THEN
               ierr = 1
               message = "Solvent groups must not share atoms"
               RETURN
            END IF
            used(atom) = .TRUE.
            IF (labels(atom) /= labels(groups(1)%atom(item)) .OR. &
                ABS(masses(atom) - masses(groups(1)%atom(item))) > &
                1.0E-10_dp*MAX(masses(atom), masses(groups(1)%atom(item)))) THEN
               ierr = 1
               message = "All 3D-2PT solvent groups must use the same ordered atom labels and masses"
               RETURN
            END IF
         END DO
      END DO
      molecular_mass = SUM(masses(groups(1)%atom))
   END SUBROUTINE check_homogeneous_groups

   SUBROUTINE molecule_motion(position_frame, velocity, masses, group, center, translation, rotation, &
                              inertia_values, principal_axes, linear, ierr, message)
      TYPE(frame_type), INTENT(IN)                        :: position_frame
      REAL(dp), INTENT(IN)                               :: velocity(:, :), masses(:)
      TYPE(group_type), INTENT(IN)                        :: group
      REAL(dp), INTENT(OUT)                              :: center(3), translation(3), rotation(3), &
                                                            inertia_values(3), principal_axes(3, 3)
      LOGICAL, INTENT(OUT)                               :: linear
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      INTEGER                                             :: atom, axis, first_atom, item
      LOGICAL                                             :: ok
      REAL(dp)                                           :: angular_momentum(3), center_offset(3), &
                                                            center_velocity(3), displacement(3), &
                                                            group_mass, inertia(3, 3), relative_velocity(3)
      REAL(dp), ALLOCATABLE                              :: eigenvalues(:), eigenvectors(:, :), relative(:, :)

      ierr = 0
      message = ""
      group_mass = SUM(masses(group%atom))
      ALLOCATE (relative(3, SIZE(group%atom)))
      relative = 0.0_dp
      first_atom = group%atom(1)
      DO item = 2, SIZE(group%atom)
         displacement = position_frame%value(:, group%atom(item)) - position_frame%value(:, first_atom)
         IF (position_frame%cell%valid) THEN
            CALL minimum_image(position_frame%cell, displacement, relative(:, item), ok)
            IF (.NOT. ok) THEN
               ierr = 1
               message = "Singular cell while reconstructing a solvent molecule"
               RETURN
            END IF
         ELSE
            relative(:, item) = displacement
         END IF
      END DO

      center_offset = 0.0_dp
      center_velocity = 0.0_dp
      DO item = 1, SIZE(group%atom)
         atom = group%atom(item)
         center_offset = center_offset + masses(atom)*relative(:, item)
         center_velocity = center_velocity + masses(atom)*velocity(:, atom)
      END DO
      center_offset = center_offset/group_mass
      center_velocity = center_velocity/group_mass
      center = position_frame%value(:, first_atom) + center_offset
      DO item = 1, SIZE(group%atom)
         relative(:, item) = relative(:, item) - center_offset
      END DO

      inertia = 0.0_dp
      angular_momentum = 0.0_dp
      DO item = 1, SIZE(group%atom)
         atom = group%atom(item)
         relative_velocity = velocity(:, atom) - center_velocity
         inertia = inertia + masses(atom)*(DOT_PRODUCT(relative(:, item), relative(:, item))*identity3() - &
                                           outer_product(relative(:, item), relative(:, item)))
         angular_momentum = angular_momentum + masses(atom)*cross(relative(:, item), relative_velocity)
      END DO
      CALL diagonalize_symmetric(inertia, eigenvalues, eigenvectors, ierr)
      IF (ierr /= 0) THEN
         message = "Inertia-tensor diagonalization did not converge"
         RETURN
      END IF
      inertia_values = eigenvalues(3:1:-1)
      principal_axes = eigenvectors(:, 3:1:-1)
      linear = SIZE(group%atom) == 2 .OR. &
               inertia_values(1) <= 1.0E-8_dp*MAX(inertia_values(3), TINY(1.0_dp))
      translation = SQRT(group_mass)*center_velocity
      rotation = 0.0_dp
      IF (linear) THEN
         rotation = angular_momentum/SQRT(0.5_dp*(inertia_values(2) + inertia_values(3)))
      ELSE
         DO axis = 1, 3
            IF (inertia_values(axis) > 1.0E-10_dp*MAX(inertia_values(3), 1.0_dp)) THEN
               rotation(axis) = DOT_PRODUCT(principal_axes(:, axis), angular_momentum)/ &
                                SQRT(inertia_values(axis))
            END IF
         END DO
      END IF
   END SUBROUTINE molecule_motion

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

      product = [a(2)*b(3) - a(3)*b(2), &
                 a(3)*b(1) - a(1)*b(3), &
                 a(1)*b(2) - a(2)*b(1)]
   END FUNCTION cross

END MODULE trajana_molecular_motion
