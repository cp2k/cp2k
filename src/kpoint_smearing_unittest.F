!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM kpoint_smearing_unittest
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_smearing_edge_status
   USE scf_control_types,               ONLY: smear_type

   IMPLICIT NONE

   LOGICAL                                            :: first_fractional, has_weight, last_occupied
   LOGICAL, DIMENSION(2)                              :: last_occupied_spin
   REAL(KIND=dp), DIMENSION(2)                        :: weight
   REAL(KIND=dp), DIMENSION(3, 2, 2)                  :: occupation
   TYPE(smear_type)                                   :: smear

   smear%do_smear = .TRUE.
   smear%eps_fermi_dirac = 1.0E-5_dp
   occupation = 0.0_dp
   occupation(1, :, :) = 1.0_dp
   weight = [1.0_dp, 0.0_dp]

   occupation(1, 2, 1) = 0.5_dp
   occupation(3, 2, 2) = 0.5_dp
   CALL kpoint_smearing_edge_status(occupation, weight, smear, 2, has_weight, &
                                    first_fractional, last_occupied, last_occupied_spin)
   IF (.NOT. has_weight) ERROR STOP "Positive K-point weight was not detected"
   IF (first_fractional) ERROR STOP "Zero-weight first-band occupation was included"
   IF (last_occupied) ERROR STOP "Zero-weight last-band occupation was included"
   IF (ANY(last_occupied_spin)) ERROR STOP "Zero-weight spin edge was included"

   weight(2) = 0.25_dp
   CALL kpoint_smearing_edge_status(occupation, weight, smear, 2, has_weight, &
                                    first_fractional, last_occupied, last_occupied_spin)
   IF (.NOT. first_fractional) ERROR STOP "Fractional first band was not detected"
   IF (.NOT. last_occupied) ERROR STOP "Occupied last band was not detected"
   IF (last_occupied_spin(1)) ERROR STOP "Wrong occupied spin edge was reported"
   IF (.NOT. last_occupied_spin(2)) ERROR STOP "Occupied spin edge was not reported"

   occupation = 0.0_dp
   occupation(1, :, 1) = 2.0_dp
   CALL kpoint_smearing_edge_status(occupation(:, :, 1:1), weight, smear, 1, has_weight, &
                                    first_fractional, last_occupied)
   IF (first_fractional) ERROR STOP "Restricted full occupation was classified as fractional"

   weight = 0.0_dp
   CALL kpoint_smearing_edge_status(occupation, weight, smear, 2, has_weight, &
                                    first_fractional, last_occupied, last_occupied_spin)
   IF (has_weight .OR. first_fractional .OR. last_occupied) THEN
      ERROR STOP "Zero-weight K-point set produced an edge status"
   END IF

   smear%do_smear = .FALSE.
   weight = 1.0_dp
   CALL kpoint_smearing_edge_status(occupation, weight, smear, 2, has_weight, &
                                    first_fractional, last_occupied, last_occupied_spin)
   IF (has_weight .OR. first_fractional .OR. last_occupied) THEN
      ERROR STOP "Disabled smearing produced an edge status"
   END IF

END PROGRAM kpoint_smearing_unittest
