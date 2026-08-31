!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM kpoint_smearing_unittest
   USE bibliography,                    ONLY: add_all_references
   USE input_constants,                 ONLY: smear_fermi_dirac,&
                                              smear_gaussian,&
                                              smear_mp,&
                                              smear_mv
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_smearing_edge_status
   USE qs_ot,                           ONLY: qs_ot_fixed_n_energy_hessian
   USE reference_manager,               ONLY: remove_all_references
   USE scf_control_types,               ONLY: smear_type
   USE smearing_utils,                  ONLY: Smearkp,&
                                              smearing_response_weight

   IMPLICIT NONE

   LOGICAL                                            :: first_fractional, has_weight, last_occupied
   LOGICAL, DIMENSION(2)                              :: last_occupied_spin
   REAL(KIND=dp), DIMENSION(2)                        :: weight
   REAL(KIND=dp), DIMENSION(3, 2, 2)                  :: occupation
   TYPE(smear_type)                                   :: smear

   CALL add_all_references()
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

   CALL test_fixed_n_smearing_response()
   CALL remove_all_references()

CONTAINS

! **************************************************************************************************
!> \brief Check the fixed-N occupation response for every supported K-point smearing method.
! **************************************************************************************************
   SUBROUTINE test_fixed_n_smearing_response()
      INTEGER, PARAMETER                                 :: nkpoint = 2, nstate = 5
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp, maxocc = 2.0_dp, &
                                                            sigma = 0.15_dp, TARGET = 4.3_dp

      INTEGER                                            :: ikpoint, imethod
      INTEGER, DIMENSION(4)                              :: methods
      REAL(KIND=dp)                                      :: error, kTS, mu, response_sum
      REAL(KIND=dp), DIMENSION(nkpoint)                  :: kpoint_weight
      REAL(KIND=dp), DIMENSION(nstate*nkpoint)           :: fd_response, flat_direction, &
                                                            predicted_response
      REAL(KIND=dp), &
         DIMENSION(nstate*nkpoint, nstate*nkpoint)       :: hessian
      REAL(KIND=dp), DIMENSION(nstate, nkpoint) :: direction, eigenvalue, eigenvalue_minus, &
         eigenvalue_plus, occupation, occupation_minus, occupation_plus, response

      methods(:) = [smear_fermi_dirac, smear_gaussian, smear_mp, smear_mv]
      kpoint_weight(:) = [0.4_dp, 0.6_dp]
      eigenvalue(:, 1) = [-0.38_dp, -0.14_dp, 0.02_dp, 0.23_dp, 0.52_dp]
      eigenvalue(:, 2) = [-0.31_dp, -0.08_dp, 0.07_dp, 0.28_dp, 0.61_dp]
      direction(:, 1) = [0.13_dp, -0.07_dp, 0.11_dp, -0.04_dp, 0.09_dp]
      direction(:, 2) = [-0.08_dp, 0.05_dp, -0.12_dp, 0.06_dp, -0.03_dp]
      flat_direction(:) = RESHAPE(direction, [nstate*nkpoint])

      DO imethod = 1, SIZE(methods)
         CALL Smearkp(occupation, mu, kTS, eigenvalue, TARGET, kpoint_weight, sigma, maxocc, &
                      methods(imethod))
         DO ikpoint = 1, nkpoint
            CALL smearing_response_weight( &
               response(:, ikpoint), occupation(:, ikpoint), eigenvalue(:, ikpoint), mu, sigma, &
               maxocc, nstate, methods(imethod))
            response(:, ikpoint) = kpoint_weight(ikpoint)*response(:, ikpoint)
         END DO
         response_sum = SUM(response)
         IF (ABS(response_sum) <= EPSILON(response_sum)) THEN
            ERROR STOP "Singular fixed-N smearing response in unit test"
         END IF
         IF ((methods(imethod) == smear_mp .OR. methods(imethod) == smear_mv) .AND. &
             .NOT. ANY(response < 0.0_dp)) THEN
            ERROR STOP "Signed smearing-response branch was not exercised"
         END IF
         CALL qs_ot_fixed_n_energy_hessian( &
            RESHAPE(response, [nstate*nkpoint]), response_sum, hessian)

         eigenvalue_plus(:, :) = eigenvalue + fd_step*direction
         eigenvalue_minus(:, :) = eigenvalue - fd_step*direction
         CALL Smearkp(occupation_plus, mu, kTS, eigenvalue_plus, TARGET, kpoint_weight, sigma, &
                      maxocc, methods(imethod))
         CALL Smearkp(occupation_minus, mu, kTS, eigenvalue_minus, TARGET, kpoint_weight, sigma, &
                      maxocc, methods(imethod))
         DO ikpoint = 1, nkpoint
            fd_response((ikpoint - 1)*nstate + 1:ikpoint*nstate) = &
               kpoint_weight(ikpoint)*(occupation_plus(:, ikpoint) - &
                                       occupation_minus(:, ikpoint))/(2.0_dp*fd_step)
         END DO
         predicted_response(:) = -MATMUL(hessian, flat_direction)
         error = MAXVAL(ABS(fd_response - predicted_response))
         error = MAX(error, ABS(SUM(fd_response)))
         error = MAX(error, MAXVAL(ABS(MATMUL(hessian, SPREAD(1.0_dp, 1, nstate*nkpoint)))))
         IF (error > 2.0E-6_dp) ERROR STOP "Fixed-N smearing response finite difference failed"
      END DO

   END SUBROUTINE test_fixed_n_smearing_response

END PROGRAM kpoint_smearing_unittest
