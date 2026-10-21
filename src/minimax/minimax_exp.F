!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate the minimax coefficients in order to
!>        approximate 1/x as a sum over exponential functions
!>        1/x ~ SUM_{i}^{K} w_i EXP(-a_i * x) for x belonging to [1:Rc].
!>
!>        Coefficients for K = 1,...,53 are obtained by local interpolation in
!>        log(Rc) from a sparse subset of the Braess-Hackbusch minimax data. The
!>        sparse representation is generated and certified against the tabulated
!>        minimax-error envelope by tools/minimax_tools/minimax_to_fortran_source.py.
!>        The same continuous backend supplies the extended ranges used by GW grids.
!> \par History
!>      02.2016 created [Patrick Seewald]
! **************************************************************************************************

MODULE minimax_exp
   USE cp_log_handling,                 ONLY: cp_to_string
   USE kinds,                           ONLY: dp
   USE minimax_exp_data,                ONLY: &
        R_max, R_mm, R_ref, err_ref, get_minimax_coeff_low, interp_degree, k_max, k_p, ref_p
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'minimax_exp'

   PUBLIC :: get_exp_minimax_coeff, get_exp_minimax_coeff_gw, validate_exp_minimax, &
             check_exp_minimax_range

   ! Sparse coefficient and reference-error data are provided by minimax_exp_data.

CONTAINS

! **************************************************************************************************
!> \brief Check that an exponential minimax approximation is available.
!> \param k Number of exponential terms.
!> \param Rc Minimax range.
!> \param ierr Zero when supported, -1 for unsupported k, 1 for invalid Rc.
! **************************************************************************************************
   SUBROUTINE check_exp_minimax_range(k, Rc, ierr)
      INTEGER, INTENT(IN)                                :: k
      REAL(KIND=dp), INTENT(IN)                          :: Rc
      INTEGER, INTENT(OUT)                               :: ierr

      ierr = 0
      IF (k < 1 .OR. k > k_max) THEN
         ierr = -1
      ELSE IF (Rc < 1.0_dp) THEN
         ierr = 1
      END IF

   END SUBROUTINE check_exp_minimax_range

! **************************************************************************************************
!> \brief Sample the numerical error of an exponential minimax approximation.
!> \param Rc Minimax range.
!> \param aw Exponents followed by weights.
!> \return Maximum sampled absolute error.
! **************************************************************************************************
   FUNCTION get_minimax_numerical_error(Rc, aw) RESULT(numerr)
      REAL(KIND=dp), INTENT(IN)                          :: Rc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: aw
      REAL(KIND=dp)                                      :: numerr

      INTEGER, PARAMETER                                 :: npoints = 100

      INTEGER                                            :: j, k, np
      REAL(KIND=dp)                                      :: mm, x

      k = SIZE(aw)/2
      numerr = 0.0_dp
      DO np = 1, npoints
         x = 1.0_dp + (Rc - 1.0_dp)*(npoints - np)/(npoints - 1)
         mm = 0.0_dp
         DO j = 1, k
            mm = mm + aw(j + k)*EXP(-aw(j)*x)
         END DO
         numerr = MAX(numerr, ABS(mm - 1.0_dp/x))
      END DO

   END FUNCTION get_minimax_numerical_error

! **************************************************************************************************
!> \brief Get exponential minimax coefficients for arbitrary Rc.
!> \param k Number of exponential terms.
!> \param Rc Minimax range.
!> \param aw Exponents followed by weights.
!> \param mm_error Numerical approximation error at Rc.
! **************************************************************************************************
   SUBROUTINE get_exp_minimax_coeff(k, Rc, aw, mm_error)
      INTEGER, INTENT(IN)                                :: k
      REAL(KIND=dp), INTENT(IN)                          :: Rc
      REAL(KIND=dp), DIMENSION(2*k), INTENT(OUT)         :: aw
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: mm_error

      IF (k < 1 .OR. k > k_max) THEN
         CPABORT("No minimax approximations available for k = "//cp_to_string(k))
      END IF
      IF (Rc < 1.0_dp) THEN
         CPABORT("Minimax range Rc must be greater than or equal to one")
      END IF

      CALL get_minimax_coeff_data(k, Rc, aw, mm_error)

   END SUBROUTINE get_exp_minimax_coeff

! **************************************************************************************************
!> \brief Interpolate sparse exponential-minimax coefficients for arbitrary Rc.
!> \param k Number of exponential terms.
!> \param Rc Minimax range.
!> \param aw Exponents followed by weights.
!> \param mm_error Numerical approximation error at Rc.
! **************************************************************************************************
   SUBROUTINE get_minimax_coeff_data(k, Rc, aw, mm_error)
      INTEGER, INTENT(IN)                                :: k
      REAL(KIND=dp), INTENT(IN)                          :: Rc
      REAL(KIND=dp), DIMENSION(2*k), INTENT(OUT)         :: aw
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: mm_error

      INTEGER                                            :: degree, i, i_first, i_last, i_right, j, &
                                                            n_full, n_interp, n_knots, start
      REAL(KIND=dp)                                      :: error_l, error_r, weight, x
      REAL(KIND=dp), DIMENSION(2*k)                      :: aw_i, aw_l, aw_r
      REAL(KIND=dp), DIMENSION(8)                        :: nodes

      CPASSERT(k >= 1 .AND. k <= k_max)
      i_first = k_p(k)
      i_last = k_p(k + 1) - 1
      n_knots = i_last - i_first + 1

      IF (Rc <= R_mm(i_first)) THEN
         CALL get_minimax_coeff_low(i_first, aw)
      ELSE IF (Rc >= R_mm(i_last)) THEN
         CALL get_minimax_coeff_low(i_last, aw)
      ELSE IF (n_knots == 2) THEN
         CALL get_minimax_coeff_low(i_first, aw_l)
         CALL get_minimax_coeff_low(i_last, aw_r)
         error_l = get_minimax_numerical_error(Rc, aw_l)
         error_r = get_minimax_numerical_error(Rc, aw_r)
         IF (error_r <= error_l) THEN
            aw = aw_r
         ELSE
            aw = aw_l
         END IF
      ELSE
         i_right = i_first + 1
         DO WHILE (Rc > R_mm(i_right))
            i_right = i_right + 1
         END DO

         degree = interp_degree(k)
         n_full = degree + 1
         IF (n_knots <= n_full) THEN
            n_interp = n_knots
            start = i_first
         ELSE IF (i_right == i_first + 1) THEN
            n_interp = degree
            start = i_first
         ELSE IF (i_right == i_last) THEN
            n_interp = degree
            start = i_last - n_interp + 1
         ELSE
            n_interp = n_full
            start = i_right - n_interp/2
         END IF
         start = MAX(i_first, MIN(start, i_last - n_interp + 1))

         x = LOG(Rc)
         DO i = 1, n_interp
            nodes(i) = LOG(R_mm(start + i - 1))
         END DO

         aw = 0.0_dp
         DO i = 1, n_interp
            weight = 1.0_dp
            DO j = 1, n_interp
               IF (i /= j) weight = weight*(x - nodes(j))/(nodes(i) - nodes(j))
            END DO
            CALL get_minimax_coeff_low(start + i - 1, aw_i)
            aw = aw + weight*aw_i
         END DO
      END IF

      IF (PRESENT(mm_error)) mm_error = get_minimax_numerical_error(Rc, aw)

   END SUBROUTINE get_minimax_coeff_data

! **************************************************************************************************
!> \brief Find the first tabulated reference range greater than or equal to Rc.
!> \param k Number of exponential terms.
!> \param Rc Minimax range.
!> \param i_ref Global reference-data index.
! **************************************************************************************************
   SUBROUTINE get_reference_index(k, Rc, i_ref)
      INTEGER, INTENT(IN)                                :: k
      REAL(KIND=dp), INTENT(IN)                          :: Rc
      INTEGER, INTENT(OUT)                               :: i_ref

      INTEGER                                            :: i_first, i_last

      CPASSERT(k >= 1 .AND. k <= k_max)
      i_first = ref_p(k)
      i_last = ref_p(k + 1) - 1

      IF (Rc <= R_ref(i_first)) THEN
         i_ref = i_first
      ELSE IF (Rc >= R_ref(i_last)) THEN
         i_ref = i_last
      ELSE
         i_ref = i_first + 1
         DO WHILE (Rc > R_ref(i_ref))
            i_ref = i_ref + 1
         END DO
      END IF

   END SUBROUTINE get_reference_index

! **************************************************************************************************
!> \brief Get exponential-minimax coefficients for GW time/frequency grids.
!>        Above the historical lower cutoff, evaluate the continuous minimax
!>        interpolation directly at E_range. Below it, preserve the legacy GW
!>        rescaling policy.
!> \param k Number of exponential terms; supported values are 26, 28, 30, 32, 34.
!> \param E_range Energy range represented by the grid.
!> \param aw Exponents followed by weights.
! **************************************************************************************************
   SUBROUTINE get_exp_minimax_coeff_gw(k, E_range, aw)
      INTEGER, INTENT(IN)                                :: k
      REAL(KIND=dp), INTENT(IN)                          :: E_range
      REAL(KIND=dp), DIMENSION(2*k), INTENT(OUT)         :: aw

      REAL(KIND=dp)                                      :: E_ratio, R_min, Rc

      SELECT CASE (k)
      CASE (26)
         R_min = 799.0_dp
      CASE (28)
         R_min = 1545.0_dp
      CASE (30)
         R_min = 2906.0_dp
      CASE (32)
         R_min = 4862.0_dp
      CASE (34)
         R_min = 9649.0_dp
      CASE DEFAULT
         CPABORT("The grid size you choose is not available.")
      END SELECT

      E_ratio = 1.0_dp
      Rc = E_range
      IF (E_range < R_min) THEN
         E_ratio = R_min/E_range
         Rc = R_min
      END IF

      CALL get_minimax_coeff_data(k, Rc, aw)
      CALL rescale_grid(aw, E_ratio)

   END SUBROUTINE get_exp_minimax_coeff_gw

! **************************************************************************************************
!> \brief Preserve the historical low-range scaling used by GW minimax grids.
!> \param aw Exponents followed by weights.
!> \param E_ratio Ratio between the minimum supported range and requested range.
! **************************************************************************************************
   PURE SUBROUTINE rescale_grid(aw, E_ratio)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: aw
      REAL(KIND=dp), INTENT(INOUT)                       :: E_ratio

      IF (E_ratio > 1.5_dp) THEN
         E_ratio = E_ratio/1.5_dp
      ELSE
         E_ratio = 1.0_dp
      END IF
      aw(:) = aw(:)*E_ratio

   END SUBROUTINE rescale_grid
! **************************************************************************************************
!> \brief Unit test checking exponential-minimax errors against the tabulated reference envelope.
!> \param n_R Number of arbitrary Rc values to be tested.
!> \param iw Output unit.
! **************************************************************************************************
   SUBROUTINE validate_exp_minimax(n_R, iw)
      INTEGER, INTENT(IN)                                :: n_R, iw

      INTEGER                                            :: i_R, i_ref, k
      LOGICAL                                            :: do_exit
      REAL(KIND=dp)                                      :: dR, mm_error, R, ref_error
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: aw

      IF (iw > 0) THEN
         WRITE (iw, '(//T2,A)') &
            "Unit tests for minimax 1/x ~ SUM_{i}^{K} w_i EXP(-a_i * x) for x belonging to [1:Rc]"
         WRITE (iw, '(T2,84("*"))')
         WRITE (iw, '(/T2,A)') &
            "1) checking numerical error against tabulated error at tabulated values Rc"
         WRITE (iw, '(/T2,A)') &
            "K, Rc, num. error, ref. error, rel. diff (num. error - ref. error)/(ref. error)"
         WRITE (iw, '(T2,54("-"))')
      END IF

      DO k = 1, k_max
         ALLOCATE (aw(2*k))
         DO i_ref = ref_p(k), ref_p(k + 1) - 1
            R = R_ref(i_ref)
            CALL get_exp_minimax_coeff(k, R, aw, mm_error)
            ref_error = err_ref(i_ref)
            IF (iw > 0) WRITE (iw, '(T2,I3, ES10.1, ES12.3, ES12.3, ES12.3)') &
               k, R, mm_error, ref_error, (mm_error - ref_error)/ref_error
            IF (mm_error > ref_error*1.05_dp + 1.0E-15_dp) THEN
               CPABORT("Test 1 failed: numerical error is larger than tabulated error")
            END IF
         END DO
         DEALLOCATE (aw)
      END DO

      IF (iw > 0) THEN
         WRITE (iw, '(T2,54("-"))')
         WRITE (iw, '(/T2,A)') "Test 1 OK"
      END IF
      IF (n_R <= 0) RETURN

      IF (iw > 0) THEN
         WRITE (iw, '(/T2,A)') &
            "2) checking numerical error against tabulated error at arbitrary values Rc"
         WRITE (iw, '(/T2,A)') &
            "K, Rc, num. error, ref. error, rel. diff (num. error - ref. error)/(ref. error)"
         WRITE (iw, '(T2,54("-"))')
      END IF

      dR = R_max**(1.0_dp/n_R)
      DO k = 1, k_max
         ALLOCATE (aw(2*k))
         do_exit = .FALSE.
         DO i_R = 1, n_R
            R = dR**i_R
            CALL get_exp_minimax_coeff(k, R, aw, mm_error)
            CALL get_reference_index(k, R, i_ref)
            IF (R > R_ref(i_ref)) THEN
               R = R_max
               do_exit = .TRUE.
            END IF
            ref_error = err_ref(i_ref)
            IF (iw > 0) WRITE (iw, '(T2,I3, ES10.1, ES12.3, ES12.3, ES12.3)') &
               k, R, mm_error, ref_error, (mm_error - ref_error)/ref_error
            IF (mm_error > ref_error*1.05_dp + 1.0E-15_dp) THEN
               CPABORT("Test 2 failed: numerical error is larger than tabulated error")
            END IF
            IF (do_exit) EXIT
         END DO
         DEALLOCATE (aw)
      END DO
      IF (iw > 0) THEN
         WRITE (iw, '(T2,54("-"))')
         WRITE (iw, '(/T2,A)') "Test 2 OK"
      END IF

   END SUBROUTINE validate_exp_minimax

END MODULE minimax_exp
