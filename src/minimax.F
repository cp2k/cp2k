!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate the minimax coefficients in order to
!>        approximate 1/x as a sum over exponential functions
!>        1/x ~ SUM_{i}^{K} w_i EXP(-a_i * x) for x belonging to [1:Rc].
!>
!>        This module is an extension of original minimax module minimax_k15
!>        (up to K = 15) to provide minimax approximations for larger
!>        ranges Rc (up to K = 53).
!>
!>        k53 implementation is based on directly tabulated coefficients from
!>        D. Braess and W. Hackbusch, IMA Journal of Numerical Analysis 25.4 (2005): 685-697
!>        http://www.mis.mpg.de/scicomp/EXP_SUM/1_x
!>
!>        Note: Due to discrete Rc values, the k53 implementation does not yield
!>        optimal approximations for arbitrary Rc. If optimal minimax coefficients
!>        are needed, the minimax_k15 module should be extended by interpolating
!>        k53 coefficients.
!> \par History
!>      02.2016 created [Patrick Seewald]
! **************************************************************************************************

MODULE minimax
   USE cp_log_handling,                 ONLY: cp_to_string
   USE kinds,                           ONLY: dp
   USE minimax_k15,                     ONLY: check_range_k15,&
                                              get_minimax_coeff_k15,&
                                              get_minimax_numerical_error
   USE minimax_k53,                     ONLY: R_max,&
                                              R_mm,&
                                              err_mm,&
                                              get_minimax_coeff_low,&
                                              k_max,&
                                              k_mm,&
                                              k_p,&
                                              n_approx
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'minimax'

   INTEGER, PARAMETER :: mm_k15 = 0, mm_k53 = 1

   PUBLIC :: get_minimax_coeff, minimax_validate, check_range

   ! Imported from minimax_k53:

   ! Number of tabulated minimax approximations:
   ! INTEGER, PARAMETER :: n_approx

   ! Number of K values:
   ! INTEGER, PARAMETER :: n_k

   ! Maximum K value:
   ! INTEGER, PARAMETER :: k_max

   ! Maximum range Rc:
   ! REAL(KIND=dp), PARAMETER :: R_max

   ! Values of K:
   ! INTEGER, PARAMETER, DIMENSION(n_approx) :: k_mm

   ! Values of Rc:
   ! REAL(KIND=dp), PARAMETER, DIMENSION(n_approx) :: R_mm

   ! Values of minimax error:
   ! REAL(KIND=dp), PARAMETER, DIMENSION(n_approx) :: err_mm

   ! Note: the coefficients (k_mm, R_mm, err_mm) are sorted w.r.t. 1) k_mm, 2) R_mm

   ! Given the ith value of K, k_p(i) points to the first minimax
   ! approximation with K terms:
   ! INTEGER, PARAMETER, DIMENSION(n_k+1) :: k_p

   ! Minimax coefficients aw of the ith minimax approximation:
   ! SUBROUTINE get_minimax_coeff_low(i, aw)

CONTAINS

! **************************************************************************************************
!> \brief Check that a minimax approximation is available for given input k, Rc.
!>        ierr ==  0: everything ok
!>        ierr ==  1: Rc too small
!>        ierr == -1: k too large
!> \param k ...
!> \param Rc ...
!> \param ierr ...
!> \note: ierr ==  1 is not a fatal error since get_minimax_coeff will return
!>        k53 minimax coefficients with smallest possible range.
! **************************************************************************************************
   SUBROUTINE check_range(k, Rc, ierr)
      INTEGER, INTENT(IN)                                :: k
      REAL(KIND=dp), INTENT(IN)                          :: Rc
      INTEGER, INTENT(OUT)                               :: ierr

      ierr = 0
      IF (k .LE. 15) THEN
         CALL check_range_k15(k, Rc, ierr)
      ELSE
         IF (k .GT. k_max) ierr = -1
      ENDIF

   END SUBROUTINE check_range

! **************************************************************************************************
!> \brief Get best minimax approximation for given input parameters. Automatically
!>        chooses the most exact set of minimax coefficients (k15 or k53) for
!>        given k, Rc.
!> \param k Number of minimax terms
!> \param Rc Minimax range
!> \param aw The a_i and w_i coefficient are stored in aw such that the first 1:K
!>        elements correspond to a_i and the K+1:2k correspond to w_i.
!> \param mm_error Numerical error of minimax approximation for given k, Rc
!> \param which_coeffs Whether the coefficients returned have been generated from
!>        k15 or k53 coefficients (mm_k15 or mm_k53).
! **************************************************************************************************
   SUBROUTINE get_minimax_coeff(k, Rc, aw, mm_error, which_coeffs)
      INTEGER, INTENT(IN)                                :: k
      REAL(KIND=dp), INTENT(IN)                          :: Rc
      REAL(KIND=dp), DIMENSION(2*k), INTENT(INOUT)       :: aw
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: mm_error
      INTEGER, INTENT(OUT), OPTIONAL                     :: which_coeffs

      INTEGER                                            :: ierr

      IF (k .LE. 15) THEN
         CALL check_range_k15(k, Rc, ierr)
         IF (ierr .EQ. 1) THEN ! Rc too small for k15 coeffs --> use k53
            CALL get_minimax_coeff_k53(k, Rc, aw, mm_error)
            IF (PRESENT(which_coeffs)) which_coeffs = mm_k53
         ELSE
            CPASSERT(ierr .EQ. 0)
            CALL get_minimax_coeff_k15(k, Rc, aw, mm_error)
            IF (PRESENT(which_coeffs)) which_coeffs = mm_k15
         ENDIF
      ELSEIF (k .LE. 53) THEN
         CALL get_minimax_coeff_k53(k, Rc, aw, mm_error)
         IF (PRESENT(which_coeffs)) which_coeffs = mm_k53
      ELSE
         CPABORT("No minimax approximations available for k = "//cp_to_string(k))
      ENDIF
   END SUBROUTINE get_minimax_coeff

! **************************************************************************************************
!> \brief Get minimax coefficients: k53 implementation (coefficients up to k=53 terms).
!>        All a_i and w_i for a set of discrete values Rc, k are tabulated and
!>        the most accurate coefficients for given input k, Rc are returned.
!> \param k ...
!> \param Rc ...
!> \param aw ...
!> \param mm_error ...
! **************************************************************************************************
   SUBROUTINE get_minimax_coeff_k53(k, Rc, aw, mm_error)
      INTEGER, INTENT(IN)                                :: k
      REAL(KIND=dp), INTENT(IN)                          :: Rc
      REAL(KIND=dp), DIMENSION(2*k), INTENT(INOUT)       :: aw
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: mm_error

      INTEGER                                            :: i_mm

      CALL get_best_minimax_approx_k53(k, Rc, i_mm)
      CALL get_minimax_coeff_low(i_mm, aw)
      IF (PRESENT(mm_error)) mm_error = get_minimax_numerical_error(Rc, aw)

   END SUBROUTINE get_minimax_coeff_k53

! **************************************************************************************************
!> \brief find minimax approx. with k terms that is most accurate for range Rc.
!> \param k ...
!> \param Rc ...
!> \param i_mm ...
!> \param ge_Rc Whether the tabulated range of the returned minimax approximation
!>              must be strictly greater than or equal to Rc. Default is .FALSE.
! **************************************************************************************************
   SUBROUTINE get_best_minimax_approx_k53(k, Rc, i_mm, ge_Rc)
      INTEGER, INTENT(IN)                                :: k
      REAL(KIND=dp), INTENT(IN)                          :: Rc
      INTEGER, INTENT(OUT)                               :: i_mm
      LOGICAL, OPTIONAL                                  :: ge_Rc

      INTEGER                                            :: i, i_k, i_l, i_r
      REAL(KIND=dp)                                      :: error_l, error_r, R_k_max, R_k_min
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: aw

      CPASSERT(k .LE. k_max)

      ! find k pointer and smallest and largest R_mm value for this k
      i_k = 1
      DO WHILE (k_mm(k_p(i_k)) .LT. k)
         i_k = i_k+1
      ENDDO
      CPASSERT(k_mm(k_p(i_k)) .EQ. k)

      R_k_min = R_mm(k_p(i_k))
      R_k_max = R_mm(k_p(i_k+1)-1)

      IF (Rc .GE. R_k_max) THEN
         i_mm = k_p(i_k+1)-1 ! pointer to largest Rc for current k
      ELSE IF (Rc .LE. R_k_min) THEN
         i_mm = k_p(i_k) ! pointer to smallest Rc for current k
      ELSE
         i = k_p(i_k)
         DO WHILE (Rc .GT. R_mm(i))
            i = i+1
         ENDDO
         i_r = i ! pointer to closest R_mm >= Rc
         i_l = i-1 ! pointer to closest R_mm < Rc

         IF (PRESENT(ge_Rc)) THEN
            IF (ge_Rc) THEN
               i_mm = i_r
               RETURN
            ENDIF
         ENDIF

         ALLOCATE (aw(2*k_mm(i_r)))
         CALL get_minimax_coeff_low(i_r, aw)
         error_l = get_minimax_numerical_error(Rc, aw)
         DEALLOCATE (aw)
         ALLOCATE (aw(2*k_mm(i_l)))
         CALL get_minimax_coeff_low(i_l, aw)
         error_r = get_minimax_numerical_error(Rc, aw)
         DEALLOCATE (aw)
         i_mm = MERGE(i_r, i_l, error_l .LE. error_r)
      ENDIF

   END SUBROUTINE get_best_minimax_approx_k53

! **************************************************************************************************
!> \brief Unit test checking that numerical error of minimax approximations
!>        generated using any k15 or k53 coefficients is consistent with
!>        tabulated error.
!> \param n_R Number of Rc values to be tested.
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE minimax_validate(n_R, iw)
      INTEGER, INTENT(IN)                                :: n_R, iw

      INTEGER                                            :: i_mm, i_R, ierr, k, which_coeffs
      LOGICAL                                            :: do_exit, verbose
      REAL(KIND=dp)                                      :: dR, mm_error, R, ref_error
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: aw

      verbose = .FALSE.

      IF (iw > 0) THEN
         WRITE (iw, '(//T2,A)') &
            "Unit tests for minimax 1/x ~ SUM_{i}^{K} w_i EXP(-a_i * x) for x belonging to [1:Rc]"
         WRITE (iw, '(T2,84("*"))')
      ENDIF

      IF (iw > 0) THEN
         WRITE (iw, '(/T2,A)') &
            "1) checking numerical error against tabulated error at tabulated values Rc"
         WRITE (iw, '(/T2,A)') &
            "which coeffs, K, Rc, num. error, ref. error, rel. diff (num. error - ref. error)/(ref. error)"
         WRITE (iw, '(T2,54("-"))')
      ENDIF
      DO i_mm = 1, n_approx
         R = R_mm(i_mm)
         k = k_mm(i_mm)
         CALL check_range(k, R, ierr)
         IF (ierr .EQ. 0) THEN
            ALLOCATE (aw(2*k))
            CALL get_minimax_coeff(k, R, aw, mm_error, which_coeffs)
            ref_error = err_mm(i_mm)
            DEALLOCATE (aw)
            IF (iw > 0) WRITE (iw, '(T2,A4, I3, ES10.1, ES12.3, ES12.3, ES12.3)') &
               MERGE("k15", "k53", which_coeffs .EQ. mm_k15), k, R, &
               mm_error, ref_error, (mm_error-ref_error)/ref_error
            CPASSERT(mm_error .LE. ref_error*1.05_dp+1.0E-15_dp)
         ELSE
            IF (iw > 0) WRITE (iw, '(T2,A4, I3, ES10.1, 3X, A)') "k15", k, R, "missing"
         ENDIF

         IF (k .LE. 15) THEN
            ALLOCATE (aw(2*k))
            CALL get_minimax_coeff_k53(k, R, aw, mm_error)
            ref_error = err_mm(i_mm)
            DEALLOCATE (aw)
            IF (iw > 0) WRITE (iw, '(T2,A4,I3, ES10.1, ES12.3, ES12.3, ES12.3)') &
               "k53", k, R, mm_error, ref_error, (mm_error-ref_error)/ref_error
            IF (mm_error .GT. ref_error*1.05_dp+1.0E-15_dp) THEN
               CPABORT("Test 1 failed: numerical error is larger than tabulated error")
            ENDIF
         ENDIF
      ENDDO

      IF (iw > 0 .AND. n_R .GT. 0) THEN
         WRITE (iw, '(T2,54("-"))')
         WRITE (iw, '(/T2,A)') "Test 1 OK"
         WRITE (iw, '(/T2,A)') &
            "2) checking numerical error against tabulated error at arbitrary values Rc"
         IF (verbose) THEN
            WRITE (iw, '(/T2,A)') &
               "which coeffs, K, Rc, num. error, ref. error, rel. diff (num. error - ref. error)/(ref. error)"
            WRITE (iw, '(T2,54("-"))')
         ENDIF
      ENDIF
      dR = R_max**(1.0_dp/n_R)

      DO k = 1, k_max
         ALLOCATE (aw(2*k))
         do_exit = .FALSE.
         DO i_R = 1, n_R
            R = dR**i_R
            CALL get_minimax_coeff(k, R, aw, mm_error, which_coeffs)
            CALL get_best_minimax_approx_k53(k, R, i_mm, ge_Rc=.TRUE.)
            IF (R .GT. R_mm(i_mm)) THEN
               R = R_max
               do_exit = .TRUE.
            ENDIF
            ref_error = err_mm(i_mm)
            IF (verbose .AND. iw > 0) WRITE (iw, '(T2, A4, I3, ES10.1, ES12.3, ES12.3, ES12.3)') &
               MERGE("k15", "k53", which_coeffs .EQ. mm_k15), k, R, &
               mm_error, ref_error, (mm_error-ref_error)/ref_error
            IF (mm_error .GT. ref_error*1.05_dp+1.0E-15_dp) THEN
               CPABORT("Test 2 failed: numerical error is larger than tabulated error")
            ENDIF
            IF (do_exit) EXIT
         ENDDO
         DEALLOCATE (aw)
      ENDDO
      IF (iw > 0) THEN
         IF (verbose) WRITE (iw, '(T2,54("-"))')
         WRITE (iw, '(/T2,A)') "Test 2 OK"
      ENDIF
   END SUBROUTINE minimax_validate

END MODULE minimax
