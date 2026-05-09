!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief CP2K-side tblite-compatible SCC Broyden mixer.
! **************************************************************************************************
MODULE tblite_scc_mixer
#if defined(__TBLITE)
   USE mctc_env, ONLY: error_type, &
                       fatal_error, &
                       wp
   USE tblite_lapack, ONLY: getrf, &
                            getrs
   USE tblite_scf_mixer_type, ONLY: mixer_type
#endif

   IMPLICIT NONE

   PRIVATE

#if defined(__TBLITE)
   PUBLIC :: new_cp2k_tblite_mixer

   TYPE, EXTENDS(mixer_type) :: cp2k_tblite_broyden_mixer_type
      INTEGER                                            :: idif = 0
      INTEGER                                            :: iget = 0
      INTEGER                                            :: iset = 0
      INTEGER                                            :: iter = 0
      INTEGER                                            :: memory = 0
      INTEGER                                            :: ndim = 0
      REAL(wp)                                           :: damp = 0.0_wp
      REAL(wp)                                           :: max_weight = 0.0_wp
      REAL(wp)                                           :: min_weight = 0.0_wp
      REAL(wp)                                           :: omega0 = 0.0_wp
      REAL(wp)                                           :: weight_factor = 0.0_wp
      REAL(wp), ALLOCATABLE, DIMENSION(:, :)             :: a, df, u
      REAL(wp), ALLOCATABLE, DIMENSION(:)                :: dq, dqlast, omega, q_in, qlast_in
   CONTAINS
      PROCEDURE                                          :: diff_1d => cp2k_tblite_mixer_diff_1d
      PROCEDURE                                          :: get_1d => cp2k_tblite_mixer_get_1d
      PROCEDURE                                          :: get_error => cp2k_tblite_mixer_get_error
      PROCEDURE                                          :: next => cp2k_tblite_mixer_next
      PROCEDURE                                          :: set_1d => cp2k_tblite_mixer_set_1d
   END TYPE cp2k_tblite_broyden_mixer_type
#endif

CONTAINS

#if defined(__TBLITE)
! **************************************************************************************************
!> \brief Create a CP2K-side tblite-compatible Broyden mixer.
!> \param mixer mixer instance
!> \param memory Broyden history length
!> \param ndim number of SCC variables
!> \param damping first-step damping
!> \param omega0 Broyden regularization weight
!> \param min_weight minimum dynamic Broyden weight
!> \param max_weight maximum dynamic Broyden weight
!> \param weight_factor residual-to-weight scaling factor
! **************************************************************************************************
   SUBROUTINE new_cp2k_tblite_mixer(mixer, memory, ndim, damping, omega0, min_weight, max_weight, &
                                    weight_factor)
      CLASS(mixer_type), ALLOCATABLE, INTENT(OUT)        :: mixer
      INTEGER, INTENT(IN)                                :: memory, ndim
      REAL(wp), INTENT(IN)                               :: damping, omega0, min_weight, max_weight, &
                                                            weight_factor

      TYPE(cp2k_tblite_broyden_mixer_type), ALLOCATABLE  :: broyden

      ALLOCATE (broyden)
      broyden%ndim = ndim
      broyden%memory = memory
      broyden%damp = damping
      broyden%omega0 = omega0
      broyden%min_weight = min_weight
      broyden%max_weight = max_weight
      broyden%weight_factor = weight_factor
      ALLOCATE (broyden%a(memory, memory))
      ALLOCATE (broyden%df(ndim, memory))
      ALLOCATE (broyden%u(ndim, memory))
      ALLOCATE (broyden%dq(ndim))
      ALLOCATE (broyden%dqlast(ndim))
      ALLOCATE (broyden%omega(memory))
      ALLOCATE (broyden%q_in(ndim))
      ALLOCATE (broyden%qlast_in(ndim))
      CALL MOVE_ALLOC(broyden, mixer)

   END SUBROUTINE new_cp2k_tblite_mixer

! **************************************************************************************************
!> \brief Set input SCC variables.
!> \param self mixer
!> \param qvec SCC-variable vector
! **************************************************************************************************
   SUBROUTINE cp2k_tblite_mixer_set_1d(self, qvec)
      CLASS(cp2k_tblite_broyden_mixer_type), INTENT(INOUT) :: self
      REAL(wp), DIMENSION(:), INTENT(IN)                  :: qvec

      self%q_in(self%iset + 1:self%iset + SIZE(qvec)) = qvec
      self%iset = self%iset + SIZE(qvec)

   END SUBROUTINE cp2k_tblite_mixer_set_1d

! **************************************************************************************************
!> \brief Set SCC residual.
!> \param self mixer
!> \param qvec output SCC-variable vector before mixing
! **************************************************************************************************
   SUBROUTINE cp2k_tblite_mixer_diff_1d(self, qvec)
      CLASS(cp2k_tblite_broyden_mixer_type), INTENT(INOUT) :: self
      REAL(wp), DIMENSION(:), INTENT(IN)                  :: qvec

      self%dq(self%idif + 1:self%idif + SIZE(qvec)) = &
         qvec - self%q_in(self%idif + 1:self%idif + SIZE(qvec))
      self%idif = self%idif + SIZE(qvec)

   END SUBROUTINE cp2k_tblite_mixer_diff_1d

! **************************************************************************************************
!> \brief Apply the next Broyden update.
!> \param self mixer
!> \param error error handling
! **************************************************************************************************
   SUBROUTINE cp2k_tblite_mixer_next(self, error)
      CLASS(cp2k_tblite_broyden_mixer_type), INTENT(INOUT) :: self
      TYPE(error_type), ALLOCATABLE, INTENT(OUT)          :: error

      INTEGER                                            :: info

      self%iset = 0
      self%idif = 0
      self%iget = 0
      self%iter = self%iter + 1
      CALL cp2k_tblite_broyden(self%ndim, self%q_in, self%qlast_in, self%dq, self%dqlast, &
                               self%iter, self%memory, self%damp, self%omega0, &
                               self%min_weight, self%max_weight, self%weight_factor, &
                               self%omega, self%df, self%u, self%a, info)
      IF (info /= 0) CALL fatal_error(error, "Broyden mixing failed to obtain next iteration")

   END SUBROUTINE cp2k_tblite_mixer_next

! **************************************************************************************************
!> \brief Get mixed SCC variables.
!> \param self mixer
!> \param qvec mixed SCC-variable vector
! **************************************************************************************************
   SUBROUTINE cp2k_tblite_mixer_get_1d(self, qvec)
      CLASS(cp2k_tblite_broyden_mixer_type), INTENT(INOUT) :: self
      REAL(wp), DIMENSION(:), INTENT(OUT)                 :: qvec

      qvec(:) = self%q_in(self%iget + 1:self%iget + SIZE(qvec))
      self%iget = self%iget + SIZE(qvec)

   END SUBROUTINE cp2k_tblite_mixer_get_1d

! **************************************************************************************************
!> \brief Get RMS SCC residual.
!> \param self mixer
!> \return residual
! **************************************************************************************************
   PURE FUNCTION cp2k_tblite_mixer_get_error(self) RESULT(error)
      CLASS(cp2k_tblite_broyden_mixer_type), INTENT(IN)  :: self
      REAL(wp)                                           :: error

      error = SQRT(SUM(self%dq**2)/SIZE(self%dq))

   END FUNCTION cp2k_tblite_mixer_get_error

! **************************************************************************************************
!> \brief Modified Broyden update matching tblite's default algorithm, with explicit constants.
!> \param n ...
!> \param q ...
!> \param qlast ...
!> \param dq ...
!> \param dqlast ...
!> \param iter ...
!> \param memory ...
!> \param alpha ...
!> \param omega0 ...
!> \param minw ...
!> \param maxw ...
!> \param wfac ...
!> \param omega ...
!> \param df ...
!> \param u ...
!> \param a ...
!> \param info ...
! **************************************************************************************************
   SUBROUTINE cp2k_tblite_broyden(n, q, qlast, dq, dqlast, iter, memory, alpha, omega0, minw, &
                                  maxw, wfac, omega, df, u, a, info)
      INTEGER, INTENT(IN)                                :: n
      REAL(wp), DIMENSION(n), INTENT(INOUT)              :: q, qlast
      REAL(wp), DIMENSION(n), INTENT(IN)                 :: dq
      REAL(wp), DIMENSION(n), INTENT(INOUT)              :: dqlast
      INTEGER, INTENT(IN)                                :: iter, memory
      REAL(wp), INTENT(IN)                               :: alpha, omega0, minw, maxw, wfac
      REAL(wp), DIMENSION(memory), INTENT(INOUT)         :: omega
      REAL(wp), DIMENSION(n, memory), INTENT(INOUT)      :: df, u
      REAL(wp), DIMENSION(memory, memory), INTENT(INOUT) :: a
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: i, it1, itn, j
      REAL(wp)                                           :: inv
      REAL(wp), ALLOCATABLE, DIMENSION(:, :)             :: beta, c

      info = 0
      itn = iter - 1
      it1 = MOD(itn - 1, memory) + 1

      IF (iter == 1) THEN
         dqlast(:) = dq
         qlast(:) = q
         q(:) = q + alpha*dq
         RETURN
      END IF

      ALLOCATE (beta(MIN(memory, itn), MIN(memory, itn)), c(MIN(memory, itn), 1))

      omega(it1) = SQRT(DOT_PRODUCT(dq, dq))
      IF (omega(it1) > (wfac/maxw)) THEN
         omega(it1) = wfac/omega(it1)
      ELSE
         omega(it1) = maxw
      END IF
      IF (omega(it1) < minw) omega(it1) = minw

      df(:, it1) = dq - dqlast
      inv = MAX(SQRT(DOT_PRODUCT(df(:, it1), df(:, it1))), EPSILON(1.0_wp))
      inv = 1.0_wp/inv
      df(:, it1) = inv*df(:, it1)

      DO j = MAX(1, itn - memory + 1), itn
         i = MOD(j - 1, memory) + 1
         a(i, it1) = DOT_PRODUCT(df(:, i), df(:, it1))
         a(it1, i) = a(i, it1)
         c(i, 1) = omega(i)*DOT_PRODUCT(df(:, i), dq)
      END DO

      DO j = MAX(1, itn - memory + 1), itn
         i = MOD(j - 1, memory) + 1
         beta(:it1, i) = omega(:it1)*omega(i)*a(:it1, i)
         beta(i, i) = beta(i, i) + omega0*omega0
      END DO

      CALL cp2k_tblite_lineq(beta, c, info)
      IF (info /= 0) RETURN

      u(:, it1) = alpha*df(:, it1) + inv*(q - qlast)
      dqlast(:) = dq
      qlast(:) = q
      q(:) = q + alpha*dq

      DO j = MAX(1, itn - memory + 1), itn
         i = MOD(j - 1, memory) + 1
         q(:) = q - omega(i)*c(i, 1)*u(:, i)
      END DO

   END SUBROUTINE cp2k_tblite_broyden

! **************************************************************************************************
!> \brief Solve the small Broyden linear system.
!> \param a ...
!> \param c ...
!> \param info ...
! **************************************************************************************************
   SUBROUTINE cp2k_tblite_lineq(a, c, info)
      REAL(wp), DIMENSION(:, :), INTENT(INOUT)           :: a, c
      INTEGER, INTENT(OUT)                               :: info

      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ipiv

      ALLOCATE (ipiv(SIZE(a, 1)))
      CALL getrf(a, ipiv, info)
      IF (info == 0) CALL getrs(a, c, ipiv, info, trans="t")

   END SUBROUTINE cp2k_tblite_lineq
#endif

END MODULE tblite_scc_mixer
