!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE xc_gauxc_gradient_policy
   IMPLICIT NONE
   PRIVATE

#if defined(__GAUXC) && defined(__GAUXC_ONEDFT_GRADIENT_FIX)
   LOGICAL, PARAMETER :: gauxc_has_onedft_gradient_fix = .TRUE.
#else
   LOGICAL, PARAMETER :: gauxc_has_onedft_gradient_fix = .FALSE.
#endif

   PUBLIC :: gauxc_has_onedft_gradient_fix, gauxc_resolve_gradient_policy

CONTAINS

! **************************************************************************************************
!> \brief Select molecular XC gradients without changing the energy/VXC runtime.
!> \param need_gradient Whether forces or a molecular virial are requested
!> \param gapw_method Whether the calculation uses GAPW or GAPW_XC
!> \param max_l Largest angular momentum in the basis
!> \param is_model Whether a OneDFT model is used
!> \param nranks Size of the force-evaluation communicator
!> \param energy_self Whether energy/VXC uses replicated rank-local runtimes
!> \param runtime Requested AUTO, SELF or MPI gradient runtime
!> \param force_fd Force finite differences for diagnostics
!> \param use_fd Selected finite-difference route
!> \param use_mpi Selected distributed analytical route
!> \param use_self Whether a separate rank-local gradient runtime is required
!> \param error Zero on success; nonzero for an unsupported explicit request
! **************************************************************************************************
   PURE SUBROUTINE gauxc_resolve_gradient_policy(need_gradient, gapw_method, max_l, is_model, nranks, &
                                                 energy_self, runtime, force_fd, use_fd, use_mpi, use_self, error)
      LOGICAL, INTENT(IN)                                :: need_gradient, gapw_method
      INTEGER, INTENT(IN)                                :: max_l
      LOGICAL, INTENT(IN)                                :: is_model
      INTEGER, INTENT(IN)                                :: nranks
      LOGICAL, INTENT(IN)                                :: energy_self
      CHARACTER(LEN=*), INTENT(IN)                       :: runtime
      LOGICAL, INTENT(IN)                                :: force_fd
      LOGICAL, INTENT(OUT)                               :: use_fd, use_mpi, use_self
      INTEGER, INTENT(OUT)                               :: error

      use_fd = .FALSE.
      use_mpi = .FALSE.
      use_self = .FALSE.
      error = 0
      IF (.NOT. need_gradient) RETURN

      ! PR #222 fixes OneDFT, not the separate conventional-functional path.
      use_fd = force_fd .OR. (gapw_method .AND. max_l > 3 .AND. &
                              .NOT. (is_model .AND. gauxc_has_onedft_gradient_fix))
      IF (.NOT. is_model .OR. nranks <= 1) RETURN

      SELECT CASE (runtime)
      CASE ("AUTO")
         use_mpi = gauxc_has_onedft_gradient_fix .AND. .NOT. energy_self .AND. .NOT. use_fd
         use_self = .NOT. energy_self .AND. .NOT. use_mpi .AND. .NOT. use_fd
      CASE ("SELF")
         use_self = .NOT. energy_self .AND. .NOT. use_fd
      CASE ("MPI")
         IF (.NOT. gauxc_has_onedft_gradient_fix) THEN
            error = 1
         ELSE IF (energy_self) THEN
            error = 2
         ELSE IF (use_fd) THEN
            error = 3
         ELSE
            use_mpi = .TRUE.
         END IF
      CASE DEFAULT
         error = 4
      END SELECT
   END SUBROUTINE gauxc_resolve_gradient_policy

END MODULE xc_gauxc_gradient_policy
