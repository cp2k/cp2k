!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM gauxc_gradient_policy_test
   USE xc_gauxc_gradient_policy,        ONLY: gauxc_has_onedft_gradient_fix,&
                                              gauxc_resolve_gradient_policy

   IMPLICIT NONE
   LOGICAL :: fd, mpi, self, fixed
   INTEGER :: error

   fixed = gauxc_has_onedft_gradient_fix
   CALL resolve(.TRUE., .TRUE., 4, .TRUE., 2, .FALSE., "AUTO", .FALSE.)
   CALL check(fd .EQV. .NOT. fixed, "high-l fallback")
   CALL check(mpi .EQV. fixed, "high-l MPI")
   CALL check(.NOT. self .AND. error == 0, "no unused gradient runtime")

   CALL resolve(.TRUE., .TRUE., 3, .TRUE., 2, .FALSE., "AUTO", .FALSE.)
   CALL check(.NOT. fd .AND. (mpi .EQV. fixed) .AND. (self .NEQV. fixed), "low-l AUTO")
   CALL resolve(.TRUE., .TRUE., 4, .TRUE., 1, .FALSE., "AUTO", .FALSE.)
   CALL check((fd .NEQV. fixed) .AND. .NOT. mpi .AND. .NOT. self, "serial high-l")
   CALL resolve(.FALSE., .TRUE., 4, .TRUE., 2, .FALSE., "AUTO", .TRUE.)
   CALL check(.NOT. fd .AND. .NOT. mpi .AND. .NOT. self, "energy-only unchanged")

   CALL resolve(.TRUE., .TRUE., 4, .FALSE., 2, .FALSE., "AUTO", .FALSE.)
   CALL check(fd .AND. .NOT. mpi .AND. .NOT. self, "conventional guard unchanged")
   CALL resolve(.TRUE., .FALSE., 1, .TRUE., 2, .FALSE., "AUTO", .TRUE.)
   CALL check(fd .AND. .NOT. mpi .AND. .NOT. self, "explicit FD for low-l GPW")
   CALL resolve(.TRUE., .TRUE., 4, .TRUE., 2, .FALSE., "SELF", .FALSE.)
   CALL check((fd .NEQV. fixed) .AND. .NOT. mpi .AND. (self .EQV. fixed), "explicit SELF")
   CALL resolve(.TRUE., .TRUE., 4, .TRUE., 2, .TRUE., "AUTO", .FALSE.)
   CALL check(.NOT. mpi .AND. .NOT. self .AND. error == 0, "open-shell energy safeguard")

   CALL resolve(.TRUE., .TRUE., 4, .TRUE., 2, .FALSE., "MPI", .FALSE.)
   CALL check(error == MERGE(0, 1, fixed), "explicit MPI requires fix")
   CALL check(mpi .EQV. fixed, "explicit MPI selection")
   CALL resolve(.TRUE., .TRUE., 4, .TRUE., 2, .TRUE., "MPI", .FALSE.)
   CALL check(error == MERGE(2, 1, fixed), "explicit MPI cannot bypass energy SELF")
   CALL resolve(.TRUE., .TRUE., 4, .TRUE., 2, .FALSE., "MPI", .TRUE.)
   CALL check(error == MERGE(3, 1, fixed), "explicit MPI cannot silently select FD")
   CALL resolve(.TRUE., .TRUE., 4, .TRUE., 2, .FALSE., "INVALID", .FALSE.)
   CALL check(error == 4, "invalid runtime")
   PRINT *, "Gradient policy checks passed; verified fix: ", fixed
CONTAINS
! **************************************************************************************************
!> \brief Resolve one policy case into the shared test result variables.
!> \param need_gradient Whether a gradient is requested
!> \param gapw Whether GAPW is used
!> \param max_l Largest basis angular momentum
!> \param model Whether a OneDFT model is used
!> \param nranks Number of MPI ranks
!> \param energy_self Whether energy/VXC uses a rank-local runtime
!> \param runtime Requested gradient runtime
!> \param force_fd Explicit diagnostic finite-difference request
! **************************************************************************************************
   SUBROUTINE resolve(need_gradient, gapw, max_l, model, nranks, energy_self, runtime, force_fd)
      LOGICAL, INTENT(IN)                                :: need_gradient, gapw
      INTEGER, INTENT(IN)                                :: max_l
      LOGICAL, INTENT(IN)                                :: model
      INTEGER, INTENT(IN)                                :: nranks
      LOGICAL, INTENT(IN)                                :: energy_self
      CHARACTER(LEN=*), INTENT(IN)                       :: runtime
      LOGICAL, INTENT(IN)                                :: force_fd

      CALL gauxc_resolve_gradient_policy(need_gradient, gapw, max_l, model, nranks, energy_self, &
                                         runtime, force_fd, fd, mpi, self, error)
   END SUBROUTINE resolve

! **************************************************************************************************
!> \brief Stop the test with the failing case label.
!> \param condition Assertion to check
!> \param label Description of the assertion
! **************************************************************************************************
   SUBROUTINE check(condition, label)
      LOGICAL, INTENT(IN)                                :: condition
      CHARACTER(LEN=*), INTENT(IN)                       :: label

      IF (.NOT. condition) THEN
         PRINT *, "FAILED: ", label
         STOP 1
      END IF
   END SUBROUTINE check
END PROGRAM gauxc_gradient_policy_test
