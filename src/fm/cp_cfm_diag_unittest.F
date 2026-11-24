!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM cp_cfm_diag_unittest
   USE ISO_C_BINDING,                   ONLY: C_ASSOCIATED
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create, cp_blacs_env_release, cp_blacs_env_type
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig, cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create, cp_cfm_get_submatrix, cp_cfm_release, &
                                              cp_cfm_set_submatrix, cp_cfm_type
   USE cp_fm_diag,                      ONLY: FM_DIAG_TYPE_CUSOLVER, FM_DIAG_TYPE_DEFAULT, FM_DIAG_TYPE_SCALAPACK, cusolver_n_min, &
                                              diag_type, direct_generalized_diagonalization
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create, cp_fm_struct_release, cp_fm_struct_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: cp2k_is_parallel, mp_para_env_type, mp_world_finalize, mp_world_init
#if defined(__CUSOLVERMP)
   USE offload_api,                     ONLY: offload_get_device_count, offload_set_chosen_device
#endif

   IMPLICIT NONE

   TYPE(mp_para_env_type), POINTER :: para_env

   ALLOCATE (para_env)
   CALL mp_world_init(para_env)
   diag_type = FM_DIAG_TYPE_DEFAULT
#if defined(__CUSOLVERMP)
   IF (offload_get_device_count() < 1) ERROR STOP "No visible device for cuSOLVERMp"
   CALL offload_set_chosen_device(MOD(para_env%mepos, offload_get_device_count()))
#endif
   direct_generalized_diagonalization = .FALSE.
   cusolver_n_min = 64
   CALL check_solve(cusolver_n_min - 1, .FALSE.)
   CALL check_solve(cusolver_n_min, .FALSE.)
   CALL check_solve(cusolver_n_min, .TRUE.)
   direct_generalized_diagonalization = .TRUE.
   CALL check_solve(cusolver_n_min, .TRUE., cusolver_n_min - 3)
   CALL check_solve(cusolver_n_min, .TRUE.)
   CALL mp_world_finalize()
   DEALLOCATE (para_env)

CONTAINS

! **************************************************************************************************
!> \brief Check the selected backend and the returned eigenpairs on a fresh BLACS grid.
!> \param n matrix order
!> \param generalized solve the generalized eigenproblem
!> \param nvec number of requested eigenvectors, all by default
! **************************************************************************************************
   SUBROUTINE check_solve(n, generalized, nvec)
      INTEGER, INTENT(IN) :: n
      LOGICAL, INTENT(IN) :: generalized
      INTEGER, INTENT(IN), OPTIONAL :: nvec

      REAL(KIND=dp), PARAMETER :: tol = 1.0E-11_dp
      COMPLEX(KIND=dp), DIMENSION(n, n) :: h, s, sz, z
      INTEGER :: i, nev
      LOGICAL :: expect_gpu, upper_only
      REAL(KIND=dp), DIMENSION(n) :: eval
      TYPE(cp_blacs_env_type), POINTER :: grid
      TYPE(cp_cfm_type) :: a, b, vec, work
      TYPE(cp_fm_struct_type), POINTER :: layout, vector_layout

      NULLIFY (grid, layout, vector_layout)
      nev = n
      IF (PRESENT(nvec)) nev = nvec
      CALL cp_blacs_env_create(grid, para_env)
      CALL cp_fm_struct_create(layout, nrow_global=n, ncol_global=n, context=grid, &
                               nrow_block=8, ncol_block=8, force_block=.TRUE.)
      CALL cp_cfm_create(a, layout)
      IF (cp2k_is_parallel .AND. .NOT. generalized) THEN
         CALL cp_fm_struct_create(vector_layout, template_fmstruct=layout, &
                                  local_leading_dimension=layout%local_leading_dimension + 2, force_block=.TRUE.)
         CALL cp_cfm_create(vec, vector_layout)
         CALL cp_fm_struct_release(vector_layout)
      ELSE
         CALL cp_cfm_create(vec, layout)
      END IF
      h = 0.0_dp
      s = 0.0_dp
      DO i = 1, n
         h(i, i) = 2.0_dp + REAL(i, dp)/n
         s(i, i) = 1.0_dp
      END DO
      DO i = 1, n - 1
         h(i, i + 1) = (0.25_dp, 0.5_dp)
         h(i + 1, i) = CONJG(h(i, i + 1))
      END DO
      IF (generalized) THEN
         DO i = 1, n
            s(i, i) = 1.0_dp + REAL(i, dp)/n
         END DO
         DO i = 1, n - 1
            s(i, i + 1) = (0.04_dp, -0.07_dp)
            s(i + 1, i) = CONJG(s(i, i + 1))
         END DO
         IF (cp2k_is_parallel .AND. nev < n) THEN
            CALL cp_fm_struct_create(vector_layout, template_fmstruct=layout, &
                                     local_leading_dimension=layout%local_leading_dimension + 2, force_block=.TRUE.)
            CALL cp_cfm_create(b, vector_layout)
            CALL cp_cfm_create(work, vector_layout)
            CALL cp_fm_struct_release(vector_layout)
         ELSE
            CALL cp_cfm_create(b, layout)
            CALL cp_cfm_create(work, layout)
         END IF
         upper_only = direct_generalized_diagonalization .AND. &
                      (diag_type == FM_DIAG_TYPE_CUSOLVER .OR. diag_type == FM_DIAG_TYPE_SCALAPACK)
         z = h
         sz = s
         IF (upper_only) THEN
            ! The lower halves are not inputs to a native solve. Keep h and s for the residual.
            DO i = 1, n - 1
               z(i + 1:n, i) = (9.0_dp, -3.0_dp)
               sz(i + 1:n, i) = (-6.0_dp, 2.0_dp)
            END DO
         END IF
         CALL cp_cfm_set_submatrix(a, z)
         CALL cp_cfm_set_submatrix(b, sz)
         CALL cp_cfm_geeig(a, b, vec, eval(:nev), work)
         CALL cp_cfm_release(work)
         CALL cp_cfm_release(b)
      ELSE
         ! Only the upper triangle is part of the standard solver's input.
         z = h
         DO i = 1, n - 1
            z(i + 1:n, i) = (9.0_dp, -3.0_dp)
         END DO
         CALL cp_cfm_set_submatrix(a, z)
         CALL cp_cfm_heevd(a, vec, eval)
      END IF

      expect_gpu = diag_type == FM_DIAG_TYPE_CUSOLVER .AND. n >= cusolver_n_min
      IF (expect_gpu .AND. para_env%num_pe > 1) expect_gpu = grid%cusolver_devices_ok == 1
      IF (C_ASSOCIATED(grid%cusolver_context) .NEQV. expect_gpu) ERROR STOP "Wrong complex standard eigensolver"
      CALL cp_cfm_get_submatrix(vec, z(:, :nev))
      sz(:, :nev) = MATMUL(s, z(:, :nev))
      IF (ANY(eval(2:nev) < eval(:nev - 1))) ERROR STOP "Eigenvalues are not ordered"
      IF (.NOT. ALL(ABS(MATMUL(h, z(:, :nev)) - sz(:, :nev)*SPREAD(eval(:nev), 1, n)) < tol)) ERROR STOP "Large eigenpair residual"
      s(:nev, :nev) = MATMUL(CONJG(TRANSPOSE(z(:, :nev))), sz(:, :nev))
      DO i = 1, nev
         s(i, i) = s(i, i) - 1.0_dp
      END DO
      IF (.NOT. ALL(ABS(s(:nev, :nev)) < tol)) ERROR STOP "Eigenvectors are not S-orthonormal"

      CALL cp_cfm_release(vec)
      CALL cp_cfm_release(a)
      CALL cp_fm_struct_release(layout)
      CALL cp_blacs_env_release(grid)
   END SUBROUTINE check_solve
END PROGRAM cp_cfm_diag_unittest
