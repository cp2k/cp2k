!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Test mesh-wide trial selection, repair and MO-gauge covariance.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
PROGRAM wannier90_projections_unittest
#if defined (__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: IEEE_ALL, ieee_get_halting_mode, ieee_set_halting_mode
#endif
   USE kinds, ONLY: dp
   USE message_passing, ONLY: mp_comm_type, mp_para_env_release, mp_para_env_type, mp_world_finalize, mp_world_init
   USE wannier90_projections, ONLY: wannier90_hybrid_candidates, wannier90_projection_rank, &
                                    wannier90_select_projections, wannier90_window_mask
#include "./base/base_uses.f90"

   IMPLICIT NONE

   COMPLEX(KIND=dp), ALLOCATABLE :: candidates(:, :, :), hybrids(:, :, :)
   COMPLEX(KIND=dp) :: amn(2, 2, 2), reference(2, 2, 2), repaired(2, 2, 3), rotation(2, 2), quartet(2, 4)
   INTEGER, PARAMETER :: groups(4, 2) = RESHAPE([1, 2, 3, 4, 4, 1, 2, 3], [4, 2])
   INTEGER :: anchor, ik, nfail, pivots(2)
   LOGICAL :: allowed(4)
   REAL(KIND=dp) :: angle, min_svalue
   TYPE(mp_comm_type) :: comm
   TYPE(mp_para_env_type), POINTER :: para_env
#if defined (__HAS_IEEE_EXCEPTIONS)
   COMPLEX(KIND=dp) :: rectangular(4, 3, 2)
   LOGICAL :: saved_halt(SIZE(IEEE_ALL)), actual_halt(SIZE(IEEE_ALL))
#endif

   CALL mp_world_init(comm)
   ALLOCATE (para_env)
   CALL para_env%from_dup(comm)
   nfail = 0
#if defined (__HAS_IEEE_EXCEPTIONS)
   CALL ieee_get_halting_mode(IEEE_ALL, saved_halt)
#endif
   IF (para_env%is_source()) THEN
      ALLOCATE (candidates(2, 3, 2))
      candidates = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      ! The strongest stacked column vanishes at the second k-point.
      candidates(1, 1, 1) = 5.0_dp
      candidates(2, 2, :) = 2.0_dp
      candidates(1, 3, :) = 0.4_dp
   ELSE
      ALLOCATE (candidates(0, 0, 0))
   END IF
   CALL wannier90_select_projections(candidates, para_env, amn, pivots, min_svalue, anchor_kpoint=0)
   IF (.NOT. (min_svalue <= 1.0E-14_dp)) nfail = nfail + 1
   CALL wannier90_select_projections(candidates, para_env, amn, pivots, min_svalue)
   IF (.NOT. (ABS(min_svalue - 0.4_dp) <= 1.0E-14_dp) .OR. ANY(pivots /= [2, 3])) nfail = nfail + 1
   reference = amn
   DO ik = 1, 2
      angle = 0.31_dp*REAL(ik, KIND=dp)
      rotation = CMPLX(0.0_dp, SIN(angle), KIND=dp)
      rotation(1, 1) = COS(angle)
      rotation(2, 2) = COS(angle)
      IF (para_env%is_source()) candidates(:, :, ik) = MATMUL(rotation, candidates(:, :, ik))
      reference(:, :, ik) = MATMUL(rotation, reference(:, :, ik))
   END DO
   CALL wannier90_select_projections(candidates, para_env, amn, pivots, min_svalue)
   IF (.NOT. (ABS(min_svalue - 0.4_dp) <= 1.0E-14_dp)) nfail = nfail + 1
   IF (ANY(pivots /= [2, 3]) .OR. .NOT. ALL(ABS(amn - reference) <= 1.0E-14_dp)) nfail = nfail + 1
   IF (para_env%is_source()) candidates = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
   CALL wannier90_select_projections(candidates, para_env, amn, pivots, min_svalue)
   IF (.NOT. (min_svalue == 0.0_dp)) nfail = nfail + 1
   IF (para_env%is_source()) THEN
      ! Numerically tied columns retain the reproducible original-index ordering.
      candidates(1, 1, :) = 1.0_dp
      candidates(1, 2, :) = 1.0_dp + 1.0E-14_dp
      candidates(2, 3, :) = 1.0_dp
   END IF
   CALL wannier90_select_projections(candidates, para_env, amn, pivots, min_svalue)
   IF (ANY(pivots /= [1, 3])) nfail = nfail + 1
   DEALLOCATE (candidates)

   IF (para_env%is_source()) THEN
      ALLOCATE (candidates(2, 4, 3))
      candidates = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      candidates(1, 1, 1) = 10.0_dp
      candidates(1, 1, 3) = 10.0_dp
      candidates(2, 2, 2) = 10.0_dp
      candidates(2, 2, 3) = 10.0_dp
      candidates(1, 3, :) = 1.0_dp
      candidates(2, 4, :) = 1.0_dp
   ELSE
      ALLOCATE (candidates(0, 0, 0))
   END IF
   ! Every fixed anchor is rank-deficient somewhere; automatic selection must repair it.
   DO anchor = 0, 3
      CALL wannier90_select_projections(candidates, para_env, repaired, pivots, min_svalue, anchor)
      IF (.NOT. (min_svalue <= 1.0E-14_dp)) nfail = nfail + 1
   END DO
   CALL wannier90_select_projections(candidates, para_env, repaired, pivots, min_svalue)
   IF (.NOT. (ABS(min_svalue - 1.0_dp) <= 1.0E-14_dp) .OR. ANY(pivots /= [3, 4])) nfail = nfail + 1
   CALL wannier90_hybrid_candidates(candidates, groups, hybrids)
   IF (para_env%is_source()) THEN
      CPASSERT(ALL(SHAPE(hybrids) == [2, 8, 3]))
      DO ik = 1, 3
         DO anchor = 1, SIZE(groups, 2)
            quartet = hybrids(:, 4*anchor - 3:4*anchor, ik)
            IF (.NOT. ALL(ABS(0.5_dp*SUM(quartet, DIM=2) - candidates(:, groups(1, anchor), ik)) &
                          <= 1.0E-14_dp)) nfail = nfail + 1
            IF (.NOT. ALL(ABS(MATMUL(quartet, CONJG(TRANSPOSE(quartet))) - &
                              MATMUL(candidates(:, :, ik), CONJG(TRANSPOSE(candidates(:, :, ik))))) &
                          <= 1.0E-14_dp)) nfail = nfail + 1
         END DO
      END DO
   ELSE
      IF (SIZE(hybrids) /= 0) nfail = nfail + 1
   END IF
   DEALLOCATE (hybrids, candidates)

   CALL wannier90_window_mask([-2.0_dp, -1.0_dp, 1.0_dp, 2.0_dp], 2, [-1.5_dp, 1.5_dp, 0.0_dp, 0.0_dp], &
                              [.TRUE., .TRUE., .FALSE., .FALSE.], allowed)
   IF (ANY(allowed .NEQV. [.FALSE., .TRUE., .TRUE., .FALSE.])) nfail = nfail + 1
   CALL wannier90_window_mask([-2.0_dp, -1.0_dp, 1.0_dp, 2.0_dp], 2, [-3.0_dp, 3.0_dp, 0.0_dp, 0.0_dp], &
                              [.TRUE., .TRUE., .FALSE., .TRUE.], allowed)
   IF (ANY(allowed .NEQV. [.TRUE., .TRUE., .FALSE., .FALSE.])) nfail = nfail + 1
#if defined (__HAS_IEEE_EXCEPTIONS)
   rectangular = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
   DO ik = 1, 3
      rectangular(ik, ik, :) = REAL(ik, KIND=dp)
      rectangular(4, ik, :) = CMPLX(0.2_dp*ik, 0.1_dp, KIND=dp)
   END DO
   min_svalue = wannier90_projection_rank(rectangular)
   CALL ieee_get_halting_mode(IEEE_ALL, actual_halt)
   IF (ANY(actual_halt .NEQV. saved_halt)) nfail = nfail + 1
   CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
   IF (.NOT. (ABS(wannier90_projection_rank(rectangular) - min_svalue) <= 1.0E-14_dp)) nfail = nfail + 1
   CALL ieee_get_halting_mode(IEEE_ALL, actual_halt)
   IF (ANY(actual_halt)) nfail = nfail + 1
   CALL ieee_set_halting_mode(IEEE_ALL, saved_halt)
#endif
   CALL comm%sum(nfail)
   CPASSERT(nfail == 0)
   CALL mp_para_env_release(para_env)
   CALL mp_world_finalize()
END PROGRAM wannier90_projections_unittest
