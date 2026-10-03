!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Test mesh-wide rank validation and unitary covariance of AO trial selection.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
PROGRAM wannier90_projections_unittest
#if defined (__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: IEEE_ALL, IEEE_DIVIDE_BY_ZERO, IEEE_INVALID, IEEE_OVERFLOW, &
                              ieee_get_halting_mode, ieee_set_halting_mode, ieee_support_halting
#endif
   USE kinds, ONLY: dp
   USE message_passing, ONLY: mp_comm_type, &
                              mp_para_env_release, &
                              mp_para_env_type, &
                              mp_world_finalize, &
                              mp_world_init
   USE wannier90_projections, ONLY: wannier90_hybrid_candidates, &
                                    wannier90_mix_projections, &
                                    wannier90_projection_rank, &
                                    wannier90_select_projections
#include "./base/base_uses.f90"

   IMPLICIT NONE

   COMPLEX(KIND=dp), ALLOCATABLE :: candidates(:, :, :)
   COMPLEX(KIND=dp), ALLOCATABLE :: hybrids(:, :, :), rotated_hybrids(:, :, :)
   COMPLEX(KIND=dp) :: quartet(2, 4)
   COMPLEX(KIND=dp) :: amn(2, 2, 2), reference(2, 2, 2), rotation(2, 2, 2)
   COMPLEX(KIND=dp) :: repaired(2, 2, 3)
   COMPLEX(KIND=dp) :: mixed(4, 3, 2), mixed_reference(4, 3, 2), single(4, 1, 2), gauge(4, 4, 2)
   INTEGER :: anchor, ik, nfail, pivots(2), reference_pivots(2)
   REAL(KIND=dp) :: angle, min_svalue
   TYPE(mp_comm_type) :: comm
   TYPE(mp_para_env_type), POINTER :: para_env
#if defined (__HAS_IEEE_EXCEPTIONS)
   LOGICAL :: saved_halt(SIZE(IEEE_ALL)), expected_halt(SIZE(IEEE_ALL)), actual_halt(SIZE(IEEE_ALL))
#endif

   CALL mp_world_init(comm)
   ALLOCATE (para_env)
   CALL para_env%from_dup(comm)
   nfail = 0
#if defined (__HAS_IEEE_EXCEPTIONS)
   CALL ieee_get_halting_mode(IEEE_ALL, saved_halt)
   IF (ieee_support_halting(IEEE_INVALID)) CALL ieee_set_halting_mode(IEEE_INVALID, .TRUE.)
   IF (ieee_support_halting(IEEE_DIVIDE_BY_ZERO)) CALL ieee_set_halting_mode(IEEE_DIVIDE_BY_ZERO, .TRUE.)
   IF (ieee_support_halting(IEEE_OVERFLOW)) CALL ieee_set_halting_mode(IEEE_OVERFLOW, .TRUE.)
   CALL ieee_get_halting_mode(IEEE_ALL, expected_halt)
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
   IF (min_svalue > 1.0e-14_dp) nfail = nfail + 1
   CALL wannier90_select_projections(candidates, para_env, amn, pivots, min_svalue)
   IF (ABS(min_svalue - 0.4_dp) > 1.0e-14_dp) nfail = nfail + 1
   IF (.NOT. ALL(pivots == [2, 3])) nfail = nfail + 1
   reference = amn
   reference_pivots = pivots

   DO ik = 1, 2
      angle = 0.31_dp*REAL(ik, KIND=dp)
      rotation(1, 1, ik) = COS(angle)
      rotation(2, 2, ik) = COS(angle)
      rotation(1, 2, ik) = CMPLX(0.0_dp, SIN(angle), KIND=dp)
      rotation(2, 1, ik) = CMPLX(0.0_dp, SIN(angle), KIND=dp)
      IF (para_env%is_source()) candidates(:, :, ik) = MATMUL(rotation(:, :, ik), candidates(:, :, ik))
   END DO
   CALL wannier90_select_projections(candidates, para_env, amn, pivots, min_svalue)
   IF (ABS(min_svalue - 0.4_dp) > 1.0e-14_dp) nfail = nfail + 1
   IF (.NOT. ALL(pivots == reference_pivots)) nfail = nfail + 1
   DO ik = 1, 2
      IF (MAXVAL(ABS(amn(:, :, ik) - MATMUL(rotation(:, :, ik), reference(:, :, ik)))) > 1.0e-14_dp) THEN
         nfail = nfail + 1
      END IF
   END DO
   IF (para_env%is_source()) candidates(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
   CALL wannier90_select_projections(candidates, para_env, amn, pivots, min_svalue)
   IF (min_svalue /= 0.0_dp) nfail = nfail + 1
   DEALLOCATE (candidates)
   IF (para_env%is_source()) THEN
      ALLOCATE (candidates(2, 5, 2))
      DO ik = 1, 2
         DO anchor = 1, 5
            candidates(1, anchor, ik) = CMPLX(0.1_dp*anchor, 0.2_dp*ik, KIND=dp)
            candidates(2, anchor, ik) = CMPLX(-0.3_dp*anchor, 0.05_dp*ik*anchor, KIND=dp)
         END DO
      END DO
   ELSE
      ALLOCATE (candidates(0, 0, 0))
   END IF
   CALL wannier90_hybrid_candidates(candidates, RESHAPE([1, 2, 3, 4], [4, 1]), hybrids)
   IF (para_env%is_source()) THEN
      IF (.NOT. ALL(SHAPE(hybrids) == [2, 5, 2])) nfail = nfail + 1
      IF (MAXVAL(ABS(hybrids(:, 1, :) - candidates(:, 5, :))) > 1.0e-14_dp) nfail = nfail + 1
      DO ik = 1, 2
         ! A full tetrahedral quartet is a unitary change of the four AO trials.
         quartet = hybrids(:, 2:5, ik)
         IF (MAXVAL(ABS(MATMUL(quartet, CONJG(TRANSPOSE(quartet))) - &
                        MATMUL(candidates(:, 1:4, ik), CONJG(TRANSPOSE(candidates(:, 1:4, ik)))))) &
             > 1.0e-14_dp) nfail = nfail + 1
         candidates(:, :, ik) = MATMUL(rotation(:, :, ik), candidates(:, :, ik))
      END DO
   END IF
   CALL wannier90_hybrid_candidates(candidates, RESHAPE([1, 2, 3, 4], [4, 1]), rotated_hybrids)
   IF (para_env%is_source()) THEN
      DO ik = 1, 2
         IF (MAXVAL(ABS(rotated_hybrids(:, :, ik) - MATMUL(rotation(:, :, ik), hybrids(:, :, ik)))) &
             > 1.0e-14_dp) nfail = nfail + 1
      END DO
   ELSE
      IF (SIZE(hybrids) /= 0 .OR. SIZE(rotated_hybrids) /= 0) nfail = nfail + 1
   END IF
   DEALLOCATE (hybrids, rotated_hybrids)
   CALL wannier90_hybrid_candidates(candidates, RESHAPE([1, 2, 3, 4, 5, 2, 3, 4], [4, 2]), hybrids)
   IF (para_env%is_source()) THEN
      IF (.NOT. ALL(SHAPE(hybrids) == [2, 8, 2])) nfail = nfail + 1
      DO ik = 1, 2
         IF (MAXVAL(ABS(SUM(hybrids(:, 1:4, ik), DIM=2)*0.5_dp - candidates(:, 1, ik))) &
             > 1.0e-14_dp) nfail = nfail + 1
         IF (MAXVAL(ABS(SUM(hybrids(:, 5:8, ik), DIM=2)*0.5_dp - candidates(:, 5, ik))) &
             > 1.0e-14_dp) nfail = nfail + 1
      END DO
   END IF
   DEALLOCATE (hybrids)
   CALL wannier90_hybrid_candidates(candidates, RESHAPE([1, 2, 3, 4], [4, 0]), hybrids)
   IF (para_env%is_source()) THEN
      IF (MAXVAL(ABS(hybrids - candidates)) > 1.0e-14_dp) nfail = nfail + 1
   END IF
   DEALLOCATE (hybrids, candidates)
   IF (para_env%is_source()) THEN
      ALLOCATE (candidates(2, 4, 3))
      candidates = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      ! Every QR anchor chooses at least one column that vanishes elsewhere.
      ! Exchanging both of those columns recovers the fixed full-rank pair.
      candidates(1, 1, 1) = 10.0_dp
      candidates(1, 1, 3) = 10.0_dp
      candidates(2, 2, 2) = 10.0_dp
      candidates(2, 2, 3) = 10.0_dp
      candidates(1, 3, :) = 1.0_dp
      candidates(2, 4, :) = 1.0_dp
   ELSE
      ALLOCATE (candidates(0, 0, 0))
   END IF
   DO anchor = 0, 3
      CALL wannier90_select_projections(candidates, para_env, repaired, pivots, min_svalue, anchor)
      IF (min_svalue > 1.0e-14_dp) nfail = nfail + 1
   END DO
   CALL wannier90_select_projections(candidates, para_env, repaired, pivots, min_svalue)
   IF (ABS(min_svalue - 1.0_dp) > 1.0e-14_dp) nfail = nfail + 1
   IF (.NOT. ALL(pivots == [3, 4])) nfail = nfail + 1
   DEALLOCATE (candidates)
   IF (para_env%is_source()) THEN
      ALLOCATE (candidates(2, 3, 2))
      candidates = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      candidates(1, 1, :) = 1.0_dp
      candidates(1, 2, :) = 1.0_dp + 1.0e-14_dp
      candidates(2, 3, :) = 1.0_dp
   ELSE
      ALLOCATE (candidates(0, 0, 0))
   END IF
   CALL wannier90_select_projections(candidates, para_env, amn, pivots, min_svalue)
   IF (.NOT. ALL(pivots == [1, 3])) nfail = nfail + 1
   IF (para_env%is_source()) THEN
      DO ik = 1, 2
         candidates(:, :, ik) = MATMUL(rotation(:, :, ik), candidates(:, :, ik))
      END DO
   END IF
   CALL wannier90_select_projections(candidates, para_env, amn, pivots, min_svalue)
   IF (.NOT. ALL(pivots == [1, 3])) nfail = nfail + 1
   DEALLOCATE (candidates)
   gauge = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
   gauge(1:2, 1:2, :) = rotation
   gauge(3, 3, :) = 1.0_dp
   gauge(4, 4, :) = 1.0_dp
   DO ik = 1, 2
      DO anchor = 1, 3
         mixed_reference(1, anchor, ik) = CMPLX(0.1_dp*anchor, 0.2_dp*ik, KIND=dp)
         mixed_reference(2, anchor, ik) = CMPLX(-0.3_dp*anchor, 0.05_dp*ik*anchor, KIND=dp)
         mixed_reference(3, anchor, ik) = CMPLX(0.2_dp*anchor**2, 0.3_dp*ik, KIND=dp)
         mixed_reference(4, anchor, ik) = CMPLX(0.7_dp, -0.1_dp*ik*anchor, KIND=dp)
      END DO
   END DO
   min_svalue = wannier90_projection_rank(mixed_reference)
#if defined (__HAS_IEEE_EXCEPTIONS)
   CALL ieee_get_halting_mode(IEEE_ALL, actual_halt)
   IF (ANY(actual_halt .NEQV. expected_halt)) nfail = nfail + 1
   ! Also check that a caller with traps disabled is left unchanged.
   CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
   IF (ABS(wannier90_projection_rank(mixed_reference) - min_svalue) > 1.0e-14_dp) nfail = nfail + 1
   CALL ieee_get_halting_mode(IEEE_ALL, actual_halt)
   IF (ANY(actual_halt)) nfail = nfail + 1
   CALL ieee_set_halting_mode(IEEE_ALL, expected_halt)
#endif
   mixed = mixed_reference
   CALL wannier90_mix_projections(mixed, 0.53_dp)
   IF (ABS(wannier90_projection_rank(mixed) - min_svalue) > 1.0e-14_dp) nfail = nfail + 1
   DO ik = 1, 2
      IF (MAXVAL(ABS(MATMUL(mixed(:, :, ik), CONJG(TRANSPOSE(mixed(:, :, ik)))) - &
                     MATMUL(mixed_reference(:, :, ik), CONJG(TRANSPOSE(mixed_reference(:, :, ik)))))) &
          > 1.0e-14_dp) nfail = nfail + 1
      mixed_reference(:, :, ik) = MATMUL(gauge(:, :, ik), mixed_reference(:, :, ik))
   END DO
   CALL wannier90_mix_projections(mixed_reference, 0.53_dp)
   DO ik = 1, 2
      IF (MAXVAL(ABS(mixed_reference(:, :, ik) - MATMUL(gauge(:, :, ik), mixed(:, :, ik)))) &
          > 1.0e-14_dp) nfail = nfail + 1
   END DO
   single = mixed(:, 1:1, :)
   CALL wannier90_mix_projections(single, -0.53_dp)
   IF (MAXVAL(ABS(single - mixed(:, 1:1, :))) /= 0.0_dp) nfail = nfail + 1
#if defined (__HAS_IEEE_EXCEPTIONS)
   CALL ieee_get_halting_mode(IEEE_ALL, actual_halt)
   IF (ANY(actual_halt .NEQV. expected_halt)) nfail = nfail + 1
   CALL ieee_set_halting_mode(IEEE_ALL, saved_halt)
#endif
   CALL comm%sum(nfail)
   CPASSERT(nfail == 0)
   CALL mp_para_env_release(para_env)
   CALL mp_world_finalize()
END PROGRAM wannier90_projections_unittest
