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
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE wannier90_projections,           ONLY: wannier90_select_projections

   IMPLICIT NONE

   COMPLEX(KIND=dp), ALLOCATABLE :: candidates(:, :, :)
   COMPLEX(KIND=dp) :: amn(2, 2, 2), reference(2, 2, 2), rotation(2, 2, 2)
   COMPLEX(KIND=dp) :: repaired(2, 2, 3)
   INTEGER :: anchor, ik, nfail, pivots(2), reference_pivots(2)
   REAL(KIND=dp) :: angle, min_svalue
   TYPE(mp_comm_type) :: comm
   TYPE(mp_para_env_type), POINTER :: para_env

   CALL mp_world_init(comm)
   ALLOCATE (para_env)
   CALL para_env%from_dup(comm)
   nfail = 0
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
   CALL comm%sum(nfail)
   IF (para_env%is_source()) WRITE (*, '(A,I0)') 'Wannier90 projection test failures: ', nfail
   CALL mp_para_env_release(para_env)
   CALL mp_world_finalize()
   IF (nfail /= 0) STOP 1
END PROGRAM wannier90_projections_unittest
