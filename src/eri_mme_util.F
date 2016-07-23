!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Some utility methods used in different contexts.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE eri_mme_util

   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: G_abs_min, R_abs_min
CONTAINS
! **************************************************************************************************
!> \brief Find minimum length of R vectors, for a general (not necessarily
!>        orthorhombic) cell.
!> \param hmat ...
!> \retval R_m ...
! **************************************************************************************************
   FUNCTION R_abs_min(hmat) RESULT(R_m)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat
      REAL(KIND=dp)                                      :: R_m

      INTEGER                                            :: sx, sy, sz
      INTEGER, DIMENSION(3)                              :: sxyz
      REAL(KIND=dp)                                      :: R_sq
      REAL(KIND=dp), DIMENSION(3)                        :: R

      R_m = 0.0_dp

      DO sx = -1, 1
      DO sy = -1, 1
      DO sz = -1, 1
         IF (.NOT. (sx == 0 .AND. sy == 0 .AND. sz == 0)) THEN
            sxyz = [sx, sy, sz]
            R = MATMUL(hmat, sxyz)
            R_sq = R(1)**2+R(2)**2+R(3)**2
            IF (R_sq < R_m .OR. R_m < EPSILON(R_m)) R_m = R_sq
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      R_m = SQRT(R_m)

   END FUNCTION R_abs_min

! **************************************************************************************************
!> \brief Find minimum length of G vectors, for a general (not necessarily
!>        orthorhombic) cell.
!> \param h_inv ...
!> \retval G_m ...
! **************************************************************************************************
   FUNCTION G_abs_min(h_inv) RESULT(G_m)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: h_inv
      REAL(KIND=dp)                                      :: G_m

      INTEGER                                            :: gx, gy, gz
      INTEGER, DIMENSION(3)                              :: gxyz
      REAL(KIND=dp)                                      :: G_sq
      REAL(KIND=dp), DIMENSION(3)                        :: G
      REAL(KIND=dp), DIMENSION(3, 3)                     :: H

      H = twopi*TRANSPOSE(h_inv)
      G_m = 0.0_dp

      DO gx = -1, 1
      DO gy = -1, 1
      DO gz = -1, 1
         IF (.NOT. (gx == 0 .AND. gy == 0 .AND. gz == 0)) THEN
            gxyz = [gx, gy, gz]
            G = MATMUL(H, gxyz)
            G_sq = G(1)**2+G(2)**2+G(3)**2
            IF (G_sq < G_m .OR. G_m < EPSILON(G_m)) G_m = G_sq
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      G_m = SQRT(G_m)

   END FUNCTION G_abs_min

END MODULE eri_mme_util
