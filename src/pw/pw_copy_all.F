!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \note
!>      Copy of pw types using an optimal match strategy
!> \par History
!>      JGH (06-May-2021) : pw_copy routine for complex match
!> \author JGH
! **************************************************************************************************
MODULE pw_copy_all
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_types,                        ONLY: pw_c1d_gs_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pw_copy_match

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_copy_all'
   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief copy a pw type variable
!> \param pw1 ...
!> \param pw2 ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE pw_copy_match(pw1, pw2)
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: pw1
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: pw2

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: cc
      INTEGER                                            :: group_size, ig1, ig2, ip, jg2, me, ng1, &
                                                            ng2, ngm, penow
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ngr
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: g_hat
      INTEGER, DIMENSION(3)                              :: k1, k2
      TYPE(mp_comm_type)                                 :: group
      TYPE(pw_grid_type), POINTER                        :: pg1, pg2

      ng1 = SIZE(pw1%array)
      ng2 = SIZE(pw2%array)

      pg1 => pw1%pw_grid
      pg2 => pw2%pw_grid

      group = pg1%para%group
      group_size = pg1%para%group%num_pe
      me = pg1%para%group%mepos
      ALLOCATE (ngr(group_size))
      ngr = 0
      ngr(me + 1) = pg1%ngpts_cut_local
      CALL group%sum(ngr)
      ngm = MAXVAL(ngr)
      ALLOCATE (cc(ngm))
      cc(1:ng1) = pw1%array(1:ng1)
      cc(ng1 + 1:ngm) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

      ALLOCATE (g_hat(3, ngm))
      g_hat = 0
      g_hat(1:3, 1:ng1) = pg1%g_hat(1:3, 1:ng1)

      DO ip = 1, group_size
         penow = me - ip + 1
         IF (penow < 0) penow = penow + group_size

         DO ig1 = 1, ngr(penow + 1)
            k1(1:3) = g_hat(1:3, ig1)
            jg2 = 0
            DO ig2 = 1, ng2
               k2(1:3) = pg2%g_hat(1:3, ig2)
               IF (SUM(ABS(k1 - k2)) == 0) THEN
                  jg2 = ig2
                  EXIT
               END IF
            END DO
            IF (jg2 /= 0) pw2%array(jg2) = cc(ig1)
         END DO
         IF (ip /= group_size) THEN
            CALL group%shift(cc)
            CALL group%shift(g_hat)
         END IF

      END DO

      DEALLOCATE (ngr, cc, g_hat)

   END SUBROUTINE pw_copy_match

END MODULE pw_copy_all

