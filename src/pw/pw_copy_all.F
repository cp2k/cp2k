!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
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
   USE message_passing,                 ONLY: mp_shift,&
                                              mp_sum
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              RECIPROCALSPACE,&
                                              pw_type
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
      TYPE(pw_type), INTENT(IN)                          :: pw1
      TYPE(pw_type), INTENT(INOUT)                       :: pw2

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: cc
      INTEGER                                            :: group, group_size, ig1, ig2, ip, jg2, &
                                                            me, ng1, ng2, ngm, penow
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ngr
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: g_hat
      INTEGER, DIMENSION(3)                              :: k1, k2
      TYPE(pw_grid_type), POINTER                        :: pg1, pg2

      IF (pw1%in_use == COMPLEXDATA1D .AND. &
          pw2%in_use == COMPLEXDATA1D .AND. &
          pw1%in_space == RECIPROCALSPACE) THEN

         ng1 = SIZE(pw1%cc)
         ng2 = SIZE(pw2%cc)

         pg1 => pw1%pw_grid
         pg2 => pw2%pw_grid

         group = pg1%para%group
         group_size = pg1%para%group_size
         me = pg1%para%my_pos
         ALLOCATE (ngr(group_size))
         ngr = 0
         ngr(me + 1) = pg1%ngpts_cut_local
         CALL mp_sum(ngr, group)
         ngm = MAXVAL(ngr)
         ALLOCATE (cc(ngm))
         cc(1:ng1) = pw1%cc(1:ng1)
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
               IF (jg2 /= 0) pw2%cc(jg2) = cc(ig1)
            END DO
            IF (ip /= group_size) THEN
               CALL mp_shift(cc, group)
               CALL mp_shift(g_hat, group)
            END IF

         END DO

         DEALLOCATE (ngr, cc, g_hat)

      ELSE
         CPABORT("No suitable data field")
      END IF

      pw2%in_space = RECIPROCALSPACE

   END SUBROUTINE pw_copy_match

END MODULE pw_copy_all

