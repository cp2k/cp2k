!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!> \author Jan Wilhelm
!> \date 05.2024
! **************************************************************************************************
MODULE gw_kp_to_real_space_and_back
   USE cp_cfm_types,                    ONLY: cp_cfm_type
   USE cp_fm_types,                     ONLY: cp_fm_set_all,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE mathconstants,                   ONLY: gaussi,&
                                              twopi,&
                                              z_one,&
                                              z_zero
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_kp_to_real_space_and_back'

   PUBLIC :: fm_trafo_rs_to_ikp, trafo_rs_to_ikp, trafo_ikp_to_rs, fm_add_ikp_to_rs, &
             add_ikp_to_all_rs

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param cfm_ikp ...
!> \param fm_rs ...
!> \param kpoints ...
!> \param ikp ...
! **************************************************************************************************
   SUBROUTINE fm_trafo_rs_to_ikp(cfm_ikp, fm_rs, kpoints, ikp)
      TYPE(cp_cfm_type)                                  :: cfm_ikp
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_rs
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ikp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_trafo_rs_to_ikp'

      INTEGER                                            :: handle, img, nimages, nimages_fm_rs

      CALL timeset(routineN, handle)

      nimages = SIZE(kpoints%index_to_cell, 1)
      nimages_fm_rs = SIZE(fm_rs)

      CPASSERT(nimages == nimages_fm_rs)

      cfm_ikp%local_data(:, :) = z_zero
      DO img = 1, nimages

         CALL add_rs_to_ikp(fm_rs(img)%local_data, cfm_ikp%local_data, kpoints%index_to_cell, &
                            kpoints%xkp(1:3, ikp), img)

      END DO

      CALL timestop(handle)

   END SUBROUTINE fm_trafo_rs_to_ikp

! **************************************************************************************************
!> \brief ...
!> \param array_rs ...
!> \param array_kp ...
!> \param index_to_cell ...
!> \param xkp ...
! **************************************************************************************************
   SUBROUTINE trafo_rs_to_ikp(array_rs, array_kp, index_to_cell, xkp)
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: array_rs
      COMPLEX(KIND=dp), DIMENSION(:, :)                  :: array_kp
      INTEGER, DIMENSION(:, :)                           :: index_to_cell
      REAL(KIND=dp)                                      :: xkp(3)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'trafo_rs_to_ikp'

      INTEGER                                            :: handle, i_cell, nimages

      CALL timeset(routineN, handle)

      nimages = SIZE(index_to_cell, 1)

      CPASSERT(nimages == SIZE(array_rs, 3))

      array_kp(:, :) = 0.0_dp
      DO i_cell = 1, nimages

         CALL add_rs_to_ikp(array_rs(:, :, i_cell), array_kp, index_to_cell, xkp, i_cell)

      END DO

      CALL timestop(handle)

   END SUBROUTINE trafo_rs_to_ikp

! **************************************************************************************************
!> \brief ...
!> \param array_rs ...
!> \param array_kp ...
!> \param index_to_cell ...
!> \param xkp ...
!> \param i_cell ...
! **************************************************************************************************
   SUBROUTINE add_rs_to_ikp(array_rs, array_kp, index_to_cell, xkp, i_cell)
      REAL(KIND=dp), DIMENSION(:, :)                     :: array_rs
      COMPLEX(KIND=dp), DIMENSION(:, :)                  :: array_kp
      INTEGER, DIMENSION(:, :)                           :: index_to_cell
      REAL(KIND=dp)                                      :: xkp(3)
      INTEGER                                            :: i_cell

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'add_rs_to_ikp'

      COMPLEX(KIND=dp)                                   :: expikR
      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: arg

      CALL timeset(routineN, handle)

      arg = REAL(index_to_cell(i_cell, 1), dp)*xkp(1) + &
            REAL(index_to_cell(i_cell, 2), dp)*xkp(2) + &
            REAL(index_to_cell(i_cell, 3), dp)*xkp(3)

      expikR = z_one*COS(twopi*arg) + gaussi*SIN(twopi*arg)

      array_kp(:, :) = array_kp(:, :) + expikR*array_rs(:, :)

      CALL timestop(handle)

   END SUBROUTINE add_rs_to_ikp

! **************************************************************************************************
!> \brief ...
!> \param array_kp ...
!> \param array_rs ...
!> \param cell ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE trafo_ikp_to_rs(array_kp, array_rs, cell, kpoints)
      COMPLEX(KIND=dp), DIMENSION(:, :, :)               :: array_kp
      REAL(KIND=dp), DIMENSION(:, :)                     :: array_rs
      INTEGER                                            :: cell(3)
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'trafo_ikp_to_rs'

      INTEGER                                            :: handle, ikp

      CALL timeset(routineN, handle)

      CPASSERT(kpoints%nkp == SIZE(array_kp, 3))

      array_rs(:, :) = 0.0_dp

      DO ikp = 1, kpoints%nkp

         CALL add_ikp_to_rs(array_kp(:, :, ikp), array_rs, cell, kpoints, ikp)

      END DO

      CALL timestop(handle)

   END SUBROUTINE trafo_ikp_to_rs

! **************************************************************************************************
!> \brief ...
!> \param cfm_ikp ...
!> \param fm_rs ...
!> \param kpoints ...
!> \param ikp ...
! **************************************************************************************************
   SUBROUTINE fm_add_ikp_to_rs(cfm_ikp, fm_rs, kpoints, ikp)
      TYPE(cp_cfm_type)                                  :: cfm_ikp
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_rs
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ikp

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fm_add_ikp_to_rs'

      INTEGER                                            :: handle, img, nimages, nimages_fm_rs
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell

      CALL timeset(routineN, handle)

      nimages = SIZE(kpoints%index_to_cell, 1)
      nimages_fm_rs = SIZE(fm_rs)

      CPASSERT(nimages == nimages_fm_rs)

      ALLOCATE (index_to_cell(nimages, 3))
      index_to_cell(1:nimages, 1:3) = kpoints%index_to_cell(1:nimages, 1:3)

      DO img = 1, nimages

         IF (ikp == 1) CALL cp_fm_set_all(fm_rs(img), 0.0_dp)

         CALL add_ikp_to_rs(cfm_ikp%local_data(:, :), fm_rs(img)%local_data, &
                            index_to_cell(img, 1:3), kpoints, ikp)

      END DO

      CALL timestop(handle)

   END SUBROUTINE fm_add_ikp_to_rs

! **************************************************************************************************
!> \brief ...
!> \param array_kp ...
!> \param array_rs ...
!> \param kpoints ...
!> \param ikp ...
!> \param index_to_cell_ext ...
! **************************************************************************************************
   SUBROUTINE add_ikp_to_all_rs(array_kp, array_rs, kpoints, ikp, index_to_cell_ext)
      COMPLEX(KIND=dp), DIMENSION(:, :)                  :: array_kp
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: array_rs
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ikp
      INTEGER, DIMENSION(:, :), OPTIONAL, POINTER        :: index_to_cell_ext

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'add_ikp_to_all_rs'

      INTEGER                                            :: cell(3), handle, img, nimages
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell

      CALL timeset(routineN, handle)

      IF (PRESENT(index_to_cell_ext)) THEN
         index_to_cell => index_to_cell_ext
      ELSE
         index_to_cell => kpoints%index_to_cell
      END IF

      nimages = SIZE(index_to_cell, 1)
      CPASSERT(SIZE(array_rs, 3) == nimages)
      DO img = 1, nimages

         cell(1:3) = index_to_cell(img, 1:3)

         CALL add_ikp_to_rs(array_kp, array_rs(:, :, img), cell, kpoints, ikp)

      END DO

      CALL timestop(handle)

   END SUBROUTINE add_ikp_to_all_rs

! **************************************************************************************************
!> \brief ...
!> \param array_kp ...
!> \param array_rs ...
!> \param cell ...
!> \param kpoints ...
!> \param ikp ...
! **************************************************************************************************
   SUBROUTINE add_ikp_to_rs(array_kp, array_rs, cell, kpoints, ikp)
      COMPLEX(KIND=dp), DIMENSION(:, :)                  :: array_kp
      REAL(KIND=dp), DIMENSION(:, :)                     :: array_rs
      INTEGER                                            :: cell(3)
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ikp

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'add_ikp_to_rs'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: arg, im, re

      CALL timeset(routineN, handle)

      arg = REAL(cell(1), dp)*kpoints%xkp(1, ikp) + &
            REAL(cell(2), dp)*kpoints%xkp(2, ikp) + &
            REAL(cell(3), dp)*kpoints%xkp(3, ikp)

      re = COS(twopi*arg)*kpoints%wkp(ikp)
      im = SIN(twopi*arg)*kpoints%wkp(ikp)

      array_rs(:, :) = array_rs(:, :) + re*REAL(array_kp(:, :)) + im*AIMAG(array_kp(:, :))

      CALL timestop(handle)

   END SUBROUTINE add_ikp_to_rs

END MODULE gw_kp_to_real_space_and_back
