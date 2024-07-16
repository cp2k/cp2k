!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Definition of the type to handle the 1/R^3 residual integral part
!> \author Teodoro Laino [tlaino] - 12.2008
! **************************************************************************************************
MODULE semi_empirical_expns3_types

   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

! **************************************************************************************************
!> \brief 1/R^3 expansion type
!> \author Teodoro Laino [tlaino] - 12.2008
! **************************************************************************************************
   TYPE semi_empirical_expns3_type
      REAL(KIND=dp)                                  :: core_core = -1.0_dp
      REAL(KIND=dp), DIMENSION(9)                    :: e1b = -1.0_dp, e2a = -1.0_dp
      REAL(KIND=dp), DIMENSION(81)                   :: w = -1.0_dp
   END TYPE semi_empirical_expns3_type

   ! *** Global parameters ***
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_expns3_types'

   PUBLIC :: semi_empirical_expns3_create, &
             semi_empirical_expns3_release, &
             semi_empirical_expns3_type

CONTAINS

! **************************************************************************************************
!> \brief Allocate semi-empirical 1/R^3 expansion type
!> \param expns3 ...
!> \author Teodoro Laino [tlaino] - 12.2008
! **************************************************************************************************
   SUBROUTINE semi_empirical_expns3_create(expns3)
      TYPE(semi_empirical_expns3_type), INTENT(OUT)      :: expns3

      expns3%core_core = 0.0_dp
      expns3%e1b = 0.0_dp
      expns3%e2a = 0.0_dp
      expns3%w = 0.0_dp
   END SUBROUTINE semi_empirical_expns3_create

! **************************************************************************************************
!> \brief Deallocate the semi-empirical type
!> \param expns3 ...
!> \author Teodoro Laino [tlaino] - 12.2008
! **************************************************************************************************
   SUBROUTINE semi_empirical_expns3_release(expns3)
      TYPE(semi_empirical_expns3_type), INTENT(INOUT)    :: expns3

      MARK_USED(expns3)

   END SUBROUTINE semi_empirical_expns3_release

END MODULE semi_empirical_expns3_types
