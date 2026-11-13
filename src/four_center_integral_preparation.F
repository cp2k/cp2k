!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Gaussian pair geometry for four-center integral screening.
! **************************************************************************************************
MODULE four_center_integral_preparation
   USE four_center_integral_types,      ONLY: eri4_pgf_image_type,&
                                              eri4_pgf_pair_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: min_image_gap_sq, prepare_eri4_pair

CONTAINS

! **************************************************************************************************
!> \brief Prepare one Gaussian pair, retaining the caller's image topology and screening metadata.
!> \param ra first center
!> \param rb second center in the selected image
!> \param zeta first Gaussian exponent
!> \param zetb second Gaussian exponent
!> \param pair exponent data; image-list ownership is unchanged
!> \param image product center for this image
! **************************************************************************************************
   SUBROUTINE prepare_eri4_pair(ra, rb, zeta, zetb, pair, image)
      REAL(dp), INTENT(IN)                               :: ra(3), rb(3), zeta, zetb
      TYPE(eri4_pgf_pair_type), INTENT(INOUT)            :: pair
      TYPE(eri4_pgf_image_type), INTENT(INOUT)           :: image

      pair%zeta = zeta
      pair%zetb = zetb
      pair%zetapzetb = zeta + zetb
      pair%ZetaInv = 1.0_dp/pair%zetapzetb
      image%ra = ra
      image%rb = rb
      image%P = (zeta*ra + zetb*rb)*pair%ZetaInv
   END SUBROUTINE prepare_eri4_pair

! **************************************************************************************************
!> \brief Squared-distance gap to the next image for an orthorhombic cell.
!> \param r minimum-image displacement
!> \param hmat cell matrix
!> \return smallest squared-distance increment to an adjacent image
! **************************************************************************************************
   PURE FUNCTION min_image_gap_sq(r, hmat) RESULT(gap)
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: r
      REAL(dp), DIMENSION(3, 3), INTENT(IN)              :: hmat
      REAL(dp)                                           :: gap

      REAL(dp), DIMENSION(3)                             :: cell_len

      cell_len = [ABS(hmat(1, 1)), ABS(hmat(2, 2)), ABS(hmat(3, 3))]
      gap = MINVAL(cell_len*(cell_len - 2.0_dp*ABS(r)))

   END FUNCTION min_image_gap_sq

END MODULE four_center_integral_preparation
