!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Primitive Gaussian pairs and products for four-center integrals.
!>        The fixed radial workspace capacity follows the compiled legacy backend.
! **************************************************************************************************
MODULE four_center_integral_types
   USE kinds,                           ONLY: dp
   USE libint_wrapper,                  ONLY: eri4_fundamentals_size => prim_data_f_size,&
                                              eri4_pgf_product_type => eri_primitive_data_type

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: eri4_cell_type, eri4_screen_coeff_type, eri4_pgf_image_type, &
             eri4_pgf_pair_type, eri4_pgf_product_type, eri4_fundamentals_size

! **************************************************************************************************
   TYPE eri4_cell_type
      REAL(dp)                                 :: cell(3) = 0.0_dp
      REAL(dp)                                 :: cell_r(3) = 0.0_dp
   END TYPE eri4_cell_type

! **************************************************************************************************
   TYPE eri4_screen_coeff_type
      REAL(dp)                                 :: x(2) = 0.0_dp
   END TYPE eri4_screen_coeff_type

! **************************************************************************************************
   TYPE eri4_pgf_image_type
      REAL(dp)                                 :: ra(3) = 0.0_dp, rb(3) = 0.0_dp
      REAL(dp)                                 :: rab2 = 0.0_dp
      REAL(dp)                                 :: S1234 = 0.0_dp
      REAL(dp)                                 :: P(3) = 0.0_dp
      REAL(dp)                                 :: R = 0.0_dp
      REAL(dp)                                 :: pgf_max = 0.0_dp
      REAL(dp), DIMENSION(3)                   :: bcell = 0.0_dp
   END TYPE eri4_pgf_image_type

! **************************************************************************************************
   TYPE eri4_pgf_pair_type
      TYPE(eri4_pgf_image_type), DIMENSION(:), POINTER &
         :: image_list => NULL()
      INTEGER                                  :: nimages = 0
      REAL(dp)                                 :: zetapzetb = 0.0_dp
      REAL(dp)                                 :: ZetaInv = 0.0_dp
      REAL(dp)                                 :: zeta = 0.0_dp, zetb = 0.0_dp
      INTEGER                                  :: ipgf = 0, jpgf = 0
   END TYPE eri4_pgf_pair_type

END MODULE four_center_integral_types
