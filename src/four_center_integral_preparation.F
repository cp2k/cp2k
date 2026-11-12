!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Pair geometry for screening and external radial seeds for the TShPSC recurrence.
! **************************************************************************************************
MODULE four_center_integral_preparation
   USE four_center_integral_types,      ONLY: eri4_pgf_image_type,&
                                              eri4_pgf_pair_type
   USE input_constants,                 ONLY: do_potential_TShPSC
   USE integral_library_types,          ONLY: coulomb_operator_type
   USE kinds,                           ONLY: dp
   USE libint_wrapper,                  ONLY: tshpsc_fundamentals_size,&
                                              tshpsc_primitive_data_type
   USE mathconstants,                   ONLY: pi
   USE r12_fundamentals,                ONLY: compute_tshpsc_fundamentals
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: min_image_gap_sq, prepare_eri4_pair, prepare_eri4_product, prepare_eri4_primitive

CONTAINS

! **************************************************************************************************
!> \brief Prepare one Gaussian pair, retaining the caller's image topology and screening metadata.
!> \param ra first center
!> \param rb second center in the selected image
!> \param zeta first Gaussian exponent
!> \param zetb second Gaussian exponent
!> \param rab2 squared center distance, reusable across primitive pairs
!> \param pair exponent data; image-list ownership is unchanged
!> \param image product center and logarithmic Gaussian attenuation for this image
! **************************************************************************************************
   SUBROUTINE prepare_eri4_pair(ra, rb, zeta, zetb, rab2, pair, image)
      REAL(dp), INTENT(IN)                               :: ra(3), rb(3), zeta, zetb, rab2
      TYPE(eri4_pgf_pair_type), INTENT(INOUT)            :: pair
      TYPE(eri4_pgf_image_type), INTENT(INOUT)           :: image

      pair%zeta = zeta
      pair%zetb = zetb
      pair%zetapzetb = zeta + zetb
      pair%ZetaInv = 1.0_dp/pair%zetapzetb
      image%ra = ra
      image%rb = rb
      image%S1234 = -zeta*zetb*pair%ZetaInv*rab2
      image%P = (zeta*ra + zetb*rb)*pair%ZetaInv
   END SUBROUTINE prepare_eri4_pair

! **************************************************************************************************
!> \brief Prepare the image-independent quartet data once per combination of Gaussian pairs.
!> \param pair_a bra Gaussian pair
!> \param pair_b ket Gaussian pair
!> \param m_max maximum radial order, including any derivative increment
!> \param product reusable template; centers and radial seeds are filled for each selected image
! **************************************************************************************************
   SUBROUTINE prepare_eri4_product(pair_a, pair_b, m_max, product)
      TYPE(eri4_pgf_pair_type), INTENT(IN)               :: pair_a, pair_b
      INTEGER, INTENT(IN)                                :: m_max
      TYPE(tshpsc_primitive_data_type), INTENT(INOUT)    :: product

      CPASSERT(m_max < tshpsc_fundamentals_size)

      product%prepared_max_order = m_max
      product%exponents = [pair_a%zeta, pair_a%zetb, pair_b%zeta, pair_b%zetb]
      product%ZetaInv = pair_a%ZetaInv
      product%EtaInv = pair_b%ZetaInv
      product%ZetapEtaInv = pair_a%zetapzetb + pair_b%zetapzetb
      product%ZetapEtaInv = 1.0_dp/product%ZetapEtaInv
      product%Rho = pair_a%zetapzetb*pair_b%zetapzetb*product%ZetapEtaInv
   END SUBROUTINE prepare_eri4_product

! **************************************************************************************************
!> \brief Complete a prepared quartet for one selected image, including normalized radial seeds.
!> \param product image-independent quartet template
!> \param image_a bra pair geometry
!> \param image_b ket pair geometry after the caller's image translation
!> \param attenuation exponential of the two pair logarithmic attenuations, reusable across translations
!> \param pq2 squared product-center distance, already available from cutoff screening
!> \param potential_parameter radial operator parameters
!> \param params complete primitive in original center order, ready for backend permutation
! **************************************************************************************************
   SUBROUTINE prepare_eri4_primitive(product, image_a, image_b, attenuation, pq2, potential_parameter, params)
      TYPE(tshpsc_primitive_data_type), INTENT(IN)       :: product
      TYPE(eri4_pgf_image_type), INTENT(IN)              :: image_a, image_b
      REAL(dp), INTENT(IN)                               :: attenuation, pq2
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter
      TYPE(tshpsc_primitive_data_type), INTENT(INOUT)    :: params

      INTEGER                                            :: i, m_max
      REAL(dp)                                           :: normalization, prefactor, T

      CPASSERT(potential_parameter%potential_type == do_potential_TShPSC)
      m_max = product%prepared_max_order
      params%prepared_max_order = m_max
      params%exponents = product%exponents
      params%ZetaInv = product%ZetaInv
      params%EtaInv = product%EtaInv
      params%ZetapEtaInv = product%ZetapEtaInv
      params%Rho = product%Rho
      params%centers(:, 1) = image_a%ra
      params%centers(:, 2) = image_a%rb
      params%centers(:, 3) = image_b%ra
      params%centers(:, 4) = image_b%rb
      params%P = image_a%P
      params%Q = image_b%P

      params%W = ((product%exponents(1) + product%exponents(2))*params%P + &
                  (product%exponents(3) + product%exponents(4))*params%Q)*params%ZetapEtaInv
      T = params%Rho*pq2
      CALL compute_tshpsc_fundamentals(params%Fm, normalization, m_max, params%Rho, T, &
                                       potential_parameter%cutoff_radius)

      prefactor = normalization*attenuation*SQRT((Pi*params%ZetapEtaInv)**3)
      DO i = 1, m_max + 1
         params%Fm(i) = params%Fm(i)*prefactor
      END DO
   END SUBROUTINE prepare_eri4_primitive

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
