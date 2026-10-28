!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Gaussian-product geometry and radial seeds shared by four-center evaluation and screening.
! **************************************************************************************************
MODULE four_center_integral_preparation
   USE four_center_integral_types,      ONLY: eri4_fundamentals_size,&
                                              eri4_pgf_image_type,&
                                              eri4_pgf_pair_type,&
                                              eri4_pgf_product_type
   USE input_constants,                 ONLY: do_potential_id
   USE integral_library_types,          ONLY: coulomb_operator_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE r12_fundamentals,                ONLY: compute_r12_fundamentals
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: prepare_eri4_pair, prepare_eri4_product, prepare_eri4_primitive

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
      image%rab2 = rab2
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
      TYPE(eri4_pgf_product_type), INTENT(INOUT)         :: product

      CPASSERT(m_max < eri4_fundamentals_size)

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
      TYPE(eri4_pgf_product_type), INTENT(IN)            :: product
      TYPE(eri4_pgf_image_type), INTENT(IN)              :: image_a, image_b
      REAL(dp), INTENT(IN)                               :: attenuation, pq2
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter
      TYPE(eri4_pgf_product_type), INTENT(INOUT)         :: params

      INTEGER                                            :: i, m_max
      REAL(dp) :: den, fundamentals_work(eri4_fundamentals_size), G(3), normalization, num, &
         prefactor, ssss, T

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

      ASSOCIATE (A => image_a%ra, B => image_a%rb, C => image_b%ra, D => image_b%rb, &
                 zeta_a => product%exponents(1), zeta_b => product%exponents(2), &
                 zeta_c => product%exponents(3), zeta_d => product%exponents(4))
         params%W = ((zeta_a + zeta_b)*params%P + (zeta_c + zeta_d)*params%Q)*params%ZetapEtaInv
         T = params%Rho*pq2
         IF (potential_parameter%potential_type == do_potential_id) THEN
            num = zeta_a*zeta_b
            den = zeta_a + zeta_b
            ssss = -num/den*SUM((A - B)**2)
            num = den*zeta_c
            den = den + zeta_c
            ssss = ssss - num/den*SUM((params%P - C)**2)
            G = (zeta_a*A + zeta_b*B + zeta_c*C)/den
            num = den*zeta_d
            den = den + zeta_d
            ssss = ssss - num/den*SUM((G - D)**2)
            params%Fm(:) = EXP(ssss)
            ! The ID seed already contains the four-Gaussian attenuation.
            normalization = 1.0_dp
            IF (attenuation > EPSILON(0.0_dp)) normalization = 1.0_dp/attenuation
         ELSE
            CALL compute_r12_fundamentals(params%Fm, normalization, m_max, potential_parameter%potential_type, &
                                          params%Rho, T, fundamentals_work, &
                                          omega=potential_parameter%omega, &
                                          cutoff_radius=potential_parameter%cutoff_radius, &
                                          scale_coulomb=potential_parameter%scale_coulomb, &
                                          scale_longrange=potential_parameter%scale_longrange, &
                                          scale_gaussian=potential_parameter%scale_gaussian)
         END IF
      END ASSOCIATE

      prefactor = normalization*attenuation*SQRT((Pi*params%ZetapEtaInv)**3)
      DO i = 1, m_max + 1
         params%Fm(i) = params%Fm(i)*prefactor
      END DO
   END SUBROUTINE prepare_eri4_primitive

END MODULE four_center_integral_preparation
