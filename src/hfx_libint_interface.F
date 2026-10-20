!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to the Libint-Library
!> \par History
!>      11.2006 created [Manuel Guidon]
!>      11.2019 Fixed potential_id initial values (A. Bussy)
!> \author Manuel Guidon
! **************************************************************************************************
MODULE hfx_libint_interface

   USE cell_types,                      ONLY: cell_type,&
                                              real_to_scaled
   USE gamma,                           ONLY: fgamma => fgamma_0
   USE hfx_contraction_methods,         ONLY: contract
   USE hfx_types,                       ONLY: hfx_pgf_product_list,&
                                              hfx_potential_type
   USE input_constants,                 ONLY: &
        do_potential_coulomb, do_potential_gaussian, do_potential_id, do_potential_long, &
        do_potential_mix_cl, do_potential_mix_cl_trunc, do_potential_mix_lg, do_potential_short, &
        do_potential_truncated
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE libint_wrapper,                  ONLY: cp_libint_get_derivs,&
                                              cp_libint_get_eris,&
                                              cp_libint_set_params_eri,&
                                              cp_libint_set_params_eri_deriv,&
                                              cp_libint_set_params_eri_screen,&
                                              cp_libint_t,&
                                              get_ssss_f_val,&
                                              prim_data_f_size
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: nco
   USE t_c_g0,                          ONLY: t_c_g0_n
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: evaluate_eri, &
             evaluate_deriv_eri, &
             evaluate_eri_screen

!***

CONTAINS

! **************************************************************************************************
!> \brief Return the canonical quartet order required by the Libint ERI kernels.
!> \param la angular momentum of the original first center
!> \param lb angular momentum of the original second center
!> \param lc angular momentum of the original third center
!> \param ld angular momentum of the original fourth center
!> \param permutation original center index for each canonical center
! **************************************************************************************************
   PURE SUBROUTINE get_eri_permutation(la, lb, lc, ld, permutation)
      INTEGER, INTENT(IN)                                :: la, lb, lc, ld
      INTEGER, INTENT(OUT)                               :: permutation(4)

      permutation = [1, 2, 3, 4]
      IF (la < lb) permutation(1:2) = [2, 1]
      IF (lc < ld) permutation(3:4) = [4, 3]
      IF (la + lb > lc + ld) permutation = permutation([3, 4, 1, 2])
   END SUBROUTINE get_eri_permutation

! **************************************************************************************************
!> \brief Restore original Cartesian axes and center labels from a canonical ERI block.
!> \param source canonical block, with the fourth Cartesian index varying fastest
!> \param destination original block, with the first Cartesian index varying fastest
!> \param ncart Cartesian shell sizes in original center order
!> \param permutation original center index for each canonical center
!> \param ncomponents one integral component or twelve center-derivative components
!> \param nplanes independent planes, one for integrals/forces or three for virial
! **************************************************************************************************
   PURE SUBROUTINE reorder_eri_block(source, destination, ncart, permutation, ncomponents, nplanes)
      INTEGER, INTENT(IN)                                :: ncart(4), permutation(4), ncomponents, &
                                                            nplanes
      REAL(dp), INTENT(OUT) :: destination(PRODUCT(ncart), ncomponents, nplanes)
      REAL(dp), INTENT(IN) :: source(PRODUCT(ncart), ncomponents, nplanes)

      INTEGER                                            :: component, extent(4), i, isource, j, k, &
                                                            l, output_component, p1, p2, p3, p4, &
                                                            plane, stride(4)

      stride = [1, ncart(1), ncart(1)*ncart(2), ncart(1)*ncart(2)*ncart(3)]
      stride = stride(permutation)
      extent = ncart(permutation)
      ! Read the Libint block contiguously and restore the original Fortran layout with strides.
      DO plane = 1, nplanes
         DO component = 1, ncomponents
            output_component = component
            IF (ncomponents == 12) THEN
               output_component = 3*(permutation((component - 1)/3 + 1) - 1) + MOD(component - 1, 3) + 1
            END IF
            isource = 0
            p1 = 1
            DO i = 1, extent(1)
               p2 = p1
               DO j = 1, extent(2)
                  p3 = p2
                  DO k = 1, extent(3)
                     p4 = p3
                     DO l = 1, extent(4)
                        isource = isource + 1
                        destination(p4, output_component, plane) = source(isource, component, plane)
                        p4 = p4 + stride(4)
                     END DO
                     p3 = p3 + stride(3)
                  END DO
                  p2 = p2 + stride(2)
               END DO
               p1 = p1 + stride(1)
            END DO
         END DO
      END DO
   END SUBROUTINE reorder_eri_block

! **************************************************************************************************
!> \brief Fill data structure used in libint
!> \param A ...
!> \param B ...
!> \param C ...
!> \param D ...
!> \param Zeta_A ...
!> \param Zeta_B ...
!> \param Zeta_C ...
!> \param Zeta_D ...
!> \param m_max ...
!> \param potential_parameter ...
!> \param libint ...
!> \param R11 ...
!> \param R22 ...
!> \par History
!>      03.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE build_quartet_data_screen(A, B, C, D, Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max, &
                                        potential_parameter, libint, R11, R22)

      REAL(KIND=dp)                                      :: A(3), B(3), C(3), D(3)
      REAL(KIND=dp), INTENT(IN)                          :: Zeta_A, Zeta_B, Zeta_C, Zeta_D
      INTEGER, INTENT(IN)                                :: m_max
      TYPE(hfx_potential_type)                           :: potential_parameter
      TYPE(cp_libint_t)                                  :: libint
      REAL(dp)                                           :: R11, R22

      INTEGER                                            :: i
      LOGICAL                                            :: use_gamma
      REAL(KIND=dp) :: AB(3), AB2, CD(3), CD2, den, Eta, EtaInv, factor, G(3), num, omega2, &
         omega_corr, omega_corr2, P(3), PQ(3), PQ2, Q(3), R, R1, R2, Rho, RhoInv, S1234, ssss, T, &
         tmp, W(3), Zeta, ZetaInv, ZetapEtaInv
      REAL(KIND=dp), DIMENSION(prim_data_f_size)         :: F, Fm

      Zeta = Zeta_A + Zeta_B
      ZetaInv = 1.0_dp/Zeta
      Eta = Zeta_C + Zeta_D
      EtaInv = 1.0_dp/Eta
      ZetapEtaInv = Zeta + Eta
      ZetapEtaInv = 1.0_dp/ZetapEtaInv
      Rho = Zeta*Eta*ZetapEtaInv
      RhoInv = 1.0_dp/Rho

      DO i = 1, 3
         P(i) = (Zeta_A*A(i) + Zeta_B*B(i))*ZetaInv
         Q(i) = (Zeta_C*C(i) + Zeta_D*D(i))*EtaInv
         AB(i) = A(i) - B(i)
         CD(i) = C(i) - D(i)
         PQ(i) = P(i) - Q(i)
         W(i) = (Zeta*P(i) + Eta*Q(i))*ZetapEtaInv
      END DO

      AB2 = DOT_PRODUCT(AB, AB)
      CD2 = DOT_PRODUCT(CD, CD)
      PQ2 = DOT_PRODUCT(PQ, PQ)

      S1234 = EXP((-Zeta_A*Zeta_B*ZetaInv*AB2) + (-Zeta_C*Zeta_D*EtaInv*CD2))
      T = Rho*PQ2

      SELECT CASE (potential_parameter%potential_type)
      CASE (do_potential_truncated)
         R = potential_parameter%cutoff_radius*SQRT(rho)
         R1 = R11
         R2 = R22
         IF (PQ2 > (R1 + R2 + potential_parameter%cutoff_radius)**2) THEN
            RETURN
         END IF
         CALL t_c_g0_n(F(1), use_gamma, R, T, m_max)
         IF (use_gamma) CALL fgamma(m_max, T, F(1))
         factor = 2.0_dp*Pi*RhoInv
      CASE (do_potential_coulomb)
         CALL fgamma(m_max, T, F(1))
         factor = 2.0_dp*Pi*RhoInv
      CASE (do_potential_short)
         R = potential_parameter%cutoff_radius*SQRT(rho)
         R1 = R11
         R2 = R22
         IF (PQ2 > (R1 + R2 + potential_parameter%cutoff_radius)**2) THEN
            RETURN
         END IF
         CALL fgamma(m_max, T, F(1))
         omega2 = potential_parameter%omega**2
         omega_corr2 = omega2/(omega2 + Rho)
         omega_corr = SQRT(omega_corr2)
         T = T*omega_corr2
         CALL fgamma(m_max, T, Fm)
         tmp = -omega_corr
         DO i = 1, m_max + 1
            F(i) = F(i) + Fm(i)*tmp
            tmp = tmp*omega_corr2
         END DO
         factor = 2.0_dp*Pi*RhoInv
      CASE (do_potential_long)
         omega2 = potential_parameter%omega**2
         omega_corr2 = omega2/(omega2 + Rho)
         omega_corr = SQRT(omega_corr2)
         T = T*omega_corr2
         CALL fgamma(m_max, T, F(1))
         tmp = omega_corr
         DO i = 1, m_max + 1
            F(i) = F(i)*tmp
            tmp = tmp*omega_corr2
         END DO
         factor = 2.0_dp*Pi*RhoInv
      CASE (do_potential_mix_cl)
         CALL fgamma(m_max, T, F(1))
         omega2 = potential_parameter%omega**2
         omega_corr2 = omega2/(omega2 + Rho)
         omega_corr = SQRT(omega_corr2)
         T = T*omega_corr2
         CALL fgamma(m_max, T, Fm)
         tmp = omega_corr
         DO i = 1, m_max + 1
            F(i) = F(i)*potential_parameter%scale_coulomb + Fm(i)*tmp*potential_parameter%scale_longrange
            tmp = tmp*omega_corr2
         END DO
         factor = 2.0_dp*Pi*RhoInv
      CASE (do_potential_mix_cl_trunc)

         ! truncated
         R = potential_parameter%cutoff_radius*SQRT(rho)
         R1 = R11
         R2 = R22
         IF (PQ2 > (R1 + R2 + potential_parameter%cutoff_radius)**2) THEN
            RETURN
         END IF
         CALL t_c_g0_n(F(1), use_gamma, R, T, m_max)
         IF (use_gamma) CALL fgamma(m_max, T, F(1))

         ! Coulomb
         CALL fgamma(m_max, T, Fm)

         DO i = 1, m_max + 1
            F(i) = F(i)*(potential_parameter%scale_coulomb + potential_parameter%scale_longrange) - &
                   Fm(i)*potential_parameter%scale_longrange
         END DO

         ! longrange
         omega2 = potential_parameter%omega**2
         omega_corr2 = omega2/(omega2 + Rho)
         omega_corr = SQRT(omega_corr2)
         T = T*omega_corr2
         CALL fgamma(m_max, T, Fm)
         tmp = omega_corr
         DO i = 1, m_max + 1
            F(i) = F(i) + Fm(i)*tmp*potential_parameter%scale_longrange
            tmp = tmp*omega_corr2
         END DO
         factor = 2.0_dp*Pi*RhoInv

      CASE (do_potential_gaussian)
         omega2 = potential_parameter%omega**2
         T = -omega2*T/(Rho + omega2)
         tmp = 1.0_dp
         DO i = 1, m_max + 1
            F(i) = EXP(T)*tmp
            tmp = tmp*omega2/(Rho + omega2)
         END DO
         factor = (Pi/(Rho + omega2))**(1.5_dp)
      CASE (do_potential_mix_lg)
         omega2 = potential_parameter%omega**2
         omega_corr2 = omega2/(omega2 + Rho)
         omega_corr = SQRT(omega_corr2)
         T = T*omega_corr2
         CALL fgamma(m_max, T, Fm)
         tmp = omega_corr*2.0_dp*Pi*RhoInv*potential_parameter%scale_longrange
         DO i = 1, m_max + 1
            Fm(i) = Fm(i)*tmp
            tmp = tmp*omega_corr2
         END DO
         T = Rho*PQ2
         T = -omega2*T/(Rho + omega2)
         tmp = (Pi/(Rho + omega2))**(1.5_dp)*potential_parameter%scale_gaussian
         DO i = 1, m_max + 1
            F(i) = EXP(T)*tmp + Fm(i)
            tmp = tmp*omega2/(Rho + omega2)
         END DO
         factor = 1.0_dp
      CASE (do_potential_id)
         ssss = -Zeta_A*Zeta_B*ZetaInv*AB2

         num = (Zeta_A + Zeta_B)*Zeta_C
         den = Zeta_A + Zeta_B + Zeta_C
         ssss = ssss - num/den*SUM((P - C)**2)

         G(:) = (Zeta_A*A(:) + Zeta_B*B(:) + Zeta_C*C(:))/den
         num = den*Zeta_D
         den = den + Zeta_D
         ssss = ssss - num/den*SUM((G - D)**2)

         F(:) = EXP(ssss)
         factor = 1.0_dp
         IF (S1234 > EPSILON(0.0_dp)) factor = 1.0_dp/S1234
      END SELECT

      tmp = (Pi*ZetapEtaInv)**3
      factor = factor*S1234*SQRT(tmp)

      DO i = 1, m_max + 1
         F(i) = F(i)*factor
      END DO

      CALL cp_libint_set_params_eri_screen(libint, A, B, C, D, P, Q, W, ZetaInv, EtaInv, ZetapEtaInv, Rho, m_max, F)

   END SUBROUTINE build_quartet_data_screen

! **************************************************************************************************
!> \brief Evaluates derivatives of electron repulsion integrals for a primitive quartet
!> \param deriv ...
!> \param nproducts ...
!> \param pgf_product_list ...
!> \param n_a ...
!> \param n_b ...
!> \param n_c ...
!> \param n_d ...
!> \param ncoa ...
!> \param ncob ...
!> \param ncoc ...
!> \param ncod ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param nsgfc ...
!> \param nsgfd ...
!> \param primitives ...
!> \param max_contraction ...
!> \param tmp_max_all ...
!> \param eps_schwarz ...
!> \param neris ...
!> \param Zeta_A ...
!> \param Zeta_B ...
!> \param Zeta_C ...
!> \param Zeta_D ...
!> \param ZetaInv ...
!> \param EtaInv ...
!> \param s_offset_a ...
!> \param s_offset_b ...
!> \param s_offset_c ...
!> \param s_offset_d ...
!> \param nl_a ...
!> \param nl_b ...
!> \param nl_c ...
!> \param nl_d ...
!> \param nsoa ...
!> \param nsob ...
!> \param nsoc ...
!> \param nsod ...
!> \param sphi_a ...
!> \param sphi_b ...
!> \param sphi_c ...
!> \param sphi_d ...
!> \param work ...
!> \param work2 ...
!> \param work_forces ...
!> \param buffer1 ...
!> \param buffer2 ...
!> \param primitives_tmp ...
!> \param use_virial ...
!> \param work_virial ...
!> \param work2_virial ...
!> \param primitives_tmp_virial ...
!> \param primitives_virial ...
!> \param cell ...
!> \param tmp_max_all_virial ...
!> \par History
!>      03.2007 created [Manuel Guidon]
!>      08.2007 refactured permutation part [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE evaluate_deriv_eri(deriv, nproducts, pgf_product_list, &
                                 n_a, n_b, n_c, n_d, &
                                 ncoa, ncob, ncoc, ncod, &
                                 nsgfa, nsgfb, nsgfc, nsgfd, &
                                 primitives, max_contraction, tmp_max_all, &
                                 eps_schwarz, neris, &
                                 Zeta_A, Zeta_B, Zeta_C, Zeta_D, ZetaInv, EtaInv, &
                                 s_offset_a, s_offset_b, s_offset_c, s_offset_d, &
                                 nl_a, nl_b, nl_c, nl_d, nsoa, nsob, nsoc, nsod, &
                                 sphi_a, sphi_b, sphi_c, sphi_d, &
                                 work, work2, work_forces, buffer1, buffer2, primitives_tmp, &
                                 use_virial, work_virial, work2_virial, primitives_tmp_virial, &
                                 primitives_virial, cell, tmp_max_all_virial)

      TYPE(cp_libint_t)                                  :: deriv
      INTEGER, INTENT(IN)                                :: nproducts
      TYPE(hfx_pgf_product_list), DIMENSION(nproducts)   :: pgf_product_list
      INTEGER, INTENT(IN)                                :: n_a, n_b, n_c, n_d, ncoa, ncob, ncoc, &
                                                            ncod, nsgfa, nsgfb, nsgfc, nsgfd
      REAL(dp), &
         DIMENSION(nsgfa, nsgfb, nsgfc, nsgfd, 12)       :: primitives
      REAL(dp)                                           :: max_contraction, tmp_max_all, eps_schwarz
      INTEGER(int_8)                                     :: neris
      REAL(dp), INTENT(IN)                               :: Zeta_A, Zeta_B, Zeta_C, Zeta_D, ZetaInv, &
                                                            EtaInv
      INTEGER                                            :: s_offset_a, s_offset_b, s_offset_c, &
                                                            s_offset_d, nl_a, nl_b, nl_c, nl_d, &
                                                            nsoa, nsob, nsoc, nsod
      REAL(dp), DIMENSION(ncoa, nsoa*nl_a)               :: sphi_a
      REAL(dp), DIMENSION(ncob, nsob*nl_b)               :: sphi_b
      REAL(dp), DIMENSION(ncoc, nsoc*nl_c)               :: sphi_c
      REAL(dp), DIMENSION(ncod, nsod*nl_d)               :: sphi_d
      REAL(dp) :: work(ncoa*ncob*ncoc*ncod, 12), work2(ncoa, ncob, ncoc, ncod, 12), &
         work_forces(ncoa*ncob*ncoc*ncod, 12)
      REAL(dp), DIMENSION(ncoa*ncob*ncoc*ncod)           :: buffer1, buffer2
      REAL(dp), DIMENSION(nsoa*nl_a, nsob*nl_b, nsoc*&
         nl_c, nsod*nl_d)                                :: primitives_tmp
      LOGICAL, INTENT(IN)                                :: use_virial
      REAL(dp) :: work_virial(ncoa*ncob*ncoc*ncod, 12, 3), &
         work2_virial(ncoa, ncob, ncoc, ncod, 12, 3)
      REAL(dp), DIMENSION(nsoa*nl_a, nsob*nl_b, nsoc*&
         nl_c, nsod*nl_d)                                :: primitives_tmp_virial
      REAL(dp), &
         DIMENSION(nsgfa, nsgfb, nsgfc, nsgfd, 12, 3)    :: primitives_virial
      TYPE(cell_type), POINTER                           :: cell
      REAL(dp)                                           :: tmp_max_all_virial

      INTEGER                                            :: a_mysize(1), angular_momenta(4), i, j, &
                                                            k, m, m_max, mysize, n, pair_order(2), &
                                                            permutation(4)
      LOGICAL                                            :: direct_contract
      REAL(dp) :: centers(3, 4), exponents(4), inverse_exponents(2), product_centers(3, 2), &
         scoord(12), tmp_max, tmp_max_virial

      m_max = n_a + n_b + n_c + n_d
      m_max = m_max + 1
      mysize = ncoa*ncob*ncoc*ncod
      a_mysize = mysize

      work = 0.0_dp

      IF (use_virial) THEN
         work_virial = 0.0_dp
      END IF

      CALL get_eri_permutation(n_a, n_b, n_c, n_d, permutation)
      angular_momenta = [n_a, n_b, n_c, n_d]
      angular_momenta = angular_momenta(permutation)
      exponents = [Zeta_A, Zeta_B, Zeta_C, Zeta_D]
      exponents = exponents(permutation)
      pair_order = [1, 2]
      IF (permutation(1) > 2) pair_order = [2, 1]
      inverse_exponents = [ZetaInv, EtaInv]
      inverse_exponents = inverse_exponents(pair_order)
      DO i = 1, nproducts
         centers(:, 1) = pgf_product_list(i)%ra
         centers(:, 2) = pgf_product_list(i)%rb
         centers(:, 3) = pgf_product_list(i)%rc
         centers(:, 4) = pgf_product_list(i)%rd
         product_centers(:, 1) = pgf_product_list(i)%P
         product_centers(:, 2) = pgf_product_list(i)%Q
         CALL cp_libint_set_params_eri_deriv(deriv, centers(:, permutation(1)), centers(:, permutation(2)), &
                                             centers(:, permutation(3)), centers(:, permutation(4)), &
                                             product_centers(:, pair_order(1)), product_centers(:, pair_order(2)), &
                                             pgf_product_list(i)%W, exponents(1), exponents(2), exponents(3), exponents(4), &
                                             inverse_exponents(1), inverse_exponents(2), pgf_product_list(i)%ZetapEtaInv, &
                                             pgf_product_list(i)%Rho, m_max, pgf_product_list(i)%Fm)
         CALL cp_libint_get_derivs(angular_momenta(4), angular_momenta(3), angular_momenta(2), &
                                   angular_momenta(1), deriv, work_forces, a_mysize)
         DO k = 1, 12
            DO j = 1, mysize
               work(j, k) = work(j, k) + work_forces(j, k)
            END DO
         END DO
         neris = neris + 12*mysize
         IF (use_virial) THEN
            DO j = 1, 4
               CALL real_to_scaled(scoord(3*j - 2:3*j), centers(:, permutation(j)), cell)
            END DO
            DO k = 1, 12
               DO j = 1, mysize
                  DO m = 1, 3
                     work_virial(j, k, m) = work_virial(j, k, m) + work_forces(j, k)*scoord(INT((k - 1)/3)*3 + m)
                  END DO
               END DO
            END DO
         END IF
      END DO

      DO n = 1, 12
         tmp_max = 0.0_dp
         DO i = 1, mysize
            tmp_max = MAX(tmp_max, ABS(work(i, n)))
         END DO
         tmp_max = tmp_max*max_contraction
         tmp_max_all = MAX(tmp_max_all, tmp_max)
      END DO
      CALL reorder_eri_block(work, work2, [ncoa, ncob, ncoc, ncod], permutation, 12, 1)
      IF (use_virial) THEN
         DO n = 1, 12
            tmp_max_virial = 0.0_dp
            DO i = 1, mysize
               tmp_max_virial = MAX(tmp_max_virial, &
                                    ABS(work_virial(i, n, 1)), ABS(work_virial(i, n, 2)), ABS(work_virial(i, n, 3)))
            END DO
            tmp_max_virial = tmp_max_virial*max_contraction
            tmp_max_all_virial = MAX(tmp_max_all_virial, tmp_max_virial)
         END DO
         CALL reorder_eri_block(work_virial, work2_virial, [ncoa, ncob, ncoc, ncod], permutation, 12, 3)
      END IF

      IF (.NOT. use_virial) THEN
         IF (tmp_max_all < eps_schwarz) RETURN
      END IF

      direct_contract = (s_offset_a == 0 .AND. s_offset_b == 0 .AND. &
                         s_offset_c == 0 .AND. s_offset_d == 0 .AND. &
                         nsgfa == nsoa*nl_a .AND. nsgfb == nsob*nl_b .AND. &
                         nsgfc == nsoc*nl_c .AND. nsgfd == nsod*nl_d)

      IF (tmp_max_all >= eps_schwarz) THEN
         IF (direct_contract) THEN
            DO i = 1, 12
               CALL contract(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod, &
                             n_a, n_b, n_c, n_d, nl_a, nl_b, nl_c, nl_d, work2(1, 1, 1, 1, i), &
                             sphi_a, &
                             sphi_b, &
                             sphi_c, &
                             sphi_d, &
                             primitives(1, 1, 1, 1, i), &
                             buffer1, buffer2)
            END DO
         ELSE
            DO i = 1, 12
               primitives_tmp(:, :, :, :) = 0.0_dp
               CALL contract(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod, &
                             n_a, n_b, n_c, n_d, nl_a, nl_b, nl_c, nl_d, work2(1, 1, 1, 1, i), &
                             sphi_a, &
                             sphi_b, &
                             sphi_c, &
                             sphi_d, &
                             primitives_tmp(1, 1, 1, 1), &
                             buffer1, buffer2)

               primitives(s_offset_a + 1:s_offset_a + nsoa*nl_a, &
                          s_offset_b + 1:s_offset_b + nsob*nl_b, &
                          s_offset_c + 1:s_offset_c + nsoc*nl_c, &
                          s_offset_d + 1:s_offset_d + nsod*nl_d, i) = &
                  primitives(s_offset_a + 1:s_offset_a + nsoa*nl_a, &
                             s_offset_b + 1:s_offset_b + nsob*nl_b, &
                             s_offset_c + 1:s_offset_c + nsoc*nl_c, &
                             s_offset_d + 1:s_offset_d + nsod*nl_d, i) + primitives_tmp(:, :, :, :)
            END DO
         END IF
      END IF

      IF (use_virial .AND. tmp_max_all_virial >= eps_schwarz) THEN
         IF (direct_contract) THEN
            DO i = 1, 12
               DO m = 1, 3
                  CALL contract(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod, &
                                n_a, n_b, n_c, n_d, nl_a, nl_b, nl_c, nl_d, work2_virial(1, 1, 1, 1, i, m), &
                                sphi_a, &
                                sphi_b, &
                                sphi_c, &
                                sphi_d, &
                                primitives_virial(1, 1, 1, 1, i, m), &
                                buffer1, buffer2)
               END DO
            END DO
         ELSE
            DO i = 1, 12
               DO m = 1, 3
                  primitives_tmp_virial(:, :, :, :) = 0.0_dp
                  CALL contract(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod, &
                                n_a, n_b, n_c, n_d, nl_a, nl_b, nl_c, nl_d, work2_virial(1, 1, 1, 1, i, m), &
                                sphi_a, &
                                sphi_b, &
                                sphi_c, &
                                sphi_d, &
                                primitives_tmp_virial(1, 1, 1, 1), &
                                buffer1, buffer2)

                  primitives_virial(s_offset_a + 1:s_offset_a + nsoa*nl_a, &
                                    s_offset_b + 1:s_offset_b + nsob*nl_b, &
                                    s_offset_c + 1:s_offset_c + nsoc*nl_c, &
                                    s_offset_d + 1:s_offset_d + nsod*nl_d, i, m) = &
                     primitives_virial(s_offset_a + 1:s_offset_a + nsoa*nl_a, &
                                       s_offset_b + 1:s_offset_b + nsob*nl_b, &
                                       s_offset_c + 1:s_offset_c + nsoc*nl_c, &
                                       s_offset_d + 1:s_offset_d + nsod*nl_d, i, m) + primitives_tmp_virial(:, :, :, :)
               END DO
            END DO
         END IF
      END IF

   END SUBROUTINE evaluate_deriv_eri

! **************************************************************************************************
!> \brief Evaluates max-abs values of  electron repulsion integrals for a primitive quartet
!> \param libint ...
!> \param A ...
!> \param B ...
!> \param C ...
!> \param D ...
!> \param Zeta_A ...
!> \param Zeta_B ...
!> \param Zeta_C ...
!> \param Zeta_D ...
!> \param n_a ...
!> \param n_b ...
!> \param n_c ...
!> \param n_d ...
!> \param max_val ...
!> \param potential_parameter ...
!> \param R1 ...
!> \param R2 ...
!> \param p_work ...
!> \par History
!>      03.2007 created [Manuel Guidon]
!>      08.2007 refactured permutation part [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE evaluate_eri_screen(libint, A, B, C, D, Zeta_A, Zeta_B, Zeta_C, Zeta_D, &
                                  n_a, n_b, n_c, n_d, &
                                  max_val, potential_parameter, R1, R2, &
                                  p_work)

      TYPE(cp_libint_t)                                  :: libint
      REAL(dp), INTENT(IN)                               :: A(3), B(3), C(3), D(3), Zeta_A, Zeta_B, &
                                                            Zeta_C, Zeta_D
      INTEGER, INTENT(IN)                                :: n_a, n_b, n_c, n_d
      REAL(dp), INTENT(INOUT)                            :: max_val
      TYPE(hfx_potential_type)                           :: potential_parameter
      REAL(dp)                                           :: R1, R2
      REAL(dp), DIMENSION(:), POINTER                    :: p_work

      INTEGER                                            :: a_mysize(1), angular_momenta(4), i, &
                                                            m_max, mysize, permutation(4)
      REAL(dp)                                           :: centers(3, 4), exponents(4)

      m_max = n_a + n_b + n_c + n_d
      mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)
      a_mysize = mysize

      IF (m_max /= 0) THEN
         CALL get_eri_permutation(n_a, n_b, n_c, n_d, permutation)
         angular_momenta = [n_a, n_b, n_c, n_d]
         angular_momenta = angular_momenta(permutation)
         centers(:, 1) = A
         centers(:, 2) = B
         centers(:, 3) = C
         centers(:, 4) = D
         exponents = [Zeta_A, Zeta_B, Zeta_C, Zeta_D]
         exponents = exponents(permutation)
         CALL build_quartet_data_screen(centers(:, permutation(1)), centers(:, permutation(2)), &
                                        centers(:, permutation(3)), centers(:, permutation(4)), &
                                        exponents(1), exponents(2), exponents(3), exponents(4), m_max, &
                                        potential_parameter, libint, R1, R2)
         CALL cp_libint_get_eris(angular_momenta(4), angular_momenta(3), angular_momenta(2), &
                                 angular_momenta(1), libint, p_work, a_mysize)
         DO i = 1, mysize
            max_val = MAX(max_val, ABS(p_work(i)))
         END DO
      ELSE
         CALL build_quartet_data_screen(A, B, C, D, Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max, &
                                        potential_parameter, libint, R1, R2)
         max_val = ABS(get_ssss_f_val(libint))
      END IF

   END SUBROUTINE evaluate_eri_screen

! **************************************************************************************************
!> \brief Evaluate electron repulsion integrals for a primitive quartet
!> \param libint ...
!> \param nproducts ...
!> \param pgf_product_list ...
!> \param n_a ...
!> \param n_b ...
!> \param n_c ...
!> \param n_d ...
!> \param ncoa ...
!> \param ncob ...
!> \param ncoc ...
!> \param ncod ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param nsgfc ...
!> \param nsgfd ...
!> \param primitives ...
!> \param max_contraction ...
!> \param tmp_max ...
!> \param eps_schwarz ...
!> \param neris ...
!> \param ZetaInv ...
!> \param EtaInv ...
!> \param s_offset_a ...
!> \param s_offset_b ...
!> \param s_offset_c ...
!> \param s_offset_d ...
!> \param nl_a ...
!> \param nl_b ...
!> \param nl_c ...
!> \param nl_d ...
!> \param nsoa ...
!> \param nsob ...
!> \param nsoc ...
!> \param nsod ...
!> \param sphi_a ...
!> \param sphi_b ...
!> \param sphi_c ...
!> \param sphi_d ...
!> \param work ...
!> \param work2 ...
!> \param buffer1 ...
!> \param buffer2 ...
!> \param primitives_tmp ...
!> \param p_work ...
!> \par History
!>      11.2006 created [Manuel Guidon]
!>      08.2007 refactured permutation part [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE evaluate_eri(libint, nproducts, pgf_product_list, &
                           n_a, n_b, n_c, n_d, &
                           ncoa, ncob, ncoc, ncod, &
                           nsgfa, nsgfb, nsgfc, nsgfd, &
                           primitives, max_contraction, tmp_max, &
                           eps_schwarz, neris, &
                           ZetaInv, EtaInv, &
                           s_offset_a, s_offset_b, s_offset_c, s_offset_d, &
                           nl_a, nl_b, nl_c, nl_d, nsoa, nsob, nsoc, nsod, &
                           sphi_a, sphi_b, sphi_c, sphi_d, work, work2, buffer1, buffer2, &
                           primitives_tmp, p_work)

      TYPE(cp_libint_t)                                  :: libint
      INTEGER, INTENT(IN)                                :: nproducts
      TYPE(hfx_pgf_product_list), DIMENSION(nproducts)   :: pgf_product_list
      INTEGER, INTENT(IN)                                :: n_a, n_b, n_c, n_d, ncoa, ncob, ncoc, &
                                                            ncod, nsgfa, nsgfb, nsgfc, nsgfd
      REAL(dp), DIMENSION(nsgfa, nsgfb, nsgfc, nsgfd)    :: primitives
      REAL(dp)                                           :: max_contraction, tmp_max, eps_schwarz
      INTEGER(int_8)                                     :: neris
      REAL(dp), INTENT(IN)                               :: ZetaInv, EtaInv
      INTEGER                                            :: s_offset_a, s_offset_b, s_offset_c, &
                                                            s_offset_d, nl_a, nl_b, nl_c, nl_d, &
                                                            nsoa, nsob, nsoc, nsod
      REAL(dp), DIMENSION(ncoa, nsoa*nl_a), INTENT(IN)   :: sphi_a
      REAL(dp), DIMENSION(ncob, nsob*nl_b), INTENT(IN)   :: sphi_b
      REAL(dp), DIMENSION(ncoc, nsoc*nl_c), INTENT(IN)   :: sphi_c
      REAL(dp), DIMENSION(ncod, nsod*nl_d), INTENT(IN)   :: sphi_d
      REAL(dp), TARGET                                   :: work(ncoa*ncob*ncoc*ncod)
      REAL(dp)                                           :: work2(ncoa, ncob, ncoc, ncod)
      REAL(dp), DIMENSION(ncoa*ncob*ncoc*ncod)           :: buffer1, buffer2
      REAL(dp), DIMENSION(nsoa*nl_a, nsob*nl_b, nsoc*&
         nl_c, nsod*nl_d)                                :: primitives_tmp
      REAL(dp), DIMENSION(:), POINTER                    :: p_work

      INTEGER                                            :: a_mysize(1), angular_momenta(4), i, &
                                                            m_max, mysize, pair_order(2), &
                                                            permutation(4)
      LOGICAL                                            :: direct_contract
      REAL(dp)                                           :: centers(3, 4), inverse_exponents(2), &
                                                            product_centers(3, 2)
      REAL(dp), DIMENSION(:), POINTER                    :: eri_src

      IF (nproducts == 0) RETURN

      m_max = n_a + n_b + n_c + n_d
      mysize = ncoa*ncob*ncoc*ncod
      a_mysize = mysize

      IF (m_max /= 0) THEN
         IF (nproducts > 1) work(1:mysize) = 0.0_dp
         CALL get_eri_permutation(n_a, n_b, n_c, n_d, permutation)
         angular_momenta = [n_a, n_b, n_c, n_d]
         angular_momenta = angular_momenta(permutation)
         pair_order = [1, 2]
         IF (permutation(1) > 2) pair_order = [2, 1]
         inverse_exponents = [ZetaInv, EtaInv]
         inverse_exponents = inverse_exponents(pair_order)
         DO i = 1, nproducts
            centers(:, 1) = pgf_product_list(i)%ra
            centers(:, 2) = pgf_product_list(i)%rb
            centers(:, 3) = pgf_product_list(i)%rc
            centers(:, 4) = pgf_product_list(i)%rd
            product_centers(:, 1) = pgf_product_list(i)%P
            product_centers(:, 2) = pgf_product_list(i)%Q
            CALL cp_libint_set_params_eri(libint, centers(:, permutation(1)), centers(:, permutation(2)), &
                                          centers(:, permutation(3)), centers(:, permutation(4)), &
                                          inverse_exponents(1), inverse_exponents(2), pgf_product_list(i)%ZetapEtaInv, &
                                          pgf_product_list(i)%Rho, product_centers(:, pair_order(1)), &
                                          product_centers(:, pair_order(2)), pgf_product_list(i)%W, m_max, pgf_product_list(i)%Fm)
            CALL cp_libint_get_eris(angular_momenta(4), angular_momenta(3), angular_momenta(2), &
                                    angular_momenta(1), libint, p_work, a_mysize)
            IF (nproducts > 1) work(1:mysize) = work(1:mysize) + p_work(1:mysize)
            neris = neris + mysize
         END DO
         IF (nproducts == 1) THEN
            eri_src => p_work
         ELSE
            eri_src => work
         END IF
         DO i = 1, mysize
            tmp_max = MAX(tmp_max, ABS(eri_src(i)))
         END DO
         tmp_max = tmp_max*max_contraction
         IF (tmp_max < eps_schwarz) RETURN
         CALL reorder_eri_block(eri_src, work2, [ncoa, ncob, ncoc, ncod], permutation, 1, 1)
      ELSE
         work(1) = 0.0_dp
         DO i = 1, nproducts
            work(1) = work(1) + pgf_product_list(i)%Fm(1)
            neris = neris + mysize
         END DO
         work2(1, 1, 1, 1) = work(1)
         tmp_max = max_contraction*ABS(work(1))
         IF (tmp_max < eps_schwarz) RETURN
      END IF

      direct_contract = (s_offset_a == 0 .AND. s_offset_b == 0 .AND. &
                         s_offset_c == 0 .AND. s_offset_d == 0 .AND. &
                         nsgfa == nsoa*nl_a .AND. nsgfb == nsob*nl_b .AND. &
                         nsgfc == nsoc*nl_c .AND. nsgfd == nsod*nl_d)
      IF (direct_contract) THEN
         CALL contract(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod, &
                       n_a, n_b, n_c, n_d, nl_a, nl_b, nl_c, nl_d, work2, &
                       sphi_a, &
                       sphi_b, &
                       sphi_c, &
                       sphi_d, &
                       primitives, &
                       buffer1, buffer2)
      ELSE
         primitives_tmp = 0.0_dp

         CALL contract(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod, &
                       n_a, n_b, n_c, n_d, nl_a, nl_b, nl_c, nl_d, work2, &
                       sphi_a, &
                       sphi_b, &
                       sphi_c, &
                       sphi_d, &
                       primitives_tmp, &
                       buffer1, buffer2)

         primitives(s_offset_a + 1:s_offset_a + nsoa*nl_a, &
                    s_offset_b + 1:s_offset_b + nsob*nl_b, &
                    s_offset_c + 1:s_offset_c + nsoc*nl_c, &
                    s_offset_d + 1:s_offset_d + nsod*nl_d) = &
            primitives(s_offset_a + 1:s_offset_a + nsoa*nl_a, &
                       s_offset_b + 1:s_offset_b + nsob*nl_b, &
                       s_offset_c + 1:s_offset_c + nsoc*nl_c, &
                       s_offset_d + 1:s_offset_d + nsod*nl_d) + primitives_tmp(:, :, :, :)
      END IF

   END SUBROUTINE evaluate_eri

END MODULE hfx_libint_interface
