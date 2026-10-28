!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Four-center integral evaluation with an encapsulated Libint backend.
!> \par History
!>      11.2006 created [Manuel Guidon]
!>      11.2019 Fixed potential_id initial values (A. Bussy)
!> \author Manuel Guidon
! **************************************************************************************************
MODULE four_center_integral_interface

   USE cell_types,                      ONLY: cell_type,&
                                              real_to_scaled
   USE four_center_integral_preparation,ONLY: prepare_eri4_pair,&
                                              prepare_eri4_primitive,&
                                              prepare_eri4_product
   USE four_center_integral_types,      ONLY: eri4_pgf_image_type,&
                                              eri4_pgf_pair_type,&
                                              eri4_pgf_product_type
   USE hfx_contraction_methods,         ONLY: contract
   USE input_constants,                 ONLY: do_potential_mix_cl_trunc,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE integral_library_types,          ONLY: coulomb_operator_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE libint_wrapper,                  ONLY: &
        cp_libint_cleanup_eri, cp_libint_cleanup_eri1, cp_libint_get_derivs, cp_libint_get_eris, &
        cp_libint_init_eri, cp_libint_init_eri1, cp_libint_set_contrdepth, &
        cp_libint_set_params_eri, cp_libint_set_params_eri_deriv, cp_libint_t
   USE orbital_pointers,                ONLY: nco
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: four_center_integral_init, four_center_integral_cleanup, &
             evaluate_eri, &
             evaluate_deriv_eri, &
             evaluate_eri_screen

!> \brief Context owning the energy and first-derivative backend workspaces.
!> \note Intrinsic copies share backend storage and must not be cleaned up independently.
   TYPE, PUBLIC :: four_center_integral_context_type
      PRIVATE
      TYPE(cp_libint_t) :: eri, deriv
   END TYPE four_center_integral_context_type

CONTAINS

! **************************************************************************************************
!> \brief Initialize the four-center energy and first-derivative workspaces.
!> \param context owning context; global backend initialization is provided by init_cp2k
!> \param max_am maximum angular momentum on any center
! **************************************************************************************************
   SUBROUTINE four_center_integral_init(context, max_am)
      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context
      INTEGER, INTENT(IN)                                :: max_am

      CALL cp_libint_init_eri(context%eri, max_am)
      CALL cp_libint_init_eri1(context%deriv, max_am)
      CALL cp_libint_set_contrdepth(context%eri, 1)
      CALL cp_libint_set_contrdepth(context%deriv, 1)
   END SUBROUTINE four_center_integral_init

! **************************************************************************************************
!> \brief Release the workspaces of an owning four-center context.
!> \param context original initialized owner; borrowed copies must not be passed here
! **************************************************************************************************
   SUBROUTINE four_center_integral_cleanup(context)
      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context

      CALL cp_libint_cleanup_eri(context%eri)
      CALL cp_libint_cleanup_eri1(context%deriv)
   END SUBROUTINE four_center_integral_cleanup

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
!> \brief Evaluates derivatives of electron repulsion integrals for a primitive quartet
!> \param context four-center integral context
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
   SUBROUTINE evaluate_deriv_eri(context, nproducts, pgf_product_list, &
                                 n_a, n_b, n_c, n_d, &
                                 ncoa, ncob, ncoc, ncod, &
                                 nsgfa, nsgfb, nsgfc, nsgfd, &
                                 primitives, max_contraction, tmp_max_all, &
                                 eps_schwarz, neris, &
                                 s_offset_a, s_offset_b, s_offset_c, s_offset_d, &
                                 nl_a, nl_b, nl_c, nl_d, nsoa, nsob, nsoc, nsod, &
                                 sphi_a, sphi_b, sphi_c, sphi_d, &
                                 work, work2, work_forces, buffer1, buffer2, primitives_tmp, &
                                 use_virial, work_virial, work2_virial, primitives_tmp_virial, &
                                 primitives_virial, cell, tmp_max_all_virial)

      TYPE(four_center_integral_context_type)            :: context
      INTEGER, INTENT(IN)                                :: nproducts
      TYPE(eri4_pgf_product_type), DIMENSION(nproducts)  :: pgf_product_list
      INTEGER, INTENT(IN)                                :: n_a, n_b, n_c, n_d, ncoa, ncob, ncoc, &
                                                            ncod, nsgfa, nsgfb, nsgfc, nsgfd
      REAL(dp), &
         DIMENSION(nsgfa, nsgfb, nsgfc, nsgfd, 12)       :: primitives
      REAL(dp)                                           :: max_contraction, tmp_max_all, eps_schwarz
      INTEGER(int_8)                                     :: neris
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
                                                            k, m, m_max, mysize, n, permutation(4)
      LOGICAL                                            :: direct_contract
      REAL(dp)                                           :: scoord(12), tmp_max, tmp_max_virial

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

      DO i = 1, nproducts
         CALL cp_libint_set_params_eri_deriv(context%deriv, pgf_product_list(i), permutation, m_max)
         CALL cp_libint_get_derivs(angular_momenta(4), angular_momenta(3), angular_momenta(2), &
                                   angular_momenta(1), context%deriv, work_forces, a_mysize)
         DO k = 1, 12
            DO j = 1, mysize
               work(j, k) = work(j, k) + work_forces(j, k)
            END DO
         END DO
         neris = neris + 12*mysize
         IF (use_virial) THEN
            DO j = 1, 4
               CALL real_to_scaled(scoord(3*j - 2:3*j), pgf_product_list(i)%centers(:, permutation(j)), cell)
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
!> \param context four-center integral context
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
!> \par History
!>      03.2007 created [Manuel Guidon]
!>      08.2007 refactured permutation part [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE evaluate_eri_screen(context, A, B, C, D, Zeta_A, Zeta_B, Zeta_C, Zeta_D, &
                                  n_a, n_b, n_c, n_d, &
                                  max_val, potential_parameter, R1, R2)

      TYPE(four_center_integral_context_type)            :: context
      REAL(dp), INTENT(IN)                               :: A(3), B(3), C(3), D(3), Zeta_A, Zeta_B, &
                                                            Zeta_C, Zeta_D
      INTEGER, INTENT(IN)                                :: n_a, n_b, n_c, n_d
      REAL(dp), INTENT(INOUT)                            :: max_val
      TYPE(coulomb_operator_type)                        :: potential_parameter
      REAL(dp)                                           :: R1, R2

      INTEGER                                            :: a_mysize(1), angular_momenta(4), i, &
                                                            m_max, mysize, permutation(4)
      REAL(dp)                                           :: centers(3, 4), exponents(4), PQ(3), PQ2, &
                                                            S1234
      REAL(dp), DIMENSION(:), POINTER                    :: p_work
      TYPE(eri4_pgf_image_type)                          :: image_a, image_b
      TYPE(eri4_pgf_pair_type)                           :: pair_a, pair_b
      TYPE(eri4_pgf_product_type)                        :: params, product

      m_max = n_a + n_b + n_c + n_d
      mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)
      a_mysize = mysize

      CALL get_eri_permutation(n_a, n_b, n_c, n_d, permutation)
      angular_momenta = [n_a, n_b, n_c, n_d]
      angular_momenta = angular_momenta(permutation)
      centers(:, 1) = A
      centers(:, 2) = B
      centers(:, 3) = C
      centers(:, 4) = D
      centers = centers(:, permutation)
      exponents = [Zeta_A, Zeta_B, Zeta_C, Zeta_D]
      exponents = exponents(permutation)

      PQ = centers(:, 1) - centers(:, 2)
      CALL prepare_eri4_pair(centers(:, 1), centers(:, 2), exponents(1), exponents(2), &
                             DOT_PRODUCT(PQ, PQ), pair_a, image_a)
      PQ = centers(:, 3) - centers(:, 4)
      CALL prepare_eri4_pair(centers(:, 3), centers(:, 4), exponents(3), exponents(4), &
                             DOT_PRODUCT(PQ, PQ), pair_b, image_b)
      PQ = image_a%P - image_b%P
      PQ2 = DOT_PRODUCT(PQ, PQ)
      SELECT CASE (potential_parameter%potential_type)
      CASE (do_potential_truncated, do_potential_short, do_potential_mix_cl_trunc)
         IF (PQ2 > (R1 + R2 + potential_parameter%cutoff_radius)**2) RETURN
      END SELECT

      CALL prepare_eri4_product(pair_a, pair_b, m_max, product)
      S1234 = EXP(image_a%S1234 + image_b%S1234)
      CALL prepare_eri4_primitive(product, image_a, image_b, S1234, PQ2, potential_parameter, params)
      IF (m_max /= 0) THEN
         CALL cp_libint_set_params_eri(context%eri, params, [1, 2, 3, 4], m_max)
         CALL cp_libint_get_eris(angular_momenta(4), angular_momenta(3), angular_momenta(2), &
                                 angular_momenta(1), context%eri, p_work, a_mysize)
         DO i = 1, mysize
            max_val = MAX(max_val, ABS(p_work(i)))
         END DO
      ELSE
         max_val = ABS(params%Fm(1))
      END IF

   END SUBROUTINE evaluate_eri_screen

! **************************************************************************************************
!> \brief Evaluate electron repulsion integrals for a primitive quartet
!> \param context four-center integral context
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
!> \par History
!>      11.2006 created [Manuel Guidon]
!>      08.2007 refactured permutation part [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE evaluate_eri(context, nproducts, pgf_product_list, &
                           n_a, n_b, n_c, n_d, &
                           ncoa, ncob, ncoc, ncod, &
                           nsgfa, nsgfb, nsgfc, nsgfd, &
                           primitives, max_contraction, tmp_max, &
                           eps_schwarz, neris, &
                           s_offset_a, s_offset_b, s_offset_c, s_offset_d, &
                           nl_a, nl_b, nl_c, nl_d, nsoa, nsob, nsoc, nsod, &
                           sphi_a, sphi_b, sphi_c, sphi_d, work, work2, buffer1, buffer2, &
                           primitives_tmp)

      TYPE(four_center_integral_context_type)            :: context
      INTEGER, INTENT(IN)                                :: nproducts
      TYPE(eri4_pgf_product_type), DIMENSION(nproducts)  :: pgf_product_list
      INTEGER, INTENT(IN)                                :: n_a, n_b, n_c, n_d, ncoa, ncob, ncoc, &
                                                            ncod, nsgfa, nsgfb, nsgfc, nsgfd
      REAL(dp), DIMENSION(nsgfa, nsgfb, nsgfc, nsgfd)    :: primitives
      REAL(dp)                                           :: max_contraction, tmp_max, eps_schwarz
      INTEGER(int_8)                                     :: neris
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

      INTEGER                                            :: a_mysize(1), angular_momenta(4), i, &
                                                            m_max, mysize, permutation(4)
      LOGICAL                                            :: direct_contract
      REAL(dp), DIMENSION(:), POINTER                    :: eri_src, p_work

      IF (nproducts == 0) RETURN

      m_max = n_a + n_b + n_c + n_d
      mysize = ncoa*ncob*ncoc*ncod
      a_mysize = mysize

      IF (m_max /= 0) THEN
         IF (nproducts > 1) work(1:mysize) = 0.0_dp
         CALL get_eri_permutation(n_a, n_b, n_c, n_d, permutation)
         angular_momenta = [n_a, n_b, n_c, n_d]
         angular_momenta = angular_momenta(permutation)

         DO i = 1, nproducts
            CALL cp_libint_set_params_eri(context%eri, pgf_product_list(i), permutation, m_max)
            CALL cp_libint_get_eris(angular_momenta(4), angular_momenta(3), angular_momenta(2), &
                                    angular_momenta(1), context%eri, p_work, a_mysize)
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

END MODULE four_center_integral_interface
