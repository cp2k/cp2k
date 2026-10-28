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
                                              pbc,&
                                              real_to_scaled
   USE four_center_integral_preparation,ONLY: min_image_gap_sq,&
                                              prepare_eri4_pair,&
                                              prepare_eri4_primitive,&
                                              prepare_eri4_product
   USE four_center_integral_types,      ONLY: eri4_cell_type,&
                                              eri4_pgf_image_type,&
                                              eri4_pgf_pair_type
   USE hfx_contraction_methods,         ONLY: contract
   USE input_constants,                 ONLY: do_potential_TShPSC,&
                                              do_potential_mix_cl_trunc,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE integral_library_types,          ONLY: coulomb_operator_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE libint_wrapper,                  ONLY: &
        cp_libint_cleanup_eri, cp_libint_cleanup_eri1, cp_libint_get_derivs, cp_libint_get_eris, &
        cp_libint_init_eri, cp_libint_init_eri1, cp_libint_set_contrdepth, &
        cp_libint_set_params_eri, cp_libint_set_params_eri_deriv, cp_libint_t, &
        eri_primitive_data_type, get_eri_permutation, reorder_eri_block
   USE orbital_pointers,                ONLY: nco
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: four_center_integral_init, four_center_integral_cleanup, &
             prepare_four_center_batch, evaluate_eri, &
             evaluate_deriv_eri, &
             evaluate_eri_screen

!> \brief Context owning the energy and first-derivative backend workspaces.
!> \note Each thread owns its context. Pass it directly without intrinsic copying.
   TYPE, PUBLIC :: four_center_integral_context_type
      PRIVATE
      TYPE(cp_libint_t) :: eri, deriv
      TYPE(eri_primitive_data_type), ALLOCATABLE :: products(:)
      INTEGER :: nproducts = 0
   END TYPE four_center_integral_context_type

CONTAINS

! **************************************************************************************************
!> \brief Select image-resolved quartets and prepare their backend cache before angular traversal.
!> \param context owning per-thread geometry and backend workspace
!> \param list1 bra primitive pair with image topology and screening bounds
!> \param list2 ket primitive pair with image topology and screening bounds
!> \param nproducts number of retained image quartets
!> \param log10_pmax logarithmic density bound
!> \param log10_eps_schwarz logarithmic integral screening threshold
!> \param neighbor_cells periodic image translations
!> \param cell simulation cell
!> \param potential_parameter physical radial operator
!> \param m_max maximum radial order including any derivative increment
!> \param do_periodic whether to enumerate periodic translations
! **************************************************************************************************
   SUBROUTINE prepare_four_center_batch(context, list1, list2, nproducts, &
                                        log10_pmax, log10_eps_schwarz, neighbor_cells, &
                                        cell, potential_parameter, m_max, do_periodic)

      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context
      TYPE(eri4_pgf_pair_type)                           :: list1, list2
      INTEGER, INTENT(OUT)                               :: nproducts
      REAL(dp), INTENT(IN)                               :: log10_pmax, log10_eps_schwarz
      TYPE(eri4_cell_type), DIMENSION(:), POINTER        :: neighbor_cells
      TYPE(cell_type), POINTER                           :: cell
      TYPE(coulomb_operator_type)                        :: potential_parameter
      INTEGER, INTENT(IN)                                :: m_max
      LOGICAL, INTENT(IN)                                :: do_periodic

      REAL(dp), PARAMETER                                :: r_cutoff_eps = 1.0E-12_dp

      INTEGER                                            :: i, j, k, k_max, nimages1, nimages2, &
                                                            tmp_i4
      LOGICAL                                            :: has_r_cutoff
      REAL(dp) :: C11(3), EtaInv, min_img_delta, P(3), pgf_max_1, pgf_max_2, PQ(3), pq2, Q(3), R1, &
         R2, r_cutoff_sq, r_cutoff_sq_guard, rc(3), rcut_eff, rd(3), rpq2, S1234, S1234a, S1234b, &
         shift(3), temp(3), temp_CC(3), temp_DD(3), tmp_D(3), Zeta_C, Zeta_D
      TYPE(eri4_pgf_image_type)                          :: shifted_image
      TYPE(eri_primitive_data_type)                      :: product
      TYPE(eri_primitive_data_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: tmp_product_list

      nimages1 = list1%nimages
      nimages2 = list2%nimages
      nproducts = 0
      context%nproducts = 0
      IF (nimages1 == 0 .OR. nimages2 == 0) RETURN
      CALL prepare_eri4_product(list1, list2, m_max, product)
      EtaInv = list2%ZetaInv
      Zeta_C = list2%zeta
      Zeta_D = list2%zetb
      temp_CC = 0.0_dp
      temp_DD = 0.0_dp
      r_cutoff_sq = 0.0_dp
      r_cutoff_sq_guard = 0.0_dp
      SELECT CASE (potential_parameter%potential_type)
      CASE (do_potential_truncated, do_potential_short, do_potential_mix_cl_trunc)
         has_r_cutoff = .TRUE.
         rcut_eff = potential_parameter%cutoff_radius
      CASE (do_potential_TShPSC)
         has_r_cutoff = .TRUE.
         rcut_eff = potential_parameter%cutoff_radius*2.0_dp
      CASE DEFAULT
         has_r_cutoff = .FALSE.
         rcut_eff = 0.0_dp
      END SELECT
      DO i = 1, nimages1
         P = list1%image_list(i)%P
         R1 = list1%image_list(i)%R
         S1234a = list1%image_list(i)%S1234
         pgf_max_1 = list1%image_list(i)%pgf_max
         DO j = 1, nimages2
            pgf_max_2 = list2%image_list(j)%pgf_max
            IF (pgf_max_1 + pgf_max_2 + log10_pmax < log10_eps_schwarz) CYCLE
            Q = list2%image_list(j)%P
            R2 = list2%image_list(j)%R
            S1234b = list2%image_list(j)%S1234
            rc = list2%image_list(j)%ra
            rd = list2%image_list(j)%rb

            IF (has_r_cutoff) THEN
               r_cutoff_sq = (R1 + R2 + rcut_eff)**2
               r_cutoff_sq_guard = r_cutoff_sq*(1.0_dp + r_cutoff_eps)
            END IF
            k_max = SIZE(neighbor_cells)
            IF (do_periodic) THEN
               temp = P - Q
               PQ = pbc(temp, cell)
               IF (has_r_cutoff .AND. cell%orthorhombic .AND. ALL(cell%perd == 1) .AND. k_max >= 1) THEN
                  IF (ALL(neighbor_cells(1)%cell_r == 0.0_dp)) THEN
                     pq2 = PQ(1)**2 + PQ(2)**2 + PQ(3)**2
                     IF (pq2 > r_cutoff_sq_guard) CYCLE
                     min_img_delta = min_image_gap_sq(PQ, cell%hmat)
                     IF (min_img_delta >= 0.0_dp .AND. pq2 + min_img_delta > r_cutoff_sq_guard) k_max = 1
                  END IF
               END IF
               shift = -PQ + temp
               temp_CC = rc + shift
               temp_DD = rd + shift
            END IF

            S1234 = EXP(S1234a + S1234b)

            DO k = 1, k_max
               IF (do_periodic) THEN
                  IF (has_r_cutoff .AND. k_max > 1) THEN
                     IF ((PQ(1) - neighbor_cells(k)%cell_r(1))**2 + &
                         (PQ(2) - neighbor_cells(k)%cell_r(2))**2 + &
                         (PQ(3) - neighbor_cells(k)%cell_r(3))**2 > r_cutoff_sq_guard) CYCLE
                  END IF
                  C11 = temp_CC + neighbor_cells(k)%cell_r(:)
                  tmp_D = temp_DD + neighbor_cells(k)%cell_r(:)
               ELSE
                  C11 = rc
                  tmp_D = rd
               END IF
               Q = (Zeta_C*C11 + Zeta_D*tmp_D)*EtaInv
               rpq2 = (P(1) - Q(1))**2 + (P(2) - Q(2))**2 + (P(3) - Q(3))**2
               IF (has_r_cutoff) THEN
                  IF (rpq2 > r_cutoff_sq) CYCLE
               END IF
               nproducts = nproducts + 1

               ! Retain the enlarged buffer in this thread's owning context.
               IF (nproducts > SIZE(context%products)) THEN
                  tmp_i4 = (3*nproducts + 1)/2
                  ALLOCATE (tmp_product_list(tmp_i4))
                  tmp_product_list(1:SIZE(context%products)) = context%products
                  CALL MOVE_ALLOC(tmp_product_list, context%products)
               END IF

               shifted_image%ra = C11
               shifted_image%rb = tmp_D
               shifted_image%P = Q
               CALL prepare_eri4_primitive(product, list1%image_list(i), shifted_image, S1234, rpq2, potential_parameter, &
                                           context%products(nproducts))
            END DO
         END DO
      END DO

      context%nproducts = nproducts
   END SUBROUTINE prepare_four_center_batch

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
      ALLOCATE (context%products(128))
   END SUBROUTINE four_center_integral_init

! **************************************************************************************************
!> \brief Release the workspaces of an owning four-center context.
!> \param context original initialized owner
! **************************************************************************************************
   SUBROUTINE four_center_integral_cleanup(context)
      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context

      CALL cp_libint_cleanup_eri(context%eri)
      CALL cp_libint_cleanup_eri1(context%deriv)
      DEALLOCATE (context%products)
      context%nproducts = 0
   END SUBROUTINE four_center_integral_cleanup

! **************************************************************************************************
!> \brief Evaluates derivatives of electron repulsion integrals for a primitive quartet
!> \param context four-center integral context
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
   SUBROUTINE evaluate_deriv_eri(context, &
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

      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context
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

      INTEGER                                            :: angular_momenta(4), i, j, k, m, mysize, &
                                                            permutation(4)
      LOGICAL                                            :: direct_contract
      REAL(dp)                                           :: scoord(12)

      mysize = ncoa*ncob*ncoc*ncod

      work = 0.0_dp

      IF (use_virial) THEN
         work_virial = 0.0_dp
      END IF

      CALL get_eri_permutation(n_a, n_b, n_c, n_d, permutation)
      angular_momenta = [n_a, n_b, n_c, n_d]
      angular_momenta = angular_momenta(permutation)

      DO i = 1, context%nproducts
         CALL cp_libint_set_params_eri_deriv(context%deriv, context%products(i), &
                                             permutation, SUM(angular_momenta) + 1)
         CALL cp_libint_get_derivs(angular_momenta(4), angular_momenta(3), angular_momenta(2), &
                                   angular_momenta(1), context%deriv, work_forces, [mysize])
         DO k = 1, 12
            DO j = 1, mysize
               work(j, k) = work(j, k) + work_forces(j, k)
            END DO
         END DO
         neris = neris + 12*mysize
         IF (use_virial) THEN
            DO j = 1, 4
               CALL real_to_scaled(scoord(3*j - 2:3*j), context%products(i)%centers(:, permutation(j)), cell)
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

      tmp_max_all = MAX(tmp_max_all, max_contraction*MAXVAL(ABS(work)))
      IF (use_virial) THEN
         tmp_max_all_virial = MAX(tmp_max_all_virial, max_contraction*MAXVAL(ABS(work_virial)))
      END IF

      direct_contract = (s_offset_a == 0 .AND. s_offset_b == 0 .AND. &
                         s_offset_c == 0 .AND. s_offset_d == 0 .AND. &
                         nsgfa == nsoa*nl_a .AND. nsgfb == nsob*nl_b .AND. &
                         nsgfc == nsoc*nl_c .AND. nsgfd == nsod*nl_d)

      IF (tmp_max_all >= eps_schwarz) THEN
         CALL reorder_eri_block(work, work2, [ncoa, ncob, ncoc, ncod], permutation, 12, 1)
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
         CALL reorder_eri_block(work_virial, work2_virial, [ncoa, ncob, ncoc, ncod], permutation, 12, 3)
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

      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context
      REAL(dp), INTENT(IN)                               :: A(3), B(3), C(3), D(3), Zeta_A, Zeta_B, &
                                                            Zeta_C, Zeta_D
      INTEGER, INTENT(IN)                                :: n_a, n_b, n_c, n_d
      REAL(dp), INTENT(INOUT)                            :: max_val
      TYPE(coulomb_operator_type)                        :: potential_parameter
      REAL(dp)                                           :: R1, R2

      INTEGER                                            :: angular_momenta(4), i, m_max, mysize, &
                                                            permutation(4)
      REAL(dp)                                           :: centers(3, 4), exponents(4), PQ(3), PQ2, &
                                                            S1234
      REAL(dp), POINTER                                  :: p_work(:)
      TYPE(eri4_pgf_image_type)                          :: image_a, image_b
      TYPE(eri4_pgf_pair_type)                           :: pair_a, pair_b
      TYPE(eri_primitive_data_type)                      :: product

      m_max = n_a + n_b + n_c + n_d
      mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)

      angular_momenta = [n_a, n_b, n_c, n_d]
      centers(:, 1) = A
      centers(:, 2) = B
      centers(:, 3) = C
      centers(:, 4) = D
      exponents = [Zeta_A, Zeta_B, Zeta_C, Zeta_D]

      CALL get_eri_permutation(n_a, n_b, n_c, n_d, permutation)
      angular_momenta = angular_momenta(permutation)
      centers = centers(:, permutation)
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
      CALL prepare_eri4_primitive(product, image_a, image_b, S1234, PQ2, potential_parameter, context%products(1))
      IF (m_max == 0) THEN
         max_val = ABS(context%products(1)%Fm(1))
      ELSE
         CALL cp_libint_set_params_eri(context%eri, context%products(1), [1, 2, 3, 4], m_max)
         CALL cp_libint_get_eris(angular_momenta(4), angular_momenta(3), angular_momenta(2), &
                                 angular_momenta(1), context%eri, p_work, [mysize])
         DO i = 1, mysize
            max_val = MAX(max_val, ABS(p_work(i)))
         END DO
      END IF

   END SUBROUTINE evaluate_eri_screen

! **************************************************************************************************
!> \brief Evaluate electron repulsion integrals for a primitive quartet
!> \param context four-center integral context
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
!> \param work2 ...
!> \param buffer1 ...
!> \param buffer2 ...
!> \param primitives_tmp ...
!> \par History
!>      11.2006 created [Manuel Guidon]
!>      08.2007 refactured permutation part [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE evaluate_eri(context, &
                           n_a, n_b, n_c, n_d, &
                           ncoa, ncob, ncoc, ncod, &
                           nsgfa, nsgfb, nsgfc, nsgfd, &
                           primitives, max_contraction, tmp_max, &
                           eps_schwarz, neris, &
                           s_offset_a, s_offset_b, s_offset_c, s_offset_d, &
                           nl_a, nl_b, nl_c, nl_d, nsoa, nsob, nsoc, nsod, &
                           sphi_a, sphi_b, sphi_c, sphi_d, work2, buffer1, buffer2, &
                           primitives_tmp)

      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context
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
      REAL(dp)                                           :: work2(ncoa, ncob, ncoc, ncod)
      REAL(dp), DIMENSION(ncoa*ncob*ncoc*ncod), TARGET   :: buffer1
      REAL(dp), DIMENSION(ncoa*ncob*ncoc*ncod)           :: buffer2
      REAL(dp), DIMENSION(nsoa*nl_a, nsob*nl_b, nsoc*&
         nl_c, nsod*nl_d)                                :: primitives_tmp

      INTEGER                                            :: angular_momenta(4), i, m_max, mysize, &
                                                            permutation(4)
      LOGICAL                                            :: direct_contract
      REAL(dp), POINTER                                  :: p_work(:)

      IF (context%nproducts == 0) RETURN
      mysize = ncoa*ncob*ncoc*ncod
      m_max = n_a + n_b + n_c + n_d
      CALL get_eri_permutation(n_a, n_b, n_c, n_d, permutation)
      angular_momenta = [n_a, n_b, n_c, n_d]
      angular_momenta = angular_momenta(permutation)
      IF (m_max /= 0) THEN
         IF (context%nproducts > 1) buffer1 = 0.0_dp
         DO i = 1, context%nproducts
            CALL cp_libint_set_params_eri(context%eri, context%products(i), permutation, m_max)
            CALL cp_libint_get_eris(angular_momenta(4), angular_momenta(3), angular_momenta(2), &
                                    angular_momenta(1), context%eri, p_work, [mysize])
            IF (context%nproducts > 1) buffer1 = buffer1 + p_work(1:mysize)
         END DO
         IF (context%nproducts > 1) p_work => buffer1
         tmp_max = MAX(tmp_max, MAXVAL(ABS(p_work)))
      ELSE
         work2 = 0.0_dp
         DO i = 1, context%nproducts
            work2(1, 1, 1, 1) = work2(1, 1, 1, 1) + context%products(i)%Fm(1)
         END DO
         tmp_max = MAX(tmp_max, ABS(work2(1, 1, 1, 1)))
      END IF
      neris = neris + INT(context%nproducts, int_8)*mysize
      tmp_max = tmp_max*max_contraction
      IF (tmp_max < eps_schwarz) RETURN
      IF (m_max /= 0) CALL reorder_eri_block(p_work, work2, [ncoa, ncob, ncoc, ncod], permutation, 1, 1)

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
