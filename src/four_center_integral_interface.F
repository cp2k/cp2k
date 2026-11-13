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
                                              prepare_eri4_pair
   USE four_center_integral_types,      ONLY: eri4_cell_type,&
                                              eri4_pgf_image_type,&
                                              eri4_pgf_pair_type
   USE hfx_contraction_methods,         ONLY: contract
   USE input_constants,                 ONLY: do_potential_mix_cl_trunc,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE integral_library_types,          ONLY: coulomb_operator_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE libint_engine,                   ONLY: cp_libint_engine_context_type,&
                                              cp_libint_engine_evaluate,&
                                              cp_libint_engine_init,&
                                              cp_libint_engine_release
   USE orbital_pointers,                ONLY: nco
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: four_center_integral_init, four_center_integral_cleanup, &
             prepare_four_center_batch, evaluate_eri, &
             evaluate_deriv_eri, &
             evaluate_eri_screen

!> \brief Shell centers and exponents for a selected image quartet.
   TYPE :: eri4_shell_product_type
      REAL(dp) :: centers(3, 4) = 0.0_dp, exponents(4) = 0.0_dp
   END TYPE eri4_shell_product_type

!> \brief Context owning the energy and first-derivative backend workspaces.
!> \note Each thread owns its context. Pass it directly without intrinsic copying.
   TYPE, PUBLIC :: four_center_integral_context_type
      PRIVATE
      TYPE(cp_libint_engine_context_type) :: eri = cp_libint_engine_context_type(), &
                                             deriv = cp_libint_engine_context_type()
      TYPE(eri4_shell_product_type), ALLOCATABLE :: products(:)
      TYPE(coulomb_operator_type) :: potential = coulomb_operator_type()
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
!> \param do_periodic whether to enumerate periodic translations
! **************************************************************************************************
   SUBROUTINE prepare_four_center_batch(context, list1, list2, nproducts, &
                                        log10_pmax, log10_eps_schwarz, neighbor_cells, &
                                        cell, potential_parameter, do_periodic)

      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context
      TYPE(eri4_pgf_pair_type)                           :: list1, list2
      INTEGER, INTENT(OUT)                               :: nproducts
      REAL(dp), INTENT(IN)                               :: log10_pmax, log10_eps_schwarz
      TYPE(eri4_cell_type), DIMENSION(:), POINTER        :: neighbor_cells
      TYPE(cell_type), POINTER                           :: cell
      TYPE(coulomb_operator_type)                        :: potential_parameter
      LOGICAL, INTENT(IN)                                :: do_periodic

      REAL(dp), PARAMETER                                :: r_cutoff_eps = 1.0E-12_dp

      INTEGER                                            :: i, j, k, k_max, nimages1, nimages2, &
                                                            tmp_i4
      LOGICAL                                            :: has_r_cutoff
      REAL(dp) :: C11(3), EtaInv, min_img_delta, P(3), pgf_max_1, pgf_max_2, PQ(3), pq2, Q(3), R1, &
         R2, r_cutoff_sq, r_cutoff_sq_guard, rc(3), rcut_eff, rd(3), rpq2, shift(3), temp(3), &
         temp_CC(3), temp_DD(3), tmp_D(3), Zeta_C, Zeta_D
      TYPE(eri4_shell_product_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: tmp_product_list

      nimages1 = list1%nimages
      nimages2 = list2%nimages
      nproducts = 0
      context%nproducts = 0
      context%potential = potential_parameter
      IF (nimages1 == 0 .OR. nimages2 == 0) RETURN
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
      CASE DEFAULT
         has_r_cutoff = .FALSE.
         rcut_eff = 0.0_dp
      END SELECT
      DO i = 1, nimages1
         P = list1%image_list(i)%P
         R1 = list1%image_list(i)%R
         pgf_max_1 = list1%image_list(i)%pgf_max
         DO j = 1, nimages2
            pgf_max_2 = list2%image_list(j)%pgf_max
            IF (pgf_max_1 + pgf_max_2 + log10_pmax < log10_eps_schwarz) CYCLE
            Q = list2%image_list(j)%P
            R2 = list2%image_list(j)%R
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

               context%products(nproducts)%centers(:, 1) = list1%image_list(i)%ra
               context%products(nproducts)%centers(:, 2) = list1%image_list(i)%rb
               context%products(nproducts)%centers(:, 3) = C11
               context%products(nproducts)%centers(:, 4) = tmp_D
               context%products(nproducts)%exponents = [list1%zeta, list1%zetb, list2%zeta, list2%zetb]
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

      CALL cp_libint_engine_init(context%eri, 4, max_am, 0)
      CALL cp_libint_engine_init(context%deriv, 4, max_am, 1)
      ALLOCATE (context%products(128))
   END SUBROUTINE four_center_integral_init

! **************************************************************************************************
!> \brief Release the workspaces of an owning four-center context.
!> \param context original initialized owner
! **************************************************************************************************
   SUBROUTINE four_center_integral_cleanup(context)
      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context

      CALL cp_libint_engine_release(context%eri)
      CALL cp_libint_engine_release(context%deriv)
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
      REAL(dp), POINTER                                  :: engine_values(:, :)

      mysize = ncoa*ncob*ncoc*ncod

      work = 0.0_dp

      IF (use_virial) THEN
         work_virial = 0.0_dp
      END IF

      CALL get_eri_permutation(n_a, n_b, n_c, n_d, permutation)
      angular_momenta = [n_a, n_b, n_c, n_d]
      angular_momenta = angular_momenta(permutation)

      DO i = 1, context%nproducts
         CALL cp_libint_engine_evaluate(context%deriv, angular_momenta, &
                                        context%products(i)%centers(:, permutation), [1, 1, 1, 1], &
                                        RESHAPE(context%products(i)%exponents(permutation), [1, 4]), &
                                        RESHAPE([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [1, 4]), context%potential, engine_values)
         work_forces = engine_values
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
      REAL(dp)                                           :: centers(3, 4), exponents(4), PQ(3), PQ2
      REAL(dp), POINTER                                  :: engine_values(:, :), p_work(:)
      TYPE(eri4_pgf_image_type)                          :: image_a, image_b
      TYPE(eri4_pgf_pair_type)                           :: pair_a, pair_b

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

      CALL prepare_eri4_pair(centers(:, 1), centers(:, 2), exponents(1), exponents(2), &
                             pair_a, image_a)
      CALL prepare_eri4_pair(centers(:, 3), centers(:, 4), exponents(3), exponents(4), &
                             pair_b, image_b)
      PQ = image_a%P - image_b%P
      PQ2 = DOT_PRODUCT(PQ, PQ)
      SELECT CASE (potential_parameter%potential_type)
      CASE (do_potential_truncated, do_potential_short, do_potential_mix_cl_trunc)
         IF (PQ2 > (R1 + R2 + potential_parameter%cutoff_radius)**2) RETURN
      END SELECT

      CALL cp_libint_engine_evaluate(context%eri, angular_momenta, centers, [1, 1, 1, 1], &
                                     RESHAPE(exponents, [1, 4]), &
                                     RESHAPE([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [1, 4]), potential_parameter, engine_values)
      p_work => engine_values(:, 1)
      IF (m_max == 0) THEN
         max_val = ABS(p_work(1))
      ELSE
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

      INTEGER                                            :: angular_momenta(4), i, mysize, &
                                                            permutation(4)
      LOGICAL                                            :: direct_contract
      REAL(dp), POINTER                                  :: engine_values(:, :), p_work(:)

      IF (context%nproducts == 0) RETURN
      mysize = ncoa*ncob*ncoc*ncod
      CALL get_eri_permutation(n_a, n_b, n_c, n_d, permutation)
      angular_momenta = [n_a, n_b, n_c, n_d]
      angular_momenta = angular_momenta(permutation)
      IF (context%nproducts > 1) buffer1 = 0.0_dp
      DO i = 1, context%nproducts
         CALL cp_libint_engine_evaluate(context%eri, angular_momenta, &
                                        context%products(i)%centers(:, permutation), [1, 1, 1, 1], &
                                        RESHAPE(context%products(i)%exponents(permutation), [1, 4]), &
                                        RESHAPE([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [1, 4]), context%potential, engine_values)
         p_work => engine_values(:, 1)
         IF (context%nproducts > 1) buffer1 = buffer1 + p_work(1:mysize)
      END DO
      IF (context%nproducts > 1) p_work => buffer1
      tmp_max = MAX(tmp_max, MAXVAL(ABS(p_work)))
      neris = neris + INT(context%nproducts, int_8)*mysize
      tmp_max = tmp_max*max_contraction
      IF (tmp_max < eps_schwarz) RETURN
      CALL reorder_eri_block(p_work, work2, [ncoa, ncob, ncoc, ncod], permutation, 1, 1)

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

! **************************************************************************************************
!> \brief Select the canonical center order for four-center integral evaluation.
!> \param la original first-center angular momentum
!> \param lb original second-center angular momentum
!> \param lc original third-center angular momentum
!> \param ld original fourth-center angular momentum
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
!> \brief Restore original Cartesian axes and derivative center labels.
!> \param source canonical block with the fourth Cartesian index varying fastest
!> \param destination original block with the first Cartesian index varying fastest
!> \param ncart Cartesian shell sizes in original center order
!> \param permutation original center index for each canonical center
!> \param ncomponents one integral component or twelve center derivatives
!> \param nplanes number of independent planes
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

END MODULE four_center_integral_interface
