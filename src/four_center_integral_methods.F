!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Prepare and evaluate screened four-center Gaussian integrals and derivatives.
! **************************************************************************************************
MODULE four_center_integral_methods
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE four_center_integral_interface,  ONLY: evaluate_deriv_eri,&
                                              evaluate_eri,&
                                              four_center_integral_context_type
   USE four_center_integral_preparation,ONLY: prepare_eri4_pair,&
                                              prepare_eri4_primitive,&
                                              prepare_eri4_product
   USE four_center_integral_types,      ONLY: eri4_cell_type,&
                                              eri4_pgf_image_type,&
                                              eri4_pgf_pair_type,&
                                              eri4_pgf_product_type,&
                                              eri4_screen_coeff_type
   USE input_constants,                 ONLY: do_potential_TShPSC,&
                                              do_potential_mix_cl_trunc,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE integral_library_types,          ONLY: coulomb_operator_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE orbital_pointers,                ONLY: nco,&
                                              nso
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: coulomb4, forces4, min_image_gap_sq, pgf_product_list_size

   ! Initial capacity estimate, shared by the per-thread product buffers.
   INTEGER, SAVE :: pgf_product_list_size = 128

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param list1 ...
!> \param list2 ...
!> \param product_list ...
!> \param nproducts ...
!> \param log10_pmax ...
!> \param log10_eps_schwarz ...
!> \param neighbor_cells ...
!> \param cell ...
!> \param potential_parameter ...
!> \param m_max ...
!> \param do_periodic ...
! **************************************************************************************************
   SUBROUTINE build_pgf_product_list(list1, list2, product_list, nproducts, &
                                     log10_pmax, log10_eps_schwarz, neighbor_cells, &
                                     cell, potential_parameter, m_max, do_periodic)

      TYPE(eri4_pgf_pair_type)                           :: list1, list2
      TYPE(eri4_pgf_product_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(INOUT)                     :: product_list
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
      TYPE(eri4_pgf_product_type)                        :: product
      TYPE(eri4_pgf_product_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: tmp_product_list

      nimages1 = list1%nimages
      nimages2 = list2%nimages
      nproducts = 0
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

               ! allocate size as needed,
               ! updating the global size estimate to make this a rare event in longer simulations
               IF (nproducts > SIZE(product_list)) THEN
!$OMP              ATOMIC READ
                  tmp_i4 = pgf_product_list_size
                  tmp_i4 = MAX(pgf_product_list_size, (3*nproducts + 1)/2)
!$OMP              ATOMIC WRITE
                  pgf_product_list_size = tmp_i4
                  ALLOCATE (tmp_product_list(SIZE(product_list)))
                  tmp_product_list(:) = product_list
                  DEALLOCATE (product_list)
                  ALLOCATE (product_list(tmp_i4))
                  product_list(1:SIZE(tmp_product_list)) = tmp_product_list
                  DEALLOCATE (tmp_product_list)
               END IF

               shifted_image%ra = C11
               shifted_image%rb = tmp_D
               shifted_image%P = Q
               CALL prepare_eri4_primitive(product, list1%image_list(i), shifted_image, S1234, rpq2, &
                                           potential_parameter, product_list(nproducts))
            END DO
         END DO
      END DO

   END SUBROUTINE build_pgf_product_list

! **************************************************************************************************
!> \brief ...
!> \param npgfa ...
!> \param npgfb ...
!> \param list ...
!> \param zeta ...
!> \param zetb ...
!> \param screen1 ...
!> \param screen2 ...
!> \param pgf ...
!> \param R_pgf ...
!> \param log10_pmax ...
!> \param log10_eps_schwarz ...
!> \param ra ...
!> \param rb ...
!> \param nelements ...
!> \param neighbor_cells ...
!> \param nimages ...
!> \param do_periodic ...
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE build_pair_list_pgf(npgfa, npgfb, list, zeta, zetb, screen1, screen2, pgf, R_pgf, &
                                  log10_pmax, log10_eps_schwarz, ra, rb, nelements, &
                                  neighbor_cells, nimages, do_periodic, cell)
      INTEGER, INTENT(IN)                                :: npgfa, npgfb
      TYPE(eri4_pgf_pair_type), DIMENSION(npgfa*npgfb)   :: list
      REAL(dp), DIMENSION(1:npgfa), INTENT(IN)           :: zeta
      REAL(dp), DIMENSION(1:npgfb), INTENT(IN)           :: zetb
      REAL(dp), INTENT(IN)                               :: screen1(2), screen2(2)
      TYPE(eri4_screen_coeff_type), DIMENSION(:, :), &
         POINTER                                         :: pgf, R_pgf
      REAL(dp), INTENT(IN)                               :: log10_pmax, log10_eps_schwarz, ra(3), &
                                                            rb(3)
      INTEGER, INTENT(OUT)                               :: nelements
      TYPE(eri4_cell_type), DIMENSION(:), POINTER        :: neighbor_cells
      INTEGER                                            :: nimages(npgfa*npgfb)
      LOGICAL, INTENT(IN)                                :: do_periodic
      TYPE(cell_type), POINTER                           :: cell

      REAL(dp), PARAMETER                                :: min_img_eps = 1.0E-12_dp

      INTEGER                                            :: element_counter, i, i_max, ipgf, j, jpgf
      REAL(dp)                                           :: AB(3), im_B(3), min_img_delta, pgf_max, &
                                                            rab2

      i_max = SIZE(neighbor_cells)
      IF (do_periodic) THEN
         IF (cell%orthorhombic .AND. screen1(1) <= 0.0_dp .AND. i_max > 1) THEN
            IF (ALL(cell%perd == 1) .AND. ALL(neighbor_cells(1)%cell_r == 0.0_dp)) THEN
               AB = ra - rb
               rab2 = AB(1)**2 + AB(2)**2 + AB(3)**2
               min_img_delta = min_image_gap_sq(AB, cell%hmat)
               IF (min_img_delta >= 0.0_dp) THEN
                  IF (screen1(1)*rab2 + screen1(2) + screen2(2) + log10_pmax < log10_eps_schwarz) THEN
                     nelements = 0
                     RETURN
                  END IF
                  IF (screen1(1)*(rab2 + min_img_delta)*(1.0_dp - min_img_eps) + &
                      screen1(2) + screen2(2) + log10_pmax < log10_eps_schwarz) i_max = 1
               END IF
            END IF
         END IF
      ELSE
         i_max = MIN(i_max, 1)
      END IF
      nimages = 0
      ! ** inner loop may never be reached
      nelements = npgfa*npgfb
      DO i = 1, i_max
         IF (do_periodic) THEN
            im_B = rb + neighbor_cells(i)%cell_r(:)
         ELSE
            im_B = rb
         END IF
         AB = ra - im_B
         rab2 = AB(1)**2 + AB(2)**2 + AB(3)**2
         IF (screen1(1)*rab2 + screen1(2) + screen2(2) + log10_pmax < log10_eps_schwarz) CYCLE
         element_counter = 0
         DO ipgf = 1, npgfa
            DO jpgf = 1, npgfb
               element_counter = element_counter + 1
               pgf_max = pgf(jpgf, ipgf)%x(1)*rab2 + pgf(jpgf, ipgf)%x(2)
               IF (pgf_max + screen2(2) + log10_pmax < log10_eps_schwarz) THEN
                  CYCLE
               END IF
               nimages(element_counter) = nimages(element_counter) + 1
               list(element_counter)%image_list(nimages(element_counter))%pgf_max = pgf_max
               CALL prepare_eri4_pair(ra, im_B, zeta(ipgf), zetb(jpgf), rab2, &
                                      list(element_counter), list(element_counter)%image_list(nimages(element_counter)))
               IF (nimages(element_counter) == 1) THEN
                  list(element_counter)%ipgf = ipgf
                  list(element_counter)%jpgf = jpgf
               END IF

               list(element_counter)%image_list(nimages(element_counter))%R = &
                  MAX(0.0_dp, R_pgf(jpgf, ipgf)%x(1)*rab2 + R_pgf(jpgf, ipgf)%x(2))
               list(element_counter)%image_list(nimages(element_counter))%bcell = neighbor_cells(i)%cell
            END DO
         END DO
         nelements = MAX(nelements, element_counter)
      END DO
      DO j = 1, nelements
         list(j)%nimages = nimages(j)
      END DO
      ! ** Remove unused elements

      element_counter = 0
      DO j = 1, nelements
         IF (list(j)%nimages == 0) CYCLE
         element_counter = element_counter + 1
         list(element_counter)%nimages = list(j)%nimages
         list(element_counter)%zetapzetb = list(j)%zetapzetb
         list(element_counter)%ZetaInv = list(j)%ZetaInv
         list(element_counter)%zeta = list(j)%zeta
         list(element_counter)%zetb = list(j)%zetb
         list(element_counter)%ipgf = list(j)%ipgf
         list(element_counter)%jpgf = list(j)%jpgf
         DO i = 1, list(j)%nimages
            list(element_counter)%image_list(i) = list(j)%image_list(i)
         END DO
      END DO

      nelements = element_counter

   END SUBROUTINE build_pair_list_pgf

! **************************************************************************************************
!> \brief Smallest increase of |r|**2 when r is moved to any other periodic image of an
!>        orthorhombic, fully periodic cell, i.e. min over n /= 0 of |r - hmat*n|**2 - |r|**2.
!>        Exact if r is a minimum-image vector (|r(i)| <= |hmat(i,i)|/2); negative otherwise.
!> \param r displacement vector
!> \param hmat cell matrix; its diagonal entries are signed cell-vector components
!> \return the squared-distance gap to the nearest other periodic image
! **************************************************************************************************
   PURE FUNCTION min_image_gap_sq(r, hmat) RESULT(gap)
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: r
      REAL(dp), DIMENSION(3, 3), INTENT(IN)              :: hmat
      REAL(dp)                                           :: gap

      REAL(dp), DIMENSION(3)                             :: cell_len

      cell_len = [ABS(hmat(1, 1)), ABS(hmat(2, 2)), ABS(hmat(3, 3))]
      gap = MINVAL(cell_len*(cell_len - 2.0_dp*ABS(r)))

   END FUNCTION min_image_gap_sq

! **************************************************************************************************
!> \brief calculates two-electron integrals of a quartet/shell using the library
!>      lib_int in the periodic case
!> \param context four-center integral context
!> \param ra ...
!> \param rb ...
!> \param rc ...
!> \param rd ...
!> \param npgfa ...
!> \param npgfb ...
!> \param npgfc ...
!> \param npgfd ...
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param lc_min ...
!> \param lc_max ...
!> \param ld_min ...
!> \param ld_max ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param nsgfc ...
!> \param nsgfd ...
!> \param sphi_a_u1 ...
!> \param sphi_a_u2 ...
!> \param sphi_a_u3 ...
!> \param sphi_b_u1 ...
!> \param sphi_b_u2 ...
!> \param sphi_b_u3 ...
!> \param sphi_c_u1 ...
!> \param sphi_c_u2 ...
!> \param sphi_c_u3 ...
!> \param sphi_d_u1 ...
!> \param sphi_d_u2 ...
!> \param sphi_d_u3 ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param zetd ...
!> \param primitive_integrals array of primitive_integrals
!> \param potential_parameter contains info for libint
!> \param neighbor_cells Periodic images
!> \param screen1 set based coefficients for near field screening
!> \param screen2 set based coefficients for near field screening
!> \param eps_schwarz threshold
!> \param max_contraction_val maximum multiplication factor for cart -> sph
!> \param cart_estimate maximum calculated integral value
!> \param cell cell
!> \param neris_tmp counter for calculated cart integrals
!> \param log10_pmax logarithm of initial p matrix max element
!> \param log10_eps_schwarz log of threshold
!> \param R1_pgf coefficients for radii of product distribution function
!> \param R2_pgf coefficients for radii of product distribution function
!> \param pgf1 schwarz coefficients pgf basid
!> \param pgf2 schwarz coefficients pgf basid
!> \param pgf_list_ij ...
!> \param pgf_list_kl ...
!> \param pgf_product_list ...
!> \param nsgfl_a ...
!> \param nsgfl_b ...
!> \param nsgfl_c ...
!> \param nsgfl_d ...
!> \param sphi_a_ext ...
!> \param sphi_b_ext ...
!> \param sphi_c_ext ...
!> \param sphi_d_ext ...
!> \param ee_work ...
!> \param ee_work2 ...
!> \param ee_buffer1 ...
!> \param ee_buffer2 ...
!> \param ee_primitives_tmp ...
!> \param nimages ...
!> \param do_periodic ...
!> \par History
!>      11.2006 created [Manuel Guidon]
!>      02.2009 completely rewritten screening part [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE coulomb4(context, ra, rb, rc, rd, npgfa, npgfb, npgfc, npgfd, &
                       la_min, la_max, lb_min, lb_max, &
                       lc_min, lc_max, ld_min, ld_max, nsgfa, nsgfb, nsgfc, nsgfd, &
                       sphi_a_u1, sphi_a_u2, sphi_a_u3, &
                       sphi_b_u1, sphi_b_u2, sphi_b_u3, &
                       sphi_c_u1, sphi_c_u2, sphi_c_u3, &
                       sphi_d_u1, sphi_d_u2, sphi_d_u3, &
                       zeta, zetb, zetc, zetd, &
                       primitive_integrals, &
                       potential_parameter, neighbor_cells, &
                       screen1, screen2, eps_schwarz, max_contraction_val, &
                       cart_estimate, cell, neris_tmp, log10_pmax, &
                       log10_eps_schwarz, R1_pgf, R2_pgf, pgf1, pgf2, &
                       pgf_list_ij, pgf_list_kl, &
                       pgf_product_list, &
                       nsgfl_a, nsgfl_b, nsgfl_c, &
                       nsgfl_d, &
                       sphi_a_ext, sphi_b_ext, sphi_c_ext, sphi_d_ext, &
                       ee_work, ee_work2, ee_buffer1, ee_buffer2, ee_primitives_tmp, &
                       nimages, do_periodic)

      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context
      REAL(dp), INTENT(IN)                               :: ra(3), rb(3), rc(3), rd(3)
      INTEGER, INTENT(IN) :: npgfa, npgfb, npgfc, npgfd, la_min, la_max, lb_min, lb_max, lc_min, &
         lc_max, ld_min, ld_max, nsgfa, nsgfb, nsgfc, nsgfd, sphi_a_u1, sphi_a_u2, sphi_a_u3, &
         sphi_b_u1, sphi_b_u2, sphi_b_u3, sphi_c_u1, sphi_c_u2, sphi_c_u3, sphi_d_u1, sphi_d_u2, &
         sphi_d_u3
      REAL(dp), DIMENSION(1:npgfa), INTENT(IN)           :: zeta
      REAL(dp), DIMENSION(1:npgfb), INTENT(IN)           :: zetb
      REAL(dp), DIMENSION(1:npgfc), INTENT(IN)           :: zetc
      REAL(dp), DIMENSION(1:npgfd), INTENT(IN)           :: zetd
      REAL(dp), DIMENSION(nsgfa, nsgfb, nsgfc, nsgfd)    :: primitive_integrals
      TYPE(coulomb_operator_type)                        :: potential_parameter
      TYPE(eri4_cell_type), DIMENSION(:), POINTER        :: neighbor_cells
      REAL(dp), INTENT(IN)                               :: screen1(2), screen2(2), eps_schwarz, &
                                                            max_contraction_val
      REAL(dp)                                           :: cart_estimate
      TYPE(cell_type), POINTER                           :: cell
      INTEGER(int_8)                                     :: neris_tmp
      REAL(dp), INTENT(IN)                               :: log10_pmax, log10_eps_schwarz
      TYPE(eri4_screen_coeff_type), DIMENSION(:, :), &
         POINTER                                         :: R1_pgf, R2_pgf, pgf1, pgf2
      TYPE(eri4_pgf_pair_type), DIMENSION(*)             :: pgf_list_ij, pgf_list_kl
      TYPE(eri4_pgf_product_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(INOUT)                     :: pgf_product_list
      INTEGER, DIMENSION(0:), INTENT(IN)                 :: nsgfl_a, nsgfl_b, nsgfl_c, nsgfl_d
      REAL(dp), INTENT(IN) :: sphi_a_ext(sphi_a_u1, sphi_a_u2, sphi_a_u3), &
         sphi_b_ext(sphi_b_u1, sphi_b_u2, sphi_b_u3), sphi_c_ext(sphi_c_u1, sphi_c_u2, sphi_c_u3), &
         sphi_d_ext(sphi_d_u1, sphi_d_u2, sphi_d_u3)
      REAL(dp), DIMENSION(*)                             :: ee_work, ee_work2, ee_buffer1, &
                                                            ee_buffer2, ee_primitives_tmp
      INTEGER, DIMENSION(*)                              :: nimages
      LOGICAL, INTENT(IN)                                :: do_periodic

      INTEGER :: ipgf, jpgf, kpgf, la, lb, lc, ld, list_ij, list_kl, lpgf, max_l, ncoa, ncob, &
         ncoc, ncod, nelements_ij, nelements_kl, nproducts, nsgfla, nsgflb, nsgflc, nsgfld, nsoa, &
         nsob, nsoc, nsod, s_offset_a, s_offset_b, s_offset_c, s_offset_d
      REAL(dp)                                           :: tmp_max

      cart_estimate = 0.0_dp
      neris_tmp = 0
      primitive_integrals = 0.0_dp

      CALL build_pair_list_pgf(npgfa, npgfb, pgf_list_ij, zeta, zetb, screen1, screen2, &
                               pgf1, R1_pgf, log10_pmax, log10_eps_schwarz, ra, rb, &
                               nelements_ij, &
                               neighbor_cells, nimages, do_periodic, cell)
      IF (nelements_ij == 0) RETURN
      CALL build_pair_list_pgf(npgfc, npgfd, pgf_list_kl, zetc, zetd, screen2, screen1, &
                               pgf2, R2_pgf, log10_pmax, log10_eps_schwarz, rc, rd, &
                               nelements_kl, &
                               neighbor_cells, nimages, do_periodic, cell)
      IF (nelements_kl == 0) RETURN

      max_l = la_max + lb_max + lc_max + ld_max

      DO list_ij = 1, nelements_ij
         ipgf = pgf_list_ij(list_ij)%ipgf
         jpgf = pgf_list_ij(list_ij)%jpgf

         DO list_kl = 1, nelements_kl
            kpgf = pgf_list_kl(list_kl)%ipgf
            lpgf = pgf_list_kl(list_kl)%jpgf

            CALL build_pgf_product_list(pgf_list_ij(list_ij), pgf_list_kl(list_kl), pgf_product_list, &
                                        nproducts, log10_pmax, log10_eps_schwarz, neighbor_cells, cell, &
                                        potential_parameter, max_l, do_periodic)
            IF (nproducts == 0) CYCLE

            s_offset_a = 0
            DO la = la_min, la_max
               s_offset_b = 0
               ncoa = nco(la)
               nsgfla = nsgfl_a(la)
               nsoa = nso(la)

               DO lb = lb_min, lb_max
                  s_offset_c = 0
                  ncob = nco(lb)
                  nsgflb = nsgfl_b(lb)
                  nsob = nso(lb)

                  DO lc = lc_min, lc_max
                     s_offset_d = 0
                     ncoc = nco(lc)
                     nsgflc = nsgfl_c(lc)
                     nsoc = nso(lc)

                     DO ld = ld_min, ld_max
                        ncod = nco(ld)
                        nsgfld = nsgfl_d(ld)
                        nsod = nso(ld)

                        tmp_max = 0.0_dp
                        CALL evaluate_eri(context, nproducts, pgf_product_list, &
                                          la, lb, lc, ld, &
                                          ncoa, ncob, ncoc, ncod, &
                                          nsgfa, nsgfb, nsgfc, nsgfd, &
                                          primitive_integrals, &
                                          max_contraction_val, tmp_max, eps_schwarz, &
                                          neris_tmp, &
                                          s_offset_a, s_offset_b, s_offset_c, s_offset_d, &
                                          nsgfla, nsgflb, nsgflc, nsgfld, nsoa, nsob, nsoc, nsod, &
                                          sphi_a_ext(1, la + 1, ipgf), &
                                          sphi_b_ext(1, lb + 1, jpgf), &
                                          sphi_c_ext(1, lc + 1, kpgf), &
                                          sphi_d_ext(1, ld + 1, lpgf), &
                                          ee_work, ee_work2, ee_buffer1, ee_buffer2, ee_primitives_tmp)
                        cart_estimate = MAX(tmp_max, cart_estimate)
                        s_offset_d = s_offset_d + nsod*nsgfld
                     END DO !ld
                     s_offset_c = s_offset_c + nsoc*nsgflc
                  END DO !lc
                  s_offset_b = s_offset_b + nsob*nsgflb
               END DO !lb
               s_offset_a = s_offset_a + nsoa*nsgfla
            END DO !la
         END DO
      END DO

   END SUBROUTINE coulomb4

! **************************************************************************************************
!> \brief ...
!> \param context four-center integral context
!> \param ra ...
!> \param rb ...
!> \param rc ...
!> \param rd ...
!> \param npgfa ...
!> \param npgfb ...
!> \param npgfc ...
!> \param npgfd ...
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param lc_min ...
!> \param lc_max ...
!> \param ld_min ...
!> \param ld_max ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param nsgfc ...
!> \param nsgfd ...
!> \param sphi_a_u1 ...
!> \param sphi_a_u2 ...
!> \param sphi_a_u3 ...
!> \param sphi_b_u1 ...
!> \param sphi_b_u2 ...
!> \param sphi_b_u3 ...
!> \param sphi_c_u1 ...
!> \param sphi_c_u2 ...
!> \param sphi_c_u3 ...
!> \param sphi_d_u1 ...
!> \param sphi_d_u2 ...
!> \param sphi_d_u3 ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param zetd ...
!> \param primitive_integrals ...
!> \param potential_parameter ...
!> \param neighbor_cells ...
!> \param screen1 ...
!> \param screen2 ...
!> \param eps_schwarz ...
!> \param max_contraction_val ...
!> \param cart_estimate ...
!> \param cell ...
!> \param neris_tmp ...
!> \param log10_pmax ...
!> \param log10_eps_schwarz ...
!> \param R1_pgf ...
!> \param R2_pgf ...
!> \param pgf1 ...
!> \param pgf2 ...
!> \param pgf_list_ij ...
!> \param pgf_list_kl ...
!> \param pgf_product_list ...
!> \param nsgfl_a ...
!> \param nsgfl_b ...
!> \param nsgfl_c ...
!> \param nsgfl_d ...
!> \param sphi_a_ext ...
!> \param sphi_b_ext ...
!> \param sphi_c_ext ...
!> \param sphi_d_ext ...
!> \param ede_work ...
!> \param ede_work2 ...
!> \param ede_work_forces ...
!> \param ede_buffer1 ...
!> \param ede_buffer2 ...
!> \param ede_primitives_tmp ...
!> \param nimages ...
!> \param do_periodic ...
!> \param use_virial ...
!> \param ede_work_virial ...
!> \param ede_work2_virial ...
!> \param ede_primitives_tmp_virial ...
!> \param primitive_integrals_virial ...
!> \param cart_estimate_virial ...
! **************************************************************************************************
   SUBROUTINE forces4(context, ra, rb, rc, rd, npgfa, npgfb, npgfc, npgfd, &
                      la_min, la_max, lb_min, lb_max, &
                      lc_min, lc_max, ld_min, ld_max, nsgfa, nsgfb, nsgfc, nsgfd, &
                      sphi_a_u1, sphi_a_u2, sphi_a_u3, &
                      sphi_b_u1, sphi_b_u2, sphi_b_u3, &
                      sphi_c_u1, sphi_c_u2, sphi_c_u3, &
                      sphi_d_u1, sphi_d_u2, sphi_d_u3, &
                      zeta, zetb, zetc, zetd, &
                      primitive_integrals, &
                      potential_parameter, neighbor_cells, &
                      screen1, screen2, eps_schwarz, max_contraction_val, &
                      cart_estimate, cell, neris_tmp, log10_pmax, &
                      log10_eps_schwarz, R1_pgf, R2_pgf, pgf1, pgf2, &
                      pgf_list_ij, pgf_list_kl, &
                      pgf_product_list, &
                      nsgfl_a, nsgfl_b, nsgfl_c, &
                      nsgfl_d, &
                      sphi_a_ext, sphi_b_ext, sphi_c_ext, sphi_d_ext, &
                      ede_work, ede_work2, ede_work_forces, &
                      ede_buffer1, ede_buffer2, ede_primitives_tmp, &
                      nimages, do_periodic, use_virial, ede_work_virial, &
                      ede_work2_virial, ede_primitives_tmp_virial, primitive_integrals_virial, &
                      cart_estimate_virial)

      TYPE(four_center_integral_context_type), &
         INTENT(INOUT)                                   :: context
      REAL(dp), INTENT(IN)                               :: ra(3), rb(3), rc(3), rd(3)
      INTEGER, INTENT(IN) :: npgfa, npgfb, npgfc, npgfd, la_min, la_max, lb_min, lb_max, lc_min, &
         lc_max, ld_min, ld_max, nsgfa, nsgfb, nsgfc, nsgfd, sphi_a_u1, sphi_a_u2, sphi_a_u3, &
         sphi_b_u1, sphi_b_u2, sphi_b_u3, sphi_c_u1, sphi_c_u2, sphi_c_u3, sphi_d_u1, sphi_d_u2, &
         sphi_d_u3
      REAL(dp), DIMENSION(1:npgfa), INTENT(IN)           :: zeta
      REAL(dp), DIMENSION(1:npgfb), INTENT(IN)           :: zetb
      REAL(dp), DIMENSION(1:npgfc), INTENT(IN)           :: zetc
      REAL(dp), DIMENSION(1:npgfd), INTENT(IN)           :: zetd
      REAL(dp), &
         DIMENSION(nsgfa, nsgfb, nsgfc, nsgfd, 12)       :: primitive_integrals
      TYPE(coulomb_operator_type)                        :: potential_parameter
      TYPE(eri4_cell_type), DIMENSION(:), POINTER        :: neighbor_cells
      REAL(dp), INTENT(IN)                               :: screen1(2), screen2(2), eps_schwarz, &
                                                            max_contraction_val
      REAL(dp)                                           :: cart_estimate
      TYPE(cell_type), POINTER                           :: cell
      INTEGER(int_8)                                     :: neris_tmp
      REAL(dp), INTENT(IN)                               :: log10_pmax, log10_eps_schwarz
      TYPE(eri4_screen_coeff_type), DIMENSION(:, :), &
         POINTER                                         :: R1_pgf, R2_pgf, pgf1, pgf2
      TYPE(eri4_pgf_pair_type), DIMENSION(*)             :: pgf_list_ij, pgf_list_kl
      TYPE(eri4_pgf_product_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(INOUT)                     :: pgf_product_list
      INTEGER, DIMENSION(0:), INTENT(IN)                 :: nsgfl_a, nsgfl_b, nsgfl_c, nsgfl_d
      REAL(dp), INTENT(IN) :: sphi_a_ext(sphi_a_u1, sphi_a_u2, sphi_a_u3), &
         sphi_b_ext(sphi_b_u1, sphi_b_u2, sphi_b_u3), sphi_c_ext(sphi_c_u1, sphi_c_u2, sphi_c_u3), &
         sphi_d_ext(sphi_d_u1, sphi_d_u2, sphi_d_u3)
      REAL(dp), DIMENSION(*)                             :: ede_work, ede_work2, ede_work_forces, &
                                                            ede_buffer1, ede_buffer2, &
                                                            ede_primitives_tmp
      INTEGER, DIMENSION(*)                              :: nimages
      LOGICAL, INTENT(IN)                                :: do_periodic, use_virial
      REAL(dp), DIMENSION(*)                             :: ede_work_virial, ede_work2_virial, &
                                                            ede_primitives_tmp_virial
      REAL(dp), &
         DIMENSION(nsgfa, nsgfb, nsgfc, nsgfd, 12, 3)    :: primitive_integrals_virial
      REAL(dp)                                           :: cart_estimate_virial

      INTEGER :: ipgf, jpgf, kpgf, la, lb, lc, ld, list_ij, list_kl, lpgf, max_l, ncoa, ncob, &
         ncoc, ncod, nelements_ij, nelements_kl, nproducts, nsgfla, nsgflb, nsgflc, nsgfld, nsoa, &
         nsob, nsoc, nsod, s_offset_a, s_offset_b, s_offset_c, s_offset_d
      REAL(dp)                                           :: tmp_max, tmp_max_virial

      cart_estimate = 0.0_dp
      primitive_integrals = 0.0_dp
      IF (use_virial) THEN
         primitive_integrals_virial = 0.0_dp
         cart_estimate_virial = 0.0_dp
      END IF
      neris_tmp = 0_int_8

      CALL build_pair_list_pgf(npgfa, npgfb, pgf_list_ij, zeta, zetb, screen1, screen2, &
                               pgf1, R1_pgf, log10_pmax, log10_eps_schwarz, ra, rb, &
                               nelements_ij, &
                               neighbor_cells, nimages, do_periodic, cell)
      IF (nelements_ij == 0) RETURN
      CALL build_pair_list_pgf(npgfc, npgfd, pgf_list_kl, zetc, zetd, screen2, screen1, &
                               pgf2, R2_pgf, log10_pmax, log10_eps_schwarz, rc, rd, &
                               nelements_kl, &
                               neighbor_cells, nimages, do_periodic, cell)
      IF (nelements_kl == 0) RETURN

      max_l = la_max + lb_max + lc_max + ld_max + 1
      DO list_ij = 1, nelements_ij
         ipgf = pgf_list_ij(list_ij)%ipgf
         jpgf = pgf_list_ij(list_ij)%jpgf

         DO list_kl = 1, nelements_kl
            kpgf = pgf_list_kl(list_kl)%ipgf
            lpgf = pgf_list_kl(list_kl)%jpgf

            CALL build_pgf_product_list(pgf_list_ij(list_ij), pgf_list_kl(list_kl), pgf_product_list, &
                                        nproducts, log10_pmax, log10_eps_schwarz, neighbor_cells, cell, &
                                        potential_parameter, max_l, do_periodic)
            IF (nproducts == 0) CYCLE
            s_offset_a = 0
            DO la = la_min, la_max
               s_offset_b = 0
               ncoa = nco(la)
               nsgfla = nsgfl_a(la)
               nsoa = nso(la)

               DO lb = lb_min, lb_max
                  s_offset_c = 0
                  ncob = nco(lb)
                  nsgflb = nsgfl_b(lb)
                  nsob = nso(lb)

                  DO lc = lc_min, lc_max
                     s_offset_d = 0
                     ncoc = nco(lc)
                     nsgflc = nsgfl_c(lc)
                     nsoc = nso(lc)

                     DO ld = ld_min, ld_max
                        ncod = nco(ld)
                        nsgfld = nsgfl_d(ld)
                        nsod = nso(ld)
                        tmp_max = 0.0_dp
                        tmp_max_virial = 0.0_dp
                        CALL evaluate_deriv_eri(context, nproducts, pgf_product_list, &
                                                la, lb, lc, ld, &
                                                ncoa, ncob, ncoc, ncod, &
                                                nsgfa, nsgfb, nsgfc, nsgfd, &
                                                primitive_integrals, &
                                                max_contraction_val, tmp_max, eps_schwarz, neris_tmp, &
                                                s_offset_a, s_offset_b, s_offset_c, s_offset_d, &
                                                nsgfla, nsgflb, nsgflc, nsgfld, nsoa, nsob, nsoc, nsod, &
                                                sphi_a_ext(1, la + 1, ipgf), &
                                                sphi_b_ext(1, lb + 1, jpgf), &
                                                sphi_c_ext(1, lc + 1, kpgf), &
                                                sphi_d_ext(1, ld + 1, lpgf), &
                                                ede_work, ede_work2, ede_work_forces, &
                                                ede_buffer1, ede_buffer2, ede_primitives_tmp, use_virial, &
                                                ede_work_virial, ede_work2_virial, ede_primitives_tmp_virial, &
                                                primitive_integrals_virial, cell, tmp_max_virial)
                        cart_estimate = MAX(tmp_max, cart_estimate)
                        IF (use_virial) cart_estimate_virial = MAX(tmp_max_virial, cart_estimate_virial)
                        s_offset_d = s_offset_d + nsod*nsgfld
                     END DO !ld
                     s_offset_c = s_offset_c + nsoc*nsgflc
                  END DO !lc
                  s_offset_b = s_offset_b + nsob*nsgflb
               END DO !lb
               s_offset_a = s_offset_a + nsoa*nsgfla
            END DO !la
         END DO
      END DO

   END SUBROUTINE forces4

END MODULE four_center_integral_methods
