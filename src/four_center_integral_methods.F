!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Prepare screened primitive Gaussian pairs and four-center products.
! **************************************************************************************************
MODULE four_center_integral_methods
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE four_center_integral_types,      ONLY: eri4_cell_type,&
                                              eri4_fundamentals_size,&
                                              eri4_pgf_pair_type,&
                                              eri4_pgf_product_type,&
                                              eri4_screen_coeff_type
   USE input_constants,                 ONLY: do_potential_TShPSC,&
                                              do_potential_id,&
                                              do_potential_mix_cl_trunc,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE integral_library_types,          ONLY: coulomb_operator_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE r12_fundamentals,                ONLY: compute_r12_fundamentals
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: build_pair_list_pgf, build_pgf_product_list, min_image_gap_sq, pgf_product_list_size

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

      INTEGER                                            :: i, j, k, k_max, l, nimages1, nimages2, &
                                                            tmp_i4
      LOGICAL                                            :: has_r_cutoff
      REAL(dp) :: C11(3), den, Eta, EtaInv, fundamentals_work(eri4_fundamentals_size), G(3), &
         min_img_delta, normalization, num, P(3), pgf_max_1, pgf_max_2, PQ(3), pq2, prefactor, &
         Q(3), R1, R2, r_cutoff_sq, r_cutoff_sq_guard, ra(3), rb(3), rc(3), rcut_eff, rd(3), Rho, &
         RhoInv, rpq2, S1234, S1234a, S1234b, shift(3), ssss, T, temp(3), temp_CC(3), temp_DD(3), &
         tmp_D(3), W(3), Zeta1, Zeta_C, Zeta_D, ZetapEtaInv
      TYPE(eri4_pgf_product_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: tmp_product_list

      CPASSERT(m_max < eri4_fundamentals_size)

      nimages1 = list1%nimages
      nimages2 = list2%nimages
      nproducts = 0
      Zeta1 = list1%zetapzetb
      Eta = list2%zetapzetb
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
         ra = list1%image_list(i)%ra
         rb = list1%image_list(i)%rb
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

            ZetapEtaInv = Zeta1 + Eta
            ZetapEtaInv = 1.0_dp/ZetapEtaInv
            Rho = Zeta1*Eta*ZetapEtaInv
            RhoInv = 1.0_dp/Rho
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

               T = Rho*rpq2
               IF (potential_parameter%potential_type == do_potential_id) THEN
                  num = list1%zeta*list1%zetb
                  den = list1%zeta + list1%zetb
                  ssss = -num/den*SUM((ra - rb)**2)

                  num = den*Zeta_C
                  den = den + Zeta_C
                  ssss = ssss - num/den*SUM((P - rc)**2)

                  G(:) = (list1%zeta*ra(:) + list1%zetb*rb(:) + Zeta_C*rc(:))/den
                  num = den*Zeta_D
                  den = den + Zeta_D
                  ssss = ssss - num/den*SUM((G - rd)**2)

                  product_list(nproducts)%Fm(:) = EXP(ssss)
                  ! The ID seed already contains the four-Gaussian attenuation.
                  normalization = 1.0_dp
                  IF (S1234 > EPSILON(0.0_dp)) normalization = 1.0_dp/S1234
               ELSE
                  CALL compute_r12_fundamentals(product_list(nproducts)%Fm, normalization, m_max, &
                                                potential_parameter%potential_type, Rho, T, fundamentals_work, &
                                                omega=potential_parameter%omega, &
                                                cutoff_radius=potential_parameter%cutoff_radius, &
                                                scale_coulomb=potential_parameter%scale_coulomb, &
                                                scale_longrange=potential_parameter%scale_longrange, &
                                                scale_gaussian=potential_parameter%scale_gaussian)
               END IF

               prefactor = normalization*S1234*SQRT((Pi*ZetapEtaInv)**3)

               DO l = 1, m_max + 1
                  product_list(nproducts)%Fm(l) = product_list(nproducts)%Fm(l)*prefactor
               END DO

               W = (Zeta1*P + Eta*Q)*ZetapEtaInv
               product_list(nproducts)%ra = ra
               product_list(nproducts)%rb = rb
               product_list(nproducts)%rc = C11
               product_list(nproducts)%rd = tmp_D
               product_list(nproducts)%ZetapEtaInv = ZetapEtaInv
               product_list(nproducts)%Rho = Rho
               product_list(nproducts)%RhoInv = RhoInv
               product_list(nproducts)%P = P
               product_list(nproducts)%Q = Q
               product_list(nproducts)%W = W
               product_list(nproducts)%AB = ra - rb
               product_list(nproducts)%CD = C11 - tmp_D
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
                                                            rab2, Zeta1, Zeta_A, Zeta_B, ZetaInv

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
               list(element_counter)%image_list(nimages(element_counter))%ra = ra
               list(element_counter)%image_list(nimages(element_counter))%rb = im_B
               list(element_counter)%image_list(nimages(element_counter))%rab2 = rab2

               Zeta_A = zeta(ipgf)
               Zeta_B = zetb(jpgf)
               Zeta1 = Zeta_A + Zeta_B
               ZetaInv = 1.0_dp/Zeta1

               IF (nimages(element_counter) == 1) THEN
                  list(element_counter)%ipgf = ipgf
                  list(element_counter)%jpgf = jpgf
                  list(element_counter)%zetaInv = ZetaInv
                  list(element_counter)%zetapzetb = Zeta1
                  list(element_counter)%zeta = Zeta_A
                  list(element_counter)%zetb = Zeta_B
               END IF

               list(element_counter)%image_list(nimages(element_counter))%S1234 = (-Zeta_A*Zeta_B*ZetaInv*rab2)
               list(element_counter)%image_list(nimages(element_counter))%P = (Zeta_A*ra + Zeta_B*im_B)*ZetaInv
               list(element_counter)%image_list(nimages(element_counter))%R = &
                  MAX(0.0_dp, R_pgf(jpgf, ipgf)%x(1)*rab2 + R_pgf(jpgf, ipgf)%x(2))
               list(element_counter)%image_list(nimages(element_counter))%ra = ra
               list(element_counter)%image_list(nimages(element_counter))%rb = im_B
               list(element_counter)%image_list(nimages(element_counter))%rab2 = rab2
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

END MODULE four_center_integral_methods
