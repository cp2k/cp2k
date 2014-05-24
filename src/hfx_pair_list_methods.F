!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines for optimizing load balance between processes in HFX calculations
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE hfx_pair_list_methods
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE gamma,                           ONLY: fgamma=>fgamma_0
  USE hfx_libint_wrapper_types,        ONLY: prim_data_f_size
  USE hfx_types,                       ONLY: &
       hfx_basis_type, hfx_block_range_type, hfx_cell_type, hfx_pgf_list, &
       hfx_pgf_product_list, hfx_potential_type, hfx_screen_coeff_type, &
       pair_list_type, pair_set_list_type
  USE input_constants,                 ONLY: &
       do_hfx_potential_coulomb, do_hfx_potential_gaussian, &
       do_hfx_potential_id, do_hfx_potential_long, do_hfx_potential_mix_cl, &
       do_hfx_potential_mix_cl_trunc, do_hfx_potential_mix_lg, &
       do_hfx_potential_short, do_hfx_potential_truncated, &
       do_mp2_potential_TShPSC
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE mp2_types,                       ONLY: pair_list_type_mp2
  USE particle_types,                  ONLY: particle_type
  USE t_c_g0,                          ONLY: t_c_g0_n
  USE t_sh_p_s_c,                      ONLY: trunc_CS_poly_n20
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC build_pair_list, &
         build_pair_list_mp2, &
         build_pair_list_pgf, &
         build_pgf_product_list, &
         build_atomic_pair_list

!***

  CONTAINS

  SUBROUTINE build_pgf_product_list(nneighbors3, list1, list2, product_list, nproducts, &
                                    log10_pmax, log10_eps_schwarz, neighbor_cells,&
                                    cell, potential_parameter, m_max, do_periodic)

    INTEGER                                  :: nneighbors3
    TYPE(hfx_pgf_list)                       :: list1, list2
    TYPE(hfx_pgf_product_list), &
      DIMENSION(nneighbors3), INTENT(OUT)    :: product_list
    INTEGER, INTENT(OUT)                     :: nproducts
    REAL(dp), INTENT(IN)                     :: log10_pmax, log10_eps_schwarz
    TYPE(hfx_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    TYPE(cell_type), POINTER                 :: cell
    TYPE(hfx_potential_type)                 :: potential_parameter
    INTEGER, INTENT(IN)                      :: m_max
    LOGICAL, INTENT(IN)                      :: do_periodic

    INTEGER                                  :: i, j, k, l, nimages1, nimages2
    LOGICAL                                  :: use_gamma
    REAL(dp) :: C11(3), Eta, EtaInv, factor, Fm(prim_data_f_size), omega2, &
      omega_corr, omega_corr2, P(3), pgf_max_1, pgf_max_2, PQ(3), Q(3), R, &
      R1, R2, ra(3), rb(3), rc(3), rd(3), Rho, RhoInv, rpq2, S1234, S1234a, &
      S1234b, shift(3), T, temp(3), temp_CC(3), temp_DD(3), tmp, tmp_D(3), &
      W(3), Zeta1, Zeta_C, Zeta_D, ZetapEtaInv

    nimages1 = list1%nimages
    nimages2 = list2%nimages
    nproducts = 0
    Zeta1 = list1%zetapzetb
    Eta = list2%zetapzetb
    EtaInv = list2%ZetaInv
    Zeta_C = list2%zeta
    Zeta_D = list2%zetb
    DO i = 1,nimages1
      P = list1%image_list(i)%P
      R1 = list1%image_list(i)%R
      S1234a = list1%image_list(i)%S1234
      pgf_max_1 = list1%image_list(i)%pgf_max
      ra = list1%image_list(i)%ra
      rb = list1%image_list(i)%rb
      DO j = 1,nimages2
        pgf_max_2 = list2%image_list(j)%pgf_max
        IF(pgf_max_1 + pgf_max_2 + log10_pmax < log10_eps_schwarz ) CYCLE
        Q = list2%image_list(j)%P
        R2 = list2%image_list(j)%R
        S1234b = list2%image_list(j)%S1234
        rc = list2%image_list(j)%ra
        rd = list2%image_list(j)%rb

        ZetapEtaInv = Zeta1+Eta
        ZetapEtaInv = 1.0_dp/ZetapEtaInv
        Rho = Zeta1*Eta*ZetapEtaInv
        RhoInv = 1.0_dp/Rho
        S1234 = EXP(S1234a+S1234b)
        IF( do_periodic ) THEN
          temp = P-Q
          PQ = pbc(temp,cell)
          shift = -PQ + temp
          temp_CC = rc+shift
          temp_DD = rd+shift
        END IF

        DO k = 1,SIZE(neighbor_cells)
          IF( do_periodic )  THEN
            C11 = temp_CC + neighbor_cells(k)%cell_r(:)
            tmp_D = temp_DD + neighbor_cells(k)%cell_r(:)
          ELSE
            C11 = rc
            tmp_D = rd
          END IF
          Q = (Zeta_C*C11 + Zeta_D*tmp_D)*EtaInv
          rpq2 = (P(1)-Q(1))**2+(P(2)-Q(2))**2+(P(3)-Q(3))**2
          IF( potential_parameter%potential_type == do_hfx_potential_truncated .OR. &
              potential_parameter%potential_type == do_hfx_potential_short .OR. &
              potential_parameter%potential_type == do_hfx_potential_mix_cl_trunc) THEN
            IF( rpq2 > (R1 + R2 + potential_parameter%cutoff_radius)**2) CYCLE
          END IF
          IF(potential_parameter%potential_type == do_mp2_potential_TShPSC) THEN
            IF( rpq2 > (R1 + R2 + potential_parameter%cutoff_radius*2.0_dp)**2) CYCLE
          END IF
          nproducts = nproducts + 1

          T = Rho * rpq2
          SELECT CASE(potential_parameter%potential_type)
          CASE( do_hfx_potential_truncated )
            R = potential_parameter%cutoff_radius * SQRT(Rho)
            CALL t_c_g0_n(product_list(nproducts)%Fm(1),use_gamma,R,T,m_max)
            IF( use_gamma ) CALL fgamma(m_max, T, product_list(nproducts)%Fm(1))
            factor = 2.0_dp * Pi * RhoInv
          CASE( do_mp2_potential_TShPSC )
            R = potential_parameter%cutoff_radius * SQRT(Rho)
            product_list(nproducts)%Fm=0.0_dp
            CALL trunc_CS_poly_n20(product_list(nproducts)%Fm(1),R,T,m_max)
            factor = 2.0_dp * Pi * RhoInv
          CASE(do_hfx_potential_coulomb)
            CALL fgamma(m_max,T,product_list(nproducts)%Fm(1))
            factor = 2.0_dp*Pi*RhoInv
          CASE(do_hfx_potential_short)
            CALL fgamma(m_max,T,product_list(nproducts)%Fm(1))
            omega2 = potential_parameter%omega**2
            omega_corr2 = omega2/(omega2+Rho)
            omega_corr = SQRT(omega_corr2)
            T = T*omega_corr2
            CALL fgamma(m_max,T,Fm)
            tmp = - omega_corr
            DO l=1,m_max+1
              product_list(nproducts)%Fm(l)=product_list(nproducts)%Fm(l) + Fm(l)*tmp
              tmp = tmp * omega_corr2
            END DO
            factor = 2.0_dp*Pi*RhoInv
          CASE(do_hfx_potential_long)
            omega2 = potential_parameter%omega**2
            omega_corr2 = omega2/(omega2+Rho)
            omega_corr = SQRT(omega_corr2)
            T = T*omega_corr2
            CALL fgamma(m_max,T,product_list(nproducts)%Fm(1))
            tmp = omega_corr
            DO l=1,m_max+1
              product_list(nproducts)%Fm(l)= product_list(nproducts)%Fm(l)*tmp
              tmp = tmp * omega_corr2
            END DO
            factor = 2.0_dp*Pi*RhoInv
          CASE(do_hfx_potential_mix_cl)
            CALL fgamma(m_max,T,product_list(nproducts)%Fm(1))
            omega2 = potential_parameter%omega**2
            omega_corr2 = omega2/(omega2+Rho)
            omega_corr = SQRT(omega_corr2)
            T = T*omega_corr2
            CALL fgamma(m_max,T,Fm)
            tmp = omega_corr
            DO l=1,m_max+1
              product_list(nproducts)%Fm(l) = &
                       product_list(nproducts)%Fm(l)*potential_parameter%scale_coulomb &
                       + Fm(l)*tmp*potential_parameter%scale_longrange
              tmp = tmp * omega_corr2
            END DO
            factor = 2.0_dp*Pi*RhoInv
          CASE(do_hfx_potential_mix_cl_trunc)

            ! truncated
            R = potential_parameter%cutoff_radius * SQRT(rho)
            CALL t_c_g0_n(product_list(nproducts)%Fm(1),use_gamma,R,T,m_max)
            IF( use_gamma ) CALL fgamma(m_max, T, product_list(nproducts)%Fm(1))

            ! Coulomb
            CALL fgamma(m_max,T,Fm)

            DO l=1,m_max+1
               product_list(nproducts)%Fm(l)=product_list(nproducts)%Fm(l)*&
                           (potential_parameter%scale_coulomb+potential_parameter%scale_longrange)- &
                         Fm(l)*potential_parameter%scale_longrange
            ENDDO

            ! longrange
            omega2 = potential_parameter%omega**2
            omega_corr2 = omega2/(omega2+Rho)
            omega_corr = SQRT(omega_corr2)
            T = T*omega_corr2
            CALL fgamma(m_max,T,Fm)
            tmp = omega_corr
            DO l=1,m_max+1
              product_list(nproducts)%Fm(l)= product_list(nproducts)%Fm(l) + Fm(l)*tmp*potential_parameter%scale_longrange
              tmp = tmp * omega_corr2
            END DO
            factor = 2.0_dp*Pi*RhoInv

          CASE(do_hfx_potential_gaussian)
            omega2 = potential_parameter%omega**2
            T = -omega2*T/(Rho+omega2)
            tmp = 1.0_dp
            DO l=1,m_max+1
              product_list(nproducts)%Fm(l) = EXP(T) * tmp
              tmp = tmp * omega2/(Rho+omega2)
            END DO
            factor = (Pi/(Rho+omega2))**(1.5_dp)
          CASE(do_hfx_potential_mix_lg)
            omega2 = potential_parameter%omega**2
            omega_corr2 = omega2/(omega2+Rho)
            omega_corr = SQRT(omega_corr2)
            T = T*omega_corr2
            CALL fgamma(m_max,T,Fm)
            tmp = omega_corr * 2.0_dp*Pi*RhoInv * potential_parameter%scale_longrange
            DO l=1,m_max+1
              Fm(l)= Fm(l)*tmp
              tmp = tmp * omega_corr2
            END DO
            T = Rho*rpq2
            T = -omega2*T/(Rho+omega2)
            tmp = (Pi/(Rho+omega2))**(1.5_dp) * potential_parameter%scale_gaussian
            DO l=1,m_max+1
              product_list(nproducts)%Fm(l) = EXP(T) * tmp + Fm(l)
              tmp = tmp * omega2/(Rho+omega2)
            END DO
            factor = 1.0_dp
          CASE(do_hfx_potential_id)
            product_list(nproducts)%Fm(1) = (Pi*RhoInv)**(1.5_dp)
            DO l=2,m_max+1
             product_list(nproducts)%Fm(l)= 0.0_dp
            END DO
            factor = 1.0_dp
          END SELECT

          tmp = (Pi*ZetapEtaInv)**3
          factor = factor*S1234*SQRT(tmp)

          DO l=1,m_max+1
            product_list(nproducts)%Fm(l)=product_list(nproducts)%Fm(l)*factor
          END DO

          W = (Zeta1*P+Eta*Q)*ZetapEtaInv
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
          product_list(nproducts)%AB = ra-rb
          product_list(nproducts)%CD = C11-tmp_D
        END DO
      END DO
    END DO

  END SUBROUTINE build_pgf_product_list


  SUBROUTINE build_pair_list_pgf(npgfa, npgfb,list, zeta, zetb, screen1, screen2, pgf, R_pgf,&
                                 log10_pmax, log10_eps_schwarz, ra, rb, nelements,&
                                 neighbor_cells, nimages, do_periodic)
    INTEGER, INTENT(IN)                      :: npgfa, npgfb
    TYPE(hfx_pgf_list), &
      DIMENSION(npgfa*npgfb)                 :: list
    REAL(dp), DIMENSION(1:npgfa), INTENT(IN) :: zeta
    REAL(dp), DIMENSION(1:npgfb), INTENT(IN) :: zetb
    REAL(dp), INTENT(IN)                     :: screen1(2), screen2(2)
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: pgf, R_pgf
    REAL(dp), INTENT(IN)                     :: log10_pmax, &
                                                log10_eps_schwarz, ra(3), &
                                                rb(3)
    INTEGER, INTENT(OUT)                     :: nelements
    TYPE(hfx_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    INTEGER                                  :: nimages(npgfa*npgfb)
    LOGICAL, INTENT(IN)                      :: do_periodic

    INTEGER                                  :: element_counter, i, ipgf, j, &
                                                jpgf
    REAL(dp)                                 :: AB(3), im_B(3), pgf_max, &
                                                rab2, Zeta1, Zeta_A, Zeta_B, &
                                                ZetaInv

    nimages = 0
    ! ** inner loop may never be reached
    nelements = npgfa*npgfb
    DO i = 1,SIZE(neighbor_cells)
      IF( do_periodic ) THEN
        im_B = rb + neighbor_cells(i)%cell_r(:)
      ELSE
        im_B = rb
      END IF
      AB = ra - im_B
      rab2 = AB(1)**2+AB(2)**2+AB(3)**2
      IF(screen1(1)*rab2+screen1(2) + screen2(2) + log10_pmax < log10_eps_schwarz) CYCLE
      element_counter = 0
      DO ipgf = 1,npgfa
        DO jpgf = 1,npgfb
          element_counter = element_counter + 1
          pgf_max = pgf(jpgf,ipgf)%x(1)*rab2+pgf(jpgf,ipgf)%x(2)
          IF(pgf_max + screen2(2) + log10_pmax < log10_eps_schwarz) THEN
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

          IF( nimages(element_counter) == 1 ) THEN
            list(element_counter)%ipgf = ipgf
            list(element_counter)%jpgf = jpgf
            list(element_counter)%zetaInv = ZetaInv
            list(element_counter)%zetapzetb = Zeta1
            list(element_counter)%zeta = Zeta_A
            list(element_counter)%zetb = Zeta_B
          END IF

          list(element_counter)%image_list(nimages(element_counter))%S1234 =  (-Zeta_A*Zeta_B*ZetaInv*rab2)
          list(element_counter)%image_list(nimages(element_counter))%P =  (Zeta_A*ra + Zeta_B*im_B)*ZetaInv
          list(element_counter)%image_list(nimages(element_counter))%R = &
                 MAX(0.0_dp, R_pgf(jpgf,ipgf)%x(1)*rab2 + R_pgf(jpgf,ipgf)%x(2))
          list(element_counter)%image_list(nimages(element_counter))%ra = ra
          list(element_counter)%image_list(nimages(element_counter))%rb = im_B
          list(element_counter)%image_list(nimages(element_counter))%rab2 = rab2
        END DO
      END DO
      nelements = MAX(nelements,element_counter)
    END DO
    DO j = 1,nelements
      list(j)%nimages = nimages(j)
    END DO
    ! ** Remove unused elements

    element_counter = 0
    DO j=1,nelements
      IF( list(j)%nimages == 0 ) CYCLE
      element_counter = element_counter + 1
      list(element_counter)%nimages = list(j)%nimages
      list(element_counter)%zetapzetb = list(j)%zetapzetb
      list(element_counter)%ZetaInv = list(j)%ZetaInv
      list(element_counter)%zeta = list(j)%zeta
      list(element_counter)%zetb = list(j)%zetb
      list(element_counter)%ipgf = list(j)%ipgf
      list(element_counter)%jpgf = list(j)%jpgf
      DO i = 1,list(j)%nimages
        list(element_counter)%image_list(i) = list(j)%image_list(i)
      END DO
    END DO

    nelements = element_counter

  END SUBROUTINE build_pair_list_pgf

  SUBROUTINE build_pair_list(natom, list,set_list,i_start,i_end,j_start,j_end,kind_of,basis_parameter,particle_set, &
                             do_periodic,coeffs_set,coeffs_kind,coeffs_kind_max0,log10_eps_schwarz,cell,&
                             pmax_blocks, atomic_pair_list)

    INTEGER, INTENT(IN)                      :: natom
    TYPE(pair_list_type), INTENT(OUT)        :: list
    TYPE(pair_set_list_type), DIMENSION(*), &
      INTENT(OUT)                            :: set_list
    INTEGER, INTENT(IN)                      :: i_start, i_end, j_start, j_end
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    LOGICAL, INTENT(IN)                      :: do_periodic
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :), POINTER         :: coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :)                        :: coeffs_kind
    REAL(KIND=dp), INTENT(IN)                :: coeffs_kind_max0, &
                                                log10_eps_schwarz
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp)                                 :: pmax_blocks
    LOGICAL, DIMENSION(natom, natom)         :: atomic_pair_list

    INTEGER                                  :: iatom, ikind, iset, jatom, &
                                                jkind, jset, n_element, &
                                                nset_ij, nseta, nsetb
    REAL(KIND=dp)                            :: rab2
    REAL(KIND=dp), DIMENSION(3)              :: B11, pbc_B, ra, rb, temp

    n_element=0
    nset_ij=0

    DO iatom=i_start,i_end
      DO jatom=j_start,j_end
        IF(atomic_pair_list(jatom,iatom) .EQV. .FALSE. ) CYCLE

        ikind = kind_of(iatom)
        nseta =  basis_parameter(ikind)%nset
        ra=particle_set(iatom)%r(:)

        IF(jatom<iatom) CYCLE
        jkind = kind_of(jatom)
        nsetb =  basis_parameter(jkind)%nset
        rb=particle_set(jatom)%r(:)

        IF( do_periodic ) THEN
          temp = rb -ra
          pbc_B = pbc(temp,cell)
          B11 = ra + pbc_B
          rab2 = (ra(1)-B11(1))**2+(ra(2)-B11(2))**2+(ra(3)-B11(3))**2
        ELSE
          rab2 = (ra(1)-rb(1))**2+(ra(2)-rb(2))**2+(ra(3)-rb(3))**2
          B11 = rb ! ra - rb
        END IF
        IF( (coeffs_kind(jkind,ikind)%x(1)*rab2+&
             coeffs_kind(jkind,ikind)%x(2) ) + coeffs_kind_max0 + pmax_blocks  < log10_eps_schwarz ) CYCLE

        n_element=n_element+1
        list%elements(n_element)%pair(1)=iatom
        list%elements(n_element)%pair(2)=jatom
        list%elements(n_element)%kind_pair(1)=ikind
        list%elements(n_element)%kind_pair(2)=jkind
        list%elements(n_element)%r1   =ra
        list%elements(n_element)%r2   =B11
        list%elements(n_element)%dist2=rab2
        ! build a list of guaranteed overlapping sets
        list%elements(n_element)%set_bounds(1)=nset_ij+1
        DO iset = 1, nseta
          DO jset = 1,nsetb
            IF (coeffs_set(jset,iset,jkind,ikind)%x(1)*rab2 + coeffs_set(jset,iset,jkind,ikind)%x(2) + &
                coeffs_kind_max0 + pmax_blocks < log10_eps_schwarz ) CYCLE
            nset_ij=nset_ij+1
            set_list(nset_ij)%pair(1)=iset
            set_list(nset_ij)%pair(2)=jset
          END DO
        END DO
        list%elements(n_element)%set_bounds(2)=nset_ij
      END DO
    END DO

    list%n_element=n_element

  END SUBROUTINE build_pair_list

  SUBROUTINE build_atomic_pair_list(natom, atomic_pair_list, kind_of,basis_parameter,particle_set, &
                                   do_periodic,coeffs_kind,coeffs_kind_max0,log10_eps_schwarz,cell,&
                                   blocks)
    INTEGER, INTENT(IN)                      :: natom
    LOGICAL, DIMENSION(natom, natom)         :: atomic_pair_list
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    LOGICAL, INTENT(IN)                      :: do_periodic
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :)                        :: coeffs_kind
    REAL(KIND=dp), INTENT(IN)                :: coeffs_kind_max0, &
                                                log10_eps_schwarz
    TYPE(cell_type), POINTER                 :: cell
    TYPE(hfx_block_range_type), &
      DIMENSION(:), POINTER                  :: blocks

    INTEGER :: iatom, iatom_end, iatom_start, iblock, ikind, jatom, &
      jatom_end, jatom_start, jblock, jkind, nseta, nsetb
    REAL(KIND=dp)                            :: rab2
    REAL(KIND=dp), DIMENSION(3)              :: B11, pbc_B, ra, rb, temp

    atomic_pair_list = .FALSE.

    DO iblock = 1,SIZE(blocks)
      iatom_start = blocks(iblock)%istart
      iatom_end = blocks(iblock)%iend
      DO jblock = 1,SIZE(blocks)
        jatom_start = blocks(jblock)%istart
        jatom_end = blocks(jblock)%iend

        DO iatom = iatom_start, iatom_end
          ikind = kind_of(iatom)
          nseta =  basis_parameter(ikind)%nset
          ra=particle_set(iatom)%r(:)
          DO jatom = jatom_start, jatom_end
            IF(jatom<iatom) CYCLE
            jkind = kind_of(jatom)
            nsetb =  basis_parameter(jkind)%nset
            rb=particle_set(jatom)%r(:)

            IF( do_periodic ) THEN
              temp = rb -ra
              pbc_B = pbc(temp,cell)
              B11 = ra + pbc_B
              rab2 = (ra(1)-B11(1))**2+(ra(2)-B11(2))**2+(ra(3)-B11(3))**2
            ELSE
              rab2 = (ra(1)-rb(1))**2+(ra(2)-rb(2))**2+(ra(3)-rb(3))**2
              B11 = rb ! ra - rb
            END IF
            IF( (coeffs_kind(jkind,ikind)%x(1)*rab2+&
                 coeffs_kind(jkind,ikind)%x(2) ) + coeffs_kind_max0  < log10_eps_schwarz ) CYCLE

            atomic_pair_list(jatom,iatom) = .TRUE.
            atomic_pair_list(iatom,jatom) = .TRUE.
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE build_atomic_pair_list

  SUBROUTINE build_pair_list_mp2(natom, list,set_list,i_start,i_end,j_start,j_end,kind_of,basis_parameter,particle_set, &
                             do_periodic,coeffs_set,coeffs_kind,coeffs_kind_max0,log10_eps_schwarz,cell,&
                             pmax_blocks, atomic_pair_list, skip_atom_symmetry)

    INTEGER, INTENT(IN)                      :: natom
    TYPE(pair_list_type_mp2)                 :: list
    TYPE(pair_set_list_type), DIMENSION(*), &
      INTENT(OUT)                            :: set_list
    INTEGER, INTENT(IN)                      :: i_start, i_end, j_start, j_end
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    LOGICAL, INTENT(IN)                      :: do_periodic
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :), POINTER         :: coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :)                        :: coeffs_kind
    REAL(KIND=dp), INTENT(IN)                :: coeffs_kind_max0, &
                                                log10_eps_schwarz
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp)                                 :: pmax_blocks
    LOGICAL, DIMENSION(natom, natom)         :: atomic_pair_list
    LOGICAL, OPTIONAL                        :: skip_atom_symmetry

    INTEGER                                  :: iatom, ikind, iset, jatom, &
                                                jkind, jset, n_element, &
                                                nset_ij, nseta, nsetb
    LOGICAL                                  :: my_skip_atom_symmetry
    REAL(KIND=dp)                            :: rab2
    REAL(KIND=dp), DIMENSION(3)              :: B11, pbc_B, ra, rb, temp

    n_element=0
    nset_ij=0

    my_skip_atom_symmetry=.FALSE.
    IF (PRESENT(skip_atom_symmetry)) my_skip_atom_symmetry=skip_atom_symmetry

    DO iatom=i_start,i_end
      DO jatom=j_start,j_end
        IF(atomic_pair_list(jatom,iatom) .EQV. .FALSE. ) CYCLE

        ikind = kind_of(iatom)
        nseta =  basis_parameter(ikind)%nset
        ra=particle_set(iatom)%r(:)

        IF(jatom<iatom .AND. (.NOT.my_skip_atom_symmetry)) CYCLE
        jkind = kind_of(jatom)
        nsetb =  basis_parameter(jkind)%nset
        rb=particle_set(jatom)%r(:)

        IF( do_periodic ) THEN
          temp = rb -ra
          pbc_B = pbc(temp,cell)
          B11 = ra + pbc_B
          rab2 = (ra(1)-B11(1))**2+(ra(2)-B11(2))**2+(ra(3)-B11(3))**2
        ELSE
          rab2 = (ra(1)-rb(1))**2+(ra(2)-rb(2))**2+(ra(3)-rb(3))**2
          B11 = rb ! ra - rb
        END IF
        IF( (coeffs_kind(jkind,ikind)%x(1)*rab2+&
             coeffs_kind(jkind,ikind)%x(2) ) + coeffs_kind_max0 + pmax_blocks  < log10_eps_schwarz ) CYCLE

        n_element=n_element+1
        list%elements(n_element)%pair(1)=iatom
        list%elements(n_element)%pair(2)=jatom
        list%elements(n_element)%kind_pair(1)=ikind
        list%elements(n_element)%kind_pair(2)=jkind
        list%elements(n_element)%r1   =ra
        list%elements(n_element)%r2   =B11
        list%elements(n_element)%dist2=rab2
        ! build a list of guaranteed overlapping sets
        list%elements(n_element)%set_bounds(1)=nset_ij+1
        DO iset = 1, nseta
          DO jset = 1,nsetb
            IF (coeffs_set(jset,iset,jkind,ikind)%x(1)*rab2 + coeffs_set(jset,iset,jkind,ikind)%x(2) + &
                coeffs_kind_max0 + pmax_blocks < log10_eps_schwarz ) CYCLE
            nset_ij=nset_ij+1
            set_list(nset_ij)%pair(1)=iset
            set_list(nset_ij)%pair(2)=jset
          END DO
        END DO
        list%elements(n_element)%set_bounds(2)=nset_ij
      END DO
    END DO

    list%n_element=n_element

  END SUBROUTINE build_pair_list_mp2

END MODULE hfx_pair_list_methods
