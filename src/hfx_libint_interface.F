!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface to the Libint-Library
!> \note
!>      IMPORTANT NOTE : this file currently is for a libint configured for
!>                       LIBINT_MAX_AM 5 and LIBINT_MAX_AM1 4
!> \par History
!>      11.2006 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE hfx_libint_interface

  USE cell_types,                      ONLY: cell_type,&
                                             real_to_scaled
  USE gamma,                           ONLY: fgamma=>fgamma_0
  USE hfx_contraction_methods,         ONLY: contract
  USE hfx_libint_wrapper,              ONLY: get_derivs,&
                                             get_eris
  USE hfx_libint_wrapper_types,        ONLY: lib_deriv,&
                                             lib_int,&
                                             prim_data,&
                                             prim_data_f_size
  USE hfx_types,                       ONLY: hfx_pgf_product_list,&
                                             hfx_potential_type
  USE input_constants,                 ONLY: &
       do_hfx_potential_coulomb, do_hfx_potential_gaussian, &
       do_hfx_potential_id, do_hfx_potential_long, do_hfx_potential_mix_cl, &
       do_hfx_potential_mix_cl_trunc, do_hfx_potential_mix_lg, &
       do_hfx_potential_short, do_hfx_potential_truncated
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: nco
  USE t_c_g0,                          ONLY: t_c_g0_n
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: evaluate_eri,&
         evaluate_deriv_eri,&
         evaluate_eri_screen

  INTEGER, DIMENSION(12), PARAMETER :: full_perm1 = (/1,2,3,4,5,6,7,8,9,10,11,12/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm2 = (/4,5,6,1,2,3,7,8,9,10,11,12/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm3 = (/1,2,3,4,5,6,10,11,12,7,8,9/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm4 = (/4,5,6,1,2,3,10,11,12,7,8,9/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm5 = (/7,8,9,10,11,12,1,2,3,4,5,6/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm6 = (/7,8,9,10,11,12,4,5,6,1,2,3/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm7 = (/10,11,12,7,8,9,1,2,3,4,5,6/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm8 = (/10,11,12,7,8,9,4,5,6,1,2,3/)


!***

  CONTAINS


! *****************************************************************************
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
!> \param prim ...
!> \param R11 ...
!> \param R22 ...
!> \par History
!>      03.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE build_quartet_data_screen(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max, &
                                       potential_parameter, prim, R11,R22)

    REAL(KIND=dp)                            :: A(3), B(3), C(3), D(3)
    REAL(KIND=dp), INTENT(IN)                :: Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: m_max
    TYPE(hfx_potential_type)                 :: potential_parameter
    TYPE(prim_data)                          :: prim
    REAL(dp)                                 :: R11, R22

    INTEGER                                  :: i
    LOGICAL                                  :: use_gamma
    REAL(KIND=dp) :: AB(3), AB2, CD(3), CD2, Eta, EtaInv, factor, omega2, &
      omega_corr, omega_corr2, P(3), PQ(3), PQ2, Q(3), R, R1, R2, Rho, &
      RhoInv, S1234, T, tmp, W(3), Zeta, ZetaInv, ZetapEtaInv
    REAL(KIND=dp), &
      DIMENSION(prim_data_f_size)            :: Fm

    Zeta = Zeta_A + Zeta_B
    ZetaInv = 1.0_dp/Zeta
    Eta  = Zeta_C + Zeta_D
    EtaInv = 1.0_dp/Eta
    ZetapEtaInv = Zeta+Eta
    ZetapEtaInv = 1.0_dp/ZetapEtaInv
    Rho  = Zeta*Eta*ZetapEtaInv
    RhoInv = 1.0_dp/Rho

    DO i=1,3
      P(i) = (Zeta_A*A(i) + Zeta_B*B(i))*ZetaInv
      Q(i) = (Zeta_C*C(i) + Zeta_D*D(i))*EtaInv
      AB(i) = A(i)-B(i)
      CD(i) = C(i)-D(i)
      PQ(i) = P(i)-Q(i)
      W(i) = (Zeta*P(i) + Eta*Q(i))*ZetapEtaInv
    END DO

    AB2 = DOT_PRODUCT(AB,AB)
    CD2 = DOT_PRODUCT(CD,CD)
    PQ2 = DOT_PRODUCT(PQ,PQ)

    S1234= EXP((-Zeta_A*Zeta_B*ZetaInv*AB2)+(-Zeta_C*Zeta_D*EtaInv*CD2))
    T = Rho*PQ2


    SELECT CASE(potential_parameter%potential_type)
      CASE(do_hfx_potential_truncated)
        R = potential_parameter%cutoff_radius * SQRT(rho)
        R1 = R11
        R2 = R22
        IF( PQ2 > (R1+R2+potential_parameter%cutoff_radius)**2 ) THEN
          RETURN
        END IF
        CALL t_c_g0_n(prim%F(1),use_gamma,R,T,m_max)
        IF( use_gamma ) CALL fgamma(m_max, T, prim%F(1))
        factor = 2.0_dp * Pi * RhoInv
      CASE(do_hfx_potential_coulomb)
        CALL fgamma(m_max,T,prim%F(1))
        factor = 2.0_dp*Pi*RhoInv
      CASE(do_hfx_potential_short)
        R = potential_parameter%cutoff_radius * SQRT(rho)
        R1 = R11
        R2 = R22
        IF( PQ2 > (R1+R2+potential_parameter%cutoff_radius)**2 ) THEN
          RETURN
        END IF
        CALL fgamma(m_max,T,prim%F(1))
        omega2 = potential_parameter%omega**2
        omega_corr2 = omega2/(omega2+Rho)
        omega_corr = SQRT(omega_corr2)
        T = T*omega_corr2
        CALL fgamma(m_max,T,Fm)
        tmp = - omega_corr
        DO i=1,m_max+1
          prim%F(i)=prim%F(i) + Fm(i)*tmp
          tmp = tmp * omega_corr2
        END DO
        factor = 2.0_dp*Pi*RhoInv
      CASE(do_hfx_potential_long)
        omega2 = potential_parameter%omega**2
        omega_corr2 = omega2/(omega2+Rho)
        omega_corr = SQRT(omega_corr2)
        T = T*omega_corr2
        CALL fgamma(m_max,T,prim%F(1))
        tmp = omega_corr
        DO i=1,m_max+1
          prim%F(i)= prim%F(i)*tmp
          tmp = tmp * omega_corr2
        END DO
        factor = 2.0_dp*Pi*RhoInv
      CASE(do_hfx_potential_mix_cl)
        CALL fgamma(m_max,T,prim%F(1))
        omega2 = potential_parameter%omega**2
        omega_corr2 = omega2/(omega2+Rho)
        omega_corr = SQRT(omega_corr2)
        T = T*omega_corr2
        CALL fgamma(m_max,T,Fm)
        tmp = omega_corr
        DO i=1,m_max+1
          prim%F(i)= prim%F(i)*potential_parameter%scale_coulomb + Fm(i)*tmp*potential_parameter%scale_longrange
          tmp = tmp * omega_corr2
        END DO
        factor = 2.0_dp*Pi*RhoInv
      CASE(do_hfx_potential_mix_cl_trunc)

        ! truncated
        R = potential_parameter%cutoff_radius * SQRT(rho)
        R1 = R11
        R2 = R22
        IF( PQ2 > (R1+R2+potential_parameter%cutoff_radius)**2 ) THEN
          RETURN
        END IF
        CALL t_c_g0_n(prim%F(1),use_gamma,R,T,m_max)
        IF( use_gamma ) CALL fgamma(m_max, T, prim%F(1))

        ! Coulomb
        CALL fgamma(m_max,T,Fm)

        DO i=1,m_max+1
           prim%F(i)=prim%F(i)*(potential_parameter%scale_coulomb+potential_parameter%scale_longrange)- &
                     Fm(i)*potential_parameter%scale_longrange
        ENDDO

        ! longrange
        omega2 = potential_parameter%omega**2
        omega_corr2 = omega2/(omega2+Rho)
        omega_corr = SQRT(omega_corr2)
        T = T*omega_corr2
        CALL fgamma(m_max,T,Fm)
        tmp = omega_corr
        DO i=1,m_max+1
          prim%F(i)= prim%F(i) + Fm(i)*tmp*potential_parameter%scale_longrange
          tmp = tmp * omega_corr2
        END DO
        factor = 2.0_dp*Pi*RhoInv

      CASE(do_hfx_potential_gaussian)
        omega2 = potential_parameter%omega**2
        T = -omega2*T/(Rho+omega2)
        tmp = 1.0_dp
        DO i=1,m_max+1
          prim%F(i) = EXP(T) * tmp
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
        DO i=1,m_max+1
          Fm(i)= Fm(i)*tmp
          tmp = tmp * omega_corr2
        END DO
        T = Rho*PQ2
        T = -omega2*T/(Rho+omega2)
        tmp = (Pi/(Rho+omega2))**(1.5_dp) * potential_parameter%scale_gaussian
        DO i=1,m_max+1
          prim%F(i) = EXP(T) * tmp + Fm(i)
          tmp = tmp * omega2/(Rho+omega2)
        END DO
        factor = 1.0_dp
      CASE(do_hfx_potential_id)
        prim%F(1) = (Pi*RhoInv)**(1.5_dp)
        DO i=2,m_max+1
         prim%F(i)= 0.0_dp
        END DO
        factor = 1.0_dp
    END SELECT

    tmp    = (Pi*ZetapEtaInv)**3
    factor = factor*S1234*SQRT(tmp)

    DO i=1,m_max+1
       prim%F(i)=prim%F(i)*factor
    ENDDO
    prim%U(1:3,1) = P-A
    prim%U(1:3,3) = Q-C
    prim%U(1:3,5) = W-P
    prim%U(1:3,6) = W-Q
    prim%oo2z      = 0.5_dp*ZetaInv
    prim%oo2n      = 0.5_dp*EtaInv
    prim%oo2zn     = 0.5_dp*ZetapEtaInv
    prim%poz       = Rho*ZetaInv
    prim%pon       = Rho*EtaInv
    prim%oo2p      = 0.5_dp*RhoInv
  END SUBROUTINE build_quartet_data_screen

! *****************************************************************************
!> \brief Fill data structure used in libderiv
!> \param prim ...
!> \param A ...
!> \param B ...
!> \param C ...
!> \param D ...
!> \param Zeta_A ...
!> \param Zeta_B ...
!> \param Zeta_C ...
!> \param Zeta_D ...
!> \param ZetaInv ...
!> \param EtaInv ...
!> \param ZetapEtaInv ...
!> \param Rho ...
!> \param RhoInv ...
!> \param P ...
!> \param Q ...
!> \param W ...
!> \par History
!>      03.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE build_deriv_data(prim, A, B, C, D,&
                              Zeta_A, Zeta_B, Zeta_C, Zeta_D,&
                              ZetaInv,EtaInv,ZetapEtaInv,Rho,RhoInv,&
                              P,Q,W)

    TYPE(prim_data)                          :: prim
    REAL(KIND=dp), INTENT(IN)                :: A(3), B(3), C(3), D(3)
    REAL(dp), INTENT(IN)                     :: Zeta_A, Zeta_B, Zeta_C, &
                                                Zeta_D, ZetaInv, EtaInv, &
                                                ZetapEtaInv, Rho, RhoInv, &
                                                P(3), Q(3), W(3)

    prim%U(1:3,1) = P-A
    prim%U(1:3,2) = P-B
    prim%U(1:3,3) = Q-C
    prim%U(1:3,4) = Q-D
    prim%U(1:3,5) = W-P
    prim%U(1:3,6) = W-Q
    prim%twozeta_a = 2.0_dp*Zeta_A
    prim%twozeta_b = 2.0_dp*Zeta_B
    prim%twozeta_c = 2.0_dp*Zeta_C
    prim%twozeta_d = 2.0_dp*Zeta_D
    prim%oo2z      = 0.5_dp*ZetaInv
    prim%oo2n      = 0.5_dp*EtaInv
    prim%oo2zn     = 0.5_dp*ZetapEtaInv
    prim%poz       = Rho*ZetaInv
    prim%pon       = Rho*EtaInv
    prim%oo2p      = 0.5_dp*RhoInv

  END SUBROUTINE build_deriv_data

! *****************************************************************************
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
! *****************************************************************************
  SUBROUTINE evaluate_deriv_eri(deriv, nproducts, pgf_product_list, &
                                n_a,n_b,n_c,n_d,&
                                ncoa,ncob,ncoc,ncod,&
                                nsgfa, nsgfb, nsgfc, nsgfd,&
                                primitives, max_contraction, tmp_max_all, &
                                eps_schwarz, neris, &
                                Zeta_A, Zeta_B, Zeta_C, Zeta_D, ZetaInv,EtaInv, &
                                s_offset_a, s_offset_b, s_offset_c, s_offset_d,&
                                nl_a, nl_b, nl_c, nl_d, nsoa, nsob, nsoc, nsod,&
                                sphi_a, sphi_b, sphi_c, sphi_d,&
                                work,work2,work_forces,buffer1,buffer2,primitives_tmp,&
                                use_virial, work_virial, work2_virial, primitives_tmp_virial,&
                                primitives_virial, cell, tmp_max_all_virial)

    TYPE(lib_deriv)                          :: deriv
    INTEGER, INTENT(IN)                      :: nproducts
    TYPE(hfx_pgf_product_list), &
      DIMENSION(nproducts)                   :: pgf_product_list
    INTEGER, INTENT(IN)                      :: n_a, n_b, n_c, n_d, ncoa, &
                                                ncob, ncoc, ncod, nsgfa, &
                                                nsgfb, nsgfc, nsgfd
    REAL(dp), DIMENSION(nsgfa, nsgfb, nsgfc,&
      nsgfd, 12)                             :: primitives
    REAL(dp)                                 :: max_contraction, tmp_max_all, &
                                                eps_schwarz
    INTEGER(int_8)                           :: neris
    REAL(dp), INTENT(IN)                     :: Zeta_A, Zeta_B, Zeta_C, &
                                                Zeta_D, ZetaInv, EtaInv
    INTEGER                                  :: s_offset_a, s_offset_b, &
                                                s_offset_c, s_offset_d, nl_a, &
                                                nl_b, nl_c, nl_d, nsoa, nsob, &
                                                nsoc, nsod
    REAL(dp), DIMENSION(ncoa, nsoa*nl_a)     :: sphi_a
    REAL(dp), DIMENSION(ncob, nsob*nl_b)     :: sphi_b
    REAL(dp), DIMENSION(ncoc, nsoc*nl_c)     :: sphi_c
    REAL(dp), DIMENSION(ncod, nsod*nl_d)     :: sphi_d
    REAL(dp) :: work(ncoa*ncob*ncoc*ncod,12), work2(ncoa,ncob,ncoc,ncod,12), &
      work_forces(ncoa*ncob*ncoc*ncod,12)
    REAL(dp), DIMENSION(ncoa*ncob*ncoc*ncod) :: buffer1, buffer2
    REAL(dp), DIMENSION(nsoa*nl_a, nsob*nl_b&
      , nsoc*nl_c, nsod*nl_d)                :: primitives_tmp
    LOGICAL, INTENT(IN)                      :: use_virial
    REAL(dp) :: work_virial(ncoa*ncob*ncoc*ncod,12,3), &
      work2_virial(ncoa,ncob,ncoc,ncod,12,3)
    REAL(dp), DIMENSION(nsoa*nl_a, nsob*nl_b&
      , nsoc*nl_c, nsod*nl_d)                :: primitives_tmp_virial
    REAL(dp), DIMENSION(nsgfa, nsgfb, nsgfc,&
      nsgfd, 12, 3)                          :: primitives_virial
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp)                                 :: tmp_max_all_virial

    INTEGER                                  :: a_mysize(1), i, j, k, l, m, &
                                                m_max, mysize, n, p1, p2, p3, &
                                                p4, perm_case
    REAL(dp) :: A(3), AB(3), B(3), C(3), CD(3), D(3), P(3), Q(3), Rho, &
      RhoInv, scoord(12), tmp_max, tmp_max_virial, W(3), ZetapEtaInv
    TYPE(prim_data), TARGET                  :: prim

    m_max = n_a+n_b+n_c+n_d
    m_max = m_max + 1
    mysize = ncoa*ncob*ncoc*ncod
    a_mysize = mysize

    work = 0.0_dp
    work2 = 0.0_dp

    IF( use_virial ) THEN
      work_virial = 0.0_dp
      work2_virial = 0.0_dp
    END IF


    perm_case = 1
    IF(n_a<n_b) THEN
      perm_case = perm_case + 1
    END IF
    IF(n_c<n_d) THEN
      perm_case = perm_case + 2
    END IF
    IF(n_a+n_b > n_c+n_d) THEN
      perm_case = perm_case + 4
    END IF

    SELECT CASE(perm_case)
      CASE(1)
        DO i = 1,nproducts
          A = pgf_product_list(i)%ra
          B = pgf_product_list(i)%rb
          C = pgf_product_list(i)%rc
          D = pgf_product_list(i)%rd
          Rho = pgf_product_list(i)%Rho
          RhoInv = pgf_product_list(i)%RhoInv
          P = pgf_product_list(i)%P
          Q = pgf_product_list(i)%Q
          W = pgf_product_list(i)%W
          AB = pgf_product_list(i)%AB
          CD = pgf_product_list(i)%CD
          ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
          prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

          CALL build_deriv_data(prim, A, B, C, D,&
                                 Zeta_A, Zeta_B, Zeta_C, Zeta_D, &
                                 ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv,&
                                 P, Q, W)
          deriv%AB=AB!A-B
          deriv%CD=CD!C-D
          CALL get_derivs(n_d, n_c, n_b, n_a, deriv, prim, work_forces, a_mysize)
          DO k=4,6
            DO j=1,mysize
              work_forces(j,k) = - 1.0_dp* (work_forces(j,k-3) + &
                                            work_forces(j,k+3) + &
                                            work_forces(j,k+6))
            END DO
          END DO
          DO k=1,12
            DO j=1,mysize
              work(j,k) = work(j,k) + work_forces(j,k)
            END DO
          END DO
          neris = neris + 12 * mysize
          IF( use_virial ) THEN
            CALL real_to_scaled(scoord(1:3),A,cell)
            CALL real_to_scaled(scoord(4:6),B,cell)
            CALL real_to_scaled(scoord(7:9),C,cell)
            CALL real_to_scaled(scoord(10:12),D,cell)
            DO k=1,12
              DO j=1,mysize
                DO m=1,3
                  work_virial(j,k,m) = work_virial(j,k,m) + work_forces(j,k)*scoord(INT((k-1)/3)*3+m)
                END DO
              END DO
            END DO
          END IF
        END DO

        DO n = 1,12
          tmp_max = 0.0_dp
          DO i = 1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i,n)))
          END DO
          tmp_max = tmp_max * max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)

          DO i = 1,ncoa
            p1 = (i-1) *ncob
            DO j = 1,ncob
              p2 = (p1 + j-1)*ncoc
              DO k = 1,ncoc
                p3 = (p2 + k-1)*ncod
                DO l = 1,ncod
                  p4 = p3 + l
                  work2(i,j,k,l,full_perm1(n)) = work(p4,n)
                END DO
              END DO
            END DO
          END DO
        END DO
        IF( use_virial ) THEN
          DO n=1,12
            tmp_max_virial = 0.0_dp
            DO i = 1,mysize
              tmp_max_virial = MAX(tmp_max_virial, &
                               ABS(work_virial(i,n,1)), ABS(work_virial(i,n,2)), ABS(work_virial(i,n,3)))
            END DO
            tmp_max_virial = tmp_max_virial * max_contraction
            tmp_max_all_virial = MAX(tmp_max_all_virial,tmp_max_virial)

            DO i = 1,ncoa
              p1 = (i-1) *ncob
              DO j = 1,ncob
                p2 = (p1 + j-1)*ncoc
                DO k = 1,ncoc
                  p3 = (p2 + k-1)*ncod
                  DO l = 1,ncod
                    p4 = p3 + l
                    work2_virial(i,j,k,l,full_perm1(n),1:3) = work_virial(p4,n,1:3)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END IF
      CASE(2)
        DO i = 1,nproducts
          A = pgf_product_list(i)%ra
          B = pgf_product_list(i)%rb
          C = pgf_product_list(i)%rc
          D = pgf_product_list(i)%rd
          Rho = pgf_product_list(i)%Rho
          RhoInv = pgf_product_list(i)%RhoInv
          P = pgf_product_list(i)%P
          Q = pgf_product_list(i)%Q
          W = pgf_product_list(i)%W
          AB = pgf_product_list(i)%AB
          CD = pgf_product_list(i)%CD
          ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
          prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

          CALL build_deriv_data(prim, B, A, C, D,&
                                 Zeta_B, Zeta_A, Zeta_C, Zeta_D,&
                                 ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv,&
                                 P, Q, W)
          deriv%AB=-AB!B-A
          deriv%CD=CD!C-D
          CALL get_derivs(n_d, n_c, n_a, n_b, deriv, prim, work_forces, a_mysize)
          DO k=4,6
            DO j=1,mysize
              work_forces(j,k) = - 1.0_dp* (work_forces(j,k-3) + &
                                            work_forces(j,k+3) + &
                                            work_forces(j,k+6) )
            ENDDO
          END DO
          DO k=1,12
            DO j=1,mysize
              work(j,k) = work(j,k) + work_forces(j,k)
            END DO
          END DO
          neris = neris + 12 * mysize
          IF( use_virial ) THEN
            CALL real_to_scaled(scoord(1:3),B,cell)
            CALL real_to_scaled(scoord(4:6),A,cell)
            CALL real_to_scaled(scoord(7:9),C,cell)
            CALL real_to_scaled(scoord(10:12),D,cell)
            DO k=1,12
              DO j=1,mysize
                DO m=1,3
                  work_virial(j,k,m) = work_virial(j,k,m) + work_forces(j,k)*scoord(INT((k-1)/3)*3+m)
                END DO
              END DO
            END DO
          END IF

        END DO
        DO n = 1,12
          tmp_max = 0.0_dp
          DO i = 1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i,n)))
          END DO
          tmp_max = tmp_max * max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)

          DO j = 1,ncob
            p1 = (j-1)*ncoa
            DO i = 1,ncoa
              p2 = (p1 + i-1)*ncoc
              DO k = 1,ncoc
                p3 = (p2 + k-1)*ncod
                DO l = 1,ncod
                  p4 = p3 + l
                  work2(i,j,k,l,full_perm2(n)) = work(p4,n)
                END DO
              END DO
            END DO
          END DO
        END DO
        IF( use_virial ) THEN
          DO n = 1,12
            tmp_max_virial = 0.0_dp
            DO i = 1,mysize
              tmp_max_virial = MAX(tmp_max_virial, &
                               ABS(work_virial(i,n,1)), ABS(work_virial(i,n,2)), ABS(work_virial(i,n,3)))
            END DO
            tmp_max_virial = tmp_max_virial * max_contraction
            tmp_max_all_virial = MAX(tmp_max_all_virial,tmp_max_virial)

            DO j = 1,ncob
              p1 = (j-1)*ncoa
              DO i = 1,ncoa
                p2 = (p1 + i-1)*ncoc
                DO k = 1,ncoc
                  p3 = (p2 + k-1)*ncod
                  DO l = 1,ncod
                    p4 = p3 + l
                    work2_virial(i,j,k,l,full_perm2(n),1:3) = work_virial(p4,n,1:3)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END IF
      CASE(3)
        DO i = 1,nproducts
          A = pgf_product_list(i)%ra
          B = pgf_product_list(i)%rb
          C = pgf_product_list(i)%rc
          D = pgf_product_list(i)%rd
          Rho = pgf_product_list(i)%Rho
          RhoInv = pgf_product_list(i)%RhoInv
          P = pgf_product_list(i)%P
          Q = pgf_product_list(i)%Q
          W = pgf_product_list(i)%W
          AB = pgf_product_list(i)%AB
          CD = pgf_product_list(i)%CD
          ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
          prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

          CALL build_deriv_data(prim, A, B, D, C,&
                                 Zeta_A, Zeta_B, Zeta_D, Zeta_C,&
                                 ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv,&
                                 P, Q, W)
          deriv%AB=AB!A-B
          deriv%CD=-CD!D-C
          CALL get_derivs(n_c, n_d, n_b, n_a, deriv, prim, work_forces, a_mysize)
          DO k=4,6
            DO j=1,mysize
              work_forces(j,k) = - 1.0_dp* (work_forces(j,k-3) + &
                                            work_forces(j,k+3) + &
                                            work_forces(j,k+6) )
            END DO
          END DO
          DO k=1,12
            DO j=1,mysize
              work(j,k) = work(j,k) + work_forces(j,k)
            END DO
          END DO
          neris = neris + 12 * mysize
          IF( use_virial ) THEN
            CALL real_to_scaled(scoord(1:3),A,cell)
            CALL real_to_scaled(scoord(4:6),B,cell)
            CALL real_to_scaled(scoord(7:9),D,cell)
            CALL real_to_scaled(scoord(10:12),C,cell)
            DO k=1,12
              DO j=1,mysize
                DO m=1,3
                  work_virial(j,k,m) = work_virial(j,k,m) + work_forces(j,k)*scoord(INT((k-1)/3)*3+m)
                END DO
              END DO
            END DO
          END IF

        END DO
        DO n = 1,12
          tmp_max = 0.0_dp
          DO i = 1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i,n)))
          END DO
          tmp_max = tmp_max * max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)

          DO i = 1,ncoa
            p1 = (i-1)*ncob
            DO j = 1,ncob
              p2 = (p1 + j-1)*ncod
              DO l = 1,ncod
                p3 = (p2 + l-1) * ncoc
                DO k = 1,ncoc
                  p4 = p3+k
                  work2(i,j,k,l,full_perm3(n)) = work(p4,n)
                END DO
              END DO
            END DO
          END DO
        END DO
        IF( use_virial ) THEN
          DO n = 1,12
            tmp_max_virial = 0.0_dp
            DO i = 1,mysize
              tmp_max_virial = MAX(tmp_max_virial, &
                               ABS(work_virial(i,n,1)), ABS(work_virial(i,n,2)), ABS(work_virial(i,n,3)))
            END DO
            tmp_max_virial = tmp_max_virial * max_contraction
            tmp_max_all_virial = MAX(tmp_max_all_virial,tmp_max_virial)

            DO i = 1,ncoa
              p1 = (i-1)*ncob
              DO j = 1,ncob
                p2 = (p1 + j-1)*ncod
                DO l = 1,ncod
                  p3 = (p2 + l-1) * ncoc
                  DO k = 1,ncoc
                    p4 = p3+k
                    work2_virial(i,j,k,l,full_perm3(n),1:3) = work_virial(p4,n,1:3)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END IF
      CASE(4)
        DO i = 1,nproducts
          A = pgf_product_list(i)%ra
          B = pgf_product_list(i)%rb
          C = pgf_product_list(i)%rc
          D = pgf_product_list(i)%rd
          Rho = pgf_product_list(i)%Rho
          RhoInv = pgf_product_list(i)%RhoInv
          P = pgf_product_list(i)%P
          Q = pgf_product_list(i)%Q
          W = pgf_product_list(i)%W
          AB = pgf_product_list(i)%AB
          CD = pgf_product_list(i)%CD
          ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
          prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

          CALL build_deriv_data(prim, B, A, D, C,&
                                 Zeta_B, Zeta_A, Zeta_D, Zeta_C,&
                                 ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv,&
                                 P, Q, W)
          deriv%AB=-AB!B-A
          deriv%CD=-CD!D-C
          CALL get_derivs(n_c, n_d, n_a, n_b, deriv, prim, work_forces, a_mysize)
          DO k=4,6
            DO j=1,mysize
              work_forces(j,k) = - 1.0_dp* (work_forces(j,k-3) + &
                                            work_forces(j,k+3) + &
                                            work_forces(j,k+6) )
            END DO
          END DO
          DO k=1,12
            DO j=1,mysize
              work(j,k) = work(j,k) + work_forces(j,k)
            END DO
          END DO
          neris = neris + 12 * mysize
          IF( use_virial ) THEN
            CALL real_to_scaled(scoord(1:3),B,cell)
            CALL real_to_scaled(scoord(4:6),A,cell)
            CALL real_to_scaled(scoord(7:9),D,cell)
            CALL real_to_scaled(scoord(10:12),C,cell)
            DO k=1,12
              DO j=1,mysize
                DO m=1,3
                  work_virial(j,k,m) = work_virial(j,k,m) + work_forces(j,k)*scoord(INT((k-1)/3)*3+m)
                END DO
              END DO
            END DO
          END IF

        END DO
        DO n = 1,12
          tmp_max = 0.0_dp
          DO i = 1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i,n)))
          END DO
          tmp_max = tmp_max * max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)

          DO j = 1,ncob
            p1 = (j-1)*ncoa
            DO i = 1,ncoa
              p2 = (p1 + i-1)*ncod
              DO l = 1,ncod
                p3 = (p2 + l-1)*ncoc
                DO k = 1,ncoc
                  p4 = p3 + k
                  work2(i,j,k,l,full_perm4(n)) = work(p4,n)
                END DO
              END DO
            END DO
          END DO
        END DO
        IF( use_virial ) THEN
          DO n = 1,12
            tmp_max_virial = 0.0_dp
            DO i = 1,mysize
              tmp_max_virial = MAX(tmp_max_virial, &
                               ABS(work_virial(i,n,1)), ABS(work_virial(i,n,2)), ABS(work_virial(i,n,3)))
            END DO
            tmp_max_virial = tmp_max_virial * max_contraction
            tmp_max_all_virial = MAX(tmp_max_all_virial,tmp_max_virial)

            DO j = 1,ncob
              p1 = (j-1)*ncoa
              DO i = 1,ncoa
                p2 = (p1 + i-1)*ncod
                DO l = 1,ncod
                  p3 = (p2 + l-1)*ncoc
                  DO k = 1,ncoc
                    p4 = p3 + k
                    work2_virial(i,j,k,l,full_perm4(n),1:3) = work_virial(p4,n,1:3)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END IF
      CASE(5)
        DO i = 1,nproducts
          A = pgf_product_list(i)%ra
          B = pgf_product_list(i)%rb
          C = pgf_product_list(i)%rc
          D = pgf_product_list(i)%rd
          Rho = pgf_product_list(i)%Rho
          RhoInv = pgf_product_list(i)%RhoInv
          P = pgf_product_list(i)%P
          Q = pgf_product_list(i)%Q
          W = pgf_product_list(i)%W
          AB = pgf_product_list(i)%AB
          CD = pgf_product_list(i)%CD
          ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
          prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

          CALL build_deriv_data(prim, C, D, A, B,&
                                 Zeta_C ,Zeta_D, Zeta_A, Zeta_B,&
                                 EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv,&
                                 Q, P, W)
          deriv%AB=CD!C-D
          deriv%CD=AB!A-B
          CALL get_derivs(n_b, n_a, n_d, n_c, deriv, prim, work_forces, a_mysize)
          DO k=4,6
            DO j=1,mysize
              work_forces(j,k) = - 1.0_dp* (work_forces(j,k-3) + &
                                            work_forces(j,k+3) + &
                                            work_forces(j,k+6) )
            END DO
          END DO
          DO k=1,12
            DO j=1,mysize
              work(j,k) = work(j,k) + work_forces(j,k)
            END DO
          END DO
          neris = neris + 12 * mysize
          IF( use_virial ) THEN
            CALL real_to_scaled(scoord(1:3),C,cell)
            CALL real_to_scaled(scoord(4:6),D,cell)
            CALL real_to_scaled(scoord(7:9),A,cell)
            CALL real_to_scaled(scoord(10:12),B,cell)
            DO k=1,12
              DO j=1,mysize
                DO m=1,3
                  work_virial(j,k,m) = work_virial(j,k,m) + work_forces(j,k)*scoord(INT((k-1)/3)*3+m)
                END DO
              END DO
            END DO
          END IF

        END DO
        DO n = 1,12
          tmp_max = 0.0_dp
          DO i = 1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i,n)))
          END DO
          tmp_max = tmp_max * max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)

          DO k = 1,ncoc
            p1 = (k-1)*ncod
            DO l = 1,ncod
              p2 = (p1 + l-1)*ncoa
              DO i = 1,ncoa
                p3 = (p2 + i-1)*ncob
                DO j = 1,ncob
                  p4 = p3+j
                  work2(i,j,k,l,full_perm5(n)) = work(p4,n)
                END DO
              END DO
            END DO
          END DO
        END DO
        IF( use_virial ) THEN
          DO n = 1,12
            tmp_max_virial = 0.0_dp
            DO i = 1,mysize
              tmp_max_virial = MAX(tmp_max_virial, &
                               ABS(work_virial(i,n,1)), ABS(work_virial(i,n,2)), ABS(work_virial(i,n,3)))
            END DO
            tmp_max_virial = tmp_max_virial * max_contraction
            tmp_max_all_virial = MAX(tmp_max_all_virial,tmp_max_virial)

            DO k = 1,ncoc
              p1 = (k-1)*ncod
              DO l = 1,ncod
                p2 = (p1 + l-1)*ncoa
                DO i = 1,ncoa
                  p3 = (p2 + i-1)*ncob
                  DO j = 1,ncob
                    p4 = p3+j
                    work2_virial(i,j,k,l,full_perm5(n),1:3) = work_virial(p4,n,1:3)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END IF
      CASE(6)
        DO i = 1,nproducts
          A = pgf_product_list(i)%ra
          B = pgf_product_list(i)%rb
          C = pgf_product_list(i)%rc
          D = pgf_product_list(i)%rd
          Rho = pgf_product_list(i)%Rho
          RhoInv = pgf_product_list(i)%RhoInv
          P = pgf_product_list(i)%P
          Q = pgf_product_list(i)%Q
          W = pgf_product_list(i)%W
          AB = pgf_product_list(i)%AB
          CD = pgf_product_list(i)%CD
          ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
          prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

          CALL build_deriv_data(prim, C, D, B, A, &
                                 Zeta_C, Zeta_D, Zeta_B, Zeta_A,&
                                 EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv,&
                                 Q, P, W)
          deriv%AB=CD!C-D
          deriv%CD=-AB!B-A
          CALL get_derivs(n_a, n_b, n_d, n_c, deriv, prim, work_forces, a_mysize)
          DO k=4,6
            DO j=1,mysize
              work_forces(j,k) = - 1.0_dp* (work_forces(j,k-3) + &
                                            work_forces(j,k+3) + &
                                            work_forces(j,k+6) )
            END DO
          END DO
          DO k=1,12
            DO j=1,mysize
              work(j,k) = work(j,k) + work_forces(j,k)
            END DO
          END DO
          neris = neris + 12 * mysize
          IF( use_virial ) THEN
            CALL real_to_scaled(scoord(1:3),C,cell)
            CALL real_to_scaled(scoord(4:6),D,cell)
            CALL real_to_scaled(scoord(7:9),B,cell)
            CALL real_to_scaled(scoord(10:12),A,cell)
            DO k=1,12
              DO j=1,mysize
                DO m=1,3
                  work_virial(j,k,m) = work_virial(j,k,m) + work_forces(j,k)*scoord(INT((k-1)/3)*3+m)
                END DO
              END DO
            END DO
          END IF

        END DO
        DO n = 1,12
          tmp_max = 0.0_dp
          DO i = 1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i,n)))
          END DO
          tmp_max = tmp_max * max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)

          DO k = 1,ncoc
            p1 = (k-1)*ncod
            DO l = 1,ncod
              p2 = (p1 + l-1)*ncob
              DO j = 1,ncob
                p3 = (p2 + j-1)*ncoa
                DO i = 1,ncoa
                  p4 = p3+i
                  work2(i,j,k,l,full_perm6(n)) = work(p4,n)
                END DO
              END DO
            END DO
          END DO
        END DO
        IF( use_virial ) THEN
          DO n = 1,12
            tmp_max_virial = 0.0_dp
            DO i = 1,mysize
              tmp_max_virial = MAX(tmp_max_virial, &
                               ABS(work_virial(i,n,1)), ABS(work_virial(i,n,2)), ABS(work_virial(i,n,3)))
            END DO
            tmp_max_virial = tmp_max_virial * max_contraction
            tmp_max_all_virial = MAX(tmp_max_all_virial,tmp_max_virial)

            DO k = 1,ncoc
              p1 = (k-1)*ncod
              DO l = 1,ncod
                p2 = (p1 + l-1)*ncob
                DO j = 1,ncob
                  p3 = (p2 + j-1)*ncoa
                  DO i = 1,ncoa
                    p4 = p3+i
                    work2_virial(i,j,k,l,full_perm6(n),1:3) = work_virial(p4,n,1:3)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END IF
      CASE(7)
        DO i = 1,nproducts
          A = pgf_product_list(i)%ra
          B = pgf_product_list(i)%rb
          C = pgf_product_list(i)%rc
          D = pgf_product_list(i)%rd
          Rho = pgf_product_list(i)%Rho
          RhoInv = pgf_product_list(i)%RhoInv
          P = pgf_product_list(i)%P
          Q = pgf_product_list(i)%Q
          W = pgf_product_list(i)%W
          AB = pgf_product_list(i)%AB
          CD = pgf_product_list(i)%CD
          ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
          prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

          CALL build_deriv_data(prim, D, C, A, B,&
                                 Zeta_D, Zeta_C, Zeta_A, Zeta_B,&
                                 EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv,&
                                 Q, P, W)
          deriv%AB=-CD!D-C
          deriv%CD=AB!A-B
          CALL get_derivs(n_b, n_a, n_c, n_d, deriv, prim, work_forces, a_mysize)
          DO k=4,6
            DO j=1,mysize
              work_forces(j,k) = - 1.0_dp* (work_forces(j,k-3) + &
                                            work_forces(j,k+3) + &
                                            work_forces(j,k+6) )
            END DO
          END DO
          DO k=1,12
            DO j=1,mysize
              work(j,k) = work(j,k) + work_forces(j,k)
            END DO
          END DO
          neris = neris + 12 * mysize
          IF( use_virial ) THEN
            CALL real_to_scaled(scoord(1:3),D,cell)
            CALL real_to_scaled(scoord(4:6),C,cell)
            CALL real_to_scaled(scoord(7:9),A,cell)
            CALL real_to_scaled(scoord(10:12),B,cell)
            DO k=1,12
              DO j=1,mysize
                DO m=1,3
                  work_virial(j,k,m) = work_virial(j,k,m) + work_forces(j,k)*scoord(INT((k-1)/3)*3+m)
                END DO
              END DO
            END DO
          END IF

        END DO
        DO n = 1,12
          tmp_max = 0.0_dp
          DO i = 1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i,n)))
          END DO
          tmp_max = tmp_max * max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)

          DO l = 1,ncod
            p1 = (l-1)*ncoc
            DO k = 1,ncoc
              p2 = (p1 + k-1) * ncoa
              DO i = 1,ncoa
                p3 = (p2 + i-1) *ncob
                DO j = 1,ncob
                  p4 = p3+j
                  work2(i,j,k,l,full_perm7(n)) = work(p4,n)
                END DO
              END DO
            END DO
          END DO
        END DO
        IF( use_virial ) THEN
          DO n = 1,12
            tmp_max_virial = 0.0_dp
            DO i = 1,mysize
              tmp_max_virial = MAX(tmp_max_virial, &
                               ABS(work_virial(i,n,1)), ABS(work_virial(i,n,2)), ABS(work_virial(i,n,3)))
            END DO
            tmp_max_virial = tmp_max_virial * max_contraction
            tmp_max_all_virial = MAX(tmp_max_all_virial,tmp_max_virial)

            DO l = 1,ncod
              p1 = (l-1)*ncoc
              DO k = 1,ncoc
                p2 = (p1 + k-1) * ncoa
                DO i = 1,ncoa
                  p3 = (p2 + i-1) *ncob
                  DO j = 1,ncob
                    p4 = p3+j
                    work2_virial(i,j,k,l,full_perm7(n),1:3) = work_virial(p4,n,1:3)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END IF
      CASE(8)
        DO i = 1,nproducts
          A = pgf_product_list(i)%ra
          B = pgf_product_list(i)%rb
          C = pgf_product_list(i)%rc
          D = pgf_product_list(i)%rd
          Rho = pgf_product_list(i)%Rho
          RhoInv = pgf_product_list(i)%RhoInv
          P = pgf_product_list(i)%P
          Q = pgf_product_list(i)%Q
          W = pgf_product_list(i)%W
          AB = pgf_product_list(i)%AB
          CD = pgf_product_list(i)%CD
          ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
          prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

          CALL build_deriv_data(prim, D, C, B, A,&
                                 Zeta_D, Zeta_C, Zeta_B, Zeta_A,&
                                 EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv,&
                                 Q, P, W)
          deriv%AB=-CD!D-C
          deriv%CD=-AB!B-A
          CALL get_derivs(n_a, n_b, n_c, n_d, deriv, prim, work_forces, a_mysize)
          DO k=4,6
            DO j=1,mysize
              work_forces(j,k) = - 1.0_dp* (work_forces(j,k-3) + &
                                            work_forces(j,k+3) + &
                                            work_forces(j,k+6) )
            END DO
          END DO
          DO k=1,12
            DO j=1,mysize
              work(j,k) = work(j,k) + work_forces(j,k)
            END DO
          END DO
          neris = neris + 12 * mysize
          IF( use_virial ) THEN
            CALL real_to_scaled(scoord(1:3),D,cell)
            CALL real_to_scaled(scoord(4:6),C,cell)
            CALL real_to_scaled(scoord(7:9),B,cell)
            CALL real_to_scaled(scoord(10:12),A,cell)
            DO k=1,12
              DO j=1,mysize
                DO m=1,3
                  work_virial(j,k,m) = work_virial(j,k,m) + work_forces(j,k)*scoord(INT((k-1)/3)*3+m)
                END DO
              END DO
            END DO
          END IF

        END DO
        DO n = 1,12
          tmp_max = 0.0_dp
          DO i = 1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i,n)))
          END DO
          tmp_max = tmp_max * max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)

          DO l = 1,ncod
            p1 = (l-1)*ncoc
            DO k = 1,ncoc
              p2 = (p1 + k-1) * ncob
              DO j = 1,ncob
                p3 = (p2 + j-1) * ncoa
                DO i = 1,ncoa
                  p4 = p3 + i
                  work2(i,j,k,l,full_perm8(n)) = work(p4,n)
                END DO
              END DO
            END DO
          END DO
        END DO
        IF( use_virial ) THEN
          DO n = 1,12
            tmp_max_virial = 0.0_dp
            DO i = 1,mysize
              tmp_max_virial = MAX(tmp_max_virial, &
                               ABS(work_virial(i,n,1)), ABS(work_virial(i,n,2)), ABS(work_virial(i,n,3)))
            END DO
            tmp_max_virial = tmp_max_virial * max_contraction
            tmp_max_all_virial = MAX(tmp_max_all_virial,tmp_max_virial)

            DO l = 1,ncod
              p1 = (l-1)*ncoc
              DO k = 1,ncoc
                p2 = (p1 + k-1) * ncob
                DO j = 1,ncob
                  p3 = (p2 + j-1) * ncoa
                  DO i = 1,ncoa
                    p4 = p3 + i
                    work2_virial(i,j,k,l,full_perm8(n),1:3) = work_virial(p4,n,1:3)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END IF
    END SELECT

    IF( .NOT. use_virial) THEN
      IF( tmp_max_all < eps_schwarz ) RETURN
    END IF

    IF( tmp_max_all >= eps_schwarz ) THEN
      DO i = 1,12
        primitives_tmp(:,:,:,:) = 0.0_dp
        CALL contract(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod,&
                      n_a, n_b, n_c, n_d, nl_a, nl_b, nl_c, nl_d, work2(1,1,1,1,i),&
                      sphi_a, &
                      sphi_b, &
                      sphi_c, &
                      sphi_d, &
                      primitives_tmp(1,1,1,1), &
                      buffer1, buffer2)

        primitives(s_offset_a+1:s_offset_a+nsoa*nl_a,&
                   s_offset_b+1:s_offset_b+nsob*nl_b,&
                   s_offset_c+1:s_offset_c+nsoc*nl_c,&
                   s_offset_d+1:s_offset_d+nsod*nl_d,i) = &
        primitives(s_offset_a+1:s_offset_a+nsoa*nl_a,&
                   s_offset_b+1:s_offset_b+nsob*nl_b,&
                   s_offset_c+1:s_offset_c+nsoc*nl_c,&
                   s_offset_d+1:s_offset_d+nsod*nl_d,i) + primitives_tmp(:,:,:,:)
      END DO
    END IF

    IF( use_virial .AND. tmp_max_all_virial >= eps_schwarz ) THEN
      DO i = 1,12
        DO m=1,3
          primitives_tmp_virial(:,:,:,:) = 0.0_dp
          CALL contract(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod,&
                        n_a, n_b, n_c, n_d, nl_a, nl_b, nl_c, nl_d, work2_virial(1,1,1,1,i,m),&
                        sphi_a, &
                        sphi_b, &
                        sphi_c, &
                        sphi_d, &
                        primitives_tmp_virial(1,1,1,1), &
                        buffer1, buffer2)

          primitives_virial(s_offset_a+1:s_offset_a+nsoa*nl_a,&
                     s_offset_b+1:s_offset_b+nsob*nl_b,&
                     s_offset_c+1:s_offset_c+nsoc*nl_c,&
                     s_offset_d+1:s_offset_d+nsod*nl_d,i,m) = &
          primitives_virial(s_offset_a+1:s_offset_a+nsoa*nl_a,&
                     s_offset_b+1:s_offset_b+nsob*nl_b,&
                     s_offset_c+1:s_offset_c+nsoc*nl_c,&
                     s_offset_d+1:s_offset_d+nsod*nl_d,i,m) + primitives_tmp_virial(:,:,:,:)
        END DO
      END DO
    END IF

  END SUBROUTINE evaluate_deriv_eri

! *****************************************************************************
!> \brief Evaluates max-abs values of  electron repulsion integrals for a primitive quartet
!> \param lib ...
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
! *****************************************************************************
  SUBROUTINE evaluate_eri_screen(lib,A,B,C,D,Zeta_A,Zeta_B,Zeta_C,Zeta_D,&
                                 n_a,n_b,n_c,n_d,&
                                 max_val, potential_parameter, R1, R2,&
                                 p_work)

    TYPE(lib_int)                            :: lib
    REAL(dp), INTENT(IN)                     :: A(3), B(3), C(3), D(3), &
                                                Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: n_a, n_b, n_c, n_d
    REAL(dp), INTENT(INOUT)                  :: max_val
    TYPE(hfx_potential_type)                 :: potential_parameter
    REAL(dp)                                 :: R1, R2
    REAL(dp), DIMENSION(:), POINTER          :: p_work

    INTEGER                                  :: a_mysize(1), i, m_max, &
                                                mysize, perm_case
    TYPE(prim_data), TARGET                  :: prim

    m_max = n_a+n_b+n_c+n_d
    mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)
    a_mysize = mysize

    IF(m_max/=0) THEN
      perm_case = 1
      IF(n_a<n_b) THEN
        perm_case = perm_case + 1
      END IF
      IF(n_c<n_d) THEN
        perm_case = perm_case + 2
      END IF
      IF(n_a+n_b > n_c+n_d) THEN
        perm_case = perm_case + 4
      END IF

      SELECT CASE(perm_case)
        CASE(1)
          CALL build_quartet_data_screen(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                            potential_parameter, prim, R1, R2)
          lib%AB=A-B
          lib%CD=C-D
          CALL get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(2)
          CALL build_quartet_data_screen(B,A,C,D,Zeta_B, Zeta_A, Zeta_C, Zeta_D, m_max,&
                            potential_parameter, prim, R1,R2)
          lib%AB=B-A
          lib%CD=C-D
          CALL get_eris(n_d, n_c, n_a, n_b, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(3)
          CALL build_quartet_data_screen(A,B,D,C,Zeta_A, Zeta_B, Zeta_D, Zeta_C, m_max,&
                            potential_parameter, prim, R1,R2)
          lib%AB=A-B
          lib%CD=D-C
          CALL get_eris(n_c, n_d, n_b, n_a, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(4)
          CALL build_quartet_data_screen(B,A,D,C,Zeta_B, Zeta_A, Zeta_D, Zeta_C, m_max,&
                            potential_parameter, prim, R1,R2)
          lib%AB=B-A
          lib%CD=D-C
          CALL get_eris(n_c, n_d, n_a, n_b, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(5)
          CALL build_quartet_data_screen(C,D,A,B,Zeta_C, Zeta_D, Zeta_A, Zeta_B, m_max,&
                            potential_parameter, prim, R1,R2)
          lib%AB=C-D
          lib%CD=A-B
          CALL get_eris(n_b, n_a, n_d, n_c, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(6)
          CALL build_quartet_data_screen(C,D,B,A,Zeta_C, Zeta_D, Zeta_B, Zeta_A, m_max,&
                            potential_parameter, prim, R1,R2)
          lib%AB=C-D
          lib%CD=B-A
          CALL get_eris(n_a, n_b, n_d, n_c, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(7)
          CALL build_quartet_data_screen(D,C,A,B,Zeta_D, Zeta_C, Zeta_A, Zeta_B, m_max,&
                            potential_parameter, prim, R1,R2)
          lib%AB=D-C
          lib%CD=A-B
          CALL get_eris(n_b, n_a, n_c, n_d, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(8)
          CALL build_quartet_data_screen(D,C,B,A,Zeta_D, Zeta_C, Zeta_B, Zeta_A, m_max,&
                            potential_parameter, prim, R1,R2)
          lib%AB=D-C
          lib%CD=B-A
          CALL get_eris(n_a, n_b, n_c, n_d, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
      END SELECT
    ELSE
      CALL build_quartet_data_screen(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                              potential_parameter, prim, R1,R2)
      max_val = ABS(prim%F(1))
    END IF

  END SUBROUTINE evaluate_eri_screen


! *****************************************************************************
!> \brief Evaluate electron repulsion integrals for a primitive quartet
!> \param lib ...
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
! *****************************************************************************

  SUBROUTINE evaluate_eri(lib, nproducts, pgf_product_list, &
                          n_a,n_b,n_c,n_d,&
                          ncoa,ncob,ncoc,ncod,&
                          nsgfa, nsgfb, nsgfc, nsgfd,&
                          primitives, max_contraction, tmp_max, &
                          eps_schwarz, neris, &
                          ZetaInv,EtaInv, &
                          s_offset_a, s_offset_b, s_offset_c, s_offset_d,&
                          nl_a, nl_b, nl_c, nl_d, nsoa, nsob, nsoc, nsod,&
                          sphi_a, sphi_b, sphi_c, sphi_d, work,work2,buffer1,buffer2, &
                          primitives_tmp, p_work)

    TYPE(lib_int)                            :: lib
    INTEGER, INTENT(IN)                      :: nproducts
    TYPE(hfx_pgf_product_list), &
      DIMENSION(nproducts)                   :: pgf_product_list
    INTEGER, INTENT(IN)                      :: n_a, n_b, n_c, n_d, ncoa, &
                                                ncob, ncoc, ncod, nsgfa, &
                                                nsgfb, nsgfc, nsgfd
    REAL(dp), &
      DIMENSION(nsgfa, nsgfb, nsgfc, nsgfd)  :: primitives
    REAL(dp)                                 :: max_contraction, tmp_max, &
                                                eps_schwarz
    INTEGER(int_8)                           :: neris
    REAL(dp), INTENT(IN)                     :: ZetaInv, EtaInv
    INTEGER                                  :: s_offset_a, s_offset_b, &
                                                s_offset_c, s_offset_d, nl_a, &
                                                nl_b, nl_c, nl_d, nsoa, nsob, &
                                                nsoc, nsod
    REAL(dp), DIMENSION(ncoa, nsoa*nl_a), &
      INTENT(IN)                             :: sphi_a
    REAL(dp), DIMENSION(ncob, nsob*nl_b), &
      INTENT(IN)                             :: sphi_b
    REAL(dp), DIMENSION(ncoc, nsoc*nl_c), &
      INTENT(IN)                             :: sphi_c
    REAL(dp), DIMENSION(ncod, nsod*nl_d), &
      INTENT(IN)                             :: sphi_d
    REAL(dp)                                 :: work(ncoa*ncob*ncoc*ncod), &
                                                work2(ncoa,ncob,ncoc,ncod)
    REAL(dp), DIMENSION(ncoa*ncob*ncoc*ncod) :: buffer1, buffer2
    REAL(dp), DIMENSION(nsoa*nl_a, nsob*nl_b&
      , nsoc*nl_c, nsod*nl_d)                :: primitives_tmp
    REAL(dp), DIMENSION(:), POINTER          :: p_work

    INTEGER                                  :: a_mysize(1), i, j, k, l, &
                                                m_max, mysize, p1, p2, p3, &
                                                p4, perm_case
    REAL(dp)                                 :: A(3), AB(3), B(3), C(3), &
                                                CD(3), D(3), P(3), Q(3), Rho, &
                                                RhoInv, W(3), ZetapEtaInv
    TYPE(prim_data), TARGET                  :: prim

    m_max = n_a+n_b+n_c+n_d
    mysize = ncoa*ncob*ncoc*ncod
    a_mysize = mysize

    work = 0.0_dp
    IF(m_max/=0) THEN
      perm_case = 1
      IF(n_a<n_b) THEN
        perm_case = perm_case + 1
      END IF
      IF(n_c<n_d) THEN
        perm_case = perm_case + 2
      END IF
      IF(n_a+n_b > n_c+n_d) THEN
        perm_case = perm_case + 4
      END IF
      SELECT CASE(perm_case)
        CASE(1)
          DO i = 1,nproducts
            A = pgf_product_list(i)%ra
            B = pgf_product_list(i)%rb
            C = pgf_product_list(i)%rc
            D = pgf_product_list(i)%rd
            Rho = pgf_product_list(i)%Rho
            RhoInv = pgf_product_list(i)%RhoInv
            P = pgf_product_list(i)%P
            Q = pgf_product_list(i)%Q
            W = pgf_product_list(i)%W
            AB = pgf_product_list(i)%AB
            CD = pgf_product_list(i)%CD
            ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
            prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

            CALL build_quartet_data(prim, A, C,&
                                    ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv,&
                                    P, Q, W)
            lib%AB=AB!A-B
            lib%CD=CD!C-D
            CALL get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
            work(1:mysize) = work(1:mysize) + p_work(1:mysize)
            neris = neris + mysize
          END DO
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) THEN
            RETURN
          END IF

          DO i = 1,ncoa
            p1 = (i-1) *ncob
            DO j = 1,ncob
              p2 = (p1 + j-1)*ncoc
              DO k = 1,ncoc
                p3 = (p2 + k-1)*ncod
                DO l = 1,ncod
                  p4 = p3 + l
                  work2(i,j,k,l) = work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(2)
          DO i = 1,nproducts
            A = pgf_product_list(i)%ra
            B = pgf_product_list(i)%rb
            C = pgf_product_list(i)%rc
            D = pgf_product_list(i)%rd
            Rho = pgf_product_list(i)%Rho
            RhoInv = pgf_product_list(i)%RhoInv
            P = pgf_product_list(i)%P
            Q = pgf_product_list(i)%Q
            W = pgf_product_list(i)%W
            AB = pgf_product_list(i)%AB
            CD = pgf_product_list(i)%CD
            ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
            prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

            CALL build_quartet_data(prim, B, C,&
                                    ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv,&
                                    P,Q,W)
            lib%AB=-AB!B-A
            lib%CD=CD!C-D
            CALL get_eris(n_d, n_c, n_a, n_b, lib, prim, p_work, a_mysize)
            work(1:mysize) = work(1:mysize) + p_work(1:mysize)
            neris = neris + mysize
          END DO
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) THEN
            RETURN
          END IF

          DO j = 1,ncob
            p1 = (j-1)*ncoa
            DO i = 1,ncoa
              p2 = (p1 + i-1)*ncoc
              DO k = 1,ncoc
                p3 = (p2 + k-1)*ncod
                DO l = 1,ncod
                  p4 = p3 + l
                  work2(i,j,k,l) = work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(3)
          DO i = 1,nproducts
            A = pgf_product_list(i)%ra
            B = pgf_product_list(i)%rb
            C = pgf_product_list(i)%rc
            D = pgf_product_list(i)%rd
            Rho = pgf_product_list(i)%Rho
            RhoInv = pgf_product_list(i)%RhoInv
            P = pgf_product_list(i)%P
            Q = pgf_product_list(i)%Q
            W = pgf_product_list(i)%W
            AB = pgf_product_list(i)%AB
            CD = pgf_product_list(i)%CD
            ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
            prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

            CALL build_quartet_data(prim, A, D,&
                                    ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv,&
                                    P,Q,W)
            lib%AB=AB!A-B
            lib%CD=-CD!D-C
            CALL get_eris(n_c, n_d, n_b, n_a, lib, prim, p_work, a_mysize)
            work(1:mysize) = work(1:mysize) + p_work(1:mysize)
            neris = neris + mysize
          END DO
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) THEN
            RETURN
          END IF

          DO i = 1,ncoa
            p1 = (i-1)*ncob
            DO j = 1,ncob
              p2 = (p1 + j-1)*ncod
              DO l = 1,ncod
                p3 = (p2 + l-1) * ncoc
                DO k = 1,ncoc
                  p4 = p3+k
                  work2(i,j,k,l) = work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(4)
          DO i = 1,nproducts
            A = pgf_product_list(i)%ra
            B = pgf_product_list(i)%rb
            C = pgf_product_list(i)%rc
            D = pgf_product_list(i)%rd
            Rho = pgf_product_list(i)%Rho
            RhoInv = pgf_product_list(i)%RhoInv
            P = pgf_product_list(i)%P
            Q = pgf_product_list(i)%Q
            W = pgf_product_list(i)%W
            AB = pgf_product_list(i)%AB
            CD = pgf_product_list(i)%CD
            ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
            prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

            CALL build_quartet_data(prim, B, D,&
                                    ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv,&
                                    P, Q, W)
            lib%AB=-AB!B-A
            lib%CD=-CD!D-C
            CALL get_eris(n_c, n_d, n_a, n_b, lib, prim, p_work, a_mysize)
            work(1:mysize) = work(1:mysize) + p_work(1:mysize)
            neris = neris + mysize
          END DO
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) THEN
            RETURN
          END IF

          DO j = 1,ncob
            p1 = (j-1)*ncoa
            DO i = 1,ncoa
              p2 = (p1 + i-1)*ncod
              DO l = 1,ncod
                p3 = (p2 + l-1)*ncoc
                DO k = 1,ncoc
                  p4 = p3 + k
                  work2(i,j,k,l) = work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(5)
          DO i = 1,nproducts
            A = pgf_product_list(i)%ra
            B = pgf_product_list(i)%rb
            C = pgf_product_list(i)%rc
            D = pgf_product_list(i)%rd
            Rho = pgf_product_list(i)%Rho
            RhoInv = pgf_product_list(i)%RhoInv
            P = pgf_product_list(i)%P
            Q = pgf_product_list(i)%Q
            W = pgf_product_list(i)%W
            AB = pgf_product_list(i)%AB
            CD = pgf_product_list(i)%CD
            ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
            prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

            CALL build_quartet_data(prim, C, A,&
                                    EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv,&
                                    Q, P, W)
            lib%AB=CD!C-D
            lib%CD=AB!A-B
            CALL get_eris(n_b, n_a, n_d, n_c, lib, prim, p_work, a_mysize)
            work(1:mysize) = work(1:mysize) + p_work(1:mysize)
            neris = neris + mysize
          END DO
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) THEN
            RETURN
          END IF

          DO k = 1,ncoc
            p1 = (k-1)*ncod
            DO l = 1,ncod
              p2 = (p1 + l-1)*ncoa
              DO i = 1,ncoa
                p3 = (p2 + i-1)*ncob
                DO j = 1,ncob
                  p4 = p3+j
                  work2(i,j,k,l) = work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(6)
          DO i = 1,nproducts
            A = pgf_product_list(i)%ra
            B = pgf_product_list(i)%rb
            C = pgf_product_list(i)%rc
            D = pgf_product_list(i)%rd
            Rho = pgf_product_list(i)%Rho
            RhoInv = pgf_product_list(i)%RhoInv
            P = pgf_product_list(i)%P
            Q = pgf_product_list(i)%Q
            W = pgf_product_list(i)%W
            AB = pgf_product_list(i)%AB
            CD = pgf_product_list(i)%CD
            ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
            prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

            CALL build_quartet_data(prim, C, B, &
                                    EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv,&
                                    Q, P, W)
            lib%AB=CD!C-D
            lib%CD=-AB!B-A
            CALL get_eris(n_a, n_b, n_d, n_c, lib, prim, p_work, a_mysize)
            work(1:mysize) = work(1:mysize) + p_work(1:mysize)
            neris = neris + mysize
          END DO
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) THEN
            RETURN
          END IF

          DO k = 1,ncoc
            p1 = (k-1)*ncod
            DO l = 1,ncod
              p2 = (p1 + l-1)*ncob
              DO j = 1,ncob
                p3 = (p2 + j-1)*ncoa
                DO i = 1,ncoa
                  p4 = p3+i
                  work2(i,j,k,l) = work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(7)
          DO i = 1,nproducts
            A = pgf_product_list(i)%ra
            B = pgf_product_list(i)%rb
            C = pgf_product_list(i)%rc
            D = pgf_product_list(i)%rd
            Rho = pgf_product_list(i)%Rho
            RhoInv = pgf_product_list(i)%RhoInv
            P = pgf_product_list(i)%P
            Q = pgf_product_list(i)%Q
            W = pgf_product_list(i)%W
            AB = pgf_product_list(i)%AB
            CD = pgf_product_list(i)%CD
            ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
            prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

            CALL build_quartet_data(prim, D, A,&
                                    EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv,&
                                    Q, P, W)
            lib%AB=-CD!D-C
            lib%CD=AB!A-B
            CALL get_eris(n_b, n_a, n_c, n_d, lib, prim, p_work, a_mysize)
            work(1:mysize) = work(1:mysize) + p_work(1:mysize)
            neris = neris + mysize
          END DO
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) THEN
            RETURN
          END IF

          DO l = 1,ncod
            p1 = (l-1)*ncoc
            DO k = 1,ncoc
              p2 = (p1 + k-1) * ncoa
              DO i = 1,ncoa
                p3 = (p2 + i-1) *ncob
                DO j = 1,ncob
                  p4 = p3+j
                  work2(i,j,k,l) = work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(8)
          DO i = 1,nproducts
            A = pgf_product_list(i)%ra
            B = pgf_product_list(i)%rb
            C = pgf_product_list(i)%rc
            D = pgf_product_list(i)%rd
            Rho = pgf_product_list(i)%Rho
            RhoInv = pgf_product_list(i)%RhoInv
            P = pgf_product_list(i)%P
            Q = pgf_product_list(i)%Q
            W = pgf_product_list(i)%W
            AB = pgf_product_list(i)%AB
            CD = pgf_product_list(i)%CD
            ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
            prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

            CALL build_quartet_data(prim, D, B,&
                                    EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv,&
                                    Q, P, W)
            lib%AB=-CD!D-C
            lib%CD=-AB!B-A
            CALL get_eris(n_a, n_b, n_c, n_d, lib, prim, p_work, a_mysize)
            work(1:mysize) = work(1:mysize) + p_work(1:mysize)
            neris = neris + mysize
          END DO
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) THEN
            RETURN
          END IF

          DO l = 1,ncod
            p1 = (l-1)*ncoc
            DO k = 1,ncoc
              p2 = (p1 + k-1) * ncob
              DO j = 1,ncob
                p3 = (p2 + j-1) * ncoa
                DO i = 1,ncoa
                  p4 = p3 + i
                  work2(i,j,k,l) = work(p4)
                END DO
              END DO
            END DO
          END DO
      END SELECT
    ELSE
      DO i = 1,nproducts
        A = pgf_product_list(i)%ra
        B = pgf_product_list(i)%rb
        C = pgf_product_list(i)%rc
        D = pgf_product_list(i)%rd
        Rho = pgf_product_list(i)%Rho
        RhoInv = pgf_product_list(i)%RhoInv
        P = pgf_product_list(i)%P
        Q = pgf_product_list(i)%Q
        W = pgf_product_list(i)%W
        ZetapEtaInv = pgf_product_list(i)%ZetapEtaInv
        prim%F(1:m_max+1) = pgf_product_list(i)%Fm(1:m_max+1)

        CALL build_quartet_data(prim, A, C,&
                                ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv,&
                                P, Q, W)
        work(1) = work(1) + prim%F(1)
        neris = neris + mysize
      END DO
      work2(1,1,1,1) = work(1)
      tmp_max = max_contraction*ABS(work(1))
      IF( tmp_max < eps_schwarz ) RETURN
    END IF

    IF( tmp_max < eps_schwarz ) RETURN
      primitives_tmp = 0.0_dp

      CALL contract(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod,&
                    n_a, n_b, n_c, n_d, nl_a, nl_b, nl_c, nl_d, work2,&
                    sphi_a, &
                    sphi_b, &
                    sphi_c, &
                    sphi_d, &
                    primitives_tmp, &
                    buffer1, buffer2)

      primitives(s_offset_a+1:s_offset_a+nsoa*nl_a,&
                 s_offset_b+1:s_offset_b+nsob*nl_b,&
                 s_offset_c+1:s_offset_c+nsoc*nl_c,&
                 s_offset_d+1:s_offset_d+nsod*nl_d) = &
      primitives(s_offset_a+1:s_offset_a+nsoa*nl_a,&
                 s_offset_b+1:s_offset_b+nsob*nl_b,&
                 s_offset_c+1:s_offset_c+nsoc*nl_c,&
                 s_offset_d+1:s_offset_d+nsod*nl_d) + primitives_tmp(:,:,:,:)

  END SUBROUTINE evaluate_eri

! *****************************************************************************
!> \brief Fill data structure used in libint
!> \param prim ...
!> \param A ...
!> \param C ...
!> \param ZetaInv ...
!> \param EtaInv ...
!> \param ZetapEtaInv ...
!> \param Rho ...
!> \param RhoInv ...
!> \param P ...
!> \param Q ...
!> \param W ...
!> \par History
!>      03.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************

  SUBROUTINE build_quartet_data(prim, A, C,&
                                ZetaInv,EtaInv,ZetapEtaInv,Rho,RhoInv,&
                                P,Q,W)

    TYPE(prim_data)                          :: prim
    REAL(KIND=dp), INTENT(IN)                :: A(3), C(3)
    REAL(dp), INTENT(IN)                     :: ZetaInv, EtaInv, ZetapEtaInv, &
                                                Rho, RhoInv, P(3), Q(3), W(3)

    prim%U(1:3,1) = P-A
    prim%U(1:3,3) = Q-C
    prim%U(1:3,5) = W-P
    prim%U(1:3,6) = W-Q
    prim%oo2z      = 0.5_dp*ZetaInv
    prim%oo2n      = 0.5_dp*EtaInv
    prim%oo2zn     = 0.5_dp*ZetapEtaInv
    prim%poz       = Rho*ZetaInv
    prim%pon       = Rho*EtaInv
    prim%oo2p      = 0.5_dp*RhoInv
  END SUBROUTINE build_quartet_data

END MODULE hfx_libint_interface
