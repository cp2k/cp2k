!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!***** cp2k/libint_interface [1.0] *
!!
!!   NAME
!!     libint_interface
!!
!!   FUNCTION
!!     Interface to the Libint-Library 
!!
!!   NOTES
!!     IMPORTANT NOTE : this file currently is for a libint configured for
!!                      LIBINT_MAX_AM 5 and LIBINT_MAX_AM1 4
!!
!!   AUTHOR
!!     Manuel Guidon 
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!
!****************************************************************************
MODULE libint_interface

  USE gamma,                           ONLY: fgamma
  USE input_constants,                 ONLY: do_hf_potential_coulomb,&
                                             do_hf_potential_long,&
                                             do_hf_potential_short
  USE kinds,                           ONLY: dp
  USE libint_wrapper,                  ONLY: get_derivs,&
                                             get_eris,&
                                             lib_deriv,&
                                             lib_int,&
                                             prim_data
  USE mathconstants
  USE orbital_pointers
  USE x_hf,                            ONLY: x_hf_potential_type
  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC evaluate_eri,&
         evaluate_deriv_eri,&
         evaluate_eri_screen, pbc

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'libint_interface'

  INTEGER, DIMENSION(12), PARAMETER :: full_perm1 = (/1,2,3,4,5,6,7,8,9,10,11,12/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm2 = (/4,5,6,1,2,3,7,8,9,10,11,12/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm3 = (/1,2,3,4,5,6,10,11,12,7,8,9/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm4 = (/4,5,6,1,2,3,10,11,12,7,8,9/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm5 = (/7,8,9,10,11,12,1,2,3,4,5,6/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm6 = (/7,8,9,10,11,12,4,5,6,1,2,3/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm7 = (/10,11,12,7,8,9,1,2,3,4,5,6/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm8 = (/10,11,12,7,8,9,4,5,6,1,2,3/)


  INTEGER                            :: global_m  = -HUGE(global_m)  
  REAL(KIND=dp), DIMENSION(17)       :: global_Fm = 0.0_dp
  REAL(KIND=dp)                      :: global_T  = 0.0_dp

!***
!****************************************************************************

  CONTAINS

!!***** libint_interface/evaluate_eri [1.0] *
!!
!!   NAME
!!     evaluate_eri
!!
!!   FUNCTION
!!     Evaluate electron repulsion integrals for a primitive quartet
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!     08.2007 refactured permutation part [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE evaluate_eri(lib,A,B,C,D,Zeta_A,Zeta_B,Zeta_C,Zeta_D,&
                          n_a,n_b,n_c,n_d,&
                          ncoa,ncob,ncoc,ncod,&
                          offset_a,offset_b,offset_c,offset_d,&
                          primitives,potential_parameter, max_contraction, tmp_max, eps_schwarz)
    
    TYPE(lib_int)                            :: lib
    REAL(dp), INTENT(IN)                     :: A(3), B(3), C(3), D(3), &
                                                Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: n_a, n_b, n_c, n_d, ncoa, &
                                                ncob, ncoc, ncod, offset_a, &
                                                offset_b, offset_c, offset_d
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod)      :: primitives
    TYPE(x_hf_potential_type)                :: potential_parameter
    REAL(dp)                                 :: max_contraction, tmp_max, &
                                                eps_schwarz

    INTEGER                                  :: a_mysize(1), i, j, k, l, &
                                                m_max, mysize, p1, p2, p3, &
                                                perm_case
    REAL(dp), DIMENSION(:), POINTER          :: p_work
    TYPE(prim_data), TARGET                  :: prim

    m_max = n_a+n_b+n_c+n_d
    mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)
    a_mysize = mysize

    !! IF we use a Compiler that does not support whether ISO_C_BINDING nor
    !! CRAY POINTERS, we have to ALLOCATE work array here
    !! See comment on __UGLY_BUT_FAST in libint_wrapper

#if defined (__XT3)
#if defined (__UGLY_BUT_FAST)
    !! do_nothing
#else
    ALLOCATE(p_work(mysize))
#endif
#endif

#if defined (__PGI)
#if defined (__UGLY_BUT_FAST)
    !! do_nothing
#else
    ALLOCATE(p_work(mysize))
#endif
#endif


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
          CALL build_quartet_data(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                            potential_parameter, prim)
          lib%AB=A-B
          lib%CD=C-D
          CALL get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(p_work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) RETURN
          DO i = 1,nco(n_a)
            p1 = (i-1) *nco(n_b)
            DO j = 1,nco(n_b)
              p2 = (p1 + j-1)*nco(n_c)
              DO k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                DO l = 1,nco(n_d)
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)=&
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)+p_work(p3+l)
                END DO
              END DO
            END DO
          END DO
        CASE(2)
          CALL build_quartet_data(B,A,C,D,Zeta_B, Zeta_A, Zeta_C, Zeta_D, m_max,&
                            potential_parameter, prim)
          lib%AB=B-A
          lib%CD=C-D
          CALL get_eris(n_d, n_c, n_a, n_b, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(p_work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) RETURN

          DO j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            DO i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_c)
              DO k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                DO l = 1,nco(n_d)
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)=&
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)+p_work(p3+l)
                END DO
              END DO
            END DO
          END DO
        CASE(3)
          CALL build_quartet_data(A,B,D,C,Zeta_A, Zeta_B, Zeta_D, Zeta_C, m_max,&
                            potential_parameter, prim)
          lib%AB=A-B
          lib%CD=D-C
          CALL get_eris(n_c, n_d, n_b, n_a, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(p_work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) RETURN

          DO i = 1,nco(n_a)
            p1 = (i-1)*nco(n_b)
            DO j = 1,nco(n_b)
              p2 = (p1 + j-1)*nco(n_d)
              DO l = 1,nco(n_d)
                p3 = (p2 + l-1) * nco(n_c)
                DO k = 1,nco(n_c)
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)=&
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)+p_work(p3+k)
                END DO
              END DO
            END DO
          END DO
        CASE(4)
          CALL build_quartet_data(B,A,D,C,Zeta_B, Zeta_A, Zeta_D, Zeta_C, m_max,&
                            potential_parameter, prim)
          lib%AB=B-A
          lib%CD=D-C
          CALL get_eris(n_c, n_d, n_a, n_b, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(p_work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) RETURN

          DO j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            DO i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_d)
              DO l = 1,nco(n_d)
                p3 = (p2 + l-1)*nco(n_c)
                DO k = 1,nco(n_c)
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)=&
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)+p_work(p3+k)
                END DO
              END DO
            END DO
          END DO
        CASE(5)
          CALL build_quartet_data(C,D,A,B,Zeta_C, Zeta_D, Zeta_A, Zeta_B, m_max,&
                            potential_parameter, prim)
          lib%AB=C-D
          lib%CD=A-B
          CALL get_eris(n_b, n_a, n_d, n_c, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(p_work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) RETURN

          DO k = 1,nco(n_c)
            p1 = (k-1)*nco(n_d)
            DO l = 1,nco(n_d)
              p2 = (p1 + l-1)*nco(n_a)
              DO i = 1,nco(n_a)
                p3 = (p2 + i-1)*nco(n_b)
                DO j = 1,nco(n_b)
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)=&
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)+p_work(p3+j)
                END DO
              END DO
            END DO
          END DO
        CASE(6)
          CALL build_quartet_data(C,D,B,A,Zeta_C, Zeta_D, Zeta_B, Zeta_A, m_max,&
                            potential_parameter, prim)
          lib%AB=C-D
          lib%CD=B-A
          CALL get_eris(n_a, n_b, n_d, n_c, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(p_work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) RETURN

          DO k = 1,nco(n_c)
            p1 = (k-1)*nco(n_d)
            DO l = 1,nco(n_d)
              p2 = (p1 + l-1)*nco(n_b)
              DO j = 1,nco(n_b)
                p3 = (p2 + j-1)*nco(n_a)
                DO i = 1,nco(n_a)
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)=&
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)+p_work(p3+i)
                END DO
              END DO
            END DO
          END DO
        CASE(7)
          CALL build_quartet_data(D,C,A,B,Zeta_D, Zeta_C, Zeta_A, Zeta_B, m_max,&
                            potential_parameter, prim)
          lib%AB=D-C
          lib%CD=A-B
          CALL get_eris(n_b, n_a, n_c, n_d, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(p_work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) RETURN

          DO l = 1,nco(n_d)
            p1 = (l-1)*nco(n_c)
            DO k = 1,nco(n_c)
              p2 = (p1 + k-1) * nco(n_a)
              DO i = 1,nco(n_a)
                p3 = (p2 + i-1) *nco(n_b)
                DO j = 1,nco(n_b)
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)=&
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)+p_work(p3+j)
                END DO
              END DO
            END DO
          END DO
        CASE(8)
          CALL build_quartet_data(D,C,B,A,Zeta_D, Zeta_C, Zeta_B, Zeta_A, m_max,&
                            potential_parameter, prim)
          lib%AB=D-C
          lib%CD=B-A
          CALL get_eris(n_a, n_b, n_c, n_d, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(p_work(i)))
          END DO
          tmp_max = tmp_max*max_contraction
          IF(tmp_max<eps_schwarz) RETURN

          DO l = 1,nco(n_d)
            p1 = (l-1)*nco(n_c)
            DO k = 1,nco(n_c)
              p2 = (p1 + k-1) * nco(n_b)
              DO j = 1,nco(n_b)
                p3 = (p2 + j-1) * nco(n_a)
                DO i = 1,nco(n_a)
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)=&
                  primitives(offset_a+i,&
                             offset_b+j,&
                             offset_c+k,&
                             offset_d+l)+p_work(p3+i)
                END DO
              END DO
            END DO
          END DO
      END SELECT
    ELSE
      CALL build_quartet_data(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                            potential_parameter, prim)
      primitives(offset_a+1,offset_b+1,offset_c+1,offset_d+1) = &
        primitives(offset_a+1,offset_b+1,offset_c+1,offset_d+1)+prim%F(1)
      tmp_max = max_contraction*ABS(prim%F(1))
    END IF


    !! IF we use a Compiler that does not support whether ISO_C_BINDING nor
    !! CRAY POINTERS, we have to DEALLOCATE work array here
    !! See comment on __UGLY_BUT_FAST in libint_wrapper

#if defined (__XT3)
#if defined (__UGLY_BUT_FAST)
    !! do_nothing
#else
    DEALLOCATE(p_work)
#endif
#endif

#if defined (__PGI)
#if defined (__UGLY_BUT_FAST)
    !! do_nothing
#else
    DEALLOCATE(p_work)
#endif
#endif

  END SUBROUTINE evaluate_eri

!***************************************************************************

!!***** libint_interface/build_quartet_data [1.0] *
!!
!!   NAME
!!     build_quartet_data
!!
!!   FUNCTION
!!     Fill data structure used in libint
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     03.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE build_quartet_data(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max, &
                              potential_parameter, prim)
     
    REAL(KIND=dp)                            :: A(3), B(3), C(3), D(3)
    REAL(KIND=dp), INTENT(IN)                :: Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: m_max
    TYPE(x_hf_potential_type)                :: potential_parameter
    TYPE(prim_data)                          :: prim

    INTEGER                                  :: i
    REAL(KIND=dp) :: AB(3), AB2, CD(3), CD2, Eta, EtaInv, factor, omega2, &
      omega_corr, omega_corr2, P(3), PQ(3), PQ2, Q(3), Rho, RhoInv, S1234, T, &
      tmp, W(3), Zeta, ZetaInv, ZetapEtaInv
    REAL(KIND=dp), DIMENSION(17)             :: Fm

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
      CASE(do_hf_potential_coulomb)
        CALL fgamma(m_max,T,prim%F)
      CASE(do_hf_potential_short)
        CALL fgamma(m_max,T,prim%F)
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
      CASE(do_hf_potential_long)
        omega2 = potential_parameter%omega**2
        omega_corr2 = omega2/(omega2+Rho)
        omega_corr = SQRT(omega_corr2)
        T = T*omega_corr2
        CALL fgamma(m_max,T,prim%F)
        tmp = omega_corr
        DO i=1,m_max+1
          prim%F(i)= prim%F(i)*tmp
          tmp = tmp * omega_corr2
        END DO
    END SELECT

    factor = 2.0_dp*Pi*RhoInv 
    tmp    = (Pi*ZetapEtaInv)**3
    factor = factor*S1234*SQRT(tmp)

    DO i=1,m_max+1
       prim%F(i)=prim%F(i)*factor
    ENDDO
    prim%U(1:3,1) = P-A
!    prim%U(1:3,2) = P-B                 !Not used in libint
    prim%U(1:3,3) = Q-C
!    prim%U(1:3,4) = Q-D                 !Not used in libint
    prim%U(1:3,5) = W-P
    prim%U(1:3,6) = W-Q
!    prim%twozeta_a = 0.0_dp                      !Not used in libint
!    prim%twozeta_b = 0.0_dp                      !Not used in libint
!    prim%twozeta_c = 0.0_dp                      !Not used in libint 
!    prim%twozeta_d = 0.0_dp                      !Not used in libint
    prim%oo2z      = 0.5_dp*ZetaInv
    prim%oo2n      = 0.5_dp*EtaInv
    prim%oo2zn     = 0.5_dp*ZetapEtaInv
    prim%poz       = Rho*ZetaInv
    prim%pon       = Rho*EtaInv
    prim%oo2p      = 0.5_dp*RhoInv
!    prim%ss_r12_ss = 0.0_dp                      !Not used in libint,libderiv
  END SUBROUTINE build_quartet_data
!***************************************************************************

!!***** libint_interface/build_deriv_data [1.0] *
!!
!!   NAME
!!     build_derive_data
!!
!!   FUNCTION
!!     Fill data structure used in libderiv
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     03.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE build_deriv_data(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                            potential_parameter, prim)
     
    REAL(KIND=dp)                            :: A(3), B(3), C(3), D(3)
    REAL(KIND=dp), INTENT(IN)                :: Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: m_max
    TYPE(x_hf_potential_type)                :: potential_parameter
    TYPE(prim_data)                          :: prim

    INTEGER                                  :: i
    REAL(KIND=dp) :: AB(3), AB2, C1(3), CD(3), CD2, D1(3), Eta, EtaInv, &
      factor, omega2, omega_corr, omega_corr2, P(3), PQ(3), PQ2, Q(3), Rho, &
      RhoInv, S1234, T, tmp, W(3), Zeta, ZetaInv, ZetapEtaInv
    REAL(KIND=dp), DIMENSION(17)             :: Fm

    Zeta = Zeta_A + Zeta_B
    ZetaInv = 1.0_dp/Zeta
    Eta  = Zeta_C + Zeta_D
    EtaInv = 1.0_dp /Eta
    ZetapEtaInv = Zeta+Eta
    ZetapEtaInv = 1.0_dp/ZetapEtaInv
    Rho  = Zeta*Eta*ZetapEtaInv
    RhoInv = 1.0_dp/Rho

    DO i=1,3
      C1(i) = C(i) 
      D1(i) = D(i)
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
      CASE(do_hf_potential_coulomb)
        CALL fgamma(m_max,T,prim%F)
      CASE(do_hf_potential_short)
        CALL fgamma(m_max,T,prim%F)
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
      CASE(do_hf_potential_long)
        omega2 = potential_parameter%omega**2
        omega_corr2 = omega2/(omega2+Rho)
        omega_corr = SQRT(omega_corr2)
        T = T*omega_corr2
        CALL fgamma(m_max,T,prim%F)
        tmp = omega_corr
        DO i=1,m_max+1
          prim%F(i)= prim%F(i)*tmp
          tmp = tmp * omega_corr2
        END DO
    END SELECT

    factor = 2.0_dp*Pi*RhoInv
    tmp = (Pi*ZetapEtaInv)**3
    factor = factor*S1234*SQRT(tmp)

    DO i=1,m_max+1
       prim%F(i)=prim%F(i)*factor
    ENDDO
    prim%U(:,1) = P-A
    prim%U(:,2) = P-B
    prim%U(:,3) = Q-C
    prim%U(:,4) = Q-D
    prim%U(:,5) = W-P
    prim%U(:,6) = W-Q
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
!    prim%ss_r12_ss = 0.0_dp                      !Not used in libint,libderiv

    C = C1
    D = D1
  END SUBROUTINE build_deriv_data
!***************************************************************************

!!***** libint_interface/evaluate_deriv_eri [1.0] *
!!
!!   NAME
!!     evaluate_deriv_eri
!!
!!   FUNCTION
!!     Evaluates derivatives of electron repulsion integrals for a primitive quartet
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     03.2007 created [Manuel Guidon]
!!     08.2007 refactured permutation part [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE evaluate_deriv_eri(deriv,A,B,C,D,Zeta_A,Zeta_B,Zeta_C,Zeta_D,&
                                n_a,n_b,n_c,n_d,work_forces,&
                                ncoa, ncob, ncoc, ncod,&
                                primitive_forces,&
                                offset_a, offset_b, offset_c, offset_d,&
                                potential_parameter, max_contraction, tmp_max_all, eps_schwarz)
    
    TYPE(lib_deriv)                          :: deriv
    REAL(dp), INTENT(IN)                     :: A(3), B(3), C(3), D(3), &
                                                Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: n_a, n_b, n_c, n_d
    REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*&
      nco(n_c)*nco(n_d), 12)                 :: work_forces
    INTEGER, INTENT(IN)                      :: ncoa, ncob, ncoc, ncod
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod, 12)  :: primitive_forces
    INTEGER, INTENT(IN)                      :: offset_a, offset_b, offset_c, &
                                                offset_d
    TYPE(x_hf_potential_type)                :: potential_parameter
    REAL(dp)                                 :: max_contraction, tmp_max_all, &
                                                eps_schwarz

    INTEGER                                  :: a_mysize(1), i, j, k, l, &
                                                m_max, mysize, n, p1, p2, p3, &
                                                perm_case
    REAL(dp)                                 :: tmp_max
    TYPE(prim_data), TARGET                  :: prim

!permutation of configuration

    perm_case = 1
    IF(n_a<n_b) THEN
      perm_case = perm_case + 1
    END IF

    IF(n_c<n_d) THEN
      perm_case = perm_case + 2
    END IF
 
    IF( n_a+n_b > n_c+n_d) THEN
      perm_case = perm_case + 4
    END IF

   
    m_max = n_a+n_b+n_c+n_d
    m_max = m_max + 1
    mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)
    a_mysize=mysize
    SELECT CASE(perm_case)
      CASE(1)
        CALL build_deriv_data(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                          potential_parameter, prim)
        deriv%AB=A-B
        deriv%CD=C-D
        CALL get_derivs(n_d, n_c, n_b, n_a, deriv, prim, work_forces, a_mysize)
        DO k=4,6
          DO i=1,mysize
             work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                                  work_forces(i,k+3) + &
                                                  work_forces(i,k+6) )
          ENDDO
        END DO
        DO n=1,12
          tmp_max = 0.0_dp
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
          END DO
          tmp_max = tmp_max*max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)
          IF(tmp_max<eps_schwarz) CYCLE

          DO i = 1,nco(n_a)
            p1 = (i-1) *nco(n_b)
            DO j = 1,nco(n_b)
              p2 = (p1 + j-1)*nco(n_c)
              DO k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                DO l = 1,nco(n_d)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm1(n))=&
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm1(n))+&
                                   work_forces(p3+l,n)
                END DO
              END DO
            END DO
          END DO
        END DO
      CASE(2)
        CALL build_deriv_data(B,A,C,D,Zeta_B, Zeta_A, Zeta_C, Zeta_D, m_max,&
                           potential_parameter, prim)
        deriv%AB=B-A
        deriv%CD=C-D
        CALL get_derivs(n_d, n_c, n_a, n_b, deriv, prim, work_forces, a_mysize)

        DO k=4,6
          DO i=1,mysize
             work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                                  work_forces(i,k+3) + &
                                                  work_forces(i,k+6) )
          ENDDO
        END DO
        DO n=1,12
          tmp_max = 0.0_dp
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
          END DO
          tmp_max = tmp_max*max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)
          IF(tmp_max<eps_schwarz) CYCLE

          DO j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            DO i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_c)
              DO k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                DO l = 1,nco(n_d)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm2(n))=&
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm2(n))+&
                                   work_forces(p3+l,n)
                END DO
              END DO
            END DO
          END DO
        END DO
     CASE(3)
        CALL build_deriv_data(A,B,D,C,Zeta_A, Zeta_B, Zeta_D, Zeta_C, m_max,&
                          potential_parameter, prim)
        deriv%AB=A-B
        deriv%CD=D-C
        CALL get_derivs(n_c, n_d, n_b, n_a, deriv, prim, work_forces, a_mysize)

        DO k=4,6
           DO i=1,mysize
              work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                               work_forces(i,k+3) + &
                                               work_forces(i,k+6) )
           ENDDO
        END DO
        DO n=1,12
          tmp_max = 0.0_dp
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
          END DO
          tmp_max = tmp_max*max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)
          IF(tmp_max<eps_schwarz) CYCLE

          DO i = 1,nco(n_a)
            p1 = (i-1)*nco(n_b)
            DO j = 1,nco(n_b)
              p2 = (p1 + j-1)*nco(n_d)
              DO l = 1,nco(n_d)
                p3 = (p2 + l-1) * nco(n_c)
                DO k = 1,nco(n_c)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm3(n))=&
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm3(n))+&
                                   work_forces(p3+k,n)
                END DO
              END DO
            END DO
          END DO
        END DO
      CASE(4)
        CALL build_deriv_data(B,A,D,C,Zeta_B, Zeta_A, Zeta_D, Zeta_C, m_max,&
                          potential_parameter, prim)
        deriv%AB=B-A
        deriv%CD=D-C
        CALL get_derivs(n_c, n_d, n_a, n_b, deriv, prim, work_forces, a_mysize)

        DO k=4,6
           DO i=1,mysize
              work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                               work_forces(i,k+3) + &
                                               work_forces(i,k+6) )
           ENDDO
        END DO
        DO n=1,12          
          tmp_max = 0.0_dp
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
          END DO
          tmp_max = tmp_max*max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)
          IF(tmp_max<eps_schwarz) CYCLE

          DO j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            DO i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_d)
              DO l = 1,nco(n_d)
                p3 = (p2 + l-1)*nco(n_c)
                DO k = 1,nco(n_c)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm4(n))=&
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm4(n))+&
                                   work_forces(p3+k,n)
                END DO
              END DO
            END DO
          END DO
        END DO
      CASE(5)
        CALL build_deriv_data(C,D,A,B,Zeta_C, Zeta_D, Zeta_A, Zeta_B, m_max,&
                          potential_parameter, prim)
        deriv%AB=C-D
        deriv%CD=A-B
        CALL get_derivs(n_b, n_a, n_d, n_c, deriv, prim, work_forces, a_mysize)

        DO k=4,6
           DO i=1,mysize
              work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                               work_forces(i,k+3) + &
                                               work_forces(i,k+6) )
           ENDDO
        END DO
        DO n=1,12
          tmp_max = 0.0_dp
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
          END DO
          tmp_max = tmp_max*max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)
          IF(tmp_max<eps_schwarz) CYCLE

          DO k = 1,nco(n_c)
            p1 = (k-1)*nco(n_d)
            DO l = 1,nco(n_d)
              p2 = (p1 + l-1)*nco(n_a)
              DO i = 1,nco(n_a)
                p3 = (p2 + i-1)*nco(n_b)
                DO j = 1,nco(n_b)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm5(n))=&
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm5(n))+&
                                   work_forces(p3+j,n)
                END DO
              END DO
            END DO
          END DO
        END DO
      CASE(6)
        CALL build_deriv_data(C,D,B,A,Zeta_C, Zeta_D, Zeta_B, Zeta_A, m_max,&
                          potential_parameter, prim)
        deriv%AB=C-D
        deriv%CD=B-A
        CALL get_derivs(n_a, n_b, n_d, n_c, deriv, prim, work_forces, a_mysize)

        DO k=4,6
          DO i=1,mysize
             work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                               work_forces(i,k+3) + &
                                               work_forces(i,k+6) )
          ENDDO
        END DO
        DO n=1,12
          tmp_max = 0.0_dp
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
          END DO
          tmp_max = tmp_max*max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)
          IF(tmp_max<eps_schwarz) CYCLE

          DO k = 1,nco(n_c)
            p1 = (k-1)*nco(n_d)
            DO l = 1,nco(n_d)
              p2 = (p1 + l-1)*nco(n_b)
              DO j = 1,nco(n_b)
                p3 = (p2 + j-1)*nco(n_a)
                DO i = 1,nco(n_a)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm6(n))=&
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm6(n))+&
                                   work_forces(p3+i,n)
                END DO
              END DO
            END DO
          END DO
        END DO
      CASE(7)
        CALL build_deriv_data(D,C,A,B,Zeta_D, Zeta_C, Zeta_A, Zeta_B, m_max,&
                          potential_parameter, prim)
        deriv%AB=D-C
        deriv%CD=A-B
        CALL get_derivs(n_b, n_a, n_c, n_d, deriv, prim, work_forces, a_mysize)

        DO k=4,6
           DO i=1,mysize
              work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                               work_forces(i,k+3) + &
                                               work_forces(i,k+6) )
           ENDDO
        END DO
        DO n=1,12
          tmp_max = 0.0_dp
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
          END DO
          tmp_max = tmp_max*max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)
          IF(tmp_max<eps_schwarz) CYCLE

          DO l = 1,nco(n_d)
            p1 = (l-1)*nco(n_c)
            DO k = 1,nco(n_c)
              p2 = (p1 + k-1) * nco(n_a)
              DO i = 1,nco(n_a)
                p3 = (p2 + i-1) *nco(n_b)
                DO j = 1,nco(n_b)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm7(n))=&
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm7(n))+&
                                   work_forces(p3+j,n)
                END DO
              END DO
            END DO
          END DO
        END DO
      CASE(8)
        CALL build_deriv_data(D,C,B,A,Zeta_D, Zeta_C, Zeta_B, Zeta_A, m_max,&
                          potential_parameter, prim)
        deriv%AB=D-C
        deriv%CD=B-A
        CALL get_derivs(n_a, n_b, n_c, n_d, deriv, prim, work_forces, a_mysize)

        DO k=4,6
           DO i=1,mysize
              work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                               work_forces(i,k+3) + &
                                               work_forces(i,k+6) )
           ENDDO
        END DO
        DO n=1,12
          tmp_max = 0.0_dp
          DO i=1,mysize
            tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
          END DO
          tmp_max = tmp_max*max_contraction
          tmp_max_all = MAX(tmp_max_all, tmp_max)
          IF(tmp_max<eps_schwarz) CYCLE

          DO l = 1,nco(n_d)
            p1 = (l-1)*nco(n_c)
            DO k = 1,nco(n_c)
              p2 = (p1 + k-1) * nco(n_b)
              DO j = 1,nco(n_b)
                p3 = (p2 + j-1) * nco(n_a)
                DO i = 1,nco(n_a)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm8(n))=&
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm8(n))+&
                                   work_forces(p3+i,n)
                END DO
              END DO
            END DO
          END DO
        END DO
    END SELECT
  END SUBROUTINE evaluate_deriv_eri

!***************************************************************************

!!***** libint_interface/evaluate_eri_screen [1.0] *
!!
!!   NAME
!!     evaluate_eri_screen
!!
!!   FUNCTION
!!     Evaluates max-abs values of  electron repulsion integrals for a primitive quartet
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     03.2007 created [Manuel Guidon]
!!     08.2007 refactured permutation part [Manuel Guidon]                  
!!
!!****************************************************************************

  SUBROUTINE evaluate_eri_screen(lib,A,B,C,D,Zeta_A,Zeta_B,Zeta_C,Zeta_D,&
                                 n_a,n_b,n_c,n_d,&
                                 max_val, potential_parameter)
    
    TYPE(lib_int)                            :: lib
    REAL(dp), INTENT(IN)                     :: A(3), B(3), C(3), D(3), &
                                                Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: n_a, n_b, n_c, n_d
    REAL(dp), INTENT(INOUT)                  :: max_val
    TYPE(x_hf_potential_type)                :: potential_parameter

    INTEGER                                  :: a_mysize(1), i, m_max, &
                                                mysize, perm_case
    REAL(dp), DIMENSION(:), POINTER          :: p_work
    TYPE(prim_data), TARGET                  :: prim

!permutation of configuration

    m_max = n_a+n_b+n_c+n_d
    mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)
    a_mysize = mysize

    !! IF we use a Compiler that does not support whether ISO_C_BINDING nor
    !! CRAY POINTERS, we have to ALLOCATE work array here
    !! See comment on __UGLY_BUT_FAST in libint_wrapper

#if defined (__XT3)
#if defined (__UGLY_BUT_FAST)
    !! do_nothing
#else
    ALLOCATE(p_work(mysize))
#endif
#endif

#if defined (__PGI)
#if defined (__UGLY_BUT_FAST)
    !! do_nothing
#else
    ALLOCATE(p_work(mysize))
#endif
#endif
    
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
          CALL build_quartet_data(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                            potential_parameter, prim)
          lib%AB=A-B
          lib%CD=C-D
          CALL get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(2)
          CALL build_quartet_data(B,A,C,D,Zeta_B, Zeta_A, Zeta_C, Zeta_D, m_max,&
                            potential_parameter, prim)
          lib%AB=B-A
          lib%CD=C-D
          CALL get_eris(n_d, n_c, n_a, n_b, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(3)
          CALL build_quartet_data(A,B,D,C,Zeta_A, Zeta_B, Zeta_D, Zeta_C, m_max,&
                            potential_parameter, prim)
          lib%AB=A-B
          lib%CD=D-C
          CALL get_eris(n_c, n_d, n_b, n_a, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(4)
          CALL build_quartet_data(B,A,D,C,Zeta_B, Zeta_A, Zeta_D, Zeta_C, m_max,&
                            potential_parameter, prim)
          lib%AB=B-A
          lib%CD=D-C
          CALL get_eris(n_c, n_d, n_a, n_b, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(5)
          CALL build_quartet_data(C,D,A,B,Zeta_C, Zeta_D, Zeta_A, Zeta_B, m_max,&
                            potential_parameter, prim)
          lib%AB=C-D
          lib%CD=A-B
          CALL get_eris(n_b, n_a, n_d, n_c, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(6)
          CALL build_quartet_data(C,D,B,A,Zeta_C, Zeta_D, Zeta_B, Zeta_A, m_max,&
                            potential_parameter, prim)
          lib%AB=C-D
          lib%CD=B-A
          CALL get_eris(n_a, n_b, n_d, n_c, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(7)
          CALL build_quartet_data(D,C,A,B,Zeta_D, Zeta_C, Zeta_A, Zeta_B, m_max,&
                            potential_parameter, prim)
          lib%AB=D-C
          lib%CD=A-B
          CALL get_eris(n_b, n_a, n_c, n_d, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(8)
          CALL build_quartet_data(D,C,B,A,Zeta_D, Zeta_C, Zeta_B, Zeta_A, m_max,&
                            potential_parameter, prim)
          lib%AB=D-C
          lib%CD=B-A
          CALL get_eris(n_a, n_b, n_c, n_d, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
      END SELECT
    ELSE
      CALL build_quartet_data(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                              potential_parameter, prim)
      max_val = ABS(prim%F(1))
    END IF

    !! IF we use a Compiler that does not support whether ISO_C_BINDING nor
    !! CRAY POINTERS, we have to DEALLOCATE work array here
    !! See comment on __UGLY_BUT_FAST in libint_wrapper

#if defined (__XT3)
#if defined(__UGLY_BUT_FAST)
    !! do nothing
#else
    DEALLOCATE(p_work)
#endif
#endif

#if defined (__PGI)
#if defined(__UGLY_BUT_FAST)
    !! do nothing
#else
    DEALLOCATE(p_work)
#endif
#endif


  END SUBROUTINE evaluate_eri_screen
!***************************************************************************

!!***** libint_interface/pbc [1.0] *
!!
!!   NAME
!!     pbc
!!
!!   FUNCTION
!!     auxiliary function to apply minimum image convention
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     06.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE pbc(x, box_size, box_size_inv, retval)
    REAL(dp), INTENT(IN)                     :: x(3), box_size(3), box_size_inv(3)
    REAL(dp), INTENT(OUT)                    :: retval(3)

    IF( 2.0_dp*ABS(x(1))<box_size(1)) THEN
      retval(1) = x(1)
    ELSE
      retval(1)=x(1) - box_size(1)*ANINT(x(1)*box_size_inv(1))
    END IF
    IF( 2.0_dp*ABS(x(2))<box_size(2)) THEN
      retval(2) = x(2)
    ELSE
      retval(2)=x(2) - box_size(2)*ANINT(x(2)*box_size_inv(2))
    END IF
    IF( 2.0_dp*ABS(x(3))<box_size(3)) THEN
      retval(3) = x(3)
    ELSE
      retval(3)=x(3) - box_size(3)*ANINT(x(3)*box_size_inv(3))
    END IF
  END SUBROUTINE pbc
!***************************************************************************

END MODULE libint_interface
