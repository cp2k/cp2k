!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Integral GKS scheme
! *****************************************************************************
MODULE semi_empirical_int_gks

  USE dg_rho0_types,                   ONLY: dg_rho0_type
  USE dg_types,                        ONLY: dg_get,&
                                             dg_type
  USE erf_fn,                          ONLY: erfc
  USE f77_blas
  USE input_constants,                 ONLY: do_multipole_quadrupole
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi,&
                                             oorootpi
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_pool_types,                   ONLY: pw_pool_type
  USE semi_empirical_int_arrays,       ONLY: rij_threshold
  USE semi_empirical_mpole_types,      ONLY: semi_empirical_mpole_type
  USE semi_empirical_types,            ONLY: se_int_control_type,&
                                             se_taper_type,&
                                             semi_empirical_type,&
                                             setup_se_int_control_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_int_gks'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module=.FALSE.

  PUBLIC :: corecore_gks, rotnuc_gks, drotnuc_gks, rotint_gks, drotint_gks

CONTAINS

  SUBROUTINE rotnuc_gks (sepi,sepj,rij,e1b,e2a,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(45), INTENT(OUT), &
      OPTIONAL                               :: e1b, e2a
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER, PARAMETER                       :: Nb = 4, Nm = Nb+Nb*(Nb-1)/2
    INTEGER, DIMENSION(Nb, Nb), PARAMETER :: &
      M = RESHAPE( (/1,2,3,4,2,5,8,10,3,8,6,9,4,10,9,7/) , (/4,4/) )

    INTEGER                                  :: i, mu, nu
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(kind=dp), DIMENSION(Nm, Nm)         :: Coul

    rab= -rij

    IF(se_int_control%do_ewald_gks) THEN
       CALL makeCoulE(rab,sepi,sepj,Coul,se_int_control,error)
    ELSE
       CALL makeCoul(rab,sepi,sepj,Coul,se_int_control,error)
    END IF

    i = 0
    DO mu = 1, sepi%natorb
       DO nu = 1, mu
          i = i + 1
          e1b(i)= -Coul(M(mu,nu),1)*sepj%zeff
       END DO
    END DO

    i = 0
    DO mu = 1, sepj%natorb
       DO nu = 1, mu
          i = i + 1
          e2a(i)= -Coul(1,M(mu,nu))*sepi%zeff
       END DO
    END DO

  END SUBROUTINE rotnuc_gks

  SUBROUTINE rotint_gks (sepi,sepj,rij,w,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(2025), INTENT(OUT), &
      OPTIONAL                               :: w
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER, PARAMETER                       :: Nb = 4, Nm = Nb+Nb*(Nb-1)/2
    INTEGER, DIMENSION(Nb, Nb), PARAMETER :: &
      M = RESHAPE( (/1,2,3,4,2,5,8,10,3,8,6,9,4,10,9,7/) , (/4,4/) )

    INTEGER                                  :: i, lam, mu, nu, sig
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(kind=dp), DIMENSION(Nm, Nm)         :: Coul

    rab= -rij

    IF(se_int_control%do_ewald_gks) THEN
       CALL makeCoulE(rab,sepi,sepj,Coul,se_int_control,error)
    ELSE
       CALL makeCoul(rab,sepi,sepj,Coul,se_int_control,error)
    END IF

    i = 0
    DO mu = 1, sepi%natorb
       DO nu = 1, mu
          DO lam = 1, sepj%natorb
             DO sig = 1, lam
                i = i + 1
                w(i)= Coul(M(mu,nu),M(lam,sig))
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE rotint_gks

  SUBROUTINE drotnuc_gks(sepi,sepj,rij,de1b,de2a,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(3, 45), &
      INTENT(OUT), OPTIONAL                  :: de1b, de2a
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER, PARAMETER                       :: Nb = 4, Nm = Nb+Nb*(Nb-1)/2
    INTEGER, DIMENSION(Nb, Nb), PARAMETER :: &
      M = RESHAPE( (/1,2,3,4,2,5,8,10,3,8,6,9,4,10,9,7/) , (/4,4/) )

    INTEGER                                  :: i, mu, nu
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(kind=dp), DIMENSION(Nm, Nm, 3)      :: dCoul

    rab= -rij

    IF(se_int_control%do_ewald_gks) THEN
       CALL makedCoulE(rab,sepi,sepj,dCoul,se_int_control,error)
    ELSE
       CALL makedCoul(rab,sepi,sepj,dCoul,se_int_control,error)
    END IF

    i = 0
    DO mu = 1, sepi%natorb
       DO nu = 1, mu
          i = i + 1
          de1b(:,i)=  dCoul(M(mu,nu),1,:)*sepj%zeff
       END DO
    END DO

    i = 0
    DO mu = 1, sepj%natorb
       DO nu = 1, mu
          i = i + 1
          de2a(:,i)=  dCoul(1,M(mu,nu),:)*sepi%zeff
       END DO
    END DO

  END SUBROUTINE drotnuc_gks

  SUBROUTINE drotint_gks(sepi,sepj,rij,dw,se_int_control,se_taper, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(3, 2025), &
      INTENT(OUT)                            :: dw
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER, PARAMETER                       :: Nb = 4, Nm = Nb+Nb*(Nb-1)/2
    INTEGER, DIMENSION(Nb, Nb), PARAMETER :: &
      M = RESHAPE( (/1,2,3,4,2,5,8,10,3,8,6,9,4,10,9,7/) , (/4,4/) )

    INTEGER                                  :: i, lam, mu, nu, sig
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(kind=dp), DIMENSION(Nm, Nm, 3)      :: dCoul

    rab= -rij

    IF(se_int_control%do_ewald_gks) THEN
       CALL makedCoulE(rab,sepi,sepj,dCoul,se_int_control,error)
    ELSE
       CALL makedCoul(rab,sepi,sepj,dCoul,se_int_control,error)
    END IF

    i = 0
    DO mu = 1, sepi%natorb
       DO nu = 1, mu
          DO lam = 1, sepj%natorb
             DO sig = 1, lam
                i = i + 1
                dw(:,i)=  -dCoul(M(mu,nu),M(lam,sig),:)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE drotint_gks

  SUBROUTINE makeCoul(RAB,sepi,sepj,Coul,se_int_control,error)
    REAL(kind=dp), DIMENSION(3)              :: RAB
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    INTEGER, PARAMETER                       :: Nb = 4, Nm = Nb+Nb*(Nb-1)/2
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    REAL(kind=dp), DIMENSION(Nm, Nm)         :: Coul

    INTEGER, DIMENSION(Nb, Nb), PARAMETER :: &
      M = RESHAPE( (/1,2,3,4,2,5,8,10,3,8,6,9,4,10,9,7/) , (/4,4/) )

    INTEGER                                  :: i, imA, imB, j, jint, k1, k2, &
                                                k3, k4, size_1c_int
    LOGICAL                                  :: shortrange
    REAL(kind=dp) :: a2, ACOULA, ACOULB, d1, d1f(3), d2, d2f(3,3), d3, &
      d3f(3,3,3), d4, d4f(3,3,3,3), f, rr, tmp, w, w0, w1, w2, w3, w4, w5
    REAL(kind=dp), DIMENSION(3)              :: v
    REAL(kind=dp), DIMENSION(Nm)             :: M0A, M0B
    REAL(kind=dp), DIMENSION(Nm, 3)          :: M1A, M1B
    REAL(kind=dp), DIMENSION(Nm, 3, 3)       :: M2A, M2B
    TYPE(semi_empirical_mpole_type), POINTER :: mpole

    shortrange = se_int_control%shortrange

    Coul(:,:)=   0.0_dp

    size_1c_int = SIZE(sepi%w_mpole)
    DO jint = 1, size_1c_int
       mpole => sepi%w_mpole(jint)%mpole
       i=mpole%indi
       j=mpole%indj
       M0A(M(i,j))=     -mpole%c
       M1A(M(i,j),:)=   -mpole%d(:)
       M2A(M(i,j),:,:)= -mpole%qc(:,:)/3.0_dp
    ENDDO
    ACOULA= sepi%acoul

    size_1c_int = SIZE(sepj%w_mpole)
    DO jint = 1, size_1c_int
       mpole => sepj%w_mpole(jint)%mpole
       i=mpole%indi
       j=mpole%indj
       M0B(M(i,j))=     -mpole%c
       M1B(M(i,j),:)=   -mpole%d(:)
       M2B(M(i,j),:,:)= -mpole%qc(:,:)/3.0_dp
    ENDDO
    ACOULB= sepj%acoul

    v(:)=RAB(:)
    rr=SQRT(DOT_PRODUCT(v,v))

    a2=0.5_dp*(1.0_dp/ACOULA+1.0_dp/ACOULB)
    w0= a2*rr
    w=  EXP(-w0)
    w1= (1.0_dp+0.5_dp*w0)
    w2= (w1+0.5_dp*w0+0.5_dp*w0**2)
    w3= (w2+w0**3/6.0_dp)
    w4= (w3+w0**4/30.0_dp)
    w5= (w3+8.0_dp*w0**4/210.0_dp+w0**5/210.0_dp)

    IF(shortrange)THEN
       f=            (-w*w1)/rr
       d1=   -1.0_dp*(-w*w2)/rr**3
       d2=    3.0_dp*(-w*w3)/rr**5
       d3=  -15.0_dp*(-w*w4)/rr**7
       d4=  105.0_dp*(-w*w5)/rr**9
    ELSE
       f=            (1.0_dp-w*w1)/rr
       d1=   -1.0_dp*(1.0_dp-w*w2)/rr**3
       d2=    3.0_dp*(1.0_dp-w*w3)/rr**5
       d3=  -15.0_dp*(1.0_dp-w*w4)/rr**7
       d4=  105.0_dp*(1.0_dp-w*w5)/rr**9
    ENDIF

    d1f(:)=       0.0_dp
    d2f(:,:)=     0.0_dp
    d3f(:,:,:)=   0.0_dp
    d4f(:,:,:,:)= 0.0_dp

    DO k1=1,3
       d1f(k1)=d1f(k1)+v(k1)*d1
    ENDDO
    DO k1=1,3
       DO k2=1,3
          d2f(k1,k2)=d2f(k1,k2)+v(k1)*v(k2)*d2 
       ENDDO
       d2f(k1,k1)=d2f(k1,k1)+ d1 
    ENDDO
    DO k1=1,3
       DO k2=1,3
          DO k3=1,3
             d3f(k1,k2,k3)=d3f(k1,k2,k3)+v(k1)*v(k2)*v(k3)*d3 
          ENDDO
          w=v(k1)*d2
          d3f(k1,k2,k2)=d3f(k1,k2,k2)+w
          d3f(k2,k1,k2)=d3f(k2,k1,k2)+w
          d3f(k2,k2,k1)=d3f(k2,k2,k1)+w
       ENDDO
    ENDDO
    DO k1=1,3
       DO k2=1,3
          DO k3=1,3
             DO k4=1,3
                d4f(k1,k2,k3,k4)=d4f(k1,k2,k3,k4)+v(k1)*v(k2)*v(k3)*v(k4)*d4 
             ENDDO
             w=v(k1)*v(k2)*d3
             d4f(k1,k2,k3,k3)=d4f(k1,k2,k3,k3)+w
             d4f(k1,k3,k2,k3)=d4f(k1,k3,k2,k3)+w
             d4f(k3,k1,k2,k3)=d4f(k3,k1,k2,k3)+w
             d4f(k1,k3,k3,k2)=d4f(k1,k3,k3,k2)+w
             d4f(k3,k1,k3,k2)=d4f(k3,k1,k3,k2)+w
             d4f(k3,k3,k1,k2)=d4f(k3,k3,k1,k2)+w
          ENDDO
          d4f(k1,k1,k2,k2)=d4f(k1,k1,k2,k2)+d2
          d4f(k1,k2,k1,k2)=d4f(k1,k2,k1,k2)+d2
          d4f(k1,k2,k2,k1)=d4f(k1,k2,k2,k1)+d2
       ENDDO
    ENDDO

    DO imA=1,Nm
       DO imB=1,Nm

          w=   M0A(imA)*M0B(imB)*f
          DO k1=1,3
             w=w+(  M1A(imA,k1)*M0B(imB)-M0A(imA)*M1B(imB,k1) )*d1f(k1)
          ENDDO
          DO k1=1,3
             DO k2=1,3
                w=w+(  M2A(imA,k1,k2)*M0B(imB)-M1A(imA,k1)*M1B(imB,k2)+M0A(imA)*M2B(imB,k1,k2) )*d2f(k1,k2)
             ENDDO
          ENDDO
          DO k1=1,3
             DO k2=1,3
                DO k3=1,3
                   w=w+( -M2A(imA,k1,k2)*M1B(imB,k3)+M1A(imA,k1)*M2B(imB,k2,k3) )*d3f(k1,k2,k3)
                ENDDO
             ENDDO
          ENDDO

          DO k1=1,3
             DO k2=1,3
                DO k3=1,3
                   DO k4=1,3
                      w=w+ M2A(imA,k1,k2)*M2B(imB,k3,k4)*d4f(k1,k2,k3,k4)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

          IF( imA==1 .AND. imB==1 )THEN
             tmp=M2A(imA,1,1)*M2B(imB,1,1)  
             tmp=tmp*a2**4
             w=w -tmp*EXP(-w0)/rr
          ENDIF

          IF( imA==1 .AND. (imB>=5 .AND. imB<=7) )THEN
             tmp=M2A(imA,1,1)*M2B(imB,imB-4,imB-4)/3.0_dp
             tmp=tmp*a2**4
             w=w -tmp*EXP(-w0)/rr
          ENDIF

          IF( (imA>=5 .AND. imA<=7) .AND. imB==1 )THEN
             tmp=M2A(imA,imA-4,imA-4)*M2B(imB,1,1)/3.0_dp  
             tmp=tmp*a2**4
             w=w -tmp*EXP(-w0)/rr
          ENDIF

          IF( (imA>=5 .AND. imA<=7) .AND. (imB>=5 .AND. imB<=7) )THEN
             tmp=M2A(imA,imA-4,imA-4)*M2B(imB,imB-4,imB-4)/9.0_dp
             tmp=tmp*a2**4
             w=w -tmp*EXP(-w0)/rr
          ENDIF

          Coul(imA,imB)=Coul(imA,imB)+w

       ENDDO
    ENDDO

  END SUBROUTINE makeCoul

  SUBROUTINE makedCoul(RAB,sepi,sepj,dCoul,se_int_control,error)
    REAL(kind=dp), DIMENSION(3)              :: RAB
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    INTEGER, PARAMETER                       :: Nb = 4, Nm = Nb+Nb*(Nb-1)/2
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    REAL(kind=dp), DIMENSION(Nm, Nm, 3)      :: dCoul

    INTEGER, DIMENSION(Nb, Nb), PARAMETER :: &
      M = RESHAPE( (/1,2,3,4,2,5,8,10,3,8,6,9,4,10,9,7/) , (/4,4/) )

    INTEGER                                  :: i, imA, imB, j, jint, k1, k2, &
                                                k3, k4, k5, size_1c_int
    LOGICAL                                  :: shortrange
    REAL(kind=dp) :: a2, ACOULA, ACOULB, d1, d1f(3), d2, d2f(3,3), d3, &
      d3f(3,3,3), d4, d4f(3,3,3,3), d5, d5f(3,3,3,3,3), f, rr, tmp, w, w0, &
      w1, w2, w3, w4, w5, w6
    REAL(kind=dp), DIMENSION(3)              :: v
    REAL(kind=dp), DIMENSION(Nm)             :: M0A, M0B
    REAL(kind=dp), DIMENSION(Nm, 3)          :: M1A, M1B
    REAL(kind=dp), DIMENSION(Nm, 3, 3)       :: M2A, M2B
    TYPE(semi_empirical_mpole_type), POINTER :: mpole

    shortrange = se_int_control%shortrange

    dCoul(:,:,:)=   0.0_dp

    size_1c_int = SIZE(sepi%w_mpole)
    DO jint = 1, size_1c_int
       mpole => sepi%w_mpole(jint)%mpole
       i=mpole%indi
       j=mpole%indj
       M0A(M(i,j))=     -mpole%c
       M1A(M(i,j),:)=   -mpole%d(:)
       M2A(M(i,j),:,:)= -mpole%qc(:,:)/3.0_dp
    ENDDO
    ACOULA= sepi%acoul

    size_1c_int = SIZE(sepj%w_mpole)
    DO jint = 1, size_1c_int
       mpole => sepj%w_mpole(jint)%mpole
       i=mpole%indi
       j=mpole%indj
       M0B(M(i,j))=     -mpole%c
       M1B(M(i,j),:)=   -mpole%d(:)
       M2B(M(i,j),:,:)= -mpole%qc(:,:)/3.0_dp
    ENDDO
    ACOULB= sepj%acoul

    v(:)=RAB(:)
    rr=SQRT(DOT_PRODUCT(v,v))

    a2=0.5_dp*(1.0_dp/ACOULA+1.0_dp/ACOULB)
    w0= a2*rr
    w=  EXP(-w0)
    w1= (1.0_dp+0.5_dp*w0)
    w2= (w1+0.5_dp*w0+0.5_dp*w0**2)
    w3= (w2+w0**3/6.0_dp)
    w4= (w3+w0**4/30.0_dp)
    w5= (w3+4.0_dp*w0**4/105.0_dp+w0**5/210.0_dp)
    w6= (w3+15.0_dp*w0**4/378.0_dp+2.0_dp*w0**5/315.0_dp+w0**6/1890.0_dp)

    IF(shortrange)THEN
       f=            (-w*w1)/rr
       d1=   -1.0_dp*(-w*w2)/rr**3
       d2=    3.0_dp*(-w*w3)/rr**5
       d3=  -15.0_dp*(-w*w4)/rr**7
       d4=  105.0_dp*(-w*w5)/rr**9
       d5= -945.0_dp*(-w*w6)/rr**11
    ELSE
       f=            (1.0_dp-w*w1)/rr
       d1=   -1.0_dp*(1.0_dp-w*w2)/rr**3
       d2=    3.0_dp*(1.0_dp-w*w3)/rr**5
       d3=  -15.0_dp*(1.0_dp-w*w4)/rr**7
       d4=  105.0_dp*(1.0_dp-w*w5)/rr**9
       d5= -945.0_dp*(1.0_dp-w*w6)/rr**11
    ENDIF

    d1f(:)=         0.0_dp
    d2f(:,:)=       0.0_dp
    d3f(:,:,:)=     0.0_dp
    d4f(:,:,:,:)=   0.0_dp
    d5f(:,:,:,:,:)= 0.0_dp

    DO k1=1,3
       d1f(k1)=d1f(k1)+v(k1)*d1
    ENDDO
    DO k1=1,3
       DO k2=1,3
          d2f(k1,k2)=d2f(k1,k2)+v(k1)*v(k2)*d2 
       ENDDO
       d2f(k1,k1)=d2f(k1,k1)+ d1 
    ENDDO
    DO k1=1,3
       DO k2=1,3
          DO k3=1,3
             d3f(k1,k2,k3)=d3f(k1,k2,k3)+v(k1)*v(k2)*v(k3)*d3 
          ENDDO
          w=v(k1)*d2
          d3f(k1,k2,k2)=d3f(k1,k2,k2)+w
          d3f(k2,k1,k2)=d3f(k2,k1,k2)+w
          d3f(k2,k2,k1)=d3f(k2,k2,k1)+w
       ENDDO
    ENDDO
    DO k1=1,3
       DO k2=1,3
          DO k3=1,3
             DO k4=1,3
                d4f(k1,k2,k3,k4)=d4f(k1,k2,k3,k4)+v(k1)*v(k2)*v(k3)*v(k4)*d4 
             ENDDO
             w=v(k1)*v(k2)*d3
             d4f(k1,k2,k3,k3)=d4f(k1,k2,k3,k3)+w
             d4f(k1,k3,k2,k3)=d4f(k1,k3,k2,k3)+w
             d4f(k3,k1,k2,k3)=d4f(k3,k1,k2,k3)+w
             d4f(k1,k3,k3,k2)=d4f(k1,k3,k3,k2)+w
             d4f(k3,k1,k3,k2)=d4f(k3,k1,k3,k2)+w
             d4f(k3,k3,k1,k2)=d4f(k3,k3,k1,k2)+w
          ENDDO
          d4f(k1,k1,k2,k2)=d4f(k1,k1,k2,k2)+d2
          d4f(k1,k2,k1,k2)=d4f(k1,k2,k1,k2)+d2
          d4f(k1,k2,k2,k1)=d4f(k1,k2,k2,k1)+d2
       ENDDO
    ENDDO
    DO k1=1,3
       DO k2=1,3
          DO k3=1,3
             DO k4=1,3
                DO k5=1,3
                   d5f(k1,k2,k3,k4,k5)=d5f(k1,k2,k3,k4,k5)+v(k1)*v(k2)*v(k3)*v(k4)*v(k5)*d5 
                ENDDO
                w=v(k1)*v(k2)*v(k3)*d4
                d5f(k1,k2,k3,k4,k4)=d5f(k1,k2,k3,k4,k4)+w
                d5f(k1,k2,k4,k3,k4)=d5f(k1,k2,k4,k3,k4)+w
                d5f(k1,k4,k2,k3,k4)=d5f(k1,k4,k2,k3,k4)+w
                d5f(k4,k1,k2,k3,k4)=d5f(k4,k1,k2,k3,k4)+w
                d5f(k1,k2,k4,k4,k3)=d5f(k1,k2,k4,k4,k3)+w
                d5f(k1,k4,k2,k4,k3)=d5f(k1,k4,k2,k4,k3)+w
                d5f(k4,k1,k2,k4,k3)=d5f(k4,k1,k2,k4,k3)+w
                d5f(k1,k4,k4,k2,k3)=d5f(k1,k4,k4,k2,k3)+w
                d5f(k4,k1,k4,k2,k3)=d5f(k4,k1,k4,k2,k3)+w
                d5f(k4,k4,k1,k2,k3)=d5f(k4,k4,k1,k2,k3)+w
             ENDDO
             w=v(k1)*d3
             d5f(k1,k2,k2,k3,k3)=d5f(k1,k2,k2,k3,k3)+w
             d5f(k1,k2,k3,k2,k3)=d5f(k1,k2,k3,k2,k3)+w
             d5f(k1,k2,k3,k3,k2)=d5f(k1,k2,k3,k3,k2)+w
             d5f(k2,k1,k2,k3,k3)=d5f(k2,k1,k2,k3,k3)+w
             d5f(k2,k1,k3,k2,k3)=d5f(k2,k1,k3,k2,k3)+w
             d5f(k2,k1,k3,k3,k2)=d5f(k2,k1,k3,k3,k2)+w
             d5f(k2,k2,k1,k3,k3)=d5f(k2,k2,k1,k3,k3)+w
             d5f(k2,k3,k1,k2,k3)=d5f(k2,k3,k1,k2,k3)+w
             d5f(k2,k3,k1,k3,k2)=d5f(k2,k3,k1,k3,k2)+w
             d5f(k2,k2,k3,k1,k3)=d5f(k2,k2,k3,k1,k3)+w
             d5f(k2,k3,k2,k1,k3)=d5f(k2,k3,k2,k1,k3)+w
             d5f(k2,k3,k3,k1,k2)=d5f(k2,k3,k3,k1,k2)+w
             d5f(k2,k2,k3,k3,k1)=d5f(k2,k2,k3,k3,k1)+w
             d5f(k2,k3,k2,k3,k1)=d5f(k2,k3,k2,k3,k1)+w
             d5f(k2,k3,k3,k2,k1)=d5f(k2,k3,k3,k2,k1)+w
          ENDDO
       ENDDO
    ENDDO

    DO k5=1,3
       DO imA=1,Nm
          DO imB=1,Nm

             w=   M0A(imA)*M0B(imB)*d1f(k5)
             DO k1=1,3
                w=w+(  M1A(imA,k1)*M0B(imB)-M0A(imA)*M1B(imB,k1) )*d2f(k1,k5)
             ENDDO
             DO k1=1,3
                DO k2=1,3
                   w=w+(  M2A(imA,k1,k2)*M0B(imB)-M1A(imA,k1)*M1B(imB,k2)+M0A(imA)*M2B(imB,k1,k2) )*d3f(k1,k2,k5)
                ENDDO
             ENDDO
             DO k1=1,3
                DO k2=1,3
                   DO k3=1,3
                      w=w+( -M2A(imA,k1,k2)*M1B(imB,k3)+M1A(imA,k1)*M2B(imB,k2,k3) )*d4f(k1,k2,k3,k5)
                   ENDDO
                ENDDO
             ENDDO

             DO k1=1,3
                DO k2=1,3
                   DO k3=1,3
                      DO k4=1,3
                         w=w+ M2A(imA,k1,k2)*M2B(imB,k3,k4)*d5f(k1,k2,k3,k4,k5)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO

             IF( imA==1 .AND. imB==1 )THEN
                tmp=M2A(imA,1,1)*M2B(imB,1,1)  
                tmp=tmp*(1.0_dp+w0)*a2**4
                w=w +v(k5)*tmp*EXP(-w0)/rr**3
             ENDIF

             IF( imA==1 .AND. (imB>=5 .AND. imB<=7) )THEN
                tmp=M2A(imA,1,1)*M2B(imB,imB-4,imB-4)/3.0_dp
                tmp=tmp*(1.0_dp+w0)*a2**4
                w=w +v(k5)*tmp*EXP(-w0)/rr**3
             ENDIF

             IF( (imA>=5 .AND. imA<=7) .AND. imB==1 )THEN
                tmp=M2A(imA,imA-4,imA-4)*M2B(imB,1,1)/3.0_dp  
                tmp=tmp*(1.0_dp+w0)*a2**4
                w=w +v(k5)*tmp*EXP(-w0)/rr**3
             ENDIF

             IF( (imA>=5 .AND. imA<=7) .AND. (imB>=5 .AND. imB<=7) )THEN
                tmp=M2A(imA,imA-4,imA-4)*M2B(imB,imB-4,imB-4)/9.0_dp
                tmp=tmp*(1.0_dp+w0)*a2**4
                w=w +v(k5)*tmp*EXP(-w0)/rr**3
             ENDIF

             dCoul(imA,imB,k5)=dCoul(imA,imB,k5)+w

          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE makedCoul

  SUBROUTINE corecore_gks (sepi,sepj,rijv,enuc,denuc,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rijv
    REAL(dp), INTENT(OUT), OPTIONAL          :: enuc
    REAL(dp), DIMENSION(3), INTENT(OUT), &
      OPTIONAL                               :: denuc
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'corecore_gks', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: Nb = 4, Nm = Nb+Nb*(Nb-1)/2

    LOGICAL                                  :: failure, l_denuc, l_enuc
    REAL(dp)                                 :: alpi, alpj, dscale, rij, &
                                                scale, zz
    REAL(kind=dp), DIMENSION(Nm, Nm)         :: Coul, CoulE
    REAL(kind=dp), DIMENSION(Nm, Nm, 3)      :: dCoul, dCoulE
    TYPE(se_int_control_type)                :: se_int_control_off

    failure = .FALSE.

    rij=DOT_PRODUCT(rijv,rijv)

    l_enuc = PRESENT(enuc)
    l_denuc= PRESENT(denuc)
    IF ((rij > rij_threshold).AND.(l_enuc.OR.l_denuc)) THEN

       rij  = SQRT(rij)

       IF(se_int_control%do_ewald_gks) THEN
          CALL setup_se_int_control_type(se_int_control_off, shortrange=.FALSE.,&
               do_ewald_r3=.FALSE., integral_screening=se_int_control%integral_screening,&
               max_multipole=do_multipole_quadrupole)
          CALL makeCoul(rijv,sepi,sepj,Coul,se_int_control_off,error)
          IF (l_denuc) CALL makedCoul(rijv,sepi,sepj,dCoul,se_int_control_off,error)
          CALL makeCoulE(rijv,sepi,sepj,CoulE,se_int_control,error)
          IF (l_denuc) CALL makedCoulE(rijv,sepi,sepj,dCoulE,se_int_control,error)
       ELSE
          CALL makeCoul(rijv,sepi,sepj,Coul,se_int_control,error)
          CoulE = Coul
          IF (l_denuc) CALL makedCoul(rijv,sepi,sepj,dCoul,se_int_control,error)
          IF (l_denuc) dCoulE = dCoul
       END IF

       scale = 0.0_dp
       dscale= 0.0_dp
       zz = sepi%zeff*sepj%zeff
       alpi = sepi%alp
       alpj = sepj%alp
       scale = EXP(-alpi*rij)+EXP(-alpj*rij)
       IF (l_enuc) THEN
          enuc=zz*CoulE(1,1)+scale*zz*Coul(1,1)
       END IF
       IF (l_denuc) THEN
          dscale= -alpi*EXP(-alpi*rij)-alpj*EXP(-alpj*rij)
          denuc(:)= zz*dCoulE(1,1,:)+dscale*(rijv(:)/rij)*zz*Coul(1,1)+scale*zz*dCoul(1,1,:)
       END IF

    ENDIF
  END SUBROUTINE corecore_gks

  SUBROUTINE makeCoulE(RAB,sepi,sepj,Coul,se_int_control,error)
    REAL(KIND=dp), DIMENSION(3)              :: RAB
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    INTEGER, PARAMETER                       :: Nb = 4, Nm = Nb+Nb*(Nb-1)/2
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    REAL(KIND=dp), DIMENSION(Nm, Nm)         :: Coul

    INTEGER, DIMENSION(Nb, Nb), PARAMETER :: &
      M = RESHAPE( (/1,2,3,4,2,5,8,10,3,8,6,9,4,10,9,7/) , (/4,4/) )

    INTEGER                                  :: gpt, i, imA, imB, j, jint, &
                                                k1, k2, k3, k4, lp, mp, np, &
                                                size_1c_int
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    REAL(KIND=dp) :: a2, ACOULA, ACOULB, alpha, cc, d1, d1f(3), d2, d2f(3,3), &
      d3, d3f(3,3,3), d4, d4f(3,3,3,3), f, ff, kr, kr2, r1, r2, r3, r5, r7, &
      r9, rr, ss, w, w0, w1, w2, w3, w4, w5
    REAL(KIND=dp), DIMENSION(3)              :: kk, v
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: rho0
    REAL(KIND=dp), DIMENSION(Nm)             :: M0A, M0B
    REAL(KIND=dp), DIMENSION(Nm, 3)          :: M1A, M1B
    REAL(KIND=dp), DIMENSION(Nm, 3, 3)       :: M2A, M2B
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(dg_type), POINTER                   :: dg
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(semi_empirical_mpole_type), POINTER :: mpole

    alpha = se_int_control%ewald_gks%alpha
    pw_pool => se_int_control%ewald_gks%pw_pool
    dg => se_int_control%ewald_gks%dg
    CALL dg_get (dg, dg_rho0=dg_rho0)
    rho0    => dg_rho0%density%pw%cr3d
    pw_grid => pw_pool%pw_grid
    bds     => pw_grid%bounds

    Coul(:,:)=   0.0_dp

    size_1c_int = SIZE(sepi%w_mpole)
    DO jint = 1, size_1c_int
       mpole => sepi%w_mpole(jint)%mpole
       i=mpole%indi
       j=mpole%indj
       M0A(M(i,j))=     -mpole%c
       M1A(M(i,j),:)=   -mpole%d(:)
       M2A(M(i,j),:,:)= -mpole%qc(:,:)/3.0_dp
    ENDDO
    ACOULA= sepi%acoul

    size_1c_int = SIZE(sepj%w_mpole)
    DO jint = 1, size_1c_int
       mpole => sepj%w_mpole(jint)%mpole
       i=mpole%indi
       j=mpole%indj
       M0B(M(i,j))=     -mpole%c
       M1B(M(i,j),:)=   -mpole%d(:)
       M2B(M(i,j),:,:)= -mpole%qc(:,:)/3.0_dp
    ENDDO
    ACOULB= sepj%acoul

    v(:)=RAB(:)
    rr=SQRT(DOT_PRODUCT(v(:),v(:)))

    r1=1.0_dp/rr
    r2=r1*r1
    r3=r2*r1
    r5=r3*r2
    r7=r5*r2
    r9=r7*r2

    a2=0.5_dp*(1.0_dp/ACOULA+1.0_dp/ACOULB)

    w0= a2*rr
    w=  EXP(-w0)
    w1= (1.0_dp+0.5_dp*w0)
    w2= (w1+0.5_dp*w0+0.5_dp*w0**2)
    w3= (w2+w0**3/6.0_dp)
    w4= (w3+w0**4/30.0_dp)
    w5= (w3+8.0_dp*w0**4/210.0_dp+w0**5/210.0_dp)

    f=            (1.0_dp-w*w1)*r1
    d1=   -1.0_dp*(1.0_dp-w*w2)*r3
    d2=    3.0_dp*(1.0_dp-w*w3)*r5
    d3=  -15.0_dp*(1.0_dp-w*w4)*r7
    d4=  105.0_dp*(1.0_dp-w*w5)*r9


    kr=  alpha*rr
    kr2= kr*kr
    w0=  1.0_dp-erfc(kr)
    w1=  2.0_dp*oorootpi*EXP(-kr2)
    w2=  w1*kr


    f=f          -w0*r1
    d1=d1+  (-w2+w0)*r3
    d2=d2+  (w2*(3.0_dp+kr2*2.0_dp)-3.0_dp*w0)*r5
    d3=d3+  (-w2*(15.0_dp+kr2*(10.0_dp+kr2*4.0_dp))+15.0_dp*w0)*r7
    d4=d4+  (w2*(105.0_dp+kr2*(70.0_dp+kr2*(28.0_dp+kr2*8.0_dp)))-105.0_dp*w0)*r9


    d1f(:)=       0.0_dp
    d2f(:,:)=     0.0_dp
    d3f(:,:,:)=   0.0_dp
    d4f(:,:,:,:)= 0.0_dp

    DO k1=1,3
       d1f(k1)=d1f(k1)+v(k1)*d1
    ENDDO
    DO k1=1,3
       DO k2=1,3
          d2f(k1,k2)=d2f(k1,k2)+v(k1)*v(k2)*d2 
       ENDDO
       d2f(k1,k1)=d2f(k1,k1)+ d1 
    ENDDO
    DO k1=1,3
       DO k2=1,3
          DO k3=1,3
             d3f(k1,k2,k3)=d3f(k1,k2,k3)+v(k1)*v(k2)*v(k3)*d3 
          ENDDO
          w=v(k1)*d2
          d3f(k1,k2,k2)=d3f(k1,k2,k2)+w
          d3f(k2,k1,k2)=d3f(k2,k1,k2)+w
          d3f(k2,k2,k1)=d3f(k2,k2,k1)+w
       ENDDO
    ENDDO
    DO k1=1,3
       DO k2=1,3
          DO k3=1,3
             DO k4=1,3
                d4f(k1,k2,k3,k4)=d4f(k1,k2,k3,k4)+v(k1)*v(k2)*v(k3)*v(k4)*d4 
             ENDDO
             w=v(k1)*v(k2)*d3
             d4f(k1,k2,k3,k3)=d4f(k1,k2,k3,k3)+w
             d4f(k1,k3,k2,k3)=d4f(k1,k3,k2,k3)+w
             d4f(k3,k1,k2,k3)=d4f(k3,k1,k2,k3)+w
             d4f(k1,k3,k3,k2)=d4f(k1,k3,k3,k2)+w
             d4f(k3,k1,k3,k2)=d4f(k3,k1,k3,k2)+w
             d4f(k3,k3,k1,k2)=d4f(k3,k3,k1,k2)+w
          ENDDO
          d4f(k1,k1,k2,k2)=d4f(k1,k1,k2,k2)+d2
          d4f(k1,k2,k1,k2)=d4f(k1,k2,k1,k2)+d2
          d4f(k1,k2,k2,k1)=d4f(k1,k2,k2,k1)+d2
       ENDDO
    ENDDO

    DO imA=1,Nm
       DO imB=1,Nm

          w=   M0A(imA)*M0B(imB)*f
          DO k1=1,3
             w=w+(  M1A(imA,k1)*M0B(imB)-M0A(imA)*M1B(imB,k1) )*d1f(k1)
          ENDDO
          DO k1=1,3
             DO k2=1,3
                w=w+(  M2A(imA,k1,k2)*M0B(imB)-M1A(imA,k1)*M1B(imB,k2)+M0A(imA)*M2B(imB,k1,k2) )*d2f(k1,k2)
             ENDDO
          ENDDO
          DO k1=1,3
             DO k2=1,3
                DO k3=1,3
                   w=w+( -M2A(imA,k1,k2)*M1B(imB,k3)+M1A(imA,k1)*M2B(imB,k2,k3) )*d3f(k1,k2,k3)
                ENDDO
             ENDDO
          ENDDO

          DO k1=1,3
             DO k2=1,3
                DO k3=1,3
                   DO k4=1,3
                      w=w+ M2A(imA,k1,k2)*M2B(imB,k3,k4)*d4f(k1,k2,k3,k4)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

          Coul(imA,imB)=Coul(imA,imB)+w

       ENDDO
    ENDDO

    v(:)=RAB(:)

    f=            0.0_dp
    d1f(:)=       0.0_dp
    d2f(:,:)=     0.0_dp
    d3f(:,:,:)=   0.0_dp
    d4f(:,:,:,:)= 0.0_dp

    DO gpt = 1, pw_grid%ngpts_cut
       lp = pw_grid%mapl%pos(pw_grid%g_hat(1,gpt))
       mp = pw_grid%mapm%pos(pw_grid%g_hat(2,gpt))
       np = pw_grid%mapn%pos(pw_grid%g_hat(3,gpt))

       lp = lp + bds(1,1)
       mp = mp + bds(1,2)
       np = np + bds(1,3)

       IF (pw_grid%gsq(gpt) == 0.0_dp) CYCLE
       kk(:)=pw_grid%g(:,gpt)
       ff = 2.0_dp * fourpi * rho0(lp,mp,np)**2 * pw_grid%vol / pw_grid%gsq(gpt)

       kr=DOT_PRODUCT(kk(:),v(:))
       cc=COS(kr)
       ss=SIN(kr)

       f=f+cc*ff
       DO k1=1,3
          d1f(k1)=d1f(k1)-kk(k1)*ss*ff
       ENDDO
       DO k1=1,3
          DO k2=1,3
             d2f(k1,k2)=d2f(k1,k2)-kk(k1)*kk(k2)*cc*ff
          ENDDO
       ENDDO
       DO k1=1,3
          DO k2=1,3
             DO k3=1,3
                d3f(k1,k2,k3)=d3f(k1,k2,k3)+kk(k1)*kk(k2)*kk(k3)*ss*ff
             ENDDO
          ENDDO
       ENDDO
       DO k1=1,3
          DO k2=1,3
             DO k3=1,3
                DO k4=1,3
                   d4f(k1,k2,k3,k4)=d4f(k1,k2,k3,k4)+kk(k1)*kk(k2)*kk(k3)*kk(k4)*cc*ff
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    ENDDO

    DO imA=1,Nm
       DO imB=1,Nm

          w=   M0A(imA)*M0B(imB)*f
          DO k1=1,3
             w=w+(  M1A(imA,k1)*M0B(imB)-M0A(imA)*M1B(imB,k1) )*d1f(k1)
          ENDDO
          DO k1=1,3
             DO k2=1,3
                w=w+(  M2A(imA,k1,k2)*M0B(imB)-M1A(imA,k1)*M1B(imB,k2)+M0A(imA)*M2B(imB,k1,k2) )*d2f(k1,k2)
             ENDDO
          ENDDO
          DO k1=1,3
             DO k2=1,3
                DO k3=1,3
                   w=w+( -M2A(imA,k1,k2)*M1B(imB,k3)+M1A(imA,k1)*M2B(imB,k2,k3) )*d3f(k1,k2,k3)
                ENDDO
             ENDDO
          ENDDO

          DO k1=1,3
             DO k2=1,3
                DO k3=1,3
                   DO k4=1,3
                      w=w+ M2A(imA,k1,k2)*M2B(imB,k3,k4)*d4f(k1,k2,k3,k4)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

          Coul(imA,imB)=Coul(imA,imB)+w

       ENDDO
    ENDDO

    DO imA=1,Nm
       DO imB=1,Nm
          w= -M0A(imA)*M0B(imB)*0.25_dp*fourpi/(pw_grid%vol*alpha**2)
          Coul(imA,imB)=Coul(imA,imB)+w
       ENDDO
    ENDDO

  END SUBROUTINE makeCoulE

  SUBROUTINE makedCoulE(RAB,sepi,sepj,dCoul,se_int_control,error)
    REAL(KIND=dp), DIMENSION(3)              :: RAB
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    INTEGER, PARAMETER                       :: Nb = 4, Nm = Nb+Nb*(Nb-1)/2
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    REAL(KIND=dp), DIMENSION(Nm, Nm, 3)      :: dCoul

    INTEGER, DIMENSION(Nb, Nb), PARAMETER :: &
      M = RESHAPE( (/1,2,3,4,2,5,8,10,3,8,6,9,4,10,9,7/) , (/4,4/) )

    INTEGER                                  :: gpt, i, imA, imB, j, jint, &
                                                k1, k2, k3, k4, k5, lp, mp, &
                                                np, size_1c_int
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    REAL(KIND=dp) :: a2, ACOULA, ACOULB, alpha, cc, d1, d1f(3), d2, d2f(3,3), &
      d3, d3f(3,3,3), d4, d4f(3,3,3,3), d5, d5f(3,3,3,3,3), f, ff, kr, kr2, &
      r1, r11, r2, r3, r5, r7, r9, rr, ss, w, w0, w1, w2, w3, w4, w5, w6
    REAL(KIND=dp), DIMENSION(3)              :: kk, v
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: rho0
    REAL(kind=dp), DIMENSION(Nm)             :: M0A, M0B
    REAL(kind=dp), DIMENSION(Nm, 3)          :: M1A, M1B
    REAL(kind=dp), DIMENSION(Nm, 3, 3)       :: M2A, M2B
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(dg_type), POINTER                   :: dg
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(semi_empirical_mpole_type), POINTER :: mpole

    alpha = se_int_control%ewald_gks%alpha
    pw_pool => se_int_control%ewald_gks%pw_pool
    dg => se_int_control%ewald_gks%dg
    CALL dg_get (dg, dg_rho0=dg_rho0)
    rho0    => dg_rho0%density%pw%cr3d
    pw_grid => pw_pool%pw_grid
    bds     => pw_grid%bounds

    dCoul(:,:,:)=   0.0_dp

    size_1c_int = SIZE(sepi%w_mpole)
    DO jint = 1, size_1c_int
       mpole => sepi%w_mpole(jint)%mpole
       i=mpole%indi
       j=mpole%indj
       M0A(M(i,j))=     -mpole%c
       M1A(M(i,j),:)=   -mpole%d(:)
       M2A(M(i,j),:,:)= -mpole%qc(:,:)/3.0_dp
    ENDDO
    ACOULA= sepi%acoul

    size_1c_int = SIZE(sepj%w_mpole)
    DO jint = 1, size_1c_int
       mpole => sepj%w_mpole(jint)%mpole
       i=mpole%indi
       j=mpole%indj
       M0B(M(i,j))=     -mpole%c
       M1B(M(i,j),:)=   -mpole%d(:)
       M2B(M(i,j),:,:)= -mpole%qc(:,:)/3.0_dp
    ENDDO
    ACOULB= sepj%acoul

    v(:)=RAB(:)
    rr=SQRT(DOT_PRODUCT(v(:),v(:)))

    a2=0.5_dp*(1.0_dp/ACOULA+1.0_dp/ACOULB)

    r1=1.0_dp/rr
    r2=r1*r1
    r3=r2*r1
    r5=r3*r2
    r7=r5*r2
    r9=r7*r2
    r11=r9*r2

    w0= a2*rr
    w=  EXP(-w0)
    w1= (1.0_dp+0.5_dp*w0)
    w2= (w1+0.5_dp*w0+0.5_dp*w0**2)
    w3= (w2+w0**3/6.0_dp)
    w4= (w3+w0**4/30.0_dp)
    w5= (w3+8.0_dp*w0**4/210.0_dp+w0**5/210.0_dp)
    w6= (w3+5.0_dp*w0**4/126.0_dp+2.0_dp*w0**5/315.0_dp+w0**6/1890.0_dp)

    f=                (1.0_dp-w*w1)*r1
    d1=   -1.0_dp*(1.0_dp-w*w2)*r3
    d2=    3.0_dp*(1.0_dp-w*w3)*r5
    d3=  -15.0_dp*(1.0_dp-w*w4)*r7
    d4=  105.0_dp*(1.0_dp-w*w5)*r9
    d5= -945.0_dp*(1.0_dp-w*w6)*r11

    kr=  alpha*rr
    kr2= kr*kr
    w0=  1.0_dp-erfc(kr)
    w1=  2.0_dp*oorootpi*EXP(-kr2)
    w2=  w1*kr

    f=f          -w0*r1
    d1=d1+  (-w2+w0)*r3
    d2=d2+  (w2*(3.0_dp+kr2*2.0_dp)-3.0_dp*w0)*r5
    d3=d3+  (-w2*(15.0_dp+kr2*(10.0_dp+kr2*4.0_dp))+15.0_dp*w0)*r7
    d4=d4+  (w2*(105.0_dp+kr2*(70.0_dp+kr2*(28.0_dp+kr2*8.0_dp)))-105.0_dp*w0)*r9
    d5=d5+  (-w2*(945.0_dp+kr2*(630.0_dp+kr2*(252.0_dp+kr2*(72.0_dp+kr2*16.0_dp))))+945.0_dp*w0)*r11 

    d1f(:)=         0.0_dp
    d2f(:,:)=       0.0_dp
    d3f(:,:,:)=     0.0_dp
    d4f(:,:,:,:)=   0.0_dp
    d5f(:,:,:,:,:)= 0.0_dp

    DO k1=1,3
       d1f(k1)=d1f(k1)+v(k1)*d1
    ENDDO
    DO k1=1,3
       DO k2=1,3
          d2f(k1,k2)=d2f(k1,k2)+v(k1)*v(k2)*d2 
       ENDDO
       d2f(k1,k1)=d2f(k1,k1)+ d1 
    ENDDO
    DO k1=1,3
       DO k2=1,3
          DO k3=1,3
             d3f(k1,k2,k3)=d3f(k1,k2,k3)+v(k1)*v(k2)*v(k3)*d3 
          ENDDO
          w=v(k1)*d2
          d3f(k1,k2,k2)=d3f(k1,k2,k2)+w
          d3f(k2,k1,k2)=d3f(k2,k1,k2)+w
          d3f(k2,k2,k1)=d3f(k2,k2,k1)+w
       ENDDO
    ENDDO
    DO k1=1,3
       DO k2=1,3
          DO k3=1,3
             DO k4=1,3
                d4f(k1,k2,k3,k4)=d4f(k1,k2,k3,k4)+v(k1)*v(k2)*v(k3)*v(k4)*d4 
             ENDDO
             w=v(k1)*v(k2)*d3
             d4f(k1,k2,k3,k3)=d4f(k1,k2,k3,k3)+w
             d4f(k1,k3,k2,k3)=d4f(k1,k3,k2,k3)+w
             d4f(k3,k1,k2,k3)=d4f(k3,k1,k2,k3)+w
             d4f(k1,k3,k3,k2)=d4f(k1,k3,k3,k2)+w
             d4f(k3,k1,k3,k2)=d4f(k3,k1,k3,k2)+w
             d4f(k3,k3,k1,k2)=d4f(k3,k3,k1,k2)+w
          ENDDO
          d4f(k1,k1,k2,k2)=d4f(k1,k1,k2,k2)+d2
          d4f(k1,k2,k1,k2)=d4f(k1,k2,k1,k2)+d2
          d4f(k1,k2,k2,k1)=d4f(k1,k2,k2,k1)+d2
       ENDDO
    ENDDO
    DO k1=1,3
       DO k2=1,3
          DO k3=1,3
             DO k4=1,3
                DO k5=1,3
                   d5f(k1,k2,k3,k4,k5)=d5f(k1,k2,k3,k4,k5)+v(k1)*v(k2)*v(k3)*v(k4)*v(k5)*d5 
                ENDDO
                w=v(k1)*v(k2)*v(k3)*d4
                d5f(k1,k2,k3,k4,k4)=d5f(k1,k2,k3,k4,k4)+w
                d5f(k1,k2,k4,k3,k4)=d5f(k1,k2,k4,k3,k4)+w
                d5f(k1,k4,k2,k3,k4)=d5f(k1,k4,k2,k3,k4)+w
                d5f(k4,k1,k2,k3,k4)=d5f(k4,k1,k2,k3,k4)+w
                d5f(k1,k2,k4,k4,k3)=d5f(k1,k2,k4,k4,k3)+w
                d5f(k1,k4,k2,k4,k3)=d5f(k1,k4,k2,k4,k3)+w
                d5f(k4,k1,k2,k4,k3)=d5f(k4,k1,k2,k4,k3)+w
                d5f(k1,k4,k4,k2,k3)=d5f(k1,k4,k4,k2,k3)+w
                d5f(k4,k1,k4,k2,k3)=d5f(k4,k1,k4,k2,k3)+w
                d5f(k4,k4,k1,k2,k3)=d5f(k4,k4,k1,k2,k3)+w
             ENDDO
             w=v(k1)*d3
             d5f(k1,k2,k2,k3,k3)=d5f(k1,k2,k2,k3,k3)+w
             d5f(k1,k2,k3,k2,k3)=d5f(k1,k2,k3,k2,k3)+w
             d5f(k1,k2,k3,k3,k2)=d5f(k1,k2,k3,k3,k2)+w
             d5f(k2,k1,k2,k3,k3)=d5f(k2,k1,k2,k3,k3)+w
             d5f(k2,k1,k3,k2,k3)=d5f(k2,k1,k3,k2,k3)+w
             d5f(k2,k1,k3,k3,k2)=d5f(k2,k1,k3,k3,k2)+w
             d5f(k2,k2,k1,k3,k3)=d5f(k2,k2,k1,k3,k3)+w
             d5f(k2,k3,k1,k2,k3)=d5f(k2,k3,k1,k2,k3)+w
             d5f(k2,k3,k1,k3,k2)=d5f(k2,k3,k1,k3,k2)+w
             d5f(k2,k2,k3,k1,k3)=d5f(k2,k2,k3,k1,k3)+w
             d5f(k2,k3,k2,k1,k3)=d5f(k2,k3,k2,k1,k3)+w
             d5f(k2,k3,k3,k1,k2)=d5f(k2,k3,k3,k1,k2)+w
             d5f(k2,k2,k3,k3,k1)=d5f(k2,k2,k3,k3,k1)+w
             d5f(k2,k3,k2,k3,k1)=d5f(k2,k3,k2,k3,k1)+w
             d5f(k2,k3,k3,k2,k1)=d5f(k2,k3,k3,k2,k1)+w
          ENDDO
       ENDDO
    ENDDO

    DO k5=1,3
       DO imA=1,Nm
          DO imB=1,Nm

             w=   M0A(imA)*M0B(imB)*d1f(k5)
             DO k1=1,3
                w=w+(  M1A(imA,k1)*M0B(imB)-M0A(imA)*M1B(imB,k1) )*d2f(k1,k5)
             ENDDO
             DO k1=1,3
                DO k2=1,3
                   w=w+(  M2A(imA,k1,k2)*M0B(imB)-M1A(imA,k1)*M1B(imB,k2)+M0A(imA)*M2B(imB,k1,k2) )*d3f(k1,k2,k5)
                ENDDO
             ENDDO
             DO k1=1,3
                DO k2=1,3
                   DO k3=1,3
                      w=w+( -M2A(imA,k1,k2)*M1B(imB,k3)+M1A(imA,k1)*M2B(imB,k2,k3) )*d4f(k1,k2,k3,k5)
                   ENDDO
                ENDDO
             ENDDO

             DO k1=1,3
                DO k2=1,3
                   DO k3=1,3
                      DO k4=1,3
                         w=w+ M2A(imA,k1,k2)*M2B(imB,k3,k4)*d5f(k1,k2,k3,k4,k5)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO

             dCoul(imA,imB,k5)=dCoul(imA,imB,k5)+w

          ENDDO
       ENDDO
    ENDDO

    v(:)=RAB(:)

    f=            0.0_dp
    d1f(:)=       0.0_dp
    d2f(:,:)=     0.0_dp
    d3f(:,:,:)=   0.0_dp
    d4f(:,:,:,:)= 0.0_dp
    d5f(:,:,:,:,:)= 0.0_dp

    DO gpt = 1, pw_grid%ngpts_cut
       lp = pw_grid%mapl%pos(pw_grid%g_hat(1,gpt))
       mp = pw_grid%mapm%pos(pw_grid%g_hat(2,gpt))
       np = pw_grid%mapn%pos(pw_grid%g_hat(3,gpt))

       lp = lp + bds(1,1)
       mp = mp + bds(1,2)
       np = np + bds(1,3)

       IF (pw_grid%gsq(gpt) == 0.0_dp) CYCLE
       kk(:)=pw_grid%g(:,gpt)
       ff = 2.0_dp * fourpi * rho0(lp,mp,np)**2 * pw_grid%vol / pw_grid%gsq(gpt)

       kr=DOT_PRODUCT(kk(:),v(:))
       cc=COS(kr)
       ss=SIN(kr)

       f=f+cc*ff
       DO k1=1,3
          d1f(k1)=d1f(k1)-kk(k1)*ss*ff
       ENDDO
       DO k1=1,3
          DO k2=1,3
             d2f(k1,k2)=d2f(k1,k2)-kk(k1)*kk(k2)*cc*ff
          ENDDO
       ENDDO
       DO k1=1,3
          DO k2=1,3
             DO k3=1,3
                d3f(k1,k2,k3)=d3f(k1,k2,k3)+kk(k1)*kk(k2)*kk(k3)*ss*ff
             ENDDO
          ENDDO
       ENDDO
       DO k1=1,3
          DO k2=1,3
             DO k3=1,3
                DO k4=1,3
                   d4f(k1,k2,k3,k4)=d4f(k1,k2,k3,k4)+kk(k1)*kk(k2)*kk(k3)*kk(k4)*cc*ff
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       DO k1=1,3
          DO k2=1,3
             DO k3=1,3
                DO k4=1,3
                   DO k5=1,3
                      d5f(k1,k2,k3,k4,k5)=d5f(k1,k2,k3,k4,k5)-kk(k1)*kk(k2)*kk(k3)*kk(k4)*kk(k5)*ss*ff
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    ENDDO

    DO k5=1,3
       DO imA=1,Nm
          DO imB=1,Nm

             w=   M0A(imA)*M0B(imB)*d1f(k5)
             DO k1=1,3
                w=w+(  M1A(imA,k1)*M0B(imB)-M0A(imA)*M1B(imB,k1) )*d2f(k1,k5)
             ENDDO
             DO k1=1,3
                DO k2=1,3
                   w=w+(  M2A(imA,k1,k2)*M0B(imB)-M1A(imA,k1)*M1B(imB,k2)+M0A(imA)*M2B(imB,k1,k2) )*d3f(k1,k2,k5)
                ENDDO
             ENDDO
             DO k1=1,3
                DO k2=1,3
                   DO k3=1,3
                      w=w+( -M2A(imA,k1,k2)*M1B(imB,k3)+M1A(imA,k1)*M2B(imB,k2,k3) )*d4f(k1,k2,k3,k5)
                   ENDDO
                ENDDO
             ENDDO

             DO k1=1,3
                DO k2=1,3
                   DO k3=1,3
                      DO k4=1,3
                         w=w+ M2A(imA,k1,k2)*M2B(imB,k3,k4)*d5f(k1,k2,k3,k4,k5)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO

             dCoul(imA,imB,k5)=dCoul(imA,imB,k5)+w

          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE makedCoulE

END MODULE semi_empirical_int_gks
