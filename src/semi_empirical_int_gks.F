!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Integral GKS scheme: The order of the integrals in makeCoul reflects
!>        the standard order by MOPAC
!> \par History
!>      Teodoro Laino [tlaino] - 04.2009 : Adapted size arrays to d-orbitals and
!>                               get rid of the alternative ordering. Using the
!>                               CP2K one.
!>      Teodoro Laino [tlaino] - 04.2009 : Skip nullification (speed-up)
!>      Teodoro Laino [tlaino] - 04.2009 : Speed-up due to fortran arrays order
!>                               optimization and collection of common pieces of
!>                               code
! *****************************************************************************
MODULE semi_empirical_int_gks

  USE dg_rho0_types,                   ONLY: dg_rho0_type
  USE dg_types,                        ONLY: dg_get,&
                                             dg_type
  USE erf_fn,                          ONLY: erfc
  USE f77_blas
  USE input_constants,                 ONLY: do_multipole_none
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi,&
                                             oorootpi
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_pool_types,                   ONLY: pw_pool_type
  USE semi_empirical_int_arrays,       ONLY: indexb,&
                                             rij_threshold
  USE semi_empirical_mpole_types,      ONLY: semi_empirical_mpole_type
  USE semi_empirical_types,            ONLY: se_int_control_type,&
                                             se_taper_type,&
                                             semi_empirical_type,&
                                             setup_se_int_control_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_int_gks'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module=.FALSE.

  PUBLIC :: corecore_gks, rotnuc_gks, drotnuc_gks, rotint_gks, drotint_gks

CONTAINS

! *****************************************************************************
!> \brief Computes the electron-nuclei integrals
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param e1b ...
!> \param e2a ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE rotnuc_gks (sepi,sepj,rij,e1b,e2a,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(45), INTENT(OUT), &
      OPTIONAL                               :: e1b, e2a
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotnuc_gks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, mu, nu
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(kind=dp), DIMENSION(45, 45)         :: Coul

    rab= -rij

    IF(se_int_control%do_ewald_gks) THEN
       IF ( DOT_PRODUCT(rij,rij) > rij_threshold) THEN
         CALL makeCoulE(rab,sepi,sepj,Coul,se_int_control,error)
       ELSE
         CALL makeCoulE0(sepi,Coul,se_int_control,error)
       END IF
    ELSE
       CALL makeCoul(rab,sepi,sepj,Coul,se_int_control,error)
    END IF

    i = 0
    DO mu = 1, sepi%natorb
       DO nu = 1, mu
          i = i + 1
          e1b(i)= -Coul(i,1)*sepj%zeff
       END DO
    END DO

    i = 0
    DO mu = 1, sepj%natorb
       DO nu = 1, mu
          i = i + 1
          e2a(i)= -Coul(1,i)*sepi%zeff
       END DO
    END DO

  END SUBROUTINE rotnuc_gks

! *****************************************************************************
!> \brief Computes the electron-electron integrals
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param w ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE rotint_gks (sepi,sepj,rij,w,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(2025), INTENT(OUT), &
      OPTIONAL                               :: w
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotint_gks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ind1, ind2, lam, mu, nu, &
                                                sig
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(kind=dp), DIMENSION(45, 45)         :: Coul

    rab= -rij

    IF(se_int_control%do_ewald_gks) THEN
       IF ( DOT_PRODUCT(rij,rij) > rij_threshold) THEN
         CALL makeCoulE(rab,sepi,sepj,Coul,se_int_control,error)
       ELSE
         CALL makeCoulE0(sepi,Coul,se_int_control,error)
       END IF
    ELSE
       CALL makeCoul(rab,sepi,sepj,Coul,se_int_control,error)
    END IF

    i    = 0
    ind1 = 0
    DO mu = 1, sepi%natorb
       DO nu = 1, mu
          ind1 = ind1 + 1
          ind2 = 0
          DO lam = 1, sepj%natorb
             DO sig = 1, lam
                i    = i + 1
                ind2 = ind2 + 1
                w(i) = Coul(ind1,ind2)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE rotint_gks

! *****************************************************************************
!> \brief Computes the derivatives of the electron-nuclei integrals
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param de1b ...
!> \param de2a ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE drotnuc_gks(sepi,sepj,rij,de1b,de2a,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(3, 45), &
      INTENT(OUT), OPTIONAL                  :: de1b, de2a
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'drotnuc_gks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, mu, nu
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(kind=dp), DIMENSION(3, 45, 45)      :: dCoul

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
          de1b(1,i)=  dCoul(1,i,1)*sepj%zeff
          de1b(2,i)=  dCoul(2,i,1)*sepj%zeff
          de1b(3,i)=  dCoul(3,i,1)*sepj%zeff
       END DO
    END DO

    i = 0
    DO mu = 1, sepj%natorb
       DO nu = 1, mu
          i = i + 1
          de2a(1,i)=  dCoul(1,1,i)*sepi%zeff
          de2a(2,i)=  dCoul(2,1,i)*sepi%zeff
          de2a(3,i)=  dCoul(3,1,i)*sepi%zeff
       END DO
    END DO

  END SUBROUTINE drotnuc_gks

! *****************************************************************************
!> \brief Computes the derivatives of the electron-electron integrals
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param dw ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE drotint_gks(sepi,sepj,rij,dw,se_int_control,se_taper, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(3, 2025), &
      INTENT(OUT)                            :: dw
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'drotint_gks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ind1, ind2, lam, mu, nu, &
                                                sig
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(kind=dp), DIMENSION(3, 45, 45)      :: dCoul

    rab= -rij

    IF(se_int_control%do_ewald_gks) THEN
       CALL makedCoulE(rab,sepi,sepj,dCoul,se_int_control,error)
    ELSE
       CALL makedCoul(rab,sepi,sepj,dCoul,se_int_control,error)
    END IF

    i    = 0
    ind1 = 0
    DO mu = 1, sepi%natorb
       DO nu = 1, mu
          ind1 = ind1 + 1
          ind2 = 0
          DO lam = 1, sepj%natorb
             DO sig = 1, lam
                i = i + 1
                ind2 = ind2 + 1
                dw(1,i)=  -dCoul(1,ind1,ind2)
                dw(2,i)=  -dCoul(2,ind1,ind2)
                dw(3,i)=  -dCoul(3,ind1,ind2)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE drotint_gks

! *****************************************************************************
!> \brief Computes the primitives of the integrals (non-periodic case)
!> \param RAB ...
!> \param sepi ...
!> \param sepj ...
!> \param Coul ...
!> \param se_int_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE makeCoul(RAB,sepi,sepj,Coul,se_int_control,error)
    REAL(kind=dp), DIMENSION(3)              :: RAB
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(kind=dp), DIMENSION(45, 45), &
      INTENT(OUT)                            :: Coul
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'makeCoul', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iA, iB, imA, imB, jA, jB, k1, &
                                                k2, k3, k4
    LOGICAL                                  :: shortrange
    REAL(kind=dp)                            :: a2, ACOULA, ACOULB, d1, &
                                                d1f(3), d2, d2f(3,3), d3, &
                                                d3f(3,3,3), d4, d4f(3,3,3,3), &
                                                f, rr, w, w0, w1, w2, w3, w4, &
                                                w5
    REAL(kind=dp), DIMENSION(3)              :: v
    REAL(kind=dp), DIMENSION(3, 3, 45)       :: M2A, M2B
    REAL(kind=dp), DIMENSION(3, 45)          :: M1A, M1B
    REAL(kind=dp), DIMENSION(45)             :: M0A, M0B

    shortrange = se_int_control%shortrange
    CALL get_se_slater_multipole(sepi, M0A, M1A, M2A, ACOULA)
    CALL get_se_slater_multipole(sepj, M0B, M1B, M2B, ACOULB)

    v(1)=RAB(1)
    v(2)=RAB(2)
    v(3)=RAB(3)
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

    CALL build_d_tensor_gks(d1f, d2f, d3f, d4f, v=v, d1=d1, d2=d2, d3=d3, d4=d4)

    imA = 0
    DO iA = 1, sepi%natorb
       DO jA = 1, iA
          imA = imA + 1

          imB = 0
          DO iB = 1, sepj%natorb
             DO jB = 1, iB
                imB = imB + 1

                w=   M0A(imA)*M0B(imB)*f
                DO k1=1,3
                   w=w+(  M1A(k1,imA)*M0B(imB)-M0A(imA)*M1B(k1,imB) )*d1f(k1)
                ENDDO
                DO k2=1,3
                   DO k1=1,3
                      w=w+(  M2A(k1,k2,imA)*M0B(imB)-M1A(k1,imA)*M1B(k2,imB)+M0A(imA)*M2B(k1,k2,imB) )*d2f(k1,k2)
                   ENDDO
                ENDDO
                DO k3=1,3
                   DO k2=1,3
                      DO k1=1,3
                         w=w+( -M2A(k1,k2,imA)*M1B(k3,imB)+M1A(k1,imA)*M2B(k2,k3,imB) )*d3f(k1,k2,k3)
                      ENDDO
                   ENDDO
                ENDDO

                DO k4=1,3
                   DO k3=1,3
                      DO k2=1,3
                         DO k1=1,3
                            w=w+ M2A(k1,k2,imA)*M2B(k3,k4,imB)*d4f(k1,k2,k3,k4)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO

                Coul(imA,imB)=w
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE makeCoul

! *****************************************************************************
!> \brief Computes the derivatives of the primitives of the integrals
!>        (non-periodic case)
!> \param RAB ...
!> \param sepi ...
!> \param sepj ...
!> \param dCoul ...
!> \param se_int_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE makedCoul(RAB,sepi,sepj,dCoul,se_int_control,error)
    REAL(kind=dp), DIMENSION(3)              :: RAB
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(kind=dp), DIMENSION(3, 45, 45), &
      INTENT(OUT)                            :: dCoul
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'makedCoul', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iA, iB, imA, imB, jA, jB, k1, &
                                                k2, k3, k4
    LOGICAL                                  :: shortrange
    REAL(kind=dp) :: a2, ACOULA, ACOULB, d1, d1f(3), d2, d2f(3,3), d3, &
      d3f(3,3,3), d4, d4f(3,3,3,3), d5, d5f(3,3,3,3,3), f, rr, tmp, w, w0, &
      w1, w2, w3, w4, w5, w6
    REAL(kind=dp), DIMENSION(3)              :: v, wv
    REAL(kind=dp), DIMENSION(3, 3, 45)       :: M2A, M2B
    REAL(kind=dp), DIMENSION(3, 45)          :: M1A, M1B
    REAL(kind=dp), DIMENSION(45)             :: M0A, M0B

    shortrange = se_int_control%shortrange
    CALL get_se_slater_multipole(sepi, M0A, M1A, M2A, ACOULA)
    CALL get_se_slater_multipole(sepj, M0B, M1B, M2B, ACOULB)

    v(1)=RAB(1)
    v(2)=RAB(2)
    v(3)=RAB(3)
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

    CALL build_d_tensor_gks(d1f, d2f, d3f, d4f, d5f, v, d1, d2, d3, d4, d5)

    imA = 0
    DO iA = 1, sepi%natorb
       DO jA = 1, iA
          imA = imA + 1

          imB = 0
          DO iB = 1, sepj%natorb
             DO jB = 1, iB
                imB = imB + 1

                tmp  = M0A(imA)*M0B(imB)
                wv(1)= tmp*d1f(1)
                wv(2)= tmp*d1f(2)
                wv(3)= tmp*d1f(3)
                DO k1=1,3
                   tmp = M1A(k1,imA)*M0B(imB)-M0A(imA)*M1B(k1,imB)
                   wv(1)=wv(1)+tmp*d2f(1,k1)
                   wv(2)=wv(2)+tmp*d2f(2,k1)
                   wv(3)=wv(3)+tmp*d2f(3,k1)
                ENDDO
                DO k2=1,3
                   DO k1=1,3
                      tmp = M2A(k1,k2,imA)*M0B(imB)-M1A(k1,imA)*M1B(k2,imB)+M0A(imA)*M2B(k1,k2,imB)
                      wv(1)=wv(1)+tmp*d3f(1,k1,k2)
                      wv(2)=wv(2)+tmp*d3f(2,k1,k2)
                      wv(3)=wv(3)+tmp*d3f(3,k1,k2)
                   ENDDO
                ENDDO
                DO k3=1,3
                   DO k2=1,3
                      DO k1=1,3
                         tmp = -M2A(k1,k2,imA)*M1B(k3,imB)+M1A(k1,imA)*M2B(k2,k3,imB)
                         wv(1)=wv(1)+tmp*d4f(1,k1,k2,k3)
                         wv(2)=wv(2)+tmp*d4f(2,k1,k2,k3)
                         wv(3)=wv(3)+tmp*d4f(3,k1,k2,k3)
                      ENDDO
                   ENDDO
                ENDDO

                DO k4=1,3
                   DO k3=1,3
                      DO k2=1,3
                         DO k1=1,3
                            tmp = M2A(k1,k2,imA)*M2B(k3,k4,imB)
                            wv(1)=wv(1)+ tmp*d5f(1,k1,k2,k3,k4)
                            wv(2)=wv(2)+ tmp*d5f(2,k1,k2,k3,k4)
                            wv(3)=wv(3)+ tmp*d5f(3,k1,k2,k3,k4)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO

                dCoul(1,imA,imB)=wv(1)
                dCoul(2,imA,imB)=wv(2)
                dCoul(3,imA,imB)=wv(3)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE makedCoul

! *****************************************************************************
!> \brief Computes nuclei-nuclei interactions
!> \param sepi ...
!> \param sepj ...
!> \param rijv ...
!> \param enuc ...
!> \param denuc ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
! *****************************************************************************
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

    LOGICAL                                  :: failure, l_denuc, l_enuc
    REAL(dp)                                 :: alpi, alpj, dscale, rij, &
                                                scale, zz
    REAL(kind=dp), DIMENSION(3, 45, 45)      :: dCoul, dCoulE
    REAL(kind=dp), DIMENSION(45, 45)         :: Coul, CoulE
    TYPE(se_int_control_type)                :: se_int_control_off

    failure = .FALSE.

    rij=DOT_PRODUCT(rijv,rijv)

    l_enuc = PRESENT(enuc)
    l_denuc= PRESENT(denuc)
    IF ((rij > rij_threshold).AND.(l_enuc.OR.l_denuc)) THEN

       rij  = SQRT(rij)

       IF(se_int_control%shortrange) THEN
          CALL setup_se_int_control_type(se_int_control_off, shortrange=.FALSE., do_ewald_r3=.FALSE.,&
               do_ewald_gks=.FALSE., integral_screening=se_int_control%integral_screening,&
               max_multipole=do_multipole_none, pc_coulomb_int=.FALSE.)
          CALL makeCoul(rijv,sepi,sepj,Coul,se_int_control_off,error)
          IF (l_denuc) CALL makedCoul(rijv,sepi,sepj,dCoul,se_int_control_off,error)
          IF (se_int_control%do_ewald_gks) THEN
             CALL makeCoulE(rijv,sepi,sepj,CoulE,se_int_control,error)
             IF (l_denuc) CALL makedCoulE(rijv,sepi,sepj,dCoulE,se_int_control,error)
          ELSE
             CALL makeCoul(rijv,sepi,sepj,CoulE,se_int_control,error)
             IF (l_denuc) CALL makedCoul(rijv,sepi,sepj,dCoulE,se_int_control,error)
          END IF
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
          denuc(1)= zz*dCoulE(1,1,1)+dscale*(rijv(1)/rij)*zz*Coul(1,1)+scale*zz*dCoul(1,1,1)
          denuc(2)= zz*dCoulE(2,1,1)+dscale*(rijv(2)/rij)*zz*Coul(1,1)+scale*zz*dCoul(2,1,1)
          denuc(3)= zz*dCoulE(3,1,1)+dscale*(rijv(3)/rij)*zz*Coul(1,1)+scale*zz*dCoul(3,1,1)
       END IF

    ELSE

       IF (se_int_control%do_ewald_gks) THEN
         zz = sepi%zeff*sepi%zeff
         CALL makeCoulE0(sepi,CoulE,se_int_control,error)
         IF (l_enuc) THEN
           enuc=zz*CoulE(1,1)
         END IF
         IF (l_denuc) THEN
           denuc=0._dp
         END IF
       END IF

    ENDIF
  END SUBROUTINE corecore_gks

! *****************************************************************************
!> \brief Computes the primitives of the integrals (periodic case)
!> \param RAB ...
!> \param sepi ...
!> \param sepj ...
!> \param Coul ...
!> \param se_int_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE makeCoulE(RAB,sepi,sepj,Coul,se_int_control,error)
    REAL(KIND=dp), DIMENSION(3)              :: RAB
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(OUT)                            :: Coul
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'makeCoulE', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: gpt, imA, imB, k1, k2, k3, &
                                                k4, lp, mp, np
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    REAL(KIND=dp) :: a2, ACOULA, ACOULB, alpha, cc, d1, d1f(3), d2, d2f(3,3), &
      d3, d3f(3,3,3), d4, d4f(3,3,3,3), f, ff, kr, kr2, r1, r2, r3, r5, r7, &
      r9, rr, ss, w, w0, w1, w2, w3, w4, w5
    REAL(KIND=dp), DIMENSION(3)              :: kk, v
    REAL(KIND=dp), DIMENSION(3, 3, 45)       :: M2A, M2B
    REAL(KIND=dp), DIMENSION(3, 45)          :: M1A, M1B
    REAL(KIND=dp), DIMENSION(45)             :: M0A, M0B
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: rho0
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(dg_type), POINTER                   :: dg
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_pool_type), POINTER              :: pw_pool

    alpha = se_int_control%ewald_gks%alpha
    pw_pool => se_int_control%ewald_gks%pw_pool
    dg => se_int_control%ewald_gks%dg
    CALL dg_get (dg, dg_rho0=dg_rho0)
    rho0    => dg_rho0%density%pw%cr3d
    pw_grid => pw_pool%pw_grid
    bds     => pw_grid%bounds

    CALL get_se_slater_multipole(sepi, M0A, M1A, M2A, ACOULA)
    CALL get_se_slater_multipole(sepj, M0B, M1B, M2B, ACOULB)

    v(1)=RAB(1)
    v(2)=RAB(2)
    v(3)=RAB(3)
    rr=SQRT(DOT_PRODUCT(v,v))

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

    CALL build_d_tensor_gks(d1f, d2f, d3f, d4f, v=v, d1=d1, d2=d2, d3=d3, d4=d4)

    DO imA=1, (sepi%natorb*(sepi%natorb+1))/2
       DO imB=1, (sepj%natorb*(sepj%natorb+1))/2

          w=   M0A(imA)*M0B(imB)*f
          DO k1=1,3
             w=w+(  M1A(k1,imA)*M0B(imB)-M0A(imA)*M1B(k1,imB) )*d1f(k1)
          ENDDO
          DO k2=1,3
             DO k1=1,3
                w=w+(  M2A(k1,k2,imA)*M0B(imB)-M1A(k1,imA)*M1B(k2,imB)+M0A(imA)*M2B(k1,k2,imB) )*d2f(k1,k2)
             ENDDO
          ENDDO
          DO k3=1,3
             DO k2=1,3
                DO k1=1,3
                   w=w+( -M2A(k1,k2,imA)*M1B(k3,imB)+M1A(k1,imA)*M2B(k2,k3,imB) )*d3f(k1,k2,k3)
                ENDDO
             ENDDO
          ENDDO

          DO k4=1,3
             DO k3=1,3
                DO k2=1,3
                   DO k1=1,3
                      w=w+ M2A(k1,k2,imA)*M2B(k3,k4,imB)*d4f(k1,k2,k3,k4)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

          Coul(imA,imB)=w
       ENDDO
    ENDDO

    v(1)=RAB(1)
    v(2)=RAB(2)
    v(3)=RAB(3)

    f=    0.0_dp
    d1f = 0.0_dp
    d2f = 0.0_dp
    d3f = 0.0_dp
    d4f = 0.0_dp

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

       kr=DOT_PRODUCT(kk,v)
       cc=COS(kr)
       ss=SIN(kr)

       f=f+cc*ff
       DO k1=1,3
          d1f(k1)=d1f(k1)-kk(k1)*ss*ff
       ENDDO
       DO k2=1,3
          DO k1=1,3
             d2f(k1,k2)=d2f(k1,k2)-kk(k1)*kk(k2)*cc*ff
          ENDDO
       ENDDO
       DO k3=1,3
          DO k2=1,3
             DO k1=1,3
                d3f(k1,k2,k3)=d3f(k1,k2,k3)+kk(k1)*kk(k2)*kk(k3)*ss*ff
             ENDDO
          ENDDO
       ENDDO
       DO k4=1,3
          DO k3=1,3
             DO k2=1,3
                DO k1=1,3
                   d4f(k1,k2,k3,k4)=d4f(k1,k2,k3,k4)+kk(k1)*kk(k2)*kk(k3)*kk(k4)*cc*ff
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    ENDDO

    DO imA=1, (sepi%natorb*(sepi%natorb+1))/2
       DO imB=1, (sepj%natorb*(sepj%natorb+1))/2

          w=   M0A(imA)*M0B(imB)*f
          DO k1=1,3
             w=w+(  M1A(k1,imA)*M0B(imB)-M0A(imA)*M1B(k1,imB) )*d1f(k1)
          ENDDO
          DO k2=1,3
             DO k1=1,3
                w=w+(  M2A(k1,k2,imA)*M0B(imB)-M1A(k1,imA)*M1B(k2,imB)+M0A(imA)*M2B(k1,k2,imB) )*d2f(k1,k2)
             ENDDO
          ENDDO
          DO k3=1,3
             DO k2=1,3
                DO k1=1,3
                   w=w+( -M2A(k1,k2,imA)*M1B(k3,imB)+M1A(k1,imA)*M2B(k2,k3,imB) )*d3f(k1,k2,k3)
                ENDDO
             ENDDO
          ENDDO

          DO k4=1,3
             DO k3=1,3
                DO k2=1,3
                   DO k1=1,3
                      w=w+ M2A(k1,k2,imA)*M2B(k3,k4,imB)*d4f(k1,k2,k3,k4)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

          Coul(imA,imB)=Coul(imA,imB)+w

       ENDDO
    ENDDO

    DO imA=1, (sepi%natorb*(sepi%natorb+1))/2
       DO imB=1, (sepj%natorb*(sepj%natorb+1))/2
          w= -M0A(imA)*M0B(imB)*0.25_dp*fourpi/(pw_grid%vol*alpha**2)
          Coul(imA,imB)=Coul(imA,imB)+w
       ENDDO
    ENDDO

  END SUBROUTINE makeCoulE

! *****************************************************************************
!> \brief Computes the derivatives of the primitives of the integrals
!>        (periodic case)
!> \param RAB ...
!> \param sepi ...
!> \param sepj ...
!> \param dCoul ...
!> \param se_int_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE makedCoulE(RAB,sepi,sepj,dCoul,se_int_control,error)
    REAL(KIND=dp), DIMENSION(3)              :: RAB
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3, 45, 45), &
      INTENT(OUT)                            :: dCoul
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'makedCoulE', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: gpt, imA, imB, k1, k2, k3, &
                                                k4, k5, lp, mp, np
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    REAL(KIND=dp) :: a2, ACOULA, ACOULB, alpha, cc, d1, d1f(3), d2, d2f(3,3), &
      d3, d3f(3,3,3), d4, d4f(3,3,3,3), d5, d5f(3,3,3,3,3), f, ff, kr, kr2, &
      r1, r11, r2, r3, r5, r7, r9, rr, ss, tmp, w, w0, w1, w2, w3, w4, w5, w6
    REAL(KIND=dp), DIMENSION(3)              :: kk, v, wv
    REAL(kind=dp), DIMENSION(3, 3, 45)       :: M2A, M2B
    REAL(kind=dp), DIMENSION(3, 45)          :: M1A, M1B
    REAL(kind=dp), DIMENSION(45)             :: M0A, M0B
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: rho0
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(dg_type), POINTER                   :: dg
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_pool_type), POINTER              :: pw_pool

    alpha = se_int_control%ewald_gks%alpha
    pw_pool => se_int_control%ewald_gks%pw_pool
    dg => se_int_control%ewald_gks%dg
    CALL dg_get (dg, dg_rho0=dg_rho0)
    rho0    => dg_rho0%density%pw%cr3d
    pw_grid => pw_pool%pw_grid
    bds     => pw_grid%bounds

    CALL get_se_slater_multipole(sepi, M0A, M1A, M2A, ACOULA)
    CALL get_se_slater_multipole(sepj, M0B, M1B, M2B, ACOULB)

    v(1)=RAB(1)
    v(2)=RAB(2)
    v(3)=RAB(3)
    rr=SQRT(DOT_PRODUCT(v,v))

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

    CALL build_d_tensor_gks(d1f, d2f, d3f, d4f, d5f, v, d1, d2, d3, d4, d5)

    DO imA=1, (sepi%natorb*(sepi%natorb+1))/2
       DO imB=1, (sepj%natorb*(sepj%natorb+1))/2

          tmp = M0A(imA)*M0B(imB)
          wv(1)=   tmp*d1f(1)
          wv(2)=   tmp*d1f(2)
          wv(3)=   tmp*d1f(3)

          DO k1=1,3
             tmp = M1A(k1,imA)*M0B(imB)-M0A(imA)*M1B(k1,imB)
             wv(1)=wv(1)+tmp*d2f(1,k1)
             wv(2)=wv(2)+tmp*d2f(2,k1)
             wv(3)=wv(3)+tmp*d2f(3,k1)
          ENDDO
          DO k2=1,3
             DO k1=1,3
                tmp = M2A(k1,k2,imA)*M0B(imB)-M1A(k1,imA)*M1B(k2,imB)+M0A(imA)*M2B(k1,k2,imB)
                wv(1)=wv(1)+tmp*d3f(1,k1,k2)
                wv(2)=wv(2)+tmp*d3f(2,k1,k2)
                wv(3)=wv(3)+tmp*d3f(3,k1,k2)
             ENDDO
          ENDDO
          DO k3=1,3
             DO k2=1,3
                DO k1=1,3
                   tmp = -M2A(k1,k2,imA)*M1B(k3,imB)+M1A(k1,imA)*M2B(k2,k3,imB)
                   wv(1)=wv(1)+tmp*d4f(1,k1,k2,k3)
                   wv(2)=wv(2)+tmp*d4f(2,k1,k2,k3)
                   wv(3)=wv(3)+tmp*d4f(3,k1,k2,k3)
                ENDDO
             ENDDO
          ENDDO

          DO k4=1,3
             DO k3=1,3
                DO k2=1,3
                   DO k1=1,3
                      tmp = M2A(k1,k2,imA)*M2B(k3,k4,imB)
                      wv(1)=wv(1)+ tmp*d5f(1,k1,k2,k3,k4)
                      wv(2)=wv(2)+ tmp*d5f(2,k1,k2,k3,k4)
                      wv(3)=wv(3)+ tmp*d5f(3,k1,k2,k3,k4)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

          dCoul(1,imA,imB)=wv(1)
          dCoul(2,imA,imB)=wv(2)
          dCoul(3,imA,imB)=wv(3)
       ENDDO
    ENDDO

    v(1)=RAB(1)
    v(2)=RAB(2)
    v(3)=RAB(3)

    f=    0.0_dp
    d1f = 0.0_dp
    d2f = 0.0_dp
    d3f = 0.0_dp
    d4f = 0.0_dp
    d5f = 0.0_dp

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

       kr=DOT_PRODUCT(kk,v)
       cc=COS(kr)
       ss=SIN(kr)

       f=f+cc*ff
       DO k1=1,3
          d1f(k1)=d1f(k1)-kk(k1)*ss*ff
       ENDDO
       DO k2=1,3
          DO k1=1,3
             d2f(k1,k2)=d2f(k1,k2)-kk(k1)*kk(k2)*cc*ff
          ENDDO
       ENDDO
       DO k3=1,3
          DO k2=1,3
             DO k1=1,3
                d3f(k1,k2,k3)=d3f(k1,k2,k3)+kk(k1)*kk(k2)*kk(k3)*ss*ff
             ENDDO
          ENDDO
       ENDDO
       DO k4=1,3
          DO k3=1,3
             DO k2=1,3
                DO k1=1,3
                   d4f(k1,k2,k3,k4)=d4f(k1,k2,k3,k4)+kk(k1)*kk(k2)*kk(k3)*kk(k4)*cc*ff
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       DO k5=1,3
          DO k4=1,3
             DO k3=1,3
                DO k2=1,3
                   DO k1=1,3
                      d5f(k1,k2,k3,k4,k5)=d5f(k1,k2,k3,k4,k5)-kk(k1)*kk(k2)*kk(k3)*kk(k4)*kk(k5)*ss*ff
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DO imA=1, (sepi%natorb*(sepi%natorb+1))/2
       DO imB=1, (sepj%natorb*(sepj%natorb+1))/2
          tmp = M0A(imA)*M0B(imB)
          wv(1)=   tmp*d1f(1)
          wv(2)=   tmp*d1f(2)
          wv(3)=   tmp*d1f(3)
          DO k1=1,3
             tmp = M1A(k1,imA)*M0B(imB)-M0A(imA)*M1B(k1,imB)
             wv(1)=wv(1)+tmp*d2f(1,k1)
             wv(2)=wv(2)+tmp*d2f(2,k1)
             wv(3)=wv(3)+tmp*d2f(3,k1)
          ENDDO
          DO k2=1,3
             DO k1=1,3
                tmp = M2A(k1,k2,imA)*M0B(imB)-M1A(k1,imA)*M1B(k2,imB)+M0A(imA)*M2B(k1,k2,imB)
                wv(1)=wv(1)+tmp*d3f(1,k1,k2)
                wv(2)=wv(2)+tmp*d3f(2,k1,k2)
                wv(3)=wv(3)+tmp*d3f(3,k1,k2)
             ENDDO
          ENDDO
          DO k3=1,3
             DO k2=1,3
                DO k1=1,3
                   tmp = -M2A(k1,k2,imA)*M1B(k3,imB)+M1A(k1,imA)*M2B(k2,k3,imB)
                   wv(1)=wv(1)+tmp*d4f(1,k1,k2,k3)
                   wv(2)=wv(2)+tmp*d4f(2,k1,k2,k3)
                   wv(3)=wv(3)+tmp*d4f(3,k1,k2,k3)
                ENDDO
             ENDDO
          ENDDO

          DO k4=1,3
             DO k3=1,3
                DO k2=1,3
                   DO k1=1,3
                      tmp =  M2A(k1,k2,imA)*M2B(k3,k4,imB)
                      wv(1)=wv(1)+tmp*d5f(1,k1,k2,k3,k4)
                      wv(2)=wv(2)+tmp*d5f(2,k1,k2,k3,k4)
                      wv(3)=wv(3)+tmp*d5f(3,k1,k2,k3,k4)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

          dCoul(1,imA,imB)=dCoul(1,imA,imB)+wv(1)
          dCoul(2,imA,imB)=dCoul(2,imA,imB)+wv(2)
          dCoul(3,imA,imB)=dCoul(3,imA,imB)+wv(3)
       ENDDO
    ENDDO

  END SUBROUTINE makedCoulE

! *****************************************************************************
!> \brief Builds the tensor for the evaluation of the integrals with the
!>        cartesian multipoles
!> \param d1f ...
!> \param d2f ...
!> \param d3f ...
!> \param d4f ...
!> \param d5f ...
!> \param v ...
!> \param d1 ...
!> \param d2 ...
!> \param d3 ...
!> \param d4 ...
!> \param d5 ...
! *****************************************************************************
  SUBROUTINE  build_d_tensor_gks(d1f, d2f, d3f, d4f, d5f, v, d1, d2, d3, d4, d5)
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: d1f
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT)                            :: d2f
    REAL(KIND=dp), DIMENSION(3, 3, 3), &
      INTENT(OUT)                            :: d3f
    REAL(KIND=dp), DIMENSION(3, 3, 3, 3), &
      INTENT(OUT)                            :: d4f
    REAL(KIND=dp), &
      DIMENSION(3, 3, 3, 3, 3), &
      INTENT(OUT), OPTIONAL                  :: d5f
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: v
    REAL(KIND=dp), INTENT(IN)                :: d1, d2, d3, d4
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: d5

    INTEGER                                  :: k1, k2, k3, k4, k5
    REAL(KIND=dp)                            :: w

    d1f = 0.0_dp
    d2f = 0.0_dp
    d3f = 0.0_dp
    d4f = 0.0_dp
    DO k1=1,3
       d1f(k1)=d1f(k1)+v(k1)*d1
    ENDDO
    DO k1=1,3
       DO k2=1,3
          d2f(k2,k1)=d2f(k2,k1)+v(k1)*v(k2)*d2
       ENDDO
       d2f(k1,k1)=d2f(k1,k1)+ d1
    ENDDO
    DO k1=1,3
       DO k2=1,3
          DO k3=1,3
             d3f(k3,k2,k1)=d3f(k3,k2,k1)+v(k1)*v(k2)*v(k3)*d3
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
                d4f(k4,k3,k2,k1)=d4f(k4,k3,k2,k1)+v(k1)*v(k2)*v(k3)*v(k4)*d4
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
    IF (PRESENT(d5f).AND.PRESENT(d5)) THEN
       d5f = 0.0_dp

       DO k1=1,3
          DO k2=1,3
             DO k3=1,3
                DO k4=1,3
                   DO k5=1,3
                      d5f(k5,k4,k3,k2,k1)=d5f(k5,k4,k3,k2,k1)+v(k1)*v(k2)*v(k3)*v(k4)*v(k5)*d5
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
    END IF
  END SUBROUTINE build_d_tensor_gks

! *****************************************************************************
!> \brief ...
!> \param sepi ...
!> \param Coul ...
!> \param se_int_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE makeCoulE0(sepi,Coul,se_int_control,error)
    TYPE(semi_empirical_type), POINTER       :: sepi
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(OUT)                            :: Coul
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'makeCoulE0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: gpt, imA, imB, k1, k2, k3, &
                                                k4, lp, mp, np
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    REAL(KIND=dp)                            :: alpha, d2f(3,3), &
                                                d4f(3,3,3,3), f, ff, w
    REAL(KIND=dp), DIMENSION(3)              :: kk
    REAL(KIND=dp), DIMENSION(3, 3, 45)       :: M2A
    REAL(KIND=dp), DIMENSION(3, 45)          :: M1A
    REAL(KIND=dp), DIMENSION(45)             :: M0A
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: rho0
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(dg_type), POINTER                   :: dg
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_pool_type), POINTER              :: pw_pool

    alpha = se_int_control%ewald_gks%alpha
    pw_pool => se_int_control%ewald_gks%pw_pool
    dg => se_int_control%ewald_gks%dg
    CALL dg_get (dg, dg_rho0=dg_rho0)
    rho0    => dg_rho0%density%pw%cr3d
    pw_grid => pw_pool%pw_grid
    bds     => pw_grid%bounds

    CALL get_se_slater_multipole(sepi, M0A, M1A, M2A)

    f=    0.0_dp
    d2f = 0.0_dp
    d4f = 0.0_dp

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

       f=f+ff
       DO k2=1,3
          DO k1=1,3
             d2f(k1,k2)=d2f(k1,k2)-kk(k1)*kk(k2)*ff
          ENDDO
       ENDDO
       DO k4=1,3
          DO k3=1,3
             DO k2=1,3
                DO k1=1,3
                   d4f(k1,k2,k3,k4)=d4f(k1,k2,k3,k4)+kk(k1)*kk(k2)*kk(k3)*kk(k4)*ff
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    ENDDO

    DO imA=1, (sepi%natorb*(sepi%natorb+1))/2
       DO imB=1, (sepi%natorb*(sepi%natorb+1))/2

          w=   M0A(imA)*M0A(imB)*f
          DO k2=1,3
             DO k1=1,3
                w=w+(  M2A(k1,k2,imA)*M0A(imB)-M1A(k1,imA)*M1A(k2,imB)+M0A(imA)*M2A(k1,k2,imB) )*d2f(k1,k2)
             ENDDO
          ENDDO

          DO k4=1,3
             DO k3=1,3
                DO k2=1,3
                   DO k1=1,3
                      w=w+ M2A(k1,k2,imA)*M2A(k3,k4,imB)*d4f(k1,k2,k3,k4)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

          Coul(imA,imB)=w

       ENDDO
    ENDDO

    DO imA=1, (sepi%natorb*(sepi%natorb+1))/2
       DO imB=1, (sepi%natorb*(sepi%natorb+1))/2
          w= -M0A(imA)*M0A(imB)*0.25_dp*fourpi/(pw_grid%vol*alpha**2)
          Coul(imA,imB)=Coul(imA,imB)+w
       ENDDO
    ENDDO

    DO imA=1, (sepi%natorb*(sepi%natorb+1))/2
       DO imB=1, (sepi%natorb*(sepi%natorb+1))/2

          w= M0A(imA)*M0A(imB)
          Coul(imA,imB)=Coul(imA,imB)-2.0_dp*alpha*oorootpi*w

          w= 0.0_dp
          DO k1=1,3
            w=w+ M1A(k1,imA)*M1A(k1,imB)
            w=w- M0A(imA)*M2A(k1,k1,imB)
            w=w- M2A(k1,k1,imA)*M0A(imB)
          ENDDO
          Coul(imA,imB)=Coul(imA,imB)-4.0_dp*alpha**3*oorootpi*w/3.0_dp

          w= 0.0_dp
          DO k2=1,3
             DO k1=1,3
                w=w+ 2.0_dp*M2A(k1,k2,imA)*M2A(k1,k2,imB)
                w=w+        M2A(k1,k1,imA)*M2A(k2,k2,imB)
             ENDDO
          ENDDO
          Coul(imA,imB)=Coul(imA,imB)-8.0_dp*alpha**5*oorootpi*w/5.0_dp
       ENDDO
    ENDDO
  END SUBROUTINE makeCoulE0

! *****************************************************************************
!> \brief Retrieves the multipole for the Slater integral evaluation
!> \param sepi ...
!> \param M0 ...
!> \param M1 ...
!> \param M2 ...
!> \param ACOUL ...
! *****************************************************************************
  SUBROUTINE get_se_slater_multipole(sepi, M0, M1, M2, ACOUL)
    TYPE(semi_empirical_type), POINTER       :: sepi
    REAL(kind=dp), DIMENSION(45), &
      INTENT(OUT)                            :: M0
    REAL(kind=dp), DIMENSION(3, 45), &
      INTENT(OUT)                            :: M1
    REAL(kind=dp), DIMENSION(3, 3, 45), &
      INTENT(OUT)                            :: M2
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: ACOUL

    INTEGER                                  :: i, j, jint, size_1c_int
    TYPE(semi_empirical_mpole_type), POINTER :: mpole

    NULLIFY(mpole)
    size_1c_int = SIZE(sepi%w_mpole)
    DO jint = 1, size_1c_int
       mpole => sepi%w_mpole(jint)%mpole
       i=mpole%indi
       j=mpole%indj
       M0(indexb(i,j))    = -mpole%cs

       M1(1,indexb(i,j))  = -mpole%ds(1)
       M1(2,indexb(i,j))  = -mpole%ds(2)
       M1(3,indexb(i,j))  = -mpole%ds(3)

       M2(1,1,indexb(i,j))= -mpole%qq(1,1)/3.0_dp
       M2(2,1,indexb(i,j))= -mpole%qq(2,1)/3.0_dp
       M2(3,1,indexb(i,j))= -mpole%qq(3,1)/3.0_dp

       M2(1,2,indexb(i,j))= -mpole%qq(1,2)/3.0_dp
       M2(2,2,indexb(i,j))= -mpole%qq(2,2)/3.0_dp
       M2(3,2,indexb(i,j))= -mpole%qq(3,2)/3.0_dp

       M2(1,3,indexb(i,j))= -mpole%qq(1,3)/3.0_dp
       M2(2,3,indexb(i,j))= -mpole%qq(2,3)/3.0_dp
       M2(3,3,indexb(i,j))= -mpole%qq(3,3)/3.0_dp
    ENDDO
    IF ( PRESENT (ACOUL) )  ACOUL = sepi%acoul
  END SUBROUTINE get_se_slater_multipole

END MODULE semi_empirical_int_gks
