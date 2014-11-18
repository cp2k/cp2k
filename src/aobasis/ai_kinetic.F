!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the kinetic energy integrals over Cartesian
!>      Gaussian-type functions.
!>
!>      [a|T|b] = [a|-nabla**2/2|b]
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      - Derivatives added (10.05.2002,MK)
!>      - Fully refactored (07.07.2014,JGH)
!> \author Matthias Krack (31.07.2000)
! *****************************************************************************
MODULE ai_kinetic
  USE ai_os_rr,                        ONLY: os_rr_ovlp
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_kinetic'

! *** Public subroutines ***

  PUBLIC :: kinetic

CONTAINS

! *****************************************************************************
!> \brief   Calculation of the two-center kinetic energy integrals [a|T|b] over
!>          Cartesian Gaussian-type functions. 
!> \param la_max Maximum L of basis on A
!> \param la_min Minimum L of basis on A
!> \param npgfa  Number of primitive functions in set of basis on A
!> \param rpgfa  Range of functions on A (used for prescreening)
!> \param zeta   Exponents of basis on center A
!> \param lb_max Maximum L of basis on A
!> \param lb_min Minimum L of basis on A
!> \param npgfb  Number of primitive functions in set of basis on B
!> \param rpgfb  Range of functions on B (used for prescreening)
!> \param zetb   Exponents of basis on center B
!> \param rab    Distance vector between centers A and B
!> \param kab    Kinetic energy integrals, optional
!> \param dab    First derivatives of Kinetic energy integrals, optional
!> \param error  CP2K error reporting
!> \date    07.07.2014
!> \author  JGH
! *****************************************************************************
  SUBROUTINE kinetic(la_max,la_min,npgfa,rpgfa,zeta,&
                     lb_max,lb_min,npgfb,rpgfb,zetb,&
                     rab,kab,dab,error)
    INTEGER, INTENT(IN)                      :: la_max, la_min, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa, zeta
    INTEGER, INTENT(IN)                      :: lb_max, lb_min, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb, zetb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: kab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT), OPTIONAL                :: dab
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kinetic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ax, ay, az, bx, by, bz, coa, &
                                                cob, ia, ib, idx, idy, idz, &
                                                ipgf, istat, jpgf, la, lb, &
                                                ldrr, lma, lmb, ma, mb, na, &
                                                nb, ofa, ofb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, b, dsx, dsy, dsz, dtx, &
                                                dty, dtz, f0, rab2, tab, xhi, &
                                                zet
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: rr, tt
    REAL(KIND=dp), DIMENSION(3)              :: rap, rbp

! Distance of the centers a and b

    rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
    tab = SQRT(rab2)

    ! Maximum l for auxiliary integrals
    IF (PRESENT(kab)) THEN
       lma = la_max + 1
       lmb = lb_max + 1
    END IF
    IF (PRESENT(dab)) THEN
       lma = la_max + 2
       lmb = lb_max + 1
       idx = coset(1,0,0) - coset(0,0,0)
       idy = coset(0,1,0) - coset(0,0,0)
       idz = coset(0,0,1) - coset(0,0,0)
    END IF
    ldrr = MAX(lma,lmb) + 1

    ! Allocate space for auxiliary integrals
    ALLOCATE(rr(0:ldrr-1,0:ldrr-1,3),tt(0:ldrr-1,0:ldrr-1,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! Number of integrals, check size of arrays
    ofa = ncoset(la_min-1)
    ofb = ncoset(lb_min-1)
    na = ncoset(la_max) - ofa
    nb = ncoset(lb_max) - ofb
    IF (PRESENT(kab)) THEN
       CPPrecondition((SIZE(kab,1) >= na*npgfa),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(kab,2) >= nb*npgfb),cp_failure_level,routineP,error,failure)
    END IF
    IF (PRESENT(dab)) THEN
       CPPrecondition((SIZE(dab,1) >= na*npgfa),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(dab,2) >= nb*npgfb),cp_failure_level,routineP,error,failure)
       CPPrecondition((SIZE(dab,3) >= 3),cp_failure_level,routineP,error,failure)
    END IF

    ! Loops over all pairs of primitive Gaussian-type functions
    ma = 0
    DO ipgf=1,npgfa
       mb = 0
       DO jpgf=1,npgfb
          ! Distance Screening
          IF(rpgfa(ipgf) + rpgfb(jpgf) < tab) THEN
             IF (PRESENT(kab)) kab(ma+1:ma+na,mb+1:mb+nb) = 0.0_dp
             IF (PRESENT(dab)) dab(ma+1:ma+na,mb+1:mb+nb,1:3) = 0.0_dp
             mb = mb + nb
             CYCLE
          ENDIF

          ! Calculate some prefactors
          a = zeta(ipgf)
          b = zetb(jpgf)
          zet =  a+b
          xhi =  a*b/zet
          rap =  b*rab/zet
          rbp = -a*rab/zet

          ! [s|s] integral
          f0 = 0.5_dp*(pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

          ! Calculate the reccurence relation, overlap
          CALL os_rr_ovlp(rap,lma,rbp,lmb,zet,ldrr,rr)

          ! kinetic energy auxiliary integrals, overlap of [da/dx|db/dx]
          DO la=0,lma-1
             DO lb=0,lmb-1
                tt(la,lb,1) = 4.0_dp*a*b*rr(la+1,lb+1,1)
                tt(la,lb,2) = 4.0_dp*a*b*rr(la+1,lb+1,2)
                tt(la,lb,3) = 4.0_dp*a*b*rr(la+1,lb+1,3)
                IF(la > 0 .AND. lb > 0) THEN
                   tt(la,lb,1) = tt(la,lb,1) + REAL(la*lb,dp)*rr(la-1,lb-1,1)
                   tt(la,lb,2) = tt(la,lb,2) + REAL(la*lb,dp)*rr(la-1,lb-1,2)
                   tt(la,lb,3) = tt(la,lb,3) + REAL(la*lb,dp)*rr(la-1,lb-1,3)
                END IF
                IF(la > 0) THEN
                   tt(la,lb,1) = tt(la,lb,1) - 2.0_dp*REAL(la,dp)*b*rr(la-1,lb+1,1)
                   tt(la,lb,2) = tt(la,lb,2) - 2.0_dp*REAL(la,dp)*b*rr(la-1,lb+1,2)
                   tt(la,lb,3) = tt(la,lb,3) - 2.0_dp*REAL(la,dp)*b*rr(la-1,lb+1,3)
                END IF
                IF(lb > 0) THEN
                   tt(la,lb,1) = tt(la,lb,1) - 2.0_dp*REAL(lb,dp)*a*rr(la+1,lb-1,1)
                   tt(la,lb,2) = tt(la,lb,2) - 2.0_dp*REAL(lb,dp)*a*rr(la+1,lb-1,2)
                   tt(la,lb,3) = tt(la,lb,3) - 2.0_dp*REAL(lb,dp)*a*rr(la+1,lb-1,3)
                END IF
             END DO
          END DO

          DO lb = lb_min,lb_max
          DO bx = 0,lb
          DO by = 0,lb-bx
             bz = lb - bx - by
             cob = coset(bx,by,bz) - ofb
             ib = mb + cob
             DO la = la_min,la_max
             DO ax = 0,la
             DO ay = 0,la-ax
                az = la - ax - ay
                coa = coset(ax,ay,az) - ofa
                ia = ma + coa
                ! integrals
                IF (PRESENT(kab)) THEN
                   kab(ia,ib) = f0 * ( tt(ax,bx,1) * rr(ay,by,2) * rr(az,bz,3) +&
                      rr(ax,bx,1) * tt(ay,by,2) * rr(az,bz,3) + &
                      rr(ax,bx,1) * rr(ay,by,2) * tt(az,bz,3) )
                END IF
                ! first derivatives
                IF (PRESENT(dab)) THEN
                   ! dx
                   dsx = 2.0_dp * a * rr(ax+1,bx,1)
                   IF(ax > 0) dsx = dsx - REAL(ax,dp) * rr(ax-1,bx,1)
                   dtx = 2.0_dp * a * tt(ax+1,bx,1)
                   IF(ax > 0) dtx = dtx - REAL(ax,dp) * tt(ax-1,bx,1)
                   dab(ia,ib,idx) = dtx*rr(ay,by,2)*rr(az,bz,3) + &
                                    dsx*(tt(ay,by,2)*rr(az,bz,3)+rr(ay,by,2)*tt(az,bz,3))
                   ! dy
                   dsy = 2.0_dp * a * rr(ay+1,by,2)
                   IF(ay > 0) dsy = dsy - REAL(ay,dp) * rr(ay-1,by,2)
                   dty = 2.0_dp * a * tt(ay+1,by,2)
                   IF(ay > 0) dty = dty - REAL(ay,dp) * tt(ay-1,by,2)
                   dab(ia,ib,idy) = dty*rr(ax,bx,1)*rr(az,bz,3) + &
                                    dsy*(tt(ax,bx,1)*rr(az,bz,3)+rr(ax,bx,1)*tt(az,bz,3))
                   ! dz
                   dsz = 2.0_dp * a * rr(az+1,bz,3)
                   IF(az > 0) dsz = dsz - REAL(az,dp) * rr(az-1,bz,3)
                   dtz = 2.0_dp * a * tt(az+1,bz,3)
                   IF(az > 0) dtz = dtz - REAL(az,dp) * tt(az-1,bz,3)
                   dab(ia,ib,idz) = dtz*rr(ax,bx,1)*rr(ay,by,2) + &
                                    dsz*(tt(ax,bx,1)*rr(ay,by,2)+rr(ax,bx,1)*tt(ay,by,2))
                   ! scale
                   dab(ia,ib,1:3) = f0 * dab(ia,ib,1:3)
                END IF
                !
             ENDDO
             ENDDO
             ENDDO !la
          ENDDO
          ENDDO
          ENDDO !lb

          mb = mb + nb
       END DO
       ma = ma + na
    END DO

    DEALLOCATE(rr,tt,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE kinetic

END MODULE ai_kinetic
