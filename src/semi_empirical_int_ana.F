!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !

!-----------------------------------------------------------------------------!

!!****** cp2k/semi_empirical_int_ana [1.0] *
!!
!!   NAME
!!     semi_empirical_int_ana
!!
!!   FUNCTION
!!     Analytical derivatives of Integrals for semi-empirical methods
!!
!!   AUTHOR
!!     Teodoro Laino - Zurich University 04.2007 [tlaino]
!!
!!   MODIFICATION HISTORY
!!     23.11.2007 jhu   short range version of integrals     
!!     
!!***
!******************************************************************************

MODULE semi_empirical_int_ana

  USE kinds,                           ONLY: dp
  USE semi_empirical_integrals,        ONLY: &
       al, drotint, drotnuc, nucint, r0, rotint, rotnuc, taper, &
       taper_fn_init, terep
  USE semi_empirical_types,            ONLY: semi_empirical_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_int_ana'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module=.FALSE.
  PUBLIC :: rotnuc_ana, rotint_ana

  INTEGER, DIMENSION(84), PARAMETER ::      map_full = (/4,7,5,8,6,10,&
       14,17,15,18,16,20,24,27,25,28,26,30,31,61,32,62,33,63,34,67,35,68,36,70,&
       37,64,38,65,39,69,40,66,41,71,42,72,43,73,44,77,45,78,46,80,47,74,48,75,&
       49,79,50,76,51,91,52,92,53,93,54,97,55,98,56,100,57,94,58,95,59,99,60,96,&
       84,87,85,88,86,90/)

  ! *****************************************************************************

CONTAINS

!!****f* semi_empirical_int_ana/rotnuc_ana [1.0] *
!!
!!   NAME
!!     rotnuc_ana
!!
!!   FUNCTION
!!     Computes analytical gradients for semiempirical integrals
!!
!!   NOTES
!!     Analytical version of the MOPAC rotnuc routine
!!
!!   INPUTS
!!     on input  sepi   = Atomic parameters of first atom
!!               sepj   = Atomic parameters of second atom
!!               rijv   = Coordinate vector i -> j 
!!               itype  = 1  (MNDO)
!!                      = 2  (AM1)
!!                      = 3  (PM3)
!!                      = 4  (PDG)
!!   
!!   on output e1b,e2a= Array of electron-nuclear attraction integrals,
!!                      e1b = Electron on atom ni attracting nucleus of nj.
!!                      e2a = Electron on atom nj attracting nucleus of ni.
!!                      enuc= nuclear-nuclear repulsion term.
!!   
!!                     de1b = derivative of e1b term
!!                     de2a = derivative of e2a term
!!                     denuc= derivative of nuclear-nuclear repulsion term. 
!!
!!   AUTHOR
!!     Teodoro Laino - Zurich University
!!
!!   MODIFICATION HISTORY
!!     04.2007 created [tlaino]
!!
!!*** **********************************************************************

  RECURSIVE SUBROUTINE rotnuc_ana (sepi,sepj,rijv,itype,enuc,e1b,e2a,denuc,de1b,de2a,shortrange)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rijv
    INTEGER, INTENT(IN)                      :: itype
    REAL(dp), INTENT(OUT), OPTIONAL          :: enuc
    REAL(dp), DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: e1b, e2a, denuc
    REAL(dp), DIMENSION(:, :), INTENT(OUT), &
      OPTIONAL                               :: de1b, de2a
    LOGICAL, INTENT(IN)                      :: shortrange

    INTEGER                                  :: i, ig, j, nt
    LOGICAL                                  :: invert, l_de1b, l_de2a, &
                                                l_denuc, l_e1b, l_e2a, &
                                                l_enuc, lgrad, si, sj
    REAL(dp) :: a, alpi, alpj, apdg, ax, cppp1, cppp2, cpps1, cpps2, csp1, &
      csp2, css1, css2, dai, daj, dax, dbi, dbj, delta, denuc_loc, dqcorr, &
      dscale, dssss, dzz, enuc2, enuc_loc, pai, paj, pbi, pbj, qcorr, rij, &
      scale, ssss, xtmp, xx11, xx21, xx22, xx31, xx32, xx33, yy11, yy21, &
      yy22, yyzz11, yyzz21, yyzz22, zaf, zbf, zz, zz11, zz21, zz22, zz31, &
      zz32, zz33
    REAL(dp), DIMENSION(3) :: da, dcppp1, dcppp2, dcpps1, dcpps2, dcsp1, &
      dcsp2, dcss1, dcss2, drij, dx1, dx2, dx3, dxx11, dxx21, dxx22, dxx31, &
      dxx32, dxx33, dyy11, dyy21, dyy22, dyyzz11, dyyzz21, dyyzz22, dzz11, &
      dzz21, dzz22, dzz31, dzz32, dzz33, x, y, z
    REAL(dp), DIMENSION(3, 3)                :: dx, dy, dz
    REAL(dp), DIMENSION(4)                   :: fni1, fni2, fni3, fnj1, fnj2, &
                                                fnj3
    REAL(dp), DIMENSION(4, 2)                :: ccore, dccore
    REAL(KIND=dp), DIMENSION(10)             :: e1b2, e2a2
    REAL(KIND=dp), DIMENSION(10, 3)          :: de1b2, de2a2
    REAL(KIND=dp), DIMENSION(3)              :: denuc2

    x(:)=-rijv(:)
    rij=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
    ! Initialization
    l_enuc = PRESENT(enuc)
    l_e1b  = PRESENT(e1b)
    l_e2a  = PRESENT(e2a)
    l_de1b = PRESENT(de1b)
    l_de2a = PRESENT(de2a)
    l_denuc= PRESENT(denuc)
    lgrad  = l_de1b.OR.l_de2a
    ! Zeros all arrays
    IF (l_e1b) THEN
       e1b =0._dp
    END IF
    IF (l_e2a) THEN
       e2a =0._dp
    END IF
    IF (l_enuc) THEN
       enuc=0._dp
    END IF
    IF (l_de1b) THEN
       de1b =0._dp
    END IF
    IF (l_de2a) THEN
       de2a =0._dp
    END IF
    IF (l_denuc) THEN
       denuc=0._dp
    END IF
    IF (rij > 0.00002_dp) THEN
       ! Compute Integrals in diatomic frame opportunely inverted
       rij  = SQRT(rij)
       a    = 1._dp/rij
       IF (lgrad.OR.l_denuc) THEN
          CALL dnucint_ana(sepi,sepj,rij,ssss,ccore,dssss,dccore,itype,shortrange)
          ! Debug Procedure.. Check valifity of analytical gradients of nucint
          IF (debug_this_module) CALL check_dnucint_ana(sepi,sepj,rij,dssss,dccore,itype,shortrange)
       ELSE
          CALL nucint(sepi,sepj,rij,ssss,ccore,itype,shortrange)
       END IF
       x(1) = x(1)*a
       x(2) = x(2)*a
       x(3) = x(3)*a
       ! Possibly Invert Frame
       IF (ABS(x(3)) > 0.99999999_dp) THEN
          ! In order to avoid divergence just change Z axes into the Y axes
          ! all quantities are rotational/invertion invariant..
          invert = .TRUE.
          xtmp=x(3)
          x(3)=x(2)
          x(2)=xtmp
       ELSE
          invert = .FALSE.
       END IF
       IF (lgrad.OR.l_denuc) THEN
          drij(:) = -x(:)
          da  = -a**2*drij
          dx1 = -(/1.0_dp,0.0_dp,0.0_dp/)
          dx2 = -(/0.0_dp,1.0_dp,0.0_dp/)
          dx3 = -(/0.0_dp,0.0_dp,1.0_dp/)
          dx(1,:) = dx1*a+(x(1)/a)*da
          dx(2,:) = dx2*a+(x(2)/a)*da
          dx(3,:) = dx3*a+(x(3)/a)*da
       END IF
       z(3)=SQRT(1._dp-x(3)*x(3))
       a=1._dp/z(3)
       y(1)=-a*x(2)*SIGN(1._dp,x(1))
       y(2)=ABS(a*x(1))
       y(3)=0._dp
       z(1)=-a*x(1)*x(3)
       z(2)=-a*x(2)*x(3)
       ! Analytical Gradients
       IF (lgrad) THEN
          dz(3,:) = -a*x(3)*dx(3,:)
          da      = -a**2*dz(3,:)
          dy(1,:) = -da*x(2)*SIGN(1._dp,x(1))-a*dx(2,:)*SIGN(1._dp,x(1))
          dy(2,:) = SIGN(1._dp,a*x(1))*(da*x(1)+a*dx(1,:))
          dy(3,:) = 0.0_dp
          dz(1,:) =-da*x(1)*x(3)-a*dx(1,:)*x(3)-a*x(1)*dx(3,:)
          dz(2,:) =-da*x(2)*x(3)-a*dx(2,:)*x(3)-a*x(2)*dx(3,:)
       END IF
       si = (sepi%natorb > 1)
       sj = (sepj%natorb > 1)
       IF ( si .OR. sj ) THEN
          xx11 = x(1)*x(1)
          xx21 = x(2)*x(1)
          xx22 = x(2)*x(2)
          xx31 = x(3)*x(1)
          xx32 = x(3)*x(2)
          xx33 = x(3)*x(3)
          yy11 = y(1)*y(1)
          yy21 = y(2)*y(1)
          yy22 = y(2)*y(2)
          zz11 = z(1)*z(1)
          zz21 = z(2)*z(1)
          zz22 = z(2)*z(2)
          zz31 = z(3)*z(1)
          zz32 = z(3)*z(2)
          zz33 = z(3)*z(3)
          yyzz11 = yy11+zz11
          yyzz21 = yy21+zz21
          yyzz22 = yy22+zz22
          ! Analytical Gradients
          IF (lgrad) THEN
             dxx11 = dx(1,:)*x(1)+x(1)*dx(1,:)
             dxx21 = dx(2,:)*x(1)+x(2)*dx(1,:)
             dxx22 = dx(2,:)*x(2)+x(2)*dx(2,:)
             dxx31 = dx(3,:)*x(1)+x(3)*dx(1,:)
             dxx32 = dx(3,:)*x(2)+x(3)*dx(2,:)
             dxx33 = dx(3,:)*x(3)+x(3)*dx(3,:)
             dyy11 = dy(1,:)*y(1)+y(1)*dy(1,:)
             dyy21 = dy(2,:)*y(1)+y(2)*dy(1,:)
             dyy22 = dy(2,:)*y(2)+y(2)*dy(2,:)
             dzz11 = dz(1,:)*z(1)+z(1)*dz(1,:)
             dzz21 = dz(2,:)*z(1)+z(2)*dz(1,:)
             dzz22 = dz(2,:)*z(2)+z(2)*dz(2,:)
             dzz31 = dz(3,:)*z(1)+z(3)*dz(1,:)
             dzz32 = dz(3,:)*z(2)+z(3)*dz(2,:)
             dzz33 = dz(3,:)*z(3)+z(3)*dz(3,:)
             dyyzz11 = dyy11+dzz11
             dyyzz21 = dyy21+dzz21
             dyyzz22 = dyy22+dzz22
          END IF
       ENDIF
       IF (l_e1b.OR.l_de1b) THEN
          css1  = ccore(1,1)
          csp1  = ccore(2,1)
          cpps1 = ccore(3,1)
          cppp1 = ccore(4,1)
       END IF
       IF (l_e1b) THEN
          e1b(1)=-css1
          IF(sepi%natorb == 4) THEN
             e1b(2) = -csp1 *x(1)
             e1b(3) = -cpps1*xx11-cppp1*yyzz11
             e1b(4) = -csp1 *x(2)
             e1b(5) = -cpps1*xx21-cppp1*yyzz21
             e1b(6) = -cpps1*xx22-cppp1*yyzz22
             e1b(7) = -csp1 *x(3)
             e1b(8) = -cpps1*xx31-cppp1*zz31
             e1b(9) = -cpps1*xx32-cppp1*zz32
             e1b(10)= -cpps1*xx33-cppp1*zz33
          END IF
          IF (invert) CALL invert_integral(e1b,3)
       END IF
       IF (l_e2a.OR.l_de2a) THEN
          css2  = ccore(1,2)
          csp2  = ccore(2,2)
          cpps2 = ccore(3,2)
          cppp2 = ccore(4,2)
       END IF
       IF (l_e2a) THEN
          e2a(1)=-css2
          IF(sepj%natorb == 4) THEN
             e2a(2) = -csp2 *x(1)
             e2a(3) = -cpps2*xx11-cppp2*yyzz11
             e2a(4) = -csp2 *x(2)
             e2a(5) = -cpps2*xx21-cppp2*yyzz21
             e2a(6) = -cpps2*xx22-cppp2*yyzz22
             e2a(7) = -csp2 *x(3)
             e2a(8) = -cpps2*xx31-cppp2*zz31
             e2a(9) = -cpps2*xx32-cppp2*zz32
             e2a(10)= -cpps2*xx33-cppp2*zz33
          END IF
          IF (invert) CALL invert_integral(e2a,3)
       END IF
       ! Analytical Gradients
       IF (lgrad) THEN
          IF (l_de1b) THEN
             dcss1  = dccore(1,1)*drij
             dcsp1  = dccore(2,1)*drij
             dcpps1 = dccore(3,1)*drij
             dcppp1 = dccore(4,1)*drij
             de1b(1,:)=-dcss1
             IF(sepi%natorb == 4) THEN
                de1b(2,:) = -dcsp1 *x(1)-csp1 *dx(1,:)
                de1b(3,:) = -dcpps1*xx11-cpps1*dxx11-dcppp1*yyzz11-cppp1*dyyzz11
                de1b(4,:) = -dcsp1 *x(2)-csp1 *dx(2,:)
                de1b(5,:) = -dcpps1*xx21-cpps1*dxx21-dcppp1*yyzz21-cppp1*dyyzz21
                de1b(6,:) = -dcpps1*xx22-cpps1*dxx22-dcppp1*yyzz22-cppp1*dyyzz22
                de1b(7,:) = -dcsp1 *x(3)-csp1 *dx(3,:)
                de1b(8,:) = -dcpps1*xx31-cpps1*dxx31-dcppp1*zz31-cppp1*dzz31
                de1b(9,:) = -dcpps1*xx32-cpps1*dxx32-dcppp1*zz32-cppp1*dzz32
                de1b(10,:)= -dcpps1*xx33-cpps1*dxx33-dcppp1*zz33-cppp1*dzz33
             END IF
             IF (invert) CALL invert_derivative(de1b,3)
          END IF
          IF (l_de2a) THEN
             dcss2  = dccore(1,2)*drij
             dcsp2  = dccore(2,2)*drij
             dcpps2 = dccore(3,2)*drij
             dcppp2 = dccore(4,2)*drij
             de2a(1,:)=-dcss2
             IF(sepj%natorb == 4) THEN
                de2a(2,:) = -dcsp2 *x(1)-csp2 *dx(1,:)
                de2a(3,:) = -dcpps2*xx11-cpps2*dxx11-dcppp2*yyzz11-cppp2*dyyzz11
                de2a(4,:) = -dcsp2 *x(2)-csp2 *dx(2,:)
                de2a(5,:) = -dcpps2*xx21-cpps2*dxx21-dcppp2*yyzz21-cppp2*dyyzz21
                de2a(6,:) = -dcpps2*xx22-cpps2*dxx22-dcppp2*yyzz22-cppp2*dyyzz22
                de2a(7,:) = -dcsp2 *x(3)-csp2 *dx(3,:)
                de2a(8,:) = -dcpps2*xx31-cpps2*dxx31-dcppp2*zz31-cppp2*dzz31
                de2a(9,:) = -dcpps2*xx32-cpps2*dxx32-dcppp2*zz32-cppp2*dzz32
                de2a(10,:)= -dcpps2*xx33-cpps2*dxx33-dcppp2*zz33-cppp2*dzz33
             END IF
             IF (invert) CALL invert_derivative(de2a,3)
          END IF
       END IF
       ! ------------------------------------
       ! Specifically for nuclear interaction
       ! ------------------------------------ 
       IF (l_enuc.OR.l_denuc) THEN
          alpi = sepi%alp
          alpj = sepj%alp
          scale = EXP(-alpi*rij)+EXP(-alpj*rij)
          IF (l_denuc) THEN
             dscale= -alpi* EXP(-alpi*rij)-alpj*EXP(-alpj*rij)
          END IF
          nt=sepi%z+sepj%z
          IF(nt == 8 .OR. nt == 9) THEN
             IF(sepi%z == 7.OR.sepi%z == 8) THEN
                scale=scale+(0.5291772083_dp*rij-1._dp)*EXP(-alpi*rij)
                IF (l_denuc) THEN
                   dscale=dscale+0.5291772083_dp*EXP(-alpi*rij)-&
                          (0.5291772083_dp*rij-1._dp)*alpi*EXP(-alpi*rij)
                END IF
             END IF
             IF(sepj%z == 7.OR.sepj%z == 8) THEN
                scale=scale+(0.5291772083_dp*rij-1._dp)*EXP(-alpj*rij)
                IF (l_denuc) THEN
                   dscale=dscale+0.5291772083_dp*EXP(-alpj*rij)-&
                        (0.5291772083_dp*rij-1._dp)*alpj*EXP(-alpj*rij)
                END IF
             END IF
          ENDIF
          zz = sepi%zeff*sepj%zeff
          enuc_loc = zz*ssss
          IF (l_denuc) THEN
             denuc_loc= zz*dssss
             dscale=SIGN(1.0_dp,scale*enuc_loc)*(dscale*enuc_loc+scale*denuc_loc)
             dzz=-zz/rij**2
          END IF
          scale=ABS(scale*enuc_loc)
          zz=zz/rij
          IF(itype == 2 .OR. itype == 3 .OR. itype == 4) THEN
             IF(itype == 2 .AND. sepi%z == 5) THEN
                !special case AM1 Boron
                SELECT CASE (sepj%z)
                CASE DEFAULT
                   nt=1
                CASE (1)
                   nt=2
                CASE (6)
                   nt=3
                CASE (9,17,35,53)
                   nt=4
                END SELECT
                fni1(:) = sepi%bfn1(:,nt)
                fni2(:) = sepi%bfn2(:,nt)
                fni3(:) = sepi%bfn3(:,nt)
             ELSE
                fni1(:) = sepi%fn1(:)
                fni2(:) = sepi%fn2(:)
                fni3(:) = sepi%fn3(:)
             END IF
             IF(itype == 2 .AND. sepj%z == 5) THEN
                !special case AM1 Boron
                SELECT CASE (sepi%z)
                CASE DEFAULT
                   nt=1
                CASE (1)
                   nt=2
                CASE (6)
                   nt=3
                CASE (9,17,35,53)
                   nt=4
                END SELECT
                fnj1(:) = sepj%bfn1(:,nt)
                fnj2(:) = sepj%bfn2(:,nt)
                fnj3(:) = sepj%bfn3(:,nt)
             ELSE
                fnj1(:) = sepj%fn1(:)
                fnj2(:) = sepj%fn2(:)
                fnj3(:) = sepj%fn3(:)
             END IF
             ! AM1/PM3/PDG correction to nuclear repulsion
             DO ig=1,SIZE(fni1)
                IF(ABS(fni1(ig)) > 0._dp) THEN
                   ax = fni2(ig)*(rij-fni3(ig))**2
                   IF(ax <= 25._dp) THEN
                      scale=scale +zz*fni1(ig)*EXP(-ax)
                      IF (l_denuc) THEN
                         dax = fni2(ig)*2.0_dp*(rij-fni3(ig))
                         dscale = dscale + dzz*fni1(ig)*EXP(-ax)-dax*zz*fni1(ig)*EXP(-ax)
                      END IF
                   ENDIF
                ENDIF
                IF(ABS(fnj1(ig)) > 0._dp) THEN
                   ax = fnj2(ig)*(rij-fnj3(ig))**2
                   IF(ax <= 25._dp) THEN
                      scale=scale +zz*fnj1(ig)*EXP(-ax)
                      IF (l_denuc) THEN
                         dax = fnj2(ig)*2.0_dp*(rij-fnj3(ig))
                         dscale = dscale + dzz*fnj1(ig)*EXP(-ax)-dax*zz*fnj1(ig)*EXP(-ax)
                      END IF
                   ENDIF
                ENDIF
             END DO
          ENDIF
          IF(itype == 4) THEN
             ! PDDG function 
             zaf=sepi%zeff/nt
             zbf=sepj%zeff/nt
             pai=sepi%pre(1)
             pbi=sepi%pre(2)
             paj=sepj%pre(1)
             pbj=sepj%pre(2)
             dai=sepi%d(1)
             dbi=sepi%d(2)
             daj=sepj%d(1)
             dbj=sepj%d(2)
             apdg=10._dp*0.5291772083_dp*0.5291772083_dp
             qcorr = (zaf*pai+zbf*paj)*EXP(-apdg*(rij-dai-daj)**2)+ &
                     (zaf*pai+zbf*pbj)*EXP(-apdg*(rij-dai-dbj)**2)+ &
                     (zaf*pbi+zbf*paj)*EXP(-apdg*(rij-dbi-daj)**2)+ &
                     (zaf*pbi+zbf*pbj)*EXP(-apdg*(rij-dbi-dbj)**2)
             IF (l_denuc) THEN
                dqcorr = (zaf*pai+zbf*paj)*EXP(-apdg*(rij-dai-daj)**2)*(-2.0_dp*apdg*(rij-dai-daj))+ &
                         (zaf*pai+zbf*pbj)*EXP(-apdg*(rij-dai-dbj)**2)*(-2.0_dp*apdg*(rij-dai-dbj))+ &
                         (zaf*pbi+zbf*paj)*EXP(-apdg*(rij-dbi-daj)**2)*(-2.0_dp*apdg*(rij-dbi-daj))+ &
                         (zaf*pbi+zbf*pbj)*EXP(-apdg*(rij-dbi-dbj)**2)*(-2.0_dp*apdg*(rij-dbi-dbj))
             END IF
          ELSEIF (itype==0) THEN
             qcorr = 0.0_dp
             scale = 0.0_dp
             dscale= 0.0_dp
             dqcorr= 0.0_dp
          ELSE
             qcorr = 0.0_dp
             dqcorr= 0.0_dp
          END IF
          IF (l_enuc) THEN
             enuc=enuc_loc+scale+qcorr
          END IF
          IF (l_denuc) THEN
             denuc = (denuc_loc + dscale + dqcorr)*drij
             IF (invert) THEN
                xtmp=denuc(2)
                denuc(2)=denuc(3)
                denuc(3)=xtmp
             END IF
          END IF
       ENDIF
       !
       ! Debug for integral values and derivatives
       !
       IF (debug_this_module) THEN
          ! Check value of integrals
          CALL rotnuc (sepi,sepj,rijv,e1b2,e2a2,enuc2,itype,shortrange)
          IF (l_e1b) THEN
             DO J = 1, 10
                IF (ABS(e1b(j))>1.0E-6_dp) THEN
                   IF (ABS(e1b2(j)-e1b(j))>1.0E-7) THEN
                      WRITE(*,*)"e1b"
                      WRITE(*,'(10F12.6)')e1b
                      WRITE(*,'(10F12.6)')e1b2
                      WRITE(*,*)invert,SUM(ABS(e1b2-e1b)),e1b2-e1b
                      STOP
                   END IF
                END IF
             END DO
          END IF
          IF (l_e2a) THEN
             DO J = 1, 10
                IF (ABS(e2a(j))>1.0E-6_dp) THEN
                   IF (ABS(e2a2(j)-e2a(j))>1.0E-7) THEN
                      WRITE(*,*)"e2a"
                      WRITE(*,'(10F12.6)')e2a
                      WRITE(*,'(10F12.6)')e2a2
                      WRITE(*,*)invert,SUM(ABS(e2a2-e2a)),e2a2-e2a
                      STOP
                   END IF
                END IF
             END DO
          END IF
          IF (l_enuc) THEN
             IF (ABS(enuc-enuc2)>1.0E-7) THEN
                WRITE(*,*)invert,enuc-enuc2
             END IF
          END IF
          ! Check derivatives
          IF (lgrad.OR.l_denuc) THEN
             ! Numerical derivatives are obviosly a big problem..
             ! First of all let's decide if the value we get for delta is compatible
             ! with a reasonable value of the integral.. (compatible if the value of the
             ! integral is greater than 1.0E-6)
             delta = 1.0E-5_dp
             CALL drotnuc(sepi,sepj,rijv,de1b2,de2a2,denuc2,itype,delta=delta,shortrange=shortrange)
             CALL rotnuc_ana(sepi,sepj,rijv,itype,e1b=e1b2,e2a=e2a2,shortrange=shortrange)
             IF (l_de1b) THEN
                DO i = 1, 3
                   DO j = 1, 10
                      IF ((ABS(e1b2(j))>1.0E-6_dp).AND.(ABS(de1b(j,i))>delta*10)) THEN
                         IF (ABS((de1b2(j,i)-de1b(j,i))/de1b(j,i))*100.0_dp>1.0_dp) THEN
                            WRITE(*,*)"check de1b",i,j,ABS(de1b2(j,i)-de1b(j,i))/de1b(j,i)*100.0_dp
                            WRITE(*,'(10F12.6)')de1b
                            WRITE(*,*)
                            WRITE(*,'(10F12.6)')de1b2
                            WRITE(*,*)invert,SUM(ABS(de1b2-de1b))
                            WRITE(*,'(10F12.6)')de1b2-de1b
                            STOP
                         END IF
                      END IF
                   END DO
                END DO
             END IF
             IF (l_de2a) THEN
                DO i = 1, 3
                   DO j = 1, 10
                      IF ((ABS(e2a2(j))>1.0E-6_dp).AND.(ABS(de2a(j,i))>delta*10)) THEN
                         IF (ABS((de2a2(j,i)-de2a(j,i))/de2a(j,i))*100.0_dp>1.0_dp) THEN
                            WRITE(*,*)"check de2a",i,j,ABS(de2a2(j,i)-de2a(j,i))/de2a(j,i)*100.0_dp
                            WRITE(*,'(10F12.6)')de2a
                            WRITE(*,*)
                            WRITE(*,'(10F12.6)')de2a2
                            WRITE(*,*)invert,SUM(ABS(de2a2-de2a))
                            WRITE(*,'(10F12.6)')de2a2-de2a
                            STOP
                         END IF
                      END IF
                   END DO
                END DO
             END IF
             IF (l_denuc) THEN
                DO i = 1, 3
                   IF (ABS(denuc2(i))>delta*10) THEN
                      IF (ABS((denuc(i)-denuc2(i))/denuc(i))*100.0_dp>1.0_dp) THEN
                         WRITE(*,*)"check denuc",i,SUM(ABS(denuc-denuc2))
                         WRITE(*,*)invert,denuc-denuc2
                         WRITE(*,'(3F12.6)')denuc
                         WRITE(*,'(3F12.6)')denuc2
                         STOP
                      END IF
                   END IF
                END DO
             END IF
          END IF
       END IF
    END IF
  END SUBROUTINE rotnuc_ana

!!****f* semi_empirical_int_ana/invert_integral [1.0] *
!!
!!   NAME
!!     invert_integral
!!
!!   FUNCTION
!!     Exploits inversion symmetry to avoid divergence
!!
!!   NOTES
!!    
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino - Zurich University
!!
!!   MODIFICATION HISTORY
!!     04.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE invert_integral(array, np)
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: array
    INTEGER, INTENT(IN)                      :: np

    INTEGER                                  :: i, i1, i2
    INTEGER, DIMENSION(2, 42)                :: map
    REAL(KIND=dp)                            :: tmp

    map = RESHAPE(map_full,(/2,42/))
    DO i = 1, np
       i1= map(1,i)
       i2= map(2,i)

       tmp=array(i1)
       array(i1)=array(i2)
       array(i2)=tmp
    END DO

  END SUBROUTINE invert_integral

!!****f* semi_empirical_int_ana/invert_derivative [1.0] *
!!
!!   NAME
!!     invert_derivative
!!
!!   FUNCTION
!!     Exploits inversion symmetry to avoid divergence
!!
!!   NOTES
!!    
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino - Zurich University
!!
!!   MODIFICATION HISTORY
!!     04.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE invert_derivative(array, np)
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: array
    INTEGER, INTENT(IN)                      :: np

    INTEGER                                  :: i, i1, i2, j, m
    INTEGER, DIMENSION(2, 42)                :: map
    REAL(KIND=dp)                            :: tmp

    map = RESHAPE(map_full,(/2,42/))
    DO j=1,3
       DO i =1, np
          i1= map(1,i)
          i2= map(2,i)

          tmp=array(i1,j)
          array(i1,j)=array(i2,j)
          array(i2,j)=tmp
       END DO
    END DO

    DO m=1,SIZE(array,1)
       tmp=array(m,2)
       array(m,2)=array(m,3)
       array(m,3)=tmp
    ENDDO
  END SUBROUTINE invert_derivative

!!****f* semi_empirical_int_ana/check_dnucint_ana [1.0] *
!!
!!   NAME
!!     check_dnucint_ana
!!
!!   FUNCTION
!!     Check Numerical Vs Analytical
!!
!!   NOTES
!!     Debug routine
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino - Zurich University
!!
!!   MODIFICATION HISTORY
!!     04.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE check_dnucint_ana (sepi,sepj,r,dssss,dccore,itype,shortrange)

    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: r, dssss
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: dccore
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: shortrange

    INTEGER                                  :: i, j
    REAL(dp)                                 :: delta, nssss, od, rn, ssssm, &
                                                ssssp
    REAL(dp), DIMENSION(4, 2)                :: ccorem, ccorep, nccore

    delta = 1.0E-8_dp
    od = 0.5_dp/delta
    rn = r + delta
    CALL nucint(sepi,sepj,rn,ssssp,ccorep,itype,shortrange)
    rn = r - delta
    CALL nucint(sepi,sepj,rn,ssssm,ccorem,itype,shortrange)
    nccore = od * (ccorep - ccorem)
    nssss  = od * (ssssp - ssssm)
    ! check
    DO i = 1, 2
       DO j = 1, 4
          IF      ((ABS(nccore(j,i))>EPSILON(0.0_dp)).AND.(ABS(dccore(j,i))>EPSILON(0.0_dp))) THEN
             IF ((nccore(j,i)-dccore(j,i))/nccore(j,i)*100._dp>0.1_dp) THEN
                WRITE(*,*) j,i,ABS(nccore(j,i)-dccore(j,i))/nccore(j,i)*100._dp
                STOP
             END IF
          ELSEIF  ((ABS(nccore(j,i))<EPSILON(0.0_dp)).AND.(ABS(dccore(j,i))>EPSILON(0.0_dp))) THEN
             WRITE(*,*) j,i,ABS(nccore(j,i)-dccore(j,i))/nccore(j,i)*100._dp
             STOP
          ELSEIF  ((ABS(nccore(j,i))>EPSILON(0.0_dp)).AND.(ABS(dccore(j,i))<EPSILON(0.0_dp))) THEN
             WRITE(*,*) j,i,ABS(nccore(j,i)-dccore(j,i))/nccore(j,i)*100._dp
             STOP
          END IF
       END DO
    END DO
    IF (ABS(nssss)>EPSILON(0.0_dp)) THEN
       IF ((nssss-dssss)/nssss*100._dp>0.1_dp) THEN
          WRITE(*,*) nssss, dssss, (nssss-dssss)/nssss*100._dp
          STOP
       END IF
    END IF
  END SUBROUTINE check_dnucint_ana

!!****f* semi_empirical_int_ana/dnucint_ana [1.0] *
!!
!!   NAME
!!     dnucint_ana
!!
!!   FUNCTION
!!     Calculates the nuclear attraction integrals.
!!
!!   NOTES
!!     Analytical version - Analytical evaluation of gradients
!!     Teodoro Laino - Zurich University 04.2007
!!     routine adapted from mopac7 (repp)
!!     vector version written by Ernest R. Davidson, Indiana University
!!
!!   INPUTS
!!      on input rij     = interatomic distance
!!               sepi    = paramters of atom i
!!               sepj    = paramters of atom j
!!               itype   = type of semi_empirical model
!!                         extension to the original routine to compute qm/mm integrals
!!     on output core    = derivative of 4 X 2 array of electron-core attraction integrals
!!               ssss    = derivative of (ssss) integral
!!                         derivatives are intended w.r.t. rij
!!        The storage of the nuclear attraction integrals  core(kl/ij) iS
!!        (SS/)=1,   (SO/)=2,   (OO/)=3,   (PP/)=4
!!        where ij=1 if the orbitals centred on atom i,  =2 if on atom j.
!!   AUTHOR
!!     Teodoro Laino - Zurich University
!!
!!   MODIFICATION HISTORY
!!     04.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE dnucint_ana ( sepi, sepj, rij, ssss, core, dssss, dcore, itype, shortrange )
    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), INTENT(OUT)                    :: ssss
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: core
    REAL(dp), INTENT(OUT)                    :: dssss
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: dcore
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: shortrange

    LOGICAL                                  :: si, sj
    REAL(dp) :: ade, adi, adj, aed, aee, aeq, ami, amj, aqe, aqi, aqj, da, &
      db, ddi, ddj, ddze, dedz, dee, deqxx, deqzz, dft, dqxxe, dqzze, dze, &
      edz, ee, eqxx, eqzz, fac, ft, qa, qb, qqi, qqj, qxxe, qzze, r, rsq, &
      xxx, zi, zj
    REAL(dp), DIMENSION(14)                  :: arg, darg, dri, dsqr, ri, sqr

    ri    = 0._dp
    dri   = 0._dp
    core  = 0._dp
    dcore = 0._dp
    r=rij
    si = (sepi%natorb >= 3)
    sj = (sepj%natorb >= 3)
    zi = sepi%zeff
    zj = sepj%zeff
    ddi = sepi%dd
    ddj = sepj%dd
    qqi = sepi%qq
    qqj = sepj%qq

    IF ((.NOT.si) .AND. (.NOT.sj)) THEN
       !
       !  hydrogen - hydrogen  (SS/SS)
       !
       ami = sepi%am
       amj = sepj%am
       aee = 0.0_dp
       IF (itype/=0) THEN
          aee = 0.5_dp/ami + 0.5_dp/amj
          aee = aee * aee
       END IF
       ri(1) = 1._dp/SQRT(r*r+aee)
       core(1,1) = zj*ri(1)
       core(1,2) = zi*ri(1)       
       !
       fac   = -r/(r*r+aee)
       dri(1) = ri(1)*fac
       dcore(1,1) = zj*dri(1)
       dcore(1,2) = zi*dri(1)

    ELSE IF (si .AND. (.NOT.sj)) THEN
       !
       !  heavy atom - hydrogen
       !
       ami = sepi%am
       adi = sepi%ad
       aqi = sepi%aq
       amj = sepj%am
       da=ddi
       qa=qqi * 2._dp
       aee = 0.0_dp
       ade = 0.0_dp
       aqe = 0.0_dp
       IF (itype /=0) THEN
          aee = 0.5_dp/ami + 0.5_dp/amj
          ade = 0.5_dp/adi + 0.5_dp/amj
          aqe = 0.5_dp/aqi + 0.5_dp/amj
          aee = aee * aee
          ade = ade * ade
          aqe = aqe * aqe
       END IF
       rsq = r*r
       arg(1) = rsq + aee
       darg(1)= 2.0_dp *r 
       xxx = r+da
       arg(2) = xxx*xxx + ade
       darg(2)= 2.0_dp *xxx
       xxx = r-da
       arg(3) = xxx*xxx + ade
       darg(3)= 2.0_dp *xxx
       xxx = r+qa
       arg(4) = xxx*xxx + aqe
       darg(4)= 2.0_dp *xxx
       xxx = r-qa
       arg(5) = xxx*xxx + aqe
       darg(5)= 2.0_dp *xxx
       arg(6) = rsq + aqe
       darg(6)= 2.0_dp *r
       arg(7) = arg(6) + qa*qa
       darg(7)= darg(6)
       !
       sqr(1:7)  = SQRT(arg(1:7))
       dsqr(1:7) = 0.5_dp/SQRT(arg(1:7))*darg(1:7)
       !
       ee = 1._dp/sqr(1)
       ri(1) = ee
       ri(2) = 0.5_dp/sqr(2) - 0.5_dp/sqr(3)
       ri(3) = ee + 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
       ri(4) = ee + 0.5_dp/sqr(7) - 0.5_dp/sqr(6)
       core(1,1) = zj*ri(1)
       core(1,2) = zi*ri(1)
       core(2,1) = zj*ri(2)
       core(3,1) = zj*ri(3)
       core(4,1) = zj*ri(4)
       !
       dee = -1._dp/(sqr(1)**2)*dsqr(1)
       dri(1) = dee
       dri(2) = -0.5_dp/(sqr(2)**2)*dsqr(2) + 0.5_dp/(sqr(3)**2)*dsqr(3)
       dri(3) = dee - 0.25_dp/(sqr(4)**2)*dsqr(4) - 0.25_dp/(sqr(5)**2)*dsqr(5) + 0.5_dp/(sqr(6)**2)*dsqr(6)
       dri(4) = dee - 0.5_dp/(sqr(7)**2)*dsqr(7) + 0.5_dp/(sqr(6)**2)*dsqr(6)
       dcore(1,1) = zj*dri(1)
       dcore(1,2) = zi*dri(1)
       dcore(2,1) = zj*dri(2)
       dcore(3,1) = zj*dri(3)
       dcore(4,1) = zj*dri(4)

    ELSE IF ((.NOT.si).AND.sj) THEN
       !
       !  hydrogen - heavy atom
       !
       ami = sepi%am
       amj = sepj%am
       adj = sepj%ad
       aqj = sepj%aq
       db=ddj
       qb=qqj * 2._dp
       aee = 0.0_dp
       aed = 0.0_dp
       aeq = 0.0_dp
       IF (itype /=0) THEN
          aee = 0.5_dp/ami + 0.5_dp/amj
          aed = 0.5_dp/ami + 0.5_dp/adj
          aeq = 0.5_dp/ami + 0.5_dp/aqj
          aee = aee * aee
          aed = aed * aed
          aeq = aeq * aeq
       END IF
       rsq = r*r
       arg(1) = rsq + aee
       darg(1)= 2.0_dp *r 
       xxx = r-db
       arg(2) = xxx*xxx + aed
       darg(2)= 2.0_dp *xxx
       xxx = r+db
       arg(3) = xxx*xxx + aed
       darg(3)= 2.0_dp *xxx
       xxx = r-qb
       arg(4) = xxx*xxx + aeq
       darg(4)= 2.0_dp *xxx
       xxx = r+qb
       arg(5) = xxx*xxx + aeq
       darg(5)= 2.0_dp *xxx
       arg(6) = rsq + aeq
       darg(6)= 2.0_dp *r
       arg(7) = arg(6) + qb*qb
       darg(7)= darg(6)
       sqr(1:7) = SQRT(arg(1:7))
       dsqr(1:7) = 0.5_dp/SQRT(arg(1:7))*darg(1:7)
       !
       ee = 1._dp/sqr(1)
       ri(1) = ee
       ri(5) = 0.5_dp/sqr(2)  - 0.5_dp/sqr(3)
       ri(11) = ee + 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
       ri(12) = ee + 0.5_dp/sqr(7) - 0.5_dp/sqr(6)
       core(1,1) = zj*ri(1)
       core(1,2) = zi*ri(1)
       core(2,2) = zi*ri(5)
       core(3,2) = zi*ri(11)
       core(4,2) = zi*ri(12)       
       !
       dee = -1._dp/(sqr(1)**2)*dsqr(1)
       dri(1) = dee
       dri(5) = -0.5_dp/(sqr(2)**2)*dsqr(2) + 0.5_dp/(sqr(3)**2)*dsqr(3)
       dri(11) = dee - 0.25_dp/(sqr(4)**2)*dsqr(4) - 0.25_dp/(sqr(5)**2)*dsqr(5) + 0.5_dp/(sqr(6)**2)*dsqr(6)
       dri(12) = dee - 0.5_dp/(sqr(7)**2)*dsqr(7) + 0.5_dp/(sqr(6)**2)*dsqr(6)
       dcore(1,1) = zj*dri(1)
       dcore(1,2) = zi*dri(1)
       dcore(2,2) = zi*dri(5)
       dcore(3,2) = zi*dri(11)
       dcore(4,2) = zi*dri(12)

    ELSE
       !
       !  heavy atom - heavy atom
       !
       !  define charge separations.
       da=ddi
       db=ddj
       qa=qqi * 2._dp
       qb=qqj * 2._dp

       ami = sepi%am
       amj = sepj%am
       adi = sepi%ad
       adj = sepj%ad
       aqi = sepi%aq
       aqj = sepj%aq
       ade = 0.0_dp
       aee = 0.0_dp
       aqe = 0.0_dp
       aed = 0.0_dp
       aeq = 0.0_dp
       IF (itype /=0) THEN
          aee = 0.5_dp/ami + 0.5_dp/amj
          ade = 0.5_dp/adi + 0.5_dp/amj
          aqe = 0.5_dp/aqi + 0.5_dp/amj
          aed = 0.5_dp/ami + 0.5_dp/adj
          aeq = 0.5_dp/ami + 0.5_dp/aqj
          ade = ade * ade
          aee = aee * aee
          aqe = aqe * aqe
          aed = aed * aed
          aeq = aeq * aeq
       END IF
       rsq = r * r
       arg(1) = rsq + aee
       darg(1)= 2.0_dp *r 
       xxx = r + da
       arg(2) = xxx * xxx + ade
       darg(2)= 2.0_dp *xxx
       xxx = r - da
       arg(3) = xxx*xxx + ade
       darg(3)= 2.0_dp *xxx
       xxx = r - qa
       arg(4) = xxx*xxx + aqe
       darg(4)= 2.0_dp *xxx
       xxx = r + qa
       arg(5) = xxx*xxx + aqe
       darg(5)= 2.0_dp *xxx
       arg(6) = rsq + aqe
       darg(6)= 2.0_dp *r
       arg(7) = arg(6) + qa*qa
       darg(7)= darg(6)
       xxx = r-db
       arg(8) = xxx*xxx + aed
       darg(8)= 2.0_dp *xxx
       xxx = r+db
       arg(9) = xxx*xxx + aed
       darg(9)= 2.0_dp *xxx
       xxx = r - qb
       arg(10) = xxx*xxx + aeq
       darg(10)= 2.0_dp *xxx
       xxx = r + qb
       arg(11) = xxx*xxx + aeq
       darg(11)= 2.0_dp *xxx
       arg(12) = rsq + aeq
       darg(12)= 2.0_dp *r
       arg(13) = arg(12) + qb*qb
       darg(13)= darg(12)
       sqr(1:13) = SQRT(arg(1:13))
       dsqr(1:13) = 0.5_dp/SQRT(arg(1:13))*darg(1:13)
       !
       ee = 1._dp/sqr(1)
       dze = -0.5_dp/sqr(2) + 0.5_dp/sqr(3)
       qzze = 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
       qxxe = 0.5_dp/sqr(7) - 0.5_dp/sqr(6)
       edz = - 0.5_dp/sqr(8) + 0.5_dp/sqr(9)
       eqzz  = 0.25_dp/sqr(10) + 0.25_dp/sqr(11) - 0.5_dp/sqr(12)
       eqxx  = 0.5_dp/sqr(13) - 0.5_dp/sqr(12)
       ri(1) = ee
       ri(2) = -dze
       ri(3) = ee + qzze
       ri(4) = ee + qxxe
       ri(5) = -edz
       ri(11) =  ee + eqzz
       ri(12) =  ee + eqxx
       core(1,1) = zj*ri(1)
       core(2,1) = zj*ri(2)
       core(3,1) = zj*ri(3)
       core(4,1) = zj*ri(4)
       core(1,2) = zi*ri(1)
       core(2,2) = zi*ri(5)
       core(3,2) = zi*ri(11)
       core(4,2) = zi*ri(12)
       !
       dee = -1._dp/(sqr(1)**2)*dsqr(1)
       ddze =  0.5_dp/(sqr(2)**2)*dsqr(2) - 0.5_dp/(sqr(3)**2)*dsqr(3)
       dqzze = -0.25_dp/(sqr(4)**2)*dsqr(4) - 0.25_dp/(sqr(5)**2)*dsqr(5) + 0.5_dp/(sqr(6)**2)*dsqr(6)
       dqxxe = -0.5_dp/(sqr(7)**2)*dsqr(7) + 0.5_dp/(sqr(6)**2)*dsqr(6)
       dedz =  0.5_dp/(sqr(8)**2)*dsqr(8) - 0.5_dp/(sqr(9)**2)*dsqr(9)
       deqzz  = -0.25_dp/(sqr(10)**2)*dsqr(10) - 0.25_dp/(sqr(11)**2)*dsqr(11) + 0.5_dp/(sqr(12)**2)*dsqr(12)
       deqxx  = -0.5_dp/(sqr(13)**2)*dsqr(13) + 0.5_dp/(sqr(12)**2)*dsqr(12)
       dri(1) = dee
       dri(2) = -ddze
       dri(3) = dee + dqzze
       dri(4) = dee + dqxxe
       dri(5) = -dedz
       dri(11) =  dee + deqzz
       dri(12) =  dee + deqxx
       dcore(1,1) = zj*dri(1)
       dcore(2,1) = zj*dri(2)
       dcore(3,1) = zj*dri(3)
       dcore(4,1) = zj*dri(4)
       dcore(1,2) = zi*dri(1)
       dcore(2,2) = zi*dri(5)
       dcore(3,2) = zi*dri(11)
       dcore(4,2) = zi*dri(12)
    END IF
    IF ( shortrange ) THEN
      IF ((.NOT.si) .AND. (.NOT.sj)) THEN
         !
         !  hydrogen - hydrogen  (SS/SS)
         !
         ri(1) = 1._dp/SQRT(r*r+aee)
         core(1,1) = core(1,1) - zj*ri(1)
         core(1,2) = core(1,2) - zi*ri(1)       
         !
         fac   = -r/(r*r+aee)
         dri(1) = ri(1)*fac
         dcore(1,1) = dcore(1,1) - zj*dri(1)
         dcore(1,2) = dcore(1,2) - zi*dri(1)
  
      ELSE IF (si .AND. (.NOT.sj)) THEN
         !
         !  heavy atom - hydrogen
         !
         da=ddi
         qa=qqi * 2._dp
         rsq = r*r
         arg(1) = rsq
         darg(1)= 2.0_dp *r 
         xxx = r+da
         arg(2) = xxx*xxx
         darg(2)= 2.0_dp *xxx
         xxx = r-da
         arg(3) = xxx*xxx
         darg(3)= 2.0_dp *xxx
         xxx = r+qa
         arg(4) = xxx*xxx
         darg(4)= 2.0_dp *xxx
         xxx = r-qa
         arg(5) = xxx*xxx
         darg(5)= 2.0_dp *xxx
         arg(6) = rsq
         darg(6)= 2.0_dp *r
         arg(7) = arg(6) + qa*qa
         darg(7)= darg(6)
         !
         sqr(1:7)  = SQRT(arg(1:7))
         dsqr(1:7) = 0.5_dp/SQRT(arg(1:7))*darg(1:7)
         !
         ee = 1._dp/sqr(1)
         ri(1) = ee
         ri(2) = 0.5_dp/sqr(2) - 0.5_dp/sqr(3)
         ri(3) = ee + 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
         ri(4) = ee + 0.5_dp/sqr(7) - 0.5_dp/sqr(6)
         core(1,1) = core(1,1) - zj*ri(1)
         core(1,2) = core(1,2) - zi*ri(1)
         core(2,1) = core(2,1) - zj*ri(2)
         core(3,1) = core(3,1) - zj*ri(3)
         core(4,1) = core(4,1) - zj*ri(4)
         !
         dee = -1._dp/(sqr(1)**2)*dsqr(1)
         dri(1) = dee
         dri(2) = -0.5_dp/(sqr(2)**2)*dsqr(2) + 0.5_dp/(sqr(3)**2)*dsqr(3)
         dri(3) = dee - 0.25_dp/(sqr(4)**2)*dsqr(4) - 0.25_dp/(sqr(5)**2)*dsqr(5) + 0.5_dp/(sqr(6)**2)*dsqr(6)
         dri(4) = dee - 0.5_dp/(sqr(7)**2)*dsqr(7) + 0.5_dp/(sqr(6)**2)*dsqr(6)
         dcore(1,1) = dcore(1,1) - zj*dri(1)
         dcore(1,2) = dcore(1,2) - zi*dri(1)
         dcore(2,1) = dcore(2,1) - zj*dri(2)
         dcore(3,1) = dcore(3,1) - zj*dri(3)
         dcore(4,1) = dcore(4,1) - zj*dri(4)
  
      ELSE IF ((.NOT.si).AND.sj) THEN
         !
         !  hydrogen - heavy atom
         !
         db=ddj
         qb=qqj * 2._dp
         rsq = r*r
         arg(1) = rsq
         darg(1)= 2.0_dp *r 
         xxx = r-db
         arg(2) = xxx*xxx
         darg(2)= 2.0_dp *xxx
         xxx = r+db
         arg(3) = xxx*xxx
         darg(3)= 2.0_dp *xxx
         xxx = r-qb
         arg(4) = xxx*xxx
         darg(4)= 2.0_dp *xxx
         xxx = r+qb
         arg(5) = xxx*xxx
         darg(5)= 2.0_dp *xxx
         arg(6) = rsq
         darg(6)= 2.0_dp *r
         arg(7) = arg(6) + qb*qb
         darg(7)= darg(6)
         sqr(1:7) = SQRT(arg(1:7))
         dsqr(1:7) = 0.5_dp/SQRT(arg(1:7))*darg(1:7)
         !
         ee = 1._dp/sqr(1)
         ri(1) = ee
         ri(5) = 0.5_dp/sqr(2)  - 0.5_dp/sqr(3)
         ri(11) = ee + 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
         ri(12) = ee + 0.5_dp/sqr(7) - 0.5_dp/sqr(6)
         core(1,1) = core(1,1) - zj*ri(1)
         core(1,2) = core(1,2) - zi*ri(1)
         core(2,2) = core(2,2) - zi*ri(5)
         core(3,2) = core(3,2) - zi*ri(11)
         core(4,2) = core(4,2) - zi*ri(12)       
         !
         dee = -1._dp/(sqr(1)**2)*dsqr(1)
         dri(1) = dee
         dri(5) = -0.5_dp/(sqr(2)**2)*dsqr(2) + 0.5_dp/(sqr(3)**2)*dsqr(3)
         dri(11) = dee - 0.25_dp/(sqr(4)**2)*dsqr(4) - 0.25_dp/(sqr(5)**2)*dsqr(5) + 0.5_dp/(sqr(6)**2)*dsqr(6)
         dri(12) = dee - 0.5_dp/(sqr(7)**2)*dsqr(7) + 0.5_dp/(sqr(6)**2)*dsqr(6)
         dcore(1,1) = dcore(1,1) - zj*dri(1)
         dcore(1,2) = dcore(1,2) - zi*dri(1)
         dcore(2,2) = dcore(2,2) - zi*dri(5)
         dcore(3,2) = dcore(3,2) - zi*dri(11)
         dcore(4,2) = dcore(4,2) - zi*dri(12)
  
      ELSE
         !
         !  heavy atom - heavy atom
         !
         !  define charge separations.
         da=ddi
         db=ddj
         qa=qqi * 2._dp
         qb=qqj * 2._dp
  
         rsq = r * r
         arg(1) = rsq
         darg(1)= 2.0_dp *r 
         xxx = r + da
         arg(2) = xxx * xxx
         darg(2)= 2.0_dp *xxx
         xxx = r - da
         arg(3) = xxx*xxx
         darg(3)= 2.0_dp *xxx
         xxx = r - qa
         arg(4) = xxx*xxx
         darg(4)= 2.0_dp *xxx
         xxx = r + qa
         arg(5) = xxx*xxx
         darg(5)= 2.0_dp *xxx
         arg(6) = rsq
         darg(6)= 2.0_dp *r
         arg(7) = arg(6) + qa*qa
         darg(7)= darg(6)
         xxx = r-db
         arg(8) = xxx*xxx
         darg(8)= 2.0_dp *xxx
         xxx = r+db
         arg(9) = xxx*xxx
         darg(9)= 2.0_dp *xxx
         xxx = r - qb
         arg(10) = xxx*xxx
         darg(10)= 2.0_dp *xxx
         xxx = r + qb
         arg(11) = xxx*xxx
         darg(11)= 2.0_dp *xxx
         arg(12) = rsq
         darg(12)= 2.0_dp *r
         arg(13) = arg(12) + qb*qb
         darg(13)= darg(12)
         sqr(1:13) = SQRT(arg(1:13))
         dsqr(1:13) = 0.5_dp/SQRT(arg(1:13))*darg(1:13)
         !
         ee = 1._dp/sqr(1)
         dze = -0.5_dp/sqr(2) + 0.5_dp/sqr(3)
         qzze = 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
         qxxe = 0.5_dp/sqr(7) - 0.5_dp/sqr(6)
         edz = - 0.5_dp/sqr(8) + 0.5_dp/sqr(9)
         eqzz  = 0.25_dp/sqr(10) + 0.25_dp/sqr(11) - 0.5_dp/sqr(12)
         eqxx  = 0.5_dp/sqr(13) - 0.5_dp/sqr(12)
         ri(1) = ee
         ri(2) = -dze
         ri(3) = ee + qzze
         ri(4) = ee + qxxe
         ri(5) = -edz
         ri(11) =  ee + eqzz
         ri(12) =  ee + eqxx
         core(1,1) = core(1,1) - zj*ri(1)
         core(2,1) = core(2,1) - zj*ri(2)
         core(3,1) = core(3,1) - zj*ri(3)
         core(4,1) = core(4,1) - zj*ri(4)
         core(1,2) = core(1,2) - zi*ri(1)
         core(2,2) = core(2,2) - zi*ri(5)
         core(3,2) = core(3,2) - zi*ri(11)
         core(4,2) = core(4,2) - zi*ri(12)
         !
         dee = -1._dp/(sqr(1)**2)*dsqr(1)
         ddze =  0.5_dp/(sqr(2)**2)*dsqr(2) - 0.5_dp/(sqr(3)**2)*dsqr(3)
         dqzze = -0.25_dp/(sqr(4)**2)*dsqr(4) - 0.25_dp/(sqr(5)**2)*dsqr(5) + 0.5_dp/(sqr(6)**2)*dsqr(6)
         dqxxe = -0.5_dp/(sqr(7)**2)*dsqr(7) + 0.5_dp/(sqr(6)**2)*dsqr(6)
         dedz =  0.5_dp/(sqr(8)**2)*dsqr(8) - 0.5_dp/(sqr(9)**2)*dsqr(9)
         deqzz  = -0.25_dp/(sqr(10)**2)*dsqr(10) - 0.25_dp/(sqr(11)**2)*dsqr(11) + 0.5_dp/(sqr(12)**2)*dsqr(12)
         deqxx  = -0.5_dp/(sqr(13)**2)*dsqr(13) + 0.5_dp/(sqr(12)**2)*dsqr(12)
         dri(1) = dee
         dri(2) = -ddze
         dri(3) = dee + dqzze
         dri(4) = dee + dqxxe
         dri(5) = -dedz
         dri(11) =  dee + deqzz
         dri(12) =  dee + deqxx
         dcore(1,1) = dcore(1,1) - zj*dri(1)
         dcore(2,1) = dcore(2,1) - zj*dri(2)
         dcore(3,1) = dcore(3,1) - zj*dri(3)
         dcore(4,1) = dcore(4,1) - zj*dri(4)
         dcore(1,2) = dcore(1,2) - zi*dri(1)
         dcore(2,2) = dcore(2,2) - zi*dri(5)
         dcore(3,2) = dcore(3,2) - zi*dri(11)
         dcore(4,2) = dcore(4,2) - zi*dri(12)
      END IF
    END IF
    ssss = ri(1)
    ! Tapering function
    ft = 1.0_dp
    dft= 0.0_dp
    IF (itype/=0) THEN
       ft = taper ( rij )
       dft= dtaper_ana ( rij )
    END IF
    ! If no orbitals are present on the atom let's zero the integral contribution
    IF (sepi%natorb == -1) core(:,1) = 0.0_dp
    IF (sepj%natorb == -1) core(:,2) = 0.0_dp
    ! First derivatives..
    dssss      = dri(1)
    dcore(:,:) = ft*dcore(:,:)+dft*core(:,:)
    dssss      = ft*dssss+dft*ssss
    ! Then values of functions
    core(:,:)  = ft*core(:,:)
    ssss       = ft*ssss
    ! If no orbitals are present on the atom let's zero the integral contribution
    IF (sepi%natorb == -1) dcore(:,1) = 0.0_dp
    IF (sepj%natorb == -1) dcore(:,2) = 0.0_dp

  END SUBROUTINE dnucint_ana

  FUNCTION dtaper_ana ( rij ) RESULT(dft)
    REAL(KIND=dp), INTENT(IN)                :: rij
    REAL(KIND=dp)                            :: dft

    REAL(KIND=dp)                            :: dr, ep

    dft = 0.0_dp
    IF ( taper_fn_init ) THEN
       dr = -al*(rij-r0)
       IF ( dr < 20._dp ) THEN
          ep = EXP(dr)
          dft = 1._dp/(1._dp+ep)**2 * ep * (-al)
       END IF
    END IF

  END FUNCTION dtaper_ana

!!****f* semi_empirical_int_ana/rotint_ana [1.0] *
!!
!!   NAME
!!     rotint_ana
!!
!!   FUNCTION
!!     calculates the derivative of the two-particle interactions
!!
!!   NOTES
!!     Analytical version - Analytical evaluation of gradients
!!     Teodoro Laino - Zurich University 04.2007
!!     routine adapted from mopac7 (repp)
!!     vector version written by Ernest R. Davidson, Indiana University
!!
!!   INPUTS
!!    on input  sepi   = Atomic parameters of first atom
!!              sepj   = Atomic parameters of second atom
!!              rijv   = Coordinate vector i -> j 
!!    
!!    on output w      = Array of two-electron repulsion integrals.
!!
!!   AUTHOR
!!     Teodoro Laino - Zurich University
!!
!!   MODIFICATION HISTORY
!!     04.2007 created [tlaino]
!!
!!*** **********************************************************************
  RECURSIVE SUBROUTINE rotint_ana (sepi,sepj,rijv,w,dw,shortrange)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rijv
    REAL(dp), DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: w
    REAL(dp), DIMENSION(:, :), INTENT(OUT), &
      OPTIONAL                               :: dw
    LOGICAL, INTENT(IN)                      :: shortrange

    INTEGER                                  :: i, j
    LOGICAL                                  :: invert, l_w, lgrad, si, sj
    REAL(dp) :: a, delta, rij, xtmp, xx11, xx21, xx22, xx31, xx32, xx33, &
      xy11, xy21, xy22, xy31, xy32, xz11, xz21, xz22, xz31, xz32, xz33, yy11, &
      yy21, yy22, yyzz11, yyzz21, yyzz22, yz11, yz21, yz22, yz31, yz32, zz11, &
      zz21, zz22, zz31, zz32, zz33
    REAL(dp), DIMENSION(100)                 :: w2
    REAL(dp), DIMENSION(100, 3)              :: dw2
    REAL(dp), DIMENSION(22)                  :: dri1, ri
    REAL(dp), DIMENSION(22, 3)               :: dri
    REAL(dp), DIMENSION(3) :: da, drij, dx1, dx2, dx3, dxx11, dxx21, dxx22, &
      dxx31, dxx32, dxx33, dxy11, dxy21, dxy22, dxy31, dxy32, dxz11, dxz21, &
      dxz22, dxz31, dxz32, dxz33, dyy11, dyy21, dyy22, dyyzz11, dyyzz21, &
      dyyzz22, dyz11, dyz21, dyz22, dyz31, dyz32, dzz11, dzz21, dzz22, dzz31, &
      dzz32, dzz33, x, y, z
    REAL(dp), DIMENSION(3, 3)                :: dx, dy, dz

    l_w  = PRESENT(w)
    lgrad = PRESENT(dw)
    IF (l_w)    w = 0.0_dp
    IF (lgrad) dw = 0.0_dp
    x(:)=-rijv(:)
    rij=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
    IF (rij < 0.00002_dp) THEN
       ! SMALL RIJ CASE
       ! w is zero 
    ELSEIF (l_w.OR.lgrad) THEN
       ! The repulsion integrals over molecular frame (w) are stored in the
       ! order in which they will later be used.  ie.  (i,j/k,l) where
       ! j.le.i  and  l.le.k     and l varies most rapidly and i least
       ! rapidly.  (anti-normal computer storage)
       rij = SQRT(rij)
       ! Compute integrals in diatomic frame as well their derivatives (if requested)
       IF (lgrad) THEN
          CALL dterep_ana(sepi,sepj,rij,ri,dri1,shortrange)
          IF (debug_this_module) THEN
             CALL check_dterep_ana(sepi, sepj, rij, ri, dri1, shortrange)
          END IF
       ELSE
          CALL terep(sepi,sepj,rij,ri,shortrange)
       END IF
       a=1._dp/rij
       x(1) = x(1)*a
       x(2) = x(2)*a
       x(3) = x(3)*a
       ! Possibly Invert Frame
       IF (ABS(x(3)) > 0.99999999_dp) THEN
          ! In order to avoid divergence just change Z axes into the Y axes
          ! all quantities are rotational/invertion invariant..
          invert = .TRUE.
          xtmp=x(3)
          x(3)=x(2)
          x(2)=xtmp
       ELSE
          invert = .FALSE.
       END IF
       IF (lgrad) THEN
          drij(:) = -x(:)
          DO i = 1, 22
             dri(i,:) = dri1(i)*drij
          END DO
          da  = -a**2*drij
          dx1 = -(/1.0_dp,0.0_dp,0.0_dp/)
          dx2 = -(/0.0_dp,1.0_dp,0.0_dp/)
          dx3 = -(/0.0_dp,0.0_dp,1.0_dp/)
          dx(1,:) = dx1*a+(x(1)/a)*da
          dx(2,:) = dx2*a+(x(2)/a)*da
          dx(3,:) = dx3*a+(x(3)/a)*da
       END IF
       z(3)=SQRT(1._dp-x(3)*x(3))
       a=1._dp/z(3)
       y(1)=-a*x(2)*SIGN(1._dp,x(1))
       y(2)=ABS(a*x(1))
       y(3)=0._dp
       z(1)=-a*x(1)*x(3)
       z(2)=-a*x(2)*x(3)
       ! Analytical Gradients
       IF (lgrad) THEN
          dz(3,:) = -a*x(3)*dx(3,:)
          da      = -a**2*dz(3,:)
          dy(1,:) = -da*x(2)*SIGN(1._dp,x(1))-a*dx(2,:)*SIGN(1._dp,x(1))
          dy(2,:) = SIGN(1._dp,a*x(1))*(da*x(1)+a*dx(1,:))
          dy(3,:) = 0.0_dp
          dz(1,:) =-da*x(1)*x(3)-a*dx(1,:)*x(3)-a*x(1)*dx(3,:)
          dz(2,:) =-da*x(2)*x(3)-a*dx(2,:)*x(3)-a*x(2)*dx(3,:)
       END IF
       si = (sepi%natorb > 1)
       sJ = (sepj%natorb > 1)
       IF ( si .OR. sj ) THEN
          xx11 = x(1)*x(1)
          xx21 = x(2)*x(1)
          xx22 = x(2)*x(2)
          xx31 = x(3)*x(1)
          xx32 = x(3)*x(2)
          xx33 = x(3)*x(3)
          yy11 = y(1)*y(1)
          YY21 = Y(2)*Y(1)
          yy22 = y(2)*y(2)
          zz11 = z(1)*z(1)
          zz21 = z(2)*z(1)
          zz22 = z(2)*z(2)
          zz31 = z(3)*z(1)
          zz32 = z(3)*z(2)
          zz33 = z(3)*z(3)
          yyzz11 = yy11+zz11
          yyzz21 = yy21+zz21
          yyzz22 = yy22+zz22
          xy11 = 2._dp*x(1)*y(1)
          xy21 =      x(1)*y(2)+x(2)*y(1)
          xy22 = 2._dp*x(2)*y(2)
          xy31 =      x(3)*y(1)
          xy32 =      x(3)*y(2)
          xz11 = 2._dp*x(1)*z(1)
          xz21 =      x(1)*z(2)+x(2)*z(1)
          xz22 = 2._dp*x(2)*z(2)
          xz31 =      x(1)*z(3)+x(3)*z(1)
          xz32 =      x(2)*z(3)+x(3)*z(2)
          xz33 = 2._dp*x(3)*z(3)
          yz11 = 2._dp*y(1)*z(1)
          yz21 =      y(1)*z(2)+y(2)*z(1)
          yz22 = 2._dp*y(2)*z(2)
          yz31 =      y(1)*z(3)
          yz32 =      y(2)*z(3)
          ! Analytical Gradients
          IF (lgrad) THEN
             dxx11 = dx(1,:)*x(1)+x(1)*dx(1,:)
             dxx21 = dx(2,:)*x(1)+x(2)*dx(1,:)
             dxx22 = dx(2,:)*x(2)+x(2)*dx(2,:)
             dxx31 = dx(3,:)*x(1)+x(3)*dx(1,:)
             dxx32 = dx(3,:)*x(2)+x(3)*dx(2,:)
             dxx33 = dx(3,:)*x(3)+x(3)*dx(3,:)
             dyy11 = dy(1,:)*y(1)+y(1)*dy(1,:)
             dyy21 = dy(2,:)*y(1)+y(2)*dy(1,:)
             dyy22 = dy(2,:)*y(2)+y(2)*dy(2,:)
             dzz11 = dz(1,:)*z(1)+z(1)*dz(1,:)
             dzz21 = dz(2,:)*z(1)+z(2)*dz(1,:)
             dzz22 = dz(2,:)*z(2)+z(2)*dz(2,:)
             dzz31 = dz(3,:)*z(1)+z(3)*dz(1,:)
             dzz32 = dz(3,:)*z(2)+z(3)*dz(2,:)
             dzz33 = dz(3,:)*z(3)+z(3)*dz(3,:)
             dyyzz11 = dyy11+dzz11
             dyyzz21 = dyy21+dzz21
             dyyzz22 = dyy22+dzz22
             dxy11 = 2._dp*dx(1,:)*y(1)+2._dp*x(1)*dy(1,:)
             dxy21 =      dx(1,:)*y(2)+x(1)*dy(2,:)+dx(2,:)*y(1)+x(2)*dy(1,:)
             dxy22 = 2._dp*dx(2,:)*y(2)+2._dp*x(2)*dy(2,:)
             dxy31 =      dx(3,:)*y(1)+x(3)*dy(1,:)
             dxy32 =      dx(3,:)*y(2)+x(3)*dy(2,:)
             dxz11 = 2._dp*dx(1,:)*z(1)+2._dp*x(1)*dz(1,:)
             dxz21 =      dx(1,:)*z(2)+x(1)*dz(2,:)+dx(2,:)*z(1)+x(2)*dz(1,:)
             dxz22 = 2._dp*dx(2,:)*z(2)+2._dp*x(2)*dz(2,:)
             dxz31 =      dx(1,:)*z(3)+x(1)*dz(3,:)+dx(3,:)*z(1)+x(3)*dz(1,:)
             dxz32 =      dx(2,:)*z(3)+x(2)*dz(3,:)+dx(3,:)*z(2)+x(3)*dz(2,:)
             dxz33 = 2._dp*dx(3,:)*z(3)+2._dp*x(3)*dz(3,:)
             dyz11 = 2._dp*dy(1,:)*z(1)+2._dp*y(1)*dz(1,:)
             dyz21 =      dy(1,:)*z(2)+y(1)*dz(2,:)+dy(2,:)*z(1)+y(2)*dz(1,:)
             dyz22 = 2._dp*dy(2,:)*z(2)+2._dp*y(2)*dz(2,:)
             dyz31 =      dy(1,:)*z(3)+y(1)*dz(3,:)
             dyz32 =      dy(2,:)*z(3)+y(2)*dz(3,:)
          END IF
       ENDIF
       IF (l_w) THEN
          w(:)=0._dp
          !(s s/s s)
          w(1)=ri(1)
          IF (sj) THEN
             !(s s/px s)
             w(2)=ri(5)*x(1)
             !(s s/px px)
             w(3)=ri(11)*xx11+ri(12)*yyzz11
             !(s s/py s)
             w(4)=ri(5)*x(2)
             !(s s/py px)
             w(5)=ri(11)*xx21+ri(12)*yyzz21
             !(s s/py py)
             w(6)=ri(11)*xx22+ri(12)*yyzz22
             !(s s/pz s)
             w(7)=ri(5)*x(3)
             !(s s/pz px)
             w(8)=ri(11)*xx31+ri(12)*zz31
             !(s s/pz py)
             w(9)=ri(11)*xx32+ri(12)*zz32
             !(s s/pz pz)
             w(10)=ri(11)*xx33+ri(12)*zz33
          END IF

          IF (si) THEN
             !(px s/s s)
             w(11)=ri(2)*x(1)
             IF (sj) THEN
                !(px s/px s)
                w(12)=ri(6)*xx11+ri(7)*yyzz11
                !(px s/px px)
                w(13)=x(1)*(ri(13)*xx11+ri(14)*yyzz11) &
                     +ri(15)*(y(1)*xy11+z(1)*xz11)
                !(px s/py s)
                w(14)=ri(6)*xx21+ri(7)*yyzz21
                !(px s/py px)
                w(15)=x(1)*(ri(13)*xx21+ri(14)*yyzz21) &
                     +ri(15)*(y(1)*xy21+z(1)*xz21)
                !(px s/py py)
                w(16)=x(1)*(ri(13)*xx22+ri(14)*yyzz22) &
                     +ri(15)*(y(1)*xy22+z(1)*xz22)
                !(px s/pz s)
                w(17)=ri(6)*xx31+ri(7)*zz31
                !(px s/pz px)
                w(18)=x(1)*(ri(13)*xx31+ri(14)*zz31) &
                     +ri(15)*(y(1)*xy31+z(1)*xz31)
                !(px s/pz py)
                w(19)=x(1)*(ri(13)*xx32+ri(14)*zz32) &
                     +ri(15)*(y(1)*xy32+z(1)*xz32)
                !(px s/pz pz)
                w(20)=x(1)*(ri(13)*xx33+ri(14)*zz33) &
                     +ri(15)*(          z(1)*xz33)
                !(px px/s s)
                w(21)=ri(3)*xx11+ri(4)*yyzz11
                !(px px/px s)
                w(22)=x(1)*(ri(8)*xx11+ri(9)*yyzz11) &
                     +ri(10)*(y(1)*xy11+z(1)*xz11)
                !(px px/px px)
                w(23) = &
                     (ri(16)*xx11+ri(17)*yyzz11)*xx11+ri(18)*xx11*yyzz11 &
                     +ri(19)*(yy11*yy11+zz11*zz11) &
                     +ri(20)*(xy11*xy11+xz11*xz11) &
                     +ri(21)*(yy11*zz11+zz11*yy11) &
                     +ri(22)*yz11*yz11
                !(px px/py s)
                w(24)=x(2)*(ri(8)*xx11+ri(9)*yyzz11) &
                     +ri(10)*(y(2)*xy11+z(2)*xz11)
                !(px px/py px)
                w(25) = &
                     (ri(16)*xx11+ri(17)*yyzz11)*xx21+ri(18)*xx11*yyzz21 &
                     +ri(19)*(yy11*yy21+zz11*zz21) &
                     +ri(20)*(xy11*xy21+xz11*xz21) &
                     +ri(21)*(yy11*zz21+zz11*yy21) &
                     +ri(22)*yz11*yz21
                !(px px/py py)
                w(26) = &
                     (ri(16)*xx11+ri(17)*yyzz11)*xx22+ri(18)*xx11*yyzz22 &
                     +ri(19)*(yy11*yy22+zz11*zz22) &
                     +ri(20)*(xy11*xy22+xz11*xz22) &
                     +ri(21)*(yy11*zz22+zz11*yy22) &
                     +ri(22)*yz11*yz22
                !(px px/pz s)
                w(27)=x(3)*(ri(8)*xx11+ri(9)*yyzz11) &
                     +ri(10)*(         +z(3)*xz11)
                !(px px/pz px)
                w(28) = &
                     (ri(16)*xx11+ri(17)*yyzz11)*xx31 &
                     +(ri(18)*xx11+ri(19)*zz11+ri(21)*yy11)*zz31 &
                     +ri(20)*(xy11*xy31+xz11*xz31) &
                     +ri(22)*yz11*yz31
                !(px px/pz py)
                w(29) = &
                     (ri(16)*xx11+ri(17)*yyzz11)*xx32 &
                     +(ri(18)*xx11+ri(19)*zz11+ri(21)*yy11)*zz32 &
                     +ri(20)*(xy11*xy32+xz11*xz32) &
                     +ri(22)*yz11*yz32
                !(px px/pz pz)
                w(30) = &
                     (ri(16)*xx11+ri(17)*yyzz11)*xx33 &
                     +(ri(18)*xx11+ri(19)*zz11+ri(21)*yy11)*zz33 &
                     +ri(20)*xz11*xz33
                !(py s/s s)
                w(31)=ri(2)*x(2)
                !(py s/px s)
                w(32)=ri(6)*xx21+ri(7)*yyzz21
                !(py s/px px)
                w(33)=x(2)*(ri(13)*xx11+ri(14)*yyzz11) &
                     +ri(15)*(y(2)*xy11+z(2)*xz11)
                !(py s/py s)
                w(34)=ri(6)*xx22+ri(7)*yyzz22
                !(py s/py px)
                w(35)=x(2)*(ri(13)*xx21+ri(14)*yyzz21) &
                     +ri(15)*(y(2)*xy21+z(2)*xz21)
                !(py s/py py)
                w(36)=x(2)*(ri(13)*xx22+ri(14)*yyzz22) &
                     +ri(15)*(y(2)*xy22+z(2)*xz22)
                !(py s/pz s)
                w(37)=ri(6)*xx32+ri(7)*zz32
                !(py s/pz px)
                w(38)=x(2)*(ri(13)*xx31+ri(14)*zz31) &
                     +ri(15)*(y(2)*xy31+z(2)*xz31)
                !(py s/pz py)
                w(39)=x(2)*(ri(13)*xx32+ri(14)*zz32) &
                     +ri(15)*(y(2)*xy32+z(2)*xz32)
                !(py s/pz pz)
                w(40)=x(2)*(ri(13)*xx33+ri(14)*zz33) &
                     +ri(15)*(         +z(2)*xz33)
                !(py px/s s)
                w(41)=ri(3)*xx21+ri(4)*yyzz21
                !(py px/px s)
                w(42)=x(1)*(ri(8)*xx21+ri(9)*yyzz21) &
                     +ri(10)*(y(1)*xy21+z(1)*xz21)
                !(py px/px px)
                w(43) = &
                     (ri(16)*xx21+ri(17)*yyzz21)*xx11+ri(18)*xx21*yyzz11 &
                     +ri(19)*(yy21*yy11+zz21*zz11) &
                     +ri(20)*(xy21*xy11+xz21*xz11) &
                     +ri(21)*(yy21*zz11+zz21*yy11) &
                     +ri(22)*yz21*yz11
                !(py px/py s)
                w(44)=x(2)*(ri(8)*xx21+ri(9)*yyzz21) &
                     +ri(10)*(y(2)*xy21+z(2)*xz21)
                !(py px/py px)
                w(45) = &
                     (ri(16)*xx21+ri(17)*yyzz21)*xx21+ri(18)*xx21*yyzz21 &
                     +ri(19)*(yy21*yy21+zz21*zz21) &
                     +ri(20)*(xy21*xy21+xz21*xz21) &
                     +ri(21)*(yy21*zz21+zz21*yy21) &
                     +ri(22)*yz21*yz21
                !(py px/py py)
                w(46) = &
                     (ri(16)*xx21+ri(17)*yyzz21)*xx22+ri(18)*xx21*yyzz22 &
                     +ri(19)*(yy21*yy22+zz21*zz22) &
                     +ri(20)*(xy21*xy22+xz21*xz22) &
                     +ri(21)*(yy21*zz22+zz21*yy22) &
                     +ri(22)*yz21*yz22
                !(py px/pz s)
                w(47)=x(3)*(ri(8)*xx21+ri(9)*yyzz21) &
                     +ri(10)*(         +z(3)*xz21)
                !(py px/pz px)
                w(48) = &
                     (ri(16)*xx21+ri(17)*yyzz21)*xx31 &
                     +(ri(18)*xx21+ri(19)*zz21+ri(21)*yy21)*zz31 &
                     +ri(20)*(xy21*xy31+xz21*xz31) &
                     +ri(22)*yz21*yz31
                !(py px/pz py)
                w(49) = &
                     (ri(16)*xx21+ri(17)*yyzz21)*xx32 &
                     +(ri(18)*xx21+ri(19)*zz21+ri(21)*yy21)*zz32 &
                     +ri(20)*(xy21*xy32+xz21*xz32) &
                     +ri(22)*yz21*yz32
                !(py px/pz pz)
                w(50) = &
                     (ri(16)*xx21+ri(17)*yyzz21)*xx33 &
                     +(ri(18)*xx21+ri(19)*zz21+ri(21)*yy21)*zz33 &
                     +ri(20)*xz21*xz33
                !(py py/s s)
                w(51)=ri(3)*xx22+ri(4)*yyzz22
                !(py py/px s)
                w(52)=x(1)*(ri(8)*xx22+ri(9)*yyzz22) &
                     +ri(10)*(y(1)*xy22+z(1)*xz22)
                !(py py/px px)
                w(53) = &
                     (ri(16)*xx22+ri(17)*yyzz22)*xx11+ri(18)*xx22*yyzz11 &
                     +ri(19)*(yy22*yy11+zz22*zz11) &
                     +ri(20)*(xy22*xy11+xz22*xz11) &
                     +ri(21)*(yy22*zz11+zz22*yy11) &
                     +ri(22)*yz22*yz11
                !(py py/py s)
                w(54)=x(2)*(ri(8)*xx22+ri(9)*yyzz22) &
                     +ri(10)*(y(2)*xy22+z(2)*xz22)
                !(py py/py px)
                w(55) = &
                     (ri(16)*xx22+ri(17)*yyzz22)*xx21+ri(18)*xx22*yyzz21 &
                     +ri(19)*(yy22*yy21+zz22*zz21) &
                     +ri(20)*(xy22*xy21+xz22*xz21) &
                     +ri(21)*(yy22*zz21+zz22*yy21) &
                     +ri(22)*yz22*yz21
                !(py py/py py)
                w(56) = &
                     (ri(16)*xx22+ri(17)*yyzz22)*xx22+ri(18)*xx22*yyzz22 &
                     +ri(19)*(yy22*yy22+zz22*zz22) &
                     +ri(20)*(xy22*xy22+xz22*xz22) &
                     +ri(21)*(yy22*zz22+zz22*yy22) &
                     +ri(22)*yz22*yz22
                !(py py/pz s)
                w(57)=x(3)*(ri(8)*xx22+ri(9)*yyzz22) &
                     +ri(10)*(         +z(3)*xz22)
                !(py py/pz px)
                w(58) = &
                     (ri(16)*xx22+ri(17)*yyzz22)*xx31 &
                     +(ri(18)*xx22+ri(19)*zz22+ri(21)*yy22)*zz31 &
                     +ri(20)*(xy22*xy31+xz22*xz31) &
                     +ri(22)*yz22*yz31
                !(py py/pz py)
                w(59) = &
                     (ri(16)*xx22+ri(17)*yyzz22)*xx32 &
                     +(ri(18)*xx22+ri(19)*zz22+ri(21)*yy22)*zz32 &
                     +ri(20)*(xy22*xy32+xz22*xz32) &
                     +ri(22)*yz22*yz32
                !(py py/pz pz)
                w(60) = &
                     (ri(16)*xx22+ri(17)*yyzz22)*xx33 &
                     +(ri(18)*xx22+ri(19)*zz22+ri(21)*yy22)*zz33 &
                     +ri(20)*xz22*xz33
                !(pz s/ss)
                w(61)=ri(2)*x(3)
                !(pz s/px s)
                w(62)=ri(6)*xx31+ri(7)*zz31
                !(pz s/px px)
                w(63)=x(3)*(ri(13)*xx11+ri(14)*yyzz11) &
                     +ri(15)*(         +z(3)*xz11)
                !(pz s/py s)
                w(64)=ri(6)*xx32+ri(7)*zz32
                !(pz s/py px)
                w(65)=x(3)*(ri(13)*xx21+ri(14)*yyzz21) &
                     +ri(15)*(         +z(3)*xz21)
                !(pz s/py py)
                w(66)=x(3)*(ri(13)*xx22+ri(14)*yyzz22) &
                     +ri(15)*(         +z(3)*xz22)
                !(pz s/pz s)
                w(67)=ri(6)*xx33+ri(7)*zz33
                !(pz s/pz px)
                w(68)=x(3)*(ri(13)*xx31+ri(14)*zz31) &
                     +ri(15)*(         +z(3)*xz31)
                !(pz s/pz py)
                w(69)=x(3)*(ri(13)*xx32+ri(14)*zz32) &
                     +ri(15)*(         +z(3)*xz32)
                !(pz s/pz pz)
                w(70)=x(3)*(ri(13)*xx33+ri(14)*zz33) &
                     +ri(15)*(         +z(3)*xz33)
                !(pz px/s s)
                w(71)=ri(3)*xx31+ri(4)*zz31
                !(pz px/px s)
                w(72)=x(1)*(ri(8)*xx31+ri(9)*zz31) &
                     +ri(10)*(y(1)*xy31+z(1)*xz31)
                !(pz px/px px)
                w(73) = &
                     (ri(16)*xx31+ri(17)*zz31)*xx11+ri(18)*xx31*yyzz11 &
                     +ri(19)*zz31*zz11 &
                     +ri(20)*(xy31*xy11+xz31*xz11) &
                     +ri(21)*zz31*yy11 &
                     +ri(22)*yz31*yz11
                !(pz px/py s)
                w(74)=x(2)*(ri(8)*xx31+ri(9)*zz31) &
                     +ri(10)*(y(2)*xy31+z(2)*xz31)
                !(pz px/py px)
                w(75) = &
                     (ri(16)*xx31+ri(17)*zz31)*xx21+ri(18)*xx31*yyzz21 &
                     +ri(19)*zz31*zz21 &
                     +ri(20)*(xy31*xy21+xz31*xz21) &
                     +ri(21)*zz31*yy21 &
                     +ri(22)*yz31*yz21
                !(pz px/py py)
                w(76) = &
                     (ri(16)*xx31+ri(17)*zz31)*xx22+ri(18)*xx31*yyzz22 &
                     +ri(19)*zz31*zz22 &
                     +ri(20)*(xy31*xy22+xz31*xz22) &
                     +ri(21)*zz31*yy22 &
                     +ri(22)*yz31*yz22
                !(pz px/pz s)
                w(77)=x(3)*(ri(8)*xx31+ri(9)*zz31) &
                     +ri(10)*(         +z(3)*xz31)
                !(pz px/pz px)
                w(78) = &
                     (ri(16)*xx31+ri(17)*zz31)*xx31 &
                     +(ri(18)*xx31+ri(19)*zz31)*zz31 &
                     +ri(20)*(xy31*xy31+xz31*xz31) &
                     +ri(22)*yz31*yz31
                !(pz px/pz py)
                w(79) = &
                     (ri(16)*xx31+ri(17)*zz31)*xx32 &
                     +(ri(18)*xx31+ri(19)*zz31)*zz32 &
                     +ri(20)*(xy31*xy32+xz31*xz32) &
                     +ri(22)*yz31*yz32
                !(pz px/pz pz)
                w(80) = &
                     (ri(16)*xx31+ri(17)*zz31)*xx33 &
                     +(ri(18)*xx31+ri(19)*zz31)*zz33 &
                     +ri(20)*xz31*xz33
                !(pz py/s s)
                w(81)=ri(3)*xx32+ri(4)*zz32
                !(pz py/px s)
                w(82)=x(1)*(ri(8)*xx32+ri(9)*zz32) &
                     +ri(10)*(y(1)*xy32+z(1)*xz32)
                !(pz py/px px)
                w(83) = &
                     (ri(16)*xx32+ri(17)*zz32)*xx11+ri(18)*xx32*yyzz11 &
                     +ri(19)*zz32*zz11 &
                     +ri(20)*(xy32*xy11+xz32*xz11) &
                     +ri(21)*zz32*yy11 &
                     +ri(22)*yz32*yz11
                !(pz py/py s)
                w(84)=x(2)*(ri(8)*xx32+ri(9)*zz32) &
                     +ri(10)*(y(2)*xy32+z(2)*xz32)
                !(pz py/py px)
                w(85) = &
                     (ri(16)*xx32+ri(17)*zz32)*xx21+ri(18)*xx32*yyzz21 &
                     +ri(19)*zz32*zz21 &
                     +ri(20)*(xy32*xy21+xz32*xz21) &
                     +ri(21)*zz32*yy21 &
                     +ri(22)*yz32*yz21
                !(pz py/py py)
                w(86) = &
                     (ri(16)*xx32+ri(17)*zz32)*xx22+ri(18)*xx32*yyzz22 &
                     +ri(19)*zz32*zz22 &
                     +ri(20)*(xy32*xy22+xz32*xz22) &
                     +ri(21)*zz32*yy22 &
                     +ri(22)*yz32*yz22
                !(pz py/pz s)
                w(87)=x(3)*(ri(8)*xx32+ri(9)*zz32) &
                     +ri(10)*(         +z(3)*xz32)
                !(pz py/pz px)
                w(88) = &
                     (ri(16)*xx32+ri(17)*zz32)*xx31 &
                     +(ri(18)*xx32+ri(19)*zz32)*zz31 &
                     +ri(20)*(xy32*xy31+xz32*xz31) &
                     +ri(22)*yz32*yz31
                !(pz py/pz py)
                w(89) = &
                     (ri(16)*xx32+ri(17)*zz32)*xx32 &
                     +(ri(18)*xx32+ri(19)*zz32)*zz32 &
                     +ri(20)*(xy32*xy32+xz32*xz32) &
                     +ri(22)*yz32*yz32
                !(pz py/pz pz)
                w(90) = &
                     (ri(16)*xx32+ri(17)*zz32)*xx33 &
                     +(ri(18)*xx32+ri(19)*zz32)*zz33 &
                     +ri(20)*xz32*xz33
                !(pz pz/s s)
                w(91)=ri(3)*xx33+ri(4)*zz33
                !(pz pz/px s)
                w(92)=x(1)*(ri(8)*xx33+ri(9)*zz33) &
                     +ri(10)*(          z(1)*xz33)
                !(pz pz/px px)
                w(93) = &
                     (ri(16)*xx33+ri(17)*zz33)*xx11+ri(18)*xx33*yyzz11 &
                     +ri(19)*zz33*zz11 &
                     +ri(20)*xz33*xz11 &
                     +ri(21)*zz33*yy11
                !(pz pz/py s)
                w(94)=x(2)*(ri(8)*xx33+ri(9)*zz33) &
                     +ri(10)*(         +z(2)*xz33)
                !(pz pz/py px)
                w(95) = &
                     (ri(16)*xx33+ri(17)*zz33)*xx21+ri(18)*xx33*yyzz21 &
                     +ri(19)*zz33*zz21 &
                     +ri(20)*xz33*xz21 &
                     +ri(21)*zz33*yy21
                !(pz pz/py py)
                w(96) = &
                     (ri(16)*xx33+ri(17)*zz33)*xx22+ri(18)*xx33*yyzz22 &
                     +ri(19)*zz33*zz22 &
                     +ri(20)*xz33*xz22 &
                     +ri(21)*zz33*yy22
                !(pz pz/pz s)
                w(97)=x(3)*(ri(8)*xx33+ri(9)*zz33) &
                     +ri(10)*(         +z(3)*xz33)
                !(pz pz/pz px)
                w(98) = &
                     (ri(16)*xx33+ri(17)*zz33)*xx31 &
                     +(ri(18)*xx33+ri(19)*zz33)*zz31 &
                     +ri(20)*xz33*xz31
                !(pz pz/pz py)
                w(99) = &
                     (ri(16)*xx33+ri(17)*zz33)*xx32 &
                     +(ri(18)*xx33+ri(19)*zz33)*zz32 &
                     +ri(20)*xz33*xz32
                !(pz pz/pz pz)
                w(100) = &
                     (ri(16)*xx33+ri(17)*zz33)*xx33 &
                     +(ri(18)*xx33+ri(19)*zz33)*zz33 &
                     +ri(20)*xz33*xz33
             ELSE
                !(px s/s s)
                w(2)=ri(2)*x(1)
                !(px px/s s)
                w(3)=ri(3)*xx11+ri(4)*yyzz11
                !(py s/s s)
                w(4)=ri(2)*x(2)
                !(py px/s s)
                w(5)=ri(3)*xx21+ri(4)*yyzz21
                !(py py/s s)
                w(6)=ri(3)*xx22+ri(4)*yyzz22
                !(pz s/ss)
                w(7)=ri(2)*x(3)
                !(pz px/s s)
                w(8)=ri(3)*xx31+ri(4)*zz31
                !(pz py/s s)
                w(9)=ri(3)*xx32+ri(4)*zz32
                !(pz pz/s s)
                w(10)=ri(3)*xx33+ri(4)*zz33
             END IF
          END IF
          IF (invert) CALL invert_integral(w,42)
          IF (debug_this_module) THEN
             ! Check value of integrals
             w2=0.0_dp
             CALL rotint (sepi,sepj,rijv,w2,shortrange)
             DO J = 1, 100
                IF (ABS(w(j))>1.0E-6_dp) THEN
                   IF (ABS(w2(j)-w(j))>1.0E-7) THEN
                      WRITE(*,*)"check integral W",j,ABS(w2(j)-w(j))
                      WRITE(*,'(10F12.6)')w
                      WRITE(*,*)
                      WRITE(*,'(10F12.6)')w2
                      WRITE(*,*)
                      WRITE(*,'(10F12.6)')w2-w
                      WRITE(*,*)
                      WRITE(*,*)invert
                      STOP
                   END IF
                END IF
             END DO
          END IF
       END IF
       ! Gradients if requested
       IF (lgrad) THEN
          !(s s/s s)
          dw(1,:)=dri(1,:)
          IF (sj) THEN
             !(s s/px s)
             dw(2,:)=dri(5,:)* x(1)+ ri(5)*dx(1,:)
             !(s s/px px)
             dw(3,:)=dri(11,:)* xx11+ ri(11)*dxx11+&
                   dri(12,:)* yyzz11+ri(12)*dyyzz11
             !(s s/py s)
             dw(4,:)=dri(5,:)* x(2)+ ri(5)*dx(2,:)
             !(s s/py px)
             dw(5,:)=dri(11,:)* xx21+ri(11)*dxx21+&
                   dri(12,:)* yyzz21+ri(12)*dyyzz21
             !(s s/py py)
             dw(6,:)=dri(11,:)* xx22+ri(11)*dxx22+&
                   dri(12,:)* yyzz22+ri(12)*dyyzz22
             !(s s/pz s)
             dw(7,:)=dri(5,:)* x(3)+ri(5)*dx(3,:)
             !(s s/pz px)
             dw(8,:)=dri(11,:)* xx31+ri(11)*dxx31+&
                   dri(12,:)* zz31+ri(12)*dzz31
             !(s s/pz py)
             dw(9,:)=dri(11,:)* xx32+ri(11)*dxx32+&
                   dri(12,:)* zz32+ri(12)*dzz32
             !(s s/pz pz)
             dw(10,:)=dri(11,:)* xx33+ri(11)*dxx33+&
                    dri(12,:)* zz33+ri(12)*dzz33
          END IF

          IF (si) THEN
             !(px s/s s)
             dw(11,:)=dri(2,:)* x(1)+ri(2)*dx(1,:)
             IF (sj) THEN
                !(px s/px s)
                dw(12,:)=dri(6,:)* xx11+ri(6)*dxx11+&
                       dri(7,:)* yyzz11+ri(7)*dyyzz11
                !(px s/px px)
                dw(13,:)=dx(1,:)* (ri(13)*xx11+ri(14)*yyzz11)+&
                        x(1)*(dri(13,:)* xx11+ri(13)*dxx11+dri(14,:)* yyzz11+ri(14)*dyyzz11)+&
                       dri(15,:)* (y(1)*xy11+z(1)*xz11)+&
                        ri(15)*(dy(1,:)* xy11+y(1)*dxy11+dz(1,:)* xz11+z(1)*dxz11)
                !(px s/py s)
                dw(14,:)=dri(6,:)* xx21+ri(6)*dxx21+&
                       dri(7,:)* yyzz21+ri(7)*dyyzz21
                !(px s/py px)
                dw(15,:)=dx(1,:)*(ri(13)*xx21+ri(14)*yyzz21)+&
                        x(1)*(dri(13,:)*xx21+ri(13)*dxx21+dri(14,:)*yyzz21+ri(14)*dyyzz21)+&
                       dri(15,:)*(y(1)*xy21+z(1)*xz21)+&
                        ri(15)*(dy(1,:)*xy21+y(1)*dxy21+dz(1,:)*xz21+z(1)*dxz21)
                !(px s/py py)
                dw(16,:)=dx(1,:)*(ri(13)*xx22+ri(14)*yyzz22)+&
                        x(1)*(dri(13,:)*xx22+ri(13)*dxx22+dri(14,:)*yyzz22+ri(14)*dyyzz22)+&
                       dri(15,:)*(y(1)*xy22+z(1)*xz22)+&
                        ri(15)*(dy(1,:)*xy22+y(1)*dxy22+dz(1,:)*xz22+z(1)*dxz22)
                !(px s/pz s)
                dw(17,:)=dri(6,:)*xx31+ri(6)*dxx31+dri(7,:)*zz31+ri(7)*dzz31
                !(px s/pz px)
                dw(18,:)=dx(1,:)*(ri(13)*xx31+ri(14)*zz31)+&
                          x(1)*(dri(13,:)*xx31+ri(13)*dxx31+dri(14,:)*zz31+ri(14)*dzz31) &
                        +dri(15,:)*(y(1)*xy31+z(1)*xz31)+&
                          ri(15)*(dy(1,:)*xy31+y(1)*dxy31+dz(1,:)*xz31+z(1)*dxz31)
                !(px s/pz py)
                dw(19,:)=dx(1,:)*(ri(13)*xx32+ri(14)*zz32)+&
                          x(1)*(dri(13,:)*xx32+ri(13)*dxx32+dri(14,:)*zz32+ri(14)*dzz32) &
                        +dri(15,:)*(y(1)*xy32+z(1)*xz32)+&
                        ri(15)*(dy(1,:)*xy32+y(1)*dxy32+dz(1,:)*xz32+z(1)*dxz32)
                !(px s/pz pz)
                dw(20,:)=dx(1,:)*(ri(13)*xx33+ri(14)*zz33)+&
                          x(1)*(dri(13,:)*xx33+ri(13)*dxx33+dri(14,:)*zz33+ri(14)*dzz33) &
                          +dri(15,:)*(z(1)*xz33)+ri(15)*(dz(1,:)*xz33+z(1)*dxz33)
                !(px px/s s)
                dw(21,:)=dri(3,:)*xx11+ri(3)*dxx11+dri(4,:)*yyzz11+ri(4)*dyyzz11
                !(px px/px s)
                dw(22,:)=dx(1,:)*(ri(8)*xx11+ri(9)*yyzz11)+&
                     x(1)*(dri(8,:)*xx11+ri(8)*dxx11+dri(9,:)*yyzz11+ri(9)*dyyzz11) &
                     +dri(10,:)*(y(1)*xy11+z(1)*xz11)+ri(10)*(dy(1,:)*xy11+y(1)*dxy11+dz(1,:)*xz11+z(1)*dxz11)
                !(px px/px px)
                dw(23,:) = &
                     (dri(16,:)*xx11+ri(16)*dxx11+dri(17,:)*yyzz11+ri(17)*dyyzz11)*xx11+(ri(16)*xx11+ri(17)*yyzz11)*dxx11&
                     +dri(18,:)*xx11*yyzz11+ri(18)*dxx11*yyzz11+ri(18)*xx11*dyyzz11 &
                     +dri(19,:)*(yy11*yy11+zz11*zz11)+ri(19)*(2.0_dp*dyy11*yy11+2.0_dp*dzz11*zz11) &
                     +dri(20,:)*(xy11*xy11+xz11*xz11)+ri(20)*(2.0_dp*dxy11*xy11+2.0_dp*dxz11*xz11) &
                     +dri(21,:)*(yy11*zz11+zz11*yy11)+ri(21)*(2.0_dp*dyy11*zz11+2.0_dp*dzz11*yy11) &
                     +dri(22,:)*yz11*yz11+2.0_dp*ri(22)*dyz11*yz11
                !(px px/py s)
                dw(24,:)=dx(2,:)*(ri(8)*xx11+ri(9)*yyzz11)+x(2)*(dri(8,:)*xx11+ri(8)*dxx11+dri(9,:)*yyzz11+ri(9)*dyyzz11) &
                     +dri(10,:)*(y(2)*xy11+z(2)*xz11)+ri(10)*(dy(2,:)*xy11+y(2)*dxy11+dz(2,:)*xz11+z(2)*dxz11)
                !(px px/py px)
                dw(25,:) = &
                     (dri(16,:)*xx11+ri(16)*dxx11+dri(17,:)*yyzz11+ri(17)*dyyzz11)*xx21+(ri(16)*xx11+ri(17)*yyzz11)*dxx21&
                     +dri(18,:)*xx11*yyzz21+ri(18)*dxx11*yyzz21+ri(18)*xx11*dyyzz21 &
                     +dri(19,:)*(yy11*yy21+zz11*zz21)+ri(19)*(dyy11*yy21+yy11*dyy21+dzz11*zz21+zz11*dzz21) &
                     +dri(20,:)*(xy11*xy21+xz11*xz21)+ri(20)*(dxy11*xy21+xy11*dxy21+dxz11*xz21+xz11*dxz21) &
                     +dri(21,:)*(yy11*zz21+zz11*yy21)+ri(21)*(dyy11*zz21+yy11*dzz21+dzz11*yy21+zz11*dyy21) &
                     +dri(22,:)*yz11*yz21+ri(22)*dyz11*yz21+ri(22)*yz11*dyz21
                !(px px/py py)
                dw(26,:) = &
                     (dri(16,:)*xx11+ri(16)*dxx11+dri(17,:)*yyzz11+ri(17)*dyyzz11)*xx22+(ri(16)*xx11+ri(17)*yyzz11)*dxx22&
                     +dri(18,:)*xx11*yyzz22+ri(18)*dxx11*yyzz22+ri(18)*xx11*dyyzz22 &
                     +dri(19,:)*(yy11*yy22+zz11*zz22)+ri(19)*(dyy11*yy22+yy11*dyy22+dzz11*zz22+zz11*dzz22) &
                     +dri(20,:)*(xy11*xy22+xz11*xz22)+ri(20)*(dxy11*xy22+xy11*dxy22+dxz11*xz22+xz11*dxz22) &
                     +dri(21,:)*(yy11*zz22+zz11*yy22)+ri(21)*(dyy11*zz22+yy11*dzz22+dzz11*yy22+zz11*dyy22) &
                     +dri(22,:)*yz11*yz22+ri(22)*dyz11*yz22+ri(22)*yz11*dyz22
                !(px px/pz s)
                dw(27,:)=dx(3,:)*(ri(8)*xx11+ri(9)*yyzz11)+x(3)*(dri(8,:)*xx11+ri(8)*dxx11+dri(9,:)*yyzz11+ri(9)*dyyzz11) &
                     +dri(10,:)*(z(3)*xz11)+ri(10)*(dz(3,:)*xz11+z(3)*dxz11)
                !(px px/pz px)
                dw(28,:) = &
                     (dri(16,:)*xx11+ri(16)*dxx11+dri(17,:)*yyzz11+ri(17)*dyyzz11)*xx31+(ri(16)*xx11+ri(17)*yyzz11)*dxx31 &
                     +(dri(18,:)*xx11+ri(18)*dxx11+dri(19,:)*zz11+ri(19)*dzz11+dri(21,:)*yy11+ri(21)*dyy11)*zz31+&
                     (ri(18)*xx11+ri(19)*zz11+ri(21)*yy11)*dzz31 &
                     +dri(20,:)*(xy11*xy31+xz11*xz31)+ri(20)*(dxy11*xy31+xy11*dxy31+dxz11*xz31+xz11*dxz31) &
                     +dri(22,:)*yz11*yz31+ri(22)*dyz11*yz31+ri(22)*yz11*dyz31
                !(px px/pz py)
                dw(29,:) = &
                     (dri(16,:)*xx11+ri(16)*dxx11+dri(17,:)*yyzz11+ri(17)*dyyzz11)*xx32+(ri(16)*xx11+ri(17)*yyzz11)*dxx32 &
                     +(dri(18,:)*xx11+ri(18)*dxx11+dri(19,:)*zz11+ri(19)*dzz11+dri(21,:)*yy11+ri(21)*dyy11)*zz32+&
                     (ri(18)*xx11+ri(19)*zz11+ri(21)*yy11)*dzz32 &
                     +dri(20,:)*(xy11*xy32+xz11*xz32)+ri(20)*(dxy11*xy32+xy11*dxy32+dxz11*xz32+xz11*dxz32) &
                     +dri(22,:)*yz11*yz32+ri(22)*dyz11*yz32+ri(22)*yz11*dyz32
                !(px px/pz pz)
                dw(30,:) = &
                     (dri(16,:)*xx11+ri(16)*dxx11+dri(17,:)*yyzz11+ri(17)*dyyzz11)*xx33+(ri(16)*xx11+ri(17)*yyzz11)*dxx33 &
                     +(dri(18,:)*xx11+ri(18)*dxx11+dri(19,:)*zz11+ri(19)*dzz11+dri(21,:)*yy11+ri(21)*dyy11)*zz33+&
                     (ri(18)*xx11+ri(19)*zz11+ri(21)*yy11)*dzz33 &
                     +dri(20,:)*xz11*xz33+ri(20)*dxz11*xz33+ri(20)*xz11*dxz33
                !(py s/s s)
                dw(31,:)=dri(2,:)*x(2)+ri(2)*dx(2,:)
                !(py s/px s)
                dw(32,:)=dri(6,:)*xx21+ri(6)*dxx21+&
                     dri(7,:)*yyzz21+ri(7)*dyyzz21
                !(py s/px px)
                dw(33,:)=dx(2,:)*(ri(13)*xx11+ri(14)*yyzz11)+&
                     x(2)*(dri(13,:)*xx11+ri(13)*dxx11+dri(14,:)*yyzz11+ri(14)*dyyzz11) &
                     +dri(15,:)*(y(2)*xy11+z(2)*xz11)+&
                     ri(15)*(dy(2,:)*xy11+y(2)*dxy11+dz(2,:)*xz11+z(2)*dxz11)
                !(py s/py s)
                dw(34,:)=dri(6,:)*xx22+ri(6)*dxx22+&
                     dri(7,:)*yyzz22+ri(7)*dyyzz22
                !(py s/py px)
                dw(35,:)=dx(2,:)*(ri(13)*xx21+ri(14)*yyzz21)+x(2)*(dri(13,:)*xx21+ri(13)*dxx21+dri(14,:)*yyzz21+ri(14)*dyyzz21) &
                     +dri(15,:)*(y(2)*xy21+z(2)*xz21)+ri(15)*(dy(2,:)*xy21+y(2)*dxy21+dz(2,:)*xz21+z(2)*dxz21)
                !(py s/py py)
                dw(36,:)=dx(2,:)*(ri(13)*xx22+ri(14)*yyzz22)+x(2)*(dri(13,:)*xx22+ri(13)*dxx22+dri(14,:)*yyzz22+ri(14)*dyyzz22) &
                     +dri(15,:)*(y(2)*xy22+z(2)*xz22)+ri(15)*(dy(2,:)*xy22+y(2)*dxy22+dz(2,:)*xz22+z(2)*dxz22)
                !(py s/pz s)
                dw(37,:)=dri(6,:)*xx32+ri(6)*dxx32+&
                     dri(7,:)*zz32+ri(7)*dzz32
                !(py s/pz px)
                dw(38,:)=dx(2,:)*(ri(13)*xx31+ri(14)*zz31)+x(2)*(dri(13,:)*xx31+ri(13)*dxx31+dri(14,:)*zz31+ri(14)*dzz31) &
                     +dri(15,:)*(y(2)*xy31+z(2)*xz31)+ri(15)*(dy(2,:)*xy31+y(2)*dxy31+dz(2,:)*xz31+z(2)*dxz31)
                !(py s/pz py)
                dw(39,:)=dx(2,:)*(ri(13)*xx32+ri(14)*zz32)+x(2)*(dri(13,:)*xx32+ri(13)*dxx32+dri(14,:)*zz32+ri(14)*dzz32) &
                     +dri(15,:)*(y(2)*xy32+z(2)*xz32)+ri(15)*(dy(2,:)*xy32+y(2)*dxy32+dz(2,:)*xz32+z(2)*dxz32)
                !(py s/pz pz)
                dw(40,:)=dx(2,:)*(ri(13)*xx33+ri(14)*zz33)+x(2)*(dri(13,:)*xx33+ri(13)*dxx33+dri(14,:)*zz33+ri(14)*dzz33) &
                     +dri(15,:)*(z(2)*xz33)+ri(15)*(dz(2,:)*xz33+z(2)*dxz33)
                !(py px/s s)
                dw(41,:)=dri(3,:)*xx21+ri(3)*dxx21+&
                     dri(4,:)*yyzz21+ri(4)*dyyzz21
                !(py px/px s)
                dw(42,:)=dx(1,:)*(ri(8)*xx21+ri(9)*yyzz21)+x(1)*(dri(8,:)*xx21+ri(8)*dxx21+dri(9,:)*yyzz21+ri(9)*dyyzz21) &
                     +dri(10,:)*(y(1)*xy21+z(1)*xz21)+ri(10)*(dy(1,:)*xy21+y(1)*dxy21+dz(1,:)*xz21+z(1)*dxz21)
                !(py px/px px)
                dw(43,:) = &
                     (dri(16,:)*xx21+ri(16)*dxx21+dri(17,:)*yyzz21+ri(17)*dyyzz21)*xx11+&
                     (ri(16)*xx21+ri(17)*yyzz21)*dxx11&
                     +dri(18,:)*xx21*yyzz11+ri(18)*dxx21*yyzz11+ri(18)*xx21*dyyzz11 &
                     +dri(19,:)*(yy21*yy11+zz21*zz11)+ri(19)*(dyy21*yy11+yy21*dyy11+dzz21*zz11+zz21*dzz11) &
                     +dri(20,:)*(xy21*xy11+xz21*xz11)+ri(20)*(dxy21*xy11+xy21*dxy11+dxz21*xz11+xz21*dxz11) &
                     +dri(21,:)*(yy21*zz11+zz21*yy11)+ri(21)*(dyy21*zz11+yy21*dzz11+dzz21*yy11+zz21*dyy11) &
                     +dri(22,:)*yz21*yz11+ri(22)*dyz21*yz11+ri(22)*yz21*dyz11
                !(py px/py s)
                dw(44,:)=dx(2,:)*(ri(8)*xx21+ri(9)*yyzz21)+x(2)*(dri(8,:)*xx21+ri(8)*dxx21+dri(9,:)*yyzz21+ri(9)*dyyzz21) &
                     +dri(10,:)*(y(2)*xy21+z(2)*xz21)+ri(10)*(dy(2,:)*xy21+y(2)*dxy21+dz(2,:)*xz21+z(2)*dxz21)
                !(py px/py px)
                dw(45,:) = &
                     (dri(16,:)*xx21+ri(16)*dxx21+dri(17,:)*yyzz21+ri(17)*dyyzz21)*xx21+(ri(16)*xx21+ri(17)*yyzz21)*dxx21&
                     +dri(18,:)*xx21*yyzz21+ri(18)*dxx21*yyzz21+ri(18)*xx21*dyyzz21 &
                     +dri(19,:)*(yy21*yy21+zz21*zz21)+ri(19)*(dyy21*yy21+yy21*dyy21+dzz21*zz21+zz21*dzz21) &
                     +dri(20,:)*(xy21*xy21+xz21*xz21)+ri(20)*(dxy21*xy21+xy21*dxy21+dxz21*xz21+xz21*dxz21) &
                     +dri(21,:)*(yy21*zz21+zz21*yy21)+ri(21)*(dyy21*zz21+yy21*dzz21+dzz21*yy21+zz21*dyy21) &
                     +dri(22,:)*yz21*yz21 +ri(22)*dyz21*yz21 +ri(22)*yz21*dyz21
                !(py px/py py)
                dw(46,:) = &
                     (dri(16,:)*xx21+ri(16)*dxx21+dri(17,:)*yyzz21+ri(17)*dyyzz21)*xx22+(ri(16)*xx21+ri(17)*yyzz21)*dxx22&
                     +dri(18,:)*xx21*yyzz22+ri(18)*dxx21*yyzz22+ri(18)*xx21*dyyzz22 &
                     +dri(19,:)*(yy21*yy22+zz21*zz22)+ri(19)*(dyy21*yy22+yy21*dyy22+dzz21*zz22+zz21*dzz22) &
                     +dri(20,:)*(xy21*xy22+xz21*xz22)+ri(20)*(dxy21*xy22+xy21*dxy22+dxz21*xz22+xz21*dxz22) &
                     +dri(21,:)*(yy21*zz22+zz21*yy22)+ri(21)*(dyy21*zz22+yy21*dzz22+dzz21*yy22+zz21*dyy22) &
                     +dri(22,:)*yz21*yz22+ri(22)*dyz21*yz22+ri(22)*yz21*dyz22
                !(py px/pz s)
                dw(47,:)=dx(3,:)*(ri(8)*xx21+ri(9)*yyzz21)+x(3)*(dri(8,:)*xx21+ri(8)*dxx21+dri(9,:)*yyzz21+ri(9)*dyyzz21) &
                     +dri(10,:)*(z(3)*xz21)+ri(10)*(dz(3,:)*xz21+z(3)*dxz21)
                !(py px/pz px)
                dw(48,:) = &
                     (dri(16,:)*xx21+ri(16)*dxx21+dri(17,:)*yyzz21+ri(17)*dyyzz21)*xx31+(ri(16)*xx21+ri(17)*yyzz21)*dxx31 &
                     +(dri(18,:)*xx21+ri(18)*dxx21+dri(19,:)*zz21+ri(19)*dzz21+dri(21,:)*yy21+ri(21)*dyy21)*zz31 &
                     +(ri(18)*xx21+ri(19)*zz21+ri(21)*yy21)*dzz31 &
                     +dri(20,:)*(xy21*xy31+xz21*xz31)+ri(20)*(dxy21*xy31+xy21*dxy31+dxz21*xz31+xz21*dxz31) &
                     +dri(22,:)*yz21*yz31 +ri(22)*dyz21*yz31 +ri(22)*yz21*dyz31
                !(py px/pz py)
                dw(49,:) = &
                     (dri(16,:)*xx21+ri(16)*dxx21+dri(17,:)*yyzz21+ri(17)*dyyzz21)*xx32+(ri(16)*xx21+ri(17)*yyzz21)*dxx32 &
                     +(dri(18,:)*xx21+ri(18)*dxx21+dri(19,:)*zz21+ri(19)*dzz21+dri(21,:)*yy21+ri(21)*dyy21)*zz32&
                     +(ri(18)*xx21+ri(19)*zz21+ri(21)*yy21)*dzz32 &
                     +dri(20,:)*(xy21*xy32+xz21*xz32)+ri(20)*(dxy21*xy32+xy21*dxy32+dxz21*xz32+xz21*dxz32) &
                     +dri(22,:)*yz21*yz32+ri(22)*dyz21*yz32+ri(22)*yz21*dyz32
                !(py px/pz pz)
                dw(50,:) = &
                     (dri(16,:)*xx21+ri(16)*dxx21+dri(17,:)*yyzz21+ri(17)*dyyzz21)*xx33&
                     +(ri(16)*xx21+ri(17)*yyzz21)*dxx33 &
                     +(dri(18,:)*xx21+ri(18)*dxx21+dri(19,:)*zz21+ri(19)*dzz21+dri(21,:)*yy21+ri(21)*dyy21)*zz33&
                     +(ri(18)*xx21+ri(19)*zz21+ri(21)*yy21)*dzz33 &
                     +dri(20,:)*xz21*xz33+ri(20)*dxz21*xz33+ri(20)*xz21*dxz33
                !(py py/s s)
                dw(51,:)=dri(3,:)*xx22+ri(3)*dxx22+&
                     dri(4,:)*yyzz22+ri(4)*dyyzz22
                !(py py/px s)
                dw(52,:)=dx(1,:)*(ri(8)*xx22+ri(9)*yyzz22)+x(1)*(dri(8,:)*xx22+ri(8)*dxx22+dri(9,:)*yyzz22+ri(9)*dyyzz22) &
                     +dri(10,:)*(y(1)*xy22+z(1)*xz22)+ri(10)*(dy(1,:)*xy22+y(1)*dxy22+dz(1,:)*xz22+z(1)*dxz22)
                !(py py/px px)
                dw(53,:) = &
                     (dri(16,:)*xx22+ri(16)*dxx22+dri(17,:)*yyzz22+ri(17)*dyyzz22)*xx11&
                     +(ri(16)*xx22+ri(17)*yyzz22)*dxx11&
                     +dri(18,:)*xx22*yyzz11+ri(18)*dxx22*yyzz11+ri(18)*xx22*dyyzz11 &
                     +dri(19,:)*(yy22*yy11+zz22*zz11)+ri(19)*(dyy22*yy11+yy22*dyy11+dzz22*zz11+zz22*dzz11) &
                     +dri(20,:)*(xy22*xy11+xz22*xz11)+ri(20)*(dxy22*xy11+xy22*dxy11+dxz22*xz11+xz22*dxz11) &
                     +dri(21,:)*(yy22*zz11+zz22*yy11)+ri(21)*(dyy22*zz11+yy22*dzz11+dzz22*yy11+zz22*dyy11) &
                     +dri(22,:)*yz22*yz11+ri(22)*dyz22*yz11+ri(22)*yz22*dyz11
                !(py py/py s)
                dw(54,:)=dx(2,:)*(ri(8)*xx22+ri(9)*yyzz22)+x(2)*(dri(8,:)*xx22+ri(8)*dxx22+dri(9,:)*yyzz22+ri(9)*dyyzz22) &
                     +dri(10,:)*(y(2)*xy22+z(2)*xz22)+ri(10)*(dy(2,:)*xy22+y(2)*dxy22+dz(2,:)*xz22+z(2)*dxz22)
                !(py py/py px)
                dw(55,:) = &
                     (dri(16,:)*xx22+ri(16)*dxx22+dri(17,:)*yyzz22+ri(17)*dyyzz22)*xx21&
                     +(ri(16)*xx22+ri(17)*yyzz22)*dxx21&
                     +dri(18,:)*xx22*yyzz21+ri(18)*dxx22*yyzz21+ri(18)*xx22*dyyzz21 &
                     +dri(19,:)*(yy22*yy21+zz22*zz21)+ri(19)*(dyy22*yy21+yy22*dyy21+dzz22*zz21+zz22*dzz21) &
                     +dri(20,:)*(xy22*xy21+xz22*xz21)+ri(20)*(dxy22*xy21+xy22*dxy21+dxz22*xz21+xz22*dxz21) &
                     +dri(21,:)*(yy22*zz21+zz22*yy21)+ri(21)*(dyy22*zz21+yy22*dzz21+dzz22*yy21+zz22*dyy21) &
                     +dri(22,:)*yz22*yz21+ri(22)*dyz22*yz21+ri(22)*yz22*dyz21
                !(py py/py py)
                dw(56,:) = &
                     (dri(16,:)*xx22+ri(16)*dxx22+dri(17,:)*yyzz22+ri(17)*dyyzz22)*xx22&
                     +(ri(16)*xx22+ri(17)*yyzz22)*dxx22&
                     +dri(18,:)*xx22*yyzz22+ri(18)*dxx22*yyzz22+ri(18)*xx22*dyyzz22 &
                     +dri(19,:)*(yy22*yy22+zz22*zz22)+ri(19)*(dyy22*yy22+yy22*dyy22+dzz22*zz22+zz22*dzz22) &
                     +dri(20,:)*(xy22*xy22+xz22*xz22)+ri(20)*(dxy22*xy22+xy22*dxy22+dxz22*xz22+xz22*dxz22) &
                     +dri(21,:)*(yy22*zz22+zz22*yy22)+ri(21)*(dyy22*zz22+yy22*dzz22+dzz22*yy22+zz22*dyy22) &
                     +dri(22,:)*yz22*yz22+ri(22)*dyz22*yz22+ri(22)*yz22*dyz22
                !(py py/pz s)
                dw(57,:)=dx(3,:)*(ri(8)*xx22+ri(9)*yyzz22)+x(3)*(dri(8,:)*xx22+ri(8)*dxx22+dri(9,:)*yyzz22+ri(9)*dyyzz22) &
                     +dri(10,:)*(z(3)*xz22)+ri(10)*(dz(3,:)*xz22+z(3)*dxz22)
                !(py py/pz px)
                dw(58,:) = &
                     (dri(16,:)*xx22+ri(16)*dxx22+dri(17,:)*yyzz22+ri(17)*dyyzz22)*xx31&
                     +(ri(16)*xx22+ri(17)*yyzz22)*dxx31 &
                     +(dri(18,:)*xx22+ri(18)*dxx22+dri(19,:)*zz22+ri(19)*dzz22+dri(21,:)*yy22+ri(21)*dyy22)*zz31&
                     +(ri(18)*xx22+ri(19)*zz22+ri(21)*yy22)*dzz31 &
                     +dri(20,:)*(xy22*xy31+xz22*xz31)+ri(20)*(dxy22*xy31+xy22*dxy31+dxz22*xz31+xz22*dxz31) &
                     +dri(22,:)*yz22*yz31+ri(22)*dyz22*yz31+ri(22)*yz22*dyz31
                !(py py/pz py)
                dw(59,:) = &
                     (dri(16,:)*xx22+ri(16)*dxx22+dri(17,:)*yyzz22+ri(17)*dyyzz22)*xx32&
                     +(ri(16)*xx22+ri(17)*yyzz22)*dxx32 &
                     +(dri(18,:)*xx22+ri(18)*dxx22+dri(19,:)*zz22+ri(19)*dzz22+dri(21,:)*yy22+ri(21)*dyy22)*zz32&
                     +(ri(18)*xx22+ri(19)*zz22+ri(21)*yy22)*dzz32 &
                     +dri(20,:)*(xy22*xy32+xz22*xz32)+ri(20)*(dxy22*xy32+xy22*dxy32+dxz22*xz32+xz22*dxz32) &
                     +dri(22,:)*yz22*yz32+ri(22)*dyz22*yz32+ri(22)*yz22*dyz32
                !(py py/pz pz)
                dw(60,:) = &
                     (dri(16,:)*xx22+ri(16)*dxx22+dri(17,:)*yyzz22+ri(17)*dyyzz22)*xx33&
                     +(ri(16)*xx22+ri(17)*yyzz22)*dxx33 &
                     +(dri(18,:)*xx22+ri(18)*dxx22+dri(19,:)*zz22+ri(19)*dzz22+dri(21,:)*yy22+ri(21)*dyy22)*zz33&
                     +(ri(18)*xx22+ri(19)*zz22+ri(21)*yy22)*dzz33 &
                     +dri(20,:)*xz22*xz33+ri(20)*dxz22*xz33+ri(20)*xz22*dxz33
                !(pz s/ss)
                dw(61,:)=dri(2,:)*x(3)+ri(2)*dx(3,:)
                !(pz s/px s)
                dw(62,:)=dri(6,:)*xx31+ri(6)*dxx31+&
                     dri(7,:)*zz31+ri(7)*dzz31
                !(pz s/px px)
                dw(63,:)=dx(3,:)*(ri(13)*xx11+ri(14)*yyzz11)&
                     +x(3)*(dri(13,:)*xx11+ri(13)*dxx11+dri(14,:)*yyzz11+ri(14)*dyyzz11) &
                     +dri(15,:)*(z(3)*xz11)+ri(15)*(dz(3,:)*xz11+z(3)*dxz11)
                !(pz s/py s)
                dw(64,:)=dri(6,:)*xx32+ri(6)*dxx32+&
                     dri(7,:)*zz32+ri(7)*dzz32
                !(pz s/py px)
                dw(65,:)=dx(3,:)*(ri(13)*xx21+ri(14)*yyzz21)&
                     +x(3)*(dri(13,:)*xx21+ri(13)*dxx21+dri(14,:)*yyzz21+ri(14)*dyyzz21) &
                     +dri(15,:)*(z(3)*xz21)+ri(15)*(dz(3,:)*xz21+z(3)*dxz21)
                !(pz s/py py)
                dw(66,:)=dx(3,:)*(ri(13)*xx22+ri(14)*yyzz22)&
                     +x(3)*(dri(13,:)*xx22+ri(13)*dxx22+dri(14,:)*yyzz22+ri(14)*dyyzz22) &
                     +dri(15,:)*(z(3)*xz22)+ri(15)*(dz(3,:)*xz22+z(3)*dxz22)
                !(pz s/pz s)
                dw(67,:)=dri(6,:)*xx33+ri(6)*dxx33+&
                     dri(7,:)*zz33+ri(7)*dzz33
                !(pz s/pz px)
                dw(68,:)=dx(3,:)*(ri(13)*xx31+ri(14)*zz31)&
                     +x(3)*(dri(13,:)*xx31+ri(13)*dxx31+dri(14,:)*zz31+ri(14)*dzz31) &
                     +dri(15,:)*(z(3)*xz31)+ri(15)*(dz(3,:)*xz31+z(3)*dxz31)
                !(pz s/pz py)
                dw(69,:)=dx(3,:)*(ri(13)*xx32+ri(14)*zz32)&
                     +x(3)*(dri(13,:)*xx32+ri(13)*dxx32+dri(14,:)*zz32+ri(14)*dzz32) &
                     +dri(15,:)*(z(3)*xz32)+ri(15)*(dz(3,:)*xz32+z(3)*dxz32)
                !(pz s/pz pz)
                dw(70,:)=dx(3,:)*(ri(13)*xx33+ri(14)*zz33)&
                     +x(3)*(dri(13,:)*xx33+ri(13)*dxx33+dri(14,:)*zz33+ri(14)*dzz33) &
                     +dri(15,:)*(z(3)*xz33)+ri(15)*(dz(3,:)*xz33+z(3)*dxz33)
                !(pz px/s s)
                dw(71,:)=dri(3,:)*xx31+ri(3)*dxx31+&
                     dri(4,:)*zz31+ri(4)*dzz31
                !(pz px/px s)
                dw(72,:)=dx(1,:)*(ri(8)*xx31+ri(9)*zz31)&
                     +x(1)*(dri(8,:)*xx31+ri(8)*dxx31+dri(9,:)*zz31+ri(9)*dzz31) &
                     +dri(10,:)*(y(1)*xy31+z(1)*xz31)+ri(10)*(dy(1,:)*xy31+y(1)*dxy31+dz(1,:)*xz31+z(1)*dxz31)
                !(pz px/px px)
                dw(73,:) = &
                     (dri(16,:)*xx31+ri(16)*dxx31+dri(17,:)*zz31+ri(17)*dzz31)*xx11&
                     +(ri(16)*xx31+ri(17)*zz31)*dxx11&
                     +dri(18,:)*xx31*yyzz11+ri(18)*dxx31*yyzz11+ri(18)*xx31*dyyzz11 &
                     +dri(19,:)*zz31*zz11+ri(19)*dzz31*zz11+ri(19)*zz31*dzz11 &
                     +dri(20,:)*(xy31*xy11+xz31*xz11)+ri(20)*(dxy31*xy11+xy31*dxy11+dxz31*xz11+xz31*dxz11) &
                     +dri(21,:)*zz31*yy11+ri(21)*dzz31*yy11+ri(21)*zz31*dyy11 &
                     +dri(22,:)*yz31*yz11+ri(22)*dyz31*yz11+ri(22)*yz31*dyz11
                !(pz px/py s)
                dw(74,:)=dx(2,:)*(ri(8)*xx31+ri(9)*zz31)+x(2)*(dri(8,:)*xx31+ri(8)*dxx31+dri(9,:)*zz31+ri(9)*dzz31) &
                     +dri(10,:)*(y(2)*xy31+z(2)*xz31)+ri(10)*(dy(2,:)*xy31+y(2)*dxy31+dz(2,:)*xz31+z(2)*dxz31)
                !(pz px/py px)
                dw(75,:) = &
                     (dri(16,:)*xx31+ri(16)*dxx31+dri(17,:)*zz31+ri(17)*dzz31)*xx21&
                     +(ri(16)*xx31+ri(17)*zz31)*dxx21&
                     +dri(18,:)*xx31*yyzz21+ri(18)*dxx31*yyzz21+ri(18)*xx31*dyyzz21 &
                     +dri(19,:)*zz31*zz21+ri(19)*dzz31*zz21+ri(19)*zz31*dzz21 &
                     +dri(20,:)*(xy31*xy21+xz31*xz21)+ri(20)*(dxy31*xy21+xy31*dxy21+dxz31*xz21+xz31*dxz21) &
                     +dri(21,:)*zz31*yy21+ri(21)*dzz31*yy21+ri(21)*zz31*dyy21 &
                     +dri(22,:)*yz31*yz21+ri(22)*dyz31*yz21+ri(22)*yz31*dyz21
                !(pz px/py py)
                dw(76,:) = &
                     (dri(16,:)*xx31+ri(16)*dxx31+dri(17,:)*zz31+ri(17)*dzz31)*xx22&
                     +(ri(16)*xx31+ri(17)*zz31)*dxx22&
                     +dri(18,:)*xx31*yyzz22+ri(18)*dxx31*yyzz22+ri(18)*xx31*dyyzz22 &
                     +dri(19,:)*zz31*zz22+ri(19)*dzz31*zz22+ri(19)*zz31*dzz22 &
                     +dri(20,:)*(xy31*xy22+xz31*xz22)+ri(20)*(dxy31*xy22+xy31*dxy22+dxz31*xz22+xz31*dxz22) &
                     +dri(21,:)*zz31*yy22+ri(21)*dzz31*yy22+ri(21)*zz31*dyy22 &
                     +dri(22,:)*yz31*yz22+ri(22)*dyz31*yz22+ri(22)*yz31*dyz22
                !(pz px/pz s)
                dw(77,:)=dx(3,:)*(ri(8)*xx31+ri(9)*zz31)+x(3)*(dri(8,:)*xx31+ri(8)*dxx31+dri(9,:)*zz31+ri(9)*dzz31) &
                     +dri(10,:)*(z(3)*xz31)+ri(10)*(dz(3,:)*xz31+z(3)*dxz31)
                !(pz px/pz px)
                dw(78,:) = &
                     (dri(16,:)*xx31+ri(16)*dxx31+dri(17,:)*zz31+ri(17)*dzz31)*xx31&
                     +(ri(16)*xx31+ri(17)*zz31)*dxx31 &
                     +(dri(18,:)*xx31+ri(18)*dxx31+dri(19,:)*zz31+ri(19)*dzz31)*zz31&
                     +(ri(18)*xx31+ri(19)*zz31)*dzz31 &
                     +dri(20,:)*(xy31*xy31+xz31*xz31)+ri(20)*(dxy31*xy31+xy31*dxy31+dxz31*xz31+xz31*dxz31) &
                     +dri(22,:)*yz31*yz31+ri(22)*dyz31*yz31+ri(22)*yz31*dyz31
                !(pz px/pz py)
                dw(79,:) = &
                     (dri(16,:)*xx31+ri(16)*dxx31+dri(17,:)*zz31+ri(17)*dzz31)*xx32&
                     +(ri(16)*xx31+ri(17)*zz31)*dxx32 &
                     +(dri(18,:)*xx31+ri(18)*dxx31+dri(19,:)*zz31+ri(19)*dzz31)*zz32&
                     +(ri(18)*xx31+ri(19)*zz31)*dzz32 &
                     +dri(20,:)*(xy31*xy32+xz31*xz32)+ri(20)*(dxy31*xy32+xy31*dxy32+dxz31*xz32+xz31*dxz32) &
                     +dri(22,:)*yz31*yz32+ri(22)*dyz31*yz32+ri(22)*yz31*dyz32
                !(pz px/pz pz)
                dw(80,:) = &
                     (dri(16,:)*xx31+ri(16)*dxx31+dri(17,:)*zz31+ri(17)*dzz31)*xx33&
                     +(ri(16)*xx31+ri(17)*zz31)*dxx33 &
                     +(dri(18,:)*xx31+ri(18)*dxx31+dri(19,:)*zz31+ri(19)*dzz31)*zz33&
                     +(ri(18)*xx31+ri(19)*zz31)*dzz33 &
                     +dri(20,:)*xz31*xz33+ri(20)*dxz31*xz33+ri(20)*xz31*dxz33
                !(pz py/s s)
                dw(81,:)=dri(3,:)*xx32+ri(3)*dxx32+&
                     dri(4,:)*zz32+ri(4)*dzz32
                !(pz py/px s)
                dw(82,:)=dx(1,:)*(ri(8)*xx32+ri(9)*zz32)&
                     +x(1)*(dri(8,:)*xx32+ri(8)*dxx32+dri(9,:)*zz32+ri(9)*dzz32) &
                     +dri(10,:)*(y(1)*xy32+z(1)*xz32)+ri(10)*(dy(1,:)*xy32+y(1)*dxy32+dz(1,:)*xz32+z(1)*dxz32)
                !(pz py/px px)
                dw(83,:) = &
                     (dri(16,:)*xx32+ri(16)*dxx32+dri(17,:)*zz32+ri(17)*dzz32)*xx11&
                     +(ri(16)*xx32+ri(17)*zz32)*dxx11&
                     +dri(18,:)*xx32*yyzz11+ri(18)*dxx32*yyzz11 +ri(18)*xx32*dyyzz11 &
                     +dri(19,:)*zz32*zz11+ri(19)*dzz32*zz11+ri(19)*zz32*dzz11 &
                     +dri(20,:)*(xy32*xy11+xz32*xz11)+ri(20)*(dxy32*xy11+xy32*dxy11+dxz32*xz11+xz32*dxz11) &
                     +dri(21,:)*zz32*yy11+ri(21)*dzz32*yy11+ri(21)*zz32*dyy11 &
                     +dri(22,:)*yz32*yz11+ri(22)*dyz32*yz11+ri(22)*yz32*dyz11
                !(pz py/py s)
                dw(84,:)=dx(2,:)*(ri(8)*xx32+ri(9)*zz32)&
                     +x(2)*(dri(8,:)*xx32+ri(8)*dxx32+dri(9,:)*zz32+ri(9)*dzz32)&
                     +dri(10,:)*(y(2)*xy32+z(2)*xz32)&
                     +ri(10)*(dy(2,:)*xy32+y(2)*dxy32+dz(2,:)*xz32+z(2)*dxz32)
                !(pz py/py px)
                dw(85,:) = &
                     (dri(16,:)*xx32+ri(16)*dxx32+dri(17,:)*zz32+ri(17)*dzz32)*xx21&
                     +(ri(16)*xx32+ri(17)*zz32)*dxx21&
                     +dri(18,:)*xx32*yyzz21+ri(18)*dxx32*yyzz21+ri(18)*xx32*dyyzz21 &
                     +dri(19,:)*zz32*zz21+ri(19)*dzz32*zz21+ri(19)*zz32*dzz21 &
                     +dri(20,:)*(xy32*xy21+xz32*xz21)+ri(20)*(dxy32*xy21+xy32*dxy21+dxz32*xz21+xz32*dxz21) &
                     +dri(21,:)*zz32*yy21+ri(21)*dzz32*yy21+ri(21)*zz32*dyy21 &
                     +dri(22,:)*yz32*yz21+ri(22)*dyz32*yz21+ri(22)*yz32*dyz21
                !(pz py/py py)
                dw(86,:) = &
                     (dri(16,:)*xx32+ri(16)*dxx32+dri(17,:)*zz32+ri(17)*dzz32)*xx22&
                     +(ri(16)*xx32+ri(17)*zz32)*dxx22&
                     +dri(18,:)*xx32*yyzz22+ri(18)*dxx32*yyzz22+ri(18)*xx32*dyyzz22 &
                     +dri(19,:)*zz32*zz22+ri(19)*dzz32*zz22+ri(19)*zz32*dzz22 &
                     +dri(20,:)*(xy32*xy22+xz32*xz22)+ri(20)*(dxy32*xy22+xy32*dxy22+dxz32*xz22+xz32*dxz22) &
                     +dri(21,:)*zz32*yy22+ri(21)*dzz32*yy22+ri(21)*zz32*dyy22 &
                     +dri(22,:)*yz32*yz22+ri(22)*dyz32*yz22+ri(22)*yz32*dyz22
                !(pz py/pz s)
                dw(87,:)=dx(3,:)*(ri(8)*xx32+ri(9)*zz32)+x(3)*(dri(8,:)*xx32+ri(8)*dxx32+dri(9,:)*zz32+ri(9)*dzz32) &
                     +dri(10,:)*(z(3)*xz32)+ri(10)*(dz(3,:)*xz32+z(3)*dxz32)
                !(pz py/pz px)
                dw(88,:) = &
                     (dri(16,:)*xx32+ri(16)*dxx32+dri(17,:)*zz32+ri(17)*dzz32)*xx31&
                     +(ri(16)*xx32+ri(17)*zz32)*dxx31 &
                     +(dri(18,:)*xx32+ri(18)*dxx32+dri(19,:)*zz32+ri(19)*dzz32)*zz31&
                     +(ri(18)*xx32+ri(19)*zz32)*dzz31 &
                     +dri(20,:)*(xy32*xy31+xz32*xz31)+ri(20)*(dxy32*xy31+xy32*dxy31+dxz32*xz31+xz32*dxz31) &
                     +dri(22,:)*yz32*yz31+ri(22)*dyz32*yz31+ri(22)*yz32*dyz31
                !(pz py/pz py)
                dw(89,:) = &
                     (dri(16,:)*xx32+ri(16)*dxx32+dri(17,:)*zz32+ri(17)*dzz32)*xx32&
                     +(ri(16)*xx32+ri(17)*zz32)*dxx32 &
                     +(dri(18,:)*xx32+ri(18)*dxx32+dri(19,:)*zz32+ri(19)*dzz32)*zz32&
                     +(ri(18)*xx32+ri(19)*zz32)*dzz32 &
                     +dri(20,:)*(xy32*xy32+xz32*xz32)+ri(20)*(dxy32*xy32+xy32*dxy32+dxz32*xz32+xz32*dxz32) &
                     +dri(22,:)*yz32*yz32+ri(22)*dyz32*yz32+ri(22)*yz32*dyz32
                !(pz py/pz pz)
                dw(90,:) = &
                     (dri(16,:)*xx32+ri(16)*dxx32+dri(17,:)*zz32+ri(17)*dzz32)*xx33&
                     +(ri(16)*xx32+ri(17)*zz32)*dxx33 &
                     +(dri(18,:)*xx32+ri(18)*dxx32+dri(19,:)*zz32+ri(19)*dzz32)*zz33&
                     +(ri(18)*xx32+ri(19)*zz32)*dzz33 &
                     +dri(20,:)*xz32*xz33+ri(20)*dxz32*xz33+ri(20)*xz32*dxz33
                !(pz pz/s s)
                dw(91,:)=dri(3,:)*xx33+ri(3)*dxx33+&
                     dri(4,:)*zz33+ri(4)*dzz33
                !(pz pz/px s)
                dw(92,:)=dx(1,:)*(ri(8)*xx33+ri(9)*zz33)&
                     +x(1)*(dri(8,:)*xx33+ri(8)*dxx33+dri(9,:)*zz33+ri(9)*dzz33) &
                     +dri(10,:)*(z(1)*xz33)+ri(10)*(dz(1,:)*xz33+z(1)*dxz33)
                !(pz pz/px px)
                dw(93,:) = &
                     (dri(16,:)*xx33+ri(16)*dxx33+dri(17,:)*zz33+ri(17)*dzz33)*xx11&
                     +(ri(16)*xx33+ri(17)*zz33)*dxx11&
                     +dri(18,:)*xx33*yyzz11+ri(18)*dxx33*yyzz11+ri(18)*xx33*dyyzz11 &
                     +dri(19,:)*zz33*zz11+ri(19)*dzz33*zz11+ri(19)*zz33*dzz11 &
                     +dri(20,:)*xz33*xz11+ri(20)*dxz33*xz11+ri(20)*xz33*dxz11 &
                     +dri(21,:)*zz33*yy11+ri(21)*dzz33*yy11+ri(21)*zz33*dyy11
                !(pz pz/py s)
                dw(94,:)=dx(2,:)*(ri(8)*xx33+ri(9)*zz33)&
                     +x(2)*(dri(8,:)*xx33+ri(8)*dxx33+dri(9,:)*zz33+ri(9)*dzz33) &
                     +dri(10,:)*(z(2)*xz33)+ri(10)*(dz(2,:)*xz33+z(2)*dxz33)
                !(pz pz/py px)
                dw(95,:) = &
                     (dri(16,:)*xx33+ri(16)*dxx33+dri(17,:)*zz33+ri(17)*dzz33)*xx21&
                     +(ri(16)*xx33+ri(17)*zz33)*dxx21&
                     +dri(18,:)*xx33*yyzz21+ri(18)*dxx33*yyzz21+ri(18)*xx33*dyyzz21&
                     +dri(19,:)*zz33*zz21+ri(19)*dzz33*zz21+ri(19)*zz33*dzz21 &
                     +dri(20,:)*xz33*xz21+ri(20)*dxz33*xz21+ri(20)*xz33*dxz21 &
                     +dri(21,:)*zz33*yy21+ri(21)*dzz33*yy21+ri(21)*zz33*dyy21
                !(pz pz/py py)
                dw(96,:) = &
                     (dri(16,:)*xx33+ri(16)*dxx33+dri(17,:)*zz33+ri(17)*dzz33)*xx22&
                     +(ri(16)*xx33+ri(17)*zz33)*dxx22&
                     +dri(18,:)*xx33*yyzz22+ri(18)*dxx33*yyzz22+ri(18)*xx33*dyyzz22 &
                     +dri(19,:)*zz33*zz22  +ri(19)*dzz33*zz22  +ri(19)*zz33*dzz22   &
                     +dri(20,:)*xz33*xz22  +ri(20)*dxz33*xz22  +ri(20)*xz33*dxz22   &
                     +dri(21,:)*zz33*yy22  +ri(21)*dzz33*yy22  +ri(21)*zz33*dyy22   
                !(pz pz/pz s)
                dw(97,:)=dx(3,:)*(ri(8)*xx33+ri(9)*zz33)&
                     +x(3)*(dri(8,:)*xx33+ri(8)*dxx33+dri(9,:)*zz33+ri(9)*dzz33) &
                     +dri(10,:)*(z(3)*xz33)+ri(10)*(dz(3,:)*xz33+z(3)*dxz33)
                !(pz pz/pz px)
                dw(98,:) = &
                     (dri(16,:)*xx33+ri(16)*dxx33+dri(17,:)*zz33+ri(17)*dzz33)*xx31&
                     +(ri(16)*xx33+ri(17)*zz33)*dxx31 &
                     +(dri(18,:)*xx33+ri(18)*dxx33+dri(19,:)*zz33+ri(19)*dzz33)*zz31&
                     +(ri(18)*xx33+ri(19)*zz33)*dzz31 &
                     +dri(20,:)*xz33*xz31+ri(20)*dxz33*xz31+ri(20)*xz33*dxz31
                !(pz pz/pz py)
                dw(99,:) = &
                     (dri(16,:)*xx33+ri(16)*dxx33+dri(17,:)*zz33+ri(17)*dzz33)*xx32&
                     +(ri(16)*xx33+ri(17)*zz33)*dxx32 &
                     +(dri(18,:)*xx33+ri(18)*dxx33+dri(19,:)*zz33+ri(19)*dzz33)*zz32&
                     +(ri(18)*xx33+ri(19)*zz33)*dzz32 &
                     +dri(20,:)*xz33*xz32+ri(20)*dxz33*xz32+ri(20)*xz33*dxz32
                !(pz pz/pz pz)
                dw(100,:) = &
                     (dri(16,:)*xx33+ri(16)*dxx33+dri(17,:)*zz33+ri(17)*dzz33)*xx33&
                     +(ri(16)*xx33+ri(17)*zz33)*dxx33 &
                     +(dri(18,:)*xx33+ri(18)*dxx33+dri(19,:)*zz33+ri(19)*dzz33)*zz33&
                     +(ri(18)*xx33+ri(19)*zz33)*dzz33 &
                     +dri(20,:)*xz33*xz33+ri(20)*dxz33*xz33+ri(20)*xz33*dxz33
             ELSE
                !(px s/s s)
                dw(2,:)=dri(2,:)*x(1)+ri(2)*dx(1,:)
                !(px px/s s)
                dw(3,:)=dri(3,:)*xx11+ri(3)*dxx11+&
                     dri(4,:)*yyzz11+ ri(4)*dyyzz11
                !(py s/s s)
                dw(4,:)=dri(2,:)*x(2)+ri(2)*dx(2,:)
                !(py px/s s)
                dw(5,:)=dri(3,:)*xx21+ri(3)*dxx21+&
                     dri(4,:)*yyzz21+ri(4)*dyyzz21
                !(py py/s s)
                dw(6,:)=dri(3,:)*xx22+ri(3)*dxx22+&
                     dri(4,:)*yyzz22+ri(4)*dyyzz22
                !(pz s/ss)
                dw(7,:)=dri(2,:)*x(3)+ri(2)*dx(3,:)
                !(pz px/s s)
                dw(8,:)=dri(3,:)*xx31+ri(3)*dxx31+&
                     dri(4,:)*zz31+ri(4)*dzz31
                !(pz py/s s)
                dw(9,:)=dri(3,:)*xx32+ri(3)*dxx32+&
                     dri(4,:)*zz32+ri(4)*dzz32
                !(pz pz/s s)
                dw(10,:)=dri(3,:)*xx33+ri(3)*dxx33+&
                     dri(4,:)*zz33+ri(4)*dzz33
             END IF
          END IF
          IF (invert) CALL invert_derivative(dw,42)
          IF (debug_this_module) THEN
             ! Check derivatives
             ! Numerical derivatives are obviosly a big problem..
             ! First of all let's decide if the value we get for delta is compatible
             ! with a reasonable value of the integral.. (compatible if the value of the
             ! integral is greater than 1.0E-6)
             delta = 1.0E-5_dp
             CALL drotint(sepi,sepj,rijv,dw2,delta=delta,shortrange=shortrange)
             CALL rotint_ana(sepi,sepj,rijv,w2,shortrange=shortrange)
             DO i = 1, 3
                DO j = 1, 100
                   IF ((ABS(w2(j))>1.0E-6_dp).AND.(ABS(dw2(j,i))>delta*10)) THEN
                      IF (ABS((dw2(j,i)-dw(j,i))/dw(j,i))*100.0_dp>1.0_dp) THEN
                         WRITE(*,*)"check de1b",i,j,ABS((dw2(j,i)-dw(j,i))/dw(j,i))*100.0_dp
                         WRITE(*,'(10F12.6)')dw
                         WRITE(*,*)
                         WRITE(*,'(10F12.6)')dw2
                         WRITE(*,*)
                         WRITE(*,'(10F12.6)')dw2-dw
                         WRITE(*,*)invert
                         STOP
                      END IF
                   END IF
                END DO
             END DO
          END IF
       END IF
    ENDIF

  END SUBROUTINE rotint_ana

!!****f* semi_empirical_int_ana/dterep_ana [1.0] *
!!
!!   NAME
!!     dterep_ana
!!
!!   FUNCTION
!!     Calculates the derivative pf two-electron repulsion integrals and the
!!     nuclear attraction integrals w.r.t. |r|
!!
!!   NOTES
!!     Analytical version - Analytical evaluation of gradients
!!     Teodoro Laino - Zurich University 04.2007
!!     routine adapted from mopac7 (repp)
!!     vector version written by Ernest R. Davidson, Indiana University
!!
!!   INPUTS
!!     on input rij     = interatomic distance
!!              sepi    = paramters of atom i
!!              sepj    = paramters of atom j
!!    on output ri      = array of two-electron repulsion integrals
!!             dri      = array of derivatives of two-electron repulsion integrals
!!       The two-centre repulsion integrals (over local coordinates) are
!!       stored as follows (where p-sigma = O,  and p-pi = P and P* )
!!       (SS/SS)=1,   (SO/SS)=2,   (OO/SS)=3,   (PP/SS)=4,   (SS/OS)=5,
!!       (SO/SO)=6,   (SP/SP)=7,   (OO/SO)=8,   (PP/SO)=9,   (PO/SP)=10,
!!       (SS/OO)=11,  (SS/PP)=12,  (SO/OO)=13,  (SO/PP)=14,  (SP/OP)=15,
!!       (OO/OO)=16,  (PP/OO)=17,  (OO/PP)=18,  (PP/PP)=19,  (PO/PO)=20,
!!       (PP/P*P*)=21,   (P*P/P*P)=22.
!!
!!   AUTHOR
!!     Teodoro Laino - Zurich University
!!
!!   MODIFICATION HISTORY
!!     04.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE dterep_ana ( sepi, sepj, rij, ri, dri, shortrange)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: ri, dri
    LOGICAL, INTENT(IN)                      :: shortrange

    LOGICAL                                  :: si, sj
    REAL(dp) :: ade, adi, adj, adq, aed, aee, aeq, ami, amj, aqd, aqe, aqi, &
      aqj, aqq, axx, da, db, ddi, ddj, ddxdx, ddxqxz, ddzdz, ddze, ddzqxx, &
      ddzqzz, dedz, dee, deqxx, deqzz, dft, dqxxdz, dqxxe, dqxxqxx, dqxxqyy, &
      dqxxqzz, dqxzdx, dqxzqxz, dqzzdz, dqzze, dqzzqxx, dqzzqzz, drsq, dwww, &
      dxdx, dxqxz, dxxx, dyyy, dzdz, dze, dzqxx, dzqzz, dzzz, edz, ee, eqxx, &
      eqzz, fac, ft, qa, qb, qqi, qqj, qxxdz, qxxe, qxxqxx, qxxqyy, qxxqzz, &
      qxzdx, qxzqxz, qzzdz, qzze, qzzqxx, qzzqzz, r, rsq, www, xxx, yyy, zi, &
      zj, zzz
    REAL(dp), DIMENSION(22)                  :: ris
    REAL(dp), DIMENSION(72)                  :: arg, darg, dsqr, sqr

    ri  = 0._dp
    dri = 0._dp
    r=rij
    si = (sepi%natorb >= 3)
    sj = (sepj%natorb >= 3)
    zi = sepi%zeff
    zj = sepj%zeff
    ddi = sepi%dd
    ddj = sepj%dd
    qqi = sepi%qq
    qqj = sepj%qq
    IF ((.NOT.si) .AND. (.NOT.sj)) THEN
       !
       !  hydrogen - hydrogen  (SS/SS)
       !
       ami = sepi%am
       amj = sepj%am
       aee = 0.5_dp/ami + 0.5_dp/amj
       aee = aee * aee
       ri(1) = 1._dp/SQRT(r*r+aee)
       !
       fac   = - r/(r*r+aee)
       dri(1) = fac*ri(1)
    ELSE IF (si .AND. (.NOT.sj)) THEN
       !
       !  heavy atom - hydrogen
       !
       ami = sepi%am
       adi = sepi%ad
       aqi = sepi%aq
       amj = sepj%am
       aee = 0.5_dp/ami + 0.5_dp/amj
       aee = aee * aee
       da=ddi
       qa=qqi * 2._dp
       ade = 0.5_dp/adi + 0.5_dp/amj
       ade = ade * ade
       aqe = 0.5_dp/aqi + 0.5_dp/amj
       aqe = aqe * aqe
       rsq = r*r
       drsq= 2.0_dp*r
       arg(1) = rsq + aee
       darg(1)= drsq
       xxx = r+da
       arg(2) = xxx*xxx + ade
       darg(2)= 2.0_dp*xxx
       xxx = r-da
       arg(3) = xxx*xxx + ade
       darg(3)= 2.0_dp*xxx
       xxx = r+qa
       arg(4) = xxx*xxx + aqe
       darg(4)= 2.0_dp*xxx
       xxx = r-qa
       arg(5) = xxx*xxx + aqe
       darg(5)= 2.0_dp*xxx
       arg(6) = rsq + aqe
       darg(6)= drsq
       arg(7) = arg(6) + qa*qa
       darg(7)= darg(6)
       sqr(1:7) = SQRT(arg(1:7))
       dsqr(1:7) = -(0.5_dp*(1.0_dp/sqr(1:7))**3)*darg(1:7)
       ee = 1._dp/sqr(1)
       dee= dsqr(1)
       ri(1) = ee
       ri(2) = 0.5_dp/sqr(2) - 0.5_dp/sqr(3)
       ri(3) = ee + 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
       ri(4) = ee + 0.5_dp/sqr(7) - 0.5_dp/sqr(6)

       dri(1) = dee
       dri(2) = 0.5_dp*dsqr(2) - 0.5_dp*dsqr(3)
       dri(3) = dee + 0.25_dp*dsqr(4) + 0.25_dp*dsqr(5) - 0.5_dp*dsqr(6)
       dri(4) = dee + 0.5_dp*dsqr(7) - 0.5_dp*dsqr(6)
 
    ELSE IF ((.NOT.si).AND.sj) THEN
       !
       !  hydrogen - heavy atom
       !
       ami = sepi%am
       amj = sepj%am
       adj = sepj%ad
       aqj = sepj%aq
       aee = 0.5_dp/ami + 0.5_dp/amj
       aee = aee * aee
       db=ddj
       qb=qqj * 2._dp
       aed = 0.5_dp/ami + 0.5_dp/adj
       aed = aed * aed
       aeq = 0.5_dp/ami + 0.5_dp/aqj
       aeq = aeq * aeq
       rsq = r*r
       drsq= 2.0_dp*r
       arg(1) = rsq + aee
       darg(1)= drsq
       xxx = r-db
       arg(2) = xxx*xxx + aed
       darg(2)= 2.0_dp*xxx
       xxx = r+db
       arg(3) = xxx*xxx + aed
       darg(3)= 2.0_dp*xxx
       xxx = r-qb
       arg(4) = xxx*xxx + aeq
       darg(4)= 2.0_dp*xxx
       xxx = r+qb
       arg(5) = xxx*xxx + aeq
       darg(5)= 2.0_dp*xxx
       arg(6) = rsq + aeq
       darg(6)= drsq
       arg(7) = arg(6) + qb*qb
       darg(7)= darg(6)
       sqr(1:7) = SQRT(arg(1:7))
       dsqr(1:7) = -(0.5_dp*(1.0_dp/sqr(1:7))**3)*darg(1:7)
       ee = 1._dp/sqr(1)
       dee= dsqr(1)
       ri(1) = ee
       ri(5) = 0.5_dp/sqr(2)  - 0.5_dp/sqr(3)
       ri(11) = ee + 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
       ri(12) = ee + 0.5_dp/sqr(7) - 0.5_dp/sqr(6)

       dri(1) = dee
       dri(5) = 0.5_dp*dsqr(2)  - 0.5_dp*dsqr(3)
       dri(11) = dee + 0.25_dp*dsqr(4) + 0.25_dp*dsqr(5) - 0.5_dp*dsqr(6)
       dri(12) = dee + 0.5_dp*dsqr(7) - 0.5_dp*dsqr(6)
 
    ELSE
       !
       !  heavy atom - heavy atom
       !
       !  define charge separations.
       da=ddi
       db=ddj
       qa=qqi * 2._dp
       qb=qqj * 2._dp
 
       ami = sepi%am
       amj = sepj%am
       adi = sepi%ad
       adj = sepj%ad
       aqi = sepi%aq
       aqj = sepj%aq

       aee = 0.5_dp/ami + 0.5_dp/amj
       aee = aee * aee
       ade = 0.5_dp/adi + 0.5_dp/amj
       ade = ade * ade
       aqe = 0.5_dp/aqi + 0.5_dp/amj
       aqe = aqe * aqe
       aed = 0.5_dp/ami + 0.5_dp/adj
       aed = aed * aed
       aeq = 0.5_dp/ami + 0.5_dp/aqj
       aeq = aeq * aeq
       axx = 0.5_dp/adi + 0.5_dp/adj
       axx = axx * axx
       adq = 0.5_dp/adi + 0.5_dp/aqj
       adq = adq * adq
       aqd = 0.5_dp/aqi + 0.5_dp/adj
       aqd = aqd * aqd
       aqq = 0.5_dp/aqi + 0.5_dp/aqj
       aqq = aqq * aqq
       rsq = r * r
       drsq= 2.0_dp*r
       arg(1) = rsq + aee
       darg(1)= drsq
       xxx = r + da
       arg(2) = xxx * xxx + ade
       darg(2) = 2.0_dp* xxx
       xxx = r - da
       arg(3) = xxx*xxx + ade
       darg(3) = 2.0_dp* xxx
       xxx = r - qa
       arg(4) = xxx*xxx + aqe
       darg(4) = 2.0_dp* xxx
       xxx = r + qa
       arg(5) = xxx*xxx + aqe
       darg(5) = 2.0_dp* xxx
       arg(6) = rsq + aqe
       darg(6)= drsq
       arg(7) = arg(6) + qa*qa
       darg(7)= darg(6)
       xxx = r-db
       arg(8) = xxx*xxx + aed
       darg(8) = 2.0_dp* xxx
       xxx = r+db
       arg(9) = xxx*xxx + aed
       darg(9) = 2.0_dp* xxx
       xxx = r - qb
       arg(10) = xxx*xxx + aeq
       darg(10)= 2.0_dp* xxx
       xxx = r + qb
       arg(11) = xxx*xxx + aeq
       darg(11)= 2.0_dp* xxx
       arg(12) = rsq + aeq
       darg(12)= drsq
       arg(13) = arg(12) + qb*qb
       darg(13)= darg(12)
       xxx = da-db
       arg(14) = rsq + axx + xxx*xxx
       darg(14)= drsq
       xxx = da+db
       arg(15) = rsq + axx + xxx*xxx
       darg(15)= drsq
       xxx = r + da - db
       arg(16) = xxx*xxx + axx
       darg(16)= 2.0_dp* xxx
       xxx = r - da + db
       arg(17) = xxx*xxx + axx
       darg(17)= 2.0_dp* xxx
       xxx = r - da - db
       arg(18) = xxx*xxx + axx
       darg(18)= 2.0_dp* xxx
       xxx = r + da + db
       arg(19) = xxx*xxx + axx
       darg(19)= 2.0_dp* xxx
       xxx = r + da
       arg(20) = xxx*xxx + adq
       darg(20)= 2.0_dp* xxx
       arg(21) = arg(20) + qb*qb
       darg(21)= darg(20)
       xxx = r - da
       arg(22) = xxx*xxx + adq
       darg(22)= 2.0_dp* xxx
       arg(23) = arg(22) + qb*qb
       darg(23)= darg(22)
       xxx = r - db
       arg(24) = xxx*xxx + aqd
       darg(24)= 2.0_dp* xxx
       arg(25) = arg(24) + qa*qa
       darg(25)= darg(24)
       xxx = r + db
       arg(26) = xxx*xxx + aqd
       darg(26)= 2.0_dp* xxx
       arg(27) = arg(26) + qa*qa
       darg(27)= darg(26)
       xxx = r + da - qb
       arg(28) = xxx*xxx + adq
       darg(28)= 2.0_dp*xxx
       xxx = r - da - qb
       arg(29) = xxx*xxx + adq
       darg(29)= 2.0_dp* xxx
       xxx = r + da + qb
       arg(30) = xxx*xxx + adq
       darg(30)= 2.0_dp* xxx
       xxx = r - da + qb
       arg(31) = xxx*xxx + adq
       darg(31)= 2.0_dp* xxx
       xxx = r + qa - db
       arg(32) = xxx*xxx + aqd
       darg(32)= 2.0_dp* xxx
       xxx = r + qa + db
       arg(33) = xxx*xxx + aqd
       darg(33)= 2.0_dp* xxx
       xxx = r - qa - db
       arg(34) = xxx*xxx + aqd
       darg(34)= 2.0_dp* xxx
       xxx = r - qa + db
       arg(35) = xxx*xxx + aqd
       darg(35)= 2.0_dp* xxx
       arg(36) = rsq + aqq
       darg(36)= drsq
       xxx = qa - qb
       arg(37) = arg(36) + xxx*xxx
       darg(37)= darg(36)
       xxx = qa + qb
       arg(38) = arg(36) + xxx*xxx
       darg(38)= darg(36)
       arg(39) = arg(36) + qa*qa
       darg(39)= darg(36)
       arg(40) = arg(36) + qb*qb
       darg(40)= darg(36)
       arg(41) = arg(39) + qb*qb
       darg(41)= darg(39)
       xxx = r - qb
       arg(42) = xxx*xxx + aqq
       darg(42)= 2.0_dp*xxx
       arg(43) = arg(42) + qa*qa
       darg(43)= darg(42)
       xxx = r + qb
       arg(44) = xxx*xxx + aqq
       darg(44)= 2.0_dp*xxx
       arg(45) = arg(44) + qa*qa
       darg(45)= darg(44)
       xxx = r + qa
       arg(46) = xxx*xxx + aqq
       darg(46)= 2.0_dp*xxx
       arg(47) = arg(46) + qb*qb
       darg(47)= darg(46)
       xxx = r - qa
       arg(48) = xxx*xxx + aqq
       darg(48)= 2.0_dp*xxx
       arg(49) = arg(48) + qb*qb
       darg(49)= darg(48)
       xxx = r + qa - qb
       arg(50) = xxx*xxx + aqq
       darg(50)= 2.0_dp*xxx
       xxx = r + qa + qb
       arg(51) = xxx*xxx + aqq
       darg(51)= 2.0_dp*xxx
       xxx = r - qa - qb
       arg(52) = xxx*xxx + aqq
       darg(52)= 2.0_dp*xxx
       xxx = r - qa + qb
       arg(53) = xxx*xxx + aqq
       darg(53)= 2.0_dp*xxx
       qa=sepi%qq
       qb=sepj%qq
       xxx = da - qb
       dxxx= 0.0_dp
       xxx = xxx*xxx
       yyy = r - qb
       dyyy= 2.0_dp*yyy
       yyy = yyy*yyy
       zzz = da + qb
       dzzz= 0.0_dp
       zzz = zzz*zzz
       www = r + qb
       dwww= 2.0_dp*www
       www = www*www
       arg(54) =  xxx +  yyy + adq
       darg(54)= dxxx + dyyy 
       arg(55) =  xxx +  www + adq
       darg(55)= dxxx + dwww
       arg(56) =  zzz +  yyy + adq
       darg(56)= dzzz + dyyy
       arg(57) =  zzz +  www + adq
       darg(57)= dzzz + dwww
       xxx = qa - db
       dxxx= 0.0_dp
       xxx = xxx*xxx
       yyy = qa + db
       dyyy= 0.0_dp
       yyy = yyy*yyy
       zzz = r + qa
       dzzz= 2.0_dp*zzz
       zzz = zzz*zzz
       www = r - qa
       dwww= 2.0_dp*www
       www = www*www
       arg(58) =  zzz +  xxx + aqd
       darg(58)= dzzz + dxxx
       arg(59) =  www +  xxx + aqd
       darg(59)= dwww + dxxx
       arg(60) =  zzz +  yyy + aqd
       darg(60)= dzzz + dyyy
       arg(61) =  www +  yyy + aqd
       darg(61)= dwww + dyyy
       xxx = qa - qb
       xxx = xxx*xxx
       arg(62) =  arg(36) + 2._dp*xxx
       darg(62)= darg(36)
       yyy = qa + qb
       yyy = yyy*yyy
       arg(63) =  arg(36) + 2._dp*yyy
       darg(63)= darg(36)
       arg(64) =  arg(36) + 2._dp*(qa*qa+qb*qb)
       darg(64)= darg(36)
       zzz = r + qa - qb
       dzzz= 2.0_dp*zzz
       zzz = zzz*zzz
       arg(65) = zzz + xxx + aqq
       darg(65)= dzzz
       arg(66) = zzz + yyy + aqq
       darg(66)= dzzz
       zzz = r + qa + qb
       dzzz= 2.0_dp*zzz
       zzz = zzz*zzz
       arg(67) = zzz + xxx + aqq
       darg(67)= dzzz
       arg(68) = zzz + yyy + aqq
       darg(68)= dzzz
       zzz = r - qa - qb
       dzzz= 2.0_dp*zzz
       zzz = zzz*zzz
       arg(69) = zzz + xxx + aqq
       darg(69)= dzzz
       arg(70) = zzz + yyy + aqq
       darg(70)= dzzz
       zzz = r - qa + qb
       dzzz= 2.0_dp*zzz
       zzz = zzz*zzz
       arg(71) = zzz + xxx + aqq
       darg(71)= dzzz
       arg(72) = zzz + yyy + aqq
       darg(72)= dzzz
       sqr(1:72) = SQRT(arg(1:72))
       dsqr(1:72)= -(0.5_dp*(1/sqr(1:72))**3)*darg(1:72)

       ee     =  1._dp/sqr(1)
       dze    = -0.5_dp/sqr(2) + 0.5_dp/sqr(3)
       qzze   =  0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
       qxxe   =  0.5_dp/sqr(7) - 0.5_dp/sqr(6)
       edz    = -0.5_dp/sqr(8) + 0.5_dp/sqr(9)
       eqzz   =  0.25_dp/sqr(10) + 0.25_dp/sqr(11) - 0.5_dp/sqr(12)
       eqxx   =  0.5_dp/sqr(13) - 0.5_dp/sqr(12)
       dxdx   =  0.5_dp/sqr(14) - 0.5_dp/sqr(15)
       dzdz   =  0.25_dp/sqr(16) + 0.25_dp/sqr(17) - 0.25_dp/sqr(18) - 0.25_dp/sqr(19)
       dzqxx  =  0.25_dp/sqr(20) - 0.25_dp/sqr(21) - 0.25_dp/sqr(22) + 0.25_dp/sqr(23)
       qxxdz  =  0.25_dp/sqr(24) - 0.25_dp/sqr(25) - 0.25_dp/sqr(26) + 0.25_dp/sqr(27)
       dzqzz  = -0.125_dp/sqr(28) + 0.125_dp/sqr(29) - 0.125_dp/sqr(30) + 0.125_dp/sqr(31) &
               - 0.25_dp/sqr(22) + 0.25_dp/sqr(20)
       qzzdz  = -0.125_dp/sqr(32) + 0.125_dp/sqr(33) - 0.125_dp/sqr(34) + 0.125_dp/sqr(35) &
               + 0.25_dp/sqr(24) - 0.25_dp/sqr(26)
       qxxqxx =  0.125_dp/sqr(37) + 0.125_dp/sqr(38) - 0.25_dp/sqr(39) - 0.25_dp/sqr(40) &
               + 0.25_dp/sqr(36)
       qxxqyy =  0.25_dp/sqr(41) - 0.25_dp/sqr(39) - 0.25_dp/sqr(40) + 0.25_dp/sqr(36)
       qxxqzz =  0.125_dp/sqr(43) + 0.125_dp/sqr(45) - 0.125_dp/sqr(42) - 0.125_dp/sqr(44) &
               - 0.25_dp/sqr(39) + 0.25_dp/sqr(36)
       qzzqxx =  0.125_dp/sqr(47) + 0.125_dp/sqr(49) - 0.125_dp/sqr(46) - 0.125_dp/sqr(48) &
               - 0.25_dp/sqr(40) + 0.25_dp/sqr(36)
       qzzqzz =  0.0625_dp/sqr(50) + 0.0625_dp/sqr(51) + 0.0625_dp/sqr(52) + 0.0625_dp/sqr(53) &
               - 0.125_dp/sqr(48) - 0.125_dp/sqr(46) - 0.125_dp/sqr(42) - 0.125_dp/sqr(44) + &
                 0.25_dp/sqr(36)
       dxqxz  = -0.25_dp/sqr(54) + 0.25_dp/sqr(55) + 0.25_dp/sqr(56) - 0.25_dp/sqr(57)
       qxzdx  = -0.25_dp/sqr(58) + 0.25_dp/sqr(59) + 0.25_dp/sqr(60) - 0.25_dp/sqr(61)
       qxzqxz =  0.125_dp/sqr(65) - 0.125_dp/sqr(67) - 0.125_dp/sqr(69) + 0.125_dp/sqr(71) &
               - 0.125_dp/sqr(66) + 0.125_dp/sqr(68) + 0.125_dp/sqr(70) - 0.125_dp/sqr(72)

       dee     =  1._dp*dsqr(1)
       ddze    = -0.5_dp*dsqr(2) + 0.5_dp*dsqr(3)
       dqzze   =  0.25_dp*dsqr(4) + 0.25_dp*dsqr(5) - 0.5_dp*dsqr(6)
       dqxxe   =  0.5_dp*dsqr(7) - 0.5_dp*dsqr(6)
       dedz    = -0.5_dp*dsqr(8) + 0.5_dp*dsqr(9)
       deqzz   =  0.25_dp*dsqr(10) + 0.25_dp*dsqr(11) - 0.5_dp*dsqr(12)
       deqxx   =  0.5_dp*dsqr(13) - 0.5_dp*dsqr(12)
       ddxdx   =  0.5_dp*dsqr(14) - 0.5_dp*dsqr(15)
       ddzdz   =  0.25_dp*dsqr(16) + 0.25_dp*dsqr(17) - 0.25_dp*dsqr(18) - 0.25_dp*dsqr(19)
       ddzqxx  =  0.25_dp*dsqr(20) - 0.25_dp*dsqr(21) - 0.25_dp*dsqr(22) + 0.25_dp*dsqr(23)
       dqxxdz  =  0.25_dp*dsqr(24) - 0.25_dp*dsqr(25) - 0.25_dp*dsqr(26) + 0.25_dp*dsqr(27)
       ddzqzz  = -0.125_dp*dsqr(28) + 0.125_dp*dsqr(29) - 0.125_dp*dsqr(30) + 0.125_dp*dsqr(31) &
                - 0.25_dp*dsqr(22) + 0.25_dp*dsqr(20)
       dqzzdz  = -0.125_dp*dsqr(32) + 0.125_dp*dsqr(33) - 0.125_dp*dsqr(34) + 0.125_dp*dsqr(35) &
                + 0.25_dp*dsqr(24) - 0.25_dp*dsqr(26)
       dqxxqxx =  0.125_dp*dsqr(37) + 0.125_dp*dsqr(38) - 0.25_dp*dsqr(39) - 0.25_dp*dsqr(40) &
                + 0.25_dp*dsqr(36)
       dqxxqyy =  0.25_dp*dsqr(41) - 0.25_dp*dsqr(39) - 0.25_dp*dsqr(40) + 0.25_dp*dsqr(36)
       dqxxqzz =  0.125_dp*dsqr(43) + 0.125_dp*dsqr(45) - 0.125_dp*dsqr(42) - 0.125_dp*dsqr(44) &
                - 0.25_dp*dsqr(39) + 0.25_dp*dsqr(36)
       dqzzqxx =  0.125_dp*dsqr(47) + 0.125_dp*dsqr(49) - 0.125_dp*dsqr(46) - 0.125_dp*dsqr(48) &
                - 0.25_dp*dsqr(40) + 0.25_dp*dsqr(36)
       dqzzqzz =  0.0625_dp*dsqr(50) + 0.0625_dp*dsqr(51) + 0.0625_dp*dsqr(52) + 0.0625_dp*dsqr(53) &
                - 0.125_dp*dsqr(48) - 0.125_dp*dsqr(46) - 0.125_dp*dsqr(42) - 0.125_dp*dsqr(44) + &
                  0.25_dp*dsqr(36)
       ddxqxz  = -0.25_dp*dsqr(54) + 0.25_dp*dsqr(55) + 0.25_dp*dsqr(56) - 0.25_dp*dsqr(57)
       dqxzdx  = -0.25_dp*dsqr(58) + 0.25_dp*dsqr(59) + 0.25_dp*dsqr(60) - 0.25_dp*dsqr(61)
       dqxzqxz =  0.125_dp*dsqr(65) - 0.125_dp*dsqr(67) - 0.125_dp*dsqr(69) + 0.125_dp*dsqr(71) &
               - 0.125_dp*dsqr(66) + 0.125_dp*dsqr(68) + 0.125_dp*dsqr(70) - 0.125_dp*dsqr(72)

       ri(1) = ee
       ri(2) = -dze
       ri(3) = ee + qzze
       ri(4) = ee + qxxe
       ri(5) = -edz
       ri(6) = dzdz
       ri(7) = dxdx
       ri(8) = -edz -qzzdz
       ri(9) = -edz -qxxdz
       ri(10) = -qxzdx
       ri(11) =  ee + eqzz
       ri(12) =  ee + eqxx
       ri(13) = -dze -dzqzz
       ri(14) = -dze -dzqxx
       ri(15) = -dxqxz
       ri(16) = ee +eqzz +qzze +qzzqzz
       ri(17) = ee +eqzz +qxxe +qxxqzz
       ri(18) = ee +eqxx +qzze +qzzqxx
       ri(19) = ee +eqxx +qxxe +qxxqxx
       ri(20) = qxzqxz
       ri(21) = ee +eqxx +qxxe +qxxqyy
       ri(22) = 0.5_dp * (qxxqxx -qxxqyy)

       dri(1) = dee
       dri(2) = -ddze
       dri(3) = dee + dqzze
       dri(4) = dee + dqxxe
       dri(5) = -dedz
       dri(6) = ddzdz
       dri(7) = ddxdx
       dri(8) = -dedz -dqzzdz
       dri(9) = -dedz -dqxxdz
       dri(10) = -dqxzdx
       dri(11) =  dee + deqzz
       dri(12) =  dee + deqxx
       dri(13) = -ddze -ddzqzz
       dri(14) = -ddze -ddzqxx
       dri(15) = -ddxqxz
       dri(16) = dee +deqzz +dqzze +dqzzqzz
       dri(17) = dee +deqzz +dqxxe +dqxxqzz
       dri(18) = dee +deqxx +dqzze +dqzzqxx
       dri(19) = dee +deqxx +dqxxe +dqxxqxx
       dri(20) = dqxzqxz
       dri(21) = dee +deqxx +dqxxe +dqxxqyy
       dri(22) = 0.5_dp * (dqxxqxx -dqxxqyy)
 
    END IF

    IF ( shortrange ) THEN
      IF ((.NOT.si) .AND. (.NOT.sj)) THEN
         !
         !  hydrogen - hydrogen  (SS/SS)
         !
         ris(1) = 1._dp/r
         ri(1) = ri(1) - ris(1)
         !
         fac   = -1._dp/r
         dri(1) = dri(1) - fac*ris(1)
      ELSE IF (si .AND. (.NOT.sj)) THEN
         !
         !  heavy atom - hydrogen
         !
         da=ddi
         qa=qqi * 2._dp
         rsq = r*r
         drsq= 2.0_dp*r
         arg(1) = rsq
         darg(1)= drsq
         xxx = r+da
         arg(2) = xxx*xxx
         darg(2)= 2.0_dp*xxx
         xxx = r-da
         arg(3) = xxx*xxx
         darg(3)= 2.0_dp*xxx
         xxx = r+qa
         arg(4) = xxx*xxx
         darg(4)= 2.0_dp*xxx
         xxx = r-qa
         arg(5) = xxx*xxx
         darg(5)= 2.0_dp*xxx
         arg(6) = rsq
         darg(6)= drsq
         arg(7) = arg(6) + qa*qa
         darg(7)= darg(6)
         sqr(1:7) = SQRT(arg(1:7))
         dsqr(1:7) = -(0.5_dp*(1.0_dp/sqr(1:7))**3)*darg(1:7)
         ee = 1._dp/sqr(1)
         dee= dsqr(1)
         ris(1) = ee
         ris(2) = 0.5_dp/sqr(2) - 0.5_dp/sqr(3)
         ris(3) = ee + 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
         ris(4) = ee + 0.5_dp/sqr(7) - 0.5_dp/sqr(6)
         ri(1) = ri(1) - ris(1)
         ri(2) = ri(2) - ris(2)
         ri(3) = ri(3) - ris(3)
         ri(4) = ri(4) - ris(4)
  
         dri(1) = dri(1) - (dee)
         dri(2) = dri(2) - (0.5_dp*dsqr(2) - 0.5_dp*dsqr(3))
         dri(3) = dri(3) - (dee + 0.25_dp*dsqr(4) + 0.25_dp*dsqr(5) - 0.5_dp*dsqr(6))
         dri(4) = dri(4) - (dee + 0.5_dp*dsqr(7) - 0.5_dp*dsqr(6))
   
      ELSE IF ((.NOT.si).AND.sj) THEN
         !
         !  hydrogen - heavy atom
         !
         db=ddj
         qb=qqj * 2._dp
         rsq = r*r
         drsq= 2.0_dp*r
         arg(1) = rsq
         darg(1)= drsq
         xxx = r-db
         arg(2) = xxx*xxx
         darg(2)= 2.0_dp*xxx
         xxx = r+db
         arg(3) = xxx*xxx
         darg(3)= 2.0_dp*xxx
         xxx = r-qb
         arg(4) = xxx*xxx
         darg(4)= 2.0_dp*xxx
         xxx = r+qb
         arg(5) = xxx*xxx
         darg(5)= 2.0_dp*xxx
         arg(6) = rsq
         darg(6)= drsq
         arg(7) = arg(6) + qb*qb
         darg(7)= darg(6)
         sqr(1:7) = SQRT(arg(1:7))
         dsqr(1:7) = -(0.5_dp*(1.0_dp/sqr(1:7))**3)*darg(1:7)
         ee = 1._dp/sqr(1)
         dee= dsqr(1)
         ris(1) = ee
         ris(5) = 0.5_dp/sqr(2)  - 0.5_dp/sqr(3)
         ris(11) = ee + 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
         ris(12) = ee + 0.5_dp/sqr(7) - 0.5_dp/sqr(6)
         ri(1) = ri(1) - ris(1)
         ri(5) = ri(5) - ris(5)
         ri(11) = ri(11) - ris(11)
         ri(12) = ri(12) - ris(12)
  
         dri(1) = dri(1) - (dee)
         dri(5) = dri(5) - (0.5_dp*dsqr(2)  - 0.5_dp*dsqr(3))
         dri(11) = dri(11) - (dee + 0.25_dp*dsqr(4) + 0.25_dp*dsqr(5) - 0.5_dp*dsqr(6))
         dri(12) = dri(12) - (dee + 0.5_dp*dsqr(7) - 0.5_dp*dsqr(6))
   
      ELSE
         !
         !  heavy atom - heavy atom
         !
         !  define charge separations.
         da=ddi
         db=ddj
         qa=qqi * 2._dp
         qb=qqj * 2._dp
   
         rsq = r * r
         drsq= 2.0_dp*r
         arg(1) = rsq
         darg(1)= drsq
         xxx = r + da
         arg(2) = xxx * xxx
         darg(2) = 2.0_dp* xxx
         xxx = r - da
         arg(3) = xxx*xxx
         darg(3) = 2.0_dp* xxx
         xxx = r - qa
         arg(4) = xxx*xxx
         darg(4) = 2.0_dp* xxx
         xxx = r + qa
         arg(5) = xxx*xxx
         darg(5) = 2.0_dp* xxx
         arg(6) = rsq
         darg(6)= drsq
         arg(7) = arg(6) + qa*qa
         darg(7)= darg(6)
         xxx = r-db
         arg(8) = xxx*xxx
         darg(8) = 2.0_dp* xxx
         xxx = r+db
         arg(9) = xxx*xxx
         darg(9) = 2.0_dp* xxx
         xxx = r - qb
         arg(10) = xxx*xxx
         darg(10)= 2.0_dp* xxx
         xxx = r + qb
         arg(11) = xxx*xxx
         darg(11)= 2.0_dp* xxx
         arg(12) = rsq
         darg(12)= drsq
         arg(13) = arg(12) + qb*qb
         darg(13)= darg(12)
         xxx = da-db
         arg(14) = rsq + xxx*xxx
         darg(14)= drsq
         xxx = da+db
         arg(15) = rsq + xxx*xxx
         darg(15)= drsq
         xxx = r + da - db
         arg(16) = xxx*xxx
         darg(16)= 2.0_dp* xxx
         xxx = r - da + db
         arg(17) = xxx*xxx
         darg(17)= 2.0_dp* xxx
         xxx = r - da - db
         arg(18) = xxx*xxx
         darg(18)= 2.0_dp* xxx
         xxx = r + da + db
         arg(19) = xxx*xxx
         darg(19)= 2.0_dp* xxx
         xxx = r + da
         arg(20) = xxx*xxx
         darg(20)= 2.0_dp* xxx
         arg(21) = arg(20) + qb*qb
         darg(21)= darg(20)
         xxx = r - da
         arg(22) = xxx*xxx
         darg(22)= 2.0_dp* xxx
         arg(23) = arg(22) + qb*qb
         darg(23)= darg(22)
         xxx = r - db
         arg(24) = xxx*xxx
         darg(24)= 2.0_dp* xxx
         arg(25) = arg(24) + qa*qa
         darg(25)= darg(24)
         xxx = r + db
         arg(26) = xxx*xxx
         darg(26)= 2.0_dp* xxx
         arg(27) = arg(26) + qa*qa
         darg(27)= darg(26)
         xxx = r + da - qb
         arg(28) = xxx*xxx
         darg(28)= 2.0_dp*xxx
         xxx = r - da - qb
         arg(29) = xxx*xxx
         darg(29)= 2.0_dp* xxx
         xxx = r + da + qb
         arg(30) = xxx*xxx
         darg(30)= 2.0_dp* xxx
         xxx = r - da + qb
         arg(31) = xxx*xxx
         darg(31)= 2.0_dp* xxx
         xxx = r + qa - db
         arg(32) = xxx*xxx
         darg(32)= 2.0_dp* xxx
         xxx = r + qa + db
         arg(33) = xxx*xxx
         darg(33)= 2.0_dp* xxx
         xxx = r - qa - db
         arg(34) = xxx*xxx
         darg(34)= 2.0_dp* xxx
         xxx = r - qa + db
         arg(35) = xxx*xxx
         darg(35)= 2.0_dp* xxx
         arg(36) = rsq
         darg(36)= drsq
         xxx = qa - qb
         arg(37) = arg(36) + xxx*xxx
         darg(37)= darg(36)
         xxx = qa + qb
         arg(38) = arg(36) + xxx*xxx
         darg(38)= darg(36)
         arg(39) = arg(36) + qa*qa
         darg(39)= darg(36)
         arg(40) = arg(36) + qb*qb
         darg(40)= darg(36)
         arg(41) = arg(39) + qb*qb
         darg(41)= darg(39)
         xxx = r - qb
         arg(42) = xxx*xxx
         darg(42)= 2.0_dp*xxx
         arg(43) = arg(42) + qa*qa
         darg(43)= darg(42)
         xxx = r + qb
         arg(44) = xxx*xxx
         darg(44)= 2.0_dp*xxx
         arg(45) = arg(44) + qa*qa
         darg(45)= darg(44)
         xxx = r + qa
         arg(46) = xxx*xxx 
         darg(46)= 2.0_dp*xxx
         arg(47) = arg(46) + qb*qb
         darg(47)= darg(46)
         xxx = r - qa
         arg(48) = xxx*xxx
         darg(48)= 2.0_dp*xxx
         arg(49) = arg(48) + qb*qb
         darg(49)= darg(48)
         xxx = r + qa - qb
         arg(50) = xxx*xxx
         darg(50)= 2.0_dp*xxx
         xxx = r + qa + qb
         arg(51) = xxx*xxx
         darg(51)= 2.0_dp*xxx
         xxx = r - qa - qb
         arg(52) = xxx*xxx
         darg(52)= 2.0_dp*xxx
         xxx = r - qa + qb
         arg(53) = xxx*xxx
         darg(53)= 2.0_dp*xxx
         qa=sepi%qq
         qb=sepj%qq
         xxx = da - qb
         dxxx= 0.0_dp
         xxx = xxx*xxx
         yyy = r - qb
         dyyy= 2.0_dp*yyy
         yyy = yyy*yyy
         zzz = da + qb
         dzzz= 0.0_dp
         zzz = zzz*zzz
         www = r + qb
         dwww= 2.0_dp*www
         www = www*www
         arg(54) =  xxx +  yyy
         darg(54)= dxxx + dyyy 
         arg(55) =  xxx +  www
         darg(55)= dxxx + dwww
         arg(56) =  zzz +  yyy
         darg(56)= dzzz + dyyy
         arg(57) =  zzz +  www
         darg(57)= dzzz + dwww
         xxx = qa - db
         dxxx= 0.0_dp
         xxx = xxx*xxx
         yyy = qa + db
         dyyy= 0.0_dp
         yyy = yyy*yyy
         zzz = r + qa
         dzzz= 2.0_dp*zzz
         zzz = zzz*zzz
         www = r - qa
         dwww= 2.0_dp*www
         www = www*www
         arg(58) =  zzz +  xxx
         darg(58)= dzzz + dxxx
         arg(59) =  www +  xxx
         darg(59)= dwww + dxxx
         arg(60) =  zzz +  yyy
         darg(60)= dzzz + dyyy
         arg(61) =  www +  yyy
         darg(61)= dwww + dyyy
         xxx = qa - qb
         xxx = xxx*xxx
         arg(62) =  arg(36) + 2._dp*xxx
         darg(62)= darg(36)
         yyy = qa + qb
         yyy = yyy*yyy
         arg(63) =  arg(36) + 2._dp*yyy
         darg(63)= darg(36)
         arg(64) =  arg(36) + 2._dp*(qa*qa+qb*qb)
         darg(64)= darg(36)
         zzz = r + qa - qb
         dzzz= 2.0_dp*zzz
         zzz = zzz*zzz
         arg(65) = zzz + xxx
         darg(65)= dzzz
         arg(66) = zzz + yyy
         darg(66)= dzzz
         zzz = r + qa + qb
         dzzz= 2.0_dp*zzz
         zzz = zzz*zzz
         arg(67) = zzz + xxx
         darg(67)= dzzz
         arg(68) = zzz + yyy
         darg(68)= dzzz
         zzz = r - qa - qb
         dzzz= 2.0_dp*zzz
         zzz = zzz*zzz
         arg(69) = zzz + xxx
         darg(69)= dzzz
         arg(70) = zzz + yyy
         darg(70)= dzzz
         zzz = r - qa + qb
         dzzz= 2.0_dp*zzz
         zzz = zzz*zzz
         arg(71) = zzz + xxx
         darg(71)= dzzz
         arg(72) = zzz + yyy
         darg(72)= dzzz
         sqr(1:72) = SQRT(arg(1:72))
         dsqr(1:72)= -(0.5_dp*(1/sqr(1:72))**3)*darg(1:72)
  
         ee     =  1._dp/sqr(1)
         dze    = -0.5_dp/sqr(2) + 0.5_dp/sqr(3)
         qzze   =  0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
         qxxe   =  0.5_dp/sqr(7) - 0.5_dp/sqr(6)
         edz    = -0.5_dp/sqr(8) + 0.5_dp/sqr(9)
         eqzz   =  0.25_dp/sqr(10) + 0.25_dp/sqr(11) - 0.5_dp/sqr(12)
         eqxx   =  0.5_dp/sqr(13) - 0.5_dp/sqr(12)
         dxdx   =  0.5_dp/sqr(14) - 0.5_dp/sqr(15)
         dzdz   =  0.25_dp/sqr(16) + 0.25_dp/sqr(17) - 0.25_dp/sqr(18) - 0.25_dp/sqr(19)
         dzqxx  =  0.25_dp/sqr(20) - 0.25_dp/sqr(21) - 0.25_dp/sqr(22) + 0.25_dp/sqr(23)
         qxxdz  =  0.25_dp/sqr(24) - 0.25_dp/sqr(25) - 0.25_dp/sqr(26) + 0.25_dp/sqr(27)
         dzqzz  = -0.125_dp/sqr(28) + 0.125_dp/sqr(29) - 0.125_dp/sqr(30) + 0.125_dp/sqr(31) &
                 - 0.25_dp/sqr(22) + 0.25_dp/sqr(20)
         qzzdz  = -0.125_dp/sqr(32) + 0.125_dp/sqr(33) - 0.125_dp/sqr(34) + 0.125_dp/sqr(35) &
                 + 0.25_dp/sqr(24) - 0.25_dp/sqr(26)
         qxxqxx =  0.125_dp/sqr(37) + 0.125_dp/sqr(38) - 0.25_dp/sqr(39) - 0.25_dp/sqr(40) &
                 + 0.25_dp/sqr(36)
         qxxqyy =  0.25_dp/sqr(41) - 0.25_dp/sqr(39) - 0.25_dp/sqr(40) + 0.25_dp/sqr(36)
         qxxqzz =  0.125_dp/sqr(43) + 0.125_dp/sqr(45) - 0.125_dp/sqr(42) - 0.125_dp/sqr(44) &
                 - 0.25_dp/sqr(39) + 0.25_dp/sqr(36)
         qzzqxx =  0.125_dp/sqr(47) + 0.125_dp/sqr(49) - 0.125_dp/sqr(46) - 0.125_dp/sqr(48) &
                 - 0.25_dp/sqr(40) + 0.25_dp/sqr(36)
         qzzqzz =  0.0625_dp/sqr(50) + 0.0625_dp/sqr(51) + 0.0625_dp/sqr(52) + 0.0625_dp/sqr(53) &
                 - 0.125_dp/sqr(48) - 0.125_dp/sqr(46) - 0.125_dp/sqr(42) - 0.125_dp/sqr(44) + &
                   0.25_dp/sqr(36)
         dxqxz  = -0.25_dp/sqr(54) + 0.25_dp/sqr(55) + 0.25_dp/sqr(56) - 0.25_dp/sqr(57)
         qxzdx  = -0.25_dp/sqr(58) + 0.25_dp/sqr(59) + 0.25_dp/sqr(60) - 0.25_dp/sqr(61)
         qxzqxz =  0.125_dp/sqr(65) - 0.125_dp/sqr(67) - 0.125_dp/sqr(69) + 0.125_dp/sqr(71) &
                 - 0.125_dp/sqr(66) + 0.125_dp/sqr(68) + 0.125_dp/sqr(70) - 0.125_dp/sqr(72)
  
         dee     =  1._dp*dsqr(1)
         ddze    = -0.5_dp*dsqr(2) + 0.5_dp*dsqr(3)
         dqzze   =  0.25_dp*dsqr(4) + 0.25_dp*dsqr(5) - 0.5_dp*dsqr(6)
         dqxxe   =  0.5_dp*dsqr(7) - 0.5_dp*dsqr(6)
         dedz    = -0.5_dp*dsqr(8) + 0.5_dp*dsqr(9)
         deqzz   =  0.25_dp*dsqr(10) + 0.25_dp*dsqr(11) - 0.5_dp*dsqr(12)
         deqxx   =  0.5_dp*dsqr(13) - 0.5_dp*dsqr(12)
         ddxdx   =  0.5_dp*dsqr(14) - 0.5_dp*dsqr(15)
         ddzdz   =  0.25_dp*dsqr(16) + 0.25_dp*dsqr(17) - 0.25_dp*dsqr(18) - 0.25_dp*dsqr(19)
         ddzqxx  =  0.25_dp*dsqr(20) - 0.25_dp*dsqr(21) - 0.25_dp*dsqr(22) + 0.25_dp*dsqr(23)
         dqxxdz  =  0.25_dp*dsqr(24) - 0.25_dp*dsqr(25) - 0.25_dp*dsqr(26) + 0.25_dp*dsqr(27)
         ddzqzz  = -0.125_dp*dsqr(28) + 0.125_dp*dsqr(29) - 0.125_dp*dsqr(30) + 0.125_dp*dsqr(31) &
                  - 0.25_dp*dsqr(22) + 0.25_dp*dsqr(20)
         dqzzdz  = -0.125_dp*dsqr(32) + 0.125_dp*dsqr(33) - 0.125_dp*dsqr(34) + 0.125_dp*dsqr(35) &
                  + 0.25_dp*dsqr(24) - 0.25_dp*dsqr(26)
         dqxxqxx =  0.125_dp*dsqr(37) + 0.125_dp*dsqr(38) - 0.25_dp*dsqr(39) - 0.25_dp*dsqr(40) &
                  + 0.25_dp*dsqr(36)
         dqxxqyy =  0.25_dp*dsqr(41) - 0.25_dp*dsqr(39) - 0.25_dp*dsqr(40) + 0.25_dp*dsqr(36)
         dqxxqzz =  0.125_dp*dsqr(43) + 0.125_dp*dsqr(45) - 0.125_dp*dsqr(42) - 0.125_dp*dsqr(44) &
                  - 0.25_dp*dsqr(39) + 0.25_dp*dsqr(36)
         dqzzqxx =  0.125_dp*dsqr(47) + 0.125_dp*dsqr(49) - 0.125_dp*dsqr(46) - 0.125_dp*dsqr(48) &
                  - 0.25_dp*dsqr(40) + 0.25_dp*dsqr(36)
         dqzzqzz =  0.0625_dp*dsqr(50) + 0.0625_dp*dsqr(51) + 0.0625_dp*dsqr(52) + 0.0625_dp*dsqr(53) &
                  - 0.125_dp*dsqr(48) - 0.125_dp*dsqr(46) - 0.125_dp*dsqr(42) - 0.125_dp*dsqr(44) + &
                    0.25_dp*dsqr(36)
         ddxqxz  = -0.25_dp*dsqr(54) + 0.25_dp*dsqr(55) + 0.25_dp*dsqr(56) - 0.25_dp*dsqr(57)
         dqxzdx  = -0.25_dp*dsqr(58) + 0.25_dp*dsqr(59) + 0.25_dp*dsqr(60) - 0.25_dp*dsqr(61)
         dqxzqxz =  0.125_dp*dsqr(65) - 0.125_dp*dsqr(67) - 0.125_dp*dsqr(69) + 0.125_dp*dsqr(71) &
                 - 0.125_dp*dsqr(66) + 0.125_dp*dsqr(68) + 0.125_dp*dsqr(70) - 0.125_dp*dsqr(72)
  
         ris(1) = ee
         ris(2) = -dze
         ris(3) = ee + qzze
         ris(4) = ee + qxxe
         ris(5) = -edz
         ris(6) = dzdz
         ris(7) = dxdx
         ris(8) = -edz -qzzdz
         ris(9) = -edz -qxxdz
         ris(10) = -qxzdx
         ris(11) =  ee + eqzz
         ris(12) =  ee + eqxx
         ris(13) = -dze -dzqzz
         ris(14) = -dze -dzqxx
         ris(15) = -dxqxz
         ris(16) = ee +eqzz +qzze +qzzqzz
         ris(17) = ee +eqzz +qxxe +qxxqzz
         ris(18) = ee +eqxx +qzze +qzzqxx
         ris(19) = ee +eqxx +qxxe +qxxqxx
         ris(20) = qxzqxz
         ris(21) = ee +eqxx +qxxe +qxxqyy
         ris(22) = 0.5_dp * (qxxqxx -qxxqyy)
  
         ri(1) = ri(1) - ris(1)
         ri(2) = ri(2) - ris(2)
         ri(3) = ri(3) - ris(3)
         ri(4) = ri(4) - ris(4)
         ri(5) = ri(5) - ris(5)
         ri(6) = ri(6) - ris(6)
         ri(7) = ri(7) - ris(7)
         ri(8) = ri(8) - ris(8)
         ri(9) = ri(9) - ris(9)
         ri(10) = ri(10) - ris(10)
         ri(11) = ri(11) - ris(11)
         ri(12) = ri(12) - ris(12)
         ri(13) = ri(13) - ris(13)
         ri(14) = ri(14) - ris(14)
         ri(15) = ri(15) - ris(15)
         ri(16) = ri(16) - ris(16)
         ri(17) = ri(17) - ris(17)
         ri(18) = ri(18) - ris(18)
         ri(19) = ri(19) - ris(19)
         ri(20) = ri(20) - ris(20)
         ri(21) = ri(21) - ris(21)
         ri(22) = ri(22) - ris(22)
  
         dri(1) = dri(1) - (dee)
         dri(2) = dri(1) - (-ddze)
         dri(3) = dri(1) - (dee + dqzze)
         dri(4) = dri(1) - (dee + dqxxe)
         dri(5) = dri(1) - (-dedz)
         dri(6) = dri(1) - (ddzdz)
         dri(7) = dri(1) - (ddxdx)
         dri(8) = dri(1) - (-dedz -dqzzdz)
         dri(9) = dri(1) - (-dedz -dqxxdz)
         dri(10) = dri(1) - (-dqxzdx)
         dri(11) = dri(1) - ( dee + deqzz)
         dri(12) = dri(1) - ( dee + deqxx)
         dri(13) = dri(1) - (-ddze -ddzqzz)
         dri(14) = dri(1) - (-ddze -ddzqxx)
         dri(15) = dri(1) - (-ddxqxz)
         dri(16) = dri(1) - (dee +deqzz +dqzze +dqzzqzz)
         dri(17) = dri(1) - (dee +deqzz +dqxxe +dqxxqzz)
         dri(18) = dri(1) - (dee +deqxx +dqzze +dqzzqxx)
         dri(19) = dri(1) - (dee +deqxx +dqxxe +dqxxqxx)
         dri(20) = dri(1) - (dqxzqxz)
         dri(21) = dri(1) - (dee +deqxx +dqxxe +dqxxqyy)
         dri(22) = dri(1) - (0.5_dp * (dqxxqxx -dqxxqyy))
   
      END IF
  
    END IF

    ! Tapering function
    ft = taper  (rij)
    dft= dtaper_ana (rij)
    ri(:) = ft*ri(:)
    dri(:) = dft*ri(:)+ft*dri(:)

  END SUBROUTINE dterep_ana

!!****f* semi_empirical_int_ana/check_dterep_ana [1.0] *
!!
!!   NAME
!!     check_dterep_ana
!!
!!   FUNCTION
!!     Check Numerical Vs Analytical
!!
!!   NOTES
!!     Debug routine
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino - Zurich University
!!
!!   MODIFICATION HISTORY
!!     04.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE check_dterep_ana (sepi,sepj,r,ri,dri,shortrange)

    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: r
    REAL(dp), DIMENSION(:), INTENT(IN)       :: ri, dri
    LOGICAL, INTENT(IN)                      :: shortrange

    INTEGER                                  :: j
    REAL(dp)                                 :: delta, od, rn
    REAL(dp), DIMENSION(22)                  :: nri, ri0, rim, rip

    delta = 1.0E-8_dp
    od = 0.5_dp/delta
    rn = r
    CALL terep(sepi,sepj,rn,ri0,shortrange)
    rn = r + delta
    CALL terep(sepi,sepj,rn,rip,shortrange)
    rn = r - delta
    CALL terep(sepi,sepj,rn,rim,shortrange)
    nri = od * (rip - rim)
    ! check
    DO j = 1, 22
       IF  (ABS(ri(j)-ri0(j))>EPSILON(0.0_dp)) THEN
          WRITE(*,*)"Error in value of the integral.."
       END IF
       IF      ((ABS(nri(j))>EPSILON(0.0_dp)).AND.(ABS(dri(j))>EPSILON(0.0_dp))) THEN
          IF ((nri(j)-dri(j))/nri(j)*100._dp>0.1_dp) THEN
             WRITE(*,*)"Error in the value of the derivative::",j,ABS(nri(j)-dri(j))/nri(j)*100._dp
             STOP
          END IF
       ELSEIF  ((ABS(nri(j))<EPSILON(0.0_dp)).AND.(ABS(dri(j))>EPSILON(0.0_dp))) THEN
          WRITE(*,*) j,ABS(nri(j)-dri(j))/nri(j)*100._dp
          STOP
       ELSEIF  ((ABS(nri(j))>EPSILON(0.0_dp)).AND.(ABS(dri(j))<EPSILON(0.0_dp))) THEN
          WRITE(*,*) j,ABS(nri(j)-dri(j))/nri(j)*100._dp
          STOP
       END IF
    END DO

  END SUBROUTINE check_dterep_ana

END MODULE semi_empirical_int_ana
