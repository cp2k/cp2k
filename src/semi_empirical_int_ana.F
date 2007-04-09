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
!!     Analytical derivatives of Integrals for semi-empiric methods
!!
!!   AUTHOR
!!     Teodoro Laino 04.2007
!!
!!   MODIFICATION HISTORY
!!     
!!***
!******************************************************************************

MODULE semi_empirical_int_ana

  USE kinds,                           ONLY: dp
  USE semi_empirical_integrals,        ONLY: al,&
                                             drotnuc,&
                                             nucint,&
                                             r0,&
                                             rotnuc,&
                                             taper,&
                                             taper_fn_init
  USE semi_empirical_types,            ONLY: semi_empirical_type

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_int_ana'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module=.FALSE.
  PUBLIC :: rotnuc_ana

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

  RECURSIVE SUBROUTINE rotnuc_ana (sepi,sepj,rijv,itype,enuc,e1b,e2a,denuc,de1b,de2a)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rijv
    INTEGER, INTENT(IN)                      :: itype
    REAL(dp), INTENT(OUT), OPTIONAL          :: enuc
    REAL(dp), DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: e1b, e2a, denuc
    REAL(dp), DIMENSION(:, :), INTENT(OUT), &
      OPTIONAL                               :: de1b, de2a

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
    IF (rij > 0.00002_dp) THEN
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
       ! Compute Integrals in diatomic frame opportunely inverted
       rij  = SQRT(rij)
       a    = 1._dp/rij
       IF (lgrad.OR.l_denuc) THEN
          CALL dnucint_ana(sepi,sepj,rij,ssss,ccore,dssss,dccore,itype)
          ! Debug Procedure.. Check valifity of analytical gradients of nucint
          IF (debug_this_module) CALL check_dnucint_ana(sepi,sepj,rij,dssss,dccore,itype)
       ELSE
          CALL nucint(sepi,sepj,rij,ssss,ccore,itype)
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
          IF (invert) CALL invert_integral(e1b)
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
          IF (invert) CALL invert_integral(e2a)
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
             IF (invert) CALL invert_derivative(de1b)
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
             IF (invert) CALL invert_derivative(de2a)
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
          scale=ABS(scale*enuc_loc)
          IF (l_denuc) THEN
             denuc_loc= zz*dssss
             dscale=SIGN(1.0_dp,scale*enuc_loc)*(dscale*enuc_loc+scale*denuc_loc)
             dzz=-zz/rij**2
          END IF
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
          CALL rotnuc (sepi,sepj,rijv,e1b2,e2a2,enuc2,itype)
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
             CALL drotnuc(sepi,sepj,rijv,de1b2,de2a2,denuc2,itype,delta=delta)
             CALL rotnuc_ana(sepi,sepj,rijv,itype,e1b=e1b2,e2a=e2a2)
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
  SUBROUTINE invert_integral(array)
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: array

    REAL(KIND=dp)                            :: tmp

    tmp=array(4)
    array(4)=array(7)
    array(7)=tmp
    
    tmp=array(5)
    array(5)=array(8)
    array(8)=tmp
    
    tmp=array(6)
    array(6)=array(10)
    array(10)=tmp

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
  SUBROUTINE invert_derivative(array)
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: array

    INTEGER                                  :: j, m
    REAL(KIND=dp)                            :: tmp

    DO j=1,3
       tmp=array(4,j)
       array(4,j)=array(7,j)
       array(7,j)=tmp
       
       tmp=array(5,j)
       array(5,j)=array(8,j)
       array(8,j)=tmp
       
       tmp=array(6,j)
       array(6,j)=array(10,j)
       array(10,j)=tmp
    END DO

    DO m=1,10
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
  SUBROUTINE check_dnucint_ana (sepi,sepj,r,dssss,dccore,itype)

    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: r, dssss
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: dccore
    INTEGER, INTENT(IN)                      :: itype

    INTEGER                                  :: i, j
    REAL(dp)                                 :: delta, nssss, od, rn, ssssm, &
                                                ssssp
    REAL(dp), DIMENSION(4, 2)                :: ccorem, ccorep, nccore

    delta = 1.0E-8_dp
    od = 0.5_dp/delta
    rn = r + delta
    CALL nucint(sepi,sepj,rn,ssssp,ccorep,itype)
    rn = r - delta
    CALL nucint(sepi,sepj,rn,ssssm,ccorem,itype)
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
  SUBROUTINE dnucint_ana ( sepi, sepj, rij, ssss, core, dssss, dcore, itype )
    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), INTENT(OUT)                    :: ssss
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: core
    REAL(dp), INTENT(OUT)                    :: dssss
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: dcore
    INTEGER, INTENT(IN)                      :: itype

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
       dri(1) = 1._dp/SQRT(r*r+aee)*fac
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
    ssss = ri(1)
    ! Tapering function
    ft = taper ( rij )
    core(:,:) = ft*core(:,:)
    ssss = ft*ssss
    ! If no orbitals are present on the atom let's zero the integral contribution
    IF (sepi%natorb == -1) core(:,1) = 0.0_dp
    IF (sepj%natorb == -1) core(:,2) = 0.0_dp
    !
    dssss = dri(1)
    ! Tapering function
    dft= dtaper_ana ( rij )
    dcore(:,:) = ft*dcore(:,:)+dft*core(:,:)
    dssss = ft*dssss+dft*ssss
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

END MODULE semi_empirical_int_ana
