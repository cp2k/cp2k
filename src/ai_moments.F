!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005 CP2K developers group                                  !
!-----------------------------------------------------------------------------!
!!****** cp2k/ai_moments [1.0] *
!!
!!   NAME
!!     ai_moments
!!
!!   FUNCTION
!!     Calculation of the moment integrals over Cartesian Gaussian-type
!!     functions.
!!
!!   AUTHOR
!!     J. Hutter (16.02.2005)
!!
!!   LITERATURE
!!     S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!******************************************************************************

MODULE ai_moments

! *****************************************************************************

! ax,ay,az  : Angular momentum index numbers of orbital a.
! bx,by,bz  : Angular momentum index numbers of orbital b.
! coset     : Cartesian orbital set pointer.
! dab       : Distance between the atomic centers a and b.
! l{a,b}    : Angular momentum quantum number of shell a or b.
! l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
! l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
! rac       : Distance vector between the atomic center a and reference point c.
! rbc       : Distance vector between the atomic center b and reference point c.
! rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
! zet{a,b}  : Exponents of the Gaussian-type functions a or b.
! zetp      : Reciprocal of the sum of the exponents of orbital a and b.

! *****************************************************************************

  USE ai_derivatives,                  ONLY: dabdr
  USE ai_overlap_new,                  ONLY: overlap
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,&
                                             indco,&
                                             ncoset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cossin, moment, diffop, contract_cossin

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE contract_cossin(cos_block, sin_block,&
                             iatom,ncoa,nsgfa,sgfa,sphi_a,ldsa,&
                             jatom,ncob,nsgfb,sgfb,sphi_b,ldsb,&
                             cosab,sinab,ldab,work,ldwork)

    REAL(dp), DIMENSION(:,:), POINTER       :: cos_block, sin_block
    INTEGER, INTENT(IN)                     :: iatom,ncoa,nsgfa,sgfa
    REAL(dp), DIMENSION(:,:), INTENT(IN)    :: sphi_a
    INTEGER, INTENT(IN)                     :: ldsa,jatom,ncob,nsgfb,sgfb
    REAL(dp), DIMENSION(:,:), INTENT(IN)    :: sphi_b
    REAL(dp), DIMENSION(:,:), INTENT(IN)    :: cosab, sinab
    INTEGER, INTENT(IN)                     :: ldsb,ldab
    REAL(dp), DIMENSION(:,:)                :: work
    INTEGER, INTENT(IN)                     :: ldwork


       ! Calculate cosine
       CALL dgemm("N","N",ncoa,nsgfb,ncob,&
                   1.0_dp,cosab(1,1),ldab,&
                   sphi_b(1,sgfb),ldsb,&
                   0.0_dp,work(1,1),ldwork)       
   
       IF(iatom<=jatom) THEN
         CALL dgemm("T","N",nsgfa,nsgfb,ncoa,&
                    1.0_dp,sphi_a(1,sgfa),ldsa,&
                    work(1,1),ldwork,&
                    1.0_dp,cos_block(sgfa,sgfb),&
                    SIZE(cos_block,1))
       ELSE
         CALL dgemm("T","N",nsgfb,nsgfa,ncoa,&
                    1.0_dp,work(1,1),ldwork,&
                    sphi_a(1,sgfa),ldsa,&
                    1.0_dp,cos_block(sgfb,sgfa),&
                    SIZE(cos_block,1))
       END IF

       ! Calculate sine
       CALL dgemm("N","N",ncoa,nsgfb,ncob,&
                   1.0_dp,sinab(1,1),ldab,&
                   sphi_b(1,sgfb),ldsb,&
                   0.0_dp,work(1,1),ldwork)       

       IF(iatom<=jatom) THEN
         CALL dgemm("T","N",nsgfa,nsgfb,ncoa,&
                    1.0_dp,sphi_a(1,sgfa),ldsa,&
                    work(1,1),ldwork,&
                    1.0_dp,sin_block(sgfa,sgfb),&
                    SIZE(sin_block,1))
       ELSE
         CALL dgemm("T","N",nsgfb,nsgfa,ncoa,&
                    1.0_dp,work(1,1),ldwork,&
                    sphi_a(1,sgfa),ldsa,&
                    1.0_dp,sin_block(sgfb,sgfa),&
                    SIZE(sin_block,1))
       END IF

    
  END SUBROUTINE contract_cossin

! *****************************************************************************

  SUBROUTINE cossin(la_max,npgfa,zeta,rpgfa,la_min,&
                    lb_max,npgfb,zetb,rpgfb,lb_min,&
                    rac,rbc,kvec,cosab,sinab)

    INTEGER, INTENT(IN)                      :: la_max,la_min,npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta, rpgfa
    INTEGER, INTENT(IN)                      :: lb_max,lb_min,npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb, rpgfb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rac, rbc
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: kvec
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: cosab,sinab

    INTEGER                                  :: ax, ay, az, bx, by, bz, i, &
                                                ipgf, j, jpgf, la, la_start, &
                                                lb, na, nb
    REAL(KIND=dp)                            :: dab, f0, f1, f2, f3, fx, fy, &
                                                fz, zetp, k2, kdp, rab2, s
    REAL(KIND=dp), DIMENSION(3)              :: rab, rap, rbp
    REAL(KIND=dp), DIMENSION(ncoset(la_max),&
      ncoset(lb_max))                        :: sc, ss

!   ---------------------------------------------------------------------------

    rab = rbc - rac
    rab2= SUM( rab**2 )
    dab = SQRT(rab2)
    k2  = kvec(1)*kvec(1) + kvec(2)*kvec(2) + kvec(3)*kvec(3)

!   *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0

    DO ipgf=1,npgfa

      nb = 0

      DO jpgf=1,npgfb

!       *** Screening ***

        IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
          DO j=nb+1,nb+ncoset(lb_max)
            DO i=na+1,na+ncoset(la_max)
              cosab(i,j) = 0.0_dp
              sinab(i,j) = 0.0_dp
            END DO
          END DO
          nb = nb + ncoset(lb_max)
          CYCLE
        END IF

!       *** Calculate some prefactors ***

        zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))

        f0 = (pi*zetp)**1.5_dp
        f1 = zetb(jpgf)*zetp
        f2 = 0.5_dp*zetp

        kdp = zetp * DOT_PRODUCT ( kvec, zeta(ipgf)*rac + zetb(jpgf)*rbc ) 


!       *** Calculate the basic two-center cos/sin integral [s|cos/sin|s] ***

        s       = f0*EXP(-zeta(ipgf)*f1*rab2)*EXP(-0.25_dp*k2*zetp)
        sc(1,1) = s*COS(kdp)
        ss(1,1) = s*SIN(kdp)

!       *** Recurrence steps: [s|O|s] -> [a|O|b] ***

        IF (la_max > 0) THEN

!         *** Vertical recurrence steps: [s|O|s] -> [a|O|s] ***

          rap(:) = f1*rab(:)

!         *** [p|O|s] = (Pi - Ai)*[s|O|s] +[s|dO|s]  (i = x,y,z) ***

          sc(2,1) = rap(1)*sc(1,1) - f2* kvec(1)*ss(1,1)
          sc(3,1) = rap(2)*sc(1,1) - f2* kvec(2)*ss(1,1)
          sc(4,1) = rap(3)*sc(1,1) - f2* kvec(3)*ss(1,1)
          ss(2,1) = rap(1)*ss(1,1) + f2* kvec(1)*sc(1,1)
          ss(3,1) = rap(2)*ss(1,1) + f2* kvec(2)*sc(1,1)
          ss(4,1) = rap(3)*ss(1,1) + f2* kvec(3)*sc(1,1)

!         *** [a|O|s] = (Pi - Ai)*[a-1i|O|s] + f2*Ni(a-1i)*[a-2i|s] ***
!         ***           + [a-1i|dO|s]                               ***

          DO la=2,la_max

!           *** Increase the angular momentum component z of function a ***

            sc(coset(0,0,la),1) = rap(3)*sc(coset(0,0,la-1),1) +&
                                 f2*REAL(la-1,dp)*sc(coset(0,0,la-2),1) -&
                            f2 *     kvec(3)*ss(coset(0,0,la-1),1)
            ss(coset(0,0,la),1) = rap(3)*ss(coset(0,0,la-1),1) +&
                                 f2*REAL(la-1,dp)*ss(coset(0,0,la-2),1) +&
                            f2 *     kvec(3)*sc(coset(0,0,la-1),1)

!           *** Increase the angular momentum component y of function a ***

            az = la - 1
            sc(coset(0,1,az),1) = rap(2)*sc(coset(0,0,az),1) -&
                            f2 *      kvec(2)*ss(coset(0,0,az),1)
            ss(coset(0,1,az),1) = rap(2)*ss(coset(0,0,az),1) +&
                            f2 *      kvec(2)*sc(coset(0,0,az),1)

            DO ay=2,la
              az = la - ay
              sc(coset(0,ay,az),1) = rap(2)*sc(coset(0,ay-1,az),1) +&
                                     f2*REAL(ay-1,dp)*sc(coset(0,ay-2,az),1) -&
                            f2 *         kvec(2)*ss(coset(0,ay-1,az),1)
              ss(coset(0,ay,az),1) = rap(2)*ss(coset(0,ay-1,az),1) +&
                                     f2*REAL(ay-1,dp)*ss(coset(0,ay-2,az),1) +&
                            f2 *         kvec(2)*sc(coset(0,ay-1,az),1)
            END DO

!           *** Increase the angular momentum component x of function a ***

            DO ay=0,la-1
              az = la - 1 - ay
              sc(coset(1,ay,az),1) = rap(1)*sc(coset(0,ay,az),1) -&
                            f2 *         kvec(1)*ss(coset(0,ay,az),1)
              ss(coset(1,ay,az),1) = rap(1)*ss(coset(0,ay,az),1) +&
                            f2 *         kvec(1)*sc(coset(0,ay,az),1)
            END DO

            DO ax=2,la
              f3 = f2*REAL(ax-1,dp)
              DO ay=0,la-ax
                az = la - ax - ay
                sc(coset(ax,ay,az),1) = rap(1)*sc(coset(ax-1,ay,az),1) +&
                                        f3*sc(coset(ax-2,ay,az),1) -&
                              f2 *          kvec(1)*ss(coset(ax-1,ay,az),1)
                ss(coset(ax,ay,az),1) = rap(1)*ss(coset(ax-1,ay,az),1) +&
                                        f3*ss(coset(ax-2,ay,az),1) +&
                              f2 *          kvec(1)*sc(coset(ax-1,ay,az),1)
              END DO
            END DO

          END DO

!         *** Recurrence steps: [a|O|s] -> [a|O|b] ***

          IF (lb_max > 0) THEN

            DO j=2,ncoset(lb_max)
              DO i=1,ncoset(la_max)
                sc(i,j) = 0.0_dp
                ss(i,j) = 0.0_dp
              END DO
            END DO

!           *** Horizontal recurrence steps ***

            rbp(:) = rap(:) - rab(:)

!           *** [a|O|p] = [a+1i|O|s] - (Bi - Ai)*[a|O|s] ***

            IF (lb_max == 1) THEN
              la_start = la_min
            ELSE
              la_start = MAX(0,la_min-1)
            END IF

            DO la=la_start,la_max-1
              DO ax=0,la
                DO ay=0,la-ax
                  az = la - ax - ay
                  sc(coset(ax,ay,az),2) = sc(coset(ax+1,ay,az),1) -&
                                         rab(1)*sc(coset(ax,ay,az),1)
                  sc(coset(ax,ay,az),3) = sc(coset(ax,ay+1,az),1) -&
                                         rab(2)*sc(coset(ax,ay,az),1)
                  sc(coset(ax,ay,az),4) = sc(coset(ax,ay,az+1),1) -&
                                         rab(3)*sc(coset(ax,ay,az),1)
                  ss(coset(ax,ay,az),2) = ss(coset(ax+1,ay,az),1) -&
                                         rab(1)*ss(coset(ax,ay,az),1)
                  ss(coset(ax,ay,az),3) = ss(coset(ax,ay+1,az),1) -&
                                         rab(2)*ss(coset(ax,ay,az),1)
                  ss(coset(ax,ay,az),4) = ss(coset(ax,ay,az+1),1) -&
                                         rab(3)*ss(coset(ax,ay,az),1)
                END DO
              END DO
            END DO

!           *** Vertical recurrence step ***

!           *** [a|O|p] = (Pi - Bi)*[a|O|s] + f2*Ni(a)*[a-1i|O|s] ***
!           ***           + [a|dO|s]                              ***

            DO ax=0,la_max
              fx = f2*REAL(ax,dp)
              DO ay=0,la_max-ax
                fy = f2*REAL(ay,dp)
                az = la_max - ax - ay
                fz = f2*REAL(az,dp)
                IF (ax == 0) THEN
                  sc(coset(ax,ay,az),2) = rbp(1)*sc(coset(ax,ay,az),1) -&
                                 f2 *         kvec(1)*ss(coset(ax,ay,az),1)
                  ss(coset(ax,ay,az),2) = rbp(1)*ss(coset(ax,ay,az),1) +&
                                 f2 *         kvec(1)*sc(coset(ax,ay,az),1)
                ELSE
                  sc(coset(ax,ay,az),2) = rbp(1)*sc(coset(ax,ay,az),1) +&
                                          fx*sc(coset(ax-1,ay,az),1) -&
                                 f2 *         kvec(1)*ss(coset(ax,ay,az),1)
                  ss(coset(ax,ay,az),2) = rbp(1)*ss(coset(ax,ay,az),1) +&
                                          fx*ss(coset(ax-1,ay,az),1) +&
                                 f2 *         kvec(1)*sc(coset(ax,ay,az),1)
                END IF
                IF (ay == 0) THEN
                  sc(coset(ax,ay,az),3) = rbp(2)*sc(coset(ax,ay,az),1) -&
                                  f2 *        kvec(2)*ss(coset(ax,ay,az),1)
                  ss(coset(ax,ay,az),3) = rbp(2)*ss(coset(ax,ay,az),1) +&
                                 f2 *         kvec(2)*sc(coset(ax,ay,az),1)
                ELSE
                  sc(coset(ax,ay,az),3) = rbp(2)*sc(coset(ax,ay,az),1) +&
                                          fy*sc(coset(ax,ay-1,az),1) -&
                                  f2 *        kvec(2)*ss(coset(ax,ay,az),1)
                  ss(coset(ax,ay,az),3) = rbp(2)*ss(coset(ax,ay,az),1) +&
                                          fy*ss(coset(ax,ay-1,az),1) +&
                                  f2 *        kvec(2)*sc(coset(ax,ay,az),1)
                END IF
                IF (az == 0) THEN
                  sc(coset(ax,ay,az),4) = rbp(3)*sc(coset(ax,ay,az),1) -&
                                  f2 *        kvec(3)*ss(coset(ax,ay,az),1)
                  ss(coset(ax,ay,az),4) = rbp(3)*ss(coset(ax,ay,az),1) +&
                                  f2 *        kvec(3)*sc(coset(ax,ay,az),1)
                ELSE
                  sc(coset(ax,ay,az),4) = rbp(3)*sc(coset(ax,ay,az),1) +&
                                          fz*sc(coset(ax,ay,az-1),1) -&
                                  f2 *        kvec(3)*ss(coset(ax,ay,az),1)
                  ss(coset(ax,ay,az),4) = rbp(3)*ss(coset(ax,ay,az),1) +&
                                          fz*ss(coset(ax,ay,az-1),1) +&
                                  f2 *        kvec(3)*sc(coset(ax,ay,az),1)
                END IF
              END DO
            END DO

!           *** Recurrence steps: [a|O|p] -> [a|O|b] ***

            DO lb=2,lb_max

!             *** Horizontal recurrence steps ***

!             *** [a|O|b] = [a+1i|O|b-1i] - (Bi - Ai)*[a|O|b-1i] ***

              IF (lb == lb_max) THEN
                la_start = la_min
              ELSE
                la_start = MAX(0,la_min-1)
              END IF

              DO la=la_start,la_max-1
                DO ax=0,la
                  DO ay=0,la-ax
                    az = la - ax - ay

!                   *** Shift of angular momentum component z from a to b ***

                    sc(coset(ax,ay,az),coset(0,0,lb)) =&
                      sc(coset(ax,ay,az+1),coset(0,0,lb-1)) -&
                      rab(3)*sc(coset(ax,ay,az),coset(0,0,lb-1))
                    ss(coset(ax,ay,az),coset(0,0,lb)) =&
                      ss(coset(ax,ay,az+1),coset(0,0,lb-1)) -&
                      rab(3)*ss(coset(ax,ay,az),coset(0,0,lb-1))

!                   *** Shift of angular momentum component y from a to b ***

                    DO by=1,lb
                      bz = lb - by
                      sc(coset(ax,ay,az),coset(0,by,bz)) =&
                        sc(coset(ax,ay+1,az),coset(0,by-1,bz)) -&
                        rab(2)*sc(coset(ax,ay,az),coset(0,by-1,bz))
                      ss(coset(ax,ay,az),coset(0,by,bz)) =&
                        ss(coset(ax,ay+1,az),coset(0,by-1,bz)) -&
                        rab(2)*ss(coset(ax,ay,az),coset(0,by-1,bz))
                    END DO

!                   *** Shift of angular momentum component x from a to b ***

                    DO bx=1,lb
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        sc(coset(ax,ay,az),coset(bx,by,bz)) =&
                          sc(coset(ax+1,ay,az),coset(bx-1,by,bz)) -&
                          rab(1)*sc(coset(ax,ay,az),coset(bx-1,by,bz))
                        ss(coset(ax,ay,az),coset(bx,by,bz)) =&
                          ss(coset(ax+1,ay,az),coset(bx-1,by,bz)) -&
                          rab(1)*ss(coset(ax,ay,az),coset(bx-1,by,bz))
                      END DO
                    END DO

                  END DO
                END DO
              END DO

!             *** Vertical recurrence step ***

!             *** [a|O|b] = (Pi - Bi)*[a|O|b-1i] + f2*Ni(a)*[a-1i|O|b-1i] + ***
!             ***           f2*Ni(b-1i)*[a|O|b-2i] + [a|dO|b-1i]            ***

              DO ax=0,la_max
                fx = f2*REAL(ax,dp)
                DO ay=0,la_max-ax
                  fy = f2*REAL(ay,dp)
                  az = la_max - ax - ay
                  fz = f2*REAL(az,dp)

!                 *** Increase the angular momentum component z of function b ***

                  f3 = f2*REAL(lb-1,dp)

                  IF (az == 0) THEN
                    sc(coset(ax,ay,az),coset(0,0,lb)) =&
                      rbp(3)*sc(coset(ax,ay,az),coset(0,0,lb-1)) +&
                      f3*sc(coset(ax,ay,az),coset(0,0,lb-2)) -&
               f2 *       kvec(3)*ss(coset(ax,ay,az),coset(0,0,lb-1))
                    ss(coset(ax,ay,az),coset(0,0,lb)) =&
                      rbp(3)*ss(coset(ax,ay,az),coset(0,0,lb-1)) +&
                      f3*ss(coset(ax,ay,az),coset(0,0,lb-2)) +&
               f2 *       kvec(3)*sc(coset(ax,ay,az),coset(0,0,lb-1))
                  ELSE
                    sc(coset(ax,ay,az),coset(0,0,lb)) =&
                      rbp(3)*sc(coset(ax,ay,az),coset(0,0,lb-1)) +&
                      fz*sc(coset(ax,ay,az-1),coset(0,0,lb-1)) +&
                      f3*sc(coset(ax,ay,az),coset(0,0,lb-2)) -&
               f2 *       kvec(3)*ss(coset(ax,ay,az),coset(0,0,lb-1))
                    ss(coset(ax,ay,az),coset(0,0,lb)) =&
                      rbp(3)*ss(coset(ax,ay,az),coset(0,0,lb-1)) +&
                      fz*ss(coset(ax,ay,az-1),coset(0,0,lb-1)) +&
                      f3*ss(coset(ax,ay,az),coset(0,0,lb-2)) +&
               f2 *       kvec(3)*sc(coset(ax,ay,az),coset(0,0,lb-1))
                  END IF

!                 *** Increase the angular momentum component y of function b ***

                  IF (ay == 0) THEN
                    bz = lb - 1
                    sc(coset(ax,ay,az),coset(0,1,bz)) =&
                      rbp(2)*sc(coset(ax,ay,az),coset(0,0,bz)) -&
                f2 *      kvec(2)*ss(coset(ax,ay,az),coset(0,0,bz))
                    ss(coset(ax,ay,az),coset(0,1,bz)) =&
                      rbp(2)*ss(coset(ax,ay,az),coset(0,0,bz)) +&
                f2 *      kvec(2)*sc(coset(ax,ay,az),coset(0,0,bz))
                    DO by=2,lb
                      bz = lb - by
                      f3 = f2*REAL(by-1,dp)
                      sc(coset(ax,ay,az),coset(0,by,bz)) =&
                        rbp(2)*sc(coset(ax,ay,az),coset(0,by-1,bz)) +&
                        f3*sc(coset(ax,ay,az),coset(0,by-2,bz)) -&
                 f2 *       kvec(2)*ss(coset(ax,ay,az),coset(0,by-1,bz))
                      ss(coset(ax,ay,az),coset(0,by,bz)) =&
                        rbp(2)*ss(coset(ax,ay,az),coset(0,by-1,bz)) +&
                        f3*ss(coset(ax,ay,az),coset(0,by-2,bz)) +&
                  f2 *      kvec(2)*sc(coset(ax,ay,az),coset(0,by-1,bz))
                    END DO
                  ELSE
                    bz = lb - 1
                    sc(coset(ax,ay,az),coset(0,1,bz)) =&
                      rbp(2)*sc(coset(ax,ay,az),coset(0,0,bz)) +&
                      fy*sc(coset(ax,ay-1,az),coset(0,0,bz)) -&
                 f2 *     kvec(2)*ss(coset(ax,ay,az),coset(0,0,bz))
                    ss(coset(ax,ay,az),coset(0,1,bz)) =&
                      rbp(2)*ss(coset(ax,ay,az),coset(0,0,bz)) +&
                      fy*ss(coset(ax,ay-1,az),coset(0,0,bz)) +&
                 f2 *     kvec(2)*sc(coset(ax,ay,az),coset(0,0,bz))
                    DO by=2,lb
                      bz = lb - by
                      f3 = f2*REAL(by-1,dp)
                      sc(coset(ax,ay,az),coset(0,by,bz)) =&
                        rbp(2)*sc(coset(ax,ay,az),coset(0,by-1,bz)) +&
                        fy*sc(coset(ax,ay-1,az),coset(0,by-1,bz)) +&
                        f3*sc(coset(ax,ay,az),coset(0,by-2,bz)) -&
                  f2 *      kvec(2)*ss(coset(ax,ay,az),coset(0,by-1,bz))
                      ss(coset(ax,ay,az),coset(0,by,bz)) =&
                        rbp(2)*ss(coset(ax,ay,az),coset(0,by-1,bz)) +&
                        fy*ss(coset(ax,ay-1,az),coset(0,by-1,bz)) +&
                        f3*ss(coset(ax,ay,az),coset(0,by-2,bz)) +&
                  f2 *      kvec(2)*sc(coset(ax,ay,az),coset(0,by-1,bz))
                    END DO
                  END IF

!                 *** Increase the angular momentum component x of function b ***

                  IF (ax == 0) THEN
                    DO by=0,lb-1
                      bz = lb - 1 - by
                      sc(coset(ax,ay,az),coset(1,by,bz)) =&
                        rbp(1)*sc(coset(ax,ay,az),coset(0,by,bz)) -&
                   f2 *     kvec(1)*ss(coset(ax,ay,az),coset(0,by,bz))
                      ss(coset(ax,ay,az),coset(1,by,bz)) =&
                        rbp(1)*ss(coset(ax,ay,az),coset(0,by,bz)) +&
                   f2 *     kvec(1)*sc(coset(ax,ay,az),coset(0,by,bz))
                    END DO
                    DO bx=2,lb
                      f3 = f2*REAL(bx-1,dp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        sc(coset(ax,ay,az),coset(bx,by,bz)) =&
                          rbp(1)*sc(coset(ax,ay,az),&
                                   coset(bx-1,by,bz)) +&
                          f3*sc(coset(ax,ay,az),coset(bx-2,by,bz)) -&
                     f2 *     kvec(1)*ss(coset(ax,ay,az),coset(bx-1,by,bz))
                        ss(coset(ax,ay,az),coset(bx,by,bz)) =&
                          rbp(1)*ss(coset(ax,ay,az),&
                                   coset(bx-1,by,bz)) +&
                          f3*ss(coset(ax,ay,az),coset(bx-2,by,bz)) +&
                     f2 *     kvec(1)*sc(coset(ax,ay,az),coset(bx-1,by,bz))
                      END DO
                    END DO
                  ELSE
                    DO by=0,lb-1
                      bz = lb - 1 - by
                      sc(coset(ax,ay,az),coset(1,by,bz)) =&
                        rbp(1)*sc(coset(ax,ay,az),coset(0,by,bz)) +&
                        fx*sc(coset(ax-1,ay,az),coset(0,by,bz)) -&
                     f2 *   kvec(1)*ss(coset(ax,ay,az),coset(0,by,bz))
                      ss(coset(ax,ay,az),coset(1,by,bz)) =&
                        rbp(1)*ss(coset(ax,ay,az),coset(0,by,bz)) +&
                        fx*ss(coset(ax-1,ay,az),coset(0,by,bz)) +&
                     f2 *   kvec(1)*sc(coset(ax,ay,az),coset(0,by,bz))
                    END DO
                    DO bx=2,lb
                      f3 = f2*REAL(bx-1,dp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        sc(coset(ax,ay,az),coset(bx,by,bz)) =&
                          rbp(1)*sc(coset(ax,ay,az),&
                                   coset(bx-1,by,bz)) +&
                          fx*sc(coset(ax-1,ay,az),coset(bx-1,by,bz)) +&
                          f3*sc(coset(ax,ay,az),coset(bx-2,by,bz)) -&
                      f2 *    kvec(1)*ss(coset(ax,ay,az),coset(bx-1,by,bz))
                        ss(coset(ax,ay,az),coset(bx,by,bz)) =&
                          rbp(1)*ss(coset(ax,ay,az),&
                                   coset(bx-1,by,bz)) +&
                          fx*ss(coset(ax-1,ay,az),coset(bx-1,by,bz)) +&
                          f3*ss(coset(ax,ay,az),coset(bx-2,by,bz)) +&
                       f2 *   kvec(1)*sc(coset(ax,ay,az),coset(bx-1,by,bz))
                      END DO
                    END DO
                  END IF

                END DO
              END DO

            END DO

          END IF

        ELSE

          IF (lb_max > 0) THEN

!           *** Vertical recurrence steps: [s|O|s] -> [s|O|b] ***

            rbp(:) = (f1 - 1.0_dp)*rab(:)

!           *** [s|O|p] = (Pi - Bi)*[s|O|s] + [s|dO|s] ***

            sc(1,2) = rbp(1)*sc(1,1) -f2 * kvec(1)*ss(1,1)
            sc(1,3) = rbp(2)*sc(1,1) -f2 * kvec(2)*ss(1,1)
            sc(1,4) = rbp(3)*sc(1,1) -f2 * kvec(3)*ss(1,1)
            ss(1,2) = rbp(1)*ss(1,1) +f2 * kvec(1)*sc(1,1)
            ss(1,3) = rbp(2)*ss(1,1) +f2 * kvec(2)*sc(1,1)
            ss(1,4) = rbp(3)*ss(1,1) +f2 * kvec(3)*sc(1,1)

!           *** [s|O|b] = (Pi - Bi)*[s|O|b-1i] + f2*Ni(b-1i)*[s|O|b-2i] ***
!           ***           + [s|dO|b-1i]                                 ***

            DO lb=2,lb_max

!             *** Increase the angular momentum component z of function b ***

              sc(1,coset(0,0,lb)) = rbp(3)*sc(1,coset(0,0,lb-1)) +&
                                    f2*REAL(lb-1,dp)*sc(1,coset(0,0,lb-2)) -&
                           f2 *         kvec(3)*ss(1,coset(0,0,lb-1))
              ss(1,coset(0,0,lb)) = rbp(3)*ss(1,coset(0,0,lb-1)) +&
                                    f2*REAL(lb-1,dp)*ss(1,coset(0,0,lb-2)) +&
                           f2 *         kvec(3)*sc(1,coset(0,0,lb-1))

!             *** Increase the angular momentum component y of function b ***

              bz = lb - 1
              sc(1,coset(0,1,bz)) = rbp(2)*sc(1,coset(0,0,bz)) -&
                           f2 *         kvec(2)*ss(1,coset(0,0,bz))
              ss(1,coset(0,1,bz)) = rbp(2)*ss(1,coset(0,0,bz)) +&
                           f2 *         kvec(2)*sc(1,coset(0,0,bz))

              DO by=2,lb
                bz = lb - by
                sc(1,coset(0,by,bz)) = rbp(2)*sc(1,coset(0,by-1,bz)) +&
                                       f2*REAL(by-1,dp)*sc(1,coset(0,by-2,bz)) -&
                           f2 *            kvec(2)*ss(1,coset(0,by-1,bz))
                ss(1,coset(0,by,bz)) = rbp(2)*ss(1,coset(0,by-1,bz)) +&
                                       f2*REAL(by-1,dp)*ss(1,coset(0,by-2,bz)) +&
                           f2 *            kvec(2)*sc(1,coset(0,by-1,bz))
              END DO

!             *** Increase the angular momentum component x of function b ***

              DO by=0,lb-1
                bz = lb - 1 - by
                sc(1,coset(1,by,bz)) = rbp(1)*sc(1,coset(0,by,bz)) -&
                              f2 *         kvec(1)*ss(1,coset(0,by,bz))
                ss(1,coset(1,by,bz)) = rbp(1)*ss(1,coset(0,by,bz)) +&
                              f2 *         kvec(1)*sc(1,coset(0,by,bz))
              END DO

              DO bx=2,lb
                f3 = f2*REAL(bx-1,dp)
                DO by=0,lb-bx
                  bz = lb - bx - by
                  sc(1,coset(bx,by,bz)) = rbp(1)*sc(1,coset(bx-1,by,bz)) +&
                                          f3*sc(1,coset(bx-2,by,bz)) -&
                              f2 *            kvec(1)*ss(1,coset(bx-1,by,bz))
                  ss(1,coset(bx,by,bz)) = rbp(1)*ss(1,coset(bx-1,by,bz)) +&
                                          f3*ss(1,coset(bx-2,by,bz)) +&
                               f2 *           kvec(1)*sc(1,coset(bx-1,by,bz))
                END DO
              END DO

            END DO

          END IF

        END IF

        DO j=ncoset(lb_min-1)+1,ncoset(lb_max)
          DO i=ncoset(la_min-1)+1,ncoset(la_max)
            cosab(na+i,nb+j) = sc(i,j)
            sinab(na+i,nb+j) = ss(i,j)
          END DO
        END DO

        nb = nb + ncoset(lb_max)

      END DO

      na = na + ncoset(la_max)

    END DO

  END SUBROUTINE cossin

! *****************************************************************************

  SUBROUTINE moment(la_max,npgfa,zeta,rpgfa,la_min,&
                    lb_max,npgfb,zetb,rpgfb,lb_min,&
                    lc_max,rac,rbc,mab)

    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta, rpgfa
    INTEGER, INTENT(IN)                      :: la_min, lb_max, npgfb, lc_max
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb, rpgfb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rac, rbc
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(OUT)                            :: mab

    INTEGER                                  :: ax, ay, az, bx, by, bz, i, &
                                                ipgf, j, jpgf, la, la_start, &
                                                lb, na, nb, k, l, lx, ly, lz, &
                                                lx1, ly1, lz1, l1, l2, ni
    REAL(KIND=dp)                            :: dab, f0, f1, f2, f3, fx, fy, &
                                                fz, zetp, rab2, f2x, f2y, f2z
    REAL(KIND=dp), DIMENSION(3)              :: rab, rap, rbp, rpc
    REAL(KIND=dp), DIMENSION(ncoset(la_max),&
      ncoset(lb_max),ncoset(lc_max))         :: s

!   ---------------------------------------------------------------------------
    rab = rac - rbc
    rab2= SUM( rab**2 )
    dab = SQRT(rab2)

!   *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0

    DO ipgf=1,npgfa

      nb = 0

      DO jpgf=1,npgfb

!       *** Screening ***

        IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
          DO k=1, ncoset(lc_max)-1
            DO j=nb+1,nb+ncoset(lb_max)
              DO i=na+1,na+ncoset(la_max)
                mab(i,j,k) = 0.0_dp
              END DO
            END DO
          END DO
          nb = nb + ncoset(lb_max)
          CYCLE
        END IF

!       *** Calculate some prefactors ***

        zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))

        f0 = (pi*zetp)**1.5_dp
        f1 = zetb(jpgf)*zetp
        f2 = 0.5_dp*zetp

!       *** Calculate the basic two-center moment integral [s|M|s] ***

        rpc = zetp*(zeta(ipgf)*rac+zetb(jpgf)*rbc)
        s(1,1,1) = f0*EXP(-zeta(ipgf)*f1*rab2)
        DO l=2, ncoset(lc_max)
          lx = indco(1,l)
          ly = indco(2,l)
          lz = indco(3,l)
          l2 = 0
          IF ( lz > 0 ) THEN
            l1 = coset(lx,ly,lz-1)
            IF ( lz > 1 ) l2 = coset(lx,ly,lz-2)
            ni = lz-1
            i  = 3
          ELSE IF ( ly > 0 ) THEN
            l1 = coset(lx,ly-1,lz)
            IF ( ly > 1 ) l2 = coset(lx,ly-2,lz)
            ni = ly-1
            i  = 2
          ELSE IF ( lx > 0 ) THEN
            l1 = coset(lx-1,ly,lz)
            IF ( lx > 1 ) l2 = coset(lx-2,ly,lz)
            ni = lx-1
            i  = 1
          END IF
          s(1,1,l) = rpc(i)*s(1,1,l1)
          IF ( l2 > 0 ) s(1,1,l) = s(1,1,l) + f2*REAL(ni,dp)*s(1,1,l2)
        END DO

!       *** Recurrence steps: [s|M|s] -> [a|M|b] ***

        DO l = 1, ncoset(lc_max)

          lx = indco(1,l)
          ly = indco(2,l)
          lz = indco(3,l)
          IF ( lx > 0 ) THEN 
            lx1 = coset(lx-1,ly,lz)
          ELSE
            lx1 = -1
          END IF
          IF ( ly > 0 ) THEN 
            ly1 = coset(lx,ly-1,lz)
          ELSE
            ly1 = -1
          END IF
          IF ( lz > 0 ) THEN 
            lz1 = coset(lx,ly,lz-1)
          ELSE
            lz1 = -1
          END IF
          f2x = f2*REAL(lx,dp)
          f2y = f2*REAL(ly,dp)
          f2z = f2*REAL(lz,dp)

          IF (la_max > 0) THEN

!           *** Vertical recurrence steps: [s|M|s] -> [a|M|s] ***

            rap(:) = -f1*rab(:)

!           *** [p|M|s] = (Pi - Ai)*[s|M|s] + f2*Ni(m-1i)[s|M-1i|s] ***

            s(2,1,l) = rap(1)*s(1,1,l)
            s(3,1,l) = rap(2)*s(1,1,l)
            s(4,1,l) = rap(3)*s(1,1,l)
            IF ( lx1 < 0 ) s(2,1,l) = s(2,1,l) + f2x*s(1,1,lx1)
            IF ( ly1 < 0 ) s(3,1,l) = s(3,1,l) + f2y*s(1,1,ly1)
            IF ( lz1 < 0 ) s(4,1,l) = s(4,1,l) + f2z*s(1,1,lz1)

!           *** [a|M|s] = (Pi - Ai)*[a-1i|M|s] + f2*Ni(a-1i)*[a-2i|M|s] ***
!           ***           + f2*Ni(m-1i)*[a-1i|M-1i|s]                   ***

            DO la=2,la_max

!             *** Increase the angular momentum component z of function a ***

              s(coset(0,0,la),1,l) = rap(3)*s(coset(0,0,la-1),1,l) +&
                                 f2*REAL(la-1,dp)*s(coset(0,0,la-2),1,l)
              IF ( lz1 > 0 ) s(coset(0,0,la),1,l) = s(coset(0,0,la),1,l) + &
                             f2z*s(coset(0,0,la-1),1,lz1)

!             *** Increase the angular momentum component y of function a ***

              az = la - 1
              s(coset(0,1,az),1,l) = rap(2)*s(coset(0,0,az),1,l)
              IF ( ly1 > 0 ) s(coset(0,1,az),1,l) = s(coset(0,1,az),1,l) + &
                             f2y*s(coset(0,0,az),1,ly1)

              DO ay=2,la
                az = la - ay
                s(coset(0,ay,az),1,l) = rap(2)*s(coset(0,ay-1,az),1,l) +&
                                       f2*REAL(ay-1,dp)*s(coset(0,ay-2,az),1,l)
                IF ( ly1 > 0 ) s(coset(0,ay,az),1,l) = s(coset(0,ay,az),1,l) + &
                             f2y*s(coset(0,ay-1,az),1,ly1)
              END DO

!             *** Increase the angular momentum component x of function a ***

              DO ay=0,la-1
                az = la - 1 - ay
                s(coset(1,ay,az),1,l) = rap(1)*s(coset(0,ay,az),1,l)
                IF ( lx1 > 0 ) s(coset(1,ay,az),1,l) = s(coset(1,ay,az),1,l) + &
                             f2x*s(coset(0,ay,az),1,lx1)
              END DO

              DO ax=2,la
                f3 = f2*REAL(ax-1,dp)
                DO ay=0,la-ax
                  az = la - ax - ay
                  s(coset(ax,ay,az),1,l) = rap(1)*s(coset(ax-1,ay,az),1,l) +&
                                          f3*s(coset(ax-2,ay,az),1,l) 
                  IF ( lx1 > 0 ) s(coset(ax,ay,az),1,l) = s(coset(ax,ay,az),1,l) + &
                                 f2x*s(coset(ax-1,ay,az),1,lx1)
                END DO
              END DO
  
            END DO

!           *** Recurrence steps: [a|M|s] -> [a|M|b] ***

            IF (lb_max > 0) THEN

              DO j=2,ncoset(lb_max)
                DO i=1,ncoset(la_max)
                  s(i,j,l) = 0.0_dp
                END DO
              END DO

!             *** Horizontal recurrence steps ***

              rbp(:) = rap(:) - rab(:)

!             *** [a|M|p] = [a+1i|M|s] - (Bi - Ai)*[a|M|s] ***

              IF (lb_max == 1) THEN
                la_start = la_min
              ELSE
                la_start = MAX(0,la_min-1)
              END IF

              DO la=la_start,la_max-1
                DO ax=0,la
                  DO ay=0,la-ax
                    az = la - ax - ay
                    s(coset(ax,ay,az),2,l) = s(coset(ax+1,ay,az),1,l) -&
                                           rab(1)*s(coset(ax,ay,az),1,l)
                    s(coset(ax,ay,az),3,l) = s(coset(ax,ay+1,az),1,l) -&
                                           rab(2)*s(coset(ax,ay,az),1,l)
                    s(coset(ax,ay,az),4,l) = s(coset(ax,ay,az+1),1,l) -&
                                           rab(3)*s(coset(ax,ay,az),1,l)
                  END DO
                END DO
              END DO

!             *** Vertical recurrence step ***

!             *** [a|M|p] = (Pi - Bi)*[a|M|s] + f2*Ni(a)*[a-1i|M|s] ***
!             ***           + f2*Ni(m)*[a|M-1i|s]                   ***

              DO ax=0,la_max
                fx = f2*REAL(ax,dp)
                DO ay=0,la_max-ax
                  fy = f2*REAL(ay,dp)
                  az = la_max - ax - ay
                  fz = f2*REAL(az,dp)
                  IF (ax == 0) THEN
                    s(coset(ax,ay,az),2,l) = rbp(1)*s(coset(ax,ay,az),1,l)
                  ELSE
                    s(coset(ax,ay,az),2,l) = rbp(1)*s(coset(ax,ay,az),1,l) +&
                                            fx*s(coset(ax-1,ay,az),1,l)
                  END IF
                  IF (lx1 > 0) s(coset(ax,ay,az),2,l) = s(coset(ax,ay,az),2,l) +&
                        f2x*s(coset(ax,ay,az),1,lx1)
                  IF (ay == 0) THEN
                    s(coset(ax,ay,az),3,l) = rbp(2)*s(coset(ax,ay,az),1,l)
                  ELSE
                    s(coset(ax,ay,az),3,l) = rbp(2)*s(coset(ax,ay,az),1,l) +&
                                            fy*s(coset(ax,ay-1,az),1,l)
                  END IF
                  IF (ly1 > 0) s(coset(ax,ay,az),3,l) = s(coset(ax,ay,az),3,l) +&
                        f2y*s(coset(ax,ay,az),1,ly1)
                  IF (az == 0) THEN
                    s(coset(ax,ay,az),4,l) = rbp(3)*s(coset(ax,ay,az),1,l)
                  ELSE
                    s(coset(ax,ay,az),4,l) = rbp(3)*s(coset(ax,ay,az),1,l) +&
                                            fz*s(coset(ax,ay,az-1),1,l)
                  END IF
                  IF (lz1 > 0) s(coset(ax,ay,az),4,l) = s(coset(ax,ay,az),4,l) +&
                        f2z*s(coset(ax,ay,az),1,lz1)
                END DO
              END DO

!             *** Recurrence steps: [a|M|p] -> [a|M|b] ***

              DO lb=2,lb_max

!               *** Horizontal recurrence steps ***

!               *** [a|M|b] = [a+1i|M|b-1i] - (Bi - Ai)*[a|M|b-1i] ***

                IF (lb == lb_max) THEN
                  la_start = la_min
                ELSE
                  la_start = MAX(0,la_min-1)
                END IF

                DO la=la_start,la_max-1
                  DO ax=0,la
                    DO ay=0,la-ax
                      az = la - ax - ay

!                     *** Shift of angular momentum component z from a to b ***

                      s(coset(ax,ay,az),coset(0,0,lb),l) =&
                        s(coset(ax,ay,az+1),coset(0,0,lb-1),l) -&
                        rab(3)*s(coset(ax,ay,az),coset(0,0,lb-1),l)

!                     *** Shift of angular momentum component y from a to b ***

                      DO by=1,lb
                        bz = lb - by
                        s(coset(ax,ay,az),coset(0,by,bz),l) =&
                          s(coset(ax,ay+1,az),coset(0,by-1,bz),l) -&
                          rab(2)*s(coset(ax,ay,az),coset(0,by-1,bz),l)
                      END DO

!                     *** Shift of angular momentum component x from a to b ***

                      DO bx=1,lb
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          s(coset(ax,ay,az),coset(bx,by,bz),l) =&
                            s(coset(ax+1,ay,az),coset(bx-1,by,bz),l) -&
                            rab(1)*s(coset(ax,ay,az),coset(bx-1,by,bz),l)
                        END DO
                      END DO

                    END DO
                  END DO
                END DO

!               *** Vertical recurrence step ***

!               *** [a|M|b] = (Pi - Bi)*[a|M|b-1i] + f2*Ni(a)*[a-1i|M|b-1i] + ***
!               ***           f2*Ni(b-1i)*[a|M|b-2i] + f2*Ni(m)[a|M-1i|b-1i]  ***

                DO ax=0,la_max
                  fx = f2*REAL(ax,dp)
                  DO ay=0,la_max-ax
                    fy = f2*REAL(ay,dp)
                    az = la_max - ax - ay
                    fz = f2*REAL(az,dp)

!                   *** Shift of angular momentum component z from a to b ***

                    f3 = f2*REAL(lb-1,dp)

                    IF (az == 0) THEN
                      s(coset(ax,ay,az),coset(0,0,lb),l) =&
                        rbp(3)*s(coset(ax,ay,az),coset(0,0,lb-1),l) +&
                        f3*s(coset(ax,ay,az),coset(0,0,lb-2),l)
                    ELSE
                      s(coset(ax,ay,az),coset(0,0,lb),l) =&
                        rbp(3)*s(coset(ax,ay,az),coset(0,0,lb-1),l) +&
                        fz*s(coset(ax,ay,az-1),coset(0,0,lb-1),l) +&
                        f3*s(coset(ax,ay,az),coset(0,0,lb-2),l)
                    END IF
                    IF (lz1 > 0) s(coset(ax,ay,az),coset(0,0,lb),l) =&
                                 s(coset(ax,ay,az),coset(0,0,lb),l) +&
                                 f2z*s(coset(ax,ay,az),coset(0,0,lb-1),lz1)

!                   *** Shift of angular momentum component y from a to b ***

                    IF (ay == 0) THEN
                      bz = lb - 1
                      s(coset(ax,ay,az),coset(0,1,bz),l) =&
                        rbp(2)*s(coset(ax,ay,az),coset(0,0,bz),l)
                      IF (ly1 > 0) s(coset(ax,ay,az),coset(0,1,bz),l) =&
                                 s(coset(ax,ay,az),coset(0,1,bz),l) +&
                                 f2y*s(coset(ax,ay,az),coset(0,0,bz),ly1)
                      DO by=2,lb
                        bz = lb - by
                        f3 = f2*REAL(by-1,dp)
                        s(coset(ax,ay,az),coset(0,by,bz),l) =&
                          rbp(2)*s(coset(ax,ay,az),coset(0,by-1,bz),l) +&
                          f3*s(coset(ax,ay,az),coset(0,by-2,bz),l)
                        IF (ly1 > 0) s(coset(ax,ay,az),coset(0,by,bz),l) =&
                                 s(coset(ax,ay,az),coset(0,by,bz),l) +&
                                 f2y*s(coset(ax,ay,az),coset(0,by-1,bz),ly1)
                      END DO
                    ELSE
                      bz = lb - 1
                      s(coset(ax,ay,az),coset(0,1,bz),l) =&
                        rbp(2)*s(coset(ax,ay,az),coset(0,0,bz),l) +&
                        fy*s(coset(ax,ay-1,az),coset(0,0,bz),l)
                      IF (ly1 > 0) s(coset(ax,ay,az),coset(0,1,bz),l) =&
                                 s(coset(ax,ay,az),coset(0,1,bz),l) +&
                                 f2y*s(coset(ax,ay,az),coset(0,0,bz),ly1)
                      DO by=2,lb
                        bz = lb - by
                        f3 = f2*REAL(by-1,dp)
                        s(coset(ax,ay,az),coset(0,by,bz),l) =&
                          rbp(2)*s(coset(ax,ay,az),coset(0,by-1,bz),l) +&
                          fy*s(coset(ax,ay-1,az),coset(0,by-1,bz),l) +&
                          f3*s(coset(ax,ay,az),coset(0,by-2,bz),l)
                        IF (ly1 > 0) s(coset(ax,ay,az),coset(0,by,bz),l) =&
                                 s(coset(ax,ay,az),coset(0,by,bz),l) +&
                                 f2y*s(coset(ax,ay,az),coset(0,by-1,bz),ly1)
                      END DO
                    END IF

!                   *** Shift of angular momentum component x from a to b ***

                    IF (ax == 0) THEN
                      DO by=0,lb-1
                        bz = lb - 1 - by
                        s(coset(ax,ay,az),coset(1,by,bz),l) =&
                          rbp(1)*s(coset(ax,ay,az),coset(0,by,bz),l)
                        IF (lx1 > 0) s(coset(ax,ay,az),coset(1,by,bz),l) =&
                                 s(coset(ax,ay,az),coset(1,by,bz),l) +&
                                 f2x*s(coset(ax,ay,az),coset(0,by,bz),lx1)
                      END DO
                      DO bx=2,lb
                        f3 = f2*REAL(bx-1,dp)
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          s(coset(ax,ay,az),coset(bx,by,bz),l) =&
                            rbp(1)*s(coset(ax,ay,az),coset(bx-1,by,bz),l) +&
                            f3*s(coset(ax,ay,az),coset(bx-2,by,bz),l)
                          IF (lx1 > 0) s(coset(ax,ay,az),coset(bx,by,bz),l) =&
                                 s(coset(ax,ay,az),coset(bx,by,bz),l) +&
                                 f2x*s(coset(ax,ay,az),coset(bx-1,by,bz),lx1)
                        END DO
                      END DO
                    ELSE
                      DO by=0,lb-1
                        bz = lb - 1 - by
                        s(coset(ax,ay,az),coset(1,by,bz),l) =&
                          rbp(1)*s(coset(ax,ay,az),coset(0,by,bz),l) +&
                          fx*s(coset(ax-1,ay,az),coset(0,by,bz),l)
                        IF (lx1 > 0) s(coset(ax,ay,az),coset(1,by,bz),l) =&
                                 s(coset(ax,ay,az),coset(1,by,bz),l) +&
                                 f2x*s(coset(ax,ay,az),coset(0,by,bz),lx1)
                      END DO
                      DO bx=2,lb
                        f3 = f2*REAL(bx-1,dp)
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          s(coset(ax,ay,az),coset(bx,by,bz),l) =&
                            rbp(1)*s(coset(ax,ay,az),coset(bx-1,by,bz),l) +&
                            fx*s(coset(ax-1,ay,az),coset(bx-1,by,bz),l) +&
                            f3*s(coset(ax,ay,az),coset(bx-2,by,bz),l)
                          IF (lx1 > 0) s(coset(ax,ay,az),coset(bx,by,bz),l) =&
                                 s(coset(ax,ay,az),coset(bx,by,bz),l) +&
                                 f2x*s(coset(ax,ay,az),coset(bx-1,by,bz),lx1)
                        END DO
                      END DO
                    END IF

                  END DO
                END DO

              END DO

            END IF

          ELSE

            IF (lb_max > 0) THEN

!             *** Vertical recurrence steps: [s|M|s] -> [s|M|b] ***

              rbp(:) = (f1 - 1.0_dp)*rab(:)

!             *** [s|M|p] = (Pi - Bi)*[s|M|s] + f2*Ni(m)*[s|M-1i|s] ***

              s(1,2,l) = rbp(1)*s(1,1,l)
              s(1,3,l) = rbp(2)*s(1,1,l)
              s(1,4,l) = rbp(3)*s(1,1,l)
              IF ( lx1 < 0 ) s(1,2,l) = s(1,2,l) + f2x*s(1,1,lx1)
              IF ( ly1 < 0 ) s(1,3,l) = s(1,3,l) + f2y*s(1,1,ly1)
              IF ( lz1 < 0 ) s(1,4,l) = s(1,4,l) + f2z*s(1,1,lz1)

!             *** [s|M|b] = (Pi - Bi)*[s|M|b-1i] + f2*Ni(b-1i)*[s|M|b-2i] ***
!             ***           + f2*Ni(m)*[s|M-1i|b-1i]                      ***

              DO lb=2,lb_max

!               *** Increase the angular momentum component z of function b ***

                s(1,coset(0,0,lb),l) = rbp(3)*s(1,coset(0,0,lb-1),l) +&
                                      f2*REAL(lb-1,dp)*s(1,coset(0,0,lb-2),l)
                IF (lz1 > 0) s(1,coset(0,0,lb),l) = s(1,coset(0,0,lb),l) +&
                             f2z*s(1,coset(0,0,lb-1),lz1)

!               *** Increase the angular momentum component y of function b ***

                bz = lb - 1
                s(1,coset(0,1,bz),l) = rbp(2)*s(1,coset(0,0,bz),l)
                IF (ly1 > 0) s(1,coset(0,1,bz),l) = s(1,coset(0,1,bz),l) +&
                             f2y*s(1,coset(0,0,bz),ly1)

              DO by=2,lb
                bz = lb - by
                s(1,coset(0,by,bz),l) = rbp(2)*s(1,coset(0,by-1,bz),l) +&
                                       f2*REAL(by-1,dp)*s(1,coset(0,by-2,bz),l)
                IF (ly1 > 0) s(1,coset(0,by,bz),l) = s(1,coset(0,by,bz),l) +&
                             f2y*s(1,coset(0,by-1,bz),ly1)
              END DO

!             *** Increase the angular momentum component x of function b ***

              DO by=0,lb-1
                bz = lb - 1 - by
                s(1,coset(1,by,bz),l) = rbp(1)*s(1,coset(0,by,bz),l)
                IF (lx1 > 0) s(1,coset(1,by,bz),l) = s(1,coset(1,by,bz),l) +&
                             f2x*s(1,coset(0,by,bz),lx1)
              END DO

              DO bx=2,lb
                f3 = f2*REAL(bx-1,dp)
                DO by=0,lb-bx
                  bz = lb - bx - by
                  s(1,coset(bx,by,bz),l) = rbp(1)*s(1,coset(bx-1,by,bz),l) +&
                                          f3*s(1,coset(bx-2,by,bz),l)
                  IF (lx1 > 0) s(1,coset(bx,by,bz),l) = s(bx,coset(1,by,bz),l) +&
                               f2x*s(1,coset(bx-1,by,bz),lx1)
                END DO
              END DO

            END DO

          END IF

        END IF

        END DO

        DO k=2,ncoset(lc_max)
          DO j=1,ncoset(lb_max)
            DO i=1,ncoset(la_max)
              mab(na+i,nb+j,k-1) = s(i,j,k)
            END DO
          END DO
        END DO

        nb = nb + ncoset(lb_max)

      END DO

      na = na + ncoset(la_max)

    END DO

  END SUBROUTINE moment

! *****************************************************************************

  SUBROUTINE diffop(la_max,npgfa,zeta,rpgfa,la_min,&
                    lb_max,npgfb,zetb,rpgfb,lb_min,&
                    rab,difab)

    INTEGER, INTENT(IN)                      :: la_max,npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta, rpgfa
    INTEGER, INTENT(IN)                      :: la_min,npgfb,lb_max
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb, rpgfb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(OUT)                            :: difab

    INTEGER                                  :: istat, lda_min, ldb_min, lmax, lds
    REAL(KIND=dp)                            :: dab, rab2
    REAL(KIND=dp), DIMENSION(npgfa*ncoset(la_max+1),&
      npgfb*ncoset(lb_max+1))           :: sab, sabdx, sabdy, sabdz
    REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: s

!   ---------------------------------------------------------------------------
    rab2= SUM( rab**2 )
    dab = SQRT(rab2)

    lda_min = MAX ( 0, la_min-1 )
    ldb_min = MAX ( 0, lb_min-1 )
    lmax = MAX(la_max+1,lb_max+1)
    lds = ncoset(lmax+1)
    ALLOCATE(s(lds,lds,ncoset(1)),STAT=istat)
    
    CALL overlap(la_max+1,lda_min,npgfa,rpgfa,zeta,&
                 lb_max+1,ldb_min,npgfb,rpgfb,zetb,&
                 rab,dab,sab,0,.FALSE.,s,lds)

    CALL dabdr(la_max,npgfa,zeta,rpgfa,la_min,&
               lb_max,npgfb,zetb,rpgfb,lb_min,&
               dab,sab,difab(:,:,1),difab(:,:,2),difab(:,:,3))

    DEALLOCATE(s,STAT=istat)

  END SUBROUTINE diffop

! *****************************************************************************

END MODULE ai_moments
