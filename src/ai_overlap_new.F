!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/ai_overlap [1.0] *
!!
!!   NAME
!!     ai_overlap
!!
!!   FUNCTION
!!     Calculation of the overlap integrals over Cartesian Gaussian-type
!!     functions.
!!
!!   AUTHOR
!!     Matthias Krack (08.10.1999)
!!
!!   LITERATURE
!!     S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!!
!!   MODIFICATION HISTORY
!!     - Derivatives added (02.05.2002,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE ai_overlap_new

! *****************************************************************************

! ax,ay,az  : Angular momentum index numbers of orbital a.
! bx,by,bz  : Angular momentum index numbers of orbital b.
! coset     : Cartesian orbital set pointer.
! dab       : Distance between the atomic centers a and b.
! l{a,b}    : Angular momentum quantum number of shell a or b.
! l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
! l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
! rab       : Distance vector between the atomic centers a and b.
! rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
! sab       : Shell set of overlap integrals.
! zet{a,b}  : Exponents of the Gaussian-type functions a or b.
! zetp      : Reciprocal of the sum of the exponents of orbital a and b.

! *****************************************************************************

  USE kinds,                           ONLY: dp,dp_size
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,ncoset
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "ai_overlap_new"

! *** Public subroutines ***

  PUBLIC :: overlap

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE overlap(la_max_set,zeta,rpgfa,la_min_set,&
                     lb_max_set,zetb,rpgfb,lb_min_set,&
                     rab,dab,sab,da_max_set,return_derivatives,&
                     sdab,pab,force_a)

!   Purpose: Calculation of the two-center overlap integrals [a|b] over
!            Cartesian Gaussian-type functions.

!   History: - Creation (19.09.2000,MK)

!   ***************************************************************************

    REAL(KIND = dp), INTENT(IN)                       :: dab
    INTEGER, INTENT(IN)                               :: da_max_set,&
                                                         la_max_set,&
                                                         la_min_set,&
                                                         lb_max_set,&
                                                         lb_min_set
    REAL(KIND = dp), DIMENSION(3), INTENT(IN)         :: rab
    REAL(KIND = dp), DIMENSION(:), INTENT(IN)         :: rpgfa,rpgfb,zeta,zetb
    REAL(KIND = dp), DIMENSION(:,:), INTENT(INOUT)    :: sab
    LOGICAL, INTENT(IN)                               :: return_derivatives

    REAL(KIND = dp), DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: sdab
    REAL(KIND = dp), DIMENSION(:,:), OPTIONAL, INTENT(IN)      :: pab
    REAL(KIND = dp), DIMENSION(3), OPTIONAL, INTENT(OUT)       :: force_a

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "overlap"

!   *** Local variables ***

    REAL(KIND = dp) :: f0,f1,f2,f3,fax,fay,faz,ftz,zetp
    INTEGER         :: ax,ay,az,bx,by,bz,cda,cdax,cday,cdaz,coa,&
                       coamx,coamy,coamz,coapx,coapy,coapz,cob,&
                       cobm2x,cobm2y,cobm2z,cobmx,cobmy,cobmz,&
                       da,da_max,dax,day,daz,i,ipgf,istat,&
                       j,jk,jpgf,jstart,k,la,la_max,la_min,la_start,&
                       lb,lb_max,lb_min,na,nb,nda
    LOGICAL         :: calculate_force_a

    REAL(KIND = dp), DIMENSION(3) :: rap,rbp

    REAL(KIND = dp), DIMENSION(:,:,:), ALLOCATABLE :: s

!   ---------------------------------------------------------------------------

    IF (PRESENT(pab).AND.PRESENT(force_a)) THEN
      calculate_force_a = .TRUE.
      force_a(:) = 0.0_dp
    ELSE
      calculate_force_a = .FALSE.
    END IF

    IF (PRESENT(sdab).OR.calculate_force_a) THEN
      IF (da_max_set == 0) THEN
        da_max = 1
        la_max = la_max_set + 1
        la_min = MAX(0,la_min_set-1)
      ELSE
        da_max = da_max_set
        la_max = la_max_set + da_max_set + 1
        la_min = MAX(0,la_min_set-da_max_set-1)
      END IF
    ELSE
      da_max = da_max_set
      la_max = la_max_set + da_max_set
      la_min = MAX(0,la_min_set-da_max_set)
    END IF

    lb_max = lb_max_set
    lb_min = lb_min_set

    ALLOCATE (s(ncoset(la_max),ncoset(lb_max),ncoset(da_max)),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "s",ncoset(la_max)*&
                                         ncoset(lb_max)*&
                                         ncoset(da_max)*dp_size)

!   *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0
    nda = 0

    DO ipgf=1,SIZE(zeta)

      nb = 0

      DO jpgf=1,SIZE(zetb)

!       *** Screening ***

        IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
          DO k=1,ncoset(da_max_set)
            jstart = (k - 1)*SIZE(sab,1)
            DO j=jstart+nb+1,jstart+nb+ncoset(lb_max_set)
              DO i=na+1,na+ncoset(la_max_set)
                sab(i,j) = 0.0_dp
              END DO
            END DO
          END DO
          nb = nb + ncoset(lb_max_set)
          CYCLE
        END IF

!       *** Calculate some prefactors ***

        zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))

        f0 = (pi*zetp)**1.5_dp
        f1 = zetb(jpgf)*zetp
        f2 = 0.5_dp*zetp

!       *** Calculate the basic two-center overlap integral [s|s] ***

        s(1,1,1) = f0*EXP(-zeta(ipgf)*f1*dab*dab)

!       *** Recurrence steps: [s|s] -> [a|b] ***

        IF (la_max > 0) THEN

!         *** Vertical recurrence steps: [s|s] -> [a|s] ***

          rap(:) = f1*rab(:)

!         *** [p|s] = (Pi - Ai)*[s|s]  (i = x,y,z) ***

          s(2,1,1) = rap(1)*s(1,1,1)
          s(3,1,1) = rap(2)*s(1,1,1)
          s(4,1,1) = rap(3)*s(1,1,1)

!         *** [a|s] = (Pi - Ai)*[a-1i|s] + f2*Ni(a-1i)*[a-2i|s] ***

          DO la=2,la_max

!           *** Increase the angular momentum component z of a ***

            s(coset(0,0,la),1,1) = rap(3)*s(coset(0,0,la-1),1,1) +&
                                   f2*REAL(la-1,dp)*s(coset(0,0,la-2),1,1)

!           *** Increase the angular momentum component y of a ***

            az = la - 1
            s(coset(0,1,az),1,1) = rap(2)*s(coset(0,0,az),1,1)
            DO ay=2,la
              az = la - ay
              s(coset(0,ay,az),1,1) = rap(2)*s(coset(0,ay-1,az),1,1) +&
                                      f2*REAL(ay-1,dp)*s(coset(0,ay-2,az),1,1)
            END DO

!           *** Increase the angular momentum component x of a ***

            DO ay=0,la-1
              az = la - 1 - ay
              s(coset(1,ay,az),1,1) = rap(1)*s(coset(0,ay,az),1,1)
            END DO
            DO ax=2,la
              f3 = f2*REAL(ax-1,dp)
              DO ay=0,la-ax
                az = la - ax - ay
                s(coset(ax,ay,az),1,1) = rap(1)*s(coset(ax-1,ay,az),1,1) +&
                                         f3*s(coset(ax-2,ay,az),1,1)
              END DO
            END DO

          END DO

!         *** Recurrence steps: [a|s] -> [a|b] ***

          IF (lb_max > 0) THEN

            DO j=2,ncoset(lb_max)
              DO i=1,ncoset(la_min)
                s(i,j,1) = 0.0_dp
              END DO
            END DO

!           *** Horizontal recurrence steps ***

            rbp(:) = rap(:) - rab(:)

!           *** [a|p] = [a+1i|s] - (Bi - Ai)*[a|s] ***

            IF (lb_max == 1) THEN
              la_start = la_min
            ELSE
              la_start = MAX(0,la_min-1)
            END IF

            DO la=la_start,la_max-1
              DO ax=0,la
                DO ay=0,la-ax
                  az = la - ax - ay
                  coa = coset(ax,ay,az)
                  coapx = coset(ax+1,ay,az)
                  coapy = coset(ax,ay+1,az)
                  coapz = coset(ax,ay,az+1)
                  s(coa,2,1) = s(coapx,1,1) - rab(1)*s(coa,1,1)
                  s(coa,3,1) = s(coapy,1,1) - rab(2)*s(coa,1,1)
                  s(coa,4,1) = s(coapz,1,1) - rab(3)*s(coa,1,1)
                END DO
              END DO
            END DO

!           *** Vertical recurrence step ***

!           *** [a|p] = (Pi - Bi)*[a|s] + f2*Ni(a)*[a-1i|s] ***

            DO ax=0,la_max
              fax = f2*REAL(ax,dp)
              DO ay=0,la_max-ax
                fay = f2*REAL(ay,dp)
                az = la_max - ax - ay
                faz = f2*REAL(az,dp)
                coa = coset(ax,ay,az)
                coamx = coset(ax-1,ay,az)
                coamy = coset(ax,ay-1,az)
                coamz = coset(ax,ay,az-1)
                s(coa,2,1) = rbp(1)*s(coa,1,1) + fax*s(coamx,1,1)
                s(coa,3,1) = rbp(2)*s(coa,1,1) + fay*s(coamy,1,1)
                s(coa,4,1) = rbp(3)*s(coa,1,1) + faz*s(coamz,1,1)
              END DO
            END DO

!           *** Recurrence steps: [a|p] -> [a|b] ***

            DO lb=2,lb_max

!             *** Horizontal recurrence steps ***

!             *** [a|b] = [a+1i|b-1i] - (Bi - Ai)*[a|b-1i] ***

              IF (lb == lb_max) THEN
                la_start = la_min
              ELSE
                la_start = MAX(0,la_min-1)
              END IF

              DO la=la_start,la_max-1
                DO ax=0,la
                  DO ay=0,la-ax
                    az = la - ax - ay
                    coa = coset(ax,ay,az)
                    coapx = coset(ax+1,ay,az)
                    coapy = coset(ax,ay+1,az)
                    coapz = coset(ax,ay,az+1)

!                   *** Shift of angular momentum component z from a to b ***

                    cob = coset(0,0,lb)
                    cobmz = coset(0,0,lb-1)
                    s(coa,cob,1) = s(coapz,cobmz,1) - rab(3)*s(coa,cobmz,1)

!                   *** Shift of angular momentum component y from a to b ***

                    DO by=1,lb
                      bz = lb - by
                      cob = coset(0,by,bz)
                      cobmy = coset(0,by-1,bz)
                      s(coa,cob,1) = s(coapy,cobmy,1) - rab(2)*s(coa,cobmy,1)
                    END DO

!                   *** Shift of angular momentum component x from a to b ***

                    DO bx=1,lb
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        cob = coset(bx,by,bz)
                        cobmx = coset(bx-1,by,bz)
                        s(coa,cob,1) = s(coapx,cobmx,1) - rab(1)*s(coa,cobmx,1)
                      END DO
                    END DO

                  END DO
                END DO
              END DO

!             *** Vertical recurrence step ***

!             *** [a|b] = (Pi - Bi)*[a|b-1i] + f2*Ni(a)*[a-1i|b-1i] + ***
!             ***         f2*Ni(b-1i)*[a|b-2i]                        ***

              DO ax=0,la_max
                fax = f2*REAL(ax,dp)
                DO ay=0,la_max-ax
                  fay = f2*REAL(ay,dp)
                  az = la_max - ax - ay
                  faz = f2*REAL(az,dp)
                  coa = coset(ax,ay,az)
                  coamx = coset(ax-1,ay,az)
                  coamy = coset(ax,ay-1,az)
                  coamz = coset(ax,ay,az-1)

!                 *** Increase the angular momentum component z of b ***

                  f3 = f2*REAL(lb-1,dp)
                  cob = coset(0,0,lb)
                  cobmz = coset(0,0,lb-1)
                  cobm2z = coset(0,0,lb-2)
                  s(coa,cob,1) = rbp(3)*s(coa,cobmz,1) +&
                                 faz*s(coamz,cobmz,1) +&
                                 f3*s(coa,cobm2z,1)

!                 *** Increase the angular momentum component y of b ***

                  bz = lb - 1
                  cob = coset(0,1,bz)
                  cobmy = coset(0,0,bz)
                  s(coa,cob,1) = rbp(2)*s(coa,cobmy,1) +&
                                 fay*s(coamy,cobmy,1)
                  DO by=2,lb
                    bz = lb - by
                    f3 = f2*REAL(by-1,dp)
                    cob = coset(0,by,bz)
                    cobmy = coset(0,by-1,bz)
                    cobm2y = coset(0,by-2,bz)
                    s(coa,cob,1) = rbp(2)*s(coa,cobmy,1) +&
                                   fay*s(coamy,cobmy,1) +&
                                   f3*s(coa,cobm2y,1)
                  END DO

!                 *** Increase the angular momentum component x of b ***

                  DO by=0,lb-1
                    bz = lb - 1 - by
                    cob = coset(1,by,bz)
                    cobmx = coset(0,by,bz)
                    s(coa,cob,1) = rbp(1)*s(coa,cobmx,1) +&
                                   fax*s(coamx,cobmx,1)
                  END DO
                  DO bx=2,lb
                    f3 = f2*REAL(bx-1,dp)
                    DO by=0,lb-bx
                      bz = lb - bx - by
                      cob = coset(bx,by,bz)
                      cobmx = coset(bx-1,by,bz)
                      cobm2x = coset(bx-2,by,bz)
                      s(coa,cob,1) = rbp(1)*s(coa,cobmx,1) +&
                                     fax*s(coamx,cobmx,1) +&
                                     f3*s(coa,cobm2x,1)
                    END DO
                  END DO

                END DO
              END DO

            END DO

          END IF

        ELSE

          IF (lb_max > 0) THEN

!           *** Vertical recurrence steps: [s|s] -> [s|b] ***

            rbp(:) = (f1 - 1.0_dp)*rab(:)

!           *** [s|p] = (Pi - Bi)*[s|s] ***

            s(1,2,1) = rbp(1)*s(1,1,1)
            s(1,3,1) = rbp(2)*s(1,1,1)
            s(1,4,1) = rbp(3)*s(1,1,1)

!           *** [s|b] = (Pi - Bi)*[s|b-1i] + f2*Ni(b-1i)*[s|b-2i] ***

            DO lb=2,lb_max

!             *** Increase the angular momentum component z of b ***

              s(1,coset(0,0,lb),1) = rbp(3)*s(1,coset(0,0,lb-1),1) +&
                                     f2*REAL(lb-1,dp)*s(1,coset(0,0,lb-2),1)

!             *** Increase the angular momentum component y of b ***

              bz = lb - 1
              s(1,coset(0,1,bz),1) = rbp(2)*s(1,coset(0,0,bz),1)
              DO by=2,lb
                bz = lb - by
                s(1,coset(0,by,bz),1) = rbp(2)*s(1,coset(0,by-1,bz),1) +&
                                        f2*REAL(by-1,dp)*s(1,coset(0,by-2,bz),1)
              END DO

!             *** Increase the angular momentum component x of b ***

              DO by=0,lb-1
                bz = lb - 1 - by
                s(1,coset(1,by,bz),1) = rbp(1)*s(1,coset(0,by,bz),1)
              END DO
              DO bx=2,lb
                f3 = f2*REAL(bx-1,dp)
                DO by=0,lb-bx
                  bz = lb - bx - by
                  s(1,coset(bx,by,bz),1) = rbp(1)*s(1,coset(bx-1,by,bz),1) +&
                                           f3*s(1,coset(bx-2,by,bz),1)
                END DO
              END DO

            END DO

          END IF

        END IF

!       *** Store the primitive overlap integrals ***

        DO j=ncoset(lb_min_set-1)+1,ncoset(lb_max_set)
          DO i=ncoset(la_min_set-1)+1,ncoset(la_max_set)
            sab(na+i,nb+j) = s(i,j,1)
          END DO
        END DO

!       *** Calculate the requested derivatives with respect  ***
!       *** to the nuclear coordinates of the atomic center a ***

        DO da=0,da_max-1
          ftz = 2.0_dp*zeta(ipgf)
          DO dax=0,da
            DO day=0,da-dax
              daz = da - dax - day
              cda = coset(dax,day,daz)
              cdax = coset(dax+1,day,daz)
              cday = coset(dax,day+1,daz)
              cdaz = coset(dax,day,daz+1)

!             *** [da/dAi|b] = 2*zeta*[a+1i|b] - Ni(a)[a-1i|b] ***

              DO la=0,la_max-da-1
                DO ax=0,la
                  fax = REAL(ax,dp)
                  DO ay=0,la-ax
                    fay = REAL(ay,dp)
                    az = la - ax - ay
                    faz = REAL(az,dp)
                    coa = coset(ax,ay,az)
                    coamx = coset(ax-1,ay,az)
                    coamy = coset(ax,ay-1,az)
                    coamz = coset(ax,ay,az-1)
                    coapx = coset(ax+1,ay,az)
                    coapy = coset(ax,ay+1,az)
                    coapz = coset(ax,ay,az+1)
                    DO lb=0,lb_max_set
                      DO bx=0,lb
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          cob = coset(bx,by,bz)
                          s(coa,cob,cdax) = ftz*s(coapx,cob,cda) -&
                                            fax*s(coamx,cob,cda)
                          s(coa,cob,cday) = ftz*s(coapy,cob,cda) -&
                                            fay*s(coamy,cob,cda)
                          s(coa,cob,cdaz) = ftz*s(coapz,cob,cda) -&
                                            faz*s(coamz,cob,cda)
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO

            END DO
          END DO
        END DO

!       *** Return all the calculated derivatives of the ***
!       *** primitive overlap integrals, if requested    ***

        IF (return_derivatives) THEN
          DO k=2,ncoset(da_max_set)
            jstart = (k - 1)*SIZE(sab,1)
            DO j=ncoset(lb_min_set-1)+1,ncoset(lb_max_set)
              jk = jstart + j
              DO i=ncoset(la_min_set-1)+1,ncoset(la_max_set)
                sab(na+i,nb+jk) = s(i,j,k)
              END DO
            END DO
          END DO
        END IF

!       *** Calculate the force contribution for the atomic center a ***

        IF (calculate_force_a) THEN
          DO k=1,3
            DO j=ncoset(lb_min_set-1)+1,ncoset(lb_max_set)
              DO i=ncoset(la_min_set-1)+1,ncoset(la_max_set)
                force_a(k) = force_a(k) + pab(na+i,nb+j)*s(i,j,k+1)
              END DO
            END DO
          END DO
        END IF

!       *** Store the first derivatives of the primitive overlap integrals ***
!       *** which are used as auxiliary integrals for the calculation of   ***
!       *** the kinetic energy integrals if requested                      ***

        IF (PRESENT(sdab)) THEN
          sdab(nda+1,nb+1,1) = s(1,1,1)
          DO k=2,4
            DO j=1,ncoset(lb_max_set)
              DO i=1,ncoset(la_max-1)
                sdab(nda+i,nb+j,k) = s(i,j,k)
              END DO
            END DO
          END DO
        END IF

        nb = nb + ncoset(lb_max_set)

      END DO

      na = na + ncoset(la_max_set)
      nda = nda + ncoset(la_max-1)

    END DO

    DEALLOCATE (s,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"s")

  END SUBROUTINE overlap

! *****************************************************************************

END MODULE ai_overlap_new
