
MODULE ai_derivatives

! Purpose: Calculate the first derivative of an integral block.

! Literature: S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)

! History: - Creation (05.10.2000, Matthias Krack)

! *****************************************************************************

! ax,ay,az  : Angular momentum index numbers of orbital a.
! bx,by,bz  : Angular momentum index numbers of orbital b.
! coset     : Cartesian orbital set pointer.
! l{a,b}    : Angular momentum quantum number of shell a or b.
! l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
! l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
! ncoset    : Number of orbitals in a Cartesian orbital set.
! npgf{a,b} : Degree of contraction of shell a or b.
! rab       : Distance vector between the atomic centers a and b.
! rab2      : Square of the distance between the atomic centers a and b.
! rac       : Distance vector between the atomic centers a and c.
! rac2      : Square of the distance between the atomic centers a and c.
! rbc       : Distance vector between the atomic centers b and c.
! rbc2      : Square of the distance between the atomic centers b and c.
! rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
! zet{a,b}  : Exponents of the Gaussian-type functions a or b.
! zetp      : Reciprocal of the sum of the exponents of orbital a and b.

! *****************************************************************************

  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: dabdr

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE dabdr(la_max,npgfa,zeta,rpgfa,la_min,&
                   lb_max,npgfb,zetb,rpgfb,lb_min,&
                   dab,ab,dabdx,dabdy,dabdz)

!   Purpose: Calculate the first derivative of an integral block.

!   History: - Creation (05.10.2000, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta, rpgfa
    INTEGER, INTENT(IN)                      :: la_min, npgfb,lb_max
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb, rpgfb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), INTENT(IN)                :: dab
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: ab
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: dabdx, dabdy, dabdz

    INTEGER                                  :: ax, ay, az, bx, by, bz, coa, &
                                                coamx, coamy, coamz, coapx, &
                                                coapy, coapz, cob, i, ipgf, &
                                                j, jpgf, la, lb, na, nb, nda
    REAL(KIND=dp)                            :: fa, fx, fy, fz

!   ---------------------------------------------------------------------------
!   *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0
    nda = 0

    dabdx=0.0_dp
    dabdy=0.0_dp
    dabdz=0.0_dp

    DO ipgf=1,npgfa

      fa = 2.0_dp*zeta(ipgf)

      nb = 0

      DO jpgf=1,npgfb

!       *** Screening ***

        IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
          DO j=nb+ncoset(lb_min-1)+1,nb+ncoset(lb_max)
            DO i=na+ncoset(la_min-1)+1,na+ncoset(la_max)
              dabdx(i,j) = 0.0_dp
              dabdy(i,j) = 0.0_dp
              dabdz(i,j) = 0.0_dp
            END DO
          END DO
          nb = nb + ncoset(lb_max)
          CYCLE
        END IF

!       *** [da/dRi|O|b] = 2*zeta*[a+1i|O|b] - Ni(a)[a-1i|O|b] ***

        DO la=MAX(0,la_min),la_max

          IF (la == 0) THEN

            coa = na + 1
            coapx = nda + 2
            coapy = nda + 3
            coapz = nda + 4

            DO lb=lb_min,lb_max
              DO bx=0,lb
                DO by=0,lb-bx
                  bz = lb - bx - by
                  cob = nb + coset(bx,by,bz)
                  dabdx(coa,cob) = fa*ab(coapx,cob)
                  dabdy(coa,cob) = fa*ab(coapy,cob)
                  dabdz(coa,cob) = fa*ab(coapz,cob)
                END DO
              END DO
            END DO

          ELSE

            DO ax=0,la
              DO ay=0,la-ax
                az = la - ax - ay

                coa = na + coset(ax,ay,az)
                coamx = nda + coset(MAX(0,ax-1),ay,az)
                coamy = nda + coset(ax,MAX(0,ay-1),az)
                coamz = nda + coset(ax,ay,MAX(0,az-1))
                coapx = nda + coset(ax+1,ay,az)
                coapy = nda + coset(ax,ay+1,az)
                coapz = nda + coset(ax,ay,az+1)

                fx = REAL(ax,dp)
                fy = REAL(ay,dp)
                fz = REAL(az,dp)

                DO lb=lb_min,lb_max
                  DO bx=0,lb
                    DO by=0,lb-bx
                      bz = lb - bx - by
                      cob = nb + coset(bx,by,bz)
                      dabdx(coa,cob) = fa*ab(coapx,cob) - fx*ab(coamx,cob)
                      dabdy(coa,cob) = fa*ab(coapy,cob) - fy*ab(coamy,cob)
                      dabdz(coa,cob) = fa*ab(coapz,cob) - fz*ab(coamz,cob)
                    END DO
                  END DO
                END DO

              END DO
            END DO

          END IF

        END DO

        nb = nb + ncoset(lb_max)

      END DO

      na = na + ncoset(la_max)
      nda = nda + ncoset(la_max+1)

    END DO

  END SUBROUTINE dabdr

! *****************************************************************************

END MODULE ai_derivatives


