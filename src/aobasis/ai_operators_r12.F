!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of integrals over Cartesian Gaussian-type functions for different r12
!>        operators: 1/r12, erf(omega*r12/r12), erfc(omega*r12/r12), exp(-omega*r12^2)/r12 and
!>                   exp(-omega*r12^2)
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!>      R. Ahlrichs, PCCP, 8, 3072 (2006)
!> \par History
!>      none
!> \par Parameters
!>       - ax,ay,az    : Angular momentum index numbers of orbital a.
!>       - cx,cy,cz    : Angular momentum index numbers of orbital c.
!>       - coset       : Cartesian orbital set pointer.
!>       - dac         : Distance between the atomic centers a and c.
!>       - l{a,c}      : Angular momentum quantum number of shell a or c.
!>       - l{a,c}_max  : Maximum angular momentum quantum number of shell a or c.
!>       - l{a,c}_min  : Minimum angular momentum quantum number of shell a or c.
!>       - ncoset      : Number of orbitals in a Cartesian orbital set.
!>       - npgf{a,c}   : Degree of contraction of shell a or c.
!>       - rac         : Distance vector between the atomic centers a and c.
!>       - rac2        : Square of the distance between the atomic centers a and c.
!>       - zet{a,c}    : Exponents of the Gaussian-type functions a or c.
!>       - zetp        : Reciprocal of the sum of the exponents of orbital a and b.
!>       - zetw        : Reciprocal of the sum of the exponents of orbital a and c.
!>       - omega       : Parameter in the operator
!> \author Dorothea Golze (05.2016)
! **************************************************************************************************
MODULE ai_operators_r12

   USE gamma,                           ONLY: fgamma=>fgamma_0
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fac,&
                                              pi
   USE orbital_pointers,                ONLY: coset,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_operators_r12'
   PRIVATE

   ! *** Public subroutines ***

   PUBLIC :: operator2, cps_coulomb2, cps_verf2, cps_verfc2, cps_vgauss2, cps_gauss2, ab_sint_os

   ABSTRACT INTERFACE
! **************************************************************************************************
!> \brief Interface for the calculation of integrals over s-functions and the s-type auxiliary
!>        integrals using the Obara-Saika (OS) scheme
!> \param v ...
!> \param nmax ...
!> \param zetp ...
!> \param zetq ...
!> \param zetw ...
!> \param rho ...
!> \param rac2 ...
!> \param omega ...
! **************************************************************************************************
      SUBROUTINE ab_sint_os(v, nmax, zetp, zetq, zetw, rho, rac2, omega)
       USE kinds,                           ONLY: dp
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: v
      INTEGER, INTENT(IN)                                :: nmax
      REAL(KIND=dp), INTENT(IN)                          :: zetp, zetq, zetw, rho, rac2, omega

      END SUBROUTINE ab_sint_os
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Calculation of the primitive two-center integrals over Cartesian Gaussian-type
!>        functions for different r12 operators.
!> \param cps_operator2 procedure pointer for the respective operator. The integrals evaluation
!>        differs only in the evaluation of the cartesian primitive s (cps) integrals [s|O(r12)|s]
!>        and auxiliary integrals [s|O(r12)|s]^n. This pointer selects the correct routine.
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param la_min ...
!> \param lc_max ...
!> \param npgfc ...
!> \param zetc ...
!> \param lc_min ...
!> \param omega ...
!> \param rac ...
!> \param rac2 ...
!> \param vac matrix storing the integrals
!> \param v temporary work array
!> \param maxder maximal derivative
!> \param vac_plus matrix storing the integrals for highler l-quantum numbers; used to
!>        construct the derivatives
! **************************************************************************************************

   SUBROUTINE operator2(cps_operator2, la_max, npgfa, zeta, la_min, lc_max, npgfc, zetc, lc_min, &
                        omega, rac, rac2, vac, v, maxder, vac_plus)
      PROCEDURE(ab_sint_os), POINTER                     :: cps_operator2
      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta
      INTEGER, INTENT(IN)                                :: la_min, lc_max, npgfc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetc
      INTEGER, INTENT(IN)                                :: lc_min
      REAL(KIND=dp), INTENT(IN)                          :: omega
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac
      REAL(KIND=dp), INTENT(IN)                          :: rac2
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: vac
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: v
      INTEGER, INTENT(IN), OPTIONAL                      :: maxder
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL           :: vac_plus

      CHARACTER(len=*), PARAMETER :: routineN = 'operator2', routineP = moduleN//':'//routineN

      INTEGER                                            :: ax, ay, az, coc, cocx, cocy, cocz, cx, &
                                                            cy, cz, i, ipgf, j, jpgf, la, lc, &
                                                            maxder_local, n, na, nap, nc, ncp, &
                                                            nmax, handle
      REAL(KIND=dp)                                      :: dac, f1, f2, f3, f4, f5, f6, fcx, &
                                                            fcy, fcz, rho, zetp, zetq, zetw
      REAL(KIND=dp), DIMENSION(3)                        :: raw, rcw

      CALL timeset(routineN, handle)

      v = 0.0_dp

      maxder_local = 0
      IF (PRESENT(maxder)) THEN
         maxder_local = maxder
         vac_plus = 0.0_dp
      END IF

      nmax = la_max+lc_max+1

      ! *** Calculate the distance of the centers a and c ***

      dac = SQRT(rac2)

      ! *** Loop over all pairs of primitive Gaussian-type functions ***

      na = 0
      nap = 0

      DO ipgf = 1, npgfa

         nc = 0
         ncp = 0

         DO jpgf = 1, npgfc

            ! *** Calculate some prefactors ***

            zetp = 1.0_dp/zeta(ipgf)
            zetq = 1.0_dp/zetc(jpgf)
            zetw = 1.0_dp/(zeta(ipgf)+zetc(jpgf))

            rho = zeta(ipgf)*zetc(jpgf)*zetw

            ! *** Calculate the basic two-center integrals [s||s]{n} ***

            CALL cps_operator2(v, nmax, zetp, zetq, zetw, rho, rac2, omega)

            ! *** Vertical recurrence steps: [s||s] -> [s||c] ***

            IF (lc_max > 0) THEN

               f1 = 0.5_dp*zetq
               f2 = -rho*zetq

               rcw(:) = -zeta(ipgf)*zetw*rac(:)

               ! *** [s||p]{n} = (Wi - Ci)*[s||s]{n+1}  (i = x,y,z) ***

               DO n = 1, nmax-1
                  v(1, 2, n) = rcw(1)*v(1, 1, n+1)
                  v(1, 3, n) = rcw(2)*v(1, 1, n+1)
                  v(1, 4, n) = rcw(3)*v(1, 1, n+1)
               END DO

               ! ** [s||c]{n} = (Wi - Ci)*[s||c-1i]{n+1} +     ***
               ! **             f1*Ni(c-1i)*(   [s||c-2i]{n} + ***
               ! **                          f2*[s||c-2i]{n+1} ***

               DO lc = 2, lc_max

                  DO n = 1, nmax-lc

                     ! **** Increase the angular momentum component z of c ***

                     v(1, coset(0, 0, lc), n) = &
                        rcw(3)*v(1, coset(0, 0, lc-1), n+1)+ &
                        f1*REAL(lc-1, dp)*(v(1, coset(0, 0, lc-2), n)+ &
                                           f2*v(1, coset(0, 0, lc-2), n+1))

                     ! *** Increase the angular momentum component y of c ***

                     cz = lc-1
                     v(1, coset(0, 1, cz), n) = rcw(2)*v(1, coset(0, 0, cz), n+1)

                     DO cy = 2, lc
                        cz = lc-cy
                        v(1, coset(0, cy, cz), n) = &
                           rcw(2)*v(1, coset(0, cy-1, cz), n+1)+ &
                           f1*REAL(cy-1, dp)*(v(1, coset(0, cy-2, cz), n)+ &
                                              f2*v(1, coset(0, cy-2, cz), n+1))
                     END DO

                     ! *** Increase the angular momentum component x of c ***

                     DO cy = 0, lc-1
                        cz = lc-1-cy
                        v(1, coset(1, cy, cz), n) = rcw(1)*v(1, coset(0, cy, cz), n+1)
                     END DO

                     DO cx = 2, lc
                        f6 = f1*REAL(cx-1, dp)
                        DO cy = 0, lc-cx
                           cz = lc-cx-cy
                           v(1, coset(cx, cy, cz), n) = &
                              rcw(1)*v(1, coset(cx-1, cy, cz), n+1)+ &
                              f6*(v(1, coset(cx-2, cy, cz), n)+ &
                                  f2*v(1, coset(cx-2, cy, cz), n+1))
                        END DO
                     END DO

                  END DO

               END DO

            END IF

            ! *** Vertical recurrence steps: [s||c] -> [a||c] ***

            IF (la_max > 0) THEN

               f3 = 0.5_dp*zetp
               f4 = -rho*zetp
               f5 = 0.5_dp*zetw

               raw(:) = zetc(jpgf)*zetw*rac(:)

               ! *** [p||s]{n} = (Wi - Ai)*[s||s]{n+1}  (i = x,y,z) ***

               DO n = 1, nmax-1
                  v(2, 1, n) = raw(1)*v(1, 1, n+1)
                  v(3, 1, n) = raw(2)*v(1, 1, n+1)
                  v(4, 1, n) = raw(3)*v(1, 1, n+1)
               END DO

               ! *** [a||s]{n} = (Wi - Ai)*[a-1i||s]{n+1} +      ***
               ! ***             f3*Ni(a-1i)*(   [a-2i||s]{n} +  ***
               ! ***                          f4*[a-2i||s]{n+1}) ***

               DO la = 2, la_max

                  DO n = 1, nmax-la

                     ! *** Increase the angular momentum component z of a ***

                     v(coset(0, 0, la), 1, n) = &
                        raw(3)*v(coset(0, 0, la-1), 1, n+1)+ &
                        f3*REAL(la-1, dp)*(v(coset(0, 0, la-2), 1, n)+ &
                                           f4*v(coset(0, 0, la-2), 1, n+1))

                     ! *** Increase the angular momentum component y of a ***

                     az = la-1
                     v(coset(0, 1, az), 1, n) = raw(2)*v(coset(0, 0, az), 1, n+1)

                     DO ay = 2, la
                        az = la-ay
                        v(coset(0, ay, az), 1, n) = &
                           raw(2)*v(coset(0, ay-1, az), 1, n+1)+ &
                           f3*REAL(ay-1, dp)*(v(coset(0, ay-2, az), 1, n)+ &
                                              f4*v(coset(0, ay-2, az), 1, n+1))
                     END DO

                     ! *** Increase the angular momentum component x of a ***

                     DO ay = 0, la-1
                        az = la-1-ay
                        v(coset(1, ay, az), 1, n) = raw(1)*v(coset(0, ay, az), 1, n+1)
                     END DO

                     DO ax = 2, la
                        f6 = f3*REAL(ax-1, dp)
                        DO ay = 0, la-ax
                           az = la-ax-ay
                           v(coset(ax, ay, az), 1, n) = &
                              raw(1)*v(coset(ax-1, ay, az), 1, n+1)+ &
                              f6*(v(coset(ax-2, ay, az), 1, n)+ &
                                  f4*v(coset(ax-2, ay, az), 1, n+1))
                        END DO
                     END DO

                  END DO

               END DO

               DO lc = 1, lc_max

                  DO cx = 0, lc
                     DO cy = 0, lc-cx
                        cz = lc-cx-cy

                        coc = coset(cx, cy, cz)
                        cocx = coset(MAX(0, cx-1), cy, cz)
                        cocy = coset(cx, MAX(0, cy-1), cz)
                        cocz = coset(cx, cy, MAX(0, cz-1))

                        fcx = f5*REAL(cx, dp)
                        fcy = f5*REAL(cy, dp)
                        fcz = f5*REAL(cz, dp)

                        ! *** [p||c]{n} = (Wi - Ai)*[s||c]{n+1} + ***
                        ! ***             f5*Ni(c)*[s||c-1i]{n+1} ***

                        DO n = 1, nmax-1-lc
                           v(2, coc, n) = raw(1)*v(1, coc, n+1)+fcx*v(1, cocx, n+1)
                           v(3, coc, n) = raw(2)*v(1, coc, n+1)+fcy*v(1, cocy, n+1)
                           v(4, coc, n) = raw(3)*v(1, coc, n+1)+fcz*v(1, cocz, n+1)
                        END DO

                        ! *** [a||c]{n} = (Wi - Ai)*[a-1i||c]{n+1} +        ***
                        ! ***             f3*Ni(a-1i)*(   [a-2i||c]{n} +    ***
                        ! ***                          f4*[a-2i||c]{n+1}) + ***
                        ! ***             f5*Ni(c)*[a-1i||c-1i]{n+1}        ***

                        DO la = 2, la_max

                           DO n = 1, nmax-la-lc

                              ! *** Increase the angular momentum component z of a ***

                              v(coset(0, 0, la), coc, n) = &
                                 raw(3)*v(coset(0, 0, la-1), coc, n+1)+ &
                                 f3*REAL(la-1, dp)*(v(coset(0, 0, la-2), coc, n)+ &
                                                    f4*v(coset(0, 0, la-2), coc, n+1))+ &
                                 fcz*v(coset(0, 0, la-1), cocz, n+1)

                              ! *** Increase the angular momentum component y of a ***

                              az = la-1
                              v(coset(0, 1, az), coc, n) = &
                                 raw(2)*v(coset(0, 0, az), coc, n+1)+ &
                                 fcy*v(coset(0, 0, az), cocy, n+1)

                              DO ay = 2, la
                                 az = la-ay
                                 v(coset(0, ay, az), coc, n) = &
                                    raw(2)*v(coset(0, ay-1, az), coc, n+1)+ &
                                    f3*REAL(ay-1, dp)*(v(coset(0, ay-2, az), coc, n)+ &
                                                       f4*v(coset(0, ay-2, az), coc, n+1))+ &
                                    fcy*v(coset(0, ay-1, az), cocy, n+1)
                              END DO

                              ! *** Increase the angular momentum component x of a ***

                              DO ay = 0, la-1
                                 az = la-1-ay
                                 v(coset(1, ay, az), coc, n) = &
                                    raw(1)*v(coset(0, ay, az), coc, n+1)+ &
                                    fcx*v(coset(0, ay, az), cocx, n+1)
                              END DO

                              DO ax = 2, la
                                 f6 = f3*REAL(ax-1, dp)
                                 DO ay = 0, la-ax
                                    az = la-ax-ay
                                    v(coset(ax, ay, az), coc, n) = &
                                       raw(1)*v(coset(ax-1, ay, az), coc, n+1)+ &
                                       f6*(v(coset(ax-2, ay, az), coc, n)+ &
                                           f4*v(coset(ax-2, ay, az), coc, n+1))+ &
                                       fcx*v(coset(ax-1, ay, az), cocx, n+1)
                                 END DO
                              END DO

                           END DO

                        END DO

                     END DO
                  END DO

               END DO

            END IF

            DO j = ncoset(lc_min-1)+1, ncoset(lc_max-maxder_local)
               DO i = ncoset(la_min-1)+1, ncoset(la_max-maxder_local)
                  vac(na+i, nc+j) = v(i, j, 1)
               END DO
            END DO

            IF (PRESENT(maxder)) THEN
               DO j = 1, ncoset(lc_max)
                  DO i = 1, ncoset(la_max)
                     vac_plus(nap+i, ncp+j) = v(i, j, 1)
                  END DO
               END DO
            END IF

            nc = nc+ncoset(lc_max-maxder_local)
            ncp = ncp+ncoset(lc_max)

         END DO

         na = na+ncoset(la_max-maxder_local)
         nap = nap+ncoset(la_max)

      END DO

      CALL timestop(handle)

   END SUBROUTINE operator2

! **************************************************************************************************
!> \brief Calculation of Coulomb integrals for s-function, i.e, [s|1/r12|s], and the auxiliary
!>        integrals [s|1/r12|s]^n
!> \param v matrix storing the integrals
!> \param nmax maximal n in the auxiliary integrals [s|1/r12|s]^n
!> \param zetp = 1/zeta
!> \param zetq = 1/zetc
!> \param zetw = 1/(zeta+zetc)
!> \param rho = zeta*zetc*zetw
!> \param rac2 square distance between center A and C, |Ra-Rc|^2
!> \param omega this parameter is actually not used, but included for the sake of the abstract
!>        interface
! **************************************************************************************************
   SUBROUTINE cps_coulomb2(v, nmax, zetp, zetq, zetw, rho, rac2, omega)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: v
      INTEGER, INTENT(IN)                                :: nmax
      REAL(KIND=dp), INTENT(IN)                          :: zetp, zetq, zetw, rho, rac2, omega

      CHARACTER(len=*), PARAMETER :: routineN = 'cps_coulomb2', routineP = moduleN//':'//routineN

      INTEGER                                            :: n
      REAL(KIND=dp)                                      :: dummy, f0, t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: f

      ALLOCATE (f(0:nmax))
      dummy = omega
      f0 = 2.0_dp*SQRT(pi**5*zetw)*zetp*zetq

      ! *** Calculate the incomplete Gamma/Boys function ***

      t = rho*rac2
      CALL fgamma(nmax-1, t, f)

      ! *** Calculate the basic two-center integrals [s||s]{n} ***

      DO n = 1, nmax
         v(1, 1, n) = f0*f(n-1)
      END DO

      DEALLOCATE (f)
   END SUBROUTINE cps_coulomb2

! **************************************************************************************************
!> \brief Calculation of verf integrals for s-function, i.e, [s|erf(omega*r12)/r12|s], and the
!>        auxiliary integrals [s|erf(omega*r12)/r12|s]^n
!> \param v matrix storing the integrals
!> \param nmax maximal n in the auxiliary integrals [s|erf(omega*r12)/r12|s]^n
!> \param zetp = 1/zeta
!> \param zetq = 1/zetc
!> \param zetw = 1/(zeta+zetc)
!> \param rho = zeta*zetc*zetw
!> \param rac2 square distance between center A and C, |Ra-Rc|^2
!> \param omega parameter in the operator
! **************************************************************************************************
   SUBROUTINE cps_verf2(v, nmax, zetp, zetq, zetw, rho, rac2, omega)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: v
      INTEGER, INTENT(IN)                                :: nmax
      REAL(KIND=dp), INTENT(IN)                          :: zetp, zetq, zetw, rho, rac2, omega

      CHARACTER(len=*), PARAMETER :: routineN = 'cps_verf2', routineP = moduleN//':'//routineN

      INTEGER                                            :: n
      REAL(KIND=dp)                                      :: arg, comega, f0, t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: f

      ALLOCATE (f(0:nmax))
      comega = omega**2/(omega**2+rho)
      f0 = 2.0_dp*SQRT(pi**5*zetw*comega)*zetp*zetq

      ! *** Calculate the incomplete Gamma/Boys function ***

      t = rho*rac2
      arg = comega*t
      CALL fgamma(nmax-1, arg, f)

      ! *** Calculate the basic two-center integrals [s||s]{n} ***

      DO n = 1, nmax
         v(1, 1, n) = f0*f(n-1)*comega**(n-1)
      END DO

      DEALLOCATE (f)

   END SUBROUTINE cps_verf2

! **************************************************************************************************
!> \brief Calculation of verfc integrals for s-function, i.e, [s|erfc(omega*r12)/r12|s], and
!>        the auxiliary integrals [s|erfc(omega*r12)/r12|s]^n
!> \param v matrix storing the integrals
!> \param nmax maximal n in the auxiliary integrals [s|erfc(omega*r12)/r12|s]^n
!> \param zetp = 1/zeta
!> \param zetq = 1/zetc
!> \param zetw = 1/(zeta+zetc)
!> \param rho = zeta*zetc*zetw
!> \param rac2 square distance between center A and C, |Ra-Rc|^2
!> \param omega parameter in the operator
! **************************************************************************************************
   SUBROUTINE cps_verfc2(v, nmax, zetp, zetq, zetw, rho, rac2, omega)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: v
      INTEGER, INTENT(IN)                                :: nmax
      REAL(KIND=dp), INTENT(IN)                          :: zetp, zetq, zetw, rho, rac2, omega

      CHARACTER(len=*), PARAMETER :: routineN = 'cps_verfc2', routineP = moduleN//':'//routineN

      INTEGER                                            :: n
      REAL(KIND=dp)                                      :: argerf, comega, f0, t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: fv, fverf

      ALLOCATE (fv(0:nmax), fverf(0:nmax))
      comega = omega**2/(omega**2+rho)
      f0 = 2.0_dp*SQRT(pi**5*zetw)*zetp*zetq

      ! *** Calculate the incomplete Gamma/Boys function ***

      t = rho*rac2
      argerf = comega*t

      CALL fgamma(nmax-1, t, fv)
      CALL fgamma(nmax-1, argerf, fverf)

      ! *** Calculate the basic two-center integrals [s||s]{n} ***

      DO n = 1, nmax
         v(1, 1, n) = f0*(fv(n-1)-SQRT(comega)*comega**(n-1)*fverf(n-1))
      END DO

      DEALLOCATE (fv, fverf)

   END SUBROUTINE cps_verfc2

! **************************************************************************************************
!> \brief Calculation of vgauss integrals for s-function, i.e, [s|exp(-omega*r12^2)/r12|s], and
!>        the auxiliary integrals [s|exp(-omega*r12^2)/r12|s]
!> \param v matrix storing the integrals
!> \param nmax maximal n in the auxiliary integrals [s|exp(-omega*r12^2)/r12|s]
!> \param zetp = 1/zeta
!> \param zetq = 1/zetc
!> \param zetw = 1/(zeta+zetc)
!> \param rho = zeta*zetc*zetw
!> \param rac2 square distance between center A and C, |Ra-Rc|^2
!> \param omega parameter in the operator
! **************************************************************************************************
   SUBROUTINE cps_vgauss2(v, nmax, zetp, zetq, zetw, rho, rac2, omega)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: v
      INTEGER, INTENT(IN)                                :: nmax
      REAL(KIND=dp), INTENT(IN)                          :: zetp, zetq, zetw, rho, rac2, omega

      CHARACTER(len=*), PARAMETER :: routineN = 'cps_vgauss2', routineP = moduleN//':'//routineN

      INTEGER                                            :: j, n
      REAL(KIND=dp)                                      :: arg, dummy, eta, expT, f0, fsign, t, tau
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: f

      ALLOCATE (f(0:nmax))

      dummy = zetp
      dummy = zetq
      eta = rho/(rho+omega)
      tau = omega/(rho+omega)

      ! *** Calculate the incomplete Gamma/Boys function ***

      t = rho*rac2
      arg = eta*t

      CALL fgamma(nmax-1, arg, f)

      expT = EXP(-omega/(omega+rho)*t)
      f0 = 2.0_dp*SQRT(pi**5*zetw**3)/(rho+omega)*expT

      ! *** Calculate the basic two-center integrals [s||s]{n} ***
      v(1, 1, 1:nmax) = 0.0_dp
      DO n = 1, nmax
         fsign = (-1.0_dp)**(n-1)
         DO j = 0, n-1
            v(1, 1, n) = v(1, 1, n)+f0*fsign* &
                         fac(n-1)/fac(n-j-1)/fac(j)*(-tau)**(n-j-1)*(-eta)**j*f(j)
         ENDDO
      ENDDO

      DEALLOCATE (f)

   END SUBROUTINE cps_vgauss2

! **************************************************************************************************
!> \brief Calculation of gauss integrals for s-function, i.e, [s|exp(-omega*r12^2)|s], and
!>        the auxiliary integrals [s|exp(-omega*r12^2)|s]
!> \param v matrix storing the integrals
!> \param nmax maximal n in the auxiliary integrals [s|exp(-omega*r12^2)|s]
!> \param zetp = 1/zeta
!> \param zetq = 1/zetc
!> \param zetw = 1/(zeta+zetc)
!> \param rho = zeta*zetc*zetw
!> \param rac2 square distance between center A and C, |Ra-Rc|^2
!> \param omega parameter in the operator
! **************************************************************************************************
   SUBROUTINE cps_gauss2(v, nmax, zetp, zetq, zetw, rho, rac2, omega)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: v
      INTEGER, INTENT(IN)                                :: nmax
      REAL(KIND=dp), INTENT(IN)                          :: zetp, zetq, zetw, rho, rac2, omega

      CHARACTER(len=*), PARAMETER :: routineN = 'cps_gauss2', routineP = moduleN//':'//routineN

      INTEGER                                            :: n
      REAL(KIND=dp)                                      :: dummy, expT, f0, t, tau
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: f

      ALLOCATE (f(0:nmax))

      dummy = zetp
      dummy = zetq
      tau = omega/(rho+omega)
      t = rho*rac2
      expT = EXP(-tau*t)
      f0 = pi**3*SQRT(zetw**3/(rho+omega)**3)*expT

      ! *** Calculate the basic two-center integrals [s||s]{n} ***

      DO n = 1, nmax
         v(1, 1, n) = f0*tau**(n-1)
      END DO

      DEALLOCATE (f)

   END SUBROUTINE cps_gauss2

END MODULE ai_operators_r12
