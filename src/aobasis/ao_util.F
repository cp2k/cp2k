!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief All kind of helpful little routines
!> \par History
!>      none
!> \author CJM & JGH
! **************************************************************************************************
MODULE ao_util

   USE kinds,                           ONLY: dp,&
                                              sp
   USE mathconstants,                   ONLY: &
        dfac, fac, fourpi, oorootpi, rootpi, sqrt105, sqrt15, sqrt2, sqrt21, sqrt3, sqrt35, sqrt5, &
        sqrt7, sqrthalf
   USE orbital_pointers,                ONLY: coset,&
                                              indco,&
                                              nco,&
                                              ncoset,&
                                              nso,&
                                              nsoset
   USE orbital_transformation_matrices, ONLY: orbtramat
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ao_util'

   !MK sqrt* constants moved to mathconstants

   REAL(KIND=dp), PARAMETER :: s_root1o4pi = 0.5_dp*oorootpi
   REAL(KIND=dp), PARAMETER :: root4pi = 2.0_dp*rootpi
   REAL(KIND=dp), PARAMETER :: s_root3o4pi = sqrt3*s_root1o4pi
   REAL(KIND=dp), PARAMETER :: root4pio3 = root4pi/sqrt3
   REAL(KIND=dp), PARAMETER :: root4pio5 = root4pi/sqrt5
   REAL(KIND=dp), PARAMETER :: s_root15o4pi = sqrt15*s_root1o4pi
   REAL(KIND=dp), PARAMETER :: root4pio15 = root4pi/sqrt15
   REAL(KIND=dp), PARAMETER :: s_root105o4pi = sqrt105*s_root1o4pi
   REAL(KIND=dp), PARAMETER :: root4pio105 = root4pi/sqrt105
   REAL(KIND=dp), PARAMETER :: s_root1o16pi = 0.25_dp*oorootpi
   REAL(KIND=dp), PARAMETER :: root16pi = 4.0_dp*rootpi
   REAL(KIND=dp), PARAMETER :: s_root5o16pi = sqrt5*s_root1o16pi
   REAL(KIND=dp), PARAMETER :: root16pio5 = root16pi/sqrt5
   REAL(KIND=dp), PARAMETER :: s_2root5o16pi = 2.0_dp*s_root5o16pi
   REAL(KIND=dp), PARAMETER :: root16pio5o2 = root16pio5*0.5_dp
   REAL(KIND=dp), PARAMETER :: s_3root5o16pi = 3.0_dp*s_root5o16pi
   REAL(KIND=dp), PARAMETER :: root16pio5o3 = root16pio5/3.0_dp
   REAL(KIND=dp), PARAMETER :: s_18root5o16pi = 18.0_dp*s_root5o16pi
   REAL(KIND=dp), PARAMETER :: root16pio5o18 = root16pio5/18.0_dp
   REAL(KIND=dp), PARAMETER :: s_2root7o16pi = 2.0_dp*sqrt7*s_root1o16pi
   REAL(KIND=dp), PARAMETER :: root16pio7o2 = root16pi/sqrt7*0.5_dp
   REAL(KIND=dp), PARAMETER :: s_3root7o16pi = 3.0_dp*sqrt7*s_root1o16pi
   REAL(KIND=dp), PARAMETER :: root16pio7o3 = root16pi/sqrt7/3.0_dp
   REAL(KIND=dp), PARAMETER :: s_root15o16pi = sqrt15*s_root1o16pi
   REAL(KIND=dp), PARAMETER :: root16pio15 = root16pi/sqrt15
   REAL(KIND=dp), PARAMETER :: s_3root35o16pi = sqrt5*s_3root7o16pi
   REAL(KIND=dp), PARAMETER :: root16pio35o3 = root16pio7o3/sqrt5
   REAL(KIND=dp), PARAMETER :: s_root105o16pi = 0.5_dp*s_root105o4pi
   REAL(KIND=dp), PARAMETER :: root16pio105 = root4pio105*2.0_dp
   REAL(KIND=dp), PARAMETER :: s_root1o32pi = 0.25_dp*sqrthalf*oorootpi
   REAL(KIND=dp), PARAMETER :: root32pi = root16pi*sqrt2
   REAL(KIND=dp), PARAMETER :: s_3root5o32pi = 3.0_dp*sqrt5*s_root1o32pi
   REAL(KIND=dp), PARAMETER :: root32pio5o3 = root32pi/sqrt5/3.0_dp
   REAL(KIND=dp), PARAMETER :: s_9root5o32pi = 9.0_dp*sqrt5*s_root1o32pi
   REAL(KIND=dp), PARAMETER :: root32pio5o9 = root32pi/sqrt5/9.0_dp
   REAL(KIND=dp), PARAMETER :: s_12root5o32pi = 12.0_dp*sqrt5*s_root1o32pi
   REAL(KIND=dp), PARAMETER :: root32pio5o12 = root32pi/sqrt5/12.0_dp
   REAL(KIND=dp), PARAMETER :: s_root21o32pi = sqrt21*s_root1o32pi
   REAL(KIND=dp), PARAMETER :: root32pio21 = root32pi/sqrt21
   REAL(KIND=dp), PARAMETER :: s_4root21o32pi = 4.0_dp*s_root21o32pi
   REAL(KIND=dp), PARAMETER :: root32pio21o4 = root32pio21/4.0_dp
   REAL(KIND=dp), PARAMETER :: s_root35o32pi = sqrt35*s_root1o32pi
   REAL(KIND=dp), PARAMETER :: root32pio35 = root32pi/sqrt35
   REAL(KIND=dp), PARAMETER :: s_3root35o32pi = 3.0_dp*s_root35o32pi
   REAL(KIND=dp), PARAMETER :: s_9root35o32pi = 9.0_dp*s_root35o32pi
   REAL(KIND=dp), PARAMETER :: s_18root35o32pi = 18.0_dp*s_root35o32pi
   REAL(KIND=dp), PARAMETER :: s_root1o64pi = 0.125_dp*oorootpi
   REAL(KIND=dp), PARAMETER :: s_3root5o64pi = 3.0_dp*sqrt5*s_root1o64pi
   REAL(KIND=dp), PARAMETER :: s_18root5o64pi = 18.0_dp*sqrt5*s_root1o64pi
   REAL(KIND=dp), PARAMETER :: s_root1o256pi = 0.0625_dp*oorootpi
   REAL(KIND=dp), PARAMETER :: s_3root1o256pi = 3.0_dp*s_root1o256pi
   REAL(KIND=dp), PARAMETER :: s_9root1o256pi = 9.0_dp*s_root1o256pi
   REAL(KIND=dp), PARAMETER :: s_18root1o256pi = 18.0_dp*s_root1o256pi
   REAL(KIND=dp), PARAMETER :: s_24root1o256pi = 24.0_dp*s_root1o256pi
   REAL(KIND=dp), PARAMETER :: s_72root1o256pi = 72.0_dp*s_root1o256pi
   REAL(KIND=dp), PARAMETER :: s_3root35o256pi = 3.0_dp*sqrt35*s_root1o256pi
   REAL(KIND=dp), PARAMETER :: s_18root35o256pi = 18.0_dp*sqrt35*s_root1o256pi

   ! *** Public subroutines ***

   PUBLIC :: exp_radius, &
             exp_radius_very_extended, &
             gauss_exponent, &
             gaussint_sph, &
             trace_r_AxB

CONTAINS

! **************************************************************************************************
!> \brief  The exponent of a primitive Gaussian function for a given radius
!>         and threshold is calculated.
!> \param l ...
!> \param radius ...
!> \param threshold ...
!> \param prefactor ...
!> \return ...
!> \date   07.03.1999
!> \par Variables
!>      - exponent : Exponent of the primitive Gaussian function.
!>      - l        : Angular momentum quantum number l.
!>      - prefactor: Prefactor of the Gaussian function (e.g. a contraction
!>                   coefficient).
!>      - radius   : Calculated radius of the Gaussian function.
!>      - threshold: Threshold for radius.
!> \author MK
!> \version 1.0
! **************************************************************************************************
   FUNCTION gauss_exponent(l, radius, threshold, prefactor) RESULT(exponent)
      INTEGER, INTENT(IN)                                :: l
      REAL(KIND=dp), INTENT(IN)                          :: radius, threshold, prefactor
      REAL(KIND=dp)                                      :: exponent

      exponent = 0.0_dp

      IF (radius < 1.0E-6_dp) RETURN
      IF (threshold < 1.0E-12_dp) RETURN

      exponent = LOG(ABS(prefactor)*radius**l/threshold)/radius**2

   END FUNCTION gauss_exponent

! **************************************************************************************************
!> \brief  The radius of a primitive Gaussian function for a given threshold
!>         is calculated.
!>               g(r) = prefactor*r**l*exp(-alpha*r**2) - threshold = 0
!> \param l Angular momentum quantum number l.
!> \param alpha Exponent of the primitive Gaussian function.
!> \param threshold Threshold for function g(r).
!> \param prefactor Prefactor of the Gaussian function (e.g. a contraction
!>                     coefficient).
!> \param epsabs Absolute tolerance (radius)
!> \param epsrel Relative tolerance (radius)
!> \param rlow Optional lower bound radius, e.g., when searching maximum radius
!> \return Calculated radius of the Gaussian function
!> \par Variables
!>        - g       : function g(r)
!>        - prec_exp: single/double precision exponential function
!>        - itermax : Maximum number of iterations
!>        - contract: Golden Section Search (GSS): [0.38, 0.62]
!>                    Equidistant sampling: [0.2, 0.4, 0.6, 0.8]
!>                    Bisection (original approach): [0.5]
!>                    Array size may trigger vectorization
! **************************************************************************************************
   FUNCTION exp_radius(l, alpha, threshold, prefactor, epsabs, epsrel, rlow) RESULT(radius)
      INTEGER, INTENT(IN)                                :: l
      REAL(KIND=dp), INTENT(IN)                          :: alpha, threshold, prefactor
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: epsabs, epsrel, rlow
      REAL(KIND=dp)                                      :: radius

      CHARACTER(len=*), PARAMETER :: routineN = 'exp_radius', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: itermax = 5000, prec_exp = sp
      REAL(KIND=dp), PARAMETER                           :: contract(*) = (/0.38, 0.62/), &
                                                            mineps = 1.0E-12, next = 2.0, &
                                                            step = 1.0

      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: a, d, dr, eps, r, rd, t
      REAL(KIND=dp), DIMENSION(SIZE(contract))           :: g, s

      IF (l .LT. 0) THEN
         CPABORT("The angular momentum quantum number is negative")
      ENDIF

      IF (alpha .EQ. 0.0_dp) THEN
         CPABORT("The Gaussian function exponent is zero")
      ELSE
         a = ABS(alpha)
      ENDIF

      IF (threshold .NE. 0.0_dp) THEN
         t = ABS(threshold)
      ELSE
         CPABORT("The requested threshold is zero")
      ENDIF

      radius = 0.0_dp
      IF (PRESENT(rlow)) radius = rlow
      IF (prefactor .EQ. 0.0_dp) RETURN

      ! MAX: facilitate early exit
      r = MAX(SQRT(0.5_dp*REAL(l, dp)/a), radius)

      d = ABS(prefactor); g(1) = d
      IF (l .NE. 0) THEN
         g(1) = g(1)*EXP(REAL(-a*r*r, KIND=prec_exp))*r**l
      END IF
      ! original approach may return radius=0
      ! with g(r) != g(radius)
      !radius = r
      IF (g(1) .LT. t) RETURN ! early exit

      radius = r*next + step
      DO i = 1, itermax
         g(1) = d*EXP(REAL(-a*radius*radius, KIND=prec_exp))*radius**l
         IF (g(1) .LT. t) EXIT
         r = radius; radius = r*next + step
      ENDDO

      ! consider absolute and relative accuracy (interval width)
      IF (PRESENT(epsabs)) THEN
         eps = epsabs
      ELSE IF (.NOT. PRESENT(epsrel)) THEN
         eps = mineps
      ELSE
         eps = HUGE(eps)
      ENDIF
      IF (PRESENT(epsrel)) eps = MIN(eps, epsrel*r)

      dr = 0.0_dp
      DO i = i + 1, itermax
         rd = radius - r
         ! check if finished or no further progress
         IF ((rd .LT. eps) .OR. (rd .EQ. dr)) RETURN
         s = r + rd*contract ! interval contraction
         g = d*EXP(REAL(-a*s*s, KIND=prec_exp))*s**l
         DO j = 1, SIZE(contract)
            IF (g(j) .LT. t) THEN
               radius = s(j) ! interval [r, sj)
               EXIT
            ELSE
               r = s(j) ! interval [sj, radius)
            ENDIF
         ENDDO
         dr = rd
      ENDDO
      IF (i .GE. itermax) THEN
         CPABORT("Maximum number of iterations reached")
      ENDIF

   END FUNCTION exp_radius

! **************************************************************************************************
!> \brief computes the radius of the Gaussian outside of which it is smaller
!>      than eps
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param pab ...
!> \param o1 ...
!> \param o2 ...
!> \param ra ...
!> \param rb ...
!> \param rp ...
!> \param zetp ...
!> \param eps ...
!> \param prefactor ...
!> \param cutoff ...
!> \param epsabs ...
!> \return ...
!> \par History
!>      03.2007 new version that assumes that the Gaussians origante from spherical
!>              Gaussians
!> \note
!>      can optionally screen by the maximum element of the pab block
! **************************************************************************************************
   FUNCTION exp_radius_very_extended(la_min, la_max, lb_min, lb_max, pab, o1, o2, ra, rb, rp, &
                                     zetp, eps, prefactor, cutoff, epsabs) RESULT(radius)

      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: pab
      INTEGER, OPTIONAL                                  :: o1, o2
      REAL(KIND=dp), INTENT(IN)                          :: ra(3), rb(3), rp(3), zetp, eps, &
                                                            prefactor, cutoff
      REAL(KIND=dp), OPTIONAL                            :: epsabs
      REAL(KIND=dp)                                      :: radius

      INTEGER                                            :: i, ico, j, jco, la(3), lb(3), lxa, lxb, &
                                                            lya, lyb, lza, lzb
      REAL(KIND=dp)                                      :: bini, binj, coef(0:20), epsin_local, &
                                                            polycoef(0:60), prefactor_local, &
                                                            rad_a, rad_b, s1, s2

! get the local prefactor, we'll now use the largest density matrix element of the block to screen

      epsin_local = 1.0E-2_dp
      IF (PRESENT(epsabs)) epsin_local = epsabs

      IF (PRESENT(pab)) THEN
         prefactor_local = cutoff
         DO lxa = 0, la_max
         DO lxb = 0, lb_max
            DO lya = 0, la_max - lxa
            DO lyb = 0, lb_max - lxb
               DO lza = MAX(la_min - lxa - lya, 0), la_max - lxa - lya
               DO lzb = MAX(lb_min - lxb - lyb, 0), lb_max - lxb - lyb
                  la = (/lxa, lya, lza/)
                  lb = (/lxb, lyb, lzb/)
                  ico = coset(lxa, lya, lza)
                  jco = coset(lxb, lyb, lzb)
                  prefactor_local = MAX(ABS(pab(o1 + ico, o2 + jco)), prefactor_local)
               ENDDO
               ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         prefactor_local = prefactor*prefactor_local
      ELSE
         prefactor_local = prefactor*MAX(1.0_dp, cutoff)
      ENDIF

      !
      ! assumes that we can compute the radius for the case where
      ! the Gaussians a and b are both on the z - axis, but at the same
      ! distance as the original a and b
      !
      rad_a = SQRT(SUM((ra - rp)**2))
      rad_b = SQRT(SUM((rb - rp)**2))

      polycoef(0:la_max + lb_max) = 0.0_dp
      DO lxa = 0, la_max
      DO lxb = 0, lb_max
         coef(0:la_max + lb_max) = 0.0_dp
         bini = 1.0_dp
         s1 = 1.0_dp
         DO i = 0, lxa
            binj = 1.0_dp
            s2 = 1.0_dp
            DO j = 0, lxb
               coef(lxa + lxb - i - j) = coef(lxa + lxb - i - j) + bini*binj*s1*s2
               binj = (binj*(lxb - j))/(j + 1)
               s2 = s2*(rad_b)
            ENDDO
            bini = (bini*(lxa - i))/(i + 1)
            s1 = s1*(rad_a)
         ENDDO
         DO i = 0, lxa + lxb
            polycoef(i) = MAX(polycoef(i), coef(i))
         ENDDO
      ENDDO
      ENDDO

      polycoef(0:la_max + lb_max) = polycoef(0:la_max + lb_max)*prefactor_local
      radius = 0.0_dp
      DO i = 0, la_max + lb_max
         radius = MAX(radius, exp_radius(i, zetp, eps, polycoef(i), epsin_local, rlow=radius))
      ENDDO

   END FUNCTION exp_radius_very_extended

! **************************************************************************************************
!> \brief ...
!> \param alpha ...
!> \param l ...
!> \return ...
! **************************************************************************************************
   FUNCTION gaussint_sph(alpha, l)

      !  calculates the radial integral over a spherical Gaussian
      !  of the form
      !     r**(2+l) * exp(-alpha * r**2)
      !

      REAL(dp), INTENT(IN)                               :: alpha
      INTEGER, INTENT(IN)                                :: l
      REAL(dp)                                           :: gaussint_sph

      IF ((l/2)*2 == l) THEN
         !even l:
         gaussint_sph = ROOTPI*0.5_dp**(l/2 + 2)*dfac(l + 1) &
                        /SQRT(alpha)**(l + 3)
      ELSE
         !odd l:
         gaussint_sph = 0.5_dp*fac((l + 1)/2)/SQRT(alpha)**(l + 3)
      ENDIF

   END FUNCTION gaussint_sph

! **************************************************************************************************
!> \brief ...
!> \param A ...
!> \param lda ...
!> \param B ...
!> \param ldb ...
!> \param m ...
!> \param n ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION trace_r_AxB(A, lda, B, ldb, m, n)

      INTEGER, INTENT(in)                                :: lda
      REAL(dp), INTENT(in)                               :: A(lda, *)
      INTEGER, INTENT(in)                                :: ldb
      REAL(dp), INTENT(in)                               :: B(ldb, *)
      INTEGER, INTENT(in)                                :: m, n
      REAL(dp)                                           :: trace_r_AxB

      INTEGER                                            :: i1, i2, imod, mminus3
      REAL(dp)                                           :: t1, t2, t3, t4

      t1 = 0._dp
      t2 = 0._dp
      t3 = 0._dp
      t4 = 0._dp
      imod = MODULO(m, 4)
      SELECT CASE (imod)
      CASE (0)
         DO i2 = 1, n
            DO i1 = 1, m, 4
               t1 = t1 + A(i1, i2)*B(i1, i2)
               t2 = t2 + A(i1 + 1, i2)*B(i1 + 1, i2)
               t3 = t3 + A(i1 + 2, i2)*B(i1 + 2, i2)
               t4 = t4 + A(i1 + 3, i2)*B(i1 + 3, i2)
            ENDDO
         ENDDO
      CASE (1)
         mminus3 = m - 3
         DO i2 = 1, n
            DO i1 = 1, mminus3, 4
               t1 = t1 + A(i1, i2)*B(i1, i2)
               t2 = t2 + A(i1 + 1, i2)*B(i1 + 1, i2)
               t3 = t3 + A(i1 + 2, i2)*B(i1 + 2, i2)
               t4 = t4 + A(i1 + 3, i2)*B(i1 + 3, i2)
            ENDDO
            t1 = t1 + A(m, i2)*B(m, i2)
         ENDDO
      CASE (2)
         mminus3 = m - 3
         DO i2 = 1, n
            DO i1 = 1, mminus3, 4
               t1 = t1 + A(i1, i2)*B(i1, i2)
               t2 = t2 + A(i1 + 1, i2)*B(i1 + 1, i2)
               t3 = t3 + A(i1 + 2, i2)*B(i1 + 2, i2)
               t4 = t4 + A(i1 + 3, i2)*B(i1 + 3, i2)
            ENDDO
            t1 = t1 + A(m - 1, i2)*B(m - 1, i2)
            t2 = t2 + A(m, i2)*B(m, i2)
         ENDDO
      CASE (3)
         mminus3 = m - 3
         DO i2 = 1, n
            DO i1 = 1, mminus3, 4
               t1 = t1 + A(i1, i2)*B(i1, i2)
               t2 = t2 + A(i1 + 1, i2)*B(i1 + 1, i2)
               t3 = t3 + A(i1 + 2, i2)*B(i1 + 2, i2)
               t4 = t4 + A(i1 + 3, i2)*B(i1 + 3, i2)
            ENDDO
            t1 = t1 + A(m - 2, i2)*B(m - 2, i2)
            t2 = t2 + A(m - 1, i2)*B(m - 1, i2)
            t3 = t3 + A(m, i2)*B(m, i2)
         ENDDO
      END SELECT
      trace_r_AxB = t1 + t2 + t3 + t4

   END FUNCTION trace_r_AxB

! **************************************************************************************************
!> \brief ...
!> \param CPC_co ...
!> \param CPC_so ...
!> \param maxl ...
!> \param lm1 ...
!> \param lm2 ...
! **************************************************************************************************
   SUBROUTINE transform_c2s(CPC_co, CPC_so, maxl, lm1, lm2)

      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: CPC_co
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: CPC_so
      INTEGER                                            :: maxl, lm1, lm2

      CHARACTER(len=*), PARAMETER :: routineN = 'transform_c2s', routineP = moduleN//':'//routineN

      INTEGER                                            :: ic1, ic2, is1, is2, l, lx, ly, lz
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: work

      ALLOCATE (work(ncoset(maxl), nsoset(maxl)))
      CPC_so = 0.0_dp
      work = 0.0_dp

      !    DO l = 0,maxl
      !      ncgf = nco(l)
      !      nsgf = nso(l)
      !      ico = ncoset(l-1) + 1
      !      iso = nsoset(l-1) + 1
      !
      !      CALL dgemm("N","T",ncoset(maxl),nsgf,ncgf,&
      !                 1.0_dp,CPC_co(1,ico),ncoset(maxl),&
      !                 orbtramat(l)%c2s(1,1),nsgf,&
      !                 0.0_dp,work(1,iso),ncoset(maxl))
      !    ENDDO
      !
      !    DO l = 0,maxl
      !      ncgf = nco(l)
      !      nsgf = nso(l)
      !      ico = ncoset(l-1) + 1
      !      iso = nsoset(l-1) + 1
      !
      !      CALL dgemm("N","N",nsgf,nsoset(maxl),ncgf,&
      !                 1.0_dp,orbtramat(l)%c2s(1,1),nsgf,&
      !                 work(ico,1),ncoset(maxl),&
      !                 0.0_dp,CPC_so(iso,1),nsoset(maxl))
      !    ENDDO

      !  do l = 0,maxl
      !     do is1 = 1,nso(l)
      !      write(*,'(A,2I3,10f10.5)') 'or ', l, is1, orbtramat(l)%c2s(is1,1:nco(l))
      !     enddo
      !  enddo
      !
      !  stop

      CPC_so = 0.0_dp
      work = 0.0_dp
      DO ic1 = 1, ncoset(lm1)
         DO l = 0, lm2
            DO is2 = 1, nso(l)
               DO ic2 = 1, nco(l)
                  lx = indco(1, ic2 + ncoset(l - 1))
                  ly = indco(2, ic2 + ncoset(l - 1))
                  lz = indco(3, ic2 + ncoset(l - 1))
                  work(ic1, is2 + nsoset(l - 1)) = &
                     work(ic1, is2 + nsoset(l - 1)) + &
                     CPC_co(ic1, ic2 + ncoset(l - 1))* &
                     orbtramat(l)%c2s(is2, ic2)* &
                     SQRT(fourpi/dfac(2*l + 1)* &
                          dfac(2*lx - 1)*dfac(2*ly - 1)*dfac(2*lz - 1))
                  !         write(*,*) 'dfac 1', dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DO is2 = 1, nsoset(lm2)
         DO l = 0, lm1
            DO is1 = 1, nso(l)
               DO ic1 = 1, nco(l)
                  lx = indco(1, ic1 + ncoset(l - 1))
                  ly = indco(2, ic1 + ncoset(l - 1))
                  lz = indco(3, ic1 + ncoset(l - 1))
                  CPC_so(is1 + nsoset(l - 1), is2) = &
                     CPC_so(is1 + nsoset(l - 1), is2) + &
                     work(ic1 + ncoset(l - 1), is2)* &
                     orbtramat(l)%c2s(is1, ic1)* &
                     SQRT(fourpi/dfac(2*l + 1)* &
                          dfac(2*lx - 1)*dfac(2*ly - 1)*dfac(2*lz - 1))
                  !          write(*,*) 'dfac 2', dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE (work)

   END SUBROUTINE transform_c2s

! **************************************************************************************************
!> \brief ...
!> \param CPC_co ...
!> \param CPC_so ...
!> \param maxl ...
! **************************************************************************************************
   SUBROUTINE transform_c2s_new(CPC_co, CPC_so, maxl)

      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: CPC_co
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: CPC_so
      INTEGER                                            :: maxl

      CHARACTER(len=*), PARAMETER :: routineN = 'transform_c2s_new', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iso, ldc, lds
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: work

      ALLOCATE (work(ncoset(maxl), ncoset(maxl), 1))
      ldc = ncoset(maxl)
      lds = nsoset(maxl)
      CPC_so = 0.0_dp
      work = 0.0_dp

      work(1:ldc, 1:ldc, 1) = CPC_co(1:ldc, 1:ldc)

      CALL cart2sph_mat(work, ldc, ldc, 1, maxl, maxl)

      DO iso = 1, nsoset(maxl)
         CPC_so(1:nsoset(maxl), iso) = work(1:nsoset(maxl), iso, 1)
      END DO

      DEALLOCATE (work)

   END SUBROUTINE transform_c2s_new

! **************************************************************************************************
!> \brief ...
!> \param matso ...
!> \param matco ...
!> \param maxl ...
!> \param lm1 ...
!> \param lm2 ...
! **************************************************************************************************
   SUBROUTINE transform_s2c(matso, matco, maxl, lm1, lm2)

      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: matso
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: matco
      INTEGER, INTENT(IN)                                :: maxl, lm1, lm2

      CHARACTER(len=*), PARAMETER :: routineN = 'transform_s2c', routineP = moduleN//':'//routineN

      INTEGER                                            :: ic1, ic2, ico, is1, is2, iso, l, lx, ly, &
                                                            lz, nc1, nc2, ns1, ns2
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: matsc

      ALLOCATE (matsc(nsoset(maxl), ncoset(maxl)))
      matco = 0.0_dp
      matsc = 0.0_dp

      ns1 = nsoset(lm1)
      ns2 = nsoset(lm2)
      nc1 = ncoset(lm1)
      nc2 = ncoset(lm2)

      !   A = matso (ns1 x ns2)
      !   C = Tt x A x T = matco (nc1 x nc2)
      !   B = A x T      = matsc (ns1 x nc2)

      !   Calculate B
      DO is1 = 1, ns1
         DO l = 0, lm2
            DO ico = 1, nco(l)
               ic2 = ncoset(l - 1) + ico
               lx = indco(1, ico + ncoset(l - 1))
               ly = indco(2, ico + ncoset(l - 1))
               lz = indco(3, ico + ncoset(l - 1))
               DO iso = 1, nso(l)
                  is2 = nsoset(l - 1) + iso
                  matsc(is1, ic2) = matsc(is1, ic2) + &
                                    matso(is1, is2)*orbtramat(l)%s2c(iso, ico)* &
                                    SQRT((fourpi)/dfac(2*l + 1)* &
                                         dfac(2*lx - 1)*dfac(2*ly - 1)*dfac(2*lz - 1))
               ENDDO ! iso
            ENDDO ! ico
         ENDDO ! l
      ENDDO ! is1

      !   Calculate C
      DO ic2 = 1, nc2
         DO l = 0, lm1
            DO ico = 1, nco(l)
               ic1 = ncoset(l - 1) + ico
               lx = indco(1, ico + ncoset(l - 1))
               ly = indco(2, ico + ncoset(l - 1))
               lz = indco(3, ico + ncoset(l - 1))
               DO iso = 1, nso(l)
                  is1 = nsoset(l - 1) + iso
                  matco(ic1, ic2) = matco(ic1, ic2) + &
                                    matsc(is1, ic2)*orbtramat(l)%s2c(iso, ico)* &
                                    SQRT((fourpi)/dfac(2*l + 1)* &
                                         dfac(2*lx - 1)*dfac(2*ly - 1)*dfac(2*lz - 1))
               ENDDO ! iso
            ENDDO ! ico
         ENDDO ! l
      ENDDO ! ic2

      DEALLOCATE (matsc)

   END SUBROUTINE transform_s2c
   !------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief ...
!> \param mat ...
!> \param ld_mat ...
!> \param sd_mat ...
!> \param n ...
!> \param lmax1 ...
!> \param lmax2 ...
! **************************************************************************************************
   SUBROUTINE sph2cart_mat(mat, ld_mat, sd_mat, n, lmax1, lmax2)
      INTEGER                                            :: ld_mat, sd_mat, n
      REAL(dp)                                           :: mat(ld_mat, sd_mat, n)
      INTEGER                                            :: lmax1, lmax2

      CHARACTER(len=*), PARAMETER :: routineN = 'sph2cart_mat', routineP = moduleN//':'//routineN

      INTEGER                                            :: idx_l1, idx_l2
      REAL(dp) :: mat_aux(ncoset(lmax1), nsoset(lmax2), n)

      IF (lmax1 > 3 .OR. lmax2 > 3) THEN
         CPABORT("l>3 not implemented")
      ENDIF
      DO idx_l2 = 1, nsoset(lmax2)

         mat_aux(1, idx_l2, :) &
            = s_root1o4pi*mat(1, idx_l2, :)

         IF (lmax1 == 0) CYCLE

         mat_aux(2, idx_l2, :) &
            = -s_root3o4pi*mat(4, idx_l2, :)

         mat_aux(3, idx_l2, :) &
            = -s_root3o4pi*mat(2, idx_l2, :)

         mat_aux(4, idx_l2, :) &
            = s_root3o4pi*mat(3, idx_l2, :)

         IF (lmax1 == 1) CYCLE

         mat_aux(5, idx_l2, :) &
            = -s_root5o16pi*mat(7, idx_l2, :) &
              + s_root15o16pi*mat(9, idx_l2, :)

         mat_aux(6, idx_l2, :) &
            = s_root15o4pi*mat(5, idx_l2, :)

         mat_aux(7, idx_l2, :) &
            = -s_root5o16pi*mat(7, idx_l2, :) &
              - s_root15o16pi*mat(9, idx_l2, :)

         mat_aux(8, idx_l2, :) &
            = -s_root15o4pi*mat(8, idx_l2, :)

         mat_aux(9, idx_l2, :) &
            = -s_root15o4pi*mat(6, idx_l2, :)

         mat_aux(10, idx_l2, :) &
            = s_2root5o16pi*mat(7, idx_l2, :)

         IF (lmax1 == 2) CYCLE

         mat_aux(11, idx_l2, :) &
            = -s_root21o32pi*mat(14, idx_l2, :) &
              + s_root35o32pi*mat(16, idx_l2, :)

         mat_aux(12, idx_l2, :) &
            = -s_3root35o32pi*mat(10, idx_l2, :) &
              + s_root21o32pi*mat(12, idx_l2, :)

         mat_aux(13, idx_l2, :) &
            = s_root21o32pi*mat(14, idx_l2, :) &
              + s_3root35o32pi*mat(16, idx_l2, :)

         mat_aux(14, idx_l2, :) &
            = s_root35o32pi*mat(10, idx_l2, :) &
              + s_root21o32pi*mat(12, idx_l2, :)

         mat_aux(15, idx_l2, :) &
            = -s_3root7o16pi*mat(13, idx_l2, :) &
              + s_root105o16pi*mat(15, idx_l2, :)

         mat_aux(16, idx_l2, :) &
            = s_root105o4pi*mat(11, idx_l2, :)

         mat_aux(17, idx_l2, :) &
            = -s_3root7o16pi*mat(13, idx_l2, :) &
              - s_root105o16pi*mat(15, idx_l2, :)

         mat_aux(18, idx_l2, :) &
            = -s_4root21o32pi*mat(14, idx_l2, :)

         mat_aux(19, idx_l2, :) &
            = -s_4root21o32pi*mat(12, idx_l2, :)

         mat_aux(20, idx_l2, :) &
            = s_2root7o16pi*mat(13, idx_l2, :)

      ENDDO

      DO idx_l1 = 1, ncoset(lmax1)

         mat(idx_l1, 1, :) &
            = s_root1o4pi*mat_aux(idx_l1, 1, :)

         IF (lmax2 == 0) CYCLE

         mat(idx_l1, 2, :) &
            = -s_root3o4pi*mat_aux(idx_l1, 4, :)

         mat(idx_l1, 3, :) &
            = -s_root3o4pi*mat_aux(idx_l1, 2, :)

         mat(idx_l1, 4, :) &
            = s_root3o4pi*mat_aux(idx_l1, 3, :)

         IF (lmax2 == 1) CYCLE

         mat(idx_l1, 5, :) &
            = -s_root5o16pi*mat_aux(idx_l1, 7, :) &
              + s_root15o16pi*mat_aux(idx_l1, 9, :)

         mat(idx_l1, 6, :) &
            = s_root15o4pi*mat_aux(idx_l1, 5, :)

         mat(idx_l1, 7, :) &
            = -s_root5o16pi*mat_aux(idx_l1, 7, :) &
              - s_root15o16pi*mat_aux(idx_l1, 9, :)

         mat(idx_l1, 8, :) &
            = -s_root15o4pi*mat_aux(idx_l1, 8, :)

         mat(idx_l1, 9, :) &
            = -s_root15o4pi*mat_aux(idx_l1, 6, :)

         mat(idx_l1, 10, :) &
            = s_2root5o16pi*mat_aux(idx_l1, 7, :)

         IF (lmax2 == 2) CYCLE

         mat(idx_l1, 11, :) &
            = -s_root21o32pi*mat_aux(idx_l1, 14, :) &
              + s_root35o32pi*mat_aux(idx_l1, 16, :)

         mat(idx_l1, 12, :) &
            = -s_3root35o32pi*mat_aux(idx_l1, 10, :) &
              + s_root21o32pi*mat_aux(idx_l1, 12, :)

         mat(idx_l1, 13, :) &
            = s_root21o32pi*mat_aux(idx_l1, 14, :) &
              + s_3root35o32pi*mat_aux(idx_l1, 16, :)

         mat(idx_l1, 14, :) &
            = s_root35o32pi*mat_aux(idx_l1, 10, :) &
              + s_root21o32pi*mat_aux(idx_l1, 12, :)

         mat(idx_l1, 15, :) &
            = -s_3root7o16pi*mat_aux(idx_l1, 13, :) &
              + s_root105o16pi*mat_aux(idx_l1, 15, :)

         mat(idx_l1, 16, :) &
            = s_root105o4pi*mat_aux(idx_l1, 11, :)

         mat(idx_l1, 17, :) &
            = -s_3root7o16pi*mat_aux(idx_l1, 13, :) &
              - s_root105o16pi*mat_aux(idx_l1, 15, :)

         mat(idx_l1, 18, :) &
            = -s_4root21o32pi*mat_aux(idx_l1, 14, :)

         mat(idx_l1, 19, :) &
            = -s_4root21o32pi*mat_aux(idx_l1, 12, :)

         mat(idx_l1, 20, :) &
            = s_2root7o16pi*mat_aux(idx_l1, 13, :)

      ENDDO

   END SUBROUTINE sph2cart_mat

! **************************************************************************************************
!> \brief ...
!> \param mat ...
!> \param ld_mat ...
!> \param sd_mat ...
!> \param n ...
!> \param lmax1 ...
!> \param lmax2 ...
! **************************************************************************************************
   SUBROUTINE cart2sph_mat(mat, ld_mat, sd_mat, n, lmax1, lmax2)

      !in:
      INTEGER                                            :: ld_mat, sd_mat, n
      REAL(dp)                                           :: mat(ld_mat, sd_mat, n)
      INTEGER                                            :: lmax1, lmax2

      CHARACTER(len=*), PARAMETER :: routineN = 'cart2sph_mat', routineP = moduleN//':'//routineN

      INTEGER                                            :: idx_l2, idx_lm1
      REAL(dp) :: mat_aux(nsoset(lmax1), ncoset(lmax2), n)

      IF ((lmax1 > 3) .OR. (lmax2 > 3)) THEN
         CPABORT("l > 3 not implemented")
      ENDIF

      DO idx_l2 = 1, ncoset(lmax2)

         mat_aux(1, idx_l2, :) &
            = root4pi*mat(1, idx_l2, :)

         IF (lmax1 == 0) CYCLE

         mat_aux(2, idx_l2, :) &
            = -root4pio3*mat(3, idx_l2, :)

         mat_aux(3, idx_l2, :) &
            = root4pio3*mat(4, idx_l2, :)

         mat_aux(4, idx_l2, :) &
            = -root4pio3*mat(2, idx_l2, :)

         IF (lmax1 == 1) CYCLE

         mat_aux(5, idx_l2, :) &
            = root4pio15*mat(6, idx_l2, :)

         mat_aux(6, idx_l2, :) &
            = -root4pio15*mat(9, idx_l2, :)

         mat_aux(7, idx_l2, :) &
            = -0.5_dp*root4pio5*mat(5, idx_l2, :) &
              - 0.5_dp*root4pio5*mat(7, idx_l2, :) &
              + root4pio5*mat(10, idx_l2, :)

         mat_aux(8, idx_l2, :) &
            = -root4pio15*mat(8, idx_l2, :)

         mat_aux(9, idx_l2, :) &
            = 0.5_dp*root4pio15*mat(5, idx_l2, :) &
              - 0.5_dp*root4pio15*mat(7, idx_l2, :)

         IF (lmax1 == 2) CYCLE

         mat_aux(10, idx_l2, :) &
            = -s_3root35o32pi*mat(12, idx_l2, :) &
              + s_root35o32pi*mat(14, idx_l2, :)

         mat_aux(11, idx_l2, :) &
            = s_root105o4pi*mat(16, idx_l2, :)

         mat_aux(12, idx_l2, :) &
            = s_root21o32pi*mat(12, idx_l2, :) &
              + s_root21o32pi*mat(14, idx_l2, :) &
              - s_4root21o32pi*mat(19, idx_l2, :)

         mat_aux(13, idx_l2, :) &
            = -s_3root7o16pi*mat(15, idx_l2, :) &
              - s_3root7o16pi*mat(17, idx_l2, :) &
              + s_2root7o16pi*mat(20, idx_l2, :)

         mat_aux(14, idx_l2, :) &
            = s_root21o32pi*mat(11, idx_l2, :) &
              + s_root21o32pi*mat(13, idx_l2, :) &
              - s_4root21o32pi*mat(18, idx_l2, :)

         mat_aux(15, idx_l2, :) &
            = s_root105o16pi*mat(15, idx_l2, :) &
              - s_root105o16pi*mat(17, idx_l2, :)

         mat_aux(16, idx_l2, :) &
            = -s_root35o32pi*mat(11, idx_l2, :) &
              + s_3root35o32pi*mat(13, idx_l2, :)

      ENDDO

      DO idx_lm1 = 1, nsoset(lmax1)

         mat(idx_lm1, 1, :) &
            = root4pi*mat_aux(idx_lm1, 1, :)

         IF (lmax2 == 0) CYCLE

         mat(idx_lm1, 2, :) &
            = -root4pio3*mat_aux(idx_lm1, 3, :)

         mat(idx_lm1, 3, :) &
            = root4pio3*mat_aux(idx_lm1, 4, :)

         mat(idx_lm1, 4, :) &
            = -root4pio3*mat_aux(idx_lm1, 2, :)

         IF (lmax2 == 1) CYCLE

         mat(idx_lm1, 5, :) &
            = root4pio15*mat_aux(idx_lm1, 6, :)

         mat(idx_lm1, 6, :) &
            = -root4pio15*mat_aux(idx_lm1, 9, :)

         mat(idx_lm1, 7, :) &
            = -0.5_dp*root4pio5*mat_aux(idx_lm1, 5, :) &
              - 0.5_dp*root4pio5*mat_aux(idx_lm1, 7, :) &
              + root4pio5*mat_aux(idx_lm1, 10, :)

         mat(idx_lm1, 8, :) &
            = -root4pio15*mat_aux(idx_lm1, 8, :)

         mat(idx_lm1, 9, :) &
            = 0.5_dp*root4pio15*mat_aux(idx_lm1, 5, :) &
              - 0.5_dp*root4pio15*mat_aux(idx_lm1, 7, :)

         IF (lmax2 == 2) CYCLE

         mat(idx_lm1, 10, :) &
            = -s_3root35o32pi*mat_aux(idx_lm1, 12, :) &
              + s_root35o32pi*mat_aux(idx_lm1, 14, :)

         mat(idx_lm1, 11, :) &
            = s_root105o4pi*mat_aux(idx_lm1, 16, :)

         mat(idx_lm1, 12, :) &
            = s_root21o32pi*mat_aux(idx_lm1, 12, :) &
              + s_root21o32pi*mat_aux(idx_lm1, 14, :) &
              - s_4root21o32pi*mat_aux(idx_lm1, 19, :)

         mat(idx_lm1, 13, :) &
            = -s_3root7o16pi*mat_aux(idx_lm1, 15, :) &
              - s_3root7o16pi*mat_aux(idx_lm1, 17, :) &
              + s_2root7o16pi*mat_aux(idx_lm1, 20, :)

         mat(idx_lm1, 14, :) &
            = s_root21o32pi*mat_aux(idx_lm1, 11, :) &
              + s_root21o32pi*mat_aux(idx_lm1, 13, :) &
              - s_4root21o32pi*mat_aux(idx_lm1, 18, :)

         mat(idx_lm1, 15, :) &
            = s_root105o16pi*mat_aux(idx_lm1, 15, :) &
              - s_root105o16pi*mat_aux(idx_lm1, 17, :)

         mat(idx_lm1, 16, :) &
            = -s_root35o32pi*mat_aux(idx_lm1, 11, :) &
              + s_3root35o32pi*mat_aux(idx_lm1, 13, :)

      ENDDO

   END SUBROUTINE cart2sph_mat

END MODULE ao_util

