!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of the moment integrals over Cartesian Gaussian-type
!>      functions.
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      none
!> \author J. Hutter (16.02.2005)
! **************************************************************************************************
MODULE ai_moments

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

   USE ai_cartesian_operators, ONLY: apply_cartesian_1d, &
                                     apply_cartesian_operator, &
                                     cartesian_polynomial_1d
   USE ai_os_rr, ONLY: os_rr_ovlp
   USE kinds, ONLY: dp
   USE mathconstants, ONLY: pi
   USE orbital_pointers, ONLY: indco, &
                               ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_moments'

   PUBLIC :: expikr, moment, diff_momop, contract_expikr, dipole_force

CONTAINS

! **************************************************************************************************
!> \brief Contract complex primitive integrals into real cosine and sine AO blocks.
!> \param cos_block ...
!> \param sin_block ...
!> \param iatom ...
!> \param ncoa ...
!> \param nsgfa ...
!> \param sgfa ...
!> \param sphi_a ...
!> \param ldsa ...
!> \param jatom ...
!> \param ncob ...
!> \param nsgfb ...
!> \param sgfb ...
!> \param sphi_b ...
!> \param ldsb ...
!> \param expab complex primitive plane-wave integrals
!> \param work ...
!> \param ldwork ...
!> \param ordered retain the ordered atom pair instead of mapping to an upper block
! **************************************************************************************************
   SUBROUTINE contract_expikr(cos_block, sin_block, &
                              iatom, ncoa, nsgfa, sgfa, sphi_a, ldsa, &
                              jatom, ncob, nsgfb, sgfb, sphi_b, ldsb, &
                              expab, work, ldwork, ordered)

      REAL(dp), DIMENSION(:, :), POINTER                 :: cos_block, sin_block
      INTEGER, INTENT(IN)                                :: iatom, ncoa, nsgfa, sgfa
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: sphi_a
      INTEGER, INTENT(IN)                                :: ldsa, jatom, ncob, nsgfb, sgfb
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: sphi_b
      INTEGER, INTENT(IN)                                :: ldsb
      COMPLEX(dp), DIMENSION(:, :), INTENT(IN)           :: expab
      REAL(dp), DIMENSION(:, :)                          :: work
      INTEGER, INTENT(IN)                                :: ldwork
      LOGICAL, INTENT(IN), OPTIONAL                      :: ordered

      LOGICAL                                            :: direct_block
      REAL(dp), ALLOCATABLE                              :: primitive(:, :)

      direct_block = iatom <= jatom
      IF (PRESENT(ordered)) direct_block = direct_block .OR. ordered

      ALLOCATE (primitive(ncoa, ncob))
      #:for component in ['cos', 'sin']
         #:if component == 'cos'
            primitive(:, :) = REAL(expab(1:ncoa, 1:ncob), KIND=dp)
         #:else
            primitive(:, :) = AIMAG(expab(1:ncoa, 1:ncob))
         #:endif
         CALL dgemm("N", "N", ncoa, nsgfb, ncob, &
                    1.0_dp, primitive(1, 1), ncoa, &
                    sphi_b(1, sgfb), ldsb, &
                    0.0_dp, work(1, 1), ldwork)

         IF (direct_block) THEN
            CALL dgemm("T", "N", nsgfa, nsgfb, ncoa, &
                       1.0_dp, sphi_a(1, sgfa), ldsa, &
                       work(1, 1), ldwork, &
                       1.0_dp, ${component}$_block(sgfa, sgfb), &
                       SIZE(${component}$_block, 1))
         ELSE
            CALL dgemm("T", "N", nsgfb, nsgfa, ncoa, &
                       1.0_dp, work(1, 1), ldwork, &
                       sphi_a(1, sgfa), ldsa, &
                       1.0_dp, ${component}$_block(sgfb, sgfa), &
                       SIZE(${component}$_block, 1))
         END IF
      #:endfor

   END SUBROUTINE contract_expikr

! **************************************************************************************************
!> \brief Evaluate primitive exp(i k.r) integrals through the complex overlap recurrence.
!> \param la_max_set highest bra shell
!> \param npgfa number of bra primitives
!> \param zeta bra Gaussian exponents
!> \param rpgfa bra primitive radii
!> \param la_min_set lowest bra shell
!> \param lb_max highest ket shell
!> \param npgfb number of ket primitives
!> \param zetb ket Gaussian exponents
!> \param rpgfb ket primitive radii
!> \param lb_min lowest ket shell
!> \param rac bra-center coordinates relative to the phase origin
!> \param rbc ket-center coordinates relative to the phase origin
!> \param kvec plane-wave vector
!> \param expab primitive integrals in the full Cartesian layout, zero below the requested shells
!> \param dexpab optional bra-center derivatives with the plane-wave operator held fixed
!> \param dexpba optional ket-center derivatives in the same bra/ket layout
!> \param momentab optional coordinate-weighted integrals relative to the phase origin
!>        Derivative and moment outputs cover the full Cartesian range.
! **************************************************************************************************
   SUBROUTINE expikr(la_max_set, npgfa, zeta, rpgfa, la_min_set, &
                     lb_max, npgfb, zetb, rpgfb, lb_min, rac, rbc, kvec, expab, dexpab, dexpba, momentab)
      INTEGER, INTENT(IN)                                :: la_max_set, npgfa
      REAL(KIND=dp), INTENT(IN)                          :: zeta(:), rpgfa(:)
      INTEGER, INTENT(IN)                                :: la_min_set, lb_max, npgfb
      REAL(KIND=dp), INTENT(IN)                          :: zetb(:), rpgfb(:)
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), INTENT(IN)                          :: rac(3), rbc(3), kvec(3)
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: expab(:, :)
      COMPLEX(KIND=dp), INTENT(INOUT), OPTIONAL, TARGET  :: dexpab(:, :, :), dexpba(:, :, :), &
                                                            momentab(:, :, :)

      COMPLEX(KIND=dp)                                   :: prefactor, zrap(3), zrbp(3)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: opr_rr(:, :, :), rr(:, :, :)
      COMPLEX(KIND=dp), POINTER                          :: opr_values(:, :, :)
      INTEGER                                            :: avec(3), bvec(3), coa, cob, ia, ib, &
                                                            idir, ipgf, jdir, jpgf, kdir, ldrr, &
                                                            lma, lmb, ma, mb, na, nb, op
      LOGICAL                                            :: derivative, right
      REAL(KIND=dp)                                      :: dab, exponent, f0, f1, f2, k2, kdp, &
                                                            origin(3), rab(3), zet, zetp

      rab = rbc - rac
      dab = NORM2(rab)
      k2 = SUM(kvec**2)
      lma = la_max_set
      lmb = lb_max
      IF (PRESENT(dexpab) .OR. PRESENT(momentab)) lma = lma + 1
      IF (PRESENT(dexpba)) lmb = lmb + 1
      ldrr = MAX(lma, lmb) + 1
      ALLOCATE (rr(0:ldrr - 1, 0:ldrr - 1, 3))
      IF (PRESENT(dexpab) .OR. PRESENT(dexpba) .OR. PRESENT(momentab)) THEN
         ALLOCATE (opr_rr(0:la_max_set, 0:lb_max, 3))
      END IF

      na = ncoset(la_max_set)
      nb = ncoset(lb_max)
      expab(1:na*npgfa, 1:nb*npgfb) = 0.0_dp
      IF (PRESENT(dexpab)) dexpab(1:na*npgfa, 1:nb*npgfb, 1:3) = 0.0_dp
      IF (PRESENT(dexpba)) dexpba(1:na*npgfa, 1:nb*npgfb, 1:3) = 0.0_dp
      IF (PRESENT(momentab)) momentab(1:na*npgfa, 1:nb*npgfb, 1:3) = 0.0_dp

      DO ipgf = 1, npgfa
         ma = (ipgf - 1)*na
         DO jpgf = 1, npgfb
            IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) CYCLE
            mb = (jpgf - 1)*nb
            zet = zeta(ipgf) + zetb(jpgf)
            zetp = 1.0_dp/zet
            f0 = (pi*zetp)**1.5_dp
            f1 = zetb(jpgf)*zetp
            f2 = 0.5_dp*zetp
            kdp = zetp*DOT_PRODUCT(kvec, zeta(ipgf)*rac + zetb(jpgf)*rbc)
            prefactor = f0*EXP(-zeta(ipgf)*f1*SUM(rab**2))*EXP(-0.25_dp*k2*zetp)* &
                        CMPLX(COS(kdp), SIN(kdp), KIND=dp)

            ! The normalized 1D tables carry no three-dimensional prefactor.
            zrap = CMPLX(f1*rab, f2*kvec, KIND=dp)
            zrbp = zrap - rab
            CALL os_rr_ovlp(zrap, lma, zrbp, lmb, zet, ldrr, rr)
            DO cob = ncoset(lb_min - 1) + 1, nb
               bvec = indco(:, cob)
               DO coa = ncoset(la_min_set - 1) + 1, na
                  avec = indco(:, coa)
                  expab(ma + coa, mb + cob) = prefactor*rr(avec(1), bvec(1), 1)* &
                                              rr(avec(2), bvec(2), 2)*rr(avec(3), bvec(3), 3)
               END DO
            END DO

            DO op = 1, 3
               SELECT CASE (op)
               CASE (1)
                  IF (.NOT. PRESENT(dexpab)) CYCLE
                  opr_values => dexpab
                  exponent = zeta(ipgf)
               CASE (2)
                  IF (.NOT. PRESENT(dexpba)) CYCLE
                  opr_values => dexpba
                  exponent = zetb(jpgf)
               CASE (3)
                  IF (.NOT. PRESENT(momentab)) CYCLE
                  opr_values => momentab
                  exponent = 0.0_dp
               END SELECT
               right = op == 2
               derivative = op /= 3
               origin = 0.0_dp
               IF (.NOT. derivative) origin = rac
               CALL apply_cartesian_1d(rr, opr_rr, la_max_set, lb_max, right, derivative, exponent, origin)
               DO cob = 1, nb
                  bvec = indco(:, cob)
                  ib = mb + cob
                  DO coa = 1, na
                     avec = indco(:, coa)
                     ia = ma + coa
                     DO idir = 1, 3
                        jdir = MOD(idir, 3) + 1
                        kdir = MOD(idir + 1, 3) + 1
                        opr_values(ia, ib, idir) = prefactor*opr_rr(avec(idir), bvec(idir), idir)* &
                                                   rr(avec(jdir), bvec(jdir), jdir)*rr(avec(kdir), bvec(kdir), kdir)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE expikr

! **************************************************************************************************
!> \brief ...
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lc_max ...
!> \param rac ...
!> \param rbc ...
!> \param mab ...
! **************************************************************************************************
   SUBROUTINE moment(la_max, npgfa, zeta, rpgfa, la_min, &
                     lb_max, npgfb, zetb, rpgfb, lc_max, rac, rbc, mab)
      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lc_max
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: mab

      INTEGER                                            :: ax, ay, az, bx, by, bz, coa, cob, i, ia, &
                                                            ib, imom, ipgf, jpgf, ldrr, lx, ly, &
                                                            lz, ma, mb, powers(0:lc_max)
      REAL(KIND=dp)                                      :: a, b, dab, f0, rab(3), rab2, rap(3), &
                                                            rbp(3), zet
      REAL(KIND=dp), ALLOCATABLE                         :: polynomial(:, :, :, :), rr(:, :, :)

      MARK_USED(la_min)
      rab = rbc - rac
      dab = NORM2(rab)
      rab2 = dab**2
      ldrr = MAX(la_max + lc_max, lb_max) + 1
      ALLOCATE (rr(0:ldrr - 1, 0:ldrr - 1, 3), polynomial(0:la_max, 0:lb_max, 0:lc_max, 3))
      powers = [(i, i=0, lc_max)]
      mab(:npgfa*ncoset(la_max), :npgfb*ncoset(lb_max), :ncoset(lc_max) - 1) = 0.0_dp
      DO ipgf = 1, npgfa
         ma = (ipgf - 1)*ncoset(la_max)
         DO jpgf = 1, npgfb
            IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) CYCLE
            mb = (jpgf - 1)*ncoset(lb_max)
            a = zeta(ipgf)
            b = zetb(jpgf)
            zet = a + b
            rap = b*rab/zet
            rbp = -a*rab/zet
            f0 = (pi/zet)**1.5_dp*EXP(-a*b*rab2/zet)
            CALL os_rr_ovlp(rap, la_max + lc_max, rbp, lb_max, zet, ldrr, rr)
            CALL cartesian_polynomial_1d(rr, polynomial, la_max, lb_max, powers, rac)

            DO cob = 1, ncoset(lb_max)
               bx = indco(1, cob)
               by = indco(2, cob)
               bz = indco(3, cob)
               ib = mb + cob
               DO coa = 1, ncoset(la_max)
                  ax = indco(1, coa)
                  ay = indco(2, coa)
                  az = indco(3, coa)
                  ia = ma + coa
                  ! Preserve the coset ordering and omission of the zeroth moment.
                  DO imom = 1, ncoset(lc_max) - 1
                     lx = indco(1, imom + 1)
                     ly = indco(2, imom + 1)
                     lz = indco(3, imom + 1)
                     mab(ia, ib, imom) = f0*polynomial(ax, bx, lx, 1)* &
                                         polynomial(ay, by, ly, 2)*polynomial(az, bz, lz, 3)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE moment

! **************************************************************************************************
!> \brief This returns the derivative of the moment integrals [a|\mu|b].
!>       By default, it differentiates the primitive on the right:
!>       [a|\mu|d/dR_bi] =  2*zetb*[a|\mu|b+1i] - Ni(b)[a|\mu|b-1i]
!>       A weighted derivative combines the left and right primitive derivatives
!>       using deltaR for the corresponding atom centers.
!>       order indicates the max order of the moment operator to be calculated
!>       1: dipole
!>       2: quadrupole
!>       ...
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param order ...
!> \param rac ...
!> \param rbc ...
!> \param difmab ...
!> \param mab_ext optional precomputed moments with zero screened primitive blocks, including raised shells
!> \param deltaR optional weights for the left and right primitive derivatives
!> \param lambda optional atom selector for the factor (iatom == lambda) - (jatom == lambda)
!> \param iatom atom associated with the left basis function
!> \param jatom atom associated with the right basis function
!> \note
! **************************************************************************************************
   SUBROUTINE diff_momop(la_max, npgfa, zeta, rpgfa, la_min, &
                         lb_max, npgfb, zetb, rpgfb, lb_min, &
                         order, rac, rbc, difmab, mab_ext, deltaR, lambda, iatom, jatom)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lb_min, order
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(OUT)  :: difmab
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: mab_ext
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: deltaR
      INTEGER, INTENT(IN), OPTIONAL                      :: lambda, iatom, jatom

      INTEGER                                            :: ider, imom, lda, lda_min, ldb, ldb_min
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: difmab_tmp
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: mab

      lda_min = MAX(0, la_min - 1)
      ldb_min = MAX(0, lb_min - 1)
      lda = ncoset(la_max)*npgfa
      ldb = ncoset(lb_max)*npgfb
      ALLOCATE (difmab_tmp(lda, ldb, 3))

      IF (PRESENT(mab_ext)) THEN
         mab => mab_ext
      ELSE
         ALLOCATE (mab(npgfa*ncoset(la_max + 1), npgfb*ncoset(lb_max + 1), &
                       ncoset(order) - 1))
         mab = 0.0_dp
!     *** Calculate the primitive overlap integrals ***
         CALL moment(la_max + 1, npgfa, zeta, rpgfa, lda_min, &
                     lb_max + 1, npgfb, zetb, rpgfb, &
                     order, rac, rbc, mab)

      END IF
      DO imom = 1, ncoset(order) - 1
         CALL apply_cartesian_operator(la_max, 0, npgfa, lb_max, 0, npgfb, mab(:, :, imom), &
                                       difmab_tmp, right=.TRUE., derivative=.TRUE., exponents=zetb)
         IF (PRESENT(deltaR)) THEN
            CPASSERT(ASSOCIATED(deltaR))
            CPASSERT(PRESENT(iatom) .AND. PRESENT(jatom))
            DO ider = 1, 3
               difmab(1:lda, 1:ldb, imom, ider) = &
                  difmab_tmp(1:lda, 1:ldb, ider)*deltaR(ider, jatom)
            END DO

            CALL apply_cartesian_operator(la_max, 0, npgfa, lb_max, 0, npgfb, mab(:, :, imom), &
                                          difmab_tmp, right=.FALSE., derivative=.TRUE., exponents=zeta)
            DO ider = 1, 3
               difmab(1:lda, 1:ldb, imom, ider) = difmab(1:lda, 1:ldb, imom, ider) &
                                                  + difmab_tmp(1:lda, 1:ldb, ider)*deltaR(ider, iatom)
            END DO
         ELSE
            difmab(1:lda, 1:ldb, imom, :) = difmab_tmp(1:lda, 1:ldb, :)
         END IF
      END DO

      IF (PRESENT(lambda)) THEN
         CPASSERT(.NOT. PRESENT(deltaR))
         CPASSERT(PRESENT(iatom) .AND. PRESENT(jatom))
         IF (iatom == lambda .AND. jatom == lambda) THEN
            difmab = 0.0_dp
         ELSE IF (iatom == lambda) THEN
            ! The right-hand derivative is selected with a positive sign.
         ELSE IF (jatom == lambda) THEN
            difmab = -difmab
         ELSE
            difmab = 0.0_dp
         END IF
      END IF

      IF (PRESENT(mab_ext)) THEN
         NULLIFY (mab)
      ELSE
         DEALLOCATE (mab)
      END IF
      DEALLOCATE (difmab_tmp)

   END SUBROUTINE diff_momop

! **************************************************************************************************
!> \brief This returns the derivative of the dipole integrals [a|x|b], with respect
!>       to the position of the primitive on the left and right, i.e.
!>       [da/dR_ai|\mu|b] =  2*zeta*[a+1i|\mu|b] - Ni(a)[a-1i|\mu|b]
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param order ...
!> \param rac ...
!> \param rbc ...
!> \param pab ...
!> \param forcea ...
!> \param forceb ...
!> \note
! **************************************************************************************************
   SUBROUTINE dipole_force(la_max, npgfa, zeta, rpgfa, la_min, &
                           lb_max, npgfb, zetb, rpgfb, lb_min, &
                           order, rac, rbc, pab, forcea, forceb)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lb_min, order
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: pab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: forcea, forceb

      INTEGER                                            :: i, imom, ipgf, j, jpgf, lda, lda_min, &
                                                            ldb, ldb_min, na, nb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: difmab, mab

      CPASSERT(order == 1)
      MARK_USED(order)

      lda_min = MAX(0, la_min - 1)
      ldb_min = MAX(0, lb_min - 1)
      lda = ncoset(la_max)*npgfa
      ldb = ncoset(lb_max)*npgfb
      ALLOCATE (difmab(lda, ldb, 3))
      ALLOCATE (mab(npgfa*ncoset(la_max + 1), npgfb*ncoset(lb_max + 1), 3))
      mab = 0.0_dp
      CALL moment(la_max + 1, npgfa, zeta, rpgfa, lda_min, &
                  lb_max + 1, npgfb, zetb, rpgfb, 1, rac, rbc, mab)

      DO imom = 1, 3
         difmab = 0.0_dp
         CALL apply_cartesian_operator(la_max, 0, npgfa, lb_max, 0, npgfb, mab(:, :, imom), &
                                       difmab, right=.TRUE., derivative=.TRUE., exponents=zetb)
         na = 0
         DO ipgf = 1, npgfa
            nb = 0
            DO jpgf = 1, npgfb
               DO j = nb + ncoset(lb_min - 1) + 1, nb + ncoset(lb_max)
                  DO i = na + ncoset(la_min - 1) + 1, na + ncoset(la_max)
                     forceb(imom, 1) = forceb(imom, 1) + pab(i, j)*difmab(i, j, 1)
                     forceb(imom, 2) = forceb(imom, 2) + pab(i, j)*difmab(i, j, 2)
                     forceb(imom, 3) = forceb(imom, 3) + pab(i, j)*difmab(i, j, 3)
                  END DO
               END DO
               nb = nb + ncoset(lb_max)
            END DO
            na = na + ncoset(la_max)
         END DO

         difmab = 0.0_dp
         CALL apply_cartesian_operator(la_max, 0, npgfa, lb_max, 0, npgfb, mab(:, :, imom), &
                                       difmab, right=.FALSE., derivative=.TRUE., exponents=zeta)
         na = 0
         DO ipgf = 1, npgfa
            nb = 0
            DO jpgf = 1, npgfb
               DO j = nb + ncoset(lb_min - 1) + 1, nb + ncoset(lb_max)
                  DO i = na + ncoset(la_min - 1) + 1, na + ncoset(la_max)
                     forcea(imom, 1) = forcea(imom, 1) + pab(i, j)*difmab(i, j, 1)
                     forcea(imom, 2) = forcea(imom, 2) + pab(i, j)*difmab(i, j, 2)
                     forcea(imom, 3) = forcea(imom, 3) + pab(i, j)*difmab(i, j, 3)
                  END DO
               END DO
               nb = nb + ncoset(lb_max)
            END DO
            na = na + ncoset(la_max)
         END DO
      END DO

      DEALLOCATE (mab, difmab)

   END SUBROUTINE dipole_force

END MODULE ai_moments
