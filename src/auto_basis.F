!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Automatic generation of auxiliary basis sets of different kind
!> \author  JGH
!>
!> <b>Modification history:</b>
!> - 11.2017 creation [JGH]
! **************************************************************************************************
MODULE auto_basis
   USE aux_basis_set,                   ONLY: create_aux_basis
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: dfac,&
                                              fac,&
                                              pi,&
                                              rootpi
   USE periodic_table,                  ONLY: get_ptable_info
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'auto_basis'

   PUBLIC :: create_ri_aux_basis_set, create_aux_fit_basis_set, create_lri_aux_basis_set

CONTAINS

! **************************************************************************************************
!> \brief Create a RI_AUX basis set using some heuristics
!> \param ri_aux_basis_set ...
!> \param qs_kind ...
!> \param basis_cntrl ...
!> \date    01.11.2017
!> \author  JGH
! **************************************************************************************************
   SUBROUTINE create_ri_aux_basis_set(ri_aux_basis_set, qs_kind, basis_cntrl)
      TYPE(gto_basis_set_type), POINTER                  :: ri_aux_basis_set
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind
      INTEGER, INTENT(IN)                                :: basis_cntrl

      CHARACTER(len=*), PARAMETER :: routineN = 'create_ri_aux_basis_set', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=default_string_length)               :: bsname
      INTEGER                                            :: i, j, jj, l, laux, linc, lmax, lval, lx, &
                                                            nsets, nx, z
      INTEGER, DIMENSION(0:18)                           :: nval
      INTEGER, DIMENSION(0:9, 1:20)                      :: nl
      INTEGER, DIMENSION(1:3)                            :: ls1, ls2, npgf
      INTEGER, DIMENSION(:), POINTER                     :: econf
      REAL(KIND=dp)                                      :: xv, zval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: zet
      REAL(KIND=dp), DIMENSION(0:18) :: bv = (/1.8_dp, 2.0_dp, 2.2_dp, 2.2_dp, 2.3_dp, 3.0_dp, &
         3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, &
         3.0_dp, 3.0_dp/), bval, fv = (/20.0_dp, 4.0_dp, 4.0_dp, 3.5_dp, 2.5_dp, 2.0_dp, 2.0_dp, &
         2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, &
         2.0_dp/), peff, pend, pmax, pmin
      REAL(KIND=dp), DIMENSION(0:9)                      :: zeff, zmax, zmin
      REAL(KIND=dp), DIMENSION(3)                        :: amax, amin, bmin
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      !
      CPABORT("Automatic basis set generation not activated")
      !
      CPASSERT(.NOT. ASSOCIATED(ri_aux_basis_set))
      NULLIFY (orb_basis_set)
      CALL get_qs_kind(qs_kind, basis_set=orb_basis_set, basis_type="ORB")
      IF (ASSOCIATED(orb_basis_set)) THEN
         CALL get_basis_keyfigures(orb_basis_set, lmax, zmin, zmax, zeff)
         CALL get_basis_products(lmax, zmin, zmax, zeff, pmin, pmax, peff)
         CALL get_qs_kind(qs_kind, zeff=zval, elec_conf=econf, element_symbol=element_symbol)
         CALL get_ptable_info(element_symbol, ielement=z)
         lval = 0
         DO l = 0, MAXVAL(UBOUND(econf))
            IF (econf(l) > 0) lval = l
         END DO
         IF (SUM(econf) /= NINT(zval)) THEN
            CPWARN("Valence charge and electron configuration not consistent")
         END IF
         pend = 0.0_dp
         linc = 1
         IF (z > 18) linc = 2
         SELECT CASE (basis_cntrl)
         CASE (0)
            laux = MAX(2*lval, lmax+linc)
         CASE (1)
            laux = MAX(2*lval, lmax+linc)
         CASE (2)
            laux = MAX(2*lval, lmax+linc+1)
         CASE (3)
            laux = MAX(2*lmax, lmax+linc+2)
         CASE DEFAULT
            CPABORT("Invalid value of control variable")
         END SELECT
         !
         DO l = 2*lmax, laux
            xv = peff(2*lmax)
            pmin(l) = xv
            pmax(l) = xv
            peff(l) = xv
            pend(l) = xv
         END DO
         !
         DO l = 0, laux
            IF (l <= 2*lval) THEN
               pend(l) = MIN(fv(l)*peff(l), pmax(l))
               bval(l) = 1.8_dp
            ELSE
               pend(l) = peff(l)
               bval(l) = bv(l)
            END IF
            xv = LOG(pend(l)/pmin(l))/LOG(bval(l))
            nval(l) = MAX(NINT(xv)+1, 0)
         END DO
         ! first set include valence only
         nsets = 1
         ls1(1) = 0
         ls2(1) = lval
         DO l = lval+1, laux
            IF (nval(l) < nval(lval)-1) EXIT
            ls2(1) = l
         END DO
         ! second set up to 2*lval
         IF (laux > ls2(1)) THEN
            IF (lval == 0 .OR. 2*lval <= ls2(1)+1) THEN
               nsets = 2
               ls1(2) = ls2(1)+1
               ls2(2) = laux
            ELSE
               nsets = 2
               ls1(2) = ls2(1)+1
               ls2(2) = MIN(2*lval, laux)
               lx = ls2(2)
               DO l = lx+1, laux
                  IF (nval(l) < nval(lx)-1) EXIT
                  ls2(2) = l
               END DO
               IF (laux > ls2(2)) THEN
                  nsets = 3
                  ls1(3) = ls2(2)+1
                  ls2(3) = laux
               END IF
            END IF
         END IF
         !
         amax = 0.0
         amin = HUGE(0.0_dp)
         bmin = HUGE(0.0_dp)
         DO i = 1, nsets
            DO j = ls1(i), ls2(i)
               amax(i) = MAX(amax(i), pend(j))
               amin(i) = MIN(amin(i), pmin(j))
               bmin(i) = MIN(bmin(i), bval(j))
            END DO
            xv = LOG(amax(i)/amin(i))/LOG(bmin(i))
            npgf(i) = MAX(NINT(xv)+1, 0)
         END DO
         nx = MAXVAL(npgf(1:nsets))
         ALLOCATE (zet(nx, nsets))
         zet = 0.0_dp
         nl = 0
         DO i = 1, nsets
            DO j = 1, npgf(i)
               jj = npgf(i)-j+1
               zet(jj, i) = amin(i)*bmin(i)**(j-1)
            END DO
            DO l = ls1(i), ls2(i)
               nl(l, i) = nval(l)
            END DO
         END DO
         bsname = TRIM(element_symbol)//"-RI-AUX-"//TRIM(orb_basis_set%name)
         !
         CALL create_aux_basis(ri_aux_basis_set, bsname, nsets, ls1, ls2, nl, npgf, zet)
         !
         DEALLOCATE (zet)
      END IF

   END SUBROUTINE create_ri_aux_basis_set
! **************************************************************************************************
!> \brief Create a LRI_AUX basis set using some heuristics
!> \param lri_aux_basis_set ...
!> \param qs_kind ...
!> \param basis_cntrl ...
!> \param exact_1c_terms ...
!> \date    01.11.2017
!> \author  JGH
! **************************************************************************************************
   SUBROUTINE create_lri_aux_basis_set(lri_aux_basis_set, qs_kind, basis_cntrl, exact_1c_terms)
      TYPE(gto_basis_set_type), POINTER                  :: lri_aux_basis_set
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind
      INTEGER, INTENT(IN)                                :: basis_cntrl
      LOGICAL, INTENT(IN)                                :: exact_1c_terms

      CHARACTER(len=*), PARAMETER :: routineN = 'create_lri_aux_basis_set', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=default_string_length)               :: bsname
      INTEGER                                            :: i, j, jj, l, laux, linc, lmax, lval, lx, &
                                                            nsets, nx, z
      INTEGER, DIMENSION(0:18)                           :: nval
      INTEGER, DIMENSION(0:9, 1:20)                      :: nl
      INTEGER, DIMENSION(1:3)                            :: ls1, ls2, npgf
      INTEGER, DIMENSION(:), POINTER                     :: econf
      REAL(KIND=dp)                                      :: xv, zval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: zet
      REAL(KIND=dp), DIMENSION(0:18) :: bv = (/1.8_dp, 2.0_dp, 2.2_dp, 2.2_dp, 2.3_dp, 3.0_dp, &
         3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, &
         3.0_dp, 3.0_dp/), bval, fv = (/20.0_dp, 4.0_dp, 4.0_dp, 3.5_dp, 2.5_dp, 2.0_dp, 2.0_dp, &
         2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, &
         2.0_dp/), peff, pend, pmax, pmin
      REAL(KIND=dp), DIMENSION(0:9)                      :: zeff, zmax, zmin
      REAL(KIND=dp), DIMENSION(3)                        :: amax, amin, bmin
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      !
      CPABORT("Automatic basis set generation not activated")
      !
      CPASSERT(.NOT. exact_1c_terms)
      CPASSERT(.NOT. ASSOCIATED(lri_aux_basis_set))
      NULLIFY (orb_basis_set)
      CALL get_qs_kind(qs_kind, basis_set=orb_basis_set, basis_type="ORB")
      IF (ASSOCIATED(orb_basis_set)) THEN
         CALL get_basis_keyfigures(orb_basis_set, lmax, zmin, zmax, zeff)
         CALL get_basis_products(lmax, zmin, zmax, zeff, pmin, pmax, peff)
         CALL get_qs_kind(qs_kind, zeff=zval, elec_conf=econf, element_symbol=element_symbol)
         CALL get_ptable_info(element_symbol, ielement=z)
         lval = 0
         DO l = 0, MAXVAL(UBOUND(econf))
            IF (econf(l) > 0) lval = l
         END DO
         IF (SUM(econf) /= NINT(zval)) THEN
            CPWARN("Valence charge and electron configuration not consistent")
         END IF
         pend = 0.0_dp
         linc = 1
         IF (z > 18) linc = 2
         SELECT CASE (basis_cntrl)
         CASE (0)
            laux = MAX(2*lval, lmax+linc)
         CASE (1)
            laux = MAX(2*lval, lmax+linc)
            laux = MAX(laux, 4)
         CASE (2)
            laux = MAX(2*lval, lmax+linc+1)
            laux = MAX(laux, 4)
         CASE (3)
            laux = MAX(2*lmax, lmax+linc+2)
            laux = MAX(laux, 5)
         CASE DEFAULT
            CPABORT("Invalid value of control variable")
         END SELECT
         !
         DO l = 2*lmax, laux
            xv = peff(2*lmax)
            pmin(l) = xv
            pmax(l) = xv
            peff(l) = xv
            pend(l) = xv
         END DO
         !
         DO l = 0, laux
            IF (l <= 2*lval) THEN
               pend(l) = MIN(fv(l)*peff(l), pmax(l))
               bval(l) = 1.8_dp
            ELSE
               pend(l) = peff(l)
               bval(l) = bv(l)
            END IF
            xv = LOG(pend(l)/pmin(l))/LOG(bval(l))
            nval(l) = MAX(NINT(xv)+1, 0)
         END DO
         ! first set include valence only
         nsets = 1
         ls1(1) = 0
         ls2(1) = lval
         DO l = lval+1, laux
            IF (nval(l) < nval(lval)-1) EXIT
            ls2(1) = l
         END DO
         ! second set up to 2*lval
         IF (laux > ls2(1)) THEN
            IF (lval == 0 .OR. 2*lval <= ls2(1)+1) THEN
               nsets = 2
               ls1(2) = ls2(1)+1
               ls2(2) = laux
            ELSE
               nsets = 2
               ls1(2) = ls2(1)+1
               ls2(2) = MIN(2*lval, laux)
               lx = ls2(2)
               DO l = lx+1, laux
                  IF (nval(l) < nval(lx)-1) EXIT
                  ls2(2) = l
               END DO
               IF (laux > ls2(2)) THEN
                  nsets = 3
                  ls1(3) = ls2(2)+1
                  ls2(3) = laux
               END IF
            END IF
         END IF
         !
         amax = 0.0
         amin = HUGE(0.0_dp)
         bmin = HUGE(0.0_dp)
         DO i = 1, nsets
            DO j = ls1(i), ls2(i)
               amax(i) = MAX(amax(i), pend(j))
               amin(i) = MIN(amin(i), pmin(j))
               bmin(i) = MIN(bmin(i), bval(j))
            END DO
            xv = LOG(amax(i)/amin(i))/LOG(bmin(i))
            npgf(i) = MAX(NINT(xv)+1, 0)
         END DO
         nx = MAXVAL(npgf(1:nsets))
         ALLOCATE (zet(nx, nsets))
         zet = 0.0_dp
         nl = 0
         DO i = 1, nsets
            DO j = 1, npgf(i)
               jj = npgf(i)-j+1
               zet(jj, i) = amin(i)*bmin(i)**(j-1)
            END DO
            DO l = ls1(i), ls2(i)
               nl(l, i) = nval(l)
            END DO
         END DO
         bsname = TRIM(element_symbol)//"-LRI-AUX-"//TRIM(orb_basis_set%name)
         !
         CALL create_aux_basis(lri_aux_basis_set, bsname, nsets, ls1, ls2, nl, npgf, zet)
         !
         DEALLOCATE (zet)
      END IF

   END SUBROUTINE create_lri_aux_basis_set
! **************************************************************************************************
!> \brief Create a AUX_FIT basis set using some heuristics
!> \param aux_fit_basis ...
!> \param qs_kind ...
!> \param basis_cntrl ...
!> \date    01.11.2017
!> \author  JGH
! **************************************************************************************************
   SUBROUTINE create_aux_fit_basis_set(aux_fit_basis, qs_kind, basis_cntrl)
      TYPE(gto_basis_set_type), POINTER                  :: aux_fit_basis
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind
      INTEGER, INTENT(IN)                                :: basis_cntrl

      CHARACTER(len=*), PARAMETER :: routineN = 'create_aux_fit_basis_set', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=default_string_length)               :: bsname
      INTEGER                                            :: i, ii, iset, l, laux, lmax, lval, &
                                                            maxpgf, np, nsets, nx, z
      INTEGER, DIMENSION(0:9)                            :: nfun, nval
      INTEGER, DIMENSION(0:9, 1:20)                      :: nl
      INTEGER, DIMENSION(1:20)                           :: lset, npgf
      INTEGER, DIMENSION(:), POINTER                     :: econf
      REAL(KIND=dp)                                      :: acmin, amax, amin, bmin, zval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: zet
      REAL(KIND=dp), DIMENSION(0:9)                      :: zeff, zmax, zmin
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      !
      CPABORT("Automatic basis set generation not activated")
      !
      CPASSERT(.NOT. ASSOCIATED(aux_fit_basis))
      NULLIFY (orb_basis_set)
      CALL get_qs_kind(qs_kind, basis_set=orb_basis_set, basis_type="ORB")
      IF (ASSOCIATED(orb_basis_set)) THEN
         CALL get_basis_keyfigures(orb_basis_set, lmax, zmin, zmax, zeff)
         CALL get_qs_kind(qs_kind, zeff=zval, elec_conf=econf, element_symbol=element_symbol)
         CALL get_ptable_info(element_symbol, ielement=z)
         lval = 0
         DO l = 0, MAXVAL(UBOUND(econf))
            IF (econf(l) > 0) lval = l
         END DO
         IF (SUM(econf) /= NINT(zval)) THEN
            CPWARN("Valence charge and electron configuration not consistent")
         END IF
         nval = 0
         DO l = 0, lval
            nx = econf(l)
            DO
               IF (nx > 0) THEN
                  nval(l) = nval(l)+1
                  nx = nx-2*(2*l+1)
               ELSE
                  EXIT
               END IF
            END DO
         END DO
         nfun = nval
         SELECT CASE (basis_cntrl)
         CASE (0)
            laux = lval
            DO l = 0, lval
               nfun(l) = nfun(l)+1
            END DO
         CASE (1)
            laux = MIN(lval+1, lmax)
            DO l = 0, lval
               nfun(l) = nfun(l)+1
            END DO
            IF (laux > lval) nfun(laux) = 1
         CASE (2)
            laux = MIN(lval+1, lmax)
            DO l = 0, lval
               nfun(l) = nfun(l)+2
            END DO
            IF (laux > lval) nfun(laux) = 1
         CASE (3)
            laux = MIN(lval+2, lmax)
            DO l = 0, lval
               nfun(l) = nfun(l)+3
            END DO
            IF (laux > lval) nfun(lval+1) = 2
            IF (laux > lval+1) nfun(laux) = 1
         CASE DEFAULT
            CPABORT("Invalid value of control variable")
         END SELECT
         nsets = 0
         maxpgf = 0
         DO l = 0, lval
            IF (nval(l) > 1) THEN
               nsets = nsets+2
               maxpgf = MAX(maxpgf, nval(l), 3)
            ELSEIF (nval(l) == 1) THEN
               nsets = nsets+1
               maxpgf = MAX(maxpgf, 3)
            END IF
            DO i = nval(l)+1, nfun(l)
               maxpgf = MAX(maxpgf, 1)
               nsets = nsets+1
            END DO
         END DO
         DO l = lval+1, laux
            maxpgf = MAX(maxpgf, 1)
            nsets = nsets+nfun(l)
         END DO
         !
         ALLOCATE (zet(maxpgf, nsets))
         zet = 0.0_dp
         nl = 0
         iset = 0
         DO l = 0, lval
            IF (nval(l) > 1) THEN
               iset = iset+1
               lset(iset) = l
               nx = MAX(nval(l), 3)
               npgf(iset) = nx
               nl(l, iset) = nval(l)-1
               np = nx+2
               amax = zmax(l)
               amin = MIN(10._dp*zmin(l), zeff(l))
               bmin = (amax/amin)**(1._dp/(np-1))
               bmin = MIN(bmin, 8.0_dp)
               acmin = amax/(bmin)**(np-3)
               DO i = 1, np
                  zet(i, iset) = amax/bmin**(np-1)
               END DO
               ! new set
               iset = iset+1
               lset(iset) = l
               nx = 2
               npgf(iset) = nx
               nl(l, iset) = 1
               np = nx+(nfun(l)-nval(l))
               amax = acmin
               amin = MAX(zmin(l), 0.10_dp+l*0.025_dp)
               bmin = (amax/amin)**(1.0_dp/(np-1))
               bmin = MIN(bmin, 6.0_dp)
               amax = amin*bmin**(np-1)
               zet(1, iset) = amax
               zet(2, iset) = amax/bmin
            ELSEIF (nval(l) == 1) THEN
               iset = iset+1
               lset(iset) = l
               nx = 3
               npgf(iset) = nx
               nl(l, iset) = 1
               np = nx+(nfun(l)-nval(l))
               amax = zmax(l)
               amin = MAX(zmin(l), 0.10_dp+l*0.025_dp)
               amin = MIN(amin, zeff(l))
               bmin = (amax/amin)**(1.0_dp/(np-1))
               bmin = MIN(bmin, 6.0_dp)
               amax = amin*bmin**(np-1)
               zet(1, iset) = amax
               zet(2, iset) = amax/bmin
               zet(3, iset) = amax/bmin/bmin
            ELSE
               nx = 0
               np = nx+(nfun(l)-nval(l))
               IF (np == 1) THEN
                  amin = zeff(l)
                  bmin = 2.00_dp
               ELSEIF (np == 2) THEN
                  amin = zeff(l)
                  bmin = zmax(l)/amin
                  bmin = MIN(4.0_dp, bmin)
               ELSE
                  bmin = (zmax(l)/zeff(l))**(1._dp/(np-2))
                  bmin = MIN(4.0_dp, bmin)
                  amin = zeff(l)/bmin
               END IF
            END IF
            DO i = nfun(l), nval(l)+1, -1
               iset = iset+1
               npgf(iset) = 1
               lset(iset) = l
               nl(l, iset) = 1
               ii = i-nval(l)-1
               zet(1, iset) = amin*bmin**ii
            END DO
         END DO
         DO l = lval+1, laux
            acmin = MAX(zeff(l), 0.10_dp+l*0.025_dp)
            np = nfun(l)
            bmin = 2.25_dp
            ii = -(np-1)/2
            amin = acmin*bmin**ii
            IF (amin < zmin(l)) THEN
               CPASSERT(ii /= 0)
               bmin = (acmin/zmin(l))**(1._dp/REAL(-ii, dp))
            END IF
            ii = np/2
            amax = acmin*bmin**ii
            IF (amax > zmax(l)) THEN
               CPASSERT(ii /= 0)
               bmin = (zmax(l)/acmin)**(1._dp/REAL(ii, dp))
            END IF
            DO i = nfun(l), 1, -1
               iset = iset+1
               npgf(iset) = 1
               lset(iset) = l
               nl(l, iset) = 1
               ii = -(np-1)/2+i-1
               zet(1, iset) = acmin*bmin**ii
            END DO
         END DO
         !
         bsname = TRIM(element_symbol)//"-AUX-FIT-"//TRIM(orb_basis_set%name)
         !
         CALL create_aux_basis(aux_fit_basis, bsname, nsets, lset, lset, nl, npgf, zet)
         !
         DEALLOCATE (zet)
         !
      END IF

   END SUBROUTINE create_aux_fit_basis_set
! **************************************************************************************************
!> \brief ...
!> \param basis_set ...
!> \param lmax ...
!> \param zmin ...
!> \param zmax ...
!> \param zeff ...
! **************************************************************************************************
   SUBROUTINE get_basis_keyfigures(basis_set, lmax, zmin, zmax, zeff)
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      INTEGER, INTENT(OUT)                               :: lmax
      REAL(KIND=dp), DIMENSION(0:9), INTENT(OUT)         :: zmin, zmax, zeff

      CHARACTER(len=*), PARAMETER :: routineN = 'get_basis_keyfigures', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, ipgf, iset, ishell, j, l, nset
      INTEGER, DIMENSION(:), POINTER                     :: lm, npgf, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: lshell
      REAL(KIND=dp)                                      :: aeff, gcca, gccb, kval, rexp, rint, rno, &
                                                            zeta
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: gcc

      CALL get_gto_basis_set(gto_basis_set=basis_set, &
                             nset=nset, &
                             nshell=nshell, &
                             npgf=npgf, &
                             l=lshell, &
                             lmax=lm, &
                             zet=zet, &
                             gcc=gcc)

      lmax = MAXVAL(lm)
      CPASSERT(lmax <= 9)

      zmax = 0.0_dp
      zmin = HUGE(0.0_dp)
      zeff = 0.0_dp

      DO iset = 1, nset
         ! zmin zmax
         DO ipgf = 1, npgf(iset)
            DO ishell = 1, nshell(iset)
               l = lshell(ishell, iset)
               zeta = zet(ipgf, iset)
               zmax(l) = MAX(zmax(l), zeta)
               zmin(l) = MIN(zmin(l), zeta)
            END DO
         END DO
         ! zeff
         DO ishell = 1, nshell(iset)
            l = lshell(ishell, iset)
            kval = fac(l+1)**2*2._dp**(2*l+1)/fac(2*l+2)
            rexp = 0.0_dp
            rno = 0.0_dp
            DO i = 1, npgf(iset)
               gcca = gcc(i, ishell, iset)
               DO j = 1, npgf(iset)
                  zeta = zet(i, iset)+zet(j, iset)
                  gccb = gcc(j, ishell, iset)
                  rint = 0.5_dp*fac(l+1)/zeta**(l+2)
                  rexp = rexp+gcca*gccb*rint
                  rint = rootpi*0.5_dp**(l+2)*dfac(2*l+1)/zeta**(l+1.5_dp)
                  rno = rno+gcca*gccb*rint
               END DO
            END DO
            rexp = rexp/rno
            aeff = (fac(l+1)/dfac(2*l+1))**2*2._dp**(2*l+1)/(pi*rexp**2)
            zeff(l) = MAX(zeff(l), aeff)
         END DO
      END DO

   END SUBROUTINE get_basis_keyfigures

! **************************************************************************************************
!> \brief ...
!> \param lmax ...
!> \param zmin ...
!> \param zmax ...
!> \param zeff ...
!> \param pmin ...
!> \param pmax ...
!> \param peff ...
! **************************************************************************************************
   SUBROUTINE get_basis_products(lmax, zmin, zmax, zeff, pmin, pmax, peff)
      INTEGER, INTENT(IN)                                :: lmax
      REAL(KIND=dp), DIMENSION(0:9), INTENT(IN)          :: zmin, zmax, zeff
      REAL(KIND=dp), DIMENSION(0:18), INTENT(OUT)        :: pmin, pmax, peff

      CHARACTER(len=*), PARAMETER :: routineN = 'get_basis_products', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: l1, l2, la

      pmin = HUGE(0.0_dp)
      pmax = 0.0_dp
      peff = 0.0_dp

      DO l1 = 0, lmax
         DO l2 = l1, lmax
            DO la = l2-l1, l2+l1
               pmax(la) = MAX(pmax(la), zmax(l1)+zmax(l2))
               pmin(la) = MIN(pmin(la), zmin(l1)+zmin(l2))
               peff(la) = MAX(peff(la), zeff(l1)+zeff(l2))
            END DO
         END DO
      END DO

   END SUBROUTINE get_basis_products
! **************************************************************************************************

END MODULE auto_basis
