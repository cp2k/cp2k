!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
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
                                              gto_basis_set_type,&
                                              sort_gto_basis_set
   USE bibliography,                    ONLY: Stoychev2016,&
                                              cite_reference
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: dfac,&
                                              fac,&
                                              gamma1,&
                                              pi,&
                                              rootpi
   USE orbital_pointers,                ONLY: init_orbital_pointers
   USE periodic_table,                  ONLY: get_ptable_info
   USE powell,                          ONLY: opt_state_type,&
                                              powell_optimize
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
!> \param basis_type ...
!> \param basis_sort ...
!> \date    01.11.2017
!> \author  JGH
! **************************************************************************************************
   SUBROUTINE create_ri_aux_basis_set(ri_aux_basis_set, qs_kind, basis_cntrl, basis_type, basis_sort)
      TYPE(gto_basis_set_type), POINTER                  :: ri_aux_basis_set
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind
      INTEGER, INTENT(IN)                                :: basis_cntrl
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type
      INTEGER, INTENT(IN), OPTIONAL                      :: basis_sort

      CHARACTER(len=*), PARAMETER :: routineN = 'create_ri_aux_basis_set'

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
      REAL(KIND=dp), DIMENSION(0:18)                     :: bv, bval, fv, peff, pend, pmax, pmin
      REAL(KIND=dp), DIMENSION(0:9)                      :: zeff, zmax, zmin
      REAL(KIND=dp), DIMENSION(3)                        :: amax, amin, bmin
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      !
      CALL cite_reference(Stoychev2016)
      !
      bv(0:18) = (/1.8_dp, 2.0_dp, 2.2_dp, 2.2_dp, 2.3_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, &
                   3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp/)
      fv(0:18) = (/20.0_dp, 4.0_dp, 4.0_dp, 3.5_dp, 2.5_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, &
                   2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp/)
      !
      CPASSERT(.NOT. ASSOCIATED(ri_aux_basis_set))
      NULLIFY (orb_basis_set)
      IF (.NOT. PRESENT(basis_type)) THEN
         CALL get_qs_kind(qs_kind, basis_set=orb_basis_set, basis_type="ORB")
      ELSE
         CALL get_qs_kind(qs_kind, basis_set=orb_basis_set, basis_type=basis_type)
      ENDIF
      IF (ASSOCIATED(orb_basis_set)) THEN
         CALL get_basis_keyfigures(orb_basis_set, lmax, zmin, zmax, zeff)
         !Note: RI basis coud require lmax up to 2*orb_lmax. This ensures that all orbital pointers
         !      are properly initialized before building the basis
         CALL init_orbital_pointers(2*lmax)
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
            laux = MAX(2*lval, lmax + linc)
         CASE (1)
            laux = MAX(2*lval, lmax + linc)
         CASE (2)
            laux = MAX(2*lval, lmax + linc + 1)
         CASE (3)
            laux = MAX(2*lmax, lmax + linc + 2)
         CASE DEFAULT
            CPABORT("Invalid value of control variable")
         END SELECT
         !
         DO l = 2*lmax + 1, laux
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
            xv = LOG(pend(l)/pmin(l))/LOG(bval(l)) + 1.e-10_dp
            nval(l) = MAX(CEILING(xv), 0)
         END DO
         ! first set include valence only
         nsets = 1
         ls1(1) = 0
         ls2(1) = lval
         DO l = lval + 1, laux
            IF (nval(l) < nval(lval) - 1) EXIT
            ls2(1) = l
         END DO
         ! second set up to 2*lval
         IF (laux > ls2(1)) THEN
            IF (lval == 0 .OR. 2*lval <= ls2(1) + 1) THEN
               nsets = 2
               ls1(2) = ls2(1) + 1
               ls2(2) = laux
            ELSE
               nsets = 2
               ls1(2) = ls2(1) + 1
               ls2(2) = MIN(2*lval, laux)
               lx = ls2(2)
               DO l = lx + 1, laux
                  IF (nval(l) < nval(lx) - 1) EXIT
                  ls2(2) = l
               END DO
               IF (laux > ls2(2)) THEN
                  nsets = 3
                  ls1(3) = ls2(2) + 1
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
            xv = LOG(amax(i)/amin(i))/LOG(bmin(i)) + 1.e-10_dp
            npgf(i) = MAX(CEILING(xv), 0)
         END DO
         nx = MAXVAL(npgf(1:nsets))
         ALLOCATE (zet(nx, nsets))
         zet = 0.0_dp
         nl = 0
         DO i = 1, nsets
            DO j = 1, npgf(i)
               jj = npgf(i) - j + 1
               zet(jj, i) = amin(i)*bmin(i)**(j - 1)
            END DO
            DO l = ls1(i), ls2(i)
               nl(l, i) = nval(l)
            END DO
         END DO
         bsname = TRIM(element_symbol)//"-RI-AUX-"//TRIM(orb_basis_set%name)
         !
         CALL create_aux_basis(ri_aux_basis_set, bsname, nsets, ls1, ls2, nl, npgf, zet)

         DEALLOCATE (zet)

         IF (PRESENT(basis_sort)) THEN
            CALL sort_gto_basis_set(ri_aux_basis_set, basis_sort)
         ENDIF

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

      CHARACTER(len=*), PARAMETER :: routineN = 'create_lri_aux_basis_set'

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=default_string_length)               :: bsname
      INTEGER                                            :: i, j, l, laux, linc, lm, lmax, lval, n1, &
                                                            n2, nsets, z
      INTEGER, DIMENSION(0:18)                           :: nval
      INTEGER, DIMENSION(0:9, 1:50)                      :: nl
      INTEGER, DIMENSION(1:50)                           :: ls1, ls2, npgf
      INTEGER, DIMENSION(:), POINTER                     :: econf
      REAL(KIND=dp)                                      :: xv, zval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: zet
      REAL(KIND=dp), DIMENSION(0:18)                     :: bval, peff, pend, pmax, pmin
      REAL(KIND=dp), DIMENSION(0:9)                      :: zeff, zmax, zmin
      REAL(KIND=dp), DIMENSION(4)                        :: bv, bx
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      !
      bv(1:4) = (/2.00_dp, 1.90_dp, 1.80_dp, 1.80_dp/)
      bx(1:4) = (/2.60_dp, 2.40_dp, 2.20_dp, 2.20_dp/)
      !
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
         !
         pend = 0.0_dp
         linc = 1
         IF (z > 18) linc = 2
         SELECT CASE (basis_cntrl)
         CASE (0)
            laux = MAX(2*lval, lmax + linc)
            laux = MIN(laux, 2 + linc)
         CASE (1)
            laux = MAX(2*lval, lmax + linc)
            laux = MIN(laux, 3 + linc)
         CASE (2)
            laux = MAX(2*lval, lmax + linc + 1)
            laux = MIN(laux, 4 + linc)
         CASE (3)
            laux = MAX(2*lval, lmax + linc + 1)
            laux = MIN(laux, 4 + linc)
         CASE DEFAULT
            CPABORT("Invalid value of control variable")
         END SELECT
         !
         DO l = 2*lmax + 1, laux
            pmin(l) = pmin(2*lmax)
            pmax(l) = pmax(2*lmax)
            peff(l) = peff(2*lmax)
         END DO
         !
         IF (exact_1c_terms) THEN
            DO l = 0, laux
               IF (l <= lval + 1) THEN
                  pend(l) = zmax(l) + 1.0_dp
                  bval(l) = bv(basis_cntrl + 1)
               ELSE
                  pend(l) = 2.0_dp*peff(l)
                  bval(l) = bx(basis_cntrl + 1)
               END IF
               pmin(l) = zmin(l)
               xv = LOG(pend(l)/pmin(l))/LOG(bval(l)) + 1.e-10_dp
               nval(l) = MAX(CEILING(xv), 0)
               bval(l) = (pend(l)/pmin(l))**(1._dp/nval(l))
            END DO
         ELSE
            DO l = 0, laux
               IF (l <= lval + 1) THEN
                  pend(l) = pmax(l)
                  bval(l) = bv(basis_cntrl + 1)
                  pmin(l) = zmin(l)
               ELSE
                  pend(l) = 4.0_dp*peff(l)
                  bval(l) = bx(basis_cntrl + 1)
               END IF
               xv = LOG(pend(l)/pmin(l))/LOG(bval(l)) + 1.e-10_dp
               nval(l) = MAX(CEILING(xv), 0)
               bval(l) = (pend(l)/pmin(l))**(1._dp/nval(l))
            END DO
         END IF
         !
         lm = MIN(2*lval, 3)
         n1 = MAXVAL(nval(0:lm))
         n2 = MAXVAL(nval(lm + 1:laux))
         nsets = n1 + n2
         ALLOCATE (zet(1, nsets))
         zet = 0.0_dp
         nl = 0
         j = MAXVAL(MAXLOC(nval(0:lm)))
         DO i = 1, n1
            ls1(i) = 0
            ls2(i) = lm
            npgf(i) = 1
            zet(1, i) = pmin(j)*bval(j)**(i - 1)
            DO l = 0, lm
               nl(l, i) = 1
            END DO
         END DO
         j = lm + 1
         DO i = n1 + 1, nsets
            ls1(i) = lm + 1
            ls2(i) = laux
            npgf(i) = 1
            zet(1, i) = pmin(j)*bval(j)**(i - n1 - 1)
            DO l = lm + 1, laux
               nl(l, i) = 1
            END DO
         END DO
         !
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

      CHARACTER(len=*), PARAMETER :: routineN = 'create_aux_fit_basis_set'

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=default_string_length)               :: bsname
      INTEGER                                            :: i, iset, l, laux, lmax, lval, maxpgf, &
                                                            mx, nf, np, nsets, nx, z
      INTEGER, DIMENSION(0:9)                            :: nfun, nval
      INTEGER, DIMENSION(0:9, 1:20)                      :: nl
      INTEGER, DIMENSION(1:20)                           :: lset, npgf
      INTEGER, DIMENSION(:), POINTER                     :: econf
      REAL(KIND=dp)                                      :: amet, z1, z2, zvel
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: zval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gcval, zet
      REAL(KIND=dp), DIMENSION(0:9)                      :: zeff, zmax, zmin
      REAL(KIND=dp), DIMENSION(25)                       :: afit, xval
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(opt_state_type)                               :: ostate

      !
      CPABORT("Automatic basis set generation not activated")
      !
      CPASSERT(.NOT. ASSOCIATED(aux_fit_basis))
      NULLIFY (orb_basis_set)
      CALL get_qs_kind(qs_kind, basis_set=orb_basis_set, basis_type="ORB")
      IF (ASSOCIATED(orb_basis_set)) THEN
         CALL get_basis_keyfigures(orb_basis_set, lmax, zmin, zmax, zeff)
         CALL get_qs_kind(qs_kind, zeff=zvel, elec_conf=econf, element_symbol=element_symbol)
         CALL get_ptable_info(element_symbol, ielement=z)
         lval = 0
         DO l = 0, MAXVAL(UBOUND(econf))
            IF (econf(l) > 0) lval = l
         END DO
         IF (SUM(econf) /= NINT(zvel)) THEN
            CPWARN("Valence charge and electron configuration not consistent")
         END IF
         nval = 0
         DO l = 0, lval
            nx = econf(l)
            DO
               IF (nx > 0) THEN
                  nval(l) = nval(l) + 1
                  nx = nx - 2*(2*l + 1)
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
               nfun(l) = nfun(l) + 1
            END DO
         CASE (1)
            laux = MIN(lval + 1, lmax)
            DO l = 0, lval
               nfun(l) = nfun(l) + 1
            END DO
            IF (laux > lval) nfun(laux) = 1
         CASE (2)
            laux = MIN(lval + 1, lmax)
            DO l = 0, lval
               nfun(l) = nfun(l) + 2
            END DO
            IF (laux > lval) nfun(laux) = 1
         CASE (3)
            laux = MIN(lval + 2, lmax)
            DO l = 0, lval
               nfun(l) = nfun(l) + 3
            END DO
            IF (laux > lval) nfun(lval + 1) = 2
            IF (laux > lval + 1) nfun(laux) = 1
         CASE DEFAULT
            CPABORT("Invalid value of control variable")
         END SELECT
         !
         nsets = 0
         maxpgf = 0
         DO l = 0, lval
            z1 = MAX(zmin(l), 0.10_dp + l*0.025_dp)
            z2 = zmax(l)
            mx = CEILING(LOG(z2/z1)/LOG(5.0_dp))
            IF (nval(l) > 1) THEN
               nsets = nsets + 2
               maxpgf = MAX(maxpgf, nval(l), mx, 3)
            ELSEIF (nval(l) == 1) THEN
               nsets = nsets + 1
               maxpgf = MAX(maxpgf, mx, 3)
            END IF
            DO i = nval(l) + 1, nfun(l)
               maxpgf = MAX(maxpgf, mx, 1)
               nsets = nsets + 1
            END DO
         END DO
         DO l = lval + 1, laux
            maxpgf = MAX(maxpgf, 1)
            nsets = nsets + nfun(l)
         END DO
         !
         ALLOCATE (zet(maxpgf, nsets))
         zet = 0.0_dp
         nl = 0
         iset = 0
         DO l = 0, laux
            amet = 2.50_dp
            ! optimize exponensts
            z1 = MAX(zmin(l), 0.20_dp + l*0.025_dp)
            z2 = zmax(l)
            mx = CEILING(LOG(z2/z1)/LOG(4.0_dp))
            IF (nval(l) > 1) THEN
               nx = MAX(nfun(l), mx, 3)
            ELSE IF (nval(l) == 1) THEN
               nx = MAX(nfun(l), mx, 2)
            ELSE
               nx = nfun(l)
            END IF
            CPASSERT(nx <= 25)
            ! Powell optimizer
            CALL get_basis_functions(orb_basis_set, l, np, nf, zval, gcval)
            ostate%nf = 0
            ostate%nvar = nx
            xval(:) = 5.00_dp
            xval(1) = z1
            ostate%rhoend = 1.e-8_dp
            ostate%rhobeg = 5.e-1_dp
            ostate%maxfun = 1000
            ostate%iprint = 0
            ostate%unit = 0
            ostate%state = 0
            DO
               IF (ostate%state == 2) THEN
                  afit(1) = xval(1)
                  DO i = 2, nx
                     afit(i) = afit(i - 1)*xval(i)
                  END DO
                  CALL overlap_maximum(l, np, nf, zval, gcval, nx, afit, amet, ostate%f)
                  CALL neb_potential(xval, ostate%nvar, ostate%f)
               END IF
               IF (ostate%state == -1) EXIT
               CALL powell_optimize(ostate%nvar, xval, ostate)
            END DO
            ostate%state = 8
            CALL powell_optimize(ostate%nvar, xval, ostate)
            afit(1) = xval(1)
            DO i = 2, nx
               afit(i) = afit(i - 1)*xval(i)
            END DO
            DEALLOCATE (zval, gcval)
            !
            IF (nval(l) > 1) THEN
               ! split set
               iset = iset + 1
               lset(iset) = l
               npgf(iset) = nx - 1
               nl(l, iset) = nval(l) - 1
               zet(1:nx - 1, iset) = afit(1:nx - 1)
               ! new set
               iset = iset + 1
               lset(iset) = l
               npgf(iset) = nx - 1
               nl(l, iset) = 1
               zet(1:nx - 1, iset) = afit(2:nx)
               DO i = 1, nfun(l) - 2
                  iset = iset + 1
                  lset(iset) = l
                  npgf(iset) = 1
                  zet(1, iset) = afit(nx - i + 1)
                  nl(l, iset) = 1
               END DO
            ELSEIF (nval(l) == 1) THEN
               iset = iset + 1
               lset(iset) = l
               npgf(iset) = nx
               zet(1:nx, iset) = afit(1:nx)
               nl(l, iset) = 1
               DO i = 1, nfun(l) - 1
                  iset = iset + 1
                  lset(iset) = l
                  npgf(iset) = 1
                  zet(1, iset) = afit(nx - i + 1)
                  nl(l, iset) = 1
               END DO
            ELSE
               DO i = 1, nfun(l)
                  iset = iset + 1
                  lset(iset) = l
                  npgf(iset) = 1
                  zet(1, iset) = afit(i)
                  nl(l, iset) = 1
               END DO
            END IF
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

      CHARACTER(len=*), PARAMETER :: routineN = 'get_basis_keyfigures'

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
            kval = fac(l + 1)**2*2._dp**(2*l + 1)/fac(2*l + 2)
            rexp = 0.0_dp
            rno = 0.0_dp
            DO i = 1, npgf(iset)
               gcca = gcc(i, ishell, iset)
               DO j = 1, npgf(iset)
                  zeta = zet(i, iset) + zet(j, iset)
                  gccb = gcc(j, ishell, iset)
                  rint = 0.5_dp*fac(l + 1)/zeta**(l + 2)
                  rexp = rexp + gcca*gccb*rint
                  rint = rootpi*0.5_dp**(l + 2)*dfac(2*l + 1)/zeta**(l + 1.5_dp)
                  rno = rno + gcca*gccb*rint
               END DO
            END DO
            rexp = rexp/rno
            aeff = (fac(l + 1)/dfac(2*l + 1))**2*2._dp**(2*l + 1)/(pi*rexp**2)
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

      CHARACTER(len=*), PARAMETER :: routineN = 'get_basis_products'

      INTEGER                                            :: l1, l2, la

      pmin = HUGE(0.0_dp)
      pmax = 0.0_dp
      peff = 0.0_dp

      DO l1 = 0, lmax
         DO l2 = l1, lmax
            DO la = l2 - l1, l2 + l1
               pmax(la) = MAX(pmax(la), zmax(l1) + zmax(l2))
               pmin(la) = MIN(pmin(la), zmin(l1) + zmin(l2))
               peff(la) = MAX(peff(la), zeff(l1) + zeff(l2))
            END DO
         END DO
      END DO

   END SUBROUTINE get_basis_products
! **************************************************************************************************
!> \brief ...
!> \param lm ...
!> \param npgf ...
!> \param nfun ...
!> \param zet ...
!> \param gcc ...
!> \param nfit ...
!> \param afit ...
!> \param amet ...
!> \param eval ...
! **************************************************************************************************
   SUBROUTINE overlap_maximum(lm, npgf, nfun, zet, gcc, nfit, afit, amet, eval)
      INTEGER, INTENT(IN)                                :: lm, npgf, nfun
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zet
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: gcc
      INTEGER, INTENT(IN)                                :: nfit
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: afit
      REAL(KIND=dp), INTENT(IN)                          :: amet
      REAL(KIND=dp), INTENT(OUT)                         :: eval

      CHARACTER(len=*), PARAMETER                        :: routineN = 'overlap_maximum'

      INTEGER                                            :: i, ia, ib, info
      REAL(KIND=dp)                                      :: fij, fxij, intab, p, xij
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: fx, tx, x2, xx

      ! SUM_i(fi M fi)
      fij = 0.0_dp
      DO ia = 1, npgf
         DO ib = 1, npgf
            p = zet(ia) + zet(ib) + amet
            intab = 0.5_dp/p**(lm + 1.5_dp)*gamma1(lm + 1)
            DO i = 1, nfun
               fij = fij + gcc(ia, i)*gcc(ib, i)*intab
            END DO
         END DO
      END DO

      !Integrals (fi M xj)
      ALLOCATE (fx(nfit, nfun), tx(nfit, nfun))
      fx = 0.0_dp
      DO ia = 1, npgf
         DO ib = 1, nfit
            p = zet(ia) + afit(ib) + amet
            intab = 0.5_dp/p**(lm + 1.5_dp)*gamma1(lm + 1)
            DO i = 1, nfun
               fx(ib, i) = fx(ib, i) + gcc(ia, i)*intab
            END DO
         END DO
      END DO

      !Integrals (xi M xj)
      ALLOCATE (xx(nfit, nfit), x2(nfit, nfit))
      DO ia = 1, nfit
         DO ib = 1, nfit
            p = afit(ia) + afit(ib) + amet
            xx(ia, ib) = 0.5_dp/p**(lm + 1.5_dp)*gamma1(lm + 1)
         END DO
      END DO

      !Solve for tab
      tx(1:nfit, 1:nfun) = fx(1:nfit, 1:nfun)
      x2(1:nfit, 1:nfit) = xx(1:nfit, 1:nfit)
      CALL DPOSV("U", nfit, nfun, x2, nfit, tx, nfit, info)
      IF (info == 0) THEN
         ! value t*xx*t
         xij = 0.0_dp
         DO i = 1, nfun
            xij = xij + DOT_PRODUCT(tx(:, i), MATMUL(xx, tx(:, i)))
         END DO
         ! value t*fx
         fxij = 0.0_dp
         DO i = 1, nfun
            fxij = fxij + DOT_PRODUCT(tx(:, i), fx(:, i))
         END DO
         !
         eval = fij - 2.0_dp*fxij + xij
      ELSE
         ! error in solving for max overlap
         eval = 1.0e10_dp
      END IF

      DEALLOCATE (fx, xx, x2, tx)

   END SUBROUTINE overlap_maximum
! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param n ...
!> \param eval ...
! **************************************************************************************************
   SUBROUTINE neb_potential(x, n, eval)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: x
      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), INTENT(INOUT)                       :: eval

      CHARACTER(len=*), PARAMETER                        :: routineN = 'neb_potential'

      INTEGER                                            :: i

      DO i = 2, n
         IF (x(i) < 1.5_dp) THEN
            eval = eval + 10.0_dp*(1.5_dp - x(i))**2
         END IF
      END DO

   END SUBROUTINE neb_potential
! **************************************************************************************************
!> \brief ...
!> \param basis_set ...
!> \param lin ...
!> \param np ...
!> \param nf ...
!> \param zval ...
!> \param gcval ...
! **************************************************************************************************
   SUBROUTINE get_basis_functions(basis_set, lin, np, nf, zval, gcval)
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      INTEGER, INTENT(IN)                                :: lin
      INTEGER, INTENT(OUT)                               :: np, nf
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: zval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gcval

      CHARACTER(len=*), PARAMETER :: routineN = 'get_basis_functions'

      INTEGER                                            :: iset, ishell, j1, j2, jf, jp, l, nset
      INTEGER, DIMENSION(:), POINTER                     :: lm, npgf, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: lshell
      LOGICAL                                            :: toadd
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

      np = 0
      nf = 0
      DO iset = 1, nset
         toadd = .TRUE.
         DO ishell = 1, nshell(iset)
            l = lshell(ishell, iset)
            IF (l == lin) THEN
               nf = nf + 1
               IF (toadd) THEN
                  np = np + npgf(iset)
                  toadd = .FALSE.
               END IF
            END IF
         END DO
      END DO
      ALLOCATE (zval(np), gcval(np, nf))
      zval = 0.0_dp
      gcval = 0.0_dp
      !
      jp = 0
      jf = 0
      DO iset = 1, nset
         toadd = .TRUE.
         DO ishell = 1, nshell(iset)
            l = lshell(ishell, iset)
            IF (l == lin) THEN
               jf = jf + 1
               IF (toadd) THEN
                  j1 = jp + 1
                  j2 = jp + npgf(iset)
                  zval(j1:j2) = zet(1:npgf(iset), iset)
                  jp = jp + npgf(iset)
                  toadd = .FALSE.
               END IF
               gcval(j1:j2, jf) = gcc(1:npgf(iset), ishell, iset)
            END IF
         END DO
      END DO

   END SUBROUTINE get_basis_functions

END MODULE auto_basis
