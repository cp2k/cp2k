!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE atom_sgp
   USE ai_onecenter,                    ONLY: sg_overlap
   USE atom_types,                      ONLY: &
        atom_basis_gridrep, atom_basis_type, atom_ecppot_type, atom_p_type, atom_type, &
        create_opmat, opmat_type, release_atom_basis, release_opmat
   USE atom_upf,                        ONLY: atom_upfpot_type
   USE atom_utils,                      ONLY: integrate_grid,&
                                              numpot_matrix
   USE erf_fn,                          ONLY: erf
   USE input_constants,                 ONLY: ecp_pseudo,&
                                              gaussian,&
                                              gth_pseudo,&
                                              no_pseudo,&
                                              upf_pseudo
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: dfac,&
                                              rootpi
   USE mathlib,                         ONLY: diamat_all,&
                                              get_pseudo_inverse_diag
#include "./base/base_uses.f90"

   IMPLICIT NONE

   TYPE sgp_potential
      LOGICAL                                  :: has_nonlocal
      INTEGER                                  :: n_nonlocal
      INTEGER                                  :: lmax
      LOGICAL, DIMENSION(0:5)                  :: is_nonlocal = .FALSE.
      REAL(KIND=dp), DIMENSION(:), POINTER     :: a_nonlocal => Null()
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_nonlocal => Null()
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER  :: c_nonlocal => Null()
      LOGICAL                                  :: has_local
      INTEGER                                  :: n_local
      REAL(KIND=dp)                            :: zval
      REAL(KIND=dp)                            :: ac_local
      REAL(KIND=dp), DIMENSION(:), POINTER     :: a_local => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER     :: c_local => Null()
      LOGICAL                                  :: has_nlcc
      INTEGER                                  :: n_nlcc
      REAL(KIND=dp), DIMENSION(:), POINTER     :: a_nlcc => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER     :: c_nlcc => Null()
   END TYPE

   PRIVATE
   PUBLIC  :: atom_sgp_construction

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_sgp'

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param atom_info ...
!> \param input_section ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE atom_sgp_construction(atom_info, input_section, iw)

      TYPE(atom_p_type), DIMENSION(:, :), POINTER        :: atom_info
      TYPE(section_vals_type), POINTER                   :: input_section
      INTEGER, INTENT(IN)                                :: iw

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_sgp_construction', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, n, ppot_type
      LOGICAL                                            :: do_transform, explicit, is_ecp, is_upf
      REAL(KIND=dp)                                      :: errcc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cgauss
      TYPE(atom_ecppot_type), POINTER                    :: ecp_pot
      TYPE(atom_type), POINTER                           :: atom_ref
      TYPE(atom_upfpot_type), POINTER                    :: upf_pot
      TYPE(opmat_type), POINTER                          :: core, hnl, score, shnl
      TYPE(sgp_potential)                                :: sgp_pot

      CALL section_vals_get(input_section, explicit=explicit)

      IF (iw > 0) WRITE (iw, '(/," ",79("*"),/,T24,A,/," ",79("*"))') "SEPARABLE GAUSSIAN PSEUDOPOTENTIAL"

      atom_ref => atom_info(1, 1)%atom

      ppot_type = atom_ref%potential%ppot_type
      SELECT CASE (ppot_type)
      CASE (gth_pseudo)
         IF (iw > 0) WRITE (iw, '(" GTH Pseudopotential is already in SGP form. ")')
         do_transform = .FALSE.
      CASE (ecp_pseudo)
         do_transform = .TRUE.
         is_ecp = .TRUE.
         is_upf = .FALSE.
         ecp_pot => atom_ref%potential%ecp_pot
      CASE (upf_pseudo)
         do_transform = .TRUE.
         is_ecp = .FALSE.
         is_upf = .TRUE.
         upf_pot => atom_ref%potential%upf_pot
      CASE (no_pseudo)
         IF (iw > 0) WRITE (iw, '(" No Pseudopotential available for transformation. ")')
         do_transform = .FALSE.
      CASE DEFAULT
         CPABORT("")
      END SELECT

      ! generate the transformed potentials
      IF (do_transform) THEN
         IF (is_ecp) THEN
            CALL ecp_sgp_constr(ecp_pot, sgp_pot, atom_ref%basis)
         ELSEIF (is_upf) THEN
            CALL upf_sgp_constr(upf_pot, sgp_pot, atom_ref%basis)
         ELSE
            CPABORT("")
         END IF
      END IF

      ! Check the result
      IF (do_transform) THEN
         NULLIFY (core, hnl)
         CALL create_opmat(core, atom_ref%basis%nbas)
         CALL create_opmat(hnl, atom_ref%basis%nbas, 5)
         NULLIFY (score, shnl)
         CALL create_opmat(score, atom_ref%basis%nbas)
         CALL create_opmat(shnl, atom_ref%basis%nbas, 5)
         !
         IF (is_ecp) THEN
            CALL ecpints(hnl%op, atom_ref%basis, ecp_pot)
         ELSEIF (is_upf) THEN
            CALL upfints(core%op, hnl%op, atom_ref%basis, upf_pot, sgp_pot%ac_local)
         ELSE
            CPABORT("")
         END IF
         !
         CALL sgpints(score%op, shnl%op, atom_ref%basis, sgp_pot)
         !
         IF (sgp_pot%has_local) THEN
            errcc = MAXVAL(ABS(core%op-score%op))
            IF (iw > 0) THEN
               WRITE (iw, '(" Local part of pseudopotential")')
               WRITE (iw, '(" Number of basis functions ",T77,i4)') sgp_pot%n_local
               WRITE (iw, '(" Max. abs. error of matrix element ",T69,f12.8)') errcc
            END IF
         END IF
         IF (sgp_pot%has_nonlocal) THEN
            errcc = MAXVAL(ABS(hnl%op-shnl%op))
            IF (iw > 0) THEN
               WRITE (iw, '(" Nonlocal part of pseudopotential")')
               WRITE (iw, '(" Max. l-quantum number",T77,i4)') sgp_pot%lmax
               WRITE (iw, '(" Number of projector basis functions ",T77,i4)') sgp_pot%n_nonlocal
               WRITE (iw, '(" Max. abs. error of matrix element ",T69,f12.8)') errcc
            END IF
         END IF
         IF (sgp_pot%has_nlcc) THEN
            IF (is_upf) THEN
               n = SIZE(upf_pot%r)
               ALLOCATE (cgauss(n))
               cgauss = 0.0_dp
               DO i = 1, sgp_pot%n_nlcc
                  cgauss(:) = cgauss(:)+sgp_pot%c_nlcc(i)*EXP(-sgp_pot%a_nlcc(i)*upf_pot%r(:)**2)
               END DO
               errcc = SUM((cgauss(:)-upf_pot%rho_nlcc(:))**2*upf_pot%r(:)**2*upf_pot%rab(:))
               errcc = SQRT(errcc/REAL(n, KIND=dp))
               DEALLOCATE (cgauss)
            ELSE
               CPABORT("")
            END IF
            IF (iw > 0) THEN
               WRITE (iw, '(" Non-linear core correction: core density")')
               WRITE (iw, '(" Number of basis functions ",T77,i4)') sgp_pot%n_nlcc
               WRITE (iw, '(" RMS error of core density ",T69,f12.8)') errcc
            END IF
         END IF
         !
         CALL release_opmat(score)
         CALL release_opmat(shnl)
         CALL release_opmat(core)
         CALL release_opmat(hnl)
      END IF

      CALL sgp_release(sgp_pot)

      IF (iw > 0) WRITE (iw, '(" ",79("*"))')

   END SUBROUTINE atom_sgp_construction

! **************************************************************************************************
!> \brief ...
!> \param ecp_pot ...
!> \param sgp_pot ...
!> \param basis ...
! **************************************************************************************************
   SUBROUTINE ecp_sgp_constr(ecp_pot, sgp_pot, basis)

      TYPE(atom_ecppot_type), POINTER                    :: ecp_pot
      TYPE(sgp_potential)                                :: sgp_pot
      TYPE(atom_basis_type)                              :: basis

      CHARACTER(len=*), PARAMETER :: routineN = 'ecp_sgp_constr', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, ia, ir, j, k, l, n, na, nl, nr
      REAL(KIND=dp)                                      :: alpha, eee, ei
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: al, cl, cpot, pgauss
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cmat, qmat, score, sinv, smat, tmat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rad

      sgp_pot%has_nlcc = .FALSE.
      sgp_pot%n_nlcc = 0
      sgp_pot%has_local = .FALSE.
      sgp_pot%n_local = 0

      ! transform semilocal potential into a separable local form
      sgp_pot%has_nonlocal = .TRUE.
      !
      nl = 12
      !
      sgp_pot%n_nonlocal = nl
      sgp_pot%lmax = ecp_pot%lmax
      ALLOCATE (sgp_pot%a_nonlocal(nl))
      ALLOCATE (sgp_pot%h_nonlocal(nl, 0:ecp_pot%lmax))
      ALLOCATE (sgp_pot%c_nonlocal(nl, nl, 0:ecp_pot%lmax))
      !
      ALLOCATE (al(nl), cl(nl))
      ALLOCATE (smat(nl, nl), sinv(nl, nl))
      ALLOCATE (tmat(nl, nl), cmat(nl, nl))
      al = 0.0_dp
      DO ir = 1, nl
         al(ir) = 80.0_dp*0.60_dp**(ir-1)
      END DO
      !
      sgp_pot%a_nonlocal(1:nl) = al(1:nl)
      !
      nr = basis%grid%nr
      rad => basis%grid%rad
      ALLOCATE (cpot(nr), pgauss(nr))
      DO l = 0, ecp_pot%lmax
         na = basis%nbas(l)
         ALLOCATE (score(na, na), qmat(na, nl))
         cpot = 0._dp
         DO k = 1, ecp_pot%npot(l)
            n = ecp_pot%nrpot(k, l)
            alpha = ecp_pot%bpot(k, l)
            cpot(:) = cpot+ecp_pot%apot(k, l)*rad**(n-2)*EXP(-alpha*rad**2)
         END DO
         DO i = 1, na
            DO j = i, na
               score(i, j) = integrate_grid(cpot, basis%bf(:, i, l), basis%bf(:, j, l), basis%grid)
               score(j, i) = score(i, j)
            END DO
         END DO
         ! overlap basis with projectors
         DO i = 1, nl
            pgauss(:) = EXP(-al(i)*rad(:)**2)*rad(:)**l
            eee = rootpi/(2._dp**(l+2)*dfac(2*l+1))/(2._dp*al(i))**(l+1.5_dp)
            pgauss(:) = pgauss(:)/SQRT(eee)
            DO ia = 1, na
               qmat(ia, i) = SUM(basis%bf(:, ia, l)*pgauss(:)*basis%grid%wr(:))
            END DO
         END DO
         ! tmat = qmat * score * qmat
         tmat(1:nl, 1:nl) = MATMUL(TRANSPOSE(qmat(1:na, 1:nl)), MATMUL(score(1:na, 1:na), qmat(1:na, 1:nl)))
         smat(1:nl, 1:nl) = MATMUL(TRANSPOSE(qmat(1:na, 1:nl)), qmat(1:na, 1:nl))
         CALL get_pseudo_inverse_diag(smat(1:nl, 1:nl), sinv(1:nl, 1:nl), 1.e-10_dp)
         cmat(1:nl, 1:nl) = MATMUL(sinv(1:nl, 1:nl), MATMUL(tmat(1:nl, 1:nl), sinv(1:nl, 1:nl)))
         cmat(1:nl, 1:nl) = (cmat(1:nl, 1:nl)+TRANSPOSE(cmat(1:nl, 1:nl)))*0.5_dp
         CALL diamat_all(cmat(1:nl, 1:nl), cl(1:nl), .TRUE.)
         !
         ! get back unnormalized Gaussians
         DO i = 1, nl
            ei = rootpi/(2._dp**(l+2)*dfac(2*l+1))/(2._dp*al(i))**(l+1.5_dp)
            cmat(i, 1:nl) = cmat(i, 1:nl)/SQRT(ei)
         END DO
         sgp_pot%h_nonlocal(1:nl, l) = cl(1:nl)
         sgp_pot%c_nonlocal(1:nl, 1:nl, l) = cmat(1:nl, 1:nl)
         sgp_pot%is_nonlocal(l) = .TRUE.
         !
         DEALLOCATE (score, qmat)
      END DO
      DEALLOCATE (cpot, pgauss)
      DEALLOCATE (al, cl, smat, sinv, tmat, cmat)

   END SUBROUTINE ecp_sgp_constr

! **************************************************************************************************
!> \brief ...
!> \param upf_pot ...
!> \param sgp_pot ...
!> \param basis ...
! **************************************************************************************************
   SUBROUTINE upf_sgp_constr(upf_pot, sgp_pot, basis)

      TYPE(atom_upfpot_type), POINTER                    :: upf_pot
      TYPE(sgp_potential)                                :: sgp_pot
      TYPE(atom_basis_type)                              :: basis

      CHARACTER(len=*), PARAMETER :: routineN = 'upf_sgp_constr', routineP = moduleN//':'//routineN

      CHARACTER(len=4)                                   :: ptype
      INTEGER                                            :: ia, ib, ipa, ipb, ir, ix, la, lb, na, &
                                                            nl, np, nr
      LOGICAL                                            :: nl_trans
      REAL(KIND=dp)                                      :: cpa, cpb, eee, ei, errcc, errloc, rc, &
                                                            zval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: al, ccharge, cgauss, cl, pgauss, pproa, &
                                                            pprob, tv, vgauss, vloc, ww
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cmat, qmat, rmat, score, sinv, smat, tmat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: gmat
      TYPE(atom_basis_type)                              :: gbasis

      IF (upf_pot%is_ultrasoft .OR. upf_pot%is_paw .OR. upf_pot%is_coulomb) THEN
         sgp_pot%has_nonlocal = .FALSE.
         sgp_pot%n_nonlocal = 0
         sgp_pot%has_local = .FALSE.
         sgp_pot%n_local = 0
         sgp_pot%has_nlcc = .FALSE.
         sgp_pot%n_nlcc = 0
         RETURN
      END IF

      ! radial grid
      nr = SIZE(upf_pot%r)
      ! weights for integration
      ALLOCATE (ww(nr))
      ww(:) = upf_pot%r(:)**2*upf_pot%rab(:)

      ! start with local potential
      ! fit local potential to Gaussian form
      ALLOCATE (vloc(nr), vgauss(nr))
      ! smearing of core charge
      zval = upf_pot%zion
      sgp_pot%ac_local = 0.5_dp
      sgp_pot%zval = zval
      DO ir = 1, nr
         IF (upf_pot%r(ir) < 1.e-12_dp) THEN
            vgauss(ir) = -2.0_dp*zval/rootpi/sgp_pot%ac_local
         ELSE
            rc = upf_pot%r(ir)/sgp_pot%ac_local
            vgauss(ir) = -zval/upf_pot%r(ir)*erf(rc)
         END IF
      END DO
      vloc(:) = upf_pot%vlocal(:)-vgauss(:)
      !
      nl = 8
      ALLOCATE (al(nl), cl(nl), tv(nl))
      ALLOCATE (smat(nl, nl), sinv(nl, nl))
      al = 0.0_dp
      cl = 0.0_dp
      DO ir = 1, nl
         al(ir) = 10.0_dp*0.60_dp**(ir-1)
      END DO
      !
      smat = 0.0_dp
      tv = 0.0_dp
      CALL atom_basis_gridrep(basis, gbasis, upf_pot%r, upf_pot%rab)
      DO la = 0, UBOUND(basis%nbas, 1)
         na = basis%nbas(la)
         ALLOCATE (rmat(na, na), gmat(na, na, nl))
         DO ia = 1, na
            DO ib = ia, na
               rmat(ia, ib) = SUM(gbasis%bf(:, ia, la)*gbasis%bf(:, ib, la)*vloc(:)*ww(:))
               rmat(ib, ia) = rmat(ia, ib)
            END DO
         END DO
         DO ir = 1, nl
            vgauss(:) = EXP(-al(ir)*upf_pot%r(:)**2)
            DO ia = 1, na
               DO ib = ia, na
                  gmat(ia, ib, ir) = SUM(gbasis%bf(:, ia, la)*gbasis%bf(:, ib, la)*vgauss(:)*ww(:))
                  gmat(ib, ia, ir) = gmat(ia, ib, ir)
               END DO
            END DO
         END DO
         DO ir = 1, nl
            tv(ir) = tv(ir)+SUM(rmat(:, :)*gmat(:, :, ir))
            DO ix = ir, nl
               smat(ir, ix) = smat(ir, ix)+SUM(gmat(:, :, ix)*gmat(:, :, ir))
               smat(ix, ir) = smat(ir, ix)
            END DO
         END DO
         DEALLOCATE (rmat, gmat)
      END DO
      CALL get_pseudo_inverse_diag(smat(1:nl, 1:nl), sinv(1:nl, 1:nl), 1.e-10_dp)
      cl(1:nl) = MATMUL(sinv(1:nl, 1:nl), tv(1:nl))
      vgauss = 0.0_dp
      DO ir = 1, nl
         vgauss(:) = vgauss(:)+cl(ir)*EXP(-al(ir)*upf_pot%r(:)**2)
      END DO
      errloc = SUM((vgauss-vloc)**2*ww)
      ALLOCATE (sgp_pot%a_local(nl), sgp_pot%c_local(nl))
      sgp_pot%n_local = nl
      sgp_pot%a_local(1:nl) = al(1:nl)
      sgp_pot%c_local(1:nl) = cl(1:nl)
      DEALLOCATE (vloc, vgauss)
      DEALLOCATE (al, cl, tv, smat, sinv)
      CALL release_atom_basis(gbasis)
      !
      ptype = ADJUSTL(TRIM(upf_pot%pseudo_type))
      IF (ptype(1:2) == "NC" .OR. ptype(1:2) == "US") THEN
         nl_trans = .FALSE.
      ELSE IF (ptype(1:2) == "SL") THEN
         nl_trans = .TRUE.
      ELSE
         CPABORT("Pseudopotential type: ["//ADJUSTL(TRIM(ptype))//"] not known")
      END IF
      sgp_pot%has_local = .TRUE.
      sgp_pot%has_nonlocal = .TRUE.

      ! Non-local pseudopotential in Gaussian form
      IF (nl_trans) THEN
         ! semi local pseudopotential
         ! fit to nonlocal form
         ! get basis representation on UPF grid
         nl = 8
         !
         sgp_pot%n_nonlocal = nl
         sgp_pot%lmax = upf_pot%l_max
         ALLOCATE (sgp_pot%a_nonlocal(nl))
         ALLOCATE (sgp_pot%h_nonlocal(nl, 0:upf_pot%l_max))
         ALLOCATE (sgp_pot%c_nonlocal(nl, nl, 0:upf_pot%l_max))
         !
         ALLOCATE (al(nl), cl(nl))
         ALLOCATE (smat(nl, nl), sinv(nl, nl))
         ALLOCATE (tmat(nl, nl), cmat(nl, nl))
         al = 0.0_dp
         DO ir = 1, nl
            al(ir) = 10.0_dp*0.60_dp**(ir-1)
         END DO
         !
         sgp_pot%a_nonlocal(1:nl) = al(1:nl)
         !
         CALL atom_basis_gridrep(basis, gbasis, upf_pot%r, upf_pot%rab)
         ALLOCATE (pgauss(nr), vloc(nr))
         DO la = 0, upf_pot%l_max
            IF (la == upf_pot%l_local) CYCLE
            sgp_pot%is_nonlocal(la) = .TRUE.
            na = gbasis%nbas(la)
            ALLOCATE (score(na, na), qmat(na, nl))
            ! Reference matrix
            vloc(:) = upf_pot%vsemi(:, la+1)-upf_pot%vlocal(:)
            DO ia = 1, na
               DO ib = ia, na
                  score(ia, ib) = SUM(vloc(:)*gbasis%bf(:, ia, la)*gbasis%bf(:, ib, la)*ww(:))
                  score(ib, ia) = score(ia, ib)
               END DO
            END DO
            ! overlap basis with projectors
            DO ir = 1, nl
               pgauss(:) = EXP(-al(ir)*upf_pot%r(:)**2)*upf_pot%r(:)**la
               eee = rootpi/(2._dp**(la+2)*dfac(2*la+1))/(2._dp*al(ir))**(la+1.5_dp)
               pgauss(:) = pgauss(:)/SQRT(eee)
               DO ia = 1, na
                  qmat(ia, ir) = SUM(gbasis%bf(:, ia, la)*pgauss(:)*ww)
               END DO
            END DO
            ! tmat = qmat * score * qmat
            tmat(1:nl, 1:nl) = MATMUL(TRANSPOSE(qmat(1:na, 1:nl)), MATMUL(score(1:na, 1:na), qmat(1:na, 1:nl)))
            smat(1:nl, 1:nl) = MATMUL(TRANSPOSE(qmat(1:na, 1:nl)), qmat(1:na, 1:nl))
            CALL get_pseudo_inverse_diag(smat(1:nl, 1:nl), sinv(1:nl, 1:nl), 1.e-10_dp)
            cmat(1:nl, 1:nl) = MATMUL(sinv(1:nl, 1:nl), MATMUL(tmat(1:nl, 1:nl), sinv(1:nl, 1:nl)))
            cmat(1:nl, 1:nl) = (cmat(1:nl, 1:nl)+TRANSPOSE(cmat(1:nl, 1:nl)))*0.5_dp
            CALL diamat_all(cmat(1:nl, 1:nl), cl(1:nl), .TRUE.)
            !
            ! get back unnormalized Gaussians
            DO ir = 1, nl
               ei = rootpi/(2._dp**(la+2)*dfac(2*la+1))/(2._dp*al(ir))**(la+1.5_dp)
               cmat(ir, 1:nl) = cmat(ir, 1:nl)/SQRT(ei)
            END DO
            sgp_pot%h_nonlocal(1:nl, la) = cl(1:nl)
            sgp_pot%c_nonlocal(1:nl, 1:nl, la) = cmat(1:nl, 1:nl)
            sgp_pot%is_nonlocal(la) = .TRUE.
            DEALLOCATE (score, qmat)
         END DO
         CALL release_atom_basis(gbasis)
         DEALLOCATE (pgauss, vloc)
         DEALLOCATE (al, cl, smat, sinv, tmat, cmat)
      ELSE
         ! non local pseudopotential
         ALLOCATE (pproa(nr), pprob(nr), pgauss(nr))
         np = upf_pot%number_of_proj
         nl = 8
         ALLOCATE (al(nl), cl(nl))
         ALLOCATE (smat(nl, nl), sinv(nl, nl))
         ALLOCATE (tmat(nl, nl), cmat(nl, nl))
         al = 0.0_dp
         cl = 0.0_dp
         DO ir = 1, nl
            al(ir) = 10.0_dp*0.60_dp**(ir-1)
         END DO
         !
         sgp_pot%lmax = MAXVAL(upf_pot%lbeta(:))
         sgp_pot%n_nonlocal = nl
         ALLOCATE (sgp_pot%a_nonlocal(nl))
         ALLOCATE (sgp_pot%h_nonlocal(nl, 0:sgp_pot%lmax))
         ALLOCATE (sgp_pot%c_nonlocal(nl, nl, 0:sgp_pot%lmax))
         !
         sgp_pot%a_nonlocal(1:nl) = al(1:nl)
         !
         CALL atom_basis_gridrep(basis, gbasis, upf_pot%r, upf_pot%rab)
         DO la = 0, sgp_pot%lmax
            sgp_pot%is_nonlocal(la) = .TRUE.
            na = gbasis%nbas(la)
            ALLOCATE (score(na, na), qmat(na, nl))
            ! Reference matrix
            score = 0.0_dp
            DO ipa = 1, np
               lb = upf_pot%lbeta(ipa)
               IF (la /= lb) CYCLE
               pproa(:) = upf_pot%beta(:, ipa)
               DO ipb = 1, np
                  lb = upf_pot%lbeta(ipb)
                  IF (la /= lb) CYCLE
                  pprob(:) = upf_pot%beta(:, ipb)
                  eee = upf_pot%dion(ipa, ipb)
                  DO ia = 1, na
                     cpa = SUM(pproa(:)*gbasis%bf(:, ia, la)*ww(:))
                     DO ib = ia, na
                        cpb = SUM(pprob(:)*gbasis%bf(:, ib, la)*ww(:))
                        score(ia, ib) = score(ia, ib)+cpa*eee*cpb
                        score(ib, ia) = score(ia, ib)
                     END DO
                  END DO
               END DO
            END DO
            ! overlap basis with projectors
            DO ir = 1, nl
               pgauss(:) = EXP(-al(ir)*upf_pot%r(:)**2)*upf_pot%r(:)**la
               eee = rootpi/(2._dp**(la+2)*dfac(2*la+1))/(2._dp*al(ir))**(la+1.5_dp)
               pgauss(:) = pgauss(:)/SQRT(eee)
               DO ia = 1, na
                  qmat(ia, ir) = SUM(gbasis%bf(:, ia, la)*pgauss(:)*ww)
               END DO
            END DO
            ! tmat = qmat * score * qmat
            tmat(1:nl, 1:nl) = MATMUL(TRANSPOSE(qmat(1:na, 1:nl)), MATMUL(score(1:na, 1:na), qmat(1:na, 1:nl)))
            smat(1:nl, 1:nl) = MATMUL(TRANSPOSE(qmat(1:na, 1:nl)), qmat(1:na, 1:nl))
            CALL get_pseudo_inverse_diag(smat(1:nl, 1:nl), sinv(1:nl, 1:nl), 1.e-10_dp)
            cmat(1:nl, 1:nl) = MATMUL(sinv(1:nl, 1:nl), MATMUL(tmat(1:nl, 1:nl), sinv(1:nl, 1:nl)))
            cmat(1:nl, 1:nl) = (cmat(1:nl, 1:nl)+TRANSPOSE(cmat(1:nl, 1:nl)))*0.5_dp
            CALL diamat_all(cmat(1:nl, 1:nl), cl(1:nl), .TRUE.)
            !
            ! get back unnormalized Gaussians
            DO ir = 1, nl
               ei = rootpi/(2._dp**(la+2)*dfac(2*la+1))/(2._dp*al(ir))**(la+1.5_dp)
               cmat(ir, 1:nl) = cmat(ir, 1:nl)/SQRT(ei)
            END DO
            sgp_pot%h_nonlocal(1:nl, la) = cl(1:nl)
            sgp_pot%c_nonlocal(1:nl, 1:nl, la) = cmat(1:nl, 1:nl)
            sgp_pot%is_nonlocal(la) = .TRUE.
            DEALLOCATE (score, qmat)
         END DO
         CALL release_atom_basis(gbasis)
         DEALLOCATE (pgauss, pproa, pprob)
         DEALLOCATE (al, cl, smat, sinv, tmat, cmat)
      ENDIF

      IF (upf_pot%core_correction) THEN
         sgp_pot%has_nlcc = .TRUE.
      ELSE
         sgp_pot%has_nlcc = .FALSE.
         sgp_pot%n_nlcc = 0
      END IF

      ! fit core charge to Gaussian form
      IF (sgp_pot%has_nlcc) THEN
         ALLOCATE (ccharge(nr), cgauss(nr))
         ccharge(:) = upf_pot%rho_nlcc(:)
         nl = 8
         ALLOCATE (al(nl), cl(nl), tv(nl))
         ALLOCATE (smat(nl, nl), sinv(nl, nl))
         al = 0.0_dp
         cl = 0.0_dp
         DO ir = 1, nl
            al(ir) = 10.0_dp*0.6_dp**(ir-1)
         END DO
         ! calculate integrals
         smat = 0.0_dp
         sinv = 0.0_dp
         tv = 0.0_dp
         CALL sg_overlap(smat(1:nl, 1:nl), 0, al(1:nl), al(1:nl))
         DO ir = 1, nl
            cgauss(:) = EXP(-al(ir)*upf_pot%r(:)**2)
            tv(ir) = SUM(cgauss*ccharge*ww)
         END DO
         CALL get_pseudo_inverse_diag(smat(1:nl, 1:nl), sinv(1:nl, 1:nl), 1.e-10_dp)
         cl(1:nl) = MATMUL(sinv(1:nl, 1:nl), tv(1:nl))
         cgauss = 0.0_dp
         DO ir = 1, nl
            cgauss(:) = cgauss(:)+cl(ir)*EXP(-al(ir)*upf_pot%r(:)**2)
         END DO
         errcc = SUM((cgauss-ccharge)**2*ww)
         ALLOCATE (sgp_pot%a_local(nl), sgp_pot%c_local(nl))
         sgp_pot%n_nlcc = nl
         sgp_pot%a_nlcc(1:nl) = al(1:nl)
         sgp_pot%c_nlcc(1:nl) = cl(1:nl)
         DEALLOCATE (ccharge, cgauss)
         DEALLOCATE (al, cl, tv, smat, sinv)
      END IF

      DEALLOCATE (ww)

   END SUBROUTINE upf_sgp_constr

! **************************************************************************************************
!> \brief ...
!> \param sgp_pot ...
! **************************************************************************************************
   SUBROUTINE sgp_release(sgp_pot)

      TYPE(sgp_potential)                                :: sgp_pot

      CHARACTER(len=*), PARAMETER :: routineN = 'sgp_release', routineP = moduleN//':'//routineN

      IF (ASSOCIATED(sgp_pot%a_nonlocal)) DEALLOCATE (sgp_pot%a_nonlocal)
      IF (ASSOCIATED(sgp_pot%h_nonlocal)) DEALLOCATE (sgp_pot%h_nonlocal)
      IF (ASSOCIATED(sgp_pot%c_nonlocal)) DEALLOCATE (sgp_pot%c_nonlocal)

      IF (ASSOCIATED(sgp_pot%a_local)) DEALLOCATE (sgp_pot%a_local)
      IF (ASSOCIATED(sgp_pot%c_local)) DEALLOCATE (sgp_pot%c_local)

      IF (ASSOCIATED(sgp_pot%a_nlcc)) DEALLOCATE (sgp_pot%a_nlcc)
      IF (ASSOCIATED(sgp_pot%c_nlcc)) DEALLOCATE (sgp_pot%c_nlcc)

   END SUBROUTINE sgp_release

! **************************************************************************************************
!> \brief ...
!> \param core ...
!> \param hnl ...
!> \param basis ...
!> \param upf_pot ...
!> \param ac_local ...
! **************************************************************************************************
   SUBROUTINE upfints(core, hnl, basis, upf_pot, ac_local)
      REAL(KIND=dp), DIMENSION(:, :, 0:)                 :: core, hnl
      TYPE(atom_basis_type), INTENT(INOUT)               :: basis
      TYPE(atom_upfpot_type)                             :: upf_pot
      REAL(KIND=dp), INTENT(IN)                          :: ac_local

      CHARACTER(len=*), PARAMETER :: routineN = 'upfints', routineP = moduleN//':'//routineN

      CHARACTER(len=4)                                   :: ptype
      INTEGER                                            :: i, j, k1, k2, la, lb, m, n
      REAL(KIND=dp)                                      :: rc, zval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: spot
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: spmat
      TYPE(atom_basis_type)                              :: gbasis

      ! get basis representation on UPF grid
      CALL atom_basis_gridrep(basis, gbasis, upf_pot%r, upf_pot%rab)

      ! local pseudopotential
      core = 0._dp
      n = SIZE(upf_pot%r)
      ALLOCATE (spot(n))
      spot(:) = upf_pot%vlocal(:)
      zval = upf_pot%zion
      DO i = 1, n
         IF (upf_pot%r(i) < 1.e-12_dp) THEN
            spot(i) = spot(i)+2.0_dp*zval/rootpi/ac_local
         ELSE
            rc = upf_pot%r(i)/ac_local
            spot(i) = spot(i)+zval/upf_pot%r(i)*erf(rc)
         END IF
      END DO

      CALL numpot_matrix(core, spot, gbasis, 0)
      DEALLOCATE (spot)

      hnl = 0._dp
      ptype = ADJUSTL(TRIM(upf_pot%pseudo_type))
      IF (ptype(1:2) == "NC" .OR. ptype(1:2) == "US") THEN
         ! non local pseudopotential
         n = MAXVAL(gbasis%nbas(:))
         m = upf_pot%number_of_proj
         ALLOCATE (spmat(n, m))
         spmat = 0.0_dp
         DO i = 1, m
            la = upf_pot%lbeta(i)
            DO j = 1, gbasis%nbas(la)
               spmat(j, i) = integrate_grid(upf_pot%beta(:, i), gbasis%bf(:, j, la), gbasis%grid)
            END DO
         END DO
         DO i = 1, m
            la = upf_pot%lbeta(i)
            DO j = 1, m
               lb = upf_pot%lbeta(j)
               IF (la == lb) THEN
                  DO k1 = 1, gbasis%nbas(la)
                     DO k2 = 1, gbasis%nbas(la)
                        hnl(k1, k2, la) = hnl(k1, k2, la)+spmat(k1, i)*upf_pot%dion(i, j)*spmat(k2, j)
                     END DO
                  END DO
               END IF
            END DO
         END DO
         DEALLOCATE (spmat)
      ELSE IF (ptype(1:2) == "SL") THEN
         ! semi local pseudopotential
         DO la = 0, upf_pot%l_max
            IF (la == upf_pot%l_local) CYCLE
            m = SIZE(upf_pot%vsemi(:, la+1))
            ALLOCATE (spot(m))
            spot(:) = upf_pot%vsemi(:, la+1)-upf_pot%vlocal(:)
            n = basis%nbas(la)
            DO i = 1, n
               DO j = i, n
                  hnl(i, j, la) = hnl(i, j, la)+ &
                                  integrate_grid(spot(:), &
                                                 gbasis%bf(:, i, la), gbasis%bf(:, j, la), gbasis%grid)
                  hnl(j, i, la) = hnl(i, j, la)
               END DO
            END DO
            DEALLOCATE (spot)
         END DO
      ELSE
         CPABORT("Pseudopotential type: ["//ADJUSTL(TRIM(ptype))//"] not known")
      END IF

      ! release basis representation on UPF grid
      CALL release_atom_basis(gbasis)

   END SUBROUTINE upfints

! **************************************************************************************************
!> \brief ...
!> \param hnl ...
!> \param basis ...
!> \param ecp_pot ...
! **************************************************************************************************
   SUBROUTINE ecpints(hnl, basis, ecp_pot)
      REAL(KIND=dp), DIMENSION(:, :, 0:)                 :: hnl
      TYPE(atom_basis_type), INTENT(INOUT)               :: basis
      TYPE(atom_ecppot_type)                             :: ecp_pot

      CHARACTER(len=*), PARAMETER :: routineN = 'ecpints', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, k, l, m, n
      REAL(KIND=dp)                                      :: alpha
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cpot
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rad

      rad => basis%grid%rad
      m = basis%grid%nr
      ALLOCATE (cpot(1:m))

      ! non local pseudopotential
      hnl = 0.0_dp
      DO l = 0, ecp_pot%lmax
         cpot = 0._dp
         DO k = 1, ecp_pot%npot(l)
            n = ecp_pot%nrpot(k, l)
            alpha = ecp_pot%bpot(k, l)
            cpot(:) = cpot(:)+ecp_pot%apot(k, l)*rad(:)**(n-2)*EXP(-alpha*rad(:)**2)
         END DO
         DO i = 1, basis%nbas(l)
            DO j = i, basis%nbas(l)
               hnl(i, j, l) = integrate_grid(cpot, basis%bf(:, i, l), basis%bf(:, j, l), basis%grid)
               hnl(j, i, l) = hnl(i, j, l)
            END DO
         END DO
      END DO
      DEALLOCATE (cpot)

   END SUBROUTINE ecpints

! **************************************************************************************************
!> \brief ...
!> \param core ...
!> \param hnl ...
!> \param basis ...
!> \param sgp_pot ...
! **************************************************************************************************
   SUBROUTINE sgpints(core, hnl, basis, sgp_pot)
      REAL(KIND=dp), DIMENSION(:, :, 0:)                 :: core, hnl
      TYPE(atom_basis_type), INTENT(INOUT)               :: basis
      TYPE(sgp_potential)                                :: sgp_pot

      CHARACTER(len=*), PARAMETER :: routineN = 'sgpints', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, ia, j, l, m, n, na
      REAL(KIND=dp)                                      :: a, c, zval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cpot, pgauss
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: qmat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rad

      rad => basis%grid%rad
      m = basis%grid%nr

      ! local pseudopotential
      ALLOCATE (cpot(m))
      IF (sgp_pot%has_local) THEN
         zval = sgp_pot%zval
         core = 0._dp
         cpot = 0.0_dp
         DO i = 1, sgp_pot%n_local
            cpot(:) = cpot(:)+sgp_pot%c_local(i)*EXP(-sgp_pot%a_local(i)*rad(:)**2)
         END DO
         CALL numpot_matrix(core, cpot, basis, 0)
      END IF
      DEALLOCATE (cpot)

      ! nonlocal pseudopotential
      IF (sgp_pot%has_nonlocal) THEN
         hnl = 0.0_dp
         ALLOCATE (pgauss(1:m))
         n = sgp_pot%n_nonlocal
         !
         DO l = 0, sgp_pot%lmax
            CPASSERT(l <= UBOUND(basis%nbas, 1))
            IF (.NOT. sgp_pot%is_nonlocal(l)) CYCLE
            ! overlap (a|p)
            na = basis%nbas(l)
            ALLOCATE (qmat(na, n))
            DO i = 1, n
               pgauss(:) = 0.0_dp
               DO j = 1, n
                  a = sgp_pot%a_nonlocal(j)
                  c = sgp_pot%c_nonlocal(j, i, l)
                  pgauss(:) = pgauss(:)+c*EXP(-a*rad(:)**2)*rad(:)**l
               END DO
               DO ia = 1, na
                  qmat(ia, i) = SUM(basis%bf(:, ia, l)*pgauss(:)*basis%grid%wr(:))
               END DO
            END DO
            DO i = 1, na
               DO j = i, na
                  DO ia = 1, n
                     hnl(i, j, l) = hnl(i, j, l)+qmat(i, ia)*qmat(j, ia)*sgp_pot%h_nonlocal(ia, l)
                  END DO
                  hnl(j, i, l) = hnl(i, j, l)
               END DO
            END DO
            DEALLOCATE (qmat)
         END DO
         DEALLOCATE (pgauss)
      END IF

   END SUBROUTINE sgpints

END MODULE atom_sgp
