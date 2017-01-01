!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE atom_admm_methods
   USE atom_operators,                  ONLY: atom_basis_projection_overlap,&
                                              atom_int_release,&
                                              atom_int_setup
   USE atom_output,                     ONLY: atom_print_basis
   USE atom_types,                      ONLY: &
        atom_basis_type, atom_integrals, atom_orbitals, atom_p_type, atom_type, create_atom_orbs, &
        init_atom_basis, release_atom_basis, release_atom_orbs
   USE atom_utils,                      ONLY: atom_consistent_method,&
                                              atom_denmat,&
                                              atom_trace,&
                                              eeri_contract
   USE atom_xc,                         ONLY: atom_dpot_lda,&
                                              atom_vxc_lda,&
                                              atom_vxc_lsd
   USE input_constants,                 ONLY: do_rhf_atom,&
                                              do_rks_atom,&
                                              do_rohf_atom,&
                                              do_uhf_atom,&
                                              do_uks_atom,&
                                              xc_funct_no_shortcut
   USE input_section_types,             ONLY: section_vals_duplicate,&
                                              section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals2,&
                                              section_vals_release,&
                                              section_vals_remove_values,&
                                              section_vals_type,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invmat_symm,&
                                              jacobi
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC  :: atom_admm

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_admm_methods'

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param atom_info ...
!> \param admm_section ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE atom_admm(atom_info, admm_section, iw)

      TYPE(atom_p_type), DIMENSION(:, :), POINTER        :: atom_info
      TYPE(section_vals_type), POINTER                   :: admm_section
      INTEGER, INTENT(IN)                                :: iw

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_admm', routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: btyp
      INTEGER                                            :: i, ifun, j, l, m, maxl, mo, n, na, nb, &
                                                            zval
      LOGICAL                                            :: pp_calc, rhf
      REAL(KIND=dp) :: admm1_k_energy, admm2_k_energy, admmq_k_energy, dfexc_admm1, dfexc_admm2, &
         dfexc_admmq, dxc, dxk, el1, el2, elq, elref, fexc_optx_admm1, fexc_optx_admm2, &
         fexc_optx_admmq, fexc_optx_ref, fexc_pbex_admm1, fexc_pbex_admm2, fexc_pbex_admmq, &
         fexc_pbex_ref, ref_energy, xsi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: lmat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :) :: admm1_k, admm2_k, admm_xcmat, admm_xcmata, &
         admm_xcmatb, admmq_k, ovlap, ref_k, ref_xcmat, ref_xcmata, ref_xcmatb, sinv, siref, tmat
      TYPE(atom_basis_type), POINTER                     :: admm_basis, ref_basis
      TYPE(atom_integrals), POINTER                      :: admm_int, ref_int
      TYPE(atom_orbitals), POINTER                       :: admm1_orbs, admm2_orbs, admmq_orbs, &
                                                            ref_orbs
      TYPE(atom_type), POINTER                           :: atom
      TYPE(section_vals_type), POINTER                   :: basis_section, xc_fun, xc_fun_section, &
                                                            xc_optx, xc_pbex, xc_section

      IF (iw > 0) THEN
         WRITE (iw, '(/,T2,A)') &
            '!-----------------------------------------------------------------------------!'
         WRITE (iw, '(T30,A)') "Analysis of ADMM approximations"
         WRITE (iw, '(T2,A)') &
            '!-----------------------------------------------------------------------------!'
      END IF

      ! setup an xc section
      NULLIFY (xc_section, xc_pbex, xc_optx)
      CALL section_vals_duplicate(atom_info(1, 1)%atom%xc_section, xc_section)
      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      ! Overwrite possible shortcut
      CALL section_vals_val_set(xc_fun_section, "_SECTION_PARAMETERS_", &
                                i_val=xc_funct_no_shortcut)
      ifun = 0
      DO
         ifun = ifun+1
         xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
         IF (.NOT. ASSOCIATED(xc_fun)) EXIT
         CALL section_vals_remove_values(xc_fun)
      END DO
      ! PBEX
      CALL section_vals_duplicate(xc_section, xc_pbex)
      xc_fun_section => section_vals_get_subs_vals(xc_pbex, "XC_FUNCTIONAL")
      CALL section_vals_val_set(xc_fun_section, "PBE%_SECTION_PARAMETERS_", l_val=.TRUE.)
      CALL section_vals_val_set(xc_fun_section, "PBE%SCALE_X", r_val=1.0_dp)
      CALL section_vals_val_set(xc_fun_section, "PBE%SCALE_C", r_val=0.0_dp)
      ! OPTX
      CALL section_vals_duplicate(xc_section, xc_optx)
      xc_fun_section => section_vals_get_subs_vals(xc_optx, "XC_FUNCTIONAL")
      CALL section_vals_val_set(xc_fun_section, "OPTX%_SECTION_PARAMETERS_", l_val=.TRUE.)
      CALL section_vals_val_set(xc_fun_section, "OPTX%SCALE_X", r_val=1.0_dp)

      ! ADMM basis set
      zval = atom_info(1, 1)%atom%z
      pp_calc = atom_info(1, 1)%atom%pp_calc
      btyp = "AE"
      IF (pp_calc) btyp = "PP"
      ALLOCATE (admm_basis)
      basis_section => section_vals_get_subs_vals(admm_section, "ADMM_BASIS")
      NULLIFY (admm_basis%grid)
      CALL init_atom_basis(admm_basis, basis_section, zval, btyp)
      IF (iw > 0) THEN
         CALL atom_print_basis(admm_basis, iw, " ADMM Basis")
      END IF
      ! reference basis set
      ref_basis => atom_info(1, 1)%atom%basis
      ! integrals
      ALLOCATE (ref_int, admm_int)
      CALL atom_int_setup(ref_int, ref_basis, eri_exchange=.TRUE.)
      CALL atom_int_setup(admm_int, admm_basis, eri_exchange=.TRUE.)
      DO l = 0, 3
         IF (admm_int%n(l) /= admm_int%nne(l)) THEN
            IF (iw > 0) WRITE (iw, *) "ADMM Basis is linear dependent ", l, admm_int%n(l), admm_int%nne(l)
            CPABORT("ADMM basis is linear dependent")
         END IF
      END DO
      ! mixed overlap
      na = MAXVAL(admm_basis%nbas(:))
      nb = MAXVAL(ref_basis%nbas(:))
      ALLOCATE (ovlap(1:na, 1:nb, 0:3))
      CALL atom_basis_projection_overlap(ovlap, admm_basis, ref_basis)
      ! Inverse of ADMM overlap matrix
      ALLOCATE (sinv(1:na, 1:na, 0:3))
      DO l = 0, 3
         n = admm_basis%nbas(l)
         IF (n < 1) CYCLE
         sinv(1:n, 1:n, l) = admm_int%ovlp(1:n, 1:n, l)
         CALL invmat_symm(sinv(1:n, 1:n, l))
      END DO
      ! ADMM transformation matrix
      ALLOCATE (tmat(1:na, 1:nb, 0:3))
      DO l = 0, 3
         n = admm_basis%nbas(l)
         m = ref_basis%nbas(l)
         IF (n < 1 .OR. m < 1) CYCLE
         tmat(1:n, 1:m, l) = MATMUL(sinv(1:n, 1:n, l), ovlap(1:n, 1:m, l))
      END DO
      ! Inverse of REF overlap matrix
      ALLOCATE (siref(1:nb, 1:nb, 0:3))
      DO l = 0, 3
         n = ref_basis%nbas(l)
         IF (n < 1) CYCLE
         siref(1:n, 1:n, l) = ref_int%ovlp(1:n, 1:n, l)
         CALL invmat_symm(siref(1:n, 1:n, l))
      END DO

      DO i = 1, SIZE(atom_info, 1)
         DO j = 1, SIZE(atom_info, 2)
            atom => atom_info(i, j)%atom
            IF (atom_consistent_method(atom%method_type, atom%state%multiplicity)) THEN
               ref_orbs => atom%orbitals
               ALLOCATE (ref_k(1:nb, 1:nb, 0:3))
               SELECT CASE (atom%method_type)
               CASE (do_rks_atom, do_rhf_atom)
                  ! restricted
                  rhf = .TRUE.
                  ref_k = 0.0_dp
                  CALL eeri_contract(ref_k, ref_int%eeri, ref_orbs%pmat, ref_basis%nbas)
                  ref_energy = 0.5_dp*atom_trace(ref_k, ref_orbs%pmat)
               CASE (do_uks_atom, do_uhf_atom)
                  ! unrestricted
                  rhf = .FALSE.
                  ref_k = 0.0_dp
                  CALL eeri_contract(ref_k, ref_int%eeri, ref_orbs%pmata, ref_basis%nbas)
                  ref_energy = atom_trace(ref_k, ref_orbs%pmata)
                  ref_k = 0.0_dp
                  CALL eeri_contract(ref_k, ref_int%eeri, ref_orbs%pmatb, ref_basis%nbas)
                  ref_energy = ref_energy+atom_trace(ref_k, ref_orbs%pmatb)
               CASE (do_rohf_atom)
                  CPABORT("ADMM not available")
               CASE DEFAULT
                  CPABORT("ADMM not available")
               END SELECT
               DEALLOCATE (ref_k)
               ! reference number of electrons
               elref = atom_trace(ref_int%ovlp, ref_orbs%pmat)
               ! admm orbitals and density matrices
               mo = MAXVAL(atom%state%maxn_calc)
               NULLIFY (admm1_orbs, admm2_orbs, admmq_orbs)
               CALL create_atom_orbs(admm1_orbs, na, mo)
               CALL create_atom_orbs(admm2_orbs, na, mo)
               CALL create_atom_orbs(admmq_orbs, na, mo)
               ALLOCATE (lmat(1:mo, 1:mo))
               ALLOCATE (admm1_k(1:na, 1:na, 0:3))
               ALLOCATE (admm2_k(1:na, 1:na, 0:3))
               ALLOCATE (admmq_k(1:na, 1:na, 0:3))
               IF (rhf) THEN
                  DO l = 0, 3
                     n = admm_basis%nbas(l)
                     m = ref_basis%nbas(l)
                     mo = atom%state%maxn_calc(l)
                     IF (n < 1 .OR. m < 1 .OR. mo < 1) CYCLE
                     admm2_orbs%wfn(1:n, 1:mo, l) = MATMUL(tmat(1:n, 1:m, l), ref_orbs%wfn(1:m, 1:mo, l))
                     CALL lowdin_matrix(admm2_orbs%wfn(1:n, 1:mo, l), lmat(1:mo, 1:mo), admm_int%ovlp(1:n, 1:n, l))
                     admm1_orbs%wfn(1:n, 1:mo, l) = MATMUL(admm2_orbs%wfn(1:n, 1:mo, l), lmat(1:mo, 1:mo))
                  END DO
                  CALL atom_denmat(admm1_orbs%pmat, admm1_orbs%wfn, admm_basis%nbas, atom%state%occupation, &
                                   atom%state%maxl_occ, atom%state%maxn_occ)
                  CALL atom_denmat(admm2_orbs%pmat, admm2_orbs%wfn, admm_basis%nbas, atom%state%occupation, &
                                   atom%state%maxl_occ, atom%state%maxn_occ)
                  el1 = atom_trace(admm_int%ovlp, admm1_orbs%pmat)
                  el2 = atom_trace(admm_int%ovlp, admm2_orbs%pmat)
                  xsi = elref/el2
                  admmq_orbs%pmat = xsi*admm2_orbs%pmat
                  elq = atom_trace(admm_int%ovlp, admmq_orbs%pmat)
                  admmq_orbs%wfn = SQRT(xsi)*admm2_orbs%wfn
                  !
                  admm1_k = 0.0_dp
                  CALL eeri_contract(admm1_k, admm_int%eeri, admm1_orbs%pmat, admm_basis%nbas)
                  admm1_k_energy = 0.5_dp*atom_trace(admm1_k, admm1_orbs%pmat)
                  admm2_k = 0.0_dp
                  CALL eeri_contract(admm2_k, admm_int%eeri, admm2_orbs%pmat, admm_basis%nbas)
                  admm2_k_energy = 0.5_dp*atom_trace(admm2_k, admm2_orbs%pmat)
                  admmq_k = 0.0_dp
                  CALL eeri_contract(admmq_k, admm_int%eeri, admmq_orbs%pmat, admm_basis%nbas)
                  admmq_k_energy = 0.5_dp*atom_trace(admmq_k, admmq_orbs%pmat)
               ELSE
                  DO l = 0, 3
                     n = admm_basis%nbas(l)
                     m = ref_basis%nbas(l)
                     mo = atom%state%maxn_calc(l)
                     IF (n < 1 .OR. m < 1 .OR. mo < 1) CYCLE
                     admm2_orbs%wfna(1:n, 1:mo, l) = MATMUL(tmat(1:n, 1:m, l), ref_orbs%wfna(1:m, 1:mo, l))
                     CALL lowdin_matrix(admm2_orbs%wfna(1:n, 1:mo, l), lmat, admm_int%ovlp(1:n, 1:n, l))
                     admm1_orbs%wfna(1:n, 1:mo, l) = MATMUL(admm2_orbs%wfna(1:n, 1:mo, l), lmat(1:mo, 1:mo))
                     admm2_orbs%wfnb(1:n, 1:mo, l) = MATMUL(tmat(1:n, 1:m, l), ref_orbs%wfnb(1:m, 1:mo, l))
                     CALL lowdin_matrix(admm2_orbs%wfnb(1:n, 1:mo, l), lmat, admm_int%ovlp(1:n, 1:n, l))
                     admm1_orbs%wfnb(1:n, 1:mo, l) = MATMUL(admm2_orbs%wfnb(1:n, 1:mo, l), lmat(1:mo, 1:mo))
                  END DO
                  CALL atom_denmat(admm1_orbs%pmata, admm1_orbs%wfna, admm_basis%nbas, atom%state%occa, &
                                   atom%state%maxl_occ, atom%state%maxn_occ)
                  CALL atom_denmat(admm1_orbs%pmatb, admm1_orbs%wfnb, admm_basis%nbas, atom%state%occb, &
                                   atom%state%maxl_occ, atom%state%maxn_occ)
                  admm1_orbs%pmat = admm1_orbs%pmata+admm1_orbs%pmatb
                  CALL atom_denmat(admm2_orbs%pmata, admm2_orbs%wfna, admm_basis%nbas, atom%state%occa, &
                                   atom%state%maxl_occ, atom%state%maxn_occ)
                  CALL atom_denmat(admm2_orbs%pmatb, admm2_orbs%wfnb, admm_basis%nbas, atom%state%occb, &
                                   atom%state%maxl_occ, atom%state%maxn_occ)
                  admm2_orbs%pmat = admm2_orbs%pmata+admm2_orbs%pmatb
                  elref = atom_trace(ref_int%ovlp, ref_orbs%pmata)
                  el2 = atom_trace(admm_int%ovlp, admm2_orbs%pmata)
                  xsi = elref/el2
                  admmq_orbs%pmata = xsi*admm2_orbs%pmata
                  admmq_orbs%wfna = SQRT(xsi)*admm2_orbs%wfna
                  elref = atom_trace(ref_int%ovlp, ref_orbs%pmatb)
                  el2 = atom_trace(admm_int%ovlp, admm2_orbs%pmatb)
                  xsi = elref/el2
                  admmq_orbs%pmatb = xsi*admm2_orbs%pmatb
                  admmq_orbs%wfnb = SQRT(xsi)*admm2_orbs%wfnb
                  admmq_orbs%pmat = admmq_orbs%pmata+admmq_orbs%pmatb
                  el1 = atom_trace(admm_int%ovlp, admm1_orbs%pmat)
                  el2 = atom_trace(admm_int%ovlp, admm2_orbs%pmat)
                  elq = atom_trace(admm_int%ovlp, admmq_orbs%pmat)
                  elref = atom_trace(ref_int%ovlp, ref_orbs%pmat)
                  !
                  admm1_k = 0.0_dp
                  CALL eeri_contract(admm1_k, admm_int%eeri, admm1_orbs%pmata, admm_basis%nbas)
                  admm1_k_energy = atom_trace(admm1_k, admm1_orbs%pmata)
                  admm1_k = 0.0_dp
                  CALL eeri_contract(admm1_k, admm_int%eeri, admm1_orbs%pmatb, admm_basis%nbas)
                  admm1_k_energy = admm1_k_energy+atom_trace(admm1_k, admm1_orbs%pmatb)
                  admm2_k = 0.0_dp
                  CALL eeri_contract(admm2_k, admm_int%eeri, admm2_orbs%pmata, admm_basis%nbas)
                  admm2_k_energy = atom_trace(admm2_k, admm2_orbs%pmata)
                  admm2_k = 0.0_dp
                  CALL eeri_contract(admm2_k, admm_int%eeri, admm2_orbs%pmatb, admm_basis%nbas)
                  admm2_k_energy = admm2_k_energy+atom_trace(admm2_k, admm2_orbs%pmatb)
                  admmq_k = 0.0_dp
                  CALL eeri_contract(admmq_k, admm_int%eeri, admmq_orbs%pmata, admm_basis%nbas)
                  admmq_k_energy = atom_trace(admmq_k, admmq_orbs%pmata)
                  admmq_k = 0.0_dp
                  CALL eeri_contract(admmq_k, admm_int%eeri, admmq_orbs%pmatb, admm_basis%nbas)
                  admmq_k_energy = admmq_k_energy+atom_trace(admmq_k, admmq_orbs%pmatb)
               END IF
               DEALLOCATE (lmat)
               !
               ! ADMM correction terms
               maxl = atom%state%maxl_occ
               IF (rhf) THEN
                  ALLOCATE (ref_xcmat(1:nb, 1:nb, 0:3))
                  ALLOCATE (admm_xcmat(1:na, 1:na, 0:3))
                  ! PBEX
                  CALL atom_vxc_lda(ref_basis, ref_orbs%pmat, maxl, xc_pbex, fexc_pbex_ref, ref_xcmat)
                  CALL atom_vxc_lda(admm_basis, admm1_orbs%pmat, maxl, xc_pbex, fexc_pbex_admm1, admm_xcmat)
                  CALL atom_vxc_lda(admm_basis, admm2_orbs%pmat, maxl, xc_pbex, fexc_pbex_admm2, admm_xcmat)
                  CALL atom_vxc_lda(admm_basis, admmq_orbs%pmat, maxl, xc_pbex, fexc_pbex_admmq, admm_xcmat)
                  ! OPTX
                  CALL atom_vxc_lda(ref_basis, ref_orbs%pmat, maxl, xc_optx, fexc_optx_ref, ref_xcmat)
                  CALL atom_vxc_lda(admm_basis, admm1_orbs%pmat, maxl, xc_optx, fexc_optx_admm1, admm_xcmat)
                  CALL atom_vxc_lda(admm_basis, admm2_orbs%pmat, maxl, xc_optx, fexc_optx_admm2, admm_xcmat)
                  CALL atom_vxc_lda(admm_basis, admmq_orbs%pmat, maxl, xc_optx, fexc_optx_admmq, admm_xcmat)
                  ! LINX
                  CALL atom_dpot_lda(ref_basis, ref_orbs%pmat, admm_basis, admm1_orbs%pmat, &
                                     maxl, "LINX", dfexc_admm1)
                  CALL atom_dpot_lda(ref_basis, ref_orbs%pmat, admm_basis, admm2_orbs%pmat, &
                                     maxl, "LINX", dfexc_admm2)
                  CALL atom_dpot_lda(ref_basis, ref_orbs%pmat, admm_basis, admmq_orbs%pmat, &
                                     maxl, "LINX", dfexc_admmq)
                  DEALLOCATE (ref_xcmat, admm_xcmat)
               ELSE
                  ALLOCATE (ref_xcmata(1:nb, 1:nb, 0:3))
                  ALLOCATE (ref_xcmatb(1:nb, 1:nb, 0:3))
                  ALLOCATE (admm_xcmata(1:na, 1:na, 0:3))
                  ALLOCATE (admm_xcmatb(1:na, 1:na, 0:3))
                  ! PBEX
                  CALL atom_vxc_lsd(ref_basis, ref_orbs%pmata, ref_orbs%pmatb, maxl, xc_pbex, fexc_pbex_ref, &
                                    ref_xcmata, ref_xcmatb)
                  CALL atom_vxc_lsd(admm_basis, admm1_orbs%pmata, admm1_orbs%pmatb, maxl, xc_pbex, &
                                    fexc_pbex_admm1, admm_xcmata, admm_xcmatb)
                  CALL atom_vxc_lsd(admm_basis, admm2_orbs%pmata, admm2_orbs%pmatb, maxl, xc_pbex, &
                                    fexc_pbex_admm2, admm_xcmata, admm_xcmatb)
                  CALL atom_vxc_lsd(admm_basis, admmq_orbs%pmata, admmq_orbs%pmatb, maxl, xc_pbex, &
                                    fexc_pbex_admmq, admm_xcmata, admm_xcmatb)
                  CALL atom_vxc_lsd(ref_basis, ref_orbs%pmata, ref_orbs%pmatb, maxl, xc_optx, fexc_optx_ref, &
                                    ref_xcmata, ref_xcmatb)
                  CALL atom_vxc_lsd(admm_basis, admm1_orbs%pmata, admm1_orbs%pmatb, maxl, xc_optx, &
                                    fexc_optx_admm1, admm_xcmata, admm_xcmatb)
                  CALL atom_vxc_lsd(admm_basis, admm2_orbs%pmata, admm2_orbs%pmatb, maxl, xc_optx, &
                                    fexc_optx_admm2, admm_xcmata, admm_xcmatb)
                  CALL atom_vxc_lsd(admm_basis, admmq_orbs%pmata, admmq_orbs%pmatb, maxl, xc_optx, &
                                    fexc_optx_admmq, admm_xcmata, admm_xcmatb)
                  DEALLOCATE (ref_xcmata, ref_xcmatb, admm_xcmata, admm_xcmatb)
               END IF
               !

               IF (iw > 0) THEN
                  WRITE (iw, "(/,A,I3,T48,A,I3)") " Electronic Structure Setting: ", i, &
                     " Electronic Structure Method: ", j
                  WRITE (iw, "(' Norm of ADMM Basis projection ',T61,F20.10)") el2/elref
                  WRITE (iw, "(' Reference Exchange Energy [Hartree]',T61,F20.10)") ref_energy
                  ! ADMM1
                  dxk = ref_energy-admm1_k_energy
                  WRITE (iw, "(A,F20.10,T60,A,F13.10)") " ADMM1 METHOD: Energy ", admm1_k_energy, &
                     " Error: ", dxk
                  dxc = fexc_pbex_ref-fexc_pbex_admm1
                  WRITE (iw, "(T10,A,F12.6,F12.3,'%',T60,A,F13.10)") "PBEX Correction  ", dxc, dxc/dxk*100._dp, &
                     " Error: ", dxk-dxc
                  dxc = fexc_optx_ref-fexc_optx_admm1
                  WRITE (iw, "(T10,A,F12.6,F12.3,'%',T60,A,F13.10)") "OPTX Correction  ", dxc, dxc/dxk*100._dp, &
                     " Error: ", dxk-dxc
                  dxc = dfexc_admm1
                  WRITE (iw, "(T10,A,F12.6,F12.3,'%',T60,A,F13.10)") "LINX Correction  ", dxc, dxc/dxk*100._dp, &
                     " Error: ", dxk-dxc
                  ! ADMM2
                  dxk = ref_energy-admm2_k_energy
                  WRITE (iw, "(A,F20.10,T60,A,F13.10)") " ADMM2 METHOD: Energy ", admm2_k_energy, &
                     " Error: ", dxk
                  dxc = fexc_pbex_ref-fexc_pbex_admm2
                  WRITE (iw, "(T10,A,F12.6,F12.3,'%',T60,A,F13.10)") "PBEX Correction  ", dxc, dxc/dxk*100._dp, &
                     " Error: ", dxk-dxc
                  dxc = fexc_optx_ref-fexc_optx_admm2
                  WRITE (iw, "(T10,A,F12.6,F12.3,'%',T60,A,F13.10)") "OPTX Correction  ", dxc, dxc/dxk*100._dp, &
                     " Error: ", dxk-dxc
                  dxc = dfexc_admm2
                  WRITE (iw, "(T10,A,F12.6,F12.3,'%',T60,A,F13.10)") "LINX Correction  ", dxc, dxc/dxk*100._dp, &
                     " Error: ", dxk-dxc
                  ! ADMMQ
                  dxk = ref_energy-admmq_k_energy
                  WRITE (iw, "(A,F20.10,T60,A,F13.10)") " ADMMQ METHOD: Energy ", admmq_k_energy, &
                     " Error: ", dxk
                  dxc = fexc_pbex_ref-fexc_pbex_admmq
                  WRITE (iw, "(T10,A,F12.6,F12.3,'%',T60,A,F13.10)") "PBEX Correction  ", dxc, dxc/dxk*100._dp, &
                     " Error: ", dxk-dxc
                  dxc = fexc_optx_ref-fexc_optx_admmq
                  WRITE (iw, "(T10,A,F12.6,F12.3,'%',T60,A,F13.10)") "OPTX Correction  ", dxc, dxc/dxk*100._dp, &
                     " Error: ", dxk-dxc
                  dxc = dfexc_admmq
                  WRITE (iw, "(T10,A,F12.6,F12.3,'%',T60,A,F13.10)") "LINX Correction  ", dxc, dxc/dxk*100._dp, &
                     " Error: ", dxk-dxc
                  ! ADMMS
                  dxk = ref_energy-admmq_k_energy
                  WRITE (iw, "(A,F20.10,T60,A,F13.10)") " ADMMS METHOD: Energy ", admmq_k_energy, &
                     " Error: ", dxk
                  dxc = fexc_pbex_ref-fexc_pbex_admmq*xsi**(2._dp/3._dp)
                  WRITE (iw, "(T10,A,F12.6,F12.3,'%',T60,A,F13.10)") "PBEX Correction  ", dxc, dxc/dxk*100._dp, &
                     " Error: ", dxk-dxc
                  dxc = fexc_optx_ref-fexc_optx_admmq*xsi**(2._dp/3._dp)
                  WRITE (iw, "(T10,A,F12.6,F12.3,'%',T60,A,F13.10)") "OPTX Correction  ", dxc, dxc/dxk*100._dp, &
                     " Error: ", dxk-dxc
               END IF
               !
               DEALLOCATE (admm1_k, admm2_k, admmq_k)
               !
               CALL release_atom_orbs(admm1_orbs)
               CALL release_atom_orbs(admm2_orbs)
               CALL release_atom_orbs(admmq_orbs)
            END IF
         END DO
      END DO

      ! clean up
      CALL atom_int_release(ref_int)
      CALL atom_int_release(admm_int)
      CALL release_atom_basis(admm_basis)
      DEALLOCATE (ref_int, admm_int, admm_basis)
      DEALLOCATE (ovlap, sinv, tmat, siref)

      CALL section_vals_release(xc_pbex)
      CALL section_vals_release(xc_optx)
      CALL section_vals_release(xc_section)

      IF (iw > 0) THEN
         WRITE (iw, '(/,T2,A)') &
            '!------------------------------End of ADMM analysis---------------------------!'
      END IF

   END SUBROUTINE atom_admm

! **************************************************************************************************
!> \brief ...
!> \param wfn ...
!> \param lmat ...
!> \param ovlp ...
! **************************************************************************************************
   SUBROUTINE lowdin_matrix(wfn, lmat, ovlp)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: wfn
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: lmat
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ovlp

      INTEGER                                            :: i, j, k, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: w
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: vmat

      n = SIZE(wfn, 2)
      IF (n > 0) THEN
         lmat = MATMUL(TRANSPOSE(wfn), MATMUL(ovlp, wfn))
         ALLOCATE (w(n), vmat(n, n))
         CALL jacobi(lmat(1:n, 1:n), w(1:n), vmat(1:n, 1:n))
         w(1:n) = 1.0_dp/SQRT(w(1:n))
         DO i = 1, n
            DO j = 1, n
               lmat(i, j) = 0.0_dp
               DO k = 1, n
                  lmat(i, j) = lmat(i, j)+vmat(i, k)*w(k)*vmat(j, k)
               END DO
            END DO
         END DO
         DEALLOCATE (vmat, w)
      END IF

   END SUBROUTINE lowdin_matrix
! **************************************************************************************************

END MODULE atom_admm_methods
