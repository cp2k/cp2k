!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief calculate the orbitals for a given atomic kind type
! **************************************************************************************************
MODULE atom_kind_orbitals
   USE ai_onecenter,                    ONLY: sg_erfc
   USE atom_electronic_structure,       ONLY: calculate_atom
   USE atom_fit,                        ONLY: atom_fit_density
   USE atom_operators,                  ONLY: atom_int_release,&
                                              atom_int_setup,&
                                              atom_ppint_release,&
                                              atom_ppint_setup,&
                                              atom_relint_release,&
                                              atom_relint_setup
   USE atom_types,                      ONLY: &
        CGTO_BASIS, Clementi_geobas, GTO_BASIS, atom_basis_type, atom_gthpot_type, atom_integrals, &
        atom_orbitals, atom_potential_type, atom_sgppot_type, atom_type, create_atom_orbs, &
        create_atom_type, lmat, release_atom_basis, release_atom_potential, release_atom_type, &
        set_atom
   USE atom_utils,                      ONLY: atom_density,&
                                              get_maxl_occ,&
                                              get_maxn_occ
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE external_potential_types,        ONLY: all_potential_type,&
                                              get_potential,&
                                              gth_potential_type,&
                                              sgp_potential_type
   USE input_constants,                 ONLY: &
        barrier_conf, do_analytic, do_dkh0_atom, do_dkh1_atom, do_dkh2_atom, do_dkh3_atom, &
        do_gapw_log, do_nonrel_atom, do_numeric, do_rks_atom, do_sczoramp_atom, do_uks_atom, &
        do_zoramp_atom, gth_pseudo, no_pseudo, poly_conf, rel_dkh, rel_none, rel_sczora_mp, &
        rel_zora, rel_zora_full, rel_zora_mp, sgp_pseudo
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: dfac,&
                                              pi
   USE periodic_table,                  ONLY: ptable
   USE physcon,                         ONLY: bohr
   USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                              create_grid_atom,&
                                              grid_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              init_atom_electronic_state,&
                                              qs_kind_type,&
                                              set_pseudo_state
   USE rel_control_types,               ONLY: rel_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_kind_orbitals'

   INTEGER, PARAMETER                                 :: nua = 40, nup = 20
   REAL(KIND=dp), DIMENSION(nua), PARAMETER :: ugbs = (/0.007299_dp, 0.013705_dp, 0.025733_dp, &
                                        0.048316_dp, 0.090718_dp, 0.170333_dp, 0.319819_dp, 0.600496_dp, 1.127497_dp, 2.117000_dp, &
                                                 3.974902_dp, 7.463317_dp, 14.013204_dp, 26.311339_dp, 49.402449_dp, 92.758561_dp, &
                                                      174.164456_dp, 327.013024_dp, 614.003114_dp, 1152.858743_dp, 2164.619772_dp, &
                                                4064.312984_dp, 7631.197056_dp, 14328.416324_dp, 26903.186074_dp, 50513.706789_dp, &
                                         94845.070265_dp, 178082.107320_dp, 334368.848683_dp, 627814.487663_dp, 1178791.123851_dp, &
                                                      2213310.684886_dp, 4155735.557141_dp, 7802853.046713_dp, 14650719.428954_dp, &
                                                  27508345.793637_dp, 51649961.080194_dp, 96978513.342764_dp, 182087882.613702_dp, &
                                                        341890134.751331_dp/)

   PUBLIC :: calculate_atomic_orbitals, calculate_atomic_density, &
             calculate_atomic_relkin, gth_potential_conversion

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param atomic_kind ...
!> \param qs_kind ...
!> \param agrid ...
!> \param iunit ...
!> \param pmat ...
!> \param fmat ...
!> \param density ...
!> \param wavefunction ...
!> \param wfninfo ...
!> \param confine ...
!> \param xc_section ...
!> \param nocc ...
! **************************************************************************************************
   SUBROUTINE calculate_atomic_orbitals(atomic_kind, qs_kind, agrid, iunit, pmat, fmat, &
                                        density, wavefunction, wfninfo, confine, xc_section, nocc)
      TYPE(atomic_kind_type), INTENT(IN)                 :: atomic_kind
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind
      TYPE(grid_atom_type), OPTIONAL                     :: agrid
      INTEGER, INTENT(IN), OPTIONAL                      :: iunit
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: pmat, fmat
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: density
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: wavefunction, wfninfo
      LOGICAL, INTENT(IN), OPTIONAL                      :: confine
      TYPE(section_vals_type), OPTIONAL, POINTER         :: xc_section
      INTEGER, DIMENSION(:), OPTIONAL                    :: nocc

      INTEGER                                            :: i, ii, ipgf, j, k, k1, k2, l, ll, m, mb, &
                                                            mo, ngp, nj, nr, ns, nset, nsgf, nu, &
                                                            quadtype, z
      INTEGER, DIMENSION(0:lmat, 10)                     :: ncalc, ncore, nelem
      INTEGER, DIMENSION(0:lmat, 100)                    :: set_index, shell_index
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf, last_sgf, ls
      LOGICAL                                            :: has_basis, uks
      REAL(KIND=dp)                                      :: al, ear, ok, rk, scal, zeff
      REAL(KIND=dp), DIMENSION(0:lmat, 10, 2)            :: edelta
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: gcc
      TYPE(all_potential_type), POINTER                  :: all_potential
      TYPE(atom_basis_type), POINTER                     :: basis
      TYPE(atom_integrals), POINTER                      :: integrals
      TYPE(atom_orbitals), POINTER                       :: orbitals
      TYPE(atom_potential_type), POINTER                 :: potential
      TYPE(atom_type), POINTER                           :: atom
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      NULLIFY (atom)
      CALL create_atom_type(atom)

      IF (PRESENT(xc_section)) THEN
         atom%xc_section => xc_section
      ELSE
         NULLIFY (atom%xc_section)
      END IF

      CALL get_atomic_kind(atomic_kind, z=z)
      NULLIFY (all_potential, gth_potential, sgp_potential, orb_basis_set)
      CALL get_qs_kind(qs_kind, zeff=zeff, &
                       basis_set=orb_basis_set, &
                       all_potential=all_potential, &
                       gth_potential=gth_potential, &
                       sgp_potential=sgp_potential)

      IF (ASSOCIATED(orb_basis_set)) THEN
         has_basis = .TRUE.
      ELSE
         has_basis = .FALSE.
      END IF

      atom%z = z
      CALL set_atom(atom, &
                    pp_calc=(ASSOCIATED(gth_potential) .OR. ASSOCIATED(sgp_potential)), &
                    do_zmp=.FALSE., &
                    doread=.FALSE., &
                    read_vxc=.FALSE., &
                    relativistic=do_nonrel_atom, &
                    coulomb_integral_type=do_numeric, &
                    exchange_integral_type=do_numeric)

      ALLOCATE (potential, basis, integrals)

      IF (PRESENT(confine)) THEN
         potential%confinement = confine
      ELSE
         IF (ASSOCIATED(gth_potential) .OR. ASSOCIATED(sgp_potential)) THEN
            potential%confinement = .TRUE.
         ELSE
            potential%confinement = .FALSE.
         END IF
      END IF
      potential%conf_type = poly_conf
      potential%acon = 0.1_dp
      potential%rcon = 2.0_dp*ptable(z)%vdw_radius*bohr
      potential%scon = 2.0_dp

      IF (ASSOCIATED(gth_potential)) THEN
         potential%ppot_type = gth_pseudo
         CALL get_potential(gth_potential, zeff=zeff)
         CALL gth_potential_conversion(gth_potential, potential%gth_pot)
         CALL set_atom(atom, zcore=NINT(zeff), potential=potential)
      ELSE IF (ASSOCIATED(sgp_potential)) THEN
         potential%ppot_type = sgp_pseudo
         CALL sgp_potential_conversion(sgp_potential, potential%sgp_pot)
         potential%sgp_pot%symbol = ptable(z)%symbol
         CALL get_potential(sgp_potential, zeff=zeff)
         CALL set_atom(atom, zcore=NINT(zeff), potential=potential)
      ELSE
         potential%ppot_type = no_pseudo
         CALL set_atom(atom, zcore=z, potential=potential)
      END IF

      NULLIFY (grid)
      IF (PRESENT(agrid)) THEN
         ngp = agrid%nr
         quadtype = agrid%quadrature
      ELSE
         ngp = 400
         quadtype = do_gapw_log
      END IF
      CALL allocate_grid_atom(grid)
      CALL create_grid_atom(grid, ngp, 1, 1, 0, quadtype)
      grid%nr = ngp
      basis%grid => grid

      NULLIFY (basis%am, basis%cm, basis%as, basis%ns, basis%bf, basis%dbf, basis%ddbf)

      IF (has_basis) THEN
         ! fill in the basis data structures
         basis%basis_type = CGTO_BASIS
         basis%eps_eig = 1.e-12_dp
         CALL get_gto_basis_set(orb_basis_set, &
                                nset=nset, nshell=nshell, npgf=npgf, lmin=lmin, lmax=lmax, &
                                l=ls, nsgf=nsgf, zet=zet, gcc=gcc, &
                                first_sgf=first_sgf, last_sgf=last_sgf)
         set_index = 0
         shell_index = 0
         basis%nprim = 0
         basis%nbas = 0
         DO i = 1, nset
            DO j = lmin(i), MIN(lmax(i), lmat)
               basis%nprim(j) = basis%nprim(j) + npgf(i)
            END DO
            DO j = 1, nshell(i)
               l = ls(j, i)
               IF (l <= lmat) THEN
                  basis%nbas(l) = basis%nbas(l) + 1
                  k = basis%nbas(l)
                  CPASSERT(k <= 100)
                  set_index(l, k) = i
                  shell_index(l, k) = j
               END IF
            END DO
         END DO

         nj = MAXVAL(basis%nprim)
         ns = MAXVAL(basis%nbas)
         ALLOCATE (basis%am(nj, 0:lmat))
         basis%am = 0._dp
         ALLOCATE (basis%cm(nj, ns, 0:lmat))
         basis%cm = 0._dp
         DO j = 0, lmat
            nj = 0
            ns = 0
            DO i = 1, nset
               IF (j >= lmin(i) .AND. j <= lmax(i)) THEN
                  DO ipgf = 1, npgf(i)
                     basis%am(nj + ipgf, j) = zet(ipgf, i)
                  END DO
                  DO ii = 1, nshell(i)
                     IF (ls(ii, i) == j) THEN
                        ns = ns + 1
                        DO ipgf = 1, npgf(i)
                           basis%cm(nj + ipgf, ns, j) = gcc(ipgf, ii, i)
                        END DO
                     END IF
                  END DO
                  nj = nj + npgf(i)
               END IF
            END DO
         END DO
      ELSE
         ! use default basis
         IF (ASSOCIATED(gth_potential) .OR. ASSOCIATED(sgp_potential)) THEN
            nu = nup
         ELSE
            nu = nua
         END IF
         basis%geometrical = .FALSE.
         basis%aval = 0._dp
         basis%cval = 0._dp
         basis%start = 0
         basis%eps_eig = 1.e-12_dp

         basis%basis_type = GTO_BASIS
         basis%nbas = nu
         basis%nprim = nu
         ALLOCATE (basis%am(nu, 0:lmat))
         DO i = 0, lmat
            basis%am(1:nu, i) = ugbs(1:nu)
         END DO
      END IF

      ! initialize basis function on a radial grid
      nr = basis%grid%nr
      m = MAXVAL(basis%nbas)
      ALLOCATE (basis%bf(nr, m, 0:lmat))
      ALLOCATE (basis%dbf(nr, m, 0:lmat))
      ALLOCATE (basis%ddbf(nr, m, 0:lmat))

      basis%bf = 0._dp
      basis%dbf = 0._dp
      basis%ddbf = 0._dp
      DO l = 0, lmat
         DO i = 1, basis%nprim(l)
            al = basis%am(i, l)
            IF (basis%basis_type == GTO_BASIS) THEN
               DO k = 1, nr
                  rk = basis%grid%rad(k)
                  ear = EXP(-al*basis%grid%rad(k)**2)
                  basis%bf(k, i, l) = rk**l*ear
                  basis%dbf(k, i, l) = (REAL(l, dp)*rk**(l - 1) - 2._dp*al*rk**(l + 1))*ear
                  basis%ddbf(k, i, l) = (REAL(l*(l - 1), dp)*rk**(l - 2) - &
                                         2._dp*al*REAL(2*l + 1, dp)*rk**(l) + 4._dp*al*rk**(l + 2))*ear
               END DO
            ELSEIF (basis%basis_type == CGTO_BASIS) THEN
               DO k = 1, nr
                  rk = basis%grid%rad(k)
                  ear = EXP(-al*basis%grid%rad(k)**2)
                  DO j = 1, basis%nbas(l)
                     basis%bf(k, j, l) = basis%bf(k, j, l) + rk**l*ear*basis%cm(i, j, l)
                     basis%dbf(k, j, l) = basis%dbf(k, j, l) &
                                          + (REAL(l, dp)*rk**(l - 1) - 2._dp*al*rk**(l + 1))*ear*basis%cm(i, j, l)
                     basis%ddbf(k, j, l) = basis%ddbf(k, j, l) + &
                                           (REAL(l*(l - 1), dp)*rk**(l - 2) - 2._dp*al*REAL(2*l + 1, dp)*rk**(l) + &
                                            4._dp*al*rk**(l + 2))*ear*basis%cm(i, j, l)
                  END DO
               END DO
            ELSE
               CPABORT('Atom basis type?')
            END IF
         END DO
      END DO

      CALL set_atom(atom, basis=basis)

      ! optimization defaults
      atom%optimization%damping = 0.2_dp
      atom%optimization%eps_scf = 1.e-6_dp
      atom%optimization%eps_diis = 100._dp
      atom%optimization%max_iter = 50
      atom%optimization%n_diis = 5

      ! set up the electronic state
      CALL init_atom_electronic_state(atomic_kind=atomic_kind, &
                                      qs_kind=qs_kind, &
                                      ncalc=ncalc, &
                                      ncore=ncore, &
                                      nelem=nelem, &
                                      edelta=edelta)

      ! restricted or unrestricted?
      IF (SUM(ABS(edelta)) > 0.0_dp) THEN
         uks = .TRUE.
         CALL set_atom(atom, method_type=do_uks_atom)
      ELSE
         uks = .FALSE.
         CALL set_atom(atom, method_type=do_rks_atom)
      END IF

      ALLOCATE (atom%state)

      atom%state%core = 0._dp
      atom%state%core(0:lmat, 1:7) = REAL(ncore(0:lmat, 1:7), dp)
      atom%state%occ = 0._dp
      IF (uks) THEN
         atom%state%occ(0:lmat, 1:7) = REAL(ncalc(0:lmat, 1:7), dp) + &
                                       edelta(0:lmat, 1:7, 1) + edelta(0:lmat, 1:7, 2)
      ELSE
         atom%state%occ(0:lmat, 1:7) = REAL(ncalc(0:lmat, 1:7), dp)
      END IF
      atom%state%occupation = 0._dp
      DO l = 0, lmat
         k = 0
         DO i = 1, 7
            IF (ncalc(l, i) > 0) THEN
               k = k + 1
               IF (uks) THEN
                  atom%state%occupation(l, k) = REAL(ncalc(l, i), dp) + &
                                                edelta(l, i, 1) + edelta(l, i, 2)
                  atom%state%occa(l, k) = 0.5_dp*REAL(ncalc(l, i), dp) + edelta(l, i, 1)
                  atom%state%occb(l, k) = 0.5_dp*REAL(ncalc(l, i), dp) + edelta(l, i, 2)
               ELSE
                  atom%state%occupation(l, k) = REAL(ncalc(l, i), dp)
               END IF
            END IF
         END DO
         ok = REAL(2*l + 1, KIND=dp)
         IF (uks) THEN
            DO i = 1, 7
               atom%state%occ(l, i) = MIN(atom%state%occ(l, i), 2.0_dp*ok)
               atom%state%occa(l, i) = MIN(atom%state%occa(l, i), ok)
               atom%state%occb(l, i) = MIN(atom%state%occb(l, i), ok)
               atom%state%occupation(l, i) = atom%state%occa(l, i) + atom%state%occb(l, i)
            END DO
         ELSE
            DO i = 1, 7
               atom%state%occ(l, i) = MIN(atom%state%occ(l, i), 2.0_dp*ok)
               atom%state%occupation(l, i) = MIN(atom%state%occupation(l, i), 2.0_dp*ok)
            END DO
         END IF
      END DO
      IF (uks) THEN
         atom%state%multiplicity = NINT(ABS(SUM(atom%state%occa - atom%state%occb)) + 1)
      ELSE
         atom%state%multiplicity = -1
      END IF

      atom%state%maxl_occ = get_maxl_occ(atom%state%occupation)
      atom%state%maxn_occ = get_maxn_occ(atom%state%occupation)
      atom%state%maxl_calc = atom%state%maxl_occ
      atom%state%maxn_calc = atom%state%maxn_occ

      ! total number of occupied orbitals
      IF (PRESENT(nocc)) THEN
         nocc = 0
         DO l = 0, lmat
            DO k = 1, 7
               IF (uks) THEN
                  IF (atom%state%occa(l, k) > 0.0_dp) THEN
                     nocc(1) = nocc(1) + 2*l + 1
                  END IF
                  IF (atom%state%occb(l, k) > 0.0_dp) THEN
                     nocc(2) = nocc(2) + 2*l + 1
                  END IF
               ELSE
                  IF (atom%state%occupation(l, k) > 0.0_dp) THEN
                     nocc(1) = nocc(1) + 2*l + 1
                     nocc(2) = nocc(2) + 2*l + 1
                  END IF
               END IF
            END DO
         END DO
      END IF

      ! calculate integrals
      ! general integrals
      CALL atom_int_setup(integrals, basis, potential=atom%potential, &
                          eri_coulomb=(atom%coulomb_integral_type == do_analytic), &
                          eri_exchange=(atom%exchange_integral_type == do_analytic))
      ! potential
      CALL atom_ppint_setup(integrals, basis, potential=atom%potential)
      ! relativistic correction terms
      NULLIFY (integrals%tzora, integrals%hdkh)
      CALL atom_relint_setup(integrals, basis, atom%relativistic, zcore=REAL(atom%zcore, dp))
      CALL set_atom(atom, integrals=integrals)

      NULLIFY (orbitals)
      mo = MAXVAL(atom%state%maxn_calc)
      mb = MAXVAL(atom%basis%nbas)
      CALL create_atom_orbs(orbitals, mb, mo)
      CALL set_atom(atom, orbitals=orbitals)

      IF (PRESENT(iunit)) THEN
         CALL calculate_atom(atom, iunit)
      ELSE
         CALL calculate_atom(atom, -1)
      END IF

      IF (PRESENT(pmat)) THEN
         ! recover density matrix in CP2K/GPW order and normalization
         IF (ASSOCIATED(pmat)) THEN
            DEALLOCATE (pmat)
         END IF
         ALLOCATE (pmat(nsgf, nsgf, 2))
         pmat = 0._dp
         DO l = 0, lmat
            ll = 2*l
            DO k1 = 1, atom%basis%nbas(l)
               DO k2 = 1, atom%basis%nbas(l)
                  scal = SQRT(atom%integrals%ovlp(k1, k1, l)*atom%integrals%ovlp(k2, k2, l))/REAL(2*l + 1, KIND=dp)
                  i = first_sgf(shell_index(l, k1), set_index(l, k1))
                  j = first_sgf(shell_index(l, k2), set_index(l, k2))
                  IF (uks) THEN
                     DO m = 0, ll
                        pmat(i + m, j + m, 1) = atom%orbitals%pmata(k1, k2, l)*scal
                        pmat(i + m, j + m, 2) = atom%orbitals%pmatb(k1, k2, l)*scal
                     END DO
                  ELSE
                     DO m = 0, ll
                        pmat(i + m, j + m, 1) = atom%orbitals%pmat(k1, k2, l)*scal
                     END DO
                  END IF
               END DO
            END DO
         ENDDO
         IF (uks) THEN
            pmat(:, :, 1) = pmat(:, :, 1) + pmat(:, :, 2)
            pmat(:, :, 2) = pmat(:, :, 1) - 2.0_dp*pmat(:, :, 2)
         END IF
      END IF

      IF (PRESENT(fmat)) THEN
         ! recover fock matrix in CP2K/GPW order.
         ! Caution: Normalization is not take care of, so it's probably weird.
         IF (uks) CPABORT("calculate_atomic_orbitals: only RKS is implemented")
         IF (ASSOCIATED(fmat)) CPABORT("fmat already associated")
         IF (.NOT. ASSOCIATED(atom%fmat)) CPABORT("atom%fmat not associated")
         ALLOCATE (fmat(nsgf, nsgf, 1))
         fmat = 0.0_dp
         DO l = 0, lmat
            ll = 2*l
            DO k1 = 1, atom%basis%nbas(l)
            DO k2 = 1, atom%basis%nbas(l)
               scal = SQRT(atom%integrals%ovlp(k1, k1, l)*atom%integrals%ovlp(k2, k2, l))
               i = first_sgf(shell_index(l, k1), set_index(l, k1))
               j = first_sgf(shell_index(l, k2), set_index(l, k2))
               DO m = 0, ll
                  fmat(i + m, j + m, 1) = atom%fmat%op(k1, k2, l)/scal
               END DO
            END DO
            END DO
         ENDDO
      END IF

      IF (PRESENT(density)) THEN
         IF (ASSOCIATED(density)) DEALLOCATE (density)
         ALLOCATE (density(nr))
         CALL atom_density(density, atom%orbitals%pmat, atom%basis, atom%state%maxl_occ)
      END IF

      IF (PRESENT(wavefunction)) THEN
         CPASSERT(PRESENT(wfninfo))
         IF (ASSOCIATED(wavefunction)) DEALLOCATE (wavefunction)
         IF (ASSOCIATED(wfninfo)) DEALLOCATE (wfninfo)
         mo = SUM(atom%state%maxn_occ)
         ALLOCATE (wavefunction(nr, mo), wfninfo(2, mo))
         wavefunction = 0.0_dp
         ii = 0
         DO l = 0, lmat
            DO i = 1, atom%state%maxn_occ(l)
               IF (atom%state%occupation(l, i) > 0.0_dp) THEN
                  ii = ii + 1
                  wfninfo(1, ii) = atom%state%occupation(l, i)
                  wfninfo(2, ii) = REAL(l, dp)
                  DO j = 1, atom%basis%nbas(l)
                     wavefunction(:, ii) = wavefunction(:, ii) + &
                                           atom%orbitals%wfn(j, i, l)*basis%bf(:, j, l)
                  END DO
               END IF
            END DO
         END DO
         CPASSERT(mo == ii)
      END IF

      ! clean up
      CALL atom_int_release(integrals)
      CALL atom_ppint_release(integrals)
      CALL atom_relint_release(integrals)
      CALL release_atom_basis(basis)
      CALL release_atom_potential(potential)
      CALL release_atom_type(atom)

      DEALLOCATE (potential, basis, integrals)

   END SUBROUTINE calculate_atomic_orbitals

! **************************************************************************************************
!> \brief ...
!> \param density ...
!> \param atomic_kind ...
!> \param qs_kind ...
!> \param ngto ...
!> \param iunit ...
!> \param allelectron ...
!> \param confine ...
! **************************************************************************************************
   SUBROUTINE calculate_atomic_density(density, atomic_kind, qs_kind, ngto, iunit, allelectron, confine)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: density
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      INTEGER, INTENT(IN)                                :: ngto
      INTEGER, INTENT(IN), OPTIONAL                      :: iunit
      LOGICAL, INTENT(IN), OPTIONAL                      :: allelectron, confine

      INTEGER, PARAMETER                                 :: num_gto = 40

      INTEGER                                            :: i, ii, k, l, ll, m, mb, mo, ngp, nn, nr, &
                                                            quadtype, relativistic, z
      INTEGER, DIMENSION(0:lmat)                         :: starti
      INTEGER, DIMENSION(0:lmat, 10)                     :: ncalc, ncore, nelem
      INTEGER, DIMENSION(:), POINTER                     :: econf
      REAL(KIND=dp)                                      :: al, aval, cc, cval, ear, rk, xx, zeff
      REAL(KIND=dp), DIMENSION(num_gto+2)                :: results
      TYPE(all_potential_type), POINTER                  :: all_potential
      TYPE(atom_basis_type), POINTER                     :: basis
      TYPE(atom_integrals), POINTER                      :: integrals
      TYPE(atom_orbitals), POINTER                       :: orbitals
      TYPE(atom_potential_type), POINTER                 :: potential
      TYPE(atom_type), POINTER                           :: atom
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      NULLIFY (atom)
      CALL create_atom_type(atom)

      CALL get_atomic_kind(atomic_kind, z=z)
      NULLIFY (all_potential, gth_potential)
      CALL get_qs_kind(qs_kind, zeff=zeff, &
                       all_potential=all_potential, &
                       gth_potential=gth_potential, &
                       sgp_potential=sgp_potential)

      IF (PRESENT(allelectron)) THEN
         IF (allelectron) THEN
            NULLIFY (gth_potential)
            zeff = z
         END IF
      END IF

      CPASSERT(ngto <= num_gto)

      IF (ASSOCIATED(gth_potential) .OR. ASSOCIATED(sgp_potential)) THEN
         ! PP calculation are non-relativistic
         relativistic = do_nonrel_atom
      ELSE
         ! AE calculations use DKH2
         relativistic = do_dkh2_atom
      END IF

      atom%z = z
      CALL set_atom(atom, &
                    pp_calc=(ASSOCIATED(gth_potential) .OR. ASSOCIATED(sgp_potential)), &
                    method_type=do_rks_atom, &
                    relativistic=relativistic, &
                    coulomb_integral_type=do_numeric, &
                    exchange_integral_type=do_numeric)

      ALLOCATE (potential, basis, integrals)

      IF (PRESENT(confine)) THEN
         potential%confinement = confine
      ELSE
         IF (ASSOCIATED(gth_potential) .OR. ASSOCIATED(sgp_potential)) THEN
            potential%confinement = .TRUE.
         ELSE
            potential%confinement = .FALSE.
         END IF
      END IF
      potential%conf_type = barrier_conf
      potential%acon = 200._dp
      potential%rcon = 4.0_dp
      potential%scon = 8.0_dp

      IF (ASSOCIATED(gth_potential)) THEN
         potential%ppot_type = gth_pseudo
         CALL get_potential(gth_potential, zeff=zeff)
         CALL gth_potential_conversion(gth_potential, potential%gth_pot)
         CALL set_atom(atom, zcore=NINT(zeff), potential=potential)
      ELSE IF (ASSOCIATED(sgp_potential)) THEN
         potential%ppot_type = sgp_pseudo
         CALL sgp_potential_conversion(sgp_potential, potential%sgp_pot)
         potential%sgp_pot%symbol = ptable(z)%symbol
         CALL get_potential(sgp_potential, zeff=zeff)
         CALL set_atom(atom, zcore=NINT(zeff), potential=potential)
      ELSE
         potential%ppot_type = no_pseudo
         CALL set_atom(atom, zcore=z, potential=potential)
      END IF

      ! atomic grid
      NULLIFY (grid)
      ngp = 400
      quadtype = do_gapw_log
      CALL allocate_grid_atom(grid)
      CALL create_grid_atom(grid, ngp, 1, 1, 0, quadtype)
      grid%nr = ngp
      basis%grid => grid

      NULLIFY (basis%am, basis%cm, basis%as, basis%ns, basis%bf, basis%dbf, basis%ddbf)

      ! fill in the basis data structures
      basis%eps_eig = 1.e-12_dp
      basis%basis_type = GTO_BASIS
      CALL Clementi_geobas(z, cval, aval, basis%nbas, starti)
      basis%nprim = basis%nbas
      m = MAXVAL(basis%nbas)
      ALLOCATE (basis%am(m, 0:lmat))
      basis%am = 0._dp
      DO l = 0, lmat
         DO i = 1, basis%nbas(l)
            ll = i - 1 + starti(l)
            basis%am(i, l) = aval*cval**(ll)
         END DO
      END DO

      basis%geometrical = .TRUE.
      basis%aval = aval
      basis%cval = cval
      basis%start = starti

      ! initialize basis function on a radial grid
      nr = basis%grid%nr
      m = MAXVAL(basis%nbas)
      ALLOCATE (basis%bf(nr, m, 0:lmat))
      ALLOCATE (basis%dbf(nr, m, 0:lmat))
      ALLOCATE (basis%ddbf(nr, m, 0:lmat))
      basis%bf = 0._dp
      basis%dbf = 0._dp
      basis%ddbf = 0._dp
      DO l = 0, lmat
         DO i = 1, basis%nbas(l)
            al = basis%am(i, l)
            DO k = 1, nr
               rk = basis%grid%rad(k)
               ear = EXP(-al*basis%grid%rad(k)**2)
               basis%bf(k, i, l) = rk**l*ear
               basis%dbf(k, i, l) = (REAL(l, dp)*rk**(l - 1) - 2._dp*al*rk**(l + 1))*ear
               basis%ddbf(k, i, l) = (REAL(l*(l - 1), dp)*rk**(l - 2) - &
                                      2._dp*al*REAL(2*l + 1, dp)*rk**(l) + 4._dp*al*rk**(l + 2))*ear
            END DO
         END DO
      END DO

      CALL set_atom(atom, basis=basis)

      ! optimization defaults
      atom%optimization%damping = 0.2_dp
      atom%optimization%eps_scf = 1.e-6_dp
      atom%optimization%eps_diis = 100._dp
      atom%optimization%max_iter = 50
      atom%optimization%n_diis = 5

      nelem = 0
      ncore = 0
      ncalc = 0
      IF (ASSOCIATED(gth_potential)) THEN
         CALL get_potential(gth_potential, elec_conf=econf)
         CALL set_pseudo_state(econf, z, ncalc, ncore, nelem)
      ELSE IF (ASSOCIATED(sgp_potential)) THEN
         CALL get_potential(sgp_potential, elec_conf=econf)
         CALL set_pseudo_state(econf, z, ncalc, ncore, nelem)
      ELSE
         DO l = 0, MIN(lmat, UBOUND(ptable(z)%e_conv, 1))
            ll = 2*(2*l + 1)
            nn = ptable(z)%e_conv(l)
            ii = 0
            DO
               ii = ii + 1
               IF (nn <= ll) THEN
                  nelem(l, ii) = nn
                  EXIT
               ELSE
                  nelem(l, ii) = ll
                  nn = nn - ll
               END IF
            END DO
         END DO
         ncalc = nelem - ncore
      END IF

      IF (qs_kind%ghost .OR. qs_kind%floating) THEN
         nelem = 0
         ncore = 0
         ncalc = 0
      END IF

      ALLOCATE (atom%state)

      atom%state%core = 0._dp
      atom%state%core(0:lmat, 1:7) = REAL(ncore(0:lmat, 1:7), dp)
      atom%state%occ = 0._dp
      atom%state%occ(0:lmat, 1:7) = REAL(ncalc(0:lmat, 1:7), dp)
      atom%state%occupation = 0._dp
      atom%state%multiplicity = -1
      DO l = 0, lmat
         k = 0
         DO i = 1, 7
            IF (ncalc(l, i) > 0) THEN
               k = k + 1
               atom%state%occupation(l, k) = REAL(ncalc(l, i), dp)
            END IF
         END DO
      END DO

      atom%state%maxl_occ = get_maxl_occ(atom%state%occupation)
      atom%state%maxn_occ = get_maxn_occ(atom%state%occupation)
      atom%state%maxl_calc = atom%state%maxl_occ
      atom%state%maxn_calc = atom%state%maxn_occ

      ! calculate integrals
      ! general integrals
      CALL atom_int_setup(integrals, basis, potential=atom%potential, &
                          eri_coulomb=(atom%coulomb_integral_type == do_analytic), &
                          eri_exchange=(atom%exchange_integral_type == do_analytic))
      ! potential
      CALL atom_ppint_setup(integrals, basis, potential=atom%potential)
      ! relativistic correction terms
      NULLIFY (integrals%tzora, integrals%hdkh)
      CALL atom_relint_setup(integrals, basis, atom%relativistic, zcore=REAL(atom%zcore, dp))
      CALL set_atom(atom, integrals=integrals)

      NULLIFY (orbitals)
      mo = MAXVAL(atom%state%maxn_calc)
      mb = MAXVAL(atom%basis%nbas)
      CALL create_atom_orbs(orbitals, mb, mo)
      CALL set_atom(atom, orbitals=orbitals)

      IF (PRESENT(iunit)) THEN
         CALL calculate_atom(atom, iunit)
         CALL atom_fit_density(atom, ngto, 0, iunit, results=results)
      ELSE
         CALL calculate_atom(atom, -1)
         CALL atom_fit_density(atom, ngto, 0, -1, results=results)
      END IF

      xx = results(1)
      cc = results(2)
      DO i = 1, ngto
         density(i, 1) = xx*cc**i
         density(i, 2) = results(2 + i)
      END DO

      ! clean up
      CALL atom_int_release(integrals)
      CALL atom_ppint_release(integrals)
      CALL atom_relint_release(integrals)
      CALL release_atom_basis(basis)
      CALL release_atom_potential(potential)
      CALL release_atom_type(atom)

      DEALLOCATE (potential, basis, integrals)

   END SUBROUTINE calculate_atomic_density

! **************************************************************************************************
!> \brief ...
!> \param atomic_kind ...
!> \param qs_kind ...
!> \param rel_control ...
!> \param rtmat ...
! **************************************************************************************************
   SUBROUTINE calculate_atomic_relkin(atomic_kind, qs_kind, rel_control, rtmat)
      TYPE(atomic_kind_type), INTENT(IN)                 :: atomic_kind
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind
      TYPE(rel_control_type), POINTER                    :: rel_control
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rtmat

      INTEGER                                            :: i, ii, ipgf, j, k, k1, k2, l, ll, m, n, &
                                                            ngp, nj, nn, nr, ns, nset, nsgf, &
                                                            quadtype, relativistic, z
      INTEGER, DIMENSION(0:lmat, 10)                     :: ncalc, ncore, nelem
      INTEGER, DIMENSION(0:lmat, 100)                    :: set_index, shell_index
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf, last_sgf, ls
      REAL(KIND=dp)                                      :: al, alpha, ear, prefac, rk, zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: omat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: gcc
      TYPE(all_potential_type), POINTER                  :: all_potential
      TYPE(atom_basis_type), POINTER                     :: basis
      TYPE(atom_integrals), POINTER                      :: integrals
      TYPE(atom_potential_type), POINTER                 :: potential
      TYPE(atom_type), POINTER                           :: atom
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      IF (rel_control%rel_method == rel_none) RETURN

      NULLIFY (all_potential, orb_basis_set)
      CALL get_qs_kind(qs_kind, basis_set=orb_basis_set, all_potential=all_potential)

      CPASSERT(ASSOCIATED(orb_basis_set))

      IF (ASSOCIATED(all_potential)) THEN
         ! only all electron atoms will get the relativistic correction

         CALL get_atomic_kind(atomic_kind, z=z)
         CALL get_qs_kind(qs_kind, zeff=zeff)
         NULLIFY (atom)
         CALL create_atom_type(atom)
         NULLIFY (atom%xc_section)
         NULLIFY (atom%orbitals)
         atom%z = z
         alpha = SQRT(all_potential%alpha_core_charge)

         ! set the method flag
         SELECT CASE (rel_control%rel_method)
         CASE DEFAULT
            CPABORT("")
         CASE (rel_dkh)
            SELECT CASE (rel_control%rel_DKH_order)
            CASE DEFAULT
               CPABORT("")
            CASE (0)
               relativistic = do_dkh0_atom
            CASE (1)
               relativistic = do_dkh1_atom
            CASE (2)
               relativistic = do_dkh2_atom
            CASE (3)
               relativistic = do_dkh3_atom
            END SELECT
         CASE (rel_zora)
            SELECT CASE (rel_control%rel_zora_type)
            CASE DEFAULT
               CPABORT("")
            CASE (rel_zora_full)
               CPABORT("")
            CASE (rel_zora_mp)
               relativistic = do_zoramp_atom
            CASE (rel_sczora_mp)
               relativistic = do_sczoramp_atom
            END SELECT
         END SELECT

         CALL set_atom(atom, &
                       pp_calc=.FALSE., &
                       method_type=do_rks_atom, &
                       relativistic=relativistic, &
                       coulomb_integral_type=do_numeric, &
                       exchange_integral_type=do_numeric)

         ALLOCATE (potential, basis, integrals)

         potential%ppot_type = no_pseudo
         CALL set_atom(atom, zcore=z, potential=potential)

         CALL get_gto_basis_set(orb_basis_set, &
                                nset=nset, nshell=nshell, npgf=npgf, lmin=lmin, lmax=lmax, l=ls, nsgf=nsgf, zet=zet, gcc=gcc, &
                                first_sgf=first_sgf, last_sgf=last_sgf)

         NULLIFY (grid)
         ngp = 400
         quadtype = do_gapw_log
         CALL allocate_grid_atom(grid)
         CALL create_grid_atom(grid, ngp, 1, 1, 0, quadtype)
         grid%nr = ngp
         basis%grid => grid

         NULLIFY (basis%am, basis%cm, basis%as, basis%ns, basis%bf, basis%dbf, basis%ddbf)
         basis%basis_type = CGTO_BASIS
         basis%eps_eig = 1.e-12_dp

         ! fill in the basis data structures
         set_index = 0
         shell_index = 0
         basis%nprim = 0
         basis%nbas = 0
         DO i = 1, nset
            DO j = lmin(i), MIN(lmax(i), lmat)
               basis%nprim(j) = basis%nprim(j) + npgf(i)
            END DO
            DO j = 1, nshell(i)
               l = ls(j, i)
               IF (l <= lmat) THEN
                  basis%nbas(l) = basis%nbas(l) + 1
                  k = basis%nbas(l)
                  CPASSERT(k <= 100)
                  set_index(l, k) = i
                  shell_index(l, k) = j
               END IF
            END DO
         END DO

         nj = MAXVAL(basis%nprim)
         ns = MAXVAL(basis%nbas)
         ALLOCATE (basis%am(nj, 0:lmat))
         basis%am = 0._dp
         ALLOCATE (basis%cm(nj, ns, 0:lmat))
         basis%cm = 0._dp
         DO j = 0, lmat
            nj = 0
            ns = 0
            DO i = 1, nset
               IF (j >= lmin(i) .AND. j <= lmax(i)) THEN
                  DO ipgf = 1, npgf(i)
                     basis%am(nj + ipgf, j) = zet(ipgf, i)
                  END DO
                  DO ii = 1, nshell(i)
                     IF (ls(ii, i) == j) THEN
                        ns = ns + 1
                        DO ipgf = 1, npgf(i)
                           basis%cm(nj + ipgf, ns, j) = gcc(ipgf, ii, i)
                        END DO
                     END IF
                  END DO
                  nj = nj + npgf(i)
               END IF
            END DO
         END DO

         ! Normalization as used in the atomic code
         ! We have to undo the Quickstep normalization
         DO j = 0, lmat
            prefac = 2.0_dp*SQRT(pi/dfac(2*j + 1))
            DO ipgf = 1, basis%nprim(j)
               DO ii = 1, basis%nbas(j)
                  basis%cm(ipgf, ii, j) = prefac*basis%cm(ipgf, ii, j)
               END DO
            END DO
         END DO

         ! initialize basis function on a radial grid
         nr = basis%grid%nr
         m = MAXVAL(basis%nbas)
         ALLOCATE (basis%bf(nr, m, 0:lmat))
         ALLOCATE (basis%dbf(nr, m, 0:lmat))
         ALLOCATE (basis%ddbf(nr, m, 0:lmat))

         basis%bf = 0._dp
         basis%dbf = 0._dp
         basis%ddbf = 0._dp
         DO l = 0, lmat
            DO i = 1, basis%nprim(l)
               al = basis%am(i, l)
               DO k = 1, nr
                  rk = basis%grid%rad(k)
                  ear = EXP(-al*basis%grid%rad(k)**2)
                  DO j = 1, basis%nbas(l)
                     basis%bf(k, j, l) = basis%bf(k, j, l) + rk**l*ear*basis%cm(i, j, l)
                     basis%dbf(k, j, l) = basis%dbf(k, j, l) &
                                          + (REAL(l, dp)*rk**(l - 1) - 2._dp*al*rk**(l + 1))*ear*basis%cm(i, j, l)
                     basis%ddbf(k, j, l) = basis%ddbf(k, j, l) + &
                                           (REAL(l*(l - 1), dp)*rk**(l - 2) - 2._dp*al*REAL(2*l + 1, dp)* &
                                            rk**(l) + 4._dp*al*rk**(l + 2))*ear*basis%cm(i, j, l)
                  END DO
               END DO
            END DO
         END DO

         CALL set_atom(atom, basis=basis)

         ! optimization defaults
         atom%optimization%damping = 0.2_dp
         atom%optimization%eps_scf = 1.e-6_dp
         atom%optimization%eps_diis = 100._dp
         atom%optimization%max_iter = 50
         atom%optimization%n_diis = 5

         ! electronic state
         nelem = 0
         ncore = 0
         ncalc = 0
         DO l = 0, MIN(lmat, UBOUND(ptable(z)%e_conv, 1))
            ll = 2*(2*l + 1)
            nn = ptable(z)%e_conv(l)
            ii = 0
            DO
               ii = ii + 1
               IF (nn <= ll) THEN
                  nelem(l, ii) = nn
                  EXIT
               ELSE
                  nelem(l, ii) = ll
                  nn = nn - ll
               END IF
            END DO
         END DO
         ncalc = nelem - ncore

         IF (qs_kind%ghost .OR. qs_kind%floating) THEN
            nelem = 0
            ncore = 0
            ncalc = 0
         END IF

         ALLOCATE (atom%state)

         atom%state%core = 0._dp
         atom%state%core(0:lmat, 1:7) = REAL(ncore(0:lmat, 1:7), dp)
         atom%state%occ = 0._dp
         atom%state%occ(0:lmat, 1:7) = REAL(ncalc(0:lmat, 1:7), dp)
         atom%state%occupation = 0._dp
         atom%state%multiplicity = -1
         DO l = 0, lmat
            k = 0
            DO i = 1, 7
               IF (ncalc(l, i) > 0) THEN
                  k = k + 1
                  atom%state%occupation(l, k) = REAL(ncalc(l, i), dp)
               END IF
            END DO
         END DO

         atom%state%maxl_occ = get_maxl_occ(atom%state%occupation)
         atom%state%maxn_occ = get_maxn_occ(atom%state%occupation)
         atom%state%maxl_calc = atom%state%maxl_occ
         atom%state%maxn_calc = atom%state%maxn_occ

         ! calculate integrals
         ! general integrals
         CALL atom_int_setup(integrals, basis)
         ! potential
         CALL atom_ppint_setup(integrals, basis, potential=atom%potential)
         ! relativistic correction terms
         NULLIFY (integrals%tzora, integrals%hdkh)
         CALL atom_relint_setup(integrals, basis, atom%relativistic, zcore=REAL(atom%zcore, dp), &
                                alpha=alpha)
         CALL set_atom(atom, integrals=integrals)

         ! for DKH we need erfc integrals to correct non-relativistic
         integrals%core = 0.0_dp
         DO l = 0, lmat
            n = integrals%n(l)
            m = basis%nprim(l)
            ALLOCATE (omat(m, m))

            CALL sg_erfc(omat(1:m, 1:m), l, alpha, basis%am(1:m, l), basis%am(1:m, l))
            integrals%core(1:n, 1:n, l) = MATMUL(TRANSPOSE(basis%cm(1:m, 1:n, l)), &
                                                 MATMUL(omat(1:m, 1:m), basis%cm(1:m, 1:n, l)))

            DEALLOCATE (omat)
         END DO

         ! recover relativistic kinetic matrix in CP2K/GPW order and normalization
         IF (ASSOCIATED(rtmat)) THEN
            DEALLOCATE (rtmat)
         END IF
         ALLOCATE (rtmat(nsgf, nsgf))
         rtmat = 0._dp
         DO l = 0, lmat
            ll = 2*l
            DO k1 = 1, basis%nbas(l)
               DO k2 = 1, basis%nbas(l)
                  i = first_sgf(shell_index(l, k1), set_index(l, k1))
                  j = first_sgf(shell_index(l, k2), set_index(l, k2))
                  SELECT CASE (atom%relativistic)
                  CASE DEFAULT
                     CPABORT("")
                  CASE (do_zoramp_atom, do_sczoramp_atom)
                     DO m = 0, ll
                        rtmat(i + m, j + m) = integrals%tzora(k1, k2, l)
                     END DO
                  CASE (do_dkh0_atom, do_dkh1_atom, do_dkh2_atom, do_dkh3_atom)
                     DO m = 0, ll
                        rtmat(i + m, j + m) = integrals%hdkh(k1, k2, l) - integrals%kin(k1, k2, l) + &
                                              atom%zcore*integrals%core(k1, k2, l)
                     END DO
                  END SELECT
               END DO
            END DO
         ENDDO
         DO k1 = 1, nsgf
            DO k2 = k1, nsgf
               rtmat(k1, k2) = 0.5_dp*(rtmat(k1, k2) + rtmat(k2, k1))
               rtmat(k2, k1) = rtmat(k1, k2)
            END DO
         END DO

         ! clean up
         CALL atom_int_release(integrals)
         CALL atom_ppint_release(integrals)
         CALL atom_relint_release(integrals)
         CALL release_atom_basis(basis)
         CALL release_atom_potential(potential)
         CALL release_atom_type(atom)

         DEALLOCATE (potential, basis, integrals)

      ELSE

         IF (ASSOCIATED(rtmat)) THEN
            DEALLOCATE (rtmat)
         END IF
         NULLIFY (rtmat)

      END IF

   END SUBROUTINE calculate_atomic_relkin

! **************************************************************************************************
!> \brief ...
!> \param gth_potential ...
!> \param gth_atompot ...
! **************************************************************************************************
   SUBROUTINE gth_potential_conversion(gth_potential, gth_atompot)
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(atom_gthpot_type)                             :: gth_atompot

      INTEGER                                            :: i, j, l, lm, n, ne, nexp_lpot, nexp_lsd, &
                                                            nexp_nlcc
      INTEGER, DIMENSION(:), POINTER                     :: nct_lpot, nct_lsd, nct_nlcc, nppnl, &
                                                            ppeconf
      LOGICAL                                            :: lpot_present, lsd_present, nlcc_present
      REAL(KIND=dp)                                      :: ac, zeff
      REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha_lpot, alpha_lsd, alpha_nlcc, ap, ce
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cval_lpot, cval_lsd, cval_nlcc
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: hp

      CALL get_potential(gth_potential, &
                         zeff=zeff, &
                         elec_conf=ppeconf, &
                         alpha_core_charge=ac, &
                         nexp_ppl=ne, &
                         cexp_ppl=ce, &
                         lppnl=lm, &
                         nprj_ppnl=nppnl, &
                         alpha_ppnl=ap, &
                         hprj_ppnl=hp)

      gth_atompot%zion = zeff
      gth_atompot%rc = SQRT(0.5_dp/ac)
      gth_atompot%ncl = ne
      gth_atompot%cl(:) = 0._dp
      IF (ac > 0._dp) THEN
         DO i = 1, ne
            gth_atompot%cl(i) = ce(i)/(2._dp*ac)**(i - 1)
         END DO
      END IF
      !extended type
      gth_atompot%lpotextended = .FALSE.
      gth_atompot%lsdpot = .FALSE.
      gth_atompot%nlcc = .FALSE.
      gth_atompot%nexp_lpot = 0
      gth_atompot%nexp_lsd = 0
      gth_atompot%nexp_nlcc = 0
      CALL get_potential(gth_potential, &
                         lpot_present=lpot_present, &
                         lsd_present=lsd_present, &
                         nlcc_present=nlcc_present)
      IF (lpot_present) THEN
         CALL get_potential(gth_potential, &
                            nexp_lpot=nexp_lpot, &
                            alpha_lpot=alpha_lpot, &
                            nct_lpot=nct_lpot, &
                            cval_lpot=cval_lpot)
         gth_atompot%lpotextended = .TRUE.
         gth_atompot%nexp_lpot = nexp_lpot
         gth_atompot%alpha_lpot(1:nexp_lpot) = SQRT(0.5_dp/alpha_lpot(1:nexp_lpot))
         gth_atompot%nct_lpot(1:nexp_lpot) = nct_lpot(1:nexp_lpot)
         DO j = 1, nexp_lpot
            ac = alpha_lpot(j)
            DO i = 1, 4
               gth_atompot%cval_lpot(i, j) = cval_lpot(i, j)/(2._dp*ac)**(i - 1)
            END DO
         END DO
      END IF
      IF (lsd_present) THEN
         CALL get_potential(gth_potential, &
                            nexp_lsd=nexp_lsd, &
                            alpha_lsd=alpha_lsd, &
                            nct_lsd=nct_lsd, &
                            cval_lsd=cval_lsd)
         gth_atompot%lsdpot = .TRUE.
         gth_atompot%nexp_lsd = nexp_lsd
         gth_atompot%alpha_lsd(1:nexp_lsd) = SQRT(0.5_dp/alpha_lsd(1:nexp_lsd))
         gth_atompot%nct_lsd(1:nexp_lsd) = nct_lsd(1:nexp_lsd)
         DO j = 1, nexp_lpot
            ac = alpha_lsd(j)
            DO i = 1, 4
               gth_atompot%cval_lsd(i, j) = cval_lsd(i, j)/(2._dp*ac)**(i - 1)
            END DO
         END DO
      END IF

      ! nonlocal part
      gth_atompot%nl(:) = 0
      gth_atompot%rcnl(:) = 0._dp
      gth_atompot%hnl(:, :, :) = 0._dp
      DO l = 0, lm
         n = nppnl(l)
         gth_atompot%nl(l) = n
         gth_atompot%rcnl(l) = SQRT(0.5_dp/ap(l))
         gth_atompot%hnl(1:n, 1:n, l) = hp(1:n, 1:n, l)
      END DO

      IF (nlcc_present) THEN
         CALL get_potential(gth_potential, &
                            nexp_nlcc=nexp_nlcc, &
                            alpha_nlcc=alpha_nlcc, &
                            nct_nlcc=nct_nlcc, &
                            cval_nlcc=cval_nlcc)
         gth_atompot%nlcc = .TRUE.
         gth_atompot%nexp_nlcc = nexp_nlcc
         gth_atompot%alpha_nlcc(1:nexp_nlcc) = alpha_nlcc(1:nexp_nlcc)
         gth_atompot%nct_nlcc(1:nexp_nlcc) = nct_nlcc(1:nexp_nlcc)
         gth_atompot%cval_nlcc(1:4, 1:nexp_nlcc) = cval_nlcc(1:4, 1:nexp_nlcc)
      END IF

   END SUBROUTINE gth_potential_conversion

! **************************************************************************************************
!> \brief ...
!> \param sgp_potential ...
!> \param sgp_atompot ...
! **************************************************************************************************
   SUBROUTINE sgp_potential_conversion(sgp_potential, sgp_atompot)
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential
      TYPE(atom_sgppot_type)                             :: sgp_atompot

      INTEGER                                            :: lm, n
      INTEGER, DIMENSION(:), POINTER                     :: ppeconf
      LOGICAL                                            :: nlcc_present
      REAL(KIND=dp)                                      :: ac, zeff
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ap, ce
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hhp
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: ccp

      CALL get_potential(sgp_potential, &
                         name=sgp_atompot%pname, &
                         zeff=zeff, &
                         elec_conf=ppeconf, &
                         alpha_core_charge=ac)
      sgp_atompot%zion = zeff
      sgp_atompot%ac_local = ac
      sgp_atompot%econf(0:3) = ppeconf(0:3)
      CALL get_potential(sgp_potential, lmax=lm, &
                         is_nonlocal=sgp_atompot%is_nonlocal, &
                         n_nonlocal=n, a_nonlocal=ap, h_nonlocal=hhp, c_nonlocal=ccp)
      ! nonlocal
      sgp_atompot%has_nonlocal = ANY(sgp_atompot%is_nonlocal)
      sgp_atompot%lmax = lm
      IF (sgp_atompot%has_nonlocal) THEN
         CPASSERT(n <= SIZE(sgp_atompot%a_nonlocal))
         sgp_atompot%n_nonlocal = n
         sgp_atompot%a_nonlocal(1:n) = ap(1:n)
         sgp_atompot%h_nonlocal(1:n, 0:lm) = hhp(1:n, 0:lm)
         sgp_atompot%c_nonlocal(1:n, 1:n, 0:lm) = ccp(1:n, 1:n, 0:lm)
      END IF
      ! local
      CALL get_potential(sgp_potential, n_local=n, a_local=ap, c_local=ce)
      CPASSERT(n <= SIZE(sgp_atompot%a_local))
      sgp_atompot%n_local = n
      sgp_atompot%a_local(1:n) = ap(1:n)
      sgp_atompot%c_local(1:n) = ce(1:n)
      ! NLCC
      CALL get_potential(sgp_potential, has_nlcc=nlcc_present, &
                         n_nlcc=n, a_nlcc=ap, c_nlcc=ce)
      IF (nlcc_present) THEN
         sgp_atompot%has_nlcc = .TRUE.
         CPASSERT(n <= SIZE(sgp_atompot%a_nlcc))
         sgp_atompot%n_nlcc = n
         sgp_atompot%a_nlcc(1:n) = ap(1:n)
         sgp_atompot%c_nlcc(1:n) = ce(1:n)
      ELSE
         sgp_atompot%has_nlcc = .FALSE.
      END IF

   END SUBROUTINE sgp_potential_conversion
! **************************************************************************************************

END MODULE atom_kind_orbitals
