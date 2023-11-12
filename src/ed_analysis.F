!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate Energy Decomposition analysis
!> \par History
!>      07.2023 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE ed_analysis
   USE admm_types,                      ONLY: admm_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: Eriksen2020,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_diag, cp_fm_get_element, cp_fm_get_info, cp_fm_get_submatrix, &
        cp_fm_release, cp_fm_set_all, cp_fm_set_element, cp_fm_set_submatrix, cp_fm_to_fm, &
        cp_fm_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_distribution_type, dbcsr_dot, dbcsr_get_info, &
        dbcsr_p_type, dbcsr_release, dbcsr_reserve_diag_blocks, dbcsr_scale, dbcsr_set, &
        dbcsr_type, dbcsr_type_symmetric
   USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals
   USE hartree_local_types,             ONLY: ecoul_1center_type
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
   USE iao_analysis,                    ONLY: iao_calculate_dmat,&
                                              iao_charges,&
                                              iao_wfn_analysis,&
                                              print_center_spread
   USE iao_types,                       ONLY: iao_env_type,&
                                              iao_set_default
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE,&
                                              pw_type
   USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                              calculate_ecore_self
   USE qs_dispersion_pairpot,           ONLY: calculate_dispersion_pairpot
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_gcp_method,                   ONLY: calculate_gcp_pairpot
   USE qs_gcp_types,                    ONLY: qs_gcp_type
   USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                              integrate_v_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_atom,                      ONLY: update_ks_atom
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_local_rho_types,              ONLY: local_rho_type
   USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              duplicate_mo_set,&
                                              get_mo_set,&
                                              mo_set_type
   USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace
   USE qs_rho_atom_types,               ONLY: rho_atom_type,&
                                              zero_rho_atom_integrals
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_xc_density
   USE qs_vxc_atom,                     ONLY: calculate_vxc_atom
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ed_analysis'

   PUBLIC ::  edmf_analysis

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param input_section ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE edmf_analysis(qs_env, input_section, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input_section
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=*), PARAMETER                        :: routineN = 'edmf_analysis'

      INTEGER                                            :: handle, iatom, ikind, iorb, ispin, jorb, &
                                                            nao, natom, nimages, nkind, no, norb, &
                                                            nref, nspin
      INTEGER, DIMENSION(2)                              :: nocc
      INTEGER, DIMENSION(:), POINTER                     :: refbas_blk_sizes
      LOGICAL                                            :: detailed_ener, do_hfx, explicit, gapw, &
                                                            gapw_xc, ref_orb_canonical, &
                                                            skip_localize, uniform_occupation, uocc
      REAL(KIND=dp)                                      :: ateps, checksum, e1, e2, e_pot, ecc, &
                                                            egcp, ehfx, ekts, evdw, focc, &
                                                            sum_energy
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ate1h, ate1xc, atecc, ateks, atener, &
                                                            odiag
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atdet, mcharge, mweight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: bcenter
      REAL(KIND=dp), DIMENSION(:), POINTER               :: occupation_numbers
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: cvec, cvec2, smo
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: c_iao_coef, ciao, fij_mat, orb_weight, &
                                                            rotmat
      TYPE(cp_fm_type), POINTER                          :: cloc, moref
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:)      :: exc_mat, ks_mat, vhxc_mat
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_h, matrix_hfx, matrix_ks, &
                                                            matrix_s, matrix_t
      TYPE(dbcsr_type)                                   :: dkmat, dmat
      TYPE(dbcsr_type), POINTER                          :: smat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: ref_basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: refbasis
      TYPE(iao_env_type)                                 :: iao_env
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos, mos_loc
      TYPE(mp_comm_type)                                 :: group
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_gcp_type), POINTER                         :: gcp_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input

      CALL section_vals_get(input_section, explicit=explicit)
      IF (.NOT. explicit) RETURN

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="(/,4(T2,A))") &
            "!-----------------------------------------------------------------------------!", &
            "!                        ENERGY DECOMPOSITION ANALYSIS                        !", &
            "!                    Janus J Eriksen, JCP 153 214109 (2020)                   !", &
            "!-----------------------------------------------------------------------------!"
      END IF
      CALL cite_reference(Eriksen2020)

      ! input options
      CALL section_vals_val_get(input_section, "REFERENCE_ORB_CANONICAL", l_val=ref_orb_canonical)
      CALL section_vals_val_get(input_section, "DETAILED_ENERGY", l_val=detailed_ener)
      CALL section_vals_val_get(input_section, "SKIP_LOCALIZATION", l_val=skip_localize)
      ! k-points?
      CALL get_qs_env(qs_env, dft_control=dft_control)
      nimages = dft_control%nimages
      IF (nimages > 1) THEN
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T2,A)") &
               "K-Points: MAO's determined and analyzed using Gamma-Point only."
         END IF
      END IF
      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc
      CALL get_qs_env(qs_env, input=input)
      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(hfx_sections, explicit=do_hfx)

      CALL get_qs_env(qs_env, mos=mos)
      nspin = dft_control%nspins
      ALLOCATE (mos_loc(nspin))

      ! do we have a uniform occupation
      uocc = .TRUE.
      DO ispin = 1, nspin
         CALL get_mo_set(mos(ispin), uniform_occupation=uniform_occupation)
         IF (.NOT. uniform_occupation) uocc = .FALSE.
      END DO
      IF (unit_nr > 0) THEN
         IF (uocc) THEN
            WRITE (UNIT=unit_nr, FMT="(T2,A)") &
               "MO's have a uniform occupation pattern"
         ELSE
            WRITE (UNIT=unit_nr, FMT="(T2,A)") &
               "MO's have varying occupations"
         END IF
      END IF

      ! perform IAO analysis
      CALL iao_set_default(iao_env)
      iao_env%do_iao = .TRUE.
      iao_env%do_charges = .TRUE.
      IF (skip_localize) THEN
         iao_env%do_bondorbitals = .FALSE.
         iao_env%do_center = .FALSE.
      ELSE
         iao_env%do_bondorbitals = .TRUE.
         iao_env%do_center = .TRUE.
      END IF
      iao_env%eps_occ = 1.0E-4_dp
      CALL get_qs_env(qs_env, cell=cell)
      iao_env%pos_periodic = .NOT. ALL(cell%perd == 0)
      no = 0
      nocc = 0
      DO ispin = 1, nspin
         CALL duplicate_mo_set(mos_loc(ispin), mos(ispin))
         CALL get_mo_set(mos_loc(ispin), nmo=norb)
         no = MAX(no, norb)
         nocc(ispin) = norb
         IF (ref_orb_canonical) THEN
            CALL get_mo_set(mos_loc(ispin), mo_coeff=moref)
            CALL get_qs_env(qs_env, matrix_ks=matrix_ks)
            CALL calculate_subspace_eigenvalues(moref, matrix_ks(ispin)%matrix, &
                                                do_rotation=.TRUE.)
         END IF
      END DO
      ALLOCATE (bcenter(5, no, nspin))
      bcenter = 0.0_dp
      CALL iao_wfn_analysis(qs_env, iao_env, unit_nr, &
                            c_iao_coef=c_iao_coef, mos=mos_loc, bond_centers=bcenter)

      ! output bond centers
      IF (iao_env%do_bondorbitals .AND. iao_env%do_center) THEN
         CALL print_center_spread(bcenter, nocc, iounit=unit_nr)
      END IF

      ! Calculate orbital rotation matrix
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      smat => matrix_s(1)%matrix
      ALLOCATE (rotmat(nspin))
      DO ispin = 1, nspin
         CALL get_mo_set(mos_loc(ispin), mo_coeff=cloc)
         CALL get_mo_set(mos(ispin), mo_coeff=moref)
         CALL cp_fm_get_info(cloc, nrow_global=nao, ncol_global=norb)
         CALL cp_fm_create(smo, cloc%matrix_struct)
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, nrow_global=norb, ncol_global=norb, &
                                  template_fmstruct=cloc%matrix_struct)
         CALL cp_fm_create(rotmat(ispin), fm_struct)
         CALL cp_fm_struct_release(fm_struct)
         ! ROTMAT = Cref(T)*S*Cloc
         CALL cp_dbcsr_sm_fm_multiply(smat, cloc, smo, ncol=norb)
         CALL parallel_gemm('T', 'N', norb, norb, nao, 1.0_dp, moref, &
                            smo, 0.0_dp, rotmat(ispin))
         CALL cp_fm_release(smo)
      END DO

      ! calculate occupation matrix
      IF (.NOT. uocc) THEN
         ALLOCATE (fij_mat(nspin))
         DO ispin = 1, nspin
            CALL cp_fm_create(fij_mat(ispin), rotmat(ispin)%matrix_struct)
            CALL cp_fm_set_all(fij_mat(ispin), 0.0_dp)
            CALL cp_fm_create(smo, rotmat(ispin)%matrix_struct)
            ! fii = f
            CALL get_mo_set(mos(ispin), nmo=norb, occupation_numbers=occupation_numbers)
            DO iorb = 1, norb
               CALL cp_fm_set_element(fij_mat(ispin), iorb, iorb, occupation_numbers(iorb))
            END DO
            ! fij = U(T)*f*U
            CALL parallel_gemm('N', 'N', norb, norb, norb, 1.0_dp, fij_mat(ispin), &
                               rotmat(ispin), 0.0_dp, smo)
            CALL parallel_gemm('T', 'N', norb, norb, norb, 1.0_dp, rotmat(ispin), &
                               smo, 0.0_dp, fij_mat(ispin))
            CALL cp_fm_release(smo)
         END DO
      END IF

      ! localized orbitals in IAO basis => CIAO
      ALLOCATE (ciao(nspin))
      DO ispin = 1, nspin
         CALL get_mo_set(mos_loc(ispin), mo_coeff=cloc)
         CALL cp_fm_get_info(cloc, nrow_global=nao, ncol_global=norb)
         CALL cp_fm_get_info(c_iao_coef(ispin), ncol_global=nref)
         CALL cp_fm_create(smo, cloc%matrix_struct)
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nref, &
                                  template_fmstruct=cloc%matrix_struct)
         CALL cp_fm_create(ciao(ispin), fm_struct)
         CALL cp_fm_struct_release(fm_struct)
         ! CIAO = A(T)*S*C
         CALL cp_dbcsr_sm_fm_multiply(smat, cloc, smo, ncol=norb)
         CALL parallel_gemm('T', 'N', nref, norb, nao, 1.0_dp, c_iao_coef(ispin), &
                            smo, 0.0_dp, ciao(ispin))
         CALL cp_fm_release(smo)
      END DO

      ! Reference basis set
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      nkind = SIZE(qs_kind_set)
      ALLOCATE (ref_basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         NULLIFY (ref_basis_set_list(ikind)%gto_basis_set)
         NULLIFY (refbasis)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=refbasis, basis_type="MIN")
         IF (ASSOCIATED(refbasis)) ref_basis_set_list(ikind)%gto_basis_set => refbasis
      END DO
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, natom=natom)
      ALLOCATE (refbas_blk_sizes(natom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=refbas_blk_sizes, &
                            basis=ref_basis_set_list)

      ! Atomic orbital weights
      ALLOCATE (orb_weight(nspin))
      ALLOCATE (mcharge(natom, 1))
      DO ispin = 1, nspin
         CALL get_mo_set(mos_loc(ispin), mo_coeff=cloc)
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, nrow_global=natom, &
                                  template_fmstruct=cloc%matrix_struct)
         CALL cp_fm_create(orb_weight(ispin), fm_struct)
         CALL cp_fm_struct_release(fm_struct)
         CALL cp_fm_set_all(orb_weight(ispin), 0.0_dp)
      END DO
      !
      CALL dbcsr_get_info(smat, distribution=dbcsr_dist)
      CALL dbcsr_create(matrix=dmat, name="DMAT", &
                        dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric, &
                        row_blk_size=refbas_blk_sizes, col_blk_size=refbas_blk_sizes, &
                        nze=0)
      CALL dbcsr_reserve_diag_blocks(dmat)
      !
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, ncol_global=1, &
                               template_fmstruct=ciao(1)%matrix_struct)
      CALL cp_fm_create(cvec, fm_struct)
      CALL cp_fm_create(cvec2, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      !
      DO ispin = 1, nspin
         CALL get_mo_set(mos_loc(ispin), &
                         mo_coeff=cloc, nmo=norb, &
                         occupation_numbers=occupation_numbers)
         IF (uocc) THEN
            ! uniform occupation
            DO iorb = 1, norb
               CALL cp_fm_to_fm(ciao(ispin), cvec, ncol=1, source_start=iorb, target_start=1)
               focc = occupation_numbers(iorb)
               CALL dbcsr_set(dmat, 0.0_dp)
               CALL cp_dbcsr_plus_fm_fm_t(dmat, cvec, ncol=1, alpha=focc)
               CALL iao_charges(dmat, mcharge(:, 1))
               CALL cp_fm_set_submatrix(orb_weight(ispin), mcharge, start_row=1, &
                                        start_col=iorb, n_rows=natom, n_cols=1)
               checksum = SUM(mcharge(:, 1))
               IF (ABS(focc - checksum) > 1.E-6_dp) THEN
                  CALL cp_warn(__LOCATION__, "Sum of atomic orbital weights is incorrect")
                  IF (unit_nr > 0) THEN
                     WRITE (UNIT=unit_nr, FMT="(T2,A,F10.6,T40,A,F10.6)") &
                        "Orbital occupation:", focc, &
                        "Sum of atomic orbital weights:", checksum
                  END IF
               END IF
            END DO
         ELSE
            ! non-diagonal occupation matrix
            ALLOCATE (odiag(norb))
            CALL cp_fm_get_diag(fij_mat(ispin), odiag)
            DO iorb = 1, norb
               IF (odiag(iorb) < 1.E-8_dp) CYCLE
               CALL dbcsr_set(dmat, 0.0_dp)
               CALL cp_fm_to_fm(ciao(ispin), cvec, ncol=1, source_start=iorb, target_start=1)
               DO jorb = 1, norb
                  CALL cp_fm_get_element(fij_mat(ispin), iorb, jorb, focc)
                  IF (focc < 1.E-12_dp) CYCLE
                  CALL cp_fm_to_fm(ciao(ispin), cvec2, ncol=1, source_start=jorb, target_start=1)
                  CALL cp_dbcsr_plus_fm_fm_t(dmat, cvec, cvec2, 1, alpha=focc, symmetry_mode=1)
               END DO
               CALL iao_charges(dmat, mcharge(:, 1))
               checksum = SUM(mcharge(:, 1))
               focc = odiag(iorb)
               IF (ABS(focc - checksum) > 1.E-6_dp) THEN
                  CALL cp_warn(__LOCATION__, "Sum of atomic orbital weights is incorrect")
                  IF (unit_nr > 0) THEN
                     WRITE (UNIT=unit_nr, FMT="(T2,A,F10.6,T40,A,F10.6)") &
                        "Orbital occupation:", focc, &
                        "Sum of atomic orbital weights:", checksum
                  END IF
               END IF
               mcharge(:, 1) = mcharge(:, 1)/focc
               CALL cp_fm_set_submatrix(orb_weight(ispin), mcharge, start_row=1, &
                                        start_col=iorb, n_rows=natom, n_cols=1)
            END DO
            DEALLOCATE (odiag)
         END IF
      END DO
      DEALLOCATE (mcharge)
      CALL dbcsr_release(dmat)
      CALL cp_fm_release(cvec)
      CALL cp_fm_release(cvec2)

      ! energy arrays
      ALLOCATE (atener(natom), ateks(natom), atecc(natom))
      ALLOCATE (ate1xc(natom), ate1h(natom))
      atener = 0.0_dp
      ateks = 0.0_dp
      atecc = 0.0_dp
      ate1xc = 0.0_dp
      ate1h = 0.0_dp
      IF (detailed_ener) THEN
         ALLOCATE (atdet(natom, 10))
         atdet = 0.0_dp
      END IF
      ! atom dependent density matrix
      CALL dbcsr_create(dkmat, template=smat)
      CALL dbcsr_copy(dkmat, smat)
      CALL dbcsr_set(dkmat, 0.0_dp)
      ! KS matrix + correction
      CALL get_qs_env(qs_env, matrix_h=matrix_h, matrix_ks=matrix_ks, kinetic=matrix_t)
      ALLOCATE (ks_mat(nspin), vhxc_mat(nspin), exc_mat(1))
      DO ispin = 1, nspin
         ALLOCATE (ks_mat(ispin)%matrix)
         CALL dbcsr_create(ks_mat(ispin)%matrix, template=matrix_h(1)%matrix)
         CALL dbcsr_copy(ks_mat(ispin)%matrix, matrix_h(1)%matrix)
         CALL dbcsr_add(ks_mat(ispin)%matrix, matrix_ks(ispin)%matrix, 1.0_dp, 1.0_dp)
         CALL dbcsr_scale(ks_mat(ispin)%matrix, 0.5_dp)
         !
         ALLOCATE (vhxc_mat(ispin)%matrix)
         CALL dbcsr_create(vhxc_mat(ispin)%matrix, template=smat)
         CALL dbcsr_copy(vhxc_mat(ispin)%matrix, smat)
         CALL dbcsr_set(vhxc_mat(ispin)%matrix, 0.0_dp)
      END DO
      ALLOCATE (exc_mat(1)%matrix)
      CALL dbcsr_create(exc_mat(1)%matrix, template=smat)
      CALL dbcsr_copy(exc_mat(1)%matrix, smat)
      CALL dbcsr_set(exc_mat(1)%matrix, 0.0_dp)
      !
      CALL vhxc_correction(qs_env, vhxc_mat, exc_mat, atecc, ate1xc, ate1h)
      !
      DO ispin = 1, nspin
         CALL dbcsr_add(ks_mat(ispin)%matrix, vhxc_mat(ispin)%matrix, 1.0_dp, 1.0_dp)
      END DO
      !
      IF (detailed_ener .AND. do_hfx) THEN
         ALLOCATE (matrix_hfx(nspin))
         DO ispin = 1, nspin
            ALLOCATE (matrix_hfx(nspin)%matrix)
            CALL dbcsr_create(matrix_hfx(ispin)%matrix, template=smat)
            CALL dbcsr_copy(matrix_hfx(ispin)%matrix, smat)
            CALL dbcsr_set(matrix_hfx(ispin)%matrix, 0.0_dp)
         END DO
         CALL get_hfx_ks_matrix(qs_env, matrix_hfx)
      END IF
      ! Loop over spins and atoms
      DO ispin = 1, nspin
         CALL get_mo_set(mos_loc(ispin), mo_coeff=cloc, nmo=norb)
         ALLOCATE (mweight(1, norb))
         DO iatom = 1, natom
            CALL cp_fm_get_submatrix(orb_weight(ispin), mweight, start_row=iatom, &
                                     start_col=1, n_rows=1, n_cols=norb)
            IF (uocc) THEN
               CALL iao_calculate_dmat(cloc, dkmat, mweight(1, :), .FALSE.)
            ELSE
               CALL iao_calculate_dmat(cloc, dkmat, mweight(1, :), fij_mat(ispin))
            END IF
            CALL dbcsr_dot(dkmat, ks_mat(ispin)%matrix, ecc)
            ateks(iatom) = ateks(iatom) + ecc
            !
            IF (detailed_ener) THEN
               CALL dbcsr_dot(dkmat, matrix_t(1)%matrix, e1)
               atdet(iatom, 1) = atdet(iatom, 1) + e1
               CALL dbcsr_dot(dkmat, matrix_h(1)%matrix, e2)
               atdet(iatom, 2) = atdet(iatom, 2) + (e2 - e1)
               IF (do_hfx) THEN
                  CALL dbcsr_scale(matrix_hfx(ispin)%matrix, 0.5_dp)
                  CALL dbcsr_dot(dkmat, matrix_hfx(ispin)%matrix, ehfx)
                  atdet(iatom, 5) = atdet(iatom, 5) + ehfx
               ELSE
                  ehfx = 0.0_dp
               END IF
               CALL dbcsr_dot(dkmat, exc_mat(1)%matrix, e1)
               atdet(iatom, 3) = atdet(iatom, 3) + ecc - e2 - e1 - ehfx
               atdet(iatom, 4) = atdet(iatom, 4) + e1
            END IF
            !
         END DO
         DEALLOCATE (mweight)
      END DO
      !
      CALL dbcsr_release(dkmat)
      DO ispin = 1, nspin
         CALL dbcsr_release(ks_mat(ispin)%matrix)
         CALL dbcsr_release(vhxc_mat(ispin)%matrix)
         DEALLOCATE (ks_mat(ispin)%matrix, vhxc_mat(ispin)%matrix)
         CALL deallocate_mo_set(mos_loc(ispin))
      END DO
      CALL dbcsr_release(exc_mat(1)%matrix)
      DEALLOCATE (exc_mat(1)%matrix)
      DEALLOCATE (ks_mat, vhxc_mat, exc_mat)
      DEALLOCATE (mos_loc)
      DEALLOCATE (refbas_blk_sizes)
      DEALLOCATE (ref_basis_set_list)
      IF (detailed_ener .AND. do_hfx) THEN
         DO ispin = 1, nspin
            CALL dbcsr_release(matrix_hfx(ispin)%matrix)
            DEALLOCATE (matrix_hfx(ispin)%matrix)
         END DO
         DEALLOCATE (matrix_hfx)
      END IF
      CALL cp_fm_release(orb_weight)
      CALL cp_fm_release(ciao)
      CALL cp_fm_release(rotmat)
      CALL cp_fm_release(c_iao_coef)
      IF (.NOT. uocc) THEN
         CALL cp_fm_release(fij_mat)
      END IF

      CALL get_qs_env(qs_env, para_env=para_env)
      group = para_env
      ! KS energy
      atener(1:natom) = ateks(1:natom)
      ! core energy corrections
      CALL group%sum(atecc)
      atener(1:natom) = atener(1:natom) + atecc(1:natom)
      IF (detailed_ener) atdet(1:natom, 6) = atdet(1:natom, 6) + atecc(1:natom)
      ! one center terms (GAPW)
      CALL group%sum(ate1xc)
      CALL group%sum(ate1h)
      atener(1:natom) = atener(1:natom) + ate1xc(1:natom) + ate1h(1:natom)
      IF (detailed_ener) THEN
         IF (gapw .OR. gapw_xc) atdet(1:natom, 8) = atdet(1:natom, 8) + ate1xc(1:natom)
         IF (gapw) atdet(1:natom, 9) = atdet(1:natom, 9) + ate1h(1:natom)
      END IF
      ! core correction
      atecc(1:natom) = 0.0_dp
      CALL calculate_ecore_overlap(qs_env, para_env, .FALSE., atecc=atecc)
      CALL group%sum(atecc)
      atener(1:natom) = atener(1:natom) + atecc(1:natom)
      IF (detailed_ener) atdet(1:natom, 7) = atdet(1:natom, 7) + atecc(1:natom)
      atecc(1:natom) = 0.0_dp
      CALL calculate_ecore_self(qs_env, atecc=atecc)
      CALL group%sum(atecc)
      atener(1:natom) = atener(1:natom) + atecc(1:natom)
      IF (detailed_ener) atdet(1:natom, 7) = atdet(1:natom, 7) + atecc(1:natom)
      ! vdW pair-potentials
      atecc(1:natom) = 0.0_dp
      CALL get_qs_env(qs_env=qs_env, dispersion_env=dispersion_env)
      CALL calculate_dispersion_pairpot(qs_env, dispersion_env, evdw, .FALSE., atevdw=atecc)
      ! Pair potential gCP energy
      CALL get_qs_env(qs_env=qs_env, gcp_env=gcp_env)
      IF (ASSOCIATED(gcp_env)) THEN
         CALL calculate_gcp_pairpot(qs_env, gcp_env, egcp, .FALSE., ategcp=atecc)
      END IF
      CALL group%sum(atecc)
      atener(1:natom) = atener(1:natom) + atecc(1:natom)
      IF (detailed_ener) atdet(1:natom, 10) = atdet(1:natom, 10) + atecc(1:natom)
      ! distribute the entropic energy
      CALL get_qs_env(qs_env, energy=energy)
      ekts = energy%kts/REAL(natom, KIND=dp)
      atener(1:natom) = atener(1:natom) + ekts

      IF (detailed_ener) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, FMT="(/,T2,A)") "Detailed Atomic Energy Components "
            CALL write_atdet(unit_nr, atdet(:, 1), "Kinetic Energy")
            CALL write_atdet(unit_nr, atdet(:, 2), "Short-Range Core and/or PP Energy")
            CALL write_atdet(unit_nr, atdet(:, 3), "Hartree Energy (long-ranged part)")
            CALL write_atdet(unit_nr, atdet(:, 5), "Exact Exchange Energy")
            CALL write_atdet(unit_nr, atdet(:, 4), "Exchange-Correlation Energy")
            CALL write_atdet(unit_nr, atdet(:, 6), "Atomic Core Hartree: Int(nc V(n+nc) dx")
            CALL write_atdet(unit_nr, atdet(:, 7), "Core Self Energy Corrections")
            IF (gapw) THEN
               CALL write_atdet(unit_nr, atdet(:, 9), "GAPW: 1-center Hartree Energy")
            END IF
            IF (gapw .OR. gapw_xc) THEN
               CALL write_atdet(unit_nr, atdet(:, 8), "GAPW: 1-center XC Energy")
            END IF
            IF (ABS(evdw) > 1.E-14_dp) THEN
               CALL write_atdet(unit_nr, atdet(:, 10), "vdW Pairpotential Energy")
            END IF
            DO iatom = 1, natom
               atecc(iatom) = SUM(atdet(iatom, 1:10)) - atener(iatom)
            END DO
            CALL write_atdet(unit_nr, atecc(:), "Missing Energy")
         END IF
      END IF

      IF (unit_nr > 0) THEN
         e_pot = energy%total
         ateps = 1.E-6_dp
         CALL write_atener(unit_nr, particle_set, atener, "Atomic Energy Decomposition")
         sum_energy = SUM(atener(:))
         checksum = ABS(e_pot - sum_energy)
         WRITE (UNIT=unit_nr, FMT="(/,(T2,A,T56,F25.13))") &
            "Potential energy (Atomic):", sum_energy, &
            "Potential energy (Total) :", e_pot, &
            "Difference               :", checksum
         CPASSERT((checksum < ateps*ABS(e_pot)))
      END IF

      IF (detailed_ener) THEN
         DEALLOCATE (atdet)
      END IF

      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="(/,T2,A)") &
            "!--------------------------- END OF ED ANALYSIS ------------------------------!"
      END IF
      DEALLOCATE (bcenter)
      DEALLOCATE (atener, ateks, atecc, ate1xc, ate1h)

      CALL timestop(handle)

   END SUBROUTINE edmf_analysis

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param vhxc_mat ...
!> \param exc_mat ...
!> \param atecc ...
!> \param ate1xc ...
!> \param ate1h ...
! **************************************************************************************************
   SUBROUTINE vhxc_correction(qs_env, vhxc_mat, exc_mat, atecc, ate1xc, ate1h)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: vhxc_mat, exc_mat
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: atecc, ate1xc, ate1h

      CHARACTER(len=*), PARAMETER                        :: routineN = 'vhxc_correction'

      INTEGER                                            :: handle, iatom, ispin, natom, nspins
      LOGICAL                                            :: do_admm_corr, gapw, gapw_xc
      REAL(KIND=dp)                                      :: eh1, exc1
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ecoul_1center_type), DIMENSION(:), POINTER    :: ecoul_1c
      TYPE(local_rho_type), POINTER                      :: local_rho_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: xc_den
      TYPE(pw_type), ALLOCATABLE, DIMENSION(:)           :: vtau, vxc
      TYPE(pw_type), POINTER                             :: v_hartree_rspace
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      TYPE(section_vals_type), POINTER                   :: xc_fun_section, xc_section
      TYPE(xc_rho_cflags_type)                           :: needs

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, ks_env=ks_env, dft_control=dft_control, pw_env=pw_env)

      nspins = dft_control%nspins
      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc
      do_admm_corr = .FALSE.
      IF (dft_control%do_admm) THEN
         IF (qs_env%admm_env%aux_exch_func /= do_admm_aux_exch_func_none) do_admm_corr = .TRUE.
      END IF
      IF (do_admm_corr) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         xc_section => admm_env%xc_section_aux
      ELSE
         xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      END IF
      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      needs = xc_functionals_get_needs(xc_fun_section, (nspins == 2), .TRUE.)

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL auxbas_pw_pool%create_pw(xc_den, &
                                    use_data=REALDATA3D, in_space=REALSPACE)
      ALLOCATE (vxc(nspins))
      DO ispin = 1, nspins
         CALL auxbas_pw_pool%create_pw(vxc(ispin), &
                                       use_data=REALDATA3D, in_space=REALSPACE)
      END DO
      IF (needs%tau .OR. needs%tau_spin) THEN
         ALLOCATE (vtau(nspins))
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%create_pw(vtau(ispin), &
                                          use_data=REALDATA3D, in_space=REALSPACE)
         END DO
      END IF

      ! Nuclear charge correction
      CALL get_qs_env(qs_env, v_hartree_rspace=v_hartree_rspace)
      CALL integrate_v_core_rspace(v_hartree_rspace, qs_env, atecc=atecc)
      IF (gapw .OR. gapw_xc) THEN
         CALL get_qs_env(qs_env=qs_env, local_rho_set=local_rho_set, &
                         rho_atom_set=rho_atom_set, ecoul_1c=ecoul_1c, &
                         natom=natom, para_env=para_env)
         CALL zero_rho_atom_integrals(rho_atom_set)
         CALL calculate_vxc_atom(qs_env, .FALSE., exc1)
         IF (gapw) THEN
            CALL Vh_1c_gg_integrals(qs_env, eh1, ecoul_1c, local_rho_set, para_env, tddft=.FALSE.)
            CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)
            CALL integrate_vhg0_rspace(qs_env, v_hartree_rspace, para_env, calculate_forces=.FALSE., &
                                       local_rho_set=local_rho_set, atener=atecc)
         END IF
      END IF

      IF (gapw_xc) THEN
         CALL get_qs_env(qs_env, rho_xc=rho_struct, dispersion_env=dispersion_env)
      ELSE
         CALL get_qs_env(qs_env, rho=rho_struct, dispersion_env=dispersion_env)
      END IF
      !
      CPASSERT(.NOT. do_admm_corr)
      !
      IF (needs%tau .OR. needs%tau_spin) THEN
         CALL qs_xc_density(ks_env, rho_struct, xc_section, dispersion_env=dispersion_env, &
                            xc_den=xc_den, vxc=vxc, vtau=vtau)
      ELSE
         CALL qs_xc_density(ks_env, rho_struct, xc_section, dispersion_env=dispersion_env, &
                            xc_den=xc_den, vxc=vxc)
      END IF
      DO ispin = 1, nspins
         CALL pw_scale(vxc(ispin), -0.5_dp)
         CALL pw_axpy(xc_den, vxc(ispin))
         CALL pw_scale(vxc(ispin), vxc(ispin)%pw_grid%dvol)
         CALL integrate_v_rspace(qs_env=qs_env, v_rspace=vxc(ispin), &
                                 hmat=vhxc_mat(ispin), calculate_forces=.FALSE., &
                                 gapw=(gapw .OR. gapw_xc))
         IF (needs%tau .OR. needs%tau_spin) THEN
            CALL pw_scale(vtau(ispin), -0.5_dp*vtau(ispin)%pw_grid%dvol)
            CALL integrate_v_rspace(qs_env=qs_env, v_rspace=vtau(ispin), &
                                    hmat=vhxc_mat(ispin), calculate_forces=.FALSE., &
                                    compute_tau=.TRUE., gapw=(gapw .OR. gapw_xc))
         END IF
      END DO
      CALL pw_scale(xc_den, xc_den%pw_grid%dvol)
      CALL integrate_v_rspace(qs_env=qs_env, v_rspace=xc_den, &
                              hmat=exc_mat(1), calculate_forces=.FALSE., &
                              gapw=(gapw .OR. gapw_xc))

      IF (gapw .OR. gapw_xc) THEN
         ! remove one-center potential matrix part
         CALL qs_rho_get(rho_struct, rho_ao_kp=matrix_p)
         CALL update_ks_atom(qs_env, vhxc_mat, matrix_p, forces=.FALSE., kscale=-0.5_dp)
         !
         DO iatom = 1, natom
            ate1xc(iatom) = ate1xc(iatom) + &
                            rho_atom_set(iatom)%exc_h - rho_atom_set(iatom)%exc_s
         END DO
         IF (gapw) THEN
            DO iatom = 1, natom
               ate1h(iatom) = ate1h(iatom) + &
                              ecoul_1c(iatom)%ecoul_1_h - ecoul_1c(iatom)%ecoul_1_s + &
                              ecoul_1c(iatom)%ecoul_1_z - ecoul_1c(iatom)%ecoul_1_0
            END DO
         END IF
      END IF

      CALL auxbas_pw_pool%give_back_pw(xc_den)
      DO ispin = 1, nspins
         CALL auxbas_pw_pool%give_back_pw(vxc(ispin))
      END DO
      IF (needs%tau .OR. needs%tau_spin) THEN
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%give_back_pw(vtau(ispin))
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE vhxc_correction

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix_hfx ...
! **************************************************************************************************
   SUBROUTINE get_hfx_ks_matrix(qs_env, matrix_hfx)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hfx

      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_hfx_ks_matrix'

      INTEGER                                            :: handle
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_ao=matrix_p)
      CALL tddft_hfx_matrix(matrix_hfx, matrix_p, qs_env, update_energy=.FALSE., &
                            recalc_integrals=.FALSE.)

      CALL timestop(handle)

   END SUBROUTINE get_hfx_ks_matrix
! **************************************************************************************************
!> \brief Write the atomic coordinates, types, and energies
!> \param iounit ...
!> \param particle_set ...
!> \param atener ...
!> \param label ...
!> \date    05.06.2023
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE write_atener(iounit, particle_set, atener, label)

      INTEGER, INTENT(IN)                                :: iounit
      TYPE(particle_type), DIMENSION(:)                  :: particle_set
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: atener
      CHARACTER(LEN=*), INTENT(IN)                       :: label

      INTEGER                                            :: i, natom

      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") TRIM(label)
         WRITE (UNIT=iounit, FMT="(T4,A,T30,A,T42,A,T54,A,T69,A)") &
            "Atom  Element", "X", "Y", "Z", "Energy[a.u.]"
         natom = SIZE(atener)
         DO i = 1, natom
            WRITE (UNIT=iounit, FMT="(I6,T12,A2,T24,3F12.6,F21.10)") i, &
               TRIM(ADJUSTL(particle_set(i)%atomic_kind%element_symbol)), &
               particle_set(i)%r(1:3)*angstrom, atener(i)
         END DO
         WRITE (UNIT=iounit, FMT="(A)") ""
      END IF

   END SUBROUTINE write_atener

! **************************************************************************************************
!> \brief Write the one component of atomic energies
!> \param iounit ...
!> \param atener ...
!> \param label ...
!> \date    18.08.2023
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE write_atdet(iounit, atener, label)
      INTEGER, INTENT(IN)                                :: iounit
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: atener
      CHARACTER(LEN=*), INTENT(IN)                       :: label

      INTEGER                                            :: natom

      IF (iounit > 0) THEN
         natom = SIZE(atener)
         WRITE (UNIT=iounit, FMT="(/,T2,A)") TRIM(label)
         WRITE (UNIT=iounit, FMT="(5G16.8)") atener(1:natom)
      END IF

   END SUBROUTINE write_atdet

END MODULE ed_analysis
