!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Analytical band velocities and effective mass tensors for periodic LCAO
!>        calculations from k-derivatives of the converged real-space Hamiltonian and
!>        overlap matrices, evaluated with zero all-rank real-space matrix replication via
!>        distributed rskp_transform, cp_cfm_type eigensolves, and streaming parallel_gemm.
! **************************************************************************************************
MODULE qs_effective_mass
   USE admm_types,                      ONLY: admm_type
   USE bibliography,                    ONLY: LuttingerKohn1955
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale,&
                                              cp_cfm_scale_and_add_fm,&
                                              cp_cfm_transpose
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig,&
                                              cp_cfm_geeig_canon,&
                                              cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_no_symmetry,&
                                              dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_fm_copy_plan_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE effective_mass_methods,          ONLY: &
        carrier_tol, degenerate_range, displace_fractional, kato_rellich_reduce, principal_masses, &
        refine_bz_bounds, spectral_newton_step, sylvester_sqrt_tile, wrap_fractional
   USE input_constants,                 ONLY: cholesky_off,&
                                              do_admm_purify_none,&
                                              do_method_dftb,&
                                              do_method_gapw,&
                                              do_method_gapw_xc,&
                                              do_method_gpw,&
                                              do_method_xtb,&
                                              plus_u_tensorial
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_transform,                ONLY: rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE mathconstants,                   ONLY: gaussi,&
                                              sqrt3,&
                                              sqrthalf,&
                                              twopi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: a_fine,&
                                              angstrom,&
                                              c_light,&
                                              evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: get_neighbor_list_set_p
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE reference_manager,               ONLY: cite_reference
   USE scf_control_types,               ONLY: scf_control_type
   USE string_utilities,                ONLY: uppercase
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_effective_mass'

   !> All 13 unique cubic directions (axes, face/body diagonals) covering +-26 lattice directions.
   REAL(dp), PARAMETER, PRIVATE :: invr3 = 1.0_dp/sqrt3
   REAL(dp), PARAMETER, PRIVATE :: refine_probe_dirs(3, 13) = RESHAPE([ &
                                                                      1.0_dp, 0.0_dp, 0.0_dp, &
                                                                      0.0_dp, 1.0_dp, 0.0_dp, &
                                                                      0.0_dp, 0.0_dp, 1.0_dp, &
                                                                      sqrthalf, sqrthalf, 0.0_dp, &
                                                                      sqrthalf, -sqrthalf, 0.0_dp, &
                                                                      sqrthalf, 0.0_dp, sqrthalf, &
                                                                      sqrthalf, 0.0_dp, -sqrthalf, &
                                                                      0.0_dp, sqrthalf, sqrthalf, &
                                                                      0.0_dp, sqrthalf, -sqrthalf, &
                                                                      invr3, invr3, invr3, &
                                                                      invr3, invr3, -invr3, &
                                                                      invr3, -invr3, invr3, &
                                                                      -invr3, invr3, invr3], [3, 13])

   !> Indices into refine_probe_dirs for the 7 Luttinger-Kohn analysis directions.
   INTEGER, PARAMETER, PRIVATE :: lk_dir_idx(7) = [1, 2, 3, 4, 6, 8, 10]

   !> Refinement parameters and state machine constants.
   REAL(dp), PARAMETER, PRIVATE :: refine_grad_tol = 1.0E-7_dp, &
                                   refine_polish_tol = 1.0E-10_dp, &
                                   refine_escape_curv = -1.0E-6_dp, &
                                   refine_decrease_tol = 1.0E-12_dp, &
                                   refine_alpha_floor = 1.0_dp/64.0_dp
   INTEGER, PARAMETER, PRIVATE  :: refine_max_iter = 25, &
                                   refine_stepping = 1, &
                                   refine_polishing = 2, &
                                   refine_converged = 3, &
                                   refine_stalled = 4

   !> Absolute fractional tolerance below which two seed k-points are the same point.
   REAL(dp), PARAMETER, PRIVATE :: seed_dedup_tol = 1.0E-9_dp

   !> Minimum overlap eigenvalue for full-rank S(k)^(1/2) when Cholesky decomposition is enabled.
   REAL(dp), PARAMETER, PRIVATE :: lowdin_spd_tol = 1.0E-12_dp

   !> Fixed carrier-type labels of the principal_masses classification constants.
   CHARACTER(LEN=13), PARAMETER, PRIVATE :: carrier_labels(0:3) = &
                                            [CHARACTER(LEN=13) :: 'flat', 'electron-like', 'hole-like', 'saddle point']

   !> Persistent world-distributed workspace and shared handles for one effective-mass call.
   TYPE, PRIVATE :: emass_work_type
      INTEGER                                            :: iw = 0
      INTEGER                                            :: nao = 0
      INTEGER                                            :: nlowdin = 0
      INTEGER                                            :: nspin = 0
      LOGICAL                                            :: has_lowdin = .FALSE.
      REAL(dp)                                           :: eps_degeneracy = 0.0_dp
      REAL(dp), ALLOCATABLE                              :: evals(:)
      REAL(dp), ALLOCATABLE                              :: lowdin_sigma(:)
      TYPE(cell_type), POINTER                           :: cell => NULL()
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env => NULL()
      TYPE(cp_cfm_type)                                  :: cfm_frame
      TYPE(cp_cfm_type)                                  :: cfm_h
      TYPE(cp_cfm_type)                                  :: cfm_s
      TYPE(cp_cfm_type)                                  :: cfm_u
      TYPE(cp_cfm_type)                                  :: cfm_work
      TYPE(cp_cfm_type)                                  :: lowdin_up
      TYPE(cp_cfm_type)                                  :: lowdin_w0
      TYPE(cp_fm_struct_type), POINTER                   :: struct_nm => NULL()
      TYPE(cp_fm_struct_type), POINTER                   :: struct_nn => NULL()
      TYPE(cp_fm_type)                                   :: fm_work
      TYPE(dbcsr_fm_copy_plan_type)                      :: plan_im
      TYPE(dbcsr_fm_copy_plan_type)                      :: plan_re
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s => NULL()
      TYPE(dbcsr_type)                                   :: cmatrix
      TYPE(dbcsr_type)                                   :: rmatrix
      TYPE(kpoint_type), POINTER                         :: kpoints => NULL()
      TYPE(mp_para_env_type), POINTER                    :: para_env => NULL()
      TYPE(qs_scf_env_type), POINTER                     :: scf_env => NULL()
      TYPE(scf_control_type), POINTER                    :: scf_control => NULL()
   END TYPE emass_work_type

   PUBLIC :: calculate_effective_mass, check_effective_mass_preconditions

CONTAINS

! **************************************************************************************************
!> \brief Fail fast on operators whose post-SCF k-space Hamiltonian is incomplete or unsupported
!>        when EFFECTIVE_MASS is active.
!> \param qs_env Quickstep environment
! **************************************************************************************************
   SUBROUTINE check_effective_mass_preconditions(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      LOGICAL                                            :: do_hfx, do_kpoints, explicit
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_im
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: em_input, hfx_section, input_section

      NULLIFY (admm_env, dft_control, matrix_ks_im, em_input, hfx_section, input_section)
      CALL get_qs_env(qs_env, input=input_section, do_kpoints=do_kpoints, dft_control=dft_control, &
                      admm_env=admm_env, matrix_ks_im_kp=matrix_ks_im)
      em_input => section_vals_get_subs_vals(input_section, "DFT%PRINT%BAND_STRUCTURE%EFFECTIVE_MASS")
      CALL section_vals_get(em_input, explicit=explicit)
      IF (.NOT. explicit) RETURN

      IF (dft_control%roks .OR. dft_control%low_spin_roks) THEN
         CPABORT('EFFECTIVE_MASS does not support ROKS.')
      END IF
      IF (ASSOCIATED(matrix_ks_im)) THEN
         CPABORT('EFFECTIVE_MASS does not support complex KS matrices.')
      END IF
      IF (dft_control%dft_plus_u .AND. dft_control%plus_u_method_id == plus_u_tensorial) THEN
         CPABORT('EFFECTIVE_MASS does not support PLUS_U_METHOD TENSORIAL.')
      END IF
      IF (dft_control%do_admm) THEN
         IF (.NOT. dft_control%do_admm_dm .AND. &
             admm_env%purification_method /= do_admm_purify_none) THEN
            CPABORT('EFFECTIVE_MASS requires ADMM purification NONE/NONE_DM.')
         END IF
      END IF

      IF (.NOT. do_kpoints) THEN
         SELECT CASE (dft_control%qs_control%method_id)
         CASE (do_method_gpw, do_method_gapw, do_method_gapw_xc)
         CASE DEFAULT
            CPABORT('Gamma-point EFFECTIVE_MASS supports GPW/GAPW/GAPW_XC only.')
         END SELECT
         do_hfx = .FALSE.
         hfx_section => section_vals_get_subs_vals(input_section, "DFT%XC%HF")
         CALL section_vals_get(hfx_section, explicit=do_hfx)
         IF (do_hfx) THEN
            CPABORT('Gamma-point EFFECTIVE_MASS rejects HFX; use KPOINTS.')
         END IF
         IF (dft_control%do_admm) THEN
            CPABORT('Gamma-point EFFECTIVE_MASS rejects ADMM; use KPOINTS.')
         END IF
         IF (dft_control%dft_plus_u) THEN
            CPABORT('Gamma-point EFFECTIVE_MASS rejects DFT+U; use KPOINTS.')
         END IF
      ELSE
         SELECT CASE (dft_control%qs_control%method_id)
         CASE (do_method_gpw, do_method_gapw, do_method_gapw_xc, do_method_dftb, do_method_xtb)
         CASE DEFAULT
            CPABORT('EFFECTIVE_MASS supports GPW/GAPW/GAPW_XC/DFTB/xTB only.')
         END SELECT
      END IF
   END SUBROUTINE check_effective_mass_preconditions

! **************************************************************************************************
!> \brief Driver of the &BAND_STRUCTURE&EFFECTIVE_MASS print section: assemble the target k-point
!>        list and evaluate band velocities and effective masses via world-distributed Bloch
!>        transforms, ScaLAPACK eigensolves, and streaming GEMMs.
!> \param qs_env Converged k-point shaped Quickstep environment
!> \param from_gamma True when qs_env was reconstructed from a Gamma-point SCF
!> \note All ranks participate in lockstep in the distributed matrix transforms, eigensolves,
!>       and narrow-block projections; only the four nao x nao cp_cfm_type matrices and one
!>       DBCSR->FM work buffer are persistent, avoiding any all-rank real-space replication.
!>       AUTO_EXTREMA seeds band-edge candidates from the SCF k-points plus the
!>       EXTREMA_SEARCH_MESH (time-reversal partners skipped) by diagonalizing the complete AO
!>       spectrum; REFINE_EXTREMA then polishes the best VBM/CBM seeds off the mesh with
!>       analytical spectral Newton steps until the band velocity vanishes.
! **************************************************************************************************
   SUBROUTINE calculate_effective_mass(qs_env, from_gamma)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: from_gamma

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_effective_mass'

      CHARACTER(LEN=default_string_length)               :: ustr
      INTEGER :: handle, i1, i2, i3, ik, ik_cbm, ik_vbm, ispin, iw, n_occ(2), nao, nc, nelectron, &
         nexplicit, nkp_scf, nseed, nspin, nval, search_mesh(3)
      INTEGER, DIMENSION(:), POINTER                     :: band_indices, mesh_vals
      LOGICAL                                            :: auto_extrema, explicit, refine_extrema
      REAL(dp)                                           :: cart_hmat(3, 3), energy_window, &
                                                            eps_degeneracy, kfrac(3), maxocc
      REAL(dp), ALLOCATABLE                              :: ev_homo(:, :), ev_lumo(:, :), &
                                                            kseeds(:, :), xkp_explicit(:, :)
      REAL(dp), DIMENSION(:), POINTER                    :: rvals
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(emass_work_type)                              :: work
      TYPE(kpoint_type), POINTER                         :: scf
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(section_vals_type), POINTER                   :: em_input, input_section

      CALL timeset(routineN, handle)

      NULLIFY (band_indices, rvals, cell, dft_control, mos, matrix_ks, scf, em_input, input_section)
      CALL get_qs_env(qs_env, dft_control=dft_control, cell=cell, input=input_section, &
                      matrix_ks_kp=matrix_ks, kpoints=scf, mos=mos)
      em_input => section_vals_get_subs_vals(input_section, "DFT%PRINT%BAND_STRUCTURE%EFFECTIVE_MASS")
      CALL section_vals_get(em_input, explicit=explicit)
      IF (.NOT. explicit) THEN
         CALL timestop(handle)
         RETURN
      END IF
      CALL cite_reference(LuttingerKohn1955)

      CALL section_vals_val_get(em_input, "AUTO_EXTREMA", l_val=auto_extrema)
      CALL section_vals_val_get(em_input, "REFINE_EXTREMA", l_val=refine_extrema)
      CALL section_vals_val_get(em_input, "EXTREMA_SEARCH_MESH", i_vals=mesh_vals)
      CPASSERT(SIZE(mesh_vals) == 3)
      search_mesh(1:3) = mesh_vals(1:3)
      CALL section_vals_val_get(em_input, "EPS_DEGENERACY", r_val=eps_degeneracy)
      CALL section_vals_val_get(em_input, "ENERGY_WINDOW", r_val=energy_window)
      CALL section_vals_val_get(em_input, "BAND_INDICES", n_rep_val=nval)
      IF (nval > 0) CALL section_vals_val_get(em_input, "BAND_INDICES", i_vals=band_indices)
      CALL section_vals_val_get(em_input, "UNITS", c_val=ustr)
      CALL uppercase(ustr)

      nspin = SIZE(matrix_ks, 1)
      CALL get_mo_set(mo_set=mos(1), nao=nao, nelectron=nelectron, maxocc=maxocc)
      n_occ(1) = NINT(REAL(nelectron, dp)/maxocc)
      IF (nspin == 2) THEN
         CALL get_mo_set(mo_set=mos(2), nelectron=nelectron, maxocc=maxocc)
         n_occ(2) = NINT(REAL(nelectron, dp)/maxocc)
      END IF

      IF (eps_degeneracy <= 0.0_dp) THEN
         CPABORT('EPS_DEGENERACY must be positive')
      END IF
      IF (auto_extrema) THEN
         IF (ANY(n_occ(1:nspin) < 1) .OR. ANY(n_occ(1:nspin) >= nao)) THEN
            CPABORT('AUTO_EXTREMA requires occupied and empty bands in every spin channel')
         END IF
         IF (ANY(search_mesh < 1)) THEN
            CPABORT('EXTREMA_SEARCH_MESH entries must be positive')
         END IF
      END IF

      ! Fractional coordinates of the explicitly requested k-points, in the cell frame of
      ! the operators (input_hmat when canonicalized): kfrac_a = (kcart . a_a)/(2 pi).
      CALL section_vals_val_get(em_input, "KPOINT", n_rep_val=nexplicit)
      IF (nexplicit > 0) THEN
         IF (ASSOCIATED(band_indices)) THEN
            IF (SIZE(band_indices) < 1) THEN
               CPABORT('BAND_INDICES must not be empty')
            END IF
            IF (ANY(band_indices < 1) .OR. ANY(band_indices > nao)) THEN
               CPABORT('BAND_INDICES entry outside the AO spectrum')
            END IF
         ELSE
            IF (energy_window < 0.0_dp) THEN
               CPABORT('ENERGY_WINDOW must be non-negative')
            END IF
            IF (ANY(n_occ(1:nspin) < 1) .OR. ANY(n_occ(1:nspin) >= nao)) THEN
               CPABORT('Band-edge selection needs occupied and empty bands; use BAND_INDICES')
            END IF
         END IF
      END IF

      iw = cp_logger_get_default_io_unit()
      IF (.NOT. auto_extrema .AND. nexplicit == 0) THEN
         IF (iw > 0) THEN
            WRITE (iw, '(T2,A)') 'EFFECTIVE_MASS| Nothing to evaluate (AUTO_EXTREMA F and no KPOINT).'
         END IF
         CALL timestop(handle)
         RETURN
      END IF

      cart_hmat(:, :) = cell%hmat(:, :)
      IF (cell%input_cell_canonicalized) cart_hmat(:, :) = cell%input_hmat(:, :)
      ALLOCATE (xkp_explicit(3, nexplicit))
      DO ik = 1, nexplicit
         CALL section_vals_val_get(em_input, "KPOINT", i_rep_val=ik, r_vals=rvals)
         CPASSERT(SIZE(rvals) == 3)
         SELECT CASE (ustr)
         CASE ("B_VECTOR")
            xkp_explicit(1:3, ik) = rvals(1:3)
         CASE ("CART_ANGSTROM")
            xkp_explicit(1:3, ik) = MATMUL(rvals(1:3), cart_hmat)*angstrom/twopi
         CASE ("CART_BOHR")
            xkp_explicit(1:3, ik) = MATMUL(rvals(1:3), cart_hmat)/twopi
         CASE DEFAULT
            CPABORT("Unknown unit <"//TRIM(ustr)//"> specified for k-point definition")
         END SELECT
      END DO

      CALL emass_work_create(work, qs_env, nao, nspin, eps_degeneracy, iw)

      IF (iw > 0) THEN
         WRITE (iw, '(T2,A)') ' '
         WRITE (iw, '(T2,A)') REPEAT('-', 79)
         WRITE (iw, '(T2,A)') ' EFFECTIVE MASS CALCULATION'
         WRITE (iw, '(T2,A)') ' Analytical k-derivatives in the full overlap metric,'// &
            ' hbar = m_e = 1 (Luttinger-Kohn k.p)'
         WRITE (iw, '(T2,A,ES12.6,A)') ' Convention: m* > 0 electron, m* < 0 hole | 1 a.u. vel = ', &
            a_fine*c_light, ' m/s'
         WRITE (iw, '(T2,A)') REPEAT('-', 79)
         IF (cell%input_cell_canonicalized) THEN
            WRITE (iw, '(T2,A)') 'EFFECTIVE_MASS| Note: Cartesian k, velocities, and tensors'// &
               ' are in CP2K canonical cell frame.'
         END IF
      END IF

      IF (auto_extrema) THEN
         ! Non-periodic cell directions carry no dispersion; do not sample them.
         WHERE (cell%perd(1:3) == 0) search_mesh(1:3) = 1
         IF (from_gamma) THEN
            nkp_scf = 1
         ELSE
            CALL get_kpoint_info(scf, nkp=nkp_scf)
         END IF
         ! Seed candidates: the SCF k-points plus a uniform search mesh; duplicate and
         ! redundant time-reversal partners are skipped because the spectra of the
         ! real-space operators are inversion symmetric.
         ALLOCATE (kseeds(3, nkp_scf + search_mesh(1)*search_mesh(2)*search_mesh(3)))
         nseed = 0
         DO ik = 1, nkp_scf
            kfrac(1:3) = 0.0_dp
            IF (.NOT. from_gamma) kfrac(1:3) = scf%xkp(1:3, ik)
            CALL add_seed(wrap_fractional(kfrac), nseed, kseeds)
         END DO
         DO i3 = 0, search_mesh(3) - 1
            DO i2 = 0, search_mesh(2) - 1
               DO i1 = 0, search_mesh(1) - 1
                  CALL add_seed(wrap_fractional([REAL(i1, dp)/REAL(search_mesh(1), dp), &
                                                 REAL(i2, dp)/REAL(search_mesh(2), dp), &
                                                 REAL(i3, dp)/REAL(search_mesh(3), dp)]), &
                                nseed, kseeds)
               END DO
            END DO
         END DO
         ALLOCATE (ev_homo(nseed, nspin), ev_lumo(nseed, nspin))
         DO ik = 1, nseed
            DO ispin = 1, nspin
               CALL solve_at(work, kseeds(1:3, ik), ispin, work%evals, nc)
               IF (n_occ(ispin) + 1 > nc) THEN
                  CPABORT('AUTO_EXTREMA target band exceeds retained basis states')
               END IF
               ev_homo(ik, ispin) = work%evals(n_occ(ispin))
               ev_lumo(ik, ispin) = work%evals(n_occ(ispin) + 1)
            END DO
         END DO
         IF (iw > 0) THEN
            IF (dft_control%smear) THEN
               WRITE (iw, '(T2,A)') 'EFFECTIVE_MASS| Warning: smearing is enabled; occupations'// &
                  ' are fractional and band edges are not well defined'
            END IF
         END IF
         DO ispin = 1, nspin
            IF (iw > 0) THEN
               IF (MAXVAL(ev_homo(1:nseed, ispin)) >= MINVAL(ev_lumo(1:nseed, ispin))) THEN
                  WRITE (iw, '(T2,A,I0,A)') 'EFFECTIVE_MASS| Warning, spin ', ispin, &
                     ': occupied and empty bands overlap on the seed mesh;'// &
                     ' band-edge localization is unreliable'
               END IF
            END IF
            ik_vbm = MAXLOC(ev_homo(1:nseed, ispin), 1)
            ik_cbm = MINLOC(ev_lumo(1:nseed, ispin), 1)
            IF (refine_extrema) THEN
               CALL refine_edge(work, kseeds(1:3, ik_vbm), ispin, n_occ(ispin), -1, 'VBM')
               CALL refine_edge(work, kseeds(1:3, ik_cbm), ispin, n_occ(ispin) + 1, 1, 'CBM')
            ELSE
               CALL solve_at(work, kseeds(1:3, ik_vbm), ispin, work%evals, nc)
               CALL report_bands(work, kseeds(1:3, ik_vbm), ispin, nc, [n_occ(ispin)], 'VBM (auto)')
               CALL solve_at(work, kseeds(1:3, ik_cbm), ispin, work%evals, nc)
               CALL report_bands(work, kseeds(1:3, ik_cbm), ispin, nc, [n_occ(ispin) + 1], 'CBM (auto)')
            END IF
         END DO
         DEALLOCATE (ev_homo, ev_lumo, kseeds)
      END IF

      DO ik = 1, nexplicit
         CALL report_explicit(work, xkp_explicit(1:3, ik), n_occ, band_indices, energy_window)
      END DO

      IF (iw > 0) THEN
         WRITE (iw, '(T2,A)') ' '
      END IF

      CALL emass_work_release(work)
      DEALLOCATE (xkp_explicit)

      CALL timestop(handle)
   END SUBROUTINE calculate_effective_mass

! **************************************************************************************************
!> \brief Allocate persistent world-distributed matrices and DBCSR workspace once per call.
!> \param work Effective-mass workspace to allocate and initialize
!> \param qs_env Converged k-point shaped Quickstep environment
!> \param nao Number of spherical Gaussian basis functions
!> \param nspin Number of spin channels (1 or 2)
!> \param eps_degeneracy Degeneracy energy threshold (Ha)
!> \param iw Output unit (> 0 on the ionode)
! **************************************************************************************************
   SUBROUTINE emass_work_create(work, qs_env, nao, nspin, eps_degeneracy, iw)
      TYPE(emass_work_type), INTENT(OUT)                 :: work
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: nao, nspin
      REAL(dp), INTENT(IN)                               :: eps_degeneracy
      INTEGER, INTENT(IN)                                :: iw

      LOGICAL                                            :: do_symmetric

      work%iw = iw
      work%nao = nao
      work%nspin = nspin
      work%eps_degeneracy = eps_degeneracy

      CALL get_qs_env(qs_env, cell=work%cell, para_env=work%para_env, &
                      blacs_env=work%blacs_env, matrix_ks_kp=work%matrix_ks, &
                      matrix_s_kp=work%matrix_s, kpoints=work%kpoints, &
                      scf_env=work%scf_env, scf_control=work%scf_control)

      CALL cp_fm_struct_create(work%struct_nn, para_env=work%para_env, context=work%blacs_env, &
                               nrow_global=nao, ncol_global=nao)
      CALL cp_fm_create(work%fm_work, work%struct_nn)
      CALL cp_cfm_create(work%cfm_h, work%struct_nn)
      CALL cp_cfm_create(work%cfm_s, work%struct_nn)
      CALL cp_cfm_create(work%cfm_frame, work%struct_nn)
      CALL cp_cfm_create(work%cfm_work, work%struct_nn)

      work%has_lowdin = .FALSE.
      work%nlowdin = 0
      IF (ASSOCIATED(work%kpoints)) THEN
         IF (ALLOCATED(work%kpoints%lowdin_v) .AND. ALLOCATED(work%kpoints%lowdin_orbitals)) THEN
            IF (SIZE(work%kpoints%lowdin_orbitals) > 0) THEN
               work%has_lowdin = ANY(work%kpoints%lowdin_v /= 0.0_dp)
            END IF
         END IF
      END IF
      IF (work%has_lowdin) THEN
         work%nlowdin = SIZE(work%kpoints%lowdin_orbitals)
         ALLOCATE (work%lowdin_sigma(nao))
         ! Match struct_nn block sizes (force_block) for ScaLAPACK pzgemm compatibility.
         CALL cp_fm_struct_create(work%struct_nm, para_env=work%para_env, context=work%blacs_env, &
                                  nrow_global=nao, ncol_global=work%nlowdin, &
                                  nrow_block=work%struct_nn%nrow_block, &
                                  ncol_block=work%struct_nn%ncol_block, &
                                  square_blocks=.TRUE., force_block=.TRUE.)
         CALL cp_cfm_create(work%cfm_u, work%struct_nn)
         CALL cp_cfm_create(work%lowdin_up, work%struct_nm)
         CALL cp_cfm_create(work%lowdin_w0, work%struct_nm)
      END IF

      CALL get_neighbor_list_set_p(neighbor_list_sets=work%kpoints%sab_nl, symmetric=do_symmetric)
      CALL dbcsr_create(work%rmatrix, template=work%matrix_s(1, 1)%matrix, &
                        matrix_type=MERGE(dbcsr_type_symmetric, dbcsr_type_no_symmetry, do_symmetric))
      CALL dbcsr_create(work%cmatrix, template=work%matrix_s(1, 1)%matrix, &
                        matrix_type=MERGE(dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, do_symmetric))
      CALL cp_dbcsr_alloc_block_from_nbl(work%rmatrix, work%kpoints%sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(work%cmatrix, work%kpoints%sab_nl)

      ALLOCATE (work%evals(nao))
   END SUBROUTINE emass_work_create

! **************************************************************************************************
!> \brief Release persistent world-distributed matrices, copy plans, and spectrum scratch.
!> \param work Effective-mass workspace to release
! **************************************************************************************************
   SUBROUTINE emass_work_release(work)
      TYPE(emass_work_type), INTENT(INOUT)               :: work

      ! Reset copy plans before releasing the target FM and its structure.
      work%plan_re = dbcsr_fm_copy_plan_type()
      work%plan_im = dbcsr_fm_copy_plan_type()
      CALL cp_fm_release(work%fm_work)
      CALL cp_cfm_release(work%cfm_h)
      CALL cp_cfm_release(work%cfm_s)
      CALL cp_cfm_release(work%cfm_frame)
      CALL cp_cfm_release(work%cfm_work)
      IF (work%has_lowdin) THEN
         CALL cp_cfm_release(work%cfm_u)
         CALL cp_cfm_release(work%lowdin_up)
         CALL cp_cfm_release(work%lowdin_w0)
         CALL cp_fm_struct_release(work%struct_nm)
         IF (ALLOCATED(work%lowdin_sigma)) DEALLOCATE (work%lowdin_sigma)
         work%has_lowdin = .FALSE.
         work%nlowdin = 0
      END IF
      CALL cp_fm_struct_release(work%struct_nn)
      CALL dbcsr_release(work%rmatrix)
      CALL dbcsr_release(work%cmatrix)
      IF (ALLOCATED(work%evals)) DEALLOCATE (work%evals)
      NULLIFY (work%cell, work%blacs_env, work%struct_nm, work%struct_nn, &
               work%matrix_ks, work%matrix_s, work%kpoints, work%para_env, &
               work%scf_env, work%scf_control)
   END SUBROUTINE emass_work_release

! **************************************************************************************************
!> \brief Bloch Fourier transform (or its k-derivative) of a real-space DBCSR matrix into cfm.
!> \param work Shared effective-mass workspace
!> \param rsmat Real-space DBCSR matrix array
!> \param row_idx Spin/row index into rsmat (1 for overlap)
!> \param kfrac Fractional k-point coordinates
!> \param cfm World-distributed complex output matrix
!> \param deriv1 Optional first Cartesian k-derivative direction (1..3)
!> \param deriv2 Optional second Cartesian k-derivative direction (1..3)
! **************************************************************************************************
   SUBROUTINE transform_to_cfm(work, rsmat, row_idx, kfrac, cfm, deriv1, deriv2)
      TYPE(emass_work_type), INTENT(INOUT)               :: work
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat
      INTEGER, INTENT(IN)                                :: row_idx
      REAL(dp), INTENT(IN)                               :: kfrac(3)
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: cfm
      INTEGER, INTENT(IN), OPTIONAL                      :: deriv1, deriv2

      CALL dbcsr_set(work%rmatrix, 0.0_dp)
      CALL dbcsr_set(work%cmatrix, 0.0_dp)
      CALL rskp_transform(rmatrix=work%rmatrix, cmatrix=work%cmatrix, rsmat=rsmat, &
                          ispin=row_idx, xkp=kfrac, cell_to_index=work%kpoints%cell_to_index, &
                          sab_nl=work%kpoints%sab_nl, hmat=work%cell%hmat, &
                          deriv_direction=deriv1, deriv_direction2=deriv2)
      CALL copy_dbcsr_to_fm(work%rmatrix, work%fm_work, plan=work%plan_re)
      CALL cp_fm_to_cfm(msourcer=work%fm_work, mtarget=cfm)
      CALL copy_dbcsr_to_fm(work%cmatrix, work%fm_work, plan=work%plan_im)
      CALL cp_cfm_scale_and_add_fm(z_one, cfm, gaussi, work%fm_work)
   END SUBROUTINE transform_to_cfm

! **************************************************************************************************
!> \brief Build H(k) and S(k) distributed and solve the generalized eigenproblem into work%cfm_frame.
!> \param work Shared effective-mass workspace
!> \param kfrac Fractional k-point coordinates
!> \param ispin Spin channel
!> \param evals_out Generalized eigenvalues (broadcast identically across all ranks)
!> \param nc_out Number of retained linear-dependency-free basis states (<= nao)
! **************************************************************************************************
   SUBROUTINE solve_at(work, kfrac, ispin, evals_out, nc_out)
      TYPE(emass_work_type), INTENT(INOUT)               :: work
      REAL(dp), INTENT(IN)                               :: kfrac(3)
      INTEGER, INTENT(IN)                                :: ispin
      REAL(dp), INTENT(OUT)                              :: evals_out(:)
      INTEGER, INTENT(OUT)                               :: nc_out

      INTEGER                                            :: ir
      REAL(dp)                                           :: s_tol
      TYPE(cp_cfm_type)                                  :: lowdin_b0

      CALL transform_to_cfm(work, work%matrix_ks, ispin, kfrac, work%cfm_h)
      CALL transform_to_cfm(work, work%matrix_s, 1, kfrac, work%cfm_s)
      IF (work%has_lowdin) THEN
         CALL cp_cfm_heevd(matrix=work%cfm_s, eigenvectors=work%cfm_u, eigenvalues=work%lowdin_sigma)
         ! Synchronize source-rank values so 1-ulp ScaLAPACK drift cannot desynchronize ranks.
         CALL work%para_env%bcast(work%lowdin_sigma, work%para_env%source)
         s_tol = MERGE(work%scf_control%eps_eigval, lowdin_spd_tol, &
                       work%scf_env%cholesky_method == cholesky_off)
         IF (ANY(work%lowdin_sigma(1:work%nao) <= s_tol)) THEN
            CPABORT('Lowdin EFFECTIVE_MASS requires positive-definite overlap.')
         END IF
         work%lowdin_sigma(1:work%nao) = SQRT(work%lowdin_sigma(1:work%nao))

         CALL cp_cfm_transpose(work%cfm_u, 'C', work%cfm_work)
         DO ir = 1, work%nlowdin
            CALL cp_cfm_to_cfm(work%cfm_work, work%lowdin_up, ncol=1, &
                               source_start=work%kpoints%lowdin_orbitals(ir), target_start=ir)
         END DO
         CALL cp_cfm_create(lowdin_b0, work%struct_nm)
         CALL cp_cfm_column_scale(work%cfm_work, &
                                  CMPLX(work%lowdin_sigma(1:work%nao), 0.0_dp, KIND=dp), &
                                  source=work%cfm_u)
         CALL parallel_gemm("N", "N", work%nao, work%nlowdin, work%nao, z_one, &
                            work%cfm_work, work%lowdin_up, z_zero, lowdin_b0)
         CALL cp_cfm_column_scale(work%lowdin_w0, &
                                  CMPLX(work%kpoints%lowdin_v(1:work%nlowdin, ispin), 0.0_dp, KIND=dp), &
                                  source=lowdin_b0)
         CALL parallel_gemm("N", "C", work%nao, work%nao, work%nlowdin, z_one, &
                            work%lowdin_w0, lowdin_b0, z_one, work%cfm_h)
         CALL cp_cfm_release(lowdin_b0)

         CALL cp_cfm_column_scale(work%cfm_s, &
                                  CMPLX(1.0_dp/work%lowdin_sigma(1:work%nao), 0.0_dp, KIND=dp), &
                                  source=work%cfm_u)
         CALL parallel_gemm("N", "N", work%nao, work%nao, work%nao, z_one, &
                            work%cfm_h, work%cfm_s, z_zero, work%cfm_work)
         CALL parallel_gemm("C", "N", work%nao, work%nao, work%nao, z_one, &
                            work%cfm_s, work%cfm_work, z_zero, work%cfm_h)
         CALL cp_cfm_heevd(matrix=work%cfm_h, eigenvectors=work%cfm_work, eigenvalues=evals_out)
         CALL parallel_gemm("N", "N", work%nao, work%nao, work%nao, z_one, &
                            work%cfm_s, work%cfm_work, z_zero, work%cfm_frame)
         nc_out = work%nao
      ELSE IF (work%scf_env%cholesky_method == cholesky_off) THEN
         CALL cp_cfm_geeig_canon(work%cfm_h, work%cfm_s, work%cfm_frame, evals_out, work%cfm_work, &
                                 work%scf_control%eps_eigval, nmo_retained=nc_out)
      ELSE
         CALL cp_cfm_geeig(work%cfm_h, work%cfm_s, work%cfm_frame, evals_out, work%cfm_work)
         nc_out = work%nao
      END IF
      ! Synchronize source-rank values so 1-ulp eigensolver drift cannot desynchronize collective k-steps.
      CALL work%para_env%bcast(evals_out, work%para_env%source)
      CALL work%para_env%bcast(nc_out, work%para_env%source)
   END SUBROUTINE solve_at

! **************************************************************************************************
!> \brief Solve Sylvester S^(1/2) Frechet derivative on cfm_tilde, project cfm_b_out = U * F~ * U_p,
!>        and accumulate cfm_b_out * W0^H + W0 * cfm_b_out^H into work%cfm_h.
!> \param work Shared effective-mass workspace
!> \param cfm_tilde Overlap derivative transformed into the overlap eigenbasis (modified in place)
!> \param lowdin_tmp Scratch matrix of shape work%struct_nm
!> \param cfm_b_out Projected Lowdin orbital derivative of shape work%struct_nm
! **************************************************************************************************
   SUBROUTINE lowdin_sylvester_term(work, cfm_tilde, lowdin_tmp, cfm_b_out)
      TYPE(emass_work_type), INTENT(INOUT)               :: work
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: cfm_tilde, lowdin_tmp, cfm_b_out

      INTEGER                                            :: ncol_s, nrow_s
      INTEGER, DIMENSION(:), POINTER                     :: cidx_s, ridx_s

      CALL cp_cfm_get_info(cfm_tilde, nrow_local=nrow_s, ncol_local=ncol_s, &
                           row_indices=ridx_s, col_indices=cidx_s)
      IF (nrow_s > 0 .AND. ncol_s > 0) THEN
         CALL sylvester_sqrt_tile(cfm_tilde%local_data(1:nrow_s, 1:ncol_s), &
                                  work%lowdin_sigma(ridx_s(1:nrow_s)), &
                                  work%lowdin_sigma(cidx_s(1:ncol_s)))
      END IF
      CALL parallel_gemm("N", "N", work%nao, work%nlowdin, work%nao, z_one, &
                         cfm_tilde, work%lowdin_up, z_zero, lowdin_tmp)
      CALL parallel_gemm("N", "N", work%nao, work%nlowdin, work%nao, z_one, &
                         work%cfm_u, lowdin_tmp, z_zero, cfm_b_out)
      CALL parallel_gemm("N", "C", work%nao, work%nao, work%nlowdin, z_one, &
                         cfm_b_out, work%lowdin_w0, z_one, work%cfm_h)
      CALL parallel_gemm("N", "C", work%nao, work%nao, work%nlowdin, z_one, &
                         work%lowdin_w0, cfm_b_out, z_one, work%cfm_h)
   END SUBROUTINE lowdin_sylvester_term

! **************************************************************************************************
!> \brief Extract group velocity and symmetric 3x3 inverse mass tensor for one non-degenerate band.
!> \param a_sub First-derivative diagonal blocks
!> \param b_sub Effective second-order k.p blocks H^(ab) - E S^(ab) + SOS^(ab) + SOS^(ba)
!> \param s_sub Overlap first-derivative diagonal blocks
!> \param isub Subspace index inside the batched arrays
!> \param vel Group velocity (a.u.)
!> \param inv Symmetric inverse effective-mass tensor (1/m_e)
! **************************************************************************************************
   PURE SUBROUTINE extract_nondeg_tensor(a_sub, b_sub, s_sub, isub, vel, inv)
      COMPLEX(dp), INTENT(IN)                            :: a_sub(:, :, :, :), b_sub(:, :, :, :, :), &
                                                            s_sub(:, :, :, :)
      INTEGER, INTENT(IN)                                :: isub
      REAL(dp), INTENT(OUT)                              :: vel(3), inv(3, 3)

      INTEGER                                            :: a, b
      REAL(dp)                                           :: s1(3)

      DO a = 1, 3
         vel(a) = REAL(a_sub(1, 1, a, isub), dp)
         s1(a) = REAL(s_sub(1, 1, a, isub), dp)
      END DO
      DO a = 1, 3
         DO b = a, 3
            inv(a, b) = REAL(b_sub(1, 1, a, b, isub), dp) - vel(a)*s1(b) - vel(b)*s1(a)
            inv(b, a) = inv(a, b)
         END DO
      END DO
   END SUBROUTINE extract_nondeg_tensor

! **************************************************************************************************
!> \brief Unified distributed projection of first-order and effective second-order k.p blocks
!>        for one or more disjoint subspaces [sub_lo(k) .. sub_hi(k)] at the solved frame.
!> \param work Shared effective-mass workspace containing cfm_frame and evals
!> \param kfrac Fractional k-point coordinates matching work%cfm_frame
!> \param ispin Spin channel
!> \param nc Number of retained states in work%cfm_frame
!> \param sub_lo Lower band indices of disjoint subspaces
!> \param sub_hi Upper band indices of disjoint subspaces
!> \param a_sub Projected velocity operators [H^(a) - e0 S^(a)]_DD (d_max, d_max, 3, nsub)
!> \param b_sub Symmetric second-order k.p blocks (d_max, d_max, 3, 3, nsub)
!> \param s_sub Projected overlap first derivatives S^(a)_DD (d_max, d_max, 3, nsub)
! **************************************************************************************************
   SUBROUTINE analyze_subspaces(work, kfrac, ispin, nc, sub_lo, sub_hi, a_sub, b_sub, s_sub)
      TYPE(emass_work_type), INTENT(INOUT)               :: work
      REAL(dp), INTENT(IN)                               :: kfrac(3)
      INTEGER, INTENT(IN)                                :: ispin, nc, sub_lo(:), sub_hi(:)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: a_sub(:, :, :, :), &
                                                            b_sub(:, :, :, :, :), s_sub(:, :, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: g_loc(:, :), sos_blk(:, :)
      INTEGER                                            :: a, b, col_hi, col_lo, d, d_max, ihi, &
                                                            ilo, iloc, irow_sub, isub, m, ncol_l, &
                                                            ncol_sub, nrow_l, nsub
      INTEGER, ALLOCATABLE                               :: offsets(:)
      INTEGER, DIMENSION(:), POINTER                     :: cidx, ridx
      REAL(dp)                                           :: e0
      REAL(dp), ALLOCATABLE                              :: e0_sub(:)
      TYPE(cp_cfm_type)                                  :: c_blk, cfm_ftilde(3), h1_proj(3), &
                                                            lowdin_ba(3), lowdin_tmp1, &
                                                            lowdin_tmp2, s1_proj(3), w_blk
      TYPE(cp_fm_struct_type), POINTER                   :: struct_sub

      nsub = SIZE(sub_lo)
      CPASSERT(nsub >= 1 .AND. SIZE(sub_hi) == nsub)
      ALLOCATE (offsets(nsub), e0_sub(nsub))
      ncol_sub = 0
      d_max = 0
      DO isub = 1, nsub
         d = sub_hi(isub) - sub_lo(isub) + 1
         CPASSERT(d >= 1)
         offsets(isub) = ncol_sub
         e0_sub(isub) = SUM(work%evals(sub_lo(isub):sub_hi(isub)))/REAL(d, dp)
         ncol_sub = ncol_sub + d
         d_max = MAX(d_max, d)
      END DO

      NULLIFY (struct_sub)
      ! Keep row_block identical to struct_nn for pzgemm compatibility; place all columns on pcol 0.
      CALL cp_fm_struct_create(struct_sub, para_env=work%para_env, context=work%blacs_env, &
                               nrow_global=work%nao, ncol_global=ncol_sub, &
                               nrow_block=work%struct_nn%nrow_block, ncol_block=ncol_sub, &
                               square_blocks=.FALSE., force_block=.TRUE.)
      CALL cp_cfm_create(c_blk, struct_sub)
      CALL cp_cfm_create(w_blk, struct_sub)
      DO a = 1, 3
         CALL cp_cfm_create(h1_proj(a), struct_sub, set_zero=.TRUE.)
         CALL cp_cfm_create(s1_proj(a), struct_sub, set_zero=.TRUE.)
      END DO
      IF (work%has_lowdin) THEN
         DO a = 1, 3
            CALL cp_cfm_create(cfm_ftilde(a), work%struct_nn)
            CALL cp_cfm_create(lowdin_ba(a), work%struct_nm)
         END DO
         CALL cp_cfm_create(lowdin_tmp1, work%struct_nm)
         CALL cp_cfm_create(lowdin_tmp2, work%struct_nm)
      END IF

      DO isub = 1, nsub
         d = sub_hi(isub) - sub_lo(isub) + 1
         CALL cp_cfm_to_cfm(work%cfm_frame, c_blk, ncol=d, &
                            source_start=sub_lo(isub), target_start=offsets(isub) + 1)
      END DO

      DO a = 1, 3
         CALL transform_to_cfm(work, work%matrix_s, 1, kfrac, work%cfm_s, deriv1=a)
         CALL parallel_gemm("N", "N", work%nao, ncol_sub, work%nao, z_one, work%cfm_s, c_blk, &
                            z_zero, w_blk)
         CALL parallel_gemm("C", "N", nc, ncol_sub, work%nao, z_one, work%cfm_frame, w_blk, &
                            z_zero, s1_proj(a))

         CALL transform_to_cfm(work, work%matrix_ks, ispin, kfrac, work%cfm_h, deriv1=a)
         IF (work%has_lowdin) THEN
            CALL parallel_gemm("N", "N", work%nao, work%nao, work%nao, z_one, &
                               work%cfm_s, work%cfm_u, z_zero, work%cfm_work)
            CALL parallel_gemm("C", "N", work%nao, work%nao, work%nao, z_one, &
                               work%cfm_u, work%cfm_work, z_zero, cfm_ftilde(a))
            CALL lowdin_sylvester_term(work, cfm_ftilde(a), lowdin_tmp1, lowdin_ba(a))
         END IF
         CALL parallel_gemm("N", "N", work%nao, ncol_sub, work%nao, z_one, work%cfm_h, c_blk, &
                            z_zero, w_blk)
         CALL parallel_gemm("C", "N", nc, ncol_sub, work%nao, z_one, work%cfm_frame, w_blk, &
                            z_zero, h1_proj(a))
      END DO

      ALLOCATE (a_sub(d_max, d_max, 3, nsub), b_sub(d_max, d_max, 3, 3, nsub), &
                s_sub(d_max, d_max, 3, nsub))
      a_sub(:, :, :, :) = z_zero
      b_sub(:, :, :, :, :) = z_zero
      s_sub(:, :, :, :) = z_zero

      CALL cp_cfm_get_info(c_blk, nrow_local=nrow_l, ncol_local=ncol_l, &
                           row_indices=ridx, col_indices=cidx)
      IF (ncol_l > 0) THEN
         CPASSERT(ncol_l == ncol_sub .AND. cidx(1) == 1 .AND. cidx(ncol_sub) == ncol_sub)
      END IF
      IF (nrow_l > 0 .AND. ncol_l == ncol_sub) THEN
         DO isub = 1, nsub
            ilo = sub_lo(isub)
            ihi = sub_hi(isub)
            d = ihi - ilo + 1
            col_lo = offsets(isub) + 1
            col_hi = offsets(isub) + d
            e0 = e0_sub(isub)
            DO a = 1, 3
               DO iloc = 1, nrow_l
                  m = ridx(iloc)
                  IF (m < 1 .OR. m > nc) CYCLE
                  h1_proj(a)%local_data(iloc, col_lo:col_hi) = &
                     h1_proj(a)%local_data(iloc, col_lo:col_hi) - &
                     e0*s1_proj(a)%local_data(iloc, col_lo:col_hi)
                  IF (m >= ilo .AND. m <= ihi) THEN
                     irow_sub = m - ilo + 1
                     a_sub(irow_sub, 1:d, a, isub) = a_sub(irow_sub, 1:d, a, isub) + &
                                                     h1_proj(a)%local_data(iloc, col_lo:col_hi)
                     s_sub(irow_sub, 1:d, a, isub) = s_sub(irow_sub, 1:d, a, isub) + &
                                                     s1_proj(a)%local_data(iloc, col_lo:col_hi)
                  END IF
               END DO
            END DO

            ALLOCATE (g_loc(nrow_l, d), sos_blk(d, d))
            DO b = 1, 3
               g_loc(1:nrow_l, 1:d) = z_zero
               DO iloc = 1, nrow_l
                  m = ridx(iloc)
                  IF (m < 1 .OR. m > nc .OR. (m >= ilo .AND. m <= ihi)) CYCLE
                  g_loc(iloc, 1:d) = h1_proj(b)%local_data(iloc, col_lo:col_hi)/(e0 - work%evals(m))
               END DO
               DO a = 1, b
                  sos_blk(1:d, 1:d) = MATMUL(CONJG(TRANSPOSE( &
                                                   h1_proj(a)%local_data(1:nrow_l, col_lo:col_hi))), &
                                             g_loc(1:nrow_l, 1:d))
                  b_sub(1:d, 1:d, a, b, isub) = b_sub(1:d, 1:d, a, b, isub) + &
                                                sos_blk(1:d, 1:d) + CONJG(TRANSPOSE(sos_blk(1:d, 1:d)))
               END DO
            END DO
            DEALLOCATE (g_loc, sos_blk)
         END DO
      END IF

      DO a = 1, 3
         DO b = a, 3
            CALL transform_to_cfm(work, work%matrix_s, 1, kfrac, work%cfm_s, deriv1=a, deriv2=b)
            CALL parallel_gemm("N", "N", work%nao, ncol_sub, work%nao, z_one, work%cfm_s, c_blk, &
                               z_zero, w_blk)
            CALL accumulate_second_deriv(a, b, is_overlap=.TRUE.)

            CALL transform_to_cfm(work, work%matrix_ks, ispin, kfrac, work%cfm_h, deriv1=a, deriv2=b)
            IF (work%has_lowdin) THEN
               CALL parallel_gemm("N", "N", work%nao, work%nao, work%nao, z_one, &
                                  work%cfm_s, work%cfm_u, z_zero, work%cfm_work)
               CALL parallel_gemm("C", "N", work%nao, work%nao, work%nao, z_one, &
                                  work%cfm_u, work%cfm_work, z_zero, work%cfm_s)
               IF (a == b) THEN
                  CALL parallel_gemm("N", "N", work%nao, work%nao, work%nao, -2.0_dp*z_one, &
                                     cfm_ftilde(a), cfm_ftilde(a), z_one, work%cfm_s)
               ELSE
                  CALL parallel_gemm("N", "N", work%nao, work%nao, work%nao, -z_one, &
                                     cfm_ftilde(a), cfm_ftilde(b), z_one, work%cfm_s)
                  CALL parallel_gemm("N", "N", work%nao, work%nao, work%nao, -z_one, &
                                     cfm_ftilde(b), cfm_ftilde(a), z_one, work%cfm_s)
               END IF
               CALL lowdin_sylvester_term(work, work%cfm_s, lowdin_tmp1, lowdin_tmp2)
               CALL cp_cfm_column_scale(lowdin_tmp1, &
                                        CMPLX(work%kpoints%lowdin_v(1:work%nlowdin, ispin), 0.0_dp, KIND=dp), &
                                        source=lowdin_ba(a))
               CALL parallel_gemm("N", "C", work%nao, work%nao, work%nlowdin, z_one, &
                                  lowdin_tmp1, lowdin_ba(b), z_one, work%cfm_h)
               CALL parallel_gemm("N", "C", work%nao, work%nao, work%nlowdin, z_one, &
                                  lowdin_ba(b), lowdin_tmp1, z_one, work%cfm_h)
            END IF
            CALL parallel_gemm("N", "N", work%nao, ncol_sub, work%nao, z_one, work%cfm_h, c_blk, &
                               z_zero, w_blk)
            CALL accumulate_second_deriv(a, b, is_overlap=.FALSE.)
         END DO
      END DO

      DO a = 1, 3
         DO b = a + 1, 3
            b_sub(1:d_max, 1:d_max, b, a, 1:nsub) = b_sub(1:d_max, 1:d_max, a, b, 1:nsub)
         END DO
      END DO

      CALL work%para_env%sum(a_sub)
      CALL work%para_env%sum(s_sub)
      DO isub = 1, nsub
         CALL work%para_env%sum(b_sub(:, :, :, :, isub))
      END DO

      DEALLOCATE (offsets, e0_sub)
      IF (work%has_lowdin) THEN
         DO a = 1, 3
            CALL cp_cfm_release(cfm_ftilde(a))
            CALL cp_cfm_release(lowdin_ba(a))
         END DO
         CALL cp_cfm_release(lowdin_tmp1)
         CALL cp_cfm_release(lowdin_tmp2)
      END IF
      CALL cp_cfm_release(c_blk)
      CALL cp_cfm_release(w_blk)
      DO a = 1, 3
         CALL cp_cfm_release(h1_proj(a))
         CALL cp_cfm_release(s1_proj(a))
      END DO
      CALL cp_fm_struct_release(struct_sub)

   CONTAINS

! **************************************************************************************************
!> \brief Accumulate local C_blk^H W_blk for S^(a_dir, b_dir) or H^(a_dir, b_dir) into b_sub.
!> \param a_dir First Cartesian derivative direction (1..3)
!> \param b_dir Second Cartesian derivative direction (1..3)
!> \param is_overlap True when accumulating -E_0 S^(a,b), false for +H^(a,b)
! **************************************************************************************************
      SUBROUTINE accumulate_second_deriv(a_dir, b_dir, is_overlap)
      INTEGER, INTENT(IN)                                :: a_dir, b_dir
      LOGICAL, INTENT(IN)                                :: is_overlap

      INTEGER                                            :: c_hi, c_lo, d_s, ksub
      REAL(dp)                                           :: wt

         IF (nrow_l <= 0 .OR. ncol_l /= ncol_sub) RETURN
         DO ksub = 1, nsub
            d_s = sub_hi(ksub) - sub_lo(ksub) + 1
            c_lo = offsets(ksub) + 1
            c_hi = offsets(ksub) + d_s
            wt = MERGE(-e0_sub(ksub), 1.0_dp, is_overlap)
            b_sub(1:d_s, 1:d_s, a_dir, b_dir, ksub) = &
               b_sub(1:d_s, 1:d_s, a_dir, b_dir, ksub) + &
               wt*MATMUL(CONJG(TRANSPOSE(c_blk%local_data(1:nrow_l, c_lo:c_hi))), &
                         w_blk%local_data(1:nrow_l, c_lo:c_hi))
         END DO
      END SUBROUTINE accumulate_second_deriv
   END SUBROUTINE analyze_subspaces

! **************************************************************************************************
!> \brief Append a seed k-point unless it or its time-reversal partner is already kept.
!> \param kf Wrapped fractional seed coordinates
!> \param nseed Running count of unique kept seeds
!> \param kseeds Seed k-point array (3, max_seeds)
! **************************************************************************************************
   SUBROUTINE add_seed(kf, nseed, kseeds)
      REAL(dp), INTENT(IN)                               :: kf(3)
      INTEGER, INTENT(INOUT)                             :: nseed
      REAL(dp), INTENT(INOUT)                            :: kseeds(:, :)

      INTEGER                                            :: s

      DO s = 1, nseed
         IF (ALL(ABS(wrap_fractional(kseeds(1:3, s) - kf)) <= seed_dedup_tol)) RETURN
         IF (ALL(ABS(wrap_fractional(kseeds(1:3, s) + kf)) <= seed_dedup_tol)) RETURN
      END DO
      nseed = nseed + 1
      kseeds(1:3, nseed) = kf(1:3)
   END SUBROUTINE add_seed

! **************************************************************************************************
!> \brief Refine one auto band edge off the seed mesh in lockstep across all ranks.
!> \param work Shared effective-mass workspace
!> \param kseed Fractional seed coordinates of the best seed point
!> \param ispin Spin channel
!> \param iband Band index of the edge (n_occ for the VBM, n_occ + 1 for the CBM)
!> \param sense +1 minimizes the band energy (CBM), -1 maximizes it (VBM)
!> \param label 'VBM' or 'CBM' for the printed status line and header
! **************************************************************************************************
   SUBROUTINE refine_edge(work, kseed, ispin, iband, sense, label)
      TYPE(emass_work_type), INTENT(INOUT)               :: work
      REAL(dp), INTENT(IN)                               :: kseed(3)
      INTEGER, INTENT(IN)                                :: ispin, iband, sense
      CHARACTER(LEN=*), INTENT(IN)                       :: label

      INTEGER                                            :: carrier, i, iterations, n_escape, &
                                                            nc_loc, state
      LOGICAL                                            :: frame_valid, is_deg, perd(3)
      REAL(dp) :: alpha, curv_dirs(3, 3), curv_eigs(3), escape_dirs(3, 3), f_cur, f_probe, &
                  f_trial, grad(3), hess(3, 3), inv_cur(3, 3), kf(3), kf_probe(3), kf_trial(3), p(3), &
                  probe_step, sensef, step_cap, vel_cur(3)
      REAL(dp), ALLOCATABLE                              :: trial_evals(:)

      ALLOCATE (trial_evals(work%nao))
      sensef = REAL(sense, dp)
      perd = (work%cell%perd(1:3) /= 0)
      CALL refine_bz_bounds(work%cell%hmat, perd, step_cap, probe_step)
      kf = wrap_fractional(kseed)
      iterations = 0
      frame_valid = .FALSE.

      CALL eval_refine_point(kf, nc_loc, is_deg, vel_cur, inv_cur)
      f_cur = sensef*work%evals(iband)

      IF (.NOT. ANY(perd)) THEN
         state = refine_converged
      ELSE
         state = refine_stepping
      END IF

      DO WHILE (state == refine_stepping)
         IF (.NOT. is_deg .AND. MAXVAL(ABS(vel_cur)) >= refine_grad_tol) THEN
            IF (iterations >= refine_max_iter) THEN
               state = refine_stalled
               EXIT
            END IF
            grad = sensef*vel_cur
            hess = sensef*inv_cur
            CALL spectral_newton_step(grad, hess, p, step_cap=step_cap)
            ! Synchronize step vector so 1-ulp eigensolver drift cannot split rank trajectories.
            CALL work%para_env%bcast(p, work%para_env%source)
            alpha = 1.0_dp
            DO WHILE (alpha >= refine_alpha_floor)
               CALL eval_trial_energy(kf, alpha*p, kf_trial, f_trial)
               IF (f_trial < f_cur) EXIT
               alpha = 0.5_dp*alpha
            END DO
            IF (alpha < refine_alpha_floor) THEN
               state = refine_stalled
               EXIT
            END IF
            kf = kf_trial
            iterations = iterations + 1
            CALL eval_refine_point(kf, nc_loc, is_deg, vel_cur, inv_cur)
            f_cur = sensef*work%evals(iband)
            CYCLE
         END IF

         ! Stationary or degenerate: probe for a descending direction.
         f_probe = f_cur
         kf_probe = kf
         IF (.NOT. is_deg) THEN
            hess = sensef*inv_cur
            CALL principal_masses(hess, curv_eigs, curv_dirs, carrier)
            CALL work%para_env%bcast(curv_eigs, work%para_env%source)
            CALL work%para_env%bcast(curv_dirs, work%para_env%source)
            n_escape = 0
            DO i = 1, 3
               IF (curv_eigs(i) < refine_escape_curv) THEN
                  n_escape = n_escape + 1
                  escape_dirs(1:3, n_escape) = curv_dirs(1:3, i)
               END IF
            END DO
            IF (n_escape > 0) THEN
               CALL probe_directions(kf, escape_dirs(1:3, 1:n_escape), kf_probe, f_probe)
            END IF
         ELSE
            CALL probe_directions(kf, refine_probe_dirs, kf_probe, f_probe)
         END IF
         IF (f_probe >= f_cur - refine_decrease_tol*MAX(1.0_dp, ABS(f_cur))) THEN
            IF (.NOT. is_deg .AND. MAXVAL(ABS(vel_cur)) > refine_polish_tol) THEN
               state = refine_polishing
               EXIT
            END IF
            state = refine_converged
            EXIT
         END IF
         IF (iterations >= refine_max_iter) THEN
            state = refine_stalled
            EXIT
         END IF
         kf = kf_probe
         iterations = iterations + 1
         CALL eval_refine_point(kf, nc_loc, is_deg, vel_cur, inv_cur)
         f_cur = sensef*work%evals(iband)
      END DO

      IF (state == refine_polishing) THEN
         grad = sensef*vel_cur
         hess = sensef*inv_cur
         CALL spectral_newton_step(grad, hess, p, step_cap=step_cap)
         CALL work%para_env%bcast(p, work%para_env%source)
         CALL eval_trial_energy(kf, p, kf_trial, f_trial)
         IF (f_trial <= f_cur) THEN
            kf = kf_trial
            iterations = iterations + 1
            CALL eval_refine_point(kf, nc_loc, is_deg, vel_cur, inv_cur)
         END IF
         state = refine_converged
      END IF

      IF (.NOT. frame_valid) THEN
         CALL solve_at(work, kf, ispin, work%evals, nc_loc)
      END IF
      DEALLOCATE (trial_evals)

      IF (work%iw > 0) THEN
         IF (state == refine_converged) THEN
            WRITE (work%iw, '(T2,A,A,A,I0,A)') 'EFFECTIVE_MASS| ', label, ' refinement: converged in ', &
               iterations, ' steps'
         ELSE
            WRITE (work%iw, '(T2,A,A,A,I0,A)') 'EFFECTIVE_MASS| ', label, &
               ' refinement: best point after ', iterations, ' steps'
         END IF
      END IF
      CALL report_bands(work, kf, ispin, nc_loc, [iband], TRIM(label)//' (auto)')

   CONTAINS

! **************************************************************************************************
!> \brief Evaluate oriented band energy at kf_base + step_cart without overwriting work%evals.
!> \param kf_base Base fractional wave vector
!> \param step_cart Cartesian k-displacement in 1/bohr
!> \param kf_trial Wrapped fractional trial wave vector
!> \param f_trial Oriented band energy at kf_trial
! **************************************************************************************************
      SUBROUTINE eval_trial_energy(kf_base, step_cart, kf_trial, f_trial)
      REAL(dp), INTENT(IN)                               :: kf_base(3), step_cart(3)
      REAL(dp), INTENT(OUT)                              :: kf_trial(3), f_trial

         kf_trial = displace_fractional(kf_base, step_cart, work%cell%hmat, perd)
         CALL solve_at(work, kf_trial, ispin, trial_evals, nc_loc)
         frame_valid = .FALSE.
         IF (iband > nc_loc) CPABORT('EFFECTIVE_MASS target band exceeds retained basis states.')
         f_trial = sensef*trial_evals(iband)
      END SUBROUTINE eval_trial_energy

! **************************************************************************************************
!> \brief Probe +/- probe_step along each direction column in dirs and record the best descent.
!> \param kf_base Base fractional wave vector
!> \param dirs Normalized Cartesian probe directions (3, ndirs)
!> \param kf_best Best fractional probe point seen
!> \param f_best Lowest oriented energy seen among sufficient-descent probes
! **************************************************************************************************
      SUBROUTINE probe_directions(kf_base, dirs, kf_best, f_best)
      REAL(dp), INTENT(IN)                               :: kf_base(3), dirs(:, :)
      REAL(dp), INTENT(INOUT)                            :: kf_best(3), f_best

      INTEGER                                            :: idir, sgn
      REAL(dp)                                           :: f_t, f_thresh, kf_t(3)

         f_thresh = f_cur - refine_decrease_tol*MAX(1.0_dp, ABS(f_cur))
         DO idir = 1, SIZE(dirs, 2)
            DO sgn = -1, 1, 2
               CALL eval_trial_energy(kf_base, REAL(sgn, dp)*probe_step*dirs(1:3, idir), kf_t, f_t)
               IF (f_t < f_thresh .AND. f_t < f_best) THEN
                  f_best = f_t
                  kf_best = kf_t
               END IF
            END DO
         END DO
      END SUBROUTINE probe_directions

! **************************************************************************************************
!> \brief Solve at kf_pt and obtain velocity/curvature when the target band is non-degenerate.
!> \param kf_pt Fractional wave vector
!> \param nc_p Number of retained states
!> \param deg_p True if iband is degenerate
!> \param vel_p Velocity of iband
!> \param inv_p Inverse mass tensor of iband
! **************************************************************************************************
      SUBROUTINE eval_refine_point(kf_pt, nc_p, deg_p, vel_p, inv_p)
      REAL(dp), INTENT(IN)                               :: kf_pt(3)
      INTEGER, INTENT(OUT)                               :: nc_p
      LOGICAL, INTENT(OUT)                               :: deg_p
      REAL(dp), INTENT(OUT)                              :: vel_p(3), inv_p(3, 3)

      COMPLEX(dp), ALLOCATABLE                           :: a_s(:, :, :, :), b_s(:, :, :, :, :), &
                                                            s_s(:, :, :, :)
      INTEGER                                            :: ihi_r, ilo_r

         CALL solve_at(work, kf_pt, ispin, work%evals, nc_p)
         frame_valid = .TRUE.
         IF (iband > nc_p) CPABORT('EFFECTIVE_MASS target band exceeds retained basis states.')
         CALL degenerate_range(work%evals(1:nc_p), iband, work%eps_degeneracy, ilo_r, ihi_r)
         deg_p = (ihi_r > ilo_r)
         IF (.NOT. deg_p) THEN
            CALL analyze_subspaces(work, kf_pt, ispin, nc_p, [iband], [iband], a_s, b_s, s_s)
            CALL extract_nondeg_tensor(a_s, b_s, s_s, 1, vel_p, inv_p)
            DEALLOCATE (a_s, b_s, s_s)
         ELSE
            vel_p = 0.0_dp
            inv_p = 0.0_dp
         END IF
      END SUBROUTINE eval_refine_point
   END SUBROUTINE refine_edge

! **************************************************************************************************
!> \brief Select and report the bands of one explicitly requested k-point.
!> \param work Shared effective-mass workspace
!> \param kfrac Fractional k-point coordinates
!> \param n_occ Occupied band counts per spin channel
!> \param band_indices Optional explicit BAND_INDICES list (null if unset)
!> \param energy_window Energy window around HOMO/LUMO when band_indices is unset
! **************************************************************************************************
   SUBROUTINE report_explicit(work, kfrac, n_occ, band_indices, energy_window)
      TYPE(emass_work_type), INTENT(INOUT)               :: work
      REAL(dp), INTENT(IN)                               :: kfrac(3)
      INTEGER, INTENT(IN)                                :: n_occ(2)
      INTEGER, DIMENSION(:), POINTER                     :: band_indices
      REAL(dp), INTENT(IN)                               :: energy_window

      INTEGER                                            :: iband, ispin, nc_loc, nselected
      INTEGER, ALLOCATABLE                               :: selected(:)

      DO ispin = 1, work%nspin
         CALL solve_at(work, kfrac, ispin, work%evals, nc_loc)
         IF (ASSOCIATED(band_indices)) THEN
            IF (ANY(band_indices > nc_loc)) THEN
               CPABORT('BAND_INDICES entry exceeds retained basis states')
            END IF
            CALL report_bands(work, kfrac, ispin, nc_loc, band_indices, 'explicit')
         ELSE
            IF (n_occ(ispin) + 1 > nc_loc) THEN
               CPABORT('Band-edge selection exceeds retained basis states')
            END IF
            IF (energy_window <= 0.0_dp) THEN
               CALL report_bands(work, kfrac, ispin, nc_loc, [n_occ(ispin), n_occ(ispin) + 1], 'explicit')
            ELSE
               ALLOCATE (selected(nc_loc))
               nselected = 0
               DO iband = 1, nc_loc
                  IF (work%evals(iband) >= work%evals(n_occ(ispin)) - energy_window .AND. &
                      work%evals(iband) <= work%evals(n_occ(ispin) + 1) + energy_window) THEN
                     nselected = nselected + 1
                     selected(nselected) = iband
                  END IF
               END DO
               CALL report_bands(work, kfrac, ispin, nc_loc, selected(1:nselected), 'explicit')
               DEALLOCATE (selected)
            END IF
         END IF
      END DO
   END SUBROUTINE report_explicit

! **************************************************************************************************
!> \brief Format one effective mass 1/inv_m, returning right-aligned INF below the flatness tol.
!> \param inv_m Inverse effective mass eigenvalue (m_e/m*)
!> \param tol Flat-band threshold
!> \return Fixed-width 14-character formatted token
! **************************************************************************************************
   ELEMENTAL FUNCTION format_mass(inv_m, tol) RESULT(str)
      REAL(dp), INTENT(IN)                               :: inv_m, tol
      CHARACTER(LEN=14)                                  :: str

      IF (ABS(inv_m) > tol) THEN
         WRITE (str, '(ES14.6)') 1.0_dp/inv_m
      ELSE
         str = '           INF'
      END IF
   END FUNCTION format_mass

! **************************************************************************************************
!> \brief Evaluate and print velocities, mass tensors, and Luttinger-Kohn directional masses
!>        for the requested bands at (kfrac, ispin); requires work%cfm_frame at (kfrac, ispin).
!> \param work Shared effective-mass workspace
!> \param kfrac Fractional k-point coordinates of the analysis
!> \param ispin Spin channel of the analysis
!> \param nc Number of retained states in work%cfm_frame
!> \param bands Band indices to report
!> \param label Target label used in the header
! **************************************************************************************************
   SUBROUTINE report_bands(work, kfrac, ispin, nc, bands, label)
      TYPE(emass_work_type), INTENT(INOUT)               :: work
      REAL(dp), INTENT(IN)                               :: kfrac(3)
      INTEGER, INTENT(IN)                                :: ispin, nc, bands(:)
      CHARACTER(LEN=*), INTENT(IN)                       :: label

      CHARACTER(LEN=14)                                  :: mass_str(3)
      CHARACTER(LEN=14), ALLOCATABLE                     :: mass_tok(:)
      COMPLEX(dp), ALLOCATABLE                           :: a_sub(:, :, :, :), aop(:, :), &
                                                            b_sub(:, :, :, :, :), bop(:, :), &
                                                            s_sub(:, :, :, :), sop(:, :)
      INTEGER                                            :: a, b, carrier, d, ib, iband, idir, ihi, &
                                                            ilo, isub, nsub
      INTEGER, ALLOCATABLE                               :: sub_hi(:), sub_lo(:)
      LOGICAL                                            :: is_stationary, linear_split
      REAL(dp)                                           :: dirs(3, 3), e0, inv(3, 3), inv_eigs(3), &
                                                            inv_tol, kcart(3), qhat(3), vel(3)
      REAL(dp), ALLOCATABLE                              :: inv_mass_eigs(:), velocity_eigs(:)

      IF (work%iw > 0) THEN
         kcart(1:3) = twopi*MATMUL(kfrac(1:3), work%cell%h_inv)
         WRITE (work%iw, '(T2,A,A,A,3(1X,F12.6),A,I0)') 'EFFECTIVE_MASS| ', TRIM(label), &
            ' k-point (fractional):', kfrac(1:3), '  spin ', ispin
         WRITE (work%iw, '(T2,A,3(1X,F12.6))') 'EFFECTIVE_MASS| Cartesian k [1/bohr]:        ', kcart
      END IF
      IF (SIZE(bands) == 0) RETURN

      ALLOCATE (sub_lo(SIZE(bands)), sub_hi(SIZE(bands)))
      nsub = 0
      DO ib = 1, SIZE(bands)
         iband = bands(ib)
         CPASSERT(iband >= 1 .AND. iband <= nc)
         CALL degenerate_range(work%evals(1:nc), iband, work%eps_degeneracy, ilo, ihi)
         IF (.NOT. ANY(sub_lo(1:nsub) == ilo)) THEN
            nsub = nsub + 1
            sub_lo(nsub) = ilo
            sub_hi(nsub) = ihi
         END IF
      END DO

      CALL analyze_subspaces(work, kfrac, ispin, nc, sub_lo(1:nsub), sub_hi(1:nsub), &
                             a_sub, b_sub, s_sub)

      DO isub = 1, nsub
         ilo = sub_lo(isub)
         ihi = sub_hi(isub)
         d = ihi - ilo + 1
         IF (d == 1) THEN
            CALL extract_nondeg_tensor(a_sub, b_sub, s_sub, isub, vel, inv)
            IF (work%iw > 0) THEN
               CALL principal_masses(inv, inv_eigs, dirs, carrier)
               inv_tol = carrier_tol*MAX(1.0_dp, MAXVAL(ABS(inv_eigs)))
               mass_str(1:3) = format_mass(inv_eigs(1:3), inv_tol)
               WRITE (work%iw, '(T2,A,1X,I6,2(1X,F14.8))') 'EFFECTIVE_MASS| Band, energy [Ha,eV]:', ilo, &
                  work%evals(ilo), work%evals(ilo)*evolt
               WRITE (work%iw, '(T2,A,*(1X,F13.7))') 'EFFECTIVE_MASS| Velocity [a.u.]:        ', vel
               WRITE (work%iw, '(T2,A)') 'EFFECTIVE_MASS| Inverse mass tensor (m_e/m*), rows x,y,z:'
               WRITE (work%iw, '(T2,A,*(1X,F13.7))') 'EFFECTIVE_MASS|   ', inv(1:3, 1)
               WRITE (work%iw, '(T2,A,*(1X,F13.7))') 'EFFECTIVE_MASS|   ', inv(1:3, 2)
               WRITE (work%iw, '(T2,A,*(1X,F13.7))') 'EFFECTIVE_MASS|   ', inv(1:3, 3)
               WRITE (work%iw, '(T2,A,*(1X,F13.7))') 'EFFECTIVE_MASS| Principal inv masses (m_e/m*):', inv_eigs
               WRITE (work%iw, '(T2,A,3(1X,A14),2X,A)') 'EFFECTIVE_MASS| Principal masses (m*/m_e):', &
                  mass_str, TRIM(carrier_labels(carrier))
               DO idir = 1, 3
                  WRITE (work%iw, '(T2,A,I2,A,3(1X,F10.6))') 'EFFECTIVE_MASS| Principal direction ', idir, &
                     ' (x,y,z):', dirs(1:3, idir)
               END DO
            END IF
         ELSE
            e0 = SUM(work%evals(ilo:ihi))/REAL(d, dp)
            IF (work%iw > 0) THEN
               WRITE (work%iw, '(T2,A,I0,A,I0,A,2(1X,F14.8))') 'EFFECTIVE_MASS| Bands ', ilo, '-', ihi, &
                  ', centroid energy [Ha,eV]:', e0, e0*evolt
               WRITE (work%iw, '(T2,A,I0,A)') 'EFFECTIVE_MASS| Degenerate subspace of size ', d, &
                  ' within EPS_DEGENERACY: a unique 3x3 Hessian does not exist.'
               WRITE (work%iw, '(T2,A)') 'EFFECTIVE_MASS| Luttinger-Kohn directional inverse masses (m_e/m*):'
            END IF
            ALLOCATE (aop(d, d), bop(d, d), sop(d, d), mass_tok(d))
            DO idir = 1, SIZE(lk_dir_idx)
               qhat(1:3) = refine_probe_dirs(1:3, lk_dir_idx(idir))
               aop(:, :) = z_zero
               bop(:, :) = z_zero
               sop(:, :) = z_zero
               DO a = 1, 3
                  aop(:, :) = aop + qhat(a)*a_sub(1:d, 1:d, a, isub)
                  sop(:, :) = sop + qhat(a)*s_sub(1:d, 1:d, a, isub)
                  DO b = 1, 3
                     bop(:, :) = bop + qhat(a)*qhat(b)*b_sub(1:d, 1:d, a, b, isub)
                  END DO
               END DO
               CALL kato_rellich_reduce(aop, bop, sop, velocity_eigs, inv_mass_eigs, linear_split)
               IF (work%iw > 0) THEN
                  is_stationary = (.NOT. linear_split) .AND. (MAXVAL(ABS(velocity_eigs)) <= refine_grad_tol)
                  WRITE (work%iw, '(T2,A,3(1X,F7.3),A,L1,A)') 'EFFECTIVE_MASS| Cartesian dir (x,y,z):', &
                     qhat, '  [linear_split=', linear_split, ']'
                  WRITE (work%iw, '(T2,A,*(1X,F13.7))') 'EFFECTIVE_MASS|   Velocities [a.u.]:      ', &
                     velocity_eigs
                  IF (is_stationary) THEN
                     WRITE (work%iw, '(T2,A,*(1X,F13.7))') 'EFFECTIVE_MASS|   Dir. inv masses (m_e/m*):', &
                        inv_mass_eigs
                     inv_tol = carrier_tol*MAX(1.0_dp, MAXVAL(ABS(inv_mass_eigs)))
                     mass_tok(1:d) = format_mass(inv_mass_eigs(1:d), inv_tol)
                     WRITE (work%iw, '(T2,A,*(1X,A14))') 'EFFECTIVE_MASS|   Dir. masses (m*/m_e):   ', &
                        mass_tok
                  ELSE
                     WRITE (work%iw, '(T2,A,*(1X,F13.7))') &
                        'EFFECTIVE_MASS|   Dir. curvature d2E/dk2 [a.u.]:', inv_mass_eigs
                  END IF
               END IF
               DEALLOCATE (velocity_eigs, inv_mass_eigs)
            END DO
            DEALLOCATE (aop, bop, sop, mass_tok)
         END IF
      END DO

      DEALLOCATE (sub_lo, sub_hi, a_sub, b_sub, s_sub)
   END SUBROUTINE report_bands

END MODULE qs_effective_mass
