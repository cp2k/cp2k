!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to Wannier90 code
!> \par History
!>      06.2016 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE qs_wannier90
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_deallocate_matrix, dbcsr_desymmetrize, dbcsr_p_type, dbcsr_set, &
        dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_get_element,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_methods,                  ONLY: kpoint_env_initialize,&
                                              kpoint_init_cell_index,&
                                              kpoint_initialize,&
                                              kpoint_initialize_mo_set,&
                                              kpoint_initialize_mos,&
                                              rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kind_rotmat_type,&
                                              kpoint_create,&
                                              kpoint_env_type,&
                                              kpoint_release,&
                                              kpoint_sym_type,&
                                              kpoint_type
   USE machine,                         ONLY: m_timestamp,&
                                              timestamp_length
   USE mathconstants,                   ONLY: twopi
   USE mathlib,                         ONLY: diag_complex
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_release,&
                                              qs_environment_type
   USE qs_gamma2kp,                     ONLY: create_kp_from_gamma
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_moments,                      ONLY: build_berry_kpoint_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_scf_diagonalization,          ONLY: do_general_diag_kp
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type
   USE wannier90,                       ONLY: wannier_setup
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_wannier90'
   INTEGER, PARAMETER, PRIVATE :: w90_kpoints_mp_grid = 0, &
                                  w90_kpoints_scf = 1

   TYPE berry_matrix_type
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER      :: sinmat => NULL(), cosmat => NULL()
   END TYPE berry_matrix_type

   PUBLIC :: wannier90_interface, prepare_wannier90_scf_mos

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param input ...
!> \param logger ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE wannier90_interface(input, logger, qs_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'wannier90_interface'

      INTEGER                                            :: handle, iw
      LOGICAL                                            :: explicit
      TYPE(section_vals_type), POINTER                   :: w_input

      !--------------------------------------------------------------------------------------------!

      CALL timeset(routineN, handle)
      w_input => section_vals_get_subs_vals(section_vals=input, &
                                            subsection_name="DFT%PRINT%WANNIER90")
      CALL section_vals_get(w_input, explicit=explicit)
      IF (explicit) THEN

         iw = cp_logger_get_default_io_unit(logger)

         IF (iw > 0) THEN
            WRITE (iw, '(/,T2,A)') &
               '!-----------------------------------------------------------------------------!'
            WRITE (iw, '(T32,A)') "Interface to Wannier90"
            WRITE (iw, '(T2,A)') &
               '!-----------------------------------------------------------------------------!'
         END IF

         CALL wannier90_files(qs_env, w_input, iw)

         IF (iw > 0) THEN
            WRITE (iw, '(/,T2,A)') &
               '!--------------------------------End of Wannier90-----------------------------!'
         END IF
      END IF
      CALL timestop(handle)

   END SUBROUTINE wannier90_interface

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param input ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE wannier90_files(qs_env, input, iw)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input
      INTEGER, INTENT(IN)                                :: iw

      INTEGER, PARAMETER                                 :: num_nnmax = 12

      CHARACTER(len=2)                                   :: asym
      CHARACTER(len=20), ALLOCATABLE, DIMENSION(:)       :: atom_symbols
      CHARACTER(len=default_string_length)               :: filename, input_kp_scheme, reuse_reason, &
                                                            seed_name
      CHARACTER(LEN=timestamp_length)                    :: timestamp
      INTEGER :: aligned_degenerate_blocks, aligned_degenerate_max_size, i, i_rep, ib, ib1, ib2, &
         ibs, ik, ik2, ikk, ikpgr, ispin, iunit, ix, iy, iz, k, kpoints_source, n_rep, nadd, nao, &
         nbs, nexcl, nkp, nmo, nntot, nspins, num_atoms, num_bands, num_bands_tot, num_kpts, &
         num_wann
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: exclude_bands
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: nblist, nnlist
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: nncell
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(3)                              :: input_nkp_grid, mp_grid
      INTEGER, DIMENSION(:), POINTER                     :: invals
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL :: diis_step, do_kpoints, full_mesh_diagonalized, gamma_only, input_full_grid, &
         input_gamma_centered, input_kpoint_symmetry, mp_grid_explicit, mp_grid_valid, my_kpgrp, &
         mygrp, reuse_scf_mos, reused_scf_mos, spinors, use_bloch_phases, validate_reuse_ok, &
         validate_reuse_scf_mos
      REAL(KIND=dp) :: aligned_degenerate_min_svalue, cmmn, gauge_arg, gauge_imag, gauge_real, &
         gauge_tmp, ksign, reuse_candidate_deviation, reuse_candidate_metric_deviation, &
         reuse_candidate_min_svalue, reuse_candidate_residual, rmmn, &
         validation_eigenvalue_deviation, validation_min_svalue, validation_subspace_deviation, &
         wkp_ref
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atoms_cart, b_latt, kpt_latt
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: reference_eigenvalues
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: reference_mo_imag, reference_mo_real
      REAL(KIND=dp), DIMENSION(3)                        :: bvec, input_kp_shift, phase_center
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_inv, real_lattice, recip_lattice
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, wkp, wkp_source
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp, xkp_source
      TYPE(berry_matrix_type), DIMENSION(:), POINTER     :: berry_matrix
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: fmk1_cfm, fmk2_cfm, mmn_cfm, omat_cfm, &
                                                            tmp_cfm
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_ao, matrix_struct_mmn, &
                                                            matrix_struct_work
      TYPE(cp_fm_type)                                   :: mat_imag, mat_real, mmn_imag, mmn_real
      TYPE(cp_fm_type), DIMENSION(2)                     :: fmk1, fmk2
      TYPE(cp_fm_type), POINTER                          :: fmdummy, fmi, fmr
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      TYPE(dbcsr_type), POINTER                          :: cmatrix, cmatrix_full, rmatrix, &
                                                            rmatrix_full
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoint, qs_kpoint
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_environment_type), POINTER                 :: qs_env_kp
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      !--------------------------------------------------------------------------------------------!

      ! add code for exclude_bands and projectors

      ! generate all arrays needed for the setup call
      CALL section_vals_val_get(input, "SEED_NAME", c_val=seed_name)
      CALL section_vals_val_get(input, "MP_GRID", i_vals=invals, explicit=mp_grid_explicit)
      CALL section_vals_val_get(input, "KPOINTS_SOURCE", i_val=kpoints_source)
      CALL section_vals_val_get(input, "WANNIER_FUNCTIONS", i_val=num_wann)
      CALL section_vals_val_get(input, "ADDED_MOS", i_val=nadd)
      CALL section_vals_val_get(input, "REUSE_SCF_MOS", l_val=reuse_scf_mos)
      CALL section_vals_val_get(input, "VALIDATE_REUSE_SCF_MOS", l_val=validate_reuse_scf_mos)
      CALL section_vals_val_get(input, "USE_BLOCH_PHASES", l_val=use_bloch_phases)
      reuse_scf_mos = reuse_scf_mos .AND. kpoints_source == w90_kpoints_scf
      validate_reuse_scf_mos = validate_reuse_scf_mos .AND. reuse_scf_mos
      mp_grid(1:3) = invals(1:3)
      ! excluded bands
      CALL section_vals_val_get(input, "EXCLUDE_BANDS", n_rep_val=n_rep)
      nexcl = 0
      DO i_rep = 1, n_rep
         CALL section_vals_val_get(input, "EXCLUDE_BANDS", i_rep_val=i_rep, i_vals=invals)
         nexcl = nexcl + SIZE(invals)
      END DO
      IF (nexcl > 0) THEN
         ALLOCATE (exclude_bands(nexcl))
         nexcl = 0
         DO i_rep = 1, n_rep
            CALL section_vals_val_get(input, "EXCLUDE_BANDS", i_rep_val=i_rep, i_vals=invals)
            exclude_bands(nexcl + 1:nexcl + SIZE(invals)) = invals(:)
            nexcl = nexcl + SIZE(invals)
         END DO
      END IF
      !
      ! lattice -> Angstrom
      CALL get_qs_env(qs_env, cell=cell)
      CALL get_cell(cell, h=real_lattice, h_inv=h_inv)
      ! k-points
      CALL get_qs_env(qs_env, particle_set=particle_set)
      CALL get_qs_env(qs_env, para_env=para_env)
      phase_center = 0.0_dp
      DO i = 1, SIZE(particle_set)
         phase_center(1:3) = phase_center(1:3) + MATMUL(h_inv, particle_set(i)%r)
      END DO
      phase_center(1:3) = phase_center(1:3)/REAL(SIZE(particle_set), KIND=dp)
      phase_center(1:3) = phase_center(1:3) - FLOOR(phase_center(1:3))
      recip_lattice(1:3, 1:3) = h_inv(1:3, 1:3)
      real_lattice(1:3, 1:3) = angstrom*real_lattice(1:3, 1:3)
      recip_lattice(1:3, 1:3) = (twopi/angstrom)*TRANSPOSE(recip_lattice(1:3, 1:3))
      NULLIFY (kpoint, qs_kpoint, xkp, wkp, xkp_source, wkp_source)
      CALL get_qs_env(qs_env, do_kpoints=do_kpoints, kpoints=qs_kpoint)
      input_kpoint_symmetry = .FALSE.
      input_full_grid = .FALSE.
      input_kp_scheme = ""
      IF (do_kpoints .AND. ASSOCIATED(qs_kpoint)) THEN
         CALL get_kpoint_info(qs_kpoint, kp_scheme=input_kp_scheme, nkp_grid=input_nkp_grid, &
                              kp_shift=input_kp_shift, symmetry=input_kpoint_symmetry, &
                              full_grid=input_full_grid, gamma_centered=input_gamma_centered, &
                              nkp=nkp, xkp=xkp, wkp=wkp)
      END IF
      CALL kpoint_create(kpoint)

      SELECT CASE (kpoints_source)
      CASE (w90_kpoints_mp_grid)
         num_kpts = mp_grid(1)*mp_grid(2)*mp_grid(3)
         ALLOCATE (kpt_latt(3, num_kpts))
         kpoint%kp_scheme = "MONKHORST-PACK"
         kpoint%symmetry = .FALSE.
         kpoint%nkp_grid(1:3) = mp_grid(1:3)
         kpoint%verbose = .FALSE.
         kpoint%full_grid = .TRUE.
         kpoint%eps_geo = 1.0e-6_dp
         kpoint%use_real_wfn = .FALSE.
         kpoint%parallel_group_size = para_env%num_pe
         i = 0
         DO ix = 0, mp_grid(1) - 1
            DO iy = 0, mp_grid(2) - 1
               DO iz = 0, mp_grid(3) - 1
                  i = i + 1
                  kpt_latt(1, i) = REAL(ix, KIND=dp)/REAL(mp_grid(1), KIND=dp)
                  kpt_latt(2, i) = REAL(iy, KIND=dp)/REAL(mp_grid(2), KIND=dp)
                  kpt_latt(3, i) = REAL(iz, KIND=dp)/REAL(mp_grid(3), KIND=dp)
               END DO
            END DO
         END DO
         kpoint%nkp = num_kpts
         ALLOCATE (kpoint%xkp(3, num_kpts), kpoint%wkp(num_kpts))
         kpoint%wkp(:) = 1._dp/REAL(num_kpts, KIND=dp)
         DO i = 1, num_kpts
            kpoint%xkp(1:3, i) = (angstrom/twopi)*MATMUL(recip_lattice, kpt_latt(:, i))
         END DO

      CASE (w90_kpoints_scf)
         IF (.NOT. do_kpoints .OR. .NOT. ASSOCIATED(qs_kpoint)) THEN
            CPABORT("WANNIER90%KPOINTS_SOURCE SCF requires an active DFT%KPOINTS section.")
         END IF
         SELECT CASE (TRIM(input_kp_scheme))
         CASE ("GAMMA")
            mp_grid(:) = 1
            num_kpts = 1
            ALLOCATE (kpt_latt(3, num_kpts))
            kpt_latt(1:3, 1) = 0.0_dp
            kpoint%kp_scheme = "GAMMA"
            kpoint%symmetry = .FALSE.
            kpoint%verbose = .FALSE.
            kpoint%full_grid = .TRUE.
            kpoint%eps_geo = 1.0e-6_dp
            kpoint%use_real_wfn = .FALSE.
            kpoint%parallel_group_size = para_env%num_pe
            kpoint%nkp = num_kpts
            ALLOCATE (kpoint%xkp(3, num_kpts), kpoint%wkp(num_kpts))
            kpoint%xkp(1:3, 1) = 0.0_dp
            kpoint%wkp(1) = 1.0_dp

         CASE ("MONKHORST-PACK", "MACDONALD")
            mp_grid(1:3) = input_nkp_grid(1:3)
            kpoint%kp_scheme = input_kp_scheme
            kpoint%symmetry = .FALSE.
            kpoint%nkp_grid(1:3) = input_nkp_grid(1:3)
            kpoint%kp_shift(1:3) = input_kp_shift(1:3)
            kpoint%gamma_centered = input_gamma_centered
            kpoint%verbose = .FALSE.
            kpoint%full_grid = .TRUE.
            kpoint%eps_geo = 1.0e-6_dp
            kpoint%use_real_wfn = .FALSE.
            kpoint%parallel_group_size = para_env%num_pe
            CALL kpoint_initialize(kpoint, particle_set, cell)
            num_kpts = kpoint%nkp
            ALLOCATE (kpt_latt(3, num_kpts))
            kpt_latt(1:3, 1:num_kpts) = kpoint%xkp(1:3, 1:num_kpts)
            IF (input_kpoint_symmetry .AND. .NOT. input_full_grid .AND. iw > 0) THEN
               WRITE (iw, '(T2,A)') &
                  "WANNIER90| SCF k-points are symmetry-reduced; regenerating the full SCF mesh."
               IF (reuse_scf_mos) THEN
                  WRITE (iw, '(T2,A)') &
                     "WANNIER90| CP2K will try to reconstruct the full-mesh MOs from the SCF orbitals."
               ELSE
                  WRITE (iw, '(T2,A)') &
                     "WANNIER90| The full exported mesh is diagonalized for the Wannier90 files."
               END IF
            END IF

         CASE ("GENERAL")
            IF (ASSOCIATED(qs_kpoint%xkp_input)) THEN
               xkp_source => qs_kpoint%xkp_input
               wkp_source => qs_kpoint%wkp_input
            ELSE
               xkp_source => xkp
               wkp_source => wkp
            END IF
            IF (.NOT. ASSOCIATED(xkp_source) .OR. .NOT. ASSOCIATED(wkp_source)) THEN
               CPABORT("Could not access the SCF GENERAL k-point set for the Wannier90 export.")
            END IF
            num_kpts = SIZE(wkp_source)
            ALLOCATE (kpt_latt(3, num_kpts))
            kpt_latt(1:3, 1:num_kpts) = xkp_source(1:3, 1:num_kpts)
            IF (mp_grid_explicit) THEN
               IF (mp_grid(1)*mp_grid(2)*mp_grid(3) /= num_kpts) THEN
                  CPABORT("WANNIER90%MP_GRID must contain exactly as many points as the SCF GENERAL mesh.")
               END IF
            ELSE
               CALL infer_wannier_mp_grid(kpt_latt, mp_grid, mp_grid_valid)
               IF (.NOT. mp_grid_valid) THEN
                  CPABORT("Could not infer WANNIER90%MP_GRID from the SCF GENERAL mesh.")
               END IF
            END IF
            wkp_ref = 1.0_dp/REAL(num_kpts, KIND=dp)
            DO i = 1, num_kpts
               IF (ABS(wkp_source(i) - wkp_ref) > 1.0e-10_dp) THEN
                  CPABORT("WANNIER90%KPOINTS_SOURCE SCF requires equally weighted GENERAL k-points.")
               END IF
            END DO
            kpoint%kp_scheme = "GENERAL"
            kpoint%symmetry = .FALSE.
            kpoint%nkp_grid(1:3) = mp_grid(1:3)
            kpoint%verbose = .FALSE.
            kpoint%full_grid = .TRUE.
            kpoint%eps_geo = 1.0e-6_dp
            kpoint%use_real_wfn = .FALSE.
            kpoint%parallel_group_size = para_env%num_pe
            kpoint%nkp = num_kpts
            ALLOCATE (kpoint%xkp(3, num_kpts), kpoint%wkp(num_kpts))
            kpoint%xkp(1:3, 1:num_kpts) = xkp_source(1:3, 1:num_kpts)
            kpoint%wkp(1:num_kpts) = wkp_ref
            IF (input_kpoint_symmetry .AND. .NOT. input_full_grid .AND. iw > 0) THEN
               WRITE (iw, '(T2,A)') &
                  "WANNIER90| SCF k-points are symmetry-reduced; using the full input GENERAL mesh."
               IF (reuse_scf_mos) THEN
                  WRITE (iw, '(T2,A)') &
                     "WANNIER90| CP2K will try to reconstruct the full-mesh MOs from the SCF orbitals."
               ELSE
                  WRITE (iw, '(T2,A)') &
                     "WANNIER90| The full exported mesh is diagonalized for the Wannier90 files."
               END IF
            END IF

         CASE DEFAULT
            CPABORT("WANNIER90%KPOINTS_SOURCE SCF does not support this DFT%KPOINTS scheme.")
         END SELECT
      CASE DEFAULT
         CPABORT("Unknown WANNIER90%KPOINTS_SOURCE setting.")
      END SELECT
      ! number of bands in calculation
      CALL get_qs_env(qs_env, mos=mos)
      CALL get_mo_set(mo_set=mos(1), nao=nao, nmo=num_bands_tot)
      num_bands_tot = MIN(nao, num_bands_tot + nadd)
      num_bands = num_bands_tot
      IF (use_bloch_phases .AND. num_wann /= num_bands) THEN
         CPABORT("WANNIER90%USE_BLOCH_PHASES requires WANNIER_FUNCTIONS to match the number of bands.")
      END IF
      num_atoms = SIZE(particle_set)
      ALLOCATE (atoms_cart(3, num_atoms))
      ALLOCATE (atom_symbols(num_atoms))
      DO i = 1, num_atoms
         atoms_cart(1:3, i) = particle_set(i)%r(1:3)
         CALL get_atomic_kind(particle_set(i)%atomic_kind, element_symbol=asym)
         atom_symbols(i) = asym
      END DO
      gamma_only = .FALSE.
      spinors = .FALSE.
      ! output
      ALLOCATE (nnlist(num_kpts, num_nnmax))
      ALLOCATE (nncell(3, num_kpts, num_nnmax))
      nnlist(:, :) = 0
      nncell(:, :, :) = 0
      nntot = 0

      IF (iw > 0) THEN
         ! setup
         CALL wannier_setup(mp_grid, num_kpts, real_lattice, recip_lattice, &
                            kpt_latt, nntot, nnlist, nncell, iw)
      END IF

      CALL get_qs_env(qs_env, para_env=para_env)
      CALL para_env%sum(nntot)
      CALL para_env%sum(nnlist)
      CALL para_env%sum(nncell)

      IF (para_env%is_source()) THEN
         ! Write the Wannier90 input file "seed_name.win"
         WRITE (filename, '(A,A)') TRIM(seed_name), ".win"
         CALL open_file(filename, unit_number=iunit, file_status="UNKNOWN", file_action="WRITE")
         !
         CALL m_timestamp(timestamp)
         WRITE (iunit, "(A)") "! Wannier90 input file generated by CP2K "
         WRITE (iunit, "(A,/)") "! Creation date "//timestamp
         !
         WRITE (iunit, "(A,I5)") "num_wann     = ", num_wann
         IF (num_bands /= num_wann .OR. use_bloch_phases) THEN
            WRITE (iunit, "(A,I5)") "num_bands    = ", num_bands
         END IF
         IF (use_bloch_phases) THEN
            ! Keep the external Wannier90 projection matrix fully defined for
            ! complete-band Bloch-phase subspaces by writing explicit identity projections.
            WRITE (iunit, "(A)") "! CP2K writes identity projections for Bloch-phase complete subspaces."
         END IF
         WRITE (iunit, "(/,A,/)") "length_unit  = bohr "
         WRITE (iunit, "(/,A,/)") "! System"
         WRITE (iunit, "(/,A)") "begin unit_cell_cart"
         WRITE (iunit, "(A)") "bohr"
         DO i = 1, 3
            WRITE (iunit, "(3F12.6)") cell%hmat(i, 1:3)
         END DO
         WRITE (iunit, "(A,/)") "end unit_cell_cart"
         WRITE (iunit, "(/,A)") "begin atoms_cart"
         WRITE (iunit, "(A)") "bohr"
         DO i = 1, num_atoms
            WRITE (iunit, "(A,3F15.10)") atom_symbols(i), atoms_cart(1:3, i)
         END DO
         WRITE (iunit, "(A,/)") "end atoms_cart"
         WRITE (iunit, "(/,A,/)") "! Kpoints"
         WRITE (iunit, "(/,A,3I6/)") "mp_grid      = ", mp_grid(1:3)
         WRITE (iunit, "(A)") "begin kpoints"
         DO i = 1, num_kpts
            WRITE (iunit, "(3F12.6)") kpt_latt(1:3, i)
         END DO
         WRITE (iunit, "(A)") "end kpoints"
         CALL close_file(iunit)
         IF (use_bloch_phases) THEN
            WRITE (filename, '(A,A)') TRIM(seed_name), ".amn"
            CALL open_file(filename, unit_number=iunit, file_status="UNKNOWN", file_action="WRITE")
            WRITE (iunit, "(A)") "! Wannier90 identity projections generated by CP2K"
            WRITE (iunit, "(3I8)") num_bands, num_kpts, num_wann
            DO ik = 1, num_kpts
               DO ib2 = 1, num_wann
                  DO ib1 = 1, num_bands
                     IF (ib1 == ib2) THEN
                        WRITE (iunit, "(3I8,2E30.14)") ib1, ib2, ik, 1.0_dp, 0.0_dp
                     ELSE
                        WRITE (iunit, "(3I8,2E30.14)") ib1, ib2, ik, 0.0_dp, 0.0_dp
                     END IF
                  END DO
               END DO
            END DO
            CALL close_file(iunit)
         END IF
      ELSE
         iunit = -1
      END IF

      ! calculate bands
      NULLIFY (qs_env_kp)
      IF (kpoints_source == w90_kpoints_mp_grid .AND. input_kpoint_symmetry .AND. iw > 0) THEN
         WRITE (iw, '(T2,A)') &
            "WANNIER90| Atomic k-point symmetry from the SCF calculation is not reused."
         WRITE (iw, '(T2,A)') &
            "WANNIER90| A full Monkhorst-Pack grid is generated for the Wannier90 interface."
      END IF
      IF (do_kpoints) THEN
         ! we already do kpoints
         qs_env_kp => qs_env
      ELSE
         ! we start from gamma point only
         ALLOCATE (qs_env_kp)
         CALL create_kp_from_gamma(qs_env, qs_env_kp)
      END IF
      IF (iw > 0) THEN
         WRITE (unit=iw, FMT="(/,T2,A)") "Start K-Point Calculation ..."
      END IF
      CALL get_qs_env(qs_env=qs_env_kp, para_env=para_env, blacs_env=blacs_env)
      CALL kpoint_env_initialize(kpoint, para_env, blacs_env)
      CALL kpoint_initialize_mos(kpoint, mos, nadd)
      CALL kpoint_initialize_mo_set(kpoint)
      !
      CALL get_qs_env(qs_env=qs_env_kp, sab_orb=sab_nl, dft_control=dft_control)
      CALL kpoint_init_cell_index(kpoint, sab_nl, para_env, dft_control%nimages)
      !
      CALL get_qs_env(qs_env=qs_env_kp, matrix_ks_kp=matrix_ks, matrix_s_kp=matrix_s, &
                      scf_env=scf_env, scf_control=scf_control)
      full_mesh_diagonalized = .FALSE.
      reused_scf_mos = .FALSE.
      reuse_reason = ""
      aligned_degenerate_blocks = 0
      aligned_degenerate_max_size = 0
      aligned_degenerate_min_svalue = 0.0_dp
      IF (reuse_scf_mos) THEN
         CALL get_kpoint_info(kpoint=kpoint, cell_to_index=cell_to_index)
         CALL do_general_diag_kp(matrix_ks, matrix_s, qs_kpoint, scf_env, scf_control, .FALSE., &
                                 diis_step)
         IF (validate_reuse_scf_mos) THEN
            IF (iw > 0) THEN
               WRITE (iw, '(T2,A)') &
                  "WANNIER90| Validating SCF MO reuse against a full-mesh diagonalization reference."
            END IF
            CALL do_general_diag_kp(matrix_ks, matrix_s, kpoint, scf_env, scf_control, .FALSE., diis_step)
            full_mesh_diagonalized = .TRUE.
            nspins = dft_control%nspins
            CALL save_wannier90_mo_snapshot(kpoint, nspins, para_env, reference_mo_real, &
                                            reference_mo_imag, reference_eigenvalues)
            CALL diagnose_wannier90_scf_reuse_candidates(kpoint, qs_kpoint, matrix_s, matrix_ks, &
                                                         cell_to_index, sab_nl, para_env, iw, &
                                                         reuse_candidate_deviation, &
                                                         reuse_candidate_min_svalue, &
                                                         reuse_candidate_metric_deviation, &
                                                         reuse_candidate_residual)
            IF (iw > 0 .AND. reuse_candidate_deviation < 1.0e100_dp) THEN
               WRITE (iw, '(T2,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') &
                  "WANNIER90| Best atom/AO candidate subspace deviation ", &
                  reuse_candidate_deviation, ", minimum singular value ", &
                  reuse_candidate_min_svalue, ", max metric deviation ", &
                  reuse_candidate_metric_deviation, ", max residual ", reuse_candidate_residual
            END IF
         END IF
         CALL prepare_wannier90_scf_mos(kpoint, qs_kpoint, matrix_s, matrix_ks, cell_to_index, &
                                        sab_nl, para_env, reused_scf_mos, reuse_reason, &
                                        aligned_degenerate_blocks, aligned_degenerate_max_size, &
                                        aligned_degenerate_min_svalue)
         IF (validate_reuse_scf_mos) THEN
            IF (reused_scf_mos) THEN
               CALL validate_wannier90_reused_mos(kpoint, matrix_s, cell_to_index, sab_nl, &
                                                  para_env, reference_mo_real, reference_mo_imag, &
                                                  reference_eigenvalues, validate_reuse_ok, &
                                                  validation_subspace_deviation, validation_min_svalue, &
                                                  validation_eigenvalue_deviation)
               IF (iw > 0) THEN
                  WRITE (iw, '(T2,A,ES10.3,A,ES10.3,A,ES10.3)') &
                     "WANNIER90| Reused MO validation: subspace deviation ", &
                     validation_subspace_deviation, ", minimum singular value ", &
                     validation_min_svalue, ", eigenvalue deviation ", &
                     validation_eigenvalue_deviation
               END IF
               IF (.NOT. validate_reuse_ok) THEN
                  reused_scf_mos = .FALSE.
                  WRITE (reuse_reason, "(A,ES10.3,A,ES10.3)") &
                     "validation failed: dS=", &
                     validation_subspace_deviation, ", dE=", validation_eigenvalue_deviation
               END IF
            END IF
            IF (.NOT. reused_scf_mos) THEN
               CALL restore_wannier90_mo_snapshot(kpoint, reference_mo_real, reference_mo_imag, &
                                                  reference_eigenvalues)
            END IF
         END IF
         IF (iw > 0) THEN
            IF (reused_scf_mos) THEN
               WRITE (iw, '(T2,A)') &
                  "WANNIER90| Reused SCF MO coefficients for the Wannier90 full k-point mesh."
               IF (use_bloch_phases) THEN
                  WRITE (iw, '(T2,A)') &
                     "WANNIER90| Wrote identity projections for Bloch-phase complete band subspaces."
                  WRITE (iw, '(T2,A,3F10.6)') &
                     "WANNIER90| Applied Bloch phase gauge to reused overlaps around fractional center", &
                     phase_center(1:3)
               END IF
               IF (aligned_degenerate_blocks > 0) THEN
                  WRITE (iw, '(T2,A,I0,A,I0,A,ES10.3)') &
                     "WANNIER90| Ritz-stabilized ", aligned_degenerate_blocks, &
                     " degenerate SCF MO subspace(s) with S(k),H(k); largest block has ", &
                     aligned_degenerate_max_size, " band(s), min metric eigenvalue ", &
                     aligned_degenerate_min_svalue
               END IF
            ELSE
               WRITE (iw, '(T2,A,A)') &
                  "WANNIER90| Could not reuse SCF MOs: ", TRIM(reuse_reason)
               WRITE (iw, '(T2,A)') &
                  "WANNIER90| Falling back to full-mesh diagonalization for the Wannier90 files."
            END IF
         END IF
      END IF
      IF (.NOT. reused_scf_mos .AND. .NOT. full_mesh_diagonalized) THEN
         CALL do_general_diag_kp(matrix_ks, matrix_s, kpoint, scf_env, scf_control, .FALSE., diis_step)
      END IF
      IF (ALLOCATED(reference_mo_real)) DEALLOCATE (reference_mo_real)
      IF (ALLOCATED(reference_mo_imag)) DEALLOCATE (reference_mo_imag)
      IF (ALLOCATED(reference_eigenvalues)) DEALLOCATE (reference_eigenvalues)
      !
      IF (iw > 0) THEN
         WRITE (iw, '(T69,A)') "... Finished"
      END IF
      !
      ! Calculate and print Overlaps
      !
      IF (para_env%is_source()) THEN
         WRITE (filename, '(A,A)') TRIM(seed_name), ".mmn"
         CALL open_file(filename, unit_number=iunit, file_status="UNKNOWN", file_action="WRITE")
         CALL m_timestamp(timestamp)
         WRITE (iunit, "(A)") "! Wannier90 file generated by CP2K "//timestamp
         WRITE (iunit, "(3I8)") num_bands, num_kpts, nntot
      ELSE
         iunit = -1
      END IF
      ! create a list of unique b vectors and a table of pointers
      ! nblist(ik,i) -> +/- b_latt(1:3,x)
      ALLOCATE (nblist(num_kpts, nntot))
      ALLOCATE (b_latt(3, num_kpts*nntot))
      nblist(:, :) = 0
      nbs = 0
      DO ik = 1, num_kpts
         DO i = 1, nntot
            bvec(1:3) = kpt_latt(1:3, nnlist(ik, i)) - kpt_latt(1:3, ik) + nncell(1:3, ik, i)
            ibs = 0
            DO k = 1, nbs
               IF (SUM(ABS(bvec(1:3) - b_latt(1:3, k))) < 1.e-6_dp) THEN
                  ibs = k
                  EXIT
               END IF
               IF (SUM(ABS(bvec(1:3) + b_latt(1:3, k))) < 1.e-6_dp) THEN
                  ibs = -k
                  EXIT
               END IF
            END DO
            IF (ibs /= 0) THEN
               ! old lattice vector
               nblist(ik, i) = ibs
            ELSE
               ! new lattice vector
               nbs = nbs + 1
               b_latt(1:3, nbs) = bvec(1:3)
               nblist(ik, i) = nbs
            END IF
         END DO
      END DO
      ! calculate all the operator matrices (a|bvec|b)
      ALLOCATE (berry_matrix(nbs))
      DO i = 1, nbs
         NULLIFY (berry_matrix(i)%cosmat)
         NULLIFY (berry_matrix(i)%sinmat)
         bvec(1:3) = twopi*MATMUL(TRANSPOSE(cell%h_inv(1:3, 1:3)), b_latt(1:3, i))
         CALL build_berry_kpoint_matrix(qs_env_kp, berry_matrix(i)%cosmat, &
                                        berry_matrix(i)%sinmat, bvec)
      END DO
      ! work matrices for MOs (all group)
      kp => kpoint%kp_env(1)%kpoint_env
      CALL get_mo_set(kp%mos(1, 1), nmo=nmo)
      NULLIFY (matrix_struct_ao, matrix_struct_work)
      CALL cp_fm_struct_create(matrix_struct_work, nrow_global=nao, &
                               ncol_global=nmo, &
                               para_env=para_env, &
                               context=blacs_env)
      DO i = 1, 2
         CALL cp_fm_create(fmk1(i), matrix_struct_work)
         CALL cp_fm_create(fmk2(i), matrix_struct_work)
      END DO
      CALL cp_cfm_create(fmk1_cfm, matrix_struct_work)
      CALL cp_cfm_create(fmk2_cfm, matrix_struct_work)
      CALL cp_cfm_create(tmp_cfm, matrix_struct_work)
      CALL cp_fm_struct_create(matrix_struct_ao, nrow_global=nao, &
                               ncol_global=nao, &
                               para_env=para_env, &
                               context=blacs_env)
      CALL cp_fm_create(mat_real, matrix_struct_ao)
      CALL cp_fm_create(mat_imag, matrix_struct_ao)
      CALL cp_cfm_create(omat_cfm, matrix_struct_ao)
      ! work matrices for Mmn(k,b) integrals
      NULLIFY (matrix_struct_mmn)
      CALL cp_fm_struct_create(matrix_struct_mmn, nrow_global=nmo, &
                               ncol_global=nmo, &
                               para_env=para_env, &
                               context=blacs_env)
      CALL cp_fm_create(mmn_real, matrix_struct_mmn)
      CALL cp_fm_create(mmn_imag, matrix_struct_mmn)
      CALL cp_cfm_create(mmn_cfm, matrix_struct_mmn)
      ! allocate some work matrices
      ALLOCATE (rmatrix, cmatrix, rmatrix_full, cmatrix_full)
      CALL dbcsr_create(rmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(cmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL dbcsr_create(rmatrix_full, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(cmatrix_full, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(cmatrix, sab_nl)
      !
      CALL get_kpoint_info(kpoint=kpoint, cell_to_index=cell_to_index)
      NULLIFY (fmdummy)
      nspins = dft_control%nspins
      DO ispin = 1, nspins
         ! loop over all k-points
         DO ik = 1, num_kpts
            ! get the MO coefficients for this k-point
            my_kpgrp = (ik >= kpoint%kp_range(1) .AND. ik <= kpoint%kp_range(2))
            IF (my_kpgrp) THEN
               ikk = ik - kpoint%kp_range(1) + 1
               kp => kpoint%kp_env(ikk)%kpoint_env
               CPASSERT(SIZE(kp%mos, 1) == 2)
               fmr => kp%mos(1, ispin)%mo_coeff
               fmi => kp%mos(2, ispin)%mo_coeff
               CALL cp_fm_copy_general(fmr, fmk1(1), para_env)
               CALL cp_fm_copy_general(fmi, fmk1(2), para_env)
            ELSE
               NULLIFY (fmr, fmi, kp)
               CALL cp_fm_copy_general(fmdummy, fmk1(1), para_env)
               CALL cp_fm_copy_general(fmdummy, fmk1(2), para_env)
            END IF
            CALL cp_fm_to_cfm(fmk1(1), fmk1(2), fmk1_cfm)
            ! loop over all connected neighbors
            DO i = 1, nntot
               ! get the MO coefficients for the connected k-point
               ik2 = nnlist(ik, i)
               mygrp = (ik2 >= kpoint%kp_range(1) .AND. ik2 <= kpoint%kp_range(2))
               IF (mygrp) THEN
                  ikk = ik2 - kpoint%kp_range(1) + 1
                  kp => kpoint%kp_env(ikk)%kpoint_env
                  CPASSERT(SIZE(kp%mos, 1) == 2)
                  fmr => kp%mos(1, ispin)%mo_coeff
                  fmi => kp%mos(2, ispin)%mo_coeff
                  CALL cp_fm_copy_general(fmr, fmk2(1), para_env)
                  CALL cp_fm_copy_general(fmi, fmk2(2), para_env)
               ELSE
                  NULLIFY (fmr, fmi, kp)
                  CALL cp_fm_copy_general(fmdummy, fmk2(1), para_env)
                  CALL cp_fm_copy_general(fmdummy, fmk2(2), para_env)
               END IF
               CALL cp_fm_to_cfm(fmk2(1), fmk2(2), fmk2_cfm)
               !
               ! transfer realspace overlaps to connected k-point
               ibs = nblist(ik, i)
               ksign = SIGN(1.0_dp, REAL(ibs, KIND=dp))
               ibs = ABS(ibs)
               CALL dbcsr_set(rmatrix, 0.0_dp)
               CALL dbcsr_set(cmatrix, 0.0_dp)
               CALL rskp_transform(rmatrix, cmatrix, rsmat=berry_matrix(ibs)%cosmat, ispin=1, &
                                   xkp=kpoint%xkp(1:3, ik2), cell_to_index=cell_to_index, sab_nl=sab_nl, &
                                   is_complex=.FALSE., rs_sign=ksign)
               CALL rskp_transform(cmatrix, rmatrix, rsmat=berry_matrix(ibs)%sinmat, ispin=1, &
                                   xkp=kpoint%xkp(1:3, ik2), cell_to_index=cell_to_index, sab_nl=sab_nl, &
                                   is_complex=.TRUE., rs_sign=ksign)
               !
               ! calculate M_(mn)^(k,b) = C(k)^H O(k,b) C(k+b)
               CALL dbcsr_desymmetrize(rmatrix, rmatrix_full)
               CALL dbcsr_desymmetrize(cmatrix, cmatrix_full)
               CALL copy_dbcsr_to_fm(rmatrix_full, mat_real)
               CALL copy_dbcsr_to_fm(cmatrix_full, mat_imag)
               CALL cp_fm_to_cfm(mat_real, mat_imag, omat_cfm)
               CALL cp_cfm_gemm("N", "N", nao, nmo, nao, CMPLX(1.0_dp, 0.0_dp, KIND=dp), &
                                omat_cfm, fmk2_cfm, CMPLX(0.0_dp, 0.0_dp, KIND=dp), tmp_cfm)
               CALL cp_cfm_gemm("C", "N", nmo, nmo, nao, CMPLX(1.0_dp, 0.0_dp, KIND=dp), &
                                fmk1_cfm, tmp_cfm, CMPLX(0.0_dp, 0.0_dp, KIND=dp), mmn_cfm)
               CALL cp_cfm_to_fm(mmn_cfm, mmn_real, mmn_imag)
               !
               ! write to output file
               IF (reused_scf_mos .AND. use_bloch_phases) THEN
                  ! Reused SCF MOs need the same global Bloch gauge in every overlap block.
                  gauge_arg = twopi*DOT_PRODUCT(kpoint%xkp(1:3, ik2) - kpoint%xkp(1:3, ik), &
                                                phase_center(1:3))
                  gauge_real = COS(gauge_arg)
                  gauge_imag = SIN(gauge_arg)
               ELSE
                  gauge_real = 1.0_dp
                  gauge_imag = 0.0_dp
               END IF
               IF (para_env%is_source()) THEN
                  WRITE (iunit, "(2I8,3I5)") ik, ik2, nncell(1:3, ik, i)
               END IF
               DO ib2 = 1, nmo
                  DO ib1 = 1, nmo
                     CALL cp_fm_get_element(mmn_real, ib1, ib2, rmmn)
                     CALL cp_fm_get_element(mmn_imag, ib1, ib2, cmmn)
                     gauge_tmp = gauge_real*rmmn - gauge_imag*cmmn
                     cmmn = gauge_imag*rmmn + gauge_real*cmmn
                     rmmn = gauge_tmp
                     IF (para_env%is_source()) THEN
                        WRITE (iunit, "(2E30.14)") rmmn, cmmn
                     END IF
                  END DO
               END DO
               !
            END DO
         END DO
      END DO
      DO i = 1, nbs
         CALL dbcsr_deallocate_matrix_set(berry_matrix(i)%cosmat)
         CALL dbcsr_deallocate_matrix_set(berry_matrix(i)%sinmat)
      END DO
      DEALLOCATE (berry_matrix)
      CALL cp_fm_struct_release(matrix_struct_work)
      DO i = 1, 2
         CALL cp_fm_release(fmk1(i))
         CALL cp_fm_release(fmk2(i))
      END DO
      CALL cp_cfm_release(fmk1_cfm)
      CALL cp_cfm_release(fmk2_cfm)
      CALL cp_cfm_release(tmp_cfm)
      CALL cp_fm_struct_release(matrix_struct_ao)
      CALL cp_fm_release(mat_real)
      CALL cp_fm_release(mat_imag)
      CALL cp_cfm_release(omat_cfm)
      CALL cp_fm_struct_release(matrix_struct_mmn)
      CALL cp_fm_release(mmn_real)
      CALL cp_fm_release(mmn_imag)
      CALL cp_cfm_release(mmn_cfm)
      CALL dbcsr_deallocate_matrix(rmatrix)
      CALL dbcsr_deallocate_matrix(cmatrix)
      CALL dbcsr_deallocate_matrix(rmatrix_full)
      CALL dbcsr_deallocate_matrix(cmatrix_full)
      !
      IF (para_env%is_source()) THEN
         CALL close_file(iunit)
      END IF
      !
      ! Calculate and print Projections
      !
      ! Print eigenvalues
      nspins = dft_control%nspins
      kp => kpoint%kp_env(1)%kpoint_env
      CALL get_mo_set(kp%mos(1, 1), nmo=nmo)
      ALLOCATE (eigval(nmo))
      CALL get_kpoint_info(kpoint, nkp=nkp, kp_range=kp_range, xkp=xkp)
      IF (para_env%is_source()) THEN
         WRITE (filename, '(A,A)') TRIM(seed_name), ".eig"
         CALL open_file(filename, unit_number=iunit, file_status="UNKNOWN", file_action="WRITE")
      ELSE
         iunit = -1
      END IF
      !
      DO ik = 1, nkp
         my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
         DO ispin = 1, nspins
            IF (my_kpgrp) THEN
               ikpgr = ik - kp_range(1) + 1
               kp => kpoint%kp_env(ikpgr)%kpoint_env
               CALL get_mo_set(kp%mos(1, ispin), eigenvalues=eigenvalues)
               eigval(1:nmo) = eigenvalues(1:nmo)
            ELSE
               eigval(1:nmo) = 0.0_dp
            END IF
            CALL kpoint%para_env_inter_kp%sum(eigval)
            eigval(1:nmo) = eigval(1:nmo)*evolt
            ! output
            IF (iunit > 0) THEN
               DO ib = 1, nmo
                  WRITE (iunit, "(2I8,F24.14)") ib, ik, eigval(ib)
               END DO
            END IF
         END DO
      END DO
      IF (para_env%is_source()) THEN
         CALL close_file(iunit)
      END IF
      !
      ! clean up
      DEALLOCATE (kpt_latt, atoms_cart, atom_symbols, eigval)
      DEALLOCATE (nnlist, nncell)
      DEALLOCATE (nblist, b_latt)
      IF (nexcl > 0) THEN
         DEALLOCATE (exclude_bands)
      END IF
      IF (do_kpoints) THEN
         NULLIFY (qs_env_kp)
      ELSE
         CALL qs_env_release(qs_env_kp)
         DEALLOCATE (qs_env_kp)
         NULLIFY (qs_env_kp)
      END IF

      CALL kpoint_release(kpoint)

   END SUBROUTINE wannier90_files

! **************************************************************************************************
!> \brief Reconstruct a full Wannier90 k-point MO set from the SCF k-point MOs.
!> \param kpoint full Wannier90 export k-point object
!> \param qs_kpoint SCF k-point object
!> \param matrix_s real-space overlap matrix
!> \param matrix_ks real-space Kohn-Sham matrix
!> \param cell_to_index real-space cell index table
!> \param sab_nl overlap neighbor list
!> \param para_env global parallel environment
!> \param success true if all full-mesh MOs were reconstructed
!> \param reason diagnostic message when reconstruction is not possible
!> \param aligned_degenerate_blocks number of aligned degenerate MO blocks
!> \param aligned_degenerate_max_size largest aligned degenerate MO block
!> \param aligned_degenerate_min_svalue smallest S(k)-metric subspace singular value
! **************************************************************************************************
   SUBROUTINE prepare_wannier90_scf_mos(kpoint, qs_kpoint, matrix_s, matrix_ks, cell_to_index, &
                                        sab_nl, para_env, success, reason, aligned_degenerate_blocks, &
                                        aligned_degenerate_max_size, &
                                        aligned_degenerate_min_svalue)
      TYPE(kpoint_type), POINTER                         :: kpoint, qs_kpoint
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s, matrix_ks
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(OUT)                               :: success
      CHARACTER(LEN=*), INTENT(OUT)                      :: reason
      INTEGER, INTENT(OUT)                               :: aligned_degenerate_blocks, &
                                                            aligned_degenerate_max_size
      REAL(KIND=dp), INTENT(OUT)                         :: aligned_degenerate_min_svalue

      CHARACTER(LEN=default_string_length)               :: best_reason, candidate_reason
      INTEGER :: aligned_blocks, aligned_max_size, candidate_aligned_blocks, &
         candidate_aligned_max_size, ik, ikpgr, ikred, ispin, isym_try, min_gap_band, &
         min_gap_kpoint, min_gap_spin, nao, nao_src, nmo, nmo_src, nspins, nsymmetry, &
         num_candidates
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: source_kpoint, sym_index
      INTEGER, DIMENSION(2)                              :: kp_range, source_kp_range
      LOGICAL                                            :: my_kpgrp, my_source_kpgrp, ok, &
                                                            source_window
      REAL(KIND=dp) :: aligned_min_svalue, band_gap, best_residual, candidate_residual, &
         degenerate_band_tol, local_min_band_gap, min_band_gap, source_owner_count, &
         source_window_min_svalue
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues_buffer, occupation_buffer, &
                                                            source_eigenvalues_buffer, &
                                                            source_occupation_buffer
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_source, matrix_struct_work
      TYPE(cp_fm_type)                                   :: dst_imag, dst_imag_full, dst_real, &
                                                            dst_real_full, src_imag, &
                                                            src_imag_full, src_real, src_real_full
      TYPE(cp_fm_type), POINTER                          :: dst_fmi, dst_fmr, src_fmi, src_fmr
      TYPE(kpoint_env_type), POINTER                     :: kp, kp_source
      TYPE(kpoint_sym_type), POINTER                     :: kpsym

      success = .FALSE.
      reason = ""
      aligned_degenerate_blocks = 0
      aligned_degenerate_max_size = 0
      aligned_degenerate_min_svalue = HUGE(1.0_dp)
      NULLIFY (matrix_struct_source, matrix_struct_work, src_fmr, src_fmi, dst_fmr, dst_fmi)

      IF (.NOT. ASSOCIATED(kpoint)) THEN
         reason = "internal Wannier90 k-point object is not available"
         RETURN
      END IF
      IF (.NOT. ASSOCIATED(qs_kpoint)) THEN
         reason = "SCF k-point object is not available"
         RETURN
      END IF
      IF (.NOT. ASSOCIATED(kpoint%kp_env) .OR. .NOT. ASSOCIATED(qs_kpoint%kp_env)) THEN
         reason = "k-point MO environments are not initialized"
         RETURN
      END IF
      IF (.NOT. ASSOCIATED(kpoint%blacs_env)) THEN
         reason = "Wannier90 k-point BLACS environment is not initialized"
         RETURN
      END IF

      CALL build_wannier90_scf_mapping(kpoint, qs_kpoint, source_kpoint, sym_index, ok, reason)
      IF (.NOT. ok) RETURN
      nsymmetry = COUNT(sym_index > 0)

      kp => kpoint%kp_env(1)%kpoint_env
      nspins = SIZE(kp%mos, 2)
      IF (SIZE(kp%mos, 1) < 2) THEN
         reason = "Wannier90 export k-point MOs are not complex-valued"
         DEALLOCATE (source_kpoint, sym_index)
         RETURN
      END IF
      CALL get_mo_set(kp%mos(1, 1), nao=nao, nmo=nmo)

      kp_source => qs_kpoint%kp_env(1)%kpoint_env
      IF (SIZE(kp_source%mos, 1) < 2) THEN
         reason = "SCF MOs are real-valued; complex symmetry phases cannot be reconstructed"
         DEALLOCATE (source_kpoint, sym_index)
         RETURN
      END IF
      CALL get_mo_set(kp_source%mos(1, 1), nao=nao_src, nmo=nmo_src)
      CALL para_env%max(nao_src)
      CALL para_env%max(nmo_src)
      IF (nao_src /= nao) THEN
         reason = "SCF and Wannier90 MO bases have different AO dimensions"
         DEALLOCATE (source_kpoint, sym_index)
         RETURN
      END IF
      IF (nmo_src < nmo) THEN
         reason = "SCF MO set has fewer bands than the Wannier90 export"
         DEALLOCATE (source_kpoint, sym_index)
         RETURN
      END IF
      source_window = nmo_src > nmo
      degenerate_band_tol = 1.0e-8_dp
      CALL get_kpoint_info(qs_kpoint, kp_range=source_kp_range)
      IF (source_kp_range(1) /= 1 .OR. source_kp_range(2) /= qs_kpoint%nkp) THEN
         reason = "SCF k-point symmetry data are distributed over k-point parallel groups"
         DEALLOCATE (source_kpoint, sym_index)
         RETURN
      END IF
      ! Positive symmetry entries require atom/AO rotations and Bloch phases. Degenerate subspaces
      ! fully contained in the exported band window are aligned below; only guard when the Wannier90
      ! window cuts through a degenerate SCF manifold at the upper band edge.
      IF (nsymmetry > 0 .AND. nmo_src > nmo) THEN
         local_min_band_gap = HUGE(1.0_dp)
         min_gap_band = nmo
         min_gap_kpoint = 0
         min_gap_spin = 0
         DO ikred = source_kp_range(1), source_kp_range(2)
            ikpgr = ikred - source_kp_range(1) + 1
            kp_source => qs_kpoint%kp_env(ikpgr)%kpoint_env
            DO ispin = 1, nspins
               CALL get_mo_set(kp_source%mos(1, ispin), eigenvalues=eigenvalues)
               band_gap = ABS(eigenvalues(nmo + 1) - eigenvalues(nmo))
               IF (band_gap < local_min_band_gap) THEN
                  local_min_band_gap = band_gap
                  min_gap_band = nmo
                  min_gap_kpoint = ikred
                  min_gap_spin = ispin
               END IF
            END DO
         END DO
         min_band_gap = local_min_band_gap
         CALL para_env%min(min_band_gap)
         IF (ABS(local_min_band_gap - min_band_gap) > degenerate_band_tol*EPSILON(1.0_dp)) THEN
            min_gap_kpoint = 0
            min_gap_spin = 0
         END IF
         CALL para_env%max(min_gap_kpoint)
         CALL para_env%max(min_gap_spin)
         CALL para_env%max(min_gap_band)
         IF (min_band_gap < degenerate_band_tol) THEN
            WRITE (reason, "(A,ES9.2,A,I0,A,I0,A,I0)") &
               "degenerate atom/AO W90 reuse guarded: edge gap=", min_band_gap, ", k=", &
               min_gap_kpoint, ", s=", min_gap_spin, ", nband=", min_gap_band
            DEALLOCATE (source_kpoint, sym_index)
            RETURN
         END IF
      END IF
      blacs_env => kpoint%blacs_env
      CALL cp_fm_struct_create(matrix_struct_work, nrow_global=nao, ncol_global=nmo, &
                               para_env=para_env, context=blacs_env)
      CALL cp_fm_create(src_real, matrix_struct_work)
      CALL cp_fm_create(src_imag, matrix_struct_work)
      CALL cp_fm_create(dst_real, matrix_struct_work)
      CALL cp_fm_create(dst_imag, matrix_struct_work)
      IF (source_window) THEN
         CALL cp_fm_struct_create(matrix_struct_source, nrow_global=nao, ncol_global=nmo_src, &
                                  para_env=para_env, context=blacs_env)
         CALL cp_fm_create(src_real_full, matrix_struct_source)
         CALL cp_fm_create(src_imag_full, matrix_struct_source)
         CALL cp_fm_create(dst_real_full, matrix_struct_source)
         CALL cp_fm_create(dst_imag_full, matrix_struct_source)
      END IF
      ALLOCATE (eigenvalues_buffer(nmo), occupation_buffer(nmo), &
                source_eigenvalues_buffer(nmo_src), source_occupation_buffer(nmo_src))

      CALL get_kpoint_info(kpoint, kp_range=kp_range)
      CALL get_kpoint_info(qs_kpoint, kp_range=source_kp_range)

      DO ik = 1, kpoint%nkp
         ikred = source_kpoint(ik)
         my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
         my_source_kpgrp = (ikred >= source_kp_range(1) .AND. ikred <= source_kp_range(2))
         DO ispin = 1, nspins
            source_eigenvalues_buffer(1:nmo_src) = 0.0_dp
            source_occupation_buffer(1:nmo_src) = 0.0_dp
            IF (my_source_kpgrp) THEN
               ikpgr = ikred - source_kp_range(1) + 1
               kp_source => qs_kpoint%kp_env(ikpgr)%kpoint_env
               src_fmr => kp_source%mos(1, ispin)%mo_coeff
               src_fmi => kp_source%mos(2, ispin)%mo_coeff
               CALL get_mo_set(kp_source%mos(1, ispin), eigenvalues=eigenvalues, &
                               occupation_numbers=occupation)
               source_eigenvalues_buffer(1:nmo_src) = eigenvalues(1:nmo_src)
               source_occupation_buffer(1:nmo_src) = occupation(1:nmo_src)
            ELSE
               NULLIFY (src_fmr, src_fmi)
            END IF
            IF (my_source_kpgrp) THEN
               source_owner_count = 1.0_dp
            ELSE
               source_owner_count = 0.0_dp
            END IF
            CALL para_env%sum(source_owner_count)
            CALL para_env%sum(source_eigenvalues_buffer)
            CALL para_env%sum(source_occupation_buffer)
            IF (source_owner_count > 0.0_dp) THEN
               source_eigenvalues_buffer(1:nmo_src) = &
                  source_eigenvalues_buffer(1:nmo_src)/source_owner_count
               source_occupation_buffer(1:nmo_src) = &
                  source_occupation_buffer(1:nmo_src)/source_owner_count
            END IF
            eigenvalues_buffer(1:nmo) = source_eigenvalues_buffer(1:nmo)
            occupation_buffer(1:nmo) = source_occupation_buffer(1:nmo)
            IF (source_window) THEN
               CALL cp_fm_copy_general(src_fmr, src_real_full, para_env)
               CALL cp_fm_copy_general(src_fmi, src_imag_full, para_env)
               CALL copy_wannier90_mo_window(src_real_full, src_real, nmo)
               CALL copy_wannier90_mo_window(src_imag_full, src_imag, nmo)
            ELSE
               CALL cp_fm_copy_general(src_fmr, src_real, para_env)
               CALL cp_fm_copy_general(src_fmi, src_imag, para_env)
            END IF

            ok = .FALSE.
            reason = ""
            aligned_blocks = 0
            aligned_max_size = 0
            aligned_min_svalue = 0.0_dp
            IF (sym_index(ik) > 0) THEN
               kpsym => qs_kpoint%kp_sym(ikred)%kpoint_sym
               IF (ASSOCIATED(kpsym)) THEN
                  best_reason = ""
                  best_residual = HUGE(1.0_dp)
                  num_candidates = 0
                  ! Little-group operations can reach the same target k-point; keep the first valid eigenspace.
                  DO isym_try = 1, kpsym%nwred
                     IF (.NOT. same_periodic_kpoint(kpoint%xkp(1:3, ik), &
                                                    kpsym%xkp(1:3, isym_try))) CYCLE
                     num_candidates = num_candidates + 1
                     CALL transform_wannier90_scf_mo(src_real, src_imag, dst_real, dst_imag, &
                                                     qs_kpoint, ikred, isym_try, para_env, ok, &
                                                     candidate_reason)
                     IF (.NOT. ok) THEN
                        reason = candidate_reason
                        CYCLE
                     END IF
                     CALL ritz_stabilize_wannier90_subspace(dst_real, dst_imag, matrix_s, matrix_ks, &
                                                            kpoint%xkp(1:3, ik), cell_to_index, &
                                                            sab_nl, ispin, eigenvalues_buffer, &
                                                            degenerate_band_tol, ok, candidate_reason, &
                                                            candidate_aligned_blocks, &
                                                            candidate_aligned_max_size, &
                                                            aligned_min_svalue, candidate_residual)
                     IF (candidate_residual < best_residual) THEN
                        best_residual = candidate_residual
                        best_reason = candidate_reason
                     END IF
                     IF (.NOT. ok) THEN
                        IF (source_window) THEN
                           CALL transform_wannier90_scf_mo(src_real_full, src_imag_full, &
                                                           dst_real_full, dst_imag_full, qs_kpoint, &
                                                           ikred, isym_try, para_env, ok, candidate_reason)
                           IF (ok) THEN
                              CALL ritz_reconstruct_wannier90_window(dst_real_full, dst_imag_full, &
                                                                     dst_real, dst_imag, matrix_s, &
                                                                     matrix_ks, kpoint%xkp(1:3, ik), &
                                                                     cell_to_index, sab_nl, ispin, &
                                                                     eigenvalues_buffer, nmo, ok, &
                                                                     candidate_reason, source_window_min_svalue, &
                                                                     candidate_residual)
                              IF (candidate_residual < best_residual) THEN
                                 best_residual = candidate_residual
                                 best_reason = candidate_reason
                              END IF
                           END IF
                        ELSE
                           CALL ritz_reconstruct_wannier90_window(dst_real, dst_imag, dst_real, &
                                                                  dst_imag, matrix_s, matrix_ks, &
                                                                  kpoint%xkp(1:3, ik), cell_to_index, &
                                                                  sab_nl, ispin, eigenvalues_buffer, nmo, &
                                                                  ok, candidate_reason, source_window_min_svalue, &
                                                                  candidate_residual)
                           IF (candidate_residual < best_residual) THEN
                              best_residual = candidate_residual
                              best_reason = candidate_reason
                           END IF
                        END IF
                     END IF
                     IF (ok) THEN
                        aligned_blocks = candidate_aligned_blocks
                        aligned_max_size = candidate_aligned_max_size
                        sym_index(ik) = isym_try
                        EXIT
                     END IF
                     reason = candidate_reason
                  END DO
                  IF (.NOT. ok .AND. num_candidates > 0 .AND. best_residual < HUGE(1.0_dp)) THEN
                     WRITE (reason, "(A,I0,A,ES9.2,A,I0,A,A32)") &
                        "atom/AO W90 guarded: best/", num_candidates, "=", best_residual, &
                        " k=", ik, " ", TRIM(best_reason)
                  ELSEIF (.NOT. ok .AND. num_candidates == 0) THEN
                     reason = "no matching SCF symmetry operation candidate"
                  END IF
               ELSE
                  reason = "SCF k-point symmetry operation is not available"
               END IF
            ELSE
               CALL transform_wannier90_scf_mo(src_real, src_imag, dst_real, dst_imag, qs_kpoint, &
                                               ikred, sym_index(ik), para_env, ok, reason)
            END IF
            IF (ok .AND. sym_index(ik) <= 0) THEN
               ! Even a direct k-point copy must be a closed H(k),S(k) subspace. This catches
               ! incomplete degenerate band windows before they can be exported to Wannier90.
               CALL ritz_stabilize_wannier90_subspace(dst_real, dst_imag, matrix_s, matrix_ks, &
                                                      kpoint%xkp(1:3, ik), cell_to_index, sab_nl, &
                                                      ispin, eigenvalues_buffer, degenerate_band_tol, &
                                                      ok, reason, aligned_blocks, aligned_max_size, &
                                                      aligned_min_svalue, candidate_residual)
               IF (.NOT. ok .AND. source_window) THEN
                  CALL transform_wannier90_scf_mo(src_real_full, src_imag_full, dst_real_full, &
                                                  dst_imag_full, qs_kpoint, ikred, sym_index(ik), &
                                                  para_env, ok, reason)
                  IF (ok) THEN
                     CALL ritz_reconstruct_wannier90_window(dst_real_full, dst_imag_full, dst_real, &
                                                            dst_imag, matrix_s, matrix_ks, &
                                                            kpoint%xkp(1:3, ik), cell_to_index, &
                                                            sab_nl, ispin, eigenvalues_buffer, nmo, ok, &
                                                            reason, source_window_min_svalue, &
                                                            candidate_residual)
                  END IF
               ELSEIF (.NOT. ok) THEN
                  CALL ritz_reconstruct_wannier90_window(dst_real, dst_imag, dst_real, dst_imag, &
                                                         matrix_s, matrix_ks, kpoint%xkp(1:3, ik), &
                                                         cell_to_index, sab_nl, ispin, &
                                                         eigenvalues_buffer, nmo, ok, reason, &
                                                         source_window_min_svalue, candidate_residual)
               END IF
            END IF
            IF (.NOT. ok) THEN
               CALL cp_fm_release(src_real)
               CALL cp_fm_release(src_imag)
               CALL cp_fm_release(dst_real)
               CALL cp_fm_release(dst_imag)
               CALL cp_fm_struct_release(matrix_struct_work)
               IF (source_window) THEN
                  CALL cp_fm_release(src_real_full)
                  CALL cp_fm_release(src_imag_full)
                  CALL cp_fm_release(dst_real_full)
                  CALL cp_fm_release(dst_imag_full)
                  CALL cp_fm_struct_release(matrix_struct_source)
               END IF
               DEALLOCATE (source_kpoint, sym_index, eigenvalues_buffer, occupation_buffer, &
                           source_eigenvalues_buffer, source_occupation_buffer)
               RETURN
            END IF
            IF (sym_index(ik) /= 0) THEN
               aligned_degenerate_blocks = aligned_degenerate_blocks + aligned_blocks
               aligned_degenerate_max_size = MAX(aligned_degenerate_max_size, aligned_max_size)
               IF (aligned_blocks > 0) THEN
                  aligned_degenerate_min_svalue = MIN(aligned_degenerate_min_svalue, aligned_min_svalue)
               END IF
            END IF

            IF (my_kpgrp) THEN
               ikpgr = ik - kp_range(1) + 1
               kp => kpoint%kp_env(ikpgr)%kpoint_env
               dst_fmr => kp%mos(1, ispin)%mo_coeff
               dst_fmi => kp%mos(2, ispin)%mo_coeff
               CALL get_mo_set(kp%mos(1, ispin), eigenvalues=eigenvalues, &
                               occupation_numbers=occupation)
               eigenvalues(1:nmo) = eigenvalues_buffer(1:nmo)
               occupation(1:nmo) = occupation_buffer(1:nmo)
               CALL get_mo_set(kp%mos(2, ispin), eigenvalues=eigenvalues, &
                               occupation_numbers=occupation)
               IF (ASSOCIATED(eigenvalues)) eigenvalues(1:nmo) = eigenvalues_buffer(1:nmo)
               IF (ASSOCIATED(occupation)) occupation(1:nmo) = occupation_buffer(1:nmo)
            ELSE
               NULLIFY (dst_fmr, dst_fmi)
            END IF
            CALL cp_fm_copy_general(dst_real, dst_fmr, para_env)
            CALL cp_fm_copy_general(dst_imag, dst_fmi, para_env)
         END DO
      END DO

      CALL cp_fm_release(src_real)
      CALL cp_fm_release(src_imag)
      CALL cp_fm_release(dst_real)
      CALL cp_fm_release(dst_imag)
      CALL cp_fm_struct_release(matrix_struct_work)
      IF (source_window) THEN
         CALL cp_fm_release(src_real_full)
         CALL cp_fm_release(src_imag_full)
         CALL cp_fm_release(dst_real_full)
         CALL cp_fm_release(dst_imag_full)
         CALL cp_fm_struct_release(matrix_struct_source)
      END IF
      DEALLOCATE (source_kpoint, sym_index, eigenvalues_buffer, occupation_buffer, &
                  source_eigenvalues_buffer, source_occupation_buffer)
      IF (aligned_degenerate_blocks == 0) aligned_degenerate_min_svalue = 0.0_dp
      success = .TRUE.

   END SUBROUTINE prepare_wannier90_scf_mos

! **************************************************************************************************
!> \brief Save a full-mesh Wannier90 MO reference on all ranks for diagnostic validation.
!> \param kpoint full Wannier90 export k-point object
!> \param nspins number of spin channels
!> \param para_env global parallel environment
!> \param mo_real real MO coefficients, indexed as AO, MO, k-point, spin
!> \param mo_imag imaginary MO coefficients, indexed as AO, MO, k-point, spin
!> \param eigenvalue_snapshot MO eigenvalues, indexed as MO, k-point, spin
! **************************************************************************************************
   SUBROUTINE save_wannier90_mo_snapshot(kpoint, nspins, para_env, mo_real, mo_imag, &
                                         eigenvalue_snapshot)
      TYPE(kpoint_type), POINTER                         :: kpoint
      INTEGER, INTENT(IN)                                :: nspins
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :), INTENT(OUT)              :: mo_real, mo_imag
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: eigenvalue_snapshot

      INTEGER                                            :: ik, ikpgr, ispin, nao, nkp, nmo
      INTEGER, DIMENSION(2)                              :: kp_range
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: owner_weight
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(cp_fm_type), POINTER                          :: fmi, fmr
      TYPE(kpoint_env_type), POINTER                     :: kp

      CALL get_kpoint_info(kpoint, nkp=nkp, kp_range=kp_range)
      kp => kpoint%kp_env(1)%kpoint_env
      CALL get_mo_set(kp%mos(1, 1), nao=nao, nmo=nmo)
      ALLOCATE (mo_real(nao, nmo, nkp, nspins), mo_imag(nao, nmo, nkp, nspins), &
                eigenvalue_snapshot(nmo, nkp, nspins), owner_weight(nkp, nspins))
      mo_real(:, :, :, :) = 0.0_dp
      mo_imag(:, :, :, :) = 0.0_dp
      eigenvalue_snapshot(:, :, :) = 0.0_dp
      owner_weight(:, :) = 0.0_dp
      DO ik = kp_range(1), kp_range(2)
         ikpgr = ik - kp_range(1) + 1
         kp => kpoint%kp_env(ikpgr)%kpoint_env
         DO ispin = 1, nspins
            fmr => kp%mos(1, ispin)%mo_coeff
            fmi => kp%mos(2, ispin)%mo_coeff
            CALL cp_fm_get_submatrix(fmr, mo_real(:, :, ik, ispin))
            CALL cp_fm_get_submatrix(fmi, mo_imag(:, :, ik, ispin))
            CALL get_mo_set(kp%mos(1, ispin), eigenvalues=eigenvalues)
            eigenvalue_snapshot(1:nmo, ik, ispin) = eigenvalues(1:nmo)
            owner_weight(ik, ispin) = 1.0_dp
         END DO
      END DO
      CALL para_env%sum(mo_real)
      CALL para_env%sum(mo_imag)
      CALL para_env%sum(eigenvalue_snapshot)
      CALL para_env%sum(owner_weight)
      DO ik = 1, nkp
         DO ispin = 1, nspins
            IF (owner_weight(ik, ispin) > 0.0_dp) THEN
               mo_real(:, :, ik, ispin) = mo_real(:, :, ik, ispin)/owner_weight(ik, ispin)
               mo_imag(:, :, ik, ispin) = mo_imag(:, :, ik, ispin)/owner_weight(ik, ispin)
               eigenvalue_snapshot(:, ik, ispin) = &
                  eigenvalue_snapshot(:, ik, ispin)/owner_weight(ik, ispin)
            END IF
         END DO
      END DO
      DEALLOCATE (owner_weight)

   END SUBROUTINE save_wannier90_mo_snapshot

! **************************************************************************************************
!> \brief Restore a full-mesh Wannier90 MO reference after a failed diagnostic reuse attempt.
!> \param kpoint full Wannier90 export k-point object
!> \param mo_real real MO coefficient snapshot
!> \param mo_imag imaginary MO coefficient snapshot
!> \param eigenvalue_snapshot MO eigenvalue snapshot
! **************************************************************************************************
   SUBROUTINE restore_wannier90_mo_snapshot(kpoint, mo_real, mo_imag, eigenvalue_snapshot)
      TYPE(kpoint_type), POINTER                         :: kpoint
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: mo_real, mo_imag
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: eigenvalue_snapshot

      INTEGER                                            :: ik, ikpgr, ispin, nmo, nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(cp_fm_type), POINTER                          :: fmi, fmr
      TYPE(kpoint_env_type), POINTER                     :: kp

      CALL get_kpoint_info(kpoint, kp_range=kp_range)
      nmo = SIZE(eigenvalue_snapshot, 1)
      nspins = SIZE(eigenvalue_snapshot, 3)
      DO ik = kp_range(1), kp_range(2)
         ikpgr = ik - kp_range(1) + 1
         kp => kpoint%kp_env(ikpgr)%kpoint_env
         DO ispin = 1, nspins
            fmr => kp%mos(1, ispin)%mo_coeff
            fmi => kp%mos(2, ispin)%mo_coeff
            CALL cp_fm_set_submatrix(fmr, mo_real(:, :, ik, ispin))
            CALL cp_fm_set_submatrix(fmi, mo_imag(:, :, ik, ispin))
            CALL get_mo_set(kp%mos(1, ispin), eigenvalues=eigenvalues)
            eigenvalues(1:nmo) = eigenvalue_snapshot(1:nmo, ik, ispin)
            CALL get_mo_set(kp%mos(2, ispin), eigenvalues=eigenvalues)
            IF (ASSOCIATED(eigenvalues)) eigenvalues(1:nmo) = eigenvalue_snapshot(1:nmo, ik, ispin)
         END DO
      END DO

   END SUBROUTINE restore_wannier90_mo_snapshot

! **************************************************************************************************
!> \brief Validate current Wannier90 MOs against a saved full-mesh diagonalization reference.
!> \param kpoint full Wannier90 export k-point object
!> \param matrix_s real-space overlap matrix
!> \param cell_to_index real-space cell index table
!> \param sab_nl overlap neighbor list
!> \param para_env global parallel environment
!> \param reference_mo_real real MO coefficient reference
!> \param reference_mo_imag imaginary MO coefficient reference
!> \param reference_eigenvalues MO eigenvalue reference
!> \param success true if the reconstructed MOs match the reference subspaces
!> \param max_subspace_deviation largest deviation of S(k)-metric singular values from one
!> \param min_svalue smallest S(k)-metric singular value
!> \param max_eigenvalue_deviation largest eigenvalue deviation
! **************************************************************************************************
   SUBROUTINE validate_wannier90_reused_mos(kpoint, matrix_s, cell_to_index, sab_nl, para_env, &
                                            reference_mo_real, reference_mo_imag, reference_eigenvalues, &
                                            success, max_subspace_deviation, min_svalue, &
                                            max_eigenvalue_deviation)
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: reference_mo_real, reference_mo_imag
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: reference_eigenvalues
      LOGICAL, INTENT(OUT)                               :: success
      REAL(KIND=dp), INTENT(OUT)                         :: max_subspace_deviation, min_svalue, &
                                                            max_eigenvalue_deviation

      REAL(KIND=dp), PARAMETER                           :: eigenvalue_tol = 1.0e-8_dp, &
                                                            subspace_tol = 1.0e-4_dp

      INTEGER                                            :: ik, ikpgr, ispin, nao, nkp, nmo, nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      LOGICAL                                            :: my_kpgrp, ok
      REAL(KIND=dp)                                      :: candidate_deviation, candidate_svalue, &
                                                            owner_count
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalue_buffer
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_work
      TYPE(cp_fm_type)                                   :: cand_imag, cand_real, ref_imag, ref_real
      TYPE(cp_fm_type), POINTER                          :: fmi, fmr
      TYPE(kpoint_env_type), POINTER                     :: kp

      success = .FALSE.
      max_subspace_deviation = 0.0_dp
      min_svalue = HUGE(1.0_dp)
      max_eigenvalue_deviation = 0.0_dp
      NULLIFY (matrix_struct_work, fmr, fmi)

      CALL get_kpoint_info(kpoint, nkp=nkp, kp_range=kp_range)
      kp => kpoint%kp_env(1)%kpoint_env
      nspins = SIZE(kp%mos, 2)
      CALL get_mo_set(kp%mos(1, 1), nao=nao, nmo=nmo)
      CPASSERT(SIZE(reference_mo_real, 1) == nao)
      CPASSERT(SIZE(reference_mo_real, 2) == nmo)
      CPASSERT(SIZE(reference_mo_real, 3) == nkp)
      CPASSERT(SIZE(reference_mo_real, 4) == nspins)

      CALL cp_fm_get_info(kp%mos(1, 1)%mo_coeff, matrix_struct=matrix_struct_work)
      CALL cp_fm_create(ref_real, matrix_struct_work)
      CALL cp_fm_create(ref_imag, matrix_struct_work)
      CALL cp_fm_create(cand_real, matrix_struct_work)
      CALL cp_fm_create(cand_imag, matrix_struct_work)
      ALLOCATE (eigenvalue_buffer(nmo))

      DO ispin = 1, nspins
         DO ik = 1, nkp
            CALL cp_fm_set_submatrix(ref_real, reference_mo_real(:, :, ik, ispin))
            CALL cp_fm_set_submatrix(ref_imag, reference_mo_imag(:, :, ik, ispin))
            my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
            IF (my_kpgrp) THEN
               ikpgr = ik - kp_range(1) + 1
               kp => kpoint%kp_env(ikpgr)%kpoint_env
               fmr => kp%mos(1, ispin)%mo_coeff
               fmi => kp%mos(2, ispin)%mo_coeff
               CALL get_mo_set(kp%mos(1, ispin), eigenvalues=eigenvalues)
               eigenvalue_buffer(1:nmo) = eigenvalues(1:nmo)
            ELSE
               NULLIFY (fmr, fmi)
               eigenvalue_buffer(1:nmo) = 0.0_dp
            END IF
            CALL cp_fm_copy_general(fmr, cand_real, para_env)
            CALL cp_fm_copy_general(fmi, cand_imag, para_env)
            IF (my_kpgrp) THEN
               owner_count = 1.0_dp
            ELSE
               owner_count = 0.0_dp
            END IF
            CALL para_env%sum(owner_count)
            CALL para_env%sum(eigenvalue_buffer)
            IF (owner_count > 0.0_dp) eigenvalue_buffer(1:nmo) = eigenvalue_buffer(1:nmo)/owner_count
            max_eigenvalue_deviation = MAX(max_eigenvalue_deviation, &
                                           MAXVAL(ABS(eigenvalue_buffer(1:nmo) - &
                                                      reference_eigenvalues(1:nmo, ik, ispin))))
            CALL measure_wannier90_subspace_error(ref_real, ref_imag, cand_real, cand_imag, matrix_s, &
                                                  kpoint%xkp(1:3, ik), cell_to_index, sab_nl, para_env, &
                                                  ok, candidate_deviation, candidate_svalue)
            IF (.NOT. ok) THEN
               max_subspace_deviation = HUGE(1.0_dp)
            ELSE
               max_subspace_deviation = MAX(max_subspace_deviation, candidate_deviation)
               min_svalue = MIN(min_svalue, candidate_svalue)
            END IF
         END DO
      END DO
      CALL para_env%max(max_subspace_deviation)
      CALL para_env%min(min_svalue)
      CALL para_env%max(max_eigenvalue_deviation)
      success = max_subspace_deviation < subspace_tol .AND. max_eigenvalue_deviation < eigenvalue_tol

      DEALLOCATE (eigenvalue_buffer)
      CALL cp_fm_release(ref_real)
      CALL cp_fm_release(ref_imag)
      CALL cp_fm_release(cand_real)
      CALL cp_fm_release(cand_imag)

   END SUBROUTINE validate_wannier90_reused_mos

! **************************************************************************************************
!> \brief Compare atom/AO reuse candidates directly to the full-mesh reference MOs.
!> \param kpoint full Wannier90 export k-point object holding reference MOs
!> \param qs_kpoint SCF k-point object
!> \param matrix_s real-space overlap matrix
!> \param matrix_ks real-space Kohn-Sham matrix
!> \param cell_to_index real-space cell index table
!> \param sab_nl overlap neighbor list
!> \param para_env global parallel environment
!> \param iw output unit
!> \param max_subspace_deviation largest best-candidate subspace deviation
!> \param min_svalue smallest best-candidate singular value
!> \param max_metric_deviation largest S(k)-metric deviation of a candidate
!> \param max_residual largest H(k),S(k) eigen-residual of a candidate
! **************************************************************************************************
   SUBROUTINE diagnose_wannier90_scf_reuse_candidates(kpoint, qs_kpoint, matrix_s, matrix_ks, &
                                                      cell_to_index, sab_nl, para_env, iw, &
                                                      max_subspace_deviation, min_svalue, &
                                                      max_metric_deviation, max_residual)
      TYPE(kpoint_type), POINTER                         :: kpoint, qs_kpoint
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s, matrix_ks
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: iw
      REAL(KIND=dp), INTENT(OUT)                         :: max_subspace_deviation, min_svalue, &
                                                            max_metric_deviation, max_residual

      REAL(KIND=dp), PARAMETER                           :: print_tol = 1.0e-4_dp, &
                                                            residual_print_tol = 1.0e-3_dp

      CHARACTER(LEN=default_string_length)               :: reason
      INTEGER                                            :: ik, ikpgr, ikred, ispin, isym_try, nao, &
                                                            nao_src, nkp, nmo, nmo_src, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: source_kpoint, sym_index
      INTEGER, DIMENSION(2)                              :: kp_range, source_kp_range
      LOGICAL                                            :: my_kpgrp, ok, source_window
      REAL(KIND=dp) :: best_deviation, best_metric_deviation, best_residual, best_svalue, &
         candidate_deviation, candidate_metric_deviation, candidate_metric_min, &
         candidate_residual, candidate_svalue, owner_count, ref_metric_deviation, ref_metric_min, &
         ref_residual
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues_buffer
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_source, matrix_struct_work
      TYPE(cp_fm_type)                                   :: dst_imag, dst_real, ref_imag, ref_real, &
                                                            src_imag, src_imag_full, src_real, &
                                                            src_real_full
      TYPE(cp_fm_type), POINTER                          :: fmi, fmr, src_fmi, src_fmr
      TYPE(kpoint_env_type), POINTER                     :: kp, kp_source
      TYPE(kpoint_sym_type), POINTER                     :: kpsym

      max_subspace_deviation = HUGE(1.0_dp)
      min_svalue = 0.0_dp
      max_metric_deviation = HUGE(1.0_dp)
      max_residual = HUGE(1.0_dp)
      NULLIFY (matrix_struct_source, matrix_struct_work, fmi, fmr, src_fmi, src_fmr)

      CALL build_wannier90_scf_mapping(kpoint, qs_kpoint, source_kpoint, sym_index, ok, reason)
      IF (.NOT. ok) RETURN
      kp => kpoint%kp_env(1)%kpoint_env
      kp_source => qs_kpoint%kp_env(1)%kpoint_env
      IF (SIZE(kp%mos, 1) < 2 .OR. SIZE(kp_source%mos, 1) < 2) THEN
         DEALLOCATE (source_kpoint, sym_index)
         RETURN
      END IF
      nspins = SIZE(kp%mos, 2)
      CALL get_mo_set(kp%mos(1, 1), nao=nao, nmo=nmo)
      CALL get_mo_set(kp_source%mos(1, 1), nao=nao_src, nmo=nmo_src)
      CALL para_env%max(nao_src)
      CALL para_env%max(nmo_src)
      IF (nao_src /= nao .OR. nmo_src < nmo) THEN
         DEALLOCATE (source_kpoint, sym_index)
         RETURN
      END IF
      source_window = nmo_src > nmo
      CALL get_kpoint_info(kpoint, nkp=nkp, kp_range=kp_range)
      CALL get_kpoint_info(qs_kpoint, kp_range=source_kp_range)
      IF (source_kp_range(1) /= 1 .OR. source_kp_range(2) /= qs_kpoint%nkp) THEN
         DEALLOCATE (source_kpoint, sym_index)
         RETURN
      END IF

      CALL cp_fm_struct_create(matrix_struct_work, nrow_global=nao, ncol_global=nmo, &
                               para_env=para_env, context=kpoint%blacs_env)
      CALL cp_fm_create(ref_real, matrix_struct_work)
      CALL cp_fm_create(ref_imag, matrix_struct_work)
      CALL cp_fm_create(src_real, matrix_struct_work)
      CALL cp_fm_create(src_imag, matrix_struct_work)
      CALL cp_fm_create(dst_real, matrix_struct_work)
      CALL cp_fm_create(dst_imag, matrix_struct_work)
      ALLOCATE (eigenvalues_buffer(nmo))
      IF (source_window) THEN
         CALL cp_fm_struct_create(matrix_struct_source, nrow_global=nao, ncol_global=nmo_src, &
                                  para_env=para_env, context=kpoint%blacs_env)
         CALL cp_fm_create(src_real_full, matrix_struct_source)
         CALL cp_fm_create(src_imag_full, matrix_struct_source)
      END IF

      max_subspace_deviation = 0.0_dp
      min_svalue = HUGE(1.0_dp)
      max_metric_deviation = 0.0_dp
      max_residual = 0.0_dp
      DO ispin = 1, nspins
         DO ik = 1, nkp
            my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
            IF (my_kpgrp) THEN
               ikpgr = ik - kp_range(1) + 1
               kp => kpoint%kp_env(ikpgr)%kpoint_env
               fmr => kp%mos(1, ispin)%mo_coeff
               fmi => kp%mos(2, ispin)%mo_coeff
            ELSE
               NULLIFY (fmr, fmi)
            END IF
            CALL cp_fm_copy_general(fmr, ref_real, para_env)
            CALL cp_fm_copy_general(fmi, ref_imag, para_env)

            ikred = source_kpoint(ik)
            my_kpgrp = (ikred >= source_kp_range(1) .AND. ikred <= source_kp_range(2))
            IF (my_kpgrp) THEN
               ikpgr = ikred - source_kp_range(1) + 1
               kp_source => qs_kpoint%kp_env(ikpgr)%kpoint_env
               src_fmr => kp_source%mos(1, ispin)%mo_coeff
               src_fmi => kp_source%mos(2, ispin)%mo_coeff
               CALL get_mo_set(kp_source%mos(1, ispin), eigenvalues=eigenvalues)
               eigenvalues_buffer(1:nmo) = eigenvalues(1:nmo)
            ELSE
               NULLIFY (src_fmr, src_fmi)
               eigenvalues_buffer(1:nmo) = 0.0_dp
            END IF
            IF (my_kpgrp) THEN
               owner_count = 1.0_dp
            ELSE
               owner_count = 0.0_dp
            END IF
            CALL para_env%sum(owner_count)
            CALL para_env%sum(eigenvalues_buffer)
            IF (owner_count > 0.0_dp) eigenvalues_buffer(1:nmo) = eigenvalues_buffer(1:nmo)/owner_count
            CALL measure_wannier90_eigenspace_quality(ref_real, ref_imag, matrix_s, matrix_ks, &
                                                      kpoint%xkp(1:3, ik), cell_to_index, sab_nl, &
                                                      para_env, ispin, eigenvalues_buffer, ok, &
                                                      ref_metric_deviation, ref_metric_min, ref_residual)
            IF (ok .AND. para_env%is_source() .AND. iw > 0 .AND. &
                (ref_metric_deviation > print_tol .OR. ref_residual > residual_print_tol)) THEN
               WRITE (iw, '(T2,A,I0,A,ES10.3,A,ES10.3,A,ES10.3)') &
                  "WANNIER90| reference k=", ik, " dM=", ref_metric_deviation, &
                  " smin=", ref_metric_min, " resid=", ref_residual
            END IF
            IF (source_window) THEN
               CALL cp_fm_copy_general(src_fmr, src_real_full, para_env)
               CALL cp_fm_copy_general(src_fmi, src_imag_full, para_env)
               CALL copy_wannier90_mo_window(src_real_full, src_real, nmo)
               CALL copy_wannier90_mo_window(src_imag_full, src_imag, nmo)
            ELSE
               CALL cp_fm_copy_general(src_fmr, src_real, para_env)
               CALL cp_fm_copy_general(src_fmi, src_imag, para_env)
            END IF

            best_deviation = HUGE(1.0_dp)
            best_metric_deviation = 0.0_dp
            best_residual = 0.0_dp
            best_svalue = 0.0_dp
            IF (sym_index(ik) <= 0) THEN
               CALL transform_wannier90_scf_mo(src_real, src_imag, dst_real, dst_imag, qs_kpoint, &
                                               ikred, sym_index(ik), para_env, ok, reason)
               IF (ok) THEN
                  CALL measure_wannier90_subspace_error(ref_real, ref_imag, dst_real, dst_imag, matrix_s, &
                                                        kpoint%xkp(1:3, ik), cell_to_index, sab_nl, &
                                                        para_env, ok, candidate_deviation, candidate_svalue)
                  IF (ok) THEN
                     CALL measure_wannier90_eigenspace_quality(dst_real, dst_imag, matrix_s, matrix_ks, &
                                                               kpoint%xkp(1:3, ik), cell_to_index, sab_nl, &
                                                               para_env, ispin, eigenvalues_buffer, ok, &
                                                               candidate_metric_deviation, &
                                                               candidate_metric_min, candidate_residual)
                  END IF
                  IF (ok) THEN
                     best_deviation = candidate_deviation
                     best_metric_deviation = candidate_metric_deviation
                     best_residual = candidate_residual
                     best_svalue = candidate_svalue
                     IF (para_env%is_source() .AND. iw > 0 .AND. &
                         (candidate_deviation > print_tol .OR. candidate_metric_deviation > print_tol .OR. &
                          candidate_residual > residual_print_tol)) THEN
                        WRITE (iw, '(T2,A,I0,A,I0,A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') &
                           "WANNIER90| reuse candidate k=", ik, " src=", ikred, " sym=", &
                           sym_index(ik), " dRef=", candidate_deviation, " dM=", &
                           candidate_metric_deviation, " smin=", candidate_metric_min, &
                           " resid=", candidate_residual
                     END IF
                  END IF
               END IF
            ELSEIF (ASSOCIATED(qs_kpoint%kp_sym)) THEN
               kpsym => qs_kpoint%kp_sym(ikred)%kpoint_sym
               IF (ASSOCIATED(kpsym)) THEN
                  DO isym_try = 1, kpsym%nwred
                     IF (.NOT. same_periodic_kpoint(kpoint%xkp(1:3, ik), &
                                                    kpsym%xkp(1:3, isym_try))) CYCLE
                     CALL transform_wannier90_scf_mo(src_real, src_imag, dst_real, dst_imag, &
                                                     qs_kpoint, ikred, isym_try, para_env, ok, reason)
                     IF (.NOT. ok) CYCLE
                     CALL measure_wannier90_subspace_error(ref_real, ref_imag, dst_real, dst_imag, &
                                                           matrix_s, kpoint%xkp(1:3, ik), cell_to_index, &
                                                           sab_nl, para_env, ok, candidate_deviation, &
                                                           candidate_svalue)
                     IF (ok) THEN
                        CALL measure_wannier90_eigenspace_quality(dst_real, dst_imag, matrix_s, matrix_ks, &
                                                                  kpoint%xkp(1:3, ik), cell_to_index, &
                                                                  sab_nl, para_env, ispin, &
                                                                  eigenvalues_buffer, ok, &
                                                                  candidate_metric_deviation, &
                                                                  candidate_metric_min, candidate_residual)
                     END IF
                     IF (ok .AND. candidate_deviation < best_deviation) THEN
                        best_deviation = candidate_deviation
                        best_metric_deviation = candidate_metric_deviation
                        best_residual = candidate_residual
                        best_svalue = candidate_svalue
                     END IF
                  END DO
               END IF
            END IF
            IF (best_deviation < HUGE(1.0_dp)) THEN
               max_subspace_deviation = MAX(max_subspace_deviation, best_deviation)
               min_svalue = MIN(min_svalue, best_svalue)
               max_metric_deviation = MAX(max_metric_deviation, best_metric_deviation)
               max_residual = MAX(max_residual, best_residual)
            END IF
         END DO
      END DO
      CALL para_env%max(max_subspace_deviation)
      CALL para_env%min(min_svalue)
      CALL para_env%max(max_metric_deviation)
      CALL para_env%max(max_residual)

      IF (source_window) THEN
         CALL cp_fm_release(src_real_full)
         CALL cp_fm_release(src_imag_full)
         CALL cp_fm_struct_release(matrix_struct_source)
      END IF
      CALL cp_fm_release(ref_real)
      CALL cp_fm_release(ref_imag)
      CALL cp_fm_release(src_real)
      CALL cp_fm_release(src_imag)
      CALL cp_fm_release(dst_real)
      CALL cp_fm_release(dst_imag)
      CALL cp_fm_struct_release(matrix_struct_work)
      DEALLOCATE (eigenvalues_buffer)
      DEALLOCATE (source_kpoint, sym_index)

   END SUBROUTINE diagnose_wannier90_scf_reuse_candidates

! **************************************************************************************************
!> \brief Measure the S(k)-metric distance between two Wannier90 MO subspaces.
!> \param ref_real real part of reference MO coefficients
!> \param ref_imag imaginary part of reference MO coefficients
!> \param cand_real real part of candidate MO coefficients
!> \param cand_imag imaginary part of candidate MO coefficients
!> \param matrix_s real-space overlap matrix
!> \param xkp target k-point coordinate
!> \param cell_to_index real-space cell index table
!> \param sab_nl overlap neighbor list
!> \param para_env global parallel environment
!> \param success true if the metric comparison was performed
!> \param max_subspace_deviation largest deviation of singular values from one
!> \param min_svalue smallest singular value of C_ref^+ S(k) C_candidate
! **************************************************************************************************
   SUBROUTINE measure_wannier90_subspace_error(ref_real, ref_imag, cand_real, cand_imag, matrix_s, &
                                               xkp, cell_to_index, sab_nl, para_env, success, &
                                               max_subspace_deviation, min_svalue)
      TYPE(cp_fm_type), INTENT(IN)                       :: ref_real, ref_imag, cand_real, cand_imag
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(OUT)                               :: success
      REAL(KIND=dp), INTENT(OUT)                         :: max_subspace_deviation, min_svalue

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: metric_projected, metric_vectors, &
                                                            overlap, ref_coeff, s_cand
      INTEGER                                            :: ib, nao, nmo, nmo_candidate
      REAL(KIND=dp)                                      :: singular_value
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: metric_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ref_i, ref_r, s_cand_i, s_cand_r
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_metric
      TYPE(cp_fm_type)                                   :: s_cand_imag, s_cand_real

      success = .FALSE.
      max_subspace_deviation = HUGE(1.0_dp)
      min_svalue = 0.0_dp
      NULLIFY (matrix_struct_metric)

      CALL cp_fm_get_info(ref_real, nrow_global=nao, ncol_global=nmo, &
                          matrix_struct=matrix_struct_metric)
      CALL cp_fm_get_info(cand_real, ncol_global=nmo_candidate)
      IF (nmo_candidate /= nmo) RETURN

      CALL cp_fm_create(s_cand_real, matrix_struct_metric)
      CALL cp_fm_create(s_cand_imag, matrix_struct_metric)
      CALL apply_wannier90_kp_matrix(matrix_s, 1, xkp, cell_to_index, sab_nl, &
                                     cand_real, cand_imag, s_cand_real, s_cand_imag)

      ALLOCATE (ref_r(nao, nmo), ref_i(nao, nmo), s_cand_r(nao, nmo), s_cand_i(nao, nmo))
      CALL cp_fm_get_submatrix(ref_real, ref_r)
      CALL cp_fm_get_submatrix(ref_imag, ref_i)
      CALL cp_fm_get_submatrix(s_cand_real, s_cand_r)
      CALL cp_fm_get_submatrix(s_cand_imag, s_cand_i)

      ALLOCATE (ref_coeff(nao, nmo), s_cand(nao, nmo), overlap(nmo, nmo), &
                metric_projected(nmo, nmo), metric_vectors(nmo, nmo), metric_values(nmo))
      ref_coeff(:, :) = CMPLX(ref_r, ref_i, KIND=dp)
      s_cand(:, :) = CMPLX(s_cand_r, s_cand_i, KIND=dp)
      overlap(:, :) = MATMUL(CONJG(TRANSPOSE(ref_coeff)), s_cand)
      metric_projected(:, :) = MATMUL(CONJG(TRANSPOSE(overlap)), overlap)
      metric_projected(:, :) = 0.5_dp*(metric_projected + CONJG(TRANSPOSE(metric_projected)))
      CALL diag_complex(metric_projected, metric_vectors, metric_values)

      min_svalue = HUGE(1.0_dp)
      max_subspace_deviation = 0.0_dp
      DO ib = 1, nmo
         singular_value = SQRT(MAX(metric_values(ib), 0.0_dp))
         min_svalue = MIN(min_svalue, singular_value)
         max_subspace_deviation = MAX(max_subspace_deviation, ABS(singular_value - 1.0_dp))
      END DO
      CALL para_env%max(max_subspace_deviation)
      CALL para_env%min(min_svalue)
      success = .TRUE.

      DEALLOCATE (ref_coeff, s_cand, overlap, metric_projected, metric_vectors, metric_values)
      DEALLOCATE (ref_r, ref_i, s_cand_r, s_cand_i)
      CALL cp_fm_release(s_cand_real)
      CALL cp_fm_release(s_cand_imag)

   END SUBROUTINE measure_wannier90_subspace_error

! **************************************************************************************************
!> \brief Measure whether transformed MOs are an H(k),S(k) invariant eigenspace.
!> \param cand_real real part of candidate MO coefficients
!> \param cand_imag imaginary part of candidate MO coefficients
!> \param matrix_s real-space overlap matrix
!> \param matrix_ks real-space Kohn-Sham matrix
!> \param xkp target k-point coordinate
!> \param cell_to_index real-space cell index table
!> \param sab_nl overlap neighbor list
!> \param para_env global parallel environment
!> \param ispin spin index
!> \param eigenvalues source MO eigenvalues corresponding to the candidate columns
!> \param success true if the metric and residual checks were performed
!> \param metric_deviation largest deviation of eigenvalues of C^+ S(k) C from one
!> \param min_metric_eigenvalue smallest eigenvalue of C^+ S(k) C
!> \param residual_norm largest element of H(k) C - S(k) C eps
! **************************************************************************************************
   SUBROUTINE measure_wannier90_eigenspace_quality(cand_real, cand_imag, matrix_s, matrix_ks, xkp, &
                                                   cell_to_index, sab_nl, para_env, ispin, eigenvalues, &
                                                   success, metric_deviation, min_metric_eigenvalue, &
                                                   residual_norm)
      TYPE(cp_fm_type), INTENT(IN)                       :: cand_real, cand_imag
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s, matrix_ks
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eigenvalues
      LOGICAL, INTENT(OUT)                               :: success
      REAL(KIND=dp), INTENT(OUT)                         :: metric_deviation, min_metric_eigenvalue, &
                                                            residual_norm

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: cand_coeff, h_coeff, metric_vectors, &
                                                            residual_block, s_coeff, s_projected
      INTEGER                                            :: ib, nao, nmo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: metric_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cand_i, cand_r, h_coeff_i, h_coeff_r, &
                                                            s_coeff_i, s_coeff_r
      TYPE(cp_cfm_type)                                  :: cand_cfm, metric_cfm, s_cfm
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_metric, &
                                                            matrix_struct_projected
      TYPE(cp_fm_type)                                   :: h_cand_imag, h_cand_real, s_cand_imag, &
                                                            s_cand_real, tmp_fm

      success = .FALSE.
      metric_deviation = HUGE(1.0_dp)
      min_metric_eigenvalue = 0.0_dp
      residual_norm = HUGE(1.0_dp)
      NULLIFY (matrix_struct_metric, matrix_struct_projected)

      CALL cp_fm_get_info(cand_real, nrow_global=nao, ncol_global=nmo, &
                          matrix_struct=matrix_struct_metric)
      IF (SIZE(eigenvalues) < nmo) RETURN

      CALL cp_fm_create(s_cand_real, matrix_struct_metric)
      CALL cp_fm_create(s_cand_imag, matrix_struct_metric)
      CALL cp_fm_create(h_cand_real, matrix_struct_metric)
      CALL cp_fm_create(h_cand_imag, matrix_struct_metric)
      CALL cp_fm_create(tmp_fm, matrix_struct_metric)
      CALL cp_cfm_create(cand_cfm, matrix_struct_metric)
      CALL cp_cfm_create(s_cfm, matrix_struct_metric)

      CALL apply_wannier90_kp_matrix(matrix_s, 1, xkp, cell_to_index, sab_nl, &
                                     cand_real, cand_imag, s_cand_real, s_cand_imag)
      CALL apply_wannier90_kp_matrix(matrix_ks, ispin, xkp, cell_to_index, sab_nl, &
                                     cand_real, cand_imag, h_cand_real, h_cand_imag)

      ALLOCATE (cand_r(nao, nmo), cand_i(nao, nmo), s_coeff_r(nao, nmo), &
                s_coeff_i(nao, nmo), h_coeff_r(nao, nmo), h_coeff_i(nao, nmo))
      CALL cp_fm_get_submatrix(cand_real, cand_r)
      CALL cp_fm_get_submatrix(cand_imag, cand_i)
      CALL cp_fm_get_submatrix(s_cand_real, s_coeff_r)
      CALL cp_fm_get_submatrix(s_cand_imag, s_coeff_i)
      CALL cp_fm_get_submatrix(h_cand_real, h_coeff_r)
      CALL cp_fm_get_submatrix(h_cand_imag, h_coeff_i)

      ALLOCATE (cand_coeff(nao, nmo), h_coeff(nao, nmo), metric_vectors(nmo, nmo), &
                residual_block(nao, nmo), s_coeff(nao, nmo), s_projected(nmo, nmo), &
                metric_values(nmo))
      cand_coeff(:, :) = CMPLX(cand_r, cand_i, KIND=dp)
      s_coeff(:, :) = CMPLX(s_coeff_r, s_coeff_i, KIND=dp)
      h_coeff(:, :) = CMPLX(h_coeff_r, h_coeff_i, KIND=dp)
      CALL cp_fm_struct_create(matrix_struct_projected, nrow_global=nmo, ncol_global=nmo, &
                               para_env=matrix_struct_metric%para_env, &
                               context=matrix_struct_metric%context)
      CALL cp_cfm_create(metric_cfm, matrix_struct_projected)
      CALL cp_fm_to_cfm(cand_real, cand_imag, cand_cfm)
      CALL cp_fm_to_cfm(s_cand_real, s_cand_imag, s_cfm)
      CALL cp_cfm_gemm("C", "N", nmo, nmo, nao, CMPLX(1.0_dp, 0.0_dp, KIND=dp), cand_cfm, &
                       s_cfm, CMPLX(0.0_dp, 0.0_dp, KIND=dp), metric_cfm)
      CALL cp_cfm_get_submatrix(metric_cfm, s_projected)
      s_projected(:, :) = 0.5_dp*(s_projected + CONJG(TRANSPOSE(s_projected)))
      CALL diag_complex(s_projected, metric_vectors, metric_values)
      metric_deviation = MAXVAL(ABS(metric_values - 1.0_dp))
      min_metric_eigenvalue = MINVAL(metric_values)

      residual_block(:, :) = h_coeff
      DO ib = 1, nmo
         residual_block(:, ib) = residual_block(:, ib) - eigenvalues(ib)*s_coeff(:, ib)
      END DO
      residual_norm = MAXVAL(ABS(residual_block))
      CALL para_env%max(metric_deviation)
      CALL para_env%min(min_metric_eigenvalue)
      CALL para_env%max(residual_norm)
      success = .TRUE.

      DEALLOCATE (cand_coeff, h_coeff, metric_vectors, residual_block, s_coeff, s_projected, &
                  metric_values)
      DEALLOCATE (cand_r, cand_i, s_coeff_r, s_coeff_i, h_coeff_r, h_coeff_i)
      CALL cp_fm_release(s_cand_real)
      CALL cp_fm_release(s_cand_imag)
      CALL cp_fm_release(h_cand_real)
      CALL cp_fm_release(h_cand_imag)
      CALL cp_fm_release(tmp_fm)
      CALL cp_cfm_release(cand_cfm)
      CALL cp_cfm_release(s_cfm)
      CALL cp_cfm_release(metric_cfm)
      CALL cp_fm_struct_release(matrix_struct_projected)

   END SUBROUTINE measure_wannier90_eigenspace_quality

! **************************************************************************************************
!> \brief Copy the leading MO columns from a larger SCF MO matrix into the Wannier90 export window.
!> \param source source MO coefficient matrix
!> \param destination destination MO coefficient matrix
!> \param ncol number of columns to copy
! **************************************************************************************************
   SUBROUTINE copy_wannier90_mo_window(source, destination, ncol)
      TYPE(cp_fm_type), INTENT(IN)                       :: source, destination
      INTEGER, INTENT(IN)                                :: ncol

      INTEGER                                            :: ncol_source, nrow
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: destination_buffer, source_buffer

      CALL cp_fm_get_info(source, nrow_global=nrow, ncol_global=ncol_source)
      CPASSERT(ncol_source >= ncol)
      ALLOCATE (source_buffer(nrow, ncol_source), destination_buffer(nrow, ncol))
      CALL cp_fm_get_submatrix(source, source_buffer)
      destination_buffer(1:nrow, 1:ncol) = source_buffer(1:nrow, 1:ncol)
      CALL cp_fm_set_submatrix(destination, destination_buffer)
      DEALLOCATE (source_buffer, destination_buffer)

   END SUBROUTINE copy_wannier90_mo_window

! **************************************************************************************************
!> \brief Apply a complex k-point matrix to a complex MO coefficient matrix.
!> \param rsmat real-space matrix images
!> \param ispin spin index for rsmat
!> \param xkp target k-point coordinate
!> \param cell_to_index real-space cell index table
!> \param sab_nl overlap neighbor list
!> \param coeff_real real part of input MO coefficients
!> \param coeff_imag imaginary part of input MO coefficients
!> \param result_real real part of matrix-vector product
!> \param result_imag imaginary part of matrix-vector product
! **************************************************************************************************
   SUBROUTINE apply_wannier90_kp_matrix(rsmat, ispin, xkp, cell_to_index, sab_nl, &
                                        coeff_real, coeff_imag, result_real, result_imag)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(cp_fm_type), INTENT(IN)                       :: coeff_real, coeff_imag, result_real, &
                                                            result_imag

      INTEGER                                            :: nao, ncol
      TYPE(cp_cfm_type)                                  :: coeff_cfm, kmat_cfm, result_cfm
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_ao, matrix_struct_coeff
      TYPE(cp_fm_type)                                   :: mat_imag, mat_real
      TYPE(dbcsr_type), POINTER                          :: kmat_imag, kmat_imag_full, kmat_real, &
                                                            kmat_real_full

      NULLIFY (matrix_struct_ao, matrix_struct_coeff, kmat_imag, kmat_imag_full, kmat_real, &
               kmat_real_full)

      CALL cp_fm_get_info(coeff_real, nrow_global=nao, ncol_global=ncol, &
                          matrix_struct=matrix_struct_coeff)

      ALLOCATE (kmat_real, kmat_imag, kmat_real_full, kmat_imag_full)
      CALL dbcsr_create(kmat_real, template=rsmat(ispin, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(kmat_imag, template=rsmat(ispin, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL dbcsr_create(kmat_real_full, template=rsmat(ispin, 1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(kmat_imag_full, template=rsmat(ispin, 1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_alloc_block_from_nbl(kmat_real, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(kmat_imag, sab_nl)
      CALL dbcsr_set(kmat_real, 0.0_dp)
      CALL dbcsr_set(kmat_imag, 0.0_dp)
      CALL rskp_transform(kmat_real, kmat_imag, rsmat=rsmat, ispin=ispin, &
                          xkp=xkp, cell_to_index=cell_to_index, sab_nl=sab_nl)
      CALL dbcsr_desymmetrize(kmat_real, kmat_real_full)
      CALL dbcsr_desymmetrize(kmat_imag, kmat_imag_full)

      CALL cp_fm_struct_create(matrix_struct_ao, nrow_global=nao, ncol_global=nao, &
                               para_env=matrix_struct_coeff%para_env, &
                               context=matrix_struct_coeff%context)
      CALL cp_fm_create(mat_real, matrix_struct_ao)
      CALL cp_fm_create(mat_imag, matrix_struct_ao)
      CALL copy_dbcsr_to_fm(kmat_real_full, mat_real)
      CALL copy_dbcsr_to_fm(kmat_imag_full, mat_imag)

      CALL cp_cfm_create(kmat_cfm, matrix_struct_ao)
      CALL cp_cfm_create(coeff_cfm, matrix_struct_coeff)
      CALL cp_cfm_create(result_cfm, matrix_struct_coeff)
      CALL cp_fm_to_cfm(mat_real, mat_imag, kmat_cfm)
      CALL cp_fm_to_cfm(coeff_real, coeff_imag, coeff_cfm)
      CALL cp_cfm_gemm("N", "N", nao, ncol, nao, CMPLX(1.0_dp, 0.0_dp, KIND=dp), kmat_cfm, &
                       coeff_cfm, CMPLX(0.0_dp, 0.0_dp, KIND=dp), result_cfm)
      CALL cp_cfm_to_fm(result_cfm, result_real, result_imag)

      CALL cp_fm_release(mat_real)
      CALL cp_fm_release(mat_imag)
      CALL cp_cfm_release(kmat_cfm)
      CALL cp_cfm_release(coeff_cfm)
      CALL cp_cfm_release(result_cfm)
      CALL cp_fm_struct_release(matrix_struct_ao)
      CALL dbcsr_deallocate_matrix(kmat_real)
      CALL dbcsr_deallocate_matrix(kmat_imag)
      CALL dbcsr_deallocate_matrix(kmat_real_full)
      CALL dbcsr_deallocate_matrix(kmat_imag_full)

   END SUBROUTINE apply_wannier90_kp_matrix

! **************************************************************************************************
!> \brief Rayleigh-Ritz stabilize a symmetry-reconstructed Wannier90 MO subspace.
!> \param dst_real real part of transformed MO coefficients
!> \param dst_imag imaginary part of transformed MO coefficients
!> \param matrix_s real-space overlap matrix
!> \param matrix_ks real-space Kohn-Sham matrix
!> \param xkp target k-point coordinate
!> \param cell_to_index real-space cell index table
!> \param sab_nl overlap neighbor list
!> \param ispin spin index
!> \param eigenvalues Ritz eigenvalues of the stabilized subspace
!> \param degenerate_band_tol degeneracy threshold
!> \param success true if the subspace was stabilized
!> \param reason diagnostic message
!> \param aligned_blocks number of stabilized subspaces
!> \param aligned_max_size largest stabilized subspace
!> \param aligned_min_svalue smallest S(k)-metric eigenvalue
!> \param max_residual largest Ritz residual
! **************************************************************************************************
   SUBROUTINE ritz_stabilize_wannier90_subspace(dst_real, dst_imag, matrix_s, matrix_ks, &
                                                xkp, cell_to_index, sab_nl, ispin, eigenvalues, &
                                                degenerate_band_tol, success, reason, aligned_blocks, &
                                                aligned_max_size, aligned_min_svalue, max_residual)
      TYPE(cp_fm_type), INTENT(IN)                       :: dst_real, dst_imag
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s, matrix_ks
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: eigenvalues
      REAL(KIND=dp), INTENT(IN)                          :: degenerate_band_tol
      LOGICAL, INTENT(OUT)                               :: success
      CHARACTER(LEN=*), INTENT(OUT)                      :: reason
      INTEGER, INTENT(OUT)                               :: aligned_blocks, aligned_max_size
      REAL(KIND=dp), INTENT(OUT)                         :: aligned_min_svalue, max_residual

      REAL(KIND=dp), PARAMETER                           :: residual_tol = 1.0e-2_dp

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: block_coeff, h_block, h_coeff, &
         h_projected, h_projected_work, metric_vectors, residual_block, ritz_vectors, s_block, &
         s_coeff, s_projected, stabilized
      INTEGER                                            :: block_first, block_last, block_size, ib, &
                                                            nao, nmo
      REAL(KIND=dp)                                      :: metric_deviation, norm_value, &
                                                            residual_norm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: metric_values, ritz_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dst_i, dst_r, h_coeff_i, h_coeff_r, &
                                                            s_coeff_i, s_coeff_r
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_metric
      TYPE(cp_fm_type)                                   :: h_dst_imag, h_dst_real, s_dst_imag, &
                                                            s_dst_real, tmp_fm

      success = .FALSE.
      reason = ""
      aligned_blocks = 0
      aligned_max_size = 0
      aligned_min_svalue = HUGE(1.0_dp)
      max_residual = 0.0_dp

      NULLIFY (matrix_struct_metric)
      CALL cp_fm_get_info(dst_real, nrow_global=nao, ncol_global=nmo, &
                          matrix_struct=matrix_struct_metric)
      IF (SIZE(eigenvalues) < nmo) THEN
         reason = "not enough eigenvalues for Wannier90 Ritz subspace stabilization"
         RETURN
      END IF

      CALL cp_fm_create(s_dst_real, matrix_struct_metric)
      CALL cp_fm_create(s_dst_imag, matrix_struct_metric)
      CALL cp_fm_create(h_dst_real, matrix_struct_metric)
      CALL cp_fm_create(h_dst_imag, matrix_struct_metric)
      CALL cp_fm_create(tmp_fm, matrix_struct_metric)

      CALL apply_wannier90_kp_matrix(matrix_s, 1, xkp, cell_to_index, sab_nl, &
                                     dst_real, dst_imag, s_dst_real, s_dst_imag)
      CALL apply_wannier90_kp_matrix(matrix_ks, ispin, xkp, cell_to_index, sab_nl, &
                                     dst_real, dst_imag, h_dst_real, h_dst_imag)

      ALLOCATE (dst_r(nao, nmo), dst_i(nao, nmo), s_coeff_r(nao, nmo), s_coeff_i(nao, nmo), &
                h_coeff_r(nao, nmo), h_coeff_i(nao, nmo))
      CALL cp_fm_get_submatrix(dst_real, dst_r)
      CALL cp_fm_get_submatrix(dst_imag, dst_i)
      CALL cp_fm_get_submatrix(s_dst_real, s_coeff_r)
      CALL cp_fm_get_submatrix(s_dst_imag, s_coeff_i)
      CALL cp_fm_get_submatrix(h_dst_real, h_coeff_r)
      CALL cp_fm_get_submatrix(h_dst_imag, h_coeff_i)

      ALLOCATE (s_coeff(nao, nmo), h_coeff(nao, nmo))
      s_coeff(:, :) = CMPLX(s_coeff_r, s_coeff_i, KIND=dp)
      h_coeff(:, :) = CMPLX(h_coeff_r, h_coeff_i, KIND=dp)

      block_first = 1
      DO WHILE (block_first <= nmo)
         block_last = block_first
         DO WHILE (block_last < nmo)
            IF (ABS(eigenvalues(block_last + 1) - eigenvalues(block_last)) >= degenerate_band_tol) EXIT
            block_last = block_last + 1
         END DO
         block_size = block_last - block_first + 1
         IF (block_size > 1) THEN
            ! The atom/AO operation fixes the subspace, while the little-group gauge inside an
            ! exactly degenerate manifold is arbitrary. Stabilize only that manifold and verify
            ! that it is an invariant H(k),S(k) subspace before exporting it to Wannier90.
            ALLOCATE (block_coeff(nao, block_size), h_block(nao, block_size), &
                      h_projected(block_size, block_size), h_projected_work(block_size, block_size), &
                      metric_vectors(block_size, block_size), residual_block(nao, block_size), &
                      ritz_vectors(block_size, block_size), s_block(nao, block_size), &
                      s_projected(block_size, block_size), stabilized(nao, block_size), &
                      metric_values(block_size), ritz_values(block_size))
            block_coeff(:, :) = CMPLX(dst_r(:, block_first:block_last), &
                                      dst_i(:, block_first:block_last), KIND=dp)
            s_block(:, :) = s_coeff(:, block_first:block_last)
            h_block(:, :) = h_coeff(:, block_first:block_last)
            s_projected(:, :) = MATMUL(CONJG(TRANSPOSE(block_coeff)), s_block)
            h_projected(:, :) = MATMUL(CONJG(TRANSPOSE(block_coeff)), h_block)
            s_projected(:, :) = 0.5_dp*(s_projected + CONJG(TRANSPOSE(s_projected)))
            h_projected(:, :) = 0.5_dp*(h_projected + CONJG(TRANSPOSE(h_projected)))

            CALL diag_complex(s_projected, metric_vectors, metric_values)
            aligned_min_svalue = MIN(aligned_min_svalue, MINVAL(metric_values))
            metric_deviation = MAXVAL(ABS(metric_values - 1.0_dp))
            IF (MINVAL(metric_values) < 1.0e-10_dp) THEN
               WRITE (reason, "(A,I0,A,ES9.2,A,ES9.2)") &
                  "singular metric blk=", block_first, " smin=", MINVAL(metric_values), &
                  " dS=", metric_deviation
               max_residual = HUGE(1.0_dp)
               DEALLOCATE (block_coeff, h_block, h_projected, h_projected_work, metric_vectors, &
                           residual_block, ritz_vectors, s_block, s_projected, stabilized, &
                           metric_values, ritz_values)
               DEALLOCATE (s_coeff, h_coeff, dst_r, dst_i, s_coeff_r, s_coeff_i, h_coeff_r, &
                           h_coeff_i)
               CALL cp_fm_release(s_dst_real)
               CALL cp_fm_release(s_dst_imag)
               CALL cp_fm_release(h_dst_real)
               CALL cp_fm_release(h_dst_imag)
               CALL cp_fm_release(tmp_fm)
               RETURN
            END IF

            DO ib = 1, block_size
               metric_vectors(:, ib) = metric_vectors(:, ib)/SQRT(metric_values(ib))
            END DO
            h_projected_work(:, :) = MATMUL(h_projected, metric_vectors)
            h_projected(:, :) = MATMUL(CONJG(TRANSPOSE(metric_vectors)), h_projected_work)
            h_projected(:, :) = 0.5_dp*(h_projected + CONJG(TRANSPOSE(h_projected)))
            CALL diag_complex(h_projected, ritz_vectors, ritz_values)
            h_projected_work(:, :) = MATMUL(metric_vectors, ritz_vectors)
            ritz_vectors(:, :) = h_projected_work
            stabilized(:, :) = MATMUL(block_coeff, ritz_vectors)
            residual_block(:, :) = MATMUL(h_block, ritz_vectors)
            h_block(:, :) = residual_block
            residual_block(:, :) = MATMUL(s_block, ritz_vectors)
            s_block(:, :) = residual_block
            DO ib = 1, block_size
               norm_value = SQRT(ABS(REAL(DOT_PRODUCT(stabilized(:, ib), s_block(:, ib)), KIND=dp)))
               IF (norm_value > EPSILON(1.0_dp)) THEN
                  stabilized(:, ib) = stabilized(:, ib)/norm_value
                  h_block(:, ib) = h_block(:, ib)/norm_value
                  s_block(:, ib) = s_block(:, ib)/norm_value
               END IF
            END DO
            residual_block(:, :) = h_block
            DO ib = 1, block_size
               residual_block(:, ib) = residual_block(:, ib) - &
                                       eigenvalues(block_first + ib - 1)*s_block(:, ib)
            END DO
            residual_norm = MAXVAL(ABS(residual_block))
            max_residual = MAX(max_residual, residual_norm)
            IF (residual_norm > residual_tol) THEN
               WRITE (reason, "(A,I0,A,ES9.2)") &
                  "blk=", block_first, " dS=", metric_deviation
               DEALLOCATE (block_coeff, h_block, h_projected, h_projected_work, metric_vectors, &
                           residual_block, ritz_vectors, s_block, s_projected, stabilized, &
                           metric_values, ritz_values)
               DEALLOCATE (s_coeff, h_coeff, dst_r, dst_i, s_coeff_r, s_coeff_i, h_coeff_r, &
                           h_coeff_i)
               CALL cp_fm_release(s_dst_real)
               CALL cp_fm_release(s_dst_imag)
               CALL cp_fm_release(h_dst_real)
               CALL cp_fm_release(h_dst_imag)
               CALL cp_fm_release(tmp_fm)
               RETURN
            END IF

            dst_r(:, block_first:block_last) = REAL(stabilized, KIND=dp)
            dst_i(:, block_first:block_last) = AIMAG(stabilized)
            h_coeff(:, block_first:block_last) = h_block
            s_coeff(:, block_first:block_last) = s_block
            aligned_blocks = aligned_blocks + 1
            aligned_max_size = MAX(aligned_max_size, block_size)
            DEALLOCATE (block_coeff, h_block, h_projected, h_projected_work, metric_vectors, &
                        residual_block, ritz_vectors, s_block, s_projected, stabilized, metric_values, &
                        ritz_values)
         END IF
         block_first = block_last + 1
      END DO

      DO ib = 1, nmo
         residual_norm = MAXVAL(ABS(h_coeff(:, ib) - eigenvalues(ib)*s_coeff(:, ib)))
         max_residual = MAX(max_residual, residual_norm)
      END DO
      IF (max_residual > residual_tol) THEN
         WRITE (reason, "(A,ES10.3)") &
            "atom/AO W90 reuse guarded: Ritz residual=", max_residual
         DEALLOCATE (s_coeff, h_coeff, dst_r, dst_i, s_coeff_r, s_coeff_i, h_coeff_r, h_coeff_i)
         CALL cp_fm_release(s_dst_real)
         CALL cp_fm_release(s_dst_imag)
         CALL cp_fm_release(h_dst_real)
         CALL cp_fm_release(h_dst_imag)
         CALL cp_fm_release(tmp_fm)
         RETURN
      END IF

      CALL cp_fm_set_submatrix(dst_real, dst_r)
      CALL cp_fm_set_submatrix(dst_imag, dst_i)

      IF (aligned_blocks == 0) aligned_min_svalue = 0.0_dp
      success = .TRUE.

      DEALLOCATE (s_coeff, h_coeff, dst_r, dst_i, s_coeff_r, s_coeff_i, h_coeff_r, h_coeff_i)
      CALL cp_fm_release(s_dst_real)
      CALL cp_fm_release(s_dst_imag)
      CALL cp_fm_release(h_dst_real)
      CALL cp_fm_release(h_dst_imag)
      CALL cp_fm_release(tmp_fm)

   END SUBROUTINE ritz_stabilize_wannier90_subspace

! **************************************************************************************************
!> \brief Reconstruct a Wannier90 export window from a larger symmetry-transformed SCF MO space.
!> \param src_real real part of the transformed source MO window
!> \param src_imag imaginary part of the transformed source MO window
!> \param dst_real real part of the exported reconstructed MO coefficients
!> \param dst_imag imaginary part of the exported reconstructed MO coefficients
!> \param matrix_s real-space overlap matrix
!> \param matrix_ks real-space Kohn-Sham matrix
!> \param xkp target k-point coordinate
!> \param cell_to_index real-space cell index table
!> \param sab_nl overlap neighbor list
!> \param ispin spin index
!> \param eigenvalues reconstructed target eigenvalues for the exported window
!> \param nmo_export number of MOs to export
!> \param success true if the reconstructed window is an invariant H(k),S(k) subspace
!> \param reason diagnostic message
!> \param min_svalue smallest S(k)-metric eigenvalue in the source window
!> \param max_residual largest target Ritz residual
! **************************************************************************************************
   SUBROUTINE ritz_reconstruct_wannier90_window(src_real, src_imag, dst_real, dst_imag, matrix_s, &
                                                matrix_ks, xkp, cell_to_index, sab_nl, ispin, &
                                                eigenvalues, nmo_export, success, reason, min_svalue, &
                                                max_residual)
      TYPE(cp_fm_type), INTENT(IN)                       :: src_real, src_imag, dst_real, dst_imag
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s, matrix_ks
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: eigenvalues
      INTEGER, INTENT(IN)                                :: nmo_export
      LOGICAL, INTENT(OUT)                               :: success
      CHARACTER(LEN=*), INTENT(OUT)                      :: reason
      REAL(KIND=dp), INTENT(OUT)                         :: min_svalue, max_residual

      REAL(KIND=dp), PARAMETER                           :: eigenvalue_tol = 1.0e-6_dp, &
                                                            residual_tol = 1.0e-7_dp

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: coeff_work, h_coeff, h_projected, &
         h_projected_work, metric_vectors, residual_block, ritz_vectors, s_coeff, s_projected, &
         source_coeff, stabilized
      INTEGER                                            :: ib, nao, nmo_source
      REAL(KIND=dp)                                      :: max_eigenvalue_shift, metric_deviation, &
                                                            norm_value
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: metric_values, ritz_values, &
                                                            source_eigenvalues
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dst_i, dst_r, h_coeff_i, h_coeff_r, &
                                                            s_coeff_i, s_coeff_r, src_i, src_r
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_metric
      TYPE(cp_fm_type)                                   :: h_src_imag, h_src_real, s_src_imag, &
                                                            s_src_real, tmp_fm

      success = .FALSE.
      reason = ""
      min_svalue = HUGE(1.0_dp)
      max_residual = HUGE(1.0_dp)

      NULLIFY (matrix_struct_metric)
      CALL cp_fm_get_info(src_real, nrow_global=nao, ncol_global=nmo_source, &
                          matrix_struct=matrix_struct_metric)
      IF (nmo_export > nmo_source) THEN
         reason = "Wannier90 export window is larger than the transformed SCF MO space"
         RETURN
      END IF
      IF (SIZE(eigenvalues) < nmo_export) THEN
         reason = "not enough eigenvalue storage for Wannier90 source-window reconstruction"
         RETURN
      END IF

      CALL cp_fm_create(s_src_real, matrix_struct_metric)
      CALL cp_fm_create(s_src_imag, matrix_struct_metric)
      CALL cp_fm_create(h_src_real, matrix_struct_metric)
      CALL cp_fm_create(h_src_imag, matrix_struct_metric)
      CALL cp_fm_create(tmp_fm, matrix_struct_metric)

      CALL apply_wannier90_kp_matrix(matrix_s, 1, xkp, cell_to_index, sab_nl, &
                                     src_real, src_imag, s_src_real, s_src_imag)
      CALL apply_wannier90_kp_matrix(matrix_ks, ispin, xkp, cell_to_index, sab_nl, &
                                     src_real, src_imag, h_src_real, h_src_imag)

      ALLOCATE (src_r(nao, nmo_source), src_i(nao, nmo_source), &
                s_coeff_r(nao, nmo_source), s_coeff_i(nao, nmo_source), &
                h_coeff_r(nao, nmo_source), h_coeff_i(nao, nmo_source), &
                dst_r(nao, nmo_export), dst_i(nao, nmo_export))
      CALL cp_fm_get_submatrix(src_real, src_r)
      CALL cp_fm_get_submatrix(src_imag, src_i)
      CALL cp_fm_get_submatrix(s_src_real, s_coeff_r)
      CALL cp_fm_get_submatrix(s_src_imag, s_coeff_i)
      CALL cp_fm_get_submatrix(h_src_real, h_coeff_r)
      CALL cp_fm_get_submatrix(h_src_imag, h_coeff_i)

      ALLOCATE (source_coeff(nao, nmo_source), s_coeff(nao, nmo_source), &
                h_coeff(nao, nmo_source), h_projected(nmo_source, nmo_source), &
                h_projected_work(nmo_source, nmo_source), metric_vectors(nmo_source, nmo_source), &
                residual_block(nao, nmo_export), ritz_vectors(nmo_source, nmo_source), &
                s_projected(nmo_source, nmo_source), stabilized(nao, nmo_export), &
                coeff_work(nao, nmo_source), metric_values(nmo_source), ritz_values(nmo_source), &
                source_eigenvalues(nmo_export))
      source_eigenvalues(1:nmo_export) = eigenvalues(1:nmo_export)
      source_coeff(:, :) = CMPLX(src_r, src_i, KIND=dp)
      s_coeff(:, :) = CMPLX(s_coeff_r, s_coeff_i, KIND=dp)
      h_coeff(:, :) = CMPLX(h_coeff_r, h_coeff_i, KIND=dp)
      s_projected(:, :) = MATMUL(CONJG(TRANSPOSE(source_coeff)), s_coeff)
      h_projected(:, :) = MATMUL(CONJG(TRANSPOSE(source_coeff)), h_coeff)
      s_projected(:, :) = 0.5_dp*(s_projected + CONJG(TRANSPOSE(s_projected)))
      h_projected(:, :) = 0.5_dp*(h_projected + CONJG(TRANSPOSE(h_projected)))

      CALL diag_complex(s_projected, metric_vectors, metric_values)
      min_svalue = MINVAL(metric_values)
      metric_deviation = MAXVAL(ABS(metric_values - 1.0_dp))
      IF (min_svalue < 1.0e-10_dp) THEN
         WRITE (reason, "(A,ES9.2,A,ES9.2)") &
            "singular expanded metric smin=", min_svalue, " dS=", metric_deviation
         GOTO 100
      END IF

      DO ib = 1, nmo_source
         metric_vectors(:, ib) = metric_vectors(:, ib)/SQRT(metric_values(ib))
      END DO
      h_projected_work(:, :) = MATMUL(h_projected, metric_vectors)
      h_projected(:, :) = MATMUL(CONJG(TRANSPOSE(metric_vectors)), h_projected_work)
      h_projected(:, :) = 0.5_dp*(h_projected + CONJG(TRANSPOSE(h_projected)))
      CALL diag_complex(h_projected, ritz_vectors, ritz_values)
      h_projected_work(:, :) = MATMUL(metric_vectors, ritz_vectors)
      ritz_vectors(:, :) = h_projected_work
      stabilized(:, :) = MATMUL(source_coeff, ritz_vectors(:, 1:nmo_export))
      coeff_work(:, :) = MATMUL(h_coeff, ritz_vectors)
      h_coeff(:, :) = coeff_work
      coeff_work(:, :) = MATMUL(s_coeff, ritz_vectors)
      s_coeff(:, :) = coeff_work
      DO ib = 1, nmo_export
         norm_value = SQRT(ABS(REAL(DOT_PRODUCT(stabilized(:, ib), s_coeff(:, ib)), KIND=dp)))
         IF (norm_value > EPSILON(1.0_dp)) THEN
            stabilized(:, ib) = stabilized(:, ib)/norm_value
            h_coeff(:, ib) = h_coeff(:, ib)/norm_value
            s_coeff(:, ib) = s_coeff(:, ib)/norm_value
         END IF
      END DO
      residual_block(:, :) = h_coeff(:, 1:nmo_export)
      DO ib = 1, nmo_export
         residual_block(:, ib) = residual_block(:, ib) - ritz_values(ib)*s_coeff(:, ib)
      END DO
      max_residual = MAXVAL(ABS(residual_block))
      IF (max_residual > residual_tol) THEN
         WRITE (reason, "(A,ES9.2)") &
            "expanded dS=", metric_deviation
         GOTO 100
      END IF
      max_eigenvalue_shift = MAXVAL(ABS(ritz_values(1:nmo_export) - source_eigenvalues(1:nmo_export)))
      IF (max_eigenvalue_shift > eigenvalue_tol) THEN
         WRITE (reason, "(A,ES9.2)") &
            "expanded dS=", metric_deviation
         GOTO 100
      END IF

      dst_r(:, :) = REAL(stabilized, KIND=dp)
      dst_i(:, :) = AIMAG(stabilized)
      CALL cp_fm_set_submatrix(dst_real, dst_r)
      CALL cp_fm_set_submatrix(dst_imag, dst_i)
      success = .TRUE.

100   CONTINUE
      DEALLOCATE (source_coeff, s_coeff, h_coeff, h_projected, h_projected_work, metric_vectors, &
                  residual_block, ritz_vectors, s_projected, stabilized, coeff_work, metric_values, &
                  ritz_values, source_eigenvalues)
      DEALLOCATE (src_r, src_i, s_coeff_r, s_coeff_i, h_coeff_r, h_coeff_i, dst_r, dst_i)
      CALL cp_fm_release(s_src_real)
      CALL cp_fm_release(s_src_imag)
      CALL cp_fm_release(h_src_real)
      CALL cp_fm_release(h_src_imag)
      CALL cp_fm_release(tmp_fm)

   END SUBROUTINE ritz_reconstruct_wannier90_window

! **************************************************************************************************
!> \brief Map the full Wannier90 mesh to SCF representative k-points and symmetry operations.
!> \param kpoint full Wannier90 export k-point object
!> \param qs_kpoint SCF k-point object
!> \param source_kpoint source representative index for each full k-point
!> \param sym_index symmetry entry in source kp_sym; 0 direct, -1 time reversal only
!> \param success true if every full k-point was mapped
!> \param reason diagnostic message
! **************************************************************************************************
   SUBROUTINE build_wannier90_scf_mapping(kpoint, qs_kpoint, source_kpoint, sym_index, success, &
                                          reason)
      TYPE(kpoint_type), POINTER                         :: kpoint, qs_kpoint
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: source_kpoint, sym_index
      LOGICAL, INTENT(OUT)                               :: success
      CHARACTER(LEN=*), INTENT(OUT)                      :: reason

      INTEGER                                            :: ik, ikred, imatch, isym, nfull
      TYPE(kpoint_sym_type), POINTER                     :: kpsym

      success = .FALSE.
      reason = ""
      nfull = kpoint%nkp
      ALLOCATE (source_kpoint(nfull), sym_index(nfull))
      source_kpoint(:) = 0
      sym_index(:) = 0

      DO ikred = 1, qs_kpoint%nkp
         imatch = find_matching_kpoint(kpoint%xkp, qs_kpoint%xkp(1:3, ikred))
         IF (imatch > 0 .AND. source_kpoint(imatch) == 0) THEN
            source_kpoint(imatch) = ikred
            sym_index(imatch) = 0
         END IF
      END DO

      ! Prefer pure time-reversal partners before general atom/AO symmetry operations.
      DO ikred = 1, qs_kpoint%nkp
         imatch = find_matching_kpoint(kpoint%xkp, -qs_kpoint%xkp(1:3, ikred))
         IF (imatch > 0 .AND. source_kpoint(imatch) == 0) THEN
            source_kpoint(imatch) = ikred
            sym_index(imatch) = -1
         END IF
      END DO

      IF (ASSOCIATED(qs_kpoint%kp_sym)) THEN
         DO ikred = 1, qs_kpoint%nkp
            kpsym => qs_kpoint%kp_sym(ikred)%kpoint_sym
            IF (.NOT. ASSOCIATED(kpsym)) CYCLE
            IF (.NOT. kpsym%apply_symmetry) CYCLE
            DO isym = 1, kpsym%nwred
               imatch = find_matching_kpoint(kpoint%xkp, kpsym%xkp(1:3, isym))
               IF (imatch > 0 .AND. source_kpoint(imatch) == 0) THEN
                  source_kpoint(imatch) = ikred
                  sym_index(imatch) = isym
               END IF
            END DO
         END DO
      END IF

      DO ik = 1, nfull
         IF (source_kpoint(ik) == 0) THEN
            reason = "not all full-mesh k-points are represented by the SCF symmetry orbits"
            RETURN
         END IF
      END DO
      success = .TRUE.

   END SUBROUTINE build_wannier90_scf_mapping

! **************************************************************************************************
!> \brief Find a fractional k-point in a periodic mesh.
!> \param xkp_mesh mesh coordinates
!> \param xkp_search coordinate to find
!> \return matching index, or zero when no match is found
! **************************************************************************************************
   INTEGER FUNCTION find_matching_kpoint(xkp_mesh, xkp_search) RESULT(ik_match)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp_mesh
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp_search

      INTEGER                                            :: ik

      ik_match = 0
      DO ik = 1, SIZE(xkp_mesh, 2)
         IF (same_periodic_kpoint(xkp_mesh(1:3, ik), xkp_search)) THEN
            ik_match = ik
            RETURN
         END IF
      END DO

   END FUNCTION find_matching_kpoint

! **************************************************************************************************
!> \brief Compare two fractional k-points modulo reciprocal lattice vectors.
!> \param xkp_a first k-point
!> \param xkp_b second k-point
!> \return true if the k-points are equivalent
! **************************************************************************************************
   LOGICAL FUNCTION same_periodic_kpoint(xkp_a, xkp_b) RESULT(same)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp_a, xkp_b

      REAL(KIND=dp), DIMENSION(3)                        :: diff

      diff(1:3) = xkp_a(1:3) - xkp_b(1:3)
      diff(1:3) = diff(1:3) - REAL(NINT(diff(1:3)), KIND=dp)
      same = SUM(ABS(diff(1:3))) < 1.0e-8_dp

   END FUNCTION same_periodic_kpoint

! **************************************************************************************************
!> \brief Transform one SCF MO coefficient matrix to an equivalent full-mesh k-point.
!> \param src_real real part of source MO coefficients
!> \param src_imag imaginary part of source MO coefficients
!> \param dst_real real part of transformed MO coefficients
!> \param dst_imag imaginary part of transformed MO coefficients
!> \param qs_kpoint SCF k-point object containing symmetry operations
!> \param ikred representative k-point index
!> \param isym symmetry operation index; 0 direct copy, -1 time reversal
!> \param para_env global parallel environment
!> \param success true if the transformation was performed
!> \param reason diagnostic message
! **************************************************************************************************
   SUBROUTINE transform_wannier90_scf_mo(src_real, src_imag, dst_real, dst_imag, qs_kpoint, ikred, &
                                         isym, para_env, success, reason)
      TYPE(cp_fm_type), INTENT(IN)                       :: src_real, src_imag, dst_real, dst_imag
      TYPE(kpoint_type), POINTER                         :: qs_kpoint
      INTEGER, INTENT(IN)                                :: ikred, isym
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(OUT)                               :: success
      CHARACTER(LEN=*), INTENT(OUT)                      :: reason

      INTEGER :: iao, iatom, ikind, imo, irow, irow_source, irow_target, irow_work, nao, natom, &
         nmo, rot_slot, source_atom, source_dim, target_atom, target_dim
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_size, ao_start
      LOGICAL                                            :: reverse_phase, time_reversal
      REAL(KIND=dp)                                      :: arg, coeff_imag, coeff_real, coskl, &
                                                            phase_imag, phase_real, sinkl
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dst_i, dst_r, src_i, src_r
      REAL(KIND=dp), DIMENSION(3)                        :: xkp_phase
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rotmat
      TYPE(kind_rotmat_type), DIMENSION(:), POINTER      :: kind_rot
      TYPE(kpoint_sym_type), POINTER                     :: kpsym

      success = .FALSE.
      reason = ""
      IF (isym == 0) THEN
         CALL cp_fm_copy_general(src_real, dst_real, para_env)
         CALL cp_fm_copy_general(src_imag, dst_imag, para_env)
         success = .TRUE.
         RETURN
      END IF

      CALL cp_fm_get_info(src_real, nrow_global=nao, ncol_global=nmo)
      ALLOCATE (src_r(nao, nmo), src_i(nao, nmo), dst_r(nao, nmo), dst_i(nao, nmo))
      CALL cp_fm_get_submatrix(src_real, src_r)
      CALL cp_fm_get_submatrix(src_imag, src_i)
      dst_r(:, :) = 0.0_dp
      dst_i(:, :) = 0.0_dp

      IF (isym == -1) THEN
         dst_r(1:nao, 1:nmo) = src_r(1:nao, 1:nmo)
         dst_i(1:nao, 1:nmo) = -src_i(1:nao, 1:nmo)
         CALL cp_fm_set_submatrix(dst_real, dst_r)
         CALL cp_fm_set_submatrix(dst_imag, dst_i)
         DEALLOCATE (src_r, src_i, dst_r, dst_i)
         success = .TRUE.
         RETURN
      END IF

      IF (.NOT. ASSOCIATED(qs_kpoint%kp_sym)) THEN
         reason = "SCF symmetry operation data are not available"
         DEALLOCATE (src_r, src_i, dst_r, dst_i)
         RETURN
      END IF
      kpsym => qs_kpoint%kp_sym(ikred)%kpoint_sym
      IF (.NOT. ASSOCIATED(kpsym)) THEN
         reason = "SCF k-point symmetry operation is not available"
         DEALLOCATE (src_r, src_i, dst_r, dst_i)
         RETURN
      END IF
      IF (isym > kpsym%nwred) THEN
         reason = "SCF k-point symmetry operation index is out of range"
         DEALLOCATE (src_r, src_i, dst_r, dst_i)
         RETURN
      END IF
      IF (.NOT. ASSOCIATED(qs_kpoint%atype) .OR. .NOT. ASSOCIATED(qs_kpoint%kind_rotmat)) THEN
         reason = "SCF atom mappings or basis rotation matrices are not available"
         DEALLOCATE (src_r, src_i, dst_r, dst_i)
         RETURN
      END IF

      rot_slot = find_wannier90_rotation_slot(qs_kpoint, kpsym%rotp(isym))
      IF (rot_slot == 0) THEN
         reason = "could not match the SCF symmetry operation to a basis rotation"
         DEALLOCATE (src_r, src_i, dst_r, dst_i)
         RETURN
      END IF
      kind_rot => qs_kpoint%kind_rotmat(rot_slot, :)
      natom = SIZE(qs_kpoint%atype)
      ALLOCATE (ao_start(natom), ao_size(natom))
      irow = 1
      DO iatom = 1, natom
         ikind = qs_kpoint%atype(iatom)
         IF (.NOT. ASSOCIATED(kind_rot(ikind)%rmat)) THEN
            reason = "a required basis rotation matrix is not available"
            DEALLOCATE (src_r, src_i, dst_r, dst_i, ao_start, ao_size)
            RETURN
         END IF
         ao_start(iatom) = irow
         ao_size(iatom) = SIZE(kind_rot(ikind)%rmat, 2)
         irow = irow + ao_size(iatom)
      END DO
      IF (irow - 1 /= nao) THEN
         reason = "atom-resolved AO dimensions do not match the MO coefficient matrix"
         DEALLOCATE (src_r, src_i, dst_r, dst_i, ao_start, ao_size)
         RETURN
      END IF

      time_reversal = kpsym%rotp(isym) < 0
      reverse_phase = qs_kpoint%gamma_centered .AND. ANY(MOD(qs_kpoint%nkp_grid, 2) == 0)
      xkp_phase(1:3) = kpsym%xkp(1:3, isym)
      DO iatom = 1, natom
         source_atom = iatom
         target_atom = kpsym%f0(iatom, isym)
         ikind = qs_kpoint%atype(source_atom)
         rotmat => kind_rot(ikind)%rmat
         source_dim = ao_size(source_atom)
         target_dim = ao_size(target_atom)
         IF (SIZE(rotmat, 1) /= target_dim .OR. SIZE(rotmat, 2) /= source_dim) THEN
            reason = "basis rotation dimensions do not match the atom/AO symmetry transform"
            DEALLOCATE (src_r, src_i, dst_r, dst_i, ao_start, ao_size)
            RETURN
         END IF
         arg = REAL(kpsym%fcell(1, source_atom, isym), KIND=dp)*xkp_phase(1) + &
               REAL(kpsym%fcell(2, source_atom, isym), KIND=dp)*xkp_phase(2) + &
               REAL(kpsym%fcell(3, source_atom, isym), KIND=dp)*xkp_phase(3)
         IF (ASSOCIATED(kpsym%kgphase)) THEN
            arg = arg + kpsym%kgphase(source_atom, isym)
         END IF
         IF (reverse_phase) arg = -arg
         coskl = COS(twopi*arg)
         sinkl = SIN(twopi*arg)
         DO imo = 1, nmo
            DO irow_work = 1, source_dim
               irow_source = ao_start(source_atom) + irow_work - 1
               coeff_real = src_r(irow_source, imo)
               coeff_imag = src_i(irow_source, imo)
               IF (time_reversal) coeff_imag = -coeff_imag
               phase_real = coskl*coeff_real - sinkl*coeff_imag
               phase_imag = sinkl*coeff_real + coskl*coeff_imag
               DO iao = 1, target_dim
                  irow_target = ao_start(target_atom) + iao - 1
                  dst_r(irow_target, imo) = dst_r(irow_target, imo) + &
                                            rotmat(iao, irow_work)*phase_real
                  dst_i(irow_target, imo) = dst_i(irow_target, imo) + &
                                            rotmat(iao, irow_work)*phase_imag
               END DO
            END DO
         END DO
      END DO

      CALL cp_fm_set_submatrix(dst_real, dst_r)
      CALL cp_fm_set_submatrix(dst_imag, dst_i)
      DEALLOCATE (src_r, src_i, dst_r, dst_i, ao_start, ao_size)
      success = .TRUE.

   END SUBROUTINE transform_wannier90_scf_mo

! **************************************************************************************************
!> \brief Locate the basis-rotation slot corresponding to a signed k-point symmetry operation.
!> \param qs_kpoint SCF k-point object
!> \param rotp signed operation identifier
!> \return rotation slot, or zero when no slot matches
! **************************************************************************************************
   INTEGER FUNCTION find_wannier90_rotation_slot(qs_kpoint, rotp) RESULT(rot_slot)
      TYPE(kpoint_type), POINTER                         :: qs_kpoint
      INTEGER, INTENT(IN)                                :: rotp

      INTEGER                                            :: irot, rot_abs

      rot_slot = 0
      rot_abs = ABS(rotp)
      IF (.NOT. ASSOCIATED(qs_kpoint%ibrot)) RETURN
      DO irot = 1, SIZE(qs_kpoint%ibrot)
         IF (rot_abs == qs_kpoint%ibrot(irot)) THEN
            rot_slot = irot
            RETURN
         END IF
      END DO

   END FUNCTION find_wannier90_rotation_slot

! **************************************************************************************************
!> \brief Infer a tensor-product Wannier90 mesh from explicit fractional k-point coordinates.
!> \param kpt_latt explicit k-point coordinates in reciprocal-lattice units
!> \param mp_grid inferred mesh dimensions
!> \param valid true if the coordinate set is compatible with a tensor-product mesh
! **************************************************************************************************
   SUBROUTINE infer_wannier_mp_grid(kpt_latt, mp_grid, valid)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: kpt_latt
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: mp_grid
      LOGICAL, INTENT(OUT)                               :: valid

      INTEGER                                            :: coord_id, i, idim, idx, n_unique, &
                                                            num_kpts, stride, unique_id
      LOGICAL                                            :: known
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: seen
      REAL(KIND=dp)                                      :: coord
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: unique_coord

      num_kpts = SIZE(kpt_latt, 2)
      mp_grid(:) = 0
      ALLOCATE (unique_coord(3, num_kpts))
      DO idim = 1, 3
         n_unique = 0
         DO i = 1, num_kpts
            coord = kpt_latt(idim, i) - FLOOR(kpt_latt(idim, i))
            IF (ABS(coord - 1.0_dp) < 1.0e-8_dp) coord = 0.0_dp
            known = .FALSE.
            DO unique_id = 1, n_unique
               IF (ABS(unique_coord(idim, unique_id) - coord) < 1.0e-8_dp) THEN
                  known = .TRUE.
                  EXIT
               END IF
            END DO
            IF (.NOT. known) THEN
               n_unique = n_unique + 1
               unique_coord(idim, n_unique) = coord
            END IF
         END DO
         mp_grid(idim) = n_unique
      END DO
      valid = (mp_grid(1)*mp_grid(2)*mp_grid(3) == num_kpts)
      IF (valid) THEN
         ALLOCATE (seen(num_kpts))
         seen(:) = .FALSE.
         DO i = 1, num_kpts
            idx = 1
            stride = 1
            DO idim = 1, 3
               coord = kpt_latt(idim, i) - FLOOR(kpt_latt(idim, i))
               IF (ABS(coord - 1.0_dp) < 1.0e-8_dp) coord = 0.0_dp
               coord_id = 0
               DO unique_id = 1, mp_grid(idim)
                  IF (ABS(unique_coord(idim, unique_id) - coord) < 1.0e-8_dp) THEN
                     coord_id = unique_id
                     EXIT
                  END IF
               END DO
               CPASSERT(coord_id > 0)
               idx = idx + (coord_id - 1)*stride
               stride = stride*mp_grid(idim)
            END DO
            IF (seen(idx)) valid = .FALSE.
            seen(idx) = .TRUE.
         END DO
         valid = valid .AND. ALL(seen)
         DEALLOCATE (seen)
      END IF
      DEALLOCATE (unique_coord)

   END SUBROUTINE infer_wannier_mp_grid

END MODULE qs_wannier90
