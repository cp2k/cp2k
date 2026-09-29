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
   USE bibliography,                    ONLY: Gresch2017,&
                                              Soluyanov2011,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: &
        cp_cfm_copy_general, cp_cfm_create, cp_cfm_get_element, cp_cfm_get_info, &
        cp_cfm_get_submatrix, cp_cfm_release, cp_cfm_set_submatrix, cp_cfm_to_cfm, cp_cfm_type, &
        cp_fm_to_cfm
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_desymmetrize, dbcsr_p_type, &
        dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, &
        dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE kpoint_methods,                  ONLY: kpoint_env_initialize,&
                                              kpoint_init_cell_index,&
                                              kpoint_initialize,&
                                              kpoint_initialize_mo_set,&
                                              kpoint_initialize_mos,&
                                              rskp_transform
   USE kpoint_mo_symmetry_methods,      ONLY: kpoint_same_periodic,&
                                              kpoint_transform_scf_mo
   USE kpoint_types,                    ONLY: get_kpoint_info,&
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
   USE parallel_gemm_api,               ONLY: parallel_gemm
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
   USE soc_pseudopotential_methods,     ONLY: V_SOC_xyz_from_pseudopotential
   USE topology_inversion,              ONLY: gaussian_inversion_action
   USE topology_state_io,               ONLY: topology_state_begin,&
                                              topology_state_point
   USE topology_symmetry,               ONLY: inversion_representation
   USE topology_tqc,                    ONLY: inversion_ebr_signature,&
                                              inversion_indicators
   USE topology_wilson,                 ONLY: chern_from_wcc,&
                                              surface_resolved,&
                                              wcc_distance,&
                                              wilson_spectrum,&
                                              wilson_step,&
                                              z2_from_wcc
   USE wannier90,                       ONLY: wannier_setup
   USE wannier90_nnkpts,                ONLY: read_wannier90_nnkpts
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_wannier90'
   INTEGER, PARAMETER, PRIVATE :: w90_kpoints_mp_grid = 0, &
                                  w90_kpoints_scf = 1, w90_kpoints_nnkp = 2, w90_kpoints_wilson = 3, &
                                  w90_kpoints_trim = 4

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

      INTEGER                                            :: handle, i, ichern, iw, iz2, &
                                                            max_refinement, old_chern, old_z2, &
                                                            refinement, source
      LOGICAL                                            :: converged, explicit, require_chern, &
                                                            require_z2
      REAL(KIND=dp)                                      :: error, movement, tolerance
      REAL(KIND=dp), ALLOCATABLE                         :: centres(:, :), previous(:, :)
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

         CALL section_vals_val_get(w_input, "KPOINTS_SOURCE", i_val=source)
         CALL section_vals_val_get(w_input, "WILSON_MAX_REFINEMENT", i_val=max_refinement)
         CALL section_vals_val_get(w_input, "WILSON_TOL", r_val=tolerance)
         CALL section_vals_val_get(w_input, "Z2", l_val=require_z2)
         CALL section_vals_val_get(w_input, "CHERN", l_val=require_chern)
         IF (tolerance <= 0.0_dp) CPABORT("WILSON_TOL must be positive.")
         IF (source == w90_kpoints_wilson) THEN
            IF (max_refinement < 1 .OR. max_refinement > 10) THEN
               CPABORT("WILSON_MAX_REFINEMENT must be between 1 and 10.")
            END IF
            converged = .FALSE.
            old_z2 = -1
            old_chern = HUGE(0)
            DO refinement = 0, max_refinement
               CALL wannier90_files(qs_env, w_input, iw, refinement, centres, iz2, ichern)
               IF (refinement > 0) THEN
                  error = 0.0_dp
                  DO i = 1, SIZE(previous, 2)
                     error = MAX(error, wcc_distance(previous(:, i), centres(:, 2*i - 1)))
                  END DO
                  movement = 0.0_dp
                  DO i = 2, SIZE(centres, 2)
                     movement = MAX(movement, wcc_distance(centres(:, i - 1), centres(:, i)))
                  END DO
                  IF (iw > 0) WRITE (iw, '(T2,A,I0,A,ES12.4,A,ES12.4)') &
                     "TOPOLOGY| Refinement ", refinement, ": WCC change ", error, ", transverse step ", movement
                  converged = error < tolerance .AND. movement < 0.1_dp .AND. iz2 == old_z2
                  IF (require_z2) converged = converged .AND. iz2 >= 0 .AND. surface_resolved(centres)
                  IF (require_chern) converged = converged .AND. ichern /= HUGE(0) .AND. ichern == old_chern
                  IF (converged) EXIT
               END IF
               CALL MOVE_ALLOC(centres, previous)
               old_z2 = iz2
               old_chern = ichern
            END DO
            IF (.NOT. converged) THEN
               CPABORT("Wilson surface unconverged; increase WILSON_MAX_REFINEMENT or mesh.")
            END IF
            IF (iw > 0) THEN
               WRITE (iw, '(T2,A)') "TOPOLOGY| Wilson surface sampling converged."
               IF (iz2 >= 0) WRITE (iw, '(T2,A,I0)') "TOPOLOGY| Converged Z2 invariant: ", iz2
               IF (require_chern) WRITE (iw, '(T2,A,I0)') "TOPOLOGY| Converged first Chern number: ", ichern
            END IF
         ELSE
            CALL wannier90_files(qs_env, w_input, iw, 0, centres, iz2, ichern)
         END IF

         IF (iw > 0) THEN
            WRITE (iw, '(/,T2,A)') &
               '!--------------------------------End of Wannier90-----------------------------!'
         END IF
      END IF
      CALL timestop(handle)

   END SUBROUTINE wannier90_interface

! **************************************************************************************************
!> \brief Copies a k-point CFM into the all-k-point Wannier work distribution.
!> \param source source CFM; it may be unassociated on ranks outside its k-point group
!> \param destination destination CFM
!> \param para_env communicator spanning the source and destination distributions
! **************************************************************************************************
   SUBROUTINE wannier90_copy_cfm_general(source, destination, para_env)
      TYPE(cp_cfm_type), INTENT(IN), POINTER             :: source
      TYPE(cp_cfm_type), INTENT(IN)                      :: destination
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env

      TYPE(cp_cfm_type)                                  :: empty_source

      IF (ASSOCIATED(source)) THEN
         CALL cp_cfm_copy_general(source, destination, para_env)
      ELSE
         CALL cp_cfm_copy_general(empty_source, destination, para_env)
      END IF

   END SUBROUTINE wannier90_copy_cfm_general

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param input ...
!> \param iw ...
!> \param refinement number of joint mesh doublings
!> \param wcc_out Wilson centres for each closed loop
!> \param z2_value candidate Z2 parity, or -1 when not requested
!> \param chern_value candidate first Chern number, or HUGE(0) when unavailable
! **************************************************************************************************
   SUBROUTINE wannier90_files(qs_env, input, iw, refinement, wcc_out, z2_value, chern_value)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input
      INTEGER, INTENT(IN)                                :: iw, refinement
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: wcc_out(:, :)
      INTEGER, INTENT(OUT)                               :: z2_value, chern_value

      INTEGER, PARAMETER                                 :: num_nnmax = 12

      CHARACTER(len=2)                                   :: asym
      CHARACTER(len=20), ALLOCATABLE, DIMENSION(:)       :: atom_symbols
      CHARACTER(len=default_path_length)                 :: nnkp_file
      CHARACTER(len=default_string_length)               :: filename, input_kp_scheme, reuse_reason, &
                                                            seed_name
      CHARACTER(LEN=timestamp_length)                    :: timestamp
      COMPLEX(KIND=dp)                                   :: zmmn
      COMPLEX(KIND=dp), ALLOCATABLE :: export_coeff(:, :), export_scalar(:, :), link_matrix(:, :), &
         parity_metric(:, :), parity_phase(:), scalar_overlap(:, :), soc_h(:, :), soc_u(:, :), &
         soc_xyz(:, :, :), spinor_coeff(:, :, :), wilson_product(:, :, :)
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: reference_mo
      INTEGER :: aligned_degenerate_blocks, aligned_degenerate_max_size, axis, base_mesh(2), &
         counts(2), first_point, i, i_rep, ib, ib1, ib2, ibs, ik, ik2, ikk, ikpgr, iloop, ipoint, &
         ispin, iunit, ix, iy, iz, jpar, k, kpoints_source, loop_direction(3), n_rep, nadd, nao, &
         nberry_images, nbs, nelectron, nexcl, nkp, nloop, nmo, nntot, npoint, nscalar, nspins, &
         num_atoms, num_bands, num_bands_tot, num_kpts, num_wann, spin_channel, state_components, &
         state_unit, status, strong, tqc_dimension, trim_id, weak(3), z4
      INTEGER, ALLOCATABLE                               :: ebr_coefficients(:, :), loop_index(:), &
                                                            parity_map(:), parity_odd(:)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: band_map, exclude_bands
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: nblist, nnlist
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: nncell
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(3)                              :: input_nkp_grid, mp_grid
      INTEGER, DIMENSION(:), POINTER                     :: invals
      INTEGER, DIMENSION(:, :, :), POINTER               :: berry_cell_index, cell_to_index
      LOGICAL :: diis_step, do_chern, do_kpoints, do_parity, do_soc, do_tqc, do_wilson, do_z2, &
         export_state, full_mesh_diagonalized, gamma_only, input_full_grid, input_gamma_centered, &
         input_kpoint_symmetry, lowest_bands, mp_grid_explicit, mp_grid_valid, my_kpgrp, mygrp, &
         nonnegative_atomic, ordered_berry, require_global_gap, reuse_scf_mos, reused_scf_mos, &
         signed_atomic, spinors, time_reversal, use_bloch_phases, validate_reuse_ok, &
         validate_reuse_scf_mos
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: keep_band
      REAL(KIND=dp) :: aligned_degenerate_min_svalue, berry_phase, chern_winding, cmmn, &
         conduction_min, direct_gap, gap_tol, gauge_arg, gauge_imag, gauge_real, gauge_tmp, ksign, &
         link_sv, loop_origin(3), pair_gap, parity_checks(4), parity_energy_error, &
         parity_energy_tolerance, parity_error, parity_gap, parity_origin(3), parity_tolerance, &
         reuse_candidate_deviation, reuse_candidate_metric_deviation, reuse_candidate_min_svalue, &
         reuse_candidate_residual, rmmn, transverse(3), valence_max, &
         validation_eigenvalue_deviation, validation_min_svalue, validation_subspace_deviation, &
         wkp_ref
      REAL(KIND=dp), ALLOCATABLE                         :: loop_sv(:), parity_energies(:), &
                                                            scalar_values(:), spinor_values(:, :)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atoms_cart, b_latt, kpt_latt
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: reference_eigenvalues
      REAL(KIND=dp), DIMENSION(3)                        :: bvec, input_kp_shift, phase_center
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_inv, real_lattice, recip_lattice
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, wkp, wkp_source
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp, xkp_source
      REAL(KIND=dp), POINTER                             :: rvals(:)
      TYPE(berry_matrix_type), DIMENSION(:), POINTER     :: berry_matrix
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: fmk1_cfm, fmk2_cfm, mmn_cfm, omat_cfm, &
                                                            tmp_cfm
      TYPE(cp_cfm_type), POINTER                         :: cmo_coeff_source
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_ao, matrix_struct_mmn, &
                                                            matrix_struct_work
      TYPE(cp_fm_type)                                   :: mat_imag, mat_real
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s, soc_matrices
      TYPE(dbcsr_type), POINTER                          :: cmatrix, cmatrix_full, loop_imag, &
                                                            loop_real, rmatrix, rmatrix_full
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: berry_kpoint, kpoint, qs_kpoint
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: overlap_nl, sab_nl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_environment_type), POINTER                 :: qs_env_kp
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      !--------------------------------------------------------------------------------------------!

      ! generate all arrays needed for the setup call
      CALL section_vals_val_get(input, "SEED_NAME", c_val=seed_name)
      CALL section_vals_val_get(input, "MP_GRID", i_vals=invals, explicit=mp_grid_explicit)
      CALL section_vals_val_get(input, "KPOINTS_SOURCE", i_val=kpoints_source)
      ordered_berry = kpoints_source >= w90_kpoints_nnkp
      CALL section_vals_val_get(input, "NNKP_FILE", c_val=nnkp_file)
      CALL section_vals_val_get(input, "SPIN_CHANNEL", i_val=spin_channel)
      CALL section_vals_val_get(input, "WILSON_LOOP", l_val=do_wilson)
      CALL section_vals_val_get(input, "Z2", l_val=do_z2)
      CALL section_vals_val_get(input, "CHERN", l_val=do_chern)
      CALL section_vals_val_get(input, "REQUIRE_GLOBAL_GAP", l_val=require_global_gap)
      CALL section_vals_val_get(input, "STATE_EXPORT", l_val=export_state)
      CALL section_vals_val_get(input, "PARITY", l_val=do_parity)
      CALL section_vals_val_get(input, "INVERSION_TQC", l_val=do_tqc)
      CALL section_vals_val_get(input, "TQC_DIMENSION", i_val=tqc_dimension)
      CALL section_vals_val_get(input, "PARITY_ORIGIN", r_vals=rvals)
      parity_origin = rvals
      CALL section_vals_val_get(input, "PARITY_TOLERANCE", r_val=parity_tolerance)
      CALL section_vals_val_get(input, "PARITY_ENERGY_TOL", r_val=parity_energy_tolerance)
      do_parity = do_parity .OR. do_tqc .OR. kpoints_source == w90_kpoints_trim
      IF (do_parity .AND. .NOT. ordered_berry) THEN
         CPABORT("PARITY requires explicit NNKP, WILSON or TRIM points.")
      END IF
      IF (tqc_dimension < 2 .OR. tqc_dimension > 3) CPABORT("TQC_DIMENSION must be 2 or 3.")
      IF (parity_tolerance <= 0.0_dp .OR. parity_energy_tolerance <= 0.0_dp) THEN
         CPABORT("Invalid parity tolerance.")
      END IF
      IF (export_state .AND. .NOT. ordered_berry) THEN
         CPABORT("STATE_EXPORT requires explicit NNKP, WILSON or TRIM points.")
      END IF
      CALL section_vals_val_get(input, "TIME_REVERSAL", l_val=time_reversal)
      CALL section_vals_val_get(input, "SOC", l_val=do_soc)
      IF (do_tqc .AND. (.NOT. do_soc .OR. .NOT. time_reversal)) THEN
         CPABORT("INVERSION_TQC requires SOC T and TIME_REVERSAL T.")
      END IF
      CALL section_vals_val_get(input, "WILSON_GAP_TOL", r_val=gap_tol)
      IF (gap_tol <= 0.0_dp) CPABORT("WILSON_GAP_TOL must be positive.")
      z2_value = -1
      chern_value = HUGE(0)
      do_wilson = do_wilson .OR. kpoints_source == w90_kpoints_wilson
      IF (do_wilson .AND. kpoints_source == w90_kpoints_trim) THEN
         CPABORT("TRIM points are not Wilson loops; use KPOINTS_SOURCE WILSON.")
      END IF
      IF (do_wilson) CALL cite_reference(Gresch2017)
      IF (do_z2) CALL cite_reference(Soluyanov2011)
      IF ((do_wilson .OR. do_soc) .AND. kpoints_source < w90_kpoints_nnkp) THEN
         CPABORT("Native Wilson/SOC requires KPOINTS_SOURCE NNKP or WILSON.")
      END IF
      IF (do_z2 .AND. (kpoints_source /= w90_kpoints_wilson .OR. .NOT. do_soc .OR. .NOT. time_reversal)) THEN
         CPABORT("Z2 requires KPOINTS_SOURCE WILSON, SOC T, and TIME_REVERSAL T.")
      END IF
      IF (do_chern .AND. (kpoints_source /= w90_kpoints_wilson .OR. do_z2)) THEN
         CPABORT("CHERN requires a full WILSON surface and cannot be combined with Z2.")
      END IF
      IF (require_global_gap .AND. .NOT. do_wilson) THEN
         CPABORT("REQUIRE_GLOBAL_GAP requires Wilson analysis.")
      END IF
      CALL get_qs_env(qs_env, dft_control=dft_control)
      IF (spin_channel < 1 .OR. spin_channel > dft_control%nspins) THEN
         CPABORT("WANNIER90%SPIN_CHANNEL is not available in this calculation.")
      END IF
      IF (do_soc .AND. dft_control%nspins /= 1) THEN
         CPABORT("WANNIER90 SOC currently requires a restricted SCF.")
      END IF
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
      NULLIFY (cmo_coeff_source, kpoint, qs_kpoint, xkp, wkp, xkp_source, wkp_source)
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
      CASE (w90_kpoints_nnkp, w90_kpoints_wilson, w90_kpoints_trim)
         IF (kpoints_source == w90_kpoints_nnkp) THEN
            CALL read_wannier90_nnkpts(nnkp_file, real_lattice, recip_lattice, kpt_latt, nnlist, nncell)
         ELSE IF (kpoints_source == w90_kpoints_trim) THEN
            num_kpts = 2**tqc_dimension
            ALLOCATE (kpt_latt(3, num_kpts), nnlist(num_kpts, 1), nncell(3, num_kpts, 1))
            kpt_latt = 0.0_dp
            nncell = 0
            DO i = 1, num_kpts
               DO axis = 1, tqc_dimension
                  IF (BTEST(i - 1, axis - 1)) kpt_latt(axis, i) = 0.5_dp
               END DO
               nnlist(i, 1) = i
            END DO
         ELSE
            CALL section_vals_val_get(input, "WILSON_MESH", i_vals=invals)
            base_mesh = invals
            IF (base_mesh(1) < 2 .OR. base_mesh(2) < 2) THEN
               CPABORT("WILSON_MESH entries must be at least two.")
            END IF
            npoint = base_mesh(1)*2**refinement
            nloop = (base_mesh(2) - 1)*2**refinement + 1
            CALL section_vals_val_get(input, "WILSON_DIRECTION", i_vals=invals)
            loop_direction = invals
            IF (ALL(loop_direction == 0)) CPABORT("WILSON_DIRECTION must be nonzero.")
            CALL section_vals_val_get(input, "WILSON_ORIGIN", r_vals=rvals)
            loop_origin = rvals
            CALL section_vals_val_get(input, "WILSON_TRANSVERSE", r_vals=rvals)
            transverse = rvals
            bvec = [loop_direction(2)*transverse(3) - loop_direction(3)*transverse(2), &
                    loop_direction(3)*transverse(1) - loop_direction(1)*transverse(3), &
                    loop_direction(1)*transverse(2) - loop_direction(2)*transverse(1)]
            IF (SUM(bvec**2) < 1.e-20_dp) THEN
               CPABORT("Wilson loop and transverse vectors must be linearly independent.")
            END IF
            IF (do_chern .AND. MAXVAL(ABS(transverse - NINT(transverse))) > 1.e-10_dp) THEN
               CPABORT("CHERN requires integer WILSON_TRANSVERSE: a closed full surface, not a half-plane.")
            END IF
            IF (do_z2) THEN
               IF (MAXVAL(ABS(2*loop_origin - NINT(2*loop_origin))) > 1.e-10_dp .OR. &
                   MAXVAL(ABS(2*transverse - NINT(2*transverse))) > 1.e-10_dp) THEN
                  CPABORT("Z2 surface boundaries must pass through time-reversal-invariant momenta.")
               END IF
               ! Initially restrict native Z2 to standard half-planes, avoiding multiple coverings.
               IF (SUM(ABS(loop_direction)) /= 1 .OR. &
                   ABS(SUM(ABS(transverse)) - 0.5_dp) > 1.e-10_dp .OR. &
                   ABS(DOT_PRODUCT(REAL(loop_direction, dp), transverse)) > 1.e-10_dp) THEN
                  CPABORT("Native Z2 requires distinct coordinate axes with windings one and one half.")
               END IF
            END IF
            num_kpts = npoint*nloop
            ALLOCATE (kpt_latt(3, num_kpts), nnlist(num_kpts, 1), nncell(3, num_kpts, 1))
            nncell = 0
            DO iloop = 1, nloop
               first_point = (iloop - 1)*npoint + 1
               DO ipoint = 1, npoint
                  i = first_point + ipoint - 1
                  kpt_latt(:, i) = loop_origin + REAL(iloop - 1, dp)/REAL(nloop - 1, dp)*transverse + &
                                   REAL(ipoint - 1, dp)/REAL(npoint, dp)*loop_direction
                  nnlist(i, 1) = i + 1
               END DO
               nnlist(i, 1) = first_point
               nncell(:, i, 1) = loop_direction
            END DO
         END IF
         num_kpts = SIZE(kpt_latt, 2)
         nntot = SIZE(nnlist, 2)
         kpoint%kp_scheme = "GENERAL"
         kpoint%symmetry = .FALSE.
         kpoint%verbose = .FALSE.
         kpoint%full_grid = .TRUE.
         kpoint%eps_geo = 1.0e-6_dp
         kpoint%use_real_wfn = .FALSE.
         kpoint%parallel_group_size = para_env%num_pe
         kpoint%nkp = num_kpts
         ALLOCATE (kpoint%xkp(3, num_kpts), kpoint%wkp(num_kpts))
         kpoint%xkp = kpt_latt
         ! Export weights do not enter the converged SCF density.
         kpoint%wkp = 1.0_dp/REAL(num_kpts, KIND=dp)
         IF (iw > 0) WRITE (iw, '(T2,A,I0,A,I0)') &
            "WANNIER90| Explicit points: ", num_kpts, ", neighbours per point: ", nntot
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
      CALL get_mo_set(mo_set=mos(1), nao=nao, nmo=num_bands_tot, nelectron=nelectron)
      num_bands_tot = MIN(nao, num_bands_tot + nadd)
      nscalar = num_bands_tot
      IF (do_soc) num_bands_tot = 2*nscalar
      ALLOCATE (keep_band(num_bands_tot))
      keep_band = .TRUE.
      DO i = 1, nexcl
         ib = exclude_bands(i)
         IF (ib < 1 .OR. ib > num_bands_tot) CPABORT("WANNIER90%EXCLUDE_BANDS: index out of range.")
         IF (.NOT. keep_band(ib)) CPABORT("WANNIER90%EXCLUDE_BANDS: duplicate band index.")
         keep_band(ib) = .FALSE.
      END DO
      num_bands = COUNT(keep_band)
      IF (num_bands == 0) CPABORT("WANNIER90: no bands left after EXCLUDE_BANDS.")
      ALLOCATE (band_map(num_bands))
      band_map(:) = PACK([(i, i=1, num_bands_tot)], keep_band)
      DEALLOCATE (keep_band)
      lowest_bands = num_bands < num_bands_tot .AND. ALL(band_map == [(i, i=1, num_bands)])
      IF (require_global_gap .AND. .NOT. lowest_bands) THEN
         CPABORT("REQUIRE_GLOBAL_GAP needs a lowest-band prefix and an excluded band above it.")
      END IF
      IF (do_z2) THEN
         IF (num_bands /= nelectron) THEN
            CPABORT("Z2 requires one occupied spinor per electron; adjust EXCLUDE_BANDS.")
         END IF
         IF (num_bands == num_bands_tot .OR. MOD(num_bands, 2) /= 0) THEN
            CPABORT("Z2 needs an even occupied spinor subspace and at least one excluded conduction band.")
         END IF
         IF (ANY(band_map /= [(i, i=1, num_bands)])) THEN
            CPABORT("Z2 requires the lowest occupied spinor bands.")
         END IF
      END IF
      IF (do_wilson) THEN
         IF (nntot /= 1) CPABORT("WILSON_LOOP requires exactly one directed neighbour per point.")
         ALLOCATE (loop_index(num_kpts))
         i = 1
         nloop = 0
         DO WHILE (i <= num_kpts)
            first_point = i
            nloop = nloop + 1
            DO
               loop_index(i) = nloop
               IF (nnlist(i, 1) == first_point) EXIT
               IF (nnlist(i, 1) /= i + 1 .OR. i == num_kpts) THEN
                  CPABORT("WILSON_LOOP requires contiguous ordered closed loops.")
               END IF
               i = i + 1
            END DO
            i = i + 1
         END DO
         ALLOCATE (wilson_product(num_bands, num_bands, nloop), loop_sv(nloop))
         wilson_product = CMPLX(0.0_dp, 0.0_dp, dp)
         DO i = 1, num_bands
            wilson_product(i, i, :) = CMPLX(1.0_dp, 0.0_dp, dp)
         END DO
         loop_sv = 1.0_dp
         ALLOCATE (wcc_out(num_bands, nloop))
      ELSE
         ALLOCATE (wcc_out(0, 0))
      END IF
      ALLOCATE (link_matrix(num_bands, num_bands))
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
      IF (kpoints_source < w90_kpoints_nnkp) THEN
         ALLOCATE (nnlist(num_kpts, num_nnmax))
         ALLOCATE (nncell(3, num_kpts, num_nnmax))
         nnlist(:, :) = 0
         nncell(:, :, :) = 0
         nntot = 0
         IF (iw > 0) THEN
            CALL wannier_setup(mp_grid, num_kpts, real_lattice, recip_lattice, &
                               kpt_latt, nntot, nnlist, nncell, iw)
         END IF
         CALL para_env%sum(nntot)
         CALL para_env%sum(nnlist)
         CALL para_env%sum(nncell)
      END IF

      CALL get_qs_env(qs_env, para_env=para_env)

      IF (para_env%is_source() .AND. kpoints_source < w90_kpoints_nnkp) THEN
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
            CALL save_wannier90_mo_snapshot(kpoint, nspins, para_env, reference_mo, &
                                            reference_eigenvalues)
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
                                                  para_env, reference_mo, &
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
               CALL restore_wannier90_mo_snapshot(kpoint, reference_mo, &
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
      IF (ALLOCATED(reference_mo)) DEALLOCATE (reference_mo)
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
      overlap_nl => sab_nl
      NULLIFY (berry_kpoint)
      IF (ordered_berry) THEN
         CALL get_qs_env(qs_env_kp, sab_all=overlap_nl)
         CALL kpoint_create(berry_kpoint)
         CALL kpoint_init_cell_index(berry_kpoint, overlap_nl, para_env, nberry_images)
         CALL get_kpoint_info(berry_kpoint, cell_to_index=berry_cell_index)
      END IF
      IF (.NOT. ASSOCIATED(overlap_nl)) CPABORT("Explicit overlaps require k-point neighbour lists.")
      ALLOCATE (berry_matrix(nbs))
      DO i = 1, nbs
         NULLIFY (berry_matrix(i)%cosmat)
         NULLIFY (berry_matrix(i)%sinmat)
         bvec(1:3) = twopi*MATMUL(TRANSPOSE(cell%h_inv(1:3, 1:3)), b_latt(1:3, i))
         IF (ordered_berry) bvec = -bvec
         CALL build_berry_kpoint_matrix(qs_env_kp, berry_matrix(i)%cosmat, &
                                        berry_matrix(i)%sinmat, bvec, ordered=ordered_berry, ordered_kpoints=berry_kpoint)
      END DO
      ! work matrices for MOs (all group)
      kp => kpoint%kp_env(1)%kpoint_env
      CALL get_mo_set(kp%mos(1), nmo=nmo)
      IF (nmo /= nscalar) CPABORT("WANNIER90: unexpected orbital count in export.")
      NULLIFY (matrix_struct_ao, matrix_struct_work)
      CALL cp_fm_struct_create(matrix_struct_work, nrow_global=nao, &
                               ncol_global=nmo, &
                               para_env=para_env, &
                               context=blacs_env)
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
      IF (ordered_berry) THEN
         ALLOCATE (loop_real, loop_imag)
         CALL dbcsr_create(loop_real, template=matrix_s(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(loop_imag, template=matrix_s(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_alloc_block_from_nbl(loop_real, overlap_nl)
         CALL cp_dbcsr_alloc_block_from_nbl(loop_imag, overlap_nl)
         CALL cp_dbcsr_alloc_block_from_nbl(rmatrix_full, overlap_nl)
         CALL cp_dbcsr_alloc_block_from_nbl(cmatrix_full, overlap_nl)
      END IF
      !
      CALL get_kpoint_info(kpoint=kpoint, cell_to_index=cell_to_index)
      nspins = dft_control%nspins
      IF (do_soc) THEN
         ! Second variation in the scalar KS eigenbasis. Retain selected spinors at each k.
         CPASSERT(ALL(kpoint%kp_range == [1, num_kpts]))
         NULLIFY (soc_matrices)
         CALL V_SOC_xyz_from_pseudopotential(qs_env_kp, soc_matrices)
         ALLOCATE (soc_xyz(nmo, nmo, 3), soc_h(2*nmo, 2*nmo), soc_u(2*nmo, 2*nmo))
         ALLOCATE (scalar_values(2*nmo), spinor_values(2*nmo, num_kpts))
         ALLOCATE (spinor_coeff(2*nmo, num_bands, num_kpts), scalar_overlap(nmo, nmo))
         DO ik = 1, num_kpts
            kp => kpoint%kp_env(ik)%kpoint_env
            CALL cp_cfm_copy_general(kp%mos(1)%cmo_coeff, fmk1_cfm, para_env)
            CALL get_mo_set(kp%mos(1), eigenvalues=eigenvalues)
            scalar_values(1:nmo) = eigenvalues(1:nmo)
            scalar_values(nmo + 1:) = eigenvalues(1:nmo)
            DO axis = 1, 3
               CALL dbcsr_set(rmatrix, 0.0_dp)
               CALL dbcsr_set(cmatrix, 0.0_dp)
               CALL rskp_transform(cmatrix, rmatrix, rsmat=soc_matrices, ispin=axis, &
                                   xkp=kpoint%xkp(:, ik), cell_to_index=cell_to_index, sab_nl=sab_nl)
               CALL dbcsr_scale(rmatrix, -1.0_dp)
               CALL dbcsr_desymmetrize(rmatrix, rmatrix_full)
               CALL dbcsr_desymmetrize(cmatrix, cmatrix_full)
               CALL copy_dbcsr_to_fm(rmatrix_full, mat_real)
               CALL copy_dbcsr_to_fm(cmatrix_full, mat_imag)
               CALL cp_fm_to_cfm(mat_real, mat_imag, omat_cfm)
               CALL parallel_gemm("N", "N", nao, nmo, nao, CMPLX(1.0_dp, 0.0_dp, dp), &
                                  omat_cfm, fmk1_cfm, CMPLX(0.0_dp, 0.0_dp, dp), tmp_cfm)
               CALL parallel_gemm("C", "N", nmo, nmo, nao, CMPLX(1.0_dp, 0.0_dp, dp), &
                                  fmk1_cfm, tmp_cfm, CMPLX(0.0_dp, 0.0_dp, dp), mmn_cfm)
               CALL cp_cfm_get_submatrix(mmn_cfm, soc_xyz(:, :, axis))
            END DO
            soc_h(1:nmo, 1:nmo) = soc_xyz(:, :, 3)
            soc_h(nmo + 1:, nmo + 1:) = -soc_xyz(:, :, 3)
            ! Match H_KS_spinor_kp and CP2K's real-spherical angular-momentum convention.
            soc_h(1:nmo, nmo + 1:) = soc_xyz(:, :, 1) + CMPLX(0.0_dp, 1.0_dp, dp)*soc_xyz(:, :, 2)
            soc_h(nmo + 1:, 1:nmo) = soc_xyz(:, :, 1) - CMPLX(0.0_dp, 1.0_dp, dp)*soc_xyz(:, :, 2)
            DO ib = 1, 2*nmo
               soc_h(ib, ib) = soc_h(ib, ib) + scalar_values(ib)
            END DO
            IF (MAXVAL(ABS(soc_h - CONJG(TRANSPOSE(soc_h)))) > 1.e-8_dp) THEN
               IF (iw > 0) WRITE (iw, '(T2,A,I0,A,ES18.10)') &
                  "TOPOLOGY| SOC Hermiticity error at point ", ik, ": ", &
                  MAXVAL(ABS(soc_h - CONJG(TRANSPOSE(soc_h))))
               CPABORT("WANNIER90 SOC: non-Hermitian second-variational Hamiltonian.")
            END IF
            CALL diag_complex(soc_h, soc_u, spinor_values(:, ik))
            spinor_coeff(:, :, ik) = soc_u(:, band_map)
         END DO
         CALL dbcsr_deallocate_matrix_set(soc_matrices)
         DEALLOCATE (soc_xyz, soc_h, soc_u, scalar_values)
         IF (iw > 0) WRITE (iw, '(T2,A,I0)') "TOPOLOGY| Scalar bands in SOC second variation: ", nmo
      END IF
      DO ispin = spin_channel, spin_channel
         ! loop over all k-points
         DO ik = 1, num_kpts
            ! get the MO coefficients for this k-point
            my_kpgrp = (ik >= kpoint%kp_range(1) .AND. ik <= kpoint%kp_range(2))
            IF (my_kpgrp) THEN
               ikk = ik - kpoint%kp_range(1) + 1
               kp => kpoint%kp_env(ikk)%kpoint_env
               cmo_coeff_source => kp%mos(ispin)%cmo_coeff
            ELSE
               NULLIFY (kp)
               NULLIFY (cmo_coeff_source)
            END IF
            CALL wannier90_copy_cfm_general(cmo_coeff_source, fmk1_cfm, para_env)
            ! loop over all connected neighbors
            DO i = 1, nntot
               ! get the MO coefficients for the connected k-point
               ik2 = nnlist(ik, i)
               mygrp = (ik2 >= kpoint%kp_range(1) .AND. ik2 <= kpoint%kp_range(2))
               IF (mygrp) THEN
                  ikk = ik2 - kpoint%kp_range(1) + 1
                  kp => kpoint%kp_env(ikk)%kpoint_env
                  cmo_coeff_source => kp%mos(ispin)%cmo_coeff
               ELSE
                  NULLIFY (kp)
                  NULLIFY (cmo_coeff_source)
               END IF
               CALL wannier90_copy_cfm_general(cmo_coeff_source, fmk2_cfm, para_env)
               !
               ! transfer realspace overlaps to connected k-point
               ibs = nblist(ik, i)
               ksign = SIGN(1.0_dp, REAL(ibs, KIND=dp))
               ibs = ABS(ibs)
               IF (ordered_berry) THEN
                  ! The cross-k AO operator is not Hermitian. Retain ordered atom pairs
                  ! and use exp(i*k'.R) [cos(b.r) - i*sin(b.r)] without symmetry shortcuts.
                  CALL dbcsr_set(rmatrix_full, 0.0_dp)
                  CALL dbcsr_set(cmatrix_full, 0.0_dp)
                  CALL dbcsr_set(loop_real, 0.0_dp)
                  CALL dbcsr_set(loop_imag, 0.0_dp)
                  CALL rskp_transform(rmatrix_full, cmatrix_full, berry_matrix(ibs)%cosmat, 1, &
                                      kpoint%xkp(:, ik2), berry_cell_index, overlap_nl)
                  CALL rskp_transform(loop_real, loop_imag, berry_matrix(ibs)%sinmat, 1, &
                                      kpoint%xkp(:, ik2), berry_cell_index, overlap_nl, rs_sign=ksign)
                  CALL dbcsr_add(rmatrix_full, loop_imag, 1.0_dp, -1.0_dp)
                  CALL dbcsr_add(cmatrix_full, loop_real, 1.0_dp, 1.0_dp)
               ELSE
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
               END IF
               CALL copy_dbcsr_to_fm(rmatrix_full, mat_real)
               CALL copy_dbcsr_to_fm(cmatrix_full, mat_imag)
               CALL cp_fm_to_cfm(mat_real, mat_imag, omat_cfm)
               CALL parallel_gemm("N", "N", nao, nmo, nao, CMPLX(1.0_dp, 0.0_dp, KIND=dp), &
                                  omat_cfm, fmk2_cfm, CMPLX(0.0_dp, 0.0_dp, KIND=dp), tmp_cfm)
               CALL parallel_gemm("C", "N", nmo, nmo, nao, CMPLX(1.0_dp, 0.0_dp, KIND=dp), &
                                  fmk1_cfm, tmp_cfm, CMPLX(0.0_dp, 0.0_dp, KIND=dp), mmn_cfm)
               IF (do_soc) THEN
                  CALL cp_cfm_get_submatrix(mmn_cfm, scalar_overlap)
                  link_matrix(:, :) = MATMUL(CONJG(TRANSPOSE(spinor_coeff(1:nmo, :, ik))), &
                                             MATMUL(scalar_overlap, spinor_coeff(1:nmo, :, ik2))) + &
                                      MATMUL(CONJG(TRANSPOSE(spinor_coeff(nmo + 1:, :, ik))), &
                                             MATMUL(scalar_overlap, spinor_coeff(nmo + 1:, :, ik2)))
               END IF
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
               DO ib2 = 1, num_bands
                  DO ib1 = 1, num_bands
                     IF (do_soc) THEN
                        rmmn = REAL(link_matrix(ib1, ib2), dp)
                        cmmn = AIMAG(link_matrix(ib1, ib2))
                     ELSE
                        CALL cp_cfm_get_element(mmn_cfm, band_map(ib1), band_map(ib2), zmmn)
                        rmmn = REAL(zmmn, KIND=dp)
                        cmmn = AIMAG(zmmn)
                     END IF
                     gauge_tmp = gauge_real*rmmn - gauge_imag*cmmn
                     cmmn = gauge_imag*rmmn + gauge_real*cmmn
                     rmmn = gauge_tmp
                     link_matrix(ib1, ib2) = CMPLX(rmmn, cmmn, dp)
                     IF (para_env%is_source()) THEN
                        WRITE (iunit, "(2E30.14)") rmmn, cmmn
                     END IF
                  END DO
               END DO
               IF (do_wilson) THEN
                  iloop = loop_index(ik)
                  CALL wilson_step(wilson_product(:, :, iloop), link_matrix, link_sv, status, 1.e-10_dp)
                  IF (status /= 0) THEN
                     CPABORT("Wilson link singular or SVD failed; refine sampling and check subspace.")
                  END IF
                  loop_sv(iloop) = MIN(loop_sv(iloop), link_sv)
               END IF
               !
            END DO
         END DO
      END DO
      ! Optional full-state snapshots, before the distributed work matrices are released.
      IF (export_state .OR. do_parity) THEN
         state_components = 1
         IF (do_soc) state_components = 2
         ALLOCATE (export_scalar(nao, nscalar), export_coeff(nao*state_components, num_bands))
         state_unit = -1
         IF (export_state) CALL topology_state_begin(qs_env_kp, TRIM(seed_name)//".topology", nao, num_bands, num_kpts, &
                                                     state_components, num_bands_tot, spin_channel, band_map, state_unit)
         IF (do_parity) THEN
            CPASSERT(ALL(kpoint%kp_range == [1, num_kpts]))
            ALLOCATE (parity_metric(nao, nao), parity_phase(nao), parity_map(nao), &
                      parity_odd(2**tqc_dimension), ebr_coefficients(2, 2**tqc_dimension), &
                      parity_energies(num_bands_tot))
            parity_odd = -1
            IF (do_tqc) THEN
               IF (num_bands /= nelectron .OR. MOD(num_bands, 2) /= 0 .OR. &
                   ANY(band_map /= [(i, i=1, num_bands)])) THEN
                  CPABORT("INVERSION_TQC requires the full occupied spinor band prefix.")
               END IF
            END IF
         END IF
         DO ik = 1, num_kpts
            kp => kpoint%kp_env(ik)%kpoint_env
            CALL cp_cfm_copy_general(kp%mos(spin_channel)%cmo_coeff, fmk1_cfm, para_env)
            CALL cp_cfm_get_submatrix(fmk1_cfm, export_scalar)
            IF (do_soc) THEN
               export_coeff(:nao, :) = MATMUL(export_scalar, spinor_coeff(:nscalar, :, ik))
               export_coeff(nao + 1:, :) = MATMUL(export_scalar, spinor_coeff(nscalar + 1:, :, ik))
               CALL topology_state_point(state_unit, ik, kpt_latt(:, ik), spinor_values(:, ik), export_coeff)
               IF (do_parity) parity_energies(:) = spinor_values(:, ik)
            ELSE
               export_coeff(:, :) = export_scalar(:, band_map)
               CALL get_mo_set(kp%mos(spin_channel), eigenvalues=eigenvalues)
               CALL topology_state_point(state_unit, ik, kpt_latt(:, ik), eigenvalues(:nscalar), export_coeff)
               IF (do_parity) parity_energies(:) = eigenvalues(:nscalar)
            END IF
            IF (.NOT. do_parity) CYCLE
            parity_gap = HUGE(1.0_dp)
            DO jpar = 1, num_bands_tot
               IF (ANY(band_map == jpar)) CYCLE
               parity_gap = MIN(parity_gap, MINVAL(ABS(parity_energies(band_map) - parity_energies(jpar))))
            END DO
            IF (num_bands == num_bands_tot .OR. parity_gap <= gap_tol) THEN
               CPABORT("PARITY requires an isolated selected subspace and computed excluded bands.")
            END IF
            IF (MAXVAL(ABS(2*kpt_latt(:, ik) - NINT(2*kpt_latt(:, ik)))) > 1.e-8_dp) CYCLE
            CALL gaussian_inversion_action(qs_env_kp, kpt_latt(:, ik), parity_origin, parity_map, &
                                           parity_phase, parity_tolerance, status)
            IF (status /= 0) THEN
               CPABORT("PARITY: inversion does not preserve geometry and atomic kinds.")
            END IF
            CALL dbcsr_set(rmatrix, 0.0_dp)
            CALL dbcsr_set(cmatrix, 0.0_dp)
            CALL rskp_transform(rmatrix, cmatrix, rsmat=matrix_s, ispin=1, xkp=kpoint%xkp(:, ik), &
                                cell_to_index=cell_to_index, sab_nl=sab_nl)
            CALL dbcsr_desymmetrize(rmatrix, rmatrix_full)
            CALL dbcsr_desymmetrize(cmatrix, cmatrix_full)
            CALL copy_dbcsr_to_fm(rmatrix_full, mat_real)
            CALL copy_dbcsr_to_fm(cmatrix_full, mat_imag)
            CALL cp_fm_to_cfm(mat_real, mat_imag, omat_cfm)
            CALL cp_cfm_get_submatrix(omat_cfm, parity_metric)
            CALL inversion_representation(parity_metric, export_coeff, parity_map, parity_phase, &
                                          parity_energies(band_map), do_soc .AND. time_reversal, parity_tolerance, &
                                          parity_energy_tolerance, counts, parity_error, parity_energy_error, status, parity_checks)
            IF (iw > 0) THEN
               WRITE (iw, '(T2,A,3F10.5,A,2I8,A,3ES14.5)') 'TQC| TRIM ', kpt_latt(:, ik), &
                  ' even/odd states ', counts, ' residuals/gap ', parity_error, parity_energy_error, parity_gap
               WRITE (iw, '(T2,A,4ES14.5)') 'TQC| Metric/normalization/inversion/TR residuals: ', parity_checks
            END IF
            IF (status /= 0) THEN
               CPABORT("PARITY: metric, subspace or symmetry representation check failed.")
            END IF
            IF (tqc_dimension == 2 .AND. MODULO(NINT(2*kpt_latt(3, ik)), 2) /= 0) CYCLE
            trim_id = 1
            DO axis = 1, tqc_dimension
               trim_id = trim_id + 2**(axis - 1)*MODULO(NINT(2*kpt_latt(axis, ik)), 2)
            END DO
            IF (parity_odd(trim_id) >= 0 .AND. parity_odd(trim_id) /= counts(2)) THEN
               CPABORT("PARITY: repeated equivalent TRIM have inconsistent multiplicities.")
            END IF
            parity_odd(trim_id) = counts(2)
         END DO
         IF (state_unit > 0) CALL close_file(state_unit)
         IF (do_parity) THEN
            IF (ALL(parity_odd < 0)) THEN
               CPABORT("PARITY requires a TRIM in the requested dimension.")
            END IF
         END IF
         IF (do_tqc) THEN
            IF (ANY(parity_odd < 0)) THEN
               CPABORT("INVERSION_TQC requires all four (2D) or eight (3D) TRIM.")
            END IF
            IF (ANY(MOD(parity_odd, 2) /= 0)) THEN
               CPABORT("INVERSION_TQC requires whole odd-parity Kramers pairs.")
            END IF
            parity_odd(:) = parity_odd/2
            CALL inversion_indicators(parity_odd, num_bands/2, tqc_dimension, strong, weak, z4, status)
            IF (status /= 0) CPABORT("Invalid inversion indicator data.")
            CALL inversion_ebr_signature(parity_odd, num_bands/2, signed_atomic, nonnegative_atomic, &
                                         ebr_coefficients, status)
            IF (status /= 0) CPABORT("Invalid inversion EBR data.")
            IF (iw > 0) THEN
               WRITE (iw, '(T2,A,I0)') 'TQC| Inversion-subgroup dimension: ', tqc_dimension
               WRITE (iw, '(T2,A,I0)') 'TQC| Fu-Kane parity index: ', strong
               IF (tqc_dimension == 3) THEN
                  WRITE (iw, '(T2,A,3I3)') 'TQC| Weak parity indices: ', weak
                  WRITE (iw, '(T2,A,I0)') 'TQC| Inversion Z4 (sum odd pairs mod 4): ', z4
               END IF
               WRITE (iw, '(T2,A,L1)') 'TQC| Signed atomic signature: ', signed_atomic
               WRITE (iw, '(T2,A,L1)') 'TQC| Nonnegative atomic signature: ', nonnegative_atomic
               IF (nonnegative_atomic) THEN
                  DO i = 1, SIZE(parity_odd)
                     WRITE (iw, '(T2,A,I0,A,2I8)') 'TQC| Center bit index ', i - 1, &
                        ' even/odd EBR pairs ', ebr_coefficients(:, i)
                  END DO
                  WRITE (iw, '(T2,A)') 'TQC| Atomic-compatible symmetry signature; not a proof of trivial topology.'
               ELSE IF (signed_atomic) THEN
                  WRITE (iw, '(T2,A)') 'TQC| Nonnegative EBR obstruction; '// &
                     'exclude hidden stable topology before calling it fragile.'
               ELSE
                  WRITE (iw, '(T2,A)') 'TQC| Stable inversion-symmetry indicator obstruction.'
               END IF
               WRITE (iw, '(T2,A)') 'TQC| TRIM checks do not establish a bulk gap; verify band isolation over the BZ.'
               WRITE (iw, '(T2,A)') 'TQC| Converge the scalar-state basis used for second-variational SOC.'
            END IF
         END IF
      END IF
      DO i = 1, nbs
         CALL dbcsr_deallocate_matrix_set(berry_matrix(i)%cosmat)
         CALL dbcsr_deallocate_matrix_set(berry_matrix(i)%sinmat)
      END DO
      DEALLOCATE (berry_matrix)
      CALL cp_fm_struct_release(matrix_struct_work)
      CALL cp_cfm_release(fmk1_cfm)
      CALL cp_cfm_release(fmk2_cfm)
      CALL cp_cfm_release(tmp_cfm)
      CALL cp_fm_struct_release(matrix_struct_ao)
      CALL cp_fm_release(mat_real)
      CALL cp_fm_release(mat_imag)
      CALL cp_cfm_release(omat_cfm)
      CALL cp_fm_struct_release(matrix_struct_mmn)
      CALL cp_cfm_release(mmn_cfm)
      CALL dbcsr_deallocate_matrix(rmatrix)
      CALL dbcsr_deallocate_matrix(cmatrix)
      CALL dbcsr_deallocate_matrix(rmatrix_full)
      CALL dbcsr_deallocate_matrix(cmatrix_full)
      IF (ordered_berry) THEN
         CALL dbcsr_deallocate_matrix(loop_real)
         CALL dbcsr_deallocate_matrix(loop_imag)
      END IF
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
      CALL get_mo_set(kp%mos(1), nmo=nmo)
      ALLOCATE (eigval(num_bands_tot))
      direct_gap = HUGE(1.0_dp)
      valence_max = -HUGE(1.0_dp)
      conduction_min = HUGE(1.0_dp)
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
         DO ispin = spin_channel, spin_channel
            IF (do_soc) THEN
               eigval(:) = spinor_values(:, ik)
            ELSE IF (my_kpgrp) THEN
               ikpgr = ik - kp_range(1) + 1
               kp => kpoint%kp_env(ikpgr)%kpoint_env
               CALL get_mo_set(kp%mos(ispin), eigenvalues=eigenvalues)
               eigval(1:nmo) = eigenvalues(1:nmo)
            ELSE
               eigval(1:nmo) = 0.0_dp
            END IF
            IF (.NOT. do_soc) CALL kpoint%para_env_inter_kp%sum(eigval)
            IF (do_wilson) THEN
               IF (lowest_bands) THEN
                  valence_max = MAX(valence_max, eigval(num_bands))
                  conduction_min = MIN(conduction_min, eigval(num_bands + 1))
               END IF
               DO ib = 1, num_bands_tot - 1
                  IF (ANY(band_map == ib) .NEQV. ANY(band_map == ib + 1)) THEN
                     pair_gap = eigval(ib + 1) - eigval(ib)
                     direct_gap = MIN(direct_gap, pair_gap)
                  END IF
               END DO
            END IF
            eigval(:) = eigval*evolt
            ! output
            IF (iunit > 0) THEN
               DO ib = 1, num_bands
                  WRITE (iunit, "(2I8,F24.14)") ib, ik, eigval(band_map(ib))
               END DO
            END IF
         END DO
      END DO
      IF (para_env%is_source()) THEN
         CALL close_file(iunit)
      END IF
      !
      IF (do_wilson) THEN
         IF (iw > 0 .AND. num_bands < num_bands_tot) WRITE (iw, '(T2,A,ES18.10)') &
            "TOPOLOGY| Minimum sampled subspace gap [eV]: ", direct_gap*evolt
         IF (direct_gap < gap_tol) CPABORT("Wilson subspace is not isolated on the sampled points.")
         IF (lowest_bands) THEN
            IF (iw > 0) THEN
               WRITE (iw, '(T2,A,ES18.10)') &
                  "TOPOLOGY| Maximum selected-band energy [eV]: ", valence_max*evolt
               WRITE (iw, '(T2,A,ES18.10)') &
                  "TOPOLOGY| Minimum excluded-band energy [eV]: ", conduction_min*evolt
               WRITE (iw, '(T2,A,ES18.10)') &
                  "TOPOLOGY| Sampled indirect gap [eV]: ", (conduction_min - valence_max)*evolt
               IF (conduction_min - valence_max <= gap_tol) WRITE (iw, '(T2,A)') &
                  "TOPOLOGY| Isolated band subspace, but no resolved common spectral gap on these samples."
            END IF
            IF (require_global_gap .AND. conduction_min - valence_max <= gap_tol) THEN
               CPABORT("No positive sampled indirect gap above the selected bands.")
            END IF
         END IF
         DO iloop = 1, nloop
            CALL wilson_spectrum(wilson_product(:, :, iloop), wcc_out(:, iloop), berry_phase, status)
            IF (status /= 0) CPABORT("Wilson eigenvalue calculation failed.")
            IF (iw > 0) WRITE (iw, '(T2,A,I0,A,F18.12,A,ES12.4)') &
               "TOPOLOGY| Loop ", iloop, " Berry phase [rad]: ", berry_phase, " min singular value: ", loop_sv(iloop)
         END DO
         IF (para_env%is_source()) THEN
            CALL open_file(TRIM(seed_name)//".wilson", unit_number=iunit, file_status="UNKNOWN", file_action="WRITE")
            WRITE (iunit, '(A)') "# loop, minimum link singular value, sorted WCC (arg(lambda)/(2*pi))"
            DO iloop = 1, nloop
               WRITE (iunit, '(I8,*(1X,ES24.16))') iloop, loop_sv(iloop), wcc_out(:, iloop)
            END DO
            CALL close_file(iunit)
         END IF
         IF (do_z2) THEN
            CALL z2_from_wcc(wcc_out, z2_value, status, 1.e-5_dp)
            IF (status /= 0 .AND. iw > 0) WRITE (iw, '(T2,A,I0)') &
               "TOPOLOGY| Kramers/crossing check requests refinement, status: ", status
            IF (iw > 0) WRITE (iw, '(T2,A,I0)') "TOPOLOGY| Z2 candidate (before sampling convergence): ", z2_value
         END IF
         IF (do_chern) THEN
            CALL chern_from_wcc(wcc_out, chern_value, chern_winding, status, 1.e-5_dp)
            IF (status /= 0) THEN
               chern_value = HUGE(0)
               IF (iw > 0) WRITE (iw, '(T2,A,I0)') &
                  "TOPOLOGY| Chern closure/winding check requests refinement, status: ", status
            ELSE IF (iw > 0) THEN
               WRITE (iw, '(T2,A,I0,A,F18.12)') &
                  "TOPOLOGY| First Chern candidate: ", chern_value, ", winding: ", chern_winding
            END IF
         END IF
      END IF
      ! clean up
      IF (ordered_berry) CALL kpoint_release(berry_kpoint)
      DEALLOCATE (kpt_latt, atoms_cart, atom_symbols, eigval)
      DEALLOCATE (nnlist, nncell)
      DEALLOCATE (nblist, b_latt)
      DEALLOCATE (band_map)
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
      TYPE(cp_cfm_type)                                  :: dst_cfm_full, src_cfm, src_cfm_full
      TYPE(cp_cfm_type), POINTER                         :: destination_cfm, source_cfm
      TYPE(cp_cfm_type), TARGET                          :: dst_cfm, empty_destination
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_source, matrix_struct_work
      TYPE(kpoint_env_type), POINTER                     :: kp, kp_source
      TYPE(kpoint_sym_type), POINTER                     :: kpsym

      success = .FALSE.
      reason = ""
      aligned_degenerate_blocks = 0
      aligned_degenerate_max_size = 0
      aligned_degenerate_min_svalue = HUGE(1.0_dp)
      NULLIFY (matrix_struct_source, matrix_struct_work)

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
      nspins = SIZE(kp%mos)
      IF (.NOT. ASSOCIATED(kp%mos(1)%cmo_coeff)) THEN
         reason = "Wannier90 export k-point MOs are not complex-valued"
         DEALLOCATE (source_kpoint, sym_index)
         RETURN
      END IF
      CALL get_mo_set(kp%mos(1), nao=nao, nmo=nmo)

      kp_source => qs_kpoint%kp_env(1)%kpoint_env
      IF (.NOT. ASSOCIATED(kp_source%mos(1)%cmo_coeff)) THEN
         reason = "SCF MOs are real-valued; complex symmetry phases cannot be reconstructed"
         DEALLOCATE (source_kpoint, sym_index)
         RETURN
      END IF
      CALL get_mo_set(kp_source%mos(1), nao=nao_src, nmo=nmo_src)
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
               CALL get_mo_set(kp_source%mos(ispin), eigenvalues=eigenvalues)
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
      CALL cp_cfm_create(src_cfm, matrix_struct_work)
      CALL cp_cfm_create(dst_cfm, matrix_struct_work)
      IF (source_window) THEN
         CALL cp_fm_struct_create(matrix_struct_source, nrow_global=nao, ncol_global=nmo_src, &
                                  para_env=para_env, context=blacs_env)
         CALL cp_cfm_create(src_cfm_full, matrix_struct_source)
         CALL cp_cfm_create(dst_cfm_full, matrix_struct_source)
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
            NULLIFY (source_cfm)
            IF (my_source_kpgrp) THEN
               ikpgr = ikred - source_kp_range(1) + 1
               kp_source => qs_kpoint%kp_env(ikpgr)%kpoint_env
               source_cfm => kp_source%mos(ispin)%cmo_coeff
               CALL get_mo_set(kp_source%mos(ispin), eigenvalues=eigenvalues, &
                               occupation_numbers=occupation)
               source_eigenvalues_buffer(1:nmo_src) = eigenvalues(1:nmo_src)
               source_occupation_buffer(1:nmo_src) = occupation(1:nmo_src)
            END IF
            IF (source_window) THEN
               CALL wannier90_copy_cfm_general(source_cfm, src_cfm_full, para_env)
            ELSE
               CALL wannier90_copy_cfm_general(source_cfm, src_cfm, para_env)
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
               CALL cp_cfm_to_cfm(src_cfm_full, src_cfm, nmo)
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
                     IF (.NOT. kpoint_same_periodic(kpoint%xkp(1:3, ik), &
                                                    kpsym%xkp(1:3, isym_try))) CYCLE
                     num_candidates = num_candidates + 1
                     CALL kpoint_transform_scf_mo(src_cfm, dst_cfm, &
                                                  qs_kpoint, ikred, isym_try, ok, &
                                                  candidate_reason)
                     IF (.NOT. ok) THEN
                        reason = candidate_reason
                        CYCLE
                     END IF
                     CALL ritz_stabilize_wannier90_subspace(dst_cfm, matrix_s, matrix_ks, &
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
                           CALL kpoint_transform_scf_mo(src_cfm_full, dst_cfm_full, qs_kpoint, &
                                                        ikred, isym_try, ok, candidate_reason)
                           IF (ok) THEN
                              CALL ritz_reconstruct_wannier90_window(dst_cfm_full, dst_cfm, matrix_s, &
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
                           CALL ritz_reconstruct_wannier90_window(dst_cfm, dst_cfm, matrix_s, matrix_ks, &
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
                  ELSE IF (.NOT. ok .AND. num_candidates == 0) THEN
                     reason = "no matching SCF symmetry operation candidate"
                  END IF
               ELSE
                  reason = "SCF k-point symmetry operation is not available"
               END IF
            ELSE
               CALL kpoint_transform_scf_mo(src_cfm, dst_cfm, qs_kpoint, &
                                            ikred, sym_index(ik), ok, reason)
            END IF
            IF (ok .AND. sym_index(ik) <= 0) THEN
               ! Even a direct k-point copy must be a closed H(k),S(k) subspace. This catches
               ! incomplete degenerate band windows before they can be exported to Wannier90.
               CALL ritz_stabilize_wannier90_subspace(dst_cfm, matrix_s, matrix_ks, &
                                                      kpoint%xkp(1:3, ik), cell_to_index, sab_nl, &
                                                      ispin, eigenvalues_buffer, degenerate_band_tol, &
                                                      ok, reason, aligned_blocks, aligned_max_size, &
                                                      aligned_min_svalue, candidate_residual)
               IF (.NOT. ok .AND. source_window) THEN
                  CALL kpoint_transform_scf_mo(src_cfm_full, dst_cfm_full, qs_kpoint, &
                                               ikred, sym_index(ik), ok, reason)
                  IF (ok) THEN
                     CALL ritz_reconstruct_wannier90_window(dst_cfm_full, dst_cfm, matrix_s, matrix_ks, &
                                                            kpoint%xkp(1:3, ik), cell_to_index, &
                                                            sab_nl, ispin, eigenvalues_buffer, nmo, ok, &
                                                            reason, source_window_min_svalue, &
                                                            candidate_residual)
                  END IF
               ELSE IF (.NOT. ok) THEN
                  CALL ritz_reconstruct_wannier90_window(dst_cfm, dst_cfm, matrix_s, matrix_ks, &
                                                         kpoint%xkp(1:3, ik), &
                                                         cell_to_index, sab_nl, ispin, &
                                                         eigenvalues_buffer, nmo, ok, reason, &
                                                         source_window_min_svalue, candidate_residual)
               END IF
            END IF
            IF (.NOT. ok) THEN
               CALL cp_cfm_release(src_cfm)
               CALL cp_cfm_release(dst_cfm)
               CALL cp_fm_struct_release(matrix_struct_work)
               IF (source_window) THEN
                  CALL cp_cfm_release(src_cfm_full)
                  CALL cp_cfm_release(dst_cfm_full)
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
               CALL get_mo_set(kp%mos(ispin), eigenvalues=eigenvalues, &
                               occupation_numbers=occupation)
               eigenvalues(1:nmo) = eigenvalues_buffer(1:nmo)
               occupation(1:nmo) = occupation_buffer(1:nmo)
               destination_cfm => kp%mos(ispin)%cmo_coeff
            ELSE
               destination_cfm => empty_destination
            END IF
            CALL cp_cfm_copy_general(dst_cfm, destination_cfm, para_env)
         END DO
      END DO

      CALL cp_cfm_release(src_cfm)
      CALL cp_cfm_release(dst_cfm)
      CALL cp_fm_struct_release(matrix_struct_work)
      IF (source_window) THEN
         CALL cp_cfm_release(src_cfm_full)
         CALL cp_cfm_release(dst_cfm_full)
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
!> \param mo MO coefficients, indexed as AO, MO, k-point, spin
!> \param eigenvalue_snapshot MO eigenvalues, indexed as MO, k-point, spin
! **************************************************************************************************
   SUBROUTINE save_wannier90_mo_snapshot(kpoint, nspins, para_env, mo, &
                                         eigenvalue_snapshot)
      TYPE(kpoint_type), POINTER                         :: kpoint
      INTEGER, INTENT(IN)                                :: nspins
      TYPE(mp_para_env_type), POINTER                    :: para_env
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :), INTENT(OUT)              :: mo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: eigenvalue_snapshot

      INTEGER                                            :: ik, ikpgr, ispin, nao, nkp, nmo
      INTEGER, DIMENSION(2)                              :: kp_range
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: owner_weight
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(kpoint_env_type), POINTER                     :: kp

      CALL get_kpoint_info(kpoint, nkp=nkp, kp_range=kp_range)
      kp => kpoint%kp_env(1)%kpoint_env
      CALL get_mo_set(kp%mos(1), nao=nao, nmo=nmo)
      ALLOCATE (mo(nao, nmo, nkp, nspins), &
                eigenvalue_snapshot(nmo, nkp, nspins), owner_weight(nkp, nspins))
      mo(:, :, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      eigenvalue_snapshot(:, :, :) = 0.0_dp
      owner_weight(:, :) = 0.0_dp
      DO ik = kp_range(1), kp_range(2)
         ikpgr = ik - kp_range(1) + 1
         kp => kpoint%kp_env(ikpgr)%kpoint_env
         DO ispin = 1, nspins
            CALL cp_cfm_get_submatrix(kp%mos(ispin)%cmo_coeff, mo(:, :, ik, ispin))
            CALL get_mo_set(kp%mos(ispin), eigenvalues=eigenvalues)
            eigenvalue_snapshot(1:nmo, ik, ispin) = eigenvalues(1:nmo)
            owner_weight(ik, ispin) = 1.0_dp
         END DO
      END DO
      CALL para_env%sum(mo)
      CALL para_env%sum(eigenvalue_snapshot)
      CALL para_env%sum(owner_weight)
      DO ik = 1, nkp
         DO ispin = 1, nspins
            IF (owner_weight(ik, ispin) > 0.0_dp) THEN
               mo(:, :, ik, ispin) = mo(:, :, ik, ispin)/owner_weight(ik, ispin)
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
!> \param mo complex MO coefficient snapshot
!> \param eigenvalue_snapshot MO eigenvalue snapshot
! **************************************************************************************************
   SUBROUTINE restore_wannier90_mo_snapshot(kpoint, mo, eigenvalue_snapshot)
      TYPE(kpoint_type), POINTER                         :: kpoint
      COMPLEX(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(IN)                                      :: mo
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: eigenvalue_snapshot

      INTEGER                                            :: ik, ikpgr, ispin, nmo, nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(kpoint_env_type), POINTER                     :: kp

      CALL get_kpoint_info(kpoint, kp_range=kp_range)
      nmo = SIZE(eigenvalue_snapshot, 1)
      nspins = SIZE(eigenvalue_snapshot, 3)
      DO ik = kp_range(1), kp_range(2)
         ikpgr = ik - kp_range(1) + 1
         kp => kpoint%kp_env(ikpgr)%kpoint_env
         DO ispin = 1, nspins
            CALL cp_cfm_set_submatrix(kp%mos(ispin)%cmo_coeff, mo(:, :, ik, ispin))
            CALL get_mo_set(kp%mos(ispin), eigenvalues=eigenvalues)
            eigenvalues(1:nmo) = eigenvalue_snapshot(1:nmo, ik, ispin)
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
!> \param reference_mo complex MO coefficient reference
!> \param reference_eigenvalues MO eigenvalue reference
!> \param success true if the reconstructed MOs match the reference subspaces
!> \param max_subspace_deviation largest deviation of S(k)-metric singular values from one
!> \param min_svalue smallest S(k)-metric singular value
!> \param max_eigenvalue_deviation largest eigenvalue deviation
! **************************************************************************************************
   SUBROUTINE validate_wannier90_reused_mos(kpoint, matrix_s, cell_to_index, sab_nl, para_env, &
                                            reference_mo, reference_eigenvalues, &
                                            success, max_subspace_deviation, min_svalue, &
                                            max_eigenvalue_deviation)
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(mp_para_env_type), POINTER                    :: para_env
      COMPLEX(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(IN)                                      :: reference_mo
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
      TYPE(cp_cfm_type)                                  :: cand_cfm, ref_cfm
      TYPE(cp_cfm_type), POINTER                         :: candidate_source
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_work
      TYPE(kpoint_env_type), POINTER                     :: kp

      success = .FALSE.
      max_subspace_deviation = 0.0_dp
      min_svalue = HUGE(1.0_dp)
      max_eigenvalue_deviation = 0.0_dp
      NULLIFY (matrix_struct_work)

      CALL get_kpoint_info(kpoint, nkp=nkp, kp_range=kp_range)
      kp => kpoint%kp_env(1)%kpoint_env
      nspins = SIZE(kp%mos)
      CALL get_mo_set(kp%mos(1), nao=nao, nmo=nmo)
      CPASSERT(SIZE(reference_mo, 1) == nao)
      CPASSERT(SIZE(reference_mo, 2) == nmo)
      CPASSERT(SIZE(reference_mo, 3) == nkp)
      CPASSERT(SIZE(reference_mo, 4) == nspins)

      matrix_struct_work => kp%mos(1)%cmo_coeff%matrix_struct
      CALL cp_cfm_create(ref_cfm, matrix_struct_work)
      CALL cp_cfm_create(cand_cfm, matrix_struct_work)
      ALLOCATE (eigenvalue_buffer(nmo))

      DO ispin = 1, nspins
         DO ik = 1, nkp
            CALL cp_cfm_set_submatrix(ref_cfm, reference_mo(:, :, ik, ispin))
            my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
            NULLIFY (candidate_source)
            IF (my_kpgrp) THEN
               ikpgr = ik - kp_range(1) + 1
               kp => kpoint%kp_env(ikpgr)%kpoint_env
               candidate_source => kp%mos(ispin)%cmo_coeff
               CALL get_mo_set(kp%mos(ispin), eigenvalues=eigenvalues)
               eigenvalue_buffer(1:nmo) = eigenvalues(1:nmo)
            ELSE
               eigenvalue_buffer(1:nmo) = 0.0_dp
            END IF
            CALL wannier90_copy_cfm_general(candidate_source, cand_cfm, para_env)
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
            CALL measure_wannier90_subspace_error(ref_cfm, cand_cfm, matrix_s, &
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
      CALL cp_cfm_release(ref_cfm)
      CALL cp_cfm_release(cand_cfm)

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
      TYPE(cp_cfm_type)                                  :: dst_cfm, ref_cfm, src_cfm, src_cfm_full
      TYPE(cp_cfm_type), POINTER                         :: source_cfm
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_source, matrix_struct_work
      TYPE(kpoint_env_type), POINTER                     :: kp, kp_source
      TYPE(kpoint_sym_type), POINTER                     :: kpsym

      max_subspace_deviation = HUGE(1.0_dp)
      min_svalue = 0.0_dp
      max_metric_deviation = HUGE(1.0_dp)
      max_residual = HUGE(1.0_dp)
      NULLIFY (matrix_struct_source, matrix_struct_work)

      CALL build_wannier90_scf_mapping(kpoint, qs_kpoint, source_kpoint, sym_index, ok, reason)
      IF (.NOT. ok) RETURN
      kp => kpoint%kp_env(1)%kpoint_env
      kp_source => qs_kpoint%kp_env(1)%kpoint_env
      IF (.NOT. ASSOCIATED(kp%mos(1)%cmo_coeff) .OR. &
          .NOT. ASSOCIATED(kp_source%mos(1)%cmo_coeff)) THEN
         DEALLOCATE (source_kpoint, sym_index)
         RETURN
      END IF
      nspins = SIZE(kp%mos)
      CALL get_mo_set(kp%mos(1), nao=nao, nmo=nmo)
      CALL get_mo_set(kp_source%mos(1), nao=nao_src, nmo=nmo_src)
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
      CALL cp_cfm_create(ref_cfm, matrix_struct_work)
      CALL cp_cfm_create(src_cfm, matrix_struct_work)
      CALL cp_cfm_create(dst_cfm, matrix_struct_work)
      ALLOCATE (eigenvalues_buffer(nmo))
      IF (source_window) THEN
         CALL cp_fm_struct_create(matrix_struct_source, nrow_global=nao, ncol_global=nmo_src, &
                                  para_env=para_env, context=kpoint%blacs_env)
         CALL cp_cfm_create(src_cfm_full, matrix_struct_source)
      END IF

      max_subspace_deviation = 0.0_dp
      min_svalue = HUGE(1.0_dp)
      max_metric_deviation = 0.0_dp
      max_residual = 0.0_dp
      DO ispin = 1, nspins
         DO ik = 1, nkp
            my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
            NULLIFY (source_cfm)
            IF (my_kpgrp) THEN
               ikpgr = ik - kp_range(1) + 1
               kp => kpoint%kp_env(ikpgr)%kpoint_env
               source_cfm => kp%mos(ispin)%cmo_coeff
            END IF
            CALL wannier90_copy_cfm_general(source_cfm, ref_cfm, para_env)

            ikred = source_kpoint(ik)
            my_kpgrp = (ikred >= source_kp_range(1) .AND. ikred <= source_kp_range(2))
            NULLIFY (source_cfm)
            IF (my_kpgrp) THEN
               ikpgr = ikred - source_kp_range(1) + 1
               kp_source => qs_kpoint%kp_env(ikpgr)%kpoint_env
               source_cfm => kp_source%mos(ispin)%cmo_coeff
               CALL get_mo_set(kp_source%mos(ispin), eigenvalues=eigenvalues)
               eigenvalues_buffer(1:nmo) = eigenvalues(1:nmo)
            ELSE
               eigenvalues_buffer(1:nmo) = 0.0_dp
            END IF
            IF (source_window) THEN
               CALL wannier90_copy_cfm_general(source_cfm, src_cfm_full, para_env)
            ELSE
               CALL wannier90_copy_cfm_general(source_cfm, src_cfm, para_env)
            END IF
            IF (my_kpgrp) THEN
               owner_count = 1.0_dp
            ELSE
               owner_count = 0.0_dp
            END IF
            CALL para_env%sum(owner_count)
            CALL para_env%sum(eigenvalues_buffer)
            IF (owner_count > 0.0_dp) eigenvalues_buffer(1:nmo) = eigenvalues_buffer(1:nmo)/owner_count
            CALL measure_wannier90_eigenspace_quality(ref_cfm, matrix_s, matrix_ks, &
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
               CALL cp_cfm_to_cfm(src_cfm_full, src_cfm, nmo)
            END IF

            best_deviation = HUGE(1.0_dp)
            best_metric_deviation = 0.0_dp
            best_residual = 0.0_dp
            best_svalue = 0.0_dp
            IF (sym_index(ik) <= 0) THEN
               CALL kpoint_transform_scf_mo(src_cfm, dst_cfm, qs_kpoint, &
                                            ikred, sym_index(ik), ok, reason)
               IF (ok) THEN
                  CALL measure_wannier90_subspace_error(ref_cfm, dst_cfm, matrix_s, &
                                                        kpoint%xkp(1:3, ik), cell_to_index, sab_nl, &
                                                        para_env, ok, candidate_deviation, candidate_svalue)
                  IF (ok) THEN
                     CALL measure_wannier90_eigenspace_quality(dst_cfm, matrix_s, matrix_ks, &
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
            ELSE IF (ASSOCIATED(qs_kpoint%kp_sym)) THEN
               kpsym => qs_kpoint%kp_sym(ikred)%kpoint_sym
               IF (ASSOCIATED(kpsym)) THEN
                  DO isym_try = 1, kpsym%nwred
                     IF (.NOT. kpoint_same_periodic(kpoint%xkp(1:3, ik), &
                                                    kpsym%xkp(1:3, isym_try))) CYCLE
                     CALL kpoint_transform_scf_mo(src_cfm, dst_cfm, &
                                                  qs_kpoint, ikred, isym_try, ok, reason)
                     IF (.NOT. ok) CYCLE
                     CALL measure_wannier90_subspace_error(ref_cfm, dst_cfm, &
                                                           matrix_s, kpoint%xkp(1:3, ik), cell_to_index, &
                                                           sab_nl, para_env, ok, candidate_deviation, &
                                                           candidate_svalue)
                     IF (ok) THEN
                        CALL measure_wannier90_eigenspace_quality(dst_cfm, matrix_s, matrix_ks, &
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
         CALL cp_cfm_release(src_cfm_full)
         CALL cp_fm_struct_release(matrix_struct_source)
      END IF
      CALL cp_cfm_release(ref_cfm)
      CALL cp_cfm_release(src_cfm)
      CALL cp_cfm_release(dst_cfm)
      CALL cp_fm_struct_release(matrix_struct_work)
      DEALLOCATE (eigenvalues_buffer)
      DEALLOCATE (source_kpoint, sym_index)

   END SUBROUTINE diagnose_wannier90_scf_reuse_candidates

! **************************************************************************************************
!> \brief Measure the S(k)-metric distance between two Wannier90 MO subspaces.
!> \param ref_cfm reference MO coefficients
!> \param cand_cfm candidate MO coefficients
!> \param matrix_s real-space overlap matrix
!> \param xkp target k-point coordinate
!> \param cell_to_index real-space cell index table
!> \param sab_nl overlap neighbor list
!> \param para_env global parallel environment
!> \param success true if the metric comparison was performed
!> \param max_subspace_deviation largest deviation of singular values from one
!> \param min_svalue smallest singular value of C_ref^+ S(k) C_candidate
! **************************************************************************************************
   SUBROUTINE measure_wannier90_subspace_error(ref_cfm, cand_cfm, matrix_s, &
                                               xkp, cell_to_index, sab_nl, para_env, success, &
                                               max_subspace_deviation, min_svalue)
      TYPE(cp_cfm_type), INTENT(IN)                      :: ref_cfm, cand_cfm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(OUT)                               :: success
      REAL(KIND=dp), INTENT(OUT)                         :: max_subspace_deviation, min_svalue

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: metric_projected, metric_vectors, &
                                                            overlap, ref_coeff, s_cand_coeff
      INTEGER                                            :: ib, nao, nmo, nmo_candidate
      REAL(KIND=dp)                                      :: singular_value
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: metric_values
      TYPE(cp_cfm_type)                                  :: s_cand
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_metric

      success = .FALSE.
      max_subspace_deviation = HUGE(1.0_dp)
      min_svalue = 0.0_dp
      NULLIFY (matrix_struct_metric)

      CALL cp_cfm_get_info(ref_cfm, nrow_global=nao, ncol_global=nmo, &
                           matrix_struct=matrix_struct_metric)
      CALL cp_cfm_get_info(cand_cfm, ncol_global=nmo_candidate)
      IF (nmo_candidate /= nmo) RETURN

      CALL cp_cfm_create(s_cand, matrix_struct_metric)
      CALL apply_wannier90_kp_matrix(matrix_s, 1, xkp, cell_to_index, sab_nl, &
                                     cand_cfm, s_cand)

      ALLOCATE (ref_coeff(nao, nmo), s_cand_coeff(nao, nmo))
      CALL cp_cfm_get_submatrix(ref_cfm, ref_coeff)
      CALL cp_cfm_get_submatrix(s_cand, s_cand_coeff)

      ALLOCATE (overlap(nmo, nmo), &
                metric_projected(nmo, nmo), metric_vectors(nmo, nmo), metric_values(nmo))
      overlap(:, :) = MATMUL(CONJG(TRANSPOSE(ref_coeff)), s_cand_coeff)
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

      DEALLOCATE (ref_coeff, s_cand_coeff, overlap, metric_projected, metric_vectors, metric_values)
      CALL cp_cfm_release(s_cand)

   END SUBROUTINE measure_wannier90_subspace_error

! **************************************************************************************************
!> \brief Measure whether transformed MOs are an H(k),S(k) invariant eigenspace.
!> \param cand_cfm candidate MO coefficients
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
   SUBROUTINE measure_wannier90_eigenspace_quality(cand_cfm, matrix_s, matrix_ks, xkp, &
                                                   cell_to_index, sab_nl, para_env, ispin, eigenvalues, &
                                                   success, metric_deviation, min_metric_eigenvalue, &
                                                   residual_norm)
      TYPE(cp_cfm_type), INTENT(IN)                      :: cand_cfm
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

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: h_coeff, metric_vectors, residual_block, &
                                                            s_coeff, s_projected
      INTEGER                                            :: ib, nao, nmo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: metric_values
      TYPE(cp_cfm_type)                                  :: h_cand, metric_cfm, s_cand
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_metric, &
                                                            matrix_struct_projected

      success = .FALSE.
      metric_deviation = HUGE(1.0_dp)
      min_metric_eigenvalue = 0.0_dp
      residual_norm = HUGE(1.0_dp)
      NULLIFY (matrix_struct_metric, matrix_struct_projected)

      CALL cp_cfm_get_info(cand_cfm, nrow_global=nao, ncol_global=nmo, &
                           matrix_struct=matrix_struct_metric)
      IF (SIZE(eigenvalues) < nmo) RETURN

      CALL cp_cfm_create(s_cand, matrix_struct_metric)
      CALL cp_cfm_create(h_cand, matrix_struct_metric)

      CALL apply_wannier90_kp_matrix(matrix_s, 1, xkp, cell_to_index, sab_nl, &
                                     cand_cfm, s_cand)
      CALL apply_wannier90_kp_matrix(matrix_ks, ispin, xkp, cell_to_index, sab_nl, &
                                     cand_cfm, h_cand)

      ALLOCATE (h_coeff(nao, nmo), metric_vectors(nmo, nmo), &
                residual_block(nao, nmo), s_coeff(nao, nmo), s_projected(nmo, nmo), &
                metric_values(nmo))
      CALL cp_cfm_get_submatrix(s_cand, s_coeff)
      CALL cp_cfm_get_submatrix(h_cand, h_coeff)
      CALL cp_fm_struct_create(matrix_struct_projected, nrow_global=nmo, ncol_global=nmo, &
                               para_env=matrix_struct_metric%para_env, &
                               context=matrix_struct_metric%context)
      CALL cp_cfm_create(metric_cfm, matrix_struct_projected)
      CALL parallel_gemm("C", "N", nmo, nmo, nao, CMPLX(1.0_dp, 0.0_dp, KIND=dp), cand_cfm, &
                         s_cand, CMPLX(0.0_dp, 0.0_dp, KIND=dp), metric_cfm)
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

      DEALLOCATE (h_coeff, metric_vectors, residual_block, s_coeff, s_projected, metric_values)
      CALL cp_cfm_release(s_cand)
      CALL cp_cfm_release(h_cand)
      CALL cp_cfm_release(metric_cfm)
      CALL cp_fm_struct_release(matrix_struct_projected)

   END SUBROUTINE measure_wannier90_eigenspace_quality

! **************************************************************************************************
!> \brief Apply a complex k-point matrix to a complex MO coefficient matrix.
!> \param rsmat real-space matrix images
!> \param ispin spin index for rsmat
!> \param xkp target k-point coordinate
!> \param cell_to_index real-space cell index table
!> \param sab_nl overlap neighbor list
!> \param coeff input MO coefficients
!> \param RESULT output of applying the k-point matrix to the MO coefficients
! **************************************************************************************************
   SUBROUTINE apply_wannier90_kp_matrix(rsmat, ispin, xkp, cell_to_index, sab_nl, &
                                        coeff, RESULT)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(cp_cfm_type), INTENT(IN)                      :: coeff
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: result

      INTEGER                                            :: nao, ncol
      TYPE(cp_cfm_type)                                  :: kmat_cfm
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_ao, matrix_struct_coeff
      TYPE(cp_fm_type)                                   :: mat_imag, mat_real
      TYPE(dbcsr_type), POINTER                          :: kmat_imag, kmat_imag_full, kmat_real, &
                                                            kmat_real_full

      NULLIFY (matrix_struct_ao, matrix_struct_coeff, kmat_imag, kmat_imag_full, kmat_real, &
               kmat_real_full)

      CALL cp_cfm_get_info(coeff, nrow_global=nao, ncol_global=ncol, &
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
      CALL cp_fm_to_cfm(mat_real, mat_imag, kmat_cfm)
      CALL parallel_gemm("N", "N", nao, ncol, nao, CMPLX(1.0_dp, 0.0_dp, KIND=dp), kmat_cfm, &
                         coeff, CMPLX(0.0_dp, 0.0_dp, KIND=dp), RESULT)

      CALL cp_fm_release(mat_real)
      CALL cp_fm_release(mat_imag)
      CALL cp_cfm_release(kmat_cfm)
      CALL cp_fm_struct_release(matrix_struct_ao)
      CALL dbcsr_deallocate_matrix(kmat_real)
      CALL dbcsr_deallocate_matrix(kmat_imag)
      CALL dbcsr_deallocate_matrix(kmat_real_full)
      CALL dbcsr_deallocate_matrix(kmat_imag_full)

   END SUBROUTINE apply_wannier90_kp_matrix

! **************************************************************************************************
!> \brief Rayleigh-Ritz stabilize a symmetry-reconstructed Wannier90 MO subspace.
!> \param dst_cfm transformed MO coefficients
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
   SUBROUTINE ritz_stabilize_wannier90_subspace(dst_cfm, matrix_s, matrix_ks, &
                                                xkp, cell_to_index, sab_nl, ispin, eigenvalues, &
                                                degenerate_band_tol, success, reason, aligned_blocks, &
                                                aligned_max_size, aligned_min_svalue, max_residual)
      TYPE(cp_cfm_type), INTENT(IN)                      :: dst_cfm
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

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: block_coeff, dst_coeff, h_block, h_coeff, &
         h_projected, h_projected_work, metric_vectors, residual_block, ritz_vectors, s_block, &
         s_coeff, s_projected, stabilized
      INTEGER                                            :: block_first, block_last, block_size, ib, &
                                                            nao, nmo
      REAL(KIND=dp)                                      :: metric_deviation, norm_value, &
                                                            residual_norm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: metric_values, ritz_values
      TYPE(cp_cfm_type)                                  :: h_dst, s_dst
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_metric

      success = .FALSE.
      reason = ""
      aligned_blocks = 0
      aligned_max_size = 0
      aligned_min_svalue = HUGE(1.0_dp)
      max_residual = 0.0_dp

      NULLIFY (matrix_struct_metric)
      CALL cp_cfm_get_info(dst_cfm, nrow_global=nao, ncol_global=nmo, &
                           matrix_struct=matrix_struct_metric)
      IF (SIZE(eigenvalues) < nmo) THEN
         reason = "not enough eigenvalues for Wannier90 Ritz subspace stabilization"
         RETURN
      END IF

      CALL cp_cfm_create(s_dst, matrix_struct_metric)
      CALL cp_cfm_create(h_dst, matrix_struct_metric)

      stabilize_subspace: BLOCK
         CALL apply_wannier90_kp_matrix(matrix_s, 1, xkp, cell_to_index, sab_nl, &
                                        dst_cfm, s_dst)
         CALL apply_wannier90_kp_matrix(matrix_ks, ispin, xkp, cell_to_index, sab_nl, &
                                        dst_cfm, h_dst)

         ALLOCATE (dst_coeff(nao, nmo), s_coeff(nao, nmo), h_coeff(nao, nmo))
         CALL cp_cfm_get_submatrix(dst_cfm, dst_coeff)
         CALL cp_cfm_get_submatrix(s_dst, s_coeff)
         CALL cp_cfm_get_submatrix(h_dst, h_coeff)

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
               block_coeff(:, :) = dst_coeff(:, block_first:block_last)
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
                  EXIT stabilize_subspace
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
                  EXIT stabilize_subspace
               END IF

               dst_coeff(:, block_first:block_last) = stabilized
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
            EXIT stabilize_subspace
         END IF

         CALL cp_cfm_set_submatrix(dst_cfm, dst_coeff)

         IF (aligned_blocks == 0) aligned_min_svalue = 0.0_dp
         success = .TRUE.
      END BLOCK stabilize_subspace

      CALL cp_cfm_release(s_dst)
      CALL cp_cfm_release(h_dst)

   END SUBROUTINE ritz_stabilize_wannier90_subspace

! **************************************************************************************************
!> \brief Reconstruct a Wannier90 export window from a larger symmetry-transformed SCF MO space.
!> \param src_cfm transformed source MO window
!> \param dst_cfm exported reconstructed MO coefficients
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
   SUBROUTINE ritz_reconstruct_wannier90_window(src_cfm, dst_cfm, matrix_s, &
                                                matrix_ks, xkp, cell_to_index, sab_nl, ispin, &
                                                eigenvalues, nmo_export, success, reason, min_svalue, &
                                                max_residual)
      TYPE(cp_cfm_type), INTENT(IN)                      :: src_cfm, dst_cfm
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
      TYPE(cp_cfm_type)                                  :: h_src, s_src
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_metric

      success = .FALSE.
      reason = ""
      min_svalue = HUGE(1.0_dp)
      max_residual = HUGE(1.0_dp)

      NULLIFY (matrix_struct_metric)
      CALL cp_cfm_get_info(src_cfm, nrow_global=nao, ncol_global=nmo_source, &
                           matrix_struct=matrix_struct_metric)
      IF (nmo_export > nmo_source) THEN
         reason = "Wannier90 export window is larger than the transformed SCF MO space"
         RETURN
      END IF
      IF (SIZE(eigenvalues) < nmo_export) THEN
         reason = "not enough eigenvalue storage for Wannier90 source-window reconstruction"
         RETURN
      END IF

      CALL cp_cfm_create(s_src, matrix_struct_metric)
      CALL cp_cfm_create(h_src, matrix_struct_metric)

      CALL apply_wannier90_kp_matrix(matrix_s, 1, xkp, cell_to_index, sab_nl, &
                                     src_cfm, s_src)
      CALL apply_wannier90_kp_matrix(matrix_ks, ispin, xkp, cell_to_index, sab_nl, &
                                     src_cfm, h_src)

      ALLOCATE (source_coeff(nao, nmo_source), s_coeff(nao, nmo_source), &
                h_coeff(nao, nmo_source), h_projected(nmo_source, nmo_source), &
                h_projected_work(nmo_source, nmo_source), metric_vectors(nmo_source, nmo_source), &
                residual_block(nao, nmo_export), ritz_vectors(nmo_source, nmo_source), &
                s_projected(nmo_source, nmo_source), stabilized(nao, nmo_export), &
                coeff_work(nao, nmo_source), metric_values(nmo_source), ritz_values(nmo_source), &
                source_eigenvalues(nmo_export))
      source_eigenvalues(1:nmo_export) = eigenvalues(1:nmo_export)
      CALL cp_cfm_get_submatrix(src_cfm, source_coeff)
      CALL cp_cfm_get_submatrix(s_src, s_coeff)
      CALL cp_cfm_get_submatrix(h_src, h_coeff)
      s_projected(:, :) = MATMUL(CONJG(TRANSPOSE(source_coeff)), s_coeff)
      h_projected(:, :) = MATMUL(CONJG(TRANSPOSE(source_coeff)), h_coeff)
      s_projected(:, :) = 0.5_dp*(s_projected + CONJG(TRANSPOSE(s_projected)))
      h_projected(:, :) = 0.5_dp*(h_projected + CONJG(TRANSPOSE(h_projected)))

      reconstruct_window: BLOCK
         CALL diag_complex(s_projected, metric_vectors, metric_values)
         min_svalue = MINVAL(metric_values)
         metric_deviation = MAXVAL(ABS(metric_values - 1.0_dp))
         IF (min_svalue < 1.0e-10_dp) THEN
            WRITE (reason, "(A,ES9.2,A,ES9.2)") &
               "singular expanded metric smin=", min_svalue, " dS=", metric_deviation
            EXIT reconstruct_window
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
            EXIT reconstruct_window
         END IF
         max_eigenvalue_shift = MAXVAL(ABS(ritz_values(1:nmo_export) - source_eigenvalues(1:nmo_export)))
         IF (max_eigenvalue_shift > eigenvalue_tol) THEN
            WRITE (reason, "(A,ES9.2)") &
               "expanded dS=", metric_deviation
            EXIT reconstruct_window
         END IF

         CALL cp_cfm_set_submatrix(dst_cfm, stabilized)
         success = .TRUE.

      END BLOCK reconstruct_window

      DEALLOCATE (source_coeff, s_coeff, h_coeff, h_projected, h_projected_work, metric_vectors, &
                  residual_block, ritz_vectors, s_projected, stabilized, coeff_work, metric_values, &
                  ritz_values, source_eigenvalues)
      CALL cp_cfm_release(s_src)
      CALL cp_cfm_release(h_src)

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
         IF (kpoint_same_periodic(xkp_mesh(1:3, ik), xkp_search)) THEN
            ik_match = ik
            RETURN
         END IF
      END DO

   END FUNCTION find_matching_kpoint

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
