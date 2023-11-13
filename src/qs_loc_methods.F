!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Driver for the localization that should be general
!>      for all the methods available and all the definition of the
!>      spread functional
!>      Write centers, spread and cubes only if required and for the
!>      selected states
!>      The localized functions are copied in the standard mos array
!>      for the next use
!> \par History
!>      01.2008 Teodoro Laino [tlaino] - University of Zurich
!>        - Merging the two localization codes and updating to new structures
!> \author MI (04.2005)
! **************************************************************************************************
MODULE qs_loc_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              deallocate_atomic_kind_set,&
                                              set_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_rot_cols,&
                                              cp_fm_rot_rows,&
                                              cp_fm_schur_product
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pool_create_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_element,&
                                              cp_fm_get_info,&
                                              cp_fm_maxabsval,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_iter_string,&
                                              cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE input_constants,                 ONLY: &
        do_loc_crazy, do_loc_direct, do_loc_gapo, do_loc_jacobi, do_loc_l1_norm_sd, do_loc_none, &
        do_loc_scdm, dump_dcd, dump_dcd_aligned_cell, dump_xmol
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_ivals,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              sp
   USE machine,                         ONLY: m_flush
   USE mathconstants,                   ONLY: pi,&
                                              twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE motion_utils,                    ONLY: get_output_format
   USE orbital_pointers,                ONLY: ncoset
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_methods,                ONLY: get_particle_set,&
                                              write_particle_coordinates
   USE particle_types,                  ONLY: allocate_particle_set,&
                                              deallocate_particle_set,&
                                              particle_type
   USE physcon,                         ONLY: angstrom
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_loc_types,                    ONLY: get_qs_loc_env,&
                                              qs_loc_env_type
   USE qs_localization_methods,         ONLY: approx_l1_norm_sd,&
                                              crazy_rotations,&
                                              direct_mini,&
                                              jacobi_cg_edf_ls,&
                                              jacobi_rotations,&
                                              rotate_orbitals,&
                                              scdm_qrfact,&
                                              zij_matrix
   USE qs_matrix_pools,                 ONLY: mpools_get
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE string_utilities,                ONLY: xstring
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_loc_methods'

! *** Public ***
   PUBLIC :: qs_print_cubes, centers_spreads_berry, centers_second_moments_loc, &
             centers_second_moments_berry, jacobi_rotation_pipek, optimize_loc_berry, &
             optimize_loc_pipek

CONTAINS

! **************************************************************************************************
!> \brief Calculate and optimize the spread functional as calculated from
!>       the selected mos  and by the definition using the berry phase
!>       as given by silvestrelli et al
!>       If required the centers and the spreads for each mos selected
!>       are calculated from z_ij and printed to file.
!>       The centers files is appended if already exists
!> \param method indicates localization algorithm
!> \param qs_loc_env new environment for the localization calculations
!> \param vectors selected mos to be localized
!> \param op_sm_set sparse matrices containing the integrals of the kind <mi e{iGr} nu>
!> \param zij_fm_set set of full matrix of size nmoloc x nmoloc, will contain the z_ij numbers
!>                    as defined by Silvestrelli et al
!> \param para_env ...
!> \param cell ...
!> \param weights ...
!> \param ispin ...
!> \param print_loc_section ...
!> \param restricted ...
!> \param nextra ...
!> \param nmo ...
!> \param vectors_2 ...
!> \param guess_mos ...
!> \par History
!>      04.2005 created [MI]
!> \author MI
!> \note
!>       This definition need the use of complex numbers, therefore the
!>       optimization routines are specific for this case
!>       The file for the centers and the spreads have a xyz format
! **************************************************************************************************
   SUBROUTINE optimize_loc_berry(method, qs_loc_env, vectors, op_sm_set, &
                                 zij_fm_set, para_env, cell, weights, ispin, print_loc_section, &
                                 restricted, nextra, nmo, vectors_2, guess_mos)

      INTEGER, INTENT(IN)                                :: method
      TYPE(qs_loc_env_type), POINTER                     :: qs_loc_env
      TYPE(cp_fm_type), INTENT(IN)                       :: vectors
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: op_sm_set
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: zij_fm_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cell_type), POINTER                           :: cell
      REAL(dp), DIMENSION(:)                             :: weights
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(section_vals_type), POINTER                   :: print_loc_section
      INTEGER                                            :: restricted
      INTEGER, INTENT(IN), OPTIONAL                      :: nextra, nmo
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: vectors_2, guess_mos

      CHARACTER(len=*), PARAMETER :: routineN = 'optimize_loc_berry'

      INTEGER                                            :: handle, max_iter, nao, nmoloc, out_each, &
                                                            output_unit, sweeps
      LOGICAL                                            :: converged, crazy_use_diag, &
                                                            do_jacobi_refinement, my_do_mixed
      REAL(dp)                                           :: crazy_scale, eps_localization, &
                                                            max_crazy_angle, start_time, &
                                                            target_time
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      output_unit = cp_print_key_unit_nr(logger, print_loc_section, "PROGRAM_RUN_INFO", &
                                         extension=".locInfo")

      ! get rows and cols of the input
      CALL cp_fm_get_info(vectors, nrow_global=nao, ncol_global=nmoloc)

      CALL zij_matrix(vectors, op_sm_set, zij_fm_set)

      max_iter = qs_loc_env%localized_wfn_control%max_iter
      max_crazy_angle = qs_loc_env%localized_wfn_control%max_crazy_angle
      crazy_use_diag = qs_loc_env%localized_wfn_control%crazy_use_diag
      crazy_scale = qs_loc_env%localized_wfn_control%crazy_scale
      eps_localization = qs_loc_env%localized_wfn_control%eps_localization
      out_each = qs_loc_env%localized_wfn_control%out_each
      target_time = qs_loc_env%target_time
      start_time = qs_loc_env%start_time
      do_jacobi_refinement = qs_loc_env%localized_wfn_control%jacobi_refinement
      my_do_mixed = qs_loc_env%localized_wfn_control%do_mixed

      CALL centers_spreads_berry(qs_loc_env, zij_fm_set, nmoloc, cell, weights, &
                                 ispin, print_loc_section, only_initial_out=.TRUE.)

      SELECT CASE (method)
      CASE (do_loc_jacobi)
         CALL jacobi_rotations(weights, zij_fm_set, vectors, para_env, max_iter=max_iter, &
                               eps_localization=eps_localization, sweeps=sweeps, &
                               out_each=out_each, target_time=target_time, start_time=start_time, &
                               restricted=restricted)
      CASE (do_loc_gapo)
         IF (my_do_mixed) THEN
            IF (nextra > 0) THEN
               IF (PRESENT(guess_mos)) THEN
                  CALL jacobi_cg_edf_ls(para_env, weights, zij_fm_set, vectors, max_iter, &
                                        eps_localization, sweeps, out_each, nextra, &
                                        qs_loc_env%localized_wfn_control%do_cg_po, &
                                        nmo=nmo, vectors_2=vectors_2, mos_guess=guess_mos)
               ELSE
                  CALL jacobi_cg_edf_ls(para_env, weights, zij_fm_set, vectors, max_iter, &
                                        eps_localization, sweeps, out_each, nextra, &
                                        qs_loc_env%localized_wfn_control%do_cg_po, &
                                        nmo=nmo, vectors_2=vectors_2)
               END IF
            ELSE
               CALL jacobi_cg_edf_ls(para_env, weights, zij_fm_set, vectors, max_iter, &
                                     eps_localization, sweeps, out_each, 0, &
                                     qs_loc_env%localized_wfn_control%do_cg_po)
            END IF
         ELSE
            CPABORT("GAPO works only with STATES MIXED")
         END IF
      CASE (do_loc_scdm)
         ! Decomposition
         CALL scdm_qrfact(vectors)
         ! Calculation of Zij
         CALL zij_matrix(vectors, op_sm_set, zij_fm_set)
         IF (do_jacobi_refinement) THEN
            ! Intermediate spread and centers
            CALL centers_spreads_berry(qs_loc_env, zij_fm_set, nmoloc, cell, weights, &
                                       ispin, print_loc_section, only_initial_out=.TRUE.)
            CALL jacobi_rotations(weights, zij_fm_set, vectors, para_env, max_iter=max_iter, &
                                  eps_localization=eps_localization, sweeps=sweeps, &
                                  out_each=out_each, target_time=target_time, start_time=start_time, &
                                  restricted=restricted)
         END IF
      CASE (do_loc_crazy)
         CALL crazy_rotations(weights, zij_fm_set, vectors, max_iter=max_iter, max_crazy_angle=max_crazy_angle, &
                              crazy_scale=crazy_scale, crazy_use_diag=crazy_use_diag, &
                              eps_localization=eps_localization, iterations=sweeps, converged=converged)
         ! Possibly fallback to jacobi if the crazy rotation fails
         IF (.NOT. converged) THEN
            IF (qs_loc_env%localized_wfn_control%jacobi_fallback) THEN
               IF (output_unit > 0) WRITE (output_unit, "(T4,A,I6,/,T4,A)") &
                  " Crazy Wannier localization not converged after ", sweeps, &
                  " iterations, switching to jacobi rotations"
               CALL jacobi_rotations(weights, zij_fm_set, vectors, para_env, max_iter=max_iter, &
                                     eps_localization=eps_localization, sweeps=sweeps, &
                                     out_each=out_each, target_time=target_time, start_time=start_time, &
                                     restricted=restricted)
            ELSE
               IF (output_unit > 0) WRITE (output_unit, "(T4,A,I6,/,T4,A)") &
                  " Crazy Wannier localization not converged after ", sweeps, &
                  " iterations, and jacobi_fallback switched off "
            END IF
         END IF
      CASE (do_loc_direct)
         CALL direct_mini(weights, zij_fm_set, vectors, max_iter=max_iter, &
                          eps_localization=eps_localization, iterations=sweeps)
      CASE (do_loc_l1_norm_sd)
         IF (.NOT. cell%orthorhombic) THEN
            CPABORT("Non-orthorhombic cell with the selected method NYI")
         ELSE
            CALL approx_l1_norm_sd(vectors, max_iter, eps_localization, converged, sweeps)
            ! here we need to set zij for the computation of the centers and spreads
            CALL zij_matrix(vectors, op_sm_set, zij_fm_set)
         END IF
      CASE (do_loc_none)
         IF (output_unit > 0) THEN
            WRITE (output_unit, '(T4,A,I6,A)') " No MOS localization applied "
         END IF
      CASE DEFAULT
         CPABORT("Unknown localization method")
      END SELECT
      IF (output_unit > 0) THEN
         IF (sweeps <= max_iter) WRITE (output_unit, '(T4,A,I3,A,I6,A)') " Localization  for spin ", ispin, &
            " converged in ", sweeps, " iterations"
      END IF

      CALL centers_spreads_berry(qs_loc_env, zij_fm_set, nmoloc, cell, weights, &
                                 ispin, print_loc_section)

      CALL cp_print_key_finished_output(output_unit, logger, print_loc_section, "PROGRAM_RUN_INFO")

      CALL timestop(handle)

   END SUBROUTINE optimize_loc_berry

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param method ...
!> \param qs_loc_env ...
!> \param vectors ...
!> \param zij_fm_set ...
!> \param ispin ...
!> \param print_loc_section ...
!> \par History
!>      04.2005 created [MI]
!> \author MI
! **************************************************************************************************
   SUBROUTINE optimize_loc_pipek(qs_env, method, qs_loc_env, vectors, zij_fm_set, &
                                 ispin, print_loc_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: method
      TYPE(qs_loc_env_type), POINTER                     :: qs_loc_env
      TYPE(cp_fm_type), INTENT(IN)                       :: vectors
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: zij_fm_set
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(section_vals_type), POINTER                   :: print_loc_section

      CHARACTER(len=*), PARAMETER :: routineN = 'optimize_loc_pipek'

      INTEGER                                            :: handle, iatom, isgf, ldz, nao, natom, &
                                                            ncol, nmoloc, output_unit, sweeps
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, last_sgf, nsgf
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_ao_fm_pools
      TYPE(cp_fm_type)                                   :: opvec
      TYPE(cp_fm_type), POINTER                          :: ov_fm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      output_unit = cp_print_key_unit_nr(logger, print_loc_section, "PROGRAM_RUN_INFO", &
                                         extension=".locInfo")

      NULLIFY (particle_set)
      ! get rows and cols of the input
      CALL cp_fm_get_info(vectors, nrow_global=nao, ncol_global=nmoloc)

      ! replicate the input kind of matrix
      CALL cp_fm_create(opvec, vectors%matrix_struct)
      CALL cp_fm_set_all(opvec, 0.0_dp)

      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set)

      natom = SIZE(particle_set, 1)
      ALLOCATE (first_sgf(natom))
      ALLOCATE (last_sgf(natom))
      ALLOCATE (nsgf(natom))

      !   construction of
      CALL get_particle_set(particle_set, qs_kind_set, &
                            first_sgf=first_sgf, last_sgf=last_sgf, nsgf=nsgf)

      !   Copy the overlap sparse matrix in a full matrix
      CALL mpools_get(qs_env%mpools, ao_ao_fm_pools=ao_ao_fm_pools)
      ALLOCATE (ov_fm)
      CALL fm_pool_create_fm(ao_ao_fm_pools(1)%pool, ov_fm, name=" ")
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, ov_fm)

      !   Compute zij here
      DO iatom = 1, natom
         CALL cp_fm_set_all(zij_fm_set(iatom, 1), 0.0_dp)
         CALL cp_fm_get_info(zij_fm_set(iatom, 1), ncol_global=ldz)
         isgf = first_sgf(iatom)
         ncol = nsgf(iatom)

         ! multiply fmxfm, using only part of the ao : Ct x S
         CALL parallel_gemm('N', 'N', nao, nmoloc, nao, 1.0_dp, ov_fm, vectors, 0.0_dp, opvec, &
                            a_first_col=1, a_first_row=1, b_first_col=1, b_first_row=1)

         CALL parallel_gemm('T', 'N', nmoloc, nmoloc, ncol, 0.5_dp, vectors, opvec, &
                            0.0_dp, zij_fm_set(iatom, 1), &
                            a_first_col=1, a_first_row=isgf, b_first_col=1, b_first_row=isgf)

         CALL parallel_gemm('N', 'N', nao, nmoloc, ncol, 1.0_dp, ov_fm, vectors, 0.0_dp, opvec, &
                            a_first_col=isgf, a_first_row=1, b_first_col=1, b_first_row=isgf)

         CALL parallel_gemm('T', 'N', nmoloc, nmoloc, nao, 0.5_dp, vectors, opvec, &
                            1.0_dp, zij_fm_set(iatom, 1), &
                            a_first_col=1, a_first_row=1, b_first_col=1, b_first_row=1)

      END DO ! iatom

      !   And now perform the optimization and rotate the orbitals
      SELECT CASE (method)
      CASE (do_loc_jacobi)
         CALL jacobi_rotation_pipek(zij_fm_set, vectors, sweeps)
      CASE (do_loc_gapo)
         CPABORT("GAPO and Pipek not implemented.")
      CASE (do_loc_crazy)
         CPABORT("Crazy and Pipek not implemented.")
      CASE (do_loc_l1_norm_sd)
         CPABORT("L1 norm and Pipek not implemented.")
      CASE (do_loc_direct)
         CPABORT("Direct and Pipek not implemented.")
      CASE (do_loc_none)
         IF (output_unit > 0) WRITE (output_unit, '(A,I6,A)') " No MOS localization applied "
      CASE DEFAULT
         CPABORT("Unknown localization method")
      END SELECT

      IF (output_unit > 0) WRITE (output_unit, '(A,I3,A,I6,A)') " Localization  for spin ", ispin, &
         " converged in ", sweeps, " iterations"

      CALL centers_spreads_pipek(qs_loc_env, zij_fm_set, particle_set, ispin, &
                                 print_loc_section)

      DEALLOCATE (first_sgf, last_sgf, nsgf)

      CALL cp_fm_release(opvec)
      CALL cp_print_key_finished_output(output_unit, logger, print_loc_section, "PROGRAM_RUN_INFO")

      CALL timestop(handle)

   END SUBROUTINE optimize_loc_pipek

! **************************************************************************************************
!> \brief 2by2 rotation for the pipek operator
!>       in this case the z_ij numbers are reals
!> \param zij_fm_set ...
!> \param vectors ...
!> \param sweeps ...
!> \par History
!>       05-2005 created
!> \author MI
! **************************************************************************************************
   SUBROUTINE jacobi_rotation_pipek(zij_fm_set, vectors, sweeps)

      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: zij_fm_set
      TYPE(cp_fm_type), INTENT(IN)                       :: vectors
      INTEGER                                            :: sweeps

      INTEGER                                            :: iatom, istate, jstate, natom, nstate
      REAL(KIND=dp)                                      :: aij, bij, ct, mii, mij, mjj, ratio, st, &
                                                            theta, tolerance
      TYPE(cp_fm_type)                                   :: rmat

      CALL cp_fm_create(rmat, zij_fm_set(1, 1)%matrix_struct)
      CALL cp_fm_set_all(rmat, 0.0_dp, 1.0_dp)

      CALL cp_fm_get_info(rmat, nrow_global=nstate)
      tolerance = 1.0e10_dp
      sweeps = 0
      natom = SIZE(zij_fm_set, 1)
      ! do jacobi sweeps until converged
      DO WHILE (tolerance >= 1.0e-4_dp)
         sweeps = sweeps + 1
         DO istate = 1, nstate
            DO jstate = istate + 1, nstate
               aij = 0.0_dp
               bij = 0.0_dp
               DO iatom = 1, natom
                  CALL cp_fm_get_element(zij_fm_set(iatom, 1), istate, istate, mii)
                  CALL cp_fm_get_element(zij_fm_set(iatom, 1), istate, jstate, mij)
                  CALL cp_fm_get_element(zij_fm_set(iatom, 1), jstate, jstate, mjj)
                  aij = aij + mij*(mii - mjj)
                  bij = bij + mij*mij - 0.25_dp*(mii - mjj)*(mii - mjj)
               END DO
               IF (ABS(bij) > 1.E-10_dp) THEN
                  ratio = -aij/bij
                  theta = 0.25_dp*ATAN(ratio)
               ELSE
                  bij = 0.0_dp
                  theta = 0.0_dp
               END IF
               ! Check max or min
               ! To minimize the spread
               IF (theta > pi*0.5_dp) THEN
                  theta = theta - pi*0.25_dp
               ELSE IF (theta < -pi*0.5_dp) THEN
                  theta = theta + pi*0.25_dp
               END IF

               ct = COS(theta)
               st = SIN(theta)

               CALL cp_fm_rot_cols(rmat, istate, jstate, ct, st)

               DO iatom = 1, natom
                  CALL cp_fm_rot_cols(zij_fm_set(iatom, 1), istate, jstate, ct, st)
                  CALL cp_fm_rot_rows(zij_fm_set(iatom, 1), istate, jstate, ct, st)
               END DO
            END DO
         END DO
         CALL check_tolerance_real(zij_fm_set, tolerance)
      END DO

      CALL rotate_orbitals(rmat, vectors)
      CALL cp_fm_release(rmat)

   END SUBROUTINE jacobi_rotation_pipek

! **************************************************************************************************
!> \brief ...
!> \param zij_fm_set ...
!> \param tolerance ...
!> \par History
!>      04.2005 created [MI]
!> \author MI
! **************************************************************************************************
   SUBROUTINE check_tolerance_real(zij_fm_set, tolerance)

      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: zij_fm_set
      REAL(dp), INTENT(OUT)                              :: tolerance

      INTEGER                                            :: iatom, istate, jstate, natom, &
                                                            ncol_local, nrow_global, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(dp)                                           :: grad_ij, zii, zij, zjj
      REAL(dp), DIMENSION(:, :), POINTER                 :: diag
      TYPE(cp_fm_type)                                   :: force

      CALL cp_fm_create(force, zij_fm_set(1, 1)%matrix_struct)
      CALL cp_fm_set_all(force, 0._dp)

      NULLIFY (diag, col_indices, row_indices)
      natom = SIZE(zij_fm_set, 1)
      CALL cp_fm_get_info(zij_fm_set(1, 1), nrow_local=nrow_local, &
                          ncol_local=ncol_local, nrow_global=nrow_global, &
                          row_indices=row_indices, col_indices=col_indices)
      ALLOCATE (diag(nrow_global, natom))

      DO iatom = 1, natom
         DO istate = 1, nrow_global
            CALL cp_fm_get_element(zij_fm_set(iatom, 1), istate, istate, diag(istate, iatom))
         END DO
      END DO

      DO istate = 1, nrow_local
         DO jstate = 1, ncol_local
            grad_ij = 0.0_dp
            DO iatom = 1, natom
               zii = diag(row_indices(istate), iatom)
               zjj = diag(col_indices(jstate), iatom)
               zij = zij_fm_set(iatom, 1)%local_data(istate, jstate)
               grad_ij = grad_ij + 4.0_dp*zij*(zjj - zii)
            END DO
            force%local_data(istate, jstate) = grad_ij
         END DO
      END DO

      DEALLOCATE (diag)

      CALL cp_fm_maxabsval(force, tolerance)
      CALL cp_fm_release(force)

   END SUBROUTINE check_tolerance_real
! **************************************************************************************************
!> \brief ...
!> \param qs_loc_env ...
!> \param zij ...
!> \param nmoloc ...
!> \param cell ...
!> \param weights ...
!> \param ispin ...
!> \param print_loc_section ...
!> \param only_initial_out ...
!> \par History
!>      04.2005 created [MI]
!> \author MI
! **************************************************************************************************
   SUBROUTINE centers_spreads_berry(qs_loc_env, zij, nmoloc, cell, weights, ispin, &
                                    print_loc_section, only_initial_out)
      TYPE(qs_loc_env_type), POINTER                     :: qs_loc_env
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: zij
      INTEGER, INTENT(IN)                                :: nmoloc
      TYPE(cell_type), POINTER                           :: cell
      REAL(dp), DIMENSION(:)                             :: weights
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(section_vals_type), POINTER                   :: print_loc_section
      LOGICAL, INTENT(IN), OPTIONAL                      :: only_initial_out

      CHARACTER(len=default_path_length)                 :: file_tmp, iter
      COMPLEX(KIND=dp)                                   :: z
      INTEGER                                            :: idir, istate, jdir, nstates, &
                                                            output_unit, unit_out_s
      LOGICAL                                            :: my_only_init
      REAL(dp)                                           :: avg_spread_ii, spread_i, spread_ii, &
                                                            sum_spread_i, sum_spread_ii
      REAL(dp), DIMENSION(3)                             :: c, c2, cpbc
      REAL(dp), DIMENSION(:, :), POINTER                 :: centers
      REAL(KIND=dp)                                      :: imagpart, realpart
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (centers, logger, print_key)
      logger => cp_get_default_logger()
      my_only_init = .FALSE.
      IF (PRESENT(only_initial_out)) my_only_init = only_initial_out

      file_tmp = TRIM(qs_loc_env%tag_mo)//"_spreads_s"//TRIM(ADJUSTL(cp_to_string(ispin)))
      output_unit = cp_print_key_unit_nr(logger, print_loc_section, "PROGRAM_RUN_INFO", &
                                         extension=".locInfo")
      unit_out_s = cp_print_key_unit_nr(logger, print_loc_section, "WANNIER_SPREADS", &
                                        middle_name=file_tmp, extension=".data")
      iter = cp_iter_string(logger%iter_info)
      IF (unit_out_s > 0 .AND. .NOT. my_only_init) WRITE (unit_out_s, '(i6,/,A)') nmoloc, TRIM(iter)

      CALL cp_fm_get_info(zij(1, 1), nrow_global=nstates)
      CPASSERT(nstates >= nmoloc)

      centers => qs_loc_env%localized_wfn_control%centers_set(ispin)%array
      CPASSERT(ASSOCIATED(centers))
      CPASSERT(SIZE(centers, 2) == nmoloc)
      sum_spread_i = 0.0_dp
      sum_spread_ii = 0.0_dp
      avg_spread_ii = 0.0_dp
      DO istate = 1, nmoloc
         c = 0.0_dp
         c2 = 0.0_dp
         spread_i = 0.0_dp
         spread_ii = 0.0_dp
         DO jdir = 1, SIZE(zij, 2)
            CALL cp_fm_get_element(zij(1, jdir), istate, istate, realpart)
            CALL cp_fm_get_element(zij(2, jdir), istate, istate, imagpart)
            z = CMPLX(realpart, imagpart, dp)
            spread_i = spread_i - weights(jdir)* &
                       LOG(realpart*realpart + imagpart*imagpart)/twopi/twopi
            spread_ii = spread_ii + weights(jdir)* &
                        (1.0_dp - (realpart*realpart + imagpart*imagpart))/twopi/twopi
            IF (jdir < 4) THEN
               DO idir = 1, 3
                  c(idir) = c(idir) + &
                            (cell%hmat(idir, jdir)/twopi)*AIMAG(LOG(z))
               END DO
            END IF
         END DO
         cpbc = pbc(c, cell)
         centers(1:3, istate) = cpbc(1:3)
         centers(4, istate) = spread_i
         centers(5, istate) = spread_ii
         sum_spread_i = sum_spread_i + centers(4, istate)
         sum_spread_ii = sum_spread_ii + centers(5, istate)
         IF (unit_out_s > 0 .AND. .NOT. my_only_init) WRITE (unit_out_s, '(I6,2F16.8)') istate, centers(4:5, istate)
      END DO
      avg_spread_ii = sum_spread_ii/REAL(nmoloc, KIND=dp)

      ! Print of wannier centers
      print_key => section_vals_get_subs_vals(print_loc_section, "WANNIER_CENTERS")
      IF (.NOT. my_only_init) CALL print_wannier_centers(qs_loc_env, print_key, centers, logger, ispin)

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(T4, A, 2x, 2A26,/,T23, A28)') " Spread Functional ", "sum_in -w_i ln(|z_in|^2)", &
            "sum_in w_i(1-|z_in|^2)", "sum_in w_i(1-|z_in|^2)/n"
         IF (my_only_init) THEN
            WRITE (output_unit, '(T4,A,T38,2F20.10,/,T38,F20.10)') " Initial Spread (Berry) : ", &
               sum_spread_i, sum_spread_ii, avg_spread_ii
         ELSE
            WRITE (output_unit, '(T4,A,T38,2F20.10,/,T38,F20.10)') " Total Spread (Berry) : ", &
               sum_spread_i, sum_spread_ii, avg_spread_ii
         END IF
      END IF

      IF (unit_out_s > 0 .AND. .NOT. my_only_init) WRITE (unit_out_s, '(A,2F16.10)') " Total ", sum_spread_i, sum_spread_ii

      CALL cp_print_key_finished_output(unit_out_s, logger, print_loc_section, "WANNIER_SPREADS")
      CALL cp_print_key_finished_output(output_unit, logger, print_loc_section, "PROGRAM_RUN_INFO")

   END SUBROUTINE centers_spreads_berry

! **************************************************************************************************
!> \brief define and print the centers and spread
!>       when the pipek operator is used
!> \param qs_loc_env ...
!> \param zij_fm_set matrix elements that define the populations on atoms
!> \param particle_set ...
!> \param ispin spin 1 or 2
!> \param print_loc_section ...
!> \par History
!>      05.2005 created [MI]
! **************************************************************************************************
   SUBROUTINE centers_spreads_pipek(qs_loc_env, zij_fm_set, particle_set, ispin, &
                                    print_loc_section)

      TYPE(qs_loc_env_type), POINTER                     :: qs_loc_env
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: zij_fm_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(section_vals_type), POINTER                   :: print_loc_section

      CHARACTER(len=default_path_length)                 :: file_tmp, iter
      INTEGER                                            :: iatom, istate, natom, nstate, unit_out_s
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_state
      REAL(dp)                                           :: r(3)
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: Qii, ziimax
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: diag
      REAL(dp), DIMENSION(:, :), POINTER                 :: centers
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (centers, logger, print_key)
      logger => cp_get_default_logger()

      CALL cp_fm_get_info(zij_fm_set(1, 1), nrow_global=nstate)
      natom = SIZE(zij_fm_set, 1)

      centers => qs_loc_env%localized_wfn_control%centers_set(ispin)%array
      CPASSERT(ASSOCIATED(centers))
      CPASSERT(SIZE(centers, 2) == nstate)

      file_tmp = TRIM(qs_loc_env%tag_mo)//"_spreads_s"//TRIM(ADJUSTL(cp_to_string(ispin)))
      unit_out_s = cp_print_key_unit_nr(logger, print_loc_section, "WANNIER_SPREADS", &
                                        middle_name=file_tmp, extension=".data", log_filename=.FALSE.)
      iter = cp_iter_string(logger%iter_info)
      IF (unit_out_s > 0) WRITE (unit_out_s, '(i6,/,A)') TRIM(iter)

      ALLOCATE (atom_of_state(nstate))
      atom_of_state = 0

      ALLOCATE (diag(nstate, natom))
      diag = 0.0_dp

      DO iatom = 1, natom
         DO istate = 1, nstate
            CALL cp_fm_get_element(zij_fm_set(iatom, 1), istate, istate, diag(istate, iatom))
         END DO
      END DO

      ALLOCATE (Qii(nstate), ziimax(nstate))
      ziimax = 0.0_dp
      Qii = 0.0_dp

      DO iatom = 1, natom
         DO istate = 1, nstate
            Qii(istate) = Qii(istate) + diag(istate, iatom)*diag(istate, iatom)
            IF (ABS(diag(istate, iatom)) > ziimax(istate)) THEN
               ziimax(istate) = ABS(diag(istate, iatom))
               atom_of_state(istate) = iatom
            END IF
         END DO
      END DO

      DO istate = 1, nstate
         iatom = atom_of_state(istate)
         r(1:3) = particle_set(iatom)%r(1:3)
         centers(1:3, istate) = r(1:3)
         centers(4, istate) = 1.0_dp/Qii(istate)
         IF (unit_out_s > 0) WRITE (unit_out_s, '(I6,F16.8)') istate, angstrom*centers(4, istate)
      END DO

      ! Print the wannier centers
      print_key => section_vals_get_subs_vals(print_loc_section, "WANNIER_CENTERS")
      CALL print_wannier_centers(qs_loc_env, print_key, centers, logger, ispin)

      CALL cp_print_key_finished_output(unit_out_s, logger, print_loc_section, "WANNIER_SPREADS")

      DEALLOCATE (Qii, ziimax, atom_of_state, diag)

   END SUBROUTINE centers_spreads_pipek

! **************************************************************************************************
!> \brief write the cube files for a set of selected states
!> \param qs_env ...
!> \param mo_coeff set mos from which the states to be printed are extracted
!> \param nstates number of states to be printed
!> \param state_list list of the indexes of the states to be printed
!> \param centers centers and spread, all=0 if they hva not been calculated
!> \param print_key ...
!> \param root initial part of the cube file names
!> \param ispin ...
!> \param idir ...
!> \param state0 ...
!> \param file_position ...
!> \par History
!>      08.2005 created [MI]
!> \author MI
!> \note
!>      This routine should not be in this module
! **************************************************************************************************
   SUBROUTINE qs_print_cubes(qs_env, mo_coeff, nstates, state_list, centers, &
                             print_key, root, ispin, idir, state0, file_position)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: mo_coeff
      INTEGER, INTENT(IN)                                :: nstates
      INTEGER, DIMENSION(:), POINTER                     :: state_list
      REAL(dp), DIMENSION(:, :), POINTER                 :: centers
      TYPE(section_vals_type), POINTER                   :: print_key
      CHARACTER(LEN=*)                                   :: root
      INTEGER, INTENT(IN), OPTIONAL                      :: ispin, idir
      INTEGER, OPTIONAL                                  :: state0
      CHARACTER(LEN=default_string_length), OPTIONAL     :: file_position

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_print_cubes'
      CHARACTER, DIMENSION(3), PARAMETER                 :: labels = (/'x', 'y', 'z'/)

      CHARACTER                                          :: label
      CHARACTER(LEN=default_path_length)                 :: file_tmp, filename, my_pos
      CHARACTER(LEN=default_string_length)               :: title
      INTEGER                                            :: handle, ia, ie, istate, ivector, &
                                                            my_ispin, my_state0, unit_out_c
      LOGICAL                                            :: add_idir, add_spin, mpi_io
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: wf_g, wf_r
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CALL timeset(routineN, handle)
      NULLIFY (auxbas_pw_pool, pw_env, logger)
      logger => cp_get_default_logger()

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, subsys=subsys)
      CALL qs_subsys_get(subsys, particles=particles)

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      CALL auxbas_pw_pool%create_pw(wf_r, &
                                    use_data=REALDATA3D, &
                                    in_space=REALSPACE)
      CALL auxbas_pw_pool%create_pw(wf_g, &
                                    use_data=COMPLEXDATA1D, &
                                    in_space=RECIPROCALSPACE)

      my_state0 = 0
      IF (PRESENT(state0)) my_state0 = state0

      add_spin = .FALSE.
      my_ispin = 1
      IF (PRESENT(ispin)) THEN
         add_spin = .TRUE.
         my_ispin = ispin
      END IF
      add_idir = .FALSE.
      IF (PRESENT(idir)) THEN
         add_idir = .TRUE.
         label = labels(idir)
      END IF

      my_pos = "REWIND"
      IF (PRESENT(file_position)) THEN
         my_pos = file_position
      END IF

      DO istate = 1, nstates
         ivector = state_list(istate) - my_state0
         CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, cell=cell, &
                         dft_control=dft_control, particle_set=particle_set, pw_env=pw_env)

         CALL calculate_wavefunction(mo_coeff, ivector, wf_r, wf_g, atomic_kind_set, &
                                     qs_kind_set, cell, dft_control, particle_set, pw_env)

         ! Formatting the middle part of the name
         ivector = state_list(istate)
         CALL xstring(root, ia, ie)
         IF (add_idir) THEN
            filename = root(ia:ie)//"_"//label//"_w"//TRIM(ADJUSTL(cp_to_string(ivector)))
         ELSE
            filename = root(ia:ie)//"_w"//TRIM(ADJUSTL(cp_to_string(ivector)))
         END IF
         IF (add_spin) THEN
            file_tmp = filename
            CALL xstring(file_tmp, ia, ie)
            filename = file_tmp(ia:ie)//"_s"//TRIM(ADJUSTL(cp_to_string(ispin)))
         END IF

         ! Using the print_key tools to open the file with the proper name
         mpi_io = .TRUE.
         unit_out_c = cp_print_key_unit_nr(logger, print_key, "", middle_name=filename, &
                                           extension=".cube", file_position=my_pos, log_filename=.FALSE., &
                                           mpi_io=mpi_io)
         IF (SIZE(centers, 1) == 6) THEN
            WRITE (title, '(A7,I5.5,A2,I1.1,A1,6F10.4)') "WFN ", ivector, "_s", my_ispin, " ", &
               centers(1:3, istate)*angstrom, centers(4:6, istate)*angstrom
         ELSE
            WRITE (title, '(A7,I5.5,A2,I1.1,A1,3F10.4)') "WFN ", ivector, "_s", my_ispin, " ", &
               centers(1:3, istate)*angstrom
         END IF
         CALL cp_pw_to_cube(wf_r, unit_out_c, title, &
                            particles=particles, &
                            stride=section_get_ivals(print_key, "STRIDE"), mpi_io=mpi_io)
         CALL cp_print_key_finished_output(unit_out_c, logger, print_key, "", mpi_io=mpi_io)
      END DO

      CALL auxbas_pw_pool%give_back_pw(wf_r)
      CALL auxbas_pw_pool%give_back_pw(wf_g)
      CALL timestop(handle)
   END SUBROUTINE qs_print_cubes

! **************************************************************************************************
!> \brief Prints wannier centers
!> \param qs_loc_env ...
!> \param print_key ...
!> \param center ...
!> \param logger ...
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE print_wannier_centers(qs_loc_env, print_key, center, logger, ispin)
      TYPE(qs_loc_env_type), POINTER                     :: qs_loc_env
      TYPE(section_vals_type), POINTER                   :: print_key
      REAL(KIND=dp), INTENT(IN)                          :: center(:, :)
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(default_string_length)                   :: iter, unit_str
      CHARACTER(LEN=default_string_length)               :: my_ext, my_form
      INTEGER                                            :: iunit, l, nstates
      LOGICAL                                            :: first_time, init_traj
      REAL(KIND=dp)                                      :: unit_conv

      nstates = SIZE(center, 2)
      my_form = "formatted"
      my_ext = ".data"
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key, first_time=first_time) &
                , cp_p_file)) THEN
         ! Find out if we want to print IONS+CENTERS or ONLY CENTERS..
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key, "/IONS+CENTERS") &
                   , cp_p_file)) THEN
            CALL get_output_format(print_key, my_form=my_form, my_ext=my_ext)
         END IF
         IF (first_time .OR. (.NOT. qs_loc_env%first_time)) THEN
            iunit = cp_print_key_unit_nr(logger, print_key, "", extension=my_ext, file_form=my_form, &
                                         middle_name=TRIM(qs_loc_env%tag_mo)//"_centers_s"//TRIM(ADJUSTL(cp_to_string(ispin))), &
                                         log_filename=.FALSE., on_file=.TRUE., is_new_file=init_traj)
            IF (iunit > 0) THEN
               ! Gather units of measure for output (if available)
               CALL section_vals_val_get(print_key, "UNIT", c_val=unit_str)
               unit_conv = cp_unit_from_cp2k(1.0_dp, TRIM(unit_str))

               IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key, "/IONS+CENTERS"), cp_p_file)) THEN
                  ! Different possible formats
                  CALL print_wannier_traj(qs_loc_env, print_key, center, iunit, init_traj, unit_conv)
               ELSE
                  ! Default print format
                  iter = cp_iter_string(logger%iter_info)
                  WRITE (iunit, '(i6,/,a)') nstates, TRIM(iter)
                  DO l = 1, nstates
                     WRITE (iunit, '(A,4F16.8)') "X", unit_conv*center(1:4, l)
                  END DO
               END IF
            END IF
            CALL cp_print_key_finished_output(iunit, logger, print_key, on_file=.TRUE.)
         END IF
      END IF
   END SUBROUTINE print_wannier_centers

! **************************************************************************************************
!> \brief computes spread and centers
!> \param qs_loc_env ...
!> \param print_key ...
!> \param center ...
!> \param iunit ...
!> \param init_traj ...
!> \param unit_conv ...
! **************************************************************************************************
   SUBROUTINE print_wannier_traj(qs_loc_env, print_key, center, iunit, init_traj, unit_conv)
      TYPE(qs_loc_env_type), POINTER                     :: qs_loc_env
      TYPE(section_vals_type), POINTER                   :: print_key
      REAL(KIND=dp), INTENT(IN)                          :: center(:, :)
      INTEGER, INTENT(IN)                                :: iunit
      LOGICAL, INTENT(IN)                                :: init_traj
      REAL(KIND=dp), INTENT(IN)                          :: unit_conv

      CHARACTER(len=default_string_length)               :: iter, remark1, remark2, title
      INTEGER                                            :: i, iskip, natom, ntot, outformat
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (particle_set, atomic_kind_set, atomic_kind, logger)
      logger => cp_get_default_logger()
      natom = SIZE(qs_loc_env%particle_set)
      ntot = natom + SIZE(center, 2)
      CALL allocate_particle_set(particle_set, ntot)
      ALLOCATE (atomic_kind_set(1))
      atomic_kind => atomic_kind_set(1)
      CALL set_atomic_kind(atomic_kind=atomic_kind, kind_number=0, &
                           name="X", element_symbol="X", mass=0.0_dp)
      ! Particles
      DO i = 1, natom
         particle_set(i)%atomic_kind => qs_loc_env%particle_set(i)%atomic_kind
         particle_set(i)%r = pbc(qs_loc_env%particle_set(i)%r, qs_loc_env%cell)
      END DO
      ! Wannier Centers
      DO i = natom + 1, ntot
         particle_set(i)%atomic_kind => atomic_kind
         particle_set(i)%r = pbc(center(1:3, i - natom), qs_loc_env%cell)
      END DO
      ! Dump the structure
      CALL section_vals_val_get(print_key, "FORMAT", i_val=outformat)

      ! Header file
      SELECT CASE (outformat)
      CASE (dump_dcd, dump_dcd_aligned_cell)
         IF (init_traj) THEN
            !Lets write the header for the coordinate dcd
            ! Note (TL) : even the new DCD format is unfortunately too poor
            !             for our capabilities.. for example here the printing
            !             of the geometry could be nested inside several iteration
            !             levels.. this cannot be exactly reproduce with DCD.
            !             Just as a compromise let's pick-up the value of the MD iteration
            !             level. In any case this is not any sensible information for the standard..
            iskip = section_get_ival(print_key, "EACH%MD")
            WRITE (iunit) "CORD", 0, -1, iskip, &
               0, 0, 0, 0, 0, 0, REAL(0, KIND=sp), 1, 0, 0, 0, 0, 0, 0, 0, 0, 24
            remark1 = "REMARK FILETYPE CORD DCD GENERATED BY CP2K"
            remark2 = "REMARK Support new DCD format with cell information"
            WRITE (iunit) 2, remark1, remark2
            WRITE (iunit) ntot
            CALL m_flush(iunit)
         END IF
      CASE (dump_xmol)
         iter = cp_iter_string(logger%iter_info)
         WRITE (UNIT=title, FMT="(A)") " Particles+Wannier centers. Iteration:"//TRIM(iter)
      CASE DEFAULT
         title = ""
      END SELECT
      CALL write_particle_coordinates(particle_set, iunit, outformat, "POS", title, qs_loc_env%cell, &
                                      unit_conv=unit_conv)
      CALL m_flush(iunit)
      CALL deallocate_particle_set(particle_set)
      CALL deallocate_atomic_kind_set(atomic_kind_set)
   END SUBROUTINE print_wannier_traj

! **************************************************************************************************
!> \brief Compute the second moments of the centers using the local (non-periodic) pos operators
!> \param qs_env ...
!> \param qs_loc_env ...
!> \param print_loc_section ...
!> \param ispin ...
!> \par History
!>      07.2020 created [H. Elgabarty]
!> \author Hossam Elgabarty
! **************************************************************************************************
   SUBROUTINE centers_second_moments_loc(qs_env, qs_loc_env, print_loc_section, ispin)

      ! I might not actually need the qs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_loc_env_type), POINTER                     :: qs_loc_env
      TYPE(section_vals_type), POINTER                   :: print_loc_section
      INTEGER, INTENT(IN)                                :: ispin

      INTEGER, PARAMETER                                 :: norder = 2
      LOGICAL, PARAMETER :: debug_this_subroutine = .FALSE.

      CHARACTER(len=default_path_length)                 :: file_tmp
      INTEGER                                            :: imoment, istate, ncol_global, nm, &
                                                            nmoloc, nrow_global, output_unit, &
                                                            output_unit_sm
      REAL(dp), DIMENSION(:, :), POINTER                 :: centers
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: moment_set
      REAL(KIND=dp), DIMENSION(3)                        :: rcc
      REAL(KIND=dp), DIMENSION(6)                        :: cov
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: momv, mvector, omvector
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: moloc_coeff
      TYPE(cp_fm_type), POINTER                          :: mo_localized
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, moments
      TYPE(mp_para_env_type), POINTER                    :: para_env

      logger => cp_get_default_logger()

      output_unit = cp_print_key_unit_nr(logger, print_loc_section, "PROGRAM_RUN_INFO", &
                                         extension=".locInfo")

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(/,T2,A)') &
            '!-----------------------------------------------------------------------------!'
         WRITE (output_unit, '(T12,A)') "Computing second moments of Wannier functions..."
         WRITE (output_unit, '(T2,A)') &
            '!-----------------------------------------------------------------------------!'
      END IF

      file_tmp = "_r_loc_cov_matrix_"//TRIM(ADJUSTL(cp_to_string(ispin)))
      output_unit_sm = cp_print_key_unit_nr(logger, print_loc_section, "WANNIER_SPREADS", &
                                            middle_name=file_tmp, extension=".data")

      CALL get_qs_loc_env(qs_loc_env=qs_loc_env, moloc_coeff=moloc_coeff)
      centers => qs_loc_env%localized_wfn_control%centers_set(ispin)%array

      CALL get_qs_env(qs_env, matrix_s=matrix_s)

      nm = ncoset(norder) - 1

      !CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, cell=cell)
      !particle_set => qs_loc_env%particle_set
      para_env => qs_loc_env%para_env

      CALL cp_fm_get_info(moloc_coeff(ispin), ncol_global=nmoloc)
      ALLOCATE (moment_set(nm, nmoloc))
      moment_set = 0.0_dp

      mo_localized => moloc_coeff(ispin)

      DO istate = 1, nmoloc
         rcc = centers(1:3, istate)
         !rcc = 0.0_dp

         ALLOCATE (moments(nm))
         DO imoment = 1, nm
            ALLOCATE (moments(imoment)%matrix)
            CALL dbcsr_copy(moments(imoment)%matrix, matrix_s(ispin)%matrix, 'MOM MAT')
            CALL dbcsr_set(moments(imoment)%matrix, 0.0_dp)
         END DO
         !

         CALL build_local_moment_matrix(qs_env, moments, norder, rcc)

         CALL cp_fm_get_info(mo_localized, ncol_global=ncol_global, nrow_global=nrow_global)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nrow_global, ncol_global=1, &
                                  para_env=mo_localized%matrix_struct%para_env, &
                                  context=mo_localized%matrix_struct%context)
         CALL cp_fm_create(mvector, fm_struct, name="mvector")
         CALL cp_fm_create(omvector, fm_struct, name="omvector")
         CALL cp_fm_create(momv, fm_struct, name="omvector")
         CALL cp_fm_struct_release(fm_struct)

         !cp_fm_to_fm_columns(msource, mtarget, ncol, source_start,target_start)
         CALL cp_fm_to_fm(mo_localized, mvector, 1, istate, 1)

         DO imoment = 1, nm
            CALL cp_dbcsr_sm_fm_multiply(moments(imoment)%matrix, mvector, omvector, 1)
            CALL cp_fm_schur_product(mvector, omvector, momv)
            !moment_set(imoment, istate) = moment_set(imoment, istate) + SUM(momv%local_data)
            moment_set(imoment, istate) = SUM(momv%local_data)
         END DO
         !
         CALL cp_fm_release(mvector)
         CALL cp_fm_release(omvector)
         CALL cp_fm_release(momv)

         DO imoment = 1, nm
            CALL dbcsr_deallocate_matrix(moments(imoment)%matrix)
         END DO
         DEALLOCATE (moments)
      END DO

      CALL para_env%sum(moment_set)

      IF (output_unit_sm > 0) THEN
         WRITE (UNIT=output_unit_sm, FMT='(A,T13,A,I1)') "#", &
            "SECOND MOMENTS OF WANNIER CENTERS FOR SPIN ", ispin
         WRITE (UNIT=output_unit_sm, FMT='(A,T29,A2,5(14x,A2))') "#", "XX", "XY", "XZ", "YY", "YZ", "ZZ"
         DO istate = 1, nmoloc
            cov = 0.0_dp
            cov(1) = moment_set(4, istate) - moment_set(1, istate)*moment_set(1, istate)
            cov(2) = moment_set(5, istate) - moment_set(1, istate)*moment_set(2, istate)
            cov(3) = moment_set(6, istate) - moment_set(1, istate)*moment_set(3, istate)
            cov(4) = moment_set(7, istate) - moment_set(2, istate)*moment_set(2, istate)
            cov(5) = moment_set(8, istate) - moment_set(2, istate)*moment_set(3, istate)
            cov(6) = moment_set(9, istate) - moment_set(3, istate)*moment_set(3, istate)
            WRITE (UNIT=output_unit_sm, FMT='(T4,A,I6,6(T20,6F16.8))') "Center:", istate, cov(1:6)
            IF (debug_this_subroutine) THEN
               WRITE (UNIT=output_unit_sm, FMT='(T4,A,I5,9(T20,9F12.6))') "Center:", istate, moment_set(1:, istate)
            END IF
         END DO
      END IF
      CALL cp_print_key_finished_output(output_unit_sm, logger, print_loc_section, "WANNIER_SPREADS")
      CALL cp_print_key_finished_output(output_unit, logger, print_loc_section, "PROGRAM_RUN_INFO")

      DEALLOCATE (moment_set)

   END SUBROUTINE centers_second_moments_loc

! **************************************************************************************************
!> \brief Compute the second moments of the centers using a periodic quadrupole operator
!> \brief c.f. Wheeler et al. PRB 100 245135 2019, and Kang et al. PRB 100 245134 2019
!> \brief The calculation is based on a Taylor expansion of the exp(i k_i r_i r_j k_j) operator to
!> \brief to first order in the cosine and the sine parts
!> \param qs_env ...
!> \param qs_loc_env ...
!> \param print_loc_section ...
!> \param ispin ...
!> \par History
!>      08.2020 created [H. Elgabarty]
!> \author Hossam Elgabarty
! **************************************************************************************************
   SUBROUTINE centers_second_moments_berry(qs_env, qs_loc_env, print_loc_section, ispin)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_loc_env_type), POINTER                     :: qs_loc_env
      TYPE(section_vals_type), POINTER                   :: print_loc_section
      INTEGER, INTENT(IN)                                :: ispin

      INTEGER, PARAMETER                                 :: taylor_order = 1
      LOGICAL, PARAMETER :: debug_this_subroutine = .FALSE.

      CHARACTER(len=default_path_length)                 :: file_tmp
      COMPLEX(KIND=dp)                                   :: z
      INTEGER                                            :: icov, imoment, istate, ncol_global, nm, &
                                                            nmoloc, norder, nrow_global, &
                                                            output_unit, output_unit_sm
      REAL(dp), DIMENSION(:, :), POINTER                 :: centers
      REAL(KIND=dp)                                      :: imagpart, Lx, Ly, Lz, realpart
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: moment_set
      REAL(KIND=dp), DIMENSION(3)                        :: rcc
      REAL(KIND=dp), DIMENSION(6)                        :: cov
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: momv, mvector, omvector
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: moloc_coeff
      TYPE(cp_fm_type), POINTER                          :: mo_localized
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, moments
      TYPE(mp_para_env_type), POINTER                    :: para_env

! two pointers, one to each spin channel's coeff. matrix (nao x nmoloc)

      logger => cp_get_default_logger()

      output_unit = cp_print_key_unit_nr(logger, print_loc_section, "PROGRAM_RUN_INFO", &
                                         extension=".locInfo")

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(/,T2,A)') &
            '!-----------------------------------------------------------------------------!'
         WRITE (output_unit, '(T12,A)') "Computing second moments of Wannier functions..."
         WRITE (output_unit, '(T2,A)') &
            '!-----------------------------------------------------------------------------!'
      END IF

      file_tmp = "_r_berry_cov_matrix_"//TRIM(ADJUSTL(cp_to_string(ispin)))
      output_unit_sm = cp_print_key_unit_nr(logger, print_loc_section, "WANNIER_SPREADS", &
                                            middle_name=file_tmp, extension=".data")

      norder = 2*(2*taylor_order + 1)
      nm = (6 + 11*norder + 6*norder**2 + norder**3)/6 - 1

      NULLIFY (cell)
      CALL get_qs_loc_env(qs_loc_env=qs_loc_env, moloc_coeff=moloc_coeff, cell=cell)
      Lx = cell%hmat(1, 1)
      Ly = cell%hmat(2, 2)
      Lz = cell%hmat(3, 3)

      centers => qs_loc_env%localized_wfn_control%centers_set(ispin)%array

      para_env => qs_loc_env%para_env

      CALL get_qs_env(qs_env, matrix_s=matrix_s)

      CALL cp_fm_get_info(moloc_coeff(ispin), ncol_global=nmoloc)

      ALLOCATE (moment_set(nm, nmoloc))
      moment_set = 0.0_dp

      mo_localized => moloc_coeff(ispin)

      DO istate = 1, nmoloc
         rcc = centers(1:3, istate)
         !rcc = 0.0_dp

         ALLOCATE (moments(nm))
         DO imoment = 1, nm
            ALLOCATE (moments(imoment)%matrix)
            CALL dbcsr_copy(moments(imoment)%matrix, matrix_s(ispin)%matrix, 'MOM MAT')
            CALL dbcsr_set(moments(imoment)%matrix, 0.0_dp)
         END DO
         !

         CALL build_local_moment_matrix(qs_env, moments, norder, rcc)

         CALL cp_fm_get_info(mo_localized, ncol_global=ncol_global, nrow_global=nrow_global)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nrow_global, ncol_global=1, &
                                  para_env=mo_localized%matrix_struct%para_env, &
                                  context=mo_localized%matrix_struct%context)
         CALL cp_fm_create(mvector, fm_struct, name="mvector")
         CALL cp_fm_create(omvector, fm_struct, name="omvector")
         CALL cp_fm_create(momv, fm_struct, name="omvector")
         CALL cp_fm_struct_release(fm_struct)

         ! cp_fm_to_fm_columns(msource, mtarget, ncol, source_start,target_start)
         CALL cp_fm_to_fm(mo_localized, mvector, 1, istate, 1)

         DO imoment = 1, nm
            CALL cp_dbcsr_sm_fm_multiply(moments(imoment)%matrix, mvector, omvector, 1)
            CALL cp_fm_schur_product(mvector, omvector, momv)
            moment_set(imoment, istate) = SUM(momv%local_data)
         END DO
         !
         CALL cp_fm_release(mvector)
         CALL cp_fm_release(omvector)
         CALL cp_fm_release(momv)

         DO imoment = 1, nm
            CALL dbcsr_deallocate_matrix(moments(imoment)%matrix)
         END DO
         DEALLOCATE (moments)
      END DO

      CALL para_env%sum(moment_set)

      IF (output_unit_sm > 0) THEN
         WRITE (UNIT=output_unit_sm, FMT='(A,T13,A,I1)') "#", &
            "SECOND MOMENTS OF WANNIER CENTERS FOR SPIN ", ispin
         WRITE (UNIT=output_unit_sm, FMT='(A,T29,A2,5(14x,A2))') "#", "XX", "XY", "XZ", "YY", "YZ", "ZZ"
         DO istate = 1, nmoloc
            cov = 0.0_dp
            DO icov = 1, 6
               realpart = 0.0_dp
               imagpart = 0.0_dp
               z = CMPLX(realpart, imagpart, dp)
               SELECT CASE (icov)

                 !! XX
               CASE (1)
                  realpart = 1.0 - 0.5*twopi*twopi/Lx/Lx/Lx/Lx &
                     *moment_set(20, istate)
                  imagpart = twopi/Lx/Lx*moment_set(4, istate) &
                             - twopi*twopi*twopi/Lx/Lx/Lx/Lx/Lx/Lx/6.0*moment_set(56, istate)

                  ! third order
                  !realpart = realpart + twopi*twopi*twopi*twopi/Lx/Lx/Lx/Lx/Lx/Lx/Lx/Lx/24.0 * moment_set(120, istate)
                  !imagpart = imagpart + twopi*twopi*twopi*twopi*twopi/Lx/Lx/Lx/Lx/Lx/Lx/Lx/Lx/Lx/Lx/120 &
                  !     * moment_set(220, istate)

                  z = CMPLX(realpart, imagpart, dp)
                  cov(1) = Lx*Lx/twopi*AIMAG(LOG(z)) - Lx*Lx/twopi/twopi*moment_set(1, istate)*moment_set(1, istate)

                 !! XY
               CASE (2)
                  realpart = 1.0 - 0.5*twopi*twopi/Lx/Ly/Lx/Ly &
                     *moment_set(23, istate)
                  imagpart = twopi/Lx/Ly*moment_set(5, istate) &
                             - twopi*twopi*twopi/Lx/Lx/Lx/Ly/Ly/Ly/6.0*moment_set(62, istate)

                  ! third order
                  !realpart = realpart + twopi*twopi*twopi*twopi/Lx/Lx/Lx/Lx/Ly/Ly/Ly/Ly/24.0 * moment_set(130, istate)
                  !imagpart = imagpart + twopi*twopi*twopi*twopi*twopi/Lx/Lx/Lx/Lx/Lx/Ly/Ly/Ly/Ly/Ly/120 &
                  !     * moment_set(235, istate)

                  z = CMPLX(realpart, imagpart, dp)
                  cov(2) = Lx*Ly/twopi*AIMAG(LOG(z)) - Lx*Ly/twopi/twopi*moment_set(1, istate)*moment_set(2, istate)

                 !! XZ
               CASE (3)
                  realpart = 1.0 - 0.5*twopi*twopi/Lx/Lz/Lx/Lz &
                     *moment_set(25, istate)
                  imagpart = twopi/Lx/Lz*moment_set(6, istate) &
                             - twopi*twopi*twopi/Lx/Lx/Lx/Lz/Lz/Lz/6.0*moment_set(65, istate)

                  ! third order
                  !realpart = realpart + twopi*twopi*twopi*twopi/Lx/Lx/Lx/Lx/Lz/Lz/Lz/Lz/24.0 * moment_set(134, istate)
                  !imagpart = imagpart + twopi*twopi*twopi*twopi*twopi/Lx/Lx/Lx/Lx/Lx/Ly/Ly/Ly/Ly/Ly/120 &
                  !     * moment_set(240, istate)

                  z = CMPLX(realpart, imagpart, dp)
                  cov(3) = Lx*Lz/twopi*AIMAG(LOG(z)) - Lx*Lz/twopi/twopi*moment_set(1, istate)*moment_set(3, istate)

                 !! YY
               CASE (4)
                  realpart = 1.0 - 0.5*twopi*twopi/Ly/Ly/Ly/Ly &
                     *moment_set(30, istate)
                  imagpart = twopi/Ly/Ly*moment_set(7, istate) &
                             - twopi*twopi*twopi/Ly/Ly/Ly/Ly/Ly/Ly/6.0*moment_set(77, istate)

                  ! third order
                  !realpart = realpart + twopi*twopi*twopi*twopi/Ly/Ly/Ly/Ly/Ly/Ly/Ly/Ly/24.0 * moment_set(156, istate)
                  !imagpart = imagpart + twopi*twopi*twopi*twopi*twopi/Ly/Ly/Ly/Ly/Ly/Ly/Ly/Ly/Ly/Ly/120 &
                  !     * moment_set(275, istate)

                  z = CMPLX(realpart, imagpart, dp)
                  cov(4) = Ly*Ly/twopi*AIMAG(LOG(z)) - Ly*Ly/twopi/twopi*moment_set(2, istate)*moment_set(2, istate)

                 !! YZ
               CASE (5)
                  realpart = 1.0 - 0.5*twopi*twopi/Ly/Lz/Ly/Lz &
                     *moment_set(32, istate)
                  imagpart = twopi/Ly/Lz*moment_set(8, istate) &
                             - twopi*twopi*twopi/Ly/Ly/Ly/Lz/Lz/Lz/6.0*moment_set(80, istate)

                  ! third order
                  !realpart = realpart + twopi*twopi*twopi*twopi/Lz/Lz/Lz/Lz/Ly/Ly/Ly/Ly/24.0 * moment_set(160, istate)
                  !imagpart = imagpart + twopi*twopi*twopi*twopi*twopi/Lz/Lz/Lz/Lz/Lz/Ly/Ly/Ly/Ly/Ly/120 &
                  !     * moment_set(280, istate)

                  z = CMPLX(realpart, imagpart, dp)
                  cov(5) = Ly*Lz/twopi*AIMAG(LOG(z)) - Ly*Lz/twopi/twopi*moment_set(2, istate)*moment_set(3, istate)

                 !! ZZ
               CASE (6)
                  realpart = 1.0 - 0.5*twopi*twopi/Lz/Lz/Lz/Lz &
                     *moment_set(34, istate)
                  imagpart = twopi/Lz/Lz*moment_set(9, istate) &
                             - twopi*twopi*twopi/Lz/Lz/Lz/Lz/Lz/Lz/6.0*moment_set(83, istate)

                  ! third order
                  !realpart = realpart + twopi*twopi*twopi*twopi/Lz/Lz/Lz/Lz/Lz/Lz/Lz/Lz/24.0 * moment_set(164, istate)
                  !imagpart = imagpart + twopi*twopi*twopi*twopi*twopi/Lz/Lz/Lz/Lz/Lz/Lz/Lz/Lz/Lz/Lz/120 &
                  !     * moment_set(285, istate)

                  z = CMPLX(realpart, imagpart, dp)
                  cov(6) = Lz*Lz/twopi*AIMAG(LOG(z)) - Lz*Lz/twopi/twopi*moment_set(3, istate)*moment_set(3, istate)

               END SELECT
            END DO
            WRITE (UNIT=output_unit_sm, FMT='(T4,A,I6,6(T20,6F16.8))') "Center:", istate, cov(1:6)
            IF (debug_this_subroutine) THEN
               WRITE (UNIT=output_unit_sm, FMT='(T4,A,I5,9(T20,9F12.6))') "Center:", istate, moment_set(1:, istate)
            END IF
         END DO
      END IF
      CALL cp_print_key_finished_output(output_unit_sm, logger, print_loc_section, "WANNIER_SPREADS")
      CALL cp_print_key_finished_output(output_unit, logger, print_loc_section, "PROGRAM_RUN_INFO")

      DEALLOCATE (moment_set)

   END SUBROUTINE centers_second_moments_berry

END MODULE qs_loc_methods
