!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Post-SCF finite-system spectral localizers and analytic AO position export.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE qs_spectral_localizer
   USE bibliography,                    ONLY: Dixon2023,&
                                              Doll2025,&
                                              Loring2015,&
                                              Wong2026
   USE cell_types,                      ONLY: cell_type,&
                                              use_perd_none,&
                                              use_perd_xy,&
                                              use_perd_xyz,&
                                              use_perd_xz,&
                                              use_perd_yz
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: do_method_gapw,&
                                              do_method_gpw
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_finite_ao,                    ONLY: finite_ao_positions,&
                                              finite_ao_soc
   USE qs_localizer_torus,              ONLY: localizer_torus_operators
   USE qs_moments,                      ONLY: build_berry_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE reference_manager,               ONLY: cite_reference
   USE spectral_localizer,              ONLY: localizer_3d_matrix_evaluate,&
                                              localizer_matrix_evaluate
   USE spectral_localizer_bloch,        ONLY: localizer_flatten_bloch
   USE spectral_localizer_dbcsr,        ONLY: localizer_ao_operators,&
                                              localizer_dbcsr_3d_build,&
                                              localizer_dbcsr_build,&
                                              localizer_dbcsr_dense,&
                                              localizer_dbcsr_release,&
                                              localizer_matrix_type
   USE spectral_localizer_flatten,      ONLY: localizer_flatten_checked
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_spectral_localizer'
   PUBLIC :: qs_scf_post_localizer

CONTAINS

! **************************************************************************************************
!> \brief Export finite AO positions and optionally evaluate localizers after SCF.
!> \param qs_env Converged Quickstep environment
! **************************************************************************************************
   SUBROUTINE qs_scf_post_localizer(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_localizer'
      CHARACTER(len=1), PARAMETER                        :: axes(3) = ['X', 'Y', 'Z']
      INTEGER, PARAMETER :: plane_periodicity(3) = [use_perd_xy, use_perd_xz, use_perd_yz]

      COMPLEX(KIND=dp), ALLOCATABLE                      :: local_dense(:, :), metric_dense(:, :)
      INTEGER :: after, axis, directions(3), flatten_handle, formulation, handle, ie, ik, &
         invariant, ip, iw, max_ao, n, nao, ndim, nposition, output_unit, plane, poisson_periodic, &
         spin_channel, status, topology
      INTEGER, POINTER                                   :: grid(:)
      LOGICAL :: active, aii, do_kpoints, do_soc, explicit_scale, flatten, omit_headers, &
         periodic_localizer, print_position, print_soc, time_reversal, torus
      REAL(dp) :: build_energy, electronic_gap(2), eps_electronic_gap, eps_flattening, &
         flatten_residuals(4), flatten_scale, local_point(3), physical_ratio
      REAL(KIND=dp)                                      :: eps_gap, eps_metric, gap, &
                                                            matrix_tolerance, metric_ratio
      REAL(KIND=dp), POINTER                             :: energies(:), kappas(:), point(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), POINTER                        :: cosines(:), image_matrices(:, :), &
                                                            matrix_ks(:), matrix_s(:), &
                                                            positions(:), soc_matrices(:, :)
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(localizer_matrix_type)                        :: localizer, metric, physical_h(1), &
                                                            physical_metric
      TYPE(localizer_matrix_type), POINTER               :: flat_hamiltonian
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), POINTER            :: sab_orb(:)
      TYPE(section_vals_type), POINTER                   :: input, section

      CALL timeset(routineN, handle)
      NULLIFY (input, section, cell, dft_control, para_env, sab_orb, matrix_ks, matrix_s, &
               positions, soc_matrices, energies, kappas, point, cosines, flat_hamiltonian)
      CALL get_qs_env(qs_env, input=input)
      section => section_vals_get_subs_vals(input, 'PROPERTIES%SPECTRAL_LOCALIZER')
      CALL section_vals_val_get(section, '_SECTION_PARAMETERS_', l_val=active)
      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)
      print_position = BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                                        'DFT%PRINT%AO_MATRICES/POSITION'), cp_p_file)
      print_soc = BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                                   'DFT%PRINT%AO_MATRICES/SOC'), cp_p_file)
      IF (.NOT. active .AND. .NOT. print_position .AND. .NOT. print_soc) THEN
         CALL timestop(handle)
         RETURN
      END IF
      CALL get_qs_env(qs_env, cell=cell, do_kpoints=do_kpoints, dft_control=dft_control, &
                      matrix_s_kp=image_matrices, para_env=para_env, sab_orb=sab_orb)
      matrix_s => image_matrices(:, 1)
      CALL section_vals_val_get(section, 'FORMULATION', i_val=formulation)
      periodic_localizer = active .AND. formulation == 1
      IF ((print_position .OR. print_soc) .AND. (ANY(cell%perd /= 0) .OR. do_kpoints)) THEN
         CPABORT('AO_MATRICES POSITION/SOC requires finite systems without KPOINTS.')
      END IF
      IF (do_kpoints .AND. .NOT. periodic_localizer) THEN
         CPABORT('Explicit KPOINTS require FORMULATION PERIODIC for localizers.')
      END IF
      IF (.NOT. periodic_localizer .AND. ANY(cell%perd /= 0)) THEN
         CPABORT('AO position/localizer requires PERIODIC NONE without KPOINTS.')
      END IF
      CALL section_vals_val_get(section, 'DIMENSION', i_val=ndim)
      IF (ndim /= 2 .AND. ndim /= 3) CPABORT('SPECTRAL_LOCALIZER DIMENSION must be 2 or 3.')
      CALL section_vals_val_get(section, 'PLANE', i_val=plane, explicit=explicit_scale)
      IF (active .AND. ndim == 3 .AND. explicit_scale) THEN
         CPABORT('PLANE is only valid for DIMENSION 2.')
      END IF
      CALL section_vals_val_get(section, 'MP_GRID', i_vals=grid, explicit=explicit_scale)
      IF (ANY(grid < 1)) CPABORT('SPECTRAL_LOCALIZER MP_GRID must be positive.')
      IF (.NOT. periodic_localizer .AND. explicit_scale) THEN
         CPABORT('SPECTRAL_LOCALIZER MP_GRID requires FORMULATION PERIODIC.')
      END IF
      IF (ANY(grid /= 1 .AND. cell%perd == 0)) THEN
         CPABORT('Localizer MP_GRID must be one in nonperiodic directions.')
      END IF
      torus = periodic_localizer .AND. (do_kpoints .OR. ANY(grid /= 1))
      directions = [1, 2, 3]
      IF (ndim == 2) THEN
         SELECT CASE (plane)
         CASE (1)
            directions(:2) = [1, 2]
         CASE (2)
            directions(:2) = [1, 3]
         CASE (3)
            directions(:2) = [2, 3]
         END SELECT
      END IF
      CALL section_vals_val_get(input, 'DFT%POISSON%PERIODIC', i_val=poisson_periodic)
      IF (periodic_localizer) THEN
         IF (ndim == 3) THEN
            IF (ANY(cell%perd == 0) .OR. poisson_periodic /= use_perd_xyz) THEN
               CPABORT('PERIODIC DIMENSION 3 requires CELL and POISSON PERIODIC XYZ.')
            END IF
         ELSE
            IF (COUNT(cell%perd /= 0) /= 2 .OR. ANY(cell%perd(directions(:2)) == 0)) THEN
               CPABORT('PERIODIC localizer requires two periodic CELL axes matching PLANE.')
            END IF
            IF (poisson_periodic /= plane_periodicity(plane)) THEN
               CPABORT('PERIODIC localizer requires matching CELL and POISSON periodicity.')
            END IF
         END IF
         CALL section_vals_val_get(section, 'KAPPA', explicit=explicit_scale)
         IF (explicit_scale) CPABORT('Use ETA, not KAPPA, for FORMULATION PERIODIC.')
      ELSE IF (active .AND. poisson_periodic /= use_perd_none) THEN
         CPABORT('SPECTRAL_LOCALIZER requires POISSON PERIODIC NONE.')
      END IF
      IF (active .AND. .NOT. periodic_localizer) THEN
         CALL section_vals_val_get(section, 'ETA', explicit=explicit_scale)
         IF (explicit_scale) CPABORT('Use KAPPA, not ETA, for FORMULATION FINITE.')
      END IF
      SELECT CASE (dft_control%qs_control%method_id)
      CASE (do_method_gpw, do_method_gapw)
      CASE DEFAULT
         CPABORT('AO position/localizer currently supports GPW and GAPW only.')
      END SELECT
      IF (active) THEN
         IF (.NOT. qs_env%scf_convergence_available .OR. .NOT. qs_env%scf_converged) THEN
            CPABORT('SPECTRAL_LOCALIZER requires a converged SCF calculation.')
         END IF
         CALL cite_reference(Dixon2023)
         IF (periodic_localizer) CALL cite_reference(Doll2025)
         CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)
         CALL section_vals_val_get(section, 'MAX_AO', i_val=max_ao)
         IF (max_ao < 1 .OR. nao > max_ao) THEN
            CPABORT('SPECTRAL_LOCALIZER exceeds the MAX_AO limit.')
         END IF
         DO axis = 1, 3
            IF (grid(axis) > max_ao/nao) THEN
               CPABORT('SPECTRAL_LOCALIZER torus exceeds the MAX_AO limit.')
            END IF
            nao = nao*grid(axis)
         END DO
      END IF

      IF (periodic_localizer .AND. .NOT. torus) THEN
         CALL periodic_ao_moments(qs_env, directions(:ndim), cosines, positions)
      ELSE IF (.NOT. periodic_localizer .AND. (active .OR. print_position)) THEN
         CALL finite_ao_positions(qs_env, positions)
      END IF
      IF (print_position) THEN
         CALL section_vals_val_get(input, 'DFT%PRINT%AO_MATRICES%NDIGITS', i_val=after)
         CALL section_vals_val_get(input, 'DFT%PRINT%AO_MATRICES%OMIT_HEADERS', l_val=omit_headers)
         after = MIN(MAX(after, 1), 16)
         DO axis = 1, 3
            iw = cp_print_key_unit_nr(logger, input, 'DFT%PRINT%AO_MATRICES/POSITION', &
                                      extension='.Log', middle_name=axes(axis))
            IF (iw > 0 .AND. .NOT. omit_headers) THEN
               WRITE (iw, '(A)') 'AO POSITION '//axes(axis)//' [bohr]; coordinate origin; covariant AO matrix'
            END IF
            CALL cp_dbcsr_write_sparse_matrix(positions(axis)%matrix, 4, after, qs_env, para_env, &
                                              output_unit=iw, omit_headers=omit_headers)
            CALL cp_print_key_finished_output(iw, logger, input, 'DFT%PRINT%AO_MATRICES/POSITION')
         END DO
      END IF
      IF (print_soc) THEN
         CALL finite_ao_soc(qs_env, soc_matrices)
         CALL section_vals_val_get(input, 'DFT%PRINT%AO_MATRICES%NDIGITS', i_val=after)
         CALL section_vals_val_get(input, 'DFT%PRINT%AO_MATRICES%OMIT_HEADERS', l_val=omit_headers)
         after = MIN(MAX(after, 1), 16)
         DO axis = 1, 3
            iw = cp_print_key_unit_nr(logger, input, 'DFT%PRINT%AO_MATRICES/SOC', &
                                      extension='.Log', middle_name='SOC_'//axes(axis))
            IF (iw > 0 .AND. .NOT. omit_headers) THEN
               WRITE (iw, '(A)') 'AO SOC '//axes(axis)//' [hartree]; real antisymmetric component; multiply by i'
            END IF
            CALL cp_dbcsr_write_sparse_matrix(soc_matrices(axis, 1)%matrix, 4, after, qs_env, para_env, &
                                              output_unit=iw, omit_headers=omit_headers)
            CALL cp_print_key_finished_output(iw, logger, input, 'DFT%PRINT%AO_MATRICES/SOC')
         END DO
      END IF
      IF (.NOT. active) THEN
         IF (ASSOCIATED(positions)) CALL dbcsr_deallocate_matrix_set(positions)
         IF (ASSOCIATED(soc_matrices)) CALL dbcsr_deallocate_matrix_set(soc_matrices)
         CALL timestop(handle)
         RETURN
      END IF

      CALL section_vals_val_get(section, 'INVARIANT', i_val=topology)
      aii = topology == 1
      IF (ndim == 3 .AND. .NOT. aii) CPABORT('DIMENSION 3 requires INVARIANT Z2.')
      IF (aii .AND. ndim == 2) CALL cite_reference(Wong2026)
      IF (ndim == 3) CALL cite_reference(Loring2015)
      CALL section_vals_val_get(section, 'SOC', l_val=do_soc)
      CALL section_vals_val_get(section, 'TIME_REVERSAL', l_val=time_reversal)
      CALL section_vals_val_get(section, 'SPIN_CHANNEL', i_val=spin_channel)
      IF (aii .AND. (.NOT. do_soc .OR. .NOT. time_reversal)) THEN
         CPABORT('SPECTRAL_LOCALIZER INVARIANT Z2 requires SOC T and TIME_REVERSAL T.')
      END IF
      CALL get_qs_env(qs_env, matrix_ks_kp=image_matrices)
      matrix_ks => image_matrices(:, 1)
      IF (spin_channel < 1 .OR. spin_channel > SIZE(matrix_ks)) THEN
         CPABORT('SPECTRAL_LOCALIZER SPIN_CHANNEL is unavailable.')
      END IF
      IF (do_soc .AND. (SIZE(matrix_ks) /= 1 .OR. spin_channel /= 1)) THEN
         CPABORT('SPECTRAL_LOCALIZER SOC requires a restricted SCF and SPIN_CHANNEL 1.')
      END IF
      CALL section_vals_val_get(section, 'ENERGY', r_vals=energies)
      IF (periodic_localizer) THEN
         CALL section_vals_val_get(section, 'ETA', r_vals=kappas)
      ELSE
         CALL section_vals_val_get(section, 'KAPPA', r_vals=kappas)
      END IF
      IF (.NOT. ALL(ieee_is_finite(energies)) .OR. .NOT. ALL(ieee_is_finite(kappas))) THEN
         CPABORT('SPECTRAL_LOCALIZER requires finite ENERGY and KAPPA/ETA.')
      END IF
      IF (ANY(kappas <= 0.0_dp)) THEN
         CPABORT('SPECTRAL_LOCALIZER requires positive KAPPA/ETA.')
      END IF
      CALL section_vals_val_get(section, 'EPS_GAP', r_val=eps_gap)
      CALL section_vals_val_get(section, 'EPS_METRIC', r_val=eps_metric)
      CALL section_vals_val_get(section, 'MATRIX_TOLERANCE', r_val=matrix_tolerance)
      IF (.NOT. ALL(ieee_is_finite([eps_gap, eps_metric, matrix_tolerance]))) THEN
         CPABORT('SPECTRAL_LOCALIZER requires finite numerical tolerances.')
      END IF
      IF (eps_gap < 0.0_dp .OR. eps_metric <= 0.0_dp .OR. eps_metric >= 1.0_dp .OR. &
          matrix_tolerance <= 0.0_dp .OR. matrix_tolerance >= 1.0_dp) THEN
         CPABORT('SPECTRAL_LOCALIZER: invalid EPS_GAP, EPS_METRIC or MATRIX_TOLERANCE.')
      END IF
      CALL section_vals_val_get(section, 'SPECTRAL_FLATTENING', l_val=flatten)
      CALL section_vals_val_get(section, 'FLATTENING_SCALE', r_val=flatten_scale)
      CALL section_vals_val_get(section, 'EPS_FLATTENING', r_val=eps_flattening)
      CALL section_vals_val_get(section, 'EPS_ELECTRONIC_GAP', r_val=eps_electronic_gap)
      IF (flatten) THEN
         IF (.NOT. ALL(ieee_is_finite([flatten_scale, eps_flattening, eps_electronic_gap]))) THEN
            CPABORT('SPECTRAL_FLATTENING requires finite scales and tolerances.')
         END IF
         IF (flatten_scale <= 0.0_dp .OR. eps_flattening <= 0.0_dp .OR. eps_flattening >= 1.0_dp .OR. &
             eps_electronic_gap < 0.0_dp) THEN
            CPABORT('SPECTRAL_FLATTENING: invalid scale or tolerance.')
         END IF
      END IF
      CALL section_vals_val_get(section, 'POSITION', n_rep_val=nposition)
      n = nao
      IF (do_soc) n = 2*nao
      IF (torus) THEN
         NULLIFY (matrix_ks, matrix_s)
         CALL localizer_torus_operators(qs_env, grid, directions(:ndim), do_soc, matrix_ks, matrix_s, &
                                        cosines, positions, soc_matrices)
      ELSE IF (do_soc) THEN
         IF (.NOT. ASSOCIATED(soc_matrices)) CALL finite_ao_soc(qs_env, soc_matrices)
      END IF
      IF (output_unit > 0) THEN
         IF (periodic_localizer) THEN
            WRITE (output_unit, '(/,T2,A)') 'SPECTRAL_LOCALIZER| Periodic Born-von-Karman DBCSR localizer'
            WRITE (output_unit, '(T2,A,3I6)') 'SPECTRAL_LOCALIZER| Full analysis mesh: ', grid
            WRITE (output_unit, '(T2,A)') 'SPECTRAL_LOCALIZER| Analytic cos(G.r)/sin(G.r); no finite-coordinate substitution.'
            WRITE (output_unit, '(T2,A)') 'SPECTRAL_LOCALIZER| Finite AO torus diagnostic; converge basis, supercell and ETA.'
         ELSE
            WRITE (output_unit, '(/,T2,A)') 'SPECTRAL_LOCALIZER| Finite-system DBCSR localizer'
         END IF
         IF (ndim == 3) THEN
            WRITE (output_unit, '(T2,A)') 'SPECTRAL_LOCALIZER| LAPACK real LU/QR sign; independent generalized chiral gap.'
         ELSE
            WRITE (output_unit, '(T2,A)') 'SPECTRAL_LOCALIZER| LAPACK LDL inertia with eigenvalue validation.'
         END IF
         WRITE (output_unit, '(T2,A,I0)') 'SPECTRAL_LOCALIZER| Localizer dimension: ', MERGE(4*n, 2*n, ndim == 3)
         IF (ndim == 3) THEN
            WRITE (output_unit, '(T2,A)') 'SPECTRAL_LOCALIZER| Three-dimensional AII; chiral determinant sign.'
         ELSE
            WRITE (output_unit, '(T2,A,2A)') 'SPECTRAL_LOCALIZER| Oriented plane: ', axes(directions(:2))
         END IF
         IF (.NOT. periodic_localizer) THEN
            WRITE (output_unit, '(T2,A)') 'SPECTRAL_LOCALIZER| Check basis, volume and KAPPA convergence.'
         END IF
         IF (.NOT. do_soc) WRITE (output_unit, '(T2,A)') &
            'SPECTRAL_LOCALIZER| Real scalar Hamiltonian: no nonzero class A Chern index is expected.'
      END IF
      IF (flatten) THEN
         ALLOCATE (flat_hamiltonian)
         IF (do_soc) THEN
            CALL localizer_ao_operators(matrix_ks(spin_channel)%matrix, matrix_s(1)%matrix, positions(:0), &
                                        physical_h, physical_metric, soc_matrices)
         ELSE
            CALL localizer_ao_operators(matrix_ks(spin_channel)%matrix, matrix_s(1)%matrix, positions(:0), &
                                        physical_h, physical_metric)
         END IF
         IF (output_unit > 0) THEN
            WRITE (output_unit, '(T2,A,ES23.14)') 'SPECTRAL_LOCALIZER| Flattening scale [hartree]: ', flatten_scale
            WRITE (output_unit, '(T2,A)') 'SPECTRAL_LOCALIZER| Gaps and indices refer to the flattened Hamiltonian.'
            IF (periodic_localizer) THEN
               WRITE (output_unit, '(T2,A)') 'SPECTRAL_LOCALIZER| Flattening in complete primitive Bloch blocks.'
               WRITE (output_unit, '(T2,A)') &
                  'SPECTRAL_LOCALIZER| Residuals: translation, metric orthonormality, eigenproblem, Hermiticity.'
            END IF
         END IF
      END IF
      DO ie = 1, SIZE(energies)
         build_energy = energies(ie)
         IF (flatten) THEN
            CALL timeset('localizer_spectral_flattening', flatten_handle)
            IF (periodic_localizer) THEN
               CALL localizer_flatten_bloch(physical_h(1), physical_metric, grid, MERGE(2, 1, do_soc), &
                                            energies(ie), flatten_scale, eps_flattening, eps_metric, eps_electronic_gap, &
                                            para_env, flat_hamiltonian, electronic_gap(1), physical_ratio, &
                                            flatten_residuals, status)
               electronic_gap(2) = electronic_gap(1)
            ELSE
               CALL localizer_flatten_checked(physical_h(1), physical_metric, energies(ie), flatten_scale, eps_flattening, &
                                              eps_metric, eps_electronic_gap, para_env, flat_hamiltonian, &
                                              electronic_gap, physical_ratio, flatten_residuals, status)
            END IF
            CALL timestop(flatten_handle)
            IF (output_unit > 0) THEN
               WRITE (output_unit, '(T2,A,ES23.14)') 'SPECTRAL_LOCALIZER| Flattening query [hartree]: ', energies(ie)
               WRITE (output_unit, '(T2,A,2ES23.14)') &
                  'SPECTRAL_LOCALIZER| Original electronic gap bracket [hartree]: ', electronic_gap
               WRITE (output_unit, '(T2,A,ES23.14)') 'SPECTRAL_LOCALIZER| Original metric ratio bound: ', physical_ratio
               WRITE (output_unit, '(T2,A,4ES14.5)') 'SPECTRAL_LOCALIZER| Flattening residuals: ', flatten_residuals
            END IF
            SELECT CASE (status)
            CASE (0)
            CASE (1)
               CPABORT('SPECTRAL_FLATTENING: electronic gap unresolved; no index.')
            CASE (-2)
               CPABORT('SPECTRAL_FLATTENING: invalid original AO metric.')
            CASE (-6)
               CPABORT('SPECTRAL_FLATTENING: sign iteration or structure failed.')
            CASE (-7)
               CPABORT('SPECTRAL_FLATTENING: broken torus translations.')
            CASE DEFAULT
               CPABORT('SPECTRAL_FLATTENING: input or gap factorization failed.')
            END SELECT
            build_energy = 0.0_dp
         END IF
         DO ip = 1, nposition
            CALL section_vals_val_get(section, 'POSITION', i_rep_val=ip, r_vals=point)
            IF (.NOT. ALL(ieee_is_finite(point))) THEN
               CPABORT('SPECTRAL_LOCALIZER requires finite POSITION.')
            END IF
            local_point(:ndim) = point(directions(:ndim))
            IF (periodic_localizer) THEN
               DO axis = 1, ndim
                  local_point(axis) = twopi*DOT_PRODUCT(cell%h_inv(directions(axis), :), point)/ &
                                      REAL(grid(directions(axis)), dp)
               END DO
            END IF
            DO ik = 1, SIZE(kappas)
               IF (ndim == 3) THEN
                  CALL localizer_dbcsr_3d_build(matrix_ks(spin_channel)%matrix, matrix_s(1)%matrix, positions, &
                                                local_point, build_energy, kappas(ik), localizer, metric, soc_matrices, &
                                                shifted_hamiltonian=flat_hamiltonian, periodic_cosines=cosines)
               ELSE IF (do_soc) THEN
                  CALL localizer_dbcsr_build(matrix_ks(spin_channel)%matrix, matrix_s(1)%matrix, &
                                             positions(directions(1))%matrix, positions(directions(2))%matrix, &
                                             local_point(:2), build_energy, kappas(ik), localizer, metric, soc_matrices, cosines, &
                                             shifted_hamiltonian=flat_hamiltonian)
               ELSE
                  CALL localizer_dbcsr_build(matrix_ks(spin_channel)%matrix, matrix_s(1)%matrix, &
                                             positions(directions(1))%matrix, positions(directions(2))%matrix, &
                                             local_point(:2), build_energy, kappas(ik), localizer, metric, &
                                             periodic_cosines=cosines, shifted_hamiltonian=flat_hamiltonian)
               END IF
               CALL localizer_dbcsr_dense(localizer, para_env, local_dense)
               CALL localizer_dbcsr_dense(metric, para_env, metric_dense)
               IF (para_env%is_source()) THEN
                  IF (ndim == 3) THEN
                     CALL localizer_3d_matrix_evaluate(local_dense, metric_dense(:n, :n), &
                                                       matrix_tolerance, eps_metric, eps_gap, invariant, gap, &
                                                       metric_ratio, status)
                  ELSE
                     CALL localizer_matrix_evaluate(local_dense, metric_dense(:n, :n), aii, &
                                                    matrix_tolerance, eps_metric, eps_gap, invariant, gap, metric_ratio, status)
                  END IF
                  IF (status < 0 .AND. output_unit > 0 .AND. ndim == 2) THEN
                     WRITE (output_unit, '(T2,A,I0,2ES14.5)') &
                        'SPECTRAL_LOCALIZER| Failure status and L/B Hermiticity residuals: ', status, &
                        MAXVAL(ABS(local_dense - CONJG(TRANSPOSE(local_dense)))), &
                        MAXVAL(ABS(metric_dense - CONJG(TRANSPOSE(metric_dense))))
                  END IF
               END IF
               DEALLOCATE (local_dense, metric_dense)
               CALL localizer_dbcsr_release(localizer)
               CALL localizer_dbcsr_release(metric)
               CALL para_env%bcast(status, para_env%source)
               IF (status < 0 .OR. status > 1) THEN
                  SELECT CASE (status)
                  CASE (-2)
                     CPABORT('SPECTRAL_LOCALIZER: invalid AO metric; check basis/EPS_METRIC.')
                  CASE (-3)
                     CPABORT('SPECTRAL_LOCALIZER: broken physical TR/AII structure.')
                  CASE DEFAULT
                     CPABORT('SPECTRAL_LOCALIZER: invalid data or factorization failure.')
                  END SELECT
               END IF
               CALL para_env%bcast(invariant, para_env%source)
               CALL para_env%bcast(gap, para_env%source)
               CALL para_env%bcast(metric_ratio, para_env%source)
               IF (output_unit > 0) THEN
                  WRITE (output_unit, '(T2,A,3ES23.14)') 'SPECTRAL_LOCALIZER| Position [bohr]: ', point
                  WRITE (output_unit, '(T2,A,ES23.14)') 'SPECTRAL_LOCALIZER| Energy [hartree]: ', energies(ie)
                  IF (periodic_localizer) THEN
                     WRITE (output_unit, '(T2,A,ES23.14)') 'SPECTRAL_LOCALIZER| Eta [hartree]: ', kappas(ik)
                  ELSE
                     WRITE (output_unit, '(T2,A,ES23.14)') 'SPECTRAL_LOCALIZER| Kappa [hartree/bohr]: ', kappas(ik)
                  END IF
                  WRITE (output_unit, '(T2,A,ES23.14)') 'SPECTRAL_LOCALIZER| Metric ratio: ', metric_ratio
                  WRITE (output_unit, '(T2,A,ES23.14)') 'SPECTRAL_LOCALIZER| Gap [hartree]: ', gap
                  IF (status == 1) THEN
                     WRITE (output_unit, '(T2,A)') 'SPECTRAL_LOCALIZER| UNRESOLVED: gap or factorization below accuracy; no index.'
                  ELSE IF (aii) THEN
                     WRITE (output_unit, '(T2,A,I0)') 'SPECTRAL_LOCALIZER| Z2 index: ', invariant
                  ELSE
                     WRITE (output_unit, '(T2,A,I0)') 'SPECTRAL_LOCALIZER| Chern index: ', invariant
                  END IF
               END IF
            END DO
         END DO
         IF (flatten) CALL localizer_dbcsr_release(flat_hamiltonian)
      END DO
      IF (flatten) THEN
         CALL localizer_dbcsr_release(physical_h(1))
         CALL localizer_dbcsr_release(physical_metric)
         DEALLOCATE (flat_hamiltonian)
      END IF
      CALL dbcsr_deallocate_matrix_set(positions)
      IF (ASSOCIATED(cosines)) CALL dbcsr_deallocate_matrix_set(cosines)
      IF (ASSOCIATED(soc_matrices)) CALL dbcsr_deallocate_matrix_set(soc_matrices)
      IF (torus) THEN
         CALL dbcsr_deallocate_matrix_set(matrix_ks)
         CALL dbcsr_deallocate_matrix_set(matrix_s)
      END IF
      CALL timestop(handle)
   END SUBROUTINE qs_scf_post_localizer

! **************************************************************************************************
!> \brief Analytic torus coordinates in the reciprocal frame of a Gamma supercell.
!> \param qs_env Quickstep environment
!> \param directions Two or three periodic lattice directions
!> \param cosines Newly allocated cosine matrices, indexed by oriented direction
!> \param sines Newly allocated sine matrices, indexed by lattice direction
! **************************************************************************************************
   SUBROUTINE periodic_ao_moments(qs_env, directions, cosines, sines)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: directions(:)
      TYPE(dbcsr_p_type), POINTER                        :: cosines(:), sines(:)

      INTEGER                                            :: j
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), POINTER                        :: matrix_s(:)
      TYPE(neighbor_list_set_p_type), POINTER            :: sab_orb(:)

      CALL get_qs_env(qs_env, cell=cell, matrix_s=matrix_s, sab_orb=sab_orb)
      ALLOCATE (cosines(SIZE(directions)), sines(3))
      DO j = 1, 3
         ALLOCATE (sines(j)%matrix)
         CALL dbcsr_create(sines(j)%matrix, template=matrix_s(1)%matrix, name='TORUS SIN')
         CALL cp_dbcsr_alloc_block_from_nbl(sines(j)%matrix, sab_orb)
         CALL dbcsr_set(sines(j)%matrix, 0.0_dp)
      END DO
      DO j = 1, SIZE(directions)
         ALLOCATE (cosines(j)%matrix)
         CALL dbcsr_create(cosines(j)%matrix, template=matrix_s(1)%matrix, name='TORUS COS')
         CALL cp_dbcsr_alloc_block_from_nbl(cosines(j)%matrix, sab_orb)
         CALL dbcsr_set(cosines(j)%matrix, 0.0_dp)
         CALL build_berry_moment_matrix(qs_env, cosines(j)%matrix, sines(directions(j))%matrix, &
                                        twopi*cell%h_inv(directions(j), :))
      END DO
   END SUBROUTINE periodic_ao_moments

END MODULE qs_spectral_localizer
