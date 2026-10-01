!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Finite or periodic quadratic pseudospectrum analysis after SCF.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE qs_quadratic_pseudospectrum
   USE bibliography,                    ONLY: Bairnsfather2026,&
                                              Cerjan2023
   USE cell_types,                      ONLY: cell_transform_input_cartesian,&
                                              cell_transform_input_reciprocal,&
                                              cell_type,&
                                              use_perd_none
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_info,&
                                              dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
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
   USE qs_ao_translations,              ONLY: ao_translation_matrices
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_finite_ao,                    ONLY: finite_ao_positions,&
                                              finite_ao_soc
   USE qs_localizer_torus,              ONLY: localizer_torus_operators,&
                                              torus_translation_operators
   USE quadratic_pseudospectrum,        ONLY: quadratic_dense,&
                                              quadratic_result_type
   USE reference_manager,               ONLY: cite_reference
   USE spectral_localizer_dbcsr,        ONLY: localizer_ao_operators,&
                                              localizer_dbcsr_check,&
                                              localizer_dbcsr_dense,&
                                              localizer_dbcsr_release,&
                                              localizer_matrix_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_quadratic_pseudospectrum'
   PUBLIC :: qs_scf_post_quadratic
CONTAINS
! **************************************************************************************************
!> \brief Evaluate scans of quadratic gaps and localized joint approximate eigenstates.
!> \param qs_env Converged Quickstep environment
! **************************************************************************************************
   SUBROUTINE qs_scf_post_quadratic(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_quadratic'

      CHARACTER(len=256)                                 :: matrix_error
      COMPLEX(dp), ALLOCATABLE                           :: complex_shift(:), dense_metric(:, :), &
                                                            dense_operators(:, :, :), &
                                                            temporary(:, :)
      INTEGER :: ao, axis, cell_periodic, formulation, handle, i, ie, ik, ip, iw, j, max_ao, &
         max_iter, max_space, n, nao, ndim, nobs, npositions, nstates, nt, output_unit, periodic, &
         spin, status
      INTEGER, ALLOCATABLE                               :: directions(:)
      INTEGER, POINTER                                   :: grid(:)
      LOGICAL :: active, do_kpoints, do_soc, explicit_grid, explicit_value, periodic_analysis, &
         torus_analysis, translation_analysis, valid
      LOGICAL, ALLOCATABLE                               :: unitary(:)
      REAL(dp)                                           :: eps_eigen, eps_metric, &
                                                            hermiticity_error, &
                                                            max_hermiticity_error, query(3)
      REAL(dp), ALLOCATABLE                              :: angles(:), radii(:), shift(:), &
                                                            translations(:, :), weights(:)
      REAL(dp), POINTER                                  :: energies(:), kappas(:), point(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), POINTER                        :: cosines(:), images_h(:, :), &
                                                            images_s(:, :), matrix_ks(:), &
                                                            matrix_s(:), positions(:), sines(:), &
                                                            soc(:, :)
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(localizer_matrix_type)                        :: metric
      TYPE(localizer_matrix_type), ALLOCATABLE           :: operators(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(quadratic_result_type)                        :: result
      TYPE(section_vals_type), POINTER                   :: input, section

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, input=input)
      section => section_vals_get_subs_vals(input, 'PROPERTIES%QUADRATIC_PSEUDOSPECTRUM')
      CALL section_vals_val_get(section, '_SECTION_PARAMETERS_', l_val=active)
      IF (.NOT. active) THEN
         CALL timestop(handle)
         RETURN
      END IF
      IF (.NOT. qs_env%scf_convergence_available .OR. .NOT. qs_env%scf_converged) THEN
         CPABORT('QUADRATIC_PSEUDOSPECTRUM requires converged SCF.')
      END IF
      CALL get_qs_env(qs_env, matrix_s_kp=images_s, matrix_ks_kp=images_h, para_env=para_env, &
                      cell=cell, dft_control=dft_control, do_kpoints=do_kpoints)
      matrix_s => images_s(:, 1)
      matrix_ks => images_h(:, 1)
      SELECT CASE (dft_control%qs_control%method_id)
      CASE (do_method_gpw, do_method_gapw)
      CASE DEFAULT
         CPABORT('QUADRATIC_PSEUDOSPECTRUM supports GPW/GAPW only.')
      END SELECT
      CALL section_vals_val_get(section, 'FORMULATION', i_val=formulation)
      periodic_analysis = formulation == 1
      translation_analysis = formulation == 2
      torus_analysis = periodic_analysis .OR. (translation_analysis .AND. ANY(cell%perd /= 0))
      CALL section_vals_val_get(section, 'MP_GRID', i_vals=grid, explicit=explicit_grid)
      IF (ANY(grid < 1)) CPABORT('Quadratic MP_GRID must be positive.')
      IF (.NOT. torus_analysis .AND. explicit_grid) THEN
         CPABORT('Quadratic MP_GRID requires a periodic property analysis.')
      END IF
      IF (ANY(grid /= 1 .AND. cell%perd == 0)) THEN
         CPABORT('Quadratic MP_GRID must be one in open directions.')
      END IF
      CALL section_vals_val_get(input, 'DFT%POISSON%PERIODIC', i_val=periodic)
      IF (torus_analysis) THEN
         CALL section_vals_val_get(input, 'SUBSYS%CELL%PERIODIC', i_val=cell_periodic)
         IF (cell_periodic /= periodic .OR. COUNT(cell%perd /= 0) == 0) THEN
            CPABORT('Quadratic torus needs matching, nonempty CELL/POISSON PERIODIC.')
         END IF
         directions = PACK([1, 2, 3], cell%perd /= 0)
         ndim = SIZE(directions)
         IF (periodic_analysis) THEN
            nobs = 1 + 2*ndim
            ALLOCATE (radii(ndim), angles(ndim))
            DO i = 1, ndim
               axis = directions(i)
               radii(i) = REAL(grid(axis), dp)/(twopi*SQRT(SUM(cell%h_inv(axis, :)**2)))
            END DO
         END IF
      ELSE
         IF (periodic /= use_perd_none .OR. ANY(cell%perd /= 0) .OR. do_kpoints) THEN
            CPABORT('Finite quadratic analysis requires PERIODIC NONE without KPOINTS.')
         END IF
         nobs = 4
      END IF
      IF (translation_analysis) THEN
         CALL section_vals_val_get(section, 'TRANSLATION_VECTOR', n_rep_val=nt)
         IF (nt < 1) THEN
            CPABORT('Quadratic TRANSLATION requires at least one TRANSLATION_VECTOR.')
         END IF
         nobs = 1 + nt
         ALLOCATE (translations(3, nt), complex_shift(nobs), unitary(nobs))
         unitary = .TRUE.
         unitary(1) = .FALSE.
         DO i = 1, nt
            CALL section_vals_val_get(section, 'TRANSLATION_VECTOR', i_rep_val=i, r_vals=point)
            IF (.NOT. ALL(ieee_is_finite(point))) THEN
               CPABORT('Quadratic TRANSLATION_VECTOR must be finite.')
            END IF
            translations(:, i) = point
            CALL cell_transform_input_cartesian(cell, translations(:, i))
         END DO
         CALL section_vals_val_get(section, 'POSITION', explicit=explicit_value)
         IF (explicit_value) THEN
            CPABORT('Quadratic TRANSLATION uses WAVE_VECTOR, not POSITION.')
         END IF
         CALL section_vals_val_get(section, 'KAPPA', explicit=explicit_value)
         IF (explicit_value) THEN
            CPABORT('Quadratic TRANSLATION uses TRANSLATION_SCALE, not KAPPA.')
         END IF
      ELSE
         CALL section_vals_val_get(section, 'TRANSLATION_VECTOR', explicit=explicit_value)
         IF (explicit_value) THEN
            CPABORT('TRANSLATION_VECTOR requires quadratic FORMULATION TRANSLATION.')
         END IF
         CALL section_vals_val_get(section, 'WAVE_VECTOR', explicit=explicit_value)
         IF (explicit_value) THEN
            CPABORT('WAVE_VECTOR requires quadratic FORMULATION TRANSLATION.')
         END IF
         CALL section_vals_val_get(section, 'TRANSLATION_SCALE', explicit=explicit_value)
         IF (explicit_value) THEN
            CPABORT('TRANSLATION_SCALE requires quadratic FORMULATION TRANSLATION.')
         END IF
      END IF
      ALLOCATE (operators(nobs), shift(nobs), weights(nobs))
      CALL section_vals_val_get(section, 'SPIN_CHANNEL', i_val=spin)
      IF (spin < 1 .OR. spin > SIZE(matrix_ks)) THEN
         CPABORT('Quadratic SPIN_CHANNEL is unavailable.')
      END IF
      CALL section_vals_val_get(section, 'SOC', l_val=do_soc)
      IF (do_soc .AND. spin /= 1) CPABORT('Quadratic SOC requires SPIN_CHANNEL 1.')
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)
      CALL section_vals_val_get(section, 'MAX_AO', i_val=max_ao)
      IF (nao > max_ao .OR. max_ao < 1) CPABORT('QUADRATIC_PSEUDOSPECTRUM exceeds MAX_AO.')
      DO axis = 1, 3
         IF (grid(axis) > max_ao/nao) CPABORT('Quadratic full AO torus exceeds MAX_AO.')
         nao = nao*grid(axis)
      END DO
      n = nao
      IF (do_soc) n = 2*nao
      CALL section_vals_val_get(section, 'NSTATES', i_val=nstates)
      max_space = n
      max_iter = 1
      CALL section_vals_val_get(section, 'EPS_EIGEN', r_val=eps_eigen)
      CALL section_vals_val_get(section, 'EPS_METRIC', r_val=eps_metric)
      CALL section_vals_val_get(section, 'ENERGY', r_vals=energies)
      CALL section_vals_val_get(section, 'KAPPA', r_vals=kappas)
      CALL section_vals_val_get(section, 'POSITION', n_rep_val=npositions)
      IF (translation_analysis) THEN
         CALL section_vals_val_get(section, 'TRANSLATION_SCALE', r_vals=kappas)
         CALL section_vals_val_get(section, 'WAVE_VECTOR', n_rep_val=npositions)
      END IF
      IF (nstates < 1 .OR. nstates > n .OR. max_iter < 1) THEN
         CPABORT('Invalid quadratic state/iteration limits.')
      END IF
      IF (.NOT. ALL(ieee_is_finite([energies, kappas, eps_eigen, eps_metric]))) THEN
         CPABORT('Quadratic analysis requires finite parameters.')
      END IF
      IF (ANY(kappas <= 0.0_dp) .OR. eps_eigen <= 0.0_dp .OR. eps_metric <= 0.0_dp .OR. eps_metric >= 1.0_dp) THEN
         CPABORT('Invalid quadratic scales or tolerances.')
      END IF
      CALL cite_reference(Cerjan2023)
      IF (translation_analysis) CALL cite_reference(Bairnsfather2026)
      NULLIFY (positions, soc)
      IF (periodic_analysis) THEN
         CALL localizer_torus_operators(qs_env, grid, directions, do_soc, matrix_ks, matrix_s, cosines, sines, soc)
         ALLOCATE (positions(2*ndim))
         DO i = 1, ndim
            positions(2*i - 1)%matrix => cosines(i)%matrix
            positions(2*i)%matrix => sines(directions(i))%matrix
         END DO
      ELSE IF (translation_analysis) THEN
         IF (torus_analysis) THEN
            CALL torus_translation_operators(qs_env, grid, translations, do_soc, matrix_ks, matrix_s, positions, soc)
         ELSE
            CALL ao_translation_matrices(qs_env, translations, positions)
         END IF
      ELSE
         CALL finite_ao_positions(qs_env, positions)
      END IF
      IF (do_soc .AND. .NOT. torus_analysis) THEN
         CALL finite_ao_soc(qs_env, soc)
      END IF
      IF (do_soc) THEN
         CALL localizer_ao_operators(matrix_ks(spin)%matrix, matrix_s(1)%matrix, positions, operators, metric, soc)
         CALL dbcsr_deallocate_matrix_set(soc)
      ELSE
         CALL localizer_ao_operators(matrix_ks(spin)%matrix, matrix_s(1)%matrix, positions, operators, metric)
      END IF
      IF (periodic_analysis) THEN
         DEALLOCATE (positions)
         CALL dbcsr_deallocate_matrix_set(cosines)
         CALL dbcsr_deallocate_matrix_set(sines)
      ELSE
         CALL dbcsr_deallocate_matrix_set(positions)
      END IF
      IF (torus_analysis) THEN
         CALL dbcsr_deallocate_matrix_set(matrix_ks)
         CALL dbcsr_deallocate_matrix_set(matrix_s)
      END IF
      CALL localizer_dbcsr_check(metric, para_env, 1.e-10_dp, valid, max_hermiticity_error)
      IF (.NOT. valid) CPABORT('Quadratic AO metric is nonfinite or not Hermitian.')
      DO i = 1, nobs
         IF (translation_analysis .AND. i > 1) CYCLE
         CALL localizer_dbcsr_check(operators(i), para_env, 1.e-10_dp, valid, hermiticity_error)
         max_hermiticity_error = MAX(max_hermiticity_error, hermiticity_error)
         IF (.NOT. valid) THEN
            WRITE (matrix_error, '(A,I0,A,ES12.4,A)') 'Quadratic AO operator ', i, &
               ' is nonfinite or not Hermitian (relative error ', hermiticity_error, &
               '). Check QS%EPS_PPNL/EPS_PGF_ORB and integration thresholds; no symmetrization is applied.'
            CPABORT(TRIM(matrix_error))
         END IF
      END DO
      ALLOCATE (dense_operators(n, n, nobs))
      CALL localizer_dbcsr_dense(metric, para_env, dense_metric)
      DO i = 1, nobs
         CALL localizer_dbcsr_dense(operators(i), para_env, temporary)
         dense_operators(:, :, i) = temporary
         DEALLOCATE (temporary)
      END DO
      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)
      iw = cp_print_key_unit_nr(logger, section, 'STATES', extension='.states')
      IF (iw > 0) WRITE (iw, '(A)') '# Quadratic eigenstates: S-normalized AO coefficients, spin-major for SOC'
      IF (iw > 0 .AND. torus_analysis) WRITE (iw, '(A,3I8)') '# Cell-major AO torus within each spin block; grid ', grid
      IF (output_unit > 0) THEN
         WRITE (output_unit, '(/,T2,A,I0)') 'QUADRATIC_PSEUDOSPECTRUM| AO/spinor dimension: ', n
         WRITE (output_unit, '(T2,A,ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Relative Hermiticity error: ', max_hermiticity_error
         WRITE (output_unit, '(T2,A)') 'QUADRATIC_PSEUDOSPECTRUM| Dense finite-basis reference'
         WRITE (output_unit, '(T2,A)') 'QUADRATIC_PSEUDOSPECTRUM| Joint approximate states, not a topology index.'
         IF (translation_analysis) THEN
            WRITE (output_unit, '(T2,A)') 'QUADRATIC_PSEUDOSPECTRUM| Continuum translations U(a) psi(r) = psi(r+a).'
            WRITE (output_unit, '(T2,A)') 'QUADRATIC_PSEUDOSPECTRUM| Full translation norms include loss outside the AO span.'
            WRITE (output_unit, '(T2,A)') 'QUADRATIC_PSEUDOSPECTRUM| Energy residual still uses the projected AO Hamiltonian.'
            DO i = 1, nt
               WRITE (output_unit, '(T2,A,I0,A,3ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Translation ', i, &
                  ' internal [bohr]: ', translations(:, i)
            END DO
         END IF
         IF (torus_analysis) THEN
            WRITE (output_unit, '(T2,A,3I6)') 'QUADRATIC_PSEUDOSPECTRUM| Full analysis mesh: ', grid
            WRITE (output_unit, '(T2,A,3I4)') 'QUADRATIC_PSEUDOSPECTRUM| Periodic lattice directions: ', directions
            IF (translation_analysis) THEN
               WRITE (output_unit, '(T2,A)') 'QUADRATIC_PSEUDOSPECTRUM| Periodized translations on the full AO torus.'
            ELSE
               WRITE (output_unit, '(T2,A)') 'QUADRATIC_PSEUDOSPECTRUM| Torus chord metric; open directions are not localized.'
            END IF
         END IF
      END IF
      DO ip = 1, npositions
         IF (translation_analysis) THEN
            CALL section_vals_val_get(section, 'WAVE_VECTOR', i_rep_val=ip, r_vals=point)
         ELSE
            CALL section_vals_val_get(section, 'POSITION', i_rep_val=ip, r_vals=point)
         END IF
         IF (.NOT. ALL(ieee_is_finite(point))) THEN
            CPABORT('Quadratic POSITION or WAVE_VECTOR must be finite.')
         END IF
         IF (translation_analysis) THEN
            query = point
            CALL cell_transform_input_reciprocal(cell, query)
            DO i = 1, nt
               complex_shift(i + 1) = EXP(CMPLX(0.0_dp, DOT_PRODUCT(query, translations(:, i)), dp))
            END DO
         ELSE IF (periodic_analysis) THEN
            query = point
            CALL cell_transform_input_cartesian(cell, query)
            DO i = 1, ndim
               axis = directions(i)
               angles(i) = twopi*MODULO(DOT_PRODUCT(cell%h_inv(axis, :), query)/REAL(grid(axis), dp), 1.0_dp)
               shift(2*i) = COS(angles(i))
               shift(2*i + 1) = SIN(angles(i))
            END DO
         ELSE
            shift(2:) = point
         END IF
         DO ik = 1, SIZE(kappas)
            DO ie = 1, SIZE(energies)
               shift(1) = energies(ie)
               weights = kappas(ik)
               weights(1) = 1.0_dp
               IF (periodic_analysis) THEN
                  DO i = 1, ndim
                     weights(2*i:2*i + 1) = kappas(ik)*radii(i)
                  END DO
               END IF
               IF (para_env%is_source()) THEN
                  IF (translation_analysis) THEN
                     complex_shift(1) = energies(ie)
                     CALL quadratic_dense(dense_operators, dense_metric, complex_shift, weights, nstates, .FALSE., max_space, &
                                          max_iter, eps_eigen, eps_metric, RESULT, status, unitary)
                  ELSE
                     CALL quadratic_dense(dense_operators, dense_metric, shift, weights, nstates, .FALSE., max_space, &
                                          max_iter, eps_eigen, eps_metric, RESULT, status)
                  END IF
               END IF
               CALL para_env%bcast(status, para_env%source)
               IF (status /= 0) THEN
                  IF (status == 1) THEN
                     IF (output_unit > 0) THEN
                        WRITE (output_unit, '(T2,A,I0,A,ES23.14)') &
                           'QUADRATIC_PSEUDOSPECTRUM| Unconverged after ', result%iterations, &
                           ' iterations; maximum residual: ', MAXVAL(result%residuals)
                     END IF
                     CPABORT('Quadratic eigensolver unconverged.')
                  END IF
                  CPABORT('Quadratic analysis failed numerical checks.')
               END IF
               IF (output_unit > 0) THEN
                  IF (translation_analysis) THEN
                     WRITE (output_unit, '(T2,A,3ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Wave vector [bohr^-1]: ', point
                     WRITE (output_unit, '(T2,A,ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Translation scale [hartree]: ', kappas(ik)
                  ELSE
                     WRITE (output_unit, '(T2,A,3ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Position [bohr]: ', point
                     WRITE (output_unit, '(T2,A,ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Kappa [hartree/bohr]: ', kappas(ik)
                  END IF
                  WRITE (output_unit, '(T2,A,ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Energy [hartree]: ', energies(ie)
                  WRITE (output_unit, '(T2,A,I0)') 'QUADRATIC_PSEUDOSPECTRUM| Ritz iterations: ', result%iterations
                  DO j = 1, nstates
                     WRITE (output_unit, '(T2,A,I0)') 'QUADRATIC_PSEUDOSPECTRUM| State: ', j
                     WRITE (output_unit, '(T2,A,ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Gap [hartree]: ', &
                        SQRT(result%eigenvalues(j))
                     WRITE (output_unit, '(T2,A,ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Energy residual [hartree]: ', &
                        result%errors(1, j)
                     IF (translation_analysis) THEN
                        DO i = 2, nobs
                           WRITE (output_unit, '(T2,A,I0,A,3ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Translation residual ', i - 1, &
                              ' full/projected/leakage: ', result%errors(i, j), result%projected_errors(i, j), result%leakage(i, j)
                          WRITE (output_unit, '(T2,A,I0,A,2ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Translation expectation ', i - 1, &
                              ' real/imag: ', result%complex_means(i, j)
                        END DO
                     ELSE IF (periodic_analysis) THEN
                        DO i = 1, ndim
                           WRITE (output_unit, '(T2,A,I1,A,ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Chord residual axis ', &
                              directions(i), ' [bohr]: ', radii(i)*SQRT(SUM(result%errors(2*i:2*i + 1, j)**2))
                           WRITE (output_unit, '(T2,A,I1,A,2ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Mean cos/sin axis ', &
                              directions(i), ': ', result%means(2*i:2*i + 1, j)
                        END DO
                     ELSE
                        WRITE (output_unit, '(T2,A,3ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Position residuals [bohr]: ', &
                           result%errors(2:4, j)
                     END IF
                     WRITE (output_unit, '(T2,A,ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Eigenpair residual [hartree^2]: ', &
                        result%residuals(j)
                     WRITE (output_unit, '(T2,A,ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Mean energy [hartree]: ', result%means(1, j)
                     IF (.NOT. periodic_analysis .AND. .NOT. translation_analysis) THEN
                        WRITE (output_unit, '(T2,A,3ES23.14)') 'QUADRATIC_PSEUDOSPECTRUM| Mean position [bohr]: ', &
                           result%means(2:4, j)
                     END IF
                  END DO
               END IF
               IF (iw > 0) THEN
                  IF (translation_analysis) THEN
                     WRITE (iw, '(A,3I8,A,5ES24.15)') '# Query ', ip, ie, ik, ' E/kx/ky/kz/scale ', energies(ie), point, kappas(ik)
                  ELSE
                     WRITE (iw, '(A,3I8,A,5ES24.15)') '# Query ', ip, ie, ik, ' E/x/y/z/kappa ', energies(ie), point, kappas(ik)
                  END IF
                  DO j = 1, nstates
                     WRITE (iw, '(A,I8,A,ES24.15)') '# State ', j, ' lambda [hartree^2] ', result%eigenvalues(j)
                     DO ao = 1, n
                        WRITE (iw, '(I10,2ES25.16)') ao, REAL(result%states(ao, j), dp), AIMAG(result%states(ao, j))
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END DO
      CALL cp_print_key_finished_output(iw, logger, section, 'STATES')
      DO i = 1, nobs
         CALL localizer_dbcsr_release(operators(i))
      END DO
      CALL localizer_dbcsr_release(metric)
      CALL timestop(handle)
   END SUBROUTINE qs_scf_post_quadratic
END MODULE qs_quadratic_pseudospectrum
