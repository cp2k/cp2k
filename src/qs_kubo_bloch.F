!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Full-mesh scalar/spinor Bloch currents from CP2K image operators and analytic AO moments.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE qs_kubo_bloch
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_info,&
                                              dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE input_constants,                 ONLY: do_method_gapw,&
                                              do_method_gpw
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_k_r_trafo_simple,         ONLY: replicate_rs_matrices,&
                                              rs_to_kp
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_release,&
                                              kpoint_type
   USE kubo_bloch,                      ONLY: bloch_current,&
                                              bloch_current_from_frame
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_release,&
                                              qs_environment_type
   USE qs_finite_ao,                    ONLY: finite_ao_soc
   USE qs_gamma2kp,                     ONLY: create_kp_from_gamma
   USE qs_mo_types,                     ONLY: mo_set_type
   USE qs_moments,                      ONLY: build_local_moment_matrix_rs_img
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_property_symmetry,            ONLY: property_symmetry_frame,&
                                              property_symmetry_plan
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kubo_bloch'
   PUBLIC :: setup_bloch_transport
CONTAINS
! **************************************************************************************************
!> \brief Build a full Gamma-centered property mesh, independent of SCF reduction.
!> \param qs_env Converged GPW/GAPW environment
!> \param grid Full property mesh; nonperiodic directions must be one
!> \param do_soc Add complete GTH SOC after a restricted SCF
!> \param projection Cartesian projection onto periodic transport directions
!> \param max_memory_mb Safety bound on replicated array memory per rank, in MiB
!> \param energies Complete spectrum, flattened (spin,k) channels
!> \param current Anti-Hermitian -i*velocity per channel
!> \param occupations Maximum occupations including normalized k weights
!> \param symmetry Reconstruct symmetry-equivalent property eigenframes
!> \param backend Geometric symmetry backend
!> \param eps_symmetry Geometric tolerance
!> \param eps_operators Target frame and current validation tolerance
!> \note This dense reference stores image operators on every rank; eigensystems are distributed by k.
! **************************************************************************************************
   SUBROUTINE setup_bloch_transport(qs_env, grid, do_soc, projection, max_memory_mb, energies, current, occupations, &
                                    symmetry, backend, eps_symmetry, eps_operators)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: grid(3)
      LOGICAL, INTENT(IN)                                :: do_soc
      REAL(dp), INTENT(IN)                               :: projection(3, 3)
      INTEGER, INTENT(IN)                                :: max_memory_mb
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: energies(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: current(:, :, :, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: occupations(:)
      LOGICAL, INTENT(IN)                                :: symmetry
      INTEGER, INTENT(IN)                                :: backend
      REAL(dp), INTENT(IN)                               :: eps_symmetry, eps_operators

      COMPLEX(dp), ALLOCATABLE                           :: a(:, :), cur(:, :, :), d(:, :, :), &
                                                            dh(:, :, :), dsk(:, :, :), h(:, :), &
                                                            s(:, :), v(:, :, :), source(:, :), TARGET(:, :), &
                                                            source_current(:, :, :), expected(:, :, :)
      INTEGER                                            :: axis, channel, failed, idir, ik, n, nao, &
                                                            ni, nk, nspin, spin, status, ir, nr, iw, &
                                                            trial, ntry, operation_try, alternatives
      INTEGER, ALLOCATABLE                              :: representative(:), operation(:)
      INTEGER(int_8)                                     :: bytes
      LOGICAL                                            :: do_kpoints, antiunitary
      REAL(dp)                                           :: k(3), rotation(3, 3), residuals(2), errors(3), &
                                                            candidate_errors(3), best_errors(3)
      REAL(dp), ALLOCATABLE                              :: mesh(:, :), representative_k(:, :), source_energies(:)
      REAL(dp), ALLOCATABLE                              :: ds(:, :, :, :), hs(:, :, :, :), &
                                                            ss(:, :, :, :), vs(:, :, :, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), POINTER                        :: dr(:, :), hr(:, :), sr(:, :), vr(:, :)
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: images, scf, property_points
      TYPE(mo_set_type), POINTER                         :: mos(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), POINTER            :: all_pairs(:)
      TYPE(qs_environment_type), POINTER                 :: qp

      CALL get_qs_env(qs_env, do_kpoints=do_kpoints, dft_control=dft_control)
      SELECT CASE (dft_control%qs_control%method_id)
      CASE (do_method_gpw, do_method_gapw)
      CASE DEFAULT
         CPABORT('Kubo BLOCH currently supports GPW/GAPW only.')
      END SELECT
      NULLIFY (qp, images, dr, vr, property_points)
      IF (do_kpoints) THEN
         qp => qs_env
      ELSE
         CALL create_kp_from_gamma(qs_env, qp)
      END IF
      CALL get_qs_env(qp, cell=cell, para_env=para_env, matrix_ks_kp=hr, matrix_s_kp=sr, &
                      sab_all=all_pairs, kpoints=scf, mos=mos)
      CALL dbcsr_get_info(sr(1, 1)%matrix, nfullrows_total=nao)
      n = nao
      IF (do_soc) n = 2*nao
      nspin = SIZE(hr, 1)
      nk = PRODUCT(grid)
      CALL kpoint_create(images)
      CALL kpoint_init_cell_index(images, all_pairs, para_env, ni)
      bytes = 8_int_8*nao*nao*ni*(nspin + 4 + MERGE(3, 0, do_soc)) + &
              16_int_8*n*n*(3_int_8*nspin*nk + 48_int_8)
      IF (max_memory_mb < 1 .OR. bytes > INT(max_memory_mb, int_8)*1024*1024) THEN
         CPABORT('Kubo BLOCH arrays exceed MAX_MEMORY_MB per MPI rank.')
      END IF
      ALLOCATE (dr(3, ni))
      CALL build_local_moment_matrix_rs_img(qp, dr)
      ALLOCATE (hs(nspin, nao, nao, ni), ss(1, nao, nao, ni), ds(3, nao, nao, ni), SOURCE=0.0_dp)
      CALL replicate_rs_matrices(hr, scf, hs, images%cell_to_index)
      CALL replicate_rs_matrices(sr, scf, ss, images%cell_to_index)
      CALL replicate_rs_matrices(dr, images, ds, images%cell_to_index)
      CALL dbcsr_deallocate_matrix_set(dr)
      IF (do_soc) THEN
         CALL finite_ao_soc(qp, vr)
         ALLOCATE (vs(3, nao, nao, ni), SOURCE=0.0_dp)
         CALL replicate_rs_matrices(vr, scf, vs, images%cell_to_index)
         CALL dbcsr_deallocate_matrix_set(vr)
      END IF
      ALLOCATE (energies(n, nspin*nk), current(n, n, 3, nspin*nk), occupations(nspin*nk))
      ALLOCATE (h(n, n), s(n, n), dh(n, n, 3), dsk(n, n, 3), d(n, n, 3), cur(n, n, 3), &
                a(nao, nao), v(nao, nao, 3))
      ALLOCATE (mesh(3, nk), source(n, n), TARGET(n, n), source_energies(n), &
                source_current(n, n, 3), expected(n, n, 3))
      DO ik = 1, nk
         mesh(:, ik) = REAL([MOD(ik - 1, grid(1)), MOD((ik - 1)/grid(1), grid(2)), &
                             (ik - 1)/(grid(1)*grid(2))], dp)/REAL(grid, dp)
      END DO
      IF (symmetry) THEN
         CALL property_symmetry_plan(qp, mesh, grid, projection, backend, eps_symmetry, &
                                     property_points, representative, operation)
         nr = property_points%nkp
         ALLOCATE (representative_k(3, nr))
         representative_k(:, :) = property_points%xkp
      ELSE
         nr = nk
         ALLOCATE (representative_k(3, nr))
         representative_k(:, :) = mesh
         ALLOCATE (representative(nk), operation(nk))
         representative(:) = [(ik, ik=1, nk)]
         operation = 0
      END IF
      energies = 0.0_dp
      current = 0.0_dp
      failed = 0
      errors = 0.0_dp
      alternatives = 0
      DO spin = 1, nspin
         DO ik = 1, nk
            channel = ik + (spin - 1)*nk
            occupations(channel) = mos(spin)%maxocc/REAL(nk, dp)
            IF (do_soc) occupations(channel) = 1.0_dp/REAL(nk, dp)
         END DO
         DO ir = 1, nr
            IF (MOD(ir - 1, para_env%num_pe) /= para_env%mepos) CYCLE
            k = representative_k(:, ir)
            CALL operators_at_k()
            CALL bloch_current(h, s, dh, dsk, d, source_energies, source_current, status, source)
            IF (status /= 0) THEN
               failed = MAX(failed, ABS(status))
               CYCLE
            END IF
            DO ik = 1, nk
               IF (representative(ik) /= ir) CYCLE
               channel = ik + (spin - 1)*nk
               energies(:, channel) = source_energies
               cur(:, :, :) = source_current
               IF (operation(ik) /= 0) THEN
                  k = mesh(:, ik)
                  CALL operators_at_k()
                  ntry = property_points%kp_sym(ir)%kpoint_sym%nwred
                  best_errors = HUGE(1.0_dp)
                  DO trial = 0, ntry
                     operation_try = operation(ik)
                     IF (trial > 0) operation_try = trial
                     IF (trial > 0 .AND. operation_try == operation(ik)) CYCLE
                     CALL property_symmetry_frame(property_points, ir, operation_try, k, source, do_soc, TARGET, &
                                                  rotation, antiunitary, status)
                     IF (status /= 0) CYCLE
                     IF (MAXVAL(ABS(MATMUL(rotation, MATMUL(projection, TRANSPOSE(rotation))) - projection)) &
                         > 1.e-9_dp) CYCLE
                     CALL bloch_current_from_frame(h, s, dh, dsk, d, source_energies, TARGET, cur, &
                                                   eps_operators, residuals, status)
                     candidate_errors = [residuals, 0.0_dp]
                     IF (status == 0) THEN
                        expected = 0.0_dp
                        DO axis = 1, 3
                           DO idir = 1, 3
                              IF (antiunitary) THEN
                                 expected(:, :, axis) = expected(:, :, axis) + &
                                                        rotation(axis, idir)*CONJG(source_current(:, :, idir))
                              ELSE
                                 expected(:, :, axis) = expected(:, :, axis) + rotation(axis, idir)*source_current(:, :, idir)
                              END IF
                           END DO
                        END DO
                        candidate_errors(3) = MAXVAL(ABS(cur - expected))/MAX(1.0_dp, MAXVAL(ABS(cur)))
                     END IF
                     IF (MAXVAL(candidate_errors) < MAXVAL(best_errors)) best_errors = candidate_errors
                     IF (status /= 0 .OR. MAXVAL(candidate_errors) > eps_operators) CYCLE
                     best_errors = candidate_errors
                     IF (trial > 0) alternatives = alternatives + 1
                     EXIT
                  END DO
                  errors = MAX(errors, best_errors)
                  IF (trial > ntry) THEN
                     failed = MAX(failed, 5)
                     CYCLE
                  END IF
               END IF
               DO axis = 1, 3
                  DO idir = 1, 3
                     current(:, :, axis, channel) = current(:, :, axis, channel) + projection(axis, idir)*cur(:, :, idir)
                  END DO
               END DO
            END DO
         END DO
      END DO
      CALL para_env%max(failed)
      CALL para_env%max(errors)
      CALL para_env%sum(alternatives)
      iw = cp_logger_get_default_io_unit()
      IF (symmetry .AND. iw > 0) THEN
         WRITE (iw, '(T2,A,I0,A,I0)') 'KUBO_TRANSPORT| Property symmetry representatives: ', nr, ' / ', nk
         WRITE (iw, '(T2,A,3ES14.5)') 'KUBO_TRANSPORT| Symmetry errors S/H/current: ', errors
         WRITE (iw, '(T2,A,I0)') 'KUBO_TRANSPORT| Alternative validated orbit operations: ', alternatives
      END IF
      IF (failed /= 0) THEN
         IF (symmetry) THEN
            CPABORT('Bloch Kubo property symmetry failed target metric/Hamiltonian/current covariance.')
         ELSE
            CPABORT('Bloch Kubo failed metric, Hermiticity or AO-connection validation.')
         END IF
      END IF
      CALL para_env%sum(energies)
      CALL para_env%sum(current)
      CALL kpoint_release(images)
      IF (ASSOCIATED(property_points)) CALL kpoint_release(property_points)
      IF (.NOT. do_kpoints) THEN
         CALL qs_env_release(qp)
         DEALLOCATE (qp)
      END IF
   CONTAINS
! **************************************************************************************************
!> \brief Fourier-transform the frozen image operators at the current (k,spin).
! **************************************************************************************************
      SUBROUTINE operators_at_k()
      INTEGER                                            :: component, direction

         CALL rs_to_kp(hs(spin, :, :, :), a, images%index_to_cell, k)
         CALL spin_diagonal(a, h)
         IF (do_soc) THEN
            DO component = 1, 3
               CALL rs_to_kp(vs(component, :, :, :), v(:, :, component), images%index_to_cell, k)
            END DO
            CALL add_soc(h, v)
         END IF
         CALL rs_to_kp(ss(1, :, :, :), a, images%index_to_cell, k)
         CALL spin_diagonal(a, s)
         DO component = 1, 3
            CALL rs_to_kp(hs(spin, :, :, :), a, images%index_to_cell, k, component, cell%hmat)
            CALL spin_diagonal(a, dh(:, :, component))
            IF (do_soc) THEN
               DO direction = 1, 3
                  CALL rs_to_kp(vs(direction, :, :, :), v(:, :, direction), images%index_to_cell, k, component, cell%hmat)
               END DO
               CALL add_soc(dh(:, :, component), v)
            END IF
            CALL rs_to_kp(ss(1, :, :, :), a, images%index_to_cell, k, component, cell%hmat)
            CALL spin_diagonal(a, dsk(:, :, component))
            CALL rs_to_kp(ds(component, :, :, :), a, images%index_to_cell, k)
            CALL spin_diagonal(a, d(:, :, component))
         END DO
      END SUBROUTINE operators_at_k
   END SUBROUTINE setup_bloch_transport
! **************************************************************************************************
!> \brief Duplicate a scalar operator in the physical spin basis.
!> \param a Scalar operator
!> \param out Scalar or doubled spinor operator
! **************************************************************************************************
   SUBROUTINE spin_diagonal(a, out)
      COMPLEX(dp), INTENT(IN)                            :: a(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: out(:, :)

      INTEGER                                            :: n

      n = SIZE(a, 1)
      out = 0.0_dp
      out(:n, :n) = a
      IF (SIZE(out, 1) == 2*n) out(n + 1:, n + 1:) = a
   END SUBROUTINE spin_diagonal
! **************************************************************************************************
!> \brief Add i*sigma times the real-spherical-harmonic SOC components, including their derivatives.
!> \param h Spinor Hamiltonian or derivative
!> \param v Fourier-transformed real antisymmetric SOC components (without the factor i)
! **************************************************************************************************
   SUBROUTINE add_soc(h, v)
      COMPLEX(dp), INTENT(INOUT)                         :: h(:, :)
      COMPLEX(dp), INTENT(IN)                            :: v(:, :, :)

      COMPLEX(dp), PARAMETER                             :: iu = (0.0_dp, 1.0_dp)

      INTEGER                                            :: n

      n = SIZE(v, 1)
      h(:n, :n) = h(:n, :n) + iu*v(:, :, 3)
      h(n + 1:, n + 1:) = h(n + 1:, n + 1:) - iu*v(:, :, 3)
      h(:n, n + 1:) = h(:n, n + 1:) + iu*v(:, :, 1) - v(:, :, 2)
      h(n + 1:, :n) = h(n + 1:, :n) + iu*v(:, :, 1) + v(:, :, 2)
   END SUBROUTINE add_soc
END MODULE qs_kubo_bloch
