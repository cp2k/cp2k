!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Finite-window AO Chern marker after a converged Quickstep calculation.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE qs_real_space_chern
   USE bibliography,                    ONLY: BiancoResta2011,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type,&
                                              use_perd_none
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_info,&
                                              dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: do_method_gapw,&
                                              do_method_gpw
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE kubo_projected,                  ONLY: projected_eigensystem
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_finite_ao,                    ONLY: finite_ao_positions,&
                                              finite_ao_soc,&
                                              finite_ao_window
   USE spectral_localizer_dbcsr,        ONLY: localizer_ao_operators,&
                                              localizer_dbcsr_dense,&
                                              localizer_dbcsr_release,&
                                              localizer_matrix_type
   USE topology_real_space,             ONLY: finite_chern_marker
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_real_space_chern'
   PUBLIC :: qs_scf_post_real_space_chern
CONTAINS
! **************************************************************************************************
!> \brief Integrate a physical Gaussian window and evaluate its metric-correct finite Chern marker.
!> \param qs_env Converged isolated GPW/GAPW environment
! **************************************************************************************************
   SUBROUTINE qs_scf_post_real_space_chern(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_real_space_chern'
      CHARACTER(len=2), PARAMETER :: planes(6) = ['XY', 'YX', 'XZ', 'ZX', 'YZ', 'ZY']
      INTEGER, PARAMETER :: axes(2, 6) = RESHAPE([1, 2, 2, 1, 1, 3, 3, 1, 2, 3, 3, 2], [2, 6])

      CHARACTER(len=256)                                 :: error
      COMPLEX(dp), ALLOCATABLE                           :: c(:, :), gram(:, :), h(:, :), s(:, :), &
                                                            temporary(:, :), w(:, :), x(:, :, :)
      INTEGER                                            :: handle, i, max_ao, n, nao, nocc, &
                                                            output_unit, periodic, plane, spin, &
                                                            status
      LOGICAL                                            :: active, do_kpoints, do_soc, &
                                                            explicit_value
      REAL(dp)                                           :: alpha(3), area, center(3), energy, &
                                                            eps_gap, gap, marker, weight, width(2)
      REAL(dp), ALLOCATABLE                              :: spectrum(:)
      REAL(dp), POINTER                                  :: values(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type)                                 :: observables(3)
      TYPE(dbcsr_p_type), POINTER                        :: images_h(:, :), images_s(:, :), &
                                                            positions(:), soc(:, :), window(:)
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(localizer_matrix_type)                        :: metric, operators(4)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: input, section

      CALL get_qs_env(qs_env, input=input)
      section => section_vals_get_subs_vals(input, 'PROPERTIES%REAL_SPACE_CHERN')
      CALL section_vals_val_get(section, '_SECTION_PARAMETERS_', l_val=active)
      IF (.NOT. active) RETURN
      CALL timeset(routineN, handle)
      IF (.NOT. qs_env%scf_convergence_available .OR. .NOT. qs_env%scf_converged) THEN
         CPABORT('REAL_SPACE_CHERN requires converged SCF.')
      END IF
      CALL get_qs_env(qs_env, cell=cell, do_kpoints=do_kpoints, dft_control=dft_control, &
                      matrix_s_kp=images_s, matrix_ks_kp=images_h, para_env=para_env)
      CALL section_vals_val_get(input, 'DFT%POISSON%PERIODIC', i_val=periodic)
      IF (periodic /= use_perd_none .OR. ANY(cell%perd /= 0) .OR. do_kpoints) THEN
         CPABORT('REAL_SPACE_CHERN requires CELL/POISSON PERIODIC NONE without KPOINTS.')
      END IF
      SELECT CASE (dft_control%qs_control%method_id)
      CASE (do_method_gpw, do_method_gapw)
      CASE DEFAULT
         CPABORT('REAL_SPACE_CHERN supports GPW/GAPW only.')
      END SELECT
      CALL section_vals_val_get(section, 'PLANE', i_val=plane)
      CALL section_vals_val_get(section, 'WINDOW_CENTER', r_vals=values)
      center = values
      CALL section_vals_val_get(section, 'WINDOW_SIGMA', explicit=explicit_value)
      IF (.NOT. explicit_value) THEN
         CPABORT('REAL_SPACE_CHERN requires explicit WINDOW_SIGMA.')
      END IF
      CALL section_vals_val_get(section, 'WINDOW_SIGMA', r_vals=values)
      width = values
      IF (.NOT. ALL(ieee_is_finite(center)) .OR. .NOT. ALL(ieee_is_finite(width))) THEN
         CPABORT('REAL_SPACE_CHERN window parameters must be finite.')
      END IF
      IF (ANY(width <= SQRT(TINY(1.0_dp))) .OR. ANY(width >= SQRT(HUGE(1.0_dp)/twopi))) THEN
         CPABORT('REAL_SPACE_CHERN WINDOW_SIGMA must be positive and numerically representable.')
      END IF
      alpha = 0.0_dp
      alpha(axes(:, plane)) = (0.5_dp/width)/width
      area = twopi*width(1)*width(2)
      CALL section_vals_val_get(section, 'ENERGY', explicit=explicit_value)
      IF (.NOT. explicit_value) THEN
         CPABORT('REAL_SPACE_CHERN requires explicit ENERGY in a resolved gap.')
      END IF
      CALL section_vals_val_get(section, 'ENERGY', r_val=energy)
      CALL section_vals_val_get(section, 'EPS_GAP', r_val=eps_gap)
      IF (.NOT. ieee_is_finite(energy) .OR. .NOT. ieee_is_finite(eps_gap)) THEN
         CPABORT('REAL_SPACE_CHERN ENERGY/EPS_GAP must be finite.')
      END IF
      IF (eps_gap < 0.0_dp) CPABORT('REAL_SPACE_CHERN EPS_GAP must be nonnegative.')
      CALL section_vals_val_get(section, 'SPIN_CHANNEL', i_val=spin)
      CALL section_vals_val_get(section, 'SOC', l_val=do_soc)
      IF (spin < 1 .OR. spin > SIZE(images_h, 1)) THEN
         CPABORT('REAL_SPACE_CHERN SPIN_CHANNEL is unavailable.')
      END IF
      IF (do_soc .AND. spin /= 1) CPABORT('REAL_SPACE_CHERN SOC requires restricted channel 1.')
      CALL section_vals_val_get(section, 'MAX_AO', i_val=max_ao)
      CALL dbcsr_get_info(images_s(1, 1)%matrix, nfullrows_total=nao)
      IF (max_ao < 1 .OR. nao > max_ao) CPABORT('REAL_SPACE_CHERN exceeds MAX_AO.')
      n = nao
      IF (do_soc) n = 2*nao
      CALL finite_ao_positions(qs_env, positions)
      ALLOCATE (window(1))
      CALL finite_ao_window(qs_env, center, alpha, window(1))
      observables(1:2) = positions(axes(:, plane))
      observables(3) = window(1)
      IF (do_soc) THEN
         CALL finite_ao_soc(qs_env, soc)
         CALL localizer_ao_operators(images_h(spin, 1)%matrix, images_s(1, 1)%matrix, &
                                     observables, operators, metric, soc)
         CALL dbcsr_deallocate_matrix_set(soc)
      ELSE
         CALL localizer_ao_operators(images_h(spin, 1)%matrix, images_s(1, 1)%matrix, &
                                     observables, operators, metric)
      END IF
      CALL dbcsr_deallocate_matrix_set(positions)
      CALL dbcsr_deallocate_matrix_set(window)
      CALL localizer_dbcsr_dense(metric, para_env, s)
      CALL localizer_dbcsr_dense(operators(1), para_env, h)
      ALLOCATE (x(n, n, 2), c(n, n), spectrum(n))
      DO i = 1, 2
         CALL localizer_dbcsr_dense(operators(i + 1), para_env, temporary)
         x(:, :, i) = temporary
         DEALLOCATE (temporary)
      END DO
      CALL localizer_dbcsr_dense(operators(4), para_env, w)
      status = 0
      IF (para_env%is_source()) THEN
         CALL projected_eigensystem(h, s, spectrum, c, status)
         IF (status == 0) THEN
            gap = MINVAL(ABS(spectrum - energy))
            IF (gap <= MAX(eps_gap, 64*EPSILON(1.0_dp)*MAX(1.0_dp, ABS(energy), MAXVAL(ABS(spectrum))))) THEN
               status = -5
            ELSE
               nocc = COUNT(spectrum < energy)
               CALL finite_chern_marker(c(:, :nocc), s, x, w, area, marker, status)
               IF (status == 0) THEN
                  gram = MATMUL(CONJG(TRANSPOSE(c(:, :nocc))), MATMUL(w, c(:, :nocc)))
                  weight = 0.0_dp
                  DO i = 1, nocc
                     weight = weight + REAL(gram(i, i), dp)
                  END DO
               END IF
            END IF
         END IF
      END IF
      CALL para_env%bcast(status)
      IF (status /= 0) THEN
         WRITE (error, '(A,I0,A)') 'REAL_SPACE_CHERN rejected its matrices/projector (status ', status, &
            '). -5: ENERGY intersects a level; -4: window range; -3: normalization; -2: overlap.'
         CPABORT(TRIM(error))
      END IF
      CALL para_env%bcast(nocc)
      CALL para_env%bcast(gap)
      CALL para_env%bcast(marker)
      CALL para_env%bcast(weight)
      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)
      IF (output_unit > 0) THEN
         WRITE (output_unit, '(/,T2,A)') 'REAL_SPACE_CHERN| Finite Gaussian window, not a rounded bulk/Z2 index.'
         WRITE (output_unit, '(T2,A,A)') 'REAL_SPACE_CHERN| Oriented Cartesian plane: ', planes(plane)
         WRITE (output_unit, '(T2,A,3ES23.14)') 'REAL_SPACE_CHERN| Window center [bohr]: ', center
         WRITE (output_unit, '(T2,A,2ES23.14)') 'REAL_SPACE_CHERN| Window sigma [bohr]: ', width
         WRITE (output_unit, '(T2,A,ES23.14)') 'REAL_SPACE_CHERN| Window area [bohr^2]: ', area
         WRITE (output_unit, '(T2,A,ES23.14)') 'REAL_SPACE_CHERN| Projector energy [hartree]: ', energy
         WRITE (output_unit, '(T2,A,ES23.14)') 'REAL_SPACE_CHERN| Spectral distance [hartree]: ', gap
         WRITE (output_unit, '(T2,A,I0,A,I0)') 'REAL_SPACE_CHERN| Occupied states: ', nocc, ' / ', n
         IF (do_soc) THEN
            WRITE (output_unit, '(T2,A)') 'REAL_SPACE_CHERN| Post-SCF GTH SOC, both spinor components.'
         ELSE
            WRITE (output_unit, '(T2,A,I0)') 'REAL_SPACE_CHERN| Scalar spin channel, no degeneracy factor: ', spin
         END IF
         WRITE (output_unit, '(T2,A,ES23.14)') 'REAL_SPACE_CHERN| Window occupied weight: ', weight
         WRITE (output_unit, '(T2,A,ES23.14)') 'REAL_SPACE_CHERN| Marker: ', marker
      END IF
      CALL cite_reference(BiancoResta2011)
      DO i = 1, 4
         CALL localizer_dbcsr_release(operators(i))
      END DO
      CALL localizer_dbcsr_release(metric)
      CALL timestop(handle)
   END SUBROUTINE qs_scf_post_real_space_chern
END MODULE qs_real_space_chern
