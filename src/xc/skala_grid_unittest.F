!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM skala_grid_unittest
   USE OMP_LIB,                         ONLY: omp_get_max_threads
   USE cell_methods,                    ONLY: cell_create
   USE cell_types,                      ONLY: cell_release,&
                                              cell_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_cp2k_xc,                   ONLY: create_xc_section
   USE input_section_types,             ONLY: section_release,&
                                              section_type,&
                                              section_vals_create,&
                                              section_vals_release,&
                                              section_vals_type,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_comm_world
   USE particle_types,                  ONLY: particle_type
   USE skala_gpw_features,              ONLY: skala_grid_layout_type
   USE skala_gpw_functional,            ONLY: skala_grid_eval,&
                                              skala_point_eval

   IMPLICIT NONE

   TYPE(mp_comm_type), TARGET :: comm
   TYPE(skala_grid_layout_type) :: grid
   TYPE(cell_type), POINTER :: cell => NULL()
   TYPE(particle_type), POINTER :: atoms(:) => NULL()
   TYPE(section_type), POINTER :: definition => NULL()
   TYPE(section_vals_type), POINTER :: xc => NULL()
   CHARACTER(len=default_path_length) :: model
   CHARACTER(len=32) :: model_mode, device_arg
   LOGICAL :: cuda_points, spinor_test, spinor_precision_test, point_fp64
   INTEGER :: ns, partition, i, s, first_z, last_z, nr, component, step, mode, caller_threads
   REAL(dp) :: hmat(3, 3), energy, reference, ep, em, derivative, h, error, maximum_error, &
               step_size, response_step_size, routing_tolerance, derivative_tolerance
   REAL(dp), ALLOCATABLE :: rho(:, :), grad(:, :, :), tau(:, :), drho(:, :), dgrad(:, :, :), dtau(:, :), &
                            vrho(:, :), vgrad(:, :, :), vtau(:, :), ref_vrho(:, :), ref_vgrad(:, :, :), ref_vtau(:, :)

   CALL get_command_argument(1, model)
   IF (LEN_TRIM(model) == 0) ERROR STOP "Pass a model exported by create_skala_grid_test_model.py"
   IF (LEN_TRIM(model) > default_string_length) THEN
      ERROR STOP "Use MODEL SKALA and set SKALA_MODEL for a long model path"
   END IF
   CALL get_command_argument(2, model_mode)
   spinor_precision_test = model_mode == "spinor-precision"
   spinor_test = model_mode == "spinor" .OR. model_mode == "spinor-float32" .OR. spinor_precision_test
   IF (model_mode == "spinor") model_mode = "geometry"
   IF (model_mode == "spinor-float32" .OR. spinor_precision_test) model_mode = "skala-float32"
   point_fp64 = TRIM(model_mode) == "skala-float32"
   CALL get_command_argument(3, device_arg)
   cuda_points = TRIM(device_arg) == "gpu"
   IF (LEN_TRIM(device_arg) > 0 .AND. .NOT. cuda_points) ERROR STOP "Expected gpu or no device argument"
   step_size = 1.0e-3_dp
   response_step_size = step_size
   routing_tolerance = 1.0e-11_dp
   derivative_tolerance = 1.0e-8_dp
   SELECT CASE (TRIM(model_mode))
   CASE ("", "geometry")
      CONTINUE
   CASE ("skala-float32")
      step_size = 1.0e-2_dp
      ! Differencing Float32 gradients needs a longer stencil than energy differences.
      response_step_size = 4.0e-2_dp
      routing_tolerance = 1.0e-5_dp
      derivative_tolerance = 1.0e-4_dp
   CASE DEFAULT
      ERROR STOP "Unknown model mode (expected geometry, skala-float32 or no mode for the analytic fixture)"
   END SELECT
   CALL cp2k_init()
   caller_threads = omp_get_max_threads()
   comm = mp_comm_world
   CALL comm%init()
   CALL create_xc_section(definition)
   CALL section_vals_create(xc, definition)
   CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%MODEL", c_val=TRIM(model))
   CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID", l_val=.TRUE.)
   hmat = RESHAPE([6.0_dp, 0.0_dp, 0.0_dp, 0.4_dp, 5.0_dp, 0.0_dp, 0.2_dp, 0.3_dp, 7.0_dp], [3, 3])
   CALL cell_create(cell, hmat=hmat, periodic=[1, 1, 1])
   ALLOCATE (atoms(2))
   atoms(1)%r = [0.8_dp, 0.6_dp, 1.1_dp]
   atoms(2)%r = [3.5_dp, 3.0_dp, 4.2_dp]
   grid%group => comm
   grid%npts = [4, 4, 3]
   IF (model_mode == "skala-float32") grid%npts = [8, 8, 7]
   grid%bounds(1, :) = 0
   grid%bounds(2, :) = grid%npts - 1
   first_z = grid%npts(3)*comm%mepos/comm%num_pe
   last_z = grid%npts(3)*(comm%mepos + 1)/comm%num_pe - 1
   grid%bounds_local = grid%bounds
   grid%bounds_local(:, 3) = [first_z, last_z]
   nr = PRODUCT(grid%npts(1:2))*(last_z - first_z + 1)
   grid%ngpts_local = nr
   grid%dvol = cell%deth/REAL(PRODUCT(grid%npts), dp)
   DO i = 1, 3
      grid%dh(:, i) = hmat(:, i)/REAL(grid%npts(i), dp)
   END DO
   maximum_error = 0.0_dp
   DO ns = 1, 2
      IF (cuda_points .OR. spinor_test) EXIT
      ALLOCATE (rho(nr, ns), grad(nr, 3, ns), tau(nr, ns), drho(nr, ns), dgrad(nr, 3, ns), dtau(nr, ns), &
                vrho(nr, ns), vgrad(nr, 3, ns), vtau(nr, ns), &
                ref_vrho(nr, ns), ref_vgrad(nr, 3, ns), ref_vtau(nr, ns))
      DO s = 1, ns
         DO i = 1, nr
            rho(i, s) = 0.1_dp + 0.03_dp*SIN(REAL(i + 16*first_z + s, dp))
            tau(i, s) = 0.04_dp + 0.01_dp*COS(REAL(i + 16*first_z + 2*s, dp))
            grad(i, :, s) = 0.01_dp*SIN(REAL(i + 16*first_z + s + [1, 2, 3], dp))
            drho(i, s) = 0.01_dp*COS(REAL(i + 16*first_z + s, dp))
            dtau(i, s) = 0.02_dp*SIN(REAL(i + 16*first_z + 2*s, dp))
            dgrad(i, :, s) = 0.01_dp*COS(REAL(i + 16*first_z + s + [1, 2, 3], dp))
         END DO
      END DO
      IF (model_mode == "skala-float32") CALL atomic_fields()
      DO partition = 1, 2
         CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID_ATOM_PARTITION", i_val=partition)
         DO mode = 0, 1
            CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID_ATOM_CHUNKS", l_val=(mode == 1))
            CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID_ATOM_CHUNK_MAX_ROWS", i_val=1)
            CALL evaluate(rho, grad, tau, energy)
            IF (model_mode == "skala-float32" .AND. energy >= -1.0e-6_dp) THEN
               ERROR STOP "The Skala field fixture must have a nontrivial XC energy"
            END IF
            IF (mode == 0) THEN
               reference = energy
               ref_vrho(:, :) = vrho
               ref_vgrad(:, :, :) = vgrad
               ref_vtau(:, :) = vtau
            ELSE
               error = ABS(energy - reference)
               IF (nr > 0) error = MAX(error, MAXVAL(ABS(vrho - ref_vrho)), &
                                       MAXVAL(ABS(vgrad - ref_vgrad)), MAXVAL(ABS(vtau - ref_vtau)))
               CALL comm%max(error)
               IF (error > routing_tolerance) ERROR STOP "Atom routing changes the energy or primitive derivatives"
            END IF
            DO component = 1, 3
               SELECT CASE (component)
               CASE (1)
                  derivative = SUM(ref_vrho*drho)*grid%dvol
               CASE (2)
                  derivative = SUM(ref_vgrad*dgrad)*grid%dvol
               CASE (3)
                  derivative = SUM(ref_vtau*dtau)*grid%dvol
               END SELECT
               CALL comm%sum(derivative)
               DO step = 1, 2
                  h = step_size/REAL(step, dp)
                  SELECT CASE (component)
                  CASE (1)
                     CALL evaluate(rho + h*drho, grad, tau, ep)
                     CALL evaluate(rho - h*drho, grad, tau, em)
                  CASE (2)
                     CALL evaluate(rho, grad + h*dgrad, tau, ep)
                     CALL evaluate(rho, grad - h*dgrad, tau, em)
                  CASE (3)
                     CALL evaluate(rho, grad, tau + h*dtau, ep)
                     CALL evaluate(rho, grad, tau - h*dtau, em)
                  END SELECT
                  error = ABS((ep - em)/(2.0_dp*h) - derivative)
                  maximum_error = MAX(maximum_error, error)
                  IF (error > derivative_tolerance) THEN
                     WRITE (*, '(A,4I3,3ES18.8)') "FD failure (spin,partition,mode,component,step,AD,FD): ", &
                        ns, partition, mode, component, h, derivative, (ep - em)/(2.0_dp*h)
                     ERROR STOP "Primitive-field derivative fails finite differences"
                  END IF
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE (rho, grad, tau, drho, dgrad, dtau, vrho, vgrad, vtau, ref_vrho, ref_vgrad, ref_vtau)
   END DO
   CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID_USE_CUDA", l_val=cuda_points)
   IF (spinor_test) THEN
      CALL noncollinear_points()
   ELSE
      CALL arbitrary_points()
   END IF
   CALL comm%max(maximum_error)
   IF (omp_get_max_threads() /= caller_threads) ERROR STOP "Skala evaluators changed caller OpenMP thread count"
   IF (comm%mepos == 0) WRITE (*, '(A,ES12.4)') "Maximum field/geometry derivative error: ", maximum_error
   DEALLOCATE (atoms)
   CALL cell_release(cell)
   CALL section_vals_release(xc)
   CALL section_release(definition)
   CALL cp2k_finalize()
CONTAINS
! **************************************************************************************************
!> \brief Check spin rotations, collinear limits and joint Hessian actions.
! **************************************************************************************************
   SUBROUTINE noncollinear_points()
      INTEGER, PARAMETER                                 :: block_size = 16

      INTEGER                                            :: a, angular, axis, b, first, ih, j, na, &
                                                            nrow, row, source
      INTEGER, ALLOCATABLE                               :: sizes(:)
      REAL(dp) :: ad, angle, charge, displacement(3), e, err, expected, local_energy, mean(4), &
         mean_point(3), moment(4), normal(3), other, phi, pi, radius, rotation(3, 3), step_h, tol, &
         total_descriptor, total_weight, z
      REAL(dp), ALLOCATABLE :: centres(:, :), col_g(:, :, :), col_r(:, :), col_t(:, :), &
         col_vg(:, :, :), col_vr(:, :), col_vt(:, :), dc(:, :), fc(:, :), fd(:, :), &
         field_grad(:, :, :), hc(:, :), hv(:, :), minus_c(:, :), minus_v(:, :), plus_c(:, :), &
         plus_v(:, :), reference_c(:, :), reference_v(:, :), rotated(:, :), tangent(:, :), &
         v(:, :), vc(:, :), y(:, :)

      na = 0
      DO a = 1, 2
         IF (MOD(a - 1, comm%num_pe) == comm%mepos) na = na + 1
      END DO
      nrow = block_size*na
      ! Packed columns are rho(4), grad(3,4), tau(4), coordinates(3), energy/descriptor weights.
      ALLOCATE (sizes(na), y(nrow, 25), tangent(nrow, 25), centres(3, na), dc(3, na), &
                v(nrow, 25), vc(3, na), hv(nrow, 25), hc(3, na), &
                plus_v(nrow, 25), plus_c(3, na), minus_v(nrow, 25), minus_c(3, na), &
                fd(nrow, 25), fc(3, na), rotated(nrow, 25), reference_v(nrow, 25), reference_c(3, na))
      sizes = block_size
      pi = ACOS(-1.0_dp)
      y = 0.0_dp
      b = 0
      row = 0
      DO a = 1, 2
         IF (MOD(a - 1, comm%num_pe) /= comm%mepos) CYCLE
         b = b + 1
         centres(:, b) = atoms(a)%r
         dc(:, b) = 0.02_dp*COS(REAL(a + [1, 2, 3], dp))
         DO j = 1, block_size
            radius = 0.4_dp + 0.5_dp*REAL((j - 1)/8, dp)
            z = 1.0_dp - 2.0_dp*(REAL(MOD(j - 1, 8), dp) + 0.5_dp)/8.0_dp
            phi = REAL(j, dp)*pi*(3.0_dp - SQRT(5.0_dp))
            row = row + 1
            y(row, 21:23) = centres(:, b) + radius*[SQRT(1.0_dp - z*z)*COS(phi), &
                                                    SQRT(1.0_dp - z*z)*SIN(phi), z]
            y(row, 25) = 4.0_dp*pi*radius**2*0.5_dp/8.0_dp
            y(row, 24) = y(row, 25)*(0.6_dp + 0.1_dp*COS(phi))
            IF (j == 1) y(row, 24) = 0.0_dp
            DO source = 1, 2
               displacement = y(row, 21:23) - atoms(source)%r
               charge = EXP(-2.0_dp*SQRT(SUM(displacement**2)))/pi
               moment = [1.0_dp, 0.24_dp, 0.08_dp, 0.0_dp]
               IF (source == 2) moment = [1.0_dp, -0.12_dp, 0.0_dp, 0.20_dp]
               y(row, 1:4) = y(row, 1:4) + charge*moment
               y(row, 17:20) = y(row, 17:20) + 0.5_dp*charge*moment
               DO axis = 1, 3
                  y(row, axis + [4, 7, 10, 13]) = y(row, axis + [4, 7, 10, 13]) - &
                                                  2.0_dp*charge*moment*displacement(axis)/SQRT(SUM(displacement**2))
               END DO
            END DO
            tangent(row, 1:4) = 0.02_dp*y(row, 1)*COS(REAL(j + [1, 2, 3, 4], dp))
            tangent(row, 5:16) = 0.02_dp*y(row, 1)*SIN(REAL(j + [(axis, axis=1, 12)], dp))
            tangent(row, 17:20) = 0.02_dp*y(row, 17)*SIN(REAL(j + [1, 2, 3, 4], dp))
            tangent(row, 21:23) = 0.01_dp*COS(REAL(j + [1, 2, 3], dp))
            tangent(row, 24:25) = 0.02_dp*y(row, 25)*[SIN(REAL(j, dp)), COS(REAL(j, dp))]
         END DO
      END DO
      angular = 26
      tol = 1.0e-8_dp
      IF (model_mode == "skala-float32") THEN
         angular = MERGE(26, 194, spinor_precision_test)
         tol = 5.0e-5_dp
      END IF
      CALL spinor_evaluate(y, centres, sizes, angular, e, v, vc)
      IF (model_mode /= "skala-float32" .OR. spinor_precision_test) THEN
         ! Independent Cartesian polynomial, including a nonlocal spin-vector descriptor.
         expected = 0.0_dp
         reference_v = 0.0_dp
         reference_c = 0.0_dp
         first = 1
         DO b = 1, na
            total_weight = SUM(y(first:first + block_size - 1, 24))
            total_descriptor = SUM(y(first:first + block_size - 1, 25))
            mean = 0.0_dp
            mean_point = 0.0_dp
            DO row = first, first + block_size - 1
               mean = mean + y(row, 25)*y(row, 1:4)/total_descriptor
               mean_point = mean_point + y(row, 25)*(y(row, 21:23) - centres(:, b))/total_descriptor
            END DO
            DO row = first, first + block_size - 1
               displacement = y(row, 21:23) - centres(:, b)
               local_energy = 0.5_dp*(0.7_dp*SUM(y(row, 1:4)**2) + 0.3_dp*SUM(y(row, 17:20)**2) + &
                                      0.2_dp*SUM(y(row, 1:4)*y(row, 17:20))) + 0.05_dp*SUM(y(row, 5:16)**2) + &
                              0.05_dp*mean(1)**2 + 0.02_dp*SUM(mean(2:4)**2) + &
                              0.04_dp*SUM(displacement**2) + (0.03_dp + 0.02_dp*mean(1))*SUM(mean_point**2)
               expected = expected + y(row, 24)*local_energy
               reference_v(row, 1:4) = y(row, 24)*(0.7_dp*y(row, 1:4) + 0.1_dp*y(row, 17:20)) + &
                                       total_weight*y(row, 25)/total_descriptor* &
                                       [0.1_dp*mean(1) + 0.02_dp*SUM(mean_point**2), 0.04_dp*mean(2:4)]
               reference_v(row, 5:16) = 0.1_dp*y(row, 24)*y(row, 5:16)
               reference_v(row, 17:20) = y(row, 24)*(0.3_dp*y(row, 17:20) + 0.1_dp*y(row, 1:4))
               reference_v(row, 21:23) = 0.08_dp*y(row, 24)*displacement + &
                                         (0.06_dp + 0.04_dp*mean(1))*total_weight*y(row, 25)*mean_point/total_descriptor
               reference_c(:, b) = reference_c(:, b) - reference_v(row, 21:23)
               reference_v(row, 24) = local_energy
               reference_v(row, 25) = total_weight/total_descriptor*( &
                                      (0.1_dp*mean(1) + 0.02_dp*SUM(mean_point**2))*(y(row, 1) - mean(1)) + &
                                      0.04_dp*SUM(mean(2:4)*(y(row, 2:4) - mean(2:4))) + &
                                      (0.06_dp + 0.04_dp*mean(1))*SUM(mean_point*(displacement - mean_point)))
            END DO
            first = first + block_size
         END DO
         CALL comm%sum(expected)
         err = MAX(ABS(e - expected), MAXVAL(ABS(v - reference_v)), MAXVAL(ABS(vc - reference_c)))
         CALL comm%max(err)
         IF (err > 1.0e-11_dp) ERROR STOP "Spin-projected analytic energy or covectors are incorrect"
         ALLOCATE (field_grad(nrow, 3, 4))
         CALL skala_point_eval(comm, TRANSPOSE(y(:, 21:23)), y(:, 24), y(:, 25), sizes, centres, &
                               y(:, 1:4), RESHAPE(y(:, 5:16), SHAPE(field_grad)), y(:, 17:20), xc, &
                               other, plus_v(:, 1:4), field_grad, plus_v(:, 17:20), spin_angular_points=angular)
         plus_v(:, 5:16) = RESHAPE(field_grad, [nrow, 12])
         err = MAX(ABS(e - other), MAXVAL(ABS(plus_v(:, 1:20) - v(:, 1:20))))
         CALL comm%max(err)
         IF (err > 1.0e-11_dp) ERROR STOP "Omitting geometry outputs changes spin-projected fields"
      END IF
      normal = [1.0_dp, 2.0_dp, 3.0_dp]/SQRT(14.0_dp)
      angle = 0.73_dp
      DO axis = 1, 3
         rotation(:, axis) = (1.0_dp - COS(angle))*normal*normal(axis)
         rotation(axis, axis) = rotation(axis, axis) + COS(angle)
      END DO
      rotation = rotation + SIN(angle)*RESHAPE([0.0_dp, normal(3), -normal(2), &
                                                -normal(3), 0.0_dp, normal(1), normal(2), -normal(1), 0.0_dp], [3, 3])
      rotated(:, :) = y
      reference_v(:, :) = v
      rotated(:, 2:4) = MATMUL(y(:, 2:4), TRANSPOSE(rotation))
      reference_v(:, 2:4) = MATMUL(v(:, 2:4), TRANSPOSE(rotation))
      rotated(:, 18:20) = MATMUL(y(:, 18:20), TRANSPOSE(rotation))
      reference_v(:, 18:20) = MATMUL(v(:, 18:20), TRANSPOSE(rotation))
      DO axis = 1, 3
         rotated(:, axis + [7, 10, 13]) = MATMUL(y(:, axis + [7, 10, 13]), TRANSPOSE(rotation))
         reference_v(:, axis + [7, 10, 13]) = MATMUL(v(:, axis + [7, 10, 13]), TRANSPOSE(rotation))
      END DO
      CALL spinor_evaluate(rotated, centres, sizes, angular, other, plus_v, plus_c)
      err = MAX(ABS(e - other), MAXVAL(ABS(plus_v - reference_v)), MAXVAL(ABS(plus_c - vc)))
      CALL comm%max(err)
      IF (comm%mepos == 0) WRITE (*, '(A,ES12.4)') "Spin-rotation energy/covector error: ", err
      IF (err > tol) ERROR STOP "Spin-projected energy or covectors violate rotation invariance"
      ad = SUM(v*tangent) + SUM(vc*dc)
      CALL comm%sum(ad)
      CALL spinor_evaluate(y, centres, sizes, angular, other, hv, hc, tangent, dc)
      fd = 0.0_dp
      fc = 0.0_dp
      DO ih = 1, 2
         step_h = response_step_size/REAL(ih, dp)
         CALL spinor_evaluate(y + step_h*tangent, centres + step_h*dc, sizes, angular, &
                              ep, plus_v, plus_c)
         CALL spinor_evaluate(y - step_h*tangent, centres - step_h*dc, sizes, angular, &
                              em, minus_v, minus_c)
         IF (ih == 1) THEN
            fd(:, :) = -(plus_v - minus_v)/(6.0_dp*step_h)
            fc(:, :) = -(plus_c - minus_c)/(6.0_dp*step_h)
         ELSE
            fd(:, :) = fd + 2.0_dp*(plus_v - minus_v)/(3.0_dp*step_h)
            fc(:, :) = fc + 2.0_dp*(plus_c - minus_c)/(3.0_dp*step_h)
            err = ABS((ep - em)/(2.0_dp*step_h) - ad)
            IF (err > tol) ERROR STOP "Spin-projected joint covectors fail energy finite differences"
         END IF
      END DO
      err = MAX(MAXVAL(ABS(fd - hv)), MAXVAL(ABS(fc - hc)))
      CALL comm%max(err)
      maximum_error = MAX(maximum_error, err)
      IF (comm%mepos == 0) WRITE (*, '(A,ES12.4)') "Spin-projected joint Hessian error: ", err
      IF (err > tol) ERROR STOP "Spin-projected joint Hessian fails Richardson finite differences"
      IF (cuda_points) THEN
         CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID_USE_CUDA", l_val=.FALSE.)
         CALL spinor_evaluate(y, centres, sizes, angular, other, plus_v, plus_c)
         err = MAX(ABS(e - other), MAXVAL(ABS(v - plus_v)), MAXVAL(ABS(vc - plus_c)))
         CALL spinor_evaluate(y, centres, sizes, angular, other, plus_v, plus_c, tangent, dc)
         err = MAX(err, MAXVAL(ABS(hv - plus_v)), MAXVAL(ABS(hc - plus_c)))
         CALL comm%max(err)
         CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID_USE_CUDA", l_val=.TRUE.)
         IF (comm%mepos == 0) WRITE (*, '(A,ES12.4)') "CPU/GPU spin-projected covector/Hessian difference: ", err
         IF (err > tol) ERROR STOP "CPU/GPU spin-projected derivatives disagree"
      END IF
      ! A common collinear axis must recover the existing two-spin functional.
      ALLOCATE (col_r(nrow, 2), col_g(nrow, 3, 2), col_t(nrow, 2), &
                col_vr(nrow, 2), col_vg(nrow, 3, 2), col_vt(nrow, 2))
      DO ih = 1, 2
         rotated(:, :) = y
         rotated(:, 3:4) = 0.0_dp
         rotated(:, 11:16) = 0.0_dp
         rotated(:, 19:20) = 0.0_dp
         IF (ih == 2) THEN
            rotated(:, 2) = 0.0_dp
            rotated(:, 8:10) = 0.0_dp
            rotated(:, 18) = 0.0_dp
         END IF
         DO j = 1, 2
            col_r(:, j) = 0.5_dp*(rotated(:, 1) + (3 - 2*j)*rotated(:, 2))
            col_g(:, :, j) = 0.5_dp*(rotated(:, 5:7) + (3 - 2*j)*rotated(:, 8:10))
            col_t(:, j) = 0.5_dp*(rotated(:, 17) + (3 - 2*j)*rotated(:, 18))
         END DO
         CALL skala_point_eval(comm, TRANSPOSE(y(:, 21:23)), y(:, 24), y(:, 25), sizes, centres, &
                               col_r, col_g, col_t, xc, expected, col_vr, col_vg, col_vt, &
                               model_fp64=.TRUE.)
         CALL spinor_evaluate(rotated, centres, sizes, angular, other, plus_v, plus_c)
         reference_v = 0.0_dp
         reference_v(:, 1) = 0.5_dp*SUM(col_vr, DIM=2)
         reference_v(:, 2) = 0.5_dp*(col_vr(:, 1) - col_vr(:, 2))
         reference_v(:, 5:7) = 0.5_dp*SUM(col_vg, DIM=3)
         reference_v(:, 8:10) = 0.5_dp*(col_vg(:, :, 1) - col_vg(:, :, 2))
         reference_v(:, 17) = 0.5_dp*SUM(col_vt, DIM=2)
         reference_v(:, 18) = 0.5_dp*(col_vt(:, 1) - col_vt(:, 2))
         err = MAX(ABS(other - expected), MAXVAL(ABS(plus_v(:, 1:20) - reference_v(:, 1:20))))
         CALL comm%max(err)
         IF (err > MERGE(1.0e-11_dp, tol, spinor_precision_test)) THEN
            ERROR STOP "Spin projection fails the collinear or unpolarized limit"
         END IF
      END DO
      IF (model_mode == "skala-float32") THEN
         ! Saturated neural activations must retain finite second/third derivatives.
         BLOCK
            INTEGER, ALLOCATABLE :: saturation_sizes(:)
            REAL(dp), ALLOCATABLE :: fields(:, :), field_direction(:, :), covectors(:, :), centre_covectors(:, :)

            ALLOCATE (saturation_sizes(na), fields(64*na, 25), field_direction(64*na, 25), &
                      covectors(64*na, 25), centre_covectors(3, na))
            saturation_sizes = 64
            fields = 0.0_dp
            moment = [1.0_dp, 0.24_dp, 0.08_dp, 0.20_dp]
            DO b = 1, na
               DO j = 1, 64
                  row = (b - 1)*64 + j
                  fields(row, 1:4) = 0.001_dp*moment
                  fields(row, 17:20) = 0.01_dp*moment
                  fields(row, 21:23) = centres(:, b) + [-1.304_dp + 2.0_dp*MOD(j - 1, 4), &
                                                        -2.439_dp + 2.25_dp*MOD((j - 1)/4, 4), &
                                                        -3.59_dp + 2.5_dp*((j - 1)/16)]
                  fields(row, 24:25) = 11.25_dp
               END DO
            END DO
            field_direction = 0.0_dp
            field_direction(:, 1:20) = 0.01_dp*fields(:, 1:20)
            CALL spinor_evaluate(fields, centres, saturation_sizes, 6, other, covectors, centre_covectors)
            CALL spinor_evaluate(fields, centres, saturation_sizes, 6, other, covectors, centre_covectors, &
                                 field_direction, 0.0_dp*centres)
         END BLOCK
      END IF
      IF (comm%mepos == 0) WRITE (*, '(A)') "Noncollinear point energy, covector and response checks passed."
   END SUBROUTINE noncollinear_points

! **************************************************************************************************
!> \brief Evaluate spin-projected point energy, covectors and optional Hessian actions.
!> \param values Primitive spinor fields and quadrature data.
!> \param centres Cartesian atom centres, in bohr.
!> \param sizes Point count of each atom block.
!> \param angular Positive Lebedev spin-direction point count.
!> \param e XC energy.
!> \param covectors Primitive-field and quadrature covectors.
!> \param centre_covectors Atom-block-centre covectors.
!> \param direction Primitive-field and quadrature direction.
!> \param centre_direction Atom-block-centre direction.
! **************************************************************************************************
   SUBROUTINE spinor_evaluate(values, centres, sizes, angular, e, covectors, centre_covectors, direction, centre_direction)
      REAL(dp), INTENT(IN)                               :: values(:, :), centres(:, :)
      INTEGER, INTENT(IN)                                :: sizes(:), angular
      REAL(dp), INTENT(OUT)                              :: e, covectors(:, :), &
                                                            centre_covectors(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: direction(:, :), centre_direction(:, :)

      REAL(dp), ALLOCATABLE                              :: grad_out(:, :, :), point_out(:, :)

      ALLOCATE (grad_out(SIZE(values, 1), 3, 4), point_out(3, SIZE(values, 1)))
      IF (PRESENT(direction)) THEN
         CALL skala_point_eval(comm, TRANSPOSE(values(:, 21:23)), values(:, 24), values(:, 25), sizes, centres, &
                               values(:, 1:4), RESHAPE(values(:, 5:16), SHAPE(grad_out)), values(:, 17:20), xc, &
                               e, covectors(:, 1:4), grad_out, covectors(:, 17:20), point_out, centre_covectors, &
                               covectors(:, 24), covectors(:, 25), direction(:, 1:4), &
                               RESHAPE(direction(:, 5:16), SHAPE(grad_out)), direction(:, 17:20), &
                               TRANSPOSE(direction(:, 21:23)), centre_direction, direction(:, 24), direction(:, 25), angular)
      ELSE
         CALL skala_point_eval(comm, TRANSPOSE(values(:, 21:23)), values(:, 24), values(:, 25), sizes, centres, &
                               values(:, 1:4), RESHAPE(values(:, 5:16), SHAPE(grad_out)), values(:, 17:20), xc, &
                               e, covectors(:, 1:4), grad_out, covectors(:, 17:20), point_out, centre_covectors, &
                               covectors(:, 24), covectors(:, 25), spin_angular_points=angular)
      END IF
      covectors(:, 5:16) = RESHAPE(grad_out, SHAPE(covectors(:, 5:16)))
      covectors(:, 21:23) = TRANSPOSE(point_out)
      IF (.NOT. ieee_is_finite(e) .OR. .NOT. ALL(ieee_is_finite(covectors)) .OR. &
          .NOT. ALL(ieee_is_finite(centre_covectors))) ERROR STOP "Nonfinite spin-projected energy or derivatives"
   END SUBROUTINE spinor_evaluate

! **************************************************************************************************
!> \brief Check point fields, quadrature derivatives and CPU/GPU agreement.
! **************************************************************************************************
   SUBROUTINE arbitrary_points()
      INTEGER, PARAMETER                                 :: block_sizes(2) = [160, 192], nang = 16, &
                                                            total_rows = SUM(block_sizes)

      INTEGER                                            :: a, axis, b, c, distribution, first, &
                                                            global_row, ia, ih, ir, na_local, &
                                                            nr_local, nspin, row, source, spin, &
                                                            sx, sy, sz
      INTEGER, ALLOCATABLE                               :: ids(:), sizes(:)
      REAL(dp) :: ad(7), centre_covector(3), charge, descriptor_derivative, direction(3), &
         displacement(3), dw, e, energy_density, err, expected, mean, mean_point(3), minus, &
         occupation, phi, pi, plus, point(3), point_covector(3), radius, ref_e, rho_sum, &
         spin_scale, step_h, sum_descriptor, sum_weights, w, weight_covector, z
      REAL(dp), ALLOCATABLE :: centres(:, :), dc(:, :), dd(:), de(:), descriptor_weights(:), &
         dg(:, :, :), dr(:, :), dt(:, :), dx(:, :), energy_weights(:), g(:, :, :), &
         gathered(:, :, :), hc(:, :), hd(:), he(:), hg(:, :, :), hr(:, :), ht(:, :), hx(:, :), &
         pc(:, :), pd(:), pe(:), pg(:, :, :), plus_c(:, :), plus_d(:), plus_e(:), plus_g(:, :, :), &
         plus_r(:, :), plus_t(:, :), plus_x(:, :), points(:, :), pr(:, :), pt(:, :), px(:, :), &
         r(:, :), reference_grad(:, :, :), t(:, :)

      ! Match SIRIUS FP64 point callbacks, including their joint Hessian actions.
      IF (point_fp64) THEN
         step_size = 1.0e-3_dp
         response_step_size = step_size
         routing_tolerance = 1.0e-11_dp
         derivative_tolerance = 1.0e-8_dp
      END IF
      pi = ACOS(-1.0_dp)
      DO nspin = 1, 2
         ALLOCATE (gathered(total_rows, 13, nspin), reference_grad(total_rows, 13, nspin))
         DO distribution = 0, 1
            na_local = 0
            nr_local = 0
            DO a = 1, 2
               IF (distribution == 0 .AND. comm%mepos /= 0) CYCLE
               IF (distribution == 1 .AND. MOD(a - 1, comm%num_pe) /= comm%mepos) CYCLE
               na_local = na_local + 1
               nr_local = nr_local + block_sizes(a)
            END DO
            ALLOCATE (sizes(na_local), centres(3, na_local), points(3, nr_local), ids(nr_local), &
                      energy_weights(nr_local), descriptor_weights(nr_local), &
                      r(nr_local, nspin), g(nr_local, 3, nspin), t(nr_local, nspin), &
                      dr(nr_local, nspin), dg(nr_local, 3, nspin), dt(nr_local, nspin), &
                      pr(nr_local, nspin), pg(nr_local, 3, nspin), pt(nr_local, nspin), &
                      dx(3, nr_local), dc(3, na_local), de(nr_local), dd(nr_local), &
                      hx(3, nr_local), hc(3, na_local), he(nr_local), hd(nr_local), &
                      plus_x(3, nr_local), plus_c(3, na_local), plus_e(nr_local), plus_d(nr_local), &
                      px(3, nr_local), pc(3, na_local), pe(nr_local), pd(nr_local), &
                      hr(nr_local, nspin), hg(nr_local, 3, nspin), ht(nr_local, nspin), &
                      plus_r(nr_local, nspin), plus_g(nr_local, 3, nspin), plus_t(nr_local, nspin))
            b = 0
            row = 0
            DO a = 1, 2
               IF (distribution == 0 .AND. comm%mepos /= 0) CYCLE
               IF (distribution == 1 .AND. MOD(a - 1, comm%num_pe) /= comm%mepos) CYCLE
               b = b + 1
               sizes(b) = block_sizes(a)
               centres(:, b) = atoms(a)%r
               dc(:, b) = 0.03_dp*COS(REAL(a + [1, 2, 3], dp))
               DO ir = 1, block_sizes(a)/nang
                  radius = 0.12_dp + (REAL(ir, dp) - 0.5_dp)*0.24_dp
                  DO ia = 1, nang
                     z = 1.0_dp - 2.0_dp*(REAL(ia, dp) - 0.5_dp)/REAL(nang, dp)
                     phi = REAL(ia, dp)*pi*(3.0_dp - SQRT(5.0_dp))
                     direction = [SQRT(1.0_dp - z*z)*COS(phi), SQRT(1.0_dp - z*z)*SIN(phi), z]
                     row = row + 1
                     global_row = SUM(block_sizes(1:a - 1)) + (ir - 1)*nang + ia
                     ids(row) = global_row
                     point = atoms(a)%r + radius*direction
                     points(:, row) = point
                     descriptor_weights(row) = 4.0_dp*pi*radius**2*0.24_dp/REAL(nang, dp)
                     energy_weights(row) = descriptor_weights(row)*(0.6_dp + 0.1_dp*COS(phi))
                     ! This row still participates in the nonlocal descriptor.
                     IF (ir == 1 .AND. ia == 1) energy_weights(row) = 0.0_dp
                     dx(:, row) = 0.01_dp*SIN(REAL(global_row + [1, 2, 3], dp))
                     de(row) = 0.02_dp*descriptor_weights(row)*COS(REAL(global_row, dp))
                     dd(row) = 0.03_dp*descriptor_weights(row)*SIN(REAL(global_row, dp))
                     r(row, :) = 0.0_dp
                     g(row, :, :) = 0.0_dp
                     t(row, :) = 0.0_dp
                     DO source = 1, 2
                        DO sz = -1, 1
                           DO sy = -1, 1
                              DO sx = -1, 1
                                 displacement = point - atoms(source)%r - MATMUL(hmat, REAL([sx, sy, sz], dp))
                                 charge = EXP(-2.0_dp*SQRT(SUM(displacement**2)))/pi
                                 DO spin = 1, nspin
                                    occupation = 1.0_dp
                                    IF (nspin == 2) occupation = 0.5_dp + (3 - 2*spin)*(0.1_dp + 0.05_dp*source)
                                    r(row, spin) = r(row, spin) + occupation*charge
                                    t(row, spin) = t(row, spin) + 0.5_dp*occupation*charge
                                    g(row, :, spin) = g(row, :, spin) - &
                                                      2.0_dp*occupation*charge*displacement/SQRT(SUM(displacement**2))
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                     DO spin = 1, nspin
                        dr(row, spin) = 0.03_dp*r(row, spin)*COS(REAL(global_row + spin, dp))
                        dt(row, spin) = 0.03_dp*t(row, spin)*SIN(REAL(global_row + 2*spin, dp))
                        DO axis = 1, 3
                           dg(row, axis, spin) = 0.03_dp*r(row, spin)*COS(REAL(global_row + spin + axis, dp))
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            IF (cuda_points) THEN
               CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID_USE_CUDA", l_val=.FALSE.)
               CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                     r, g, t, xc, plus, plus_r, plus_g, plus_t, plus_x, plus_c, plus_e, plus_d, &
                                     model_fp64=point_fp64)
               CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID_USE_CUDA", l_val=.TRUE.)
            END IF
            CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                  r, g, t, xc, e, pr, pg, pt, px, pc, pe, pd, model_fp64=point_fp64)
            IF (cuda_points) THEN
               err = ABS(e - plus)
               IF (nr_local > 0) err = MAX(err, MAXVAL(ABS(pr - plus_r)), MAXVAL(ABS(pg - plus_g)), &
                                           MAXVAL(ABS(pt - plus_t)), MAXVAL(ABS(px - plus_x)), &
                                           MAXVAL(ABS(pc - plus_c)), MAXVAL(ABS(pe - plus_e)), MAXVAL(ABS(pd - plus_d)))
               CALL comm%max(err)
               IF (comm%mepos == 0) WRITE (*, '(A,ES12.4)') "CPU/GPU point energy/covector difference: ", err
               IF (err > routing_tolerance) ERROR STOP "CPU/GPU point energy or covectors disagree"
            END IF
            IF (model_mode == "skala-float32" .AND. e >= -1.0e-6_dp) THEN
               ERROR STOP "The arbitrary-point fixture must have a nontrivial Skala energy"
            END IF
            IF (model_mode /= "skala-float32") THEN
               expected = 0.0_dp
               err = 0.0_dp
               descriptor_derivative = 0.0_dp
               first = 1
               spin_scale = 1.0_dp
               IF (nspin == 1) spin_scale = 0.5_dp
               DO b = 1, na_local
                  mean = 0.0_dp
                  mean_point = 0.0_dp
                  centre_covector = 0.0_dp
                  sum_weights = SUM(energy_weights(first:first + sizes(b) - 1))
                  sum_descriptor = SUM(descriptor_weights(first:first + sizes(b) - 1))
                  DO row = first, first + sizes(b) - 1
                     mean = mean + SUM(r(row, :))*descriptor_weights(row)/sum_descriptor
                     IF (model_mode == "geometry") THEN
                        mean_point = mean_point + (points(:, row) - centres(:, b))*descriptor_weights(row)/sum_descriptor
                     END IF
                  END DO
                  DO row = first, first + sizes(b) - 1
                     w = energy_weights(row)
                     dw = descriptor_weights(row)
                     energy_density = spin_scale*(0.7_dp*SUM(r(row, :)**2) + &
                                                  0.3_dp*SUM(t(row, :)**2) + 0.2_dp*SUM(r(row, :)*t(row, :)) + &
                                                  0.1_dp*SUM(g(row, :, :)**2)) + 0.05_dp*mean**2
                     weight_covector = 0.1_dp*mean*sum_weights*(SUM(r(row, :)) - mean)/sum_descriptor
                     point_covector = 0.0_dp
                     IF (model_mode == "geometry") THEN
                        displacement = points(:, row) - centres(:, b)
                        energy_density = energy_density + 0.04_dp*SUM(displacement**2) + &
                                         (0.03_dp + 0.02_dp*mean)*SUM(mean_point**2)
                        point_covector = 0.08_dp*w*displacement + &
                                         (0.06_dp + 0.04_dp*mean)*sum_weights*dw*mean_point/sum_descriptor
                        weight_covector = weight_covector + &
                                          (0.06_dp + 0.04_dp*mean)*sum_weights* &
                                          SUM(mean_point*(displacement - mean_point))/sum_descriptor + &
                                          0.02_dp*sum_weights*SUM(mean_point**2)*(SUM(r(row, :)) - mean)/sum_descriptor
                     END IF
                     expected = expected + w*energy_density
                     centre_covector = centre_covector - point_covector
                     err = MAX(err, ABS(pe(row) - energy_density), ABS(pd(row) - weight_covector), &
                               MAXVAL(ABS(px(:, row) - point_covector)))
                     rho_sum = 0.1_dp*mean*sum_weights*dw/sum_descriptor
                     IF (model_mode == "geometry") THEN
                        rho_sum = rho_sum + 0.02_dp*SUM(mean_point**2)*sum_weights*dw/sum_descriptor
                     END IF
                     DO spin = 1, nspin
                        err = MAX(err, ABS(pr(row, spin) - &
                                           (w*spin_scale*(1.4_dp*r(row, spin) + 0.2_dp*t(row, spin)) + rho_sum)))
                        err = MAX(err, ABS(pt(row, spin) - &
                                           w*spin_scale*(0.6_dp*t(row, spin) + 0.2_dp*r(row, spin))))
                        err = MAX(err, MAXVAL(ABS(pg(row, :, spin) - w*spin_scale*0.2_dp*g(row, :, spin))))
                     END DO
                  END DO
                  err = MAX(err, MAXVAL(ABS(pc(:, b) - centre_covector)))
                  descriptor_derivative = MAX(descriptor_derivative, MAXVAL(ABS(pr(first, :))))
                  first = first + sizes(b)
               END DO
               CALL comm%sum(expected)
               CALL comm%max(err)
               CALL comm%max(descriptor_derivative)
               IF (MAX(ABS(e - expected), err) > 1.0e-11_dp) THEN
                  ERROR STOP "Arbitrary-point analytic energy or raw covector is incorrect"
               END IF
               IF (descriptor_derivative < 1.0e-10_dp) THEN
                  ERROR STOP "Zero-energy-weight rows lost their descriptor derivative"
               END IF
            END IF
            gathered = 0.0_dp
            DO row = 1, nr_local
               gathered(ids(row), 1, :) = pr(row, :)
               gathered(ids(row), 2:4, :) = pg(row, :, :)
               gathered(ids(row), 5, :) = pt(row, :)
               gathered(ids(row), 6:8, 1) = px(:, row)
               gathered(ids(row), 9, 1) = pe(row)
               gathered(ids(row), 10, 1) = pd(row)
            END DO
            first = 1
            DO b = 1, na_local
               gathered(ids(first), 11:13, 1) = pc(:, b)
               first = first + sizes(b)
            END DO
            CALL comm%sum(gathered)
            IF (distribution == 0) THEN
               ref_e = e
               reference_grad(:, :, :) = gathered
            ELSE
               IF (MAX(ABS(e - ref_e), MAXVAL(ABS(gathered - reference_grad))) > routing_tolerance) THEN
                  ERROR STOP "Redistributing complete atom blocks changed Skala energy or derivatives"
               END IF
            END IF
            point = SUM(px, DIM=2) + SUM(pc, DIM=2)
            CALL comm%sum(point)
            IF (MAXVAL(ABS(point)) > derivative_tolerance) THEN
               ERROR STOP "Point and block-centre derivatives violate translation invariance"
            END IF
            ad = [SUM(pr*dr), SUM(pg*dg), SUM(pt*dt), SUM(px*dx), SUM(pc*dc), SUM(pe*de), SUM(pd*dd)]
            CALL comm%sum(ad)
            CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                  r, g, t, xc, e, hr, hg, ht, hx, hc, he, hd, &
                                  response_density=dr, response_gradient=dg, response_kinetic=dt, model_fp64=point_fp64)
            DO ih = 1, 2
               step_h = response_step_size/REAL(ih, dp)
               CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                     r + step_h*dr, g + step_h*dg, t + step_h*dt, xc, plus, &
                                     plus_r, plus_g, plus_t, plus_x, plus_c, plus_e, plus_d, model_fp64=point_fp64)
               CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                     r - step_h*dr, g - step_h*dg, t - step_h*dt, xc, minus, &
                                     pr, pg, pt, px, pc, pe, pd, model_fp64=point_fp64)
               err = 0.0_dp
               IF (nr_local > 0) err = MAX(MAXVAL(ABS(hr - (plus_r - pr)/(2.0_dp*step_h))), &
                                           MAXVAL(ABS(hg - (plus_g - pg)/(2.0_dp*step_h))), &
                                           MAXVAL(ABS(ht - (plus_t - pt)/(2.0_dp*step_h))), &
                                           MAXVAL(ABS(hx - (plus_x - px)/(2.0_dp*step_h))), &
                                           MAXVAL(ABS(hc - (plus_c - pc)/(2.0_dp*step_h))), &
                                           MAXVAL(ABS(he - (plus_e - pe)/(2.0_dp*step_h))), &
                                           MAXVAL(ABS(hd - (plus_d - pd)/(2.0_dp*step_h))))
               CALL comm%max(err)
               IF (err > derivative_tolerance) THEN
                  WRITE (*, '(A,3I3,ES18.8)') "Point Hessian failure (spin,distribution,step,error): ", &
                     nspin, distribution, ih, err
                  IF (nr_local > 0) THEN
                     WRITE (*, '(A,3ES18.8)') "Density/gradient/kinetic Hessian errors: ", &
                        MAXVAL(ABS(hr - (plus_r - pr)/(2.0_dp*step_h))), &
                        MAXVAL(ABS(hg - (plus_g - pg)/(2.0_dp*step_h))), &
                        MAXVAL(ABS(ht - (plus_t - pt)/(2.0_dp*step_h)))
                     WRITE (*, '(A,3ES18.8)') "Density/gradient/kinetic Hessian magnitudes: ", &
                        MAXVAL(ABS(hr)), MAXVAL(ABS(hg)), MAXVAL(ABS(ht))
                  END IF
                  ERROR STOP "Arbitrary-point Hessian-vector product fails finite differences"
               END IF
               maximum_error = MAX(maximum_error, err)
            END DO
            CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                  r, g, t, xc, e, hr, hg, ht, hx, hc, he, hd, &
                                  response_density=0.0_dp*dr, response_gradient=0.0_dp*dg, response_kinetic=0.0_dp*dt, &
                                  model_fp64=point_fp64)
            err = 0.0_dp
            IF (nr_local > 0) err = MAX(MAXVAL(ABS(hr)), MAXVAL(ABS(hg)), MAXVAL(ABS(ht)), &
                                        MAXVAL(ABS(hx)), MAXVAL(ABS(hc)), MAXVAL(ABS(he)), MAXVAL(ABS(hd)))
            CALL comm%max(err)
            IF (err /= 0.0_dp) ERROR STOP "A zero direction must give zero field and geometry Hessian responses"
            IF (cuda_points) THEN
               CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID_USE_CUDA", l_val=.FALSE.)
               CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                     r, g, t, xc, plus, plus_r, plus_g, plus_t, plus_x, plus_c, plus_e, plus_d, &
                                     dr, dg, dt, dx, dc, de, dd, model_fp64=point_fp64)
               CALL section_vals_val_set(xc, "XC_FUNCTIONAL%GAUXC%NATIVE_GRID_USE_CUDA", l_val=.TRUE.)
            END IF
            CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                  r, g, t, xc, e, hr, hg, ht, hx, hc, he, hd, &
                                  dr, dg, dt, dx, dc, de, dd, model_fp64=point_fp64)
            IF (cuda_points) THEN
               err = 0.0_dp
               IF (nr_local > 0) err = MAX(MAXVAL(ABS(hr - plus_r)), MAXVAL(ABS(hg - plus_g)), &
                                           MAXVAL(ABS(ht - plus_t)), MAXVAL(ABS(hx - plus_x)), &
                                           MAXVAL(ABS(hc - plus_c)), MAXVAL(ABS(he - plus_e)), MAXVAL(ABS(hd - plus_d)))
               CALL comm%max(err)
               IF (comm%mepos == 0) WRITE (*, '(A,ES12.4)') "CPU/GPU joint point Hessian difference: ", err
               IF (err > routing_tolerance) ERROR STOP "CPU/GPU joint point Hessians disagree"
            END IF
            DO ih = 1, 2
               step_h = response_step_size/REAL(ih, dp)
               CALL skala_point_eval(comm, points + step_h*dx, energy_weights + step_h*de, &
                                     descriptor_weights + step_h*dd, sizes, centres + step_h*dc, &
                                     r + step_h*dr, g + step_h*dg, t + step_h*dt, xc, plus, &
                                     plus_r, plus_g, plus_t, plus_x, plus_c, plus_e, plus_d, model_fp64=point_fp64)
               CALL skala_point_eval(comm, points - step_h*dx, energy_weights - step_h*de, &
                                     descriptor_weights - step_h*dd, sizes, centres - step_h*dc, &
                                     r - step_h*dr, g - step_h*dg, t - step_h*dt, xc, minus, &
                                     pr, pg, pt, px, pc, pe, pd, model_fp64=point_fp64)
               err = 0.0_dp
               IF (nr_local > 0) err = MAX(MAXVAL(ABS(hr - (plus_r - pr)/(2.0_dp*step_h))), &
                                           MAXVAL(ABS(hg - (plus_g - pg)/(2.0_dp*step_h))), &
                                           MAXVAL(ABS(ht - (plus_t - pt)/(2.0_dp*step_h))), &
                                           MAXVAL(ABS(hx - (plus_x - px)/(2.0_dp*step_h))), &
                                           MAXVAL(ABS(hc - (plus_c - pc)/(2.0_dp*step_h))), &
                                           MAXVAL(ABS(he - (plus_e - pe)/(2.0_dp*step_h))), &
                                           MAXVAL(ABS(hd - (plus_d - pd)/(2.0_dp*step_h))))
               CALL comm%max(err)
               IF (err > derivative_tolerance) THEN
                  ERROR STOP "Joint field/geometry Hessian-vector product fails finite differences"
               END IF
               maximum_error = MAX(maximum_error, err)
            END DO
            ! Check the two mixed blocks independently, including the scalar-spin chain rule.
            CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                  r, g, t, xc, e, hr, hg, ht, hx, hc, he, hd, &
                                  response_density=dr, response_gradient=dg, response_kinetic=dt, model_fp64=point_fp64)
            expected = SUM(hx*dx) + SUM(hc*dc) + SUM(he*de) + SUM(hd*dd)
            CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                  r, g, t, xc, e, hr, hg, ht, &
                                  response_density=0.0_dp*dr, response_gradient=0.0_dp*dg, response_kinetic=0.0_dp*dt, &
                                  response_points=dx, response_atom_coords=dc, response_energy_weights=de, &
                                  response_descriptor_weights=dd, model_fp64=point_fp64)
            err = SUM(hr*dr) + SUM(hg*dg) + SUM(ht*dt) - expected
            CALL comm%sum(err)
            IF (ABS(err) > derivative_tolerance) ERROR STOP "Mixed field/geometry Hessian is not symmetric"
            DO c = 1, 7
               DO ih = 1, 2
                  step_h = step_size/REAL(ih, dp)
                  SELECT CASE (c)
                  CASE (1)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                           r + step_h*dr, g, t, xc, plus, pr, pg, pt, model_fp64=point_fp64)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                           r - step_h*dr, g, t, xc, minus, pr, pg, pt, model_fp64=point_fp64)
                  CASE (2)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                           r, g + step_h*dg, t, xc, plus, pr, pg, pt, model_fp64=point_fp64)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                           r, g - step_h*dg, t, xc, minus, pr, pg, pt, model_fp64=point_fp64)
                  CASE (3)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                           r, g, t + step_h*dt, xc, plus, pr, pg, pt, model_fp64=point_fp64)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                           r, g, t - step_h*dt, xc, minus, pr, pg, pt, model_fp64=point_fp64)
                  CASE (4)
                     CALL skala_point_eval(comm, points + step_h*dx, energy_weights, descriptor_weights, sizes, centres, &
                                           r, g, t, xc, plus, pr, pg, pt, model_fp64=point_fp64)
                     CALL skala_point_eval(comm, points - step_h*dx, energy_weights, descriptor_weights, sizes, centres, &
                                           r, g, t, xc, minus, pr, pg, pt, model_fp64=point_fp64)
                  CASE (5)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres + step_h*dc, &
                                           r, g, t, xc, plus, pr, pg, pt, model_fp64=point_fp64)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres - step_h*dc, &
                                           r, g, t, xc, minus, pr, pg, pt, model_fp64=point_fp64)
                  CASE (6)
                     CALL skala_point_eval(comm, points, energy_weights + step_h*de, descriptor_weights, sizes, centres, &
                                           r, g, t, xc, plus, pr, pg, pt, model_fp64=point_fp64)
                     CALL skala_point_eval(comm, points, energy_weights - step_h*de, descriptor_weights, sizes, centres, &
                                           r, g, t, xc, minus, pr, pg, pt, model_fp64=point_fp64)
                  CASE (7)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights + step_h*dd, sizes, centres, &
                                           r, g, t, xc, plus, pr, pg, pt, model_fp64=point_fp64)
                     CALL skala_point_eval(comm, points, energy_weights, descriptor_weights - step_h*dd, sizes, centres, &
                                           r, g, t, xc, minus, pr, pg, pt, model_fp64=point_fp64)
                  END SELECT
                  err = ABS((plus - minus)/(2.0_dp*step_h) - ad(c))
                  maximum_error = MAX(maximum_error, err)
                  IF (err > derivative_tolerance) THEN
                     WRITE (*, '(A,4I3,3ES18.8)') "Point FD failure (spin,distribution,component,step,AD,FD): ", &
                        nspin, distribution, c, ih, step_h, ad(c), (plus - minus)/(2.0_dp*step_h)
                     ERROR STOP "Arbitrary-point derivative fails finite differences"
                  END IF
               END DO
            END DO
            energy_weights = 0.0_dp
            CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                  r, g, t, xc, e, pr, pg, pt, px, pc, pe, pd, model_fp64=point_fp64)
            err = ABS(e)
            ! The energy-weight derivative itself need not vanish at zero weights.
            IF (nr_local > 0) err = MAX(err, MAXVAL(ABS(pr)), MAXVAL(ABS(pg)), MAXVAL(ABS(pt)), &
                                        MAXVAL(ABS(px)), MAXVAL(ABS(pc)), MAXVAL(ABS(pd)))
            CALL comm%max(err)
            IF (err > 1.0e-12_dp) ERROR STOP "Zero energy weights must give zero energy, field and geometry covectors"
            CALL skala_point_eval(comm, points, energy_weights, descriptor_weights, sizes, centres, &
                                  r, g, t, xc, e, hr, hg, ht, &
                                  response_density=dr, response_gradient=dg, response_kinetic=dt, model_fp64=point_fp64)
            err = 0.0_dp
            IF (nr_local > 0) err = MAX(MAXVAL(ABS(hr)), MAXVAL(ABS(hg)), MAXVAL(ABS(ht)))
            CALL comm%max(err)
            IF (err > 1.0e-12_dp) ERROR STOP "Zero energy weights must give zero Hessian-vector product"
            DEALLOCATE (sizes, centres, points, ids, energy_weights, descriptor_weights, &
                        r, g, t, dr, dg, dt, pr, pg, pt, dx, dc, de, dd, px, pc, pe, pd, &
                        hr, hg, ht, hx, hc, he, hd, plus_r, plus_g, plus_t, plus_x, plus_c, plus_e, plus_d)
         END DO
         DEALLOCATE (gathered, reference_grad)
      END DO
      IF (comm%mepos == 0) WRITE (*, '(A)') "Arbitrary-point atom-block energy and derivative checks passed."
   END SUBROUTINE arbitrary_points

! **************************************************************************************************
!> \brief Check direct and reconstructed atom-centred primitive fields.
! **************************************************************************************************
   SUBROUTINE atomic_fields()
      INTEGER                                            :: a, ix, iy, iz, row, spin, sx, sy, sz
      REAL(dp)                                           :: charge, displacement(3), occupation, &
                                                            point(3), radius

      ! Occupied hydrogenic orbitals and periodic images give mutually consistent primitives.
      rho = 0.0_dp
      grad = 0.0_dp
      tau = 0.0_dp
      row = 0
      DO iz = first_z, last_z
         DO iy = 0, grid%npts(2) - 1
            DO ix = 0, grid%npts(1) - 1
               row = row + 1
               point = MATMUL(grid%dh, REAL([ix, iy, iz], dp))
               DO a = 1, SIZE(atoms)
                  DO sz = -1, 1
                     DO sy = -1, 1
                        DO sx = -1, 1
                           displacement = point - atoms(a)%r - MATMUL(hmat, REAL([sx, sy, sz], dp))
                           radius = SQRT(SUM(displacement**2))
                           charge = EXP(-2.0_dp*radius)/ACOS(-1.0_dp)
                           DO spin = 1, ns
                              occupation = 1.0_dp
                              IF (ns == 2) occupation = 0.5_dp + (3 - 2*spin)*(0.1_dp + 0.05_dp*a)
                              rho(row, spin) = rho(row, spin) + occupation*charge
                              tau(row, spin) = tau(row, spin) + 0.5_dp*occupation*charge
                              grad(row, :, spin) = grad(row, :, spin) - &
                                                   2.0_dp*occupation*charge*displacement/radius
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
               DO spin = 1, ns
                  drho(row, spin) = 0.01_dp*rho(row, spin)*COS(REAL(row + 64*first_z + spin, dp))
                  dtau(row, spin) = 0.02_dp*tau(row, spin)*SIN(REAL(row + 64*first_z + spin, dp))
                  dgrad(row, :, spin) = 0.01_dp*grad(row, :, spin)* &
                                        COS(REAL(row + 64*first_z + spin + [1, 2, 3], dp))
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE atomic_fields

! **************************************************************************************************
!> \brief Evaluate regular-grid XC energy and covectors.
!> \param r Spin densities.
!> \param g Cartesian spin-density gradients.
!> \param t Spin kinetic-energy densities.
!> \param e XC energy.
! **************************************************************************************************
   SUBROUTINE evaluate(r, g, t, e)
      REAL(dp), INTENT(IN)                               :: r(:, :), g(:, :, :), t(:, :)
      REAL(dp), INTENT(OUT)                              :: e

      CALL skala_grid_eval(grid, r, g, t, xc, atoms, cell, e, vrho, vgrad, vtau)
      IF (.NOT. ieee_is_finite(e) .OR. .NOT. ALL(ieee_is_finite(vrho)) .OR. &
          .NOT. ALL(ieee_is_finite(vgrad)) .OR. .NOT. ALL(ieee_is_finite(vtau))) THEN
         ERROR STOP "Nonfinite model energy or primitive derivative"
      END IF
   END SUBROUTINE evaluate
END PROGRAM skala_grid_unittest
