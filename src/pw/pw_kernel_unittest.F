!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM pw_kernel_unittest
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_comm_world
   USE mt_util,                         ONLY: mt2d_kernel_strain
   USE pw_grid_types,                   ONLY: FULLSPACE,&
                                              HALFSPACE,&
                                              pw_grid_type
   USE pw_grids,                        ONLY: pw_grid_create,&
                                              pw_grid_release
   USE pw_methods,                      ONLY: do_accurate_sum,&
                                              do_standard_sum,&
                                              pw_integral_ab
   USE pw_poisson_methods,              ONLY: pw_poisson_set,&
                                              pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_mt,&
                                              pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r1d_gs_type

   IMPLICIT NONE
   REAL(dp) :: max_error

   CALL cp2k_init()
   max_error = 0.0_dp
   CALL check_pw_integral_and_poisson()
   CALL check_mt2d_kernel()
   PRINT '(A,ES16.7)', 'PW integral and MT2D kernel maximum scaled error: ', max_error
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Check PW multiplicity and volume normalization, and auxiliary-density Poisson routes.
! **************************************************************************************************
   SUBROUTINE check_pw_integral_and_poisson()
      REAL(dp), PARAMETER :: step = 5.e-4_dp, &
         strain_steps(4) = [step, -step, 2.0_dp*step, -2.0_dp*step]

      COMPLEX(dp)                                        :: a_mode, b_mode
      INTEGER                                            :: axis, ig, sign_step, span, sumtype, test
      REAL(dp) :: cell(3, 3), derivative, energy, energy_dv, energy_v, energy_v_dv, expected, &
         expected_real, kernel, reference, shifted_g(3), strain_energy(4), stress(3, 3), &
         stress_dv(3, 3), stress_only(3, 3), stress_v(3, 3), stress_v_dv(3, 3)
      TYPE(mp_comm_type)                                 :: comm
      TYPE(pw_c1d_gs_type)                               :: a, b, gradient(3), gradient_v_dv(3), &
                                                            potential, potential_v_dv
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_poisson_type)                              :: poisson
      TYPE(pw_r1d_gs_type)                               :: a_real, b_real

      CALL comm%set_handle(mp_comm_world%get_handle())
      ! Nonunit volume and an appreciable truncation factor exp(-pi).
      cell = 0.0_dp
      DO axis = 1, 3
         cell(axis, axis) = 7.0_dp**(1.0_dp/3.0_dp)
      END DO
      a_mode = CMPLX(0.4_dp, -0.2_dp, dp)
      b_mode = CMPLX(-0.7_dp, 0.9_dp, dp)
      NULLIFY (grid)
      DO span = FULLSPACE, HALFSPACE, -1
         ! Reciprocal inputs and outputs need no FFT. The production grid still
         ! supplies the distributed G0 ownership and half/full-space mappings.
         CALL pw_grid_create(grid, comm, cell, grid_span=span, npts=[5, 5, 5], odd=.TRUE., fft_usage=.FALSE.)
         CALL a%create(grid)
         CALL b%create(grid)
         CALL a_real%create(grid)
         CALL b_real%create(grid)
         CALL potential%create(grid)
         CALL potential_v_dv%create(grid)
         DO axis = 1, 3
            CALL gradient(axis)%create(grid)
            CALL gradient_v_dv(axis)%create(grid)
         END DO
         ALLOCATE (poisson%pw_pools(1))
         CALL pw_pool_create(poisson%pw_pools(1)%pool, grid)
         poisson%pw_level = 1
         poisson%parameters%solver = pw_poisson_mt
         poisson%parameters%periodic = [1, 1, 0]
         poisson%parameters%mt_alpha = 1.0_dp
         CALL pw_poisson_set(poisson, cell_hmat=cell, use_level=1)
         DO test = 1, 2
            ! Constant G0 only and G0 plus a Hermitian mode pair.
            a%array = 0.0_dp
            b%array = 0.0_dp
            expected = 0.0_dp
            expected_real = 0.0_dp
            reference = 0.0_dp
            IF (grid%have_g0) THEN
               a%array(1) = 2.0_dp
               b%array(1) = 3.0_dp
            END IF
            expected = 6.0_dp
            expected_real = 6.0_dp
            reference = 6.0_dp*mt2d_kernel_reference([0.0_dp, 0.0_dp, 0.0_dp], 3, cell(3, 3))
            IF (test /= 1) THEN
               DO ig = 1, grid%ngpts_cut_local
                  IF (ANY(grid%g_hat(2:3, ig) /= 0)) CYCLE
                  IF (grid%g_hat(1, ig) == 1) THEN
                     a%array(ig) = a_mode
                     b%array(ig) = b_mode
                  ELSE IF (grid%g_hat(1, ig) == -1) THEN
                     a%array(ig) = CONJG(a_mode)
                     b%array(ig) = CONJG(b_mode)
                  END IF
               END DO
               expected = expected + 2.0_dp*REAL(CONJG(a_mode)*b_mode, dp)
               expected_real = expected_real + 2.0_dp*REAL(a_mode, dp)*REAL(b_mode, dp)
               kernel = mt2d_kernel_reference([2.0_dp*pi/cell(1, 1), 0.0_dp, 0.0_dp], 3, cell(3, 3))
               reference = reference + 2.0_dp*REAL(CONJG(a_mode)*b_mode, dp)*kernel
            END IF
            a_real%array = REAL(a%array, dp)
            b_real%array = REAL(b%array, dp)
            ! G0 must have weight one, the mode pair weight two, and both must
            ! receive V=7 for an integral. just_sum retains the unscaled oracle.
            DO sumtype = do_accurate_sum, do_standard_sum
               CALL check_error(pw_integral_ab(a, b, sumtype=sumtype), 7.0_dp*expected, 3.e-14_dp)
               CALL check_error(pw_integral_ab(a, b, sumtype=sumtype, just_sum=.TRUE.), expected, 3.e-14_dp)
               CALL check_error(pw_integral_ab(a_real, b_real, sumtype=sumtype), 7.0_dp*expected_real, 3.e-14_dp)
               CALL check_error(pw_integral_ab(a_real, b_real, sumtype=sumtype, just_sum=.TRUE.), &
                                expected_real, 3.e-14_dp)
               CALL check_error(pw_integral_ab(a, b_real, sumtype=sumtype), 7.0_dp*expected_real, 3.e-14_dp)
               CALL check_error(pw_integral_ab(a_real, b, sumtype=sumtype), 7.0_dp*expected_real, 3.e-14_dp)
            END DO
            CALL pw_poisson_solve(poisson, a, energy, h_stress=stress, aux_density=b)
            CALL pw_poisson_solve(poisson, a, energy_dv, dvhartree=gradient, h_stress=stress_dv, aux_density=b)
            CALL pw_poisson_solve(poisson, a, energy_v, potential, h_stress=stress_v, aux_density=b)
            CALL pw_poisson_solve(poisson, a, energy_v_dv, potential_v_dv, dvhartree=gradient_v_dv, &
                                  h_stress=stress_v_dv, aux_density=b)
            CALL pw_poisson_solve(poisson, a, h_stress=stress_only, aux_density=b)
            CALL check_error(energy, 0.5_dp*7.0_dp*reference, 3.e-13_dp)
            CALL check_error(energy_dv, energy, 3.e-13_dp)
            CALL check_error(energy_v, energy, 3.e-13_dp)
            CALL check_error(energy_v_dv, energy, 3.e-13_dp)
            CALL check_error(MAXVAL(ABS(stress_dv - stress)), 0.0_dp, 3.e-11_dp)
            CALL check_error(MAXVAL(ABS(stress_v - stress)), 0.0_dp, 3.e-11_dp)
            CALL check_error(MAXVAL(ABS(stress_v_dv - stress)), 0.0_dp, 3.e-11_dp)
            CALL check_error(MAXVAL(ABS(stress_only - stress)), 0.0_dp, 3.e-11_dp)
            IF (test == 2) THEN
               ! Fixed Fourier charges Q=V*rho give E(epsilon)=V/2/(1+epsilon)
               ! times the strained kernel sum. This absolute production-virial
               ! oracle includes both G0/volume and the nonzero-mode kernel term.
               DO sign_step = 1, 4
                  shifted_g = [2.0_dp*pi/cell(1, 1)/(1.0_dp + strain_steps(sign_step)), 0.0_dp, 0.0_dp]
                  reference = 6.0_dp*mt2d_kernel_reference([0.0_dp, 0.0_dp, 0.0_dp], 3, cell(3, 3)) &
                              + 2.0_dp*REAL(CONJG(a_mode)*b_mode, dp)*mt2d_kernel_reference(shifted_g, 3, cell(3, 3))
                  strain_energy(sign_step) = 0.5_dp*7.0_dp*reference/(1.0_dp + strain_steps(sign_step))
               END DO
               derivative = (8.0_dp*(strain_energy(1) - strain_energy(2)) &
                             - strain_energy(3) + strain_energy(4))/(12.0_dp*step)
               CALL check_error(stress(1, 1) + derivative, 0.0_dp, 3.e-9_dp)
            END IF
            DO ig = 1, grid%ngpts_cut_local
               kernel = mt2d_kernel_reference(grid%g(:, ig), 3, cell(3, 3))
               CALL check_error(ABS(potential%array(ig) - kernel*a%array(ig)), 0.0_dp, 3.e-13_dp)
               CALL check_error(ABS(potential_v_dv%array(ig) - kernel*a%array(ig)), 0.0_dp, 3.e-13_dp)
               DO axis = 1, 3
                  CALL check_error(ABS(gradient(axis)%array(ig) &
                                       - CMPLX(0.0_dp, grid%g(axis, ig), dp)*kernel*a%array(ig)), 0.0_dp, 3.e-13_dp)
                  CALL check_error(ABS(gradient_v_dv(axis)%array(ig) &
                                       - CMPLX(0.0_dp, grid%g(axis, ig), dp)*kernel*a%array(ig)), 0.0_dp, 3.e-13_dp)
               END DO
            END DO
         END DO
         CALL poisson%release()
         CALL a%release()
         CALL b%release()
         CALL a_real%release()
         CALL b_real%release()
         CALL potential%release()
         CALL potential_v_dv%release()
         DO axis = 1, 3
            CALL gradient(axis)%release()
            CALL gradient_v_dv(axis)%release()
         END DO
         CALL pw_grid_release(grid)
      END DO
   END SUBROUTINE check_pw_integral_and_poisson

! **************************************************************************************************
!> \brief Compare the MT2D kernel and in-plane virial with independent strained-cell energies.
! **************************************************************************************************
   SUBROUTINE check_mt2d_kernel()
      REAL(dp), PARAMETER :: step = 3.e-4_dp, shifts(4) = [step, -step, 2.0_dp*step, -2.0_dp*step]

      INTEGER                                            :: alpha, beta, ip, iz, jp, miller(3), &
                                                            mode, modes(3, 5), sign_step
      REAL(dp)                                           :: derivative, g(3), h_inv(3, 3), kernel, &
                                                            lengths(3), shifted_g(3), &
                                                            shifted_kernel(4), strain(3, 3)

      ! Rotate the Cartesian vacuum axis. A shear within the periodic plane is
      ! permitted; only normal/periodic mixing is excluded by the MT2D kernel.
      DO iz = 1, 3
         ip = MOD(iz, 3) + 1
         jp = MOD(iz + 1, 3) + 1
         lengths = [9.0_dp, 13.0_dp, 17.0_dp]
         h_inv = 0.0_dp
         DO alpha = 1, 3
            h_inv(alpha, alpha) = 1.0_dp/lengths(alpha)
         END DO
         h_inv(ip, jp) = -1.4_dp/(lengths(ip)*lengths(jp))
         modes = 0
         ! G=0, both odd and even normal-only modes, and several modes with
         ! appreciable exponential truncation reject the periodic K**2 formula.
         modes(iz, 2:3) = [1, 2]
         modes(ip, 4) = 1
         modes(:, 5) = [1, 2, 1]
         DO mode = 1, SIZE(modes, 2)
            miller = modes(:, mode)
            g = 2.0_dp*pi*MATMUL(TRANSPOSE(h_inv), REAL(miller, dp))
            CALL mt2d_kernel_strain(g, iz, lengths(iz), kernel, strain)
            CALL check_error(kernel, mt2d_kernel_reference(g, iz, lengths(iz)), 3.e-14_dp)
            CALL check_error(MAXVAL(ABS(strain(iz, :))), 0.0_dp, 0.0_dp)
            CALL check_error(MAXVAL(ABS(strain(:, iz))), 0.0_dp, 0.0_dp)
            IF (mode <= 3) CALL check_error(MAXVAL(ABS(strain)), 0.0_dp, 0.0_dp)
            DO alpha = 1, 3
               IF (alpha == iz) CYCLE
               DO beta = 1, 3
                  IF (beta == iz) CYCLE
                  DO sign_step = 1, 4
                     ! g'=(I+epsilon_ab)^(-T) g, evaluated independently of dK.
                     shifted_g = g
                     IF (alpha == beta) THEN
                        shifted_g(beta) = g(beta)/(1.0_dp + shifts(sign_step))
                     ELSE
                        shifted_g(beta) = g(beta) - shifts(sign_step)*g(alpha)
                     END IF
                     shifted_kernel(sign_step) = mt2d_kernel_reference(shifted_g, iz, lengths(iz))
                  END DO
                  derivative = (8.0_dp*(shifted_kernel(1) - shifted_kernel(2)) &
                                - shifted_kernel(3) + shifted_kernel(4))/(12.0_dp*step)
                  CALL check_error(strain(alpha, beta), derivative, 3.e-9_dp)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE check_mt2d_kernel

! **************************************************************************************************
!> \brief Evaluate the original MT2D influence kernel independently of its strain derivative.
!> \param g ...
!> \param iz ...
!> \param slab_size ...
!> \return ...
! **************************************************************************************************
   FUNCTION mt2d_kernel_reference(g, iz, slab_size) RESULT(kernel)
      REAL(dp), INTENT(IN)                               :: g(3)
      INTEGER, INTENT(IN)                                :: iz
      REAL(dp), INTENT(IN)                               :: slab_size
      REAL(dp)                                           :: kernel

      REAL(dp)                                           :: g2, g_parallel(3), radius

      g2 = SUM(g**2)
      IF (g2 == 0.0_dp) THEN
         kernel = pi*slab_size**2/2.0_dp
      ELSE
         g_parallel = g
         g_parallel(iz) = 0.0_dp
         radius = slab_size/2.0_dp
         kernel = fourpi/g2 - fourpi/g2*COS(g(iz)*radius)*EXP(-SQRT(SUM(g_parallel**2))*radius)
      END IF
   END FUNCTION mt2d_kernel_reference

! **************************************************************************************************
!> \brief Require the scaled error to be finite and within the supplied tolerance.
!> \param value ...
!> \param reference ...
!> \param tolerance ...
! **************************************************************************************************
   SUBROUTINE check_error(value, reference, tolerance)
      REAL(dp), INTENT(IN)                               :: value, reference, tolerance

      REAL(dp)                                           :: error

      error = ABS(value - reference)/MAX(1.0_dp, ABS(reference))
      max_error = MAX(max_error, error)
      IF (.NOT. (error <= tolerance)) THEN
         PRINT *, 'Coulomb check: value, reference, error, tolerance:', value, reference, error, tolerance
         ERROR STOP 'PW integral or MT2D kernel check failed'
      END IF
   END SUBROUTINE check_error

END PROGRAM pw_kernel_unittest
