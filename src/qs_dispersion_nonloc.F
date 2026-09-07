!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of non local dispersion functionals
!>  Some routines adapted from:
!>  Copyright (C) 2001-2009 Quantum ESPRESSO group
!>  Copyright (C) 2009 Brian Kolb, Timo Thonhauser - Wake Forest University
!>  This file is distributed under the terms of the
!>  GNU General Public License. See the file `License'
!>  in the root directory of the present distribution,
!>  or http://www.gnu.org/copyleft/gpl.txt .
!> \author JGH
! **************************************************************************************************
MODULE qs_dispersion_nonloc
   USE bibliography,                    ONLY: Dion2004,&
                                              Romanperez2009,&
                                              Sabatini2013,&
                                              cite_reference
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE input_constants,                 ONLY: vdw_nl_DRSLL,&
                                              vdw_nl_LMKLL,&
                                              vdw_nl_RVV10,&
                                              xc_vdw_fun_nonloc
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE mathconstants,                   ONLY: pi,&
                                              rootpi
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_grid_types,                   ONLY: HALFSPACE,&
                                              pw_grid_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_derive,&
                                              pw_transfer
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   REAL(KIND=dp), PARAMETER :: epsr = 1.e-12_dp

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_nonloc'

   PUBLIC :: qs_dispersion_nonloc_init, calculate_dispersion_nonloc

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param dispersion_env ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE qs_dispersion_nonloc_init(dispersion_env, para_env)
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_dispersion_nonloc_init'

      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: funit, handle, ipair, itable, nqs, &
                                                            nr_points, vdw_type

      CALL timeset(routineN, handle)

      SELECT CASE (dispersion_env%nl_type)
      CASE DEFAULT
         CPABORT("Unknown vdW-DF functional")
      CASE (vdw_nl_DRSLL, vdw_nl_LMKLL)
         CALL cite_reference(Dion2004)
      CASE (vdw_nl_RVV10)
         CALL cite_reference(Sabatini2013)
      END SELECT
      CALL cite_reference(RomanPerez2009)

      vdw_type = dispersion_env%type
      SELECT CASE (vdw_type)
      CASE DEFAULT
         ! do nothing
      CASE (xc_vdw_fun_nonloc)
         ! setup information on non local functionals
         filename = dispersion_env%kernel_file_name
         IF (para_env%is_source()) THEN
            ! Read the kernel information from file "filename"
            CALL open_file(file_name=filename, unit_number=funit, file_form="FORMATTED")
            READ (funit, *) nqs, nr_points
            READ (funit, *) dispersion_env%r_max
         END IF
         CALL para_env%bcast(nqs)
         CALL para_env%bcast(nr_points)
         CALL para_env%bcast(dispersion_env%r_max)
         ALLOCATE (dispersion_env%q_mesh(nqs), dispersion_env%kernel_table(nqs*(nqs + 1)/2, 0:nr_points, 2))
         dispersion_env%nqs = nqs
         dispersion_env%nr_points = nr_points
         IF (para_env%is_source()) THEN
            !! Read in the values of the q points used to generate this kernel
            READ (funit, "(1p, 4e23.14)") dispersion_env%q_mesh
            ! The file stores kernel values followed by second derivatives, both in
            ! lower-triangular pair order. Keep this immutable table in packed form.
            DO itable = 1, 2
               DO ipair = 1, nqs*(nqs + 1)/2
                  READ (funit, "(1p, 4e23.14)") dispersion_env%kernel_table(ipair, 0:nr_points, itable)
               END DO
            END DO
            CALL close_file(unit_number=funit)
         END IF
         CALL para_env%bcast(dispersion_env%q_mesh)
         CALL para_env%bcast(dispersion_env%kernel_table)
         ! 2nd derivates for interpolation
         ALLOCATE (dispersion_env%d2y_dx2(nqs, nqs))
         CALL initialize_spline_interpolation(dispersion_env%q_mesh, dispersion_env%d2y_dx2)
         !
         dispersion_env%q_cut = dispersion_env%q_mesh(nqs)
         dispersion_env%q_min = dispersion_env%q_mesh(1)
         dispersion_env%dk = 2.0_dp*pi/dispersion_env%r_max

      END SELECT

      CALL timestop(handle)

   END SUBROUTINE qs_dispersion_nonloc_init

! **************************************************************************************************
!> \brief Calculates the non-local vdW functional using the method of Soler
!>        For spin polarized cases we use E(a,b) = E(a+b), i.e. total density
!> \param vxc_rho ...
!> \param rho_r ...
!> \param rho_g ...
!> \param edispersion ...
!> \param dispersion_env ...
!> \param energy_only ...
!> \param pw_pool ...
!> \param xc_pw_pool ...
!> \param para_env ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE calculate_dispersion_nonloc(vxc_rho, rho_r, rho_g, edispersion, &
                                          dispersion_env, energy_only, pw_pool, xc_pw_pool, para_env, virial)
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: vxc_rho, rho_r
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      REAL(KIND=dp), INTENT(OUT)                         :: edispersion
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      LOGICAL, INTENT(IN)                                :: energy_only
      TYPE(pw_pool_type), POINTER                        :: pw_pool, xc_pw_pool
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(virial_type), OPTIONAL, POINTER               :: virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_dispersion_nonloc'
      INTEGER, DIMENSION(3, 3), PARAMETER :: nd = RESHAPE([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])

      INTEGER                                            :: handle, handle_fft, i, i_grid, idir, &
                                                            ispin, nl_type, np, nspin, p, q, r, s
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: q_low
      INTEGER, DIMENSION(1:3)                            :: hi, lo, n
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: b_value, beta, Ec_nl, sumnp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dq0_dgradrho, dq0_drho, hpot, q0, rho, &
                                                            theta_scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: drho, spline_coeff, u_contract
      REAL(KIND=dp), CONTIGUOUS, POINTER                 :: tmp_1d(:), vxc_1d(:)
      TYPE(pw_c1d_gs_type)                               :: div_g, rho_tot_g, tmp_g, vxc_g
      TYPE(pw_c1d_gs_type), ALLOCATABLE, DIMENSION(:)    :: thetas_g
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type)                               :: tmp_r, vxc_r

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(rho_r))
      CPASSERT(ASSOCIATED(rho_g))
      CPASSERT(ASSOCIATED(pw_pool))

      IF (PRESENT(virial)) THEN
         use_virial = virial%pv_calculate .AND. (.NOT. virial%pv_numer)
      ELSE
         use_virial = .FALSE.
      END IF
      IF (use_virial) THEN
         CPASSERT(.NOT. energy_only)
      END IF
      IF (.NOT. energy_only) THEN
         CPASSERT(ASSOCIATED(vxc_rho))
      END IF

      nl_type = dispersion_env%nl_type

      b_value = dispersion_env%b_value
      beta = 0.03125_dp*(3.0_dp/(b_value**2.0_dp))**0.75_dp
      nspin = SIZE(rho_r)

      ! temporary arrays for FFT
      CALL pw_pool%create_pw(tmp_g)
      CALL pw_pool%create_pw(tmp_r)

      ! Sum the spin densities on the vdW grid before transforming or differentiating.
      CALL pw_pool%create_pw(rho_tot_g)
      CALL pw_transfer(rho_g(1), rho_tot_g)
      DO ispin = 2, nspin
         CALL pw_transfer(rho_g(ispin), tmp_g)
         CALL pw_axpy(tmp_g, rho_tot_g, 1._dp)
      END DO
      CALL pw_transfer(rho_tot_g, tmp_r)

      np = SIZE(tmp_r%array)
      tmp_1d(1:np) => tmp_r%array
      ALLOCATE (rho(np), drho(np, 3))
      DO i = 1, 3
         lo(i) = LBOUND(tmp_r%array, i)
         hi(i) = UBOUND(tmp_r%array, i)
         n(i) = hi(i) - lo(i) + 1
      END DO
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             SHARED(n, lo, rho, tmp_r) PRIVATE(s) COLLAPSE(3)
      DO r = 0, n(3) - 1
         DO q = 0, n(2) - 1
            DO p = 0, n(1) - 1
               s = r*n(2)*n(1) + q*n(1) + p + 1
               rho(s) = tmp_r%array(p + lo(1), q + lo(2), r + lo(3))
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
      DO idir = 1, 3
         CALL pw_transfer(rho_tot_g, tmp_g)
         CALL pw_derive(tmp_g, nd(:, idir))
         CALL pw_transfer(tmp_g, tmp_r)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             SHARED(idir, n, lo, drho, tmp_r) PRIVATE(s) COLLAPSE(3)
         DO r = 0, n(3) - 1
            DO q = 0, n(2) - 1
               DO p = 0, n(1) - 1
                  s = r*n(2)*n(1) + q*n(1) + p + 1
                  drho(s, idir) = tmp_r%array(p + lo(1), q + lo(2), r + lo(3))
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
      END DO
      CALL pw_pool%give_back_pw(rho_tot_g)

      !! ---------------------------------------------------------------------------------
      !! Find the value of q0 for all assigned grid points.  q is defined in equations
      !! 11 and 12 of DION and q0 is the saturated version of q defined in equation
      !! 5 of SOLER.  This routine also returns the derivatives of the q0s with respect
      !! to the charge-density and the gradient of the charge-density.  These are needed
      !! for the potential calculated below.
      !! ---------------------------------------------------------------------------------

      IF (energy_only) THEN
         ALLOCATE (q0(np))
         SELECT CASE (nl_type)
         CASE DEFAULT
            CPABORT("Unknown vdW-DF functional")
         CASE (vdw_nl_DRSLL, vdw_nl_LMKLL)
            CALL get_q0_on_grid_eo_vdw(rho, drho, q0, dispersion_env)
         CASE (vdw_nl_RVV10)
            CALL get_q0_on_grid_eo_rvv10(rho, drho, q0, dispersion_env)
         END SELECT
      ELSE
         ALLOCATE (q0(np), dq0_drho(np), dq0_dgradrho(np))
         SELECT CASE (nl_type)
         CASE DEFAULT
            CPABORT("Unknown vdW-DF functional")
         CASE (vdw_nl_DRSLL, vdw_nl_LMKLL)
            CALL get_q0_on_grid_vdw(rho, drho, q0, dq0_drho, dq0_dgradrho, dispersion_env)
         CASE (vdw_nl_RVV10)
            CALL get_q0_on_grid_rvv10(rho, drho, q0, dq0_drho, dq0_dgradrho, dispersion_env)
         END SELECT
      END IF

      ! Generate one theta channel directly in the FFT workspace. Only reciprocal-space
      ! channels must coexist for the convolution; no real-space np-by-nqs array is needed.
      ALLOCATE (q_low(np), spline_coeff(np, 4), theta_scale(np))
      CALL prepare_splines(q0, rho, dispersion_env, q_low, spline_coeff, theta_scale)
      ALLOCATE (thetas_g(dispersion_env%nqs))
      CALL timeset("vdW_theta_forward", handle_fft)
      DO i = 1, dispersion_env%nqs
         CALL build_theta(i, q_low, spline_coeff, theta_scale, dispersion_env, tmp_1d)
         CALL pw_pool%create_pw(thetas_g(i))
         CALL pw_transfer(tmp_r, thetas_g(i))
      END DO
      CALL timestop(handle_fft)
      DEALLOCATE (spline_coeff, theta_scale)
      grid => thetas_g(1)%pw_grid
      !! ---------------------------------------------------------------------------------------------
      !! Carry out the integration in equation 8 of SOLER.  This also turns the thetas array into the
      !! precursor to the u_i(k) array which is inverse fourier transformed to get the u_i(r) functions
      !! of SOLER equation 11.  Add the energy we find to the output variable etxc.
      !! --------------------------------------------------------------------------------------------------
      sumnp = np
      CALL para_env%sum(sumnp)
      IF (use_virial) THEN
         ! calculates kernel contribution to stress
         CALL vdW_energy(thetas_g, dispersion_env, Ec_nl, energy_only, virial)
         SELECT CASE (nl_type)
         CASE (vdw_nl_RVV10)
            Ec_nl = 0.5_dp*Ec_nl + beta*SUM(rho(:))*grid%vol/sumnp
         END SELECT
         ! calculates energy contribution to stress
         ! potential contribution to stress is calculated together with other potentials (Hxc)
         DO idir = 1, 3
            virial%pv_xc(idir, idir) = virial%pv_xc(idir, idir) + Ec_nl
         END DO
      ELSE
         CALL vdW_energy(thetas_g, dispersion_env, Ec_nl, energy_only)
         SELECT CASE (nl_type)
         CASE (vdw_nl_RVV10)
            Ec_nl = 0.5_dp*Ec_nl + beta*SUM(rho(:))*grid%vol/sumnp
         END SELECT
      END IF
      CALL para_env%sum(Ec_nl)
      IF (nl_type == vdw_nl_RVV10) Ec_nl = Ec_nl*dispersion_env%scale_rvv10
      edispersion = Ec_nl

      IF (energy_only) THEN
         DEALLOCATE (q0, q_low)
      ELSE
         ! Accumulate the two spline contractions and their endpoint values as each
         ! inverse FFT finishes. Keep the original potential-side knot convention.
         ALLOCATE (u_contract(np, 4), hpot(np))
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(np, q_low, q0, dispersion_env, u_contract) PRIVATE(s)
         DO i_grid = 1, np
            s = q_low(i_grid)
            IF (s > 0 .AND. s < dispersion_env%nqs - 1) THEN
               IF (q0(i_grid) == dispersion_env%q_mesh(s + 1)) q_low(i_grid) = s + 1
            END IF
            u_contract(i_grid, :) = 0.0_dp
         END DO
!$OMP END PARALLEL DO
         CALL timeset("vdW_theta_inverse", handle_fft)
         DO i = 1, dispersion_env%nqs
            CALL pw_transfer(thetas_g(i), tmp_r)
            CALL accumulate_potential(i, q_low, dispersion_env, tmp_1d, u_contract)
         END DO
         CALL timestop(handle_fft)

         ! Write the local potential directly into its PW object.
         CALL pw_pool%create_pw(vxc_r)
         vxc_1d(1:np) => vxc_r%array
         IF (use_virial) THEN
            grid => tmp_g%pw_grid
            CALL get_potential(q0, dq0_drho, dq0_dgradrho, rho, q_low, u_contract, vxc_1d, hpot, &
                               dispersion_env, drho, grid%dvol, virial)
         ELSE
            CALL get_potential(q0, dq0_drho, dq0_dgradrho, rho, q_low, u_contract, vxc_1d, hpot, &
                               dispersion_env)
         END IF
         DEALLOCATE (u_contract, q_low, q0, dq0_drho, dq0_dgradrho)
         SELECT CASE (nl_type)
         CASE (vdw_nl_RVV10)
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(np, vxc_1d, hpot, beta, dispersion_env)
            DO i_grid = 1, np
               vxc_1d(i_grid) = (0.5_dp*vxc_1d(i_grid) + beta)*dispersion_env%scale_rvv10
               hpot(i_grid) = 0.5_dp*dispersion_env%scale_rvv10*hpot(i_grid)
            END DO
!$OMP END PARALLEL DO
         END SELECT
         NULLIFY (vxc_1d)
         ! Sum the derivatives before the inverse FFT. Keep the real-space projection
         ! before transferring the potential to a possibly different XC grid.
         CALL pw_pool%create_pw(div_g)
         DO idir = 1, 3
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             SHARED(n, lo, tmp_r, hpot, drho, idir) PRIVATE(s) COLLAPSE(3)
            DO r = 0, n(3) - 1
               DO q = 0, n(2) - 1
                  DO p = 0, n(1) - 1
                     s = r*n(2)*n(1) + q*n(1) + p + 1
                     tmp_r%array(p + lo(1), q + lo(2), r + lo(3)) = hpot(s)*drho(s, idir)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
            CALL pw_transfer(tmp_r, tmp_g)
            CALL pw_derive(tmp_g, nd(:, idir))
            IF (idir == 1) THEN
               CALL pw_transfer(tmp_g, div_g)
            ELSE
               CALL pw_axpy(tmp_g, div_g, 1._dp)
            END IF
         END DO
         CALL pw_transfer(div_g, tmp_r)
         CALL pw_pool%give_back_pw(div_g)
         CALL pw_axpy(tmp_r, vxc_r, -1._dp)
         CALL pw_transfer(vxc_r, tmp_g)
         CALL pw_pool%give_back_pw(vxc_r)
         CALL xc_pw_pool%create_pw(vxc_r)
         CALL xc_pw_pool%create_pw(vxc_g)
         CALL pw_transfer(tmp_g, vxc_g)
         CALL pw_transfer(vxc_g, vxc_r)
         DO ispin = 1, nspin
            CALL pw_axpy(vxc_r, vxc_rho(ispin), 1._dp)
         END DO
         CALL xc_pw_pool%give_back_pw(vxc_r)
         CALL xc_pw_pool%give_back_pw(vxc_g)
      END IF

      NULLIFY (tmp_1d)

      DO i = 1, dispersion_env%nqs
         CALL pw_pool%give_back_pw(thetas_g(i))
      END DO
      CALL pw_pool%give_back_pw(tmp_r)
      CALL pw_pool%give_back_pw(tmp_g)

      DEALLOCATE (rho, drho, thetas_g)

      CALL timestop(handle)

   END SUBROUTINE calculate_dispersion_nonloc

! **************************************************************************************************
!> \brief This routine carries out the integration of equation 8 of SOLER.  It returns the non-local
!> exchange-correlation energy and the u_alpha(k) arrays used to find the u_alpha(r) arrays via
!> equations 11 and 12 in SOLER.
!> energy contribution to stress is added in qs_force
!> \param thetas_g ...
!> \param dispersion_env ...
!> \param vdW_xc_energy ...
!> \param energy_only ...
!> \param virial ...
!> \par History
!>    OpenMP added: Aug 2016  MTucker
! **************************************************************************************************
   SUBROUTINE vdW_energy(thetas_g, dispersion_env, vdW_xc_energy, energy_only, virial)
      TYPE(pw_c1d_gs_type), DIMENSION(:), INTENT(IN)     :: thetas_g
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      REAL(KIND=dp), INTENT(OUT)                         :: vdW_xc_energy
      LOGICAL, INTENT(IN)                                :: energy_only
      TYPE(virial_type), OPTIONAL, POINTER               :: virial

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'vdW_energy'

      INTEGER                                            :: handle, ig, iq, l, m, nl_type, nqs, &
                                                            q1_i, q2_i
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: g, g2, g2_last, g_multiplier, gm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: theta_im, theta_re, u_im, u_re
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dkernel_of_dk, kernel_of_k
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_thread
      TYPE(pw_grid_type), POINTER                        :: grid

      CALL timeset(routineN, handle)
      nqs = dispersion_env%nqs

      use_virial = PRESENT(virial)
      virial_thread(:, :) = 0.0_dp ! always initialize to avoid floating point exceptions in OMP REDUCTION

      vdW_xc_energy = 0._dp
      grid => thetas_g(1)%pw_grid

      IF (grid%grid_span == HALFSPACE) THEN
         g_multiplier = 2._dp
      ELSE
         g_multiplier = 1._dp
      END IF

      nl_type = dispersion_env%nl_type

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          SHARED(nqs, energy_only, grid, dispersion_env, use_virial, thetas_g, &
!$OMP                 g_multiplier, nl_type) &
!$OMP          PRIVATE(g2_last, kernel_of_k, dkernel_of_dk, theta_re, theta_im, &
!$OMP                  g2, g, iq, q2_i, u_re, u_im, q1_i, gm, l, m) &
!$OMP          REDUCTION(+:vdW_xc_energy, virial_thread)

      g2_last = HUGE(0._dp)

      ALLOCATE (kernel_of_k(nqs, nqs))
      IF (use_virial) ALLOCATE (dkernel_of_dk(nqs, nqs))
      ALLOCATE (theta_re(nqs), theta_im(nqs), u_re(nqs), u_im(nqs))

!$OMP DO
      DO ig = 1, grid%ngpts_cut_local
         g2 = grid%gsq(ig)
         IF (ABS(g2 - g2_last) > 1.e-10) THEN
            g2_last = g2
            g = SQRT(g2)
            IF (use_virial) THEN
               CALL interpolate_kernel(g, kernel_of_k, dispersion_env, dispersion_env%kernel_table, dkernel_of_dk)
            ELSE
               CALL interpolate_kernel(g, kernel_of_k, dispersion_env, dispersion_env%kernel_table)
            END IF
         END IF
         ! Save all inputs before in-place output. Vectorize over output channels;
         ! the input-channel accumulation order remains unchanged for every output.
         DO iq = 1, nqs
            theta_re(iq) = REAL(thetas_g(iq)%array(ig), KIND=dp)
            theta_im(iq) = AIMAG(thetas_g(iq)%array(ig))
         END DO
         u_re(:) = 0.0_dp
         u_im(:) = 0.0_dp
         DO q1_i = 1, nqs
!$OMP SIMD
            DO q2_i = 1, nqs
               u_re(q2_i) = u_re(q2_i) + kernel_of_k(q2_i, q1_i)*theta_re(q1_i)
               u_im(q2_i) = u_im(q2_i) + kernel_of_k(q2_i, q1_i)*theta_im(q1_i)
            END DO
!$OMP END SIMD
         END DO
         DO q2_i = 1, nqs
            IF (ig < grid%first_gne0) THEN
               vdW_xc_energy = vdW_xc_energy + (u_re(q2_i)*theta_re(q2_i) + u_im(q2_i)*theta_im(q2_i))
            ELSE
               vdW_xc_energy = vdW_xc_energy &
                               + g_multiplier*(u_re(q2_i)*theta_re(q2_i) + u_im(q2_i)*theta_im(q2_i))
            END IF
            IF (.NOT. energy_only) thetas_g(q2_i)%array(ig) = CMPLX(u_re(q2_i), u_im(q2_i), KIND=dp)
         END DO

         IF (use_virial .AND. ig >= grid%first_gne0) THEN
            ! Reduce over all channel pairs before assembling the stress tensor.
            gm = 0.0_dp
            DO q2_i = 1, nqs
!$OMP SIMD REDUCTION(+:gm)
               DO q1_i = 1, nqs
                  gm = gm + dkernel_of_dk(q1_i, q2_i) &
                       *(theta_re(q1_i)*theta_re(q2_i) + theta_im(q1_i)*theta_im(q2_i))
               END DO
!$OMP END SIMD
            END DO
            gm = 0.5_dp*g_multiplier*grid%vol*gm
            IF (nl_type == vdw_nl_RVV10) gm = 0.5_dp*gm
            DO l = 1, 3
               DO m = 1, l
                  virial_thread(l, m) = virial_thread(l, m) - gm*(grid%g(l, ig)*grid%g(m, ig))/g
               END DO
            END DO
         END IF
      END DO
!$OMP END DO

      DEALLOCATE (theta_re, theta_im, u_re, u_im, kernel_of_k)
      IF (use_virial) DEALLOCATE (dkernel_of_dk)

!$OMP END PARALLEL

      vdW_xc_energy = vdW_xc_energy*grid%vol*0.5_dp

      IF (use_virial) THEN
         DO l = 1, 3
            DO m = 1, (l - 1)
               virial%pv_xc(l, m) = virial%pv_xc(l, m) + virial_thread(l, m)
               virial%pv_xc(m, l) = virial%pv_xc(l, m)
            END DO
            m = l
            virial%pv_xc(l, m) = virial%pv_xc(l, m) + virial_thread(l, m)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE vdW_energy

! **************************************************************************************************
!> \brief This routine finds the non-local correlation contribution to the potential
!> (i.e. the derivative of the non-local piece of the energy with respect to
!> density) given in SOLER equation 10.  The u_alpha(k) functions were found
!> while calculating the energy. Their spline contractions are accumulated during the inverse FFTs.
!> Most of the required derivatives were calculated in the "get_q0_on_grid"
!> routine, but the derivative of the interpolation polynomials, P_alpha(q),
!> (SOLER equation 3) with respect to q is interpolated here, along with the
!> polynomials themselves.
!> \param q0 ...
!> \param dq0_drho ...
!> \param dq0_dgradrho ...
!> \param total_rho ...
!> \param q_low_grid Lower spline endpoint at each grid point
!> \param u_contract Second-derivative contractions and endpoint values
!> \param potential ...
!> \param h_prefactor ...
!> \param dispersion_env ...
!> \param drho ...
!> \param dvol ...
!> \param virial ...
!> \par History
!> OpenMP added: Aug 2016  MTucker
! **************************************************************************************************
   SUBROUTINE get_potential(q0, dq0_drho, dq0_dgradrho, total_rho, q_low_grid, u_contract, potential, h_prefactor, &
                            dispersion_env, drho, dvol, virial)

      REAL(dp), DIMENSION(:), INTENT(in)                 :: q0, dq0_drho, dq0_dgradrho, total_rho
      INTEGER, DIMENSION(:), INTENT(IN)                  :: q_low_grid
      REAL(dp), DIMENSION(:, :), INTENT(in)              :: u_contract
      REAL(dp), DIMENSION(:), INTENT(out)                :: potential, h_prefactor
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      REAL(dp), DIMENSION(:, :), INTENT(in), OPTIONAL    :: drho
      REAL(dp), INTENT(IN), OPTIONAL                     :: dvol
      TYPE(virial_type), OPTIONAL, POINTER               :: virial

      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_potential'

      INTEGER                                            :: handle, i_grid, l, m, nl_type, nqs, &
                                                            q_hi, q_low
      LOGICAL                                            :: use_virial
      REAL(dp)                                           :: a, b, b_value, c, const, d, dq, dq_6, e, &
                                                            f, prefactor, tmp_1_2, tmp_1_4, &
                                                            tmp_3_4, u_dp_dq0, u_p
      REAL(dp), DIMENSION(3, 3)                          :: virial_thread
      REAL(dp), DIMENSION(:), POINTER                    :: q_mesh

      CALL timeset(routineN, handle)

      use_virial = PRESENT(virial)
      CPASSERT(.NOT. use_virial .OR. PRESENT(drho))
      CPASSERT(.NOT. use_virial .OR. PRESENT(dvol))

      virial_thread(:, :) = 0.0_dp ! always initialize to avoid floating point exceptions in OMP REDUCTION
      b_value = dispersion_env%b_value
      const = 1.0_dp/(3.0_dp*b_value**(3.0_dp/2.0_dp)*pi**(5.0_dp/4.0_dp))

      q_mesh => dispersion_env%q_mesh
      nqs = dispersion_env%nqs
      nl_type = dispersion_env%nl_type

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          SHARED(nqs, u_contract, q_low_grid, q_mesh, q0, nl_type, potential, h_prefactor, &
!$OMP                 dq0_drho, dq0_dgradrho, total_rho, const, use_virial, drho, dvol, virial) &
!$OMP          PRIVATE(q_low, q_hi, dq, dq_6, A, b, c, d, e, f, u_p, u_dp_dq0, &
!$OMP                  prefactor, l, m, tmp_1_2, tmp_1_4, tmp_3_4) &
!$OMP          REDUCTION(+:virial_thread)

!$OMP DO
      DO i_grid = 1, SIZE(q0)
         potential(i_grid) = 0.0_dp
         h_prefactor(i_grid) = 0.0_dp
         IF (nl_type == vdw_nl_RVV10 .AND. total_rho(i_grid) <= epsr) CYCLE
         q_low = q_low_grid(i_grid)
         q_hi = q_low + 1

         dq = q_mesh(q_hi) - q_mesh(q_low)
         dq_6 = dq/6.0_dp

         a = (q_mesh(q_hi) - q0(i_grid))/dq
         b = (q0(i_grid) - q_mesh(q_low))/dq
         c = (a**3 - a)*dq*dq_6
         d = (b**3 - b)*dq*dq_6
         e = (3.0_dp*a**2 - 1.0_dp)*dq_6
         f = (3.0_dp*b**2 - 1.0_dp)*dq_6

         u_p = a*u_contract(i_grid, 3) + b*u_contract(i_grid, 4) &
               + c*u_contract(i_grid, 1) + d*u_contract(i_grid, 2)
         u_dp_dq0 = (u_contract(i_grid, 4) - u_contract(i_grid, 3))/dq &
                    - e*u_contract(i_grid, 1) + f*u_contract(i_grid, 2)

         !! The first term in equation 13 of SOLER
         SELECT CASE (nl_type)
         CASE DEFAULT
            CPABORT("Unknown vdW-DF functional")
         CASE (vdw_nl_DRSLL, vdw_nl_LMKLL)
            potential(i_grid) = u_p + u_dp_dq0*dq0_drho(i_grid)
            prefactor = u_dp_dq0*dq0_dgradrho(i_grid)
         CASE (vdw_nl_RVV10)
            tmp_1_2 = SQRT(total_rho(i_grid))
            tmp_1_4 = SQRT(tmp_1_2)
            tmp_3_4 = tmp_1_4*tmp_1_4*tmp_1_4
            potential(i_grid) = const*0.75_dp/tmp_1_4*u_p + const*tmp_3_4*u_dp_dq0*dq0_drho(i_grid)
            prefactor = const*tmp_3_4*u_dp_dq0*dq0_dgradrho(i_grid)
         END SELECT
         IF (q0(i_grid) /= q_mesh(nqs)) THEN
            h_prefactor(i_grid) = prefactor
         END IF

         ! The saturation guard applies only to h_prefactor, not to the virial.
         IF (use_virial .AND. ABS(prefactor) > 0.0_dp) THEN
            IF (nl_type == vdw_nl_RVV10) prefactor = 0.5_dp*prefactor
            prefactor = prefactor*dvol
            DO l = 1, 3
               DO m = 1, l
                  virial_thread(l, m) = virial_thread(l, m) - prefactor*drho(i_grid, l)*drho(i_grid, m)
               END DO
            END DO
         END IF
      END DO ! i_grid = 1, SIZE(q0)
!$OMP END DO

!$OMP END PARALLEL

      IF (use_virial) THEN
         DO l = 1, 3
            DO m = 1, (l - 1)
               virial%pv_xc(l, m) = virial%pv_xc(l, m) + virial_thread(l, m)
               virial%pv_xc(m, l) = virial%pv_xc(l, m)
            END DO
            m = l
            virial%pv_xc(l, m) = virial%pv_xc(l, m) + virial_thread(l, m)
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE get_potential

! **************************************************************************************************
!> \brief calculates exponent = sum(from i=1 to hi, ((alpha)**i)/i) ) without <<< calling power >>>
!> \param hi = upper index for sum
!> \param alpha ...
!> \param exponent = output value
!> \par History
!>     Created:  MTucker, Aug 2016
! **************************************************************************************************
   ELEMENTAL SUBROUTINE calculate_exponent(hi, alpha, exponent)
      INTEGER, INTENT(in)                                :: hi
      REAL(dp), INTENT(in)                               :: alpha
      REAL(dp), INTENT(out)                              :: exponent

      INTEGER                                            :: i
      REAL(dp)                                           :: multiplier

      multiplier = alpha
      exponent = alpha

      DO i = 2, hi
         multiplier = multiplier*alpha
         exponent = exponent + (multiplier/i)
      END DO
   END SUBROUTINE calculate_exponent

! **************************************************************************************************
!> \brief calculate exponent = sum(from i=1 to hi, ((alpha)**i)/i) )  without calling power
!>        also calculates derivative using similar series
!> \param hi = upper index for sum
!> \param alpha ...
!> \param exponent = output value
!> \param derivative ...
!> \par History
!>     Created:  MTucker, Aug 2016
! **************************************************************************************************
   ELEMENTAL SUBROUTINE calculate_exponent_derivative(hi, alpha, exponent, derivative)
      INTEGER, INTENT(in)                                :: hi
      REAL(dp), INTENT(in)                               :: alpha
      REAL(dp), INTENT(out)                              :: exponent, derivative

      INTEGER                                            :: i
      REAL(dp)                                           :: multiplier

      derivative = 0.0d0
      multiplier = 1.0d0
      exponent = 0.0d0

      DO i = 1, hi
         derivative = derivative + multiplier
         multiplier = multiplier*alpha
         exponent = exponent + (multiplier/i)
      END DO
   END SUBROUTINE calculate_exponent_derivative

   !! This routine first calculates the q value defined in (DION equations 11 and 12), then
   !! saturates it according to (SOLER equation 5).
! **************************************************************************************************
!> \brief This routine first calculates the q value defined in (DION equations 11 and 12), then
!> saturates it according to (SOLER equation 5).
!> \param total_rho ...
!> \param gradient_rho ...
!> \param q0 ...
!> \param dq0_drho ...
!> \param dq0_dgradrho ...
!> \param dispersion_env ...
! **************************************************************************************************
   SUBROUTINE get_q0_on_grid_vdw(total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho, dispersion_env)
      !!
      !! more specifically it calculates the following
      !!
      !!     q0(ir) = q0 as defined above
      !!     dq0_drho(ir) = total_rho * d q0 /d rho
      !!     dq0_dgradrho = total_rho / |gradient_rho| * d q0 / d |gradient_rho|
      !!
      REAL(dp), INTENT(IN)                               :: total_rho(:), gradient_rho(:, :)
      REAL(dp), INTENT(OUT)                              :: q0(:), dq0_drho(:), dq0_dgradrho(:)
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env

      INTEGER, PARAMETER                                 :: m_cut = 12
      REAL(dp), PARAMETER :: LDA_A = 0.031091_dp, LDA_a1 = 0.2137_dp, LDA_b1 = 7.5957_dp, &
         LDA_b2 = 3.5876_dp, LDA_b3 = 1.6382_dp, LDA_b4 = 0.49294_dp

      INTEGER                                            :: i_grid
      REAL(dp)                                           :: dq0_dq, exponent, gradient_correction, &
                                                            kF, LDA_1, LDA_2, q, q__q_cut, q_cut, &
                                                            q_min, r_s, sqrt_r_s, Z_ab

      q_cut = dispersion_env%q_cut
      q_min = dispersion_env%q_min
      SELECT CASE (dispersion_env%nl_type)
      CASE DEFAULT
         CPABORT("Unknown vdW-DF functional")
      CASE (vdw_nl_DRSLL)
         Z_ab = -0.8491_dp
      CASE (vdw_nl_LMKLL)
         Z_ab = -1.887_dp
      END SELECT

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             SHARED(total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho, q_cut, q_min, Z_ab) &
!$OMP             PRIVATE(dq0_dq, exponent, gradient_correction, kF, LDA_1, LDA_2, q, q__q_cut, r_s, sqrt_r_s) &
!$OMP             SCHEDULE(STATIC)
      DO i_grid = 1, SIZE(total_rho)
         q0(i_grid) = q_cut
         dq0_drho(i_grid) = 0.0_dp
         dq0_dgradrho(i_grid) = 0.0_dp

         !! This prevents numerical problems.  If the charge density is negative (an
         !! unphysical situation), we simply treat it as very small.  In that case,
         !! q0 will be very large and will be saturated.  For a saturated q0 the derivative
         !! dq0_dq will be 0 so we set q0 = q_cut and dq0_drho = dq0_dgradrho = 0 and go on
         !! to the next point.
         !! ------------------------------------------------------------------------------------
         IF (total_rho(i_grid) < epsr) CYCLE
         !! ------------------------------------------------------------------------------------
         !! Calculate some intermediate values needed to find q
         !! ------------------------------------------------------------------------------------
         kF = (3.0_dp*pi*pi*total_rho(i_grid))**(1.0_dp/3.0_dp)
         r_s = (3.0_dp/(4.0_dp*pi*total_rho(i_grid)))**(1.0_dp/3.0_dp)
         sqrt_r_s = SQRT(r_s)

         gradient_correction = -Z_ab/(36.0_dp*kF*total_rho(i_grid)**2) &
                               *(gradient_rho(i_grid, 1)**2 + gradient_rho(i_grid, 2)**2 + gradient_rho(i_grid, 3)**2)

         LDA_1 = 8.0_dp*pi/3.0_dp*(LDA_A*(1.0_dp + LDA_a1*r_s))
         LDA_2 = 2.0_dp*LDA_A*(LDA_b1*sqrt_r_s + LDA_b2*r_s + LDA_b3*r_s*sqrt_r_s + LDA_b4*r_s*r_s)
         !! ---------------------------------------------------------------
         !! This is the q value defined in equations 11 and 12 of DION
         !! ---------------------------------------------------------------
         q = kF + LDA_1*LOG(1.0_dp + 1.0_dp/LDA_2) + gradient_correction
         !! ---------------------------------------------------------------
         !! Here, we calculate q0 by saturating q according to equation 5 of SOLER.  Also, we find
         !! the derivative dq0_dq needed for the derivatives dq0_drho and dq0_dgradrh0 discussed below.
         !! ---------------------------------------------------------------------------------------
         q__q_cut = q/q_cut
         CALL calculate_exponent_derivative(m_cut, q__q_cut, exponent, dq0_dq)
         q0(i_grid) = q_cut*(1.0_dp - EXP(-exponent))
         dq0_dq = dq0_dq*EXP(-exponent)
         !! ---------------------------------------------------------------------------------------
         !! This is to handle a case with q0 too small.  We simply set it to the smallest q value in
         !! out q_mesh.  Hopefully this doesn't get used often (ever)
         !! ---------------------------------------------------------------------------------------
         IF (q0(i_grid) < q_min) THEN
            q0(i_grid) = q_min
         END IF
         !! ---------------------------------------------------------------------------------------
         !! Here we find derivatives.  These are actually the density times the derivative of q0 with respect
         !! to rho and gradient_rho.  The density factor comes in since we are really differentiating
         !! theta = (rho)*P(q0) with respect to density (or its gradient) which will be
         !! dtheta_drho = P(q0) + dP_dq0 * [rho * dq0_dq * dq_drho]   and
         !! dtheta_dgradient_rho =  dP_dq0  * [rho * dq0_dq * dq_dgradient_rho]
         !! The parts in square brackets are what is calculated here.  The dP_dq0 term will be interpolated
         !! later.  There should actually be a factor of the magnitude of the gradient in the gradient_rho derivative
         !! but that cancels out when we differentiate the magnitude of the gradient with respect to a particular
         !! component.
         !! ------------------------------------------------------------------------------------------------

         dq0_drho(i_grid) = dq0_dq*(kF/3.0_dp - 7.0_dp/3.0_dp*gradient_correction &
                                    - 8.0_dp*pi/9.0_dp*LDA_A*LDA_a1*r_s*LOG(1.0_dp + 1.0_dp/LDA_2) &
                                    + LDA_1/(LDA_2*(1.0_dp + LDA_2)) &
                                    *(2.0_dp*LDA_A*(LDA_b1/6.0_dp*sqrt_r_s + LDA_b2/3.0_dp*r_s + LDA_b3/2.0_dp*r_s*sqrt_r_s &
                                                    + 2.0_dp*LDA_b4/3.0_dp*r_s**2)))

         dq0_dgradrho(i_grid) = total_rho(i_grid)*dq0_dq*2.0_dp*(-Z_ab)/(36.0_dp*kF*total_rho(i_grid)**2)

      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE get_q0_on_grid_vdw

! **************************************************************************************************
!> \brief ...
!> \param total_rho ...
!> \param gradient_rho ...
!> \param q0 ...
!> \param dq0_drho ...
!> \param dq0_dgradrho ...
!> \param dispersion_env ...
! **************************************************************************************************
   SUBROUTINE get_q0_on_grid_rvv10(total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho, dispersion_env)
      !!
      !! more specifically it calculates the following
      !!
      !!     q0(ir) = q0 as defined above
      !!     dq0_drho(ir) = total_rho * d q0 /d rho
      !!     dq0_dgradrho = total_rho / |gradient_rho| * d q0 / d |gradient_rho|
      !!
      REAL(dp), INTENT(IN)                               :: total_rho(:), gradient_rho(:, :)
      REAL(dp), INTENT(OUT)                              :: q0(:), dq0_drho(:), dq0_dgradrho(:)
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env

      INTEGER, PARAMETER                                 :: m_cut = 12

      INTEGER                                            :: i_grid
      REAL(dp)                                           :: b_value, C_value, dk_dn, dq0_dq, dw0_dn, &
                                                            exponent, gmod2, k, mod_grad, q, &
                                                            q__q_cut, q_cut, q_min, w0, wg2, wp2

      q_cut = dispersion_env%q_cut
      q_min = dispersion_env%q_min
      b_value = dispersion_env%b_value
      C_value = dispersion_env%c_value

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             SHARED(total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho, q_cut, q_min, b_value, C_value) &
!$OMP             PRIVATE(dk_dn, dq0_dq, dw0_dn, exponent, gmod2, k, mod_grad, q, q__q_cut, w0, wg2, wp2) &
!$OMP             SCHEDULE(STATIC)
      DO i_grid = 1, SIZE(total_rho)
         q0(i_grid) = q_cut
         dq0_drho(i_grid) = 0.0_dp
         dq0_dgradrho(i_grid) = 0.0_dp

         gmod2 = gradient_rho(i_grid, 1)**2 + gradient_rho(i_grid, 2)**2 + gradient_rho(i_grid, 3)**2

         !if (total_rho(i_grid) > epsr .and. gmod2 > epsr) cycle
         IF (total_rho(i_grid) > epsr) THEN

            !! Calculate some intermediate values needed to find q
            !! ------------------------------------------------------------------------------------
            mod_grad = SQRT(gmod2)

            wp2 = 16.0_dp*pi*total_rho(i_grid)
            wg2 = 4_dp*C_value*(mod_grad/total_rho(i_grid))**4

            k = b_value*3.0_dp*pi*((total_rho(i_grid)/(9.0_dp*pi))**(1.0_dp/6.0_dp))
            w0 = SQRT(wg2 + wp2/3.0_dp)

            q = w0/k

            !! Here, we calculate q0 by saturating q according
            !! ---------------------------------------------------------------------------------------
            q__q_cut = q/q_cut
            CALL calculate_exponent_derivative(m_cut, q__q_cut, exponent, dq0_dq)
            q0(i_grid) = q_cut*(1.0_dp - EXP(-exponent))
            dq0_dq = dq0_dq*EXP(-exponent)

            !! ---------------------------------------------------------------------------------------
            IF (q0(i_grid) < q_min) THEN
               q0(i_grid) = q_min
            END IF

            !!---------------------------------Final values---------------------------------
            dw0_dn = 1.0_dp/(2.0_dp*w0)*(16.0_dp/3.0_dp*pi - 4.0_dp*wg2/total_rho(i_grid))
            dk_dn = k/(6.0_dp*total_rho(i_grid))

            dq0_drho(i_grid) = dq0_dq*1.0_dp/(k**2.0)*(dw0_dn*k - dk_dn*w0)
            ! wg2 is proportional to |gradient_rho|**4, so this limit is zero.
            IF (gmod2 > 0.0_dp) THEN
               dq0_dgradrho(i_grid) = dq0_dq*1.0_dp/(2.0_dp*k*w0)*4.0_dp*wg2/gmod2
            END IF
         END IF

      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE get_q0_on_grid_rvv10

! **************************************************************************************************
!> \brief ...
!> \param total_rho ...
!> \param gradient_rho ...
!> \param q0 ...
!> \param dispersion_env ...
! **************************************************************************************************
   SUBROUTINE get_q0_on_grid_eo_vdw(total_rho, gradient_rho, q0, dispersion_env)

      REAL(dp), INTENT(IN)                               :: total_rho(:), gradient_rho(:, :)
      REAL(dp), INTENT(OUT)                              :: q0(:)
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env

      INTEGER, PARAMETER                                 :: m_cut = 12
      REAL(dp), PARAMETER :: LDA_A = 0.031091_dp, LDA_a1 = 0.2137_dp, LDA_b1 = 7.5957_dp, &
         LDA_b2 = 3.5876_dp, LDA_b3 = 1.6382_dp, LDA_b4 = 0.49294_dp

      INTEGER                                            :: i_grid
      REAL(dp)                                           :: exponent, gradient_correction, kF, &
                                                            LDA_1, LDA_2, q, q__q_cut, q_cut, &
                                                            q_min, r_s, sqrt_r_s, Z_ab

      q_cut = dispersion_env%q_cut
      q_min = dispersion_env%q_min
      SELECT CASE (dispersion_env%nl_type)
      CASE DEFAULT
         CPABORT("Unknown vdW-DF functional")
      CASE (vdw_nl_DRSLL)
         Z_ab = -0.8491_dp
      CASE (vdw_nl_LMKLL)
         Z_ab = -1.887_dp
      END SELECT

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             SHARED(total_rho, gradient_rho, q0, q_cut, q_min, Z_ab) &
!$OMP             PRIVATE(exponent, gradient_correction, kF, LDA_1, LDA_2, q, q__q_cut, r_s, sqrt_r_s) &
!$OMP             SCHEDULE(STATIC)
      DO i_grid = 1, SIZE(total_rho)
         q0(i_grid) = q_cut
         !! This prevents numerical problems.  If the charge density is negative (an
         !! unphysical situation), we simply treat it as very small.  In that case,
         !! q0 will be very large and will be saturated.  For a saturated q0 the derivative
         !! dq0_dq will be 0 so we set q0 = q_cut and dq0_drho = dq0_dgradrho = 0 and go on
         !! to the next point.
         !! ------------------------------------------------------------------------------------
         IF (total_rho(i_grid) < epsr) CYCLE
         !! ------------------------------------------------------------------------------------
         !! Calculate some intermediate values needed to find q
         !! ------------------------------------------------------------------------------------
         kF = (3.0_dp*pi*pi*total_rho(i_grid))**(1.0_dp/3.0_dp)
         r_s = (3.0_dp/(4.0_dp*pi*total_rho(i_grid)))**(1.0_dp/3.0_dp)
         sqrt_r_s = SQRT(r_s)

         gradient_correction = -Z_ab/(36.0_dp*kF*total_rho(i_grid)**2) &
                               *(gradient_rho(i_grid, 1)**2 + gradient_rho(i_grid, 2)**2 + gradient_rho(i_grid, 3)**2)

         LDA_1 = 8.0_dp*pi/3.0_dp*(LDA_A*(1.0_dp + LDA_a1*r_s))
         LDA_2 = 2.0_dp*LDA_A*(LDA_b1*sqrt_r_s + LDA_b2*r_s + LDA_b3*r_s*sqrt_r_s + LDA_b4*r_s*r_s)
         !! ------------------------------------------------------------------------------------
         !! This is the q value defined in equations 11 and 12 of DION
         !! ---------------------------------------------------------------
         q = kF + LDA_1*LOG(1.0_dp + 1.0_dp/LDA_2) + gradient_correction

         !! ---------------------------------------------------------------
         !! Here, we calculate q0 by saturating q according to equation 5 of SOLER.  Also, we find
         !! the derivative dq0_dq needed for the derivatives dq0_drho and dq0_dgradrh0 discussed below.
         !! ---------------------------------------------------------------------------------------
         q__q_cut = q/q_cut
         CALL calculate_exponent(m_cut, q__q_cut, exponent)
         q0(i_grid) = q_cut*(1.0_dp - EXP(-exponent))

         !! ---------------------------------------------------------------------------------------
         !! This is to handle a case with q0 too small.  We simply set it to the smallest q value in
         !! out q_mesh.  Hopefully this doesn't get used often (ever)
         !! ---------------------------------------------------------------------------------------
         IF (q0(i_grid) < q_min) THEN
            q0(i_grid) = q_min
         END IF
      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE get_q0_on_grid_eo_vdw

! **************************************************************************************************
!> \brief ...
!> \param total_rho ...
!> \param gradient_rho ...
!> \param q0 ...
!> \param dispersion_env ...
! **************************************************************************************************
   SUBROUTINE get_q0_on_grid_eo_rvv10(total_rho, gradient_rho, q0, dispersion_env)

      REAL(dp), INTENT(IN)                               :: total_rho(:), gradient_rho(:, :)
      REAL(dp), INTENT(OUT)                              :: q0(:)
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env

      INTEGER, PARAMETER                                 :: m_cut = 12

      INTEGER                                            :: i_grid
      REAL(dp)                                           :: b_value, C_value, exponent, gmod2, k, q, &
                                                            q__q_cut, q_cut, q_min, w0, wg2, wp2

      q_cut = dispersion_env%q_cut
      q_min = dispersion_env%q_min
      b_value = dispersion_env%b_value
      C_value = dispersion_env%c_value

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             SHARED(total_rho, gradient_rho, q0, q_cut, q_min, b_value, C_value) &
!$OMP             PRIVATE(exponent, gmod2, k, q, q__q_cut, w0, wg2, wp2) &
!$OMP             SCHEDULE(STATIC)
      DO i_grid = 1, SIZE(total_rho)
         q0(i_grid) = q_cut

         gmod2 = gradient_rho(i_grid, 1)**2 + gradient_rho(i_grid, 2)**2 + gradient_rho(i_grid, 3)**2

         !if (total_rho(i_grid) > epsr .and. gmod2 > epsr) cycle
         IF (total_rho(i_grid) > epsr) THEN

            !! Calculate some intermediate values needed to find q
            !! ------------------------------------------------------------------------------------
            wp2 = 16.0_dp*pi*total_rho(i_grid)
            wg2 = 4_dp*C_value*(gmod2*gmod2)/(total_rho(i_grid)**4)

            k = b_value*3.0_dp*pi*((total_rho(i_grid)/(9.0_dp*pi))**(1.0_dp/6.0_dp))
            w0 = SQRT(wg2 + wp2/3.0_dp)

            q = w0/k

            !! Here, we calculate q0 by saturating q according
            !! ---------------------------------------------------------------------------------------
            q__q_cut = q/q_cut
            CALL calculate_exponent(m_cut, q__q_cut, exponent)
            q0(i_grid) = q_cut*(1.0_dp - EXP(-exponent))

            IF (q0(i_grid) < q_min) THEN
               q0(i_grid) = q_min
            END IF

         END IF

      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE get_q0_on_grid_eo_rvv10

! **************************************************************************************************
!> \brief Prepare the spline intervals, weights and density factors for streamed theta construction.
!> \param q0 Saturated interpolation points
!> \param rho Total density
!> \param dispersion_env Non-local functional parameters
!> \param q_low Lower spline endpoint; zero denotes an inactive rVV10 point
!> \param coeff Cubic spline coefficients a, b, c, d
!> \param theta_scale Density factor multiplying the interpolated spline
! **************************************************************************************************
   SUBROUTINE prepare_splines(q0, rho, dispersion_env, q_low, coeff, theta_scale)
      REAL(dp), INTENT(IN)                               :: q0(:), rho(:)
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      INTEGER, INTENT(OUT)                               :: q_low(:)
      REAL(dp), INTENT(OUT)                              :: coeff(:, :), theta_scale(:)

      INTEGER                                            :: i, j, lower, nqs, upper
      LOGICAL                                            :: rvv10
      REAL(dp)                                           :: a, b, const, dx, dx2_6
      REAL(dp), POINTER                                  :: q_mesh(:)

      q_mesh => dispersion_env%q_mesh
      nqs = dispersion_env%nqs
      CPASSERT(nqs >= 2)
      rvv10 = dispersion_env%nl_type == vdw_nl_RVV10
      const = 1.0_dp/(3.0_dp*rootpi*dispersion_env%b_value**1.5_dp)/(pi**0.75_dp)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             SHARED(q0, rho, q_mesh, nqs, rvv10, const, q_low, coeff, theta_scale) &
!$OMP             PRIVATE(j, lower, upper, A, b, dx, dx2_6) SCHEDULE(STATIC)
      DO i = 1, SIZE(q0)
         IF (rvv10 .AND. rho(i) <= epsr) THEN
            q_low(i) = 0
            coeff(i, :) = 0.0_dp
            theta_scale(i) = 0.0_dp
            CYCLE
         END IF
         lower = 1
         upper = nqs
         DO WHILE (upper - lower > 1)
            j = (upper + lower)/2
            IF (q0(i) > q_mesh(j)) THEN
               lower = j
            ELSE
               upper = j
            END IF
         END DO
         q_low(i) = lower
         dx = q_mesh(upper) - q_mesh(lower)
         dx2_6 = dx*dx/6.0_dp
         a = (q_mesh(upper) - q0(i))/dx
         b = (q0(i) - q_mesh(lower))/dx
         coeff(i, 1) = a
         coeff(i, 2) = b
         coeff(i, 3) = (a**3 - a)*dx2_6
         coeff(i, 4) = (b**3 - b)*dx2_6
         IF (rvv10) THEN
            theta_scale(i) = const*rho(i)**0.75_dp
         ELSE
            theta_scale(i) = rho(i)
         END IF
      END DO
!$OMP END PARALLEL DO
   END SUBROUTINE prepare_splines

! **************************************************************************************************
!> \brief Generate one real-space theta channel directly in the existing FFT workspace.
!> \param iq Channel index
!> \param q_low Lower spline endpoint
!> \param coeff Spline coefficients
!> \param theta_scale Density factor
!> \param dispersion_env Non-local functional parameters
!> \param theta FFT input, viewed as a contiguous one-dimensional array
! **************************************************************************************************
   SUBROUTINE build_theta(iq, q_low, coeff, theta_scale, dispersion_env, theta)
      INTEGER, INTENT(IN)                                :: iq, q_low(:)
      REAL(dp), INTENT(IN)                               :: coeff(:, :), theta_scale(:)
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      REAL(dp), INTENT(OUT)                              :: theta(:)

      INTEGER                                            :: i, lower
      REAL(dp)                                           :: p
      REAL(dp), POINTER                                  :: d2y(:, :)

      d2y => dispersion_env%d2y_dx2
!$OMP PARALLEL DO SIMD DEFAULT(NONE) &
!$OMP             SHARED(iq, q_low, coeff, theta_scale, d2y, theta) PRIVATE(lower, p) SCHEDULE(STATIC)
      DO i = 1, SIZE(theta)
         lower = q_low(i)
         theta(i) = 0.0_dp
         IF (lower == 0) CYCLE
         p = coeff(i, 1)*MERGE(1.0_dp, 0.0_dp, iq == lower) &
             + coeff(i, 2)*MERGE(1.0_dp, 0.0_dp, iq == lower + 1) &
             + (coeff(i, 3)*d2y(iq, lower) + coeff(i, 4)*d2y(iq, lower + 1))
         theta(i) = p*theta_scale(i)
      END DO
!$OMP END PARALLEL DO SIMD
   END SUBROUTINE build_theta

! **************************************************************************************************
!> \brief Accumulate one inverse-transformed channel without storing a real-space channel matrix.
!> \param iq Channel index
!> \param q_low Lower spline endpoint, using the potential-side knot convention
!> \param dispersion_env Non-local functional parameters
!> \param u Inverse FFT output
!> \param u_contract Sums against the two second-derivative columns, then the two endpoint values
! **************************************************************************************************
   SUBROUTINE accumulate_potential(iq, q_low, dispersion_env, u, u_contract)
      INTEGER, INTENT(IN)                                :: iq, q_low(:)
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      REAL(dp), INTENT(IN)                               :: u(:)
      REAL(dp), INTENT(INOUT)                            :: u_contract(:, :)

      INTEGER                                            :: i, lower
      REAL(dp), POINTER                                  :: d2y(:, :)

      d2y => dispersion_env%d2y_dx2
!$OMP PARALLEL DO SIMD DEFAULT(NONE) &
!$OMP             SHARED(iq, q_low, d2y, u, u_contract) PRIVATE(lower) SCHEDULE(STATIC)
      DO i = 1, SIZE(u)
         lower = q_low(i)
         IF (lower == 0) CYCLE
         u_contract(i, 1) = u_contract(i, 1) + u(i)*d2y(iq, lower)
         u_contract(i, 2) = u_contract(i, 2) + u(i)*d2y(iq, lower + 1)
         IF (iq == lower) u_contract(i, 3) = u(i)
         IF (iq == lower + 1) u_contract(i, 4) = u(i)
      END DO
!$OMP END PARALLEL DO SIMD
   END SUBROUTINE accumulate_potential

! **************************************************************************************************
!> \brief This routine is modeled after an algorithm from "Numerical Recipes in C" by Cambridge
!> University Press, pages 96-97.  It was adapted for Fortran and for the problem at hand.
!> \param x ...
!> \param d2y_dx2 ...
!> \par History
!>     OpenMP added: Aug 2016  MTucker
! **************************************************************************************************
   SUBROUTINE initialize_spline_interpolation(x, d2y_dx2)

      REAL(dp), INTENT(in)                               :: x(:)
      REAL(dp), INTENT(inout)                            :: d2y_dx2(:, :)

      INTEGER                                            :: index, Nx, P_i
      REAL(dp)                                           :: temp1, temp2
      REAL(dp), ALLOCATABLE                              :: temp_array(:), y(:)

      Nx = SIZE(x)

!$OMP PARALLEL DEFAULT( NONE )                  &
!$OMP           SHARED( x, d2y_dx2, Nx )        &
!$OMP          PRIVATE( temp_array, y           &
!$OMP                 , index, temp1, temp2     &
!$OMP                 )

      ALLOCATE (temp_array(Nx), y(Nx))

!$OMP DO
      DO P_i = 1, Nx
         !! In the Soler method, the polynomials that are interpolated are Kronecker delta functions
         !! at a particular q point.  So, we set all y values to 0 except the one corresponding to
         !! the particular function P_i.
         !! ----------------------------------------------------------------------------------------
         y = 0.0_dp
         y(P_i) = 1.0_dp
         !! ----------------------------------------------------------------------------------------

         d2y_dx2(P_i, 1) = 0.0_dp
         temp_array(1) = 0.0_dp
         DO index = 2, Nx - 1
            temp1 = (x(index) - x(index - 1))/(x(index + 1) - x(index - 1))
            temp2 = temp1*d2y_dx2(P_i, index - 1) + 2.0_dp
            d2y_dx2(P_i, index) = (temp1 - 1.0_dp)/temp2
            temp_array(index) = (y(index + 1) - y(index))/(x(index + 1) - x(index)) &
                                - (y(index) - y(index - 1))/(x(index) - x(index - 1))
            temp_array(index) = (6.0_dp*temp_array(index)/(x(index + 1) - x(index - 1)) &
                                 - temp1*temp_array(index - 1))/temp2
         END DO
         d2y_dx2(P_i, Nx) = 0.0_dp
         DO index = Nx - 1, 1, -1
            d2y_dx2(P_i, index) = d2y_dx2(P_i, index)*d2y_dx2(P_i, index + 1) + temp_array(index)
         END DO
      END DO
!$OMP END DO

      DEALLOCATE (temp_array, y)
!$OMP END PARALLEL

   END SUBROUTINE initialize_spline_interpolation

! **************************************************************************************************
!> \brief Interpolate the symmetric kernel and optionally its radial derivative from packed tables.
!> \param k Reciprocal-vector length
!> \param kernel_of_k Interpolated kernel matrix
!> \param dispersion_env Radial mesh and channel count
!> \param kernel_table Packed channel pairs, radial points, and value/second derivative
!> \param dkernel_of_dk Optional derivative matrix for the kernel virial
! **************************************************************************************************
   SUBROUTINE interpolate_kernel(k, kernel_of_k, dispersion_env, kernel_table, dkernel_of_dk)
      REAL(dp), INTENT(IN)                               :: k
      REAL(dp), INTENT(OUT)                              :: kernel_of_k(:, :)
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      REAL(dp), INTENT(IN)                               :: kernel_table(:, 0:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: dkernel_of_dk(:, :)

      INTEGER                                            :: ipair, k_i, q1_i, q2_i
      LOGICAL                                            :: on_mesh
      REAL(dp)                                           :: a, b, c, d, da, db, dc, dd, dk, dk_6, &
                                                            value

      dk = dispersion_env%dk
      CPASSERT(k < dispersion_env%nr_points*dk)
      k_i = INT(k/dk)
      on_mesh = MOD(k, dk) == 0.0_dp
      a = (dk*(k_i + 1.0_dp) - k)/dk
      b = (k - dk*k_i)/dk
      c = (a**3 - a)*(dk*dk/6.0_dp)
      d = (b**3 - b)*(dk*dk/6.0_dp)
      DO q1_i = 1, dispersion_env%nqs
!$OMP SIMD PRIVATE(ipair, value)
         DO q2_i = 1, q1_i
            ipair = q1_i*(q1_i - 1)/2 + q2_i
            IF (on_mesh) THEN
               value = kernel_table(ipair, k_i, 1)
            ELSE
               value = a*kernel_table(ipair, k_i, 1) + b*kernel_table(ipair, k_i + 1, 1) &
                       + (c*kernel_table(ipair, k_i, 2) + d*kernel_table(ipair, k_i + 1, 2))
            END IF
            kernel_of_k(q2_i, q1_i) = value
            kernel_of_k(q1_i, q2_i) = value
         END DO
!$OMP END SIMD
      END DO
      IF (PRESENT(dkernel_of_dk)) THEN
         dk_6 = dk/6.0_dp
         da = -1.0_dp/dk
         db = 1.0_dp/dk
         dc = -(3*a**2 - 1.0_dp)*dk_6
         dd = (3*b**2 - 1.0_dp)*dk_6
         DO q1_i = 1, dispersion_env%nqs
!$OMP SIMD PRIVATE(ipair, value)
            DO q2_i = 1, q1_i
               ipair = q1_i*(q1_i - 1)/2 + q2_i
               value = da*kernel_table(ipair, k_i, 1) + db*kernel_table(ipair, k_i + 1, 1) &
                       + dc*kernel_table(ipair, k_i, 2) + dd*kernel_table(ipair, k_i + 1, 2)
               dkernel_of_dk(q2_i, q1_i) = value
               dkernel_of_dk(q1_i, q2_i) = value
            END DO
!$OMP END SIMD
         END DO
      END IF
   END SUBROUTINE interpolate_kernel

! **************************************************************************************************

END MODULE qs_dispersion_nonloc
