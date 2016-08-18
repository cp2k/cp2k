!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
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
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE input_constants,                 ONLY: vdw_nl_DRSLL,&
                                              vdw_nl_LMKLL,&
                                              vdw_nl_RVV10,&
                                              xc_vdw_fun_nonloc
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: pi
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_sum
   USE pw_grid_types,                   ONLY: HALFSPACE,&
                                              pw_grid_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_derive,&
                                              pw_transfer
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type,&
                                              pw_type
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
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_dispersion_nonloc_init', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: filename
      INTEGER                                            :: funit, handle, nqs, nr_points, q1_i, &
                                                            q2_i, vdw_type

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
         IF (para_env%source == para_env%mepos) THEN
            ! Read the kernel information from file "filename"
            CALL open_file(file_name=filename, unit_number=funit, file_form="FORMATTED")
            READ (funit, *) nqs, nr_points
            READ (funit, *) dispersion_env%r_max
         END IF
         CALL mp_bcast(nqs, para_env%source, para_env%group)
         CALL mp_bcast(nr_points, para_env%source, para_env%group)
         CALL mp_bcast(dispersion_env%r_max, para_env%source, para_env%group)
         ALLOCATE (dispersion_env%q_mesh(nqs), dispersion_env%kernel(0:nr_points, nqs, nqs), &
                   dispersion_env%d2phi_dk2(0:nr_points, nqs, nqs))
         dispersion_env%nqs = nqs
         dispersion_env%nr_points = nr_points
         IF (para_env%source == para_env%mepos) THEN
            !! Read in the values of the q points used to generate this kernel
            READ (funit, "(1p, 4e23.14)") dispersion_env%q_mesh
            !! For each pair of q values, read in the function phi_q1_q2(k).
            !! That is, the fourier transformed kernel function assuming q1 and q2
            !! for all the values of r used.
            DO q1_i = 1, nqs
               DO q2_i = 1, q1_i
                  READ (funit, "(1p, 4e23.14)") dispersion_env%kernel(0:nr_points, q1_i, q2_i)
                  dispersion_env%kernel(0:nr_points, q2_i, q1_i) = dispersion_env%kernel(0:nr_points, q1_i, q2_i)
               END DO
            END DO
            !! Again, for each pair of q values (q1 and q2), read in the value
            !! of the second derivative of the above mentiond Fourier transformed
            !! kernel function phi_alpha_beta(k).  These are used for spline
            !! interpolation of the Fourier transformed kernel.
            DO q1_i = 1, nqs
               DO q2_i = 1, q1_i
                  READ (funit, "(1p, 4e23.14)") dispersion_env%d2phi_dk2(0:nr_points, q1_i, q2_i)
                  dispersion_env%d2phi_dk2(0:nr_points, q2_i, q1_i) = dispersion_env%d2phi_dk2(0:nr_points, q1_i, q2_i)
               END DO
            END DO
            CALL close_file(unit_number=funit)
         END IF
         CALL mp_bcast(dispersion_env%q_mesh, para_env%source, para_env%group)
         CALL mp_bcast(dispersion_env%kernel, para_env%source, para_env%group)
         CALL mp_bcast(dispersion_env%d2phi_dk2, para_env%source, para_env%group)
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
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rho, rho_r, rho_g
      REAL(KIND=dp), INTENT(OUT)                         :: edispersion
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      LOGICAL, INTENT(IN)                                :: energy_only
      TYPE(pw_pool_type), POINTER                        :: pw_pool, xc_pw_pool
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(virial_type), OPTIONAL, POINTER               :: virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_dispersion_nonloc', &
         routineP = moduleN//':'//routineN
      INTEGER, DIMENSION(3, 3), PARAMETER :: nd = RESHAPE((/1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3, 3/))

      INTEGER                                            :: handle, i, i_grid, idir, ispin, n1, n2, &
                                                            n3, np, nspin
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: b_value, beta, Ec_nl, sumnp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dq0_dgradrho, dq0_drho, hpot, potential, &
                                                            q0, rho
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: drho, thetas
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_p_type), ALLOCATABLE, DIMENSION(:)         :: thetas_g
      TYPE(pw_p_type), ALLOCATABLE, DIMENSION(:, :)      :: drho_r
      TYPE(pw_type), POINTER                             :: tmp_g, tmp_r, vxc_g, vxc_r

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(rho_r))
      CPASSERT(ASSOCIATED(rho_g))
      CPASSERT(ASSOCIATED(pw_pool))

      IF (PRESENT(virial)) THEN
         use_virial = virial%pv_calculate .AND. (.NOT. virial%pv_numer)
      ELSE
         use_virial = .FALSE.
      ENDIF
      IF (use_virial) THEN
         CPASSERT(.NOT. energy_only)
      END IF
      IF (.NOT. energy_only) THEN
         CPASSERT(ASSOCIATED(vxc_rho))
      END IF

      b_value = dispersion_env%b_value
      beta = 0.03125_dp*(3.0_dp/(b_value**2.0_dp))**0.75_dp
      nspin = SIZE(rho_r)

      ! temporary arrays for FFT
      CALL pw_pool_create_pw(pw_pool, tmp_g, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(pw_pool, tmp_r, use_data=REALDATA3D, in_space=REALSPACE)

      ! get density derivatives
      ALLOCATE (drho_r(3, nspin))
      DO ispin = 1, nspin
         DO idir = 1, 3
            NULLIFY (drho_r(idir, ispin)%pw)
            CALL pw_pool_create_pw(pw_pool, drho_r(idir, ispin)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_transfer(rho_g(ispin)%pw, tmp_g)
            CALL pw_derive(tmp_g, nd(:, idir))
            CALL pw_transfer(tmp_g, drho_r(idir, ispin)%pw)
         END DO
      END DO

      np = SIZE(tmp_r%cr3d)
      ALLOCATE (rho(np), drho(np, 3))
      rho = 0.0_dp
      drho = 0.0_dp
      DO ispin = 1, nspin
         CPASSERT(rho_r(ispin)%pw%in_use == REALDATA3D)
         CALL pw_transfer(rho_g(ispin)%pw, tmp_g)
         CALL pw_transfer(tmp_g, tmp_r)
         rho(:) = rho(:)+RESHAPE(tmp_r%cr3d, (/np/))
         drho(:, 1) = drho(:, 1)+RESHAPE(drho_r(1, ispin)%pw%cr3d, (/np/))
         drho(:, 2) = drho(:, 2)+RESHAPE(drho_r(2, ispin)%pw%cr3d, (/np/))
         drho(:, 3) = drho(:, 3)+RESHAPE(drho_r(3, ispin)%pw%cr3d, (/np/))
      END DO
      !! ---------------------------------------------------------------------------------
      !! Find the value of q0 for all assigned grid points.  q is defined in equations
      !! 11 and 12 of DION and q0 is the saturated version of q defined in equation
      !! 5 of SOLER.  This routine also returns the derivatives of the q0s with respect
      !! to the charge-density and the gradient of the charge-density.  These are needed
      !! for the potential calculated below.
      !! ---------------------------------------------------------------------------------

      IF (energy_only) THEN
         ALLOCATE (q0(np))
         SELECT CASE (dispersion_env%nl_type)
         CASE DEFAULT
            CPABORT("Unknown vdW-DF functional")
         CASE (vdw_nl_DRSLL, vdw_nl_LMKLL)
            CALL get_q0_on_grid_eo_vdw(rho, drho, q0, dispersion_env)
         CASE (vdw_nl_RVV10)
            CALL get_q0_on_grid_eo_rvv10(rho, drho, q0, dispersion_env)
         END SELECT
      ELSE
         ALLOCATE (q0(np), dq0_drho(np), dq0_dgradrho(np))
         SELECT CASE (dispersion_env%nl_type)
         CASE DEFAULT
            CPABORT("Unknown vdW-DF functional")
         CASE (vdw_nl_DRSLL, vdw_nl_LMKLL)
            CALL get_q0_on_grid_vdw(rho, drho, q0, dq0_drho, dq0_dgradrho, dispersion_env)
         CASE (vdw_nl_RVV10)
            CALL get_q0_on_grid_rvv10(rho, drho, q0, dq0_drho, dq0_dgradrho, dispersion_env)
         END SELECT
      END IF

      !! ---------------------------------------------------------------------------------------------
      !! Here we allocate and calculate the theta functions appearing in equations 8-12 of SOLER.
      !! They are defined as rho*P_i(q0(rho, gradient_rho)) for vdW and vdW2 or
      !! constant*rho^(3/4)*P_i(q0(rho, gradient_rho)) for rVV10 where P_i is a polynomial that
      !! interpolates a Kronecker delta function at the point q_i (taken from the q_mesh) and q0 is
      !! the saturated version of q.
      !! q is defined in equations 11 and 12 of DION and the saturation proceedure is defined in equation 5
      !! of soler.  This is the biggest memory consumer in the method since the thetas array is
      !! (total # of FFT points)*Nqs complex numbers.  In a parallel run, each processor will hold the
      !! values of all the theta functions on just the points assigned to it.
      !! --------------------------------------------------------------------------------------------------
      !! thetas are stored in reciprocal space as  theta_i(k) because this is the way they are used later
      !! for the convolution (equation 11 of SOLER).
      !! --------------------------------------------------------------------------------------------------
      ALLOCATE (thetas(np, dispersion_env%nqs))
      !! Interpolate the P_i polynomials defined in equation 3 in SOLER for the particular
      !! q0 values we have.
      CALL spline_interpolation(dispersion_env%q_mesh, dispersion_env%d2y_dx2, q0, thetas)

      !! Form the thetas where theta is defined as rho*p_i(q0) for vdW and vdW2 or
      !! constant*rho^(3/4)*p_i(q0) for rVV10
      !! ------------------------------------------------------------------------------------
      SELECT CASE (dispersion_env%nl_type)
      CASE DEFAULT
         CPABORT("Unknown vdW-DF functional")
      CASE (vdw_nl_DRSLL, vdw_nl_LMKLL)
         DO i = 1, dispersion_env%nqs
            thetas(:, i) = thetas(:, i)*rho(:)
         END DO
      CASE (vdw_nl_RVV10)
         DO i_grid = 1, np
            IF (rho(i_grid) > epsr) THEN
               DO i = 1, dispersion_env%nqs
                  thetas(i_grid, i) = thetas(i_grid, i)*(1.0_dp/(3.0_dp*SQRT(pi)* &
                                                                 (b_value**(3.0_dp/2.0_dp))))*(rho(i_grid)/pi)**(3.0_dp/4.0_dp)
               END DO
            ELSE
               thetas(i_grid, :) = 0.0_dp
            END IF
         END DO
      END SELECT
      !! ------------------------------------------------------------------------------------
      !! Get thetas in reciprocal space.
      n1 = SIZE(tmp_r%cr3d, 1)
      n2 = SIZE(tmp_r%cr3d, 2)
      n3 = SIZE(tmp_r%cr3d, 3)
      ALLOCATE (thetas_g(dispersion_env%nqs))
      DO i = 1, dispersion_env%nqs
         tmp_r%cr3d = RESHAPE(thetas(:, i), (/n1, n2, n3/))
         NULLIFY (thetas_g(i)%pw)
         CALL pw_pool_create_pw(pw_pool, thetas_g(i)%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         CALL pw_transfer(tmp_r, thetas_g(i)%pw)
      END DO
      grid => thetas_g(1)%pw%pw_grid
      !! ---------------------------------------------------------------------------------------------
      !! Carry out the integration in equation 8 of SOLER.  This also turns the thetas array into the
      !! precursor to the u_i(k) array which is inverse fourier transformed to get the u_i(r) functions
      !! of SOLER equation 11.  Add the energy we find to the output variable etxc.
      !! --------------------------------------------------------------------------------------------------
      sumnp = np
      CALL mp_sum(sumnp, para_env%group)
      IF (use_virial) THEN
         ! calculates kernel contribution to stress
         CALL vdW_energy(thetas_g, dispersion_env, Ec_nl, energy_only, virial)
         SELECT CASE (dispersion_env%nl_type)
         CASE (vdw_nl_RVV10)
            Ec_nl = 0.5_dp*Ec_nl+beta*SUM(rho(:))*grid%vol/sumnp
         END SELECT
         ! calculates energy contribution to stress
         ! potential contribution to stress is calculated together with other potentials (Hxc)
         DO idir = 1, 3
            virial%pv_xc(idir, idir) = virial%pv_xc(idir, idir)+Ec_nl
         END DO
      ELSE
         CALL vdW_energy(thetas_g, dispersion_env, Ec_nl, energy_only)
         SELECT CASE (dispersion_env%nl_type)
         CASE (vdw_nl_RVV10)
            Ec_nl = 0.5_dp*Ec_nl+beta*SUM(rho(:))*grid%vol/sumnp
         END SELECT
      END IF
      CALL mp_sum(Ec_nl, para_env%group)
      edispersion = Ec_nl

      IF (.NOT. energy_only) THEN
         !! ----------------------------------------------------------------------------
         !! Inverse Fourier transform the u_i(k) to get the u_i(r) of SOLER equation 11.
         !!-----------------------------------------------------------------------------
         DO i = 1, dispersion_env%nqs
            CALL pw_transfer(thetas_g(i)%pw, tmp_r)
            thetas(:, i) = RESHAPE(tmp_r%cr3d, (/np/))
         END DO
         !! -------------------------------------------------------------------------
         !! Here we allocate the array to hold the potential. This is calculated via
         !! equation 10 of SOLER, using the u_i(r) calculated from equations 11 and
         !! 12 of SOLER.  Each processor allocates the array to be the size of the
         !! full grid  because, as can be seen in SOLER equation 9, processors need
         !! to access grid points outside their allocated regions.
         !! -------------------------------------------------------------------------
         ALLOCATE (potential(np), hpot(np))
         IF (use_virial) THEN
            ! calculates gradient contribution to stress
            grid => tmp_g%pw_grid
            CALL get_potential(q0, dq0_drho, dq0_dgradrho, rho, thetas, potential, hpot, &
                               dispersion_env, drho, grid%dvol, virial)
         ELSE
            CALL get_potential(q0, dq0_drho, dq0_dgradrho, rho, thetas, potential, hpot, &
                               dispersion_env)
         END IF
         SELECT CASE (dispersion_env%nl_type)
         CASE (vdw_nl_RVV10)
            potential(:) = 0.5_dp*potential(:)+beta
            hpot(:) = 0.5_dp*hpot(:)
         END SELECT
         CALL pw_pool_create_pw(pw_pool, vxc_r, use_data=REALDATA3D, in_space=REALSPACE)
         vxc_r%cr3d = RESHAPE(potential, (/n1, n2, n3/))
         DO idir = 1, 3
            tmp_r%cr3d = 0.0_dp
            DO ispin = 1, nspin
               tmp_r%cr3d = tmp_r%cr3d+RESHAPE(hpot, (/n1, n2, n3/))*drho_r(idir, ispin)%pw%cr3d
            END DO
            CALL pw_transfer(tmp_r, tmp_g)
            CALL pw_derive(tmp_g, nd(:, idir))
            CALL pw_transfer(tmp_g, tmp_r)
            CALL pw_axpy(tmp_r, vxc_r, -1._dp)
         END DO
         CALL pw_transfer(vxc_r, tmp_g)
         CALL pw_pool_give_back_pw(pw_pool, vxc_r)
         CALL pw_pool_create_pw(xc_pw_pool, vxc_r, use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_pool_create_pw(xc_pw_pool, vxc_g, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         CALL pw_transfer(tmp_g, vxc_g)
         CALL pw_transfer(vxc_g, vxc_r)
         DO ispin = 1, nspin
            CALL pw_axpy(vxc_r, vxc_rho(ispin)%pw, 1._dp)
         END DO
         CALL pw_pool_give_back_pw(xc_pw_pool, vxc_r)
         CALL pw_pool_give_back_pw(xc_pw_pool, vxc_g)
      END IF

      DEALLOCATE (thetas)

      DO i = 1, dispersion_env%nqs
         CALL pw_pool_give_back_pw(pw_pool, thetas_g(i)%pw)
      END DO
      DO ispin = 1, nspin
         DO idir = 1, 3
            CALL pw_pool_give_back_pw(pw_pool, drho_r(idir, ispin)%pw)
         END DO
      END DO
      CALL pw_pool_give_back_pw(pw_pool, tmp_r)
      CALL pw_pool_give_back_pw(pw_pool, tmp_g)

      DEALLOCATE (rho, drho, drho_r, thetas_g)
      IF (energy_only) THEN
         DEALLOCATE (q0)
      ELSE
         DEALLOCATE (q0, dq0_drho, dq0_dgradrho)
      END IF

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
      TYPE(pw_p_type), ALLOCATABLE, DIMENSION(:)         :: thetas_g
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      REAL(KIND=dp), INTENT(OUT)                         :: vdW_xc_energy
      LOGICAL, INTENT(IN)                                :: energy_only
      TYPE(virial_type), OPTIONAL, POINTER               :: virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'vdW_energy', routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: uu
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: theta
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: u_vdw(:, :)
      INTEGER                                            :: handle, ig, iq, l, m, nl_type, nqs, &
                                                            q1_i, q2_i
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: g, g2, g2_last, g_multiplier, gm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dkernel_of_dk, kernel_of_k
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_thread
      TYPE(pw_grid_type), POINTER                        :: grid

      CALL timeset(routineN, handle)
      nqs = dispersion_env%nqs

      IF (PRESENT(virial)) THEN
         use_virial = .TRUE.
         virial_thread(:, :) = 0.0_dp
      ELSE
         use_virial = .FALSE.
      END IF

      vdW_xc_energy = 0._dp
      grid => thetas_g(1)%pw%pw_grid

      IF (grid%grid_span == HALFSPACE) THEN
         g_multiplier = 2._dp
      ELSE
         g_multiplier = 1._dp
      END IF

      nl_type = dispersion_env%nl_type

      IF (.NOT. energy_only) THEN
         ALLOCATE (u_vdW(grid%ngpts_cut_local, nqs))
         u_vdW(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      END IF

!$OMP PARALLEL DEFAULT( NONE )                                             &
!$OMP           SHARED( nqs, energy_only, grid, dispersion_env             &
!$OMP                 , use_virial, thetas_g, g_multiplier, nl_type        &
!$OMP                 , u_vdW                                              &
!$OMP                 )                                                    &
!$OMP          PRIVATE( g2_last, kernel_of_k, dkernel_of_dk, theta         &
!$OMP                 , g2, g, iq                                          &
!$OMP                 , q2_i, uu, q1_i, gm, l, m                           &
!$OMP                 )                                                    &
!$OMP        REDUCTION( +: vdW_xc_energy, virial_thread                    &
!$OMP                 )

      g2_last = HUGE(0._dp)

      ALLOCATE (kernel_of_k(nqs, nqs), dkernel_of_dk(nqs, nqs))
      ALLOCATE (theta(nqs))

!$OMP DO
      DO ig = 1, grid%ngpts_cut_local
         g2 = grid%gsq(ig)
         IF (ABS(g2-g2_last) > 1.e-10) THEN
            g2_last = g2
            g = SQRT(g2)
            CALL interpolate_kernel(g, kernel_of_k, dispersion_env)
            IF (use_virial) CALL interpolate_dkernel_dk(g, dkernel_of_dk, dispersion_env)
         END IF
         DO iq = 1, nqs
            theta(iq) = thetas_g(iq)%pw%cc(ig)
         END DO
         DO q2_i = 1, nqs
            uu = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            DO q1_i = 1, nqs
               uu = uu+kernel_of_k(q2_i, q1_i)*theta(q1_i)
            END DO
            IF (ig < grid%first_gne0) THEN
               vdW_xc_energy = vdW_xc_energy+REAL(uu*CONJG(theta(q2_i)), KIND=dp)
            ELSE
               vdW_xc_energy = vdW_xc_energy+g_multiplier*REAL(uu*CONJG(theta(q2_i)), KIND=dp)
            END IF
            IF (.NOT. energy_only) u_vdW(ig, q2_i) = uu

            IF (use_virial .AND. ig >= grid%first_gne0) THEN
               DO q1_i = 1, nqs
                  gm = 0.5_dp*g_multiplier*grid%vol*dkernel_of_dk(q1_i, q2_i)*REAL(theta(q1_i)*CONJG(theta(q2_i)), KIND=dp)
                  IF (nl_type == vdw_nl_RVV10) THEN
                     gm = 0.5_dp*gm
                  END IF
                  DO l = 1, 3
                     DO m = 1, l
                        virial_thread(l, m) = virial_thread(l, m)-gm*(grid%g(l, ig)*grid%g(m, ig))/g
                     END DO
                  END DO
               END DO
            END IF

         END DO
      END DO
!$OMP END DO

      DEALLOCATE (theta)
      DEALLOCATE (kernel_of_k, dkernel_of_dk)

      IF (.NOT. energy_only) THEN
!$OMP DO
         DO ig = 1, grid%ngpts_cut_local
            DO iq = 1, nqs
               thetas_g(iq)%pw%cc(ig) = u_vdW(ig, iq)
            END DO
         END DO
!$OMP END DO
      END IF

!$OMP END PARALLEL

      IF (.NOT. energy_only) THEN
         DEALLOCATE (u_vdW)
      END IF

      vdW_xc_energy = vdW_xc_energy*grid%vol*0.5_dp

      IF (use_virial) THEN
         DO l = 1, 3
            DO m = 1, (l-1)
               virial%pv_xc(l, m) = virial%pv_xc(l, m)+virial_thread(l, m)
               virial%pv_xc(m, l) = virial%pv_xc(l, m)
            END DO
            m = l
            virial%pv_xc(l, m) = virial%pv_xc(l, m)+virial_thread(l, m)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE vdW_energy

! **************************************************************************************************
!> \brief This routine finds the non-local correlation contribution to the potential
!> (i.e. the derivative of the non-local piece of the energy with respect to
!> density) given in SOLER equation 10.  The u_alpha(k) functions were found
!> while calculating the energy. They are passed in as the matrix u_vdW.
!> Most of the required derivatives were calculated in the "get_q0_on_grid"
!> routine, but the derivative of the interpolation polynomials, P_alpha(q),
!> (SOLER equation 3) with respect to q is interpolated here, along with the
!> polynomials themselves.
!> \param q0 ...
!> \param dq0_drho ...
!> \param dq0_dgradrho ...
!> \param total_rho ...
!> \param u_vdW ...
!> \param potential ...
!> \param h_prefactor ...
!> \param dispersion_env ...
!> \param drho ...
!> \param dvol ...
!> \param virial ...
!> OpenMP added: Aug 2016  MTucker
! **************************************************************************************************
   SUBROUTINE get_potential(q0, dq0_drho, dq0_dgradrho, total_rho, u_vdW, potential, h_prefactor, &
                            dispersion_env, drho, dvol, virial)

      REAL(dp), DIMENSION(:), INTENT(in)                 :: q0, dq0_drho, dq0_dgradrho, total_rho
      REAL(dp), DIMENSION(:, :), INTENT(in)              :: u_vdW
      REAL(dp), DIMENSION(:), INTENT(out)                :: potential, h_prefactor
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      REAL(dp), DIMENSION(:, :), INTENT(in), OPTIONAL    :: drho
      REAL(dp), INTENT(IN), OPTIONAL                     :: dvol
      TYPE(virial_type), OPTIONAL, POINTER               :: virial

      CHARACTER(len=*), PARAMETER :: routineN = 'get_potential', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_grid, l, m, nl_type, nqs, P_i, &
                                                            q, q_hi, q_low
      LOGICAL                                            :: use_virial
      REAL(dp)                                           :: a, b, b_value, c, const, d, dP_dq0, dq, &
                                                            dq_6, e, f, P, prefactor, tmp_1_2, &
                                                            tmp_1_4, tmp_3_4
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: y
      REAL(dp), DIMENSION(:), POINTER                    :: q_mesh
      REAL(dp), DIMENSION(:, :), POINTER                 :: d2y_dx2
      REAL(dp), DIMENSION(3, 3)                          :: virial_thread

      CALL timeset(routineN, handle)

      IF (PRESENT(virial)) THEN
         use_virial = .TRUE.
         virial_thread(:, :) = 0.0_dp
         CPASSERT(PRESENT(drho))
         CPASSERT(PRESENT(dvol))
      ELSE
         use_virial = .FALSE.
      END IF

      b_value = dispersion_env%b_value
      const = 1.0_dp/(3.0_dp*b_value**(3.0_dp/2.0_dp)*pi**(5.0_dp/4.0_dp))
      potential = 0.0_dp
      h_prefactor = 0.0_dp

      d2y_dx2 => dispersion_env%d2y_dx2
      q_mesh => dispersion_env%q_mesh
      nqs = dispersion_env%nqs
      nl_type = dispersion_env%nl_type

!$OMP PARALLEL DEFAULT( NONE )                                         &
!$OMP           SHARED( nqs, u_vdW, q_mesh, q0, d2y_dx2, nl_type       &
!$OMP                 , potential, h_prefactor                         &
!$OMP                 , dq0_drho, dq0_dgradrho, total_rho, const       &
!$OMP                 , use_virial, drho, dvol, virial                 &
!$OMP                 )                                                &
!$OMP          PRIVATE( y                                              &
!$OMP                 , q_low, q_hi, q, dq, dq_6, A, b, c, d, e, f     &
!$OMP                 , P_i, dP_dq0, P, prefactor, l, m                &
!$OMP                 , tmp_1_2, tmp_1_4, tmp_3_4                      &
!$OMP                 )                                                &
!$OMP        REDUCTION(+: virial_thread                                &
!$OMP                 )

      ALLOCATE (y(nqs))

!$OMP DO
      DO i_grid = 1, SIZE(u_vdW, 1)
         q_low = 1
         q_hi = nqs
         ! Figure out which bin our value of q0 is in in the q_mesh
         DO WHILE ((q_hi-q_low) > 1)
            q = INT((q_hi+q_low)/2)
            IF (q_mesh(q) > q0(i_grid)) THEN
               q_hi = q
            ELSE
               q_low = q
            END IF
         END DO
         IF (q_hi == q_low) CPABORT("get_potential:  qhi == qlow")

         dq = q_mesh(q_hi)-q_mesh(q_low)
         dq_6 = dq/6.0_dp

         a = (q_mesh(q_hi)-q0(i_grid))/dq
         b = (q0(i_grid)-q_mesh(q_low))/dq
         c = (a**3-a)*dq*dq_6
         d = (b**3-b)*dq*dq_6
         e = (3.0_dp*a**2-1.0_dp)*dq_6
         f = (3.0_dp*b**2-1.0_dp)*dq_6

         DO P_i = 1, nqs
            y = 0.0_dp
            y(P_i) = 1.0_dp
            dP_dq0 = (y(q_hi)-y(q_low))/dq-e*d2y_dx2(P_i, q_low)+f*d2y_dx2(P_i, q_hi)
            P = a*y(q_low)+b*y(q_hi)+c*d2y_dx2(P_i, q_low)+d*d2y_dx2(P_i, q_hi)
            !! The first term in equation 13 of SOLER
            SELECT CASE (nl_type)
            CASE DEFAULT
               CPABORT("Unknown vdW-DF functional")
            CASE (vdw_nl_DRSLL, vdw_nl_LMKLL)
               potential(i_grid) = potential(i_grid)+u_vdW(i_grid, P_i)*(P+dP_dq0*dq0_drho(i_grid))
               prefactor = u_vdW(i_grid, P_i)*dP_dq0*dq0_dgradrho(i_grid)
            CASE (vdw_nl_RVV10)
               IF (total_rho(i_grid) > epsr) THEN
                  tmp_1_2 = SQRT(total_rho(i_grid))
                  tmp_1_4 = SQRT(tmp_1_2) ! == total_rho(i_grid)**(1.0_dp/4.0_dp)
                  tmp_3_4 = tmp_1_4*tmp_1_4*tmp_1_4 ! == total_rho(i_grid)**(3.0_dp/4.0_dp)
                  potential(i_grid) = potential(i_grid) &
                                      +u_vdW(i_grid, P_i)*(const*0.75_dp/tmp_1_4*P &
                                                           +const*tmp_3_4*dP_dq0*dq0_drho(i_grid))
                  prefactor = u_vdW(i_grid, P_i)*const*tmp_3_4*dP_dq0*dq0_dgradrho(i_grid)
               ELSE
                  prefactor = 0.0_dp
               ENDIF
            END SELECT
            IF (q0(i_grid) .NE. q_mesh(nqs)) THEN
               h_prefactor(i_grid) = h_prefactor(i_grid)+prefactor
            END IF

            IF (use_virial .AND. ABS(prefactor) > 0.0_dp) THEN
               SELECT CASE (nl_type)
               CASE DEFAULT
                  ! do nothing
               CASE (vdw_nl_RVV10)
                  prefactor = 0.5_dp*prefactor
               END SELECT
               prefactor = prefactor*dvol
               DO l = 1, 3
                  DO m = 1, l
                     virial_thread(l, m) = virial_thread(l, m)-prefactor*drho(i_grid, l)*drho(i_grid, m)
                  ENDDO
               ENDDO
            END IF

         END DO ! P_i = 1, nqs
      END DO ! i_grid = 1, SIZE(u_vdW, 1)
!$OMP END DO

      DEALLOCATE (y)
!$OMP END PARALLEL

      IF (use_virial) THEN
         DO l = 1, 3
            DO m = 1, (l-1)
               virial%pv_xc(l, m) = virial%pv_xc(l, m)+virial_thread(l, m)
               virial%pv_xc(m, l) = virial%pv_xc(l, m)
            ENDDO
            m = l
            virial%pv_xc(l, m) = virial%pv_xc(l, m)+virial_thread(l, m)
         ENDDO
      END IF

      CALL timestop(handle)
   END SUBROUTINE get_potential

   !! This routine first calculates the q value defined in (DION equations 11 and 12), then
   !! saturates it according to (SOLER equation 5).
! **************************************************************************************************
!> \brief ...
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

      CHARACTER(len=*), PARAMETER :: routineN = 'get_q0_on_grid_vdw', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: m_cut = 12
      REAL(dp), PARAMETER :: LDA_A = 0.031091_dp, LDA_a1 = 0.2137_dp, LDA_b1 = 7.5957_dp, &
         LDA_b2 = 3.5876_dp, LDA_b3 = 1.6382_dp, LDA_b4 = 0.49294_dp

      INTEGER                                            :: i_grid, index
      REAL(dp)                                           :: dq0_dq, exponent, gradient_correction, &
                                                            kF, LDA_1, LDA_2, q, q_cut, q_min, &
                                                            r_s, sqrt_r_s, Z_ab

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

      ! initialize q0-related arrays ...
      q0(:) = q_cut
      dq0_drho(:) = 0.0_dp
      dq0_dgradrho(:) = 0.0_dp

      DO i_grid = 1, SIZE(total_rho)

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
                               *(gradient_rho(i_grid, 1)**2+gradient_rho(i_grid, 2)**2+gradient_rho(i_grid, 3)**2)

         LDA_1 = 8.0_dp*pi/3.0_dp*(LDA_A*(1.0_dp+LDA_a1*r_s))
         LDA_2 = 2.0_dp*LDA_A*(LDA_b1*sqrt_r_s+LDA_b2*r_s+LDA_b3*r_s*sqrt_r_s+LDA_b4*r_s*r_s)
         !! ---------------------------------------------------------------
         !! This is the q value defined in equations 11 and 12 of DION
         !! ---------------------------------------------------------------
         q = kF+LDA_1*LOG(1.0_dp+1.0_dp/LDA_2)+gradient_correction
         !! ---------------------------------------------------------------
         !! Here, we calculate q0 by saturating q according to equation 5 of SOLER.  Also, we find
         !! the derivative dq0_dq needed for the derivatives dq0_drho and dq0_dgradrh0 discussed below.
         !! ---------------------------------------------------------------------------------------
         exponent = 0.0_dp
         dq0_dq = 0.0_dp
         DO index = 1, m_cut
            exponent = exponent+((q/q_cut)**index)/index
            dq0_dq = dq0_dq+((q/q_cut)**(index-1))
         END DO
         q0(i_grid) = q_cut*(1.0_dp-EXP(-exponent))
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

         dq0_drho(i_grid) = dq0_dq*(kF/3.0_dp-7.0_dp/3.0_dp*gradient_correction &
                                    -8.0_dp*pi/9.0_dp*LDA_A*LDA_a1*r_s*LOG(1.0_dp+1.0_dp/LDA_2) &
                                    +LDA_1/(LDA_2*(1.0_dp+LDA_2)) &
                                    *(2.0_dp*LDA_A*(LDA_b1/6.0_dp*sqrt_r_s+LDA_b2/3.0_dp*r_s+LDA_b3/2.0_dp*r_s*sqrt_r_s &
                                                    +2.0_dp*LDA_b4/3.0_dp*r_s**2)))

         dq0_dgradrho(i_grid) = total_rho(i_grid)*dq0_dq*2.0_dp*(-Z_ab)/(36.0_dp*kF*total_rho(i_grid)**2)

      END DO

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

      INTEGER                                            :: i_grid, index
      REAL(dp)                                           :: b_value, C_value, dk_dn, dq0_dq, dw0_dn, &
                                                            exponent, gmod2, k, mod_grad, q, &
                                                            q_cut, q_min, w0, wg2, wp2

      q_cut = dispersion_env%q_cut
      q_min = dispersion_env%q_min
      b_value = dispersion_env%b_value
      C_value = dispersion_env%c_value

      ! initialize q0-related arrays ...
      q0(:) = q_cut
      dq0_drho(:) = 0.0_dp
      dq0_dgradrho(:) = 0.0_dp

      DO i_grid = 1, SIZE(total_rho)

         gmod2 = gradient_rho(i_grid, 1)**2+gradient_rho(i_grid, 2)**2+gradient_rho(i_grid, 3)**2

         !if (total_rho(i_grid) > epsr .and. gmod2 > epsr) cycle
         IF (total_rho(i_grid) > epsr) THEN

            !! Calculate some intermediate values needed to find q
            !! ------------------------------------------------------------------------------------
            mod_grad = SQRT(gmod2)

            wp2 = 16.0_dp*pi*total_rho(i_grid)
            wg2 = 4_dp*C_value*(mod_grad/total_rho(i_grid))**4

            k = b_value*3.0_dp*pi*((total_rho(i_grid)/(9.0_dp*pi))**(1.0_dp/6.0_dp))
            w0 = SQRT(wg2+wp2/3.0_dp)

            q = w0/k

            !! Here, we calculate q0 by saturating q according
            !! ---------------------------------------------------------------------------------------

            exponent = 0.0_dp
            dq0_dq = 0.0_dp

            DO index = 1, m_cut

               exponent = exponent+((q/q_cut)**index)/index
               dq0_dq = dq0_dq+((q/q_cut)**(index-1))

            END DO

            q0(i_grid) = q_cut*(1.0_dp-EXP(-exponent))
            dq0_dq = dq0_dq*EXP(-exponent)

            !! ---------------------------------------------------------------------------------------

            IF (q0(i_grid) < q_min) THEN
               q0(i_grid) = q_min
            END IF

            !!---------------------------------Final values---------------------------------

            dw0_dn = 1.0_dp/(2.0_dp*w0)*(16.0_dp/3.0_dp*pi-4.0_dp*wg2/total_rho(i_grid))
            dk_dn = k/(6.0_dp*total_rho(i_grid))

            dq0_drho(i_grid) = dq0_dq*1.0_dp/(k**2.0)*(dw0_dn*k-dk_dn*w0)
            dq0_dgradrho(i_grid) = dq0_dq*1.0_dp/(2.0_dp*k*w0)*4.0_dp*wg2/(mod_grad**2)

         ENDIF

      END DO

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

      CHARACTER(len=*), PARAMETER :: routineN = 'get_q0_on_grid_eo_vdw', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: m_cut = 12
      REAL(dp), PARAMETER :: LDA_A = 0.031091_dp, LDA_a1 = 0.2137_dp, LDA_b1 = 7.5957_dp, &
         LDA_b2 = 3.5876_dp, LDA_b3 = 1.6382_dp, LDA_b4 = 0.49294_dp

      INTEGER                                            :: i_grid, index
      REAL(dp)                                           :: exponent, gradient_correction, kF, &
                                                            LDA_1, LDA_2, q, q_cut, q_min, r_s, &
                                                            sqrt_r_s, Z_ab

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

      ! initialize q0-related arrays ...
      q0(:) = q_cut
      DO i_grid = 1, SIZE(total_rho)
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
                               *(gradient_rho(i_grid, 1)**2+gradient_rho(i_grid, 2)**2+gradient_rho(i_grid, 3)**2)

         LDA_1 = 8.0_dp*pi/3.0_dp*(LDA_A*(1.0_dp+LDA_a1*r_s))
         LDA_2 = 2.0_dp*LDA_A*(LDA_b1*sqrt_r_s+LDA_b2*r_s+LDA_b3*r_s*sqrt_r_s+LDA_b4*r_s*r_s)
         !! ------------------------------------------------------------------------------------
         !! This is the q value defined in equations 11 and 12 of DION
         !! ---------------------------------------------------------------

         q = kF+LDA_1*LOG(1.0_dp+1.0_dp/LDA_2)+gradient_correction

         !! ---------------------------------------------------------------
         !! Here, we calculate q0 by saturating q according to equation 5 of SOLER.  Also, we find
         !! the derivative dq0_dq needed for the derivatives dq0_drho and dq0_dgradrh0 discussed below.
         !! ---------------------------------------------------------------------------------------
         exponent = 0.0_dp
         DO index = 1, m_cut
            exponent = exponent+((q/q_cut)**index)/index
         END DO
         q0(i_grid) = q_cut*(1.0_dp-EXP(-exponent))
         !! ---------------------------------------------------------------------------------------
         !! This is to handle a case with q0 too small.  We simply set it to the smallest q value in
         !! out q_mesh.  Hopefully this doesn't get used often (ever)
         !! ---------------------------------------------------------------------------------------
         IF (q0(i_grid) < q_min) THEN
            q0(i_grid) = q_min
         END IF
      END DO

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

      INTEGER                                            :: i_grid, index
      REAL(dp)                                           :: b_value, C_value, exponent, gmod2, k, &
                                                            mod_grad, q, q_cut, q_min, w0, wg2, wp2

      q_cut = dispersion_env%q_cut
      q_min = dispersion_env%q_min
      b_value = dispersion_env%b_value
      C_value = dispersion_env%c_value

      ! initialize q0-related arrays ...
      q0(:) = q_cut

      DO i_grid = 1, SIZE(total_rho)

         gmod2 = gradient_rho(i_grid, 1)**2+gradient_rho(i_grid, 2)**2+gradient_rho(i_grid, 3)**2

         !if (total_rho(i_grid) > epsr .and. gmod2 > epsr) cycle
         IF (total_rho(i_grid) > epsr) THEN

            !! Calculate some intermediate values needed to find q
            !! ------------------------------------------------------------------------------------
            mod_grad = SQRT(gmod2)

            wp2 = 16.0_dp*pi*total_rho(i_grid)
            wg2 = 4_dp*C_value*(mod_grad/total_rho(i_grid))**4

            k = b_value*3.0_dp*pi*((total_rho(i_grid)/(9.0_dp*pi))**(1.0_dp/6.0_dp))
            w0 = SQRT(wg2+wp2/3.0_dp)

            q = w0/k

            !! Here, we calculate q0 by saturating q according
            !! ---------------------------------------------------------------------------------------

            exponent = 0.0_dp

            DO index = 1, m_cut

               exponent = exponent+((q/q_cut)**index)/index

            END DO

            q0(i_grid) = q_cut*(1.0_dp-EXP(-exponent))

            IF (q0(i_grid) < q_min) THEN
               q0(i_grid) = q_min
            END IF

         ENDIF

      END DO

   END SUBROUTINE get_q0_on_grid_eo_rvv10

! **************************************************************************************************
!> \brief This routine is modeled after an algorithm from "Numerical Recipes in C" by Cambridge University
!> press, page 97.  It was adapted for Fortran, of course and for the problem at hand, in that it finds
!> the bin a particular x value is in and then loops over all the P_i functions so we only have to find
!> the bin once.
!> \param x ...
!> \param d2y_dx2 ...
!> \param evaluation_points ...
!> \param values ...
!> \par History
!>     OpenMP added: Aug 2016  MTucker
! **************************************************************************************************
   SUBROUTINE spline_interpolation(x, d2y_dx2, evaluation_points, values)
      REAL(dp), INTENT(in)                               :: x(:), d2y_dx2(:, :), evaluation_points(:)
      REAL(dp), INTENT(out)                              :: values(:, :)

      INTEGER                                            :: i_grid, index, lower_bound, &
                                                            Ngrid_points, Nx, P_i, upper_bound
      REAL(dp)                                           :: a, b, c, d, dx, const
      REAL(dp), ALLOCATABLE                              :: y(:)

      Nx = SIZE(x)
      Ngrid_points = SIZE(evaluation_points)

!$OMP PARALLEL DEFAULT( NONE )                                         &
!$OMP           SHARED( x, d2y_dx2, evaluation_points, values          &
!$OMP                 , Nx, Ngrid_points                               &
!$OMP                 )                                                &
!$OMP          PRIVATE( y, lower_bound, upper_bound, index             &
!$OMP                 , dx, const, a, b, c, d, P_i                     &
!$OMP                 )

      ALLOCATE (y(Nx))

!$OMP DO
      DO i_grid = 1, Ngrid_points
         lower_bound = 1
         upper_bound = Nx
         DO WHILE ((upper_bound-lower_bound) > 1)
            index = (upper_bound+lower_bound)/2
            IF (evaluation_points(i_grid) > x(index)) THEN
               lower_bound = index
            ELSE
               upper_bound = index
            END IF
         END DO

         dx = x(upper_bound)-x(lower_bound)
         const = dx*dx/6.0_dp

         a = (x(upper_bound)-evaluation_points(i_grid))/dx
         b = (evaluation_points(i_grid)-x(lower_bound))/dx
         c = (a**3-a)*const
         d = (b**3-b)*const

         DO P_i = 1, Nx
            y = 0
            y(P_i) = 1
            values(i_grid, P_i) = a*y(lower_bound) + b*y(upper_bound) &
                                + (c*d2y_dx2(P_i, lower_bound) + d*d2y_dx2(P_i, upper_bound))
         END DO
      END DO
!$OMP END DO

      DEALLOCATE (y)
!$OMP END PARALLEL
   END SUBROUTINE spline_interpolation

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
         !! In the Soler method, the polynomials that are interpolated are Kronecker delta funcions
         !! at a particular q point.  So, we set all y values to 0 except the one corresponding to
         !! the particular function P_i.
         !! ----------------------------------------------------------------------------------------
         y = 0.0_dp
         y(P_i) = 1.0_dp
         !! ----------------------------------------------------------------------------------------

         d2y_dx2(P_i, 1) = 0.0_dp
         temp_array(1) = 0.0_dp
         DO index = 2, Nx-1
            temp1 = (x(index)-x(index-1))/(x(index+1)-x(index-1))
            temp2 = temp1*d2y_dx2(P_i, index-1)+2.0_dp
            d2y_dx2(P_i, index) = (temp1-1.0_dp)/temp2
            temp_array(index) = (y(index+1)-y(index))/(x(index+1)-x(index)) &
                                -(y(index)-y(index-1))/(x(index)-x(index-1))
            temp_array(index) = (6.0_dp*temp_array(index)/(x(index+1)-x(index-1)) &
                                 -temp1*temp_array(index-1))/temp2
         END DO
         d2y_dx2(P_i, Nx) = 0.0_dp
         DO index = Nx-1, 1, -1
            d2y_dx2(P_i, index) = d2y_dx2(P_i, index)*d2y_dx2(P_i, index+1)+temp_array(index)
         END DO
      END DO
!$OMP END DO

      DEALLOCATE (temp_array, y)
!$OMP END PARALLEL

   END SUBROUTINE initialize_spline_interpolation

   ! **************************************************************************************************
   !! This routine is modeled after an algorithm from "Numerical Recipes in C" by Cambridge
   !! University Press, page 97.  Adapted for Fortran and the problem at hand.  This function is used to
   !! find the Phi_alpha_beta needed for equations 8 and 11 of SOLER.
! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param kernel_of_k ...
!> \param dispersion_env ...
!> \par History
!>     Optimised when adding OpenMP to vdw_energy (which calls this routine): Aug 2016 MTucker
! **************************************************************************************************
   SUBROUTINE interpolate_kernel(k, kernel_of_k, dispersion_env)
      REAL(dp), INTENT(in)                               :: k
      REAL(dp), INTENT(out)                              :: kernel_of_k(:, :)
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env

      CHARACTER(len=*), PARAMETER :: routineN = 'interpolate_kernel', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: k_i, Nr_points, q1_i, q2_i
      REAL(dp)                                           :: A, B, C, const, D, dk, r_max
      REAL(dp), DIMENSION(:, :, :), POINTER              :: d2phi_dk2, kernel

      Nr_points = dispersion_env%nr_points
      r_max = dispersion_env%r_max
      dk = dispersion_env%dk
      kernel => dispersion_env%kernel
      d2phi_dk2 => dispersion_env%d2phi_dk2

      !! Check to make sure that the kernel table we have is capable of dealing with this
      !! value of k.  If k is larger than Nr_points*2*pi/r_max then we can't perform the
      !! interpolation.  In that case, a kernel file should be generated with a larger number
      !! of radial points.
      !! -------------------------------------------------------------------------------------
      CPASSERT(k < Nr_points*dk)
      !! -------------------------------------------------------------------------------------
      kernel_of_k = 0.0_dp
      !! This integer division figures out which bin k is in since the kernel
      !! is set on a uniform grid.
      k_i = INT(k/dk)

      !! Test to see if we are trying to interpolate a k that is one of the actual
      !! function points we have.  The value is just the value of the function in that
      !! case.
      !! ----------------------------------------------------------------------------------------
      IF (MOD(k, dk) == 0) THEN
         DO q1_i = 1, dispersion_env%Nqs
            DO q2_i = 1, q1_i
               kernel_of_k(q1_i, q2_i) = kernel(k_i, q1_i, q2_i)
               kernel_of_k(q2_i, q1_i) = kernel(k_i, q2_i, q1_i)
            END DO
         END DO
         RETURN
      END IF
      !! ----------------------------------------------------------------------------------------
      !! If we are not on a function point then we carry out the interpolation
      !! ----------------------------------------------------------------------------------------
      const = dk*dk/6.0_dp
      A = (dk*(k_i+1.0_dp)-k)/dk
      B = (k-dk*k_i)/dk
      C = (A**3-A)*const
      D = (B**3-B)*const
      DO q1_i = 1, dispersion_env%Nqs
         DO q2_i = 1, q1_i
            kernel_of_k(q1_i, q2_i) = A*kernel(k_i, q1_i, q2_i)+B*kernel(k_i+1, q1_i, q2_i) &
                                      +(C*d2phi_dk2(k_i, q1_i, q2_i)+D*d2phi_dk2(k_i+1, q1_i, q2_i))
            kernel_of_k(q2_i, q1_i) = kernel_of_k(q1_i, q2_i)
         END DO
      END DO

   END SUBROUTINE interpolate_kernel

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param dkernel_of_dk ...
!> \param dispersion_env ...
!> \par History
!>     Optimised when adding OpenMP to vdw_energy (which calls this routine): Aug 2016 MTucker
! **************************************************************************************************
   SUBROUTINE interpolate_dkernel_dk(k, dkernel_of_dk, dispersion_env)
      REAL(dp), INTENT(in)                               :: k
      REAL(dp), INTENT(out)                              :: dkernel_of_dk(:, :)
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env

      CHARACTER(len=*), PARAMETER :: routineN = 'interpolate_dkernel_dk', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: k_i, Nr_points, q1_i, q2_i
      REAL(dp)                                           :: A, B, dAdk, dBdk, dCdk, dDdk, dk, dk_6
      REAL(dp), DIMENSION(:, :, :), POINTER              :: d2phi_dk2, kernel

      Nr_points = dispersion_env%nr_points
      dk = dispersion_env%dk
      kernel => dispersion_env%kernel
      d2phi_dk2 => dispersion_env%d2phi_dk2
      dk_6 = dk/6.0_dp

      CPASSERT(k < Nr_points*dk)

      dkernel_of_dk = 0.0_dp
      k_i = INT(k/dk)

      A = (dk*(k_i+1.0_dp)-k)/dk
      B = (k-dk*k_i)/dk
      dAdk = -1.0_dp/dk
      dBdk = 1.0_dp/dk
      dCdk = -(3*A**2-1.0_dp)*dk_6
      dDdk = (3*B**2-1.0_dp)*dk_6
      DO q1_i = 1, dispersion_env%Nqs
         DO q2_i = 1, q1_i
            dkernel_of_dk(q1_i, q2_i) = dAdk*kernel(k_i, q1_i, q2_i)+dBdk*kernel(k_i+1, q1_i, q2_i) &
                                        +dCdk*d2phi_dk2(k_i, q1_i, q2_i)+dDdk*d2phi_dk2(k_i+1, q1_i, q2_i)
            dkernel_of_dk(q2_i, q1_i) = dkernel_of_dk(q1_i, q2_i)
         END DO
      END DO

   END SUBROUTINE interpolate_dkernel_dk

! **************************************************************************************************

END MODULE qs_dispersion_nonloc
