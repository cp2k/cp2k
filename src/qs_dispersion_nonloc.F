!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of non local dispersion functionals
!>  Some routines adapted from:
!>  Copyright (C) 2001-2009 Quantum ESPRESSO group
!>  Copyright (C) 2009 Brian Kolb, Timo Thonhauser - Wake Forest University
!>  This file is distributed under the terms of the
!>  GNU General Public License. See the file `License'
!>  in the root directory of the present distribution,
!>  or http://www.gnu.org/copyleft/gpl.txt .
!> \author JGH
! *****************************************************************************
MODULE qs_dispersion_nonloc

  USE bibliography,                    ONLY: Dion2004,&
                                             Romanperez2009,&
                                             Sabatini2013,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: vdw_nl_DRSLL,&
                                             vdw_nl_LMKLL,&
                                             vdw_nl_RVV10,&
                                             xc_vdw_fun_nonloc
  USE input_section_types,             ONLY: section_vals_type
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
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  REAL(KIND=dp), PARAMETER :: epsr = 1.e-12_dp

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_nonloc'

  PUBLIC :: qs_dispersion_nonloc_init, calculate_dispersion_nonloc

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE qs_dispersion_nonloc_init(dispersion_env,nl_section,para_env,error)
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: nl_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_dispersion_nonloc_init', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: filename
    INTEGER                                  :: funit, handle, nqs, &
                                                nr_points, q1_i, q2_i, stat, &
                                                vdw_type
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    SELECT CASE (dispersion_env%nl_type)
       CASE DEFAULT
          STOP 'Unknown vdW-DF functional'
       CASE (vdw_nl_DRSLL,vdw_nl_LMKLL)
          CALL cite_reference(Dion2004)
       CASE (vdw_nl_RVV10)
          CALL cite_reference(Sabatini2013)
    END SELECT
    CALL cite_reference(RomanPerez2009)

    failure = .FALSE.

    vdw_type = dispersion_env%type
    SELECT CASE (vdw_type)
       CASE DEFAULT
          ! do nothing
       CASE (xc_vdw_fun_nonloc)
          ! setup information on non local functionals
          filename = dispersion_env%kernel_file_name
          IF(para_env%source==para_env%mepos) THEN
            ! Read the kernel information from file "filename"
            CALL open_file(file_name=filename,unit_number=funit,file_form="FORMATTED")
            READ(funit,*) nqs, nr_points
            READ(funit,*) dispersion_env%r_max
          END IF
          CALL mp_bcast(nqs,para_env%source,para_env%group)
          CALL mp_bcast(nr_points,para_env%source,para_env%group)
          CALL mp_bcast(dispersion_env%r_max,para_env%source,para_env%group)
          ALLOCATE(dispersion_env%q_mesh(nqs),dispersion_env%kernel(0:nr_points,nqs,nqs),&
                   dispersion_env%d2phi_dk2(0:nr_points,nqs,nqs),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          dispersion_env%nqs = nqs
          dispersion_env%nr_points = nr_points
          IF(para_env%source==para_env%mepos) THEN
            !! Read in the values of the q points used to generate this kernel
            READ(funit,"(1p4e23.14)") dispersion_env%q_mesh
            !! For each pair of q values, read in the function phi_q1_q2(k).
            !! That is, the fourier transformed kernel function assuming q1 and q2
            !! for all the values of r used.
            DO q1_i = 1, nqs
               DO q2_i = 1, q1_i
                  READ(funit,"(1p4e23.14)") dispersion_env%kernel(0:nr_points, q1_i, q2_i)
                  dispersion_env%kernel(0:nr_points, q2_i, q1_i) = dispersion_env%kernel(0:nr_points, q1_i, q2_i)
               END DO
            END DO
            !! Again, for each pair of q values (q1 and q2), read in the value
            !! of the second derivative of the above mentiond Fourier transformed
            !! kernel function phi_alpha_beta(k).  These are used for spline
            !! interpolation of the Fourier transformed kernel.
            DO q1_i = 1, nqs
               DO q2_i = 1, q1_i
                  READ(funit,"(1p4e23.14)") dispersion_env%d2phi_dk2(0:nr_points, q1_i, q2_i)
                  dispersion_env%d2phi_dk2(0:nr_points, q2_i, q1_i) = dispersion_env%d2phi_dk2(0:nr_points, q1_i, q2_i)
               END DO
            END DO
            CALL close_file(unit_number=funit)
          END IF
          CALL mp_bcast(dispersion_env%q_mesh,para_env%source,para_env%group)
          CALL mp_bcast(dispersion_env%kernel,para_env%source,para_env%group)
          CALL mp_bcast(dispersion_env%d2phi_dk2,para_env%source,para_env%group)
          ! 2nd derivates for interpolation
          ALLOCATE(dispersion_env%d2y_dx2(nqs,nqs),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL initialize_spline_interpolation (dispersion_env%q_mesh, dispersion_env%d2y_dx2)
          !
          dispersion_env%q_cut = dispersion_env%q_mesh(nqs)
          dispersion_env%q_min = dispersion_env%q_mesh(1)
          dispersion_env%dk = 2.0_dp*pi/dispersion_env%r_max

    END SELECT

    CALL timestop(handle)

  END SUBROUTINE qs_dispersion_nonloc_init

! *****************************************************************************

  SUBROUTINE calculate_dispersion_nonloc(qs_env,vxc_rho,rho_r,rho_g,energy_only,&
                                         cell,pw_pool,xc_pw_pool,virial,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: vxc_rho, rho_r, rho_g
    LOGICAL, INTENT(IN)                      :: energy_only
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool, xc_pw_pool
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_dispersion_nonloc', &
      routineP = moduleN//':'//routineN
    INTEGER, DIMENSION(3, 3), PARAMETER :: &
      nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))

    INTEGER                                  :: handle, i, i_grid, idir, &
                                                ispin, n1, n2, n3, np, nspin, &
                                                stat
    LOGICAL                                  :: failure, use_virial
    REAL(KIND=dp)                            :: b_value, beta, Ec_nl, sumnp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dq0_dgradrho, dq0_drho, hpot, &
                                                potential, q0, rho
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: drho, thetas
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_grid_type), POINTER              :: grid
    TYPE(pw_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: thetas_g
    TYPE(pw_p_type), ALLOCATABLE, &
      DIMENSION(:, :)                        :: drho_r
    TYPE(pw_type), POINTER                   :: tmp_g, tmp_r, vxc_g, vxc_r
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env
    TYPE(qs_energy_type), POINTER            :: energy

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(rho_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho_g),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)

    ! only spin nonpolarized
    nspin=1

    IF (PRESENT(virial)) THEN
       use_virial = virial%pv_calculate.AND.(.NOT.virial%pv_numer)
    ELSE
       use_virial = .FALSE.
    ENDIF
    IF(use_virial) THEN
       CPPrecondition(.NOT.energy_only,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. energy_only) THEN
       CPPrecondition(ASSOCIATED(vxc_rho),cp_failure_level,routineP,error,failure)
    END IF
    CALL get_qs_env(qs_env, dispersion_env=dispersion_env, energy=energy, &
                    para_env=para_env, error=error)

    b_value = dispersion_env%b_value
    beta = 0.03125_dp * (3.0_dp / (b_value**2.0_dp) )**0.75_dp

    ! tempory arrays for FFT
    CALL pw_pool_create_pw(pw_pool, tmp_g, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, error=error)
    CALL pw_pool_create_pw(pw_pool, tmp_r, use_data=REALDATA3D, in_space=REALSPACE, error=error)

    ! get density derivatives
    ALLOCATE ( drho_r(3,nspin), STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ispin=1,nspin
       DO idir=1,3
          NULLIFY(drho_r(idir,ispin)%pw)
          CALL pw_pool_create_pw(pw_pool,drho_r(idir,ispin)%pw, &
               use_data=REALDATA3D, in_space=REALSPACE, error=error)
          CALL pw_transfer ( rho_g(ispin)%pw, tmp_g, error=error)
          CALL pw_derive ( tmp_g, nd(:,idir) ,error=error)
          CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw, error=error)
       END DO
    END DO

    ispin = 1
    CPPrecondition(rho_r(ispin)%pw%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
    CALL pw_transfer ( rho_g(ispin)%pw, tmp_g, error=error)
    CALL pw_transfer ( tmp_g, tmp_r, error=error)
    np = SIZE(tmp_r%cr3d)
    ALLOCATE ( rho(np), drho(np,3), STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    rho(:) = RESHAPE(tmp_r%cr3d,(/np/))
    drho(:,1) = RESHAPE(drho_r(1,ispin)%pw%cr3d,(/np/))
    drho(:,2) = RESHAPE(drho_r(2,ispin)%pw%cr3d,(/np/))
    drho(:,3) = RESHAPE(drho_r(3,ispin)%pw%cr3d,(/np/))
    !! ---------------------------------------------------------------------------------
    !! Find the value of q0 for all assigned grid points.  q is defined in equations
    !! 11 and 12 of DION and q0 is the saturated version of q defined in equation
    !! 5 of SOLER.  This routine also returns the derivatives of the q0s with respect
    !! to the charge-density and the gradient of the charge-density.  These are needed
    !! for the potential calculated below.
    !! ---------------------------------------------------------------------------------

    IF (energy_only) THEN
       ALLOCATE ( q0(np), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       SELECT CASE (dispersion_env%nl_type)
          CASE DEFAULT
             STOP 'Unknown vdW-DF functional'
          CASE (vdw_nl_DRSLL,vdw_nl_LMKLL)
             CALL get_q0_on_grid_eo_vdw (rho,drho,q0,dispersion_env)
          CASE (vdw_nl_RVV10)
             CALL get_q0_on_grid_eo_rvv10 (rho,drho,q0,dispersion_env)
       END SELECT
    ELSE
       ALLOCATE ( q0(np), dq0_drho(np), dq0_dgradrho(np), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       SELECT CASE (dispersion_env%nl_type)
          CASE DEFAULT
             STOP 'Unknown vdW-DF functional'
          CASE (vdw_nl_DRSLL,vdw_nl_LMKLL)
             CALL get_q0_on_grid_vdw (rho,drho,q0,dq0_drho,dq0_dgradrho,dispersion_env)
          CASE (vdw_nl_RVV10)
             CALL get_q0_on_grid_rvv10 (rho,drho,q0,dq0_drho,dq0_dgradrho,dispersion_env)
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
    ALLOCATE ( thetas(np,dispersion_env%nqs), STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !! Interpolate the P_i polynomials defined in equation 3 in SOLER for the particular
    !! q0 values we have.
    CALL spline_interpolation(dispersion_env%q_mesh, dispersion_env%d2y_dx2, q0, thetas)

    !! Form the thetas where theta is defined as rho*p_i(q0) for vdW and vdW2 or
    !! constant*rho^(3/4)*p_i(q0) for rVV10
    !! ------------------------------------------------------------------------------------
    SELECT CASE (dispersion_env%nl_type)
       CASE DEFAULT
          STOP 'Unknown vdW-DF functional'
       CASE (vdw_nl_DRSLL,vdw_nl_LMKLL)
          DO i=1,dispersion_env%nqs
             thetas(:,i) = thetas(:,i) * rho(:)
          END DO
       CASE (vdw_nl_RVV10)
          DO i_grid=1, np
             IF (rho(i_grid) > epsr) THEN 
                DO i=1,dispersion_env%nqs
                   thetas(i_grid,i) = thetas(i_grid,i) * (1.0_dp / (3.0_dp * SQRT(pi) * &
                   ( b_value**(3.0_dp/2.0_dp) ) ) ) * (rho(i_grid) / pi)**(3.0_dp/4.0_dp)
                END DO
             ELSE
                thetas(i_grid,:) = 0.0_dp
             END IF
          END DO
    END SELECT
    !! ------------------------------------------------------------------------------------
    !! Get thetas in reciprocal space.
    n1 = SIZE(tmp_r%cr3d,1)
    n2 = SIZE(tmp_r%cr3d,2)
    n3 = SIZE(tmp_r%cr3d,3)
    ALLOCATE ( thetas_g(dispersion_env%nqs), STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i=1,dispersion_env%nqs
       tmp_r%cr3d = RESHAPE(thetas(:,i),(/n1,n2,n3/))
       NULLIFY(thetas_g(i)%pw)
       CALL pw_pool_create_pw(pw_pool, thetas_g(i)%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, error=error)
       CALL pw_transfer ( tmp_r, thetas_g(i)%pw ,error=error)
    END DO
    grid => thetas_g(1)%pw%pw_grid
    !! ---------------------------------------------------------------------------------------------
    !! Carry out the integration in equation 8 of SOLER.  This also turns the thetas array into the 
    !! precursor to the u_i(k) array which is inverse fourier transformed to get the u_i(r) functions
    !! of SOLER equation 11.  Add the energy we find to the output variable etxc.
    !! --------------------------------------------------------------------------------------------------
    sumnp = np
    CALL mp_sum(sumnp,para_env%group)
    IF (use_virial) THEN
       ! calculates kernel contribution to stress
       CALL vdW_energy(thetas_g, dispersion_env, Ec_nl, energy_only, virial, error=error)
       SELECT CASE (dispersion_env%nl_type)
          CASE (vdw_nl_RVV10)
             Ec_nl = 0.5_dp*Ec_nl + beta*SUM(rho(:))*grid%vol/sumnp
       END SELECT
       ! calculates energy contribution to stress
       ! potential contribution to stress is calculated together with other potentials (Hxc)
       DO idir=1,3
          virial%pv_xc(idir,idir) = virial%pv_xc(idir,idir) + Ec_nl
       END DO
    ELSE
       CALL vdW_energy(thetas_g, dispersion_env, Ec_nl, energy_only, error=error)
       SELECT CASE (dispersion_env%nl_type)
          CASE (vdw_nl_RVV10)
             Ec_nl = 0.5_dp*Ec_nl + beta*SUM(rho(:))*grid%vol/sumnp
       END SELECT
    END IF
    CALL mp_sum(Ec_nl,para_env%group)
    energy%dispersion = Ec_nl

    IF (.NOT.energy_only) THEN
       !! ----------------------------------------------------------------------------
       !! Inverse Fourier transform the u_i(k) to get the u_i(r) of SOLER equation 11.
       !!-----------------------------------------------------------------------------
       DO i=1,dispersion_env%nqs
          CALL pw_transfer ( thetas_g(i)%pw, tmp_r ,error=error)
          thetas(:,i) = RESHAPE(tmp_r%cr3d,(/np/))
       END DO
       !! -------------------------------------------------------------------------
       !! Here we allocate the array to hold the potential. This is calculated via 
       !! equation 10 of SOLER, using the u_i(r) calculated from equations 11 and 
       !! 12 of SOLER.  Each processor allocates the array to be the size of the 
       !! full grid  because, as can be seen in SOLER equation 9, processors need
       !! to access grid points outside their allocated regions.
       !! -------------------------------------------------------------------------
       ALLOCATE ( potential(np), hpot(np), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF(use_virial) THEN
          ! calculates gradient contribution to stress
          grid => tmp_g%pw_grid
          CALL get_potential(q0, dq0_drho, dq0_dgradrho, rho, thetas, potential, hpot, &
                             dispersion_env, drho, grid%dvol, virial, error=error)
       ELSE
          CALL get_potential(q0, dq0_drho, dq0_dgradrho, rho, thetas, potential, hpot, &
                             dispersion_env, error=error)
       END IF
       SELECT CASE (dispersion_env%nl_type)
          CASE (vdw_nl_RVV10)
             potential = 0.5_dp*potential + beta
             hpot = 0.5_dp*hpot
       END SELECT
       CALL pw_pool_create_pw(pw_pool, vxc_r, use_data=REALDATA3D, in_space=REALSPACE, error=error)
       vxc_r%cr3d = RESHAPE(potential,(/n1,n2,n3/))
       DO idir=1,3
          tmp_r%cr3d = RESHAPE(hpot,(/n1,n2,n3/)) * drho_r(idir,ispin)%pw%cr3d
          CALL pw_transfer ( tmp_r, tmp_g ,error=error)
          CALL pw_derive ( tmp_g, nd(:,idir) ,error=error)
          CALL pw_transfer ( tmp_g, tmp_r ,error=error)
          CALL pw_axpy(tmp_r,vxc_r,-1._dp,error)
       END DO
       CALL pw_transfer ( vxc_r, tmp_g ,error=error)
       CALL pw_pool_give_back_pw(pw_pool, vxc_r, error=error)
       CALL pw_pool_create_pw(xc_pw_pool, vxc_r, use_data=REALDATA3D, in_space=REALSPACE, error=error)
       CALL pw_pool_create_pw(xc_pw_pool, vxc_g, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, error=error)
       CALL pw_transfer ( tmp_g, vxc_g ,error=error)
       CALL pw_transfer ( vxc_g, vxc_r ,error=error)
       CALL pw_axpy(vxc_r,vxc_rho(ispin)%pw,1._dp,error)
       CALL pw_pool_give_back_pw(xc_pw_pool, vxc_r, error=error)
       CALL pw_pool_give_back_pw(xc_pw_pool, vxc_g, error=error)
    END IF

    DEALLOCATE ( thetas, STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO i=1,dispersion_env%nqs
       CALL pw_pool_give_back_pw(pw_pool, thetas_g(i)%pw, error=error)
    END DO
    DO ispin=1,nspin
       DO idir=1,3
          CALL pw_pool_give_back_pw(pw_pool, drho_r(idir,ispin)%pw, error=error)
       END DO
    END DO
    CALL pw_pool_give_back_pw(pw_pool, tmp_r, error=error)
    CALL pw_pool_give_back_pw(pw_pool, tmp_g, error=error)

    DEALLOCATE ( rho, drho, drho_r, thetas_g, STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (energy_only) THEN
       DEALLOCATE ( q0, STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
       DEALLOCATE ( q0, dq0_drho, dq0_dgradrho, STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE calculate_dispersion_nonloc

!! #############################################################################################
!! This routine carries out the integration of equation 8 of SOLER.  It returns the non-local 
!! exchange-correlation energy and the u_alpha(k) arrays used to find the u_alpha(r) arrays via 
!! equations 11 and 12 in SOLER.
!! energy contribution to stress is added in qs_force

  SUBROUTINE vdW_energy(thetas_g, dispersion_env, vdW_xc_energy, energy_only, virial, error)
    TYPE(pw_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: thetas_g
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env
    REAL(KIND=dp), INTENT(OUT)               :: vdW_xc_energy
    LOGICAL, INTENT(IN)                      :: energy_only
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'vdW_energy', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: uu
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: theta, theta_g, thetam
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: u_vdw(:,:)
    INTEGER                                  :: handle, ig, igs, iq, l, m, &
                                                nqs, q1_i, q2_i, stat
    LOGICAL                                  :: failure, use_virial
    REAL(KIND=dp)                            :: g, g2, g2_last, g_multiplier, &
                                                gm
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: dkernel_of_dk, kernel_of_k
    TYPE(pw_grid_type), POINTER              :: grid

    failure = .FALSE.
    CALL timeset(routineN,handle)
    nqs = dispersion_env%nqs
    ALLOCATE(kernel_of_k(nqs,nqs),dkernel_of_dk(nqs,nqs),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(theta(nqs),thetam(nqs),theta_g(nqs),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(virial)) THEN
       use_virial = .TRUE.
    ELSE
       use_virial = .FALSE.
    END IF

    vdW_xc_energy = 0._dp
    g2_last = HUGE(0._dp)
    grid => thetas_g(1)%pw%pw_grid

    IF (.NOT.energy_only) THEN
       ALLOCATE(u_vdW(grid%ngpts_cut_local,nqs),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       u_vdW(:,:) = CMPLX(0.0_dp,0.0_dp,KIND=dp)
    END IF

    IF (grid%grid_span == HALFSPACE) THEN
       g_multiplier = 2._dp
    ELSE
       g_multiplier = 1._dp
    END IF
    igs = 0
    DO ig = 1, grid%ngpts_cut_local
       g2 = grid%gsq(ig)
       IF (ABS(g2 - g2_last) > 1.e-10) THEN
          g2_last = g2
          igs = igs + 1
          g = SQRT(g2)
          CALL interpolate_kernel(g,kernel_of_k,dispersion_env,error)
          IF (use_virial) CALL interpolate_dkernel_dk(g,dkernel_of_dk,dispersion_env,error)
       END IF
       DO iq=1,nqs
          theta(iq) = thetas_g(iq)%pw%cc(ig)
       END DO
       DO q2_i = 1, Nqs
          uu = CMPLX(0.0_dp,0.0_dp,KIND=dp)
          DO q1_i = 1, Nqs
             uu  = uu + kernel_of_k(q2_i,q1_i)*theta(q1_i)
          END DO
          IF (ig < grid%first_gne0) THEN
            vdW_xc_energy = vdW_xc_energy + (uu*CONJG(theta(q2_i)))
          ELSE
            vdW_xc_energy = vdW_xc_energy + G_multiplier * (uu*CONJG(theta(q2_i)))
          END IF
          IF (.NOT.energy_only) u_vdW(ig,q2_i)  = uu
          IF (use_virial .AND. ig >=grid%first_gne0) THEN
             DO q1_i = 1, Nqs
                 gm = 0.5_dp * G_multiplier * grid%vol*theta(q1_i)*dkernel_of_dk(q1_i,q2_i)*CONJG(theta(q2_i))
                 SELECT CASE (dispersion_env%nl_type)
                    CASE (vdw_nl_RVV10)
                       gm = 0.5_dp*gm
                 END SELECT
                 DO l = 1, 3
                     DO m = 1, l
                        virial%pv_xc(l,m) = virial%pv_xc(l,m) - gm * (grid%g(l,ig)*grid%g(m,ig)) / g
                        virial%pv_xc(m,l) = virial%pv_xc(l,m)
                     END DO
                 END DO
             END DO
          END IF
       END DO
    END DO
    vdW_xc_energy = vdW_xc_energy  * grid%vol * 0.5_dp

    IF (.NOT.energy_only) THEN
       DO iq=1,nqs
          DO ig = 1, grid%ngpts_cut_local
             thetas_g(iq)%pw%cc(ig) = u_vdW(ig,iq)
          END DO
       END DO
       DEALLOCATE(u_vdW,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(theta,thetam,theta_g,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(kernel_of_k,dkernel_of_dk,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE vdW_energy

! *****************************************************************************
!! This routine finds the non-local correlation contribution to the potential 
!! (i.e. the derivative of the non-local piece of the energy with respect to 
!! density) given in SOLER equation 10.  The u_alpha(k) functions were found
!! while calculating the energy. They are passed in as the matrix u_vdW.  
!! Most of the required derivatives were calculated in the "get_q0_on_grid" 
!! routine, but the derivative of the interpolation polynomials, P_alpha(q),
!! (SOLER equation 3) with respect to q is interpolated here, along with the 
!! polynomials themselves.
  SUBROUTINE get_potential(q0, dq0_drho, dq0_dgradrho, total_rho, u_vdW, potential, h_prefactor, &
                         dispersion_env, drho, dvol, virial, error)

    REAL(dp), DIMENSION(:), INTENT(in)       :: q0, dq0_drho, dq0_dgradrho, &
                                                total_rho
    REAL(dp), DIMENSION(:, :), INTENT(in)    :: u_vdW
    REAL(dp), DIMENSION(:), INTENT(inout)    :: potential, h_prefactor
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env
    REAL(dp), DIMENSION(:, :), OPTIONAL      :: drho
    REAL(dp), INTENT(IN), OPTIONAL           :: dvol
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_grid, l, m, nqs, &
                                                P_i, q, q_hi, q_low, stat
    LOGICAL                                  :: failure, use_virial
    REAL(dp)                                 :: a, b, b_value, c, const, d, &
                                                dP_dq0, dq, e, f, P, prefactor
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: y
    REAL(dp), DIMENSION(:), POINTER          :: q_mesh
    REAL(dp), DIMENSION(:, :), POINTER       :: d2y_dx2

    failure = .FALSE.
    CALL timeset(routineN,handle)

    IF (PRESENT(virial)) THEN
       use_virial = .TRUE.
       CPPrecondition(PRESENT(drho),cp_failure_level,routineP,error,failure)
       CPPrecondition(PRESENT(dvol),cp_failure_level,routineP,error,failure)
    ELSE
       use_virial = .FALSE.
    END IF
  
    b_value = dispersion_env%b_value
    const = 1.0_dp / (3.0_dp * b_value**(3.0_dp/2.0_dp) * pi**(5.0_dp/4.0_dp) )
    potential = 0.0_dp
    h_prefactor   = 0.0_dp

    d2y_dx2 => dispersion_env%d2y_dx2
    q_mesh => dispersion_env%q_mesh
    nqs = dispersion_env%nqs
    ALLOCATE(y(nqs),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO i_grid = 1,SIZE(u_vdw,1)
       q_low = 1
       q_hi = Nqs
       ! Figure out which bin our value of q0 is in in the q_mesh
       DO WHILE ( (q_hi - q_low) > 1)
          q = INT((q_hi + q_low)/2)
          IF (q_mesh(q) > q0(i_grid)) THEN
             q_hi = q
          ELSE
             q_low = q
          END IF
       END DO
       IF (q_hi == q_low) STOP 'get_potential:  qhi == qlow'

       dq = q_mesh(q_hi) - q_mesh(q_low)

       a = (q_mesh(q_hi) - q0(i_grid))/dq
       b = (q0(i_grid) - q_mesh(q_low))/dq
       c = (a**3 - a)*dq**2/6.0_dp
       d = (b**3 - b)*dq**2/6.0_dp
       e = (3.0_dp*a**2 - 1.0_dp)*dq/6.0_dp
       f = (3.0_dp*b**2 - 1.0_dp)*dq/6.0_dp

       DO P_i = 1, Nqs
          y = 0.0_dp
          y(P_i) = 1.0_dp
          dP_dq0 = (y(q_hi) - y(q_low))/dq - e*d2y_dx2(P_i,q_low) + f*d2y_dx2(P_i,q_hi)
          P = a*y(q_low) + b*y(q_hi) + c*d2y_dx2(P_i,q_low) + d*d2y_dx2(P_i,q_hi)
          !! The first term in equation 13 of SOLER
          SELECT CASE (dispersion_env%nl_type)
             CASE DEFAULT
                STOP 'Unknown vdW-DF functional'
             CASE (vdw_nl_DRSLL,vdw_nl_LMKLL)
                potential(i_grid) = potential(i_grid) + u_vdW(i_grid,P_i)* (P + dP_dq0 * dq0_drho(i_grid))
                prefactor = u_vdW(i_grid,P_i) * dP_dq0 * dq0_dgradrho(i_grid)
             CASE (vdw_nl_RVV10)
                IF (total_rho(i_grid) > epsr) THEN
                   potential(i_grid) = potential(i_grid) + u_vdW(i_grid,P_i)* &
                      (const * (3.0_dp/4.0_dp) / (total_rho(i_grid)**(1.0_dp/4.0_dp))  * P + &
                      const * total_rho(i_grid)**(3.0_dp/4.0_dp) *  dP_dq0 * dq0_drho(i_grid))
                   prefactor = u_vdW(i_grid,P_i) * const * total_rho(i_grid)**(3.0_dp/4.0_dp) * &
                      dP_dq0 * dq0_dgradrho(i_grid)
                ELSE
                   prefactor = 0.0_dp
                ENDIF
          END SELECT
          IF (q0(i_grid) .NE. q_mesh(Nqs)) THEN
             h_prefactor(i_grid) = h_prefactor(i_grid) + prefactor
          END IF
          IF (use_virial .AND. ABS(prefactor) > 0.0_dp) THEN
             SELECT CASE (dispersion_env%nl_type)
                CASE DEFAULT
                   ! do nothing
                CASE (vdw_nl_RVV10)
                   prefactor = 0.5_dp*prefactor
             END SELECT
             prefactor = prefactor*dvol
             DO l = 1, 3
                DO m = 1, l
                   virial%pv_xc(l,m) = virial%pv_xc(l,m) - prefactor*drho(i_grid,l)*drho(i_grid,m)
                   virial%pv_xc(m,l) = virial%pv_xc(l,m)
                ENDDO
             ENDDO
          END IF
       END DO
    END DO
    DEALLOCATE(y,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE get_potential

  !! This routine first calculates the q value defined in (DION equations 11 and 12), then
  !! saturates it according to (SOLER equation 5).  
  SUBROUTINE get_q0_on_grid_vdw (total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho, dispersion_env)
  !!
  !! more specifically it calculates the following
  !!
  !!     q0(ir) = q0 as defined above
  !!     dq0_drho(ir) = total_rho * d q0 /d rho 
  !!     dq0_dgradrho = total_rho / |gradient_rho| * d q0 / d |gradient_rho|
  !!
    REAL(dp), INTENT(IN)                     :: total_rho(:), &
                                                gradient_rho(:,:)
    REAL(dp), INTENT(OUT)                    :: q0(:), dq0_drho(:), &
                                                dq0_dgradrho(:)
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env

    INTEGER, PARAMETER :: m_cut = 12                                
    REAL(dp), PARAMETER :: LDA_A = 0.031091_dp, LDA_a1 = 0.2137_dp, &
      LDA_b1 = 7.5957_dp, LDA_b2 = 3.5876_dp    , LDA_b3 = 1.6382_dp, &
      LDA_b4 = 0.49294_dp

    INTEGER                                  :: i_grid, index
    REAL(dp)                                 :: dq0_dq, exponent, &
                                                gradient_correction, kF, &
                                                LDA_1, LDA_2, q, q_cut, &
                                                q_min, r_s, sqrt_r_s, Z_ab

    q_cut = dispersion_env%q_cut
    q_min = dispersion_env%q_min
    SELECT CASE (dispersion_env%nl_type)
       CASE DEFAULT
          STOP 'Unknown vdW-DF functional'
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
           * (gradient_rho(i_grid,1)**2+gradient_rho(i_grid,2)**2+gradient_rho(i_grid,3)**2)
 
       LDA_1 =  8.0_dp*pi/3.0_dp*(LDA_A*(1.0_dp+LDA_a1*r_s))
       LDA_2 =  2.0_dp*LDA_A * (LDA_b1*sqrt_r_s + LDA_b2*r_s + LDA_b3*r_s*sqrt_r_s + LDA_b4*r_s*r_s)
       !! ---------------------------------------------------------------
       !! This is the q value defined in equations 11 and 12 of DION
       !! ---------------------------------------------------------------
       q = kF + LDA_1 * LOG(1.0_dp+1.0_dp/LDA_2) + gradient_correction
       !! ---------------------------------------------------------------
       !! Here, we calculate q0 by saturating q according to equation 5 of SOLER.  Also, we find 
       !! the derivative dq0_dq needed for the derivatives dq0_drho and dq0_dgradrh0 discussed below.
       !! ---------------------------------------------------------------------------------------
       exponent = 0.0_dp
       dq0_dq = 0.0_dp
       DO index = 1, m_cut
          exponent = exponent + ( (q/q_cut)**index)/index
          dq0_dq = dq0_dq + ( (q/q_cut)**(index-1))
       END DO
       q0(i_grid) = q_cut*(1.0_dp - EXP(-exponent))
       dq0_dq = dq0_dq * EXP(-exponent)
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
  
       dq0_drho(i_grid) = dq0_dq * (kF/3.0_dp - 7.0_dp/3.0_dp*gradient_correction &
            - 8.0_dp*pi/9.0_dp * LDA_A*LDA_a1*r_s*LOG(1.0_dp+1.0_dp/LDA_2) &
            + LDA_1/(LDA_2*(1.0_dp + LDA_2)) &
            * (2.0_dp*LDA_A*(LDA_b1/6.0_dp*sqrt_r_s + LDA_b2/3.0_dp*r_s + LDA_b3/2.0_dp*r_s*sqrt_r_s &
            + 2.0_dp*LDA_b4/3.0_dp*r_s**2)))
  
       dq0_dgradrho(i_grid) = total_rho(i_grid) * dq0_dq * 2.0_dp * (-Z_ab)/(36.0_dp*kF*total_rho(i_grid)**2)
  
    END DO
 
  END SUBROUTINE get_q0_on_grid_vdw

  SUBROUTINE get_q0_on_grid_rvv10 (total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho, dispersion_env)
  !!
  !! more specifically it calculates the following
  !!
  !!     q0(ir) = q0 as defined above
  !!     dq0_drho(ir) = total_rho * d q0 /d rho 
  !!     dq0_dgradrho = total_rho / |gradient_rho| * d q0 / d |gradient_rho|
  !!
    REAL(dp), INTENT(IN)                     :: total_rho(:), &
                                                gradient_rho(:,:)
    REAL(dp), INTENT(OUT)                    :: q0(:), dq0_drho(:), &
                                                dq0_dgradrho(:)
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env

    INTEGER, PARAMETER                       :: m_cut = 12

    INTEGER                                  :: i_grid, index
    REAL(dp)                                 :: b_value, C_value, dk_dn, &
                                                dq0_dq, dw0_dn, exponent, &
                                                gmod2, k, mod_grad, q, q_cut, &
                                                q_min, w0, wg2, wp2

    q_cut = dispersion_env%q_cut
    q_min = dispersion_env%q_min
    b_value = dispersion_env%b_value
    C_value = dispersion_env%c_value

    ! initialize q0-related arrays ... 
    q0(:) = q_cut
    dq0_drho(:) = 0.0_dp
    dq0_dgradrho(:) = 0.0_dp
  
    DO i_grid = 1, SIZE(total_rho)
  
     gmod2 = gradient_rho(i_grid,1)**2+gradient_rho(i_grid,2)**2+gradient_rho(i_grid,3)**2
 
     !if (total_rho(i_grid) > epsr .and. gmod2 > epsr) cycle
     IF (total_rho(i_grid) > epsr) THEN
     
       !! Calculate some intermediate values needed to find q
       !! ------------------------------------------------------------------------------------
       mod_grad = SQRT(gmod2)

       wp2= 16.0_dp*pi*total_rho(i_grid)
       wg2 = 4_dp*C_value * (mod_grad/total_rho(i_grid))**4

       k = b_value*3.0_dp* pi* ((total_rho(i_grid)/(9.0_dp*pi))**(1.0_dp/6.0_dp))
       w0 = SQRT( wg2 + wp2/3.0_dp  )

       q = w0 / k 
     
       !! Here, we calculate q0 by saturating q according 
       !! ---------------------------------------------------------------------------------------

       exponent = 0.0_dp
       dq0_dq = 0.0_dp
     
       DO index = 1, m_cut
        
          exponent = exponent + ( (q/q_cut)**index)/index
          dq0_dq = dq0_dq + ( (q/q_cut)**(index-1))
        
       END DO

       q0(i_grid) = q_cut*(1.0_dp - EXP(-exponent))
       dq0_dq = dq0_dq * EXP(-exponent)

       !! ---------------------------------------------------------------------------------------

       IF (q0(i_grid) < q_min) THEN
         q0(i_grid) = q_min
       END IF

       !!---------------------------------Final values---------------------------------  

       dw0_dn = 1.0_dp / (2.0_dp * w0 ) * ( 16.0_dp/3.0_dp * pi  - 4.0_dp * wg2 / total_rho(i_grid) )                
       dk_dn = k / ( 6.0_dp * total_rho(i_grid) )

       dq0_drho(i_grid) = dq0_dq * 1.0_dp / (k**2.0) * (dw0_dn * k - dk_dn * w0 )   
       dq0_dgradrho(i_grid) = dq0_dq * 1.0_dp / ( 2.0_dp * k * w0 ) * 4.0_dp * wg2 / (mod_grad**2)

     ENDIF  

    END DO
 
  END SUBROUTINE get_q0_on_grid_rvv10

  SUBROUTINE get_q0_on_grid_eo_vdw (total_rho, gradient_rho, q0, dispersion_env)

    REAL(dp), INTENT(IN)                     :: total_rho(:), &
                                                gradient_rho(:,:)
    REAL(dp), INTENT(OUT)                    :: q0(:)
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env

    INTEGER, PARAMETER :: m_cut = 12                                
    REAL(dp), PARAMETER :: LDA_A = 0.031091_dp, LDA_a1 = 0.2137_dp, &
      LDA_b1 = 7.5957_dp, LDA_b2 = 3.5876_dp    , LDA_b3 = 1.6382_dp, &
      LDA_b4 = 0.49294_dp

    INTEGER                                  :: i_grid, index
    REAL(dp)                                 :: exponent, &
                                                gradient_correction, kF, &
                                                LDA_1, LDA_2, q, q_cut, &
                                                q_min, r_s, sqrt_r_s, Z_ab

    q_cut = dispersion_env%q_cut
    q_min = dispersion_env%q_min
    SELECT CASE (dispersion_env%nl_type)
       CASE DEFAULT
          STOP 'Unknown vdW-DF functional'
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
           * (gradient_rho(i_grid,1)**2+gradient_rho(i_grid,2)**2+gradient_rho(i_grid,3)**2)
 
       LDA_1 =  8.0_dp*pi/3.0_dp*(LDA_A*(1.0_dp+LDA_a1*r_s))
       LDA_2 =  2.0_dp*LDA_A * (LDA_b1*sqrt_r_s + LDA_b2*r_s + LDA_b3*r_s*sqrt_r_s + LDA_b4*r_s*r_s)
       !! ------------------------------------------------------------------------------------
       !! This is the q value defined in equations 11 and 12 of DION
       !! ---------------------------------------------------------------
 
       q = kF + LDA_1 * LOG(1.0_dp+1.0_dp/LDA_2) + gradient_correction
 
       !! ---------------------------------------------------------------
       !! Here, we calculate q0 by saturating q according to equation 5 of SOLER.  Also, we find 
       !! the derivative dq0_dq needed for the derivatives dq0_drho and dq0_dgradrh0 discussed below.
       !! ---------------------------------------------------------------------------------------
       exponent = 0.0_dp
       DO index = 1, m_cut
          exponent = exponent + ( (q/q_cut)**index)/index
       END DO
       q0(i_grid) = q_cut*(1.0_dp - EXP(-exponent))
       !! ---------------------------------------------------------------------------------------
       !! This is to handle a case with q0 too small.  We simply set it to the smallest q value in
       !! out q_mesh.  Hopefully this doesn't get used often (ever)
       !! ---------------------------------------------------------------------------------------
       IF (q0(i_grid) < q_min) THEN
          q0(i_grid) = q_min
       END IF
    END DO
 
  END SUBROUTINE get_q0_on_grid_eo_vdw


  SUBROUTINE get_q0_on_grid_eo_rvv10 (total_rho, gradient_rho, q0, dispersion_env)

    REAL(dp), INTENT(IN)                     :: total_rho(:), &
                                                gradient_rho(:,:)
    REAL(dp), INTENT(OUT)                    :: q0(:)
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env

    INTEGER, PARAMETER :: m_cut = 12                                

    INTEGER                                  :: i_grid, index
    REAL(dp)                                 :: b_value, C_value, exponent, &
                                                gmod2, k, mod_grad, q, q_cut, &
                                                q_min, w0, wg2, wp2

    q_cut = dispersion_env%q_cut
    q_min = dispersion_env%q_min
    b_value = dispersion_env%b_value
    C_value = dispersion_env%c_value

    ! initialize q0-related arrays ... 
    q0(:) = q_cut
  
    DO i_grid = 1, SIZE(total_rho)
  
     gmod2 = gradient_rho(i_grid,1)**2+gradient_rho(i_grid,2)**2+gradient_rho(i_grid,3)**2
 
     !if (total_rho(i_grid) > epsr .and. gmod2 > epsr) cycle
     IF (total_rho(i_grid) > epsr) THEN
     
       !! Calculate some intermediate values needed to find q
       !! ------------------------------------------------------------------------------------
       mod_grad = SQRT(gmod2)

       wp2= 16.0_dp*pi*total_rho(i_grid)
       wg2 = 4_dp*C_value * (mod_grad/total_rho(i_grid))**4

       k = b_value*3.0_dp* pi* ((total_rho(i_grid)/(9.0_dp*pi))**(1.0_dp/6.0_dp))
       w0 = SQRT( wg2 + wp2/3.0_dp  )

       q = w0 / k 
     
       !! Here, we calculate q0 by saturating q according 
       !! ---------------------------------------------------------------------------------------

       exponent = 0.0_dp
     
       DO index = 1, m_cut
        
          exponent = exponent + ( (q/q_cut)**index)/index
        
       END DO

       q0(i_grid) = q_cut*(1.0_dp - EXP(-exponent))

       IF (q0(i_grid) < q_min) THEN
         q0(i_grid) = q_min
       END IF

     ENDIF  

    END DO
 
  END SUBROUTINE get_q0_on_grid_eo_rvv10

  !! This routine is modeled after an algorithm from "Numerical Recipes in C" by Cambridge University
  !! press, page 97.  It was adapted for Fortran, of course and for the problem at hand, in that it finds
  !! the bin a particular x value is in and then loops over all the P_i functions so we only have to find
  !! the bin once.
  SUBROUTINE spline_interpolation (x, d2y_dx2, evaluation_points, values)

    REAL(dp), INTENT(in)                     :: x(:), d2y_dx2(:,:), &
                                                evaluation_points(:)
    REAL(dp), INTENT(inout)                  :: values(:,:)

    INTEGER                                  :: i_grid, index, lower_bound, &
                                                Ngrid_points, Nx, P_i, &
                                                upper_bound
    REAL(dp)                                 :: a, b, c, d, dx
    REAL(dp), ALLOCATABLE                    :: y(:)

    Nx = SIZE(x)
    Ngrid_points = SIZE(evaluation_points)

    !! Allocate the temporary array
    ALLOCATE( y(Nx) )

    DO i_grid=1, Ngrid_points
       lower_bound = 1
       upper_bound = Nx
       DO WHILE ( (upper_bound - lower_bound) > 1 )
          index = (upper_bound+lower_bound)/2
          IF ( evaluation_points(i_grid) > x(index) ) THEN
             lower_bound = index
          ELSE
             upper_bound = index
          END IF
       END DO

       dx = x(upper_bound)-x(lower_bound)

       a = (x(upper_bound) - evaluation_points(i_grid))/dx
       b = (evaluation_points(i_grid) - x(lower_bound))/dx
       c = ((a**3-a)*dx**2)/6.0_dp
       d = ((b**3-b)*dx**2)/6.0_dp

       DO P_i = 1, Nx
          y = 0
          y(P_i) = 1
          values(i_grid, P_i) = a*y(lower_bound) + b*y(upper_bound) &
             + (c*d2y_dx2(P_i,lower_bound) + d*d2y_dx2(P_i, upper_bound))
       END DO
    END DO

    DEALLOCATE( y )

  END SUBROUTINE spline_interpolation

  !! This routine is modeled after an algorithm from "Numerical Recipes in C" by Cambridge
  !! University Press, pages 96-97.  It was adapted for Fortran and for the problem at hand.
  SUBROUTINE initialize_spline_interpolation (x, d2y_dx2)

    REAL(dp), INTENT(in)                     :: x(:)
    REAL(dp), INTENT(inout)                  :: d2y_dx2(:,:)

    INTEGER                                  :: index, Nx, P_i
    REAL(dp)                                 :: temp1, temp2
    REAL(dp), ALLOCATABLE                    :: temp_array(:), y(:)

    Nx = SIZE(x)
    ALLOCATE( temp_array(Nx), y(Nx) )

    DO P_i=1, Nx
       !! In the Soler method, the polynomicals that are interpolated are Kronecker delta funcions
       !! at a particular q point.  So, we set all y values to 0 except the one corresponding to 
       !! the particular function P_i.
       !! ----------------------------------------------------------------------------------------
       y = 0.0_dp
       y(P_i) = 1.0_dp
       !! ----------------------------------------------------------------------------------------
  
       d2y_dx2(P_i,1) = 0.0_dp
       temp_array(1) = 0.0_dp
       DO index = 2, Nx-1
          temp1 = (x(index)-x(index-1))/(x(index+1)-x(index-1))
          temp2 = temp1 * d2y_dx2(P_i,index-1) + 2.0_dp
          d2y_dx2(P_i,index) = (temp1-1.0_dp)/temp2
          temp_array(index) = (y(index+1)-y(index))/(x(index+1)-x(index)) &
               - (y(index)-y(index-1))/(x(index)-x(index-1))
          temp_array(index) = (6.0_dp*temp_array(index)/(x(index+1)-x(index-1)) &
               - temp1*temp_array(index-1))/temp2
       END DO
       d2y_dx2(P_i,Nx) = 0.0_dp
       DO index=Nx-1, 1, -1
          d2y_dx2(P_i,index) = d2y_dx2(P_i,index) * d2y_dx2(P_i,index+1) + temp_array(index)
       END DO
    END DO

    DEALLOCATE( temp_array, y)

  END SUBROUTINE initialize_spline_interpolation

  ! *****************************************************************************
  !! This routine is modeled after an algorithm from "Numerical Recipes in C" by Cambridge
  !! University Press, page 97.  Adapted for Fortran and the problem at hand.  This function is used to 
  !! find the Phi_alpha_beta needed for equations 8 and 11 of SOLER.
  SUBROUTINE interpolate_kernel(k, kernel_of_k, dispersion_env,error)

    REAL(dp), INTENT(in)                     :: k
    REAL(dp), INTENT(inout)                  :: kernel_of_k(:,:)
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'interpolate_kernel', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k_i, Nr_points, q1_i, q2_i
    LOGICAL                                  :: failure
    REAL(dp)                                 :: A, B, C, D, dk, r_max
    REAL(dp), DIMENSION(:, :, :), POINTER    :: d2phi_dk2, kernel

    failure = .FALSE.
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
    CPPrecondition(k < Nr_points*dk,cp_failure_level,routineP,error,failure)
    !! -------------------------------------------------------------------------------------
    kernel_of_k = 0.0_dp
    !! This integer division figures out which bin k is in since the kernel
    !! is set on a uniform grid.
    k_i = INT(k/dk)

    !! Test to see if we are trying to interpolate a k that is one of the actual
    !! function points we have.  The value is just the value of the function in that
    !! case.
    !! ----------------------------------------------------------------------------------------
    IF (MOD(k,dk) == 0) THEN
       DO q1_i = 1, dispersion_env%Nqs
          DO q2_i = 1, q1_i
             kernel_of_k(q1_i, q2_i) = kernel(k_i,q1_i, q2_i)
             kernel_of_k(q2_i, q1_i) = kernel(k_i,q2_i, q1_i)
          END DO
       END DO
       RETURN
    END IF
    !! ----------------------------------------------------------------------------------------
    !! If we are not on a function point then we carry out the interpolation
    !! ----------------------------------------------------------------------------------------
    A = (dk*(k_i+1.0_dp) - k)/dk
    B = (k - dk*k_i)/dk
    C = (A**3-A)*dk**2/6.0_dp
    D = (B**3-B)*dk**2/6.0_dp
    DO q1_i = 1, dispersion_env%Nqs
       DO q2_i = 1, q1_i
          kernel_of_k(q1_i, q2_i) = A*kernel(k_i, q1_i, q2_i) + B*kernel(k_i+1, q1_i, q2_i) &
               +(C*d2phi_dk2(k_i, q1_i, q2_i) + D*d2phi_dk2(k_i+1, q1_i, q2_i))
          kernel_of_k(q2_i, q1_i) = kernel_of_k(q1_i, q2_i)
       END DO
    END DO

  END SUBROUTINE interpolate_kernel

  SUBROUTINE interpolate_dkernel_dk(k, dkernel_of_dk, dispersion_env, error)
    REAL(dp), INTENT(in)                     :: k
    REAL(dp), INTENT(inout)                  :: dkernel_of_dk(:,:)
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'interpolate_dkernel_dk', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k_i, Nr_points, q1_i, q2_i
    LOGICAL                                  :: failure
    REAL(dp)                                 :: A, B, dAdk, dBdk, dCdk, dDdk, &
                                                dk, r_max
    REAL(dp), DIMENSION(:, :, :), POINTER    :: d2phi_dk2, kernel

    failure = .FALSE.
    Nr_points = dispersion_env%nr_points
    r_max = dispersion_env%r_max
    dk = dispersion_env%dk
    kernel => dispersion_env%kernel
    d2phi_dk2 => dispersion_env%d2phi_dk2

    CPPrecondition(k < Nr_points*dk,cp_failure_level,routineP,error,failure)

    dkernel_of_dk = 0.0_dp
    k_i = INT(k/dk)
 
    A = (dk*(k_i+1.0_dp) - k)/dk
    B = (k - dk*k_i)/dk
    dAdk = -1.0_dp/dk
    dBdk = 1.0_dp/dk
    dCdk = -((3*A**2 -1.0_dp)/6.0_dp)*dk
    dDdk = ((3*B**2 -1.0_dp)/6.0_dp)*dk
    DO q1_i = 1, dispersion_env%Nqs
       DO q2_i = 1, q1_i
          dkernel_of_dk(q1_i, q2_i) = dAdk*kernel(k_i, q1_i, q2_i) + dBdk*kernel(k_i+1, q1_i, q2_i) &
               + dCdk*d2phi_dk2(k_i, q1_i, q2_i) + dDdk*d2phi_dk2(k_i+1, q1_i, q2_i)
          dkernel_of_dk(q2_i, q1_i) = dkernel_of_dk(q1_i, q2_i)
       END DO
    END DO

  END SUBROUTINE interpolate_Dkernel_Dk

! *****************************************************************************

END  MODULE qs_dispersion_nonloc
