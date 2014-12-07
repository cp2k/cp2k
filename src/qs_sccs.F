!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Self-consistent continuum solvation (SCCS) model implementation
!> \author  Matthias Krack (MK)
!> \version 1.0
!> \par Literature:
!>      - J.-L. Fattebert and F. Gygi,
!>        Density functional theory for efficient ab initio molecular dynamics
!>        simulations in solution, J. Comput. Chem. 23, 662-666 (2002)
!>      - O. Andreussi, I. Dabo, and N. Marzari,
!>        Revised self-consistent continuum solvation in electronic-structure
!>        calculations, J. Chem. Phys. 136, 064102-20 (2012)
!> \par History:
!>      - Creation (10.10.2013,MK)
!>      - Derivatives using finite differences (26.11.2013,MK)
!>      - Cube file dump of the dielectric function (19.12.2013,MK)
!>      - Cube file dump of the polarisation potential (20.12.2013,MK)
!>      - Calculation of volume and surface of the cavity (21.12.2013,MK)
!>      - Functional derivative of the cavitation energy (28.12.2013,MK)
! *****************************************************************************

MODULE qs_sccs

  USE cp_control_types,                ONLY: dft_control_type,&
                                             sccs_control_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr,&
                                             low_print_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
  USE cp_realspace_grid_init,          ONLY: init_input_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE input_constants,                 ONLY: sccs_andreussi,&
                                             sccs_derivative_cd3,&
                                             sccs_derivative_cd5,&
                                             sccs_derivative_cd7,&
                                             sccs_derivative_fft,&
                                             sccs_fattebert_gygi
  USE input_section_types,             ONLY: section_get_ivals,&
                                             section_get_lval,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8
  USE mathconstants,                   ONLY: twopi
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sum
  USE particle_list_types,             ONLY: particle_list_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_derive,&
                                             pw_integral_ab,&
                                             pw_scale,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE realspace_grid_types,            ONLY: &
       pw2rs, realspace_grid_desc_type, realspace_grid_input_type, &
       realspace_grid_p_type, realspace_grid_type, rs2pw, rs_grid_create, &
       rs_grid_create_descriptor, rs_grid_release, &
       rs_grid_release_descriptor, rs_grid_zero, rs_pw_transfer
  USE termination,                     ONLY: print_warning,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_sccs'

  PUBLIC :: sccs

CONTAINS

! *****************************************************************************
!> \brief   Self-consistent continuum solvation (SCCS) model implementation
!> \param qs_env ...
!> \param rho_tot_gspace ...
!> \param v_hartree_gspace ...
!> \param v_sccs ...
!> \param h_stress ...
!> \param error ...
!> \par History:
!>      - Creation (10.10.2013,MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************

  SUBROUTINE sccs(qs_env,rho_tot_gspace,v_hartree_gspace,v_sccs,h_stress,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_type), POINTER                   :: rho_tot_gspace, &
                                                v_hartree_gspace, v_sccs
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT), OPTIONAL                  :: h_stress
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'sccs', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: tol = 1.0E-12_dp

    CHARACTER(LEN=4*default_string_length)   :: message
    CHARACTER(LEN=default_string_length)     :: cube_path, filename, &
                                                my_pos_cube, print_path
    INTEGER                                  :: cube_unit, di, dj, handle, i, &
                                                ispin, iter, j, k, nspin, &
                                                output_unit, print_level
    INTEGER(KIND=int_8)                      :: ngpts
    INTEGER, DIMENSION(3)                    :: lb, ub
    LOGICAL                                  :: append_cube, &
                                                calculate_stress_tensor, &
                                                failure, should_output
    REAL(KIND=dp) :: cavity_surface, cavity_volume, cell_volume, dphi2, dvol, &
      eps0, f, norm_drho2, polarisation_charge, rho_delta, rho_delta_avg, &
      rho_delta_max, rho_iter_new, tot_rho_elec, tot_rho_solute
    REAL(KIND=dp), DIMENSION(3)              :: abc
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: cp_subsys
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(3)            :: dln_eps_elec, dphi_tot, &
                                                drho_elec
    TYPE(pw_p_type), DIMENSION(3, 3)         :: d2rho_elec
    TYPE(pw_p_type), DIMENSION(5)            :: work_r3d
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_r
    TYPE(pw_p_type), POINTER                 :: rho_r_sccs
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_type), POINTER :: deps_elec, eps_elec, ln_eps_elec, &
      norm_drho_elec, phi_pol, phi_solute, rho_elec, rho_iter_old, &
      rho_solute, rho_tot, rho_tot_zero, theta
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(sccs_control_type), POINTER         :: sccs_control
    TYPE(section_vals_type), POINTER         :: input

    failure = .FALSE.

    CALL timeset(routineN,handle)

    NULLIFY (auxbas_pw_pool)
    NULLIFY (cp_subsys)
    NULLIFY (deps_elec)
    NULLIFY (dft_control)
    NULLIFY (energy)
    NULLIFY (eps_elec)
    NULLIFY (input)
    NULLIFY (ln_eps_elec)
    NULLIFY (logger)
    NULLIFY (norm_drho_elec)
    NULLIFY (para_env)
    NULLIFY (particles)
    NULLIFY (phi_pol)
    NULLIFY (phi_solute)
    NULLIFY (poisson_env)
    NULLIFY (pw_env)
    NULLIFY (pw_pools)
    NULLIFY (rho)
    NULLIFY (rho_elec)
    NULLIFY (rho_iter_old)
    NULLIFY (rho_solute)
    NULLIFY (rho_tot)
    NULLIFY (rho_tot_zero)
    NULLIFY (sccs_control)
    NULLIFY (scf_env)
    NULLIFY (theta)

    ! Load data from Quickstep environment
    CALL get_qs_env(qs_env=qs_env,&
                    cp_subsys=cp_subsys,&
                    dft_control=dft_control,&
                    energy=energy,&
                    input=input,&
                    para_env=para_env,&
                    pw_env=pw_env,&
                    rho=rho,&
                    scf_env=scf_env,&
                    error=error)
    CALL cp_subsys_get(cp_subsys,particles=particles,error=error)

    sccs_control => dft_control%sccs_control

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho_tot_gspace),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(v_hartree_gspace),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(v_sccs),cp_failure_level,routineP,error,failure)

    IF (PRESENT(h_stress)) THEN
      calculate_stress_tensor = .TRUE.
      h_stress(:,:) = 0.0_dp
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "No stress tensor is implemented for SCCS models yet",para_env)
    ELSE
       calculate_stress_tensor = .FALSE.
    END IF

    ! Get access to the PW grid pool
    CALL pw_env_get(pw_env,&
                    auxbas_pw_pool=auxbas_pw_pool,&
                    pw_pools=pw_pools,&
                    poisson_env=poisson_env,&
                    error=error)

    CALL pw_zero(v_sccs,error=error)

    ! Calculate no SCCS contribution, if the requested SCF convergence threshold is not reached yet
    IF (.NOT.sccs_control%sccs_activated) THEN
      IF (sccs_control%eps_scf > 0.0_dp) THEN
        IF ((scf_env%iter_delta > sccs_control%eps_scf).OR.&
            ((qs_env%scf_env%outer_scf%iter_count == 0).AND.&
             (qs_env%scf_env%iter_count <= 1))) THEN
          CALL pw_poisson_solve(poisson_env=poisson_env,&
                                density=rho_tot_gspace,&
                                ehartree=energy%hartree,&
                                vhartree=v_hartree_gspace,&
                                error=error)
          energy%sccs_pol = 0.0_dp
          energy%sccs_cav = 0.0_dp
          energy%sccs_dis = 0.0_dp
          energy%sccs_rep = 0.0_dp
          CALL timestop(handle)
          RETURN
        END IF
      END IF
      sccs_control%sccs_activated = .TRUE.
    END IF

    nspin = dft_control%nspins

    ! Manage print output control
    logger => cp_error_get_logger(error)
    print_level = logger%iter_info%print_level
    print_path = "DFT%PRINT%SCCS"
    should_output = (BTEST(cp_print_key_should_output(logger%iter_info,input,&
                           TRIM(print_path),error=error),cp_p_file))
    output_unit = cp_print_key_unit_nr(logger,input,TRIM(print_path),&
                                       extension=".sccs",&
                                       ignore_should_output=should_output,&
                                       log_filename=.FALSE.,&
                                       error=error)

    ! Get work storage for the 3d grids in r-space
    DO i=1,SIZE(work_r3d)
      NULLIFY (work_r3d(i)%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool,&
                             work_r3d(i)%pw,&
                             use_data=REALDATA3D,&
                             in_space=REALSPACE,&
                             error=error)
    END DO

    ! Assign total electronic density in r-space
    rho_elec => work_r3d(1)%pw

    ! Retrieve grid parameters
    ngpts = rho_elec%pw_grid%ngpts
    dvol = rho_elec%pw_grid%dvol
    cell_volume = rho_elec%pw_grid%vol
    abc(1:3) = REAL(rho_elec%pw_grid%npts(1:3),KIND=dp)*rho_elec%pw_grid%dr(1:3)
    lb(1:3) = rho_elec%pw_grid%bounds_local(1,1:3)
    ub(1:3) = rho_elec%pw_grid%bounds_local(2,1:3)

    ! Get rho
    CALL qs_rho_get(rho_struct=rho,&
                    rho_r=rho_r,&
                    rho_r_sccs=rho_r_sccs,&
                    error=error)

    ! Retrieve the last rho_iter from the previous SCCS cycle if available
    CPPrecondition(ASSOCIATED(rho_r_sccs),cp_failure_level,routineP,error,failure)
    rho_iter_old => rho_r_sccs%pw

    ! Retrieve the total electronic density in r-space
    CALL pw_copy(rho_r(1)%pw,rho_elec,error=error)
    DO ispin=2,nspin
      CALL pw_axpy(rho_r(ispin)%pw,rho_elec,error=error)
    END DO
    tot_rho_elec = accurate_sum(rho_elec%cr3d)*dvol
    CALL mp_sum(tot_rho_elec,para_env%group)

    ! Calculate the dielectric (smoothed) function of rho_elec in r-space
    eps_elec => work_r3d(2)%pw
    deps_elec => work_r3d(3)%pw
    eps0 = sccs_control%epsilon_solvent
    SELECT CASE (sccs_control%method_id)
    CASE (sccs_andreussi)
      CALL andreussi(rho_elec,eps_elec,deps_elec,eps0,sccs_control%rho_max,&
                     sccs_control%rho_min)
    CASE (sccs_fattebert_gygi)
      CALL fattebert_gygi(rho_elec,eps_elec,deps_elec,eps0,sccs_control%beta,&
                          sccs_control%rho_zero)
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Invalid method specified for SCCS model",&
                        para_env)
    END SELECT

    ! Optional output of the dielectric function in cube file format
    filename = "DIELECTRIC_FUNCTION"
    cube_path = TRIM(print_path)//"%"//TRIM(filename)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,input,TRIM(cube_path),error=error),&
              cp_p_file)) THEN
      append_cube = section_get_lval(input,TRIM(cube_path)//"%APPEND",error=error)
      my_pos_cube="REWIND"
      IF (append_cube) my_pos_cube="APPEND"
      cube_unit = cp_print_key_unit_nr(logger,input,TRIM(cube_path),&
                                       extension=".cube",middle_name=TRIM(filename),&
                                       file_position=my_pos_cube,log_filename=.FALSE.,&
                                       error=error)
      IF (output_unit > 0) THEN
        INQUIRE (UNIT=cube_unit,NAME=filename)
        WRITE (UNIT=output_unit,FMT="(/,T2,A,/,/,T2,A)")&
         "The dielectric function is written in cube file format to the file:",TRIM(filename)
      END IF
      CALL cp_pw_to_cube(eps_elec,cube_unit,TRIM(filename),particles=particles,&
                         stride=section_get_ivals(input,TRIM(cube_path)//"%STRIDE",error=error),&
                         error=error)
      CALL cp_print_key_finished_output(cube_unit,logger,input,TRIM(cube_path),error=error)
    END IF

    ! Calculate the (quantum) volume and surface of the solute cavity
    cavity_surface = 0.0_dp
    cavity_volume = 0.0_dp

    IF (should_output.AND.(ABS(eps0 - 1.0_dp) > tol)) THEN

      ! Initialise the switching function theta
      theta => work_r3d(4)%pw
      CALL pw_zero(theta,error=error)

      ! Calculate the (quantum) volume of the solute cavity
      f = 1.0_dp/(eps0 - 1.0_dp)
      !$omp parallel do default(none) &
      !$omp             private(i,j,k) &
      !$omp             shared(eps0,eps_elec,f,lb,theta,ub)
      DO k=lb(3),ub(3)
        DO j=lb(2),ub(2)
          DO i=lb(1),ub(1)
            theta%cr3d(i,j,k) = f*(eps0 - eps_elec%cr3d(i,j,k))
          END DO
        END DO
      END DO
      !$omp end parallel do
      cavity_volume = accurate_sum(theta%cr3d)*dvol
      CALL mp_sum(cavity_volume,para_env%group)

      ! Calculate the derivative of the electronic density in r-space
      ! TODO: Could be retrieved from the QS environment
      DO i=1,3
        NULLIFY (drho_elec(i)%pw)
        CALL pw_pool_create_pw(auxbas_pw_pool,&
                               drho_elec(i)%pw,&
                               use_data=REALDATA3D,&
                               in_space=REALSPACE,&
                               error=error)
      END DO
      CALL derive(rho_elec,drho_elec,sccs_derivative_fft,pw_env,input,para_env,error)

      ! Calculate the norm of the gradient of the electronic density in r-space
      norm_drho_elec => work_r3d(5)%pw
      !$omp parallel do default(none) &
      !$omp             private(i,j,k) &
      !$omp             shared(drho_elec,lb,norm_drho_elec,ub)
      DO k=lb(3),ub(3)
        DO j=lb(2),ub(2)
          DO i=lb(1),ub(1)
            norm_drho_elec%cr3d(i,j,k) = SQRT(drho_elec(1)%pw%cr3d(i,j,k)*&
                                              drho_elec(1)%pw%cr3d(i,j,k) +&
                                              drho_elec(2)%pw%cr3d(i,j,k)*&
                                              drho_elec(2)%pw%cr3d(i,j,k) +&
                                              drho_elec(3)%pw%cr3d(i,j,k)*&
                                              drho_elec(3)%pw%cr3d(i,j,k))
          END DO
        END DO
      END DO
      !$omp end parallel do

      ! Optional output of the norm of the density gradient in cube file format
      filename = "DENSITY_GRADIENT"
      cube_path = TRIM(print_path)//"%"//TRIM(filename)
      IF (BTEST(cp_print_key_should_output(logger%iter_info,input,TRIM(cube_path),error=error),&
                cp_p_file)) THEN
        append_cube = section_get_lval(input,TRIM(cube_path)//"%APPEND",error=error)
        my_pos_cube="REWIND"
        IF (append_cube) my_pos_cube="APPEND"
        cube_unit = cp_print_key_unit_nr(logger,input,TRIM(cube_path),&
                                         extension=".cube",middle_name=TRIM(filename),&
                                         file_position=my_pos_cube,log_filename=.FALSE.,&
                                         error=error)
        IF (output_unit > 0) THEN
          INQUIRE (UNIT=cube_unit,NAME=filename)
          WRITE (UNIT=output_unit,FMT="(/,T2,A,/,/,T2,A)")&
           "The norm of the density gradient is written in cube file format to the file:",TRIM(filename)
        END IF
        CALL cp_pw_to_cube(norm_drho_elec,cube_unit,TRIM(filename),particles=particles,&
                           stride=section_get_ivals(input,TRIM(cube_path)//"%STRIDE",error=error),&
                           error=error)
        CALL cp_print_key_finished_output(cube_unit,logger,input,TRIM(cube_path),error=error)
      END IF

      ! Calculate the (quantum) surface of the solute cavity
      SELECT CASE (sccs_control%method_id)
      CASE (sccs_andreussi)
        CALL surface_andreussi(rho_elec,norm_drho_elec,drho_elec,theta,eps0,&
                               sccs_control%rho_max,sccs_control%rho_min,&
                               sccs_control%delta_rho)
      CASE (sccs_fattebert_gygi)
        CALL surface_fattebert_gygi(rho_elec,norm_drho_elec,drho_elec,theta,eps0,&
                                    sccs_control%beta,sccs_control%rho_zero,&
                                    sccs_control%delta_rho)
      CASE DEFAULT
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "Invalid method specified for SCCS model",&
                          para_env)
      END SELECT
      cavity_surface = accurate_sum(theta%cr3d)*dvol
      CALL mp_sum(cavity_surface,para_env%group)

      IF (sccs_control%gamma_solvent > 0.0_dp) THEN

        CALL print_warning(routineN,moduleN,__LINE__,&
                           "Experimental feature requested: The contribution from the "//&
                           "cavitation energy is added as an extra potential term to the "//&
                           "Kohn-Sham potential",para_env)

        ! Calculate the second derivatives of the electronic density in r-space
        DO di=1,3
          DO dj=1,3
            NULLIFY (d2rho_elec(dj,di)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool,&
                                   d2rho_elec(dj,di)%pw,&
                                   use_data=REALDATA3D,&
                                   in_space=REALSPACE,&
                                   error=error)
          END DO
          CALL derive(drho_elec(di)%pw,d2rho_elec(:,di),sccs_derivative_fft,pw_env,&
                      input,para_env,error)
        END DO

        ! Calculate the contribution of the cavitation energy to the Kohn-Sham potential
        !$omp parallel do default(none) &
        !$omp             private(di,dj,i,j,k,norm_drho2) &
        !$omp             shared(d2rho_elec,drho_elec,lb,norm_drho_elec,sccs_control) &
        !$omp             shared(theta,ub,v_sccs)
        DO k=lb(3),ub(3)
          DO j=lb(2),ub(2)
            DO i=lb(1),ub(1)
              norm_drho2 = norm_drho_elec%cr3d(i,j,k)*norm_drho_elec%cr3d(i,j,k)
              DO di=1,3
                DO dj=1,3
                  v_sccs%cr3d(i,j,k) = sccs_control%gamma_solvent*theta%cr3d(i,j,k)*&
                                       (drho_elec(di)%pw%cr3d(i,j,k)*&
                                        drho_elec(dj)%pw%cr3d(i,j,k)*&
                                        d2rho_elec(di,dj)%pw%cr3d(i,j,k)/norm_drho2 -&
                                        d2rho_elec(di,di)%pw%cr3d(i,j,k))/norm_drho2
                END DO
              END DO
            END DO
          END DO
        END DO
        !$omp end parallel do
        CALL pw_scale(v_sccs,dvol,error=error)

      END IF

      ! Release storage
      NULLIFY (theta)
      DO i=1,3
        CALL pw_pool_give_back_pw(auxbas_pw_pool,drho_elec(i)%pw,error=error)
      END DO

    END IF

    ! Retrieve the total charge density (core + elec) of the solute in r-space
    rho_solute => work_r3d(4)%pw
    CALL pw_zero(rho_solute,error=error)
    CALL pw_transfer(rho_tot_gspace,rho_solute,error=error)
    tot_rho_solute = accurate_sum(rho_solute%cr3d)*dvol
    CALL mp_sum(tot_rho_solute,para_env%group)

    ! Reassign work storage to rho_tot_zero, because rho_elec is no longer needed
    rho_tot_zero => rho_elec
    NULLIFY (rho_elec)

    ! Build the initial (rho_iter = 0) total charge density (solute plus polarisation) in r-space
    ! eps_elec <- ln(eps_elec)
    !$omp parallel do default(none) &
    !$omp             private(i,j,k) &
    !$omp             shared(eps_elec,lb,message,output_unit,para_env,ub) &
    !$omp             shared(rho_solute,rho_tot_zero)
    DO k=lb(3),ub(3)
      DO j=lb(2),ub(2)
        DO i=lb(1),ub(1)
          IF (eps_elec%cr3d(i,j,k) < 1.0E-12_dp) THEN
            WRITE (UNIT=message,FMT="(A,ES12.3,A,3(I0,A))")&
             "SCCS| Invalid dielectric function value ",eps_elec%cr3d(i,j,k),&
             " encountered at grid point (",i,",",j,",",k,")"
            CALL stop_program(routineN,moduleN,__LINE__,message,para_env)
          END IF
          rho_tot_zero%cr3d(i,j,k) = rho_solute%cr3d(i,j,k)/eps_elec%cr3d(i,j,k)
          eps_elec%cr3d(i,j,k) = LOG(eps_elec%cr3d(i,j,k))
        END DO
      END DO
    END DO
    !$omp end parallel do

    ! Reassign pointers due to new content
    ln_eps_elec => eps_elec
    NULLIFY (eps_elec)

    ! Build the derivative of ln_eps_elec
    DO i=1,3
      NULLIFY (dln_eps_elec(i)%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool,&
                             dln_eps_elec(i)%pw,&
                             use_data=REALDATA3D,&
                             in_space=REALSPACE,&
                             error=error)
      CALL pw_zero(dln_eps_elec(i)%pw,error=error)
    END DO
    CALL derive(ln_eps_elec,dln_eps_elec,sccs_control%derivative_method,pw_env,input,para_env,error)

    ! Print header for the SCCS cycle
    IF (should_output.AND.(output_unit > 0)) THEN
      IF (print_level > low_print_level) THEN
        WRITE (UNIT=output_unit,FMT="(/,(T3,A,T61,F20.10))")&
         "SCCS| Total electronic charge density ",-tot_rho_elec,&
         "SCCS| Total charge density (solute)   ",-tot_rho_solute
        WRITE (UNIT=output_unit,FMT="(T3,A,T56,F25.3)")&
         "SCCS| Surface of the solute cavity [bohr^2]",cavity_surface,&
         "SCCS|                              [angstrom^2]",&
         cp_unit_from_cp2k(cavity_surface,"angstrom^2",error=error),&
         "SCCS| Volume of the solute cavity  [bohr^3]",cavity_volume,&
         "SCCS|                              [angstrom^3]",&
         cp_unit_from_cp2k(cavity_volume,"angstrom^3",error=error),&
         "SCCS| Volume of the cell           [bohr^3]",cell_volume,&
         "SCCS|                              [angstrom^3]",&
         cp_unit_from_cp2k(cell_volume,"angstrom^3",error=error)
        WRITE (UNIT=output_unit,FMT="(T3,A)")&
         "SCCS|",&
         "SCCS|   Step    Average residual    Maximum residual"
      END IF
    END IF

    ! Get storage for the derivative of the total potential (field) in r-space
    DO i=1,3
      NULLIFY (dphi_tot(i)%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool,&
                             dphi_tot(i)%pw,&
                             use_data=REALDATA3D,&
                             in_space=REALSPACE,&
                             error=error)
    END DO

    ! Reassign work storage to rho_tot, because ln_eps_elec is no longer needed
    rho_tot => ln_eps_elec
    NULLIFY (ln_eps_elec)

    ! Initialise the total electronic density in r-space rho_tot with rho_tot_zero + rho_iter_zero
    CALL pw_copy(rho_tot_zero,rho_tot,error=error)
    CALL pw_axpy(rho_iter_old,rho_tot,error=error)

    ! Main SCCS iteration loop
    iter = 0

    iter_loop: DO

      ! Increment iteration counter
      iter = iter + 1

      ! Check if the requested maximum number of SCCS iterations is reached
      IF (iter > sccs_control%max_iter) THEN
        IF (output_unit > 0) THEN
          WRITE (UNIT=output_unit,FMT="(T3,A,/,T3,A,I0,A)")&
           "SCCS| Maximum number of SCCS iterations reached",&
           "SCCS| Iteration cycle did not converge in ",sccs_control%max_iter," steps"
        END IF
        EXIT iter_loop
      END IF

      ! Calculate derivative of the current total potential in r-space
      CALL pw_poisson_solve(poisson_env=poisson_env,&
                            density=rho_tot,&
                            dvhartree=dphi_tot,&
                            error=error)

      ! Update total charge density (solute plus polarisation) in r-space
      ! based on the iterated polarisation charge density
      f = 0.5_dp/twopi
      rho_delta_avg = 0.0_dp
      rho_delta_max = 0.0_dp
      !$omp parallel do default(none) &
      !$omp             private(i,j,k,rho_delta,rho_iter_new) &
      !$omp             shared(dln_eps_elec,dphi_tot,f,lb,rho_iter_old,ub) &
      !$omp             shared(rho_tot,rho_tot_zero,sccs_control) &
      !$omp             reduction(+:rho_delta_avg) &
      !$omp             reduction(MAX:rho_delta_max)
      DO k=lb(3),ub(3)
        DO j=lb(2),ub(2)
          DO i=lb(1),ub(1)
            rho_iter_new = (dln_eps_elec(1)%pw%cr3d(i,j,k)*dphi_tot(1)%pw%cr3d(i,j,k) +&
                            dln_eps_elec(2)%pw%cr3d(i,j,k)*dphi_tot(2)%pw%cr3d(i,j,k) +&
                            dln_eps_elec(3)%pw%cr3d(i,j,k)*dphi_tot(3)%pw%cr3d(i,j,k))*f
            rho_iter_new = rho_iter_old%cr3d(i,j,k) +&
                           sccs_control%mixing*(rho_iter_new - rho_iter_old%cr3d(i,j,k))
            rho_delta = ABS(rho_iter_new - rho_iter_old%cr3d(i,j,k))
            rho_delta_max = MAX(rho_delta,rho_delta_max)
            rho_delta_avg = rho_delta_avg + rho_delta
            rho_tot%cr3d(i,j,k) = rho_tot_zero%cr3d(i,j,k) + rho_iter_new
            rho_iter_old%cr3d(i,j,k) = rho_iter_new
          END DO
        END DO
      END DO
      !$omp end parallel do

      CALL mp_sum(rho_delta_avg,para_env%group)
      rho_delta_avg = rho_delta_avg/REAL(ngpts,KIND=dp)
      CALL mp_max(rho_delta_max,para_env%group)

      IF (should_output.AND.(output_unit > 0)) THEN
        IF (print_level > low_print_level) THEN
          IF ((ABS(rho_delta_avg) < 1.0E-8_dp).OR.&
              (ABS(rho_delta_avg) >= 1.0E5_dp)) THEN
            WRITE (UNIT=output_unit,FMT="(T3,A,I6,4X,ES16.4,4X,ES16.4)")&
             "SCCS| ",iter,rho_delta_avg,rho_delta_max
          ELSE
            WRITE (UNIT=output_unit,FMT="(T3,A,I6,4X,F16.8,4X,F16.8)")&
             "SCCS| ",iter,rho_delta_avg,rho_delta_max
          END IF
        END IF
      END IF

      ! Check if the SCCS iteration cycle is converged to the requested tolerance
      IF (rho_delta_max <= sccs_control%eps_sccs) THEN
        IF (should_output.AND.(output_unit > 0)) THEN
          WRITE (UNIT=output_unit,FMT="(T3,A,I0,A)")&
           "SCCS| Iteration cycle converged in ",iter," steps"
        END IF
        EXIT iter_loop
      END IF

    END DO iter_loop

    ! Calculate the total Hartree energy, potential, and its derivatives of
    ! the solute and the implicit solvent
    CALL pw_transfer(rho_tot,rho_tot_gspace,error=error)
    IF (calculate_stress_tensor) THEN
      CALL pw_poisson_solve(poisson_env=poisson_env,&
                            density=rho_tot_gspace,&
                            ehartree=energy%sccs_hartree,&
                            vhartree=v_hartree_gspace,&
                            dvhartree=dphi_tot,&
                            h_stress=h_stress,&
                            error=error)
    ELSE
      CALL pw_poisson_solve(poisson_env=poisson_env,&
                            density=rho_tot_gspace,&
                            ehartree=energy%sccs_hartree,&
                            vhartree=v_hartree_gspace,&
                            dvhartree=dphi_tot,&
                            error=error)
    END IF
    phi_pol => work_r3d(5)%pw
    CALL pw_transfer(v_hartree_gspace,phi_pol,error=error)

    ! Calculate the Hartree energy and potential of the solute only
    phi_solute => rho_tot_zero
    CALL pw_zero(phi_solute,error=error)
    NULLIFY (rho_tot_zero)
    CALL pw_poisson_solve(poisson_env=poisson_env,&
                          density=rho_solute,&
                          ehartree=energy%hartree,&
                          vhartree=phi_solute,&
                          error=error)

    ! Calculate the polarisation potential
    ! phi_pol = phi_tot - phi_solute
    CALL pw_axpy(phi_solute,phi_pol,alpha=-1.0_dp,error=error)

    ! Calculate the polarisation charge
    ! rho_pol = rho_tot - rho_solute
    CALL pw_axpy(rho_solute,rho_tot,alpha=-1.0_dp,error=error)
    polarisation_charge = accurate_sum(rho_tot%cr3d)*dvol
    CALL mp_sum(polarisation_charge,para_env%group)

    filename = "POLARISATION_POTENTIAL"
    cube_path = TRIM(print_path)//"%"//TRIM(filename)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,input,TRIM(cube_path),error=error),&
              cp_p_file)) THEN
      append_cube = section_get_lval(input,TRIM(cube_path)//"%APPEND",error=error)
      my_pos_cube="REWIND"
      IF (append_cube) my_pos_cube="APPEND"
      cube_unit = cp_print_key_unit_nr(logger,input,TRIM(cube_path),&
                                       extension=".cube",middle_name=TRIM(filename),&
                                       file_position=my_pos_cube,log_filename=.FALSE.,&
                                       error=error)
      IF (output_unit > 0) THEN
        INQUIRE (UNIT=cube_unit,NAME=filename)
        WRITE (UNIT=output_unit,FMT="(/,T2,A,/,/,T2,A)")&
         "The SCCS polarisation potential is written in cube file format to the file:",TRIM(filename)
      END IF
      CALL cp_pw_to_cube(phi_pol,cube_unit,TRIM(filename),particles=particles,&
                         stride=section_get_ivals(input,TRIM(cube_path)//"%STRIDE",error=error),&
                         error=error)
      CALL cp_print_key_finished_output(cube_unit,logger,input,TRIM(cube_path),error=error)
    END IF

    ! Calculate SCCS polarisation energy
    energy%sccs_pol = 0.5_dp*pw_integral_ab(rho_solute,phi_pol,error=error)

    ! Calculate the Makov-Payne energy correction for charged systems with PBC
    ! Madelung energy of a simple cubic lattice of point charges immersed in neutralising jellium
    energy%sccs_mpc = -0.5_dp*2.8373_dp*tot_rho_solute**2/MINVAL(abc(1:3))

    ! Calculate additional solvation terms
    energy%sccs_cav = sccs_control%gamma_solvent*cavity_surface
    energy%sccs_dis = sccs_control%beta_solvent*cavity_volume
    energy%sccs_rep = sccs_control%alpha_solvent*cavity_surface

    IF (should_output.AND.(output_unit > 0)) THEN
      WRITE (UNIT=output_unit,FMT="(T3,A)")&
       "SCCS|"
      WRITE (UNIT=output_unit,FMT="(T3,A,T56,F25.12)")&
       "SCCS| Polarisation charge",polarisation_charge
      WRITE (UNIT=output_unit,FMT="(T3,A,T56,F25.12)")&
       "SCCS| Hartree energy of the solute only [Hartree]",energy%hartree,&
       "SCCS| Hartree energy of solute and solvent [Hartree]",energy%sccs_hartree
      WRITE (UNIT=output_unit,FMT="(T3,A,T56,F25.12,/,T3,A,T61,F20.3)")&
       "SCCS| Polarisation energy    [Hartree]",energy%sccs_pol,&
       "SCCS|                        [kcal/mol]",&
       cp_unit_from_cp2k(energy%sccs_pol,"kcalmol",error=error),&
       "SCCS| Makov-Payne correction [Hartree]",energy%sccs_mpc,&
       "SCCS|                        [kcal/mol]",&
       cp_unit_from_cp2k(energy%sccs_mpc,"kcalmol",error=error),&
       "SCCS| Cavitation energy      [Hartree]",energy%sccs_cav,&
       "SCCS|                        [kcal/mol]",&
       cp_unit_from_cp2k(energy%sccs_cav,"kcalmol",error=error),&
       "SCCS| Dispersion free energy [Hartree]",energy%sccs_dis,&
       "SCCS|                        [kcal/mol]",&
       cp_unit_from_cp2k(energy%sccs_dis,"kcalmol",error=error),&
       "SCCS| Repulsion free energy  [Hartree]",energy%sccs_rep,&
       "SCCS|                        [kcal/mol]",&
       cp_unit_from_cp2k(energy%sccs_rep,"kcalmol",error=error)
    END IF

    ! Calculate SCCS contribution to the Kohn-Sham potential
    f = -0.25_dp*dvol/twopi
    !$omp parallel do default(none) &
    !$omp             private(dphi2,i,j,k) &
    !$omp             shared(calculate_stress_tensor,f,deps_elec,dphi_tot) &
    !$omp             shared(dvol,lb,ub,v_sccs)
    DO k=lb(3),ub(3)
      DO j=lb(2),ub(2)
        DO i=lb(1),ub(1)
          dphi2 = dphi_tot(1)%pw%cr3d(i,j,k)*dphi_tot(1)%pw%cr3d(i,j,k) +&
                  dphi_tot(2)%pw%cr3d(i,j,k)*dphi_tot(2)%pw%cr3d(i,j,k) +&
                  dphi_tot(3)%pw%cr3d(i,j,k)*dphi_tot(3)%pw%cr3d(i,j,k)
          v_sccs%cr3d(i,j,k) =  v_sccs%cr3d(i,j,k) + f*deps_elec%cr3d(i,j,k)*dphi2
        END DO
      END DO
    END DO
    !$omp end parallel do

    ! Release work storage
    NULLIFY (phi_solute)
    NULLIFY (rho_tot)
    NULLIFY (deps_elec)
    NULLIFY (rho_solute)
    DO i=1,SIZE(work_r3d)
      CALL pw_pool_give_back_pw(auxbas_pw_pool,work_r3d(i)%pw,error=error)
    END DO
    DO i=1,3
      CALL pw_pool_give_back_pw(auxbas_pw_pool,dln_eps_elec(i)%pw,error=error)
      CALL pw_pool_give_back_pw(auxbas_pw_pool,dphi_tot(i)%pw,error=error)
    END DO

    ! Release the SCCS printout environment
    CALL cp_print_key_finished_output(output_unit,logger,input,TRIM(print_path),&
                                      ignore_should_output=should_output,&
                                      error=error)

    CALL timestop(handle)

  END SUBROUTINE sccs

! *****************************************************************************
!> \brief      Calculate the smoothed dielectric function of Andreussi et al.
!> \param rho_elec ...
!> \param eps_elec ...
!> \param deps_elec ...
!> \param epsilon_solvent ...
!> \param rho_max ...
!> \param rho_min ...
!> \par History:
!>      - Creation (16.10.2013,MK)
!>      - Finite difference of isosurfaces implemented (21.12.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.1
! *****************************************************************************
  SUBROUTINE andreussi(rho_elec,eps_elec,deps_elec,epsilon_solvent,rho_max,&
                       rho_min)

    TYPE(pw_type), POINTER                   :: rho_elec, eps_elec, deps_elec
    REAL(KIND=dp), INTENT(IN)                :: epsilon_solvent, rho_max, &
                                                rho_min

    CHARACTER(LEN=*), PARAMETER :: routineN = 'andreussi', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: tol = 1.0E-12_dp

    INTEGER                                  :: handle, i, j, k
    INTEGER, DIMENSION(3)                    :: lb, ub
    REAL(KIND=dp)                            :: diff, dq, dt, f, ln_rho_max, &
                                                ln_rho_min, q, rho, t, x, y

    CALL timeset(routineN,handle)

    f = LOG(epsilon_solvent)/twopi
    diff = rho_max - rho_min
    IF (diff > tol) THEN
      ln_rho_max = LOG(rho_max)
      ln_rho_min = LOG(rho_min)
      q = twopi/(ln_rho_max - ln_rho_min)
      dq = -f*q
    END IF

    lb(1:3) = rho_elec%pw_grid%bounds_local(1,1:3)
    ub(1:3) = rho_elec%pw_grid%bounds_local(2,1:3)

    ! Calculate the dielectric function and its derivative
    !$omp parallel do default(none) &
    !$omp             private(dt,i,j,k,rho,t,x,y) &
    !$omp             shared(deps_elec,diff,dq,eps_elec,epsilon_solvent,f,lb,ub) &
    !$omp             shared(ln_rho_max,rho_elec,q,rho_max,rho_min)
    DO k=lb(3),ub(3)
      DO j=lb(2),ub(2)
        DO i=lb(1),ub(1)
          rho = rho_elec%cr3d(i,j,k)
          IF (rho < rho_min) THEN
            eps_elec%cr3d(i,j,k) = epsilon_solvent
            deps_elec%cr3d(i,j,k) = 0.0_dp
          ELSE IF (rho <= rho_max) THEN
            IF (diff > tol) THEN
              x = LOG(rho)
              y = q*(ln_rho_max - x)
              t = f*(y - SIN(y))
              eps_elec%cr3d(i,j,k) = EXP(t)
              dt = dq*(1.0_dp - COS(y))
              deps_elec%cr3d(i,j,k) = eps_elec%cr3d(i,j,k)*dt/rho
            ELSE
              eps_elec%cr3d(i,j,k) = 1.0_dp
              deps_elec%cr3d(i,j,k) = 0.0_dp
            END IF
          ELSE
            eps_elec%cr3d(i,j,k) = 1.0_dp
            deps_elec%cr3d(i,j,k) = 0.0_dp
          END IF
        END DO
      END DO
    END DO
    !$omp end parallel do

    CALL timestop(handle)

  END SUBROUTINE andreussi

! *****************************************************************************
!> \brief      Calculate the smoothed dielectric function of Fattebert and Gygi
!> \param rho_elec ...
!> \param eps_elec ...
!> \param deps_elec ...
!> \param epsilon_solvent ...
!> \param beta ...
!> \param rho_zero ...
!> \par History:
!>      - Creation (15.10.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! *****************************************************************************
  SUBROUTINE fattebert_gygi(rho_elec,eps_elec,deps_elec,epsilon_solvent,beta,&
                            rho_zero)

    TYPE(pw_type), POINTER                   :: rho_elec, eps_elec, deps_elec
    REAL(KIND=dp), INTENT(IN)                :: epsilon_solvent, beta, &
                                                rho_zero

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fattebert_gygi', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: tol = 1.0E-12_dp

    INTEGER                                  :: handle, i, j, k
    INTEGER, DIMENSION(3)                    :: lb, ub
    REAL(KIND=dp)                            :: df, f, p, q, rho, s, t, &
                                                twobeta

    CALL timeset(routineN,handle)

    df = (1.0_dp - epsilon_solvent)/rho_zero
    f = 0.5_dp*(epsilon_solvent - 1.0_dp)
    q = 1.0_dp/rho_zero
    twobeta = 2.0_dp*beta

    lb(1:3) = rho_elec%pw_grid%bounds_local(1,1:3)
    ub(1:3) = rho_elec%pw_grid%bounds_local(2,1:3)

    ! Calculate the smoothed dielectric function and its derivative
    !$omp parallel do default(none) &
    !$omp             private(i,j,k,p,rho,s,t) &
    !$omp             shared(df,deps_elec,eps_elec,epsilon_solvent,f,lb,ub) &
    !$omp             shared(q,rho_elec,twobeta)
    DO k=lb(3),ub(3)
      DO j=lb(2),ub(2)
        DO i=lb(1),ub(1)
          rho = rho_elec%cr3d(i,j,k)
          IF (rho < tol) THEN
            eps_elec%cr3d(i,j,k) = epsilon_solvent
            deps_elec%cr3d(i,j,k) = 0.0_dp
          ELSE
            s = rho*q
            p = s**twobeta
            t = 1.0_dp/(1.0_dp + p)
            eps_elec%cr3d(i,j,k) = 1.0_dp + f*(1.0_dp + (1.0_dp - p)*t)
            deps_elec%cr3d(i,j,k) = df*twobeta*t*t*p/s
          END IF
        END DO
      END DO
    END DO
    !$omp end parallel do

    CALL timestop(handle)

  END SUBROUTINE fattebert_gygi

! *****************************************************************************
!> \brief      Build the numerical derivative of a function on realspace grid
!> \param f ...
!> \param df ...
!> \param method ...
!> \param pw_env ...
!> \param input ...
!> \param para_env ...
!> \param error ...
!> \par History:
!>      - Creation (15.11.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! *****************************************************************************
  SUBROUTINE derive(f,df,method,pw_env,input,para_env,error)

    TYPE(pw_type), POINTER                   :: f
    TYPE(pw_p_type), DIMENSION(3), &
      INTENT(OUT)                            :: df
    INTEGER, INTENT(IN)                      :: method
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'derive', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: border_points, handle, i, j, k
    INTEGER, DIMENSION(3)                    :: lb, n, ub
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: h
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: drdx, drdy, drdz, r
    TYPE(pw_p_type), DIMENSION(2)            :: work_g1d
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    TYPE(realspace_grid_input_type)          :: input_settings
    TYPE(realspace_grid_p_type), &
      DIMENSION(3)                           :: drs_grid
    TYPE(realspace_grid_type), POINTER       :: rs_grid
    TYPE(section_vals_type), POINTER         :: rs_grid_section

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(f),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)

    ! Perform method specific setup
    SELECT CASE (method)
    CASE (sccs_derivative_cd3,sccs_derivative_cd5,sccs_derivative_cd7)
      NULLIFY (rs_desc)
      NULLIFY (rs_grid)
      rs_grid_section => section_vals_get_subs_vals(input,"DFT%MGRID%RS_GRID",error=error)
      SELECT CASE (method)
      CASE (sccs_derivative_cd3)
        border_points = 1
      CASE (sccs_derivative_cd5)
        border_points = 2
      CASE (sccs_derivative_cd7)
        border_points = 3
      END SELECT
      CALL init_input_type(input_settings,2*border_points+1,rs_grid_section,&
                           1,(/-1,-1,-1/),error)
      CALL rs_grid_create_descriptor(rs_desc,f%pw_grid,input_settings,&
                                     border_points=border_points,error=error)
      CALL rs_grid_create(rs_grid,rs_desc,error=error)
!MK      CALL rs_grid_print(rs_grid,6,error=error)
      CALL rs_pw_transfer(rs_grid,f,pw2rs,error=error)
      DO i=1,3
        CALL rs_grid_create(drs_grid(i)%rs_grid,rs_desc,error=error)
        CALL rs_grid_zero(drs_grid(i)%rs_grid)
      END DO
      lb(1:3) = rs_grid%lb_real(1:3)
      ub(1:3) = rs_grid%ub_real(1:3)
      r => rs_grid%r
      drdx => drs_grid(1)%rs_grid%r
      drdy => drs_grid(2)%rs_grid%r
      drdz => drs_grid(3)%rs_grid%r
    CASE (sccs_derivative_fft)
      lb(1:3) = f%pw_grid%bounds_local(1,1:3)
      ub(1:3) = f%pw_grid%bounds_local(2,1:3)
      NULLIFY (auxbas_pw_pool)
      CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
      ! Get work storage for the 1d grids in g-space (derivative calculation)
      DO i=1,SIZE(work_g1d)
        NULLIFY (work_g1d(i)%pw)
        CALL pw_pool_create_pw(auxbas_pw_pool,&
                               work_g1d(i)%pw,&
                               use_data=COMPLEXDATA1D,&
                               in_space=RECIPROCALSPACE,&
                               error=error)
      END DO
    END SELECT

    ! Calculate the derivatives
    SELECT CASE (method)
    CASE (sccs_derivative_cd3)
      ! 3-point stencil central differences
      h(1:3) = 2.0_dp*f%pw_grid%dr(1:3)
      !$omp parallel do default(none) &
      !$omp             private(i,j,k) &
      !$omp             shared(drdx,drdy,drdz,h,lb,r,ub)
      DO k=lb(3),ub(3)
        DO j=lb(2),ub(2)
          DO i=lb(1),ub(1)
            drdx(i,j,k) = (r(i+1,j,k) - r(i-1,j,k))/h(1)
            drdy(i,j,k) = (r(i,j+1,k) - r(i,j-1,k))/h(2)
            drdz(i,j,k) = (r(i,j,k+1) - r(i,j,k-1))/h(3)
          END DO
        END DO
      END DO
      !$omp end parallel do
    CASE (sccs_derivative_cd5)
      ! 5-point stencil central differences
      h(1:3) = 12.0_dp*f%pw_grid%dr(1:3)
      !$omp parallel do default(none) &
      !$omp             private(i,j,k) &
      !$omp             shared(drdx,drdy,drdz,h,lb,r,ub)
      DO k=lb(3),ub(3)
        DO j=lb(2),ub(2)
          DO i=lb(1),ub(1)
            drdx(i,j,k) = (r(i-2,j,k) - r(i+2,j,k) + 8.0_dp*(r(i+1,j,k) - r(i-1,j,k)))/h(1)
            drdy(i,j,k) = (r(i,j-2,k) - r(i,j+2,k) + 8.0_dp*(r(i,j+1,k) - r(i,j-1,k)))/h(2)
            drdz(i,j,k) = (r(i,j,k-2) - r(i,j,k+2) + 8.0_dp*(r(i,j,k+1) - r(i,j,k-1)))/h(3)
          END DO
        END DO
      END DO
      !$omp end parallel do
    CASE (sccs_derivative_cd7)
      ! 7-point stencil central differences
      h(1:3) = 60.0_dp*f%pw_grid%dr(1:3)
      !$omp parallel do default(none) &
      !$omp             private(i,j,k) &
      !$omp             shared(drdx,drdy,drdz,h,lb,r,ub)
      DO k=lb(3),ub(3)
        DO j=lb(2),ub(2)
          DO i=lb(1),ub(1)
            drdx(i,j,k) = (r(i+3,j,k) - r(i-3,j,k) + 9.0_dp*(r(i-2,j,k) - r(i+2,j,k)) +&
                           45.0_dp*(r(i+1,j,k) - r(i-1,j,k)))/h(1)
            drdy(i,j,k) = (r(i,j+3,k) - r(i,j-3,k) + 9.0_dp*(r(i,j-2,k) - r(i,j+2,k)) +&
                           45.0_dp*(r(i,j+1,k) - r(i,j-1,k)))/h(2)
            drdz(i,j,k) = (r(i,j,k+3) - r(i,j,k-3) + 9.0_dp*(r(i,j,k-2) - r(i,j,k+2)) +&
                           45.0_dp*(r(i,j,k+1) - r(i,j,k-1)))/h(3)
          END DO
        END DO
      END DO
      !$omp end parallel do
    CASE (sccs_derivative_fft)
      ! FFT
      CALL pw_transfer(f,work_g1d(1)%pw,error=error)
      DO i=1,3
        n(:) = 0
        n(i) = 1
        CALL pw_copy(work_g1d(1)%pw,work_g1d(2)%pw,error=error)
        CALL pw_derive(work_g1d(2)%pw,n(:),error=error)
        CALL pw_transfer(work_g1d(2)%pw,df(i)%pw,error=error)
      END DO
    CASE DEFAULT
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Invalid derivative method for SCCS specified",&
                        para_env)
    END SELECT

    ! Perform method specific cleanup
    SELECT CASE (method)
    CASE (sccs_derivative_cd3,sccs_derivative_cd5,sccs_derivative_cd7)
      CALL rs_grid_release(rs_grid,error=error)
      DO i=1,3
        CALL rs_pw_transfer(drs_grid(i)%rs_grid,df(i)%pw,rs2pw,error=error)
        CALL rs_grid_release(drs_grid(i)%rs_grid,error=error)
      END DO
      CALL rs_grid_release_descriptor(rs_desc,error=error)
    CASE (sccs_derivative_fft)
      DO i=1,SIZE(work_g1d)
        CALL pw_pool_give_back_pw(auxbas_pw_pool,work_g1d(i)%pw,error=error)
      END DO
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE derive

! *****************************************************************************
!> \brief      Calculate the finite difference between two isosurfaces of the
!>             electronic density. The smoothed dielectric function of
!>             Andreussi et al. is used as switching function eventually
!>             defining the quantum volume and surface of the cavity.
!> \param rho_elec ...
!> \param norm_drho_elec ...
!> \param drho_elec ...
!> \param dtheta ...
!> \param epsilon_solvent ...
!> \param rho_max ...
!> \param rho_min ...
!> \param delta_rho ...
!> \par History:
!>      - Creation (21.12.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! *****************************************************************************
  SUBROUTINE surface_andreussi(rho_elec,norm_drho_elec,drho_elec,dtheta,&
                               epsilon_solvent,rho_max,rho_min,delta_rho)

    TYPE(pw_type), POINTER                   :: rho_elec, norm_drho_elec
    TYPE(pw_p_type), DIMENSION(3), &
      INTENT(IN)                             :: drho_elec
    TYPE(pw_type), POINTER                   :: dtheta
    REAL(KIND=dp), INTENT(IN)                :: epsilon_solvent, rho_max, &
                                                rho_min, delta_rho

    CHARACTER(LEN=*), PARAMETER :: routineN = 'surface_andreussi', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: tol = 1.0E-12_dp

    INTEGER                                  :: handle, i, j, k, l
    INTEGER, DIMENSION(3)                    :: lb, ub
    REAL(KIND=dp)                            :: diff, e, eps_elec, f, &
                                                ln_rho_max, ln_rho_min, q, &
                                                rho, t, x, y
    REAL(KIND=dp), DIMENSION(2)              :: theta

    CALL timeset(routineN,handle)

    e = epsilon_solvent - 1.0_dp
    f = LOG(epsilon_solvent)/twopi
    diff = rho_max - rho_min
    IF (diff > tol) THEN
      ln_rho_max = LOG(rho_max)
      ln_rho_min = LOG(rho_min)
      q = twopi/(ln_rho_max - ln_rho_min)
    END IF

    lb(1:3) = rho_elec%pw_grid%bounds_local(1,1:3)
    ub(1:3) = rho_elec%pw_grid%bounds_local(2,1:3)

    ! Calculate finite difference between two isosurfaces
    !$omp parallel do default(none) &
    !$omp             private(eps_elec,i,j,k,l,rho,t,theta,x,y) &
    !$omp             shared(delta_rho,diff,drho_elec,dtheta,e,epsilon_solvent,f,lb) &
    !$omp             shared(ln_rho_max,norm_drho_elec,rho_elec,q,rho_max,rho_min,ub)
    DO k=lb(3),ub(3)
      DO j=lb(2),ub(2)
        DO i=lb(1),ub(1)
          DO l=1,2
            rho = rho_elec%cr3d(i,j,k) + (REAL(l,KIND=dp) - 1.5_dp)*delta_rho
            IF (rho < rho_min) THEN
              eps_elec = epsilon_solvent
            ELSE IF (rho <= rho_max) THEN
              IF (diff > tol) THEN
                x = LOG(rho)
                y = q*(ln_rho_max - x)
                t = f*(y - SIN(y))
                eps_elec = EXP(t)
              ELSE
                eps_elec = 1.0_dp
              END IF
            ELSE
              eps_elec = 1.0_dp
            END IF
            theta(l) = (epsilon_solvent - eps_elec)/e
          END DO
          dtheta%cr3d(i,j,k) = (theta(2) - theta(1))*norm_drho_elec%cr3d(i,j,k)/delta_rho
        END DO
      END DO
    END DO
    !$omp end parallel do

    CALL timestop(handle)

  END SUBROUTINE surface_andreussi

! *****************************************************************************
!> \brief      Calculate the finite difference between two isosurfaces of the
!>             the electronic density. The smoothed dielectric function of
!>             Fattebert and Gygi is used as switching function eventually
!>             defining the quantum volume and surface of the cavity.
!> \param rho_elec ...
!> \param norm_drho_elec ...
!> \param drho_elec ...
!> \param dtheta ...
!> \param epsilon_solvent ...
!> \param beta ...
!> \param rho_zero ...
!> \param delta_rho ...
!> \par History:
!>      - Creation (21.12.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! *****************************************************************************
  SUBROUTINE surface_fattebert_gygi(rho_elec,norm_drho_elec,drho_elec,dtheta,&
                                    epsilon_solvent,beta,rho_zero,delta_rho)

    TYPE(pw_type), POINTER                   :: rho_elec, norm_drho_elec
    TYPE(pw_p_type), DIMENSION(3), &
      INTENT(IN)                             :: drho_elec
    TYPE(pw_type), POINTER                   :: dtheta
    REAL(KIND=dp), INTENT(IN)                :: epsilon_solvent, beta, &
                                                rho_zero, delta_rho

    CHARACTER(LEN=*), PARAMETER :: routineN = 'surface_fattebert_gygi', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: tol = 1.0E-12_dp

    INTEGER                                  :: handle, i, j, k, l
    INTEGER, DIMENSION(3)                    :: lb, ub
    REAL(KIND=dp)                            :: e, eps_elec, f, p, q, rho, s, &
                                                t, twobeta
    REAL(KIND=dp), DIMENSION(2)              :: theta

    CALL timeset(routineN,handle)

    e = epsilon_solvent - 1.0_dp
    f = 0.5_dp*e
    q = 1.0_dp/rho_zero
    twobeta = 2.0_dp*beta

    lb(1:3) = rho_elec%pw_grid%bounds_local(1,1:3)
    ub(1:3) = rho_elec%pw_grid%bounds_local(2,1:3)

    ! Calculate finite difference between two isosurfaces
    !$omp parallel do default(none) &
    !$omp             private(eps_elec,i,j,k,l,p,rho,s,t,theta) &
    !$omp             shared(delta_rho,drho_elec,dtheta,e,epsilon_solvent,f,lb) &
    !$omp             shared(norm_drho_elec,q,rho_elec,twobeta,ub)
    DO k=lb(3),ub(3)
      DO j=lb(2),ub(2)
        DO i=lb(1),ub(1)
          DO l=1,2
            rho = rho_elec%cr3d(i,j,k) + (REAL(l,KIND=dp) - 1.5_dp)*delta_rho
            IF (rho < tol) THEN
              eps_elec = epsilon_solvent
            ELSE
              s = rho*q
              p = s**twobeta
              t = 1.0_dp/(1.0_dp + p)
              eps_elec = 1.0_dp + f*(1.0_dp + (1.0_dp - p)*t)
            END IF
            theta(l) = (epsilon_solvent - eps_elec)/e
          END DO
          dtheta%cr3d(i,j,k) = (theta(2) - theta(1))*norm_drho_elec%cr3d(i,j,k)/delta_rho
        END DO
      END DO
    END DO
    !$omp end parallel do

    CALL timestop(handle)

  END SUBROUTINE surface_fattebert_gygi

END MODULE qs_sccs
