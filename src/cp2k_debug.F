!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Debug energy and derivatives w.r.t. finite differences
!> \note
!>      Use INTERPOLATION USE_GUESS, in order to perform force and energy
!>      calculations with the same density. This is not compulsory when iterating
!>      to selfconsistency, but essential in the non-selfconsistent case [08.2005,TdK].
!> \par History
!>      12.2004 created [tlaino]
!>      08.2005 consistent_energies option added, to allow FD calculations
!>              with the correct energies in the non-selfconsistent case, but
!>              keep in mind, that the QS energies and forces are then NOT
!>              consistent to each other [TdK].
!>      08.2005 In case the Harris functional is used, consistent_energies is
!>              et to .FALSE., otherwise the QS energies are spuriously used [TdK].
!>      01.2015 Remove Harris functional option
!>      - Revised (20.11.2013,MK)
!> \author Teodoro Laino
! *****************************************************************************
MODULE cp2k_debug
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force,&
                                             force_env_calc_num_pressure
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             use_qs_force
  USE input_constants,                 ONLY: do_stress_analytical,&
                                             do_stress_diagonal_anal,&
                                             do_stress_diagonal_numer,&
                                             do_stress_numerical
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE particle_methods,                ONLY: write_fist_particle_coordinates,&
                                             write_qs_particle_coordinates
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE termination,                     ONLY: print_warning,&
                                             stop_program
  USE virial_types,                    ONLY: cp_virial,&
                                             virial_create,&
                                             virial_release,&
                                             virial_set,&
                                             virial_type,&
                                             zero_virial
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp2k_debug'
  REAL(KIND=dp), PRIVATE, PARAMETER :: maxerr = 5.0_dp

  PUBLIC :: cp2k_debug_energy_and_forces

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param force_env ...
! *****************************************************************************
  SUBROUTINE cp2k_debug_energy_and_forces(force_env)

    TYPE(force_env_type), POINTER            :: force_env

    CHARACTER(len=*), PARAMETER :: routineN = 'cp2k_debug_energy_and_forces', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=3*default_string_length)   :: message
    INTEGER                                  :: i, ip, iw, j, k, np, stat, &
                                                stress_tensor
    LOGICAL                                  :: check_failed, debug_forces, &
                                                debug_stress_tensor, failure, &
                                                skip, stop_on_mismatch
    REAL(KIND=dp)                            :: dx, eps_no_error_check, &
                                                std_value
    REAL(KIND=dp), DIMENSION(2)              :: numer_energy
    REAL(KIND=dp), DIMENSION(3)              :: err, my_maxerr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: analyt_forces, numer_forces
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(section_vals_type), POINTER         :: root_section, subsys_section
    TYPE(virial_type), POINTER               :: virial, virial_analytical, &
                                                virial_numerical

    failure = .FALSE.

    NULLIFY (analyt_forces, numer_forces, subsys, particles)
    
    root_section => force_env%root_section


    CALL force_env_get(force_env,para_env=para_env,subsys=subsys,cell=cell)
    subsys_section => section_vals_get_subs_vals(force_env%force_env_section, &
                                                 "SUBSYS")

    CALL section_vals_val_get(root_section,"DEBUG%DEBUG_STRESS_TENSOR",&
                              l_val=debug_stress_tensor)
    CALL section_vals_val_get(root_section,"DEBUG%DEBUG_FORCES",&
                              l_val=debug_forces)
    CALL section_vals_val_get(root_section,"DEBUG%DX",&
                              r_val=dx)
    dx = ABS(dx)
    CALL section_vals_val_get(root_section,"DEBUG%EPS_NO_ERROR_CHECK",&
                              r_val=eps_no_error_check)
    eps_no_error_check = MAX(eps_no_error_check,EPSILON(0.0_dp))
    CALL section_vals_val_get(root_section,"DEBUG%STOP_ON_MISMATCH",&
                              l_val=stop_on_mismatch)

    logger => cp_get_default_logger()
    iw = cp_print_key_unit_nr(logger,root_section,"DEBUG%PROGRAM_RUN_INFO",&
                              extension=".log")
    IF (debug_stress_tensor) THEN
       ! To debug stress tensor the stress tensor calculation must be
       ! first enabled..
       CALL section_vals_val_get(force_env%force_env_section,"STRESS_TENSOR",&
                                 i_val=stress_tensor)
       skip = .FALSE.
       SELECT CASE (stress_tensor)
       CASE (do_stress_analytical, do_stress_diagonal_anal)
          ! OK
       CASE (do_stress_numerical, do_stress_diagonal_numer)
          ! Nothing to check
          message = "Numerical stress tensor was requested in the FORCE_EVAL section. "//&
                    "Nothing to debug!"
          CALL print_warning(routineN,moduleN,__LINE__,TRIM(message),para_env)
          skip = .TRUE.
       CASE DEFAULT
          message = "Stress tensor calculation was not enabled in the FORCE_EVAL section. "//&
                    "Nothing to debug!"
          CALL print_warning(routineN,moduleN,__LINE__,TRIM(message),para_env)
          skip = .TRUE.
       END SELECT

       IF (.NOT.skip) THEN

          ! Compute the analytical stress tensor
          CALL cp_subsys_get(subsys,virial=virial)
          CALL virial_set(virial,pv_numer=.FALSE.)
          CALL force_env_calc_energy_force(force_env,&
                                           calc_force=.TRUE.,&
                                           calc_stress_tensor=.TRUE.)

          ! Retrieve the analytical virial
          CALL virial_create(virial_analytical)
          CALL zero_virial(virial_analytical)
          CALL cp_virial(virial,virial_analytical)

          ! Debug stress tensor (numerical vs analytical)
          CALL virial_set(virial,pv_numer=.TRUE.)
          CALL force_env_calc_num_pressure(force_env,dx=dx)

          ! Retrieve the numerical virial
          CALL cp_subsys_get(subsys,virial=virial)
          CALL virial_create(virial_numerical)
          CALL zero_virial(virial_numerical)
          CALL cp_virial(virial,virial_numerical)

          ! Print results
          IF (iw > 0) THEN
             WRITE (UNIT=iw,FMT="((T2,A))")&
              "DEBUG| Numerical pv_virial"
             WRITE (UNIT=iw,FMT="((T2,A,3F16.10))")&
              ("DEBUG|",virial_numerical%pv_virial(i,1:3),i=1,3)
             WRITE (UNIT=iw,FMT="(/,(T2,A))")&
              "DEBUG| Analytical pv_virial"
             WRITE (UNIT=iw,FMT="((T2,A,3F16.10))")&
              ("DEBUG|",virial_analytical%pv_virial(i,1:3),i=1,3)
             WRITE (UNIT=iw,FMT="(/,(T2,A))")&
              "DEBUG| Difference pv_virial"
             WRITE (UNIT=iw,FMT="((T2,A,3F16.10))")&
              ("DEBUG|",virial_numerical%pv_virial(i,1:3) - virial_analytical%pv_virial(i,1:3),i=1,3)
             WRITE (UNIT=iw,FMT="(T2,A,T40,F16.10,/)")&
              "DEBUG| Sum of differences:",&
              SUM(ABS(virial_numerical%pv_virial(:,:) - virial_analytical%pv_virial(:,:)))
          END IF

          ! Check and abort on failure
          check_failed = .FALSE.
          DO i=1,3
             err(:) = 0.0_dp
             DO k=1,3
                IF (ABS(virial_analytical%pv_virial(i,k)) >= eps_no_error_check) THEN
                   err(k) = 100.0_dp*(virial_numerical%pv_virial(i,k) - virial_analytical%pv_virial(i,k))/&
                            virial_analytical%pv_virial(i,k)
                END IF
             END DO
             IF (ANY(ABS(err(1:3)) > maxerr)) check_failed = .TRUE.
          END DO
          CALL virial_release(virial_analytical)
          CALL virial_release(virial_numerical)
          IF (check_failed) THEN
             message = "A mismatch between the analytical and the numerical "//&
                       "stress tensor has been detected. Check the implementation "//&
                       "of the stress tensor"
             IF (stop_on_mismatch) THEN
                CALL stop_program(routineN,moduleN,__LINE__,TRIM(message),para_env)
             ELSE
                CALL print_warning(routineN,moduleN,__LINE__,TRIM(message),para_env)
             END IF
          END IF
       END IF
    END IF

    IF (debug_forces) THEN
       ! Debug forces (numerical vs analytical)
       particles => subsys%particles%els
       SELECT CASE (force_env%in_use)
       CASE (use_qs_force)
          CALL get_qs_env(force_env%qs_env,qs_kind_set=qs_kind_set)
          CALL write_qs_particle_coordinates(particles,qs_kind_set,subsys_section,"DEBUG")
       CASE DEFAULT
          CALL write_fist_particle_coordinates(particles,subsys_section,charges=Null())
       END SELECT
       ! First evaluate energy and forces
       CALL force_env_calc_energy_force(force_env,&
                                        calc_force=.TRUE.,&
                                        calc_stress_tensor=.FALSE.)
       ! Copy forces in array and start the numerical calculation
       IF (ASSOCIATED(analyt_forces)) DEALLOCATE(analyt_forces)
       np = subsys%particles%n_els
       ALLOCATE (analyt_forces(np,3),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,failure)
       DO ip=1,np
          analyt_forces(ip,1:3) = particles(ip)%f(1:3)
       END DO
       ! Loop on atoms and coordinates
       IF (ASSOCIATED(numer_forces)) DEALLOCATE (numer_forces)
       ALLOCATE (numer_forces(subsys%particles%n_els,3),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,failure)
       Atom: DO ip=1,np
          Coord: DO k=1,3
             numer_energy = 0.0_dp
             std_value = particles(ip)%r(k)
             DO j=1,2
                particles(ip)%r(k) = std_value - (-1.0_dp)**j*dx
                SELECT CASE (force_env%in_use)
                CASE (use_qs_force)
                   CALL get_qs_env(force_env%qs_env,qs_kind_set=qs_kind_set)
                   CALL write_qs_particle_coordinates(particles,qs_kind_set,subsys_section,"DEBUG")
                CASE DEFAULT
                   CALL write_fist_particle_coordinates(particles,subsys_section,charges=Null())
                END SELECT
                ! Compute energy
                CALL force_env_calc_energy_force(force_env,&
                                                 calc_force=.FALSE.,&
                                                 calc_stress_tensor=.FALSE.,&
                                                 consistent_energies=.TRUE.)
                CALL force_env_get(force_env,potential_energy=numer_energy(j))
             END DO
             particles(ip)%r(k) = std_value
             numer_forces(ip,k) = -0.5_dp*(numer_energy(1) - numer_energy(2))/dx
             IF (iw > 0) THEN
                WRITE (UNIT=iw,FMT="(T2,A,T17,A,F7.4,A,T34,A,F7.4,A,T52,A,T68,A)")&
                 "DEBUG| Atom","E("//ACHAR(119+k)//" +",dx,")",&
                               "E("//ACHAR(119+k)//" -",dx,")",&
                               "f(numerical)","f(analytical)"
                WRITE (UNIT=iw,FMT="(T2,A,I5,4(1X,F16.8))")&
                 "DEBUG|",ip,numer_energy(1:2),numer_forces(ip,k),analyt_forces(ip,k)
             END IF
          END DO Coord
          ! Check analytical forces using the numerical forces as reference
          my_maxerr = maxerr
          err(1:3) = 0.0_dp
          DO k=1,3
             ! Calculate percentage but ignore very small force values
             IF (ABS(analyt_forces(ip,k)) >= eps_no_error_check) THEN
                err(k) = 100.0_dp*(numer_forces(ip,k) - analyt_forces(ip,k))/analyt_forces(ip,k)
             END IF
             ! Increase error tolerance for small force values
             IF (ABS(analyt_forces(ip,k)) <= 0.0001_dp) my_maxerr(k) = 5.0_dp*my_maxerr(k)
          END DO
          IF (iw > 0) THEN
             WRITE (UNIT=iw,FMT="(/,T2,A)")&
              "DEBUG| Atom  Coordinate   f(numerical)   f(analytical)   Difference   Error [%]"
             DO k=1,3
                WRITE (UNIT=iw,FMT="(T2,A,I5,T19,A1,T26,F14.8,T42,F14.8,T57,F12.8,T71,F10.2)")&
                 "DEBUG|",ip,ACHAR(119+k),numer_forces(ip,k),analyt_forces(ip,k),&
                 analyt_forces(ip,k) - numer_forces(ip,k),err(k)
             END DO
          END IF
          IF (ANY(ABS(err(1:3)) > my_maxerr(1:3))) THEN
             message = "A mismatch between analytical and numerical forces "//&
                       "has been detected. Check the implementation of the "//&
                       "analytical force calculation"
             IF (stop_on_mismatch) THEN
                CALL stop_program(routineN,moduleN,__LINE__,TRIM(message),para_env)
             ELSE
                CALL print_warning(routineN,moduleN,__LINE__,TRIM(message),para_env)
             END IF
          END IF
       END DO Atom
       ! Print summary
       IF (iw > 0) THEN
          WRITE (UNIT=iw,FMT="(/,(T2,A))")&
           "DEBUG|======================== BEGIN OF SUMMARY ===============================",&
           "DEBUG| Atom  Coordinate   f(numerical)   f(analytical)   Difference   Error [%]"
          DO ip=1,np
             err(1:3) = 0.0_dp
             DO k=1,3
                IF (ABS(analyt_forces(ip,k)) >= eps_no_error_check) THEN
                   err(k) = 100_dp*(numer_forces(ip,k) - analyt_forces(ip,k))/analyt_forces(ip,k)
                   WRITE (UNIT=iw,FMT="(T2,A,I5,T19,A1,T26,F14.8,T42,F14.8,T57,F12.8,T71,F10.2)")&
                    "DEBUG|",ip,ACHAR(119+k),numer_forces(ip,k),analyt_forces(ip,k),&
                    analyt_forces(ip,k) - numer_forces(ip,k),err(k)
                ELSE
                   WRITE (UNIT=iw,FMT="(T2,A,I5,T19,A1,T26,F14.8,T42,F14.8,T57,F12.8)")&
                    "DEBUG|",ip,ACHAR(119+k),numer_forces(ip,k),analyt_forces(ip,k),&
                    analyt_forces(ip,k) - numer_forces(ip,k)
                END IF
             END DO
          END DO
          WRITE (UNIT=iw,FMT="(T2,A)")&
           "DEBUG|======================== END OF SUMMARY ================================="
       END IF
       ! Release work storage
       IF (ASSOCIATED(analyt_forces)) DEALLOCATE(analyt_forces)
       IF (ASSOCIATED(numer_forces)) DEALLOCATE(numer_forces)
    END IF
    CALL cp_print_key_finished_output(iw,logger,root_section,&
                                      "DEBUG%PROGRAM_RUN_INFO")

  END SUBROUTINE cp2k_debug_energy_and_forces

END MODULE cp2k_debug
