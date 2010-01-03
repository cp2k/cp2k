!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
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
!> \author Teodoro Laino
! *****************************************************************************
MODULE cp2k_debug
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_calc_energy_force,&
                                             force_env_calc_num_pressure
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_constants,                 ONLY: do_stress_analytical,&
                                             do_stress_diagonal_anal,&
                                             do_stress_diagonal_numer,&
                                             do_stress_numerical
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type,&
                                             write_qs_particle_coordinates
  USE virial_types,                    ONLY: cp_virial,&
                                             virial_create,&
                                             virial_release,&
                                             virial_type,&
                                             zero_virial
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp2k_debug'
  REAL(KIND=dp), PRIVATE, PARAMETER :: MaxErr = 5.0_dp

  PUBLIC :: cp2k_debug_energy_and_forces

CONTAINS

! *****************************************************************************
  SUBROUTINE cp2k_debug_energy_and_forces(force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp2k_debug_energy_and_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ip, iseq, iw, j, k, stat, &
                                                stress_tensor
    LOGICAL                                  :: debug_forces, &
                                                debug_stress_tensor, failure, &
                                                harris_flag, skip
    REAL(kind=dp)                            :: dx, std_value
    REAL(kind=dp), DIMENSION(2)              :: numer_energy
    REAL(kind=dp), DIMENSION(3)              :: Err, my_maxerr
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: analyt_forces, numer_forces
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(section_vals_type), POINTER         :: harris_section, root_section, &
                                                subsys_section
    TYPE(virial_type), POINTER               :: virial, virial_analytical, &
                                                virial_numerical

    failure = .FALSE.
    NULLIFY(analyt_forces, numer_forces, subsys, particles)
    root_section=>force_env%root_section
    IF (.NOT.failure) THEN
       harris_section => section_vals_get_subs_vals(force_env%force_env_section, &
            "DFT%QS%HARRIS",error=error)
       subsys_section => section_vals_get_subs_vals(force_env%force_env_section, &
            "SUBSYS",error=error)
       CALL section_vals_val_get(harris_section, "ACTIVATE", &
            l_val=harris_flag, error=error)

       CALL section_vals_val_get(root_section,"DEBUG%DEBUG_STRESS_TENSOR",&
            l_val=debug_stress_tensor, error=error)
       CALL section_vals_val_get(root_section,"DEBUG%DEBUG_FORCES",&
            l_val=debug_forces, error=error)
       CALL section_vals_val_get(root_section,"DEBUG%DX",&
            r_val=dx, error=error)

       logger => cp_error_get_logger(error)
       iw=cp_print_key_unit_nr(logger,root_section,"DEBUG%PROGRAM_RUN_INFO",&
            extension=".log",error=error)
       IF (debug_stress_tensor) THEN
          ! To debug stress tensor the stress tensor calculation must be
          ! first enabled..
          CALL section_vals_val_get(force_env%force_env_section,"STRESS_TENSOR",&
               i_val=stress_tensor,error=error)
          skip = .FALSE.
          SELECT CASE (stress_tensor)
          CASE (do_stress_analytical, do_stress_diagonal_anal)
             ! OK..
          CASE (do_stress_numerical, do_stress_diagonal_numer)
             ! Nothing to check..
             CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                  "Numerical Stress Tensor was requested in the FORCE_EVAL section. "//&
                  "Nothing to debug!"//&
CPSourceFileRef,&
             only_ionode=.TRUE.)
             skip = .TRUE.
          CASE DEFAULT
             CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                  "Stress Tensor calculation was not enabled in the FORCE_EVAL section. "//&
                  "Nothing to debug!"//&
CPSourceFileRef,&
                  only_ionode=.TRUE.)
             skip = .TRUE.
          END SELECT

          IF (.NOT.skip) THEN
             ! Compute the analytical stress tensor
             CALL force_env_calc_energy_force(force_env,calc_force=.TRUE.,&
                  error=error)
             
             ! Retrieve the analytical virial
             CALL force_env_get(force_env, virial=virial, error=error)
             CALL virial_create(virial_analytical, error=error)
             CALL zero_virial(virial_analytical)
             CALL cp_virial(virial, virial_analytical)
             
             ! Debug Stress Tensor (numerical Vs Analytical)
             CALL force_env_calc_num_pressure(force_env,dx=dx,error=error)
             
             ! Retrieve the numerical virial
             CALL force_env_get(force_env, virial=virial, error=error)
             CALL virial_create(virial_numerical, error=error)
             CALL zero_virial(virial_numerical)
             CALL cp_virial(virial, virial_numerical)
             
             ! Print results
             IF (iw>0) THEN
                WRITE (UNIT=iw,FMT="(/,'  *** DEBUG PV_VIRIAL - START ***',/)")
                WRITE (UNIT=iw,FMT="((T34,A))") "Numerical pv_virial"
                WRITE (UNIT=iw,FMT="((T34,3F16.10))") (virial_numerical%pv_virial(i,1:3), i=1,3)
                WRITE (UNIT=iw,FMT="(/)")
                WRITE (UNIT=iw,FMT="((T34,A))") "Analytical pv_virial"
                WRITE (UNIT=iw,FMT="((T34,3F16.10))") (virial_analytical%pv_virial(i,1:3), i=1,3)
                WRITE (UNIT=iw,FMT="(/)")
                WRITE (UNIT=iw,FMT="((T34,A))") "Difference pv_virial"
                WRITE (UNIT=iw,FMT="((T34,3F16.10))")&
                     (virial_numerical%pv_virial(i,1:3)-virial_analytical%pv_virial(i,1:3), i=1,3)
                WRITE (UNIT=iw,FMT="(/,'  *** DEBUG PV_VIRIAL - END   ***',/)")
             END IF
             
             ! Checkd and abort (in case..)
             DO i = 1,3
                Err = 0.0_dp
                DO K = 1, 3
                   IF (ABS(virial_numerical%pv_virial(i,K)) >= 1.0E-6_dp) THEN
                      Err(K) = (virial_numerical%pv_virial(i,K)-&
                           virial_analytical%pv_virial(i,K))/virial_numerical%pv_virial(i,K)*100.0_dp
                   END IF
                END DO
                CPPostcondition(ABS(Err(1))<=MaxErr,cp_failure_level,routineP,error,failure)
                CPPostcondition(ABS(Err(2))<=MaxErr,cp_failure_level,routineP,error,failure)
                CPPostcondition(ABS(Err(3))<=MaxErr,cp_failure_level,routineP,error,failure)
             END DO
             CALL virial_release(virial_analytical, error=error)
             CALL virial_release(virial_numerical, error=error)
          END IF
       END IF
       
       IF (debug_forces) THEN
          ! Debug Forces (numerical Vs Analytical)
          CALL force_env_get(force_env,subsys=subsys,cell=cell,error=error)
          particles => subsys%particles%els
          CALL write_qs_particle_coordinates(particles,cell,subsys_section,"DEBUG",error)
          !
          ! First evaluate energy and forces...
          !
          CALL force_env_calc_energy_force(force_env,calc_force=.TRUE.,&
               error=error)
          !
          ! Copy forces in array and start the numerical calculation
          !
          IF (ASSOCIATED(analyt_forces)) DEALLOCATE(analyt_forces)
          iseq = subsys%particles%n_els
          ALLOCATE( analyt_forces(iseq,3), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          iseq = 0
          DO ip = 1, subsys%particles%n_els
             iseq = iseq + 1
             analyt_forces(iseq,:) = particles(ip)%f
          END DO
          !
          ! Loop on atoms and coordinates
          !
          iseq = 0
          IF (ASSOCIATED( numer_forces)) DEALLOCATE( numer_forces)
          ALLOCATE(  numer_forces(subsys%particles%n_els,3), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          Atom: DO ip = 1, subsys%particles%n_els
             iseq = iseq + 1
             Coord: DO k = 1, 3
                numer_energy = 0.0_dp
                std_value = particles(ip)%r(k)
                DO j = 1, 2
                   particles(ip)%r(k) = std_value - (-1.0_dp)**j * Dx
                   CALL write_qs_particle_coordinates(particles,cell,subsys_section,"DEBUG",error)
                   ! Compute forces..
                   CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,&
                        consistent_energies=(.NOT.harris_flag), error=error)
                   CALL force_env_get (force_env, potential_energy=numer_energy(j), error=error)
                END DO
                particles(ip)%r(k) = std_value
                numer_forces(ip,k) = - (numer_energy(1) - numer_energy(2) ) / (2.0_dp*Dx)
                IF (iw>0) THEN
                   WRITE(iw,'(A,I5,A,I5,4F15.9)')"ATOM NUMBER ::",iseq," COORD:",k,numer_energy,&
                        numer_forces(ip,k), analyt_forces(iseq,k)
                ENDIF
             END DO Coord
             !
             ! Check analytical forces Vs numerical forces
             !
             IF (iw>0) THEN
                WRITE(iw,'(A,I5,6F15.9)')"ATOM NUMBER ::",iseq,analyt_forces(iseq,:), numer_forces(ip,:)
             ENDIF
          END DO Atom
          IF (iw>0) THEN
             WRITE(iw,'(A,I5)')"SUMMARY DEBUG :: "
          ENDIF
          iseq = iseq - subsys%particles%n_els
          DO ip = 1, subsys%particles%n_els
             iseq = iseq + 1
             IF (iw>0) THEN
                WRITE(iw,'(A,I5,9F12.6)')"ATOM NUMBER ::",iseq,analyt_forces(iseq,:), numer_forces(ip,:),&
                     analyt_forces(iseq,:) - numer_forces(ip,:)
             ENDIF
          END DO
          !
          ! Runtime check...
          !
          DO ip =  1, subsys%particles%n_els
             Err = 0.0_dp
             DO K = 1, 3
                IF (ABS(numer_forces(ip,K)) >= 1.0E-6_dp) THEN
                   Err(K) = (analyt_forces(ip,K)-numer_forces(ip,K))/numer_forces(ip,K)*100.0_dp
                END IF
             END DO
             IF (iw>0) THEN
100             FORMAT(I5,F15.9," ( ",F7.2," ) ",F15.9," ( ",F7.2," ) ",F15.9," ( ",F7.2," ) ")
                WRITE(iw,100)ip,Analyt_Forces(Ip,1),Err(1),&
                     Analyt_Forces(Ip,2),Err(2),&
                     Analyt_Forces(Ip,3),Err(3)
             ENDIF
             my_MaxErr = MaxErr
             IF (ABS(Analyt_Forces(Ip,1)) <= 0.0001_dp ) my_MaxErr(1)=my_MaxErr(1)*5.0_dp
             IF (ABS(Analyt_Forces(Ip,2)) <= 0.0001_dp ) my_MaxErr(2)=my_MaxErr(2)*5.0_dp
             IF (ABS(Analyt_Forces(Ip,3)) <= 0.0001_dp ) my_MaxErr(3)=my_MaxErr(3)*5.0_dp
             CPPostcondition(ABS(Err(1))<=my_MaxErr(1),cp_failure_level,routineP,error,failure)
             CPPostcondition(ABS(Err(2))<=my_MaxErr(2),cp_failure_level,routineP,error,failure)
             CPPostcondition(ABS(Err(3))<=my_MaxErr(3),cp_failure_level,routineP,error,failure)
          END DO
          
          IF (ASSOCIATED(analyt_forces)) DEALLOCATE(analyt_forces)
          IF (ASSOCIATED( numer_forces)) DEALLOCATE( numer_forces)
       END IF
       CALL cp_print_key_finished_output(iw,logger,root_section,&
            "DEBUG%PROGRAM_RUN_INFO", error=error)
    END IF
  END SUBROUTINE cp2k_debug_energy_and_forces

END MODULE cp2k_debug
