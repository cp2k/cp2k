!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****h* cp2k/qs_debug [1.0] *
!!
!!   NAME
!!     qs_debug
!!
!!   FUNCTION
!!     Debug energy and derivatives w.r.t. finite differences
!!
!!   NOTES
!!     Use INTERPOLATION USE_GUESS, in order to perform force and energy
!!     calculations with the same density. This is not compulsory when iterating
!!     to selfconsistency, but essential in the non-selfconsistent case [08.2005,TdK].
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     12.2004 created [tlaino]
!!     08.2005 consistent_energies option added, to allow FD calculations
!!             with the correct energies in the non-selfconsistent case, but
!!             keep in mind, that the QS energies and forces are then NOT
!!             consistent to each other [TdK].
!!     08.2005 In case the Harris functional is used, consistent_energies is
!!             et to .FALSE., otherwise the QS energies are spuriously used [TdK].
!!
!!
!!   SOURCE
!****************************************************************************
MODULE qs_debug
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qs_numerical_pressure,           ONLY: qs_calc_numerical_pressure
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_debug'

  PUBLIC :: qs_debug_energy_and_forces
  REAL(KIND=dp), PRIVATE, PARAMETER :: MaxErr = 5.0_dp

CONTAINS

  SUBROUTINE qs_debug_energy_and_forces(force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_debug_energy_and_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip, iseq, iw, j, k, stat
    LOGICAL                                  :: failure, harris_flag, &
                                                num_pressure
    REAL(kind=dp)                            :: dx, std_value
    REAL(kind=dp), DIMENSION(2)              :: numer_energy
    REAL(kind=dp), DIMENSION(3)              :: Err, my_maxerr
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: analyt_forces, numer_forces
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(section_vals_type), POINTER         :: harris_section, root_section

    failure = .FALSE.
    NULLIFY(analyt_forces, numer_forces, subsys, particles)
    para_env=>force_env%para_env
    root_section=>force_env%root_section
    IF (.NOT.failure) THEN
       harris_section => section_vals_get_subs_vals(force_env%force_env_section, &
            "DFT%QS%HARRIS",error=error)
       CALL section_vals_val_get(harris_section, "ACTIVATE", &
            l_val=harris_flag, error=error)

       CALL section_vals_val_get(root_section,"DEBUG%NUMERICAL_PRESSURE",&
            l_val=num_pressure, error=error)
       CALL section_vals_val_get(root_section,"DEBUG%DX",&
            r_val=dx, error=error)
       IF (num_pressure) CALL qs_calc_numerical_pressure(force_env,error=error)
       logger => cp_error_get_logger(error)
       iw=cp_print_key_unit_nr(logger,root_section,"GLOBAL%PROGRAM_RUN_INFO",&
            extension=".log",error=error)
       !
       ! First evaluate energy and forces...
       !
       CALL force_env_calc_energy_force(force_env,calc_force=.TRUE.,&
            error=error)
       !
       ! Copy forces in array and start the numerical calculation
       !
       CALL force_env_get(force_env,subsys=subsys,error=error)
       IF (ASSOCIATED(analyt_forces)) DEALLOCATE(analyt_forces)
       iseq = subsys%particles%n_els
       ALLOCATE( analyt_forces(iseq,3), stat=stat)
       iseq = 0
       particles => subsys%particles%els
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
       particles => subsys%particles%els
       Atom: DO ip = 1, subsys%particles%n_els
          iseq = iseq + 1
          Coord: DO k = 1, 3
             numer_energy = 0.0_dp
             std_value = particles(ip)%r(k)
             DO j = 1, 2
                particles(ip)%r(k) = std_value - (-1.0_dp)**j * Dx
                IF (.NOT. harris_flag) THEN
                   CALL force_env_calc_energy_force(force_env,calc_force=.TRUE.,&
                        consistent_energies=.TRUE., error=error)
                ELSE
                   CALL force_env_calc_energy_force(force_env,calc_force=.TRUE.,&
                        consistent_energies=.FALSE., error=error)
                END IF
                CALL force_env_get (force_env, potential_energy=numer_energy(j), error=error)
             END DO
             particles(ip)%r(k) = std_value
             numer_forces(ip,k) = - (numer_energy(1) - numer_energy(2) ) / (2.0_dp*Dx)
             IF (para_env%ionode) THEN
                WRITE(iw,'(A,I5,A,I5,4F15.9)')"ATOM NUMBER ::",iseq," COORD:",k,numer_energy,&
                     numer_forces(ip,k), analyt_forces(iseq,k)
             ENDIF
          END DO Coord
          !
          ! Check analytical forces Vs numerical forces
          !
          IF (para_env%ionode) THEN
             WRITE(iw,'(A,I5,6F15.9)')"ATOM NUMBER ::",iseq,analyt_forces(iseq,:), numer_forces(ip,:)
          ENDIF
       END DO Atom
       IF (para_env%ionode) THEN
          WRITE(iw,'(A,I5)')"SUMMARY DEBUG :: "
       ENDIF
       iseq = iseq - subsys%particles%n_els
       DO ip = 1, subsys%particles%n_els
          iseq = iseq + 1
          IF (para_env%ionode) THEN
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
          IF (para_env%ionode) THEN
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
       CALL cp_print_key_finished_output(iw,logger,root_section,&
            "GLOBAL%PROGRAM_RUN_INFO", error=error)
    END IF
100 FORMAT(I5,F15.9," ( ",F7.2," ) ",F15.9," ( ",F7.2," ) ",F15.9," ( ",F7.2," ) ")
  END SUBROUTINE qs_debug_energy_and_forces


END MODULE qs_debug
