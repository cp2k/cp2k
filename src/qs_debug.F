!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_debug [1.0] *
!!
!!   NAME
!!     qs_debug
!!
!!   FUNCTION
!!     Debug energy and derivatives w.r.t. finite differences
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     12.2004 created [tlaino]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_debug
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_assertion_failed
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level 
  USE global_types,                    ONLY: global_environment_type,&
                                             globenv_create,&
                                             globenv_retain,&
                                             globenv_release
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_type, &
                                             force_env_release,&
                                             force_env_set, force_env_get
  USE cp_subsystem_types,              ONLY: cp_subsys_retain,&
                                             cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_type, &
                                             force_env_release,&
                                             force_env_set, force_env_get
  USE kinds,                           only: dp
  USE particle_types,             ONLY: particle_type

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_debug'

  PUBLIC :: qs_debug_energy_and_forces
  REAL(KIND=dp), PRIVATE, PARAMETER :: MaxErr = 1.0_dp

CONTAINS

  SUBROUTINE qs_debug_energy_and_forces(force_env, globenv, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER  :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error    
    ! Local
    CHARACTER(len=*), PARAMETER :: routineN = 'qs_debug_energy_and_forces', &
      routineP = moduleN//':'//routineN    
    INTEGER :: iw, stat, iseq, isubsys, ip, k, j
    LOGICAL :: failure
    REAL(kind=dp) :: std_value
    REAL(kind=dp), PARAMETER :: Dx=0.001_dp
    REAL(kind=dp), DIMENSION(2) :: numer_energy
    REAL(kind=dp), DIMENSION(3) :: Err, my_maxerr
    REAL(kind=dp), DIMENSION(:,:), POINTER     :: analyt_forces, numer_forces
    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys
    TYPE(particle_type), POINTER, DIMENSION(:) :: particles

    failure = .FALSE.
    NULLIFY(analyt_forces, numer_forces, subsys, particles)
    IF (.NOT.failure) THEN
       iw = globenv%scr
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
       iseq = 0
       DO isubsys=1, SIZE(subsys)
          iseq = iseq + subsys(isubsys)%subsys%particles%n_els
       END DO
       ALLOCATE( analyt_forces(iseq,3), stat=stat)
       iseq = 0
       DO isubsys=1, SIZE(subsys)
          particles => subsys(isubsys)%subsys%particles%els
          DO ip = 1, subsys(isubsys)%subsys%particles%n_els
             iseq = iseq + 1
             analyt_forces(iseq,:) = particles(ip)%f
          END DO
       END DO
       !
       ! Loop on atoms and coordinates
       !
       iseq = 0
       DO isubsys=1, SIZE(subsys)
          IF (ASSOCIATED( numer_forces)) DEALLOCATE( numer_forces)
          ALLOCATE(  numer_forces(subsys(isubsys)%subsys%particles%n_els,3), stat=stat)          
          particles => subsys(isubsys)%subsys%particles%els
          Atom: DO ip = 1, subsys(isubsys)%subsys%particles%n_els
             iseq = iseq + 1
             Coord: DO k = 1, 3
                numer_energy = 0.0_dp
                std_value = particles(ip)%r(k)
                DO j = 1, 2 
                   particles(ip)%r(k) = std_value - (-1.0_dp)**j * Dx                   
                   CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,&
                        error=error)
                   CALL force_env_get (force_env, potential_energy=numer_energy(j), error=error)
                END DO
                particles(ip)%r(k) = std_value
                numer_forces(ip,k) = - (numer_energy(1) - numer_energy(2) ) / (2.0_dp*Dx)
                WRITE(iw,'(A,I5,A,I5,4F15.9)')"ATOM NUMBER ::",iseq," COORD:",k,numer_energy,&
                     numer_forces(ip,k), analyt_forces(iseq,k)
             END DO Coord
             !
             ! Check analytical forces Vs numerical forces
             !             
             WRITE(iw,'(A,I5,6F15.9)')"ATOM NUMBER ::",iseq,analyt_forces(iseq,:), numer_forces(ip,:)
          END DO Atom
          WRITE(iw,'(A,I5)')"SUMMARY DEBUG :: SUBSYS NUMBER :",isubsys
          iseq = iseq - subsys(isubsys)%subsys%particles%n_els
          DO ip = 1, subsys(isubsys)%subsys%particles%n_els
             iseq = iseq + 1
             WRITE(iw,'(A,I5,9F12.6)')"ATOM NUMBER ::",iseq,analyt_forces(iseq,:), numer_forces(ip,:),&
                  analyt_forces(iseq,:) - numer_forces(ip,:)   
          END DO
          !
          ! Runtime check...
          !
          DO ip =  1, subsys(isubsys)%subsys%particles%n_els
             Err = 0.0_dp
             DO K = 1, 3
                IF (ABS(numer_forces(ip,K)) >= 1.0E-6_dp) THEN 
                   Err(K) = (analyt_forces(ip,K)-numer_forces(ip,K))/numer_forces(ip,K)*100.0_dp
                END IF
             END DO
             WRITE(*,100)ip,Analyt_Forces(Ip,1),Err(1),&
                            Analyt_Forces(Ip,2),Err(2),&
                            Analyt_Forces(Ip,3),Err(3)
             my_MaxErr = MaxErr
             IF (ABS(Analyt_Forces(Ip,1)) <= 0.0001_dp ) my_MaxErr(1)=my_MaxErr(1)*5.0_dp
             IF (ABS(Analyt_Forces(Ip,2)) <= 0.0001_dp ) my_MaxErr(1)=my_MaxErr(2)*5.0_dp
             IF (ABS(Analyt_Forces(Ip,3)) <= 0.0001_dp ) my_MaxErr(1)=my_MaxErr(3)*5.0_dp
             WRITE(*,*)my_MaxErr
             CPPostcondition(ABS(Err(1))<=my_MaxErr(1),cp_failure_level,routineP,error,failure)
             CPPostcondition(ABS(Err(2))<=my_MaxErr(2),cp_failure_level,routineP,error,failure)
             CPPostcondition(ABS(Err(3))<=my_MaxErr(3),cp_failure_level,routineP,error,failure)
          END DO
       END DO
       IF (ASSOCIATED(analyt_forces)) DEALLOCATE(analyt_forces)
       IF (ASSOCIATED( numer_forces)) DEALLOCATE( numer_forces)
    END IF
100 FORMAT(I5,F15.9," ( ",F7.2," ) ",F15.9," ( ",F7.2," ) ",F15.9," ( ",F7.2," ) ")
  END SUBROUTINE qs_debug_energy_and_forces


END MODULE qs_debug
