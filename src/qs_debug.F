!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
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
  USE global_types,                    ONLY: global_environment_type
  USE cp_subsystem_types,              ONLY: cp_subsys_retain,&
                                             cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_type, &
                                             force_env_release,&
                                             force_env_set, force_env_get
  USE input_section_types,             ONLY: section_vals_get_subs_vals, &
                                             section_vals_val_get, &
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
#include "cp_common_uses.h"

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
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error    
    ! Local
    CHARACTER(len=*), PARAMETER :: routineN = 'qs_debug_energy_and_forces', &
      routineP = moduleN//':'//routineN    
    INTEGER :: iw, stat, iseq, isubsys, ip, k, j
    LOGICAL :: failure, harris_flag
    REAL(kind=dp) :: std_value
    REAL(kind=dp), PARAMETER :: Dx=0.001_dp
    REAL(kind=dp), DIMENSION(2) :: numer_energy
    REAL(kind=dp), DIMENSION(3) :: Err, my_maxerr
    REAL(kind=dp), DIMENSION(:,:), POINTER     :: analyt_forces, numer_forces
    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys
    TYPE(particle_type), POINTER, DIMENSION(:) :: particles
    TYPE(section_vals_type), POINTER           :: harris_section

    failure = .FALSE.
    NULLIFY(analyt_forces, numer_forces, subsys, particles)
    IF (.NOT.failure) THEN
       harris_section => section_vals_get_subs_vals(globenv%input_file, &
                                                    "FORCE_EVAL%DFT%QS%HARRIS")
       CALL section_vals_val_get(harris_section, "ACTIVATE", &
                                 l_val=harris_flag, error=error)

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
                   IF (.NOT. harris_flag) THEN
                     CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,&
                                                      consistent_energies=.TRUE., error=error)
                   ELSE
                     CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,&
                                                      consistent_energies=.FALSE., error=error)
                   END IF
                   CALL force_env_get (force_env, potential_energy=numer_energy(j), error=error)
                END DO
                particles(ip)%r(k) = std_value
                numer_forces(ip,k) = - (numer_energy(1) - numer_energy(2) ) / (2.0_dp*Dx)
                IF (globenv%ionode) THEN
                   WRITE(iw,'(A,I5,A,I5,4F15.9)')"ATOM NUMBER ::",iseq," COORD:",k,numer_energy,&
                     numer_forces(ip,k), analyt_forces(iseq,k)
                ENDIF
             END DO Coord
             !
             ! Check analytical forces Vs numerical forces
             !             
             IF (globenv%ionode) THEN
                WRITE(iw,'(A,I5,6F15.9)')"ATOM NUMBER ::",iseq,analyt_forces(iseq,:), numer_forces(ip,:)
             ENDIF
          END DO Atom
          IF (globenv%ionode) THEN
              WRITE(iw,'(A,I5)')"SUMMARY DEBUG :: SUBSYS NUMBER :",isubsys
          ENDIF
          iseq = iseq - subsys(isubsys)%subsys%particles%n_els
          DO ip = 1, subsys(isubsys)%subsys%particles%n_els
             iseq = iseq + 1 
             IF (globenv%ionode) THEN
                 WRITE(iw,'(A,I5,9F12.6)')"ATOM NUMBER ::",iseq,analyt_forces(iseq,:), numer_forces(ip,:),&
                      analyt_forces(iseq,:) - numer_forces(ip,:)   
             ENDIF
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
             IF (globenv%ionode) THEN
                 WRITE(iw,100)ip,Analyt_Forces(Ip,1),Err(1),&
                            Analyt_Forces(Ip,2),Err(2),&
                            Analyt_Forces(Ip,3),Err(3)
             ENDIF
             my_MaxErr = MaxErr
             IF (ABS(Analyt_Forces(Ip,1)) <= 0.0001_dp ) my_MaxErr(1)=my_MaxErr(1)*5.0_dp
             IF (ABS(Analyt_Forces(Ip,2)) <= 0.0001_dp ) my_MaxErr(1)=my_MaxErr(2)*5.0_dp
             IF (ABS(Analyt_Forces(Ip,3)) <= 0.0001_dp ) my_MaxErr(1)=my_MaxErr(3)*5.0_dp
             ! ? does the whole of this routine actually work in parallel 
             IF (globenv%ionode) THEN
                WRITE(iw,*)my_MaxErr
             ENDIF
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
