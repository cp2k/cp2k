!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/qs_numerical_pressure [1.0] *
!!
!!   NAME
!!     qs_numerical_pressure
!!
!!   FUNCTION
!!     Calculate stress tensor and pressure numerically
!!
!!   AUTHOR
!!     JCS
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [JCS]
!!             
!!
!!   SOURCE
!****************************************************************************
MODULE qs_numerical_pressure
  USE global_types,                    ONLY: global_environment_type
  USE cp_subsystem_types,              ONLY: cp_subsys_retain,&
                                             cp_subsys_get,&
                                             cp_subsystem_p_type
#include "cp_common_uses.h"
  
  USE qs_energy,                       ONLY: qs_energies
  USE force_env_types,                 ONLY: force_env_type, &
                                             force_env_release,&
                                             force_env_set, force_env_get
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type

  USE cell_types,                 ONLY: get_hinv, cell_type, get_cell,&
                                             cell_create, cell_clone,&
                                             cell_release
  USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env, &
                                             qs_env_setup
  USE qs_environment_types,            ONLY: get_qs_env,set_qs_env
  USE pw_env_types,                    ONLY: pw_env_type,pw_env_get
  USE qs_ks_types,                     ONLY: qs_ks_release, qs_ks_env_type
  USE qs_ks_methods,                   ONLY: qs_ks_create
  USE virial_types,                    ONLY: virial_type
  USE input_constants,                 ONLY: qs_debug_run
     

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_numerical_pressure'

  PUBLIC :: qs_calc_numerical_pressure

  CONTAINS

  SUBROUTINE qs_calc_numerical_pressure(force_env, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error    
    ! Local
    CHARACTER(len=*), PARAMETER :: routineN = 'qs_calc_numerical_pressure', &
      routineP = moduleN//':'//routineN    
    INTEGER :: iw, stat, iseq, isubsys, ip, k, j, l, n
    LOGICAL :: failure 
    REAL(kind=dp) :: numer_pressure = 0.0_dp
    REAL(kind=dp), PARAMETER :: Dx=0.001_dp
    REAL(kind=dp), DIMENSION(2) :: numer_energy
    REAL(kind=dp), DIMENSION(3) :: numer_stress  = 0.0_dp
    REAL(kind=dp), DIMENSION(:,:), POINTER     :: pos_ref
    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys
    TYPE(particle_type), POINTER, DIMENSION(:) :: particles
    TYPE(cell_type), POINTER                   :: cell, cell_ref
    TYPE(cell_type), POINTER                   :: cell_local
    TYPE(qs_ks_env_type), POINTER              :: new_ks_env
    TYPE(virial_type), POINTER                 :: virial
    TYPE(global_environment_type), POINTER     :: globenv


    NULLIFY(cell_local)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       !
       ! Copy atomic positions into pos_ref
       !
       CALL force_env_get(force_env,subsys=subsys,globenv=globenv,error=error)       
       iw = globenv%scr
 
       IF (globenv%ionode) THEN
       WRITE ( iw, '( /A,A/ )' ) ' **************************** ', &
           'NUMERICAL STRESS ********************************'
       END IF

!       CALL qs_energies(force_env%qs_env,globenv,&
!                        consistent_energies=.TRUE.,box_changed=.FALSE.)
       iseq = 0
       DO isubsys=1, SIZE(subsys)
          iseq = iseq + subsys(isubsys)%subsys%particles%n_els
       END DO
       ALLOCATE( pos_ref(iseq,3), stat=stat)
       iseq = 0
       DO isubsys=1, SIZE(subsys)
         particles => subsys(isubsys)%subsys%particles%els
         DO ip = 1, subsys(isubsys)%subsys%particles%n_els
            iseq = iseq + 1
            pos_ref(iseq,:) = particles(ip)%r
         END DO
       END DO
       CALL force_env_get(force_env,cell=cell,cell_ref=cell_ref,&
                          virial=virial, error=error)
       CALL cell_create(cell_local)
       CALL cell_clone(cell,cell_local)
       !
       ! First change box.
       !

       DO ip = 1, 1
         DO k = 1, 2
           cell%hmat(ip,ip) = cell_local%hmat(ip,ip)  - (-1.0_dp)**k * Dx
           CALL get_hinv( cell )
           iseq = 0
           !          
           ! Then scale positions.   
           !
           DO  isubsys = 1, SIZE(subsys)
             particles => subsys(isubsys)%subsys%particles%els
           DO j = 1, subsys(isubsys)%subsys%particles%n_els
             iseq = iseq + 1
             DO l = 1, 3
             IF (l == ip) THEN
                   particles(j)%r(l) = pos_ref(iseq,l) * cell%hmat(l,l) / cell_local%hmat(l,l)
                 ELSE
                   particles(j)%r(l) = pos_ref(iseq,l)
                 END IF
              END DO
             END DO
         END DO
           !
           ! Since the box has changed, rebuild grids, i.e. pw_env and ks_env
           !
           force_env%box_changed = .TRUE.
           force_env%qs_env%use_ref_cell=.TRUE.
           CALL qs_energies(force_env%qs_env,globenv,&
                            consistent_energies=.TRUE.,box_changed=.TRUE.,error=error)
           numer_energy(k) = force_env%qs_env%energy%total
! reset cell
           cell%hmat(ip,ip) = cell_local%hmat(ip,ip)
         END DO
       CALL get_hinv( cell )
         numer_stress(ip) = (numer_energy(1) - numer_energy(2) ) / (2.0_dp*Dx)                
         virial%pv_virial(ip,ip) = (-1.0_dp)*numer_stress(ip) * cell_local%hmat(ip,ip)
         IF (globenv%run_type_id.EQ.qs_debug_run) THEN
           IF (globenv%ionode) THEN
              WRITE(iw,'(/A/A,I1,A,I1,A,F15.9)')"----- NUMERICAL STRESS -----",&
                "pv_virial(",ip,",",ip,"): ",virial%pv_virial(ip,ip)
           END IF
         END IF
       END DO
       !
       ! Reset postions and rebuild environment.
       ! 
       iseq = 0
       DO isubsys=1, SIZE(subsys)
         particles => subsys(isubsys)%subsys%particles%els
         DO ip = 1, subsys(isubsys)%subsys%particles%n_els
            iseq = iseq + 1
            particles(ip)%r = pos_ref(iseq,:)
         END DO
       END DO
       CALL get_hinv ( cell )
       force_env%box_changed = .TRUE.
       force_env%qs_env%use_ref_cell=.TRUE.
!       CALL qs_env_rebuild_pw_env(force_env%qs_env,error=error)
!       CALL qs_ks_create ( new_ks_env, force_env%qs_env, error=error )
!       CALL set_qs_env(force_env%qs_env, ks_env=new_ks_env,error=error)
!       CALL qs_ks_release ( new_ks_env, error=error )
!       force_env%box_changed = .FALSE.
!      CALL qs_energies(force_env%qs_env,globenv,&
!                       consistent_energies=.TRUE.,box_changed=.TRUE.)

       IF (globenv%ionode) THEN
       WRITE ( iw, '( /A,A/ )' ) ' **************************** ', &
           'NUMERICAL STRESS END *****************************'
       END IF

       IF (ASSOCIATED( pos_ref )) DEALLOCATE( pos_ref )
       IF (ASSOCIATED( cell_local )) CALL cell_release( cell_local )
    END IF
  END SUBROUTINE qs_calc_numerical_pressure


END MODULE qs_numerical_pressure
