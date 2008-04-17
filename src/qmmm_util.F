!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      09.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE qmmm_util
  USE cell_types,                      ONLY: cell_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             use_qmmm
  USE input_constants,                 ONLY: do_qmmm_wall_quadratic,&
                                             do_qmmm_wall_reflective
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: fist_subsys,&
                                             qs_subsys
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_util'
  PUBLIC :: apply_qmmm_walls_reflective,&
            apply_qmmm_walls,&
            apply_qmmm_translate

CONTAINS

! *****************************************************************************
!> \brief Apply QM quadratic walls in order to avoid QM atoms escaping from 
!>      the QM Box
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      02.2008 created 
!> \author Benjamin G Levine
! *****************************************************************************
  SUBROUTINE apply_qmmm_walls(force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_qmmm_walls', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iwall_type
    LOGICAL                                  :: explicit, failure
    TYPE(section_vals_type), POINTER         :: walls_section

    failure         = .FALSE.
    walls_section => section_vals_get_subs_vals(force_env%root_section,&
        "FORCE_EVAL%QMMM%WALLS",error=error)
    CALL section_vals_get(walls_section, explicit=explicit, error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(walls_section,"TYPE",i_val=iwall_type,error=error)
       SELECT CASE(iwall_type)
       CASE(do_qmmm_wall_quadratic)
          CALL apply_qmmm_walls_quadratic(force_env, walls_section, error)
       CASE(do_qmmm_wall_reflective)
          ! Do nothing.. reflective walls are applied directly in the integrator
       END SELECT
    ENDIF

  END SUBROUTINE apply_qmmm_walls

! *****************************************************************************
!> \brief Apply reflective QM walls in order to avoid QM atoms escaping from 
!>      the QM Box 
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      08.2007 created [tlaino] - Zurich University
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE apply_qmmm_walls_reflective(force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_qmmm_walls_reflective', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip, iwall_type, qm_index
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index
    LOGICAL                                  :: explicit, failure, is_x(2), &
                                                is_y(2), is_z(2)
    REAL(KIND=dp)                            :: skin
    REAL(KIND=dp), DIMENSION(3)              :: coord, qm_cell_diag
    TYPE(cell_type), POINTER                 :: mm_cell, qm_cell
    TYPE(cp_subsys_type), POINTER            :: subsys_mm, subsys_qm
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm
    TYPE(section_vals_type), POINTER         :: walls_section

    failure         = .FALSE.
    SELECT CASE(force_env%in_use)
    CASE(use_qmmm)
       NULLIFY(subsys_mm, subsys_qm, qm_atom_index,particles_mm, qm_cell, mm_cell,&
            walls_section)
       walls_section => section_vals_get_subs_vals(force_env%root_section,"FORCE_EVAL%QMMM%WALLS",error=error)
       CALL section_vals_get(walls_section, explicit=explicit, error=error)
       IF (explicit) THEN
          CALL section_vals_val_get(walls_section,"WALL_SKIN",r_val=skin,error=error)
          CALL section_vals_val_get(walls_section,"TYPE",i_val=iwall_type,error=error)
       ELSE
          skin = 0.0_dp
       END IF
       CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
       CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(force_env%qmmm_env),cp_failure_level,routineP,error,failure)
       CPPrecondition(force_env%qmmm_env%ref_count>0,cp_failure_level,routineP,error,failure)
   
       CALL force_env_get(force_env%sub_force_env(fist_subsys)%force_env,&
                          cell=mm_cell,subsys=subsys_mm,error=error)
       CALL force_env_get(force_env%sub_force_env(qs_subsys)%force_env,&
                          cell=qm_cell,subsys=subsys_qm,error=error)
       qm_atom_index   => force_env%qmmm_env%qm_atom_index
       CPPrecondition(ASSOCIATED(qm_atom_index),cp_failure_level,routineP,error,failure)
   
       qm_cell_diag = (/qm_cell%hmat(1,1),&
                        qm_cell%hmat(2,2),&
                        qm_cell%hmat(3,3)/)
       particles_mm => subsys_mm%particles%els
       DO ip=1,SIZE(qm_atom_index)
          qm_index = qm_atom_index(ip)
          coord = particles_mm(qm_index)%r
          IF (ANY(coord<skin).OR.ANY(coord>(qm_cell_diag-skin))) THEN
             IF (explicit) THEN
                IF (iwall_type==do_qmmm_wall_reflective) THEN
                   ! Apply Walls
                   is_x(1) = (coord(1)<skin)
                   is_x(2) = (coord(1)>(qm_cell_diag(1)-skin))
                   is_y(1) = (coord(2)<skin)                  
                   is_y(2) = (coord(2)>(qm_cell_diag(2)-skin))
                   is_z(1) = (coord(3)<skin)                  
                   is_z(2) = (coord(3)>(qm_cell_diag(3)-skin))
                   IF (ANY(is_x)) THEN
                      ! X coordinate
                      IF      (is_x(1)) THEN
                         particles_mm(qm_index)%v(1) =  ABS(particles_mm(qm_index)%v(1))
                      ELSE IF (is_x(2)) THEN
                         particles_mm(qm_index)%v(1) =  -ABS(particles_mm(qm_index)%v(1))
                      END IF
                   END IF
                   IF (ANY(is_y)) THEN
                      ! Y coordinate
                      IF      (is_y(1)) THEN
                         particles_mm(qm_index)%v(2) =  ABS(particles_mm(qm_index)%v(2))
                      ELSE IF (is_y(2)) THEN
                         particles_mm(qm_index)%v(2) =  -ABS(particles_mm(qm_index)%v(2))
                      END IF
                   END IF
                   IF (ANY(is_z)) THEN
                      ! Z coordinate
                      IF      (is_z(1)) THEN
                         particles_mm(qm_index)%v(3) =  ABS(particles_mm(qm_index)%v(3))
                      ELSE IF (is_z(2)) THEN
                         particles_mm(qm_index)%v(3) =  -ABS(particles_mm(qm_index)%v(3))
                      END IF
                   END IF
                ENDIF
             ELSE
                ! Otherwise print a warning and continue crossing cp2k's finger..
                CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                "One or few QM atoms are within the SKIN of the quantum box. Check your run "//&
                "and you may possibly consider: the activation of the QMMM WALLS "//&
                "around the QM box, switching ON the centering of the QM box or increase "//&
                "the size of the QM cell. CP2K CONTINUE but results could be meaningless. "//&
                CPSourceFileRef,&
                only_ionode=.TRUE.)
             END IF
          END IF 
       END DO 
    END SELECT

  END SUBROUTINE apply_qmmm_walls_reflective

! *****************************************************************************
!> \brief Apply QM quadratic walls in order to avoid QM atoms escaping from 
!>      the QM Box
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      02.2008 created 
!> \author Benjamin G Levine
! *****************************************************************************
  SUBROUTINE apply_qmmm_walls_quadratic(force_env, walls_section, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: walls_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_qmmm_walls_quadratic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip, qm_index
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index
    LOGICAL                                  :: failure, is_x(2), is_y(2), &
                                                is_z(2)
    REAL(KIND=dp)                            :: k, skin, wallenergy, wallforce
    REAL(KIND=dp), DIMENSION(3)              :: coord, qm_cell_diag
    TYPE(cell_type), POINTER                 :: mm_cell, qm_cell
    TYPE(cp_subsys_type), POINTER            :: subsys_mm, subsys_qm
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm
    TYPE(qs_energy_type), POINTER            :: energy

    CALL section_vals_val_get(walls_section,"WALL_SKIN",r_val=skin,error=error)
    CALL section_vals_val_get(walls_section,"K",r_val=k,error=error)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env%qmmm_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(force_env%qmmm_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CALL force_env_get(force_env%sub_force_env(fist_subsys)%force_env,&
                       cell=mm_cell,subsys=subsys_mm,error=error)
    CALL force_env_get(force_env%sub_force_env(qs_subsys)%force_env,&
                       cell=qm_cell,subsys=subsys_qm,error=error)
    qm_atom_index   => force_env%qmmm_env%qm_atom_index
    CPPrecondition(ASSOCIATED(qm_atom_index),cp_failure_level,routineP,error,failure)
  
    qm_cell_diag = (/qm_cell%hmat(1,1),&
                     qm_cell%hmat(2,2),&
                     qm_cell%hmat(3,3)/)
    particles_mm => subsys_mm%particles%els
    wallenergy=0.0_dp
    DO ip=1,SIZE(qm_atom_index)
       qm_index = qm_atom_index(ip)
       coord = particles_mm(qm_index)%r
       IF (ANY(coord<skin).OR.ANY(coord>(qm_cell_diag-skin))) THEN
          is_x(1) = (coord(1)<skin)
          is_x(2) = (coord(1)>(qm_cell_diag(1)-skin))
          is_y(1) = (coord(2)<skin)                  
          is_y(2) = (coord(2)>(qm_cell_diag(2)-skin))
          is_z(1) = (coord(3)<skin)                  
          is_z(2) = (coord(3)>(qm_cell_diag(3)-skin))
          IF (is_x(1)) THEN
            wallforce=2.0_dp*k*(skin-coord(1))
            particles_mm(qm_index)%f(1)=particles_mm(qm_index)%f(1)+&
                  wallforce
            wallenergy=wallenergy+wallforce*(skin-coord(1))*0.5_dp
          ENDIF
          IF (is_x(2)) THEN
            wallforce=2.0_dp*k*((qm_cell_diag(1)-skin)-coord(1))
            particles_mm(qm_index)%f(1)=particles_mm(qm_index)%f(1)+&
                  wallforce
            wallenergy=wallenergy+wallforce*((qm_cell_diag(1)-skin)-&
                  coord(1))*0.5_dp
          ENDIF
          IF (is_y(1)) THEN
            wallforce=2.0_dp*k*(skin-coord(2))
            particles_mm(qm_index)%f(2)=particles_mm(qm_index)%f(2)+&
                  wallforce
            wallenergy=wallenergy+wallforce*(skin-coord(2))*0.5_dp
          ENDIF
          IF (is_y(2)) THEN
            wallforce=2.0_dp*k*((qm_cell_diag(2)-skin)-coord(2))
            particles_mm(qm_index)%f(2)=particles_mm(qm_index)%f(2)+&
                  wallforce
            wallenergy=wallenergy+wallforce*((qm_cell_diag(2)-skin)-&
                  coord(2))*0.5_dp
          ENDIF
          IF (is_z(1)) THEN
            wallforce=2.0_dp*k*(skin-coord(3))
            particles_mm(qm_index)%f(3)=particles_mm(qm_index)%f(3)+&
                  wallforce
            wallenergy=wallenergy+wallforce*(skin-coord(3))*0.5_dp
          ENDIF
          IF (is_z(2)) THEN
            wallforce=2.0_dp*k*((qm_cell_diag(3)-skin)-coord(3))
            particles_mm(qm_index)%f(3)=particles_mm(qm_index)%f(3)+&
                  wallforce
            wallenergy=wallenergy+wallforce*((qm_cell_diag(3)-skin)-&
                  coord(3))*0.5_dp
          ENDIF
       ENDIF
    ENDDO

    CALL get_qs_env(qs_env=force_env%sub_force_env(qs_subsys)%&
      force_env%qs_env, energy=energy,error=error)
    energy%total = energy%total + wallenergy

  END SUBROUTINE

! *****************************************************************************
!> \brief Apply translation to the full system in order to center the QM
!>      system into the QM box
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      08.2007 created [tlaino] - Zurich University
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE apply_qmmm_translate(force_env,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_qmmm_translate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: max_coord, min_coord, transl_v
    TYPE(cell_type), POINTER                 :: mm_cell, qm_cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys_mm, subsys_qm
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm, particles_qm

    SELECT CASE(force_env%in_use)
    CASE(use_qmmm)
       min_coord       =  HUGE(0.0_dp)
       max_coord       = -HUGE(0.0_dp)
       failure         = .FALSE.
       logger => cp_error_get_logger(error)
       NULLIFY(subsys_mm, subsys_qm, qm_atom_index,particles_mm,particles_qm, qm_cell, mm_cell)

       CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
       CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(force_env%qmmm_env),cp_failure_level,routineP,error,failure)
       CPPrecondition(force_env%qmmm_env%ref_count>0,cp_failure_level,routineP,error,failure)

       CALL force_env_get(force_env%sub_force_env(fist_subsys)%force_env,&
                          cell=mm_cell,subsys=subsys_mm,error=error)
       CALL force_env_get(force_env%sub_force_env(qs_subsys)%force_env,&
                          cell=qm_cell,subsys=subsys_qm,error=error)
       qm_atom_index   => force_env%qmmm_env%qm_atom_index
       CPPrecondition(ASSOCIATED(qm_atom_index),cp_failure_level,routineP,error,failure)

       particles_qm => subsys_qm%particles%els
       particles_mm => subsys_mm%particles%els
       DO ip=1,SIZE(qm_atom_index)
          min_coord=MIN(min_coord,particles_mm(qm_atom_index(ip))%r)
          max_coord=MAX(max_coord,particles_mm(qm_atom_index(ip))%r)
       END DO
       IF (.NOT.force_env%qmmm_env%center_qm_subsys0) force_env%qmmm_env%do_translate = .FALSE.
       IF (force_env%qmmm_env%do_translate) THEN
          !
          ! The first time we always translate all the system in order
          ! to centre the QM system in the box.
          !
          transl_v = (max_coord + min_coord) / 2.0_dp
          transl_v(1) = transl_v(1) - qm_cell%hmat(1,1)/2.0_dp
          transl_v(2) = transl_v(2) - qm_cell%hmat(2,2)/2.0_dp
          transl_v(3) = transl_v(3) - qm_cell%hmat(3,3)/2.0_dp

          IF (ANY(force_env%qmmm_env%utrasl /= 1.0_dp)) THEN
             transl_v = REAL( FLOOR(transl_v/force_env%qmmm_env%utrasl),KIND=dp) *&
                  force_env%qmmm_env%utrasl
          END IF
          force_env%qmmm_env%transl_v = force_env%qmmm_env%transl_v + transl_v
          particles_mm => subsys_mm%particles%els
          DO ip=1,subsys_mm%particles%n_els
             particles_mm(ip)%r = particles_mm(ip)%r - transl_v
          END DO
          IF (logger%para_env%mepos==logger%para_env%source) &
               WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.FALSE.),fmt='(/1X,A)')&
               " Translating the system in order to center the QM fragment in the QM box."
          IF (.NOT.force_env%qmmm_env%center_qm_subsys) force_env%qmmm_env%do_translate = .FALSE.
       END IF
       particles_mm => subsys_mm%particles%els
       DO ip=1,SIZE(qm_atom_index)
          particles_qm(ip)%r=particles_mm(qm_atom_index(ip))%r
       END DO
    END SELECT

  END SUBROUTINE apply_qmmm_translate

END MODULE qmmm_util
