!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      09.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE qmmm_util
  USE cell_methods,                    ONLY: cell_copy,&
                                             cell_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             use_qmmm
  USE input_constants,                 ONLY: do_qmmm_wall_none,&
                                             do_qmmm_wall_quadratic,&
                                             do_qmmm_wall_reflective
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE particle_types,                  ONLY: particle_type,&
                                             write_fist_particle_coordinates,&
                                             write_qs_particle_coordinates
  USE qmmm_types,                      ONLY: fist_subsys,&
                                             force_mixing_core_subsys,&
                                             force_mixing_extended_subsys,&
                                             primary_subsys,&
                                             qs_subsys
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_util'
  PUBLIC :: apply_qmmm_walls_reflective,&
            apply_qmmm_walls,&
            apply_qmmm_translate,&
            apply_qmmm_translate_low,&
            apply_qmmm_wrap,&
            apply_qmmm_unwrap,&
            spherical_cutoff_factor,&
            qmmm_force_mixing_active

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
    REAL(KIND=dp), DIMENSION(3)              :: coord, qm_cell_diag, skin
    REAL(KIND=dp), DIMENSION(:), POINTER     :: list
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
          NULLIFY(list)
          CALL section_vals_val_get(walls_section,"WALL_SKIN",r_vals=list,error=error)
          CALL section_vals_val_get(walls_section,"TYPE",i_val=iwall_type,error=error)
          skin(:) = list(:)
       ELSE
          ![NB]
          iwall_type=do_qmmm_wall_reflective
          skin(:) = 0.0_dp
       END IF
       CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
       CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(force_env%qmmm_env),cp_failure_level,routineP,error,failure)
       CPPrecondition(force_env%qmmm_env%ref_count>0,cp_failure_level,routineP,error,failure)

       IF (force_env%qmmm_env%do_force_mixing) THEN
          CALL cp_assert(iwall_type == do_qmmm_wall_none, &
            cp_warning_level,cp_assertion_failed,routineP,&
            "Reflective walls for QM/MM are not implemented (or useful) when "//&
            "force mixing is active.  Skipping!"//CPSourceFileRef)
          RETURN
       END IF
   
       CALL force_env_get(force_env%sub_force_env(primary_subsys)%force_env%sub_force_env(fist_subsys)%force_env,&
                          cell=mm_cell,subsys=subsys_mm,error=error)
       CALL force_env_get(force_env%sub_force_env(primary_subsys)%force_env%sub_force_env(qs_subsys)%force_env,&
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
                   is_x(1) = (coord(1)<skin(1))
                   is_x(2) = (coord(1)>(qm_cell_diag(1)-skin(1)))
                   is_y(1) = (coord(2)<skin(2))
                   is_y(2) = (coord(2)>(qm_cell_diag(2)-skin(2)))
                   is_z(1) = (coord(3)<skin(3))
                   is_z(2) = (coord(3)>(qm_cell_diag(3)-skin(3)))
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
    REAL(KIND=dp)                            :: k, wallenergy, wallforce
    REAL(KIND=dp), DIMENSION(3)              :: coord, qm_cell_diag, skin
    REAL(KIND=dp), DIMENSION(:), POINTER     :: list
    TYPE(cell_type), POINTER                 :: mm_cell, qm_cell
    TYPE(cp_subsys_type), POINTER            :: subsys_mm, subsys_qm
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm
    TYPE(qs_energy_type), POINTER            :: energy

    NULLIFY(list)
    CALL section_vals_val_get(walls_section,"WALL_SKIN",r_vals=list,error=error)
    CALL section_vals_val_get(walls_section,"K",r_val=k,error=error)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env%qmmm_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(force_env%qmmm_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (force_env%qmmm_env%do_force_mixing) THEN
       CALL cp_assert(.FALSE., &
          cp_warning_level,cp_assertion_failed,routineP,&
          "Quadratic walls for QM/MM are not implemented (or useful), when "//&
          "force mixing is active.  Skipping!"//CPSourceFileRef)
       RETURN
    END IF

    CALL force_env_get(force_env%sub_force_env(fist_subsys)%force_env,&
                       cell=mm_cell,subsys=subsys_mm,error=error)
    CALL force_env_get(force_env%sub_force_env(qs_subsys)%force_env,&
                       cell=qm_cell,subsys=subsys_qm,error=error)

    qm_atom_index   => force_env%qmmm_env%qm_atom_index
    CPPrecondition(ASSOCIATED(qm_atom_index),cp_failure_level,routineP,error,failure)

    skin(:) = list(:)

    qm_cell_diag = (/qm_cell%hmat(1,1),&
                     qm_cell%hmat(2,2),&
                     qm_cell%hmat(3,3)/)
    particles_mm => subsys_mm%particles%els
    wallenergy=0.0_dp
    DO ip=1,SIZE(qm_atom_index)
       qm_index = qm_atom_index(ip)
       coord = particles_mm(qm_index)%r
       IF (ANY(coord<skin).OR.ANY(coord>(qm_cell_diag-skin))) THEN
          is_x(1) = (coord(1)<skin(1))
          is_x(2) = (coord(1)>(qm_cell_diag(1)-skin(1)))
          is_y(1) = (coord(2)<skin(2))
          is_y(2) = (coord(2)>(qm_cell_diag(2)-skin(2)))
          is_z(1) = (coord(3)<skin(3))
          is_z(2) = (coord(3)>(qm_cell_diag(3)-skin(3)))
          IF (is_x(1)) THEN
            wallforce=2.0_dp*k*(skin(1)-coord(1))
            particles_mm(qm_index)%f(1)=particles_mm(qm_index)%f(1)+&
                  wallforce
            wallenergy=wallenergy+wallforce*(skin(1)-coord(1))*0.5_dp
          ENDIF
          IF (is_x(2)) THEN
            wallforce=2.0_dp*k*((qm_cell_diag(1)-skin(1))-coord(1))
            particles_mm(qm_index)%f(1)=particles_mm(qm_index)%f(1)+&
                  wallforce
            wallenergy=wallenergy+wallforce*((qm_cell_diag(1)-skin(1))-&
                  coord(1))*0.5_dp
          ENDIF
          IF (is_y(1)) THEN
            wallforce=2.0_dp*k*(skin(2)-coord(2))
            particles_mm(qm_index)%f(2)=particles_mm(qm_index)%f(2)+&
                  wallforce
            wallenergy=wallenergy+wallforce*(skin(2)-coord(2))*0.5_dp
          ENDIF
          IF (is_y(2)) THEN
            wallforce=2.0_dp*k*((qm_cell_diag(2)-skin(2))-coord(2))
            particles_mm(qm_index)%f(2)=particles_mm(qm_index)%f(2)+&
                  wallforce
            wallenergy=wallenergy+wallforce*((qm_cell_diag(2)-skin(2))-&
                  coord(2))*0.5_dp
          ENDIF
          IF (is_z(1)) THEN
            wallforce=2.0_dp*k*(skin(3)-coord(3))
            particles_mm(qm_index)%f(3)=particles_mm(qm_index)%f(3)+&
                  wallforce
            wallenergy=wallenergy+wallforce*(skin(3)-coord(3))*0.5_dp
          ENDIF
          IF (is_z(2)) THEN
            wallforce=2.0_dp*k*((qm_cell_diag(3)-skin(3))-coord(3))
            particles_mm(qm_index)%f(3)=particles_mm(qm_index)%f(3)+&
                  wallforce
            wallenergy=wallenergy+wallforce*((qm_cell_diag(3)-skin(3))-&
                  coord(3))*0.5_dp
          ENDIF
       ENDIF
    ENDDO

    CALL get_qs_env(qs_env=force_env%sub_force_env(qs_subsys)%force_env%qs_env, energy=energy,error=error)
    energy%total = energy%total + wallenergy

  END SUBROUTINE apply_qmmm_walls_quadratic

  ! wrap positions (with mm periodicity)
  SUBROUTINE apply_qmmm_wrap(subsys_mm, mm_cell, subsys_qm, qm_atom_index, saved_pos, error)
    TYPE(cp_subsys_type), POINTER            :: subsys_mm
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys_qm
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: qm_atom_index
    REAL(dp), ALLOCATABLE                    :: saved_pos(:,:)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_qmmm_wrap', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i_dim, ip
    REAL(dp)                                 :: r_lat(3)

    ALLOCATE(saved_pos(3,subsys_mm%particles%n_els))

    DO ip=1, subsys_mm%particles%n_els
       saved_pos(1:3,ip) = subsys_mm%particles%els(ip)%r(1:3)
       r_lat = MATMUL(mm_cell%h_inv, subsys_mm%particles%els(ip)%r)
       DO i_dim=1, 3
          IF (mm_cell%perd(i_dim) /= 1) THEN
            r_lat(i_dim) = 0.0_dp
          END IF
       END DO
       subsys_mm%particles%els(ip)%r = subsys_mm%particles%els(ip)%r - MATMUL(mm_cell%hmat, FLOOR(r_lat))
    END DO

    IF (PRESENT(subsys_qm) .AND. PRESENT(qm_atom_index)) THEN
       DO ip=1,SIZE(qm_atom_index)
          subsys_qm%particles%els(ip)%r=subsys_mm%particles%els(qm_atom_index(ip))%r
       END DO
    END IF
  END SUBROUTINE apply_qmmm_wrap

  SUBROUTINE apply_qmmm_unwrap(subsys_mm, subsys_qm, qm_atom_index, saved_pos, error)
    TYPE(cp_subsys_type), POINTER            :: subsys_mm
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys_qm
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: qm_atom_index
    REAL(dp), ALLOCATABLE                    :: saved_pos(:,:)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_qmmm_unwrap', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip

    DO ip=1, subsys_mm%particles%n_els
       subsys_mm%particles%els(ip)%r(1:3) = saved_pos(1:3,ip)
    END DO

    IF (PRESENT(subsys_qm) .AND. PRESENT(qm_atom_index)) THEN
       DO ip=1,SIZE(qm_atom_index)
          subsys_qm%particles%els(ip)%r=subsys_mm%particles%els(qm_atom_index(ip))%r
       END DO
    END IF

    DEALLOCATE(saved_pos)
  END SUBROUTINE apply_qmmm_unwrap

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
    REAL(dp), POINTER                        :: charges(:)
    TYPE(cell_type), POINTER                 :: cell_core, cell_extended
    TYPE(cp_subsys_type), POINTER            :: subsys_core, subsys_extended, &
                                                subsys_fist, subsys_qs
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_core, &
                                                particles_extended, &
                                                particles_mm, particles_qm
    TYPE(section_vals_type), POINTER         :: subsys_section

    SELECT CASE(force_env%in_use)

    CASE(use_qmmm)
      IF (SIZE(force_env%sub_force_env) == 1) THEN
        CALL apply_qmmm_translate_low(force_env%sub_force_env(primary_subsys)%force_env,error)
      ELSE IF (SIZE(force_env%sub_force_env) == 2) THEN
        ! want to center extended, and make core consistent with that
        CALL apply_qmmm_translate_low(force_env%sub_force_env(force_mixing_extended_subsys)%force_env,error)

        ! translate core fist particles
        CALL force_env_get(force_env%sub_force_env(force_mixing_extended_subsys)%force_env,&
                          subsys=subsys_extended,cell=cell_extended,error=error)
        CALL force_env_get(force_env%sub_force_env(force_mixing_core_subsys)%force_env,&
                          subsys=subsys_core,cell=cell_core,error=error)
        particles_extended => subsys_extended%particles%els
        particles_core => subsys_core%particles%els
        DO ip=1,SIZE(particles_extended)
           particles_core(ip)%r = particles_extended(ip)%r
        END DO
        CALL cell_copy(cell_extended, cell_core, error)

        ! also update core QM particles (like in apply_qmmm_translate_low)
        CALL force_env_get(force_env%sub_force_env(force_mixing_core_subsys)%force_env%sub_force_env(fist_subsys)%force_env,&
                         subsys=subsys_fist,error=error)
        CALL force_env_get(force_env%sub_force_env(force_mixing_core_subsys)%force_env%sub_force_env(qs_subsys)%force_env,&
                         subsys=subsys_qs,error=error)
        qm_atom_index   => force_env%sub_force_env(force_mixing_core_subsys)%force_env%qmmm_env%qm_atom_index
        particles_qm => subsys_qs%particles%els
        particles_mm => subsys_fist%particles%els
        CPPrecondition(ASSOCIATED(qm_atom_index),cp_failure_level,routineP,error,failure)
        DO ip=1,SIZE(qm_atom_index)
           particles_qm(ip)%r=particles_mm(qm_atom_index(ip))%r
        END DO

        IF (debug_this_module) THEN
          CALL force_env_get(force_env%sub_force_env(force_mixing_core_subsys)%force_env%sub_force_env(qs_subsys)%force_env,&
                            subsys=subsys_core,cell=cell_core,error=error)
          particles_core => subsys_core%particles%els
          subsys_section => section_vals_get_subs_vals(force_env%force_env_section, &
                                                              "SUBSYS",error=error)
          CALL write_qs_particle_coordinates(particles_core,subsys_section,"QM/MM core calc first QM, then MM (0 charges)",error)
          CALL force_env_get(force_env%sub_force_env(force_mixing_core_subsys)%force_env%sub_force_env(fist_subsys)%force_env,&
                            subsys=subsys_core,cell=cell_core,error=error)
          particles_core => subsys_core%particles%els
          ALLOCATE(charges(SIZE(particles_core)))
          charges = 0.0_dp
          CALL write_fist_particle_coordinates(particles_core,subsys_section,charges,error)
          DEALLOCATE(charges)
        ENDIF

      ELSE
         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
      ENDIF
    END SELECT
  END SUBROUTINE apply_qmmm_translate

! *****************************************************************************
!> \brief Apply translation to the full system in order to center the QM
!>      system into the QM box
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      08.2007 created [tlaino] - Zurich University
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE apply_qmmm_translate_low(force_env,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_qmmm_translate_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bigger_ip, i_dim, ip, max_ip, &
                                                min_ip, smaller_ip, tmp_ip, &
                                                unit_nr
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index
    LOGICAL                                  :: failure
    LOGICAL, ALLOCATABLE                     :: avoid(:)
    REAL(DP) :: bigger_lat_dv, center_p(3), lat_dv, lat_dv3(3), lat_min(3), &
      lat_p(3), max_coord_lat(3), min_coord_lat(3), smaller_lat_dv
    REAL(DP), POINTER                        :: charges(:)
    REAL(KIND=dp), DIMENSION(3)              :: max_coord, min_coord, transl_v
    TYPE(cell_type), POINTER                 :: mm_cell, qm_cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys_mm, subsys_qm
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm, particles_qm
    TYPE(section_vals_type), POINTER         :: subsys_section

    SELECT CASE(force_env%in_use)
    CASE(use_qmmm)
       failure         = .FALSE.
       logger => cp_error_get_logger(error)
       NULLIFY(subsys_mm, subsys_qm, qm_atom_index, particles_mm, particles_qm,&
               subsys_section, qm_cell, mm_cell)

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
       IF (.NOT.force_env%qmmm_env%center_qm_subsys0) force_env%qmmm_env%do_translate = .FALSE.
       IF (force_env%qmmm_env%do_translate) THEN
          IF (.NOT. force_env%qmmm_env%center_qm_subsys_pbc_aware) THEN
             ! naive coordinate based min-max
             min_coord       =  HUGE(0.0_dp)
             max_coord       = -HUGE(0.0_dp)
             DO ip=1,SIZE(qm_atom_index)
                min_coord=MIN(min_coord,particles_mm(qm_atom_index(ip))%r)
                max_coord=MAX(max_coord,particles_mm(qm_atom_index(ip))%r)
             END DO
          ELSE
             !! periodic based min max (uses complex number based mean)
             center_p = qmmm_pbc_aware_mean(particles_mm, mm_cell, qm_atom_index)
             ALLOCATE(avoid(SIZE(qm_atom_index)))
             DO i_dim=1, 3
                IF (mm_cell%perd(i_dim) /= 1) THEN
                  ! find absolute min and max positions (along i_dim direction) in lattice coordinates
                  min_coord_lat(i_dim) = HUGE(0.0_dp)
                  max_coord_lat(i_dim) = -HUGE(0.0_dp)
                  DO ip=1, SIZE(qm_atom_index)
                    lat_p = MATMUL(mm_cell%h_inv, particles_mm(qm_atom_index(ip))%r)
                    min_coord_lat(i_dim)=MIN(lat_p(i_dim), min_coord_lat(i_dim))
                    max_coord_lat(i_dim)=MAX(lat_p(i_dim), max_coord_lat(i_dim))
                  END DO
                ELSE
                  ! find min_ip closest to (pbc-aware) mean pos
                  avoid = .FALSE.
                  min_ip = qmmm_find_closest(particles_mm, mm_cell, qm_atom_index, avoid, center_p, i_dim, 0)
                  avoid(min_ip) = .TRUE.
                  ! find max_ip closest to min_ip
                  max_ip = qmmm_find_closest(particles_mm, mm_cell, qm_atom_index, avoid, &
                           particles_mm(qm_atom_index(min_ip))%r, i_dim, 0, lat_dv)
                  avoid(max_ip) = .TRUE.
                  ! switch min and max if necessary
                  IF (lat_dv < 0.0) THEN
                    tmp_ip = min_ip
                    min_ip = max_ip
                    max_ip = tmp_ip
                  ENDIF
                  ! loop over all other atoms
                  DO WHILE (.NOT. ALL(avoid))
                    ! find smaller below min, bigger after max
                    smaller_ip = qmmm_find_closest(particles_mm, mm_cell, qm_atom_index,&
                                 avoid, particles_mm(qm_atom_index(min_ip))%r, i_dim, -1, smaller_lat_dv)
                    bigger_ip = qmmm_find_closest(particles_mm, mm_cell, qm_atom_index,&
                                 avoid, particles_mm(qm_atom_index(max_ip))%r, i_dim, 1, bigger_lat_dv)
                    ! move min or max, not both
                    IF (ABS(smaller_lat_dv) < ABS(bigger_lat_dv)) THEN
                       min_ip = smaller_ip
                       avoid(min_ip) = .TRUE.
                    ELSE
                       max_ip = bigger_ip
                       avoid(max_ip) = .TRUE.
                    END IF
                  END DO
                  ! find min and max coordinates in lattice positions (i_dim ! only)
                  lat_dv3 = qmmm_lat_dv(mm_cell, particles_mm(qm_atom_index(min_ip))%r, particles_mm(qm_atom_index(max_ip))%r)
                  IF (lat_dv3(i_dim) < 0.0_dp) lat_dv3(i_dim) = lat_dv3(i_dim) + 1.0_dp
                  lat_min = MATMUL(mm_cell%h_inv, particles_mm(qm_atom_index(min_ip))%r)
                  min_coord_lat(i_dim) = lat_min(i_dim)
                  max_coord_lat(i_dim) = lat_min(i_dim) + lat_dv3(i_dim)
                END IF ! periodic
             END DO ! i_dim
             ! min and max coordinates from lattice positions to Cartesian
             min_coord = MATMUL(mm_cell%hmat, min_coord_lat)
             max_coord = MATMUL(mm_cell%hmat, max_coord_lat)
             DEALLOCATE(avoid)
          END IF ! pbc aware center
          transl_v = (max_coord + min_coord) / 2.0_dp

          !
          ! The first time we always translate all the system in order
          ! to centre the QM system in the box.
          !
          transl_v(:) = transl_v(:) - SUM(qm_cell%hmat, 2)/2.0_dp

          IF (ANY(force_env%qmmm_env%utrasl /= 1.0_dp)) THEN
             transl_v = REAL( FLOOR(transl_v/force_env%qmmm_env%utrasl),KIND=dp) *&
                  force_env%qmmm_env%utrasl
          END IF
          force_env%qmmm_env%transl_v = force_env%qmmm_env%transl_v + transl_v
          particles_mm => subsys_mm%particles%els
          DO ip=1,subsys_mm%particles%n_els
             particles_mm(ip)%r = particles_mm(ip)%r - transl_v
          END DO
          unit_nr=cp_logger_get_default_io_unit(logger)
          IF (unit_nr>0) WRITE (unit=unit_nr,fmt='(/1X,A)')&
               " Translating the system in order to center the QM fragment in the QM box."
          IF (.NOT.force_env%qmmm_env%center_qm_subsys) force_env%qmmm_env%do_translate = .FALSE.
       END IF
       particles_mm => subsys_mm%particles%els
       DO ip=1,SIZE(qm_atom_index)
          particles_qm(ip)%r=particles_mm(qm_atom_index(ip))%r
       END DO

       subsys_section => section_vals_get_subs_vals(force_env%force_env_section, &
                                                    "SUBSYS",error=error)
       CALL write_qs_particle_coordinates(particles_qm,subsys_section,"QM/MM first QM, then MM (0 charges)",error)
       ALLOCATE(charges(SIZE(particles_mm)))
       charges = 0.0_dp
       CALL write_fist_particle_coordinates(particles_mm,subsys_section,charges,error)
       DEALLOCATE(charges)

    END SELECT

  END SUBROUTINE apply_qmmm_translate_low

  ! pbc-aware mean QM atom position
  FUNCTION qmmm_pbc_aware_mean(particles_mm, mm_cell, qm_atom_index)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm
    TYPE(cell_type), POINTER                 :: mm_cell
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index
    REAL(dp)                                 :: qmmm_pbc_aware_mean(3)

    COMPLEX(dp)                              :: mean_z(3)
    INTEGER                                  :: ip

      mean_z = 0.0_dp
      DO ip=1,SIZE(qm_atom_index)
        mean_z = mean_z + EXP(CMPLX(0.0_dp,1.0_dp)*2.0*pi*MATMUL(mm_cell%h_inv, particles_mm(qm_atom_index(ip))%r))
      END DO
      mean_z = mean_z / ABS(mean_z)
      qmmm_pbc_aware_mean = MATMUL(mm_cell%hmat, REAL(LOG(mean_z)/(CMPLX(0.0_dp,1.0_dp)*2.0_dp*pi), dp))
  END FUNCTION

  ! minimum image lattice coordinates difference vector
  FUNCTION qmmm_lat_dv(mm_cell, p1, p2)
    TYPE(cell_type), POINTER                 :: mm_cell
    REAL(dp)                                 :: p1(3), p2(3), qmmm_lat_dv(3)

    REAL(dp)                                 :: lat_v1(3), lat_v2(3)

      lat_v1 = MATMUL(mm_cell%h_inv, p1)
      lat_v2 = MATMUL(mm_cell%h_inv, p2)

      qmmm_lat_dv = lat_v2-lat_v1
      qmmm_lat_dv = qmmm_lat_dv - FLOOR(qmmm_lat_dv)
  END FUNCTION qmmm_lat_dv

  ! find closest QM particle, in position/negative direction if dir is 1 or -1,
  ! respectively
  FUNCTION qmmm_find_closest(particles_mm, mm_cell, qm_atom_index, avoid, p, i_dim, dir, closest_dv) RESULT (closest_ip)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm
    TYPE(cell_type), POINTER                 :: mm_cell
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index
    LOGICAL                                  :: avoid(:)
    REAL(dp)                                 :: p(3)
    INTEGER                                  :: i_dim, dir
    REAL(dp), OPTIONAL                       :: closest_dv
    INTEGER                                  :: closest_ip

    INTEGER                                  :: ip, shift
    REAL(dp)                                 :: lat_dv3(3), lat_dv_shifted, &
                                                my_closest_dv

      closest_ip = -1
      my_closest_dv=HUGE(0.0)
      DO ip=1, SIZE(qm_atom_index)
        IF (avoid(ip)) CYCLE
        lat_dv3 = qmmm_lat_dv(mm_cell, p, particles_mm(qm_atom_index(ip))%r)
        DO shift=-1, 1
           lat_dv_shifted = lat_dv3(i_dim)+shift*1.0_dp
           IF (ABS(lat_dv_shifted) < ABS(my_closest_dv) .AND. (dir*lat_dv_shifted >= 0.0)) THEN
              my_closest_dv = lat_dv_shifted
              closest_ip = ip
           ENDIF
        END DO
      END DO

      IF (PRESENT(closest_dv)) THEN
        closest_dv = my_closest_dv
      ENDIF

  END FUNCTION qmmm_find_closest


! *****************************************************************************
!> \brief Computes a spherical cutoff factor for the QMMM interactions
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2008 created
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE spherical_cutoff_factor(spherical_cutoff, rij, factor, error)
    REAL(KIND=dp), DIMENSION(2), INTENT(IN)  :: spherical_cutoff
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    REAL(KIND=dp), INTENT(OUT)               :: factor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spherical_cutoff_factor', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: r, r0

    r   = SQRT(DOT_PRODUCT(rij,rij))
    r0  = spherical_cutoff(1)-20.0_dp*spherical_cutoff(2)
    factor = 0.5_dp*(1.0_dp-TANH((r-r0)/spherical_cutoff(2)))

  END SUBROUTINE spherical_cutoff_factor

  RECURSIVE FUNCTION qmmm_force_mixing_active(force_env, error) RESULT(active)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: active

    INTEGER                                  :: iforce_eval, nforce_eval
    TYPE(section_vals_type), POINTER         :: qmmm_force_mixing_section

    active = .FALSE. 
    qmmm_force_mixing_section => section_vals_get_subs_vals(force_env%force_env_section,"QMMM%FORCE_MIXING",&
                                 can_return_null=.TRUE.,error=error)
    IF (ASSOCIATED(qmmm_force_mixing_section)) CALL section_vals_get(qmmm_force_mixing_section,&
                                                    explicit=active,error=error)

    IF (active) RETURN

    ! top level wasn't a QMMM with FORCE_MIXING, so loop over sub force envs, e.g. for mixed force env
    nforce_eval = SIZE(force_env%sub_force_env)
    DO iforce_eval=1, nforce_eval
      ! if the current force env is QMMM, then any sub force envs are dups of it generated by FORCE_MIXING, so explicitly avoid recursing into those
      IF (force_env%in_use /= use_qmmm) THEN 
         active = qmmm_force_mixing_active(force_env%sub_force_env(iforce_eval)%force_env, error)
         IF (active) RETURN
      ENDIF
    END DO

  END FUNCTION qmmm_force_mixing_active

END MODULE qmmm_util
