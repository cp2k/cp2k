!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Update a QM/MM calculations with force mixing
!> \par History
!>      5.2004 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qmmmx_update
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_restart_force_eval,        ONLY: update_force_eval
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_release,&
                                             section_vals_type
  USE qmmm_create,                     ONLY: qmmm_env_create
  USE qmmm_types,                      ONLY: qmmm_env_get
  USE qmmmx_types,                     ONLY: qmmmx_env_release,&
                                             qmmmx_env_type
  USE qmmmx_util,                      ONLY: setup_force_mixing_qmmm_sections,&
                                             update_force_mixing_labels
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmmx_update'

  PUBLIC :: qmmmx_update_force_env

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param force_env ...
!> \param root_section ...
! *****************************************************************************
  SUBROUTINE qmmmx_update_force_env(force_env, root_section)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_update_force_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: force_mixing_active, &
                                                labels_changed
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds, new_atomic_kinds
    TYPE(cp_subsys_type), POINTER            :: subsys, subsys_new
    TYPE(distribution_1d_type), POINTER      :: local_particles, &
                                                new_local_particles
    TYPE(qmmmx_env_type), POINTER            :: new_qmmmx_env
    TYPE(section_vals_type), POINTER :: qmmm_core_section, &
      qmmm_extended_Section, qmmm_force_mixing, qmmm_section, subsys_section

! check everything for not null, because sometimes (e.g. metadynamics in parallel) it happens

    IF (.NOT. ASSOCIATED(force_env)) RETURN
    IF (.NOT. ASSOCIATED(force_env%force_env_section)) RETURN
    ! these two should never happen, because the sections exist, but just in case...
    qmmm_section => section_vals_get_subs_vals(force_env%force_env_section,"QMMM",can_return_null=.TRUE.)
    IF (.NOT. ASSOCIATED(qmmm_section)) RETURN
    qmmm_force_mixing => section_vals_get_subs_vals(qmmm_section,"FORCE_MIXING",can_return_null=.TRUE.)
    IF (.NOT. ASSOCIATED(qmmm_force_mixing)) RETURN
    CALL section_vals_get(qmmm_force_mixing,explicit=force_mixing_active)
    IF(.NOT.force_mixing_active) RETURN
    IF(.NOT.ASSOCIATED(force_env%qmmmx_env)) STOP "qmmmx_update_force_env: force_env%qmmmx not associated"

    CALL force_env_get(force_env, subsys=subsys)
    CALL update_force_mixing_labels(subsys, qmmm_section, labels_changed=labels_changed)
    IF (.NOT. labels_changed) RETURN
    CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
            "Adaptive force-mixing labels changed, rebuilding QM/MM calculations! "//&
            CPSourceFileRef,only_ionode=.TRUE.)

    CALL update_force_eval(force_env, root_section, .FALSE.)

    ALLOCATE(new_qmmmx_env)

    ! using CUR_INDICES and CUR_LABELS, create appropriate QM_KIND sections for two QM/MM calculations
    CALL setup_force_mixing_qmmm_sections(subsys, qmmm_section, qmmm_core_section, qmmm_extended_section)

    subsys_section => section_vals_get_subs_vals(force_env%force_env_section,"SUBSYS")
    ![ADAPT] no sure about use_motion_section
    CALL qmmm_env_create(new_qmmmx_env%core,&
      force_env%root_section, force_env%para_env, force_env%globenv,&
      force_env%force_env_section, qmmm_core_section, subsys_section, use_motion_section=.TRUE., &
      prev_subsys=subsys, ignore_outside_box=.TRUE.)
    CALL qmmm_env_create(new_qmmmx_env%ext, &
      force_env%root_section, force_env%para_env, force_env%globenv,&
      force_env%force_env_section, qmmm_extended_section, subsys_section, use_motion_section=.TRUE., &
      prev_subsys=subsys, ignore_outside_box=.TRUE.)

    ! [NB] need to copy wiener process data, since it's not recreated when
    ! fist subsys is recreated by qmmm_env_create
    CALL qmmm_env_get(force_env%qmmmx_env%core, subsys=subsys)
    CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, local_particles=local_particles)
    CALL qmmm_env_get(new_qmmmx_env%core, subsys=subsys_new)
    CALL cp_subsys_get(subsys_new, atomic_kinds=new_atomic_kinds, local_particles=new_local_particles)
    IF (ASSOCIATED(local_particles%local_particle_set)) THEN
       CALL copy_wiener_process(atomic_kinds, local_particles, new_atomic_kinds, new_local_particles)
    ENDIF

    CALL qmmm_env_get(force_env%qmmmx_env%ext, subsys=subsys)
    CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, local_particles=local_particles)
    CALL qmmm_env_get(new_qmmmx_env%ext, subsys=subsys_new)
    CALL cp_subsys_get(subsys_new, atomic_kinds=new_atomic_kinds, local_particles=new_local_particles)
    IF (ASSOCIATED(local_particles%local_particle_set)) THEN
       CALL copy_wiener_process(atomic_kinds, local_particles, new_atomic_kinds, new_local_particles)
    ENDIF

    CALL section_vals_release(qmmm_core_section)
    CALL section_vals_release(qmmm_extended_section)

    ! release old qmmmx_env and point to new one
    CALL qmmmx_env_release(force_env%qmmmx_env)
    force_env%qmmmx_env => new_qmmmx_env

  END SUBROUTINE qmmmx_update_force_env


! *****************************************************************************
!> \brief ...
!> \param from_local_particle_kinds ...
!> \param from_local_particles ...
!> \param to_local_particle_kinds ...
!> \param to_local_particles ...
! *****************************************************************************
  SUBROUTINE copy_wiener_process(from_local_particle_kinds, from_local_particles,&
                                 to_local_particle_kinds, to_local_particles)
    TYPE(atomic_kind_list_type), POINTER     :: from_local_particle_kinds
    TYPE(distribution_1d_type), POINTER      :: from_local_particles
    TYPE(atomic_kind_list_type), POINTER     :: to_local_particle_kinds
    TYPE(distribution_1d_type), POINTER      :: to_local_particles

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_wiener_process', &
      routineP = moduleN//':'//routineN

    INTEGER :: from_iparticle_kind, from_iparticle_local(1), &
      from_nparticle_kind, from_nparticle_local, handle, to_iparticle_global, &
      to_iparticle_kind, to_iparticle_local, to_nparticle_kind, &
      to_nparticle_local, tot_from_nparticle_local, tot_to_nparticle_local
    LOGICAL                                  :: failure, found_it

    CALL timeset(routineN,handle)
    CPPrecondition (ASSOCIATED(from_local_particles),cp_failure_level,routineP,failure)
    CPPrecondition (ASSOCIATED(to_local_particles),cp_failure_level,routineP,failure)

    IF (.NOT. ASSOCIATED(from_local_particles%local_particle_set)) RETURN
    CPPrecondition (.NOT. ASSOCIATED(to_local_particles%local_particle_set),cp_failure_level,routineP,failure)

    from_nparticle_kind = from_local_particle_kinds%n_els
    to_nparticle_kind = to_local_particle_kinds%n_els

    ! make sure total number of particles hasn't changed, even if particle kinds have
    tot_from_nparticle_local = 0
    DO from_iparticle_kind=1, from_nparticle_kind
      tot_from_nparticle_local = tot_from_nparticle_local + from_local_particles%n_el(from_iparticle_kind)
    END DO
    tot_to_nparticle_local = 0
    DO to_iparticle_kind=1, to_nparticle_kind
      tot_to_nparticle_local = tot_to_nparticle_local + to_local_particles%n_el(to_iparticle_kind)
    END DO
    CPPrecondition (tot_from_nparticle_local == tot_to_nparticle_local, cp_failure_level,routineP,failure)

    ALLOCATE(to_local_particles%local_particle_set(to_nparticle_kind))
    DO to_iparticle_kind=1, to_nparticle_kind

       to_nparticle_local = to_local_particles%n_el(to_iparticle_kind)
       ALLOCATE(to_local_particles%local_particle_set(to_iparticle_kind)%rng(to_nparticle_local))

       DO to_iparticle_local=1, to_nparticle_local
          to_iparticle_global = to_local_particles%list(to_iparticle_kind)%array(to_iparticle_local)
          ALLOCATE(to_local_particles%local_particle_set(to_iparticle_kind)%rng(to_iparticle_local)%stream)

          found_it = .FALSE.
          ! find the matching kind/index where this particle was before
          DO from_iparticle_kind=1, from_nparticle_kind
             from_nparticle_local = from_local_particles%n_el(from_iparticle_kind)
             IF (MINVAL(ABS( from_local_particles%list(from_iparticle_kind)%array(1:from_nparticle_local) -&
                 to_iparticle_global)) == 0) THEN
               from_iparticle_local = MINLOC(ABS( &
                   from_local_particles%list(from_iparticle_kind)%array(1:from_nparticle_local) - to_iparticle_global))
               to_local_particles%local_particle_set(to_iparticle_kind)%rng(to_iparticle_local)%stream = &
                  from_local_particles%local_particle_set(from_iparticle_kind)%rng(from_iparticle_local(1))%stream
               found_it = .TRUE.
               EXIT
             ENDIF
          END DO
          CPPostcondition (found_it, cp_failure_level,routineP,failure)

       END DO ! to_iparticle_local

    END DO ! to_iparticle_kind
    CALL timestop(handle)

  END SUBROUTINE copy_wiener_process

END MODULE qmmmx_update
