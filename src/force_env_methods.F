!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2003  CP2K developers group                            !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/force_env_types *
!!
!!   NAME
!!     force_env_types
!!
!!   FUNCTION
!!     Interface for the force calculations
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     cjm, FEB-20-2001: pass variable box_ref
!!     cjm, SEPT-12-2002: major reorganization
!!     fawzi, APR-12-2003: introduced force_env (based on the work by CJM&JGH)
!!     fawzi, NOV-3-2004: reorganized interface for f77 interface
!!
!!   SOURCE
!******************************************************************************

MODULE force_env_methods
  USE realspace_grid_types,            ONLY: rs_pw_to_cube
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_output_handling, ONLY: cp_should_output
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level,&
                                             cp_logger_get_default_unit_nr
  USE cp_para_env,                     ONLY: cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_retain,&
                                             cp_subsys_get,&
                                             cp_subsystem_p_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_retain,&
                                             dyn_coeff_set_p_type
  USE fist_environment_types,          ONLY: fist_environment_type,&
                                             fist_env_release,&
                                             fist_env_retain,&
                                             get_fist_env
  USE fist_force,                      ONLY: fist_force_control
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_read_variables,&
                                             force_env_release,&
                                             force_env_retain,&
                                             force_env_set,&
                                             force_env_type,&
                                             force_env_write_variables,&
                                             use_fist_force,&
                                             use_kg_force,&
                                             use_pol_force,&
                                             use_qs_force,&
                                             use_qmmm,&
                                             use_list_force
  USE global_types,                    ONLY: global_environment_type,&
       globenv_retain, globenv_release
  USE kg_energy,                       ONLY: kg_energies
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_env_retain,&
                                             kg_environment_type
  USE kg_force,                        ONLY: kg_force_control
  USE kinds,                           ONLY: dp
  USE metadynamics,                    ONLY: metadyn_add_forces
  USE metadynamics_types,              ONLY: meta_env_retain,&
                                             meta_env_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qmmm_electrostatic_methods,      ONLY: qmmm_el_coupling
  USE qmmm_methods,                    ONLY: qmmm_forces
  USE qmmm_types,                      ONLY: qmmm_env_qm_type,&
                                             qmmm_env_qm_retain,&
                                             qmmm_links_type
  USE qs_energy,                       ONLY: qs_energies
  USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_retain,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force,                        ONLY: qs_forces
  USE simulation_cell,                 ONLY: cell_type
  USE virial_types,                    ONLY: virial_create, &
                                             virial_type,&
                                             virial_retain
  USE qs_ks_qmmm_methods,              ONLY: ks_qmmm_env_rebuild
  USE qmmm_links_methods,              ONLY: qmmm_link_Imomm_coord,&
                                             qmmm_link_Imomm_forces
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_describe_html, section_describe_index_html, section_vals_type,&
       section_vals_get_subs_vals, section_vals_get, section_vals_val_get

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='force_env_methods'

  PUBLIC :: force_env_create, &
       force_env_calc_energy_force, &
       force_env_set_cell

  INTEGER, SAVE, PRIVATE :: last_force_env_id=0
  
!!***
!******************************************************************************

CONTAINS

!!****** force_env_types/force_env_calc_energy_force [1.0] *
!!
!!   NAME
!!     force_env_calc_energy_force
!!
!!   FUNCTION
!!     Interface routine for force and energy calculations
!!
!!   ARGUMENTS
!!     - force_env: the force_env of which you want the energy and forces
!!     - calc_force: if false the forces *might* be left unchanged
!!       or be unvalid, no guarantee on them is done.Defaults to true
!!
!!   AUTHOR
!!     CJM & fawzi
!!
!!   SOURCE
!******************************************************************************

RECURSIVE SUBROUTINE force_env_calc_energy_force ( force_env, &
     calc_force, error)

    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(IN), OPTIONAL            :: calc_force
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_calc_energy_force', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_box_change,&
         need_f
    INTEGER                                  :: i

!------------------------------------------------------------------------------

  failure=.FALSE.
  need_f=.TRUE.
  IF (PRESENT(calc_force)) need_f=calc_force

  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     my_box_change=force_env%box_changed
     
     SELECT CASE ( force_env%in_use )
     CASE ( use_fist_force )
        CALL fist_force_control( force_env%fist_env,force_env%virial, &
                                 force_env%globenv, my_box_change )
    CASE ( use_qs_force )
       IF (.not.need_f) THEN
          CALL qs_energies(force_env%qs_env,force_env%globenv)
       ELSE
          CALL qs_forces(force_env%qs_env,force_env%globenv)
       END IF
    CASE ( use_pol_force )
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    CASE ( use_kg_force )
       IF (.not.need_f) THEN
          CALL kg_energies ( force_env%kg_env, force_env%globenv )
       ELSE
          CALL kg_force_control ( force_env%kg_env, force_env%globenv )
       END IF
    CASE ( use_qmmm )
      CALL qmmm_energy_and_forces(force_env,need_f,error)
    CASE ( use_list_force )
!HAF: we could split this out into a separate subroutine, but for now:
!HAF: atm we assume all sub-force_env's share the same para_env, just like QS
      CPPrecondition(ASSOCIATED(force_env%force_env_list),cp_failure_level, routineP,error,failure)
      IF (.NOT.failure) THEN
         DO i = 1, SIZE (force_env%force_env_list)
            CPPrecondition(ASSOCIATED(force_env%force_env_list(i)%force_env), cp_failure_level,routineP,error,failure)
            CALL force_env_calc_energy_force( &
                 force_env%force_env_list(i)%force_env, need_f,error)
         END DO
      END IF
    CASE default
      CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
  END SELECT

  IF (ASSOCIATED(force_env%meta_env)) THEN
     CALL metadyn_add_forces(force_env%meta_env,force_env, error=error)
  END IF
  force_env%box_changed=.false.
END IF
END SUBROUTINE force_env_calc_energy_force
!******************************************************************************

!!****f* force_env_types/force_env_create [1.0] *
!!
!!   NAME
!!     force_env_create
!!
!!   FUNCTION
!!     creates and initializes a force environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - force_env: the force env to create
!!     - fist_env, kg_env, qs_env: exactly one of these should be
!!       associated, the one that is active
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE force_env_create(force_env,globenv,force_env_list,fist_env,&
     qs_env,kg_env,meta_env,sub_force_env,qmmm_env,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(force_env_p_type), DIMENSION(:), OPTIONAL, &
      POINTER                                :: force_env_list
    TYPE(fist_environment_type), OPTIONAL, &
      POINTER                                :: fist_env
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    TYPE(global_environment_type), &
      POINTER                                :: globenv
    TYPE(meta_env_type), OPTIONAL, POINTER   :: meta_env
!!HAF: ARGHHHHH now we have the following new from the main branch.... :
!!HAF I guess I can delete my version force_env_list now...
!!HAF why do I even bother doing anything here???? *sigh*
    TYPE(force_env_p_type), DIMENSION(:), POINTER, OPTIONAL :: sub_force_env
    TYPE(qmmm_env_qm_type), POINTER, OPTIONAL   :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat,isubsys
    INTEGER                                  :: i,j,k
    LOGICAL                                  :: failure

    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys
    TYPE(dyn_coeff_set_p_type), DIMENSION(:),POINTER :: sub_dyn_coeff


  failure=.FALSE.

  ALLOCATE ( force_env, stat=stat )
  CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
  IF (.NOT. failure) THEN
     NULLIFY ( force_env%subsys, force_env%fist_env, &
               force_env%kg_env, force_env%qs_env,   &
               force_env%para_env, force_env%globenv, &
               force_env%dyn_coeff_set, force_env%force_env_list, &
               force_env%meta_env, force_env%sub_force_env, &
               force_env%qmmm_env )
     last_force_env_id=last_force_env_id+1
     force_env%id_nr=last_force_env_id
     force_env%ref_count=1
     force_env%in_use=0
     force_env%box_changed=.false.

     force_env%globenv => globenv
     CALL globenv_retain(force_env%globenv,error=error)
     force_env%para_env => globenv%para_env
     CALL cp_para_env_retain(force_env%para_env, error=error)
     IF (PRESENT(force_env_list)) THEN
        IF (ASSOCIATED(force_env_list)) THEN
           CPPrecondition(force_env%in_use==0,cp_failure_level,routineP,error,failure)
           force_env%in_use=use_list_force
           force_env%force_env_list => force_env_list
           DO i = 1, SIZE(force_env_list)
              IF (ASSOCIATED(force_env_list(i)%force_env)) THEN
                 CALL force_env_retain (force_env_list(i)%force_env, &
                                        error = error )
              END IF
           END DO
        END IF
     END IF
     IF (PRESENT(fist_env)) THEN
        IF (ASSOCIATED(fist_env)) THEN
           CPPrecondition(force_env%in_use==0,cp_failure_level,routineP,error,failure)
           force_env%in_use=use_fist_force
           force_env%fist_env => fist_env
           CALL fist_env_retain ( fist_env, error = error )
           CALL virial_create ( force_env % virial, virial_in_use = .TRUE. )
        END IF
     END IF
     IF (PRESENT(kg_env)) THEN
        IF (ASSOCIATED(kg_env)) THEN
           CPPrecondition(force_env%in_use==0,cp_failure_level,routineP,error,failure)
           force_env%in_use=use_kg_force
           force_env%kg_env => kg_env
           CALL kg_env_retain(kg_env,error=error)
           CALL virial_create ( force_env%virial, virial_in_use = .FALSE. )
        END IF
     END IF
     IF (PRESENT(qs_env)) THEN
        IF (ASSOCIATED(qs_env)) THEN
           CPPrecondition(force_env%in_use==0,cp_failure_level,routineP,error,failure)
           force_env%in_use=use_qs_force
           force_env%qs_env => qs_env
           CALL qs_env_retain(qs_env,error=error)
           CALL virial_create ( force_env%virial, virial_in_use = .FALSE. )
        END IF
     END IF
     IF (PRESENT(qmmm_env)) THEN
        CPPrecondition(present(sub_force_env),cp_failure_level,routineP,error,failure)
        force_env%in_use=use_qmmm
        force_env%qmmm_env => qmmm_env
        CALL qmmm_env_qm_retain(qmmm_env,error=error)
        force_env%virial => sub_force_env(1)%force_env%virial
        CALL virial_retain(force_env%virial,error=error)
     END IF
     CPPostcondition(force_env%in_use/=0,cp_failure_level,routineP,error,failure)

     IF (PRESENT(sub_force_env)) THEN
        force_env%sub_force_env => sub_force_env
     END IF

     IF (PRESENT(meta_env)) THEN
        force_env%meta_env => meta_env
        CALL meta_env_retain(meta_env,error=error)
     END IF
     SELECT CASE(force_env%in_use)
     CASE(use_list_force)
        j = 0;
! get total number of subsystems
        DO i = 1, SIZE(force_env%force_env_list)
           CPPrecondition(ASSOCIATED(force_env%force_env_list(i)%force_env), cp_failure_level,routineP,error,failure)
           CALL force_env_get ( force_env % force_env_list (i) % force_env, &
                                 subsys = subsys, error = error )
           j = j + SIZE(subsys)
        END DO
        ALLOCATE ( force_env % subsys ( j ), stat = stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
! assign subsystems
        k = 0
        DO i = 1, SIZE(force_env%force_env_list)
           CALL force_env_get ( force_env % force_env_list (i) % force_env, &
                            subsys = subsys, error = error ) 
           DO j = 1, SIZE(subsys)
              force_env % subsys ( k + j ) = subsys( j )
              CALL cp_subsys_retain ( force_env % subsys ( k + j ) % subsys, &
                                       error = error ) 
           END DO
           k = k + SIZE(subsys)
        END DO
! same procedure for dyn_coeff_set:
        j = 0;
        DO i = 1, SIZE(force_env%force_env_list)
           CALL force_env_get ( force_env % force_env_list (i) % force_env, &
                                dyn_coeff_set = sub_dyn_coeff, error = error )
           j = j + SIZE(sub_dyn_coeff)
        END DO
        ALLOCATE ( force_env % dyn_coeff_set ( j ), stat = stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        k = 0
        DO i = 1, SIZE(force_env%force_env_list)
           CALL force_env_get ( force_env % force_env_list (i) % force_env,&
                                dyn_coeff_set = sub_dyn_coeff, error = error )
           DO j = 1, SIZE(sub_dyn_coeff)
              force_env % dyn_coeff_set ( k + j ) = sub_dyn_coeff(j)
              CALL dyn_coeff_set_retain( &
                     force_env%dyn_coeff_set(k+j)%dyn_coeff_set,error=error)
           END DO
           k = k + SIZE(sub_dyn_coeff)
        END DO

     CASE(use_fist_force)
        ALLOCATE ( force_env % subsys ( 1 ), stat = stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        CALL get_fist_env ( force_env % fist_env, &
                            subsys=force_env % subsys ( 1 ) % subsys, &
                            error = error )
        CALL cp_subsys_retain ( force_env % subsys ( 1 ) % subsys, &
                                error = error )

     CASE(use_kg_force)
        ALLOCATE(force_env%subsys(1), stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        NULLIFY(force_env%subsys(1)%subsys)
        CALL get_kg_env(force_env%kg_env,&
             subsys=force_env%subsys(1)%subsys,&
             error=error)
        CALL cp_subsys_retain(force_env%subsys(1)%subsys,&
             error=error)
        IF(kg_env%dft_control%qs_control%polarization) THEN
           ALLOCATE(force_env%dyn_coeff_set(1), stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           NULLIFY(force_env%dyn_coeff_set(1)%dyn_coeff_set)
           CALL get_kg_env(force_env%kg_env,&
                    dyn_coeff_set=force_env%dyn_coeff_set(1)%dyn_coeff_set,&
                    error=error)
           CALL dyn_coeff_set_retain(force_env%dyn_coeff_set(1)%dyn_coeff_set,&
                    error=error)
        END IF

     CASE(use_qs_force)
        ALLOCATE(force_env%subsys(1), stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        NULLIFY(force_env%subsys(1)%subsys)

        CALL get_qs_env(force_env%qs_env,&
             subsys=force_env%subsys(1)%subsys,&
             error=error)
        CALL cp_subsys_retain(force_env%subsys(1)%subsys,&
             error=error)
     CASE(use_qmmm)
        subsys => force_env%sub_force_env(1)%force_env%subsys
        ALLOCATE(force_env%subsys(size(subsys)), stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        
        DO isubsys=1,SIZE(subsys)
           force_env%subsys(isubsys)%subsys => subsys(isubsys)%subsys
           CALL cp_subsys_retain(subsys(isubsys)%subsys,&
                error=error)
        END DO
     CASE default
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     END SELECT
     
  END IF

END SUBROUTINE force_env_create

!***************************************************************************

!!****f* force_env_types/force_env_set_cell [1.0] *
!!
!!   NAME
!!     force_env_set_cell
!!
!!   FUNCTION
!!     changes the cell of the system
!!
!!   NOTES
!!     for future use, each method should make the changes that it needs
!!     to cope with the new cell
!!
!!   INPUTS
!!     - force_env: the force environment where the cell should be changed
!!     - cell: the new cell
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2003 created [fawzi]
!!
!!*** **********************************************************************
RECURSIVE SUBROUTINE force_env_set_cell(force_env, cell, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_set_cell', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: i

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)
  CPPrecondition(cell%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     SELECT CASE(force_env%in_use)
     CASE (use_qs_force)
        CPPrecondition(ASSOCIATED(force_env%qs_env),cp_failure_level,routineP,error,failure)
        CALL set_qs_env(force_env%qs_env,&
             cell=cell,&
             error=error)
        CALL qs_env_rebuild_pw_env(force_env%qs_env,error=error) !avoidable?
     CASE (use_list_force)
        DO i = 1, SIZE(force_env%force_env_list)
           CALL force_env_set_cell(force_env%force_env_list(i)%force_env, &
                cell, error)
        END DO
     CASE (use_fist_force)
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     CASE (use_kg_force)
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"unknown in_use flag value "//&
             CPSourceFileRef,&
             error,failure)
     END SELECT
     force_env%box_changed=.true.
  END IF
END SUBROUTINE force_env_set_cell
!***************************************************************************

!!****f* force_env_methods/qmmm_energy_and_forces [1.0] *
!!
!!   NAME
!!     qmmm_forces
!!
!!   FUNCTION
!!     calculates the qm/mm energy and forces
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - calc_force: if also the forces should be calculated
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [fawzi]
!!
!!*** **********************************************************************
RECURSIVE SUBROUTINE qmmm_energy_and_forces(force_env,calc_force,error)
  TYPE(force_env_type), POINTER            :: force_env
  LOGICAL, INTENT(IN), OPTIONAL            :: calc_force
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure,need_f,qmmm_link, qmmm_link_imomm
  character(len=*), parameter :: routineN='qmmm_energy_and_forces',&
       routineP=moduleN//':'//routineN
  TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys_qm,subsys_mm
  INTEGER, DIMENSION(:), POINTER :: qm_atom_index, qm_subsys_index
  INTEGER :: ip, isubsys
  INTEGER, PARAMETER :: fist_subsys = 1, qs_subsys = 2
  TYPE(particle_type), DIMENSION(:), POINTER :: particles_qm, particles_mm
  REAL(KIND=dp), DIMENSION(3) :: min_coord,max_coord,transl_v
  TYPE(cell_type), POINTER :: qm_cell, mm_cell
  REAL(KIND=dp) :: energy_mm, energy_qm
  TYPE(cp_logger_type), pointer :: logger
  TYPE(section_vals_type), POINTER :: print_key, force_env_section, qmmm_section
  TYPE(qmmm_links_type), POINTER :: qmmm_links
  LOGICAL, SAVE :: do_translate=.TRUE.
  

  min_coord       =  HUGE(0.0_dp)
  max_coord       = -HUGE(0.0_dp)
  failure         = .FALSE.
  need_f          = .TRUE.
  qmmm_link       = .FALSE.
  qmmm_link_imomm = .FALSE.
  logger => cp_error_get_logger(error)
  IF (PRESENT(calc_force)) need_f=calc_force
  NULLIFY(subsys_mm, subsys_qm, qm_atom_index,particles_mm,particles_qm, qm_cell, mm_cell)
  NULLIFY(print_key, force_env_section, qmmm_section)
  force_env_section => &
       section_vals_get_subs_vals(force_env%sub_force_env(qs_subsys)%force_env%globenv%input_file,&
       "FORCE_EVAL",error=error)
  qmmm_section => section_vals_get_subs_vals(force_env_section,"QMMM",error=error)

  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(force_env%qmmm_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%qmmm_env%ref_count>0,cp_failure_level,routineP,error,failure)

  CALL force_env_get(force_env%sub_force_env(fist_subsys)%force_env,&
       cell=mm_cell,subsys=subsys_mm,error=error)
  CALL force_env_get(force_env%sub_force_env(qs_subsys)%force_env,&
       cell=qm_cell,subsys=subsys_qm,error=error)
  qm_atom_index   => force_env%qmmm_env%qm_atom_index
  qm_subsys_index => force_env%qmmm_env%qm_subsys_index
  qmmm_link       =  force_env%qmmm_env%qmmm_link
  qmmm_links      => force_env%qmmm_env%qmmm_links
  IF (qmmm_link) THEN
     CPPrecondition(ASSOCIATED(qmmm_links),cp_failure_level,routineP,error,failure)
     qmmm_link_imomm = (SIZE(qmmm_links%imomm) /= 0)
  END IF
  CPPrecondition(ASSOCIATED(qm_atom_index),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(qm_subsys_index),cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(subsys_qm)==1,cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(subsys_mm)==1,cp_failure_level,routineP,error,failure)

  particles_qm => subsys_qm(1)%subsys%particles%els
  isubsys=-1
  DO ip=1,SIZE(qm_atom_index)
     IF (isubsys/=qm_subsys_index(ip)) THEN
        isubsys=qm_subsys_index(ip)
        particles_mm => subsys_mm(isubsys)%subsys%particles%els
     END IF
     min_coord=MIN(min_coord,particles_mm(qm_atom_index(ip))%r)
     max_coord=MAX(max_coord,particles_mm(qm_atom_index(ip))%r)
  END DO
  IF (do_translate) THEN
     !
     ! The first time we always translate all the system in order 
     ! to centre the QM system in the box. For DEBUG run_type after
     ! the first translation nomore translation have to be done !!
     !
     transl_v = (max_coord + min_coord) / 2.0_dp
     transl_v(1) = transl_v(1) - qm_cell%hmat(1,1)/2.0_dp
     transl_v(2) = transl_v(2) - qm_cell%hmat(2,2)/2.0_dp
     transl_v(3) = transl_v(3) - qm_cell%hmat(3,3)/2.0_dp

     IF (logger%para_env%mepos==logger%para_env%source) &
          WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.FALSE.),fmt='(A)')&
          " Translating the system in order to center the QM fragment in the QM box."

     DO isubsys=1,SIZE(subsys_mm)     
        particles_mm => subsys_mm(isubsys)%subsys%particles%els
        DO ip=1,subsys_mm(isubsys)%subsys%particles%n_els
           particles_mm(ip)%r = &
                particles_mm(ip)%r - transl_v
        END DO
     END DO
     IF (.NOT.force_env%qmmm_env%center_qm_subsys) do_translate = .FALSE.
  END IF
  !  CALL force_env_write_variables(force_env, 6, .true., error)
  isubsys=-1
  DO ip=1,SIZE(qm_atom_index)
     IF (isubsys/=qm_subsys_index(ip)) THEN
        isubsys=qm_subsys_index(ip)
        particles_mm => subsys_mm(isubsys)%subsys%particles%els
     END IF
     particles_qm(ip)%r=particles_mm(qm_atom_index(ip))%r
  END DO
  ! If present QM/MM links (just IMOMM) correct the position of the qm-link atom
  IF (qmmm_link_imomm) CALL qmmm_link_Imomm_coord(qmmm_links, particles_qm, qm_atom_index, error)
  
  ! Initialize ks_qmmm_env
  CALL ks_qmmm_env_rebuild(qs_env=force_env%sub_force_env(qs_subsys)%force_env%qs_env,&
       qmmm_env=force_env%qmmm_env,error=error)

  ! Compute the short range QM/MM Electrostatic Potential
  CALL qmmm_el_coupling(globenv=force_env%sub_force_env(qs_subsys)%force_env%globenv,&
       qs_env=force_env%sub_force_env(qs_subsys)%force_env%qs_env,&
       qmmm_env=force_env%qmmm_env,&
       mm_particles=particles_mm,&
       mm_cell=mm_cell,&
       error=error)

  ! Fist
  CALL force_env_calc_energy_force(force_env%sub_force_env(fist_subsys)%force_env,&
       calc_force=need_f,error=error)
  !
  ! Print Out information on fist energy calculation...
  !
  CALL force_env_get(force_env%sub_force_env(fist_subsys)%force_env,&
                     potential_energy=energy_mm,&
                     error=error)
  print_key => section_vals_get_subs_vals(qmmm_section,"fist_derivatives",&
       error=error)
  IF (cp_should_output(logger%iter_info,print_key,error=error).AND.&
       logger%para_env%mepos==logger%para_env%source) THEN
     WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.FALSE.),&
          fmt='(/1X,A,F15.9)')"Energy after Fist calculation: ",energy_mm
     IF (need_f) THEN
        WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.FALSE.),&
             fmt='(/1X,A)')"Derivatives on all atoms after Fist calculation: "
        DO ip=1,SIZE(particles_mm)
           WRITE(unit=cp_logger_get_default_unit_nr(logger,local=.FALSE.),&
                fmt='(1X,3F15.9)')particles_mm(ip)%f
        END DO
     END IF
  END IF
  !
  ! QS
  CALL force_env_calc_energy_force(force_env%sub_force_env(qs_subsys)%force_env,&
       calc_force=need_f,error=error)
  !
  ! Print Out information on QS energy calculation...
  !
  CALL force_env_get(force_env%sub_force_env(qs_subsys)%force_env,&
                     potential_energy=energy_qm,&
                     error=error)

  print_key => section_vals_get_subs_vals(qmmm_section,"qs_derivatives",&
       error=error)
  IF (cp_should_output(logger%iter_info,print_key,error=error).AND.&
       logger%para_env%mepos==logger%para_env%source) THEN
     WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.FALSE.),&  
          fmt='(/1X,A,F15.9)')"Energy after QS calculation: ",energy_qm
     IF (need_f) THEN
        WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.FALSE.),&  
             fmt='(/1X,A)')"Derivatives on qm atoms after QS calculation: "
        DO ip=1,SIZE(particles_qm)
           WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.FALSE.),&
                fmt='(1X,3F15.9)')particles_qm(ip)%f
        END DO
     END IF
  END IF
  !
  ! QM/MM Interaction Potential forces
  !
  CALL qmmm_forces(force_env%sub_force_env(qs_subsys)%force_env%qs_env,&
       force_env%qmmm_env,particles_mm,&
       mm_cell=mm_cell,&
       calc_force=need_f,error=error)
  !
  ! Print Out information on QS energy calculation...
  !
  CALL force_env_get(force_env%sub_force_env(qs_subsys)%force_env,&
                     potential_energy=energy_qm,&
                     error=error)

  IF (need_f) THEN
     isubsys=-1
     ! If present QM/MM links (just IMOMM) correct the position of the qm-link atom
     IF (qmmm_link_imomm) CALL qmmm_link_Imomm_forces(qmmm_links,particles_qm,qm_atom_index,error)
     DO ip=1,SIZE(qm_atom_index)
        IF (isubsys/=qm_subsys_index(ip)) THEN
           isubsys=qm_subsys_index(ip)
           particles_mm => subsys_mm(isubsys)%subsys%particles%els
        END IF
        particles_mm(qm_atom_index(ip))%f=particles_mm(qm_atom_index(ip))%f+&
             particles_qm(ip)%f
     END DO
  END IF

  print_key => section_vals_get_subs_vals(qmmm_section,"qmmm_derivatives",&
       error=error)
  IF (cp_should_output(logger%iter_info,print_key,error=error).AND.&
       logger%para_env%mepos==logger%para_env%source) THEN
     WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.FALSE.),&  
          fmt='(/1X,A,F15.9)')"Energy after QMMM calculation: ",energy_qm
     IF (need_f) THEN
        WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.FALSE.),&
             fmt='(/1X,A)')"Derivatives on all atoms after QMMM calculation: "
        DO ip=1,SIZE(particles_mm)
           WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.FALSE.),&
                fmt='(1X,3F15.9)')particles_mm(ip)%f
        END DO
     END IF
  END IF

END SUBROUTINE qmmm_energy_and_forces
!***************************************************************************

END MODULE force_env_methods
