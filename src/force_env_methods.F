!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2003  CP2K developers group                            !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/force_env_types [1.0] *
!!
!!   NAME
!!     force_env_types
!!
!!   FUNCTION
!!     Interface for the force calculations
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     cjm, FEB-20-2001: pass variable box_ref
!!     cjm, SEPT-12-2002: major reorganization
!!     fawzi, APR-12-2003: introduced force_env
!!
!!   SOURCE
!******************************************************************************

MODULE force_env_methods
!  USE fist_environment_types,          ONLY: fist_environment_type,&
!                                             replica_environment_type
!  USE fist_force,                      ONLY: fist_force_control
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_retain
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_retain
  USE kinds,                           ONLY: dbl,&
                                             dp
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
                                             use_qmmm
  USE global_types,                    ONLY: global_environment_type
  USE kg_energy,                       ONLY: kg_energies
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_env_retain,&
                                             kg_environment_type
  USE kg_force,                        ONLY: kg_force_control
  USE metadynamics,                    ONLY: metadyn_add_forces
  USE metadynamics_types,              ONLY: meta_env_retain,&
                                             meta_env_type
  USE qs_energy,                       ONLY: qs_energies
  USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_retain,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force,                        ONLY: qs_forces
  USE simulation_cell,                 ONLY: cell_type
  USE qmmm_types, ONLY: qmmm_env_type, qmmm_env_retain
  USE particle_types, only: particle_type
  USE qmmm_electrostatic_methods, ONLY: qmmm_el_coupling

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='force_env_methods'

  PUBLIC :: force_env_create, &
       force_env_calculate_force, force_env_calculate_energy,&
       force_env_set_cell

  INTEGER, SAVE, PRIVATE :: last_force_env_id=0
  
!!***
!******************************************************************************

CONTAINS

!!****** force_env_types/force_env_calculate_force [1.0] *
!!
!!   NAME
!!     force_env_calculate_force
!!
!!   FUNCTION
!!     Interface routine for force calculations
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     04.2003 adapted to force_env [fawzi]
!!
!!   SOURCE
!******************************************************************************

  RECURSIVE SUBROUTINE force_env_calculate_force ( force_env, box_change, error )
  

    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(IN), OPTIONAL            :: box_change
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_calculate_force', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_box_change

! dbg TYPE(fist_environment_type), POINTER :: fist_env
!------------------------------------------------------------------------------

  failure=.FALSE.
  my_box_change=.FALSE.
  IF (PRESENT(box_change)) my_box_change=box_change
!  NULLIFY(fist_env)

  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)

  SELECT CASE ( force_env%in_use )
   !dbg CASE ( use_fist_force )
   !dbg   fist_env => force_env % fist_env
   !   SELECT CASE ( force_env % fist_env % pimd )
   !   CASE ( .TRUE. )
   !     beads = SIZE ( fist_env % rep_env )
   !     e_pot = 0._dbl
   !     DO i = 1, beads
   !       CALL fist_force_control ( fist_env % rep_env ( i ), &
   !                                 fist_env % ewald_param, &
   !                                 fist_env % potparm, fist_env % thermo, &
   !                                 my_box_change, force_env%globenv )
   !       e_pot = e_pot + fist_env % thermo % pot
   !     ENDDO
   !     fist_env % thermo % pot = e_pot
   !   CASE ( .FALSE. )
!
!          CALL fist_force_control ( fist_env % rep_env ( 1 ), &
!                                    fist_env % ewald_param, &
!                                    fist_env % potparm, fist_env % thermo, &
!                                    my_box_change, force_env%globenv )
!      END SELECT 
    CASE ( use_qs_force )
      CALL qs_forces(force_env%qs_env,force_env%globenv,&
           neighbor_lists_from_scratch=.TRUE.)
    CASE ( use_pol_force )
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    CASE ( use_kg_force )
      CALL kg_force_control ( force_env%kg_env, force_env%globenv )
   CASE (use_qmmm)
      CALL qmmm_forces(force_env,box_change,error)
   CASE default
      CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
  END SELECT

  IF (ASSOCIATED(force_env%meta_env)) THEN
     CALL metadyn_add_forces(force_env%meta_env,force_env, error=error)
  END IF

END SUBROUTINE force_env_calculate_force
!******************************************************************************

!!****f* force_env_types/force_env_calculate_energy [1.0] *
!!
!!   NAME
!!     force_env_calculate_energy
!!
!!   FUNCTION
!!     Is guaranteed to updates only the energy of the actual configuration
!!     (but it *might* also update the force)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - force_env: the en
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!*** **********************************************************************
RECURSIVE SUBROUTINE force_env_calculate_energy(force_env, box_change, error)
    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, OPTIONAL                        :: box_change
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_calculate_energy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_box_change

!dbg  TYPE(fist_environment_type), POINTER :: fist_env
!------------------------------------------------------------------------------

  failure=.FALSE.
  my_box_change=.FALSE.
  IF (PRESENT(box_change)) my_box_change=box_change
!dbg  NULLIFY(fist_env)

  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)

  SELECT CASE ( force_env%in_use )
!    CASE ( use_fist_force ) ! NOT optimized!
!      fist_env => force_env%fist_env
!      SELECT CASE ( fist_env % pimd )
!      CASE ( .TRUE. )
!        beads = SIZE ( fist_env % rep_env )
!        e_pot = 0._dbl
!        DO i = 1, beads
!          CALL fist_force_control ( fist_env % rep_env ( i ), &
!                                    fist_env % ewald_param, &
!                                    fist_env % potparm, fist_env % thermo, &
!                                    my_box_change, force_env%globenv )
!          e_pot = e_pot + fist_env % thermo % pot
!        ENDDO
!        fist_env % thermo % pot = e_pot
!      CASE ( .FALSE. )
!
!          CALL fist_force_control ( fist_env % rep_env ( 1 ), &
!                                    fist_env % ewald_param, &
!                                    fist_env % potparm, fist_env % thermo, &
!                                    my_box_change, force_env%globenv )
!      END SELECT 
    CASE ( use_qs_force )
      CALL qs_energies(force_env%qs_env,force_env%globenv,&
           neighbor_lists_from_scratch=.TRUE.)
    CASE ( use_pol_force )
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    CASE ( use_kg_force ) 
      CALL kg_energies ( force_env%kg_env, force_env%globenv )
   CASE (use_qmmm)
      call qmmm_energy(force_env, box_change, error=error)
  END SELECT
  IF (ASSOCIATED(force_env%meta_env)) THEN
     CALL metadyn_add_forces(force_env%meta_env,force_env, error=error)
  END IF
  
END SUBROUTINE force_env_calculate_energy
!***************************************************************************


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
!SUBROUTINE force_env_create(force_env,fist_env,qs_env,kg_env,para_env,&
SUBROUTINE force_env_create(force_env,qs_env,kg_env,para_env,&
     globenv,meta_env,sub_force_env,qmmm_env,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(global_environment_type), &
      INTENT(inout), OPTIONAL, TARGET        :: globenv
    TYPE(meta_env_type), OPTIONAL, POINTER   :: meta_env
    TYPE(force_env_p_type), DIMENSION(:), POINTER, OPTIONAL :: sub_force_env
    TYPE(qmmm_env_type), POINTER, OPTIONAL   :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat,isubsys
    LOGICAL                                  :: failure
    TYPE(cp_subsystem_p_type), DIMENSION(:), pointer :: subsys

!  TYPE(fist_environment_type), POINTER, OPTIONAL :: fist_env
!  TYPE(replica_environment_type), POINTER :: rep_att
!  NULLIFY(rep_att,atomic_kinds, particles, local_molecules)
!  NULLIFY(atomic_kinds, particles)

  failure=.FALSE.

  ALLOCATE(force_env, stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  IF (.NOT. failure) THEN
!     NULLIFY(force_env%subsys, force_env%fist_env, &
!          force_env%kg_env, force_env%qs_env, force_env%para_env,&
!          force_env%globenv)
     NULLIFY(force_env%subsys, force_env%dyn_coeff_set, &
          force_env%kg_env, force_env%qs_env, force_env%para_env,&
          force_env%globenv, force_env%meta_env, force_env%sub_force_env,&
          force_env%qmmm_env)
     last_force_env_id=last_force_env_id+1
     force_env%id_nr=last_force_env_id
     force_env%ref_count=1
     force_env%in_use=0

     IF (PRESENT(globenv)) THEN
        force_env%globenv => globenv
        force_env%para_env => globenv%para_env
     ELSE
        CPAssert(PRESENT(para_env),cp_failure_level,routineP,error,failure)
     END IF
     IF (PRESENT(para_env)) force_env%para_env => para_env
     CALL cp_para_env_retain(force_env%para_env, error=error)

!     IF (PRESENT(fist_env)) THEN
!        IF (ASSOCIATED(fist_env)) THEN
!           force_env%in_use=use_fist_force
!           force_env%fist_env => fist_env
!        END IF
!     END IF
     IF (PRESENT(kg_env)) THEN
        IF (ASSOCIATED(kg_env)) THEN
           CPPrecondition(force_env%in_use==0,cp_failure_level,routineP,error,failure)
           force_env%in_use=use_kg_force
           force_env%kg_env => kg_env
           CALL kg_env_retain(kg_env,error=error)
        END IF
     END IF
     IF (PRESENT(qs_env)) THEN
        CPPrecondition(PRESENT(globenv),cp_failure_level,routineP,error,failure)

        IF (ASSOCIATED(qs_env)) THEN
           CPPrecondition(force_env%in_use==0,cp_failure_level,routineP,error,failure)
           force_env%in_use=use_qs_force
           force_env%qs_env => qs_env
           CALL qs_env_retain(qs_env,error=error)
        END IF
     END IF
     IF (PRESENT(qmmm_env)) THEN
        CPPrecondition(present(sub_force_env),cp_failure_level,routineP,error,failure)
        force_env%in_use=use_qmmm
        force_env%qmmm_env => qmmm_env
        CALL qmmm_env_retain(qmmm_env,error=error)
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
!     CASE(use_fist_force)
!        CPPrecondition(PRESENT(globenv),cp_failure_level,routineP,error,failure)
!
!        n_rep=SIZE(force_env%fist_env%rep_env)
!        ALLOCATE(force_env%subsys(n_rep),stat=stat)
!        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!
!        DO i=1,n_rep
!           NULLIFY(force_env%subsys(i)%subsys)
!           CALL cp_subsys_create(force_env%subsys(i)%subsys,&
!                para_env=force_env % para_env, error=error)
!
!           rep_att => force_env%fist_env%rep_env(i)
!
!           CALL atomic_kind_list_create(atomic_kinds,&
!                els_ptr=rep_att%atomic_kind_set, owns_els=.FALSE., error=error)
!           CALL particle_list_create(particles,&
!                els_ptr=rep_att%part, owns_els=.FALSE., error=error)
!           CALL mol_struct_list_create(local_molecules,&
!                els_ptr=rep_att%molecule, error=error)
!           CALL cp_subsys_set(force_env%subsys(i)%subsys,&
!                local_particles=rep_att%local_particles,&
!                atomic_kinds=atomic_kinds, particles=particles,&
!                pnode=rep_att%pnode,local_molecules=local_molecules,&
!                error=error)
!           CALL atomic_kind_list_release(atomic_kinds, error=error)
!           CALL particle_list_release(particles, error=error)
!           CALL mol_struct_list_release(local_molecules, error=error)
!        END DO

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
SUBROUTINE force_env_set_cell(force_env, cell, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_set_cell', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
!     CASE (use_fist_force)
!        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     CASE (use_kg_force)
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"unknown in_use flag value "//&
             CPSourceFileRef,&
             error,failure)
     END SELECT
  END IF
END SUBROUTINE force_env_set_cell
!***************************************************************************

!!****f* force_env_methods/qmmm_forces [1.0] *
!!
!!   NAME
!!     qmmm_forces
!!
!!   FUNCTION
!!     calculates the qm/mm forces
!!
!!   NOTES
!!     -
!!
!!   INPUTS
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
RECURSIVE SUBROUTINE qmmm_forces(force_env,box_change,error)
  TYPE(force_env_type), POINTER            :: force_env
  LOGICAL, INTENT(IN), OPTIONAL            :: box_change
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='qmmm_forces',&
        routineP=moduleN//':'//routineN
  TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys_qm,subsys_mm
  INTEGER, DIMENSION(:), POINTER :: qm_atom_index, qm_subsys_index
  INTEGER :: ip, isubsys
  TYPE(particle_type), DIMENSION(:), POINTER :: particles_qm, particles_mm
  REAL(dbl), DIMENSION(3) :: min_coord,max_coord,transl_v
  TYPE(cell_type), POINTER :: cell_loc

  failure=.false.
  NULLIFY(subsys_mm, subsys_qm, qm_atom_index,particles_mm,particles_qm)

  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(force_env%qmmm_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%qmmm_env%ref_count>0,cp_failure_level,routineP,error,failure)
  ! transfer position to qm part
  CALL force_env_get(force_env%sub_force_env(1)%force_env,&
       subsys=subsys_mm,error=error)
  CALL force_env_get(force_env%sub_force_env(2)%force_env,&
       cell=cell_loc,subsys=subsys_qm,error=error)
  qm_atom_index => force_env%qmmm_env%qm_atom_index
  qm_subsys_index => force_env%qmmm_env%qm_subsys_index
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
  WRITE(*,*)max_coord, min_coord
  transl_v = (max_coord + min_coord) / 2.0_dp
  transl_v(1) = transl_v(1) - cell_loc%hmat(1,1)/2.0_dp
  transl_v(2) = transl_v(2) - cell_loc%hmat(2,2)/2.0_dp
  transl_v(3) = transl_v(3) - cell_loc%hmat(3,3)/2.0_dp
  
  DO isubsys=1,SIZE(subsys_mm)     
     particles_mm => subsys_mm(isubsys)%subsys%particles%els
     DO ip=1,subsys_mm(isubsys)%subsys%particles%n_els
        particles_mm(ip)%r = &
             particles_mm(ip)%r - transl_v
     END DO
  END DO
  WRITE(*,'(A)')'PIPPO'
  CALL force_env_write_variables(force_env, 6, .true., error)
  isubsys=-1
  DO ip=1,SIZE(qm_atom_index)
     IF (isubsys/=qm_subsys_index(ip)) THEN
        isubsys=qm_subsys_index(ip)
        particles_mm => subsys_mm(isubsys)%subsys%particles%els
     END IF
     particles_qm(ip)%r=particles_mm(qm_atom_index(ip))%r
  END DO

  ! Initialize ks_qmmm_env
  CALL init_ks_qmmm_env(qs_env=force_env%sub_force_env(2)%force_env%qs_env,error=error)

  ! Compute the QM/MM Electrostatic Potential
  CALL qmmm_el_coupling(force_env%sub_force_env(2)%force_env%qs_env,&
                        force_env%qmmm_env,particles_mm,error)
!TL
! Fist: isubsys = 1
!  isubsys = 1
!  CALL force_env_calculate_force(force_env%sub_force_env(isubsys)%force_env,&
!       box_change=box_change,error=error)
  isubsys = 2
  CALL force_env_calculate_force(force_env%sub_force_env(isubsys)%force_env,&
       box_change=box_change,error=error)
! TL
  isubsys=-1
  DO ip=1,SIZE(qm_atom_index)
     IF (isubsys/=qm_subsys_index(ip)) THEN
        isubsys=qm_subsys_index(ip)
        particles_mm => subsys_mm(isubsys)%subsys%particles%els
     END IF
     particles_mm(qm_atom_index(ip))%f=particles_qm(ip)%f
  END DO  
END SUBROUTINE qmmm_forces
!***************************************************************************


!!****f* force_env_methods/qmmm_energy [1.0] *
!!
!!   NAME
!!     qmmm_energy
!!
!!   FUNCTION
!!     calculates the qm/mm energy
!!
!!   NOTES
!!     -
!!
!!   INPUTS
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
RECURSIVE SUBROUTINE qmmm_energy(force_env,box_change,error)
  TYPE(force_env_type), POINTER            :: force_env
  LOGICAL, INTENT(IN), OPTIONAL            :: box_change
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='qmmm_energy',&
        routineP=moduleN//':'//routineN
  TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys_qm,subsys_mm
  INTEGER, DIMENSION(:), POINTER :: qm_atom_index, qm_subsys_index
  TYPE(particle_type), DIMENSION(:), POINTER :: particles_qm, particles_mm
  INTEGER :: ip,isubsys
  REAL(dbl), DIMENSION(3) :: min_coord,max_coord,transl_v
  TYPE(cell_type), POINTER :: cell_loc

  failure=.false.
  NULLIFY(subsys_mm, subsys_qm, qm_atom_index,particles_mm,particles_qm)

  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(force_env%qmmm_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%qmmm_env%ref_count>0,cp_failure_level,routineP,error,failure)
  ! transfer position to qm part
  CALL force_env_get(force_env%sub_force_env(1)%force_env,&
       subsys=subsys_mm,error=error)
  CALL force_env_get(force_env%sub_force_env(2)%force_env,&
       cell=cell_loc,subsys=subsys_qm,error=error)
  qm_atom_index => force_env%qmmm_env%qm_atom_index
  qm_subsys_index => force_env%qmmm_env%qm_subsys_index
  CPPrecondition(ASSOCIATED(qm_atom_index),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(qm_subsys_index),cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(subsys_qm)==1,cp_failure_level,routineP,error,failure)
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
  WRITE(*,'(A)')'PIPPO1'
  CALL force_env_write_variables(force_env, 6, .true., error)
  WRITE(*,*)max_coord, min_coord
  transl_v = (max_coord + min_coord) / 2.0_dp
  transl_v(1) = transl_v(1) - cell_loc%hmat(1,1)/2.0_dp
  transl_v(2) = transl_v(2) - cell_loc%hmat(2,2)/2.0_dp
  transl_v(3) = transl_v(3) - cell_loc%hmat(3,3)/2.0_dp
  
  DO isubsys=1,SIZE(subsys_mm)     
     particles_mm => subsys_mm(isubsys)%subsys%particles%els
     DO ip=1,subsys_mm(isubsys)%subsys%particles%n_els
        particles_mm(ip)%r = &
             particles_mm(ip)%r - transl_v
     END DO
  END DO
  WRITE(*,'(A)')'PIPPO2'
  CALL force_env_write_variables(force_env, 6, .true., error)
  isubsys=-1
  DO ip=1,SIZE(qm_atom_index)
     IF (isubsys/=qm_subsys_index(ip)) THEN
        isubsys=qm_subsys_index(ip)
        particles_mm => subsys_mm(isubsys)%subsys%particles%els
     END IF
     particles_qm(ip)%r=particles_mm(qm_atom_index(ip))%r
  END DO

  ! Initialize ks_qmmm_env
  CALL init_ks_qmmm_env(qs_env=force_env%sub_force_env(2)%force_env%qs_env,error=error)

  ! Compute the QM/MM Electrostatic Potential
  CALL qmmm_el_coupling(force_env%sub_force_env(2)%force_env%qs_env,&
                        force_env%qmmm_env,particles_mm,error)
!TL
! Fist: isubsys = 1
!  isubsys = 1
!  CALL force_env_calculate_energy(force_env%sub_force_env(isubsys)%force_env,&
!       box_change=box_change,error=error) 
  isubsys = 2
  CALL force_env_calculate_energy(force_env%sub_force_env(isubsys)%force_env,&
       box_change=box_change,error=error) 
!TL
END SUBROUTINE qmmm_energy
!***************************************************************************


!!****f* force_env_methods/init_ks_qmmm_env [1.0] *
!!
!!   NAME
!!     init_ks_qmmm_env
!!
!!   FUNCTION
!!     Initialize the ks_qmmm_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
SUBROUTINE init_ks_qmmm_env(qs_env,error)
  USE qs_ks_methods, ONLY:     qs_ks_qmmm_create
  USE qs_ks_qmmm_types, ONLY:  qs_ks_qmmm_release,&
                               qs_ks_qmmm_env_type
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error
  TYPE(qs_environment_type), OPTIONAL, &
       POINTER                                :: qs_env 
  TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env

  NULLIFY(ks_qmmm_env)
  CALL get_qs_env(qs_env=qs_env,&
       ks_qmmm_env=ks_qmmm_env)
 
!   *** allocate the ks_qmmm env if not allocated yet!**
    IF (.not.ASSOCIATED(ks_qmmm_env)) THEN
       CALL qs_ks_qmmm_create(ks_qmmm_env=ks_qmmm_env,qs_env=qs_env,error=error)
       CALL set_qs_env(qs_env=qs_env, ks_qmmm_env=ks_qmmm_env)
       CALL qs_ks_qmmm_release(ks_qmmm_env=ks_qmmm_env,error=error)
    END IF
END SUBROUTINE init_ks_qmmm_env
!***************************************************************************

END MODULE force_env_methods
