!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2003  CP2K developers group                            !
!-----------------------------------------------------------------------------!
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

MODULE force_env_types
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_env,                     ONLY: cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_release,&
                                             cp_subsystem_p_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_p_type,&
                                             dyn_coeff_set_release
  USE eip_environment_types,           ONLY: eip_env_get,&
                                             eip_env_release,&
                                             eip_environment_type
  USE fist_energy_types,               ONLY: fist_energy_type
  USE fist_environment_types,          ONLY: fist_env_release,&
                                             fist_environment_type,&
                                             get_fist_env
  USE global_types,                    ONLY: global_environment_type,&
                                             globenv_release
  USE input_section_types,             ONLY: section_vals_type
  USE kg_energy_types,                 ONLY: kg_energy_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_env_release,&
                                             kg_environment_type
  USE kinds,                           ONLY: dp
  USE metadynamics_types,              ONLY: meta_env_release,&
                                             meta_env_retain,&
                                             meta_env_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_release,&
                                             qmmm_env_qm_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_release,&
                                             qs_environment_type
  USE tbmd_energy_types,               ONLY: tbmd_energy_type
  USE tbmd_environment_types,          ONLY: tbmd_env_get,&
                                             tbmd_env_release,&
                                             tbmd_environment_type
  USE virial_types,                    ONLY: virial_release,&
                                             virial_retain,&
                                             virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_env_types'

  INTEGER, PARAMETER, PUBLIC :: use_fist_force=501, use_qs_force=502,&
       use_kg_force=503, use_pol_force=504, use_qmmm=505, &
       use_tbmd_force=506, use_eip_force=507
  PUBLIC :: force_env_type, force_env_p_type
  PUBLIC :: force_env_retain, force_env_release, force_env_set,&
       force_env_get
  
!!***
!******************************************************************************

!!****s* force_env_types/force_env_type [1.0] *
!!
!!   NAME
!!     force_env_type
!!
!!   FUNCTION
!!     wrapper to abstract the force evaluation of the various methods
!!
!!   NOTES
!!     as always direct manipulation of these attributes can have very
!!     bad effects. In this case it can be quite bad and the variables 
!!     might not be up to date. You are warned, use only the get method...
!!
!!   ATTRIBUTES
!!     - id_nr: identification number
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - in_use: which method is in use
!!     - fist_env: the fist environment (allocated only if fist is in use)
!!     - kg_env: the kg environment (allocated only if kg is in use)
!!     - qs_env: qs_env (activated only if quickstep is in use)
!!     - tbmd_env: The TBMD environment (of course only allocated when used)
!!     - globenv: the globenv to have the input that generated this force_env
!!     - para_env: the parallel environment that contains all the parallel
!!       environment of the fragments
!!     - meta_env: the metadynamics environment, allocated if there is
!!       metadynamics
!!     - box_changed: if the cell has canged since the last evaluation
!!       of the forces
!!
!!     read-only attributes (get them *only* through force_env_get):
!!     - subsys: the fragments that build up the actual system.
!!       The subsystems (and thus their para_envs) are supposed to be ordered in
!!       some way, so that looping with i=1,.. and doing blocking mpi
!!       communication in the actual subsystem on all processors should not
!!       deadlock.
!!     - cell: the cell of the actual system
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!     07.2003 tried to adapt to multiple mpi groups
!!
!!   SOURCE
  !***************************************************************************
  TYPE force_env_type
!     PRIVATE
     INTEGER :: id_nr,ref_count,in_use
     LOGICAL :: box_changed
     TYPE ( fist_environment_type ), POINTER :: fist_env 
     TYPE ( kg_environment_type ), POINTER :: kg_env
     TYPE ( meta_env_type ), POINTER :: meta_env
     TYPE ( qs_environment_type ), POINTER :: qs_env
     TYPE ( tbmd_environment_type ), POINTER :: tbmd_env
     TYPE ( eip_environment_type ), POINTER  :: eip_env
     TYPE ( cp_subsystem_p_type ), DIMENSION(:), POINTER :: subsys
     TYPE ( dyn_coeff_set_p_type ), DIMENSION(:),POINTER :: dyn_coeff_set
     TYPE ( global_environment_type ), POINTER :: globenv
     TYPE ( cp_para_env_type ), POINTER :: para_env
     TYPE ( force_env_p_type ), DIMENSION(:), POINTER :: sub_force_env
     TYPE ( qmmm_env_qm_type ), POINTER :: qmmm_env
     TYPE ( virial_type  ), POINTER :: virial
  END TYPE force_env_type

!!****s* force_env_types/force_env_p_type [1.0] *
!!
!!   NAME
!!     force_env_p_type
!!
!!   FUNCTION
!!     allows for the creation of an array of force_env
!!
!!   NOTES
!!     added by MJM for MC swap moves
!!
!!   ATTRIBUTES
!!     - force_env: a force environment (see above)
!!
!!   AUTHOR
!!     MJM
!!
!!***************************************************************************

  TYPE force_env_p_type
    TYPE ( force_env_type ), POINTER :: force_env
  END TYPE force_env_p_type
!!***
!****************************************************************************

CONTAINS

!******************************************************************************
!!****f* force_env_types/force_env_retain [1.0] *
!!
!!   NAME
!!     force_env_retain
!!
!!   FUNCTION
!!     retains the given force env
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - force_env: the force environment to retain
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
SUBROUTINE force_env_retain(force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(force_env%ref_count>0,cp_failure_level,routineP,error)
     force_env%ref_count=force_env%ref_count+1
  END IF
END SUBROUTINE force_env_retain
!***************************************************************************

!!****f* force_env_types/force_env_release [1.0] *
!!
!!   NAME
!!     force_env_release
!!
!!   FUNCTION
!!     releases the given force env
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - force_env: the force environment to release
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
RECURSIVE SUBROUTINE force_env_release(force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(force_env)) THEN
     CPPreconditionNoFail(force_env%ref_count>0,cp_failure_level,routineP,error)
     force_env%ref_count=force_env%ref_count-1
     IF (force_env%ref_count==0) THEN
        force_env%ref_count=1
        SELECT CASE ( force_env%in_use )
        CASE ( use_fist_force )
          CALL fist_env_release(force_env%fist_env,error=error)
        CASE ( use_qs_force )
          CALL qs_env_release(force_env%qs_env,error=error)
        CASE ( use_kg_force )
          CALL kg_env_release(force_env%kg_env,error=error)
        CASE ( use_tbmd_force )
          CALL tbmd_env_release(force_env%tbmd_env, error=error)
        CASE ( use_eip_force )
          CALL eip_env_release(force_env%eip_env, error=error)
        END SELECT
        IF (ASSOCIATED(force_env%subsys)) THEN
           DO i=1,SIZE(force_env%subsys)
              CALL cp_subsys_release(force_env%subsys(i)%subsys, error=error)
           END DO
           DEALLOCATE(force_env%subsys, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(force_env%dyn_coeff_set)) THEN
           DO i=1,SIZE(force_env%dyn_coeff_set)
              CALL dyn_coeff_set_release(force_env%dyn_coeff_set(i)%dyn_coeff_set,error=error)
           END DO
           DEALLOCATE(force_env%dyn_coeff_set, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(force_env%sub_force_env)) THEN
           DO i=1,SIZE(force_env%sub_force_env)
              CALL force_env_release(force_env%sub_force_env(i)%force_env,&
                   error=error)
           END DO
           DEALLOCATE(force_env%sub_force_env,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        CALL cp_para_env_release(force_env%para_env,error=error)
        CALL globenv_release(force_env%globenv)
        ! not deallocated
        CPAssert(.NOT.ASSOCIATED(force_env%kg_env),cp_warning_level,routineP,error,failure)
        CPAssert(.NOT.ASSOCIATED(force_env%fist_env),cp_warning_level,routineP,error,failure)
        CALL meta_env_release(force_env%meta_env,error=error)
        CALL qmmm_env_qm_release(force_env%qmmm_env,error=error)
        CALL virial_release(force_env%virial,error=error)
        force_env%ref_count=0
        DEALLOCATE(force_env,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(force_env)
END SUBROUTINE force_env_release
!***************************************************************************

!!****f* force_env_types/force_env_get [1.0] *
!!
!!   NAME
!!     force_env_get
!!
!!   FUNCTION
!!     returns various attributes about the force environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - force_env: the force environment you what informations about
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!     for the orther arguments see the attributes of force_env_type
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
RECURSIVE SUBROUTINE force_env_get(force_env, in_use, &
     fist_env, qs_env, kg_env, meta_env, dyn_coeff_set, subsys, para_env, &
     potential_energy, kinetic_energy, cell, cell_ref, sub_force_env, qmmm_env,&
     tbmd_env, eip_env, virial, box_changed, use_ref_cell, globenv, input, error)
    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, INTENT(out), OPTIONAL           :: in_use
    TYPE(fist_environment_type), OPTIONAL, &
      POINTER                                :: fist_env
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    TYPE(meta_env_type), OPTIONAL, POINTER   :: meta_env
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: dyn_coeff_set
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: subsys
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: potential_energy, &
                                                kinetic_energy
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    TYPE(force_env_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: sub_force_env
    TYPE(qmmm_env_qm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(tbmd_environment_type), OPTIONAL, &
      POINTER                                :: tbmd_env
    TYPE(eip_environment_type), OPTIONAL, &
      POINTER                                :: eip_env
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    LOGICAL, INTENT(out), OPTIONAL           :: box_changed, use_ref_cell
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: eip_kinetic_energy, &
                                                eip_potential_energy, kin, &
                                                penergy_mm, penergy_qm, pot
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(fist_energy_type), POINTER          :: thermo
    TYPE(kg_energy_type), POINTER            :: kg_energy
    TYPE(qs_energy_type), POINTER            :: qs_energy
    TYPE(qs_environment_type), POINTER       :: my_qs_env
    TYPE(tbmd_energy_type), POINTER          :: tbmd_energy

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN

     SELECT CASE(force_env%in_use)
     CASE (use_qs_force)
        CPPrecondition(ASSOCIATED(force_env%qs_env),cp_failure_level,routineP,error,failure)
        CPPreconditionNoFail(.NOT.PRESENT(fist_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(kg_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(tbmd_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(eip_env),cp_warning_level,routineP,error)
        CALL get_qs_env(force_env%qs_env,&
                        energy=qs_energy,&
                        cell=cell,&
                        cell_ref=cell_ref,&
                        use_ref_cell=use_ref_cell,&
                        input=input,&
                        error=error)
        IF (PRESENT(potential_energy)) potential_energy = qs_energy%total
        CPPrecondition(.NOT.PRESENT(kinetic_energy),cp_failure_level,routineP,error,failure)
     CASE (use_fist_force)
        CPPrecondition(ASSOCIATED(force_env%fist_env),cp_failure_level,routineP,error,failure)
        CPPreconditionNoFail(.NOT.PRESENT(qs_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(kg_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(tbmd_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(eip_env),cp_warning_level,routineP,error)
        CPPrecondition(.NOT.PRESENT(input),cp_failure_level,routineP,error,failure)
        CALL get_fist_env(force_env%fist_env,&
                          thermo=thermo,&
                          cell=cell,&
                          cell_ref=cell_ref,&
                          error=error)
        IF (PRESENT(potential_energy)) potential_energy = thermo%pot
        IF (PRESENT(kinetic_energy)) kinetic_energy = thermo%kin
        IF (PRESENT(use_ref_cell)) use_ref_cell = .FALSE.
     CASE (use_tbmd_force)
        CPPrecondition(ASSOCIATED(force_env%tbmd_env), cp_failure_level, routineP, error, failure)
        CPPreconditionNoFail(.NOT. PRESENT(qs_env), cp_warning_level, routineP, error)
        CPPreconditionNoFail(.NOT. PRESENT(kg_env), cp_warning_level, routineP, error)
        CPPreconditionNoFail(.NOT. PRESENT(fist_env), cp_warning_level, routineP, error)
        CPPreconditionNoFail(.NOT.PRESENT(eip_env), cp_warning_level, routineP, error)
        CALL tbmd_env_get(force_env%tbmd_env,&
                          tbmd_energy=tbmd_energy,&
                          cell=cell, &
                          cell_ref=cell_ref,&
                          use_ref_cell=use_ref_cell,&
                          error=error)
        IF (PRESENT(potential_energy)) THEN
          potential_energy = tbmd_energy%E_bandstructure + tbmd_energy%E_repulsion &
                                                         + tbmd_energy%E_correction
        END IF
        IF (PRESENT(kinetic_energy)) kinetic_energy = tbmd_energy%E_kinetic
     CASE (use_eip_force)
        CPPrecondition(ASSOCIATED(force_env%eip_env), cp_failure_level, routineP, error, failure)
        CPPreconditionNoFail(.NOT. PRESENT(qs_env), cp_warning_level, routineP, error)
        CPPreconditionNoFail(.NOT. PRESENT(kg_env), cp_warning_level, routineP, error)
        CPPreconditionNoFail(.NOT. PRESENT(fist_env), cp_warning_level, routineP, error)
        CPPreconditionNoFail(.NOT.PRESENT(tbmd_env), cp_warning_level, routineP, error)
        CALL eip_env_get(force_env%eip_env,&
                         cell=cell,&
                         cell_ref=cell_ref, &
                         use_ref_cell=use_ref_cell,&
                         eip_potential_energy=eip_potential_energy,&
                         eip_kinetic_energy=eip_kinetic_energy,&
                         error=error)
        IF (PRESENT(potential_energy)) THEN
          potential_energy = eip_potential_energy
        END IF
        IF (PRESENT(kinetic_energy)) kinetic_energy = eip_kinetic_energy
        CPPrecondition(.NOT.PRESENT(kinetic_energy),cp_failure_level,routineP,error,failure)
     CASE (use_kg_force)
        CPPrecondition(ASSOCIATED(force_env%kg_env),cp_failure_level,routineP,error,failure)
        ! do not check ?
        CPPreconditionNoFail(.NOT.PRESENT(qs_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(fist_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(tbmd_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(eip_env),cp_warning_level,routineP,error)
        NULLIFY(dft_control)
        CALL get_kg_env(force_env%kg_env,&
                        energy=kg_energy,&
                        cell=cell,&
                        cell_ref=cell_ref,&
                        use_ref_cell=use_ref_cell,&
                        dft_control=dft_control,&
                        error=error)
        IF(dft_control % qs_control % method == "KG_GPW") THEN
           NULLIFY(my_qs_env,qs_energy)
           CALL get_kg_env(force_env%kg_env,&
                          sub_qs_env=my_qs_env)
           CALL get_qs_env(qs_env=my_qs_env,&
                          energy=qs_energy)
           IF (PRESENT(potential_energy)) potential_energy=qs_energy%total
           NULLIFY(my_qs_env,qs_energy)
        ELSE
          IF (PRESENT(potential_energy)) potential_energy=kg_energy%total        
          IF (PRESENT(kinetic_energy)) kinetic_energy=kg_energy%kin
        END IF
     CASE (use_qmmm)
          IF (PRESENT(cell)) THEN
               CALL force_env_get(force_env%sub_force_env(1)%force_env,&
                                  cell=cell,&
                                  error=error)
          ENDIF
          IF (PRESENT(cell_ref)) THEN
               CALL force_env_get(force_env%sub_force_env(1)%force_env,&
                                  cell_ref=cell_ref,&
                                  error=error)
          ENDIF
          IF (PRESENT(use_ref_cell)) use_ref_cell = .FALSE.
          IF (PRESENT(kinetic_energy)) THEN
             CALL force_env_get(force_env%sub_force_env(1)%force_env,&
                                kinetic_energy=kinetic_energy,&
                                error=error)
          END IF
          IF (PRESENT(potential_energy)) THEN
               CALL force_env_get(force_env%sub_force_env(1)%force_env,&
                                  potential_energy=penergy_mm,&
                                  error=error)
               CALL force_env_get(force_env%sub_force_env(2)%force_env,&
                                  potential_energy=penergy_qm,&
                                  error=error)
               potential_energy = penergy_qm+penergy_mm
          ENDIF
     CASE DEFAULT
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"unknown in_use flag value "//&
             CPSourceFileRef,&
             error,failure)
     END SELECT

     IF (PRESENT(in_use)) in_use=force_env%in_use
     IF (PRESENT(fist_env)) THEN
        fist_env => force_env%fist_env
     END IF
     IF (PRESENT(qs_env)) THEN
        qs_env => force_env%qs_env
     END IF
     IF (PRESENT(kg_env)) THEN
        kg_env => force_env%kg_env
     END IF
     IF (PRESENT(tbmd_env)) THEN
        tbmd_env => force_env%tbmd_env
     END IF
     IF (PRESENT(eip_env)) THEN
        eip_env => force_env%eip_env
     END IF
     IF (PRESENT(dyn_coeff_set)) dyn_coeff_set => force_env%dyn_coeff_set
     IF (PRESENT(subsys)) subsys => force_env%subsys
     IF (PRESENT(para_env)) para_env => force_env%para_env
     IF (ASSOCIATED(force_env%meta_env)) THEN
       IF (PRESENT(potential_energy)) THEN
         potential_energy=potential_energy+force_env%meta_env%cv_energy
       END IF
       IF (PRESENT(kinetic_energy)) THEN
         kinetic_energy=kinetic_energy+force_env%meta_env%ekin_s
       END IF
     END IF
     IF (PRESENT(meta_env)) meta_env => force_env%meta_env
     IF (PRESENT(sub_force_env)) sub_force_env => force_env%sub_force_env
     IF (PRESENT(qmmm_env)) qmmm_env => force_env%qmmm_env
     IF ( PRESENT ( virial ) )   virial => force_env % virial
     IF (PRESENT(box_changed)) box_changed = force_env%box_changed
     IF (PRESENT(globenv)) globenv => force_env%globenv
  END IF
END SUBROUTINE force_env_get
!***************************************************************************

!
!!****f* force_env_types/force_env_set [1.0] *
!!
!!   NAME
!!     force_env_set
!!
!!   FUNCTION
!!     changes some attributes of the force_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - force_env: the force environment where the cell should be changed
!!     - meta_env: the new meta environment
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
SUBROUTINE force_env_set(force_env, meta_env, virial, box_changed, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(meta_env_type), OPTIONAL, POINTER   :: meta_env
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    LOGICAL, INTENT(in), OPTIONAL            :: box_changed
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     IF (PRESENT(meta_env)) THEN
        IF (ASSOCIATED(meta_env)) THEN
           CALL meta_env_retain(meta_env,error=error)
        END IF
        CALL meta_env_release(force_env%meta_env,error=error)
        force_env%meta_env => meta_env
     END IF
     IF (PRESENT(virial)) THEN
        IF (ASSOCIATED(virial)) THEN
           CALL virial_retain(virial,error=error)
        END IF
        CALL virial_release(force_env%virial,error=error)
        force_env%virial => virial
     END IF
     IF (PRESENT(box_changed)) force_env%box_changed=box_changed
  END IF
END SUBROUTINE force_env_set
!***************************************************************************


END MODULE force_env_types

!******************************************************************************
