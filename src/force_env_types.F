!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface for the force calculations
!> \par History
!>      cjm, FEB-20-2001: pass variable box_ref
!>      cjm, SEPT-12-2002: major reorganization
!>      fawzi, APR-12-2003: introduced force_env
!>      cjm, FEB-27-2006: no more box_change
!> \author CJM & JGH
! *****************************************************************************
MODULE force_env_types
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_env,                     ONLY: cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_types,                 ONLY: cp_result_release,&
                                             cp_result_retain,&
                                             cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_release,&
                                             cp_subsys_type
  USE cpot_types,                      ONLY: cpot_release,&
                                             cpot_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_release,&
                                             dyn_coeff_set_type
  USE eip_environment_types,           ONLY: eip_env_get,&
                                             eip_env_release,&
                                             eip_environment_type
  USE ep_types,                        ONLY: ep_env_release,&
                                             ep_env_type
  USE extended_system_types,           ONLY: extended_energy_type
  USE f77_blas
  USE fist_energy_types,               ONLY: fist_energy_type
  USE fist_environment_types,          ONLY: fist_env_release,&
                                             fist_environment_type,&
                                             get_fist_env
  USE fp_types,                        ONLY: fp_env_release,&
                                             fp_env_retain,&
                                             fp_type
  USE global_types,                    ONLY: global_environment_type,&
                                             globenv_release
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kg_energy_types,                 ONLY: kg_energy_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_env_release,&
                                             kg_environment_type
  USE kinds,                           ONLY: dp
  USE metadynamics_types,              ONLY: meta_env_release,&
                                             meta_env_retain,&
                                             meta_env_type
  USE mixed_energy_types,              ONLY: mixed_energy_type
  USE mixed_environment_types,         ONLY: get_mixed_env,&
                                             mixed_env_release,&
                                             mixed_environment_type
  USE particle_list_types,             ONLY: particle_list_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_release,&
                                             qmmm_env_qm_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_release,&
                                             qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_release,&
                                             virial_retain,&
                                             virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_env_types'

  INTEGER, PARAMETER, PUBLIC :: use_fist_force = 501, &
                                use_qs_force   = 502, &
                                use_kg_force   = 503, &
                                use_qmmm       = 504, &
                                use_eip_force  = 505, &
                                use_ep_force   = 506, &
                                use_mixed_force= 507

  CHARACTER(LEN=10), DIMENSION(501:507), PARAMETER, PUBLIC :: &
       use_prog_name = (/&
       "FIST  ", &
       "QS    ", &
       "KG    ", &
       "QMMM  ", &
       "EIP   ", &
       "EP    ", &
       "MIXED "/)

  PUBLIC :: force_env_type,&
            force_env_p_type
  PUBLIC :: force_env_retain,&
            force_env_release,&
            force_env_set,&
            force_env_get,&
            force_env_pos_get,&
            force_env_get_natom,&
            multiple_fe_list

! *****************************************************************************
!> \brief wrapper to abstract the force evaluation of the various methods
!> \param id_nr identification number
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param in_use which method is in use
!> \param fist_env the fist environment (allocated only if fist is in use)
!> \param kg_env the kg environment (allocated only if kg is in use)
!> \param qs_env qs_env (activated only if quickstep is in use)
!> \param globenv the globenv to have the input that generated this force_env
!> \param para_env the parallel environment that contains all the parallel
!>        environment of the fragments
!> \param meta_env the metadynamics environment, allocated if there is
!>        metadynamics
!> \param fp_env the flexible partitioning environment
!>      read-only attributes (get them *only* through force_env_get):
!> \param subsys the fragments that build up the actual system.
!> \param ep_env environment for an ep calculation
!> \param cell the cell of the actual system
!> \note
!>      as always direct manipulation of these attributes can have very
!>      bad effects. In this case it can be quite bad and the variables
!>      might not be up to date. You are warned, use only the get method...
!> \par History
!>      04.2003 created [fawzi]
!>      07.2003 tried to adapt to multiple mpi groups
!> \author fawzi
! *****************************************************************************
  TYPE force_env_type
!     PRIVATE
     INTEGER :: id_nr,ref_count,in_use,method_name_id
     TYPE ( fist_environment_type ), POINTER             :: fist_env
     TYPE ( kg_environment_type ), POINTER               :: kg_env
     TYPE ( meta_env_type ), POINTER                     :: meta_env
     TYPE ( fp_type ), POINTER                           :: fp_env
     TYPE ( qs_environment_type ), POINTER               :: qs_env
     TYPE ( eip_environment_type ), POINTER              :: eip_env
     TYPE ( cp_subsys_type ), POINTER                    :: subsys
     TYPE ( dyn_coeff_set_type ), POINTER                :: dyn_coeff_set
     TYPE ( global_environment_type ), POINTER           :: globenv
     TYPE ( cp_para_env_type ), POINTER                  :: para_env
     TYPE ( force_env_p_type ), DIMENSION(:), POINTER    :: sub_force_env
     TYPE ( qmmm_env_qm_type ), POINTER                  :: qmmm_env
     TYPE ( virial_type  ), POINTER                      :: virial
     TYPE ( ep_env_type ), POINTER                       :: ep_env
     TYPE ( mixed_environment_type ), POINTER            :: mixed_env     
     TYPE ( section_vals_type), POINTER                  :: force_env_section
     TYPE ( section_vals_type), POINTER                  :: root_section
     TYPE ( extended_energy_type), POINTER               :: extended_energies
     TYPE ( cpot_type), POINTER                          :: cpot_env
     TYPE ( cp_result_type), POINTER                     :: results 
  END TYPE force_env_type

! *****************************************************************************
!> \brief allows for the creation of an array of force_env
!> \param force_env a force environment (see above)
!> \note
!>      added by MJM for MC swap moves
!> \author MJM
! *****************************************************************************
  TYPE force_env_p_type
    TYPE ( force_env_type ), POINTER :: force_env
  END TYPE force_env_p_type

CONTAINS

! *****************************************************************************
!> \brief retains the given force env
!> \param force_env the force environment to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      see doc/ReferenceCounting.html
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE force_env_retain(force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief releases the given force env
!> \param force_env the force environment to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      see doc/ReferenceCounting.html
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  RECURSIVE SUBROUTINE force_env_release(force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, my_group, stat
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: my_error

    failure=.FALSE.
    IF (ASSOCIATED(force_env)) THEN
       CPPreconditionNoFail(force_env%ref_count>0,cp_failure_level,routineP,error)
       force_env%ref_count=force_env%ref_count-1
       IF (force_env%ref_count==0) THEN
          force_env%ref_count=1
          ! Deallocate Subsys
          IF (ASSOCIATED(force_env%subsys)) THEN
             CALL cp_subsys_release(force_env%subsys, error=error)
          END IF
          ! Deallocate KG structures
          IF (ASSOCIATED(force_env%dyn_coeff_set)) THEN
             CALL dyn_coeff_set_release(force_env%dyn_coeff_set,error=error)
          END IF
          ! Deallocate SUB_FORCE_ENV
          IF (ASSOCIATED(force_env%sub_force_env)) THEN
             DO i=1,SIZE(force_env%sub_force_env)
                IF (.NOT.ASSOCIATED(force_env%sub_force_env(i)%force_env)) CYCLE
                ! Use the proper error to deallocate..
                my_error=error
                IF (force_env%in_use==use_mixed_force) THEN
                   my_group = force_env%mixed_env%group_distribution(force_env%para_env%mepos)
                   my_error = force_env%mixed_env%sub_error(my_group+1)
                END IF
                CALL force_env_release(force_env%sub_force_env(i)%force_env,&
                     error=my_error)
                CALL cp_error_check(my_error, failure)
             END DO
             DEALLOCATE(force_env%sub_force_env,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          ! Deallocate Energy type
          IF (ASSOCIATED(force_env%extended_energies)) THEN
             DEALLOCATE(force_env%extended_energies,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF

          SELECT CASE ( force_env%in_use )
          CASE ( use_fist_force )
             CALL fist_env_release(force_env%fist_env,error=error)
          CASE ( use_qs_force )
             CALL qs_env_release(force_env%qs_env,error=error)
          CASE ( use_kg_force )
             CALL kg_env_release(force_env%kg_env,error=error)
          CASE ( use_eip_force )
             CALL eip_env_release(force_env%eip_env, error=error)
          CASE ( use_ep_force )
             CALL ep_env_release(force_env%ep_env, error=error)
          CASE (use_mixed_force)
             CALL mixed_env_release(force_env%mixed_env,error=error)
          END SELECT
          CALL cp_result_release(results=force_env%results,error=error)
          CALL cpot_release(force_env%cpot_env,error=error)
          CALL globenv_release(force_env%globenv,force_env%para_env,error=error)
          CALL cp_para_env_release(force_env%para_env,error=error)
          ! Not deallocated
          CPAssert(.NOT.ASSOCIATED(force_env%kg_env),cp_warning_level,routineP,error,failure)
          CPAssert(.NOT.ASSOCIATED(force_env%fist_env),cp_warning_level,routineP,error,failure)
          CPAssert(.NOT.ASSOCIATED(force_env%qs_env),cp_warning_level,routineP,error,failure)
          CPAssert(.NOT.ASSOCIATED(force_env%eip_env),cp_warning_level,routineP,error,failure)
          CPAssert(.NOT.ASSOCIATED(force_env%ep_env),cp_warning_level,routineP,error,failure)
          CPAssert(.NOT.ASSOCIATED(force_env%mixed_env),cp_warning_level,routineP,error,failure)
          CALL meta_env_release(force_env%meta_env,error=error)
          CALL fp_env_release(force_env%fp_env,error=error)
          CALL qmmm_env_qm_release(force_env%qmmm_env,error=error)
          CALL virial_release(force_env%virial,error=error)
          CALL section_vals_release(force_env%force_env_section,error=error)
          CALL section_vals_release(force_env%root_section,error=error)
          force_env%ref_count=0
          DEALLOCATE(force_env,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(force_env)
  END SUBROUTINE force_env_release

! *****************************************************************************
!> \brief returns various attributes about the force environment
!> \param force_env the force environment you what informations about
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>      for the orther arguments see the attributes of force_env_type
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
RECURSIVE SUBROUTINE force_env_get(force_env, in_use, &
     fist_env, qs_env, kg_env, meta_env, fp_env, dyn_coeff_set, subsys, para_env, &
     potential_energy, kinetic_energy, harmonic_shell, kinetic_shell,&
     cell, cell_ref, sub_force_env, qmmm_env, &
     eip_env, virial, use_ref_cell, globenv, input, ep_env, force_env_section, &
     method_name_id, root_section, mixed_env,results, error)
    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, INTENT(out), OPTIONAL           :: in_use
    TYPE(fist_environment_type), OPTIONAL, &
      POINTER                                :: fist_env
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    TYPE(meta_env_type), OPTIONAL, POINTER   :: meta_env
    TYPE(fp_type), OPTIONAL, POINTER         :: fp_env
    TYPE(dyn_coeff_set_type), OPTIONAL, &
      POINTER                                :: dyn_coeff_set
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: potential_energy, &
                                                kinetic_energy, &
                                                harmonic_shell, kinetic_shell
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    TYPE(force_env_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: sub_force_env
    TYPE(qmmm_env_qm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(eip_environment_type), OPTIONAL, &
      POINTER                                :: eip_env
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    LOGICAL, INTENT(out), OPTIONAL           :: use_ref_cell
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
    TYPE(ep_env_type), OPTIONAL, POINTER     :: ep_env
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: force_env_section
    INTEGER, INTENT(out), OPTIONAL           :: method_name_id
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: root_section
    TYPE(mixed_environment_type), OPTIONAL, &
      POINTER                                :: mixed_env
    TYPE(cp_result_type), OPTIONAL, POINTER  :: results
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: eip_kinetic_energy, &
                                                eip_potential_energy, &
                                                penergy_mm, penergy_qm
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(fist_energy_type), POINTER          :: thermo
    TYPE(kg_energy_type), POINTER            :: kg_energy
    TYPE(mixed_energy_type), POINTER         :: mixed_energy
    TYPE(qs_energy_type), POINTER            :: qs_energy
    TYPE(qs_environment_type), POINTER       :: my_qs_env

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN

     SELECT CASE(force_env%in_use)
     CASE(use_ep_force)
        CPPrecondition(ASSOCIATED(force_env%ep_env),cp_failure_level,routineP,error,failure)
        CALL get_qs_env(force_env%ep_env%main_qs_env,&
                        cell=cell,&
                        cell_ref=cell_ref,&
                        use_ref_cell=use_ref_cell,&
                        input=input,&
                        error=error)
        IF (PRESENT(potential_energy)) potential_energy = force_env%ep_env%energy%e_tot
        CPPrecondition(.NOT.PRESENT(kinetic_energy),cp_failure_level,routineP,error,failure)
     CASE (use_qs_force)
        CPPrecondition(ASSOCIATED(force_env%qs_env),cp_failure_level,routineP,error,failure)
        CPPreconditionNoFail(.NOT.PRESENT(fist_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(kg_env),cp_warning_level,routineP,error)
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
        IF (PRESENT(kinetic_shell)) kinetic_shell = thermo%kin_shell
        IF (PRESENT(harmonic_shell)) harmonic_shell = thermo%harm_shell
     CASE (use_eip_force)
        CPPrecondition(ASSOCIATED(force_env%eip_env), cp_failure_level, routineP, error, failure)
        CPPreconditionNoFail(.NOT. PRESENT(qs_env), cp_warning_level, routineP, error)
        CPPreconditionNoFail(.NOT. PRESENT(kg_env), cp_warning_level, routineP, error)
        CPPreconditionNoFail(.NOT. PRESENT(fist_env), cp_warning_level, routineP, error)
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
                          sub_qs_env=my_qs_env,error=error)
           CALL get_qs_env(qs_env=my_qs_env,&
                          energy=qs_energy,error=error)
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
     CASE (use_mixed_force)
        CPPrecondition(ASSOCIATED(force_env%mixed_env),cp_failure_level,routineP,error,failure)
        CPPreconditionNoFail(.NOT.PRESENT(qs_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(kg_env),cp_warning_level,routineP,error)
        CPPreconditionNoFail(.NOT.PRESENT(eip_env),cp_warning_level,routineP,error)
        CPPrecondition(.NOT.PRESENT(input),cp_failure_level,routineP,error,failure)
        CALL get_mixed_env(force_env%mixed_env,&
                          mixed_energy=mixed_energy,&
                          cell=cell,&
                          cell_ref=cell_ref,&
                          error=error)
        IF (PRESENT(potential_energy)) potential_energy = mixed_energy%pot
        IF (PRESENT(kinetic_energy)) kinetic_energy = mixed_energy%kin
        IF (PRESENT(use_ref_cell)) use_ref_cell = .FALSE.
     CASE DEFAULT
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"unknown in_use flag value "//&
             CPSourceFileRef,&
             error,failure)
     END SELECT

     IF (PRESENT(force_env_section)) force_env_section => force_env%force_env_section
     IF (PRESENT(in_use)) in_use=force_env%in_use
     IF (PRESENT(method_name_id))method_name_id=force_env%method_name_id
     IF (PRESENT(fist_env)) THEN
        fist_env => force_env%fist_env
     END IF
     IF (PRESENT(qs_env)) THEN
        qs_env => force_env%qs_env
     END IF
     IF (PRESENT(kg_env)) THEN
        kg_env => force_env%kg_env
     END IF
     IF (PRESENT(eip_env)) THEN
        eip_env => force_env%eip_env
     END IF
     IF (PRESENT(results)) results => force_env%results     
     IF (PRESENT(dyn_coeff_set)) dyn_coeff_set => force_env%dyn_coeff_set
     IF (PRESENT(subsys)) subsys => force_env%subsys
     IF (PRESENT(para_env)) para_env => force_env%para_env
     ! adjust the total energy for the metadynamics
     IF (ASSOCIATED(force_env%meta_env)) THEN
       IF (PRESENT(potential_energy)) THEN
         potential_energy=potential_energy+force_env%meta_env%cv_energy
       END IF
       IF (PRESENT(kinetic_energy)) THEN
         kinetic_energy=kinetic_energy+force_env%meta_env%ekin_s
       END IF
     END IF
     ! adjust the total energy for the flexible partitioning
     IF (ASSOCIATED(force_env%fp_env) .AND. PRESENT(potential_energy)) THEN
        IF (force_env%fp_env%use_fp) THEN
             potential_energy=potential_energy+force_env%fp_env%energy
        ENDIF
     ENDIF
     ! extended energies
     IF (ASSOCIATED(force_env%extended_energies) .AND. PRESENT(potential_energy)) THEN
        potential_energy=potential_energy+ force_env%extended_energies%energy_g3x3_restraint+&
                                           force_env%extended_energies%energy_fixd_restraint+&
                                           force_env%extended_energies%energy_g4x6_restraint+&
                                           force_env%extended_energies%energy_colv_restraint
     END IF
     IF (PRESENT(fp_env)) fp_env => force_env%fp_env
     IF (PRESENT(meta_env)) meta_env => force_env%meta_env
     IF (PRESENT(sub_force_env)) sub_force_env => force_env%sub_force_env
     IF (PRESENT(qmmm_env))      qmmm_env      => force_env%qmmm_env
     IF (PRESENT(mixed_env))     mixed_env     => force_env%mixed_env
     IF (PRESENT(virial))        virial        => force_env%virial
     IF (PRESENT(globenv))       globenv       => force_env%globenv
     IF (PRESENT(root_section))  root_section  => force_env%root_section
     IF (PRESENT(ep_env))        ep_env        => force_env%ep_env
  END IF
END SUBROUTINE force_env_get

! *****************************************************************************
!> \brief changes some attributes of the force_env
!> \param force_env the force environment where the cell should be changed
!> \param meta_env the new meta environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE force_env_set(force_env, meta_env,fp_env, virial, force_env_section,&
     method_name_id, extended_energies,results, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(meta_env_type), OPTIONAL, POINTER   :: meta_env
    TYPE(fp_type), OPTIONAL, POINTER         :: fp_env
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: force_env_section
    INTEGER, OPTIONAL                        :: method_name_id
    TYPE(extended_energy_type), OPTIONAL, &
      POINTER                                :: extended_energies
    TYPE(cp_result_type), OPTIONAL, POINTER  :: results
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     IF (PRESENT(extended_energies)) THEN
        IF (ASSOCIATED(force_env%extended_energies)) THEN
           DEALLOCATE(force_env%extended_energies,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        force_env%extended_energies => extended_energies
     END IF
     IF (PRESENT(meta_env)) THEN
        IF (ASSOCIATED(meta_env)) THEN
           CALL meta_env_retain(meta_env,error=error)
        END IF
        CALL meta_env_release(force_env%meta_env,error=error)
        force_env%meta_env => meta_env
     END IF
     IF (PRESENT(fp_env)) THEN
        CALL fp_env_retain(fp_env,error=error)
        CALL fp_env_release(force_env%fp_env,error=error)
        force_env%fp_env => fp_env
     END IF
     IF (PRESENT(virial)) THEN
        IF (ASSOCIATED(virial)) THEN
           CALL virial_retain(virial,error=error)
        END IF
        CALL virial_release(force_env%virial,error=error)
        force_env%virial => virial
     END IF
     IF (PRESENT(force_env_section)) THEN
        IF (ASSOCIATED(force_env_section)) THEN
           CALL section_vals_retain(force_env_section,error=error)
           CALL section_vals_release(force_env%force_env_section,error=error)
           force_env%force_env_section => force_env_section
        END IF
     END IF
     IF(PRESENT(results))THEN
        IF(ASSOCIATED(results))THEN
           CALL cp_result_retain(results,error)
        END IF
        CALL cp_result_release(force_env%results,error)
        force_env%results => results
     END IF
     IF (PRESENT(method_name_id))THEN
        force_env%method_name_id=method_name_id
     END IF
  END IF
END SUBROUTINE force_env_set

! *****************************************************************************
!> \brief returns the position in a dimension(*) array
!> \param force_env the force_env you want to get the positions
!> \param pos the array of the positions
!> \param n_el the number of positions (3*natom), just to check
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE force_env_pos_get(force_env,pos,n_el,error)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(kind=dp), DIMENSION(*), INTENT(out) :: pos
    INTEGER, INTENT(in)                      :: n_el
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_pos_get', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idir, ii, ip, my_n_el
    LOGICAL                                  :: failure
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles

  failure=.FALSE.

  CALL timeset(routineN,handle)
  CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     IF (.NOT.failure) THEN
        CALL force_env_get(force_env,&
             subsys=subsys,error=error)
        CALL cp_error_check(error,failure)
     END IF
     IF (.NOT.failure) THEN
        CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
     END IF
     IF (.NOT.failure) THEN
        particles => subsys%particles
        my_n_el=3*particles%n_els
        CALL cp_assert(my_n_el==n_el,cp_failure_level,cp_assertion_failed,&
                       routineP,"wrong pos size ("//cp_to_string(n_el)&
                       //" vs" // &
                       cp_to_string(my_n_el) // &
 CPSourceFileRef,error,failure)
     END IF
     IF (.NOT.failure) THEN
        ii=0
        particles => subsys%particles
        DO ip=1,particles%n_els
           DO idir=1,3
              ii=ii+1
              pos(ii)=particles%els(ip)%r(idir)
           END DO
        END DO
     END IF
  END IF
  CALL timestop(handle)
END SUBROUTINE force_env_pos_get

! *****************************************************************************
!> \brief returns the number of atoms
!> \param force_env the force_env you what information about
!> \param natom the numer of atoms
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
FUNCTION force_env_get_natom(force_env,error) RESULT(n_atom)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: n_atom

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_get_natom', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles

  failure=.FALSE.

  n_atom=0
  IF (.NOT.failure) THEN
     NULLIFY(subsys)
     CALL force_env_get(force_env,&
          subsys=subsys,error=error)
     CALL cp_error_check(error,failure)
     IF (.NOT.failure) THEN
        CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
     END IF
     IF (.NOT.failure) THEN
        particles => subsys%particles
        n_atom=particles%n_els
     END IF
  END IF
END FUNCTION force_env_get_natom

! *****************************************************************************
!> \brief returns the order of the multiple force_env
!> \param force_env the force_env you what information about
!> \param natom the numer of atoms
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
SUBROUTINE multiple_fe_list(force_env_sections, root_section, i_force_eval, nforce_eval, error)
    TYPE(section_vals_type), POINTER         :: force_env_sections, &
                                                root_section
    INTEGER, DIMENSION(:), POINTER           :: i_force_eval
    INTEGER                                  :: nforce_eval
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'multiple_fe_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iforce_eval, main_force_eval, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: my_i_force_eval
    LOGICAL                                  :: failure

! Let's treat the case of Multiple force_eval

  CALL section_vals_get(force_env_sections, n_repetition=nforce_eval, error=error)
  CALL section_vals_val_get(root_section,"MULTIPLE_FORCE_EVALS%FORCE_EVAL_ORDER",&
       i_vals=my_i_force_eval,ignore_required=.TRUE.,error=error)
  ALLOCATE(i_force_eval(nforce_eval),stat=stat)
  CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
  IF (nforce_eval>0) THEN
     IF (nforce_eval==SIZE(my_i_force_eval)) THEN
        i_force_eval                = my_i_force_eval
     ELSE
        ! The difference in the amount of defined force_env MUST be one..
        CPPostcondition(nforce_eval-SIZE(my_i_force_eval)==1,cp_fatal_level,routineP,error,failure)
        DO iforce_eval = 1, nforce_eval
           IF (ANY(my_i_force_eval==iforce_eval)) CYCLE
           main_force_eval = iforce_eval
           EXIT
        END DO
        i_force_eval(1)             = main_force_eval
        i_force_eval(2:nforce_eval) = my_i_force_eval
     END IF
  END IF
END SUBROUTINE multiple_fe_list

END MODULE force_env_types

