!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface for the force calculations
!> \par History
!>      cjm, FEB-20-2001: pass variable box_ref
!>      cjm, SEPT-12-2002: major reorganization
!>      fawzi, APR-12-2003: introduced force_env
!>      cjm, FEB-27-2006: no more box_change
!>      MK, Nov. 2010: new interfaces added and others were updated
!> \author CJM & JGH
! *****************************************************************************
MODULE force_env_types
  USE atprop_types,                    ONLY: atprop_type
  USE cell_types,                      ONLY: cell_type
  USE cp_para_env,                     ONLY: cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_types,                 ONLY: cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_type,&
                                             pack_subsys_particles
  USE eip_environment_types,           ONLY: eip_env_get,&
                                             eip_env_release,&
                                             eip_environment_type
  USE ep_f77,                          ONLY: ep_f77_get_energy,&
                                             ep_f77_get_main_qs_env,&
                                             ep_f77_release
  USE ewald_environment_types,         ONLY: ewald_environment_type
  USE ewald_pw_methods,                ONLY: ewald_pw_grid_change
  USE ewald_pw_types,                  ONLY: ewald_pw_type
  USE f77_blas
  USE fist_energy_types,               ONLY: fist_energy_type
  USE fist_environment_types,          ONLY: fist_env_get,&
                                             fist_env_release,&
                                             fist_env_set,&
                                             fist_environment_type
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
  USE kinds,                           ONLY: dp
  USE metadynamics_types,              ONLY: meta_env_release,&
                                             meta_env_retain,&
                                             meta_env_type
  USE mixed_energy_types,              ONLY: mixed_energy_type
  USE mixed_environment_types,         ONLY: get_mixed_env,&
                                             mixed_env_release,&
                                             mixed_environment_type
  USE qmmm_types,                      ONLY: fist_subsys,&
                                             force_mixing_extended_subsys,&
                                             primary_subsys,&
                                             qmmm_env_qm_release,&
                                             qmmm_env_qm_type,&
                                             qs_subsys
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_release,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_methods,                   ONLY: qs_ks_create
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_env_types'

  INTEGER, PARAMETER, PUBLIC :: use_fist_force = 501, &
                                use_qs_force   = 502, &
                                use_qmmm       = 503, &
                                use_eip_force  = 504, &
                                use_ep_force   = 505, &
                                use_mixed_force= 506

  CHARACTER(LEN=10), DIMENSION(501:506), PARAMETER, PUBLIC :: &
       use_prog_name = (/&
       "FIST  ", &
       "QS    ", &
       "QMMM  ", &
       "EIP   ", &
       "EP    ", &
       "MIXED "/)

  PUBLIC :: force_env_type,&
            force_env_p_type

  PUBLIC :: force_env_retain,&
            force_env_release,&
            force_env_get,&
            force_env_get_natom,&
            force_env_get_nparticle,&
            force_env_get_frc,&
            force_env_get_pos,&
            force_env_get_vel,&
            force_env_set,&
            force_env_set_cell,&
            multiple_fe_list

! *****************************************************************************
!> \brief wrapper to abstract the force evaluation of the various methods
!> \param id_nr identification number
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param in_use which method is in use
!> \param fist_env the fist environment (allocated only if fist is in use)
!> \param qs_env qs_env (activated only if quickstep is in use)
!> \param globenv the globenv to have the input that generated this force_env
!> \param para_env the parallel environment that contains all the parallel
!>        environment of the fragments
!> \param meta_env the metadynamics environment, allocated if there is
!>        metadynamics
!> \param fp_env the flexible partitioning environment
!>      read-only attributes (get them *only* through force_env_get):
!> \param subsys the fragments that build up the actual system.
!> \param ep_env_id environment for an ep calculation
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
     INTEGER :: id_nr,ref_count,in_use,method_name_id
     REAL ( KIND=dp )                                    :: additional_potential
     TYPE ( fist_environment_type ), POINTER             :: fist_env
     TYPE ( meta_env_type ), POINTER                     :: meta_env
     TYPE ( fp_type ), POINTER                           :: fp_env
     TYPE ( qs_environment_type ), POINTER               :: qs_env
     TYPE ( eip_environment_type ), POINTER              :: eip_env
     TYPE ( cp_subsys_type ), POINTER                    :: subsys
     TYPE ( global_environment_type ), POINTER           :: globenv
     TYPE ( cp_para_env_type ), POINTER                  :: para_env
     TYPE ( force_env_p_type ), DIMENSION(:), POINTER    :: sub_force_env
     TYPE ( qmmm_env_qm_type ), POINTER                  :: qmmm_env
     INTEGER, POINTER                                    :: ep_env_id
     TYPE ( mixed_environment_type ), POINTER            :: mixed_env
     TYPE ( section_vals_type), POINTER                  :: force_env_section
     TYPE ( section_vals_type), POINTER                  :: root_section
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
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
!> \note
!>      see doc/ReferenceCounting.html
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
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
!> \note
!>      see doc/ReferenceCounting.html
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

          SELECT CASE ( force_env%in_use )
          CASE ( use_fist_force )
             CALL fist_env_release(force_env%fist_env,error=error)
          CASE ( use_qs_force )
             CALL qs_env_release(force_env%qs_env,error=error)
          CASE ( use_eip_force )
             CALL eip_env_release(force_env%eip_env, error=error)
          CASE ( use_ep_force )
             CALL ep_f77_release(force_env%ep_env_id, error=error)
          CASE (use_mixed_force)
             CALL mixed_env_release(force_env%mixed_env,error=error)
          END SELECT
          CALL globenv_release(force_env%globenv,error=error)
          CALL cp_para_env_release(force_env%para_env,error=error)
          ! Not deallocated
          CPAssert(.NOT.ASSOCIATED(force_env%fist_env),cp_warning_level,routineP,error,failure)
          CPAssert(.NOT.ASSOCIATED(force_env%qs_env),cp_warning_level,routineP,error,failure)
          CPAssert(.NOT.ASSOCIATED(force_env%eip_env),cp_warning_level,routineP,error,failure)
          CPAssert(.NOT.ASSOCIATED(force_env%ep_env_id),cp_warning_level,routineP,error,failure)
          CPAssert(.NOT.ASSOCIATED(force_env%mixed_env),cp_warning_level,routineP,error,failure)
          CALL meta_env_release(force_env%meta_env,error=error)
          CALL fp_env_release(force_env%fp_env,error=error)
          CALL qmmm_env_qm_release(force_env%qmmm_env,error=error)
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
!> \param in_use ...
!> \param fist_env ...
!> \param qs_env ...
!> \param meta_env ...
!> \param fp_env ...
!> \param atprop ...
!> \param subsys ...
!> \param para_env ...
!> \param potential_energy ...
!> \param additional_potential ...
!> \param kinetic_energy ...
!> \param harmonic_shell ...
!> \param kinetic_shell ...
!> \param cell ...
!> \param cell_ref ...
!> \param sub_force_env ...
!> \param qmmm_env ...
!> \param eip_env ...
!> \param virial ...
!> \param use_ref_cell ...
!> \param globenv ...
!> \param input ...
!> \param ep_env_id ...
!> \param force_env_section ...
!> \param method_name_id ...
!> \param root_section ...
!> \param mixed_env ...
!> \param results ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>      for the orther arguments see the attributes of force_env_type
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  RECURSIVE SUBROUTINE force_env_get(force_env, in_use, fist_env, qs_env, &
       meta_env, fp_env, atprop, subsys, para_env, potential_energy, additional_potential, &
       kinetic_energy, harmonic_shell, kinetic_shell, cell, cell_ref, sub_force_env,&
       qmmm_env, eip_env, virial, use_ref_cell, globenv, input, ep_env_id, force_env_section, &
       method_name_id, root_section, mixed_env, results, error)
    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, INTENT(out), OPTIONAL           :: in_use
    TYPE(fist_environment_type), OPTIONAL, &
      POINTER                                :: fist_env
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(meta_env_type), OPTIONAL, POINTER   :: meta_env
    TYPE(fp_type), OPTIONAL, POINTER         :: fp_env
    TYPE(atprop_type), OPTIONAL, POINTER     :: atprop
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: potential_energy, &
      additional_potential, kinetic_energy, harmonic_shell, kinetic_shell
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
    INTEGER, OPTIONAL, POINTER               :: ep_env_id
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

    INTEGER                                  :: cur_subsys
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: eip_kinetic_energy, &
                                                eip_potential_energy, &
                                                penergy_mm, penergy_qm
    TYPE(fist_energy_type), POINTER          :: thermo
    TYPE(mixed_energy_type), POINTER         :: mixed_energy
    TYPE(qs_energy_type), POINTER            :: qs_energy
    TYPE(qs_environment_type), POINTER       :: ep_main_qs_env

    NULLIFY(ep_main_qs_env)

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       SELECT CASE(force_env%in_use)
       CASE(use_ep_force)
          CPPrecondition(ASSOCIATED(force_env%ep_env_id),cp_failure_level,routineP,error,failure)
          CALL ep_f77_get_main_qs_env(force_env%ep_env_id, ep_main_qs_env)
          CALL get_qs_env(ep_main_qs_env,&
                          cell=cell,&
                          cell_ref=cell_ref,&
                          use_ref_cell=use_ref_cell,&
                          input=input,&
                          error=error)
          NULLIFY(ep_main_qs_env)
          IF (PRESENT(potential_energy)) potential_energy = ep_f77_get_energy(force_env%ep_env_id)
          CPPrecondition(.NOT.PRESENT(kinetic_energy),cp_failure_level,routineP,error,failure)
       CASE (use_qs_force)
          CPPrecondition(ASSOCIATED(force_env%qs_env),cp_failure_level,routineP,error,failure)
          CPPreconditionNoFail(.NOT.PRESENT(fist_env),cp_warning_level,routineP,error)
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
!          CPPreconditionNoFail(.NOT.PRESENT(qs_env),cp_warning_level,routineP,error)
!          CPPreconditionNoFail(.NOT.PRESENT(eip_env),cp_warning_level,routineP,error)
          CPPrecondition(.NOT.PRESENT(input),cp_failure_level,routineP,error,failure)
          CALL fist_env_get(force_env%fist_env,&
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
       CASE (use_qmmm)
            IF (PRESENT(cell)) THEN
                 CALL force_env_get(force_env%sub_force_env(primary_subsys)%force_env,&
                                    cell=cell,&
                                    error=error)
            ENDIF
            IF (PRESENT(cell_ref)) THEN
                 CALL force_env_get(force_env%sub_force_env(primary_subsys)%force_env,&
                                    cell_ref=cell_ref,&
                                    error=error)
            ENDIF
            IF (PRESENT(use_ref_cell)) use_ref_cell = .FALSE.
            IF (PRESENT(kinetic_energy)) THEN
               CALL force_env_get(force_env%sub_force_env(primary_subsys)%force_env,&
                                  kinetic_energy=kinetic_energy,&
                                  error=error)
            END IF
            IF (PRESENT(potential_energy)) THEN
                 ! get the underlying energies from primary subsys.  This is the only subsys
                 ! for conventional QM/MM, and force-mixing knows to put relevant energy there.
                 IF (SIZE(force_env%sub_force_env) == 1) THEN
                   cur_subsys = primary_subsys
                 ELSE IF (force_env%qmmm_env%do_force_mixing .AND. SIZE(force_env%sub_force_env) == 2) THEN ! force mixing
                   cur_subsys = force_mixing_extended_subsys
                 ELSE
                    CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                         routineP,"bad number of sub_force_envs (or do_force_mixing setting) getting qmmm potential energy "//&
                         CPSourceFileRef,&
                         error,failure)
                 ENDIF
                 CALL force_env_get(force_env%sub_force_env(cur_subsys)%force_env%sub_force_env(fist_subsys)%force_env,&
                                    potential_energy=penergy_mm,&
                                    error=error)
                 CALL force_env_get(force_env%sub_force_env(cur_subsys)%force_env%sub_force_env(qs_subsys)%force_env,&
                                    potential_energy=penergy_qm,&
                                    error=error)
                 potential_energy = penergy_qm+penergy_mm
            ENDIF
       CASE (use_mixed_force)
          CPPrecondition(ASSOCIATED(force_env%mixed_env),cp_failure_level,routineP,error,failure)
!          CPPreconditionNoFail(.NOT.PRESENT(qs_env),cp_warning_level,routineP,error)
!          CPPreconditionNoFail(.NOT.PRESENT(eip_env),cp_warning_level,routineP,error)
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
       IF (PRESENT(eip_env)) THEN
          eip_env => force_env%eip_env
       END IF
       IF (PRESENT(results)) &
         CALL cp_subsys_get(force_env%subsys, results=results, error=error)
       IF (PRESENT(subsys)) subsys => force_env%subsys
       IF (PRESENT(para_env)) para_env => force_env%para_env
       ! adjust the total energy for the metadynamics
       IF (ASSOCIATED(force_env%meta_env)) THEN
         IF (PRESENT(potential_energy)) THEN
           potential_energy=potential_energy + &
                            force_env%meta_env%epot_s + &
                            force_env%meta_env%epot_walls + &
                            force_env%meta_env%hills_env%energy
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
       IF (PRESENT(potential_energy)) THEN
          potential_energy = potential_energy + force_env%additional_potential
       END IF
       IF (PRESENT(additional_potential)) THEN
          additional_potential = force_env%additional_potential
       END IF
       IF (PRESENT(fp_env)) fp_env => force_env%fp_env
       IF (PRESENT(meta_env)) meta_env => force_env%meta_env
       IF (PRESENT(sub_force_env)) sub_force_env => force_env%sub_force_env
       IF (PRESENT(qmmm_env))      qmmm_env      => force_env%qmmm_env
       IF (PRESENT(mixed_env))     mixed_env     => force_env%mixed_env
       IF (PRESENT(atprop)) &
          CALL cp_subsys_get(force_env%subsys, atprop=atprop, error=error)
       IF (PRESENT(virial)) &
          CALL cp_subsys_get(force_env%subsys, virial=virial, error=error)
       IF (PRESENT(globenv))       globenv       => force_env%globenv
       IF (PRESENT(root_section))  root_section  => force_env%root_section
       IF (PRESENT(ep_env_id))     ep_env_id        => force_env%ep_env_id
    END IF

  END SUBROUTINE force_env_get

! *****************************************************************************
!> \brief returns the number of atoms
!> \param force_env the force_env you what information about
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval n_atom the number of atoms
!> \date   22.11.2010 updated (MK)
!> \author fawzi
! *****************************************************************************
  FUNCTION force_env_get_natom(force_env,error) RESULT(n_atom)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: n_atom

    CHARACTER(LEN=*), PARAMETER :: routineN = 'force_env_get_natom', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_subsys_type), POINTER            :: subsys

    failure = .FALSE.
    n_atom = 0
    NULLIFY (subsys)
    CALL force_env_get(force_env,subsys=subsys,error=error)
    CALL cp_error_check(error,failure)
    IF (.NOT.failure) THEN
       CALL cp_subsys_get(subsys,natom=n_atom,error=error)
    END IF

  END FUNCTION force_env_get_natom

! *****************************************************************************
!> \brief returns the number of particles in a force environment
!> \param force_env the force_env you what information about
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval n_particle the number of particles
!> \date   22.11.2010 (MK)
!> \author Matthias Krack
! *****************************************************************************
  FUNCTION force_env_get_nparticle(force_env,error) RESULT(n_particle)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: n_particle

    CHARACTER(LEN=*), PARAMETER :: routineN = 'force_env_get_nparticle', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_subsys_type), POINTER            :: subsys

    failure = .FALSE.
    n_particle = 0
    NULLIFY (subsys)
    CALL force_env_get(force_env,subsys=subsys,error=error)
    CALL cp_error_check(error,failure)
    IF (.NOT.failure) THEN
      CALL cp_subsys_get(subsys,nparticle=n_particle,error=error)
    END IF

  END FUNCTION force_env_get_nparticle

! *****************************************************************************
!> \brief returns the particle forces in a dimension(*) array
!> \param force_env the force_env you want to get the forces
!> \param frc the array of the forces
!> \param n ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \date   22.11.2010 Creation
!> \author Matthias Krack
! *****************************************************************************
  SUBROUTINE force_env_get_frc(force_env,frc,n,error)

    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(*), INTENT(OUT) :: frc
    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'force_env_get_frc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_subsys_type), POINTER            :: subsys

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CPPrecondition((force_env%ref_count > 0),cp_failure_level,routineP,error,failure)
    CALL force_env_get(force_env,subsys=subsys,error=error)
    CALL cp_error_check(error,failure)
    IF (.NOT.failure) THEN
       CALL pack_subsys_particles(subsys=subsys,f=frc(1:n),error=error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE force_env_get_frc

! *****************************************************************************
!> \brief returns the particle positions in a dimension(*) array
!> \param force_env the force_env you want to get the positions
!> \param pos the array of the positions
!> \param n ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \date   22.11.2010 updated (MK)
!> \author fawzi
! *****************************************************************************
  SUBROUTINE force_env_get_pos(force_env,pos,n,error)

    TYPE(force_env_type), POINTER            :: force_env
    REAL(kind=dp), DIMENSION(*), INTENT(OUT) :: pos
    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'force_env_get_pos', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_subsys_type), POINTER            :: subsys

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CPPrecondition((force_env%ref_count > 0),cp_failure_level,routineP,error,failure)
    CALL force_env_get(force_env,subsys=subsys,error=error)
    CALL cp_error_check(error,failure)
    IF (.NOT.failure) THEN
       CALL pack_subsys_particles(subsys=subsys,r=pos(1:n),error=error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE force_env_get_pos

! *****************************************************************************
!> \brief returns the particle velocities in a dimension(*) array
!> \param force_env the force_env you want to get the velocities
!> \param vel the array of the velocities
!> \param n ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \date   22.11.2010 Creation (MK)
!> \author Matthias Krack
! *****************************************************************************
  SUBROUTINE force_env_get_vel(force_env,vel,n,error)

    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(*), INTENT(OUT) :: vel
    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'force_env_get_vel', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_subsys_type), POINTER            :: subsys

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CPPrecondition((force_env%ref_count > 0),cp_failure_level,routineP,error,failure)
    CALL force_env_get(force_env,subsys=subsys,error=error)
    CALL cp_error_check(error,failure)
    IF (.NOT.failure) THEN
       CALL pack_subsys_particles(subsys=subsys,v=vel(1:n),error=error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE force_env_get_vel

! *****************************************************************************
!> \brief changes some attributes of the force_env
!> \param force_env the force environment where the cell should be changed
!> \param meta_env the new meta environment
!> \param fp_env ...
!> \param force_env_section ...
!> \param method_name_id ...
!> \param additional_potential ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param results ... UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par History
!>      09.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE force_env_set(force_env, meta_env,fp_env, force_env_section,&
     method_name_id, additional_potential, error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(meta_env_type), OPTIONAL, POINTER   :: meta_env
    TYPE(fp_type), OPTIONAL, POINTER         :: fp_env
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: force_env_section
    INTEGER, OPTIONAL                        :: method_name_id
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: additional_potential
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_env_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
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
       IF (PRESENT(fp_env)) THEN
          CALL fp_env_retain(fp_env,error=error)
          CALL fp_env_release(force_env%fp_env,error=error)
          force_env%fp_env => fp_env
       END IF
       IF (PRESENT(force_env_section)) THEN
          IF (ASSOCIATED(force_env_section)) THEN
             CALL section_vals_retain(force_env_section,error=error)
             CALL section_vals_release(force_env%force_env_section,error=error)
             force_env%force_env_section => force_env_section
          END IF
       END IF
       IF (PRESENT(additional_potential)) THEN
          force_env%additional_potential = additional_potential
       END IF
       IF (PRESENT(method_name_id)) THEN
          force_env%method_name_id=method_name_id
       END IF
    END IF

  END SUBROUTINE force_env_set

! *****************************************************************************
!> \brief changes the cell of the system
!> \param force_env the force environment where the cell should be changed
!> \param cell the new cell
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2003 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      for future use, each method should make the changes that it needs
!>      to cope with the new cell
! *****************************************************************************
  RECURSIVE SUBROUTINE force_env_set_cell(force_env, cell, error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'force_env_set_cell', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iforce_eval, my_group, &
                                                nforce_eval
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: my_error
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(qs_environment_type), POINTER       :: ep_main_qs_env
    TYPE(qs_ks_env_type), POINTER            :: new_ks_env

    NULLIFY(ep_main_qs_env)

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(force_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)
    CPPrecondition(cell%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       SELECT CASE(force_env%in_use)
       CASE (use_qs_force)
          NULLIFY ( new_ks_env )
          CPPrecondition(ASSOCIATED(force_env%qs_env),cp_failure_level,routineP,error,failure)
          CALL set_qs_env(force_env%qs_env,cell=cell,error=error)
          CALL qs_env_rebuild_pw_env( force_env % qs_env, error = error)
          CALL qs_ks_create ( new_ks_env,  force_env % qs_env ,error = error )
          CALL set_qs_env( force_env % qs_env, ks_env=new_ks_env ,error=error)
          CALL qs_ks_release ( new_ks_env ,error=error)
       CASE (use_ep_force)
          NULLIFY ( new_ks_env )
          CPPrecondition(ASSOCIATED(force_env%ep_env_id),cp_failure_level,routineP,error,failure)
          CALL ep_f77_get_main_qs_env(force_env%ep_env_id, ep_main_qs_env)
          CALL set_qs_env(ep_main_qs_env,cell=cell,error=error)
          CALL qs_env_rebuild_pw_env(ep_main_qs_env, error=error)
          CALL qs_ks_create ( new_ks_env,  force_env % qs_env ,error=error)
          CALL set_qs_env(ep_main_qs_env, ks_env=new_ks_env ,error=error)
          CALL qs_ks_release ( new_ks_env ,error=error)
          NULLIFY(ep_main_qs_env)
       CASE (use_fist_force)
          CPPrecondition(ASSOCIATED(force_env%fist_env),cp_failure_level,routineP,error,failure)
          CALL fist_env_get ( force_env%fist_env, ewald_pw = ewald_pw,  &
                              ewald_env = ewald_env ,error=error)
          CALL fist_env_set(force_env%fist_env, cell=cell, error=error)
          CALL ewald_pw_grid_change ( ewald_pw, ewald_env, cell, error )
       CASE (use_eip_force)
          CPAssert(.FALSE., cp_failure_level, routineP, error, failure)
       CASE (use_mixed_force)
          CPPrecondition(ASSOCIATED(force_env%mixed_env),cp_failure_level,routineP,error,failure)
          nforce_eval = SIZE(force_env%sub_force_env)
          my_group = force_env%mixed_env%group_distribution(force_env%para_env%mepos)
          my_error = force_env%mixed_env%sub_error(my_group+1)
          DO iforce_eval = 1, nforce_eval
             IF (ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) THEN
                CALL force_env_set_cell(force_env=force_env%sub_force_env(iforce_eval)%force_env,&
                     cell=cell, error=my_error)
             END IF
          END DO
       CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown in_use flag value "//&
CPSourceFileRef,&
               error,failure)
       END SELECT
    END IF

  END SUBROUTINE force_env_set_cell

! *****************************************************************************
!> \brief returns the order of the multiple force_env
!> \param force_env_sections ...
!> \param root_section ...
!> \param i_force_eval ...
!> \param nforce_eval ...
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

    failure = .FALSE.
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

