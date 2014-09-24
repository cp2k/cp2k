!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE input_restart_force_eval
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cell_methods,                    ONLY: cell_type,&
                                             real_to_scaled
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_linked_list_val,              ONLY: cp_sll_val_create,&
                                             cp_sll_val_get_length,&
                                             cp_sll_val_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             multiple_fe_list,&
                                             use_qmmm
  USE input_constants,                 ONLY: do_coord_cp2k,&
                                             ehrenfest,&
                                             mol_dyn_run,&
                                             mon_car_run,&
                                             pint_run,&
                                             use_rt_restart
  USE input_cp2k_restarts_util,        ONLY: section_velocity_val_set
  USE input_restart_rng,               ONLY: section_rng_val_set
  USE input_section_types,             ONLY: &
       section_get_keyword_index, section_type, section_vals_add_values, &
       section_vals_get, section_vals_get_subs_vals, &
       section_vals_get_subs_vals3, section_vals_remove_values, &
       section_vals_type, section_vals_val_get, section_vals_val_set, &
       section_vals_val_unset
  USE input_val_types,                 ONLY: val_create,&
                                             val_release,&
                                             val_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_sum
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: get_molecule_kind
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE multipole_types,                 ONLY: multipole_type
  USE parallel_rng_types,              ONLY: dump_rng_stream,&
                                             rng_record_length
  USE particle_list_types,             ONLY: particle_list_type
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE string_utilities,                ONLY: string_to_ascii
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_restart_force_eval'

  PUBLIC ::           update_force_eval, update_subsys, section_rng_val_set

CONTAINS

! *****************************************************************************
!> \brief Updates the force_eval section of the input file
!> \param force_env ...
!> \param root_section ...
!> \param write_binary_restart_file ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_force_eval(force_env,root_section,&
                               write_binary_restart_file,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    LOGICAL, INTENT(IN)                      :: write_binary_restart_file
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_force_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iforce_eval, myid, nforce_eval
    INTEGER, DIMENSION(:), POINTER           :: i_force_eval
    LOGICAL                                  :: failure, multiple_subsys, &
                                                skip_vel_section
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(section_vals_type), POINTER         :: force_env_sections, &
                                                qmmm_section, rng_section, &
                                                subsys_section
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.
    NULLIFY(rng_section, subsys_section, virial, subsys)
    ! If it's not a dynamical run don't update the velocity section
    CALL section_vals_val_get(root_section,"GLOBAL%RUN_TYPE",i_val=myid,error=error)
    skip_vel_section =((myid /= mol_dyn_run).AND.&
                       (myid /= mon_car_run).AND.&
                       (myid /= pint_run).AND.&
                       (myid /= ehrenfest))

    ! Go on updatig the force_env_section
    force_env_sections => section_vals_get_subs_vals(root_section,"FORCE_EVAL",error=error)
    CALL multiple_fe_list(force_env_sections, root_section, i_force_eval, nforce_eval, error)
    ! The update of the input MUST be realized only on the main force_eval
    ! All the others will be left not updated because there is no real need to update them...
    iforce_eval=1
    subsys_section => section_vals_get_subs_vals3(force_env_sections,"SUBSYS",&
         i_rep_section=i_force_eval(iforce_eval),error=error)
    CALL update_subsys(subsys_section,force_env,skip_vel_section,&
                       write_binary_restart_file,error)

    rng_section    => section_vals_get_subs_vals(subsys_section,"RNG_INIT",error=error)
    CALL update_rng_particle(rng_section,force_env,error)

    qmmm_section   => section_vals_get_subs_vals3(force_env_sections,"QMMM",&
         i_rep_section=i_force_eval(iforce_eval),error=error)
    CALL update_qmmm(qmmm_section,force_env,error)

    ! And now update only the cells of all other force_eval
    ! This is to make consistent for cell variable runs..
    IF (nforce_eval>1) THEN
       CALL force_env_get(force_env, subsys=subsys, error=error)
       CALL cp_subsys_get(subsys, virial=virial, error=error)
       CALL section_vals_val_get(root_section,"MULTIPLE_FORCE_EVALS%MULTIPLE_SUBSYS",l_val=multiple_subsys,error=error)
       IF (virial%pv_availability.AND.multiple_subsys) THEN
          DO iforce_eval = 2, nforce_eval
             subsys_section => section_vals_get_subs_vals3(force_env_sections,"SUBSYS",&
                  i_rep_section=i_force_eval(iforce_eval),error=error)
             CALL update_cell_section(force_env, subsys_section, error)
          END DO
       END IF
    END IF

    IF(myid == ehrenfest)CALL section_vals_val_set(root_section,"FORCE_EVAL%DFT%REAL_TIME_PROPAGATION%INITIAL_WFN",&
                                                   i_val=use_rt_restart,error=error)
    DEALLOCATE(i_force_eval)

  END SUBROUTINE update_force_eval

! *****************************************************************************
!> \brief Updates the qmmm section if needed
!> \param qmmm_section ...
!> \param force_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2007 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_qmmm(qmmm_section,force_env,error)
    TYPE(section_vals_type), POINTER         :: qmmm_section
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_qmmm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work

    failure = .FALSE.
    SELECT CASE(force_env%in_use)
    CASE(use_qmmm)
       CALL section_vals_get(qmmm_section, explicit=explicit, error=error)
       CPPostcondition(explicit,cp_failure_level,routineP,error,failure)

       ALLOCATE(work(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       work = force_env%qmmm_env%transl_v
       CALL section_vals_val_set(qmmm_section,"INITIAL_TRANSLATION_VECTOR",r_vals_ptr=work,error=error)
    END SELECT

  END SUBROUTINE update_qmmm

! *****************************************************************************
!> \brief Updates the rng section of the input file
!>      Write current status of the parallel random number generator (RNG)
!> \param rng_section ...
!> \param force_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_rng_particle(rng_section,force_env,error)

    TYPE(section_vals_type), POINTER         :: rng_section
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_rng_particle', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=rng_record_length)         :: rng_record
    INTEGER :: iparticle, iparticle_kind, iparticle_local, istat, nparticle, &
      nparticle_kind, nparticle_local
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: ascii
    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_list_type), POINTER        :: particles

    CALL force_env_get(force_env, subsys=subsys, para_env=para_env, error=error)
    CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, local_particles=local_particles,&
         particles=particles, error=error)

    IF (ASSOCIATED(local_particles%local_particle_set)) THEN
       nparticle_kind = atomic_kinds%n_els
       nparticle = particles%n_els

       ALLOCATE (ascii(rng_record_length,nparticle),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ascii = 0

       DO iparticle=1,nparticle
          DO iparticle_kind=1,nparticle_kind
             nparticle_local = local_particles%n_el(iparticle_kind)
             DO iparticle_local=1,nparticle_local
                IF (iparticle == local_particles%list(iparticle_kind)%array(iparticle_local)) THEN
                   CALL dump_rng_stream(rng_stream=local_particles%local_particle_set(iparticle_kind)%&
                         rng(iparticle_local)%stream,&
                        rng_record=rng_record, error=error)
                   CALL string_to_ascii(rng_record,ascii(:,iparticle))
                END IF
             END DO
          END DO
       END DO

       CALL mp_sum(ascii,para_env%group)

       CALL section_rng_val_set(rng_section=rng_section,nsize=nparticle,ascii=ascii,&
            error=error)

       DEALLOCATE (ascii,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE update_rng_particle

! *****************************************************************************
!> \brief Updates the subsys section of the input file
!> \param subsys_section ...
!> \param force_env ...
!> \param skip_vel_section ...
!> \param write_binary_restart_file ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_subsys(subsys_section,force_env,skip_vel_section,&
                           write_binary_restart_file,error)
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(IN)                      :: skip_vel_section, &
                                                write_binary_restart_file
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_subsys', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: coord_file_name, unit_str
    INTEGER                                  :: coord_file_format, handle, &
                                                output_unit, stat
    INTEGER, DIMENSION(:), POINTER           :: multiple_unit_cell
    LOGICAL                                  :: failure, scale
    REAL(KIND=dp)                            :: conv_factor
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(multipole_type), POINTER            :: multipoles
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(section_vals_type), POINTER         :: work_section

    failure = .FALSE.
    NULLIFY (work_section, core_particles, particles, shell_particles, &
             subsys, cell, para_env, multipoles)
    CALL timeset(routineN,handle)
    CALL force_env_get(force_env, subsys=subsys, cell=cell, para_env=para_env, &
                       error=error)

    CALL cp_subsys_get(subsys, particles=particles, molecules_new=molecules,&
         shell_particles=shell_particles, core_particles=core_particles,&
         multipoles=multipoles, error=error)

    ! Remove the multiple_unit_cell from the input structure.. at this point we have
    ! already all the information we need..
    ALLOCATE(multiple_unit_cell(3),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    multiple_unit_cell = 1
    CALL section_vals_val_set(subsys_section, "TOPOLOGY%MULTIPLE_UNIT_CELL",&
         i_vals_ptr=multiple_unit_cell, error=error)

    ! Coordinates and Velocities
    work_section => section_vals_get_subs_vals(subsys_section,"COORD",error=error)
    CALL section_vals_val_get(work_section,"UNIT",c_val=unit_str,ignore_required=.TRUE.,error=error)
    CALL section_vals_val_get(work_section,"SCALED",l_val=scale,ignore_required=.TRUE.,error=error)
    conv_factor = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
    IF (.NOT.write_binary_restart_file) THEN
       CALL section_vals_val_get(subsys_section,"TOPOLOGY%COORD_FILE_FORMAT",&
                                 i_val=coord_file_format,&
                                 error=error)
       IF (coord_file_format == do_coord_cp2k) THEN
          CALL section_vals_val_get(subsys_section,"TOPOLOGY%COORD_FILE_NAME",&
                                    c_val=coord_file_name,&
                                    error=error)
          output_unit = 0
          IF (para_env%ionode) THEN
             CALL open_file(file_name=TRIM(ADJUSTL(coord_file_name)),&
                            file_action="READWRITE",&
                            file_form="FORMATTED",&
                            file_position="REWIND",&
                            file_status="UNKNOWN",&
                            unit_number=output_unit)
            CALL dump_coordinates_cp2k(particles,molecules,cell,conv_factor,&
                                       output_unit=output_unit,&
                                       core_or_shell=.FALSE.,&
                                       scaled_coordinates=scale,&
                                       error=error)
            CALL close_file(unit_number=output_unit)
          END IF
       ELSE
          CALL section_coord_val_set(work_section, particles, molecules, conv_factor, scale,&
                                     cell, error=error)
       END IF
    END IF
    CALL section_vals_val_set(subsys_section,"TOPOLOGY%NUMBER_OF_ATOMS",i_val=particles%n_els,error=error)
    work_section => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
    IF (.NOT.skip_vel_section) THEN
       IF (.NOT.write_binary_restart_file) THEN
          CALL section_velocity_val_set(work_section, particles, conv_factor=1.0_dp, error=error)
       END IF
    ELSE
       CALL section_vals_remove_values(work_section, error)
    END IF

    ! Write restart input for core-shell model
    IF (.NOT.write_binary_restart_file) THEN
       IF (ASSOCIATED(shell_particles)) THEN
          work_section => section_vals_get_subs_vals(subsys_section,"SHELL_COORD",error=error)
          CALL section_vals_val_get(work_section,"UNIT",c_val=unit_str,ignore_required=.TRUE.,error=error)
          CALL section_vals_val_get(work_section,"SCALED",l_val=scale,ignore_required=.TRUE.,error=error)
          conv_factor = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
          CALL section_coord_val_set(work_section, shell_particles, molecules, &
                                     conv_factor, scale, cell, shell=.TRUE., error=error)
          IF (.NOT.skip_vel_section) THEN
             work_section => section_vals_get_subs_vals(subsys_section,"SHELL_VELOCITY",error=error)
             CALL section_velocity_val_set(work_section, shell_particles, conv_factor=1.0_dp, error=error)
          END IF
       END IF
       IF (ASSOCIATED(core_particles)) THEN
          work_section => section_vals_get_subs_vals(subsys_section,"CORE_COORD",error=error)
          CALL section_vals_val_get(work_section,"UNIT",c_val=unit_str,ignore_required=.TRUE.,error=error)
          CALL section_vals_val_get(work_section,"SCALED",l_val=scale,ignore_required=.TRUE.,error=error)
          conv_factor = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
          CALL section_coord_val_set(work_section, core_particles, molecules, &
                                     conv_factor, scale, cell, shell=.TRUE., error=error)
          IF (.NOT.skip_vel_section) THEN
             work_section => section_vals_get_subs_vals(subsys_section,"CORE_VELOCITY",error=error)
             CALL section_velocity_val_set(work_section, core_particles, conv_factor=1.0_dp, error=error)
          END IF
       END IF
    END IF

    ! Updating cell info
    CALL update_cell_section(force_env, subsys_section, error)
    ! Updating multipoles
    IF (ASSOCIATED(multipoles)) THEN
       work_section => section_vals_get_subs_vals(subsys_section,"MULTIPOLES",error=error)
       DO
          IF (SIZE(work_section%values,2)==1) EXIT
          CALL section_vals_add_values(work_section,error=error)
       END DO
       IF (ASSOCIATED(multipoles%dipoles)) THEN
          work_section => section_vals_get_subs_vals(subsys_section,"MULTIPOLES%DIPOLES",error=error)
          CALL update_dipoles_section(multipoles%dipoles, work_section, error)
       END IF
       IF (ASSOCIATED(multipoles%quadrupoles)) THEN
          work_section => section_vals_get_subs_vals(subsys_section,"MULTIPOLES%QUADRUPOLES",error=error)
          CALL update_quadrupoles_section(multipoles%quadrupoles, work_section, error)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE update_subsys

! *****************************************************************************
!> \brief routine to update cell section.. splitted from the main update subsys_section
!> \param force_env ...
!> \param subsys_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2007 created [teo] Zurich University
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_cell_section(force_env, subsys_section, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_cell_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    INTEGER, DIMENSION(:), POINTER           :: iwork
    LOGICAL                                  :: failure, use_ref_cell
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(section_vals_type), POINTER         :: work_section

    failure = .FALSE.
    NULLIFY(work, work_section, cell, cell_ref, iwork)
    CALL timeset(routineN,handle)
    ! Handle case for parallel runs
    IF (ASSOCIATED(force_env)) THEN
       CALL force_env_get(force_env,&
                          cell=cell,&
                          cell_ref=cell_ref,&
                          use_ref_cell=use_ref_cell,&
                          error=error)
    END IF

    ! CELL
    work_section => section_vals_get_subs_vals(subsys_section,"CELL",error=error)
    ALLOCATE(work(3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! CELL VECTORS - A
    work(1:3) = cell%hmat(1:3,1)
    CALL section_vals_val_set(work_section,"A",r_vals_ptr=work,error=error)
    ALLOCATE(work(3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! CELL VECTORS - B
    work(1:3) = cell%hmat(1:3,2)
    CALL section_vals_val_set(work_section,"B",r_vals_ptr=work,error=error)
    ALLOCATE(work(3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! CELL VECTORS - C
    work(1:3) = cell%hmat(1:3,3)
    CALL section_vals_val_set(work_section,"C",r_vals_ptr=work,error=error)
    ! MULTIPLE_UNIT_CELL
    ALLOCATE(iwork(3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    iwork = 1
    CALL section_vals_val_set(work_section,"MULTIPLE_UNIT_CELL",i_vals_ptr=iwork,error=error)
    ! Unset unused or misleading information
    CALL section_vals_val_unset(work_section,"ABC",error=error)
    CALL section_vals_val_unset(work_section,"ALPHA_BETA_GAMMA",error=error)
    ! CELL_REF
    IF (use_ref_cell.AND.ASSOCIATED(cell_ref)) THEN
       work_section => section_vals_get_subs_vals(subsys_section,"CELL%CELL_REF",error=error)
       ALLOCATE(work(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! CELL VECTORS - A
       work(1:3) = cell_ref%hmat(1:3,1)
       CALL section_vals_val_set(work_section,"A",r_vals_ptr=work,error=error)
       ALLOCATE(work(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! CELL VECTORS - B
       work(1:3) = cell_ref%hmat(1:3,2)
       CALL section_vals_val_set(work_section,"B",r_vals_ptr=work,error=error)
       ALLOCATE(work(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! CELL VECTORS - C
       work(1:3) = cell_ref%hmat(1:3,3)
       CALL section_vals_val_set(work_section,"C",r_vals_ptr=work,error=error)
       ! MULTIPLE_UNIT_CELL
       ALLOCATE(iwork(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       iwork = 1
       CALL section_vals_val_set(work_section,"MULTIPLE_UNIT_CELL",i_vals_ptr=iwork,error=error)
       ! Unset unused or misleading information
       CALL section_vals_val_unset(work_section,"ABC",error=error)
       CALL section_vals_val_unset(work_section,"ALPHA_BETA_GAMMA",error=error)
    END IF
    CALL timestop(handle)
  END SUBROUTINE update_cell_section

! *****************************************************************************
!> \brief routine to dump coordinates.. fast implementation
!> \param coord_section ...
!> \param particles ...
!> \param molecules ...
!> \param conv_factor ...
!> \param scale ...
!> \param cell ...
!> \param shell ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE section_coord_val_set(coord_section, particles, molecules, conv_factor,&
       scale, cell, shell, error)
    TYPE(section_vals_type), POINTER         :: coord_section
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(mol_new_list_type), POINTER         :: molecules
    REAL(KIND=dp)                            :: conv_factor
    LOGICAL, INTENT(IN)                      :: scale
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL, INTENT(IN), OPTIONAL            :: shell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'section_coord_val_set', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2*default_string_length)   :: line
    CHARACTER(LEN=default_string_length)     :: my_tag, name
    INTEGER                                  :: handle, ik, imol, irk, &
                                                last_atom, Nlist
    LOGICAL                                  :: failure, ldum, &
                                                molname_generated, my_shell
    REAL(KIND=dp), DIMENSION(3)              :: r, s
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(molecule_type), POINTER             :: molecule_now
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    my_shell = .FALSE.
    IF(PRESENT(shell)) my_shell = shell
    CPPrecondition(ASSOCIATED(coord_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(coord_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => coord_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)

    DO
       IF (SIZE(coord_section%values,2)==1) EXIT
       CALL section_vals_add_values(coord_section,error=error)
    END DO
    vals  => coord_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    imol=0
    last_atom=0
    DO irk=1,particles%n_els
       CALL get_atomic_kind(particles%els(irk)%atomic_kind,name=name)
       IF (my_shell) THEN
         s = particles%els(irk)%r
         IF (scale) THEN
            r = s
            CALL real_to_scaled(s,r,cell)
         ELSE
            s = s*conv_factor
         END IF
         WRITE (UNIT=line,FMT="(A,3(1X,ES25.16),1X,I0)")&
           TRIM(name),s(1:3),particles%els(irk)%atom_index
         CALL val_create(my_val,lc_val=line,error=error)
         IF (Nlist /= 0) THEN
            IF (irk==1) THEN
               new_pos => vals
            ELSE
               new_pos => new_pos%rest
            END IF
            old_val => new_pos%first_el
            CALL val_release(old_val,error=error)
            new_pos%first_el => my_val
         ELSE
            IF (irk==1) THEN
               NULLIFY(new_pos)
               CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
               vals => new_pos
            ELSE
               NULLIFY(new_pos%rest)
               CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
               new_pos => new_pos%rest
            END IF
         END IF
         NULLIFY(my_val)
       ELSE
         IF(last_atom<irk) THEN
           imol=imol+1
           molecule_now => molecules%els(imol)
           CALL get_molecule(molecule_now,last_atom=last_atom)
           CALL get_molecule_kind(molecule_now%molecule_kind,molname_generated=molname_generated,&
                name=my_tag)
           IF (molname_generated) my_tag=""
         END IF
         ldum = qmmm_ff_precond_only_qm(my_tag)
         ldum = qmmm_ff_precond_only_qm(name)
         s = particles%els(irk)%r
         IF (scale) THEN
            r = s
            CALL real_to_scaled(s,r,cell)
         ELSE
            s = s*conv_factor
         END IF
         IF (LEN_TRIM(my_tag) > 0) THEN
            WRITE (UNIT=line,FMT="(A,3(1X,ES25.16),1X,A,1X,I0)")&
              TRIM(name),s(1:3),TRIM(my_tag),imol
         ELSE
            WRITE (UNIT=line,FMT="(A,3(1X,ES25.16))")&
              TRIM(name),s(1:3)
         END IF
         CALL val_create(my_val,lc_val=line,error=error)

         IF (Nlist /= 0) THEN
            IF (irk==1) THEN
               new_pos => vals
            ELSE
               new_pos => new_pos%rest
            END IF
            old_val => new_pos%first_el
            CALL val_release(old_val,error=error)
            new_pos%first_el => my_val
         ELSE
            IF (irk==1) THEN
               NULLIFY(new_pos)
               CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
               vals => new_pos
            ELSE
               NULLIFY(new_pos%rest)
               CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
               new_pos => new_pos%rest
            END IF
         END IF
         NULLIFY(my_val)
       END IF
    END DO
    coord_section%values(ik,1)%list => vals

    CALL timestop(handle)
  END SUBROUTINE section_coord_val_set

! *****************************************************************************
!> \brief routine to dump dipoles.. fast implementation
!> \param dipoles ...
!> \param dipoles_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2007 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_dipoles_section(dipoles, dipoles_section, error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dipoles
    TYPE(section_vals_type), POINTER         :: dipoles_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_dipoles_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ik, irk, Nlist, &
                                                nloop, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    CPPrecondition(ASSOCIATED(dipoles_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(dipoles_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => dipoles_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)

    ! At least one of the two arguments must be present..
    nloop = SIZE(dipoles,2)
    DO
       IF (SIZE(dipoles_section%values,2)==1) EXIT
       CALL section_vals_add_values(dipoles_section,error=error)
    END DO
    vals  => dipoles_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    DO irk=1,nloop
       ALLOCATE(work(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Always stored in A.U.
       work = dipoles(1:3,irk)
       CALL val_create(my_val,r_vals_ptr=work,error=error)

       IF (Nlist /= 0) THEN
          IF (irk==1) THEN
             new_pos => vals
          ELSE
             new_pos => new_pos%rest
          END IF
          old_val => new_pos%first_el
          CALL val_release(old_val,error=error)
          new_pos%first_el => my_val
       ELSE
          IF (irk==1) THEN
             NULLIFY(new_pos)
             CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
             vals => new_pos
          ELSE
             NULLIFY(new_pos%rest)
             CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
             new_pos => new_pos%rest
          END IF
       END IF
       NULLIFY(my_val)
    END DO
    dipoles_section%values(ik,1)%list => vals

    CALL timestop(handle)
  END SUBROUTINE update_dipoles_section

! *****************************************************************************
!> \brief routine to dump quadrupoles.. fast implementation
!> \param quadrupoles ...
!> \param quadrupoles_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2007 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_quadrupoles_section(quadrupoles, quadrupoles_section, error)

    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: quadrupoles
    TYPE(section_vals_type), POINTER         :: quadrupoles_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_quadrupoles_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ik, ind, irk, j, &
                                                Nlist, nloop, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    CPPrecondition(ASSOCIATED(quadrupoles_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(quadrupoles_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => quadrupoles_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)

    ! At least one of the two arguments must be present..
    nloop = SIZE(quadrupoles,2)
    DO
       IF (SIZE(quadrupoles_section%values,2)==1) EXIT
       CALL section_vals_add_values(quadrupoles_section,error=error)
    END DO
    vals  => quadrupoles_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    DO irk=1,nloop
       ALLOCATE(work(6),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Always stored in A.U.
       ind = 0
       DO i = 1, 3
          DO j = i, 3
             ind  = ind + 1
             work(ind) = quadrupoles(j,i,irk)
          END DO
       END DO
       CALL val_create(my_val,r_vals_ptr=work,error=error)

       IF (Nlist /= 0) THEN
          IF (irk==1) THEN
             new_pos => vals
          ELSE
             new_pos => new_pos%rest
          END IF
          old_val => new_pos%first_el
          CALL val_release(old_val,error=error)
          new_pos%first_el => my_val
       ELSE
          IF (irk==1) THEN
             NULLIFY(new_pos)
             CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
             vals => new_pos
          ELSE
             NULLIFY(new_pos%rest)
             CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
             new_pos => new_pos%rest
          END IF
       END IF
       NULLIFY(my_val)
    END DO
    quadrupoles_section%values(ik,1)%list => vals

    CALL timestop(handle)
  END SUBROUTINE update_quadrupoles_section

! *****************************************************************************
!> \brief   Dump atomic, core, or shell coordinates to a file in CP2K &COORD
!>          section format
!> \param particles ...
!> \param molecules ...
!> \param cell ...
!> \param conv_factor ...
!> \param output_unit ...
!> \param core_or_shell ...
!> \param scaled_coordinates ...
!> \param error ...
!> \par History
!>      07.02.2011 (Creation, MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE dump_coordinates_cp2k(particles,molecules,cell,conv_factor,&
                                   output_unit,core_or_shell,&
                                   scaled_coordinates,error)

    TYPE(particle_list_type), POINTER        :: particles
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(IN)                :: conv_factor
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(IN)                      :: core_or_shell, &
                                                scaled_coordinates
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dump_coordinates_cp2k', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: kind_name, tag
    INTEGER                                  :: handle, imolecule, iparticle, &
                                                last_atom
    LOGICAL                                  :: failure, ldum, &
                                                molname_generated
    REAL(KIND=dp), DIMENSION(3)              :: r, s
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(particles),cp_failure_level,routineP,error,failure)
    IF (.NOT.core_or_shell) THEN
       CPPrecondition(ASSOCIATED(molecules),cp_failure_level,routineP,error,failure)
    END IF
    IF (scaled_coordinates) THEN
       CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)
    END IF

    kind_name = ""
    tag = ""
    imolecule = 0
    last_atom = 0
    DO iparticle=1,particles%n_els
       CALL get_atomic_kind(particles%els(iparticle)%atomic_kind,name=kind_name)
       IF (.NOT.core_or_shell) THEN
          IF (iparticle > last_atom) THEN
             imolecule = imolecule + 1
             molecule => molecules%els(imolecule)
             CALL get_molecule(molecule,last_atom=last_atom)
             CALL get_molecule_kind(molecule%molecule_kind,&
                                    molname_generated=molname_generated,&
                                    name=tag)
             IF (molname_generated) tag = ""
          END IF
          ldum = qmmm_ff_precond_only_qm(tag)
          ldum = qmmm_ff_precond_only_qm(kind_name)
       END IF
       IF (scaled_coordinates) THEN
          CALL real_to_scaled(s,particles%els(iparticle)%r,cell)
          r(1:3) = s(1:3)
       ELSE
          r(1:3) = particles%els(iparticle)%r(1:3)*conv_factor
       END IF
       IF (core_or_shell) THEN
          WRITE (UNIT=output_unit,FMT="(A,3(1X,ES25.16),1X,I0)")&
            TRIM(ADJUSTL(kind_name)),r(1:3),particles%els(iparticle)%atom_index
       ELSE
          IF (LEN_TRIM(tag) > 0) THEN
             WRITE (UNIT=output_unit,FMT="(A,3(1X,ES25.16),1X,A,1X,I0)")&
               TRIM(ADJUSTL(kind_name)),r(1:3),TRIM(tag),imolecule
          ELSE
             WRITE (UNIT=output_unit,FMT="(A,3(1X,ES25.16))")&
               TRIM(ADJUSTL(kind_name)),r(1:3)
          END IF
       END IF
    END DO

    CALL timestop(handle)

  END SUBROUTINE dump_coordinates_cp2k

END MODULE input_restart_force_eval
