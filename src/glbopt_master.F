!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Master's routines for global optimization
!> \author Ole Schuett
! *****************************************************************************
MODULE glbopt_master
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             deallocate_atomic_kind_set
  USE colvar_types,                    ONLY: colvar_p_release,&
                                             colvar_p_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE exclusion_types,                 ONLY: exclusion_release,&
                                             exclusion_type
  USE glbopt_mincrawl,                 ONLY: mincrawl_finalize,&
                                             mincrawl_init,&
                                             mincrawl_steer,&
                                             mincrawl_type
  USE glbopt_minhop,                   ONLY: minhop_finalize,&
                                             minhop_init,&
                                             minhop_steer,&
                                             minhop_type
  USE input_constants,                 ONLY: dump_xmol,&
                                             glbopt_do_mincrawl,&
                                             glbopt_do_minhop
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8
  USE molecule_kind_types,             ONLY: deallocate_molecule_kind_set,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: deallocate_global_constraint,&
                                             deallocate_molecule_set,&
                                             global_constraint_type,&
                                             molecule_type
  USE particle_types,                  ONLY: deallocate_particle_set,&
                                             particle_type,&
                                             write_particle_coordinates
  USE swarm_message,                   ONLY: swarm_message_get,&
                                             swarm_message_type
  USE topology,                        ONLY: topology_control
#include "./common/cp_common_uses.f90"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_master'


 PUBLIC :: glbopt_master_type
 PUBLIC :: glbopt_master_init, glbopt_master_finalize
 PUBLIC :: glbopt_master_steer


 TYPE glbopt_master_type
   PRIVATE
   REAL(KIND=dp)                                       :: E_lowest = HUGE(1.0_dp)
   REAL(KIND=dp)                                       :: E_target = TINY(1.0_dp)
   INTEGER                                             :: iw = 0
   INTEGER                                             :: progress_traj_unit =0
   INTEGER(int_8)                                      :: total_md_steps = 0
   INTEGER(int_8)                                      :: total_gopt_steps = 0
   INTEGER(int_8)                                      :: count_reports = 0
   INTEGER                                             :: method
   TYPE(minhop_type),POINTER                           :: minhop
   TYPE(mincrawl_type), POINTER                        :: mincrawl
   INTEGER                                             :: i_iteration = 0
   TYPE(atomic_kind_type), DIMENSION(:), POINTER       :: atomic_kind_set => Null()
   TYPE(particle_type), DIMENSION(:), POINTER          :: particle_set => Null()
   TYPE(section_vals_type), POINTER                    :: glbopt_section => Null()
   TYPE(cp_error_type)                                 :: error
 END TYPE glbopt_master_type

 CONTAINS


! *****************************************************************************
!> \brief Initializes the master of the global optimizer
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE glbopt_master_init(this, para_env, root_section, input_path, n_walkers, iw, error)
    TYPE(glbopt_master_type), INTENT(INOUT)  :: this
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: root_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(IN)                      :: n_walkers, iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)

    this%iw = iw
    this%error = error

    this%glbopt_section => section_vals_get_subs_vals(root_section, "SWARM%GLOBAL_OPT", error=error)
    CALL section_vals_retain(this%glbopt_section,error=error)

    CALL section_vals_val_get(this%glbopt_section,"E_TARGET", r_val=this%E_target, error=error)
    CALL section_vals_val_get(this%glbopt_section,"METHOD", i_val=this%method, error=error)

    CALL glbopt_read_particle_set(this, para_env, root_section, iw, error)

    logger => cp_error_get_logger(error)
    this%progress_traj_unit = cp_print_key_unit_nr(logger, &
      this%glbopt_section,"PROGRESS_TRAJECTORY", &
      middle_name="progress", extension=".xyz", &
      file_action="WRITE", file_position="REWIND", error=error)

    SELECT CASE (this%method)
       CASE(glbopt_do_minhop)
         ALLOCATE(this%minhop)
         CALL minhop_init(this%minhop, this%glbopt_section, n_walkers, iw, error)
       CASE(glbopt_do_mincrawl)
         ALLOCATE(this%mincrawl)
         CALL mincrawl_init(this%mincrawl, this%glbopt_section, n_walkers, iw, this%particle_set, error)
       CASE DEFAULT
        STOP "Unkown glbopt_method"
    END SELECT
 END SUBROUTINE glbopt_master_init


! *****************************************************************************
!> \brief Helper-routine for glbopt_master_init, reads part of SUBSYS-section
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE glbopt_read_particle_set(this, para_env, root_section, iw, error)
    TYPE(glbopt_master_type), INTENT(INOUT)  :: this
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: root_section
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(colvar_p_type), DIMENSION(:), &
      POINTER                                :: colvar_p
    TYPE(exclusion_type), DIMENSION(:), &
      POINTER                                :: exclusions
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section

    NULLIFY(atomic_kind_set, particle_set, molecule_kind_set, molecule_set)
    NULLIFY(colvar_p, gci, exclusions, force_env_section, subsys_section)

    force_env_section => section_vals_get_subs_vals(root_section, "FORCE_EVAL", error=error)
    subsys_section => section_vals_get_subs_vals(root_section, "FORCE_EVAL%SUBSYS", error=error)


    CALL topology_control(atomic_kind_set, &
                          particle_set, &
                          molecule_kind_set, &
                          molecule_set, &
                          colvar_p, &
                          gci, &
                          root_section=root_section, &
                          para_env=para_env, &
                          force_env_section=force_env_section, &
                          subsys_section=subsys_section, &
                          use_motion_section=.FALSE., &
                          exclusions=exclusions, &
                          error=error)

    !This is the only thing we need to write trajectories.
    this%atomic_kind_set => atomic_kind_set
    this%particle_set => particle_set
    CALL exclusion_release(exclusions, error)
    CALL deallocate_molecule_set(molecule_set, error)
    CALL deallocate_molecule_kind_set(molecule_kind_set, error)
    CALL deallocate_global_constraint(gci, error)
    CALL colvar_p_release(colvar_p, error)

  END SUBROUTINE glbopt_read_particle_set


! *****************************************************************************
!> \brief Central steering routine of global optimizer
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE glbopt_master_steer(this, report, cmd, should_stop)
    TYPE(glbopt_master_type), INTENT(INOUT)  :: this
    TYPE(swarm_message_type)                 :: report, cmd
    LOGICAL, INTENT(INOUT)                   :: should_stop

    CALL progress_report(this, report)

    SELECT CASE (this%method)
      CASE(glbopt_do_minhop)
          CALL minhop_steer(this%minhop, report, cmd)
      CASE(glbopt_do_mincrawl)
          CALL mincrawl_steer(this%mincrawl, report, cmd)
      CASE DEFAULT
       STOP "Unkown glbopt_method"
    END SELECT

    IF(this%E_lowest<this%E_target) THEN
       IF(this%iw>0)  WRITE (this%iw,"(A,I8,A)" )  &
          " GLBOPT| Reached E_pot < E_target after ", this%i_iteration," iterations. Quitting."
       should_stop = .TRUE.
    END IF
 END SUBROUTINE glbopt_master_steer


! *****************************************************************************
!> \brief Helper routine for glbopt_master_steer(), updates stats, etc.
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE progress_report(this, report)
    TYPE(glbopt_master_type), INTENT(INOUT)  :: this
    TYPE(swarm_message_type)                 :: report

    CHARACTER(len=default_string_length)     :: status
    INTEGER                                  :: gopt_steps, md_steps, &
                                                report_worker_id
    REAL(KIND=dp)                            :: report_Epot

    this%i_iteration = this%i_iteration + 1

    CALL swarm_message_get(report, "worker_id", report_worker_id)
    CALL swarm_message_get(report, "status", status)

    IF(TRIM(status) == "ok") THEN
       CALL swarm_message_get(report, "Epot", report_Epot)
       CALL swarm_message_get(report, "md_steps", md_steps)
       CALL swarm_message_get(report, "gopt_steps", gopt_steps)
       this%total_md_steps   = this%total_md_steps + md_steps
       this%total_gopt_steps = this%total_gopt_steps + gopt_steps
       this%count_reports    = this%count_reports + 1

       IF(report_Epot < this%E_lowest) THEN
          this%E_lowest = report_Epot
          CALL write_progress_traj(this, report)
       ENDIF

       IF(this%iw>0) THEN
          WRITE (this%iw,'(A,46X,I8)')&
             " GLBOPT| Reporting worker ", report_worker_id
          WRITE (this%iw,'(A,20X,E15.8)')&
             " GLBOPT| Reported potential energy [Hartree] ", report_Epot
          WRITE (this%iw,'(A,13X,E15.8)')&
             " GLBOPT| Lowest reported potential energy [Hartree] ", this%E_lowest
          WRITE (this%iw,'(A,T71,F10.1)')&
             " GLBOPT| Average number of MD steps", REAL(this%total_md_steps,KIND=dp)/this%count_reports
          WRITE (this%iw,'(A,T71,F10.1)')&
             " GLBOPT| Average number of Geo-Opt steps", REAL(this%total_gopt_steps,KIND=dp)/this%count_reports
       END IF
    END IF
 END SUBROUTINE progress_report


! *****************************************************************************
!> \brief Helper routine for progress_report(), write frames to trajectory.
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE write_progress_traj(this, report)
    TYPE(glbopt_master_type), INTENT(INOUT)  :: this
    TYPE(swarm_message_type), INTENT(IN)     :: report

    CHARACTER(len=default_string_length)     :: title, unit_str
    INTEGER                                  :: report_worker_id
    REAL(KIND=dp)                            :: report_Epot, unit_conv
    REAL(KIND=dp), DIMENSION(:), POINTER     :: report_positions

    NULLIFY(report_positions)

    IF(this%progress_traj_unit <= 0) RETURN

    CALL swarm_message_get(report, "worker_id", report_worker_id)
    CALL swarm_message_get(report, "positions", report_positions)
    CALL swarm_message_get(report, "Epot", report_Epot)

    WRITE (title, '(A,I8,A,I5,A,F20.10)') 'i = ',this%i_iteration, &
      ' worker_id = ',report_worker_id, ' Epot = ', report_Epot

    !get the conversion factor for the length unit
    CALL section_vals_val_get(this%glbopt_section,"PROGRESS_TRAJECTORY%UNIT",&
            c_val=unit_str,error=this%error)
    unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=this%error)
    CALL write_particle_coordinates(this%particle_set,&
                                    iunit=this%progress_traj_unit, &
                                    output_format=dump_xmol, &
                                    content="POS", &
                                    title=TRIM(title), &
                                    array=report_positions,&
                                    unit_conv=unit_conv,&
                                    error=this%error)
     DEALLOCATE(report_positions)
 END SUBROUTINE write_progress_traj


! *****************************************************************************
!> \brief Finalized the master of the global optimizer
!> \author Ole
! *****************************************************************************
 SUBROUTINE glbopt_master_finalize(this)
    TYPE(glbopt_master_type), INTENT(INOUT)  :: this

    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)

    SELECT CASE (this%method)
       CASE(glbopt_do_minhop)
         CALL minhop_finalize(this%minhop)
         DEALLOCATE(this%minhop)
       CASE(glbopt_do_mincrawl)
         CALL mincrawl_finalize(this%mincrawl)
         DEALLOCATE(this%mincrawl)
       CASE DEFAULT
        STOP "Unkown glbopt_method"
    END SELECT

    logger => cp_error_get_logger(this%error)
    CALL cp_print_key_finished_output(this%progress_traj_unit, logger, &
      this%glbopt_section, "PROGRESS_TRAJECTORY", error=this%error)

    CALL section_vals_release(this%glbopt_section, error=this%error)
    CALL deallocate_particle_set(this%particle_set, this%error)
    CALL deallocate_atomic_kind_set(this%atomic_kind_set, this%error)

 END SUBROUTINE glbopt_master_finalize

END MODULE glbopt_master

