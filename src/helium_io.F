!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  I/O subroutines for helium
!> \author Lukasz Walewski
!> \date   2009-06-08
! *****************************************************************************
MODULE helium_io

  USE cell_types,                      ONLY: get_cell
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE helium_common,                   ONLY: helium_cycle_number,&
                                             helium_path_length,&
                                             helium_pbc,&
                                             helium_rotate
  USE helium_types,                    ONLY: e_id_interact,&
                                             e_id_kinetic,&
                                             e_id_potential,&
                                             e_id_thermo,&
                                             e_id_total,&
                                             e_id_virial,&
                                             helium_solvent_type
  USE input_constants,                 ONLY: fmt_id_pdb,&
                                             fmt_id_xyz,&
                                             helium_cell_shape_cube,&
                                             helium_cell_shape_octahedron
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_set,&
                                             section_vals_write
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE message_passing,                 ONLY: mp_gather
  USE parallel_rng_types,              ONLY: get_rng_stream
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_io'

  INTEGER, PUBLIC, PARAMETER :: wrt_mode_current = 1, &
                                wrt_mode_average = 2, &
                                wrt_mode_sampling = 3

  PUBLIC :: helium_update
  PUBLIC :: helium_write_restart
  PUBLIC :: helium_write_rdf
  PUBLIC :: helium_write_line
  PUBLIC :: helium_write_setup
  PUBLIC :: helium_write_energy
  PUBLIC :: helium_write_wnumber

  CONTAINS

  ! ***************************************************************************
  !> \brief  Update the helium state.
  !> \author Lukasz Walewski
  !> \date   2009-11-12
  !> \descr  Transfer the current helium state from the runtime environment
  !>         to the input structure, so that it can be used for I/O, etc.
  !> \note   This routine does message passing - use with care.
  ! ***************************************************************************
  SUBROUTINE helium_update( helium, error )
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_update', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: msgstr, stmp, time_unit
    INTEGER                                  :: handle, itmp
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: time_start, time_stop, &
                                                time_used
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)
    time_start = m_walltime()

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    ! choose the print key
    NULLIFY(print_key)
    IF ( helium%solute_present ) THEN
      print_key => section_vals_get_subs_vals( helium%input, &
                   "MOTION%PINT%PRINT%RESTART", error=error)
    ELSE
      print_key => section_vals_get_subs_vals( helium%input, &
                   "MOTION%PINT%HELIUM%PRINT%RESTART", error=error)
    END IF

    ! decide whether to update or not
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
        basis_section=print_key,error=error),cp_p_file) ) THEN

      IF ( .NOT. helium%solute_present ) THEN
        ! update iteration number
        itmp = logger%iter_info%iteration(2)
        CALL section_vals_val_set( &
          helium%input, &
          "MOTION%PINT%ITERATION", &
          i_val=itmp, &
          error=error )
        ! else - PINT will do that
      END IF

      ! save coordinates
      CALL helium_coord_save( helium, error )

      ! save permutation state
      CALL helium_perm_save( helium, error )

      ! save RNG state
      CALL helium_rng_state_save( helium, error )

      IF (helium%solute_present) THEN
        ! save forces on the solute
        CALL helium_force_save( helium, error )
      END IF

      time_stop = m_walltime()
      time_used = time_stop - time_start
      time_unit = "sec"
      IF ( time_used .GT. 60 ) THEN
        time_used = time_used / 60.0_dp
        time_unit = "min"
      END IF
      IF ( time_used .GT. 60 ) THEN
        time_used = time_used / 60.0_dp
        time_unit = "hours"
      END IF
      stmp = ""
      WRITE(stmp, '(F20.1)') time_used
      msgstr = "Runtime data updated in " // TRIM(ADJUSTL(stmp))
      msgstr = TRIM(ADJUSTL(msgstr)) // " " // TRIM(ADJUSTL(time_unit)) // "."
      CALL helium_write_line(msgstr,error)

    END IF

    CALL timestop(handle)

    RETURN
  END SUBROUTINE helium_update

  ! ***************************************************************************
  !> \brief  Save coordinates from the runtime environment to the input tree.
  !> \author Lukasz Walewski
  !> \date   2009-11-09
  ! ***************************************************************************
  SUBROUTINE helium_coord_save( helium, error )
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_coord_save', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: msglen, status
    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:), POINTER     :: msg, msg_gather
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    CALL helium_rotate(helium,-helium%relrot, error)

    ! allocate the buffer to be passed and fill it with local coords at each proc
    NULLIFY(msg)
    msglen = SIZE(helium%pos)
    ALLOCATE(msg(msglen),STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)
    msg(:) = PACK( helium%pos, .TRUE. )

    ! allocate the buffer for message passing
    NULLIFY(msg_gather)
    msglen = SIZE(helium%pos) * logger%para_env%num_pe
    ALLOCATE(msg_gather(msglen),STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)

    ! pass the message from all processors to logger%para_env%source
    msg_gather(:) = 0.0_dp
    CALL mp_gather(msg,msg_gather,&
         logger%para_env%source,logger%para_env%group)

    ! update coordinates in the global input structure
    CALL section_vals_val_set(helium%input, &
         "MOTION%PINT%HELIUM%COORD%_DEFAULT_KEYWORD_", &
         r_vals_ptr=msg_gather,error=error)

    ! NULLIFY, but do not DEALLOCATE! - a new pointer to this array is silently
    ! assigned in section_vals_val_set - this memory will be used later on!
    ! "The val becomes the owner of the array" - from section_vals_val_set docu
    NULLIFY(msg_gather)

    ! DEALLOCATE since this array is only used locally
    DEALLOCATE(msg,STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)

    CALL helium_rotate(helium, helium%relrot, error)

    RETURN
  END SUBROUTINE helium_coord_save

  ! ***************************************************************************
  !> \brief  Save forces from the runtime environment to the input tree.
  !> \author Lukasz Walewski
  !> \date   2009-11-10
  !> \descr  The average forces on the solute at the time of restart update
  !>         are equal on all processors, since we average them over
  !>         replicated helium environments.
  ! ***************************************************************************
  SUBROUTINE helium_force_save( helium, error )
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_force_save', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: err_str, stmp
    INTEGER                                  :: msglen, reqlen, status
    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:), POINTER     :: msg
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    ! check that the number of values match the current runtime
    reqlen = helium%solute_atoms * helium%solute_beads * 3
    msglen = SIZE(helium%force_avrg)
    err_str = "Invalid size of HELIUM%FORCE: received '"
    stmp = ""
    WRITE(stmp,*) msglen
    err_str = TRIM(ADJUSTL(err_str)) // &
              TRIM(ADJUSTL(stmp)) // "' but expected '"
    stmp = ""
    WRITE(stmp,*) reqlen
    err_str = TRIM(ADJUSTL(err_str)) // &
              TRIM(ADJUSTL(stmp)) // "'."
    CALL cp_assert(msglen==reqlen,cp_failure_level,&
         cp_assertion_failed, routineP, err_str)

    ! allocate the buffer to be saved and fill it with forces
    NULLIFY(msg)
    ALLOCATE(msg(msglen),STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)
    msg(:) = PACK( helium%force_avrg, .TRUE. )

    ! update forces in the global input structure
    CALL section_vals_val_set(helium%input, &
         "MOTION%PINT%HELIUM%FORCE%_DEFAULT_KEYWORD_", &
         r_vals_ptr=msg,error=error)

    ! NULLIFY, but do not DEALLOCATE! - a new pointer to this array is silently
    ! assigned in section_vals_val_set - this memeory will be used later on!
    ! "The val becomes the owner of the array" - from section_vals_val_set docu
    NULLIFY(msg)

    RETURN
  END SUBROUTINE helium_force_save

  ! ***************************************************************************
  !> \brief  Save the current permutation state in the input structure.
  !> \author Lukasz Walewski
  !> \date   2009-11-05
  !> \descr  Gather the permutation array from all processors and transfer
  !>         collected data to the input structure. The inverse permutation
  !>         array is not saved since it can be easily recalculated given the
  !>         permutation state.
  ! ***************************************************************************
  SUBROUTINE helium_perm_save( helium, error )
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_perm_save', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: msglen, status
    INTEGER, DIMENSION(:), POINTER           :: msg_gather
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    ! allocate the buffer for message passing
    NULLIFY(msg_gather)
    msglen = SIZE(helium%permutation) * logger%para_env%num_pe
    ALLOCATE(msg_gather(msglen),STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)

    ! pass the message from all processors to logger%para_env%source
    msg_gather(:) = 0
    CALL mp_gather(helium%permutation,msg_gather,&
         logger%para_env%source,logger%para_env%group)

    ! update permutation state in the global input structure
    CALL section_vals_val_set(helium%input, &
         "MOTION%PINT%HELIUM%PERM%_DEFAULT_KEYWORD_", &
         i_vals_ptr=msg_gather,error=error)

    ! NULLIFY, but do not DEALLOCATE! - a new pointer to this array is silently
    ! assigned in section_vals_val_set - this memory will be used later on!
    ! "The val becomes the owner of the array" - from section_vals_val_set docu
    NULLIFY(msg_gather)

    RETURN
  END SUBROUTINE helium_perm_save

  ! ***************************************************************************
  !> \brief  Save RNG state to the input structure.
  !> \author Lukasz Walewski
  !> \date   2009-11-04
  ! ***************************************************************************
  SUBROUTINE helium_rng_state_save( helium, error )
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_rng_state_save', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: msglen, offset, status
    LOGICAL                                  :: failure, lbf
    REAL(kind=dp)                            :: bf, bu
    REAL(kind=dp), DIMENSION(3, 2)           :: bg, cg, ig
    REAL(kind=dp), DIMENSION(40)             :: msg
    REAL(kind=dp), DIMENSION(:), POINTER     :: msg_gather
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    ! pack RNG state on each processor to the local array
    CALL get_rng_stream(helium%rng_stream_uniform,bg=bg,cg=cg,ig=ig,&
         buffer=bu,buffer_filled=lbf,error=error)
    offset = 0
    msg(offset+1:offset+6)   = PACK( bg, .TRUE. )
    msg(offset+7:offset+12)  = PACK( cg, .TRUE. )
    msg(offset+13:offset+18) = PACK( ig, .TRUE. )
    IF ( lbf ) THEN
      bf = 1.0_dp
    ELSE
      bf = -1.0_dp
    END IF
    msg(offset+19) = bf
    msg(offset+20) = bu
    CALL get_rng_stream(helium%rng_stream_gaussian,bg=bg,cg=cg,ig=ig,&
         buffer=bu,buffer_filled=lbf,error=error)
    offset = 20
    msg(offset+1:offset+6)   = PACK( bg, .TRUE. )
    msg(offset+7:offset+12)  = PACK( cg, .TRUE. )
    msg(offset+13:offset+18) = PACK( ig, .TRUE. )
    IF ( lbf ) THEN
      bf = 1.0_dp
    ELSE
      bf = -1.0_dp
    END IF
    msg(offset+19) = bf
    msg(offset+20) = bu

    ! Gather RNG state (in msg_gather vector) from all processors at
    ! logger%para_env%source
    NULLIFY(msg_gather)
    msglen = SIZE(msg)*logger%para_env%num_pe
    ALLOCATE(msg_gather(msglen),STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)
    msg_gather(:) = 0.0_dp
    CALL mp_gather(msg,msg_gather,logger%para_env%source,logger%para_env%group)

    ! update the RNG state in the global input structure
    CALL section_vals_val_set(helium%input, &
         "MOTION%PINT%HELIUM%RNG_STATE%_DEFAULT_KEYWORD_", &
         r_vals_ptr=msg_gather,error=error)

    ! NULLIFY, but do not DEALLOCATE! - a new pointer to this array is silently
    ! assigned in section_vals_val_set - this memeory will be used later on!
    ! "The val becomes the owner of the array" - from section_vals_val_set docu
    NULLIFY(msg_gather)

    RETURN
  END SUBROUTINE helium_rng_state_save

  ! ***************************************************************************
  !> \brief  Update and write out the helium restart data
  !> \author hforbert
  !> \par    History
  !>         2009-07-03 added FORCE section with av. helium forces [lwalewski]
  ! ***************************************************************************
  SUBROUTINE helium_write_restart(helium,error)

    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_restart', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(helium%input),cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    ! If we are simulating helium only - actually write the restart file,
    ! otherwise the helium data will be written to the main restart file,
    ! created higher up the call tree (i.e. in pint_write_restart)
    IF ( helium%solute_present ) RETURN

    NULLIFY(print_key)
    print_key => section_vals_get_subs_vals( helium%input, &
                 "MOTION%PINT%HELIUM%PRINT%RESTART", error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
      basis_section=print_key,error=error),cp_p_file) ) THEN
      unit_nr = cp_print_key_unit_nr(logger,print_key,&
        middle_name="helium",extension=".restart",&
        do_backup=.TRUE.,error=error)
      ! cp_print_key_unit_nr returns -1 on nodes other than logger%para_env%ionode
      IF (unit_nr>0) THEN
        CALL section_vals_write(helium%input,unit_nr,hide_root=.TRUE.,&
          error=error)
      END IF
      CALL cp_print_key_finished_output(unit_nr,logger,print_key,&
        error=error)
      CALL helium_write_line("Restart data written to file.",error)
    END IF

    RETURN
  END SUBROUTINE helium_write_restart

  ! ***************************************************************************
  !> \brief  Write helium parameters to the output unit.
  !> \author Lukasz Walewski
  !> \date   2009-06-03
  ! ***************************************************************************
  SUBROUTINE helium_write_setup( helium, error )

    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_setup', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: my_label, stmp, unit_str
    INTEGER                                  :: i, itmp, j, unit_nr
    INTEGER(kind=int_8)                      :: i8tmp
    LOGICAL                                  :: failure, first
    REAL(kind=dp)                            :: rtmp, v1, v2, v3
    REAL(kind=dp), DIMENSION(3)              :: my_abc
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    my_label = "HELIUM| "

    IF (logger%para_env%ionode) THEN
      unit_nr = cp_logger_get_default_unit_nr(logger)

      WRITE(unit_nr,*)
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Number of helium environments:     ", helium%num_proc

      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Number of solvent atoms:           ", helium%atoms
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Number of solvent beads:           ", helium%beads
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Total number of solvent particles: ", helium%atoms*helium%beads

      unit_str = "angstrom^-3"
      rtmp = cp_unit_from_cp2k(helium%density, &
        unit_str, error=error)
      WRITE(unit_nr,'(T2,A,F12.6)') TRIM(my_label)//" Density   ["// &
        TRIM(unit_str)//"]:", rtmp

      unit_str = "angstrom"
      rtmp = cp_unit_from_cp2k(helium%cell_size, &
        unit_str, error=error)
      WRITE(unit_nr,'(T2,A,F12.6)') TRIM(my_label)//" Cell size ["// &
        TRIM(unit_str)//"]:   ", rtmp

      IF ( helium%periodic ) THEN
        IF ( helium%cell_shape .EQ. helium_cell_shape_cube ) THEN
          CALL helium_write_line("PBC cell shape: CUBE.",error)
        ELSE IF ( helium%cell_shape .EQ. helium_cell_shape_octahedron ) THEN
          CALL helium_write_line("PBC cell shape: TRUNCATED OCTAHEDRON.",error)
        ELSE
          CALL helium_write_line("*** Warning: unknown cell shape.",error)
        END IF
      ELSE
        CALL helium_write_line("PBC turned off.",error)
      END IF

      ! first step gets incremented during first iteration
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " First MC step                      :", helium%first_step + 1
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Last MC step                       :", helium%last_step
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Total number of MC steps           :", helium%num_steps
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Number of outer MC trials per step :", helium%iter_rot
      WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
        " Number of inner MC trials per step :", helium%iter_norot
      i8tmp = helium%iter_rot
      i8tmp = i8tmp * helium%iter_norot
      stmp = ""
      WRITE(stmp, *) i8tmp
      WRITE(unit_nr,'(T2,A)') TRIM(my_label)//&
        " Total number of MC trials per step : " // TRIM(ADJUSTL(stmp))
      i8tmp = helium%num_steps
      i8tmp = i8tmp * helium%iter_rot
      i8tmp = i8tmp * helium%iter_norot
      stmp = ""
      WRITE(stmp, *) i8tmp
      WRITE(unit_nr,'(T2,A)') TRIM(my_label) //&
        " Total number of MC trials          : " // TRIM(ADJUSTL(stmp))

      IF (helium%solute_present) THEN
        WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
          " Number of solute atoms:            ", helium%solute_atoms
        WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
          " Number of solute beads:            ", helium%solute_beads
        WRITE(unit_nr,'(T2,A,1X,I0)') TRIM(my_label)//&
          " Total number of solute particles:  ", helium%solute_atoms*&
          helium%solute_beads
        DO i = 1, 3
          WRITE(unit_nr,'(T2,A,1X,I0,A,I0,A)', ADVANCE='NO') TRIM(my_label)//&
            " Solute atom type: ", i,&
            "(",helium%solute_number(i),")"
          first = .TRUE.
          DO j = 1, helium%solute_number(i)
            IF (first) THEN
              WRITE(unit_nr,'(A)', ADVANCE='NO')", indices: "
              first = .FALSE.
            ELSE
              WRITE(unit_nr,'(A)', ADVANCE='NO') ", "
            END IF
            WRITE(unit_nr,'(I0)', ADVANCE='NO') helium%solute_index(i,j)
          END DO
          WRITE(unit_nr,'(1X)')
        END DO
        CALL get_cell(helium%solute_cell, abc=my_abc)
        unit_str = "angstrom"
        v1 = cp_unit_from_cp2k(my_abc(1), unit_str, error=error)
        v2 = cp_unit_from_cp2k(my_abc(2), unit_str, error=error)
        v3 = cp_unit_from_cp2k(my_abc(3), unit_str, error=error)
        WRITE(unit_nr,'(T2,A,F12.6,1X,F12.6,1X,F12.6)') &
          TRIM(my_label)//" Solute cell size ["// &
          TRIM(unit_str)//"]:   ", v1, v2, v3
      ELSE
        WRITE(unit_nr,'(T2,A)') TRIM(my_label)//" Solute is NOT present"
      END IF
    END IF

    ! radial distribution function related settings
    rtmp = cp_unit_from_cp2k(helium%rdf_delr, "angstrom", error=error)
    WRITE(stmp, '(1X,F12.6)') rtmp
    CALL helium_write_line("RDF| delr [angstrom]: "//TRIM(stmp),error)
    rtmp = cp_unit_from_cp2k(helium%rdf_maxr, "angstrom", error=error)
    WRITE(stmp, '(1X,F12.6)') rtmp
    CALL helium_write_line("RDF| maxr [angstrom]: "//TRIM(stmp),error)
    itmp = helium%rdf_nbin
    WRITE(stmp, '(I6)') itmp
    CALL helium_write_line("RDF| nbin           : "//TRIM(stmp),error)

    CALL helium_write_line("",error)

    RETURN
  END SUBROUTINE helium_write_setup

  ! ***************************************************************************
  !> \brief  Writes out a line of text to the default output unit.
  !> \author Lukasz Walewski
  !> \date   2009-07-10
  ! ***************************************************************************
  SUBROUTINE helium_write_line(line,error)

    CHARACTER(len=*), INTENT(IN)             :: line
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_line', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: my_label
    INTEGER                                  :: unit_nr
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    my_label = "HELIUM|"

    IF (logger%para_env%ionode) THEN
      unit_nr = cp_logger_get_default_unit_nr(logger)
      WRITE(unit_nr,'(T2,A)') TRIM(my_label)//" "//TRIM(line)
    END IF

    RETURN
  END SUBROUTINE helium_write_line

  ! ***************************************************************************
  !> \brief  Writes out helium energies according to HELIUM%PRINT%ENER_INFO
  !> \author Lukasz Walewski
  !> \date   2009-06-08
  ! ***************************************************************************
  SUBROUTINE helium_write_energy( helium, error )

    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iteration, unit_nr
    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:), POINTER     :: my_energy
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(helium%input),cp_failure_level,routineP,error,failure)

    my_energy => helium%energy_avrg

    NULLIFY(print_key,logger)
    print_key => section_vals_get_subs_vals(helium%input, &
      "MOTION%PINT%HELIUM%PRINT%ENER_INFO", error=error)
    logger => cp_error_get_logger(error)
    iteration = logger%iter_info%iteration(2)
    IF ( BTEST(cp_print_key_should_output(iteration_info=logger%iter_info,&
      basis_section=print_key,error=error),cp_p_file) ) THEN
      unit_nr=cp_print_key_unit_nr(logger,print_key, &
        middle_name="helium-energy",extension=".dat",error=error )
      ! cp_print_key_unit_nr returns -1 on nodes other than logger%para_env%ionode
      IF (unit_nr>0) THEN
        WRITE (unit_nr,"(I8,1X,F8.1,6(1X,F20.9))") &
          iteration, &
          helium%n_accepted_permutations, &
          my_energy(e_id_potential), &
          my_energy(e_id_kinetic), &
          my_energy(e_id_thermo), &
          my_energy(e_id_virial), &
          my_energy(e_id_interact), &
          my_energy(e_id_total)
        CALL m_flush(unit_nr)
      END IF
      CALL cp_print_key_finished_output(unit_nr,logger,print_key,error=error)
    END IF

    RETURN
  END SUBROUTINE helium_write_energy

  ! ***************************************************************************
  !> \brief  Writes out helium winding number according to HELIUM%PRINT%WNUMBER
  !> \author Lukasz Walewski
  !> \date   2009-10-19
  ! ***************************************************************************
  SUBROUTINE helium_write_wnumber( helium, error )

    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_wnumber', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iteration, unit_nr
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: wsquared
    REAL(kind=dp), DIMENSION(3)              :: my_wnumber
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(helium%input),cp_failure_level,routineP,error,failure)

    my_wnumber(:) = helium%wnumber_avrg(:)
    wsquared = my_wnumber(1)**2 + my_wnumber(2)**2 + my_wnumber(3)**2

    NULLIFY(print_key,logger)
    print_key => section_vals_get_subs_vals(helium%input, &
      "MOTION%PINT%HELIUM%PRINT%WNUMBER", error=error)
    logger => cp_error_get_logger(error)
    iteration = logger%iter_info%iteration(2)
    IF ( BTEST(cp_print_key_should_output(iteration_info=logger%iter_info,&
      basis_section=print_key,error=error),cp_p_file) ) THEN
      unit_nr=cp_print_key_unit_nr(logger,print_key, &
        middle_name="helium-wnumber",extension=".dat",error=error )
      IF (unit_nr>0) THEN
        WRITE (unit_nr,"(I8,4(1X,F20.9))") &
          iteration, wsquared, my_wnumber(1), my_wnumber(2), my_wnumber(3)
        CALL m_flush(unit_nr)
      END IF
      CALL cp_print_key_finished_output(unit_nr,logger,print_key,error=error)
    END IF

    RETURN
  END SUBROUTINE helium_write_wnumber

  ! ***************************************************************************
  !> \brief  Writes helium configuration according to HELIUM%PRINT%COORDINATES
  !> \author Lukasz Walewski
  !> \date   2009-07-16
  !> \note   particle_types->write_particle_coordinates is of no use here,
  !>         since it does not support atom connectivity information needed
  !>         for helium paths
  ! ***************************************************************************
  SUBROUTINE helium_write_coordinates( helium, r0, fmt_id, middle_name, error )

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(kind=dp), DIMENSION(3), INTENT(IN)  :: r0
    INTEGER, INTENT(IN)                      :: fmt_id
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: middle_name
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_coordinates', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: fmt_string, my_middle_name, &
                                                stmp
    INTEGER                                  :: ia, ib, ib1, ib2, ic, icycle, &
                                                tmp1, tmp2, unit_nr
    LOGICAL                                  :: are_connected, failure, ltmp, &
                                                should_output
    REAL(kind=dp)                            :: xtmp, ytmp, ztmp
    REAL(kind=dp), DIMENSION(3)              :: r1, r2
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(helium%input),cp_failure_level,routineP,error,failure)

    IF (PRESENT(middle_name)) THEN
      my_middle_name = middle_name
    ELSE
      my_middle_name = "helium-pos"
    END IF

    NULLIFY(logger,print_key)
    logger => cp_error_get_logger(error)

    print_key => section_vals_get_subs_vals(helium%input, &
      "MOTION%PINT%HELIUM%PRINT%COORDINATES", error=error)
    should_output = BTEST(cp_print_key_should_output( &
      iteration_info=logger%iter_info,&
      basis_section=print_key,error=error),cp_p_file)
    IF (should_output) THEN

      ! prepare the coordinates for output (print one image of the periodic
      ! space that is centered around r0)
      DO ia = 1, helium%atoms
        DO ib = 1, helium%beads
          r1(:) = helium%pos(:,ia,ib) - r0(:)
          r2(:) = helium%pos(:,ia,ib) - r0(:)
          CALL helium_pbc( helium, r2 )
          ltmp = .FALSE.
          DO ic = 1, 3
            IF ( ABS(r1(ic)-r2(ic)) .GT. 100.0_dp*EPSILON(0.0_dp) ) THEN
              ltmp = .TRUE.
              CYCLE
                END IF
          END DO
          IF ( ltmp ) THEN
            helium%work(:,ia,ib) = r0(:) + r2(:)
          ELSE
            helium%work(:,ia,ib) = helium%pos(:,ia,ib)
          END IF
        END DO
      END DO

!TODO write out each He replica (at the moment only coords from rank = 0 are printed)
      IF (logger%para_env%ionode) THEN

        SELECT CASE (fmt_id)

          CASE (fmt_id_xyz)
            unit_nr=cp_print_key_unit_nr(logger,print_key, &
              middle_name=my_middle_name,extension=".xyz",error=error)
            fmt_string = "('He',3(1X,F20.10))"
            WRITE(unit_nr,'(I0)') helium%atoms*helium%beads
            WRITE(unit_nr,'(A)') "He coordinates"
            DO ia = 1, helium%atoms
              DO ib = 1, helium%beads
                xtmp = helium%work(1,ia,ib)
                xtmp = cp_unit_from_cp2k(xtmp, "angstrom", error=error)
                ytmp = helium%work(2,ia,ib)
                ytmp = cp_unit_from_cp2k(ytmp, "angstrom", error=error)
                ztmp = helium%work(3,ia,ib)
                ztmp = cp_unit_from_cp2k(ztmp, "angstrom", error=error)
                WRITE(unit_nr,fmt_string) xtmp, ytmp, ztmp
              END DO
            END DO
            CALL m_flush(unit_nr)
            CALL cp_print_key_finished_output(unit_nr,logger,print_key,&
              error=error)

          CASE (fmt_id_pdb)
            unit_nr=cp_print_key_unit_nr(logger,print_key, &
              middle_name=my_middle_name,extension=".pdb",error=error)
            ! write out the unit cell parameters
            fmt_string = "(A6,3F9.3,3F7.2,1X,A11,1X,I3)"
            xtmp = helium%cell_size
            xtmp = cp_unit_from_cp2k(xtmp, "angstrom", error=error)
            IF ( helium%cell_shape .EQ. helium_cell_shape_octahedron ) THEN
              stmp = "O          "
            ELSE
              stmp = "C          "
            END IF
            WRITE(unit_nr,fmt_string) "CRYST1", &
              xtmp, xtmp, xtmp, &
              90.0_dp, 90.0_dp, 90.0_dp, &
              stmp, helium%beads
            ! write out the bead coordinates
            fmt_string = &
              "(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)"
            DO ia = 1, helium%atoms
              icycle = helium_cycle_number(helium,ia)
              DO ib = 1, helium%beads
                xtmp = helium%work(1,ia,ib)
                xtmp = cp_unit_from_cp2k(xtmp, "angstrom", error=error)
                ytmp = helium%work(2,ia,ib)
                ytmp = cp_unit_from_cp2k(ytmp, "angstrom", error=error)
                ztmp = helium%work(3,ia,ib)
                ztmp = cp_unit_from_cp2k(ztmp, "angstrom", error=error)
                WRITE(unit_nr,fmt_string) "ATOM  ", &
                  (ia-1)*helium%beads+ib, &
                  " He ", " ", "   ", "X", &
                  icycle, &
                  " ", &
                  xtmp, ytmp, ztmp, &
                  1.0_dp, 0.0_dp, "HE", "  "
              END DO
            END DO

            ! write out the bead connectivity information
            DO ia = 1, helium%atoms

              ! write connectivity records for this atom only if the path
              ! it belongs to is longer than 1.
              IF ( helium_path_length(helium,ia) .LE. 1 ) THEN
                CYCLE
              END IF

              DO ib = 1, helium%beads-1
                ! check wheather the consecutive beads belong to the same box
                r1(:) = helium%work(:,ia,ib) - helium%work(:,ia,ib+1)
                r2(:) = r1(:)
                CALL helium_pbc( helium, r2 )
                are_connected = .TRUE.
                DO ic = 1, 3
                  IF ( ABS(r1(ic)-r2(ic)) .GT. 100.0_dp*EPSILON(0.0_dp) ) THEN
                        ! if the distance betw ib and ib+1 changes upon applying
                        ! PBC do not connect them
                        are_connected = .FALSE.
                        CYCLE
                      END IF
                    END DO
                IF ( are_connected ) THEN
                  tmp1 = (ia-1)*helium%beads+ib
                  tmp2 = (ia-1)*helium%beads+ib+1
                  ! smaller value has to go first
                  IF (tmp1 .LT. tmp2) THEN
                    ib1 = tmp1
                    ib2 = tmp2
                  ELSE
                    ib1 = tmp2
                    ib2 = tmp1
                  END IF
                  WRITE(unit_nr,'(A6,2I5)') "CONECT", ib1, ib2
                END IF
              END DO
              ! last bead of atom <ia> connects to the first bead
              ! of the next atom in the permutation cycle
              r1(:) = helium%work(:,ia,helium%beads) - helium%work(:,helium%permutation(ia),1)
              r2(:) = r1(:)
              CALL helium_pbc( helium, r2 )
              are_connected = .TRUE.
              DO ic = 1, 3
                IF ( ABS(r1(ic)-r2(ic)) .GT. 100.0_dp*EPSILON(0.0_dp) ) THEN
                      ! if the distance betw ib and ib+1 changes upon applying
                      ! PBC do not connect them
                      are_connected = .FALSE.
                      CYCLE
                    END IF
                  END DO
              IF ( are_connected ) THEN
                tmp1 = ia*helium%beads
                tmp2 = (helium%permutation(ia)-1)*helium%beads+1
                IF (tmp1 .LT. tmp2) THEN
                  ib1 = tmp1
                  ib2 = tmp2
                ELSE
                  ib1 = tmp2
                  ib2 = tmp1
                END IF
                WRITE(unit_nr,'(A6,2I5)') "CONECT", ib1, ib2
              END IF
            END DO
            WRITE(unit_nr,'(A)') "END"
            CALL m_flush(unit_nr)
            CALL cp_print_key_finished_output(unit_nr,logger,print_key,&
              error=error)

          CASE DEFAULT
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
            routineP,"Unknown file format id ("//&
            TRIM(ADJUSTL(cp_to_string(fmt_id)))//")",&
            error,failure)
        END SELECT

      END IF
    END IF

    RETURN
  END SUBROUTINE helium_write_coordinates

  ! ***************************************************************************
  !> \brief  Write helium RDF according to HELIUM%PRINT%RDF
  !> \author Lukasz Walewski
  !> \date   2009-07-23
  ! ***************************************************************************
  SUBROUTINE helium_write_rdf( helium, error )

    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_write_rdf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, unit_nr
    LOGICAL                                  :: failure, is_new, should_output
    REAL(kind=dp)                            :: rtmp
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(helium%input),cp_failure_level,routineP,error,failure)

    NULLIFY(logger,print_key)
    logger => cp_error_get_logger(error)
    print_key => section_vals_get_subs_vals(helium%input, &
      "MOTION%PINT%HELIUM%PRINT%RDF", error=error)
    should_output = BTEST(cp_print_key_should_output( &
      iteration_info=logger%iter_info,&
      basis_section=print_key,error=error),cp_p_file)
    IF (should_output) THEN
      IF (logger%para_env%ionode) THEN
        unit_nr=cp_print_key_unit_nr(logger,print_key, &
          middle_name="helium-hst",extension=".dat",&
          is_new_file=is_new, error=error)
        IF (.NOT. is_new) THEN
          WRITE(unit_nr,'(A1)') "&"
        END IF
        DO i = 1, helium%rdf_nbin
          rtmp = ( REAL(i) - 0.5_dp ) * helium%rdf_delr
          rtmp = cp_unit_from_cp2k(rtmp, "angstrom", error=error)
          WRITE(unit_nr,'(2F20.10)') rtmp,&
            helium%rdf_avrg(i)
        END DO
        CALL m_flush(unit_nr)
        CALL cp_print_key_finished_output(unit_nr,logger,print_key,&
          error=error)
      END IF
    END IF

    RETURN
  END SUBROUTINE helium_write_rdf

END MODULE helium_io
