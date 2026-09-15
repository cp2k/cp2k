!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_event_analysis
   USE trajana_cell_source,             ONLY: cell_source_type
   USE trajana_command_line,            ONLY: fail,&
                                              get_integer_option,&
                                              get_option,&
                                              get_real_option
   USE trajana_frame_controls,          ONLY: frame_selected
   USE trajana_geometry_analysis,       ONLY: action_type,&
                                              evaluate_action,&
                                              read_actions,&
                                              validate_actions
   USE trajana_kinds,                   ONLY: dp
   USE trajana_text_utils,              ONLY: lower_case
   USE trajana_trajectory_io,           ONLY: close_output,&
                                              open_output,&
                                              xyz_reader_type
   USE trajana_trajectory_types,        ONLY: frame_type

   IMPLICIT NONE
   PRIVATE

   TYPE :: event_state_type
      LOGICAL :: active = .FALSE.
      LOGICAL :: seen = .FALSE.
      INTEGER :: active_samples = 0
      INTEGER :: event_count = 0
      INTEGER :: open_events = 0
      INTEGER :: start_frame = 0
      INTEGER :: last_frame = 0
      INTEGER :: extremum_frame = 0
      REAL(dp) :: start_ps = 0.0_dp
      REAL(dp) :: last_ps = 0.0_dp
      REAL(dp) :: extremum = 0.0_dp
      REAL(dp) :: extremum_ps = 0.0_dp
      REAL(dp) :: first_entry_ps = -1.0_dp
      REAL(dp) :: total_residence_ps = 0.0_dp
      REAL(dp) :: longest_residence_ps = 0.0_dp
      REAL(dp) :: trajectory_extremum = 0.0_dp
   END TYPE event_state_type

   PUBLIC :: run_events

CONTAINS

   SUBROUTINE run_events()
      CHARACTER(LEN=512)                                 :: message
      CHARACTER(LEN=:), ALLOCATABLE                      :: action_path, cell_path, cell_text, &
                                                            direction, input_path, output_path, &
                                                            periodic_text, summary_path
      INTEGER                                            :: expected_atoms, first, frame_count, &
                                                            ierr, item, last, output_unit, stride, &
                                                            summary_unit
      LOGICAL                                            :: above, eof, found, summary_requested
      REAL(dp)                                           :: dt_fs, effective_dt_ps, entry_value, &
                                                            exit_value, time_ps
      REAL(dp), ALLOCATABLE                              :: result(:)
      TYPE(action_type), ALLOCATABLE                     :: actions(:)
      TYPE(cell_source_type)                             :: cells
      TYPE(event_state_type), ALLOCATABLE                :: states(:)
      TYPE(frame_type)                                   :: frame
      TYPE(xyz_reader_type)                              :: trajectory

      CALL get_option("--input", input_path, found, "-")
      CALL get_option("--output", output_path, found, "-")
      CALL get_option("--summary", summary_path, summary_requested)
      CALL get_option("--actions", action_path, found, "trajana.in")
      CALL get_option("--cell", cell_text, found, "")
      CALL get_option("--cell-file", cell_path, found, "")
      CALL get_option("--periodic", periodic_text, found, "XYZ")
      CALL get_option("--direction", direction, found, "above")
      CALL get_real_option("--dt-fs", dt_fs, -1.0_dp)
      CALL get_real_option("--entry", entry_value, HUGE(1.0_dp))
      CALL get_real_option("--exit", exit_value, HUGE(1.0_dp))
      CALL get_integer_option("--first", first, 1)
      CALL get_integer_option("--last", last, HUGE(1))
      CALL get_integer_option("--stride", stride, 1)

      IF (dt_fs <= 0.0_dp) CALL fail("events requires a positive --dt-fs")
      IF (entry_value > 0.5_dp*HUGE(1.0_dp)) CALL fail("events requires --entry")
      IF (exit_value > 0.5_dp*HUGE(1.0_dp)) CALL fail("events requires --exit")
      IF (first < 1 .OR. last < first .OR. stride < 1) CALL fail("Invalid frame range")
      direction = lower_case(direction)
      SELECT CASE (direction)
      CASE ("above")
         above = .TRUE.
         IF (exit_value >= entry_value) CALL fail("For --direction above, --exit must be smaller than --entry")
      CASE ("below")
         above = .FALSE.
         IF (exit_value <= entry_value) CALL fail("For --direction below, --exit must be larger than --entry")
      CASE DEFAULT
         CALL fail("--direction expects above or below")
      END SELECT
      effective_dt_ps = dt_fs*REAL(stride, dp)/1000.0_dp

      CALL read_actions(action_path, actions, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL cells%configure(cell_text, cell_path, periodic_text, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL trajectory%open_file(input_path, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL open_output(output_path, output_unit, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      WRITE (output_unit, "(A)") "# Hysteretic events in geometric trajectory observables"
      WRITE (output_unit, "(A,A)") "# direction: ", TRIM(direction)
      WRITE (output_unit, "(A,ES24.16)") "# entry: ", entry_value
      WRITE (output_unit, "(A,ES24.16)") "# exit: ", exit_value
      WRITE (output_unit, "(A,ES24.16)") "# selected_frame_time_step_ps: ", effective_dt_ps
      WRITE (output_unit, "(A)") &
         "# action event start_frame end_frame start_ps end_ps residence_ps extremum "// &
         "extremum_frame extremum_ps returned_past_exit"

      ALLOCATE (RESULT(SIZE(actions)), states(SIZE(actions)))
      expected_atoms = -1
      frame_count = 0
      DO
         CALL trajectory%read_frame(frame, eof, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         IF (eof) EXIT
         CALL cells%attach(frame, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         IF (frame%number > last) EXIT
         IF (.NOT. frame_selected(frame%number, first, last, stride)) CYCLE
         IF (expected_atoms < 0) THEN
            expected_atoms = frame%natoms
            CALL validate_actions(actions, expected_atoms, ierr, message)
            IF (ierr /= 0) CALL fail(message)
         END IF
         IF (frame%natoms /= expected_atoms) CALL fail("The atom count changes along the trajectory")

         frame_count = frame_count + 1
         time_ps = REAL(frame_count - 1, dp)*effective_dt_ps
         DO item = 1, SIZE(actions)
            CALL evaluate_action(actions(item), frame, RESULT(item), ierr)
            IF (ierr /= 0) CALL fail("Degenerate geometry or singular cell in frame")
            CALL update_state(states(item), RESULT(item), frame%number, time_ps, effective_dt_ps, &
                              entry_value, exit_value, above, item, output_unit)
         END DO
      END DO
      IF (frame_count == 0) CALL fail("No trajectory frames were selected")
      DO item = 1, SIZE(actions)
         IF (states(item)%active) THEN
            states(item)%open_events = states(item)%open_events + 1
            CALL write_event(states(item), item, effective_dt_ps, .FALSE., output_unit)
            states(item)%active = .FALSE.
         END IF
      END DO

      CALL close_output(output_unit)
      CALL trajectory%close_file()
      CALL cells%close()
      IF (summary_requested) THEN
         CALL open_output(summary_path, summary_unit, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         CALL write_summary(summary_unit, actions, states, direction, entry_value, exit_value, frame_count, &
                            effective_dt_ps)
         CALL close_output(summary_unit)
      END IF
   END SUBROUTINE run_events

   SUBROUTINE update_state(state, value, frame, time_ps, dt_ps, entry_value, exit_value, &
                           above, action, output_unit)
      TYPE(event_state_type), INTENT(INOUT)              :: state
      REAL(dp), INTENT(IN)                               :: value
      INTEGER, INTENT(IN)                                :: frame
      REAL(dp), INTENT(IN)                               :: time_ps, dt_ps, entry_value, exit_value
      LOGICAL, INTENT(IN)                                :: above
      INTEGER, INTENT(IN)                                :: action, output_unit

      LOGICAL                                            :: entered, exited

      IF (.NOT. state%seen) THEN
         state%trajectory_extremum = value
         state%seen = .TRUE.
      ELSE IF ((above .AND. value > state%trajectory_extremum) .OR. &
               (.NOT. above .AND. value < state%trajectory_extremum)) THEN
         state%trajectory_extremum = value
      END IF
      IF (above) THEN
         entered = value >= entry_value
         exited = value <= exit_value
      ELSE
         entered = value <= entry_value
         exited = value >= exit_value
      END IF

      IF (.NOT. state%active) THEN
         IF (entered) CALL start_event(state, value, frame, time_ps)
      ELSE IF (exited) THEN
         CALL write_event(state, action, dt_ps, .TRUE., output_unit)
         state%active = .FALSE.
      ELSE
         state%active_samples = state%active_samples + 1
         state%last_frame = frame
         state%last_ps = time_ps
         IF ((above .AND. value > state%extremum) .OR. (.NOT. above .AND. value < state%extremum)) THEN
            state%extremum = value
            state%extremum_frame = frame
            state%extremum_ps = time_ps
         END IF
      END IF
   END SUBROUTINE update_state

   SUBROUTINE start_event(state, value, frame, time_ps)
      TYPE(event_state_type), INTENT(INOUT)              :: state
      REAL(dp), INTENT(IN)                               :: value
      INTEGER, INTENT(IN)                                :: frame
      REAL(dp), INTENT(IN)                               :: time_ps

      state%active = .TRUE.
      state%active_samples = 1
      state%event_count = state%event_count + 1
      state%start_frame = frame
      state%last_frame = frame
      state%extremum_frame = frame
      state%start_ps = time_ps
      state%last_ps = time_ps
      state%extremum = value
      state%extremum_ps = time_ps
      IF (state%first_entry_ps < 0.0_dp) state%first_entry_ps = time_ps
   END SUBROUTINE start_event

   SUBROUTINE write_event(state, action, dt_ps, returned, unit)
      TYPE(event_state_type), INTENT(INOUT)              :: state
      INTEGER, INTENT(IN)                                :: action
      REAL(dp), INTENT(IN)                               :: dt_ps
      LOGICAL, INTENT(IN)                                :: returned
      INTEGER, INTENT(IN)                                :: unit

      REAL(dp)                                           :: residence_ps

      residence_ps = REAL(state%active_samples, dp)*dt_ps
      state%total_residence_ps = state%total_residence_ps + residence_ps
      state%longest_residence_ps = MAX(state%longest_residence_ps, residence_ps)
      WRITE (unit, "(2I8,2I12,4ES24.16,I12,ES24.16,I4)") action, state%event_count, &
         state%start_frame, state%last_frame, state%start_ps, state%last_ps, residence_ps, state%extremum, &
         state%extremum_frame, state%extremum_ps, MERGE(1, 0, returned)
   END SUBROUTINE write_event

   SUBROUTINE write_summary(unit, actions, states, direction, entry_value, exit_value, frames, dt_ps)
      INTEGER, INTENT(IN)                                :: unit
      TYPE(action_type), INTENT(IN)                      :: actions(:)
      TYPE(event_state_type), INTENT(IN)                 :: states(:)
      CHARACTER(LEN=*), INTENT(IN)                       :: direction
      REAL(dp), INTENT(IN)                               :: entry_value, exit_value
      INTEGER, INTENT(IN)                                :: frames
      REAL(dp), INTENT(IN)                               :: dt_ps

      INTEGER                                            :: item

      WRITE (unit, "(A)") "# Summary of hysteretic events"
      WRITE (unit, "(A,I0)") "# selected_frames: ", frames
      WRITE (unit, "(A,ES24.16)") "# selected_frame_time_step_ps: ", dt_ps
      WRITE (unit, "(A)") &
         "# action kind atom1 atom2 atom3 atom4 direction entry exit events open_events first_entry_ps "// &
         "total_residence_ps longest_residence_ps trajectory_extremum"
      DO item = 1, SIZE(actions)
         WRITE (unit, "(I8,1X,A,4I8,1X,A,2ES24.16,2I8,4ES24.16)") item, &
            TRIM(action_name(actions(item)%kind)), actions(item)%atom, TRIM(direction), entry_value, exit_value, &
            states(item)%event_count, states(item)%open_events, states(item)%first_entry_ps, &
            states(item)%total_residence_ps, states(item)%longest_residence_ps, states(item)%trajectory_extremum
      END DO
   END SUBROUTINE write_summary

   FUNCTION action_name(kind) RESULT(name)
      INTEGER, INTENT(IN)                                :: kind
      CHARACTER(LEN=8)                                   :: name

      SELECT CASE (kind)
      CASE (1)
         name = "distance"
      CASE (2)
         name = "angle"
      CASE (3)
         name = "torsion"
      CASE DEFAULT
         name = "unknown"
      END SELECT
   END FUNCTION action_name

END MODULE trajana_event_analysis
