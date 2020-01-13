!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interactive shell of CP2K
!> \note
!>       sample of a simple runner that uses the f77_interface
!>       it can be used to connect c programs, communicating through std-in/ std-out
!>
!>       positions are in angstrom, energies in evolt
!>
!>       commands:
!>       load filename:   loads the filename, returns the env_id, or -1 in case of error
!>       natom [env_id]:  returns the number of atoms in the environment env_id
!>                        (defaults to the last loaded)
!>       setpos [env_id]: sets the positions of the atoms, should be followed
!>                        by natom*3 (on a line) and then all the positions [angstrom]
!>       getpos [env_id]: gets the positions of the atoms, returns
!>                        natom*3 (on a line) and then all the positions [angstrom]
!>       calcE [env_id]:  calculate the energy and returns it (in eV)
!>       calcEF [env_id]: calculate the energy and forces and returns it,
!>                        first the energy on a line (in eV), then the natom*3 (on a line)
!>                        and finally all the values (in eV/angstrom)
!> \par History
!>      2019: Complete refactoring (Ole Schuett)
!>
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE cp2k_shell
   USE ISO_FORTRAN_ENV,                 ONLY: IOSTAT_END
   USE cp2k_info,                       ONLY: compile_arch,&
                                              compile_date,&
                                              compile_host,&
                                              compile_revision,&
                                              cp2k_home,&
                                              cp2k_version,&
                                              id_cp2k_version,&
                                              print_cp2k_license
   USE cp2k_runs,                       ONLY: run_input
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE f77_interface,                   ONLY: &
        calc_energy_force, calc_force, create_force_env, destroy_force_env, get_cell, get_energy, &
        get_force, get_natom, get_pos, get_stress_tensor, set_cell, set_pos
   USE input_constants,                 ONLY: id_development_version
   USE input_cp2k_read,                 ONLY: empty_initial_variables
   USE input_section_types,             ONLY: section_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE machine,                         ONLY: m_chdir,&
                                              m_flush,&
                                              m_getcwd,&
                                              m_getlog,&
                                              m_getpid,&
                                              m_hostnm
   USE message_passing,                 ONLY: mp_bcast
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE string_utilities,                ONLY: uppercase
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Queried by ASE. Increase version after bug-fixing or behavior changes.
   CHARACTER(LEN=*), PARAMETER                           :: CP2K_SHELL_VERSION = "4.0"

   TYPE cp2k_shell_type
      REAL(dp)                                           :: pos_fact = 1.0_dp
      REAL(dp)                                           :: e_fact = 1.0_dp
      LOGICAL                                            :: harsh = .FALSE.
      TYPE(cp_para_env_type), POINTER                    :: para_env => Null()
      CHARACTER(LEN=5)                                   :: units = "au"
      INTEGER                                            :: env_id = -1
      INTEGER                                            :: iw
   END TYPE cp2k_shell_type

   PUBLIC :: launch_cp2k_shell

CONTAINS

! **************************************************************************************************
!> \brief Launch the interactive CP2K shell.
!> \param input_declaration ...
! **************************************************************************************************
   SUBROUTINE launch_cp2k_shell(input_declaration)
      TYPE(section_type), POINTER                        :: input_declaration

      CHARACTER(LEN=default_string_length)               :: arg1, arg2, cmd
      TYPE(cp2k_shell_type)                              :: shell
      TYPE(cp_logger_type), POINTER                      :: logger

      logger => cp_get_default_logger()
      shell%para_env => logger%para_env
      shell%iw = cp_logger_get_default_io_unit()

      DO
         IF (.NOT. parse_next_line(shell, cmd, arg1, arg2)) EXIT

         ! dispatch command
         SELECT CASE (cmd)
         CASE ('HELP')
            CALL help_command(shell)
         CASE ('INFO', 'INFORMATION', 'LICENSE')
            CALL info_license_command(shell)
         CASE ('VERSION')
            CALL version_command(shell)
         CASE ('WRITE_FILE')
            CALL write_file_command(shell)
         CASE ('LAST_ENV_ID')
            CALL get_last_env_id(shell)
         CASE ('BG_LOAD', 'BGLOAD')
            CALL bg_load_command(shell, input_declaration, arg1)
         CASE ('LOAD')
            CALL load_command(shell, input_declaration, arg1, arg2)
         CASE ('DESTROY')
            CALL destroy_force_env_command(shell, arg1)
         CASE ('NATOM', 'N_ATOM')
            CALL get_natom_command(shell, arg1)
         CASE ('SETPOS', 'SET_POS')
            CALL set_pos_command(shell, arg1)
         CASE ('SETCELL', 'SET_CELL')
            CALL set_cell_command(shell, arg1)
         CASE ('GETCELL', 'GET_CELL')
            CALL get_cell_command(shell, arg1)
         CASE ('GETSTRESS', 'GET_STRESS')
            CALL get_stress_command(shell, arg1)
         CASE ('GETPOS', 'GET_POS')
            CALL get_pos_command(shell, arg1)
         CASE ('GETE', 'GET_E')
            CALL get_energy_command(shell, arg1)
         CASE ('EVALE', 'EVAL_E')
            CALL eval_energy_command(shell, arg1)
         CASE ('CALCE', 'CALC_E')
            CALL calc_energy_command(shell, arg1)
         CASE ('EVALEF', 'EVAL_EF')
            CALL eval_energy_force_command(shell, arg1)
         CASE ('GETF', 'GET_F')
            CALL get_forces_command(shell, arg1)
         CASE ('CALCEF', 'CALC_EF')
            CALL calc_energy_forces_command(shell, arg1)
         CASE ('RUN')
            CALL run_command(shell, input_declaration, arg1, arg2)
         CASE ('UNITS_EVA', 'UNITS_EV_A')
            CALL set_units_ev_a(shell)
         CASE ('UNITS_AU')
            CALL set_units_au(shell)
         CASE ('UNITS')
            CALL get_units(shell)
         CASE ('HARSH')
            shell%harsh = .TRUE.
         CASE ('PERMISSIVE')
            shell%harsh = .FALSE.
         CASE ('CD', 'CHDIR')
            CALL set_pwd_command(shell, arg1)
         CASE ('PWD', 'CWD')
            CALL get_pwd_command(shell)
         CASE ('EXIT')
            IF (shell%iw > 0) WRITE (shell%iw, '(a)') '* EXIT'
            EXIT
         CASE default
            CALL print_error('unknown command: '//cmd, shell)
         END SELECT
      END DO

   END SUBROUTINE launch_cp2k_shell

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param cmd ...
!> \param arg1 ...
!> \param arg2 ...
!> \return ...
! **************************************************************************************************
   FUNCTION parse_next_line(shell, cmd, arg1, arg2) RESULT(success)
      TYPE(cp2k_shell_type)                              :: shell
      CHARACTER(LEN=*), INTENT(out)                      :: cmd, arg1, arg2
      LOGICAL                                            :: success

      CHARACTER(LEN=10*default_string_length)            :: line
      INTEGER                                            :: i, iostat

      success = .TRUE.
      IF (shell%iw > 0) THEN
         WRITE (shell%iw, '("* READY")')
         CALL m_flush(shell%iw)
         READ (*, '(a)', iostat=iostat) line
         IF (iostat /= 0) THEN
            IF (iostat == IOSTAT_END) THEN
               WRITE (shell%iw, '(a)') '* EOF'
            END IF
            success = .FALSE. ! EOF
         END IF
      END IF
      CALL mp_bcast(success, shell%para_env%source, shell%para_env%group)
      IF (.NOT. success) RETURN
      CALL mp_bcast(line, shell%para_env%source, shell%para_env%group)

      ! extract command
      line = TRIM(line)
      DO i = 1, LEN_TRIM(line)
         IF (line(i:i) == ' ') EXIT
      ENDDO
      cmd = line(1:i)
      CALL uppercase(cmd)
      line = ADJUSTL(line(i:)) ! shift

      ! extract first arg
      DO i = 1, LEN_TRIM(line)
         IF (line(i:i) == ' ') EXIT
      ENDDO
      arg1 = line(1:i)
      line = ADJUSTL(line(i:)) ! shift

      ! extract second arg
      DO i = 1, LEN_TRIM(line)
         IF (line(i:i) == ' ') EXIT
      ENDDO
      arg2 = line(1:i)

      ! ignore remaining line
   END FUNCTION parse_next_line

! **************************************************************************************************
!> \brief Falls be env_id unchagned if not provided
!> \param str ...
!> \param shell ...
!> \return ...
! **************************************************************************************************
   FUNCTION parse_env_id(str, shell) RESULT(success)
      CHARACTER(LEN=*), INTENT(in)                       :: str
      TYPE(cp2k_shell_type)                              :: shell
      LOGICAL                                            :: success

      INTEGER                                            :: iostat

      success = .TRUE.
      IF (LEN_TRIM(str) > 0) THEN
         READ (str, *, iostat=iostat) shell%env_id
         IF (iostat /= 0) THEN
            shell%env_id = -1
            success = .FALSE.
            CALL print_error("parse_env_id failed", shell)
         ENDIF
      ELSE IF (shell%env_id < 1) THEN
         CALL print_error("last env_id not set", shell)
         success = .FALSE.
      ENDIF
      ! fallback: reuse last env_id
   END FUNCTION parse_env_id

! **************************************************************************************************
!> \brief ...
!> \param condition ...
!> \param message ...
!> \param shell ...
!> \return ...
! **************************************************************************************************
   FUNCTION my_assert(condition, message, shell) RESULT(success)
      LOGICAL, INTENT(in)                                :: condition
      CHARACTER(LEN=*), INTENT(in)                       :: message
      TYPE(cp2k_shell_type)                              :: shell
      LOGICAL                                            :: success

      success = condition
      IF (.NOT. success) THEN
         CALL print_error(message, shell)
      ENDIF
   END FUNCTION my_assert

! **************************************************************************************************
!> \brief ...
!> \param message ...
!> \param shell ...
! **************************************************************************************************
   SUBROUTINE print_error(message, shell)
      CHARACTER(LEN=*), INTENT(in)                       :: message
      TYPE(cp2k_shell_type)                              :: shell

      IF (shell%harsh) CPABORT(message)

      IF (shell%iw > 0) THEN
         WRITE (shell%iw, '("* ERROR ",a)') message
      END IF
   END SUBROUTINE print_error

! **************************************************************************************************
!> \brief ...
!> \param shell ...
! **************************************************************************************************
   SUBROUTINE help_command(shell)
      TYPE(cp2k_shell_type)                              :: shell

      IF (shell%iw > 0) THEN
         WRITE (shell%iw, *) 'Commands'
         WRITE (shell%iw, *) ' '
         WRITE (shell%iw, *) ' If there is [env_id] it means that an optional env_id can be given,'
         WRITE (shell%iw, *) ' if none is given it defaults to the last env_id loaded'
         WRITE (shell%iw, *) ' All commands are case insensitive.'
         WRITE (shell%iw, *) ' '
         WRITE (shell%iw, *) ' INFO: returns some information about cp2k.'
         WRITE (shell%iw, *) ' VERSION: returns shell version. (queried by ASE to assert features & bugfixes)'
         WRITE (shell%iw, *) ' WRITE_FILE: Writes content to a file (allows for using ASE over ssh).'
         WRITE (shell%iw, *) ' LOAD <inp-filename> [out-filename]: loads the filename, returns the env_id, or -1 in case of error'
         WRITE (shell%iw, *) '   out-filename is optional and defaults to <inp-filename>.out'
         WRITE (shell%iw, *) '   use "__STD_OUT__" for printing to the screen'
         WRITE (shell%iw, *) ' BG_LOAD <filename>: loads the filename, without returning the env_id'
         WRITE (shell%iw, *) ' LAST_ENV_ID: returns the env_id of the last environment loaded'
         WRITE (shell%iw, *) ' DESTROY [env_id]: destroys the given environment (last and default env'
         WRITE (shell%iw, *) '   might become invalid)'
         WRITE (shell%iw, *) ' NATOM [env_id]: returns the number of atoms in the environment env_id'
         WRITE (shell%iw, *) ' SET_POS [env_id]: sets the positions of the atoms, should be followed'
         WRITE (shell%iw, *) '   by natom*3 (on a line) and then all the positions. Returns the max'
         WRITE (shell%iw, *) '   change of the coordinates (useful to avoid extra calculations).'
         WRITE (shell%iw, *) ' SET_CELL [env_id]: sets the cell, should be followed by 9 numbers'
         WRITE (shell%iw, *) ' GET_CELL [env_id]: gets the cell vectors'
         WRITE (shell%iw, *) ' GET_STRESS [env_id]: gets the stress tensor of the last calculation on env_id'
         WRITE (shell%iw, *) ' GET_POS [env_id]: gets the positions of the atoms, returns'
         WRITE (shell%iw, *) '   natom*3 (on a line) and then all the positions then "* END" '
         WRITE (shell%iw, *) '   (alone on a line)'
         WRITE (shell%iw, *) ' GET_E [env_id]: gets the energy of the last calculation on env_id'
         WRITE (shell%iw, *) ' GET_F [env_id]: gets the forces on the atoms,of the last calculation on '
         WRITE (shell%iw, *) '   env_id, if only the energy was calculated the content is undefined. Returns'
         WRITE (shell%iw, *) '   natom*3 (on a line) and then all the forces then "* END" (alone on'
         WRITE (shell%iw, *) '   a line)'
         WRITE (shell%iw, *) ' CALC_E [env_id]: calculate the energy and returns it'
         WRITE (shell%iw, *) ' EVAL_E [env_id]: calculate the energy (without returning it)'
         WRITE (shell%iw, *) ' CALC_EF [env_id]: calculate energy and forces and returns them,'
         WRITE (shell%iw, *) '   first the energy on a line, then the natom*3 (on a line)'
         WRITE (shell%iw, *) '   and finally all the values and "* END" (alone on a line)'
         WRITE (shell%iw, *) ' EVAL_EF [env_id]: calculate the energy and forces (without returning them)'
         WRITE (shell%iw, *) ' RUN <inp-filename> <out-filename>: run the given input file'
         WRITE (shell%iw, *) ' HARSH: stops on any error'
         WRITE (shell%iw, *) ' PERMISSIVE: stops only on serious errors'
         WRITE (shell%iw, *) ' UNITS: returns the units used for energy and position'
         WRITE (shell%iw, *) ' UNITS_EV_A: sets the units to electron volt (energy)  and Angstrom (positions)'
         WRITE (shell%iw, *) ' UNITS_AU: sets the units atomic units'
         WRITE (shell%iw, *) ' CD <dir>: change working directory'
         WRITE (shell%iw, *) ' PWD: print working directory'
         WRITE (shell%iw, *) ' EXIT: Quit the shell'
         WRITE (shell%iw, *) ' HELP: writes the present help'
         CALL m_flush(shell%iw)
      END IF
   END SUBROUTINE help_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
! **************************************************************************************************
   SUBROUTINE info_license_command(shell)
      TYPE(cp2k_shell_type)                              :: shell

      CHARACTER(LEN=default_path_length)                 :: cwd
      CHARACTER(LEN=default_string_length)               :: dev_flag, host_name, user_name
      INTEGER                                            :: pid

      IF (shell%iw > 0) THEN
         CALL m_getcwd(cwd)
         CALL m_getpid(pid)
         CALL m_getlog(user_name)
         CALL m_hostnm(host_name)
         dev_flag = ""
         IF (id_cp2k_version == id_development_version) dev_flag = " (Development Version)"
         WRITE (UNIT=shell%iw, FMT="(A,A)") &
            " PROGRAM STARTED ON ", TRIM(host_name)
         WRITE (UNIT=shell%iw, FMT="(A,A)") &
            " PROGRAM STARTED BY ", TRIM(user_name)
         WRITE (UNIT=shell%iw, FMT="(A,i10)") &
            " PROGRAM PROCESS ID ", pid
         WRITE (UNIT=shell%iw, FMT="(A,A)") &
            " PROGRAM STARTED IN ", TRIM(cwd)
         WRITE (UNIT=shell%iw, FMT="(/,T2,A,T31,A50)") &
            "CP2K| version string: ", &
            ADJUSTR(TRIM(cp2k_version)//TRIM(dev_flag))
         WRITE (UNIT=shell%iw, FMT="(T2,A,T41,A40)") &
            "CP2K| source code revision number:", &
            ADJUSTR(compile_revision)
         WRITE (UNIT=shell%iw, FMT="(T2,A,T41,A40)") &
            "CP2K| is freely available from ", &
            ADJUSTR(TRIM(cp2k_home))
         WRITE (UNIT=shell%iw, FMT="(T2,A,T31,A50)") &
            "CP2K| Program compiled at", &
            ADJUSTR(compile_date(1:MIN(50, LEN(compile_date))))
         WRITE (UNIT=shell%iw, FMT="(T2,A,T31,A50)") &
            "CP2K| Program compiled on", &
            ADJUSTR(compile_host(1:MIN(50, LEN(compile_host))))
         WRITE (UNIT=shell%iw, FMT="(T2,A,T31,A50)") &
            "CP2K| Program compiled for", &
            ADJUSTR(compile_arch(1:MIN(50, LEN(compile_arch))))

         CALL print_cp2k_license(shell%iw)
         CALL m_flush(shell%iw)
      END IF

   END SUBROUTINE info_license_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
! **************************************************************************************************
   SUBROUTINE version_command(shell)
      TYPE(cp2k_shell_type)                              :: shell

      IF (shell%iw > 0) THEN
         WRITE (shell%iw, '(a,a)') "CP2K Shell Version: ", CP2K_SHELL_VERSION
         CALL m_flush(shell%iw)
      END IF
   END SUBROUTINE version_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
! **************************************************************************************************
   SUBROUTINE write_file_command(shell)
      TYPE(cp2k_shell_type)                              :: shell

      CHARACTER(LEN=10*default_string_length)            :: line
      CHARACTER(LEN=default_path_length)                 :: out_filename
      INTEGER                                            :: file_unit, i, iostat, n_lines

      IF (shell%iw > 0) THEN
         READ (*, '(a)', iostat=iostat) out_filename
         IF (iostat /= 0) CPABORT('WRITE_FILE bad filename')
         READ (*, *, iostat=iostat) n_lines
         IF (iostat /= 0) CPABORT('WRITE_FILE bad number of lines')
         CALL open_file(file_name=TRIM(out_filename), unit_number=file_unit, &
                        file_status="UNKNOWN", file_form="FORMATTED", file_action="WRITE")
         DO i = 1, n_lines
            READ (*, '(a)', iostat=iostat) line
            IF (iostat /= 0) CPABORT('WRITE_FILE read error')
            WRITE (file_unit, '(a)', iostat=iostat) TRIM(line)
            IF (iostat /= 0) CPABORT('WRITE_FILE write error')
         END DO
         READ (*, '(a)', iostat=iostat) line
         IF (iostat /= 0) CPABORT('WRITE_FILE read error')
         IF (TRIM(line) /= "*END") CPABORT('WRITE_FILE bad end delimiter')
         CALL close_file(unit_number=file_unit)
      END IF
   END SUBROUTINE write_file_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
! **************************************************************************************************
   SUBROUTINE get_last_env_id(shell)
      TYPE(cp2k_shell_type)                              :: shell

      IF (shell%iw > 0) THEN
         WRITE (shell%iw, '(i10)') shell%env_id
         CALL m_flush(shell%iw)
      END IF
   END SUBROUTINE get_last_env_id

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param input_declaration ...
!> \param arg1 ...
! **************************************************************************************************
   SUBROUTINE bg_load_command(shell, input_declaration, arg1)
      TYPE(cp2k_shell_type)                              :: shell
      TYPE(section_type), POINTER                        :: input_declaration
      CHARACTER(LEN=*)                                   :: arg1

      INTEGER                                            :: ierr

      IF (.NOT. my_assert(LEN_TRIM(arg1) > 0, "file argument missing", shell)) RETURN
      CALL create_force_env(new_env_id=shell%env_id, &
                            input_declaration=input_declaration, &
                            input_path=TRIM(arg1), &
                            output_path=TRIM(arg1)//'.out', &
                            owns_out_unit=.TRUE., ierr=ierr)
      IF (ierr /= 0) THEN
         shell%env_id = -1
         CALL print_error("create_force_env failed", shell)
      ENDIF
   END SUBROUTINE bg_load_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param input_declaration ...
!> \param arg1 ...
!> \param arg2 ...
! **************************************************************************************************
   SUBROUTINE load_command(shell, input_declaration, arg1, arg2)
      TYPE(cp2k_shell_type)                              :: shell
      TYPE(section_type), POINTER                        :: input_declaration
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1, arg2

      CHARACTER(LEN=default_path_length)                 :: inp_filename, out_filename
      INTEGER                                            :: ierr

      IF (.NOT. my_assert(LEN_TRIM(arg1) > 0, "file argument missing", shell)) RETURN
      inp_filename = arg1
      out_filename = TRIM(inp_filename)//'.out'
      IF (LEN_TRIM(arg2) > 0) out_filename = arg2
      CALL create_force_env(new_env_id=shell%env_id, &
                            input_declaration=input_declaration, &
                            input_path=inp_filename, &
                            output_path=out_filename, &
                            owns_out_unit=.TRUE., ierr=ierr)
      IF (ierr /= 0) THEN
         shell%env_id = -1
         CALL print_error("create_force_env failed", shell)
      ELSE IF (shell%iw > 0) THEN
         WRITE (shell%iw, '(i10)') shell%env_id
         CALL m_flush(shell%iw)
      END IF
   END SUBROUTINE load_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param arg1 ...
! **************************************************************************************************
   SUBROUTINE destroy_force_env_command(shell, arg1)
      TYPE(cp2k_shell_type)                              :: shell
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1

      INTEGER                                            :: ierr

      IF (.NOT. parse_env_id(arg1, shell)) RETURN
      CALL destroy_force_env(shell%env_id, ierr)
      shell%env_id = -1
      IF (ierr /= 0) CALL print_error('destroy_force_env failed', shell)
   END SUBROUTINE destroy_force_env_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param arg1 ...
! **************************************************************************************************
   SUBROUTINE get_natom_command(shell, arg1)
      TYPE(cp2k_shell_type)                              :: shell
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1

      INTEGER                                            :: ierr, iostat, n_atom

      IF (.NOT. parse_env_id(arg1, shell)) RETURN
      CALL get_natom(shell%env_id, n_atom, ierr)
      IF (.NOT. my_assert(ierr == 0, 'get_natom failed', shell)) RETURN
      IF (shell%iw > 0) THEN
         WRITE (shell%iw, '(i10)', iostat=iostat) n_atom
         CALL m_flush(shell%iw)
      END IF
   END SUBROUTINE get_natom_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param arg1 ...
! **************************************************************************************************
   SUBROUTINE set_pos_command(shell, arg1)
      TYPE(cp2k_shell_type)                              :: shell
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1

      CHARACTER(LEN=10*default_string_length)            :: line
      INTEGER                                            :: i, ierr, iostat, n_atom
      REAL(KIND=dp)                                      :: max_change
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: old_pos, pos

      IF (.NOT. parse_env_id(arg1, shell)) RETURN
      CALL get_natom(shell%env_id, n_atom, ierr)
      IF (.NOT. my_assert(ierr == 0, 'get_natom failed', shell)) RETURN
      ALLOCATE (pos(3*n_atom), old_pos(3*n_atom))
      IF (shell%iw > 0) THEN
         READ (*, *, iostat=iostat) n_atom
         IF (.NOT. my_assert(iostat == 0, 'setpos read n_atom failed', shell)) RETURN
         IF (.NOT. my_assert(n_atom == SIZE(pos), 'setpos invalid number of atoms', shell)) RETURN
         READ (*, *, iostat=iostat) pos
         IF (.NOT. my_assert(iostat == 0, 'setpos read coords failed', shell)) RETURN
         pos(:) = pos(:)/shell%pos_fact
         READ (*, '(a)', iostat=iostat) line
         IF (.NOT. my_assert(iostat == 0, 'setpos read END failed', shell)) RETURN
         CALL uppercase(line)
         IF (.NOT. my_assert(TRIM(line) == '*END', 'missing *END', shell)) RETURN
      END IF
      CALL mp_bcast(pos, shell%para_env%source, shell%para_env%group)
      CALL get_pos(shell%env_id, old_pos, n_el=3*n_atom, ierr=ierr)
      IF (.NOT. my_assert(ierr == 0, 'get_pos error', shell)) RETURN
      CALL set_pos(shell%env_id, new_pos=pos, n_el=3*n_atom, ierr=ierr)
      IF (.NOT. my_assert(ierr == 0, 'set_pos error', shell)) RETURN
      max_change = 0.0_dp
      DO i = 1, SIZE(pos)
         max_change = MAX(max_change, ABS(pos(i) - old_pos(i)))
      END DO
      DEALLOCATE (pos, old_pos)
      IF (shell%iw > 0) THEN
         WRITE (shell%iw, '(ES22.13)') max_change*shell%pos_fact
         CALL m_flush(shell%iw)
      END IF
   END SUBROUTINE set_pos_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param arg1 ...
! **************************************************************************************************
   SUBROUTINE set_cell_command(shell, arg1)
      TYPE(cp2k_shell_type)                              :: shell
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1

      INTEGER                                            :: ierr, iostat
      REAL(KIND=dp), DIMENSION(3, 3)                     :: cell

      IF (.NOT. parse_env_id(arg1, shell)) RETURN
      IF (shell%iw > 0) THEN
         READ (*, *, iostat=iostat) cell
         IF (.NOT. my_assert(iostat == 0, 'setcell read failed', shell)) RETURN
         cell(:, :) = cell(:, :)/shell%pos_fact
      END IF
      CALL mp_bcast(cell, shell%para_env%source, shell%para_env%group)
      CALL set_cell(shell%env_id, new_cell=cell, ierr=ierr)
      IF (.NOT. my_assert(ierr == 0, 'set_cell failed', shell)) RETURN
   END SUBROUTINE set_cell_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param arg1 ...
! **************************************************************************************************
   SUBROUTINE get_cell_command(shell, arg1)
      TYPE(cp2k_shell_type)                              :: shell
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1

      INTEGER                                            :: ierr
      REAL(KIND=dp), DIMENSION(3, 3)                     :: cell

      IF (.NOT. parse_env_id(arg1, shell)) RETURN
      CALL get_cell(shell%env_id, cell=cell, ierr=ierr)
      IF (.NOT. my_assert(ierr == 0, 'get_cell failed', shell)) RETURN
      cell(:, :) = cell(:, :)*shell%pos_fact
      IF (shell%iw > 0) THEN
         WRITE (shell%iw, '(9ES22.13)') cell
         CALL m_flush(shell%iw)
      END IF
   END SUBROUTINE get_cell_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param arg1 ...
! **************************************************************************************************
   SUBROUTINE get_stress_command(shell, arg1)
      TYPE(cp2k_shell_type)                              :: shell
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1

      INTEGER                                            :: ierr
      REAL(KIND=dp), DIMENSION(3, 3)                     :: stress_tensor

      IF (.NOT. parse_env_id(arg1, shell)) RETURN
      CALL get_stress_tensor(shell%env_id, stress_tensor=stress_tensor, ierr=ierr)
      IF (.NOT. my_assert(ierr == 0, 'get_stress_tensor failed', shell)) RETURN
      stress_tensor(:, :) = stress_tensor(:, :)*(shell%e_fact/shell%pos_fact**3)
      IF (shell%iw > 0) THEN
         WRITE (shell%iw, '(9ES22.13)') stress_tensor
         CALL m_flush(shell%iw)
      END IF
   END SUBROUTINE get_stress_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param arg1 ...
! **************************************************************************************************
   SUBROUTINE get_pos_command(shell, arg1)
      TYPE(cp2k_shell_type)                              :: shell
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1

      INTEGER                                            :: ierr, n_atom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: pos

      IF (.NOT. parse_env_id(arg1, shell)) RETURN
      CALL get_natom(shell%env_id, n_atom, ierr)
      IF (.NOT. my_assert(ierr == 0, 'get_natom failed', shell)) RETURN
      ALLOCATE (pos(3*n_atom))
      CALL get_pos(shell%env_id, pos=pos, n_el=3*n_atom, ierr=ierr)
      IF (.NOT. my_assert(ierr == 0, 'get_pos failed', shell)) RETURN
      IF (shell%iw > 0) THEN
         WRITE (shell%iw, '(i10)') 3*n_atom
         WRITE (shell%iw, '(3ES22.13)') pos(:)*shell%pos_fact
         WRITE (shell%iw, '(a)') "* END"
         CALL m_flush(shell%iw)
      END IF
      DEALLOCATE (pos)
   END SUBROUTINE get_pos_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param arg1 ...
! **************************************************************************************************
   SUBROUTINE get_energy_command(shell, arg1)
      TYPE(cp2k_shell_type)                              :: shell
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1

      INTEGER                                            :: ierr
      REAL(KIND=dp)                                      :: e_pot

      IF (.NOT. parse_env_id(arg1, shell)) RETURN
      CALL get_energy(shell%env_id, e_pot, ierr)
      IF (.NOT. my_assert(ierr == 0, 'get_energy failed', shell)) RETURN
      IF (shell%iw > 0) THEN
         WRITE (shell%iw, '(ES22.13)') e_pot*shell%e_fact
         CALL m_flush(shell%iw)
      END IF
   END SUBROUTINE get_energy_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param arg1 ...
! **************************************************************************************************
   SUBROUTINE eval_energy_command(shell, arg1)
      TYPE(cp2k_shell_type)                              :: shell
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1

      INTEGER                                            :: ierr

      IF (.NOT. parse_env_id(arg1, shell)) RETURN
      CALL calc_energy_force(shell%env_id, calc_force=.FALSE., ierr=ierr)
      IF (ierr /= 0) CALL print_error('calc_energy_force failed', shell)
   END SUBROUTINE eval_energy_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param arg1 ...
! **************************************************************************************************
   SUBROUTINE calc_energy_command(shell, arg1)
      TYPE(cp2k_shell_type)                              :: shell
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1

      INTEGER                                            :: ierr
      REAL(KIND=dp)                                      :: e_pot

      IF (.NOT. parse_env_id(arg1, shell)) RETURN
      CALL calc_energy_force(shell%env_id, calc_force=.FALSE., ierr=ierr)
      IF (.NOT. my_assert(ierr == 0, 'calc_energy_force failed', shell)) RETURN
      CALL get_energy(shell%env_id, e_pot, ierr)
      IF (.NOT. my_assert(ierr == 0, 'get_energy failed', shell)) RETURN
      IF (shell%iw > 0) THEN
         WRITE (shell%iw, '(ES22.13)') e_pot*shell%e_fact
         CALL m_flush(shell%iw)
      END IF
   END SUBROUTINE calc_energy_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param arg1 ...
! **************************************************************************************************
   SUBROUTINE eval_energy_force_command(shell, arg1)
      TYPE(cp2k_shell_type)                              :: shell
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1

      INTEGER                                            :: ierr

      IF (.NOT. parse_env_id(arg1, shell)) RETURN
      CALL calc_energy_force(shell%env_id, calc_force=.TRUE., ierr=ierr)
      IF (ierr /= 0) CALL print_error('calc_energy_force failed', shell)
   END SUBROUTINE eval_energy_force_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param arg1 ...
! **************************************************************************************************
   SUBROUTINE get_forces_command(shell, arg1)
      TYPE(cp2k_shell_type)                              :: shell
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1

      INTEGER                                            :: ierr, n_atom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: forces

      IF (.NOT. parse_env_id(arg1, shell)) RETURN
      CALL get_natom(shell%env_id, n_atom, ierr)
      IF (.NOT. my_assert(ierr == 0, 'get_natom failed', shell)) RETURN
      ALLOCATE (forces(3*n_atom))
      CALL get_force(shell%env_id, frc=forces, n_el=3*n_atom, ierr=ierr)
      IF (.NOT. my_assert(ierr == 0, 'get_force failed', shell)) RETURN
      forces(:) = forces(:)*(shell%e_fact/shell%pos_fact)
      IF (shell%iw > 0) THEN
         WRITE (shell%iw, '(i10)') 3*n_atom
         WRITE (shell%iw, '(3ES22.13)') forces
         WRITE (shell%iw, '("* END")')
         CALL m_flush(shell%iw)
      END IF
      DEALLOCATE (forces)
   END SUBROUTINE get_forces_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param arg1 ...
! **************************************************************************************************
   SUBROUTINE calc_energy_forces_command(shell, arg1)
      TYPE(cp2k_shell_type)                              :: shell
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1

      INTEGER                                            :: ierr, n_atom
      REAL(KIND=dp)                                      :: e_pot
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: forces

      IF (.NOT. parse_env_id(arg1, shell)) RETURN
      CALL calc_energy_force(shell%env_id, calc_force=.TRUE., ierr=ierr)
      IF (.NOT. my_assert(ierr == 0, 'calc_energy_force failed', shell)) RETURN
      CALL get_energy(shell%env_id, e_pot, ierr)
      IF (.NOT. my_assert(ierr == 0, 'get_energy failed', shell)) RETURN
      CALL get_natom(shell%env_id, n_atom, ierr)
      IF (.NOT. my_assert(ierr == 0, 'get_natom failed', shell)) RETURN
      ALLOCATE (forces(3*n_atom))
      CALL get_force(shell%env_id, frc=forces, n_el=3*n_atom, ierr=ierr)
      IF (.NOT. my_assert(ierr == 0, 'get_energy failed', shell)) RETURN
      IF (shell%iw > 0) THEN
         WRITE (shell%iw, '(ES22.13)') e_pot*shell%e_fact
         WRITE (shell%iw, '(i10)') 3*n_atom
         WRITE (shell%iw, '(3ES22.13)') forces*(shell%e_fact/shell%pos_fact)
         WRITE (shell%iw, '("* END")')
         CALL m_flush(shell%iw)
      END IF
      DEALLOCATE (forces)
   END SUBROUTINE calc_energy_forces_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param input_declaration ...
!> \param arg1 ...
!> \param arg2 ...
! **************************************************************************************************
   SUBROUTINE run_command(shell, input_declaration, arg1, arg2)
      TYPE(cp2k_shell_type)                              :: shell
      TYPE(section_type), POINTER                        :: input_declaration
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1, arg2

      IF (.NOT. my_assert(LEN_TRIM(arg1) > 0, "input-file argument missing", shell)) RETURN
      IF (.NOT. my_assert(LEN_TRIM(arg2) > 0, "input-file argument missing", shell)) RETURN
      CALL run_input(input_declaration, arg1, arg2, empty_initial_variables)
   END SUBROUTINE run_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
! **************************************************************************************************
   SUBROUTINE set_units_ev_a(shell)
      TYPE(cp2k_shell_type)                              :: shell

      shell%e_fact = evolt
      shell%pos_fact = angstrom
      shell%units = 'eV_A'
   END SUBROUTINE set_units_ev_a

! **************************************************************************************************
!> \brief ...
!> \param shell ...
! **************************************************************************************************
   SUBROUTINE set_units_au(shell)
      TYPE(cp2k_shell_type)                              :: shell

      shell%e_fact = 1.0_dp
      shell%pos_fact = 1.0_dp
      shell%units = 'au'
   END SUBROUTINE set_units_au

! **************************************************************************************************
!> \brief ...
!> \param shell ...
! **************************************************************************************************
   SUBROUTINE get_units(shell)
      TYPE(cp2k_shell_type)                              :: shell

      IF (shell%iw > 0) THEN
         WRITE (shell%iw, '(a)') TRIM(shell%units)
         CALL m_flush(shell%iw)
      END IF
   END SUBROUTINE get_units

! **************************************************************************************************
!> \brief ...
!> \param shell ...
!> \param arg1 ...
! **************************************************************************************************
   SUBROUTINE set_pwd_command(shell, arg1)
      TYPE(cp2k_shell_type)                              :: shell
      CHARACTER(LEN=*), INTENT(IN)                       :: arg1

      INTEGER                                            :: ierr

      IF (.NOT. my_assert(LEN_TRIM(arg1) > 0, 'missing directory', shell)) RETURN
      CALL m_chdir(arg1, ierr)
      IF (ierr /= 0) CALL print_error('changing directory failed', shell)
   END SUBROUTINE set_pwd_command

! **************************************************************************************************
!> \brief ...
!> \param shell ...
! **************************************************************************************************
   SUBROUTINE get_pwd_command(shell)
      TYPE(cp2k_shell_type)                              :: shell

      CHARACTER(LEN=default_path_length)                 :: cwd

      IF (shell%iw > 0) THEN
         CALL m_getcwd(cwd)
         WRITE (shell%iw, '(a)') TRIM(cwd)
      END IF
   END SUBROUTINE get_pwd_command

END MODULE cp2k_shell
