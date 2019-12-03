!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Main program of CP2K
!> \par Copyright
!>    CP2K: A general program to perform molecular dynamics simulations
!>    Copyright (C) 2000, 2001, 2002, 2003  CP2K developers group
!>    Copyright (C) 2004, 2005, 2006, 2007  CP2K developers group
!>    Copyright (C) 2008, 2009, 2010, 2011  CP2K developers group
!>    Copyright (C) 2012, 2013, 2014, 2015  CP2K developers group
!>    Copyright (C) 2016                    CP2K developers group
!> \par
!>    This program is free software; you can redistribute it and/or modify
!>    it under the terms of the GNU General Public License as published by
!>    the Free Software Foundation; either version 2 of the License, or
!>    (at your option) any later version.
!> \par
!>    This program is distributed in the hope that it will be useful,
!>    but WITHOUT ANY WARRANTY; without even the implied warranty of
!>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!>    GNU General Public License for more details.
!> \par
!>    You should have received a copy of the GNU General Public License
!>    along with this program; if not, write to the Free Software
!>    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!>    Boston, MA  02110-1301, USA.
!> \par
!>    See also https://www.fsf.org/licensing/licenses/gpl.html
!> \par
!>    CP2K, including its sources and pointers to the authors
!>    can be found at  https://www.cp2k.org/
!> \note
!>       should be kept as lean as possible.
!>       see cp2k_run for more comments
!> \author Joost VandeVondele
! **************************************************************************************************
PROGRAM cp2k

   USE cp2k_info,                       ONLY: compile_revision,&
                                              cp2k_flags,&
                                              cp2k_version,&
                                              id_cp2k_version,&
                                              print_cp2k_license
   USE cp2k_runs,                       ONLY: run_input,&
                                              write_xml_file
   USE cp2k_shell,                      ONLY: launch_cp2k_shell
   USE f77_interface,                   ONLY: check_input,&
                                              default_para_env,&
                                              finalize_cp2k,&
                                              init_cp2k
   USE input_constants,                 ONLY: id_development_version
   USE input_cp2k,                      ONLY: create_cp2k_root_section
   USE input_section_types,             ONLY: section_release,&
                                              section_type
   USE iso_fortran_env,                 ONLY: compiler_options,&
                                              compiler_version
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length
   USE machine,                         ONLY: default_output_unit
#include "../base/base_uses.f90"

   IMPLICIT NONE

   CHARACTER(LEN=default_path_length)   :: input_file_name, output_file_name, &
                                           arg_att, command
   CHARACTER(LEN=default_string_length) :: dev_flag
   CHARACTER(LEN=:), ALLOCATABLE        :: compiler_options_string
   INTEGER                              :: output_unit, l, i
   INTEGER                              :: ierr, i_arg
   LOGICAL                              :: check, usage, echo_input, command_line_error
   LOGICAL                              :: run_it, force_run, has_input, xml, print_version, print_license, shell_mode
   TYPE(section_type), POINTER          :: input_declaration

   NULLIFY (input_declaration)

   ! output goes to the screen by default
   output_unit = default_output_unit

   ! set default behaviour for the command line switches
   check = .FALSE.
   usage = .FALSE.
   echo_input = .FALSE.
   has_input = .FALSE.
   run_it = .TRUE.
   shell_mode = .FALSE.
   force_run = .FALSE.
   print_version = .FALSE.
   print_license = .FALSE.
   command_line_error = .FALSE.
   xml = .FALSE.
   input_file_name = "Missing input file name" ! no default
   output_file_name = "__STD_OUT__" ! by default we go to std_out

   ! get command and strip path
   CALL GET_COMMAND(command, status=ierr)
   CPASSERT(ierr == 0)
   l = LEN_TRIM(command)
   DO i = l, 1, -1
      IF (command(i:i) == "/" .OR. command(i:i) == "\") EXIT
   END DO
   command = command(i + 1:l)

   ! check if binary was invoked as cp2k_shell
   IF (command(1:10) == "cp2k_shell") THEN
      shell_mode = .TRUE.
      run_it = .FALSE.
   ELSE IF (COMMAND_ARGUMENT_COUNT() < 1) THEN
      WRITE (output_unit, "(/,T2,A)") &
         "ERROR: At least one command line argument must be specified"
      command_line_error = .TRUE.
   END IF

   i_arg = 0
   arg_loop: DO WHILE (i_arg < COMMAND_ARGUMENT_COUNT())
      i_arg = i_arg + 1
      CALL GET_COMMAND_ARGUMENT(i_arg, arg_att, status=ierr)
      CPASSERT(ierr == 0)
      SELECT CASE (arg_att)
      CASE ("--check", "-c")
         check = .TRUE.
         run_it = .FALSE.
         echo_input = .FALSE.
      CASE ("--echo", "-e")
         check = .TRUE.
         run_it = .FALSE.
         echo_input = .TRUE.
      CASE ("-v", "--version")
         print_version = .TRUE.
         run_it = .FALSE.
      CASE ("--license")
         print_license = .TRUE.
         run_it = .FALSE.
      CASE ("--run", "-r")
         force_run = .TRUE.
      CASE ("--shell", "-s")
         shell_mode = .TRUE.
         run_it = .FALSE.
      CASE ("-help", "--help", "-h")
         usage = .TRUE.
         run_it = .FALSE.
      CASE ("-i")
         i_arg = i_arg + 1
         CALL GET_COMMAND_ARGUMENT(i_arg, arg_att, status=ierr)
         CPASSERT(ierr == 0)
         ! argument does not start with a - it is an filename
         IF (.NOT. arg_att(1:1) == "-") THEN
            input_file_name = arg_att
            has_input = .TRUE.
         ELSE
            WRITE (output_unit, "(/,T2,A)") &
               "ERROR: The input file name"//TRIM(arg_att)//" starts with -"
            command_line_error = .TRUE.
            EXIT arg_loop
         ENDIF
      CASE ("-o")
         i_arg = i_arg + 1
         CALL GET_COMMAND_ARGUMENT(i_arg, arg_att, status=ierr)
         CPASSERT(ierr == 0)
         ! argument does not start with a - it is an filename
         IF (.NOT. arg_att(1:1) == "-") THEN
            output_file_name = arg_att
         ELSE
            WRITE (output_unit, "(/,T2,A)") &
               "ERROR: The output file name"//TRIM(arg_att)//" starts with -"
            command_line_error = .TRUE.
            EXIT arg_loop
         ENDIF
      CASE ("--xml")
         xml = .TRUE.
         run_it = .FALSE.
      CASE default
         ! if the last argument does not start with a - it is an input filename
         !MK in order to digest the additional flags of mpirun
         IF ((.NOT. has_input) .AND. &
             (i_arg == COMMAND_ARGUMENT_COUNT()) .AND. &
             (.NOT. arg_att(1:1) == "-")) THEN
            input_file_name = arg_att
            has_input = .TRUE.
         ELSEIF (has_input .AND. &
                 (.NOT. arg_att(1:1) == "-")) THEN
            WRITE (output_unit, "(/,T2,A)") &
               "Error: Tried to specify two input files"
            command_line_error = .TRUE.
            EXIT arg_loop
         END IF
      END SELECT
   END DO arg_loop
   IF ((run_it .OR. force_run .OR. check .OR. echo_input) .AND. &
       (.NOT. has_input) .AND. (.NOT. command_line_error)) THEN
      WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
         "ERROR: An input file name is required"
      command_line_error = .TRUE.
   END IF

   CALL init_cp2k(init_mpi=.TRUE., ierr=ierr)

   IF (ierr == 0) THEN
      ! some first info concerning how to run CP2K

      IF (usage .OR. command_line_error) THEN
         IF (default_para_env%ionode) THEN
            l = LEN_TRIM(command)
            WRITE (UNIT=output_unit, FMT="(/,(T2,A))") &
               TRIM(command)//" [-c|--check] [-e|--echo] [-h|--help]", &
               REPEAT(" ", l)//" [-i] <input_file>", &
               REPEAT(" ", l)//" [-mpi-mapping|--mpi-mapping] <method>", &
               REPEAT(" ", l)//" [-o] <output_file>", &
               REPEAT(" ", l)//" [-r|-run] [-s|--shell] [--xml]"
            WRITE (UNIT=output_unit, FMT="(/,T2,A,/,/,T2,A,/,/,T2,A,/,/,(T3,A))") &
               "starts the CP2K program, see <https://www.cp2k.org/>", &
               "The easiest way is "//TRIM(command)//" <input_file>", &
               "The following options can be used:", &
               "-i <input_file>   : provides an input file name, if it is the last", &
               "                    argument, the -i flag is not needed", &
               "-o <output_file>  : provides an output file name [default: screen]"
            WRITE (UNIT=output_unit, FMT="(/,T2,A,/,/,(T3,A))") &
               "These switches skip the simulation, unless [-r|-run] is specified:", &
               "--check, -c       : performs a syntax check of the <input_file>", &
               "--echo, -e        : echoes the <input_file>, and make all defaults explicit", &
               "                    The input is also checked, but only a failure is reported", &
               "--help, -h        : writes this message", &
               "--license         : prints the CP2K license", &
               "--mpi-mapping     : applies a given MPI reordering to CP2K", &
               "--run, -r         : forces a CP2K run regardless of other specified flags", &
               "--shell, -s       : start interactive shell mode", &
               "--version, -v     : prints the CP2K version and the revision number", &
               "--xml             : dumps the whole CP2K input structure as a XML file", &
               "                    xml2htm generates a HTML manual from this XML file", &
               ""
         END IF
      END IF

      IF (.NOT. command_line_error) THEN

         ! write the version string
         IF (print_version) THEN
            IF (default_para_env%ionode) THEN
               dev_flag = ""
               IF (id_cp2k_version == id_development_version) &
                  dev_flag = " (Development Version)"
               WRITE (output_unit, "(T2,A)") &
                  cp2k_version//TRIM(dev_flag), &
                  "Source code revision "//TRIM(compile_revision), &
                  TRIM(cp2k_flags())
               compiler_options_string = compiler_options()
               WRITE (output_unit, "(T2,A,A)") "compiler: ", compiler_version()
               WRITE (output_unit, "(T2,A)") "compiler options:"
               DO i = 0, (LEN(compiler_options_string) - 1)/68
                  WRITE (output_unit, "(T4,A)") &
                     compiler_options_string(i*68 + 1:MIN(LEN(compiler_options_string), (i + 1)*68))
               END DO
            ENDIF
         ENDIF

         ! write the license
         IF (print_license) THEN
            IF (default_para_env%ionode) THEN
               CALL print_cp2k_license(output_unit)
            ENDIF
         ENDIF

         IF (xml) THEN
            IF (default_para_env%ionode) THEN
               CALL write_xml_file()
            END IF
         END IF

         CALL create_cp2k_root_section(input_declaration)

         IF (check) THEN
            CALL check_input(input_declaration, input_file_name, output_file_name, &
                             echo_input=echo_input, ierr=ierr)
            IF (default_para_env%ionode) THEN
               IF (ierr == 0) THEN
                  IF (.NOT. echo_input) THEN
                     WRITE (output_unit, "(A)") "SUCCESS, the input could be parsed correctly."
                     WRITE (output_unit, "(A)") "         This is does not guarantee that this input is meaningful"
                     WRITE (output_unit, "(A)") "         or will run successfully"
                  ENDIF
               ELSE
                  WRITE (output_unit, "(A)") "ERROR,   the input could *NOT* be parsed correctly."
                  WRITE (output_unit, "(A)") "         Please, check and correct it"
               END IF
            END IF
         ENDIF

         IF (shell_mode) THEN
            CALL launch_cp2k_shell(input_declaration)
         ENDIF

         IF (run_it .OR. force_run) THEN
            CALL run_input(input_declaration, input_file_name, output_file_name)
         END IF

         CALL section_release(input_declaration)
      ENDIF
   ELSE
      WRITE (UNIT=output_unit, FMT="(/,A)") "initial setup (MPI ?) error"
   END IF

   ! and the final cleanup
   CALL finalize_cp2k(finalize_mpi=.TRUE., ierr=ierr)
   CPASSERT(ierr == 0)

END PROGRAM
