!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief starts the program
!> \par Copyright
!>    CP2K: A general program to perform molecular dynamics simulations         !
!>    Copyright (C) 2000, 2001, 2002, 2003  CP2K developers group               !
!>    Copyright (C) 2004, 2005, 2006, 2007  CP2K developers group               !
!>    Copyright (C) 2008, 2009, 2010        CP2K developers group               !
!> \par
!>    This program is free software; you can redistribute it and/or modify      !
!>    it under the terms of the GNU General Public License as published by      !
!>    the Free Software Foundation; either version 2 of the License, or         !
!>    (at your option) any later version.                                       !
!> \par
!>    This program is distributed in the hope that it will be useful,           !
!>    but WITHOUT ANY WARRANTY; without even the implied warranty of            !
!>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             !
!>    GNU General Public License for more details.                              !
!> \par
!>    You should have received a copy of the GNU General Public License         !
!>    along with this program; if not, write to the Free Software               !
!>    Foundation, Inc., 51 Franklin Street, Fifth Floor,                        !
!>    Boston, MA  02110-1301, USA.                                              !
!> \par
!>    See also http://www.fsf.org/licensing/licenses/gpl.html                   !
!> \par
!>    CP2K, including its sources and pointers to the authors                   !
!>    can be found at  http://cp2k.berlios.de/                                  !
!> \note
!>       should be kept as lean as possible.
!>       see cp2k_run for more comments  
!> \author Joost VandeVondele !
! *****************************************************************************
PROGRAM cp2k
  USE cp2k_info,         ONLY : print_cp2k_license,&
                                cp2k_version,&
                                id_cp2k_version, get_runtime_info
  USE cp2k_runs,         ONLY : write_cp2k_html_manual, &
                                write_xml_file, run_input
  USE input_constants,   ONLY : id_development_version
  USE machine,           ONLY : m_iargc,m_getarg, default_output_unit
  USE message_passing,   ONLY : mp_sync, mp_abort
  USE kinds,             ONLY : default_path_length,&
                                default_string_length
  USE f77_interface,     ONLY : init_cp2k, finalize_cp2k, &
                                check_input, default_para_env
  USE f77_blas
#include "cp_common_uses.h"
  IMPLICIT NONE
#if defined(__USE_PAT)
#include "pat_apif.h"
#endif

  CHARACTER(LEN=default_path_length)   :: input_file_name,output_file_name,arg_att
  CHARACTER(LEN=default_string_length) :: dev_flag
  INTEGER                              :: output_unit,istat
  INTEGER                              :: ierr, n_arg, i_arg
  LOGICAL                              :: check, manual, usage, echo_input, command_line_error,&
                                          permissive
  LOGICAL                              :: run_it, force_run, has_input, xml, print_version, print_license
  TYPE(cp_error_type)                  :: error 

  ! output goes to the screen by default
  output_unit=default_output_unit

  ! get runtime information
  CALL get_runtime_info

  ! initialize with a default error (a proper error will be initialized later)
  CALL cp_error_init(error)

  ! do the initial setup of cp2k, provides a default_para_env
  CALL init_cp2k(init_mpi=.TRUE.,ierr=ierr)

#if defined(__USE_PAT)
  CALL PAT_sampling_state(PAT_STATE_OFF,istat)
#endif

  IF (ierr==0) THEN
#if defined(__HMD)
     DO i=0,default_para_env%num_pe-1
        IF (default_para_env%mepos==i) THEN
           PRINT *, "***** proc ",i
           CALL _dump_allocated(%VAL(0))
        END IF
        CALL mp_sync(default_para_env%group)
     END DO
#endif

     ! set default behaviour for the command line switches
     check      =.FALSE.
     manual     =.FALSE.
     usage      =.FALSE.
     echo_input =.FALSE.
     has_input  =.FALSE.
     run_it     =.TRUE.
     force_run  =.FALSE.
     permissive =.FALSE.
     print_version =.FALSE.
     print_license =.FALSE.
     command_line_error=.FALSE.
     xml = .FALSE.
     input_file_name  = "Missing input file name"   ! no default
     output_file_name = "__STD_OUT__"               ! by default we go to std_out

     ! parsing of command line arguments
     n_arg=m_iargc()

     ! this is an error condition, provide usage information
     IF (n_arg<1) THEN
         WRITE (output_unit,"(a)") ""
         WRITE (output_unit,"(a)") "ERROR: At least one command line argument must be specified"
         command_line_error=.TRUE.
     ENDIF

     i_arg=0
arg_loop: DO WHILE (i_arg<n_arg)
        i_arg=i_arg+1
        CALL m_getarg(i_arg,arg_att)
        SELECT CASE(arg_att)
        CASE("--permissive-echo")
           permissive=.TRUE.
        CASE("--check","-c")
           check=.TRUE.
           run_it=.FALSE.
           echo_input=.FALSE.
        CASE("--echo","-e")
           check=.TRUE.
           run_it=.FALSE.
           echo_input=.TRUE.
        CASE("--html-manual")
           manual=.TRUE.
           run_it=.FALSE.
        CASE("-v","--version")
           print_version=.TRUE.
           run_it=.FALSE.
        CASE("--license")
           print_license=.TRUE.
           run_it=.FALSE.
        CASE("--run","-r")
           force_run=.TRUE.
        CASE("-help","--help","-h")
           usage=.TRUE.
           run_it=.FALSE.
        CASE("-i")
           i_arg=i_arg+1
           CALL m_getarg(i_arg,arg_att)
           ! argument does not start with a - it is an filename
           IF (.NOT. arg_att(1:1)=="-") THEN
              input_file_name=arg_att
              has_input=.TRUE.
           ELSE
              WRITE (output_unit,"(a)") ""
              WRITE (output_unit,"(a)") "ERROR: the input file name"//TRIM(arg_att)//" starts with -"
              command_line_error=.TRUE.
              EXIT arg_loop
           ENDIF
        CASE("-o")
           i_arg=i_arg+1
           CALL m_getarg(i_arg,arg_att)
           ! argument does not start with a - it is an filename
           IF (.NOT. arg_att(1:1)=="-") THEN
              output_file_name=arg_att
           ELSE
              WRITE (output_unit,"(a)") ""
              WRITE (output_unit,"(a)") "ERROR: the output file name"//TRIM(arg_att)//" starts with -"
              command_line_error=.TRUE.
              EXIT arg_loop
           ENDIF
        CASE("--xml")
           xml = .TRUE.
           run_it = .FALSE.
        CASE default
           ! if the last argument does not start with a - it is an input filename
           !MK in order to digest the additional flags of mpirun
           IF ((.NOT.has_input).AND.&
               (i_arg == n_arg).AND.&
               (.NOT.arg_att(1:1) == "-")) THEN
              input_file_name=arg_att
              has_input=.TRUE.
           END IF
        END SELECT
     END DO arg_loop
     IF (permissive) THEN
        run_it=.FALSE.
        echo_input=.TRUE.
        check=.TRUE.
     END IF
     IF ((run_it.OR.force_run.OR.check.OR.echo_input).AND.&
         (.NOT.has_input).AND.(.NOT.command_line_error)) THEN
       WRITE (UNIT=output_unit,FMT="(/,A)")&
         "ERROR: an input file name is required"
       command_line_error = .TRUE.
     END IF

     ! some first info concerning how to run CP2K 

     IF (usage.OR.command_line_error) THEN
       IF (default_para_env%source == default_para_env%mepos) THEN
         CALL m_getarg(0,arg_att) ! actually this is non-standard
         WRITE (UNIT=output_unit,FMT="(/,A,4(/,/,T5,A),/,(T5,A))")&
           TRIM(arg_att)//" [-c|--check] [-e|--echo] [-h|--help]"//&
                          " [--html-manual] [--xml] [-r|-run]"//&
                          " [-o ouput_file] [--permissive-echo]"//&
                          " [[-i] input_file]",&
           "Starts the CP2K program, see <http://cp2k.berlios.de/>",&
           "The easiest way is "//TRIM(arg_att)//" input_file",&
           "The following options can be used:",&
           "-o output_file : provides an output file name [default: screen]",&
           "-i input_file  : provides an input file name, if it is the last",&
           "                 argument, the -i is not needed"
         WRITE (UNIT=output_unit,FMT="(/,T5,A,/,/,(T5,A))")&
           "These switches skip the simulation, unless [-r|-run] is specified:",&
           "--run, -r         : force a run after any of the following",&
           "--check, -c       : perform a syntax check of [input_file]",&
           "--echo, -e        : echos [input_file], and makes all defaults explicit",&
           "                    The input is also checked, but only failure is reported",&
           "--version,-v      : print the CP2K version number",&
           "--license         : print the CP2K license",&
           "--help, -h        : writes this message",&
           "--html-manual     : writes a HTML reference manual of the CP2K input file",&
           "                    in the current directory. index.html is a good starting",&
           "                    point for browsing",&
           "--permissive-echo : ignores unknown keywords and sections in the input",&
           "                    and echoes it, you cannot run, you have to run the",&
           "                    dumped input",&
           "--xml             : dumps the whole CP2K input structure as a XML file",&
           "                    xml2htm generates a HTML manual from this XML file"
       END IF
     END IF

     IF (.NOT. command_line_error) THEN

        ! write the version string
        IF (print_version) THEN
           IF (default_para_env%source==default_para_env%mepos) THEN
              dev_flag = ""
              IF (id_cp2k_version==id_development_version) dev_flag=" (Development Version)"
              WRITE(output_unit,'(T2,A)') cp2k_version//TRIM(dev_flag)
           ENDIF
        ENDIF

        ! write the license
        IF (print_license) THEN
           IF (default_para_env%source==default_para_env%mepos) THEN
              CALL print_cp2k_license(output_unit)
           ENDIF
        ENDIF

        ! if a manual is needed
        IF (manual) THEN
           IF (default_para_env%mepos==default_para_env%source) THEN
              CALL write_cp2k_html_manual(error=error)
           END IF
        ENDIF

        IF (xml) THEN
          IF (default_para_env%mepos == default_para_env%source) THEN
            CALL write_xml_file(error=error)
          END IF
        END IF

        IF (check) THEN
           CALL check_input(input_file_name,output_file_name,echo_input=echo_input,&
                permissive=permissive,ierr=ierr)
           IF (default_para_env%source==default_para_env%mepos) THEN
              IF (ierr==0) THEN
                IF (.NOT. echo_input) THEN
                   WRITE (output_unit,"(a)") "SUCCESS, the input could be parsed correctly."
                   WRITE (output_unit,"(a)") "         This is does not guarantee that this input is meaningful"
                   WRITE (output_unit,"(a)") "         or will run successfully"
                ENDIF
              ELSE
                 WRITE (output_unit,"(a)") "ERROR,   the input could *NOT* be parsed correctly."
                 WRITE (output_unit,"(a)") "         Please, check and correct it"
              END IF
           END IF
        ENDIF

        IF (run_it .OR. force_run) THEN
           CALL run_input(input_file_name,output_file_name,ierr=ierr)
           IF (ierr.NE.0) THEN
              CALL mp_abort()
           ENDIF
        END IF
     ENDIF

#if defined(__HMD)
     DO i=0,default_para_env%num_pe-1
        IF (default_para_env%mepos==i) THEN
           PRINT *, "***** proc ",i
           CALL _dump_allocated_delta(%VAL(0))
        END IF
        CALL mp_sync(default_para_env%group)
     END DO
#endif
  ELSE
     WRITE(UNIT=output_unit,FMT="(/,A)") "initial setup (MPI ?) error" 
  END IF

#if defined(__USE_PAT)
  CALL PAT_sampling_state(PAT_STATE_ON,istat)
#endif

  ! and the final cleanup
  CALL finalize_cp2k(finalize_mpi=.TRUE.,ierr=ierr)

END PROGRAM
