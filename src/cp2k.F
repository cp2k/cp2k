!!****** cp2k/cp2k [1.0] *
!!
!!   NAME
!!     cp2k, The main program
!!
!!   COPYRIGHT
!!I---------------------------------------------------------------------------I
!!I                                                                           I
!!I CP2K: A general program to perform molecular dynamics simulations         I
!!I Copyright (C) 2000,2001,2002,2003,2004  CP2K developers group             I
!!I                                                                           I
!!I This program is free software; you can redistribute it and/or modify      I
!!I it under the terms of the GNU General Public License as published by      I
!!I the Free Software Foundation; either version 2 of the License, or         I
!!I (at your option) any later version.                                       I
!!I                                                                           I
!!I This program is distributed in the hope that it will be useful,           I
!!I but WITHOUT ANY WARRANTY; without even the implied warranty of            I
!!I MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             I
!!I GNU General Public License for more details.                              I
!!I                                                                           I
!!I You should have received a copy of the GNU General Public License         I
!!I along with this program; if not, write to the Free Software               I
!!I Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   FUNCTION
!!      starts the program
!!
!!   AUTHOR
!!      Joost VandeVondele 
!!
!!   NOTES
!!      should be kept as lean as possible.
!!      see cp2k_run for more comments
!!
!------------------------------------------------------------------------------
PROGRAM cp2k
  USE cp2k_runs,         ONLY : write_cp2k_html_manual
  USE machine,           ONLY : m_iargc,m_getarg
  USE message_passing,   ONLY : mp_world_init, &
       mp_world_finalize, mp_sync
  USE kinds,             ONLY : default_string_length
  USE cp_para_types, only: cp_para_env_type
  USE cp_para_env, ONLY: cp_para_env_create, cp_para_env_release
  USE f77_interface, ONLY: init_cp2k, finalize_cp2k, run_input, &
       check_input, default_para_env
  IMPLICIT NONE

  CHARACTER(LEN=default_string_length) :: input_file_name,arg_att
  INTEGER                              :: output_unit
  INTEGER                              :: ierr, i, n_arg
  LOGICAL                              :: check, manual, usage, echo_input
  ! this is the only place that should need unit 6 all others will get 
  ! something reasonable from the cp_logger
  INTEGER, PARAMETER                   :: default_output_unit = 6 
  TYPE(cp_para_env_type), pointer      :: para_env

  ! output goes to the screen by default
  output_unit=default_output_unit

  ! do the initial setup of cp2k
  CALL init_cp2k(init_mpi=.TRUE.,ierr=ierr)

  IF (ierr==0) then
     para_env => default_para_env
#if defined(__HMD)
    DO i=0,para_env%num_pe-1
       IF (para_env%mepos==i) THEN
          PRINT *, "***** proc ",i
          CALL _dump_allocated(%VAL(0))
       END IF
       CALL mp_sync(para_env%group)
    END DO
#endif
     ! get the input file name and command line flags
     check=.FALSE.
     manual=.false.
     usage=.false.
     echo_input=.false.
     input_file_name = "input"
     n_arg=0
     DO i=1,m_iargc()
        CALL m_getarg(i,arg_att)
        SELECT CASE(arg_att)
        CASE("-check","--check-input")
           check=.TRUE.
           echo_input=.false.
        CASE("-CHECK","--CHECK-INPUT")
           check=.TRUE.
           echo_input=.true.
        CASE("--html-manual")
           manual=.true.
        CASE("-help","--help","-h")
           usage=.true.
        CASE default
           n_arg=n_arg+1
           IF (n_arg==1) THEN
              input_file_name=arg_att
           END IF
        END SELECT
     END DO
     !     IF (n_arg>1) usage=.true. ! check ?
     IF (manual) THEN
        IF (para_env%mepos==para_env%source) THEN
           CALL write_cp2k_html_manual()
        END IF
     ELSE IF (usage) THEN
        IF (para_env%source==para_env%mepos) THEN
           CALL m_getarg(0,arg_att)
           WRITE (output_unit,"(a)") trim(arg_att)
           WRITE (output_unit,"('           [--html-manual] [-check|-CHECK|--check-input|--CHECK-INPUT]')")
           WRITE (output_unit,"('           [-h|-help|--help] input_file',/)")
           WRITE (output_unit,"('runs the cp2k simulation program see <http://cp2k.berlios.de>.')")
           WRITE (output_unit,"(' --html-manual',/,"//&
                "'   writes a html manual of the input in the current directory')")
           WRITE (output_unit,"(' -check -CHECK --check-input --CHECK-INPUT',/,"//&
                "'   just checks the input without running the program',/"//&
                "'   The upcase versions write the input with all the defaults made explicit')")
           WRITE (output_unit,"(' -h -help --help',/,"//&
                "'   writes out this message')")
        END IF
     ELSE IF (check) then
        CALL check_input(input_file_name,"__STD_OUT__",echo_input=echo_input,&
             ierr=ierr)
        IF (para_env%source==para_env%mepos) THEN
           IF (ierr==0) THEN
              WRITE (output_unit,&
                   "(/,'SUCCESS, the input could be interpreted correctly')")
              WRITE (output_unit,&
                   "('         This does not imply the the the input is meaningful or guarantees')")
              WRITE (output_unit,&
                   "('         a succesful run')")
           ELSE
              WRITE (output_unit,&
                   "(/,'ERROR,   the input could NOT be interpreted correctly')")
              WRITE (output_unit,&
                   "('         please check and correct it')")
           END IF
        END IF
     ELSE
        CALL run_input(input_file_name,"__STD_OUT__",ierr=ierr)
     END IF
  ELSE
     WRITE(output_unit,"(/,'initial setup (mpi) error')")
  ENDIF
 #if defined(__HMD)
    DO i=0,para_env%num_pe-1
       IF (para_env%mepos==i) THEN
          PRINT *, "***** proc ",i
          CALL _dump_allocated_delta(%VAL(0))
       END IF
       CALL mp_sync(para_env%group)
    END DO
#endif
 
  ! and the final cleanup
  CALL finalize_cp2k(finalize_mpi=.TRUE.,ierr=ierr)

END PROGRAM
