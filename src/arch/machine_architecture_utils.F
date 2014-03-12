!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group group!
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Machine Hardware support for cp2k
!           Pack and unpack for machine representation
!> \author  Christiane Pousa Ribeiro
!> \date    2011-08-08
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-08-23
! *****************************************************************************
MODULE machine_architecture_utils
#if !defined (__HAS_NO_ISO_C_BINDING)
  USE ISO_C_BINDING
#endif
  USE dummy_c_bindings

 USE machine_architecture_types
 USE ma_kinds

 IMPLICIT NONE

 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'machine_architecture_utils'

 PUBLIC :: ma_pack_threads, ma_unpack_threads
 PUBLIC :: ma_create_process, ma_set_threads

 PUBLIC :: ascii_to_string, string_to_ascii
 PUBLIC :: integer_to_string

 CONTAINS

! *****************************************************************************
!> \brief   Create the information structure of a process
!> \var     proc_info is the processor structure
!> \var     mp_env is the comm information for MPI
!> \var     nr_threads is the number of threads of a processor
! *****************************************************************************
 SUBROUTINE ma_create_process(proc_info, mp_env, nr_threads, core, node)
    TYPE(ma_process), INTENT(OUT)            :: proc_info
    TYPE(ma_mp_type), INTENT(IN)             :: mp_env
    INTEGER, INTENT(IN)                      :: nr_threads, core, node

 proc_info%mp_info%myproc = mp_env%myproc
 proc_info%mp_info%numproc = mp_env%numproc
 proc_info%mp_info%mp_group = mp_env%mp_group
 proc_info%nr_threads = nr_threads
 proc_info%core = core
 proc_info%node = node
 proc_info%mp_info%myid = mp_env%myid
 proc_info%mp_info%all_proc = mp_env%all_proc
 ALLOCATE (proc_info%threads_info(nr_threads))

 END SUBROUTINE ma_create_process

! *****************************************************************************
!> \brief   Set the information of the group of threads of a process
!> \var     proc_info is the processor structure
!> \var     thread_info is the thread information
!> \var     id is the OMP id of the thread
! *****************************************************************************
 SUBROUTINE ma_set_threads(proc_info, thread_info, id)
    TYPE(ma_process)                         :: proc_info
    TYPE(thread_inf), INTENT(IN)             :: thread_info
    INTEGER                                  :: id

 proc_info%threads_info(id)%core = thread_info%core
 proc_info%threads_info(id)%node = thread_info%node
 proc_info%threads_info(id)%id_omp = thread_info%id_omp
 proc_info%threads_info(id)%id_real = thread_info%id_real

 END SUBROUTINE ma_set_threads

! *****************************************************************************
!> \brief   Pack components of a machine architecture into a single vector
!  thread is a vector that must be allocated first
! *****************************************************************************
 SUBROUTINE ma_pack_threads(proc_info, thread_send)
    TYPE(ma_process), INTENT(IN)             :: proc_info
    INTEGER, DIMENSION(:, :), INTENT(OUT), &
      OPTIONAL                               :: thread_send

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ma_pack_threads', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: it

 DO it=1, proc_info%nr_threads
  thread_send(1,it) =  proc_info%threads_info(it)%core
  thread_send(2,it) =  proc_info%threads_info(it)%node
  thread_send(3,it) =  proc_info%threads_info(it)%id_omp
  thread_send(4,it) =  proc_info%threads_info(it)%id_real
 END DO

 END SUBROUTINE ma_pack_threads

! *****************************************************************************
!> \brief   Unpack components of a machine architecture into the thread
!>          structre
! *****************************************************************************
 SUBROUTINE ma_unpack_threads(proc_info, thread_recv, thread_map)
    TYPE(ma_process), INTENT(IN)             :: proc_info
    INTEGER, DIMENSION(:, :, :), &
      INTENT(IN), OPTIONAL                   :: thread_recv
    TYPE(thread_inf), DIMENSION(:, :)        :: thread_map

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ma_unpack_threads', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: from, it

 DO from=1, proc_info%mp_info%numproc
  DO it=1, proc_info%nr_threads
   thread_map(from,it)%core =  thread_recv(1,it,from)
   thread_map(from,it)%node =  thread_recv(2,it,from)
   thread_map(from,it)%id_omp =  thread_recv(3,it,from)
   thread_map(from,it)%id_real =  thread_recv(4,it,from)
  END DO
 END DO
 END SUBROUTINE ma_unpack_threads

 ! *****************************************************************************
 ! COPY from string_utils - Convert a sequence of integer numbers (ASCII code)
 ! to a string.
 ! *****************************************************************************
   SUBROUTINE ascii_to_string(nascii,string)

    INTEGER, DIMENSION(:), INTENT(IN)        :: nascii
    CHARACTER(LEN=*), INTENT(OUT)            :: string

    INTEGER                                  :: i

   string = ""

    DO i=1,MIN(LEN(string),SIZE(nascii))
      IF ((nascii(i) >= 0).AND.(nascii(i) <= 127)) THEN
             string(i:i) = CHAR(nascii(i))
      ELSE
          string(i:i) = " "
      END IF
    END DO

END SUBROUTINE ascii_to_string

! *****************************************************************************
! COPY from string_utils - Convert a string to sequence of integer numbers.
! *****************************************************************************
  SUBROUTINE string_to_ascii(string,nascii)

    CHARACTER(LEN=*), INTENT(IN)             :: string
    INTEGER, DIMENSION(:), INTENT(OUT)       :: nascii

    INTEGER                                  :: i

  nascii(:) = 0

  DO i=1,MIN(LEN(string),SIZE(nascii))
         nascii(i) = ICHAR(string(i:i))
  END DO

END SUBROUTINE string_to_ascii

! *****************************************************************************
! COPY from string_utils - Convert a integer to a string.
! *****************************************************************************
SUBROUTINE integer_to_string(inumber,string)
    INTEGER, INTENT(IN)                      :: inumber
    CHARACTER(LEN=*), INTENT(OUT)            :: string

 WRITE (UNIT=string,FMT='(I0)') inumber

END SUBROUTINE integer_to_string

END MODULE machine_architecture_utils
