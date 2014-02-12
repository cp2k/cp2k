!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      - m_flush added (12.06.2002,MK)
!>      - print_memory changed (24.09.2002,MK)
!>      - m_memory implemented (28.04.2011,IB)
!> \author APSI & JGH
! *****************************************************************************
  USE ISO_C_BINDING
  USE f77_blas
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: m_cputime, m_flush, m_memory, &
            m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
            m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c,m_mov, m_memory_details, &
            m_procrun

  INTERFACE m_loc_r
     MODULE PROCEDURE m_loc_r1,m_loc_r2,m_loc_r3,m_loc_r4
  END INTERFACE

  INTERFACE m_loc_c
     MODULE PROCEDURE m_loc_c1,m_loc_c2,m_loc_c3,m_loc_c4
  END INTERFACE

  TYPE, BIND(C) :: timeval
     INTEGER(C_LONG) :: tv_sec, tv_usec
  END TYPE

  TYPE, BIND(C) :: rusage
     TYPE(timeval) :: ru_utime, ru_stime
     INTEGER(C_LONG) :: ru_maxrss, ru_ixrss, ru_idrss, ru_isrss, ru_minflt, ru_majflt, &
                     ru_nswap, ru_inblock, ru_oublock, ru_msgsnd, ru_msgrcv, &
                     ru_nsignals, ru_nvcsw, ru_nivcsw
  END TYPE

  INTERFACE
    INTEGER(C_INT) FUNCTION getrusage (who, rusage) BIND(C, name='getrusage')
    USE ISO_C_BINDING
      INTEGER(C_INT), VALUE :: who
      TYPE(C_PTR),    VALUE :: rusage
    END FUNCTION getrusage
  END INTERFACE

CONTAINS
! *** get more detailed memory info, all units are bytes.
! *** the only 'useful' option is MemLikelyFree which is an estimate of remaining memory
! *** assumed to give info like /proc/meminfo while MeMLikelyFree is the amount of
! *** memory we're likely to be able to allocate, but not necessarily in one chunk
! *** zero means not available
  SUBROUTINE m_memory_details(MemTotal,MemFree,Buffers,Cached,Slab,SReclaimable,MemLikelyFree)

    INTEGER(kind=int_8), OPTIONAL            :: MemTotal, MemFree, Buffers, &
                                                Cached, Slab, SReclaimable, &
                                                MemLikelyFree

     MemTotal=0
     MemFree=0
     Buffers=0
     Cached=0
     Slab=0
     SReclaimable=0
     MemLikelyFree=0

  END SUBROUTINE m_memory_details

! returns if a process is running on the local machine
! 1 if yes and 0 if not
INTEGER FUNCTION m_procrun(id) RESULT (run_on)
    INTEGER           ::   id, ios
    CHARACTER(len=80) ::   filename, tmp
    CHARACTER(len=8)  ::   id_s

    run_on = 0

END FUNCTION m_procrun



! *****************************************************************************
FUNCTION m_loc_r1(a) RESULT(res)
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_r1

! *****************************************************************************
FUNCTION m_loc_r2(a) RESULT(res)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(in)                             :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_r2

! *****************************************************************************
FUNCTION m_loc_r3(a) RESULT(res)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(in)                             :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_r3

! *****************************************************************************
FUNCTION m_loc_r4(a) RESULT(res)
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(in)                             :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_r4

! *****************************************************************************
FUNCTION m_loc_c1(a) RESULT(res)
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(in)                             :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_c1

! *****************************************************************************
FUNCTION m_loc_c2(a) RESULT(res)
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(in)                             :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_c2

! *****************************************************************************
FUNCTION m_loc_c3(a) RESULT(res)
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(in)                             :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_c3

! *****************************************************************************
FUNCTION m_loc_c4(a) RESULT(res)
    COMPLEX(KIND=dp), &
      DIMENSION(:, :, :, :), INTENT(in)      :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_c4

! can be used to get a nice core
! *****************************************************************************
SUBROUTINE m_abort()
   CALL abort()
END SUBROUTINE m_abort

! *****************************************************************************
FUNCTION m_iargc() RESULT (ic)
    INTEGER                                  :: ic

    INTEGER, EXTERNAL                        :: iargc

  ic = iargc()
END FUNCTION m_iargc

!!  cpu time in seconds
! *****************************************************************************
FUNCTION m_cputime() RESULT (ct)
    REAL(KIND=dp)                            :: ct

    REAL                                     :: wtl

  CALL CPU_TIME(wtl)
  ct=wtl
END FUNCTION m_cputime

! *****************************************************************************
!> \brief   Flush the output to a logical unit.
!> \author  MK
!> \date    14.10.1999
!> \version 1.0
! *****************************************************************************
  SUBROUTINE m_flush(lunit)

    INTEGER, INTENT(IN)                      :: lunit

    CALL flush_(lunit)

  END SUBROUTINE m_flush

! returns the total amount of memory [bytes] in use, if known, zero otherwise
! *****************************************************************************
  FUNCTION m_memory()
    INTEGER(KIND=int_8)                      :: m_memory

    INTEGER(C_INT)                           :: ret
    TYPE(rusage), TARGET                     :: usage

    ret = getrusage(0, C_LOC(usage))
    m_memory = usage%ru_maxrss * 1024
  END FUNCTION m_memory

! *****************************************************************************
  SUBROUTINE m_mov(source,TARGET)

    CHARACTER(LEN=*), INTENT(IN)             :: source, TARGET

    CHARACTER(LEN=LEN(source)+1)             :: ls1
    CHARACTER(LEN=LEN(TARGET)+1)             :: ls2
    INTEGER                                  :: pos, stat
    INTEGER, EXTERNAL                        :: rename

! **** system call caused problems on infiniband / Linux setup on blanc.cscs.ch
! **** due to the fork.
! **** rename is a better solution
! CHARACTER(LEN=2*default_string_length+4) :: cmd
! cmd = "mv " // source(1:LEN_TRIM(source)) // " " // TARGET(1:LEN_TRIM(TARGET))
! CALL system(cmd)

    ls1=source
    ls2=TARGET
    pos=LEN_TRIM(source)+1
    ls1(pos:)=CHAR(0)
    pos=LEN_TRIM(TARGET)+1
    ls2(pos:)=CHAR(0)
    stat=rename(ls1,ls2)
    IF (stat .NE. 0) THEN
      WRITE(6,*) "Trying to move "//TRIM(source)//" to "//TRIM(TARGET)//"."
      WRITE(6,*) "rename returned status: ",stat
      STOP "Problem moving file"
    ENDIF

  END SUBROUTINE m_mov

! *****************************************************************************
SUBROUTINE m_hostnm(hname)
    CHARACTER(len=*), INTENT(OUT)            :: hname

    INTEGER                                  :: hostnm_, ierror

  ierror = hostnm_(hname)
END SUBROUTINE m_hostnm
! *****************************************************************************
SUBROUTINE m_getcwd(curdir)
    CHARACTER(len=*), INTENT(OUT)            :: curdir

    INTEGER                                  :: getcwd_, ierror

  ierror = getcwd_(curdir)
END SUBROUTINE m_getcwd

! *****************************************************************************
SUBROUTINE m_chdir(dir,ierror)
    CHARACTER(len=*), INTENT(IN)             :: dir
    INTEGER, INTENT(OUT)                     :: ierror

    CHARACTER(LEN=1000)                      :: dir_local
    INTEGER                                  :: chdir

    dir_local=dir
    dir_local(LEN_TRIM(dir)+1:LEN_TRIM(dir)+1)=ACHAR(0)
    ierror = chdir(dir_local)
END SUBROUTINE m_chdir

! *****************************************************************************
SUBROUTINE m_getlog(user)
    CHARACTER(len=*), INTENT(OUT)            :: user

  CALL getlog_(user)
END SUBROUTINE m_getlog
! *****************************************************************************
SUBROUTINE m_getuid(uid)
    INTEGER, INTENT(OUT)                     :: uid

    INTEGER                                  :: getuid_

  uid = getuid_()
END SUBROUTINE m_getuid
! *****************************************************************************
SUBROUTINE m_getpid(pid)
    INTEGER, INTENT(OUT)                     :: pid

    INTEGER                                  :: getpid_

  pid = getpid_()
END SUBROUTINE m_getpid
! *****************************************************************************
SUBROUTINE m_getarg(i,arg)
    INTEGER, INTENT(IN)                      :: i
    CHARACTER(len=*), INTENT(OUT)            :: arg

  CALL getarg(i,arg)
END SUBROUTINE m_getarg
