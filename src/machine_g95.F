!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/machine_pgi [1.0] *
!!
!!   NAME
!!     machine_g95
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     APSI & JGH
!!
!!   MODIFICATION HISTORY
!!     - m_flush added (12.06.2002,MK)
!!     - print_memory changed (24.09.2002,MK)
!!     - adapted for g95 (29.05.2003,JVdV)
!!
!!   SOURCE
!!*****
!******************************************************************************

MODULE machine_g95

  USE kinds,                           ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: m_walltime, m_cputime, m_datum, m_flush, print_memory, &
            m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
            m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c

!!*****
!******************************************************************************

CONTAINS

FUNCTION m_loc_r(a) RESULT(res)
  REAL(KIND=dp), INTENT(in), DIMENSION(*) :: a
  integer :: res
  res=-1
END FUNCTION m_loc_r
  
FUNCTION m_loc_c(a) result(res)
  COMPLEX(KIND=dp), INTENT(in), DIMENSION(*) :: a
  integer :: res
  res=-1
END FUNCTION m_loc_c

! can be used to get a nice core
SUBROUTINE m_abort()
   CALL abort()
END SUBROUTINE m_abort

! the number of arguments of the fortran program
FUNCTION m_iargc() RESULT (ic)
    INTEGER                                  :: ic
    INTEGER                                  :: iargc
    ic = iargc()
END FUNCTION m_iargc

!******************************************************************************
!!  cpu time in seconds
FUNCTION m_cputime() RESULT (ct)
    REAL(KIND=dp)                                :: ct
    REAL                                     :: wtl
    CALL CPU_TIME(wtl)
    ct=wtl
END FUNCTION m_cputime

!******************************************************************************
!!  wall clock time in seconds
FUNCTION m_walltime() RESULT (wt)
    REAL(KIND=dp)                                :: wt

    INTEGER                                  :: count, count_rate
    CALL system_clock(count=count,count_rate=count_rate)
    wt = REAL(count,KIND=dp)/REAL(count_rate,KIND=dp)
END FUNCTION m_walltime

!******************************************************************************
!!  time and date
SUBROUTINE m_datum(cal_date)
    CHARACTER(len=*), INTENT(OUT)            :: cal_date

    character(len=10)                        :: time
    character(len=5)                         :: zone
    character(len=8)                         :: date
    integer, dimension(8)                    :: values
    call date_and_time(date=date, time=time, zone=zone, values=values)
    cal_date=date//" "//time
END SUBROUTINE m_datum

! *****************************************************************************
! flush a given unit
  SUBROUTINE m_flush(lunit)
    INTEGER, INTENT(IN)                      :: lunit
    CALL flush(lunit)
  END SUBROUTINE m_flush

! *****************************************************************************

  SUBROUTINE print_memory(checkpoint)

!   Purpose: Print the memory usage of the program at checkpoint.

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: checkpoint
    CHARACTER(LEN=200)                       :: cmdstr
    INTEGER                                  :: pid

!   ---------------------------------------------------------------------------

    CALL m_getpid(pid)

    WRITE(UNIT=cmdstr,FMT="(A,I8,A,I8,A)")&
      "echo -e ""\\n Memory size of process""",pid,""" <"//TRIM(checkpoint)//&
      ">:"" `ps -o vsz -p",pid," | tail -1` KBytes""\\c"""

    CALL system(cmdstr)

  END SUBROUTINE print_memory

! *****************************************************************************
SUBROUTINE m_hostnm(hname)
    CHARACTER(len=*), INTENT(OUT)            :: hname
    INTEGER :: hostnm,ierror
    ierror=hostnm(hname)
END SUBROUTINE m_hostnm

!******************************************************************************
SUBROUTINE m_getcwd(curdir)
    CHARACTER(len=*), INTENT(OUT)            :: curdir
    INTEGER                                  :: getcwd,ierror
    ierror = getcwd(curdir)
END SUBROUTINE m_getcwd

!******************************************************************************
SUBROUTINE m_chdir(dir,ierror)
    CHARACTER(len=*), INTENT(IN)          :: dir
    INTEGER,          INTENT(OUT)            :: ierror
    INTEGER                                  :: chdir
    ierror = chdir(dir)
END SUBROUTINE m_chdir

!******************************************************************************
SUBROUTINE m_getlog(user)
    CHARACTER(len=*), INTENT(OUT)            :: user
    CALL getlog(user)
END SUBROUTINE m_getlog

!******************************************************************************
SUBROUTINE m_getuid(uid)
    INTEGER, INTENT(OUT)                     :: uid
    INTEGER :: getuid
    uid = getuid()
END SUBROUTINE m_getuid

!******************************************************************************
SUBROUTINE m_getpid(pid)
    INTEGER, INTENT(OUT)                     :: pid
    INTEGER :: getpid
    pid = getpid()
END SUBROUTINE m_getpid

!******************************************************************************
SUBROUTINE m_getarg(i,arg)
    INTEGER, INTENT(IN)                      :: i
    CHARACTER(len=*), INTENT(OUT)            :: arg
    CALL getarg(i,arg)
END SUBROUTINE m_getarg

END MODULE machine_g95
