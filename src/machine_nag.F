!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/machine_nag [1.0] *
!!
!!   NAME
!!     machine_nag
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!     based on guessing calls from the .mod files. No manual yet.
!!
!!   SOURCE
!!*****
!******************************************************************************

MODULE machine_nag

  USE f90_unix
  USE f90_unix_dir
  USE f90_unix_env,                    ONLY: gethostname,&
                                             getlogin
  USE kinds,                           ONLY: dbl

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: m_walltime, m_cputime, m_datum, m_flush, print_memory, &
            m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, m_iargc, m_abort

!!*****
!******************************************************************************

CONTAINS


! can be used to get a nice core
SUBROUTINE m_abort()
   CALL abort()
END SUBROUTINE m_abort


FUNCTION m_iargc() RESULT (ic)
    INTEGER                                  :: ic

  ic = iargc()
END FUNCTION m_iargc

!******************************************************************************
!!  cpu time in seconds
FUNCTION m_cputime() RESULT (ct)
    REAL(dbl)                                :: ct

    REAL                                     :: ctr

  CALL CPU_TIME(ctr)
  ct=ctr
END FUNCTION m_cputime
!******************************************************************************
!!  wall clock time in seconds
FUNCTION m_walltime() RESULT (wt)
    REAL(dbl)                                :: wt

    INTEGER                                  :: count, count_rate

  CALL system_clock(count=count,count_rate=count_rate)
  wt = DBLE(count)/count_rate
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

  SUBROUTINE m_flush(lunit)

!   Purpose: Flush the output to a logical unit.

!   History: - Creation (14.10.1999,MK)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: lunit

!   ---------------------------------------------------------------------------

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

    INTEGER                                  :: len
    CHARACTER(len=60)                        :: tmpname

  len=60
  CALL gethostname(tmpname,len)
  hname=tmpname(1:len)
END SUBROUTINE m_hostnm
!******************************************************************************
SUBROUTINE m_getcwd(curdir)
    CHARACTER(len=*), INTENT(OUT)            :: curdir

    INTEGER                                  :: len
    CHARACTER(len=60)                        :: tmpdir

  len=60
  CALL getcwd(tmpdir,len)
  curdir=tmpdir(1:len)
END SUBROUTINE m_getcwd
!******************************************************************************
SUBROUTINE m_getlog(user)
    CHARACTER(len=*), INTENT(OUT)            :: user

    INTEGER                                  :: len
    CHARACTER(len=60)                        :: tmp

  len=60
  CALL getlogin(tmp,len)
  user=tmp(1:len)
END SUBROUTINE m_getlog
!******************************************************************************
SUBROUTINE m_getuid(uid)
    INTEGER, INTENT(OUT)                     :: uid

  uid = getuid()
END SUBROUTINE m_getuid
!******************************************************************************
SUBROUTINE m_getpid(pid)
    INTEGER, INTENT(OUT)                     :: pid

  pid = getpid()
END SUBROUTINE m_getpid
!******************************************************************************
SUBROUTINE m_getarg(i,arg)
    INTEGER, INTENT(IN)                      :: i
    CHARACTER(len=*), INTENT(OUT)            :: arg

  CALL getarg(i,arg)
END SUBROUTINE m_getarg
!******************************************************************************
END MODULE machine_nag
