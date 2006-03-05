!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/machine_sun [1.0] *
!!
!!   NAME
!!     machine_sun
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     APSI & JGH
!!
!!   MODIFICATION HISTORY
!!     - m_flush added (12.06.2002,MK)
!!     - print_memory changed (24.09.2002,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE machine_sun

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
    REAL(KIND=dp), DIMENSION(*), INTENT(in)  :: a
    INTEGER                                  :: res

  res=-1
END FUNCTION m_loc_r
  
FUNCTION m_loc_c(a) RESULT(res)
    COMPLEX(KIND=dp), DIMENSION(*), &
      INTENT(in)                             :: a
    INTEGER                                  :: res

  res=-1
END FUNCTION m_loc_c

! can be used to get a nice core
SUBROUTINE m_abort()
   CALL abort()
END SUBROUTINE m_abort


FUNCTION m_iargc() RESULT (ic)
    INTEGER                                  :: ic

    INTEGER, EXTERNAL                        :: iargc

  ic = iargc()
END FUNCTION m_iargc

!******************************************************************************
!!  cpu time in seconds
FUNCTION m_cputime() RESULT (ct)
    REAL(KIND=dp)                            :: ct

    REAL                                     :: cto

! REAL ( KIND = 4 ) :: etime, timearray ( 2 )
! ct = REAL ( etime ( timearray ),KIND=dp)
! ct = REAL ( timearray (1),KIND=dp)

  CALL cpu_time(cto)
  ct = cto
END FUNCTION m_cputime
!******************************************************************************
!!  wall clock time in seconds
FUNCTION m_walltime() RESULT (wt)
    REAL(KIND=dp)                            :: wt

    INTEGER*8                                :: count, count_rate

! wt = omp_get_wtime() ! omp call is fast and accurate

  CALL SYSTEM_CLOCK(count=count,count_rate=count_rate)
  wt=REAL(count,dp)/REAL(count_rate,dp)
END FUNCTION m_walltime
!******************************************************************************
!!  time and date
SUBROUTINE m_datum(cal_date)
    CHARACTER(len=*), INTENT(OUT)            :: cal_date

    INTEGER                                  :: l, le

  CALL fdate(cal_date)
  l = LEN(cal_date)
  le = LEN_TRIM(cal_date)
  cal_date(le:l) = ' '
END SUBROUTINE m_datum

! *****************************************************************************

  SUBROUTINE m_flush(lunit)

!   Purpose: Flush the output to a logical unit.

!   History: - Creation (14.10.1999,MK)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: lunit

!   ---------------------------------------------------------------------------
! there seems to be a problem with calling flush too often (MD) ....
! CALL flush(lunit)

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
     "/usr/bin/echo ""\n PID""",pid,""" <"//TRIM(checkpoint)//&
     ">:"" `ps -o vsz -p",pid," | tail -1` KBytes""\c"""

   CALL system(cmdstr)

  END SUBROUTINE print_memory

! *****************************************************************************

SUBROUTINE m_hostnm(hname)
    CHARACTER(len=*), INTENT(OUT)            :: hname

    INTEGER                                  :: hostnm, ierror

  ierror = hostnm(hname)
END SUBROUTINE m_hostnm
!******************************************************************************
SUBROUTINE m_getcwd(curdir)
    CHARACTER(len=*), INTENT(OUT)            :: curdir

    INTEGER                                  :: getcwd, ierror

  ierror = getcwd(curdir)
END SUBROUTINE m_getcwd
!******************************************************************************
SUBROUTINE m_chdir(dir,ierror)
    CHARACTER(len=*), INTENT(IN)             :: dir
    INTEGER, INTENT(OUT)                     :: ierror

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

    INTEGER                                  :: getuid

  uid = getuid()
END SUBROUTINE m_getuid
!******************************************************************************
SUBROUTINE m_getpid(pid)
    INTEGER, INTENT(OUT)                     :: pid

    INTEGER                                  :: getpid

  pid = getpid()
END SUBROUTINE m_getpid
!******************************************************************************
SUBROUTINE m_getarg(i,arg)
    INTEGER, INTENT(IN)                      :: i
    CHARACTER(len=*), INTENT(OUT)            :: arg

  CALL getarg(i,arg)
END SUBROUTINE m_getarg
!******************************************************************************
END MODULE machine_sun
