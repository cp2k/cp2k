!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/machine_pgi [1.0] *
!!
!!   NAME
!!     machine_pgi
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
!!*****
!******************************************************************************

MODULE machine_pgi

  USE kinds,                           ONLY: dbl,&
                                             sgl

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: m_walltime, m_cputime, m_datum, m_flush, print_memory, &
            m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, m_abort, m_iargc

!!*****
!******************************************************************************

CONTAINS

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
    REAL(dbl)                                :: ct

    INTEGER                                  :: mclock

  ct = mclock()*0.01_dbl
END FUNCTION m_cputime
!******************************************************************************
!!  wall clock time in seconds
FUNCTION m_walltime() RESULT (wt)
    REAL(dbl)                                :: wt

    REAL(sgl)                                :: etime

  TYPE tb_type
     SEQUENCE
     REAL (sgl) :: usrtime
     REAL (sgl) :: systime
  END TYPE tb_type
  TYPE (tb_type) etime_struct
  wt = etime(etime_struct)
END FUNCTION m_walltime
!******************************************************************************
!!  time and date
SUBROUTINE m_datum(cal_date)
    CHARACTER(len=*), INTENT(OUT)            :: cal_date

    INTEGER                                  :: l, le

  CALL fdate(cal_date)
  l = len(cal_date)
  le = len_trim(cal_date)
  cal_date(le:l) = ' '
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
END MODULE machine_pgi
