!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/machine_t3e [1.0] *
!!
!!   NAME
!!     machine_t3e
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

MODULE machine_t3e

  USE kinds,                           ONLY: dbl

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: m_walltime, m_cputime, m_datum, m_flush, print_memory, &
            m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
            m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c

!!*****
!******************************************************************************

CONTAINS


FUNCTION m_loc_r(a) RESULT(res)
  REAL(kind=dbl), INTENT(in), DIMENSION(*) :: a
  integer :: res

  res=-1
END FUNCTION m_loc_r
  
FUNCTION m_loc_c(a) result(res)
  COMPLEX(kind=dbl), INTENT(in), DIMENSION(*) :: a
  integer :: res

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
    REAL(dbl)                                :: ct

    REAL(dbl)                                :: tsecnd

  ct = tsecnd()
END FUNCTION m_cputime
!******************************************************************************
!!  wall clock time in seconds
FUNCTION m_walltime() RESULT (wt)
    REAL(dbl)                                :: wt

    REAL(dbl)                                :: time

  EXTERNAL time
  wt = time()
END FUNCTION m_walltime
!******************************************************************************
!!  time and date
SUBROUTINE m_datum(cal_date)
    CHARACTER(len=*), INTENT(OUT)            :: cal_date

    INTEGER                                  :: l, le

  WRITE (cal_date,'(a8,'' on '',a8)') clock(), date()
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
!MK CALL flush(lunit)

  END SUBROUTINE m_flush

! *****************************************************************************

  SUBROUTINE print_memory(checkpoint)

!   Purpose: Print the memory usage of the program at checkpoint.

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: checkpoint

!   ---------------------------------------------------------------------------
!MK CALL m_getpid(pid)
!MK WRITE(UNIT=cmdstr,FMT="(A,I8,A,I8,A)")&
!MK   "echo -e ""\\n Memory size of process""",pid,""" <"//TRIM(checkpoint)//&
!MK   ">:"" `ps -o vsz -p",pid," | tail -1` KBytes""\\c"""
!MK CALL ishell(cmdstr)

  END SUBROUTINE print_memory

! *****************************************************************************

SUBROUTINE m_hostnm(hname)
    CHARACTER(len=*), INTENT(OUT)            :: hname

  CALL gethost(hname)
END SUBROUTINE m_hostnm
!******************************************************************************
SUBROUTINE m_getcwd(curdir)
    CHARACTER(len=*), INTENT(OUT)            :: curdir

  CALL getcwd(curdir)
END SUBROUTINE m_getcwd
!******************************************************************************
SUBROUTINE m_chdir(dir,ierror)
    CHARACTER(len=*), INTENT(IN)          :: dir
    INTEGER,          INTENT(OUT)            :: ierror
    INTEGER                                  :: chdir
    ! ierror = chdir(dir)  ! adapt as needed
    STOP
END SUBROUTINE m_chdir
!******************************************************************************
SUBROUTINE m_getlog(user)
    CHARACTER(len=*), INTENT(OUT)            :: user

    INTEGER                                  :: ierror, ilen

  CALL pxfgetlogin(user,ilen,ierror)
END SUBROUTINE m_getlog
!******************************************************************************
SUBROUTINE m_getuid(uid)
    INTEGER, INTENT(OUT)                     :: uid

    INTEGER                                  :: ierror

  CALL pxfgetuid(uid,ierror)
END SUBROUTINE m_getuid
!******************************************************************************
SUBROUTINE m_getpid(pid)
    INTEGER, INTENT(OUT)                     :: pid

    INTEGER                                  :: ierror

  CALL pxfgetpid(pid,ierror)
END SUBROUTINE m_getpid
!******************************************************************************
SUBROUTINE m_getarg(i,arg)
    INTEGER, INTENT(IN)                      :: i
    CHARACTER(len=*), INTENT(OUT)            :: arg

    INTEGER                                  :: ierror, ilen

  CALL pxfgetarg(i,arg,ilen,ierror)
END SUBROUTINE m_getarg
!******************************************************************************
END MODULE machine_t3e
