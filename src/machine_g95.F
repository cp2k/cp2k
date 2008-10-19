!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      - m_flush added (12.06.2002,MK)
!>      - print_memory changed (24.09.2002,MK)
!>      - adapted for g95 (29.05.2003,JVdV)
!> \author APSI & JGH
! *****************************************************************************
MODULE machine_g95
  USE ISO_C_BINDING,                   ONLY: C_INT64_T
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             int_8

  IMPLICIT NONE

  INTEGER(KIND=C_INT64_T), bind(C, name='_g95_total_alloc') :: total_memory

  PRIVATE

  PUBLIC :: m_cputime, m_datum, m_flush, m_memory, &
            m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
            m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, total_memory, m_mov

CONTAINS

! *****************************************************************************
FUNCTION m_loc_r(a) RESULT(res)
    REAL(KIND=dp), DIMENSION(*), INTENT(in)  :: a
    INTEGER                                  :: res

    res=LOC(a)
END FUNCTION m_loc_r
  
! *****************************************************************************
FUNCTION m_loc_c(a) RESULT(res)
    COMPLEX(KIND=dp), DIMENSION(*), &
      INTENT(in)                             :: a
    INTEGER                                  :: res

    res=LOC(a)
END FUNCTION m_loc_c

! can be used to get a nice core
! *****************************************************************************
SUBROUTINE m_abort()
   CALL abort()
END SUBROUTINE m_abort

! the number of arguments of the fortran program
! *****************************************************************************
FUNCTION m_iargc() RESULT (ic)
    INTEGER                                  :: ic

    INTEGER                                  :: iargc

    ic = iargc()
END FUNCTION m_iargc

!!  cpu time in seconds
! *****************************************************************************
FUNCTION m_cputime() RESULT (ct)
    REAL(KIND=dp)                            :: ct

    CALL CPU_TIME(ct)
END FUNCTION m_cputime

!!  time and date
! *****************************************************************************
SUBROUTINE m_datum(cal_date)
    CHARACTER(len=*), INTENT(OUT)            :: cal_date

    CHARACTER(len=10)                        :: time
    CHARACTER(len=5)                         :: zone
    CHARACTER(len=8)                         :: date
    INTEGER, DIMENSION(8)                    :: values

    CALL DATE_AND_TIME(date=date, time=time, zone=zone, values=values)
    cal_date=date//" "//time
END SUBROUTINE m_datum

! flush a given unit
! *****************************************************************************
  SUBROUTINE m_flush(lunit)
    INTEGER, INTENT(IN)                      :: lunit

    CALL flush(lunit)
  END SUBROUTINE m_flush

! returns the total amount of memory [bytes] in use, if known, zero otherwise
! *****************************************************************************
  FUNCTION m_memory()
    INTEGER(KIND=int_8)                      :: m_memory

    m_memory=total_memory
  END FUNCTION m_memory

! *****************************************************************************
  SUBROUTINE m_mov(source,TARGET)

    CHARACTER(LEN=*), INTENT(IN)             :: source, TARGET

!   ---------------------------------------------------------------------------

    CALL rename(TRIM(source),TRIM(TARGET))

  END SUBROUTINE m_mov

! *****************************************************************************
SUBROUTINE m_hostnm(hname)
    CHARACTER(len=*), INTENT(OUT)            :: hname

    INTEGER                                  :: hostnm, ierror

    ierror=hostnm(hname)
END SUBROUTINE m_hostnm

! *****************************************************************************
SUBROUTINE m_getcwd(curdir)
    CHARACTER(len=*), INTENT(OUT)            :: curdir

    INTEGER                                  :: getcwd, ierror

    ierror = getcwd(curdir)
END SUBROUTINE m_getcwd

! *****************************************************************************
SUBROUTINE m_chdir(dir,ierror)
    CHARACTER(len=*), INTENT(IN)             :: dir
    INTEGER, INTENT(OUT)                     :: ierror

    INTEGER                                  :: chdir

    ierror = chdir(dir)
END SUBROUTINE m_chdir

! *****************************************************************************
SUBROUTINE m_getlog(user)
    CHARACTER(len=*), INTENT(OUT)            :: user

    CALL getlog(user)
END SUBROUTINE m_getlog

! *****************************************************************************
SUBROUTINE m_getuid(uid)
    INTEGER, INTENT(OUT)                     :: uid

    INTEGER                                  :: getuid

    uid = getuid()
END SUBROUTINE m_getuid

! *****************************************************************************
SUBROUTINE m_getpid(pid)
    INTEGER, INTENT(OUT)                     :: pid

    INTEGER                                  :: getpid

    pid = getpid()
END SUBROUTINE m_getpid

! *****************************************************************************
SUBROUTINE m_getarg(i,arg)
    INTEGER, INTENT(IN)                      :: i
    CHARACTER(len=*), INTENT(OUT)            :: arg

    CALL getarg(i,arg)
END SUBROUTINE m_getarg

END MODULE machine_g95

