!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/machine_irix [1.0] *
!!
!!   NAME
!!     machine_irix
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     APSI & JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE machine_irix
  
  USE kinds, ONLY : dbl, sgl
  
  PRIVATE
  PUBLIC :: m_walltime, m_cputime, m_datum
  PUBLIC :: m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!  cpu time in seconds
FUNCTION m_cputime() RESULT (ct)
  IMPLICIT NONE
  REAL (dbl) :: ct
  INTEGER :: mclock

  ct = mclock()*0.01_dbl
END FUNCTION m_cputime
!
!******************************************************************************
!!  wall clock time in seconds
FUNCTION m_walltime() RESULT (wt)
  IMPLICIT NONE
  REAL (dbl) :: wt
  REAL (sgl) :: etime
  TYPE tb_type
     SEQUENCE
     REAL (sgl) usrtime
     REAL (sgl) systime
  END TYPE tb_type
  TYPE (tb_type) etime_struct

  wt = etime(etime_struct)
END FUNCTION m_walltime
!******************************************************************************
!!  time and date
SUBROUTINE m_datum(cal_date)
  IMPLICIT NONE
  CHARACTER (len=*), INTENT (OUT) :: cal_date
  INTEGER :: l, le

  CALL fdate(cal_date)
  l = len(cal_date)
  le = len_trim(cal_date)
  cal_date(le:l) = ' '
END SUBROUTINE m_datum
!******************************************************************************
SUBROUTINE m_hostnm(hname)
  IMPLICIT NONE
  CHARACTER (len=*), INTENT (OUT) :: hname
  INTEGER :: hostnm, ierror

  ierror = hostnm(hname)
END SUBROUTINE m_hostnm
!******************************************************************************
SUBROUTINE m_getcwd(curdir)
  IMPLICIT NONE
  CHARACTER (len=*), INTENT (OUT) :: curdir
  INTEGER :: getcwd, ierror

  ierror = getcwd(curdir)
END SUBROUTINE m_getcwd
!******************************************************************************
SUBROUTINE m_getlog(user)
  IMPLICIT NONE
  CHARACTER (len=*), INTENT (OUT) :: user

  CALL getlog(user)
END SUBROUTINE m_getlog
!******************************************************************************
SUBROUTINE m_getuid(uid)
  IMPLICIT NONE
  INTEGER, INTENT (OUT) :: uid
  INTEGER :: getuid

  uid = getuid()
END SUBROUTINE m_getuid
!******************************************************************************
SUBROUTINE m_getpid(pid)
  IMPLICIT NONE
  INTEGER, INTENT (OUT) :: pid
  INTEGER :: getpid

  pid = getpid()
END SUBROUTINE m_getpid
!******************************************************************************
SUBROUTINE m_getarg(i,arg)
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: i
  CHARACTER (len=*), INTENT (OUT) :: arg

  CALL getarg(i,arg)
END SUBROUTINE m_getarg
!******************************************************************************

END MODULE machine_irix

!******************************************************************************
