!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/machine_aix [1.0] *
!!
!!   NAME
!!     machine_aix
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     APSI & JGH
!!
!!   MODIFICATION HISTORY
!!     - m_flush added (12.06.2002,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE machine_aix

  USE kinds, ONLY : dbl, sgl

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: m_walltime, m_cputime, m_datum, m_flush, print_memory, &
            m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!  cpu time in seconds
FUNCTION m_cputime() RESULT (ct)
  REAL (dbl) :: ct
  INTEGER :: mclock
  ct = mclock()*0.01_dbl
END FUNCTION m_cputime
!******************************************************************************
!!  wall clock time in seconds
FUNCTION m_walltime() RESULT (wt)
  REAL (dbl) :: wt
  REAL (sgl) :: etime_
  TYPE tb_type
     SEQUENCE
     REAL (sgl) usrtime
     REAL (sgl) systime
  END TYPE tb_type
  TYPE (tb_type) etime_struct
  wt = etime_(etime_struct)
END FUNCTION m_walltime
!******************************************************************************
!! memory usage
SUBROUTINE print_memory(subr)
  CHARACTER(len=24), INTENT(IN) :: subr
!locals:
  CHARACTER(len=100):: string
  INTEGER :: ihpstat,getpid_,my_pid
  MY_PID=GETPID_()
  WRITE(6,'(A,A,T25,A,$)') ' ***',SUBR,' the size of the program is '
  WRITE(STRING,'(A,I10,A)') 'ps -Fvsz -p ',MY_PID,&
       ' | sed -e "/VSZ/d" | sed -e "s/$/ KBytes ***/" '
  CALL SYSTEM(STRING)
END SUBROUTINE print_memory
!******************************************************************************
!!  time and date
SUBROUTINE m_datum(cal_date)
  CHARACTER (len=*), INTENT (OUT) :: cal_date
  INTEGER :: l, le
  CALL fdate_(cal_date)
  l = len(cal_date)
  le = len_trim(cal_date)
  cal_date(le:l) = ' '
END SUBROUTINE m_datum
!******************************************************************************

! *** Flush the output to a logical unit ***

  SUBROUTINE m_flush(lunit)

    INTEGER, INTENT(IN) :: lunit

    CALL flush_(lunit)

  END SUBROUTINE m_flush

!******************************************************************************
SUBROUTINE m_hostnm(hname)
  CHARACTER (len=*), INTENT (OUT) :: hname
  INTEGER :: hostnm_, ierror
  ierror = hostnm_(hname)
END SUBROUTINE m_hostnm
!******************************************************************************
SUBROUTINE m_getcwd(curdir)
  CHARACTER (len=*), INTENT (OUT) :: curdir
  INTEGER :: getcwd_, ierror
  ierror = getcwd_(curdir)
END SUBROUTINE m_getcwd
!******************************************************************************
SUBROUTINE m_getlog(user)
  CHARACTER (len=*), INTENT (OUT) :: user
  CALL getlog_(user)
END SUBROUTINE m_getlog
!******************************************************************************
SUBROUTINE m_getuid(uid)
  INTEGER, INTENT (OUT) :: uid
  INTEGER :: getuid_
  uid = getuid_()
END SUBROUTINE m_getuid
!******************************************************************************
SUBROUTINE m_getpid(pid)
  INTEGER, INTENT (OUT) :: pid
  INTEGER :: getpid_
  pid = getpid_()
END SUBROUTINE m_getpid
!******************************************************************************
SUBROUTINE m_getarg(i,arg)
  INTEGER, INTENT (IN) :: i
  CHARACTER (len=*), INTENT (OUT) :: arg
  CALL getarg(i,arg)
END SUBROUTINE m_getarg
!******************************************************************************
END MODULE machine_aix
