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
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE machine_t3e
  
  USE kinds, ONLY : dbl, sgl
  
  PRIVATE
  PUBLIC :: m_walltime, m_cputime, m_datum, print_memory
  PUBLIC :: m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg
  
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!  cpu time in seconds
FUNCTION m_cputime() RESULT (ct)
  IMPLICIT NONE
  REAL (dbl) :: ct
  REAL (dbl) :: tsecnd

  ct = tsecnd()
END FUNCTION m_cputime
!
!******************************************************************************
!!  wall clock time in seconds
FUNCTION m_walltime() RESULT (wt)
  IMPLICIT NONE
  REAL (dbl) :: wt
  REAL (dbl) :: time
  EXTERNAL time

  wt = time()
END FUNCTION m_walltime
!******************************************************************************
!! memory usage
SUBROUTINE print_memory(subr)
  IMPLICIT NONE
  CHARACTER(len=24), INTENT(IN) :: subr
!locals:
  CHARACTER(len=200):: string
  INTEGER :: ihpstat,getpid_,my_pid, ierror
  CALL pxfgetpid ( MY_PID, ierror )
!*apsi  WRITE(6,'(A,A,T25,A,$)') ' ***',SUBR,' the size of the program is '
!*apsi        WRITE(STRING,'(A,I10,A)') 'ps -o vsz -p ', MY_PID, &
!*apsi          ' | sed -e "/VSZ/d" | sed -e "s/$/ KBytes ***/" '
  WRITE(STRING,'(A,I10,3A)') 'ps -o vsz -p ', MY_PID, &
       ' | sed -e "/VSZ/d" | sed -e "s/$/ KBytes ***/" | sed -e "s/^/ ***', &
       subr, ' the size of the program is /"'
  CALL ishell(string)
END SUBROUTINE print_memory
!******************************************************************************
!!  time and date
SUBROUTINE m_datum(cal_date)
  IMPLICIT NONE
  CHARACTER (len=*), INTENT (OUT) :: cal_date
  INTEGER :: l, le

  WRITE (cal_date,'(a8,'' on '',a8)') clock(), date()
  l = len(cal_date)
  le = len_trim(cal_date)
  cal_date(le:l) = ' '
END SUBROUTINE m_datum
!******************************************************************************
SUBROUTINE m_hostnm(hname)
  IMPLICIT NONE
  CHARACTER (len=*), INTENT (OUT) :: hname

  CALL gethost(hname)
END SUBROUTINE m_hostnm
!******************************************************************************
SUBROUTINE m_getcwd(curdir)
  IMPLICIT NONE
  CHARACTER (len=*), INTENT (OUT) :: curdir

  CALL getcwd(curdir)
END SUBROUTINE m_getcwd
!******************************************************************************
SUBROUTINE m_getlog(user)
  IMPLICIT NONE
  CHARACTER (len=*), INTENT (OUT) :: user
  INTEGER :: ilen, ierror

  CALL pxfgetlogin(user,ilen,ierror)
END SUBROUTINE m_getlog
!******************************************************************************
SUBROUTINE m_getuid(uid)
  IMPLICIT NONE
  INTEGER, INTENT (OUT) :: uid
  INTEGER :: ierror

  CALL pxfgetuid(uid,ierror)
END SUBROUTINE m_getuid
!******************************************************************************
SUBROUTINE m_getpid(pid)
  IMPLICIT NONE
  INTEGER, INTENT (OUT) :: pid
  INTEGER :: ierror

  CALL pxfgetpid(pid,ierror)
END SUBROUTINE m_getpid
!******************************************************************************
SUBROUTINE m_getarg(i,arg)
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: i
  CHARACTER (len=*), INTENT (OUT) :: arg
  INTEGER :: ilen, ierror

  CALL pxfgetarg(i,arg,ilen,ierror)
END SUBROUTINE m_getarg
!******************************************************************************

END MODULE machine_t3e

!******************************************************************************

!*apsi* Since this one does not exist...
SUBROUTINE gsyev ( jobz, uplo, n, a, lda, w, work, lwork, info )
  USE stop_program, ONLY : stop_prg
  IMPLICIT NONE
  CHARACTER :: jobz, uplo
  INTEGER :: INFO, LDA, LWORK, N
  REAL ( KIND = 4 ) :: A( LDA, * ), W( * ), WORK( * )
  
  CALL stop_prg ( "gsyev", "not defined" )
END SUBROUTINE gsyev
