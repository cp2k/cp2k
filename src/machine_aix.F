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
!!     - print_memory changed (24.09.2002,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE machine_aix

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

    INTEGER, EXTERNAL                        :: iargc

  ic = iargc()
END FUNCTION m_iargc

!******************************************************************************
!!  cpu time in seconds
FUNCTION m_cputime() RESULT (ct)
    REAL(dbl)                                :: ct

    REAL                                     :: wtl

  CALL CPU_TIME(wtl)
  ct=wtl
END FUNCTION m_cputime
!******************************************************************************
!!  wall clock time in seconds
FUNCTION m_walltime() RESULT (wt)
  REAL (dbl) :: wt,omp_get_wtime
  integer count,count_rate
#if defined(__OMPWTIME)
  wt = omp_get_wtime() ! this should not work in serial, or so but seems to work, presumably xlf_r
#else
  CALL system_clock(count=count,count_rate=count_rate)
  wt = DBLE(count)/count_rate
#endif
END FUNCTION m_walltime
!******************************************************************************
!!  time and date
SUBROUTINE m_datum(cal_date)
    CHARACTER(len=*), INTENT(OUT)            :: cal_date

    INTEGER                                  :: l, le

  CALL fdate_(cal_date)
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

    CALL flush_(lunit)

  END SUBROUTINE m_flush

! *****************************************************************************

  SUBROUTINE print_memory(checkpoint)

!   Purpose: Print the memory usage of the program at checkpoint.

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: checkpoint

!   ---------------------------------------------------------------------------
! here there is definitely a problem with integer sizes. 
! they should be the C integer size, which is not really know
! and depends on the 32 (INTEGER*4) or 64 (INTEGER*8) compilation flag
! this should present a hack
! we use the fact that the default LOC kind (on IBM) is 4 or 8 in 32 or 64 bit mode

    TYPE memorytype4
       SEQUENCE
       INTEGER*4 arena,ordblks,smblks,hblks,hblkhd,usmblks,fsmblks,uordblks,fordblks,keepcost
    END TYPE memorytype4
    TYPE memorytype8
       SEQUENCE
       INTEGER*8 arena,ordblks,smblks,hblks,hblkhd,usmblks,fsmblks,uordblks,fordblks,keepcost
    END TYPE memorytype8
    TYPE(memorytype4) res4
    TYPE(memorytype8) res8,mallinfo
    EQUIVALENCE(res8,res4)
    
    res8=mallinfo()
    ! arena is the size of heap
    ! the sum is the actually allocated memory
    IF (KIND(LOC(res4%arena)).eq.4) THEN
        write(6,*) checkpoint,res4%usmblks+res4%uordblks,"/",res4%arena
    ELSE
        write(6,*) checkpoint,res8%usmblks+res8%uordblks,"/",res8%arena
    ENDIF

  END SUBROUTINE print_memory

! *****************************************************************************

SUBROUTINE m_hostnm(hname)
    CHARACTER(len=*), INTENT(OUT)            :: hname

    INTEGER                                  :: hostnm_, ierror

  ierror = hostnm_(hname)
END SUBROUTINE m_hostnm
!******************************************************************************
SUBROUTINE m_getcwd(curdir)
    CHARACTER(len=*), INTENT(OUT)            :: curdir

    INTEGER                                  :: getcwd_, ierror

  ierror = getcwd_(curdir)
END SUBROUTINE m_getcwd
!******************************************************************************
SUBROUTINE m_getlog(user)
    CHARACTER(len=*), INTENT(OUT)            :: user

  CALL getlog_(user)
END SUBROUTINE m_getlog
!******************************************************************************
SUBROUTINE m_getuid(uid)
    INTEGER, INTENT(OUT)                     :: uid

    INTEGER                                  :: getuid_

  uid = getuid_()
END SUBROUTINE m_getuid
!******************************************************************************
SUBROUTINE m_getpid(pid)
    INTEGER, INTENT(OUT)                     :: pid

    INTEGER                                  :: getpid_

  pid = getpid_()
END SUBROUTINE m_getpid
!******************************************************************************
SUBROUTINE m_getarg(i,arg)
    INTEGER, INTENT(IN)                      :: i
    CHARACTER(len=*), INTENT(OUT)            :: arg

  CALL getarg(i,arg)
END SUBROUTINE m_getarg
!******************************************************************************
END MODULE machine_aix
