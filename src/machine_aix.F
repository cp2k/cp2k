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

  USE kinds,                           ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: m_walltime, m_cputime, m_datum, m_flush, print_memory, &
            m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
            m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c

  INTERFACE m_loc_r
     MODULE PROCEDURE m_loc_r1,m_loc_r2,m_loc_r3,m_loc_r4
  END INTERFACE

  INTERFACE m_loc_c
     MODULE PROCEDURE m_loc_c1,m_loc_c2,m_loc_c3,m_loc_c4
  END INTERFACE

!!*****
!******************************************************************************

CONTAINS

FUNCTION m_loc_r1(a) RESULT(res)
    REAL(KIND=dp), DIMENSION(:), INTENT(in) :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_r1

FUNCTION m_loc_r2(a) RESULT(res)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(in)                             :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_r2

FUNCTION m_loc_r3(a) RESULT(res)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(in)                             :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_r3

FUNCTION m_loc_r4(a) RESULT(res)
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(in)                             :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_r4

FUNCTION m_loc_c1(a) RESULT(res)
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(in)                             :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_c1

FUNCTION m_loc_c2(a) RESULT(res)
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(in)                             :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_c2

FUNCTION m_loc_c3(a) RESULT(res)
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(in)                             :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_c3

FUNCTION m_loc_c4(a) RESULT(res)
    COMPLEX(KIND=dp), &
      DIMENSION(:, :, :, :), INTENT(in)      :: a
    INTEGER(kind=8)                          :: res

  res=INT(loc(a),8)
END FUNCTION m_loc_c4

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
    REAL(KIND=dp)                                :: ct

    REAL                                     :: wtl

  CALL CPU_TIME(wtl)
  ct=wtl
END FUNCTION m_cputime
!******************************************************************************
!!  wall clock time in seconds
FUNCTION m_walltime() RESULT (wt)
  REAL(KIND=dp) :: wt,omp_get_wtime
  INTEGER count,count_rate
#if defined(__OMPWTIME)
  wt = omp_get_wtime() ! this should not work in serial, or so but seems to work, presumably xlf_r
#else
  CALL SYSTEM_CLOCK(count=count,count_rate=count_rate)
  wt = REAL(count,KIND=dp)/count_rate
#endif
END FUNCTION m_walltime
!******************************************************************************
!!  time and date
SUBROUTINE m_datum(cal_date)
    CHARACTER(len=*), INTENT(OUT)            :: cal_date

    INTEGER                                  :: l, le

  CALL fdate_(cal_date)
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

    CALL flush_(lunit)

  END SUBROUTINE m_flush

! *****************************************************************************

  SUBROUTINE print_memory(checkpoint)

!   Purpose: Print the memory usage of the program at checkpoint.

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: checkpoint
#if defined(__OSX)
    !
    ! This is a hack to by-pass the problem of the xlf compiler
    ! on OS/X..
    !
    CHARACTER(LEN=200)                       :: cmdstr
    INTEGER                                  :: pid

    CALL m_getpid(pid)

    WRITE(UNIT=cmdstr,FMT="(A,I8,A,I8,A)")&
      "echo -e ""\\n Memory size of process""",pid,""" <"//TRIM(checkpoint)//&
      ">:"" `ps -o vsz -p",pid," | tail -1` KBytes""\\c"""

    CALL system(cmdstr)
#else
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
        WRITE(6,*) checkpoint,res4%usmblks+res4%uordblks,"/",res4%arena
    ELSE
        WRITE(6,*) checkpoint,res8%usmblks+res8%uordblks,"/",res8%arena
    ENDIF
#endif
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
SUBROUTINE m_chdir(dir,ierror)
    CHARACTER(len=*), INTENT(IN)             :: dir
    INTEGER, INTENT(OUT)                     :: ierror

    CHARACTER(LEN=1000)                      :: dir_local
    INTEGER                                  :: chdir

    dir_local=dir
    dir_local(LEN_TRIM(dir)+1:LEN_TRIM(dir)+1)=ACHAR(0)
    ierror = chdir(dir_local)
END SUBROUTINE m_chdir

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
