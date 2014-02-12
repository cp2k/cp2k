!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      - m_flush added (12.06.2002,MK)
!>      - print_memory changed (24.09.2002,MK)
!> \author APSI & JGH
! *****************************************************************************
  USE f77_blas
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8,&
                                             sp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: m_cputime, m_flush, m_memory, &
            m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
            m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c,m_mov, m_memory_details, &
            m_procrun

CONTAINS
! *** get more detailed memory info, all units are bytes.
! *** the only 'useful' option is MemLikelyFree which is an estimate of remaining memory
! *** assumed to give info like /proc/meminfo while MeMLikelyFree is the amount of
! *** memory we're likely to be able to allocate, but not necessarily in one chunk
! *** zero means not available
  SUBROUTINE m_memory_details(MemTotal,MemFree,Buffers,Cached,Slab,SReclaimable,MemLikelyFree)

    INTEGER(kind=int_8), OPTIONAL            :: MemTotal, MemFree, Buffers, &
                                                Cached, Slab, SReclaimable, &
                                                MemLikelyFree

     MemTotal=0
     MemFree=0
     Buffers=0
     Cached=0
     Slab=0
     SReclaimable=0
     MemLikelyFree=0

  END SUBROUTINE m_memory_details

! returns if a process is running on the local machine
! 1 if yes and 0 if not
INTEGER FUNCTION m_procrun(id) RESULT (run_on)
    INTEGER           ::   id, ios
    CHARACTER(len=80) ::   filename, tmp
    CHARACTER(len=8)  ::   id_s

    run_on = 0

END FUNCTION m_procrun


! *****************************************************************************
FUNCTION m_loc_r(a) RESULT(res)
    REAL(KIND=dp), DIMENSION(*), INTENT(in)  :: a
    INTEGER                                  :: res

  res=-1
END FUNCTION m_loc_r

! *****************************************************************************
FUNCTION m_loc_c(a) RESULT(res)
    COMPLEX(KIND=dp), DIMENSION(*), &
      INTENT(in)                             :: a
    INTEGER                                  :: res

  res=-1
END FUNCTION m_loc_c

! can be used to get a nice core
! *****************************************************************************
SUBROUTINE m_abort()
   CALL abort()
END SUBROUTINE m_abort

! *****************************************************************************
FUNCTION m_iargc() RESULT (ic)
    INTEGER                                  :: ic

    INTEGER, EXTERNAL                        :: iargc

  ic = iargc()
END FUNCTION m_iargc

!!  cpu time in seconds
! *****************************************************************************
FUNCTION m_cputime() RESULT (ct)
    REAL(KIND=dp)                            :: ct

    REAL(sp)                                 :: d, etime, tarray(2)

  EXTERNAL etime
  d=ETIME(tarray)
  ct = REAL(tarray(1),KIND=dp)
END FUNCTION m_cputime

! *****************************************************************************
!> \brief   Flush the output to a logical unit.
!> \author  MK
!> \date    14.10.1999
!> \version 1.0
! *****************************************************************************
  SUBROUTINE m_flush(lunit)
    INTEGER, INTENT(IN)                      :: lunit

    CALL flush(lunit)

  END SUBROUTINE m_flush

! returns the total amount of memory [bytes] in use, if known, zero otherwise
! *****************************************************************************
  FUNCTION m_memory()
    INTEGER(KIND=int_8)                      :: m_memory

    m_memory=0
  END FUNCTION m_memory

! *****************************************************************************
  SUBROUTINE m_mov(source,TARGET)

    CHARACTER(LEN=*), INTENT(IN)             :: source, TARGET

    CHARACTER(LEN=2*default_string_length+4) :: cmd

    cmd = "mv " // source(1:LEN_TRIM(source)) // " " // TARGET(1:LEN_TRIM(TARGET))
    CALL system(cmd)

  END SUBROUTINE m_mov

! *****************************************************************************
SUBROUTINE m_hostnm(hname)
    CHARACTER(len=*), INTENT(OUT)            :: hname

    INTEGER                                  :: hostnm, ierror

  ierror = hostnm(hname)
END SUBROUTINE m_hostnm
! *****************************************************************************
SUBROUTINE m_getcwd(curdir)
    CHARACTER(len=*), INTENT(OUT)            :: curdir

    CHARACTER(len=255)                       :: dir
    INTEGER                                  :: getcwd, ierror

  ierror = getcwd(dir)
  curdir = dir
END SUBROUTINE m_getcwd
! *****************************************************************************
SUBROUTINE m_chdir(dir,ierror)
    CHARACTER(len=*), INTENT(IN)             :: dir
    INTEGER, INTENT(OUT)                     :: ierror

! ierror = chdir(dir)

    STOP
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
