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
  USE kinds,                           ONLY: dp,&
                                             int_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: m_cputime, m_flush, m_memory, &
            m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg,&
            m_abort, m_iargc, m_chdir, m_loc_r, m_loc_c,m_mov, m_memory_details, &
            m_procrun

CONTAINS

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

#if defined(__parallel)
    REAL(KIND=dp), EXTERNAL                  :: MPI_WTIME

    ct = MPI_WTIME()
#else
    INTEGER                                  :: mclock

    ct = mclock()*0.01_dp
#endif
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

      !
      ! __NO_STATM_ACCESS can be used to disable the stuff, if getpagesize
      ! lead to linking errors or /proc/self/statm can not be opened
      !
#if defined(__NO_STATM_ACCESS) || defined (__HAS_NO_ISO_C_BINDING)
      m_memory=0
#else
      CHARACTER(LEN=80) :: DATA
      INTEGER :: iostat,i

      ! the size of a page, might not be available everywhere
      INTERFACE
       FUNCTION getpagesize() BIND(C,name="getpagesize") RESULT(RES)
         USE ISO_C_BINDING
         INTEGER(C_INT) :: RES
       END FUNCTION
      END INTERFACE

      !
      ! reading from statm
      !
      m_memory=-1
      DATA=""
      OPEN(121245,FILE="/proc/self/statm",ACTION="READ",STATUS="OLD",ACCESS="STREAM")
      DO I=1,80
         READ(121245,END=999) DATA(I:I)
      ENDDO
999   CLOSE(121245)
      DATA(I:80)=""
      READ(DATA,*,IOSTAT=iostat) m_memory
      IF (iostat.NE.0) THEN
         m_memory=0
      ELSE
         m_memory=m_memory*getpagesize()
      ENDIF
#endif

  END FUNCTION m_memory

! *** get more detailed memory info, all units are bytes.
! *** the only 'useful' option is MemLikelyFree which is an estimate of remaining memory
! *** assumed to give info like /proc/meminfo while MeMLikelyFree is the amount of
! *** memory we're likely to be able to allocate, but not necessarily in one chunk
! *** zero means not available
  SUBROUTINE m_memory_details(MemTotal,MemFree,Buffers,Cached,Slab,SReclaimable,MemLikelyFree)

     INTEGER(kind=int_8), OPTIONAL :: MemTotal,MemFree,Buffers,Cached,Slab,SReclaimable,MemLikelyFree

     INTEGER, PARAMETER :: Nbuffer=10000
     CHARACTER(LEN=Nbuffer) :: meminfo


     INTEGER :: i

     MemTotal=0
     MemFree=0
     Buffers=0
     Cached=0
     Slab=0
     SReclaimable=0
     MemLikelyFree=0
     meminfo=""

#ifndef __NO_STATM_ACCESS
     OPEN(UNIT=8123,file="/proc/meminfo",ACCESS="STREAM",ERR=901)
     i=0
     DO
       i=i+1
       IF (i>Nbuffer) EXIT
       READ(8123,END=900,ERR=900) meminfo(i:i)
     ENDDO
 900 CONTINUE
     meminfo(i:Nbuffer)=""
 901 CONTINUE
     CLOSE(8123,ERR=902)
 902 CONTINUE
     MemTotal=get_field_value_in_bytes('MemTotal:')
     MemFree=get_field_value_in_bytes('MemFree:')
     Buffers=get_field_value_in_bytes('Buffers:')
     Cached=get_field_value_in_bytes('Cached:')
     Slab=get_field_value_in_bytes('Slab:')
     SReclaimable=get_field_value_in_bytes('SReclaimable:')
     ! opinions here vary but this might work
     MemLikelyFree=MemFree+Buffers+Cached+SReclaimable


  CONTAINS
        INTEGER(int_8) FUNCTION get_field_value_in_bytes(field)
           CHARACTER(LEN=*) :: field
           INTEGER :: start
           INTEGER(KIND=int_8) :: value
           get_field_value_in_bytes=0
           start=INDEX(meminfo,field)
           IF (start.NE.0) THEN
              start=start+LEN_TRIM(field)
              IF (start.LT.Nbuffer) THEN
                 READ(meminfo(start:),*,ERR=999,END=999) value
                 ! XXXXXXX convert from Kb to bytes XXXXXXXX
                 get_field_value_in_bytes=value*1024
 999             CONTINUE
              ENDIF
           ENDIF
        END FUNCTION
#endif

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
  SUBROUTINE m_mov(source,TARGET)

    CHARACTER(LEN=*), INTENT(IN)             :: source, TARGET

    CALL rename(source(1:LEN_TRIM(source)), TARGET(1:LEN_TRIM(TARGET)))

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
