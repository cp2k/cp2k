!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Implemenation of machine interface based on Fortran 2003 and POSIX
!> \author Ole Schuett
! *****************************************************************************
  USE kinds,                           ONLY: dp, int_8, default_path_length,&
                                             default_string_length
  USE ISO_C_BINDING, ONLY: C_INT, C_NULL_CHAR, C_CHAR, C_PTR, C_NULL_PTR, C_ASSOCIATED, C_F_POINTER

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: m_flush, m_memory, &
            m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
            m_iargc, m_abort, m_chdir, m_mov, &
            m_memory_details, m_procrun

  INTEGER(KIND=int_8), PUBLIC, SAVE :: m_memory_max=0

CONTAINS

! *****************************************************************************
! can be used to get a nice core
! *****************************************************************************
  SUBROUTINE m_abort()
    INTERFACE
      SUBROUTINE abort() BIND(C,name="abort")
      END SUBROUTINE
    END INTERFACE

    CALL abort()
  END SUBROUTINE m_abort


! *****************************************************************************
! the number of arguments of the fortran program
! *****************************************************************************
  FUNCTION m_iargc() RESULT (ic)
    INTEGER                                  :: ic

    ic = COMMAND_ARGUMENT_COUNT ()
  END FUNCTION m_iargc


! *****************************************************************************
! flush a given unit
! *****************************************************************************
  SUBROUTINE m_flush(lunit)
    INTEGER, INTENT(IN)                      :: lunit

    FLUSH(lunit)
  END SUBROUTINE m_flush


! *****************************************************************************
! returns if a process is running on the local machine
! 1 if yes and 0 if not
! *****************************************************************************
  FUNCTION m_procrun(pid) RESULT (run_on)
    INTEGER, INTENT(IN)       ::   pid
    INTEGER                   ::   run_on
    INTEGER                   ::   istat

    INTERFACE
      FUNCTION kill(pid, sig) BIND(C,name="kill") RESULT(errno)
        IMPORT
        INTEGER(KIND=C_INT),VALUE                :: pid, sig
        INTEGER(KIND=C_INT)                      :: errno
      END FUNCTION
    END INTERFACE

    ! If sig is 0, then no signal is sent, but error checking is still
    ! performed; this can be used to check for the existence of a process
    ! ID or process group ID.

    istat = kill(pid=pid, sig=0)
    IF(istat == 0) THEN
       run_on = 1 ! no error, process exists
    ELSE
       run_on = 0 ! error, process probably does not exist
    ENDIF

  END FUNCTION m_procrun


! *****************************************************************************
! returns the total amount of memory [bytes] in use, if known, zero otherwise
! *****************************************************************************
  SUBROUTINE m_memory(mem)

      INTEGER(KIND=int_8), OPTIONAL, INTENT(OUT)         :: mem
      INTEGER(KIND=int_8)                      :: m1,m2,m3,mem_local

      !
      ! __NO_STATM_ACCESS can be used to disable the stuff, if getpagesize
      ! lead to linking errors or /proc/self/statm can not be opened
      !
#if defined(__NO_STATM_ACCESS)
      mem_local=0
#else
      CHARACTER(LEN=80) :: DATA
      INTEGER :: iostat,i

      ! the size of a page, might not be available everywhere
      INTERFACE
       FUNCTION getpagesize() BIND(C,name="getpagesize") RESULT(RES)
         IMPORT
         INTEGER(C_INT) :: RES
       END FUNCTION
      END INTERFACE

      !
      ! reading from statm
      !
      mem_local=-1
      DATA=""
      OPEN(121245,FILE="/proc/self/statm",ACTION="READ",STATUS="OLD",ACCESS="STREAM")
      DO I=1,80
         READ(121245,END=999) DATA(I:I)
      ENDDO
999   CLOSE(121245)
      DATA(I:80)=""
      ! m1 = total
      ! m2 = resident
      ! m3 = shared
      READ(DATA,*,IOSTAT=iostat) m1,m2,m3
      IF (iostat.NE.0) THEN
         mem_local=0
      ELSE
         mem_local=m2
#if defined(__STATM_TOTAL)
         mem_local=m1
#endif
#if defined(__STATM_RESIDENT)
         mem_local=m2
#endif
         mem_local=mem_local*getpagesize()
      ENDIF
#endif

      m_memory_max=MAX(mem_local,m_memory_max)
      IF (PRESENT(mem)) mem=mem_local

  END SUBROUTINE m_memory

! *****************************************************************************
! *** get more detailed memory info, all units are bytes.
! *** the only 'useful' option is MemLikelyFree which is an estimate of remaining memory
! *** assumed to give info like /proc/meminfo while MeMLikelyFree is the amount of
! *** memory we're likely to be able to allocate, but not necessarily in one chunk
! *** zero means not available
! *****************************************************************************
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
  END SUBROUTINE m_memory_details


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE m_mov(source,TARGET)

    CHARACTER(LEN=*), INTENT(IN)             :: source, TARGET

    INTEGER                                  :: istat

    INTERFACE
      FUNCTION unlink(path) BIND(C,name="unlink") RESULT(errno)
        IMPORT
        CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: path
        INTEGER(KIND=C_INT)                      :: errno
      END FUNCTION
    END INTERFACE

    INTERFACE
      FUNCTION rename(src, dest) BIND(C,name="rename") RESULT(errno)
        IMPORT
        CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: src, dest
        INTEGER(KIND=C_INT)                      :: errno
      END FUNCTION
    END INTERFACE

    IF (TARGET==source) THEN
       WRITE(*,*) "Warning: m_mov ",TRIM(TARGET)," equals ", TRIM(source)
       RETURN
    ENDIF

    ! first remove target (needed on windows / mingw)
    istat = unlink(TRIM(TARGET)//c_null_char)
    ! ignore istat of unlink

    ! now move
    istat = rename(TRIM(source)//c_null_char, TRIM(TARGET)//c_null_char)
    IF (istat .NE. 0) THEN
      WRITE(*,*) "Trying to move "//TRIM(source)//" to "//TRIM(TARGET)//"."
      WRITE(*,*) "rename returned status: ",istat
      WRITE(*,*) "Problem moving file"
      CALL m_abort()
    ENDIF
  END SUBROUTINE m_mov


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE m_hostnm(hname)
    CHARACTER(len=*), INTENT(OUT)            :: hname

    INTEGER                                  :: istat, i
    CHARACTER(len=default_path_length)       :: buf

    INTERFACE
      FUNCTION  gethostname(buf, buflen) BIND(C,name="gethostname") RESULT(errno)
        IMPORT
        CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: buf
        INTEGER(KIND=C_INT), VALUE               :: buflen
        INTEGER(KIND=C_INT)                      :: errno
      END FUNCTION
    END INTERFACE

    istat = gethostname(buf, LEN(buf))
    IF(istat /= 0) THEN
       WRITE (*,*) "m_hostnm failed"
       CALL m_abort()
    ENDIF
    i = INDEX(buf, c_null_char) -1
    hname = buf(1:i)
  END SUBROUTINE m_hostnm


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE m_getcwd(curdir)
    CHARACTER(len=*), INTENT(OUT)            :: curdir
    TYPE(C_PTR)                              :: stat
    INTEGER                                  :: i
    CHARACTER(len=default_path_length), TARGET  :: tmp

    INTERFACE
      FUNCTION  getcwd(buf, buflen) BIND(C,name="getcwd") RESULT(stat)
        IMPORT
        CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: buf
        INTEGER(KIND=C_INT), VALUE               :: buflen
        TYPE(C_PTR)                              :: stat
      END FUNCTION
    END INTERFACE

    stat = getcwd(tmp, LEN(tmp))
    IF(.NOT. C_ASSOCIATED(stat)) THEN
       WRITE (*,*) "m_getcwd failed"
       CALL m_abort()
    ENDIF
    i = INDEX(tmp, c_null_char) -1
    curdir = tmp(1:i)
  END SUBROUTINE m_getcwd


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE m_chdir(dir,ierror)
    CHARACTER(len=*), INTENT(IN)             :: dir
    INTEGER, INTENT(OUT)                     :: ierror

    INTERFACE
      FUNCTION chdir(path) BIND(C,name="chdir") RESULT(errno)
        IMPORT
        CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: path
        INTEGER(KIND=C_INT)                      :: errno
      END FUNCTION
    END INTERFACE

    ierror = chdir(TRIM(dir)//c_null_char)
  END SUBROUTINE m_chdir


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE m_getlog(user)
    CHARACTER(len=*), INTENT(OUT)            :: user
    INTEGER :: status

    ! on a posix system LOGNAME should be defined
    CALL get_environment_variable("LOGNAME", value=user, status=status)
    ! nope, check alternative
    IF (status/=0) &
       CALL get_environment_variable("USER", value=user, status=status)
    ! fall back
    IF (status/=0) &
       user="root ;-)"

  END SUBROUTINE m_getlog


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE m_getuid(uid)
   INTEGER, INTENT(OUT)                     :: uid

   INTERFACE
     FUNCTION getuid() BIND(C,name="getuid") RESULT(uid)
       IMPORT
       INTEGER(KIND=C_INT)              :: uid
     END FUNCTION
   END INTERFACE

   uid = getuid()
  END SUBROUTINE m_getuid


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE m_getpid(pid)
   INTEGER, INTENT(OUT)                     :: pid

   INTERFACE
     FUNCTION getpid() BIND(C,name="getpid") RESULT(pid)
       IMPORT
       INTEGER(KIND=C_INT)              :: pid
     END FUNCTION
   END INTERFACE

   pid = getpid()
  END SUBROUTINE m_getpid


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE m_getarg(i,arg)
    INTEGER, INTENT(IN)                      :: i
    CHARACTER(len=*), INTENT(OUT)            :: arg
    CHARACTER(len=1024)                      :: tmp
    INTEGER                                  :: istat

    CALL GET_COMMAND_ARGUMENT(i, tmp, status=istat)
    IF(istat /= 0) THEN
       WRITE (*,*) "m_getarg failed"
       CALL m_abort()
    ENDIF
    arg = TRIM(tmp)
  END SUBROUTINE m_getarg

