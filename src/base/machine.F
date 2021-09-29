!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Machine interface based on Fortran 2003 and POSIX
!> \par History
!>      JGH (05.07.2001) : added G95 interface
!>      - m_flush added (12.06.2002,MK)
!>      - Missing print_memory added (24.09.2002,MK)
!>      - Migrate to generic implementation based on F2003 + POSIX (2014, Ole Schuett)
!> \author APSI, JGH, Ole Schuett
! **************************************************************************************************
MODULE machine
   USE ISO_C_BINDING, ONLY: C_CHAR, &
                            C_INT, &
                            C_PTR, &
                            C_NULL_CHAR, &
                            C_ASSOCIATED
   USE ISO_FORTRAN_ENV, ONLY: input_unit, &
                              output_unit
   USE omp_lib, ONLY: omp_get_wtime
   USE kinds, ONLY: default_path_length, &
                    default_string_length, &
                    dp, &
                    int_8
#if defined(__LIBXSMM)
   USE libxsmm, ONLY: libxsmm_timer_tick, &
                      libxsmm_timer_duration, &
                      libxsmm_get_target_archid, &
                      LIBXSMM_X86_SSE4
#endif

   IMPLICIT NONE

   ! Except for some error handling code, all code should
   ! get a unit number from the print keys or from the logger, in order
   ! to guarantee correct output behavior,
   ! for example in farming or path integral runs
   ! default_input_unit should never be used
   ! but we need to know what it is, as we should not try to open it for output
   INTEGER, PUBLIC, PARAMETER                   :: default_output_unit = output_unit, &
                                                   default_input_unit = input_unit

#include "machine_cpuid.h"
   ! Enumerates the target architectures or instruction set extensions.
   ! A feature is present if within range for the respective architecture.
   ! For example, to check for MACHINE_X86_AVX the following is true:
   ! MACHINE_X86_AVX <= m_cpuid() and MACHINE_X86 >= m_cpuid().
   ! For example, to check for MACHINE_ARM_SOME the following is true:
   ! MACHINE_ARM_SOME <= m_cpuid() and MACHINE_ARM >= m_cpuid().
   INTEGER, PUBLIC, PARAMETER :: &
      MACHINE_CPU_GENERIC = CP_MACHINE_CPU_GENERIC, &
      MACHINE_X86_SSE4 = CP_MACHINE_X86_SSE4, &
      MACHINE_X86_AVX = CP_MACHINE_X86_AVX, &
      MACHINE_X86_AVX2 = CP_MACHINE_X86_AVX2, &
      MACHINE_X86_AVX512 = CP_MACHINE_X86_AVX512, &
      MACHINE_X86 = MACHINE_X86_AVX512 ! marks end of range
   ! other arch to be added as needed e.g.,
   !MACHINE_ARM_SOME    = 2000
   !MACHINE_ARM_ELSE    = 2001
   !MACHINE_ARM = MACHINE_ARM_ELSE
   !MACHINE_PWR_????    = 3000

   PRIVATE

   PUBLIC :: m_walltime, m_datum, m_hostnm, m_flush, &
             m_getcwd, m_getlog, m_getpid, m_procrun, m_abort, &
             m_chdir, m_mov, m_memory, m_memory_details, m_energy, &
             m_cpuinfo, m_cpuid_static, m_cpuid, m_cpuid_name, &
             m_get_omp_stacksize

   INTERFACE
      ! **********************************************************************************************
      !> \brief Target architecture or instruction set extension according to compiler target flags.
      !> \return cpuid according to MACHINE_* integer-parameter.
      !> \par History
      !>      04.2019 created [Hans Pabst]
      ! **********************************************************************************************
      PURE FUNCTION m_cpuid_static() BIND(C)
         IMPORT :: C_INT
         INTEGER(C_INT) :: m_cpuid_static
      END FUNCTION m_cpuid_static
   END INTERFACE

   ! Flushing is enabled by default because without it crash reports can get lost.
   ! For performance reasons it can be disabled via the input in &GLOBAL.
   LOGICAL, SAVE, PUBLIC :: flush_should_flush = .TRUE.

   INTEGER(KIND=int_8), SAVE, PUBLIC :: m_memory_max = 0

CONTAINS

! **************************************************************************************************
!> \brief flushes units if the &GLOBAL flag is set accordingly
!> \param lunit ...
!> \par History
!>      10.2008 created [Joost VandeVondele]
!> \note
!>      flushing might degrade performance significantly (30% and more)
! **************************************************************************************************
   SUBROUTINE m_flush(lunit)
      INTEGER, INTENT(IN)                                :: lunit

      IF (flush_should_flush) FLUSH (lunit)

   END SUBROUTINE

! **************************************************************************************************
!> \brief returns time from a real-time clock, protected against rolling
!>      early/easily
!> \return ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
!> \note
!>      same implementation for all machines.
!>      might still roll, if not called multiple times per count_max/count_rate
! **************************************************************************************************
   FUNCTION m_walltime() RESULT(wt)
      REAL(KIND=dp)                                      :: wt

#if defined(__LIBXSMM)
      wt = libxsmm_timer_duration(0_int_8, libxsmm_timer_tick())
#else
      wt = omp_get_wtime()
#endif
   END FUNCTION m_walltime

! **************************************************************************************************
!> \brief reads /proc/cpuinfo if it exists (i.e. Linux) to return relevant info
!> \param model_name as obtained from the 'model name' field, UNKNOWN otherwise
! **************************************************************************************************
   SUBROUTINE m_cpuinfo(model_name)
      CHARACTER(LEN=default_string_length)               :: model_name

      INTEGER, PARAMETER                                 :: bufferlen = 2048

      CHARACTER(LEN=bufferlen)                           :: buffer
      INTEGER                                            :: i, icol, iline, imod, stat

      model_name = "UNKNOWN"
      buffer = ""
      OPEN (121245, FILE="/proc/cpuinfo", ACTION="READ", STATUS="OLD", ACCESS="STREAM", IOSTAT=stat)
      IF (stat == 0) THEN
         DO i = 1, bufferlen
            READ (121245, END=999) buffer(I:I)
         ENDDO
999      CLOSE (121245)
         imod = INDEX(buffer, "model name")
         IF (imod > 0) THEN
            icol = imod - 1 + INDEX(buffer(imod:), ":")
            iline = icol - 1 + INDEX(buffer(icol:), NEW_LINE('A'))
            IF (iline == icol - 1) iline = bufferlen + 1
            model_name = buffer(icol + 1:iline - 1)
         ENDIF
      ENDIF
   END SUBROUTINE m_cpuinfo

! **************************************************************************************************
!> \brief Target architecture or instruction set extension according to CPU-check at runtime.
!> \return cpuid according to MACHINE_* integer-parameter.
!> \par History
!>      04.2019 created [Hans Pabst]
! **************************************************************************************************
   PURE FUNCTION m_cpuid() RESULT(cpuid)
      INTEGER :: cpuid
#if defined(__LIBXSMM)
      cpuid = libxsmm_get_target_archid()
      cpuid = MERGE(MIN(MACHINE_X86_SSE4 + cpuid - LIBXSMM_X86_SSE4, MACHINE_X86), &
                    MACHINE_CPU_GENERIC, LIBXSMM_X86_SSE4 .LE. cpuid)
#else
      cpuid = m_cpuid_static()
#endif
   END FUNCTION m_cpuid

! **************************************************************************************************
!> \brief Determine name of target architecture for a given CPUID.
!> \param cpuid integer value (MACHINE_*)
!> \return name or short name.
!> \par History
!>      06.2019 created [Hans Pabst]
! **************************************************************************************************
   FUNCTION m_cpuid_name(cpuid)
      INTEGER                                            :: cpuid
      CHARACTER(len=default_string_length), POINTER      :: m_cpuid_name

      CHARACTER(len=default_string_length), SAVE, TARGET :: name_generic = "generic", &
         name_unknown = "unknown", name_x86_avx = "x86_avx", name_x86_avx2 = "x86_avx2", &
         name_x86_avx512 = "x86_avx512", name_x86_sse4 = "x86_sse4"

      SELECT CASE (cpuid)
      CASE (MACHINE_CPU_GENERIC)
         m_cpuid_name => name_generic
      CASE (MACHINE_X86_SSE4)
         m_cpuid_name => name_x86_sse4
      CASE (MACHINE_X86_AVX)
         m_cpuid_name => name_x86_avx
      CASE (MACHINE_X86_AVX2)
         m_cpuid_name => name_x86_avx2
      CASE (MACHINE_X86_AVX512)
         m_cpuid_name => name_x86_avx512
      CASE DEFAULT
         m_cpuid_name => name_unknown
      END SELECT
   END FUNCTION m_cpuid_name

! **************************************************************************************************
!> \brief returns the energy used since some time in the past.
!>        The precise meaning depends on the infrastructure is available.
!>        In the cray_pm_energy case, this is the energy used by the node in kJ.
!> \return ...
!> \par History
!>      09.2013 created [Joost VandeVondele, Ole Schuett]
! **************************************************************************************************
   FUNCTION m_energy() RESULT(wt)
      REAL(KIND=dp)                            :: wt

#if defined(__CRAY_PM_ENERGY)
      wt = read_energy("/sys/cray/pm_counters/energy")
#elif defined(__CRAY_PM_ACCEL_ENERGY)
      wt = read_energy("/sys/cray/pm_counters/accel_energy")
#else
      wt = 0.0 ! fallback default
#endif

   END FUNCTION m_energy

#if defined(__CRAY_PM_ACCEL_ENERGY) || defined(__CRAY_PM_ENERGY)
! **************************************************************************************************
!> \brief reads energy values from the sys-filesystem
!> \param filename ...
!> \return ...
!> \par History
!>      09.2013 created [Joost VandeVondele, Ole Schuett]
! **************************************************************************************************
   FUNCTION read_energy(filename) RESULT(wt)
      CHARACTER(LEN=*)                                   :: filename
      REAL(KIND=dp)                                      :: wt

      CHARACTER(LEN=80)                                  :: DATA
      INTEGER                                            :: i, iostat
      INTEGER(KIND=int_8)                                :: raw

      OPEN (121245, FILE=filename, ACTION="READ", STATUS="OLD", ACCESS="STREAM")
      DO I = 1, 80
         READ (121245, END=999) DATA(I:I)
      ENDDO
999   CLOSE (121245)
      DATA(I:80) = ""
      READ (DATA, *, IOSTAT=iostat) raw
      IF (iostat .NE. 0) THEN
         wt = 0.0_dp
      ELSE
         ! convert from J to kJ
         wt = raw/1000.0_dp
      ENDIF
   END FUNCTION read_energy
#endif

! **************************************************************************************************
!> \brief returns a datum in human readable format using a standard Fortran routine
!> \param cal_date ...
!> \par History
!>      10.2009 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE m_datum(cal_date)
      CHARACTER(len=*), INTENT(OUT)                      :: cal_date

      CHARACTER(len=10)                                  :: time
      CHARACTER(len=8)                                   :: date

      CALL DATE_AND_TIME(date=date, time=time)
      cal_date = date(1:4)//"-"//date(5:6)//"-"//date(7:8)//" "//time(1:2)//":"//time(3:4)//":"//time(5:10)

   END SUBROUTINE m_datum

! **************************************************************************************************
!> \brief Can be used to get a nice core
! **************************************************************************************************
   SUBROUTINE m_abort()
      INTERFACE
         SUBROUTINE abort() BIND(C, name="abort")
         END SUBROUTINE
      END INTERFACE

      CALL abort()
   END SUBROUTINE m_abort

! **************************************************************************************************
!> \brief Returns if a process is running on the local machine
!>        1 if yes and 0 if not
!> \param pid ...
!> \return ...
! **************************************************************************************************
   FUNCTION m_procrun(pid) RESULT(run_on)
      INTEGER, INTENT(IN)       ::   pid
      INTEGER                   ::   run_on
#if defined(__MINGW)
      run_on = 0
#else
      INTEGER                   ::   istat

      INTERFACE
         FUNCTION kill(pid, sig) BIND(C, name="kill") RESULT(errno)
            IMPORT
            INTEGER(KIND=C_INT), VALUE                :: pid, sig
            INTEGER(KIND=C_INT)                      :: errno
         END FUNCTION
      END INTERFACE

      ! If sig is 0, then no signal is sent, but error checking is still
      ! performed; this can be used to check for the existence of a process
      ! ID or process group ID.

      istat = kill(pid=pid, sig=0)
      IF (istat == 0) THEN
         run_on = 1 ! no error, process exists
      ELSE
         run_on = 0 ! error, process probably does not exist
      ENDIF
#endif
   END FUNCTION m_procrun

! **************************************************************************************************
!> \brief Returns the total amount of memory [bytes] in use, if known, zero otherwise
!> \param mem ...
! **************************************************************************************************
   SUBROUTINE m_memory(mem)

      INTEGER(KIND=int_8), OPTIONAL, INTENT(OUT)         :: mem
      INTEGER(KIND=int_8)                      :: mem_local

      !
      ! __NO_STATM_ACCESS can be used to disable the stuff, if getpagesize
      ! lead to linking errors or /proc/self/statm can not be opened
      !
#if defined(__NO_STATM_ACCESS)
      mem_local = 0
#else
      INTEGER(KIND=int_8)                      :: m1, m2, m3
      CHARACTER(LEN=80) :: DATA
      INTEGER :: iostat, i

      ! the size of a page, might not be available everywhere
      INTERFACE
         FUNCTION getpagesize() BIND(C, name="getpagesize") RESULT(RES)
            IMPORT
            INTEGER(C_INT) :: RES
         END FUNCTION
      END INTERFACE

      !
      ! reading from statm
      !
      mem_local = -1
      DATA = ""
      OPEN (121245, FILE="/proc/self/statm", ACTION="READ", STATUS="OLD", ACCESS="STREAM")
      DO I = 1, 80
         READ (121245, END=999) DATA(I:I)
      ENDDO
999   CLOSE (121245)
      DATA(I:80) = ""
      ! m1 = total
      ! m2 = resident
      ! m3 = shared
      READ (DATA, *, IOSTAT=iostat) m1, m2, m3
      IF (iostat .NE. 0) THEN
         mem_local = 0
      ELSE
         mem_local = m2
#if defined(__STATM_TOTAL)
         mem_local = m1
#endif
#if defined(__STATM_RESIDENT)
         mem_local = m2
#endif
         mem_local = mem_local*getpagesize()
      ENDIF
#endif

      m_memory_max = MAX(mem_local, m_memory_max)
      IF (PRESENT(mem)) mem = mem_local

   END SUBROUTINE m_memory

! **************************************************************************************************
!> \brief get more detailed memory info, all units are bytes.
!>         the only 'useful' option is MemLikelyFree which is an estimate of remaining memory
!>         assumed to give info like /proc/meminfo while MeMLikelyFree is the amount of
!>         memory we're likely to be able to allocate, but not necessarily in one chunk
!>         zero means not available...
!> \param MemTotal ...
!> \param MemFree ...
!> \param Buffers ...
!> \param Cached ...
!> \param Slab ...
!> \param SReclaimable ...
!> \param MemLikelyFree ...
! **************************************************************************************************
   SUBROUTINE m_memory_details(MemTotal, MemFree, Buffers, Cached, Slab, SReclaimable, MemLikelyFree)

      INTEGER(kind=int_8), OPTIONAL :: MemTotal, MemFree, Buffers, Cached, Slab, SReclaimable, MemLikelyFree

      INTEGER, PARAMETER :: Nbuffer = 10000
      CHARACTER(LEN=Nbuffer) :: meminfo

      INTEGER :: i

      MemTotal = 0
      MemFree = 0
      Buffers = 0
      Cached = 0
      Slab = 0
      SReclaimable = 0
      MemLikelyFree = 0
      meminfo = ""

      OPEN (UNIT=8123, file="/proc/meminfo", ACCESS="STREAM", ERR=901)
      i = 0
      DO
         i = i + 1
         IF (i > Nbuffer) EXIT
         READ (8123, END=900, ERR=900) meminfo(i:i)
      ENDDO
900   CONTINUE
      meminfo(i:Nbuffer) = ""
901   CONTINUE
      CLOSE (8123, ERR=902)
902   CONTINUE
      MemTotal = get_field_value_in_bytes('MemTotal:')
      MemFree = get_field_value_in_bytes('MemFree:')
      Buffers = get_field_value_in_bytes('Buffers:')
      Cached = get_field_value_in_bytes('Cached:')
      Slab = get_field_value_in_bytes('Slab:')
      SReclaimable = get_field_value_in_bytes('SReclaimable:')
      ! opinions here vary but this might work
      MemLikelyFree = MemFree + Buffers + Cached + SReclaimable

   CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param field ...
!> \return ...
! **************************************************************************************************
      INTEGER(int_8) FUNCTION get_field_value_in_bytes(field)
      CHARACTER(LEN=*)                                   :: field

      INTEGER                                            :: start
      INTEGER(KIND=int_8)                                :: value

         get_field_value_in_bytes = 0
         start = INDEX(meminfo, field)
         IF (start .NE. 0) THEN
            start = start + LEN_TRIM(field)
            IF (start .LT. Nbuffer) THEN
               READ (meminfo(start:), *, ERR=999, END=999) value
               ! XXXXXXX convert from Kb to bytes XXXXXXXX
               get_field_value_in_bytes = value*1024
999            CONTINUE
            ENDIF
         ENDIF
      END FUNCTION
   END SUBROUTINE m_memory_details

! **************************************************************************************************
!> \brief ...
!> \param hname ...
! **************************************************************************************************
   SUBROUTINE m_hostnm(hname)
      CHARACTER(len=*), INTENT(OUT)            :: hname
#if defined(__MINGW)
      ! While there is a gethostname in the Windows (POSIX) API, it requires that winsocks is
      ! initialised prior to using it via WSAStartup(..), respectively cleaned up at the end via WSACleanup().
      hname = "<unknown>"
#else
      INTEGER                                  :: istat, i
      CHARACTER(len=default_path_length)       :: buf

      INTERFACE
         FUNCTION gethostname(buf, buflen) BIND(C, name="gethostname") RESULT(errno)
            IMPORT
            CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: buf
            INTEGER(KIND=C_INT), VALUE               :: buflen
            INTEGER(KIND=C_INT)                      :: errno
         END FUNCTION
      END INTERFACE

      istat = gethostname(buf, LEN(buf))
      IF (istat /= 0) THEN
         WRITE (*, *) "m_hostnm failed"
         CALL m_abort()
      ENDIF
      i = INDEX(buf, c_null_char) - 1
      hname = buf(1:i)
#endif
   END SUBROUTINE m_hostnm

! **************************************************************************************************
!> \brief ...
!> \param curdir ...
! **************************************************************************************************
   SUBROUTINE m_getcwd(curdir)
      CHARACTER(len=*), INTENT(OUT)            :: curdir
      TYPE(C_PTR)                              :: stat
      INTEGER                                  :: i
      CHARACTER(len=default_path_length), TARGET  :: tmp

      INTERFACE
         FUNCTION getcwd(buf, buflen) BIND(C, name="getcwd") RESULT(stat)
            IMPORT
            CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: buf
            INTEGER(KIND=C_INT), VALUE               :: buflen
            TYPE(C_PTR)                              :: stat
         END FUNCTION
      END INTERFACE

      stat = getcwd(tmp, LEN(tmp))
      IF (.NOT. C_ASSOCIATED(stat)) THEN
         WRITE (*, *) "m_getcwd failed"
         CALL m_abort()
      ENDIF
      i = INDEX(tmp, c_null_char) - 1
      curdir = tmp(1:i)
   END SUBROUTINE m_getcwd

! **************************************************************************************************
!> \brief ...
!> \param dir ...
!> \param ierror ...
! **************************************************************************************************
   SUBROUTINE m_chdir(dir, ierror)
      CHARACTER(len=*), INTENT(IN)             :: dir
      INTEGER, INTENT(OUT)                     :: ierror

      INTERFACE
         FUNCTION chdir(path) BIND(C, name="chdir") RESULT(errno)
            IMPORT
            CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: path
            INTEGER(KIND=C_INT)                      :: errno
         END FUNCTION
      END INTERFACE

      ierror = chdir(TRIM(dir)//c_null_char)
   END SUBROUTINE m_chdir

! **************************************************************************************************
!> \brief ...
!> \param pid ...
! **************************************************************************************************
   SUBROUTINE m_getpid(pid)
      INTEGER, INTENT(OUT)                     :: pid

      INTERFACE
         FUNCTION getpid() BIND(C, name="getpid") RESULT(pid)
            IMPORT
            INTEGER(KIND=C_INT)              :: pid
         END FUNCTION
      END INTERFACE

      pid = getpid()
   END SUBROUTINE m_getpid

! **************************************************************************************************
!> \brief ...
!> \param path ...
!> \return ...
! **************************************************************************************************
   FUNCTION m_unlink(path) RESULT(istat)

      CHARACTER(LEN=*), INTENT(IN)             :: path

      INTEGER                                  :: istat

      INTERFACE
         FUNCTION unlink(path) BIND(C, name="unlink") RESULT(errno)
            IMPORT
            CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: path
            INTEGER(KIND=C_INT)                      :: errno
         END FUNCTION
      END INTERFACE

      istat = unlink(TRIM(path)//c_null_char)
   END FUNCTION m_unlink

! **************************************************************************************************
!> \brief ...
!> \param source ...
!> \param TARGET ...
! **************************************************************************************************
   SUBROUTINE m_mov(source, TARGET)

      CHARACTER(LEN=*), INTENT(IN)             :: source, TARGET

      INTEGER                                  :: istat

      INTERFACE
         FUNCTION rename(src, dest) BIND(C, name="rename") RESULT(errno)
            IMPORT
            CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: src, dest
            INTEGER(KIND=C_INT)                      :: errno
         END FUNCTION
      END INTERFACE

      IF (TARGET == source) THEN
         WRITE (*, *) "Warning: m_mov ", TRIM(TARGET), " equals ", TRIM(source)
         RETURN
      ENDIF

      ! first remove target (needed on windows / mingw)
      istat = m_unlink(TARGET)
      ! ignore istat of unlink

      ! now move
      istat = rename(TRIM(source)//c_null_char, TRIM(TARGET)//c_null_char)
      IF (istat .NE. 0) THEN
         WRITE (*, *) "Trying to move "//TRIM(source)//" to "//TRIM(TARGET)//"."
         WRITE (*, *) "rename returned status: ", istat
         WRITE (*, *) "Problem moving file"
         CALL m_abort()
      ENDIF
   END SUBROUTINE m_mov

! **************************************************************************************************
!> \brief ...
!> \param user ...
! **************************************************************************************************
   SUBROUTINE m_getlog(user)

      CHARACTER(LEN=*), INTENT(OUT) :: user

      INTEGER                       :: istat

      ! on a posix system LOGNAME should be defined
      CALL get_environment_variable("LOGNAME", value=user, status=istat)
      ! nope, check alternative
      IF (istat /= 0) &
         CALL get_environment_variable("USER", value=user, status=istat)
      ! nope, check alternative
      IF (istat /= 0) &
         CALL get_environment_variable("USERNAME", value=user, status=istat)
      ! fall back
      IF (istat /= 0) &
         user = "<unknown>"

   END SUBROUTINE m_getlog

! **************************************************************************************************
!> \brief Retrieve environment variable OMP_STACKSIZE
!> \param omp_stacksize Value of OMP_STACKSIZE
! **************************************************************************************************
   SUBROUTINE m_get_omp_stacksize(omp_stacksize)
      CHARACTER(LEN=*), INTENT(OUT)                      :: omp_stacksize

      INTEGER                                            :: istat

      omp_stacksize = ""
      CALL get_environment_variable("OMP_STACKSIZE", value=omp_stacksize, status=istat)
      ! Fall back, if OMP_STACKSIZE is not set
      IF (istat /= 0) omp_stacksize = "default"

   END SUBROUTINE m_get_omp_stacksize

END MODULE machine
