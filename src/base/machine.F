!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      JGH (05.07.2001) : added G95 interface
!>      - m_flush added (12.06.2002,MK)
!>      - Missing print_memory added (24.09.2002,MK)
!> \author APSI & JGH
! **************************************************************************************************
MODULE machine
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
   USE ISO_FORTRAN_ENV, ONLY: input_unit, &
                              output_unit
   USE kinds, ONLY: default_string_length, &
                    dp, &
                    int_8
   USE machine_internal, ONLY: &
      m_abort, m_chdir, m_flush_internal => m_flush, m_getcwd, m_getlog, m_getpid, &
      m_hostnm, m_memory, m_memory_details, m_memory_max, m_mov, m_procrun

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads, OMP_GET_WTIME

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

   PUBLIC :: m_walltime, m_datum, m_hostnm, m_flush, m_flush_internal, &
             m_getcwd, m_getlog, m_getpid, m_procrun, m_abort, &
             m_chdir, m_mov, m_memory, m_memory_details, m_memory_max, m_energy, &
             m_cpuinfo, m_cpuid_static, m_cpuid, m_cpuid_name

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

   ! should only be set according to the state in &GLOBAL
   LOGICAL, SAVE, PUBLIC :: flush_should_flush = .FALSE.

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

      IF (flush_should_flush) CALL m_flush_internal(lunit)
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
#if defined(__LIBXSMM)
      USE libxsmm, ONLY: libxsmm_timer_tick, libxsmm_timer_duration
#endif
      REAL(KIND=dp)                                      :: wt

#if defined(__LIBXSMM)
      wt = libxsmm_timer_duration(0_int_8, libxsmm_timer_tick())
#else

      INTEGER(KIND=int_8)                                :: count
      INTEGER(KIND=int_8), SAVE                          :: count_max, count_rate, cycles = -1, &
                                                            last_count
!$    IF (.FALSE.) THEN
! count lies in [0,count_max] and increases monotonically

         IF (cycles == -1) THEN ! get parameters of system_clock and initialise
            CALL SYSTEM_CLOCK(count_rate=count_rate, count_max=count_max)
            cycles = 0
            last_count = 0
         ENDIF

         CALL SYSTEM_CLOCK(count=count)

         ! protect against non-standard cases where time might be non-monotonous,
         ! but it is unlikely that the clock cycled (e.g. underlying system clock adjustments)
         ! i.e. if count is smaller than last_count by only a small fraction of count_max,
         ! we use last_count instead
         ! if count is smaller, we assume that the clock cycled.
         IF (count < last_count) THEN
            IF (last_count - count < count_max/100) THEN
               count = last_count
            ELSE
               cycles = cycles + 1
            ENDIF
         ENDIF

         ! keep track of our history
         last_count = count

         wt = (REAL(count, KIND=dp) + REAL(cycles, KIND=dp)*(1.0_dp + REAL(count_max, KIND=dp))) &
              /REAL(count_rate, KIND=dp)
!$    ELSE
!$       wt = OMP_GET_WTIME()
!$    ENDIF
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
#if defined(__LIBXSMM)
      USE libxsmm, ONLY: libxsmm_get_target_archid, LIBXSMM_X86_SSE4
#endif
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

END MODULE machine
