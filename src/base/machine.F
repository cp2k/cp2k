!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH (05.07.2001) : added G95 interface
!>      - m_flush added (12.06.2002,MK)
!>      - Missing print_memory added (24.09.2002,MK)
!> \author APSI & JGH
! *****************************************************************************
MODULE machine
  USE ISO_FORTRAN_ENV,                 ONLY: input_unit,&
                                             output_unit
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE machine_internal,                ONLY: &
       m_abort, m_chdir, m_flush_internal=>m_flush, m_getarg, m_getcwd, &
       m_getlog, m_getpid, m_hostnm, m_iargc, m_memory, m_memory_details, &
       m_memory_max, m_mov, m_procrun

 !$ USE OMP_LIB

  IMPLICIT NONE

  ! Except for some error handling code, all code should
  ! get a unit number from the print keys or from the logger, in order
  ! to guarantee correct output behavior,
  ! for example in farming or path integral runs
  ! default_input_unit should never be used
  ! but we need to know what it is, as we should not try to open it for output
  INTEGER, PUBLIC, PARAMETER                   :: default_output_unit = output_unit, &
                                                  default_input_unit  = input_unit

  PRIVATE

  PUBLIC :: m_walltime, m_datum, m_flush, m_flush_internal,&
            m_hostnm, m_getcwd, m_getlog, m_getpid, m_getarg, m_procrun,&
            m_memory, m_iargc, m_abort, m_chdir, m_mov, m_memory_details, m_energy, m_memory_max

  ! should only be set according to the state in &GLOBAL
  LOGICAL, SAVE, PUBLIC :: flush_should_flush=.FALSE.

CONTAINS


! *****************************************************************************
!> \brief flushes units if the &GLOBAL flag is set accordingly
!> \param lunit ...
!> \par History
!>      10.2008 created [Joost VandeVondele]
!> \note
!>      flushing might degrade performance significantly (30% and more)
! *****************************************************************************
SUBROUTINE m_flush(lunit)
    INTEGER, INTENT(IN)                      :: lunit

   IF (flush_should_flush) CALL m_flush_internal(lunit)
END SUBROUTINE
! *****************************************************************************
!> \brief returns time from a real-time clock, protected against rolling
!>      early/easily
!> \retval wt ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
!> \note
!>      same implementation for all machines.
!>      might still roll, if not called multiple times per count_max/count_rate
! *****************************************************************************
FUNCTION m_walltime() RESULT (wt)
    REAL(KIND=dp)                            :: wt

    INTEGER(KIND=int_8)                      :: count
    INTEGER(KIND=int_8), SAVE                :: count_max, count_rate, &
                                                cycles = -1, last_count

    !$ IF (.FALSE.) THEN
! count lies in [0,count_max] and increases monotonically

    IF (cycles == -1) THEN ! get parameters of system_clock and initialise
        CALL SYSTEM_CLOCK(count_rate=count_rate,count_max=count_max)
        cycles = 0
        last_count = 0
    ENDIF

    CALL SYSTEM_CLOCK(count=count)

    ! protect against non-standard cases where time might be non-monotonous,
    ! but it is unlikely that the clock cycled (e.g. underlying system clock adjustments)
    ! i.e. if count is smaller than last_count by only a small fraction of count_max,
    ! we use last_count instead
    ! if count is smaller, we assume that the clock cycled.
    IF (count<last_count) THEN
       IF ( last_count-count < count_max / 100 ) THEN
          count=last_count
       ELSE
          cycles=cycles+1
       ENDIF
    ENDIF

    ! keep track of our history
    last_count=count

    wt = ( REAL(count,KIND=dp)+REAL(cycles,KIND=dp)*(1.0_dp+REAL(count_max,KIND=dp)) ) &
         / REAL(count_rate,KIND=dp)
    !$ ELSE
    !$    wt = OMP_GET_WTIME ()
    !$ ENDIF
END FUNCTION m_walltime


! *****************************************************************************
!> \brief returns the energy used since some time in the past.
!>        The precise meaning depends on the infrastructure is available.
!>        In the cray_pm_energy case, this is the energy used by the node in kJ.
!> \retval wt ...
!> \par History
!>      09.2013 created [Joost VandeVondele, Ole Schuett]
! *****************************************************************************
FUNCTION m_energy() RESULT (wt)
    REAL(KIND=dp)                            :: wt

    wt = 0.0 ! fallback default

#if defined(__CRAY_PM_ENERGY)
   wt = read_energy("/sys/cray/pm_counters/energy")
#elif defined(__CRAY_PM_ACCEL_ENERGY)
   wt = read_energy("/sys/cray/pm_counters/accel_energy")
#elif defined(__CRAY_PM_FAKE_ENERGY)
   ! reads from the current working dir
   wt = read_energy("energy")
#endif
END FUNCTION m_energy


! *****************************************************************************
!> \brief reads energy values from the sys-filesystem
!> \param filename ...
!> \retval wt ...
!> \par History
!>      09.2013 created [Joost VandeVondele, Ole Schuett]
! *****************************************************************************
FUNCTION read_energy(filename) RESULT (wt)
    CHARACTER(LEN=*)                         :: filename
    REAL(KIND=dp)                            :: wt

    CHARACTER(LEN=80)                        :: DATA
    INTEGER                                  :: i, iostat
    INTEGER(KIND=int_8)                      :: raw

    OPEN(121245,FILE=filename,ACTION="READ",STATUS="OLD",ACCESS="STREAM")
    DO I=1,80
       READ(121245,END=999) DATA(I:I)
    ENDDO
999 CLOSE(121245)
    DATA(I:80)=""
    READ(DATA,*,IOSTAT=iostat) raw
    IF (iostat.NE.0) THEN
       wt=0.0_dp
    ELSE
       ! convert from J to kJ
       wt=raw/1000.0_dp
    ENDIF
END FUNCTION read_energy


! *****************************************************************************
!> \brief returns a datum in human readable format using a standard Fortran routine
!> \param cal_date ...
!> \par History
!>      10.2009 created [Joost VandeVondele]
! *****************************************************************************
SUBROUTINE m_datum(cal_date)
    CHARACTER(len=*), INTENT(OUT)            :: cal_date

    CHARACTER(len=10)                        :: time
    CHARACTER(len=8)                         :: date

    CALL DATE_AND_TIME(date=date, time=time)
    cal_date=date(1:4)//"-"//date(5:6)//"-"//date(7:8)//" "//time(1:2)//":"//time(3:4)//":"//time(5:10)

END SUBROUTINE m_datum

END MODULE machine
