!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH (05.07.2001) : added G95 interface
!>      - m_flush added (12.06.2002,MK)
!>      - Missing print_memory added (24.09.2002,MK)
!> \author APSI & JGH
! *****************************************************************************
MODULE machine
  USE  kinds,  ONLY        : dp
#if defined ( __AIX )
  USE machine_aix, ONLY    : & 
       m_cputime, m_datum, m_flush_internal=>m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       m_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, m_mov
#elif defined ( __IRIX )
  USE machine_irix, ONLY   : &
       m_cputime, m_datum, m_flush_internal=>m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       m_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, m_mov
#elif defined ( __G95 )
  USE machine_g95, ONLY    : &
       m_cputime, m_datum, m_flush_internal=>m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       m_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, m_mov
#elif defined ( __GFORTRAN )
  USE machine_gfortran, ONLY    : &
       m_cputime, m_datum, m_flush_internal=>m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       m_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, m_mov
#elif defined ( __INTEL )
  USE machine_intel, ONLY  : &
       m_cputime, m_datum, m_flush_internal=>m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       m_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, m_mov
! must be before __PGI because pgf90 defines __PGI
#elif defined ( __XT3 )
  USE machine_xt3, ONLY    : &
       m_cputime, m_datum, m_flush_internal=>m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       m_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, m_mov
#elif defined ( __XT5 )
  USE machine_xt5, ONLY    : &
       m_cputime, m_datum, m_flush_internal=>m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       m_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, m_mov
#elif defined ( __PGI )
  USE machine_pgi, ONLY    : &
       m_cputime, m_datum, m_flush_internal=>m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       m_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, m_mov
#elif defined ( __T3E )
  USE machine_t3e, ONLY    : &
       m_cputime, m_datum, m_flush_internal=>m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       m_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, m_mov
#elif defined ( __DEC )
  USE machine_dec, ONLY    : &
       m_cputime, m_datum, m_flush_internal=>m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       m_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, m_mov
#elif defined ( __ABSOFT )
  USE machine_absoft, ONLY : &
       m_cputime, m_datum, m_flush_internal=>m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       m_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, m_mov
#elif defined ( __SUN )
  USE machine_sun, ONLY    : &
       m_cputime, m_datum, m_flush_internal=>m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       m_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, m_mov
#elif defined ( __NAG )
  USE machine_nag, ONLY    : &
       m_cputime, m_datum, m_flush_internal=>m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       m_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, m_mov
#elif defined ( __SX )
  USE machine_sx, ONLY    : &
       m_cputime, m_datum, m_flush_internal=>m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       m_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, m_mov
#endif

  IMPLICIT NONE

  ! this is the unit that is normally preconnected to stdout.
  ! this is unit 6 for all compilers we know, but in principle could be
  ! compiler dependent.
  ! however, except for some error handling code, all code should
  ! get a unit number from the print keys or from the logger, in order
  ! to guarantee correct output behavior,
  ! for example in farming or path integral runs
  ! default_input_unit should never be used
  ! but we need to know what it is, as we should not try to open it for output
  INTEGER, PUBLIC, PARAMETER                   :: default_output_unit = 6, &
                                                  default_input_unit  = 5

  PRIVATE

  PUBLIC :: m_walltime, m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       m_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c, m_mov

  ! should only be set according to the state in &GLOBAL
  LOGICAL, SAVE, PUBLIC :: flush_should_flush=.FALSE.

CONTAINS
! *****************************************************************************
!> \brief flushes units if the &GLOBAL flag is set accordingly
!> \note
!>      flushing might degrade performance significantly (30% and more)
!> \par History
!>      10.2008 created [Joost VandeVondele]
! *****************************************************************************
SUBROUTINE m_flush(lunit)
    INTEGER, INTENT(IN)                      :: lunit

   IF (flush_should_flush) CALL m_flush_internal(lunit)
END SUBROUTINE
! *****************************************************************************
!> \brief returns time from a real-time clock, protected against rolling 
!>      early/easily 
!> \note
!>      same implementation for all machines. 
!>      might still roll, if not called multiple times per count_max/count_rate
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
FUNCTION m_walltime() RESULT (wt)
    REAL(KIND=dp)                            :: wt

    INTEGER                                  :: count
    INTEGER, SAVE                            :: count_max, count_rate, &
                                                cycles = -1, last_count

!$OMP CRITICAL
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
!$OMP END CRITICAL
END FUNCTION m_walltime

END MODULE machine
