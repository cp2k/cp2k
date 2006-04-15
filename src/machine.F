!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000,2004  CP2K developers group                            !
!-----------------------------------------------------------------------------!
!!****** cp2k/machine [1.0] *
!!
!!   NAME
!!     machine
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     APSI & JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (05.07.2001) : added G95 interface
!!     - m_flush added (12.06.2002,MK)
!!     - Missing print_memory added (24.09.2002,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE machine
  USE  kinds,  ONLY        : dp
#if defined ( __AIX )
  USE machine_aix, ONLY    : & 
       m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __IRIX )
  USE machine_irix, ONLY   : &
       m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __G95 )
  USE machine_g95, ONLY    : &
       m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __INTEL )
  USE machine_intel, ONLY  : &
       m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __PGI )
  USE machine_pgi, ONLY    : &
       m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __XT3 )
  USE machine_xt3, ONLY    : &
       m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __T3E )
  USE machine_t3e, ONLY    : &
       m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __DEC )
  USE machine_dec, ONLY    : &
       m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __ABSOFT )
  USE machine_absoft, ONLY : &
       m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __SUN )
  USE machine_sun, ONLY    : &
       m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __NAG )
  USE machine_nag, ONLY    : &
       m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __SX )
  USE machine_sx, ONLY    : &
       m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: m_walltime, m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c

CONTAINS

!!****f* machine/m_walltime *
!!
!!   NAME
!!     m_walltime
!!
!!   FUNCTION
!!     returns time from a real-time clock, protected against rolling 
!!     early/easily 
!!
!!   NOTES
!!     same implementation for all machines. 
!!     might still roll, if not called multiple times per count_max/count_rate
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
FUNCTION m_walltime() RESULT (wt)
    REAL(KIND=dp)                            :: wt

    INTEGER                                  :: count, count_max, count_rate
    INTEGER, SAVE                            :: cycles = 0, last_count = -1

    ! count lies in [0,count_max] and increases monotonically
    CALL SYSTEM_CLOCK(count=count,count_rate=count_rate,count_max=count_max)

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

END FUNCTION m_walltime

END MODULE machine
!!*****
