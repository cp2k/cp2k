MODULE dbcsr_machine

#if defined ( __NAG )
  USE f90_unix
#endif

  USE dbcsr_kinds,                     ONLY: dp

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: default_output_unit=6

  PUBLIC :: m_walltime, m_flush, m_abort

CONTAINS

  SUBROUTINE m_flush(lunit)
    INTEGER, INTENT(IN)                      :: lunit

#if defined ( __G95 )
    CALL flush(lunit)
#elif defined ( __GFORTRAN )
    CALL flush(lunit)
#elif defined ( __INTEL )
    CALL flush(lunit)
#elif defined ( __XT5 )
    CALL flush(lunit)
#elif defined ( __PGI )
#elif defined ( __NAG )
    CALL flush(lunit)

#endif
  END SUBROUTINE m_flush

  SUBROUTINE m_abort()
#if defined ( __G95 )
    CALL abort()
#elif defined ( __GFORTRAN )
    CALL abort()
#elif defined ( __INTEL )
    CALL abort()
#elif defined ( __XT5 )
    CALL abort()
#elif defined ( __PGI )
    CALL abort()
#elif defined ( __NAG )
    CALL abort()
#endif
  END SUBROUTINE m_abort

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

END MODULE dbcsr_machine
