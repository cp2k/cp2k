!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief dielectric constant data type
!> \par History
!>       06.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE dielectric_types

  USE kinds,                           ONLY: dp
  USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: pw_p_type,&
                                             pw_release,&
                                             pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE 

  PRIVATE 

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dielectric_types'
 
  TYPE dielectric_parameters
    LOGICAL                            :: dielec_core_correction
    REAL(KIND=dp)                      :: eps0
    REAL(KIND=dp)                      :: rho_min
    REAL(KIND=dp)                      :: rho_max
    INTEGER                            :: derivative_method
  END TYPE dielectric_parameters

  TYPE dielectric_type
     TYPE(dielectric_parameters)       :: params
     TYPE(pw_type), POINTER            :: eps
     TYPE(pw_type), POINTER            :: deps_drho
     TYPE(pw_p_type), DIMENSION(3)     :: dln_eps
  END TYPE dielectric_type

  PUBLIC :: dielectric_type, dielectric_parameters
  PUBLIC :: dielectric_release

  INTEGER, PARAMETER, PUBLIC           :: derivative_cd3 = 1, &
                                          derivative_cd5 = 2, &
                                          derivative_cd7 = 3, &
                                          derivative_fft = 100, &
                                          derivative_fft_use_deps = 101, &
                                          derivative_fft_use_drho = 102
CONTAINS

! *****************************************************************************
!> \brief   deallocates dielectric data type
!> \param dielectric the dielectric data type to be released
!> \param pw_pool pool of the plane wave grid  
!> \param error cp2k error
!>
!> \par History
!>       06.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dielectric_release(dielectric, pw_pool, error)
    TYPE(dielectric_type), POINTER           :: dielectric
    TYPE(pw_pool_type), INTENT(IN), &
      OPTIONAL, POINTER                      :: pw_pool
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dielectric_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    LOGICAL                                  :: can_give_back

    CALL timeset(routineN,handle)
    IF (ASSOCIATED(dielectric)) THEN
       can_give_back = PRESENT(pw_pool)
       IF (can_give_back) can_give_back = ASSOCIATED(pw_pool)
       IF (can_give_back) THEN
          CALL pw_pool_give_back_pw(pw_pool, dielectric%eps,&
                                    accept_non_compatible=.TRUE., error=error)
          CALL pw_pool_give_back_pw(pw_pool, dielectric%deps_drho,&
                                    accept_non_compatible=.TRUE., error=error)
          DO i = 1, 3
             CALL pw_pool_give_back_pw(pw_pool, dielectric%dln_eps(i)%pw,&
                                    accept_non_compatible=.TRUE., error=error)
          END DO
       ELSE
          CALL pw_release(dielectric%eps, error = error)
          CALL pw_release(dielectric%deps_drho, error = error)
          DO i = 1, 3
             CALL pw_release(dielectric%dln_eps(i)%pw, error=error)
          END DO
       END IF
       DEALLOCATE(dielectric)
    END IF

    CALL timestop(handle)

  END SUBROUTINE dielectric_release

END MODULE dielectric_types
