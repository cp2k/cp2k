!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  numerical operations on real-space grid
!> \par History
!>       12.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE rs_methods

  USE kinds,                           ONLY: dp
  USE pw_types,                        ONLY: pw_p_type,&
                                             pw_type
  USE realspace_grid_types,            ONLY: &
       pw2rs, realspace_grid_desc_type, realspace_grid_p_type, &
       realspace_grid_type, rs2pw, rs_grid_create, rs_grid_release, &
       rs_grid_zero, rs_pw_transfer
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rs_methods'

  PUBLIC derive_fdm_cd3,&
         derive_fdm_cd5,&
         derive_fdm_cd7

CONTAINS 

! *****************************************************************************
!> \brief    2nd order finite difference derivative of a function on realspace grid
!> \param f  input function
!> \param df derivative of f
!> \param rs_grid real-space grid
!> \param error cp2k error
!> \par History:
!>      - Creation (15.11.2013,MK)
!>      - Refactored and moved here from qs_sccs.F (12.2014, Hossein Bani-Hashemian)
!> \author     Matthias Krack (MK)
!> \version    1.0
! *****************************************************************************
  SUBROUTINE derive_fdm_cd3(f,df,rs_grid,error)

    TYPE(pw_type), POINTER                   :: f
    TYPE(pw_p_type), DIMENSION(3), &
      INTENT(OUT)                            :: df
    TYPE(realspace_grid_type), POINTER       :: rs_grid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'derive_fdm_cd3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k
    INTEGER, DIMENSION(3)                    :: lb, ub
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: h
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: drdx, drdy, drdz, r
    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    TYPE(realspace_grid_p_type), &
      DIMENSION(3)                           :: drs_grid

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(f),cp_failure_level,routineP,error,failure)

    ! Setup
    rs_desc => rs_grid%desc
    CALL rs_pw_transfer(rs_grid,f,pw2rs,error=error)
    DO i=1,3
       CALL rs_grid_create(drs_grid(i)%rs_grid,rs_desc,error=error)
       CALL rs_grid_zero(drs_grid(i)%rs_grid)
    END DO

    lb(1:3) = rs_grid%lb_real(1:3)
    ub(1:3) = rs_grid%ub_real(1:3)
    r => rs_grid%r
    drdx => drs_grid(1)%rs_grid%r
    drdy => drs_grid(2)%rs_grid%r
    drdz => drs_grid(3)%rs_grid%r

    ! 3-point stencil central differences
    h(1:3) = 2.0_dp*f%pw_grid%dr(1:3)
    !$omp parallel do default(none) &
    !$omp             private(i,j,k) &
    !$omp             shared(drdx,drdy,drdz,h,lb,r,ub)
    DO k=lb(3),ub(3)
      DO j=lb(2),ub(2)
        DO i=lb(1),ub(1)
          drdx(i,j,k) = (r(i+1,j,k) - r(i-1,j,k))/h(1)
          drdy(i,j,k) = (r(i,j+1,k) - r(i,j-1,k))/h(2)
          drdz(i,j,k) = (r(i,j,k+1) - r(i,j,k-1))/h(3)
        END DO
      END DO
    END DO
    !$omp end parallel do

    ! Cleanup
    DO i=1,3
       CALL rs_pw_transfer(drs_grid(i)%rs_grid, df(i)%pw,rs2pw,error=error)
       CALL rs_grid_release(drs_grid(i)%rs_grid,error=error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE derive_fdm_cd3

! *****************************************************************************
!> \brief    4th order finite difference derivative of a function on realspace grid
!> \param f  input function
!> \param df derivative of f
!> \param rs_grid real-space grid
!> \param error cp2k error
!> \par History:
!>      - Creation (15.11.2013,MK)
!>      - Refactored and moved here from qs_sccs.F (12.2014, Hossein Bani-Hashemian)
!> \author     Matthias Krack (MK)
!> \version    1.0
! *****************************************************************************
  SUBROUTINE derive_fdm_cd5(f,df,rs_grid,error)

    TYPE(pw_type), POINTER                   :: f
    TYPE(pw_p_type), DIMENSION(3), &
      INTENT(OUT)                            :: df
    TYPE(realspace_grid_type), POINTER       :: rs_grid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'derive_fdm_cd5', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k
    INTEGER, DIMENSION(3)                    :: lb, ub
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: h
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: drdx, drdy, drdz, r
    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    TYPE(realspace_grid_p_type), &
      DIMENSION(3)                           :: drs_grid

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(f),cp_failure_level,routineP,error,failure)

    ! Setup
    rs_desc => rs_grid%desc
    CALL rs_pw_transfer(rs_grid,f,pw2rs,error=error)
    DO i=1,3
       CALL rs_grid_create(drs_grid(i)%rs_grid,rs_desc,error=error)
       CALL rs_grid_zero(drs_grid(i)%rs_grid)
    END DO

    lb(1:3) = rs_grid%lb_real(1:3)
    ub(1:3) = rs_grid%ub_real(1:3)
    r => rs_grid%r
    drdx => drs_grid(1)%rs_grid%r
    drdy => drs_grid(2)%rs_grid%r
    drdz => drs_grid(3)%rs_grid%r

    ! 5-point stencil central differences
    h(1:3) = 12.0_dp*f%pw_grid%dr(1:3)
    !$omp parallel do default(none) &
    !$omp             private(i,j,k) &
    !$omp             shared(drdx,drdy,drdz,h,lb,r,ub)
    DO k=lb(3),ub(3)
      DO j=lb(2),ub(2)
        DO i=lb(1),ub(1)
          drdx(i,j,k) = (r(i-2,j,k) - r(i+2,j,k) + 8.0_dp*(r(i+1,j,k) - r(i-1,j,k)))/h(1)
          drdy(i,j,k) = (r(i,j-2,k) - r(i,j+2,k) + 8.0_dp*(r(i,j+1,k) - r(i,j-1,k)))/h(2)
          drdz(i,j,k) = (r(i,j,k-2) - r(i,j,k+2) + 8.0_dp*(r(i,j,k+1) - r(i,j,k-1)))/h(3)
        END DO
      END DO
    END DO
    !$omp end parallel do

    ! Cleanup
    DO i=1,3
       CALL rs_pw_transfer(drs_grid(i)%rs_grid, df(i)%pw,rs2pw,error=error)
       CALL rs_grid_release(drs_grid(i)%rs_grid,error=error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE derive_fdm_cd5

! *****************************************************************************
!> \brief    6th order finite difference derivative of a function on realspace grid
!> \param f  input function
!> \param df derivative of f
!> \param rs_grid real-space grid
!> \param error cp2k error
!> \par History:
!>      - Creation (15.11.2013,MK)
!>      - Refactored and moved here from qs_sccs.F (12.2014, Hossein Bani-Hashemian)
!> \author     Matthias Krack (MK)
!> \version    1.0
! *****************************************************************************
  SUBROUTINE derive_fdm_cd7(f,df,rs_grid,error)

    TYPE(pw_type), POINTER                   :: f
    TYPE(pw_p_type), DIMENSION(3), &
      INTENT(OUT)                            :: df
    TYPE(realspace_grid_type), POINTER       :: rs_grid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'derive_fdm_cd7', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k
    INTEGER, DIMENSION(3)                    :: lb, ub
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: h
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: drdx, drdy, drdz, r
    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    TYPE(realspace_grid_p_type), &
      DIMENSION(3)                           :: drs_grid

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(f),cp_failure_level,routineP,error,failure)

    ! Setup
    rs_desc => rs_grid%desc
    CALL rs_pw_transfer(rs_grid,f,pw2rs,error=error)
    DO i=1,3
       CALL rs_grid_create(drs_grid(i)%rs_grid,rs_desc,error=error)
       CALL rs_grid_zero(drs_grid(i)%rs_grid)
    END DO

    lb(1:3) = rs_grid%lb_real(1:3)
    ub(1:3) = rs_grid%ub_real(1:3)
    r => rs_grid%r
    drdx => drs_grid(1)%rs_grid%r
    drdy => drs_grid(2)%rs_grid%r
    drdz => drs_grid(3)%rs_grid%r

    ! 7-point stencil central differences
    h(1:3) = 60.0_dp*f%pw_grid%dr(1:3)
    !$omp parallel do default(none) &
    !$omp             private(i,j,k) &
    !$omp             shared(drdx,drdy,drdz,h,lb,r,ub)
    DO k=lb(3),ub(3)
      DO j=lb(2),ub(2)
        DO i=lb(1),ub(1)
          drdx(i,j,k) = (r(i+3,j,k) - r(i-3,j,k) + 9.0_dp*(r(i-2,j,k) - r(i+2,j,k)) +&
                         45.0_dp*(r(i+1,j,k) - r(i-1,j,k)))/h(1)
          drdy(i,j,k) = (r(i,j+3,k) - r(i,j-3,k) + 9.0_dp*(r(i,j-2,k) - r(i,j+2,k)) +&
                         45.0_dp*(r(i,j+1,k) - r(i,j-1,k)))/h(2)
          drdz(i,j,k) = (r(i,j,k+3) - r(i,j,k-3) + 9.0_dp*(r(i,j,k-2) - r(i,j,k+2)) +&
                         45.0_dp*(r(i,j,k+1) - r(i,j,k-1)))/h(3)
        END DO
      END DO
    END DO
    !$omp end parallel do

    ! Cleanup
    DO i=1,3
       CALL rs_pw_transfer(drs_grid(i)%rs_grid, df(i)%pw,rs2pw,error=error)
       CALL rs_grid_release(drs_grid(i)%rs_grid,error=error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE derive_fdm_cd7

END MODULE rs_methods
