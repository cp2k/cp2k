!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2003  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pw_spline_utils [1.0] *
!!
!!   NAME
!!     pw_spline_utils
!!
!!   FUNCTION
!!     different utils that are useful to manipulate splines on the regular
!!     grid of a pw
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2003 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE pw_spline_utils
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE message_passing,                 ONLY: mp_sum
  USE pw_grid_types,                   ONLY: FULLSPACE
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pw_spline_utils'

    REAL ( KIND = dp ), PUBLIC, PARAMETER, DIMENSION(4) :: nn10_coeffs=&
         (/ 125._dp/216._dp, 25._dp/432._dp, 5._dp/864._dp, 1._dp/1728._dp /)
    ! factor 10 between subsequent values

  PUBLIC :: pw_spline3_interpolate_values_g, pw_spline3_evaluate_values_g,&
       pw_spline3_deriv_r, pw_spline3_deriv_g, pw_spline3_test_coeff
  PUBLIC :: pw_spline_scale_deriv, pw_nn_smear_g
  PUBLIC :: pw_spline2_interpolate_values_g, pw_spline2_evaluate_values_g,&
       pw_spline2_deriv_r, pw_spline2_deriv_g, pw_spline2_test_coeff
!***
!****************************************************************************
CONTAINS

!!****f* pw_spline_utils/pw_spline2_interpolate_values_g [1.0] *
!!
!!   NAME
!!     pw_spline2_interpolate_values_g
!!
!!   SYNOPSIS
!!     Subroutine pw_spline2_interpolate_values_g(spline_g, error)
!!       Type(pw_type), Pointer:: spline_g
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_spline2_interpolate_values_g
!!
!!   FUNCTION
!!     calculates the FFT of the coefficents of the quadratic spline that
!!     interpolates the given values
!!
!!   NOTES
!!     does not work with spherical cutoff
!!
!!   INPUTS
!!     - spline_g: on entry the FFT of the values to interpolate as cc,
!!       will contain the FFT of the coefficents of the spline
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_spline2_interpolate_values_g(spline_g,error)
    TYPE(pw_type),POINTER :: spline_g
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline2_interpolate_values_g',&
         routineP=moduleN//':'//routineN
    REAL (KIND = dp) :: c23, c_val, coeff
    REAL (KIND = dp), DIMENSION(:), ALLOCATABLE :: cosIVals,cosJVals,cosKVals
    INTEGER, DIMENSION(3) :: n_tot
    INTEGER :: ii,i, j, k, handle, stat
    INTEGER, DIMENSION (2,3) :: gbo
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    n_tot(1:3) = spline_g%pw_grid%npts (1:3)
    gbo = spline_g%pw_grid%bounds
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(spline_g),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_use==COMPLEXDATA1D,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_space==RECIPROCALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.spline_g%pw_grid%spherical,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%pw_grid%grid_span==FULLSPACE,cp_failure_level,routineP,error,failure)

    ALLOCATE(cosIVals(gbo(1,1):gbo(2,1)),cosJVals(gbo(1,2):gbo(2,2)),&
         cosKVals(gbo(1,3):gbo(2,3)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    coeff=twopi/n_tot(1)
    !$omp parallel do private(i)
    DO i=gbo(1,1),gbo(2,1)
       cosIVals(i)=COS(coeff*REAL(i,dp))
    END DO
    coeff=twopi/n_tot(2)
    !$omp parallel do private(j)
    DO j=gbo(1,2),gbo(2,2)
       cosJVals(j)=COS(coeff*REAL(j,dp))
    END DO
    coeff=twopi/n_tot(3)
    !$omp parallel do private(k)
    DO k=gbo(1,3),gbo(2,3)
       cosKVals(k)=COS(coeff*REAL(k,dp))
    END DO

    !$omp parallel do private(i,j,k,ii,coeff,c23)
    DO ii=1,SIZE(spline_g%cc)
       i=spline_g%pw_grid%g_hat(1,ii)
       j=spline_g%pw_grid%g_hat(2,ii)
       k=spline_g%pw_grid%g_hat(3,ii)

       c23=cosJVals(j)*cosKVals(k)
       coeff=64.0_dp/(cosIVals(i)*c23+&
            (cosIVals(i)*cosJVals(j)+cosIVals(i)*cosKVals(k)+c23)*3.0_dp+&
            (cosIVals(i)+cosJVals(j)+cosKVals(k))*9.0_dp+&
            27.0_dp)
       
       spline_g%cc(ii)=spline_g%cc(ii)*coeff
       
    END DO
    DEALLOCATE(cosIVals, cosJVals, cosKVals, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE pw_spline2_interpolate_values_g
!***************************************************************************

!!****f* pw_spline_utils/pw_spline3_interpolate_values_g [1.0] *
!!
!!   NAME
!!     pw_spline3_interpolate_values_g
!!
!!   SYNOPSIS
!!     Subroutine pw_spline3_interpolate_values_g(spline_g, error)
!!       Type(pw_type), Pointer:: spline_g
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_spline3_interpolate_values_g
!!
!!   FUNCTION
!!     calculates the FFT of the coefficents of the cubic spline that
!!     interpolates the given values
!!
!!   NOTES
!!     does not work with spherical cutoff
!!     stupid distribution for cos calculation, it should calculate only the
!!     needed cos, and avoid the mp_sum
!!
!!   INPUTS
!!     - spline_g: on entry the FFT of the values to interpolate as cc,
!!       will contain the FFT of the coefficents of the spline
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_spline3_interpolate_values_g(spline_g,error)
    TYPE(pw_type),POINTER :: spline_g
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline3_interpolate_values_g',&
         routineP=moduleN//':'//routineN
    REAL (KIND = dp) :: c23, c_val, coeff
    REAL (KIND = dp), DIMENSION(:), ALLOCATABLE :: cosIVals,cosJVals,cosKVals
    INTEGER, DIMENSION(3) :: n_tot
    INTEGER :: ii,i, j, k, handle, stat
    INTEGER, DIMENSION (2,3) :: gbo
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    n_tot(1:3) = spline_g%pw_grid%npts (1:3)
    gbo = spline_g%pw_grid%bounds
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(spline_g),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_use==COMPLEXDATA1D,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_space==RECIPROCALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.spline_g%pw_grid%spherical,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%pw_grid%grid_span==FULLSPACE,cp_failure_level,routineP,error,failure)

    ALLOCATE(cosIVals(gbo(1,1):gbo(2,1)),&
         cosJVals(gbo(1,2):gbo(2,2)),&
         cosKVals(gbo(1,3):gbo(2,3)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    coeff=twopi/n_tot(1)
    !$omp parallel do private(i)
    DO i=gbo(1,1),gbo(2,1)
       cosIVals(i)=COS(coeff*REAL(i,dp))
    END DO
    coeff=twopi/n_tot(2)
    !$omp parallel do private(j)
    DO j=gbo(1,2),gbo(2,2)
       cosJVals(j)=COS(coeff*REAL(j,dp))
    END DO
    coeff=twopi/n_tot(3)
    !$omp parallel do private(k)
    DO k=gbo(1,3),gbo(2,3)
       cosKVals(k)=COS(coeff*REAL(k,dp))
    END DO

    !$omp parallel do private(i,j,k,ii,coeff,c23)
    DO ii=1,SIZE(spline_g%cc)
       i=spline_g%pw_grid%g_hat(1,ii)
       j=spline_g%pw_grid%g_hat(2,ii)
       k=spline_g%pw_grid%g_hat(3,ii)
             ! no opt
!FM                coeff=1.0/((cosVal(1)*cosVal(2)*cosVal(3))/27.0_dp+&
!FM                     (cosVal(1)*cosVal(2)+cosVal(1)*cosVal(3)+&
!FM                     cosVal(2)*cosVal(3))*2.0_dp/27.0_dp+&
!FM                     (cosVal(1)+cosVal(2)+cosVal(3))*4.0_dp/27.0_dp+&
!FM                     8.0_dp/27.0_dp)
             ! opt
       c23=cosJVals(j)*cosKVals(k)
       coeff=27.0_dp/(cosIVals(i)*c23+&
            (cosIVals(i)*cosJVals(j)+cosIVals(i)*cosKVals(k)+c23)*2.0_dp+&
            (cosIVals(i)+cosJVals(j)+cosKVals(k))*4.0_dp+&
            8.0_dp)
       
       spline_g%cc(ii)=spline_g%cc(ii)*coeff
       
    END DO
    DEALLOCATE(cosIVals, cosJVals, cosKVals, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE pw_spline3_interpolate_values_g
!***************************************************************************

!!****f* pw_spline_utils/pw_spline2_evaluate_values_g [1.0] *
!!
!!   NAME
!!     pw_spline2_evaluate_values_g
!!
!!   SYNOPSIS
!!     Subroutine pw_spline2_evaluate_values_g(spline_g, error)
!!       Type(pw_type), Pointer:: spline_g
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_spline2_evaluate_values_g
!!
!!   FUNCTION
!!     calculates the FFT of the values of the quadratic spline given by the
!!     FFT of the coefficents of the spline (inverse of 
!!     pw_spline_interpolate_values_g)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - spline_g: on entry the FFT of the coefficents of the spline as cr3d,
!!       will contain the FFT of the values of the spline
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_spline2_evaluate_values_g(spline_g,error)
    TYPE(pw_type),POINTER :: spline_g
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline2_evaluate_values_g',&
         routineP=moduleN//':'//routineN
    REAL (KIND = dp) :: c23, c_val, coeff, inv64
    REAL (KIND = dp), DIMENSION(:), ALLOCATABLE :: cosIVals,cosJVals, cosKVals
    INTEGER, DIMENSION(3) :: n,n_tot
    INTEGER :: ii,i, j, k, handle, stat
    INTEGER, DIMENSION (2,3) :: bo, gbo
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    n(1:3) = spline_g%pw_grid%npts_local (1:3)
    n_tot(1:3) = spline_g%pw_grid%npts (1:3)
    bo = spline_g%pw_grid%bounds_local
    gbo = spline_g%pw_grid%bounds
    inv64=1.0_dp/64.0_dp

    CPPrecondition(ASSOCIATED(spline_g),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_use==COMPLEXDATA1D,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_space==RECIPROCALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.spline_g%pw_grid%spherical,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%pw_grid%grid_span==FULLSPACE,cp_failure_level,routineP,error,failure)

    ALLOCATE(cosIVals(gbo(1,1):gbo(2,1)),cosJVals(gbo(1,2):gbo(2,2)), &
         cosKVals(gbo(1,3):gbo(2,3)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    coeff=twopi/n_tot(1)
    !$omp parallel do private(i)
    DO i=gbo(1,1),gbo(2,1)
       cosIVals(i)=COS(coeff*REAL(i,dp))
    END DO
    coeff=twopi/n_tot(2)
    !$omp parallel do private(j)
    DO j=gbo(1,2),gbo(2,2)
       cosJVals(j)=COS(coeff*REAL(j,dp))
    END DO
    !$omp parallel do private(k)
    DO k=gbo(1,3),gbo(2,3)
       cosKVals(k)=COS(coeff*REAL(k,dp))
    END DO

    !$omp parallel do private(ii,k,j,i,coeff,c23)
    DO ii=1,SIZE(spline_g%cc)
       i=spline_g%pw_grid%g_hat(1,ii)
       j=spline_g%pw_grid%g_hat(2,ii)
       k=spline_g%pw_grid%g_hat(3,ii)

       c23=cosJVals(j)*cosKVals(k)
       coeff=(cosIVals(i)*c23+&
            (cosIVals(i)*cosJVals(j)+cosIVals(i)*cosKVals(k)+c23)*3.0_dp+&
            (cosIVals(i)+cosJVals(j)+cosKVals(k))*9.0_dp+&
            27.0_dp)*inv64
       
       spline_g%cc(ii)=spline_g%cc(ii)*coeff
    END DO
    DEALLOCATE(cosIVals, cosJVals, cosKVals, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE pw_spline2_evaluate_values_g
!***************************************************************************

!!****f* pw_spline_utils/pw_nn_smear_g [1.0] *
!!
!!   NAME
!!     pw_nn_smear_g
!!
!!   SYNOPSIS
!!     Subroutine pw_nn_smear_g(spline_g, coeffs, error)
!!       Type(pw_type), Pointer:: spline_g
!!       Real(Kind=dp), Dimension(4), Intent (IN):: coeffs
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_nn_smear_g
!!
!!   FUNCTION
!!     calculates the FFT of the values of a nearest neighbor smearing
!!
!!   NOTES
!!     does not normalize the smear to 1.
!!     with coeff=(/ 8._dp/27._dp, 2._dp/27._dp, 1._dp/54._dp, 1._dp/216._dp /)
!!     is equivalent to pw_spline3_evaluate_values_g, with
!!     coeff=(/ 27._dp/64._dp, 9._dp/128._dp, 3._dp/256._dp, 1._dp/512._dp /)
!!
!!   INPUTS
!!     - spline_g: on entry the FFT of the values to interpolate as cc3d,
!!       will contain the FFT of the smeared values
!!     - coeffs: array with the coefficent of the smearing, ordered with 
!!       the distance from the center: coeffs(1) the coeff of the central
!!       element, coeffs(2) the coeff of the 6 element with distance 1,
!!       coeff(3) the coeff of the 12 elements at distance sqrt(2),
!!       coeff(4) the coeff of the 8 elements at distance sqrt(3).
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_nn_smear_g(spline_g,coeffs,error)
    TYPE(pw_type),POINTER :: spline_g
    REAL(KIND = dp), DIMENSION(4), INTENT(in) :: coeffs
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_nn_smear_g',&
         routineP=moduleN//':'//routineN
    REAL (KIND = dp) :: c23, c_val, coeff
    REAL (KIND = dp), DIMENSION(4) :: r_coeffs
    REAL (KIND = dp), DIMENSION(:), ALLOCATABLE :: cosIVals,cosJVals, cosKVals
    INTEGER, DIMENSION(3) :: n,n_tot
    INTEGER :: ii,i, j, k, handle, stat
    INTEGER, DIMENSION (2,3) :: bo, gbo
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    n(1:3) = spline_g%pw_grid%npts_local (1:3)
    n_tot(1:3) = spline_g%pw_grid%npts (1:3)
    bo = spline_g%pw_grid%bounds_local
    gbo = spline_g%pw_grid%bounds

    r_coeffs=coeffs
    r_coeffs(2)=r_coeffs(2)*2.0_dp
    r_coeffs(3)=r_coeffs(3)*4.0_dp
    r_coeffs(4)=r_coeffs(4)*8.0_dp

    CPPrecondition(ASSOCIATED(spline_g),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_use==COMPLEXDATA1D,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_space==RECIPROCALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.spline_g%pw_grid%spherical,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%pw_grid%grid_span==FULLSPACE,cp_failure_level,routineP,error,failure)

    ALLOCATE(cosIVals(gbo(1,1):gbo(2,1)),cosJVals(gbo(1,2):gbo(2,2)), &
         cosKVals(gbo(1,3):gbo(2,3)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    coeff=twopi/n_tot(1)
    !$omp parallel do private(i)
    DO i=gbo(1,1),gbo(2,1)
       cosIVals(i)=COS(coeff*REAL(i,dp))
    END DO
    coeff=twopi/n_tot(2)
    !$omp parallel do private(j)
    DO j=gbo(1,2),gbo(2,2)
       cosJVals(j)=COS(coeff*REAL(j,dp))
    END DO
    !$omp parallel do private(k)
    DO k=gbo(1,3),gbo(2,3)
       cosKVals(k)=COS(coeff*REAL(k,dp))
    END DO

    !$omp parallel do private(ii,k,j,i,coeff,c23)
    DO ii=1,SIZE(spline_g%cc)
       i=spline_g%pw_grid%g_hat(1,ii)
       j=spline_g%pw_grid%g_hat(2,ii)
       k=spline_g%pw_grid%g_hat(3,ii)

       c23=cosJVals(j)*cosKVals(k)
       coeff=(r_coeffs(4)*cosIVals(i)*c23+&
            (cosIVals(i)*cosJVals(j)+cosIVals(i)*cosKVals(k)+c23)*r_coeffs(3)+&
            (cosIVals(i)+cosJVals(j)+cosKVals(k))*r_coeffs(2)+&
            r_coeffs(1))
       
       spline_g%cc(ii)=spline_g%cc(ii)*coeff
    END DO
    DEALLOCATE(cosIVals, cosJVals, cosKVals, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE pw_nn_smear_g
!***************************************************************************

!!****f* pw_spline_utils/pw_spline3_evaluate_values_g [1.0] *
!!
!!   NAME
!!     pw_spline3_evaluate_values_g
!!
!!   SYNOPSIS
!!     Subroutine pw_spline3_evaluate_values_g(spline_g, error)
!!       Type(pw_type), Pointer:: spline_g
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_spline3_evaluate_values_g
!!
!!   FUNCTION
!!     calculates the FFT of the values of the cubic spline given by the
!!     FFT of the coefficents of the spline (inverse of 
!!     pw_spline_interpolate_values_g)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - spline_g: on entry the FFT of the values to interpolate as cc3d,
!!       will contain the FFT of the coefficents of the spline
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_spline3_evaluate_values_g(spline_g,error)
    TYPE(pw_type),POINTER :: spline_g
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline3_evaluate_values_g',&
         routineP=moduleN//':'//routineN
    REAL (KIND = dp) :: c23, c_val, coeff, inv27
    REAL (KIND = dp), DIMENSION(:), ALLOCATABLE :: cosIVals,cosJVals, cosKVals
    INTEGER, DIMENSION(3) :: n,n_tot
    INTEGER :: ii,i, j, k, handle, stat
    INTEGER, DIMENSION (2,3) :: bo, gbo
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    n(1:3) = spline_g%pw_grid%npts_local (1:3)
    n_tot(1:3) = spline_g%pw_grid%npts (1:3)
    bo = spline_g%pw_grid%bounds_local
    gbo = spline_g%pw_grid%bounds
    inv27=1.0_dp/27.0_dp

    CPPrecondition(ASSOCIATED(spline_g),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_use==COMPLEXDATA1D,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_space==RECIPROCALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.spline_g%pw_grid%spherical,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%pw_grid%grid_span==FULLSPACE,cp_failure_level,routineP,error,failure)

    ALLOCATE(cosIVals(gbo(1,1):gbo(2,1)),cosJVals(gbo(1,2):gbo(2,2)), &
         cosKVals(gbo(1,3):gbo(2,3)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    coeff=twopi/n_tot(1)
    !$omp parallel do private(i)
    DO i=gbo(1,1),gbo(2,1)
       cosIVals(i)=COS(coeff*REAL(i,dp))
    END DO
    coeff=twopi/n_tot(2)
    !$omp parallel do private(j)
    DO j=gbo(1,2),gbo(2,2)
       cosJVals(j)=COS(coeff*REAL(j,dp))
    END DO
    coeff=twopi/n_tot(3)
    !$omp parallel do private(k)
    DO k=gbo(1,3),gbo(2,3)
       cosKVals(k)=COS(coeff*REAL(k,dp))
    END DO

    !$omp parallel do private(ii,k,j,i,coeff,c23)
    DO ii=1,SIZE(spline_g%cc)
       i=spline_g%pw_grid%g_hat(1,ii)
       j=spline_g%pw_grid%g_hat(2,ii)
       k=spline_g%pw_grid%g_hat(3,ii)
             ! no opt
!FM                coeff=((cosVal(1)*cosVal(2)*cosVal(3))/27.0_dp+&
!FM                     (cosVal(1)*cosVal(2)+cosVal(1)*cosVal(3)+&
!FM                     cosVal(2)*cosVal(3))*2.0_dp/27.0_dp+&
!FM                     (cosVal(1)+cosVal(2)+cosVal(3))*4.0_dp/27.0_dp+&
!FM                     8.0_dp/27.0_dp)
             ! opt
       c23=cosJVals(j)*cosKVals(k)
       coeff=(cosIVals(i)*c23+&
            (cosIVals(i)*cosJVals(j)+cosIVals(i)*cosKVals(k)+c23)*2.0_dp+&
            (cosIVals(i)+cosJVals(j)+cosKVals(k))*4.0_dp+&
            8.0_dp)*inv27
       
       spline_g%cc(ii)=spline_g%cc(ii)*coeff
    END DO
    DEALLOCATE(cosIVals, cosJVals, cosKVals, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE pw_spline3_evaluate_values_g
!***************************************************************************

!!****f* xc/pw_spline3_deriv_r [1.0] *
!!
!!   NAME
!!     pw_spline3_deriv_r
!!
!!   SYNOPSIS
!!     Subroutine pw_spline3_deriv_r(spline_coeff_r, deriv_val_r, idir,&
!!         error)
!!       Type(pw_type), Pointer:: spline_coeff_r, deriv_val_r
!!       Integer, Intent (IN):: idir
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_spline3_deriv_r
!!
!!   FUNCTION
!!     puts the x,y,z derivatives of the given cubic spline
!!
!!   NOTES
!!     non parallel, no multigrid test, derivatives taken assuming
!!     unit distance between gridpoints
!!
!!   ARGUMENTS
!!     - spline_coeff_r: the cofficents of the spline in real space
!!     - deriv_val_r: will contain the derivative of the spline
!!     - idir: direction of the derivative (1,2,3: x,y,z)
!!     - pw_pool: plane waves pool
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     05.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_spline3_deriv_r(spline_coeff_r,deriv_val_r, idir, error)
    TYPE(pw_type),POINTER :: spline_coeff_r, deriv_val_r
    INTEGER, INTENT(in) :: idir
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline3_deriv',&
         routineP=moduleN//':'//routineN
    REAL (KIND = dp), DIMENSION(3) :: cosVal
    REAL (KIND = dp), DIMENSION(3,3) :: deriv_coeff
    REAL (KIND = dp) :: c23, c_val, cosKVal
    REAL (KIND = dp), DIMENSION(:), ALLOCATABLE :: cosIVals,cosJVals
    INTEGER, DIMENSION(3) :: n,n_tot
    INTEGER :: ii, i, j, k, handle, stat, ij
    INTEGER, DIMENSION (2,3) :: bo
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    deriv_coeff = RESHAPE ((/ &
         1.0_dp/72.0_dp, 1.0_dp/18.0_dp, 1.0_dp/72.0_dp,&
         1.0_dp/18.0_dp, 2.0_dp/9.0_dp,  1.0_dp/18.0_dp,&
         1.0_dp/72.0_dp, 1.0_dp/18.0_dp, 1.0_dp/72.0_dp /),(/3,3/))
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(spline_coeff_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_coeff_r%in_space==REALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_coeff_r%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_val_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_val_r%in_space==REALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_val_r%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_coeff_r%pw_grid%para%group_size==0,cp_failure_level,routineP,error,failure)

    n(1:3) = spline_coeff_r%pw_grid%npts_local (1:3)
    n_tot(1:3) = spline_coeff_r%pw_grid%npts (1:3)
    bo = spline_coeff_r%pw_grid%bounds_local


    SELECT CASE(idir)
    CASE(1)
       ! x deriv
       !$omp parallel do private(k,j,i,ii,ij,c_val)
       DO k = 0,n(3)-1
          DO j = 0,n(2)-1
             DO i = 0,n(1)-1

                c_val=0.0_dp
                DO ii=-1,1
                   DO ij=-1,1
                      c_val=c_val+deriv_coeff(ii+2,ij+2)*&
                           (&
                           spline_coeff_r%cr3d &
                           (bo(1,1)+MODULO(i+1,n(1)),&
                           bo(1,2)+MODULO(j+ii,n(2)),&
                           bo(1,3)+MODULO(k+ij,n(3)))&
                           -&
                           spline_coeff_r%cr3d &
                           (bo(1,1)+MODULO(i-1,n(1)),&
                           bo(1,2)+MODULO(j+ii,n(2)),&
                           bo(1,3)+MODULO(k+ij,n(3)))&
                           )
                   END DO
                END DO
                deriv_val_r%cr3d&
                     (bo(1,1)+i,bo(1,2)+j,bo(1,3)+k)=c_val

             END DO
          END DO
       END DO
    CASE(2)
       ! y deriv
       !$omp parallel do private(k,j,i,ii,ij,c_val)
       DO k = 0,n(3)-1
          DO j = 0,n(2)-1
             DO i = 0,n(1)-1

                c_val=0.0_dp
                DO ii=-1,1
                   DO ij=-1,1
                      c_val=c_val+deriv_coeff(ii+2,ij+2)*&
                           (&
                           spline_coeff_r%cr3d &
                           (bo(1,1)+MODULO(i+ii,n(1)),&
                           bo(1,2)+MODULO(j+1,n(2)),&
                           bo(1,3)+MODULO(k+ij,n(3)))&
                           -&
                           spline_coeff_r%cr3d &
                           (bo(1,1)+MODULO(i+ii,n(1)),&
                           bo(1,2)+MODULO(j-1,n(2)),&
                           bo(1,3)+MODULO(k+ij,n(3)))&
                           )
                   END DO
                END DO
                deriv_val_r%cr3d&
                     (bo(1,1)+i,bo(1,2)+j,bo(1,3)+k)=c_val

             END DO
          END DO
       END DO
    CASE(3)
       ! z deriv
       !$omp parallel do private(k,j,i,ii,ij,c_val)
       DO k = 0,n(3)-1
          DO j = 0,n(2)-1
             DO i = 0,n(1)-1

                c_val=0.0_dp
                DO ii=-1,1
                   DO ij=-1,1
                      c_val=c_val+deriv_coeff(ii+2,ij+2)*&
                           (&
                           spline_coeff_r%cr3d &
                           (bo(1,1)+MODULO(i+ii,n(1)),&
                           bo(1,2)+MODULO(j+ij,n(2)),&
                           bo(1,3)+MODULO(k+1,n(3)))&
                           -&
                           spline_coeff_r%cr3d &
                           (bo(1,1)+MODULO(i+ii,n(1)),&
                           bo(1,2)+MODULO(j+ij,n(2)),&
                           bo(1,3)+MODULO(k-1,n(3)))&
                           )
                   END DO
                END DO
                deriv_val_r%cr3d&
                     (bo(1,1)+i,bo(1,2)+j,bo(1,3)+k)=c_val

             END DO
          END DO
       END DO
    CASE default
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
            routineP,"invalid idir="//cp_to_string(idir)//&
            CPSourceFileRef,&
            error=error,failure=failure)
    END SELECT

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE pw_spline3_deriv_r
!***************************************************************************

!!****f* xc/pw_spline2_deriv_r [1.0] *
!!
!!   NAME
!!     pw_spline2_deriv_r
!!
!!   SYNOPSIS
!!     Subroutine pw_spline2_deriv_r(spline_coeff_r, deriv_val_r, idir,&
!!         error)
!!       Type(pw_type), Pointer:: spline_coeff_r, deriv_val_r
!!       Integer, Intent (IN):: idir
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_spline2_deriv_r
!!
!!   FUNCTION
!!     puts the x,y,z derivatives of the given quadratic spline
!!
!!   NOTES
!!     non parallel, no multigrid test, derivatives taken assuming
!!     unit distance between gridpoints
!!
!!   ARGUMENTS
!!     - spline_coeff_r: the cofficents of the spline in real space
!!     - deriv_val_r: will contain the derivative of the spline
!!     - idir: direction of the derivative (1,2,3: x,y,z)
!!     - pw_pool: plane waves pool
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     05.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_spline2_deriv_r(spline_coeff_r,deriv_val_r, idir, error)
    TYPE(pw_type),POINTER :: spline_coeff_r, deriv_val_r
    INTEGER, INTENT(in) :: idir
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline2_deriv',&
         routineP=moduleN//':'//routineN
    REAL (KIND = dp), DIMENSION(3) :: cosVal
    REAL (KIND = dp), DIMENSION(3,3) :: deriv_coeff
    REAL (KIND = dp) :: c23, c_val, cosKVal
    REAL (KIND = dp), DIMENSION(:), ALLOCATABLE :: cosIVals,cosJVals
    INTEGER, DIMENSION(3) :: n,n_tot
    INTEGER :: ii, i, j, k, handle, stat, ij
    INTEGER, DIMENSION (2,3) :: bo
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    deriv_coeff = RESHAPE ((/ &
         1.0_dp/128.0_dp, 3.0_dp/64.0_dp, 1.0_dp/128.0_dp,&
         3.0_dp/64.0_dp, 9.0_dp/32.0_dp,  3.0_dp/64.0_dp,&
         1.0_dp/128.0_dp, 3.0_dp/64.0_dp, 1.0_dp/128.0_dp /),(/3,3/))
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(spline_coeff_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_coeff_r%in_space==REALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_coeff_r%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_val_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_val_r%in_space==REALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_val_r%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_coeff_r%pw_grid%para%group_size==0,cp_failure_level,routineP,error,failure)

    n(1:3) = spline_coeff_r%pw_grid%npts_local (1:3)
    n_tot(1:3) = spline_coeff_r%pw_grid%npts (1:3)
    bo = spline_coeff_r%pw_grid%bounds_local


    SELECT CASE(idir)
    CASE(1)
       ! x deriv
       !$omp parallel do private(k,j,i,ii,ij,c_val)
       DO k = 0,n(3)-1
          DO j = 0,n(2)-1
             DO i = 0,n(1)-1

                c_val=0.0_dp
                DO ii=-1,1
                   DO ij=-1,1
                      c_val=c_val+deriv_coeff(ii+2,ij+2)*&
                           (&
                           spline_coeff_r%cr3d &
                           (bo(1,1)+MODULO(i+1,n(1)),&
                           bo(1,2)+MODULO(j+ii,n(2)),&
                           bo(1,3)+MODULO(k+ij,n(3)))&
                           -&
                           spline_coeff_r%cr3d &
                           (bo(1,1)+MODULO(i-1,n(1)),&
                           bo(1,2)+MODULO(j+ii,n(2)),&
                           bo(1,3)+MODULO(k+ij,n(3)))&
                           )
                   END DO
                END DO
                deriv_val_r%cr3d&
                     (bo(1,1)+i,bo(1,2)+j,bo(1,3)+k)=c_val

             END DO
          END DO
       END DO
    CASE(2)
       ! y deriv
       !$omp parallel do private(k,j,i,ii,ij,c_val)
       DO k = 0,n(3)-1
          DO j = 0,n(2)-1
             DO i = 0,n(1)-1

                c_val=0.0_dp
                DO ii=-1,1
                   DO ij=-1,1
                      c_val=c_val+deriv_coeff(ii+2,ij+2)*&
                           (&
                           spline_coeff_r%cr3d &
                           (bo(1,1)+MODULO(i+ii,n(1)),&
                           bo(1,2)+MODULO(j+1,n(2)),&
                           bo(1,3)+MODULO(k+ij,n(3)))&
                           -&
                           spline_coeff_r%cr3d &
                           (bo(1,1)+MODULO(i+ii,n(1)),&
                           bo(1,2)+MODULO(j-1,n(2)),&
                           bo(1,3)+MODULO(k+ij,n(3)))&
                           )
                   END DO
                END DO
                deriv_val_r%cr3d&
                     (bo(1,1)+i,bo(1,2)+j,bo(1,3)+k)=c_val

             END DO
          END DO
       END DO
    CASE(3)
       ! z deriv
       !$omp parallel do private(k,j,i,ii,ij,c_val)
       DO k = 0,n(3)-1
          DO j = 0,n(2)-1
             DO i = 0,n(1)-1

                c_val=0.0_dp
                DO ii=-1,1
                   DO ij=-1,1
                      c_val=c_val+deriv_coeff(ii+2,ij+2)*&
                           (&
                           spline_coeff_r%cr3d &
                           (bo(1,1)+MODULO(i+ii,n(1)),&
                           bo(1,2)+MODULO(j+ij,n(2)),&
                           bo(1,3)+MODULO(k+1,n(3)))&
                           -&
                           spline_coeff_r%cr3d &
                           (bo(1,1)+MODULO(i+ii,n(1)),&
                           bo(1,2)+MODULO(j+ij,n(2)),&
                           bo(1,3)+MODULO(k-1,n(3)))&
                           )
                   END DO
                END DO
                deriv_val_r%cr3d&
                     (bo(1,1)+i,bo(1,2)+j,bo(1,3)+k)=c_val

             END DO
          END DO
       END DO
    CASE default
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
            routineP,"invalid idir="//cp_to_string(idir)//&
            CPSourceFileRef,&
            error=error,failure=failure)
    END SELECT

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE pw_spline2_deriv_r
!***************************************************************************

!!****f* pw_spline_utils/pw_spline_scale_deriv [1.0] *
!!
!!   NAME
!!     pw_spline_scale_deriv
!!
!!   SYNOPSIS
!!     Subroutine pw_spline_scale_deriv(deriv_vals_r, cell, transpose,&
!!         scale, error)
!!       Type(pw_p_type), Dimension(3):: deriv_vals_r
!!       Type(cell_type), Pointer:: cell
!!       Logical, Intent (IN), Optional:: transpose
!!       Real(Kind=dp), Intent (IN), Optional:: scale
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_spline_scale_deriv
!!
!!   FUNCTION
!!     rescales the derivatives from gridspacing=1 to the real derivatives
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - deriv_vals_r: an array of x,y,z derivatives
!!     - cell: the mapping to the real cell
!!     - transpose: if true applies the transpose of the map (defaults to 
!!       false)
!!     - scale: a scaling factor (defaults to 1.0)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_spline_scale_deriv(deriv_vals_r,cell,transpose,scale,error)
    TYPE(pw_p_type), DIMENSION(3) :: deriv_vals_r ! sun does not like that its components are modified
    TYPE(cell_type), POINTER :: cell
    LOGICAL, INTENT(in), OPTIONAL :: transpose
    REAL(KIND = dp), INTENT(in), OPTIONAL :: scale
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure, my_transpose
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline_scale_deriv',&
         routineP=moduleN//':'//routineN
    REAL (KIND = dp), DIMENSION(3) :: derivVal
    REAL (KIND = dp), DIMENSION(3,3) :: h_grid
    REAL(KIND = dp) :: my_scale
    INTEGER, DIMENSION(3) :: n,n_tot
    INTEGER :: i, j, k, handle, stat, idir
    INTEGER, DIMENSION (2,3) :: bo

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    my_transpose=.FALSE.
    IF (PRESENT(transpose)) my_transpose=transpose
    my_scale=1.0_dp
    IF (PRESENT(scale)) my_scale=scale
    n_tot(1:3) = deriv_vals_r(1)%pw%pw_grid%npts (1:3)
    bo = deriv_vals_r(1)%pw%pw_grid%bounds_local

    CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)
!FM    CPPrecondition(ASSOCIATED(deriv_vals_r),cp_failure_level,routineP,error,failure)

    ! map grid to real derivative
    IF (my_transpose) THEN
       DO j=1,3
          DO i=1,3
             h_grid(j,i)=my_scale*REAL(n_tot(i),dp)*cell%h_inv(i,j)      
          END DO
       END DO
    ELSE
       DO j=1,3
          DO i=1,3
             h_grid(i,j)=my_scale*REAL(n_tot(i),dp)*cell%h_inv(i,j)      
          END DO
       END DO
    END IF
    !$omp parallel do private(k,j,i,idir,derivVal)
    DO k = bo(1,3), bo(2,3)
       DO j = bo(1,2), bo(2,2)
          DO i = bo(1,1), bo(2,1)

             DO idir=1,3
                derivVal(idir)= &
                     h_grid(1,idir)*&
                     deriv_vals_r(1)%pw%cr3d(i,j,k)&
                     +h_grid(2,idir)*&
                     deriv_vals_r(2)%pw%cr3d(i,j,k)&
                     +h_grid(3,idir)*&
                     deriv_vals_r(3)%pw%cr3d(i,j,k)
             END DO
             DO idir=1,3
                deriv_vals_r(idir)%pw%cr3d(i,j,k)=&
                     derivVal(idir)
             END DO

          END DO
       END DO
    END DO

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE pw_spline_scale_deriv
!***************************************************************************

!!****f* pw_spline_utils/pw_spline3_deriv_g [1.0] *
!!
!!   NAME
!!     pw_spline3_deriv_g
!!
!!   SYNOPSIS
!!     Subroutine pw_spline3_deriv_g(spline_g, idir, error)
!!       Type(pw_type), Pointer:: spline_g
!!       Integer, Intent (IN):: idir
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_spline3_deriv_g
!!
!!   FUNCTION
!!     calculates the FFT of the values of the x,y,z (idir=1,2,3)
!!     derivative of the cubic spline
!!
!!   NOTES
!!     the distance between gridpoints is assumed to be 1
!!
!!   INPUTS
!!     - spline_g: on entry the FFT of the coefficents of the spline
!!       will contain the FFT of the derivative
!!     - idir: direction of the derivative
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_spline3_deriv_g(spline_g,idir,error)
    TYPE(pw_type),POINTER :: spline_g
    INTEGER, INTENT(in) :: idir
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline3_deriv_g',&
         routineP=moduleN//':'//routineN
    REAL (KIND = dp) :: tmp, c_val, coeff, inv9
    REAL (KIND = dp), DIMENSION(:), ALLOCATABLE :: csIVals,csJVals, csKVals
    INTEGER, DIMENSION(3) :: n,n_tot
    INTEGER :: ii,i, j, k, handle, stat
    INTEGER, DIMENSION (2,3) :: bo, gbo
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    n(1:3) = spline_g%pw_grid%npts_local (1:3)
    n_tot(1:3) = spline_g%pw_grid%npts (1:3)
    bo = spline_g%pw_grid%bounds_local
    gbo = spline_g%pw_grid%bounds
    inv9=1.0_dp/9.0_dp

    CPPrecondition(ASSOCIATED(spline_g),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_use==COMPLEXDATA1D,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_space==RECIPROCALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(.not.spline_g%pw_grid%spherical,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%pw_grid%grid_span==FULLSPACE,cp_failure_level,routineP,error,failure)

    ALLOCATE(csIVals(gbo(1,1):gbo(2,1)),&
         csJVals(gbo(1,2):gbo(2,2)),&
         csKVals(gbo(1,3):gbo(2,3)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    coeff=twopi/n_tot(1)
    IF (idir==1) THEN
       !$omp parallel do private(i)
       DO i=gbo(1,1),gbo(2,1)
          csIVals(i)=SIN(coeff*REAL(i,dp))
       END DO
    ELSE
       !$omp parallel do private(i)
       DO i=gbo(1,1),gbo(2,1)
          csIVals(i)=COS(coeff*REAL(i,dp))
       END DO
    END IF
    coeff=twopi/n_tot(2)
    IF (idir==2) THEN
       !$omp parallel do private(j)
       DO j=gbo(1,2),gbo(2,2)
          csJVals(j)=SIN(coeff*REAL(j,dp))
       END DO
    ELSE
       !$omp parallel do private(j)
       DO j=gbo(1,2),gbo(2,2)
          csJVals(j)=COS(coeff*REAL(j,dp))
       END DO
    END IF
    coeff=twopi/n_tot(3)
    IF (idir==3) THEN
       !$omp parallel do private(k)
       DO k=gbo(1,3),gbo(2,3)
          csKVals(k)=SIN(coeff*REAL(k,dp))
       END DO
    ELSE
       !$omp parallel do private(k)
       DO k=gbo(1,3),gbo(2,3)
          csKVals(k)=COS(coeff*REAL(k,dp))
       END DO
    END IF

    SELECT CASE(idir)
    CASE (1)
       ! x deriv
       !$omp parallel do private(ii,k,j,i,coeff,tmp)
       DO ii=1,SIZE(spline_g%cc)
          i=spline_g%pw_grid%g_hat(1,ii)
          j=spline_g%pw_grid%g_hat(2,ii)
          k=spline_g%pw_grid%g_hat(3,ii)
!FM                ! formula
!FM                coeff=(sinVal(1)*cosVal(2)*cosVal(3))/9.0_dp+&
!FM                     (sinVal(1)*cosVal(2)+sinVal(1)*cosVal(3))*2.0_dp/9.0_dp+&
!FM                     sinVal(1)*4.0_dp/9.0_dp
          tmp=csIVals(i)*csJVals(j)
          coeff=(tmp*csKVals(k)+&
               (tmp+csIVals(i)*csKVals(k))*2.0_dp+&
               csIVals(i)*4.0_dp)*inv9
          
          spline_g%cc(ii)=spline_g%cc(ii)*&
               CMPLX(0.0_dp,coeff,dp)
       END DO
    CASE(2)
       ! y deriv
       !$omp parallel do private(ii,k,j,i,coeff,tmp)
       DO ii=1,SIZE(spline_g%cc)
          i=spline_g%pw_grid%g_hat(1,ii)
          j=spline_g%pw_grid%g_hat(2,ii)
          k=spline_g%pw_grid%g_hat(3,ii)

          tmp=csIVals(i)*csJVals(j)
          coeff=(tmp*csKVals(k)+&
               (tmp+csJVals(j)*csKVals(k))*2.0_dp+&
               csJVals(j)*4.0_dp)*inv9
          
          spline_g%cc(ii)=spline_g%cc(ii)*&
               CMPLX(0.0_dp,coeff,dp)
       END DO
    CASE(3)
       ! z deriv
       !$omp parallel do private(ii,k,j,i,coeff,tmp)
       DO ii=1,SIZE(spline_g%cc)
          i=spline_g%pw_grid%g_hat(1,ii)
          j=spline_g%pw_grid%g_hat(2,ii)
          k=spline_g%pw_grid%g_hat(3,ii)

          tmp=csIVals(i)*csKVals(k)
          coeff=(tmp*csJVals(j)+&
               (tmp+csJVals(j)*csKVals(k))*2.0_dp+&
               csKVals(k)*4.0_dp)*inv9
          
          spline_g%cc(ii)=spline_g%cc(ii)*&
               CMPLX(0.0_dp,coeff,dp)
       END DO
    END SELECT

    DEALLOCATE(csIVals, csJVals, csKVals, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE pw_spline3_deriv_g
!***************************************************************************

!!****f* pw_spline_utils/pw_spline2_deriv_g [1.0] *
!!
!!   NAME
!!     pw_spline_deriv_g
!!
!!   SYNOPSIS
!!     Subroutine pw_spline_deriv_g(spline_g, idir, error)
!!       Type(pw_type), Pointer:: spline_g
!!       Integer, Intent (IN):: idir
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_spline_deriv_g
!!
!!   FUNCTION
!!     calculates the FFT of the values of the x,y,z (idir=1,2,3)
!!     derivative of the quadratic spline
!!
!!   NOTES
!!     the distance between gridpoints is assumed to be 1
!!
!!   INPUTS
!!     - spline_g: on entry the FFT of the coefficents of the spline
!!       will contain the FFT of the derivative
!!     - idir: direction of the derivative
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_spline2_deriv_g(spline_g,idir,error)
    TYPE(pw_type),POINTER :: spline_g
    INTEGER, INTENT(in) :: idir
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline2_deriv_g',&
         routineP=moduleN//':'//routineN
    REAL (KIND = dp) :: tmp, c_val, coeff, inv16
    REAL (KIND = dp), DIMENSION(:), ALLOCATABLE :: csIVals,csJVals, csKVals
    INTEGER, DIMENSION(3) :: n,n_tot
    INTEGER :: ii,i, j, k, handle, stat
    INTEGER, DIMENSION (2,3) :: bo
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    n(1:3) = spline_g%pw_grid%npts_local (1:3)
    n_tot(1:3) = spline_g%pw_grid%npts (1:3)
    bo = spline_g%pw_grid%bounds
    inv16=1.0_dp/16.0_dp

    CPPrecondition(ASSOCIATED(spline_g),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_use==COMPLEXDATA1D,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_space==RECIPROCALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(.not.spline_g%pw_grid%spherical,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%pw_grid%grid_span==FULLSPACE,cp_failure_level,routineP,error,failure)

    ALLOCATE(csIVals(bo(1,1):bo(2,1)),csJVals(bo(1,2):bo(2,2)),&
         csKVals(bo(1,3):bo(2,3)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    coeff=twopi/n_tot(1)
    IF (idir==1) THEN
       !$omp parallel do private(i)
       DO i=bo(1,1),bo(2,1)
          csIVals(i)=SIN(coeff*REAL(i,dp))
       END DO
    ELSE
       !$omp parallel do private(i)
       DO i=bo(1,1),bo(2,1)
          csIVals(i)=COS(coeff*REAL(i,dp))
       END DO
    END IF
    coeff=twopi/n_tot(2)
    IF (idir==2) THEN
       !$omp parallel do private(j)
       DO j=bo(1,2),bo(2,2)
          csJVals(j)=SIN(coeff*REAL(j,dp))
       END DO
    ELSE
       !$omp parallel do private(j)
       DO j=bo(1,2),bo(2,2)
          csJVals(j)=COS(coeff*REAL(j,dp))
       END DO
    END IF
    coeff=twopi/n_tot(3)
    IF (idir==3) THEN
       !$omp parallel do private(k)
       DO k=bo(1,3),bo(2,3)
          csKVals(k)=SIN(coeff*REAL(k,dp))
       END DO
    ELSE
       !$omp parallel do private(k)
       DO k=bo(1,3),bo(2,3)
          csKVals(k)=COS(coeff*REAL(k,dp))
       END DO
    END IF

    SELECT CASE(idir)
    CASE (1)
       ! x deriv
       !$omp parallel do private(ii,k,j,i,coeff,tmp)
       DO ii=1,SIZE(spline_g%cc)
          i=spline_g%pw_grid%g_hat(1,ii)
          j=spline_g%pw_grid%g_hat(2,ii)
          k=spline_g%pw_grid%g_hat(3,ii)
!FM                ! formula
!FM                coeff=(sinVal(1)*cosVal(2)*cosVal(3))/16.0_dp+&
!FM                     (sinVal(1)*cosVal(2)+sinVal(1)*cosVal(3))*3.0_dp/16.0_dp+&
!FM                     sinVal(1)*9.0_dp/16.0_dp
          tmp=csIVals(i)*csJVals(j)
          coeff=(tmp*csKVals(k)+&
               (tmp+csIVals(i)*csKVals(k))*3.0_dp+&
               csIVals(i)*9.0_dp)*inv16
          
          spline_g%cc(ii)=spline_g%cc(ii)*&
               CMPLX(0.0_dp,coeff,dp)
       END DO
    CASE(2)
       ! y deriv
       !$omp parallel do private(ii,k,j,i,coeff,tmp)
       DO ii=1,SIZE(spline_g%cc)
          i=spline_g%pw_grid%g_hat(1,ii)
          j=spline_g%pw_grid%g_hat(2,ii)
          k=spline_g%pw_grid%g_hat(3,ii)

          tmp=csIVals(i)*csJVals(j)
          coeff=(tmp*csKVals(k)+&
               (tmp+csJVals(j)*csKVals(k))*3.0_dp+&
               csJVals(j)*9.0_dp)*inv16
          
          spline_g%cc(ii)=spline_g%cc(ii)*&
               CMPLX(0.0_dp,coeff,dp)
       END DO
    CASE(3)
       ! z deriv
       !$omp parallel do private(ii,k,j,i,coeff,tmp)
       DO ii=1,SIZE(spline_g%cc)
          i=spline_g%pw_grid%g_hat(1,ii)
          j=spline_g%pw_grid%g_hat(2,ii)
          k=spline_g%pw_grid%g_hat(3,ii)

          tmp=csIVals(i)*csKVals(k)
          coeff=(tmp*csJVals(j)+&
               (tmp+csJVals(j)*csKVals(k))*3.0_dp+&
               csKVals(k)*9.0_dp)*inv16
          
          spline_g%cc(ii)=spline_g%cc(ii)*&
               CMPLX(0.0_dp,coeff,dp)
       END DO
    END SELECT

    DEALLOCATE(csIVals, csJVals, csKVals, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE pw_spline2_deriv_g
!***************************************************************************

!!****f* xc/pw_spline3_test_coeff [1.0] *
!!
!!   NAME
!!     pw_spline3_test_coeff
!!
!!   SYNOPSIS
!!     Subroutine pw_spline3_test_coeff(spline_coeff_r, spline_values_r,&
!!         error)
!!       Type(pw_type), Pointer:: spline_coeff_r, spline_values_r
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_spline3_test_coeff
!!
!!   FUNCTION
!!     checks if the given coefficents of the cubic spline really 
!!     interpolate the given values (looks at the differences on the grid)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - spline_coeff_r: the cofficents of the spline in realspace
!!     - spline_values_r: the values that the spline should interpolate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     05.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_spline3_test_coeff(spline_coeff_r, spline_values_r,error)
    TYPE(pw_type), POINTER :: spline_coeff_r, spline_values_r
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline3_test_coeff',&
         routineP=moduleN//':'//routineN
    REAL (KIND = dp), DIMENSION(3,3,3) :: spline_coeff
    REAL(KIND = dp) :: maxdiff,diff, c_val
    INTEGER :: i,j,k,ii,ij,ik
    INTEGER, DIMENSION(2,3) :: bo
    INTEGER, DIMENSION(3) :: n_tot

    failure=.FALSE.
    bo = spline_values_r%pw_grid%bounds_local
    n_tot = spline_values_r%pw_grid%npts

    spline_coeff=RESHAPE( (/&
         1._dp/(27._dp*8._dp),1._dp/(27._dp*2._dp),1._dp/(27._dp*8._dp),&
         1._dp/(27._dp*2._dp),2._dp/(27._dp),1._dp/(27._dp*2._dp),&
         1._dp/(27._dp*8._dp),1._dp/(27._dp*2._dp),1._dp/(27._dp*8._dp),&
         1._dp/(27._dp*2._dp),2._dp/(27._dp),1._dp/(27._dp*2._dp),&
         2._dp/(27._dp),8._dp/(27._dp),2._dp/(27._dp),&
         1._dp/(27._dp*2._dp),2._dp/(27._dp),1._dp/(27._dp*2._dp),&
         1._dp/(27._dp*8._dp),1._dp/(27._dp*2._dp),1._dp/(27._dp*8._dp),&
         1._dp/(27._dp*2._dp),2._dp/(27._dp),1._dp/(27._dp*2._dp),&
         1._dp/(27._dp*8._dp),1._dp/(27._dp*2._dp),1._dp/(27._dp*8._dp)&
         /),(/3,3,3/))

    CPPrecondition(ASSOCIATED(spline_coeff_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_coeff_r%in_space==REALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_coeff_r%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(spline_values_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_values_r%in_space==REALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_values_r%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_coeff_r%pw_grid%para%group_size==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       maxdiff=0.0_dp
       DO k = 0,spline_values_r%pw_grid%npts_local(3)-1
          DO j = 0,spline_values_r%pw_grid%npts_local(2)-1
             DO i = 0,spline_values_r%pw_grid%npts_local(1)-1

                c_val=0.0_dp
                DO ii=-1,1
                   DO ij=-1,1
                      DO ik=-1,1
                         c_val=c_val+spline_coeff(ii+2,ij+2,ik+2)*&
                              (&
                              spline_coeff_r%cr3d &
                              (bo(1,1)+MODULO(i+ii,n_tot(1)),&
                              bo(1,2)+MODULO(j+ij,n_tot(2)),&
                              bo(1,3)+MODULO(k+ik,n_tot(3)))&
                              )
                      END DO
                   END DO
                END DO
                diff=ABS(spline_values_r%cr3d &
                     (bo(1,1)+i,bo(1,2)+j,bo(1,3)+k)-c_val)
                IF (diff > maxdiff) THEN
                   maxdiff=diff
                   PRINT *,"maxdiff",maxdiff,i,j,k
                END IF

             END DO
          END DO
       END DO
    END IF
  END SUBROUTINE pw_spline3_test_coeff
!***************************************************************************

!!****f* xc/pw_spline2_test_coeff [1.0] *
!!
!!   NAME
!!     pw_spline_test_coeff
!!
!!   SYNOPSIS
!!     Subroutine pw_spline_test_coeff(spline_coeff_r, spline_values_r,&
!!         error)
!!       Type(pw_type), Pointer:: spline_coeff_r, spline_values_r
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_spline_test_coeff
!!
!!   FUNCTION
!!     checks if the coefficents of the quadratic spline really interpolate
!!     the given values (looks at the differences on the grid)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - spline_coeff_r: the cofficents of the spline in realspace
!!     - spline_values_r: the values that the spline should interpolate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     05.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_spline2_test_coeff(spline_coeff_r, spline_values_r,error)
    TYPE(pw_type), POINTER :: spline_coeff_r, spline_values_r
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline2_test_coeff',&
         routineP=moduleN//':'//routineN
    REAL (KIND = dp), DIMENSION(3,3,3) :: spline_coeff
    REAL(KIND = dp) :: maxdiff,diff, c_val
    INTEGER :: i,j,k,ii,ij,ik
    INTEGER, DIMENSION(2,3) :: bo
    INTEGER, DIMENSION(3) :: n_tot

    failure=.FALSE.
    bo = spline_values_r%pw_grid%bounds_local
    n_tot = spline_values_r%pw_grid%npts

    spline_coeff=RESHAPE( (/&
         1._dp/(64._dp*8._dp),3._dp/(64._dp*4._dp),1._dp/(64._dp*8._dp),&
         3._dp/(64._dp*4._dp),9._dp/(64._dp*2._dp),3._dp/(64._dp*4._dp),&
         1._dp/(64._dp*8._dp),3._dp/(64._dp*4._dp),1._dp/(64._dp*8._dp),&
         3._dp/(64._dp*4._dp),9._dp/(64._dp*2._dp),3._dp/(64._dp*4._dp),&
         9._dp/(64._dp*2._dp),27._dp/(64._dp),9._dp/(64._dp*2_dp),&
         3._dp/(64._dp*4._dp),9._dp/(64._dp*2._dp),3._dp/(64._dp*4._dp),&
         1._dp/(64._dp*8._dp),3._dp/(64._dp*4._dp),1._dp/(64._dp*8._dp),&
         3._dp/(64._dp*4._dp),9._dp/(64._dp*2._dp),3._dp/(64._dp*4._dp),&
         1._dp/(64._dp*8._dp),3._dp/(64._dp*4._dp),1._dp/(64._dp*8._dp)&
         /),(/3,3,3/))

    CPPrecondition(ASSOCIATED(spline_coeff_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_coeff_r%in_space==REALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_coeff_r%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(spline_values_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_values_r%in_space==REALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_values_r%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_coeff_r%pw_grid%para%group_size==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       maxdiff=0.0_dp
       DO k = 0,spline_values_r%pw_grid%npts_local(3)-1
          DO j = 0,spline_values_r%pw_grid%npts_local(2)-1
             DO i = 0,spline_values_r%pw_grid%npts_local(1)-1

                c_val=0.0_dp
                DO ii=-1,1
                   DO ij=-1,1
                      DO ik=-1,1
                         c_val=c_val+spline_coeff(ii+2,ij+2,ik+2)*&
                              (&
                              spline_coeff_r%cr3d &
                              (bo(1,1)+MODULO(i+ii,n_tot(1)),&
                              bo(1,2)+MODULO(j+ij,n_tot(2)),&
                              bo(1,3)+MODULO(k+ik,n_tot(3)))&
                              )
                      END DO
                   END DO
                END DO
                diff=ABS(spline_values_r%cr3d &
                     (bo(1,1)+i,bo(1,2)+j,bo(1,3)+k)-c_val)
                IF (diff > maxdiff) THEN
                   maxdiff=diff
                   PRINT *,"maxdiff",maxdiff,i,j,k
                END IF

             END DO
          END DO
       END DO
    END IF
  END SUBROUTINE pw_spline2_test_coeff
!***************************************************************************

END MODULE pw_spline_utils
