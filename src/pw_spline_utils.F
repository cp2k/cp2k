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
  USE coefficient_types,               ONLY: SQUARE
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
  USE message_passing,                 ONLY: mp_sendrecv
  USE pw_grid_types,                   ONLY: FULLSPACE
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type,&
                                             pw_prolongate,&
                                             pw_inject,&
                                             pw_zero,&
                                             pw_transfer,&
                                             pw_sumup,&
                                             pw_copy
  USE pw_pool_types, ONLY: pw_pool_type, pw_pool_create_pw,&
       pw_pool_give_back_pw
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE realspace_grid_types,            ONLY: rs_pw_to_cube
  USE pws, ONLY: pw_integral_ab, pw_integral_aa

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pw_spline_utils'

    REAL ( KIND = dp ), PUBLIC, PARAMETER, DIMENSION(4) :: nn10_coeffs=&
         (/ 125._dp/216._dp, 25._dp/432._dp, 5._dp/864._dp, 1._dp/1728._dp /),&
         spline3_coeffs=&
         (/ 8._dp/(27._dp),2._dp/(27._dp),1._dp/(27._dp*2._dp),&
         1._dp/(27._dp*8._dp) /),&
         spline2_coeffs=&
         (/ 27._dp/(64._dp),9._dp/(64._dp*2_dp),3._dp/(64._dp*4._dp),&
         1._dp/(64._dp*8._dp) /),&
         nn50_coeffs=&
         (/ 15625._dp/17576._dp, 625._dp/35152._dp, 25._dp/70304._dp,&
         1._dp/140608._dp /),&
         spl3_1d_coeffs=&
         (/ 2._dp/3._dp,23._dp/48._dp,1._dp/6._dp,1._dp/48._dp /),&
         spl3_aint_coeff=&
         (/ 46._dp/27._dp,-2._dp/(27._dp),-1._dp/(27._dp*2._dp),&
         -1._dp/(27._dp*8._dp) /),&
         spl3_precond1_coeff=&
         (/ 64._dp/3._dp,-8._dp/3._dp,-1._dp/3._dp,&
         -1._dp/24._dp /)
    
         
    REAL ( KIND = dp ), PUBLIC, PARAMETER, DIMENSION(3) :: spline3_deriv_coeffs=&
         (/ 2.0_dp/9.0_dp,  1.0_dp/18.0_dp, 1.0_dp/72.0_dp /),&
         spline2_deriv_coeffs=&
         (/ 9.0_dp/32.0_dp, 3.0_dp/64.0_dp, 1.0_dp/128.0_dp /),&
         nn10_deriv_coeffs=&
         (/ 25._dp/72._dp, 5._dp/144, 1._dp/288._dp /),&
         nn50_deriv_coeffs=&
         (/ 625._dp/1352._dp, 25._dp/2704._dp,1._dp/5408._dp /)

  PUBLIC :: pw_spline3_interpolate_values_g, pw_spline3_evaluate_values_g,&
       pw_spline3_deriv_g
  PUBLIC :: pw_spline_scale_deriv, pw_nn_smear_g
  PUBLIC :: pw_spline2_interpolate_values_g, pw_spline2_evaluate_values_g,&
       pw_spline2_deriv_g
  PUBLIC :: pw_nn_compose_r, pw_nn_smear_r, pw_nn_deriv_r, pw_inject_s3

  INTEGER, PRIVATE, PARAMETER :: method_rs=2
!***
!****************************************************************************

  !!****s* pw_spline_utils/pw_spline_precond_type [1.0] *
  !!
  !!   NAME
  !!     pw_spline_precond_type
  !!
  !!   FUNCTION
  !!     stores information for the preconditioner used to calculate the 
  !!     coeffs of splines
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   SOURCE
  !***************************************************************************
  type pw_spline_precond_type
     INTEGER :: ref_count,id_nr
  end type pw_spline_precond_type
  !!***
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
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
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
    TYPE(pw_type), POINTER                   :: spline_g
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'pw_spline2_interpolate_values_g', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: gbo
    INTEGER, DIMENSION(3)                    :: n_tot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: c23, coeff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cosIVals, cosJVals, cosKVals

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    n_tot(1:3) = spline_g%pw_grid%npts (1:3)
    gbo = spline_g%pw_grid%bounds

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
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine pw_spline3_interpolate_values_g
!!
!!   FUNCTION
!!     calculates the FFT of the coefficents of the2 cubic spline that
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
    TYPE(pw_type), POINTER                   :: spline_g
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'pw_spline3_interpolate_values_g', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: gbo
    INTEGER, DIMENSION(3)                    :: n_tot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: c23, coeff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cosIVals, cosJVals, cosKVals

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    n_tot(1:3) = spline_g%pw_grid%npts (1:3)
    gbo = spline_g%pw_grid%bounds

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
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
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
    TYPE(pw_type), POINTER                   :: spline_g
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_spline2_evaluate_values_g', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    INTEGER, DIMENSION(3)                    :: n, n_tot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: c23, coeff, inv64
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cosIVals, cosJVals, cosKVals

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
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
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
    TYPE(pw_type), POINTER                   :: spline_g
    REAL(KIND=dp), DIMENSION(4), INTENT(in)  :: coeffs
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_nn_smear_g', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    INTEGER, DIMENSION(3)                    :: n, n_tot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: c23, coeff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cosIVals, cosJVals, cosKVals
    REAL(KIND=dp), DIMENSION(4)              :: r_coeffs

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
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
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
    TYPE(pw_type), POINTER                   :: spline_g
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_spline3_evaluate_values_g', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    INTEGER, DIMENSION(3)                    :: n, n_tot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: c23, coeff, inv27
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cosIVals, cosJVals, cosKVals

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
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
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
    TYPE(pw_p_type), DIMENSION(3)            :: deriv_vals_r
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL, INTENT(in), OPTIONAL            :: transpose
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: scale
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_spline_scale_deriv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, idir, j, k
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(3)                    :: n_tot
    LOGICAL                                  :: diag, failure, my_transpose
    REAL(KIND=dp)                            :: dVal1, dVal2, dVal3, &
                                                my_scale, scalef
    REAL(KIND=dp), DIMENSION(3, 3)           :: h_grid
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: ddata, ddata2, ddata3

! sun does not like that its components are modified

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
    diag=.TRUE.
    IF (my_transpose) THEN
       DO j=1,3
          DO i=1,3
             h_grid(j,i)=my_scale*REAL(n_tot(i),dp)*cell%h_inv(i,j)
             IF (i/=j.AND.h_grid(j,i)/=0.0_dp) diag=.FALSE.
          END DO
       END DO
    ELSE
       DO j=1,3
          DO i=1,3
             h_grid(i,j)=my_scale*REAL(n_tot(i),dp)*cell%h_inv(i,j)      
             IF (i/=j.AND.h_grid(i,j)/=0.0_dp) diag=.FALSE.
          END DO
       END DO
    END IF

    IF (diag) THEN
       DO idir=1,3
          ddata => deriv_vals_r(idir)%pw%cr3d
          scalef=h_grid(idir,idir)
          CALL dscal((bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1),&
               scalef,ddata,1)
!FM          !$omp parallel do default(none) private(k,j,i) shared(ddata,scalef)
!FM          DO k = bo(1,3), bo(2,3)
!FM             DO j = bo(1,2), bo(2,2)
!FM                DO i = bo(1,1), bo(2,1)
!FM                   ddata(i,j,k)=scalef*ddata(i,j,k)
!FM                END DO
!FM             END DO
!FM          END DO
       END DO
    ELSE
       ddata => deriv_vals_r(1)%pw%cr3d
       ddata2 => deriv_vals_r(2)%pw%cr3d
       ddata3 => deriv_vals_r(3)%pw%cr3d
       !$omp parallel do default(none) private(k,j,i,dVal1,dVal2,dVal3) &
       !$omp          shared(ddata,ddata2,ddata3,h_grid,bo)
       DO k = bo(1,3), bo(2,3)
          DO j = bo(1,2), bo(2,2)
             DO i = bo(1,1), bo(2,1)

                dVal1=ddata(i,j,k)
                dVal2=ddata2(i,j,k)
                dVal3=ddata3(i,j,k)

                ddata(i,j,k)=h_grid(1,1)*dVal1+&
                     h_grid(2,1)*dVal2+ h_grid(3,1)*dVal3
                ddata2(i,j,k)=h_grid(1,2)*dVal1+&
                     h_grid(2,2)*dVal2+ h_grid(3,2)*dVal3
                ddata3(i,j,k)=h_grid(1,3)*dVal1+&
                     h_grid(2,3)*dVal2+ h_grid(3,3)*dVal3

             END DO
          END DO
       END DO
    END IF

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
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
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
    TYPE(pw_type), POINTER                   :: spline_g
    INTEGER, INTENT(in)                      :: idir
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_spline3_deriv_g', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    INTEGER, DIMENSION(3)                    :: n, n_tot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: coeff, inv9, tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: csIVals, csJVals, csKVals

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
    TYPE(pw_type), POINTER                   :: spline_g
    INTEGER, INTENT(in)                      :: idir
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_spline2_deriv_g', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(3)                    :: n, n_tot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: coeff, inv16, tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: csIVals, csJVals, csKVals

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

!!****f* pw_spine_utils/pw_compose_stripe [1.0] *
!!
!!   NAME
!!     pw_compose_stripe
!!
!!   FUNCTION
!!     applies a nearest neighbor linear operator to a stripe in x direction:
!!     out_val(i)=sum(weight(j)*in_val(i+j-1),j=0..2)
!!
!!   NOTES
!!     uses 2 read streams and 1 write stream
!!
!!   ARGUMENTS
!!     - weights: the weights of the linear operator
!!     - in_val: the argument to the operator
!!     - in_val_first: the first argument (needed to calculate out_val(1))
!!     - in_val_last: the last argument (needed to calculate out_val(n_el))
!!     - out_val: the place where the result is accumulated
!!     - n_el: the number of elements in in_v and out_v
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2004 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_compose_stripe(weights,in_val,in_val_first,in_val_last,&
     out_val,n_el)
    REAL(kind=dp), DIMENSION(0:2), &
      INTENT(in)                             :: weights
    INTEGER                                  :: n_el
    REAL(kind=dp), DIMENSION(1:n_el), &
      INTENT(inout)                          :: out_val
    REAL(kind=dp), INTENT(in)                :: in_val_last, in_val_first
    REAL(kind=dp), DIMENSION(1:n_el), &
      INTENT(in)                             :: in_val

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_compose_stripe', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: v0, v1, v2

!FM    TYPE(cp_error_type), INTENT(inout), &
!FM      OPTIONAL                               :: error
!    REAL(kind=dp), DIMENSION(0:2)            :: v

  failure=.FALSE.
  v0=in_val_first
  v1=in_val(1)
  IF (weights(1)==0.0_dp) THEN
     ! optimized version for x deriv
     DO i=1,n_el-3,3
        v2=in_val(i+1)
        out_val(i)= out_val(i)+&
             weights(0)*v0+&
             weights(2)*v2
        v0=in_val(i+2)
        out_val(i+1)= out_val(i+1)+&
             weights(0)*v1+&
             weights(2)*v0
        v1=in_val(i+3)
        out_val(i+2)= out_val(i+2)+&
             weights(0)*v2+&
             weights(2)*v1
     END DO
  ELSE
     ! generic version
     DO i=1,n_el-3,3
        v2=in_val(i+1)
        out_val(i)= out_val(i)+&
             weights(0)*v0+&
             weights(1)*v1+&
             weights(2)*v2
        v0=in_val(i+2)
        out_val(i+1)= out_val(i+1)+&
             weights(0)*v1+&
             weights(1)*v2+&
             weights(2)*v0
        v1=in_val(i+3)
        out_val(i+2)= out_val(i+2)+&
             weights(0)*v2+&
             weights(1)*v0+&
             weights(2)*v1
     END DO
  END IF
  SELECT CASE(MODULO(n_el-1,3))
  CASE(0)
     IF (n_el>0) THEN
        v2=in_val_last
        out_val(n_el)= out_val(n_el)+&
             weights(0)*v0+&
             weights(1)*v1+&
             weights(2)*v2
     END IF
  CASE(1)
     v2=in_val(n_el)
     out_val(n_el-1)= out_val(n_el-1)+&
          weights(0)*v0+&
          weights(1)*v1+&
          weights(2)*v2
     v0=in_val_last
     out_val(n_el)= out_val(n_el)+&
          weights(0)*v1+&
          weights(1)*v2+&
          weights(2)*v0
  CASE(2)
     v2=in_val(n_el-1)
     out_val(n_el-2)= out_val(n_el-2)+&
          weights(0)*v0+&
          weights(1)*v1+&
          weights(2)*v2
     v0=in_val(n_el)
     out_val(n_el-1)= out_val(n_el-1)+&
          weights(0)*v1+&
          weights(1)*v2+&
          weights(2)*v0
     v1=in_val_last
     out_val(n_el)= out_val(n_el)+&
          weights(0)*v2+&
          weights(1)*v0+&
          weights(2)*v1
  END SELECT
  
END SUBROUTINE pw_compose_stripe
!***************************************************************************

!!****f* pw_spine_utils/pw_compose_stripe_r3 [1.0] *
!!
!!   NAME
!!     pw_compose_stripe_r3
!!
!!   FUNCTION
!!     applies a linear operator using the nearest neighbors in x and y
!!     direction to a stripe in x direction:
!!     out_val(i)=sum(weight(j,0)*in_val_0(i+j-1)+
!!                    weight(j,1)*in_val_1(i+j-1)+
!!                    weight(j,2)*in_val_2(i+j-1),j=0..2)
!!
!!   NOTES
!!     if method_rs>=4 uses 4 read streams and one write stream
!!
!!   ARGUMENTS
!!     - weights: the weights of the linear operator
!!     - in_val_0: the first stripe in x direction of the argument
!!       of the operator
!!     - in_val_1: the central stripe in x direction of the argument
!!       of the operator
!!     - in_val_2: the last stripe in x direction of the argument
!!       of the operator
!!     - in_val_first: the first arguments (needed to calculate out_val(1))
!!     - in_val_last: the last arguments (needed to calculate out_val(n_el))
!!     - out_val: the place where the result is accumulated
!!     - n_el: the number of elements in in_v and out_v
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2004 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_compose_stripe_r3(weights,in_val_0,in_val_1,in_val_2,&
     in_val_first,in_val_last,out_val,n_el)
    REAL(kind=dp), DIMENSION(0:2, 0:2), &
      INTENT(in)                             :: weights
    INTEGER                                  :: n_el
    REAL(kind=dp), DIMENSION(1:n_el), &
      INTENT(inout)                          :: out_val
    REAL(kind=dp), DIMENSION(0:2), &
      INTENT(in)                             :: in_val_last, in_val_first
    REAL(kind=dp), DIMENSION(1:n_el), &
      INTENT(in)                             :: in_val_2, in_val_1, in_val_0

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_compose_stripe_r3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: xderiv, yderiv
    REAL(kind=dp)                            :: v_0_0, v_0_1, v_0_2, v_1_0, &
                                                v_1_1, v_1_2, v_2_0, v_2_1, &
                                                v_2_2

!FM    TYPE(cp_error_type), INTENT(inout), &
!FM      OPTIONAL                               :: error
!REAL(kind=dp), DIMENSION(0:2, 0:2)       :: v

  yderiv=ALL(weights(:,1)==0.0_dp)
  xderiv=ALL(weights(1,:)==0.0_dp)
!FM  IF (method_rs<3 .OR.((.NOT.yderiv).AND.method_rs<4)) THEN
!FM     CALL pw_compose_stripe(weights=weights(:,0),in_val=in_val_0,&
!FM          in_val_first=in_val_first(0),in_val_last=in_val_last(0),out_val=out_val,&
!FM          n_el=n_el)
!FM     IF (.NOT.yderiv) THEN
!FM        CALL pw_compose_stripe(weights=weights(:,1),in_val=in_val_1,&
!FM             in_val_first=in_val_first(1),in_val_last=in_val_last(1),out_val=out_val,&
!FM             n_el=n_el)
!FM     END IF
!FM     CALL pw_compose_stripe(weights=weights(:,2),in_val=in_val_2,&
!FM          in_val_first=in_val_first(2),in_val_last=in_val_last(2),out_val=out_val,&
!FM          n_el=n_el)
!FM  ELSE
     v_0_0=in_val_first(0)
     v_0_1=in_val_first(1)
     v_0_2=in_val_first(2)
     v_1_0=in_val_0(1)
     v_1_1=in_val_1(1)
     v_1_2=in_val_2(1)

     IF (xderiv) THEN
        ! optimized version
        DO i=1,n_el-3,3
           v_2_0=in_val_0(i+1)
           v_2_1=in_val_1(i+1)
           v_2_2=in_val_2(i+1)
           out_val(i)= out_val(i)+&
                weights(0,0)*v_0_0+&
                weights(2,0)*v_2_0+&
                weights(0,1)*v_0_1+&
                weights(2,1)*v_2_1+&
                weights(0,2)*v_0_2+&
                weights(2,2)*v_2_2
           v_0_0=in_val_0(i+2)
           v_0_1=in_val_1(i+2)
           v_0_2=in_val_2(i+2)
           out_val(i+1)= out_val(i+1)+&
                   weights(0,0)*v_1_0+&
                   weights(2,0)*v_0_0+&
                   weights(0,1)*v_1_1+&
                   weights(2,1)*v_0_1+&
                   weights(0,2)*v_1_2+&
                   weights(2,2)*v_0_2
           v_1_0=in_val_0(i+3)
           v_1_1=in_val_1(i+3)
           v_1_2=in_val_2(i+3)
           out_val(i+2)= out_val(i+2)+&
                weights(0,0)*v_2_0+&
                weights(2,0)*v_1_0+&
                weights(0,1)*v_2_1+&
                weights(2,1)*v_1_1+&
                weights(0,2)*v_2_2+&
                weights(2,2)*v_1_2
        END DO
     ELSEIF (yderiv) THEN
        ! optimized version
        DO i=1,n_el-3,3
           v_2_0=in_val_0(i+1)
           v_2_2=in_val_2(i+1)
           out_val(i)= out_val(i)+&
                weights(0,0)*v_0_0+&
                weights(1,0)*v_1_0+&
                weights(2,0)*v_2_0+&
                weights(0,2)*v_0_2+&
                weights(1,2)*v_1_2+&
                weights(2,2)*v_2_2
           v_0_0=in_val_0(i+2)
           v_0_2=in_val_2(i+2)
           out_val(i+1)= out_val(i+1)+&
                weights(0,0)*v_1_0+&
                weights(1,0)*v_2_0+&
                weights(2,0)*v_0_0+&
                weights(0,2)*v_1_2+&
                weights(1,2)*v_2_2+&
                weights(2,2)*v_0_2
           v_1_0=in_val_0(i+3)
           v_1_2=in_val_2(i+3)
           out_val(i+2)= out_val(i+2)+&
                weights(0,0)*v_2_0+&
                weights(1,0)*v_0_0+&
                weights(2,0)*v_1_0+&
                weights(0,2)*v_2_2+&
                weights(1,2)*v_0_2+&
                weights(2,2)*v_1_2
        END DO
     ELSE
        ! generic version
        DO i=1,n_el-3,3
           v_2_0=in_val_0(i+1)
           v_2_1=in_val_1(i+1)
           v_2_2=in_val_2(i+1)
           out_val(i)= out_val(i)+&
                weights(0,0)*v_0_0+&
                weights(1,0)*v_1_0+&
                weights(2,0)*v_2_0+&
                weights(0,1)*v_0_1+&
                weights(1,1)*v_1_1+&
                weights(2,1)*v_2_1+&
                weights(0,2)*v_0_2+&
                weights(1,2)*v_1_2+&
                weights(2,2)*v_2_2
           v_0_0=in_val_0(i+2)
           v_0_1=in_val_1(i+2)
           v_0_2=in_val_2(i+2)
           out_val(i+1)= out_val(i+1)+&
                weights(0,0)*v_1_0+&
                weights(1,0)*v_2_0+&
                weights(2,0)*v_0_0+&
                weights(0,1)*v_1_1+&
                weights(1,1)*v_2_1+&
                weights(2,1)*v_0_1+&
                weights(0,2)*v_1_2+&
                weights(1,2)*v_2_2+&
                weights(2,2)*v_0_2
           v_1_0=in_val_0(i+3)
           v_1_1=in_val_1(i+3)
           v_1_2=in_val_2(i+3)
           out_val(i+2)= out_val(i+2)+&
                weights(0,0)*v_2_0+&
                weights(1,0)*v_0_0+&
                weights(2,0)*v_1_0+&
                weights(0,1)*v_2_1+&
                weights(1,1)*v_0_1+&
                weights(2,1)*v_1_1+&
                weights(0,2)*v_2_2+&
                weights(1,2)*v_0_2+&
                weights(2,2)*v_1_2
        END DO
     END IF
     SELECT CASE(MODULO(n_el-1,3))
     CASE(0)
        IF (n_el>0) THEN
           v_2_0=in_val_last(0)
           v_2_1=in_val_last(1)
           v_2_2=in_val_last(2)
           out_val(n_el)= out_val(n_el)+&
                weights(0,0)*v_0_0+&
                weights(1,0)*v_1_0+&
                weights(2,0)*v_2_0+&
                weights(0,1)*v_0_1+&
                weights(1,1)*v_1_1+&
                weights(2,1)*v_2_1+&
                weights(0,2)*v_0_2+&
                weights(1,2)*v_1_2+&
                weights(2,2)*v_2_2
        END IF
     CASE(1)
        v_2_0=in_val_0(n_el)
        v_2_1=in_val_1(n_el)
        v_2_2=in_val_2(n_el)
        out_val(n_el-1)= out_val(n_el-1)+&
             weights(0,0)*v_0_0+&
             weights(1,0)*v_1_0+&
             weights(2,0)*v_2_0+&
             weights(0,1)*v_0_1+&
             weights(1,1)*v_1_1+&
             weights(2,1)*v_2_1+&
             weights(0,2)*v_0_2+&
             weights(1,2)*v_1_2+&
             weights(2,2)*v_2_2
        v_0_0=in_val_last(0)
        v_0_1=in_val_last(1)
        v_0_2=in_val_last(2)
        out_val(n_el)= out_val(n_el)+&
             weights(0,0)*v_1_0+&
             weights(1,0)*v_2_0+&
             weights(2,0)*v_0_0+&
             weights(0,1)*v_1_1+&
             weights(1,1)*v_2_1+&
             weights(2,1)*v_0_1+&
             weights(0,2)*v_1_2+&
             weights(1,2)*v_2_2+&
             weights(2,2)*v_0_2
     CASE(2)
        v_2_0=in_val_0(n_el-1)
        v_2_1=in_val_1(n_el-1)
        v_2_2=in_val_2(n_el-1)
        out_val(n_el-2)= out_val(n_el-2)+&
             weights(0,0)*v_0_0+&
             weights(1,0)*v_1_0+&
             weights(2,0)*v_2_0+&
             weights(0,1)*v_0_1+&
             weights(1,1)*v_1_1+&
             weights(2,1)*v_2_1+&
             weights(0,2)*v_0_2+&
             weights(1,2)*v_1_2+&
             weights(2,2)*v_2_2
        v_0_0=in_val_0(n_el)
        v_0_1=in_val_1(n_el)
        v_0_2=in_val_2(n_el)
        out_val(n_el-1)= out_val(n_el-1)+&
             weights(0,0)*v_1_0+&
             weights(1,0)*v_2_0+&
             weights(2,0)*v_0_0+&
             weights(0,1)*v_1_1+&
             weights(1,1)*v_2_1+&
             weights(2,1)*v_0_1+&
             weights(0,2)*v_1_2+&
             weights(1,2)*v_2_2+&
             weights(2,2)*v_0_2
        v_1_0=in_val_last(0)
        v_1_1=in_val_last(1)
        v_1_2=in_val_last(2)
        out_val(n_el)= out_val(n_el)+&
             weights(0,0)*v_2_0+&
             weights(1,0)*v_0_0+&
             weights(2,0)*v_1_0+&
             weights(0,1)*v_2_1+&
             weights(1,1)*v_0_1+&
             weights(2,1)*v_1_1+&
             weights(0,2)*v_2_2+&
             weights(1,2)*v_0_2+&
             weights(2,2)*v_1_2
     END SELECT
!FM  END IF
  
END SUBROUTINE pw_compose_stripe_r3
!***************************************************************************

!!****f* pw_spline_utils/pw_nn_compose_r_work [1.0] *
!!
!!   NAME
!!     pw_nn_compose_r_work
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     private routine that computes pw_nn_compose_r (it seems that without 
!!     passing arrays in this way either some compiler do a copyin/out (xlf)
!!     or by inlining suboptimal code is produced (nag)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - weights: a 3x3x3 array with the linear operator
!!     - in_val: the argument for the linear operator
!!     - out_val: place where the value of the linear oprator should be added
!!     - pw_in: pw to be able to get the needed meta data about in_val and 
!!       out_val
!!     - bo: boundaries of in_val and out_val
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pw_nn_compose_r_work(weights,in_val,out_val,pw_in,bo,pbc,error)
    REAL(kind=dp), DIMENSION(0:2, 0:2, 0:2)  :: weights
    INTEGER, DIMENSION(2, 3)                 :: bo
    TYPE(pw_type), POINTER                   :: pw_in
    REAL(kind=dp), DIMENSION(bo(1, 1):bo(2, &
      1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, &
      3)), INTENT(inout)                     :: out_val
    REAL(kind=dp), DIMENSION(bo(1, 1):bo(2, &
      1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, &
      3)), INTENT(in)                        :: in_val
    LOGICAL, OPTIONAL, intent(in)            :: pbc
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "pw_nn_compose_r_work", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: i, j, jw, k, kw, myj, myk, &
                                                stat
    INTEGER, DIMENSION(2, 3)                 :: gbo
    INTEGER, DIMENSION(3)                    :: s
    LOGICAL                                  :: failure, has_boundary, &
                                                yderiv, zderiv,my_pbc
    REAL(kind=dp)                            :: in_val_f, in_val_l
    REAL(kind=dp), DIMENSION(0:2)            :: in_val_first, in_val_last
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: l_boundary, tmp, u_boundary

     zderiv=ALL(weights(:,:,1)==0.0_dp)
     yderiv=ALL(weights(:,1,:)==0.0_dp)
     bo=pw_in%pw_grid%bounds_local
     gbo=pw_in%pw_grid%bounds
     DO i=1,3
        s(i)=bo(2,i)-bo(1,i)+1
     END DO
     my_pbc=.true.
     if (present(pbc)) my_pbc=pbc
     IF (ANY(s<1)) RETURN
     has_boundary= ANY(pw_in%pw_grid%bounds_local(:,1) /= &
          pw_in%pw_grid%bounds(:,1)).or..not.my_pbc
     IF (has_boundary) THEN
        ALLOCATE(l_boundary(bo(1,2):bo(2,2),bo(1,3):bo(2,3)),&
             u_boundary(bo(1,2):bo(2,2),bo(1,3):bo(2,3)),&
             tmp(bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        IF (my_pbc) THEN ! to do
        tmp(:,:)=pw_in%cr3d(bo(2,1),:,:)
        CALL mp_sendrecv(tmp,pw_in%pw_grid%para%pos_of_x(&
             gbo(1,1)+MODULO(bo(2,1)+1-gbo(1,1),gbo(2,1)-gbo(1,1)+1)),&
             l_boundary,pw_in%pw_grid%para%pos_of_x(&
             gbo(1,1)+MODULO(bo(1,1)-1-gbo(1,1),gbo(2,1)-gbo(1,1)+1)),&
             pw_in%pw_grid%para%group)
        tmp(:,:)=pw_in%cr3d(bo(1,1),:,:)
        CALL mp_sendrecv(tmp,pw_in%pw_grid%para%pos_of_x(&
             gbo(1,1)+MODULO(bo(1,1)-1-gbo(1,1),gbo(2,1)-gbo(1,1)+1)),&
             u_boundary,pw_in%pw_grid%para%pos_of_x(&
             gbo(1,1)+MODULO(bo(2,1)+1-gbo(1,1),gbo(2,1)-gbo(1,1)+1)),&
             pw_in%pw_grid%para%group)
        END IF
        DEALLOCATE(tmp,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

        IF (.NOT.my_pbc) THEN
           IF (gbo(1,1)==bo(1,1)) THEN
              l_boundary=0.0_dp
           END IF
           IF (gbo(2,1)==bo(2,1)) THEN
              u_boundary=0.0_dp
           END IF
        END IF

        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF

     IF (method_rs<3) THEN
        !$omp parallel do default(none) private(k,kw,myk,j,jw,myj,in_val_f,&
        !$omp     in_val_l) shared(zderiv,yderiv,bo,in_val,out_val,s,l_boundary,&
        !$omp     u_boundary,weights,has_boundary)
        DO k=0,s(3)-1
           DO kw=0,2
              myk=bo(1,3)+MODULO(k+kw-1,s(3))
              IF (zderiv.AND.kw==1.or..NOT.my_pbc.AND.myk/=bo(1,3)+k+kw-1) CYCLE
              DO j=0,s(2)-1
                 DO jw=0,2
                    myj=bo(1,2)+MODULO(j+jw-1,s(2))
                    IF (yderiv.AND.jw==1.OR..not.my_pbc.AND.myj/=bo(1,2)+j+jw-1) CYCLE
                    IF (has_boundary) THEN
                       in_val_f=l_boundary(myj,myk)
                       in_val_l=u_boundary(myj,myk)
                    ELSE
                       in_val_f=in_val(bo(2,1),myj,myk)
                       in_val_l=in_val(bo(1,1),myj,myk)
                    END IF
                    CALL pw_compose_stripe(weights=weights(:,jw,kw),&
                         in_val=in_val(:,myj,myk),&
                         in_val_first=in_val_f,in_val_last=in_val_l,&
                         out_val=out_val(:,bo(1,2)+j,bo(1,3)+k),n_el=s(1)) 
                 END DO
              END DO
           END DO
        END DO
     ELSE
        CPAssert(my_pbc,cp_failure_level,routineP,error,failure)
        !$omp parallel do default(none) private(k,kw,myk,j,myj,in_val_first,&
        !$omp     in_val_last) shared(zderiv,yderiv,bo,in_val,out_val,s,&
        !$omp     l_boundary, u_boundary,weights,has_boundary)
        DO k=0,s(3)-1
           DO j=0,s(2)-1
              DO kw=0,2
                 myk=bo(1,3)+MODULO(k+kw-1,s(3))
                 IF (zderiv.AND.kw==1) CYCLE
                 IF (has_boundary) THEN
                    DO jw=0,2
                       myj=bo(1,2)+MODULO(j+jw-1,s(2))
                       in_val_first(jw)=l_boundary(myj, myk)
                       in_val_last(jw)=u_boundary(myj, myk)
                    END DO
                 ELSE
                    DO jw=0,2
                       myj=bo(1,2)+MODULO(j+jw-1,s(2))
                       in_val_first(jw)=in_val(bo(2,1),&
                            myj,myk)
                       in_val_last(jw)=in_val(bo(1,1),&
                            myj,myk)
                    END DO
                 END IF
                 CALL pw_compose_stripe_r3(weights=weights(:,:,kw),&
                      in_val_0=in_val(:,bo(1,2)+MODULO(j-1,s(2)), myk),&
                      in_val_1=in_val(:,bo(1,2)+j, myk),&
                      in_val_2=in_val(:,bo(1,2)+MODULO(j+1,s(2)), myk),&
                      in_val_first=in_val_first,&
                      in_val_last=in_val_last,&
                      out_val=out_val(:,bo(1,2)+j, bo(1,3)+k),&
                      n_el=s(1))
              END DO
           END DO
        END DO
     END IF
     IF (has_boundary) THEN
        DEALLOCATE(l_boundary,u_boundary,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
   END SUBROUTINE pw_nn_compose_r_work
!***************************************************************************

!!****f* pw_spline_utils/pw_nn_compose_r [1.0] *
!!
!!   NAME
!!     pw_nn_compose_r
!!
!!   SYNOPSIS
!!     Subroutine pw_nn_compose_r(weights, pw_in, pw_out, error)
!!       Real(Kind=dp), Dimension(0:2,0:2,0:2):: weights
!!       Type(pw_type), Pointer:: pw_in, pw_out
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine pw_nn_compose_r
!!
!!   FUNCTION
!!     applies a nearest neighbor linear operator to a pw in real space
!!
!!   NOTES
!!     has specialized versions for derivative operator (with central values==0)
!!
!!   ARGUMENTS
!!     - weights: a 3x3x3 array with the linear operator
!!     - pw_in: the argument for the linear operator
!!     - pw_out: place where the value of the linear oprator should be added
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pw_nn_compose_r(weights, pw_in, pw_out, pbc,error)
    REAL(kind=dp), DIMENSION(0:2, 0:2, 0:2)  :: weights
    TYPE(pw_type), POINTER                   :: pw_in, pw_out
    LOGICAL, OPTIONAL, intent(in)            :: pbc
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_nn_compose_r', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

  CALL timeset(routineN,"I","",handle)
  failure=.FALSE.
  CPPrecondition(ASSOCIATED(pw_in),cp_failure_level,routineP,error,failure)
  CPPrecondition(pw_in%in_space==REALSPACE,cp_failure_level,routineP,error,failure)
  CPPrecondition(pw_in%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(pw_out),cp_failure_level,routineP,error,failure)
  CPPrecondition(pw_out%in_space==REALSPACE,cp_failure_level,routineP,error,failure)
  CPPrecondition(pw_out%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
  CALL cp_assert(ALL(pw_in%pw_grid%bounds_local(:,2:3) == pw_in%pw_grid%bounds(:,2:3)),&
       cp_failure_level,cp_assertion_failed,routineP,&
       "wrong pw distribution",error,failure)
  IF (.NOT.failure) THEN
     CALL pw_nn_compose_r_work(weights=weights,in_val=pw_in%cr3d,&
          out_val=pw_out%cr3d,pw_in=pw_in,bo=pw_in%pw_grid%bounds_local,&
          pbc=pbc,error=error)
  END IF
  CALL timestop(0.0_dp,handle)
END SUBROUTINE pw_nn_compose_r
!***************************************************************************

!!****f* pw_spline_utils/pw_nn_smear_r [1.0] *
!!
!!   NAME
!!     pw_nn_smear_r
!!
!!   SYNOPSIS
!!     Subroutine pw_nn_smear_r(pw_in, pw_out, coeffs, error)
!!       Type(pw_type), Pointer:: pw_in, pw_out
!!       Real(Kind=dp), Dimension(4), Intent (IN):: coeffs
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine pw_nn_smear_r
!!
!!   FUNCTION
!!     calculates the values of a nearest neighbor smearing
!!
!!   NOTES
!!     does not normalize the smear to 1.
!!     with coeff=(/ 8._dp/27._dp, 2._dp/27._dp, 1._dp/54._dp, 1._dp/216._dp /)
!!     is equivalent to pw_spline3_evaluate_values_g, with
!!     coeff=(/ 27._dp/64._dp, 9._dp/128._dp, 3._dp/256._dp, 1._dp/512._dp /)
!!
!!   INPUTS
!!     - pw_in: the argument for the linear operator
!!     - pw_out: place where the smeared values should be added
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
!!*** **********************************************************************
  SUBROUTINE pw_nn_smear_r(pw_in,pw_out,coeffs,pbc,error)
    TYPE(pw_type), POINTER                   :: pw_in, pw_out
    REAL(KIND=dp), DIMENSION(4), INTENT(in)  :: coeffs
    LOGICAL, OPTIONAL, intent(in)            :: pbc
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_nn_smear_g', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k
    REAL(kind=dp), &
      DIMENSION(-1:1, -1:1, -1:1)            :: weights

    DO k=-1,1
       DO j=-1,1
          DO i=-1,1
             weights(i,j,k)=coeffs(ABS(i)+ABS(j)+ABS(k)+1)
          END DO
       END DO
    END DO

    CALL pw_nn_compose_r(weights=weights,pw_in=pw_in,pw_out=pw_out,&
         pbc=pbc,error=error)
  END SUBROUTINE pw_nn_smear_r
!***************************************************************************


!!****f* pw_spline_utils/pw_nn_deriv_r [1.0] *
!!
!!   NAME
!!     pw_nn_deriv_r
!!
!!   SYNOPSIS
!!     Subroutine pw_nn_deriv_r(pw_in, pw_out, coeffs, idir, error)
!!       Type(pw_type), Pointer:: pw_in, pw_out
!!       Real(Kind=dp), Dimension(3), Intent (IN):: coeffs
!!       Integer:: idir
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine pw_nn_deriv_r
!!
!!   FUNCTION
!!     calculates a nearest neighbor central derivative.
!!     for the x dir:
!!     pw_out%cr3d(i,j,k)=( pw_in(i+1,j,k)-pw_in(i-1,j,k) )*coeff(1)+
!!            ( pw_in(i+1,j(+-)1,k)-pw_in(i-1,j(+-)1,k)+
!!              pw_in(i+1,j,k(+-)1)-pw_in(i-1,j,k(+-)1) )*coeff(2)+
!!            ( pw_in(i+1,j(+-)1,k(+-)1)-pw_in(i-1,j(+-)1,k(+-)1)+
!!              pw_in(i+1,j(+-)1,k(-+)1)-pw_in(i-1,j(+-)1,k(-+)1) )*coeff(3)
!!     periodic boundary conditions are applied
!!
!!   NOTES
!!     with coeff=(/ 2.0_dp/9.0_dp,  1.0_dp/18.0_dp, 1.0_dp/72.0_dp /)
!!     is equivalent to pw_spline3_deriv_r, with
!!     coeff=(/ 9.0_dp/32.0_dp, 3.0_dp/64.0_dp, 1.0_dp/128.0_dp /)
!!     to pw_spline2_deriv_r
!!     coeff=(/ 25._dp/72._dp, 5._dp/144, 1._dp/288._dp /)
!!
!!   INPUTS
!!     - pw_in: the argument for the linear operator
!!     - pw_out: place where the smeared values should be added
!!     - spline_g: on entry the FFT of the values to interpolate as cc3d,
!!       will contain the FFT of the smeared values
!!     - coeffs: array with the coefficent of the front (positive) plane
!!       of the central derivative, ordered with 
!!       the distance from the center: coeffs(1) the coeff of the central
!!       element, coeffs(2) the coeff of the 4 element with distance 1,
!!       coeff(3) the coeff of the 4 elements at distance sqrt(2)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
  SUBROUTINE pw_nn_deriv_r(pw_in,pw_out,coeffs,idir,error)
    TYPE(pw_type), POINTER                   :: pw_in, pw_out
    REAL(KIND=dp), DIMENSION(3), INTENT(in)  :: coeffs
    INTEGER                                  :: idir
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_nn_deriv_r', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, idirVal, j, k
    LOGICAL                                  :: failure
    REAL(kind=dp), &
      DIMENSION(-1:1, -1:1, -1:1)            :: weights

    failure=.FALSE.
    DO k=-1,1
       DO j=-1,1
          DO i=-1,1
             SELECT CASE(idir)
             CASE (1)
                idirVal=i
             CASE (2)
                idirVal=j
             CASE (3)
                idirVal=k
             CASE default
                CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                     routineP,"invalid idir ("//TRIM(cp_to_string(idir))//")",&
                     error,failure)
             END SELECT
             IF (idirVal==0) THEN
                weights(i,j,k)=0.0_dp
             ELSE
                weights(i,j,k)=REAL(idirVal,dp)*coeffs(ABS(i)+ABS(j)+ABS(k))
             END IF
          END DO
       END DO
    END DO
    
    IF (.NOT.failure) THEN
       CALL pw_nn_compose_r(weights=weights,pw_in=pw_in,pw_out=pw_out,&
            error=error)
    END IF
  END SUBROUTINE pw_nn_deriv_r
!***************************************************************************

!!****f* pw_spline_utils/add_coarse2fine [1.0] *
!!
!!   NAME
!!     add_coarse2fine
!!
!!   FUNCTION
!!     low level function that adds a coarse grid (without boundary)
!!     to a fine grid.
!!
!!     It will add to
!!
!!       fine_values(2*coarse_bounds(1,1):2*coarse_bounds(2,1),
!!                   2*coarse_bounds(1,2):2*coarse_bounds(2,2),
!!                   2*coarse_bounds(1,3):2*coarse_bounds(2,3))
!!
!!     using
!!
!!       coarse_coeffs(coarse_bounds(1,1):coarse_bounds(2,1),
!!                     coarse_bounds(1,2):coarse_bounds(2,2),
!!                     coarse_bounds(1,3):coarse_bounds(2,3))
!!
!!     composed with the weights obtained by the direct product of the
!!     1d coefficents weights:
!!
!!     for i,j,k in -3..3
!!        w(i,j,k)=weights_1d(abs(i)+1)*weights_1d(abs(j)+1)*
!!                 weights_1d(abs(k)+1)
!!
!!   NOTES
!!     enormous optimization potential
!!
!!   ARGUMENTS
!!     - coarse_coeffs: 3d array with boundary of size 1 with the values of the
!!       coefficients
!!     - fine_values: 3d array where to add the values due to the
!!       coarse coeffs
!!     - coarse_bounds: bounds of the coarse coefficents
!!     - weights_1d: the weights of the 1d smearing
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE add_coarse2fine(coarse_coeffs,fine_values,coarse_bounds,&
       weights_1d,error)
    INTEGER, DIMENSION(2,3) :: coarse_bounds
    REAL(kind=dp), DIMENSION(coarse_bounds(1,1):coarse_bounds(2,1),&
         coarse_bounds(1,2):coarse_bounds(2,2),&
         coarse_bounds(1,3):coarse_bounds(2,3)), INTENT(in) :: coarse_coeffs
    REAL(kind=dp), DIMENSION(2*coarse_bounds(1,1):2*coarse_bounds(2,1),&
                   2*coarse_bounds(1,2):2*coarse_bounds(2,2),&
                   2*coarse_bounds(1,3):2*coarse_bounds(2,3)),&
                   INTENT(inout) :: fine_values
    REAL(kind=dp), DIMENSION(4), intent(in) :: weights_1d
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='add_coarse2fine',&
         routineP=moduleN//':'//routineN
    INTEGER :: i,j,k,fi,fj,fk,ii,ij,ik
    REAL(kind=dp) :: wk,wj
    
    failure=.FALSE.
    
    IF (.NOT. failure) THEN
!!       DO k=coarse_bounds(1,3)-1,coarse_bounds(2,3)+1
!!          DO ik=-3,3
!!             fk=2*k+ik
!!             wk= weights_1d(ABS(ik)+1)
!!             IF (fk<2*coarse_bounds(1,3).OR.fk>2*coarse_bounds(2,3)) cycle
!!             DO j=coarse_bounds(1,2)-1,coarse_bounds(2,2)+1
!!                DO ij=-3,3
!!                   fj=2*j+ij
!!                   IF (fj<2*coarse_bounds(1,2).OR.fj>2*coarse_bounds(2,2)) cycle
!!                   wj=weights_1d(ABS(ij)+1)*wk
!!                   DO i=coarse_bounds(1,1)-1,coarse_bounds(2,1)+1
!!                      DO ii=-3,3
!!                         fi=2*i+ii
!!                         IF (fi<2*coarse_bounds(1,1).OR.fi>2*coarse_bounds(2,1)) cycle
!!                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
!!                              weights_1d(ABS(ii)+1)*wj*coarse_coeffs(i,j,k)
!!                      END DO
!!                   END DO
!!                END DO
!!             END DO
!!          END DO
!!       END DO

       DO k=coarse_bounds(1,3),coarse_bounds(2,3)
          DO ik=-3,3
             fk=2*k+ik
             wk= weights_1d(ABS(ik)+1)
             IF (fk<2*coarse_bounds(1,3).OR.fk>2*coarse_bounds(2,3)) cycle
             DO j=coarse_bounds(1,2),coarse_bounds(2,2)
                DO ij=-3,3
                   fj=2*j+ij
                   IF (fj<2*coarse_bounds(1,2).OR.fj>2*coarse_bounds(2,2)) cycle
                   wj=weights_1d(ABS(ij)+1)*wk
                   DO i=coarse_bounds(1,1),coarse_bounds(2,1)
                      DO ii=-3,3
                         fi=2*i+ii
                         IF (fi<2*coarse_bounds(1,1).OR.fi>2*coarse_bounds(2,1)) cycle
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              weights_1d(ABS(ii)+1)*wj*coarse_coeffs(i,j,k)
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
    END IF

  END SUBROUTINE add_coarse2fine
!***************************************************************************

!!****f* m/ pw_spline_do_precond *
!!
!!   NAME
!!     pw_spline_do_precond
!!
!!   FUNCTION
!!     applies the preconditioner to the system of equations to find the
!!     coefficents of the spline
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pw_spline_do_precond(preconditioner,in_v,out_v,error)
  TYPE(pw_spline_precond_type), pointer :: preconditioner
  TYPE(pw_type), POINTER :: in_v, out_v
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='pw_spline_do_precond',&
        routineP=moduleN//':'//routineN

  failure=.false.
  
  if (.not. failure) then
     call pw_zero(out_v)
     CALL pw_nn_smear_r(in_v,out_v,spl3_precond1_coeff,pbc=.false.,error=error)
!     CALL pw_copy(in_v,out_v)
  end if
END SUBROUTINE pw_spline_do_precond
!***************************************************************************

!!****f* pw_spline_utils/find_coeffs [1.0] *
!!
!!   NAME
!!     find_coeffs
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     -
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION find_coeffs(values,coeffs,weights,preconditioner, pool, &
     eps_r,eps_x,max_iter,error) result(res)
  TYPE(pw_type), POINTER :: values,coeffs
  type(cp_error_type), optional, intent(inout) :: error
  REAL(kind=dp), DIMENSION(4), intent(in) :: weights
  TYPE(pw_spline_precond_type), pointer :: preconditioner
  TYPE(pw_pool_type), pointer :: pool
  INTEGER, intent(in) :: max_iter
  REAL(kind=dp), INTENT(in) :: eps_r, eps_x
  logical :: res

  logical :: failure
  character(len=*), parameter :: routineN='find_coeffs',&
       routineP=moduleN//':'//routineN
  TYPE(pw_type), POINTER :: r,z,p,Ap
  REAL(kind=dp) :: alpha,beta,r_z,r_z_new,eps_x_att,eps_r_att
  INTEGER :: i,j,k,iter,iiter
  INTEGER, DIMENSION(2,3) :: bo
  failure=.false.

  res=.false.
  IF (.NOT. failure) THEN
     CALL pw_pool_create_pw(pool,r,use_data=REALDATA3D,in_space=REALSPACE,&
          error=error)
     CALL pw_pool_create_pw(pool,z,use_data=REALDATA3D,in_space=REALSPACE,&
          error=error)
     CALL pw_pool_create_pw(pool,p,use_data=REALDATA3D,in_space=REALSPACE,&
          error=error)
     CALL pw_pool_create_pw(pool,Ap,use_data=REALDATA3D,in_space=REALSPACE,&
          error=error)

     DO iiter=1,max_iter,10
        CALL pw_zero(r)
        CALL pw_nn_smear_r(pw_in=coeffs,pw_out=r,coeffs=-weights,pbc=.FALSE.,&
             error=error)
        CALL pw_sumup(values,r)
        CALL pw_spline_do_precond(preconditioner,in_v=r,out_v=z,error=error)
        CALL pw_copy(z,p)
        r_z=pw_integral_ab(r,z)

        DO iter=iiter,MIN(iiter+9,max_iter)
           call pw_zero(Ap)
           CALL pw_nn_smear_r(pw_in=p,pw_out=Ap,coeffs=weights,pbc=.false.,&
                error=error)
           alpha=r_z/pw_integral_ab(Ap,p)

           CALL pw_sumup(p,coeffs,alpha=alpha)

           eps_x_att=alpha*sqrt(pw_integral_aa(p,SQUARE)) ! try to spare if unneded?
           eps_r_att=sqrt(pw_integral_aa(r,SQUARE))
           PRINT *, iter,"eps_x=",eps_x_att,"eps_r=",eps_r_att

           IF (eps_r_att< eps_r.AND.eps_x_att< eps_x) THEN
              res=.true.
              exit
           END IF

!FM           IF (MODULO(iter,10)==0) THEN
!FM              ! rebuilds from scratch now and then
!FM              CALL pw_zero(r)
!FM              CALL pw_nn_smear_r(pw_in=coeffs,pw_out=r,coeffs=-weights,pbc=.FALSE.,&
!FM                   error=error)
!FM              CALL pw_sumup(values,r)
!FM           ELSE
              CALL pw_sumup(Ap,r,alpha=-alpha)
!FM           END IF

           CALL pw_spline_do_precond(preconditioner,in_v=r,out_v=z,error=error)

           r_z_new=pw_integral_ab(r,z)
           beta=r_z_new/r_z
           r_z=r_z_new

           bo=p%pw_grid%bounds_local
           DO k=bo(1,3),bo(2,3)
              DO j=bo(1,2),bo(2,2)
                 DO i=bo(1,1),bo(2,1)
                    p%cr3d(i,j,k)=z%cr3d(i,j,k)+beta*p%cr3d(i,j,k)
                 END DO
              END DO
           END DO

        END DO
     END DO

     CALL pw_pool_give_back_pw(pool,r,error=error)
     CALL pw_pool_give_back_pw(pool,z,error=error)
     CALL pw_pool_give_back_pw(pool,p,error=error)
     CALL pw_pool_give_back_pw(pool,Ap,error=error)
  END IF
END FUNCTION find_coeffs
!***************************************************************************

!!****f* pw_spline_utils/pw_project_s3 [1.0] *
!!
!!   NAME
!!     pw_project_s3
!!
!!   FUNCTION
!!     projects the function from a fine grid to a coarse one
!!
!!   NOTES
!!     extremely slow (but correct) version
!!
!!   ARGUMENTS
!!     - pw_fine_in: the fine grid
!!     - pw_coarse_out: the coarse grid
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pw_project_s3(pw_fine_in,pw_coarse_out,pw_fine_pool,aug_fine_pool,error)
  TYPE(pw_type), pointer :: pw_fine_in
  TYPE(pw_type), pointer :: pw_coarse_out
  TYPE(pw_pool_type), POINTER :: pw_fine_pool,aug_fine_pool
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='pw_project_s3',&
        routineP=moduleN//':'//routineN
  TYPE(pw_type), POINTER :: tmp_fine1,tmp_fine2,tmp_fine_g
  INTEGER :: i,j,k
  INTEGER, DIMENSION(2,3) :: bo

  failure=.false.
  
  IF (.NOT. failure) THEN
     NULLIFY(tmp_fine1,tmp_fine2)
     CALL pw_pool_create_pw(aug_fine_pool,tmp_fine1, use_data=REALDATA3D,&
          in_space=REALSPACE,error=error)
     call pw_zero(tmp_fine1)
!FM     CALL pw_nn_smear_r(pw_in=pw_fine_in,pw_out=tmp_fine1,&
!FM          coeffs=nn2_fine_coeff, pbc=.false.,error=error)
!FM     CALL pw_pool_create_pw(aug_fine_pool,tmp_fine2, use_data=REALDATA3D,&
!FM          in_space=REALSPACE,error=error)
!FM     call pw_zero(tmp_fine2)
!FM     CALL pw_nn_smear_r(pw_in=tmp_fine1,pw_out=tmp_fine2,&
!FM          coeffs=nn2_fine_coeff, pbc=.false.,error=error)
     CALL pw_pool_give_back_pw(aug_fine_pool,tmp_fine2,error=error)
     CALL pw_pool_create_pw(pw_fine_pool,tmp_fine2, use_data=REALDATA3D,&
          in_space=REALSPACE,error=error)
     CALL pw_zero(tmp_fine2)
     bo=pw_coarse_out%pw_grid%bounds_local
     DO k=bo(1,3),bo(2,3)
        DO j=bo(1,2),bo(2,2)
           DO i=bo(1,1),bo(2,1)
              tmp_fine2%cr3d(i,j,k)=tmp_fine1%cr3d(2*i,2*j,2*k)
           END DO
        END DO
     END DO
     CALL pw_pool_give_back_pw(aug_fine_pool,tmp_fine1,error=error)

     CALL pw_pool_create_pw(pw_fine_pool,tmp_fine_g, use_data=COMPLEXDATA1D,&
          in_space=RECIPROCALSPACE,error=error)
     CALL pw_transfer(tmp_fine2,tmp_fine_g)
     CALL pw_spline3_interpolate_values_g(tmp_fine_g,error=error)
     CALL pw_transfer(tmp_fine_g,tmp_fine2)
     bo=pw_coarse_out%pw_grid%bounds_local
     DO k=bo(1,3),bo(2,3)
        DO j=bo(1,2),bo(2,2)
           DO i=bo(1,1),bo(2,1)
              pw_coarse_out%cr3d(i,j,k)=tmp_fine2%cr3d(2*i,2*j,2*k)
           END DO
        END DO
     END DO
     CALL pw_pool_give_back_pw(pw_fine_pool,tmp_fine2,error=error)
     CALL pw_pool_give_back_pw(pw_fine_pool,tmp_fine_g,error=error)
     
  END IF
END SUBROUTINE pw_project_s3
!***************************************************************************

!!****f* pw_spline_utils/pw_inject_s3 [1.0] *
!!
!!   NAME
!!     pw_inject_s3
!!
!!   FUNCTION
!!     injects a function from a coarse grid into a fine one
!!
!!   NOTES
!!     extremely slow (but correct) version
!!
!!   ARGUMENTS
!!     - pw_fine_in: the fine grid
!!     - pw_coarse_out: the coarse grid
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE pw_inject_s3(pw_coarse_in,pw_fine_out,coarse_pool,error)
  TYPE(pw_type), pointer :: pw_coarse_in
  TYPE(pw_type), pointer :: pw_fine_out
  TYPE(pw_pool_type), POINTER :: coarse_pool
  type(cp_error_type), optional, intent(inout) :: error
  
  LOGICAL :: failure, success
  character(len=*), parameter :: routineN='pw_inject_s3',&
        routineP=moduleN//':'//routineN
  TYPE(pw_type), POINTER :: coeffs, pw_tmp ! pw_fine1,tmp_fine_g
  INTEGER :: i,j,k, handle !,side,inc
  INTEGER, DIMENSION(2,3) :: bo
  TYPE(pw_spline_precond_type), pointer :: precond

  failure=.false.
  
  CALL timeset(routineN,"I","",handle)
  IF (.NOT. failure) THEN
     NULLIFY(coeffs)!pw_fine1,tmp_fine_g)
     CALL pw_pool_create_pw(coarse_pool,coeffs, use_data=REALDATA3D,&
          in_space=REALSPACE,error=error)
     CALL pw_copy(pw_coarse_in,coeffs)
     bo=pw_coarse_in%pw_grid%bounds_local

     call pw_zero(coeffs)
     CALL pw_nn_smear_r(pw_coarse_in,coeffs,spl3_aint_coeff,pbc=.false.,error=error)
     success=find_coeffs(values=pw_coarse_in,coeffs=coeffs,&
          weights=spline3_coeffs,preconditioner=precond, pool=coarse_pool, &
          eps_r=1.e-14_dp,eps_x=1.e-14_dp,&
          max_iter=100,error=error)
     PRINT *, "success=",success

!FM     CALL pw_pool_create_pw(coarse_pool,pw_tmp, use_data=REALDATA3D,&
!FM          in_space=REALSPACE,error=error)
!FM     call pw_zero(coeffs)
!FM     call pw_zero(pw_tmp)
!FM     CALL pw_nn_smear_r(pw_coarse_in,pw_tmp,spl3_aint_coeff,pbc=.false.,error=error)
!FM     CALL pw_nn_smear_r(pw_tmp,coeffs,spl3_aint_coeff,pbc=.FALSE.,error=error)
!FM     CALL pw_pool_give_back_pw(coarse_pool,pw_tmp)
!FM     success=find_coeffs(values=pw_coarse_in,coeffs=coeffs,&
!FM          weights=spline3_coeffs,preconditioner=precond, pool=coarse_pool, &
!FM          eps_r=10._dp*EPSILON(0._dp),eps_x=10.*epsilon(0._dp),&
!FM          max_iter=2,error=error)
!FM     PRINT *, "success=",success
!FM     
!FM     call pw_zero(coeffs)
!FM     CALL pw_nn_smear_r(pw_coarse_in,coeffs,spl3_precond1_coeff,pbc=.false.,error=error)
!FM     success=find_coeffs(values=pw_coarse_in,coeffs=coeffs,&
!FM          weights=spline3_coeffs,preconditioner=precond, pool=coarse_pool, &
!FM          eps_r=10._dp*EPSILON(0._dp),eps_x=10.*epsilon(0._dp),&
!FM          max_iter=100,error=error)
!FM     PRINT *, "success=",success
!FM
!FM     CALL pw_pool_create_pw(coarse_pool,pw_tmp, use_data=REALDATA3D,&
!FM          in_space=REALSPACE,error=error)
!FM     call pw_zero(coeffs)
!FM     call pw_zero(pw_tmp)
!FM     CALL pw_nn_smear_r(pw_coarse_in,pw_tmp,spl3_precond1_coeff,pbc=.false.,error=error)
!FM     CALL pw_nn_smear_r(pw_tmp,coeffs,spl3_precond1_coeff,pbc=.FALSE.,error=error)
!FM     CALL pw_pool_give_back_pw(coarse_pool,pw_tmp)
!FM     success=find_coeffs(values=pw_coarse_in,coeffs=coeffs,&
!FM          weights=spline3_coeffs,preconditioner=precond, pool=coarse_pool, &
!FM          eps_r=10._dp*EPSILON(0._dp),eps_x=10.*epsilon(0._dp),&
!FM          max_iter=2,error=error)
!FM     PRINT *, "success=",success

     bo=pw_coarse_in%pw_grid%bounds_local
!FM     DO k=bo(1,3),bo(2,3)
!FM        DO j=bo(1,2),bo(2,2)
!FM           DO i=bo(1,1),bo(2,1)
!FM              pw_fine1%cr3d(i,j,k)=pw_coarse_in%cr3d(i,j,k)
!FM           END DO
!FM        END DO
!FM     END DO
!FM
!FM     do side=0,1
!FM        k=bo(side+1,3)
!FM        inc=2*side-1
!FM        DO j=bo(1,2),bo(2,2)
!FM           DO i=bo(1,1),bo(2,1)
!FM              pw_fine1%cr3d(i,j,k+inc)=pw_fine1%cr3d(i,j,k+inc)+&
!FM                   2._dp*pw_coarse_in%cr3d(i,j,k)-&
!FM                   pw_coarse_in%cr3d(i,j,k-inc)
!FM           END DO
!FM        END DO
!FM     end do
!FM     do side=0,1
!FM        inc=2*side-1
!FM        j=bo(side+1,2)
!FM        DO k=bo(1,3),bo(2,3)
!FM           DO i=bo(1,1),bo(2,1)
!FM              pw_fine1%cr3d(i,j+inc,k)=pw_fine1%cr3d(i,j+inc,k)+&
!FM                   2._dp*pw_coarse_in%cr3d(i,j,k)-&
!FM                   pw_coarse_in%cr3d(i,j-inc,k)
!FM           END DO
!FM        END DO
!FM     end do
!FM     do side=0,1
!FM        inc=2*side-1
!FM        i=bo(side+1,1)
!FM        DO k=bo(1,3),bo(2,3)
!FM           DO j=bo(1,2),bo(2,2)
!FM              pw_fine1%cr3d(i+inc,j,k)=pw_fine1%cr3d(i+inc,j,k)+&
!FM                   2._dp*pw_coarse_in%cr3d(i,j,k)-&
!FM                   pw_coarse_in%cr3d(i-inc,j,k)
!FM           END DO
!FM        END DO
!FM     end do
     
     OPEN(UNIT=21,FILE="coeffs.cube")
     CALL rs_pw_to_cube(coeffs,21,.TRUE.,"coeffs")
     CLOSE(21)

!FM     CALL pw_pool_create_pw(pw_fine_pool,tmp_fine_g, use_data=COMPLEXDATA1D,&
!FM          in_space=RECIPROCALSPACE,error=error)
!FM     CALL pw_transfer(pw_fine1,tmp_fine_g)
!FM     CALL pw_spline3_interpolate_values_g(tmp_fine_g,error=error)
!FM     CALL pw_transfer(tmp_fine_g,pw_fine1)
!FM     CALL pw_pool_give_back_pw(pw_fine_pool,tmp_fine_g,error=error)

     CALL add_coarse2fine(coarse_coeffs=coeffs%cr3d,&
          fine_values=pw_fine_out%cr3d,&
          coarse_bounds=bo,&
          weights_1d=spl3_1d_coeffs,&
          error=error)

     CALL pw_pool_give_back_pw(coarse_pool,coeffs,error=error)
     
  END IF
  CALL timestop(0._dp,handle)
END SUBROUTINE pw_inject_s3
!***************************************************************************

END MODULE pw_spline_utils
