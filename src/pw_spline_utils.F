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
!!     08.2004 removed spline evaluation method using more than 2 read streams
!!             (pw_compose_stripe_rs3), added linear solver based spline 
!!             inversion [fawzi]
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
                                             cp_warning_level,&
                                             cp_logger_get_default_unit_nr
  USE cp_files,                        ONLY: open_file,close_file
  USE kinds,                           ONLY: dp
  USE machine, only: m_flush
  USE mathconstants,                   ONLY: twopi
  USE message_passing,                 ONLY: mp_sendrecv, mp_isend,mp_irecv,&
       mp_waitall, mp_alltoall, mp_sync, mp_comm_compare
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
       pw_pool_give_back_pw, pw_pool_retain, pw_pool_release
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE realspace_grid_types,            ONLY: rs_pw_to_cube
  USE pws, ONLY: pw_integral_ab, pw_integral_aa
  USE input_section_types, ONLY: section_vals_type, section_vals_get_subs_vals,&
       section_vals_val_get
  USE input_keyword_types, only: keyword_type
  USE cp_output_handling, ONLY: cp_should_output, cp_iter_string, cp_iterate,&
       cp_add_iter_level, cp_rm_iter_level, cp_unitnr_for_output, &
       cp_finished_output
  USE string_utilities, only: compress
  USE input_constants, only: no_precond,precond_spl3_aint, precond_spl3_1

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
         spl3_aint_coeff=&
         (/ 46._dp/27._dp,-2._dp/(27._dp),-1._dp/(27._dp*2._dp),&
         -1._dp/(27._dp*8._dp) /),&
         spl3_precond1_coeff=&
         (/ 64._dp/3._dp,-8._dp/3._dp,-1._dp/3._dp, -1._dp/24._dp /),&
         spl3_1d_transf_coeffs=&
         (/ 2._dp/3._dp,23._dp/48._dp,1._dp/6._dp,1._dp/48._dp /)
    
         
    REAL ( KIND = dp ), PUBLIC, PARAMETER, DIMENSION(3) :: spline3_deriv_coeffs=&
         (/ 2.0_dp/9.0_dp,  1.0_dp/18.0_dp, 1.0_dp/72.0_dp /),&
         spline2_deriv_coeffs=&
         (/ 9.0_dp/32.0_dp, 3.0_dp/64.0_dp, 1.0_dp/128.0_dp /),&
         nn10_deriv_coeffs=&
         (/ 25._dp/72._dp, 5._dp/144, 1._dp/288._dp /),&
         nn50_deriv_coeffs=&
         (/ 625._dp/1352._dp, 25._dp/2704._dp,1._dp/5408._dp /),&
         spl3_1d_coeffs0=&
         (/ 1._dp/6_dp,2._dp/3._dp,1._dp/6._dp /),&
         spl3_1d_transf_border1=&
         (/ 0.517977704_dp, 0.464044595_dp, 0.17977701e-1_dp /)
    
    INTEGER, SAVE, PRIVATE :: last_precond_id=0

  PUBLIC :: pw_spline3_interpolate_values_g, pw_spline3_evaluate_values_g,&
       pw_spline3_deriv_g
  PUBLIC :: pw_spline_scale_deriv, pw_nn_smear_g
  PUBLIC :: pw_spline2_interpolate_values_g, pw_spline2_evaluate_values_g,&
       pw_spline2_deriv_g
  PUBLIC :: pw_nn_compose_r, pw_nn_smear_r, pw_nn_deriv_r, pw_inject_s3,&
       pw_project_s3, spl3_nopbc, pw_nn_compose_r_no_pbc, spl3_nopbct

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
     INTEGER :: ref_count,id_nr, kind
     REAL(kind=dp), DIMENSION(4) :: coeffs
     REAL(kind=dp), DIMENSION(3) :: coeffs_1d
     LOGICAL :: sharpen, normalize, pbc
     TYPE(pw_pool_type), pointer :: pool
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

  IF (n_el<1) RETURN
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
     v2=in_val_last
     out_val(n_el)= out_val(n_el)+&
          weights(0)*v0+&
          weights(1)*v1+&
          weights(2)*v2
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
!!     or by inlining suboptimal code is produced (nag))
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
SUBROUTINE pw_nn_compose_r_work(weights,in_val,out_val,pw_in,bo,error)
  REAL(kind=dp), DIMENSION(0:2, 0:2, 0:2)  :: weights
  INTEGER, DIMENSION(2, 3)                 :: bo
  TYPE(pw_type), POINTER                   :: pw_in
  REAL(kind=dp), DIMENSION(bo(1, 1):bo(2, &
       1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, &
       3)), INTENT(inout)                     :: out_val
  REAL(kind=dp), DIMENSION(bo(1, 1):bo(2, &
       1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, &
       3)), INTENT(in)                        :: in_val
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error

  CHARACTER(len=*), PARAMETER :: routineN = "pw_nn_compose_r_work", &
       routineP = moduleN//":"//routineN

  INTEGER                                  :: i, j, jw, k, kw, myj, myk, &
       stat
  INTEGER, DIMENSION(2, 3)                 :: gbo
  INTEGER, DIMENSION(3)                    :: s
  LOGICAL                                  :: failure, has_boundary, &
       yderiv, zderiv
  REAL(kind=dp)                            :: in_val_f, in_val_l
  REAL(kind=dp), DIMENSION(:, :), POINTER  :: l_boundary, tmp, u_boundary

  zderiv=ALL(weights(:,:,1)==0.0_dp)
  yderiv=ALL(weights(:,1,:)==0.0_dp)
  bo=pw_in%pw_grid%bounds_local
  gbo=pw_in%pw_grid%bounds
  DO i=1,3
     s(i)=bo(2,i)-bo(1,i)+1
  END DO
  IF (ANY(s<1)) RETURN
  has_boundary= ANY(pw_in%pw_grid%bounds_local(:,1) /= &
       pw_in%pw_grid%bounds(:,1))
  IF (has_boundary) THEN
     ALLOCATE(l_boundary(bo(1,2):bo(2,2),bo(1,3):bo(2,3)),&
          u_boundary(bo(1,2):bo(2,2),bo(1,3):bo(2,3)),&
          tmp(bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
     DEALLOCATE(tmp,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF

  !$omp parallel do default(none) private(k,kw,myk,j,jw,myj,in_val_f,&
  !$omp     in_val_l) shared(zderiv,yderiv,bo,in_val,out_val,s,l_boundary,&
  !$omp     u_boundary,weights,has_boundary)
  DO k=0,s(3)-1
     DO kw=0,2
        myk=bo(1,3)+MODULO(k+kw-1,s(3))
        IF (zderiv.AND.kw==1) CYCLE
        DO j=0,s(2)-1
           DO jw=0,2
              myj=bo(1,2)+MODULO(j+jw-1,s(2))
              IF (yderiv.AND.jw==1) CYCLE
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
SUBROUTINE pw_nn_compose_r(weights, pw_in, pw_out,error)
    REAL(kind=dp), DIMENSION(0:2, 0:2, 0:2)  :: weights
    TYPE(pw_type), POINTER                   :: pw_in, pw_out
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
          error=error)
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
  SUBROUTINE pw_nn_smear_r(pw_in,pw_out,coeffs,error)
    TYPE(pw_type), POINTER                   :: pw_in, pw_out
    REAL(KIND=dp), DIMENSION(4), INTENT(in)  :: coeffs
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
         error=error)
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

!!****f* pw_spline_utils/add_coarse2fine *
!!
!!   NAME
!!     add_coarse2fine
!!
!!   FUNCTION
!!     low level function that adds a coarse grid
!!     to a fine grid.
!!     If pbc is true periodic boundary conditions are applied
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
!!     coarse looping is continuos, I did not check if keeping the fine looping
!!     contiguos is better.
!!     And I ask myself yet again why, why we use x-slice distribution,
!!     z-slice distribution would be much better performancewise
!!     (and would semplify this code enormously).
!!     fine2coarse has much more understandable parallel part (build up of
!!     send/rcv sizes,... but worse if you have really a lot of processors, 
!!     probabily irrelevant because it is not critical) [fawzi].
!!
!!   ARGUMENTS
!!     - coarse_coeffs_pw: the values of the coefficients
!!     - fine_values_pw: where to add the values due to the
!!       coarse coeffs
!!     - weights_1d: the weights of the 1d smearing
!!     - w_border0: the 1d weight at the border (when pbc is false)
!!     - w_border1: the 1d weights for a point one off the border
!!       (w_border1(1) is the weight of the coefficent at the border)
!!       (used if pbc is false)
!!     - pbc: if periodic boundary conditions should be applied
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE add_coarse2fine(coarse_coeffs_pw,fine_values_pw,&
       weights_1d,w_border0,w_border1,pbc,safe_computation,error)
    TYPE(pw_type), POINTER :: coarse_coeffs_pw, fine_values_pw
    REAL(kind=dp), DIMENSION(4), INTENT(in) :: weights_1d
    REAL(kind=dp), intent(in) :: w_border0
    REAL(kind=dp), DIMENSION(3), INTENT(in) :: w_border1
    LOGICAL, intent(in) :: pbc
    LOGICAL, OPTIONAL, intent(in) :: safe_computation
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure, is_split, has_i_ubound, has_i_lbound, safe_calc
    CHARACTER(len=*), PARAMETER :: routineN='add_coarse2fine',&
         routineP=moduleN//':'//routineN
    INTEGER :: i,j,k,fi,fj,fk,ii,ij,ik,rest_b,s(3), reqs(3), nReq, f_shift(3)
    INTEGER :: rcv_tot_size, ip, send_tot_size, x, xx, x_att, p, p_old, p_lb,&
         p_ub, fi_lb, fi_ub,  n_procs, coarse_slice_size, stat, ss, sf
!FMtmp    integer ::  unit_nr
    REAL(kind=dp) :: wk,wj,wi,v0,v1,v2,v3
    REAL(kind=dp),DIMENSION(3) :: ww0, w_0
    REAL(kind=dp),DIMENSION(4) :: ww1, w_1
    INTEGER, DIMENSION(2,3) :: coarse_bo,coarse_gbo,fine_bo,fine_gbo, my_coarse_bo
    INTEGER, DIMENSION(:), POINTER :: pos_of_x
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_size, send_offset, rcv_size,&
         rcv_offset, sent_size, real_rcv_size
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: send_buf, rcv_buf
    REAL(kind=dp), DIMENSION(:,:,:), POINTER :: coarse_coeffs, fine_values
    REAL(kind=dp), DIMENSION(:,:), POINTER :: u_boundary1,u_boundary2,&
         l_boundary1, l_boundary2, tmp1,tmp2
    TYPE(cp_logger_type), pointer :: logger

    failure=.FALSE.
    safe_calc=.false.
    if (present(safe_computation)) safe_calc=safe_computation
    logger => cp_error_get_logger(error)
!tmpFM    unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    CALL mp_comm_compare(coarse_coeffs_pw%pw_grid%para%group,&
         fine_values_pw%pw_grid%para%group,ii)
    IF (ii>1) then
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF
!tmpFM    WRITE(unit_nr,*) "START add_coarse2fine"
    IF (.NOT. failure) THEN
       my_coarse_bo=coarse_coeffs_pw%pw_grid%bounds_local
       coarse_gbo=coarse_coeffs_pw%pw_grid%bounds
       fine_bo=fine_values_pw%pw_grid%bounds_local
       fine_gbo=fine_values_pw%pw_grid%bounds
       f_shift=fine_gbo(1,:)-2*coarse_gbo(1,:)
       DO j=2,3
          DO i=1,2
             coarse_bo(i,j)=FLOOR((fine_bo(i,j)-f_shift(j))/2.)
          END DO
       END DO
       coarse_bo(1,1)=floor((fine_bo(1,1)-2-f_shift(1))/2.)
       coarse_bo(2,1)=floor((fine_bo(2,1)+3-f_shift(1))/2.)
       is_split= ANY(coarse_gbo(:,1) /= my_coarse_bo(:,1))
       IF (.NOT.is_split.OR..NOT.pbc) THEN
          coarse_bo(1,1)=MAX(coarse_gbo(1,1),coarse_bo(1,1))
          coarse_bo(2,1)=MIN(coarse_gbo(2,1),coarse_bo(2,1))
       END IF
       has_i_ubound=(coarse_gbo(2,1)/=coarse_bo(2,1))
       has_i_lbound=(coarse_gbo(1,1)/=coarse_bo(1,1))

       IF (pbc) THEN
          CPPrecondition(ALL(fine_gbo(1,:)==2*coarse_gbo(1,:)+f_shift),cp_failure_level,routineP,error,failure)
          CPPrecondition(ALL(fine_gbo(2,:)==2*coarse_gbo(2,:)+1+f_shift),cp_failure_level,routineP,error,failure)
       ELSE
          CPPrecondition(ALL(fine_gbo(2,:)==2*coarse_gbo(2,:)+f_shift),cp_failure_level,routineP,error,failure)
          CPPrecondition(ALL(fine_gbo(1,:)==2*coarse_gbo(1,:)+f_shift),cp_failure_level,routineP,error,failure)
       END IF

!tmpFM       WRITE (unit_nr,*) "coarse_gbo",coarse_gbo,"coarse_bo",coarse_bo,"my_coarse_bo",my_coarse_bo,"f_shift",f_shift
       
       coarse_coeffs => coarse_coeffs_pw%cr3d
       DO i=1,3
          s(i)=coarse_gbo(2,i)-coarse_gbo(1,i)+1
       END DO

       ! *** parallel case
       IF (is_split) THEN
          coarse_slice_size=(coarse_bo(2,2)-coarse_bo(1,2)+1)*&
               (coarse_bo(2,3)-coarse_bo(1,3)+1)
          n_procs=coarse_coeffs_pw%pw_grid%para%group_size
          ALLOCATE(send_size(0:n_procs-1),send_offset(0:n_procs-1),&
               sent_size(0:n_procs-1),rcv_size(0:n_procs-1), &
               rcv_offset(0:n_procs-1),real_rcv_size(0:n_procs-1), stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

          ! ** rcv size count

          pos_of_x => coarse_coeffs_pw%pw_grid%para%pos_of_x
          p_old=pos_of_x(coarse_gbo(1,1) &
               +MODULO(coarse_bo(1,1)-coarse_gbo(1,1),s(1)))
          rcv_size=0
          DO x=coarse_bo(1,1),coarse_bo(2,1)
             p=pos_of_x(coarse_gbo(1,1)+MODULO(x-coarse_gbo(1,1),s(1)))
              rcv_size(p)=rcv_size(p)+coarse_slice_size
          END DO

          ! ** send size count

          pos_of_x => fine_values_pw%pw_grid%para%pos_of_x
          sf=fine_gbo(2,1)-fine_gbo(1,1)+1
!tmpFM          WRITE(unit_nr,*) "mepos=",logger%para_env%mepos,"pos_of_x",pos_of_x
!tmpFM          WRITE(unit_nr,*) "fine_bo",fine_values_pw%pw_grid%bounds
!tmpFM          WRITE(unit_nr,*) "fine_bo_local",fine_values_pw%pw_grid%bounds_local
          fi_lb=2*my_coarse_bo(1,1)-3+f_shift(1)
          fi_ub=2*my_coarse_bo(2,1)+3+f_shift(1)
          IF (.not.pbc) THEN
             fi_lb=MAX(fi_lb,fine_gbo(1,1))
             fi_ub=MIN(fi_ub,fine_gbo(2,1))
          ELSE
             fi_ub=MIN(fi_ub,fi_lb+sf-1)
          END IF
          p_old=pos_of_x(fine_gbo(1,1)+MODULO(fi_lb-fine_gbo(1,1),sf))
          p_lb=FLOOR((fi_lb-2-f_shift(1))/2.)
          send_size=0
          DO x=fi_lb,fi_ub
             p=pos_of_x(fine_gbo(1,1)+MODULO(x-fine_gbo(1,1),sf))
!tmpFM             WRITE(unit_nr,*) "x=",x,"p=",p
             IF (p/=p_old) THEN
                p_ub=FLOOR((x-1+3-f_shift(1))/2.)
                
                send_size(p_old)=send_size(p_old)+(MIN(p_ub,my_coarse_bo(2,1)) &
                     -MAX(p_lb,my_coarse_bo(1,1))+1)*coarse_slice_size
                
                IF (pbc) THEN
                   DO xx=p_lb,coarse_gbo(1,1)-1
                      x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
                      IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
                         send_size(p_old)=send_size(p_old)+coarse_slice_size
                      END IF
                   END DO
                   DO xx=coarse_gbo(2,1)+1,p_ub
                      x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
                      IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
                         send_size(p_old)=send_size(p_old)+coarse_slice_size
                      END IF
                   END DO
                END IF

                p_old=p
                p_lb=FLOOR((x-2-f_shift(1))/2.)
             END IF
          END DO
          p_ub=FLOOR((fi_ub+3-f_shift(1))/2.)
          
          send_size(p_old)=send_size(p_old)+(MIN(p_ub,my_coarse_bo(2,1)) &
               -MAX(p_lb,my_coarse_bo(1,1))+1)*coarse_slice_size
          
          IF (pbc) THEN
             DO xx=p_lb,coarse_gbo(1,1)-1
                x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
                IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
                   send_size(p_old)=send_size(p_old)+coarse_slice_size
                END IF
             END DO
             DO xx=coarse_gbo(2,1)+1,p_ub
                x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
                IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
                   send_size(p_old)=send_size(p_old)+coarse_slice_size
                END IF
             END DO
          END IF
!tmpFM          WRITE(unit_nr,*) "send_size",send_size
          ! ** offsets & alloc send-rcv

          send_tot_size=0
          DO ip=0,n_procs-1
             send_offset(ip)=send_tot_size
             send_tot_size=send_tot_size+send_size(ip)
          END DO
          ALLOCATE(send_buf(0:send_tot_size-1),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

!tmpFM          WRITE(unit_nr,*) "rcv_size",rcv_size
          rcv_tot_size=0
          DO ip=0,n_procs-1
             rcv_offset(ip)=rcv_tot_size
             rcv_tot_size=rcv_tot_size+rcv_size(ip)
          END DO
          call cp_assert(rcv_tot_size==&
               (coarse_bo(2,1)-coarse_bo(1,1)+1)*coarse_slice_size,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "Error calculating rcv_tot_size "//&
               CPSourceFileRef,&
               error,failure)
          ALLOCATE(rcv_buf(0:rcv_tot_size-1),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          
          ! ** fill send buffer

!-FM          pos_of_x => fine_values_pw%pw_grid%para%pos_of_x
!-FM          fi_lb=2*my_coarse_bo(1,1)-2+f_shift(1)
!-FM          fi_ub=2*my_coarse_bo(2,1)+3+f_shift(1)
!-FM          IF (.not.pbc) THEN
!-FM             fi_lb=MAX(fi_lb,fine_gbo(1,1))
!-FM             fi_ub=MIN(fi_ub,fine_gbo(2,1))
!-FM          END IF
          p_old=pos_of_x(fine_gbo(1,1)+MODULO(fi_lb-fine_gbo(1,1),sf))
          p_lb=FLOOR((fi_lb-2-f_shift(1))/2.)
          sent_size=send_offset
          ss=my_coarse_bo(2,1)-my_coarse_bo(1,1)+1
          DO x=fi_lb,fi_ub
             p=pos_of_x(fine_gbo(1,1)+MODULO(x-fine_gbo(1,1),sf))
             IF (p/=p_old) THEN
                p_ub=FLOOR((x-1+3-f_shift(1))/2.)
                
                IF (pbc) THEN
                   DO xx=p_lb,coarse_gbo(1,1)-1
                      x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),sf)
                      IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
                         CALL dcopy(coarse_slice_size,&
                              coarse_coeffs(x_att,my_coarse_bo(1,2),&
                              my_coarse_bo(1,3)),ss,send_buf(sent_size(p_old)),1)
                         sent_size(p_old)=sent_size(p_old)+coarse_slice_size
                      END IF
                   END DO
                END IF

                ii=sent_size(p_old)
                DO k=coarse_bo(1,3),coarse_bo(2,3)
                   DO j=coarse_bo(1,2),coarse_bo(2,2)
                      DO i=MAX(p_lb,my_coarse_bo(1,1)),MIN(p_ub,my_coarse_bo(2,1))
                         send_buf(ii)=coarse_coeffs(i,j,k)
                         ii=ii+1
                      END DO
                   END DO
                END DO
                sent_size(p_old)=ii
                
                IF (pbc) THEN
                   DO xx=coarse_gbo(2,1)+1,p_ub
                      x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
                      IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
                         CALL dcopy(coarse_slice_size,&
                              coarse_coeffs(x_att,my_coarse_bo(1,2),&
                              my_coarse_bo(1,3)),ss,&
                              send_buf(sent_size(p_old)),1)
                         sent_size(p_old)=sent_size(p_old)+coarse_slice_size
                      END IF
                   END DO
                END IF

                p_old=p
                p_lb=FLOOR((x-2-f_shift(1))/2.)
             END IF
          END DO
          p_ub=floor((fi_ub+3-f_shift(1))/2.)
          
          IF (pbc) THEN
             DO xx=p_lb,coarse_gbo(1,1)-1
                x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
                IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
                   CALL dcopy(coarse_slice_size,&
                        coarse_coeffs(x_att,my_coarse_bo(1,2),&
                        my_coarse_bo(1,3)),ss,send_buf(sent_size(p_old)),1)
                   sent_size(p_old)=sent_size(p_old)+coarse_slice_size
                END IF
             END DO
          END IF

          ii=sent_size(p_old)
          DO k=coarse_bo(1,3),coarse_bo(2,3)
             DO j=coarse_bo(1,2),coarse_bo(2,2)
                DO i=MAX(p_lb,my_coarse_bo(1,1)),MIN(p_ub,my_coarse_bo(2,1))
                   send_buf(ii)=coarse_coeffs(i,j,k)
                   ii=ii+1
                END DO
             END DO
          END DO
          sent_size(p_old)=ii

          IF (pbc) THEN
             DO xx=coarse_gbo(2,1)+1,p_ub
                x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
                IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
                   CALL dcopy(coarse_slice_size,&
                        coarse_coeffs(x_att,my_coarse_bo(1,2),&
                        my_coarse_bo(1,3)),ss,send_buf(sent_size(p_old)),1)
                   sent_size(p_old)=sent_size(p_old)+coarse_slice_size
                END IF
             END DO
          END IF

!tmpFM          WRITE(unit_nr,*)"sent_size",sent_size,"send_offset",send_offset
          CPPostcondition(ALL(sent_size(:n_procs-2)==send_offset(1:)),cp_failure_level,routineP,error,failure)
!tmpFM          WRITE(unit_nr,*)"ssize",sent_size(n_procs-1),"send_tot_size",send_tot_size
          CPPostcondition(sent_size(n_procs-1)==send_tot_size,cp_failure_level,routineP,error,failure)
!tmpFM          CALL m_flush(unit_nr)
          ! test send/rcv sizes
          CALL mp_alltoall(send_size,real_rcv_size,1,coarse_coeffs_pw%pw_grid%para%group)
          CPAssert(ALL(real_rcv_size==rcv_size),cp_failure_level,routineP,error,failure)
          ! all2all
          call mp_alltoall( sb=send_buf, scount=send_size, sdispl=send_offset,&
               rb=rcv_buf, rcount=rcv_size, rdispl=rcv_offset, &
               group=coarse_coeffs_pw%pw_grid%para%group )

          ! ** reorder rcv buffer
          ! (actually reordering should be needed only with pbc)

          ALLOCATE(coarse_coeffs(coarse_bo(1,1):coarse_bo(2,1),&
               coarse_bo(1,2):coarse_bo(2,2),&
               coarse_bo(1,3):coarse_bo(2,3)),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

          pos_of_x => coarse_coeffs_pw%pw_grid%para%pos_of_x
          p_old=pos_of_x(coarse_gbo(1,1) &
               +MODULO(coarse_bo(1,1)-coarse_gbo(1,1),s(1)))
          sent_size=rcv_offset
          p_lb=coarse_bo(1,1)
          DO x=coarse_bo(1,1),coarse_bo(2,1)
             p=pos_of_x(coarse_gbo(1,1)+MODULO(x-coarse_gbo(1,1),s(1)))
             IF (p/=p_old) THEN
                p_ub=x-1

                ii=sent_size(p_old)
                DO k=coarse_bo(1,3),coarse_bo(2,3)
                   DO j=coarse_bo(1,2),coarse_bo(2,2)
                      DO i=p_lb,p_ub
                         coarse_coeffs(i,j,k)=rcv_buf(ii)
                         ii=ii+1
                      END DO
                   END DO
                END DO
                sent_size(p_old)=ii

                p_lb=x
                p_old=p
             END IF
             rcv_size(p)=rcv_size(p)+coarse_slice_size
          END DO
          p_ub=x-1
          ii=sent_size(p_old)
          DO k=coarse_bo(1,3),coarse_bo(2,3)
             DO j=coarse_bo(1,2),coarse_bo(2,2)
                DO i=p_lb,p_ub
                   coarse_coeffs(i,j,k)=rcv_buf(ii)
                   ii=ii+1
                END DO
             END DO
          END DO
          sent_size(p_old)=ii

          CPPostcondition(ALL(sent_size(0:n_procs-2)==rcv_offset(1:)),cp_failure_level,routineP,error,failure)
          CPPostcondition(sent_size(n_procs-1)==rcv_tot_size,cp_failure_level,routineP,error,failure)
          
          ! dealloc
          DEALLOCATE(send_size,send_offset, rcv_size, rcv_offset, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(send_buf,rcv_buf,real_rcv_size,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       fine_values => fine_values_pw%cr3d
       w_0=(/ weights_1d(3),weights_1d(1),weights_1d(3) /)
       w_1=(/ weights_1d(4),weights_1d(2),weights_1d(2),weights_1d(4) /)

       DO k=coarse_bo(1,3),coarse_bo(2,3)
          DO ik=-3,3
             IF (pbc) THEN
                wk=weights_1d(ABS(ik)+1)
                fk=coarse_gbo(1,3)+MODULO(2*k+ik-coarse_bo(1,3),2*s(3))
             ELSE
                fk=2*k+ik+f_shift(3)
                IF (fk<=fine_bo(1,3)+1.OR.fk>=fine_bo(2,3)-1) THEN
                   IF (fk<fine_bo(1,3).OR.fk>fine_bo(2,3)) CYCLE
                   IF (fk==fine_bo(1,3).OR.fk==fine_bo(2,3)) THEN
                      IF (ik/=0) CYCLE
                      wk=w_border0
                   ELSE IF (fk==2*coarse_bo(1,3)+1+f_shift(3)) THEN
                      SELECT CASE(ik)
                      CASE(1)
                         wk=w_border1(1)
                      CASE(-1)
                         wk=w_border1(2)
                      CASE(-3)
                         wk=w_border1(3)
                      CASE default
                         CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
                         CYCLE
                      END SELECT
                   ELSE
                      SELECT CASE(ik)
                      CASE(3)
                         wk=w_border1(3)
                      CASE(1)
                         wk=w_border1(2)
                      CASE(-1)
                         wk=w_border1(1)
                      CASE default
                         CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
                         CYCLE
                      END SELECT
                   END IF
                ELSE
                   wk=weights_1d(ABS(ik)+1)
                END IF
             END IF
             DO j=coarse_bo(1,2),coarse_bo(2,2)
                DO ij=-3,3
                   IF (pbc) THEN
                      wj=weights_1d(ABS(ij)+1)
                      fj=coarse_bo(1,2)+MODULO(2*j+ij-coarse_bo(1,2),2*s(2))
                   ELSE
                      fj=2*j+ij+f_shift(2)
                      IF (fj<=fine_bo(1,2)+1.OR.fj>=fine_bo(2,2)-1) THEN
                         IF (fj<fine_bo(1,2).OR.fj>fine_bo(2,2)) CYCLE
                         IF (fj==fine_bo(1,2).OR.fj==fine_bo(2,2)) THEN
                            IF (ij/=0) CYCLE
                            wj=w_border0*wk
                         ELSE IF (fj==2*coarse_bo(1,2)+1+f_shift(2)) THEN
                            SELECT CASE(ij)
                            CASE(1)
                               wj=w_border1(1)*wk
                            CASE(-1)
                               wj=w_border1(2)*wk
                            CASE(-3)
                               wj=w_border1(3)*wk
                            CASE default
                               CYCLE
                            END SELECT
                         ELSE
                            SELECT CASE(ij)
                            CASE(-1)
                               wj=w_border1(1)*wk
                            CASE(1)
                               wj=w_border1(2)*wk
                            CASE(3)
                               wj=w_border1(3)*wk
                            CASE default
                               CYCLE
                            END SELECT
                         END IF
                      ELSE
                         wj=weights_1d(ABS(ij)+1)*wk
                      END IF
                   END IF

                   IF (fine_bo(2,1)-fine_bo(1,1)<7.or.safe_calc) THEN
                      DO i=coarse_bo(1,1),coarse_bo(2,1)
                         DO ii=-3,3
                            IF (pbc.and..not.is_split) THEN
                               wi=weights_1d(ABS(ii)+1)
                               fi=coarse_bo(1,2)+MODULO(2*i+ii-coarse_bo(1,1),2*s(1))
                            ELSE
                               fi=2*i+ii+f_shift(1)
                               IF (fi<fine_bo(1,1).OR.fi>fine_bo(2,1)) CYCLE
                               IF (.NOT.has_i_lbound.AND.fi<=fine_gbo(1,1)+1.OR. &
                                    .not.has_i_ubound.and.fi>=fine_gbo(2,1)-1) THEN
                                  IF (fi==fine_gbo(1,1).OR.fi==fine_gbo(2,1)) THEN
                                     IF (ii/=0) CYCLE
                                     wi=w_border0*wj
                                  ELSE IF (fi==fine_gbo(1,1)+1) THEN
                                     SELECT CASE(ii)
                                     CASE(1)
                                        wi=w_border1(1)*wj
                                     CASE(-1)
                                        wi=w_border1(2)*wj
                                     CASE(-3)
                                        wi=w_border1(3)*wj
                                     CASE default
                                        CYCLE
                                     END SELECT
                                  ELSE
                                     SELECT CASE(ii)
                                     CASE(-1)
                                        wi=w_border1(1)*wj
                                     CASE(1)
                                        wi=w_border1(2)*wj
                                     CASE(3)
                                        wi=w_border1(3)*wj
                                     CASE default
                                        CYCLE
                                     END SELECT
                                  END IF
                               ELSE
                                  wi=weights_1d(ABS(ii)+1)*wj
                               END IF
                            END IF
                            fine_values(fi,fj,fk)=&
                                 fine_values(fi,fj,fk)+&
                                 wi*coarse_coeffs(i,j,k)
                         END DO
                      END DO
                   ELSE
                      ww0=wj*w_0
                      ww1=wj*w_1
                      IF (pbc.AND..NOT.is_split) THEN
                         v3=coarse_coeffs(coarse_bo(2,1),j,k)
                         i=coarse_bo(1,1)
                         fi=2*i+f_shift(1)
                         v0= coarse_coeffs(i,j,k)
                         v1= coarse_coeffs(i+1,j,k)
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              ww0(1)*v3+ww0(2)*v0+ww0(3)*v1
                         v2= coarse_coeffs(i+2,j,k)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              ww1(1)*v3+ww1(2)*v0+ww1(3)*v1+ww1(4)*v2
                      ELSE IF (.NOT.has_i_lbound) THEN
                         i=coarse_bo(1,1)
                         fi=2*i+f_shift(1)
                         v0= coarse_coeffs(i,j,k)
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              w_border0*wj*v0
                         v1= coarse_coeffs(i+1,j,k)
                         v2= coarse_coeffs(i+2,j,k)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              wj*(w_border1(1)*v0+w_border1(2)*v1+&
                              w_border1(3)*v2)
                      ELSE
                         i=coarse_bo(1,1)
                         v0= coarse_coeffs(i,j,k)
                         v1= coarse_coeffs(i+1,j,k)
                         v2= coarse_coeffs(i+2,j,k)
                         fi=2*i+f_shift(1)+1
                         CALL cp_assert(fi+1==fine_bo(1,1).OR.&
                              fi+2==fine_bo(1,1),cp_failure_level,&
                              cp_assertion_failed,routineN,&
                              "unexpected start index "//&
                              TRIM(cp_to_string(coarse_bo(1,1)))//" "//&
                              TRIM(cp_to_string(fi)),error,failure)
                      END IF
                      fi=fi+1
                      IF (fi>=fine_bo(1,1)) THEN
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              ww0(1)*v0+ww0(2)*v1+&
                              ww0(3)*v2
                      ELSE
                         CPPrecondition(fi+1==fine_bo(1,1),cp_failure_level,routineP,error,failure)
                      END IF
                      DO i=coarse_bo(1,1)+3,floor((fine_bo(2,1)-f_shift(1))/2.)-3,4
                         v3=coarse_coeffs(i,j,k)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              (ww1(1)*v0+ww1(2)*v1+&
                              ww1(3)*v2+ww1(4)*v3)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              (ww0(1)*v1+ww0(2)*v2+&
                              ww0(3)*v3)
                         v0=coarse_coeffs(i+1,j,k)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              (ww1(4)*v0+ww1(1)*v1+&
                              ww1(2)*v2+ww1(3)*v3)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              (ww0(1)*v2+ww0(2)*v3+&
                              ww0(3)*v0)
                         v1=coarse_coeffs(i+2,j,k)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              (ww1(3)*v0+ww1(4)*v1+&
                              ww1(1)*v2+ww1(2)*v3)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              (ww0(1)*v3+ww0(2)*v0+&
                              ww0(3)*v1)
                         v2=coarse_coeffs(i+3,j,k)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              (ww1(2)*v0+ww1(3)*v1+&
                              ww1(4)*v2+ww1(1)*v3)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              (ww0(1)*v0+ww0(2)*v1+&
                              ww0(3)*v2)
                      END DO
                      rest_b=MODULO(floor((fine_bo(2,1)-f_shift(1))/2.)-coarse_bo(1,1)-3+1,4)
                      IF (rest_b>0) THEN
                         i=floor((fine_bo(2,1)-f_shift(1))/2.)-rest_b+1
                         v3=coarse_coeffs(i,j,k)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              (ww1(1)*v0+ww1(2)*v1+&
                              ww1(3)*v2+ww1(4)*v3)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              (ww0(1)*v1+ww0(2)*v2+&
                              ww0(3)*v3)
                         IF (rest_b>1) THEN
                            v0=coarse_coeffs(i+1,j,k)
                            fi=fi+1
                            fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                 (ww1(4)*v0+ww1(1)*v1+&
                                 ww1(2)*v2+ww1(3)*v3)
                            fi=fi+1
                            fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                 (ww0(1)*v2+ww0(2)*v3+&
                                 ww0(3)*v0)
                            IF (rest_b>2) THEN
                               v1=coarse_coeffs(i+2,j,k)
                               fi=fi+1
                               fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                    (ww1(3)*v0+ww1(4)*v1+&
                                    ww1(1)*v2+ww1(2)*v3)
                               fi=fi+1
                               fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                    (ww0(1)*v3+ww0(2)*v0+&
                                    ww0(3)*v1)
                               IF (pbc.and..not.is_split) THEN
                                  v2=coarse_coeffs(coarse_bo(1,1),j,k)
                                  fi=fi+1
                                  fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                       ww1(1)*v3+ww1(2)*v0+ww1(3)*v1+ww1(4)*v2
                                  fi=fi+1
                                  fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                       ww0(1)*v0+ww0(2)*v1+ww0(3)*v2
                                  v3=coarse_coeffs(coarse_bo(1,1)+1,j,k)
                                  fi=fi+1
                                  fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                       ww1(1)*v0+ww1(2)*v1+ww1(3)*v2+ww1(4)*v3
                               ELSE IF (has_i_ubound) THEN
                                  v2=coarse_coeffs(i+3,j,k)
                                  fi=fi+1
                                  fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                       ww1(1)*v3+ww1(2)*v0+ww1(3)*v1+ww1(4)*v2
                                  fi=fi+1
                                  fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                       ww0(1)*v0+ww0(2)*v1+ww0(3)*v2
                                  IF (fi+1==fine_bo(2,1)) THEN
                                     v3=coarse_coeffs(i+4,j,k)
                                     fi=fi+1
                                     fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                          ww1(1)*v0+ww1(2)*v1+ww1(3)*v2+ww1(4)*v3
                                  END IF
                               ELSE
                                  fi=fi+1
                                  fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                       wj*(w_border1(3)*v3+w_border1(2)*v0+&
                                       w_border1(1)*v1)
                                  fi=fi+1
                                  fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                       w_border0*wj*v1
                               END IF
                            ELSE IF (pbc.and..not.is_split) THEN
                               v1=coarse_coeffs(coarse_bo(1,1),j,k)
                               fi=fi+1
                               fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                    ww1(1)*v2+ww1(2)*v3+ww1(3)*v0+ww1(4)*v1
                               fi=fi+1
                               fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                    ww0(1)*v3+ww0(2)*v0+ww0(3)*v1
                               v2=coarse_coeffs(coarse_bo(1,1)+1,j,k)
                               fi=fi+1
                               fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                    ww1(1)*v3+ww1(2)*v0+ww1(3)*v1+ww1(4)*v2
                            ELSE IF(has_i_ubound) THEN
                               v1=coarse_coeffs(i+2,j,k)
                               fi=fi+1
                               fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                    ww1(1)*v2+ww1(2)*v3+ww1(3)*v0+ww1(4)*v1
                               fi=fi+1
                               fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                    ww0(1)*v3+ww0(2)*v0+ww0(3)*v1
                               IF (fi+1==fine_bo(2,1)) THEN
                                  v2=coarse_coeffs(i+3,j,k)
                                  fi=fi+1
                                  fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                       ww1(1)*v3+ww1(2)*v0+ww1(3)*v1+ww1(4)*v2
                               END IF
                            ELSE
                               fi=fi+1
                               fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                    wj*(w_border1(3)*v2+w_border1(2)*v3+&
                                    w_border1(1)*v0)
                               fi=fi+1
                               fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                    w_border0*wj*v0
                            END IF
                         ELSE if (pbc) then
                            v0=coarse_coeffs(coarse_bo(1,1),j,k)
                            fi=fi+1
                            fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                 ww1(1)*v1+ww1(2)*v2+ww1(3)*v3+ww1(4)*v0
                            fi=fi+1
                            fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                 ww0(1)*v2+ww0(2)*v3+ww0(3)*v0
                            v1=coarse_coeffs(coarse_bo(1,1)+1,j,k)
                            fi=fi+1
                            fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                 ww1(1)*v2+ww1(2)*v3+ww1(3)*v0+ww1(4)*v1
                         ELSE IF (has_i_ubound) THEN
                            v0=coarse_coeffs(i+1,j,k)
                            fi=fi+1
                            fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                 ww1(1)*v1+ww1(2)*v2+ww1(3)*v3+ww1(4)*v0
                            fi=fi+1
                            fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                 ww0(1)*v2+ww0(2)*v3+ww0(3)*v0
                            IF (fi+1==fine_bo(2,1)) THEN
                               v1=coarse_coeffs(i+2,j,k)
                               fi=fi+1
                               fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                    ww1(1)*v2+ww1(2)*v3+ww1(3)*v0+ww1(4)*v1
                            END IF
                         ELSE
                            fi=fi+1
                            fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                 wj*(w_border1(3)*v1+w_border1(2)*v2+&
                                 w_border1(1)*v3)
                            fi=fi+1
                            fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                 w_border0*wj*v3
                         END IF
                      ELSE if (pbc.and..not.is_split) then
                         v3=coarse_coeffs(coarse_bo(1,1),j,k)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              ww1(1)*v0+ww1(2)*v1+ww1(3)*v2+ww1(4)*v3
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              ww0(1)*v1+ww0(2)*v2+ww0(3)*v3
                         v0=coarse_coeffs(coarse_bo(1,1)+1,j,k)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              ww1(1)*v1+ww1(2)*v2+ww1(3)*v3+ww1(4)*v0
                      ELSE IF (has_i_ubound) then
                         v3=coarse_coeffs(i,j,k)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              ww1(1)*v0+ww1(2)*v1+ww1(3)*v2+ww1(4)*v3
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              ww0(1)*v1+ww0(2)*v2+ww0(3)*v3
                         IF (fi+1==fine_bo(2,1)) THEN
                            v0=coarse_coeffs(i+1,j,k)
                            fi=fi+1
                            fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                                 ww1(1)*v1+ww1(2)*v2+ww1(3)*v3+ww1(4)*v0
                         END IF
                      ELSE
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              wj*(w_border1(3)*v0+w_border1(2)*v1+&
                              w_border1(1)*v2)
                         fi=fi+1
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              w_border0*wj*v2
                      END IF
                      CPPostcondition(fi==fine_bo(2,1),cp_failure_level,routineP,error,failure)
                   END IF
                END DO
             END DO
          END DO
       END DO

       IF (is_split) THEN
          DEALLOCATE(coarse_coeffs,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
!tmpFM    WRITE(unit_nr,*) "END add_coarse2fine"
    call mp_sync ( coarse_coeffs_pw%pw_grid%para%group )
  END SUBROUTINE add_coarse2fine
!***************************************************************************



!!****f* pw_spline_utils/add_fine2coarse [1.0] *
!!
!!   NAME
!!     add_fine2coarse
!!
!!   FUNCTION
!!     low level function that adds a coarse grid (without boundary)
!!     to a fine grid.
!!
!!     It will add to
!!
!!       coarse_coeffs(coarse_bounds(1,1):coarse_bounds(2,1),
!!                     coarse_bounds(1,2):coarse_bounds(2,2),
!!                     coarse_bounds(1,3):coarse_bounds(2,3))
!!
!!     using
!!
!!       fine_values(2*coarse_bounds(1,1):2*coarse_bounds(2,1),
!!                   2*coarse_bounds(1,2):2*coarse_bounds(2,2),
!!                   2*coarse_bounds(1,3):2*coarse_bounds(2,3))
!!
!!     composed with the weights obtained by the direct product of the
!!     1d coefficents weights:
!!
!!     for i,j,k in -3..3
!!        w(i,j,k)=weights_1d(abs(i)+1)*weights_1d(abs(j)+1)*
!!                 weights_1d(abs(k)+1)
!!
!!   NOTES
!!     see coarse2fine for some relevant notes
!!
!!   ARGUMENTS
!!     - coarse_coeffs: 3d array with boundary of size 1 with the values of the
!!       coefficients
!!     - fine_values: 3d array where to add the values due to the
!!       coarse coeffs
!!     - coarse_bounds: bounds of the coarse coefficents
!!     - weights_1d: the weights of the 1d smearing
!!     - w_border0: the 1d weight at the border
!!     - w_border1: the 1d weights for a point one off the border
!!       (w_border1(1) is the weight of the coefficent at the border)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE add_fine2coarse(fine_values_pw,coarse_coeffs_pw,&
       weights_1d,w_border0,w_border1,pbc,safe_computation,error)
    TYPE(pw_type), POINTER :: coarse_coeffs_pw, fine_values_pw
    REAL(kind=dp), DIMENSION(4), INTENT(in) :: weights_1d
    REAL(kind=dp), intent(in) :: w_border0
    REAL(kind=dp), DIMENSION(3), INTENT(in) :: w_border1
    LOGICAL, intent(in) :: pbc
    LOGICAL, OPTIONAL, intent(in) :: safe_computation
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    REAL(kind=dp), DIMENSION(:,:,:), POINTER :: fine_values, coarse_coeffs
    CHARACTER(len=*), PARAMETER :: routineN = 'add_fine2coarse', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: fi, fj, fk, i, ij, ik, j, k, &
         rest_b, ii, s(3), f_shift(3)
!FMtmp    integer :: unit_nr
    INTEGER :: rcv_tot_size, ip, send_tot_size, x, xx, x_att, p, p_old, p_lb,&
         p_ub, fi_lb, fi_ub,  n_procs, coarse_slice_size, stat, ss,&
         sf
    INTEGER, DIMENSION(2,3) :: coarse_bo,coarse_gbo,fine_bo, fine_gbo, my_coarse_bo
    LOGICAL                                  :: failure, is_split, has_i_ubound, &
         has_i_lbound, local_data, safe_calc
    REAL(kind=dp)                            :: vv0, vv1, vv2, vv3,vv4,vv5,vv6,vv7, wj, wk, wi
    REAL(kind=dp), DIMENSION(3)              :: w_0, ww0
    REAL(kind=dp), DIMENSION(4)              :: w_1, ww1
    INTEGER, DIMENSION(:), POINTER :: pos_of_x
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_size, send_offset, rcv_size,&
         rcv_offset, sent_size, pp_lb, pp_ub, real_rcv_size
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: send_buf, rcv_buf
    TYPE(cp_logger_type), pointer :: logger

    failure=.FALSE.

    IF (.NOT. failure) THEN
       safe_calc=.FALSE.
       IF (PRESENT(safe_computation)) safe_calc=safe_computation
       logger => cp_error_get_logger(error)
!tmpFM       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)

       my_coarse_bo=coarse_coeffs_pw%pw_grid%bounds_local
       coarse_gbo=coarse_coeffs_pw%pw_grid%bounds
       fine_bo=fine_values_pw%pw_grid%bounds_local
       fine_gbo=fine_values_pw%pw_grid%bounds
       f_shift=fine_gbo(1,:)-2*coarse_gbo(1,:)
       is_split= ANY(coarse_gbo(:,1) /= my_coarse_bo(:,1))
       coarse_bo=my_coarse_bo
       coarse_bo(1,1)=FLOOR(REAL(fine_bo(1,1)-f_shift(1),dp)/2._dp)-1
       coarse_bo(2,1)=FLOOR(REAL(fine_bo(2,1)+1-f_shift(1),dp)/2._dp)+1
       IF (.NOT.is_split.OR..NOT.pbc) THEN
          coarse_bo(1,1)=MAX(coarse_gbo(1,1),coarse_bo(1,1))
          coarse_bo(2,1)=MIN(coarse_gbo(2,1),coarse_bo(2,1))
       END IF
       has_i_ubound=(coarse_gbo(2,1)/=coarse_bo(2,1))
       has_i_lbound=(coarse_gbo(1,1)/=coarse_bo(1,1))
!tmpFM       WRITE(unit_nr,*) "START add_fine2coarse"
!tmpFM       WRITE(unit_nr,*)"coarse_bo",coarse_bo,"my_coarse_bo",my_coarse_bo
!tmpFM       WRITE(unit_nr,*) "has_i_ubound",has_i_ubound,"has_i_lbound",has_i_lbound
!tmpFM       WRITE(unit_nr,*) "is_split",is_split
       IF (pbc) THEN
          CPPrecondition(ALL(fine_gbo(1,:)==2*coarse_gbo(1,:)+f_shift),cp_failure_level,routineP,error,failure)
          CPPrecondition(ALL(fine_gbo(2,:)==2*coarse_gbo(2,:)+f_shift+1),cp_failure_level,routineP,error,failure)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ELSE
          CPPrecondition(ALL(fine_gbo(2,:)==2*coarse_gbo(2,:)+f_shift),cp_failure_level,routineP,error,failure)
          CPPrecondition(ALL(fine_gbo(1,:)==2*coarse_gbo(1,:)+f_shift),cp_failure_level,routineP,error,failure)
       END IF
       CPPrecondition(coarse_gbo(2,1)-coarse_gbo(1,2)>1,cp_failure_level,routineP,error,failure)
       local_data=is_split ! ANY(coarse_bo/=my_coarse_bo)
       IF (local_data) THEN
          ALLOCATE(coarse_coeffs(coarse_bo(1,1):coarse_bo(2,1),&
               coarse_bo(1,2):coarse_bo(2,2),&
               coarse_bo(1,3):coarse_bo(2,3)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          coarse_coeffs=0._dp
       ELSE
          coarse_coeffs => coarse_coeffs_pw%cr3d
       END IF

       fine_values => fine_values_pw%cr3d
       w_0=(/ weights_1d(3),weights_1d(1),weights_1d(3) /)
       w_1=(/ weights_1d(4),weights_1d(2),weights_1d(2),weights_1d(4) /)

       DO i=1,3
          s(i)=coarse_gbo(2,i)-coarse_gbo(1,i)+1
       END DO
       IF (ANY(s<1)) RETURN

       DO k=coarse_bo(1,3),coarse_bo(2,3)
          DO ik=-3,3
             fk=2*k+ik+f_shift(3)
             IF (fk<=fine_bo(1,3)+1.OR.fk>=fine_bo(2,3)-1) THEN
                IF (fk<fine_bo(1,3).OR.fk>fine_bo(2,3)) CYCLE
                IF (fk==fine_bo(1,3).OR.fk==fine_bo(2,3)) THEN
                   IF (ik/=0) CYCLE
                   wk=w_border0
                ELSE IF (fk==fine_bo(1,3)+1) THEN
                   SELECT CASE(ik)
                   CASE(1)
                      wk=w_border1(1)
                   CASE(-1)
                      wk=w_border1(2)
                   CASE(-3)
                      wk=w_border1(3)
                   CASE default
                      CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
                      CYCLE
                   END SELECT
                ELSE
                   SELECT CASE(ik)
                   CASE(3)
                      wk=w_border1(3)
                   CASE(1)
                      wk=w_border1(2)
                   CASE(-1)
                      wk=w_border1(1)
                   CASE default
                      CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
                      CYCLE
                   END SELECT
                END IF
             ELSE
                wk=weights_1d(ABS(ik)+1)
             END IF
             DO j=coarse_bo(1,2),coarse_bo(2,2)
                DO ij=-3,3
                   fj=2*j+ij+f_shift(2)
                   IF (fj<=fine_bo(1,2)+1.OR.fj>=fine_bo(2,2)-1) THEN
                      IF (fj<fine_bo(1,2).OR.fj>fine_bo(2,2)) CYCLE
                      IF (fj==fine_bo(1,2).OR.fj==fine_bo(2,2)) THEN
                         IF (ij/=0) CYCLE
                         wj=w_border0*wk
                      ELSE IF (fj==fine_bo(1,2)+1) THEN
                         SELECT CASE(ij)
                         CASE(1)
                            wj=w_border1(1)*wk
                         CASE(-1)
                            wj=w_border1(2)*wk
                         CASE(-3)
                            wj=w_border1(3)*wk
                         CASE default
                            CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
                            CYCLE
                         END SELECT
                      ELSE
                         SELECT CASE(ij)
                         CASE(-1)
                            wj=w_border1(1)*wk
                         CASE(1)
                            wj=w_border1(2)*wk
                         CASE(3)
                            wj=w_border1(3)*wk
                         CASE default
                            CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
                            CYCLE
                         END SELECT
                      END IF
                   ELSE
                      wj=weights_1d(ABS(ij)+1)*wk
                   END IF

                   IF (coarse_bo(2,1)-coarse_bo(1,1)<7.or.safe_calc) THEN
                      DO i=coarse_bo(1,1),coarse_bo(2,1)
                         DO ii=-3,3
                            IF (pbc.and..not.is_split) THEN
                               wk=weights_1d(ABS(ii)+1)
                               fi=coarse_bo(1,1)+MODULO(2*i+ii-coarse_bo(1,1),2*s(1))+f_shift(1)
                            ELSE
                               fi=2*i+ii+f_shift(1)
                               IF (fi<fine_bo(1,1).OR.fi>fine_bo(2,1)) CYCLE
                               IF (.NOT.has_i_lbound.AND.fi<=fine_gbo(1,1)+1.OR. &
                                    .NOT.has_i_ubound.AND.fi>=fine_gbo(2,1)-1) THEN
                                  IF (fi==fine_bo(1,1).OR.fi==fine_bo(2,1)) THEN
                                     IF (ii/=0) CYCLE
                                     wi=w_border0*wj
                                  ELSE IF (fi==fine_bo(1,1)+1) THEN
                                     SELECT CASE(ii)
                                     CASE(1)
                                        wi=w_border1(1)*wj
                                     CASE(-1)
                                        wi=w_border1(2)*wj
                                     CASE(-3)
                                        wi=w_border1(3)*wj
                                     CASE default
                                        CYCLE
                                     END SELECT
                                  ELSE
                                     SELECT CASE(ii)
                                     CASE(-1)
                                        wi=w_border1(1)*wj
                                     CASE(1)
                                        wi=w_border1(2)*wj
                                     CASE(3)
                                        wi=w_border1(3)*wj
                                     CASE default
                                        CYCLE
                                     END SELECT
                                  END IF
                               ELSE
                                  wi=weights_1d(ABS(ii)+1)*wj
                               END IF
                            END IF
                            coarse_coeffs(i,j,k)=&
                                 coarse_coeffs(i,j,k)+&
                                 wi*fine_values(fi,fj,fk)
                         END DO
                      END DO
                   ELSE
                      ww0=wj*w_0
                      ww1=wj*w_1
                      IF (pbc.AND..NOT.is_split) THEN
                         i=coarse_bo(1,1)-1
                         vv2=fine_values(fine_bo(2,1)-2,fj,fk)
                         vv3=fine_values(fine_bo(2,1)-1,fj,fk)
                         vv4=fine_values(fine_bo(2,1),fj,fk)
                         fi=fine_bo(1,1)
                         vv5=fine_values(fi,fj,fk)
                         fi=fi+1
                         vv6=fine_values(fi,fj,fk)
                         fi=fi+1
                         vv7=fine_values(fi,fj,fk)
                         coarse_coeffs(i+1,j,k)=coarse_coeffs(i+1,j,k)&
                              +ww1(4)*vv2+ww0(3)*vv3+ww1(3)*vv4+ww0(2)*vv5+ww1(2)*vv6+ww0(1)*vv7
                         coarse_coeffs(i+2,j,k)=coarse_coeffs(i+2,j,k)&
                              +ww1(4)*vv4+ww0(3)*vv5+ww1(3)*vv6+ww0(2)*vv7
                         coarse_coeffs(i+3,j,k)=coarse_coeffs(i+3,j,k)&
                              +ww1(4)*vv6+ww0(3)*vv7
                      ELSE IF (has_i_lbound) THEN
!-FM                         i=coarse_bo(1,1)
!-FM                         fi=fine_bo(1,1)-1
!-FM                         IF (i+1==FLOOR((fine_bo(1,1)+1)/2._dp)) THEN
!-FM                            fi=fi+1
!-FM                            vv0=fine_values(fi,fj,fk)
!-FM                         END IF
!-FM                         coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)+&
!-FM                              vv0*ww0(3)*
                      ELSE
                         i=coarse_bo(1,1)
                         fi=2*i+f_shift(1)
                         vv0=fine_values(fi,fj,fk)
                         fi=fi+1
                         vv1=fine_values(fi,fj,fk)
                         fi=fi+1
                         vv2=fine_values(fi,fj,fk)
                         coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)+&
                              (vv0*w_border0+vv1*w_border1(1))*wj+vv2*ww0(1)
                         coarse_coeffs(i+1,j,k)=coarse_coeffs(i+1,j,k)+&
                              wj*w_border1(2)*vv1+ww0(2)*vv2
                         coarse_coeffs(i+2,j,k)=coarse_coeffs(i+2,j,k)+&
                              wj*w_border1(3)*vv1+ww0(3)*vv2
                      END IF
!-FM                      ubound=MIN(FLOOR((fine_bo(2,1)-f_shift(1))/2._dp),coarse_bo(2,1))
                      DO i=coarse_bo(1,1)+3,FLOOR((fine_bo(2,1)-f_shift(1))/2._dp)-3,4
                         fi=fi+1
                         vv0=fine_values(fi,fj,fk)
                         fi=fi+1
                         vv1=fine_values(fi,fj,fk)
                         fi=fi+1
                         vv2=fine_values(fi,fj,fk)
                         fi=fi+1
                         vv3=fine_values(fi,fj,fk)
                         fi=fi+1
                         vv4=fine_values(fi,fj,fk)
                         fi=fi+1
                         vv5=fine_values(fi,fj,fk)
                         fi=fi+1
                         vv6=fine_values(fi,fj,fk)
                         fi=fi+1
                         vv7=fine_values(fi,fj,fk)
                         coarse_coeffs(i-3,j,k)=coarse_coeffs(i-3,j,k)&
                              +ww1(1)*vv0
                         coarse_coeffs(i-2,j,k)=coarse_coeffs(i-2,j,k)&
                              +ww1(2)*vv0+ww0(1)*vv1+ww1(1)*vv2
                         coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                              +ww1(3)*vv0+ww0(2)*vv1+ww1(2)*vv2+ww0(1)*vv3+ww1(1)*vv4
                         coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                              +ww1(4)*vv0+ww0(3)*vv1+ww1(3)*vv2+ww0(2)*vv3+ww1(2)*vv4+ww0(1)*vv5+ww1(1)*vv6
                         coarse_coeffs(i+1,j,k)=coarse_coeffs(i+1,j,k)&
                              +ww1(4)*vv2+ww0(3)*vv3+ww1(3)*vv4+ww0(2)*vv5+ww1(2)*vv6+ww0(1)*vv7
                         coarse_coeffs(i+2,j,k)=coarse_coeffs(i+2,j,k)&
                              +ww1(4)*vv4+ww0(3)*vv5+ww1(3)*vv6+ww0(2)*vv7
                         coarse_coeffs(i+3,j,k)=coarse_coeffs(i+3,j,k)&
                              +ww1(4)*vv6+ww0(3)*vv7
                      END DO
                      call cp_assert(FLOOR((fine_bo(2,1)-f_shift(1))/2._dp)-coarse_bo(1,1)>=4,&
                           cp_failure_level,cp_assertion_failed,routineP,&
                           "FLOOR((fine_bo(2,1)-f_shift(1))/2._dp)-coarse_bo(1,1)>=4",&
                           error,failure)
                      rest_b=MODULO(FLOOR((fine_bo(2,1)-f_shift(1))/2._dp)-coarse_bo(1,1)-6,4)
                      !PRINT *,"old i=",i
                      i=FLOOR((fine_bo(2,1)-f_shift(1))/2._dp)-3-rest_b+4
                      !PRINT *, "i=",i,"rest_b=",rest_b,"fi=",fi
                      CPPrecondition(fi==(i-2)*2+f_shift(1),cp_failure_level,routineP,error,failure)
                      IF (rest_b>0) THEN
                         fi=fi+1
                         vv0=fine_values(fi,fj,fk)
                         fi=fi+1
                         vv1=fine_values(fi,fj,fk)
                         coarse_coeffs(i-3,j,k)=coarse_coeffs(i-3,j,k)&
                              +ww1(1)*vv0
                         coarse_coeffs(i-2,j,k)=coarse_coeffs(i-2,j,k)&
                              +ww1(2)*vv0+ww0(1)*vv1
                         coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                              +ww1(3)*vv0+ww0(2)*vv1
                         coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                              +ww1(4)*vv0+ww0(3)*vv1
                         IF (rest_b>1) THEN
                            fi=fi+1
                            vv2=fine_values(fi,fj,fk)
                            fi=fi+1
                            vv3=fine_values(fi,fj,fk)
                            coarse_coeffs(i-2,j,k)=coarse_coeffs(i-2,j,k)&
                                 +ww1(1)*vv2
                            coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                                 +ww1(2)*vv2+ww0(1)*vv3
                            coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                                 +ww1(3)*vv2+ww0(2)*vv3
                            coarse_coeffs(i+1,j,k)=coarse_coeffs(i+1,j,k)&
                                 +ww1(4)*vv2+ww0(3)*vv3
                            IF (rest_b>2) THEN
                               fi=fi+1
                               vv4=fine_values(fi,fj,fk)
                               fi=fi+1
                               vv5=fine_values(fi,fj,fk)
                               fi=fi+1
                               vv6=fine_values(fi,fj,fk)
                               fi=fi+1
                               vv7=fine_values(fi,fj,fk)
                               coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                                    +ww1(1)*vv4
                               coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                                    +ww1(2)*vv4+ww0(1)*vv5+wj*w_border1(3)*vv6
                               coarse_coeffs(i+1,j,k)=coarse_coeffs(i+1,j,k)&
                                    +ww1(3)*vv4+ww0(2)*vv5+wj*w_border1(2)*vv6
                               coarse_coeffs(i+2,j,k)=coarse_coeffs(i+2,j,k)&
                                    +ww1(4)*vv4+ww0(3)*vv5+wj*w_border1(1)*vv6+w_border0*wj*vv7
                            ELSE
                               fi=fi+1
                               vv4=fine_values(fi,fj,fk)
                               fi=fi+1
                               vv5=fine_values(fi,fj,fk)
                               coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                                    +wj*w_border1(3)*vv4
                               coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                                    +wj*w_border1(2)*vv4
                               coarse_coeffs(i+1,j,k)=coarse_coeffs(i+1,j,k)&
                                    +wj*(w_border1(1)*vv4+w_border0*vv5)
                            END IF
                         ELSE
                            fi=fi+1
                            vv2=fine_values(fi,fj,fk)
                            fi=fi+1
                            vv3=fine_values(fi,fj,fk)
                            coarse_coeffs(i-2,j,k)=coarse_coeffs(i-2,j,k)&
                                 +wj*w_border1(3)*vv2
                            coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                                 +wj*w_border1(2)*vv2
                            coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                                 +wj*(w_border1(1)*vv2+w_border0*vv3)
                         END IF
                      ELSE
                         fi=fi+1
                         vv0=fine_values(fi,fj,fk)
                         fi=fi+1
                         vv1=fine_values(fi,fj,fk)

                         coarse_coeffs(i-3,j,k)=coarse_coeffs(i-3,j,k)&
                              +wj*w_border1(3)*vv0
                         coarse_coeffs(i-2,j,k)=coarse_coeffs(i-2,j,k)&
                              +wj*w_border1(2)*vv0
                         coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                              +wj*(w_border1(1)*vv0+w_border0*vv1)
                      END IF
                      CPPostcondition(fi==fine_bo(2,1),cp_failure_level,routineP,error,failure)
                   END IF
                END DO
             END DO
          END DO
       END DO

       ! *** parallel case
       IF (is_split) THEN
          coarse_slice_size=(coarse_bo(2,2)-coarse_bo(1,2)+1)*&
               (coarse_bo(2,3)-coarse_bo(1,3)+1)
          n_procs=coarse_coeffs_pw%pw_grid%para%group_size
          ALLOCATE(send_size(0:n_procs-1),send_offset(0:n_procs-1),&
               sent_size(0:n_procs-1),rcv_size(0:n_procs-1), &
               rcv_offset(0:n_procs-1), pp_lb(0:n_procs-1),&
               pp_ub(0:n_procs-1),real_rcv_size(0:n_procs-1),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

          ! ** send size count

          pos_of_x => coarse_coeffs_pw%pw_grid%para%pos_of_x
          send_size=0
          DO x=coarse_bo(1,1),coarse_bo(2,1)
             p=pos_of_x(coarse_gbo(1,1)+MODULO(x-coarse_gbo(1,1),s(1)))
             send_size(p)=send_size(p)+coarse_slice_size
          END DO
!tmpFM          WRITE(unit_nr,*) "send_size",send_size,"norm",send_size/coarse_slice_size

          ! ** rcv size count

          pos_of_x => fine_values_pw%pw_grid%para%pos_of_x
          p_old=pos_of_x(fine_gbo(1,1))
          pp_lb(p_old)=fine_gbo(1,1)
          DO x=fine_gbo(1,1),fine_gbo(2,1)
             p=pos_of_x(x)
             IF (p/=p_old) THEN
                pp_ub(p_old)=x-1
                pp_lb(p)=x
                p_old=p
             END IF
          END DO
          pp_ub(p_old)=fine_gbo(2,2)
          
          DO ip=0,n_procs-1
             pp_lb(ip)=FLOOR(REAL(pp_lb(ip)-f_shift(1),dp)/2._dp)-1
             pp_ub(ip)=FLOOR(REAL(pp_ub(ip)+1-f_shift(1),dp)/2._dp)+1
             IF (.NOT.is_split.OR..NOT.pbc) THEN
                pp_lb(ip)=MAX(pp_lb(ip),coarse_gbo(1,1))
                pp_ub(ip)=MIN(pp_ub(ip),coarse_gbo(2,1))
             END IF
          END DO

          rcv_size=0
          DO ip=0,n_procs-1
             DO x=pp_lb(ip),coarse_gbo(1,1)-1
                x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
                IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
                   rcv_size(ip)=rcv_size(ip)+coarse_slice_size
                END IF
             END DO
             rcv_size(ip)=rcv_size(ip)+coarse_slice_size*&
               MAX(0,&
               MIN(pp_ub(ip),my_coarse_bo(2,1))-MAX(pp_lb(ip),my_coarse_bo(1,1))+1)
             DO x=coarse_gbo(2,1)+1,pp_ub(ip)
                x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
                IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
                   rcv_size(ip)=rcv_size(ip)+coarse_slice_size
                END IF
             END DO
          END DO
!-FM
!-FM          pos_of_x => fine_values_pw%pw_grid%para%pos_of_x
!-FM          sf=fine_gbo(2,1)-fine_gbo(1,1)+1
!-FM          WRITE(unit_nr,*) "mepos=",logger%para_env%mepos,"pos_of_x",pos_of_x
!-FM          WRITE(unit_nr,*) "fine_bo",fine_values_pw%pw_grid%bounds
!-FM          WRITE(unit_nr,*) "fine_bo_local",fine_values_pw%pw_grid%bounds_local
!-FM          fi_lb=2*my_coarse_bo(1,1)-3+f_shift(1)
!-FM          fi_ub=2*my_coarse_bo(2,1)+3+f_shift(1)
!-FM          IF (.not.pbc) THEN
!-FM             fi_lb=MAX(fi_lb,fine_gbo(1,1))
!-FM             fi_ub=MIN(fi_ub,fine_gbo(2,1))
!-FM          ELSE
!-FM             fi_ub=MIN(fi_ub,fi_lb+sf-1)
!-FM          END IF
!-FM          p_old=pos_of_x(fine_gbo(1,1)+MODULO(fi_lb-fine_gbo(1,1),sf))
!-FM          p_lb=FLOOR((fi_lb-2-f_shift(1))/2.)
!-FM          rcv_size=0
!-FM          DO x=fi_lb,fi_ub
!-FM             p=pos_of_x(fine_gbo(1,1)+MODULO(x-fine_gbo(1,1),sf))
!-FM             WRITE(unit_nr,*) "x=",x,"p=",p
!-FM             IF (p/=p_old) THEN
!-FM                p_ub=FLOOR((x-1+3-f_shift(1))/2.)
!-FM
!-FM                rcv_size(p_old)=rcv_size(p_old)+(MIN(p_ub,my_coarse_bo(2,1)) &
!-FM                     -MAX(p_lb,my_coarse_bo(1,1))+1)*coarse_slice_size
!-FM
!-FM                IF (pbc) THEN
!-FM                   DO xx=p_lb,coarse_gbo(1,1)-1
!-FM                      x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
!-FM                      IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
!-FM                         rcv_size(p_old)=rcv_size(p_old)+coarse_slice_size
!-FM                      END IF
!-FM                   END DO
!-FM                   DO xx=coarse_gbo(2,1)+1,p_ub
!-FM                      x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
!-FM                      IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
!-FM                         rcv_size(p_old)=rcv_size(p_old)+coarse_slice_size
!-FM                      END IF
!-FM                   END DO
!-FM                END IF
!-FM
!-FM                p_old=p
!-FM                p_lb=FLOOR((x-2-f_shift(1))/2.)
!-FM             END IF
!-FM          END DO
!-FM          p_ub=FLOOR((fi_ub+3-f_shift(1))/2.)
!-FM
!-FM          rcv_size(p_old)=rcv_size(p_old)+(MIN(p_ub,my_coarse_bo(2,1)) &
!-FM               -MAX(p_lb,my_coarse_bo(1,1))+1)*coarse_slice_size
!-FM
!-FM          IF (pbc) THEN
!-FM             DO xx=p_lb,coarse_gbo(1,1)-1
!-FM                x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
!-FM                IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
!-FM                   rcv_size(p_old)=rcv_size(p_old)+coarse_slice_size
!-FM                END IF
!-FM             END DO
!-FM             DO xx=coarse_gbo(2,1)+1,p_ub
!-FM                x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
!-FM                IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
!-FM                   rcv_size(p_old)=rcv_size(p_old)+coarse_slice_size
!-FM                END IF
!-FM             END DO
!-FM          END IF

!tmpFM          WRITE(unit_nr,*) "rcv_size",rcv_size,"norm",rcv_size/coarse_slice_size

          ! ** offsets & alloc send-rcv

          send_tot_size=0
          DO ip=0,n_procs-1
             send_offset(ip)=send_tot_size
             send_tot_size=send_tot_size+send_size(ip)
          END DO
          call cp_assert(send_tot_size==&
               (coarse_bo(2,1)-coarse_bo(1,1)+1)*coarse_slice_size,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "Error calculating send_tot_size "//&
               CPSourceFileRef,&
               error,failure)
          ALLOCATE(send_buf(0:send_tot_size-1),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

          rcv_tot_size=0
          DO ip=0,n_procs-1
             rcv_offset(ip)=rcv_tot_size
             rcv_tot_size=rcv_tot_size+rcv_size(ip)
          END DO
          ALLOCATE(rcv_buf(0:rcv_tot_size-1),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

          ! ** fill send buffer

          pos_of_x => coarse_coeffs_pw%pw_grid%para%pos_of_x
          p_old=pos_of_x(coarse_gbo(1,1) &
               +MODULO(coarse_bo(1,1)-coarse_gbo(1,1),s(1)))
          sent_size=send_offset
          ss=coarse_bo(2,1)-coarse_bo(1,1)+1
          DO x=coarse_bo(1,1),coarse_bo(2,1)
             p=pos_of_x(coarse_gbo(1,1)+MODULO(x-coarse_gbo(1,1),s(1)))
             CALL dcopy(coarse_slice_size,&
                  coarse_coeffs(x,coarse_bo(1,2),&
                  coarse_bo(1,3)),ss,send_buf(sent_size(p)),1)
             sent_size(p)=sent_size(p)+coarse_slice_size
          END DO

          CALL cp_assert(ALL(sent_size(0:n_procs-2)==send_offset(1:n_procs-1)),&
               cp_failure_level,cp_assertion_failed,routineP,&
               "error 1 filling send buffer "//&
               CPSourceFileRef,&
               error,failure)
          CALL cp_assert(sent_size(n_procs-1)==send_tot_size,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "error 2 filling send buffer "//&
               CPSourceFileRef,&
               error,failure)

!tmpFM          WRITE(unit_nr,*)"sent_size",sent_size,"send_offset",send_offset
          CPPostcondition(ALL(sent_size(:n_procs-2)==send_offset(1:)),cp_failure_level,routineP,error,failure)
!tmpFM          WRITE(unit_nr,*)"ssize",sent_size(n_procs-1),"send_tot_size",send_tot_size
          CPPostcondition(sent_size(n_procs-1)==send_tot_size,cp_failure_level,routineP,error,failure)
!tmpFM          CALL m_flush(unit_nr)
          ! test send/rcv sizes
          CALL mp_alltoall(send_size,real_rcv_size,1,coarse_coeffs_pw%pw_grid%para%group)
          CPAssert(ALL(real_rcv_size==rcv_size),cp_failure_level,routineP,error,failure)
          ! all2all
          call mp_alltoall( sb=send_buf, scount=send_size, sdispl=send_offset,&
               rb=rcv_buf, rcount=rcv_size, rdispl=rcv_offset, &
               group=coarse_coeffs_pw%pw_grid%para%group )

          ! ** sum & reorder rcv buffer

          sent_size=rcv_offset
          DO ip=0,n_procs-1

             DO x=pp_lb(ip),coarse_gbo(1,1)-1
                x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
                IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
                   ii=sent_size(ip)
                   DO k=coarse_bo(1,3),coarse_bo(2,3)
                      DO j=coarse_bo(1,2),coarse_bo(2,2)
                         coarse_coeffs_pw%cr3d(x_att,j,k)=coarse_coeffs_pw%cr3d(x_att,j,k)+rcv_buf(ii)
                         ii=ii+1
                      END DO
                   END DO
                   sent_size(ip)=ii
                END IF
             END DO

             ii=sent_size(ip)
             DO x_att=MAX(pp_lb(ip),my_coarse_bo(1,1)),MIN(pp_ub(ip),my_coarse_bo(2,1))
                DO k=coarse_bo(1,3),coarse_bo(2,3)
                   DO j=coarse_bo(1,2),coarse_bo(2,2)
                      coarse_coeffs_pw%cr3d(x_att,j,k)=coarse_coeffs_pw%cr3d(x_att,j,k)+rcv_buf(ii)
                      ii=ii+1
                   END DO
                END DO
             END DO
             sent_size(ip)=ii

             DO x=coarse_gbo(2,1)+1,pp_ub(ip)
                x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
                IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
                   ii=sent_size(ip)
                   DO k=coarse_bo(1,3),coarse_bo(2,3)
                      DO j=coarse_bo(1,2),coarse_bo(2,2)
                         coarse_coeffs_pw%cr3d(x_att,j,k)=coarse_coeffs_pw%cr3d(x_att,j,k)+rcv_buf(ii)
                         ii=ii+1
                      END DO
                   END DO
                   sent_size(ip)=ii
                END IF
             END DO

          END DO

          CALL cp_assert(ALL(sent_size(0:n_procs-2)==rcv_offset(1:n_procs-1)),&
               cp_failure_level,cp_assertion_failed,routineP,&
               "error 1 handling the rcv buffer "//&
               CPSourceFileRef,&
               error,failure)
          CALL cp_assert(sent_size(n_procs-1)==rcv_tot_size,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "error 2 handling the rcv buffer "//&
               CPSourceFileRef,&
               error,failure)
          
!-FM          pos_of_x => fine_values_pw%pw_grid%para%pos_of_x
!-FM          sf=fine_gbo(2,1)-fine_gbo(1,1)+1
!-FM          fi_lb=2*my_coarse_bo(1,1)-3+f_shift(1)
!-FM          fi_ub=2*my_coarse_bo(2,1)+3+f_shift(1)
!-FM          IF (.not.pbc) THEN
!-FM             fi_lb=MAX(fi_lb,fine_gbo(1,1))
!-FM             fi_ub=MIN(fi_ub,fine_gbo(2,1))
!-FM          ELSE
!-FM             fi_ub=MIN(fi_ub,fi_lb+sf-1)
!-FM          END IF
!-FM          p_old=pos_of_x(fine_gbo(1,1)+MODULO(fi_lb-fine_gbo(1,1),sf))
!-FM          p_lb=FLOOR((fi_lb-2-f_shift(1))/2.)
!-FM          sent_size=rcv_offset
!-FM          DO x=fi_lb,fi_ub
!-FM             p=pos_of_x(fine_gbo(1,1)+MODULO(x-fine_gbo(1,1),sf))
!-FM             WRITE(unit_nr,*) "x=",x,"p=",p
!-FM             IF (p/=p_old) THEN
!-FM                p_ub=FLOOR((x-1+3-f_shift(1))/2.)
!-FM
!-FM                IF (pbc) THEN
!-FM                   DO xx=p_lb,coarse_gbo(1,1)-1
!-FM                      x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
!-FM                      IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
!-FM                         ii=sent_size(p_old)
!-FM                         DO k=coarse_bo(1,3),coarse_bo(2,3)
!-FM                            DO j=coarse_bo(1,2),coarse_bo(2,2)
!-FM                               coarse_coeffs_pw%cr3d(x_att,j,k)=coarse_coeffs_pw%cr3d(x_att,j,k)+rcv_buf(ii)
!-FM                               ii=ii+1
!-FM                            END DO
!-FM                         END DO
!-FM                         sent_size(p_old)=ii
!-FM                      END IF
!-FM                   END DO
!-FM                END IF
!-FM
!-FM                ii=sent_size(p_old)
!-FM                DO k=coarse_bo(1,3),coarse_bo(2,3)
!-FM                   DO j=coarse_bo(1,2),coarse_bo(2,2)
!-FM                      DO i=MAX(p_lb,my_coarse_bo(1,1)),MIN(p_ub,my_coarse_bo(2,1))
!-FM                         coarse_coeffs_pw%cr3d(i,j,k)=coarse_coeffs_pw%cr3d(x_att,j,k)+rcv_buf(ii)
!-FM                         ii=ii+1
!-FM                      END DO
!-FM                   END DO
!-FM                END DO
!-FM                sent_size(p_old)=ii
!-FM
!-FM                IF (pbc) THEN
!-FM                   DO xx=coarse_gbo(2,1)+1,p_ub
!-FM                      x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
!-FM                      IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
!-FM                         ii=sent_size(p_old)
!-FM                         DO k=coarse_bo(1,3),coarse_bo(2,3)
!-FM                            DO j=coarse_bo(1,2),coarse_bo(2,2)
!-FM                               coarse_coeffs_pw%cr3d(x_att,j,k)=coarse_coeffs_pw%cr3d(x_att,j,k)+rcv_buf(ii)
!-FM                               ii=ii+1
!-FM                            END DO
!-FM                         END DO
!-FM                         sent_size(p_old)=ii
!-FM                      END IF
!-FM                   END DO
!-FM                END IF
!-FM
!-FM                p_old=p
!-FM                p_lb=FLOOR((x-2-f_shift(1))/2.)
!-FM             END IF
!-FM          END DO
!-FM          p_ub=FLOOR((fi_ub+3-f_shift(1))/2.)
!-FM
!-FM          IF (pbc) THEN
!-FM             DO xx=p_lb,coarse_gbo(1,1)-1
!-FM                x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
!-FM                IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
!-FM                   ii=sent_size(p_old)
!-FM                   DO k=coarse_bo(1,3),coarse_bo(2,3)
!-FM                      DO j=coarse_bo(1,2),coarse_bo(2,2)
!-FM                         coarse_coeffs_pw%cr3d(x_att,j,k)=coarse_coeffs_pw%cr3d(x_att,j,k)+rcv_buf(ii)
!-FM                         ii=ii+1
!-FM                      END DO
!-FM                   END DO
!-FM                   sent_size(p_old)=ii
!-FM                END IF
!-FM             END DO
!-FM          END IF
!-FM
!-FM          ii=sent_size(p_old)
!-FM          DO k=coarse_bo(1,3),coarse_bo(2,3)
!-FM             DO j=coarse_bo(1,2),coarse_bo(2,2)
!-FM                DO i=MAX(p_lb,my_coarse_bo(1,1)),MIN(p_ub,my_coarse_bo(2,1))
!-FM                   coarse_coeffs_pw%cr3d(i,j,k)=coarse_coeffs_pw%cr3d(x_att,j,k)+rcv_buf(ii)
!-FM                   ii=ii+1
!-FM                END DO
!-FM             END DO
!-FM          END DO
!-FM          sent_size(p_old)=ii
!-FM
!-FM          IF (pbc) THEN
!-FM             DO xx=coarse_gbo(2,1)+1,p_ub
!-FM                x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
!-FM                IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
!-FM                   ii=sent_size(p_old)
!-FM                   DO k=coarse_bo(1,3),coarse_bo(2,3)
!-FM                      DO j=coarse_bo(1,2),coarse_bo(2,2)
!-FM                         coarse_coeffs_pw%cr3d(x_att,j,k)=coarse_coeffs_pw%cr3d(x_att,j,k)+rcv_buf(ii)
!-FM                         ii=ii+1
!-FM                      END DO
!-FM                   END DO
!-FM                   sent_size(p_old)=ii
!-FM                END IF
!-FM             END DO
!-FM          END IF

          ! dealloc
          DEALLOCATE(send_size,send_offset, rcv_size, rcv_offset, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(send_buf,rcv_buf,real_rcv_size, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(pp_ub,pp_lb,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
!tmpFM       WRITE(unit_nr,*) 'is_split',is_split, ' has_i_lbound',has_i_lbound,&
!tmpFM            ' has_i_ubound',has_i_ubound
!tmpFM       WRITE(unit_nr,*) "END add_fine2coarse"
    END IF
    call mp_sync ( coarse_coeffs_pw%pw_grid%para%group )

  END SUBROUTINE add_fine2coarse
  !***************************************************************************

  !!****f* pw_spline_utils/pw_spline_precond_create [1.0] *
  !!
  !!   NAME
  !!     creates a preconditioner for pw splines
  !!
  !!   FUNCTION
  !!     -
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - preconditioner: the preconditioner to create
  !!     - precond_kind: the kind of preconditioner to use
  !!     - pool: a pool with grids of the same type as the elements to 
  !!       precondition
  !!     - pbc: if periodic boundary conditions should be applied
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE pw_spline_precond_create(preconditioner,precond_kind,&
       pool,pbc,error)
    TYPE(pw_spline_precond_type), pointer :: preconditioner
    INTEGER, intent(in) :: precond_kind
    TYPE(pw_pool_type), pointer :: pool
    LOGICAL, intent(in) :: pbc
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='pw_spline_precond_create',&
         routineP=moduleN//':'//routineN
    integer :: stat

    failure=.false.

    ALLOCATE(preconditioner,stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       last_precond_id=last_precond_id+1
       preconditioner%id_nr=last_precond_id
       preconditioner%ref_count=1
       preconditioner%kind=no_precond
       preconditioner%pool => pool
       preconditioner%pbc=pbc
       CALL pw_pool_retain(pool,error=error)
       CALL pw_spline_precond_set_kind(preconditioner,precond_kind,error=error)
    END IF
  END SUBROUTINE pw_spline_precond_create
  !***************************************************************************

  !!****f* pw_spline_utils/pw_spline_precond_set_kind *
  !!
  !!   NAME
  !!     pw_spline_precond_set_kind
  !!
  !!   FUNCTION
  !!     switches the types of precoditioner to use
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - precoditioner: the preconditioner to be changed
  !!     - precond_kind: the new kind of preconditioner to use
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE pw_spline_precond_set_kind(preconditioner,precond_kind,pbc,error)
    TYPE(pw_spline_precond_type), pointer :: preconditioner
    INTEGER, intent(in) :: precond_kind
    LOGICAL, INTENT(in), OPTIONAL :: pbc
    type(cp_error_type), optional, intent(inout) :: error

    LOGICAL :: failure,do_3d_coeff
    INTEGER :: i,j,k
    real(kind=dp) :: s
    character(len=*), parameter :: routineN='pw_spline_precond_set_kind',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(ASSOCIATED(preconditioner),cp_failure_level,routineP,error,failure)
    CPPrecondition(preconditioner%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       do_3d_coeff=.false.
       preconditioner%kind=precond_kind
       if (present(pbc)) preconditioner%pbc=pbc
       SELECT CASE(precond_kind)
       CASE(no_precond)
       CASE (precond_spl3_aint)
          preconditioner%coeffs_1d=spl3_1d_coeffs0
          preconditioner%sharpen=.true.
          preconditioner%normalize=.true.
          do_3d_coeff=.true.
       CASE (precond_spl3_1)
          preconditioner%coeffs_1d(1)=0.5_dp/3._dp**(1._dp/3._dp)
          preconditioner%coeffs_1d(2)=4._dp/3._dp**(1._dp/3._dp)
          preconditioner%coeffs_1d(3)=0.5_dp/3._dp**(1._dp/3._dp)
          preconditioner%sharpen=.true.
          preconditioner%normalize=.false.
          do_3d_coeff=.true.
       CASE default
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       IF (do_3d_coeff) THEN
          s=1._dp
          IF (preconditioner%sharpen) s=-1._dp
          preconditioner%coeffs(1)=&
               s*preconditioner%coeffs_1d(2)*&
               preconditioner%coeffs_1d(2)*&
               preconditioner%coeffs_1d(2)
          preconditioner%coeffs(2)=&
               s*preconditioner%coeffs_1d(1)*&
               preconditioner%coeffs_1d(2)*&
               preconditioner%coeffs_1d(2)
          preconditioner%coeffs(3)=&
               s*preconditioner%coeffs_1d(1)*&
               preconditioner%coeffs_1d(1)*&
               preconditioner%coeffs_1d(2)
          preconditioner%coeffs(4)=&
               s*preconditioner%coeffs_1d(1)*&
               preconditioner%coeffs_1d(1)*&
               preconditioner%coeffs_1d(1)
          IF (preconditioner%sharpen) THEN
             IF (preconditioner%normalize) THEN
                preconditioner%coeffs(1)=2._dp+&
                     preconditioner%coeffs(1)
             ELSE
                preconditioner%coeffs(1)=-preconditioner%coeffs(1)
             END IF
          END IF
       END IF
    END IF
  END SUBROUTINE pw_spline_precond_set_kind
  !***************************************************************************

  !!****f* pw_spline_utils/pw_spline_precond_retain [1.0] *
  !!
  !!   NAME
  !!     pw_spline_precond_retain
  !!
  !!   FUNCTION
  !!     retains the preconditioner
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - preconditioner: the preconditioner to retain
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE pw_spline_precond_retain(preconditioner,error)
    TYPE(pw_spline_precond_type), pointer :: preconditioner
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='pw_spline_precond_retain',&
         routineP=moduleN//':'//routineN

    failure=.false.

    CPPrecondition(ASSOCIATED(preconditioner),cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       CPPreconditionNoFail(preconditioner%ref_count>1,cp_failure_level,routineP,error)
       preconditioner%ref_count=preconditioner%ref_count+1
    end if
  END SUBROUTINE pw_spline_precond_retain
  !***************************************************************************

  !!****f* pw_spline_utils/pw_spline_precond_release [1.0] *
  !!
  !!   NAME
  !!     pw_spline_precond_release
  !!
  !!   FUNCTION
  !!     releases the preconditioner
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - preconditioner: the preconditioner to release
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE pw_spline_precond_release(preconditioner,error)
    TYPE(pw_spline_precond_type), pointer :: preconditioner
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='pw_spline_precond_release',&
         routineP=moduleN//':'//routineN
    integer :: stat

    failure=.false.

    if (associated(preconditioner)) then
       CPPreconditionNoFail(preconditioner%ref_count>0,cp_failure_level,routineP,error)
       preconditioner%ref_count=preconditioner%ref_count-1
       IF (preconditioner%ref_count==0) THEN
          CALL pw_pool_release(preconditioner%pool,error=error)
          DEALLOCATE(preconditioner,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    end if
  END SUBROUTINE pw_spline_precond_release
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
  !!     - preconditioner: the preconditioner to apply
  !!     - v_in: the grid on which the preconditioner should be applied
  !!     - v_out: place to store the preconditioner applied on v_out
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
    TYPE(pw_type), pointer :: pw_tmp

    failure=.false.

    CPPrecondition(ASSOCIATED(preconditioner),cp_failure_level,routineP,error,failure)
    CPPrecondition(preconditioner%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       SELECT CASE(preconditioner%kind)
       CASE(no_precond)
          CALL pw_copy(in_v,out_v)
       CASE (precond_spl3_aint,precond_spl3_1)
          CALL pw_zero(out_v)
          IF (preconditioner%pbc) THEN
             CALL pw_nn_smear_r(pw_in=in_v,pw_out=out_v,&
                  coeffs=preconditioner%coeffs,error=error)
          ELSE
             CALL pw_nn_compose_r_no_pbc(weights_1d=preconditioner%coeffs_1d,&
                  pw_in=in_v, pw_out=out_v, sharpen=preconditioner%sharpen,&
                  normalize=preconditioner%normalize, error=error)
          END IF
       CASE default
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
    END IF
  END SUBROUTINE pw_spline_do_precond
  !***************************************************************************

  !!****f* pw_spline_utils/find_coeffs [1.0] *
  !!
  !!   NAME
  !!     find_coeffs
  !!
  !!   FUNCTION
  !!     solves iteratively (CG) a systmes of linear equations
  !!          linOp(coeffs)=values
  !!     (for example those needed to find the coefficents of a spline)
  !!     Returns true if the it succeded to acheive the requested accuracy
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - values: the right hand side of the system
  !!     - coeff: will contain the solution of the system (and on entry
  !!       it contains the starting point)
  !!     - linOp: the linear operator to be inverted
  !!     - preconditioner: the preconditioner to apply
  !!     - pool: a pool of grids (for the temporary objects)
  !!     - eps_r: the requested precision on the residual
  !!     - eps_x: the requested precision on the solution
  !!     - max_iter: maximum number of iteration allowed
  !!     - print_section: information about when to print information
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  FUNCTION find_coeffs(values,coeffs,linOp,preconditioner, pool, &
       eps_r,eps_x,max_iter,print_section,error) RESULT(res)
    TYPE(pw_type), POINTER :: values,coeffs
    type(cp_error_type), optional, intent(inout) :: error
    TYPE(pw_spline_precond_type), pointer :: preconditioner
    TYPE(pw_pool_type), pointer :: pool
    INTEGER, intent(in) :: max_iter
    REAL(kind=dp), INTENT(in) :: eps_r, eps_x
    TYPE(section_vals_type), POINTER, optional :: print_section
    logical :: res
    INTERFACE
       SUBROUTINE linOp(pw_in,pw_out,error)
         USE pw_types, only: pw_type
         USE cp_error_handling, only: cp_error_type
         TYPE(pw_type), POINTER :: pw_in,pw_out
         TYPE(cp_error_type), INTENT(inout), optional :: error
       END SUBROUTINE linOp
    END INTERFACE

    LOGICAL :: failure, last
    character(len=*), parameter :: routineN='find_coeffs',&
         routineP=moduleN//':'//routineN
    TYPE(pw_type), POINTER :: r,z,p,Ap
    REAL(kind=dp) :: alpha,beta,r_z,r_z_new,eps_x_att,eps_r_att
    INTEGER :: i,j,k,iter,iiter,n_rlevel
    INTEGER, DIMENSION(2,3) :: bo
    TYPE(section_vals_type), pointer :: print_key
    TYPE(cp_logger_type), pointer :: logger
    failure=.false.
    last=.false.

    res=.false.
    IF (.NOT. failure) THEN
       logger => cp_error_get_logger(error)
       CALL pw_pool_create_pw(pool,r,use_data=REALDATA3D,in_space=REALSPACE,&
            error=error)
       CALL pw_pool_create_pw(pool,z,use_data=REALDATA3D,in_space=REALSPACE,&
            error=error)
       CALL pw_pool_create_pw(pool,p,use_data=REALDATA3D,in_space=REALSPACE,&
            error=error)
       CALL pw_pool_create_pw(pool,Ap,use_data=REALDATA3D,in_space=REALSPACE,&
            error=error)

       nullify(print_key)
       IF (PRESENT(print_section)) THEN
          IF (ASSOCIATED(print_section)) THEN
             print_key => section_vals_get_subs_vals(print_section,"conv_info",&
                  error=error)
          END IF
       END IF

       CALL cp_add_iter_level(logger%iter_info,level_name="findcoef",n_rlevel_new=n_rlevel,error=error)
       ext_do:DO iiter=1,max_iter,10
          CALL pw_zero(r)
          CALL linOp(pw_in=coeffs,pw_out=r,error=error)
          r%cr3d=-r%cr3d
          CALL pw_sumup(values,r)
          CALL pw_spline_do_precond(preconditioner,in_v=r,out_v=z,error=error)
          CALL pw_copy(z,p)
          r_z=pw_integral_ab(r,z)

          DO iter=iiter,MIN(iiter+9,max_iter)
             eps_r_att=sqrt(pw_integral_aa(r,SQUARE))
             IF (eps_r_att==0._dp) THEN
                eps_x_att=0._dp
                last=.true.
             ELSE
                CALL pw_zero(Ap)
                CALL linOp(pw_in=p,pw_out=Ap,error=error)
                alpha=r_z/pw_integral_ab(Ap,p)

                CALL pw_sumup(p,coeffs,alpha=alpha)

                eps_x_att=alpha*SQRT(pw_integral_aa(p,SQUARE)) ! try to spare if unneded?
                IF (eps_r_att<eps_r .AND. eps_x_att<eps_x) last=.TRUE.
             END IF
             CALL cp_iterate(logger%iter_info,last=last,error=error)
             IF (cp_should_output(logger%iter_info,print_key,error=error).AND.&
                  logger%para_env%mepos==logger%para_env%source) THEN
                WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.false.),&
                     fmt="(i6,' eps_x=',e12.4,' eps_r=',e12.4)")&
                     iter,eps_x_att,eps_r_att
             END IF
             IF (last) THEN
                res=.true.
                exit ext_do
             END IF

             CALL pw_sumup(Ap,r,alpha=-alpha)

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
       END DO ext_do
       CALL cp_rm_iter_level(logger%iter_info,n_rlevel_att=n_rlevel,error=error)

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
  SUBROUTINE pw_project_s3(pw_fine_in,pw_coarse_out,coarse_pool,param_section,&
       error)
    TYPE(pw_type), pointer :: pw_fine_in
    TYPE(pw_type), pointer :: pw_coarse_out
    TYPE(pw_pool_type), POINTER :: coarse_pool
    TYPE(section_vals_type), pointer :: param_section
    type(cp_error_type), optional, intent(inout) :: error

    LOGICAL :: failure,safe_computation
    character(len=*), parameter :: routineN='pw_project_s3',&
         routineP=moduleN//':'//routineN
    TYPE(pw_type), POINTER :: values,coeffs
    INTEGER :: i,j,k, max_iter, aint_precond, precond_kind, unit_nr
    INTEGER, DIMENSION(2,3) :: bo
    TYPE(pw_spline_precond_type), pointer :: precond
    logical :: success
    REAL(kind=dp) :: eps_r,eps_x
    INTEGER, save :: ifile=0
    CHARACTER(80) :: filename

    failure=.false.
    ifile=ifile+1

    IF (.NOT. failure) THEN
       CALL section_vals_val_get(param_section,"safe_computation", &
            l_val=safe_computation, error=error)
       CALL section_vals_val_get(param_section,"aint_precond", &
            i_val=aint_precond, error=error)
       CALL section_vals_val_get(param_section,"precond", &
            i_val=precond_kind, error=error)
       CALL section_vals_val_get(param_section,"max_iter", &
            i_val=max_iter, error=error)
       CALL section_vals_val_get(param_section,"eps_r", &
            r_val=eps_r, error=error)
       CALL section_vals_val_get(param_section,"eps_x", &
            r_val=eps_x, error=error)

       aint_precond=no_precond
       precond_kind=no_precond
       bo=pw_coarse_out%pw_grid%bounds_local
       NULLIFY(values,coeffs)
       CALL pw_pool_create_pw(coarse_pool,values, use_data=REALDATA3D,&
            in_space=REALSPACE,error=error)

       call pw_zero(values)
       call add_fine2coarse(fine_values_pw=pw_fine_in,&
            coarse_coeffs_pw=values,&
            weights_1d=spl3_1d_transf_coeffs/2._dp, w_border0=0.5_dp,&
            w_border1=spl3_1d_transf_border1/2._dp,pbc=.FALSE.,&
            safe_computation=safe_computation,error=error)
!FM     WRITE(filename,"(a,i6,a)") "p1_fine-",ifile,".cube"
!FM     CALL compress(filename,full=.true.)
!FM     CALL open_file(unit_number=unit_nr,file_name=filename,file_action="write")
!FM     call rs_pw_to_cube(pw=pw_fine_in, iunit=unit_nr, ionode=.true., title="values")
!FM     call close_file(unit_nr)
!FM     WRITE(filename,"(a,i6,a)") "p2_coarse-",ifile,".cube"
!FM     CALL compress(filename,full=.true.)
!FM     CALL open_file(unit_number=unit_nr,file_name=filename,file_action="write")
!FM     call rs_pw_to_cube(pw=values, iunit=unit_nr, ionode=.true., title="values")
!FM     call close_file(unit_nr)

       CALL pw_pool_create_pw(coarse_pool,coeffs, use_data=REALDATA3D,&
            in_space=REALSPACE,error=error)
       nullify(precond)
       call pw_spline_precond_create(precond,precond_kind=aint_precond,&
            pool=coarse_pool,pbc=.false.,error=error)
       CALL pw_spline_do_precond(precond,values,coeffs,error=error)
       CALL pw_spline_precond_set_kind(precond,precond_kind,error=error)
       success=find_coeffs(values=values,coeffs=coeffs,&
            linOp=spl3_nopbct,preconditioner=precond, pool=coarse_pool, &
            eps_r=eps_r,eps_x=eps_x,print_section=param_section,&
            max_iter=max_iter,error=error)
!FM     WRITE(filename,"(a,i6,a)") "p3_inv-",ifile,".cube"
!FM     CALL compress(filename,full=.true.)
!FM     CALL open_file(unit_number=unit_nr,file_name=filename,file_action="write")
!FM     call rs_pw_to_cube(pw=coeffs, iunit=unit_nr, ionode=.true., title="values")
!FM     call close_file(unit_nr)
!FM     CPPostconditionNoFail(success,cp_warning_level,routineP,error)
       CALL pw_spline_precond_release(precond,error=error)

       call pw_zero(pw_coarse_out)
       CALL pw_sumup(coeffs,pw_coarse_out)

       CALL pw_pool_give_back_pw(coarse_pool,values,error=error)
       CALL pw_pool_give_back_pw(coarse_pool,coeffs,error=error)
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
  SUBROUTINE pw_inject_s3(pw_coarse_in,pw_fine_out,coarse_pool,param_section,&
       error)
    TYPE(pw_type), pointer :: pw_coarse_in
    TYPE(pw_type), pointer :: pw_fine_out
    TYPE(pw_pool_type), POINTER :: coarse_pool
    TYPE(section_vals_type), pointer :: param_section
    type(cp_error_type), optional, intent(inout) :: error

    LOGICAL :: failure, success, safe_computation
    character(len=*), parameter :: routineN='pw_inject_s3',&
         routineP=moduleN//':'//routineN
    TYPE(pw_type), POINTER :: coeffs
    INTEGER :: i,j,k, handle,aint_precond,precond_kind,max_iter,unit_nr
    INTEGER, DIMENSION(2,3) :: bo
    REAL(kind=dp) :: eps_x, eps_r
    TYPE(pw_spline_precond_type), pointer :: precond
    TYPE(section_vals_type), pointer :: print_key
    TYPE(cp_logger_type), pointer :: logger
    INTEGER, save :: ifile=0
    character(len=80) :: filename

    failure=.false.

    ifile=ifile+1
    CALL timeset(routineN,"I","",handle)
    IF (.NOT. failure) THEN
       NULLIFY(coeffs)
       logger => cp_error_get_logger(error)
       CALL pw_pool_create_pw(coarse_pool,coeffs, use_data=REALDATA3D,&
            in_space=REALSPACE,error=error)
       bo=pw_coarse_in%pw_grid%bounds_local

       CALL section_vals_val_get(param_section,"safe_computation", &
            l_val=safe_computation, error=error)
       CALL section_vals_val_get(param_section,"aint_precond", &
            i_val=aint_precond, error=error)
       CALL section_vals_val_get(param_section,"precond", &
            i_val=precond_kind, error=error)
       CALL section_vals_val_get(param_section,"max_iter", &
            i_val=max_iter, error=error)
       CALL section_vals_val_get(param_section,"eps_r", &
            r_val=eps_r, error=error)
       CALL section_vals_val_get(param_section,"eps_x", &
            r_val=eps_x, error=error)

       nullify(precond)
       call pw_spline_precond_create(precond,precond_kind=aint_precond,&
            pool=coarse_pool,pbc=.FALSE.,error=error)
       CALL pw_spline_do_precond(precond,pw_coarse_in,coeffs,error=error)
       CALL pw_spline_precond_set_kind(precond,precond_kind,error=error)
       success=find_coeffs(values=pw_coarse_in,coeffs=coeffs,&
            linOp=spl3_nopbc,preconditioner=precond, pool=coarse_pool, &
            eps_r=eps_r,eps_x=eps_x,&
            max_iter=max_iter,print_section=param_section,error=error)
       CPPostconditionNoFail(success,cp_warning_level,routineP,error)
       CALL pw_spline_precond_release(precond,error=error)

       print_key => section_vals_get_subs_vals(param_section,"spl_coeffs",&
            error=error)
       IF (cp_should_output(logger%iter_info,print_key,error=error)) THEN
          unit_nr=cp_unitnr_for_output(logger=logger, print_key=print_key,&
               extension=".cube", local=.false., log_filename=.true., error=error)
          CALL rs_pw_to_cube(coeffs,unit_nr,&
               logger%para_env%source==logger%para_env%mepos,"spl_coeffs")
          CALL cp_finished_output(unit_nr,logger=logger, print_key=print_key,&
               error=error)
       END IF

       CALL add_coarse2fine(coarse_coeffs_pw=coeffs,&
            fine_values_pw=pw_fine_out,&
            weights_1d=spl3_1d_transf_coeffs,&
            w_border0=1._dp,&
            w_border1=spl3_1d_transf_border1,&
            pbc=.FALSE.,safe_computation=safe_computation,&
            error=error)
       CALL pw_pool_give_back_pw(coarse_pool,coeffs,error=error)

    END IF
    CALL timestop(0._dp,handle)
  END SUBROUTINE pw_inject_s3
  !***************************************************************************

  !!****f* pw_spline_utils/pw_nn_compose_r_no_pbc *
  !!
  !!   NAME
  !!     pw_compose_r_no_pbc
  !!
  !!   FUNCTION
  !!     adds to pw_out pw_in composed with the weights
  !!     pw_out%cr3d(i,j,k)=pw_out%cr3d(i,j,k)+sum(pw_in%cr3d(i+l,j+m,k+n)*
  !!        weights_1d(abs(l)+1)*weights_1d(abs(m)+1)*weights_1d(abs(n)+1),
  !!        l=-1..1,m=-1..1,n=-1..1)
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE pw_nn_compose_r_no_pbc(weights_1d,pw_in,pw_out,&
       sharpen,normalize,transpose,error)
    REAL(kind=dp), DIMENSION(-1:1)           :: weights_1d
    TYPE(pw_type), POINTER                   :: pw_in, pw_out
    LOGICAL, INTENT(in), OPTIONAL            :: sharpen,normalize,transpose
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                            :: error

    CHARACTER(len=*), PARAMETER :: routineN = "pw_nn_compose_r_work", &
         routineP = moduleN//":"//routineN

    INTEGER                                  :: i, j, jw, k, kw, myj, myk, &
         stat,n_els,first_index, last_index
    INTEGER, DIMENSION(2, 3)                 :: bo,gbo
    INTEGER, DIMENSION(3)                    :: s
    LOGICAL                                  :: failure, has_l_boundary, &
         has_u_boundary,is_split,my_sharpen,my_normalize, my_transpose
    REAL(kind=dp)                            :: in_val_f, in_val_l,w_k,w_j
    REAL(kind=dp), DIMENSION(-1:1)           :: w
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: l_boundary, tmp, u_boundary
    REAL(kind=dp), DIMENSION(:,:,:), POINTER :: in_val, out_val
    REAL(kind=dp), parameter :: w_border=0.0_dp

    bo=pw_in%pw_grid%bounds_local
    gbo=pw_in%pw_grid%bounds
    in_val => pw_in%cr3d
    out_val => pw_out%cr3d
    my_sharpen=.false.
    if (present(sharpen)) my_sharpen=sharpen
    my_normalize=.false.
    if (present(normalize)) my_normalize=normalize
    my_transpose=.FALSE.
    if (present(transpose)) my_transpose=transpose
    DO i=1,3
       s(i)=bo(2,i)-bo(1,i)+1
    END DO
    IF (ANY(s<1)) RETURN
    is_split= ANY(pw_in%pw_grid%bounds_local(:,1) /= &
         pw_in%pw_grid%bounds(:,1))
    has_l_boundary=(gbo(1,1)==bo(1,1))
    has_u_boundary= (gbo(2,1)==bo(2,1))
    IF (is_split) THEN
       ALLOCATE(l_boundary(bo(1,2):bo(2,2),bo(1,3):bo(2,3)),&
            u_boundary(bo(1,2):bo(2,2),bo(1,3):bo(2,3)),&
            tmp(bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
       DEALLOCATE(tmp,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

    n_els=s(1)
    IF (has_l_boundary) THEN
       n_els=n_els-1
       first_index=bo(1,1)+1
    ELSE
       first_index=bo(1,1)
    END IF
    IF (has_u_boundary) THEN
       n_els=n_els-1
       last_index=bo(2,1)-1
    ELSE
       last_index=bo(2,1)
    END IF
    !#omp parallel do default(none) private(k,kw,myk,j,jw,myj,in_val_f,&
    !#omp     in_val_l) shared(zderiv,yderiv,bo,in_val,out_val,s,l_boundary,&
    !#omp     u_boundary,weights,is_split)
    DO k=bo(1,3),bo(2,3)
       DO kw=-1,1
          myk=k+kw
          IF (my_transpose) THEN
             IF (k>=gbo(2,3)-1.OR.k<=gbo(1,3)+1) THEN
                IF (k==gbo(2,3).OR.k==gbo(1,3)) THEN
                   IF (myk<gbo(2,3).AND.myk>gbo(1,3)) THEN
                      w_k=weights_1d(kw)
                   ELSE IF (kw==0) THEN
                      w_k=1._dp
                   ELSE
                      CYCLE
                   END IF
                ELSE
                   IF (myk==gbo(2,3).OR.myk==gbo(1,3)) CYCLE
                   w_k=weights_1d(kw)
                END IF
             ELSE
                w_k=weights_1d(kw)
             END IF
          ELSE
             IF (k==gbo(2,3).OR.k==gbo(1,3)) THEN
                IF (kw/=0) CYCLE
                w_k=1._dp
             ELSE
                w_k=weights_1d(kw)
             END IF
          END IF
          DO j=bo(1,2),bo(2,2)
             DO jw=-1,1
                myj=j+jw
                IF (my_transpose) THEN
                   IF (j>=gbo(2,2)-1.OR.j<=gbo(1,2)+1) THEN
                      IF (j==gbo(2,2).OR.j==gbo(1,2)) THEN
                         IF (myj<gbo(2,2).AND.myj>gbo(1,2)) THEN
                            w_j=weights_1d(jw)*w_k
                         ELSE IF (jw==0) THEN
                            w_j=w_k
                         ELSE
                            CYCLE
                         END IF
                      ELSE
                         IF (myj==gbo(2,2).OR.myj==gbo(1,2)) CYCLE
                         w_j=w_k*weights_1d(jw)
                      END IF
                   ELSE
                      w_j=w_k*weights_1d(jw)
                   END IF
                ELSE
                   IF (j==gbo(2,2).OR.j==gbo(1,2)) THEN
                      IF (jw/=0) CYCLE
                      w_j=w_k
                   ELSE
                      w_j=w_k*weights_1d(jw)
                   END IF
                END IF

                IF (has_l_boundary) THEN
                   IF (my_transpose) THEN
                      IF (s(1)==1) THEN
                         CPAssert(.not.has_u_boundary,cp_failure_level,routineP,error,failure)
                         out_val(bo(1,1),j,k)=out_val(bo(1,1),j,k)+&
                              w_j*(in_val(bo(1,1),myj,myk)+&
                              u_boundary(myj,myk)*weights_1d(1))
                         in_val_f=0.0_dp
                      ELSE
                         out_val(bo(1,1),j,k)=out_val(bo(1,1),j,k)+&
                              w_j*(in_val(bo(1,1),myj,myk)+&
                              in_val(bo(1,1)+1,myj,myk)*weights_1d(1))
                         in_val_f=0.0_dp
                      END IF
                   ELSE
                      in_val_f=in_val(bo(1,1),myj,myk)
                      IF (my_sharpen) THEN
                         IF (kw==0.AND.jw==0) THEN
                            IF (my_normalize) THEN
                               out_val(bo(1,1),j,k)=out_val(bo(1,1),j,k)+&
                                    (2.0_dp-w_j)*in_val_f
                            ELSE
                               out_val(bo(1,1),j,k)=out_val(bo(1,1),j,k)+&
                                    in_val_f*w_j
                            END IF
                         ELSE
                            out_val(bo(1,1),j,k)=out_val(bo(1,1),j,k)-&
                                 in_val_f*w_j
                         END IF
                      ELSE
                         out_val(bo(1,1),j,k)=out_val(bo(1,1),j,k)+&
                              in_val_f*w_j
                      END IF
                   END IF
                ELSE
                   in_val_f=l_boundary(myj,myk)
                END IF
                if (has_u_boundary) then
                   IF (my_transpose) THEN
                      IF (s(1)==1) THEN
                         CPAssert(.not.has_l_boundary,cp_failure_level,routineP,error,failure)
                         out_val(bo(2,1),j,k)=out_val(bo(2,1),j,k)+&
                              w_j*(in_val(bo(2,1),myj,myk)+&
                              l_boundary(myj,myk)*weights_1d(1))
                         in_val_l=0.0_dp
                      ELSE
                         in_val_l=0.0_dp
                         out_val(bo(2,1),j,k)=out_val(bo(2,1),j,k)+&
                              w_j*(in_val(bo(2,1),myj,myk)+&
                              in_val(bo(2,1)-1,myj,myk)*weights_1d(1))
                      END IF
                   ELSE
                      in_val_l=in_val(bo(2,1),myj,myk)
                      IF (my_sharpen) THEN
                         IF (kw==0.AND.jw==0) THEN
                            IF (my_normalize) THEN
                               out_val(bo(2,1),j,k)=out_val(bo(2,1),j,k)+&
                                    in_val_l*(2._dp-w_j)
                            ELSE
                               out_val(bo(2,1),j,k)=out_val(bo(2,1),j,k)+&
                                    in_val_l*w_j
                            END IF
                         ELSE
                            out_val(bo(2,1),j,k)=out_val(bo(2,1),j,k)-&
                                 w_j*in_val_l
                         END IF
                      ELSE
                         out_val(bo(2,1),j,k)=out_val(bo(2,1),j,k)+&
                              w_j*in_val_l
                      END IF
                   END IF
                ELSE
                   in_val_l=u_boundary(myj,myk)
                END IF
                IF (last_index>=first_index) THEN
                   IF (my_transpose.AND.bo(1,1)-1==gbo(1,1)) in_val_f=0._dp
                   IF (my_transpose.AND.bo(2,1)+1==gbo(2,1)) in_val_l=0._dp
                   IF (my_sharpen) THEN
                      w=-weights_1d*w_j
                      IF (kw==0.AND.jw==0) THEN
                         IF (my_normalize) THEN
                            w(0)=w(0)+2._dp
                         ELSE
                            w(0)=-w(0)
                         END IF
                      END IF
                   ELSE
                      w=weights_1d*w_j
                   END IF
                   CALL pw_compose_stripe(weights=w,&
                        in_val=in_val(first_index:last_index,myj,myk),&
                        in_val_first=in_val_f,in_val_last=in_val_l,&
                        out_val=out_val(first_index:last_index,j,k),&
                        n_el=n_els)
                END IF
             END DO
          END DO
       END DO
    END DO

    IF (is_split) THEN
       DEALLOCATE(l_boundary,u_boundary,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

  END SUBROUTINE pw_nn_compose_r_no_pbc
  !***************************************************************************

  SUBROUTINE spl3_nopbc(pw_in,pw_out,error)
    TYPE(pw_type), POINTER :: pw_in, pw_out
    TYPE(cp_error_type), INTENT(inout), optional :: error

    CALL pw_zero(pw_out)
    CALL pw_nn_compose_r_no_pbc(weights_1d=spl3_1d_coeffs0,pw_in=pw_in,&
         pw_out=pw_out,sharpen=.FALSE.,normalize=.FALSE.,error=error)
  END SUBROUTINE spl3_nopbc

  SUBROUTINE spl3_nopbct(pw_in,pw_out,error)
    TYPE(pw_type), POINTER :: pw_in, pw_out
    TYPE(cp_error_type), INTENT(inout), optional :: error

    CALL pw_zero(pw_out)
    CALL pw_nn_compose_r_no_pbc(weights_1d=spl3_1d_coeffs0,pw_in=pw_in,&
         pw_out=pw_out,sharpen=.FALSE.,normalize=.FALSE.,transpose=.TRUE.,&
         error=error)
  END SUBROUTINE spl3_nopbct

  SUBROUTINE spl3_pbc(pw_in,pw_out,error)
    TYPE(pw_type), POINTER :: pw_in, pw_out
    TYPE(cp_error_type), INTENT(inout), optional :: error

    CALL pw_zero(pw_out)
    CALL pw_nn_smear_r(pw_in,pw_out,coeffs=spline3_coeffs,&
         error=error)
  END SUBROUTINE spl3_pbc

END MODULE pw_spline_utils
