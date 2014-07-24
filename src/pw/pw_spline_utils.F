!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief different utils that are useful to manipulate splines on the regular
!>      grid of a pw
!> \par History
!>      05.2003 created [fawzi]
!>      08.2004 removed spline evaluation method using more than 2 read streams
!>              (pw_compose_stripe_rs3), added linear solver based spline
!>              inversion [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE pw_spline_utils
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE message_passing,                 ONLY: mp_alltoall,&
                                             mp_comm_compare,&
                                             mp_sendrecv,&
                                             mp_sum
  USE pw_grid_types,                   ONLY: FULLSPACE,&
                                             PW_MODE_LOCAL
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_integral_aa_fast,&
                                             pw_integral_ab_fast,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_release,&
                                             pw_pool_retain,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             SQUARE,&
                                             pw_p_type,&
                                             pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_spline_utils'

  INTEGER, PARAMETER, PUBLIC               :: no_precond=0,&
                                              precond_spl3_aint=1,&
                                              precond_spl3_1=2,&
                                              precond_spl3_aint2=3,&
                                              precond_spl3_2=4,&
                                              precond_spl3_3=5

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

  PUBLIC :: pw_spline3_interpolate_values_g,&
            pw_spline3_deriv_g
  PUBLIC :: pw_spline_scale_deriv
  PUBLIC :: pw_spline2_interpolate_values_g,&
            pw_spline2_deriv_g
  PUBLIC :: pw_nn_smear_r, pw_nn_deriv_r,&
            spl3_nopbc, spl3_pbc, spl3_nopbct
  PUBLIC :: add_fine2coarse, add_coarse2fine
  PUBLIC :: pw_spline_precond_create,&
            pw_spline_do_precond,&
            pw_spline_precond_set_kind,&
            find_coeffs,&
            pw_spline_precond_release,&
            pw_spline_precond_type,&
            Eval_Interp_Spl3_pbc,&
            Eval_d_Interp_Spl3_pbc

!***

! *****************************************************************************
!> \brief stores information for the preconditioner used to calculate the
!>      coeffs of splines
!> \author fawzi
! *****************************************************************************
  TYPE pw_spline_precond_type
     INTEGER :: ref_count,id_nr, kind
     REAL(kind=dp), DIMENSION(4) :: coeffs
     REAL(kind=dp), DIMENSION(3) :: coeffs_1d
     LOGICAL :: sharpen, normalize, pbc, transpose
     TYPE(pw_pool_type), POINTER :: pool
  END TYPE pw_spline_precond_type

CONTAINS

! *****************************************************************************
!> \brief calculates the FFT of the coefficents of the quadratic spline that
!>      interpolates the given values
!> \param spline_g on entry the FFT of the values to interpolate as cc,
!>        will contain the FFT of the coefficents of the spline
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2003 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      does not work with spherical cutoff
! *****************************************************************************
  SUBROUTINE pw_spline2_interpolate_values_g(spline_g,error)
    TYPE(pw_type), POINTER                   :: spline_g
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'pw_spline2_interpolate_values_g', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: gbo
    INTEGER, DIMENSION(3)                    :: n_tot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: c23, coeff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cosIVals, cosJVals, cosKVals

    CALL timeset(routineN,handle)

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

    CALL timestop(handle)
  END SUBROUTINE pw_spline2_interpolate_values_g

! *****************************************************************************
!> \brief calculates the FFT of the coefficents of the2 cubic spline that
!>      interpolates the given values
!> \param spline_g on entry the FFT of the values to interpolate as cc,
!>        will contain the FFT of the coefficents of the spline
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2003 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      does not work with spherical cutoff
!>      stupid distribution for cos calculation, it should calculate only the
!>      needed cos, and avoid the mp_sum
! *****************************************************************************
  SUBROUTINE pw_spline3_interpolate_values_g(spline_g,error)
    TYPE(pw_type), POINTER                   :: spline_g
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'pw_spline3_interpolate_values_g', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: gbo
    INTEGER, DIMENSION(3)                    :: n_tot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: c23, coeff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cosIVals, cosJVals, cosKVals

    CALL timeset(routineN,handle)

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

    CALL timestop(handle)
  END SUBROUTINE pw_spline3_interpolate_values_g

! *****************************************************************************
!> \brief calculates the FFT of the values of the quadratic spline given by the
!>      FFT of the coefficents of the spline (inverse of
!>      pw_spline_interpolate_values_g)
!> \param spline_g on entry the FFT of the coefficents of the spline as cr3d,
!>        will contain the FFT of the values of the spline
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE pw_spline2_evaluate_values_g(spline_g,error)
    TYPE(pw_type), POINTER                   :: spline_g
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_spline2_evaluate_values_g', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    INTEGER, DIMENSION(3)                    :: n, n_tot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: c23, coeff, inv64
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cosIVals, cosJVals, cosKVals

    CALL timeset(routineN,handle)

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

    CALL timestop(handle)
  END SUBROUTINE pw_spline2_evaluate_values_g

! *****************************************************************************
!> \brief calculates the FFT of the values of a nearest neighbor smearing
!> \param spline_g on entry the FFT of the values to interpolate as cc3d,
!>        will contain the FFT of the smeared values
!> \param coeffs array with the coefficent of the smearing, ordered with
!>        the distance from the center: coeffs(1) the coeff of the central
!>        element, coeffs(2) the coeff of the 6 element with distance 1,
!>        coeff(3) the coeff of the 12 elements at distance sqrt(2),
!>        coeff(4) the coeff of the 8 elements at distance sqrt(3).
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2003 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      does not normalize the smear to 1.
!>      with coeff=(/ 8._dp/27._dp, 2._dp/27._dp, 1._dp/54._dp, 1._dp/216._dp /)
!>      is equivalent to pw_spline3_evaluate_values_g, with
!>      coeff=(/ 27._dp/64._dp, 9._dp/128._dp, 3._dp/256._dp, 1._dp/512._dp /)
! *****************************************************************************
  SUBROUTINE pw_nn_smear_g(spline_g,coeffs,error)
    TYPE(pw_type), POINTER                   :: spline_g
    REAL(KIND=dp), DIMENSION(4), INTENT(in)  :: coeffs
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_nn_smear_g', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    INTEGER, DIMENSION(3)                    :: n, n_tot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: c23, coeff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cosIVals, cosJVals, cosKVals
    REAL(KIND=dp), DIMENSION(4)              :: r_coeffs

    CALL timeset(routineN,handle)

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

    CALL timestop(handle)
  END SUBROUTINE pw_nn_smear_g

! *****************************************************************************
!> \brief calculates the FFT of the values of the cubic spline given by the
!>      FFT of the coefficents of the spline (inverse of
!>      pw_spline_interpolate_values_g)
!> \param spline_g on entry the FFT of the values to interpolate as cc3d,
!>        will contain the FFT of the coefficents of the spline
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE pw_spline3_evaluate_values_g(spline_g,error)
    TYPE(pw_type), POINTER                   :: spline_g
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_spline3_evaluate_values_g', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    INTEGER, DIMENSION(3)                    :: n, n_tot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: c23, coeff, inv27
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cosIVals, cosJVals, cosKVals

    CALL timeset(routineN,handle)

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

    CALL timestop(handle)
  END SUBROUTINE pw_spline3_evaluate_values_g

! *****************************************************************************
!> \brief rescales the derivatives from gridspacing=1 to the real derivatives
!> \param deriv_vals_r an array of x,y,z derivatives
!> \param transpose if true applies the transpose of the map (defaults to
!>        false)
!> \param scale a scaling factor (defaults to 1.0)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param cell the mapping to the real cell UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par History
!>      06.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE pw_spline_scale_deriv(deriv_vals_r,transpose,scale,error)
    TYPE(pw_p_type), DIMENSION(3)            :: deriv_vals_r
    LOGICAL, INTENT(in), OPTIONAL            :: transpose
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: scale
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_spline_scale_deriv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, idir, j, k
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(3)                    :: n_tot
    LOGICAL                                  :: diag, failure, my_transpose
    REAL(KIND=dp)                            :: dVal1, dVal2, dVal3, &
                                                my_scale, scalef
    REAL(KIND=dp), DIMENSION(3, 3)           :: dh_inv, h_grid
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: ddata, ddata2, ddata3

    CALL timeset(routineN,handle)

    failure=.FALSE.
    my_transpose=.FALSE.
    IF (PRESENT(transpose)) my_transpose=transpose
    my_scale=1.0_dp
    IF (PRESENT(scale)) my_scale=scale
    n_tot(1:3) = deriv_vals_r(1)%pw%pw_grid%npts (1:3)
    bo = deriv_vals_r(1)%pw%pw_grid%bounds_local
    dh_inv = deriv_vals_r(1)%pw%pw_grid%dh_inv

    ! map grid to real derivative
    diag=.TRUE.
    IF (my_transpose) THEN
       DO j=1,3
          DO i=1,3
             h_grid(j,i)=my_scale*dh_inv(i,j) ! REAL(n_tot(i),dp)*cell_h_inv(i,j)
             IF (i/=j.AND.h_grid(j,i)/=0.0_dp) diag=.FALSE.
          END DO
       END DO
    ELSE
       DO j=1,3
          DO i=1,3
             h_grid(i,j)=my_scale*dh_inv(i,j) ! REAL(n_tot(i),dp)*cell_h_inv(i,j)
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

    CALL timestop(handle)
  END SUBROUTINE pw_spline_scale_deriv

! *****************************************************************************
!> \brief calculates the FFT of the values of the x,y,z (idir=1,2,3)
!>      derivative of the cubic spline
!> \param spline_g on entry the FFT of the coefficents of the spline
!>        will contain the FFT of the derivative
!> \param idir direction of the derivative
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2003 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      the distance between gridpoints is assumed to be 1
! *****************************************************************************
  SUBROUTINE pw_spline3_deriv_g(spline_g,idir,error)
    TYPE(pw_type), POINTER                   :: spline_g
    INTEGER, INTENT(in)                      :: idir
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_spline3_deriv_g', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    INTEGER, DIMENSION(3)                    :: n, n_tot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: coeff, inv9, tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: csIVals, csJVals, csKVals

    CALL timeset(routineN,handle)

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

    CALL timestop(handle)
  END SUBROUTINE pw_spline3_deriv_g

! *****************************************************************************
!> \brief calculates the FFT of the values of the x,y,z (idir=1,2,3)
!>      derivative of the quadratic spline
!> \param spline_g on entry the FFT of the coefficents of the spline
!>        will contain the FFT of the derivative
!> \param idir direction of the derivative
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2003 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      the distance between gridpoints is assumed to be 1
! *****************************************************************************
  SUBROUTINE pw_spline2_deriv_g(spline_g,idir,error)
    TYPE(pw_type), POINTER                   :: spline_g
    INTEGER, INTENT(in)                      :: idir
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_spline2_deriv_g', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(3)                    :: n, n_tot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: coeff, inv16, tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: csIVals, csJVals, csKVals

    CALL timeset(routineN,handle)

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

    CALL timestop(handle)
  END SUBROUTINE pw_spline2_deriv_g

! *****************************************************************************
!> \brief applies a nearest neighbor linear operator to a stripe in x direction:
!>      out_val(i)=sum(weight(j)*in_val(i+j-1),j=0..2)
!> \param weights the weights of the linear operator
!> \param in_val the argument to the operator
!> \param in_val_first the first argument (needed to calculate out_val(1))
!> \param in_val_last the last argument (needed to calculate out_val(n_el))
!> \param out_val the place where the result is accumulated
!> \param n_el the number of elements in in_v and out_v
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par History
!>      04.2004 created [fawzi]
!> \author fawzi
!> \note
!>      uses 2 read streams and 1 write stream
! *****************************************************************************
SUBROUTINE pw_compose_stripe(weights,in_val,in_val_first,in_val_last,&
     out_val,n_el)
    REAL(kind=dp), DIMENSION(0:2), &
      INTENT(in)                             :: weights
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: in_val
    REAL(kind=dp), INTENT(in)                :: in_val_first, in_val_last
    REAL(kind=dp), DIMENSION(*), &
      INTENT(inout)                          :: out_val
    INTEGER                                  :: n_el

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_compose_stripe', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: v0, v1, v2

!1:n_el), &
!1:n_el), &

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
! *****************************************************************************
!> \brief ...
!> \param weights ...
!> \param in_val ...
!> \param in_val_first ...
!> \param in_val_last ...
!> \param out_val ...
!> \param first_val ...
!> \param last_val ...
!> \param myj ...
!> \param myk ...
!> \param j ...
!> \param k ...
! *****************************************************************************
SUBROUTINE pw_compose_stripe2(weights,in_val,in_val_first,in_val_last,&
     out_val,first_val,last_val,myj,myk,j,k)
    REAL(kind=dp), DIMENSION(0:2), &
      INTENT(in)                             :: weights
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: in_val
    REAL(kind=dp), INTENT(in)                :: in_val_first, in_val_last
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: out_val
    INTEGER, INTENT(in)                      :: first_val, last_val, myj, &
                                                myk, j, k

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_compose_stripe2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n_el
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: v0, v1, v2

  IF (last_val-first_val<0) RETURN
  failure=.FALSE.
  v0=in_val_first
  v1=in_val(first_val,myj,myk)
  IF (weights(1)==0.0_dp) THEN
     ! optimized version for x deriv
     DO i=first_val,last_val-3,3
        v2=in_val(i+1,myj,myk)
        out_val(i,j,k)= out_val(i,j,k)+&
             weights(0)*v0+&
             weights(2)*v2
        v0=in_val(i+2,myj,myk)
        out_val(i+1,j,k)= out_val(i+1,j,k)+&
             weights(0)*v1+&
             weights(2)*v0
        v1=in_val(i+3,myj,myk)
        out_val(i+2,j,k)= out_val(i+2,j,k)+&
             weights(0)*v2+&
             weights(2)*v1
     END DO
  ELSE
     ! generic version
     DO i=first_val,last_val-3,3
        v2=in_val(i+1,myj,myk)
        out_val(i,j,k)= out_val(i,j,k)+&
             weights(0)*v0+&
             weights(1)*v1+&
             weights(2)*v2
        v0=in_val(i+2,myj,myk)
        out_val(i+1,j,k)= out_val(i+1,j,k)+&
             weights(0)*v1+&
             weights(1)*v2+&
             weights(2)*v0
        v1=in_val(i+3,myj,myk)
        out_val(i+2,j,k)= out_val(i+2,j,k)+&
             weights(0)*v2+&
             weights(1)*v0+&
             weights(2)*v1
     END DO
  END IF
  n_el=last_val
  SELECT CASE(MODULO(last_val-first_val,3))
  CASE(0)
     v2=in_val_last
     out_val(n_el,j,k)= out_val(n_el,j,k)+&
          weights(0)*v0+&
          weights(1)*v1+&
          weights(2)*v2
  CASE(1)
     v2=in_val(n_el,myj,myk)
     out_val(n_el-1,j,k)= out_val(n_el-1,j,k)+&
          weights(0)*v0+&
          weights(1)*v1+&
          weights(2)*v2
     v0=in_val_last
     out_val(n_el,j,k)= out_val(n_el,j,k)+&
          weights(0)*v1+&
          weights(1)*v2+&
          weights(2)*v0
  CASE(2)
     v2=in_val(n_el-1,myj,myk)
     out_val(n_el-2,j,k)= out_val(n_el-2,j,k)+&
          weights(0)*v0+&
          weights(1)*v1+&
          weights(2)*v2
     v0=in_val(n_el,myj,myk)
     out_val(n_el-1,j,k)= out_val(n_el-1,j,k)+&
          weights(0)*v1+&
          weights(1)*v2+&
          weights(2)*v0
     v1=in_val_last
     out_val(n_el,j,k)= out_val(n_el,j,k)+&
          weights(0)*v2+&
          weights(1)*v0+&
          weights(2)*v1
  END SELECT

END SUBROUTINE pw_compose_stripe2
!**********

! *****************************************************************************
!> \brief private routine that computes pw_nn_compose_r (it seems that without
!>      passing arrays in this way either some compiler do a copyin/out (xlf)
!>      or by inlining suboptimal code is produced (nag))
!> \param weights a 3x3x3 array with the linear operator
!> \param in_val the argument for the linear operator
!> \param out_val place where the value of the linear oprator should be added
!> \param pw_in pw to be able to get the needed meta data about in_val and
!>        out_val
!> \param bo boundaries of in_val and out_val
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
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
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_nn_compose_r_work', &
      routineP = moduleN//':'//routineN

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

! *****************************************************************************
!> \brief applies a nearest neighbor linear operator to a pw in real space
!> \param weights a 3x3x3 array with the linear operator
!> \param pw_in the argument for the linear operator
!> \param pw_out place where the value of the linear oprator should be added
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
!> \note
!>      has specialized versions for derivative operator (with central values==0)
! *****************************************************************************
SUBROUTINE pw_nn_compose_r(weights, pw_in, pw_out,error)
    REAL(kind=dp), DIMENSION(0:2, 0:2, 0:2)  :: weights
    TYPE(pw_type), POINTER                   :: pw_in, pw_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_nn_compose_r', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

  CALL timeset(routineN,handle)
  failure=.FALSE.
  CPPrecondition(ASSOCIATED(pw_in),cp_failure_level,routineP,error,failure)
  CPPrecondition(pw_in%in_space==REALSPACE,cp_failure_level,routineP,error,failure)
  CPPrecondition(pw_in%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(pw_out),cp_failure_level,routineP,error,failure)
  CPPrecondition(pw_out%in_space==REALSPACE,cp_failure_level,routineP,error,failure)
  CPPrecondition(pw_out%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
  IF (.NOT.ALL(pw_in%pw_grid%bounds_local(:,2:3) == pw_in%pw_grid%bounds(:,2:3))) THEN
     CALL cp_assert(.FALSE.,&
          cp_failure_level,cp_assertion_failed,routineP,&
          "wrong pw distribution",error,failure)
  END IF
  IF (.NOT.failure) THEN
     CALL pw_nn_compose_r_work(weights=weights,in_val=pw_in%cr3d,&
          out_val=pw_out%cr3d,pw_in=pw_in,bo=pw_in%pw_grid%bounds_local,&
          error=error)
  END IF
  CALL timestop(handle)
END SUBROUTINE pw_nn_compose_r

! *****************************************************************************
!> \brief calculates the values of a nearest neighbor smearing
!> \param pw_in the argument for the linear operator
!> \param pw_out place where the smeared values should be added
!> \param coeffs array with the coefficent of the smearing, ordered with
!>        the distance from the center: coeffs(1) the coeff of the central
!>        element, coeffs(2) the coeff of the 6 element with distance 1,
!>        coeff(3) the coeff of the 12 elements at distance sqrt(2),
!>        coeff(4) the coeff of the 8 elements at distance sqrt(3).
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Fawzi Mohamed
!> \note
!>      does not normalize the smear to 1.
!>      with coeff=(/ 8._dp/27._dp, 2._dp/27._dp, 1._dp/54._dp, 1._dp/216._dp /)
!>      is equivalent to pw_spline3_evaluate_values_g, with
!>      coeff=(/ 27._dp/64._dp, 9._dp/128._dp, 3._dp/256._dp, 1._dp/512._dp /)
! *****************************************************************************
  SUBROUTINE pw_nn_smear_r(pw_in,pw_out,coeffs,error)
    TYPE(pw_type), POINTER                   :: pw_in, pw_out
    REAL(KIND=dp), DIMENSION(4), INTENT(in)  :: coeffs
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_nn_smear_r', &
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

! *****************************************************************************
!> \brief calculates a nearest neighbor central derivative.
!>      for the x dir:
!>      pw_out%cr3d(i,j,k)=( pw_in(i+1,j,k)-pw_in(i-1,j,k) )*coeff(1)+
!>             ( pw_in(i+1,j(+-)1,k)-pw_in(i-1,j(+-)1,k)+
!>               pw_in(i+1,j,k(+-)1)-pw_in(i-1,j,k(+-)1) )*coeff(2)+
!>             ( pw_in(i+1,j(+-)1,k(+-)1)-pw_in(i-1,j(+-)1,k(+-)1)+
!>               pw_in(i+1,j(+-)1,k(-+)1)-pw_in(i-1,j(+-)1,k(-+)1) )*coeff(3)
!>      periodic boundary conditions are applied
!> \param pw_in the argument for the linear operator
!> \param pw_out place where the smeared values should be added
!> \param coeffs array with the coefficent of the front (positive) plane
!>        of the central derivative, ordered with
!>        the distance from the center: coeffs(1) the coeff of the central
!>        element, coeffs(2) the coeff of the 4 element with distance 1,
!>        coeff(3) the coeff of the 4 elements at distance sqrt(2)
!> \param idir ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param spline_g on entry the FFT of the values to interpolate as cc3d,
!>        will contain the FFT of the smeared values UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author Fawzi Mohamed
!> \note
!>      with coeff=(/ 2.0_dp/9.0_dp,  1.0_dp/18.0_dp, 1.0_dp/72.0_dp /)
!>      is equivalent to pw_spline3_deriv_r, with
!>      coeff=(/ 9.0_dp/32.0_dp, 3.0_dp/64.0_dp, 1.0_dp/128.0_dp /)
!>      to pw_spline2_deriv_r
!>      coeff=(/ 25._dp/72._dp, 5._dp/144, 1._dp/288._dp /)
! *****************************************************************************
  SUBROUTINE pw_nn_deriv_r(pw_in,pw_out,coeffs,idir,error)
    TYPE(pw_type), POINTER                   :: pw_in, pw_out
    REAL(KIND=dp), DIMENSION(3), INTENT(in)  :: coeffs
    INTEGER                                  :: idir
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief low level function that adds a coarse grid
!>      to a fine grid.
!>      If pbc is true periodic boundary conditions are applied
!>
!>      It will add to
!>
!>        fine_values(2*coarse_bounds(1,1):2*coarse_bounds(2,1),
!>                    2*coarse_bounds(1,2):2*coarse_bounds(2,2),
!>                    2*coarse_bounds(1,3):2*coarse_bounds(2,3))
!>
!>      using
!>
!>        coarse_coeffs(coarse_bounds(1,1):coarse_bounds(2,1),
!>                      coarse_bounds(1,2):coarse_bounds(2,2),
!>                      coarse_bounds(1,3):coarse_bounds(2,3))
!>
!>      composed with the weights obtained by the direct product of the
!>      1d coefficents weights:
!>
!>      for i,j,k in -3..3
!>         w(i,j,k)=weights_1d(abs(i)+1)*weights_1d(abs(j)+1)*
!>                  weights_1d(abs(k)+1)
!> \param coarse_coeffs_pw the values of the coefficients
!> \param fine_values_pw where to add the values due to the
!>        coarse coeffs
!> \param weights_1d ...
!> \param w_border0 ...
!> \param w_border1 ...
!> \param pbc if periodic boundary conditions should be applied
!> \param safe_computation ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param w_border 0: the 1d weight at the border (when pbc is false) UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param w_border 1: the 1d weights for a point one off the border
!>        (w_border1(1) is the weight of the coefficent at the border)
!>        (used if pbc is false) UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param weights_ 1d: the weights of the 1d smearing UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author fawzi
!> \note
!>      coarse looping is continuos, I did not check if keeping the fine looping
!>      contiguos is better.
!>      And I ask myself yet again why, why we use x-slice distribution,
!>      z-slice distribution would be much better performancewise
!>      (and would semplify this code enormously).
!>      fine2coarse has much more understandable parallel part (build up of
!>      send/rcv sizes,... but worse if you have really a lot of processors,
!>      probabily irrelevant because it is not critical) [fawzi].
! *****************************************************************************
  SUBROUTINE add_coarse2fine(coarse_coeffs_pw,fine_values_pw,&
       weights_1d,w_border0,w_border1,pbc,safe_computation,error)
    TYPE(pw_type), POINTER                   :: coarse_coeffs_pw, &
                                                fine_values_pw
    REAL(kind=dp), DIMENSION(4), INTENT(in)  :: weights_1d
    REAL(kind=dp), INTENT(in)                :: w_border0
    REAL(kind=dp), DIMENSION(3), INTENT(in)  :: w_border1
    LOGICAL, INTENT(in)                      :: pbc
    LOGICAL, INTENT(in), OPTIONAL            :: safe_computation
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_coarse2fine', &
      routineP = moduleN//':'//routineN

    INTEGER :: coarse_slice_size, f_shift(3), fi, fi_lb, fi_ub, fj, fk, &
      handle, handle2, i, ii, ij, ik, ip, j, k, my_lb, my_ub, n_procs, p, &
      p_lb, p_old, p_ub, rcv_tot_size, rest_b, s(3), send_tot_size, sf, &
      shift, ss, stat, x, x_att, xx
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rcv_offset, rcv_size, &
                                                real_rcv_size, send_offset, &
                                                send_size, sent_size
    INTEGER, DIMENSION(2, 3)                 :: coarse_bo, coarse_gbo, &
                                                fine_bo, fine_gbo, &
                                                my_coarse_bo
    INTEGER, DIMENSION(:), POINTER           :: pos_of_x
    LOGICAL                                  :: failure, has_i_lbound, &
                                                has_i_ubound, is_split, &
                                                safe_calc
    REAL(kind=dp)                            :: v0, v1, v2, v3, wi, wj, wk
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: rcv_buf, send_buf
    REAL(kind=dp), DIMENSION(3)              :: w_0, ww0
    REAL(kind=dp), DIMENSION(4)              :: w_1, ww1
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: coarse_coeffs, fine_values

    CALL timeset(routineN,handle)
!    CALL timeset(routineN//"_pre",handle2)
    failure=.FALSE.
    safe_calc=.FALSE.
    IF (PRESENT(safe_computation)) safe_calc=safe_computation
    CALL mp_comm_compare(coarse_coeffs_pw%pw_grid%para%group,&
         fine_values_pw%pw_grid%para%group,ii)
    IF (ii>1) THEN
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF
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
       IF (fine_bo(1,1)<=fine_bo(2,1)) THEN
          coarse_bo(1,1)=FLOOR((fine_bo(1,1)-2-f_shift(1))/2.)
          coarse_bo(2,1)=FLOOR((fine_bo(2,1)+3-f_shift(1))/2.)
       ELSE
          coarse_bo(1,1)=coarse_gbo(2,1)
          coarse_bo(2,1)=coarse_gbo(2,1)-1
       END IF
       is_split= ANY(coarse_gbo(:,1) /= my_coarse_bo(:,1))
       IF (.NOT.is_split.OR..NOT.pbc) THEN
          coarse_bo(1,1)=MAX(coarse_gbo(1,1),coarse_bo(1,1))
          coarse_bo(2,1)=MIN(coarse_gbo(2,1),coarse_bo(2,1))
       END IF
       has_i_ubound=(fine_gbo(2,1)/=fine_bo(2,1)).or.pbc.and.is_split
       has_i_lbound=(fine_gbo(1,1)/=fine_bo(1,1)).or.pbc.and.is_split

       IF (pbc) THEN
          CPPrecondition(ALL(fine_gbo(1,:)==2*coarse_gbo(1,:)+f_shift),cp_failure_level,routineP,error,failure)
          CPPrecondition(ALL(fine_gbo(2,:)==2*coarse_gbo(2,:)+1+f_shift),cp_failure_level,routineP,error,failure)
       ELSE
          CPPrecondition(ALL(fine_gbo(2,:)==2*coarse_gbo(2,:)+f_shift),cp_failure_level,routineP,error,failure)
          CPPrecondition(ALL(fine_gbo(1,:)==2*coarse_gbo(1,:)+f_shift),cp_failure_level,routineP,error,failure)
       END IF

       coarse_coeffs => coarse_coeffs_pw%cr3d
       DO i=1,3
          s(i)=coarse_gbo(2,i)-coarse_gbo(1,i)+1
       END DO
!       CALL timestop(handle2)
       ! *** parallel case
       IF (is_split) THEN
          CALL timeset(routineN//"_comm",handle2)
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
          ! ** offsets & alloc send-rcv

          send_tot_size=0
          DO ip=0,n_procs-1
             send_offset(ip)=send_tot_size
             send_tot_size=send_tot_size+send_size(ip)
          END DO
          ALLOCATE(send_buf(0:send_tot_size-1),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

          rcv_tot_size=0
          DO ip=0,n_procs-1
             rcv_offset(ip)=rcv_tot_size
             rcv_tot_size=rcv_tot_size+rcv_size(ip)
          END DO
          IF (.NOT.rcv_tot_size==(coarse_bo(2,1)-coarse_bo(1,1)+1)*coarse_slice_size) THEN
             CALL cp_assert(.FALSE.,&
                  cp_failure_level,cp_assertion_failed,routineP,&
                  "Error calculating rcv_tot_size "//&
CPSourceFileRef,&
                  error,failure)
          END IF
          ALLOCATE(rcv_buf(0:rcv_tot_size-1),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

          ! ** fill send buffer

          p_old=pos_of_x(fine_gbo(1,1)+MODULO(fi_lb-fine_gbo(1,1),sf))
          p_lb=FLOOR((fi_lb-2-f_shift(1))/2.)
          sent_size(:)=send_offset
          ss=my_coarse_bo(2,1)-my_coarse_bo(1,1)+1
          DO x=fi_lb,fi_ub
             p=pos_of_x(fine_gbo(1,1)+MODULO(x-fine_gbo(1,1),sf))
             IF (p/=p_old) THEN
                shift=FLOOR((fine_gbo(1,1)+MODULO(x-1-fine_gbo(1,1),sf)-f_shift(1))/2._dp)-&
                     FLOOR((x-1-f_shift(1))/2._dp)
                p_ub=FLOOR((x-1+3-f_shift(1))/2._dp)

                IF (pbc) THEN
                   DO xx=p_lb+shift,coarse_gbo(1,1)-1
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
                      DO i=MAX(p_lb+shift,my_coarse_bo(1,1)),MIN(p_ub+shift,my_coarse_bo(2,1))
                         send_buf(ii)=coarse_coeffs(i,j,k)
                         ii=ii+1
                      END DO
                   END DO
                END DO
                sent_size(p_old)=ii

                IF (pbc) THEN
                   DO xx=coarse_gbo(2,1)+1,p_ub+shift
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
          shift=FLOOR((fine_gbo(1,1)+MODULO(x-1-fine_gbo(1,1),sf)-f_shift(1))/2._dp)-&
               FLOOR((x-1-f_shift(1))/2._dp)
          p_ub=FLOOR((fi_ub+3-f_shift(1))/2.)

          IF (pbc) THEN
             DO xx=p_lb+shift,coarse_gbo(1,1)-1
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
                DO i=MAX(p_lb+shift,my_coarse_bo(1,1)),MIN(p_ub+shift,my_coarse_bo(2,1))
                   send_buf(ii)=coarse_coeffs(i,j,k)
                   ii=ii+1
                END DO
             END DO
          END DO
          sent_size(p_old)=ii

          IF (pbc) THEN
             DO xx=coarse_gbo(2,1)+1,p_ub+shift
                x_att=coarse_gbo(1,1)+MODULO(xx-coarse_gbo(1,1),s(1))
                IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
                   CALL dcopy(coarse_slice_size,&
                        coarse_coeffs(x_att,my_coarse_bo(1,2),&
                        my_coarse_bo(1,3)),ss,send_buf(sent_size(p_old)),1)
                   sent_size(p_old)=sent_size(p_old)+coarse_slice_size
                END IF
             END DO
          END IF

          CPPostcondition(ALL(sent_size(:n_procs-2)==send_offset(1:)),cp_failure_level,routineP,error,failure)
          CPPostcondition(sent_size(n_procs-1)==send_tot_size,cp_failure_level,routineP,error,failure)
          ! test send/rcv sizes
          CALL mp_alltoall(send_size,real_rcv_size,1,coarse_coeffs_pw%pw_grid%para%group)
          CPAssert(ALL(real_rcv_size==rcv_size),cp_failure_level,routineP,error,failure)
          ! all2all
          CALL mp_alltoall( sb=send_buf, scount=send_size, sdispl=send_offset,&
               rb=rcv_buf, rcount=rcv_size, rdispl=rcv_offset, &
               group=coarse_coeffs_pw%pw_grid%para%group )

          ! ** reorder rcv buffer
          ! (actually reordering should be needed only with pbc)

          ALLOCATE(coarse_coeffs(coarse_bo(1,1):coarse_bo(2,1),&
               coarse_bo(1,2):coarse_bo(2,2),&
               coarse_bo(1,3):coarse_bo(2,3)),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

          my_lb=MAX(coarse_gbo(1,1),coarse_bo(1,1))
          my_ub=MIN(coarse_gbo(2,1),coarse_bo(2,1))
          pos_of_x => coarse_coeffs_pw%pw_grid%para%pos_of_x
          sent_size(:)=rcv_offset
          ss=coarse_bo(2,1)-coarse_bo(1,1)+1
          DO x=my_ub+1,coarse_bo(2,1)
             p_old=pos_of_x(coarse_gbo(1,1)+MODULO(x-coarse_gbo(1,1),s(1)))
             CALL dcopy(coarse_slice_size,&
                  rcv_buf(sent_size(p_old)),1,&
                  coarse_coeffs(x,coarse_bo(1,2),&
                  coarse_bo(1,3)),ss)
             sent_size(p_old)=sent_size(p_old)+coarse_slice_size
          END DO
          p_old=pos_of_x(coarse_gbo(1,1) &
               +MODULO(my_lb-coarse_gbo(1,1),s(1)))
          p_lb=my_lb
          DO x=my_lb,my_ub
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
          p_ub=my_ub
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
          DO x=coarse_bo(1,1),my_lb-1
             p_old=pos_of_x(coarse_gbo(1,1)+MODULO(x-coarse_gbo(1,1),s(1)))
             CALL dcopy(coarse_slice_size,&
                  rcv_buf(sent_size(p_old)),1,&
                  coarse_coeffs(x,coarse_bo(1,2),&
                  coarse_bo(1,3)),ss)
             sent_size(p_old)=sent_size(p_old)+coarse_slice_size
          END DO

          CPPostcondition(ALL(sent_size(0:n_procs-2)==rcv_offset(1:)),cp_failure_level,routineP,error,failure)
          CPPostcondition(sent_size(n_procs-1)==rcv_tot_size,cp_failure_level,routineP,error,failure)

          ! dealloc
          DEALLOCATE(send_size,send_offset, rcv_size, rcv_offset, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(send_buf,rcv_buf,real_rcv_size,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          CALL timestop(handle2)

       END IF
       fine_values => fine_values_pw%cr3d
       w_0=(/ weights_1d(3),weights_1d(1),weights_1d(3) /)
       w_1=(/ weights_1d(4),weights_1d(2),weights_1d(2),weights_1d(4) /)

       DO k=coarse_bo(1,3),coarse_bo(2,3)
          DO ik=-3,3
             IF (pbc) THEN
                wk=weights_1d(ABS(ik)+1)
                fk=fine_gbo(1,3)+MODULO(2*k+ik-fine_gbo(1,3)+f_shift(3),2*s(3))
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
                      wj=weights_1d(ABS(ij)+1)*wk
                      fj=fine_gbo(1,2)+MODULO(2*j+ij-fine_gbo(1,2)+f_shift(2),2*s(2))
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
!                      CALL timeset(routineN//"_safe",handle2)
                      DO i=coarse_bo(1,1),coarse_bo(2,1)
                         DO ii=-3,3
                            IF (pbc.and..not.is_split) THEN
                               wi=weights_1d(ABS(ii)+1)*wj
                               fi=fine_gbo(1,1)+MODULO(2*i+ii-fine_gbo(1,1)+f_shift(1),2*s(1))
                            ELSE
                               fi=2*i+ii+f_shift(1)
                               IF (fi<fine_bo(1,1).OR.fi>fine_bo(2,1)) CYCLE
                               IF (.NOT.pbc.AND.(fi<=fine_gbo(1,1)+1.OR. &
                                                 fi>=fine_gbo(2,1)-1)) THEN
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
!                      CALL timestop(handle2)
                   ELSE
!                      CALL timeset(routineN//"_core1",handle2)
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
                         IF (.NOT.(fi+1==fine_bo(1,1).OR.&
                              fi+2==fine_bo(1,1))) THEN
                            CALL cp_assert(.FALSE.,cp_failure_level,&
                                 cp_assertion_failed,routineN,&
                                 "unexpected start index "//&
                                 TRIM(cp_to_string(coarse_bo(1,1)))//" "//&
                                 TRIM(cp_to_string(fi)),error,failure)
                         END IF
                      END IF
                      fi=fi+1
                      IF (fi>=fine_bo(1,1)) THEN
                         fine_values(fi,fj,fk)=fine_values(fi,fj,fk)+&
                              ww0(1)*v0+ww0(2)*v1+&
                              ww0(3)*v2
                      ELSE
                         CPPrecondition(fi+1==fine_bo(1,1),cp_failure_level,routineP,error,failure)
                      END IF
!                      CALL timestop(handle2)
!                      CALL timeset(routineN//"_core",handle2)
                      DO i=coarse_bo(1,1)+3,FLOOR((fine_bo(2,1)-f_shift(1))/2.)-3,4
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
!                      CALL timestop(handle2)
!                      CALL timeset(routineN//"_clean",handle2)
                      rest_b=MODULO(FLOOR((fine_bo(2,1)-f_shift(1))/2.)-coarse_bo(1,1)-3+1,4)
                      IF (rest_b>0) THEN
                         i=FLOOR((fine_bo(2,1)-f_shift(1))/2.)-rest_b+1
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
                         ELSE IF (pbc.and..not.is_split) THEN
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
                      ELSE IF (pbc.and..not.is_split) THEN
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
                      ELSE IF (has_i_ubound) THEN
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
!                   CALL timestop(handle2)
                END DO
             END DO
          END DO
       END DO

       IF (is_split) THEN
          DEALLOCATE(coarse_coeffs,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE add_coarse2fine

! *****************************************************************************
!> \brief low level function that adds a coarse grid (without boundary)
!>      to a fine grid.
!>
!>      It will add to
!>
!>        coarse_coeffs(coarse_bounds(1,1):coarse_bounds(2,1),
!>                      coarse_bounds(1,2):coarse_bounds(2,2),
!>                      coarse_bounds(1,3):coarse_bounds(2,3))
!>
!>      using
!>
!>        fine_values(2*coarse_bounds(1,1):2*coarse_bounds(2,1),
!>                    2*coarse_bounds(1,2):2*coarse_bounds(2,2),
!>                    2*coarse_bounds(1,3):2*coarse_bounds(2,3))
!>
!>      composed with the weights obtained by the direct product of the
!>      1d coefficents weights:
!>
!>      for i,j,k in -3..3
!>         w(i,j,k)=weights_1d(abs(i)+1)*weights_1d(abs(j)+1)*
!>                  weights_1d(abs(k)+1)
!> \param fine_values_pw ...
!> \param coarse_coeffs_pw ...
!> \param weights_1d ...
!> \param w_border0 ...
!> \param w_border1 ...
!> \param pbc ...
!> \param safe_computation ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param coarse_bounds bounds of the coarse coefficents UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param coarse_coeffs 3d array with boundary of size 1 with the values of the
!>        coefficients UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param fine_values 3d array where to add the values due to the
!>        coarse coeffs UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param w_border 0: the 1d weight at the border UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param w_border 1: the 1d weights for a point one off the border
!>        (w_border1(1) is the weight of the coefficent at the border) UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param weights_ 1d: the weights of the 1d smearing UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author fawzi
!> \note
!>      see coarse2fine for some relevant notes
! *****************************************************************************
  SUBROUTINE add_fine2coarse(fine_values_pw,coarse_coeffs_pw,&
       weights_1d,w_border0,w_border1,pbc,safe_computation,error)
    TYPE(pw_type), POINTER                   :: fine_values_pw, &
                                                coarse_coeffs_pw
    REAL(kind=dp), DIMENSION(4), INTENT(in)  :: weights_1d
    REAL(kind=dp), INTENT(in)                :: w_border0
    REAL(kind=dp), DIMENSION(3), INTENT(in)  :: w_border1
    LOGICAL, INTENT(in)                      :: pbc
    LOGICAL, INTENT(in), OPTIONAL            :: safe_computation
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_fine2coarse', &
      routineP = moduleN//':'//routineN

    INTEGER :: coarse_slice_size, f_shift(3), fi, fj, fk, handle, handle2, i, &
      ii, ij, ik, ip, j, k, n_procs, p, p_old, rcv_tot_size, rest_b, s(3), &
      send_tot_size, ss, stat, x, x_att
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: pp_lb, pp_ub, rcv_offset, &
                                                rcv_size, real_rcv_size, &
                                                send_offset, send_size, &
                                                sent_size
    INTEGER, DIMENSION(2, 3)                 :: coarse_bo, coarse_gbo, &
                                                fine_bo, fine_gbo, &
                                                my_coarse_bo
    INTEGER, DIMENSION(:), POINTER           :: pos_of_x
    LOGICAL                                  :: failure, has_i_lbound, &
                                                has_i_ubound, is_split, &
                                                local_data, safe_calc
    REAL(kind=dp)                            :: vv0, vv1, vv2, vv3, vv4, vv5, &
                                                vv6, vv7, wi, wj, wk
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: rcv_buf, send_buf
    REAL(kind=dp), DIMENSION(3)              :: w_0, ww0
    REAL(kind=dp), DIMENSION(4)              :: w_1, ww1
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: coarse_coeffs, fine_values

    failure=.FALSE.
    CALL timeset(routineN,handle)

    IF (.NOT. failure) THEN
       safe_calc=.FALSE.
       IF (PRESENT(safe_computation)) safe_calc=safe_computation

       my_coarse_bo=coarse_coeffs_pw%pw_grid%bounds_local
       coarse_gbo=coarse_coeffs_pw%pw_grid%bounds
       fine_bo=fine_values_pw%pw_grid%bounds_local
       fine_gbo=fine_values_pw%pw_grid%bounds
       f_shift=fine_gbo(1,:)-2*coarse_gbo(1,:)
       is_split= ANY(coarse_gbo(:,1) /= my_coarse_bo(:,1))
       coarse_bo=my_coarse_bo
       IF (fine_bo(1,1)<=fine_bo(2,1)) THEN
          coarse_bo(1,1)=FLOOR(REAL(fine_bo(1,1)-f_shift(1),dp)/2._dp)-1
          coarse_bo(2,1)=FLOOR(REAL(fine_bo(2,1)+1-f_shift(1),dp)/2._dp)+1
       ELSE
          coarse_bo(1,1)=coarse_gbo(2,1)
          coarse_bo(2,1)=coarse_gbo(2,1)-1
       END IF
       IF (.NOT.is_split.OR..NOT.pbc) THEN
          coarse_bo(1,1)=MAX(coarse_gbo(1,1),coarse_bo(1,1))
          coarse_bo(2,1)=MIN(coarse_gbo(2,1),coarse_bo(2,1))
       END IF
       has_i_ubound=(fine_gbo(2,1)/=fine_bo(2,1)).or.pbc.and.is_split
       has_i_lbound=(fine_gbo(1,1)/=fine_bo(1,1)).or.pbc.and.is_split

       IF (pbc) THEN
          CPPrecondition(ALL(fine_gbo(1,:)==2*coarse_gbo(1,:)+f_shift),cp_failure_level,routineP,error,failure)
          CPPrecondition(ALL(fine_gbo(2,:)==2*coarse_gbo(2,:)+f_shift+1),cp_failure_level,routineP,error,failure)
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
             IF (pbc) THEN
                wk=weights_1d(ABS(ik)+1)
                fk=fine_gbo(1,3)+MODULO(2*k+ik-fine_gbo(1,3)+f_shift(3),2*s(3))
             ELSE
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
             END IF
             DO j=coarse_bo(1,2),coarse_bo(2,2)
                DO ij=-3,3
                   IF (pbc) THEN
                      fj=fine_gbo(1,2)+MODULO(2*j+ij-fine_gbo(1,2)+f_shift(2),&
                           2*s(2))
                      wj=weights_1d(ABS(ij)+1)*wk
                   ELSE
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
                   END IF

                   IF (coarse_bo(2,1)-coarse_bo(1,1)<7.or.safe_calc) THEN
                      DO i=coarse_bo(1,1),coarse_bo(2,1)
                         DO ii=-3,3
                            IF (pbc.AND..NOT.is_split) THEN
                               wi=weights_1d(ABS(ii)+1)*wj
                               fi=fine_gbo(1,1)+MODULO(2*i+ii-fine_gbo(1,1)+f_shift(1),2*s(1))
                            ELSE
                               fi=2*i+ii+f_shift(1)
                               IF (fi<fine_bo(1,1).OR.fi>fine_bo(2,1)) CYCLE
                               IF (((.NOT.pbc).AND.fi<=fine_gbo(1,1)+1).OR. &
                                    ((.NOT.pbc).AND.fi>=fine_gbo(2,1)-1)) THEN
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
                         i=coarse_bo(1,1)
                         fi=fine_bo(1,1)-1
                         IF (i+1==FLOOR((fine_bo(1,1)+1-f_shift(1))/2._dp)) THEN
                            fi=fi+1
                            vv0=fine_values(fi,fj,fk)
                            coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)+&

                                 vv0*ww0(3)
                            coarse_coeffs(i+1,j,k)=coarse_coeffs(i+1,j,k)+&

                                 vv0*ww0(2)
                            coarse_coeffs(i+2,j,k)=coarse_coeffs(i+2,j,k)+&

                                 vv0*ww0(1)
                         END IF
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
                      IF (.NOT.FLOOR((fine_bo(2,1)-f_shift(1))/2._dp)-coarse_bo(1,1)>=4) THEN
                         CALL cp_assert(.FALSE.,&
                              cp_failure_level,cp_assertion_failed,routineP,&
                              "FLOOR((fine_bo(2,1)-f_shift(1))/2._dp)-coarse_bo(1,1)>=4",&
                              error,failure)
                      END IF
                      rest_b=MODULO(FLOOR((fine_bo(2,1)-f_shift(1))/2._dp)-coarse_bo(1,1)-6,4)
                      i=FLOOR((fine_bo(2,1)-f_shift(1))/2._dp)-3-rest_b+4
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
                               IF (has_i_ubound) THEN
                                  IF (coarse_bo(2,1)-2==FLOOR((fine_bo(2,1)-f_shift(1))/2._dp)) THEN
                                     fi=fi+1
                                     vv0=fine_values(fi,fj,fk)
                                     coarse_coeffs(i+4,j,k)=coarse_coeffs(i+4,j,k)&
                                          +vv0*ww1(4)
                                  ELSE
                                     vv0=0._dp
                                  END IF
                                  coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                                       +ww1(2)*vv4+ww0(1)*vv5+ww1(1)*vv6
                                  coarse_coeffs(i+1,j,k)=coarse_coeffs(i+1,j,k)&
                                       +ww1(3)*vv4+ww0(2)*vv5+ww1(2)*vv6+ww0(1)*vv7+vv0*ww1(1)
                                  coarse_coeffs(i+2,j,k)=coarse_coeffs(i+2,j,k)&
                                       +ww1(4)*vv4+ww0(3)*vv5+ww1(3)*vv6+ww0(2)*vv7+vv0*ww1(2)
                                  coarse_coeffs(i+3,j,k)=coarse_coeffs(i+3,j,k)&
                                       +ww1(4)*vv6+ww0(3)*vv7+vv0*ww1(3)
                               ELSEIF (pbc.and..not.is_split) THEN
                                  fi=fi+1
                                  vv0=fine_values(fi,fj,fk)
                                  vv1=fine_values(fine_bo(1,1),fj,fk)
                                  vv2=fine_values(fine_bo(1,1)+1,fj,fk)
                                  coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                                       +ww1(2)*vv4+ww0(1)*vv5+ww1(1)*vv6
                                  coarse_coeffs(i+1,j,k)=coarse_coeffs(i+1,j,k)&
                                       +ww1(3)*vv4+ww0(2)*vv5+ww1(2)*vv6+ww0(1)*vv7+vv0*ww1(1)
                                  coarse_coeffs(i+2,j,k)=coarse_coeffs(i+2,j,k)&
                                       +ww1(4)*vv4+ww0(3)*vv5+ww1(3)*vv6+ww0(2)*vv7+vv0*ww1(2)&
                                       +vv1*ww0(1)+vv2*ww1(1)
                               ELSE
                                  coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                                       +ww1(2)*vv4+ww0(1)*vv5+wj*w_border1(3)*vv6
                                  coarse_coeffs(i+1,j,k)=coarse_coeffs(i+1,j,k)&
                                       +ww1(3)*vv4+ww0(2)*vv5+wj*w_border1(2)*vv6
                                  coarse_coeffs(i+2,j,k)=coarse_coeffs(i+2,j,k)&
                                       +ww1(4)*vv4+ww0(3)*vv5+wj*w_border1(1)*vv6+w_border0*wj*vv7
                               END IF
                            ELSE
                               fi=fi+1
                               vv4=fine_values(fi,fj,fk)
                               fi=fi+1
                               vv5=fine_values(fi,fj,fk)
                               IF (has_i_ubound) THEN
                                  IF (coarse_bo(2,1)-2==FLOOR((fine_bo(2,1)-f_shift(1))/2._dp)) THEN
                                     fi=fi+1
                                     vv6=fine_values(fi,fj,fk)
                                     coarse_coeffs(i+3,j,k)=coarse_coeffs(i+3,j,k)&
                                          +ww1(4)*vv6
                                  ELSE
                                     vv6=0._dp
                                  END IF
                                  coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                                       +ww1(1)*vv4
                                  coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                                       +ww1(2)*vv4+ww0(1)*vv5+ww1(1)*vv6
                                  coarse_coeffs(i+1,j,k)=coarse_coeffs(i+1,j,k)&
                                       +ww1(3)*vv4+ww0(2)*vv5+ww1(2)*vv6
                                  coarse_coeffs(i+2,j,k)=coarse_coeffs(i+2,j,k)&
                                       +ww1(4)*vv4+ww0(3)*vv5+ww1(3)*vv6
                               ELSEIF (pbc.AND..not.is_split) THEN
                                  fi=fi+1
                                  vv6=fine_values(fi,fj,fk)
                                  vv7=fine_values(fine_bo(1,1),fj,fk)
                                  vv0=fine_values(fine_bo(1,1)+1,fj,fk)
                                  coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                                       +ww1(1)*vv4
                                  coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                                       +ww1(4)*vv0+ww0(3)*vv1+ww1(3)*vv2+ww0(2)*vv3+ww1(2)*vv4+ww0(1)*vv5+ww1(1)*vv6
                                  coarse_coeffs(i+1,j,k)=coarse_coeffs(i+1,j,k)&
                                       +ww1(4)*vv2+ww0(3)*vv3+ww1(3)*vv4+ww0(2)*vv5+ww1(2)*vv6&
                                       +ww0(1)*vv7+ww1(1)*vv0
                               ELSE
                                  coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                                       +wj*w_border1(3)*vv4
                                  coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                                       +wj*w_border1(2)*vv4
                                  coarse_coeffs(i+1,j,k)=coarse_coeffs(i+1,j,k)&
                                       +wj*(w_border1(1)*vv4+w_border0*vv5)
                               END IF
                            END IF
                         ELSE
                            fi=fi+1
                            vv2=fine_values(fi,fj,fk)
                            fi=fi+1
                            vv3=fine_values(fi,fj,fk)
                            IF (has_i_ubound) THEN
                               IF (coarse_bo(2,1)-2==FLOOR((fine_bo(2,1)-f_shift(1))/2._dp)) THEN
                                  fi=fi+1
                                  vv4=fine_values(fi,fj,fk)
                                  coarse_coeffs(i+2,j,k)=coarse_coeffs(i+2,j,k)&
                                       +ww1(4)*vv4
                               ELSE
                                  vv4=0._dp
                               END IF
                               coarse_coeffs(i-2,j,k)=coarse_coeffs(i-2,j,k)&
                                    +ww1(1)*vv2
                               coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                                    +ww1(2)*vv2+ww0(1)*vv3+ww1(1)*vv4
                               coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                                    +ww1(3)*vv2+ww0(2)*vv3+ww1(2)*vv4
                               coarse_coeffs(i+1,j,k)=coarse_coeffs(i+1,j,k)&
                                    +ww1(4)*vv2+ww0(3)*vv3+ww1(3)*vv4
                            ELSEIF (pbc.AND..not.is_split) THEN
                               fi=fi+1
                               vv4=fine_values(fi,fj,fk)
                               vv5=fine_values(fine_bo(1,1),fj,fk)
                               vv6=fine_values(fine_bo(1,1)+1,fj,fk)
                               coarse_coeffs(i-2,j,k)=coarse_coeffs(i-2,j,k)&
                                    +ww1(1)*vv2
                               coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                                    +ww1(2)*vv2+ww0(1)*vv3+ww1(1)*vv4
                               coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                                    +ww1(3)*vv2+ww0(2)*vv3+ww1(2)*vv4+vv5*ww0(1)+ww1(1)*vv6
                            ELSE
                               coarse_coeffs(i-2,j,k)=coarse_coeffs(i-2,j,k)&
                                    +wj*w_border1(3)*vv2
                               coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                                    +wj*w_border1(2)*vv2
                               coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                                    +wj*(w_border1(1)*vv2+w_border0*vv3)
                            END IF
                         END IF
                      ELSE
                         fi=fi+1
                         vv0=fine_values(fi,fj,fk)
                         fi=fi+1
                         vv1=fine_values(fi,fj,fk)
                         IF (has_i_ubound) THEN
                            IF (coarse_bo(2,1)-2==FLOOR((fine_bo(2,1)-f_shift(1))/2._dp)) THEN
                               fi=fi+1
                               vv2=fine_values(fi,fj,fk)
                               coarse_coeffs(i+1,j,k)=coarse_coeffs(i+1,j,k)&
                                    +ww1(4)*vv2
                            ELSE
                               vv2=0._dp
                            END IF
                            coarse_coeffs(i-3,j,k)=coarse_coeffs(i-3,j,k)&
                                 +ww1(1)*vv0
                            coarse_coeffs(i-2,j,k)=coarse_coeffs(i-2,j,k)&
                                 +ww1(2)*vv0+ww0(1)*vv1+ww1(1)*vv2
                            coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                                 +ww1(3)*vv0+ww0(2)*vv1+ww1(2)*vv2
                            coarse_coeffs(i,j,k)=coarse_coeffs(i,j,k)&
                                 +ww1(4)*vv0+ww0(3)*vv1+ww1(3)*vv2
                         ELSEIF (pbc.AND..not.is_split) THEN
                            fi=fi+1
                            vv2=fine_values(fi,fj,fk)
                            vv3=fine_values(fine_bo(1,1),fk,fk)
                            vv4=fine_values(fine_bo(1,1)+1,fk,fk)
                            coarse_coeffs(i-3,j,k)=coarse_coeffs(i-3,j,k)&
                                 +ww1(1)*vv0
                            coarse_coeffs(i-2,j,k)=coarse_coeffs(i-2,j,k)&
                                 +ww1(2)*vv0+ww0(1)*vv1+ww1(1)*vv2
                            coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                                 +ww1(3)*vv0+ww0(2)*vv1+ww1(2)*vv2+ww0(1)*vv3+ww1(1)*vv4
                         ELSE
                            coarse_coeffs(i-3,j,k)=coarse_coeffs(i-3,j,k)&
                                 +wj*w_border1(3)*vv0
                            coarse_coeffs(i-2,j,k)=coarse_coeffs(i-2,j,k)&
                                 +wj*w_border1(2)*vv0
                            coarse_coeffs(i-1,j,k)=coarse_coeffs(i-1,j,k)&
                                 +wj*(w_border1(1)*vv0+w_border0*vv1)
                         END IF
                      END IF
                      CPPostcondition(fi==fine_bo(2,1),cp_failure_level,routineP,error,failure)
                   END IF
                END DO
             END DO
          END DO
       END DO

       ! *** parallel case
       IF (is_split) THEN
          CALL timeset(routineN//"_comm",handle2)
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

          ! ** rcv size count

          pos_of_x => fine_values_pw%pw_grid%para%pos_of_x
          p_old=pos_of_x(fine_gbo(1,1))
          pp_lb=fine_gbo(2,1)
          pp_ub=fine_gbo(2,1)-1
          pp_lb(p_old)=fine_gbo(1,1)
          DO x=fine_gbo(1,1),fine_gbo(2,1)
             p=pos_of_x(x)
             IF (p/=p_old) THEN
                pp_ub(p_old)=x-1
                pp_lb(p)=x
                p_old=p
             END IF
          END DO
          pp_ub(p_old)=fine_gbo(2,1)

          DO ip=0,n_procs-1
             IF (pp_lb(ip)<=pp_ub(ip)) THEN
                pp_lb(ip)=FLOOR(REAL(pp_lb(ip)-f_shift(1),dp)/2._dp)-1
                pp_ub(ip)=FLOOR(REAL(pp_ub(ip)+1-f_shift(1),dp)/2._dp)+1
             ELSE
                pp_lb(ip)=coarse_gbo(2,1)
                pp_ub(ip)=coarse_gbo(2,1)-1
             END IF
             IF (.NOT.is_split.OR..NOT.pbc) THEN
                pp_lb(ip)=MAX(pp_lb(ip),coarse_gbo(1,1))
                pp_ub(ip)=MIN(pp_ub(ip),coarse_gbo(2,1))
             END IF
          END DO

          rcv_size=0
          DO ip=0,n_procs-1
             DO x=pp_lb(ip),coarse_gbo(1,1)-1
                x_att=coarse_gbo(1,1)+MODULO(x-coarse_gbo(1,1),s(1))
                IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
                   rcv_size(ip)=rcv_size(ip)+coarse_slice_size
                END IF
             END DO
             rcv_size(ip)=rcv_size(ip)+coarse_slice_size*&
               MAX(0,&
               MIN(pp_ub(ip),my_coarse_bo(2,1))-MAX(pp_lb(ip),my_coarse_bo(1,1))+1)
             DO x=coarse_gbo(2,1)+1,pp_ub(ip)
                x_att=coarse_gbo(1,1)+MODULO(x-coarse_gbo(1,1),s(1))
                IF (x_att>=my_coarse_bo(1,1).AND.x_att<=my_coarse_bo(2,1)) THEN
                   rcv_size(ip)=rcv_size(ip)+coarse_slice_size
                END IF
             END DO
          END DO

          ! ** offsets & alloc send-rcv

          send_tot_size=0
          DO ip=0,n_procs-1
             send_offset(ip)=send_tot_size
             send_tot_size=send_tot_size+send_size(ip)
          END DO
          CALL cp_assert(send_tot_size==&
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
          sent_size(:)=send_offset
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

       IF (local_data) THEN
          DEALLOCATE(coarse_coeffs,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       ELSE
          NULLIFY(coarse_coeffs)
       END IF

          CPPostcondition(ALL(sent_size(:n_procs-2)==send_offset(1:)),cp_failure_level,routineP,error,failure)
          CPPostcondition(sent_size(n_procs-1)==send_tot_size,cp_failure_level,routineP,error,failure)
          ! test send/rcv sizes
          CALL mp_alltoall(send_size,real_rcv_size,1,coarse_coeffs_pw%pw_grid%para%group)

          CPAssert(ALL(real_rcv_size==rcv_size),cp_failure_level,routineP,error,failure)
          ! all2all
          CALL mp_alltoall( sb=send_buf, scount=send_size, sdispl=send_offset,&
               rb=rcv_buf, rcount=rcv_size, rdispl=rcv_offset, &
               group=coarse_coeffs_pw%pw_grid%para%group )

          ! ** sum & reorder rcv buffer

          sent_size(:)=rcv_offset
          DO ip=0,n_procs-1

             DO x=pp_lb(ip),coarse_gbo(1,1)-1
                x_att=coarse_gbo(1,1)+MODULO(x-coarse_gbo(1,1),s(1))
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
                x_att=coarse_gbo(1,1)+MODULO(x-coarse_gbo(1,1),s(1))
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

          ! dealloc
          DEALLOCATE(send_size,send_offset, rcv_size, rcv_offset, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(send_buf,rcv_buf,real_rcv_size, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(pp_ub,pp_lb,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          CALL timestop(handle2)
       ELSE
          CPPostcondition(.NOT.local_data,cp_failure_level,routineP,error,failure)
       END IF
    END IF

    CALL timestop(handle)
  END SUBROUTINE add_fine2coarse

! *****************************************************************************
!> \brief ...
!> \param preconditioner the preconditioner to create
!> \param precond_kind the kind of preconditioner to use
!> \param pool a pool with grids of the same type as the elements to
!>        precondition
!> \param pbc if periodic boundary conditions should be applied
!> \param transpose ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pw_spline_precond_create(preconditioner,precond_kind,&
       pool,pbc,transpose,error)
    TYPE(pw_spline_precond_type), POINTER    :: preconditioner
    INTEGER, INTENT(in)                      :: precond_kind
    TYPE(pw_pool_type), POINTER              :: pool
    LOGICAL, INTENT(in)                      :: pbc, transpose
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_spline_precond_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(preconditioner,stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       last_precond_id=last_precond_id+1
       preconditioner%id_nr=last_precond_id
       preconditioner%ref_count=1
       preconditioner%kind=no_precond
       preconditioner%pool => pool
       preconditioner%pbc=pbc
       preconditioner%transpose=transpose
       CALL pw_pool_retain(pool,error=error)
       CALL pw_spline_precond_set_kind(preconditioner,precond_kind,error=error)
    END IF
  END SUBROUTINE pw_spline_precond_create

! *****************************************************************************
!> \brief switches the types of precoditioner to use
!> \param preconditioner ...
!> \param precond_kind the new kind of preconditioner to use
!> \param pbc ...
!> \param transpose ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param precoditioner the preconditioner to be changed UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pw_spline_precond_set_kind(preconditioner,precond_kind,pbc,&
       transpose,error)
    TYPE(pw_spline_precond_type), POINTER    :: preconditioner
    INTEGER, INTENT(in)                      :: precond_kind
    LOGICAL, INTENT(in), OPTIONAL            :: pbc, transpose
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_spline_precond_set_kind', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: do_3d_coeff, failure
    REAL(kind=dp)                            :: s

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(preconditioner),cp_failure_level,routineP,error,failure)
    CPPrecondition(preconditioner%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(transpose)) preconditioner%transpose=transpose
       do_3d_coeff=.FALSE.
       preconditioner%kind=precond_kind
       IF (PRESENT(pbc)) preconditioner%pbc=pbc
       SELECT CASE(precond_kind)
       CASE(no_precond)
       CASE (precond_spl3_aint2)
          preconditioner%coeffs_1d=(/ -1.66_dp*0.25_dp,1.66_dp,-1.66_dp*0.25_dp /)
          preconditioner%sharpen=.FALSE.
          preconditioner%normalize=.FALSE.
          do_3d_coeff=.TRUE.
       CASE (precond_spl3_3)
          preconditioner%coeffs_1d(1)=-0.25_dp*1.6_dp
          preconditioner%coeffs_1d(2)=1.6_dp
          preconditioner%coeffs_1d(3)=-0.25_dp*1.6_dp
          preconditioner%sharpen=.FALSE.
          preconditioner%normalize=.FALSE.
          do_3d_coeff=.TRUE.
       CASE (precond_spl3_2)
          preconditioner%coeffs_1d(1)=-0.26_dp*1.76_dp
          preconditioner%coeffs_1d(2)=1.76_dp
          preconditioner%coeffs_1d(3)=-0.26_dp*1.76_dp
          preconditioner%sharpen=.FALSE.
          preconditioner%normalize=.FALSE.
          do_3d_coeff=.TRUE.
       CASE (precond_spl3_aint)
          preconditioner%coeffs_1d=spl3_1d_coeffs0
          preconditioner%sharpen=.TRUE.
          preconditioner%normalize=.TRUE.
          do_3d_coeff=.TRUE.
       CASE (precond_spl3_1)
          preconditioner%coeffs_1d(1)=0.5_dp/3._dp**(1._dp/3._dp)
          preconditioner%coeffs_1d(2)=4._dp/3._dp**(1._dp/3._dp)
          preconditioner%coeffs_1d(3)=0.5_dp/3._dp**(1._dp/3._dp)
          preconditioner%sharpen=.TRUE.
          preconditioner%normalize=.FALSE.
          do_3d_coeff=.TRUE.
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

! *****************************************************************************
!> \brief retains the preconditioner
!> \param preconditioner the preconditioner to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pw_spline_precond_retain(preconditioner,error)
    TYPE(pw_spline_precond_type), POINTER    :: preconditioner
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_spline_precond_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(preconditioner),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(preconditioner%ref_count>1,cp_failure_level,routineP,error)
       preconditioner%ref_count=preconditioner%ref_count+1
    END IF
  END SUBROUTINE pw_spline_precond_retain

! *****************************************************************************
!> \brief releases the preconditioner
!> \param preconditioner the preconditioner to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pw_spline_precond_release(preconditioner,error)
    TYPE(pw_spline_precond_type), POINTER    :: preconditioner
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_spline_precond_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(preconditioner)) THEN
       CPPreconditionNoFail(preconditioner%ref_count>0,cp_failure_level,routineP,error)
       preconditioner%ref_count=preconditioner%ref_count-1
       IF (preconditioner%ref_count==0) THEN
          CALL pw_pool_release(preconditioner%pool,error=error)
          DEALLOCATE(preconditioner,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
  END SUBROUTINE pw_spline_precond_release

! *****************************************************************************
!> \brief applies the preconditioner to the system of equations to find the
!>      coefficents of the spline
!> \param preconditioner the preconditioner to apply
!> \param in_v ...
!> \param out_v ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param v_in the grid on which the preconditioner should be applied UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param v_out place to store the preconditioner applied on v_out UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pw_spline_do_precond(preconditioner,in_v,out_v,error)
    TYPE(pw_spline_precond_type), POINTER    :: preconditioner
    TYPE(pw_type), POINTER                   :: in_v, out_v
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_spline_do_precond', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(preconditioner),cp_failure_level,routineP,error,failure)
    CPPrecondition(preconditioner%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       SELECT CASE(preconditioner%kind)
       CASE(no_precond)
          CALL pw_copy(in_v,out_v,error=error)
       CASE (precond_spl3_aint,precond_spl3_1)
          CALL pw_zero(out_v,error=error)
          IF (preconditioner%pbc) THEN
             CALL pw_nn_smear_r(pw_in=in_v,pw_out=out_v,&
                  coeffs=preconditioner%coeffs,error=error)
          ELSE
             CALL pw_nn_compose_r_no_pbc(weights_1d=preconditioner%coeffs_1d,&
                  pw_in=in_v, pw_out=out_v, sharpen=preconditioner%sharpen,&
                  normalize=preconditioner%normalize,&
                  transpose=preconditioner%transpose,error=error)
          END IF
       CASE(precond_spl3_3,precond_spl3_2,precond_spl3_aint2)
          CALL pw_zero(out_v,error=error)
          IF (preconditioner%pbc) THEN
             CALL pw_nn_smear_r(pw_in=in_v,pw_out=out_v,&
                  coeffs=preconditioner%coeffs,error=error)
          ELSE
             CALL pw_nn_compose_r_no_pbc(weights_1d=preconditioner%coeffs_1d,&
                  pw_in=in_v, pw_out=out_v, sharpen=preconditioner%sharpen,&
                  normalize=preconditioner%normalize, smooth_boundary=.TRUE.,&
                  transpose=preconditioner%transpose,error=error)
          END IF
       CASE default
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
    END IF
  END SUBROUTINE pw_spline_do_precond

! *****************************************************************************
!> \brief solves iteratively (CG) a systmes of linear equations
!>           linOp(coeffs)=values
!>      (for example those needed to find the coefficents of a spline)
!>      Returns true if the it succeded to acheive the requested accuracy
!> \param values the right hand side of the system
!> \param coeffs ...
!> \param linOp the linear operator to be inverted
!> \param preconditioner the preconditioner to apply
!> \param pool a pool of grids (for the temporary objects)
!> \param eps_r the requested precision on the residual
!> \param eps_x the requested precision on the solution
!> \param max_iter maximum number of iteration allowed
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param coeff will contain the solution of the system (and on entry
!>        it contains the starting point) UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param print_section information about when to print information UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION find_coeffs(values,coeffs,linOp,preconditioner, pool, &
       eps_r,eps_x,max_iter,error) RESULT(res)
    TYPE(pw_type), POINTER                   :: values, coeffs
    INTERFACE
! *****************************************************************************
       SUBROUTINE linOp(pw_in,pw_out,error)
         USE pw_types, ONLY: pw_type
         USE cp_error_handling, ONLY: cp_error_type
         TYPE(pw_type), POINTER :: pw_in,pw_out
         TYPE(cp_error_type), INTENT(inout)  :: error
       END SUBROUTINE linOp
    END INTERFACE
    TYPE(pw_spline_precond_type), POINTER    :: preconditioner
    TYPE(pw_pool_type), POINTER              :: pool
    REAL(kind=dp), INTENT(in)                :: eps_r, eps_x
    INTEGER, INTENT(in)                      :: max_iter
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'find_coeffs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iiter, iter, j, k
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure, last
    REAL(kind=dp)                            :: alpha, beta, eps_r_att, &
                                                eps_x_att, r_z, r_z_new
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(pw_type), POINTER                   :: Ap, p, r, z

    failure=.FALSE.
    last=.FALSE.

    res=.FALSE.
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

       !CALL cp_add_iter_level(logger%iter_info,level_name="SPLINE_FIND_COEFFS",error=error)
       ext_do:DO iiter=1,max_iter,10
          CALL pw_zero(r,error=error)
          CALL linOp(pw_in=coeffs,pw_out=r,error=error)
          r%cr3d=-r%cr3d
          CALL pw_axpy(values,r,error=error)
          CALL pw_spline_do_precond(preconditioner,in_v=r,out_v=z,error=error)
          CALL pw_copy(z,p,error=error)
          r_z=pw_integral_ab_fast(r,z,error=error)

          DO iter=iiter,MIN(iiter+9,max_iter)
             eps_r_att=SQRT(pw_integral_aa_fast(r,SQUARE,error=error))
             IF (eps_r_att==0._dp) THEN
                eps_x_att=0._dp
                last=.TRUE.
             ELSE
                CALL pw_zero(Ap,error=error)
                CALL linOp(pw_in=p,pw_out=Ap,error=error)
                alpha=r_z/pw_integral_ab_fast(Ap,p,error=error)

                CALL pw_axpy(p,coeffs,alpha=alpha,error=error)

                eps_x_att=alpha*SQRT(pw_integral_aa_fast(p,SQUARE,error=error)) ! try to spare if unneded?
                IF (eps_r_att<eps_r .AND. eps_x_att<eps_x) last=.TRUE.
             END IF
             !CALL cp_iterate(logger%iter_info,last=last,error=error)
             IF (last) THEN
                res=.TRUE.
                EXIT ext_do
             END IF

             CALL pw_axpy(Ap,r,alpha=-alpha,error=error)

             CALL pw_spline_do_precond(preconditioner,in_v=r,out_v=z,error=error)

             r_z_new=pw_integral_ab_fast(r,z,error=error)
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
       !CALL cp_rm_iter_level(logger%iter_info,level_name="SPLINE_FIND_COEFFS",error=error)

       CALL pw_pool_give_back_pw(pool,r,error=error)
       CALL pw_pool_give_back_pw(pool,z,error=error)
       CALL pw_pool_give_back_pw(pool,p,error=error)
       CALL pw_pool_give_back_pw(pool,Ap,error=error)
    END IF
  END FUNCTION find_coeffs


! *****************************************************************************
!> \brief adds to pw_out pw_in composed with the weights
!>      pw_out%cr3d(i,j,k)=pw_out%cr3d(i,j,k)+sum(pw_in%cr3d(i+l,j+m,k+n)*
!>         weights_1d(abs(l)+1)*weights_1d(abs(m)+1)*weights_1d(abs(n)+1),
!>         l=-1..1,m=-1..1,n=-1..1)
!> \param weights_1d ...
!> \param pw_in ...
!> \param pw_out ...
!> \param sharpen ...
!> \param normalize ...
!> \param transpose ...
!> \param smooth_boundary ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pw_nn_compose_r_no_pbc(weights_1d,pw_in,pw_out,&
       sharpen,normalize,transpose,smooth_boundary,error)
    REAL(kind=dp), DIMENSION(-1:1)           :: weights_1d
    TYPE(pw_type), POINTER                   :: pw_in, pw_out
    LOGICAL, INTENT(in), OPTIONAL            :: sharpen, normalize, &
                                                transpose, smooth_boundary
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_nn_compose_r_no_pbc', &
      routineP = moduleN//':'//routineN
    REAL(kind=dp), PARAMETER                 :: w_border = 0.0_dp

    INTEGER                                  :: first_index, i, j, jw, k, kw, &
                                                last_index, myj, myk, n_els, &
                                                stat
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    INTEGER, DIMENSION(3)                    :: s
    LOGICAL :: failure, has_l_boundary, has_u_boundary, is_split, &
      my_normalize, my_sharpen, my_smooth_boundary, my_transpose
    REAL(kind=dp)                            :: in_val_f, in_val_l, &
                                                in_val_tmp, w_j, w_k
    REAL(kind=dp), DIMENSION(-1:1)           :: w
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: l_boundary, tmp, u_boundary
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: in_val, out_val

    bo=pw_in%pw_grid%bounds_local
    gbo=pw_in%pw_grid%bounds
    in_val => pw_in%cr3d
    out_val => pw_out%cr3d
    my_sharpen=.FALSE.
    IF (PRESENT(sharpen)) my_sharpen=sharpen
    my_normalize=.FALSE.
    IF (PRESENT(normalize)) my_normalize=normalize
    my_transpose=.FALSE.
    IF (PRESENT(transpose)) my_transpose=transpose
    my_smooth_boundary=.FALSE.
    IF (PRESENT(smooth_boundary)) my_smooth_boundary=smooth_boundary
    CPAssert(.not.my_normalize.OR.my_sharpen,cp_failure_level,routineP,error,failure)
    CPAssert(.NOT.my_smooth_boundary.OR..NOT.my_sharpen,cp_failure_level,routineP,error,failure)
    DO i=1,3
       s(i)=bo(2,i)-bo(1,i)+1
    END DO
    IF (ANY(s<1)) go to 666
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
                      IF (my_smooth_boundary) THEN
                         w_k=weights_1d(kw)/weights_1d(0)
                      END IF
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
             IF (k>=gbo(2,3)-1.OR.k<=gbo(1,3)+1) THEN
                IF (k==gbo(2,3).OR.k==gbo(1,3)) THEN
                   IF (kw/=0) CYCLE
                   w_k=1._dp
                ELSE
                   IF (my_smooth_boundary.AND.((k==gbo(1,3)+1.AND.myk==gbo(1,3)).OR.&
                        (k==gbo(2,3)-1.AND.myk==gbo(2,3)))) THEN
                      w_k=weights_1d(kw)/weights_1d(0)
                   ELSE
                      w_k=weights_1d(kw)
                   END IF
                END IF
             ELSE
                w_k=weights_1d(kw)
             END IF
          END IF
          DO j=bo(1,2),bo(2,2)
             DO jw=-1,1
                myj=j+jw
                IF (j<gbo(2,2)-1.AND.j>gbo(1,2)+1) THEN
                   w_j=w_k*weights_1d(jw)
                ELSE
                   IF (my_transpose) THEN
                      IF (j==gbo(2,2).OR.j==gbo(1,2)) THEN
                         IF (myj<gbo(2,2).AND.myj>gbo(1,2)) THEN
                            w_j=weights_1d(jw)*w_k
                            IF (my_smooth_boundary) THEN
                               w_j=weights_1d(jw)/weights_1d(0)*w_k
                            END IF
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
                      IF (j==gbo(2,2).OR.j==gbo(1,2)) THEN
                         IF (jw/=0) CYCLE
                         w_j=w_k
                      ELSE IF (my_smooth_boundary.AND.((j==gbo(1,2)+1.AND.myj==gbo(1,2)).OR.&
                           (j==gbo(2,2)-1.AND.myj==gbo(2,2)))) THEN
                         w_j=w_k*weights_1d(jw)/weights_1d(0)
                      ELSE
                         w_j=w_k*weights_1d(jw)
                      END IF
                   END IF
                END IF

                IF (has_l_boundary) THEN
                   IF (my_transpose) THEN
                      IF (s(1)==1) THEN
                         CPAssert(.not.has_u_boundary,cp_failure_level,routineP,error,failure)
                         in_val_tmp=u_boundary(myj,myk)
                      ELSE
                         in_val_tmp=in_val(bo(1,1)+1,myj,myk)
                      END IF
                      IF (my_sharpen) THEN
                         IF (kw==0.AND.jw==0) THEN
                            IF (my_normalize) THEN
                               out_val(bo(1,1),j,k)=out_val(bo(1,1),j,k)+&
                                    (2.0_dp-w_j)*in_val(bo(1,1),myj,myk)-&
                                    in_val_tmp*weights_1d(1)*w_j
                            ELSE
                               out_val(bo(1,1),j,k)=out_val(bo(1,1),j,k)+&
                                    in_val(bo(1,1),myj,myk)*w_j-&
                                    in_val_tmp*weights_1d(1)*w_j
                            END IF
                         ELSE
                            out_val(bo(1,1),j,k)=out_val(bo(1,1),j,k)-&
                                 in_val(bo(1,1),myj,myk)*w_j-&
                                 in_val_tmp*weights_1d(1)*w_j
                         END IF
                      ELSE IF (my_smooth_boundary) THEN
                         out_val(bo(1,1),j,k)=out_val(bo(1,1),j,k)+&
                              w_j*(in_val(bo(1,1),myj,myk)+&
                              in_val_tmp*weights_1d(1)/weights_1d(0))
                      ELSE
                         out_val(bo(1,1),j,k)=out_val(bo(1,1),j,k)+&
                              w_j*(in_val(bo(1,1),myj,myk)+&
                              in_val_tmp*weights_1d(1))
                      END IF
                      in_val_f=0.0_dp
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
                IF (has_u_boundary) THEN
                   IF (my_transpose) THEN
                      in_val_l=in_val(bo(2,1),myj,myk)
                      IF (s(1)==1) THEN
                         CPAssert(.not.has_l_boundary,cp_failure_level,routineP,error,failure)
                         in_val_tmp=l_boundary(myj,myk)
                      ELSE
                         in_val_tmp=in_val(bo(2,1)-1,myj,myk)
                      END IF
                      IF (my_sharpen) THEN
                         IF (kw==0.AND.jw==0) THEN
                            IF (my_normalize) THEN
                               out_val(bo(2,1),j,k)=out_val(bo(2,1),j,k)+&
                                    in_val_l*(2._dp-w_j)-&
                                    in_val_tmp*weights_1d(1)*w_j
                            ELSE
                               out_val(bo(2,1),j,k)=out_val(bo(2,1),j,k)+&
                                    in_val_l*w_j-&
                                    in_val_tmp*weights_1d(1)*w_j
                            END IF
                         ELSE
                            out_val(bo(2,1),j,k)=out_val(bo(2,1),j,k)-&
                                 w_j*in_val_l-&
                                    in_val_tmp*weights_1d(1)*w_j
                         END IF
                      ELSE IF (my_smooth_boundary) THEN
                         out_val(bo(2,1),j,k)=out_val(bo(2,1),j,k)+&
                              w_j*(in_val_l+in_val_tmp*weights_1d(1)/weights_1d(0))
                      ELSE
                         out_val(bo(2,1),j,k)=out_val(bo(2,1),j,k)+&
                              w_j*(in_val_l+in_val_tmp*weights_1d(1))
                      END IF
                      in_val_l=0._dp
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
                   IF (my_transpose) THEN
                      IF (bo(1,1)-1==gbo(1,1)) THEN
                         in_val_f=0._dp
                      ELSE IF (bo(2,1)+1==gbo(2,1)) THEN
                         in_val_l=0._dp
                      END IF
                   END IF
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
                   IF (my_smooth_boundary.AND.(.NOT.my_transpose)) THEN
                      IF (gbo(1,1)+1>=bo(1,1).AND.&
                           gbo(1,1)+1<=bo(2,1).AND.gbo(2,1)-gbo(1,1)>2) THEN
                         IF (gbo(1,1)>=bo(1,1)) THEN
                            out_val(gbo(1,1)+1,j,k)=out_val(gbo(1,1)+1,j,k)+&
                                 in_val(gbo(1,1),myj,myk)*w_j*weights_1d(-1)*&
                                 (1._dp/weights_1d(0)-1._dp)
                         ELSE
                            out_val(gbo(1,1)+1,j,k)=out_val(gbo(1,1)+1,j,k)+&
                                 l_boundary(myj,myk)*w_j*weights_1d(-1)*&
                                 (1._dp/weights_1d(0)-1._dp)
                         END IF
                      END IF
                   END IF
                   CALL pw_compose_stripe(weights=w,&
                        in_val=in_val(first_index:last_index,myj,myk),&
                        in_val_first=in_val_f,in_val_last=in_val_l,&
                        out_val=out_val(first_index:last_index,j,k),&
                        n_el=n_els)
!FM                   call pw_compose_stripe2(weights=w,&
!FM                        in_val=in_val,&
!FM                        in_val_first=in_val_f,in_val_last=in_val_l,&
!FM                        out_val=out_val,&
!FM                        first_val=first_index,last_val=last_index,&
!FM                        myj=myj,myk=myk,j=j,k=k)
                   IF (my_smooth_boundary.AND.(.NOT.my_transpose))THEN
                      IF (gbo(2,1)-1>=bo(1,1).AND.&
                           gbo(2,1)-1<=bo(2,1).AND.gbo(2,1)-gbo(1,1)>2) THEN
                         IF (gbo(2,1)<=bo(2,1)) THEN
                            out_val(gbo(2,1)-1,j,k)=out_val(gbo(2,1)-1,j,k)+&
                                 in_val(gbo(2,1),myj,myk)*w_j*weights_1d(1)*&
                                 (1._dp/weights_1d(0)-1._dp)
                         ELSE
                            out_val(gbo(2,1)-1,j,k)=out_val(gbo(2,1)-1,j,k)+&
                                 u_boundary(myj,myk)*w_j*weights_1d(1)*&
                                 (1._dp/weights_1d(0)-1._dp)
                         END IF
                      END IF
                   END IF

                END IF
             END DO
          END DO
       END DO
    END DO

    IF (is_split) THEN
       DEALLOCATE(l_boundary,u_boundary,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
666 CONTINUE
  END SUBROUTINE pw_nn_compose_r_no_pbc

! *****************************************************************************
!> \brief ...
!> \param pw_in ...
!> \param pw_out ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE spl3_nopbc(pw_in,pw_out,error)
    TYPE(pw_type), POINTER                   :: pw_in, pw_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL pw_zero(pw_out,error=error)
    CALL pw_nn_compose_r_no_pbc(weights_1d=spl3_1d_coeffs0,pw_in=pw_in,&
         pw_out=pw_out,sharpen=.FALSE.,normalize=.FALSE.,error=error)
  END SUBROUTINE spl3_nopbc

! *****************************************************************************
!> \brief ...
!> \param pw_in ...
!> \param pw_out ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE spl3_nopbct(pw_in,pw_out,error)
    TYPE(pw_type), POINTER                   :: pw_in, pw_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL pw_zero(pw_out,error=error)
    CALL pw_nn_compose_r_no_pbc(weights_1d=spl3_1d_coeffs0,pw_in=pw_in,&
         pw_out=pw_out,sharpen=.FALSE.,normalize=.FALSE.,transpose=.TRUE.,&
         error=error)
  END SUBROUTINE spl3_nopbct

! *****************************************************************************
!> \brief ...
!> \param pw_in ...
!> \param pw_out ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE spl3_pbc(pw_in,pw_out,error)
    TYPE(pw_type), POINTER                   :: pw_in, pw_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL pw_zero(pw_out,error=error)
    CALL pw_nn_smear_r(pw_in,pw_out,coeffs=spline3_coeffs,&
         error=error)
  END SUBROUTINE spl3_pbc

! *****************************************************************************
!> \brief Evaluates the PBC interpolated Spline (pw) function on the generic
!>      input vector (vec)
!> \param vec ...
!> \param pw ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval val ...
!> \par History
!>      12.2007 Adapted for use with distributed grids [rdeclerck]
!> \author Teodoro Laino 12/2005 [tlaino]
!> \note
!>      Requires the Spline coefficients to be computed with PBC
! *****************************************************************************
  FUNCTION Eval_Interp_Spl3_pbc(vec,pw,error) RESULT(val)
    REAL(KIND=dp), DIMENSION(3), INTENT(in)  :: vec
    TYPE(pw_type), POINTER                   :: pw
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: val

    CHARACTER(len=*), PARAMETER :: routineN = 'Eval_Interp_Spl3_pbc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ivec(3), j, k, npts(3)
    INTEGER, DIMENSION(2, 3)                 :: bo, bo_l
    INTEGER, DIMENSION(4)                    :: ii, ij, ik
    LOGICAL                                  :: failure, my_mpsum
    REAL(KIND=dp) :: a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3, dr1, &
      dr2, dr3, e1, e2, e3, f1, f2, f3, g1, g2, g3, h1, h2, h3, p1, p2, p3, &
      q1, q2, q3, r1, r2, r3, s1, s2, s3, s4, t1, t2, t3, t4, u1, u2, u3, v1, &
      v2, v3, v4, xd1, xd2, xd3
    REAL(KIND=dp), DIMENSION(4, 4, 4)        :: box
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid

    NULLIFY(grid)
    failure=.FALSE.
    IF (.NOT.failure) THEN
       my_mpsum = (pw%pw_grid%para%mode/=PW_MODE_LOCAL)
       npts = pw%pw_grid%npts
       ivec = FLOOR(vec/pw%pw_grid%dr)
       dr1  = pw%pw_grid%dr(1)
       dr2  = pw%pw_grid%dr(2)
       dr3  = pw%pw_grid%dr(3)

       xd1  = (vec(1)/dr1)-REAL(ivec(1),kind=dp)
       xd2  = (vec(2)/dr2)-REAL(ivec(2),kind=dp)
       xd3  = (vec(3)/dr3)-REAL(ivec(3),kind=dp)
       grid => pw%cr3d(:,:,:)
       bo = pw%pw_grid%bounds
       bo_l = pw%pw_grid%bounds_local

       ik(1) = MODULO(ivec(3)-1,npts(3)) + bo(1,3)
       ik(2) = MODULO(ivec(3)  ,npts(3)) + bo(1,3)
       ik(3) = MODULO(ivec(3)+1,npts(3)) + bo(1,3)
       ik(4) = MODULO(ivec(3)+2,npts(3)) + bo(1,3)

       ij(1) = MODULO(ivec(2)-1,npts(2)) + bo(1,2)
       ij(2) = MODULO(ivec(2)  ,npts(2)) + bo(1,2)
       ij(3) = MODULO(ivec(2)+1,npts(2)) + bo(1,2)
       ij(4) = MODULO(ivec(2)+2,npts(2)) + bo(1,2)

       ii(1) = MODULO(ivec(1)-1,npts(1)) + bo(1,1)
       ii(2) = MODULO(ivec(1)  ,npts(1)) + bo(1,1)
       ii(3) = MODULO(ivec(1)+1,npts(1)) + bo(1,1)
       ii(4) = MODULO(ivec(1)+2,npts(1)) + bo(1,1)

       DO k = 1,4
          DO j = 1,4
             DO i = 1,4
                IF(&
                   ii(i) >= bo_l(1,1) .AND. &
                   ii(i) <= bo_l(2,1) .AND. &
                   ij(j) >= bo_l(1,2) .AND. &
                   ij(j) <= bo_l(2,2) .AND. &
                   ik(k) >= bo_l(1,3) .AND. &
                   ik(k) <= bo_l(2,3) &
                ) THEN
                   box(i,j,k) = grid(ii(i) + 1 - bo_l(1,1),&
                                     ij(j) + 1 - bo_l(1,2),&
                                     ik(k) + 1 - bo_l(1,3))
                ELSE
                   box(i,j,k) = 0.0_dp
                END IF
             END DO
          END DO
       END DO

       a1  = 3.0_dp + xd1
       a2  = a1*a1
       a3  = a2*a1
       b1  = 2.0_dp + xd1
       b2  = b1*b1
       b3  = b2*b1
       c1  = 1.0_dp + xd1
       c2  = c1*c1
       c3  = c2*c1
       d1  = xd1
       d2  = d1*d1
       d3  = d2*d1
       e1  = 3.0_dp + xd2
       e2  = e1*e1
       e3  = e2*e1
       f1  = 2.0_dp + xd2
       f2  = f1*f1
       f3  = f2*f1
       g1  = 1.0_dp + xd2
       g2  = g1*g1
       g3  = g2*g1
       h1  = xd2
       h2  = h1*h1
       h3  = h2*h1
       p1  = 3.0_dp + xd3
       p2  = p1*p1
       p3  = p2*p1
       q1  = 2.0_dp + xd3
       q2  = q1*q1
       q3  = q2*q1
       r1  = 1.0_dp + xd3
       r2  = r1*r1
       r3  = r2*r1
       u1  = xd3
       u2  = u1*u1
       u3  = u2*u1

       t1 =   1.0_dp/6.0_dp * (64.0_dp - 48.0_dp*a1 + 12.0_dp*a2 - a3)
       t2 = -22.0_dp/3.0_dp + 10.0_dp*b1 - 4.0_dp*b2 + 0.5_dp*b3
       t3 =   2.0_dp/3.0_dp -  2.0_dp*c1 + 2.0_dp*c2 - 0.5_dp*c3
       t4 =   1.0_dp/6.0_dp*d3
       s1 =   1.0_dp/6.0_dp * (64.0_dp - 48.0_dp*e1 + 12.0_dp*e2 - e3)
       s2 = -22.0_dp/3.0_dp + 10.0_dp*f1 - 4.0_dp*f2 + 0.5_dp*f3
       s3 =   2.0_dp/3.0_dp -  2.0_dp*g1 + 2.0_dp*g2 - 0.5_dp*g3
       s4 =   1.0_dp/6.0_dp*h3
       v1 =   1.0_dp/6.0_dp * (64.0_dp - 48.0_dp*p1 + 12.0_dp*p2 - p3)
       v2 = -22.0_dp/3.0_dp + 10.0_dp*q1 - 4.0_dp*q2 + 0.5_dp*q3
       v3 =   2.0_dp/3.0_dp -  2.0_dp*r1 + 2.0_dp*r2 - 0.5_dp*r3
       v4 =   1.0_dp/6.0_dp*u3

       val = (( box(1,1,1) * t1 + box(2,1,1) * t2 + box(3,1,1) * t3 + box(4,1,1) * t4 ) * s1  +&
              ( box(1,2,1) * t1 + box(2,2,1) * t2 + box(3,2,1) * t3 + box(4,2,1) * t4 ) * s2  +&
              ( box(1,3,1) * t1 + box(2,3,1) * t2 + box(3,3,1) * t3 + box(4,3,1) * t4 ) * s3  +&
              ( box(1,4,1) * t1 + box(2,4,1) * t2 + box(3,4,1) * t3 + box(4,4,1) * t4 ) * s4  ) * v1 +&
             (( box(1,1,2) * t1 + box(2,1,2) * t2 + box(3,1,2) * t3 + box(4,1,2) * t4 ) * s1  +&
              ( box(1,2,2) * t1 + box(2,2,2) * t2 + box(3,2,2) * t3 + box(4,2,2) * t4 ) * s2  +&
              ( box(1,3,2) * t1 + box(2,3,2) * t2 + box(3,3,2) * t3 + box(4,3,2) * t4 ) * s3  +&
              ( box(1,4,2) * t1 + box(2,4,2) * t2 + box(3,4,2) * t3 + box(4,4,2) * t4 ) * s4  ) * v2 +&
             (( box(1,1,3) * t1 + box(2,1,3) * t2 + box(3,1,3) * t3 + box(4,1,3) * t4 ) * s1  +&
              ( box(1,2,3) * t1 + box(2,2,3) * t2 + box(3,2,3) * t3 + box(4,2,3) * t4 ) * s2  +&
              ( box(1,3,3) * t1 + box(2,3,3) * t2 + box(3,3,3) * t3 + box(4,3,3) * t4 ) * s3  +&
              ( box(1,4,3) * t1 + box(2,4,3) * t2 + box(3,4,3) * t3 + box(4,4,3) * t4 ) * s4  ) * v3 +&
             (( box(1,1,4) * t1 + box(2,1,4) * t2 + box(3,1,4) * t3 + box(4,1,4) * t4 ) * s1  +&
              ( box(1,2,4) * t1 + box(2,2,4) * t2 + box(3,2,4) * t3 + box(4,2,4) * t4 ) * s2  +&
              ( box(1,3,4) * t1 + box(2,3,4) * t2 + box(3,3,4) * t3 + box(4,3,4) * t4 ) * s3  +&
              ( box(1,4,4) * t1 + box(2,4,4) * t2 + box(3,4,4) * t3 + box(4,4,4) * t4 ) * s4  ) * v4

       IF (my_mpsum) CALL mp_sum ( val, pw%pw_grid%para%group )
    END IF

  END FUNCTION Eval_Interp_Spl3_pbc

! *****************************************************************************
!> \brief Evaluates the derivatives of the PBC interpolated Spline (pw)
!>      function on the generic input vector (vec)
!> \param vec ...
!> \param pw ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval val ...
!> \par History
!>      12.2007 Adapted for use with distributed grids [rdeclerck]
!> \author Teodoro Laino 12/2005 [tlaino]
!> \note
!>      Requires the Spline coefficients to be computed with PBC
! *****************************************************************************
  FUNCTION Eval_d_Interp_Spl3_pbc(vec,pw,error) RESULT(val)
    REAL(KIND=dp), DIMENSION(3), INTENT(in)  :: vec
    TYPE(pw_type), POINTER                   :: pw
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: val(3)

    CHARACTER(len=*), PARAMETER :: routineN = 'Eval_d_Interp_Spl3_pbc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ivec(3), j, k, npts(3)
    INTEGER, DIMENSION(2, 3)                 :: bo, bo_l
    INTEGER, DIMENSION(4)                    :: ii, ij, ik
    LOGICAL                                  :: failure, my_mpsum
    REAL(KIND=dp) :: a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3, dr1, &
      dr1i, dr2, dr2i, dr3, dr3i, e1, e2, e3, f1, f2, f3, g1, g2, g3, h1, h2, &
      h3, p1, p2, p3, q1, q2, q3, r1, r2, r3, s1, s1d, s1o, s2, s2d, s2o, s3, &
      s3d, s3o, s4, s4d, s4o, t1, t1d, t1o, t2, t2d, t2o, t3, t3d, t3o, t4, &
      t4d, t4o, u1, u2, u3, v1, v1d, v1o, v2, v2d, v2o, v3, v3d, v3o, v4, &
      v4d, v4o, xd1, xd2, xd3
    REAL(KIND=dp), DIMENSION(4, 4, 4)        :: box
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid

    NULLIFY(grid)
    failure=.FALSE.
    IF (.NOT.failure) THEN
       my_mpsum = (pw%pw_grid%para%mode/=PW_MODE_LOCAL)
       npts = pw%pw_grid%npts
       ivec = FLOOR(vec/pw%pw_grid%dr)
       dr1  = pw%pw_grid%dr(1)
       dr2  = pw%pw_grid%dr(2)
       dr3  = pw%pw_grid%dr(3)
       dr1i = 1.0_dp / dr1
       dr2i = 1.0_dp / dr2
       dr3i = 1.0_dp / dr3
       xd1  = (vec(1)/dr1)-REAL(ivec(1),kind=dp)
       xd2  = (vec(2)/dr2)-REAL(ivec(2),kind=dp)
       xd3  = (vec(3)/dr3)-REAL(ivec(3),kind=dp)
       grid => pw%cr3d(:,:,:)
       bo = pw%pw_grid%bounds
       bo_l = pw%pw_grid%bounds_local

       ik(1) = MODULO(ivec(3)-1,npts(3)) + bo(1,3)
       ik(2) = MODULO(ivec(3)  ,npts(3)) + bo(1,3)
       ik(3) = MODULO(ivec(3)+1,npts(3)) + bo(1,3)
       ik(4) = MODULO(ivec(3)+2,npts(3)) + bo(1,3)

       ij(1) = MODULO(ivec(2)-1,npts(2)) + bo(1,2)
       ij(2) = MODULO(ivec(2)  ,npts(2)) + bo(1,2)
       ij(3) = MODULO(ivec(2)+1,npts(2)) + bo(1,2)
       ij(4) = MODULO(ivec(2)+2,npts(2)) + bo(1,2)

       ii(1) = MODULO(ivec(1)-1,npts(1)) + bo(1,1)
       ii(2) = MODULO(ivec(1)  ,npts(1)) + bo(1,1)
       ii(3) = MODULO(ivec(1)+1,npts(1)) + bo(1,1)
       ii(4) = MODULO(ivec(1)+2,npts(1)) + bo(1,1)

       DO k = 1,4
          DO j = 1,4
             DO i = 1,4
                IF(&
                   ii(i) >= bo_l(1,1) .AND. &
                   ii(i) <= bo_l(2,1) .AND. &
                   ij(j) >= bo_l(1,2) .AND. &
                   ij(j) <= bo_l(2,2) .AND. &
                   ik(k) >= bo_l(1,3) .AND. &
                   ik(k) <= bo_l(2,3) &
                ) THEN
                   box(i,j,k) = grid(ii(i) + 1 - bo_l(1,1),&
                                     ij(j) + 1 - bo_l(1,2),&
                                     ik(k) + 1 - bo_l(1,3))
                ELSE
                   box(i,j,k) = 0.0_dp
                END IF
             END DO
          END DO
       END DO

       a1  = 3.0_dp + xd1
       a2  = a1*a1
       a3  = a2*a1
       b1  = 2.0_dp + xd1
       b2  = b1*b1
       b3  = b2*b1
       c1  = 1.0_dp + xd1
       c2  = c1*c1
       c3  = c2*c1
       d1  = xd1
       d2  = d1*d1
       d3  = d2*d1
       e1  = 3.0_dp + xd2
       e2  = e1*e1
       e3  = e2*e1
       f1  = 2.0_dp + xd2
       f2  = f1*f1
       f3  = f2*f1
       g1  = 1.0_dp + xd2
       g2  = g1*g1
       g3  = g2*g1
       h1  = xd2
       h2  = h1*h1
       h3  = h2*h1
       p1  = 3.0_dp + xd3
       p2  = p1*p1
       p3  = p2*p1
       q1  = 2.0_dp + xd3
       q2  = q1*q1
       q3  = q2*q1
       r1  = 1.0_dp + xd3
       r2  = r1*r1
       r3  = r2*r1
       u1  = xd3
       u2  = u1*u1
       u3  = u2*u1

       t1o =   1.0_dp/6.0_dp * (64.0_dp - 48.0_dp*a1 + 12.0_dp*a2 - a3)
       t2o = -22.0_dp/3.0_dp + 10.0_dp*b1 - 4.0_dp*b2 + 0.5_dp*b3
       t3o =   2.0_dp/3.0_dp -  2.0_dp*c1 + 2.0_dp*c2 - 0.5_dp*c3
       t4o =   1.0_dp/6.0_dp*d3
       s1o =   1.0_dp/6.0_dp * (64.0_dp - 48.0_dp*e1 + 12.0_dp*e2 - e3)
       s2o = -22.0_dp/3.0_dp + 10.0_dp*f1 - 4.0_dp*f2 + 0.5_dp*f3
       s3o =   2.0_dp/3.0_dp -  2.0_dp*g1 + 2.0_dp*g2 - 0.5_dp*g3
       s4o =   1.0_dp/6.0_dp*h3
       v1o =   1.0_dp/6.0_dp * (64.0_dp - 48.0_dp*p1 + 12.0_dp*p2 - p3)
       v2o = -22.0_dp/3.0_dp + 10.0_dp*q1 - 4.0_dp*q2 + 0.5_dp*q3
       v3o =   2.0_dp/3.0_dp -  2.0_dp*r1 + 2.0_dp*r2 - 0.5_dp*r3
       v4o =   1.0_dp/6.0_dp*u3

       t1d =  -8.0_dp + 4.0_dp*a1 - 0.5_dp*a2
       t2d =  10.0_dp - 8.0_dp*b1 + 1.5_dp*b2
       t3d =  -2.0_dp + 4.0_dp*c1 - 1.5_dp*c2
       t4d =   0.5_dp*d2
       s1d =  -8.0_dp + 4.0_dp*e1 - 0.5_dp*e2
       s2d =  10.0_dp - 8.0_dp*f1 + 1.5_dp*f2
       s3d =  -2.0_dp + 4.0_dp*g1 - 1.5_dp*g2
       s4d =   0.5_dp*h2
       v1d =  -8.0_dp + 4.0_dp*p1 - 0.5_dp*p2
       v2d =  10.0_dp - 8.0_dp*q1 + 1.5_dp*q2
       v3d =  -2.0_dp + 4.0_dp*r1 - 1.5_dp*r2
       v4d =   0.5_dp*u2

       t1     = t1d*dr1i
       t2     = t2d*dr1i
       t3     = t3d*dr1i
       t4     = t4d*dr1i
       s1     = s1o
       s2     = s2o
       s3     = s3o
       s4     = s4o
       v1     = v1o
       v2     = v2o
       v3     = v3o
       v4     = v4o
       val(1) = (( box(1,1,1) * t1 + box(2,1,1) * t2 + box(3,1,1) * t3 + box(4,1,1) * t4 ) * s1  +&
                 ( box(1,2,1) * t1 + box(2,2,1) * t2 + box(3,2,1) * t3 + box(4,2,1) * t4 ) * s2  +&
                 ( box(1,3,1) * t1 + box(2,3,1) * t2 + box(3,3,1) * t3 + box(4,3,1) * t4 ) * s3  +&
                 ( box(1,4,1) * t1 + box(2,4,1) * t2 + box(3,4,1) * t3 + box(4,4,1) * t4 ) * s4  ) * v1 +&
                (( box(1,1,2) * t1 + box(2,1,2) * t2 + box(3,1,2) * t3 + box(4,1,2) * t4 ) * s1  +&
                 ( box(1,2,2) * t1 + box(2,2,2) * t2 + box(3,2,2) * t3 + box(4,2,2) * t4 ) * s2  +&
                 ( box(1,3,2) * t1 + box(2,3,2) * t2 + box(3,3,2) * t3 + box(4,3,2) * t4 ) * s3  +&
                 ( box(1,4,2) * t1 + box(2,4,2) * t2 + box(3,4,2) * t3 + box(4,4,2) * t4 ) * s4  ) * v2 +&
                (( box(1,1,3) * t1 + box(2,1,3) * t2 + box(3,1,3) * t3 + box(4,1,3) * t4 ) * s1  +&
                 ( box(1,2,3) * t1 + box(2,2,3) * t2 + box(3,2,3) * t3 + box(4,2,3) * t4 ) * s2  +&
                 ( box(1,3,3) * t1 + box(2,3,3) * t2 + box(3,3,3) * t3 + box(4,3,3) * t4 ) * s3  +&
                 ( box(1,4,3) * t1 + box(2,4,3) * t2 + box(3,4,3) * t3 + box(4,4,3) * t4 ) * s4  ) * v3 +&
                (( box(1,1,4) * t1 + box(2,1,4) * t2 + box(3,1,4) * t3 + box(4,1,4) * t4 ) * s1  +&
                 ( box(1,2,4) * t1 + box(2,2,4) * t2 + box(3,2,4) * t3 + box(4,2,4) * t4 ) * s2  +&
                 ( box(1,3,4) * t1 + box(2,3,4) * t2 + box(3,3,4) * t3 + box(4,3,4) * t4 ) * s3  +&
                 ( box(1,4,4) * t1 + box(2,4,4) * t2 + box(3,4,4) * t3 + box(4,4,4) * t4 ) * s4  ) * v4

       t1     = t1o
       t2     = t2o
       t3     = t3o
       t4     = t4o
       s1     = s1d*dr2i
       s2     = s2d*dr2i
       s3     = s3d*dr2i
       s4     = s4d*dr2i
       v1     = v1o
       v2     = v2o
       v3     = v3o
       v4     = v4o
       val(2) = (( box(1,1,1) * t1 + box(2,1,1) * t2 + box(3,1,1) * t3 + box(4,1,1) * t4 ) * s1  +&
                 ( box(1,2,1) * t1 + box(2,2,1) * t2 + box(3,2,1) * t3 + box(4,2,1) * t4 ) * s2  +&
                 ( box(1,3,1) * t1 + box(2,3,1) * t2 + box(3,3,1) * t3 + box(4,3,1) * t4 ) * s3  +&
                 ( box(1,4,1) * t1 + box(2,4,1) * t2 + box(3,4,1) * t3 + box(4,4,1) * t4 ) * s4  ) * v1 +&
                (( box(1,1,2) * t1 + box(2,1,2) * t2 + box(3,1,2) * t3 + box(4,1,2) * t4 ) * s1  +&
                 ( box(1,2,2) * t1 + box(2,2,2) * t2 + box(3,2,2) * t3 + box(4,2,2) * t4 ) * s2  +&
                 ( box(1,3,2) * t1 + box(2,3,2) * t2 + box(3,3,2) * t3 + box(4,3,2) * t4 ) * s3  +&
                 ( box(1,4,2) * t1 + box(2,4,2) * t2 + box(3,4,2) * t3 + box(4,4,2) * t4 ) * s4  ) * v2 +&
                (( box(1,1,3) * t1 + box(2,1,3) * t2 + box(3,1,3) * t3 + box(4,1,3) * t4 ) * s1  +&
                 ( box(1,2,3) * t1 + box(2,2,3) * t2 + box(3,2,3) * t3 + box(4,2,3) * t4 ) * s2  +&
                 ( box(1,3,3) * t1 + box(2,3,3) * t2 + box(3,3,3) * t3 + box(4,3,3) * t4 ) * s3  +&
                 ( box(1,4,3) * t1 + box(2,4,3) * t2 + box(3,4,3) * t3 + box(4,4,3) * t4 ) * s4  ) * v3 +&
                (( box(1,1,4) * t1 + box(2,1,4) * t2 + box(3,1,4) * t3 + box(4,1,4) * t4 ) * s1  +&
                 ( box(1,2,4) * t1 + box(2,2,4) * t2 + box(3,2,4) * t3 + box(4,2,4) * t4 ) * s2  +&
                 ( box(1,3,4) * t1 + box(2,3,4) * t2 + box(3,3,4) * t3 + box(4,3,4) * t4 ) * s3  +&
                 ( box(1,4,4) * t1 + box(2,4,4) * t2 + box(3,4,4) * t3 + box(4,4,4) * t4 ) * s4  ) * v4

       t1     = t1o
       t2     = t2o
       t3     = t3o
       t4     = t4o
       s1     = s1o
       s2     = s2o
       s3     = s3o
       s4     = s4o
       v1     = v1d*dr3i
       v2     = v2d*dr3i
       v3     = v3d*dr3i
       v4     = v4d*dr3i
       val(3) = (( box(1,1,1) * t1 + box(2,1,1) * t2 + box(3,1,1) * t3 + box(4,1,1) * t4 ) * s1  +&
                 ( box(1,2,1) * t1 + box(2,2,1) * t2 + box(3,2,1) * t3 + box(4,2,1) * t4 ) * s2  +&
                 ( box(1,3,1) * t1 + box(2,3,1) * t2 + box(3,3,1) * t3 + box(4,3,1) * t4 ) * s3  +&
                 ( box(1,4,1) * t1 + box(2,4,1) * t2 + box(3,4,1) * t3 + box(4,4,1) * t4 ) * s4  ) * v1 +&
                (( box(1,1,2) * t1 + box(2,1,2) * t2 + box(3,1,2) * t3 + box(4,1,2) * t4 ) * s1  +&
                 ( box(1,2,2) * t1 + box(2,2,2) * t2 + box(3,2,2) * t3 + box(4,2,2) * t4 ) * s2  +&
                 ( box(1,3,2) * t1 + box(2,3,2) * t2 + box(3,3,2) * t3 + box(4,3,2) * t4 ) * s3  +&
                 ( box(1,4,2) * t1 + box(2,4,2) * t2 + box(3,4,2) * t3 + box(4,4,2) * t4 ) * s4  ) * v2 +&
                (( box(1,1,3) * t1 + box(2,1,3) * t2 + box(3,1,3) * t3 + box(4,1,3) * t4 ) * s1  +&
                 ( box(1,2,3) * t1 + box(2,2,3) * t2 + box(3,2,3) * t3 + box(4,2,3) * t4 ) * s2  +&
                 ( box(1,3,3) * t1 + box(2,3,3) * t2 + box(3,3,3) * t3 + box(4,3,3) * t4 ) * s3  +&
                 ( box(1,4,3) * t1 + box(2,4,3) * t2 + box(3,4,3) * t3 + box(4,4,3) * t4 ) * s4  ) * v3 +&
                (( box(1,1,4) * t1 + box(2,1,4) * t2 + box(3,1,4) * t3 + box(4,1,4) * t4 ) * s1  +&
                 ( box(1,2,4) * t1 + box(2,2,4) * t2 + box(3,2,4) * t3 + box(4,2,4) * t4 ) * s2  +&
                 ( box(1,3,4) * t1 + box(2,3,4) * t2 + box(3,3,4) * t3 + box(4,3,4) * t4 ) * s3  +&
                 ( box(1,4,4) * t1 + box(2,4,4) * t2 + box(3,4,4) * t3 + box(4,4,4) * t4 ) * s4  ) * v4

       IF (my_mpsum) CALL mp_sum ( val, pw%pw_grid%para%group )

    END IF

  END FUNCTION Eval_d_Interp_Spl3_pbc

END MODULE pw_spline_utils
