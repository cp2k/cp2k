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
  USE cp_array_r_utils,                ONLY: cp_3d_r_output
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_check,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: wp=>dp
  USE mathconstants,                   ONLY: twopi
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_compare
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_coeff,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_give_back_coeffs,&
                                             pw_pools_init_coeffs
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             NOSPACE,&
                                             REALDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_allocate,&
                                             pw_copy,&
                                             pw_deallocate,&
                                             pw_derive,&
                                             pw_p_type,&
                                             pw_sumup,&
                                             pw_transfer,&
                                             pw_type,&
                                             pw_zero
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pw_spline_utils'

  PUBLIC :: pw_spline_interpolate_values_g, &
       pw_spline_deriv_r, pw_spline_deriv_g, pw_spline_scale_deriv,&
       pw_spline_test_coeff
!***
!****************************************************************************
CONTAINS

!!****f* pw_spline_utils/pw_spline_interpolate_values_g [1.0] *
!!
!!   NAME
!!     pw_spline_interpolate_values_g
!!
!!   SYNOPSIS
!!     Subroutine pw_spline_interpolate_values_g(spline_g, error)
!!       Type(pw_type), Pointer:: spline_g
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_spline_interpolate_values_g
!!
!!   FUNCTION
!!     calculates the FFT of the coefficents of the spline that interpolates
!!     the given values
!!
!!   NOTES
!!     the distance between gridpoints is assumed to be 1
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
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_spline_interpolate_values_g(spline_g,error)
    TYPE(pw_type),POINTER :: spline_g
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline_interpolate_values_g',&
         routineP=moduleN//':'//routineN
    REAL (kind=wp) :: c23, c_val, coeff, cosKVal
    REAL (kind=wp), DIMENSION(:), ALLOCATABLE :: cosIVals,cosJVals
    INTEGER, DIMENSION(3) :: n,n_tot
    INTEGER :: i, j, k, handle, stat
    INTEGER, DIMENSION (2,3) :: bo
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    n(1:3) = spline_g%pw_grid%npts_local (1:3)
    n_tot(1:3) = spline_g%pw_grid%npts (1:3)
    bo = spline_g%pw_grid%bounds_local

    CPPrecondition(ASSOCIATED(spline_g),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_use==COMPLEXDATA3D,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_space==RECIPROCALSPACE,cp_failure_level,routineP,error,failure)

    ALLOCATE(cosIVals(bo(1,1):bo(2,1)),cosJVals(bo(1,2):bo(2,2)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    coeff=twopi/n_tot(1)
    !$omp parallel do private(i)
    DO i=bo(1,1),bo(2,1)
       cosIVals(i)=COS(coeff*REAL(i+n_tot(1)/2,wp))
    END DO
    coeff=twopi/n_tot(2)
    !$omp parallel do private(j)
    DO j=bo(1,2),bo(2,2)
       cosJVals(j)=COS(coeff*REAL(j+n_tot(2)/2,wp))
    END DO

    !$omp parallel do private(k,j,i,cosKVal,coeff,c23)
    DO k = bo(1,3), bo(2,3)
       cosKVal=COS(twopi/n_tot(3)*REAL(k+n_tot(3)/2,wp))
       DO j = bo(1,2), bo(2,2)
          DO i = bo(1,1), bo(2,1)

             ! no opt
!FM                coeff=1.0/((cosVal(1)*cosVal(2)*cosVal(3))/27.0_wp+&
!FM                     (cosVal(1)*cosVal(2)+cosVal(1)*cosVal(3)+&
!FM                     cosVal(2)*cosVal(3))*2.0_wp/27.0_wp+&
!FM                     (cosVal(1)+cosVal(2)+cosVal(3))*4.0_wp/27.0_wp+&
!FM                     8.0_wp/27.0_wp)
             ! opt
             c23=cosJVals(j)*cosKVal
             coeff=27.0_wp/(cosIVals(i)*c23+&
                  (cosIVals(i)*cosJVals(j)+cosIVals(i)*cosKVal+c23)*2.0_wp+&
                  (cosIVals(i)+cosJVals(j)+cosKVal)*4.0_wp+&
                  8.0_wp)

             spline_g%cc3d(i,j,k)=spline_g%cc3d(i,j,k)*coeff

          END DO
       END DO
    END DO
    DEALLOCATE(cosIVals, cosJVals, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    CALL timestop(0.0_wp,handle)
  END SUBROUTINE pw_spline_interpolate_values_g
!***************************************************************************

!!****f* xc/pw_spline_deriv_r [1.0] *
!!
!!   NAME
!!     pw_spline_deriv_r
!!
!!   SYNOPSIS
!!     Subroutine pw_spline_deriv_r(spline_coeff_r, deriv_val_r, idir,&
!!         error)
!!       Type(pw_type), Pointer:: spline_coeff_r, deriv_val_r
!!       Integer, Intent (IN):: idir
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_spline_deriv_r
!!
!!   FUNCTION
!!     puts the x,y,z derivatives of the given spline
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
  SUBROUTINE pw_spline_deriv_r(spline_coeff_r,deriv_val_r, idir, error)
    TYPE(pw_type),POINTER :: spline_coeff_r, deriv_val_r
    INTEGER, INTENT(in) :: idir
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline_deriv',&
         routineP=moduleN//':'//routineN
    REAL (kind=wp), DIMENSION(3) :: cosVal
    REAL (kind=wp), DIMENSION(3,3) :: deriv_coeff
    REAL (kind=wp) :: c23, c_val, cosKVal
    REAL (kind=wp), DIMENSION(:), ALLOCATABLE :: cosIVals,cosJVals
    INTEGER, DIMENSION(3) :: n,n_tot
    INTEGER :: ii, i, j, k, handle, stat, ij
    INTEGER, DIMENSION (2,3) :: bo
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    deriv_coeff = RESHAPE ((/ &
         1.0_wp/72.0_wp, 1.0_wp/18.0_wp, 1.0_wp/72.0_wp,&
         1.0_wp/18.0_wp, 2.0_wp/9.0_wp,  1.0_wp/18.0_wp,&
         1.0_wp/72.0_wp, 1.0_wp/18.0_wp, 1.0_wp/72.0_wp /),(/3,3/))
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

                c_val=0.0_wp
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

                c_val=0.0_wp
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

                c_val=0.0_wp
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

    CALL timestop(0.0_wp,handle)
  END SUBROUTINE pw_spline_deriv_r
!***************************************************************************


!!****f* pw_spline_utils/pw_spline_scale_deriv [1.0] *
!!
!!   NAME
!!     pw_spline_scale_deriv
!!
!!   SYNOPSIS
!!     Subroutine pw_spline_scale_deriv(deriv_vals_r, cell, transpose,&
!!         scale, error)
!!       Type(pw_p_type), Dimension(3), Intent (IN):: deriv_vals_r
!!       Type(cell_type), Pointer:: cell
!!       Logical, Intent (IN), Optional:: transpose
!!       Real(Kind=wp), Intent (IN), Optional:: scale
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
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_spline_scale_deriv(deriv_vals_r,cell,transpose,scale,error)
    TYPE(pw_p_type), DIMENSION(3), INTENT(in) :: deriv_vals_r
    TYPE(cell_type), POINTER :: cell
    LOGICAL, INTENT(in), OPTIONAL :: transpose
    REAL(kind=wp), INTENT(in), OPTIONAL :: scale
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure, my_transpose
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline_scale_deriv',&
         routineP=moduleN//':'//routineN
    REAL (kind=wp), DIMENSION(3) :: derivVal
    REAL (kind=wp), DIMENSION(3,3) :: h_grid
    REAL(kind=wp) :: my_scale
    INTEGER, DIMENSION(3) :: n,n_tot
    INTEGER :: i, j, k, handle, stat, idir
    INTEGER, DIMENSION (2,3) :: bo

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    my_transpose=.FALSE.
    IF (PRESENT(transpose)) my_transpose=transpose
    my_scale=1.0_wp
    IF (PRESENT(scale)) my_scale=scale
    n_tot(1:3) = deriv_vals_r(1)%pw%pw_grid%npts (1:3)
    bo = deriv_vals_r(1)%pw%pw_grid%bounds_local

    CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)
!FM    CPPrecondition(ASSOCIATED(deriv_vals_r),cp_failure_level,routineP,error,failure)

    ! map grid to real derivative
    IF (my_transpose) THEN
       DO j=1,3
          DO i=1,3
             h_grid(j,i)=my_scale*REAL(n_tot(i),wp)*cell%h_inv(i,j)      
          END DO
       END DO
    ELSE
       DO j=1,3
          DO i=1,3
             h_grid(i,j)=my_scale*REAL(n_tot(i),wp)*cell%h_inv(i,j)      
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

    CALL timestop(0.0_wp,handle)
  END SUBROUTINE pw_spline_scale_deriv
!***************************************************************************

!!****f* pw_spline_utils/pw_spline_deriv_g [1.0] *
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
!!     derivative of the spline
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
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_spline_deriv_g(spline_g,idir,error)
    TYPE(pw_type),POINTER :: spline_g
    INTEGER, INTENT(in) :: idir
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline_deriv_g',&
         routineP=moduleN//':'//routineN
    REAL (kind=wp) :: tmp, c_val, coeff, csKVal, inv9
    REAL (kind=wp), DIMENSION(:), ALLOCATABLE :: csIVals,csJVals
    INTEGER, DIMENSION(3) :: n,n_tot
    INTEGER :: i, j, k, handle, stat
    INTEGER, DIMENSION (2,3) :: bo
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    n(1:3) = spline_g%pw_grid%npts_local (1:3)
    n_tot(1:3) = spline_g%pw_grid%npts (1:3)
    bo = spline_g%pw_grid%bounds_local
    inv9=1.0_wp/9.0_wp

    CPPrecondition(ASSOCIATED(spline_g),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_use==COMPLEXDATA3D,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_g%in_space==RECIPROCALSPACE,cp_failure_level,routineP,error,failure)

    ALLOCATE(csIVals(bo(1,1):bo(2,1)),csJVals(bo(1,2):bo(2,2)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    coeff=twopi/n_tot(1)
    IF (idir==1) THEN
       !$omp parallel do private(i)
       DO i=bo(1,1),bo(2,1)
          csIVals(i)=SIN(coeff*REAL(i+n_tot(1)/2,wp))
       END DO
    ELSE
       !$omp parallel do private(i)
       DO i=bo(1,1),bo(2,1)
          csIVals(i)=COS(coeff*REAL(i+n_tot(1)/2,wp))
       END DO
    END IF
    coeff=twopi/n_tot(2)
    IF (idir==2) THEN
       !$omp parallel do private(j)
       DO j=bo(1,2),bo(2,2)
          csJVals(j)=SIN(coeff*REAL(j+n_tot(2)/2,wp))
       END DO
    ELSE
       !$omp parallel do private(j)
       DO j=bo(1,2),bo(2,2)
          csJVals(j)=COS(coeff*REAL(j+n_tot(2)/2,wp))
       END DO
    END IF

    SELECT CASE(idir)
    CASE (1)
       ! x deriv
       !$omp parallel do private(k,j,i,csVal,coeff,c12)
       DO k = bo(1,3), bo(2,3)
          csKVal=COS(twopi/n_tot(3)*REAL(k+n_tot(3)/2,wp))
          DO j = bo(1,2), bo(2,2)
             DO i = bo(1,1), bo(2,1)

!FM                ! formula
!FM                coeff=(sinVal(1)*cosVal(2)*cosVal(3))/9.0_wp+&
!FM                     (sinVal(1)*cosVal(2)+sinVal(1)*cosVal(3))*2.0_wp/9.0_wp+&
!FM                     sinVal(1)*4.0_wp/9.0_wp

                tmp=csIVals(i)*csJVals(j)
                coeff=(tmp*csKVal+&
                     (tmp+csIVals(i)*csKVal)*2.0_wp+&
                     csIVals(i)*4.0_wp)*inv9

                spline_g%cc3d(i,j,k)=spline_g%cc3d(i,j,k)*&
                     CMPLX(0.0_wp,coeff,wp)

             END DO
          END DO
       END DO
    CASE(2)
       ! y deriv
       !$omp parallel do private(k,j,i,csVal,coeff,c12)
       DO k = bo(1,3), bo(2,3)
          csKVal=COS(twopi/n_tot(3)*REAL(k+n_tot(3)/2,wp))
          DO j = bo(1,2), bo(2,2)
             DO i = bo(1,1), bo(2,1)

                tmp=csIVals(i)*csJVals(j)
                coeff=(tmp*csKVal+&
                     (tmp+csJVals(j)*csKVal)*2.0_wp+&
                     csJVals(j)*4.0_wp)*inv9

                spline_g%cc3d(i,j,k)=spline_g%cc3d(i,j,k)*&
                     CMPLX(0.0_wp,coeff,wp)

             END DO
          END DO
       END DO
    CASE(3)
       ! z deriv
       !$omp parallel do private(k,j,i,csVal,coeff,c12)
       DO k = bo(1,3), bo(2,3)
          csKVal=SIN(twopi/n_tot(3)*REAL(k+n_tot(3)/2,wp))
          DO j = bo(1,2), bo(2,2)
             DO i = bo(1,1), bo(2,1)

                tmp=csIVals(i)*csKVal
                coeff=(tmp*csJVals(j)+&
                     (tmp+csJVals(j)*csKVal)*2.0_wp+&
                     csKVal*4.0_wp)*inv9

                spline_g%cc3d(i,j,k)=spline_g%cc3d(i,j,k)*&
                     CMPLX(0.0_wp,coeff,wp)

             END DO
          END DO
       END DO
    END SELECT

    DEALLOCATE(csIVals, csJVals, stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    CALL timestop(0.0_wp,handle)
  END SUBROUTINE pw_spline_deriv_g
!***************************************************************************

!!****f* xc/pw_spline_test_coeff [1.0] *
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
!!     checks if the coefficents really interpolate the given values
!!     (looks at the differences on the grid)
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
  SUBROUTINE pw_spline_test_coeff(spline_coeff_r, spline_values_r,error)
    TYPE(pw_type), POINTER :: spline_coeff_r, spline_values_r
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_spline_test_coeff',&
         routineP=moduleN//':'//routineN
    REAL (kind=wp), DIMENSION(3,3,3) :: spline_coeff
    REAL(kind=wp) :: maxdiff,diff, c_val
    INTEGER :: i,j,k,ii,ij,ik
    INTEGER, DIMENSION(2,3) :: bo
    INTEGER, DIMENSION(3) :: n_tot

    failure=.FALSE.
    bo = spline_values_r%pw_grid%bounds_local
    n_tot = spline_values_r%pw_grid%npts

    spline_coeff=RESHAPE( (/&
         1._wp/(27._wp*8._wp),1._wp/(27._wp*2._wp),1._wp/(27._wp*8._wp),&
         1._wp/(27._wp*2._wp),2._wp/(27._wp),1._wp/(27._wp*2._wp),&
         1._wp/(27._wp*8._wp),1._wp/(27._wp*2._wp),1._wp/(27._wp*8._wp),&
         1._wp/(27._wp*2._wp),2._wp/(27._wp),1._wp/(27._wp*2._wp),&
         2._wp/(27._wp),8._wp/(27._wp),2._wp/(27._wp),&
         1._wp/(27._wp*2._wp),2._wp/(27._wp),1._wp/(27._wp*2._wp),&
         1._wp/(27._wp*8._wp),1._wp/(27._wp*2._wp),1._wp/(27._wp*8._wp),&
         1._wp/(27._wp*2._wp),2._wp/(27._wp),1._wp/(27._wp*2._wp),&
         1._wp/(27._wp*8._wp),1._wp/(27._wp*2._wp),1._wp/(27._wp*8._wp)&
         /),(/3,3,3/))

    CPPrecondition(ASSOCIATED(spline_coeff_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_coeff_r%in_space==REALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_coeff_r%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(spline_values_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_values_r%in_space==REALSPACE,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_values_r%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
    CPPrecondition(spline_coeff_r%pw_grid%para%group_size==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       maxdiff=0.0_wp
       DO k = 0,spline_values_r%pw_grid%npts_local(3)-1
          DO j = 0,spline_values_r%pw_grid%npts_local(2)-1
             DO i = 0,spline_values_r%pw_grid%npts_local(1)-1

                c_val=0.0_wp
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
  END SUBROUTINE pw_spline_test_coeff
!***************************************************************************

END MODULE pw_spline_utils
