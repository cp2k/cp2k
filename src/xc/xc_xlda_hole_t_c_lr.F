!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates the lda exchange hole in a truncated coulomb potential.
!>        Can be used as longrange correction for truncated hfx calculations
!> \par History
!>      Manuel Guidon (12.2008)  : created
!> \author Manuel Guidon (06.2008)
! *****************************************************************************

MODULE xc_xlda_hole_t_c_lr
  
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE mathlib,                         ONLY: expint
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  PUBLIC :: xlda_hole_t_c_lr_lda_eval, xlda_hole_t_c_lr_lda_info,&
            xlda_hole_t_c_lr_lsd_eval, xlda_hole_t_c_lr_lsd_info,&
            xlda_hole_t_c_lr_lda_calc_0

  REAL(KIND=dp), PARAMETER :: A  = 1.0161144_dp, &
                              B  = -0.37170836_dp, &
                              C  = -0.077215461_dp, &
                              D  = 0.57786348_dp, &
                              E  = -0.051955731_dp

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_xlda_hole_t_c_lr'

  CONTAINS

! *****************************************************************************
!> \brief returns various information on the functional
!> \param xlda_params section selecting the various parameters for the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv controls the number of derivatives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE xlda_hole_t_c_lr_lda_info ( xlda_params, reference, shortform, needs, max_deriv, error)
    TYPE(section_vals_type), POINTER         :: xlda_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    IF ( PRESENT ( reference ) ) THEN
      reference = "{LDA version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "{LDA}"
    END IF
    IF ( PRESENT(needs) ) THEN
       needs%rho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE xlda_hole_t_c_lr_lda_info

! *****************************************************************************
!> \brief returns various information on the functional
!> \param xlda_params section selecting the various parameters for the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv controls the number of derivatives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>        12.2008 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE xlda_hole_t_c_lr_lsd_info ( xlda_params, reference, shortform, needs, max_deriv, error)
    TYPE(section_vals_type), POINTER         :: xlda_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    IF ( PRESENT ( reference ) ) THEN
      reference = "{LSD version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "{LSD}"
    END IF
    IF ( PRESENT(needs) ) THEN
       needs%rho_spin=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE xlda_hole_t_c_lr_lsd_info

! *****************************************************************************
!> \brief evaluates the truncated lda exchange hole
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param order degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param params input parameters (scaling, cutoff_radius)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xlda_hole_t_c_lr_lda_eval ( rho_set, deriv_set, order, params, error )

    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN)                      :: order
    TYPE(section_vals_type), POINTER         :: params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xlda_hole_t_c_lr_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_rho, R, sx
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0, e_rho, rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY(bo)

    CALL section_vals_val_get(params,"SCALE_X",r_val=sx,error=error)
    CALL section_vals_val_get(params,"CUTOFF_RADIUS",r_val=R,error=error)


    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

      CALL xc_rho_set_get(rho_set,rho=rho,&
           local_bounds=bo,rho_cutoff=epsilon_rho,&
           error=error)
      npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

      ! meaningful default for the arrays we don't need: let us make compiler
      ! and debugger happy...
      IF (cp_debug) THEN
        ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ELSE
        dummy=> rho
      END IF

      e_0 => dummy
      e_rho => dummy

      IF (order>=0) THEN
        deriv => xc_dset_get_derivative(deriv_set,"",&
              allocate_deriv=.TRUE., error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
      END IF
      IF (order>=1.OR.order==-1) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
              allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)
      END IF
      IF (order>1.OR.order<-1) THEN
         CALL cp_unimplemented_error(fromWhere=routineP, &
              message="derivatives bigger than 1 not implemented", &
              error=error, error_level=cp_failure_level)
      END IF

      IF ( R == 0.0_dp ) THEN
        CALL cp_unimplemented_error(fromWhere=routineP, &
             message="Cutoff_Radius 0.0 not implemented", &
              error=error, error_level=cp_failure_level)
      END IF
      CALL xlda_hole_t_c_lr_lda_calc(npoints,order,rho=rho,&
                                e_0=e_0,e_rho=e_rho,&
                                epsilon_rho=epsilon_rho,&
                                sx=sx,R=R,error=error)

      IF (cp_debug) THEN
        DEALLOCATE(dummy,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      ELSE
        NULLIFY(dummy)
      END IF
    END IF
    CALL timestop(handle)

  END SUBROUTINE xlda_hole_t_c_lr_lda_eval

! *****************************************************************************
!> \brief Call low level routine
!> \param npoints ...
!> \param order ...
!> \param rho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param epsilon_rho ...
!> \param sx ...
!> \param R ...
!> \param error ...
!> \par History
!>      12.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xlda_hole_t_c_lr_lda_calc(npoints, order, rho,  e_0, e_rho,&
                                  epsilon_rho,sx, R,&
                                  error)

    INTEGER, INTENT(in)                      :: npoints, order
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(inout)                          :: rho, e_0, e_rho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx, R
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xlda_hole_t_c_lr_lda_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(dp)                                 :: my_rho

    !$omp parallel do default(none) &
    !$omp             shared(npoints, rho, epsilon_rho, order, e_0, e_rho) &
    !$omp             shared(sx, r, error) &
    !$omp             private(ip, my_rho)

    DO ip =1,npoints
      my_rho = MAX(rho(ip),0.0_dp)
      IF(my_rho > epsilon_rho) THEN
        CALL xlda_hole_t_c_lr_lda_calc_0(order, my_rho, e_0(ip), e_rho(ip),&
                                    sx, R, error)
      END IF
    END DO

    !$omp end parallel do

  END SUBROUTINE xlda_hole_t_c_lr_lda_calc

! *****************************************************************************
!> \brief low level routine
!> \param order ...
!> \param rho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param sx ...
!> \param R ...
!> \param error ...
!> \par History
!>      12.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xlda_hole_t_c_lr_lda_calc_0(order, rho,  e_0, e_rho,&
                                    sx, R,&
                                    error)
    INTEGER, INTENT(IN)                      :: order
    REAL(KIND=dp), INTENT(IN)                :: rho
    REAL(kind=dp), INTENT(INOUT)             :: e_0, e_rho
    REAL(KIND=dp), INTENT(IN)                :: sx, R
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xlda_hole_t_c_lr_lda_calc_0', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp) :: t1, t12, t14, t15, t19, t2, t22, t23, t24, t25, t3, t32, &
      t33, t36, t4, t41, t46, t5, t6, t62, t64, t67, t68, t7, t82, t86, t9, &
      t91, t95

    IF( order >= 0 ) THEN
      t1 = rho ** 2
      t2 = t1 * pi
      t3 = 3 ** (0.1e1_dp / 0.3e1_dp)
      t4 = pi ** 2
      t5 = t4 * rho
      t6 = t5 ** (0.1e1_dp / 0.3e1_dp)
      t7 = t6 ** 2
      t9 = t3 / t7
      t12 = LOG(R * t3 * t6)
      t14 = R ** 2
      t15 = t14 ** 2
      t19 = 0.1e1_dp / D
      t22 = t3 ** 2
      t23 = t22 * t7
      t24 = D * t14 * t23
      t25 = EXP(-t24)
      t32 = 9 + 4 * A * t14 * t23
      t33 = LOG(t32)
      t36 = D ** 2
      t41 = expint(1, t24)
      t46 = 0.1e1_dp / t36
      t62 = LOG(0.2e1_dp)
      t64 = LOG(A)
      t67 = A * t12 + 0.3e1_dp / 0.2e1_dp * E * t15 * t3 * t6 * t5 * t19 * t25  &
      - A * t33 / 0.2e1_dp + E / t36 / D * t25 + A * t41 / 0.2e1_dp + E * t14  &
      * t22 * t7 * t46 * t25 + B * t19 * t25 / 0.2e1_dp + C * t46 * t25 / 0.2e1_dp  &
      + C * t14 * t22 * t7 * t19 * t25 / 0.2e1_dp + A * t62 + A * t64  &
      / 0.2e1_dp
      t68 = t9 * t67
      e_0 = e_0 + (0.2e1_dp / 0.3e1_dp * t2 * t68) * sx
    END IF
    IF( order >=1 .OR. order ==-1) THEN
      t82 = A / rho
      t86 = t4 ** 2
      t91 = A ** 2
      t95 = 0.1e1_dp / t6 * t4
      e_rho = e_rho + (0.4e1_dp / 0.3e1_dp * rho * pi * t68 - 0.4e1_dp / 0.9e1_dp * t1 * t4 * pi  &
      * t3 / t7 / t5 * t67 + 0.2e1_dp / 0.3e1_dp * t2 * t9 * (t82 / 0.3e1_dp -  &
      0.3e1_dp * E * t15 * t14 * t86 * rho * t25 - 0.4e1_dp / 0.3e1_dp * t91 * t14  &
      * t22 * t95 / t32 - t82 * t25 / 0.3e1_dp - B * t14 * t22 * t95 * t25  &
      / 0.3e1_dp - C * t15 * t3 * t6 * t4 * t25) ) * sx
    END IF


  END SUBROUTINE xlda_hole_t_c_lr_lda_calc_0

! *****************************************************************************
!> \brief evaluates the truncated lsd exchange hole. Calls the lda routine and
!>        applies spin scaling relation
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param order degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param params input parameters (scaling, cutoff_radius)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xlda_hole_t_c_lr_lsd_eval ( rho_set, deriv_set, order, params, error )

    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN)                      :: order
    TYPE(section_vals_type), POINTER         :: params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xlda_hole_t_c_lr_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_rho, R, sx
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0, e_rhoa, e_rhob, &
                                                rhoa, rhob
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY(bo)

    CALL section_vals_val_get(params,"SCALE_X",r_val=sx,error=error)
    CALL section_vals_val_get(params,"CUTOFF_RADIUS",r_val=R,error=error)


    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

      CALL xc_rho_set_get(rho_set,rhoa=rhoa, rhob=rhob,&
           local_bounds=bo,rho_cutoff=epsilon_rho,&
           error=error)
      npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

      ! meaningful default for the arrays we don't need: let us make compiler
      ! and debugger happy...
      IF (cp_debug) THEN
        ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ELSE
        dummy=> rhoa
      END IF

      e_0 => dummy
      e_rhoa => dummy
      e_rhob => dummy

      IF (order>=0) THEN
        deriv => xc_dset_get_derivative(deriv_set,"",&
              allocate_deriv=.TRUE., error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
      END IF
      IF (order>=1.OR.order==-1) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
              allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rhoa,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
              allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rhob,error=error)
      END IF
      IF (order>1.OR.order<-1) THEN
         CALL cp_unimplemented_error(fromWhere=routineP, &
              message="derivatives bigger than 2 not implemented", &
              error=error, error_level=cp_failure_level)
      END IF
      IF ( R == 0.0_dp ) THEN
        CALL cp_unimplemented_error(fromWhere=routineP, &
                       message="Cutoff_Radius 0.0 not implemented", &
                       error=error, error_level=cp_failure_level)
      END IF

      !$omp parallel default(none) &
      !$omp          shared(npoints, order, rhoa, e_0, e_rhoa, epsilon_rho) &
      !$omp          shared(sx, r, error, rhob, e_rhob)

      CALL xlda_hole_t_c_lr_lsd_calc(npoints,order,rho=rhoa,&
                                e_0=e_0,e_rho=e_rhoa,&
                                epsilon_rho=epsilon_rho,&
                                sx=sx,R=R,error=error)

      CALL xlda_hole_t_c_lr_lsd_calc(npoints,order,rho=rhob,&
                                e_0=e_0,e_rho=e_rhob,&
                                epsilon_rho=epsilon_rho,&
                                sx=sx,R=R,error=error)
      !$omp end parallel

      IF (cp_debug) THEN
        DEALLOCATE(dummy,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      ELSE
        NULLIFY(dummy)
      END IF
    END IF
    CALL timestop(handle)

  END SUBROUTINE xlda_hole_t_c_lr_lsd_eval

! *****************************************************************************
!> \brief low level routine
!> \param npoints ...
!> \param order ...
!> \param rho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param epsilon_rho ...
!> \param sx ...
!> \param R ...
!> \param error ...
!> \par History
!>      12.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xlda_hole_t_c_lr_lsd_calc(npoints, order, rho,  e_0, e_rho,&
                                  epsilon_rho,sx, R,&
                                  error)

    INTEGER, INTENT(in)                      :: npoints, order
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(inout)                          :: rho, e_0, e_rho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx, R
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xlda_hole_t_c_lr_lsd_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(dp)                                 :: e_tmp, my_rho

    !$omp do

    DO ip =1,npoints
      my_rho = 2.0_dp*MAX(rho(ip),0.0_dp)
      IF(my_rho > epsilon_rho) THEN
        e_tmp = 0.0_dp
        CALL xlda_hole_t_c_lr_lda_calc_0(order, my_rho, e_tmp , e_rho(ip),&
                                    sx, R, error)
        e_0(ip) = e_0(ip) + 0.5_dp * e_tmp
      END IF
    END DO

    !$omp end do

  END SUBROUTINE xlda_hole_t_c_lr_lsd_calc
END MODULE xc_xlda_hole_t_c_lr

