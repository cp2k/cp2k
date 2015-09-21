!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief calculates the Becke 88 exchange functional
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE xc_xbecke88
  USE bibliography,                    ONLY: Becke1988,&
                                             cite_reference
  USE cp_array_utils_r,                ONLY: cp_3d_r_p_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
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

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_xbecke88'
  REAL(kind=dp), PARAMETER :: beta=0.0042_dp

  PUBLIC :: xb88_lda_info, xb88_lsd_info, xb88_lda_eval, xb88_lsd_eval
CONTAINS

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE xb88_lda_info(reference,shortform, needs, max_deriv)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lda_info', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
       reference = "A. Becke, Phys. Rev. A 38, 3098 (1988) {LDA version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "Becke 1988 Exchange Functional (LDA)"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho=.TRUE.
       needs%rho_1_3=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=3

  END SUBROUTINE xb88_lda_info

! *****************************************************************************
!> \brief return various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE xb88_lsd_info(reference,shortform, needs, max_deriv)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lsd_info', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
       reference = "A. Becke, Phys. Rev. A 38, 3098 (1988) {LSD version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "Becke 1988 Exchange Functional (LSD)"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho_spin=.TRUE.
       needs%rho_spin_1_3=.TRUE.
       needs%norm_drho_spin=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=3

  END SUBROUTINE xb88_lsd_info

! *****************************************************************************
!> \brief evaluates the becke 88 exchange functional for lda
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param xb88_params input parameters (scaling)
!> \par History
!>      11.2003 created [fawzi]
!>      01.2007 added scaling [Manuel Guidon]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE xb88_lda_eval(rho_set,deriv_set,grad_deriv,xb88_params)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: xb88_params

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_rho, sx
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndrho, &
      e_ndrho_ndrho, e_ndrho_ndrho_ndrho, e_ndrho_ndrho_rho, e_ndrho_rho, &
      e_ndrho_rho_rho, e_rho, e_rho_rho, e_rho_rho_rho, norm_drho, rho, &
      rho_1_3
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(bo)

    CALL section_vals_val_get(xb88_params,"scale_x",r_val=sx)

    CALL cite_reference(Becke1988)

    CPASSERT(ASSOCIATED(rho_set))
    CPASSERT(rho_set%ref_count>0)
    CPASSERT(ASSOCIATED(deriv_set))
    CPASSERT(deriv_set%ref_count>0)
    CALL xc_rho_set_get(rho_set,rho_1_3=rho_1_3,rho=rho,&
         norm_drho=norm_drho,local_bounds=bo,rho_cutoff=epsilon_rho)
    npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

    ! meaningful default for the arrays we don't need: let us make compiler
    ! and debugger happy...
    IF (cp_debug) THEN
       ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)))
    ELSE
       dummy=> rho
    END IF

    e_0 => dummy
    e_rho => dummy
    e_ndrho => dummy
    e_rho_rho => dummy
    e_ndrho_rho => dummy
    e_ndrho_ndrho => dummy
    e_rho_rho_rho => dummy
    e_ndrho_rho_rho => dummy
    e_ndrho_ndrho_rho => dummy
    e_ndrho_ndrho_ndrho => dummy

    IF (grad_deriv>=0) THEN
       deriv => xc_dset_get_derivative(deriv_set,"",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_0)
    END IF
    IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
       deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho)
       deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho)
    END IF
    IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
       deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho_rho)
       deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rho)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(norm_drho)(norm_drho)", allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho)
    END IF
    IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
       deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)(rho)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(norm_drho)(rho)(rho)",allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho_rho)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(norm_drho)(norm_drho)(rho)",allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_rho)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(norm_drho)(norm_drho)(norm_drho)", allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho)
    END IF
    IF (grad_deriv>3.OR.grad_deriv<-3) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="derivatives bigger than 3 not implemented", &
            error_level=cp_failure_level)
    END IF

    !$omp parallel default(none) &
    !$omp          shared(rho, rho_1_3, norm_drho, e_0, e_rho) &
    !$omp          shared(e_ndrho, e_rho_rho, e_ndrho_rho) &
    !$omp          shared(e_ndrho_ndrho, e_rho_rho_rho) &
    !$omp          shared(e_ndrho_rho_rho, e_ndrho_ndrho_rho) &
    !$omp          shared(e_ndrho_ndrho_ndrho, grad_deriv, npoints) &
    !$omp          shared(epsilon_rho, sx)

    CALL xb88_lda_calc(rho=rho, rho_1_3=rho_1_3, norm_drho=norm_drho,&
         e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,e_rho_rho=e_rho_rho,&
         e_ndrho_rho=e_ndrho_rho, e_ndrho_ndrho=e_ndrho_ndrho, &
         e_rho_rho_rho=e_rho_rho_rho, e_ndrho_rho_rho=e_ndrho_rho_rho,&
         e_ndrho_ndrho_rho=e_ndrho_ndrho_rho,&
         e_ndrho_ndrho_ndrho=e_ndrho_ndrho_ndrho,grad_deriv=grad_deriv,&
         npoints=npoints,epsilon_rho=epsilon_rho,sx=sx)

    !$omp end parallel

    IF (cp_debug) THEN
       DEALLOCATE(dummy)
    ELSE
       NULLIFY(dummy)
    END IF
    CALL timestop(handle)
  END SUBROUTINE xb88_lda_eval

! *****************************************************************************
!> \brief evaluates the becke 88 exchange functional for lda
!> \param rho the density where you want to evaluate the functional
!> \param rho_1_3 ...
!> \param norm_drho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_rho_rho ...
!> \param e_ndrho_rho ...
!> \param e_ndrho_ndrho ...
!> \param e_rho_rho_rho ...
!> \param e_ndrho_rho_rho ...
!> \param e_ndrho_ndrho_rho ...
!> \param e_ndrho_ndrho_ndrho ...
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param npoints ...
!> \param epsilon_rho ...
!> \param sx scaling-parameter for exchange
!> \par History
!>      11.2003 created [fawzi]
!>      01.2007 added scaling [Manuel Guidon]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE xb88_lda_calc(rho, rho_1_3, norm_drho,&
       e_0,e_rho,e_ndrho,e_rho_rho,e_ndrho_rho,&
       e_ndrho_ndrho, e_rho_rho_rho, e_ndrho_rho_rho, e_ndrho_ndrho_rho,&
       e_ndrho_ndrho_ndrho,grad_deriv,npoints,epsilon_rho,sx)
    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), DIMENSION(1:npoints), INTENT(inout) :: &
      e_ndrho_ndrho_ndrho, e_ndrho_ndrho_rho, e_ndrho_rho_rho, e_rho_rho_rho, &
      e_ndrho_ndrho, e_ndrho_rho, e_rho_rho, e_ndrho, e_rho, e_0
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(in)                             :: norm_drho, rho_1_3, rho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lda_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii
    LOGICAL                                  :: failure
    REAL(kind=dp) :: c, epsilon_rho43, my_rho, my_rho_1_3, t1, t10, t100, &
      t104, t11, t12, t126, t13, t159, t164, t17, t170, t176, t18, t189, t2, &
      t21, t23, t29, t3, t31, t33, t37, t39, t40, t43, t44, t49, t5, t51, &
      t54, t6, t64, t67, t7, t72, t74, t75, t79, t83, t86, t90, t91, t98, &
      t99, x

    failure=.FALSE.
    c= 1.5_dp*(0.75_dp/Pi)**(1.0_dp/3.0_dp)
    t1 = 2 ** (0.1e1_dp / 0.3e1_dp)
    t2 = t1 ** 2
    t5 = beta * t2
    t7 = beta * t1
    epsilon_rho43=epsilon_rho**(4./3.)

    SELECT CASE(grad_deriv)
    CASE (0)

       !$omp do

       DO ii=1,npoints

          my_rho = rho(ii)
          IF (my_rho>epsilon_rho) THEN
             my_rho_1_3 = rho_1_3(ii)
             t3 = my_rho_1_3 * my_rho
             x = norm_drho(ii)/MAX(t3,epsilon_rho43)
             t6 = x ** 2
             t10 = t2 * t6 + 0.1e1_dp
             t11 = SQRT(t10)
             t12 = t1 * x + t11
             t13 = LOG(t12)
             t17 = 0.1e1_dp + 0.6e1_dp * t7 * x * t13
             t18 = 0.1e1_dp / t17
             t21 = -c - t5 * t6 * t18

             e_0(ii) = e_0(ii) + ( t2 * t3 * t21 / 0.2e1_dp ) * sx
          END IF
       END DO

       !$omp end do

    CASE (1)

       !$omp do

       DO ii=1,npoints

          my_rho = rho(ii)
          IF (my_rho > epsilon_rho) THEN
             my_rho_1_3 = rho_1_3(ii)
             t3 = my_rho_1_3 * my_rho
             x = norm_drho(ii)/t3
             t6 = x ** 2
             t10 = t2 * t6 + 0.1e1_dp
             t11 = SQRT(t10)
             t12 = t1 * x + t11
             t13 = LOG(t12)
             t17 = 0.1e1_dp + 0.6e1_dp * t7 * x * t13
             t18 = 0.1e1_dp / t17
             t21 = -c - t5 * t6 * t18

             e_0(ii) = e_0(ii) + ( t2 * t3 * t21 / 0.2e1_dp ) * sx

             t23 = t2 * my_rho_1_3
             t29 = 0.1e1_dp / t11
             t31 = t1 + t2 * x * t29
             t33 = 0.1e1_dp / t12
             t37 = 0.6e1_dp * t7 * t13 + 0.6e1_dp * t7 * x * t31 * t33
             t39 = t17 ** 2
             t40 = 0.1e1_dp / t39
             t43 = -0.2e1_dp * t5 * x * t18 + t5 * t6 * t37 * t40
             t44 = t43 * x

             e_rho(ii) = e_rho(ii)&
                  - ( 0.2e1_dp / 0.3e1_dp * t23 * t44 - 0.2e1_dp / 0.3e1_dp *&
                  t23 * t21 ) * sx
             e_ndrho(ii) = e_ndrho(ii)&
                  + ( t2 * t43 / 0.2e1_dp ) * sx
          END IF

       END DO

       !$omp end do

    CASE (-1)

       !$omp do

       DO ii=1,npoints

          my_rho = rho(ii)
          IF (my_rho>epsilon_rho) THEN
             my_rho_1_3 = rho_1_3(ii)
             t3 = my_rho_1_3 * my_rho
             x = norm_drho(ii)/t3
             t6 = x ** 2
             t10 = t2 * t6 + 0.1e1_dp
             t11 = SQRT(t10)
             t12 = t1 * x + t11
             t13 = LOG(t12)
             t17 = 0.1e1_dp + 0.6e1_dp * t7 * x * t13
             t18 = 0.1e1_dp / t17
             t21 = -c - t5 * t6 * t18

             t23 = t2 * my_rho_1_3
             t29 = 0.1e1_dp / t11
             t31 = t1 + t2 * x * t29
             t33 = 0.1e1_dp / t12
             t37 = 0.6e1_dp * t7 * t13 + 0.6e1_dp * t7 * x * t31 * t33
             t39 = t17 ** 2
             t40 = 0.1e1_dp / t39
             t43 = -0.2e1_dp * t5 * x * t18 + t5 * t6 * t37 * t40
             t44 = t43 * x

             e_rho(ii) = e_rho(ii)&
                  - ( 0.2e1_dp / 0.3e1_dp * t23 * t44 - 0.2e1_dp / 0.3e1_dp *&
                  t23 * t21 ) * sx
             e_ndrho(ii) = e_ndrho(ii)&
                  + ( t2 * t43 / 0.2e1_dp ) * sx
          END IF
       END DO

       !$omp end do

    CASE (2)

       !$omp do

       DO ii=1,npoints

          my_rho = rho(ii)
          IF (my_rho > epsilon_rho) THEN
             my_rho_1_3 = rho_1_3(ii)
             t3 = my_rho_1_3 * my_rho
             x = norm_drho(ii)/t3
             t6 = x ** 2
             t10 = t2 * t6 + 0.1e1_dp
             t11 = SQRT(t10)
             t12 = t1 * x + t11
             t13 = LOG(t12)
             t17 = 0.1e1_dp + 0.6e1_dp * t7 * x * t13
             t18 = 0.1e1_dp / t17
             t21 = -c - t5 * t6 * t18

             e_0(ii) = e_0(ii) + ( t2 * t3 * t21 / 0.2e1_dp ) * sx

             t23 = t2 * my_rho_1_3
             t29 = 0.1e1_dp / t11
             t31 = t1 + t2 * x * t29
             t33 = 0.1e1_dp / t12
             t37 = 0.6e1_dp * t7 * t13 + 0.6e1_dp * t7 * x * t31 * t33
             t39 = t17 ** 2
             t40 = 0.1e1_dp / t39
             t43 = -0.2e1_dp * t5 * x * t18 + t5 * t6 * t37 * t40
             t44 = t43 * x

             e_rho(ii) = e_rho(ii)&
                  - ( 0.2e1_dp / 0.3e1_dp * t23 * t44 - 0.2e1_dp / 0.3e1_dp *&
                  t23 * t21 ) * sx
             e_ndrho(ii) = e_ndrho(ii)&
                  + ( t2 * t43 / 0.2e1_dp ) * sx

             t49 = my_rho_1_3 ** 2
             t51 = t2 / t49
             t54 = x * t40
             t64 = 0.1e1_dp / t11 / t10
             t67 = t2 * t29 - 0.2e1_dp * t1 * t6 * t64
             t72 = t31 ** 2
             t74 = t12 ** 2
             t75 = 0.1e1_dp / t74
             t79 = 0.12e2_dp * t7 * t31 * t33 + 0.6e1_dp * t7 * x * t67 * t33 &
                  - 0.6e1_dp * t7 * x * t72 * t75
             t83 = t37 ** 2
             t86 = 0.1e1_dp / t39 / t17
             t90 = -0.2e1_dp * t5 * t18 + 0.4e1_dp * t5 * t54 * t37 + t5 * t6 *&
                  t79 * t40 - 0.2e1_dp * t5 * t6 * t83 * t86
             t91 = t90 * t6
             t98 = 0.1e1_dp / my_rho
             t99 = t2 * t98
             t100 = t90 * x
             t104 = 0.1e1_dp / t3

             e_rho_rho(ii) = e_rho_rho(ii)&
                  + ( 0.8e1_dp / 0.9e1_dp * t51 * t91 - 0.2e1_dp / 0.9e1_dp * &
                  t51 * t44 + 0.2e1_dp / 0.9e1_dp * t51 * t21 ) * sx
             e_ndrho_rho(ii) = e_ndrho_rho(ii)&
                  - ( 0.2e1_dp / 0.3e1_dp * t99 * t100 ) * sx
             e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii)&
                  + ( t2 * t90 * t104 / 0.2e1_dp ) * sx
          END IF
       END DO

       !$omp end do

    CASE (-2)

       !$omp do

       DO ii=1,npoints

          my_rho = rho(ii)
          IF (my_rho>epsilon_rho) THEN
             my_rho_1_3 = rho_1_3(ii)
             t3 = my_rho_1_3 * my_rho
             x = norm_drho(ii)/t3
             t6 = x ** 2
             t10 = t2 * t6 + 0.1e1_dp
             t11 = SQRT(t10)
             t12 = t1 * x + t11
             t13 = LOG(t12)
             t17 = 0.1e1_dp + 0.6e1_dp * t7 * x * t13
             t18 = 0.1e1_dp / t17
             t21 = -c - t5 * t6 * t18

             t23 = t2 * my_rho_1_3
             t29 = 0.1e1_dp / t11
             t31 = t1 + t2 * x * t29
             t33 = 0.1e1_dp / t12
             t37 = 0.6e1_dp * t7 * t13 + 0.6e1_dp * t7 * x * t31 * t33
             t39 = t17 ** 2
             t40 = 0.1e1_dp / t39
             t43 = -0.2e1_dp * t5 * x * t18 + t5 * t6 * t37 * t40
             t44 = t43 * x

             t49 = my_rho_1_3 ** 2
             t51 = t2 / t49
             t54 = x * t40
             t64 = 0.1e1_dp / t11 / t10
             t67 = t2 * t29 - 0.2e1_dp * t1 * t6 * t64
             t72 = t31 ** 2
             t74 = t12 ** 2
             t75 = 0.1e1_dp / t74
             t79 = 0.12e2_dp * t7 * t31 * t33 + 0.6e1_dp * t7 * x * t67 * t33 &
                  - 0.6e1_dp * t7 * x * t72 * t75
             t83 = t37 ** 2
             t86 = 0.1e1_dp / t39 / t17
             t90 = -0.2e1_dp * t5 * t18 + 0.4e1_dp * t5 * t54 * t37 + t5 * t6 *&
                  t79 * t40 - 0.2e1_dp * t5 * t6 * t83 * t86
             t91 = t90 * t6
             t98 = 0.1e1_dp / my_rho
             t99 = t2 * t98
             t100 = t90 * x
             t104 = 0.1e1_dp / t3

             e_rho_rho(ii) = e_rho_rho(ii)&
                  + ( 0.8e1_dp / 0.9e1_dp * t51 * t91 - 0.2e1_dp / 0.9e1_dp * &
                  t51 * t44 + 0.2e1_dp / 0.9e1_dp * t51 * t21 ) * sx
             e_ndrho_rho(ii) = e_ndrho_rho(ii)&
                  - ( 0.2e1_dp / 0.3e1_dp * t99 * t100 ) * sx
             e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii)&
                  + ( t2 * t90 * t104 / 0.2e1_dp ) * sx
          END IF
       END DO

       !$omp end do

    CASE default

       !$omp do

       DO ii=1,npoints

          my_rho = rho(ii)
          IF (my_rho>epsilon_rho) THEN
             my_rho_1_3 = rho_1_3(ii)
             t3 = my_rho_1_3 * my_rho
             x = norm_drho(ii)/t3
             t6 = x ** 2
             t10 = t2 * t6 + 0.1e1_dp
             t11 = SQRT(t10)
             t12 = t1 * x + t11
             t13 = LOG(t12)
             t17 = 0.1e1_dp + 0.6e1_dp * t7 * x * t13
             t18 = 0.1e1_dp / t17
             t21 = -c - t5 * t6 * t18

             IF (grad_deriv>=0) THEN
                e_0(ii) = e_0(ii) + ( t2 * t3 * t21 / 0.2e1_dp ) * sx
             END IF

             t23 = t2 * my_rho_1_3
             t29 = 0.1e1_dp / t11
             t31 = t1 + t2 * x * t29
             t33 = 0.1e1_dp / t12
             t37 = 0.6e1_dp * t7 * t13 + 0.6e1_dp * t7 * x * t31 * t33
             t39 = t17 ** 2
             t40 = 0.1e1_dp / t39
             t43 = -0.2e1_dp * t5 * x * t18 + t5 * t6 * t37 * t40
             t44 = t43 * x

             IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
                e_rho(ii) = e_rho(ii)&
                     - ( 0.2e1_dp / 0.3e1_dp * t23 * t44 - 0.2e1_dp / 0.3e1_dp *&
                     t23 * t21 ) * sx
                e_ndrho(ii) = e_ndrho(ii)&
                     + ( t2 * t43 / 0.2e1_dp ) * sx
             END IF

             t49 = my_rho_1_3 ** 2
             t51 = t2 / t49
             t54 = x * t40
             t64 = 0.1e1_dp / t11 / t10
             t67 = t2 * t29 - 0.2e1_dp * t1 * t6 * t64
             t72 = t31 ** 2
             t74 = t12 ** 2
             t75 = 0.1e1_dp / t74
             t79 = 0.12e2_dp * t7 * t31 * t33 + 0.6e1_dp * t7 * x * t67 * t33 &
                  - 0.6e1_dp * t7 * x * t72 * t75
             t83 = t37 ** 2
             t86 = 0.1e1_dp / t39 / t17
             t90 = -0.2e1_dp * t5 * t18 + 0.4e1_dp * t5 * t54 * t37 + t5 * t6 *&
                  t79 * t40 - 0.2e1_dp * t5 * t6 * t83 * t86
             t91 = t90 * t6
             t98 = 0.1e1_dp / my_rho
             t99 = t2 * t98
             t100 = t90 * x
             t104 = 0.1e1_dp / t3

             IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
                e_rho_rho(ii) = e_rho_rho(ii)&
                     + ( 0.8e1_dp / 0.9e1_dp * t51 * t91 - 0.2e1_dp / 0.9e1_dp * &
                     t51 * t44 + 0.2e1_dp / 0.9e1_dp * t51 * t21 ) * sx
                e_ndrho_rho(ii) = e_ndrho_rho(ii)&
                     - ( 0.2e1_dp / 0.3e1_dp * t99 * t100 ) * sx
                e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii)&
                     + ( t2 * t90 * t104 / 0.2e1_dp ) * sx
             END IF

             t126 = t10 ** 2
             t159 = t39 ** 2
             t164 = 0.6e1_dp * t5 * t40 * t37 - 0.12e2_dp * t5 * x * t86 * t83 &
                  + 0.6e1_dp * t5 * t54 * t79 + t5 * t6 * (0.18e2_dp * t7 * t67 &
                  * t33 - 0.18e2_dp * t7 * t72 * t75 + 0.6e1_dp * t7 * x * &
                  (-0.6e1_dp * t1 * t64 * x + 0.12e2_dp * t6 * x / t11 / t126)&
                  * t33 - 0.18e2_dp * t7 * x * t67 * t75 * t31 + 0.12e2_dp * t7 &
                  * x * t72 * t31 / t74 / t12) * t40 - 0.6e1_dp * t5 * t6 * t79 &
                  * t86 * t37 + 0.6e1_dp * t5 * t6 * t83 * t37 / t159
             t170 = 0.8e1_dp / 0.9e1_dp * t51 * t164 * t6 + 0.14e2_dp / 0.9e1_dp&
                  * t51 * t100
             t176 = t2 / t49 / my_rho
             t189 = my_rho ** 2

             IF (grad_deriv==-3.or.grad_deriv>=3) THEN
                e_rho_rho_rho(ii) = e_rho_rho_rho(ii) &
                     - ( 0.4e1_dp / 0.3e1_dp * t170 * x * t98 + 0.16e2_dp / 0.27e2_dp &
                     * t176 * t91 - 0.4e1_dp / 0.27e2_dp * t176 * t44 + &
                     0.4e1_dp / 0.27e2_dp * t176 * t21 ) * sx
                e_ndrho_rho_rho(ii) = e_ndrho_rho_rho(ii)&
                     + ( t170 * t104 ) * sx
                e_ndrho_ndrho_rho(ii) = e_ndrho_ndrho_rho(ii)&
                     + ( (-0.2e1_dp / 0.3e1_dp * t99 * t164 * x - &
                     0.2e1_dp / 0.3e1_dp * t99 * t90) * t104 ) * sx
                e_ndrho_ndrho_ndrho(ii) = e_ndrho_ndrho_ndrho(ii)&
                     + ( t2 * t164 / t49 / t189 / 0.2e1_dp ) * sx
             END IF
          END IF
       END DO

       !$omp end do

    END SELECT

  END SUBROUTINE xb88_lda_calc

! *****************************************************************************
!> \brief evaluates the becke 88 exchange functional for lsd
!> \param rho_set ...
!> \param deriv_set ...
!> \param grad_deriv ...
!> \param xb88_params ...
!> \par History
!>      11.2003 created [fawzi]
!>      01.2007 added scaling [Manuel Guidon]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE xb88_lsd_eval(rho_set,deriv_set,grad_deriv,xb88_params)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: xb88_params

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_rho, sx
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0
    TYPE(cp_3d_r_p_type), DIMENSION(2) :: e_ndrho, e_ndrho_ndrho, &
      e_ndrho_ndrho_ndrho, e_ndrho_ndrho_rho, e_ndrho_rho, e_ndrho_rho_rho, &
      e_rho, e_rho_rho, e_rho_rho_rho, norm_drho, rho, rho_1_3
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)

    CALL cite_reference(Becke1988)

    failure=.FALSE.
    NULLIFY(deriv, bo)
    DO i=1,2
       NULLIFY(norm_drho(i)%array, rho(i)%array, rho_1_3(i)%array)
    END DO

    CALL section_vals_val_get(xb88_params,"scale_x",r_val=sx)
    CPASSERT(ASSOCIATED(rho_set))
    CPASSERT(rho_set%ref_count>0)
    CPASSERT(ASSOCIATED(deriv_set))
    CPASSERT(deriv_set%ref_count>0)
    CALL xc_rho_set_get(rho_set,rhoa_1_3=rho_1_3(1)%array,&
         rhob_1_3=rho_1_3(2)%array,rhoa=rho(1)%array,&
         rhob=rho(2)%array,norm_drhoa=norm_drho(1)%array, &
         norm_drhob=norm_drho(2)%array,rho_cutoff=epsilon_rho,&
         local_bounds=bo)
    npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

    ! meaningful default for the arrays we don't need: let us make compiler
    ! and debugger happy...
    IF (cp_debug) THEN
       ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)))
    ELSE
       dummy=> rho(1)%array
    END IF

    e_0=>dummy
    DO i=1,2
       e_rho(i)%array => dummy
       e_ndrho(i)%array => dummy
       e_rho_rho(i)%array => dummy
       e_ndrho_rho(i)%array => dummy
       e_ndrho_ndrho(i)%array => dummy
       e_rho_rho_rho(i)%array => dummy
       e_ndrho_rho_rho(i)%array => dummy
       e_ndrho_ndrho_rho(i)%array => dummy
       e_ndrho_ndrho_ndrho(i)%array => dummy
    END DO

    IF (grad_deriv>=0) THEN
       deriv => xc_dset_get_derivative(deriv_set,"",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv, deriv_data=e_0)
    END IF
    IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
       deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho(1)%array)
       deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho(2)%array)
       deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho(1)%array)
       deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho(2)%array)
    END IF
    IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
       deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho_rho(1)%array)
       deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho_rho(2)%array)
       deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhoa)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho(1)%array)
       deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhob)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho(2)%array)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(norm_drhoa)(norm_drhoa)", allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho(1)%array)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(norm_drhob)(norm_drhob)", allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho(2)%array)
    END IF
    IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
       deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)(rhoa)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho(1)%array)
       deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)(rhob)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho(2)%array)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(norm_drhoa)(rhoa)(rhoa)",allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho_rho(1)%array)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(norm_drhob)(rhob)(rhob)",allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho_rho(2)%array)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(norm_drhoa)(norm_drhoa)(rhoa)",allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_rho(1)%array)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(norm_drhob)(norm_drhob)(rhob)",allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_rho(2)%array)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(norm_drhoa)(norm_drhoa)(norm_drhoa)", allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho(1)%array)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(norm_drhob)(norm_drhob)(norm_drhob)", allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho(2)%array)
    END IF
    IF (grad_deriv>3.OR.grad_deriv<-3) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="derivatives bigger than 3 not implemented", &
            error_level=cp_failure_level)
    END IF

    DO ispin=1,2

       !$omp parallel default(none) &
       !$omp          shared(rho, ispin, rho_1_3, norm_drho, e_0) &
       !$omp          shared(e_rho, e_ndrho, e_rho_rho, e_ndrho_rho) &
       !$omp          shared(e_ndrho_ndrho, e_rho_rho_rho) &
       !$omp          shared(e_ndrho_rho_rho, e_ndrho_ndrho_rho) &
       !$omp          shared(e_ndrho_ndrho_ndrho, grad_deriv, npoints) &
       !$omp          shared(epsilon_rho, sx)

       CALL xb88_lsd_calc(&
            rho_spin=rho(ispin)%array,&
            rho_1_3_spin=rho_1_3(ispin)%array,&
            norm_drho_spin=norm_drho(ispin)%array,&
            e_0=e_0,&
            e_rho_spin=e_rho(ispin)%array,&
            e_ndrho_spin=e_ndrho(ispin)%array,&
            e_rho_rho_spin=e_rho_rho(ispin)%array,&
            e_ndrho_rho_spin=e_ndrho_rho(ispin)%array,&
            e_ndrho_ndrho_spin=e_ndrho_ndrho(ispin)%array,&
            e_rho_rho_rho_spin=e_rho_rho_rho(ispin)%array,&
            e_ndrho_rho_rho_spin=e_ndrho_rho_rho(ispin)%array,&
            e_ndrho_ndrho_rho_spin=e_ndrho_ndrho_rho(ispin)%array,&
            e_ndrho_ndrho_ndrho_spin=e_ndrho_ndrho_ndrho(ispin)%array,&
            grad_deriv=grad_deriv, npoints=npoints,&
            epsilon_rho=epsilon_rho, sx=sx)

       !$omp end parallel

    END DO

    IF (cp_debug) THEN
       DEALLOCATE(dummy)
    ELSE
       NULLIFY(dummy)
    END IF
    CALL timestop(handle)

  END SUBROUTINE xb88_lsd_eval

! *****************************************************************************
!> \brief low level calculation of the becke 88 exchange functional for lsd
!> \param rho_spin alpha or beta spin density
!> \param rho_1_3_spin rho_spin**(1./3.)
!> \param norm_drho_spin || grad rho_spin ||
!> \param e_0 adds to it the local value of the functional
!> \param e_rho_spin e_*_spin: derivative of the functional wrt. to the variables
!>        named where the * is. Everything wrt. to the spin of the arguments.
!> \param e_ndrho_spin ...
!> \param e_rho_rho_spin ...
!> \param e_ndrho_rho_spin ...
!> \param e_ndrho_ndrho_spin ...
!> \param e_rho_rho_rho_spin ...
!> \param e_ndrho_rho_rho_spin ...
!> \param e_ndrho_ndrho_rho_spin ...
!> \param e_ndrho_ndrho_ndrho_spin ...
!> \param grad_deriv ...
!> \param npoints ...
!> \param epsilon_rho ...
!> \param sx scaling-parameter for exchange
!> \par History
!>      01.2004 created [fawzi]
!>      01.2007 added scaling [Manuel Guidon]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE xb88_lsd_calc(rho_spin, rho_1_3_spin, norm_drho_spin,e_0,&
       e_rho_spin,e_ndrho_spin,e_rho_rho_spin,e_ndrho_rho_spin,&
       e_ndrho_ndrho_spin, e_rho_rho_rho_spin, e_ndrho_rho_rho_spin,&
       e_ndrho_ndrho_rho_spin,&
       e_ndrho_ndrho_ndrho_spin,grad_deriv,npoints,epsilon_rho,sx)
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rho_spin, rho_1_3_spin, &
                                                norm_drho_spin
    REAL(kind=dp), DIMENSION(*), INTENT(inout) :: e_0, e_rho_spin, &
      e_ndrho_spin, e_rho_rho_spin, e_ndrho_rho_spin, e_ndrho_ndrho_spin, &
      e_rho_rho_rho_spin, e_ndrho_rho_rho_spin, e_ndrho_ndrho_rho_spin, &
      e_ndrho_ndrho_ndrho_spin
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx

    CHARACTER(len=*), PARAMETER :: routineN = 'xb88_lsd_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii
    LOGICAL                                  :: failure
    REAL(kind=dp) :: c, epsilon_rho43, my_epsilon_rho, my_rho, t1, t103, t11, &
      t12, t127, t133, t138, t14, t151, t17, t18, t2, t20, t22, t23, t27, &
      t28, t3, t30, t35, t36, t4, t42, t43, t44, t5, t51, t53, t57, t58, t59, &
      t6, t63, t64, t66, t67, t7, t75, t76, t79, t8, t87, x

    failure=.FALSE.
    c= 1.5_dp*(0.75_dp/Pi)**(1.0_dp/3.0_dp)
    my_epsilon_rho=0.5_dp*epsilon_rho
    epsilon_rho43=my_epsilon_rho**(4._dp/3._dp)


    SELECT CASE(grad_deriv)
    CASE(0)

       !$omp do

       DO ii=1,npoints
          my_rho=rho_spin(ii)
          IF (my_rho>my_epsilon_rho) THEN
             t1 = rho_1_3_spin(ii) * my_rho
             x = norm_drho_spin(ii)/t1
             t2 = x ** 2
             t3 = beta * t2
             t4 = beta * x
             t5 = t2 + 0.1e1_dp
             t6 = SQRT(t5)
             t7 = x + t6
             t8 = LOG(t7)
             t11 = 0.1e1_dp + 0.6e1_dp * t4 * t8
             t12 = 0.1e1_dp / t11
             t14 = -c - t3 * t12

             e_0(ii) = e_0(ii)&
                  + ( t1 * t14 ) * sx
          END IF
       END DO

       !$omp end do

    CASE (1)

       !$omp do

       DO ii=1,npoints
          my_rho=rho_spin(ii)
          IF (my_rho>my_epsilon_rho) THEN
             t1 = rho_1_3_spin(ii) * my_rho
             x = norm_drho_spin(ii)/t1
             t2 = x ** 2
             t3 = beta * t2
             t4 = beta * x
             t5 = t2 + 0.1e1_dp
             t6 = SQRT(t5)
             t7 = x + t6
             t8 = LOG(t7)
             t11 = 0.1e1_dp + 0.6e1_dp * t4 * t8
             t12 = 0.1e1_dp / t11
             t14 = -c - t3 * t12

             e_0(ii) = e_0(ii)&
                  + ( t1 * t14 ) * sx

             t17 = t11 ** 2
             t18 = 0.1e1_dp / t17
             t20 = 0.1e1_dp / t6
             t22 = 0.1e1_dp + t20 * x
             t23 = 0.1e1_dp / t7
             t27 = 0.6e1_dp * beta * t8 + 0.6e1_dp * t4 * t22 * t23
             t28 = t18 * t27
             t30 = -0.2e1_dp * t4 * t12 + t3 * t28

             e_rho_spin(ii) = e_rho_spin(ii)&
                  - ( 0.4e1_dp / 0.3e1_dp * rho_1_3_spin(ii) * t30 * x - &
                  0.4e1_dp / 0.3e1_dp * rho_1_3_spin(ii) * t14 ) * sx
             e_ndrho_spin(ii) = e_ndrho_spin(ii)&
                  + ( t30 ) * sx
          END IF
       END DO

       !$omp end do

    CASE default

       !$omp do

       DO ii=1,npoints
          my_rho=rho_spin(ii)
          IF (my_rho>my_epsilon_rho) THEN
             t1 = rho_1_3_spin(ii) * my_rho
             x = norm_drho_spin(ii)/t1
             t2 = x ** 2
             t3 = beta * t2
             t4 = beta * x
             t5 = t2 + 0.1e1_dp
             t6 = SQRT(t5)
             t7 = x + t6
             t8 = LOG(t7)
             t11 = 0.1e1_dp + 0.6e1_dp * t4 * t8
             t12 = 0.1e1_dp / t11
             t14 = -c - t3 * t12

             IF (grad_deriv>=0) THEN
                e_0(ii) = e_0(ii)&
                     + ( t1 * t14 ) * sx
             END IF

             t17 = t11 ** 2
             t18 = 0.1e1_dp / t17
             t20 = 0.1e1_dp / t6
             t22 = 0.1e1_dp + t20 * x
             t23 = 0.1e1_dp / t7
             t27 = 0.6e1_dp * beta * t8 + 0.6e1_dp * t4 * t22 * t23
             t28 = t18 * t27
             t30 = -0.2e1_dp * t4 * t12 + t3 * t28

             IF (grad_deriv==-1.or.grad_deriv>=1) THEN
                e_rho_spin(ii) = e_rho_spin(ii)&
                     - ( 0.4e1_dp / 0.3e1_dp * rho_1_3_spin(ii) * t30 * x - &
                     0.4e1_dp / 0.3e1_dp * rho_1_3_spin(ii) * t14 ) * sx
                e_ndrho_spin(ii) = e_ndrho_spin(ii)&
                     + ( t30 ) * sx
             END IF

             t35 = rho_1_3_spin(ii) ** 2
             t36 = 0.1e1_dp / t35
             t42 = 0.1e1_dp / t17 / t11
             t43 = t27 ** 2
             t44 = t42 * t43
             t51 = 0.1e1_dp / t6 / t5
             t53 = -t51 * t2 + t20
             t57 = t22 ** 2
             t58 = t7 ** 2
             t59 = 0.1e1_dp / t58
             t63 = 0.12e2_dp * beta * t22 * t23 + 0.6e1_dp * t4 * t53 * t23 - 0.6e1_dp * t4 * t57 * t59
             t64 = t18 * t63
             t66 = -0.2e1_dp * beta * t12 + 0.4e1_dp * t4 * t28 - 0.2e1_dp * t3 * t44 + t3 * t64
             t67 = t36 * t66
             t75 = 0.1e1_dp / my_rho
             t76 = t75 * t66
             t79 = 0.1e1_dp / t1

             IF (grad_deriv==-2.or.grad_deriv>=2) THEN
                e_rho_rho_spin(ii) = e_rho_rho_spin(ii)&
                     + ( 0.16e2_dp / 0.9e1_dp * t67 * t2 - 0.4e1_dp / 0.9e1_dp &
                     * t36 * t30 * x + 0.4e1_dp / 0.9e1_dp * t36 * t14 ) * sx
                e_ndrho_rho_spin(ii) = e_ndrho_rho_spin(ii)&
                     - ( 0.4e1_dp / 0.3e1_dp * t76 * x ) * sx
                e_ndrho_ndrho_spin(ii) = e_ndrho_ndrho_spin(ii)+&
                     ( t66 * t79 ) * sx
             END IF

             t87 = t17 ** 2
             t103 = t5 ** 2
             t127 = 0.6e1_dp * beta * t18 * t27 - 0.12e2_dp * t4 * t44 + 0.6e1_dp&
                  * t4 * t64 + 0.6e1_dp * t3 / t87 * t43 * t27 - 0.6e1_dp * t3 &
                  * t42 * t27 * t63 + t3 * t18 * (0.18e2_dp * beta * t53 * t23 &
                  - 0.18e2_dp * beta * t57 * t59 + 0.6e1_dp * t4 * (0.3e1_dp / &
                  t6 / t103 * t2 * x - 0.3e1_dp * t51 * x) * t23 - 0.18e2_dp * &
                  t4 * t53 * t59 * t22 + 0.12e2_dp * t4 * t57 * t22 / t58 / t7)
             t133 = 0.16e2_dp / 0.9e1_dp * t36 * t127 * t2 + 0.28e2_dp / 0.9e1_dp &
                  * t67 * x
             t138 = 0.1e1_dp / t35 / my_rho
             t151 = my_rho ** 2

             IF (grad_deriv==-3.or.grad_deriv>=3) THEN
                e_rho_rho_rho_spin(ii) = e_rho_rho_rho_spin(ii)&
                     - ( 0.4e1_dp / 0.3e1_dp * t133 * x * t75 + 0.32e2_dp / 0.27e2_dp &
                     * t138 * t66 * t2 - 0.8e1_dp / 0.27e2_dp * t138 * t30 * &
                     x + 0.8e1_dp / 0.27e2_dp * t138 * t14 ) * sx
                e_ndrho_rho_rho_spin(ii) = e_ndrho_rho_rho_spin(ii)&
                     + ( t133 * t79 ) * sx
                e_ndrho_ndrho_rho_spin(ii) = e_ndrho_ndrho_rho_spin(ii)&
                     + ( (-0.4e1_dp / 0.3e1_dp * t75 * t127 * x - &
                     0.4e1_dp / 0.3e1_dp * t76) * t79 ) * sx
                e_ndrho_ndrho_ndrho_spin(ii) = e_ndrho_ndrho_ndrho_spin(ii)&
                     + ( t127 / t35 / t151 ) * sx
             END IF
          END IF
       END DO

       !$omp end do

    END SELECT

  END SUBROUTINE xb88_lsd_calc

END MODULE xc_xbecke88
