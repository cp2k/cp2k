!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief calculates the lyp correlation functional
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE xc_lyp
  USE bibliography,                    ONLY: Lee1988,&
                                             cite_reference
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
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_lyp'
  REAL(kind=dp), PARAMETER, PRIVATE :: a=0.04918_dp, b=0.132_dp,&
       c=0.2533_dp,d=0.349_dp

  PUBLIC :: lyp_lda_info, lyp_lsd_info, lyp_lda_eval, lyp_lsd_eval
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
  SUBROUTINE lyp_lda_info(reference,shortform, needs, max_deriv)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_lda_info', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
       reference = "C. Lee, W. Yang, R.G. Parr, Phys. Rev. B, 37, 785 (1988) {LDA version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "Lee-Yang-Parr correlation energy functional (LDA)"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho=.TRUE.
       needs%rho_1_3=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=3

  END SUBROUTINE lyp_lda_info

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
SUBROUTINE lyp_lsd_info(reference,shortform, needs, max_deriv)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_lsd_info', &
      routineP = moduleN//':'//routineN

  IF ( PRESENT ( reference ) ) THEN
     reference = "C. Lee, W. Yang, R.G. Parr, Phys. Rev. B, 37, 785 (1988) {LSD version}"
  END IF
  IF ( PRESENT ( shortform ) ) THEN
     shortform = "Lee-Yang-Parr correlation energy functional (LSD)"
  END IF
  IF (PRESENT(needs)) THEN
     needs%rho_spin=.TRUE.
     needs%norm_drho_spin=.TRUE.
     needs%norm_drho=.TRUE.
  END IF
  IF (PRESENT(max_deriv)) max_deriv=2

END SUBROUTINE lyp_lsd_info

! *****************************************************************************
!> \brief evaluates the lyp correlation functional for lda
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param lyp_params input parameters (scaling)
!> \par History
!>      11.2003 created [fawzi]
!>      01.2007 added scaling [Manuel Guidon]
!> \author fawzi
! *****************************************************************************
SUBROUTINE lyp_lda_eval(rho_set,deriv_set,grad_deriv,lyp_params)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: lyp_params

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_norm_drho, &
                                                epsilon_rho, sc
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndrho, &
      e_ndrho_ndrho, e_ndrho_ndrho_rho, e_ndrho_rho, e_ndrho_rho_rho, e_rho, &
      e_rho_rho, e_rho_rho_rho, norm_drho, rho, rho_1_3
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)

  failure=.FALSE.
  NULLIFY(bo)

  CALL section_vals_val_get(lyp_params,"scale_c",r_val=sc)
  CALL cite_reference(Lee1988)

  CPASSERT(ASSOCIATED(rho_set))
  CPASSERT(rho_set%ref_count>0)
  CPASSERT(ASSOCIATED(deriv_set))
  CPASSERT(deriv_set%ref_count>0)
  CALL xc_rho_set_get(rho_set,rho_1_3=rho_1_3,rho=rho,&
       norm_drho=norm_drho,local_bounds=bo,rho_cutoff=epsilon_rho,&
       drho_cutoff=epsilon_norm_drho)
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
!FM        deriv => xc_dset_get_derivative(deriv_set,&
!FM             "(norm_drho)(norm_drho)(norm_drho)", allocate_deriv=.TRUE.,&
!FM             
!FM        call xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho)
  END IF
  IF (grad_deriv>3.OR.grad_deriv<-3) THEN
     CALL cp_unimplemented_error(fromWhere=routineP, &
          message="derivatives bigger than 3 not implemented", &
          error_level=cp_failure_level)
  END IF

  !$omp parallel default(none) &
  !$omp          shared(rho, rho_1_3, norm_drho, e_0, e_rho, e_ndrho) &
  !$omp          shared(e_rho_rho, e_ndrho_rho, e_ndrho_ndrho) &
  !$omp          shared(e_rho_rho_rho, e_ndrho_rho_rho) &
  !$omp          shared(e_ndrho_ndrho_rho, grad_deriv, npoints) &
  !$omp          shared(epsilon_rho,  sc)

  CALL lyp_lda_calc(rho=rho, rho_1_3=rho_1_3, norm_drho=norm_drho,&
       e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,e_rho_rho=e_rho_rho,&
       e_ndrho_rho=e_ndrho_rho, e_ndrho_ndrho=e_ndrho_ndrho, &
       e_rho_rho_rho=e_rho_rho_rho, e_ndrho_rho_rho=e_ndrho_rho_rho,&
       e_ndrho_ndrho_rho=e_ndrho_ndrho_rho,&
       grad_deriv=grad_deriv,&
       npoints=npoints,epsilon_rho=epsilon_rho,sc=sc)

  !$omp end parallel

  IF (cp_debug) THEN
     DEALLOCATE(dummy)
  ELSE
     NULLIFY(dummy)
  END IF
  CALL timestop(handle)
END SUBROUTINE lyp_lda_eval

! *****************************************************************************
!> \brief evaluates the lyp correlation functional for lda
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
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param npoints ...
!> \param epsilon_rho ...
!> \param sc scaling-parameter for correlation
!> \par History
!>      11.2003 created [fawzi]
!>      01.2007 added scaling [Manuel Guidon]
!> \author fawzi
! *****************************************************************************
SUBROUTINE lyp_lda_calc(rho, rho_1_3, norm_drho,&
     e_0,e_rho,e_ndrho,e_rho_rho,e_ndrho_rho,&
     e_ndrho_ndrho, e_rho_rho_rho, e_ndrho_rho_rho, e_ndrho_ndrho_rho,&
     grad_deriv,npoints,epsilon_rho,sc)
    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), DIMENSION(1:npoints), INTENT(inout) :: e_ndrho_ndrho_rho, &
      e_ndrho_rho_rho, e_rho_rho_rho, e_ndrho_ndrho, e_ndrho_rho, e_rho_rho, &
      e_ndrho, e_rho, e_0
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(in)                             :: norm_drho, rho_1_3, rho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sc

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_lda_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii
    REAL(kind=dp) :: cf, epsilon_rho13, my_ndrho, my_rho, my_rho_1_3, t1, &
      t102, t103, t105, t106, t11, t112, t123, t124, t127, t128, t13, t133, &
      t134, t14, t161, t165, t166, t173, t176, t184, t185, t189, t192, t193, &
      t196, t199, t2, t200, t201, t202, t203, t215, t22, t227, t228, t231, &
      t245, t246, t248, t26, t265, t278, t279, t3, t332, t37, t38, t39, t4, &
      t40, t41, t44, t45, t48, t5, t52, t6, t61, t62, t69, t7, t70, t77, t78, &
      t80, t87, t88, t89, t93, t94, t95, t98, t99

  epsilon_rho13=epsilon_rho**(1.0_dp/3.0_dp)
  cf=0.3_dp*(3._dp*pi*pi)**(2._dp/3._dp)
  SELECT CASE(grad_deriv)
  CASE (1)
     !$omp do
     DO ii=1,npoints
        my_rho = rho(ii)
        IF (my_rho>epsilon_rho) THEN
           my_rho_1_3 = rho_1_3(ii)
           my_ndrho = norm_drho(ii)
           t1 = my_rho_1_3 ** 2
           t2 = t1 * my_rho
           t3 = 0.1e1_dp / t2
           t4 = a * t3
           t5 = my_rho ** 2
           t6 = t5 * my_rho
           t7 = my_rho_1_3 * t6
           t11 = 0.1e1_dp / my_rho_1_3
           t13 = EXP(-c * t11)
           t14 = b * t13
           t22 = my_ndrho ** 2
           t26 = my_rho_1_3 * t22
           t37 = -0.72e2_dp * t7 - 0.72e2_dp * t6 * d - 0.72e2_dp * t14 *&
                  t7 * cf - 0.72e2_dp * t14 * t6 * cf * d + 0.3e1_dp * t14&
                & * t1 * t22 + 0.10e2_dp * t14 * t26 * d + 0.7e1_dp * t14 &
                &* t26 * c + 0.7e1_dp * t14 * t22 * c * d
           t38 = my_rho_1_3 + d
           t39 = t38 ** 2
           t40 = 0.1e1_dp / t39
           t41 = t37 * t40

           e_0(ii) = e_0(ii)&
                + ( t4 * t41 / 0.72e2_dp ) * sc
           t44 = 0.1e1_dp / t1 / t5
           t45 = a * t44
           t48 = my_rho_1_3 * t5
           t52 = b * c
           t61 = t13 * cf
           t62 = t61 * d
           t69 = 0.1e1_dp / t1
           t70 = t69 * t13
           t77 = 0.1e1_dp / my_rho
           t78 = t52 * t77
           t80 = t13 * t22 * d
           t87 = c ** 2
           t88 = b * t87
           t89 = t77 * t13
           t93 = my_rho_1_3 * my_rho
           t94 = 0.1e1_dp / t93
           t95 = t88 * t94
           t98 = -0.240e3_dp * t48 - 0.216e3_dp * t5 * d - 0.24e2_dp * t52&
                & * t5 * t13 * cf - 0.240e3_dp * t14 * t48 * cf -&
                  0.24e2_dp * t52 * t2 * t62 - 0.216e3_dp * t14 * t5 * cf &
                &* d + 0.10e2_dp / 0.3e1_dp * t52 * t70 * t22 + 0.2e1_dp *&
                  t14 * t11 * t22 + 0.10e2_dp / 0.3e1_dp * t78 * t80 +&
                  0.10e2_dp / 0.3e1_dp * t14 * t69 * t22 * d + 0.7e1_dp /&
                  0.3e1_dp * t88 * t89 * t22 + 0.7e1_dp / 0.3e1_dp * t95 *&
                  t80
           t99 = t98 * t40
           t102 = 0.1e1_dp / t48
           t103 = a * t102
           t105 = 0.1e1_dp / t39 / t38
           t106 = t37 * t105

           e_rho(ii) = e_rho(ii) &
                 - ( 0.5e1_dp / 0.216e3_dp * t45 * t41 - t4 * t99 / 0.72e2_dp&
                & + t103 * t106 / 0.108e3_dp) * sc

           t112 = my_rho_1_3 * my_ndrho
           t123 = 0.6e1_dp * t14 * t1 * my_ndrho + 0.20e2_dp * t14 * t112 &
                &* d + 0.14e2_dp * t14 * t112 * c + 0.14e2_dp * t14 *&
                  my_ndrho * c * d
           t124 = t123 * t40

           e_ndrho(ii) = e_ndrho(ii)&
                + ( t4 * t124 / 0.72e2_dp ) * sc
        END IF
     END DO
     !$end do
  CASE default
     !$omp do
     DO ii=1,npoints
        my_rho = rho(ii)
        IF (my_rho>epsilon_rho) THEN
           my_rho_1_3 = rho_1_3(ii)
           my_ndrho = norm_drho(ii)
           t1 = my_rho_1_3 ** 2
           t2 = t1 * my_rho
           t3 = 0.1e1_dp / t2
           t4 = a * t3
           t5 = my_rho ** 2
           t6 = t5 * my_rho
           t7 = my_rho_1_3 * t6
           t11 = 0.1e1_dp / my_rho_1_3
           t13 = EXP(-c * t11)
           t14 = b * t13
           t22 = my_ndrho ** 2
           t26 = my_rho_1_3 * t22
           t37 = -0.72e2_dp * t7 - 0.72e2_dp * t6 * d - 0.72e2_dp * t14 *&
                  t7 * cf - 0.72e2_dp * t14 * t6 * cf * d + 0.3e1_dp * t14&
                & * t1 * t22 + 0.10e2_dp * t14 * t26 * d + 0.7e1_dp * t14 &
                &* t26 * c + 0.7e1_dp * t14 * t22 * c * d
           t38 = my_rho_1_3 + d
           t39 = t38 ** 2
           t40 = 0.1e1_dp / t39
           t41 = t37 * t40

           IF (grad_deriv>=0) THEN
              e_0(ii) = e_0(ii)&
                   + ( t4 * t41 / 0.72e2_dp ) * sc
           END IF

           t44 = 0.1e1_dp / t1 / t5
           t45 = a * t44
           t48 = my_rho_1_3 * t5
           t52 = b * c
           t61 = t13 * cf
           t62 = t61 * d
           t69 = 0.1e1_dp / t1
           t70 = t69 * t13
           t77 = 0.1e1_dp / my_rho
           t78 = t52 * t77
           t80 = t13 * t22 * d
           t87 = c ** 2
           t88 = b * t87
           t89 = t77 * t13
           t93 = my_rho_1_3 * my_rho
           t94 = 0.1e1_dp / t93
           t95 = t88 * t94
           t98 = -0.240e3_dp * t48 - 0.216e3_dp * t5 * d - 0.24e2_dp * t52&
                & * t5 * t13 * cf - 0.240e3_dp * t14 * t48 * cf -&
                  0.24e2_dp * t52 * t2 * t62 - 0.216e3_dp * t14 * t5 * cf &
                &* d + 0.10e2_dp / 0.3e1_dp * t52 * t70 * t22 + 0.2e1_dp *&
                  t14 * t11 * t22 + 0.10e2_dp / 0.3e1_dp * t78 * t80 +&
                  0.10e2_dp / 0.3e1_dp * t14 * t69 * t22 * d + 0.7e1_dp /&
                  0.3e1_dp * t88 * t89 * t22 + 0.7e1_dp / 0.3e1_dp * t95 *&
                  t80
           t99 = t98 * t40
           t102 = 0.1e1_dp / t48
           t103 = a * t102
           t105 = 0.1e1_dp / t39 / t38
           t106 = t37 * t105
           t112 = my_rho_1_3 * my_ndrho
           t123 = 0.6e1_dp * t14 * t1 * my_ndrho + 0.20e2_dp * t14 * t112 &
                &* d + 0.14e2_dp * t14 * t112 * c + 0.14e2_dp * t14 *&
                  my_ndrho * c * d
           t124 = t123 * t40

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_rho(ii) = e_rho(ii) &
                   - ( 0.5e1_dp / 0.216e3_dp * t45 * t41 - t4 * t99 / 0.72e2_dp&
                   & + t103 * t106 / 0.108e3_dp ) * sc
              e_ndrho(ii) = e_ndrho(ii)&
                   + ( t4 * t124 / 0.72e2_dp ) * sc
           END IF

           t127 = 0.1e1_dp / t1 / t6
           t128 = a * t127
           t133 = 0.1e1_dp / t7
           t134 = a * t133
           t161 = t3 * t13
           t165 = 0.1e1_dp / t5
           t166 = t165 * t13
           t173 = t52 * t165
           t176 = t88 * t102
           t184 = b * t87 * c
           t185 = t102 * t13
           t189 = t184 * t44
           t192 = -0.560e3_dp * t93 - 0.432e3_dp * my_rho * d - 0.128e3_dp&
                & * t52 * my_rho * t13 * cf - 0.8e1_dp * t88 * t1 * t13 *&
                  cf - 0.560e3_dp * t14 * t93 * cf - 0.112e3_dp * t52 * t1&
                & * t62 - 0.8e1_dp * t88 * my_rho_1_3 * t62 - 0.432e3_dp *&
                  t14 * my_rho * cf * d - 0.14e2_dp / 0.9e1_dp * t52 *&
                  t161 * t22 - 0.11e2_dp / 0.9e1_dp * t88 * t166 * t22 -&
                  0.2e1_dp / 0.3e1_dp * t14 * t94 * t22 - 0.20e2_dp /&
                  0.9e1_dp * t173 * t80 - 0.2e1_dp * t176 * t80 -&
                  0.20e2_dp / 0.9e1_dp * t14 * t3 * t22 * d + 0.7e1_dp /&
                  0.9e1_dp * t184 * t185 * t22 + 0.7e1_dp / 0.9e1_dp *&
                  t189 * t80
           t193 = t192 * t40
           t196 = t98 * t105
           t199 = 0.1e1_dp / t6
           t200 = a * t199
           t201 = t39 ** 2
           t202 = 0.1e1_dp / t201
           t203 = t37 * t202
           t215 = t13 * my_ndrho * d
           t227 = 0.20e2_dp / 0.3e1_dp * t52 * t70 * my_ndrho + 0.4e1_dp *&
                  t14 * t11 * my_ndrho + 0.20e2_dp / 0.3e1_dp * t78 * t215&
                & + 0.20e2_dp / 0.3e1_dp * t14 * t69 * my_ndrho * d +&
                  0.14e2_dp / 0.3e1_dp * t88 * t89 * my_ndrho + 0.14e2_dp &
                &/ 0.3e1_dp * t95 * t215
           t228 = t227 * t40
           t231 = t123 * t105
           t245 = 0.6e1_dp * t14 * t1 + 0.20e2_dp * t14 * my_rho_1_3 * d +&
                  0.14e2_dp * t14 * my_rho_1_3 * c + 0.14e2_dp * t14 * c *&
                  d
           t246 = t245 * t40

           IF (grad_deriv>=2 .OR.grad_deriv==-2) THEN
              e_rho_rho(ii) = e_rho_rho(ii)&
                   + ( 0.5e1_dp / 0.81e2_dp * t128 * t41 - 0.5e1_dp / 0.108e3_dp&
                   & * t45 * t99 + t134 * t106 / 0.27e2_dp + t4 * t193 /&
                     0.72e2_dp - t103 * t196 / 0.54e2_dp + t200 * t203 /&
                     0.108e3_dp ) * sc
              e_ndrho_rho(ii) = e_ndrho_rho(ii)&
                   - ( 0.5e1_dp / 0.216e3_dp * t45 * t124 - t4 * t228 /&
                     0.72e2_dp + t103 * t231 / 0.108e3_dp ) * sc
              e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii)&
                   + ( t4 * t246 / 0.72e2_dp ) * sc
           END IF

           t248 = t5 ** 2
           t265 = 0.1e1_dp / t248
           t278 = t87 ** 2
           t279 = b * t278
           t332 = -0.432e3_dp * d - 0.2240e4_dp / 0.3e1_dp * my_rho_1_3 -&
                  0.74e2_dp / 0.27e2_dp * t184 * t127 * t80 + 0.100e3_dp /&
                  0.27e2_dp * t14 * t44 * t22 * d + 0.7e1_dp / 0.27e2_dp *&
                  t279 * t127 * t13 * t22 - 0.8e1_dp / 0.3e1_dp * t184 *&
                  t70 * cf - 0.2240e4_dp / 0.3e1_dp * t14 * my_rho_1_3 *&
                  cf - 0.48e2_dp * t88 * t11 * t13 * cf - 0.944e3_dp /&
                  0.3e1_dp * t52 * t61 - 0.40e2_dp * t88 * t69 * t62 -&
                  0.656e3_dp / 0.3e1_dp * t52 * t11 * t62 + 0.7e1_dp /&
                  0.27e2_dp * t279 * t265 * t80 + 0.64e2_dp / 0.27e2_dp *&
                  t52 * t44 * t13 * t22 - 0.8e1_dp / 0.3e1_dp * t184 * t77&
                & * t62 - 0.432e3_dp * t14 * cf * d + 0.52e2_dp /&
                  0.27e2_dp * t88 * t199 * t13 * t22 - 0.20e2_dp /&
                  0.9e1_dp * t184 * t133 * t13 * t22 + 0.8e1_dp / 0.9e1_dp&
                & * t14 * t102 * t22 + 0.100e3_dp / 0.27e2_dp * t52 * t199&
                & * t80 + 0.106e3_dp / 0.27e2_dp * t88 * t133 * t80

           IF (grad_deriv>=3 .OR. grad_deriv==-3) THEN
              e_rho_rho_rho(ii) = e_rho_rho_rho(ii)&
                   - ( 0.55e2_dp / 0.243e3_dp * a / t1 / t248 * t41 - 0.5e1_dp &
                   &/ 0.27e2_dp * t128 * t99 + 0.40e2_dp / 0.243e3_dp * a /&
                     my_rho_1_3 / t248 * t106 + 0.5e1_dp / 0.72e2_dp * t45 *&
                     t193 - t134 * t196 / 0.9e1_dp + 0.7e1_dp / 0.108e3_dp *&
                     a * t265 * t203 - t4 * t332 * t40 / 0.72e2_dp + t103 *&
                     t192 * t105 / 0.36e2_dp - t200 * t98 * t202 / 0.36e2_dp &
                   &+ t128 * t37 / t201 / t38 / 0.81e2_dp ) * sc
              e_ndrho_rho_rho(ii) = e_ndrho_rho_rho(ii)&
                   + ( 0.5e1_dp / 0.81e2_dp * t128 * t124 - 0.5e1_dp /&
                     0.108e3_dp * t45 * t228 + t134 * t231 / 0.27e2_dp + t4 *&
                   & (-0.28e2_dp / 0.9e1_dp * t52 * t161 * my_ndrho -&
                     0.22e2_dp / 0.9e1_dp * t88 * t166 * my_ndrho - 0.4e1_dp &
                   &/ 0.3e1_dp * t14 * t94 * my_ndrho - 0.40e2_dp / 0.9e1_dp &
                   &* t173 * t215 - 0.4e1_dp * t176 * t215 - 0.40e2_dp /&
                     0.9e1_dp * t14 * t3 * my_ndrho * d + 0.14e2_dp /&
                     0.9e1_dp * t184 * t185 * my_ndrho + 0.14e2_dp / 0.9e1_dp&
                   & * t189 * t215) * t40 / 0.72e2_dp - t103 * t227 * t105 /&
                     0.54e2_dp + t200 * t123 * t202 / 0.108e3_dp ) * sc
              e_ndrho_ndrho_rho(ii) = e_ndrho_ndrho_rho(ii)&
                   - ( 0.5e1_dp / 0.216e3_dp * t45 * t246 - t4 * (0.20e2_dp /&
                     0.3e1_dp * t52 * t70 + 0.4e1_dp * t14 * t11 + 0.20e2_dp &
                   &/ 0.3e1_dp * t52 * t89 * d + 0.20e2_dp / 0.3e1_dp * t14 *&
                     t69 * d + 0.14e2_dp / 0.3e1_dp * t88 * t89 + 0.14e2_dp /&
                     0.3e1_dp * t88 * t94 * t13 * d) * t40 / 0.72e2_dp + t103&
                   & * t245 * t105 / 0.108e3_dp ) * sc
           END IF
        END IF

!FM      t1 = my_rho_1_3 ** 2
!FM      t2 = t1 * my_rho
!FM      t3 = 0.1e1_dp / t2
!FM      t4 = a * t3
!FM      t5 = my_rho ** 2
!FM      t6 = t5 * my_rho
!FM      t7 = my_rho_1_3 * t6
!FM      t11 = 0.1e1_dp / my_rho_1_3
!FM      t13 = exp(-c * t11)
!FM      t14 = b * t13
!FM      t22 = my_ndrho ** 2
!FM      t26 = my_rho_1_3 * t22
!FM      t37 = -0.72e2_dp *( t7 - t14 *&
!FM           & t7 * cf - t6 * d * (1.0_dp+ t14 * cf)) + t14 *(t22*(0.3e1_dp &
!FM           & * t1 + 0.7e1_dp * c * d)&
!FM           + 0.10e2_dp * t26 * d + 0.7e1_dp &
!FM           &* t26 * c )
!FM      t38 = my_rho_1_3 + d
!FM      t39 = t38 ** 2
!FM      t40 = 0.1e1_dp / t39
!FM      t41 = t37 * t40
!FM
!FM     IF (grad_deriv>=0.AND.my_rho>epsilon_rho) THEN
!FM        e_0(ii) = e_0(ii)&
!FM             + t4 * t41 / 0.72e2_dp
!FM     END IF
!FM
!FM      t44 = 0.1e1_dp / (t1 * t5)
!FM      t45 = a * t44
!FM      t48 = my_rho_1_3 * t5
!FM      t52 = b * c
!FM      t61 = t13 * cf
!FM      t62 = t61 * d
!FM      t69 = 0.1e1_dp / t1
!FM      t70 = t69 * t13
!FM      t77 = 0.1e1_dp / my_rho
!FM      t78 = t52 * t77
!FM      t80 = t13 * t22 * d
!FM      t87 = c ** 2
!FM      t88 = b * t87
!FM      t89 = t77 * t13
!FM      t93 = my_rho_1_3 * my_rho
!FM      t94 = 0.1e1_dp / t93
!FM      t95 = t88 * t94
!FM      t98 = -0.216e3_dp * t5 *d -0.240e3_dp * t48(1.0_dp+ t14 * cf) -&
!FM           & 0.24e2_dp * t52 * (t2 * t62 +t13*t5*cf) &
!FM           - 0.216e3_dp * t14 * t5 * cf &
!FM           &* d + t22 *(0.10e2_dp / 0.3e1_dp * t52 * t70 + 0.2e1_dp *&
!FM           & t14 * t11 + 0.10e2_dp / 0.3e1_dp * t14 * t69 * d + 0.7e1_dp /&
!FM           & 0.3e1_dp * t88 * t89 ) +(0.10e2_dp / 0.3e1_dp * t78  +&
!FM           & 0.7e1_dp / 0.3e1_dp * t95 ) *&
!FM           & t80
!FM      t99 = t98 * t40
!FM      t102 = 0.1e1_dp / t48
!FM      t103 = a * t102
!FM      t105 = 0.1e1_dp / (t39 * t38)
!FM      t106 = t37 * t105
!FM      t112 = my_rho_1_3 * my_ndrho
!FM      t123 = t14 *(0.6e1_dp  t1 * my_ndrho + t112 * 0.20e2_dp &
!FM           &* d + 0.14e2_dp * c *(t112 + t14 *&
!FM           & my_ndrho * d))
!FM      t124 = t123 * t40
!FM
!FM     IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
!FM        e_rho(ii) = e_rho(ii) &
!FM             -0.5e1_dp / 0.216e3_dp * t45 * t41 + t4 * t99 / 0.72e2_dp&
!FM             & - t103 * t106 / 0.108e3_dp
!FM        e_ndrho(ii) = e_ndrho(ii)&
!FM             +t4 * t124 / 0.72e2_dp
!FM     END IF
!FM
!FM      t127 = 0.1e1_dp / (t1 * t6)
!FM      t128 = a * t127
!FM      t133 = 0.1e1_dp / t7
!FM      t134 = a * t133
!FM      t161 = t3 * t13
!FM      t165 = 0.1e1_dp / t5
!FM      t166 = t165 * t13
!FM      t173 = t52 * t165
!FM      t176 = t88 * t102
!FM      t184 = b * t87 * c
!FM      t185 = t102 * t13
!FM      t189 = t184 * t44
!FM      t192 = -0.560e3_dp * t93 - 0.432e3_dp * my_rho * d +cf*(&
!FM           -t13*(0.128e3_dp&
!FM           & * t52 * my_rho + 0.8e1_dp * t88 * t1 )&
!FM           & - 0.560e3_dp * t14 * t93 ) - 0.112e3_dp * t52 * t1&
!FM           & * t62 - 0.8e1_dp * t88 * my_rho_1_3 * t62 - 0.432e3_dp *&
!FM           & t14 * my_rho * cf * d - 0.14e2_dp / 0.9e1_dp * t52 *&
!FM           & t161 * t22 - 0.11e2_dp / 0.9e1_dp * t88 * t166 * t22 -&
!FM           & 0.2e1_dp / 0.3e1_dp * t14 * t94 * t22 - 0.20e2_dp /&
!FM           & 0.9e1_dp * t173 * t80 - 0.2e1_dp * t176 * t80 -&
!FM           & 0.20e2_dp / 0.9e1_dp * t14 * t3 * t22 * d + 0.7e1_dp /&
!FM           & 0.9e1_dp * t184 * t185 * t22 + 0.7e1_dp / 0.9e1_dp *&
!FM           & t189 * t80
!FM      t193 = t192 * t40
!FM      t196 = t98 * t105
!FM      t199 = 0.1e1_dp / t6
!FM      t200 = a * t199
!FM      t201 = t39 ** 2
!FM      t202 = 0.1e1_dp / t201
!FM      t203 = t37 * t202
!FM      t215 = t13 * my_ndrho * d
!FM      t227 = 0.20e2_dp / 0.3e1_dp * t52 * t70 * my_ndrho + 0.4e1_dp *&
!FM           & t14 * t11 * my_ndrho + 0.20e2_dp / 0.3e1_dp * t78 * t215&
!FM           & + 0.20e2_dp / 0.3e1_dp * t14 * t69 * my_ndrho * d +&
!FM           & 0.14e2_dp / 0.3e1_dp * t88 * t89 * my_ndrho + 0.14e2_dp &
!FM           &/ 0.3e1_dp * t95 * t215
!FM      t228 = t227 * t40
!FM      t231 = t123 * t105
!FM      t245 = 0.6e1_dp * t14 * t1 + 0.20e2_dp * t14 * my_rho_1_3 * d +&
!FM           & 0.14e2_dp * t14 * my_rho_1_3 * c + 0.14e2_dp * t14 * c *&
!FM           & d
!FM      t246 = t245 * t40
!FM
!FM     IF (grad_deriv>=2 .OR.grad_deriv==-2) THEN
!FM        e_rho_rho(ii) = e_rho_rho(ii)&
!FM             +0.5e1_dp / 0.81e2_dp * t128 * t41 - 0.5e1_dp / 0.108e3_dp&
!FM           & * t45 * t99 + t134 * t106 / 0.27e2_dp + t4 * t193 /&
!FM           & 0.72e2_dp - t103 * t196 / 0.54e2_dp + t200 * t203 /&
!FM           & 0.108e3_dp
!FM        e_ndrho_rho(ii) = e_ndrho_rho(ii)&
!FM             -0.5e1_dp / 0.216e3_dp * t45 * t124 + t4 * t228 /&
!FM           & 0.72e2_dp - t103 * t231 / 0.108e3_dp
!FM        e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii)&
!FM             +t4 * t246 / 0.72e2_dp
!FM     END IF
!FM
!FM      t248 = t5 ** 2
!FM      t265 = 0.1e1_dp / t248
!FM      t278 = t87 ** 2
!FM      t279 = b * t278
!FM      t332 = -0.432e3_dp * d - 0.2240e4_dp / 0.3e1_dp * my_rho_1_3 -&
!FM           & 0.74e2_dp / 0.27e2_dp * t184 * t127 * t80 + 0.100e3_dp /&
!FM           & 0.27e2_dp * t14 * t44 * t22 * d + 0.7e1_dp / 0.27e2_dp *&
!FM           & t279 * t127 * t13 * t22 - 0.8e1_dp / 0.3e1_dp * t184 *&
!FM           & t70 * cf - 0.2240e4_dp / 0.3e1_dp * t14 * my_rho_1_3 *&
!FM           & cf - 0.48e2_dp * t88 * t11 * t13 * cf - 0.944e3_dp /&
!FM           & 0.3e1_dp * t52 * t61 - 0.40e2_dp * t88 * t69 * t62 -&
!FM           & 0.656e3_dp / 0.3e1_dp * t52 * t11 * t62 + 0.7e1_dp /&
!FM           & 0.27e2_dp * t279 * t265 * t80 + 0.64e2_dp / 0.27e2_dp *&
!FM           & t52 * t44 * t13 * t22 - 0.8e1_dp / 0.3e1_dp * t184 * t77&
!FM           & * t62 - 0.432e3_dp * t14 * cf * d + 0.52e2_dp /&
!FM           & 0.27e2_dp * t88 * t199 * t13 * t22 - 0.20e2_dp /&
!FM           & 0.9e1_dp * t184 * t133 * t13 * t22 + 0.8e1_dp / 0.9e1_dp&
!FM           & * t14 * t102 * t22 + 0.100e3_dp / 0.27e2_dp * t52 * t199&
!FM           & * t80 + 0.106e3_dp / 0.27e2_dp * t88 * t133 * t80
!FM
!FM     IF (grad_deriv>=3 .OR. grad_deriv==-3) THEN
!FM        e_rho_rho_rho(ii) = e_rho_rho_rho(ii)&
!FM             -0.55e2_dp / 0.243e3_dp * a / t1 / t248 * t41 + 0.5e1_dp &
!FM           &/ 0.27e2_dp * t128 * t99 - 0.40e2_dp / 0.243e3_dp * a /&
!FM           & my_rho_1_3 / t248 * t106 - 0.5e1_dp / 0.72e2_dp * t45 *&
!FM           & t193 + t134 * t196 / 0.9e1_dp - 0.7e1_dp / 0.108e3_dp *&
!FM           & a * t265 * t203 + t4 * t332 * t40 / 0.72e2_dp - t103 *&
!FM           & t192 * t105 / 0.36e2_dp + t200 * t98 * t202 / 0.36e2_dp &
!FM           &- t128 * t37 / t201 / t38 / 0.81e2_dp
!FM        e_ndrho_rho_rho(ii) = e_ndrho_rho_rho(ii)&
!FM             +0.5e1_dp / 0.81e2_dp * t128 * t124 - 0.5e1_dp /&
!FM           & 0.108e3_dp * t45 * t228 + t134 * t231 / 0.27e2_dp + t4 *&
!FM           & (-0.28e2_dp / 0.9e1_dp * t52 * t161 * my_ndrho -&
!FM           & 0.22e2_dp / 0.9e1_dp * t88 * t166 * my_ndrho - 0.4e1_dp &
!FM           &/ 0.3e1_dp * t14 * t94 * my_ndrho - 0.40e2_dp / 0.9e1_dp &
!FM           &* t173 * t215 - 0.4e1_dp * t176 * t215 - 0.40e2_dp /&
!FM           & 0.9e1_dp * t14 * t3 * my_ndrho * d + 0.14e2_dp /&
!FM           & 0.9e1_dp * t184 * t185 * my_ndrho + 0.14e2_dp / 0.9e1_dp&
!FM           & * t189 * t215) * t40 / 0.72e2_dp - t103 * t227 * t105 /&
!FM           & 0.54e2_dp + t200 * t123 * t202 / 0.108e3_dp
!FM        e_ndrho_ndrho_rho(ii) = e_ndrho_ndrho_rho(ii)&
!FM              -0.5e1_dp / 0.216e3_dp * t45 * t246 + t4 * (0.20e2_dp /&
!FM           & 0.3e1_dp * t52 * t70 + 0.4e1_dp * t14 * t11 + 0.20e2_dp &
!FM           &/ 0.3e1_dp * t52 * t89 * d + 0.20e2_dp / 0.3e1_dp * t14 *&
!FM           & t69 * d + 0.14e2_dp / 0.3e1_dp * t88 * t89 + 0.14e2_dp /&
!FM           & 0.3e1_dp * t88 * t94 * t13 * d) * t40 / 0.72e2_dp - t103&
!FM           & * t245 * t105 / 0.108e3_dp
!FM     END IF

     END DO

     !$omp end do

  END SELECT
END SUBROUTINE lyp_lda_calc

! *****************************************************************************
!> \brief evaluates the becke 88 exchange functional for lsd
!> \param rho_set ...
!> \param deriv_set ...
!> \param grad_deriv ...
!> \param lyp_params input parameter (scaling)
!> \par History
!>      11.2003 created [fawzi]
!>      01.2007 added scaling [Manuel Guidon]
!> \author fawzi
! *****************************************************************************
SUBROUTINE lyp_lsd_eval(rho_set,deriv_set,grad_deriv,lyp_params)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: lyp_params

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_drho, epsilon_rho, sc
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndr, &
      e_ndr_ndr, e_ndr_ra, e_ndr_rb, e_ndra, e_ndra_ndra, e_ndra_ra, &
      e_ndra_rb, e_ndrb, e_ndrb_ndrb, e_ndrb_ra, e_ndrb_rb, e_ra, e_ra_ra, &
      e_ra_rb, e_rb, e_rb_rb, norm_drho, norm_drhoa, norm_drhob, rhoa, rhob
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
  failure=.FALSE.
  NULLIFY(deriv, bo)

  CALL section_vals_val_get(lyp_params,"scale_c",r_val=sc)
  CALL cite_reference(Lee1988)

  CPASSERT(ASSOCIATED(rho_set))
  CPASSERT(rho_set%ref_count>0)
  CPASSERT(ASSOCIATED(deriv_set))
  CPASSERT(deriv_set%ref_count>0)
  CALL xc_rho_set_get(rho_set,&
       rhoa=rhoa, rhob=rhob, norm_drhoa=norm_drhoa, &
       norm_drhob=norm_drhob, norm_drho=norm_drho, &
       rho_cutoff=epsilon_rho,&
       drho_cutoff=epsilon_drho, local_bounds=bo)
  npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

  ! meaningful default for the arrays we don't need: let us make compiler
  ! and debugger happy...
  IF (cp_debug) THEN
     ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)))
  ELSE
     dummy=> rhoa
  END IF
  e_0 => dummy
  e_ra => dummy
  e_rb => dummy
  e_ndra_ra => dummy
  e_ndra_rb => dummy
  e_ndrb_ra => dummy
  e_ndrb_rb => dummy
  e_ndr_ndr => dummy
  e_ndra_ndra => dummy
  e_ndrb_ndrb => dummy
  e_ndr => dummy
  e_ndra => dummy
  e_ndrb => dummy
  e_ra_ra => dummy
  e_ra_rb => dummy
  e_rb_rb => dummy
  e_ndr_ra => dummy
  e_ndr_rb => dummy

  IF (grad_deriv>=0) THEN
     deriv => xc_dset_get_derivative(deriv_set,"",&
          allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv, deriv_data=e_0)
  END IF
  IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
     deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
          allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ra)
     deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
          allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_rb)
     deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
          allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ndr)
     deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)",&
          allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ndra)
     deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)",&
          allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ndrb)
  END IF
  IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
     deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)",&
          allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ra_ra)
     deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhob)",&
          allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ra_rb)
     deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)",&
          allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_rb_rb)
     deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rhoa)",&
          allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ndr_ra)
     deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rhob)",&
          allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ndr_rb)
     deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhoa)",&
          allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ndra_ra)
     deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhob)",&
          allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ndra_rb)
     deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhob)",&
          allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ndrb_rb)
     deriv => xc_dset_get_derivative(deriv_set,&
          "(norm_drho)(norm_drho)", allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ndr_ndr)
     deriv => xc_dset_get_derivative(deriv_set,&
          "(norm_drhoa)(norm_drhoa)", allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ndra_ndra)
     deriv => xc_dset_get_derivative(deriv_set,&
          "(norm_drhob)(norm_drhob)", allocate_deriv=.TRUE.)
     CALL xc_derivative_get(deriv,deriv_data=e_ndrb_ndrb)
  END IF

  !$omp parallel default(none) &
  !$omp          shared(rhoa, rhob, norm_drho, norm_drhoa, norm_drhob) &
  !$omp          shared(e_0, e_ra, e_rb, e_ndra_ra, e_ndra_rb, e_ndrb_ra) &
  !$omp          shared(e_ndrb_rb, e_ndr_ndr, e_ndra_ndra, e_ndrb_ndrb) &
  !$omp          shared(e_ndr, e_ndra, e_ndrb, e_ra_ra, e_ra_rb, e_rb_rb) &
  !$omp          shared(e_ndr_ra, e_ndr_rb, grad_deriv) &
  !$omp          shared(npoints, epsilon_rho, sc)

  CALL lyp_lsd_calc(&
       rhoa=rhoa, rhob=rhob, norm_drho=norm_drho, norm_drhoa=norm_drhoa,&
       norm_drhob=norm_drhob, e_0=e_0, e_ra=e_ra, e_rb=e_rb,&
         e_ndra_ra=e_ndra_ra, e_ndra_rb=e_ndra_rb, e_ndrb_ra&
       &=e_ndrb_ra, e_ndrb_rb=e_ndrb_rb, e_ndr_ndr=e_ndr_ndr,&
       e_ndra_ndra=e_ndra_ndra, e_ndrb_ndrb=e_ndrb_ndrb, e_ndr=e_ndr,&
       e_ndra=e_ndra, e_ndrb=e_ndrb, e_ra_ra=e_ra_ra, &
       e_ra_rb=e_ra_rb, e_rb_rb=e_rb_rb, e_ndr_ra=e_ndr_ra,&
       e_ndr_rb=e_ndr_rb,&
       grad_deriv=grad_deriv, npoints=npoints, &
       epsilon_rho=epsilon_rho,sc=sc)

  !$omp end parallel

  IF (cp_debug) THEN
     DEALLOCATE(dummy)
  ELSE
     NULLIFY(dummy)
  END IF
  CALL timestop(handle)
END SUBROUTINE lyp_lsd_eval

! *****************************************************************************
!> \brief low level calculation of the becke 88 exchange functional for lsd
!> \param rhoa alpha spin density
!> \param rhob beta spin density
!> \param norm_drho || grad rhoa + grad rhob ||
!> \param norm_drhoa || grad rhoa ||
!> \param norm_drhob || grad rhob ||
!> \param e_0 adds to it the local value of the functional
!> \param e_ra e_*: derivative of the functional wrt. to the variables
!>        named where the * is.
!> \param e_rb ...
!> \param e_ndra_ra ...
!> \param e_ndra_rb ...
!> \param e_ndrb_ra ...
!> \param e_ndrb_rb ...
!> \param e_ndr_ndr ...
!> \param e_ndra_ndra ...
!> \param e_ndrb_ndrb ...
!> \param e_ndr ...
!> \param e_ndra ...
!> \param e_ndrb ...
!> \param e_ra_ra ...
!> \param e_ra_rb ...
!> \param e_rb_rb ...
!> \param e_ndr_ra ...
!> \param e_ndr_rb ...
!> \param grad_deriv ...
!> \param npoints ...
!> \param epsilon_rho ...
!> \param sc scaling parameter for correlation
!> \par History
!>      01.2004 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE lyp_lsd_calc(rhoa, rhob, norm_drho, norm_drhoa, norm_drhob,&
     e_0, e_ra, e_rb, e_ndra_ra, e_ndra_rb, e_ndrb_ra, e_ndrb_rb, e_ndr_ndr,&
     e_ndra_ndra, e_ndrb_ndrb, e_ndr,&
     e_ndra, e_ndrb, e_ra_ra, e_ra_rb, e_rb_rb, e_ndr_ra, e_ndr_rb,&
     grad_deriv,npoints,epsilon_rho,sc)
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rhoa, rhob, norm_drho, &
                                                norm_drhoa, norm_drhob
    REAL(kind=dp), DIMENSION(*), INTENT(inout) :: e_0, e_ra, e_rb, e_ndra_ra, &
      e_ndra_rb, e_ndrb_ra, e_ndrb_rb, e_ndr_ndr, e_ndra_ndra, e_ndrb_ndrb, &
      e_ndr, e_ndra, e_ndrb, e_ra_ra, e_ra_rb, e_rb_rb, e_ndr_ra, e_ndr_rb
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sc

    CHARACTER(len=*), PARAMETER :: routineN = 'lyp_lsd_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii
    LOGICAL                                  :: failure
    REAL(kind=dp) :: cf, my_ndrho, my_ndrhoa, my_ndrhob, my_rho, my_rhoa, &
      my_rhob, t1, t100, t101, t102, t103, t104, t107, t108, t109, t112, &
      t115, t118, t12, t120, t124, t129, t13, t130, t132, t135, t138, t14, &
      t140, t142, t145, t15, t153, t155, t159, t16, t164, t168, t17, t171, &
      t176, t18, t181, t182, t185, t189, t194, t195, t198, t2, t20, t202, &
      t205, t206, t21, t210, t214, t215, t218, t22, t222, t228, t23, t232, &
      t236, t238, t24, t243, t249, t25, t252, t253, t255, t257, t26, t260, &
      t265, t268, t27, t270, t29, t3, t30, t304, t31, t310, t313, t316, t319, &
      t322, t326, t329, t332, t334, t341, t35, t354
    REAL(kind=dp) :: t356, t360, t363, t37, t373, t376, t381, t39, t391, t4, &
      t40, t408, t41, t415, t419, t422, t430, t44, t445, t449, t45, t452, &
      t46, t467, t47, t48, t483, t487, t49, t490, t5, t50, t505, t515, t519, &
      t52, t520, t54, t56, t57, t6, t61, t62, t64, t66, t7, t72, t75, t78, &
      t8, t82, t85, t86, t87, t88, t90, t92, t94, t95, t98, t99

  failure=.FALSE.
  cf=0.3_dp*(3._dp*pi*pi)**(2._dp/3._dp)

  !$omp do

  DO ii=1,npoints
     my_rhoa=MAX(rhoa(ii),0.0_dp)
     my_rhob=MAX(rhob(ii),0.0_dp)
     my_rho=my_rhoa+my_rhob
     IF (my_rho>epsilon_rho) THEN
     my_ndrho=norm_drho(ii)
     my_ndrhoa=norm_drhoa(ii)
     my_ndrhob=norm_drhob(ii)

     t1 = my_rho ** (0.1e1_dp / 0.3e1_dp)
     t2 = 0.1e1_dp / t1
     t3 = d * t2
     t4 = 0.1e1_dp + t3
     t5 = 0.1e1_dp / t4
     t6 = a * t5
     t7 = my_rhoa * my_rhob
     t8 = 0.1e1_dp / my_rho
     t12 = a * b
     t13 = c * t2
     t14 = EXP(-t13)
     t15 = t12 * t14
     t16 = my_rho ** 2
     t17 = t16 * my_rho
     t18 = t1 ** 2
     t20 = 0.1e1_dp / t18 / t17
     t21 = t5 * t20
     t22 = 2 ** (0.1e1_dp / 0.3e1_dp)
     t23 = t22 ** 2
     t24 = t23 * cf
     t25 = my_rhoa ** 2
     t26 = my_rhoa ** (0.1e1_dp / 0.3e1_dp)
     t27 = t26 ** 2
     t29 = my_rhob ** 2
     t30 = my_rhob ** (0.1e1_dp / 0.3e1_dp)
     t31 = t30 ** 2
     t35 = 0.8e1_dp * t24 * (t27 * t25 + t31 * t29)
     t37 = t3 * t5
     t39 = 0.47e2_dp / 0.18e2_dp - 0.7e1_dp / 0.18e2_dp * t13 -&
            0.7e1_dp / 0.18e2_dp * t37
     t40 = my_ndrho ** 2
     t41 = t39 * t40
     t44 = 0.5e1_dp / 0.2e1_dp - t13 / 0.18e2_dp - t37 / 0.18e2_dp
     t45 = my_ndrhoa ** 2
     t46 = my_ndrhob ** 2
     t47 = t45 + t46
     t48 = t44 * t47
     t49 = t13 + t37 - 0.11e2_dp
     t50 = my_rhoa * t8
     t52 = my_rhob * t8
     t54 = t50 * t45 + t52 * t46
     t56 = t49 * t54 / 0.9e1_dp
     t57 = t35 + t41 - t48 - t56
     t61 = 0.2e1_dp / 0.3e1_dp * t16
     t62 = t61 - t25
     t64 = t61 - t29
     t66 = t7 * t57 - 0.2e1_dp / 0.3e1_dp * t16 * t40 + t62 * t46 +&
            t64 * t45

     IF (grad_deriv>=0 .AND. my_rho>epsilon_rho) THEN
        e_0(ii) = e_0(ii)&
             - ( 0.4e1_dp * t6 * t7 * t8 + t15 * t21 * t66 ) * sc
     END IF
     !--------

     t72 = t27 * my_rhoa
     t75 = t49 * t8
     t78 = 0.64e2_dp / 0.3e1_dp * t24 * t72 - t75 * t45 / 0.9e1_dp
     t82 = my_rhob * t57 + t7 * t78 - 0.2e1_dp * my_rhoa * t46
     t85 = t4 ** 2
     t86 = 0.1e1_dp / t85
     t87 = a * t86
     t88 = t87 * my_rhoa
     t90 = 0.1e1_dp / t1 / t16
     t92 = my_rhob * t90 * d
     t94 = 0.4e1_dp / 0.3e1_dp * t88 * t92
     t95 = 0.1e1_dp / t16
     t98 = 0.4e1_dp * t6 * t7 * t95
     t99 = t12 * c
     t100 = t16 ** 2
     t101 = t100 * my_rho
     t102 = 0.1e1_dp / t101
     t103 = t102 * t14
     t104 = t5 * t66
     t107 = t99 * t103 * t104 / 0.3e1_dp
     t108 = t86 * t102
     t109 = t66 * d
     t112 = t15 * t108 * t109 / 0.3e1_dp
     t115 = t5 / t18 / t100
     t118 = 0.11e2_dp / 0.3e1_dp * t15 * t115 * t66
     t120 = 0.1e1_dp / t1 / my_rho
     t124 = d ** 2
     t129 = c * t120 + d * t120 * t5 - t124 / t18 / my_rho * t86
     t130 = 0.7e1_dp / 0.54e2_dp * t129
     t132 = t129 / 0.54e2_dp
     t135 = -t129 / 0.3e1_dp
     t138 = my_rhoa * t95
     t140 = my_rhob * t95
     t142 = -t138 * t45 - t140 * t46
     t145 = t130 * t40 - t132 * t47 - t135 * t54 / 0.9e1_dp - t49 *&
            t142 / 0.9e1_dp
     t153 = t7 * t145 - 0.4e1_dp / 0.3e1_dp * my_rho * t40 +&
            0.4e1_dp / 0.3e1_dp * my_rho * t46 + 0.4e1_dp / 0.3e1_dp&
          & * my_rho * t45
     t155 = t15 * t21 * t153

     IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
        e_ra(ii)=e_ra(ii)&
             - ( 0.4e1_dp * t6 * t52 + t15 * t21 * t82 + t94 - t98 +&
               t107 + t112 - t118 + t155 ) * sc
     END IF

     t159 = t31 * my_rhob
     t164 = 0.64e2_dp / 0.3e1_dp * t24 * t159 - t75 * t46 / 0.9e1_dp
     t168 = my_rhoa * t57 + t7 * t164 - 0.2e1_dp * my_rhob * t45

     IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
        e_rb(ii)=e_rb(ii)&
             - ( 0.4e1_dp * t6 * t50 + t15 * t21 * t168 + t94 - t98 +&
               t107 + t112 - t118 + t155 ) * sc
     END IF

     t171 = t39 * my_ndrho
     t176 = 0.2e1_dp * t7 * t171 - 0.4e1_dp / 0.3e1_dp * t16 * my_ndrho

     IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
        e_ndr(ii)=e_ndr(ii)&
             - ( t15 * t21 * t176 ) * sc
     END IF

     t181 = t49 * my_rhoa
     t182 = t8 * my_ndrhoa
     t185 = -0.2e1_dp * t44 * my_ndrhoa - 0.2e1_dp / 0.9e1_dp * t181 * t182
     t189 = t7 * t185 + 0.2e1_dp * t64 * my_ndrhoa

     IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
        e_ndra(ii) = e_ndra(ii)&
             - ( t15 * t21 * t189 ) * sc
     END IF

     t194 = t49 * my_rhob
     t195 = t8 * my_ndrhob
     t198 = -0.2e1_dp * t44 * my_ndrhob - 0.2e1_dp / 0.9e1_dp * t194 * t195
     t202 = t7 * t198 + 0.2e1_dp * t62 * my_ndrhob

     IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
        e_ndrb(ii) = e_ndrb(ii)&
             - ( t15 * t21 * t202 ) * sc
     END IF
     !-------

     t205 = t100 * t16
     t206 = 0.1e1_dp / t205
     t210 = 0.26e2_dp / 0.9e1_dp * t99 * t206 * t14 * t104
     t214 = 0.2e1_dp / 0.3e1_dp * t99 * t103 * t5 * t153
     t215 = c ** 2
     t218 = 0.1e1_dp / t1 / t205
     t222 = t12 * t215 * t218 * t14 * t104 / 0.9e1_dp
     t228 = 0.2e1_dp / 0.9e1_dp * t12 * c * t218 * t14 * t86 * t109
     t232 = 0.26e2_dp / 0.9e1_dp * t15 * t86 * t206 * t109
     t236 = 0.2e1_dp / 0.3e1_dp * t15 * t108 * t153 * d
     t238 = 0.1e1_dp / t85 / t4
     t243 = 0.2e1_dp / 0.9e1_dp * t15 * t238 * t218 * t66 * t124
     t249 = 0.154e3_dp / 0.9e1_dp * t15 * t5 / t18 / t101 * t66
     t252 = 0.22e2_dp / 0.3e1_dp * t15 * t115 * t153
     t253 = t6 * t140
     t255 = t87 * t92
     t257 = c * t90
     t260 = d * t90 * t5
     t265 = t124 / t18 / t16 * t86
     t268 = 0.1e1_dp / t17
     t270 = t124 * d * t268 * t238
     t304 = t15 * t21 * (t7 * ((-0.14e2_dp / 0.81e2_dp * t257 -&
            0.14e2_dp / 0.81e2_dp * t260 + 0.7e1_dp / 0.27e2_dp *&
            t265 - 0.7e1_dp / 0.81e2_dp * t270) * t40 - (-0.2e1_dp /&
            0.81e2_dp * t257 - 0.2e1_dp / 0.81e2_dp * t260 + t265 /&
            0.27e2_dp - t270 / 0.81e2_dp) * t47 - (0.4e1_dp /&
            0.9e1_dp * t257 + 0.4e1_dp / 0.9e1_dp * t260 - 0.2e1_dp &
          &/ 0.3e1_dp * t265 + 0.2e1_dp / 0.9e1_dp * t270) * t54 /&
            0.9e1_dp - 0.2e1_dp / 0.9e1_dp * t135 * t142 - t49 *&
          & (0.2e1_dp * my_rhoa * t268 * t45 + 0.2e1_dp * my_rhob *&
            t268 * t46) / 0.9e1_dp) - 0.4e1_dp / 0.3e1_dp * t40 +&
            0.4e1_dp / 0.3e1_dp * t46 + 0.4e1_dp / 0.3e1_dp * t45)
     t310 = 0.40e2_dp / 0.9e1_dp * t88 * my_rhob / t1 / t17 * d
     t313 = my_rhob * t20
     t316 = 0.8e1_dp / 0.9e1_dp * a * t238 * my_rhoa * t313 * t124
     t319 = 0.8e1_dp * t6 * t7 * t268
     t322 = t99 * t103 * t5 * t82
     t326 = t15 * t108 * t82 * d
     t329 = t15 * t115 * t82
     t332 = t135 * t8
     t334 = t49 * t95
     t341 = t15 * t21 * (my_rhob * t145 + t7 * (-t332 * t45 /&
            0.9e1_dp + t334 * t45 / 0.9e1_dp))

     IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
        e_ra_ra(ii) = e_ra_ra(ii)&
             + ( t210 - t214 - t222 - t228 + t232 - t236 - t243 - t249 +&
               t252 + 0.8e1_dp * t253 - 0.8e1_dp / 0.3e1_dp * t255 -&
               t304 + t310 - t316 - t319 - 0.2e1_dp / 0.3e1_dp * t322 -&
               0.2e1_dp / 0.3e1_dp * t326 + 0.22e2_dp / 0.3e1_dp * t329&
             & - 0.2e1_dp * t341 - t15 * t21 * (0.2e1_dp * my_rhob *&
               t78 + 0.320e3_dp / 0.9e1_dp * t72 * my_rhob * t24 -&
               0.2e1_dp * t46) ) * sc
     END IF

     t354 = t99 * t103 * t5 * t168
     t356 = t6 * t138
     t360 = t87 * my_rhoa * t90 * d
     t363 = t15 * t115 * t168
     t373 = t15 * t21 * (my_rhoa * t145 + t7 * (-t332 * t46 /&
            0.9e1_dp + t334 * t46 / 0.9e1_dp))
     t376 = t15 * t108 * t168 * d

     IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
        e_rb_rb(ii)=e_rb_rb(ii)&
             + ( t210 - t214 - t222 - t228 + t232 - t236 - t243 - t249 +&
               t252 - t304 + t310 - t316 - t319 + 0.8e1_dp * t356 -&
               0.8e1_dp / 0.3e1_dp * t360 + 0.22e2_dp / 0.3e1_dp *&
               t363 - 0.2e1_dp * t373 - 0.2e1_dp / 0.3e1_dp * t354 -&
               0.2e1_dp / 0.3e1_dp * t376 - t15 * t21 * (0.2e1_dp *&
               my_rhoa * t164 + 0.320e3_dp / 0.9e1_dp * my_rhoa *&
               t159 * t24 - 0.2e1_dp * t45) ) * sc
     END IF

     t381 = -t354 / 0.3e1_dp + 0.4e1_dp * t356 - 0.4e1_dp / 0.3e1_dp&
          & * t360 + 0.11e2_dp / 0.3e1_dp * t363 - t373 - t341 -&
            t376 / 0.3e1_dp + t310 - 0.4e1_dp * t6 * t8 + 0.4e1_dp *&
            t253 + t210 - t214 - t222
     t391 = -t228 + t232 - t236 - t243 - t249 + t252 - 0.4e1_dp /&
            0.3e1_dp * t255 - t304 - t316 - t319 - t322 / 0.3e1_dp -&
            t326 / 0.3e1_dp + 0.11e2_dp / 0.3e1_dp * t329 - t15 *&
            t21 * (t35 + t41 - t48 - t56 + my_rhob * t164 + my_rhoa &
          &* t78)

     IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
        e_ra_rb(ii)=e_ra_rb(ii)&
             + ( t381 + t391 ) * sc
     END IF

     t408 = t12 * t14 * t5
     t415 = t99 * t103 * t5 * t176 / 0.3e1_dp
     t419 = t15 * t108 * t176 * d / 0.3e1_dp
     t422 = 0.11e2_dp / 0.3e1_dp * t15 * t115 * t176
     t430 = t15 * t21 * (0.2e1_dp * t7 * t130 * my_ndrho - 0.8e1_dp &
          &/ 0.3e1_dp * my_rho * my_ndrho)

     IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
        e_ndr_ra(ii) = e_ndr_ra(ii)&
             - ( 0.2e1_dp * t408 * t313 * t171 + t415 + t419 - t422 +&
               t430 ) * sc
        e_ndr_rb(ii) = e_ndr_rb(ii)&
             - ( 0.2e1_dp * t408 * t20 * my_rhoa * t171 + t415 + t419 -&
               t422 + t430 ) * sc
     END IF

     t445 = t99 * t103 * t5 * t189 / 0.3e1_dp
     t449 = t15 * t108 * t189 * d / 0.3e1_dp
     t452 = 0.11e2_dp / 0.3e1_dp * t15 * t115 * t189
     t467 = t15 * t21 * (t7 * (-0.2e1_dp * t132 * my_ndrhoa -&
            0.2e1_dp / 0.9e1_dp * t135 * my_rhoa * t182 + 0.2e1_dp /&
            0.9e1_dp * t181 * t95 * my_ndrhoa) + 0.8e1_dp / 0.3e1_dp&
          & * my_rho * my_ndrhoa)

     IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
        e_ndra_ra(ii)=e_ndra_ra(ii)&
             - ( t15 * t21 * (my_rhob * t185 - 0.2e1_dp / 0.9e1_dp * t7&
             & * t75 * my_ndrhoa) + t445 + t449 - t452 + t467 ) * sc
        e_ndra_rb(ii)=e_ndra_rb(ii)&
             - ( t15 * t21 * (my_rhoa * t185 - 0.4e1_dp * my_rhob *&
               my_ndrhoa) + t445 + t449 - t452 + t467 ) * sc
     END IF

     t483 = t99 * t103 * t5 * t202 / 0.3e1_dp
     t487 = t15 * t108 * t202 * d / 0.3e1_dp
     t490 = 0.11e2_dp / 0.3e1_dp * t15 * t115 * t202
     t505 = t15 * t21 * (t7 * (-0.2e1_dp * t132 * my_ndrhob -&
            0.2e1_dp / 0.9e1_dp * t135 * my_rhob * t195 + 0.2e1_dp /&
            0.9e1_dp * t194 * t95 * my_ndrhob) + 0.8e1_dp / 0.3e1_dp&
          & * my_rho * my_ndrhob)

     IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
        e_ndrb_ra(ii) =e_ndrb_ra(ii)&
             - ( t15 * t21 * (my_rhob * t198 - 0.4e1_dp * my_rhoa *&
               my_ndrhob) + t483 + t487 - t490 + t505 ) * sc
        e_ndrb_rb(ii)=e_ndrb_rb(ii)&
             - ( t15 * t21 * (my_rhoa * t198 - 0.2e1_dp / 0.9e1_dp * t7&
             & * t75 * my_ndrhob) + t483 + t487 - t490 + t505 ) * sc
     END IF

     t515 = 0.4e1_dp / 0.3e1_dp * t16

     IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
        e_ndr_ndr(ii) = e_ndr_ndr(ii)&
             - ( t15 * t21 * (0.2e1_dp * t7 * t39 - t515) ) * sc
     END IF

     t519 = t13 / 0.9e1_dp
     t520 = t37 / 0.9e1_dp

     IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
        e_ndra_ndra(ii) = e_ndra_ndra(ii)&
             - ( t15 * t21 * (t7 * (-0.5e1_dp + t519 + t520 - 0.2e1_dp &
             &/ 0.9e1_dp * t181 * t8) + t515 - 0.2e1_dp * t29) ) * sc

        e_ndrb_ndrb(ii) = e_ndrb_ndrb(ii)&
             - ( t15 * t21 * (t7 * (-0.5e1_dp + t519 + t520 - 0.2e1_dp &
             &/ 0.9e1_dp * t194 * t8) + t515 - 0.2e1_dp * t25) ) * sc
     END IF
  END IF
 END DO
 !$omp end do

END SUBROUTINE lyp_lsd_calc

END MODULE xc_lyp
