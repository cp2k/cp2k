!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2003  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/xc_xbecke88 [1.0] *
!!
!!   NAME
!!     xc_xbecke88
!!
!!   FUNCTION
!!     calculates the Becke 88 exchange functional
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE xc_xbecke88
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger, cp_a_l, cp_error_check
  use kinds, only: dp
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
  USE xc_derivative_set_types, ONLY: xc_derivative_set_type, xc_dset_get_derivative
  USE xc_derivative_types, ONLY: xc_derivative_get, xc_derivative_type
  USE xc_rho_set_types, ONLY: xc_rho_set_type, xc_rho_set_get, xc_rho_cflags_type,&
       xc_rho_cflags_setall
  USE cp_array_r_utils, only: cp_3d_r_p_type
  USE mathconstants,                   ONLY: pi
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='xc_xbecke88'
  REAL(kind=dp), PARAMETER :: beta=0.0042_dp
  
  PUBLIC :: xb88_lda_info, xb88_lsd_info, xb88_lda_eval, xb88_lsd_eval
!!***
!****************************************************************************
contains

!!****f* xc_xbecke88/xb88_lda_info [1.0] *
!!
!!   NAME
!!     xb88_lda_info
!!
!!   FUNCTION
!!     return various information on the functional
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - reference: string with the reference of the actual functional
!!     - shortform: string with the shortform of the functional name
!!     - needs: the components needed by this functional are set to
!!       true (does not set the unneeded components to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE xb88_lda_info(reference,shortform, needs, max_deriv, error)
  CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
  CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform
  TYPE(xc_rho_cflags_type), INTENT(inout), optional :: needs
  INTEGER, INTENT(out), OPTIONAL :: max_deriv
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='xb88_lda_info',&
       routineP=moduleN//':'//routineN
  
  IF ( PRESENT ( reference ) ) THEN
     reference = "A. Becke, Phys. Rev. A 38, 3098 (1988) {LDA version}"
  END IF
  IF ( PRESENT ( shortform ) ) THEN
     shortform = "Becke 1988 Exchange Functional (LDA)"
  END IF
  IF (PRESENT(needs)) THEN
     needs%rho=.true.
     needs%rho_1_3=.true.
     needs%norm_drho=.true.
  END IF
  if (present(max_deriv)) max_deriv=3

END SUBROUTINE xb88_lda_info
!***************************************************************************

!!****f* xc_xbecke88/xb88_lsd_info [1.0] *
!!
!!   NAME
!!     xb88_lsd_info
!!
!!   FUNCTION
!!     return various information on the functional
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - reference: string with the reference of the actual functional
!!     - shortform: string with the shortform of the functional name
!!     - needs: the components needed by this functional are set to
!!       true (does not set the unneeded components to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE xb88_lsd_info(reference,shortform, needs, max_deriv, error)
  CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: reference
  CHARACTER (LEN=*), INTENT(OUT), OPTIONAL :: shortform
  TYPE(xc_rho_cflags_type), INTENT(inout), optional :: needs
  INTEGER, INTENT(out), OPTIONAL :: max_deriv
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='xb88_lsd_info',&
       routineP=moduleN//':'//routineN
  
  IF ( PRESENT ( reference ) ) THEN
     reference = "A. Becke, Phys. Rev. A 38, 3098 (1988) {LSD version}"
  END IF
  IF ( PRESENT ( shortform ) ) THEN
     shortform = "Becke 1988 Exchange Functional (LSD)"
  END IF
  IF (PRESENT(needs)) THEN
     needs%rho_spin=.true.
     needs%rho_spin_1_3=.true.
     needs%norm_drho_spin=.true.
  END IF
  if (present(max_deriv)) max_deriv=3

END SUBROUTINE xb88_lsd_info
!***************************************************************************

!!****f* xc_xbecke88/xb88_lda_eval [1.0] *
!!
!!   NAME
!!     xb88_lda_eval
!!
!!   FUNCTION
!!     evaluates the becke 88 exchange functional for lda
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rho_set: the density where you want to evaluate the functional
!!     - deriv_set: place where to store the functional derivatives (they are
!!       added to the derivatives)
!!     - grad_deriv: degree of the derivative that should be evalated,
!!       if positive all the derivatives up to the given degree are evaluated,
!!       if negative only the given degree is calculated
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE xb88_lda_eval(rho_set,deriv_set,grad_deriv,error)
  TYPE(xc_rho_set_type), pointer :: rho_set
  TYPE(xc_derivative_set_type), pointer :: deriv_set
  INTEGER, INTENT(in) :: grad_deriv
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='xb88_lda_eval',&
        routineP=moduleN//':'//routineN
  INTEGER :: i,j,k,npoints
  integer, DIMENSION(:,:), pointer :: bo
  REAL(kind=dp),DIMENSION(:,:,:), POINTER :: dummy,rho, rho_1_3, norm_drho,&
       e_0,e_rho,e_ndrho,e_rho_rho,e_ndrho_rho,&
       e_ndrho_ndrho, e_rho_rho_rho, e_ndrho_rho_rho, e_ndrho_ndrho_rho,&
       e_ndrho_ndrho_ndrho
  TYPE(xc_derivative_type), POINTER :: deriv
  REAL(kind=dp) :: epsilon_rho, epsilon_norm_drho, x, t1,t2,t7,&
       t3,t4,t5,t6,t8,t9,t10,t11,t12,t16,t17,t20,t30,t32,t36,&
       t38,t39,t42,t44,t51,t58,t60,t64,t66,t67,t71,t75,t78,t82,t84,t86,&
       t87,t90,t91,t93,t94,t98,t117,t146,t151,t155,t159,t161,t168,t173,t181
  INTEGER :: stat

  failure=.false.
  nullify(bo)

  CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL xc_rho_set_get(rho_set,rho_1_3=rho_1_3,rho=rho,&
          norm_drho=norm_drho,local_bounds=bo, error=error)
     npoints=(bo(2,1)-bo(1,1))*(bo(2,2)-bo(1,2))*(bo(2,3)-bo(1,3))

     ! meaningful default for the arrays we don't need: let us make compiler
     ! and debugger happy...
     ALLOCATE(dummy(bo(1,1):bo(1,1),bo(1,2):bo(1,2),bo(1,3):bo(1,3)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     
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
             allocate_deriv=.TRUE., error=error)
        call xc_derivative_get(deriv,deriv_data=e_0,error=error)
     END IF
     IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
             allocate_deriv=.TRUE.,error=error)
        call xc_derivative_get(deriv,deriv_data=e_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
             allocate_deriv=.TRUE.,error=error)
        call xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)
     END IF
     IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
             allocate_deriv=.TRUE.,error=error)
        call xc_derivative_get(deriv,deriv_data=e_rho_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rho)",&
             allocate_deriv=.true.,error=error)
        call xc_derivative_get(deriv,deriv_data=e_ndrho_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drho)(norm_drho)", allocate_deriv=.TRUE.,error=error)
        call xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho,error=error)
     END IF
     IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)(rho)",&
             allocate_deriv=.TRUE.,error=error)
        call xc_derivative_get(deriv,deriv_data=e_rho_rho_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drho)(rho)(rho)",allocate_deriv=.true.,error=error)
        call xc_derivative_get(deriv,deriv_data=e_ndrho_rho_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drho)(norm_drho)(rho)",allocate_deriv=.true.,error=error)
        call xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drho)(norm_drho)(norm_drho)", allocate_deriv=.TRUE.,&
             error=error)
        call xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho,error=error)
     END IF
     IF (grad_deriv>3.OR.grad_deriv<-3) THEN
        call cp_unimplemented_error(fromWhere=routineP, &
             message="derivatives bigger than 3 not implemented", &
             error=error, error_level=cp_failure_level)
     END IF

     CALL xb88_lda_calc(rho=rho, rho_1_3=rho_1_3, norm_drho=norm_drho,&
          e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,e_rho_rho=e_rho_rho,&
          e_ndrho_rho=e_ndrho_rho, e_ndrho_ndrho=e_ndrho_ndrho, &
          e_rho_rho_rho=e_rho_rho_rho, e_ndrho_rho_rho=e_ndrho_rho_rho,&
          e_ndrho_ndrho_rho=e_ndrho_ndrho_rho,&
          e_ndrho_ndrho_ndrho=e_ndrho_ndrho_ndrho,grad_deriv=grad_deriv,&
          npoints=npoints,error=error)

     DEALLOCATE(dummy,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
END SUBROUTINE xb88_lda_eval
!***************************************************************************

!!****f* xc_xbecke88/xb88_lda_calc [1.0] *
!!
!!   NAME
!!     xb88_lda_eval
!!
!!   FUNCTION
!!     evaluates the becke 88 exchange functional for lda
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rho_set: the density where you want to evaluate the functional
!!     - deriv_set: place where to store the functional derivatives (they are
!!       added to the derivatives)
!!     - grad_deriv: degree of the derivative that should be evalated,
!!       if positive all the derivatives up to the given degree are evaluated,
!!       if negative only the given degree is calculated
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE xb88_lda_calc(rho, rho_1_3, norm_drho,&
     e_0,e_rho,e_ndrho,e_rho_rho,e_ndrho_rho,&
     e_ndrho_ndrho, e_rho_rho_rho, e_ndrho_rho_rho, e_ndrho_ndrho_rho,&
     e_ndrho_ndrho_ndrho,grad_deriv,npoints,error)
  REAL(kind=dp),DIMENSION(*), intent(in) :: rho, rho_1_3, norm_drho
  REAL(kind=dp), DIMENSION(*), intent(inout) :: e_0,e_rho,e_ndrho,&
       e_rho_rho,e_ndrho_rho,&
       e_ndrho_ndrho, e_rho_rho_rho, e_ndrho_rho_rho, e_ndrho_ndrho_rho,&
       e_ndrho_ndrho_ndrho
  INTEGER, INTENT(in) :: grad_deriv, npoints
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='xb88_lda_calc',&
        routineP=moduleN//':'//routineN
  INTEGER :: ii
  REAL(kind=dp) :: epsilon_rho, epsilon_norm_drho, x, t1,t2,t7,&
       t3,t4,t5,t6,t8,t9,t10,t11,t12,t16,t17,t20,t30,t32,t36,&
       t38,t39,t42,t44,t51,t58,t60,t64,t66,t67,t71,t75,t78,t82,t84,t86,&
       t87,t90,t91,t93,t94,t98,t117,t146,t151,t155,t159,t161,t168,t173,t181

  REAL(kind=dp), PARAMETER :: t_2_3=2.0_dp/3.0_dp,&
       t_4_3=4.0_dp/3.0_dp, t_1_9=1.0_dp/9.0_dp,&
       t_2_9=2.0_dp/9.0_dp, t_4_27=4.0_dp/27.0_dp


  failure=.false.
  t1 = 2.0_dp ** (0.1e1_dp / 0.3e1_dp)
  t2 = t1 ** 2
  t7 = beta * t1
  t3 = 1.5_dp*(0.75*pi)**(1.0_dp/3.0_dp) * t2
  !FM$omp parallel do defaults(none),&
  !FM$omp   shared(rho, rho_1_3, norm_drho,&
  !FM$omp     e_0,e_rho,e_ndrho,e_rho_rho,e_ndrho_rho,&
  !FM$omp     e_ndrho_ndrho, e_rho_rho_rho, e_ndrho_rho_rho, e_ndrho_ndrho_rho,&
  !FM$omp     e_ndrho_ndrho_ndrho,t1,t2,t3,t7),&
  !FM$omp   private(ii,t4,t5,t6,t8,t9,t10,t11,t12,t16,t17,t20,t30,t32,t36,&
  !FM$omp     t38,t39,t42,t44,t51,t58,t60,t64,t66,t67,t71,t75,t78,t82,t84,t86,&
  !FM$omp     t87,t90,t91,t93,t94,t98,t117,t146,t151,t155,t159,t161,t168,t173,t181)
  DO ii=1,npoints
     IF (rho(ii)> epsilon_rho .AND.&
          norm_drho(ii)>epsilon_norm_drho) THEN
        t4 = rho_1_3(ii) * rho(ii)
        x = norm_drho(ii)/t4
        t5 = beta * t2
        t6 = x ** 2
        t8 = t1 * x
        t9 = t8 + 0.1e1_dp
        t10 = SQRT(t9)
        t11 = t8 + t10
        t12 = LOG(t11)
        t16 = 0.1e1_dp + 0.6e1_dp * t7 * x * t12
        t17 = 0.1e1_dp / t16
        t20 = 0.1e1_dp - t5 * t6 * t17

        IF (grad_deriv>=0) THEN
           e_0(ii) = e_0(ii)+t3 * t4 * t20 *0.5_dp
        END IF
        
        t30 = t1 + t1 / (t10 * 0.2e1_dp)
        t32 = 0.1e1_dp / t11
        t36 = 0.6e1_dp * t7 *(t12 + x * t30 * t32)
        t38 = t16 ** 2
        t39 = 0.1e1_dp / t38
        t42 = -0.2e1_dp * t5 * x * t17 + t5 * t6 * t36 * t39
        t44 = t3 * rho_1_3(ii) * t42
        
        IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
           e_rho(ii) = e_rho(ii)+t_2_3 * (-t44 + t3 * rho_1_3(ii) * t20)
           e_ndrho(ii) = e_ndrho(ii)+t3 * t42 *0.5_dp
        END IF
        
        t51 = x * t39
        t58 = beta * x
        t60 = 0.1e1_dp / (t10 * t9)
        t64 = t30 ** 2
        t66 = t11 ** 2
        t67 = 0.1e1_dp / t66
        t71 = 0.12e2_dp * t7 * t30 * t32 - 0.3e1_dp * t58 * t60 *&
             t32 - 0.6e1_dp * t7 * x * t64 * t67
        t75 = t36 ** 2
        t78 = 0.1e1_dp / (t38 * t16)
        t82 = -0.2e1_dp * t5 * t17 + 0.4e1_dp * t5 * t51 * t36 + t5 *&
             t6 * t71 * t39 - 0.2e1_dp * t5 * t6 * t75 * t78
        t84 = t3 * rho_1_3(ii) * t82
        t86 = t_2_3*(- t84 + t44)
        t87 = 0.1e1_dp / rho(ii)
        t90 = rho_1_3(ii) ** 2
        t91 = 0.1e1_dp / t90
        t93 = t3 * t91 * t42
        t94 = 0.2e1_dp *t_1_9 * t93
        t98 = 0.1e1_dp / t4
                    
        IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
           e_rho_rho(ii) = e_rho_rho(ii)-t_4_3*t86 * t87 - t94 + &
                0.2e1_dp*t_1_9 * t3 * t91 * t20
           e_ndrho_rho(ii) = e_ndrho_rho(ii)+t86 * t98
           e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii)+t3 * t82 * t98 * 0.5_dp
        END IF

        t117 = t9 ** 2
        t146 = t38 ** 2
        t151 = t5*(0.6e1_dp * (t39 * t36+t51 * t71+ t6 *(- t71 * t78 * &
             t36+ t75 * t36 / t146)) &
             - 0.12e2_dp * x * t78 * t75 &
             + t6 * (-0.9e1_dp *&
             beta * t60 * t32 - 0.18e2_dp * t7 * t64 * t67 + 0.9e1_dp *&
             t58 /( 0.2e1_dp * t10 * t117) * t32 * t1 + 0.9e1_dp * t58 *&
             t60 * t67 * t30 + 0.12e2_dp * t7 * x * t64 * t30 &
             / (t66 * t11)) * t39)
        t155 = t_2_3 * (- t3 * rho_1_3(ii) * t151 + t84)
        t159 = t3 * t91 * t82
        t161 = -2.0_dp*t_2_3 * t155 * t87 - t_2_9 *&
             t159 + t94
        t168 = rho(ii) ** 2
        t173 = 0.1e1_dp / (t90 * rho(ii))
        t181 = 0.1e1_dp / (t90 * t168)
        
        IF (grad_deriv==-3.or.grad_deriv>=3) THEN
           e_rho_rho_rho(ii) = e_rho_rho_rho(ii) +t_4_3 * (-t161 * t87 -&
                t_2_9 * (-t159 + t93) * t87 + t86 &
                / t168 )+&
                t_4_27 *( t3 * t173 * t42 - t3 * t173 * t20)
           e_ndrho_rho_rho(ii) = e_ndrho_rho_rho(ii)+t161 * t98
           e_ndrho_ndrho_rho(ii) = e_ndrho_ndrho_rho(ii)+t155 * t181
           e_ndrho_ndrho_ndrho(ii) = e_ndrho_ndrho_ndrho(ii)+t3 * t151 * t181 *&
                0.5_dp
        END IF
     END IF
  END DO
  
END SUBROUTINE xb88_lda_calc
!***************************************************************************

!!****f* xc_xbecke88/xb88_lsd_eval [1.0] *
!!
!!   NAME
!!     xb88_lsd_eval
!!
!!   FUNCTION
!!     evaluates the becke 88 exchange functional for lsd
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
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE xb88_lsd_eval(rho_set,deriv_set,grad_deriv,error)
  TYPE(xc_rho_set_type), pointer :: rho_set
  TYPE(xc_derivative_set_type), pointer :: deriv_set
  INTEGER, INTENT(in) :: grad_deriv
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='xb88_lsd_eval',&
       routineP=moduleN//':'//routineN
  INTEGER :: i,j,k, ispin, stat, npoints
  INTEGER, DIMENSION(:,:),pointer :: bo
  REAL(kind=dp), DIMENSION(:,:,:), POINTER :: dummy,e_0
  TYPE(cp_3d_r_p_type),DIMENSION(2) :: rho, rho_1_3, norm_drho,&
       e_rho,e_ndrho,e_rho_rho,e_ndrho_rho,&
       e_ndrho_ndrho, e_rho_rho_rho, e_ndrho_rho_rho, e_ndrho_ndrho_rho,&
       e_ndrho_ndrho_ndrho
  TYPE(xc_derivative_type), POINTER :: deriv
  REAL(kind=dp) :: epsilon_rho, epsilon_norm_drho, x,c, &
       t1,t2,t3,t4,t5,t6,t7,t8,t9,t12,t13,t15,t16,t19,t20,t24,t25,t29,t30,t32,&
       t33,t41,t42,t43,t50,t54,t55,t56,t60,t61,t63,t64,t66,t67,t70,t72,t73,&
       t74,t77,t86,t102,t121,t124,t127,t129,t136,t142,t148
  REAL(kind=dp) :: t_2_3=2.0_dp/3.0_dp, t_1_9=1.0_dp/9.0_dp,t_4_9=4.0_dp/9.0_dp,&
       t_4_27=4.0_dp/27.0_dp, t_8_27=8.0_dp/27.0_dp,t_4_3=4.0_dp/3.0_dp


  failure=.false.
  NULLIFY(deriv, bo)

  CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL xc_rho_set_get(rho_set,rhoa_1_3=rho_1_3(1)%array,&
          rhob_1_3=rho_1_3(2)%array,rhoa=rho(1)%array,&
          rhob=rho(2)%array,norm_drhoa=norm_drho(1)%array, &
          norm_drhob=norm_drho(2)%array,rho_cutoff=epsilon_rho,&
          drho_cutoff=epsilon_norm_drho, local_bounds=bo, error=error)
     npoints=(bo(2,1)-bo(1,1))*(bo(2,2)-bo(1,2))*(bo(2,3)-bo(1,3))

     ! meaningful default for the arrays we don't need: let us make compiler
     ! and debugger happy...
     ALLOCATE(dummy(bo(1,1):bo(1,1),bo(1,2):bo(1,2),bo(1,3):bo(1,3)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

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
             allocate_deriv=.TRUE., error=error)
        CALL xc_derivative_get(deriv, deriv_data=e_0,error=error)
     END IF
     IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rho(1)%array,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rho(2)%array,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho(1)%array,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho(1)%array,error=error)
     END IF
     IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rho_rho(1)%array,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rho_rho(2)%array,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho(1)%array,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho(2)%array,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drhoa)(norm_drhoa)", allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho(1)%array,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drhob)(norm_drhob)", allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho(2)%array,error=error)
     END IF
     IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)(rhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho(1)%array,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)(rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho(2)%array,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drhoa)(rhoa)(rhoa)",allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho_rho(1)%array,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drhob)(rhob)(rhob)",allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho_rho(2)%array,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drhoa)(norm_drhoa)(rhoa)",allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_rho(1)%array,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drhob)(norm_drhob)(rhob)",allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_rho(2)%array,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drhoa)(norm_drhoa)(norm_drhoa)", allocate_deriv=.TRUE.,&
             error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho(1)%array,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drhob)(norm_drhob)(norm_drhob)", allocate_deriv=.TRUE.,&
             error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho(2)%array,error=error)
     END IF

     DO ispin=1,2
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
             grad_deriv=grad_deriv, npoints=npoints,error=error)
     END DO

     DEALLOCATE(dummy,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF

END SUBROUTINE xb88_lsd_eval
!***************************************************************************

!!****f* xc_xbecke88/xb88_lsd_calc [1.0] *
!!
!!   NAME
!!     xb88_lsd_calc
!!
!!   FUNCTION
!!     low level calculation of the becke 88 exchange functional for lsd
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rho_spin: alpha or beta spin density
!!     - rho_1_3_spin: rho_spin**(1./3.)
!!     - norm_drho_spin: || grad rho_spin ||
!!     - e_0: adds to it the local value of the functional
!!     - e_*_spin: derivative of the functional wrt. to the variables
!!       named where the * is. Everything wrt. to the spin of the arguments.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     01.2004 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE xb88_lsd_calc(rho_spin, rho_1_3_spin, norm_drho_spin,e_0,&
     e_rho_spin,e_ndrho_spin,e_rho_rho_spin,e_ndrho_rho_spin,&
     e_ndrho_ndrho_spin, e_rho_rho_rho_spin, e_ndrho_rho_rho_spin,&
     e_ndrho_ndrho_rho_spin,&
     e_ndrho_ndrho_ndrho_spin,grad_deriv,npoints,error)
  REAL(kind=dp),DIMENSION(*), intent(in) :: rho_spin, rho_1_3_spin, norm_drho_spin
  REAL(kind=dp), DIMENSION(*), intent(inout) :: e_0,&
       e_rho_spin,e_ndrho_spin, e_rho_rho_spin, e_ndrho_rho_spin,&
       e_ndrho_ndrho_spin, e_rho_rho_rho_spin, e_ndrho_rho_rho_spin,&
       e_ndrho_ndrho_rho_spin, e_ndrho_ndrho_ndrho_spin
  INTEGER, INTENT(in) :: grad_deriv, npoints
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='xb88_lsd_calc_123',&
       routineP=moduleN//':'//routineN
  INTEGER :: ii
  REAL(kind=dp) :: epsilon_rho, epsilon_norm_drho, x,c, &
       t1,t2,t3,t4,t5,t6,t7,t8,t9,t12,t13,t15,t16,t19,t20,t24,t25,t29,t30,t32,&
       t33,t41,t42,t43,t50,t54,t55,t56,t60,t61,t63,t64,t66,t67,t70,t72,t73,&
       t74,t77,t86,t102,t121,t124,t127,t129,t136,t142,t148
  REAL(kind=dp),PARAMETER :: t_2_3=2.0_dp/3.0_dp, t_1_9=1.0_dp/9.0_dp,&
       t_4_9=4.0_dp/9.0_dp,&
       t_4_27=4.0_dp/27.0_dp, t_8_27=8.0_dp/27.0_dp,t_4_3=4.0_dp/3.0_dp


  failure=.false.
  c=3.0_dp/2.0_dp*(3.0_dp/4.0_dp*Pi)**(1.0_dp/3.0_dp)

  IF (.NOT. failure) THEN

     c=2.0_dp/3.0_dp*(3.0_dp/4.0_dp*pi)**(1.0_dp/3.0_dp)

     !FM$omp parallel do default(none),shared(c,rho_spin, rho_1_3_spin,&
     !FM$omp norm_drho_spin,e_0,&
     !FM$omp e_rho_spin,e_ndrho_spin,e_rho_rho_spin,e_ndrho_rho_spin,&
     !FM$omp e_ndrho_ndrho_spin, e_rho_rho_rho_spin, e_ndrho_rho_rho_spin,&
     !FM$omp e_ndrho_ndrho_rho_spin,&
     !FM$omp e_ndrho_ndrho_ndrho_spin, grad_deriv),&
     !FM$omp private(ii,x,t1,t2,t3,t4,t5,t6,t7,t8,t9,t12,t13,&
     !FM$omp t15,t16,t19,t20,t24,t25,t29,t30,t32,t33,t41,t42,t43,t50,t54,t55,&
     !FM$omp t56,t60,t61,t63,t64,t66,t67,t70,t72,t73,t74,t77,t86,t102,t121,&
     !FM$omp t124,t127,t129,t136,t142,t148)
     DO ii=1,npoints
        IF (rho_spin(ii)> epsilon_rho .AND.&
             norm_drho_spin(ii)>epsilon_norm_drho) THEN
           t1 = rho_1_3_spin(ii) * rho_spin(ii)
           x=norm_drho_spin(ii)/t1
           t3 = x ** 2
           t4 = beta * t3
           t5 = beta * x
           t6 = x + 0.1D1
           t7 = SQRT(t6)
           t8 = x + t7
           t9 = LOG(t8)
           t12 = 0.1D1 + 0.6D1 * t5 * t9
           t13 = 0.1D1 / t12
           t15 = 0.1D1 - t4 * t13

           IF (grad_deriv>=0) THEN
              e_0(ii) = e_0(ii)+c * t1 * t15
           END IF

           t16 = c * rho_1_3_spin(ii)
           t19 = t12 ** 2
           t20 = 0.1D1 / t19
           t24 = 0.1D1 + 0.1D1 / (t7 * 0.2D1)
           t25 = 0.1D1 / t8
           t29 = 0.6D1 * beta * t9 + 0.6D1 * t5 * t24 * t25
           t30 = t20 * t29
           t32 = -0.2D1 * t5 * t13 + t4 * t30
           t33 = t16 * t32

           IF (grad_deriv==-1.or.grad_deriv>=1) THEN
              e_rho_spin(ii) = e_rho_spin(ii)+t_4_3*&
                   (-t33 + t16 * t15)
              e_ndrho_spin(ii) = e_ndrho_spin(ii) +&
                   c * t32
           END IF

           t41 = 0.1D1 / (t19 * t12)
           t42 = t29 ** 2
           t43 = t41 * t42
           t50 = 0.1D1 / (t7 * t6)
           t54 = t24 ** 2
           t55 = t8 ** 2
           t56 = 0.1D1 / t55
           t60 = 0.12D2 * beta * t24 * t25 - 1.5_dp * t5 * t50 * t25 -&
                0.6D1 * t5 * t54 * t56
           t61 = t20 * t60
           t63 = -0.2D1 * beta * t13 + 0.4D1 * t5 * t30 - 0.2D1 * t4 * t43 + &
                t4 * t61
           t64 = t16 * t63
           t66 = t_4_3 * (- t64 + t33)
           t67 = 0.1D1 / rho_spin(ii)
           t70 = rho_1_3_spin(ii) ** 2
           t72 = c / t70
           t73 = t72 * t32
           t74 = t_4_9 * t73
           t77 = 0.1D1 / t1

           IF (grad_deriv==-2.or.grad_deriv>=2) THEN
              e_rho_rho_spin(ii) = e_rho_rho_spin(ii)&
                   -t_4_3*t66 * t67 - t74 + t_4_9 * t72 * t15
              e_ndrho_rho_spin(ii) = e_ndrho_rho_spin(ii)+&
                   t66 * t77
              e_ndrho_ndrho_spin(ii) = e_ndrho_ndrho_spin(ii)+&
                   c * t63 * t77
           END IF

           t86 = t19 ** 2
           t102 = t6 ** 2
           t121 = 0.6D1 * beta * t20 * t29 - 0.12D2 * t5 * t43 + 0.6D1 * t5 *&
                t61 + 0.6D1 * t4 / t86 * t42 * t29 - 0.6D1 * t4 * t41 * t29 * t60 &
                + t4 * t20 * (-4.5_dp * beta * t50 * t25 - 0.18D2 * beta *&
                t54 * t56 + 2.25_dp * t5 /( t7 * t102) * t25 + 4.5_dp &
                * t5 * t50 * t56 * t24 + 0.12D2 * t5 * t54 * t24 / t55 / t8)
           t124 = t_4_3*(-t16 * t121 + t64)
           t127 = t72 * t63
           t129 = -t_4_3 * t124 * t67 -t_4_9 * t127 + t74
           t136 = rho_spin(ii) ** 2
           t142 = c / (t70 * rho_spin(ii))
           t148 = 0.1D1 / (t70 * t136)

           IF (grad_deriv==-3.or.grad_deriv>=3) THEN
              e_rho_rho_rho_spin(ii) = e_rho_rho_rho_spin(ii)&
                   -t_4_3*(t129 * t67 + t_4_9 * (- t127 + t73) * t67 + t66 / t136)&
                   + t_4_27 * t142 * t32 - t_8_27 * t142 * t15
              e_ndrho_rho_rho_spin(ii) = e_ndrho_rho_rho_spin(ii)&
                   +t129 * t77
              e_ndrho_ndrho_rho_spin(ii) = e_ndrho_ndrho_rho_spin(ii)&
                   + t124 * t148
              e_ndrho_ndrho_ndrho_spin(ii) = e_ndrho_ndrho_ndrho_spin(ii)&
                   + c * t121 * t148
           END IF
        END IF
     END DO

  END IF

END SUBROUTINE xb88_lsd_calc
!***************************************************************************

END MODULE xc_xbecke88
