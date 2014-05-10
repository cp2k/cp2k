!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief density matrix optimization using exponential transformations
!> \par History
!>       2012.05 created [Florian Schiffmann]
!> \author Florian Schiffmann
! *****************************************************************************

MODULE dm_ls_scf_curvy
  USE bibliography,                    ONLY: Shao2003,&
                                             cite_reference
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_filter, &
       cp_dbcsr_frobenius_norm, cp_dbcsr_init, cp_dbcsr_multiply, &
       cp_dbcsr_norm, cp_dbcsr_release, cp_dbcsr_scale, cp_dbcsr_set, &
       cp_dbcsr_trace, cp_dbcsr_transposed, cp_dbcsr_type
  USE dbcsr_api,                       ONLY: dbcsr_type_no_symmetry
  USE dm_ls_scf_types,                 ONLY: ls_scf_curvy_type,&
                                             ls_scf_env_type
  USE f77_blas
  USE input_constants,                 ONLY: ls_scf_line_search_3point,&
                                             ls_scf_line_search_3point_2d
  USE iterate_matrix,                  ONLY: purify_mcweeny
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE mathconstants,                   ONLY: ifac
  USE mathlib,                         ONLY: invmat
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dm_ls_scf_curvy'

  PUBLIC :: dm_ls_curvy_optimization, deallocate_curvy_data

CONTAINS

  ! *****************************************************************************
  !> \brief driver routine for Head-Gordon curvy step approach  
  !> \par History
  !>       2012.05 created [Florian Schiffmann]
  !> \author Florian Schiffmann
  ! *****************************************************************************

  SUBROUTINE dm_ls_curvy_optimization(ls_scf_env,energy,check_conv,error)
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    REAL(KIND=dp)                            :: energy
    LOGICAL                                  :: check_conv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dm_ls_curvy_optimization', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, lsstep
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    CALL cite_reference(Shao2003)

! Check that the transformation matrices into orthonormal basis are present
    failure=.FALSE.
    CPPrecondition(ls_scf_env%matrix_s_sqrt_inv%ref_count/=0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ls_scf_env%matrix_s_sqrt%ref_count/=0,cp_failure_level,routineP,error,failure)

! Upon first call initialize all matrices needed curing optimization
! In addtion transform P into orthonormal basis. Will be scaled by 0.5 in closed shell case
! Only to be done once as it will be stored and reused afterwards
! TRS4 might give nondiempotent P therefore McWeeny purification is applied on initial P

    IF(.NOT.ALLOCATED(ls_scf_env%curvy_data%matrix_dp))THEN
        CALL init_curvy(ls_scf_env%curvy_data,ls_scf_env%matrix_s,ls_scf_env%nspins,error)
        ls_scf_env%curvy_data%line_search_step=1
        
        IF(ls_scf_env%curvy_data%line_search_type==ls_scf_line_search_3point_2d)THEN
           DO i=1,ls_scf_env%nspins
              CALL cp_dbcsr_copy(ls_scf_env%curvy_data%matrix_psave(i,1),&
                                 ls_scf_env%matrix_p(i),error=error)
           END DO
        END IF
        IF(ls_scf_env%nspins==1)CALL cp_dbcsr_scale(ls_scf_env%matrix_p(1),0.5_dp,error=error) 
        CALL transform_matrix_orth(ls_scf_env%matrix_p,ls_scf_env%matrix_s_sqrt,&
                                   ls_scf_env%eps_filter,error)
        CALL purify_mcweeny(ls_scf_env%matrix_p,ls_scf_env%eps_filter,3,error)
        DO i=1,ls_scf_env%nspins
           CALL cp_dbcsr_copy(ls_scf_env%curvy_data%matrix_p(i),ls_scf_env%matrix_p(i),error=error)
        END DO
    END IF
        
    lsstep=ls_scf_env%curvy_data%line_search_step

! If new search direction has to be computed transform H into the orthnormal basis

    IF(ls_scf_env%curvy_data%line_search_step==1)&
        CALL transform_matrix_orth(ls_scf_env%matrix_ks,ls_scf_env%matrix_s_sqrt_inv,&
                                   ls_scf_env%eps_filter,error)

! Set the energies for the line search and make sure to give the correct energy back to scf_main
    ls_scf_env%curvy_data%energies(lsstep)=energy
    IF(lsstep.ne.1)energy=ls_scf_env%curvy_data%energies(1)

! start the optimization by calling the driver routine or simply combine saved P(2D line search)
    IF(lsstep.le.2)THEN
       CALL optimization_step(ls_scf_env%curvy_data,ls_scf_env,error) 
    ELSE IF(lsstep==ls_scf_env%curvy_data%line_search_type)THEN
! line_search type has the value appropriate to the number of energy calculations needed 
       CALL optimization_step(ls_scf_env%curvy_data,ls_scf_env,error)
    ELSE 
       CALL new_p_from_save(ls_scf_env%matrix_p,ls_scf_env%curvy_data%matrix_psave,lsstep,&
                            ls_scf_env%curvy_data%double_step_size,error)
       ls_scf_env%curvy_data%line_search_step=ls_scf_env%curvy_data%line_search_step+1
       CALL timestop(handle)
       RETURN
    END IF
    lsstep=ls_scf_env%curvy_data%line_search_step

! transform new density matrix back into nonorthonormal basis (again scaling might apply)

    CALL transform_matrix_orth(ls_scf_env%matrix_p,ls_scf_env%matrix_s_sqrt_inv,&
                               ls_scf_env%eps_filter,error)
    IF(ls_scf_env%nspins==1)CALL cp_dbcsr_scale(ls_scf_env%matrix_p(1),2.0_dp,error=error)

! P-matrices only need to be stored in case of 2D line search
    IF(lsstep.le.3.AND.ls_scf_env%curvy_data%line_search_type==ls_scf_line_search_3point_2d)THEN
       DO i=1,ls_scf_env%nspins
           CALL cp_dbcsr_copy(ls_scf_env%curvy_data%matrix_psave(i,lsstep),&
                              ls_scf_env%matrix_p(i),error=error)
       END DO
    END IF
    check_conv=lsstep==1   

    CALL timestop(handle)


  END SUBROUTINE dm_ls_curvy_optimization

  ! *****************************************************************************
  !> \brief low level routine for Head-Gordons curvy step approach  
  !>        computes gradients, performs a cg and line search, 
  !>        and evaluates the BCH series to obtain the new P matrix
  !> \par History
  !>       2012.05 created [Florian Schiffmann]
  !> \author Florian Schiffmann
  ! *****************************************************************************

  SUBROUTINE optimization_step(curvy_data,ls_scf_env,error)
    TYPE(ls_scf_curvy_type)                  :: curvy_data
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'optimization_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    REAL(KIND=dp)                            :: filter, step_size(2)

! Upon first line search step compute new search direction and apply CG if required

    CALL timeset(routineN,handle)

    IF(curvy_data%line_search_step==1)THEN     
       curvy_data%step_size=MAXVAL(curvy_data%step_size)
       curvy_data%step_size=MIN(MAX(0.10_dp,0.5_dp*ABS(curvy_data%step_size(1))),0.5_dp)
! Dynamic eps_filter for newton steps
       filter=MAX(ls_scf_env%eps_filter*curvy_data%min_filter,&
                  ls_scf_env%eps_filter*curvy_data%filter_factor)
       CALL compute_direction_newton(curvy_data%matrix_p, ls_scf_env%matrix_ks,&
               curvy_data%matrix_dp,filter, curvy_data%fix_shift, curvy_data%shift,&
               curvy_data%cg_numer, curvy_data%cg_denom, curvy_data%min_shift,error)
       curvy_data%filter_factor=curvy_data%scale_filter*curvy_data%filter_factor
       step_size=curvy_data%step_size
       curvy_data%BCH_saved=0
    ELSE IF(curvy_data%line_search_step==2)THEN
       step_size=curvy_data%step_size
       IF(curvy_data%energies(1)-curvy_data%energies(2).gt.0.0_dp)THEN
          curvy_data%step_size=curvy_data%step_size*2.0_dp
          curvy_data%double_step_size=.TRUE.
       ELSE
          curvy_data%step_size=curvy_data%step_size*0.5_dp
          curvy_data%double_step_size=.FALSE.
       END IF
       step_size=curvy_data%step_size
    ELSE IF(curvy_data%line_search_step==ls_scf_line_search_3point_2d)THEN
       CALL line_search_2d(curvy_data%energies,curvy_data%step_size,error)
       step_size=curvy_data%step_size
    ELSE IF(curvy_data%line_search_step==ls_scf_line_search_3point)THEN
       CALL line_search_3pnt(curvy_data%energies,curvy_data%step_size,error)
       step_size=curvy_data%step_size
    END IF


    CALL update_p_exp (curvy_data%matrix_p,ls_scf_env%matrix_p,curvy_data%matrix_dp,&
              curvy_data%matrix_BCH,ls_scf_env%eps_filter,step_size,curvy_data%BCH_saved,&
              curvy_data%n_bch_hist,error)

! line_search type has the value appropriate to the numeber of energy calculations needed 
    curvy_data%line_search_step=MOD(curvy_data%line_search_step,curvy_data%line_search_type)+1
    IF(curvy_data%line_search_step==1)THEN
       DO ispin=1,SIZE(curvy_data%matrix_p)
          CALL cp_dbcsr_copy(curvy_data%matrix_p(ispin),ls_scf_env%matrix_p(ispin),error=error)
       END DO
    END IF
    CALL timestop(handle)

  END SUBROUTINE optimization_step

  ! *****************************************************************************
  !> \brief Perform a 6pnt-2D line search for spin polarized calculations.  
  !>        Fit a 2D parabolic function to 6 points
  !> \par History
  !>       2012.05 created [Florian Schiffmann]
  !> \author Florian Schiffmann
  ! *****************************************************************************


  SUBROUTINE line_search_2d(energies,step_size,error)
    REAL(KIND=dp)                            :: energies(6), step_size(2)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'line_search_2d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: info, unit_nr
    REAL(KIND=dp)                            :: a, e_pred, param(6), s1, &
                                                s1sq, s2, s2sq, &
                                                sys_lin_eq(6,6), tmp_e, v1, v2
    TYPE(cp_logger_type), POINTER            :: logger

    logger => cp_error_get_logger(error)
    IF(energies(1)-energies(2).lt.0._dp)THEN
       tmp_e=energies(2); energies(2)=energies(3); energies(3)=tmp_e
       step_size=step_size*2.0_dp
    END IF
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF
    s1=0.5_dp*step_size(1); s2=step_size(1); s1sq=s1**2; s2sq=s2**2 
    sys_lin_eq=0.0_dp;sys_lin_eq(:,6)=1.0_dp
    sys_lin_eq(2,1)=s1sq;sys_lin_eq(2,2)=s1sq;sys_lin_eq(2,3)=s1sq;sys_lin_eq(2,4)=s1;sys_lin_eq(2,5)=s1
    sys_lin_eq(3,1)=s2sq;sys_lin_eq(3,2)=s2sq;sys_lin_eq(3,3)=s2sq;sys_lin_eq(3,4)=s2;sys_lin_eq(3,5)=s2
    sys_lin_eq(4,3)=s1sq;sys_lin_eq(4,5)=s1
    sys_lin_eq(5,1)=s1sq;sys_lin_eq(5,4)=s1
    sys_lin_eq(6,3)=s2sq;sys_lin_eq(6,5)=s2
    
    CALL invmat(sys_lin_eq,info,error)
    param=MATMUL(sys_lin_eq,energies)
    v1=(param(2)*param(4))/(2.0_dp*param(1))-param(5)
    v2=-(param(2)**2)/(2.0_dp*param(1))+2.0_dp*param(3)
    step_size(2)=v1/v2
    step_size(1)=(-param(2)*step_size(2)-param(4))/(2.0_dp*param(1))
    IF(step_size(1).lt.0.0_dp)step_size(1)=1.0_dp
    IF(step_size(2).lt.0.0_dp)step_size(2)=1.0_dp
!    step_size(1)=MIN(step_size(1),2.0_dp)
!    step_size(2)=MIN(step_size(2),2.0_dp)
    e_pred=param(1)*step_size(1)**2+param(2)*step_size(1)*step_size(2)+&
           param(3)*step_size(2)**2+param(4)*step_size(1)+param(5)*step_size(2)+param(6)
    IF(unit_nr.gt.0)WRITE(unit_nr,"(t3,a,F10.5,F10.5,A,F20.9)")&
                   " Line Search: Step Size",step_size," Predicted energy",e_pred
    e_pred=param(1)*s1**2+param(2)*s2*s1*0.0_dp+&
          param(3)*s1**2*0.0_dp+param(4)*s1+param(5)*s1*0.0_dp+param(6)
    
  END SUBROUTINE line_search_2d    

  ! *****************************************************************************
  !> \brief Perform a 3pnt line search   
  !> \par History
  !>       2012.05 created [Florian Schiffmann]
  !> \author Florian Schiffmann
  ! *****************************************************************************
  SUBROUTINE line_search_3pnt(energies,step_size,error)
    REAL(KIND=dp)                            :: energies(3), step_size(2)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'line_search_3pnt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    REAL(KIND=dp)                            :: a, b, c, e_pred, min_val, &
                                                step1, tmp, tmp_e
    TYPE(cp_logger_type), POINTER            :: logger

    logger => cp_error_get_logger(error)
    IF(energies(1)-energies(2).lt.0._dp)THEN
       tmp_e=energies(2); energies(2)=energies(3); energies(3)=tmp_e
       step_size=step_size*2.0_dp
    END IF
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF
    step1=0.5_dp*step_size(1)
    c=energies(1)
    a=(energies(3)+c-2.0_dp*energies(2))/(2.0_dp*step1**2)
    b=(energies(2)-c-a*step1**2)/step1
    IF(a.lt.1.0E-12_dp)a=-1.0E-12_dp
    min_val=-b/(2.0_dp*a)
    e_pred=a*min_val**2+b*min_val+c
    tmp=step_size(1)
    IF(e_pred.lt.energies(1).AND.e_pred.lt.energies(2))THEN
       step_size=  MAX(-1.0_dp,&
            MIN(min_val,10_dp*step_size))
    ELSE
       step_size=1.0_dp
    END IF
    e_pred=a*(step_size(1))**2+b*(step_size(1))+c
    IF(unit_nr.gt.0) THEN
      WRITE(unit_nr,"(t3,a,f16.8,a,F20.9)")"Line Search: Step Size",step_size(1)," Predicted energy",e_pred
      CALL m_flush(unit_nr)
    ENDIF
  END SUBROUTINE line_search_3pnt

  ! *****************************************************************************
  !> \brief Get a new search direction. Iterate to obtain a Newton like step 
  !>        Refine with a CG update of the search direction   
  !> \par History
  !>       2012.05 created [Florian Schiffmann]
  !> \author Florian Schiffmann
  ! *****************************************************************************

  SUBROUTINE compute_direction_newton(matrix_p,matrix_ks,matrix_dp,eps_filter,fix_shift,&
                                      curvy_shift,cg_numer,cg_denom,min_shift,error)
    TYPE(cp_dbcsr_type), DIMENSION(:)        :: matrix_p, matrix_ks, matrix_dp
    REAL(KIND=dp)                            :: eps_filter
    LOGICAL                                  :: fix_shift(2)
    REAL(KIND=dp)                            :: curvy_shift(2), cg_numer(2), &
                                                cg_denom(2), min_shift
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_direction_newton', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, ncyc, &
                                                nspin, unit_nr
    LOGICAL                                  :: at_limit
    REAL(KIND=dp)                            :: beta, conv_val, maxel, &
                                                old_conv, shift
    TYPE(cp_dbcsr_type) :: matrix_Ax, matrix_b, matrix_cg, matrix_dp_old, &
      matrix_PKs, matrix_res, matrix_tmp, matrix_tmp1
    TYPE(cp_logger_type), POINTER            :: logger

    logger => cp_error_get_logger(error)

    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF
    CALL timeset(routineN,handle)
    nspin=SIZE(matrix_p)

    CALL cp_dbcsr_init(matrix_PKs,error=error)
    CALL cp_dbcsr_create(matrix_PKs,template=matrix_dp(1),matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_Ax,error=error)
    CALL cp_dbcsr_create(matrix_Ax,template=matrix_dp(1),matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_tmp,error=error)
    CALL cp_dbcsr_create(matrix_tmp,template=matrix_dp(1),matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_tmp1,error=error)
    CALL cp_dbcsr_create(matrix_tmp1,template=matrix_dp(1),matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_res,error=error)
    CALL cp_dbcsr_create(matrix_res,template=matrix_dp(1),matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_cg,error=error)
    CALL cp_dbcsr_create(matrix_cg,template=matrix_dp(1),matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_b,error=error)
    CALL cp_dbcsr_create(matrix_b,template=matrix_dp(1),matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_dp_old,error=error)
    CALL cp_dbcsr_create(matrix_dp_old,template=matrix_dp(1),matrix_type=dbcsr_type_no_symmetry,error=error)

    DO ispin=1,nspin
       CALL cp_dbcsr_copy(matrix_dp_old,matrix_dp(ispin),error=error)

! Precompute some matrices to save work during iterations 
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_p(ispin), matrix_ks(ispin),&
                               0.0_dp, matrix_PKs, filter_eps=eps_filter,error=error)
       CALL cp_dbcsr_transposed(matrix_b,matrix_PKs,error=error)
       CALL cp_dbcsr_copy(matrix_cg,matrix_b,error=error)

! Starting CG with guess 0-matrix gives -2*gradient=[Ks*P-(Ks*P)T] for cg_matrix in second step
       CALL cp_dbcsr_add(matrix_cg,matrix_PKs,2.0_dp,-2.0_dp,error=error)

! Residual matrix in first step=cg matrix. Keep Pks for later use in CG!
       CALL cp_dbcsr_copy(matrix_res,matrix_cg,error=error)

! Precompute -FP-[FP]T which will be used throughout the CG iterations
       CALL cp_dbcsr_add(matrix_b,matrix_PKs, -1.0_dp,-1.0_dp,error=error)

! Setup some values to check convergence and safty checks for eigenvalue shifting
      CALL cp_dbcsr_norm(matrix_res,which_norm=2,norm_scalar=old_conv,error=error)
      old_conv=cp_dbcsr_frobenius_norm(matrix_res)
      shift=MIN(10.0_dp,MAX(min_shift,0.05_dp*old_conv))
      conv_val=MAX(0.010_dp*old_conv,100.0_dp*eps_filter)
      old_conv=100.0_dp
      IF(fix_shift(ispin))THEN
         shift=MAX(min_shift,MIN(10.0_dp,MAX(shift,curvy_shift(ispin)-0.5_dp*curvy_shift(ispin))))
         curvy_shift(ispin)=shift
      END IF

! Begin the real optimization loop
      CALL cp_dbcsr_set(matrix_dp(ispin),0.0_dp,error=error)
      ncyc=10
      DO i=1,ncyc

! One step to compute: -FPD-DPF-DFP-PFD (not obvious but symmetry allows for some tricks)
         CALL commutator_symm(matrix_b,matrix_cg,matrix_Ax,eps_filter,1.0_dp,error)

! Compute the missing bits 2*(FDP+PDF) (again use symmetry to compute as a commutator)
         CALL cp_dbcsr_multiply("N", "N",1.0_dp, matrix_cg,matrix_p(ispin), &
              0.0_dp, matrix_tmp, filter_eps=eps_filter,error=error)
         CALL commutator_symm(matrix_ks(ispin),matrix_tmp,matrix_tmp1,eps_filter,2.0_dp,error)
         CALL cp_dbcsr_add(matrix_Ax,matrix_tmp1,1.0_dp,1.0_dp,error=error)

! Apply the shift and hope it's enough to stabilize the CG iterations
         CALL cp_dbcsr_add(matrix_Ax,matrix_cg,1.0_dp,shift,error=error)
              
         CALL compute_cg_matrices(matrix_Ax,matrix_res,matrix_cg,matrix_dp(ispin),&
              matrix_tmp,eps_filter,at_limit,unit_nr,error)
         CALL cp_dbcsr_filter(matrix_cg,eps_filter,error=error)

! check for convergence of the newton step
         maxel=cp_dbcsr_frobenius_norm(matrix_res) 
         IF(unit_nr.gt.0) THEN
           WRITE(unit_nr,"(T3,A,F12.6)")"Convergence of Newton iteration ",maxel
           CALL m_flush(unit_nr) 
         ENDIF
         at_limit=at_limit.OR.(old_conv/maxel.lt.1.01_dp)
         old_conv=maxel
         IF(i==ncyc.AND.maxel/conv_val.gt.5.0_dp)THEN
             fix_shift(ispin)=.TRUE.
             curvy_shift(ispin)=4.0_dp*shift
         END IF
         IF(maxel.lt.conv_val.OR.at_limit)EXIT
       END DO

! Refine the Newton like search direction with a preconditioned cg update
       CALL cp_dbcsr_transposed(matrix_b,matrix_PKs, error=error)
       !compute b= -2*KsP+2*PKs=-(2*gradient)
       CALL cp_dbcsr_copy(matrix_cg,matrix_b,error=error)
       CALL cp_dbcsr_add(matrix_cg,matrix_PKs,1.0_dp,-1.0_dp,error=error)
       cg_denom(ispin)=cg_numer(ispin)
       CALL cp_dbcsr_trace(matrix_cg,matrix_dp(ispin),cg_numer(ispin),error=error)
       beta=cg_numer(ispin)/MAX(cg_denom(ispin),1.0E-6_dp)
       IF(beta.lt.1.0_dp)THEN
          beta=MAX(0.0_dp,beta)
          CALL cp_dbcsr_add(matrix_dp(ispin),matrix_dp_old,1.0_dp,beta,error=error)  
       END IF
       IF(unit_nr.gt.0)WRITE(unit_nr,"(A)")" "
    END DO

    CALL cp_dbcsr_release(matrix_PKs,error)
    CALL cp_dbcsr_release(matrix_dp_old,error)
    CALL cp_dbcsr_release(matrix_b,error)
    CALL cp_dbcsr_release(matrix_Ax,error)
    CALL cp_dbcsr_release(matrix_tmp,error)
    CALL cp_dbcsr_release(matrix_tmp1,error)
    CALL cp_dbcsr_release(matrix_b,error)
    CALL cp_dbcsr_release(matrix_res,error)
    CALL cp_dbcsr_release(matrix_cg,error)

    IF (unit_nr.gt.0) CALL m_flush(unit_nr)
    CALL timestop(handle)
  END SUBROUTINE compute_direction_newton

              

  ! *****************************************************************************
  !> \brief compute the optimal step size of the current cycle and update the 
  !>        matrices needed to solve the system of linear equations
  !> \par History
  !>       2012.05 created [Florian Schiffmann]
  !> \author Florian Schiffmann
  ! *****************************************************************************

  SUBROUTINE compute_cg_matrices(Ax,res,cg,deltp,tmp,eps_filter,at_limit,unit_nr,error)          
    TYPE(cp_dbcsr_type)                      :: Ax, res, cg, deltp, tmp
    REAL(KIND=dp)                            :: eps_filter
    LOGICAL                                  :: at_limit
    INTEGER                                  :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_cg_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, info
    REAL(KIND=dp)                            :: alpha, beta, devi(3), fac, &
                                                fac1, lin_eq(3,3), new_norm, &
                                                norm_cA, norm_rr, vec(3)

     at_limit=.FALSE.
     CALL cp_dbcsr_trace(res,res,norm_rr,error=error)
     CALL cp_dbcsr_trace(cg,Ax,norm_cA,error=error)
     lin_eq=0.0_dp
     fac=norm_rr/norm_cA
     fac1=fac
! Use a 3point line serach and a fit to a quadratic function to determine optimal step size
     DO i=1,3
        CALL cp_dbcsr_copy(tmp,res,error=error)
        CALL cp_dbcsr_add(tmp,Ax,1.0_dp,-fac,error=error)
        devi(i)=cp_dbcsr_frobenius_norm(tmp)
        lin_eq(i,:)=(/fac**2,fac,1.0_dp/)
        fac=fac1+fac1*((-1)**i)*0.5_dp
     END DO        
     CALL invmat(lin_eq,info,error) 
     vec=MATMUL(lin_eq,devi)
     alpha=-vec(2)/(2.0_dp*vec(1))
     fac=SQRT(norm_rr/(norm_cA*alpha))
!scale the previous matrices to match the step size
     CALL cp_dbcsr_scale(Ax,fac,error=error)
     CALL cp_dbcsr_scale(cg,fac,error=error)
     norm_cA=norm_cA*fac**2

! USe CG to get the new matrices
     alpha=norm_rr/norm_cA
     CALL cp_dbcsr_add(res,Ax,1.0_dp,-alpha,error=error)
     CALL cp_dbcsr_trace(res,res,new_norm,error=error)
     IF(norm_rr.lt.eps_filter*0.001_dp.OR.new_norm.lt.eps_filter*0.001_dp)THEN
        beta=0.0_dp
        at_limit=.TRUE.
     ELSE
        beta=new_norm/norm_rr
        CALL cp_dbcsr_add(deltp,cg,1.0_dp,alpha,error=error)
     END IF
     beta=new_norm/norm_rr
     CALL cp_dbcsr_add(cg,res,beta,1.0_dp,error=error)
     
  END SUBROUTINE compute_cg_matrices

  ! *****************************************************************************
  !> \brief Only for 2D line serach. Use saved P-components to construct new
  !>        test density matrix. Takes care as well, whether step_size
  !>        increased or decreased during 2nd step and combines matrices accordingly
  !> \par History
  !>       2012.05 created [Florian Schiffmann]
  !> \author Florian Schiffmann
  ! *****************************************************************************
     

  SUBROUTINE new_p_from_save(matrix_p,matrix_psave,lsstep,DOUBLE,error)
    TYPE(cp_dbcsr_type), DIMENSION(:)        :: matrix_p
    TYPE(cp_dbcsr_type), DIMENSION(:, :)     :: matrix_psave
    INTEGER                                  :: lsstep
    LOGICAL                                  :: DOUBLE
    TYPE(cp_error_type), INTENT(INOUT)       :: error

     SELECT CASE(lsstep)
     CASE(3)
       CALL cp_dbcsr_copy(matrix_p(1),matrix_psave(1,1),error=error)  
       IF(DOUBLE)THEN
          CALL cp_dbcsr_copy(matrix_p(2),matrix_psave(2,2),error=error)  
       ELSE
          CALL cp_dbcsr_copy(matrix_p(2),matrix_psave(2,3),error=error)  
       END IF   
     CASE(4)
       IF(DOUBLE)THEN
          CALL cp_dbcsr_copy(matrix_p(1),matrix_psave(1,2),error=error) 
       ELSE
          CALL cp_dbcsr_copy(matrix_p(1),matrix_psave(1,3),error=error)  
       END IF
       CALL cp_dbcsr_copy(matrix_p(2),matrix_psave(2,1),error=error)  
     CASE(5)
       CALL cp_dbcsr_copy(matrix_p(1),matrix_psave(1,1),error=error)  
       IF(DOUBLE)THEN
          CALL cp_dbcsr_copy(matrix_p(2),matrix_psave(2,3),error=error)  
       ELSE
          CALL cp_dbcsr_copy(matrix_p(2),matrix_psave(2,2),error=error)  
       END IF
     END SELECT
   
  END SUBROUTINE new_p_from_save

  ! *****************************************************************************
  !> \brief computes a commutator exploiting symmetry RES=k*[A,B]=k*[AB-(AB)T]  
  !> \param eps_filter   filtering threshold for sparse matrices 
  !> \param prefac      prefactor k in above equation 
  !> \par History
  !>       2012.05 created [Florian Schiffmann]
  !> \author Florian Schiffmann
  ! *****************************************************************************
 
  SUBROUTINE commutator_symm(a,b,res,eps_filter,prefac,error)
    TYPE(cp_dbcsr_type)                      :: a, b, res
    REAL(KIND=dp)                            :: eps_filter, prefac
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'commutator_symm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_dbcsr_type)                      :: work

     CALL timeset(routineN,handle)

     CALL cp_dbcsr_init(work,error=error)
     CALL cp_dbcsr_create(work,template=a,matrix_type=dbcsr_type_no_symmetry,error=error)

     CALL cp_dbcsr_multiply("N", "N", prefac, a, b,0.0_dp, res, filter_eps=eps_filter,error=error)
     CALL cp_dbcsr_transposed(work,res,error=error)
     CALL cp_dbcsr_add(res,work,1.0_dp,-1.0_dp,error=error)
   
     CALL cp_dbcsr_release(work,error)
   
     CALL timestop(handle)
  END SUBROUTINE commutator_symm
 
  ! *****************************************************************************
  !> \brief Use the BCH update to get the new idempotent P
  !>        Numerics don't allow for perfect idempotency, therefore a mc weeny 
  !>        step is used to make sure we stay close to the idempotent surface  
  !> \par History
  !>       2012.05 created [Florian Schiffmann]
  !> \author Florian Schiffmann
  ! *****************************************************************************
 

  SUBROUTINE update_p_exp(matrix_p_in,matrix_p_out,matrix_dp,matrix_BCH,threshold,step_size,&
                          BCH_saved,n_bch_hist,error)
    TYPE(cp_dbcsr_type), DIMENSION(:)        :: matrix_p_in, matrix_p_out, &
                                                matrix_dp
    TYPE(cp_dbcsr_type), DIMENSION(:, :)     :: matrix_BCH
    REAL(KIND=dp)                            :: threshold, step_size(2)
    INTEGER                                  :: BCH_saved(2), n_bch_hist
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_p_exp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, nsave, &
                                                nspin, unit_nr
    LOGICAL                                  :: save_BCH
    REAL(KIND=dp)                            :: frob_norm, step_fac
    TYPE(cp_dbcsr_type)                      :: matrix, matrix_tmp
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
      unit_nr=-1
    ENDIF

    CALL cp_dbcsr_init(matrix,error=error)
    CALL cp_dbcsr_create(matrix,template=matrix_p_in(1),matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_tmp,error=error)
    CALL cp_dbcsr_create(matrix_tmp,template=matrix_p_in(1),matrix_type=dbcsr_type_no_symmetry,error=error)
    nspin=SIZE(matrix_p_in)
       
    DO ispin=1,nspin
       step_fac=1.0_dp
       frob_norm=1.0_dp
       nsave=0       
 
       CALL cp_dbcsr_copy(matrix_tmp,matrix_p_in(ispin),error=error)
       CALL cp_dbcsr_copy(matrix_p_out(ispin),matrix_p_in(ispin),error=error)
! If a BCH history is used make good use of it and do a few steps as a copy and scale update of P
! else BCH_saved will be 0 and loop is skipped
       DO i=1,BCH_saved(ispin)
          step_fac=step_fac*step_size(ispin)
          CALL cp_dbcsr_copy(matrix_tmp,matrix_p_out(ispin),error=error)
          CALL cp_dbcsr_add(matrix_p_out(ispin),matrix_BCH(ispin,i),1.0_dp,ifac(i)*step_fac,error=error)
          CALL cp_dbcsr_add(matrix_tmp,matrix_p_out(ispin),1.0_dp,-1.0_dp,error=error)
          frob_norm=cp_dbcsr_frobenius_norm(matrix_tmp)
          IF(unit_nr.gt.0)WRITE(unit_nr,"(t3,a,i3,a,f16.8)")"BCH: step",i," Norm of P_old-Pnew:",frob_norm
          IF(frob_norm.lt.threshold) EXIT
       END DO
       IF(frob_norm.lt.threshold)CYCLE  
 
! If the copy and scale isn't enough compute a few more BCH steps. 20 seems high but except of the first step it will never be close
       save_BCH=BCH_saved(ispin)==0.AND.n_bch_hist.GT.0
       DO i=BCH_saved(ispin)+1,20
          step_fac=step_fac*step_size(ispin)
          !allow for a bit of matrix magic here by exploiting matrix and matrix_tmp
          !matrix_tmp is alway the previous order of the BCH series
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp, matrix_dp(ispin),&
               0.0_dp, matrix, filter_eps=threshold, error=error)

          !(anti)symmetry allows to sum the transposed instead of the full commutator, matrix becomes the latest result

          CALL cp_dbcsr_transposed(matrix_tmp,matrix,error=error)
          CALL cp_dbcsr_add(matrix,matrix_tmp,1.0_dp,1.0_dp,error=error)

          !Finally, add the new BCH order to P, but store the previous one for a convergence check
          CALL cp_dbcsr_copy(matrix_tmp,matrix_p_out(ispin),error=error)
          CALL cp_dbcsr_add(matrix_p_out(ispin),matrix,1.0_dp,ifac(i)*step_fac,error=error)
          IF(save_BCH.AND.i.le.n_bch_hist)THEN
             CALL cp_dbcsr_copy(matrix_BCH(ispin,i),matrix,error=error)
             nsave=i
          END IF

          CALL cp_dbcsr_add(matrix_tmp,matrix_p_out(ispin),1.0_dp,-1.0_dp,error=error)

          !Stop the BCH-series if two successive P's differ by less the threshold
          frob_norm=cp_dbcsr_frobenius_norm(matrix_tmp)
          IF(unit_nr.gt.0)WRITE(unit_nr,"(t3,a,i3,a,f16.8)")"BCH: step",i," Norm of P_old-Pnew:",frob_norm
          IF(frob_norm.lt.threshold) EXIT

          !Copy the latest BCH-matrix on matrix tmp, so we can cycle with all matrices in place
          CALL cp_dbcsr_copy(matrix_tmp,matrix,error=error)
          CALL cp_dbcsr_filter(matrix_tmp,threshold,error=error)
       END DO
       BCH_saved(ispin)=nsave
       IF(unit_nr.gt.0)WRITE(unit_nr,"(A)")" "
    END DO

    CALL purify_mcweeny(matrix_p_out,threshold,1,error)
    IF (unit_nr.gt.0) CALL m_flush(unit_nr)
    CALL cp_dbcsr_release (matrix_tmp,error)
    CALL cp_dbcsr_release (matrix,error)
    CALL timestop(handle)
  END SUBROUTINE update_p_exp

  ! *****************************************************************************
  !> \brief performs a tranformation of a matrix back to/into orthonormal basis  
  !>        in case of P a scaling of 0.5 has to be applied for closed shell case
  !> \param matrix       matrix to be transformed
  !> \param matrix_trafo transformation matrix
  !> \param eps_filter   filtering threshold for sparse matrices 
  !> \par History
  !>       2012.05 created [Florian Schiffmann]
  !> \author Florian Schiffmann
  ! *****************************************************************************

  SUBROUTINE transform_matrix_orth(matrix,matrix_trafo,eps_filter,error)
    TYPE(cp_dbcsr_type), DIMENSION(:)        :: matrix
    TYPE(cp_dbcsr_type)                      :: matrix_trafo
    REAL(KIND=dp)                            :: eps_filter
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'transform_matrix_orth', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    TYPE(cp_dbcsr_type)                      :: matrix_tmp, matrix_work

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_init(matrix_work,error=error)
    CALL cp_dbcsr_create(matrix_work,template=matrix(1),matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_tmp,error=error)
    CALL cp_dbcsr_create(matrix_tmp,template=matrix(1),matrix_type=dbcsr_type_no_symmetry,error=error)

    DO ispin=1,SIZE(matrix)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix(ispin), matrix_trafo,&
            0.0_dp, matrix_work, filter_eps=eps_filter, error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_trafo, matrix_work,&
            0.0_dp, matrix_tmp,filter_eps=eps_filter, error=error)
       ! symmetrize results (this is again needed to make sure everything is stable)
       CALL cp_dbcsr_transposed(matrix_work,matrix_tmp,error=error)
       CALL cp_dbcsr_add(matrix_tmp,matrix_work,0.5_dp,0.5_dp,error=error)
       CALL cp_dbcsr_copy(matrix(ispin),matrix_tmp,error=error)   
    END DO

    CALL cp_dbcsr_release(matrix_tmp,error=error)
    CALL cp_dbcsr_release(matrix_work,error=error)
    CALL timestop(handle)

  END SUBROUTINE
       
  SUBROUTINE deallocate_curvy_data(curvy_data,error)
    TYPE(ls_scf_curvy_type)                  :: curvy_data
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_curvy_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j

    CALL release_dbcsr_array(curvy_data%matrix_dp,error)
    CALL release_dbcsr_array(curvy_data%matrix_p,error)

    IF(ALLOCATED(curvy_data%matrix_psave))THEN
       DO i=1,SIZE(curvy_data%matrix_psave,1)
          DO j=1,3
             CALL cp_dbcsr_release(curvy_data%matrix_psave(i,j),error)  
          END DO
       END DO
       DEALLOCATE(curvy_data%matrix_psave)
    END IF
    IF(ALLOCATED(curvy_data%matrix_BCH))THEN
       DO i=1,SIZE(curvy_data%matrix_BCH,1)
          DO j=1,7
             CALL cp_dbcsr_release(curvy_data%matrix_BCH(i,j),error)  
          END DO
       END DO
       DEALLOCATE(curvy_data%matrix_BCH)  
    END IF
  END SUBROUTINE deallocate_curvy_data

  SUBROUTINE release_dbcsr_array(matrix,error)
    TYPE(cp_dbcsr_type), ALLOCATABLE, &
      DIMENSION(:)                           :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'release_dbcsr_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF(ALLOCATED(matrix))THEN
       DO i=1,SIZE(matrix)
          CALL cp_dbcsr_release(matrix(i),error)
       END DO
       DEALLOCATE(matrix)
    END IF
  END SUBROUTINE release_dbcsr_array



  SUBROUTINE init_curvy(curvy_data,matrix_s,nspins,error)
    TYPE(ls_scf_curvy_type)                  :: curvy_data
    TYPE(cp_dbcsr_type)                      :: matrix_s
    INTEGER                                  :: nspins
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_curvy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, j

    ALLOCATE (curvy_data%matrix_dp(nspins))
    ALLOCATE (curvy_data%matrix_p(nspins))
    DO ispin=1,nspins
       CALL cp_dbcsr_init(curvy_data%matrix_dp(ispin),error=error)
       CALL cp_dbcsr_create(curvy_data%matrix_dp(ispin),template=matrix_s,&
            matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_set(curvy_data%matrix_dp(ispin),0.0_dp,error=error)
       CALL cp_dbcsr_init(curvy_data%matrix_p(ispin),error=error)
       CALL cp_dbcsr_create(curvy_data%matrix_p(ispin),template=matrix_s,&
            matrix_type=dbcsr_type_no_symmetry,error=error)
       curvy_data%fix_shift=.FALSE.
       curvy_data%double_step_size=.TRUE.
       curvy_data%shift=1.0_dp
       curvy_data%BCH_saved=0
       curvy_data%step_size=0.60_dp
       curvy_data%cg_numer=0.00_dp
       curvy_data%cg_denom=0.00_dp
    END DO
    IF(curvy_data%line_search_type==ls_scf_line_search_3point_2d)THEN
       ALLOCATE (curvy_data%matrix_psave(nspins,3))
       DO ispin=1,nspins
          DO j=1,3
             CALL cp_dbcsr_init(curvy_data%matrix_psave(ispin,j),error=error)
             CALL cp_dbcsr_create(curvy_data%matrix_psave(ispin,j),template=matrix_s,&
                                  matrix_type=dbcsr_type_no_symmetry,error=error)
          END DO
       END DO
    END IF
    IF(curvy_data%n_bch_hist.gt.0)THEN
       ALLOCATE (curvy_data%matrix_BCH(nspins,curvy_data%n_bch_hist))
       DO ispin=1,nspins
          DO j=1,curvy_data%n_bch_hist
             CALL cp_dbcsr_init(curvy_data%matrix_BCH(ispin,j),error=error)
             CALL cp_dbcsr_create(curvy_data%matrix_BCH(ispin,j),template=matrix_s,&
                                  matrix_type=dbcsr_type_no_symmetry,error=error)
          END DO
        END DO
     END IF
       
       
  END SUBROUTINE init_curvy

END MODULE dm_ls_scf_curvy
