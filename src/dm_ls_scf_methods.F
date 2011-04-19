!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief lower level routines for linear scaling SCF
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE dm_ls_scf_methods
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_add_on_diag, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_finalize, cp_dbcsr_frobenius_norm, cp_dbcsr_init, &
       cp_dbcsr_multiply, cp_dbcsr_release, cp_dbcsr_scale, cp_dbcsr_set, &
       cp_dbcsr_trace
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE dbcsr_block_access,              ONLY: dbcsr_put_block
  USE dbcsr_data_methods,              ONLY: dbcsr_data_clear_pointer,&
                                             dbcsr_data_init,&
                                             dbcsr_data_new,&
                                             dbcsr_data_release
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_methods,                   ONLY: dbcsr_get_data_type
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_iterator,&
                                             dbcsr_type_no_symmetry
  USE dm_ls_scf_types,                 ONLY: ls_cluster_atomic,&
                                             ls_mstruct_type
  USE f77_blas
  USE input_constants,                 ONLY: ls_s_preconditioner_atomic,&
                                             ls_s_preconditioner_molecular,&
                                             ls_s_preconditioner_none
  USE iterate_matrix,                  ONLY: matrix_sign_Newton_Schulz,&
                                             matrix_sqrt_Newton_Schulz
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dm_ls_scf_methods'

  PUBLIC density_matrix_sign, density_matrix_sign_fixed_mu
  PUBLIC compute_matrix_preconditioner, apply_matrix_preconditioner

CONTAINS

! *****************************************************************************
!> \brief compute for a block positive definite matrix s (bs)
!>        the sqrt(bs) and inv(sqrt(bs))
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE compute_matrix_preconditioner(matrix_s,preconditioner_type, ls_mstruct, &
             matrix_bs_sqrt,matrix_bs_sqrt_inv,threshold,error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_s
    INTEGER                                  :: preconditioner_type
    TYPE(ls_mstruct_type)                    :: ls_mstruct
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_bs_sqrt, &
                                                matrix_bs_sqrt_inv
    REAL(KIND=dp)                            :: threshold
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'compute_matrix_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_col, iblock_row
    LOGICAL                                  :: block_needed, tr
    TYPE(cp_dbcsr_type)                      :: matrix_bs
    TYPE(dbcsr_data_obj)                     :: block_any
    TYPE(dbcsr_iterator)                     :: iter

    CALL timeset(routineN,handle)

    ! first generate a block diagonal copy of s
    CALL cp_dbcsr_init(matrix_bs,error=error)
    CALL cp_dbcsr_create(matrix_bs,template=matrix_s,error=error)

    CALL dbcsr_data_init (block_any)
    CALL dbcsr_data_new(block_any, dbcsr_get_data_type(matrix_s%matrix))

    SELECT CASE(preconditioner_type)
    CASE(ls_s_preconditioner_none)
    CASE(ls_s_preconditioner_atomic, ls_s_preconditioner_molecular)
       CALL dbcsr_iterator_start(iter, matrix_s%matrix)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          !CALL cp_dbcsr_iterator_next_block(iter, iblock_row, iblock_col, block_s, blk)
          CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, block_any, tr)

          ! do we need the block ?
          ! this depends on the preconditioner, but also the matrix clustering method employed
          ! for a clustered matrix, right now, we assume that atomic and molecular preconditioners
          ! are actually the same, and only require that the diagonal blocks (clustered) are present

          block_needed=.FALSE.

          IF (iblock_row==iblock_col) THEN
              block_needed=.TRUE.
          ELSE
              IF (preconditioner_type==ls_s_preconditioner_molecular .AND. &
                  ls_mstruct%cluster_type==ls_cluster_atomic) THEN
                 IF (ls_mstruct%atom_to_molecule(iblock_row)==ls_mstruct%atom_to_molecule(iblock_col)) block_needed=.TRUE.
              ENDIF
          ENDIF
          
          ! add it 
          IF (block_needed) THEN
             !CALL cp_dbcsr_put_block(matrix=matrix_bs,&
             !                   row=iblock_row, col=iblock_col, block=block_s)
             CALL dbcsr_put_block(matrix=matrix_bs%matrix,&
                                row=iblock_row, col=iblock_col, block=block_any)

          ENDIF

       ENDDO
       CALL dbcsr_iterator_stop(iter)
    END SELECT

    CALL dbcsr_data_clear_pointer (block_any)
    CALL dbcsr_data_release (block_any)

    CALL cp_dbcsr_finalize(matrix_bs,error=error)

    SELECT CASE(preconditioner_type)
    CASE(ls_s_preconditioner_none)
       ! for now make it a simple identity matrix
       CALL cp_dbcsr_copy(matrix_bs_sqrt,matrix_bs,error=error)
       CALL cp_dbcsr_set(matrix_bs_sqrt,0.0_dp,error=error)
       CALL cp_dbcsr_add_on_diag(matrix_bs_sqrt,1.0_dp,error=error)

       ! for now make it a simple identity matrix
       CALL cp_dbcsr_copy(matrix_bs_sqrt_inv,matrix_bs,error=error)
       CALL cp_dbcsr_set(matrix_bs_sqrt_inv,0.0_dp,error=error)
       CALL cp_dbcsr_add_on_diag(matrix_bs_sqrt_inv,1.0_dp,error=error)
    CASE(ls_s_preconditioner_atomic, ls_s_preconditioner_molecular)
       CALL cp_dbcsr_copy(matrix_bs_sqrt,matrix_bs,error=error)
       CALL cp_dbcsr_copy(matrix_bs_sqrt_inv,matrix_bs,error=error)
       ! XXXXXXXXXXX 
       ! XXXXXXXXXXX the threshold here could be done differently,
       ! XXXXXXXXXXX using eps_filter is reducing accuracy for no good reason, this is cheap
       ! XXXXXXXXXXX 
       CALL matrix_sqrt_Newton_Schulz(matrix_bs_sqrt,matrix_bs_sqrt_inv,matrix_bs,&
                                      threshold=MIN(threshold,1.0E-10_dp),error=error)
    END SELECT

    CALL cp_dbcsr_release(matrix_bs,error=error)

    CALL timestop(handle)

  END SUBROUTINE compute_matrix_preconditioner

! *****************************************************************************
!> \brief apply a preconditioner either
!>        forward (precondition)            inv(sqrt(bs)) * A * inv(sqrt(bs))
!>        backward (restore to old form)        sqrt(bs)  * A * sqrt(bs)
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE apply_matrix_preconditioner(matrix,direction,matrix_bs_sqrt,matrix_bs_sqrt_inv,error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    CHARACTER(LEN=*)                         :: direction
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_bs_sqrt, &
                                                matrix_bs_sqrt_inv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_matrix_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_type)                      :: matrix_tmp

    CALL timeset(routineN,handle)
    failure=.FALSE.
    CALL cp_dbcsr_init(matrix_tmp,error=error)
    CALL cp_dbcsr_create(matrix_tmp,template=matrix,matrix_type=dbcsr_type_no_symmetry,error=error)

    SELECT CASE(direction)
    CASE("forward")
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix, matrix_bs_sqrt_inv,&
                                0.0_dp, matrix_tmp, error=error)
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_bs_sqrt_inv, matrix_tmp,&
                                0.0_dp, matrix, error=error)
    CASE("backward")
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix, matrix_bs_sqrt,&
                                0.0_dp, matrix_tmp, error=error)
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_bs_sqrt, matrix_tmp,&
                                0.0_dp, matrix, error=error)
    CASE DEFAULT
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    CALL cp_dbcsr_release(matrix_tmp,error=error)

    CALL timestop(handle)

  END SUBROUTINE apply_matrix_preconditioner

! *****************************************************************************
!> \brief compute the density matrix with a trace that is close to nelectron.
!>        take a mu as input, and improve by bisection as needed.
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE density_matrix_sign(matrix_p,mu,fixed_mu,matrix_ks,matrix_s,matrix_s_inv,nelectron,threshold,error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_p
    REAL(KIND=dp), INTENT(INOUT)             :: mu
    LOGICAL                                  :: fixed_mu
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_ks, matrix_s, &
                                                matrix_s_inv
    INTEGER, INTENT(IN)                      :: nelectron
    REAL(KIND=dp), INTENT(IN)                :: threshold
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'density_matrix_sign', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: initial_increment = 0.01_dp

    INTEGER                                  :: handle, iter, unit_nr
    LOGICAL                                  :: has_mu_high, has_mu_low
    REAL(KIND=dp)                            :: increment, mu_high, mu_low, &
                                                trace
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    increment=initial_increment

    has_mu_low=.FALSE.
    has_mu_high=.FALSE.

    ! bisect if both bounds are known, otherwise find the bounds with a linear search
    DO iter=1,30
       IF (has_mu_low .AND. has_mu_high) THEN
          mu=(mu_low+mu_high)/2
          IF (ABS(mu_high-mu_low)<threshold) EXIT
       ENDIF

       CALL density_matrix_sign_fixed_mu(matrix_p,trace,mu,matrix_ks,matrix_s,matrix_s_inv,threshold,error)
       IF (unit_nr>0) WRITE(unit_nr,'(T2,A,I2,1X,F13.9,1X,F15.9)') &
             "Density matrix:  iter, mu, trace error: ", iter, mu, trace-nelectron
    
       ! OK, we can skip early if we are as close as possible to the exact result
       ! smaller differences should be considered 'noise' 
       IF (ABS(trace-nelectron)<0.5_dp .OR. fixed_mu) EXIT
   
       IF (trace<nelectron) THEN
          mu_low=mu
          mu=mu+increment
          has_mu_low=.TRUE.
          increment=increment*2
       ELSE
          mu_high=mu
          mu=mu-increment
          has_mu_high=.TRUE.
          increment=increment*2
       ENDIF
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE density_matrix_sign

! *****************************************************************************
!> \brief for a fixed mu, compute the corresponding density matrix and its trace
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE density_matrix_sign_fixed_mu(matrix_p,trace,mu,matrix_ks,matrix_s,matrix_s_inv,threshold,error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_p
    REAL(KIND=dp), INTENT(OUT)               :: trace
    REAL(KIND=dp), INTENT(INOUT)             :: mu
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_ks, matrix_s, &
                                                matrix_s_inv
    REAL(KIND=dp), INTENT(IN)                :: threshold
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'density_matrix_sign_fixed_mu', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, unit_nr
    REAL(KIND=dp)                            :: frob_matrix
    TYPE(cp_dbcsr_type)                      :: matrix_p_ud, matrix_sign, &
                                                matrix_sinv_ks, matrix_tmp
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! get inv(S)*H-I*mu
    CALL cp_dbcsr_init(matrix_sinv_ks,error=error)
    CALL cp_dbcsr_create(matrix_sinv_ks,template=matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_s_inv, matrix_ks,&
                           0.0_dp, matrix_sinv_ks, filter_eps=threshold,error=error)
    CALL cp_dbcsr_add_on_diag(matrix_sinv_ks,-mu,error=error)

    ! compute sign(inv(S)*H-I*mu)
    CALL cp_dbcsr_init(matrix_sign,error=error)
    CALL cp_dbcsr_create(matrix_sign,template=matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL matrix_sign_Newton_Schulz(matrix_sign,matrix_sinv_ks,threshold,error) 
    CALL cp_dbcsr_release(matrix_sinv_ks,error=error)

    ! now construct the density matrix PS=0.5*(I-sign(inv(S)H-I*mu))
    CALL cp_dbcsr_init(matrix_p_ud,error=error)
    CALL cp_dbcsr_create(matrix_p_ud,template=matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_copy(matrix_p_ud,matrix_sign,error=error)
    CALL cp_dbcsr_scale(matrix_p_ud,-0.5_dp,error=error)
    CALL cp_dbcsr_add_on_diag(matrix_p_ud,0.5_dp,error=error)
    CALL cp_dbcsr_release(matrix_sign,error=error)

    ! we now have PS, lets get its trace
    CALL cp_dbcsr_trace(matrix_p_ud,trace,error=error)

    ! we can also check it is idempotent PS*PS=PS
    CALL cp_dbcsr_init(matrix_tmp,error=error)
    CALL cp_dbcsr_create(matrix_tmp,template=matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_p_ud, matrix_p_ud,&
                           0.0_dp, matrix_tmp, filter_eps=threshold,error=error)
    CALL cp_dbcsr_add(matrix_tmp, matrix_p_ud, 1.0_dp, -1.0_dp, error=error)
    frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp)
    CALL cp_dbcsr_release(matrix_tmp,error=error)
    IF (unit_nr>0) WRITE(unit_nr,'(T2,A,F20.12)') "Deviation from idempotency: ", frob_matrix

    ! get P=PS*inv(S)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_p_ud, matrix_s_inv,&
                           0.0_dp, matrix_p, filter_eps=threshold,error=error)
    CALL cp_dbcsr_release(matrix_p_ud,error=error)

    CALL timestop(handle)

  END SUBROUTINE density_matrix_sign_fixed_mu

END MODULE dm_ls_scf_methods


