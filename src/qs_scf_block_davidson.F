!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief module that contains the algorithms to perform an itrative
!>         diagonalization by the block-Davidson approach
!>         P. Blaha, et al J. Comp. Physics, 229, (2010), 453-460
!>         \Iterative diagonalization in augmented plane wave based 
!>              methods in electronic structure calculations\
!> \par History
!>      05.2011 created [MI]
!> \author MI
! *****************************************************************************
MODULE qs_scf_block_davidson

  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_transpose,&
                                             cp_fm_triangular_invert
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_diag,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_to_fm_submat,&
                                             cp_fm_type,&
                                             cp_fm_vectorsnorm
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE preconditioner,                  ONLY: apply_full_all
  USE preconditioner_types,            ONLY: preconditioner_type
  USE qs_block_davidson_types,         ONLY: davidson_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_block_davidson'

  PUBLIC :: generate_extended_space

CONTAINS

! *****************************************************************************
  SUBROUTINE generate_extended_space(bdav_env,mo_set,matrix_h,matrix_s,output_unit,&
             preconditioner,error)

    TYPE(davidson_type)                      :: bdav_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h, matrix_s
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(preconditioner_type), OPTIONAL, &
      POINTER                                :: preconditioner
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'generate_extended_space', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, imo, istat, iter, &
                                                max_iter, nao, nmo, &
                                                nmo_converged
    LOGICAL                                  :: converged, &
                                                do_apply_preconditioner, &
                                                failure
    REAL(dp)                                 :: max_norm, min_norm, t1, t2
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: ritz_coeff
    REAL(dp), DIMENSION(:), POINTER          :: eigenvalues
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: block_mat, matrix_hc, &
                                                matrix_sc, matrix_tmp, &
                                                matrix_tmp_t, matrix_z, &
                                                mo_coeff

    failure=.FALSE.

    CALL timeset(routineN,handle)

    NULLIFY(block_mat, eigenvalues, fm_struct_tmp, matrix_hc, mo_coeff, &
         matrix_tmp, matrix_tmp_t, matrix_sc, matrix_z, mo_coeff)
    do_apply_preconditioner = .FALSE.
    IF(PRESENT(preconditioner)) do_apply_preconditioner=.TRUE.
    CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff,eigenvalues=eigenvalues,nao=nao,nmo=nmo)
    IF(do_apply_preconditioner) THEN
         max_iter =  bdav_env%max_iter
    ELSE
       max_iter = 1
    END IF

    CALL cp_fm_create(matrix_sc,mo_coeff%matrix_struct,name="sc",error=error)
    CALL cp_fm_create(matrix_hc,mo_coeff%matrix_struct,name="hc",error=error)

    t1 = m_walltime()
    IF (output_unit > 0) THEN
          WRITE(output_unit,"(T15,A,T23,A,T36,A,T49,A,T60,A,/,T8,A)")   &
                " Cycle ", " conv. MOS ", " B2MAX ", " B2MIN ", " Time",  REPEAT("-",60)
      END IF

!    bdav_env%max_iter = 10
    converged=.FALSE.
    DO iter = 1, max_iter
      ! Prepare HC and SC of same type and size of mo_coeff
      CALL cp_dbcsr_sm_fm_multiply(matrix_s,mo_coeff,matrix_sc,nmo,error=error)
      CALL cp_dbcsr_sm_fm_multiply(matrix_h,mo_coeff,matrix_hc,nmo,error=error)
   
      ! compute Ritz values 
      ALLOCATE(ritz_coeff(nmo),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ritz_coeff=0.0_dp
      CALL ritz_coefficients(bdav_env,mo_coeff,matrix_sc,matrix_hc,ritz_coeff,error=error)
   
      ! extended subspace Z = P [H - theta S]C  this ia another matrix of type and size as mo_coeff 
      CALL cp_fm_create(matrix_tmp,mo_coeff%matrix_struct,name="tmp",error=error)
      CALL cp_fm_to_fm(matrix_sc,matrix_tmp,error=error)
      CALL cp_fm_column_scale(matrix_tmp,ritz_coeff)
      CALL cp_fm_scale_and_add(-1.0_dp,matrix_tmp,1.0_dp,matrix_hc,error=error)
   
      ! Check for converged eigenvectors
      CALL cp_fm_vectorsnorm(matrix_tmp,ritz_coeff,error=error)
      nmo_converged = 0
      max_norm = 0.0_dp
      min_norm = 1.e10_dp
      DO imo = 1,nmo
        max_norm = MAX(max_norm,ritz_coeff(imo))
        min_norm = MIN(min_norm,ritz_coeff(imo))
      END DO
      DO  imo = 1,nmo
        IF(ritz_coeff(imo) <= bdav_env%eps_iter ) THEN
            nmo_converged = nmo_converged + 1
        END IF
      END DO
      IF(REAL(nmo_converged,dp)/REAL(nmo,dp)>bdav_env%conv_percent) converged=.TRUE. 
      DEALLOCATE(ritz_coeff,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
   
   
      ! calculate or retrieve (if [H0- l S0]^-1) the preconditioner P: 
      !  this matrix has the full size nao x nao 
      matrix_z => bdav_env%matrix_z
      IF(do_apply_preconditioner) THEN
        CALL apply_full_all(preconditioner,matrix_tmp,matrix_z,error=error)
      ELSE
       CALL cp_fm_to_fm(matrix_tmp,matrix_z,error=error)
      END IF
      CALL cp_fm_release(matrix_tmp,error=error)
   
      CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nmo,ncol_global=nmo, &
                               context=mo_coeff%matrix_struct%context, &
                               para_env=mo_coeff%matrix_struct%para_env,error=error)
      CALL cp_fm_create(matrix_tmp,fm_struct_tmp,name="matrix_tmp",error=error)
      CALL cp_fm_create(matrix_tmp_t,fm_struct_tmp,name="matrix_tmp_t",error=error)
      CALL cp_fm_struct_release(fm_struct_tmp,error=error)
   
      ! compute the bottom left  ZSC (top right is transpose)
      CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,matrix_z,matrix_sc,0.0_dp,matrix_tmp,error=error)
      ! set the bottom left part of S[C,Z] block matrix  ZSC
      block_mat => bdav_env%S_block_mat
      CALL cp_fm_to_fm_submat(matrix_tmp,block_mat,nmo,nmo,1,1,1+nmo,1,error=error)
      CALL cp_fm_transpose(matrix_tmp,matrix_tmp_t,error=error)
      CALL cp_fm_to_fm_submat(matrix_tmp_t,block_mat,nmo,nmo,1,1,1,1+nmo,error=error)
   
      ! compute the bottom left  ZHC (top right is transpose)
      CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,matrix_z,matrix_hc,0.0_dp,matrix_tmp,error=error)
      ! set the bottom left part of S[C,Z] block matrix  ZHC
      block_mat => bdav_env%H_block_mat
      CALL cp_fm_to_fm_submat(matrix_tmp,block_mat,nmo,nmo,1,1,1+nmo,1,error=error)
      CALL cp_fm_transpose(matrix_tmp,matrix_tmp_t,error=error)
      CALL cp_fm_to_fm_submat(matrix_tmp_t,block_mat,nmo,nmo,1,1,1,1+nmo,error=error)
      CALL cp_fm_release(matrix_tmp_t,error=error)
   
      ! (reuse matrix_sc and matrix_hc to computr HZ and SZ)
      CALL cp_dbcsr_sm_fm_multiply(matrix_s,matrix_z,matrix_sc,nmo,error=error)
      CALL cp_dbcsr_sm_fm_multiply(matrix_h,matrix_z,matrix_hc,nmo,error=error)
   
      ! compute the bottom right  ZSZ 
      CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,matrix_z,matrix_sc,0.0_dp,matrix_tmp,error=error)
      ! set the bottom right part of S[C,Z] block matrix  ZSZ
      block_mat => bdav_env%S_block_mat
      CALL cp_fm_to_fm_submat(matrix_tmp,block_mat,nmo,nmo,1,1,1+nmo,1+nmo,error=error)
   
      ! compute the bottom right  ZHZ 
      CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,matrix_z,matrix_hc,0.0_dp,matrix_tmp,error=error)
      ! set the bottom right part of S[C,Z] block matrix  ZHZ
      block_mat => bdav_env%H_block_mat
      CALL cp_fm_to_fm_submat(matrix_tmp,block_mat,nmo,nmo,1,1,1+nmo,1+nmo,error=error)
   
      ! solution of the reduced eigenvalues problem
      CALL reduce_extended_space(bdav_env,error=error)
   
      block_mat => bdav_env%H_block_vec
      CALL cp_fm_to_fm_submat(block_mat,matrix_tmp,nmo,nmo,1,1,1,1,error=error)
      CALL cp_fm_gemm('N','N',nao,nmo,nmo,1.0_dp,mo_coeff,matrix_tmp,0.0_dp,matrix_hc,error=error)
      CALL cp_fm_to_fm_submat(block_mat,matrix_tmp,nmo,nmo,1+nmo,1,1,1,error=error)
      CALL cp_fm_gemm('N','N',nao,nmo,nmo,1.0_dp,matrix_z,matrix_tmp,1.0_dp,matrix_hc,error=error)
      CALL cp_fm_release(matrix_tmp,error=error)
   
      ! in case some vector are already converged only a subset of vectors are copied in the MOS
      IF (.FALSE.) THEN
   
      ELSE
        CALL cp_fm_to_fm(matrix_hc,mo_coeff,error=error)
        eigenvalues(1:nmo) = bdav_env%evals(1:nmo)
      END IF
      t2 = m_walltime()
      IF (output_unit > 0) THEN
         WRITE(output_unit,'(T16,I5,T24,I6,T33,E12.4,2x,E12.4,T60,F8.3)')&
                 iter, nmo_converged,  max_norm, min_norm, t2-t1
      END IF
      t1=m_walltime()

      IF(converged) THEN
        IF (output_unit > 0)  WRITE(output_unit,*)  " Reached convergence in ", iter, &
             " Davidson iterations"
        EXIT
      END IF

    END DO
    CALL cp_fm_release(matrix_hc,error=error)
    CALL cp_fm_release(matrix_sc,error=error)

    CALL timestop(handle)

  END SUBROUTINE generate_extended_space

  SUBROUTINE ritz_coefficients(bdav_env,mo_coeff,matrix_sc,matrix_hc,ritz_coeff,error)

    TYPE(davidson_type)                      :: bdav_env
    TYPE(cp_fm_type), POINTER                :: mo_coeff, matrix_sc, matrix_hc
    REAL(dp), DIMENSION(:)                   :: ritz_coeff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ritz_coefficients', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, istat, nao, nmo
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: chc_diag, csc_diag
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: block_mat, matrix_tmp

    failure=.FALSE.

    CALL timeset(routineN,handle)

    NULLIFY(block_mat,fm_struct_tmp,matrix_tmp)
    CALL cp_fm_get_info(mo_coeff,nrow_global=nao,ncol_global=nmo,error=error)

    ALLOCATE(csc_diag(nmo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(chc_diag(nmo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! storage matrix of size mos x mos, only the diagonal elements are used
    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nmo,ncol_global=nmo, &
                             context=mo_coeff%matrix_struct%context, &
                             para_env=mo_coeff%matrix_struct%para_env,error=error)
    CALL cp_fm_create(matrix_tmp,fm_struct_tmp,name="matrix_tmp",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    ! since we only use diagonal elements this is a bit of a waste
    ! compute CSC
    CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,mo_coeff,matrix_sc,0.0_dp,matrix_tmp,error=error)
    CALL cp_fm_get_diag(matrix_tmp,csc_diag,error=error)
    ! set the top left part of S[C,Z] block matrix  CSC
    block_mat => bdav_env%S_block_mat
    CALL cp_fm_to_fm_submat(matrix_tmp,block_mat,nmo,nmo,1,1,1,1,error=error)

    ! compute CHC
    CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,mo_coeff,matrix_hc,0.0_dp,matrix_tmp,error=error)
    CALL cp_fm_get_diag(matrix_tmp,chc_diag,error=error)
    ! set the top left part of H[C,Z] block matrix CHC 
    block_mat => bdav_env%H_block_mat
    CALL cp_fm_to_fm_submat(matrix_tmp,block_mat,nmo,nmo,1,1,1,1,error=error)
 
    DO i=1,nmo
      IF(ABS(csc_diag(i))>EPSILON(0.0_dp)) THEN
        ritz_coeff(i) = chc_diag(i)/csc_diag(i)
      END IF
    END DO

    CALL cp_fm_release(matrix_tmp,error=error)

    CALL timestop(handle)

  END SUBROUTINE ritz_coefficients


  SUBROUTINE  reduce_extended_space(bdav_env,error)

    TYPE(davidson_type)                      :: bdav_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reduce_extended_space', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nmox2
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:), POINTER          :: evals
    TYPE(cp_fm_type), POINTER                :: H_block_mat, H_block_vec, &
                                                matrix_tmp, S_block_mat

    failure=.FALSE.

    CALL timeset(routineN,handle)
    NULLIFY(evals, H_block_mat,H_block_vec,matrix_tmp, S_block_mat)

    S_block_mat => bdav_env%S_block_mat
    matrix_tmp  => bdav_env%W_block_mat
    CALL cp_fm_get_info(S_block_mat,nrow_global=nmox2,error=error)
!
!    IF(.false.) THEN
      CALL cp_fm_cholesky_decompose(S_block_mat,error=error)
      CALL cp_fm_triangular_invert(S_block_mat,error=error)
!    ELSE
!      ndep =0
!      CALL cp_fm_power(S_block_mat,matrix_tmp,-0.5_dp,&
!                       1.E-5_dp,ndep,error=error)
!    END IF
!    stop 'invert'
!
    H_block_mat => bdav_env%H_block_mat
    CALL cp_fm_cholesky_restore(H_block_mat,nmox2,S_block_mat,matrix_tmp,&
         "MULTIPLY",pos="RIGHT",error=error)
    CALL cp_fm_cholesky_restore(matrix_tmp,nmox2,S_block_mat,H_block_mat,&
         "MULTIPLY",pos="LEFT",transa="T",error=error)
    evals => bdav_env%evals
    CALL cp_fm_syevd(H_block_mat,matrix_tmp,evals,error=error)
    H_block_vec => bdav_env%H_block_vec
    CALL cp_fm_cholesky_restore(matrix_tmp,nmox2,S_block_mat,H_block_vec,"MULTIPLY",error=error)
!     write(*,*) nmox2, evals(1:nmox2)
!     stop 'restore'
    CALL timestop(handle)


  END SUBROUTINE reduce_extended_space


END MODULE qs_scf_block_davidson
