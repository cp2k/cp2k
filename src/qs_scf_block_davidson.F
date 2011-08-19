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

  USE array_types,                     ONLY: array_i1d_obj,&
                                             array_release
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add, &
       cp_create_bl_distribution, cp_dbcsr_col_block_sizes, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_distribution, cp_dbcsr_distribution_release,&
       cp_dbcsr_get_diag, cp_dbcsr_get_info, cp_dbcsr_init_p, cp_dbcsr_multiply, &
       cp_dbcsr_norm,  cp_dbcsr_release_p, cp_dbcsr_row_block_sizes,&
       cp_dbcsr_scale_by_vector, cp_dbcsr_set
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply, copy_dbcsr_to_fm, copy_fm_to_dbcsr
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
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_to_fm_submat,&
                                             cp_fm_type,&
                                             cp_fm_vectorsnorm
  USE dbcsr_methods,                   ONLY: dbcsr_distribution_mp,&
                                             dbcsr_distribution_new,&
                                             dbcsr_distribution_row_dist,&
                                             dbcsr_mp_npcols,&
                                             dbcsr_mp_nprows
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_type_complex_default,&
                                             dbcsr_type_no_symmetry,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_default, dbcsr_norm_column
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE message_passing,                 ONLY: mp_sum
  USE preconditioner,                  ONLY: apply_preconditioner
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

  PUBLIC :: generate_extended_space, generate_extended_space_sparse

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
                                                nmo_converged, nmo_not_converged
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iconv
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

    Allocate(iconv(nmo), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ritz_coeff(nmo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

!    bdav_env%max_iter = 10
    converged=.FALSE.
    DO iter = 1, max_iter
      ! Prepare HC and SC of same type and size of mo_coeff
      CALL cp_dbcsr_sm_fm_multiply(matrix_s,mo_coeff,matrix_sc,nmo,error=error)
      CALL cp_dbcsr_sm_fm_multiply(matrix_h,mo_coeff,matrix_hc,nmo,error=error)
   
      ! compute Ritz values 
      ritz_coeff=0.0_dp
      CALL ritz_coefficients(bdav_env,mo_coeff,matrix_sc,matrix_hc,ritz_coeff,error=error)

      ! extended subspace Z = P [H - theta S]C  this ia another matrix of type and size as mo_coeff 
      CALL cp_fm_create(matrix_tmp,mo_coeff%matrix_struct,name="tmp",error=error)
      CALL cp_fm_to_fm(matrix_sc,matrix_tmp,error=error)
      CALL cp_fm_column_scale(matrix_tmp,ritz_coeff)
      CALL cp_fm_scale_and_add(-1.0_dp,matrix_tmp,1.0_dp,matrix_hc,error=error)
   
      ! Check for converged eigenvectors
      ritz_coeff=0.0_dp
      CALL cp_fm_vectorsnorm(matrix_tmp,ritz_coeff,error=error)
      nmo_converged = 0
      max_norm = 0.0_dp
      min_norm = 1.e10_dp
      DO imo = 1,nmo
        max_norm = MAX(max_norm,ritz_coeff(imo))
        min_norm = MIN(min_norm,ritz_coeff(imo))
      END DO
      iconv = 0
      DO  imo = 1,nmo
        IF(ritz_coeff(imo) <= bdav_env%eps_iter ) THEN
            nmo_converged = nmo_converged + 1
            iconv(nmo_converged)=imo
!            econv(nmo_converged)=eigenvalues(imo)
!        ELSE
!          ! keep only sequential converged vectors, starting from lowest
!          exit
        END IF
      END DO
      nmo_not_converged=nmo-nmo_converged
!dbg
!      IF(nmo_converged>0) THEN
!        write(*,*) " NMO ", nmo, nmo_converged, nmo_not_converged
!        write(*,*) "iconv",  iconv(1:nmo_converged)
!      END IF
!dbg

      IF(REAL(nmo_converged,dp)/REAL(nmo,dp)>bdav_env%conv_percent) converged=.TRUE. 
   
      ! calculate or retrieve (if [H0- l S0]^-1) the preconditioner P: 
      !  this matrix has the full size nao x nao 
      matrix_z => bdav_env%matrix_z
      IF(do_apply_preconditioner) THEN
        CALL apply_preconditioner(preconditioner,matrix_tmp,matrix_z,error=error)
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
   
      ! Project out converged vectors
!      IF(nmo_converged>0) THEN
      IF(.FALSE.) THEN
          block_mat => bdav_env%H_block_vec
          CALL cp_fm_set_all(block_mat,0.0_dp, 1.0_dp,error=error)
          CALL cp_fm_to_fm_submat(block_mat,bdav_env%S_block_mat,2*nmo,nmo_converged,1,1,1,1,error=error)
          CALL cp_fm_to_fm_submat(block_mat,bdav_env%S_block_mat,nmo_converged,2*nmo,1,1,1,1,error=error)
          block_mat => bdav_env%H_block_vec
          CALL cp_fm_set_all(block_mat,0.0_dp, 10000.0_dp,error=error)
          CALL cp_fm_to_fm_submat(block_mat,bdav_env%H_block_mat,2*nmo,nmo_converged,1,1,1,1,error=error)
          CALL cp_fm_to_fm_submat(block_mat,bdav_env%H_block_mat,nmo_converged,2*nmo,1,1,1,1,error=error)
      ELSE
        nmo_not_converged = nmo
        nmo_converged = 0
      END IF

    
      ! solution of the reduced eigenvalues problem
      CALL reduce_extended_space(bdav_env,error=error)
   
      block_mat => bdav_env%H_block_vec
      CALL cp_fm_to_fm_submat(block_mat,matrix_tmp,nmo,nmo_not_converged,1,1,1,1,error=error)
      CALL cp_fm_gemm('N','N',nao,nmo_not_converged,nmo,1.0_dp,mo_coeff,matrix_tmp,0.0_dp,matrix_hc,error=error)
      CALL cp_fm_to_fm_submat(block_mat,matrix_tmp,nmo,nmo_not_converged,1+nmo,1,1,1,error=error)
      CALL cp_fm_gemm('N','N',nao,nmo_not_converged,nmo,1.0_dp,matrix_z,matrix_tmp,1.0_dp,matrix_hc,error=error)
      CALL cp_fm_release(matrix_tmp,error=error)
   
      ! in case some vector are already converged only a subset of vectors are copied in the MOS
      IF (nmo_converged>0) THEN
        CALL cp_fm_to_fm_submat(matrix_hc,mo_coeff,nao,nmo_not_converged,1,1,1,nmo_converged+1,error=error)
        eigenvalues(nmo_converged+1:nmo) = bdav_env%evals(1:nmo_not_converged)
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
    DEALLOCATE(iconv,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ritz_coeff,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
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
!   CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,mo_coeff,matrix_sc,0.0_dp,matrix_tmp,error=error)
!    CALL cp_fm_get_diag(matrix_tmp,csc_diag,error=error)
    ! set the top left part of S[C,Z] block matrix  CSC
    block_mat => bdav_env%S_block_mat
    CALL cp_fm_set_all(block_mat,0.0_dp, 1.0_dp,error=error)
!    CALL cp_fm_to_fm_submat(matrix_tmp,block_mat,nmo,nmo,1,1,1,1,error=error)

    ! compute CHC
    CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,mo_coeff,matrix_hc,0.0_dp,matrix_tmp,error=error)
    CALL cp_fm_get_diag(matrix_tmp,chc_diag,error=error)
    ! set the top left part of H[C,Z] block matrix CHC 
    block_mat => bdav_env%H_block_mat
    CALL cp_fm_to_fm_submat(matrix_tmp,block_mat,nmo,nmo,1,1,1,1,error=error)
 
    DO i=1,nmo
!      IF(ABS(csc_diag(i))>EPSILON(0.0_dp)) THEN
        ritz_coeff(i) = chc_diag(i)!/csc_diag(i)
!      END IF
    END DO
    CALL cp_fm_release(matrix_tmp,error=error)

    CALL timestop(handle)

  END SUBROUTINE ritz_coefficients
! *****************************************************************************
  SUBROUTINE generate_extended_space_sparse(bdav_env,mo_set,matrix_h,matrix_s,output_unit,&
             preconditioner,error)

    TYPE(davidson_type)                      :: bdav_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h, matrix_s
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(preconditioner_type), OPTIONAL, &
      POINTER                                :: preconditioner
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'generate_extended_space_sparse', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, imo, istat, iter, max_iter, k,&
                                                n, nao,nmo, nmo_converged, nmo_not_converged
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iconv
    LOGICAL ::                                  converged, &
                                                do_apply_preconditioner, &
                                                failure
    REAL(dp)                                 :: max_norm, min_norm, t1, t2
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: ritz_coeff
    REAL(dp), DIMENSION(:), POINTER          :: eigenvalues

    TYPE(array_i1d_obj)                      :: col_blk_size, col_dist
    TYPE(array_i1d_obj)                      :: row_blk_size, row_dist

    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_b
    TYPE(cp_dbcsr_type), POINTER             :: matrix_hc, matrix_mm, matrix_sc, matrix_z, matrix_pz
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: block_mat, matrix_mm_fm, matrix_mm_fm_t, &
                                                matrix_nm_fm, mo_coeff, matrix_pz_fm
    TYPE(dbcsr_distribution_obj)             :: dist


    failure=.FALSE.
    CALL timeset(routineN,handle)

    do_apply_preconditioner = .FALSE.
    IF(PRESENT(preconditioner)) do_apply_preconditioner=.TRUE.

    NULLIFY(mo_coeff,mo_coeff_b, matrix_hc, matrix_sc, matrix_z, matrix_pz, matrix_mm)
    NULLIFY(block_mat, matrix_mm_fm, matrix_mm_fm_t, mo_coeff, matrix_nm_fm, matrix_pz_fm)
    CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, mo_coeff_b=mo_coeff_b,&
         eigenvalues=eigenvalues,nao=nao,nmo=nmo)
    IF(do_apply_preconditioner) THEN
         max_iter =  bdav_env%max_iter
    ELSE
       max_iter = 1
    END IF

    t1 = m_walltime()
    IF (output_unit > 0) THEN
          WRITE(output_unit,"(T15,A,T23,A,T36,A,T49,A,T60,A,/,T8,A)")   &
                " Cycle ", " conv. MOS ", " B2MAX ", " B2MIN ", " Time",  REPEAT("-",60)
    END IF


    ! Allocate array for Ritz values
    ALLOCATE(ritz_coeff(nmo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    Allocate(iconv(nmo), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! Allocate Sparse matrices: same distribution as mo_coeff_b => naoxnmo
    ! matrix_hc, matrix_sc, matrix_z, matrix_pz
    CALL cp_dbcsr_init_p(matrix_hc,error=error)
    CALL cp_dbcsr_create(matrix_hc,"matrix_hc",cp_dbcsr_distribution(mo_coeff_b),&
         dbcsr_type_no_symmetry,cp_dbcsr_row_block_sizes(mo_coeff_b),&
         cp_dbcsr_col_block_sizes(mo_coeff_b),0,0,dbcsr_type_real_default,error=error)
    CALL cp_dbcsr_init_p(matrix_sc,error=error)
    CALL cp_dbcsr_create(matrix_sc,"matrix_sc",cp_dbcsr_distribution(mo_coeff_b),&
         dbcsr_type_no_symmetry,cp_dbcsr_row_block_sizes(mo_coeff_b),&
         cp_dbcsr_col_block_sizes(mo_coeff_b),0,0,dbcsr_type_real_default,error=error)
    CALL cp_dbcsr_init_p(matrix_z,error=error)
    CALL cp_dbcsr_create(matrix_z,"matrix_z",cp_dbcsr_distribution(mo_coeff_b),&
         dbcsr_type_no_symmetry,cp_dbcsr_row_block_sizes(mo_coeff_b),&
         cp_dbcsr_col_block_sizes(mo_coeff_b),0,0,dbcsr_type_real_default,error=error)
    CALL cp_dbcsr_init_p(matrix_pz,error=error)
    CALL cp_dbcsr_create(matrix_pz,"matrix_pz",cp_dbcsr_distribution(mo_coeff_b),&
         dbcsr_type_no_symmetry,cp_dbcsr_row_block_sizes(mo_coeff_b),&
         cp_dbcsr_col_block_sizes(mo_coeff_b),0,0,dbcsr_type_real_default,error=error)

    ! Allocate Sparse matrices: nmoxnmo
    ! matrix_mm
    CALL cp_create_bl_distribution (col_dist, col_blk_size, nmo, &
          dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(mo_coeff_b))))
    CALL cp_create_bl_distribution (row_dist, row_blk_size, nmo, &
          dbcsr_mp_nprows(dbcsr_distribution_mp(cp_dbcsr_distribution(mo_coeff_b))))
    CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp(cp_dbcsr_distribution(mo_coeff_b)),&
          row_dist,col_dist)
    CALL cp_dbcsr_init_p(matrix_mm,error=error)
    CALL cp_dbcsr_create(matrix_mm,"matrix_mm",dist,dbcsr_type_no_symmetry,&
         row_blk_size,col_blk_size,0,0,dbcsr_type_real_default,error=error)
    CALL cp_dbcsr_distribution_release (dist)
    CALL array_release (col_blk_size)
    CALL array_release (col_dist)
    CALL array_release (row_blk_size)
    CALL array_release (row_dist)

    ! Allocate Full matrices: same structure as mo_coeff => naoxnmo
    ! matrix_pz_fm, matrix_nm_fm
    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nao,ncol_global=nmo, &
                             context=mo_coeff%matrix_struct%context, &
                             para_env=mo_coeff%matrix_struct%para_env,error=error)
    CALL cp_fm_create(matrix_nm_fm,fm_struct_tmp,name="matrix_nm",error=error)
    CALL cp_fm_create(matrix_pz_fm,fm_struct_tmp,name="matrix_pz",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    ! Allocate Full matrices: nmoxnmo
    ! matrix_mm_fm, matrix_mm_fm_t
    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nmo,ncol_global=nmo, &
                             context=mo_coeff%matrix_struct%context, &
                             para_env=mo_coeff%matrix_struct%para_env,error=error)
    CALL cp_fm_create(matrix_mm_fm,fm_struct_tmp,name="matrix_mm",error=error)
    CALL cp_fm_create(matrix_mm_fm_t,fm_struct_tmp,name="matrix_mm_t",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    converged=.FALSE.
    DO iter = 1, max_iter
      ! Prepare HC and SC, using mo_coeff_b (sparse), these are still sparse 

      CALL cp_dbcsr_get_info(mo_coeff_b,nfullrows_total=n,nfullcols_total=k)
      CALL cp_dbcsr_multiply('n','n',1.0_dp,matrix_h,mo_coeff_b,0.0_dp,matrix_hc,last_column=k,error=error)
      CALL cp_dbcsr_multiply('n','n',1.0_dp,matrix_s,mo_coeff_b,0.0_dp,matrix_sc,last_column=k,error=error)

      ! compute Ritz values 
      ritz_coeff=0.0_dp

      CALL cp_dbcsr_multiply('t','n',1.0_dp,mo_coeff_b,matrix_hc,0.0_dp,matrix_mm, last_column=k, error=error)
      CALL cp_dbcsr_get_diag(matrix_mm,ritz_coeff,error=error) 
      CALL mp_sum(ritz_coeff,mo_coeff%matrix_struct%para_env%group)

      ! set the top left part of H[C,Z] block matrix CHC 
      block_mat => bdav_env%H_block_mat
      !copy sparse to full
      CALL copy_dbcsr_to_fm(matrix_mm,matrix_mm_fm,error=error)
      CALL cp_fm_to_fm_submat(matrix_mm_fm,block_mat,nmo,nmo,1,1,1,1,error=error)

      ! set the top left part of S[C,Z] block matrix  CSC
      block_mat => bdav_env%S_block_mat
      CALL cp_fm_set_all(block_mat,0.0_dp, 1.0_dp,error=error)

      ! extended subspace P Z = P [H - theta S]C  this ia another matrix of type and size as mo_coeff_b 
      CALL cp_dbcsr_copy(matrix_z,matrix_sc,error=error)
      CALL cp_dbcsr_scale_by_vector(matrix_z,ritz_coeff,side='right',error=error)
      CALL cp_dbcsr_add(matrix_z,matrix_hc,-1.0_dp,1.0_dp,error=error)

      ! Check for converged eigenvectors
      ritz_coeff =0.0_dp
      CALL cp_dbcsr_norm(matrix_z,which_norm=dbcsr_norm_column,norm_vector=ritz_coeff,error=error)
      nmo_converged = 0
      max_norm = 0.0_dp
      min_norm = 1.e10_dp
      DO imo = 1,nmo
        max_norm = MAX(max_norm,ritz_coeff(imo))
        min_norm = MIN(min_norm,ritz_coeff(imo))
      END DO
      iconv = 0
      DO  imo = 1,nmo
        IF(ritz_coeff(imo) <= bdav_env%eps_iter ) THEN
            nmo_converged = nmo_converged + 1
            iconv(nmo_converged)=imo
!            econv(nmo_converged)=eigenvalues(imo)
!        ELSE
!          ! keep only sequential converged vectors, starting from lowest
!          exit
        END IF
      END DO
      nmo_not_converged=nmo-nmo_converged

!dbg
!      IF(nmo_converged>0) THEN
!        write(*,*) " NMO ", nmo, nmo_converged, nmo_not_converged
!        write(*,*) "iconv",  iconv(1:nmo_converged)
!      END IF
!dbg

      IF(REAL(nmo_converged,dp)/REAL(nmo,dp)>bdav_env%conv_percent) converged=.TRUE. 
   
   
      ! calculate or retrieve (if [H0- l S0]^-1) the preconditioner P: 
      !  this matrix has the full size nao x nao 
      IF(do_apply_preconditioner) THEN
        CALL apply_preconditioner(preconditioner,matrix_z,matrix_pz,error=error)
      ELSE
       CALL cp_dbcsr_copy(matrix_pz,matrix_z,error=error)
      END IF

!Dbg
 !     ! check preconditioned matrix
 !     ritz_coeff =0.0_dp
 !     CALL cp_dbcsr_norm(matrix_pz,which_norm=dbcsr_norm_column,norm_vector=ritz_coeff,error=error)
 !     write(*,*) 'norm pz ', ritz_coeff(:)
!    stop 'ritz'
!dbg

      ! compute the bottom left  ZSC (top right is transpose)
      CALL cp_dbcsr_multiply('t','n',1.0_dp,matrix_pz,matrix_sc,0.0_dp,matrix_mm,last_column=k,error=error)
      !  set the bottom left part of S[C,Z] block matrix  ZSC
      block_mat => bdav_env%S_block_mat
      !copy sparse to full
      CALL copy_dbcsr_to_fm(matrix_mm,matrix_mm_fm,error=error)
      CALL cp_fm_to_fm_submat(matrix_mm_fm,block_mat,nmo,nmo,1,1,1+nmo,1,error=error)
      CALL cp_fm_transpose(matrix_mm_fm,matrix_mm_fm_t,error=error)
      CALL cp_fm_to_fm_submat(matrix_mm_fm_t,block_mat,nmo,nmo,1,1,1,1+nmo,error=error)

      ! compute the bottom left  ZHC (top right is transpose)
      CALL cp_dbcsr_multiply('t','n',1.0_dp,matrix_pz,matrix_hc,0.0_dp,matrix_mm,last_column=k,error=error)
      ! set the bottom left part of S[C,Z] block matrix  ZHC
      block_mat => bdav_env%H_block_mat
      CALL copy_dbcsr_to_fm(matrix_mm,matrix_mm_fm,error=error)
      CALL cp_fm_to_fm_submat(matrix_mm_fm,block_mat,nmo,nmo,1,1,1+nmo,1,error=error)
      CALL cp_fm_transpose(matrix_mm_fm,matrix_mm_fm_t,error=error)
      CALL cp_fm_to_fm_submat(matrix_mm_fm_t,block_mat,nmo,nmo,1,1,1,1+nmo,error=error)

      ! (reuse matrix_sc and matrix_hc to computr HZ and SZ)
      CALL cp_dbcsr_get_info(matrix_pz,nfullrows_total=n,nfullcols_total=k)
      CALL cp_dbcsr_multiply('n','n',1.0_dp,matrix_h,matrix_pz,0.0_dp,matrix_hc,last_column=k,error=error)
      CALL cp_dbcsr_multiply('n','n',1.0_dp,matrix_s,matrix_pz,0.0_dp,matrix_sc,last_column=k,error=error)

      ! compute the bottom right  ZSZ 
      CALL cp_dbcsr_multiply('t','n',1.0_dp,matrix_pz,matrix_sc,0.0_dp,matrix_mm,last_column=k,error=error)
      ! set the bottom right part of S[C,Z] block matrix  ZSZ
      block_mat => bdav_env%S_block_mat
      CALL copy_dbcsr_to_fm(matrix_mm,matrix_mm_fm,error=error)
      CALL cp_fm_to_fm_submat(matrix_mm_fm,block_mat,nmo,nmo,1,1,1+nmo,1+nmo,error=error)

      ! compute the bottom right  ZHZ 
      CALL cp_dbcsr_multiply('t','n',1.0_dp,matrix_pz,matrix_hc,0.0_dp,matrix_mm,last_column=k,error=error)
      ! set the bottom right part of H[C,Z] block matrix  ZHZ
      block_mat => bdav_env%H_block_mat
      CALL copy_dbcsr_to_fm(matrix_mm,matrix_mm_fm,error=error)
      CALL cp_fm_to_fm_submat(matrix_mm_fm,block_mat,nmo,nmo,1,1,1+nmo,1+nmo,error=error)

!      IF(nmo_converged>0) THEN
      IF(.FALSE.) THEN

      ELSE
        nmo_not_converged = nmo
!        nmo_converged = 0
      END IF

      ! solution of the reduced eigenvalues problem
      CALL reduce_extended_space(bdav_env,error=error)
   
      block_mat => bdav_env%H_block_vec
      CALL cp_fm_to_fm_submat(block_mat,matrix_mm_fm,nmo,nmo_not_converged,1,1,1,1,error=error)
      CALL cp_fm_gemm('N','N',nao,nmo_not_converged,nmo,1.0_dp,mo_coeff,matrix_mm_fm,&
           0.0_dp,matrix_nm_fm,error=error)
      CALL cp_fm_to_fm_submat(block_mat,matrix_mm_fm,nmo,nmo_not_converged,1+nmo,1,1,1,error=error)
      CALL copy_dbcsr_to_fm(matrix_pz,matrix_pz_fm,error=error)
      CALL cp_fm_gemm('N','N',nao,nmo_not_converged,nmo,1.0_dp,matrix_pz_fm,matrix_mm_fm,&
           1.0_dp,matrix_nm_fm,error=error)
   
      ! in case some vector are already converged only a subset of vectors are copied in the MOS
!      IF (nmo_converged>0) THEN
      IF (.FALSE.) THEN
        CALL cp_fm_to_fm_submat(matrix_nm_fm,mo_coeff,nao,nmo_not_converged,1,1,1,nmo_converged+1,error=error)
        eigenvalues(nmo_converged+1:nmo) = bdav_env%evals(1:nmo_not_converged)
      ELSE
        CALL cp_fm_to_fm(matrix_nm_fm,mo_coeff,error=error)
        eigenvalues(1:nmo) = bdav_env%evals(1:nmo)
      END IF

      CALL copy_fm_to_dbcsr(mo_coeff,mo_coeff_b,error=error)!fm->dbcsr

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
    END DO ! iter

!dbg
 !    write(*,*) 'eval ', eigenvalues
!dbg

    DEALLOCATE(ritz_coeff,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(iconv,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL cp_dbcsr_release_p(matrix_hc, error=error)
    CALL cp_dbcsr_release_p(matrix_sc, error=error)
    CALL cp_dbcsr_release_p(matrix_z, error=error)
    CALL cp_dbcsr_release_p(matrix_pz, error=error)
    CALL cp_dbcsr_release_p(matrix_mm, error=error)

    CALL cp_fm_release(matrix_mm_fm,error=error)
    CALL cp_fm_release(matrix_mm_fm_t,error=error)
    CALL cp_fm_release(matrix_nm_fm,error=error)
    CALL cp_fm_release(matrix_pz_fm,error=error)

    CALL timestop(handle)

  END SUBROUTINE generate_extended_space_sparse

! *****************************************************************************
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
      CALL cp_fm_cholesky_decompose(S_block_mat,error=error)
      CALL cp_fm_triangular_invert(S_block_mat,error=error)
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
    CALL timestop(handle)


  END SUBROUTINE reduce_extended_space


END MODULE qs_scf_block_davidson
