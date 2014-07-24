!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief groups fairly general SCF methods, so that modules other than qs_scf can use them too
!>        split off from qs_scf to reduce dependencies
!> \par History
!>      - Joost VandeVondele (03.2006)
!>      - combine_ks_matrices added (05.04.06,MK)
!>      - second ROKS scheme added (15.04.06,MK)
!>      - MO occupation management moved (29.08.2008,MK)
! *****************************************************************************
MODULE qs_scf_methods

  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_desymmetrize, cp_dbcsr_get_block_p, &
       cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_p_type, &
       cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_symm,&
                                             cp_fm_triangular_multiply,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                             cp_fm_block_jacobi
  USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_set_element,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: cholesky_inverse,&
                                             cholesky_off,&
                                             cholesky_reduce,&
                                             cholesky_restore
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_max
  USE qs_density_mixing_types,         ONLY: mixing_storage_type
  USE qs_mo_types,                     ONLY: correct_mo_eigenvalues,&
                                             get_mo_set,&
                                             mo_set_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_methods'
  REAL(KIND=dp), PARAMETER    :: ratio = 0.25_dp

  PUBLIC :: combine_ks_matrices,&
            cp_sm_mix,&
            eigensolver,&
            eigensolver_dbcsr,&
            eigensolver_symm,&
            scf_env_density_mixing,&
            simple_eigensolver

  INTERFACE combine_ks_matrices
    MODULE PROCEDURE combine_ks_matrices_1,&
                     combine_ks_matrices_2
  END INTERFACE combine_ks_matrices

CONTAINS

! *****************************************************************************
!> \brief perform (if requested) a density mixing
!> \param p_mix_new ...
!> \param mixing_store ...
!> \param rho_ao ...
!> \param para_env ...
!> \param iter_delta ...
!> \param iter_count ...
!> \param diis ...
!> \param invert ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE scf_env_density_mixing(p_mix_new,mixing_store, rho_ao, para_env, iter_delta,&
                                    iter_count,diis,invert,error)
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p_mix_new
    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), INTENT(INOUT)             :: iter_delta
    INTEGER, INTENT(IN)                      :: iter_count
    LOGICAL, INTENT(in), OPTIONAL            :: diis, invert
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_density_mixing', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: failure, my_diis, my_invert
    REAL(KIND=dp)                            :: my_p_mix, tmp

    CALL timeset(routineN,handle)

    failure=.FALSE.
    my_diis=.FALSE.
    IF (PRESENT(diis)) my_diis=diis
    my_invert=.FALSE.
    IF (PRESENT(invert)) my_invert=invert
    my_p_mix = mixing_store%alpha
    IF (my_diis .OR. iter_count<mixing_store%nskip_mixing) THEN
        my_p_mix = 1.0_dp
    END IF

    iter_delta = 0.0_dp
    IF (.NOT.failure) THEN
       CPPrecondition(ASSOCIATED(p_mix_new),cp_failure_level,routineP,error,failure)
       DO ispin=1,SIZE(p_mix_new)
          IF (my_invert) THEN
             CPPrecondition(my_p_mix/=0.0_dp,cp_failure_level,routineP,error,failure)
             IF (my_p_mix/=1.0_dp) THEN
                CALL cp_dbcsr_add(matrix_a=p_mix_new(ispin)%matrix,&
                     alpha_scalar=1.0_dp/my_p_mix,&
                     matrix_b=rho_ao(ispin)%matrix,&
                     beta_scalar=(my_p_mix-1.0_dp)/my_p_mix,&
                     error=error)
             END IF
          ELSE
             CALL cp_sm_mix(m1=p_mix_new(ispin)%matrix,&
                  m2=rho_ao(ispin)%matrix,&
                  p_mix=my_p_mix,&
                  delta=tmp,&
                  para_env=para_env,&
                  error=error)
             iter_delta=MAX(iter_delta,tmp)
          END IF
       END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE scf_env_density_mixing

! *****************************************************************************
!> \brief   Diagonalise the Kohn-Sham matrix to get a new set of MO eigen-
!>          vectors and MO eigenvalues. ks will be modified
!> \param matrix_ks_fm ...
!> \param mo_set ...
!> \param ortho ...
!> \param work ...
!> \param cholesky_method ...
!> \param use_jacobi ...
!> \param error ...
!> \date    01.05.2001
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE eigensolver(matrix_ks_fm,mo_set,ortho,work,&
                         cholesky_method,use_jacobi,&
                         error)
    TYPE(cp_fm_type), POINTER                :: matrix_ks_fm
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: ortho, work
    INTEGER, INTENT(INout)                   :: cholesky_method
    LOGICAL, INTENT(IN)                      :: use_jacobi
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eigensolver', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homo, nao, nmo
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    CALL timeset(routineN,handle)

    NULLIFY (mo_coeff)
    NULLIFY (mo_eigenvalues)


    ! Diagonalise the Kohn-Sham matrix

    CALL get_mo_set(mo_set=mo_set,&
                    nao=nao,&
                    nmo=nmo,&
                    homo=homo, &
                    eigenvalues=mo_eigenvalues,&
                    mo_coeff=mo_coeff)

       SELECT CASE(cholesky_method)
       CASE(cholesky_reduce)
         CALL cp_fm_cholesky_reduce(matrix_ks_fm,ortho,error=error)
         CALL choose_eigv_solver(matrix_ks_fm,work,mo_eigenvalues,error=error)
         CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"SOLVE",error=error)
       CASE(cholesky_restore)
         CALL cp_fm_upper_to_full(matrix_ks_fm,work,error=error)
         CALL cp_fm_cholesky_restore(matrix_ks_fm,nao,ortho,work,&
                "SOLVE",pos="RIGHT",error=error)
         CALL cp_fm_cholesky_restore(work,nao,ortho,matrix_ks_fm,&
               "SOLVE",pos="LEFT",transa="T",error=error)
         CALL choose_eigv_solver(matrix_ks_fm,work,mo_eigenvalues,error=error)
         CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"SOLVE",error=error)
       CASE(cholesky_inverse)
         CALL cp_fm_upper_to_full(matrix_ks_fm,work,error=error)

         CALL cp_fm_triangular_multiply(ortho,matrix_ks_fm ,side="R",transpose_tr=.FALSE.,&
               invert_tr=.FALSE., uplo_tr="U",n_rows=nao,n_cols=nao,alpha=1.0_dp,error=error)
         CALL cp_fm_triangular_multiply(ortho, matrix_ks_fm,side="L",transpose_tr=.TRUE.,&
               invert_tr=.FALSE., uplo_tr="U",n_rows=nao,n_cols=nao,alpha=1.0_dp,error=error)
         CALL choose_eigv_solver(matrix_ks_fm,work,mo_eigenvalues,error=error)
         CALL cp_fm_triangular_multiply(ortho,work,side="L",transpose_tr=.FALSE.,&
              invert_tr=.FALSE., uplo_tr="U",n_rows=nao,n_cols=nmo,alpha=1.0_dp,error=error)
         CALL cp_fm_to_fm(work,mo_coeff,nmo,1,1)

       END SELECT

       IF (use_jacobi) THEN
          CALL cp_fm_to_fm(mo_coeff,ortho,error=error)
          cholesky_method = cholesky_off
       END IF

    CALL timestop(handle)

  END SUBROUTINE eigensolver

! *****************************************************************************
!> \brief ...
!> \param matrix_ks ...
!> \param matrix_ks_fm ...
!> \param mo_set ...
!> \param ortho_dbcsr ...
!> \param ksbuf1 ...
!> \param ksbuf2 ...
!> \param work ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE eigensolver_dbcsr(matrix_ks,matrix_ks_fm,mo_set,ortho_dbcsr,ksbuf1,ksbuf2,work,&
                         error)
    TYPE(cp_dbcsr_type), POINTER             :: matrix_ks
    TYPE(cp_fm_type), POINTER                :: matrix_ks_fm
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: ortho_dbcsr, ksbuf1, ksbuf2
    TYPE(cp_fm_type), POINTER                :: work
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eigensolver_dbcsr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao, nmo
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    CALL timeset(routineN,handle)

    NULLIFY (mo_coeff)
    NULLIFY (mo_eigenvalues)

    CALL get_mo_set(mo_set=mo_set,&
                    nao=nao,&
                    nmo=nmo,&
                    eigenvalues=mo_eigenvalues,&
                    mo_coeff=mo_coeff)

!    Reduce KS matrix
      CALL cp_dbcsr_desymmetrize(matrix_ks,ksbuf2,error=error)
      CALL cp_dbcsr_multiply('N','N',1.0_dp,ksbuf2,ortho_dbcsr,0.0_dp,ksbuf1,error=error)
      CALL cp_dbcsr_multiply('T','N',1.0_dp,ortho_dbcsr,ksbuf1,0.0_dp,ksbuf2,error=error)

!    Solve the eigenvalue problem
      CALL copy_dbcsr_to_fm(ksbuf2, matrix_ks_fm, error=error)
      CALL choose_eigv_solver(matrix_ks_fm,work,mo_eigenvalues,error=error)

!    Restore the eigenvector of the general eig. problem
      CALL cp_dbcsr_sm_fm_multiply(ortho_dbcsr,work,mo_coeff,nmo, error=error)



    CALL timestop(handle)

  END SUBROUTINE eigensolver_dbcsr


! *****************************************************************************
!> \brief ...
!> \param matrix_ks_fm ...
!> \param mo_set ...
!> \param ortho ...
!> \param work ...
!> \param do_level_shift ...
!> \param level_shift ...
!> \param use_jacobi ...
!> \param jacobi_threshold ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE eigensolver_symm(matrix_ks_fm,mo_set,ortho,work,do_level_shift,&
                         level_shift,use_jacobi,&
                         jacobi_threshold,error)
    TYPE(cp_fm_type), POINTER                :: matrix_ks_fm
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: ortho, work
    LOGICAL, INTENT(IN)                      :: do_level_shift
    REAL(KIND=dp), INTENT(IN)                :: level_shift
    LOGICAL, INTENT(IN)                      :: use_jacobi
    REAL(KIND=dp), INTENT(IN)                :: jacobi_threshold
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eigensolver_symm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homo, imo, nao, &
                                                nelectron, nmo
    REAL(KIND=dp)                            :: alpha
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    CALL timeset(routineN,handle)

    NULLIFY (mo_coeff)
    NULLIFY (mo_eigenvalues)


    ! Diagonalise the Kohn-Sham matrix

    CALL get_mo_set(mo_set=mo_set,&
                    nao=nao,&
                    nmo=nmo,&
                    homo=homo, &
                    nelectron=nelectron,&
                    eigenvalues=mo_eigenvalues,&
                    mo_coeff=mo_coeff)

    IF (use_jacobi) THEN

       CALL cp_fm_symm("L","U",nao,homo,1.0_dp,matrix_ks_fm,mo_coeff,0.0_dp,work,error=error)
       CALL cp_gemm("T","N",homo,nao-homo,nao,1.0_dp,work,mo_coeff,&
                       0.0_dp,matrix_ks_fm,b_first_col=homo+1,c_first_col=homo+1,error=error)

       ! Block Jacobi (pseudo-diagonalization, only one sweep)
       CALL cp_fm_block_jacobi(matrix_ks_fm,mo_coeff,mo_eigenvalues,&
                               jacobi_threshold,homo+1)

    ELSE ! full S^(-1/2) has been computed

       CALL cp_fm_symm("L","U",nao,nao,1.0_dp,matrix_ks_fm,ortho,0.0_dp,work,error=error)
       CALL cp_gemm("T","N",nao,nao,nao,1.0_dp,ortho,work,0.0_dp,matrix_ks_fm,error=error)

       IF (do_level_shift) THEN
          DO imo=homo+1,nmo
             ! now unnecessary broadcast of element
             CALL cp_fm_get_element(matrix_ks_fm,imo,imo,alpha)
             alpha=alpha+level_shift
             CALL cp_fm_set_element(matrix_ks_fm,imo,imo,alpha,error=error)
          END DO
       END IF

       CALL choose_eigv_solver(matrix_ks_fm,work,mo_eigenvalues,error=error)

       CALL cp_gemm("N","N",nao,nmo,nao,1.0_dp,ortho,work,0.0_dp,&
                       mo_coeff,error=error)

       IF (do_level_shift) THEN
         CALL correct_mo_eigenvalues(mo_set,level_shift)
       END IF

    END IF


    CALL timestop(handle)

  END SUBROUTINE eigensolver_symm


! *****************************************************************************


! *****************************************************************************
!> \brief ...
!> \param matrix_ks ...
!> \param mo_set ...
!> \param work ...
!> \param do_level_shift ...
!> \param level_shift ...
!> \param use_jacobi ...
!> \param jacobi_threshold ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE simple_eigensolver(matrix_ks,mo_set,work,do_level_shift,&
                                level_shift,use_jacobi,jacobi_threshold,error)

    TYPE(cp_fm_type), POINTER                :: matrix_ks
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: work
    LOGICAL, INTENT(IN)                      :: do_level_shift
    REAL(KIND=dp), INTENT(IN)                :: level_shift
    LOGICAL, INTENT(IN)                      :: use_jacobi
    REAL(KIND=dp), INTENT(IN)                :: jacobi_threshold
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'simple_eigensolver', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homo, imo, nao, &
                                                nelectron, nmo
    REAL(KIND=dp)                            :: alpha
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    CALL timeset(routineN,handle)

    NULLIFY (mo_coeff)
    NULLIFY (mo_eigenvalues)


    CALL get_mo_set(mo_set=mo_set,&
                    nao=nao,&
                    nmo=nmo,&
                    homo=homo, &
                    nelectron=nelectron,&
                    eigenvalues=mo_eigenvalues,&
                    mo_coeff=mo_coeff)

    IF (do_level_shift) THEN
      DO imo=homo+1,nmo
         ! now unnecessary broadcast of element
         CALL cp_fm_get_element(matrix_ks,imo,imo,alpha)
         alpha=alpha+level_shift
         CALL cp_fm_set_element(matrix_ks,imo,imo,alpha,error=error)
      END DO
    END IF

    IF ( use_jacobi ) THEN
       CALL cp_fm_symm("L","U",nao,homo,1.0_dp,matrix_ks,mo_coeff,0.0_dp,work,error=error)
       CALL cp_gemm("T","N",homo,nao-homo,nao,1.0_dp,work,mo_coeff,&
                       0.0_dp,matrix_ks,b_first_col=homo+1,c_first_col=homo+1,error=error)
       ! Block Jacobi (pseudo-diagonalization, only one sweep)
       CALL cp_fm_block_jacobi(matrix_ks,mo_coeff,mo_eigenvalues,jacobi_threshold,homo+1)
    ELSE

       CALL choose_eigv_solver(matrix_ks,work,mo_eigenvalues,error=error)

      CALL cp_fm_to_fm(work,mo_coeff,nmo,1,1)

    END IF

    IF (do_level_shift) THEN
      CALL correct_mo_eigenvalues(mo_set,level_shift)
    END IF

    CALL timestop(handle)

  END SUBROUTINE simple_eigensolver

! *****************************************************************************
!> \brief Perform a mixing of the given matrixes into the first matrix
!>      m1 = m2 + p_mix (m1-m2)
!> \param m1 ...
!> \param m2 ...
!> \param p_mix how much m1 is conserved (0: none, 1: all)
!> \param delta maximum norm of m1-m2
!> \param para_env ...
!> \param m3 ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param m 1: first (new) matrix, is modified UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param m 2: the second (old) matrix UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par History
!>      02.2003 rewamped [fawzi]
!> \author fawzi
!> \note
!>      if you what to store the result in m2 swap m1 and m2 an use
!>      (1-pmix) as pmix
!>      para_env should be removed (embedded in matrix)
! *****************************************************************************
  SUBROUTINE cp_sm_mix(m1,m2,p_mix,delta,para_env,m3,error)

    TYPE(cp_dbcsr_type), POINTER             :: m1, m2
    REAL(KIND=dp), INTENT(IN)                :: p_mix
    REAL(KIND=dp), INTENT(OUT)               :: delta
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: m3
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_mix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, i, iblock_col, &
                                                iblock_row, j
    LOGICAL                                  :: failure, found
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_delta_block, p_new_block, &
                                                p_old_block
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL timeset(routineN,handle)
    failure = .FALSE.
    delta = 0.0_dp

    CALL cp_dbcsr_iterator_start(iter, m1)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, iblock_row, iblock_col, p_new_block, blk)
       CALL cp_dbcsr_get_block_p(matrix=m2,row=iblock_row,col=iblock_col,&
            BLOCK=p_old_block,found=found)
       CPPostcondition(ASSOCIATED(p_old_block),cp_failure_level,routineP,error,failure)
       IF(PRESENT(m3)) THEN
          CALL cp_dbcsr_get_block_p(matrix=m3,row=iblock_row,col=iblock_col,&
               BLOCK=p_delta_block,found=found)
          CPPostcondition(ASSOCIATED(p_delta_block),cp_failure_level,routineP,error,failure)

          DO j=1,SIZE(p_new_block,2)
             DO i=1,SIZE(p_new_block,1)
                p_delta_block(i,j) = p_new_block(i,j) - p_old_block(i,j)
                delta = MAX(delta,ABS(p_delta_block(i,j)))
             END DO
          END DO
       ELSE
          DO j=1,SIZE(p_new_block,2)
             DO i=1,SIZE(p_new_block,1)
                p_new_block(i,j) = p_new_block(i,j) - p_old_block(i,j)
                delta = MAX(delta,ABS(p_new_block(i,j)))
                p_new_block(i,j) = p_old_block(i,j) + p_mix*p_new_block(i,j)
             END DO
          END DO
       ENDIF
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL mp_max(delta,para_env%group)

    CALL timestop(handle)

  END SUBROUTINE cp_sm_mix

! *****************************************************************************
!> \brief ...
!> \param ksa ...
!> \param ksb ...
!> \param occa ...
!> \param occb ...
!> \param roks_parameter ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE combine_ks_matrices_1(ksa,ksb,occa,occb,roks_parameter,error)

    ! Combine the alpha and beta Kohn-Sham matrices during a restricted open
    ! Kohn-Sham (ROKS) calculation
    ! On input ksa and ksb contain the alpha and beta Kohn-Sham matrices,
    ! respectively. occa and occb contain the corresponding MO occupation
    ! numbers. On output the combined ROKS operator matrix is returned in ksa.

    ! Literature: - C. C. J. Roothaan, Rev. Mod. Phys. 32, 179 (1960)
    !             - M. F. Guest and V. R. Saunders, Mol. Phys. 28(3), 819 (1974)

    TYPE(cp_fm_type), POINTER                :: ksa, ksb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: occa, occb
    REAL(KIND=dp), &
      DIMENSION(0:2, 0:2, 1:2), INTENT(IN)   :: roks_parameter
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'combine_ks_matrices_1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, icol_global, &
                                                icol_local, irow_global, &
                                                irow_local, j, ncol_local, &
                                                nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: compatible_matrices, failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fa, fb
    TYPE(cp_fm_struct_type), POINTER         :: ksa_struct, ksb_struct

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(ksa),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(ksb),cp_failure_level,routineP,error,failure)

    CALL cp_fm_get_info(matrix=ksa,&
                        matrix_struct=ksa_struct,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        local_data=fa,&
                        error=error)

    CALL cp_fm_get_info(matrix=ksb,&
                        matrix_struct=ksb_struct,&
                        local_data=fb,&
                        error=error)

    compatible_matrices = cp_fm_struct_equivalent(ksa_struct,ksb_struct,error=error)
    CPPrecondition(compatible_matrices,cp_failure_level,routineP,error,failure)

    IF (SUM(occb) == 0.0_dp) fb = 0.0_dp

    DO icol_local=1,ncol_local
      icol_global = col_indices(icol_local)
      j = INT(occa(icol_global)) + INT(occb(icol_global))
      DO irow_local=1,nrow_local
        irow_global = row_indices(irow_local)
        i = INT(occa(irow_global)) + INT(occb(irow_global))
        fa(irow_local,icol_local) =&
          roks_parameter(i,j,1)*fa(irow_local,icol_local) +&
          roks_parameter(i,j,2)*fb(irow_local,icol_local)
      END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE combine_ks_matrices_1

! *****************************************************************************
!> \brief ...
!> \param ksa ...
!> \param ksb ...
!> \param occa ...
!> \param occb ...
!> \param f ...
!> \param nalpha ...
!> \param nbeta ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE combine_ks_matrices_2(ksa,ksb,occa,occb,f,nalpha,nbeta,error)

    ! Combine the alpha and beta Kohn-Sham matrices during a restricted open
    ! Kohn-Sham (ROKS) calculation
    ! On input ksa and ksb contain the alpha and beta Kohn-Sham matrices,
    ! respectively. occa and occb contain the corresponding MO occupation
    ! numbers. On output the combined ROKS operator matrix is returned in ksa.

    ! Literature: - C. C. J. Roothaan, Rev. Mod. Phys. 32, 179 (1960)
    !             - M. Filatov and S. Shaik, Chem. Phys. Lett. 288, 689 (1998)

    TYPE(cp_fm_type), POINTER                :: ksa, ksb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: occa, occb
    REAL(KIND=dp), INTENT(IN)                :: f
    INTEGER, INTENT(IN)                      :: nalpha, nbeta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'combine_ks_matrices_2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, icol_global, &
                                                icol_local, irow_global, &
                                                irow_local, ncol_local, &
                                                nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: compatible_matrices, failure
    REAL(KIND=dp)                            :: beta, t1, t2, ta, tb
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fa, fb
    TYPE(cp_fm_struct_type), POINTER         :: ksa_struct, ksb_struct

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(ksa),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(ksb),cp_failure_level,routineP,error,failure)

    CALL cp_fm_get_info(matrix=ksa,&
                        matrix_struct=ksa_struct,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        local_data=fa,&
                        error=error)

    CALL cp_fm_get_info(matrix=ksb,&
                        matrix_struct=ksb_struct,&
                        local_data=fb,&
                        error=error)

    compatible_matrices = cp_fm_struct_equivalent(ksa_struct,ksb_struct,error=error)
    CPPrecondition(compatible_matrices,cp_failure_level,routineP,error,failure)

    beta = 1.0_dp/(1.0_dp - f)

    DO icol_local=1,ncol_local

      icol_global = col_indices(icol_local)

      DO irow_local=1,nrow_local

        irow_global = row_indices(irow_local)

        t1 = 0.5_dp*(fa(irow_local,icol_local) + fb(irow_local,icol_local))

        IF ((0 < irow_global).AND.(irow_global <= nbeta)) THEN
          IF ((0 < icol_global).AND.(icol_global <= nbeta)) THEN
            ! closed-closed
            fa(irow_local,icol_local) = t1
          ELSE IF ((nbeta < icol_global).AND.(icol_global <= nalpha)) THEN
            ! closed-open
            ta = 0.5_dp*(f - REAL(occa(icol_global),KIND=dp))/f
            tb = 0.5_dp*(f - REAL(occb(icol_global),KIND=dp))/f
            t2 = ta*fa(irow_local,icol_local) + tb*fb(irow_local,icol_local)
            fa(irow_local,icol_local) = t1 + (beta - 1.0_dp)*t2
          ELSE
            ! closed-virtual
            fa(irow_local,icol_local) = t1
          END IF
        ELSE IF ((nbeta < irow_global).AND.(irow_global <= nalpha)) THEN
          IF ((0 < irow_global).AND.(irow_global <= nbeta)) THEN
            ! open-closed
            ta = 0.5_dp*(f - REAL(occa(irow_global),KIND=dp))/f
            tb = 0.5_dp*(f - REAL(occb(irow_global),KIND=dp))/f
            t2 = ta*fa(irow_local,icol_local) + tb*fb(irow_local,icol_local)
            fa(irow_local,icol_local) = t1 + (beta - 1.0_dp)*t2
          ELSE IF ((nbeta < icol_global).AND.(icol_global <= nalpha)) THEN
            ! open-open
            ta = 0.5_dp*(f - REAL(occa(icol_global),KIND=dp))/f
            tb = 0.5_dp*(f - REAL(occb(icol_global),KIND=dp))/f
            t2 = ta*fa(irow_local,icol_local) + tb*fb(irow_local,icol_local)
            IF (irow_global == icol_global) THEN
              fa(irow_local,icol_local) = t1 - t2
            ELSE
              fa(irow_local,icol_local) = t1 - 0.5_dp*t2
            END IF
          ELSE
            ! open-virtual
            ta = 0.5_dp*(f - REAL(occa(irow_global),KIND=dp))/f
            tb = 0.5_dp*(f - REAL(occb(irow_global),KIND=dp))/f
            t2 = ta*fa(irow_local,icol_local) + tb*fb(irow_local,icol_local)
            fa(irow_local,icol_local) = t1 - t2
          END IF
        ELSE
          IF ((0 < irow_global).AND.(irow_global < nbeta)) THEN
            ! virtual-closed
            fa(irow_local,icol_local) = t1
          ELSE IF ((nbeta < icol_global).AND.(icol_global <= nalpha)) THEN
            ! virtual-open
            ta = 0.5_dp*(f - REAL(occa(icol_global),KIND=dp))/f
            tb = 0.5_dp*(f - REAL(occb(icol_global),KIND=dp))/f
            t2 = ta*fa(irow_local,icol_local) + tb*fb(irow_local,icol_local)
            fa(irow_local,icol_local) = t1 - t2
          ELSE
            ! virtual-virtual
            fa(irow_local,icol_local) = t1
          END IF
        END IF

      END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE combine_ks_matrices_2

END MODULE qs_scf_methods
