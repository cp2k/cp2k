!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
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

  USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                             dbcsr_deallocate_matrix_set,&
                                             dbcsr_from_sm,&
                                             dbcsr_get_block_node,&
                                             dbcsr_obj_p_type,&
                                             sm_from_dbcsr
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_symm,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: cp_fm_block_jacobi,&
                                             cp_fm_syevd
  USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_set_element,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_methods,                   ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_operations,                ONLY: dbcsr_add
  USE dbcsr_types,                     ONLY: dbcsr_iterator,&
                                             dbcsr_obj
  USE f77_blas
  USE input_constants,                 ONLY: cholesky_inverse,&
                                             cholesky_off,&
                                             cholesky_reduce,&
                                             cholesky_restore
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_max
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_mo_types,                     ONLY: correct_mo_eigenvalues,&
                                             get_mo_set,&
                                             mo_set_type
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_methods'
  REAL(KIND=dp), PARAMETER    :: ratio = 0.25_dp

  PUBLIC :: combine_ks_matrices,&
            cp_sm_mix,&
            eigensolver,&
            scf_env_density_mixing,&
            simple_eigensolver

  INTERFACE combine_ks_matrices
    MODULE PROCEDURE combine_ks_matrices_1,&
                     combine_ks_matrices_2
  END INTERFACE combine_ks_matrices

CONTAINS

! *****************************************************************************
!> \brief perform (if requested) a density mixing
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE scf_env_density_mixing(p_mix_new_sm,scf_env,iter_delta,&
                                    qs_env,diis,invert,error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: p_mix_new_sm
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    REAL(KIND=dp), INTENT(INOUT)             :: iter_delta
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in)                      :: diis
    LOGICAL, INTENT(in), OPTIONAL            :: invert
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_density_mixing', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: failure, my_invert
    REAL(KIND=dp)                            :: my_p_mix, tmp
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: p_mix_new, rho_ao

    CALL timeset(routineN,handle)

    NULLIFY(p_mix_new,rho_ao)
    CALL dbcsr_allocate_matrix_set( p_mix_new, SIZE(p_mix_new_sm), error )
    CALL dbcsr_allocate_matrix_set( rho_ao, SIZE(p_mix_new_sm), error )
    DO ispin=1,SIZE(p_mix_new)
       ALLOCATE(p_mix_new(ispin)%matrix,rho_ao(ispin)%matrix)!super dirty
       CALL dbcsr_from_sm(p_mix_new(ispin)%matrix, p_mix_new_sm(ispin)%matrix, error)
       CALL dbcsr_from_sm(rho_ao(ispin)%matrix, qs_env%rho%rho_ao(ispin)%matrix, error)
    ENDDO

    failure=.FALSE.
    my_invert=.FALSE.
    IF (PRESENT(invert)) my_invert=invert
    my_p_mix = scf_env%p_mix_alpha
    IF (diis .OR. scf_env%iter_count<scf_env%nskip_mixing) THEN
        my_p_mix = 1.0_dp
    END IF

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(ASSOCIATED(p_mix_new),cp_failure_level,routineP,error,failure)
       DO ispin=1,SIZE(p_mix_new)
          IF (my_invert) THEN
             CPPrecondition(my_p_mix/=0.0_dp,cp_failure_level,routineP,error,failure)
             IF (my_p_mix/=1.0_dp) THEN
                CALL dbcsr_add(matrix_a=p_mix_new(ispin)%matrix,&
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
                  para_env=qs_env%para_env,&
                  error=error)
             iter_delta=MAX(iter_delta,tmp)
          END IF
       END DO
    END IF

    DO ispin=1,SIZE(p_mix_new)
       CALL sm_from_dbcsr( p_mix_new_sm(ispin)%matrix, p_mix_new(ispin)%matrix, &
            qs_env%rho%rho_ao(ispin)%matrix%distribution_2d, error=error )
    ENDDO

    CALL dbcsr_deallocate_matrix_set( p_mix_new, error=error )
    CALL dbcsr_deallocate_matrix_set( rho_ao, error=error )

    CALL timestop(handle)

  END SUBROUTINE scf_env_density_mixing

! *****************************************************************************
!> \brief   Diagonalise the Kohn-Sham matrix to get a new set of MO eigen-
!>          vectors and MO eigenvalues. ks will be modified
!> \author  Matthias Krack
!> \date    01.05.2001
!> \version 1.0
! *****************************************************************************
  SUBROUTINE eigensolver(matrix_ks,mo_set,ortho,work,do_level_shift,&
                         level_shift,cholesky_method,use_jacobi,&
                         jacobi_threshold,error)
    TYPE(cp_fm_type), POINTER                :: matrix_ks
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: ortho, work
    LOGICAL, INTENT(IN)                      :: do_level_shift
    REAL(KIND=dp), INTENT(IN)                :: level_shift
    INTEGER, INTENT(INout)                   :: cholesky_method
    LOGICAL, INTENT(IN)                      :: use_jacobi
    REAL(KIND=dp), INTENT(IN)                :: jacobi_threshold
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eigensolver', &
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

    IF (do_level_shift) THEN
       IF (cholesky_method>cholesky_off) THEN
          CALL stop_program("eigensolver","level shift not implemented")
       ENDIF
    ENDIF

    IF (cholesky_method>cholesky_off) THEN
       SELECT CASE(cholesky_method)
       CASE(cholesky_reduce)
         CALL cp_fm_cholesky_reduce(matrix_ks,ortho,error=error)
       CASE(cholesky_restore)
         CALL cp_fm_upper_to_full(matrix_ks,work,error=error)
         CALL cp_fm_cholesky_restore(matrix_ks,nao,ortho,work,&
                "SOLVE",pos="RIGHT",error=error)
         CALL cp_fm_cholesky_restore(work,nao,ortho,matrix_ks,&
               "SOLVE",pos="LEFT",transa="T",error=error)
       CASE(cholesky_inverse)
         CALL cp_fm_upper_to_full(matrix_ks,work,error=error)
         CALL cp_fm_cholesky_restore(matrix_ks,nao,ortho,work,&
                "MULTIPLY",pos="RIGHT",error=error)
         CALL cp_fm_cholesky_restore(work,nao,ortho,matrix_ks,&
               "MULTIPLY",pos="LEFT",transa="T",error=error)
       END SELECT

       IF (use_jacobi) THEN
          CALL cp_fm_syevd(matrix_ks,work,mo_eigenvalues,error=error)
          IF(cholesky_method==cholesky_inverse) THEN
            CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"MULTIPLY",error=error)
          ELSE
            CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"SOLVE",error=error)
          END IF
          CALL cp_fm_to_fm(mo_coeff,ortho,error=error)
          cholesky_method = cholesky_off
       ELSEIF(nmo==nao) THEN
          ! In most of the cases syevd is faster than syevx, also when anly partial spectrum is needed
          CALL cp_fm_syevd(matrix_ks,work,mo_eigenvalues,error=error)
          IF(cholesky_method==cholesky_inverse) THEN
            CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"MULTIPLY",error=error)
          ELSE
            CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"SOLVE",error=error)
          END IF
       ELSE
          CALL cp_fm_syevd(matrix_ks,work,mo_eigenvalues,error=error)
!          CALL cp_fm_syevr(matrix_ks,work,mo_eigenvalues,1,nmo,error=error)
          IF(cholesky_method==cholesky_inverse) THEN
            CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"MULTIPLY",error=error)
          ELSE
            CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"SOLVE",error=error)
          END IF
       END IF

    ELSE IF (use_jacobi) THEN

       CALL cp_fm_symm("L","U",nao,homo,1.0_dp,matrix_ks,mo_coeff,0.0_dp,work,error=error)
       CALL cp_fm_gemm("T","N",homo,nao-homo,nao,1.0_dp,work,mo_coeff,&
                       0.0_dp,matrix_ks,b_first_col=homo+1,c_first_col=homo+1,error=error)

       ! Block Jacobi (pseudo-diagonalization, only one sweep)

       CALL cp_fm_block_jacobi(matrix_ks,mo_coeff,mo_eigenvalues,&
                               jacobi_threshold,homo+1)

    ELSE ! full S^(-1/2) has been computed

       CALL cp_fm_symm("L","U",nao,nao,1.0_dp,matrix_ks,ortho,0.0_dp,work,error=error)
       CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_dp,ortho,work,0.0_dp,matrix_ks,error=error)

       IF (do_level_shift) THEN
          DO imo=homo+1,nmo
             ! now unnecessary broadcast of element
             CALL cp_fm_get_element(matrix_ks,imo,imo,alpha)
             alpha=alpha+level_shift
             CALL cp_fm_set_element(matrix_ks,imo,imo,alpha,error=error)
          END DO
       END IF

       ! In most of the cases syevd is faster than syevx, also when anly partial spectrum is needed
       IF(nmo==nao) THEN
         CALL cp_fm_syevd(matrix_ks,work,mo_eigenvalues,error)
       ELSE
         CALL cp_fm_syevd(matrix_ks,work,mo_eigenvalues,error)
!         CALL cp_fm_syevr(matrix_ks,work,mo_eigenvalues,1,nmo,error=error)
       END IF

       CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,ortho,work,0.0_dp,&
                       mo_coeff,error=error)

       IF (do_level_shift) THEN
         CALL correct_mo_eigenvalues(mo_set,level_shift)
       END IF

    END IF

    CALL timestop(handle)

  END SUBROUTINE eigensolver

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

    INTEGER                                  :: handle, homo, imo, istat, &
                                                nao, nelectron, nmo
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp)                            :: alpha
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigvals
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
       CALL cp_fm_gemm("T","N",homo,nao-homo,nao,1.0_dp,work,mo_coeff,&
                       0.0_dp,matrix_ks,b_first_col=homo+1,c_first_col=homo+1,error=error)
       ! Block Jacobi (pseudo-diagonalization, only one sweep)
       CALL cp_fm_block_jacobi(matrix_ks,mo_coeff,mo_eigenvalues,jacobi_threshold,homo+1)
    ELSE
      ALLOCATE(eigvals(nao),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ! In most of the cases syevd is faster than syevx, also when anly partial spectrum is needed
      CALL cp_fm_syevd(matrix_ks,work,eigvals,error)

      mo_eigenvalues(1:nmo) = eigvals(1:nmo)
      CALL cp_fm_to_fm(work,mo_coeff,nmo,1,1)

      DEALLOCATE(eigvals,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    IF (do_level_shift) THEN
      CALL correct_mo_eigenvalues(mo_set,level_shift)
    END IF

    CALL timestop(handle)

  END SUBROUTINE simple_eigensolver

! *****************************************************************************
!> \brief Perform a mixing of the given matrixes into the first matrix
!>      m1 = m2 + p_mix (m1-m2)
!> \param m 1: first (new) matrix, is modified
!> \param m 2: the second (old) matrix
!> \param p_mix how much m1 is conserved (0: none, 1: all)
!> \param delta maximum norm of m1-m2
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      if you what to store the result in m2 swap m1 and m2 an use
!>      (1-pmix) as pmix
!>      para_env should be removed (embedded in matrix)
!> \par History
!>      02.2003 rewamped [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE cp_sm_mix(m1,m2,p_mix,delta,para_env,m3,error)

    TYPE(dbcsr_obj), POINTER                 :: m1, m2
    REAL(KIND=dp), INTENT(IN)                :: p_mix
    REAL(KIND=dp), INTENT(OUT)               :: delta
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_obj), OPTIONAL, POINTER       :: m3
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_mix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, i, iblock_col, &
                                                iblock_row, j
    LOGICAL                                  :: failure, tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_delta_block, p_new_block, &
                                                p_old_block
    TYPE(dbcsr_iterator)                     :: iter

    CALL timeset(routineN,handle)
    failure = .FALSE.
    delta = 0.0_dp

    CALL dbcsr_iterator_start(iter, m1)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, p_new_block, tr, blk)
       CALL dbcsr_get_block_node(matrix=m2,&
                                 block_row=iblock_row,&
                                 block_col=iblock_col,&
                                 BLOCK=p_old_block)
       CPPostcondition(ASSOCIATED(p_old_block),cp_failure_level,routineP,error,failure)
       IF(PRESENT(m3)) THEN
          CALL dbcsr_get_block_node(matrix=m3,&
                                    block_row=iblock_row,&
                                    block_col=iblock_col,&
                                    BLOCK=p_delta_block)
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
    CALL dbcsr_iterator_stop(iter)

    CALL mp_max(delta,para_env%group)

    CALL timestop(handle)

  END SUBROUTINE cp_sm_mix

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
