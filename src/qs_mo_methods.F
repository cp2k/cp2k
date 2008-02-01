!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief collects routines that perform operations directly related to MOs
!> \note
!>      first version : most routines imported
!> \author Joost VandeVondele (2003-08)
! *****************************************************************************
MODULE qs_mo_methods
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_syrk,&
                                             cp_fm_triangular_multiply
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_retain,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply,&
                                             cp_sm_plus_fm_fm_t
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_max
  USE qs_mo_types,                     ONLY: mo_set_p_type,&
                                             mo_set_type
  USE sparse_matrix_types,             ONLY: real_matrix_type,&
                                             set_matrix
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_mo_methods'

  INTERFACE calculate_density_matrix
    MODULE PROCEDURE calculate_density_matrix_sparse,&
                     calculate_density_matrix_blacs
  END INTERFACE

  INTERFACE calculate_subspace_eigenvalues
    MODULE PROCEDURE subspace_eigenvalues_ks
  END INTERFACE

  PUBLIC :: make_basis_simple, make_basis_cholesky, make_basis_sv, make_basis_sm, &
            calculate_density_matrix, calculate_subspace_eigenvalues, &
            calculate_orthonormality, calculate_magnitude

CONTAINS

! *****************************************************************************
!> \brief returns an S-orthonormal basis v (v^T S v ==1)
!> \param v and S
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE make_basis_sm(vmatrix,ncol,matrix_s,error)
    TYPE(cp_fm_type), POINTER                :: vmatrix
    INTEGER, INTENT(IN)                      :: ncol
    TYPE(real_matrix_type), POINTER          :: matrix_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'make_basis_sm', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, i, n, ncol_global
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: overlap_vv, svmatrix

!   ---------------------------------------------------------------------------

    IF (ncol.EQ.0) RETURN

    CALL timeset(routineN,"I","",handle)

    CALL cp_fm_get_info(matrix=vmatrix,nrow_global=n,ncol_global=ncol_global,error=error)
    IF (ncol.gt.ncol_global) CALL stop_program(routineN,"wrong ncol value")

    CALL cp_fm_create(svmatrix,vmatrix%matrix_struct,"SV",error=error)
    CALL cp_sm_fm_multiply(matrix_s,vmatrix,svmatrix,ncol,error=error)

    NULLIFY(fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol,&
                             para_env=vmatrix%matrix_struct%para_env, &
                             context=vmatrix%matrix_struct%context,error=error)
    CALL cp_fm_create(overlap_vv,fm_struct_tmp,"overlap_vv",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    CALL cp_fm_gemm('T','N',ncol,ncol,n,rone,vmatrix,svmatrix,rzero, overlap_vv,error=error)
    CALL cp_fm_cholesky_decompose(overlap_vv,error=error)
    CALL cp_fm_triangular_multiply(overlap_vv,vmatrix,n_cols=ncol,side='R',invert_tr=.TRUE.,error=error)

    CALL cp_fm_release(overlap_vv,error=error)
    CALL cp_fm_release(svmatrix,error=error)

    CALL timestop(0.0_dp,handle)

    END SUBROUTINE make_basis_sm

! *****************************************************************************
!> \brief returns an S-orthonormal basis v and the corresponding matrix S*v as well
!> \param v and S*v
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE make_basis_sv(vmatrix,ncol,svmatrix,error)

    TYPE(cp_fm_type), POINTER                :: vmatrix
    INTEGER, INTENT(IN)                      :: ncol
    TYPE(cp_fm_type), POINTER                :: svmatrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'make_basis_sv', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, i, n, ncol_global
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: overlap_vv

!   ---------------------------------------------------------------------------

    IF (ncol.EQ.0) RETURN

    CALL timeset(routineN,"I","",handle)
    NULLIFY(fm_struct_tmp)

    CALL cp_fm_get_info(matrix=vmatrix,nrow_global=n,ncol_global=ncol_global,error=error)
    IF (ncol.gt.ncol_global) CALL stop_program(routineN,"wrong ncol value")

    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol,&
                             para_env=vmatrix%matrix_struct%para_env, &
                             context=vmatrix%matrix_struct%context,error=error)
    CALL cp_fm_create(overlap_vv,fm_struct_tmp,"overlap_vv",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    CALL cp_fm_gemm('T','N',ncol,ncol,n,rone,vmatrix,svmatrix,rzero, overlap_vv,error=error)
    CALL cp_fm_cholesky_decompose(overlap_vv,error=error)
    CALL cp_fm_triangular_multiply(overlap_vv,vmatrix,n_cols=ncol,side='R',invert_tr=.TRUE.,error=error)
    CALL cp_fm_triangular_multiply(overlap_vv,svmatrix,n_cols=ncol,side='R',invert_tr=.TRUE.,error=error)

    CALL cp_fm_release(overlap_vv,error=error)

    CALL timestop(0.0_dp,handle)

    END SUBROUTINE make_basis_sv

! *****************************************************************************
!> \brief return a set of S orthonormal vectors (C^T S C == 1) where
!>      the cholesky decomposed form of S is passed as an argument
!> \param ortho cholesky decomposed S matrix
!> \note
!>      if the cholesky decomposed S matrix is not available
!>      use make_basis_sm since this is much faster than computing the
!>      cholesky decomposition of S
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE make_basis_cholesky(vmatrix,ncol,ortho,error)

    TYPE(cp_fm_type), POINTER                :: vmatrix
    INTEGER, INTENT(IN)                      :: ncol
    TYPE(cp_fm_type), POINTER                :: ortho
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'make_basis_cholesky', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, n, ncol_global
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: overlap_vv

!   ---------------------------------------------------------------------------

    IF (ncol.EQ.0) RETURN

    CALL timeset(routineN,handle)
    NULLIFY(fm_struct_tmp)

    CALL cp_fm_get_info(matrix=vmatrix,nrow_global=n,ncol_global=ncol_global,error=error)
    IF (ncol.gt.ncol_global) CALL stop_program(routineN,"wrong ncol value")

    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol,&
                             para_env=vmatrix%matrix_struct%para_env, &
                             context=vmatrix%matrix_struct%context,error=error)
    CALL cp_fm_create(overlap_vv,fm_struct_tmp,"overlap_vv",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    CALL cp_fm_triangular_multiply(ortho,vmatrix,n_cols=ncol,error=error)
    CALL cp_fm_syrk('U','T',n,rone,vmatrix,1,1,rzero,overlap_vv,error=error)
    CALL cp_fm_cholesky_decompose(overlap_vv,error=error)
    CALL cp_fm_triangular_multiply(overlap_vv,vmatrix,n_cols=ncol,side='R',invert_tr=.TRUE.,error=error)
    CALL cp_fm_triangular_multiply(ortho,vmatrix,n_cols=ncol,invert_tr=.TRUE.,error=error)

    CALL cp_fm_release(overlap_vv,error=error)

    CALL timestop(0.0_dp,handle)

    END SUBROUTINE make_basis_cholesky

! *****************************************************************************
!> \brief given a set of vectors, return an orthogonal (C^T C == 1) set
!>      spanning the same space (notice, only for cases where S==1)
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE make_basis_simple(vmatrix,ncol,error)

    TYPE(cp_fm_type), POINTER                :: vmatrix
    INTEGER, INTENT(IN)                      :: ncol
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'make_basis_simple', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, n, ncol_global
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: overlap_vv

!   ---------------------------------------------------------------------------

    IF (ncol.EQ.0) RETURN

    CALL timeset(routineN,handle)

    NULLIFY(fm_struct_tmp)

    CALL cp_fm_get_info(matrix=vmatrix,nrow_global=n,ncol_global=ncol_global,error=error)
    IF (ncol.gt.ncol_global) CALL stop_program(routineN,"wrong ncol value")

    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol,&
                             para_env=vmatrix%matrix_struct%para_env, &
                             context=vmatrix%matrix_struct%context,error=error)
    CALL cp_fm_create(overlap_vv,fm_struct_tmp,"overlap_vv",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    CALL cp_fm_gemm('T','N',ncol,ncol,n,rone,vmatrix,vmatrix,rzero, overlap_vv,error=error)
    CALL cp_fm_cholesky_decompose(overlap_vv,error=error)
    CALL cp_fm_triangular_multiply(overlap_vv,vmatrix,n_cols=ncol,side='R',invert_tr=.TRUE.,error=error)

    CALL cp_fm_release(overlap_vv,error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE make_basis_simple

! *****************************************************************************
  SUBROUTINE calculate_density_matrix_blacs(mo_set,density_matrix,error)

!   Purpose: Calculate the density matrix from the MO coeff and the
!            MO occupation numbers.

!   History: - Creation (01.05.2001,MK)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: density_matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'calculate_density_matrix_blacs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, imo

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    CALL cp_fm_syrk("U","N",mo_set%lfomo-1,mo_set%maxocc,mo_set%mo_coeff,&
                    1,1,0.0_dp,density_matrix,error=error)

!   *** Add the density contribution from the fractional occupied MOs ***

    DO imo=mo_set%lfomo,mo_set%homo
      CALL cp_fm_syrk("U","N",1,mo_set%occupation_numbers(imo),&
                      mo_set%mo_coeff,1,imo,1.0_dp,density_matrix,error=error)

    END DO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_density_matrix_blacs

! *****************************************************************************
  SUBROUTINE calculate_density_matrix_sparse(mo_set,density_matrix,error)

!   Purpose: Calculate the density matrix

!   History: - Creation (06.2002, Joost VandeVondele)
!            - Fractional occupied orbitals (MK)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(real_matrix_type), POINTER          :: density_matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'calculate_density_matrix_sparse', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: alpha
    TYPE(cp_fm_type), POINTER                :: fm_tmp

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    CALL set_matrix(density_matrix,0.0_dp)

    IF ( .NOT. mo_set%uniform_occupation ) THEN ! not all orbitals 1..homo are equally occupied
       NULLIFY(fm_tmp)
       CALL cp_fm_create(fm_tmp,mo_set%mo_coeff%matrix_struct,error=error)
       CALL cp_fm_to_fm(mo_set%mo_coeff,fm_tmp,error=error)
       CALL cp_fm_column_scale(fm_tmp,mo_set%occupation_numbers(1:mo_set%homo))
       alpha=1.0_dp

    ELSE
       alpha=mo_set%maxocc
       fm_tmp => mo_set%mo_coeff
       CALL cp_fm_retain(fm_tmp,error=error)
    ENDIF

    CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix,&
                            matrix_v=mo_set%mo_coeff,&
                            matrix_g=fm_tmp,&
                            ncol=mo_set%homo,&
                            alpha=alpha,error=error)

    CALL cp_fm_release(fm_tmp,error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_density_matrix_sparse

! *****************************************************************************
!> \brief computes ritz values of a set of orbitals given a ks_matrix
!>      rotates the orbitals into eigenstates depending on do_rotation
!>      writes the evals to the screen depending on ionode/scr
!> \param orbitals S-orthonormal orbitals
!> \param ks_matrix Kohn-Sham matrix
!> \param para_env of the Kohn-Sham matrix
!> \param evals_arg optional, filled with the evals
!> \param ionode /scr  : if present write to unit scr where ionode
!> \param do_rotation optional rotate orbitals if .TRUE. (default)
!> \param co_rotate an optional set of orbitals rotated by the same amout
!> \par History
!>      08.2004 documented and added do_rotation [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE subspace_eigenvalues_ks(orbitals,ks_matrix,evals_arg,para_env,ionode,scr,do_rotation,co_rotate,error)

!   Purpose: Calculates the eigenvalues (ritz values) for a set of orbitals
!            and rotates them to make them eigenvectors
!            if evals_arg is present its filled
!            if scr and ionode are present they go to the screen as well.

    TYPE(cp_fm_type), POINTER                :: orbitals
    TYPE(real_matrix_type), POINTER          :: ks_matrix
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: evals_arg
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN), OPTIONAL            :: ionode
    INTEGER, INTENT(IN), OPTIONAL            :: scr
    LOGICAL, INTENT(IN), OPTIONAL            :: do_rotation
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: co_rotate
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'subspace_eigenvalues_ks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, ncol_global, &
                                                nrow_global
    LOGICAL                                  :: do_rotation_local
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: e_vectors, h_block, &
                                                weighted_vectors

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    do_rotation_local=.TRUE.
    IF (PRESENT(do_rotation)) do_rotation_local=do_rotation

    NULLIFY(weighted_vectors,h_block,e_vectors,fm_struct_tmp)
    CALL cp_fm_get_info(matrix=orbitals, &
                                 ncol_global=ncol_global, &
                                 nrow_global=nrow_global,error=error)

    IF (ncol_global.GT.0) THEN

        ALLOCATE(evals(ncol_global))

        CALL cp_fm_create(weighted_vectors,orbitals%matrix_struct,"weighted_vectors",error=error)
        CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol_global,  ncol_global=ncol_global, &
                                                para_env=orbitals%matrix_struct%para_env, &
                                                context=orbitals%matrix_struct%context,error=error)
        CALL cp_fm_create(h_block,fm_struct_tmp, name="h block",error=error)
        CALL cp_fm_create(e_vectors,fm_struct_tmp, name="e vectors",error=error)
        CALL cp_fm_struct_release(fm_struct_tmp,error=error)

        ! h subblock and diag
        CALL cp_sm_fm_multiply(ks_matrix,orbitals,weighted_vectors, ncol_global,error=error)

        CALL cp_fm_gemm('T','N',ncol_global,ncol_global,nrow_global,1.0_dp, &
                        orbitals,weighted_vectors,0.0_dp,h_block,error=error)

        CALL cp_fm_syevd(h_block,e_vectors,evals,error=error)

        ! rotate the orbitals
        IF (do_rotation_local) THEN
            CALL cp_fm_gemm('N','N',nrow_global,ncol_global,ncol_global,1.0_dp, &
                        orbitals,e_vectors,0.0_dp,weighted_vectors,error=error)
            CALL cp_fm_to_fm(weighted_vectors,orbitals,error=error)
            IF (PRESENT(co_rotate)) THEN
              IF (ASSOCIATED(co_rotate)) THEN
                CALL cp_fm_gemm('N','N',nrow_global,ncol_global,ncol_global,1.0_dp, &
                        co_rotate,e_vectors,0.0_dp,weighted_vectors,error=error)
                CALL cp_fm_to_fm(weighted_vectors,co_rotate,error=error)
              ENDIF
            ENDIF
        ENDIF

        ! give output
        IF (PRESENT(evals_arg)) THEN
          evals_arg(:)=evals(:)
        ENDIF

        IF (PRESENT(ionode) .OR. PRESENT(scr)) THEN
          IF (.NOT. PRESENT(ionode)) CALL stop_program("subspace evals","IONODE?")
          IF (.NOT. PRESENT(scr)) CALL stop_program("subspace evals","SCR?")
          IF (ionode) THEN
             DO i=1,ncol_global,4
                j=MIN(3,ncol_global-i)
                SELECT CASE (j)
                CASE(3)
                 WRITE(scr,'(1X,4F16.8)') evals(i:i+j)
                CASE(2)
                 WRITE(scr,'(1X,3F16.8)') evals(i:i+j)
                CASE(1)
                 WRITE(scr,'(1X,2F16.8)') evals(i:i+j)
                CASE(0)
                 WRITE(scr,'(1X,1F16.8)') evals(i:i+j)
                END SELECT
             ENDDO
          ENDIF
        ENDIF

        CALL cp_fm_release(weighted_vectors,error=error)
        CALL cp_fm_release(h_block,error=error)
        CALL cp_fm_release(e_vectors,error=error)

        DEALLOCATE(evals)

    ENDIF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE subspace_eigenvalues_ks

! computes the effective orthonormality of a set of mos given an s-matrix
! orthonormality is the max deviation from unity of the C^T S C
! *****************************************************************************
  SUBROUTINE calculate_orthonormality(orthonormality,mo_array,matrix_s,error)
    REAL(KIND=dp)                            :: orthonormality
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: matrix_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_orthonormality', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, j, k, n, &
                                                ncol_local, nrow_local, nspin
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: alpha, max_alpha
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: overlap, svec

   NULLIFY(tmp_fm_struct,svec,overlap)

   CALL timeset(routineN,"I","",handle)

   nspin=SIZE(mo_array)
   max_alpha=0.0_dp

   DO ispin=1,nspin
      IF (PRESENT(matrix_s)) THEN
        ! get S*C
        CALL cp_fm_create(svec,mo_array(ispin)%mo_set%mo_coeff%matrix_struct,error=error)
        CALL cp_fm_get_info(mo_array(ispin)%mo_set%mo_coeff,&
             nrow_global=n,ncol_global=k,error=error)
        CALL cp_sm_fm_multiply(matrix_s,mo_array(ispin)%mo_set%mo_coeff,&
             svec,k,error=error)
        ! get C^T (S*C)
        CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=k,ncol_global=k, &
             para_env=mo_array(ispin)%mo_set%mo_coeff%matrix_struct%para_env, &
             context=mo_array(ispin)%mo_set%mo_coeff%matrix_struct%context,error=error)
        CALL cp_fm_create(overlap,tmp_fm_struct,error=error)
        CALL cp_fm_struct_release(tmp_fm_struct,error=error)
        CALL cp_fm_gemm('T','N',k,k,n,1.0_dp, mo_array(ispin)%mo_set%mo_coeff,&
             svec,0.0_dp,overlap,error=error)
        CALL cp_fm_release(svec,error=error)
      ELSE
        ! orthogonal basis C^T C
        CALL cp_fm_get_info(mo_array(ispin)%mo_set%mo_coeff,&
             nrow_global=n,ncol_global=k,error=error)
        CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=k,ncol_global=k, &
             para_env=mo_array(ispin)%mo_set%mo_coeff%matrix_struct%para_env, &
             context=mo_array(ispin)%mo_set%mo_coeff%matrix_struct%context,error=error)
        CALL cp_fm_create(overlap,tmp_fm_struct,error=error)
        CALL cp_fm_struct_release(tmp_fm_struct,error=error)
        CALL cp_fm_gemm('T','N',k,k,n,1.0_dp, mo_array(ispin)%mo_set%mo_coeff,&
             mo_array(ispin)%mo_set%mo_coeff,0.0_dp,overlap,error=error)
      ENDIF
      CALL cp_fm_get_info(overlap,nrow_local=nrow_local,ncol_local=ncol_local, &
           row_indices=row_indices,col_indices=col_indices,error=error)
      DO i=1,nrow_local
         DO j=1,ncol_local
            alpha=overlap%local_data(i,j)
            IF (row_indices(i).eq.col_indices(j)) alpha=alpha-1.0_dp
            max_alpha=MAX(max_alpha,ABS(alpha))
         ENDDO
      ENDDO
      CALL cp_fm_release(overlap,error=error)
    ENDDO
    CALL mp_max(max_alpha,mo_array(1)%mo_set%mo_coeff%matrix_struct%para_env%group)
    orthonormality=max_alpha
    ! write(6,*) "max deviation from orthonormalization ",orthonormality

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_orthonormality

! computes the minimum/maximum magnitudes of C^T C. This could be useful
! to detect problems in the case of nearly singular overlap matrices.
! in this case, we expect the ratio of min/max to be large
! this routine is only similar to mo_orthonormality if S==1
! *****************************************************************************
  SUBROUTINE calculate_magnitude(mo_array,mo_mag_min,mo_mag_max,error)
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    REAL(KIND=dp)                            :: mo_mag_min, mo_mag_max
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_magnitude', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, k, n, nspin
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: evals
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: evecs, overlap

   NULLIFY(tmp_fm_struct,overlap)

   CALL timeset(routineN,"I","",handle)

   nspin=SIZE(mo_array)
   mo_mag_min= HUGE(0.0_dp)
   mo_mag_max=-HUGE(0.0_dp)
   DO ispin=1,nspin
      CALL cp_fm_get_info(mo_array(ispin)%mo_set%mo_coeff,&
             nrow_global=n,ncol_global=k,error=error)
      ALLOCATE(evals(k))
      CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=k,ncol_global=k, &
             para_env=mo_array(ispin)%mo_set%mo_coeff%matrix_struct%para_env, &
             context=mo_array(ispin)%mo_set%mo_coeff%matrix_struct%context,error=error)
      CALL cp_fm_create(overlap,tmp_fm_struct,error=error)
      CALL cp_fm_create(evecs,tmp_fm_struct,error=error)
      CALL cp_fm_struct_release(tmp_fm_struct,error=error)
      CALL cp_fm_gemm('T','N',k,k,n,1.0_dp, mo_array(ispin)%mo_set%mo_coeff,&
             mo_array(ispin)%mo_set%mo_coeff,0.0_dp,overlap,error=error)
      CALL cp_fm_syevd(overlap,evecs,evals,error)
      mo_mag_min=MIN(MINVAL(evals),mo_mag_min)
      mo_mag_max=MAX(MAXVAL(evals),mo_mag_max)
      CALL cp_fm_release(overlap,error=error)
      CALL cp_fm_release(evecs,error=error)
      DEALLOCATE(evals)
    ENDDO
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_magnitude

END MODULE qs_mo_methods
