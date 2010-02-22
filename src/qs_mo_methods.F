!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief collects routines that perform operations directly related to MOs
!> \note
!>      first version : most routines imported
!> \author Joost VandeVondele (2003-08)
! *****************************************************************************
MODULE qs_mo_methods
  USE array_types,                     ONLY: array_i1d_obj,&
                                             array_release
  USE cp_dbcsr_diag,                   ONLY: dbcsr_syevd,&
                                             dbcsr_syevx
  USE cp_dbcsr_operations,             ONLY: &
       copy_dbcsr_to_fm, copy_fm_to_dbcsr, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_finalize, cp_dbcsr_init, cp_dbcsr_init_p, cp_dbcsr_multiply, &
       cp_dbcsr_release, cp_dbcsr_release_p, cp_dbcsr_scale_by_vector, &
       cp_dbcsr_set, dbcsr_plus_fm_fm_t, dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_syrk,&
                                             cp_fm_triangular_multiply
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
  USE cp_fm_diag,                      ONLY: cp_fm_power,&
                                             cp_fm_syevd,&
                                             cp_fm_syevx
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE dbcsr_dist_operations,           ONLY: create_bl_distribution
  USE dbcsr_methods,                   ONLY: dbcsr_distribution,&
                                             dbcsr_distribution_mp,&
                                             dbcsr_distribution_new,&
                                             dbcsr_distribution_release,&
                                             dbcsr_get_info,&
                                             dbcsr_mp_npcols,&
                                             dbcsr_mp_nprows
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_obj,&
                                             dbcsr_type_normal,&
                                             dbcsr_type_real_default
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_max
  USE qs_mo_types,                     ONLY: mo_set_p_type,&
                                             mo_set_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_mo_methods'

  PUBLIC :: make_basis_simple, make_basis_cholesky, make_basis_sv, make_basis_sm, &
            make_basis_lowdin, calculate_density_matrix, calculate_subspace_eigenvalues, &
            calculate_orthonormality, calculate_magnitude

  INTERFACE calculate_density_matrix
    MODULE PROCEDURE calculate_density_matrix_blacs
    MODULE PROCEDURE calculate_dm_sparse
  END INTERFACE

  INTERFACE calculate_subspace_eigenvalues
    MODULE PROCEDURE subspace_eigenvalues_ks_fm
    MODULE PROCEDURE subspace_eigenvalues_ks_dbcsr
  END INTERFACE

  INTERFACE make_basis_sv
     MODULE PROCEDURE make_basis_sv_fm
     MODULE PROCEDURE make_basis_sv_dbcsr
  END INTERFACE

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
    TYPE(dbcsr_obj), POINTER                 :: matrix_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'make_basis_sm', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, n, ncol_global
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: overlap_vv, svmatrix

    IF (ncol.EQ.0) RETURN

    CALL timeset(routineN,handle)

    CALL cp_fm_get_info(matrix=vmatrix,nrow_global=n,ncol_global=ncol_global,error=error)
    IF (ncol.gt.ncol_global) CALL stop_program(routineN,"wrong ncol value")

    CALL cp_fm_create(svmatrix,vmatrix%matrix_struct,"SV",error=error)
    CALL dbcsr_sm_fm_multiply(matrix_s,vmatrix,svmatrix,ncol,error=error)

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

    CALL timestop(handle)

  END SUBROUTINE make_basis_sm

! *****************************************************************************
!> \brief returns an S-orthonormal basis v and the corresponding matrix S*v as well
!> \param v and S*v
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE make_basis_sv_fm(vmatrix,ncol,svmatrix,error)

    TYPE(cp_fm_type), POINTER                :: vmatrix
    INTEGER, INTENT(IN)                      :: ncol
    TYPE(cp_fm_type), POINTER                :: svmatrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'make_basis_sv_fm', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, n, ncol_global
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: overlap_vv

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

    CALL cp_fm_gemm('T','N',ncol,ncol,n,rone,vmatrix,svmatrix,rzero, overlap_vv,error=error)
    CALL cp_fm_cholesky_decompose(overlap_vv,error=error)
    CALL cp_fm_triangular_multiply(overlap_vv,vmatrix,n_cols=ncol,side='R',invert_tr=.TRUE.,error=error)
    CALL cp_fm_triangular_multiply(overlap_vv,svmatrix,n_cols=ncol,side='R',invert_tr=.TRUE.,error=error)

    CALL cp_fm_release(overlap_vv,error=error)

    CALL timestop(handle)

    END SUBROUTINE make_basis_sv_fm

  SUBROUTINE make_basis_sv_dbcsr(vmatrix,ncol,svmatrix,para_env,blacs_env,error)

    TYPE(dbcsr_obj)                          :: vmatrix
    INTEGER, INTENT(IN)                      :: ncol
    TYPE(dbcsr_obj)                          :: svmatrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'make_basis_sv_dbcsr', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, n, ncol_global
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_svmatrix, fm_vmatrix, &
                                                overlap_vv

    IF (ncol.EQ.0) RETURN

    CALL timeset(routineN,handle)

    !CALL cp_fm_get_info(matrix=vmatrix,nrow_global=n,ncol_global=ncol_global,error=error)
    CALL dbcsr_get_info(vmatrix,nfullrows_total=n,nfullcols_total=ncol_global)
    IF (ncol.gt.ncol_global) CALL stop_program(routineN,"wrong ncol value")


    CALL cp_fm_struct_create(fm_struct_tmp,context=blacs_env,nrow_global=ncol,&
         ncol_global=ncol,para_env=para_env,error=error)
    CALL cp_fm_create(overlap_vv,fm_struct_tmp,name="fm_overlap_vv",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    CALL cp_fm_struct_create(fm_struct_tmp,context=blacs_env,nrow_global=n,&
         ncol_global=ncol_global,para_env=para_env,error=error)
    CALL cp_fm_create(fm_vmatrix,fm_struct_tmp,name="fm_vmatrix",error=error)
    CALL cp_fm_create(fm_svmatrix,fm_struct_tmp,name="fm_svmatrix",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    CALL copy_dbcsr_to_fm(vmatrix,fm_vmatrix,error=error)
    CALL copy_dbcsr_to_fm(svmatrix,fm_svmatrix,error=error)

    CALL cp_fm_gemm('T','N',ncol,ncol,n,rone,fm_vmatrix,fm_svmatrix,rzero, overlap_vv,error=error)
    CALL cp_fm_cholesky_decompose(overlap_vv,error=error)
    CALL cp_fm_triangular_multiply(overlap_vv,fm_vmatrix,n_cols=ncol,side='R',invert_tr=.TRUE.,error=error)
    CALL cp_fm_triangular_multiply(overlap_vv,fm_svmatrix,n_cols=ncol,side='R',invert_tr=.TRUE.,error=error)

    CALL copy_fm_to_dbcsr(fm_vmatrix,vmatrix, error=error)
    CALL copy_fm_to_dbcsr(fm_svmatrix,svmatrix, error=error)

    CALL cp_fm_release(overlap_vv,error=error)
    CALL cp_fm_release(fm_vmatrix,error=error)
    CALL cp_fm_release(fm_svmatrix,error=error)

    CALL timestop(handle)

    END SUBROUTINE make_basis_sv_dbcsr

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
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: overlap_vv

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

    CALL timestop(handle)

    END SUBROUTINE make_basis_cholesky

! *****************************************************************************
!> \brief return a set of S orthonormal vectors (C^T S C == 1) where
!>      a Loedwin transformation is applied to keep the rotated vectors as close
!>      as possible to the original ones
!> \param 
!> \note
!> \par History
!>      05.2009 created [MI]
! *****************************************************************************
  SUBROUTINE make_basis_lowdin(vmatrix,ncol,matrix_s,error)

    TYPE(cp_fm_type), POINTER                :: vmatrix
    INTEGER, INTENT(IN)                      :: ncol
    TYPE(dbcsr_obj), POINTER                 :: matrix_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'make_basis_lowdin', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, n, ncol_global, ndep
    REAL(dp)                                 :: threshold
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: csc, sc, work

    IF (ncol.EQ.0) RETURN

    CALL timeset(routineN,handle)
    NULLIFY(fm_struct_tmp)
    threshold = 1.0E-7_dp
    CALL cp_fm_get_info(matrix=vmatrix,nrow_global=n,ncol_global=ncol_global,error=error)
    IF (ncol.gt.ncol_global) CALL stop_program(routineN,"wrong ncol value")

    CALL cp_fm_create(sc,vmatrix%matrix_struct,"SC",error=error)
    CALL dbcsr_sm_fm_multiply(matrix_s,vmatrix,sc,ncol,error=error)

    NULLIFY(fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol,&
                             para_env=vmatrix%matrix_struct%para_env, &
                             context=vmatrix%matrix_struct%context,error=error)
    CALL cp_fm_create(csc,fm_struct_tmp,"csc",error=error)
    CALL cp_fm_create(work,fm_struct_tmp,"work",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    CALL cp_fm_gemm('T','N',ncol,ncol,n,rone,vmatrix,sc,rzero, csc,error=error)
    CALL cp_fm_power(csc,work,-0.5_dp, threshold,ndep,error=error)
    CALL cp_fm_gemm('N','N',n,ncol,ncol,rone,vmatrix,csc,rzero,sc,error=error)
    CALL cp_fm_to_fm(sc, vmatrix, ncol, 1,1)

    CALL cp_fm_release(csc,error=error)
    CALL cp_fm_release(sc,error=error)
    CALL cp_fm_release(work,error=error)

    CALL timestop(handle)

  END SUBROUTINE make_basis_lowdin

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
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: overlap_vv

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

    CALL timestop(handle)

  END SUBROUTINE make_basis_simple

! *****************************************************************************
!> \brief  Calculate the density matrix from the MO coeff and the
!>         MO occupation numbers. 
!> \author  MK
!> \date    01.05.2001
!> \version 1.0
! *****************************************************************************
  SUBROUTINE calculate_density_matrix_blacs(mo_set,density_matrix,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: density_matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'calculate_density_matrix_blacs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, imo

    CALL timeset(routineN,handle)

    CALL cp_fm_syrk("U","N",mo_set%lfomo-1,mo_set%maxocc,mo_set%mo_coeff,&
                    1,1,0.0_dp,density_matrix,error=error)

!   *** Add the density contribution from the fractional occupied MOs ***

    DO imo=mo_set%lfomo,mo_set%homo
      CALL cp_fm_syrk("U","N",1,mo_set%occupation_numbers(imo),&
                      mo_set%mo_coeff,1,imo,1.0_dp,density_matrix,error=error)

    END DO

    CALL timestop(handle)

  END SUBROUTINE calculate_density_matrix_blacs

! *****************************************************************************
!> \brief   Calculate the density matrix
!> \author  Joost VandeVondele
!> \date    06.2002
!> \par History
!>       - Fractional occupied orbitals (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE calculate_dm_sparse(mo_set,density_matrix,use_dbcsr,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(dbcsr_obj), POINTER                 :: density_matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: use_dbcsr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_dm_sparse', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: my_use_dbcsr
    REAL(KIND=dp)                            :: alpha
    TYPE(cp_fm_type), POINTER                :: fm_tmp
    TYPE(dbcsr_obj)                          :: dbcsr_tmp

    CALL timeset(routineN,handle)

    my_use_dbcsr = .FALSE.
    IF(PRESENT(use_dbcsr)) my_use_dbcsr = use_dbcsr
    IF(my_use_dbcsr) THEN
       IF(.NOT.ASSOCIATED(mo_set%mo_coeff_b)) CALL stop_program(routineN,"mo_coeff_b NOT ASSOCIATED")
    ENDIF


    CALL cp_dbcsr_set(density_matrix,0.0_dp,error=error)

    IF ( .NOT. mo_set%uniform_occupation ) THEN ! not all orbitals 1..homo are equally occupied

       IF(my_use_dbcsr) THEN
          CALL cp_dbcsr_init(dbcsr_tmp,error=error)
          CALL cp_dbcsr_copy(dbcsr_tmp,mo_set%mo_coeff_b,error=error)
          CALL cp_dbcsr_scale_by_vector(dbcsr_tmp,mo_set%occupation_numbers(1:mo_set%homo),&
               side='right',error=error)
          CALL cp_dbcsr_multiply("N", "T", 1.0_dp, mo_set%mo_coeff_b, dbcsr_tmp,&
               1.0_dp, density_matrix, retain_sparsity=.TRUE.,&
               last_k = mo_set%homo, error=error)
          CALL cp_dbcsr_release(dbcsr_tmp, error=error)
       ELSE
          NULLIFY(fm_tmp)
          CALL cp_fm_create(fm_tmp,mo_set%mo_coeff%matrix_struct,error=error)
          CALL cp_fm_to_fm(mo_set%mo_coeff,fm_tmp,error=error)
          CALL cp_fm_column_scale(fm_tmp,mo_set%occupation_numbers(1:mo_set%homo))
          alpha=1.0_dp
          CALL dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix,&
                                  matrix_v=mo_set%mo_coeff,&
                                  matrix_g=fm_tmp,&
                                  ncol=mo_set%homo,&
                                  alpha=alpha,error=error)
          CALL cp_fm_release(fm_tmp,error=error)
       ENDIF
    ELSE
       IF(my_use_dbcsr) THEN
          CALL cp_dbcsr_multiply("N", "T", mo_set%maxocc, mo_set%mo_coeff_b, mo_set%mo_coeff_b,&
               1.0_dp, density_matrix, retain_sparsity=.TRUE.,&
               last_k = mo_set%homo, error=error)
       ELSE
          alpha=mo_set%maxocc
          CALL dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix,&
                                  matrix_v=mo_set%mo_coeff,&
                                  ncol=mo_set%homo,&
                                  alpha=alpha,error=error)
       ENDIF
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE calculate_dm_sparse

! *****************************************************************************
!> \brief computes ritz values of a set of orbitals given a ks_matrix
!>      rotates the orbitals into eigenstates depending on do_rotation
!>      writes the evals to the screen depending on ionode/scr
!> \param orbitals S-orthonormal orbitals
!> \param ks_matrix Kohn-Sham matrix
!> \param para_env of the Kohn-Sham matrix
!> \param evals_arg optional, filled with the evals
!> \param ionode /scr  : if present write to unit scr where ionode
!> \param do_rotation optional rotate orbitals (default=.TRUE.)
!>        note that rotating the orbitals is slower
!> \param co_rotate an optional set of orbitals rotated by the same rotation matrix
!> \par History
!>      08.2004 documented and added do_rotation [Joost VandeVondele]
!>      09.2008 only compute eigenvalues if rotation is not needed
! *****************************************************************************
  SUBROUTINE subspace_eigenvalues_ks_fm(orbitals,ks_matrix,evals_arg,ionode,scr,&
       do_rotation,co_rotate,co_rotate_dbcsr,error)

    TYPE(cp_fm_type), POINTER                :: orbitals
    TYPE(dbcsr_obj), POINTER                 :: ks_matrix
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: evals_arg
    LOGICAL, INTENT(IN), OPTIONAL            :: ionode
    INTEGER, INTENT(IN), OPTIONAL            :: scr
    LOGICAL, INTENT(IN), OPTIONAL            :: do_rotation
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: co_rotate
    TYPE(dbcsr_obj), OPTIONAL, POINTER       :: co_rotate_dbcsr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'subspace_eigenvalues_ks_fm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, ncol_global, &
                                                nrow_global
    LOGICAL                                  :: compute_evecs, &
                                                do_rotation_local
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: e_vectors, h_block, &
                                                weighted_vectors, &
                                                weighted_vectors2

    CALL timeset(routineN,handle)

    do_rotation_local=.TRUE.
    IF (PRESENT(do_rotation)) do_rotation_local=do_rotation

    NULLIFY(weighted_vectors,weighted_vectors2,h_block,e_vectors,fm_struct_tmp)
    CALL cp_fm_get_info(matrix=orbitals, &
                                 ncol_global=ncol_global, &
                                 nrow_global=nrow_global,error=error)

    IF (do_rotation_local) THEN
       compute_evecs=.TRUE.
    ELSE
       ! this would be the logical choice if syevx computing only evals were faster than syevd computing evecs and evals.
       compute_evecs=.FALSE.
       ! this is not the case, so lets compute evecs always
       compute_evecs=.TRUE.
    ENDIF

    IF (ncol_global.GT.0) THEN

        ALLOCATE(evals(ncol_global))

        CALL cp_fm_create(weighted_vectors,orbitals%matrix_struct,"weighted_vectors",error=error)
        CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol_global,  ncol_global=ncol_global, &
                                                para_env=orbitals%matrix_struct%para_env, &
                                                context=orbitals%matrix_struct%context,error=error)
        CALL cp_fm_create(h_block,fm_struct_tmp, name="h block",error=error)
        IF (compute_evecs) THEN
           CALL cp_fm_create(e_vectors,fm_struct_tmp, name="e vectors",error=error)
        ENDIF
        CALL cp_fm_struct_release(fm_struct_tmp,error=error)

        ! h subblock and diag
        CALL dbcsr_sm_fm_multiply(ks_matrix,orbitals,weighted_vectors, ncol_global,error=error)

        CALL cp_fm_gemm('T','N',ncol_global,ncol_global,nrow_global,1.0_dp, &
                        orbitals,weighted_vectors,0.0_dp,h_block,error=error)

        ! if eigenvectors are required, go for syevd, otherwise only compute eigenvalues
        IF (compute_evecs) THEN
            CALL cp_fm_syevd(h_block,e_vectors,evals,error=error)
        ELSE
            CALL cp_fm_syevx(h_block,eigenvalues=evals,error=error)
        ENDIF

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
            IF(PRESENT(co_rotate_dbcsr)) THEN
               IF(ASSOCIATED(co_rotate_dbcsr)) THEN
                  CALL cp_fm_create(weighted_vectors2,orbitals%matrix_struct,"weighted_vectors",error=error)
                  CALL copy_dbcsr_to_fm(co_rotate_dbcsr,weighted_vectors2,error)
                  CALL cp_fm_gemm('N','N',nrow_global,ncol_global,ncol_global,1.0_dp, &
                       weighted_vectors2,e_vectors,0.0_dp,weighted_vectors,error=error)
                  CALL copy_fm_to_dbcsr(weighted_vectors,co_rotate_dbcsr,error=error)
                  CALL cp_fm_release(weighted_vectors2,error=error)
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
        IF (compute_evecs) THEN
           CALL cp_fm_release(e_vectors,error=error)
        ENDIF

        DEALLOCATE(evals)

    ENDIF

    CALL timestop(handle)

  END SUBROUTINE subspace_eigenvalues_ks_fm


  SUBROUTINE subspace_eigenvalues_ks_dbcsr(orbitals,ks_matrix,evals_arg,ionode,scr,&
       do_rotation,co_rotate,para_env,blacs_env,error)

    TYPE(dbcsr_obj), POINTER                 :: orbitals, ks_matrix
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: evals_arg
    LOGICAL, INTENT(IN), OPTIONAL            :: ionode
    INTEGER, INTENT(IN), OPTIONAL            :: scr
    LOGICAL, INTENT(IN), OPTIONAL            :: do_rotation
    TYPE(dbcsr_obj), OPTIONAL, POINTER       :: co_rotate
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'subspace_eigenvalues_ks_dbcsr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, ncol_global, &
                                                nrow_global
    LOGICAL                                  :: compute_evecs, &
                                                do_rotation_local
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: evals
    TYPE(array_i1d_obj)                      :: col_blk_size, col_dist, &
                                                row_blk_size, row_dist
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(dbcsr_obj), POINTER                 :: e_vectors, h_block, &
                                                weighted_vectors

    CALL timeset(routineN,handle)

    do_rotation_local=.TRUE.
    IF (PRESENT(do_rotation)) do_rotation_local=do_rotation

    NULLIFY(e_vectors, h_block, weighted_vectors)

    CALL dbcsr_get_info(matrix=orbitals, &
         nfullcols_total=ncol_global, &
         nfullrows_total=nrow_global)

    IF (do_rotation_local) THEN
       compute_evecs=.TRUE.
    ELSE
       ! this would be the logical choice if syevx computing only evals were faster than syevd computing evecs and evals.
       compute_evecs=.FALSE.
       ! this is not the case, so lets compute evecs always
       compute_evecs=.TRUE.
    ENDIF

    IF (ncol_global.GT.0) THEN

        ALLOCATE(evals(ncol_global))

        CALL cp_dbcsr_init_p(weighted_vectors,error=error)
        CALL cp_dbcsr_copy(weighted_vectors,orbitals,name="weighted_vectors",error=error)

        CALL create_bl_distribution (col_dist, col_blk_size, ncol_global, &
             dbcsr_mp_npcols(dbcsr_distribution_mp(dbcsr_distribution(orbitals))))
        CALL create_bl_distribution (row_dist, row_blk_size, ncol_global, &
             dbcsr_mp_nprows(dbcsr_distribution_mp(dbcsr_distribution(orbitals))))
        CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp (dbcsr_distribution(orbitals)),&
             row_dist, col_dist)
        CALL array_release (row_dist); CALL array_release (col_dist)

        !CALL cp_fm_create(weighted_vectors,orbitals%matrix_struct,"weighted_vectors",error=error)
        !CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol_global,  ncol_global=ncol_global, &
        !                                        para_env=orbitals%matrix_struct%para_env, &
        !                                        context=orbitals%matrix_struct%context,error=error)
        CALL cp_dbcsr_init_p(h_block,error=error)
        CALL cp_dbcsr_create(h_block, "h_block", dist, dbcsr_type_normal,&
             row_blk_size, col_blk_size, 0, 0, dbcsr_type_real_default, special=.TRUE.,&
             error=error)
        CALL cp_dbcsr_finalize(h_block, error=error)
        !CALL cp_fm_create(h_block,fm_struct_tmp, name="h block",error=error)
        IF (compute_evecs) THEN
           CALL cp_dbcsr_init_p(e_vectors,error=error)
           CALL cp_dbcsr_create(e_vectors, "e_vectors", dist, dbcsr_type_normal,&
                row_blk_size, col_blk_size, 0, 0, dbcsr_type_real_default, special=.TRUE.,&
                error=error)
           CALL cp_dbcsr_finalize(e_vectors, error=error)
           !CALL cp_fm_create(e_vectors,fm_struct_tmp, name="e vectors",error=error)
        ENDIF
        CALL dbcsr_distribution_release (dist)
        CALL array_release (row_blk_size);CALL array_release (col_blk_size)
        !CALL cp_fm_struct_release(fm_struct_tmp,error=error)

        ! h subblock and diag
        CALL cp_dbcsr_multiply('N','N',1.0_dp,ks_matrix,orbitals,&
             0.0_dp,weighted_vectors,last_column=ncol_global,error=error)
        !CALL dbcsr_sm_fm_multiply(ks_matrix,orbitals,weighted_vectors, ncol_global,error=error)

        CALL cp_dbcsr_multiply('T','N',1.0_dp,orbitals,weighted_vectors,0.0_dp,h_block,error=error)
        !CALL cp_fm_gemm('T','N',ncol_global,ncol_global,nrow_global,1.0_dp, &
        !                orbitals,weighted_vectors,0.0_dp,h_block,error=error)

        ! if eigenvectors are required, go for syevd, otherwise only compute eigenvalues
        IF (compute_evecs) THEN
           CALL dbcsr_syevd(h_block,e_vectors,evals,para_env=para_env,blacs_env=blacs_env,error=error)
        ELSE
           CALL dbcsr_syevx(h_block,eigenvalues=evals,para_env=para_env,blacs_env=blacs_env,error=error)
        ENDIF

        ! rotate the orbitals
        IF (do_rotation_local) THEN
           CALL cp_dbcsr_multiply('N','N',1.0_dp,orbitals,e_vectors,0.0_dp,weighted_vectors,error=error)
           !CALL cp_fm_gemm('N','N',nrow_global,ncol_global,ncol_global,1.0_dp, &
           !             orbitals,e_vectors,0.0_dp,weighted_vectors,error=error)
           CALL cp_dbcsr_copy(orbitals,weighted_vectors,error=error)
           !CALL cp_fm_to_fm(weighted_vectors,orbitals,error=error)
           IF (PRESENT(co_rotate)) THEN
              IF (ASSOCIATED(co_rotate)) THEN
                 CALL cp_dbcsr_multiply('N','N',1.0_dp,co_rotate,e_vectors,0.0_dp,weighted_vectors,error=error)
                 !CALL cp_fm_gemm('N','N',nrow_global,ncol_global,ncol_global,1.0_dp, &
                 !     co_rotate,e_vectors,0.0_dp,weighted_vectors,error=error)
                 CALL cp_dbcsr_copy(co_rotate,weighted_vectors,error=error)
                 !CALL cp_fm_to_fm(weighted_vectors,co_rotate,error=error)
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

        CALL cp_dbcsr_release_p(weighted_vectors,error=error)
        CALL cp_dbcsr_release_p(h_block,error=error)
        !CALL cp_fm_release(weighted_vectors,error=error)
        !CALL cp_fm_release(h_block,error=error)
        IF (compute_evecs) THEN
           CALL cp_dbcsr_release_p(e_vectors,error=error)
           !CALL cp_fm_release(e_vectors,error=error)
        ENDIF

        DEALLOCATE(evals)

    ENDIF

    CALL timestop(handle)

  END SUBROUTINE subspace_eigenvalues_ks_dbcsr


! computes the effective orthonormality of a set of mos given an s-matrix
! orthonormality is the max deviation from unity of the C^T S C
! *****************************************************************************
  SUBROUTINE calculate_orthonormality(orthonormality,mo_array,matrix_s,error)
    REAL(KIND=dp)                            :: orthonormality
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(dbcsr_obj), OPTIONAL, POINTER       :: matrix_s
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

   CALL timeset(routineN,handle)

   nspin=SIZE(mo_array)
   max_alpha=0.0_dp

   DO ispin=1,nspin
      IF (PRESENT(matrix_s)) THEN
        ! get S*C
        CALL cp_fm_create(svec,mo_array(ispin)%mo_set%mo_coeff%matrix_struct,error=error)
        CALL cp_fm_get_info(mo_array(ispin)%mo_set%mo_coeff,&
             nrow_global=n,ncol_global=k,error=error)
        CALL dbcsr_sm_fm_multiply(matrix_s,mo_array(ispin)%mo_set%mo_coeff,&
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

    CALL timestop(handle)

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

    INTEGER                                  :: handle, ispin, k, n, nspin
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: evals
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: evecs, overlap

   NULLIFY(tmp_fm_struct,overlap)

   CALL timeset(routineN,handle)

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
    CALL timestop(handle)

  END SUBROUTINE calculate_magnitude

END MODULE qs_mo_methods
