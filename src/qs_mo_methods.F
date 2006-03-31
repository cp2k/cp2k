!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2003 2006 CP2K developers group
!-----------------------------------------------------------------------------!
!!****h* cp2k/qs_mo_methods [1.0] *
!!
!!   NAME
!!     qs_mo_methods
!!
!!   FUNCTION
!!     collects routines that perform operations directly related to MOs
!!
!!   NOTES
!!     first version : most routines imported 
!!
!!   AUTHOR
!!     Joost VandeVondele (2003-08)
!!
! ************************************************************************************

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
!!****f* qs_mo_methods/make_basis_sm*
!!
!!   NAME
!!     make_basis_sm
!!
!!   FUNCTION
!!     returns an S-orthonormal basis v (v^T S v ==1)
!!
!!   NOTES
!!
!!   INPUTS
!!     v and S
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE make_basis_sm(vmatrix,ncol,matrix_s)
    TYPE(cp_fm_type), POINTER                :: vmatrix
    INTEGER, INTENT(IN)                      :: ncol
    TYPE(real_matrix_type), POINTER          :: matrix_s

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


    CALL cp_fm_get_info(matrix=vmatrix,nrow_global=n,ncol_global=ncol_global)
    IF (ncol.gt.ncol_global) CALL stop_program(routineN,"wrong ncol value")

    CALL cp_fm_create(svmatrix,vmatrix%matrix_struct,"SV")
    CALL cp_sm_fm_multiply(matrix_s,vmatrix,svmatrix,ncol)
    
    NULLIFY(fm_struct_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol,&
                             para_env=vmatrix%matrix_struct%para_env, & 
                             context=vmatrix%matrix_struct%context)
    CALL cp_fm_create(overlap_vv,fm_struct_tmp,"overlap_vv")
    CALL cp_fm_struct_release(fm_struct_tmp)

    CALL cp_fm_gemm('T','N',ncol,ncol,n,rone,vmatrix,svmatrix,rzero, overlap_vv)
    CALL cp_fm_cholesky_decompose(overlap_vv)
    CALL cp_fm_triangular_multiply(overlap_vv,vmatrix,n_cols=ncol,side='R',invert_tr=.TRUE.)

    CALL cp_fm_release(overlap_vv)
    CALL cp_fm_release(svmatrix)

    CALL timestop(0.0_dp,handle)

    END SUBROUTINE make_basis_sm

!!****f* qs_mo_methods/make_basis_sv*
!!
!!   NAME
!!     make_basis_sv
!!
!!   FUNCTION
!!     returns an S-orthonormal basis v and the corresponding matrix S*v as well
!!
!!   NOTES
!!
!!   INPUTS
!!     v and S*v
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE make_basis_sv(vmatrix,ncol,svmatrix)

    TYPE(cp_fm_type), POINTER                :: vmatrix
    INTEGER, INTENT(IN)                      :: ncol
    TYPE(cp_fm_type), POINTER                :: svmatrix

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

    CALL cp_fm_get_info(matrix=vmatrix,nrow_global=n,ncol_global=ncol_global)
    IF (ncol.gt.ncol_global) CALL stop_program(routineN,"wrong ncol value")
    
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol,&
                             para_env=vmatrix%matrix_struct%para_env, & 
                             context=vmatrix%matrix_struct%context)
    CALL cp_fm_create(overlap_vv,fm_struct_tmp,"overlap_vv")
    CALL cp_fm_struct_release(fm_struct_tmp)

    CALL cp_fm_gemm('T','N',ncol,ncol,n,rone,vmatrix,svmatrix,rzero, overlap_vv)
    CALL cp_fm_cholesky_decompose(overlap_vv)
    CALL cp_fm_triangular_multiply(overlap_vv,vmatrix,n_cols=ncol,side='R',invert_tr=.TRUE.)
    CALL cp_fm_triangular_multiply(overlap_vv,svmatrix,n_cols=ncol,side='R',invert_tr=.TRUE.) 

    CALL cp_fm_release(overlap_vv)

    CALL timestop(0.0_dp,handle)

    END SUBROUTINE make_basis_sv

!!****f* qs_mo_methods/make_basis_cholesky *
!!
!!   NAME
!!     make_basis_cholesky
!!
!!   FUNCTION
!!     return a set of S orthonormal vectors (C^T S C == 1) where
!!     the cholesky decomposed form of S is passed as an argument
!!
!!   NOTES
!!     if the cholesky decomposed S matrix is not available
!!     use make_basis_sm since this is much faster than computing the
!!     cholesky decomposition of S
!!
!!   INPUTS
!!    - ortho : cholesky decomposed S matrix
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE make_basis_cholesky(vmatrix,ncol,ortho)

    TYPE(cp_fm_type), POINTER                :: vmatrix
    INTEGER, INTENT(IN)                      :: ncol
    TYPE(cp_fm_type), POINTER                :: ortho

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

    CALL cp_fm_get_info(matrix=vmatrix,nrow_global=n,ncol_global=ncol_global)
    IF (ncol.gt.ncol_global) CALL stop_program(routineN,"wrong ncol value")
    
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol,&
                             para_env=vmatrix%matrix_struct%para_env, & 
                             context=vmatrix%matrix_struct%context)
    CALL cp_fm_create(overlap_vv,fm_struct_tmp,"overlap_vv")
    CALL cp_fm_struct_release(fm_struct_tmp)

    CALL cp_fm_triangular_multiply(ortho,vmatrix,n_cols=ncol)
    CALL cp_fm_syrk('U','T',n,rone,vmatrix,1,1,rzero,overlap_vv)
    CALL cp_fm_cholesky_decompose(overlap_vv)
    CALL cp_fm_triangular_multiply(overlap_vv,vmatrix,n_cols=ncol,side='R',invert_tr=.TRUE.)
    CALL cp_fm_triangular_multiply(ortho,vmatrix,n_cols=ncol,invert_tr=.TRUE.)

    CALL cp_fm_release(overlap_vv)

    CALL timestop(0.0_dp,handle)

    END SUBROUTINE make_basis_cholesky

!!****f* qs_mo_methods/make_basis_simple *
!!
!!   NAME
!!     make_basis_simple
!!
!!   FUNCTION
!!     given a set of vectors, return an orthogonal (C^T C == 1) set
!!     spanning the same space (notice, only for cases where S==1)
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE make_basis_simple(vmatrix,ncol)

    TYPE(cp_fm_type), POINTER                :: vmatrix
    INTEGER, INTENT(IN)                      :: ncol

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

    CALL cp_fm_get_info(matrix=vmatrix,nrow_global=n,ncol_global=ncol_global)
    IF (ncol.gt.ncol_global) CALL stop_program(routineN,"wrong ncol value")
    
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol,&
                             para_env=vmatrix%matrix_struct%para_env, & 
                             context=vmatrix%matrix_struct%context)
    CALL cp_fm_create(overlap_vv,fm_struct_tmp,"overlap_vv")
    CALL cp_fm_struct_release(fm_struct_tmp)

    CALL cp_fm_gemm('T','N',ncol,ncol,n,rone,vmatrix,vmatrix,rzero, overlap_vv)
    CALL cp_fm_cholesky_decompose(overlap_vv)
    CALL cp_fm_triangular_multiply(overlap_vv,vmatrix,n_cols=ncol,side='R',invert_tr=.TRUE.)

    CALL cp_fm_release(overlap_vv)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE make_basis_simple

! *****************************************************************************

  SUBROUTINE calculate_density_matrix_blacs(mo_set,density_matrix,error)

!   Purpose: Calculate the density matrix from the MO coeff and the
!            MO occupation numbers.

!   History: - Creation (01.05.2001,MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: density_matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: handle, imo

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_density_matrix_blacs","I","",handle)

    CALL cp_fm_syrk("U","N",mo_set%lfomo-1,mo_set%maxocc,mo_set%mo_coeff,&
                    1,1,0.0_dp,density_matrix)

!   *** Add the density contribution from the fractional occupied MOs ***

    DO imo=mo_set%lfomo,mo_set%homo
      CALL cp_fm_syrk("U","N",1,mo_set%occupation_numbers(imo),&
                      mo_set%mo_coeff,1,imo,1.0_dp,density_matrix)

    END DO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_density_matrix_blacs

! *****************************************************************************

  SUBROUTINE calculate_density_matrix_sparse(mo_set,density_matrix,error)

!   Purpose: Calculate the density matrix

!   History: - Creation (06.2002, Joost VandeVondele)
!            - Fractional occupied orbitals (MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(real_matrix_type), POINTER          :: density_matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: alpha
    TYPE(cp_fm_type), POINTER                :: fm_tmp

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_density_matrix_sparse","I","",handle)

    CALL set_matrix(density_matrix,0.0_dp)

    IF ( .NOT. mo_set%uniform_occupation ) THEN ! not all orbitals 1..homo are equally occupied
       NULLIFY(fm_tmp)
       CALL cp_fm_create(fm_tmp,mo_set%mo_coeff%matrix_struct)
       CALL cp_fm_to_fm(mo_set%mo_coeff,fm_tmp)
       CALL cp_fm_column_scale(fm_tmp,mo_set%occupation_numbers(1:mo_set%homo))
       alpha=1.0_dp

    ELSE
       alpha=mo_set%maxocc
       fm_tmp => mo_set%mo_coeff
       CALL cp_fm_retain(fm_tmp)
    ENDIF

    CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix,&
                            matrix_v=mo_set%mo_coeff,&
                            matrix_g=fm_tmp,&
                            ncol=mo_set%homo,&
                            alpha=alpha)

    CALL cp_fm_release(fm_tmp)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_density_matrix_sparse

!!****f* qs_mo_methods/subspace_eigenvalues_ks *
!!
!!   NAME
!!     subspace_eigenvalues_ks
!!
!!   FUNCTION
!!     computes ritz values of a set of orbitals given a ks_matrix 
!!     rotates the orbitals into eigenstates depending on do_rotation
!!     writes the evals to the screen depending on ionode/scr
!!
!!   NOTES
!!
!!   INPUTS
!!    - orbitals    : S-orthonormal orbitals
!!    - ks_matrix   : Kohn-Sham matrix
!!    - para_env    : of the Kohn-Sham matrix
!!    - evals_arg   : optional, filled with the evals
!!    - ionode/scr  : if present write to unit scr where ionode
!!    - do_rotation : optional rotate orbitals if .TRUE. (default)
!!
!!   MODIFICATION HISTORY
!!     08.2004 documented and added do_rotation [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
! *****************************************************************************

  SUBROUTINE subspace_eigenvalues_ks(orbitals,ks_matrix,evals_arg,para_env,ionode,scr,do_rotation,co_rotate)

!   Purpose: Calculates the eigenvalues (ritz values) for a set of orbitals
!            and rotates them to make them eigenvectors
!            if evals_arg is present its filled
!            if scr and ionode are present they go to the screen as well.

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER                :: orbitals
    TYPE(real_matrix_type), POINTER          :: ks_matrix
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: evals_arg
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN), OPTIONAL            :: ionode
    INTEGER, INTENT(IN), OPTIONAL            :: scr
    LOGICAL, INTENT(IN), OPTIONAL            :: do_rotation
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: co_rotate

    INTEGER                                  :: handle, i, j, ncol_global, &
                                                nrow_global
    LOGICAL                                  :: do_rotation_local
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: e_vectors, h_block, &
                                                weighted_vectors

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_subspace_eigenvalues","I","",handle)

    do_rotation_local=.TRUE.
    IF (PRESENT(do_rotation)) do_rotation_local=do_rotation

    NULLIFY(weighted_vectors,h_block,e_vectors,fm_struct_tmp)
    CALL cp_fm_get_info(matrix=orbitals, &
                                 ncol_global=ncol_global, &
                                 nrow_global=nrow_global)

    IF (ncol_global.GT.0) THEN

        ALLOCATE(evals(ncol_global))

        CALL cp_fm_create(weighted_vectors,orbitals%matrix_struct,"weighted_vectors")
        CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol_global,  ncol_global=ncol_global, &
                                                para_env=orbitals%matrix_struct%para_env, &
                                                context=orbitals%matrix_struct%context)
        CALL cp_fm_create(h_block,fm_struct_tmp, name="h block")
        CALL cp_fm_create(e_vectors,fm_struct_tmp, name="e vectors")
        CALL cp_fm_struct_release(fm_struct_tmp)

        ! h subblock and diag
        CALL cp_sm_fm_multiply(ks_matrix,orbitals,weighted_vectors, ncol_global)


        CALL cp_fm_gemm('T','N',ncol_global,ncol_global,nrow_global,1.0_dp, &
                        orbitals,weighted_vectors,0.0_dp,h_block)

        CALL cp_fm_syevd(h_block,e_vectors,evals)

        ! rotate the orbitals
        IF (do_rotation_local) THEN
            CALL cp_fm_gemm('N','N',nrow_global,ncol_global,ncol_global,1.0_dp, &
                        orbitals,e_vectors,0.0_dp,weighted_vectors)
            CALL cp_fm_to_fm(weighted_vectors,orbitals)
            IF (PRESENT(co_rotate)) THEN
                CALL cp_fm_gemm('N','N',nrow_global,ncol_global,ncol_global,1.0_dp, &
                        co_rotate,e_vectors,0.0_dp,weighted_vectors)
                CALL cp_fm_to_fm(weighted_vectors,co_rotate)
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

        CALL cp_fm_release(weighted_vectors)
        CALL cp_fm_release(h_block)
        CALL cp_fm_release(e_vectors)

        DEALLOCATE(evals)

    ENDIF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE subspace_eigenvalues_ks

! *****************************************************************************
! computes the effective orthonormality of a set of mos given an s-matrix
! orthonormality is the max deviation from unity of the C^T S C
! *****************************************************************************
  SUBROUTINE calculate_orthonormality(orthonormality,mo_array,matrix_s)
    REAL(KIND=dp)                            :: orthonormality
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: matrix_s

    INTEGER                                  :: handle, i, ispin, j, k, n, &
                                                ncol_local, nrow_local, nspin
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: alpha, max_alpha
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: overlap, svec

   NULLIFY(tmp_fm_struct,svec,overlap)

   CALL timeset("calculate_orthonormality","I","",handle)

   nspin=SIZE(mo_array)
   max_alpha=0.0_dp

   DO ispin=1,nspin
      IF (PRESENT(matrix_s)) THEN
        ! get S*C
        CALL cp_fm_create(svec,mo_array(ispin)%mo_set%mo_coeff%matrix_struct)       
        CALL cp_fm_get_info(mo_array(ispin)%mo_set%mo_coeff,&
             nrow_global=n,ncol_global=k)
        CALL cp_sm_fm_multiply(matrix_s,mo_array(ispin)%mo_set%mo_coeff,&
             svec,k)
        ! get C^T (S*C)
        CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=k,ncol_global=k, &
             para_env=mo_array(ispin)%mo_set%mo_coeff%matrix_struct%para_env, &
             context=mo_array(ispin)%mo_set%mo_coeff%matrix_struct%context)
        CALL cp_fm_create(overlap,tmp_fm_struct)       
        CALL cp_fm_struct_release(tmp_fm_struct)
        CALL cp_fm_gemm('T','N',k,k,n,1.0_dp, mo_array(ispin)%mo_set%mo_coeff,&
             svec,0.0_dp,overlap)
        CALL cp_fm_release(svec)
      ELSE
        ! orthogonal basis C^T C
        CALL cp_fm_get_info(mo_array(ispin)%mo_set%mo_coeff,&
             nrow_global=n,ncol_global=k)
        CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=k,ncol_global=k, &
             para_env=mo_array(ispin)%mo_set%mo_coeff%matrix_struct%para_env, &
             context=mo_array(ispin)%mo_set%mo_coeff%matrix_struct%context)
        CALL cp_fm_create(overlap,tmp_fm_struct)       
        CALL cp_fm_struct_release(tmp_fm_struct)
        CALL cp_fm_gemm('T','N',k,k,n,1.0_dp, mo_array(ispin)%mo_set%mo_coeff,&
             mo_array(ispin)%mo_set%mo_coeff,0.0_dp,overlap)
      ENDIF
      CALL cp_fm_get_info(overlap,nrow_local=nrow_local,ncol_local=ncol_local, &
           row_indices=row_indices,col_indices=col_indices)
      DO i=1,nrow_local
         DO j=1,ncol_local
            alpha=overlap%local_data(i,j)
            IF (row_indices(i).eq.col_indices(j)) alpha=alpha-1.0_dp
            max_alpha=MAX(max_alpha,ABS(alpha))
         ENDDO
      ENDDO
      CALL cp_fm_release(overlap)
    ENDDO
    CALL mp_max(max_alpha,mo_array(1)%mo_set%mo_coeff%matrix_struct%para_env%group)
    orthonormality=max_alpha
    ! write(6,*) "max deviation from orthonormalization ",orthonormality

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_orthonormality

! *****************************************************************************
! computes the minimum/maximum magnitudes of C^T C. This could be useful
! to detect problems in the case of nearly singular overlap matrices.
! in this case, we expect the ratio of min/max to be large
! this routine is only similar to mo_orthonormality if S==1
! *****************************************************************************
  SUBROUTINE calculate_magnitude(mo_array,mo_mag_min,mo_mag_max,error)
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    REAL(KIND=dp)                            :: mo_mag_min, mo_mag_max
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_magnitude'

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
             nrow_global=n,ncol_global=k)
      ALLOCATE(evals(k))
      CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=k,ncol_global=k, &
             para_env=mo_array(ispin)%mo_set%mo_coeff%matrix_struct%para_env, &
             context=mo_array(ispin)%mo_set%mo_coeff%matrix_struct%context)
      CALL cp_fm_create(overlap,tmp_fm_struct)
      CALL cp_fm_create(evecs,tmp_fm_struct)
      CALL cp_fm_struct_release(tmp_fm_struct)
      CALL cp_fm_gemm('T','N',k,k,n,1.0_dp, mo_array(ispin)%mo_set%mo_coeff,&
             mo_array(ispin)%mo_set%mo_coeff,0.0_dp,overlap)
      CALL cp_fm_syevd(overlap,evecs,evals,error)
      mo_mag_min=MIN(MINVAL(evals),mo_mag_min)
      mo_mag_max=MAX(MAXVAL(evals),mo_mag_max)
      CALL cp_fm_release(overlap)
      CALL cp_fm_release(evecs)
      DEALLOCATE(evals) 
    ENDDO
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_magnitude

END MODULE qs_mo_methods
