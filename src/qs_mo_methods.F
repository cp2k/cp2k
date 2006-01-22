!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2003 CP2K developers group
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
                                             cp_fm_syrk
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_init_random,&
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
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_mo_types'

  INTERFACE calculate_density_matrix
    MODULE PROCEDURE calculate_density_matrix_sparse,&
                     calculate_density_matrix_blacs
  END INTERFACE

  INTERFACE calculate_subspace_eigenvalues
    MODULE PROCEDURE subspace_eigenvalues_ks 
  END INTERFACE

  PUBLIC :: make_basis, &
            calculate_density_matrix,&
            calculate_subspace_eigenvalues, &
            calculate_orthonormality

CONTAINS


! XXXXXXXXXXXXXXXXXXXX SHOULD USE cp_fm_ instead of home-cooked scalapack

! *****************************************************************************
! orthogonalises a set of vectors, taking into account S
!
! c) if ortho is present, this is
!            1) the decomposed 'S' matrix of the system
!               so that v^T*S*v=1
!            2) S*V (no pmatrix in this case) both v, and sv will be changed
!               v->v*u^-1 sv*u^-1
!
! a) if pmatrix is not present, we find an orthogonal basis for the space the input spans
!    using S as a metric (that is we ortho U*V with U^T*U=S)
!    (a random initial basis can be obtained by cp_fm_init_random)
!
! b) if pmatrix is present, it is assumed that this is a projector of rank nmo
! and a basis for this projector is found i.e. pmatrix=v*v^T
!
!
! *****************************************************************************
  SUBROUTINE make_basis(vmatrix,ncol,ortho,otype,pmatrix,work)
! *****************************************************************************

    TYPE(cp_fm_type), POINTER                :: vmatrix
    INTEGER, INTENT(IN)                      :: ncol
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: ortho
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: otype
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: pmatrix, work

    INTEGER                                  :: handle, i, info, n
    INTEGER, DIMENSION(9)                    :: desca, desco, descp, descvv, &
                                                descw
    LOGICAL                                  :: use_cholesky, &
                                                use_pmatrix, use_simple, &
                                                use_sv
    REAL(KIND=dp), PARAMETER                 :: rone=1.0_dp, rzero=0.0_dp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, o, p, vv, w
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: overlap_vv

!   ---------------------------------------------------------------------------
    IF (ncol.EQ.0) return

    CALL timeset("make_basis","I","",handle)
    NULLIFY(fm_struct_tmp)

! check input ...

    use_pmatrix=.FALSE.
    use_sv=.FALSE.
    use_cholesky=.FALSE.
    use_simple=.FALSE.

    IF (PRESENT(pmatrix)) THEN
       IF (.NOT. PRESENT(work)) THEN
          CALL stop_program("make_basis","pmatrix needs work space !?")
       ENDIF
       use_pmatrix=.TRUE.
    ELSE
      IF (PRESENT(ortho)) THEN ! ortho can be two different things here (see comment above)
        IF (.NOT. PRESENT(otype)) THEN
            CALL stop_program("make_basis","ortho needs type !?")
        ENDIF
        SELECT CASE(otype)
        CASE ("CHOLESKY")
           use_cholesky=.TRUE.
        CASE ("SV")
           use_sv=.TRUE.
        CASE DEFAULT
           CALL stop_program("make_basis","sorry wrong type")
        END SELECT 
      ELSE
        use_simple=.TRUE.
      ENDIF
    ENDIF

    context => vmatrix%matrix_struct%context
    CALL cp_fm_get_info(matrix=vmatrix,nrow_global=n)
    IF (ncol.lt.0 .OR. ncol.gt.vmatrix%matrix_struct%ncol_global) THEN
        CALL stop_program("make_basis","wrong ncol value")
    ENDIF
    
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol,&
                             para_env=vmatrix%matrix_struct%para_env, context=context)
    CALL cp_fm_create(overlap_vv,fm_struct_tmp,"overlap_vv")
    CALL cp_fm_struct_release(fm_struct_tmp)

    a => vmatrix%local_data
    vv => overlap_vv%local_data
#if defined(__SCALAPACK)
    desca(:) = vmatrix%matrix_struct%descriptor(:)
    descvv(:) = overlap_vv%matrix_struct%descriptor(:)
#endif
    IF (use_pmatrix) THEN
       CALL cp_fm_init_random(work,ncol)
       p => pmatrix%local_data
       w => work%local_data
#if defined(__SCALAPACK)
       descp(:) = pmatrix%matrix_struct%descriptor(:)
       descw(:) = work%matrix_struct%descriptor(:)
       CALL PDSYMM('L','U',n,ncol,rone,p(1,1),1,1,descp,w(1,1),1,1,descw,rzero,a(1,1),1,1,desca)
       IF (use_cholesky) THEN
         o => ortho%local_data
         desco(:) = ortho%matrix_struct%descriptor(:)
         CALL PDTRMM('L','U','N','N',n,ncol,rone,o(1,1),1,1,desco,a(1,1),1,1,desca)
       ENDIF
#else
       ! U*P*U^T*V
       CALL DSYMM('L','U',n,ncol,rone,p(1,1),n,w(1,1),n,rzero,a(1,1),n)
       IF (use_cholesky) THEN
         o => ortho%local_data
         CALL DTRMM('L','U','N','N',n,ncol,rone,o(1,1),n,a(1,1),n)
       ENDIF
#endif
    ELSE
       ! U*V should be used
       IF (use_cholesky) THEN
         o => ortho%local_data
#if defined(__SCALAPACK)
         desco(:) = ortho%matrix_struct%descriptor(:)
         CALL PDTRMM('L','U','N','N',n,ncol,rone,o(1,1),1,1,desco,a(1,1),1,1,desca)
#else
         CALL DTRMM('L','U','N','N',n,ncol,rone,o(1,1),n,a(1,1),n)
#endif
       ENDIF
    ENDIF

    IF (use_sv .OR. use_simple) THEN
      IF (use_simple) THEN
          CALL cp_fm_gemm('T','N',ncol,ncol,n,rone,vmatrix,vmatrix,rzero, &
                                              overlap_vv)
      ELSE
          CALL cp_fm_gemm('T','N',ncol,ncol,n,rone,vmatrix,ortho,rzero, &
                                              overlap_vv)
      ENDIF
      CALL cp_fm_cholesky_decompose(overlap_vv)
      a => vmatrix%local_data
      vv => overlap_vv%local_data
#if defined(__SCALAPACK)
      desca(:) = vmatrix%matrix_struct%descriptor(:)
      descvv(:) = overlap_vv%matrix_struct%descriptor(:)
      CALL pdtrsm('R','U','N','N',n,ncol,rone,vv(1,1),1,1,descvv,a(1,1),1,1,desca)
      IF (.NOT. use_simple) THEN
          o => ortho%local_data
          desco(:) = ortho%matrix_struct%descriptor(:)
          CALL pdtrsm('R','U','N','N',n,ncol,rone,vv(1,1),1,1,descvv,o(1,1),1,1,desco)
      ENDIF
#else
      CALL dtrsm('R','U','N','N',n,ncol,rone,vv(1,1),ncol,a(1,1),n)
      IF (.NOT. use_simple) THEN
          o => ortho%local_data
          CALL dtrsm('R','U','N','N',n,ncol,rone,vv(1,1),ncol,o(1,1),n)
      ENDIF
#endif
    ELSE
      ! ortho (U*V) inner product should be one
      a => vmatrix%local_data
      vv => overlap_vv%local_data
#if defined(__SCALAPACK)
      desca(:) = vmatrix%matrix_struct%descriptor(:)
      descvv(:) = overlap_vv%matrix_struct%descriptor(:)
      CALL PDSYRK('U','T',ncol,n,rone,a(1,1),1,1,desca,rzero,vv(1,1),1,1,descvv,ncol)
      CALL PDPOTRF('U',ncol,vv(1,1),1,1,descvv,info)
      IF (info.NE.0) CALL stop_program("make_basis","DPOTRF degenerate vectors")
      CALL PDTRSM('R','U','N','N',n,ncol,rone,vv(1,1),1,1,descvv,a(1,1),1,1,desca)
#else
      CALL DSYRK('U','T',ncol,n,rone,a(1,1),n,rzero,vv(1,1),ncol)
      CALL DPOTRF('U',ncol,vv(1,1),ncol,info)
      IF (info.NE.0) CALL stop_program("make_basis","DPOTRF degenerate vectors")
      CALL DTRSM('R','U','N','N',n,ncol,rone,vv(1,1),ncol,a(1,1),n)
#endif
    ENDIF

    IF (use_cholesky) THEN
#if defined(__SCALAPACK)
      CALL PDTRSM('L','U','N','N',n,ncol,rone,o(1,1),1,1,desco,a(1,1),1,1,desca)
#else
      CALL DTRSM('L','U','N','N',n,ncol,rone,o(1,1),n,a(1,1),n)
#endif
    ENDIF

    CALL cp_fm_release(overlap_vv)

    CALL timestop(0.0_dp,handle)

    END SUBROUTINE make_basis

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
    LOGICAL, OPTIONAL, INTENT(IN)            :: ionode
    INTEGER, OPTIONAL, INTENT(IN)            :: scr
    LOGICAL, OPTIONAL, INTENT(IN)            :: do_rotation
    TYPE(cp_fm_type), POINTER, OPTIONAL      :: co_rotate

    INTEGER                                  :: handle, i, j, &
                                                ncol_global,  &
                                                nrow_global
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: e_vectors, h_block, &
                                                weighted_vectors
    LOGICAL :: do_rotation_local

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
    TYPE(real_matrix_type),POINTER, OPTIONAL :: matrix_s

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

END MODULE qs_mo_methods
