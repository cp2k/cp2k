!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2003 CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
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
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_type,&
                                             cp_to_string,&
                                             cp_warning_level
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

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_mo_types'

  INTERFACE calculate_density_matrix
    MODULE PROCEDURE calculate_density_matrix_sparse,&
                     calculate_density_matrix_blacs
  END INTERFACE

  INTERFACE calculate_subspace_eigenvalues
    MODULE PROCEDURE subspace_eigenvalues_ks ! outdated ,subspace_eigenvalues_vec
  END INTERFACE

  INTERFACE calculate_w_matrix
    MODULE PROCEDURE calculate_w_matrix_1,&
                     calculate_w_matrix_2
  END INTERFACE


  PUBLIC :: make_basis, &
            calculate_density_matrix,&
            calculate_w_matrix,&
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

    INTEGER                                  :: handle, i, info, mypcol, &
                                                myprow, n, ncol_block, npcol, &
                                                nprow, nrow_block
    INTEGER, DIMENSION(9)                    :: desca, desco, descp, descvv, &
                                                descw
    LOGICAL                                  :: found_type, use_cholesky, &
                                                use_pmatrix, use_simple, &
                                                use_sv
    REAL(KIND=dp)                            :: rone, rzero
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, o, p, vv, w
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: overlap_vv

!   ---------------------------------------------------------------------------

    CALL timeset("make_basis","I","",handle)
    NULLIFY(fm_struct_tmp)

! check input ...

    IF (PRESENT(pmatrix)) THEN
       IF (.NOT. PRESENT(work)) THEN
          CALL stop_program("make_basis","pmatrix needs work space !?")
       ENDIF
       use_pmatrix=.TRUE.
    ELSE
       use_pmatrix=.FALSE.
    ENDIF

    IF (PRESENT(ortho)) THEN
       IF (.NOT. PRESENT(otype)) THEN
           CALL stop_program("make_basis","ortho needs type !?")
       ENDIF
       found_type=.FALSE.
       use_sv=.FALSE.
       use_cholesky=.FALSE.
       IF (otype .EQ. "CHOLESKY") THEN
          found_type=.TRUE.
          use_cholesky=.TRUE.
       ENDIF
       IF (otype .EQ. "SV") THEN
          found_type=.TRUE.
          use_sv=.TRUE.
       ENDIF
       IF (.NOT. found_type) CALL stop_program("make_basis","sorry wrong type")
    ELSE
       use_simple=.TRUE.
    ENDIF

    context => vmatrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    rone=1.0_dp
    rzero=0.0_dp
    n   = vmatrix%matrix_struct%nrow_global
    CALL cp_fm_get_info(matrix=vmatrix,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block)
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

    IF (use_sv) THEN
      CALL cp_fm_gemm('T','N',ncol,ncol,n,rone,vmatrix,ortho,rzero, &
                                          overlap_vv)
      CALL cp_fm_cholesky_decompose(overlap_vv)
      a => vmatrix%local_data
      vv => overlap_vv%local_data
      o => ortho%local_data
#if defined(__SCALAPACK)
      desca(:) = vmatrix%matrix_struct%descriptor(:)
      descvv(:) = overlap_vv%matrix_struct%descriptor(:)
      desco(:) = ortho%matrix_struct%descriptor(:)
      CALL pdtrsm('R','U','N','N',n,ncol,rone,vv(1,1),1,1,descvv,a(1,1),1,1,desca)
      CALL pdtrsm('R','U','N','N',n,ncol,rone,vv(1,1),1,1,descvv,o(1,1),1,1,desco)
#else
      CALL dtrsm('R','U','N','N',n,ncol,rone,vv(1,1),ncol,a(1,1),n)
      CALL dtrsm('R','U','N','N',n,ncol,rone,vv(1,1),ncol,o(1,1),n)
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

!   Purpose: Calculate the density matrix from the MO eigenvectors and the
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

    CALL cp_fm_syrk("U","N",mo_set%lfomo-1,mo_set%maxocc,mo_set%eigenvectors,&
                    1,1,0.0_dp,density_matrix)

!   *** Add the density contribution from the fractional occupied MOs ***

    DO imo=mo_set%lfomo,mo_set%homo
      CALL cp_fm_syrk("U","N",1,mo_set%occupation_numbers(imo),&
                      mo_set%eigenvectors,1,imo,1.0_dp,density_matrix)
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

    IF ( mo_set%lfomo .LE. mo_set%homo) THEN ! not all orbitals 1..homo are equally occupied
       NULLIFY(fm_tmp)
       CALL cp_fm_create(fm_tmp,mo_set%eigenvectors%matrix_struct)
       CALL cp_fm_to_fm(mo_set%eigenvectors,fm_tmp)
       CALL cp_fm_column_scale(fm_tmp,mo_set%occupation_numbers(1:mo_set%homo))
       alpha=1.0_dp
    ELSE
       alpha=mo_set%maxocc
       fm_tmp => mo_set%eigenvectors
       CALL cp_fm_retain(fm_tmp)
    ENDIF

    CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix,&
                            matrix_v=mo_set%eigenvectors,&
                            matrix_g=fm_tmp,&
                            ncol=mo_set%homo,&
                            alpha=alpha)

    CALL cp_fm_release(fm_tmp)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_density_matrix_sparse

! *****************************************************************************

  SUBROUTINE calculate_w_matrix_1(mo_set,w_matrix)

!   Purpose: Calculate the W matrix from the MO eigenvectors, MO eigenvalues,
!            and the MO occupation numbers.

!   History: - Creation (03.03.03,MK)
!            - Modification that computes it as a full block, several times (e.g. 20) 
!              faster at the cost of some additional memory

!   ***************************************************************************

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(real_matrix_type), POINTER          :: w_matrix

    CHARACTER(LEN=*), PARAMETER :: routine_name = "calculate_w_matrix_1", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//moduleN//")"

    INTEGER                                  :: handle, imo, istat
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigocc
    TYPE(cp_fm_type), POINTER                :: weighted_vectors

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)
    NULLIFY(weighted_vectors)

    CALL set_matrix(w_matrix,0.0_dp)
    CALL cp_fm_create(weighted_vectors,mo_set%eigenvectors%matrix_struct,"weighted_vectors")
    CALL cp_fm_to_fm(mo_set%eigenvectors,weighted_vectors)

    ! scale every column with the occupation
    ALLOCATE(eigocc(mo_set%homo),stat=istat)
    IF (istat /= 0) CALL stop_memory(routine,"eigocc")
    DO imo=1,mo_set%homo
       eigocc(imo) = mo_set%eigenvalues(imo)*mo_set%occupation_numbers(imo)
    ENDDO
    CALL cp_fm_column_scale(weighted_vectors,eigocc)
    DEALLOCATE(eigocc)

    CALL cp_sm_plus_fm_fm_t(sparse_matrix=w_matrix,&
                            matrix_v=mo_set%eigenvectors,&
                            matrix_g=weighted_vectors,&
                            ncol=mo_set%homo)

    CALL cp_fm_release(weighted_vectors)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_w_matrix_1

! *****************************************************************************

  SUBROUTINE calculate_w_matrix_2(mo_set,ks_matrix,w_matrix)

!   Purpose: Calculate the W matrix from the MO eigenvectors, MO eigenvalues,
!            and the MO occupation numbers.

!   History: - Creation (17.07.2002,MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(real_matrix_type), POINTER          :: ks_matrix, w_matrix

    INTEGER                                  :: handle, ncol_block, &
                                                ncol_global, nrow_block, &
                                                nrow_global
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: h_block, weighted_vectors

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_w_matrix_2","I","",handle)

    NULLIFY(weighted_vectors,h_block,fm_struct_tmp)

    CALL cp_fm_get_info(matrix=mo_set%eigenvectors,&
                        ncol_global=ncol_global,&
                        nrow_global=nrow_global,&
                        nrow_block=nrow_block,&
                        ncol_block=ncol_block)

    CALL cp_fm_create(weighted_vectors,mo_set%eigenvectors%matrix_struct,"weighted_vectors")
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol_global,  ncol_global=ncol_global, &
                                            para_env=mo_set%eigenvectors%matrix_struct%para_env, &
                                            context=mo_set%eigenvectors%matrix_struct%context)
    CALL cp_fm_create(h_block,fm_struct_tmp, name="h block")
    CALL cp_fm_struct_release(fm_struct_tmp)


    IF ( mo_set%lfomo .LE. mo_set%homo) THEN ! not all orbitals 1..homo are equally occupied
      CALL stop_program("calculate_w_matrix_2","NYI")
    END IF

    CALL cp_sm_fm_multiply(ks_matrix,mo_set%eigenvectors,weighted_vectors, &
                           mo_set%homo,mo_set%eigenvectors%matrix_struct%para_env)
    CALL cp_fm_gemm('T','N',ncol_global,ncol_global,nrow_global,1.0_dp, &
                    mo_set%eigenvectors,weighted_vectors,0.0_dp,h_block)
    CALL cp_fm_gemm('N','N',nrow_global,ncol_global,ncol_global,1.0_dp, &
                    mo_set%eigenvectors,h_block,0.0_dp,weighted_vectors)

    ! scale every column with the occupation
    CALL cp_fm_column_scale(weighted_vectors,mo_set%occupation_numbers)

    CALL set_matrix(w_matrix,0.0_dp)
    CALL cp_sm_plus_fm_fm_t(sparse_matrix=w_matrix,&
                            matrix_v=mo_set%eigenvectors,&
                            matrix_g=weighted_vectors,&
                            ncol=mo_set%homo)

    CALL cp_fm_release(weighted_vectors)
    CALL cp_fm_release(h_block)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_w_matrix_2

! *****************************************************************************

  SUBROUTINE subspace_eigenvalues_ks(orbitals,ks_matrix,evals_arg,para_env,ctxt,ionode,scr)

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
    TYPE(cp_blacs_env_type), POINTER         :: ctxt
    LOGICAL, OPTIONAL                        :: ionode
    INTEGER, OPTIONAL                        :: scr

    INTEGER                                  :: handle, i, j, ncol_block, &
                                                ncol_global, nrow_block, &
                                                nrow_global
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: evals
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: e_vectors, h_block, &
                                                weighted_vectors

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_subspace_eigenvalues","I","",handle)

    NULLIFY(weighted_vectors,h_block,e_vectors,fm_struct_tmp)
    CALL cp_fm_get_info(matrix=orbitals, &
                                 ncol_global=ncol_global, &
                                 nrow_global=nrow_global, &
                                 nrow_block=nrow_block, &
                                 ncol_block=ncol_block)

    ALLOCATE(evals(ncol_global))

    CALL cp_fm_create(weighted_vectors,orbitals%matrix_struct,"weighted_vectors")
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol_global,  ncol_global=ncol_global, &
                                            para_env=orbitals%matrix_struct%para_env, &
                                            context=orbitals%matrix_struct%context)
    CALL cp_fm_create(h_block,fm_struct_tmp, name="h block")
    CALL cp_fm_create(e_vectors,fm_struct_tmp, name="e vectors")
    CALL cp_fm_struct_release(fm_struct_tmp)

    ! h subblock and diag
    CALL cp_sm_fm_multiply(ks_matrix,orbitals,weighted_vectors, &
                                       ncol_global,para_env=para_env)


    CALL cp_fm_gemm('T','N',ncol_global,ncol_global,nrow_global,1.0_dp, &
                    orbitals,weighted_vectors,0.0_dp,h_block)

    CALL cp_fm_syevd(h_block,e_vectors,evals)

    ! rotate the orbitals
    CALL cp_fm_gemm('N','N',nrow_global,ncol_global,ncol_global,1.0_dp, &
                    orbitals,e_vectors,0.0_dp,weighted_vectors)
    CALL cp_fm_to_fm(weighted_vectors,orbitals)

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
             WRITE(scr,'(1X,4F12.8)') evals(i:i+j)
            CASE(2)
             WRITE(scr,'(1X,3F12.8)') evals(i:i+j)
            CASE(1)
             WRITE(scr,'(1X,2F12.8)') evals(i:i+j)
            CASE(0)
             WRITE(scr,'(1X,1F12.8)') evals(i:i+j)
            END SELECT
         ENDDO
      ENDIF
    ENDIF

    CALL cp_fm_release(weighted_vectors)
    CALL cp_fm_release(h_block)
    CALL cp_fm_release(e_vectors)

    DEALLOCATE(evals)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE subspace_eigenvalues_ks

! *****************************************************************************
! computes the effective orthonormality of a set of mos given an s-matrix
! orthonormality is the max deviation from unity of the C^T S C
! *****************************************************************************
  SUBROUTINE calculate_orthonormality(mo_array,matrix_s,orthonormality)
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(real_matrix_type), POINTER          :: matrix_s
    REAL(KIND=dp)                            :: orthonormality

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
          ! get S*C
          CALL cp_fm_create(svec,mo_array(ispin)%mo_set%eigenvectors%matrix_struct)       
          CALL cp_fm_get_info(mo_array(ispin)%mo_set%eigenvectors,nrow_global=n,ncol_global=k)
          CALL cp_sm_fm_multiply(matrix_s,mo_array(ispin)%mo_set%eigenvectors, svec, &
                                      k, para_env=mo_array(ispin)%mo_set%eigenvectors%matrix_struct%para_env)
          ! get C^T (S*C)
          CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=k,ncol_global=k, &
                                   para_env=mo_array(ispin)%mo_set%eigenvectors%matrix_struct%para_env, &
                                   context=mo_array(ispin)%mo_set%eigenvectors%matrix_struct%context)
          CALL cp_fm_create(overlap,tmp_fm_struct)       
          CALL cp_fm_struct_release(tmp_fm_struct)
          CALL cp_fm_gemm('T','N',k,k,n,1.0_dp, mo_array(ispin)%mo_set%eigenvectors,svec,0.0_dp,overlap)
          CALL cp_fm_get_info(overlap,nrow_local=nrow_local,ncol_local=ncol_local, &
                              row_indices=row_indices,col_indices=col_indices)
          DO i=1,nrow_local
           DO j=1,ncol_local
              alpha=overlap%local_data(i,j)
              IF (row_indices(i).eq.col_indices(j)) alpha=alpha-1.0_dp
              max_alpha=MAX(max_alpha,ABS(alpha))
           ENDDO
          ENDDO
          CALL cp_fm_release(svec)
          CALL cp_fm_release(overlap)
    ENDDO
    CALL mp_max(max_alpha,mo_array(1)%mo_set%eigenvectors%matrix_struct%para_env%group)
    orthonormality=max_alpha
    ! write(6,*) "max deviation from orthonormalization ",orthonormality

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_orthonormality

END MODULE qs_mo_methods
