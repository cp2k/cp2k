!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/qs_blacs [1.0] *
!!
!!   NAME
!!     qs_blacs
!!
!!   FUNCTION
!!     BLACS
!!
!!   AUTHOR
!!     Matthias Krack (22.05.2001)
!!
!!   MODIFICATION HISTORY
!!     1) Changed the size work array for syevx (Joost VandeVondele 28.02.02)
!!     1) Introduced some cholesky guys         (Joost VandeVondele    05.02)
!!     2) Introduced a symmetric matrix x distributed vector multiply ( JVdV)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_blacs
! *****************************************************************************
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_from_globenv,&
                                             cp_blacs_env_release,&
                                             cp_blacs_env_retain,&
                                             finish_blacs,&
                                             start_blacs
  USE cp_error_handling,               ONLY: cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_unimplemented_error
  USE cp_fm_basic_linalg,              ONLY: cp_fm_add,&
                                             cp_fm_gemm,&
                                             cp_fm_schur_product,&
                                             cp_fm_symm,&
                                             cp_fm_syrk,&
                                             cp_fm_trace,&
                                             cp_fm_triangular_multiply,&
                                             cp_sm_fm_multiply,&
                                             cp_sm_plus_fm_fm_t
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_struct,                    ONLY: cp_f_matrix_struct_type,&
                                             cp_fmstruct_create,&
                                             cp_fmstruct_equivalent,&
                                             cp_fmstruct_release,&
                                             cp_fmstruct_retain,&
                                             optimal_blacs_col_block_size,&
                                             optimal_blacs_row_block_size
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_create2,&
                                             cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_init_random,&
                                             cp_fm_maxval,&
                                             cp_fm_release,&
                                             cp_fm_retain,&
                                             cp_fm_set_all,&
                                             cp_fm_set_element,&
                                             cp_fm_set_submatrix,&
                                             cp_full_matrix_p_type,&
                                             cp_full_matrix_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_create2,&
                                             cp_fm_vect_dealloc
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_log,&
                                             cp_logger_get_unit_nr,&
                                             cp_logger_type,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_matrix_utils,                 ONLY: cp_dealloc_ref,&
                                             cp_get,&
                                             cp_init,&
                                             cp_matrix_block_iterator,&
                                             cp_next
  USE cp_para_env,                     ONLY: cp_para_env_from_globenv,&
                                             cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             wp => dp,&
                                             wp_size => dp_size
  USE machine,                         ONLY: m_cputime
  USE mathlib,                         ONLY: symmetrize_matrix
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_gather,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sync
  USE parallel,                        ONLY: cp2k_is_parallel
  USE sparse_matrix_types,             ONLY: first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_type,&
                                             sparse_plus_rep_rept,&
                                             sparse_times_replicated
  USE string_utilities,                ONLY: compress
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN='qs_blacs'

! *** Public parameters ***
  PUBLIC :: optimal_blacs_col_block_size,&
       optimal_blacs_row_block_size

! *** Public data types ***

  PUBLIC :: cp_full_matrix_p_type,&
            cp_full_matrix_type

! *** Public subroutines ***
! at most 39 continuation lines :-)
  PUBLIC :: cp_fm_create2,&
       cp_fm_release,cp_fm_vect_dealloc,&
       cp_fm_init_random,&
       copy_blacs_to_blacs_matrix, copy_blacs_to_full_matrix,&
       copy_blacs_to_sparse_matrix, copy_sparse_to_blacs_matrix,&
       blacs_replicated_copy,&
       blacs_column_copy,&
       cp_fm_get_info,&
       write_blacs_matrix, &
       cp_fm_get_submatrix, &
       cp_fm_set_submatrix, blacs_scale_and_d,&
       copy_vector_to_blacs, &
       copy_blacs_to_vector, blacs_find_vector_distribution,&
       read_blacs_matrix, replicate_blacs_matrix,&
       blacs_set_local_block, cp_fm_vect_create2,&
       blacs_daxpy,  blacs_get_element
  
  PUBLIC :: start_blacs,&
       finish_blacs, get_blacs_info
       
  PUBLIC ::  cp_fm_add,&
            blacs_add_to_element, &
            cp_fm_gemm,&
            blacs_vector_symv, &
            cp_fm_maxval, cp_fm_set_all,&
            cp_fm_set_element, blacs_syevd,&
            blacs_syevx, cp_fm_symm,&
            cp_fm_syrk, cp_fm_trace,&
            power_blacs_matrix,&
            symmetrise_blacs_matrix,&
            add_blacs_to_block_diag_sm,  cp_fm_cholesky_decompose, &
            cp_fm_cholesky_invert,  cp_fm_cholesky_reduce, &
            cp_fm_cholesky_restore,  blacs_jacobi_davidson, &
            blacs_make_basis, &
            cp_fm_schur_product, &
            blacs_transpose,  blacs_scale_matrix, &
            cp_sm_fm_multiply, &
            cp_sm_plus_fm_fm_t, &
            cp_fm_triangular_multiply

  PUBLIC :: blacs_jacobi,&
            blacs_block_jacobi_scaled,&
            blacs_block_jacobi_classic

! *****************************************************************************

CONTAINS

! *****************************************************************************



  SUBROUTINE blacs_set_local_block(matrix,BLOCK)
! *****************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix
    REAL(wp),DIMENSION(:,:), POINTER          :: BLOCK

    REAL(wp),DIMENSION(:,:), POINTER          :: mblock
    INTEGER :: nelements

    mblock => matrix%local_data
    nelements = SIZE(mblock,1)*&
         matrix%matrix_struct%ncol_local(matrix%matrix_struct%context%mepos(2))

    CALL DCOPY(nelements,BLOCK(1,1),1,mblock(1,1),1)

  END SUBROUTINE blacs_set_local_block


! *****************************************************************************

  SUBROUTINE blacs_daxpy(matrix_a,c1,alpha,c2,error)
    TYPE(cp_full_matrix_type), POINTER       :: matrix_a
    INTEGER, INTENT(IN)                    :: c1,c2
    REAL(wp), INTENT(IN)                   :: alpha
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    INTEGER :: myprow,mypcol,n
    REAL(wp), DIMENSION(:,:), POINTER :: a
    INTEGER, DIMENSION(9) :: desca

    myprow=matrix_a%matrix_struct%context%mepos(1)
    mypcol=matrix_a%matrix_struct%context%mepos(2)
    a => matrix_a%local_data
    n = SIZE(a,1)

#if defined(__SCALAPACK)
    desca(:)=matrix_a%matrix_struct%descriptor(:)
    CALL PDAXPY(n,alpha,a(1,1),1,c1,desca,1,a(1,1),1,c2,desca,1)
#else
    CALL DAXPY(n,alpha,a(1,c1),1,a(1,c2),1)
#endif

  END SUBROUTINE blacs_daxpy

! *****************************************************************************

  SUBROUTINE blacs_add_to_element(matrix,irow_global,icol_global,alpha)

!   Purpose: Add alpha to the BLACS matrix element
!            matrix(irow_global,icol_global).

!   History: - Creation (08.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix
    REAL(wp), INTENT(IN)                      :: alpha
    INTEGER, INTENT(IN)                       :: icol_global,&
                                                 irow_global

!   *** Local variables ***

    INTEGER :: icol_local,ipcol,iprow,irow_local,mypcol,myprow,npcol,&
               nprow

    REAL(wp), DIMENSION(:,:), POINTER :: a

!   ---------------------------------------------------------------------------

    myprow=matrix%matrix_struct%context%mepos(1)
    mypcol=matrix%matrix_struct%context%mepos(2)
    nprow=matrix%matrix_struct%context%num_pe(1)
    npcol=matrix%matrix_struct%context%num_pe(2)

    a => matrix%local_data

#if defined(__SCALAPACK)
    CALL infog2l(irow_global,icol_global,matrix%matrix_struct%descriptor,&
         nprow,npcol,myprow,mypcol,&
         irow_local,icol_local,iprow,ipcol)

    IF ((iprow == myprow).AND.(ipcol == mypcol)) THEN
       a(irow_local,icol_local) = a(irow_local,icol_local) + alpha
    END IF
#else
    a(irow_global,icol_global) = a(irow_global,icol_global) + alpha
#endif

  END SUBROUTINE blacs_add_to_element
! *****************************************************************************

  SUBROUTINE blacs_scale_matrix(matrixa,scaling,kmax)
! ugly routine to perform a scaling of every column in a matrix
! ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrixa
    REAL(wp), DIMENSION(:), POINTER           :: scaling
    INTEGER, OPTIONAL :: kmax

!   *** Local variables ***

    INTEGER :: icol_global,icol_local,ipcol,iprow,irow_local,mypcol,&
               myprow,npcol,nprow,source,i,k,n

    REAL(wp), DIMENSION(:,:), POINTER :: a

!   ---------------------------------------------------------------------------

    myprow=matrixa%matrix_struct%context%mepos(1)
    mypcol=matrixa%matrix_struct%context%mepos(2)
    nprow=matrixa%matrix_struct%context%num_pe(1)
    npcol=matrixa%matrix_struct%context%num_pe(2)

    a => matrixa%local_data
    n =  SIZE(a,1)
    k =  matrixa%matrix_struct%ncol_global

#if defined(__SCALAPACK)
    IF (PRESENT(kmax)) k =  kmax

    DO icol_global=1,k
       CALL infog2l(1,icol_global,matrixa%matrix_struct%descriptor,&
            nprow,npcol,myprow,mypcol,&
            irow_local,icol_local,iprow,ipcol)
       IF ((ipcol == mypcol)) THEN
           CALL DSCAL(n,scaling(icol_global),a(1,icol_local),1)
       END IF
    ENDDO
#else
    DO i=1,k
      CALL DSCAL(n,scaling(i),a(1,i),1)
    ENDDO
#endif
  END SUBROUTINE blacs_scale_matrix

!
! decides what fraction of the vector has to be on the local pe,
! giving each pe about the same amount of data (that is nloc elements),
! in this way we can map more easily to the blacs distribution
! the size of the vector is determined from matrix_a
! columns of matrix_a should represent a vector
! notice that it is the blacs layout that defines how the vector is distributed
! and what element of the vector actually is where

  SUBROUTINE blacs_find_vector_distribution(matrix_a,nloc)
    TYPE(cp_full_matrix_type), POINTER       :: matrix_a
    INTEGER, INTENT(OUT)                   :: nloc

    TYPE(cp_blacs_env_type), POINTER       :: context
    INTEGER :: n,nloc_blacs

    n = matrix_a%matrix_struct%nrow_global
    context => matrix_a%matrix_struct%context

!!FM#if defined(__SCALAPACK)
    nloc_blacs = matrix_a%matrix_struct%nrow_local(context%mepos(1))
    nloc = nloc_blacs / context%num_pe(2)
    IF (context%mepos(2).EQ.0) THEN
         nloc = nloc_blacs - (context%num_pe(2) - 1) * nloc
         ! bug if nloc on proc 0 has a different value than on other proc?
    ENDIF
!!FM#else
!!FM    nloc   = n
!!FM#endif

  END SUBROUTINE blacs_find_vector_distribution
!
! copy a blacs vector to a distributed vector
! has obviously to be consistent with find_vector_distribution
! see notes at blacs_vector_symv
! use_col_nr specifies the column of the matrix_a to use.
! this is the global col index
!
  SUBROUTINE copy_blacs_to_vector(matrix_a,use_col_nr,v)
    TYPE(cp_full_matrix_type), POINTER       :: matrix_a
    INTEGER, INTENT(IN)                    :: use_col_nr
    REAL(wp), POINTER, DIMENSION(:)        :: v
    INTEGER :: nprow,npcol,myprow,mypcol,n,nloc_blacs,source,nloc,nloc0,i
    INTEGER :: irow_global,icol_global,irow_local,icol_local,iprow,ipcol
    INTEGER :: ndata,offset
    INTEGER, DIMENSION(9) :: desca
    TYPE(cp_blacs_env_type), POINTER       :: context

    n = matrix_a%matrix_struct%nrow_global
    context => matrix_a%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

#if defined (__SCALAPACK)
    desca(:) = matrix_a%matrix_struct%descriptor(:)
    irow_global=1
    icol_global=use_col_nr
    CALL infog2l(irow_global,icol_global,desca,nprow,npcol,myprow,mypcol,&
         irow_local,icol_local,iprow,ipcol)

    nloc_blacs = matrix_a%matrix_struct%nrow_local(context%mepos(1))
    nloc  = nloc_blacs / npcol                  ! all cols
    nloc0 = nloc_blacs - (npcol - 1) * nloc     ! col zero

    DO i=0,npcol-1

       ndata=nloc
       IF (i.EQ.0) ndata=nloc0

       IF (mypcol.EQ.ipcol) THEN ! get your part of the data
          offset=nloc0+1+(i-1)*nloc
          IF (i.EQ.0) offset=1

          IF (mypcol.EQ.i) THEN
             CALL DCOPY(ndata,matrix_a%local_data(offset,icol_local),1,&
                        v,1)
          ELSE
             CALL DGESD2D(context%group,ndata,1,   &
             matrix_a%local_data(offset,icol_local),  &
             ndata,myprow,i)
          ENDIF
       ELSE
          IF (mypcol.EQ.i) THEN
             CALL DGERV2D(context%group,ndata,1,   &
             v, &
             ndata,myprow,ipcol)
          ENDIF
       ENDIF
    ENDDO
#else
    nloc = n
    IF (matrix_a%matrix_struct%nrow_local(myprow).LT.use_col_nr) THEN
       CALL stop_program("copy_blacs_to_vector","hmm")
    ENDIF
    CALL DCOPY(nloc,matrix_a%local_data(1,use_col_nr),1,v,1)
#endif

  END SUBROUTINE copy_blacs_to_vector

!
! copy a full distributed vector into a blacs vector
! has obviously to be consistent with find_vector_distribution
! see notes at blacs_vector_symv
!
  SUBROUTINE copy_vector_to_blacs(v,matrix_a,use_col_nr)
    TYPE(cp_full_matrix_type), POINTER       :: matrix_a
    TYPE(global_environment_type)          :: globenv
    INTEGER, INTENT(IN)                    :: use_col_nr
    REAL(wp), POINTER, DIMENSION(:)        :: v
    INTEGER :: nprow,npcol,myprow,mypcol,n,nloc_blacs,source,nloc,nloc0,i
    INTEGER :: irow_global,icol_global,irow_local,icol_local,iprow,ipcol
    INTEGER :: ndata,offset
    INTEGER, DIMENSION(9) :: desca
    TYPE(cp_blacs_env_type), POINTER       :: context

    n = matrix_a%matrix_struct%nrow_global
    context => matrix_a%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    
#if defined (__SCALAPACK)
    desca(:) = matrix_a%matrix_struct%descriptor(:)
    irow_global=1
    icol_global=use_col_nr
    CALL infog2l(irow_global,icol_global,desca,nprow,npcol,myprow,mypcol,&
         irow_local,icol_local,iprow,ipcol)

    nloc_blacs = matrix_a%matrix_struct%nrow_local(myprow)
    nloc  = nloc_blacs / npcol                ! all cols
    nloc0 = nloc_blacs - (npcol - 1) * nloc   ! col zero

    DO i=0,npcol-1

       ndata=nloc
       IF (i.EQ.0) ndata=nloc0

       IF (mypcol.EQ.ipcol) THEN ! get your part of the data
          offset=nloc0+1+(i-1)*nloc
          IF (i.EQ.0) offset=1

          IF (mypcol.EQ.i) THEN
             CALL DCOPY(ndata,v,1,matrix_a%local_data(offset,icol_local),1)
          ELSE
             CALL DGERV2D(context%group,ndata,1,   &
             matrix_a%local_data(offset,icol_local),  &
             ndata,myprow,i)
          ENDIF
       ELSE
          IF (mypcol.EQ.i) THEN
             CALL DGESD2D(context%group,ndata,1,   &
             v, &
             ndata,myprow,ipcol)
          ENDIF
       ENDIF
    ENDDO
#else
    nloc = n
    IF (matrix_a%matrix_struct%ncol_local(mypcol).LT.use_col_nr) THEN
       CALL stop_program("copy_vector_to_blacs","hmm")
    ENDIF
    CALL DCOPY(nloc,v,1,matrix_a%local_data(1,use_col_nr),1)
#endif

  END SUBROUTINE copy_vector_to_blacs

! *****************************************************************************
! computes the product of an upper symmetric blacs matrix with vin
! and stores the result in vout ! vin (vout) is a distributed vector,
! (so that every pe has a part of vin, not necessarily the blacs way thus)
! it is not evident where a given vector element v_i will be
! (i.e. pe and v(xxx) ) since this is determined by the blacs ordering
! of matrix elements
!
! this routine needs a matrix_buf, i.e. a blacs matrix with at least two columns
! that will be overwritten by the results (i.e. first column vin, second column
! vout )
! *****************************************************************************

  SUBROUTINE blacs_vector_symv(matrix_a,vin,vout,matrix_buf)
! *****************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix_a,matrix_buf
    REAL(wp), DIMENSION(:), POINTER           :: vin,vout

!   *** Local variables ***

    INTEGER :: handle,n,inc
    REAL(wp), DIMENSION(:,:), POINTER :: a,b
    REAL(wp) :: alpha,beta
    INTEGER :: nprow,npcol,myprow,mypcol
    INTEGER, DIMENSION(9) :: desca,descb
    TYPE(cp_blacs_env_type), POINTER :: context

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_vector_symv","I","",handle)

    context => matrix_a%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    alpha=1.0_wp
    beta=0.0_wp
    inc=1

    a => matrix_a%local_data
    n = matrix_a%matrix_struct%nrow_global
#if defined(__SCALAPACK)
    desca(:) = matrix_a%matrix_struct%descriptor(:)
    b => matrix_buf%local_data
    descb(:) = matrix_buf%matrix_struct%descriptor(:)
    CALL copy_vector_to_blacs(vin,matrix_buf,1)
    CALL pdsymv('U',n,alpha,a(1,1),1,1,desca,b(1,1),1,1,descb,inc,beta,b(1,1),1,2,descb,inc)
    CALL copy_blacs_to_vector(matrix_buf,2,vout)
#else
    CALL dsymv('U',n,alpha,a(1,1),n,vin(1),inc,beta,vout(1),inc)
#endif

    CALL timestop(0.0_wp,handle)

    END SUBROUTINE blacs_vector_symv

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
  SUBROUTINE blacs_make_basis(vmatrix,ncol,ortho,otype,pmatrix,work)
! *****************************************************************************

    TYPE(cp_full_matrix_type), POINTER                   :: vmatrix
    TYPE(cp_full_matrix_type), POINTER, OPTIONAL         :: pmatrix,work,ortho
    INTEGER, INTENT(IN)                                :: ncol
    CHARACTER ( LEN = * ), INTENT ( IN ),OPTIONAL      :: otype

!   *** Local variables ***

    INTEGER :: handle,n,info
    REAL(wp), DIMENSION(:,:), POINTER :: a,b,vv,w,p,o
    REAL(wp) :: rone,rzero
    INTEGER :: nprow,npcol,myprow,mypcol,nrow_block,ncol_block,i
    INTEGER, DIMENSION(9) :: desca,descb,descvv,descw,descp,desco
    TYPE(cp_full_matrix_type), POINTER           :: overlap_vv
    LOGICAL :: found_type,use_pmatrix,use_cholesky,use_sv,use_simple
    TYPE(cp_blacs_env_type), POINTER :: context


!   ---------------------------------------------------------------------------

    CALL timeset("blacs_make_basis","I","",handle)

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
    rone=1.0_wp
    rzero=0.0_wp
    n   = vmatrix%matrix_struct%nrow_global
    CALL cp_fm_get_info(matrix=vmatrix,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block)

    CALL cp_fm_create2(new_matrix=overlap_vv,&
                               nrow_global=ncol,&
                               ncol_global=ncol,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="overlap_vv",&
                               context=context,&
                               para_env=vmatrix%matrix_struct%para_env)

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

    CALL timestop(0.0_wp,handle)

    END SUBROUTINE blacs_make_basis

! *****************************************************************************

  SUBROUTINE blacs_get_element(matrix,irow_global,icol_global,alpha)

!   Purpose: Get the BLACS matrix element (irow_global,icol_global).

!   History: - Creation (22.01.2002, Matthias Krack)
!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix
    REAL(wp), INTENT(OUT)                     :: alpha
    INTEGER, INTENT(IN)                       :: icol_global,&
                                                 irow_global

!   *** Local variables ***

    INTEGER :: icol_local,ipcol,iprow,irow_local,mypcol,myprow,npcol,&
               nprow

    INTEGER, DIMENSION(9) :: desca

    REAL(wp), DIMENSION(:,:), POINTER :: a
    TYPE(cp_blacs_env_type), POINTER :: context

!   ---------------------------------------------------------------------------

    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

#if defined(__SCALAPACK)

    a => matrix%local_data
    desca(:) = matrix%matrix_struct%descriptor(:)

    CALL infog2l(irow_global,icol_global,desca,nprow,npcol,myprow,mypcol,&
                 irow_local,icol_local,iprow,ipcol)

    IF ((iprow == myprow).AND.(ipcol == mypcol)) THEN
      alpha = a(irow_local,icol_local)
    ELSE
      alpha = 0.0_wp
    END IF

#else

    alpha = matrix%local_data(irow_global,icol_global)

#endif
  END SUBROUTINE blacs_get_element

! *****************************************************************************

!
! copies ncol vectors out of matrix_v starting at firstcol to the
! replicated vector rep_v (or the other way around), it is assumed that
! all PEs have all replicated data, firstcol is a blacs_block boundary
! (i.e. ncol_block*k+1)

!
! attention ! rep_v(block_size,nrows) for increased locality in matrix multiply
! (i.e. vectors are stored as rows of this matrix)

!
! direction="B2R" .or. "R2B"  Blacs -> Replicated and ...
!
  SUBROUTINE blacs_replicated_copy(matrix_v,rep_v,firstcol,ncol,nblock,direction)
    TYPE(cp_full_matrix_type) , POINTER              :: matrix_v
    REAL(wp), DIMENSION(:,:), POINTER              :: rep_v
    INTEGER, INTENT(IN)                            :: firstcol,ncol,nblock
    CHARACTER(LEN=3)                               :: direction
    !----
    REAL(wp), DIMENSION(:,:), POINTER :: blacs_block
    INTEGER                           :: i,j,nrow_global,ncol_global,ncol_block
    INTEGER                           :: myprow,mypcol,this_col,nrow_local
    INTEGER                           :: ncol_local,nrow_block,handle
    INTEGER, DIMENSION(:), POINTER                 :: row_indices,col_indices
    TYPE(cp_blacs_env_type), POINTER :: context

    CALL timeset("blacs_replicated_copy","I","",handle)

    context => matrix_v%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)

    blacs_block=>matrix_v%local_data

    CALL cp_fm_get_info(matrix=matrix_v, &
                          nrow_global=nrow_global,ncol_global=ncol_global, &
                          nrow_block =nrow_block ,ncol_block =ncol_block , &
                          nrow_local =nrow_local ,ncol_local =ncol_local , &
                          row_indices=row_indices,col_indices=col_indices)
  
#if defined(__SCALAPACK) 
    IF (direction.EQ."B2R") THEN
       CALL dcopy(nblock*nrow_global,0.0_wp,0,rep_v(1,1),1)
    ENDIF
#endif

    DO i=1,ncol_local
       this_col=col_indices(i)-firstcol+1
       IF (this_col.GE.1 .AND. this_col.LE.ncol) THEN
          IF (direction.EQ."B2R") THEN
             DO j=1,nrow_local
               rep_v(this_col,row_indices(j))=blacs_block(j,i) 
             ENDDO
          ELSE
             DO j=1,nrow_local
                blacs_block(j,i)=rep_v(this_col,row_indices(j))
             ENDDO
          ENDIF
       ENDIF
    ENDDO

#if defined(__SCALAPACK) 
    IF (direction.EQ."B2R") THEN
       CALL mp_sum(rep_v,matrix_v%matrix_struct%para_env%group)
    ENDIF
#endif

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_replicated_copy

  SUBROUTINE blacs_column_copy(msource,mtarget,ncol,source_start,target_start)
    TYPE(cp_full_matrix_type), POINTER          :: msource,mtarget
    INTEGER, INTENT(IN)                       :: ncol
    INTEGER, INTENT(IN), OPTIONAL             :: source_start,target_start

    REAL(wp), DIMENSION(:,:), POINTER         :: a,b
    INTEGER                                   :: nprow,npcol,&
                                                 myprow,mypcol,i,n
    INTEGER, DIMENSION(9)                     :: desca,descb
    TYPE(cp_blacs_env_type), POINTER          :: context
    INTEGER ss,ts
    ss=1
    ts=1
    IF (PRESENT(source_start)) ss=source_start
    IF (PRESENT(target_start)) ts=target_start

    context => msource%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    n = msource%matrix_struct%nrow_global

    a => msource%local_data
    b => mtarget%local_data

#if defined(__SCALAPACK)
    desca(:) = msource%matrix_struct%descriptor(:)
    descb(:) = mtarget%matrix_struct%descriptor(:)
    DO i=0,ncol-1
      CALL pdcopy(n,a(1,1),1,ss+i,desca,1,b(1,1),1,ts+i,descb,1)
    ENDDO
#else
    CALL dcopy(ncol*n,a(1,ss),1,b(1,ts),1)
#endif

  END SUBROUTINE blacs_column_copy

!*************************************************************************
!computes all eigenvalues and vectors of a real symmetric matrix
!should be quite a bit faster than syevx for that case
!especially in parallel with thightly clustered evals
!needs more workspace
!*************************************************************************
  SUBROUTINE blacs_syevd(matrix,eigenvectors,eigenvalues,error)

    TYPE(cp_full_matrix_type), POINTER          :: eigenvectors,matrix
    REAL(wp), DIMENSION(:), INTENT(OUT)       :: eigenvalues
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    INTEGER  :: handle,info,istat,liwork,lwork,mypcol,myprow,n

    INTEGER, DIMENSION(9) :: descm,descv
    REAL(wp), DIMENSION(:,:), POINTER :: m,v
    REAL(wp), DIMENSION(:), POINTER :: work
    INTEGER, DIMENSION(:), POINTER  :: iwork
    TYPE(cp_blacs_env_type), POINTER          :: context
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='blacs_syevd',&
         routineP=moduleN//':'//routineN
!!FM    integer :: muno1,muno2

    CALL timeset("blacs_syevd","I","",handle)
    failure=.FALSE.
   
    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    n = matrix%matrix_struct%nrow_global

    m => matrix%local_data

#if defined(__SCALAPACK)
    descm(:) = matrix%matrix_struct%descriptor(:)
    v => eigenvectors%local_data
    descv(:) = eigenvectors%matrix_struct%descriptor(:)

    liwork=7*n+8*context%num_pe(2)+2
    ALLOCATE(iwork(liwork),STAT=istat)
    IF (istat.NE.0) CALL stop_memory("blacs_syevd","iwork")
! work space query
    lwork = -1
    ALLOCATE(work(1),STAT=istat)
    IF (istat.NE.0) CALL stop_memory("blacs_syevd","work")
    CALL PDSYEVD('V','U',n,m(1,1),1,1,descm,eigenvalues(1),v(1,1),1,1,descv, &
                  work(1),lwork,iwork(1),liwork,info)
    ! look here for a PDORMTR warning :-)
    ! this routine seems to need more workspace than reported
    ! (? lapack bug ...)
    ! arbitrary additional memory  ... we give 1000 more words
    lwork = NINT(work(1)+1000)
    DEALLOCATE(work)
    ALLOCATE(work(lwork),STAT=istat)
    IF (istat.NE.0) CALL stop_memory("blacs_syevd","lwork")
    CALL PDSYEVD('V','U',n,m(1,1),1,1,descm,eigenvalues(1),v(1,1),1,1,descv, &
                  work(1),lwork,iwork(1),liwork,info)
#else
    ! for n=4 with my lapack (linux) I need to add 50 to the value of lwork.
    ! unfortunately query is not supported. I did not test n dependence[fawzi]
    lwork=1+6*n+2*n**2+50
    liwork=5*n+3
    ALLOCATE(work(lwork),STAT=istat)
    IF (istat.NE.0) CALL stop_memory("blacs_syevd","work")
    ALLOCATE(iwork(liwork),STAT=istat)
    IF (istat.NE.0) CALL stop_memory("blacs_syevd","iwork")
    CALL DSYEVD('V','U', n, m(1,1), n, eigenvalues(1), work(1), lwork, iwork(1), liwork, info)
    CALL copy_blacs_to_blacs_matrix(matrix,eigenvectors)
#endif

    IF (info.NE.0) CALL stop_program("blacs_dsyevd","problems")

    DEALLOCATE(iwork)
    DEALLOCATE(work)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_syevd

! *****************************************************************************

  SUBROUTINE blacs_syevx(matrix,eigenvectors,eigenvalues,neig,work_syevx,error)
!  alternative , needs dsygvx
!  SUBROUTINE blacs_syevx(matrix,eigenvectors,eigenvalues,neig,work_syevx,&
!                         globenv,bmatrix)

!   Purpose: Diagonalise the symmetric n by n matrix using the LAPACK library.
!            if bmatrix present, use general solver, will overwrite bmatrix
!            with its cholesky decom.

!   History: - Creation (06.06.2001, Matthias Krack)
!   History: - Mod      (05.2001, Joost VandeVondele)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: eigenvectors,matrix
!    alternative needs dsygvx
!    TYPE(cp_full_matrix_type), POINTER, OPTIONAL:: bmatrix
    REAL(wp), INTENT(IN)                      :: work_syevx
    INTEGER, INTENT(IN)                       :: neig
    REAL(wp), DIMENSION(:), INTENT(OUT)       :: eigenvalues
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "blacs_syevx",&
      routineP=moduleN//':'//routineN
    REAL(wp), PARAMETER :: abstol = 0.0_wp,&
                           orfac = -1.0_wp,&
                           vl = 0.0_wp,&
                           vu = 0.0_wp

!   *** Local variables ***

    INTEGER  :: handle,info,istat,liwork,lwork,m,mypcol,myprow,n,nb,&
                nn,np0,npcol,npe,nprow,nq0,nz,output_unit,itype
    LOGICAL  :: ionode

    INTEGER, DIMENSION(9) :: desca,descz,descb

    REAL(wp), DIMENSION(:), ALLOCATABLE :: gap,w,work
    INTEGER, DIMENSION(:), ALLOCATABLE  :: iclustr,ifail,iwork
    REAL(wp), DIMENSION(:,:), POINTER   :: a,z,b
    TYPE(cp_blacs_env_type), POINTER    :: context
    TYPE(cp_logger_type), POINTER       :: logger

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL  :: iceil,numroc
#else
    INTEGER, EXTERNAL  :: ilaenv
#endif

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_syevx","I","",handle)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source

    context =>  matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    
    n = matrix%matrix_struct%nrow_global

    ALLOCATE (w(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"w",n*wp_size)
    w(:) = 0.0_wp

    eigenvalues(:) = 0.0_wp
#if defined(__SCALAPACK)

    IF (matrix%matrix_struct%nrow_block /= matrix%matrix_struct%ncol_block) THEN
       CALL stop_program(routineP,"Invalid blocksize (no square blocks)")
    END IF

    a => matrix%local_data
    desca(:) = matrix%matrix_struct%descriptor(:)
    z => eigenvectors%local_data
    descz(:) = eigenvectors%matrix_struct%descriptor(:)

!   *** Get the optimal work storage size ***

    npe = nprow*npcol
    nb = matrix%matrix_struct%nrow_block
    nn = MAX(n,nb,2)
    np0 = numroc(nn,nb,0,0,nprow)
    nq0 = MAX(numroc(nn,nb,0,0,npcol),nb)

    lwork = 5*n + MAX(5*nn,np0*nq0) + iceil(neig,npe)*nn + 2*nb*nb +&
            INT(work_syevx*REAL((neig - 1)*n,wp))
    liwork = MAX(3*n + npe + 1,4*n,14) + 2*n

    ALLOCATE (gap(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"gap",npe*wp_size)
    gap(:) = 0.0_wp
    ALLOCATE (iclustr(2*npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"iclustr",2*npe*int_size)
    iclustr(:) = 0
    ALLOCATE (ifail(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"ifail",n*int_size)
    ifail(:) = 0
    ALLOCATE (iwork(liwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"iwork",liwork*int_size)
    iwork(:) = 0
    ALLOCATE (work(lwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"work",lwork*wp_size)
    work(:) = 0.0_wp
    ! write(6,*) "done allocate"

!   *** Diagonalise matrix ***
!    if (PRESENT(bmatrix)) then
!       itype=1
!       b => bmatrix%local_data
!       descb(:) = matrix%matrix_struct%descriptor(:)
!       CALL pdsygvx(itype,"V","I","U",n,a,1,1,desca,b,1,1,descb,vl,vu,1,neig, &
!                 abstol,m,nz,w,orfac,&
!                 z,1,1,descz,work,lwork,iwork,liwork,ifail,iclustr,gap,info)
!    else
       CALL pdsyevx("V","I","U",n,a(1,1),1,1,desca,vl,vu,1,neig,abstol,m,nz,w(1),orfac,&
                 z(1,1),1,1,descz,work(1),lwork,iwork(1),liwork,ifail(1),iclustr(1),gap,info)
!    endif
    ! write(6,*) "done pdsyevx"


!   *** Error handling ***

    IF (info /= 0) THEN
       output_unit = cp_logger_get_unit_nr(logger,&
            cp_failure_level,routineP,local=.FALSE.)
      IF (ionode) THEN
        WRITE (unit=output_unit,FMT="(/,(T3,A,T12,1X,I10))")&
          "info    = ",info,&
          "lwork   = ",lwork,&
          "liwork  = ",liwork,&
          "nz      = ",nz
        IF (info > 0) THEN
          WRITE (unit=output_unit,FMT="(/,T3,A,(T12,6(1X,I10)))")&
            "ifail   = ",ifail
          WRITE (unit=output_unit,FMT="(/,T3,A,(T12,6(1X,I10)))")&
            "iclustr = ",iclustr
          WRITE (unit=output_unit,FMT="(/,T3,A,(T12,6(1X,E10.3)))")&
            "gap     = ",gap
        END IF
      END IF
      CALL stop_program(routineP,"Error in pdsyevx")
    END IF

!   *** Release work storage ***

    DEALLOCATE (gap,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"gap")
    DEALLOCATE (iclustr,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"iclustr")
    DEALLOCATE (ifail,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"ifail")
    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"iwork")
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"work")

#else

    a => matrix%local_data
    z => eigenvectors%local_data

!   *** Get the optimal work storage size ***

    nb = MAX(ilaenv(1,"DSYTRD","U",n,-1,-1,-1),&
             ilaenv(1,"DORMTR","U",n,-1,-1,-1))

    lwork = MAX((nb + 3)*n,8*n)+n ! sun bug fix
    liwork = 5*n

    ALLOCATE (ifail(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"ifail",n*int_size)
    ifail(:) = 0
    ALLOCATE (iwork(liwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"iwork",liwork*int_size)
    iwork(:) = 0
    ALLOCATE (work(lwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"work",lwork*wp_size)
    work(:) = 0.0_wp

!   *** Diagonalise matrix ***
!   dsygvx seems not to be in all lapack versions ...
!    if (PRESENT(bmatrix)) then
!       itype=1
!       b => bmatrix%local_data
!       CALL dsygvx(itype,"V","I","U",n,a,n,b,n,vl,vu,1,neig,abstol,m,w,z,n, &
!                   work,lwork,iwork,ifail,info)
!    else
       CALL dsyevx("V","I","U",n,a(1,1),n,vl,vu,1,neig,abstol,m,w,z(1,1),n,work(1),lwork,&
                iwork(1),ifail(1),info)
!    endif

!   *** Error handling ***

    IF (info /= 0) CALL stop_program(routineP,"Error in dsyevx")

!   *** Release work storage ***

    DEALLOCATE (ifail,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"ifail")
    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"iwork")
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"work")

#endif
    eigenvalues(1:neig) = w(1:neig)
    DEALLOCATE (w,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"w")

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_syevx

! *****************************************************************************

  SUBROUTINE copy_blacs_to_blacs_matrix(source,TARGET,error)

!   Purpose: Copy BLACS matrix to a BLACS matrix of the same type.

!   History: - Creation (08.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER :: source,TARGET
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_blacs_to_blacs_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    INTEGER :: ipcol,iprow,istat,npcol,nprow,ncol_local,nrow_local,msize1,msize2
    REAL(wp), DIMENSION(:,:), POINTER :: a,b
    LOGICAL :: failure

!   ---------------------------------------------------------------------------

    failure=.TRUE.
    nprow = source%matrix_struct%context%num_pe(1)
    npcol = source%matrix_struct%context%num_pe(2)

    IF(.NOT. cp2k_is_parallel .OR.&
         cp_fmstruct_equivalent(source%matrix_struct,&
         target%matrix_struct,error=error)) THEN
      CALL cp_assert(SIZE(source%local_data,1)==SIZE(target%local_data,1).AND.&
         SIZE(source%local_data,2)==SIZE(target%local_data,2),cp_failure_level,&
         cp_assertion_failed, routine, &
         "internal error, local_data has different sizes"//&
         CPSourceFileRef,&
         error=error)
      CALL dcopy(SIZE(source%local_data,1)*SIZE(source%local_data,2),&
           source%local_data(1,1),1,target%local_data(1,1),1)
   ELSE
      CALL cp_assert(source%matrix_struct%nrow_global==&
           target%matrix_struct%nrow_global,&
           cp_failure_level,&
           cp_assertion_failed, routine, &
           "cannot copy between blacs matrixes of differen sizes"//&
           CPSourceFileRef,&
           error=error)
      CALL cp_assert(source%matrix_struct%ncol_global==&
           target%matrix_struct%ncol_global,&
           cp_failure_level,&
           cp_assertion_failed, routine, &
           "cannot copy between blacs matrixes of differen sizes"//&
           CPSourceFileRef,&
           error=error)
#ifdef __SCALAPACK
      CALL pdcopy(source%matrix_struct%nrow_global*&
           source%matrix_struct%ncol_global,&
           source%local_data(1,1),1,1,source%matrix_struct%descriptor,1,&
           target%local_data(1,1),1,1,target%matrix_struct%descriptor,1)
#else
      CPAssert(.FALSE.,cp_failure_level,routine,error,failure)
#endif
   END IF
  END SUBROUTINE copy_blacs_to_blacs_matrix

! *****************************************************************************

  SUBROUTINE copy_blacs_to_full_matrix(blacs_matrix,full_matrix,error)

!   Purpose: Copy a BLACS matrix to a full matrix.

!   History: - Creation (18.06.2001, Matthias Krack)
!   - 08.2002 changed for the transition to local_data, unoptimized 
!     (but used for debugging only...) [fawzi]
!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: blacs_matrix
    REAL(wp), DIMENSION(:,:), POINTER         :: full_matrix
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_blacs_to_full_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    INTEGER :: handle,icol_global,icol_local,ipcol,ipe,iprow,&
               irow_global,irow_local,istat,mypcol,mype,myprow,ncol_block,&
               ncol_global,ncol_local,npcol,npe,nprow,nrow_block,nrow_global,&
               nrow_local, source
    LOGICAL :: ionode

    REAL(wp), DIMENSION(:,:), POINTER :: blacs_block
    TYPE(cp_blacs_env_type), POINTER :: context

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: blacs_pnum,indxl2g

#endif
!   ---------------------------------------------------------------------------

    CALL timeset("copy_blacs_to_full_matrix","I","",handle)

    ! to do : use new io based on error
    source = blacs_matrix%matrix_struct%para_env%source
    ionode = blacs_matrix%matrix_struct%para_env%mepos==source

    context => blacs_matrix%matrix_struct%context

    nrow_global = blacs_matrix%matrix_struct%nrow_global
    ncol_global = blacs_matrix%matrix_struct%ncol_global

    IF (ionode) THEN
      ALLOCATE (full_matrix(nrow_global,ncol_global),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"full_matrix",&
                                       nrow_global*ncol_global*wp_size)
      full_matrix(:,:) = 0.0_wp
    END IF

    mype=context%my_pid
    npe=context%n_pid
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    nrow_block = blacs_matrix%matrix_struct%nrow_block
    ncol_block = blacs_matrix%matrix_struct%ncol_block

#if defined(__SCALAPACK)

    DO iprow=0,nprow-1
      DO ipcol=0,npcol-1

        ipe = blacs_pnum(context%group,iprow,ipcol)

        nrow_local = blacs_matrix%matrix_struct%nrow_local(iprow)
        ncol_local = blacs_matrix%matrix_struct%ncol_local(ipcol)

        IF (ionode) THEN

          IF (ipe /= mype) THEN
            ALLOCATE (blacs_block(nrow_local, ncol_local),&
                      STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"alloc blacs_block",&
                               nrow_local*ncol_local*wp_size)
            END IF
            CALL dgerv2d(context%group,nrow_local,ncol_local,&
                         blacs_block(1,1),nrow_local,&
                         iprow,ipcol)
          ELSE
            blacs_block => blacs_matrix%local_data
          END IF

          DO icol_local=1,ncol_local
            icol_global = indxl2g(icol_local,ncol_block,ipcol,&
                 blacs_matrix%matrix_struct%first_p_pos(2),npcol)
            DO irow_local=1,nrow_local
              irow_global = indxl2g(irow_local,nrow_block,iprow,&
                   blacs_matrix%matrix_struct%first_p_pos(1),nprow)
              full_matrix(irow_global,icol_global) = blacs_block(irow_local,&
                                                                 icol_local)
            END DO
          END DO

          IF (ipe /= mype) THEN
            DEALLOCATE (blacs_block,STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"dealloc blacs_block")
            END IF
          END IF

        ELSE

          IF (ipe == mype) THEN
            CALL dgesd2d(context%group,nrow_local,ncol_local,&
                         blacs_matrix%local_data(1,1),nrow_local,&
                         blacs_matrix%matrix_struct%first_p_pos(1),&
                         blacs_matrix%matrix_struct%first_p_pos(2))
          END IF

        END IF

        CALL blacs_barrier(context%group,"A")

      END DO
    END DO

#else

    full_matrix(:,:) = blacs_matrix%local_data(:,:)

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE copy_blacs_to_full_matrix

! *****************************************************************************

  SUBROUTINE copy_blacs_to_sparse_matrix(blacs_matrix,sparse_matrix)

!   Purpose: Copy a BLACS matrix to a sparse matrix. The BLACS matrix blocks
!            are deallocated during the copy procedure.

!   History: - Creation (06.06.2001, Matthias Krack)
!   08.2002 adapted to local_data, could be optimized
!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: blacs_matrix
    TYPE(real_matrix_type), POINTER           :: sparse_matrix

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_blacs_to_sparse_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: group,handle,iblock_col,iblock_row,icol,icol_global,&
               icol_local,ipcol,ipe,iprow,irow,irow_global,irow_local,istat,&
               jpcol,jprow,mypcol,mype,myprow,nblock_row,ncol_block,&
               ncol_local,npcol,npe,nprow,nrow_block,nrow_local

    INTEGER, DIMENSION(:), POINTER    :: first_col,first_row,last_col,last_row
    REAL(wp), DIMENSION(:,:), POINTER :: blacs_block,sparse_block
    TYPE(cp_blacs_env_type), POINTER  :: context

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: blacs_pnum,indxg2l,indxg2p

#endif
!   ---------------------------------------------------------------------------

    CALL timeset("copy_blacs_to_sparse_matrix","I","",handle)

    group = blacs_matrix%matrix_struct%para_env%group
    context => blacs_matrix%matrix_struct%context

    CALL get_matrix_info(matrix=sparse_matrix,&
                         nblock_row=nblock_row,&
                         first_row=first_row,&
                         first_col=first_col,&
                         last_row=last_row,&
                         last_col=last_col)
    mype=context%my_pid
    npe=context%n_pid
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    nrow_block = blacs_matrix%matrix_struct%nrow_block
    ncol_block = blacs_matrix%matrix_struct%ncol_block

#if defined(__SCALAPACK)

    DO iprow=0,nprow-1
      DO ipcol=0,npcol-1

        ipe = blacs_pnum(context%group,iprow,ipcol)

        IF (ipe /= mype) THEN
          nrow_local = blacs_matrix%matrix_struct%nrow_local(iprow)
          ncol_local = blacs_matrix%matrix_struct%ncol_local(ipcol)
          ALLOCATE (blacs_block(nrow_local,ncol_local),&
                    STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"alloc blacs_block",&
                             nrow_local*ncol_local*wp_size)
          END IF
        ELSE
          blacs_block => blacs_matrix%local_data
        END IF


        CALL mp_bcast(blacs_block,ipe,group)

        DO iblock_row=1,nblock_row

          block_node => first_block_node(matrix=sparse_matrix,&
                                         block_row=iblock_row)

          DO WHILE (ASSOCIATED(block_node))

            CALL get_block_node(block_node=block_node,&
                                block_col=iblock_col,&
                                BLOCK=sparse_block)

            icol = 1

            DO icol_global=first_col(iblock_col),last_col(iblock_col)

              jpcol = indxg2p(icol_global,ncol_block,mypcol,&
                   blacs_matrix%matrix_struct%first_p_pos(2),npcol)

              IF (jpcol == ipcol) THEN

                icol_local = indxg2l(icol_global,ncol_block,mypcol,&
                     blacs_matrix%matrix_struct%first_p_pos(2),npcol)

                irow = 1

                DO irow_global=first_row(iblock_row),last_row(iblock_row)

                  jprow = indxg2p(irow_global,nrow_block,myprow,&
                       blacs_matrix%matrix_struct%first_p_pos(1),nprow)

                  IF (jprow == iprow) THEN

                     irow_local = indxg2l(irow_global,nrow_block,myprow,&
                          blacs_matrix%matrix_struct%first_p_pos(1),nprow)

                    sparse_block(irow,icol) = blacs_block(irow_local,&
                                                          icol_local)

                  END IF

                  irow = irow + 1

                END DO

              END IF

              icol = icol + 1

            END DO

            block_node => next_block_node(block_node)

          END DO

        END DO

        IF (ipe /= mype) THEN
          DEALLOCATE (blacs_block,STAT=istat)
          IF (istat /= 0) THEN
             CALL stop_memory(routine,"dealloc blacs_block")
          END IF
        END IF

      END DO
    END DO

#else

    blacs_block => blacs_matrix%local_data

    DO iblock_row=1,nblock_row

      block_node => first_block_node(matrix=sparse_matrix,&
                                     block_row=iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            BLOCK=sparse_block)

        icol = 1

        DO icol_global=first_col(iblock_col),last_col(iblock_col)

          irow = 1

          DO irow_global=first_row(iblock_row),last_row(iblock_row)

            sparse_block(irow,icol) = blacs_block(irow_global,icol_global)

            irow = irow + 1

          END DO

          icol = icol + 1

        END DO

        block_node => next_block_node(block_node)

      END DO

    END DO

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE copy_blacs_to_sparse_matrix

! *****************************************************************************

  SUBROUTINE copy_sparse_to_blacs_matrix(sparse_matrix,blacs_matrix)

!   Purpose: Copy a sparse matrix to a BLACS matrix. The BLACS matrix blocks
!            are allocated during the copy procedure.

!   History: - Creation (05.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: blacs_matrix
    TYPE(real_matrix_type), POINTER           :: sparse_matrix

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_sparse_to_blacs_matrix (MODULE qs_blacs)",&
      routineN="copy_sparse_to_blacs_matrix"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: group,handle,iblock_col,iblock_row,icol,icol_global,&
               icol_local,ipcol,ipe,iprow,irow,irow_global,irow_local,istat,&
               jpcol,jprow,mypcol,mype,myprow,nblock_row,ncol_block,&
               ncol_local,npcol,npe,nprow,nrow_block,nrow_local

    INTEGER, DIMENSION(:), POINTER    :: first_col,first_row,last_col,last_row
    REAL(wp), DIMENSION(:,:), POINTER :: blacs_block,sparse_block
    TYPE(cp_blacs_env_type), POINTER  :: context

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: blacs_pnum,indxg2l,indxg2p

#endif
!   ---------------------------------------------------------------------------
    
    CALL timeset("copy_sparse_to_blacs_matrix","I","",handle)

    group = blacs_matrix%matrix_struct%para_env%group
    context => blacs_matrix%matrix_struct%context

    mype=context%my_pid
    npe=context%n_pid
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    nrow_block = blacs_matrix%matrix_struct%nrow_block
    ncol_block = blacs_matrix%matrix_struct%ncol_block

    CALL get_matrix_info(matrix=sparse_matrix,&
                         nblock_row=nblock_row,&
                         first_row=first_row,&
                         first_col=first_col,&
                         last_row=last_row,&
                         last_col=last_col)
#if defined(__SCALAPACK)

    DO iprow=0,nprow-1
      DO ipcol=0,npcol-1

        ipe = blacs_pnum(context%group,iprow,ipcol)

        IF (ipe /= mype) THEN
          nrow_local = blacs_matrix%matrix_struct%nrow_local(iprow)
          ncol_local = blacs_matrix%matrix_struct%ncol_local(ipcol)
          ALLOCATE (blacs_block(nrow_local,ncol_local),&
                    STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"alloc blacs_block",&
                             nrow_local*ncol_local*wp_size)
          END IF
        ELSE
           blacs_block => blacs_matrix%local_data
        END IF

        blacs_block(:,:) = 0.0_wp

        DO iblock_row=1,nblock_row

          block_node => first_block_node(matrix=sparse_matrix,&
                                         block_row=iblock_row)

          DO WHILE (ASSOCIATED(block_node))

            CALL get_block_node(block_node=block_node,&
                                block_col=iblock_col,&
                                BLOCK=sparse_block)

            icol = 1

            DO icol_global=first_col(iblock_col),last_col(iblock_col)

               jpcol = indxg2p(icol_global,ncol_block,mypcol,&
                    blacs_matrix%matrix_struct%first_p_pos(2),npcol)

              IF (jpcol == ipcol) THEN

                icol_local = indxg2l(icol_global,ncol_block,mypcol,&
                     blacs_matrix%matrix_struct%first_p_pos(2),npcol)

                irow = 1

                DO irow_global=first_row(iblock_row),last_row(iblock_row)

                  jprow = indxg2p(irow_global,nrow_block,myprow,&
                       blacs_matrix%matrix_struct%first_p_pos(1),nprow)

                  IF (jprow == iprow) THEN

                    irow_local = indxg2l(irow_global,nrow_block,myprow,&
                         blacs_matrix%matrix_struct%first_p_pos(1),nprow)

                    blacs_block(irow_local,icol_local) = sparse_block(irow,&
                                                                      icol)

                  END IF

                  irow = irow + 1

                END DO

              END IF

              icol = icol + 1

            END DO

            block_node => next_block_node(block_node)

          END DO

        END DO

        CALL mp_sum(blacs_block,ipe,group)

        IF (ipe /= mype) THEN
          DEALLOCATE (blacs_block,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"dealloc blacs_block")
          END IF
        END IF

      END DO
    END DO

#else

    IF (.NOT.ASSOCIATED(blacs_matrix%local_data)) THEN
      nrow_local = blacs_matrix%matrix_struct%nrow_local(myprow)
      ncol_local = blacs_matrix%matrix_struct%ncol_local(mypcol)
      ALLOCATE (blacs_matrix%local_data(nrow_local,ncol_local),&
                STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"blacs_matrix%local_data",&
                         nrow_local*ncol_local*wp_size)
      END IF
    END IF

    blacs_block => blacs_matrix%local_data

    blacs_block(:,:) = 0.0_wp

    DO iblock_row=1,nblock_row

      block_node => first_block_node(matrix=sparse_matrix,&
                                     block_row=iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            BLOCK=sparse_block)

        icol = 1

        DO icol_global=first_col(iblock_col),last_col(iblock_col)

          irow = 1

          DO irow_global=first_row(iblock_row),last_row(iblock_row)

            blacs_block(irow_global,icol_global) = sparse_block(irow,icol)

            irow = irow + 1

          END DO

          icol = icol + 1

        END DO

        block_node => next_block_node(block_node)

      END DO

    END DO

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE copy_sparse_to_blacs_matrix

! *****************************************************************************

  SUBROUTINE get_blacs_info(globenv,my_process_row,my_process_column,&
                            my_process_number,number_of_process_rows,&
                            number_of_process_columns,number_of_processes)

!   Purpose: Return informations about the specified BLACS context.

!   History: - Creation (19.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, OPTIONAL, INTENT(OUT)            :: my_process_column,&
                                                 my_process_number,&
                                                 my_process_row,&
                                                 number_of_process_columns,&
                                                 number_of_process_rows,&
                                                 number_of_processes

!   *** Local variables ***

    INTEGER :: mypcol,mype,myprow,npcol,npe,nprow
    TYPE(cp_blacs_env_type), POINTER :: context

!   ---------------------------------------------------------------------------

    context => cp_blacs_env_from_globenv(globenv)
    IF (.NOT.ASSOCIATED(context)) CALL stop_program("get_blacs_info","no blacs environment")

    mype=context%my_pid
    npe=context%n_pid
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

    IF (PRESENT(my_process_row)) my_process_row = myprow
    IF (PRESENT(my_process_column)) my_process_column = mypcol
    IF (PRESENT(my_process_number)) my_process_number = mype
    IF (PRESENT(number_of_process_rows)) number_of_process_rows = nprow
    IF (PRESENT(number_of_process_columns)) number_of_process_columns = npcol
    IF (PRESENT(number_of_processes)) number_of_processes = npe

  END SUBROUTINE get_blacs_info

! *****************************************************************************


  SUBROUTINE power_blacs_matrix(matrix,work,exponent,threshold,n_dependent,&
                                work_syevx)

!   Purpose: Raise the real symmetric n by n matrix to the power given by
!            exponent. All eigenvectors with a corresponding eigenvalue lower
!            than threshold are quenched.

!   History: - Creation (29.03.1999, Matthias Krack)
!            - Parallelised using BLACS and ScaLAPACK (06.06.2001, MK)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix,work
    REAL(wp), INTENT(IN)                      :: exponent,threshold,work_syevx
    INTEGER, INTENT(OUT)                      :: n_dependent

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE power_blacs_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    REAL(wp) :: f,p
    INTEGER  :: handle,icol_global,icol_local,ipcol,iprow,irow_global,&
                irow_local,istat,mypcol,myprow,ncol_block,ncol_global,npcol,&
                nprow,nrow_block,nrow_global

    REAL(wp), DIMENSION(:), ALLOCATABLE :: eigenvalues
    REAL(wp), DIMENSION(:,:), POINTER   :: eigenvectors
    TYPE(cp_blacs_env_type), POINTER    :: context

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: indxg2l,indxg2p

#endif
!   ---------------------------------------------------------------------------

    CALL timeset("power_blacs_matrix","I","",handle)

    
    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    n_dependent = 0
    p = 0.5_wp*exponent

    nrow_global = matrix%matrix_struct%nrow_global
    ncol_global = matrix%matrix_struct%ncol_global

    ALLOCATE (eigenvalues(ncol_global),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"eigenvalues",ncol_global*wp_size)
    eigenvalues(:) = 0.0_wp

!   *** Compute the eigenvectors and eigenvalues ***

    CALL blacs_syevx(matrix,work,eigenvalues,ncol_global,work_syevx)

#if defined(__SCALAPACK)

    nrow_block = work%matrix_struct%nrow_block
    ncol_block = work%matrix_struct%ncol_block

    eigenvectors => work%local_data

!   *** Build matrix**exponent with eigenvector quenching ***

    p = 0.5_wp*exponent

    n_dependent = 0

    DO icol_global=1,ncol_global

      IF (eigenvalues(icol_global) < threshold) THEN

        n_dependent = n_dependent + 1

        ipcol = indxg2p(icol_global,ncol_block,mypcol,&
             work%matrix_struct%first_p_pos(2),npcol)

        IF (mypcol == ipcol) THEN
          icol_local = indxg2l(icol_global,ncol_block,mypcol,&
               work%matrix_struct%first_p_pos(2),npcol)
          DO irow_global=1,nrow_global
            iprow = indxg2p(irow_global,nrow_block,myprow,&
                 work%matrix_struct%first_p_pos(1),nprow)
            IF (myprow == iprow) THEN
              irow_local = indxg2l(irow_global,nrow_block,myprow,&
                   work%matrix_struct%first_p_pos(1),nprow)
              eigenvectors(irow_local,icol_local) = 0.0_wp
            END IF
          END DO
        END IF

      ELSE

        f = eigenvalues(icol_global)**p

        ipcol = indxg2p(icol_global,ncol_block,mypcol,&
             work%matrix_struct%first_p_pos(2),npcol)

        IF (mypcol == ipcol) THEN
          icol_local = indxg2l(icol_global,ncol_block,mypcol,&
               work%matrix_struct%first_p_pos(2),npcol)
          DO irow_global=1,nrow_global
            iprow = indxg2p(irow_global,nrow_block,myprow,&
                 work%matrix_struct%first_p_pos(1),nprow)
            IF (myprow == iprow) THEN
              irow_local = indxg2l(irow_global,nrow_block,myprow,&
                   work%matrix_struct%first_p_pos(1),nprow)
              eigenvectors(irow_local,icol_local) =&
                f*eigenvectors(irow_local,icol_local)
            END IF
          END DO
        END IF

      END IF

    END DO

#else

    eigenvectors => work%local_data

!   *** Build matrix**exponent with eigenvector quenching ***

    DO icol_global=1,ncol_global

      IF (eigenvalues(icol_global) < threshold) THEN

        n_dependent = n_dependent + 1
        eigenvectors(1:nrow_global,icol_global) = 0.0_wp

      ELSE

        f = eigenvalues(icol_global)**p
        eigenvectors(1:nrow_global,icol_global) =&
          f*eigenvectors(1:nrow_global,icol_global)

      END IF

    END DO

#endif
    CALL cp_fm_syrk("U","N",ncol_global,1.0_wp,work,1,1,0.0_wp,matrix)

    DEALLOCATE (eigenvalues,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"eigenvalues")

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE power_blacs_matrix

! *****************************************************************************

  SUBROUTINE read_blacs_matrix(matrix,lunit,globenv)

!   Purpose: Read a BLACS matrix from the logical unit number "lunit".

!   History: - Creation (19.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: lunit

!   *** Local variables ***

    INTEGER :: i,j,mypcol,myprow,ncol_local,nrow_local
    TYPE(cp_blacs_env_type), POINTER :: context
!   ---------------------------------------------------------------------------

    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)

    nrow_local = matrix%matrix_struct%nrow_local(myprow)
    ncol_local = matrix%matrix_struct%ncol_local(mypcol)

    READ (UNIT=lunit) ((matrix%local_data(i,j),i=1,nrow_local),&
                                                           j=1,ncol_local)

  END SUBROUTINE read_blacs_matrix

! *****************************************************************************

  SUBROUTINE replicate_blacs_matrix(prototype_matrix,new_matrix,name)

!   Purpose: Allocate a distributed BLACS matrix using a prototype matrix.

!   History: - Creation (08.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER :: new_matrix,prototype_matrix
    CHARACTER(LEN=*), INTENT(IN)     :: name

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE replicate_blacs_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    TYPE(cp_blacs_env_type), POINTER :: context

!   ---------------------------------------------------------------------------

    CALL cp_fm_create(new_matrix,prototype_matrix%matrix_struct)
    CALL dcopy(SIZE(new_matrix%local_data,1)*&
         SIZE(new_matrix%local_data,2),&
         prototype_matrix%local_data(1,1),&
         1,new_matrix%local_data(1,1),1)

  END SUBROUTINE replicate_blacs_matrix

! *****************************************************************************

  SUBROUTINE blacs_transpose(matrix,matrixt)
! matrixt=matrix^T ! assumes no symmetry
!***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix,matrixt

!   *** Local variables ***

    INTEGER :: handle,icol_global,icol_local,ipcol,iprow,irow_global,&
               irow_local,mypcol,myprow,ncol_block,ncol_global,ncol_local,&
               npcol,nprow,nrow_block,nrow_global,nrow_local

    INTEGER, DIMENSION(9) :: desca,descc

    REAL(wp), DIMENSION(:,:), POINTER :: a,c
    INTEGER :: i,j
    TYPE(cp_blacs_env_type), POINTER :: context

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_transpose","I","",handle)

    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
   
    nrow_global = matrix%matrix_struct%nrow_global
    ncol_global = matrix%matrix_struct%ncol_global

    a => matrix%local_data
    c => matrixt%local_data

#if defined(__SCALAPACK)
    desca(:) = matrix%matrix_struct%descriptor(:)
    descc(:) = matrixt%matrix_struct%descriptor(:)
    CALL pdtran(nrow_global,ncol_global,1.0_wp,a(1,1),1,1,desca,0.0_wp,c(1,1),1,1,descc)
#else
    DO j=1,ncol_global
     DO i=1,nrow_global
        c(j,i)=a(i,j)
     ENDDO
    ENDDO
#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_transpose

! *****************************************************************************

  SUBROUTINE symmetrise_blacs_matrix(matrix,work)

!   Purpose: Symmetrise a symmetric BLACS matrix.

!   History: - Creation (12.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix,work

!   *** Local variables ***

    INTEGER :: handle,icol_global,icol_local,ipcol,iprow,irow_global,&
               irow_local,mypcol,myprow,ncol_block,ncol_global,ncol_local,&
               npcol,nprow,nrow_block,nrow_global,nrow_local

    INTEGER, DIMENSION(9) :: desca,descc

    REAL(wp), DIMENSION(:,:), POINTER :: a,c
    TYPE(cp_blacs_env_type), POINTER :: context

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: indxl2g

#endif
!   ---------------------------------------------------------------------------

    CALL timeset("symmetrise_blacs_matrix","I","",handle)

    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

#if defined(__SCALAPACK)

    nrow_global = matrix%matrix_struct%nrow_global
    ncol_global = matrix%matrix_struct%ncol_global

    nrow_block = matrix%matrix_struct%nrow_block
    ncol_block = matrix%matrix_struct%ncol_block

    nrow_local = matrix%matrix_struct%nrow_local(myprow)
    ncol_local = matrix%matrix_struct%ncol_local(mypcol)

    a => work%local_data
    desca(:) = work%matrix_struct%descriptor(:)
    c => matrix%local_data
    descc(:) = matrix%matrix_struct%descriptor(:)

    DO icol_local=1,ncol_local
      icol_global = indxl2g(icol_local,ncol_block,mypcol,&
           matrix%matrix_struct%first_p_pos(2),npcol)
      DO irow_local=1,nrow_local
        irow_global = indxl2g(irow_local,nrow_block,myprow,&
             matrix%matrix_struct%first_p_pos(1),nprow)
        IF (irow_global > icol_global) THEN
          c(irow_local,icol_local) = 0.0_wp
        ELSE IF (irow_global == icol_global) THEN
          c(irow_local,icol_local) = 0.5_wp*c(irow_local,icol_local)
        END IF
      END DO
    END DO

    a(:,:) = c(:,:)

    CALL pdtran(nrow_global,ncol_global,1.0_wp,a(1,1),1,1,desca,1.0_wp,c(1,1),1,1,descc)

#else

    a => matrix%local_data

    CALL symmetrize_matrix(a,"upper_to_lower")

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE symmetrise_blacs_matrix

! *****************************************************************************

  SUBROUTINE write_blacs_matrix(matrix,lunit,globenv,formatted)

!   Purpose: Write a BLACS matrix to the logical unit number "lunit".

!   History: - Creation (19.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN), OPTIONAL :: globenv
    INTEGER, INTENT(IN)                       :: lunit
    LOGICAL, INTENT(IN), OPTIONAL             :: formatted

!   *** Local variables ***

    INTEGER :: i,j,mypcol,myprow,ncol_local,npcol,nprow,nrow_local
    TYPE(cp_blacs_env_type), POINTER :: context

!   ---------------------------------------------------------------------------

    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)

    nrow_local = matrix%matrix_struct%nrow_local(myprow)
    ncol_local = matrix%matrix_struct%ncol_local(mypcol)

    IF (PRESENT(formatted)) THEN
       IF (formatted) THEN
         WRITE(6,*) matrix%name
         WRITE(6,*) "row ,col ,value"
         DO j=1,ncol_local
         DO i=1,nrow_local
         WRITE (UNIT=lunit,FMT='(I5,I5,E16.8)') i,j,matrix%local_data(i,j)
         ENDDO
         ENDDO
       ELSE
          WRITE (UNIT=lunit)&
            ((matrix%local_data(i,j),i=1,nrow_local),&
                                                 j=1,ncol_local)
       ENDIF
    ELSE
          WRITE (UNIT=lunit) ((matrix%local_data(i,j),i=1,nrow_local),&
                                                            j=1,ncol_local)

    ENDIF

  END SUBROUTINE write_blacs_matrix
!!***
! *****************************************************************************

!!****f* blacs/add_blacs_to_block_diag_sm [1.0] *
!!
!!   NAME
!!     add_blacs_to_block_diag_sm
!!
!!   FUNCTION
!!     add the diagonal blocs of a blacs matrix to a sparse matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     source_m: the blacs matrix whose diagonal blocks are copied
!!     target_m: the sparse matrix that will contain the diagonal
!!               blocks (the blocks should already be allocated)
!!     d_struct: the distribution of the the blocks among the processors:
!!               the number of the processor that has each block.
!!               If an associated pointer is given its contents should
!!               be rightly initialized, if not associated then
!!               after the call it will contain the actual distribution
!!               and you are responsible of deallocating it.
!!     error: variable to control error logging, stopping,...
!!            see module cp_error_handling
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE add_blacs_to_block_diag_sm(source_m, target_m, para_env, &
       d_struct, error)
    TYPE(cp_full_matrix_type), POINTER :: source_m
    TYPE(real_matrix_type), INTENT(inout), TARGET :: target_m
    TYPE(cp_para_env_type), POINTER :: para_env
    INTEGER, DIMENSION(:), POINTER, OPTIONAL :: d_struct
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    INTEGER :: handle, i, nrows, ncols, stat, j
    INTEGER, DIMENSION(:), POINTER :: my_d_struct, first_col,first_row,&
         last_col, last_row
    INTEGER :: start_proc_row,start_proc_col,n_blacs_proc, sm_block_proc,&
         nprow,npcol,myprow,mypcol,b_block_row_start,b_block_row_stop,&
         b_block_col_start,b_block_col_stop,sm_block_nr,b_block_row,&
         b_block_col, start_row_of_b, end_row_of_b, start_row_of_sm,&
         end_row_of_sm, start_col_of_b, end_col_of_b, start_col_of_sm,&
         end_col_of_sm,nblock_col
    CHARACTER(len=*), PARAMETER :: routineN='add_blacs_to_block_diag_sm',&
         routineP=moduleN//':'//routineN
    TYPE(real_matrix_type), POINTER :: target_ptr
    TYPE(cp_matrix_block_iterator) :: iterator
    TYPE(cp_blacs_env_type), POINTER :: context
    REAL(kind=wp), DIMENSION(:,:), POINTER :: block_val
    failure=.FALSE.
    target_ptr => target_m
    NULLIFY(my_d_struct)

    CALL timeset(routineN//','//moduleN,'I',"",handle)
    CPPrecondition(ASSOCIATED(target_ptr),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       context => source_m%matrix_struct%context
       CALL get_matrix_info(target_ptr, nrow=nrows, ncol=ncols,&
            nblock_col=nblock_col)
       CPAssert(ncols==nrows,cp_warning_level,routineP,error,failure)
       ncols=MIN(nrows,ncols)
       IF (PRESENT(d_struct)) my_d_struct => d_struct
       IF (.NOT.ASSOCIATED(my_d_struct)) THEN
! build distribution structure
! assume every diagonal block is on one (and just one) proc
          ALLOCATE(my_d_struct(ncols),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          my_d_struct=0
          CALL cp_init(iterator, matrix=target_m, error=error)
          DO
             IF (.NOT.cp_next(iterator)) EXIT
             CALL cp_get(iterator,block_row=i,block_col=j,error=error)
             IF (i==j) THEN
                my_d_struct(i)=para_env%mepos
             ELSE
                CALL cp_error_message(cp_warning_level,&
                     routineP,"there is an off diagonal block! in "//&
                     CPSourceFileRef,error)
! just as info (only on the local proc)
             END IF
          END DO
          CALL cp_dealloc_ref(iterator,error=error)
          CALL mp_sum(my_d_struct,para_env%group)
       END IF
    END IF
    CPAssert(ASSOCIATED(my_d_struct),cp_failure_level,routineP,error,failure)
    CPAssert(SIZE(my_d_struct)==nblock_col,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPAssert(ALL(my_d_struct<para_env%num_pe),cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
!group = global_env%group
!source = global_env%source

       CALL get_matrix_info(matrix=target_ptr,&
            first_row=first_row,&
            first_col=first_col,&
            last_row=last_row,&
            last_col=last_col)
#if (__SCALAPACK)
!MK       call blacs_pinfo(blacs_proc_nr, n_blacs_proc)
!MK       CPPrecondition(blacs_proc_nr==global_env%mepos,cp_failure_level,routineP,error,failure)
!MK       call blacs_gridinfo(source_m%context,nprow,npcol,myprow,mypcol)
!MK       start_proc_row=source_m%matrix_struct%descriptor(7)
!MK       start_proc_col=source_m%matrix_struct%descriptor(8)
!MK!    if ( (myprow >= n_proc_row) .or. (mypcol >= n_proc_col) ) then
!MK!       has_b_matrix=.false.
!MK!    else
!MK!       has_b_matrix=.true.
!MK!    end if
!MK
!MK       do sm_block_nr=1,ncols
!MK          CPPrecondition(first_row(sm_block_nr)==first_col(sm_block_nr),cp_warning_level,routineP,error,failure)
!MK! the processor that has the sm_block
!MK          sm_block_proc=d_struct(sm_block_nr)
!MK
!MK! 0 based indexing
!MK          b_block_row_start= (first_row(sm_block_nr)-1)/source_m%nrow_block
!MK          b_block_row_stop= (last_row(sm_block_nr)-1)/source_m%nrow_block
!MK          b_block_col_start= (first_col(sm_block_nr)-1)/source_m%ncol_block
!MK          b_block_col_stop= (last_col(sm_block_nr)-1)/source_m%ncol_block
!MK
!MK          do b_block_row= b_block_row_start,b_block_row_stop
!MK             do b_block_col= b_block_col_start,b_block_col_stop
!MK                blacs_block_proc_row= mod(b_block_row+start_proc_row,nprow)
!MK                blacs_block_proc_col= mod(b_block_col+start_proc_col,npcol)
!MK                blacs_proc_nr= blacs_pnum(source_matrix%context,&
!MK                     blacs_block_proc_row, blacs_block_proc_col)
!MK! the blacs proc nr should be just the mpi nr
!MK                if (blacs_block_proc_nr == sm_block_nr) then
!MK                   if (blacs_block_proc_nr == global_env%mepos) then
!MK                      start_row_of_b= first_col(sm_block_nr)-b_block_row*&
!MK                           source_m%nrow_block
!MK                      end_row_of_b= min(source_m%nrow_block,&
!MK                           last_col(sm_block_nr)-b_block_row*&
!MK                           source_m%nrow_block)
!MK                      start_row_of_sm=max(0,-start_row_of_sm)
!MK                      start_row_of_b=max(0,start_row_of_sm)
!MK                      end_row_of_sm=start_row_of_sm+end_row_of_b-&
!MK                           start_row_of_b
!MK! to do
!MK!source_m%p(blacs_block_proc_row,blacs_block_proc_col)&
!MK!     %block(start_row_of_b:end_row_of_b,&
!MK!     start_col_of_b,end_col_of_b)
!MK
!MK                   end if
!MK                end if
!MK             end do
!MK          end do
!MK       end do
!MK
#else
       DO sm_block_nr=1,nblock_col
          CALL get_block_node(target_ptr,sm_block_nr,sm_block_nr,&
               BLOCK=block_val)
          block_val=block_val+source_m%local_data &
               (first_row(sm_block_nr):last_row(sm_block_nr),&
               first_col(sm_block_nr):last_col(sm_block_nr))
       END DO
#endif
    END IF
    CALL timestop(0.0_wp,handle)
  END SUBROUTINE add_blacs_to_block_diag_sm
!***************************************************************************

!***************************************************************************
!
! this is a jacobi-davidson method.
! diagonalises matrix_a, using matrix_v as an initial guess
! tol is the final tol for the eigenvalues
!
! the method expands the subspace by selecting a bunch of vectors
! and computing (A-shift*I)^-1*V  where shift is fixed somewhere in the spectrum
!
! for the 32 water, doesn't outperform dsyevx if more than 1/6 of the vectors are needed
!
! hence has not been parallellised. Kept here for later use. might need some
! cleaning.
! additionally, the method used to expand the subspace might be not optimal,
! experimenting (or reading literature) might help.
! (mai 2002 JVdV)
!***************************************************************************

  SUBROUTINE blacs_jacobi_davidson(matrix_a,matrix_v,evals,neig,tol, &
                                                 matrix_work,globenv)

! **************************************************************************

    TYPE(cp_full_matrix_type), POINTER          :: matrix_a, matrix_v, matrix_work
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), DIMENSION(:)                    :: evals
    INTEGER, INTENT(IN)                       :: neig
    REAL(wp), INTENT(IN)                      :: tol


!   *** Local variables ***

    INTEGER :: mypcol,myprow,npcol,nprow,source,n
    INTEGER :: nrow_block,ncol_block,ncol

    INTEGER, DIMENSION(9) :: desca,descv,descwork,descvecbuf,deschh,descvh,descvv
    INTEGER, DIMENSION(:), POINTER :: ipivot,iwork
    INTEGER :: lipivot, istat, info, i

    TYPE(cp_full_matrix_type), POINTER          ::  vecbuffer,matrix_overlap_hh
    TYPE(cp_full_matrix_type), POINTER          ::  matrix_overlap_vh,matrix_overlap_vv

    REAL(wp), DIMENSION(:,:), POINTER :: a,v,work,vecbuf,overlap_hh,overlap_vh
    REAL(wp), DIMENSION(:,:), POINTER :: overlap_vv
    REAL(wp), DIMENSION(:), POINTER :: evwork,evalsl,evalslold
    REAL(wp) :: rone,rzero,dum,dum2,ddot,dnrm2,shift,gap,gapnow
    INTEGER  :: nhigh,nvec,nlow,levwork,liwork,iter,ksteps,kblock,iter_getrs
    INTEGER  :: handle,handle2,handle3,this_one,count,j,igap,idamax,split
    INTEGER  :: context
    INTEGER, DIMENSION(:), POINTER :: to_be_used,to_be_used_old
    REAL(wp) :: time_init_start,time_init_end,time_init_total
    REAL(wp) :: time_full_start,time_full_end,time_full_total
    REAL(wp) :: time_diag_start,time_diag_end,time_diag_total
    REAL(wp) :: time_getrf_start,time_getrf_end,time_getrf_total
    REAL(wp) :: time_av_start,time_av_end,time_av_total
    REAL(wp) :: time_extend_start,time_extend_end,time_extend_total
    REAL(wp) :: time_getrs_start,time_getrs_end,time_getrs_total
    REAL(wp) :: time_ortho_start,time_ortho_end,time_ortho_total
    REAL(wp) :: time_suba_start,time_suba_end,time_suba_total
    REAL(wp) :: time_rot_start,time_rot_end,time_rot_total
    INTEGER :: count_av,count_getrs
    TYPE(cp_blacs_env_type), POINTER :: ctxt

!   ---------------------------------------------------------------------------
    CALL timeset("blacs_jacobi_davidson","I","",handle)

    context = matrix_a%matrix_struct%context%group

    time_init_total=0.0_wp
    time_full_total=0.0_wp
    time_diag_total=0.0_wp
    time_getrf_total=0.0_wp
    time_av_total=0.0_wp
    time_extend_total=0.0_wp
    time_getrs_total=0.0_wp
    time_ortho_total=0.0_wp
    time_suba_total=0.0_wp
    time_rot_total=0.0_wp

    count_av=0
    count_getrs=0

    time_init_start=m_cputime()
    time_full_start=m_cputime()

    rone=1.0_wp
    rzero=0.0_wp

    source = matrix_a%matrix_struct%para_env%source
    n = matrix_a%matrix_struct%nrow_global ! matrix dimension
    ksteps=5
    nhigh=MAX(ksteps*((neig/3)/ksteps),ksteps)  ! the upper part of the spectrum
    kblock=nhigh/ksteps
    nlow=neig-nhigh          ! the lower part of the spectrum
    nvec=neig+nhigh          ! the number of vectors used is subspace diag
    IF (nhigh.LT.1 .OR. nlow.LT.1 ) CALL stop_program("jacobi_davidson","nhigh, nlow, too small")
    IF (nhigh.GE.neig/2 ) CALL stop_program("jacobi_davidson","nhigh, too large")

    levwork=1+6*nvec+2*nvec**2
    liwork=5*nvec+3
    ALLOCATE(evwork(levwork))
    ALLOCATE(iwork(liwork))
    ALLOCATE(evalsl(nvec))
    ALLOCATE(evalslold(nvec))
    ALLOCATE(to_be_used(nhigh))
    ALLOCATE(to_be_used_old(nhigh))
    evalslold=0.0_wp

    CALL cp_fm_get_info(matrix=matrix_a,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block)
    ncol = 2*nvec
    ! some buffer space
    CALL cp_fm_create2(new_matrix=vecbuffer,&
                               nrow_global=n,&
                               ncol_global=ncol,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="VEC BUFFER",&
                               globenv=globenv)
    CALL cp_fm_create2(new_matrix=matrix_overlap_vh,&
                               nrow_global=neig,&
                               ncol_global=nhigh,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="overlap_vh",&
                               globenv=globenv)
    CALL cp_fm_create2(new_matrix=matrix_overlap_hh,&
                               nrow_global=nhigh,&
                               ncol_global=nhigh,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="overlap_hh",&
                               globenv=globenv)
    CALL cp_fm_create2(new_matrix=matrix_overlap_vv,&
                               nrow_global=nvec,&
                               ncol_global=nvec,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="overlap_vv",&
                               globenv=globenv)

! put the right pointers in place
    
    ctxt => matrix_a%matrix_struct%context
    myprow=ctxt%mepos(1)
    mypcol=ctxt%mepos(2)
    nprow=ctxt%num_pe(1)
    npcol=ctxt%num_pe(2)

#if defined(__SCALAPACK)
    CALL stop_program("blacs_jacobi_davidson","is not parallel")

    a => matrix_a%local_data
    desca(:) = matrix_a%matrix_struct%descriptor(:)

    v => matrix_v%local_data
    descv(:) = matrix_v%matrix_struct%descriptor(:)

    work => matrix_work%local_data
    descwork(:) = matrix_work%matrix_struct%descriptor(:)

    vecbuf => vecbuffer%local_data
    descvecbuf(:) = vecbuffer%matrix_struct%descriptor(:)

    overlap_vh => matrix_overlap_vh%local_data
    descvh(:) = matrix_overlap_vh%matrix_struct%descriptor(:)

    overlap_hh => matrix_overlap_hh%local_data
    deschh(:) = matrix_overlap_hh%matrix_struct%descriptor(:)

    overlap_vv => matrix_overlap_vv%local_data
    descvv(:) = matrix_overlap_vv%matrix_struct%descriptor(:)
#else
    a => matrix_a%local_data
    v => matrix_v%local_data
    work => matrix_work%local_data
    vecbuf => vecbuffer%local_data
    overlap_vh => matrix_overlap_vh%local_data
    overlap_hh => matrix_overlap_hh%local_data
    overlap_vv => matrix_overlap_vv%local_data
#endif

! ****************************************************************************
! make the input symmetrical
! ****************************************************************************

    CALL symmetrise_blacs_matrix(matrix_a,matrix_work)


! ****************************************************************************
! get the A*V for the occupied orbitals
! ****************************************************************************

#if defined(__SCALAPACK)

#else
    CALL dcopy(n*neig,v,1,vecbuf,1)
#endif

    time_av_start=m_cputime()
#if defined(__SCALAPACK)

#else
    count_av=count_av+neig
    CALL dgemm('N','N',n,neig,n,rone,a,n,vecbuf(1,1),n,rzero,vecbuf(1,nvec+1),n)
#endif
    time_av_end=m_cputime()
    time_av_total=time_av_total+(time_av_end-time_av_start)

! ****************************************************************************
! get the V^T*A*V for the occupied orbitals
! ****************************************************************************
    time_suba_start=m_cputime()
#if defined(__SCALAPACK)

#else
    CALL DGEMM('T','N',neig,neig,n,rone,vecbuf(1,1),n,vecbuf(1,nvec+1),n,rzero,overlap_vv,nvec)
#endif
    time_suba_end=m_cputime()
    time_suba_total=time_suba_total+(time_suba_end-time_suba_start)

! ****************************************************************************
! diag subblock
! ****************************************************************************
    time_diag_start=m_cputime()
#if defined(__SCALAPACK)

#else
    CALL DSYEVD('V','U',neig,overlap_vv,nvec,evalslold,evwork,levwork,iwork,liwork,info)
#endif
    IF (info.NE.0) CALL stop_program("jacobi_davidson","dsyevd overlap_vv 1")
    time_diag_end=m_cputime()
    time_diag_total=time_diag_total+(time_diag_end-time_diag_start)

! ****************************************************************************
! rotate vectors
! ****************************************************************************
    time_rot_start=m_cputime()
#if defined(__SCALAPACK)

#else
    CALL DGEMM('N','N',n,neig,neig,rone,vecbuf(1,     1),n,overlap_vv,nvec,rzero,v,n)
    CALL DCOPY(neig*n,v,1,vecbuf(1,1),1)
#endif
    time_rot_end=m_cputime()
    time_rot_total=time_rot_total+(time_rot_end-time_rot_start)

    ! start with the upper part of the spectrum
    DO i=1,nhigh
        to_be_used(i)=neig-nhigh+i
    ENDDO
    to_be_used_old=to_be_used

! ****************************************************************************
! we LU decompose the matrix a-shift*I
! ****************************************************************************
    time_getrf_start=m_cputime()
    CALL copy_blacs_to_blacs_matrix(matrix_a,matrix_work)
    lipivot = MAX(n,128) ! somehow too long ...
    ALLOCATE (ipivot(lipivot),stat=istat)
    IF (istat.NE.0) CALL stop_memory("blacs_jacobi_davidson","ipivot")
    ! find some gap in the lower half of the spectrum
    gap=0.0_wp
    igap=-1
    DO i=1,neig/2
       gapnow=evalslold(i+1)-evalslold(i)
       IF (gapnow.GE.gap) THEN
          gap=gapnow
          igap=i
       ENDIF
    ENDDO
    shift=0.5_wp*(evalslold(igap+1)+evalslold(igap))
    ! shift=1.1*evalslold(1)
#if defined(__SCALAPACK)
    CALL pdgetrf(n,n,work,1,1,descwork,ipivot,info)
#else
    ! put the singularity in the middle of the spectrum
    DO i=1,n
       work(i,i)=work(i,i)-shift
    ENDDO
    CALL dgetrf(n,n,work,n,ipivot,info)
#endif
    IF (info.NE.0) CALL stop_program("blacs_jacobi_davidson","problems in getrf")
    time_getrf_end=m_cputime()
    time_getrf_total=time_getrf_total+(time_getrf_end-time_getrf_start)

    time_init_end=m_cputime()
    time_init_total=time_init_total+(time_init_end-time_init_start)

! ****************************************************************************
! start main iteration
! ****************************************************************************
    DO iter=1,30

       ! ************************************************************************
       !  extension of the basis
       ! ************************************************************************
       time_extend_start=m_cputime()
       ! copy the needed vectors in place
#if defined(__SCALAPACK)
#else
       DO i=1,nhigh
          CALL DCOPY(n,vecbuf(1,to_be_used(i)),1,vecbuf(1,neig+i),1)
       ENDDO
#endif
       ! *********************************************************************
       ! get (LU)^-1 * vectors
       ! *********************************************************************
       time_getrs_start=m_cputime()
#if defined(__SCALAPACK)

#else
          DO i=1,ksteps-1
            CALL DAXPY(n*kblock,rone,vecbuf(1,neig+i*kblock+1),1,vecbuf(1,neig+1),1)
          ENDDO
          CALL DGETRS('N',n,kblock,work,n,ipivot,vecbuf(1,neig+1),n,info)
          DO i=1,ksteps-1
            CALL DCOPY(n*kblock,vecbuf(1,neig+1),1,vecbuf(1,neig+1+i*kblock),1)
            CALL DGETRS('N',n,kblock,work,n,ipivot,vecbuf(1,neig+1),n,info)
          ENDDO
          count_getrs=count_getrs+nhigh
#endif
          IF (info.NE.0) CALL stop_program("jacobi_davidson","Dgetrs")
       time_getrs_end=m_cputime()
       time_getrs_total=time_getrs_total+(time_getrs_end-time_getrs_start)


       ! *********************************************************************
       ! keep the new vectors orthogonal, is numerically rather hard, it seems
       ! *********************************************************************
#if defined(__SCALAPACK)

#else
       time_ortho_start=m_cputime()

       ! ortho to the existing basis vectors
       CALL DGEMM('T','N',neig,nhigh,n,rone,vecbuf(1,1),n,vecbuf(1,neig+1),n,rzero,overlap_vh,neig)
       CALL DGEMM('N','N',n,nhigh,neig,-rone,vecbuf(1,1),n,overlap_vh,neig,rone,vecbuf(1,neig+1),n)
       ! and ortho again within the new set
       CALL DSYRK('U','T',nhigh,n,rone,vecbuf(1,neig+1),n,rzero,overlap_hh,nhigh)
       CALL DPOTRF('U',nhigh,overlap_hh,nhigh,info)
       IF (info.NE.0) CALL stop_program("jacobi_davidson","DPOTRF overlap_hh")
       CALL DTRSM('R','U','N','N',n,nhigh,rone,overlap_hh,nhigh,vecbuf(1,neig+1),n)
       ! ortho to the existing basis vectors
       CALL DGEMM('T','N',neig,nhigh,n,rone,vecbuf(1,1),n,vecbuf(1,neig+1),n,rzero,overlap_vh,neig)
       CALL DGEMM('N','N',n,nhigh,neig,-rone,vecbuf(1,1),n,overlap_vh,neig,rone,vecbuf(1,neig+1),n)
       ! and ortho again within the new set
       CALL DSYRK('U','T',nhigh,n,rone,vecbuf(1,neig+1),n,rzero,overlap_hh,nhigh)
       CALL DPOTRF('U',nhigh,overlap_hh,nhigh,info)
       IF (info.NE.0) CALL stop_program("jacobi_davidson","DPOTRF overlap_hh")
       CALL DTRSM('R','U','N','N',n,nhigh,rone,overlap_hh,nhigh,vecbuf(1,neig+1),n)

       time_ortho_end=m_cputime()
       time_ortho_total=time_ortho_total+(time_ortho_end-time_ortho_start)
#endif

       ! *********************************************************************
       ! compute A*vectors for the new set
       ! *********************************************************************
       time_av_start=m_cputime()
#if defined(__SCALAPACK)

#else
       CALL dgemm('N','N',n,nhigh,n,rone,a,n,vecbuf(1,neig+1),n,rzero,vecbuf(1,nvec+neig+1),n)
       count_av=count_av+nhigh
#endif
       time_av_end=m_cputime()
       time_av_total=time_av_total+(time_av_end-time_av_start)

       time_extend_end=m_cputime()
       time_extend_total=time_extend_total+(time_extend_end-time_extend_start)

       ! *********************************************************************
       ! get v^T*(A*v)
       ! *********************************************************************
       time_suba_start=m_cputime()
#if defined(__SCALAPACK)

#else
       CALL DCOPY(nvec*nvec,rzero,0,overlap_vv,1)
       DO i=1,neig
          overlap_vv(i,i)=evalslold(i)
       ENDDO
       CALL DGEMM('T','N',nvec,nhigh,n,rone,vecbuf(1,1),n,vecbuf(1,nvec+neig+1),n, &
                                                    rzero,overlap_vv(1,neig+1),nvec)
       !CALL DGEMM('T','N',nvec,nvec,n,rone,vecbuf(1,1),n,vecbuf(1,nvec+1),n, &
       !                                                  rzero,overlap_vv,nvec)
#endif
       time_suba_end=m_cputime()
       time_suba_total=time_suba_total+(time_suba_end-time_suba_start)

       ! *********************************************************************
       ! diagonalise the matrix v^T*(A*v)
       ! *********************************************************************
       time_diag_start=m_cputime()
#if defined(__SCALAPACK)

#else
       CALL DSYEVD('V','U',nvec,overlap_vv,nvec,evalsl,evwork,levwork,iwork,liwork,info)
#endif
       IF (info.NE.0) CALL stop_program("jacobi_davidson","dsyevd overlap_vv")
       time_diag_end=m_cputime()
       time_diag_total=time_diag_total+(time_diag_end-time_diag_start)

       ! *********************************************************************
       ! here we select the eigenvectors with the largest change
       ! they will be used in the next iteration
       ! *********************************************************************
       ! write(6,*) "evalsl",NINT(-log(abs(evalsl(1:neig)-evalslold(1:neig))+1E-20)/log(10.0_wp))
       evwork(1:neig)=ABS(evalsl(1:neig)-evalslold(1:neig))
       CALL sort(evwork,neig,iwork)
       WRITE(6,*) "changes",iwork(neig),evwork(neig)
       ! use split of the part that changed most, half new
       count=neig
       split=nhigh/4
       DO i=1,split
          to_be_used(i)=iwork(count)
          count=count-1
       ENDDO
       DO i=split+1,nhigh
          DO
           this_one=iwork(count)
           count=count-1
           DO j=1,neig
             IF (this_one.EQ.to_be_used_old(j)) this_one=0
           ENDDO
           IF (this_one.NE.0) EXIT
          ENDDO
          to_be_used(i)=this_one
       ENDDO
       to_be_used_old=to_be_used
       evalslold=evalsl

       ! *********************************************************************
       ! rotate V, and AV so that the lowest neig are first
       ! *********************************************************************
       time_rot_start=m_cputime()
#if defined(__SCALAPACK)

#else
       CALL DGEMM('N','N',n,neig,nvec,rone,vecbuf(1,     1),n,overlap_vv,nvec,rzero,v,n)
       CALL DCOPY(neig*n,v,1,vecbuf(1,1),1)
#endif
       time_rot_end=m_cputime()
       time_rot_total=time_rot_total+(time_rot_end-time_rot_start)

       IF (evwork(neig).LE.tol) EXIT

    ENDDO

    ! *********************************************************************
    ! make the final vectors exactly ortho
    ! *********************************************************************
    time_ortho_start=m_cputime()
#if defined(__SCALAPACK)

#else
      CALL DSYRK('U','T',neig,n,rone,v(1,1),n,rzero,overlap_vv,nvec)
      CALL DPOTRF('U',neig,overlap_vv,nvec,info)
      ! in this case the vectors are exactly degenerate
      IF (info.NE.0) CALL stop_program("jacobi_davidson","DPOTRF final wave functions")
      CALL DTRSM('R','U','N','N',n,neig,rone,overlap_vv,nvec,v(1,1),n)
#endif
    time_ortho_end=m_cputime()
    time_ortho_total=time_ortho_total+(time_ortho_end-time_ortho_start)

    evals(1:neig)=evalsl(1:neig)

    CALL cp_fm_release(vecbuffer)
    CALL cp_fm_release(matrix_overlap_vh)
    CALL cp_fm_release(matrix_overlap_hh)
    CALL cp_fm_release(matrix_overlap_vv)
    DEALLOCATE(ipivot)
    DEALLOCATE(evwork)
    DEALLOCATE(iwork)
    DEALLOCATE(evalsl)
    DEALLOCATE(evalslold)
    DEALLOCATE(to_be_used)
    DEALLOCATE(to_be_used_old)

    time_full_end=m_cputime()
    time_full_total=time_full_total+(time_full_end-time_full_start)
    WRITE(6,*) "time full    : ", time_full_total
    WRITE(6,*) "time init    : ", time_init_total
    WRITE(6,*) "time extend  : ", time_extend_total
    WRITE(6,*) "---------------"
    WRITE(6,*) "time getrf   : ", time_getrf_total
    WRITE(6,*) "time getrs   : ", time_getrs_total
    WRITE(6,*) "time av      : ", time_av_total
    WRITE(6,*) "time diag    : ", time_diag_total
    WRITE(6,*) "time rot     : ", time_rot_total
    WRITE(6,*) "time suba    : ", time_suba_total
    WRITE(6,*) "time ortho   : ", time_ortho_total
    WRITE(6,*) "---------------"
    WRITE(6,*) "# neig/nhigh/n ", neig, "/",nhigh,"/",n
    WRITE(6,*) "# A*V /t       ", count_av,"/",time_av_total/count_av
    WRITE(6,*) "# (L*U)^-1*V/t ", count_getrs,"/",time_getrs_total/count_getrs

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_jacobi_davidson

! *****************************************************************************


!!****f* qs_blacs/blacs_scale_and_d [1.0] *
!!
!!   NAME
!!     blacs_scale_and_d
!!
!!   FUNCTION
!!     scales a balcs matrix and adds a scalar to the diagonal elements.
!!     Can be used to calculate (I-blacs_matrix) in place
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - blacs_matrix: the matrix to change
!!     - scale: how much to rescale the matrix (defaults to 1.0)
!!     - diagonal_add: how much is added to the diagonal elements of
!!       the matrix
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE blacs_scale_and_d(blacs_matrix,scale, diagonal_add,error)
    TYPE(cp_full_matrix_type), POINTER :: blacs_matrix
    TYPE(cp_para_env_type), POINTER :: para_env
    REAL(kind=wp), INTENT(in), OPTIONAL :: scale, diagonal_add
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error
    
    INTEGER :: i,j,nrow_local, ncol_local,myprow,mypcol
    INTEGER, DIMENSION(:), POINTER :: row_indices, col_indices
    REAL(kind=wp) :: sc, diag
    REAL(kind=wp), DIMENSION(:,:), POINTER :: blacs_block
    TYPE(cp_blacs_env_type), POINTER :: context

    sc=1.0_wp;  diag=0.0_wp
    IF (PRESENT(scale)) sc=scale
    IF (PRESENT(diagonal_add)) diag=diagonal_add
    
    context => blacs_matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)

    blacs_block => blacs_matrix%local_data

    IF (diag==0.0_wp) THEN
      blacs_block=sc*blacs_block
    ELSE
      CALL cp_fm_get_info(matrix=blacs_matrix, &
          nrow_local =nrow_local ,ncol_local =ncol_local , &
          row_indices=row_indices,col_indices=col_indices)
      IF (sc==0.0) THEN
        CALL dcopy(SIZE(blacs_block,1)*SIZE(blacs_block,2),0.0_wp,0,blacs_block,1)
        DO j=1,ncol_local
          DO i=1,nrow_local
            IF (col_indices(j)== row_indices(i)) THEN
              blacs_block(i,j)=diag
            END IF
          END DO
        END DO
      ELSE
        DO j=1,ncol_local
          DO i=1,nrow_local
            IF (col_indices(j)== row_indices(i)) THEN
              blacs_block(i,j)=diag+sc*blacs_block(i,j)
            ELSE
              blacs_block(i,j)=sc*blacs_block(i,j)
            END IF
          END DO
        END DO
      END IF
      
    END IF
END SUBROUTINE blacs_scale_and_d

!   ***************************************************************************
!   ***************************************************************************
!   ***************************************************************************
!   ***************************************************************************

  SUBROUTINE blacs_jacobi(matrix,eigenvectors,eigenvalues,epsi,Anz_max,&
                          thresh,blocksize)

!   Purpose: Calculates Eigenvalues und Eigenvectors from full symmetric matrix

!   History: - Creation (16.09.2002, Martin Fengler) Has its origin in blacs_syevx

!   Remark: The Eigenvalues and Eigenvectors are NOT sorted! In our purpose it
!           is not necessary to sort them since they are already at their
!           correct posistions after a couple of SCF-Steps.
!           So this routine should only be used after some general eigenvalue
!           solvers worked on the same problem

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER        :: eigenvectors,matrix
    INTEGER, INTENT(IN),optional              :: blocksize   ! Gibt an wieviele Eigenwerte/vektoren berechnet werden solln
    INTEGER, INTENT(IN)                       :: Anz_max     ! Gibt die maximal Anzahl an Jacobi-Iterationen an
    REAL(wp), INTENT(IN)                      :: epsi,thresh ! epsi: Abbruchkriterium, wenn epsi kleiner als die Quadratsumme
                                                             ! der Ausserdiagonalwerte ist
                                                             ! thresh: Ein Element wird nur rotiert, wenn es groesser als thresh ist.         
    REAL(wp), DIMENSION(:), INTENT(OUT)       :: eigenvalues ! Enthaelt die Eigenwerte 1:blocksize (vom kleinsten bis zum groessten)
                                                             ! Da ab einem gewissen SCF-SChritt nicht mehr umsortiert werden muss 
!   *** Local variables ***
    INTEGER :: handle
    REAL(wp), DIMENSION(:,:), POINTER   :: A,EV         ! EV: Eigenvektoren

!   ***** LOKALE VARIABLEN *******
    REAL(wp) :: tan_theta,tau,c,s,J_pp,J_pq,J_qp,J_qq,A_off
    integer  :: I,J,Anz,q,p,bs,N
    REAL(wp), dimension (:),allocatable :: c_ip,c_iq

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_jacobi","I","",handle)
        
    N = matrix%matrix_struct%nrow_global ! Groesse der Matrix A, die diagonisiert werden soll
  
  
    allocate(c_ip(N),c_iq(N))   ! Speicher fuer den lokalen Eigenwertvektor

    eigenvalues(:) = 0.0_wp      ! Contains at the end the eigenvalues (more precise: The diagonal of A)
    A => matrix%local_data       ! Contains the Matrix to be worked on
    EV => eigenvectors%local_data ! Contains the eigenvectors up to blocksize: Rest ist Muell

!   *** START Diagonalising matrix ***
!   *** START Diagonalising matrix ***

       ! Arguments checken
       if(.not. present(blocksize)) then
         bs = N
       else      
         bs = blocksize
       endif      
      
      
       tan_theta = 0.0_wp;
       tau       = 0.0_wp;
       
       ! Q zur Einheitsmatrix machen
       EV = 0.0_wp ;
       forall (I=1:N)
         EV(I,I)=1.0_wp
       endforall
       
       ! Summe der Off-Diagonalwerte bis BLocksize ausrechnen
       A_off = 0.0_wp;
       do J=2,bs
         do I=1,(J-1)
           A_off = A_off + A(I,J)*A(I,J);
         enddo 
       enddo      
       A_off=2.0_wp * A_off;
      
      ! Anzahl der Durchgaenge
       Anz=1;
 
      do 
         if (A_off < epsi ) then
           exit
         endif  
    
   
         do q=2,bs
          do p=1,q-1
    
  
          if(abs(A(p,q))>thresh) then
 
           ! Durch das Wegrotieren erniedirgt sich die Quadratsumme der 
           ! Ausserdiagonalwerte:
          A_off = A_off - 2.0_wp*A(p,q)*A(p,q);
   
          tau = (A(q,q)-A(p,p))/(2.0_wp*A(p,q));
    
          tan_theta = sign(1.0_wp,tau)/(abs(tau)+sqrt(1.0_wp+tau*tau));
    
          ! cos theta
          c = 1.0_wp/sqrt(1.0_wp+tan_theta*tan_theta);
          s = tan_theta*c;   
     
          ! Schnelle Implementierung von: J' * A * J
 
          ! Erst wird von links dranmultiplitziert
          ! J'*A , weil die Symmetrie mir die Strides rettet

          J_pp=c;  J_pq=s;
          J_qp=-s; J_qq=c;          


          c_ip(1:bs)=A(1:bs,p)*J_pp + A(1:bs,q)*J_qp;
          c_iq(1:bs)=A(1:bs,p)*J_pq + A(1:bs,q)*J_qq;
     
          A(1:bs,p)=c_ip(1:bs);
          A(1:bs,q)=c_iq(1:bs);

    
          ! Jetzt dann von rechts
          ! A*J
             
          A(p,p) = J_pp*A(p,p) + J_qp*A(q,p);
          A(q,q) = J_pq*A(p,q) + J_qq*A(q,q);          
          
          A(p,q) = 0.0_wp;
          A(q,p) = 0.0_wp;        
          
          ! Die 4 Elemente A(q,q),A(p,q),A(q,p),A(q,q) werden unnoetig mit transferiert...soll mal latte sein
           A(p,1:bs)=A(1:bs,p)
           A(q,1:bs)=A(1:bs,q)

          
          !A(p,1:p-1)  = A(1:p-1,p)
          !A(p,p+1:q-1)= A(p+1:q-1,p)
          !A(p,q+1:bs) = A(q+1:bs,p)
          
          !A(q,1:p-1)  = A(1:p-1,q)
          !A(q,p+1:q-1)= A(p+1:q-1,q)
          !A(q,q+1:bs) = A(q+1:bs,q)

             
 
          ! Und jetzt noch die Eigenvektoren produzieren:
          ! Erst wird von rechts dranmultiplitziert
          ! Q * J
          if (Anz==1) then       
            c_ip(1:q)=EV(1:q,p)*J_pp + EV(1:q,q)*J_qp;
            c_iq(1:q)=EV(1:q,p)*J_pq + EV(1:q,q)*J_qq;
     
            EV(1:q,p)=c_ip;
            EV(1:q,q)=c_iq;
          else
            c_ip=EV(:,p)*J_pp + EV(:,q)*J_qp;
            c_iq=EV(:,p)*J_pq + EV(:,q)*J_qq;
     
            EV(:,p)=c_ip;
            EV(:,q)=c_iq;
          endif  

        endif


        if(A_off<epsi) then
          exit
        endif

 
      enddo  ! for-loop Innere
     
      if(A_off<epsi) then
       exit
      endif

    
    enddo   ! for-loop 

    Anz=Anz+1;
  
    if(Anz>Anz_max) then
      exit
    endif
  
   enddo  ! END-WHILE

  if(A_off<=epsi) then
    write(*,*) 'Die Quadrat-Summe der Ausserdiagnonalelemente ist mit unter epsi: ',A_off
  endif
 
 ! Abschliessend noch sortieren der Eigenwerte und Eigenvektoren
 ! nach der Groesse

 write(*,*) 'Anzahl der Zyklischen Durchlaeufe: ',(Anz-1)
   
 forall (I=1:blocksize)
    eigenvalues(I)=A(I,I)
 endforall


!   *** END Diagonalising matrix ***
!   *** END Diagonalising matrix ***


!   *** Release work storage ***
    DEALLOCATE (c_ip,c_iq)
    
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_jacobi

! *****************************************************************************

!   ***************************************************************************
!   ***************************************************************************
!   ***************************************************************************
!   ***************************************************************************

  SUBROUTINE blacs_block_jacobi_scaled(matrix,eigenvectors,eigenvalues,thresh,start_sec_block)
    implicit none
!
!   Purpose: Calculates Block diagonalizazion from full symmetric matrix
!   History: - Creation (07.10.2002, Martin Fengler) Has its origin in blacs_syevx
!   This Routines scales on-the-fly the matrix elements to detect significant influences on the eigenvalues
!   IT DOES ONLY ONE SWEEP!

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER        :: eigenvectors,matrix
    INTEGER, INTENT(IN)                       :: start_sec_block  ! Gibt an wo der zweite Block BEGINNT
    REAL(wp), INTENT(IN)                      :: thresh ! thresh: Ein SKALIERTES Element wird nur rotiert, wenn es grsser als thresh ist.         
    REAL(wp), DIMENSION(:), INTENT(IN)        :: eigenvalues ! Enthaelt die Eigenwerte 1:N(vom kleinsten bis zum grssten)
                                                             ! Da ab einem gewissen SCF-SChritt nicht mehr umsortiert werden muss 

!   *** Local variables ***
    INTEGER :: handle
    REAL(wp), DIMENSION(:,:), POINTER   :: A,EV         ! EV: Eigenvektoren

!   ***** LOKALE VARIABLEN *******
    REAL(wp) :: tan_theta,tau,c,s,k,max_thresh,new_thresh, temp
    integer  :: I,J,q,p,bs,N
    REAL(wp), dimension (:),allocatable :: c_ip,c_iq
      


!   ---------------------------------------------------------------------------

    CALL timeset("blacs_block_jacobi_scaled","I","",handle)
        
    N = matrix%matrix_struct%nrow_global ! Groesse der Matrix A, die diagonisiert werden soll
  
    allocate(c_ip(N),c_iq(N))   ! Speicher fuer den lokalen Eigenwertvektor

    A => matrix%local_data       ! Contains the Matrix to be worked on
    EV => eigenvectors%local_data ! Contains the eigenvectors up to blocksize: Rest ist Muell


!   *** START Diagonalising matrix ***
!   *** START Diagonalising matrix ***
      
       tan_theta = 0.0_wp;
       tau       = 0.0_wp;
       
       ! Q zur Einheitsmatrix machen
       EV = 0.0_wp ;
       forall (I=1:N)
         EV(I,I)=1.0_wp
       endforall

    k=thresh
    max_thresh=0.0_wp
    do q=start_sec_block,N
      do p=1,(start_sec_block-1)
       temp=abs(1.0_wp/(eigenvalues(p)-eigenvalues(q))*A(p,q))
       if(temp>max_thresh) then
         max_thresh=temp
       endif     
      enddo
    enddo  

    new_thresh=k*max_thresh;
       
      
    do q=start_sec_block,N
      do p=1,(start_sec_block-1)
    
          temp=abs(1.0_wp/(eigenvalues(p)-eigenvalues(q))*A(p,q))
          if(temp>new_thresh .or. temp>0.00001_wp) then
         
           ! Durch das Wegrotieren erniedirgt sich die Quadratsumme der 
           ! Ausserdiagonalwerte:
  
          tau = (A(q,q)-A(p,p))/(2.0_wp*A(p,q));
    
          tan_theta = sign(1.0_wp,tau)/(abs(tau)+sqrt(1.0_wp+tau*tau));
    
          ! cos theta
          c = 1.0_wp/sqrt(1.0_wp+tan_theta*tan_theta);
          s = tan_theta*c;   

          ! Und jetzt noch die Eigenvektoren produzieren:
          ! Q * J
            c_ip = c*EV(:,p) - s*EV(:,q);
            c_iq = s*EV(:,p) + c*EV(:,q);
     
            EV(:,p)=c_ip;
            EV(:,q)=c_iq;

         endif


      enddo  ! for-loop Innere
    enddo   ! for-loop 

!   *** END Diagonalising matrix ***
!   *** END Diagonalising matrix ***


!   *** Release work storage ***
    DEALLOCATE (c_ip,c_iq)
    
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_block_jacobi_scaled

! *****************************************************************************


  SUBROUTINE blacs_block_jacobi_classic(matrix,eigenvectors,thresh,&
                                        start_sec_block)
!
!   Purpose: Calculates Block diagonalizazion from full symmetric matrix
!   History: - Creation (07.10.2002, Martin Fengler) Has its origin in blacs_syevx
!   This Routines rotates only elements which are larger than a threshold "thresh"
!   IT DOES ONLY ONE SWEEP!

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER        :: eigenvectors,matrix
    INTEGER, INTENT(IN)                       :: start_sec_block  ! Gibt an wo der zweite Block BEGINNT
    REAL(wp), INTENT(IN)                      :: thresh ! thresh: Ein Element wird nur rotiert, wenn es groesser als thresh ist.         

!   *** Local variables ***
    INTEGER :: handle
    REAL(wp), DIMENSION(:,:), POINTER   :: A,EV         ! EV: Eigenvektoren

!   ***** LOKALE VARIABLEN *******

    REAL(wp) :: tan_theta,tau,c,s
    integer  :: I,J,q,p,bs,N
    REAL(wp), dimension (:),allocatable :: c_ip

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL  :: numroc,indxl2g
    integer :: myprow,mypcol,nprow,npcol,ictxt_loc,block_dim_row,block_dim_col,info,ev_row_block_size,iam,counter,status,ierr
    integer :: source,allgrp,mynumrows,ictxt,temp_row,temp_col,temp_nprow,temp_npcol,mype,npe,a_row,a_col,a_nb,a_mb,iglob,jglob
    INTEGER, DIMENSION(9) :: desca,descz,desc_a_block, desc_ev_loc
    REAL(wp), allocatable :: A_loc(:,:), EV_loc(:,:),diagonale(:)
    TYPE(cp_blacs_env_type), POINTER    :: context
#endif

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_block_jacobi_classic","I","",handle)
        
#if defined(__SCALAPACK)
    
    ! 'Umgebung fuer die Kohn-Sham-Matrix und die Eigenvektoren uebernehmen!'
    
    context =>  matrix%matrix_struct%context
    source =    matrix%matrix_struct%para_env%source ! globenv%ionode
    allgrp =    matrix%matrix_struct%para_env%group ! globenv%ionode

    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    
    N = matrix%matrix_struct%nrow_global

    A => matrix%local_data
    desca(:) = matrix%matrix_struct%descriptor(:)
    EV => eigenvectors%local_data
    descz(:) = eigenvectors%matrix_struct%descriptor(:)
    
    allocate(diagonale(N),STAT=status)
    if(status > 0) then
      write (*,*) 'Allok-Error fuer die Diagonale!'
    endif  
    diagonale(:) = 0.0_wp

    ! Extrahiere Diagonale der Matrix 
    ! Dazu guckt jeder in seinem lokalen Block und schict dann diese Diagonalteile samt INdexe an den Master

    a_row = size(A,1)
    a_col = size(A,2)
    a_mb  = desca(5)
    a_nb  = desca(6)

    do J=1,a_col
      do I=1,a_row
        iglob = indxl2g(I,a_mb,myprow,0,nprow)
        jglob = indxl2g(J,a_nb,mypcol,0,npcol)
        if(iglob==jglob) diagonale(iglob) = a(i,j)
      enddo
    enddo                

    ! Hole die Diagonalteile zusammen...
    
     CALL mp_sum(diagonale,allgrp)

    ! Kopiere den Block, der wegrotiert werden soll zunaechst auf den Masterprozessor, und anschliessend
    ! per Broadcast an alle Prozis
    ! ACHTUNG start_sec_block sagt aus WO der ZWEITE Block STARTET!!!
    ! Der Block wird mitsamt dem OO-Block bearbeitet
    
    block_dim_row = start_sec_block-1
    block_dim_col = N
    allocate(A_loc(block_dim_row,block_dim_col))
    
    !call blacs_get(-1,0,ictxt)
    ictxt=0
    
    CALL blacs_pinfo(mype,npe)
    call blacs_get(ictxt,10,ictxt_loc)
    
    call blacs_gridinit(ictxt_loc,'Row-major',NPROW*NPCOL,1)

    call descinit(desc_a_block,block_dim_row,block_dim_col,block_dim_row,&
                  block_dim_col,0,0,ictxt_loc,block_dim_row,info)
    
    ! Jetzt sind eigentlich nur im Master-Prozess Daten reingekommen
    call pdgemr2d(block_dim_row,block_dim_col,A,1,1,desca,A_loc,1,1,desc_a_block,ictxt)

    CALL mp_bcast (A_loc,0, allgrp )
  
    ! Da nun jeder ueber den oberen Block verfuegt, koennen wir jetzt die Eigenvektoren so umsortieren, dass
    ! man ein NN*1 Prozessgrid hat, und somit jeder Prozessor ueber ein Buendel von Zeilel verfuegt,
    ! die selbst unabhaengig modifizieren darf.

    ! Aufsetzen der Eigenvektorverteilung
    iam = mype
    ev_row_block_size = 32
    mynumrows = NUMROC(N,ev_row_block_size,iam,0,NPROW*NPCOL)
 
    allocate(EV_loc(mynumrows,N),c_ip(mynumrows))
        
    call descinit(desc_ev_loc,N,N,ev_row_block_size,N,0,0,ictxt_loc,&
                  mynumrows,info)

!   *** START Diagonalising matrix ***
!   *** START Diagonalising matrix ***
       
    ! 'Q zur Einheitsmatrix machen!'
    call blacs_gridinfo(ictxt_loc,nprow,npcol,myprow,mypcol)
    
       do J=1,N
         do I=1,mynumrows
           if(indxl2g(I,ev_row_block_size,myprow,0,nprow)==indxl2g(J,N,mypcol,0,npcol)) then
             EV_loc(I,J) = 1.0_wp ;
           else  
             EV_loc(I,J) = 0.0_wp ;
           endif
         enddo
       enddo                

        
    ! Eigentliche Blockdiagonalisierung
    do q=start_sec_block,N
      do p=1,(start_sec_block-1)
    
          if(abs(A_loc(p,q))>thresh) then         
  
          tau = (diagonale(q)-diagonale(p))/(2.0_wp*A_loc(p,q));
    
          tan_theta = sign(1.0_wp,tau)/(abs(tau)+sqrt(1.0_wp+tau*tau));
    
          ! cos theta
          c = 1.0_wp/sqrt(1.0_wp+tan_theta*tan_theta);
          s = tan_theta*c;   

          ! Und jetzt noch die Eigenvektoren produzieren:
          ! Q * J
          !  Verstaendliche Version (bevor die BLAS-Aufrufe sie ersetzt haben)
          !  c_ip = c*EV_loc(:,p) - s*EV_loc(:,q);
          !  c_iq = s*EV_loc(:,p) + c*EV_loc(:,q);
     
          !  EV(:,p)=c_ip;
          !  EV(:,q)=c_iq;
          
          call dcopy(mynumrows,EV_loc(1,q),1,c_ip,1)
          call dscal(mynumrows,-s,c_ip,1)
          call daxpy(mynumrows,c,EV_loc(1,p),1,c_ip,1)
          
          call dscal(mynumrows,c,EV_loc(1,q),1)
          call daxpy(mynumrows,s,EV_loc(1,p),1,EV_loc(1,q),1)
          
          call dcopy(mynumrows,c_ip,1,EV_loc(1,p),1)
          
         endif


      enddo  ! for-loop Innere
    enddo   ! for-loop 

!   *** END Diagonalising matrix ***
!   *** END Diagonalising matrix ***

    ! Nun muessen die Eigenvektoren wieder in die alte Verteilung zurueckverschickt werden.
    ! Verschicke EVs zurueck!'
    call pdgemr2d(N,N,EV_loc,1,1,desc_ev_loc,EV,1,1,descz,ictxt_loc)

    ! Speicher freigeben
    deallocate(A_loc,EV_loc,c_ip,diagonale)

    ! ...und tschuess!
    
#else

    N = matrix%matrix_struct%nrow_global ! Groesse der Matrix A, die bearbeitet werden soll
  
  
    allocate(c_ip(N))   ! Speicher fuer den lokalen Eigenwertvektor

    A => matrix%local_data       ! Contains the Matrix to be worked on
    EV => eigenvectors%local_data ! Contains the eigenvectors up to blocksize: Rest ist Muell

!   *** START Diagonalising matrix ***
!   *** START Diagonalising matrix ***
      
       tan_theta = 0.0_wp;
       tau       = 0.0_wp;
       
       ! Q zur Einheitsmatrix machen
       EV = 0.0_wp ;
       forall (I=1:N)
         EV(I,I)=1.0_wp
       endforall
       
          
    do q=start_sec_block,N
      do p=1,(start_sec_block-1)
    
          if(abs(A(p,q))>thresh) then         
  
          tau = (A(q,q)-A(p,p))/(2.0_wp*A(p,q));
    
          tan_theta = sign(1.0_wp,tau)/(abs(tau)+sqrt(1.0_wp+tau*tau));
    
          ! cos theta
          c = 1.0_wp/sqrt(1.0_wp+tan_theta*tan_theta);
          s = tan_theta*c;   

          ! Und jetzt noch die Eigenvektoren produzieren:
          ! Q * J
          !  Verstaendliche Version (bevor die BLAS-Aufrufe sie ersetzt haben)
          !  c_ip = c*EV(:,p) - s*EV(:,q);
          !  c_iq = s*EV(:,p) + c*EV(:,q);
     
          !  EV(:,p)=c_ip;
          !  EV(:,q)=c_iq;
          
          call dcopy(N,EV(1,q),1,c_ip,1)
          call dscal(N,-s,c_ip,1)
          call daxpy(N,c,EV(1,p),1,c_ip,1)
          
          call dscal(N,c,EV(1,q),1)
          call daxpy(N,s,EV(1,p),1,EV(1,q),1)
          
          call dcopy(N,c_ip,1,EV(1,p),1)
          
         endif


      enddo  ! for-loop Innere
    enddo   ! for-loop 

!   *** END Diagonalising matrix ***
!   *** END Diagonalising matrix ***


!   *** Release work storage ***
    DEALLOCATE (c_ip)
    
#endif    

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_block_jacobi_classic

! *****************************************************************************

END MODULE qs_blacs
