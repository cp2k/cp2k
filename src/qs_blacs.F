
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001, 2002, 2003  CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/qs_blacs [1.0] *
!!
!!   NAME
!!     qs_blacs
!!
!!   FUNCTION
!!     most of the routines in this module should not be used anymore (use cp_fm_* )
!!     the remaining ones should also go one day
!!
!!   AUTHOR
!!     Matthias Krack (22.05.2001)
!!
!!   MODIFICATION HISTORY
!!     1) Changed the size work array for syevx (Joost VandeVondele 28.02.02)
!!     1) Introduced some cholesky guys         (Joost VandeVondele    05.02)
!!     2) Introduced a symmetric matrix x distributed vector multiply ( JVdV)
!!     3) start to remove the file
!!
!!   SOURCE
!******************************************************************************

MODULE qs_blacs
! *****************************************************************************
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_unimplemented_error
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_syrk
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_equivalent,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_retain,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_init_random,&
                                             cp_fm_release,&
                                             cp_fm_retain,&
                                             cp_fm_type
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
  USE mathlib,                         ONLY: symmetrize_matrix
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_gather,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sync, &
                                             mp_alltoall
  USE parallel,                        ONLY: cp2k_is_parallel
  USE sparse_matrix_types,             ONLY: first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_type
  USE string_utilities,                ONLY: compress
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN='qs_blacs'

! *** Public subroutines ***
  PUBLIC :: copy_blacs_to_ionode_matrix,&
            blacs_replicated_copy,&
            blacs_column_copy,&
            write_blacs_matrix, &
       blacs_scale_and_d,&
       copy_vector_to_blacs, &
       copy_blacs_to_vector, &
       blacs_find_vector_distribution,&
       read_blacs_matrix, &
       blacs_daxpy,  &
       blacs_get_element

  PUBLIC :: blacs_add_to_element, &
            blacs_vector_symv, &
            symmetrise_blacs_matrix,&
            add_blacs_to_block_diag_sm,  &
            blacs_make_basis, &
            blacs_transpose,  &
            blacs_scale_matrix

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE blacs_daxpy(matrix_a,c1,alpha,c2,error)
    TYPE(cp_fm_type), POINTER       :: matrix_a
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

    TYPE(cp_fm_type), POINTER          :: matrix
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

    TYPE(cp_fm_type), POINTER          :: matrixa
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
    TYPE(cp_fm_type), POINTER       :: matrix_a
    INTEGER, INTENT(OUT)                   :: nloc

    TYPE(cp_blacs_env_type), POINTER       :: context
    INTEGER :: n,nloc_blacs

    n = matrix_a%matrix_struct%nrow_global
    context => matrix_a%matrix_struct%context

!!FM#if defined(__SCALAPACK)
    nloc_blacs = matrix_a%matrix_struct%nrow_locals(context%mepos(1))
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
    TYPE(cp_fm_type), POINTER       :: matrix_a
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

    nloc_blacs = matrix_a%matrix_struct%nrow_locals(context%mepos(1))
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
    IF (matrix_a%matrix_struct%nrow_locals(myprow).LT.use_col_nr) THEN
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
    TYPE(cp_fm_type), POINTER       :: matrix_a
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

    nloc_blacs = matrix_a%matrix_struct%nrow_locals(myprow)
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
    IF (matrix_a%matrix_struct%ncol_locals(mypcol).LT.use_col_nr) THEN
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

    TYPE(cp_fm_type), POINTER          :: matrix_a,matrix_buf
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

    TYPE(cp_fm_type), POINTER                   :: vmatrix
    TYPE(cp_fm_type), POINTER, OPTIONAL         :: pmatrix,work,ortho
    INTEGER, INTENT(IN)                                :: ncol
    CHARACTER ( LEN = * ), INTENT ( IN ),OPTIONAL      :: otype

!   *** Local variables ***

    INTEGER :: handle,n,info
    REAL(wp), DIMENSION(:,:), POINTER :: a,b,vv,w,p,o
    REAL(wp) :: rone,rzero
    INTEGER :: nprow,npcol,myprow,mypcol,nrow_block,ncol_block,i
    INTEGER, DIMENSION(9) :: desca,descb,descvv,descw,descp,desco
    TYPE(cp_fm_type), POINTER           :: overlap_vv
    LOGICAL :: found_type,use_pmatrix,use_cholesky,use_sv,use_simple
    TYPE(cp_blacs_env_type), POINTER :: context
    TYPE(cp_fm_struct_type), POINTER :: fm_struct_tmp

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_make_basis","I","",handle)
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
    rone=1.0_wp
    rzero=0.0_wp
    n   = vmatrix%matrix_struct%nrow_global
    CALL cp_fm_get_info(matrix=vmatrix,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block)

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

    CALL timestop(0.0_wp,handle)

    END SUBROUTINE blacs_make_basis

! *****************************************************************************

  SUBROUTINE blacs_get_element(matrix,irow_global,icol_global,alpha)

!   Purpose: Get the BLACS matrix element (irow_global,icol_global).

!   History: - Creation (22.01.2002, Matthias Krack)
!   ***************************************************************************

    TYPE(cp_fm_type), POINTER          :: matrix
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
    TYPE(cp_fm_type) , POINTER              :: matrix_v
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

    ! CALL timeset("blacs_replicated_copy","I","",handle)

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

    ! CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_replicated_copy

  SUBROUTINE blacs_column_copy(msource,mtarget,ncol,source_start,target_start)
    TYPE(cp_fm_type), POINTER          :: msource,mtarget
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

! *****************************************************************************

  SUBROUTINE copy_blacs_to_ionode_matrix(blacs_matrix,full_matrix,error)

!   Purpose: Copy a BLACS matrix to a full matrix on the ionode (does it work in parallel ? ? ? ? ? ? )

!   History: - Creation (18.06.2001, Matthias Krack)
!   - 08.2002 changed for the transition to local_data, unoptimized
!     (but used for debugging only...) [fawzi]
!   ***************************************************************************

    TYPE(cp_fm_type), POINTER          :: blacs_matrix
    REAL(wp), DIMENSION(:,:), POINTER         :: full_matrix
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_blacs_to_ionode_matrix (MODULE qs_blacs)"

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

    CALL timeset("copy_blacs_to_ionode_matrix","I","",handle)

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

        nrow_local = blacs_matrix%matrix_struct%nrow_locals(iprow)
        ncol_local = blacs_matrix%matrix_struct%ncol_locals(ipcol)

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

  END SUBROUTINE copy_blacs_to_ionode_matrix

! *****************************************************************************

  SUBROUTINE read_blacs_matrix(matrix,lunit,globenv)

!   Purpose: Read a BLACS matrix from the logical unit number "lunit".

!   History: - Creation (19.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: lunit

!   *** Local variables ***

    INTEGER :: i,j,mypcol,myprow,ncol_local,nrow_local
    TYPE(cp_blacs_env_type), POINTER :: context
!   ---------------------------------------------------------------------------

    context => matrix%matrix_struct%context
    myprow=context%mepos(1)
    mypcol=context%mepos(2)

    nrow_local = matrix%matrix_struct%nrow_locals(myprow)
    ncol_local = matrix%matrix_struct%ncol_locals(mypcol)

    READ (UNIT=lunit) ((matrix%local_data(i,j),i=1,nrow_local),&
                                                           j=1,ncol_local)

  END SUBROUTINE read_blacs_matrix

! *****************************************************************************
  SUBROUTINE blacs_transpose(matrix,matrixt)
! matrixt=matrix^T ! assumes no symmetry
!***************************************************************************

    TYPE(cp_fm_type), POINTER          :: matrix,matrixt

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

    TYPE(cp_fm_type), POINTER          :: matrix,work

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

    nrow_local = matrix%matrix_struct%nrow_locals(myprow)
    ncol_local = matrix%matrix_struct%ncol_locals(mypcol)

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

    TYPE(cp_fm_type), POINTER          :: matrix
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

    nrow_local = matrix%matrix_struct%nrow_locals(myprow)
    ncol_local = matrix%matrix_struct%ncol_locals(mypcol)

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
!!   SYNOPSIS
!!     Subroutine add_blacs_to_block_diag_sm(source_m, target_m, para_env,&
!!         d_struct, error)
!!       Type(cp_fm_type), Pointer:: source_m
!!       Type(real_matrix_type), Intent (INOUT), Target:: target_m
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Integer, Dimension(:), Pointer, Optional:: d_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine add_blacs_to_block_diag_sm
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
    TYPE(cp_fm_type), POINTER :: source_m
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

!!****f* qs_blacs/blacs_scale_and_d [1.0] *
!!
!!   NAME
!!     blacs_scale_and_d
!!
!!   SYNOPSIS
!!     Subroutine blacs_scale_and_d(blacs_matrix, scale, diagonal_add,&
!!         error)
!!       Type(cp_fm_type), Pointer:: blacs_matrix
!!       Real(Kind=wp), Intent (IN), Optional:: scale, diagonal_add
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine blacs_scale_and_d
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
    TYPE(cp_fm_type), POINTER :: blacs_matrix
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

END MODULE qs_blacs
