
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
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
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
            blacs_column_copy,&
            write_blacs_matrix, &
            blacs_make_basis

! *****************************************************************************

CONTAINS

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

END MODULE qs_blacs
