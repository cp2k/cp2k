!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
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
!!
!!   SOURCE
!******************************************************************************

MODULE qs_blacs

! *****************************************************************************

  USE kinds, ONLY: int_size,&
                   wp => dp,&
                   wp_size => dp_size

  USE global_types,        ONLY: global_environment_type
  USE mathlib,             ONLY: symmetrize_matrix
  USE sparse_matrix_types, ONLY: first_block_node,&
                                 get_block_node,&
                                 get_matrix_info,&
                                 next_block_node,&
                                 real_block_node_type,&
                                 real_matrix_type
  USE message_passing,     ONLY: mp_bcast,&
                                 mp_max,&
                                 mp_sum,&
                                 mp_sync
  USE string_utilities,    ONLY: compress
  USE termination,         ONLY: stop_memory,&
                                 stop_program
  USE timings,             ONLY: timeset,&
                                 timestop

  use cp_error_handling, only: cp_assert, cp_error_message, cp_error_type
  use cp_log_handling, only: cp_to_string, cp_warning_level, cp_failure_level
  use cp_matrix_utils, only: cp_init, cp_dealloc_ref, cp_next, cp_get, &
       cp_matrix_block_iterator

  IMPLICIT NONE

  PRIVATE

  character(len=*), parameter, private :: moduleN='qs_blacs'

  TYPE blacs_matrix_block_type
    PRIVATE
    INTEGER                           :: ncol_local,nrow_local
    REAL(wp), DIMENSION(:,:), POINTER :: block
  END TYPE blacs_matrix_block_type

  TYPE blacs_matrix_type
    PRIVATE
    CHARACTER(LEN=60)                                      :: name
    INTEGER                                                :: context,&
                                                              ncol_block,&
                                                              ncol_global,&
                                                              nrow_block,&
                                                              nrow_global
    INTEGER, DIMENSION(9)                                  :: descriptor
    ! when using blacs
    ! ncol_block=descriptor(6)
    ! ncol_global=descriptor(4)
    ! nrow_block=descriptor(5)
    ! nrow_global=descriptor(3)
    TYPE(blacs_matrix_block_type), DIMENSION(:,:), POINTER :: p
  END TYPE blacs_matrix_type

  TYPE blacs_matrix_p_type
    TYPE(blacs_matrix_type), POINTER :: blacs_matrix
  END TYPE blacs_matrix_p_type

! *** Public data types ***

  PUBLIC :: blacs_matrix_p_type,&
            blacs_matrix_type

! *** Public subroutines ***

  PUBLIC :: allocate_blacs_matrix,&
            blacs_add,&
            blacs_add_to_element,&
            blacs_gemm,&
            blacs_get_element,&
            blacs_maxval,&
            blacs_set_all,&
            blacs_set_element,&
            blacs_syevx,&
            blacs_symm,&
            blacs_syrk,&
            blacs_trace,&
            copy_blacs_to_blacs_matrix,&
            copy_blacs_to_full_matrix,&
            copy_blacs_to_sparse_matrix,&
            copy_sparse_to_blacs_matrix,&
            deallocate_blacs_matrix,&
            finish_blacs,&
            get_blacs_info,&
            get_blacs_matrix_info,&
            power_blacs_matrix,&
            read_blacs_matrix,&
            replicate_blacs_matrix,&
            start_blacs,&
            symmetrise_blacs_matrix,&
            write_blacs_matrix,&
            blacs_diag_mult,&
            add_blacs_to_block_diag_sm, &
            blacs_cholesky_decompose, &
            blacs_cholesky_reduce, &
            blacs_cholesky_restore

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_blacs_matrix(new_matrix,nrow_global,ncol_global,&
                                   nrow_block,ncol_block,name,context,globenv)

!   Purpose: Allocate a new distributed BLACS matrix.

!   History: - Creation (23.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: new_matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=*), INTENT(IN)              :: name
    INTEGER, INTENT(IN)                       :: context,ncol_block,&
                                                 ncol_global,nrow_block,&
                                                 nrow_global

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_blacs_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    CHARACTER(LEN=40) :: message
    INTEGER           :: group,ierror,ipcol,ipe,iprow,istat,mype,mypcol,&
                         myprow,ncol_local,npcol,npe,nprow,nrow_local,&
                         output_unit,source
    LOGICAL           :: ionode

    INTEGER, DIMENSION(:), ALLOCATABLE :: pcol,prow

#if defined(__parallel)
    INTEGER, EXTERNAL :: blacs_pnum,numroc

#endif
!   ---------------------------------------------------------------------------

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr
    source = globenv%source

    ALLOCATE (new_matrix,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"new_matrix",int_size)

    new_matrix%name = name
#if defined(__parallel)

    CALL blacs_pinfo(mype,npe)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    new_matrix%context = context
    new_matrix%nrow_global = nrow_global
    new_matrix%ncol_global = ncol_global

    new_matrix%nrow_block = MIN(nrow_block,nrow_global/nprow,&
                                ncol_block,ncol_global/npcol)
    new_matrix%ncol_block = new_matrix%nrow_block

    IF ((new_matrix%nrow_block == 0).OR.&
        (new_matrix%ncol_block == 0)) THEN
      CALL stop_program(routine,"More processes than matrix elements",globenv)
    END IF

    nrow_local = numroc(nrow_global,new_matrix%nrow_block,myprow,source,nprow)
    ncol_local = numroc(ncol_global,new_matrix%ncol_block,mypcol,source,npcol)

    ALLOCATE (prow(0:npe-1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"prow",npe*int_size)
    prow(:) = 0
    ALLOCATE (pcol(0:npe-1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pcol",npe*int_size)
    pcol(:) = 0

    prow(mype) = nrow_local
    pcol(mype) = ncol_local

    CALL mp_sum(prow,group)
    CALL mp_sum(pcol,group)

    IF (ionode) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,/,/,(T3,A,I6))")&
        "BLACS INFORMATION (BLACS matrix allocation)",&
        "Matrix name: "//TRIM(name),&
        "Number of rows of the global matrix:    ",new_matrix%nrow_global,&
        "Number of columns of the global matrix: ",new_matrix%ncol_global,&
        "Number of rows of a matrix block:       ",new_matrix%nrow_block,&
        "Number of columns of a matrix block:    ",new_matrix%ncol_block
      WRITE (UNIT=output_unit,FMT="(/,T4,A,/)")&
        "PE      block rows   block columns      rows   columns"
      WRITE (UNIT=output_unit,FMT="(I5,T16,I6,T32,I6,T42,I6,T52,I6)")&
        (ipe,prow(ipe)/new_matrix%nrow_block,pcol(ipe)/new_matrix%nrow_block,&
         prow(ipe),pcol(ipe),ipe=0,npe-1)
    END IF

    ALLOCATE (new_matrix%p(0:nprow-1,0:npcol-1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"new_matrix%p(nprow,npcol)",0)

    DO iprow=0,nprow-1
      DO ipcol=0,npcol-1
        NULLIFY (new_matrix%p(iprow,ipcol)%block)
        ipe = blacs_pnum(context,iprow,ipcol)
        new_matrix%p(iprow,ipcol)%nrow_local = prow(ipe)
        new_matrix%p(iprow,ipcol)%ncol_local = pcol(ipe)
      END DO
    END DO

    ALLOCATE (new_matrix%p(myprow,mypcol)%block(nrow_local,ncol_local),&
              STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"new_matrix%p(myprow,mypcol)%block",&
                       nrow_local*ncol_local*wp_size)
    END IF
    new_matrix%p(myprow,mypcol)%block(:,:) = 0.0_wp

    new_matrix%descriptor(:) = 0
    CALL descinit(new_matrix%descriptor,new_matrix%nrow_global,&
                  new_matrix%ncol_global,new_matrix%nrow_block,&
                  new_matrix%ncol_block,source,source,context,nrow_local,&
                  ierror)

    IF (ierror /= 0) THEN
      WRITE (UNIT=message,FMT="(A,I6)") "Error in descinit: ierror = ",ierror
      CALL compress(message)
      CALL stop_program(routine,message,globenv)
    END IF

    DEALLOCATE (prow,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"prow")
    DEALLOCATE (pcol,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pcol")

#else

    new_matrix%context = 0
    new_matrix%nrow_block = nrow_global
    new_matrix%ncol_block = ncol_global
    new_matrix%nrow_global = nrow_global
    new_matrix%ncol_global = ncol_global
    new_matrix%descriptor(:) = 0
    ALLOCATE (new_matrix%p(source:source,source:source))
    new_matrix%p(source,source)%nrow_local = nrow_global
    new_matrix%p(source,source)%ncol_local = ncol_global
    ALLOCATE (new_matrix%p(source,source)%block(nrow_global,ncol_global),&
              STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"new_matrix%p(source,source)%block",&
                       nrow_global*ncol_global*wp_size)
    END IF
    new_matrix%p(source,source)%block(:,:) = 0.0_wp

#endif
  END SUBROUTINE allocate_blacs_matrix

! *****************************************************************************

  SUBROUTINE blacs_add(alpha,matrix_a,beta,matrix_b,context,globenv)

!   Purpose: Scale and add two BLACS matrices (a <- alpha*a + beta*b).

!   History: - Creation (11.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix_a,matrix_b
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(IN)                      :: alpha,beta
    INTEGER, INTENT(IN)                       :: context

!   *** Local variables ***

    INTEGER :: handle,mypcol,myprow,npcol,nprow,source

    REAL(wp), DIMENSION(:,:), POINTER :: a,b

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_add","I","",handle)

    source = globenv%source
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

#else

    myprow = source
    mypcol = source

#endif
    a => matrix_a%p(myprow,mypcol)%block
    b => matrix_b%p(myprow,mypcol)%block

    IF (alpha == 0.0_wp) THEN
      IF (beta == 0.0_wp) THEN
        a(:,:) = 0.0_wp
      ELSE IF (beta == 1.0_wp) THEN
        a(:,:) = b(:,:)
      ELSE
        a(:,:) = beta*b(:,:)
      END IF
    ELSE IF (beta == 0.0_wp) THEN
      IF (alpha == 1.0_wp) THEN
        RETURN
      ELSE
        a(:,:) = alpha*a(:,:)
      END IF
    ELSE IF (alpha == 1.0_wp) THEN
      IF (beta == 1.0_wp) THEN
        a(:,:) = a(:,:) + b(:,:)
      ELSE
        a(:,:) = a(:,:) + beta*b(:,:)
      END IF
    ELSE IF (beta == 1.0_wp) THEN
      a(:,:) = alpha*a(:,:) + b(:,:)
    ELSE
      a(:,:) = alpha*a(:,:) + beta*b(:,:)
    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_add

! *****************************************************************************

  SUBROUTINE blacs_add_to_element(matrix,irow_global,icol_global,alpha,&
                                  context,globenv)

!   Purpose: Add alpha to the BLACS matrix element
!            matrix(irow_global,icol_global).

!   History: - Creation (08.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(IN)                      :: alpha
    INTEGER, INTENT(IN)                       :: context,icol_global,&
                                                 irow_global

!   *** Local variables ***

    INTEGER :: icol_local,ipcol,iprow,irow_local,mypcol,myprow,npcol,nprow,&
               source

    INTEGER, DIMENSION(9) :: desca

    REAL(wp), DIMENSION(:,:), POINTER :: a

!   ---------------------------------------------------------------------------

    source = globenv%source
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    a => matrix%p(myprow,mypcol)%block
    desca(:) = matrix%descriptor(:)

    CALL infog2l(irow_global,icol_global,desca,nprow,npcol,myprow,mypcol,&
                 irow_local,icol_local,iprow,ipcol)

    IF ((iprow == myprow).AND.(ipcol == mypcol)) THEN
      a(irow_local,icol_local) = a(irow_local,icol_local) + alpha
    END IF

#else

    a => matrix%p(source,source)%block

    a(irow_global,icol_global) = a(irow_global,icol_global) + alpha

#endif
  END SUBROUTINE blacs_add_to_element

! *****************************************************************************

  SUBROUTINE blacs_gemm(transa,transb,m,n,k,alpha,matrix_a,matrix_b,beta,&
                        matrix_c,context,globenv)

!   Purpose: BLACS interface to the BLAS routine dgemm.

!   History: - Creation (07.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix_a,matrix_b,matrix_c
    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=1), INTENT(IN)              :: transa,transb
    REAL(wp), INTENT(IN)                      :: alpha,beta
    INTEGER, INTENT(IN)                       :: context,k,m,n

!   *** Local variables ***

    INTEGER :: handle,lda,ldb,ldc,mypcol,myprow,npcol,nprow,source

    INTEGER, DIMENSION(9) :: desca,descb,descc

    REAL(wp), DIMENSION(:,:), POINTER :: a,b,c

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_gemm","I","",handle)

    source = globenv%source
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    a => matrix_a%p(myprow,mypcol)%block
    desca(:) = matrix_a%descriptor(:)
    b => matrix_b%p(myprow,mypcol)%block
    descb(:) = matrix_b%descriptor(:)
    c => matrix_c%p(myprow,mypcol)%block
    descc(:) = matrix_c%descriptor(:)

    CALL pdgemm(transa,transb,m,n,k,alpha,a,1,1,desca,b,1,1,descb,beta,c,1,1,&
                descc)

#else

    a => matrix_a%p(source,source)%block
    b => matrix_b%p(source,source)%block
    c => matrix_c%p(source,source)%block

    lda = matrix_a%nrow_global
    ldb = matrix_b%nrow_global
    ldc = matrix_c%nrow_global

    CALL dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_gemm

! *****************************************************************************

  SUBROUTINE blacs_get_element(matrix,irow_global,icol_global,alpha,context,&
                               globenv)

!   Purpose: Get the BLACS matrix element (irow_global,icol_global).

!   History: - Creation (22.01.2002, Matthias Krack)
!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(OUT)                     :: alpha
    INTEGER, INTENT(IN)                       :: context,&
                                                 icol_global,&
                                                 irow_global

!   *** Local variables ***

    INTEGER :: icol_local,ipcol,iprow,irow_local,mypcol,myprow,npcol,nprow,&
               source

    INTEGER, DIMENSION(9) :: desca

    REAL(wp), DIMENSION(:,:), POINTER :: a

!   ---------------------------------------------------------------------------

    source = globenv%source
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    a => matrix%p(myprow,mypcol)%block
    desca(:) = matrix%descriptor(:)

    CALL infog2l(irow_global,icol_global,desca,nprow,npcol,myprow,mypcol,&
                 irow_local,icol_local,iprow,ipcol)

    IF ((iprow == myprow).AND.(ipcol == mypcol)) THEN
      alpha = a(irow_local,icol_local)
    ELSE
      alpha = 0.0_wp
    END IF

#else

    alpha = matrix%p(source,source)%block(irow_global,icol_global)

#endif
  END SUBROUTINE blacs_get_element

! *****************************************************************************

  SUBROUTINE blacs_maxval(matrix,a_max,context,globenv)

!   Purpose: Get the maximum absolute element of a BLACS matrix.

!   History: - Creation (11.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(OUT)                     :: a_max
    INTEGER, INTENT(IN)                       :: context

!   *** Local variables ***

    INTEGER :: group,handle,mypcol,myprow,npcol,nprow,source
    REAL(kind=wp), dimension(:,:), pointer    :: my_block

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_maxval","I","",handle)

    source = globenv%source
    group = globenv%group
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

#else

    myprow = source
    mypcol = source

#endif
    my_block => matrix%p(myprow,mypcol)%block
    a_max = MAXVAL(ABS(my_block))
    
    CALL mp_max(a_max,group)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_maxval

! *****************************************************************************

  SUBROUTINE blacs_set_all(matrix,alpha,context,globenv)

!   Purpose: Set the BLACS matrix elements to alpha.

!   History: - Creation (12.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(IN)                      :: alpha
    INTEGER, INTENT(IN)                       :: context

!   *** Local variables ***

    INTEGER :: mypcol,myprow,npcol,nprow,source

!   ---------------------------------------------------------------------------

    source = globenv%source
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    matrix%p(myprow,mypcol)%block(:,:) = alpha

#else

    matrix%p(source,source)%block(:,:) = alpha

#endif
  END SUBROUTINE blacs_set_all

! *****************************************************************************

  SUBROUTINE blacs_set_element(matrix,irow_global,icol_global,alpha,context,&
                               globenv)

!   Purpose: Set the BLACS matrix element (irow_global,icol_global) to alpha.

!   History: - Creation (08.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(IN)                      :: alpha
    INTEGER, INTENT(IN)                       :: context,icol_global,&
                                                 irow_global

!   *** Local variables ***

    INTEGER :: icol_local,ipcol,iprow,irow_local,mypcol,myprow,npcol,nprow,&
               source

    INTEGER, DIMENSION(9) :: desca

    REAL(wp), DIMENSION(:,:), POINTER :: a

!   ---------------------------------------------------------------------------

    source = globenv%source
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    a => matrix%p(myprow,mypcol)%block
    desca(:) = matrix%descriptor(:)

    CALL infog2l(irow_global,icol_global,desca,nprow,npcol,myprow,mypcol,&
                 irow_local,icol_local,iprow,ipcol)

    IF ((iprow == myprow).AND.(ipcol == mypcol)) THEN
      a(irow_local,icol_local) = alpha
    END IF

#else

    matrix%p(source,source)%block(irow_global,icol_global) = alpha

#endif
  END SUBROUTINE blacs_set_element

!******************************************************************
! used to replace a symmetric positive def. matrix by its cholesky 
! decomposition
!******************************************************************

  SUBROUTINE blacs_cholesky_decompose(matrix,context,globenv)
   TYPE(blacs_matrix_type), POINTER           :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: context

    REAL(wp), DIMENSION(:,:), POINTER         :: a
    integer                                   :: info,handle
    INTEGER                                   :: output_unit,source,n
    LOGICAL                                   :: ionode
    integer                                   :: nprow,npcol,myprow,mypcol
    INTEGER, DIMENSION(9)                     :: desca

    CALL timeset("blacs_cholesky_decompose","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr
    source = globenv%source
    n = matrix%nrow_global

#if defined(__parallel)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    a => matrix%p(myprow,mypcol)%block
    desca(:) = matrix%descriptor(:)
    CALL pdpotrf('U',n,a,1,1,desca,info)

#else

    a => matrix%p(source,source)%block

    CALL dpotrf('U',n,a,n,info)

#endif

    if (info.ne.0) then
      CALL stop_program(" blacs_cholesky_decompose","Error dpotrf",globenv)
    endif


    CALL timestop(0.0_wp,handle)

  END  SUBROUTINE blacs_cholesky_decompose

!******************************************************************
! reduce a matrix pencil A,B to normal form
! B has to be cholesky decomposed with  blacs_cholesky_decompose
! before calling this routine
!******************************************************************

  SUBROUTINE blacs_cholesky_reduce(matrix,matrixb,context,globenv)
   TYPE(blacs_matrix_type), POINTER           :: matrix,matrixb
    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: context

    REAL(wp), DIMENSION(:,:), POINTER         :: a,b
    integer                                   :: info,itype,handle
    INTEGER                                   :: output_unit,source,n
    LOGICAL                                   :: ionode
    integer                                   :: nprow,npcol,myprow,mypcol
    INTEGER, DIMENSION(9)                     :: desca,descb
    REAL(wp)                                  :: scale

    CALL timeset("blacs_cholesky_reduce","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr
    source = globenv%source
    n = matrix%nrow_global
    itype =1

#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    a => matrix%p(myprow,mypcol)%block
    desca(:) = matrix%descriptor(:)
    b => matrixb%p(myprow,mypcol)%block
    descb(:) = matrixb%descriptor(:)

    CALL pdsygst(itype,'U',n,a,1,1,desca,b,1,1,descb,scale,info) 
    IF (scale .ne. 1.0_wp) THEN
       ! this is supposed to be one in current version of lapack
       ! if not, eigenvalues have to be scaled by this number
       CALL stop_program(" blacs_cholesky_reduce","scale not equal 1",globenv)
    ENDIF
#else

    a => matrix%p(source,source)%block
    b => matrixb%p(source,source)%block

    CALL dsygst(itype,'U',n,a,n,b,n,info)

#endif

    if (info.ne.0) then
      CALL stop_program(" blacs_cholesky_reduce","Error dsygst",globenv)
    endif


    CALL timestop(0.0_wp,handle)

  END  SUBROUTINE blacs_cholesky_reduce

!******************************************************************
! restore after blacs_cholesky_reduce
! matrix should be the eigenvectors
!******************************************************************
  SUBROUTINE blacs_cholesky_restore(matrix,neig,matrixb,context,globenv)
    TYPE(blacs_matrix_type), POINTER          :: matrix,matrixb
    INTEGER, INTENT(IN)                       :: neig
    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: context

    REAL(wp), DIMENSION(:,:), POINTER         :: a,b
    integer                                   :: info,itype,handle
    INTEGER                                   :: output_unit,source,n
    LOGICAL                                   :: ionode
    REAL(wp)                                  :: alfa
    integer                                   :: nprow,npcol,myprow,mypcol
    INTEGER, DIMENSION(9)                     :: desca,descb

    CALL timeset("blacs_cholesky_restore","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr
    source = globenv%source
    n = matrix%nrow_global
    itype = 1 

#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    a => matrix%p(myprow,mypcol)%block
    desca(:) = matrix%descriptor(:)
    b => matrixb%p(myprow,mypcol)%block
    descb(:) = matrixb%descriptor(:)
    alfa=1.0_wp
    CALL pdtrsm('L','U','N','N',n,neig,alfa,b,1,1,descb,a,1,1,desca)

#else

    ! notice b is the cholesky guy
    a => matrix%p(source,source)%block
    b => matrixb%p(source,source)%block
    alfa=1.0_wp
    CALL dtrsm('L','U','N','N',n,neig,alfa,b,n,a,n)

#endif

    CALL timestop(0.0_wp,handle)

  END  SUBROUTINE blacs_cholesky_restore
   

   

! *****************************************************************************

  SUBROUTINE blacs_syevx(matrix,eigenvectors,eigenvalues,neig,work_syevx,&
                         context,globenv)
!  alternative , needs dsygvx
!  SUBROUTINE blacs_syevx(matrix,eigenvectors,eigenvalues,neig,work_syevx,&
!                         context,globenv,bmatrix)

!   Purpose: Diagonalise the symmetric n by n matrix using the LAPACK library.
!            if bmatrix present, use general solver, will overwrite bmatrix
!            with its cholesky decom.

!   History: - Creation (06.06.2001, Matthias Krack)
!   History: - Mod      (05.2001, Joost VandeVondele)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: eigenvectors,matrix
!    alternative needs dsygvx
!    TYPE(blacs_matrix_type), POINTER, OPTIONAL:: bmatrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(IN)                      :: work_syevx
    INTEGER, INTENT(IN)                       :: context,neig
    REAL(wp), DIMENSION(:), INTENT(OUT)       :: eigenvalues

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE blas_syevx (MODULE qs_blacs)"
    REAL(wp), PARAMETER :: abstol = 0.0_wp,&
                           orfac = -1.0_wp,&
                           vl = 0.0_wp,&
                           vu = 0.0_wp

!   *** Local variables ***

    INTEGER  :: handle,info,istat,liwork,lwork,m,mypcol,myprow,n,nb,nn,np0,&
                npcol,npe,nprow,nq0,nz,output_unit,source,itype
    LOGICAL  :: ionode

    INTEGER, DIMENSION(9) :: desca,descz,descb

    REAL(wp), DIMENSION(:), ALLOCATABLE :: gap,w,work
    INTEGER, DIMENSION(:), ALLOCATABLE  :: iclustr,ifail,iwork
    REAL(wp), DIMENSION(:,:), POINTER   :: a,z,b

#if defined(__parallel)
    INTEGER, EXTERNAL  :: iceil,numroc
#else
    INTEGER, EXTERNAL  :: ilaenv
#endif

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_syevx","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr
    source = globenv%source

    n = matrix%nrow_global

    ALLOCATE (w(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"w",n*wp_size)
    w(:) = 0.0_wp

    eigenvalues(:) = 0.0_wp
#if defined(__parallel)

    IF (matrix%nrow_block /= matrix%ncol_block) THEN
      CALL stop_program(routine,"Invalid blocksize (no square blocks)",globenv)
    END IF

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    a => matrix%p(myprow,mypcol)%block
    desca(:) = matrix%descriptor(:)
    z => eigenvectors%p(myprow,mypcol)%block
    descz(:) = eigenvectors%descriptor(:)

!   *** Get the optimal work storage size ***

    npe = nprow*npcol
    nb = matrix%nrow_block
    nn = MAX(n,nb,2)
    np0 = numroc(nn,nb,0,0,nprow)
    nq0 = MAX(numroc(nn,nb,0,0,npcol),nb)

    lwork = 5*n + MAX(5*nn,np0*nq0) + iceil(neig,npe)*nn + 2*nb*nb +&
            INT(work_syevx*REAL((neig - 1)*n,wp))
    liwork = MAX(3*n + npe + 1,4*n,14) + 2*n

    ALLOCATE (gap(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gap",npe*wp_size)
    gap(:) = 0.0_wp
    ALLOCATE (iclustr(2*npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"iclustr",2*npe*int_size)
    iclustr(:) = 0
    ALLOCATE (ifail(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"ifail",n*int_size)
    ifail(:) = 0
    ALLOCATE (iwork(liwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"iwork",liwork*int_size)
    iwork(:) = 0
    ALLOCATE (work(lwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",lwork*wp_size)
    work(:) = 0.0_wp

!   *** Diagonalise matrix ***
!    if (PRESENT(bmatrix)) then
!       itype=1
!       b => bmatrix%p(myprow,mypcol)%block
!       descb(:) = matrix%descriptor(:)
!       CALL pdsygvx(itype,"V","I","U",n,a,1,1,desca,b,1,1,descb,vl,vu,1,neig, &
!                 abstol,m,nz,w,orfac,&
!                 z,1,1,descz,work,lwork,iwork,liwork,ifail,iclustr,gap,info)
!    else
       CALL pdsyevx("V","I","U",n,a,1,1,desca,vl,vu,1,neig,abstol,m,nz,w,orfac,&
                 z,1,1,descz,work,lwork,iwork,liwork,ifail,iclustr,gap,info)
!    endif


!   *** Error handling ***

    IF (info /= 0) THEN
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
      CALL stop_program(routine,"Error in pdsyevx",globenv)
    END IF

!   *** Release work storage ***

    DEALLOCATE (gap,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gap")
    DEALLOCATE (iclustr,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"iclustr")
    DEALLOCATE (ifail,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"ifail")
    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"iwork")
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work")

#else

    a => matrix%p(source,source)%block
    z => eigenvectors%p(source,source)%block

!   *** Get the optimal work storage size ***

    nb = MAX(ilaenv(1,"DSYTRD","U",n,-1,-1,-1),&
             ilaenv(1,"DORMTR","U",n,-1,-1,-1))

    lwork = MAX((nb + 3)*n,8*n)+n ! sun bug fix
    liwork = 5*n

    ALLOCATE (ifail(n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"ifail",n*int_size)
    ifail(:) = 0
    ALLOCATE (iwork(liwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"iwork",liwork*int_size)
    iwork(:) = 0
    ALLOCATE (work(lwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",lwork*wp_size)
    work(:) = 0.0_wp

!   *** Diagonalise matrix ***
!   dsygvx seems not to be in all lapack versions ...
!    if (PRESENT(bmatrix)) then
!       itype=1
!       b => bmatrix%p(source,source)%block
!       CALL dsygvx(itype,"V","I","U",n,a,n,b,n,vl,vu,1,neig,abstol,m,w,z,n, &
!                   work,lwork,iwork,ifail,info)
!    else
       CALL dsyevx("V","I","U",n,a,n,vl,vu,1,neig,abstol,m,w,z,n,work,lwork,&
                iwork,ifail,info)
!    endif

!   *** Error handling ***

    IF (info /= 0) CALL stop_program(routine,"Error in dsyevx",globenv)

!   *** Release work storage ***

    DEALLOCATE (ifail,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"ifail")
    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"iwork")
    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work")

#endif
    eigenvalues(1:neig) = w(1:neig)
    DEALLOCATE (w,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"w")

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_syevx

! *****************************************************************************

  SUBROUTINE blacs_symm(side,uplo,m,n,alpha,matrix_a,matrix_b,beta,matrix_c,&
                        context,globenv)

!   Purpose: BLACS interface to the BLAS routine dsymm.

!   History: - Creation (07.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix_a,matrix_b,matrix_c
    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=1), INTENT(IN)              :: side,uplo
    REAL(wp), INTENT(IN)                      :: alpha,beta
    INTEGER, INTENT(IN)                       :: context,m,n

!   *** Local variables ***

    INTEGER :: handle,lda,ldb,ldc,mypcol,myprow,npcol,nprow,source

    INTEGER, DIMENSION(9) :: desca,descb,descc

    REAL(wp), DIMENSION(:,:), POINTER :: a,b,c

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_symm","I","",handle)

    source = globenv%source
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    a => matrix_a%p(myprow,mypcol)%block
    desca(:) = matrix_a%descriptor(:)
    b => matrix_b%p(myprow,mypcol)%block
    descb(:) = matrix_b%descriptor(:)
    c => matrix_c%p(myprow,mypcol)%block
    descc(:) = matrix_c%descriptor(:)

    CALL pdsymm(side,uplo,m,n,alpha,a,1,1,desca,b,1,1,descb,beta,c,1,1,descc)

#else

    a => matrix_a%p(source,source)%block
    b => matrix_b%p(source,source)%block
    c => matrix_c%p(source,source)%block

    lda = matrix_a%nrow_global
    ldb = matrix_b%nrow_global
    ldc = matrix_c%nrow_global

    CALL dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_symm

! *****************************************************************************

  SUBROUTINE blacs_syrk(uplo,trans,k,alpha,matrix_a,ia,ja,beta,matrix_c,&
                        context,globenv)

!   Purpose: BLACS interface to the BLAS routine dsyrk.

!   History: - Creation (07.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix_a,matrix_c
    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=1), INTENT(IN)              :: trans,uplo
    REAL(wp), INTENT(IN)                      :: alpha,beta
    INTEGER, INTENT(IN)                       :: context,ia,ja,k

!   *** Local variables ***

    INTEGER :: handle,lda,ldc,mypcol,myprow,n,npcol,nprow,source

    INTEGER, DIMENSION(9) :: desca,descc

    REAL(wp), DIMENSION(:,:), POINTER :: a,c

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_syrk","I","",handle)

    source = globenv%source
    n = matrix_a%nrow_global
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    a => matrix_a%p(myprow,mypcol)%block
    desca(:) = matrix_a%descriptor(:)
    c => matrix_c%p(myprow,mypcol)%block
    descc(:) = matrix_c%descriptor(:)

    CALL pdsyrk(uplo,trans,n,k,alpha,a,ia,ja,desca,beta,c,1,1,descc)

#else

    a => matrix_a%p(source,source)%block
    c => matrix_c%p(source,source)%block

    lda = matrix_a%nrow_global
    ldc = matrix_c%nrow_global

    CALL dsyrk(uplo,trans,n,k,alpha,a(ia,ja),lda,beta,c,ldc)

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_syrk

! *****************************************************************************

  SUBROUTINE blacs_trace(matrix_a,matrix_b,trace,context,globenv)

!   Purpose: Calculate the trace of the product of two BLACS matrices.

!   History: - Creation (11.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix_a,matrix_b
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(OUT)                     :: trace
    INTEGER, INTENT(IN)                       :: context

!   *** Local variables ***

    INTEGER :: group,handle,icol_local,irow_local,mypcol,myprow,ncol_local,&
               npcol,nprow,nrow_local,source

    REAL(wp), DIMENSION(:,:), POINTER :: a,b

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_trace","I","",handle)

    group = globenv%group
    source = globenv%source
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

#else

    myprow = source
    mypcol = source

#endif
    a => matrix_a%p(myprow,mypcol)%block
    b => matrix_b%p(myprow,mypcol)%block

    nrow_local = matrix_a%p(myprow,mypcol)%nrow_local
    ncol_local = matrix_b%p(myprow,mypcol)%ncol_local

    trace = 0.0_wp

    DO icol_local=1,ncol_local
      DO irow_local=1,nrow_local
        trace = trace + a(irow_local,icol_local)*b(irow_local,icol_local)
      END DO
    END DO

    CALL mp_sum(trace,group)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_trace

! *****************************************************************************

  SUBROUTINE copy_blacs_to_blacs_matrix(source,target)

!   Purpose: Copy BLACS matrix to a BLACS matrix of the same type.

!   History: - Creation (08.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER :: source,target

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_blacs_to_blacs_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    INTEGER :: ipcol,iprow,istat,npcol,nprow,ncol_local,nrow_local

!   ---------------------------------------------------------------------------

    nprow = SIZE(source%p,1)
    npcol = SIZE(source%p,2)

    DO iprow=0,nprow-1
      DO ipcol=0,npcol-1
        IF (ASSOCIATED(source%p(iprow,ipcol)%block)) THEN
          IF (.NOT.ASSOCIATED(target%p(iprow,ipcol)%block)) THEN
            nrow_local = target%p(iprow,ipcol)%nrow_local
            ncol_local = target%p(iprow,ipcol)%ncol_local
            ALLOCATE (target%p(iprow,ipcol)%block(nrow_local,ncol_local),&
                      STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"target%p(iprow,ipcol)%block",&
                               nrow_local*ncol_local*wp_size)
            END IF
          END IF
          target%p(iprow,ipcol)%block(:,:) = source%p(iprow,ipcol)%block(:,:)
        ELSE
          IF (ASSOCIATED(target%p(iprow,ipcol)%block)) THEN
            DEALLOCATE (target%p(iprow,ipcol)%block,STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"target%p(iprow,ipcol)%block")
            END IF
          END IF
        END IF
      END DO
    END DO

  END SUBROUTINE copy_blacs_to_blacs_matrix

! *****************************************************************************

  SUBROUTINE copy_blacs_to_full_matrix(blacs_matrix,full_matrix,context,&
                                       globenv)

!   Purpose: Copy a BLACS matrix to a full matrix.

!   History: - Creation (18.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: blacs_matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: context
    REAL(wp), DIMENSION(:,:), POINTER         :: full_matrix

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_blacs_to_full_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    INTEGER :: handle,icol_global,icol_local,ipcol,ipe,iprow,irow_global,&
               irow_local,istat,mypcol,mype,myprow,ncol_block,ncol_global,&
               ncol_local,npcol,npe,nprow,nrow_block,nrow_global,nrow_local,&
               source
    LOGICAL :: ionode

    REAL(wp), DIMENSION(:,:), POINTER :: blacs_block

#if defined(__parallel)
    INTEGER, EXTERNAL :: blacs_pnum,indxl2g

#endif
!   ---------------------------------------------------------------------------

    CALL timeset("copy_blacs_to_full_matrix","I","",handle)

    ionode = globenv%ionode
    source = globenv%source

    nrow_global = blacs_matrix%nrow_global
    ncol_global = blacs_matrix%ncol_global

    IF (ionode) THEN
      ALLOCATE (full_matrix(nrow_global,ncol_global),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"full_matrix",&
                                       nrow_global*ncol_global*wp_size)
      full_matrix(:,:) = 0.0_wp
    END IF
#if defined(__parallel)

    CALL blacs_pinfo(mype,npe)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    nrow_block = blacs_matrix%nrow_block
    ncol_block = blacs_matrix%ncol_block

    DO iprow=0,nprow-1
      DO ipcol=0,npcol-1

        ipe = blacs_pnum(context,iprow,ipcol)

        nrow_local = blacs_matrix%p(iprow,ipcol)%nrow_local
        ncol_local = blacs_matrix%p(iprow,ipcol)%ncol_local

        IF (ionode) THEN

          IF (ipe /= mype) THEN
            ALLOCATE (blacs_matrix%p(iprow,ipcol)%block(nrow_local,&
                                                        ncol_local),&
                      STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"blacs_matrix%p(iprow,ipcol)%block",&
                               nrow_local*ncol_local*wp_size)
            END IF
            CALL dgerv2d(context,nrow_local,ncol_local,&
                         blacs_matrix%p(iprow,ipcol)%block,nrow_local,&
                         iprow,ipcol)
          END IF

          blacs_block => blacs_matrix%p(iprow,ipcol)%block

          DO icol_local=1,ncol_local
            icol_global = indxl2g(icol_local,ncol_block,ipcol,source,npcol)
            DO irow_local=1,nrow_local
              irow_global = indxl2g(irow_local,nrow_block,iprow,source,nprow)
              full_matrix(irow_global,icol_global) = blacs_block(irow_local,&
                                                                 icol_local)
            END DO
          END DO

          IF (ipe /= mype) THEN
            DEALLOCATE (blacs_matrix%p(iprow,ipcol)%block,STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"blacs_matrix%p(iprow,ipcol)%block")
            END IF
          END IF

        ELSE

          IF (ipe == mype) THEN
            CALL dgesd2d(context,nrow_local,ncol_local,&
                         blacs_matrix%p(iprow,ipcol)%block,nrow_local,&
                         source,source)
          END IF

        END IF

        CALL blacs_barrier(context,"A")

      END DO
    END DO

#else

    full_matrix(:,:) = blacs_matrix%p(source,source)%block(:,:)

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE copy_blacs_to_full_matrix

! *****************************************************************************

  SUBROUTINE copy_blacs_to_sparse_matrix(blacs_matrix,sparse_matrix,context,&
                                         globenv)

!   Purpose: Copy a BLACS matrix to a sparse matrix. The BLACS matrix blocks
!            are deallocated during the copy procedure.

!   History: - Creation (06.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: blacs_matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(real_matrix_type), POINTER           :: sparse_matrix
    INTEGER, INTENT(IN)                       :: context

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_blacs_to_sparse_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: group,handle,iblock_col,iblock_row,icol,icol_global,icol_local,&
               ipcol,ipe,iprow,irow,irow_global,irow_local,istat,jpcol,jprow,&
               mypcol,mype,myprow,nblock_row,ncol_block,ncol_local,npcol,npe,&
               nprow,nrow_block,nrow_local,source

    INTEGER, DIMENSION(:), POINTER    :: first_col,first_row,last_col,last_row
    REAL(wp), DIMENSION(:,:), POINTER :: blacs_block,sparse_block

#if defined(__parallel)
    INTEGER, EXTERNAL :: blacs_pnum,indxg2l,indxg2p

#endif
!   ---------------------------------------------------------------------------

    CALL timeset("copy_blacs_to_sparse_matrix","I","",handle)

    group = globenv%group
    source = globenv%source

    CALL get_matrix_info(matrix=sparse_matrix,&
                         nblock_row=nblock_row,&
                         first_row=first_row,&
                         first_col=first_col,&
                         last_row=last_row,&
                         last_col=last_col)
#if defined(__parallel)

    CALL blacs_pinfo(mype,npe)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    nrow_block = blacs_matrix%nrow_block
    ncol_block = blacs_matrix%ncol_block

    DO iprow=0,nprow-1
      DO ipcol=0,npcol-1

        ipe = blacs_pnum(context,iprow,ipcol)

        IF (ipe /= mype) THEN
          nrow_local = blacs_matrix%p(iprow,ipcol)%nrow_local
          ncol_local = blacs_matrix%p(iprow,ipcol)%ncol_local
          ALLOCATE (blacs_matrix%p(iprow,ipcol)%block(nrow_local,ncol_local),&
                    STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"blacs_matrix%p(iprow,ipcol)%block",&
                             nrow_local*ncol_local*wp_size)
          END IF
        END IF

        blacs_block => blacs_matrix%p(iprow,ipcol)%block

        CALL mp_bcast(blacs_block,ipe,group)

        DO iblock_row=1,nblock_row

          block_node => first_block_node(matrix=sparse_matrix,&
                                         block_row=iblock_row)

          DO WHILE (ASSOCIATED(block_node))

            CALL get_block_node(block_node=block_node,&
                                block_col=iblock_col,&
                                block=sparse_block)

            icol = 1

            DO icol_global=first_col(iblock_col),last_col(iblock_col)

              jpcol = indxg2p(icol_global,ncol_block,mypcol,source,npcol)

              IF (jpcol == ipcol) THEN

                icol_local = indxg2l(icol_global,ncol_block,mypcol,source,&
                                     npcol)

                irow = 1

                DO irow_global=first_row(iblock_row),last_row(iblock_row)

                  jprow = indxg2p(irow_global,nrow_block,myprow,source,nprow)

                  IF (jprow == iprow) THEN

                    irow_local = indxg2l(irow_global,nrow_block,myprow,source,&
                                         nprow)

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
          DEALLOCATE (blacs_matrix%p(iprow,ipcol)%block,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"blacs_matrix%p(iprow,ipcol)%block")
          END IF
        END IF

      END DO
    END DO

#else

    blacs_block => blacs_matrix%p(source,source)%block

    DO iblock_row=1,nblock_row

      block_node => first_block_node(matrix=sparse_matrix,&
                                     block_row=iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            block=sparse_block)

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

  SUBROUTINE copy_sparse_to_blacs_matrix(sparse_matrix,blacs_matrix,context,&
                                         globenv)

!   Purpose: Copy a sparse matrix to a BLACS matrix. The BLACS matrix blocks
!            are allocated during the copy procedure.

!   History: - Creation (05.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: blacs_matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(real_matrix_type), POINTER           :: sparse_matrix
    INTEGER, INTENT(IN)                       :: context

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_sparse_to_blacs_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: group,handle,iblock_col,iblock_row,icol,icol_global,icol_local,&
               ipcol,ipe,iprow,irow,irow_global,irow_local,istat,jpcol,jprow,&
               mypcol,mype,myprow,nblock_row,ncol_block,ncol_local,npcol,npe,&
               nprow,nrow_block,nrow_local,source

    INTEGER, DIMENSION(:), POINTER    :: first_col,first_row,last_col,last_row
    REAL(wp), DIMENSION(:,:), POINTER :: blacs_block,sparse_block

#if defined(__parallel)
    INTEGER, EXTERNAL :: blacs_pnum,indxg2l,indxg2p

#endif
!   ---------------------------------------------------------------------------

    CALL timeset("copy_sparse_to_blacs_matrix","I","",handle)

    group = globenv%group
    source = globenv%source

    CALL get_matrix_info(matrix=sparse_matrix,&
                         nblock_row=nblock_row,&
                         first_row=first_row,&
                         first_col=first_col,&
                         last_row=last_row,&
                         last_col=last_col)
#if defined(__parallel)

    CALL blacs_pinfo(mype,npe)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    nrow_block = blacs_matrix%nrow_block
    ncol_block = blacs_matrix%ncol_block

    DO iprow=0,nprow-1
      DO ipcol=0,npcol-1

        ipe = blacs_pnum(context,iprow,ipcol)

        IF (ipe /= mype) THEN
          nrow_local = blacs_matrix%p(iprow,ipcol)%nrow_local
          ncol_local = blacs_matrix%p(iprow,ipcol)%ncol_local
          ALLOCATE (blacs_matrix%p(iprow,ipcol)%block(nrow_local,ncol_local),&
                    STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"blacs_matrix%p(iprow,ipcol)%block",&
                             nrow_local*ncol_local*wp_size)
          END IF
        END IF

        blacs_block => blacs_matrix%p(iprow,ipcol)%block

        blacs_block(:,:) = 0.0_wp

        DO iblock_row=1,nblock_row

          block_node => first_block_node(matrix=sparse_matrix,&
                                         block_row=iblock_row)

          DO WHILE (ASSOCIATED(block_node))

            CALL get_block_node(block_node=block_node,&
                                block_col=iblock_col,&
                                block=sparse_block)

            icol = 1

            DO icol_global=first_col(iblock_col),last_col(iblock_col)

              jpcol = indxg2p(icol_global,ncol_block,mypcol,source,npcol)

              IF (jpcol == ipcol) THEN

                icol_local = indxg2l(icol_global,ncol_block,mypcol,source,&
                                     npcol)

                irow = 1

                DO irow_global=first_row(iblock_row),last_row(iblock_row)

                  jprow = indxg2p(irow_global,nrow_block,myprow,source,nprow)

                  IF (jprow == iprow) THEN

                    irow_local = indxg2l(irow_global,nrow_block,myprow,source,&
                                         nprow)

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
          DEALLOCATE (blacs_matrix%p(iprow,ipcol)%block,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"blacs_matrix%p(iprow,ipcol)%block")
          END IF
        END IF

      END DO
    END DO

#else

    IF (.NOT.ASSOCIATED(blacs_matrix%p(source,source)%block)) THEN
      nrow_local = blacs_matrix%p(source,source)%nrow_local
      ncol_local = blacs_matrix%p(source,source)%ncol_local
      ALLOCATE (blacs_matrix%p(source,source)%block(nrow_local,ncol_local),&
                STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"blacs_matrix%p(source,source)%block",&
                         nrow_local*ncol_local*wp_size)
      END IF
    END IF

    blacs_block => blacs_matrix%p(source,source)%block

    blacs_block(:,:) = 0.0_wp

    DO iblock_row=1,nblock_row

      block_node => first_block_node(matrix=sparse_matrix,&
                                     block_row=iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            block=sparse_block)

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

  SUBROUTINE deallocate_blacs_matrix(matrix)

!   Purpose: Deallocate a distributed BLACS matrix.

!   History: - Creation (08.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER :: matrix

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_blacs_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    INTEGER :: ipcol,iprow,istat,npcol,nprow

!   ---------------------------------------------------------------------------

    matrix%name = ""

    matrix%context = 0
    matrix%nrow_block = 0
    matrix%ncol_block = 0
    matrix%nrow_global = 0
    matrix%ncol_global = 0
    matrix%descriptor(:) = 0

    IF (ASSOCIATED(matrix%p)) THEN

      nprow = SIZE(matrix%p,1)
      npcol = SIZE(matrix%p,2)

      DO iprow=0,nprow-1
        DO ipcol=0,npcol-1
          IF (ASSOCIATED(matrix%p(iprow,ipcol)%block)) THEN
            DEALLOCATE (matrix%p(iprow,ipcol)%block,STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"matrix%p(iprow,ipcol)%block")
            END IF
          END IF
        END DO
      END DO

      DEALLOCATE (matrix%p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"matrix%p")

    END IF

    DEALLOCATE (matrix,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"matrix")

  END SUBROUTINE deallocate_blacs_matrix

! *****************************************************************************

  SUBROUTINE finish_blacs(context,globenv)

!   Purpose: Release the resources of a BLACS context.

!   History: - Creation (22.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: context

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE finish_blacs (MODULE qs_blacs)"

!   *** Local variables ***

    INTEGER :: group,ipe,istat,mype,npe,output_unit
    LOGICAL :: ionode

    INTEGER, DIMENSION(:), ALLOCATABLE :: pcon

!   ---------------------------------------------------------------------------

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr
#if defined(__parallel)

    IF (globenv%print%blacs_info) THEN

      CALL blacs_pinfo(mype,npe)

      ALLOCATE (pcon(0:npe-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcon",npe*int_size)
      pcon(:) = 0

      pcon(mype) = context

      CALL mp_sum(pcon,group)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
          "BLACS INFORMATION (BLACS finished)"
        WRITE (UNIT=output_unit,FMT="(/,T3,A,/)")&
          " PE   BLACS context"
        WRITE (UNIT=output_unit,FMT="(I5,T10,I12)")&
          (ipe,pcon(ipe),ipe=0,npe-1)
      END IF

      DEALLOCATE (pcon,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcon")

    END IF

    CALL mp_sync(group)

    CALL blacs_gridexit(context)

#endif
  END SUBROUTINE finish_blacs

! *****************************************************************************

  SUBROUTINE get_blacs_info(context,globenv,my_process_row,my_process_column,&
                            my_process_number,number_of_process_rows,&
                            number_of_process_columns,number_of_processes)

!   Purpose: Return informations about the specified BLACS context.

!   History: - Creation (19.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: context
    INTEGER, OPTIONAL, INTENT(OUT)            :: my_process_column,&
                                                 my_process_number,&
                                                 my_process_row,&
                                                 number_of_process_columns,&
                                                 number_of_process_rows,&
                                                 number_of_processes

!   *** Local variables ***

    INTEGER :: mypcol,mype,myprow,npcol,npe,nprow,source

#if defined(__parallel)
    INTEGER, EXTERNAL :: blacs_pnum

#endif
!   ---------------------------------------------------------------------------

    source = globenv%source
#if defined(__parallel)

    CALL blacs_pinfo(mype,npe)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

#else

    myprow = source
    mypcol = source
    mype = source
    nprow = 1
    npcol = 1
    npe = 1

#endif
    IF (PRESENT(my_process_row)) my_process_row = myprow
    IF (PRESENT(my_process_column)) my_process_column = mypcol
    IF (PRESENT(my_process_number)) my_process_number = mype
    IF (PRESENT(number_of_process_rows)) number_of_process_rows = nprow
    IF (PRESENT(number_of_process_columns)) number_of_process_columns = npcol
    IF (PRESENT(number_of_processes)) number_of_processes = npe

  END SUBROUTINE get_blacs_info

! *****************************************************************************

  SUBROUTINE get_blacs_matrix_info(matrix,name,nrow_global,ncol_global,&
                                   nrow_block,ncol_block)

!   Purpose: Return informations about the specified BLACS matrix.

!   History: - Creation (08.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER        :: matrix
    CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: name
    INTEGER, OPTIONAL, INTENT(OUT)          :: ncol_block,ncol_global,&
                                               nrow_block,nrow_global

!   ---------------------------------------------------------------------------

    IF (PRESENT(name)) name = matrix%name
    IF (PRESENT(nrow_global)) nrow_global = matrix%nrow_global
    IF (PRESENT(ncol_global)) ncol_global = matrix%ncol_global
    IF (PRESENT(nrow_block)) nrow_block = matrix%nrow_block
    IF (PRESENT(ncol_block)) ncol_block = matrix%ncol_block

  END SUBROUTINE get_blacs_matrix_info

! *****************************************************************************

  SUBROUTINE power_blacs_matrix(matrix,work,exponent,threshold,n_dependent,&
                                work_syevx,context,globenv)

!   Purpose: Raise the real symmetric n by n matrix to the power given by
!            exponent. All eigenvectors with a corresponding eigenvalue lower
!            than threshold are quenched.

!   History: - Creation (29.03.1999, Matthias Krack)
!            - Parallelised using BLACS and ScaLAPACK (06.06.2001, MK)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix,work
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(IN)                      :: exponent,threshold,work_syevx
    INTEGER, INTENT(IN)                       :: context
    INTEGER, INTENT(OUT)                      :: n_dependent

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE power_blacs_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    REAL(wp) :: f,p
    INTEGER  :: handle,icol_global,icol_local,ipcol,iprow,irow_global,&
                irow_local,istat,mypcol,myprow,ncol_block,ncol_global,npcol,&
                nprow,nrow_block,nrow_global,source

    REAL(wp), DIMENSION(:), ALLOCATABLE :: eigenvalues
    REAL(wp), DIMENSION(:,:), POINTER   :: eigenvectors

#if defined(__parallel)
    INTEGER, EXTERNAL :: indxg2l,indxg2p

#endif
!   ---------------------------------------------------------------------------

    CALL timeset("power_blacs_matrix","I","",handle)

    source = globenv%source
    n_dependent = 0
    p = 0.5_wp*exponent

    nrow_global = matrix%nrow_global
    ncol_global = matrix%ncol_global

    ALLOCATE (eigenvalues(ncol_global),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"eigenvalues",ncol_global*wp_size)
    eigenvalues(:) = 0.0_wp

!   *** Compute the eigenvectors and eigenvalues ***

    CALL blacs_syevx(matrix,work,eigenvalues,ncol_global,work_syevx,context,&
                     globenv)

#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    nrow_block = work%nrow_block
    ncol_block = work%ncol_block

    eigenvectors => work%p(myprow,mypcol)%block

!   *** Build matrix**exponent with eigenvector quenching ***

    p = 0.5_wp*exponent

    n_dependent = 0

    DO icol_global=1,ncol_global

      IF (eigenvalues(icol_global) < threshold) THEN

        n_dependent = n_dependent + 1

        ipcol = indxg2p(icol_global,ncol_block,mypcol,source,npcol)

        IF (mypcol == ipcol) THEN
          icol_local = indxg2l(icol_global,ncol_block,mypcol,source,npcol)
          DO irow_global=1,nrow_global
            iprow = indxg2p(irow_global,nrow_block,myprow,source,nprow)
            IF (myprow == iprow) THEN
              irow_local = indxg2l(irow_global,nrow_block,myprow,source,nprow)
              eigenvectors(irow_local,icol_local) = 0.0_wp
            END IF
          END DO
        END IF

      ELSE

        f = eigenvalues(icol_global)**p

        ipcol = indxg2p(icol_global,ncol_block,mypcol,source,npcol)

        IF (mypcol == ipcol) THEN
          icol_local = indxg2l(icol_global,ncol_block,mypcol,source,npcol)
          DO irow_global=1,nrow_global
            iprow = indxg2p(irow_global,nrow_block,myprow,source,nprow)
            IF (myprow == iprow) THEN
              irow_local = indxg2l(irow_global,nrow_block,myprow,source,nprow)
              eigenvectors(irow_local,icol_local) =&
                f*eigenvectors(irow_local,icol_local)
            END IF
          END DO
        END IF

      END IF

    END DO

#else

    eigenvectors => work%p(source,source)%block

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
    CALL blacs_syrk("U","N",ncol_global,1.0_wp,work,1,1,0.0_wp,matrix,context,&
                    globenv)

    DEALLOCATE (eigenvalues,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"eigenvalues")

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE power_blacs_matrix

! *****************************************************************************

  SUBROUTINE read_blacs_matrix(matrix,lunit,context,globenv)

!   Purpose: Read a BLACS matrix from the logical unit number "lunit".

!   History: - Creation (19.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: context,lunit

!   *** Local variables ***

    INTEGER :: i,j,mypcol,myprow,ncol_local,npcol,nprow,nrow_local,source

!   ---------------------------------------------------------------------------

    source =  globenv%source
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

#else

    myprow = source
    mypcol = source

#endif
    nrow_local = matrix%p(myprow,mypcol)%nrow_local
    ncol_local = matrix%p(myprow,mypcol)%ncol_local

    READ (UNIT=lunit) ((matrix%p(myprow,mypcol)%block(i,j),i=1,nrow_local),&
                                                           j=1,ncol_local)

  END SUBROUTINE read_blacs_matrix

! *****************************************************************************

  SUBROUTINE replicate_blacs_matrix(prototype_matrix,new_matrix,name)

!   Purpose: Allocate a distributed BLACS matrix using a prototype matrix.

!   History: - Creation (08.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER :: new_matrix,prototype_matrix
    CHARACTER(LEN=*), INTENT(IN)     :: name

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE replicate_blacs_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    INTEGER :: ipcol,iprow,istat,ncol_local,npcol,nprow,nrow_local

!   ---------------------------------------------------------------------------

    ALLOCATE (new_matrix,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"new_matrix",int_size)

    new_matrix%name = name

    new_matrix%context = prototype_matrix%context
    new_matrix%nrow_block = prototype_matrix%nrow_block
    new_matrix%ncol_block = prototype_matrix%ncol_block
    new_matrix%nrow_global = prototype_matrix%nrow_global
    new_matrix%ncol_global = prototype_matrix%ncol_global
    new_matrix%descriptor(:) = prototype_matrix%descriptor(:)

    nprow = SIZE(prototype_matrix%p,1)
    npcol = SIZE(prototype_matrix%p,2)

    ALLOCATE (new_matrix%p(0:nprow-1,0:npcol-1),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"new_matrix%p",nprow*npcol*int_size)
    END IF

    DO iprow=0,nprow-1
      DO ipcol=0,npcol-1
        NULLIFY (new_matrix%p(iprow,ipcol)%block)
        new_matrix%p(iprow,ipcol)%nrow_local =&
          prototype_matrix%p(iprow,ipcol)%nrow_local
        new_matrix%p(iprow,ipcol)%ncol_local =&
          prototype_matrix%p(iprow,ipcol)%ncol_local
        IF (ASSOCIATED(prototype_matrix%p(iprow,ipcol)%block)) THEN
          nrow_local = new_matrix%p(iprow,ipcol)%nrow_local
          ncol_local = new_matrix%p(iprow,ipcol)%ncol_local
          ALLOCATE (new_matrix%p(iprow,ipcol)%block(nrow_local,ncol_local),&
                    STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"new_matrix%p(iprow,ipcol)%block",&
                             nrow_local*ncol_local*wp_size)
          END IF
          new_matrix%p(iprow,ipcol)%block(:,:) =&
            prototype_matrix%p(iprow,ipcol)%block(:,:)
        END IF
      END DO
    END DO

  END SUBROUTINE replicate_blacs_matrix

! *****************************************************************************

  SUBROUTINE start_blacs(nprow,npcol,context,globenv)

!   Purpose: Initialize a BLACS process grid. The BLACS context is returned.

!   History: - Creation (22.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(INOUT)                    :: context,npcol,nprow

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE start_blacs (MODULE qs_blacs)"

!   *** Local variables ***

    INTEGER :: group,ipe,istat,mypcol,mype,myprow,npe,output_unit,source
    LOGICAL :: ionode

    INTEGER, DIMENSION(:), ALLOCATABLE :: pcol,pcon,prow

!   ---------------------------------------------------------------------------

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr
    source = globenv%source
    context = 0
#if defined(__parallel)

    CALL blacs_pinfo(mype,npe)
    CALL blacs_get(-1,0,context)

    IF (nprow*npcol /= npe) THEN
      DO ipe=CEILING(SQRT(REAL(npe,wp))),npe
        IF (MODULO(npe,ipe) == 0) THEN
          nprow = ipe
          npcol = npe/nprow
          EXIT
        END IF
      END DO
    END IF

    CALL blacs_gridinit(context,"Row-major",nprow,npcol)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    IF (globenv%print%blacs_info) THEN

      ALLOCATE (prow(0:npe-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"prow",npe*int_size)
      prow(:) = 0
      ALLOCATE (pcol(0:npe-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcol",npe*int_size)
      pcol(:) = 0
      ALLOCATE (pcon(0:npe-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcon",npe*int_size)
      pcon(:) = 0

      prow(mype) = myprow
      pcol(mype) = mypcol
      pcon(mype) = context

      CALL mp_sum(prow,group)
      CALL mp_sum(pcol,group)
      CALL mp_sum(pcon,group)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,(T3,A,T32,I6))")&
          "BLACS INFORMATION (BLACS started)",&
          "Number of processes:         ",nprow*npcol,&
          "Number of process rows:      ",nprow,&
          "Number of process columns:   ",npcol
        WRITE (UNIT=output_unit,FMT="(/,T3,A,/)")&
          " PE   process row   process column   BLACS context"
        WRITE (UNIT=output_unit,FMT="(I5,T14,I6,T31,I6,T41,I12)")&
          (ipe,prow(ipe),pcol(ipe),pcon(ipe),ipe=0,npe-1)
      END IF

      DEALLOCATE (prow,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"prow")
      DEALLOCATE (pcol,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcol")
      DEALLOCATE (pcon,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcon")

    END IF

#endif
  END SUBROUTINE start_blacs

! *****************************************************************************

  SUBROUTINE symmetrise_blacs_matrix(matrix,work,context,globenv)

!   Purpose: Symmetrise a symmetric BLACS matrix.

!   History: - Creation (12.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix,work
    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: context

!   *** Local variables ***

    INTEGER :: handle,icol_global,icol_local,ipcol,iprow,irow_global,&
               irow_local,mypcol,myprow,ncol_block,ncol_global,ncol_local,&
               npcol,nprow,nrow_block,nrow_global,nrow_local,source

    INTEGER, DIMENSION(9) :: desca,descc

    REAL(wp), DIMENSION(:,:), POINTER :: a,c

#if defined(__parallel)
    INTEGER, EXTERNAL :: indxl2g

#endif
!   ---------------------------------------------------------------------------

    CALL timeset("symmetrise_blacs_matrix","I","",handle)

    source = globenv%source
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    nrow_global = matrix%nrow_global
    ncol_global = matrix%ncol_global

    nrow_block = matrix%nrow_block
    ncol_block = matrix%ncol_block

    nrow_local = matrix%p(myprow,mypcol)%nrow_local
    ncol_local = matrix%p(myprow,mypcol)%ncol_local

    a => work%p(myprow,mypcol)%block
    desca(:) = work%descriptor(:)
    c => matrix%p(myprow,mypcol)%block
    descc(:) = matrix%descriptor(:)

    DO icol_local=1,ncol_local
      icol_global = indxl2g(icol_local,ncol_block,mypcol,source,npcol)
      DO irow_local=1,nrow_local
        irow_global = indxl2g(irow_local,nrow_block,myprow,source,nprow)
        IF (irow_global > icol_global) THEN
          c(irow_local,icol_local) = 0.0_wp
        ELSE IF (irow_global == icol_global) THEN
          c(irow_local,icol_local) = 0.5_wp*c(irow_local,icol_local)
        END IF
      END DO
    END DO

    a(:,:) = c(:,:)

    CALL pdtran(nrow_global,ncol_global,1.0_wp,a,1,1,desca,1.0_wp,c,1,1,descc)

#else

    a => matrix%p(source,source)%block

    CALL symmetrize_matrix(a,"upper_to_lower")

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE symmetrise_blacs_matrix

! *****************************************************************************

  SUBROUTINE write_blacs_matrix(matrix,lunit,context,globenv)

!   Purpose: Write a BLACS matrix to the logical unit number "lunit".

!   History: - Creation (19.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: context,lunit

!   *** Local variables ***

    INTEGER :: i,j,mypcol,myprow,ncol_local,npcol,nprow,nrow_local,source

!   ---------------------------------------------------------------------------

    source =  globenv%source
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

#else

    myprow = source
    mypcol = source

#endif
    nrow_local = matrix%p(myprow,mypcol)%nrow_local
    ncol_local = matrix%p(myprow,mypcol)%ncol_local

    WRITE (UNIT=lunit) ((matrix%p(myprow,mypcol)%block(i,j),i=1,nrow_local),&
                                                            j=1,ncol_local)

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
  subroutine add_blacs_to_block_diag_sm(source_m, target_m, &
       global_env, d_struct, error)
    TYPE(blacs_matrix_type), POINTER :: source_m
    type(real_matrix_type), intent(inout), target :: target_m
    type(global_environment_type), intent(in) :: global_env
    integer, dimension(:), pointer, optional :: d_struct
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: handle, i, nrows, ncols, stat, j
    integer, dimension(:), pointer :: my_d_struct, first_col,first_row,&
         last_col, last_row
    integer :: start_proc_row,start_proc_col,n_blacs_proc, sm_block_proc,&
         nprow,npcol,myprow,mypcol,b_block_row_start,b_block_row_stop,&
         b_block_col_start,b_block_col_stop,sm_block_nr,b_block_row,&
         b_block_col, start_row_of_b, end_row_of_b, start_row_of_sm,&
         end_row_of_sm, start_col_of_b, end_col_of_b, start_col_of_sm,&
         end_col_of_sm
    character(len=*), parameter :: routineN='add_blacs_to_block_diag_sm',&
         routineP=moduleN//':'//routineN
    type(real_matrix_type), pointer :: target_ptr
    type(cp_matrix_block_iterator) :: iterator
    real(kind=wp), dimension(:,:), pointer :: block_val
    failure=.false.
    target_ptr => target_m
    nullify(my_d_struct)

    call timeset(routineN//','//moduleN,'I',"",handle)
    CPPrecondition(associated(target_ptr),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       call get_matrix_info(target_ptr, nrow=nrows, ncol=ncols)
       CPAssert(ncols==nrows,cp_warning_level,routineP,error,failure)
       ncols=min(nrows,ncols)
       if (present(d_struct)) my_d_struct=d_struct
       if (.not.associated(my_d_struct)) then
! build distribution structure
! assume every diagonal block is on one (and just one) proc
          allocate(my_d_struct(ncols),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          my_d_struct=0
          call cp_init(iterator, matrix=target_m, error=error)
          do
             if (.not.cp_next(iterator)) exit
             call cp_get(iterator,block_row=i,block_col=j,error=error)
             if (i==j) then
                my_d_struct(i)=global_env%mepos
             else
                call cp_error_message(cp_warning_level,&
                     routineP,"there is an off diagonal block! in "//&
                     CPSourceFileRef,error)
! just as info (only on the local proc)
             end if
          end do
          call cp_dealloc_ref(iterator,error=error)
          call mp_sum(my_d_struct,global_env%group)
       end if
    end if
    CPAssert(associated(my_d_struct),cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       CPAssert(all(my_d_struct<global_env%num_pe),cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
!group = global_env%group
!source = global_env%source

       call get_matrix_info(matrix=target_ptr,&
            first_row=first_row,&
            first_col=first_col,&
            last_row=last_row,&
            last_col=last_col)
#if (__parallel)
!MK       call blacs_pinfo(blacs_proc_nr, n_blacs_proc)
!MK       CPPrecondition(blacs_proc_nr==global_env%mepos,cp_failure_level,routineP,error,failure)
!MK       call blacs_gridinfo(source_m%context,nprow,npcol,myprow,mypcol)
!MK       start_proc_row=source_m%descriptor(7)
!MK       start_proc_col=source_m%descriptor(8)
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
       do sm_block_nr=1,ncols
          call get_block_node(target_ptr,sm_block_nr,sm_block_nr,&
               block=block_val)
          block_val=block_val+source_m%p(global_env%source,global_env%source)%&
               block(first_row(sm_block_nr):last_row(sm_block_nr),&
               first_col(sm_block_nr):last_col(sm_block_nr))
       end do
#endif
    end if
    call timestop(0.0_wp,handle)
  end subroutine add_blacs_to_block_diag_sm
!***************************************************************************

!!****f* blacs/blacs_diag_mult [1.0] *
!!
!!   NAME
!!     blacs_diag_mult
!!
!!   FUNCTION
!!     multiplies a blacs matrix with a diagonal matrix given with a vector
!!     If the vector is smaller than the size of the matrix, then the
!!     other elements are considered 0.
!!     c = alpha diag(diagonal) b + beta c
!!
!!   NOTES
!!     if aliasing is accepted target_m and matrix might be the same object.
!!     Actually the vector could also be distributed, but for now it is global
!!
!!   INPUTS
!!     diagonal: the value of the diagonal elements
!!     matrix: the second matrix that is multiplied
!!     alpha: defaults to 1
!!     beta: defaults to 0
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
  subroutine blacs_diag_mult(diagonal,b,c,alpha,beta,global_env,error)
    real(kind=wp), dimension(:), intent(in) :: diagonal
    TYPE(blacs_matrix_type), POINTER :: b,c
    real(kind=wp), intent(in), optional :: alpha, beta
    type(global_environment_type), intent(in), target :: global_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: handle, global_row, my_p_row, my_p_col, size_diag,&
         n_p_rows, n_p_cols, i, j
    real(kind=wp) :: al, be
    character(len=*), parameter :: routineN='blacs_diag_mult',&
         routineP=moduleN//':'//routineN
    real(wp), dimension(:,:), pointer :: block_b, block_c
#if defined(__parallel)
    INTEGER, EXTERNAL :: indxg2l
    call blacs_gridinfo(b%context,n_p_rows,n_p_cols,my_p_row,my_p_col)
#else
    my_p_row=global_env% source; my_p_col=global_env% source
    n_p_rows=1; n_p_cols=1
#endif
    failure=.false.
    al=1.0_wp
    be=0.0_wp

    call timeset(routineN//','//moduleN,'I',"",handle)
    CPPrecondition(b% context==c% context,cp_failure_level,routineP,error,failure)
    CPPrecondition(b% ncol_global==c% ncol_global,cp_failure_level,routineP,error,failure)
    CPPrecondition(b% nrow_global==c% nrow_global,cp_failure_level,routineP,error,failure)
    CPPrecondition(b% ncol_block==c% ncol_block,cp_failure_level,routineP,error,failure)
    CPPrecondition(b% nrow_block==c% nrow_block,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (present(alpha)) al=alpha
       if (present(beta)) be=beta
       size_diag=size(diagonal)
       block_b => b%p(my_p_row,my_p_col)%block
       block_c => c%p(my_p_row,my_p_col)%block
       if (alpha==1.0_wp .and. beta==0.0_wp) then ! make common case fast....
#if (__parallel)
!MK! might be faster to do column by column (more contiguous)
!MK          do i=1,size_diag
!MK             global_row=indxl2g(i,descriptor(5),descriptor(my_p_row),&
!MK                  descriptor(7),n_p_rows)
!MK             if (global_row>size_diag) then
!MK                c%p(my_p_row,my_p_col)%block(i:,:)= 0.0_wp
!MK                exit
!MK             end if
!MK             c%p(my_p_row,my_p_col)%block(i,:)=&
!MK                  diagonal(global_row)*b%p(my_p_row,my_p_col)%block(i,:)
!MK          end do
#else
          do j=1,size(block_b,2)
             do i=1,min(size_diag,size(block_b,1))
                block_c(i,j)=diagonal(i)*&
                     block_b(i,j)
             end do
          end do
          if (size(block_b)>size(diagonal)) then
             block_c((size_diag+1):,:)= 0.0_wp
          end if
#endif
       else
          ! might be faster to do column by column (more contiguous)
          do i=1,size_diag
#if (__parallel)
!MK             global_row=indxl2g(i,descriptor(5),descriptor(my_p_row),&
!MK                  descriptor(7),n_p_rows)
#else
             global_row=i
#endif
             if (global_row>size_diag) then
                block_c(i,:)=&
                     al*diagonal(global_row)*block_b(i,:)+&
                     be*block_c(i,:)
                exit
             end if
             block_c(i,:)=&
                  al*diagonal(global_row)*block_b(i,:)+&
                  be*block_c(i,:)
          end do
       end if
    end if
    call timestop(0.0_wp,handle)
  end subroutine blacs_diag_mult
!***************************************************************************

END MODULE qs_blacs
