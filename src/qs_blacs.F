!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999-2002  CP2K developers group                            !
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

  USE kinds, ONLY: int_size,&
                   wp => dp,&
                   wp_size => dp_size

  USE global_types,        ONLY: global_environment_type
  USE mathlib,             ONLY: symmetrize_matrix
  USE util,                ONLY: sort
  USE sparse_matrix_types, ONLY: first_block_node,&
                                 get_block_node,&
                                 get_matrix_info,&
                                 next_block_node,&
                                 real_block_node_type,&
                                 real_matrix_type, &
                                 sparse_times_replicated, &
                                 sparse_plus_rep_rept
  USE message_passing,     ONLY: mp_bcast,&
                                 mp_max,&
                                 mp_sum,&
                                 mp_sync
  USE string_utilities,    ONLY: compress
  USE termination,         ONLY: stop_memory,&
                                 stop_program
  USE timings,             ONLY: timeset,&
                                 timestop
  USE timesl,              ONLY: cputime

  use cp_error_handling, only: cp_assert, cp_error_message, cp_error_type,&
       cp_assertion_failed
  use cp_log_handling, only: cp_to_string, cp_warning_level, cp_failure_level
  use cp_matrix_utils, only: cp_init, cp_dealloc_ref, cp_next, cp_get, &
       cp_matrix_block_iterator
  use cp_mpi_env, only: cp_mpi_env_type, cp_mpi_get_from_glob

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

  ! the default blacs block sizes
  ! consider using #ifdefs to give them the optimal values
  INTEGER, PARAMETER :: optimal_blacs_col_block_size=32
  INTEGER, PARAMETER :: optimal_blacs_row_block_size=32

! *** Public parameters ***
  PUBLIC :: optimal_blacs_col_block_size,&
       optimal_blacs_row_block_size

! *** Public data types ***

  PUBLIC :: blacs_matrix_p_type,&
            blacs_matrix_type

! *** Public subroutines ***
! at most 39 continuation lines :-)
  PUBLIC :: allocate_blacs_matrix, blacs_add,&
            blacs_add_to_element, blacs_gemm,&
            blacs_vector_symv, copy_vector_to_blacs, &
            copy_blacs_to_vector, blacs_find_vector_distribution,&
            blacs_daxpy,  blacs_get_element,&
            blacs_maxval, blacs_set_all,&
            blacs_set_element, blacs_syevd,&
            blacs_syevx, blacs_symm,&
            blacs_syrk, blacs_trace,&
            copy_blacs_to_blacs_matrix, copy_blacs_to_full_matrix,&
            copy_blacs_to_sparse_matrix, copy_sparse_to_blacs_matrix,&
            blacs_replicated_copy,  deallocate_blacs_matrix,&
            finish_blacs, get_blacs_info,&
            get_blacs_matrix_info, power_blacs_matrix,&
            read_blacs_matrix, replicate_blacs_matrix,&
            start_blacs, symmetrise_blacs_matrix,&
            write_blacs_matrix, blacs_diag_mult,&
            add_blacs_to_block_diag_sm,  blacs_cholesky_decompose, &
            blacs_cholesky_invert,  blacs_cholesky_reduce, &
            blacs_cholesky_restore,  blacs_jacobi_davidson, &
            blacs_column_copy,  blacs_make_basis, &
            blacs_init_random,  blacs_schur_product, &
            blacs_transpose,  blacs_scale_matrix, &
            blacs_set_local_block,  sparse_times_blacs, &
            sparse_plus_blacs_blacst, allocate_blacs_matrix_vect,&
            deallocate_blacs_matrix_vect, blacs_get_submatrix, &
            blacs_set_submatrix, blacs_scale_and_d

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_blacs_matrix(new_matrix,nrow_global,ncol_global,&
                                   nrow_block,ncol_block,name,globenv)

!   Purpose: Allocate a new distributed BLACS matrix.

!   History: - Creation (23.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: new_matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=*), INTENT(IN)              :: name
    INTEGER, INTENT(IN)                       :: ncol_block,&
                                                 ncol_global,nrow_block,&
                                                 nrow_global

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_blacs_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    CHARACTER(LEN=40) :: message
    INTEGER           :: context,group,ierror,ipcol,ipe,iprow,istat,mype,&
                         mypcol,myprow,ncol_local,npcol,npe,nprow,nrow_local,&
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
    context = globenv%context

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

    IF (ionode .and. .false.) THEN
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
    CALL DSCAL(nrow_local*ncol_local,0.0_wp,new_matrix%p(myprow,mypcol)%block,1)
    !new_matrix%p(myprow,mypcol)%block(:,:) = 0.0_wp

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
    CALL DSCAL(nrow_global*ncol_global,0.0_wp,new_matrix%p(source,source)%block,1)
    !new_matrix%p(source,source)%block(:,:) = 0.0_wp

#endif
  END SUBROUTINE allocate_blacs_matrix

! *****************************************************************************
!!****f* qs_blacs/allocate_blacs_matrix_vect [1.0] *
!!
!!   NAME
!!     allocate_blacs_matrix_vect
!!
!!   FUNCTION
!!     allocates an array of identical blacs matrixes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - blacs_matrixes: a pointer that will point to the array of matrixes
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
subroutine allocate_blacs_matrix_vect(blacs_matrixes,n_matrixes,&
     nrow_global,ncol_global,globenv,name,&
     nrow_block,ncol_block,error)
  type(blacs_matrix_p_type), dimension(:), pointer :: blacs_matrixes
  type(global_environment_type), intent(in):: globenv
  character(len=*), intent(in) :: name
  integer, intent(in) :: n_matrixes, nrow_global, ncol_global
  integer, intent(in), optional :: nrow_block, ncol_block
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='blacs_matrixes',&
        routineP=moduleN//':'//routineN
  integer :: i, stat, row_bl, col_bl
  failure=.false.
  row_bl=optimal_blacs_row_block_size
  col_bl=optimal_blacs_col_block_size
  if (present(nrow_block)) row_bl=nrow_block
  if (present(ncol_block)) col_bl=ncol_block
  
  allocate(blacs_matrixes(n_matrixes),stat=stat)
  CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  if (.not. failure) then
     do i=1,n_matrixes
        call allocate_blacs_matrix(blacs_matrixes(i)%blacs_matrix,&
             nrow_global, ncol_global,row_bl,col_bl,name//cp_to_string(i),&
             globenv)
     end do
  end if
end subroutine allocate_blacs_matrix_vect
!***************************************************************************

  SUBROUTINE blacs_set_local_block(matrix,block,globenv)
! *****************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp),DIMENSION(:,:), POINTER          :: block,mblock

    INTEGER :: context,nprow,npcol,myprow,mypcol,&
               nrow_local,ncol_local,source,nelements

    source=globenv%source
    context = globenv%context

#if defined(__parallel)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    nrow_local  = matrix%p(myprow,mypcol)%nrow_local
    ncol_local  = matrix%p(myprow,mypcol)%ncol_local
    nelements = nrow_local*ncol_local
    mblock => matrix%p(myprow,mypcol)%block
#else
    nrow_local  = matrix%p(source,source)%nrow_local
    ncol_local  = matrix%p(source,source)%ncol_local
    nelements = nrow_local*ncol_local
    mblock => matrix%p(source,source)%block
#endif

    CALL DCOPY(nelements,block(1,1),1,mblock(1,1),1)

  END SUBROUTINE blacs_set_local_block


! *****************************************************************************

  SUBROUTINE blacs_add(alpha,matrix_a,beta,matrix_b,globenv)

!   Purpose: Scale and add two BLACS matrices (a <- alpha*a + beta*b).

!   History: - Creation (11.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix_a,matrix_b
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(IN)                      :: alpha,beta

!   *** Local variables ***

    INTEGER :: context,handle,mypcol,myprow,npcol,nprow,source

    REAL(wp), DIMENSION(:,:), POINTER :: a,b

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_add","I","",handle)

    source = globenv%source
    context = globenv%context
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
  SUBROUTINE blacs_scale_matrix(matrixa,scaling,globenv,kmax)
! ugly routine to perform a scaling of every column in a matrix
! ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrixa
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), DIMENSION(:), POINTER           :: scaling
    INTEGER, OPTIONAL :: kmax

!   *** Local variables ***

    INTEGER :: context,icol_global,icol_local,ipcol,iprow,irow_local,mypcol,&
               myprow,npcol,nprow,source,i,k,n

    INTEGER, DIMENSION(9) :: desca
    REAL(wp), DIMENSION(:,:), POINTER :: a

!   ---------------------------------------------------------------------------

    source = globenv%source
    context = globenv%context
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    a => matrixa%p(myprow,mypcol)%block
    desca(:) = matrixa%descriptor(:)
    if (present(kmax)) then
      k =  kmax
    else
      k =  matrixa%ncol_global
    endif
    n =  matrixa%p(myprow,mypcol)%nrow_local ! mind this

    DO icol_global=1,k
       CALL infog2l(1,icol_global,desca,nprow,npcol,myprow,mypcol,&
                 irow_local,icol_local,iprow,ipcol)
       IF ((ipcol == mypcol)) THEN
           CALL DSCAL(n,scaling(icol_global),a(1,icol_local),1)
       END IF
    ENDDO


#else

    a => matrixa%p(source,source)%block
    n =  matrixa%nrow_global
    k =  matrixa%ncol_global
    do i=1,k
      CALL DSCAL(n,scaling(i),a(1,i),1)
    enddo

#endif
  END SUBROUTINE blacs_scale_matrix

! *****************************************************************************

  SUBROUTINE blacs_add_to_element(matrix,irow_global,icol_global,alpha,globenv)

!   Purpose: Add alpha to the BLACS matrix element
!            matrix(irow_global,icol_global).

!   History: - Creation (08.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(IN)                      :: alpha
    INTEGER, INTENT(IN)                       :: icol_global,&
                                                 irow_global

!   *** Local variables ***

    INTEGER :: context,icol_local,ipcol,iprow,irow_local,mypcol,myprow,npcol,&
               nprow,source

    INTEGER, DIMENSION(9) :: desca

    REAL(wp), DIMENSION(:,:), POINTER :: a

!   ---------------------------------------------------------------------------

    source = globenv%source
    context = globenv%context
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

! *****************************************************************************
  SUBROUTINE sparse_times_blacs(sparse_matrix,v_in,v_out,ncol, &
       globenv,alpha,beta,mpi_env,block_owner)
    TYPE(real_matrix_type), POINTER   :: sparse_matrix
    TYPE(blacs_matrix_type) , POINTER :: v_in
    TYPE(blacs_matrix_type) , POINTER :: v_out
    INTEGER, INTENT(IN)               :: ncol
    TYPE(global_environment_type), INTENT(IN),optional :: globenv
    type(cp_mpi_env_type), intent(in),optional :: mpi_env
    real(kind=wp), intent(in), optional :: alpha, beta
    integer, intent(in), dimension(:,:), optional :: block_owner

!   *** Local parameters ***
    INTEGER, PARAMETER          :: blocksize = 100
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE sparse_times_blacs (MODULE qs_blacs)"

!   *** Local variables ***
    REAL(wp), DIMENSION(:,:), POINTER :: rep_v_in,rep_v_out
    integer i,nrow_in,nrow_out,ncol_block,handle
    type(cp_mpi_env_type) :: mpienv
!   ---------------------------------------------------------------------------
    CALL timeset("sparse_times_blacs","I","",handle)
    if (present(mpi_env)) then
       mpienv=mpi_env
    else if (present(globenv)) then
       mpienv=cp_mpi_get_from_glob(globenv)
    else
       call cp_assert(.false.,cp_failure_level,cp_assertion_failed,routine,&
            "either globenv or mpi_env must be given")
    end if
    nrow_in=v_in%nrow_global
    nrow_out=v_out%nrow_global
    allocate(rep_v_in(blocksize,nrow_in))
    allocate(rep_v_out(blocksize,nrow_out))

    DO i=1,ncol,blocksize
       ncol_block=MIN(ncol-i+1,blocksize)
       call blacs_get_submatrix(v_in,rep_v_in,start_col=i,n_cols=ncol_block,&
            transpose_target_m=.true.,mpi_env=mpienv)
       !FM call blacs_replicated_copy(v_in,rep_v_in,i,ncol_block,blocksize,"B2R",globenv)
       call sparse_times_replicated(sparse_matrix,rep_v_in,rep_v_out,ncol_block, &
                                    blocksize,mpienv, block_owner)
       !FM call blacs_replicated_copy(v_out,rep_v_out,i,ncol_block,blocksize,"R2B",globenv)
       call blacs_set_submatrix(v_out,rep_v_out,&
            mpi_env=mpienv,&
            start_col=i,n_cols=ncol_block,transpose_new_val=.true.,&
            alpha=alpha,beta=beta)
            
    ENDDO

    deallocate(rep_v_in,rep_v_out)
    CALL timestop(0.0_wp,handle)
  END SUBROUTINE sparse_times_blacs

! *****************************************************************************
! computes sparse=sparse+alpha*v*g^T 
! only v*v^T tested
! *****************************************************************************
  SUBROUTINE sparse_plus_blacs_blacst(sparse_matrix,matrix_v,matrix_g,ncol,&
                                      globenv,first_col,alpha)
    TYPE(real_matrix_type), POINTER   :: sparse_matrix
    TYPE(blacs_matrix_type) , POINTER :: matrix_v
    TYPE(blacs_matrix_type) , POINTER, OPTIONAL :: matrix_g
    INTEGER, INTENT(IN)               :: ncol
    REAL(wp), OPTIONAL, INTENT(IN)    :: alpha
    INTEGER, OPTIONAL, INTENT(IN)     :: first_col
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***
    INTEGER, PARAMETER          :: blocksize = 100
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE sparse_plus_blacs_blacst (MODULE qs_blacs)"

!   *** Local variables ***
    REAL(wp), DIMENSION(:,:), POINTER :: rep_v,rep_g
    REAL(wp) :: prefactor
    INTEGER col1,i,nrow,ncol_block,handle
    LOGICAL g_present
!   ---------------------------------------------------------------------------
    CALL timeset("sparse_plus_blacs_blacst","I","",handle)

    g_present=.false.
    IF (PRESENT(matrix_g)) g_present=.true.

    nrow=matrix_v%nrow_global
    IF (g_present) THEN
       allocate(rep_v(blocksize,nrow))
       allocate(rep_g(blocksize,nrow))
    ELSE
       allocate(rep_v(blocksize,nrow))
       rep_g=>rep_v
    ENDIF

    IF (PRESENT(alpha)) THEN
      prefactor = alpha
    ELSE
      prefactor = 1.0_wp
    END IF

    IF (PRESENT(first_col)) THEN
      col1 = first_col
    ELSE
      col1 = 1
    END IF

    DO i=col1,ncol,blocksize
       ncol_block=MIN(ncol-i+1,blocksize)
       call blacs_replicated_copy(matrix_v,rep_v,i,ncol_block, &
                                        blocksize,"B2R",globenv)
       IF (g_present) THEN
         call blacs_replicated_copy(matrix_g,rep_g,i,ncol_block, &
                                          blocksize,"B2R",globenv)
       ENDIF
       call sparse_plus_rep_rept(sparse_matrix,rep_v,rep_g,ncol_block, &
                                    blocksize,globenv,alpha=prefactor)
    ENDDO

    IF (g_present) THEN
       deallocate(rep_v,rep_g)
    ELSE
       deallocate(rep_v)
    ENDIF

    CALL timestop(0.0_wp,handle)
  END SUBROUTINE sparse_plus_blacs_blacst
! *****************************************************************************


! *****************************************************************************

  SUBROUTINE blacs_gemm(transa,transb,m,n,k,alpha,matrix_a,matrix_b,beta,&
                        matrix_c,globenv)

!   Purpose: BLACS interface to the BLAS routine dgemm.

!   History: - Creation (07.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix_a,matrix_b,matrix_c
    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=1), INTENT(IN)              :: transa,transb
    REAL(wp), INTENT(IN)                      :: alpha,beta
    INTEGER, INTENT(IN)                       :: k,m,n

!   *** Local variables ***

    INTEGER :: context,handle,lda,ldb,ldc,mypcol,myprow,npcol,nprow,source

    INTEGER, DIMENSION(9) :: desca,descb,descc

    REAL(wp), DIMENSION(:,:), POINTER :: a,b,c
    REAL(wp) :: flops

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_gemm","I","Gflops",handle)

    source = globenv%source
    context = globenv%context
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    a => matrix_a%p(myprow,mypcol)%block
    desca(:) = matrix_a%descriptor(:)
    b => matrix_b%p(myprow,mypcol)%block
    descb(:) = matrix_b%descriptor(:)
    c => matrix_c%p(myprow,mypcol)%block
    descc(:) = matrix_c%descriptor(:)
    flops=2.0_wp*(1E-9_wp*m*n*k)/(nprow*npcol)
    CALL pdgemm(transa,transb,m,n,k,alpha,a(1,1),1,1,desca,b(1,1),1,1, &
                descb,beta,c(1,1),1,1,&
                descc)

#else

    a => matrix_a%p(source,source)%block
    b => matrix_b%p(source,source)%block
    c => matrix_c%p(source,source)%block

    lda = matrix_a%nrow_global
    ldb = matrix_b%nrow_global
    ldc = matrix_c%nrow_global

    flops=2.0_wp*(1E-9_wp*m*n*k)
    CALL dgemm(transa,transb,m,n,k,alpha,a(1,1),lda,b(1,1),ldb,beta,c(1,1),ldc)

#endif
    CALL timestop(flops,handle)

  END SUBROUTINE blacs_gemm

  SUBROUTINE blacs_daxpy(matrix_a,c1,alpha,c2,globenv)
    TYPE(blacs_matrix_type), POINTER       :: matrix_a
    TYPE(global_environment_type)          :: globenv
    INTEGER, INTENT(IN)                    :: c1,c2
    REAL(wp), INTENT(IN)                   :: alpha

    INTEGER :: context,nprow,npcol,myprow,mypcol,n,source
    REAL(wp), DIMENSION(:,:), POINTER :: a
    INTEGER, DIMENSION(9) :: desca

    n = matrix_a%nrow_global
    source = globenv%source
    context = globenv%context

#if defined(__parallel)
    desca(:)=matrix_a%descriptor(:)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    a => matrix_a%p(myprow,mypcol)%block
    CALL PDAXPY(n,alpha,a(1,1),1,c1,desca,1,a(1,1),1,c2,desca,1)
#else
    a => matrix_a%p(source,source)%block
    CALL DAXPY(n,alpha,a(1,c1),1,a(1,c2),1)
#endif

  END SUBROUTINE

!
! decides what fraction of the vector has to be on the local pe,
! giving each pe about the same amount of data (that is nloc elements),
! in this way we can map more easily to the blacs distribution
! the size of the vector is determined from matrix_a
! columns of matrix_a should represent a vector
! notice that it is the blacs layout that defines how the vector is distributed
! and what element of the vector actually is where

  SUBROUTINE blacs_find_vector_distribution(matrix_a,nloc,globenv)
    TYPE(blacs_matrix_type), POINTER       :: matrix_a
    TYPE(global_environment_type), INTENT(IN)      :: globenv
    INTEGER, INTENT(OUT)                   :: nloc
    INTEGER :: context,nprow,npcol,myprow,mypcol,n,nloc_blacs

    n = matrix_a%nrow_global
    context = globenv%context

#if defined(__parallel)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    nloc_blacs = matrix_a%p(myprow,mypcol)%nrow_local
    nloc = nloc_blacs / npcol
    if (mypcol.eq.0) then
         nloc = nloc_blacs - (npcol - 1) * nloc
    endif
#else
    nloc   = n
#endif

  END SUBROUTINE blacs_find_vector_distribution
!
! copy a blacs vector to a distributed vector
! has obviously to be consistent with find_vector_distribution
! see notes at blacs_vector_symv
! use_col_nr specifies the column of the matrix_a to use.
! this is the global col index
!
  SUBROUTINE copy_blacs_to_vector(matrix_a,use_col_nr,v,globenv)
    TYPE(blacs_matrix_type), POINTER       :: matrix_a
    TYPE(global_environment_type)          :: globenv
    INTEGER, INTENT(IN)                    :: use_col_nr
    REAL(wp), POINTER, DIMENSION(:)        :: v
    INTEGER :: context,nprow,npcol,myprow,mypcol,n,nloc_blacs,source,nloc,nloc0,i
    INTEGER :: irow_global,icol_global,irow_local,icol_local,iprow,ipcol
    INTEGER :: ndata,offset
    INTEGER, DIMENSION(9) :: desca

    n = matrix_a%nrow_global
    source = globenv%source
    context = globenv%context

#if defined (__parallel)
    desca(:) = matrix_a%descriptor(:)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    irow_global=1
    icol_global=use_col_nr
    CALL infog2l(irow_global,icol_global,desca,nprow,npcol,myprow,mypcol,&
                 irow_local,icol_local,iprow,ipcol)

    nloc_blacs = matrix_a%p(myprow,mypcol)%nrow_local
    nloc  = nloc_blacs / npcol                  ! all cols
    nloc0 = nloc_blacs - (npcol - 1) * nloc     ! col zero

    do i=0,npcol-1

       ndata=nloc
       if (i.eq.0) ndata=nloc0

       if (mypcol.eq.ipcol) then ! get your part of the data
          offset=nloc0+1+(i-1)*nloc
          if (i.eq.0) offset=1

          if (mypcol.eq.i) then
             CALL DCOPY(ndata,matrix_a%p(myprow,mypcol)%block(offset,icol_local),1,&
                        v,1)
          else
             CALL DGESD2D(context,ndata,1,   &
             matrix_a%p(myprow,mypcol)%block(offset,icol_local),  &
             ndata,myprow,i)
          endif
       else
          if (mypcol.eq.i) then
             CALL DGERV2D(context,ndata,1,   &
             v, &
             ndata,myprow,ipcol)
          endif
       endif
    enddo
#else
    nloc = n
    if (matrix_a%p(source,source)%ncol_local.lt.use_col_nr) then
       call stop_program("copy_blacs_to_vector","hmm",globenv)
    endif
    CALL DCOPY(nloc,matrix_a%p(source,source)%block(1,use_col_nr),1,v,1)
#endif

  END SUBROUTINE copy_blacs_to_vector

!
! copy a full distributed vector into a blacs vector
! has obviously to be consistent with find_vector_distribution
! see notes at blacs_vector_symv
!
  SUBROUTINE copy_vector_to_blacs(v,matrix_a,use_col_nr,globenv)
    TYPE(blacs_matrix_type), POINTER       :: matrix_a
    TYPE(global_environment_type)          :: globenv
    INTEGER, INTENT(IN)                    :: use_col_nr
    REAL(wp), POINTER, DIMENSION(:)        :: v
    INTEGER :: context,nprow,npcol,myprow,mypcol,n,nloc_blacs,source,nloc,nloc0,i
    INTEGER :: irow_global,icol_global,irow_local,icol_local,iprow,ipcol
    INTEGER :: ndata,offset
    INTEGER, DIMENSION(9) :: desca

    n = matrix_a%nrow_global
    source = globenv%source
    context = globenv%context

#if defined (__parallel)
    desca(:) = matrix_a%descriptor(:)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    irow_global=1
    icol_global=use_col_nr
    CALL infog2l(irow_global,icol_global,desca,nprow,npcol,myprow,mypcol,&
                 irow_local,icol_local,iprow,ipcol)


    nloc_blacs = matrix_a%p(myprow,mypcol)%nrow_local
    nloc  = nloc_blacs / npcol                  ! all cols
    nloc0 = nloc_blacs - (npcol - 1) * nloc     ! col zero

    do i=0,npcol-1

       ndata=nloc
       if (i.eq.0) ndata=nloc0

       if (mypcol.eq.ipcol) then ! get your part of the data
          offset=nloc0+1+(i-1)*nloc
          if (i.eq.0) offset=1

          if (mypcol.eq.i) then
             CALL DCOPY(ndata,v,1,matrix_a%p(myprow,mypcol)%block(offset,icol_local),1)
          else
             CALL DGERV2D(context,ndata,1,   &
             matrix_a%p(myprow,mypcol)%block(offset,icol_local),  &
             ndata,myprow,i)
          endif
       else
          if (mypcol.eq.i) then
             CALL DGESD2D(context,ndata,1,   &
             v, &
             ndata,myprow,ipcol)
          endif
       endif
    enddo
#else
    nloc = n
    if (matrix_a%p(source,source)%ncol_local.lt.use_col_nr) then
       call stop_program("copy_vector_to_blacs","hmm",globenv)
    endif
    CALL DCOPY(nloc,v,1,matrix_a%p(source,source)%block(1,use_col_nr),1)
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

  SUBROUTINE blacs_vector_symv(matrix_a,vin,vout,matrix_buf,globenv)
! *****************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix_a,matrix_buf
    TYPE(global_environment_type)             :: globenv
    REAL(wp), DIMENSION(:), POINTER           :: vin,vout

!   *** Local variables ***

    INTEGER :: context,handle,source,n,inc
    REAL(wp), DIMENSION(:,:), POINTER :: a,b
    REAL(wp) :: alpha,beta
    INTEGER :: nprow,npcol,myprow,mypcol
    INTEGER, DIMENSION(9) :: desca,descb


!   ---------------------------------------------------------------------------

    CALL timeset("blacs_vector_symv","I","",handle)

    source = globenv%source
    context = globenv%context
    alpha=1.0D0
    beta=0.0D0
    inc=1

    n   = matrix_a%nrow_global
#if defined(__parallel)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    a => matrix_a%p(myprow,mypcol)%block
    desca(:) = matrix_a%descriptor(:)
    b => matrix_buf%p(myprow,mypcol)%block
    descb(:) = matrix_buf%descriptor(:)
    call copy_vector_to_blacs(vin,matrix_buf,1,globenv)
    call pdsymv('U',n,alpha,a(1,1),1,1,desca,b(1,1),1,1,descb,inc,beta,b(1,1),1,2,descb,inc)
    call copy_blacs_to_vector(matrix_buf,2,vout,globenv)
#else
    a => matrix_a%p(source,source)%block
    call dsymv('U',n,alpha,a(1,1),n,vin(1),inc,beta,vout(1),inc)
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
!    (a random initial basis can be obtained by blacs_init_random)
!
! b) if pmatrix is present, it is assumed that this is a projector of rank nmo
! and a basis for this projector is found i.e. pmatrix=v*v^T
!
!
! *****************************************************************************
  SUBROUTINE blacs_make_basis(vmatrix,ncol,globenv,ortho,otype,pmatrix,work)
! *****************************************************************************

    TYPE(blacs_matrix_type), POINTER                   :: vmatrix
    TYPE(blacs_matrix_type), POINTER, OPTIONAL         :: pmatrix,work,ortho
    TYPE(global_environment_type)                      :: globenv
    INTEGER, INTENT(IN)                                :: ncol
    CHARACTER ( LEN = * ), INTENT ( IN ),OPTIONAL      :: otype

!   *** Local variables ***

    INTEGER :: context,handle,n,source,info
    REAL(wp), DIMENSION(:,:), POINTER :: a,b,vv,w,p,o
    REAL(wp) :: rone,rzero
    INTEGER :: nprow,npcol,myprow,mypcol,nrow_block,ncol_block,i
    INTEGER, DIMENSION(9) :: desca,descb,descvv,descw,descp,desco
    TYPE(blacs_matrix_type), POINTER           :: overlap_vv
    logical :: found_type,use_pmatrix,use_cholesky,use_sv,use_simple



!   ---------------------------------------------------------------------------

    CALL timeset("blacs_make_basis","I","",handle)

! check input ...

    IF (present(pmatrix)) THEN
       IF (.not. present(work)) THEN
          call stop_program("make_basis","pmatrix needs work space !?")
       ENDIF
       use_pmatrix=.true.
    ELSE
       use_pmatrix=.false.
    ENDIF

    IF (present(ortho)) THEN
       IF (.not. present(otype)) THEN
           call stop_program("make_basis","ortho needs type !?")
       ENDIF
       found_type=.false.
       use_sv=.false.
       use_cholesky=.false.
       IF (otype .eq. "CHOLESKY") THEN
          found_type=.true.
          use_cholesky=.true.
       ENDIF
       IF (otype .eq. "SV") THEN
          found_type=.true.
          use_sv=.true.
       ENDIF
       IF (.not. found_type) call stop_program("make_basis","sorry wrong type")
    ELSE
       use_simple=.true.
    ENDIF

    source = globenv%source
    context = globenv%context
    rone=1.0D0
    rzero=0.0D0
    n   = vmatrix%nrow_global
    CALL get_blacs_matrix_info(matrix=vmatrix,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block)

    CALL allocate_blacs_matrix(new_matrix=overlap_vv,&
                               nrow_global=ncol,&
                               ncol_global=ncol,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="overlap_vv",&
                               globenv=globenv)

#if defined(__parallel)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    a => vmatrix%p(myprow,mypcol)%block
    desca(:) = vmatrix%descriptor(:)
    vv => overlap_vv%p(myprow,mypcol)%block
    descvv(:) = overlap_vv%descriptor(:)
#else
    a => vmatrix%p(source,source)%block
    vv => overlap_vv%p(source,source)%block
#endif
    IF (use_pmatrix) THEN
       CALL blacs_init_random(work,ncol,globenv)
#if defined(__parallel)
       p => pmatrix%p(myprow,mypcol)%block
       descp(:) = pmatrix%descriptor(:)
       w => work%p(myprow,mypcol)%block
       descw(:) = work%descriptor(:)
       CALL PDSYMM('L','U',n,ncol,rone,p(1,1),1,1,descp,w(1,1),1,1,descw,rzero,a(1,1),1,1,desca)
       IF (use_cholesky) THEN
         o => ortho%p(myprow,mypcol)%block
         desco(:) = ortho%descriptor(:)
         CALL PDTRMM('L','U','N','N',n,ncol,rone,o(1,1),1,1,desco,a(1,1),1,1,desca)
       ENDIF
#else
       p => pmatrix%p(source,source)%block
       w => work%p(source,source)%block
       ! U*P*U^T*V
       CALL DSYMM('L','U',n,ncol,rone,p(1,1),n,w(1,1),n,rzero,a(1,1),n)
       IF (use_cholesky) THEN
         o => ortho%p(source,source)%block
         CALL DTRMM('L','U','N','N',n,ncol,rone,o(1,1),n,a(1,1),n)
       ENDIF
#endif
    ELSE
       ! U*V should be used
       IF (use_cholesky) THEN
#if defined(__parallel)
         o => ortho%p(myprow,mypcol)%block
         desco(:) = ortho%descriptor(:)
         CALL PDTRMM('L','U','N','N',n,ncol,rone,o(1,1),1,1,desco,a(1,1),1,1,desca)
#else
         o => ortho%p(source,source)%block
         CALL DTRMM('L','U','N','N',n,ncol,rone,o(1,1),n,a(1,1),n)
#endif
       ENDIF
    ENDIF

    IF (use_sv) THEN
      CALL blacs_GEMM('T','N',ncol,ncol,n,rone,vmatrix,ortho,rzero, &
                                          overlap_vv,globenv)
      CALL blacs_cholesky_decompose(overlap_vv,globenv)
#if defined(__parallel)
      CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
      a => vmatrix%p(myprow,mypcol)%block
      desca(:) = vmatrix%descriptor(:)
      vv => overlap_vv%p(myprow,mypcol)%block
      descvv(:) = overlap_vv%descriptor(:)
      o => ortho%p(myprow,mypcol)%block
      desco(:) = ortho%descriptor(:)
      CALL pdtrsm('R','U','N','N',n,ncol,rone,vv(1,1),1,1,descvv,a(1,1),1,1,desca)
      CALL pdtrsm('R','U','N','N',n,ncol,rone,vv(1,1),1,1,descvv,o(1,1),1,1,desco)
#else
      o => ortho%p(source,source)%block
      a => vmatrix%p(source,source)%block
      vv => overlap_vv%p(source,source)%block
      CALL dtrsm('R','U','N','N',n,ncol,rone,vv(1,1),ncol,a(1,1),n)
      CALL dtrsm('R','U','N','N',n,ncol,rone,vv(1,1),ncol,o(1,1),n)
#endif
    ELSE
      ! ortho (U*V) inner product should be one
#if defined(__parallel)
      CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
      a => vmatrix%p(myprow,mypcol)%block
      desca(:) = vmatrix%descriptor(:)
      vv => overlap_vv%p(myprow,mypcol)%block
      descvv(:) = overlap_vv%descriptor(:)
      CALL PDSYRK('U','T',ncol,n,rone,a(1,1),1,1,desca,rzero,vv(1,1),1,1,descvv,ncol)
      CALL PDPOTRF('U',ncol,vv(1,1),1,1,descvv,info)
      if (info.ne.0) call stop_program("make_basis","DPOTRF degenerate vectors")
      CALL PDTRSM('R','U','N','N',n,ncol,rone,vv(1,1),1,1,descvv,a(1,1),1,1,desca)
#else
      a => vmatrix%p(source,source)%block
      vv => overlap_vv%p(source,source)%block
      CALL DSYRK('U','T',ncol,n,rone,a(1,1),n,rzero,vv(1,1),ncol)
      CALL DPOTRF('U',ncol,vv(1,1),ncol,info)
      if (info.ne.0) call stop_program("make_basis","DPOTRF degenerate vectors")
      CALL DTRSM('R','U','N','N',n,ncol,rone,vv(1,1),ncol,a(1,1),n)
#endif
    ENDIF

    IF (use_cholesky) THEN
#if defined(__parallel)
      CALL PDTRSM('L','U','N','N',n,ncol,rone,o(1,1),1,1,desco,a(1,1),1,1,desca)
#else
      CALL DTRSM('L','U','N','N',n,ncol,rone,o(1,1),n,a(1,1),n)
#endif
    ENDIF

    call deallocate_blacs_matrix(overlap_vv)

    CALL timestop(0.0_wp,handle)

    END SUBROUTINE blacs_make_basis

! *****************************************************************************
! init ncol vectors of a blacs matrix with random numbers
! *****************************************************************************
    SUBROUTINE BLACS_INIT_RANDOM(matrix,ncol,globenv)

    TYPE(blacs_matrix_type), POINTER           :: matrix
    TYPE(global_environment_type)              :: globenv
    INTEGER, INTENT(IN)                        :: ncol

    INTEGER :: context,handle,n,source,i
    REAL(wp), DIMENSION(:,:), POINTER :: a
    INTEGER :: nprow,npcol,myprow,mypcol,nrow_block,ncol_block
    INTEGER, DIMENSION(9) :: desca
    INTEGER, DIMENSION(4),SAVE :: ISEED
    LOGICAL, SAVE :: FIRST = .true.
    INTEGER :: irow_global,icol_global,irow_local,icol_local,iprow,ipcol

    CALL timeset("blacs_init_random","I","",handle)

    source = globenv%source
    context = globenv%context
    n   = matrix%nrow_global

#if defined(__parallel)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    IF (FIRST) THEN
       ISEED(1)=myprow
       ISEED(2)=mypcol
       ISEED(3)=2
       ISEED(4)=1 ! last one has to be odd
       FIRST=.false.
    ENDIF
    a => matrix%p(myprow,mypcol)%block
    desca(:)=matrix%descriptor(:)
    do i=1,ncol
     CALL infog2l(1,i,desca,nprow,npcol,myprow,mypcol,&
                 irow_local,icol_local,iprow,ipcol)
     IF (ipcol == mypcol) THEN
        CALL dlarnv(1,iseed,size(a,1),a(1,icol_local))
     ENDIF
    enddo
#else
    IF (FIRST) THEN
       ISEED(1)=4
       ISEED(2)=3
       ISEED(3)=2
       ISEED(4)=1 ! last one has to be odd
       FIRST=.false.
    ENDIF
    a => matrix%p(source,source)%block
    CALL dlarnv(1,iseed,n*ncol,a(1,1))
#endif

    END SUBROUTINE

! *****************************************************************************

  SUBROUTINE blacs_get_element(matrix,irow_global,icol_global,alpha,globenv)

!   Purpose: Get the BLACS matrix element (irow_global,icol_global).

!   History: - Creation (22.01.2002, Matthias Krack)
!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(OUT)                     :: alpha
    INTEGER, INTENT(IN)                       :: icol_global,&
                                                 irow_global

!   *** Local variables ***

    INTEGER :: context,icol_local,ipcol,iprow,irow_local,mypcol,myprow,npcol,&
               nprow,source

    INTEGER, DIMENSION(9) :: desca

    REAL(wp), DIMENSION(:,:), POINTER :: a

!   ---------------------------------------------------------------------------

    source = globenv%source
    context = globenv%context
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

  SUBROUTINE blacs_maxval(matrix,a_max,globenv)

!   Purpose: Get the maximum absolute element of a BLACS matrix.

!   History: - Creation (11.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(OUT)                     :: a_max

!   *** Local variables ***

    INTEGER :: context,group,handle,mypcol,myprow,npcol,nprow,source

    REAL(wp), DIMENSION(:,:), POINTER :: my_block

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_maxval","I","",handle)

    source = globenv%source
    group = globenv%group
    context = globenv%context
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

  SUBROUTINE blacs_set_all(matrix,alpha,globenv)

!   Purpose: Set the BLACS matrix elements to alpha.

!   History: - Creation (12.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(IN)                      :: alpha

!   *** Local variables ***

    INTEGER :: context,mypcol,myprow,npcol,nprow,source

!   ---------------------------------------------------------------------------

    source = globenv%source
    context = globenv%context
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    matrix%p(myprow,mypcol)%block(:,:) = alpha

#else

    matrix%p(source,source)%block(:,:) = alpha

#endif
  END SUBROUTINE blacs_set_all

! *****************************************************************************

  SUBROUTINE blacs_set_element(matrix,irow_global,icol_global,alpha,globenv)

!   Purpose: Set the BLACS matrix element (irow_global,icol_global) to alpha.

!   History: - Creation (08.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(IN)                      :: alpha
    INTEGER, INTENT(IN)                       :: icol_global,&
                                                 irow_global

!   *** Local variables ***

    INTEGER :: context,icol_local,ipcol,iprow,irow_local,mypcol,myprow,npcol,&
               nprow,source

    INTEGER, DIMENSION(9) :: desca

    REAL(wp), DIMENSION(:,:), POINTER :: a

!   ---------------------------------------------------------------------------

    source = globenv%source
    context = globenv%context
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
  SUBROUTINE blacs_replicated_copy(matrix_v,rep_v,firstcol,ncol,nblock,direction,globenv)
    TYPE(blacs_matrix_type) , POINTER              :: matrix_v
    REAL(wp), DIMENSION(:,:), POINTER              :: rep_v
    INTEGER, INTENT(IN)                            :: firstcol,ncol,nblock
    TYPE(global_environment_type), INTENT(IN)      :: globenv
    CHARACTER(LEN=3)                               :: direction
    !----
    REAL(wp), DIMENSION(:,:), POINTER :: blacs_block
    INTEGER                           :: i,j,nrow_global,ncol_global,ncol_block
    INTEGER                           :: source,myprow,mypcol,nprow,npcol,context
    INTEGER                           :: this_col,nrow_local,ncol_local,nrow_block,handle
    INTEGER, DIMENSION(:), POINTER                 :: row_indices,col_indices

    CALL timeset("blacs_replicated_copy","I","",handle)

    source =globenv%source 
    context=globenv%context

#if defined(__parallel)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    blacs_block=>matrix_v%p(myprow,mypcol)%block
#else
    blacs_block=>matrix_v%p(source,source)%block
#endif

    CALL get_blacs_matrix_info(matrix=matrix_v, &
                          nrow_global=nrow_global,ncol_global=ncol_global, &
                          nrow_block =nrow_block ,ncol_block =ncol_block , &
                          nrow_local =nrow_local ,ncol_local =ncol_local , &
                          row_indices=row_indices,col_indices=col_indices,globenv=globenv)
  
#if defined(__parallel) 
    if (direction.eq."B2R") then
       CALL dcopy(nblock*nrow_global,0.0_wp,0,rep_v(1,1),1)
    endif
#endif

    DO i=1,ncol_local
       this_col=col_indices(i)-firstcol+1
       if (this_col.ge.1 .and. this_col.le.ncol) then
          if (direction.eq."B2R") then
             do j=1,nrow_local
               rep_v(this_col,row_indices(j))=blacs_block(j,i) 
             enddo
          else
             do j=1,nrow_local
                blacs_block(j,i)=rep_v(this_col,row_indices(j))
             enddo
          endif
       endif
    ENDDO

#if defined(__parallel) 
    if (direction.eq."B2R") then
       call mp_sum(rep_v,globenv%group)  
    endif
#endif

    deallocate(row_indices,col_indices)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_replicated_copy

  !!****f* qs_blacs/blacs_submatrix_set [1.0] *
  !!
  !!   NAME
  !!     blacs_submatrix_set
  !!
  !!   FUNCTION
  !!     sets a submatrix of a blacs matrix
  !!      blacs_matrix(start_row:start_row+n_rows,start_col:start_col+n_cols)
  !!      = alpha*op(new_values)(1:n_rows,1:n_cols)+ beta 
  !!      * blacs_matrix(start_row:start_row+n_rows,start_col:start_col+n_cols)
  !!
  !!   NOTES
  !!     optimized for full column updates and alpha=1.0, beta=0.0
  !!     col_indices, and row_indices should be cached in some structure
  !!     common to all the "similar" blacs matrixes
  !!
  !!   INPUTS
  !!     - new_values: a replicated full matrix with the new values
  !!     - blacs_matrix: the blacs to change
  !!     - mpi_env: the mpi environment
  !!     - start_row: the starting row of b_matrix (defaults to 1)
  !!     - start_col: the starting col of b_matrix (defaults to 1)
  !!     - n_rows: the number of row to change in b (defaults to
  !!       size(op(new_values),1))
  !!     - n_cols: the number of columns to change in b (defaults to
  !!       size(op(new_values),2))
  !!     - alpha: rescaling factor for the new values (defaults to 1.0)
  !!     - beta: rescaling factor for the old values (defaults to 0.0)
  !!     - transpose_a: if new_values should be transposed: if true
  !!       op(new_values)=new_values^T, else op(new_values)=new_values
  !!       (defaults to false)
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Fawzi Mohamed
  !!
  !!   MODIFICATION HISTORY
  !!     07.2002 created borrowing from Joost's blacs_replicated_copy [fawzi]
  !!
  !!*** *********************************************************************
  subroutine blacs_set_submatrix(blacs_matrix,new_values,mpi_env, start_row,&
       start_col, n_rows, n_cols, alpha, beta, transpose_new_val, error)
    type(blacs_matrix_type) , pointer :: blacs_matrix
    real(wp), dimension(:,:), pointer :: new_values
    type(cp_mpi_env_type), intent(in)      :: mpi_env
    integer, intent(in), optional :: start_row, start_col,&
         n_rows, n_cols
    real(kind=wp), intent(in), optional :: alpha, beta
    logical, intent(in), optional :: transpose_new_val
    type(cp_error_type), intent(inout), optional :: error

    real(kind=wp), dimension(:,:), pointer :: blacs_block
    real(kind=wp) :: al,be
    integer :: i0,j0,ncol,nrow,i,j,nrow_global,ncol_global,ncol_block,&
         source,myprow,mypcol,nprow,npcol,context, &
         this_col,this_row,nrow_local,ncol_local,nrow_block,handle
    integer, dimension(:), pointer :: row_indices,col_indices
    character(len=*), parameter :: routineN='blacs_set_submatrix',&
         routineP=moduleN//':'//routineN
    logical :: tr_a
    al=1.0_wp; be=0.0_wp; i0=1; j0=1; tr_a=.false.

    call timeset(routineN//','//moduleN,"I","",handle)

    if (present(alpha)) al=alpha
    if (present(beta)) be=beta
    if (present(start_row)) i0=start_row
    if (present(start_col)) j0=start_col
    if (present(transpose_new_val)) tr_a=transpose_new_val
    if (tr_a) then
       nrow=size(new_values,2)
       ncol=size(new_values,1)
    else
       nrow=size(new_values,1)
       ncol=size(new_values,2)
    end if
    if (present(n_rows)) nrow=n_rows
    if (present(n_cols)) ncol=n_cols

    source =mpi_env%source
    context=blacs_matrix%context

#if defined(__parallel)
    call blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    blacs_block => blacs_matrix%p(myprow,mypcol)%block
#else
    blacs_block => blacs_matrix%p(source,source)%block
#endif

    call get_blacs_matrix_info(matrix=blacs_matrix, &
         nrow_global=nrow_global,ncol_global=ncol_global, &
         nrow_block =nrow_block ,ncol_block =ncol_block , &
         nrow_local =nrow_local ,ncol_local =ncol_local , &
         row_indices=row_indices,col_indices=col_indices,&
         mpi_env=mpi_env)


    if (al==1.0.and.be==0.0) then
       do j=1,ncol_local
          this_col=col_indices(j)-j0+1
          if (this_col.ge.1 .and. this_col.le.ncol) then
             if (tr_a) then
                if (i0==1.and.nrow_global==nrow) then
                   do i=1,nrow_local
                      blacs_block(i,j)=new_values(this_col,row_indices(i))
                   end do
                else
                   do i=1,nrow_local
                      this_row=row_indices(i)-i0+1
                      if (this_row>=1 .and. this_row<=nrow) then
                         blacs_block(i,j)=new_values(this_col,this_row)
                      end if
                   end do
                end if
             else
                if (i0==1.and.nrow_global==nrow) then
                   do i=1,nrow_local
                      blacs_block(i,j)=new_values(row_indices(i),this_col)
                   end do
                else
                   do i=1,nrow_local
                      this_row=row_indices(i)-i0+1
                      if (this_row>=1 .and. this_row<=nrow) then
                         blacs_block(i,j)=new_values(this_row,this_col)
                      end if
                   end do
                end if
             end if
          end if
       end do
    else
       do j=1,ncol_local
          this_col=col_indices(j)-j0+1
          if (this_col.ge.1 .and. this_col.le.ncol) then
             if (tr_a) then
                do i=1,nrow_local
                   this_row=row_indices(i)-i0+1
                   if (this_row>=1 .and. this_row<=nrow) then
                      blacs_block(i,j)=al*new_values(this_col,this_row)+&
                           be*blacs_block(i,j)
                   end if
                end do
             else
                do i=1,nrow_local
                   this_row=row_indices(i)-i0+1
                   if (this_row>=1 .and. this_row<=nrow) then
                      blacs_block(i,j)=al*new_values(this_row,this_col)+&
                           be*blacs_block(i,j)
                   end if
                end do
             end if
          end if
       end do
    end if
    deallocate(row_indices,col_indices)

    call timestop(0.0_wp,handle)

  end subroutine blacs_set_submatrix

  !!****f* qs_blacs/blacs_get_submatrix [1.0] *
  !!
  !!   NAME
  !!     blacs_get_submatrix
  !!
  !!   FUNCTION
  !!     gets a submatrix of a blacs matrix
  !!      op(target_m)(1:n_rows,1:n_cols)
  !!      =blacs_matrix(start_row:start_row+n_rows,start_col:start_col+n_cols)
  !!
  !!   NOTES
  !!     optimized for full column updates. Zeros out a little too much
  !!     of target_m
  !!     col_indices, and row_indices should be cached in some structure
  !!     common to all the "similar" blacs matrixes
  !!
  !!   INPUTS
  !!     - target_m: a replicated full matrix that will contain the result
  !!     - blacs_matrix: the blacs you want to get the info from
  !!     - start_row: the starting row of b_matrix (defaults to 1)
  !!     - start_col: the starting col of b_matrix (defaults to 1)
  !!     - n_rows: the number of row to change in b (defaults to
  !!       size(op(new_values),1))
  !!     - n_cols: the number of columns to change in b (defaults to
  !!       size(op(new_values),2))
  !!     - transpose_target_m: if target_m should be transposed: if true
  !!       op(target_m)=target_m^T, else op(target_m)=target_m
  !!       (defaults to false)
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Fawzi Mohamed
  !!
  !!   MODIFICATION HISTORY
  !!     07.2002 created borrowing from Joost's blacs_replicated_copy [fawzi]
  !!
  !!*** *********************************************************************
  subroutine blacs_get_submatrix(blacs_matrix,target_m,mpi_env, start_row,&
       start_col, n_rows, n_cols, transpose_target_m, error)
    type(blacs_matrix_type) , pointer :: blacs_matrix
    real(wp), dimension(:,:), pointer :: target_m
    type(cp_mpi_env_type), intent(in)      :: mpi_env
    integer, intent(in), optional :: start_row, start_col,&
         n_rows, n_cols
    logical, intent(in), optional :: transpose_target_m
    type(cp_error_type), intent(inout), optional :: error

    real(kind=wp), dimension(:,:), pointer :: blacs_block
    integer :: i0,j0,ncol,nrow,i,j,nrow_global,ncol_global,ncol_block,&
         source,myprow,mypcol,nprow,npcol,context, &
         this_col,this_row,nrow_local,ncol_local,nrow_block,handle
    integer, dimension(:), pointer :: row_indices,col_indices
    character(len=*), parameter :: routineN='blacs_get_submatrix',&
         routineP=moduleN//':'//routineN
    logical :: tr_a
    i0=1; j0=1; tr_a=.false.

    call timeset(routineN//','//moduleN,"I","",handle)

    if (present(start_row)) i0=start_row
    if (present(start_col)) j0=start_col
    if (present(transpose_target_m)) tr_a=transpose_target_m
    if (tr_a) then
       nrow=size(target_m,2)
       ncol=size(target_m,1)
    else
       nrow=size(target_m,1)
       ncol=size(target_m,2)
    end if
    if (present(n_rows)) nrow=n_rows
    if (present(n_cols)) ncol=n_cols

    source =mpi_env%source
    context=blacs_matrix%context

#if defined(__parallel)
    ! zero-out whole target_m 
    call dcopy(size(target_m,1)*size(target_m,2),0.0_wp,0,target_m(1,1),1)
    call blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    blacs_block => blacs_matrix%p(myprow,mypcol)%block
#else
    blacs_block => blacs_matrix%p(source,source)%block
#endif

    call get_blacs_matrix_info(matrix=blacs_matrix, &
         nrow_global=nrow_global,ncol_global=ncol_global, &
         nrow_block =nrow_block ,ncol_block =ncol_block , &
         nrow_local =nrow_local ,ncol_local =ncol_local , &
         row_indices=row_indices,col_indices=col_indices,&
         mpi_env=mpi_env)


    do j=1,ncol_local
       this_col=col_indices(j)-j0+1
       if (this_col.ge.1 .and. this_col.le.ncol) then
          if (tr_a) then
             if (i0==1.and.nrow_global==nrow) then
                do i=1,nrow_local
                   target_m(this_col,row_indices(i))=blacs_block(i,j)
                end do
             else
                do i=1,nrow_local
                   this_row=row_indices(i)-i0+1
                   if (this_row>=1 .and. this_row<=nrow) then
                      target_m(this_col,this_row)=blacs_block(i,j)
                   end if
                end do
             end if
          else
             if (i0==1.and.nrow_global==nrow) then
                do i=1,nrow_local
                   target_m(row_indices(i),this_col)=blacs_block(i,j)
                end do
             else
                do i=1,nrow_local
                   this_row=row_indices(i)-i0+1
                   if (this_row>=1 .and. this_row<=nrow) then
                      target_m(this_row,this_col)=blacs_block(i,j)
                   end if
                end do
             end if
          end if
       end if
    end do
    deallocate(row_indices,col_indices)

    call mp_sum(target_m,mpi_env%group)
    
    call timestop(0.0_wp,handle)

  end subroutine blacs_get_submatrix

!******************************************************************
! used to replace the cholesky decomposition by the inverse
!******************************************************************
  SUBROUTINE blacs_cholesky_invert(matrix,globenv)
   TYPE(blacs_matrix_type), POINTER           :: matrix

    TYPE(global_environment_type), INTENT(IN) :: globenv

    REAL(wp), DIMENSION(:,:), POINTER         :: a
    integer                                   :: context,info,handle
    INTEGER                                   :: source,n
    integer                                   :: nprow,npcol,myprow,mypcol
    INTEGER, DIMENSION(9)                     :: desca

    CALL timeset("blacs_cholesky_invert","I","",handle)

    source = globenv%source
    context = globenv%context
    n = matrix%nrow_global

#if defined(__parallel)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    a => matrix%p(myprow,mypcol)%block
    desca(:) = matrix%descriptor(:)
    CALL pdpotri('U',n,a(1,1),1,1,desca,info)

#else

    a => matrix%p(source,source)%block

    CALL dpotri('U',n,a(1,1),n,info)

#endif

    if (info.ne.0) then
      CALL stop_program(" blacs_cholesky_invert","Error dpotri",globenv)
    endif


    CALL timestop(0.0_wp,handle)

  END  SUBROUTINE blacs_cholesky_invert




!******************************************************************
! used to replace a symmetric positive def. matrix by its cholesky
! decomposition
!******************************************************************

  SUBROUTINE blacs_cholesky_decompose(matrix,globenv)
   TYPE(blacs_matrix_type), POINTER           :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv

    REAL(wp), DIMENSION(:,:), POINTER         :: a
    integer                                   :: context,info,handle
    INTEGER                                   :: output_unit,source,n
    LOGICAL                                   :: ionode
    integer                                   :: nprow,npcol,myprow,mypcol
    INTEGER, DIMENSION(9)                     :: desca

    CALL timeset("blacs_cholesky_decompose","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr
    source = globenv%source
    context = globenv%context
    n = matrix%nrow_global

#if defined(__parallel)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    a => matrix%p(myprow,mypcol)%block
    desca(:) = matrix%descriptor(:)
    CALL pdpotrf('U',n,a(1,1),1,1,desca,info)

#else

    a => matrix%p(source,source)%block

    CALL dpotrf('U',n,a(1,1),n,info)

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
! A,B -> inv(U^T)*A*inv(U),1 (AX=BX -> inv(U^T)*A*inv(U)*U*X=U*X hence evecs U*X)
!******************************************************************

  SUBROUTINE blacs_cholesky_reduce(matrix,matrixb,globenv)
   TYPE(blacs_matrix_type), POINTER           :: matrix,matrixb
    TYPE(global_environment_type), INTENT(IN) :: globenv

    REAL(wp), DIMENSION(:,:), POINTER         :: a,b
    integer                                   :: context,info,itype,handle
    INTEGER                                   :: output_unit,source,n
    LOGICAL                                   :: ionode
    integer                                   :: nprow,npcol,myprow,mypcol
    INTEGER, DIMENSION(9)                     :: desca,descb
    REAL(wp)                                  :: scale

    CALL timeset("blacs_cholesky_reduce","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr
    source = globenv%source
    context = globenv%context
    n = matrix%nrow_global
    itype =1

#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    a => matrix%p(myprow,mypcol)%block
    desca(:) = matrix%descriptor(:)
    b => matrixb%p(myprow,mypcol)%block
    descb(:) = matrixb%descriptor(:)

    CALL pdsygst(itype,'U',n,a(1,1),1,1,desca,b(1,1),1,1,descb,scale,info)
    IF (scale .ne. 1.0_wp) THEN
       ! this is supposed to be one in current version of lapack
       ! if not, eigenvalues have to be scaled by this number
       CALL stop_program(" blacs_cholesky_reduce","scale not equal 1",globenv)
    ENDIF
#else

    a => matrix%p(source,source)%block
    b => matrixb%p(source,source)%block

    CALL dsygst(itype,'U',n,a(1,1),n,b(1,1),n,info)

#endif

    if (info.ne.0) then
      CALL stop_program(" blacs_cholesky_reduce","Error dsygst",globenv)
    endif


    CALL timestop(0.0_wp,handle)

  END  SUBROUTINE blacs_cholesky_reduce

!******************************************************************
!
! op can be SOLVE (out = U^-1 * in ) or MULTIPLY   (out = U * in )
!
!******************************************************************
  SUBROUTINE blacs_cholesky_restore(matrix,neig,matrixb,matrixout,op,globenv)
    TYPE(blacs_matrix_type), POINTER          :: matrix,matrixb,matrixout
    INTEGER, INTENT(IN)                       :: neig
    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER ( LEN = * ), INTENT ( IN )      :: op

    REAL(wp), DIMENSION(:,:), POINTER         :: a,b,out
    integer                                   :: context,info,itype,handle
    INTEGER                                   :: output_unit,source,n
    LOGICAL                                   :: ionode
    REAL(wp)                                  :: alpha
    integer                                   :: nprow,npcol,myprow,mypcol,i
    INTEGER, DIMENSION(9)                     :: desca,descb,descout

    CALL timeset("blacs_cholesky_restore","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr
    source = globenv%source
    context = globenv%context
    n = matrix%nrow_global
    itype = 1
    IF (op .ne. "SOLVE" .and. op .ne. "MULTIPLY") THEN
       call stop_program("blacs_cholesky_restore","wrong argument op")
    ENDIF

#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    a => matrix%p(myprow,mypcol)%block
    desca(:) = matrix%descriptor(:)
    b => matrixb%p(myprow,mypcol)%block
    descb(:) = matrixb%descriptor(:)
    out => matrixout%p(myprow,mypcol)%block
    descout(:) = matrixout%descriptor(:)
    alpha=1.0_wp
    do i=1,neig
      CALL pdcopy(n,a(1,1),1,i,desca,1,out(1,1),1,i,descout,1)
    enddo
    IF (op.eq."SOLVE") THEN
      CALL pdtrsm('L','U','N','N',n,neig,alpha,b(1,1),1,1,descb,out(1,1),1,1,descout)
    ELSE
      CALL pdtrmm('L','U','N','N',n,neig,alpha,b(1,1),1,1,descb,out(1,1),1,1,descout)
    ENDIF
#else

    ! notice b is the cholesky guy
    a => matrix%p(source,source)%block
    b => matrixb%p(source,source)%block
    out => matrixout%p(source,source)%block
    alpha=1.0_wp
    CALL dcopy(neig*n,a(1,1),1,out(1,1),1)
    IF (op.eq."SOLVE") THEN
      CALL dtrsm('L','U','N','N',n,neig,alpha,b(1,1),n,out(1,1),n)
    ELSE
      CALL dtrmm('L','U','N','N',n,neig,alpha,b(1,1),n,out(1,1),n)
    ENDIF

#endif

    CALL timestop(0.0_wp,handle)

  END  SUBROUTINE blacs_cholesky_restore
! ****
  SUBROUTINE blacs_column_copy(msource,mtarget,ncol,globenv)
    TYPE(blacs_matrix_type), POINTER          :: msource,mtarget
    INTEGER, INTENT(IN)                       :: ncol
    TYPE(global_environment_type), INTENT(IN) :: globenv

    REAL(wp), DIMENSION(:,:), POINTER         :: a,b
    integer                                   :: context,nprow,npcol,&
                                                 myprow,mypcol,i,source,n
    INTEGER, DIMENSION(9)                     :: desca,descb

    source = globenv%source
    context = globenv%context
    n = msource%nrow_global

#if defined(__parallel)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    a => msource%p(myprow,mypcol)%block
    desca(:) = msource%descriptor(:)
    b => mtarget%p(myprow,mypcol)%block
    descb(:) = mtarget%descriptor(:)
    do i=1,ncol
      CALL pdcopy(n,a(1,1),1,i,desca,1,b(1,1),1,i,descb,1)
    enddo
#else

    a => msource%p(source,source)%block
    b => mtarget%p(source,source)%block
    CALL dcopy(ncol*n,a(1,1),1,b(1,1),1)
#endif
  END SUBROUTINE blacs_column_copy

!*************************************************************************
!computes all eigenvalues and vectors of a real symmetric matrix
!should be quite a bit faster than syevx for that case
!especially in parallel with thightly clustered evals
!needs more workspace
!*************************************************************************
  SUBROUTINE blacs_syevd(matrix,eigenvectors,eigenvalues,globenv)

    TYPE(blacs_matrix_type), POINTER          :: eigenvectors,matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), DIMENSION(:), INTENT(OUT)       :: eigenvalues

    INTEGER  :: context,handle,info,istat,liwork,lwork,mypcol,myprow,n,&
                npcol,nprow,source

    INTEGER, DIMENSION(9) :: descm,descv
    REAL(wp), DIMENSION(:,:), POINTER :: m,v
    REAL(wp), DIMENSION(:), POINTER :: work
    INTEGER, DIMENSION(:), POINTER  :: iwork

    CALL timeset("blacs_syevd","I","",handle)

    source = globenv%source
    context = globenv%context
    n = matrix%nrow_global

#if defined(__parallel)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    m => matrix%p(myprow,mypcol)%block
    descm(:) = matrix%descriptor(:)
    v => eigenvectors%p(myprow,mypcol)%block
    descv(:) = eigenvectors%descriptor(:)

    liwork=7*n+8*npcol+2
    allocate(iwork(liwork),STAT=istat)
    if (istat.ne.0) CALL stop_memory("blacs_syevd","iwork")
! work space query
    lwork = -1
    allocate(work(1),STAT=istat)
    if (istat.ne.0) CALL stop_memory("blacs_syevd","work")
    CALL PDSYEVD('V','U',n,m(1,1),1,1,descm,eigenvalues(1),v(1,1),1,1,descv, &
                  work(1),lwork,iwork(1),liwork,info)
    lwork = NINT(work(1))
    deallocate(work)
    allocate(work(lwork),STAT=istat)
    if (istat.ne.0) CALL stop_memory("blacs_syevd","lwork")
    CALL PDSYEVD('V','U',n,m(1,1),1,1,descm,eigenvalues(1),v(1,1),1,1,descv, &
                  work(1),lwork,iwork(1),liwork,info)
#else
    m => matrix%p(source,source)%block
    lwork=1+6*n+2*n**2
    liwork=5*n+3
    allocate(work(lwork),STAT=istat)
    if (istat.ne.0) CALL stop_memory("blacs_syevd","work")
    allocate(iwork(liwork),STAT=istat)
    if (istat.ne.0) CALL stop_memory("blacs_syevd","iwork")

    CALL DSYEVD('V','U', n, m(1,1), n, eigenvalues(1), work(1), lwork, iwork(1), liwork, info)
    CALL copy_blacs_to_blacs_matrix(matrix,eigenvectors)
#endif

    IF (info.ne.0) CALL stop_program("blacs_dsyevd","problems")

    deallocate(iwork)
    deallocate(work)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_syevd

! *****************************************************************************

  SUBROUTINE blacs_syevx(matrix,eigenvectors,eigenvalues,neig,work_syevx,&
                         globenv)
!  alternative , needs dsygvx
!  SUBROUTINE blacs_syevx(matrix,eigenvectors,eigenvalues,neig,work_syevx,&
!                         globenv,bmatrix)

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
    INTEGER, INTENT(IN)                       :: neig
    REAL(wp), DIMENSION(:), INTENT(OUT)       :: eigenvalues

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE blas_syevx (MODULE qs_blacs)"
    REAL(wp), PARAMETER :: abstol = 0.0_wp,&
                           orfac = -1.0_wp,&
                           vl = 0.0_wp,&
                           vu = 0.0_wp

!   *** Local variables ***

    INTEGER  :: context,handle,info,istat,liwork,lwork,m,mypcol,myprow,n,nb,&
                nn,np0,npcol,npe,nprow,nq0,nz,output_unit,source,itype
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
    context = globenv%context

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
    ! write(6,*) "done allocate"

!   *** Diagonalise matrix ***
!    if (PRESENT(bmatrix)) then
!       itype=1
!       b => bmatrix%p(myprow,mypcol)%block
!       descb(:) = matrix%descriptor(:)
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
       CALL dsyevx("V","I","U",n,a(1,1),n,vl,vu,1,neig,abstol,m,w,z(1,1),n,work(1),lwork,&
                iwork(1),ifail(1),info)
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
                        globenv)

!   Purpose: BLACS interface to the BLAS routine dsymm.

!   History: - Creation (07.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix_a,matrix_b,matrix_c
    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=1), INTENT(IN)              :: side,uplo
    REAL(wp), INTENT(IN)                      :: alpha,beta
    INTEGER, INTENT(IN)                       :: m,n

!   *** Local variables ***

    INTEGER :: context,handle,lda,ldb,ldc,mypcol,myprow,npcol,nprow,source

    INTEGER, DIMENSION(9) :: desca,descb,descc

    REAL(wp), DIMENSION(:,:), POINTER :: a,b,c

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_symm","I","",handle)

    source = globenv%source
    context = globenv%context
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    a => matrix_a%p(myprow,mypcol)%block
    desca(:) = matrix_a%descriptor(:)
    b => matrix_b%p(myprow,mypcol)%block
    descb(:) = matrix_b%descriptor(:)
    c => matrix_c%p(myprow,mypcol)%block
    descc(:) = matrix_c%descriptor(:)

    CALL pdsymm(side,uplo,m,n,alpha,a(1,1),1,1,desca,b(1,1),1,1,descb,beta,c(1,1),1,1,descc)

#else

    a => matrix_a%p(source,source)%block
    b => matrix_b%p(source,source)%block
    c => matrix_c%p(source,source)%block

    lda = matrix_a%nrow_global
    ldb = matrix_b%nrow_global
    ldc = matrix_c%nrow_global

    CALL dsymm(side,uplo,m,n,alpha,a(1,1),lda,b(1,1),ldb,beta,c(1,1),ldc)

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_symm

! *****************************************************************************

  SUBROUTINE blacs_syrk(uplo,trans,k,alpha,matrix_a,ia,ja,beta,matrix_c,&
                        globenv)

!   Purpose: BLACS interface to the BLAS routine dsyrk.

!   History: - Creation (07.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix_a,matrix_c
    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=1), INTENT(IN)              :: trans,uplo
    REAL(wp), INTENT(IN)                      :: alpha,beta
    INTEGER, INTENT(IN)                       :: ia,ja,k

!   *** Local variables ***

    INTEGER :: context,handle,lda,ldc,mypcol,myprow,n,npcol,nprow,source

    INTEGER, DIMENSION(9) :: desca,descc

    REAL(wp), DIMENSION(:,:), POINTER :: a,c

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_syrk","I","",handle)

    source = globenv%source
    context = globenv%context
    n = matrix_a%nrow_global
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    a => matrix_a%p(myprow,mypcol)%block
    desca(:) = matrix_a%descriptor(:)
    c => matrix_c%p(myprow,mypcol)%block
    descc(:) = matrix_c%descriptor(:)

    CALL pdsyrk(uplo,trans,n,k,alpha,a(1,1),ia,ja,desca,beta,c(1,1),1,1,descc)

#else

    a => matrix_a%p(source,source)%block
    c => matrix_c%p(source,source)%block

    lda = matrix_a%nrow_global
    ldc = matrix_c%nrow_global

    CALL dsyrk(uplo,trans,n,k,alpha,a(ia,ja),lda,beta,c(1,1),ldc)

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_syrk

! *****************************************************************************

  SUBROUTINE blacs_schur_product(matrix_a,matrix_b,matrix_c,globenv)

!   Purpose: Calculate the schur product of two matrices

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix_a,matrix_b,matrix_c
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: context,group,handle,icol_local,irow_local,mypcol,myprow,&
               ncol_local,npcol,nprow,nrow_local,source

    REAL(wp), DIMENSION(:,:), POINTER :: a,b,c

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_schur_product","I","",handle)

    group = globenv%group
    source = globenv%source
    context = globenv%context
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

#else

    myprow = source
    mypcol = source

#endif
    a => matrix_a%p(myprow,mypcol)%block
    b => matrix_b%p(myprow,mypcol)%block
    c => matrix_c%p(myprow,mypcol)%block

    nrow_local = matrix_a%p(myprow,mypcol)%nrow_local
    ncol_local = matrix_a%p(myprow,mypcol)%ncol_local

    DO icol_local=1,ncol_local
      DO irow_local=1,nrow_local
        c(irow_local,icol_local) = a(irow_local,icol_local)*b(irow_local,icol_local)
      END DO
    END DO

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_schur_product

! *****************************************************************************

  SUBROUTINE blacs_trace(matrix_a,matrix_b,trace,globenv)

!   Purpose: Calculate the trace of the product of two BLACS matrices.

!   History: - Creation (11.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix_a,matrix_b
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(OUT)                     :: trace

!   *** Local variables ***

    INTEGER :: context,group,handle,icol_local,irow_local,mypcol,myprow,&
               ncol_local,npcol,nprow,nrow_local,source

    REAL(wp), DIMENSION(:,:), POINTER :: a,b

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_trace","I","",handle)

    group = globenv%group
    source = globenv%source
    context = globenv%context
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

    INTEGER :: ipcol,iprow,istat,npcol,nprow,ncol_local,nrow_local,msize1,msize2
    REAL(wp), DIMENSION(:,:), POINTER :: a,b

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
          msize1=SIZE(source%p(iprow,ipcol)%block,1)*SIZE(source%p(iprow,ipcol)%block,2)
          msize2=SIZE(target%p(iprow,ipcol)%block,1)*SIZE(target%p(iprow,ipcol)%block,2)
          if (msize1.ne.msize2) then
             CALL stop_program("copy_blacs_to_blacs_matrix","no identical sizes")
          endif
          a => source%p(iprow,ipcol)%block
          b => target%p(iprow,ipcol)%block
          CALL DCOPY(msize1,a(1,1),1,b(1,1),1)
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

  SUBROUTINE copy_blacs_to_full_matrix(blacs_matrix,full_matrix,globenv)

!   Purpose: Copy a BLACS matrix to a full matrix.

!   History: - Creation (18.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: blacs_matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), DIMENSION(:,:), POINTER         :: full_matrix

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_blacs_to_full_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    INTEGER :: context,handle,icol_global,icol_local,ipcol,ipe,iprow,&
               irow_global,irow_local,istat,mypcol,mype,myprow,ncol_block,&
               ncol_global,ncol_local,npcol,npe,nprow,nrow_block,nrow_global,&
               nrow_local,source
    LOGICAL :: ionode

    REAL(wp), DIMENSION(:,:), POINTER :: blacs_block

#if defined(__parallel)
    INTEGER, EXTERNAL :: blacs_pnum,indxl2g

#endif
!   ---------------------------------------------------------------------------

    CALL timeset("copy_blacs_to_full_matrix","I","",handle)

    ionode = globenv%ionode
    source = globenv%source
    context = globenv%context

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
                         blacs_matrix%p(iprow,ipcol)%block(1,1),nrow_local,&
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
                         blacs_matrix%p(iprow,ipcol)%block(1,1),nrow_local,&
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

  SUBROUTINE copy_blacs_to_sparse_matrix(blacs_matrix,sparse_matrix,globenv)

!   Purpose: Copy a BLACS matrix to a sparse matrix. The BLACS matrix blocks
!            are deallocated during the copy procedure.

!   History: - Creation (06.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: blacs_matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(real_matrix_type), POINTER           :: sparse_matrix

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_blacs_to_sparse_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: context,group,handle,iblock_col,iblock_row,icol,icol_global,&
               icol_local,ipcol,ipe,iprow,irow,irow_global,irow_local,istat,&
               jpcol,jprow,mypcol,mype,myprow,nblock_row,ncol_block,&
               ncol_local,npcol,npe,nprow,nrow_block,nrow_local,source

    INTEGER, DIMENSION(:), POINTER    :: first_col,first_row,last_col,last_row
    REAL(wp), DIMENSION(:,:), POINTER :: blacs_block,sparse_block

#if defined(__parallel)
    INTEGER, EXTERNAL :: blacs_pnum,indxg2l,indxg2p

#endif
!   ---------------------------------------------------------------------------

    CALL timeset("copy_blacs_to_sparse_matrix","I","",handle)

    group = globenv%group
    source = globenv%source
    context = globenv%context

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

  SUBROUTINE copy_sparse_to_blacs_matrix(sparse_matrix,blacs_matrix,globenv)

!   Purpose: Copy a sparse matrix to a BLACS matrix. The BLACS matrix blocks
!            are allocated during the copy procedure.

!   History: - Creation (05.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: blacs_matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(real_matrix_type), POINTER           :: sparse_matrix

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_sparse_to_blacs_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: context,group,handle,iblock_col,iblock_row,icol,icol_global,&
               icol_local,ipcol,ipe,iprow,irow,irow_global,irow_local,istat,&
               jpcol,jprow,mypcol,mype,myprow,nblock_row,ncol_block,&
               ncol_local,npcol,npe,nprow,nrow_block,nrow_local,source

    INTEGER, DIMENSION(:), POINTER    :: first_col,first_row,last_col,last_row
    REAL(wp), DIMENSION(:,:), POINTER :: blacs_block,sparse_block

#if defined(__parallel)
    INTEGER, EXTERNAL :: blacs_pnum,indxg2l,indxg2p

#endif
!   ---------------------------------------------------------------------------

    CALL timeset("copy_sparse_to_blacs_matrix","I","",handle)

    group = globenv%group
    source = globenv%source
    context = globenv%context

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

!!****f* qs_blacs/deallocate_blacs_matrix_vect [1.0] *
!!
!!   NAME
!!     deallocate_blacs_matrix_vect
!!
!!   FUNCTION
!!     deallocate an array of pointers to blacs matrixes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - blacs_matrixes: the array of matrixes to deallocate
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
subroutine deallocate_blacs_matrix_vect(blacs_matrixes, error)
  type(blacs_matrix_p_type), dimension(:), pointer :: blacs_matrixes
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='blacs_matrixes',&
        routineP=moduleN//':'//routineN
  integer :: stat, i
  failure=.false.
  
  if (associated(blacs_matrixes)) then
     do i=1,size(blacs_matrixes)
        call deallocate_blacs_matrix(blacs_matrixes(i)%blacs_matrix)
     end do
     deallocate(blacs_matrixes,stat=stat)
     CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  end if
end subroutine deallocate_blacs_matrix_vect
!***************************************************************************

  SUBROUTINE finish_blacs(globenv)

!   Purpose: Release the resources of a BLACS context.

!   History: - Creation (22.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(INOUT) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE finish_blacs (MODULE qs_blacs)"

!   *** Local variables ***

    INTEGER :: context,group,ipe,istat,mype,npe,output_unit
    LOGICAL :: ionode

    INTEGER, DIMENSION(:), ALLOCATABLE :: pcon

!   ---------------------------------------------------------------------------

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr
    context = globenv%context
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

    globenv%context = 0
    globenv%nprow = 0
    globenv%npcol = 0

#endif
  END SUBROUTINE finish_blacs

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

    INTEGER :: context,mypcol,mype,myprow,npcol,npe,nprow,source

#if defined(__parallel)
    INTEGER, EXTERNAL :: blacs_pnum

#endif
!   ---------------------------------------------------------------------------

    source = globenv%source
    context = globenv%context
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
                                   nrow_block,ncol_block,nrow_local,ncol_local,&
                                   row_indices,col_indices,globenv, mpi_env)

!   Purpose: Return informations about the specified BLACS matrix.
!   always allocates the indices vector if present.
!   indices is an array that maps the local indices to the global ones

!   History: - Creation (08.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER        :: matrix
    CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: name
    INTEGER, OPTIONAL, INTENT(OUT)          :: ncol_block,ncol_global,&
                                               nrow_block,nrow_global,&
                                               nrow_local,ncol_local
    INTEGER, OPTIONAL, DIMENSION(:), POINTER:: row_indices,col_indices
    TYPE(global_environment_type), INTENT(IN), OPTIONAL :: globenv
    type(cp_mpi_env_type), intent(in), optional :: mpi_env

    INTEGER source,i,nprow,npcol,myprow,mypcol,context
#if defined(__parallel)
    INTEGER indxl2g
#endif

!   ---------------------------------------------------------------------------

    IF (PRESENT(name)) name = matrix%name
    IF (PRESENT(nrow_global)) nrow_global = matrix%nrow_global
    IF (PRESENT(ncol_global)) ncol_global = matrix%ncol_global
    IF (PRESENT(nrow_block)) nrow_block = matrix%nrow_block
    IF (PRESENT(ncol_block)) ncol_block = matrix%ncol_block

    IF (present(globenv)) THEN
        source=globenv%source
        context=globenv%context
    ENDIF
    if (present(mpi_env)) then
       source=mpi_env%source
       context=matrix%context
    end if

    IF (present(nrow_local)) THEN
      IF (.not. present(globenv).and..not.present(mpi_env)) then
         CALL stop_program("get_blacs_matrix_info","globenv")
      END IF
#if defined(__parallel)
          CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
          nrow_local = matrix%p(myprow,mypcol)%nrow_local
#else
          nrow_local = matrix%p(source,source)%nrow_local
#endif
    ENDIF

    IF (present(ncol_local)) THEN
      IF (.not. present(globenv).and..not.present(mpi_env)) then
         CALL stop_program("get_blacs_matrix_info","globenv")
      END IF
#if defined(__parallel)
          CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
          ncol_local = matrix%p(myprow,mypcol)%ncol_local
#else
          ncol_local = matrix%p(source,source)%ncol_local
#endif
    ENDIF

    IF (present(row_indices)) THEN
      IF (.not. present(globenv).and..not.present(mpi_env)) then
         CALL stop_program("get_blacs_matrix_info","globenv")
      END IF
#if defined(__parallel)
       CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
       allocate(row_indices(matrix%p(myprow,mypcol)%nrow_local))
       do i=1,matrix%p(myprow,mypcol)%nrow_local
         row_indices(i)=indxl2g(i,matrix%nrow_block,myprow,source,nprow)
       enddo
#else
       allocate(row_indices(matrix%p(source,source)%nrow_local))
       do i=1,matrix%p(source,source)%nrow_local
         row_indices(i)=i
       enddo
#endif
    ENDIF

    IF (present(col_indices)) THEN
      IF (.not. present(globenv).and..not.present(mpi_env)) then
         CALL stop_program("get_blacs_matrix_info","globenv")
      END IF
#if defined(__parallel)
       CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
       allocate(col_indices(matrix%p(myprow,mypcol)%ncol_local))
       do i=1,matrix%p(myprow,mypcol)%ncol_local
         col_indices(i)=indxl2g(i,matrix%ncol_block,mypcol,source,npcol)
       enddo
#else
       allocate(col_indices(matrix%p(source,source)%ncol_local))
       do i=1,matrix%p(source,source)%ncol_local
         col_indices(i)=i
       enddo
#endif
    ENDIF

  END SUBROUTINE get_blacs_matrix_info

! *****************************************************************************

  SUBROUTINE power_blacs_matrix(matrix,work,exponent,threshold,n_dependent,&
                                work_syevx,globenv)

!   Purpose: Raise the real symmetric n by n matrix to the power given by
!            exponent. All eigenvectors with a corresponding eigenvalue lower
!            than threshold are quenched.

!   History: - Creation (29.03.1999, Matthias Krack)
!            - Parallelised using BLACS and ScaLAPACK (06.06.2001, MK)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix,work
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(wp), INTENT(IN)                      :: exponent,threshold,work_syevx
    INTEGER, INTENT(OUT)                      :: n_dependent

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE power_blacs_matrix (MODULE qs_blacs)"

!   *** Local variables ***

    REAL(wp) :: f,p
    INTEGER  :: context,handle,icol_global,icol_local,ipcol,iprow,irow_global,&
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
    context = globenv%context
    n_dependent = 0
    p = 0.5_wp*exponent

    nrow_global = matrix%nrow_global
    ncol_global = matrix%ncol_global

    ALLOCATE (eigenvalues(ncol_global),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"eigenvalues",ncol_global*wp_size)
    eigenvalues(:) = 0.0_wp

!   *** Compute the eigenvectors and eigenvalues ***

    CALL blacs_syevx(matrix,work,eigenvalues,ncol_global,work_syevx,globenv)

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
    CALL blacs_syrk("U","N",ncol_global,1.0_wp,work,1,1,0.0_wp,matrix,globenv)

    DEALLOCATE (eigenvalues,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"eigenvalues")

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE power_blacs_matrix

! *****************************************************************************

  SUBROUTINE read_blacs_matrix(matrix,lunit,globenv)

!   Purpose: Read a BLACS matrix from the logical unit number "lunit".

!   History: - Creation (19.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: lunit

!   *** Local variables ***

    INTEGER :: context,i,j,mypcol,myprow,ncol_local,npcol,nprow,nrow_local,&
               source

!   ---------------------------------------------------------------------------

    source =  globenv%source
    context =  globenv%context
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

  SUBROUTINE start_blacs(globenv)

!   Purpose: Initialize a BLACS process grid. The BLACS context is returned.

!   History: - Creation (22.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(INOUT) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE start_blacs (MODULE qs_blacs)"

!   *** Local variables ***

    INTEGER :: context,group,ipe,istat,mypcol,mype,myprow,npcol,npe,nprow,&
               output_unit,source
    LOGICAL :: ionode

    INTEGER, DIMENSION(:), ALLOCATABLE :: pcol,pcon,prow

!   ---------------------------------------------------------------------------

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr
    source = globenv%source
    context = 0
    npcol = 0
    nprow = 0
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
    globenv%context = context
    globenv%nprow = nprow
    globenv%npcol = npcol

  END SUBROUTINE start_blacs

! *****************************************************************************

  SUBROUTINE blacs_transpose(matrix,matrixt,globenv)
! matrixt=matrix^T ! assumes no symmetry
!***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix,matrixt
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: context,handle,icol_global,icol_local,ipcol,iprow,irow_global,&
               irow_local,mypcol,myprow,ncol_block,ncol_global,ncol_local,&
               npcol,nprow,nrow_block,nrow_global,nrow_local,source

    INTEGER, DIMENSION(9) :: desca,descc

    REAL(wp), DIMENSION(:,:), POINTER :: a,c
    INTEGER :: i,j

!   ---------------------------------------------------------------------------

    CALL timeset("blacs_transpose","I","",handle)

    source = globenv%source
    context = globenv%context
    nrow_global = matrix%nrow_global
    ncol_global = matrix%ncol_global
#if defined(__parallel)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    a => matrix%p(myprow,mypcol)%block
    desca(:) = matrix%descriptor(:)
    c => matrixt%p(myprow,mypcol)%block
    descc(:) = matrixt%descriptor(:)
    CALL pdtran(nrow_global,ncol_global,1.0_wp,a(1,1),1,1,desca,0.0_wp,c(1,1),1,1,descc)
#else
    a => matrix%p(source,source)%block
    c => matrixt%p(source,source)%block
    DO j=1,ncol_global
     DO i=1,nrow_global
        c(j,i)=a(i,j)
     ENDDO
    ENDDO
#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE blacs_transpose

! *****************************************************************************

  SUBROUTINE symmetrise_blacs_matrix(matrix,work,globenv)

!   Purpose: Symmetrise a symmetric BLACS matrix.

!   History: - Creation (12.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix,work
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: context,handle,icol_global,icol_local,ipcol,iprow,irow_global,&
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
    context = globenv%context
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

    CALL pdtran(nrow_global,ncol_global,1.0_wp,a(1,1),1,1,desca,1.0_wp,c(1,1),1,1,descc)

#else

    a => matrix%p(source,source)%block

    CALL symmetrize_matrix(a,"upper_to_lower")

#endif
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE symmetrise_blacs_matrix

! *****************************************************************************

  SUBROUTINE write_blacs_matrix(matrix,lunit,globenv,formatted)

!   Purpose: Write a BLACS matrix to the logical unit number "lunit".

!   History: - Creation (19.06.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: lunit
    LOGICAL, INTENT(IN), OPTIONAL             :: formatted

!   *** Local variables ***

    INTEGER :: context,i,j,mypcol,myprow,ncol_local,npcol,nprow,nrow_local,&
               source

!   ---------------------------------------------------------------------------

    source =  globenv%source
    context =  globenv%context
#if defined(__parallel)

    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

#else

    myprow = source
    mypcol = source

#endif
    nrow_local = matrix%p(myprow,mypcol)%nrow_local
    ncol_local = matrix%p(myprow,mypcol)%ncol_local

    IF (PRESENT(formatted)) THEN
       IF (formatted) THEN
         write(6,*) matrix%name
         write(6,*) "row ,col ,value"
         DO j=1,ncol_local
         DO i=1,nrow_local
         WRITE (UNIT=lunit,FMT='(I5,I5,F6.2)') i,j,matrix%p(myprow,mypcol)%block(i,j)
         ENDDO
         ENDDO
       ELSE
          WRITE (UNIT=lunit)&
            ((matrix%p(myprow,mypcol)%block(i,j),i=1,nrow_local),&
                                                 j=1,ncol_local)
       ENDIF
    ELSE
          WRITE (UNIT=lunit) ((matrix%p(myprow,mypcol)%block(i,j),i=1,nrow_local),&
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
         end_col_of_sm,nblock_col
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
       call get_matrix_info(target_ptr, nrow=nrows, ncol=ncols,&
            nblock_col=nblock_col)
       CPAssert(ncols==nrows,cp_warning_level,routineP,error,failure)
       ncols=min(nrows,ncols)
       if (present(d_struct)) my_d_struct => d_struct
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
    CPAssert(size(my_d_struct)==nblock_col,cp_failure_level,routineP,error,failure)
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
       do sm_block_nr=1,nblock_col
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
    nullify(block_c, block_b)
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

    TYPE(blacs_matrix_type), POINTER          :: matrix_a, matrix_v, matrix_work
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

    TYPE(blacs_matrix_type), POINTER          ::  vecbuffer,matrix_overlap_hh
    TYPE(blacs_matrix_type), POINTER          ::  matrix_overlap_vh,matrix_overlap_vv

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

!   ---------------------------------------------------------------------------
    CALL timeset("blacs_jacobi_davidson","I","",handle)

    context = globenv%context

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

    time_init_start=cputime()
    time_full_start=cputime()

    rone=1.0_wp
    rzero=0.0_wp

    source = globenv%source
    n = matrix_a%nrow_global ! matrix dimension
    ksteps=5
    nhigh=MAX(ksteps*((neig/3)/ksteps),ksteps)  ! the upper part of the spectrum
    kblock=nhigh/ksteps
    nlow=neig-nhigh          ! the lower part of the spectrum
    nvec=neig+nhigh          ! the number of vectors used is subspace diag
    if (nhigh.lt.1 .or. nlow.lt.1 ) call stop_program("jacobi_davidson","nhigh, nlow, too small")
    if (nhigh.ge.neig/2 ) call stop_program("jacobi_davidson","nhigh, too large")

    levwork=1+6*nvec+2*nvec**2
    liwork=5*nvec+3
    allocate(evwork(levwork))
    allocate(iwork(liwork))
    allocate(evalsl(nvec))
    allocate(evalslold(nvec))
    allocate(to_be_used(nhigh))
    allocate(to_be_used_old(nhigh))
    evalslold=0.0_wp

    CALL get_blacs_matrix_info(matrix=matrix_a,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block)
    ncol = 2*nvec
    ! some buffer space
    CALL allocate_blacs_matrix(new_matrix=vecbuffer,&
                               nrow_global=n,&
                               ncol_global=ncol,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="VEC BUFFER",&
                               globenv=globenv)
    CALL allocate_blacs_matrix(new_matrix=matrix_overlap_vh,&
                               nrow_global=neig,&
                               ncol_global=nhigh,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="overlap_vh",&
                               globenv=globenv)
    CALL allocate_blacs_matrix(new_matrix=matrix_overlap_hh,&
                               nrow_global=nhigh,&
                               ncol_global=nhigh,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="overlap_hh",&
                               globenv=globenv)
    CALL allocate_blacs_matrix(new_matrix=matrix_overlap_vv,&
                               nrow_global=nvec,&
                               ncol_global=nvec,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="overlap_vv",&
                               globenv=globenv)

! put the right pointers in place
#if defined(__parallel)
    call stop_program("blacs_jacobi_davidson","is not parallel")
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
    a => matrix_a%p(myprow,mypcol)%block
    desca(:) = matrix_a%descriptor(:)

    v => matrix_v%p(myprow,mypcol)%block
    descv(:) = matrix_v%descriptor(:)

    work => matrix_work%p(myprow,mypcol)%block
    descwork(:) = matrix_work%descriptor(:)

    vecbuf => vecbuffer%p(myprow,mypcol)%block
    descvecbuf(:) = vecbuffer%descriptor(:)

    overlap_vh => matrix_overlap_vh%p(myprow,mypcol)%block
    descvh(:) = matrix_overlap_vh%descriptor(:)

    overlap_hh => matrix_overlap_hh%p(myprow,mypcol)%block
    deschh(:) = matrix_overlap_hh%descriptor(:)

    overlap_vv => matrix_overlap_vv%p(myprow,mypcol)%block
    descvv(:) = matrix_overlap_vv%descriptor(:)
#else
    a => matrix_a%p(source,source)%block
    v => matrix_v%p(source,source)%block
    work => matrix_work%p(source,source)%block
    vecbuf => vecbuffer%p(source,source)%block
    overlap_vh => matrix_overlap_vh%p(source,source)%block
    overlap_hh => matrix_overlap_hh%p(source,source)%block
    overlap_vv => matrix_overlap_vv%p(source,source)%block
#endif

! ****************************************************************************
! make the input symmetrical
! ****************************************************************************

    call symmetrise_blacs_matrix(matrix_a,matrix_work,globenv)


! ****************************************************************************
! get the A*V for the occupied orbitals
! ****************************************************************************

#if defined(__parallel)

#else
    call dcopy(n*neig,v,1,vecbuf,1)
#endif

    time_av_start=cputime()
#if defined(__parallel)

#else
    count_av=count_av+neig
    CALL dgemm('N','N',n,neig,n,rone,a,n,vecbuf(1,1),n,rzero,vecbuf(1,nvec+1),n)
#endif
    time_av_end=cputime()
    time_av_total=time_av_total+(time_av_end-time_av_start)

! ****************************************************************************
! get the V^T*A*V for the occupied orbitals
! ****************************************************************************
    time_suba_start=cputime()
#if defined(__parallel)

#else
    CALL DGEMM('T','N',neig,neig,n,rone,vecbuf(1,1),n,vecbuf(1,nvec+1),n,rzero,overlap_vv,nvec)
#endif
    time_suba_end=cputime()
    time_suba_total=time_suba_total+(time_suba_end-time_suba_start)

! ****************************************************************************
! diag subblock
! ****************************************************************************
    time_diag_start=cputime()
#if defined(__parallel)

#else
    CALL DSYEVD('V','U',neig,overlap_vv,nvec,evalslold,evwork,levwork,iwork,liwork,info)
#endif
    if (info.ne.0) call stop_program("jacobi_davidson","dsyevd overlap_vv 1")
    time_diag_end=cputime()
    time_diag_total=time_diag_total+(time_diag_end-time_diag_start)

! ****************************************************************************
! rotate vectors
! ****************************************************************************
    time_rot_start=cputime()
#if defined(__parallel)

#else
    CALL DGEMM('N','N',n,neig,neig,rone,vecbuf(1,     1),n,overlap_vv,nvec,rzero,v,n)
    CALL DCOPY(neig*n,v,1,vecbuf(1,1),1)
#endif
    time_rot_end=cputime()
    time_rot_total=time_rot_total+(time_rot_end-time_rot_start)

    ! start with the upper part of the spectrum
    DO i=1,nhigh
        to_be_used(i)=neig-nhigh+i
    ENDDO
    to_be_used_old=to_be_used

! ****************************************************************************
! we LU decompose the matrix a-shift*I
! ****************************************************************************
    time_getrf_start=cputime()
    call copy_blacs_to_blacs_matrix(matrix_a,matrix_work)
    lipivot = MAX(n,128) ! somehow too long ...
    allocate (ipivot(lipivot),stat=istat)
    if (istat.ne.0) call stop_memory("blacs_jacobi_davidson","ipivot")
    ! find some gap in the lower half of the spectrum
    gap=0.0_wp
    igap=-1
    DO i=1,neig/2
       gapnow=evalslold(i+1)-evalslold(i)
       if (gapnow.ge.gap) then
          gap=gapnow
          igap=i
       endif
    ENDDO
    shift=0.5_wp*(evalslold(igap+1)+evalslold(igap))
    ! shift=1.1*evalslold(1)
#if defined(__parallel)
    CALL pdgetrf(n,n,work,1,1,descwork,ipivot,info)
#else
    ! put the singularity in the middle of the spectrum
    DO i=1,n
       work(i,i)=work(i,i)-shift
    ENDDO
    CALL dgetrf(n,n,work,n,ipivot,info)
#endif
    if (info.ne.0) call stop_program("blacs_jacobi_davidson","problems in getrf")
    time_getrf_end=cputime()
    time_getrf_total=time_getrf_total+(time_getrf_end-time_getrf_start)

    time_init_end=cputime()
    time_init_total=time_init_total+(time_init_end-time_init_start)

! ****************************************************************************
! start main iteration
! ****************************************************************************
    DO iter=1,30

       ! ************************************************************************
       !  extension of the basis
       ! ************************************************************************
       time_extend_start=cputime()
       ! copy the needed vectors in place
#if defined(__parallel)
#else
       DO i=1,nhigh
          CALL DCOPY(n,vecbuf(1,to_be_used(i)),1,vecbuf(1,neig+i),1)
       ENDDO
#endif
       ! *********************************************************************
       ! get (LU)^-1 * vectors
       ! *********************************************************************
       time_getrs_start=cputime()
#if defined(__parallel)

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
          if (info.ne.0) call stop_program("jacobi_davidson","Dgetrs")
       time_getrs_end=cputime()
       time_getrs_total=time_getrs_total+(time_getrs_end-time_getrs_start)


       ! *********************************************************************
       ! keep the new vectors orthogonal, is numerically rather hard, it seems
       ! *********************************************************************
#if defined(__parallel)

#else
       time_ortho_start=cputime()

       ! ortho to the existing basis vectors
       CALL DGEMM('T','N',neig,nhigh,n,rone,vecbuf(1,1),n,vecbuf(1,neig+1),n,rzero,overlap_vh,neig)
       CALL DGEMM('N','N',n,nhigh,neig,-rone,vecbuf(1,1),n,overlap_vh,neig,rone,vecbuf(1,neig+1),n)
       ! and ortho again within the new set
       CALL DSYRK('U','T',nhigh,n,rone,vecbuf(1,neig+1),n,rzero,overlap_hh,nhigh)
       CALL DPOTRF('U',nhigh,overlap_hh,nhigh,info)
       if (info.ne.0) call stop_program("jacobi_davidson","DPOTRF overlap_hh")
       CALL DTRSM('R','U','N','N',n,nhigh,rone,overlap_hh,nhigh,vecbuf(1,neig+1),n)
       ! ortho to the existing basis vectors
       CALL DGEMM('T','N',neig,nhigh,n,rone,vecbuf(1,1),n,vecbuf(1,neig+1),n,rzero,overlap_vh,neig)
       CALL DGEMM('N','N',n,nhigh,neig,-rone,vecbuf(1,1),n,overlap_vh,neig,rone,vecbuf(1,neig+1),n)
       ! and ortho again within the new set
       CALL DSYRK('U','T',nhigh,n,rone,vecbuf(1,neig+1),n,rzero,overlap_hh,nhigh)
       CALL DPOTRF('U',nhigh,overlap_hh,nhigh,info)
       if (info.ne.0) call stop_program("jacobi_davidson","DPOTRF overlap_hh")
       CALL DTRSM('R','U','N','N',n,nhigh,rone,overlap_hh,nhigh,vecbuf(1,neig+1),n)

       time_ortho_end=cputime()
       time_ortho_total=time_ortho_total+(time_ortho_end-time_ortho_start)
#endif

       ! *********************************************************************
       ! compute A*vectors for the new set
       ! *********************************************************************
       time_av_start=cputime()
#if defined(__parallel)

#else
       CALL dgemm('N','N',n,nhigh,n,rone,a,n,vecbuf(1,neig+1),n,rzero,vecbuf(1,nvec+neig+1),n)
       count_av=count_av+nhigh
#endif
       time_av_end=cputime()
       time_av_total=time_av_total+(time_av_end-time_av_start)

       time_extend_end=cputime()
       time_extend_total=time_extend_total+(time_extend_end-time_extend_start)

       ! *********************************************************************
       ! get v^T*(A*v)
       ! *********************************************************************
       time_suba_start=cputime()
#if defined(__parallel)

#else
       CALL DSCAL(nvec*nvec,rzero,overlap_vv,1)
       DO i=1,neig
          overlap_vv(i,i)=evalslold(i)
       ENDDO
       CALL DGEMM('T','N',nvec,nhigh,n,rone,vecbuf(1,1),n,vecbuf(1,nvec+neig+1),n, &
                                                    rzero,overlap_vv(1,neig+1),nvec)
       !CALL DGEMM('T','N',nvec,nvec,n,rone,vecbuf(1,1),n,vecbuf(1,nvec+1),n, &
       !                                                  rzero,overlap_vv,nvec)
#endif
       time_suba_end=cputime()
       time_suba_total=time_suba_total+(time_suba_end-time_suba_start)

       ! *********************************************************************
       ! diagonalise the matrix v^T*(A*v)
       ! *********************************************************************
       time_diag_start=cputime()
#if defined(__parallel)

#else
       CALL DSYEVD('V','U',nvec,overlap_vv,nvec,evalsl,evwork,levwork,iwork,liwork,info)
#endif
       if (info.ne.0) call stop_program("jacobi_davidson","dsyevd overlap_vv")
       time_diag_end=cputime()
       time_diag_total=time_diag_total+(time_diag_end-time_diag_start)

       ! *********************************************************************
       ! here we select the eigenvectors with the largest change
       ! they will be used in the next iteration
       ! *********************************************************************
       ! write(6,*) "evalsl",NINT(-log(abs(evalsl(1:neig)-evalslold(1:neig))+1E-20)/log(10.0_wp))
       evwork(1:neig)=abs(evalsl(1:neig)-evalslold(1:neig))
       call sort(evwork,neig,iwork)
       write(6,*) "changes",iwork(neig),evwork(neig)
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
             if (this_one.eq.to_be_used_old(j)) this_one=0
           ENDDO
           if (this_one.ne.0) exit
          ENDDO
          to_be_used(i)=this_one
       ENDDO
       to_be_used_old=to_be_used
       evalslold=evalsl

       ! *********************************************************************
       ! rotate V, and AV so that the lowest neig are first
       ! *********************************************************************
       time_rot_start=cputime()
#if defined(__parallel)

#else
       CALL DGEMM('N','N',n,neig,nvec,rone,vecbuf(1,     1),n,overlap_vv,nvec,rzero,v,n)
       CALL DCOPY(neig*n,v,1,vecbuf(1,1),1)
#endif
       time_rot_end=cputime()
       time_rot_total=time_rot_total+(time_rot_end-time_rot_start)

       if (evwork(neig).le.tol) EXIT

    ENDDO

    ! *********************************************************************
    ! make the final vectors exactly ortho
    ! *********************************************************************
    time_ortho_start=cputime()
#if defined(__parallel)

#else
      CALL DSYRK('U','T',neig,n,rone,v(1,1),n,rzero,overlap_vv,nvec)
      CALL DPOTRF('U',neig,overlap_vv,nvec,info)
      ! in this case the vectors are exactly degenerate
      if (info.ne.0) call stop_program("jacobi_davidson","DPOTRF final wave functions")
      CALL DTRSM('R','U','N','N',n,neig,rone,overlap_vv,nvec,v(1,1),n)
#endif
    time_ortho_end=cputime()
    time_ortho_total=time_ortho_total+(time_ortho_end-time_ortho_start)

    evals(1:neig)=evalsl(1:neig)

    CALL deallocate_blacs_matrix(vecbuffer)
    CALL deallocate_blacs_matrix(matrix_overlap_vh)
    CALL deallocate_blacs_matrix(matrix_overlap_hh)
    CALL deallocate_blacs_matrix(matrix_overlap_vv)
    deallocate(ipivot)
    deallocate(evwork)
    deallocate(iwork)
    deallocate(evalsl)
    deallocate(evalslold)
    deallocate(to_be_used)
    deallocate(to_be_used_old)

    time_full_end=cputime()
    time_full_total=time_full_total+(time_full_end-time_full_start)
    write(6,*) "time full    : ", time_full_total
    write(6,*) "time init    : ", time_init_total
    write(6,*) "time extend  : ", time_extend_total
    write(6,*) "---------------"
    write(6,*) "time getrf   : ", time_getrf_total
    write(6,*) "time getrs   : ", time_getrs_total
    write(6,*) "time av      : ", time_av_total
    write(6,*) "time diag    : ", time_diag_total
    write(6,*) "time rot     : ", time_rot_total
    write(6,*) "time suba    : ", time_suba_total
    write(6,*) "time ortho   : ", time_ortho_total
    write(6,*) "---------------"
    write(6,*) "# neig/nhigh/n ", neig, "/",nhigh,"/",n
    write(6,*) "# A*V /t       ", count_av,"/",time_av_total/count_av
    write(6,*) "# (L*U)^-1*V/t ", count_getrs,"/",time_getrs_total/count_getrs

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
subroutine blacs_scale_and_d(blacs_matrix,mpi_env,scale, diagonal_add,error)
    type(blacs_matrix_type), pointer :: blacs_matrix
    type(cp_mpi_env_type), intent(in) :: mpi_env
    real(kind=wp), intent(in), optional :: scale, diagonal_add
    type(cp_error_type), intent(inout), optional :: error
    
    integer :: i,j,nrow_local, ncol_local, nprow,npcol,myprow,mypcol
    integer, dimension(:), pointer :: row_indices, col_indices
    real(kind=wp) :: sc, diag
    real(kind=wp), dimension(:,:), pointer :: blacs_block

    sc=1.0_wp;  diag=0.0_wp
    if (present(scale)) sc=scale
    if (present(diagonal_add)) diag=diagonal_add
    
#if defined(__parallel)
    call blacs_gridinfo(blacs_matrix%context,nprow,npcol,myprow,mypcol)
    blacs_block => blacs_matrix%p(myprow,mypcol)%block
#else
    blacs_block => blacs_matrix%p(mpi_env%source,mpi_env%source)%block
#endif
    if (diag==0.0_wp) then
      blacs_block=sc*blacs_block
    else
      call get_blacs_matrix_info(matrix=blacs_matrix, &
          nrow_local =nrow_local ,ncol_local =ncol_local , &
          row_indices=row_indices,col_indices=col_indices,&
          mpi_env=mpi_env)
      if (sc==0.0) then
        call dcopy(size(blacs_block,1)*size(blacs_block,2),0.0_wp,0,blacs_block,1)
        do j=1,ncol_local
          do i=1,nrow_local
            if (col_indices(j)== row_indices(i)) then
              blacs_block(i,j)=diag
            end if
          end do
        end do
      else
        do j=1,ncol_local
          do i=1,nrow_local
            if (col_indices(j)== row_indices(i)) then
              blacs_block(i,j)=diag+sc*blacs_block(i,j)
            else
              blacs_block(i,j)=sc*blacs_block(i,j)
            end if
          end do
        end do
      end if
      
      deallocate(row_indices,col_indices)
    end if
end subroutine blacs_scale_and_d

END MODULE qs_blacs
