!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/sparse_matrix_types [1.0] *
!!
!!   NAME
!!     sparse_matrix_types
!!
!!   FUNCTION
!!     Define the matrix data types.
!!
!!   AUTHOR
!!     MK (23.06.2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE sparse_matrix_types

  USE kinds, ONLY: int_size,&
                   wp => dp,&
                   wp_size => dp_size

  USE termination, ONLY: stop_memory,&
                         stop_program
  USE timings,     ONLY: timeset,&
                         timestop
  USE global_types,        ONLY: global_environment_type
  USE message_passing,     ONLY: mp_sum

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: moduleN='sparse_matrix_types'

! *** Define a real matrix data type ***

  TYPE real_block_node_type
    PRIVATE
    TYPE(real_block_node_type), POINTER :: next_block_node
    INTEGER                             :: block_col
    REAL(wp), DIMENSION(:,:), POINTER   :: block
  END TYPE real_block_node_type

  TYPE real_block_list_type
    PRIVATE
    TYPE(real_block_node_type), POINTER :: first_block_node
  END TYPE real_block_list_type

  TYPE real_matrix_type
    PRIVATE
    TYPE(real_block_list_type), DIMENSION(:), POINTER :: block_list
    CHARACTER(LEN=60)                                 :: name
    CHARACTER(LEN=40)                                 :: symmetry
    INTEGER                                           :: nblock_col,&
                                                         nblock_row,&
                                                         ncol,nrow
    INTEGER, DIMENSION(:), POINTER                    :: first_col,&
                                                         first_row,&
                                                         last_col,&
                                                         last_row
  END TYPE real_matrix_type

  TYPE real_matrix_p_type
    TYPE(real_matrix_type), POINTER :: matrix
  END TYPE real_matrix_p_type

! *** Public data types ***

  PUBLIC :: real_block_node_type,&
            real_matrix_p_type,&
            real_matrix_type

! *** Public subroutines ***

  PUBLIC :: add_block_node,&
            add_matrices,&
            allocate_matrix,&
            allocate_matrix_set,&
            copy_matrix,&
            copy_sparse_to_full_matrix,&
            deallocate_matrix,&
            deallocate_matrix_row,&
            deallocate_matrix_set,&
            get_block_node,&
            get_matrix_info,&
            merge_matrices,&
            put_block_node,&
            replicate_matrix,&
            replicate_matrix_structure,&
            set_matrix,&
            symmetrise_diagonal_blocks,&
            transfer_matrix, &
            sparse_times_replicated, &
            sparse_plus_rep_rept, &
            remove_block_node, &
            write_blocks_maxabsval, &
            add_all_real_matrix_blocks

! *** Public functions ***

  PUBLIC :: first_block_node,&
            next_block_node,&
            find_block_node

! *****************************************************************************

  INTERFACE add_block_node
    MODULE PROCEDURE add_real_matrix_block
  END INTERFACE

  INTERFACE add_matrices
    MODULE PROCEDURE add_real_matrices
  END INTERFACE

  INTERFACE allocate_matrix
    MODULE PROCEDURE allocate_real_matrix,allocate_full_real_matrix
  END INTERFACE

  INTERFACE allocate_matrix_set
    MODULE PROCEDURE allocate_real_matrix_set
  END INTERFACE

  INTERFACE copy_matrix
    MODULE PROCEDURE copy_real_matrix
  END INTERFACE

  INTERFACE deallocate_matrix
    MODULE PROCEDURE deallocate_real_matrix
  END INTERFACE

  INTERFACE deallocate_matrix_row
    MODULE PROCEDURE deallocate_real_matrix_row
  END INTERFACE

  INTERFACE deallocate_matrix_set
    MODULE PROCEDURE deallocate_real_matrix_set
  END INTERFACE

  INTERFACE find_block_node
    MODULE PROCEDURE find_real_block_node
  END INTERFACE

  INTERFACE first_block_node
    MODULE PROCEDURE first_real_block_node
  END INTERFACE

  INTERFACE get_block_node
    MODULE PROCEDURE get_real_block_node,get_real_matrix_block
  END INTERFACE

  INTERFACE merge_matrices
    MODULE PROCEDURE merge_real_matrices
  END INTERFACE

  INTERFACE next_block_node
    MODULE PROCEDURE next_real_block_node
  END INTERFACE

  INTERFACE put_block_node
    MODULE PROCEDURE put_real_block_node,put_real_matrix_block
  END INTERFACE

  INTERFACE replicate_matrix
    MODULE PROCEDURE replicate_real_matrix
  END INTERFACE

  INTERFACE replicate_matrix_structure
    MODULE PROCEDURE replicate_real_matrix_structure
  END INTERFACE

  INTERFACE set_matrix
    MODULE PROCEDURE set_real_matrix
  END INTERFACE

  INTERFACE transfer_matrix
    MODULE PROCEDURE transfer_real_matrix
  END INTERFACE

  INTERFACE remove_block_node
     MODULE PROCEDURE remove_real_matrix_block
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE add_real_matrices(c,alpha,a,beta,b)

!   Purpose: Calculate the sum of the  sparse matrices a and b which have the
!            same structure.

!            c = alpha*a + beta*b

!   History: - Creation (10.04.2002, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER :: a,b,c
    REAL(wp), INTENT(IN)            :: alpha,beta

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
         "SUBROUTINE add_real_matrices (MODULE sparse_matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: iblock_col,iblock_row

    REAL(wp), DIMENSION(:,:), POINTER :: a_block,b_block,c_block

!   ---------------------------------------------------------------------------

!   *** Check the association status of the matrices ***

    IF (.NOT.ASSOCIATED(a)) THEN
       CALL stop_program(routine,"The matrix pointer a is not associated")
    END IF

    IF (.NOT.ASSOCIATED(b)) THEN
       CALL stop_program(routine,"The matrix pointer b is not associated")
    END IF

    IF (.NOT.ASSOCIATED(c)) THEN
       CALL stop_program(routine,"The matrix pointer c is not associated")
    END IF

!   *** Add matrix elements ***

    IF ((alpha == 1.0_wp).AND.(beta == 1.0_wp)) THEN ! make common case fast

       DO iblock_row=1,c%nblock_row

          block_node => first_block_node(c,iblock_row)

          DO WHILE (ASSOCIATED(block_node))

             CALL get_block_node(block_node=block_node,&
                                 block_col=iblock_col,&
                                 block=c_block)

             CALL get_block_node(matrix=a,&
                                 block_row=iblock_row,&
                                 block_col=iblock_col,&
                                 block=a_block)

             CALL get_block_node(matrix=b,&
                                 block_row=iblock_row,&
                                 block_col=iblock_col,&
                                 block=b_block)
             
             c_block(:,:) = a_block(:,:) + b_block(:,:)

             block_node => block_node%next_block_node

          END DO

       END DO

    ELSE

       DO iblock_row=1,c%nblock_row

          block_node => first_block_node(c,iblock_row)

          DO WHILE (ASSOCIATED(block_node))

            CALL get_block_node(block_node=block_node,&
                                block_col=iblock_col,&
                                block=c_block)

            CALL get_block_node(matrix=a,&
                                block_row=iblock_row,&
                                block_col=iblock_col,&
                                block=a_block)

            CALL get_block_node(matrix=b,&
                                block_row=iblock_row,&
                                block_col=iblock_col,&
                                block=b_block)

            c_block(:,:) = alpha*a_block(:,:) + beta*b_block(:,:)

            block_node => block_node%next_block_node

          END DO

       END DO

    END IF

  END SUBROUTINE add_real_matrices

! *****************************************************************************

  SUBROUTINE add_real_matrix_block(matrix,block_row,block_col,block)

!   Purpose: Allocate and initialize a new block node.

!   History: - Creation (28.07.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER             :: matrix
    INTEGER, INTENT(IN)                         :: block_col,block_row
    REAL(wp), DIMENSION(:,:), OPTIONAL, POINTER :: block

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE add_real_matrix_block (MODULE sparse_matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node, node_p

    INTEGER :: istat,ncol,nrow

!   ---------------------------------------------------------------------------

!   *** Calculate the block dimensions ***

    nrow = matrix%last_row(block_row) - matrix%first_row(block_row) + 1
    ncol = matrix%last_col(block_col) - matrix%first_col(block_col) + 1

    ALLOCATE (block_node,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"Block node of "//matrix%name,0)

!   *** Define the data set of the new block node ***

    block_node%block_col = block_col

    ALLOCATE (block_node%block(nrow,ncol),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"block_node%block",wp_size*nrow*ncol)
    END IF
    block_node%block(:,:) = 0.0_wp

    IF (PRESENT(block)) THEN
      IF (ASSOCIATED(block)) THEN
        IF ((SIZE(block,1) == nrow).AND.(SIZE(block,2) == ncol)) THEN
          block_node%block(:,:) = block(:,:)
        ELSE
          CALL stop_program(routine,"Incompatible block dimensions")
        END IF
      ELSE
        block => block_node%block
      END IF
    END IF

!   *** Link the new block node to the block list ***

!   *** Insert the new blocks in an ordered fashion (block_col are in ascending order)
    node_p => matrix%block_list(block_row)%first_block_node
    IF (.NOT.ASSOCIATED(node_p)) THEN
       NULLIFY(block_node%next_block_node)
       matrix%block_list(block_row)%first_block_node => block_node
    ELSE IF (node_p%block_col > block_node%block_col) THEN
       block_node%next_block_node => matrix%block_list(block_row)%first_block_node
       matrix%block_list(block_row)%first_block_node => block_node
    ELSE
       DO WHILE (ASSOCIATED(node_p%next_block_node))
          IF (node_p%next_block_node%block_col > block_node%block_col) EXIT
          node_p => node_p%next_block_node
       END DO
       block_node%next_block_node => node_p%next_block_node
       node_p%next_block_node => block_node
    END IF

  END SUBROUTINE add_real_matrix_block

!*****************************************************************************************
! writes for all local blocks the maximum absolute value of 
! the block entries to unit fileunit
!*****************************************************************************************

  SUBROUTINE write_blocks_maxabsval(sparse_matrix,fileunit)

    TYPE(real_matrix_type), POINTER   :: sparse_matrix
    INTEGER, INTENT(IN)               :: fileunit

!   *** Local parameters ***
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE write_block_magnitude (MODULE sparse_matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: iblock_col,iblock_row,icol,icol_block,irow,irow_block,istat,&
               ncol,nrow,sbncol,sbnrow,  handle

    REAL(wp), DIMENSION(:,:), POINTER :: sparse_block

    DO iblock_row=1,sparse_matrix%nblock_row

      block_node => first_block_node(sparse_matrix,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            block=sparse_block)

        write(fileunit,*) iblock_row,iblock_col,maxval(abs(sparse_block))

        block_node => next_block_node(block_node)

      END DO

    END DO

  END SUBROUTINE write_blocks_maxabsval

!******************************************************************************************
! add all real matrix blocks, but distributed over all cpus, keeps into account symmetric
!******************************************************************************************
  SUBROUTINE add_all_real_matrix_blocks(sparse_matrix,globenv)

    TYPE(real_matrix_type), POINTER   :: sparse_matrix
    TYPE(global_environment_type), INTENT(IN)      :: globenv

    integer iblock_row,iblock_col,ij
    ij=0
    DO iblock_row=1,sparse_matrix%nblock_row
      DO iblock_col=1,sparse_matrix%nblock_col
         if (sparse_matrix%symmetry.eq."symmetric") then
            if (iblock_col.lt.iblock_row) CYCLE
         endif
         ij=ij+1
         if (MOD(ij,globenv%num_pe).eq. globenv%mepos) then
            CALL add_real_matrix_block(sparse_matrix, &
                 block_row=iblock_row,block_col=iblock_col)
         endif
      END DO
    END DO

  END SUBROUTINE add_all_real_matrix_blocks

! *****************************************************************************
  SUBROUTINE allocate_real_matrix(matrix,nrow,ncol,nblock_row,nblock_col,&
                                  first_row,last_row,first_col,last_col,&
                                  matrix_name,matrix_symmetry)

!   Purpose: Allocate and initialize a real matrix at the real_matrix_type
!            level.

!   History: - Creation (16.06.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER :: matrix
    CHARACTER(LEN=*), INTENT(IN)    :: matrix_name,matrix_symmetry
    INTEGER, INTENT(IN)             :: nblock_row,nblock_col,ncol,nrow

    INTEGER, DIMENSION(nblock_row), INTENT(IN) :: first_row,last_row
    INTEGER, DIMENSION(nblock_col), INTENT(IN) :: first_col,last_col

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_real_matrix (MODULE sparse_matrix_types)"

!   *** Local variables ***

    INTEGER :: irow,istat

!   ---------------------------------------------------------------------------

!   *** Deallocate the old matrix ***

    IF (ASSOCIATED(matrix)) CALL deallocate_matrix(matrix)

!   *** Allocate a set of block lists ***

    ALLOCATE (matrix,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,matrix%name,0)

    ALLOCATE (matrix%block_list(nblock_row),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,matrix%name,0)

    matrix%name = matrix_name
    matrix%symmetry = matrix_symmetry

    matrix%nblock_row = nblock_row
    matrix%nblock_col = nblock_col

    matrix%nrow = nrow
    matrix%ncol = ncol

    ALLOCATE (matrix%first_row(nblock_row),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"matrix%first_row",int_size*nblock_row)
    END IF
    matrix%first_row(:) = first_row(:)

    ALLOCATE (matrix%last_row(nblock_row),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"matrix%last_row",int_size*nblock_row)
    END IF
    matrix%last_row(:) = last_row(:)

    ALLOCATE (matrix%first_col(nblock_col),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"matrix%first_col",int_size*nblock_col)
    END IF
    matrix%first_col(:) = first_col(:)

    ALLOCATE (matrix%last_col(nblock_col),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"matrix%last_col",int_size*nblock_col)
    END IF
    matrix%last_col(:) = last_col(:)

!   *** Initialize all block lists ***

    DO irow=1,nblock_row
      NULLIFY (matrix%block_list(irow)%first_block_node)
    END DO

  END SUBROUTINE allocate_real_matrix

! *****************************************************************************

  SUBROUTINE allocate_real_matrix_set(matrix_set,nmatrix)

!   Purpose: Allocate and initialize a real matrix set.

!   History: - Creation (13.03.2002, MK)

!   ***************************************************************************

    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: matrix_set
    INTEGER, INTENT(IN)                             :: nmatrix

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_real_matrix_set (MODULE sparse_matrix_types)"

!   *** Local variables ***

    INTEGER :: imatrix,istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(matrix_set)) CALL deallocate_matrix_set(matrix_set)

    ALLOCATE (matrix_set(nmatrix),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"matrix_set",0)

    DO imatrix=1,nmatrix
      NULLIFY (matrix_set(imatrix)%matrix)
    END DO

  END SUBROUTINE allocate_real_matrix_set

! *****************************************************************************

  SUBROUTINE allocate_full_real_matrix(matrix,nrow,ncol,matrix_name,&
                                       matrix_symmetry)

!   Purpose: Allocate and initialize a real matrix structure with only one
!            block.

!   History: - Creation (10.07.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER :: matrix
    CHARACTER(LEN=*), INTENT(IN)    :: matrix_name,matrix_symmetry
    INTEGER, INTENT(IN)             :: ncol,nrow

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_full_real_matrix (MODULE sparse_matrix_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

!   *** Deallocate the old matrix ***

    IF (ASSOCIATED(matrix)) CALL deallocate_matrix(matrix)

!   *** Allocate a set of block lists ***

    ALLOCATE (matrix,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,matrix%name,0)

    ALLOCATE (matrix%block_list(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,matrix%name,0)

    matrix%name = matrix_name
    matrix%symmetry = matrix_symmetry

    matrix%nblock_row = 1
    matrix%nblock_col = 1

    matrix%nrow = nrow
    matrix%ncol = ncol

    ALLOCATE (matrix%first_row(1),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"matrix%first_row",int_size)
    END IF
    matrix%first_row(1) = 1

    ALLOCATE (matrix%last_row(1),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"matrix%last_row",int_size)
    END IF
    matrix%last_row(1) = nrow

    ALLOCATE (matrix%first_col(1),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"matrix%first_col",int_size)
    END IF
    matrix%first_col(1) = 1

    ALLOCATE (matrix%last_col(1),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"matrix%last_col",int_size)
    END IF
    matrix%last_col(1) = ncol

    NULLIFY (matrix%block_list(1)%first_block_node)

    CALL add_real_matrix_block(matrix=matrix,&
                               block_row=1,&
                               block_col=1)

  END SUBROUTINE allocate_full_real_matrix

! *****************************************************************************

  SUBROUTINE copy_sparse_to_full_matrix(sparse_matrix,full_matrix)

!   Purpose: Copy the matrix blocks of a sparse matrix to the corresponding
!            full matrix which is allocated in this routine.

!   History: - Creation (19.06.2001, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER   :: sparse_matrix
    REAL(wp), DIMENSION(:,:), POINTER :: full_matrix

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_sparse_to_full_matrix (MODULE sparse_matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: iblock_col,iblock_row,icol,icol_block,irow,irow_block,istat,&
               ncol,nrow

    REAL(wp), DIMENSION(:,:), POINTER :: sparse_block

!   ---------------------------------------------------------------------------

!   *** Check the association status of the input matrix ***

    IF (.NOT.ASSOCIATED(sparse_matrix)) THEN
      CALL stop_program(routine,"The input matrix pointer is not associated")
    END IF

    IF (ASSOCIATED(full_matrix)) DEALLOCATE (full_matrix)

    nrow = sparse_matrix%nrow
    ncol = sparse_matrix%ncol

    ALLOCATE (full_matrix(nrow,ncol),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"full_matrix",nrow*ncol*wp_size)
    full_matrix(:,:) = 0.0_wp

!   *** Traverse all block nodes of the sparse matrix ***

    DO iblock_row=1,sparse_matrix%nblock_row

      block_node => first_block_node(sparse_matrix,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            block=sparse_block)

        icol_block = 1

        DO icol=sparse_matrix%first_col(iblock_col),&
                sparse_matrix%last_col(iblock_col)

          irow_block = 1

          DO irow=sparse_matrix%first_row(iblock_row),&
                  sparse_matrix%last_row(iblock_row)

            full_matrix(irow,icol) = sparse_block(irow_block,icol_block)

            irow_block = irow_block + 1

          END DO

          icol_block = icol_block + 1

        END DO

        block_node => next_block_node(block_node)

      END DO

    END DO

  END SUBROUTINE copy_sparse_to_full_matrix

! *****************************************************************************
! multiplies a sparse matrix times ncol replicated vectors
! (that are stored as rows ! in v_in leading dimension nblock )
! if the sparse_matrix is symmetric, off diagonal blocks will be used twice
! it is still assumed that the diagonal blocks are symmetric (and hence full)
! Joost VandeVondele july 2002
! *****************************************************************************
  SUBROUTINE sparse_times_replicated(sparse_matrix,rep_v_in,rep_v_out,ncol,nblock,globenv)
    TYPE(real_matrix_type), POINTER   :: sparse_matrix
    REAL(wp), DIMENSION(:,:), POINTER :: rep_v_in,rep_v_out
    INTEGER, INTENT(IN)               :: ncol,nblock
    TYPE(global_environment_type), INTENT(IN)      :: globenv

!   *** Local parameters ***
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE sparse_times_replicated (MODULE sparse_matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: iblock_col,iblock_row,icol,icol_block,irow,irow_block,istat,&
               nrow,sbncol,sbnrow,  handle
    REAL(wp) :: flops

    REAL(wp), DIMENSION(:,:), POINTER :: sparse_block

!   ---------------------------------------------------------------------------

!   *** Check the association status of the input matrix ***

    IF (.NOT.ASSOCIATED(sparse_matrix)) THEN
      CALL stop_program(routine,"The input matrix pointer is not associated")
    END IF

    if (sparse_matrix%symmetry .ne. "symmetric") then
      CALL stop_program(routine,"not tested "//sparse_matrix%symmetry)
    endif

    if (sparse_matrix%nrow .ne. sparse_matrix%ncol) then
      CALL stop_program(routine,"not implemented")
    endif
    if (sparse_matrix%ncol .ne. size(rep_v_in,2)) then
      CALL stop_program(routine,"n x n * n x k ?")
    endif

    CALL dcopy(sparse_matrix%ncol*nblock,0.0_wp,0,rep_v_out,1)

    CALL timeset("sparse_times_replicated_local","I","Gflops",handle)
    flops=0.0_wp

!   *** Traverse all block nodes of the sparse matrix ***

    DO iblock_row=1,sparse_matrix%nblock_row

      block_node => first_block_node(sparse_matrix,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            block=sparse_block)

        ! write(17+globenv%mepos,*) maxval(abs(sparse_block))

        sbnrow=sparse_matrix%last_row(iblock_row)-sparse_matrix%first_row(iblock_row)+1
        sbncol=sparse_matrix%last_col(iblock_col)-sparse_matrix%first_col(iblock_col)+1

        flops=flops+2.0_wp*sbnrow*sbncol*ncol
        CALL DGEMM('N','T',ncol,sbnrow,sbncol,1.0_wp, &
                           rep_v_in (1,sparse_matrix%first_col(iblock_col)), &
                           nblock, &
                           sparse_block(1,1),sbnrow, &
                    1.0_wp,rep_v_out(1,sparse_matrix%first_row(iblock_row)), &
                           nblock)

        if (iblock_col .ne. iblock_row .and. sparse_matrix%symmetry.eq."symmetric") then

           flops=flops+2.0_wp*sbnrow*sbncol*ncol
           CALL DGEMM('N','N',ncol,sbncol,sbnrow,1.0_wp, &
                           rep_v_in (1,sparse_matrix%first_row(iblock_row)), &
                           nblock, &
                           sparse_block(1,1),sbnrow, &
                    1.0_wp,rep_v_out(1,sparse_matrix%first_col(iblock_col)), &
                           nblock)
        endif

        block_node => next_block_node(block_node)

      END DO

    END DO
    flops=flops*1E-9_wp
    CALL timestop(flops,handle)

#if defined(__parallel)
    CALL mp_sum(rep_v_out,globenv%group)
#endif

  END SUBROUTINE sparse_times_replicated

! *****************************************************************************
! intended to compute sparse=sparse+v*g^T
! for all the local blocks of the sparse matrix
! same comment as sparse_times_replicated, v,g stored as rows instead of cols
! Joost VandeVondele july 2002
! *****************************************************************************
  SUBROUTINE sparse_plus_rep_rept(sparse_matrix,rep_v,rep_g,ncol,nblock,globenv)
    TYPE(real_matrix_type), POINTER   :: sparse_matrix
    REAL(wp), DIMENSION(:,:), POINTER :: rep_v,rep_g
    INTEGER, INTENT(IN)               :: ncol,nblock
    TYPE(global_environment_type), INTENT(IN)      :: globenv

!   *** Local parameters ***
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE sparse_plus_rep_rept(MODULE sparse_matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: iblock_col,iblock_row,icol,icol_block,irow,irow_block,istat,&
               nrow,sbncol,sbnrow,  handle
    REAL(wp) :: flops

    REAL(wp), DIMENSION(:,:), POINTER :: sparse_block

!   ---------------------------------------------------------------------------

!   *** Check the association status of the input matrix ***

    IF (.NOT.ASSOCIATED(sparse_matrix)) THEN
      CALL stop_program(routine,"The input matrix pointer is not associated")
    END IF

    if (sparse_matrix%nrow .ne. sparse_matrix%ncol) then
      CALL stop_program(routine,"not allowed")
    endif
    CALL timeset("sparse_plus_rep_rept","I","Gflops",handle)
    flops=0.0_wp

!   *** Traverse all block nodes of the sparse matrix ***

    DO iblock_row=1,sparse_matrix%nblock_row

      block_node => first_block_node(sparse_matrix,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            block=sparse_block)

        icol=sparse_matrix%first_col(iblock_col)
        irow=sparse_matrix%first_row(iblock_row)
        sbnrow=sparse_matrix%last_row(iblock_row)-irow+1
        sbncol=sparse_matrix%last_col(iblock_col)-icol+1

        flops=flops+2.0_wp*sbnrow*sbncol*ncol
        CALL DGEMM('T','N',sbnrow,sbncol,ncol,1.0_wp, &
                           rep_v(1,irow),nblock, &
                           rep_g(1,icol),nblock, &
                    1.0_wp,sparse_block,sbnrow)

        block_node => next_block_node(block_node)

      END DO

    END DO
    flops=flops*1E-9_wp
    CALL timestop(flops,handle)

  END SUBROUTINE sparse_plus_rep_rept

! *****************************************************************************

! *****************************************************************************

  SUBROUTINE copy_real_matrix(source,target)

!   Purpose: Copy the real matrix source to the existing real matrix target.

!   History: - Creation (07.07.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER :: source,target

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_real_matrix (MODULE sparse_matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: source_block_node,&
                                           target_block_node

    INTEGER :: first_col,first_row,handle,&
               iblock_col,iblock_row,icol,irow,istat,&
               jblock_col,jblock_row,jcol,jrow,last_col,last_row,&
               source_first_col,source_first_row,&
               source_last_col,source_last_row,&
               target_first_col,target_first_row,&
               target_last_col,target_last_row

!   ---------------------------------------------------------------------------

!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(source)) THEN
      CALL stop_program(routine,"The source matrix pointer is not associated")
    END IF

    CALL timeset("copy_real_matrix","I","",handle)

!   *** Check the association status of the target matrix ***

    IF (ASSOCIATED(target)) THEN

      DO iblock_row=1,source%nblock_row

        source_first_row = source%first_row(iblock_row)
        source_last_row = source%last_row(iblock_row)

        source_block_node => source%block_list(iblock_row)%first_block_node

        DO WHILE (ASSOCIATED(source_block_node))

          iblock_col = source_block_node%block_col

          source_first_col = source%first_col(iblock_col)
          source_last_col = source%last_col(iblock_col)

          DO jblock_row=1,target%nblock_row

            target_first_row = target%first_row(jblock_row)
            target_last_row = target%last_row(jblock_row)

            IF ((target_first_row <= source_last_row).AND.&
                (target_last_row >= source_first_row)) THEN

              target_block_node =>&
                target%block_list(jblock_row)%first_block_node

              DO WHILE (ASSOCIATED(target_block_node))

                jblock_col = target_block_node%block_col

                target_first_col = target%first_col(jblock_col)
                target_last_col = target%last_col(jblock_col)

                IF ((target_first_col <= source_last_col).AND.&
                    (target_last_col >= source_first_col)) THEN

                  first_col = MAX(source_first_col,target_first_col)
                  first_row = MAX(source_first_row,target_first_row)

                  last_col = MIN(source_last_col,target_last_col)
                  last_row = MIN(source_last_row,target_last_row)

                  jcol = first_col - target_first_col + 1

                  DO icol=first_col-source_first_col+1,&
                          last_col-source_first_col+1
                    jrow = first_row - target_first_row + 1
                    DO irow=first_row-source_first_row+1,&
                            last_row-source_first_row+1
                      target_block_node%block(jrow,jcol) =&
                        source_block_node%block(irow,icol)
                      jrow = jrow + 1
                    END DO
                    jcol = jcol + 1
                  END DO

                END IF

                target_block_node => target_block_node%next_block_node

              END DO

            END IF

          END DO

          source_block_node => source_block_node%next_block_node

        END DO

      END DO

    ELSE

      CALL stop_program(routine,&
                        "The target matrix pointer is not associated. "//&
                        "Use replicate_matrix instead of copy_matrix")

    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE copy_real_matrix

! *****************************************************************************

  SUBROUTINE deallocate_real_matrix(matrix)

!   Purpose: Deallocate a real matrix at the real_matrix_type level.

!   History: - Creation (11.07.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER :: matrix

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_real_matrix (MODULE sparse_matrix_types)"

!   *** Local variables ***

    INTEGER :: iblock_row,istat

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(matrix)) RETURN

    DO iblock_row=1,matrix%nblock_row
      CALL deallocate_matrix_row(matrix,iblock_row)
    END DO

    DEALLOCATE (matrix%block_list,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,matrix%name//" (block_list)")

    DEALLOCATE (matrix%first_row,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,matrix%name//" (first_row)")

    DEALLOCATE (matrix%last_row,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,matrix%name//" (last_row)")

    DEALLOCATE (matrix%first_col,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,matrix%name//" (first_col)")

    DEALLOCATE (matrix%last_col,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,matrix%name//" (last_col)")

    DEALLOCATE (matrix,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,matrix%name)

  END SUBROUTINE deallocate_real_matrix

! *****************************************************************************

  SUBROUTINE deallocate_real_matrix_row(matrix,block_row)

!   Purpose: Deallocate a matrix row.

!   History: - Creation (30.11.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER :: matrix
    INTEGER, INTENT(IN)             :: block_row

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_real_matrix_row (MODULE sparse_matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: current_block_node,next_block_node

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    current_block_node => matrix%block_list(block_row)%first_block_node

    DO WHILE (ASSOCIATED(current_block_node))
      next_block_node => current_block_node%next_block_node
      DEALLOCATE (current_block_node%block,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,matrix%name//" (block)")
      DEALLOCATE (current_block_node,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,matrix%name//" (block_node)")
      current_block_node => next_block_node
    END DO

    NULLIFY (matrix%block_list(block_row)%first_block_node)

  END SUBROUTINE deallocate_real_matrix_row

! *****************************************************************************

  SUBROUTINE deallocate_real_matrix_set(matrix_set)

!   Purpose: Deallocate a real matrix set.

!   History: - Creation (13.03.2002, MK)

!   ***************************************************************************

    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: matrix_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_real_matrix_set (MODULE sparse_matrix_types)"

!   *** Local variables ***

    INTEGER :: imatrix,istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(matrix_set)) THEN
      DO imatrix=1,SIZE(matrix_set)
        CALL deallocate_matrix(matrix_set(imatrix)%matrix)
      END DO
    ELSE
      CALL stop_program(routine,"The matrix_set pointer is not associated")
    END IF

    DEALLOCATE (matrix_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"matrix_set")

  END SUBROUTINE deallocate_real_matrix_set

! *****************************************************************************

  FUNCTION find_real_block_node(matrix,block_row,block_col) RESULT(block_node)

!   Purpose: Return a pointer to the requested block node.

!   History: - Creation (23.06.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER :: matrix
    INTEGER, INTENT(IN)             :: block_col,block_row

    TYPE(real_block_node_type), POINTER :: block_node

!   ---------------------------------------------------------------------------

    block_node => matrix%block_list(block_row)%first_block_node

    DO WHILE (ASSOCIATED(block_node))
      IF (block_node%block_col == block_col) EXIT
      block_node => block_node%next_block_node
    END DO

  END FUNCTION find_real_block_node

! *****************************************************************************

  FUNCTION first_real_block_node(matrix,block_row) RESULT(first_block_node)

!   Purpose: Return a pointer to the first block node of a block list.

!   History: - Creation (23.06.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER :: matrix
    INTEGER, INTENT(IN)             :: block_row

    TYPE(real_block_node_type), POINTER :: first_block_node

!   ---------------------------------------------------------------------------

    first_block_node => matrix%block_list(block_row)%first_block_node

  END FUNCTION first_real_block_node

! *****************************************************************************

  SUBROUTINE get_matrix_info(matrix,matrix_name,matrix_symmetry,&
                             nblock_row,nblock_col,nrow,ncol,&
                             first_row,last_row,first_col,last_col,&
                             nblock_allocated,nelement_allocated)

!   Purpose: Return the requested matrix information.

!   History: - Creation (10.07.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: matrix
    CHARACTER(LEN=60), OPTIONAL, INTENT(OUT) :: matrix_name
    CHARACTER(LEN=40), OPTIONAL, INTENT(OUT) :: matrix_symmetry
    INTEGER, OPTIONAL, INTENT(OUT)           :: nblock_allocated,nblock_row,&
                                                nblock_col,ncol,&
                                                nelement_allocated,nrow
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: first_col,first_row,&
                                                last_col,last_row

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: iblock_row

    REAL(wp), DIMENSION(:,:), POINTER :: block

!   ---------------------------------------------------------------------------

    IF (PRESENT(matrix_name)) matrix_name = matrix%name
    IF (PRESENT(matrix_symmetry)) matrix_symmetry = matrix%symmetry
    IF (PRESENT(nblock_row)) nblock_row = matrix%nblock_row
    IF (PRESENT(nblock_col)) nblock_col = matrix%nblock_col
    IF (PRESENT(nrow)) nrow = matrix%nrow
    IF (PRESENT(ncol)) ncol = matrix%ncol
    IF (PRESENT(first_row)) first_row => matrix%first_row
    IF (PRESENT(last_row)) last_row => matrix%last_row
    IF (PRESENT(first_col)) first_col => matrix%first_col
    IF (PRESENT(last_col)) last_col => matrix%last_col

    IF (PRESENT(nblock_allocated)) THEN
      nblock_allocated = 0
      DO iblock_row=1,matrix%nblock_row
        block_node => first_block_node(matrix,iblock_row)
        DO WHILE (ASSOCIATED(block_node))
          nblock_allocated = nblock_allocated + 1
          block_node => next_block_node(block_node)
        END DO
      END DO
    END IF

    IF (PRESENT(nelement_allocated)) THEN
      nelement_allocated = 0
      DO iblock_row=1,matrix%nblock_row
        block_node => first_block_node(matrix,iblock_row)
        DO WHILE (ASSOCIATED(block_node))
          CALL get_block_node(block_node=block_node,block=block)
          nelement_allocated = nelement_allocated + SIZE(block)
          block_node => next_block_node(block_node)
        END DO
      END DO
    END IF

  END SUBROUTINE get_matrix_info

! *****************************************************************************

  SUBROUTINE get_real_block_node(block_node,block_col,block)

!   Purpose: Get block node data set.

!   History: - Creation (28.05.2000, MK)

!   ***************************************************************************

    TYPE(real_block_node_type), INTENT(IN)      :: block_node
    INTEGER, OPTIONAL, INTENT(OUT)              :: block_col
    REAL(wp), DIMENSION(:,:), OPTIONAL, POINTER :: block

!   ---------------------------------------------------------------------------

    IF (PRESENT(block_col)) block_col = block_node%block_col
    IF (PRESENT(block)) block => block_node%block

  END SUBROUTINE get_real_block_node

! *****************************************************************************

  SUBROUTINE get_real_matrix_block(matrix,block_row,block_col,&
                                   first_row,last_row,first_col,last_col,&
                                   block_node,block)

!   Purpose: Get block node data set.

!   History: - Creation (28.05.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER               :: matrix
    TYPE(real_block_node_type), OPTIONAL, POINTER :: block_node
    INTEGER, INTENT(IN)                           :: block_col,block_row
    INTEGER, OPTIONAL, INTENT(OUT)                :: first_col,first_row,&
                                                     last_col,last_row
    REAL(wp), DIMENSION(:,:), OPTIONAL, POINTER   :: block

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: current_block_node

!   ---------------------------------------------------------------------------

    IF (PRESENT(first_row)) first_row = matrix%first_row(block_row)
    IF (PRESENT(last_row)) last_row = matrix%last_row(block_row)

    IF (PRESENT(first_col)) first_col = matrix%first_col(block_col)
    IF (PRESENT(last_col)) last_col = matrix%last_col(block_col)

    current_block_node => find_real_block_node(matrix,block_row,block_col)

    IF (ASSOCIATED(current_block_node)) THEN
      IF (PRESENT(block_node)) block_node => current_block_node
      IF (PRESENT(block)) block => current_block_node%block
    ELSE
      IF (PRESENT(block_node)) NULLIFY (block_node)
      IF (PRESENT(block)) NULLIFY (block)
    END IF

  END SUBROUTINE get_real_matrix_block

! *****************************************************************************

  SUBROUTINE merge_real_matrices(a,b)

!   Purpose: Merge two sparse matrices a and b which have the same structure,
!            i.e. all blocks of b which are not present in a will be created
!            and blocks which are present in both matrices will be added.
!            b will be removed before exit.

!   History: - Creation (16.07.2002,MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER :: a,b

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE merge_real_matrices (MODULE sparse_matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: iblock_col,iblock_row

    REAL(wp), DIMENSION(:,:), POINTER :: a_block,b_block

!   ---------------------------------------------------------------------------

!   *** Check the association status of the matrices ***

    IF (.NOT.ASSOCIATED(a)) THEN
       CALL stop_program(routine,"The matrix pointer a is not associated")
    END IF

    IF (.NOT.ASSOCIATED(b)) THEN
       CALL stop_program(routine,"The matrix pointer b is not associated")
    END IF

    DO iblock_row=1,b%nblock_row

      block_node => first_block_node(b,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            block=b_block)

        CALL get_block_node(matrix=a,&
                            block_row=iblock_row,&
                            block_col=iblock_col,&
                            block=a_block)

        IF (ASSOCIATED(a_block)) THEN
          a_block(:,:) = a_block(:,:) + b_block(:,:)
        ELSE
          CALL add_real_matrix_block(matrix=a,&
                                     block_row=iblock_row,&
                                     block_col=iblock_col,&
                                     block=b_block)
        END IF

        block_node => block_node%next_block_node

      END DO

      CALL deallocate_real_matrix_row(matrix=b,block_row=iblock_row)

    END DO

    CALL deallocate_real_matrix(b)

  END SUBROUTINE merge_real_matrices

! *****************************************************************************

  FUNCTION next_real_block_node(block_node) RESULT(next_block_node)

!   Purpose: Return a pointer to the next block node of a block list.

!   History: - Creation (28.06.2000, MK)

!   ***************************************************************************

    TYPE(real_block_node_type), POINTER :: block_node

    TYPE(real_block_node_type), POINTER :: next_block_node

!   ---------------------------------------------------------------------------
    if (associated(block_node)) then
       next_block_node => block_node%next_block_node
    else
       CALL stop_program("sparse_matrix_types:next_real_block_node",&
            "The block_node pointer a is not associated")
       nullify(next_block_node)
    end if

  END FUNCTION next_real_block_node

! *****************************************************************************

  SUBROUTINE put_real_block_node(block_node,matrix,block_row,block_col,block)

!   Purpose: Set block node data set.

!   History: - Creation (27.07.2000, MK)

!   ***************************************************************************

    TYPE(real_block_node_type), POINTER         :: block_node
    TYPE(real_matrix_type), POINTER             :: matrix
    INTEGER, INTENT(IN)                         :: block_row
    INTEGER, OPTIONAL, INTENT(IN)               :: block_col
    REAL(wp), DIMENSION(:,:), OPTIONAL, POINTER :: block

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE put_real_block_node (MODULE sparse_matrix_types)"

!   *** Local variables ***

    INTEGER :: ncol,nrow

!   ---------------------------------------------------------------------------

    IF (PRESENT(block_col)) block_node%block_col = block_col

    IF (PRESENT(block)) THEN

      nrow = matrix%last_row(block_row) - matrix%first_row(block_row) + 1
      ncol = matrix%last_col(block_col) - matrix%first_col(block_col) + 1

      IF (PRESENT(block)) THEN
        IF (ASSOCIATED(block)) THEN
          IF ((SIZE(block,1) == nrow).AND.(SIZE(block,2) == ncol)) THEN
            block_node%block(:,:) = block(:,:)
          ELSE
            CALL stop_program(routine,"Incompatible block dimensions")
          END IF
        ELSE
          CALL stop_program(routine,"Block is not associated")
        END IF
      END IF

    END IF

  END SUBROUTINE put_real_block_node

! *****************************************************************************

  SUBROUTINE put_real_matrix_block(matrix,block_row,block_col,block)

!   Purpose: Allocate and initialize a new block node or update an existing
!            block node.

!   History: - Creation (07.07.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER            :: matrix
    INTEGER, INTENT(IN)                        :: block_col,block_row
    REAL(wp), DIMENSION(:,:), OPTIONAL, POINTER:: block

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE put_real_matrix_block (MODULE sparse_matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: istat,ncol,nrow

!   ---------------------------------------------------------------------------

!   *** Check, if the requested block node already exists ***

    block_node => find_real_block_node(matrix,block_row,block_col)

    IF (ASSOCIATED(block_node)) THEN

!     *** Calculate the block dimensions ***

      nrow = matrix%last_row(block_row) - matrix%first_row(block_row) + 1
      ncol = matrix%last_col(block_col) - matrix%first_col(block_col) + 1

!     *** Update the data set of an existing block node ***

      IF (PRESENT(block)) THEN
        IF (ASSOCIATED(block)) THEN
          IF ((SIZE(block,1) == nrow).AND.(SIZE(block,2) == ncol)) THEN
            block_node%block(:,:) = block(:,:)
          ELSE
            CALL stop_program(routine,"Incompatible block dimensions")
          END IF
        ELSE
          CALL stop_program(routine,"Matrix block is not associated")
        END IF
      ELSE
        CALL stop_program(routine,&
                          "Block update is impossible, because no matrix "//&
                          "block was specified")
      END IF

    ELSE

!     *** Add a new block node ***

      IF (PRESENT(block)) THEN
        CALL add_block_node(matrix,block_row,block_col,block)
      ELSE
        CALL add_block_node(matrix,block_row,block_col)
      END IF

    END IF

  END SUBROUTINE put_real_matrix_block

! *****************************************************************************

  SUBROUTINE replicate_real_matrix(source,target,target_name,allocate_blocks)

!   Purpose: Replicate the existing matrix source. The replicated matrix is
!            target with the name target_name.
!            by default copies the data blocks,
!            if optional allocate_blocks=false then no blocks are allocated

!   History: - Creation (17.11.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER :: source,target
    CHARACTER(LEN=*), INTENT(IN)    :: target_name
    LOGICAL, OPTIONAL, INTENT(IN)   :: allocate_blocks

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE replicate_real_matrix (MODULE sparse_matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: iblock_col,iblock_row
    LOGICAL :: make_data_blocks

    REAL(wp), DIMENSION(:,:), POINTER :: block

!   ---------------------------------------------------------------------------

!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(source)) THEN
      CALL stop_program(routine,"The source matrix pointer is not associated")
    END IF

    IF (ASSOCIATED(target)) CALL deallocate_matrix(target)

!   *** Allocate a new matrix structure ***

    CALL allocate_matrix(matrix=target,&
                         nrow=source%nrow,&
                         ncol=source%ncol,&
                         nblock_row=source%nblock_row,&
                         nblock_col=source%nblock_col,&
                         first_row=source%first_row(:),&
                         last_row=source%last_row(:),&
                         first_col=source%first_col(:),&
                         last_col=source%last_col(:),&
                         matrix_name=target_name,&
                         matrix_symmetry=source%symmetry)

    IF (PRESENT(allocate_blocks)) THEN
      make_data_blocks = allocate_blocks
    ELSE
      make_data_blocks = .TRUE.
    END IF

!   *** Initialize all block nodes ***

    IF (make_data_blocks) THEN

      DO iblock_row=1,source%nblock_row

        block_node => first_block_node(source,iblock_row)

        DO WHILE (ASSOCIATED(block_node))

          CALL get_block_node(block_node=block_node,&
                              block_col=iblock_col,&
                              block=block)

          CALL add_block_node(matrix=target,&
                              block_row=iblock_row,&
                              block_col=iblock_col,&
                              block=block)

          block_node => next_block_node(block_node)

        END DO

      END DO

    ENDIF

  END SUBROUTINE replicate_real_matrix

! *****************************************************************************

  SUBROUTINE replicate_real_matrix_structure(source,target,target_name,&
                                             target_symmetry)

!   Purpose: Replicate the matrix structure of the existing matrix source. The
!            new matrix target with the name target_name has the same
!            structure.

!   History: - Creation (17.11.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER        :: source,target
    CHARACTER(LEN=*), INTENT(IN)           :: target_name
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: target_symmetry

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE replicate_real_matrix_structure (MODULE sparse_matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node
    CHARACTER(LEN=40)                   :: matrix_symmetry
    INTEGER                             :: iblock_col,iblock_row

!   ---------------------------------------------------------------------------

!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(source)) THEN
      CALL stop_program(routine,"The source matrix pointer is not associated")
    END IF

    IF (PRESENT(target_symmetry)) THEN
      matrix_symmetry = target_symmetry
    ELSE
      matrix_symmetry = source%symmetry
    END IF

!   *** Allocate a new matrix structure ***

    CALL allocate_matrix(matrix=target,&
                         nrow=source%nrow,&
                         ncol=source%ncol,&
                         nblock_row=source%nblock_row,&
                         nblock_col=source%nblock_col,&
                         first_row=source%first_row(:),&
                         last_row=source%last_row(:),&
                         first_col=source%first_col(:),&
                         last_col=source%last_col(:),&
                         matrix_name=target_name,&
                         matrix_symmetry=matrix_symmetry)

!   *** Initialize all block nodes ***

    DO iblock_row=1,source%nblock_row

      block_node => first_block_node(source,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col)

        CALL add_block_node(matrix=target,&
                            block_row=iblock_row,&
                            block_col=iblock_col)

        block_node => next_block_node(block_node)

      END DO

    END DO

  END SUBROUTINE replicate_real_matrix_structure

! *****************************************************************************

  SUBROUTINE set_real_matrix(matrix,value)

!   Purpose: Set all elements of matrix to value.

!   History: - Creation (10.04.2002, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER :: matrix
    REAL(wp), INTENT(IN)            :: value

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE set_real_matrix (MODULE sparse_matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: iblock_row

    REAL(wp), DIMENSION(:,:), POINTER :: block

!   ---------------------------------------------------------------------------

!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(matrix)) THEN
      CALL stop_program(routine,"The matrix pointer is not associated")
    END IF

!   *** Set matrix elements to value ***

    DO iblock_row=1,matrix%nblock_row

      block_node => first_block_node(matrix,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block=block)

        block(:,:) = value

        block_node => block_node%next_block_node

      END DO

    END DO

  END SUBROUTINE set_real_matrix

! *****************************************************************************

  SUBROUTINE symmetrise_diagonal_blocks(matrix)

!   Purpose: Symmetrise the diagonal blocks of matrix.

!   History: - Creation (13.06.2001, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER :: matrix

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE symmetrise_diagonal_blocks (MODULE sparse_matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: iblock_col,iblock_row,icol,irow

    REAL(wp), DIMENSION(:,:), POINTER :: block

!   ---------------------------------------------------------------------------

!   *** Check the association status of the input matrix ***

    IF (.NOT.ASSOCIATED(matrix)) THEN
      CALL stop_program(routine,"The input matrix pointer is not associated")
    END IF

!   *** Traverse all block nodes ***

    DO iblock_row=1,matrix%nblock_row

      block_node => first_block_node(matrix,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            block=block)

        IF (iblock_row == iblock_col) THEN
          DO irow=1,SIZE(block,1)
            DO icol=irow+1,SIZE(block,2)
              block(icol,irow) = block(irow,icol)
            END DO
          END DO
        END IF

        block_node => next_block_node(block_node)

      END DO

    END DO

  END SUBROUTINE symmetrise_diagonal_blocks

! *****************************************************************************

  SUBROUTINE transfer_real_matrix(source,target)

!   Purpose: Transfer source matrix to target matrix. Both matrices have the
!            same structure.

!   History: - Creation (10.04.2002, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER :: source,target

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE transfer_real_matrix (MODULE sparse_matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: iblock_col,iblock_row

    REAL(wp), DIMENSION(:,:), POINTER :: source_block,target_block

!   ---------------------------------------------------------------------------

!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(source)) THEN
      CALL stop_program(routine,"The source matrix pointer is not associated")
    END IF

    IF (.NOT.ASSOCIATED(target)) THEN
      CALL stop_program(routine,"The target matrix pointer is not associated")
    END IF

!   *** Copy matrix elements ***

    DO iblock_row=1,source%nblock_row

      block_node => first_block_node(source,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            block=source_block)

        CALL get_block_node(matrix=target,&
                            block_row=iblock_row,&
                            block_col=iblock_col,&
                            block=target_block)

        target_block(:,:) = source_block(:,:)

        block_node => block_node%next_block_node

      END DO

    END DO

  END SUBROUTINE transfer_real_matrix

! *****************************************************************************

!!****f* sparse_matrix_types/remove_real_matrix_block [1.0] *
!!
!!   NAME
!!     remove_real_matrix_block
!!
!!   FUNCTION
!!     removes a block from the matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix you remove the block from
!!     - block_row: the row of the block to remove
!!     - block_col: the column of the block to remove
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     7.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE remove_real_matrix_block(matrix,block_row, block_col, error)
    USE cp_log_handling, ONLY: cp_to_string, cp_failure_level, cp_warning_level
    USE cp_error_handling, ONLY: cp_assert, cp_error_type
    TYPE(real_matrix_type), POINTER :: matrix
    INTEGER, INTENT(in) :: block_col, block_row
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    INTEGER :: stat
    CHARACTER(len=*), PARAMETER :: routineN='remove_real_matrix_block',&
         routineP=moduleN//':'//routineN
    TYPE(real_block_node_type), POINTER :: block_node, old_block
    NULLIFY(block_node,old_block)
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       block_node => matrix%block_list(block_row)%first_block_node
       DO WHILE (ASSOCIATED(block_node))
          IF (block_node%block_col == block_col) EXIT
          old_block => block_node
          block_node => block_node%next_block_node
       END DO
       IF (ASSOCIATED(block_node)) THEN
          IF (ASSOCIATED(old_block)) THEN
             old_block%next_block_node => block_node%next_block_node
          ELSE
             matrix%block_list(block_row)%first_block_node => &
                  block_node%next_block_node
          END IF
          DEALLOCATE (block_node%block,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          DEALLOCATE (block_node,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE remove_real_matrix_block
!***************************************************************************

END MODULE sparse_matrix_types
