!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Define the matrix data types.
!> \par History
!>      none
!> \author MK (23.06.2000)
! *****************************************************************************
MODULE sparse_matrix_types
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_retain,&
                                             distribution_2d_type
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sum
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sparse_matrix_types'
  INTEGER, SAVE, PRIVATE :: last_matrix_id=0

! *** Define a real matrix data type ***

! *****************************************************************************
  TYPE real_block_node_type
    PRIVATE
    TYPE(real_block_node_type), POINTER :: next_block_node
    INTEGER                             :: block_col
    REAL(KIND = dp), DIMENSION(:,:), POINTER   :: BLOCK
  END TYPE real_block_node_type

! *****************************************************************************
  TYPE real_block_list_type
    PRIVATE
    TYPE(real_block_node_type), POINTER :: first_block_node
    INTEGER                             :: nblock_node
    TYPE(real_block_node_type), POINTER :: last_used_block_node ! not last in list
  END TYPE real_block_list_type

! *****************************************************************************
  TYPE real_matrix_type
!    PRIVATE
    TYPE(real_block_list_type), DIMENSION(:), POINTER :: block_list
    CHARACTER(LEN=80)                                 :: name
    CHARACTER(LEN=40)                                 :: symmetry
    INTEGER                                           :: nblock_col,&
                                                         nblock_row,&
                                                         ncol,nrow
    INTEGER :: id_nr, ref_count!, print_count
    INTEGER, DIMENSION(:), POINTER                    :: first_col,&
                                                         first_row,&
                                                         last_col,&
                                                         last_row
    !TYPE(cp_sm_struct_type), pointer                 :: matrix_struct
    TYPE(distribution_2d_type), POINTER               :: distribution_2d

    ! the sparsity of type of this matrix is of this kind
    ! this is not yet strictly enforced
    INTEGER                                           ::  sparsity_id
  END TYPE real_matrix_type

! *****************************************************************************
  TYPE real_matrix_p_type
    TYPE(real_matrix_type), POINTER :: matrix
  END TYPE real_matrix_p_type

! *** Public data types ***

  PUBLIC :: real_block_node_type,&
            real_matrix_p_type,&
            real_matrix_type, &
            real_block_list_type

! *** Public subroutines ***

  PUBLIC :: add_block_node,&
            add_matrices,&
            add_1d_block_node,&
            allocate_matrix,&
            allocate_matrix_set,&
            copy_matrix, transfer_matrix, & ! what is the difference ?
            checksum_matrix, &
            copy_local_sm_to_replicated_fm,&
            deallocate_matrix,&
            deallocate_matrix_set,&
            get_block_node,&
            get_matrix_diagonal,&
            get_matrix_info,&
            put_block_node,&
            replicate_matrix, replicate_matrix_structure,& ! should be put together
            scale_matrix,&
            set_matrix,&
            set_matrix_diagonal,&
            sparse_times_local, &
            sparse_plus_loc_loct, &
            remove_block_node, &
            write_blocks_maxabsval, get_real_matrix_occupation,& ! should be removed
            cp_sm_scale_and_add,&
            cp_sm_get_id_nr,&
            cp_sm_set,&
            cp_sm_sm_trace,&
            cp_sm_get_block_list,&
            cp_sm_maxabsval,&
            create_dbm_from_sm,&
            deallocate_dbm_from_sm,&
            cp_sm_add_constant, &
            mpsum_dbm

! *** Public functions ***

  PUBLIC :: first_block_node,&
            next_block_node,&
            find_block_node

  INTERFACE add_block_node
    MODULE PROCEDURE add_real_matrix_block
  END INTERFACE

  INTERFACE add_matrices
    MODULE PROCEDURE add_real_matrices
  END INTERFACE

  INTERFACE allocate_matrix
    MODULE PROCEDURE allocate_real_matrix
  END INTERFACE

  INTERFACE allocate_matrix_set
    MODULE PROCEDURE allocate_real_matrix_set
    MODULE PROCEDURE allocate_real_matrix_set_2d
  END INTERFACE

  INTERFACE copy_matrix
    MODULE PROCEDURE copy_real_matrix
  END INTERFACE

  INTERFACE checksum_matrix
     MODULE PROCEDURE checksum_real_matrix
  END INTERFACE

  INTERFACE deallocate_matrix
    MODULE PROCEDURE deallocate_real_matrix
  END INTERFACE

  INTERFACE deallocate_matrix_row
    MODULE PROCEDURE deallocate_real_matrix_row
  END INTERFACE

  INTERFACE deallocate_matrix_set
    MODULE PROCEDURE deallocate_real_matrix_set
    MODULE PROCEDURE deallocate_real_matrix_set_2d
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

  INTERFACE get_matrix_diagonal
    MODULE PROCEDURE get_real_matrix_diagonal
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

  INTERFACE scale_matrix
    MODULE PROCEDURE scale_real_matrix
  END INTERFACE

  INTERFACE set_matrix
    MODULE PROCEDURE set_real_matrix
  END INTERFACE

  INTERFACE set_matrix_diagonal
    MODULE PROCEDURE set_real_matrix_diagonal
  END INTERFACE

  INTERFACE transfer_matrix
    MODULE PROCEDURE transfer_real_matrix
  END INTERFACE

  INTERFACE remove_block_node
     MODULE PROCEDURE remove_real_matrix_block
  END INTERFACE

  INTERFACE cp_sm_sm_trace
     MODULE PROCEDURE cp_sm_sm_trace_scalar, cp_sm_sm_trace_vector
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief scales the matrix and add another matrix
!>      matrix_a=alpha*matrix_a+beta*matrix_b
!>      filters the result with the structure (sparsity) of matrix_a
!> \param matrix_a a sparse matrix (inout)
!> \param alpha scaling factor of matrix_a (defaults to 1.0)
!> \param matrix_b a sparse matrix (in)
!> \param beta scaling factor of matrix_b (defaults to 1.0)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      should be rewitten avoiding the get_block_node taht calls find_block_node
!> \par History
!>      2.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_sm_scale_and_add(matrix_a,alpha,matrix_b,beta,error)
    TYPE(real_matrix_type), POINTER          :: matrix_a
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: alpha
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: matrix_b
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: beta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_scale_and_add', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_col, iblock_row
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_alpha, my_beta
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a_block, b_block
    TYPE(real_block_node_type), POINTER      :: block_node

  CALL timeset(routineN,handle)

  NULLIFY(a_block,b_block,block_node)
  failure=.FALSE.
  my_beta=1.0_dp
  my_alpha=1.0_dp
  IF (PRESENT(alpha)) my_alpha=alpha
  IF (PRESENT(beta)) my_beta=beta

  IF (.not.PRESENT(matrix_b)) my_beta=0.0_dp
  IF (my_beta==0) THEN
     IF (my_alpha/=1.0_dp) THEN
        CALL scale_real_matrix(matrix_a,my_alpha)
     END IF
  ELSE
     CPPrecondition(ASSOCIATED(matrix_a),cp_failure_level,routineP,error,failure)
     CPPrecondition(ASSOCIATED(matrix_b),cp_failure_level,routineP,error,failure)
    IF (my_alpha == 1.0_dp) THEN ! make common case fast
       DO iblock_row=1,matrix_a%nblock_row
          block_node => first_block_node(matrix_a,iblock_row)
          DO WHILE (ASSOCIATED(block_node))

             CALL get_block_node(block_node=block_node,&
                                 block_col=iblock_col,&
                                 BLOCK=a_block)
             CPPrecondition(ASSOCIATED(a_block),cp_failure_level,routineP,error,failure)
             CALL get_block_node(matrix=matrix_b,&
                                 block_row=iblock_row,&
                                 block_col=iblock_col,&
                                 BLOCK=b_block)

             IF (ASSOCIATED(b_block)) THEN
                CALL daxpy(SIZE(a_block,1)*SIZE(a_block,2),&
                     my_beta, b_block(1,1),1,&
                     a_block(1,1),1)
             END IF

             block_node => block_node%next_block_node
          END DO
       END DO
    ELSE
       DO iblock_row=1,matrix_a%nblock_row
          block_node => first_block_node(matrix_a,iblock_row)
          DO WHILE (ASSOCIATED(block_node))
             CALL get_block_node(block_node=block_node,&
                                 block_col=iblock_col,&
                                 BLOCK=a_block)
             CPPrecondition(ASSOCIATED(a_block),cp_failure_level,routineP,error,failure)
             CALL get_block_node(matrix=matrix_b,&
                                 block_row=iblock_row,&
                                 block_col=iblock_col,&
                                 BLOCK=b_block)

             IF (ASSOCIATED(b_block)) THEN
                a_block(:,:) = my_alpha*a_block(:,:) + my_beta*b_block(:,:)
             ELSE
                CALL dscal(SIZE(a_block,1)*SIZE(a_block,2),&
                     my_alpha,a_block(1,1),1)
             END IF

             block_node => block_node%next_block_node
          END DO
       END DO
    END IF

  END IF
  CALL timestop(handle)
END SUBROUTINE cp_sm_scale_and_add

! *****************************************************************************
!> \brief  Calculate the sum of the sparse matrices a and b which have the
!>           same structure.
!>           c = alpha*a + beta*b 
!> \author  MK
!> \date    10.04.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE add_real_matrices(c,alpha,a,beta,b)

    TYPE(real_matrix_type), POINTER          :: c
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(real_matrix_type), POINTER          :: a
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(real_matrix_type), POINTER          :: b

    CHARACTER(len=*), PARAMETER :: routineN = 'add_real_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iblock_col, iblock_row
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a_block, b_block, c_block
    TYPE(real_block_node_type), POINTER      :: block_node

!   *** Check the association status of the matrices ***

    IF (.NOT.ASSOCIATED(a)) THEN
      CALL stop_program(routineP,"The matrix pointer a is not associated")
    END IF

    IF (.NOT.ASSOCIATED(b)) THEN
      CALL stop_program(routineP,"The matrix pointer b is not associated")
    END IF

    IF (.NOT.ASSOCIATED(c)) THEN
      CALL stop_program(routineP,"The matrix pointer c is not associated")
    END IF

!   *** Add matrix elements ***

    IF ((alpha == 1.0_dp).AND.(beta == 1.0_dp)) THEN

      DO iblock_row=1,c%nblock_row

        block_node => first_block_node(c,iblock_row)

        DO WHILE (ASSOCIATED(block_node))

          CALL get_block_node(block_node=block_node,&
                              block_col=iblock_col,&
                              BLOCK=c_block)

          CALL get_block_node(matrix=a,&
                              block_row=iblock_row,&
                              block_col=iblock_col,&
                              BLOCK=a_block)

          CALL get_block_node(matrix=b,&
                              block_row=iblock_row,&
                              block_col=iblock_col,&
                              BLOCK=b_block)

          IF (ASSOCIATED(a_block)) THEN
            IF (ASSOCIATED(b_block)) THEN
              c_block(:,:) = alpha*a_block(:,:) + beta*b_block(:,:)
            ELSE
              c_block(:,:) = alpha*a_block(:,:)
            END IF
          ELSE
            IF (ASSOCIATED(b_block)) THEN
              c_block(:,:) = beta*b_block(:,:)
            ELSE
              c_block(:,:) = 0.0_dp
            END IF
          END IF

          block_node => block_node%next_block_node

        END DO

      END DO

    ELSE

      DO iblock_row=1,c%nblock_row

        block_node => first_block_node(c,iblock_row)

        DO WHILE (ASSOCIATED(block_node))

          CALL get_block_node(block_node=block_node,&
                              block_col=iblock_col,&
                              BLOCK=c_block)

          CALL get_block_node(matrix=a,&
                              block_row=iblock_row,&
                              block_col=iblock_col,&
                              BLOCK=a_block)

          CALL get_block_node(matrix=b,&
                              block_row=iblock_row,&
                              block_col=iblock_col,&
                              BLOCK=b_block)

          IF (ASSOCIATED(a_block)) THEN
            IF (ASSOCIATED(b_block)) THEN
              c_block(:,:) = alpha*a_block(:,:) + beta*b_block(:,:)
            ELSE
              c_block(:,:) = alpha*a_block(:,:)
            END IF
          ELSE
            IF (ASSOCIATED(b_block)) THEN
              c_block(:,:) = beta*b_block(:,:)
            ELSE
              c_block(:,:) = 0.0_dp
            END IF
          END IF

          block_node => block_node%next_block_node

        END DO

      END DO

    END IF

  END SUBROUTINE add_real_matrices

! *****************************************************************************
!> \brief   Allocate and initialize a new block node.
!> \author  MK
!> \date    28.07.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE add_real_matrix_block(matrix,block_row,block_col,BLOCK, error)

    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row, block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: BLOCK
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_real_matrix_block', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cbc, cbr, istat, mypcol, &
                                                myprow, ncol, nrow, rbc, rbr
    LOGICAL                                  :: failure
    TYPE(real_block_node_type), POINTER      :: block_node, new_block_node

!   perform checks that we are only adding allowed blocks

    failure = .FALSE.
    IF (ASSOCIATED(matrix%distribution_2d)) THEN
       myprow = matrix%distribution_2d%blacs_env%mepos(1)
       mypcol = matrix%distribution_2d%blacs_env%mepos(2)
       rbr    = matrix%distribution_2d%row_distribution(block_row)
       rbc    = matrix%distribution_2d%row_distribution(block_col)
       cbr    = matrix%distribution_2d%col_distribution(block_row)
       cbc    = matrix%distribution_2d%col_distribution(block_col)
       SELECT CASE(matrix%symmetry)
       CASE("symmetric","antisymmetric")
         ! we allow for putting a block at i,j even if we officially only own j,i
         IF ( .NOT. ((rbr.EQ.myprow .AND. cbc.EQ.mypcol) .OR. &
                     (rbc.EQ.myprow .AND. cbr.EQ.mypcol))       ) THEN
            CPPostcondition(.FALSE.,cp_warning_level,routineP,error,failure)
         ENDIF
       CASE DEFAULT
         IF ( .NOT. ( rbr.EQ.myprow .AND. cbc.EQ.mypcol) ) THEN
            CPPostcondition(.FALSE.,cp_warning_level,routineP,error,failure)
         ENDIF
       END SELECT
    ENDIF

!   *** Calculate the block dimensions ***

    nrow = matrix%last_row(block_row) - matrix%first_row(block_row) + 1
    ncol = matrix%last_col(block_col) - matrix%first_col(block_col) + 1

    ALLOCATE (new_block_node,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"Block node of "//matrix%name,0)

!   *** Define the data set of the new block node ***

    new_block_node%block_col = block_col

    ALLOCATE (new_block_node%block(nrow,ncol),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineP,"new_block_node%block",dp_size*nrow*ncol)
    END IF

    new_block_node%block(:,:) = 0.0_dp
    IF (PRESENT(BLOCK)) THEN
      IF (ASSOCIATED(BLOCK)) THEN
        IF ((SIZE(BLOCK,1) == nrow).AND.(SIZE(BLOCK,2) == ncol)) THEN
          new_block_node%block(:,:) = BLOCK(:,:)
        ELSE
          CALL stop_program(routineP,"Incompatible block dimensions")
        END IF
      ELSE
        BLOCK => new_block_node%block
      END IF
    END IF

!   *** Link the new block node to the block list ***

    block_node => matrix%block_list(block_row)%first_block_node
    ! if the last used block node is a short cut, use it
    IF (ASSOCIATED(matrix%block_list(block_row)%last_used_block_node)) THEN
       IF (matrix%block_list(block_row)%last_used_block_node%block_col .LT. block_col) THEN
           block_node => matrix%block_list(block_row)%last_used_block_node
       ENDIF
    ENDIF

    IF (.NOT.ASSOCIATED(block_node)) THEN
      NULLIFY (new_block_node%next_block_node)
      matrix%block_list(block_row)%first_block_node => new_block_node
    ELSE IF (block_node%block_col > new_block_node%block_col) THEN
      new_block_node%next_block_node =>&
        matrix%block_list(block_row)%first_block_node
      matrix%block_list(block_row)%first_block_node => new_block_node
    ELSE
      DO WHILE (ASSOCIATED(block_node%next_block_node))
        IF (block_node%next_block_node%block_col > new_block_node%block_col) EXIT
        block_node => block_node%next_block_node
      END DO
      new_block_node%next_block_node => block_node%next_block_node
      block_node%next_block_node => new_block_node
    END IF

    ! update
    matrix%block_list(block_row)%last_used_block_node => new_block_node

!   *** Increment block counter ***

    matrix%block_list(block_row)%nblock_node =&
      matrix%block_list(block_row)%nblock_node + 1

  END SUBROUTINE add_real_matrix_block

! *****************************************************************************
!> \brief   Allocate and initialize a new block node in the oce matrix
!> \version 1.0
! *****************************************************************************
  SUBROUTINE  add_1d_block_node(matrix,block_row,block_col,BLOCK,error)

    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row, block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: BLOCK
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_1d_block_node', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cbc, istat, mypcol, myprow, &
                                                ncol, nrow, rbr
    LOGICAL                                  :: failure
    TYPE(real_block_node_type), POINTER      :: block_node, new_block_node

!   perform checks that we are only adding allowed blocks

    failure = .FALSE.
    IF (ASSOCIATED(matrix%distribution_2d)) THEN
       myprow = matrix%distribution_2d%blacs_env%mepos(1)
       mypcol = matrix%distribution_2d%blacs_env%mepos(2)
       rbr    = matrix%distribution_2d%row_distribution(block_col)
       cbc    = matrix%distribution_2d%col_distribution(block_col)
       IF ( .NOT. ( rbr.EQ.myprow .OR. cbc.EQ.mypcol ) ) THEN
          CPPostcondition(.FALSE.,cp_warning_level,routineP,error,failure)
       ENDIF
    ENDIF

!   *** Calculate the block dimensions ***

    nrow = matrix%last_row(block_row) - matrix%first_row(block_row) + 1
    ncol = matrix%last_col(block_col) - matrix%first_col(block_col) + 1

    ALLOCATE (new_block_node,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"Block node of "//matrix%name,0)

!   *** Define the data set of the new block node ***

    new_block_node%block_col = block_col

    ALLOCATE (new_block_node%block(nrow,ncol),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineP,"new_block_node%block",dp_size*nrow*ncol)
    END IF
    new_block_node%block(:,:) = 0.0_dp

    IF (PRESENT(BLOCK)) THEN
      IF (ASSOCIATED(BLOCK)) THEN
        IF ((SIZE(BLOCK,1) == nrow).AND.(SIZE(BLOCK,2) == ncol)) THEN
          new_block_node%block(:,:) = BLOCK(:,:)
        ELSE
          CALL stop_program(routineP,"Incompatible block dimensions")
        END IF
      ELSE
        BLOCK => new_block_node%block
      END IF
    END IF

!   *** Link the new block node to the block list ***

    block_node => matrix%block_list(block_row)%first_block_node

    IF (.NOT.ASSOCIATED(block_node)) THEN
      NULLIFY (new_block_node%next_block_node)
      matrix%block_list(block_row)%first_block_node => new_block_node
    ELSE IF (block_node%block_col > new_block_node%block_col) THEN
      new_block_node%next_block_node =>&
        matrix%block_list(block_row)%first_block_node
      matrix%block_list(block_row)%first_block_node => new_block_node
    ELSE
      DO WHILE (ASSOCIATED(block_node%next_block_node))
        IF (block_node%next_block_node%block_col > new_block_node%block_col) EXIT
        block_node => block_node%next_block_node
      END DO
      new_block_node%next_block_node => block_node%next_block_node
      block_node%next_block_node => new_block_node
    END IF

!   *** Increment block counter ***

    matrix%block_list(block_row)%nblock_node =&
      matrix%block_list(block_row)%nblock_node + 1

  END SUBROUTINE add_1d_block_node

! writes for all local blocks the maximum absolute value of
! the block entries to unit fileunit

! *****************************************************************************
  SUBROUTINE write_blocks_maxabsval(sparse_matrix,fileunit)

    TYPE(real_matrix_type), POINTER          :: sparse_matrix
    INTEGER, INTENT(IN)                      :: fileunit

    CHARACTER(len=*), PARAMETER :: routineN = 'write_blocks_maxabsval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iblock_col, iblock_row
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sparse_block
    TYPE(real_block_node_type), POINTER      :: block_node

    DO iblock_row=1,sparse_matrix%nblock_row

      block_node => first_block_node(sparse_matrix,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            BLOCK=sparse_block)

        WRITE(fileunit,*) iblock_row,iblock_col,MAXVAL(ABS(sparse_block))

        block_node => next_block_node(block_node)

      END DO

    END DO

  END SUBROUTINE write_blocks_maxabsval

! *****************************************************************************
!> \brief   Allocate and initialize a real matrix at the real_matrix_type level.
!>          distribution_2d if present and not NULL, implies that this matrix has
!>          a distribution_2d compatible shape
!> \author  MK
!> \date    16.06.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_real_matrix(matrix,nrow,ncol,nblock_row,nblock_col,&
       first_row,last_row,first_col,last_col,matrix_name,matrix_symmetry,&
       sparsity_id,distribution_2d,error)
    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: nrow, ncol, nblock_row, &
                                                nblock_col
    INTEGER, DIMENSION(nblock_row), &
      INTENT(IN)                             :: first_row, last_row
    INTEGER, DIMENSION(nblock_col), &
      INTENT(IN)                             :: first_col, last_col
    CHARACTER(LEN=*), INTENT(IN)             :: matrix_name, matrix_symmetry
    INTEGER, INTENT(IN)                      :: sparsity_id
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                :: distribution_2d
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_real_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: irow, istat

!   *** Deallocate the old matrix ***

    IF (ASSOCIATED(matrix)) CALL deallocate_matrix(matrix,error=error)
!   *** Allocate a set of block lists ***

    ALLOCATE (matrix,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,matrix%name,0)

    NULLIFY(matrix%distribution_2d)
    IF (PRESENT(distribution_2d)) matrix%distribution_2d=>distribution_2d
    IF (ASSOCIATED(matrix%distribution_2d)) THEN
       CALL distribution_2d_retain(matrix%distribution_2d,error=error)
    END IF

    ALLOCATE (matrix%block_list(nblock_row),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,matrix%name,0)

    last_matrix_id=last_matrix_id+1
    matrix%id_nr=last_matrix_id
    matrix%ref_count=1

    matrix%sparsity_id=sparsity_id

    matrix%name = matrix_name
    matrix%symmetry = matrix_symmetry

    matrix%nblock_row = nblock_row
    matrix%nblock_col = nblock_col

    matrix%nrow = nrow
    matrix%ncol = ncol

    ALLOCATE (matrix%first_row(nblock_row),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineP,"matrix%first_row",int_size*nblock_row)
    END IF
    matrix%first_row(:) = first_row(:)

    ALLOCATE (matrix%last_row(nblock_row),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineP,"matrix%last_row",int_size*nblock_row)
    END IF
    matrix%last_row(:) = last_row(:)

    ALLOCATE (matrix%first_col(nblock_col),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineP,"matrix%first_col",int_size*nblock_col)
    END IF
    matrix%first_col(:) = first_col(:)

    ALLOCATE (matrix%last_col(nblock_col),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineP,"matrix%last_col",int_size*nblock_col)
    END IF
    matrix%last_col(:) = last_col(:)

!   *** Initialize all block lists ***

    DO irow=1,nblock_row
      NULLIFY (matrix%block_list(irow)%first_block_node)
      matrix%block_list(irow)%nblock_node = 0
      NULLIFY (matrix%block_list(irow)%last_used_block_node)
    END DO

  END SUBROUTINE allocate_real_matrix

! *****************************************************************************
!> \brief   Allocate and initialize a real matrix set.
!> \author  MK
!> \date    13.03.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_real_matrix_set(matrix_set,nmatrix,error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_set
    INTEGER, INTENT(IN)                      :: nmatrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_real_matrix_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix, istat

    IF (ASSOCIATED(matrix_set)) CALL deallocate_matrix_set(matrix_set,error=error)

    ALLOCATE (matrix_set(nmatrix),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"matrix_set",0)

    DO imatrix=1,nmatrix
      NULLIFY (matrix_set(imatrix)%matrix)
    END DO

  END SUBROUTINE allocate_real_matrix_set

! *****************************************************************************
!> \brief   Allocate and initialize a real matrix set.
!> \author  MK
!> \date    13.03.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_real_matrix_set_2d(matrix_set,nmatrix,mmatrix,error)
    TYPE(real_matrix_p_type), &
      DIMENSION(:, :), POINTER               :: matrix_set
    INTEGER, INTENT(IN)                      :: nmatrix, mmatrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_real_matrix_set_2d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix, istat, jmatrix

    IF (ASSOCIATED(matrix_set)) CALL deallocate_matrix_set(matrix_set,error=error)

    ALLOCATE (matrix_set(nmatrix,mmatrix),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"matrix_set",0)

    DO jmatrix=1,mmatrix
      DO imatrix=1,nmatrix
        NULLIFY (matrix_set(imatrix,jmatrix)%matrix)
      END DO
    END DO

  END SUBROUTINE allocate_real_matrix_set_2d

! *****************************************************************************
!> \brief  copy to local blocks of a sparse matrix in a replicated full matrix 
!> \par Personal Opinion  
!>          according to me this routine is a bug in parallel JVDV 
!> \author  MK
!> \date    19.06.2001
!> \version 1.0
! *****************************************************************************
  SUBROUTINE copy_local_sm_to_replicated_fm(sparse_matrix,fm)

    TYPE(real_matrix_type), POINTER          :: sparse_matrix
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fm

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'copy_local_sm_to_replicated_fm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_col, &
                                                iblock_row, icol, icol_block, &
                                                irow, irow_block, istat, &
                                                ncol, nrow
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sparse_block
    TYPE(real_block_node_type), POINTER      :: block_node

!   *** Check the association status of the input matrix ***

    CALL timeset(routineN,handle)

    IF (.NOT.ASSOCIATED(sparse_matrix)) THEN
      CALL stop_program(routineN,"The input matrix pointer is not associated")
    END IF

    IF (ASSOCIATED(fm)) DEALLOCATE (fm)

    nrow = sparse_matrix%nrow
    ncol = sparse_matrix%ncol

    ALLOCATE (fm(nrow,ncol),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,"fm",nrow*ncol*dp_size)
    fm(:,:) = 0.0_dp

!   *** Traverse all block nodes of the sparse matrix ***

    DO iblock_row=1,sparse_matrix%nblock_row

      block_node => first_block_node(sparse_matrix,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            BLOCK=sparse_block)

        icol_block = 1

        DO icol=sparse_matrix%first_col(iblock_col),&
                sparse_matrix%last_col(iblock_col)

          irow_block = 1

          DO irow=sparse_matrix%first_row(iblock_row),&
                  sparse_matrix%last_row(iblock_row)

            fm(irow,icol) = sparse_block(irow_block,icol_block)

            irow_block = irow_block + 1

          END DO

          icol_block = icol_block + 1

        END DO

        block_node => next_block_node(block_node)

      END DO

    END DO
 
    CALL timestop(handle)

  END SUBROUTINE copy_local_sm_to_replicated_fm

! multiplies a sparse matrix times ncol local vectors, adding it the the output
! (that are stored as rows ! in v_in leading dimension nblock )
! if the sparse_matrix is symmetric, off diagonal blocks will be used twice
! it is still assumed that the diagonal blocks are symmetric (and hence full)
! Joost VandeVondele july 2002
  ! notice, adding a beta is incompatible with the definition of this function
! *****************************************************************************
  SUBROUTINE sparse_times_local(sparse_matrix,rep_v_in,rep_v_out,ncol,nblock, alpha)
    TYPE(real_matrix_type), POINTER          :: sparse_matrix
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rep_v_in, rep_v_out
    INTEGER, INTENT(IN)                      :: ncol, nblock
    REAL(KIND=dp), OPTIONAL                  :: alpha

    CHARACTER(len=*), PARAMETER :: routineN = 'sparse_times_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_col, &
                                                iblock_row, ithread, nthread, &
                                                sbncol, sbnrow
    LOGICAL                                  :: antisymmetric, owner, &
                                                symmetric
    REAL(KIND=dp)                            :: alpha_l, alpha_l_ji
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sparse_block
    TYPE(real_block_node_type), POINTER      :: block_node

!$  INTEGER :: omp_get_num_threads, omp_get_thread_num

    IF (PRESENT(alpha)) THEN
        alpha_l=alpha
    ELSE
        alpha_l=1.0_dp
    ENDIF

!   *** Check the association status of the input matrix ***
    CALL timeset(routineN,handle)

    IF (ncol.EQ.0) THEN
        CALL timestop(handle)
        RETURN
    ENDIF

    IF (.NOT.ASSOCIATED(sparse_matrix)) THEN
      CALL stop_program(routineP,"The input matrix pointer is not associated")
    END IF

    symmetric=.FALSE.
    antisymmetric=.FALSE.

    SELECT CASE(TRIM(sparse_matrix%symmetry))
    CASE("symmetric")
       symmetric=.TRUE.
       alpha_l_ji=alpha_l  ! we have the same sign for the multiply
    CASE("antisymmetric")
       antisymmetric=.TRUE.
       alpha_l_ji=-alpha_l  ! we have the opposite sign for the multiply
    CASE("none","no symmetry")
       ! nothing
    CASE DEFAULT
       ! possibly just OK (i.e. general matrix) just provide an empty slot
       WRITE(6,*) sparse_matrix%symmetry
       CALL stop_program(routineP,"wrong matrix symmetry specification")
    END SELECT

    IF ((symmetric .OR. antisymmetric) &
           .AND. (sparse_matrix%nrow .NE. sparse_matrix%ncol)) THEN
       CALL stop_program(routineP,"error nonsquare symmetric matrix")
    END IF
    IF (sparse_matrix%ncol .NE. SIZE(rep_v_in,2)) THEN
      CALL stop_program(routineP,"n x n * n x k ?")
    ENDIF

!   *** Traverse all block nodes of the sparse matrix ***
!   no simple parallel do because we write both the iblock_row and iblock_col

!$OMP PARALLEL PRIVATE(ithread,nthread,iblock_row,block_node,iblock_col) &
!$OMP PRIVATE(sparse_block,owner,sbnrow,sbncol)  &
!$OMP SHARED(symmetric,antisymmetric)

    ithread=0
    nthread=1
!$  nthread=omp_get_num_threads()
!$  ithread=omp_get_thread_num()
    DO iblock_row=1,sparse_matrix%nblock_row

      block_node => first_block_node(sparse_matrix,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            BLOCK=sparse_block)

        owner=.TRUE.
        IF (owner) THEN
           sbnrow=sparse_matrix%last_row(iblock_row)- &
                sparse_matrix%first_row(iblock_row)+1
           sbncol=sparse_matrix%last_col(iblock_col)- &
                sparse_matrix%first_col(iblock_col)+1
           IF (sbnrow.NE.0.AND.sbncol.NE.0) THEN
              IF (MOD(iblock_row,nthread).eq.ithread) THEN
                 CALL DGEMM('N','T',ncol,sbnrow,sbncol,alpha_l, &
                            rep_v_in (1,sparse_matrix%first_col(iblock_col)), &
                            nblock, &
                            sparse_block(1,1),sbnrow, &
                            1.0_dp,rep_v_out(1,sparse_matrix%first_row(iblock_row)), &
                            nblock)
              ENDIF

              IF (iblock_col .NE. iblock_row .AND. (symmetric .OR. antisymmetric)) THEN
                 IF (MOD(iblock_col,nthread).eq.ithread) THEN
                    CALL DGEMM('N','N',ncol,sbncol,sbnrow,alpha_l_ji, &
                         rep_v_in (1,sparse_matrix%first_row(iblock_row)), &
                         nblock, &
                         sparse_block(1,1),sbnrow, &
                         1.0_dp,rep_v_out(1,sparse_matrix%first_col(iblock_col)), &
                         nblock)
                 ENDIF
              ENDIF
           ENDIF
        END IF

        block_node => next_block_node(block_node)

      END DO

    END DO
!$OMP END PARALLEL

    CALL timestop(handle)

  END SUBROUTINE sparse_times_local

! intended to compute sparse=sparse+alpha*v*g^T
! for all the local blocks of the sparse matrix
! same comment as sparse_times_local, v,g stored as rows instead of cols
! Joost VandeVondele july 2002
! *****************************************************************************
  SUBROUTINE sparse_plus_loc_loct(sparse_matrix,rep_v,rep_g,ncol,nblock,&
                                  alpha)
    TYPE(real_matrix_type), POINTER          :: sparse_matrix
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rep_v, rep_g
    INTEGER, INTENT(IN)                      :: ncol, nblock
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: alpha

    CHARACTER(len=*), PARAMETER :: routineN = 'sparse_plus_loc_loct', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_col, &
                                                iblock_row, icol, irow, &
                                                ithread, nthread, sbncol, &
                                                sbnrow
    REAL(KIND=dp)                            :: prefactor
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sparse_block
    TYPE(real_block_node_type), POINTER      :: block_node

!$  INTEGER :: omp_get_num_threads, omp_get_thread_num

!   *** Check the association status of the input matrix ***

    IF (ncol.EQ.0) RETURN

    IF (.NOT.ASSOCIATED(sparse_matrix)) THEN
      CALL stop_program(routineP,"The input matrix pointer is not associated")
    END IF

    IF (PRESENT(alpha)) THEN
      prefactor = alpha
    ELSE
      prefactor = 1.0_dp
    END IF

    IF (sparse_matrix%nrow .NE. sparse_matrix%ncol) THEN
      CALL stop_program(routineP,"not allowed")
    ENDIF
    CALL timeset(routineN,handle)

!   *** Traverse all block nodes of the sparse matrix ***
! take into account the triangular schape of the matrix
!$OMP PARALLEL PRIVATE(ithread,nthread,iblock_row,block_node,iblock_col,sparse_block,icol,irow,sbnrow,sbncol) 
    ithread=0
    nthread=1
!$  nthread=omp_get_num_threads()
!$  ithread=omp_get_thread_num()
    DO iblock_row=1,sparse_matrix%nblock_row

      block_node => first_block_node(sparse_matrix,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            BLOCK=sparse_block)

        IF (MOD(iblock_row*sparse_matrix%nblock_row+iblock_col,nthread).eq.ithread) THEN
           icol=sparse_matrix%first_col(iblock_col)
           irow=sparse_matrix%first_row(iblock_row)
           sbnrow=sparse_matrix%last_row(iblock_row)-irow+1
           sbncol=sparse_matrix%last_col(iblock_col)-icol+1

           IF (sbnrow.NE.0 .AND. sbncol.NE.0) THEN
              CALL DGEMM('T','N',sbnrow,sbncol,ncol,prefactor,&
                                 rep_v(1,irow),nblock, &
                                 rep_g(1,icol),nblock, &
                          1.0_dp,sparse_block(1,1),sbnrow)
           ENDIF
        ENDIF

        block_node => next_block_node(block_node)

      END DO

    END DO
!$OMP END PARALLEL

    CALL timestop(handle)

  END SUBROUTINE sparse_plus_loc_loct

! *****************************************************************************
!> \brief   Copy the real matrix source to the existing real matrix target.
!> \author  MK
!> \date    07.07.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE copy_real_matrix(source,TARGET)

    TYPE(real_matrix_type), POINTER          :: source, TARGET

    CHARACTER(len=*), PARAMETER :: routineN = 'copy_real_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: first_col, first_row, handle, iblock_col, iblock_row, icol, &
      irow, jblock_col, jblock_row, jcol, jrow, last_col, last_row, &
      source_first_col, source_first_row, source_last_col, source_last_row, &
      target_first_col, target_first_row, target_last_col, target_last_row
    TYPE(real_block_node_type), POINTER      :: source_block_node, &
                                                target_block_node

!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(source)) THEN
      CALL stop_program(routineP,"The source matrix pointer is not associated")
    END IF

    CALL timeset(routineN,handle)

!   *** Check the association status of the target matrix ***

    IF (ASSOCIATED(TARGET)) THEN

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

      CALL stop_program(routineP,&
                        "The target matrix pointer is not associated. "//&
                        "Use replicate_matrix instead of copy_matrix")

    END IF

    CALL timestop(handle)

  END SUBROUTINE copy_real_matrix

! *****************************************************************************
!> \brief  Deallocate a real matrix at the real_matrix_type level. 
!> \author MK 
!> \date   11.07.2000 
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_real_matrix(matrix,error)

    TYPE(real_matrix_type), POINTER          :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_real_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iblock_row, istat

    IF (.NOT.ASSOCIATED(matrix)) RETURN

    DO iblock_row=1,matrix%nblock_row
      CALL deallocate_matrix_row(matrix,iblock_row)
    END DO

    DEALLOCATE (matrix%block_list,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,matrix%name//" (block_list)")

    DEALLOCATE (matrix%first_row,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,matrix%name//" (first_row)")

    DEALLOCATE (matrix%last_row,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,matrix%name//" (last_row)")

    DEALLOCATE (matrix%first_col,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,matrix%name//" (first_col)")

    DEALLOCATE (matrix%last_col,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,matrix%name//" (last_col)")

    CALL distribution_2d_release(matrix%distribution_2d,error=error)

    DEALLOCATE (matrix,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,matrix%name)

  END SUBROUTINE deallocate_real_matrix

! *****************************************************************************
!> \brief   Deallocate a matrix row.
!> \author  MK
!> \date    30.11.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_real_matrix_row(matrix,block_row)

    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_real_matrix_row', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    TYPE(real_block_node_type), POINTER      :: current_block_node, &
                                                next_block_node

    current_block_node => matrix%block_list(block_row)%first_block_node

    DO WHILE (ASSOCIATED(current_block_node))
      next_block_node => current_block_node%next_block_node
      DEALLOCATE (current_block_node%block,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,matrix%name//" (block)")
      DEALLOCATE (current_block_node,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,matrix%name//" (block_node)")
      current_block_node => next_block_node
    END DO

    NULLIFY (matrix%block_list(block_row)%first_block_node)
    NULLIFY (matrix%block_list(block_row)%last_used_block_node)

    matrix%block_list(block_row)%nblock_node = 0

  END SUBROUTINE deallocate_real_matrix_row

! *****************************************************************************
!> \brief   Deallocate a real matrix set.
!> \author  MK
!> \date    13.03.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_real_matrix_set(matrix_set,error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_real_matrix_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix, istat

    IF (ASSOCIATED(matrix_set)) THEN
      DO imatrix=1,SIZE(matrix_set)
        CALL deallocate_matrix(matrix_set(imatrix)%matrix,error=error)
      END DO
      DEALLOCATE (matrix_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"matrix_set")
    END IF

  END SUBROUTINE deallocate_real_matrix_set

! *****************************************************************************
!> \brief  Deallocate a real matrix set. 
!> \author MK 
!> \date   13.03.2002 
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_real_matrix_set_2d(matrix_set,error)

    TYPE(real_matrix_p_type), &
      DIMENSION(:, :), POINTER               :: matrix_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'deallocate_real_matrix_set_2d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix, istat, jmatrix

    IF (ASSOCIATED(matrix_set)) THEN
      DO jmatrix=1,SIZE(matrix_set,2)
        DO imatrix=1,SIZE(matrix_set,1)
          CALL deallocate_matrix(matrix_set(imatrix,jmatrix)%matrix,error=error)
        END DO
      END DO
      DEALLOCATE (matrix_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"matrix_set")
    END IF

  END SUBROUTINE deallocate_real_matrix_set_2d

! *****************************************************************************
!> \brief   Return a pointer to the requested block node. 
!> \author  MK
!> \date    23.06.2000
!> \version 1.0
! *****************************************************************************
  FUNCTION find_real_block_node(matrix,block_row,block_col) RESULT(block_node)

    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row, block_col
    TYPE(real_block_node_type), POINTER      :: block_node

    CHARACTER(len=*), PARAMETER :: routineN = 'find_real_block_node', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lb, ub

    ub = UBOUND(matrix%block_list,1)
    lb = LBOUND(matrix%block_list,1)

    block_node => matrix%block_list(block_row)%first_block_node

    ! if the last used block node is a short cut, use it
    IF (ASSOCIATED(matrix%block_list(block_row)%last_used_block_node)) THEN
       IF (matrix%block_list(block_row)%last_used_block_node%block_col .LE. block_col) THEN
           block_node => matrix%block_list(block_row)%last_used_block_node
       ENDIF
    ENDIF

    DO WHILE (ASSOCIATED(block_node))
      IF (block_node%block_col == block_col) EXIT
      block_node => block_node%next_block_node
    END DO

    ! update
    matrix%block_list(block_row)%last_used_block_node => block_node

  END FUNCTION find_real_block_node

! *****************************************************************************
!> \brief  Return a pointer to the first block node of a block list. 
!> \author MK 
!> \date   23.06.2000 
!> \version 1.0
! *****************************************************************************
  FUNCTION first_real_block_node(matrix,block_row) RESULT(first_block_node)

    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row
    TYPE(real_block_node_type), POINTER      :: first_block_node

    first_block_node => matrix%block_list(block_row)%first_block_node

  END FUNCTION first_real_block_node

! *****************************************************************************
!> \brief   Return the requested matrix information.
!> \author  MK
!> \date    10.07.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_matrix_info(matrix,matrix_name,matrix_symmetry,&
                             nblock_row,nblock_col,nrow,ncol,&
                             first_row,last_row,first_col,last_col,&
                             nblock_allocated,nelement_allocated,sparsity_id)

    TYPE(real_matrix_type), POINTER          :: matrix
    CHARACTER(LEN=80), INTENT(OUT), OPTIONAL :: matrix_name
    CHARACTER(LEN=40), INTENT(OUT), OPTIONAL :: matrix_symmetry
    INTEGER, INTENT(OUT), OPTIONAL           :: nblock_row, nblock_col, nrow, &
                                                ncol
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: first_row, last_row, &
                                                first_col, last_col
    INTEGER, INTENT(OUT), OPTIONAL           :: nblock_allocated, &
                                                nelement_allocated, &
                                                sparsity_id

    INTEGER                                  :: iblock_row
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: BLOCK
    TYPE(real_block_node_type), POINTER      :: block_node

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
    IF (PRESENT(sparsity_id)) sparsity_id = matrix%sparsity_id

    IF (PRESENT(nblock_allocated)) THEN
      nblock_allocated = 0
      DO iblock_row=1,matrix%nblock_row
        nblock_allocated = nblock_allocated +&
                           matrix%block_list(iblock_row)%nblock_node
      END DO
    END IF

    IF (PRESENT(nelement_allocated)) THEN
      nelement_allocated = 0
      DO iblock_row=1,matrix%nblock_row
        block_node => first_block_node(matrix,iblock_row)
        DO WHILE (ASSOCIATED(block_node))
          CALL get_block_node(block_node=block_node,BLOCK=BLOCK)
          nelement_allocated = nelement_allocated + SIZE(BLOCK)
          block_node => next_block_node(block_node)
        END DO
      END DO
    END IF

  END SUBROUTINE get_matrix_info

! *****************************************************************************
!> \brief  Get the the column and/or matrix-block for a given sparse matrix block node. 
!> \author MK 
!> \date   28.05.2000 
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_real_block_node(block_node,block_col,BLOCK)

    TYPE(real_block_node_type), INTENT(IN)   :: block_node
    INTEGER, INTENT(OUT), OPTIONAL           :: block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: BLOCK

    IF (PRESENT(block_col)) block_col = block_node%block_col
    IF (PRESENT(BLOCK)) BLOCK => block_node%block

  END SUBROUTINE get_real_block_node

! *****************************************************************************
!> \brief  Get a specific (block_row,block_col)-block from the real matrix 
!> \author  MK
!> \date   28.05.2000 
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_real_matrix_block(matrix,block_row,block_col,&
                                   first_row,last_row,first_col,last_col,&
                                   block_node,BLOCK)
    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row, block_col
    INTEGER, INTENT(OUT), OPTIONAL           :: first_row, last_row, &
                                                first_col, last_col
    TYPE(real_block_node_type), OPTIONAL, &
      POINTER                                :: block_node
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: BLOCK

    TYPE(real_block_node_type), POINTER      :: current_block_node

    IF (PRESENT(first_row)) first_row = matrix%first_row(block_row)
    IF (PRESENT(last_row)) last_row = matrix%last_row(block_row)

    IF (PRESENT(first_col)) first_col = matrix%first_col(block_col)
    IF (PRESENT(last_col)) last_col = matrix%last_col(block_col)

    current_block_node => find_real_block_node(matrix,block_row,block_col)

    IF (ASSOCIATED(current_block_node)) THEN
      IF (PRESENT(block_node)) block_node => current_block_node
      IF (PRESENT(BLOCK)) BLOCK => current_block_node%block
    ELSE
      IF (PRESENT(block_node)) NULLIFY (block_node)
      IF (PRESENT(BLOCK)) NULLIFY (BLOCK)
    END IF

  END SUBROUTINE get_real_matrix_block

! *****************************************************************************
!> \brief   Set the diagonal elements of matrix to the values of vector.
!> \author  MK
!> \date    06.11.2003
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_real_matrix_diagonal(matrix,vector)

    TYPE(real_matrix_type), POINTER          :: matrix
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: vector

    CHARACTER(len=*), PARAMETER :: routineN = 'get_real_matrix_diagonal', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_row, i, iblock_col, &
                                                iblock_row, irow
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block
    TYPE(real_block_node_type), POINTER      :: block_node

!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(matrix)) THEN
      CALL stop_program(routineP,"The matrix pointer is not associated")
    END IF

    vector(:) = 0.0_dp

    DO iblock_row=1,matrix%nblock_row

      first_row = matrix%first_row(iblock_row)

      block_node => first_block_node(matrix,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            block=block)

        IF (iblock_row == iblock_col) THEN
          DO i=1,SIZE(block,1)
            irow = first_row + i - 1
            vector(irow) = block(i,i)
          END DO
        END IF

        block_node => next_block_node(block_node)

      END DO

    END DO

  END SUBROUTINE get_real_matrix_diagonal

! *****************************************************************************
!> \brief   Return a pointer to the next block node of a block list.
!> \author  MK
!> \date    28.06.2000
!> \version 1.0
! *****************************************************************************
  FUNCTION next_real_block_node(block_node) RESULT(next_block_node)

    TYPE(real_block_node_type), POINTER      :: block_node, next_block_node

    IF (ASSOCIATED(block_node)) THEN
       next_block_node => block_node%next_block_node
    ELSE
       CALL stop_program("sparse_matrix_types:next_real_block_node",&
            "The block_node pointer a is not associated")
       NULLIFY(next_block_node)
    END IF

  END FUNCTION next_real_block_node

! *****************************************************************************
!> \brief   Set block node data set.
!> \author  MK
!> \date    27.07.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE put_real_block_node(block_node,matrix,block_row,block_col,BLOCK,add)
    TYPE(real_block_node_type), POINTER      :: block_node
    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row
    INTEGER, INTENT(IN), OPTIONAL            :: block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: BLOCK
    LOGICAL, INTENT(IN), OPTIONAL            :: add

    CHARACTER(len=*), PARAMETER :: routineN = 'put_real_block_node', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ncol, nrow
    LOGICAL                                  :: add_block

    IF (PRESENT(add)) THEN
      add_block = add
    ELSE
      add_block = .FALSE.
    END IF

    IF (PRESENT(block_col)) block_node%block_col = block_col

    IF (PRESENT(BLOCK)) THEN

      nrow = matrix%last_row(block_row) - matrix%first_row(block_row) + 1
      ncol = matrix%last_col(block_col) - matrix%first_col(block_col) + 1

      IF (PRESENT(BLOCK)) THEN
        IF (ASSOCIATED(BLOCK)) THEN
          IF ((SIZE(BLOCK,1) == nrow).AND.(SIZE(BLOCK,2) == ncol)) THEN
            IF ( add_block ) THEN
              block_node%block(:,:) = block_node%block(:,:) + BLOCK(:,:)
            ELSE
              block_node%block(:,:) = BLOCK(:,:)
            END IF
          ELSE
            CALL stop_program(routineP,"Incompatible block dimensions")
          END IF
        ELSE
          CALL stop_program(routineP,"Block is not associated")
        END IF
      END IF

    END IF

  END SUBROUTINE put_real_block_node

! *****************************************************************************
!> \brief Allocate and initialize a new block node or update an existing block node.  
!> \author  MK
!> \date    07.07.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE put_real_matrix_block(matrix,block_row,block_col,BLOCK,add,error)

    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row, block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: BLOCK
    LOGICAL, INTENT(IN), OPTIONAL            :: add
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'put_real_matrix_block', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ncol, nrow
    LOGICAL                                  :: add_block
    TYPE(real_block_node_type), POINTER      :: block_node

    IF ( PRESENT (add) ) THEN
      add_block = add
    ELSE
      add_block = .FALSE.
    ENDIF

!   *** Check, if the requested block node already exists ***

    block_node => find_real_block_node(matrix,block_row,block_col)

    IF (ASSOCIATED(block_node)) THEN

!     *** Calculate the block dimensions ***

      nrow = matrix%last_row(block_row) - matrix%first_row(block_row) + 1
      ncol = matrix%last_col(block_col) - matrix%first_col(block_col) + 1

!     *** Update the data set of an existing block node ***

      IF (PRESENT(BLOCK)) THEN
        IF (ASSOCIATED(BLOCK)) THEN
          IF ((SIZE(BLOCK,1) == nrow).AND.(SIZE(BLOCK,2) == ncol)) THEN
            IF ( add_block ) THEN
              block_node%block(:,:) = block_node%block(:,:) + BLOCK(:,:)
            ELSE
              block_node%block(:,:) = BLOCK(:,:)
            END IF
          ELSE
            CALL stop_program(routineP,"Incompatible block dimensions")
          END IF
        ELSE
          CALL stop_program(routineP,"Matrix block is not associated")
        END IF
      ELSE
        CALL stop_program(routineP,&
                          "Block update is impossible, because no matrix "//&
                          "block was specified")
      END IF

    ELSE

!     *** Add a new block node ***

      IF (PRESENT(BLOCK)) THEN
        CALL add_block_node(matrix,block_row,block_col,BLOCK,error=error)
      ELSE
        CALL add_block_node(matrix,block_row,block_col,error=error)
      END IF

    END IF

  END SUBROUTINE put_real_matrix_block

! *****************************************************************************
!> \brief  Replicate the existing matrix source. The replicated matrix is
!>         target with the name target_name.
!>         by default copies the data blocks,
!>         if optional allocate_blocks=false then no blocks are allocated  
!> \author MK 
!> \date    17.11.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE replicate_real_matrix(source,TARGET,target_name,allocate_blocks,error)

    TYPE(real_matrix_type), POINTER          :: source, TARGET
    CHARACTER(LEN=*), INTENT(IN)             :: target_name
    LOGICAL, INTENT(IN), OPTIONAL            :: allocate_blocks
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'replicate_real_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iblock_col, iblock_row
    LOGICAL                                  :: make_data_blocks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: BLOCK
    TYPE(real_block_node_type), POINTER      :: block_node

!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(source)) THEN
      CALL stop_program(routineP,"The source matrix pointer is not associated")
    END IF

    IF (ASSOCIATED(TARGET)) CALL deallocate_matrix(TARGET,error=error)

!   *** Allocate a new matrix structure ***

    CALL allocate_matrix(matrix=TARGET,&
                         nrow=source%nrow,&
                         ncol=source%ncol,&
                         nblock_row=source%nblock_row,&
                         nblock_col=source%nblock_col,&
                         first_row=source%first_row(:),&
                         last_row=source%last_row(:),&
                         first_col=source%first_col(:),&
                         last_col=source%last_col(:),&
                         matrix_name=target_name,&
                         sparsity_id=source%sparsity_id, &
                         distribution_2d=source%distribution_2d,&
                         matrix_symmetry=source%symmetry,error=error)

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
                              BLOCK=BLOCK)

          CALL add_block_node(matrix=TARGET,&
                              block_row=iblock_row,&
                              block_col=iblock_col,&
                              BLOCK=BLOCK,error=error)

          block_node => next_block_node(block_node)

        END DO

      END DO

    ENDIF

  END SUBROUTINE replicate_real_matrix

! *****************************************************************************
!> \brief  Replicate the matrix structure of the existing matrix source. The
!>          new matrix target with the name target_name has the same
!>          structure. 
!> \author  MK
!> \date    17.11.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE replicate_real_matrix_structure(source,TARGET,target_name,&
                                             target_symmetry,error)
    TYPE(real_matrix_type), POINTER          :: source, TARGET
    CHARACTER(LEN=*), INTENT(IN)             :: target_name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: target_symmetry
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'replicate_real_matrix_structure', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=40)                        :: matrix_symmetry
    INTEGER                                  :: handle, iblock_col, iblock_row
    TYPE(real_block_node_type), POINTER      :: block_node

!   *** Check the association status of the source matrix ***

    CALL timeset(routineN,handle)

    IF (.NOT.ASSOCIATED(source)) THEN
      CALL stop_program(routineP,"The source matrix pointer is not associated")
    END IF

    IF (PRESENT(target_symmetry)) THEN
      matrix_symmetry = target_symmetry
    ELSE
      matrix_symmetry = source%symmetry
    END IF

!   *** Allocate a new matrix structure ***
    CALL allocate_matrix(matrix=TARGET,&
                         nrow=source%nrow,&
                         ncol=source%ncol,&
                         nblock_row=source%nblock_row,&
                         nblock_col=source%nblock_col,&
                         first_row=source%first_row(:),&
                         last_row=source%last_row(:),&
                         first_col=source%first_col(:),&
                         last_col=source%last_col(:),&
                         matrix_name=target_name,&
                         sparsity_id=source%sparsity_id,&
                         distribution_2d=source%distribution_2d,&
                         matrix_symmetry=matrix_symmetry,error=error)

!   *** Initialize all block nodes ***

    DO iblock_row=1,source%nblock_row

      block_node => first_block_node(source,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col)

        IF(matrix_symmetry=="none" .OR.matrix_symmetry=="non symmetric")THEN
          CALL add_1d_block_node(matrix=TARGET,&
                            block_row=iblock_row,&
                            block_col=iblock_col,error=error)
        ELSE
          CALL add_block_node(matrix=TARGET,&
                            block_row=iblock_row,&
                            block_col=iblock_col,error=error)
        END IF

        block_node => next_block_node(block_node)

      END DO

    END DO

    CALL timestop(handle)

  END SUBROUTINE replicate_real_matrix_structure

! *****************************************************************************
!> \brief  Multiply the sparse matrix a with alpha.
!>            a <- alpha*a 
!> \author MK 
!> \date    17.07.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE scale_real_matrix(a,alpha)

    TYPE(real_matrix_type), POINTER          :: a
    REAL(KIND=dp), INTENT(IN)                :: alpha

    CHARACTER(len=*), PARAMETER :: routineN = 'scale_real_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iblock_row
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: BLOCK
    TYPE(real_block_node_type), POINTER      :: block_node

!   *** Check the association status of the matrices ***

    IF (.NOT.ASSOCIATED(a)) THEN
       CALL stop_program(routineP,"The matrix pointer a is not associated")
    END IF

    IF (alpha == 1.0_dp) RETURN

    IF (alpha == 0.0_dp) THEN

      CALL set_matrix(a,0.0_dp)

    ELSE

      DO iblock_row=1,a%nblock_row

        block_node => first_block_node(a,iblock_row)

        DO WHILE (ASSOCIATED(block_node))

          CALL get_block_node(block_node=block_node,&
                              BLOCK=BLOCK)

          BLOCK(:,:) = alpha*BLOCK(:,:)

          block_node => block_node%next_block_node

        END DO

      END DO

    END IF

  END SUBROUTINE scale_real_matrix

! *****************************************************************************
  FUNCTION checksum_real_matrix(a,para_env)
    TYPE(real_matrix_type), POINTER          :: a
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp)                            :: checksum_real_matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'checksum_real_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iblock_row, n
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: BLOCK
    REAL(KIND=dp), EXTERNAL                  :: DDOT
    TYPE(real_block_node_type), POINTER      :: block_node

!   *** Check the association status of the matrices ***

    checksum_real_matrix = 0.0_dp
    IF (.NOT.ASSOCIATED(a)) THEN
       CALL stop_program(routineP,"The matrix pointer a is not associated")
    END IF
    IF (.NOT.ASSOCIATED(para_env)) THEN
       CALL stop_program(routineP,"The para_env pointer a is not associated")
    END IF
    DO iblock_row=1,a%nblock_row
       block_node => first_block_node(a,iblock_row)
       DO WHILE (ASSOCIATED(block_node))
          CALL get_block_node(block_node=block_node,&
                              BLOCK=BLOCK)
          n = SIZE(BLOCK,1)*SIZE(BLOCK,2)
          checksum_real_matrix = checksum_real_matrix + DDOT(n,BLOCK(1,1),1,BLOCK(1,1),1)
          block_node => block_node%next_block_node
       END DO
    END DO
    CALL mp_sum(checksum_real_matrix,para_env%group)
  END FUNCTION checksum_real_matrix

! *****************************************************************************
!> \brief   Set all elements of matrix to value.
!> \author  MK
!> \date    10.04.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_real_matrix(matrix,value)
    TYPE(real_matrix_type), POINTER          :: matrix
    REAL(KIND=dp), INTENT(IN)                :: value

    CHARACTER(len=*), PARAMETER :: routineN = 'set_real_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iblock_row
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: BLOCK
    TYPE(real_block_node_type), POINTER      :: block_node

!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(matrix)) THEN
      CALL stop_program(routineP,"The matrix pointer is not associated")
    END IF

!   *** Set matrix elements to value ***

    DO iblock_row=1,matrix%nblock_row

      block_node => first_block_node(matrix,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            BLOCK=BLOCK)

        BLOCK(:,:) = value

        block_node => block_node%next_block_node

      END DO

    END DO

  END SUBROUTINE set_real_matrix

! *****************************************************************************
!> \brief  Set the diagonal elements of matrix to the values of vector. 
!> \author MK 
!> \date   06.11.2003 
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_real_matrix_diagonal(matrix,vector)

    TYPE(real_matrix_type), POINTER          :: matrix
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: vector

    CHARACTER(len=*), PARAMETER :: routineN = 'set_real_matrix_diagonal', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_row, i, iblock_col, &
                                                iblock_row, irow
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block
    TYPE(real_block_node_type), POINTER      :: block_node

!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(matrix)) THEN
      CALL stop_program(routineP,"The matrix pointer is not associated")
    END IF

    DO iblock_row=1,matrix%nblock_row

      first_row = matrix%first_row(iblock_row)

      block_node => first_block_node(matrix,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            block=block)

        IF (iblock_row == iblock_col) THEN
          DO i=1,SIZE(block,1)
            irow = first_row + i - 1
            block(i,i) = vector(irow)
          END DO
        END IF

        block_node => next_block_node(block_node)

      END DO

    END DO

  END SUBROUTINE set_real_matrix_diagonal

! *****************************************************************************
!> \brief Transfers source matrix to target matrix (copies the values)
!>    filtering the result through the structure (sparsity) of target matrix
!>    zeros blocks as necessary
!> \param source_matrix the matrix to copy
!> \param target_matrix the matrix that will contain the copy
!> \param error variable to control error logging, stopping,...
!>             see module cp_error_handling
!> \note
!>      assumes that both matrixes have the same distribution
!>      XXXXX very dangerous (i.e. wrong?) in parallell, if the sparsities are not the same XXXXXXX
!> \par History
!>      - 10.04.2002 Creation [MK]
!>      - 4.2003 extended for matrixes with different structures [fawzi]
!> \author Matthias Krack
! *****************************************************************************
  SUBROUTINE transfer_real_matrix(source_matrix,target_matrix ,error)

    TYPE(real_matrix_type), POINTER          :: source_matrix, target_matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'transfer_real_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_col_source, &
                                                iblock_col_target, iblock_row
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: source_block, target_block
    TYPE(real_block_node_type), POINTER      :: block_node_source, &
                                                block_node_target

!   *** Check the association status of the source matrix ***

    CALL timeset(routineN,handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(source_matrix),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(target_matrix),cp_failure_level,routineP,error,failure)

!   *** Copy matrix elements ***
    !
    ! we loop over the blocks, making use of the fact that they are ordered
    ! since the sparsity might explicitly be different, this requires checking that
    ! for a given target block, the source block is actually corresponding
    !
    ! for all rows, go over the linked list of blocks
    DO iblock_row=1,source_matrix%nblock_row

      block_node_target => first_block_node(target_matrix,iblock_row)
      block_node_source => first_block_node(source_matrix,iblock_row)

      DO WHILE (ASSOCIATED(block_node_target))

        CALL get_block_node(block_node=block_node_target,&
                            block_col=iblock_col_target,&
                            BLOCK=target_block)

        IF (ASSOCIATED(target_block)) THEN

           NULLIFY(source_block)

           DO WHILE(ASSOCIATED(block_node_source))
              CALL get_block_node(block_node=block_node_source,&
                                  block_col=iblock_col_source,&
                                  BLOCK=source_block)
              IF (iblock_col_source == iblock_col_target) EXIT
              NULLIFY(source_block)
              IF (iblock_col_source >  iblock_col_target) EXIT
              block_node_source => block_node_source%next_block_node
           ENDDO

           IF (ASSOCIATED(source_block)) THEN
              IF (SIZE(target_block,1)*SIZE(target_block,2).NE.0) THEN
                 CALL dcopy(SIZE(target_block,1)*SIZE(target_block,2),&
                      source_block(1,1),1,target_block(1,1),1)
              ENDIF
           ELSE
              IF (SIZE(target_block,1)*SIZE(target_block,2).NE.0) THEN
                 CALL dcopy(SIZE(target_block,1)*SIZE(target_block,2),&
                      0.0_dp,0,target_block(1,1),1)
              ENDIF
           END IF
        END IF

        block_node_target => block_node_target%next_block_node

      END DO

    END DO
    CALL timestop(handle)

  END SUBROUTINE transfer_real_matrix

! *****************************************************************************
!> \brief removes a block from the matrix
!> \param matrix the matrix you remove the block from
!> \param block_row the row of the block to remove
!> \param block_col the column of the block to remove
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      7.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE remove_real_matrix_block(matrix,block_row, block_col, error)
    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(in)                      :: block_row, block_col
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'remove_real_matrix_block', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(real_block_node_type), POINTER      :: block_node, old_block

    NULLIFY(block_node,old_block)
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       block_node => matrix%block_list(block_row)%first_block_node
       ! better safe than sorry
       NULLIFY(matrix%block_list(block_row)%last_used_block_node)
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

! *****************************************************************************
  SUBROUTINE get_real_matrix_occupation(matrix, block_occupation, memory_occupation, fill_symmetric)

    TYPE(real_matrix_type), POINTER          :: matrix
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: block_occupation, &
                                                memory_occupation
    LOGICAL, OPTIONAL                        :: fill_symmetric

    INTEGER :: block_counter, block_memory, diagonal_memory, i, j, &
      max_blocks, max_memory, memory_counter, n_diagonal_blocks
    TYPE(real_block_node_type), POINTER      :: node

! should symmetric matrices be fully counted or just one triangle?
! code

    IF (SIZE(matrix%block_list) /= SIZE(matrix%first_col)) THEN
       CALL stop_program("get_real_matrix_occupation", "precondition")
    END IF

    max_blocks = matrix%nblock_row * matrix%nblock_col
    max_memory = matrix%ncol * matrix%nrow
    block_counter = 0
    memory_counter = 0
    n_diagonal_blocks = 0
    diagonal_memory = 0

    DO i=1, matrix%nblock_row
       node => first_block_node(matrix, i)
       DO WHILE (ASSOCIATED(NODE))
          block_memory = SIZE(node%block,1) * SIZE(node%block,2)
          j = node%block_col
          IF (i == j) THEN
             n_diagonal_blocks = n_diagonal_blocks + 1
             diagonal_memory = diagonal_memory + block_memory
          END IF
          block_counter = block_counter + 1
          memory_counter = memory_counter + block_memory
          node => next_block_node(node)
       END DO
    END DO

    IF (TRIM(matrix%symmetry) == "symmetric") THEN
       IF (.NOT.PRESENT(fill_symmetric) .OR. fill_symmetric) THEN
          block_counter = 2*block_counter - n_diagonal_blocks
          memory_counter = 2*memory_counter - diagonal_memory
       END IF
    END IF
    IF (PRESENT(block_occupation)) block_occupation = REAL(block_counter, dp ) / REAL(max_blocks, dp )
    IF (PRESENT(memory_occupation)) memory_occupation = REAL(memory_counter, dp ) / REAL(max_memory, dp )

  END SUBROUTINE get_real_matrix_occupation

! *****************************************************************************
!> \brief returns the id of the given matrix
!> \param matrix the matrix you want info about
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
FUNCTION cp_sm_get_id_nr(matrix,error) RESULT(res)
    TYPE(real_matrix_type), POINTER          :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_get_id_nr', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
!  CPPrecondition(associated(matrix),cp_failure_level,routineP,error,failure)
  CPPrecondition(matrix%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     res=matrix%id_nr
  ELSE
    res=0
  END IF
END FUNCTION cp_sm_get_id_nr

! *****************************************************************************
!> \brief sets various attributes of the given matrix
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE cp_sm_set(matrix,name,error)
    TYPE(real_matrix_type), POINTER          :: matrix
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
  CPPrecondition(matrix%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(name)) matrix%name=name
  END IF
END SUBROUTINE cp_sm_set

! *****************************************************************************
!> \brief gets a list of indices of all allocated blocks (locally) of a
!>      given matrix
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      none
!> \author JGH (04-07-2003)
! *****************************************************************************
SUBROUTINE cp_sm_get_block_list ( matrix, list, error )
    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, DIMENSION(:, :), POINTER        :: list
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_get_block_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ib, istat, j, num_blocks
    LOGICAL                                  :: failure
    TYPE(real_block_node_type), POINTER      :: node

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
  CPPrecondition(.NOT.ASSOCIATED(list),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN

    CALL get_matrix_info ( matrix, nblock_allocated=num_blocks )
    ALLOCATE (list(2,num_blocks),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"list",2*num_blocks)

    ib = 0
    DO i=1, matrix%nblock_row
       node => first_block_node(matrix, i)
       DO WHILE (ASSOCIATED(node))
          IF ( ASSOCIATED(node%block) ) THEN
            j = node%block_col
            ib = ib + 1
            IF ( ib > num_blocks ) CALL stop_program(routineP,"too many blocks")
            list(1,ib) = i
            list(2,ib) = j
          END IF
          node => next_block_node(node)
       END DO
    END DO

  END IF

END SUBROUTINE cp_sm_get_block_list

! *****************************************************************************
!> \brief returns the maximum of the abs of the matrix entries
!> \param matrix the matrix you are interested in
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE cp_sm_maxabsval(matrix,maxabsval,para_env,error)
    TYPE(real_matrix_type), POINTER          :: matrix
    REAL(KIND=dp), INTENT(out)               :: maxabsval
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_maxabsval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iblock_row
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix_block
    TYPE(real_block_node_type), POINTER      :: block_node

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
  CPPrecondition(matrix%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)

  maxabsval=0.0_dp

  DO iblock_row=1,matrix%nblock_row

     block_node => first_block_node(matrix,iblock_row)

     DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            BLOCK=matrix_block)

        maxabsval=MAX(maxabsval,MAXVAL(ABS(matrix_block)))

        block_node => block_node%next_block_node

      END DO

    END DO

    CALL mp_max(maxabsval,para_env%group)

  END SUBROUTINE cp_sm_maxabsval

! *****************************************************************************
!> \brief returns the trace of the product of two sparse matrices
!> \param matrix the matrix you are interested in
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      uses Kahan's sum
!> \par History
!>      10.2004 created [jgh]
!> \author jgh
! *****************************************************************************
SUBROUTINE cp_sm_sm_trace_scalar(matrix1,matrix2,trace,para_env,error)
    TYPE(real_matrix_type), POINTER          :: matrix1, matrix2
    REAL(KIND=dp), INTENT(out)               :: trace
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_sm_trace_scalar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: block_col1, block_col2, &
                                                iblock_row, nblock_row
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: c, ks, sfactor, t, &
                                                trace_block, y
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix_block1, matrix_block2
    TYPE(real_block_node_type), POINTER      :: block_node1, block_node2

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(matrix1),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix2),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix1%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix2%ref_count>0,cp_failure_level,routineP,error,failure)
    IF ( PRESENT ( para_env ) ) THEN
      CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
      CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF

    trace=0.0_dp
    ks = 0.0_dp ;  t = 0.0_dp ; y = 0.0_dp ; c = 0.0_dp

    CALL get_matrix_info(matrix=matrix1,nblock_row=nblock_row)
    ! we assume symmetric matrices
    sfactor = 2._dp

    DO iblock_row=1,nblock_row

     block_node1 => first_block_node(matrix1,iblock_row)
     block_node2 => first_block_node(matrix2,iblock_row)

     DO WHILE (ASSOCIATED(block_node1) .AND. ASSOCIATED(block_node2))

        CALL get_block_node(block_node=block_node1,&
                            block_col=block_col1,BLOCK=matrix_block1)
        CALL get_block_node(block_node=block_node2,&
                            block_col=block_col2,BLOCK=matrix_block2)
        CPPostcondition(block_col1==block_col2,cp_warning_level,routineP,error,failure)

        trace_block = accurate_sum ( matrix_block1*matrix_block2 )

        IF ( iblock_row == block_col1 ) THEN
          trace_block = trace_block
        ELSE
          trace_block = sfactor*trace_block
        END IF

        ! use Kahan's sum to reduce roundoff
        y  = trace_block - c
        t  = ks + y
        c  = (t - ks) - y
        ks = t

        block_node1 => block_node1%next_block_node
        block_node2 => block_node2%next_block_node

     END DO

    END DO

    trace = ks

! If not present, assumes that you have a replicated matrix structure
    IF ( PRESENT ( para_env ) ) &
    CALL mp_sum(trace,para_env%group)

END SUBROUTINE cp_sm_sm_trace_scalar

! *****************************************************************************
!> \brief Calculates blockwise the trace of a matrix product of two sparse matrices,
!>      this is necessary if the forceone application is the calculation of the
!>      forces via F=Tr[p*(dH/dR)]
!> \note
!>      - matrix1: 1st multiplier
!>      - matrix2: 2nd multiplier
!>      - trace_vec: The vector in which the blockwise traces are written
!>      - para_env: The paralle environment
!>      - error: variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2006 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
SUBROUTINE cp_sm_sm_trace_vector(matrix1, matrix2, trace_vec, para_env, error)
    TYPE(real_matrix_type), POINTER          :: matrix1, matrix2
    REAL(KIND=dp), DIMENSION(:), POINTER     :: trace_vec
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_sm_trace_vector', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: block_col1, block_col2, &
                                                iblock_row, nblock_row
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: trace_block
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix_block1, matrix_block2
    TYPE(real_block_node_type), POINTER      :: block_node1, block_node2

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(matrix1), cp_failure_level, routineP, error, failure)
  CPPrecondition(ASSOCIATED(matrix2), cp_failure_level, routineP, error, failure)
  CPPrecondition(matrix1%ref_count>0, cp_failure_level, routineP, error, failure)
  CPPrecondition(matrix2%ref_count>0, cp_failure_level, routineP, error, failure)
  CPPrecondition(ASSOCIATED(para_env), cp_failure_level, routineP, error, failure)
  CPPrecondition(para_env%ref_count>0, cp_failure_level, routineP, error, failure)

  trace_vec(:) = 0.0_dp

  CALL get_matrix_info(matrix=matrix1, nblock_row=nblock_row)

  DO iblock_row = 1,nblock_row

     block_node1 => first_block_node(matrix1, iblock_row)
     block_node2 => first_block_node(matrix2, iblock_row)

     DO WHILE (ASSOCIATED(block_node1) .AND. ASSOCIATED(block_node2))

        CALL get_block_node(block_node=block_node1, &
                            block_col=block_col1, BLOCK=matrix_block1)
        CALL get_block_node(block_node=block_node2, &
                            block_col=block_col2, BLOCK=matrix_block2)
        CPPostcondition(block_col1==block_col2, cp_warning_level, routineP, error, failure)

        trace_block = accurate_sum(matrix_block1*matrix_block2)

        IF (iblock_row == block_col1) THEN
          trace_vec(iblock_row) = trace_vec(iblock_row) + trace_block
        ELSE
          ! we assume symmetric matrices
          trace_vec(iblock_row) = trace_vec(iblock_row) + trace_block
          trace_vec(nblock_row-iblock_row+1) = trace_vec(nblock_row-iblock_row+1) &
                                             + trace_block
        END IF

        block_node1 => block_node1%next_block_node
        block_node2 => block_node2%next_block_node

     END DO

  END DO

  CALL mp_sum(trace_vec, para_env%group)

END SUBROUTINE cp_sm_sm_trace_vector

! *****************************************************************************
!> \brief  Create a diagonal block matrix from a real sparse matrix
!>         same structure. 
!> \version 1.0
! *****************************************************************************
SUBROUTINE create_dbm_from_sm(dbm,sm,matrix_name,copy_in,para_env,error)

    TYPE(real_matrix_type), POINTER          :: dbm, sm
    CHARACTER(LEN=*), INTENT(IN)             :: matrix_name
    LOGICAL, INTENT(IN)                      :: copy_in
    TYPE(cp_para_env_type), INTENT(IN), &
      OPTIONAL                               :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_dbm_from_sm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_col, iblock_row
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block
    TYPE(real_block_node_type), POINTER      :: block_node

    CALL timeset(routineN,handle)

    CALL allocate_matrix(matrix=dbm,&
                         nblock_row=sm%nblock_row,&
                         nblock_col=sm%nblock_col,&
                         nrow=sm%nrow,&
                         ncol=sm%ncol,&
                         first_row=sm%first_row,&
                         last_row=sm%last_row,&
                         first_col=sm%first_col,&
                         last_col=sm%last_col,&
                         matrix_name=matrix_name,&
                         matrix_symmetry="symmetric",&
                         sparsity_id=1111,error=error)

    DO iblock_row = 1, sm%nblock_row
       CALL add_block_node(matrix=dbm,&
                           block_row=iblock_row,&
                           block_col=iblock_row,error=error)
    END DO
    IF ( copy_in ) THEN

      DO iblock_row=1,sm%nblock_row
        block_node => first_block_node(sm,iblock_row)
        DO WHILE (ASSOCIATED(block_node))
          NULLIFY(block)
          CALL get_block_node(block_node=block_node,&
                              block_col=iblock_col,&
                              block=block)
          IF (iblock_row == iblock_col) THEN
            CALL put_block_node(dbm,iblock_row,iblock_row,block,error=error)
          END IF
          block_node => next_block_node(block_node)
        END DO
      END DO
      IF ( PRESENT(para_env) ) THEN
        DO iblock_row=1,dbm%nblock_row
          NULLIFY(block)
          CALL get_block_node(matrix=dbm,&
                              block_row=iblock_row,block_col=iblock_row,&
                              BLOCK=block)
          CALL mp_sum(block,para_env%group)
        END DO
      END IF

    END IF

    CALL timestop(handle)

END SUBROUTINE create_dbm_from_sm
! *****************************************************************************
!> \brief  mp_sum the diagonal blocks 
!> \version 1.0
! *****************************************************************************
SUBROUTINE mpsum_dbm(dbm,para_env,error)
    TYPE(real_matrix_type), POINTER          :: dbm
    TYPE(cp_para_env_type), INTENT(IN)       :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mpsum_dbm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_row
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block

    CALL timeset(routineN,handle)
    DO iblock_row=1,dbm%nblock_row
      NULLIFY(block)
      CALL get_block_node(matrix=dbm,&
                          block_row=iblock_row,block_col=iblock_row,&
                          BLOCK=block)
      CALL mp_sum(block,para_env%group)
    END DO
    CALL timestop(handle)
END SUBROUTINE mpsum_dbm

! *****************************************************************************
!> \brief   Deallocate a diagonal block matrix
!> \version 1.0
! *****************************************************************************
SUBROUTINE deallocate_dbm_from_sm(dbm,sm,copy_out,para_env,error)
    TYPE(real_matrix_type), POINTER          :: dbm, sm
    LOGICAL, INTENT(IN)                      :: copy_out
    TYPE(cp_para_env_type), INTENT(IN), &
      OPTIONAL                               :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_dbm_from_sm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_col, iblock_row
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block
    TYPE(real_block_node_type), POINTER      :: block_node

  CALL timeset(routineN,handle)

  IF ( copy_out ) THEN
     IF ( PRESENT(para_env) ) THEN
        DO iblock_row=1,dbm%nblock_row
          NULLIFY(block)
          CALL get_block_node(matrix=dbm,&
                              block_row=iblock_row,block_col=iblock_row,&
                              BLOCK=block)
          CALL mp_sum(block,para_env%group)
        END DO
      END IF
      DO iblock_row=1,sm%nblock_row
        block_node => first_block_node(sm,iblock_row)
        DO WHILE (ASSOCIATED(block_node))
          CALL get_block_node(block_node=block_node,&
                              block_col=iblock_col)
          IF (iblock_row == iblock_col) THEN
            CALL get_block_node(matrix=dbm,&
                                block_row=iblock_row,block_col=iblock_row,&
                                BLOCK=block)
            CALL put_block_node(sm,iblock_row,iblock_row,block,.TRUE.,error=error)
          END IF
          block_node => next_block_node(block_node)
        END DO
      END DO
  END IF

  CALL deallocate_matrix(dbm,error=error)

  CALL timestop(handle)

END SUBROUTINE deallocate_dbm_from_sm

! *****************************************************************************
!> \brief Adds a constant value to each matrix element
!> \param matrix_a a sparse matrix (inout)
!> \param alpha value to be added
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2008 created [mguidon]
!> \author Manuel Guidon 
! *****************************************************************************
SUBROUTINE cp_sm_add_constant(matrix_a,alpha,error)
    TYPE(real_matrix_type), POINTER          :: matrix_a
    REAL(KIND=dp), INTENT(in)                :: alpha
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_add_constant', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iblock_row, j
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a_block
    TYPE(real_block_node_type), POINTER      :: block_node

  CALL timeset(routineN,handle)

  NULLIFY(a_block,block_node)
  failure=.FALSE.
  CPPrecondition(ASSOCIATED(matrix_a),cp_failure_level,routineP,error,failure)

  DO iblock_row=1,matrix_a%nblock_row
   block_node => first_block_node(matrix_a,iblock_row)
    DO WHILE (ASSOCIATED(block_node))

      CALL get_block_node(block_node=block_node,&
                          BLOCK=a_block)
      
      DO j=1,SIZE(a_block,2)
        DO i=1,SIZE(a_block,1)
          a_block(i,j) = a_block(i,j) + alpha
        END DO
      END DO
      block_node => block_node%next_block_node
    END DO
  END DO
  CALL timestop(handle)
END SUBROUTINE cp_sm_add_constant

END MODULE sparse_matrix_types
