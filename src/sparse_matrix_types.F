!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004 CP2K developers group                           !
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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_type,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_retain,&
                                             distribution_2d_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size,&
                                             default_string_length
  USE message_passing,                 ONLY: mp_max,mp_sum
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: moduleN='sparse_matrix_types'
  INTEGER, SAVE, PRIVATE :: last_matrix_id=0

! *** Define a real matrix data type ***

  TYPE real_block_node_type
    PRIVATE
    TYPE(real_block_node_type), POINTER :: next_block_node
    INTEGER                             :: block_col
    REAL(KIND = dp), DIMENSION(:,:), POINTER   :: BLOCK
  END TYPE real_block_node_type

  TYPE real_block_list_type
    PRIVATE
    TYPE(real_block_node_type), POINTER :: first_block_node
    INTEGER                             :: nblock_node
    TYPE(real_block_node_type), POINTER :: last_used_block_node ! not last in list
  END TYPE real_block_list_type

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
            cp_sm_maxabsval

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
    MODULE PROCEDURE allocate_real_matrix 
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

! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****f* sparse_matrix_types/cp_sm_scale_and_add [1.0] *
!!
!!   NAME
!!     cp_sm_scale_and_add
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_scale_and_add(matrix_a, alpha, matrix_b, beta,&
!!         error)
!!       Type(real_matrix_type), Pointer:: matrix_a
!!       Type(real_matrix_type), Pointer, Optional:: matrix_b
!!       Real(KIND = dp), Intent (IN), Optional:: alpha, beta
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_scale_and_add
!!
!!   FUNCTION
!!     scales the matrix and add another matrix
!!     matrix_a=alpha*matrix_a+beta*matrix_b
!!     filters the result with the structure (sparsity) of matrix_a
!!
!!   NOTES
!!     should be rewitten avoiding the get_block_node taht calls find_block_node
!!
!!   INPUTS
!!     - matrix_a: a sparse matrix (inout)
!!     - alpha: scaling factor of matrix_a (defaults to 1.0)
!!     - matrix_b: a sparse matrix (in)
!!     - beta: scaling factor of matrix_b (defaults to 1.0)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     2.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sm_scale_and_add(matrix_a,alpha,matrix_b,beta,error)
    TYPE(real_matrix_type), POINTER          :: matrix_a
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: alpha
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: matrix_b
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: beta
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_scale_and_add', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_col, iblock_row
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_alpha, my_beta
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a_block, b_block
    TYPE(real_block_node_type), POINTER      :: block_node

  CALL timeset(routineN,"I","",handle)

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
  CALL timestop(0.0_dp,handle)
END SUBROUTINE cp_sm_scale_and_add
!***************************************************************************

  SUBROUTINE add_real_matrices(c,alpha,a,beta,b)

!   Purpose: Calculate the sum of the sparse matrices a and b which have the
!            same structure.

!            c = alpha*a + beta*b

!   History: - Creation (10.04.2002, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: c
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(real_matrix_type), POINTER          :: a
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(real_matrix_type), POINTER          :: b

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE add_real_matrices (MODULE sparse_matrix_types)"

    INTEGER                                  :: iblock_col, iblock_row
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a_block, b_block, c_block
    TYPE(real_block_node_type), POINTER      :: block_node

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

  SUBROUTINE add_real_matrix_block(matrix,block_row,block_col,BLOCK,&
       error)

!   Purpose: Allocate and initialize a new block node.

!   History: - Creation (28.07.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row, block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: BLOCK
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_real_matrix_block', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cbc, cbr, istat, mypcol, &
                                                myprow, ncol, nrow, rbc, rbr
    LOGICAL                                  :: failure
    TYPE(real_block_node_type), POINTER      :: block_node, new_block_node

!   ---------------------------------------------------------------------------
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
!*****************************************************************************************

  SUBROUTINE  add_1d_block_node(matrix,block_row,block_col,BLOCK,error)

!   Purpose: Allocate and initialize a new block node in the oce matrix

! *****************************************************************************

    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row, block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: BLOCK
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_1d_block_node', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cbc, cbr, istat, mypcol, &
                                                myprow, ncol, nrow, rbc, rbr
    LOGICAL                                  :: failure
    TYPE(real_block_node_type), POINTER      :: block_node, new_block_node
!   ---------------------------------------------------------------------------
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

!*****************************************************************************************
! writes for all local blocks the maximum absolute value of 
! the block entries to unit fileunit
!*****************************************************************************************

  SUBROUTINE write_blocks_maxabsval(sparse_matrix,fileunit)

    TYPE(real_matrix_type), POINTER          :: sparse_matrix
    INTEGER, INTENT(IN)                      :: fileunit

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE write_block_magnitude (MODULE sparse_matrix_types)"

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
  SUBROUTINE allocate_real_matrix(matrix,nrow,ncol,nblock_row,nblock_col,&
                                  first_row,last_row,first_col,last_col,&
                                  matrix_name,matrix_symmetry,sparsity_id,distribution_2d)

!   Purpose: Allocate and initialize a real matrix at the real_matrix_type
!            level.

!   distribution_2d if present and not NULL, implies that this matrix has a distribution_2d compatible shape

!   History: - Creation (16.06.2000, MK)

!   ***************************************************************************

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

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE allocate_real_matrix (MODULE sparse_matrix_types)"

    INTEGER                                  :: irow, istat

!   ---------------------------------------------------------------------------
!   *** Deallocate the old matrix ***

    IF (ASSOCIATED(matrix)) CALL deallocate_matrix(matrix)

!   *** Allocate a set of block lists ***

    ALLOCATE (matrix,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,matrix%name,0)

    NULLIFY(matrix%distribution_2d) 
    IF (PRESENT(distribution_2d)) matrix%distribution_2d=>distribution_2d
    IF (ASSOCIATED(matrix%distribution_2d)) THEN
       CALL distribution_2d_retain(matrix%distribution_2d)
    END IF

    ALLOCATE (matrix%block_list(nblock_row),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,matrix%name,0)
    
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
      matrix%block_list(irow)%nblock_node = 0
      NULLIFY (matrix%block_list(irow)%last_used_block_node)
    END DO

  END SUBROUTINE allocate_real_matrix

! *****************************************************************************

  SUBROUTINE allocate_real_matrix_set(matrix_set,nmatrix)

!   Purpose: Allocate and initialize a real matrix set.

!   History: - Creation (13.03.2002, MK)

!   ***************************************************************************

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_set
    INTEGER, INTENT(IN)                      :: nmatrix

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE allocate_real_matrix_set (MODULE sparse_matrix_types)"

    INTEGER                                  :: imatrix, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(matrix_set)) CALL deallocate_matrix_set(matrix_set)

    ALLOCATE (matrix_set(nmatrix),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"matrix_set",0)

    DO imatrix=1,nmatrix
      NULLIFY (matrix_set(imatrix)%matrix)
    END DO

  END SUBROUTINE allocate_real_matrix_set
!******************************************************************
  SUBROUTINE copy_local_sm_to_replicated_fm(sparse_matrix,fm)

!   Purpose: copy to local blocks of a sparse matrix in a replicated full matrix (why ??)
!            ! according to me this routine is a bug in parallel JVDV 

!   History: - Creation (19.06.2001, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: sparse_matrix
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fm

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE copy_local_sm_to_replicated_fm (MODULE sparse_matrix_types)"

    INTEGER                                  :: iblock_col, iblock_row, icol, &
                                                icol_block, irow, irow_block, &
                                                istat, ncol, nrow
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sparse_block
    TYPE(real_block_node_type), POINTER      :: block_node

!   ---------------------------------------------------------------------------
!   *** Check the association status of the input matrix ***

    IF (.NOT.ASSOCIATED(sparse_matrix)) THEN
      CALL stop_program(routine,"The input matrix pointer is not associated")
    END IF

    IF (ASSOCIATED(fm)) DEALLOCATE (fm)

    nrow = sparse_matrix%nrow
    ncol = sparse_matrix%ncol

    ALLOCATE (fm(nrow,ncol),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"fm",nrow*ncol*dp_size)
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

  END SUBROUTINE copy_local_sm_to_replicated_fm

! *****************************************************************************
! multiplies a sparse matrix times ncol local vectors, adding it the the output
! (that are stored as rows ! in v_in leading dimension nblock )
! if the sparse_matrix is symmetric, off diagonal blocks will be used twice
! it is still assumed that the diagonal blocks are symmetric (and hence full)
! Joost VandeVondele july 2002
! *****************************************************************************
  ! notice, adding a beta is incompatible with the definition of this function
  SUBROUTINE sparse_times_local(sparse_matrix,rep_v_in,rep_v_out,ncol,nblock, alpha)
    TYPE(real_matrix_type), POINTER          :: sparse_matrix
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rep_v_in, rep_v_out
    INTEGER, INTENT(IN)                      :: ncol, nblock
    REAL(KIND=dp), OPTIONAL                  :: alpha

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE sparse_times_local (MODULE sparse_matrix_types)"

    INTEGER :: handle, iblock_col, iblock_row, ithread, nthread, &
      omp_get_num_threads, omp_get_thread_num, &
      sbncol, sbnrow
    LOGICAL                                  :: antisymmetric, owner, &
                                                symmetric
    REAL(KIND=dp)                            :: alpha_l, alpha_l_ji, flops
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sparse_block
    TYPE(real_block_node_type), POINTER      :: block_node

!   ---------------------------------------------------------------------------

    IF (ncol.EQ.0) RETURN

    IF (PRESENT(alpha)) THEN
        alpha_l=alpha
    ELSE
        alpha_l=1.0_dp
    ENDIF

!   *** Check the association status of the input matrix ***
    CALL timeset("sparse_times_local","I","Gflops",handle)

    IF (.NOT.ASSOCIATED(sparse_matrix)) THEN
      CALL stop_program(routine,"The input matrix pointer is not associated")
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
       CALL stop_program(routine,"wrong matrix symmetry specification")
    END SELECT

    IF ((symmetric .OR. antisymmetric) &
           .AND. (sparse_matrix%nrow .NE. sparse_matrix%ncol)) THEN
       CALL stop_program(routine,"error nonsquare symmetric matrix")
    END IF
    IF (sparse_matrix%ncol .NE. SIZE(rep_v_in,2)) THEN
      CALL stop_program(routine,"n x n * n x k ?")
    ENDIF

    flops=0.0_dp

!   *** Traverse all block nodes of the sparse matrix ***
!   no simple parallel do because we write both the iblock_row and iblock_col

!$OMP PARALLEL PRIVATE(ithread,nthread,iblock_row,block_node,iblock_col) &
!$OMP PRIVATE(sparse_block,owner,sbnrow,sbncol) REDUCTION(+:flops) &
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
           IF (MOD(iblock_row,nthread).eq.ithread) THEN
              flops=flops+2.0_dp*sbnrow*sbncol*ncol
              CALL DGEMM('N','T',ncol,sbnrow,sbncol,alpha_l, &
                         rep_v_in (1,sparse_matrix%first_col(iblock_col)), &
                         nblock, &
                         sparse_block(1,1),sbnrow, &
                         1.0_dp,rep_v_out(1,sparse_matrix%first_row(iblock_row)), &
                         nblock)
           ENDIF
           
           IF (iblock_col .NE. iblock_row .AND. (symmetric .OR. antisymmetric)) THEN
              IF (MOD(iblock_col,nthread).eq.ithread) THEN 
                 flops=flops+2.0_dp*sbnrow*sbncol*ncol
                 CALL DGEMM('N','N',ncol,sbncol,sbnrow,alpha_l_ji, &
                      rep_v_in (1,sparse_matrix%first_row(iblock_row)), &
                      nblock, &
                      sparse_block(1,1),sbnrow, &
                      1.0_dp,rep_v_out(1,sparse_matrix%first_col(iblock_col)), &
                      nblock)
              ENDIF 
           ENDIF
        END IF

        block_node => next_block_node(block_node)

      END DO

    END DO
!$OMP END PARALLEL

    flops=flops*1E-9_dp
    CALL timestop(flops,handle)

  END SUBROUTINE sparse_times_local

! *****************************************************************************
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

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE sparse_plus_loc_loct(MODULE sparse_matrix_types)"

    INTEGER :: handle, iblock_col, iblock_row, icol, irow, ithread, nthread, &
       omp_get_num_threads, omp_get_thread_num, sbncol, sbnrow
    REAL(KIND=dp)                            :: flops, prefactor
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sparse_block
    TYPE(real_block_node_type), POINTER      :: block_node

!   ---------------------------------------------------------------------------
!   *** Check the association status of the input matrix ***

    IF (ncol.EQ.0) RETURN

    IF (.NOT.ASSOCIATED(sparse_matrix)) THEN
      CALL stop_program(routine,"The input matrix pointer is not associated")
    END IF

    IF (PRESENT(alpha)) THEN
      prefactor = alpha
    ELSE
      prefactor = 1.0_dp
    END IF

    IF (sparse_matrix%nrow .NE. sparse_matrix%ncol) THEN
      CALL stop_program(routine,"not allowed")
    ENDIF
    CALL timeset("sparse_plus_loc_loct","I","Gflops",handle)
    flops=0.0_dp

!   *** Traverse all block nodes of the sparse matrix ***
! take into account the triangular schape of the matrix
!$OMP PARALLEL PRIVATE(ithread,nthread,iblock_row,block_node,iblock_col,sparse_block,icol,irow,sbnrow,sbncol) REDUCTION(+:flops)
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

           flops=flops+2.0_dp*sbnrow*sbncol*ncol
           CALL DGEMM('T','N',sbnrow,sbncol,ncol,prefactor,&
                              rep_v(1,irow),nblock, &
                              rep_g(1,icol),nblock, &
                       1.0_dp,sparse_block(1,1),sbnrow)
        ENDIF

        block_node => next_block_node(block_node)

      END DO

    END DO
!$OMP END PARALLEL

    flops=flops*1E-9_dp
    CALL timestop(flops,handle)

  END SUBROUTINE sparse_plus_loc_loct

! *****************************************************************************

! *****************************************************************************

  SUBROUTINE copy_real_matrix(source,TARGET)

!   Purpose: Copy the real matrix source to the existing real matrix target.

!   History: - Creation (07.07.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: source, TARGET

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE copy_real_matrix (MODULE sparse_matrix_types)"

    INTEGER :: first_col, first_row, handle, iblock_col, iblock_row, icol, &
      irow, jblock_col, jblock_row, jcol, jrow, last_col, last_row, &
      source_first_col, source_first_row, source_last_col, source_last_row, &
      target_first_col, target_first_row, target_last_col, target_last_row
    TYPE(real_block_node_type), POINTER      :: source_block_node, &
                                                target_block_node

!   ---------------------------------------------------------------------------
!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(source)) THEN
      CALL stop_program(routine,"The source matrix pointer is not associated")
    END IF

    CALL timeset("copy_real_matrix","I","",handle)

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

      CALL stop_program(routine,&
                        "The target matrix pointer is not associated. "//&
                        "Use replicate_matrix instead of copy_matrix")

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE copy_real_matrix

! *****************************************************************************

  SUBROUTINE deallocate_real_matrix(matrix)

!   Purpose: Deallocate a real matrix at the real_matrix_type level.

!   History: - Creation (11.07.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: matrix

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE deallocate_real_matrix (MODULE sparse_matrix_types)"

    INTEGER                                  :: iblock_row, istat

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

    CALL distribution_2d_release(matrix%distribution_2d)

    DEALLOCATE (matrix,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,matrix%name)

  END SUBROUTINE deallocate_real_matrix

! *****************************************************************************

  SUBROUTINE deallocate_real_matrix_row(matrix,block_row)

!   Purpose: Deallocate a matrix row.

!   History: - Creation (30.11.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE deallocate_real_matrix_row (MODULE sparse_matrix_types)"

    INTEGER                                  :: istat
    TYPE(real_block_node_type), POINTER      :: current_block_node, &
                                                next_block_node

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
    NULLIFY (matrix%block_list(block_row)%last_used_block_node)

    matrix%block_list(block_row)%nblock_node = 0

  END SUBROUTINE deallocate_real_matrix_row

! *****************************************************************************

  SUBROUTINE deallocate_real_matrix_set(matrix_set)

!   Purpose: Deallocate a real matrix set.

!   History: - Creation (13.03.2002, MK)

!   ***************************************************************************

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_set

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE deallocate_real_matrix_set (MODULE sparse_matrix_types)"

    INTEGER                                  :: imatrix, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(matrix_set)) THEN
      DO imatrix=1,SIZE(matrix_set)
        CALL deallocate_matrix(matrix_set(imatrix)%matrix)
      END DO
      DEALLOCATE (matrix_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"matrix_set")
    END IF

  END SUBROUTINE deallocate_real_matrix_set

! *****************************************************************************

  FUNCTION find_real_block_node(matrix,block_row,block_col) RESULT(block_node)

!   Purpose: Return a pointer to the requested block node.

!   History: - Creation (23.06.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row, block_col
    TYPE(real_block_node_type), POINTER      :: block_node

    CHARACTER(len=*), PARAMETER :: routineN = 'find_real_block_node', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lb, ub
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: error

!   ---------------------------------------------------------------------------

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

  FUNCTION first_real_block_node(matrix,block_row) RESULT(first_block_node)

!   Purpose: Return a pointer to the first block node of a block list.

!   History: - Creation (23.06.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row
    TYPE(real_block_node_type), POINTER      :: first_block_node

!   ---------------------------------------------------------------------------

    first_block_node => matrix%block_list(block_row)%first_block_node

  END FUNCTION first_real_block_node

! *****************************************************************************

  SUBROUTINE get_matrix_info(matrix,matrix_name,matrix_symmetry,&
                             nblock_row,nblock_col,nrow,ncol,&
                             first_row,last_row,first_col,last_col,&
                             nblock_allocated,nelement_allocated,sparsity_id)

!   Purpose: Return the requested matrix information.

!   History: - Creation (10.07.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: matrix
    CHARACTER(LEN=80), INTENT(OUT), OPTIONAL :: matrix_name
    CHARACTER(LEN=40), INTENT(OUT), OPTIONAL :: matrix_symmetry
    INTEGER, INTENT(OUT), OPTIONAL           :: nblock_row, nblock_col, nrow, &
                                                ncol
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: first_row, last_row, &
                                                first_col, last_col
    INTEGER, INTENT(OUT), OPTIONAL           :: nblock_allocated, &
                                                nelement_allocated,sparsity_id

    INTEGER                                  :: iblock_row
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: BLOCK
    TYPE(real_block_node_type), POINTER      :: block_node

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

  SUBROUTINE get_real_block_node(block_node,block_col,BLOCK)

!   Purpose: Get the the column and/or matrix-block for a given sparse matrix block node.

!   History: - Creation (28.05.2000, MK)

!   ***************************************************************************

    TYPE(real_block_node_type), INTENT(IN)   :: block_node
    INTEGER, INTENT(OUT), OPTIONAL           :: block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: BLOCK

!   ---------------------------------------------------------------------------

    IF (PRESENT(block_col)) block_col = block_node%block_col
    IF (PRESENT(BLOCK)) BLOCK => block_node%block

  END SUBROUTINE get_real_block_node

! *****************************************************************************

  SUBROUTINE get_real_matrix_block(matrix,block_row,block_col,&
                                   first_row,last_row,first_col,last_col,&
                                   block_node,BLOCK)

!   Purpose: Get a specific (block_row,block_col)-block from the real matrix

!   History: - Creation (28.05.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row, block_col
    INTEGER, INTENT(OUT), OPTIONAL           :: first_row, last_row, &
                                                first_col, last_col
    TYPE(real_block_node_type), OPTIONAL, &
      POINTER                                :: block_node
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: BLOCK

    TYPE(real_block_node_type), POINTER      :: current_block_node

!   ---------------------------------------------------------------------------

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

  SUBROUTINE get_real_matrix_diagonal(matrix,vector)

!   Purpose: Set the diagonal elements of matrix to the values of vector.

!   History: - Creation (06.11.2003,MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: matrix
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: vector

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE get_real_matrix_diagonal (MODULE sparse_matrix_types)"

    INTEGER                                  :: first_row, i, iblock_col, &
                                                iblock_row, irow
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block
    TYPE(real_block_node_type), POINTER      :: block_node

!   ---------------------------------------------------------------------------
!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(matrix)) THEN
      CALL stop_program(routine,"The matrix pointer is not associated")
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

  FUNCTION next_real_block_node(block_node) RESULT(next_block_node)

!   Purpose: Return a pointer to the next block node of a block list.

!   History: - Creation (28.06.2000, MK)

!   ***************************************************************************

    TYPE(real_block_node_type), POINTER      :: block_node, next_block_node

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(block_node)) THEN
       next_block_node => block_node%next_block_node
    ELSE
       CALL stop_program("sparse_matrix_types:next_real_block_node",&
            "The block_node pointer a is not associated")
       NULLIFY(next_block_node)
    END IF

  END FUNCTION next_real_block_node

! *****************************************************************************

  SUBROUTINE put_real_block_node(block_node,matrix,block_row,block_col,BLOCK)

!   Purpose: Set block node data set.

!   History: - Creation (27.07.2000, MK)

!   ***************************************************************************

    TYPE(real_block_node_type), POINTER      :: block_node
    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row
    INTEGER, INTENT(IN), OPTIONAL            :: block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: BLOCK

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE put_real_block_node (MODULE sparse_matrix_types)"

    INTEGER                                  :: ncol, nrow

!   ---------------------------------------------------------------------------

    IF (PRESENT(block_col)) block_node%block_col = block_col

    IF (PRESENT(BLOCK)) THEN

      nrow = matrix%last_row(block_row) - matrix%first_row(block_row) + 1
      ncol = matrix%last_col(block_col) - matrix%first_col(block_col) + 1

      IF (PRESENT(BLOCK)) THEN
        IF (ASSOCIATED(BLOCK)) THEN
          IF ((SIZE(BLOCK,1) == nrow).AND.(SIZE(BLOCK,2) == ncol)) THEN
            block_node%block(:,:) = BLOCK(:,:)
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

  SUBROUTINE put_real_matrix_block(matrix,block_row,block_col,BLOCK)

!   Purpose: Allocate and initialize a new block node or update an existing
!            block node.

!   History: - Creation (07.07.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row, block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: BLOCK

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE put_real_matrix_block (MODULE sparse_matrix_types)"

    INTEGER                                  :: ncol, nrow
    TYPE(real_block_node_type), POINTER      :: block_node

!   ---------------------------------------------------------------------------
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
            block_node%block(:,:) = BLOCK(:,:)
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

      IF (PRESENT(BLOCK)) THEN
        CALL add_block_node(matrix,block_row,block_col,BLOCK)
      ELSE
        CALL add_block_node(matrix,block_row,block_col)
      END IF

    END IF

  END SUBROUTINE put_real_matrix_block

! *****************************************************************************

  SUBROUTINE replicate_real_matrix(source,TARGET,target_name,allocate_blocks)

!   Purpose: Replicate the existing matrix source. The replicated matrix is
!            target with the name target_name.
!            by default copies the data blocks,
!            if optional allocate_blocks=false then no blocks are allocated

!   History: - Creation (17.11.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: source, TARGET
    CHARACTER(LEN=*), INTENT(IN)             :: target_name
    LOGICAL, INTENT(IN), OPTIONAL            :: allocate_blocks

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE replicate_real_matrix (MODULE sparse_matrix_types)"

    INTEGER                                  :: iblock_col, iblock_row
    LOGICAL                                  :: make_data_blocks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: BLOCK
    TYPE(real_block_node_type), POINTER      :: block_node

!   ---------------------------------------------------------------------------
!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(source)) THEN
      CALL stop_program(routine,"The source matrix pointer is not associated")
    END IF

    IF (ASSOCIATED(TARGET)) CALL deallocate_matrix(TARGET)

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
                              BLOCK=BLOCK)

          CALL add_block_node(matrix=TARGET,&
                              block_row=iblock_row,&
                              block_col=iblock_col,&
                              BLOCK=BLOCK)

          block_node => next_block_node(block_node)

        END DO

      END DO

    ENDIF

  END SUBROUTINE replicate_real_matrix

! *****************************************************************************

  SUBROUTINE replicate_real_matrix_structure(source,TARGET,target_name,&
                                             target_symmetry)

!   Purpose: Replicate the matrix structure of the existing matrix source. The
!            new matrix target with the name target_name has the same
!            structure.

!   History: - Creation (17.11.2000, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: source, TARGET
    CHARACTER(LEN=*), INTENT(IN)             :: target_name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: target_symmetry

    CHARACTER(LEN=*), PARAMETER :: routine = "replicate_real_matrix_structure"

    CHARACTER(LEN=40)                        :: matrix_symmetry
    INTEGER                                  :: handle, iblock_col, iblock_row
    TYPE(real_block_node_type), POINTER      :: block_node

!   ---------------------------------------------------------------------------
!   *** Check the association status of the source matrix ***

    CALL timeset(routine,"I","",handle)

    IF (.NOT.ASSOCIATED(source)) THEN
      CALL stop_program(routine,"The source matrix pointer is not associated")
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
                         matrix_symmetry=matrix_symmetry)

!   *** Initialize all block nodes ***

    DO iblock_row=1,source%nblock_row

      block_node => first_block_node(source,iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col)

        CALL add_block_node(matrix=TARGET,&
                            block_row=iblock_row,&
                            block_col=iblock_col)

        block_node => next_block_node(block_node)

      END DO

    END DO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE replicate_real_matrix_structure

! *****************************************************************************

  SUBROUTINE scale_real_matrix(a,alpha)

!   Purpose: Multiply the sparse matrix a with alpha.

!            a <- alpha*a

!   History: - Creation (17.07.2002,MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: a
    REAL(KIND=dp), INTENT(IN)                :: alpha

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE scale_real_matrix (MODULE sparse_matrix_types)"

    INTEGER                                  :: iblock_row
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: BLOCK
    TYPE(real_block_node_type), POINTER      :: block_node

!   ---------------------------------------------------------------------------
!   *** Check the association status of the matrices ***

    IF (.NOT.ASSOCIATED(a)) THEN
       CALL stop_program(routine,"The matrix pointer a is not associated")
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

  SUBROUTINE set_real_matrix(matrix,value)

!   Purpose: Set all elements of matrix to value.

!   History: - Creation (10.04.2002, MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: matrix
    REAL(KIND=dp), INTENT(IN)                :: value

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE set_real_matrix (MODULE sparse_matrix_types)"

    INTEGER                                  :: iblock_row
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: BLOCK
    TYPE(real_block_node_type), POINTER      :: block_node

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
                            BLOCK=BLOCK)

        BLOCK(:,:) = value

        block_node => block_node%next_block_node

      END DO

    END DO

  END SUBROUTINE set_real_matrix

! *****************************************************************************

  SUBROUTINE set_real_matrix_diagonal(matrix,vector)

!   Purpose: Set the diagonal elements of matrix to the values of vector.

!   History: - Creation (06.11.2003,MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: matrix
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: vector

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE set_real_matrix_diagonal (MODULE sparse_matrix_types)"

    INTEGER                                  :: first_row, i, iblock_col, &
                                                iblock_row, irow
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block
    TYPE(real_block_node_type), POINTER      :: block_node

!   ---------------------------------------------------------------------------
!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(matrix)) THEN
      CALL stop_program(routine,"The matrix pointer is not associated")
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

!!****f* sparse_matrix_types/transfer_real_matrix [1.0] *
!!
!!   NAME
!!     transfer_real_matrix
!!
!!   FUNCTION
!!   Transfers source matrix to target matrix (copies the values) 
!!   filtering the result through the structure (sparsity) of target matrix
!!   zeros blocks as necessary
!!
!!   NOTES
!!     assumes that both matrixes have the same distribution
!!     XXXXX very dangerous (i.e. wrong?) in parallell, if the sparsities are not the same XXXXXXX
!!
!!   INPUTS
!!     source_matrix: the matrix to copy
!!     target_matrix: the matrix that will contain the copy
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     - 10.04.2002 Creation [MK]
!!     - 4.2003 extended for matrixes with different structures [fawzi]
!!
!!*** *********************************************************************
  SUBROUTINE transfer_real_matrix(source_matrix,target_matrix ,error)

    TYPE(real_matrix_type), POINTER          :: source_matrix, target_matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'transfer_real_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle 
    INTEGER                                  :: iblock_row
    INTEGER                                  :: iblock_col_source
    INTEGER                                  :: iblock_col_target
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: source_block, target_block
    TYPE(real_block_node_type), POINTER      :: block_node_target, block_node_source

!   ---------------------------------------------------------------------------
!   *** Check the association status of the source matrix ***

    CALL timeset(routineN,"I","",handle)
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
              CALL dcopy(SIZE(target_block,1)*SIZE(target_block,2),&
                   source_block(1,1),1,target_block(1,1),1)
           ELSE
              CALL dcopy(SIZE(target_block,1)*SIZE(target_block,2),&
                   0.0_dp,0,target_block(1,1),1)
           END IF
        END IF
        
        block_node_target => block_node_target%next_block_node

      END DO

    END DO
    CALL timestop(0.0_dp,handle)

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
    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(in)                      :: block_row, block_col
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
!***************************************************************************

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

!!****f* sparse_matrix_types/cp_sm_get_id_nr [1.0] *
!!
!!   NAME
!!     cp_sm_get_id_nr
!!
!!   SYNOPSIS
!!     Function cp_sm_get_id_nr(matrix, error) Result(res)
!!       Integer:: res
!!       Type(real_matrix_type), Intent (IN), Target:: matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_sm_get_id_nr
!!
!!   FUNCTION
!!     returns the id of the given matrix
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - matrix: the matrix you want info about
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
FUNCTION cp_sm_get_id_nr(matrix,error) RESULT(res)
    TYPE(real_matrix_type), INTENT(in), &
      TARGET                                 :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
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
!***************************************************************************

!!****f* sparse_matrix_types/cp_sm_set [1.0] *
!!
!!   NAME
!!     cp_sm_set
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_set(matrix, name, error)
!!       Type(real_matrix_type), Pointer:: matrix
!!       Character(Len=*), Intent (IN), Optional:: name
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_set
!!
!!   FUNCTION
!!     sets various attributes of the given matrix
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sm_set(matrix,name,error)
    TYPE(real_matrix_type), POINTER          :: matrix
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

!***************************************************************************
!!****f* sparse_matrix_types/cp_sm_get_block_list [1.0] *
!!
!!   NAME
!!     cp_sm_get_block_list
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_get_block_list(matrix, list, error)
!!       Type(real_matrix_type), Pointer:: matrix
!!       Integer, Dimension(:,:), Pointer:: list
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_get_block_list
!!
!!   FUNCTION
!!     gets a list of indices of all allocated blocks (locally) of a
!!     given matrix
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     JGH (04-07-2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
SUBROUTINE cp_sm_get_block_list ( matrix, list, error )
    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, DIMENSION(:, :), POINTER        :: list
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

!***************************************************************************

!!****f* cp_sm_maxabsval [1.0] *
!!
!!   NAME
!!     cp_sm_maxabsval
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_maxabsval(matrix, maxabsval, para_env, error)
!!       Type(real_matrix_type), Pointer:: matrix
!!       Real(KIND = dp), Intent (OUT):: maxabsval
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_maxabsval
!!
!!   FUNCTION
!!     returns the maximum of the abs of the matrix entries
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - matrix: the matrix you are interested in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sm_maxabsval(matrix,maxabsval,para_env,error)
    TYPE(real_matrix_type), POINTER          :: matrix
    REAL(KIND=dp), INTENT(out)               :: maxabsval
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
!***************************************************************************
!!****f* cp_sm_sm_trace [1.0] *
!!
!!   NAME
!!     cp_sm_sm_trace
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_sm_trace(matrix1, matrix2, trace, para_env, error)
!!       Type(real_matrix_type), Pointer:: matrix1
!!       Type(real_matrix_type), Pointer:: matrix2
!!       Real(KIND = dp), Intent (OUT):: trace
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_sm_trace
!!
!!   FUNCTION
!!     returns the trace of the product of two sparse matrices
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - matrix: the matrix you are interested in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     jgh
!!
!!   MODIFICATION HISTORY
!!     10.2004 created [jgh]
!!
!!*** **********************************************************************
SUBROUTINE cp_sm_sm_trace(matrix1,matrix2,trace,para_env,error)
    TYPE(real_matrix_type), POINTER          :: matrix1,matrix2
    REAL(KIND=dp), INTENT(out)               :: trace
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_sm_trace',&
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iblock_row,nblock_row,&
                                                block_col1,block_col2
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix_block1
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix_block2
    TYPE(real_block_node_type), POINTER      :: block_node1,block_node2
    REAL(KIND=dp)                            :: trace_block, sfactor
    CHARACTER(len=default_string_length)     :: matrix_symmetry

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(matrix1),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(matrix2),cp_failure_level,routineP,error,failure)
  CPPrecondition(matrix1%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(matrix2%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
  
  trace=0.0_dp

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

        trace_block = SUM ( matrix_block1*matrix_block2 )

        IF ( iblock_row == block_col1 ) THEN
          trace = trace + trace_block
        ELSE
          trace = trace + sfactor*trace_block
        END IF

        block_node1 => block_node1%next_block_node
        block_node2 => block_node2%next_block_node

     END DO

  END DO

  CALL mp_sum(trace,para_env%group)

END SUBROUTINE cp_sm_sm_trace

!***************************************************************************

END MODULE sparse_matrix_types
