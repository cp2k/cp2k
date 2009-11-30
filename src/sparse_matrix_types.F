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
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE message_passing,                 ONLY: mp_sum
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
            allocate_matrix,&
            allocate_matrix_set,&
            checksum_matrix, &
            copy_local_sm_to_replicated_fm,&
            deallocate_matrix,&
            deallocate_matrix_set,&
            get_block_node,&
            get_matrix_diagonal,&
            get_matrix_info,&
            replicate_matrix, replicate_matrix_structure,& ! should be put together
            set_matrix,&
            set_matrix_diagonal,&
            sparse_times_local, &
            sparse_plus_loc_loct, &
            cp_sm_scale_and_add,&
            cp_sm_get_id_nr,&
            cp_sm_set

! *** Public functions ***

  PUBLIC :: first_block_node,&
            next_block_node,&
            find_block_node

  INTERFACE add_block_node
    MODULE PROCEDURE add_real_matrix_block_obsolete
  END INTERFACE

  INTERFACE allocate_matrix
    MODULE PROCEDURE allocate_real_matrix
  END INTERFACE

  INTERFACE allocate_matrix_set
    MODULE PROCEDURE allocate_real_matrix_set
    MODULE PROCEDURE allocate_real_matrix_set_2d
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

  INTERFACE replicate_matrix
    MODULE PROCEDURE replicate_real_matrix
  END INTERFACE

  INTERFACE replicate_matrix_structure
    MODULE PROCEDURE replicate_real_matrix_structure
  END INTERFACE

  INTERFACE set_matrix
    MODULE PROCEDURE set_real_matrix
  END INTERFACE

  INTERFACE set_matrix_diagonal
    MODULE PROCEDURE set_real_matrix_diagonal
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
!> \brief   Allocate and initialize a new block node.
!> \author  MK
!> \date    28.07.2000
!>  DONT USE THIS ROUTINE ANYMORE! WILL BE REMOVED SOON.
!>  DONT USE THIS ROUTINE ANYMORE! WILL BE REMOVED SOON.
!>  DONT USE THIS ROUTINE ANYMORE! WILL BE REMOVED SOON.
!>  DONT USE THIS ROUTINE ANYMORE! Use the new sparse library dbcsr*.
!> \version 1.0
! *****************************************************************************
  SUBROUTINE add_real_matrix_block_obsolete(matrix,block_row,block_col,BLOCK, error)

    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(IN)                      :: block_row, block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: BLOCK
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'add_real_matrix_block_obsolete', &
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
            !vw turn off warnings (sm->dbcsr)
            !vwCPPostcondition(.FALSE.,cp_warning_level,routineP,error,failure)
         ENDIF
       CASE DEFAULT
         IF ( .NOT. ( rbr.EQ.myprow .AND. cbc.EQ.mypcol) ) THEN
            !vw turn off warnings (sm->dbcsr)
            !vwCPPostcondition(.FALSE.,cp_warning_level,routineP,error,failure)
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

  END SUBROUTINE add_real_matrix_block_obsolete

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

! *****************************************************************************
!> \brief   Allocate and initialize a real matrix at the real_matrix_type level.
!>          distribution_2d if present and not NULL, implies that this matrix has
!>          a distribution_2d compatible shape
!> \author  MK
!> \date    16.06.2000
!> \version 1.1
!> \par History
!>      2009-08-27 UB Adds flag to prevent repointing the matrix pointer
! *****************************************************************************
  SUBROUTINE allocate_real_matrix(matrix,nrow,ncol,nblock_row,nblock_col,&
       first_row,last_row,first_col,last_col,matrix_name,matrix_symmetry,&
       sparsity_id,distribution_2d,error,keep_pointer)
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
    LOGICAL, INTENT(in), OPTIONAL            :: keep_pointer

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_real_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: irow, istat
    LOGICAL                                  :: existed, kp

    IF (PRESENT (keep_pointer)) THEN
       kp = .TRUE.
    ELSE
       kp = .FALSE.
    ENDIF
    existed = ASSOCIATED (matrix)
!   *** Deallocate the old matrix ***

    IF (existed) THEN
       CALL deallocate_matrix(matrix,error=error,keep_pointer=kp)
    ENDIF
    IF (.NOT. (existed .AND. kp)) THEN
       ALLOCATE (matrix,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,matrix%name,0)
    ENDIF

!   *** Allocate a set of block lists ***

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
!> \brief  Deallocate a real matrix at the real_matrix_type level. 
!> \author MK 
!> \date   11.07.2000 
!> \version 1.1
!> \par History
!>      2009-08-27 UB Adds flag to prevent deallocating pointer
! *****************************************************************************
  SUBROUTINE deallocate_real_matrix(matrix,error,keep_pointer)

    TYPE(real_matrix_type), POINTER          :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL, INTENT(in), OPTIONAL            :: keep_pointer

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

    IF (PRESENT (keep_pointer)) THEN
       IF (keep_pointer) RETURN
    ENDIF
    DEALLOCATE (matrix,STAT=istat)
    NULLIFY (matrix)
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

END MODULE sparse_matrix_types
