!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/matrix_types [1.0] *
!!
!!   NAME
!!     matrix_types
!!
!!   FUNCTION
!!     Define the matrix data types.
!!
!!   AUTHOR
!!     Matthias Krack (23.06.2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE matrix_types

! *****************************************************************************

! Index:

! SUBROUTINE add_real_matrix_block(matrix,block_row,block_col,block)
! SUBROUTINE allocate_real_matrix(matrix,nrow,ncol,nblock_row,nblock_col,&
!                                 first_row,last_row,first_col,last_col,&
!                                 matrix_name,matrix_symmetry)
! SUBROUTINE allocate_full_real_matrix(matrix,nrow,ncol,matrix_name,&
!                                      matrix_symmetry)
! SUBROUTINE copy_real_matrix(source,target)
! SUBROUTINE deallocate_real_matrix(matrix)
! SUBROUTINE deallocate_real_matrix_row(matrix,block_row)
! SUBROUTINE get_matrix_info(matrix,nrow,ncol,nblock_row,nblock_col,&
!                            matrix_name,matrix_symmetry)
! SUBROUTINE get_real_block_node(block_node,block_col,block)
! SUBROUTINE get_real_matrix_block(matrix,block_row,block_col,block_node,block)
! SUBROUTINE put_real_block_node(block_node,matrix,block_row,block_col,block)
! SUBROUTINE put_real_matrix_block(matrix,block_row,block_col,block)
! SUBROUTINE replicate_real_matrix(source,target,target_name)
! SUBROUTINE replicate_real_matrix_structure(source,target,target_name)

! FUNCTION find_real_block_node(matrix,block_row,block_col) RESULT(block_node)
! FUNCTION first_real_block_node(matrix,block_row) RESULT(first_block_node)
! FUNCTION next_real_block_node(block_node) RESULT(next_block_node)

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  IMPLICIT NONE

  PRIVATE

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

  TYPE real_matrix_set_type
    TYPE(real_matrix_type), POINTER :: matrix
  END TYPE real_matrix_set_type

! *** Public data types ***

  PUBLIC :: real_block_node_type,&
            real_matrix_set_type,&
            real_matrix_type

! *** Public subroutines ***

  PUBLIC :: add_block_node,&
            allocate_matrix,&
            copy_matrix,&
            deallocate_matrix,&
            deallocate_matrix_row,&
            get_block_node,&
            get_matrix_info,&
            put_block_node,&
            replicate_matrix,&
            replicate_matrix_structure

! *** Public functions ***

  PUBLIC :: first_block_node,&
            next_block_node

! *****************************************************************************

  INTERFACE add_block_node
    MODULE PROCEDURE add_real_matrix_block
  END INTERFACE

  INTERFACE allocate_matrix
    MODULE PROCEDURE allocate_real_matrix,allocate_full_real_matrix
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

  INTERFACE find_block_node
    MODULE PROCEDURE find_real_block_node
  END INTERFACE

  INTERFACE first_block_node
    MODULE PROCEDURE first_real_block_node
  END INTERFACE

  INTERFACE get_block_node
    MODULE PROCEDURE get_real_block_node,get_real_matrix_block
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

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE add_real_matrix_block(matrix,block_row,block_col,block)

!   Purpose: Allocate and initialize a new block node.

!   History: - Creation (28.07.2000, Matthias Krack)

!   ***************************************************************************

    USE memory_utilities, ONLY: reallocate
    USE termination,      ONLY: stop_memory,stop_program

    TYPE(real_matrix_type), POINTER             :: matrix
    INTEGER, INTENT(IN)                         :: block_col,block_row
    REAL(wp), DIMENSION(:,:), OPTIONAL, POINTER :: block

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE add_real_matrix_block (MODULE matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: istat,ncol,nrow

!   ---------------------------------------------------------------------------

!   *** Calculate the block dimensions ***

    nrow = matrix%last_row(block_row) - matrix%first_row(block_row) + 1
    ncol = matrix%last_col(block_col) - matrix%first_col(block_col) + 1

    ALLOCATE (block_node,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"Block node of "//matrix%name,0)

!   *** Define the data set of the new block node ***

    block_node%block_col = block_col

    NULLIFY (block_node%block)

    block_node%block => reallocate(block_node%block,1,nrow,1,ncol)

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

    block_node%next_block_node => matrix%block_list(block_row)%first_block_node

    matrix%block_list(block_row)%first_block_node => block_node

  END SUBROUTINE add_real_matrix_block

! *****************************************************************************

  SUBROUTINE allocate_real_matrix(matrix,nrow,ncol,nblock_row,nblock_col,&
                                  first_row,last_row,first_col,last_col,&
                                  matrix_name,matrix_symmetry)

!   Purpose: Allocate and initialize a real matrix.

!   History: - Creation (16.06.2000, Matthias Krack)

!   ***************************************************************************

    USE memory_utilities, ONLY: reallocate
    USE termination,      ONLY: stop_memory

    TYPE(real_matrix_type), POINTER :: matrix
    CHARACTER(LEN=*), INTENT(IN)    :: matrix_name,matrix_symmetry
    INTEGER, INTENT(IN)             :: nblock_row,nblock_col,ncol,nrow

    INTEGER, DIMENSION(nblock_row), INTENT(IN) :: first_row,last_row
    INTEGER, DIMENSION(nblock_col), INTENT(IN) :: first_col,last_col

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_real_matrix (MODULE matrix_types)"

!   *** Local variables ***

    INTEGER :: irow,istat

!   ---------------------------------------------------------------------------

!   *** Deallocate the old matrix ***

    IF (ASSOCIATED(matrix)) CALL deallocate_real_matrix(matrix)

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

    NULLIFY (matrix%first_row)
    matrix%first_row => reallocate(matrix%first_row,1,nblock_row)
    matrix%first_row(:) = first_row(:)

    NULLIFY (matrix%last_row)
    matrix%last_row => reallocate(matrix%last_row,1,nblock_row)
    matrix%last_row(:) = last_row(:)

    NULLIFY (matrix%first_col)
    matrix%first_col => reallocate(matrix%first_col,1,nblock_col)
    matrix%first_col(:) = first_col(:)

    NULLIFY (matrix%last_col)
    matrix%last_col => reallocate(matrix%last_col,1,nblock_col)
    matrix%last_col(:) = last_col(:)

!   *** Initialize all block lists ***

    DO irow=1,nblock_row
      NULLIFY (matrix%block_list(irow)%first_block_node)
    END DO

  END SUBROUTINE allocate_real_matrix

! *****************************************************************************

  SUBROUTINE allocate_full_real_matrix(matrix,nrow,ncol,matrix_name,&
                                       matrix_symmetry)

!   Purpose: Allocate and initialize a real matrix structure with only one
!            block.

!   History: - Creation (10.07.2000, Matthias Krack)

!   ***************************************************************************

    USE memory_utilities, ONLY: reallocate
    USE termination,      ONLY: stop_memory

    TYPE(real_matrix_type), POINTER :: matrix
    CHARACTER(LEN=*), INTENT(IN)    :: matrix_name,matrix_symmetry
    INTEGER, INTENT(IN)             :: ncol,nrow

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_full_real_matrix (MODULE matrix_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

!   *** Deallocate the old matrix ***

    IF (ASSOCIATED(matrix)) CALL deallocate_real_matrix(matrix)

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

    NULLIFY (matrix%first_row)
    matrix%first_row => reallocate(matrix%first_row,1,1)
    matrix%first_row(1) = 1

    NULLIFY (matrix%last_row)
    matrix%last_row => reallocate(matrix%last_row,1,1)
    matrix%last_row(1) = nrow

    NULLIFY (matrix%first_col)
    matrix%first_col => reallocate(matrix%first_col,1,1)
    matrix%first_col(1) = 1

    NULLIFY (matrix%last_col)
    matrix%last_col => reallocate(matrix%last_col,1,1)
    matrix%last_col(1) = ncol

    NULLIFY (matrix%block_list(1)%first_block_node)

    CALL add_real_matrix_block(matrix=matrix,&
                               block_row=1,&
                               block_col=1)

  END SUBROUTINE allocate_full_real_matrix

! *****************************************************************************

  SUBROUTINE copy_real_matrix(source,target)

!   Purpose: Copy the real matrix source to the existing real matrix target.

!   History: - Creation (07.07.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_program

    TYPE(real_matrix_type), POINTER :: source,target

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE copy_real_matrix (MODULE matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: source_block_node,&
                                           target_block_node

    INTEGER :: first_col,first_row,iblock_col,iblock_row,icol,irow,istat,&
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

  END SUBROUTINE copy_real_matrix

! *****************************************************************************

  SUBROUTINE deallocate_real_matrix(matrix)

!   Purpose: Deallocate a matrix.

!   History: - Creation (11.07.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_memory

    TYPE(real_matrix_type), POINTER :: matrix

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_real_matrix (MODULE matrix_types)"

!   *** Local variables ***

    INTEGER :: iblock_row,istat

!   ---------------------------------------------------------------------------

    DO iblock_row=1,matrix%nblock_row
      CALL deallocate_real_matrix_row(matrix,iblock_row)
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

!   History: - Creation (30.11.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_memory

    TYPE(real_matrix_type), POINTER :: matrix
    INTEGER, INTENT(IN)             :: block_row

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_real_matrix_row (MODULE matrix_types)"

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

  FUNCTION find_real_block_node(matrix,block_row,block_col) RESULT(block_node)

!   Purpose: Return a pointer to the requested block node.

!   History: - Creation (23.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER :: matrix
    INTEGER, INTENT(IN)             :: block_col,block_row

    TYPE(real_block_node_type), POINTER :: block_node

!   ---------------------------------------------------------------------------

    block_node => matrix%block_list(block_row)%first_block_node

    DO WHILE (ASSOCIATED(block_node))
      IF (block_node%block_col == block_col) RETURN
      block_node => block_node%next_block_node
    END DO

  END FUNCTION find_real_block_node

! *****************************************************************************

  FUNCTION first_real_block_node(matrix,block_row) RESULT(first_block_node)

!   Purpose: Return a pointer to the first block node of a block list.

!   History: - Creation (23.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER :: matrix
    INTEGER, INTENT(IN)             :: block_row

    TYPE(real_block_node_type), POINTER :: first_block_node

!   ---------------------------------------------------------------------------

    first_block_node => matrix%block_list(block_row)%first_block_node

  END FUNCTION first_real_block_node

! *****************************************************************************

  SUBROUTINE get_matrix_info(matrix,nrow,ncol,nblock_row,nblock_col,&
                             matrix_name,matrix_symmetry)

!   Purpose: Return the requested matrix information.

!   History: - Creation (10.07.2000, Matthias Krack)

!   ***************************************************************************

    USE memory_utilities, ONLY: reallocate

    TYPE(real_matrix_type), POINTER          :: matrix
    CHARACTER(LEN=60), OPTIONAL, INTENT(OUT) :: matrix_name
    CHARACTER(LEN=40), OPTIONAL, INTENT(OUT) :: matrix_symmetry
    INTEGER, OPTIONAL, INTENT(OUT)           :: nblock_row,nblock_col,ncol,nrow

!   ---------------------------------------------------------------------------

    IF (PRESENT(matrix_name)) matrix_name = matrix%name
    IF (PRESENT(matrix_symmetry)) matrix_symmetry = matrix%symmetry
    IF (PRESENT(nblock_row)) nblock_row = matrix%nblock_row
    IF (PRESENT(nblock_col)) nblock_col = matrix%nblock_col
    IF (PRESENT(nrow)) nrow = matrix%nrow
    IF (PRESENT(ncol)) ncol = matrix%ncol

  END SUBROUTINE get_matrix_info

! *****************************************************************************

  SUBROUTINE get_real_block_node(block_node,block_col,block)

!   Purpose: Get block node data set.

!   History: - Creation (28.05.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(real_block_node_type), INTENT(IN)      :: block_node
    INTEGER, OPTIONAL, INTENT(OUT)              :: block_col
    REAL(wp), DIMENSION(:,:), OPTIONAL, POINTER :: block

!   ---------------------------------------------------------------------------

    IF (PRESENT(block_col)) block_col = block_node%block_col
    IF (PRESENT(block)) block => block_node%block

  END SUBROUTINE get_real_block_node

! *****************************************************************************

  SUBROUTINE get_real_matrix_block(matrix,block_row,block_col,block_node,block)

!   Purpose: Get block node data set.

!   History: - Creation (28.05.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER               :: matrix
    TYPE(real_block_node_type), OPTIONAL, POINTER :: block_node
    INTEGER, INTENT(IN)                           :: block_col,block_row
    REAL(wp), DIMENSION(:,:), OPTIONAL, POINTER   :: block

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: current_block_node

!   ---------------------------------------------------------------------------

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

  FUNCTION next_real_block_node(block_node) RESULT(next_block_node)

!   Purpose: Return a pointer to the next block node of a block list.

!   History: - Creation (28.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(real_block_node_type), INTENT(IN) :: block_node

    TYPE(real_block_node_type), POINTER :: next_block_node

!   ---------------------------------------------------------------------------

    next_block_node => block_node%next_block_node

  END FUNCTION next_real_block_node

! *****************************************************************************

  SUBROUTINE put_real_block_node(block_node,matrix,block_row,block_col,block)

!   Purpose: Set block node data set.

!   History: - Creation (27.07.2000, Matthias Krack)

!   ***************************************************************************

    USE memory_utilities, ONLY: reallocate
    USE termination,      ONLY: stop_program

    TYPE(real_block_node_type), POINTER         :: block_node
    TYPE(real_matrix_type), POINTER             :: matrix
    INTEGER, INTENT(IN)                         :: block_row
    INTEGER, OPTIONAL, INTENT(IN)               :: block_col
    REAL(wp), DIMENSION(:,:), OPTIONAL, POINTER :: block

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE put_real_block_node (MODULE matrix_types)"

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

!   History: - Creation (07.07.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_program

    TYPE(real_matrix_type), POINTER            :: matrix
    INTEGER, INTENT(IN)                        :: block_col,block_row
    REAL(wp), DIMENSION(:,:), OPTIONAL, POINTER:: block

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE put_real_matrix_block (MODULE matrix_types)"

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

  SUBROUTINE replicate_real_matrix(source,target,target_name)

!   Purpose: Replicate the existing matrix source. The replicated matrix is
!            target with the name target_name.

!   History: - Creation (17.11.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_program

    TYPE(real_matrix_type), POINTER :: source,target
    CHARACTER(LEN=*), INTENT(IN)    :: target_name

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE replicate_real_matrix (MODULE matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: iblock_col,iblock_row

    REAL(wp), DIMENSION(:,:), POINTER :: block

!   ---------------------------------------------------------------------------

!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(source)) THEN
      CALL stop_program(routine,"The source matrix pointer is not associated")
    END IF

!   *** Allocate a new matrix structure ***

    CALL allocate_real_matrix(matrix=target,&
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

!   *** Initialize all block nodes ***

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

  END SUBROUTINE replicate_real_matrix

! *****************************************************************************

  SUBROUTINE replicate_real_matrix_structure(source,target,target_name)

!   Purpose: Replicate the matrix structure of the existing matrix source. The
!            new matrix target with the name target_name has the same
!            structure.

!   History: - Creation (17.11.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_program

    TYPE(real_matrix_type), POINTER :: source,target
    CHARACTER(LEN=*), INTENT(IN)    :: target_name

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE replicate_real_matrix_structure (MODULE matrix_types)"

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: iblock_col,iblock_row

!   ---------------------------------------------------------------------------

!   *** Check the association status of the source matrix ***

    IF (.NOT.ASSOCIATED(source)) THEN
      CALL stop_program(routine,"The source matrix pointer is not associated")
    END IF

!   *** Allocate a new matrix structure ***

    CALL allocate_real_matrix(matrix=target,&
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

END MODULE matrix_types
