!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief represent a complex full matrix distribued on many processors
!> \author Joost VandeVondele, based on Fawzi's cp_fm_* routines
! **************************************************************************************************
MODULE cp_cfm_types

   USE cp_blacs_calls,                  ONLY: cp_blacs_zgebr2d,&
                                              cp_blacs_zgebs2d
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_retain,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE cp_log_handling,                 ONLY: cp_to_string
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: cp2k_is_parallel,&
                                              mp_sum
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_cfm_types'

   INTEGER, SAVE, PRIVATE :: last_cfm_id_nr = 0

   PUBLIC :: cp_cfm_type, cp_cfm_p_type
   PUBLIC :: cp_cfm_create, cp_cfm_release
   PUBLIC :: cp_cfm_set_all, cp_cfm_get_info, cp_cfm_get_element
   PUBLIC :: cp_cfm_set_element, cp_cfm_to_cfm, cp_cfm_to_fm

   INTERFACE cp_cfm_to_cfm
      MODULE PROCEDURE cp_cfm_to_cfm_matrix, & ! a full matrix
         cp_cfm_to_cfm_columns ! just a number of columns
   END INTERFACE

! **************************************************************************************************
!> \brief represent a complex full matrix
!> \param name the name of the matrix, used for printing
!> \param id_nr identification number that makes this matrix unique
!> \param ref_cont reference count (see doc/ReferenceCounting.html)
!> \param print_count sequential number to make the printed files unique
!> \param matrix_struct structure of this matrix
!> \param local_data array with the data of the matrix (its contents
!>        depend on the matrix type used: in parallel runs it will be
!>        in scalapack format, in sequential, it will simply contain
!>        the matrix)
! **************************************************************************************************
   TYPE cp_cfm_type
!    PRIVATE
      CHARACTER(LEN=60) :: name
      INTEGER :: id_nr, ref_count, print_count
      TYPE(cp_fm_struct_type), POINTER :: matrix_struct
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER :: local_data
   END TYPE cp_cfm_type

! **************************************************************************************************
!> \brief just to build arrays of pointers to matrixes
!> \param matrix the pointer to the matrix
! **************************************************************************************************
   TYPE cp_cfm_p_type
      TYPE(cp_cfm_type), POINTER :: matrix
   END TYPE cp_cfm_p_type

CONTAINS

! **************************************************************************************************
!> \brief creates a new full matrix with the given structure
!> \param matrix the matrix to be created
!> \param matrix_struct the structure of matrix
!> \param name ...
!> \note
!>      preferred allocation routine
! **************************************************************************************************
   SUBROUTINE cp_cfm_create(matrix, matrix_struct, name)
      TYPE(cp_cfm_type), POINTER               :: matrix
      TYPE(cp_fm_struct_type), POINTER         :: matrix_struct
      CHARACTER(len=*), INTENT(in), OPTIONAL   :: name

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_create', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: ncol_local, npcol, nprow, &
                                                  nrow_local
      TYPE(cp_blacs_env_type), POINTER         :: context

#if defined(__parallel) && ! defined(__SCALAPACK)
      CPABORT("full matrixes need scalapack for parallel runs ")
#endif

      CPASSERT(ASSOCIATED(matrix_struct))
      ALLOCATE (matrix)
      context => matrix_struct%context
      matrix%matrix_struct => matrix_struct
      CALL cp_fm_struct_retain(matrix%matrix_struct)
      last_cfm_id_nr = last_cfm_id_nr+1
      matrix%id_nr = last_cfm_id_nr
      matrix%ref_count = 1
      matrix%print_count = 0

      nprow = context%num_pe(1)
      npcol = context%num_pe(2)
      NULLIFY (matrix%local_data)

      nrow_local = matrix_struct%local_leading_dimension
      ncol_local = MAX(1, matrix_struct%ncol_locals(context%mepos(2)))
      ALLOCATE (matrix%local_data(nrow_local, ncol_local))

      CALL dcopy(2*nrow_local*ncol_local, 0.0_dp, 0, matrix%local_data, 1)
! matrix%local_data(:,:) = 0.0_dp

      IF (PRESENT(name)) THEN
         matrix%name = name
      ELSE
         matrix%name = 'full matrix'//cp_to_string(matrix%id_nr)
      END IF
   END SUBROUTINE cp_cfm_create

! **************************************************************************************************
!> \brief retains a full matrix
!> \param matrix the matrix to retain
! **************************************************************************************************
   SUBROUTINE cp_cfm_retain(matrix)
      TYPE(cp_cfm_type), POINTER                         :: matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_retain', routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(matrix))
      CPASSERT(matrix%ref_count > 0)
      matrix%ref_count = matrix%ref_count+1
   END SUBROUTINE cp_cfm_retain

! **************************************************************************************************
!> \brief releases a full matrix
!> \param matrix the matrix to release
! **************************************************************************************************
   SUBROUTINE cp_cfm_release(matrix)
      TYPE(cp_cfm_type), POINTER                         :: matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_release', routineP = moduleN//':'//routineN

      IF (ASSOCIATED(matrix)) THEN
         CPASSERT(matrix%ref_count > 0)
         matrix%ref_count = matrix%ref_count-1
         IF (matrix%ref_count < 1) THEN
            IF (ASSOCIATED(matrix%local_data)) THEN
               DEALLOCATE (matrix%local_data)
            END IF
            matrix%name = ""
            CALL cp_fm_struct_release(matrix%matrix_struct)
            DEALLOCATE (matrix)
         END IF
      END IF
      NULLIFY (matrix)
   END SUBROUTINE cp_cfm_release

! **************************************************************************************************
!> \brief   Set the FULL matrix elements to alpha.
!> \param matrix ...
!> \param alpha ...
!> \param beta ...
!> \date    12.06.2001
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE cp_cfm_set_all(matrix, alpha, beta)
      TYPE(cp_cfm_type), POINTER                         :: matrix
      COMPLEX(KIND=dp), INTENT(IN)                       :: alpha
      COMPLEX(KIND=dp), INTENT(IN), OPTIONAL             :: beta

      INTEGER                                            :: i

      matrix%local_data(:, :) = alpha
      IF (PRESENT(beta)) THEN
         DO i = 1, matrix%matrix_struct%nrow_global
            CALL cp_cfm_set_element(matrix, i, i, beta)
         ENDDO
      ENDIF

   END SUBROUTINE cp_cfm_set_all

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param irow_global ...
!> \param icol_global ...
!> \param alpha ...
!> \par History
!>      , TCH, created
!>      always return the answer
! **************************************************************************************************
   SUBROUTINE cp_cfm_get_element(matrix, irow_global, icol_global, alpha)

      TYPE(cp_cfm_type), POINTER               :: matrix
      INTEGER, INTENT(IN)                      :: irow_global, icol_global
      COMPLEX(KIND=dp), INTENT(OUT)            :: alpha

      INTEGER                                  :: mypcol, myprow, npcol, nprow
      TYPE(cp_blacs_env_type), POINTER         :: context
#if defined(__SCALAPACK)
      COMPLEX(KIND=dp), DIMENSION(:, :), &
         POINTER                                :: a
      INTEGER                                  :: icol_local, ipcol, iprow, &
                                                  irow_local
      INTEGER, DIMENSION(9)                    :: desca
#endif

      context => matrix%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      nprow = context%num_pe(1)
      npcol = context%num_pe(2)

#if defined(__SCALAPACK)
      a => matrix%local_data
      desca(:) = matrix%matrix_struct%descriptor(:)

      CALL infog2l(irow_global, icol_global, desca, nprow, npcol, myprow, mypcol, &
                   irow_local, icol_local, iprow, ipcol)

      IF ((iprow == myprow) .AND. (ipcol == mypcol)) THEN
         alpha = a(irow_local, icol_local)
         CALL cp_blacs_ZGEBS2D(context%group, 'All', ' ', 1, 1, alpha, 1)
      ELSE
         CALL cp_blacs_ZGEBR2D(context%group, 'All', ' ', 1, 1, alpha, 1, iprow, ipcol)
      END IF
#else
      alpha = matrix%local_data(irow_global, icol_global)
#endif

   END SUBROUTINE cp_cfm_get_element

! **************************************************************************************************
!> \brief gets a submatrix of a full matrix
!>       op(target_m)(1:n_rows,1:n_cols)
!>       =fm(start_row:start_row+n_rows,start_col:start_col+n_cols)
!>      target_m is replicated on all cpus
!>      using this call is expensive
!> \param fm        : the full you want to get the info from
!> \param target_m  : a replicated full matrix that will contain the result
!> \param start_row : the starting row of b_matrix (defaults to 1)
!> \param start_col : the starting col of b_matrix (defaults to 1)
!> \param n_rows    : the number of row to change in b (defaults to
!>                    size(op(new_values),1))
!> \param n_cols    : the number of columns to change in b (defaults to
!>                    size(op(new_values),2))
!> \param transpose : if target_m should be transposed: if true
!>                    op(target_m)=target_m^T, else op(target_m)=target_m
!>                    (defaults to false)
!> \par History
!>      2016/04/05 created borrowing from Fawzi's cp_fm_get_submatrix
!> \author Lianheng Tong
!> \note
!>      optimized for full column updates. Zeros out a little too much
!>      of target_m
!>      the target_m is replicated and valid on all cpus
! **************************************************************************************************
   SUBROUTINE cp_cfm_get_submatrix(fm, &
                                   target_m, &
                                   start_row, &
                                   start_col, &
                                   n_rows, &
                                   n_cols, &
                                   transpose)
      TYPE(cp_cfm_type), POINTER                         :: fm
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(out)     :: target_m
      INTEGER, INTENT(in), OPTIONAL                      :: start_row, start_col, n_rows, n_cols
      LOGICAL, INTENT(in), OPTIONAL                      :: transpose

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_get_submatrix', &
         routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: full_block
      INTEGER                                            :: handle, i, i0, j, j0, ncol, ncol_global, &
                                                            ncol_local, nrow, nrow_global, &
                                                            nrow_local, this_col, this_row
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: tr_a
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      i0 = 1
      j0 = 1
      tr_a = .FALSE.

      CPASSERT(ASSOCIATED(fm))
      CPASSERT(fm%ref_count .GT. 0)

      IF (PRESENT(start_row)) i0 = start_row
      IF (PRESENT(start_col)) j0 = start_col
      IF (PRESENT(transpose)) tr_a = transpose
      IF (tr_a) THEN
         nrow = SIZE(target_m, 2)
         ncol = SIZE(target_m, 1)
      ELSE
         nrow = SIZE(target_m, 1)
         ncol = SIZE(target_m, 2)
      END IF
      IF (PRESENT(n_rows)) nrow = n_rows
      IF (PRESENT(n_cols)) ncol = n_cols

      para_env => fm%matrix_struct%para_env
      full_block => fm%local_data

#if defined(__SCALAPACK)
      ! zero-out whole target_m
      IF (SIZE(target_m, 1)*SIZE(target_m, 2) .NE. 0) THEN
         CALL dcopy(SIZE(target_m, 1)*SIZE(target_m, 2), 0.0_dp, 0, target_m(1, 1), 1)
      END IF
#endif
      CALL cp_cfm_get_info(matrix=fm, &
                           nrow_global=nrow_global, ncol_global=ncol_global, &
                           nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices)
      DO j = 1, ncol_local
         this_col = col_indices(j)-j0+1
         IF (this_col .GE. 1 .AND. this_col .LE. ncol) THEN
            IF (tr_a) THEN
               IF (i0 .EQ. 1 .AND. nrow_global .EQ. nrow) THEN
                  DO i = 1, nrow_local
                     target_m(this_col, row_indices(i)) = full_block(i, j)
                  END DO
               ELSE
                  DO i = 1, nrow_local
                     this_row = row_indices(i)-i0+1
                     IF (this_row .GE. 1 .AND. this_row .LE. nrow) THEN
                        target_m(this_col, this_row) = full_block(i, j)
                     END IF
                  END DO
               END IF
            ELSE
               IF (i0 .EQ. 1 .AND. nrow_global .EQ. nrow) THEN
                  DO i = 1, nrow_local
                     target_m(row_indices(i), this_col) = full_block(i, j)
                  END DO
               ELSE
                  DO i = 1, nrow_local
                     this_row = row_indices(i)-i0+1
                     IF (this_row .GE. 1 .AND. this_row .LE. nrow) THEN
                        target_m(this_row, this_col) = full_block(i, j)
                     END IF
                  END DO
               END IF
            END IF
         END IF
      END DO

      CALL mp_sum(target_m, para_env%group)

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_get_submatrix

! **************************************************************************************************
!> \brief   Set the FULL matrix element (irow_global,icol_global) to alpha.
!> \param matrix ...
!> \param irow_global ...
!> \param icol_global ...
!> \param alpha ...
!> \date    12.06.2001
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE cp_cfm_set_element(matrix, irow_global, icol_global, alpha)

      TYPE(cp_cfm_type), POINTER               :: matrix
      INTEGER, INTENT(IN)                      :: irow_global, icol_global
      COMPLEX(KIND=dp), INTENT(IN)             :: alpha

      INTEGER                                  :: mypcol, myprow, npcol, nprow
      TYPE(cp_blacs_env_type), POINTER         :: context
#if defined(__SCALAPACK)
      INTEGER                                  :: icol_local, ipcol, iprow, &
                                                  irow_local
      INTEGER, DIMENSION(9)                    :: desca
      COMPLEX(KIND=dp), DIMENSION(:, :), &
         POINTER                                :: a
#endif

      context => matrix%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      nprow = context%num_pe(1)
      npcol = context%num_pe(2)

#if defined(__SCALAPACK)

      a => matrix%local_data

      desca(:) = matrix%matrix_struct%descriptor(:)

      CALL infog2l(irow_global, icol_global, desca, nprow, npcol, myprow, mypcol, &
                   irow_local, icol_local, iprow, ipcol)

      IF ((iprow == myprow) .AND. (ipcol == mypcol)) THEN
         a(irow_local, icol_local) = alpha
      END IF

#else

      matrix%local_data(irow_global, icol_global) = alpha

#endif
   END SUBROUTINE cp_cfm_set_element

! **************************************************************************************************
!> \brief sets a submatrix of a full matrix
!>       matrix(start_row:start_row+n_rows,start_col:start_col+n_cols)
!>       = alpha*op(new_values)(1:n_rows,1:n_cols)+ beta
!>       * matrix(start_row:start_row+n_rows,start_col:start_col+n_cols)
!> \param matrix the full to change
!> \param new_values a replicated full matrix with the new values
!> \param start_row the starting row of b_matrix (defaults to 1)
!> \param start_col the starting col of b_matrix (defaults to 1)
!> \param n_rows the number of row to change in b (defaults to
!>        size(op(new_values),1))
!> \param n_cols the number of columns to change in b (defaults to
!>        size(op(new_values),2))
!> \param alpha rescaling factor for the new values (defaults to 1.0)
!> \param beta rescaling factor for the old values (defaults to 0.0)
!> \param transpose if new_values should be transposed: if true
!>        op(new_values)=new_values^T, else op(new_values)=new_values
!>        (defaults to false)
!> \par History
!>      07.2002 created borrowing from Fawzi's cp_fm_set_submatrix
!> \author Lianheng Tong
!> \note
!>      optimized for full column updates and alpha=1.0, beta=0.0
!>      the new_values need to be valid on all cpus
! **************************************************************************************************
   SUBROUTINE cp_cfm_set_submatrix(matrix, new_values, start_row, &
                                   start_col, n_rows, n_cols, alpha, beta, transpose)
      TYPE(cp_cfm_type), POINTER                         :: matrix
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(in)      :: new_values
      INTEGER, INTENT(in), OPTIONAL                      :: start_row, start_col, n_rows, n_cols
      COMPLEX(KIND=dp), INTENT(in), OPTIONAL             :: alpha, beta
      LOGICAL, INTENT(in), OPTIONAL                      :: transpose

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_set_submatrix', &
         routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: al, be
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: full_block
      INTEGER                                            :: i, i0, j, j0, ncol, ncol_block, &
                                                            ncol_global, ncol_local, nrow, &
                                                            nrow_block, nrow_global, nrow_local, &
                                                            this_col, this_row
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: tr_a

      al = (1.0_dp, 0.0_dp)
      be = (0.0_dp, 0.0_dp)
      i0 = 1
      j0 = 1
      tr_a = .FALSE.
      ! can be called too many times, making it a bit useless
      ! CALL timeset(routineN//','//moduleN,handle)

      CPASSERT(ASSOCIATED(matrix))
      CPASSERT(matrix%ref_count > 0)

      IF (PRESENT(alpha)) al = alpha
      IF (PRESENT(beta)) be = beta
      IF (PRESENT(start_row)) i0 = start_row
      IF (PRESENT(start_col)) j0 = start_col
      IF (PRESENT(transpose)) tr_a = transpose
      IF (tr_a) THEN
         nrow = SIZE(new_values, 2)
         ncol = SIZE(new_values, 1)
      ELSE
         nrow = SIZE(new_values, 1)
         ncol = SIZE(new_values, 2)
      END IF
      IF (PRESENT(n_rows)) nrow = n_rows
      IF (PRESENT(n_cols)) ncol = n_cols

      full_block => matrix%local_data

      CALL cp_cfm_get_info(matrix=matrix, &
                           nrow_global=nrow_global, &
                           ncol_global=ncol_global, &
                           nrow_block=nrow_block, &
                           ncol_block=ncol_block, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      IF (al .EQ. (1.0_dp, 0.0_dp) .AND. be .EQ. (0.0_dp, 0.0_dp)) THEN
         DO j = 1, ncol_local
            this_col = col_indices(j)-j0+1
            IF (this_col .GE. 1 .AND. this_col .LE. ncol) THEN
               IF (tr_a) THEN
                  IF (i0 .EQ. 1 .AND. nrow_global .EQ. nrow) THEN
                     DO i = 1, nrow_local
                        full_block(i, j) = new_values(this_col, row_indices(i))
                     END DO
                  ELSE
                     DO i = 1, nrow_local
                        this_row = row_indices(i)-i0+1
                        IF (this_row .GE. 1 .AND. this_row .LE. nrow) THEN
                           full_block(i, j) = new_values(this_col, this_row)
                        END IF
                     END DO
                  END IF
               ELSE
                  IF (i0 .EQ. 1 .AND. nrow_global .EQ. nrow) THEN
                     DO i = 1, nrow_local
                        full_block(i, j) = new_values(row_indices(i), this_col)
                     END DO
                  ELSE
                     DO i = 1, nrow_local
                        this_row = row_indices(i)-i0+1
                        IF (this_row .GE. 1 .AND. this_row .LE. nrow) THEN
                           full_block(i, j) = new_values(this_row, this_col)
                        END IF
                     END DO
                  END IF
               END IF
            END IF
         END DO
      ELSE
         DO j = 1, ncol_local
            this_col = col_indices(j)-j0+1
            IF (this_col .GE. 1 .AND. this_col .LE. ncol) THEN
               IF (tr_a) THEN
                  DO i = 1, nrow_local
                     this_row = row_indices(i)-i0+1
                     IF (this_row .GE. 1 .AND. this_row .LE. nrow) THEN
                        full_block(i, j) = al*new_values(this_col, this_row)+ &
                                           be*full_block(i, j)
                     END IF
                  END DO
               ELSE
                  DO i = 1, nrow_local
                     this_row = row_indices(i)-i0+1
                     IF (this_row .GE. 1 .AND. this_row .LE. nrow) THEN
                        full_block(i, j) = al*new_values(this_row, this_col)+ &
                                           be*full_block(i, j)
                     END IF
                  END DO
               END IF
            END IF
         END DO
      END IF

      ! CALL timestop(handle)

   END SUBROUTINE cp_cfm_set_submatrix

! **************************************************************************************************
!> \brief   Return informations about the specified FULL matrix.
!>          the indices vectors are cached.
!>          indices is an array that maps the local indices to the global ones
!>          local data is a little hole to access the blacs block, needed in some rare cases
!>          for efficient use
!> \param matrix ...
!> \param name ...
!> \param nrow_global ...
!> \param ncol_global ...
!> \param nrow_block ...
!> \param ncol_block ...
!> \param nrow_local ...
!> \param ncol_local ...
!> \param row_indices ...
!> \param col_indices ...
!> \param local_data ...
!> \param context ...
!> \param matrix_struct ...
!> \param para_env ...
!> \date    12.06.2001
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE cp_cfm_get_info(matrix, name, nrow_global, ncol_global, &
                              nrow_block, ncol_block, nrow_local, ncol_local, &
                              row_indices, col_indices, local_data, context, &
                              matrix_struct, para_env)
      TYPE(cp_cfm_type), POINTER        :: matrix
      CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: name
      INTEGER, OPTIONAL, INTENT(OUT)          :: ncol_block, ncol_global, &
                                                 nrow_block, nrow_global, &
                                                 nrow_local, ncol_local
      INTEGER, OPTIONAL, DIMENSION(:), POINTER   :: row_indices, col_indices
      TYPE(cp_para_env_type), POINTER, OPTIONAL :: para_env
      TYPE(cp_blacs_env_type), POINTER, OPTIONAL :: context
      TYPE(cp_fm_struct_type), POINTER, OPTIONAL :: matrix_struct
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER, OPTIONAL :: local_data

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_get_info', &
                                     routineP = moduleN//':'//routineN
      INTEGER i, nprow, npcol, myprow, mypcol
      TYPE(cp_blacs_env_type), POINTER :: ctxt
#if defined(__SCALAPACK)
      INTEGER, EXTERNAL :: indxl2g
#endif

      IF (PRESENT(name)) name = matrix%name
      IF (PRESENT(nrow_global)) nrow_global = matrix%matrix_struct%nrow_global
      IF (PRESENT(ncol_global)) ncol_global = matrix%matrix_struct%ncol_global
      IF (PRESENT(nrow_block)) nrow_block = matrix%matrix_struct%nrow_block
      IF (PRESENT(ncol_block)) ncol_block = matrix%matrix_struct%ncol_block
      IF (PRESENT(context)) context => matrix%matrix_struct%context
      IF (PRESENT(matrix_struct)) matrix_struct => matrix%matrix_struct
      IF (PRESENT(local_data)) local_data => matrix%local_data ! not hiding things anymore :-(

      ctxt => matrix%matrix_struct%context
      myprow = ctxt%mepos(1)
      mypcol = ctxt%mepos(2)
      nprow = ctxt%num_pe(1)
      npcol = ctxt%num_pe(2)

      IF (PRESENT(nrow_local)) THEN
         nrow_local = matrix%matrix_struct%nrow_locals(myprow)
      ENDIF

      IF (PRESENT(ncol_local)) THEN
         ncol_local = matrix%matrix_struct%ncol_locals(mypcol)
      ENDIF

      IF (PRESENT(row_indices)) THEN
         row_indices => matrix%matrix_struct%row_indices
         IF (.NOT. ASSOCIATED(row_indices)) THEN
            ALLOCATE (matrix%matrix_struct%row_indices &
                      (matrix%matrix_struct%nrow_locals(myprow)))
            row_indices => matrix%matrix_struct%row_indices
#ifdef __SCALAPACK
            DO i = 1, SIZE(row_indices)
               matrix%matrix_struct%row_indices(i) = &
                  indxl2g(i, matrix%matrix_struct%nrow_block, myprow, &
                          matrix%matrix_struct%first_p_pos(1), nprow)
            ENDDO
#else
            DO i = 1, SIZE(row_indices)
               row_indices(i) = i
            ENDDO
#endif
         END IF
      ENDIF

      IF (PRESENT(col_indices)) THEN
         col_indices => matrix%matrix_struct%col_indices
         IF (.NOT. ASSOCIATED(col_indices)) THEN
            ALLOCATE (matrix%matrix_struct%col_indices &
                      (matrix%matrix_struct%ncol_locals(mypcol)))
            col_indices => matrix%matrix_struct%col_indices
#ifdef __SCALAPACK
            DO i = 1, SIZE(col_indices)
               col_indices(i) = &
                  indxl2g(i, matrix%matrix_struct%ncol_block, mypcol, &
                          matrix%matrix_struct%first_p_pos(2), npcol)
            ENDDO
#else
            DO i = 1, SIZE(col_indices)
               col_indices(i) = i
            ENDDO
#endif
         END IF
      ENDIF

      IF (PRESENT(para_env)) para_env => matrix%matrix_struct%para_env

   END SUBROUTINE cp_cfm_get_info

! **************************************************************************************************
!> \brief copy one identically sized matrix in the other
!> \param source ...
!> \param destination ...
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE cp_cfm_to_cfm_matrix(source, destination)
      TYPE(cp_cfm_type), POINTER                         :: source, destination

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_to_cfm_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: npcol, nprow

      nprow = source%matrix_struct%context%num_pe(1)
      npcol = source%matrix_struct%context%num_pe(2)

      IF (.NOT. cp2k_is_parallel .OR. &
          cp_fm_struct_equivalent(source%matrix_struct, &
                                  destination%matrix_struct)) THEN
         IF (SIZE(source%local_data, 1) /= SIZE(destination%local_data, 1) .OR. &
             SIZE(source%local_data, 2) /= SIZE(destination%local_data, 2)) &
            CPABORT("internal local_data has different sizes")
         CALL dcopy(SIZE(source%local_data, 1)*SIZE(source%local_data, 2)*2, &
                    source%local_data(1, 1), 1, destination%local_data(1, 1), 1)
      ELSE
         IF (source%matrix_struct%nrow_global /= destination%matrix_struct%nrow_global) &
            CPABORT("cannot copy between full matrixes of differen sizes")
         IF (source%matrix_struct%ncol_global /= destination%matrix_struct%ncol_global) &
            CPABORT("cannot copy between full matrixes of differen sizes")
#ifdef __SCALAPACK
         CALL pzcopy(source%matrix_struct%nrow_global* &
                     source%matrix_struct%ncol_global, &
                     source%local_data(1, 1), 1, 1, source%matrix_struct%descriptor, 1, &
                     destination%local_data(1, 1), 1, 1, destination%matrix_struct%descriptor, 1)
#else
         CPABORT("")
#endif
      END IF
   END SUBROUTINE cp_cfm_to_cfm_matrix

! **************************************************************************************************
!> \brief copy just a subset of columns of a cfm to a cfm
!> \param msource ...
!> \param mtarget ...
!> \param ncol ...
!> \param source_start ...
!> \param target_start ...
! **************************************************************************************************
   SUBROUTINE cp_cfm_to_cfm_columns(msource, mtarget, ncol, source_start, &
                                    target_start)

      TYPE(cp_cfm_type), POINTER               :: msource, mtarget
      INTEGER, INTENT(IN)                      :: ncol
      INTEGER, INTENT(IN), OPTIONAL            :: source_start, target_start

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_cfm_to_cfm_columns', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: handle, n, ss, ts
      COMPLEX(KIND=dp), DIMENSION(:, :), &
         POINTER                                :: a, b
#if defined(__SCALAPACK)
      INTEGER                                  :: i
      INTEGER, DIMENSION(9)                    :: desca, descb
#endif

      CALL timeset(routineN, handle)

      ss = 1
      ts = 1

      IF (PRESENT(source_start)) ss = source_start
      IF (PRESENT(target_start)) ts = target_start

      n = msource%matrix_struct%nrow_global

      a => msource%local_data
      b => mtarget%local_data

#if defined(__SCALAPACK)
      desca(:) = msource%matrix_struct%descriptor(:)
      descb(:) = mtarget%matrix_struct%descriptor(:)
      DO i = 0, ncol-1
         CALL pzcopy(n, a(1, 1), 1, ss+i, desca, 1, b(1, 1), 1, ts+i, descb, 1)
      END DO
#else
      CALL zcopy(ncol*n, a(1, ss), 1, b(1, ts), 1)
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_to_cfm_columns

! **************************************************************************************************
!> \brief copy just a triangular matrix
!> \param msource ...
!> \param mtarget ...
!> \param uplo : 'U' for upper triangular, 'L' for lower triangular
! **************************************************************************************************
   SUBROUTINE cp_cfm_to_cfm_triangular(msource, mtarget, uplo)
      TYPE(cp_cfm_type), POINTER               :: msource, mtarget
      CHARACTER(LEN=*), INTENT(IN)             :: uplo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_cfm_to_cfm_triangular', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: handle, ncol, nrow
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER  :: aa, bb
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                    :: desca, descb
#endif

      CALL timeset(routineN, handle)

      nrow = msource%matrix_struct%nrow_global
      ncol = msource%matrix_struct%ncol_global

      aa => msource%local_data
      bb => mtarget%local_data

#if defined(__SCALAPACK)
      desca(:) = msource%matrix_struct%descriptor(:)
      descb(:) = mtarget%matrix_struct%descriptor(:)
      CALL pzlacpy(uplo, nrow, ncol, aa(1, 1), 1, 1, desca, bb(1, 1), 1, 1, descb)
#else
      CALL zlacpy(uplo, nrow, ncol, aa(1, 1), nrow, bb(1, 1), nrow)
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_to_cfm_triangular

! **************************************************************************************************
!> \brief copy real and imaginary parts of a cfm to fm
!>        works only for absolutely equivalent structures
!> \param msource  : complex matrix
!> \param mtargetr : optional, if present outputs real part of msource
!> \param mtargeti : optional, if present outputs imag part of msource
! **************************************************************************************************
   SUBROUTINE cp_cfm_to_fm(msource, mtargetr, mtargeti)

      TYPE(cp_cfm_type), POINTER                         :: msource
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: mtargetr, mtargeti

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_cfm_to_fm', routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: zmat
      INTEGER                                            :: handle
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: imat, rmat

      CALL timeset(routineN, handle)

      zmat => msource%local_data
      IF (PRESENT(mtargetr)) THEN
         rmat => mtargetr%local_data
         IF ((.NOT. cp_fm_struct_equivalent(mtargetr%matrix_struct, msource%matrix_struct)) .OR. &
             (SIZE(rmat, 1) .NE. SIZE(zmat, 1)) .OR. &
             (SIZE(rmat, 2) .NE. SIZE(zmat, 2))) THEN
            CPABORT("size of local_data of mtargetr differ to msource")
         END IF
         ! copy local data
         rmat = REAL(zmat, KIND=dp)
      ELSE
         NULLIFY (rmat)
      END IF
      IF (PRESENT(mtargeti)) THEN
         imat => mtargeti%local_data
         IF ((.NOT. cp_fm_struct_equivalent(mtargeti%matrix_struct, msource%matrix_struct)) .OR. &
             (SIZE(imat, 1) .NE. SIZE(zmat, 1)) .OR. &
             (SIZE(imat, 2) .NE. SIZE(zmat, 2))) THEN
            CPABORT("size of local_data of mtargeti differ to msource")
         END IF
         ! copy local data
         imat = REAL(AIMAG(zmat), KIND=dp)
      ELSE
         NULLIFY (imat)
      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_to_fm

! **************************************************************************************************
!> \brief copy real and imaginary in forms of fms to a cfm
!>        works only for absolutely equivalent structures
!> \param msourcer : optional, if present input for real part.
!>                   if not present, treat as zero
!> \param msourcei : optional, if present input for imag part.
!>                   if not present, treat as zero
!> \param mtarget  : outputs complex matrix
! **************************************************************************************************
   SUBROUTINE cp_fm_to_cfm(msourcer, msourcei, mtarget)
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: msourcer, msourcei
      TYPE(cp_cfm_type), POINTER                         :: mtarget

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_to_cfm', routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: zmat
      INTEGER                                            :: handle, mode
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: imat, rmat

      CALL timeset(routineN, handle)

      mode = 0
      zmat => mtarget%local_data
      IF (PRESENT(msourcer)) THEN
         rmat => msourcer%local_data
         IF ((.NOT. cp_fm_struct_equivalent(msourcer%matrix_struct, mtarget%matrix_struct)) .OR. &
             (SIZE(rmat, 1) .NE. SIZE(zmat, 1)) .OR. &
             (SIZE(rmat, 2) .NE. SIZE(zmat, 2))) THEN
            CPABORT("size of local_data of msourcer differ to mtarget")
         END IF
         mode = mode+1
      ELSE
         NULLIFY (rmat)
      END IF
      IF (PRESENT(msourcei)) THEN
         imat => msourcei%local_data
         IF ((.NOT. cp_fm_struct_equivalent(msourcei%matrix_struct, mtarget%matrix_struct)) .OR. &
             (SIZE(imat, 1) .NE. SIZE(zmat, 1)) .OR. &
             (SIZE(imat, 2) .NE. SIZE(zmat, 2))) THEN
            CPABORT("size of local_data of msourcei differ to mtarget")
         END IF
         mode = mode+2
      ELSE
         NULLIFY (imat)
      END IF
      ! copy local data
      SELECT CASE (mode)
      CASE (0)
         zmat(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      CASE (1)
         zmat(:, :) = CMPLX(rmat(:, :), 0.0_dp, KIND=dp)
      CASE (2)
         zmat(:, :) = CMPLX(0.0_dp, imat(:, :), KIND=dp)
      CASE (3)
         zmat = CMPLX(rmat, imat, KIND=dp)
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE cp_fm_to_cfm

END MODULE cp_cfm_types
