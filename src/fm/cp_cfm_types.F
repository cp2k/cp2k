!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Represents a complex full matrix distributed on many processors.
!> \author Joost VandeVondele, based on Fawzi's cp_fm_* routines
! **************************************************************************************************
MODULE cp_cfm_types
   USE cp_blacs_calls,                  ONLY: cp_blacs_zgebr2d,&
                                              cp_blacs_zgebs2d
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent,&
                                              cp_fm_struct_get,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_retain,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE cp_log_handling,                 ONLY: cp_to_string
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE message_passing,                 ONLY: &
        cp2k_is_parallel, mp_allgather, mp_any_source, mp_bcast, mp_irecv, mp_isend, mp_proc_null, &
        mp_request_null, mp_sum, mp_waitall
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_cfm_types'
   INTEGER, PARAMETER, PRIVATE :: src_tag = 3, dest_tag = 5, send_tag = 7, recv_tag = 11

   INTEGER, SAVE, PRIVATE :: last_cfm_id_nr = 0

   PUBLIC :: cp_cfm_type, cp_cfm_p_type, copy_cfm_info_type
   PUBLIC :: cp_cfm_cleanup_copy_general, &
             cp_cfm_create, &
             cp_cfm_finish_copy_general, &
             cp_cfm_get_element, &
             cp_cfm_get_info, &
             cp_cfm_get_submatrix, &
             cp_cfm_release, &
             cp_cfm_retain, &
             cp_cfm_set_all, &
             cp_cfm_set_element, &
             cp_cfm_set_submatrix, &
             cp_cfm_start_copy_general, &
             cp_cfm_to_cfm, &
             cp_cfm_to_fm, &
             cp_fm_to_cfm

   INTERFACE cp_cfm_to_cfm
      MODULE PROCEDURE cp_cfm_to_cfm_matrix, & ! a full matrix
         cp_cfm_to_cfm_columns ! just a number of columns
   END INTERFACE

! **************************************************************************************************
!> \brief Represent a complex full matrix.
!> \param name           the name of the matrix, used for printing
!> \param id_nr         identification number that makes this matrix unique
!> \param ref_cont      reference count (see doc/ReferenceCounting.html)
!> \param print_count   sequential number to make the printed files unique
!> \param matrix_struct structure of this matrix
!> \param local_data    array with the data of the matrix (its content depends on
!>                      the matrix type used: in parallel run it will be in
!>                      ScaLAPACK format, in sequential run it will simply contain the matrix)
! **************************************************************************************************
   TYPE cp_cfm_type
      CHARACTER(len=60) :: name
      INTEGER :: id_nr, ref_count, print_count
      TYPE(cp_fm_struct_type), POINTER :: matrix_struct
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER :: local_data
   END TYPE cp_cfm_type

! **************************************************************************************************
!> \brief Just to build arrays of pointers to matrices.
!> \param matrix the pointer to the matrix
! **************************************************************************************************
   TYPE cp_cfm_p_type
      TYPE(cp_cfm_type), POINTER :: matrix
   END TYPE cp_cfm_p_type

! **************************************************************************************************
!> \brief Stores the state of a copy between cp_cfm_start_copy_general
!>        and cp_cfm_finish_copy_general.
!> \par History
!>      Jan 2017  derived type 'copy_info_type' has been created [Mark T]
!>      Jan 2018  the type 'copy_info_type' has been adapted for complex matrices [Sergey Chulkov]
! **************************************************************************************************
   TYPE copy_cfm_info_type
      !> number of MPI processes that send data
      INTEGER                                     :: send_size
      !> number of locally stored rows (1) and columns (2) of the destination matrix
      INTEGER, DIMENSION(2)                       :: nlocal_recv
      !> number of rows (1) and columns (2) of the ScaLAPACK block of the source matrix
      INTEGER, DIMENSION(2)                       :: nblock_src
      !> BLACS process grid shape of the source matrix: (1) nproc_row, (2) nproc_col
      INTEGER, DIMENSION(2)                       :: src_num_pe
      !> displacements into recv_buf
      INTEGER, ALLOCATABLE, DIMENSION(:)          :: recv_disp
      !> MPI requests for non-blocking receive and send operations
      INTEGER, ALLOCATABLE, DIMENSION(:)          :: recv_request, send_request
      !> global column and row indices of locally stored elements of the destination matrix
      INTEGER, DIMENSION(:), POINTER              :: recv_col_indices, recv_row_indices
      !> rank of MPI process with BLACS coordinates (prow, pcol)
      INTEGER, ALLOCATABLE, DIMENSION(:, :)       :: src_blacs2mpi
      !> receiving and sending buffers for non-blocking MPI communication
      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:) :: recv_buf, send_buf
   END TYPE copy_cfm_info_type

CONTAINS

! **************************************************************************************************
!> \brief Creates a new full matrix with the given structure.
!> \param matrix        matrix to be created
!> \param matrix_struct structure of the matrix
!> \param name          name of the matrix
!> \note
!>      preferred allocation routine
! **************************************************************************************************
   SUBROUTINE cp_cfm_create(matrix, matrix_struct, name)
      TYPE(cp_cfm_type), POINTER                         :: matrix
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      CHARACTER(len=*), INTENT(in), OPTIONAL             :: name

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_create', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: ncol_local, npcol, nprow, nrow_local
      TYPE(cp_blacs_env_type), POINTER                   :: context

#if defined(__parallel) && ! defined(__SCALAPACK)
      CPABORT("full matrices need ScaLAPACK for parallel run")
#endif

      CPASSERT(ASSOCIATED(matrix_struct))
      ALLOCATE (matrix)
      context => matrix_struct%context
      matrix%matrix_struct => matrix_struct
      CALL cp_fm_struct_retain(matrix%matrix_struct)
      last_cfm_id_nr = last_cfm_id_nr + 1
      matrix%id_nr = last_cfm_id_nr
      matrix%ref_count = 1
      matrix%print_count = 0

      nprow = context%num_pe(1)
      npcol = context%num_pe(2)
      NULLIFY (matrix%local_data)

      nrow_local = matrix_struct%local_leading_dimension
      ncol_local = MAX(1, matrix_struct%ncol_locals(context%mepos(2)))
      ALLOCATE (matrix%local_data(nrow_local, ncol_local))

      ! do not initialise created matrix as it is up to the user to do so
      !CALL zcopy(nrow_local*ncol_local, z_zero, 0, matrix%local_data, 1)

      IF (PRESENT(name)) THEN
         matrix%name = name
      ELSE
         matrix%name = 'full matrix'//cp_to_string(matrix%id_nr)
      END IF
   END SUBROUTINE cp_cfm_create

! **************************************************************************************************
!> \brief Retains a full matrix.
!> \param matrix the matrix to retain
! **************************************************************************************************
   SUBROUTINE cp_cfm_retain(matrix)
      TYPE(cp_cfm_type), POINTER                         :: matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_retain', routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(matrix))
      CPASSERT(matrix%ref_count > 0)
      matrix%ref_count = matrix%ref_count + 1
   END SUBROUTINE cp_cfm_retain

! **************************************************************************************************
!> \brief Releases a full matrix.
!> \param matrix the matrix to release
! **************************************************************************************************
   SUBROUTINE cp_cfm_release(matrix)
      TYPE(cp_cfm_type), POINTER                         :: matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_release', routineP = moduleN//':'//routineN

      IF (ASSOCIATED(matrix)) THEN
         CPASSERT(matrix%ref_count > 0)
         matrix%ref_count = matrix%ref_count - 1
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
!> \brief Set all elements of the full matrix to alpha. Besides, set all
!>        diagonal matrix elements to beta (if given explicitly).
!> \param matrix  matrix to initialise
!> \param alpha   value of off-diagonal matrix elements
!> \param beta    value of diagonal matrix elements (equal to alpha if absent)
!> \date    12.06.2001
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE cp_cfm_set_all(matrix, alpha, beta)
      TYPE(cp_cfm_type), POINTER                         :: matrix
      COMPLEX(kind=dp), INTENT(in)                       :: alpha
      COMPLEX(kind=dp), INTENT(in), OPTIONAL             :: beta

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_set_all', routineP = moduleN//':'//routineN

      INTEGER                                            :: irow_local, nrow_local
#if defined(__SCALAPACK)
      INTEGER                                            :: icol_local, ncol_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
#endif

      CALL zcopy(SIZE(matrix%local_data), alpha, 0, matrix%local_data(1, 1), 1)

      IF (PRESENT(beta)) THEN
#if defined(__SCALAPACK)
         CALL cp_cfm_get_info(matrix, nrow_local=nrow_local, ncol_local=ncol_local, &
                              row_indices=row_indices, col_indices=col_indices)

         icol_local = 1
         irow_local = 1

         DO WHILE (irow_local <= nrow_local .AND. icol_local <= ncol_local)
            IF (row_indices(irow_local) < col_indices(icol_local)) THEN
               irow_local = irow_local + 1
            ELSE IF (row_indices(irow_local) > col_indices(icol_local)) THEN
               icol_local = icol_local + 1
            ELSE
               matrix%local_data(irow_local, icol_local) = beta
               irow_local = irow_local + 1
               icol_local = icol_local + 1
            END IF
         END DO
#else
         nrow_local = MIN(matrix%matrix_struct%nrow_global, matrix%matrix_struct%ncol_global)

         DO irow_local = 1, nrow_local
            matrix%local_data(irow_local, irow_local) = beta
         END DO
#endif
      END IF

   END SUBROUTINE cp_cfm_set_all

! **************************************************************************************************
!> \brief Get the matrix element by its global index.
!> \param matrix      full matrix
!> \param irow_global global row index
!> \param icol_global global column index
!> \param alpha       matrix element
!> \par History
!>      , TCH, created
!>      always return the answer
! **************************************************************************************************
   SUBROUTINE cp_cfm_get_element(matrix, irow_global, icol_global, alpha)
      TYPE(cp_cfm_type), POINTER                         :: matrix
      INTEGER, INTENT(in)                                :: irow_global, icol_global
      COMPLEX(kind=dp), INTENT(out)                      :: alpha

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_get_element', &
                                     routineP = moduleN//':'//routineN

#if defined(__SCALAPACK)
      INTEGER                                            :: icol_local, ipcol, iprow, irow_local, &
                                                            mypcol, myprow, npcol, nprow
      INTEGER, DIMENSION(9)                              :: desca
      TYPE(cp_blacs_env_type), POINTER                   :: context
#endif

#if defined(__SCALAPACK)
      context => matrix%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      nprow = context%num_pe(1)
      npcol = context%num_pe(2)

      desca(:) = matrix%matrix_struct%descriptor(:)

      CALL infog2l(irow_global, icol_global, desca, nprow, npcol, myprow, mypcol, &
                   irow_local, icol_local, iprow, ipcol)

      IF ((iprow == myprow) .AND. (ipcol == mypcol)) THEN
         alpha = matrix%local_data(irow_local, icol_local)
         CALL cp_blacs_ZGEBS2D(context%group, 'All', ' ', 1, 1, alpha, 1)
      ELSE
         CALL cp_blacs_ZGEBR2D(context%group, 'All', ' ', 1, 1, alpha, 1, iprow, ipcol)
      END IF

#else
      alpha = matrix%local_data(irow_global, icol_global)
#endif

   END SUBROUTINE cp_cfm_get_element

! **************************************************************************************************
!> \brief Set the matrix element (irow_global,icol_global) of the full matrix to alpha.
!> \param matrix      full matrix
!> \param irow_global global row index
!> \param icol_global global column index
!> \param alpha       value of the matrix element
!> \date    12.06.2001
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE cp_cfm_set_element(matrix, irow_global, icol_global, alpha)
      TYPE(cp_cfm_type), POINTER                         :: matrix
      INTEGER, INTENT(in)                                :: irow_global, icol_global
      COMPLEX(kind=dp), INTENT(in)                       :: alpha

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_set_element', &
                                     routineP = moduleN//':'//routineN

#if defined(__SCALAPACK)
      INTEGER                                            :: icol_local, ipcol, iprow, irow_local, &
                                                            mypcol, myprow, npcol, nprow
      INTEGER, DIMENSION(9)                              :: desca
      TYPE(cp_blacs_env_type), POINTER                   :: context
#endif

#if defined(__SCALAPACK)
      context => matrix%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      nprow = context%num_pe(1)
      npcol = context%num_pe(2)

      desca(:) = matrix%matrix_struct%descriptor(:)

      CALL infog2l(irow_global, icol_global, desca, nprow, npcol, myprow, mypcol, &
                   irow_local, icol_local, iprow, ipcol)

      IF ((iprow == myprow) .AND. (ipcol == mypcol)) THEN
         matrix%local_data(irow_local, icol_local) = alpha
      END IF

#else
      matrix%local_data(irow_global, icol_global) = alpha
#endif

   END SUBROUTINE cp_cfm_set_element

! **************************************************************************************************
!> \brief Extract a sub-matrix from the full matrix:
!>        op(target_m)(1:n_rows,1:n_cols) = fm(start_row:start_row+n_rows,start_col:start_col+n_cols).
!>        Sub-matrix 'target_m' is replicated on each CPU. Using this call is expensive.
!> \param fm          full matrix you want to get the elements from
!> \param target_m    2-D array to store the extracted sub-matrix
!> \param start_row   global row index of the matrix element target_m(1,1) (defaults to 1)
!> \param start_col   global column index of the matrix element target_m(1,1) (defaults to 1)
!> \param n_rows      number of rows to extract (defaults to size(op(target_m),1))
!> \param n_cols      number of columns to extract (defaults to size(op(target_m),2))
!> \param transpose   indicates that the extracted sub-matrix target_m should be transposed:
!>                    op(target_m) = target_m^T if .TRUE.,
!>                    op(target_m) = target_m   if .FALSE. (defaults to false)
!> \par History
!>   * 04.2016 created borrowing from Fawzi's cp_fm_get_submatrix [Lianheng Tong]
!>   * 01.2018 drop innermost conditional branching [Sergey Chulkov]
!> \author Lianheng Tong
!> \note
!>      Optimized for full column updates. The matrix target_m is replicated and valid on all CPUs.
! **************************************************************************************************
   SUBROUTINE cp_cfm_get_submatrix(fm, target_m, start_row, start_col, n_rows, n_cols, transpose)
      TYPE(cp_cfm_type), POINTER                         :: fm
      COMPLEX(kind=dp), DIMENSION(:, :), INTENT(out)     :: target_m
      INTEGER, INTENT(in), OPTIONAL                      :: start_row, start_col, n_rows, n_cols
      LOGICAL, INTENT(in), OPTIONAL                      :: transpose

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_get_submatrix', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: local_data
      INTEGER :: end_col_global, end_col_local, end_row_global, end_row_local, handle, i, j, &
         ncol_global, ncol_local, nrow_global, nrow_local, start_col_global, start_col_local, &
         start_row_global, start_row_local, this_col
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: do_zero, tr_a
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(fm))
      CPASSERT(fm%ref_count .GT. 0)

      IF (SIZE(target_m) /= 0) THEN
#if defined(__SCALAPACK)
         do_zero = .TRUE.
#else
         do_zero = .FALSE.
#endif

         tr_a = .FALSE.
         IF (PRESENT(transpose)) tr_a = transpose

         ! find out the first and last global row/column indices
         start_row_global = 1
         start_col_global = 1
         IF (PRESENT(start_row)) start_row_global = start_row
         IF (PRESENT(start_col)) start_col_global = start_col

         IF (tr_a) THEN
            end_row_global = SIZE(target_m, 2)
            end_col_global = SIZE(target_m, 1)
         ELSE
            end_row_global = SIZE(target_m, 1)
            end_col_global = SIZE(target_m, 2)
         END IF
         IF (PRESENT(n_rows)) end_row_global = n_rows
         IF (PRESENT(n_cols)) end_col_global = n_cols

         end_row_global = end_row_global + start_row_global - 1
         end_col_global = end_col_global + start_col_global - 1

         CALL cp_cfm_get_info(matrix=fm, &
                              nrow_global=nrow_global, ncol_global=ncol_global, &
                              nrow_local=nrow_local, ncol_local=ncol_local, &
                              row_indices=row_indices, col_indices=col_indices)
         IF (end_row_global > nrow_global) THEN
            end_row_global = nrow_global
            do_zero = .TRUE.
         END IF
         IF (end_col_global > ncol_global) THEN
            end_col_global = ncol_global
            do_zero = .TRUE.
         END IF

         ! find out row/column indices of locally stored matrix elements that needs to be copied.
         ! Arrays row_indices and col_indices are assumed to be sorted in ascending order
         DO start_row_local = 1, nrow_local
            IF (row_indices(start_row_local) >= start_row_global) EXIT
         END DO

         DO end_row_local = start_row_local, nrow_local
            IF (row_indices(end_row_local) > end_row_global) EXIT
         END DO
         end_row_local = end_row_local - 1

         DO start_col_local = 1, ncol_local
            IF (col_indices(start_col_local) >= start_col_global) EXIT
         END DO

         DO end_col_local = start_col_local, ncol_local
            IF (col_indices(end_col_local) > end_col_global) EXIT
         END DO
         end_col_local = end_col_local - 1

         para_env => fm%matrix_struct%para_env
         local_data => fm%local_data

         ! wipe the content of the target matrix if:
         !  * the source matrix is distributed across a number of processes, or
         !  * not all elements of the target matrix will be assigned, e.g.
         !        when the target matrix is larger then the source matrix
         IF (do_zero) &
            CALL zcopy(SIZE(target_m), z_zero, 0, target_m(1, 1), 1)

         IF (tr_a) THEN
            DO j = start_col_local, end_col_local
               this_col = col_indices(j) - start_col_global + 1
               DO i = start_row_local, end_row_local
                  target_m(this_col, row_indices(i) - start_row_global + 1) = local_data(i, j)
               END DO
            END DO
         ELSE
            DO j = start_col_local, end_col_local
               this_col = col_indices(j) - start_col_global + 1
               DO i = start_row_local, end_row_local
                  target_m(row_indices(i) - start_row_global + 1, this_col) = local_data(i, j)
               END DO
            END DO
         END IF

         CALL mp_sum(target_m, para_env%group)
      END IF

      CALL timestop(handle)
   END SUBROUTINE cp_cfm_get_submatrix

! **************************************************************************************************
!> \brief Set a sub-matrix of the full matrix:
!>       matrix(start_row:start_row+n_rows,start_col:start_col+n_cols)
!>       = alpha*op(new_values)(1:n_rows,1:n_cols) +
!>         beta*matrix(start_row:start_row+n_rows,start_col:start_col+n_cols)
!> \param matrix      full to update
!> \param new_values  replicated 2-D array that holds new elements of the updated sub-matrix
!> \param start_row   global row index of the matrix element new_values(1,1) (defaults to 1)
!> \param start_col   global column index of the matrix element new_values(1,1) (defaults to 1)
!> \param n_rows      number of rows to update (defaults to size(op(new_values),1))
!> \param n_cols      number of columns to update (defaults to size(op(new_values),2))
!> \param alpha       scale factor for the new values (defaults to (1.0,0.0))
!> \param beta        scale factor for the old values (defaults to (0.0,0.0))
!> \param transpose   indicates that the matrix new_values should be transposed:
!>                    op(new_values) = new_values^T if .TRUE.,
!>                    op(new_values) = new_values   if .FALSE. (defaults to false)
!> \par History
!>   * 04.2016 created borrowing from Fawzi's cp_fm_set_submatrix [Lianheng Tong]
!>   * 01.2018 drop innermost conditional branching [Sergey Chulkov]
!> \author Lianheng Tong
!> \note
!>      Optimized for alpha=(1.0,0.0), beta=(0.0,0.0)
!>      All matrix elements 'new_values' need to be valid on all CPUs
! **************************************************************************************************
   SUBROUTINE cp_cfm_set_submatrix(matrix, new_values, start_row, &
                                   start_col, n_rows, n_cols, alpha, beta, transpose)
      TYPE(cp_cfm_type), POINTER                         :: matrix
      COMPLEX(kind=dp), DIMENSION(:, :), INTENT(in)      :: new_values
      INTEGER, INTENT(in), OPTIONAL                      :: start_row, start_col, n_rows, n_cols
      COMPLEX(kind=dp), INTENT(in), OPTIONAL             :: alpha, beta
      LOGICAL, INTENT(in), OPTIONAL                      :: transpose

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_set_submatrix', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp)                                   :: al, be
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: local_data
      INTEGER :: end_col_global, end_col_local, end_row_global, end_row_local, handle, i, j, &
         ncol_global, ncol_local, nrow_global, nrow_local, start_col_global, start_col_local, &
         start_row_global, start_row_local, this_col
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: tr_a

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(matrix))
      CPASSERT(matrix%ref_count > 0)

      al = z_one
      be = z_zero
      IF (PRESENT(alpha)) al = alpha
      IF (PRESENT(beta)) be = beta

      ! find out the first and last global row/column indices
      start_row_global = 1
      start_col_global = 1
      IF (PRESENT(start_row)) start_row_global = start_row
      IF (PRESENT(start_col)) start_col_global = start_col

      tr_a = .FALSE.
      IF (PRESENT(transpose)) tr_a = transpose

      IF (tr_a) THEN
         end_row_global = SIZE(new_values, 2)
         end_col_global = SIZE(new_values, 1)
      ELSE
         end_row_global = SIZE(new_values, 1)
         end_col_global = SIZE(new_values, 2)
      END IF
      IF (PRESENT(n_rows)) end_row_global = n_rows
      IF (PRESENT(n_cols)) end_col_global = n_cols

      end_row_global = end_row_global + start_row_global - 1
      end_col_global = end_col_global + start_col_global - 1

      CALL cp_cfm_get_info(matrix=matrix, &
                           nrow_global=nrow_global, ncol_global=ncol_global, &
                           nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices)
      IF (end_row_global > nrow_global) end_row_global = nrow_global
      IF (end_col_global > ncol_global) end_col_global = ncol_global

      ! find out row/column indices of locally stored matrix elements that needs to be set.
      ! Arrays row_indices and col_indices are assumed to be sorted in ascending order
      DO start_row_local = 1, nrow_local
         IF (row_indices(start_row_local) >= start_row_global) EXIT
      END DO

      DO end_row_local = start_row_local, nrow_local
         IF (row_indices(end_row_local) > end_row_global) EXIT
      END DO
      end_row_local = end_row_local - 1

      DO start_col_local = 1, ncol_local
         IF (col_indices(start_col_local) >= start_col_global) EXIT
      END DO

      DO end_col_local = start_col_local, ncol_local
         IF (col_indices(end_col_local) > end_col_global) EXIT
      END DO
      end_col_local = end_col_local - 1

      local_data => matrix%local_data

      IF (al == z_one .AND. be == z_zero) THEN
         IF (tr_a) THEN
            DO j = start_col_local, end_col_local
               this_col = col_indices(j) - start_col_global + 1
               DO i = start_row_local, end_row_local
                  local_data(i, j) = new_values(this_col, row_indices(i) - start_row_global + 1)
               END DO
            END DO
         ELSE
            DO j = start_col_local, end_col_local
               this_col = col_indices(j) - start_col_global + 1
               DO i = start_row_local, end_row_local
                  local_data(i, j) = new_values(row_indices(i) - start_row_global + 1, this_col)
               END DO
            END DO
         END IF
      ELSE
         IF (tr_a) THEN
            DO j = start_col_local, end_col_local
               this_col = col_indices(j) - start_col_global + 1
               DO i = start_row_local, end_row_local
                  local_data(i, j) = al*new_values(this_col, row_indices(i) - start_row_global + 1) + &
                                     be*local_data(i, j)
               END DO
            END DO
         ELSE
            DO j = start_col_local, end_col_local
               this_col = col_indices(j) - start_col_global + 1
               DO i = start_row_local, end_row_local
                  local_data(i, j) = al*new_values(row_indices(i) - start_row_global + 1, this_col) + &
                                     be*local_data(i, j)
               END DO
            END DO
         END IF
      END IF

      CALL timestop(handle)
   END SUBROUTINE cp_cfm_set_submatrix

! **************************************************************************************************
!> \brief Returns information about a full matrix.
!> \param matrix        matrix
!> \param name          name of the matrix
!> \param nrow_global   total number of rows
!> \param ncol_global   total number of columns
!> \param nrow_block    number of rows per ScaLAPACK block
!> \param ncol_block    number of columns per ScaLAPACK block
!> \param nrow_local    number of locally stored rows
!> \param ncol_local    number of locally stored columns
!> \param row_indices   global indices of locally stored rows
!> \param col_indices   global indices of locally stored columns
!> \param local_data    locally stored matrix elements
!> \param context       BLACS context
!> \param matrix_struct matrix structure
!> \param para_env      parallel environment
!> \date    12.06.2001
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE cp_cfm_get_info(matrix, name, nrow_global, ncol_global, &
                              nrow_block, ncol_block, nrow_local, ncol_local, &
                              row_indices, col_indices, local_data, context, &
                              matrix_struct, para_env)
      TYPE(cp_cfm_type), POINTER        :: matrix
      CHARACTER(len=*), OPTIONAL, INTENT(OUT) :: name
      INTEGER, OPTIONAL, INTENT(OUT)          :: ncol_block, ncol_global, &
                                                 nrow_block, nrow_global, &
                                                 nrow_local, ncol_local
      INTEGER, OPTIONAL, DIMENSION(:), POINTER   :: row_indices, col_indices
      TYPE(cp_para_env_type), POINTER, OPTIONAL :: para_env
      TYPE(cp_blacs_env_type), POINTER, OPTIONAL :: context
      TYPE(cp_fm_struct_type), POINTER, OPTIONAL :: matrix_struct
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER, OPTIONAL :: local_data

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
#if defined(__SCALAPACK)
            DO i = 1, SIZE(row_indices)
               row_indices(i) = &
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
#if defined(__SCALAPACK)
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
!> \brief Copy content of a full matrix into another full matrix of the same size.
!> \param source      source matrix
!> \param destination destination matrix
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
         CALL zcopy(SIZE(source%local_data), source%local_data(1, 1), 1, destination%local_data(1, 1), 1)
      ELSE
         IF (source%matrix_struct%nrow_global /= destination%matrix_struct%nrow_global) &
            CPABORT("cannot copy between full matrixes of differen sizes")
         IF (source%matrix_struct%ncol_global /= destination%matrix_struct%ncol_global) &
            CPABORT("cannot copy between full matrixes of differen sizes")
#if defined(__SCALAPACK)
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
!> \brief Copy a number of sequential columns of a full matrix into another full matrix.
!> \param msource      source matrix
!> \param mtarget      destination matrix
!> \param ncol         number of columns to copy
!> \param source_start global index of the first column to copy within the source matrix
!> \param target_start global index of the first column to copy within the destination matrix
! **************************************************************************************************
   SUBROUTINE cp_cfm_to_cfm_columns(msource, mtarget, ncol, source_start, &
                                    target_start)

      TYPE(cp_cfm_type), POINTER                         :: msource, mtarget
      INTEGER, INTENT(IN)                                :: ncol
      INTEGER, INTENT(IN), OPTIONAL                      :: source_start, target_start

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_to_cfm_columns', &
                                     routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: a, b
      INTEGER                                            :: handle, n, ss, ts
#if defined(__SCALAPACK)
      INTEGER                                            :: i
      INTEGER, DIMENSION(9)                              :: desca, descb
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
      DO i = 0, ncol - 1
         CALL pzcopy(n, a(1, 1), 1, ss + i, desca, 1, b(1, 1), 1, ts + i, descb, 1)
      END DO
#else
      CALL zcopy(ncol*n, a(1, ss), 1, b(1, ts), 1)
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_to_cfm_columns

! **************************************************************************************************
!> \brief Copy just a triangular matrix.
!> \param msource source matrix
!> \param mtarget target matrix
!> \param uplo    'U' for upper triangular, 'L' for lower triangular
! **************************************************************************************************
   SUBROUTINE cp_cfm_to_cfm_triangular(msource, mtarget, uplo)
      TYPE(cp_cfm_type), POINTER                         :: msource, mtarget
      CHARACTER(len=*), INTENT(IN)                       :: uplo

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_to_cfm_triangular', &
                                     routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: aa, bb
      INTEGER                                            :: handle, ncol, nrow
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                              :: desca, descb
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
!> \brief Copy real and imaginary parts of a complex full matrix into
!>        separate real-value full matrices.
!> \param msource    complex matrix
!> \param mtargetr   (optional) real part of the source matrix
!> \param mtargeti   (optional) imaginary part of the source matrix
!> \note
!>        Matrix structures are assumed to be equivalent.
! **************************************************************************************************
   SUBROUTINE cp_cfm_to_fm(msource, mtargetr, mtargeti)

      TYPE(cp_cfm_type), POINTER                         :: msource
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: mtargetr, mtargeti

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_to_fm', routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: zmat
      INTEGER                                            :: handle
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: imat, rmat

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
         rmat = REAL(zmat, kind=dp)
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
         imat = REAL(AIMAG(zmat), kind=dp)
      ELSE
         NULLIFY (imat)
      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_to_fm

! **************************************************************************************************
!> \brief Construct a complex full matrix by taking its real and imaginary parts from
!>        two separate real-value full matrices.
!> \param msourcer   (optional) real part of the complex matrix (defaults to 0.0)
!> \param msourcei   (optional) imaginary part of the complex matrix (defaults to 0.0)
!> \param mtarget    resulting complex matrix
!> \note
!>        Matrix structures are assumed to be equivalent.
! **************************************************************************************************
   SUBROUTINE cp_fm_to_cfm(msourcer, msourcei, mtarget)
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: msourcer, msourcei
      TYPE(cp_cfm_type), POINTER                         :: mtarget

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_to_cfm', routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: zmat
      INTEGER                                            :: handle, mode
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: imat, rmat

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
         mode = mode + 1
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
         mode = mode + 2
      ELSE
         NULLIFY (imat)
      END IF
      ! copy local data
      SELECT CASE (mode)
      CASE (0)
         zmat(:, :) = z_zero
      CASE (1)
         zmat(:, :) = CMPLX(rmat(:, :), 0.0_dp, kind=dp)
      CASE (2)
         zmat(:, :) = CMPLX(0.0_dp, imat(:, :), kind=dp)
      CASE (3)
         zmat(:, :) = CMPLX(rmat(:, :), imat(:, :), kind=dp)
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE cp_fm_to_cfm

! **************************************************************************************************
!> \brief Initiate the copy operation: get distribution data, post MPI isend and irecvs.
!> \param source       input complex-valued fm matrix
!> \param destination  output complex-valued fm matrix
!> \param para_env     parallel environment corresponding to the BLACS env that covers all parts
!>                     of the input and output matrices
!> \param info         all of the data that will be needed to complete the copy operation
!> \note a slightly modified version of the subroutine cp_fm_start_copy_general() that uses
!>       allocatable arrays instead of pointers wherever possible.
! **************************************************************************************************
   SUBROUTINE cp_cfm_start_copy_general(source, destination, para_env, info)
      TYPE(cp_cfm_type), POINTER                         :: source, destination
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(copy_cfm_info_type), INTENT(out)              :: info

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_start_copy_general', &
         routineP = moduleN//':'//routineN

      INTEGER :: dest_p_i, dest_q_j, global_comm, global_rank, global_size, handle, i, j, k, &
         mpi_rank, ncol_block_dest, ncol_block_src, ncol_local_recv, ncol_local_send, ncols, &
         nrow_block_dest, nrow_block_src, nrow_local_recv, nrow_local_send, nrows, p, q, &
         recv_rank, recv_size, send_rank, send_size
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: all_ranks, dest2global, dest_p, dest_q, &
                                                            recv_count, send_count, send_disp, &
                                                            source2global, src_p, src_q
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: dest_blacs2mpi
      INTEGER, DIMENSION(2)                              :: dest_block, dest_block_tmp, dest_num_pe, &
                                                            src_block, src_block_tmp, src_num_pe
      INTEGER, DIMENSION(6)                              :: recv_req, send_req
      INTEGER, DIMENSION(:), POINTER                     :: recv_col_indices, recv_row_indices, &
                                                            send_col_indices, send_row_indices
      TYPE(cp_fm_struct_type), POINTER                   :: recv_dist, send_dist

      CALL timeset(routineN, handle)

      IF (.NOT. cp2k_is_parallel) THEN
         ! Just copy all of the matrix data into a 'send buffer', to be unpacked later
         nrow_local_send = SIZE(source%local_data, 1)
         ncol_local_send = SIZE(source%local_data, 2)
         ALLOCATE (info%send_buf(nrow_local_send*ncol_local_send))
         k = 0
         DO j = 1, ncol_local_send
            DO i = 1, nrow_local_send
               k = k + 1
               info%send_buf(k) = source%local_data(i, j)
            END DO
         END DO
      ELSE
         NULLIFY (recv_dist, send_dist)
         NULLIFY (recv_col_indices, recv_row_indices, send_col_indices, send_row_indices)

         ! The 'global' communicator contains both the source and destination decompositions
         global_comm = para_env%group
         global_size = para_env%num_pe
         global_rank = para_env%mepos

         ! The source/send decomposition and destination/recv decompositions may only exist on
         ! on a subset of the processes involved in the communication
         ! Check if the source and/or destination arguments are .not. ASSOCIATED():
         ! if so, skip the send / recv parts (since these processes do not participate in the sending/receiving distribution)
         IF (ASSOCIATED(destination)) THEN
            recv_dist => destination%matrix_struct
            recv_rank = recv_dist%para_env%mepos
         ELSE
            recv_rank = mp_proc_null
         END IF

         IF (ASSOCIATED(source)) THEN
            send_dist => source%matrix_struct
            send_rank = send_dist%para_env%mepos
         ELSE
            send_rank = mp_proc_null
         END IF

         ! Map the rank in the source/dest communicator to the global rank
         ALLOCATE (all_ranks(0:global_size - 1))

         CALL mp_allgather(send_rank, all_ranks, global_comm)
         IF (ASSOCIATED(destination)) THEN
            ALLOCATE (source2global(0:COUNT(all_ranks .NE. mp_proc_null) - 1))
            DO i = 0, global_size - 1
               IF (all_ranks(i) .NE. mp_proc_null) THEN
                  source2global(all_ranks(i)) = i
               END IF
            END DO
         END IF

         CALL mp_allgather(recv_rank, all_ranks, global_comm)
         IF (ASSOCIATED(source)) THEN
            ALLOCATE (dest2global(0:COUNT(all_ranks .NE. mp_proc_null) - 1))
            DO i = 0, global_size - 1
               IF (all_ranks(i) .NE. mp_proc_null) THEN
                  dest2global(all_ranks(i)) = i
               END IF
            END DO
         END IF
         DEALLOCATE (all_ranks)

         ! Some data from the two decompositions will be needed by all processes in the global group :
         ! process grid shape, block size, and the BLACS-to-MPI mapping

         ! The global root process will receive the data (from the root process in each decomposition)
         send_req(:) = mp_request_null
         IF (global_rank == 0) THEN
            recv_req(:) = mp_request_null
            CALL mp_irecv(src_block, mp_any_source, global_comm, recv_req(1), tag=src_tag)
            CALL mp_irecv(dest_block, mp_any_source, global_comm, recv_req(2), tag=dest_tag)
            CALL mp_irecv(src_num_pe, mp_any_source, global_comm, recv_req(3), tag=src_tag)
            CALL mp_irecv(dest_num_pe, mp_any_source, global_comm, recv_req(4), tag=dest_tag)
         END IF

         IF (ASSOCIATED(source)) THEN
            IF ((send_rank == 0)) THEN
               ! need to use separate buffers here in case this is actually global rank 0
               src_block_tmp = (/send_dist%nrow_block, send_dist%ncol_block/)
               CALL mp_isend(src_block_tmp, 0, global_comm, send_req(1), tag=src_tag)
               CALL mp_isend(send_dist%context%num_pe, 0, global_comm, send_req(2), tag=src_tag)
            END IF
         END IF

         IF (ASSOCIATED(destination)) THEN
            IF ((recv_rank == 0)) THEN
               dest_block_tmp = (/recv_dist%nrow_block, recv_dist%ncol_block/)
               CALL mp_isend(dest_block_tmp, 0, global_comm, send_req(3), tag=dest_tag)
               CALL mp_isend(recv_dist%context%num_pe, 0, global_comm, send_req(4), tag=dest_tag)
            END IF
         END IF

         IF (global_rank == 0) THEN
            CALL mp_waitall(recv_req(1:4))
            ! Now we know the process decomposition, we can allocate the arrays to hold the blacs2mpi mapping
            ALLOCATE (info%src_blacs2mpi(0:src_num_pe(1) - 1, 0:src_num_pe(2) - 1), &
                      dest_blacs2mpi(0:dest_num_pe(1) - 1, 0:dest_num_pe(2) - 1))
            CALL mp_irecv(info%src_blacs2mpi, mp_any_source, global_comm, recv_req(5), tag=src_tag)
            CALL mp_irecv(dest_blacs2mpi, mp_any_source, global_comm, recv_req(6), tag=dest_tag)
         END IF

         IF (ASSOCIATED(source)) THEN
            IF ((send_rank == 0)) THEN
               CALL mp_isend(send_dist%context%blacs2mpi(:, :), 0, global_comm, send_req(5), tag=src_tag)
            END IF
         END IF

         IF (ASSOCIATED(destination)) THEN
            IF ((recv_rank == 0)) THEN
               CALL mp_isend(recv_dist%context%blacs2mpi(:, :), 0, global_comm, send_req(6), tag=dest_tag)
            END IF
         END IF

         IF (global_rank == 0) THEN
            CALL mp_waitall(recv_req(5:6))
         END IF

         ! Finally, broadcast the data to all processes in the global communicator
         CALL mp_bcast(src_block, 0, global_comm)
         CALL mp_bcast(dest_block, 0, global_comm)
         CALL mp_bcast(src_num_pe, 0, global_comm)
         CALL mp_bcast(dest_num_pe, 0, global_comm)
         info%src_num_pe(1:2) = src_num_pe(1:2)
         info%nblock_src(1:2) = src_block(1:2)
         IF (global_rank /= 0) THEN
            ALLOCATE (info%src_blacs2mpi(0:src_num_pe(1) - 1, 0:src_num_pe(2) - 1), &
                      dest_blacs2mpi(0:dest_num_pe(1) - 1, 0:dest_num_pe(2) - 1))
         END IF
         CALL mp_bcast(info%src_blacs2mpi, 0, global_comm)
         CALL mp_bcast(dest_blacs2mpi, 0, global_comm)

         recv_size = dest_num_pe(1)*dest_num_pe(2)
         send_size = src_num_pe(1)*src_num_pe(2)
         info%send_size = send_size
         CALL mp_waitall(send_req(:))

         ! Setup is now complete, we can start the actual communication here.
         ! The order implemented here is:
         !  DEST_1
         !      compute recv sizes
         !      call mp_irecv
         !  SRC_1
         !      compute send sizes
         !      pack send buffers
         !      call mp_isend
         !  DEST_2
         !      wait for the recvs and unpack buffers (this part eventually will go into another routine to allow comms to run concurrently)
         !  SRC_2
         !      wait for the sends

         ! DEST_1
         IF (ASSOCIATED(destination)) THEN
            CALL cp_fm_struct_get(recv_dist, row_indices=recv_row_indices, &
                                  col_indices=recv_col_indices)
            info%recv_col_indices => recv_col_indices
            info%recv_row_indices => recv_row_indices
            nrow_block_src = src_block(1)
            ncol_block_src = src_block(2)
            ALLOCATE (recv_count(0:send_size - 1), info%recv_disp(0:send_size - 1), info%recv_request(0:send_size - 1))

            ! Determine the recv counts, allocate the receive buffers, call mpi_irecv for all the non-zero sized receives
            nrow_local_recv = recv_dist%nrow_locals(recv_dist%context%mepos(1))
            ncol_local_recv = recv_dist%ncol_locals(recv_dist%context%mepos(2))
            info%nlocal_recv(1) = nrow_local_recv
            info%nlocal_recv(2) = ncol_local_recv
            ! Initialise src_p, src_q arrays (sized using number of rows/cols in the receiving distribution)
            ALLOCATE (src_p(nrow_local_recv), src_q(ncol_local_recv))
            DO i = 1, nrow_local_recv
               ! For each local row we will receive, we look up its global row (in recv_row_indices),
               ! then work out which row block it comes from, and which process row that row block comes from.
               src_p(i) = MOD(((recv_row_indices(i) - 1)/nrow_block_src), src_num_pe(1))
            END DO
            DO j = 1, ncol_local_recv
               ! Similarly for the columns
               src_q(j) = MOD(((recv_col_indices(j) - 1)/ncol_block_src), src_num_pe(2))
            END DO
            ! src_p/q now contains the process row/column ID that will send data to that row/column

            DO q = 0, src_num_pe(2) - 1
               ncols = COUNT(src_q .EQ. q)
               DO p = 0, src_num_pe(1) - 1
                  nrows = COUNT(src_p .EQ. p)
                  ! Use the send_dist here as we are looking up the processes where the data comes from
                  recv_count(info%src_blacs2mpi(p, q)) = nrows*ncols
               END DO
            END DO
            DEALLOCATE (src_p, src_q)

            ! Use one long buffer (and displacements into that buffer)
            !     this prevents the need for a rectangular array where not all elements will be populated
            ALLOCATE (info%recv_buf(SUM(recv_count(:))))
            info%recv_disp(0) = 0
            DO i = 1, send_size - 1
               info%recv_disp(i) = info%recv_disp(i - 1) + recv_count(i - 1)
            END DO

            ! Issue receive calls on ranks which expect data
            DO k = 0, send_size - 1
               IF (recv_count(k) .GT. 0) THEN
                  CALL mp_irecv(info%recv_buf(info%recv_disp(k) + 1:info%recv_disp(k) + recv_count(k)), &
                                source2global(k), global_comm, info%recv_request(k))
               ELSE
                  info%recv_request(k) = mp_request_null
               END IF
            END DO
            DEALLOCATE (source2global)
         END IF ! ASSOCIATED(destination)

         ! SRC_1
         IF (ASSOCIATED(source)) THEN
            CALL cp_fm_struct_get(send_dist, row_indices=send_row_indices, &
                                  col_indices=send_col_indices)
            nrow_block_dest = dest_block(1)
            ncol_block_dest = dest_block(2)
            ALLOCATE (send_count(0:recv_size - 1), send_disp(0:recv_size - 1), info%send_request(0:recv_size - 1))

            ! Determine the send counts, allocate the send buffers
            nrow_local_send = send_dist%nrow_locals(send_dist%context%mepos(1))
            ncol_local_send = send_dist%ncol_locals(send_dist%context%mepos(2))

            ! Initialise dest_p, dest_q arrays (sized nrow_local, ncol_local)
            !   i.e. number of rows,cols in the sending distribution
            ALLOCATE (dest_p(nrow_local_send), dest_q(ncol_local_send))

            DO i = 1, nrow_local_send
               ! Use the send_dist%row_indices() here (we are looping over the local rows we will send)
               dest_p(i) = MOD(((send_row_indices(i) - 1)/nrow_block_dest), dest_num_pe(1))
            END DO
            DO j = 1, ncol_local_send
               dest_q(j) = MOD(((send_col_indices(j) - 1)/ncol_block_dest), dest_num_pe(2))
            END DO
            ! dest_p/q now contain the process row/column ID that will receive data from that row/column

            DO q = 0, dest_num_pe(2) - 1
               ncols = COUNT(dest_q .EQ. q)
               DO p = 0, dest_num_pe(1) - 1
                  nrows = COUNT(dest_p .EQ. p)
                  send_count(dest_blacs2mpi(p, q)) = nrows*ncols
               END DO
            END DO
            DEALLOCATE (dest_p, dest_q)

            ! Allocate the send buffer using send_count -- and calculate the offset into the buffer for each process
            ALLOCATE (info%send_buf(SUM(send_count(:))))
            send_disp(0) = 0
            DO k = 1, recv_size - 1
               send_disp(k) = send_disp(k - 1) + send_count(k - 1)
            END DO

            ! Loop over the smat, pack the send buffers
            send_count(:) = 0
            DO j = 1, ncol_local_send
               ! Use send_col_indices and row_indices here, as we are looking up the global row/column number of local rows.
               dest_q_j = MOD(((send_col_indices(j) - 1)/ncol_block_dest), dest_num_pe(2))
               DO i = 1, nrow_local_send
                  dest_p_i = MOD(((send_row_indices(i) - 1)/nrow_block_dest), dest_num_pe(1))
                  mpi_rank = dest_blacs2mpi(dest_p_i, dest_q_j)
                  send_count(mpi_rank) = send_count(mpi_rank) + 1
                  info%send_buf(send_disp(mpi_rank) + send_count(mpi_rank)) = source%local_data(i, j)
               END DO
            END DO

            ! For each non-zero send_count, call mpi_isend
            DO k = 0, recv_size - 1
               IF (send_count(k) .GT. 0) THEN
                  CALL mp_isend(info%send_buf(send_disp(k) + 1:send_disp(k) + send_count(k)), &
                                dest2global(k), global_comm, info%send_request(k))
               ELSE
                  info%send_request(k) = mp_request_null
               END IF
            END DO
            DEALLOCATE (send_count, send_disp, dest2global)
         END IF ! ASSOCIATED(source)
         DEALLOCATE (dest_blacs2mpi)

      END IF !IF (.NOT. cp2k_is_parallel)

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_start_copy_general

! **************************************************************************************************
!> \brief Complete the copy operation: wait for comms, unpack, clean up MPI state.
!> \param destination  output cfm matrix
!> \param info         all of the data that will be needed to complete the copy operation
!> \note a slightly modified version of the subroutine cp_fm_finish_copy_general() that uses
!>       allocatable arrays instead of pointers wherever possible.
! **************************************************************************************************
   SUBROUTINE cp_cfm_finish_copy_general(destination, info)
      TYPE(cp_cfm_type), POINTER                         :: destination
      TYPE(copy_cfm_info_type), INTENT(inout)            :: info

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_finish_copy_general', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, k, mpi_rank, ni, nj, &
                                                            src_q_j
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: recv_count, src_p_i
      INTEGER, DIMENSION(:), POINTER                     :: recv_col_indices, recv_row_indices

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(destination))

      IF (.NOT. cp2k_is_parallel) THEN
         ! Now unpack the data from the 'send buffer'
         k = 0
         DO j = 1, SIZE(destination%local_data, 2)
            DO i = 1, SIZE(destination%local_data, 1)
               k = k + 1
               destination%local_data(i, j) = info%send_buf(k)
            END DO
         END DO
         DEALLOCATE (info%send_buf)
      ELSE
         ! Set up local variables ...
         recv_col_indices => info%recv_col_indices
         recv_row_indices => info%recv_row_indices

         ! ... use the local variables to do the work
         ! DEST_2
         CALL mp_waitall(info%recv_request(:))

         nj = info%nlocal_recv(2)
         ni = info%nlocal_recv(1)
         ALLOCATE (recv_count(0:info%send_size - 1), src_p_i(ni))
         ! Loop over the rmat, filling it in with data from the recv buffers
         ! (here the block sizes, num_pes refer to the distribution of the source matrix)
         recv_count(:) = 0
         DO i = 1, ni
            src_p_i(i) = MOD(((recv_row_indices(i) - 1)/info%nblock_src(1)), info%src_num_pe(1))
         END DO

         DO j = 1, nj
            src_q_j = MOD(((recv_col_indices(j) - 1)/info%nblock_src(2)), info%src_num_pe(2))
            DO i = 1, ni
               mpi_rank = info%src_blacs2mpi(src_p_i(i), src_q_j)
               recv_count(mpi_rank) = recv_count(mpi_rank) + 1
               destination%local_data(i, j) = info%recv_buf(info%recv_disp(mpi_rank) + recv_count(mpi_rank))
            END DO
         END DO

         DEALLOCATE (recv_count, src_p_i)
         ! Invalidate the stored state
         NULLIFY (info%recv_col_indices, info%recv_row_indices)
         DEALLOCATE (info%recv_disp, info%recv_request, info%recv_buf, info%src_blacs2mpi)
      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_finish_copy_general

! **************************************************************************************************
!> \brief Complete the copy operation: wait for comms clean up MPI state.
!> \param source  cfm matrix
!> \param info    all of the data that will be needed to complete the copy operation
!> \note a slightly modified version of the subroutine cp_fm_cleanup_copy_general() that uses
!>       allocatable arrays instead of pointers wherever possible.
! **************************************************************************************************
   SUBROUTINE cp_cfm_cleanup_copy_general(source, info)
      TYPE(cp_cfm_type), POINTER                         :: source
      TYPE(copy_cfm_info_type), INTENT(inout)            :: info

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_cleanup_copy_general', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(source))
      MARK_USED(source)

      IF (.NOT. cp2k_is_parallel) THEN
         ! Don't do anything - no MPI state for the serial case
      ELSE
         ! SRC_2
         ! If this process is also in the destination decomposition, this deallocate
         ! Was already done in cp_fm_finish_copy_general
         IF (ALLOCATED(info%src_blacs2mpi)) DEALLOCATE (info%src_blacs2mpi)
         CALL mp_waitall(info%send_request(:))
         DEALLOCATE (info%send_request, info%send_buf)
      END IF

      CALL timestop(handle)
   END SUBROUTINE cp_cfm_cleanup_copy_general
END MODULE cp_cfm_types
