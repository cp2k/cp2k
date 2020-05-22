!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief represent a full matrix distributed on many processors
!> \par History
!>      3) separated structure object, removed globenv, renamed to full matrix
!>         many changes (fawzi 08.2002)
!> \author Matthias Krack (22.05.2001)
! **************************************************************************************************
MODULE cp_fm_types
   USE cp_blacs_calls,                  ONLY: cp_blacs_dgebr2d,&
                                              cp_blacs_dgebs2d,&
                                              cp_blacs_gridexit,&
                                              cp_blacs_gridinfo,&
                                              cp_blacs_gridinit
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent,&
                                              cp_fm_struct_get,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_retain,&
                                              cp_fm_struct_type,&
                                              cp_fm_struct_write_info
   USE cp_log_handling,                 ONLY: cp_to_string
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: dp,&
                                              sp
   USE message_passing,                 ONLY: &
        cp2k_is_parallel, mp_allgather, mp_any_source, mp_bcast, mp_irecv, mp_isend, mp_max, &
        mp_proc_null, mp_recv, mp_request_null, mp_send, mp_sum, mp_waitall
   USE parallel_rng_types,              ONLY: UNIFORM,&
                                              rng_stream_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_types'
   LOGICAL, PARAMETER          :: debug_this_module = .TRUE.
   INTEGER, PARAMETER :: src_tag = 3, dest_tag = 5, send_tag = 7, recv_tag = 11

   INTEGER, SAVE :: last_fm_id_nr = 0
   INTEGER, PRIVATE :: mm_type = 1

   PUBLIC :: cp_fm_type, &
             cp_fm_p_type, copy_info_type

   PUBLIC :: cp_fm_add_to_element, &
             cp_fm_create, &
             cp_fm_retain, &
             cp_fm_release, &
             cp_fm_get_info, &
             cp_fm_set_element, &
             cp_fm_get_element, &
             cp_fm_get_diag, & ! get diagonal
             cp_fm_set_all, & ! set all elements and diagonal
             cp_fm_set_submatrix, & ! set a submatrix to given values
             cp_fm_get_submatrix, & ! get a submatrix of given values
             cp_fm_init_random, &
             cp_fm_maxabsval, & ! find the maximum absolute value
             cp_fm_maxabsrownorm, & ! find the maximum of the sum of the abs of the elements of a row
             cp_fm_to_fm, & ! copy (parts of) a fm to a fm
             cp_fm_vectorsnorm, & ! compute the norm of the column-vectors
             cp_fm_vectorssum, & ! compute the sum of all elements of the column-vectors
             cp_fm_to_fm_submat, & ! copy (parts of) a fm to a fm
             cp_fm_to_fm_triangular, &
             cp_fm_copy_general, &
             cp_fm_start_copy_general, &
             cp_fm_finish_copy_general, &
             cp_fm_cleanup_copy_general, &
             cp_fm_write_unformatted, & ! writes a full matrix to an open unit
             cp_fm_write_formatted, & ! writes a full matrix to an open unit
             cp_fm_read_unformatted, & ! reads a full matrix from an open unit
             cp_fm_setup, & ! allows to set flags for fms
             cp_fm_get_mm_type, &
             cp_fm_write_info

   PUBLIC :: cp_fm_indxg2p, &
             cp_fm_indxg2l, &
             cp_fm_indxl2g

   INTERFACE cp_fm_to_fm
      MODULE PROCEDURE cp_fm_to_fm_matrix, & ! a full matrix
         cp_fm_to_fm_columns ! just a number of columns
   END INTERFACE

! **************************************************************************************************
!> \brief represent a full matrix
!> \param name the name of the matrix, used for printing
!> \param id_nr identification number that makes this matrix unique
!> \param ref_cont reference count (see doc/ReferenceCounting.html)
!> \param print_count sequential number to make the printed files unique
!> \param matrix_struct structure of this matrix
!> \param local_data array with the data of the matrix (its contents
!>        depend on the matrix type used: in parallel runs it will be
!>        in scalapack format, in sequential, it will simply contain
!>        the matrix)
!> \par History
!>      08.2002 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   TYPE cp_fm_type
!    PRIVATE
      CHARACTER(LEN=60) :: name
      INTEGER :: id_nr, ref_count, print_count
      LOGICAL :: use_sp
      TYPE(cp_fm_struct_type), POINTER :: matrix_struct
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: local_data
      REAL(KIND=sp), DIMENSION(:, :), POINTER :: local_data_sp
   END TYPE cp_fm_type

! **************************************************************************************************
!> \brief just to build arrays of pointers to matrices
!> \param matrix the pointer to the matrix
!> \par History
!>      08.2002 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   TYPE cp_fm_p_type
      TYPE(cp_fm_type), POINTER :: matrix
   END TYPE cp_fm_p_type

! **************************************************************************************************
!> \brief Stores the state of a copy between cp_fm_start_copy_general
!>        and cp_fm_finish_copy_general
!> \par History
!>      Jan 2017  [Mark T]
! **************************************************************************************************
   TYPE copy_info_type
      INTEGER :: send_size
      INTEGER, DIMENSION(2) :: nlocal_recv, nblock_src, src_num_pe ! 1->row  2->col
      INTEGER, DIMENSION(:), POINTER       :: send_request, recv_request, recv_disp
      INTEGER, DIMENSION(:), POINTER       :: recv_col_indices, recv_row_indices
      INTEGER, DIMENSION(:, :), POINTER    :: src_blacs2mpi
      REAL(KIND=dp), DIMENSION(:), POINTER :: recv_buf, send_buf
   END TYPE copy_info_type

CONTAINS

! **************************************************************************************************
!> \brief creates a new full matrix with the given structure
!> \param matrix the matrix to be created
!> \param matrix_struct the structure of matrix
!> \param name ...
!> \param use_sp ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      preferred allocation routine
! **************************************************************************************************
   SUBROUTINE cp_fm_create(matrix, matrix_struct, name, use_sp)
      TYPE(cp_fm_type), POINTER                          :: matrix
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      CHARACTER(len=*), INTENT(in), OPTIONAL             :: name
      LOGICAL, INTENT(in), OPTIONAL                      :: use_sp

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_create', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ncol_local, npcol, nprow, &
                                                            nrow_local
      TYPE(cp_blacs_env_type), POINTER                   :: context

      CALL timeset(routineN, handle)

#if defined(__parallel) && ! defined(__SCALAPACK)
      CPABORT("full matrices need scalapack for parallel runs")
#endif

      CPASSERT(ASSOCIATED(matrix_struct))
      ALLOCATE (matrix)
      context => matrix_struct%context
      matrix%matrix_struct => matrix_struct
      CALL cp_fm_struct_retain(matrix%matrix_struct)
      last_fm_id_nr = last_fm_id_nr + 1
      matrix%id_nr = last_fm_id_nr
      matrix%ref_count = 1
      matrix%print_count = 0

      matrix%use_sp = .FALSE.
      IF (PRESENT(use_sp)) matrix%use_sp = use_sp

      nprow = context%num_pe(1)
      npcol = context%num_pe(2)
      NULLIFY (matrix%local_data)
      NULLIFY (matrix%local_data_sp)

      ! OK, we allocate here at least a 1 x 1 matrix
      ! this must (and is) compatible with the descinit call
      ! in cp_fm_struct
      nrow_local = matrix_struct%local_leading_dimension
      ncol_local = MAX(1, matrix_struct%ncol_locals(context%mepos(2)))
      IF (matrix%use_sp) THEN
         ALLOCATE (matrix%local_data_sp(nrow_local, ncol_local))
      ELSE
         ALLOCATE (matrix%local_data(nrow_local, ncol_local))
      ENDIF

      ! JVDV we should remove this, as it is up to the user to zero afterwards
      IF (matrix%use_sp) THEN
         CALL scopy(nrow_local*ncol_local, 0.0_sp, 0, matrix%local_data_sp, 1)
      ELSE
         CALL dcopy(nrow_local*ncol_local, 0.0_dp, 0, matrix%local_data, 1)
      ENDIF

      IF (PRESENT(name)) THEN
         matrix%name = name
      ELSE
         matrix%name = 'full matrix'//cp_to_string(matrix%id_nr)
      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_fm_create

! **************************************************************************************************
!> \brief retains a full matrix
!> \param matrix the matrix to retain
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE cp_fm_retain(matrix)
      TYPE(cp_fm_type), POINTER                          :: matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_retain', routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(matrix))
      CPASSERT(matrix%ref_count > 0)
      matrix%ref_count = matrix%ref_count + 1

   END SUBROUTINE cp_fm_retain

! **************************************************************************************************
!> \brief releases a full matrix
!> \param matrix the matrix to release
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE cp_fm_release(matrix)
      TYPE(cp_fm_type), POINTER                          :: matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_release', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(matrix)) THEN
         CPASSERT(matrix%ref_count > 0)
         matrix%ref_count = matrix%ref_count - 1
         IF (matrix%ref_count < 1) THEN
            IF (ASSOCIATED(matrix%local_data)) THEN
               DEALLOCATE (matrix%local_data)
            END IF
            IF (ASSOCIATED(matrix%local_data_sp)) THEN
               DEALLOCATE (matrix%local_data_sp)
            END IF
            matrix%name = ""
            CALL cp_fm_struct_release(matrix%matrix_struct)
            DEALLOCATE (matrix)
         END IF
      END IF
      NULLIFY (matrix)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_release

! **************************************************************************************************
!> \brief fills a matrix with random numbers
!> \param matrix : to be initialized
!> \param ncol : numbers of cols to fill
!> \param start_col : starting at coll number
!> \author Joost VandeVondele
!> \note
!>      the value of a_ij is independent of the number of cpus
! **************************************************************************************************
   SUBROUTINE cp_fm_init_random(matrix, ncol, start_col)
      TYPE(cp_fm_type), POINTER                          :: matrix
      INTEGER, INTENT(IN), OPTIONAL                      :: ncol, start_col

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_init_random', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, icol_global, icol_local, irow_local, my_ncol, my_start_col, ncol_global, &
         ncol_local, nrow_global, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL, SAVE                                      :: FIRST = .TRUE.
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: buff
      REAL(KIND=dp), DIMENSION(3, 2), SAVE               :: seed
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: local_data
      TYPE(rng_stream_type)                              :: rng

      CALL timeset(routineN, handle)

      ! first setup a parallel rng
      IF (FIRST) THEN
         seed(:, :) = RESHAPE((/1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp/), (/3, 2/))
         FIRST = .FALSE.
      ENDIF
      ! guarantee same seed over all tasks
      CALL mp_bcast(seed, 0, matrix%matrix_struct%para_env%group)

      rng = rng_stream_type("cp_fm_init_random_stream", distribution_type=UNIFORM, &
                            extended_precision=.TRUE., seed=seed)

      CPASSERT(.NOT. matrix%use_sp)

      CALL cp_fm_get_info(matrix, nrow_global=nrow_global, ncol_global=ncol_global, &
                          nrow_local=nrow_local, ncol_local=ncol_local, &
                          local_data=local_data, &
                          row_indices=row_indices, col_indices=col_indices)

      my_start_col = 1
      IF (PRESENT(start_col)) my_start_col = start_col
      my_ncol = matrix%matrix_struct%ncol_global
      IF (PRESENT(ncol)) my_ncol = ncol

      IF (ncol_global < (my_start_col + my_ncol - 1)) &
         CPABORT("ncol_global>=(my_start_col+my_ncol-1)")

      ALLOCATE (buff(nrow_global))

      ! each global row has its own substream, in order to reach the stream for the local col,
      ! we just reset to the next substream
      ! following this, we fill the full buff with random numbers, and pick those we need
      icol_global = 0
      DO icol_local = 1, ncol_local
         CPASSERT(col_indices(icol_local) > icol_global)
         DO
            CALL rng%reset_to_next_substream()
            icol_global = icol_global + 1
            IF (icol_global == col_indices(icol_local)) EXIT
         ENDDO
         CALL rng%fill(buff)
         DO irow_local = 1, nrow_local
            local_data(irow_local, icol_local) = buff(row_indices(irow_local))
         ENDDO
      ENDDO

      DEALLOCATE (buff)

      ! store seed before deletion (unclear if this is the proper seed)
      CALL rng%get(ig=seed)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_init_random

! **************************************************************************************************
!> \brief set all elements of a matrix to the same value,
!>      and optionally the diagonal to a different one
!> \param matrix input matrix
!> \param alpha scalar used to set all elements of the matrix
!> \param beta scalar used to set diagonal of the matrix
!> \note
!>      can be used to zero a matrix
!>      can be used to create a unit matrix (I-matrix) alpha=0.0_dp beta=1.0_dp
! **************************************************************************************************
   SUBROUTINE cp_fm_set_all(matrix, alpha, beta)

      TYPE(cp_fm_type), POINTER                          :: matrix
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: beta

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_set_all', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, matrix_size, n

      CALL timeset(routineN, handle)

      IF (matrix%use_sp) THEN
         matrix_size = SIZE(matrix%local_data_sp, 1)*SIZE(matrix%local_data_sp, 2)
         CALL scopy(matrix_size, REAL(alpha, sp), 0, matrix%local_data_sp, 1)
      ELSE
         matrix_size = SIZE(matrix%local_data, 1)*SIZE(matrix%local_data, 2)
         CALL dcopy(matrix_size, alpha, 0, matrix%local_data, 1)
      ENDIF

      IF (PRESENT(beta)) THEN
         CPASSERT(.NOT. matrix%use_sp)
         n = MIN(matrix%matrix_struct%nrow_global, matrix%matrix_struct%ncol_global)
         DO i = 1, n
            CALL cp_fm_set_element(matrix, i, i, beta)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_fm_set_all

! **************************************************************************************************
!> \brief returns the diagonal elements of a fm
!> \param matrix ...
!> \param diag ...
! **************************************************************************************************
   SUBROUTINE cp_fm_get_diag(matrix, diag)

      IMPLICIT NONE

      ! arguments
      TYPE(cp_fm_type), POINTER                  :: matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: diag

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_get_diag', &
                                     routineP = moduleN//':'//routineN

      ! locals
      INTEGER :: i, nrow_global

#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9) :: desca
      TYPE(cp_blacs_env_type), POINTER :: context
      INTEGER :: icol_local, ipcol, iprow, irow_local, mypcol, myprow, npcol, &
                 nprow
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: a
      REAL(KIND=sp), DIMENSION(:, :), POINTER :: a_sp
#endif

      CALL cp_fm_get_info(matrix, nrow_global=nrow_global)

#if defined(__SCALAPACK)
      diag = 0.0_dp
      context => matrix%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      nprow = context%num_pe(1)
      npcol = context%num_pe(2)

      a => matrix%local_data
      a_sp => matrix%local_data_sp
      desca(:) = matrix%matrix_struct%descriptor(:)

      DO i = 1, nrow_global
         CALL infog2l(i, i, desca, nprow, npcol, myprow, mypcol, &
                      irow_local, icol_local, iprow, ipcol)
         IF ((iprow == myprow) .AND. (ipcol == mypcol)) THEN
            IF (matrix%use_sp) THEN
               diag(i) = REAL(a_sp(irow_local, icol_local), dp)
            ELSE
               diag(i) = a(irow_local, icol_local)
            ENDIF
         ENDIF
      END DO
#else
      DO i = 1, nrow_global
         IF (matrix%use_sp) THEN
            diag(i) = REAL(matrix%local_data_sp(i, i), dp)
         ELSE
            diag(i) = matrix%local_data(i, i)
         ENDIF
      ENDDO
#endif
      CALL mp_sum(diag, matrix%matrix_struct%para_env%group)

   END SUBROUTINE cp_fm_get_diag

! **************************************************************************************************
!> \brief returns an element of a fm
!>      this value is valid on every cpu
!>      using this call is expensive
!> \param matrix the matrix to read
!> \param irow_global the row
!> \param icol_global the col
!> \param alpha the value of matrix(irow_global, icol_global)
!> \param local true if the element is on this cpu, false otherwise
!> \note
!>      - modified semantics. now this function always returns the value
!>        previously the value was zero on cpus that didn't own the relevant
!>        part of the matrix (Joost VandeVondele, May 2003)
!>      - usage of the function should be avoided, as it is likely to rather slow
!>        using row_indices/col_indices/local_data + some smart scheme normally
!>        yields a real parallel code
! **************************************************************************************************
   SUBROUTINE cp_fm_get_element(matrix, irow_global, icol_global, alpha, local)

      IMPLICIT NONE

      ! arguments
      TYPE(cp_fm_type), POINTER          :: matrix
      REAL(KIND=dp), INTENT(OUT)                     :: alpha
      INTEGER, INTENT(IN)                       :: icol_global, &
                                                   irow_global
      LOGICAL, INTENT(OUT), OPTIONAL            :: local

      ! locals
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9) :: desca
      TYPE(cp_blacs_env_type), POINTER :: context
      INTEGER :: icol_local, ipcol, iprow, irow_local, mypcol, myprow, npcol, &
                 nprow
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: a
#endif

#if defined(__SCALAPACK)
      context => matrix%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      nprow = context%num_pe(1)
      npcol = context%num_pe(2)

      a => matrix%local_data
      desca(:) = matrix%matrix_struct%descriptor(:)

      CALL infog2l(irow_global, icol_global, desca, nprow, npcol, myprow, mypcol, &
                   irow_local, icol_local, iprow, ipcol)

      IF ((iprow == myprow) .AND. (ipcol == mypcol)) THEN
         alpha = a(irow_local, icol_local)
         CALL cp_blacs_dgebs2d(context%group, 'All', ' ', 1, 1, alpha, 1)
         IF (PRESENT(local)) local = .TRUE.
      ELSE
         CALL cp_blacs_dgebr2d(context%group, 'All', ' ', 1, 1, alpha, 1, iprow, ipcol)
         IF (PRESENT(local)) local = .FALSE.
      END IF

#else
      IF (PRESENT(local)) local = .TRUE.
      alpha = matrix%local_data(irow_global, icol_global)
#endif

   END SUBROUTINE cp_fm_get_element

! **************************************************************************************************
!> \brief sets an element of a matrix
!> \param matrix ...
!> \param irow_global ...
!> \param icol_global ...
!> \param alpha ...
!> \note
!>      we expect all cpus to have the same arguments in the call to this function
!>      (otherwise one should use local_data tricks)
! **************************************************************************************************
   SUBROUTINE cp_fm_set_element(matrix, irow_global, icol_global, alpha)
      TYPE(cp_fm_type), POINTER                :: matrix
      INTEGER, INTENT(IN)                      :: irow_global, icol_global
      REAL(KIND=dp), INTENT(IN)                :: alpha

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_set_element', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: mypcol, myprow, npcol, nprow
      TYPE(cp_blacs_env_type), POINTER         :: context
#if defined(__SCALAPACK)
      INTEGER                                  :: icol_local, ipcol, iprow, &
                                                  irow_local
      INTEGER, DIMENSION(9)                    :: desca
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
#endif

      context => matrix%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      nprow = context%num_pe(1)
      npcol = context%num_pe(2)

      CPASSERT(.NOT. matrix%use_sp)

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
   END SUBROUTINE cp_fm_set_element

! **************************************************************************************************
!> \brief sets a submatrix of a full matrix
!>       fm(start_row:start_row+n_rows,start_col:start_col+n_cols)
!>       = alpha*op(new_values)(1:n_rows,1:n_cols)+ beta
!>       * fm(start_row:start_row+n_rows,start_col:start_col+n_cols)
!> \param fm the full to change
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
!>      07.2002 created borrowing from Joost's blacs_replicated_copy [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      optimized for full column updates and alpha=1.0, beta=0.0
!>      the new_values need to be valid on all cpus
! **************************************************************************************************
   SUBROUTINE cp_fm_set_submatrix(fm, new_values, start_row, &
                                  start_col, n_rows, n_cols, alpha, beta, transpose)
      TYPE(cp_fm_type), POINTER                          :: fm
      REAL(KIND=dp), DIMENSION(:, :), INTENT(in)         :: new_values
      INTEGER, INTENT(in), OPTIONAL                      :: start_row, start_col, n_rows, n_cols
      REAL(KIND=dp), INTENT(in), OPTIONAL                :: alpha, beta
      LOGICAL, INTENT(in), OPTIONAL                      :: transpose

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_set_submatrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, i0, j, j0, ncol, ncol_block, &
                                                            ncol_global, ncol_local, nrow, &
                                                            nrow_block, nrow_global, nrow_local, &
                                                            this_col, this_row
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: tr_a
      REAL(KIND=dp)                                      :: al, be
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: full_block

      al = 1.0_dp; be = 0.0_dp; i0 = 1; j0 = 1; tr_a = .FALSE.
      ! can be called too many times, making it a bit useless
      ! CALL timeset(routineN//','//moduleN,handle)

      CPASSERT(ASSOCIATED(fm))
      CPASSERT(fm%ref_count > 0)

      CPASSERT(.NOT. fm%use_sp)

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

      full_block => fm%local_data

      CALL cp_fm_get_info(matrix=fm, &
                          nrow_global=nrow_global, ncol_global=ncol_global, &
                          nrow_block=nrow_block, ncol_block=ncol_block, &
                          nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)

      IF (al == 1.0 .AND. be == 0.0) THEN
         DO j = 1, ncol_local
            this_col = col_indices(j) - j0 + 1
            IF (this_col .GE. 1 .AND. this_col .LE. ncol) THEN
               IF (tr_a) THEN
                  IF (i0 == 1 .AND. nrow_global == nrow) THEN
                     DO i = 1, nrow_local
                        full_block(i, j) = new_values(this_col, row_indices(i))
                     END DO
                  ELSE
                     DO i = 1, nrow_local
                        this_row = row_indices(i) - i0 + 1
                        IF (this_row >= 1 .AND. this_row <= nrow) THEN
                           full_block(i, j) = new_values(this_col, this_row)
                        END IF
                     END DO
                  END IF
               ELSE
                  IF (i0 == 1 .AND. nrow_global == nrow) THEN
                     DO i = 1, nrow_local
                        full_block(i, j) = new_values(row_indices(i), this_col)
                     END DO
                  ELSE
                     DO i = 1, nrow_local
                        this_row = row_indices(i) - i0 + 1
                        IF (this_row >= 1 .AND. this_row <= nrow) THEN
                           full_block(i, j) = new_values(this_row, this_col)
                        END IF
                     END DO
                  END IF
               END IF
            END IF
         END DO
      ELSE
         DO j = 1, ncol_local
            this_col = col_indices(j) - j0 + 1
            IF (this_col .GE. 1 .AND. this_col .LE. ncol) THEN
               IF (tr_a) THEN
                  DO i = 1, nrow_local
                     this_row = row_indices(i) - i0 + 1
                     IF (this_row >= 1 .AND. this_row <= nrow) THEN
                        full_block(i, j) = al*new_values(this_col, this_row) + &
                                           be*full_block(i, j)
                     END IF
                  END DO
               ELSE
                  DO i = 1, nrow_local
                     this_row = row_indices(i) - i0 + 1
                     IF (this_row >= 1 .AND. this_row <= nrow) THEN
                        full_block(i, j) = al*new_values(this_row, this_col) + &
                                           be*full_block(i, j)
                     END IF
                  END DO
               END IF
            END IF
         END DO
      END IF

      ! CALL timestop(handle)

   END SUBROUTINE cp_fm_set_submatrix

! **************************************************************************************************
!> \brief gets a submatrix of a full matrix
!>       op(target_m)(1:n_rows,1:n_cols)
!>       =fm(start_row:start_row+n_rows,start_col:start_col+n_cols)
!>      target_m is replicated on all cpus
!>      using this call is expensive
!> \param fm the full you want to get the info from
!> \param target_m a replicated full matrix that will contain the result
!> \param start_row the starting row of b_matrix (defaults to 1)
!> \param start_col the starting col of b_matrix (defaults to 1)
!> \param n_rows the number of row to change in b (defaults to
!>        size(op(new_values),1))
!> \param n_cols the number of columns to change in b (defaults to
!>        size(op(new_values),2))
!> \param transpose if target_m should be transposed: if true
!>        op(target_m)=target_m^T, else op(target_m)=target_m
!>        (defaults to false)
!> \par History
!>      07.2002 created borrowing from Joost's blacs_replicated_copy [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      optimized for full column updates. Zeros out a little too much
!>      of target_m
!>      the target_m is replicated and valid on all cpus
! **************************************************************************************************
   SUBROUTINE cp_fm_get_submatrix(fm, target_m, start_row, &
                                  start_col, n_rows, n_cols, transpose)
      TYPE(cp_fm_type), POINTER                          :: fm
      REAL(KIND=dp), DIMENSION(:, :), INTENT(out)        :: target_m
      INTEGER, INTENT(in), OPTIONAL                      :: start_row, start_col, n_rows, n_cols
      LOGICAL, INTENT(in), OPTIONAL                      :: transpose

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_get_submatrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, i0, j, j0, ncol, ncol_global, &
                                                            ncol_local, nrow, nrow_global, &
                                                            nrow_local, this_col, this_row
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: tr_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: full_block
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      i0 = 1; j0 = 1; tr_a = .FALSE.

      CPASSERT(ASSOCIATED(fm))
      CPASSERT(fm%ref_count > 0)

      CPASSERT(.NOT. fm%use_sp)

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
      IF (SIZE(target_m, 1)*SIZE(target_m, 2) /= 0) THEN
         CALL dcopy(SIZE(target_m, 1)*SIZE(target_m, 2), 0.0_dp, 0, target_m(1, 1), 1)
      END IF
#endif

      CALL cp_fm_get_info(matrix=fm, &
                          nrow_global=nrow_global, ncol_global=ncol_global, &
                          nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)

      DO j = 1, ncol_local
         this_col = col_indices(j) - j0 + 1
         IF (this_col .GE. 1 .AND. this_col .LE. ncol) THEN
            IF (tr_a) THEN
               IF (i0 == 1 .AND. nrow_global == nrow) THEN
                  DO i = 1, nrow_local
                     target_m(this_col, row_indices(i)) = full_block(i, j)
                  END DO
               ELSE
                  DO i = 1, nrow_local
                     this_row = row_indices(i) - i0 + 1
                     IF (this_row >= 1 .AND. this_row <= nrow) THEN
                        target_m(this_col, this_row) = full_block(i, j)
                     END IF
                  END DO
               END IF
            ELSE
               IF (i0 == 1 .AND. nrow_global == nrow) THEN
                  DO i = 1, nrow_local
                     target_m(row_indices(i), this_col) = full_block(i, j)
                  END DO
               ELSE
                  DO i = 1, nrow_local
                     this_row = row_indices(i) - i0 + 1
                     IF (this_row >= 1 .AND. this_row <= nrow) THEN
                        target_m(this_row, this_col) = full_block(i, j)
                     END IF
                  END DO
               END IF
            END IF
         END IF
      END DO

      CALL mp_sum(target_m, para_env%group)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_get_submatrix

! **************************************************************************************************
!> \brief returns all kind of information about the full matrix
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
!> \param nrow_locals ...
!> \param ncol_locals ...
!> \param matrix_struct ...
!> \param para_env ...
!> \note
!>       see also cp_fm_struct for explanation
!>       - nrow_local, ncol_local, row_indices, col_indices, local_data are hooks for efficient
!>         access to the local blacs block
! **************************************************************************************************
   SUBROUTINE cp_fm_get_info(matrix, name, nrow_global, ncol_global, &
                             nrow_block, ncol_block, nrow_local, ncol_local, &
                             row_indices, col_indices, local_data, context, &
                             nrow_locals, ncol_locals, matrix_struct, para_env)

      TYPE(cp_fm_type), POINTER                          :: matrix
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL            :: name
      INTEGER, INTENT(OUT), OPTIONAL                     :: nrow_global, ncol_global, nrow_block, &
                                                            ncol_block, nrow_local, ncol_local
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: row_indices, col_indices
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: local_data
      TYPE(cp_blacs_env_type), OPTIONAL, POINTER         :: context
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: nrow_locals, ncol_locals
      TYPE(cp_fm_struct_type), OPTIONAL, POINTER         :: matrix_struct
      TYPE(cp_para_env_type), OPTIONAL, POINTER          :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_get_info', routineP = moduleN//':'//routineN

      IF (PRESENT(name)) name = matrix%name
      IF (PRESENT(context)) context => matrix%matrix_struct%context
      IF (PRESENT(matrix_struct)) matrix_struct => matrix%matrix_struct
      IF (PRESENT(local_data)) local_data => matrix%local_data ! not hiding things anymore :-(
      IF (PRESENT(row_indices)) THEN
         CALL cp_fm_struct_get(matrix%matrix_struct, row_indices=row_indices)
      ENDIF
      IF (PRESENT(col_indices)) THEN
         CALL cp_fm_struct_get(matrix%matrix_struct, col_indices=col_indices)
      ENDIF
      IF (PRESENT(nrow_locals)) THEN
         nrow_locals => matrix%matrix_struct%nrow_locals
      END IF
      IF (PRESENT(ncol_locals)) THEN
         ncol_locals => matrix%matrix_struct%ncol_locals
      END IF

      CALL cp_fm_struct_get(matrix%matrix_struct, nrow_local=nrow_local, &
                            ncol_local=ncol_local, nrow_global=nrow_global, &
                            ncol_global=ncol_global, nrow_block=nrow_block, &
                            ncol_block=ncol_block)

      IF (PRESENT(para_env)) para_env => matrix%matrix_struct%para_env

   END SUBROUTINE cp_fm_get_info

! **************************************************************************************************
!> \brief Write nicely formatted info about the FM to the given I/O unit (including the underlying FM struct)
!> \param matrix a cp_fm_type instance
!> \param io_unit the I/O unit to use for writing
! **************************************************************************************************
   SUBROUTINE cp_fm_write_info(matrix, io_unit)
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: matrix
      INTEGER, INTENT(IN)                                :: io_unit

      WRITE (io_unit, '(/,A,A12)') "CP_FM | Name:                           ", matrix%name
      CALL cp_fm_struct_write_info(matrix%matrix_struct, io_unit)
   END SUBROUTINE cp_fm_write_info

! **************************************************************************************************
!> \brief find the maximum absolute value of the matrix element
!>      maxval(abs(matrix))
!> \param matrix ...
!> \param a_max ...
!> \param ir_max ...
!> \param ic_max ...
! **************************************************************************************************
   SUBROUTINE cp_fm_maxabsval(matrix, a_max, ir_max, ic_max)
      TYPE(cp_fm_type), POINTER                          :: matrix
      REAL(KIND=dp), INTENT(OUT)                         :: a_max
      INTEGER, INTENT(OUT), OPTIONAL                     :: ir_max, ic_max

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_maxabsval', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ic_max_local, ir_max_local, &
                                                            j, mepos, ncol_local, nrow_local, &
                                                            num_pe
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ic_max_vec, ir_max_vec
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(dp)                                           :: my_max
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: a_max_vec
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: my_block
      REAL(KIND=sp), DIMENSION(:, :), POINTER            :: my_block_sp

      CALL timeset(routineN, handle)

      my_block => matrix%local_data
      my_block_sp => matrix%local_data_sp

      CALL cp_fm_get_info(matrix, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)

      IF (matrix%use_sp) THEN
         a_max = REAL(MAXVAL(ABS(my_block_sp(1:nrow_local, 1:ncol_local))), dp)
      ELSE
         a_max = MAXVAL(ABS(my_block(1:nrow_local, 1:ncol_local)))
      ENDIF

      IF (PRESENT(ir_max)) THEN
         num_pe = matrix%matrix_struct%para_env%num_pe
         mepos = matrix%matrix_struct%para_env%mepos
         ALLOCATE (ir_max_vec(0:num_pe - 1))
         ir_max_vec(0:num_pe - 1) = 0
         ALLOCATE (ic_max_vec(0:num_pe - 1))
         ic_max_vec(0:num_pe - 1) = 0
         ALLOCATE (a_max_vec(0:num_pe - 1))
         a_max_vec(0:num_pe - 1) = 0.0_dp
         my_max = 0.0_dp

         IF ((ncol_local > 0) .AND. (nrow_local > 0)) THEN
            DO i = 1, ncol_local
               DO j = 1, nrow_local
                  IF (matrix%use_sp) THEN
                     IF (ABS(my_block_sp(j, i)) > my_max) THEN
                        my_max = REAL(my_block_sp(j, i), dp)
                        ir_max_local = j
                        ic_max_local = i
                     ENDIF
                  ELSE
                     IF (ABS(my_block(j, i)) > my_max) THEN
                        my_max = my_block(j, i)
                        ir_max_local = j
                        ic_max_local = i
                     ENDIF
                  END IF
               END DO
            END DO

            a_max_vec(mepos) = my_max
            ir_max_vec(mepos) = row_indices(ir_max_local)
            ic_max_vec(mepos) = col_indices(ic_max_local)

         END IF

         CALL mp_sum(a_max_vec, matrix%matrix_struct%para_env%group)
         CALL mp_sum(ir_max_vec, matrix%matrix_struct%para_env%group)
         CALL mp_sum(ic_max_vec, matrix%matrix_struct%para_env%group)

         my_max = 0.0_dp
         DO i = 0, num_pe - 1
            IF (a_max_vec(i) > my_max) THEN
               ir_max = ir_max_vec(i)
               ic_max = ic_max_vec(i)
            END IF
         END DO

         DEALLOCATE (ir_max_vec, ic_max_vec, a_max_vec)
         CPASSERT(ic_max > 0)
         CPASSERT(ir_max > 0)

      END IF

      CALL mp_max(a_max, matrix%matrix_struct%para_env%group)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_maxabsval

! **************************************************************************************************
!> \brief find the maximum over the rows of the sum of the absolute values of the elements of a given row
!>      = || A ||_infinity
!> \param matrix ...
!> \param a_max ...
!> \note
!>      for a real symmetric matrix it holds that || A ||_2 = |lambda_max| < || A ||_infinity
!>      Hence this can be used to estimate an upper bound for the eigenvalues of a matrix
!>      http://mathworld.wolfram.com/MatrixNorm.html
!>      (but the bound is not so tight in the general case)
! **************************************************************************************************
   SUBROUTINE cp_fm_maxabsrownorm(matrix, a_max)
      TYPE(cp_fm_type), POINTER                          :: matrix
      REAL(KIND=dp), INTENT(OUT)                         :: a_max

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_maxabsrownorm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, ncol_local, nrow_global, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: values
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: my_block

      CALL timeset(routineN, handle)

      my_block => matrix%local_data

      CPASSERT(.NOT. matrix%use_sp)

      CALL cp_fm_get_info(matrix, row_indices=row_indices, nrow_global=nrow_global, &
                          nrow_local=nrow_local, ncol_local=ncol_local)

      ! the efficiency could be improved by making use of the row-col distribution of scalapack
      ALLOCATE (values(nrow_global))
      values = 0.0_dp
      DO j = 1, ncol_local
         DO i = 1, nrow_local
            values(row_indices(i)) = values(row_indices(i)) + ABS(my_block(i, j))
         ENDDO
      ENDDO
      CALL mp_sum(values, matrix%matrix_struct%para_env%group)
      a_max = MAXVAL(values)
      DEALLOCATE (values)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief find the inorm of each column norm_{j}= sqrt( \sum_{i} A_{ij}*A_{ij} )
!> \param matrix ...
!> \param norm_array ...
! **************************************************************************************************
   SUBROUTINE cp_fm_vectorsnorm(matrix, norm_array)
      TYPE(cp_fm_type), POINTER                          :: matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: norm_array

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_vectorsnorm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, ncol_global, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: my_block

      CALL timeset(routineN, handle)

      my_block => matrix%local_data

      CPASSERT(.NOT. matrix%use_sp)

      CALL cp_fm_get_info(matrix, col_indices=col_indices, ncol_global=ncol_global, &
                          nrow_local=nrow_local, ncol_local=ncol_local)

      ! the efficiency could be improved by making use of the row-col distribution of scalapack
      norm_array = 0.0_dp
      DO j = 1, ncol_local
         DO i = 1, nrow_local
            norm_array(col_indices(j)) = norm_array(col_indices(j)) + my_block(i, j)*my_block(i, j)
         ENDDO
      ENDDO
      CALL mp_sum(norm_array, matrix%matrix_struct%para_env%group)
      norm_array = SQRT(norm_array)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief summing up all the elements along the matrix's i-th index
!>        \f$ \mathrm{sum}_{j} = \sum_{i} A_{ij} \f$
!>        or
!>        \f$ \mathrm{sum}_{i} = \sum_{j} A_{ij} \f$
!> \param matrix     an input matrix A
!> \param sum_array  sums of elements in each column/row
!> \param dir ...
!> \note  forked from cp_fm_vectorsnorm() to be used with
!>        the maximum overlap method
!>        added row variation
! **************************************************************************************************
   SUBROUTINE cp_fm_vectorssum(matrix, sum_array, dir)
      TYPE(cp_fm_type), POINTER                          :: matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: sum_array
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL             :: dir

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_vectorssum', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: docol
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: my_block

      CALL timeset(routineN, handle)

      IF (PRESENT(dir)) THEN
         IF (dir == 'c' .OR. dir == 'C') THEN
            docol = .TRUE.
         ELSEIF (dir == 'r' .OR. dir == 'R') THEN
            docol = .FALSE.
         ELSE
            CPABORT('Wrong argument DIR')
         END IF
      ELSE
         docol = .TRUE.
      END IF

      my_block => matrix%local_data

      CPASSERT(.NOT. matrix%use_sp)

      CALL cp_fm_get_info(matrix, col_indices=col_indices, row_indices=row_indices, &
                          nrow_local=nrow_local, ncol_local=ncol_local)

      ! the efficiency could be improved by making use of the row-col distribution of scalapack
      sum_array(:) = 0.0_dp
      IF (docol) THEN
      DO j = 1, ncol_local
         DO i = 1, nrow_local
            sum_array(col_indices(j)) = sum_array(col_indices(j)) + my_block(i, j)
         ENDDO
      ENDDO
      ELSE
      DO j = 1, ncol_local
         DO i = 1, nrow_local
            sum_array(row_indices(i)) = sum_array(row_indices(i)) + my_block(i, j)
         ENDDO
      ENDDO
      END IF
      CALL mp_sum(sum_array, matrix%matrix_struct%para_env%group)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief copy one identically sized matrix in the other
!> \param source ...
!> \param destination ...
!> \note
!>      see also cp_fm_to_fm_columns
! **************************************************************************************************
   SUBROUTINE cp_fm_to_fm_matrix(source, destination)

      TYPE(cp_fm_type), POINTER                          :: source, destination

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_to_fm_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, npcol, nprow

      CALL timeset(routineN, handle)

      nprow = source%matrix_struct%context%num_pe(1)
      npcol = source%matrix_struct%context%num_pe(2)

      IF ((.NOT. cp2k_is_parallel) .OR. &
          cp_fm_struct_equivalent(source%matrix_struct, &
                                  destination%matrix_struct)) THEN
         IF (source%use_sp .AND. destination%use_sp) THEN
            IF (SIZE(source%local_data_sp, 1) /= SIZE(destination%local_data_sp, 1) .OR. &
                SIZE(source%local_data_sp, 2) /= SIZE(destination%local_data_sp, 2)) &
               CALL cp_abort(__LOCATION__, &
                             "Cannot copy full matrix <"//TRIM(source%name)// &
                             "> to full matrix <"//TRIM(destination%name)// &
                             ">. The local_data blocks have different sizes.")
            CALL scopy(SIZE(source%local_data_sp, 1)*SIZE(source%local_data_sp, 2), &
                       source%local_data_sp(1, 1), 1, destination%local_data_sp(1, 1), 1)
         ELSEIF (source%use_sp .AND. .NOT. destination%use_sp) THEN
            IF (SIZE(source%local_data_sp, 1) /= SIZE(destination%local_data, 1) .OR. &
                SIZE(source%local_data_sp, 2) /= SIZE(destination%local_data, 2)) &
               CALL cp_abort(__LOCATION__, &
                             "Cannot copy full matrix <"//TRIM(source%name)// &
                             "> to full matrix <"//TRIM(destination%name)// &
                             ">. The local_data blocks have different sizes.")
            destination%local_data = REAL(source%local_data_sp, dp)
         ELSEIF (.NOT. source%use_sp .AND. destination%use_sp) THEN
            IF (SIZE(source%local_data, 1) /= SIZE(destination%local_data_sp, 1) .OR. &
                SIZE(source%local_data, 2) /= SIZE(destination%local_data_sp, 2)) &
               CALL cp_abort(__LOCATION__, &
                             "Cannot copy full matrix <"//TRIM(source%name)// &
                             "> to full matrix <"//TRIM(destination%name)// &
                             ">. The local_data blocks have different sizes.")
            destination%local_data_sp = REAL(source%local_data, sp)
         ELSE
            IF (SIZE(source%local_data, 1) /= SIZE(destination%local_data, 1) .OR. &
                SIZE(source%local_data, 2) /= SIZE(destination%local_data, 2)) &
               CALL cp_abort(__LOCATION__, &
                             "Cannot copy full matrix <"//TRIM(source%name)// &
                             "> to full matrix <"//TRIM(destination%name)// &
                             ">. The local_data blocks have different sizes.")
            CALL dcopy(SIZE(source%local_data, 1)*SIZE(source%local_data, 2), &
                       source%local_data(1, 1), 1, destination%local_data(1, 1), 1)
         ENDIF
      ELSE
         CPABORT("")
      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_fm_to_fm_matrix

! **************************************************************************************************
!> \brief copy just a subset of columns of a fm to a fm
!> \param msource ...
!> \param mtarget ...
!> \param ncol ...
!> \param source_start ...
!> \param target_start ...
! **************************************************************************************************
   SUBROUTINE cp_fm_to_fm_columns(msource, mtarget, ncol, source_start, &
                                  target_start)

      TYPE(cp_fm_type), POINTER                :: msource, mtarget
      INTEGER, INTENT(IN)                      :: ncol
      INTEGER, INTENT(IN), OPTIONAL            :: source_start, target_start

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_to_fm_columns', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: handle, n, ss, ts
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b
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
      DO i = 0, ncol - 1
         CALL pdcopy(n, a(1, 1), 1, ss + i, desca, 1, b(1, 1), 1, ts + i, descb, 1)
      END DO
#else
      CALL dcopy(ncol*n, a(1, ss), 1, b(1, ts), 1)
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_to_fm_columns

! **************************************************************************************************
!> \brief copy just a triangular matrix
!> \param msource ...
!> \param mtarget ...
!> \param uplo ...
! **************************************************************************************************
   SUBROUTINE cp_fm_to_fm_triangular(msource, mtarget, uplo)

      TYPE(cp_fm_type), POINTER                :: msource, mtarget
      CHARACTER(LEN=*), INTENT(IN)             :: uplo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_to_fm_triangular', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: handle, ncol, nrow
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                    :: desca, descb
#endif

      CALL timeset(routineN, handle)

      nrow = msource%matrix_struct%nrow_global
      ncol = msource%matrix_struct%ncol_global

      a => msource%local_data
      b => mtarget%local_data

#if defined(__SCALAPACK)
      desca(:) = msource%matrix_struct%descriptor(:)
      descb(:) = mtarget%matrix_struct%descriptor(:)
      CALL pdlacpy(uplo, nrow, ncol, a(1, 1), 1, 1, desca, b(1, 1), 1, 1, descb)
#else
      CALL dlacpy(uplo, nrow, ncol, a(1, 1), nrow, b(1, 1), nrow)
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_to_fm_triangular

! **************************************************************************************************
!> \brief copy just a part ot the matrix
!> \param msource ...
!> \param mtarget ...
!> \param nrow ...
!> \param ncol ...
!> \param s_firstrow ...
!> \param s_firstcol ...
!> \param t_firstrow ...
!> \param t_firstcol ...
! **************************************************************************************************

   SUBROUTINE cp_fm_to_fm_submat(msource, mtarget, nrow, ncol, s_firstrow, s_firstcol, t_firstrow, t_firstcol)

      TYPE(cp_fm_type), POINTER                :: msource, mtarget
      INTEGER, INTENT(IN)                      :: nrow, ncol, s_firstrow, &
                                                  s_firstcol, t_firstrow, &
                                                  t_firstcol

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_to_fm_submat', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: handle, i, na, nb, ss, ts
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                    :: desca, descb
#endif

      CALL timeset(routineN, handle)

      a => msource%local_data
      b => mtarget%local_data

      na = msource%matrix_struct%nrow_global
      nb = mtarget%matrix_struct%nrow_global
!    nrow must be <= na and nb
      IF (nrow > na) &
         CPABORT("cannot copy because nrow > number of rows of source matrix")
      IF (nrow > nb) &
         CPABORT("cannot copy because nrow > number of rows of target matrix")
      na = msource%matrix_struct%ncol_global
      nb = mtarget%matrix_struct%ncol_global
!    ncol must be <= na_col and nb_col
      IF (ncol > na) &
         CPABORT("cannot copy because nrow > number of rows of source matrix")
      IF (ncol > nb) &
         CPABORT("cannot copy because nrow > number of rows of target matrix")

#if defined(__SCALAPACK)
      desca(:) = msource%matrix_struct%descriptor(:)
      descb(:) = mtarget%matrix_struct%descriptor(:)
      DO i = 0, ncol - 1
         ss = s_firstcol + i
         ts = t_firstcol + i
         CALL pdcopy(nrow, a(1, 1), s_firstrow, ss, desca, 1, b(1, 1), t_firstrow, ts, descb, 1)
      END DO
#else
      DO i = 0, ncol - 1
         ss = s_firstcol + i
         ts = t_firstcol + i
         CALL dcopy(nrow, a(s_firstrow, ss), 1, b(t_firstrow, ts), 1)
      END DO
#endif

      CALL timestop(handle)
   END SUBROUTINE cp_fm_to_fm_submat

! **************************************************************************************************
!> \brief General copy of a fm matrix to another fm matrix.
!>        Uses non-blocking MPI rather than ScaLAPACK.
!>
!> \param source          input fm matrix
!> \param destination     output fm matrix
!> \param para_env        parallel environment corresponding to the BLACS env that covers all parts
!>                        of the input and output matrices
!> \par History
!>      31-Jan-2017 : Re-implemented using non-blocking MPI [IainB, MarkT]
! **************************************************************************************************
   SUBROUTINE cp_fm_copy_general(source, destination, para_env)
      TYPE(cp_fm_type), POINTER                          :: source, destination
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_copy_general', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(copy_info_type)                               :: info

      CALL timeset(routineN, handle)

      CALL cp_fm_start_copy_general(source, destination, para_env, info)
      IF (ASSOCIATED(destination)) THEN
         CALL cp_fm_finish_copy_general(destination, info)
      END IF
      IF (ASSOCIATED(source)) THEN
         CALL cp_fm_cleanup_copy_general(source, info)
      END IF

      CALL timestop(handle)
   END SUBROUTINE cp_fm_copy_general

! **************************************************************************************************
!> \brief Initiates the copy operation: get distribution data, post MPI isend and irecvs
!> \param source input fm matrix
!> \param destination output fm matrix
!> \param para_env parallel environment corresponding to the BLACS env that covers all parts
!>                        of the input and output matrices
!> \param info all of the data that will be needed to complete the copy operation
! **************************************************************************************************
   SUBROUTINE cp_fm_start_copy_general(source, destination, para_env, info)
      TYPE(cp_fm_type), POINTER                          :: source, destination
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(copy_info_type), INTENT(OUT)                  :: info

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_start_copy_general', &
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
      INTEGER, DIMENSION(:), POINTER :: recv_col_indices, recv_disp, recv_request, &
         recv_row_indices, send_col_indices, send_request, send_row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: src_blacs2mpi
      REAL(KIND=dp), DIMENSION(:), POINTER               :: recv_buf, send_buf
      TYPE(cp_fm_struct_type), POINTER                   :: recv_dist, send_dist

      CALL timeset(routineN, handle)

      IF (.NOT. cp2k_is_parallel) THEN
         ! Just copy all of the matrix data into a 'send buffer', to be unpacked later
         nrow_local_send = SIZE(source%local_data, 1)
         ncol_local_send = SIZE(source%local_data, 2)
         ALLOCATE (send_buf(nrow_local_send*ncol_local_send))
         info%send_buf => send_buf
         k = 0
         DO j = 1, ncol_local_send
            DO i = 1, nrow_local_send
               k = k + 1
               send_buf(k) = source%local_data(i, j)
            END DO
         END DO
      ELSE
         ! assumption of double precision data is carried forward from ScaLAPACK version
         IF (ASSOCIATED(source)) THEN
            IF (source%use_sp) CPABORT("only DP kind implemented")
         END IF
         IF (ASSOCIATED(destination)) THEN
            IF (destination%use_sp) CPABORT("only DP kind implemented")
         END IF

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
            IF ((send_rank .EQ. 0)) THEN
               ! need to use separate buffers here in case this is actually global rank 0
               src_block_tmp = (/send_dist%nrow_block, send_dist%ncol_block/)
               CALL mp_isend(src_block_tmp, 0, global_comm, send_req(1), tag=src_tag)
               CALL mp_isend(send_dist%context%num_pe, 0, global_comm, send_req(2), tag=src_tag)
            END IF
         END IF

         IF (ASSOCIATED(destination)) THEN
            IF ((recv_rank .EQ. 0)) THEN
               dest_block_tmp = (/recv_dist%nrow_block, recv_dist%ncol_block/)
               CALL mp_isend(dest_block_tmp, 0, global_comm, send_req(3), tag=dest_tag)
               CALL mp_isend(recv_dist%context%num_pe, 0, global_comm, send_req(4), tag=dest_tag)
            END IF
         END IF

         IF (global_rank == 0) THEN
            CALL mp_waitall(recv_req(1:4))
            ! Now we know the process decomposition, we can allocate the arrays to hold the blacs2mpi mapping
            ALLOCATE (src_blacs2mpi(0:src_num_pe(1) - 1, 0:src_num_pe(2) - 1), &
                      dest_blacs2mpi(0:dest_num_pe(1) - 1, 0:dest_num_pe(2) - 1) &
                      )
            CALL mp_irecv(src_blacs2mpi, mp_any_source, global_comm, recv_req(5), tag=src_tag)
            CALL mp_irecv(dest_blacs2mpi, mp_any_source, global_comm, recv_req(6), tag=dest_tag)
         END IF

         IF (ASSOCIATED(source)) THEN
            IF ((send_rank .EQ. 0)) THEN
               CALL mp_isend(send_dist%context%blacs2mpi(:, :), 0, global_comm, send_req(5), tag=src_tag)
            END IF
         END IF

         IF (ASSOCIATED(destination)) THEN
            IF ((recv_rank .EQ. 0)) THEN
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
         IF (global_rank .NE. 0) THEN
            ALLOCATE (src_blacs2mpi(0:src_num_pe(1) - 1, 0:src_num_pe(2) - 1), &
                      dest_blacs2mpi(0:dest_num_pe(1) - 1, 0:dest_num_pe(2) - 1) &
                      )
         END IF
         CALL mp_bcast(src_blacs2mpi, 0, global_comm)
         CALL mp_bcast(dest_blacs2mpi, 0, global_comm)
         info%src_blacs2mpi => src_blacs2mpi

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
                                  col_indices=recv_col_indices &
                                  )
            info%recv_col_indices => recv_col_indices
            info%recv_row_indices => recv_row_indices
            nrow_block_src = src_block(1)
            ncol_block_src = src_block(2)
            ALLOCATE (recv_count(0:send_size - 1), recv_disp(0:send_size - 1), recv_request(0:send_size - 1))

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
                  recv_count(src_blacs2mpi(p, q)) = nrows*ncols
               END DO
            END DO
            DEALLOCATE (src_p, src_q)

            ! Use one long buffer (and displacements into that buffer)
            !     this prevents the need for a rectangular array where not all elements will be populated
            ALLOCATE (recv_buf(SUM(recv_count(:))))
            info%recv_buf => recv_buf
            recv_disp(0) = 0
            DO i = 1, send_size - 1
               recv_disp(i) = recv_disp(i - 1) + recv_count(i - 1)
            END DO
            info%recv_disp => recv_disp

            ! Issue receive calls on ranks which expect data
            DO k = 0, send_size - 1
               IF (recv_count(k) .GT. 0) THEN
                  CALL mp_irecv(recv_buf(recv_disp(k) + 1:recv_disp(k) + recv_count(k)), &
                                source2global(k), global_comm, recv_request(k))
               ELSE
                  recv_request(k) = mp_request_null
               END IF
            END DO
            info%recv_request => recv_request
            DEALLOCATE (source2global)
         END IF ! ASSOCIATED(destination)

         ! SRC_1
         IF (ASSOCIATED(source)) THEN
            CALL cp_fm_struct_get(send_dist, row_indices=send_row_indices, &
                                  col_indices=send_col_indices &
                                  )
            nrow_block_dest = dest_block(1)
            ncol_block_dest = dest_block(2)
            ALLOCATE (send_count(0:recv_size - 1), send_disp(0:recv_size - 1), send_request(0:recv_size - 1))

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
            ALLOCATE (send_buf(SUM(send_count(:))))
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
                  send_buf(send_disp(mpi_rank) + send_count(mpi_rank)) = source%local_data(i, j)
               END DO
            END DO
            info%send_buf => send_buf

            ! For each non-zero send_count, call mpi_isend
            DO k = 0, recv_size - 1
               IF (send_count(k) .GT. 0) THEN
                  CALL mp_isend(send_buf(send_disp(k) + 1:send_disp(k) + send_count(k)), &
                                dest2global(k), global_comm, send_request(k))
               ELSE
                  send_request(k) = mp_request_null
               END IF
            END DO
            info%send_request => send_request
            DEALLOCATE (send_count, send_disp, dest2global)
         END IF ! ASSOCIATED(source)
         DEALLOCATE (dest_blacs2mpi)

      END IF !IF (.NOT. cp2k_is_parallel)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_start_copy_general

! **************************************************************************************************
!> \brief Completes the copy operation: wait for comms, unpack, clean up MPI state
!> \param destination output fm matrix
!> \param info all of the data that will be needed to complete the copy operation
! **************************************************************************************************
   SUBROUTINE cp_fm_finish_copy_general(destination, info)
      TYPE(cp_fm_type), POINTER                          :: destination
      TYPE(copy_info_type), INTENT(INOUT)                :: info

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_finish_copy_general', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, k, mpi_rank, send_size, &
                                                            src_p_i, src_q_j
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: recv_count
      INTEGER, DIMENSION(2)                              :: nblock_src, nlocal_recv, src_num_pe
      INTEGER, DIMENSION(:), POINTER                     :: recv_col_indices, recv_disp, &
                                                            recv_request, recv_row_indices, &
                                                            send_request
      INTEGER, DIMENSION(:, :), POINTER                  :: src_blacs2mpi
      REAL(KIND=dp), DIMENSION(:), POINTER               :: recv_buf, send_buf

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
         send_size = info%send_size
         nlocal_recv(1:2) = info%nlocal_recv(:)
         nblock_src(1:2) = info%nblock_src(:)
         src_num_pe(1:2) = info%src_num_pe(:)
         recv_col_indices => info%recv_col_indices
         recv_row_indices => info%recv_row_indices
         src_blacs2mpi => info%src_blacs2mpi
         recv_request => info%recv_request
         recv_buf => info%recv_buf
         recv_disp => info%recv_disp
         send_request => info%send_request
         send_buf => info%send_buf

         ! ... use the local variables to do the work
         ! DEST_2
         CALL mp_waitall(recv_request(:))
         ALLOCATE (recv_count(0:send_size - 1))
         ! Loop over the rmat, filling it in with data from the recv buffers
         ! (here the block sizes, num_pes refer to the distribution of the source matrix)
         recv_count(:) = 0
         DO j = 1, nlocal_recv(2)
            src_q_j = MOD(((recv_col_indices(j) - 1)/nblock_src(2)), src_num_pe(2))
            DO i = 1, nlocal_recv(1)
               src_p_i = MOD(((recv_row_indices(i) - 1)/nblock_src(1)), src_num_pe(1))
               mpi_rank = src_blacs2mpi(src_p_i, src_q_j)
               recv_count(mpi_rank) = recv_count(mpi_rank) + 1
               destination%local_data(i, j) = recv_buf(recv_disp(mpi_rank) + recv_count(mpi_rank))
            END DO
         END DO
         DEALLOCATE (recv_count, recv_disp, recv_request, recv_buf, src_blacs2mpi)
         ! Invalidate the stored state
         NULLIFY (info%recv_col_indices, &
                  info%recv_row_indices, &
                  info%recv_buf, &
                  info%recv_disp, &
                  info%recv_request, &
                  info%src_blacs2mpi)

      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_fm_finish_copy_general

! **************************************************************************************************
!> \brief Completes the copy operation: wait for comms clean up MPI state
!> \param source fm matrix
!> \param info all of the data that will be needed to complete the copy operation
! **************************************************************************************************
   SUBROUTINE cp_fm_cleanup_copy_general(source, info)
      TYPE(cp_fm_type), POINTER                          :: source
      TYPE(copy_info_type), INTENT(INOUT)                :: info

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_cleanup_copy_general', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: send_request
      INTEGER, DIMENSION(:, :), POINTER                  :: src_blacs2mpi
      REAL(KIND=dp), DIMENSION(:), POINTER               :: send_buf

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(source))
      MARK_USED(source)

      IF (.NOT. cp2k_is_parallel) THEN
         ! Don't do anything - no MPI state for the serial case
      ELSE
         send_request => info%send_request
         send_buf => info%send_buf
         src_blacs2mpi => info%src_blacs2mpi

         ! SRC_2
         ! If this process is also in the destination decomposition, this deallocate
         ! Was already done in cp_fm_finish_copy_general
         IF (ASSOCIATED(src_blacs2mpi)) THEN
            DEALLOCATE (src_blacs2mpi)
         ENDIF
         CALL mp_waitall(send_request(:))
         DEALLOCATE (send_request, send_buf)
         ! Invalidate the stored state
         NULLIFY (info%send_request, info%send_buf, info%src_blacs2mpi)

      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_fm_cleanup_copy_general

! **************************************************************************************************
!> \brief General copy of a submatrix of fm matrix to  a submatrix of another fm matrix.
!>        The two matrices can have different contexts.
!>
!>        Summary of distribution routines for dense matrices
!>        The following will copy A(iA:iA+M-1,jA:jA+N-1) to B(iB:iB+M-1,jB:jB+N-1):
!>
!>        call pdgemr2d(M,N,Aloc,iA,jA,descA,Bloc,iB,jB,descB,context)
!>
!>        A process that is not a part of the context of A should set descA(2)
!>        to -1, and similarly for B.
!>
!> \param source          input fm matrix
!> \param destination     output fm matrix
!> \param nrows           number of rows of sub matrix to be copied
!> \param ncols           number of cols of sub matrix to be copied
!> \param s_firstrow      starting global row index of sub matrix in source
!> \param s_firstcol      starting global col index of sub matrix in source
!> \param d_firstrow      starting global row index of sub matrix in destination
!> \param d_firstcol      starting global col index of sub matrix in destination
!> \param global_context  process grid that covers all parts of either A or B.
! **************************************************************************************************
   SUBROUTINE cp_fm_to_fm_submat_general(source, &
                                         destination, &
                                         nrows, &
                                         ncols, &
                                         s_firstrow, &
                                         s_firstcol, &
                                         d_firstrow, &
                                         d_firstcol, &
                                         global_context)

      TYPE(cp_fm_type), POINTER                :: source, destination
      INTEGER, INTENT(IN)                      :: nrows, ncols, &
                                                  s_firstrow, s_firstcol, &
                                                  d_firstrow, d_firstcol, &
                                                  global_context

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_to_fm_submat_general', &
                                     routineP = moduleN//':'//routineN

      LOGICAL                                  :: debug
      INTEGER                                  :: handle
#if defined(__SCALAPACK)
      INTEGER, DIMENSION(9)                    :: desca, descb
      REAL(KIND=dp), DIMENSION(1, 1), TARGET    :: dummy
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: smat, dmat
#endif

      CALL timeset(routineN, handle)

      debug = debug_this_module

      IF (.NOT. cp2k_is_parallel) THEN
         CALL cp_fm_to_fm_submat(source, &
                                 destination, &
                                 nrows, &
                                 ncols, &
                                 s_firstrow, &
                                 s_firstcol, &
                                 d_firstrow, &
                                 d_firstcol)
      ELSE
#ifdef __SCALAPACK
         NULLIFY (smat, dmat)
         ! check whether source is available on this process
         IF (ASSOCIATED(source)) THEN
            desca = source%matrix_struct%descriptor
            IF (source%use_sp) CPABORT("only DP kind implemented")
            IF (nrows .GT. source%matrix_struct%nrow_global) &
               CPABORT("nrows is greater than nrow_global of source")
            IF (ncols .GT. source%matrix_struct%ncol_global) &
               CPABORT("ncols is greater than ncol_global of source")
            smat => source%local_data
         ELSE
            desca = -1
            smat => dummy
         END IF
         ! check destination is available on this process
         IF (ASSOCIATED(destination)) THEN
            descb = destination%matrix_struct%descriptor
            IF (destination%use_sp) CPABORT("only DP kind implemented")
            IF (nrows .GT. destination%matrix_struct%nrow_global) &
               CPABORT("nrows is greater than nrow_global of destination")
            IF (ncols .GT. destination%matrix_struct%ncol_global) &
               CPABORT("ncols is greater than ncol_global of destination")
            dmat => destination%local_data
         ELSE
            descb = -1
            dmat => dummy
         END IF
         ! do copy

         CALL pdgemr2d(nrows, &
                       ncols, &
                       smat, &
                       s_firstrow, &
                       s_firstcol, &
                       desca, &
                       dmat, &
                       d_firstrow, &
                       d_firstcol, &
                       descb, &
                       global_context)
#else
         MARK_USED(global_context)
         CPABORT("this subroutine only supports SCALAPACK")
#endif
      END IF

      CALL timestop(handle)

   END SUBROUTINE cp_fm_to_fm_submat_general

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param irow_global ...
!> \param icol_global ...
!> \param alpha ...
! **************************************************************************************************
   SUBROUTINE cp_fm_add_to_element(matrix, irow_global, icol_global, alpha)

      ! Add alpha to the matrix element specified by the global indices
      ! irow_global and icol_global

      ! - Creation (05.05.06,MK)

      TYPE(cp_fm_type), POINTER                :: matrix
      INTEGER, INTENT(IN)                      :: irow_global, icol_global
      REAL(KIND=dp), INTENT(IN)                :: alpha

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_add_to_element', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: mypcol, myprow, npcol, nprow
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
      TYPE(cp_blacs_env_type), POINTER         :: context
#if defined(__SCALAPACK)
      INTEGER                                  :: icol_local, ipcol, iprow, &
                                                  irow_local
      INTEGER, DIMENSION(9)                    :: desca
#endif

      CPASSERT(ASSOCIATED(matrix))

      context => matrix%matrix_struct%context

      myprow = context%mepos(1)
      mypcol = context%mepos(2)

      nprow = context%num_pe(1)
      npcol = context%num_pe(2)

      a => matrix%local_data

#if defined(__SCALAPACK)

      desca(:) = matrix%matrix_struct%descriptor(:)

      CALL infog2l(irow_global, icol_global, desca, nprow, npcol, myprow, mypcol, &
                   irow_local, icol_local, iprow, ipcol)

      IF ((iprow == myprow) .AND. (ipcol == mypcol)) THEN
         a(irow_local, icol_local) = a(irow_local, icol_local) + alpha
      END IF

#else

      a(irow_global, icol_global) = a(irow_global, icol_global) + alpha

#endif

   END SUBROUTINE cp_fm_add_to_element

! **************************************************************************************************
!> \brief ...
!> \param fm ...
!> \param unit ...
! **************************************************************************************************
   SUBROUTINE cp_fm_write_unformatted(fm, unit)
      TYPE(cp_fm_type), POINTER                :: fm
      INTEGER                                  :: unit

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_write_unformatted', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: handle, j, max_block, &
                                                  ncol_global, nrow_global
      TYPE(cp_para_env_type), POINTER          :: para_env
#if defined(__SCALAPACK)
      INTEGER                                  :: i, i_block, icol_local, &
                                                  ictxt_loc, in, info, ipcol, &
                                                  iprow, irow_local, &
                                                  mepos, mypcol, myprow, &
                                                  npcol, nprow, num_pe, rb, tag
      INTEGER, DIMENSION(9)                    :: desc
      REAL(KIND=dp), DIMENSION(:), POINTER     :: vecbuf
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: newdat
      INTEGER, EXTERNAL                        :: numroc
#endif

      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(fm, nrow_global=nrow_global, ncol_global=ncol_global, ncol_block=max_block, &
                          para_env=para_env)

#if defined(__SCALAPACK)
      num_pe = para_env%num_pe
      mepos = para_env%mepos
      rb = nrow_global
      tag = 0
      ! get a new context
      ictxt_loc = para_env%group
      CALL cp_blacs_gridinit(ictxt_loc, 'R', 1, num_pe)
      CALL cp_blacs_gridinfo(ictxt_loc, nprow, npcol, myprow, mypcol)
      CALL descinit(desc, nrow_global, ncol_global, rb, max_block, 0, 0, ictxt_loc, nrow_global, info)
      CPASSERT(info == 0)
      in = numroc(ncol_global, max_block, mypcol, 0, npcol)

      ALLOCATE (newdat(nrow_global, MAX(1, in)))

      ! do the actual scalapack to cols reordering
      CALL pdgemr2d(nrow_global, ncol_global, fm%local_data(1, 1), 1, 1, &
                    fm%matrix_struct%descriptor, &
                    newdat(1, 1), 1, 1, desc, ictxt_loc)

      ALLOCATE (vecbuf(nrow_global*max_block))
      vecbuf = HUGE(1.0_dp) ! init for valgrind

      DO i = 1, ncol_global, MAX(max_block, 1)
         i_block = MIN(max_block, ncol_global - i + 1)
         CALL infog2l(1, i, desc, nprow, npcol, myprow, mypcol, &
                      irow_local, icol_local, iprow, ipcol)
         IF (ipcol == mypcol) THEN
            DO j = 1, i_block
               vecbuf((j - 1)*nrow_global + 1:nrow_global*j) = newdat(:, icol_local + j - 1)
            END DO
         END IF

         IF (ipcol == 0) THEN
            ! do nothing
         ELSE
            IF (ipcol == mypcol) THEN
               CALL mp_send(vecbuf(:), 0, tag, para_env%group)
            ENDIF
            IF (mypcol == 0) THEN
               CALL mp_recv(vecbuf(:), ipcol, tag, para_env%group)
            ENDIF
         ENDIF

         IF (unit > 0) THEN
            DO j = 1, i_block
               WRITE (unit) vecbuf((j - 1)*nrow_global + 1:nrow_global*j)
            END DO
         END IF

      END DO
      DEALLOCATE (vecbuf)

      CALL cp_blacs_gridexit(ictxt_loc)

      DEALLOCATE (newdat)

#else

      IF (unit > 0) THEN
         DO j = 1, ncol_global
            WRITE (unit) fm%local_data(:, j)
         END DO
      END IF

#endif
      CALL timestop(handle)

   END SUBROUTINE cp_fm_write_unformatted

! **************************************************************************************************
!> \brief Write out a full matrix in plain text.
!> \param fm the matrix to be outputted
!> \param unit the unit number for I/O
!> \param header optional header
!> \param value_format ...
! **************************************************************************************************
   SUBROUTINE cp_fm_write_formatted(fm, unit, header, value_format)
      TYPE(cp_fm_type), POINTER                :: fm
      INTEGER                                  :: unit
      CHARACTER(LEN=*), OPTIONAL               :: header, value_format

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_write_formatted', &
                                     routineP = moduleN//':'//routineN

      CHARACTER(LEN=21)                        :: my_value_format
      INTEGER                                  :: handle, i, j, max_block, &
                                                  ncol_global, nrow_global, &
                                                  nrow_block
      TYPE(cp_para_env_type), POINTER          :: para_env
#if defined(__SCALAPACK)
      INTEGER                                  :: i_block, icol_local, &
                                                  ictxt_loc, in, info, ipcol, &
                                                  iprow, irow_local, &
                                                  mepos, mypcol, myprow, &
                                                  npcol, nprow, num_pe, rb, tag, k, &
                                                  icol, irow
      INTEGER, DIMENSION(9)                    :: desc
      REAL(KIND=dp), DIMENSION(:), POINTER     :: vecbuf
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: newdat
      INTEGER, EXTERNAL                        :: numroc
#endif

      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(fm, nrow_global=nrow_global, ncol_global=ncol_global, ncol_block=max_block, &
                          nrow_block=nrow_block, para_env=para_env)

      IF (PRESENT(value_format)) THEN
         CPASSERT(LEN_TRIM(ADJUSTL(value_format)) < 11)
         my_value_format = "(I10, I10, "//TRIM(ADJUSTL(value_format))//")"
      ELSE
         my_value_format = "(I10, I10, ES24.12)"
      END IF

      IF (unit > 0) THEN
         IF (PRESENT(header)) WRITE (unit, *) header
         WRITE (unit, "(A2, A8, A10, A24)") "#", "Row", "Column", ADJUSTL("Value")
      END IF

#if defined(__SCALAPACK)
      num_pe = para_env%num_pe
      mepos = para_env%mepos
      rb = nrow_global
      tag = 0
      ! get a new context
      ictxt_loc = para_env%group
      CALL cp_blacs_gridinit(ictxt_loc, 'R', 1, num_pe)
      CALL cp_blacs_gridinfo(ictxt_loc, nprow, npcol, myprow, mypcol)
      CALL descinit(desc, nrow_global, ncol_global, rb, max_block, 0, 0, ictxt_loc, nrow_global, info)
      CPASSERT(info == 0)
      in = numroc(ncol_global, max_block, mypcol, 0, npcol)

      ALLOCATE (newdat(nrow_global, MAX(1, in)))

      ! do the actual scalapack to cols reordering
      CALL pdgemr2d(nrow_global, ncol_global, fm%local_data(1, 1), 1, 1, &
                    fm%matrix_struct%descriptor, &
                    newdat(1, 1), 1, 1, desc, ictxt_loc)

      ALLOCATE (vecbuf(nrow_global*max_block))
      vecbuf = HUGE(1.0_dp) ! init for valgrind
      irow = 1
      icol = 1

      DO i = 1, ncol_global, MAX(max_block, 1)
         i_block = MIN(max_block, ncol_global - i + 1)
         CALL infog2l(1, i, desc, nprow, npcol, myprow, mypcol, &
                      irow_local, icol_local, iprow, ipcol)
         IF (ipcol == mypcol) THEN
            DO j = 1, i_block
               vecbuf((j - 1)*nrow_global + 1:nrow_global*j) = newdat(:, icol_local + j - 1)
            END DO
         END IF

         IF (ipcol == 0) THEN
            ! do nothing
         ELSE
            IF (ipcol == mypcol) THEN
               CALL mp_send(vecbuf(:), 0, tag, para_env%group)
            ENDIF
            IF (mypcol == 0) THEN
               CALL mp_recv(vecbuf(:), ipcol, tag, para_env%group)
            ENDIF
         ENDIF

         IF (unit > 0) THEN
            DO j = 1, i_block
               DO k = (j - 1)*nrow_global + 1, nrow_global*j
                  WRITE (UNIT=unit, FMT=my_value_format) irow, icol, vecbuf(k)
                  irow = irow + 1
                  IF (irow > nrow_global) THEN
                     irow = 1
                     icol = icol + 1
                  END IF
               END DO
            END DO
         END IF

      END DO
      DEALLOCATE (vecbuf)

      CALL cp_blacs_gridexit(ictxt_loc)

      DEALLOCATE (newdat)

#else

      IF (unit > 0) THEN
         DO j = 1, ncol_global
            DO i = 1, nrow_global
               WRITE (UNIT=unit, FMT=my_value_format) i, j, fm%local_data(i, j)
            END DO
         END DO
      END IF

#endif
      CALL timestop(handle)

   END SUBROUTINE cp_fm_write_formatted

! **************************************************************************************************
!> \brief ...
!> \param fm ...
!> \param unit ...
! **************************************************************************************************
   SUBROUTINE cp_fm_read_unformatted(fm, unit)
      TYPE(cp_fm_type), POINTER                :: fm
      INTEGER                                  :: unit

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_read_unformatted', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: handle, j, max_block, &
                                                  ncol_global, nrow_global
      TYPE(cp_para_env_type), POINTER          :: para_env
#if defined(__SCALAPACK)
      INTEGER                                  :: k, n_cols
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuf
#endif

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(fm, nrow_global=nrow_global, ncol_global=ncol_global, ncol_block=max_block, &
                          para_env=para_env)

#if defined(__SCALAPACK)

      ! the parallel case could be made more efficient (see cp_fm_write_unformatted)

      ALLOCATE (vecbuf(nrow_global, max_block))

      DO j = 1, ncol_global, max_block

         n_cols = MIN(max_block, ncol_global - j + 1)
         IF (para_env%mepos == 0) THEN
            DO k = 1, n_cols
               READ (unit) vecbuf(:, k)
            ENDDO
         ENDIF
         CALL mp_bcast(vecbuf, 0, para_env%group)
         CALL cp_fm_set_submatrix(fm, vecbuf, start_row=1, start_col=j, n_cols=n_cols)

      ENDDO

      DEALLOCATE (vecbuf)

#else

      DO j = 1, ncol_global
         READ (unit) fm%local_data(:, j)
      END DO

#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_read_unformatted

! **************************************************************************************************
!> \brief wrapper to scalapack function INDXG2P that computes the process
!>         coordinate which possesses the entry of a distributed matrix specified
!>         by a global index INDXGLOB.
!>
!>  Arguments
!>  =========
!>
!>  INDXGLOB  (global input) INTEGER
!>            The global index of the element.
!>
!>  NB        (global input) INTEGER
!>            Block size, size of the blocks the distributed matrix is
!>            split into.
!>
!>  IPROC     (local dummy) INTEGER
!>            Dummy argument in this case in order to unify the calling
!>            sequence of the tool-routines.
!>
!>  ISRCPROC  (global input) INTEGER
!>            The coordinate of the process that possesses the first
!>            row/column of the distributed matrix.
!>
!>  NPROCS    (global input) INTEGER
!>            The total number processes over which the matrix is
!>            distributed.
!>
!> \param INDXGLOB ...
!> \param NB ...
!> \param IPROC ...
!> \param ISRCPROC ...
!> \param NPROCS ...
!> \return ...
!> \author Mauro Del Ben [MDB] - 12.2012
! **************************************************************************************************
   FUNCTION cp_fm_indxg2p(INDXGLOB, NB, IPROC, ISRCPROC, NPROCS) RESULT(G2P)
      INTEGER, INTENT(IN)                      :: INDXGLOB, NB, IPROC, ISRCPROC, NPROCS
      INTEGER                                  :: G2P

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_indxg2p', &
                                     routineP = moduleN//':'//routineN

#if defined(__SCALAPACK)
      INTEGER, EXTERNAL :: indxg2p
#endif

#if defined(__SCALAPACK)

      G2P = indxg2p(INDXGLOB, NB, IPROC, ISRCPROC, NPROCS)

#else
      MARK_USED(indxglob)
      MARK_USED(iproc)
      MARK_USED(isrcproc)
      MARK_USED(nb)
      MARK_USED(nprocs)

      G2P = 0

#endif

   END FUNCTION cp_fm_indxg2p

! **************************************************************************************************
!> \brief wrapper to scalapack function INDXG2L that computes the local index
!>         of a distributed matrix entry pointed to by the global index INDXGLOB.
!>
!>  Arguments
!>  =========
!>
!>  INDXGLOB  (global input) INTEGER
!>            The global index of the distributed matrix entry.
!>
!>  NB        (global input) INTEGER
!>            Block size, size of the blocks the distributed matrix is
!>            split into.
!>
!>  IPROC     (local dummy) INTEGER
!>            Dummy argument in this case in order to unify the calling
!>            sequence of the tool-routines.
!>
!>  ISRCPROC  (local dummy) INTEGER
!>            Dummy argument in this case in order to unify the calling
!>            sequence of the tool-routines.
!>
!>  NPROCS    (global input) INTEGER
!>            The total number processes over which the distributed
!>            matrix is distributed.
!>
!> \param INDXGLOB ...
!> \param NB ...
!> \param IPROC ...
!> \param ISRCPROC ...
!> \param NPROCS ...
!> \return ...
!> \author Mauro Del Ben [MDB] - 12.2012
! **************************************************************************************************
   FUNCTION cp_fm_indxg2l(INDXGLOB, NB, IPROC, ISRCPROC, NPROCS) RESULT(G2L)
      INTEGER, INTENT(IN)                      :: INDXGLOB, NB, IPROC, ISRCPROC, NPROCS
      INTEGER                                  :: G2L

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_indxg2l', &
                                     routineP = moduleN//':'//routineN

#if defined(__SCALAPACK)
      INTEGER, EXTERNAL :: indxg2l
#endif

#if defined(__SCALAPACK)

      G2L = indxg2l(INDXGLOB, NB, IPROC, ISRCPROC, NPROCS)

#else
      MARK_USED(iproc)
      MARK_USED(isrcproc)
      MARK_USED(nb)
      MARK_USED(nprocs)
      MARK_USED(indxglob)

      G2L = INDXGLOB

#endif

   END FUNCTION cp_fm_indxg2l

! **************************************************************************************************
!> \brief wrapper to scalapack function INDXL2G that computes the global index
!>         of a distributed matrix entry pointed to by the local index INDXLOC
!>         of the process indicated by IPROC.
!>
!>  Arguments
!>  =========
!>
!>  INDXLOC   (global input) INTEGER
!>            The local index of the distributed matrix entry.
!>
!>  NB        (global input) INTEGER
!>            Block size, size of the blocks the distributed matrix is
!>            split into.
!>
!>  IPROC     (local input) INTEGER
!>            The coordinate of the process whose local array row or
!>            column is to be determined.
!>
!>  ISRCPROC  (global input) INTEGER
!>            The coordinate of the process that possesses the first
!>            row/column of the distributed matrix.
!>
!>  NPROCS    (global input) INTEGER
!>            The total number processes over which the distributed
!>            matrix is distributed.
!>
!> \param INDXLOC ...
!> \param NB ...
!> \param IPROC ...
!> \param ISRCPROC ...
!> \param NPROCS ...
!> \return ...
!> \author Mauro Del Ben [MDB] - 12.2012
! **************************************************************************************************
   FUNCTION cp_fm_indxl2g(INDXLOC, NB, IPROC, ISRCPROC, NPROCS) RESULT(L2G)
      INTEGER, INTENT(IN)                      :: INDXLOC, NB, IPROC, ISRCPROC, NPROCS
      INTEGER                                  :: L2G

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_indxl2g', &
                                     routineP = moduleN//':'//routineN

#if defined(__SCALAPACK)
      INTEGER, EXTERNAL :: indxl2g
#endif

#if defined(__SCALAPACK)

      L2G = indxl2g(INDXLOC, NB, IPROC, ISRCPROC, NPROCS)

#else
      MARK_USED(iproc)
      MARK_USED(isrcproc)
      MARK_USED(nb)
      MARK_USED(nprocs)

      L2G = INDXLOC

#endif

   END FUNCTION cp_fm_indxl2g

! **************************************************************************************************
!> \brief ...
!> \param mult_type ...
! **************************************************************************************************
   SUBROUTINE cp_fm_setup(mult_type)
      INTEGER                                            :: mult_type

      mm_type = mult_type
   END SUBROUTINE cp_fm_setup

! **************************************************************************************************
!> \brief ...
!> \return ...
! **************************************************************************************************
   FUNCTION cp_fm_get_mm_type() RESULT(res)
      INTEGER                                            :: res

      res = mm_type
   END FUNCTION

END MODULE cp_fm_types
