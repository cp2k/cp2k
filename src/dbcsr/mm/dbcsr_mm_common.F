!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Common variables and routines for the dbcsr matrix-matrix multiplication algorithms.
!>
!> \author  Alfio Lazzaro
!>
!> <b>Modification history:</b>
!>  - 2016-08    Code organization (Alfio Lazzaro).
! **************************************************************************************************
MODULE dbcsr_mm_common
   USE acc_stream,                      ONLY: acc_stream_type
   USE array_types,                     ONLY: array_hold
   USE dbcsr_data_methods,              ONLY: dbcsr_data_ensure_size
   USE dbcsr_methods,                   ONLY: dbcsr_get_index_memory_type,&
                                              dbcsr_nfullcols_local,&
                                              dbcsr_nfullrows_local
   USE dbcsr_mm_multrec,                ONLY: dbcsr_mm_multrec_type
   USE dbcsr_ptr_util,                  ONLY: ensure_array_size
   USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                              dbcsr_memtype_type,&
                                              dbcsr_mpi_size_limits,&
                                              dbcsr_mpi_statistics_type,&
                                              dbcsr_type
   USE dbcsr_work_operations,           ONLY: dbcsr_create
   USE kinds,                           ONLY: int_4,&
                                              int_8,&
                                              real_8
#include "../../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_common'

   TYPE dbcsr_memtype_type_p
      TYPE(dbcsr_memtype_type), POINTER :: p => Null()
      ! ensure that array-elements are on different cache lines
      INTEGER(kind=int_4), DIMENSION(64)    :: padding
   END TYPE dbcsr_memtype_type_p

   TYPE(dbcsr_memtype_type_p), DIMENSION(:), POINTER, SAVE :: memtype_product_wm => Null()

   TYPE(dbcsr_mpi_statistics_type), SAVE :: dbcsr_mpi_statistics
   INTEGER, SAVE :: num_multiplications = 0
   REAL, SAVE :: max_memory = 0

   TYPE(dbcsr_memtype_type), SAVE  :: memtype_abpanel_1, memtype_abpanel_2, &
                                      memtype_trsbuffer_1, memtype_trsbuffer_2, &
                                      memtype_mpi_buffer, memtype_mpi_product
   TYPE(acc_stream_type), SAVE     :: stream_1, stream_2
   ! ab-panels and streams are shared between all threads

   TYPE dbcsr_mm_multrec_type_p
      TYPE(dbcsr_mm_multrec_type), POINTER :: p => Null()
      ! ensure that array-elements are on different cache lines
      INTEGER(kind=int_4), DIMENSION(64)       :: padding
   END TYPE dbcsr_mm_multrec_type_p

   PUBLIC :: memtype_product_wm
   PUBLIC :: dbcsr_mpi_statistics, num_multiplications
   PUBLIC :: max_memory

   PUBLIC :: memtype_abpanel_1, memtype_abpanel_2, &
             memtype_trsbuffer_1, memtype_trsbuffer_2, &
             memtype_mpi_buffer, memtype_mpi_product
   PUBLIC :: stream_1, stream_2

   PUBLIC :: dbcsr_mm_multrec_type_p
   PUBLIC :: count_mpi_statistics
   PUBLIC :: setup_buffer_matrix
   PUBLIC :: rec_sort_index
   PUBLIC :: enumerate_blk_sizes

   PUBLIC :: product_matrix_size_guess

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mpi_statistics ...
!> \param data_size ...
!> \param size_breakdown ...
!> \param element_size_bytes ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE count_mpi_statistics(mpi_statistics, data_size, &
                                   size_breakdown, element_size_bytes)
      REAL, DIMENSION(:), INTENT(INOUT)                  :: mpi_statistics
      INTEGER, INTENT(IN)                                :: data_size
      INTEGER(KIND=int_8), DIMENSION(:, :), &
         INTENT(INOUT), OPTIONAL                         :: size_breakdown
      INTEGER, INTENT(IN), OPTIONAL                      :: element_size_bytes

      INTEGER                                            :: ilimit, nlimits
      INTEGER(KIND=int_8)                                :: data_size_bytes, llimit

      mpi_statistics(1) = mpi_statistics(1)+REAL(data_size)
      mpi_statistics(2) = MIN(mpi_statistics(2), REAL(data_size))
      mpi_statistics(3) = MAX(mpi_statistics(3), REAL(data_size))
      IF (PRESENT(size_breakdown)) THEN
         data_size_bytes = data_size
         ! change in bytes
         IF (PRESENT(element_size_bytes)) data_size_bytes = data_size_bytes*element_size_bytes
         nlimits = SIZE(dbcsr_mpi_size_limits)
         ! check for oversize messages
         IF (data_size_bytes .GT. dbcsr_mpi_size_limits(nlimits)) THEN
            size_breakdown(nlimits+1, 1) = size_breakdown(nlimits+1, 1)+1
            size_breakdown(nlimits+1, 2) = size_breakdown(nlimits+1, 2)+data_size_bytes
            RETURN
         ENDIF
         llimit = 0
         DO ilimit = 1, nlimits
            IF (data_size_bytes .GE. llimit .AND. data_size_bytes .LE. dbcsr_mpi_size_limits(ilimit)) THEN
               size_breakdown(ilimit, 1) = size_breakdown(ilimit, 1)+1
               size_breakdown(ilimit, 2) = size_breakdown(ilimit, 2)+data_size_bytes
               RETURN
            ENDIF
            llimit = dbcsr_mpi_size_limits(ilimit)
         ENDDO
      ENDIF
   END SUBROUTINE count_mpi_statistics

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param source_matrix ...
!> \param index_size ...
!> \param data_size ...
!> \param data_buffer ...
!> \param data_memory_type ...
! **************************************************************************************************
   SUBROUTINE setup_buffer_matrix(matrix, source_matrix, &
                                  index_size, data_size, data_buffer, data_memory_type)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      TYPE(dbcsr_type), INTENT(IN)                       :: source_matrix
      INTEGER, INTENT(IN), OPTIONAL                      :: index_size, data_size
      TYPE(dbcsr_data_obj), INTENT(IN), OPTIONAL         :: data_buffer
      TYPE(dbcsr_memtype_type), INTENT(IN), OPTIONAL     :: data_memory_type

      matrix = dbcsr_type()
      CALL dbcsr_create(matrix, &
                        template=source_matrix, &
                        name=TRIM("Buffer of "//TRIM(source_matrix%name)), &
                        nze=data_size, &
                        data_buffer=data_buffer, &
                        data_memory_type=data_memory_type, &
                        index_memory_type=memtype_mpi_buffer)
      IF (PRESENT(data_size)) THEN
         CALL dbcsr_data_ensure_size( &
            matrix%data_area, &
            data_size, nocopy=.TRUE.)
      ENDIF
      IF (PRESENT(index_size)) THEN
         CALL ensure_array_size( &
            matrix%index, &
            ub=index_size, nocopy=.TRUE., &
            memory_type=dbcsr_get_index_memory_type(matrix))
      ENDIF
      matrix%negate_real = source_matrix%negate_real
      matrix%negate_imaginary = source_matrix%negate_imaginary
      matrix%local_indexing = source_matrix%local_indexing
      matrix%list_indexing = source_matrix%list_indexing
      !
      IF (source_matrix%has_local_rows) THEN
         matrix%local_rows = source_matrix%local_rows
         CALL array_hold(matrix%local_rows)
         matrix%has_local_rows = .TRUE.
      ENDIF
      IF (source_matrix%has_global_rows) THEN
         matrix%global_rows = source_matrix%global_rows
         CALL array_hold(matrix%global_rows)
         matrix%has_global_rows = .TRUE.
      ENDIF
      IF (source_matrix%has_local_cols) THEN
         matrix%local_cols = source_matrix%local_cols
         CALL array_hold(matrix%local_cols)
         matrix%has_local_cols = .TRUE.
      ENDIF
      IF (source_matrix%has_global_cols) THEN
         matrix%global_cols = source_matrix%global_cols
         CALL array_hold(matrix%global_cols)
         matrix%has_global_cols = .TRUE.
      ENDIF

   END SUBROUTINE setup_buffer_matrix

! **************************************************************************************************
!> \brief Sorts index for recursing.
!> \param mi ...
!> \param mf ...
!> \param ni ...
!> \param nf ...
!> \param nele ...
!> \param a ...
!> \param d ...
!> \par History
!> - 2011-02-17 [UB] modified for use in DBCSR; reduced memory usage.
!> \author JV
!> \note Always cut longest first. On a tie cut N
! **************************************************************************************************
   RECURSIVE SUBROUTINE rec_sort_index(mi, mf, ni, nf, nele, a, d)
      INTEGER, INTENT(IN)                                :: mi, mf, ni, nf, nele
      INTEGER, DIMENSION(3, 1:nele), INTENT(inout)       :: a
      INTEGER, INTENT(IN)                                :: d

      LOGICAL, PARAMETER                                 :: dbg = .FALSE.

      INTEGER                                            :: half, M, N, nlow
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: tmp

!   ---------------------------------------------------------------------------

      IF (dbg) THEN
         WRITE (*, *) " rs", mi, mf, "/", ni, nf, "=>", nele, d
         WRITE (*, '(3(1X,I7))') a(:, 1:nele)
      ENDIF
      IF (dbg) THEN
         IF (d .GT. 20) THEN
            WRITE (*, *) a(1, -d*1000)
         ENDIF
      ENDIF
      ALLOCATE (tmp(3, nele))
      M = mf-mi+1
      N = nf-ni+1
      IF (M > N) THEN
         half = M/2
         CALL rec_split(nele, a, tmp, 1, nlow, mi, half)
         a = tmp
         DEALLOCATE (tmp)
         IF (nlow .GT. 1) THEN
            CALL rec_sort_index(mi, mi+half-1, ni, nf, nlow, a(:, 1:nlow), d+1)
         ENDIF
         IF (nele-nlow .GT. 1) THEN
            CALL rec_sort_index(mi+half, mf, ni, nf, nele-nlow, a(:, nlow+1:nele), d+1)
         ENDIF
      ELSE
         half = N/2
         CALL rec_split(nele, a, tmp, 2, nlow, ni, half)
         a = tmp
         DEALLOCATE (tmp)
         IF (nlow .GT. 1) THEN
            CALL rec_sort_index(mi, mf, ni, ni+half-1, nlow, a(:, 1:nlow), d+1)
         ENDIF
         IF (nele-nlow .GT. 1) THEN
            CALL rec_sort_index(mi, mf, ni+half, nf, nele-nlow, a(:, nlow+1:nele), d+1)
         ENDIF
      ENDIF
   END SUBROUTINE rec_sort_index

! **************************************************************************************************
!> \brief ...
!> \param nele ...
!> \param a ...
!> \param split ...
!> \param row_or_col ...
!> \param nlow ...
!> \param mi ...
!> \param half ...
! **************************************************************************************************
   SUBROUTINE rec_split(nele, a, split, row_or_col, nlow, mi, half)
      INTEGER, INTENT(IN)                                :: nele
      INTEGER, DIMENSION(3, nele), INTENT(IN)            :: a
      INTEGER, DIMENSION(3, nele), INTENT(OUT)           :: split
      INTEGER, INTENT(IN)                                :: row_or_col
      INTEGER, INTENT(OUT)                               :: nlow
      INTEGER, INTENT(IN)                                :: mi, half

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rec_split', routineP = moduleN//':'//routineN

      INTEGER                                            :: el, half_m, p_high, p_low

      half_m = mi+half-1
      p_low = 1
      p_high = nele
      DO el = 1, nele
         IF (a(row_or_col, el) <= half_m) THEN
            split(1:3, p_low) = a(1:3, el)
            p_low = p_low+1
         ELSE
            split(1:3, p_high) = a(1:3, el)
            p_high = p_high-1
         ENDIF
      ENDDO
      nlow = p_low-1
      CPASSERT(p_high .EQ. nlow)

   END SUBROUTINE rec_split

! **************************************************************************************************
!> \brief Enumerate all occuring blocksizes
!> \param blk_sizes ...
!> \param enum ...
!> \param rev_enum ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE enumerate_blk_sizes(blk_sizes, enum, rev_enum)
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes, enum, rev_enum

      CHARACTER(len=*), PARAMETER :: routineN = 'enumerate_blk_sizes', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, n

      CALL timeset(routineN, handle)

      n = MAXVAL(blk_sizes)
      ALLOCATE (enum(0:n))
      enum(:) = 0

      DO i = 1, SIZE(blk_sizes)
         enum(blk_sizes(i)) = 1
      ENDDO

      n = SUM(enum)
      ALLOCATE (rev_enum(n))

      n = 0
      DO i = 0, SIZE(enum)-1
         IF (enum(i) > 0) THEN
            n = n+1
            enum(i) = n
            rev_enum(n) = i
         END IF
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE enumerate_blk_sizes

! **************************************************************************************************
!> \brief Guess the size of the product matrix from the A and B sparsities
!> \param matrix_left ...
!> \param matrix_right ...
!> \param product_matrix ...
!> \param left_data_size ...
!> \param right_data_size ...
!> \param left_col_nimages ...
!> \param right_row_nimages ...
!> \param nthreads ...
!> \retval size_guess ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   FUNCTION product_matrix_size_guess(matrix_left, matrix_right, product_matrix, &
                                      left_data_size, right_data_size, &
                                      left_col_nimages, right_row_nimages, &
                                      nthreads) RESULT(size_guess)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_left, matrix_right, product_matrix
      INTEGER, INTENT(IN)                                :: left_data_size, right_data_size, &
                                                            left_col_nimages, right_row_nimages, &
                                                            nthreads
      INTEGER                                            :: size_guess

      INTEGER(KIND=int_8)                                :: size8
      REAL(kind=real_8)                                  :: factor, fill_guess, left_fill, right_fill

      ! First we calculate the sparsities
      size8 = INT(dbcsr_nfullrows_local(matrix_left), KIND=int_8)* &
              INT(dbcsr_nfullcols_local(matrix_left), KIND=int_8)
      size8 = MAX(1_int_8, size8)
      left_fill = (REAL(left_data_size, KIND=real_8)*REAL(left_col_nimages, KIND=real_8))/REAL(size8, KIND=real_8)
      size8 = INT(dbcsr_nfullrows_local(matrix_right), KIND=int_8)* &
              INT(dbcsr_nfullcols_local(matrix_right), KIND=int_8)
      size8 = MAX(1_int_8, size8)
      right_fill = (REAL(right_data_size, KIND=real_8)*REAL(right_row_nimages, KIND=real_8))/REAL(size8, KIND=real_8)
      size8 = INT(dbcsr_nfullrows_local(product_matrix), KIND=int_8)* &
              INT(dbcsr_nfullcols_local(product_matrix), KIND=int_8)
      size8 = MAX(1_int_8, size8)
!     factor = 7.0 ! Old guess
      factor = 2.4 ! New guess
      fill_guess = factor*MAX(left_fill, right_fill)
      fill_guess = MIN(1.0_real_8, MAX(0.0_real_8, fill_guess))
      IF (nthreads .GT. 1) THEN
         fill_guess = fill_guess*3.0_real_8/REAL(2*nthreads, KIND=real_8)
      ENDIF
      size_guess = INT(REAL(size8, KIND=real_8)*fill_guess, KIND=int_4)
   END FUNCTION product_matrix_size_guess

END MODULE dbcsr_mm_common
