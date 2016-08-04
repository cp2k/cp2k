!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
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
   USE acc_event,                       ONLY: acc_event_record,&
                                              acc_event_synchronize,&
                                              acc_stream_wait_event
   USE acc_stream,                      ONLY: acc_stream_type
   USE array_types,                     ONLY: array_data,&
                                              array_hold
   USE dbcsr_acc_operations,            ONLY: dbcsr_acc_transpose
   USE dbcsr_data_methods,              ONLY: dbcsr_data_ensure_size,&
                                              dbcsr_data_host2dev
   USE dbcsr_methods,                   ONLY: dbcsr_get_index_memory_type,&
                                              dbcsr_init
   USE dbcsr_mm_multrec,                ONLY: dbcsr_mm_multrec_type
   USE dbcsr_ptr_util,                  ONLY: ensure_array_size
   USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                              dbcsr_memtype_type,&
                                              dbcsr_mpi_size_limits,&
                                              dbcsr_mpi_statistics_type,&
                                              dbcsr_obj,&
                                              dbcsr_type_int_4
   USE dbcsr_work_operations,           ONLY: dbcsr_create
   USE kinds,                           ONLY: int_4,&
                                              int_8
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
   INTEGER, SAVE :: max_nblocks = 0
   INTEGER, SAVE :: num_multiplications = 0
   REAL, SAVE :: max_memory = 0

   TYPE(dbcsr_memtype_type), SAVE  :: memtype_abpanel_1, memtype_abpanel_2, &
                                      memtype_trsbuffer_1, memtype_trsbuffer_2, &
                                      memtype_mpi_buffer
   TYPE(acc_stream_type), SAVE     :: stream_1, stream_2
   ! ab-panels and streams are shared between all threads

   TYPE dbcsr_mm_multrec_type_p
      TYPE(dbcsr_mm_multrec_type), POINTER :: p => Null()
      ! ensure that array-elements are on different cache lines
      INTEGER(kind=int_4), DIMENSION(64)       :: padding
   END TYPE dbcsr_mm_multrec_type_p

   PUBLIC :: memtype_product_wm
   PUBLIC :: dbcsr_mpi_statistics, num_multiplications
   PUBLIC :: max_nblocks, max_memory

   PUBLIC :: memtype_abpanel_1, memtype_abpanel_2, &
             memtype_trsbuffer_1, memtype_trsbuffer_2, &
             memtype_mpi_buffer
   PUBLIC :: stream_1, stream_2

   PUBLIC :: dbcsr_mm_multrec_type_p
   PUBLIC :: count_mpi_statistics
   PUBLIC :: setup_buffer_matrix
   PUBLIC :: rec_sort_index
   PUBLIC :: acc_transpose_blocks, enumerate_blk_sizes

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
      TYPE(dbcsr_obj), INTENT(INOUT)                     :: matrix
      TYPE(dbcsr_obj), INTENT(IN)                        :: source_matrix
      INTEGER, INTENT(IN), OPTIONAL                      :: index_size, data_size
      TYPE(dbcsr_data_obj), INTENT(IN), OPTIONAL         :: data_buffer
      TYPE(dbcsr_memtype_type), INTENT(IN), OPTIONAL     :: data_memory_type

      CALL dbcsr_init(matrix)
      CALL dbcsr_create(matrix, &
                        template=source_matrix, &
                        name=TRIM("Buffer of "//TRIM(source_matrix%m%name)), &
                        nze=data_size, &
                        data_buffer=data_buffer, &
                        data_memory_type=data_memory_type, &
                        index_memory_type=memtype_mpi_buffer)
      IF (PRESENT(data_size)) THEN
         CALL dbcsr_data_ensure_size( &
            matrix%m%data_area, &
            data_size, nocopy=.TRUE.)
      ENDIF
      IF (PRESENT(index_size)) THEN
         CALL ensure_array_size( &
            matrix%m%index, &
            ub=index_size, nocopy=.TRUE., &
            memory_type=dbcsr_get_index_memory_type(matrix))
      ENDIF
      matrix%m%negate_real = source_matrix%m%negate_real
      matrix%m%negate_imaginary = source_matrix%m%negate_imaginary
      matrix%m%local_indexing = source_matrix%m%local_indexing
      matrix%m%list_indexing = source_matrix%m%list_indexing
      !
      IF (source_matrix%m%has_local_rows) THEN
         matrix%m%local_rows = source_matrix%m%local_rows
         CALL array_hold(matrix%m%local_rows)
         matrix%m%has_local_rows = .TRUE.
      ENDIF
      IF (source_matrix%m%has_global_rows) THEN
         matrix%m%global_rows = source_matrix%m%global_rows
         CALL array_hold(matrix%m%global_rows)
         matrix%m%has_global_rows = .TRUE.
      ENDIF
      IF (source_matrix%m%has_local_cols) THEN
         matrix%m%local_cols = source_matrix%m%local_cols
         CALL array_hold(matrix%m%local_cols)
         matrix%m%has_local_cols = .TRUE.
      ENDIF
      IF (source_matrix%m%has_global_cols) THEN
         matrix%m%global_cols = source_matrix%m%global_cols
         CALL array_hold(matrix%m%global_cols)
         matrix%m%has_global_cols = .TRUE.
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
!> \brief write out a stack for transposing the blocks
!> \param matrix ...
!> \param trs_stackbuf ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE acc_transpose_blocks(matrix, trs_stackbuf)
      TYPE(dbcsr_obj), INTENT(IN)                        :: matrix
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: trs_stackbuf

      CHARACTER(len=*), PARAMETER :: routineN = 'acc_transpose_blocks', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk_p, col, handle, i, m, mi, mi_max, n, &
                                                            nblks, ni, ni_max, offset, row, x
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: counters, filled, offsets
      LOGICAL                                            :: local
      INTEGER, DIMENSION(:), POINTER :: blk_index, col_blk_sizes, col_blk_sizes2enum, &
         enum2col_blk_sizes, enum2row_blk_sizes, local2global_cols, local2global_rows, &
         row_blk_sizes, row_blk_sizes2enum, trs_stack
      INTEGER, DIMENSION(3, matrix%m%nblks)              :: tmp_stack

      CALL timeset(routineN, handle)

      NULLIFY (row_blk_sizes2enum, enum2row_blk_sizes)
      NULLIFY (col_blk_sizes2enum, enum2col_blk_sizes)
      NULLIFY (local2global_rows, local2global_cols, trs_stack)

      IF (.NOT. matrix%m%list_indexing) &
         CPABORT("build_trs_stack: only list_indexing supported.")
      IF (trs_stackbuf%d%data_type /= dbcsr_type_int_4) &
         CPABORT("build_trs_stack: stac_buf has wrong datatype")
      blk_index => matrix%m%coo_l
      row_blk_sizes => array_data(matrix%m%row_blk_size)
      col_blk_sizes => array_data(matrix%m%col_blk_size)
      local = matrix%m%local_indexing
      local2global_rows => array_data(matrix%m%local_rows)
      local2global_cols => array_data(matrix%m%local_cols)
      nblks = matrix%m%nblks

      ! enumerate the blocksizes to keep the following 2D-arrays small.
      CALL enumerate_blk_sizes(row_blk_sizes, row_blk_sizes2enum, enum2row_blk_sizes)
      CALL enumerate_blk_sizes(col_blk_sizes, col_blk_sizes2enum, enum2col_blk_sizes)
      mi_max = SIZE(enum2row_blk_sizes); ni_max = SIZE(enum2col_blk_sizes)
      ALLOCATE (counters(mi_max, ni_max), offsets(mi_max, ni_max), filled(mi_max, ni_max))
      counters(:, :) = 0; offsets(:, :) = 0; filled(:, :) = 0

      ! make sure buffer from previous cannon-tick was uploaded
      CALL acc_event_synchronize(trs_stackbuf%d%acc_ready)

      CALL dbcsr_data_ensure_size(trs_stackbuf, data_size=nblks, nocopy=.TRUE.)
      trs_stack => trs_stackbuf%d%i4

      ! collect block addresses and dimensions in a temporary stack
      ! while doing so, also count number of blocks per block-dimensions
      DO i = 1, nblks
         row = blk_index(3*(i-1)+1)
         col = blk_index(3*(i-1)+2)
         blk_p = blk_index(3*(i-1)+3)
         IF (blk_p == 0) CYCLE
         IF (local) THEN
            row = local2global_rows(row)
            col = local2global_cols(col)
         ENDIF
         m = row_blk_sizes(row)
         n = col_blk_sizes(col)
         mi = row_blk_sizes2enum(m)
         ni = col_blk_sizes2enum(n)
         tmp_stack(1, i) = mi
         tmp_stack(2, i) = ni
         tmp_stack(3, i) = blk_p-1
         counters(mi, ni) = counters(mi, ni)+1
      ENDDO

      ! calculate offsets for first element of each sub-stack
      offset = 0
      DO mi = 1, mi_max
         DO ni = 1, ni_max
            offsets(mi, ni) = offset
            offset = offset+counters(mi, ni)
         ENDDO
      ENDDO

      ! write all sub-stacks into the host-pinned buffer
      DO i = 1, nblks
         mi = tmp_stack(1, i)
         ni = tmp_stack(2, i)
         blk_p = tmp_stack(3, i)
         x = offsets(mi, ni)+filled(mi, ni)+1
         trs_stack(x) = blk_p
         filled(mi, ni) = filled(mi, ni)+1
      ENDDO

      !sanity check
      DO mi = 1, mi_max
         DO ni = 1, ni_max
            IF (filled(mi, ni) /= counters(mi, ni)) &
               CPABORT("acc_transpose_blocks: bug")
         END DO
      END DO

      !transfer all stacks
      CALL dbcsr_data_host2dev(trs_stackbuf)

      ! make sure block-buffer is uploaded befor running the kernels
      CALL acc_stream_wait_event(trs_stackbuf%d%memory_type%acc_stream, matrix%m%data_area%d%acc_ready)

      ! launch kernels
      DO mi = 1, mi_max
         DO ni = 1, ni_max
            IF (counters(mi, ni) > 0) THEN
               m = enum2row_blk_sizes(mi)
               n = enum2col_blk_sizes(ni)
               CALL dbcsr_acc_transpose( &
                  trs_stack=trs_stackbuf%d%acc_devmem, &
                  offset=offsets(mi, ni), &
                  nblks=counters(mi, ni), &
                  datatype=matrix%m%data_type, &
                  buffer=matrix%m%data_area%d%acc_devmem, &
                  m=m, n=n, &
                  stream=trs_stackbuf%d%memory_type%acc_stream)
            END IF
         ENDDO
      ENDDO

      ! make sure block-buffer are not used until transpose kernels finnished
      CALL acc_event_record(trs_stackbuf%d%acc_ready, trs_stackbuf%d%memory_type%acc_stream)
      CALL acc_stream_wait_event(matrix%m%data_area%d%memory_type%acc_stream, trs_stackbuf%d%acc_ready)
      CALL acc_event_record(matrix%m%data_area%d%acc_ready, matrix%m%data_area%d%memory_type%acc_stream)

      DEALLOCATE (row_blk_sizes2enum, enum2row_blk_sizes)
      DEALLOCATE (col_blk_sizes2enum, enum2col_blk_sizes)
      CALL timestop(handle)
   END SUBROUTINE acc_transpose_blocks

! **************************************************************************************************
!> \brief Enumerate all occuring blocksizes
!> \param blk_sizes ...
!> \param enum ...
!> \param rev_enum ...
!>  \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE enumerate_blk_sizes(blk_sizes, enum, rev_enum)
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes, enum, rev_enum

      INTEGER                                            :: i, n

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
   END SUBROUTINE enumerate_blk_sizes

END MODULE dbcsr_mm_common
