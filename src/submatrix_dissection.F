!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \note
!> This module implements a modified version of the submatrix method, proposed in
!>   M. Lass, S. Mohr, H. Wiebeler, T. Kuehne, C. Plessl
!>   A Massively Parallel Algorithm for the Approximate Calculation of Inverse p-th Roots of Large Sparse Matrices
!>   Proc. Platform for Advanced Scientific Computing (PASC) Conference, ACM, 2018
!>
!> The method is extended to minimize the required data transfers and floating-point operations under the assumption that non-zero
!> blocks of the matrix are close to its diagonal.
!>
!> Submatrices can be constructed not for single columns of the matrix but for a set of w consecutive submatrices. The underlying
!> assumption is that columns next to each other have relatively similar occupation patterns, i.e., constructing a submatrix from
!> columns x and x+1 should not lead to a much bigger submatrix than contructing it only from column x.
!>
!> The construction of the submatrices requires communication between all ranks. It is crucial to implement this communication as
!> efficient as possible, i.e., data should only ever be transferred once between two ranks and message sizes need to be
!> sufficiently large to utilize the communication bandwidth. To achieve this, we communicate the required blocks for all
!> submatrices at once and copy them to large buffers before transmitting them via MPI.
!>
!> Note on multi-threading:
!> Submatrices can be constructed and processed in parallel by multiple threads. However, generate_submatrix, get_sm_ids_for_rank
!> and copy_resultcol are the only thread-safe routines in this module. All other routines involve MPI communication or operate on
!> common, non-protected data and are hence not thread-safe.
!>
!> TODO:
!> * generic types (for now only dp supported)
!> * optimization of threaded initialization
!> * sanity checks at the beginning of all methods
!> * optimize execution without MPI, i.e., if __parallel is not defined
!>
!> \author Michael Lass
! **************************************************************************************************

MODULE submatrix_dissection

   USE dbcsr_api, ONLY: dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_type, dbcsr_finalize, &
                        dbcsr_get_block_p, dbcsr_get_info, dbcsr_get_stored_coordinates, &
                        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
                        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_put_block, dbcsr_type
   USE dbcsr_mpiwrap, ONLY: mp_allgather, mp_alltoall, mp_irecv, mp_isend, mp_send, mp_wait
   USE kinds, ONLY: dp
   USE submatrix_methods, ONLY: qsort_two
   USE submatrix_types, ONLY: buffer_type, bufptr_type, intBuffer_type, set_type, setarray_type

#ifdef _OPENMP
   USE omp_lib, ONLY: omp_get_max_threads, omp_get_thread_num
#endif

#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PARAMETER   :: DEBUG = .FALSE.
   CHARACTER, PARAMETER :: endl = NEW_LINE(endl)

   TYPE, PUBLIC :: submatrix_dissection_type
      TYPE(dbcsr_type)                                :: dbcsr_mat
      TYPE(dbcsr_distribution_type)                   :: dist
      LOGICAL                                         :: initialized = .FALSE.
      INTEGER                                         :: numnodes, group, myrank, nblkcols, nblkrows, nblks, local_blocks, &
                                                         cols_per_sm, number_of_submatrices
      INTEGER, DIMENSION(:), POINTER                  :: row_blk_size, col_blk_size
      INTEGER, DIMENSION(:), ALLOCATABLE              :: coo_cols, coo_rows, coo_col_offsets, coo_cols_local, coo_rows_local, &
                                                         coo_col_offsets_local, submatrix_owners, submatrix_sizes
      TYPE(buffer_type), DIMENSION(:), ALLOCATABLE    :: recvbufs, result_sendbufs ! Indexing starts with 0 to match rank ids!
      TYPE(set_type), DIMENSION(:), ALLOCATABLE       :: result_blocks_for_rank, result_blocks_from_rank
      TYPE(bufptr_type), DIMENSION(:), ALLOCATABLE    :: coo_dptr
      TYPE(intBuffer_type), DIMENSION(:), ALLOCATABLE :: result_blocks_for_rank_buf_offsets
   CONTAINS
      PROCEDURE :: init => submatrix_dissection_init
      PROCEDURE :: final => submatrix_dissection_final
      PROCEDURE :: get_sm_ids_for_rank => submatrix_get_sm_ids_for_rank
      PROCEDURE :: generate_submatrix => submatrix_generate_sm
      PROCEDURE :: copy_resultcol => submatrix_cpy_resultcol
      PROCEDURE :: communicate_results => submatrix_communicate_results
      PROCEDURE :: get_first_relevant_col => submatrix_get_first_relevant_col
   END TYPE submatrix_dissection_type

CONTAINS

! **************************************************************************************************
!> \brief determine which columns of the submatrix are relevant for the result matrix
!> \param this - object of class submatrix_dissection_type
!> \param sm_id - id of the submatrix
!> \return first column of interest for result matrix
! **************************************************************************************************
   FUNCTION submatrix_get_first_relevant_col(this, sm_id) RESULT(sm_col_offset)
      CLASS(submatrix_dissection_type), INTENT(IN)     :: this
      INTEGER, INTENT(IN)                              :: sm_id
      INTEGER                                          :: sm_col_offset, i, j
      TYPE(set_type)                                   :: nonzero_rows

      ! TODO: Should we buffer the list of non-zero rows for each submatrix instead of recalculating it each time?
      sm_col_offset = 0
      DO i = (sm_id - 1)*this%cols_per_sm + 1, sm_id*this%cols_per_sm     ! all colums that determine submatrix sm_id
         DO j = this%coo_col_offsets(i), this%coo_col_offsets(i + 1) - 1 ! all blocks that are within this column
            CALL nonzero_rows%insert(this%coo_rows(j))
         ENDDO
      ENDDO

      sm_col_offset = 0
      DO i = 1, nonzero_rows%elements
         IF (nonzero_rows%get(i) .EQ. (sm_id - 1)*this%cols_per_sm + 1) THEN
            ! We just found the nonzero row what corresponds to the first inducing column of our submatrix
            sm_col_offset = i
            IF (DEBUG) PRINT *, "For sm", sm_id, "the relevant result columns start with column", sm_col_offset
            EXIT
         ENDIF
      ENDDO
      IF (sm_col_offset .EQ. 0) THEN
         CPABORT("Could not determine relevant result columns of submatrix")
      ENDIF

      CALL nonzero_rows%reset

   END FUNCTION submatrix_get_first_relevant_col

! **************************************************************************************************
!> \brief initialize submatrix dissection and communicate, needs to be called before constructing any submatrices.
!> \param this - object of class submatrix_dissection_type
!> \param matrix_p - dbcsr input matrix
! **************************************************************************************************
   SUBROUTINE submatrix_dissection_init(this, matrix_p) ! Should be PURE but the iterator routines are not
      CLASS(submatrix_dissection_type), INTENT(INOUT)  :: this
      TYPE(dbcsr_type), INTENT(IN)                     :: matrix_p

      INTEGER                                          :: cur_row, cur_col, cur_blk, i, j, k, l, m, l_limit_left, l_limit_right, &
                                                          bufsize, bufsize_next
      INTEGER, DIMENSION(:), ALLOCATABLE               :: blocks_per_rank, coo_dsplcmnts, num_blockids_send, num_blockids_recv
      TYPE(dbcsr_iterator_type)                        :: iter
      TYPE(set_type)                                   :: nonzero_rows
      REAL(KIND=dp)                                    :: flops_total, flops_per_rank, flops_per_sm, flops_threshold, &
                                                          flops_current, flops_remaining

      ! Indexing for the following arrays starts with 0 to match rank ids
      TYPE(set_type), DIMENSION(:), ALLOCATABLE        :: blocks_from_rank
      TYPE(buffer_type), DIMENSION(:), ALLOCATABLE     :: sendbufs
      TYPE(intBuffer_type), DIMENSION(:), ALLOCATABLE  :: blocks_for_rank

      ! Additional structures for threaded parts
      INTEGER                                          :: numthreads, mytid
      ! Indexing starts at 0 to match thread ids
      TYPE(setarray_type), DIMENSION(:), ALLOCATABLE   :: nonzero_rows_t
      TYPE(setarray_type), DIMENSION(:), ALLOCATABLE   :: result_blocks_from_rank_t, result_blocks_for_rank_t, &
                                                          blocks_from_rank_t

      LOGICAL                                          :: valid
      REAL(KIND=dp), DIMENSION(:), POINTER             :: blockp

      this%dbcsr_mat = matrix_p
      CALL dbcsr_get_info(matrix=this%dbcsr_mat, nblkcols_total=this%nblkcols, nblkrows_total=this%nblkrows, &
                          row_blk_size=this%row_blk_size, col_blk_size=this%col_blk_size, group=this%group, distribution=this%dist)
      CALL dbcsr_distribution_get(dist=this%dist, mynode=this%myrank, numnodes=this%numnodes)

      IF (this%nblkcols .NE. this%nblkrows) THEN
         CPABORT("Number of block rows and cols need to be identical")
      ENDIF

      DO i = 1, this%nblkcols
         IF (this%col_blk_size(i) .NE. this%row_blk_size(i)) THEN
            CPABORT("Block row sizes and col sizes need to be identical")
         ENDIF
      ENDDO

      ! TODO: Do some more sanity checks here, e.g., the matrix must not be stored symmetrically

      ! For the submatrix method, we need global knwoledge about which blocks are actually used. Therefore, we create a COO
      ! representation of the blocks (without their contents) on all ranks.

      ! TODO: Right now, the COO contains all blocks. Also we transfer all blocks. We could skip half of them due to the matrix
      ! being symmetric (however, we need to transpose the blocks). This can increase performance only by a factor of 2 and is
      ! therefore deferred.

      ! Determine number of locally stored blocks
      this%local_blocks = 0
      CALL dbcsr_iterator_start(iter, this%dbcsr_mat, read_only=.TRUE.)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, cur_row, cur_col, cur_blk)
         this%local_blocks = this%local_blocks + 1
      ENDDO
      CALL dbcsr_iterator_stop(iter)

      ALLOCATE (this%coo_cols_local(this%local_blocks), this%coo_rows_local(this%local_blocks), blocks_per_rank(this%numnodes), &
                coo_dsplcmnts(this%numnodes), this%coo_col_offsets_local(this%nblkcols + 1), &
                blocks_for_rank(0:this%numnodes - 1), blocks_from_rank(0:this%numnodes - 1), &
                sendbufs(0:this%numnodes - 1), this%recvbufs(0:this%numnodes - 1), this%result_sendbufs(0:this%numnodes - 1), &
                this%result_blocks_for_rank(0:this%numnodes - 1), this%result_blocks_from_rank(0:this%numnodes - 1), &
                this%result_blocks_for_rank_buf_offsets(0:this%numnodes - 1))

      i = 1
      CALL dbcsr_iterator_start(iter, this%dbcsr_mat, read_only=.TRUE.)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, cur_row, cur_col, cur_blk)
         this%coo_cols_local(i) = cur_col
         this%coo_rows_local(i) = cur_row
         i = i + 1
      ENDDO
      CALL dbcsr_iterator_stop(iter)

      IF (DEBUG) THEN
         PRINT *, this%myrank, "Rows:", this%coo_rows_local
         PRINT *, this%myrank, "Cols:", this%coo_cols_local
      ENDIF

      ! We only know how many blocks we own. What's with the other ranks?
#ifdef __parallel
      CALL mp_allgather(msgout=this%local_blocks, msgin=blocks_per_rank, gid=this%group)
#else
      blocks_per_rank(1) = this%local_blocks
#endif
      IF (DEBUG) PRINT *, this%myrank, "Blocks per rank:", blocks_per_rank
      coo_dsplcmnts(1) = 0
      DO i = 1, this%numnodes - 1
         coo_dsplcmnts(i + 1) = coo_dsplcmnts(i) + blocks_per_rank(i)
      ENDDO
      IF (DEBUG) PRINT *, this%myrank, "Displacements:", coo_dsplcmnts

      ! Get a global view on the matrix
      this%nblks = SUM(blocks_per_rank)
      ALLOCATE (this%coo_cols(this%nblks), this%coo_rows(this%nblks), this%coo_col_offsets(this%nblkcols + 1), &
                this%coo_dptr(this%nblks))
#ifdef __parallel
      CALL mp_allgather(msgout=this%coo_rows_local, msgin=this%coo_rows, gid=this%group, rcount=blocks_per_rank, &
                        rdispl=coo_dsplcmnts)
      CALL mp_allgather(msgout=this%coo_cols_local, msgin=this%coo_cols, gid=this%group, rcount=blocks_per_rank, &
                        rdispl=coo_dsplcmnts)
#else
      this%coo_rows(:) = this%coo_rows_local
      this%coo_cols(:) = this%coo_cols_local
#endif

      DEALLOCATE (blocks_per_rank, coo_dsplcmnts)

      ! Sort COO arrays according to their columns
      CALL qsort_two(this%coo_cols_local, this%coo_rows_local, 1, this%local_blocks)
      CALL qsort_two(this%coo_cols, this%coo_rows, 1, this%nblks)

      IF (DEBUG) THEN
         PRINT *, endl, this%myrank, "Global COO Cols:", this%coo_cols, endl
         PRINT *, endl, this%myrank, "Global COO Rows:", this%coo_rows, endl
         PRINT *, endl, this%myrank, "Local COO Cols:", this%coo_cols_local, endl
         PRINT *, endl, this%myrank, "Local COO Rows:", this%coo_rows_local, endl
      ENDIF

      ! Get COO array offsets for columns to accelerate lookups
      this%coo_col_offsets(this%nblkcols + 1) = this%nblks + 1
      j = 1
      DO i = 1, this%nblkcols
         DO WHILE ((j .LE. this%nblks))
            IF (this%coo_cols(j) .GE. i) EXIT
            j = j + 1
         ENDDO
         this%coo_col_offsets(i) = j
      ENDDO

      ! Same for local COO
      this%coo_col_offsets_local(this%nblkcols + 1) = this%local_blocks + 1
      j = 1
      DO i = 1, this%nblkcols
         DO WHILE ((j .LE. this%local_blocks))
            IF (this%coo_cols_local(j) .GE. i) EXIT
            j = j + 1
         ENDDO
         this%coo_col_offsets_local(i) = j
      ENDDO

      IF (DEBUG) THEN
         PRINT *, endl, this%myrank, "Global COO Offsets:", endl, this%coo_col_offsets, endl
         PRINT *, endl, this%myrank, "Local COO Offsets:", endl, this%coo_col_offsets_local, endl
      ENDIF

      ! We could combine multiple columns to generate a single submatrix. For now, we have not found a practical use-case for this
      ! so we only look at single columns for now.
      this%cols_per_sm = 1

      ! Determine sizes of all submatrices. This is required in order to assess the computational effort that is required to process
      ! the submatrices.
      this%number_of_submatrices = this%nblkcols/this%cols_per_sm
      ALLOCATE (this%submatrix_sizes(this%number_of_submatrices))
      this%submatrix_sizes = 0
      flops_total = 0.0D0
      DO i = 1, this%number_of_submatrices
         CALL nonzero_rows%reset
         DO j = (i - 1)*this%cols_per_sm + 1, i*this%cols_per_sm ! all colums that determine submatrix i
            DO k = this%coo_col_offsets(j), this%coo_col_offsets(j + 1) - 1 ! all blocks that are within this column
               CALL nonzero_rows%insert(this%coo_rows(k))
            ENDDO
         ENDDO
         DO j = 1, nonzero_rows%elements
            this%submatrix_sizes(i) = this%submatrix_sizes(i) + this%col_blk_size(nonzero_rows%get(j))
         ENDDO
         flops_total = flops_total + 2.0D0*this%submatrix_sizes(i)*this%submatrix_sizes(i)*this%submatrix_sizes(i)
      ENDDO
      IF (DEBUG) THEN
         PRINT *, "Submatrix sizes:", this%submatrix_sizes
         PRINT *, "Total flops:", flops_total
      ENDIF

      ! Create mapping from submatrices to ranks. Since submatrices can be of different sizes, we need to perform some load
      ! balancing here. For that we assume that arithmetic operations performed on the submatrices scale cubically.
      IF ((this%myrank .EQ. 0) .AND. DEBUG) &
         PRINT *, "Creating", this%number_of_submatrices, "submatrices from", this%cols_per_sm, "cols each"
      ALLOCATE (this%submatrix_owners(this%number_of_submatrices))
      flops_per_sm = flops_total/this%number_of_submatrices
      flops_per_rank = flops_total/this%numnodes
      flops_current = 0.0D0
      j = 0
      DO i = 1, this%number_of_submatrices
         this%submatrix_owners(i) = j
         flops_current = flops_current + 2.0D0*this%submatrix_sizes(i)*this%submatrix_sizes(i)*this%submatrix_sizes(i)
         flops_remaining = flops_total - flops_current
         flops_threshold = (this%numnodes - j - 1)*flops_per_rank
         IF ((j .LT. (this%numnodes - 1)) &
             .AND. ((flops_remaining .LE. flops_threshold &
                     .OR. (this%number_of_submatrices - i) .LT. (this%numnodes - j)))) THEN
            IF (DEBUG) PRINT *, "Flops for rank", j, ":", flops_current
            j = j + 1
            flops_total = flops_total - flops_current
            flops_current = 0.0D0
         ENDIF
      ENDDO
      IF (DEBUG) PRINT *, "Submatrix owners:", this%submatrix_owners

      ! Prepare data structures for multithreaded loop
#ifdef _OPENMP
      numthreads = omp_get_max_threads()
#else
      numthreads = 1
#endif
      ALLOCATE (result_blocks_from_rank_t(0:numthreads - 1), &
                result_blocks_for_rank_t(0:numthreads - 1), &
                blocks_from_rank_t(0:numthreads - 1), &
                nonzero_rows_t(0:numthreads - 1))

      ! Figure out which blocks we need to receive. Blocks are identified here as indices into our COO representation.
      ! FIXME: This currently shows bad parallel efficiency.
      !$OMP PARALLEL NUM_THREADS(numthreads) DEFAULT(NONE) &
      !$OMP          PRIVATE(i,j,k,l,m,l_limit_left,l_limit_right,cur_col,cur_row,mytid) &
      !$OMP          SHARED(result_blocks_from_rank_t,result_blocks_for_rank_t,blocks_from_rank_t,this,numthreads,nonzero_rows_t)
#ifdef _OPENMP
      mytid = omp_get_thread_num()
#else
      mytid = 0
#endif
      IF (DEBUG) PRINT *, "Hello from thread", mytid, "/", numthreads

      ALLOCATE (nonzero_rows_t(mytid)%sets(1), &
                result_blocks_from_rank_t(mytid)%sets(0:this%numnodes - 1), &
                result_blocks_for_rank_t(mytid)%sets(0:this%numnodes - 1), &
                blocks_from_rank_t(mytid)%sets(0:this%numnodes - 1))

      !$OMP DO schedule(guided)
      DO i = 1, this%number_of_submatrices
         CALL nonzero_rows_t(mytid)%sets(1)%reset
         DO j = (i - 1)*this%cols_per_sm + 1, i*this%cols_per_sm ! all colums that determine submatrix i
            DO k = this%coo_col_offsets(j), this%coo_col_offsets(j + 1) - 1 ! all blocks that are within this column
               ! This block will be required to assemble the final block matrix as it is within an inducing column for submatrix i.
               ! Figure out who stores it and insert it into the result_blocks_* sets.
               CALL dbcsr_get_stored_coordinates(matrix=this%dbcsr_mat, row=this%coo_rows(k), column=j, processor=m)
               IF (m .EQ. this%myrank) THEN
                  CALL result_blocks_from_rank_t(mytid)%sets(this%submatrix_owners(i))%insert(k)
               ENDIF
               IF (this%submatrix_owners(i) .EQ. this%myrank) THEN
                  CALL nonzero_rows_t(mytid)%sets(1)%insert(this%coo_rows(k))
                  CALL result_blocks_for_rank_t(mytid)%sets(m)%insert(k)
               ENDIF
            ENDDO
         ENDDO

         IF (this%submatrix_owners(i) .NE. this%myrank) CYCLE

         IF (DEBUG) PRINT *, "nonzero rows: ", nonzero_rows_t(mytid)%sets(1)%getall()

         ! In the following, we determine all blocks required to build the submatrix. We interpret nonzero_rows_t(mytid)(j) as
         ! column and nonzero_rows_t(mytid)(k) as row.
         DO j = 1, nonzero_rows_t(mytid)%sets(1)%elements
            cur_col = nonzero_rows_t(mytid)%sets(1)%get(j)
            l_limit_left = this%coo_col_offsets(cur_col)
            l_limit_right = this%coo_col_offsets(cur_col + 1) - 1
            ! IF (DEBUG) print *, this%coo_rows(l_limit_left:l_limit_right)
            DO k = 1, nonzero_rows_t(mytid)%sets(1)%elements
               cur_row = nonzero_rows_t(mytid)%sets(1)%get(k)
               ! IF (DEBUG) print *, "Looking for row", cur_row
               l = l_limit_left
               DO WHILE (l .LE. l_limit_right)
                  IF (this%coo_rows(l) .GE. cur_row) EXIT
                  l = l + 1
               ENDDO
               l_limit_left = l
               IF (l .LE. l_limit_right) THEN
                  IF (this%coo_rows(l) .EQ. cur_row) THEN
                     ! We found a valid block. Figure out what to do with it.
                     ! IF (DEBUG) print *, "Found row", this%coo_rows(l)
                     IF (DEBUG) PRINT *, "(", cur_row, ",", cur_col, ")", l
                     CALL dbcsr_get_stored_coordinates(matrix=this%dbcsr_mat, row=this%coo_rows(l), &
                                                       column=this%coo_cols(l), processor=m)
                     CALL blocks_from_rank_t(mytid)%sets(m)%insert(l)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL

      ! Merge partial results from threads
      DO i = 0, this%numnodes - 1
         DO j = 0, numthreads - 1
            DO k = 1, result_blocks_from_rank_t(j)%sets(i)%elements
               CALL this%result_blocks_from_rank(i)%insert(result_blocks_from_rank_t(j)%sets(i)%get(k))
            ENDDO
            CALL result_blocks_from_rank_t(j)%sets(i)%reset
            DO k = 1, result_blocks_for_rank_t(j)%sets(i)%elements
               CALL this%result_blocks_for_rank(i)%insert(result_blocks_for_rank_t(j)%sets(i)%get(k))
            ENDDO
            CALL result_blocks_for_rank_t(j)%sets(i)%reset
            DO k = 1, blocks_from_rank_t(j)%sets(i)%elements
               CALL blocks_from_rank(i)%insert(blocks_from_rank_t(j)%sets(i)%get(k))
            ENDDO
            CALL blocks_from_rank_t(j)%sets(i)%reset
         ENDDO
      ENDDO
      DO i = 0, numthreads - 1
         CALL nonzero_rows_t(i)%sets(1)%reset
         DEALLOCATE (result_blocks_from_rank_t(i)%sets, result_blocks_for_rank_t(i)%sets, blocks_from_rank_t(i)%sets, &
                     nonzero_rows_t(i)%sets)
      ENDDO
      DEALLOCATE (result_blocks_from_rank_t, result_blocks_for_rank_t, blocks_from_rank_t, nonzero_rows_t)

      ! Make other ranks aware of our needs
      ALLOCATE (num_blockids_send(0:this%numnodes - 1), num_blockids_recv(0:this%numnodes - 1))
      DO i = 0, this%numnodes - 1
         num_blockids_send(i) = blocks_from_rank(i)%elements
      ENDDO
      IF (DEBUG) PRINT *, this%myrank, "Number of block ids to send to each rank:", num_blockids_send
#ifdef __parallel
      CALL mp_alltoall(num_blockids_send, num_blockids_recv, 1, this%group)
#else
      num_blockids_recv(:) = num_blockids_send
#endif
      IF (DEBUG) PRINT *, this%myrank, "Number of block ids to receive from each rank:", num_blockids_recv
      DO i = 0, this%numnodes - 1
         CALL blocks_for_rank(i)%alloc(num_blockids_recv(i))
      ENDDO
      DEALLOCATE (num_blockids_send, num_blockids_recv)
#ifdef __parallel
      DO i = 1, this%numnodes
         k = MODULO(this%myrank - i, this%numnodes) ! rank to receive from
         IF (DEBUG) PRINT *, this%myrank, "Expecting message from", k
         CALL mp_irecv(msgout=blocks_for_rank(k)%data, source=k, comm=this%group, request=blocks_for_rank(k)%mpi_request)
      ENDDO
      DO i = 1, this%numnodes ! TODO: use isend and merge with previous loop
         k = MODULO(this%myrank + i, this%numnodes) ! rank to send to
         IF (DEBUG) PRINT *, this%myrank, "Sending message to", k, ":", blocks_from_rank(k)%getall()
         CALL mp_send(blocks_from_rank(k)%getall(), k, 0, this%group)
      ENDDO
      DO i = 0, this%numnodes - 1
         CALL mp_wait(blocks_for_rank(i)%mpi_request)
         IF (DEBUG) PRINT *, this%myrank, "Received message from", i, ":", blocks_for_rank(i)%data
      ENDDO
#else
      blocks_for_rank(0)%data = blocks_from_rank(0)%getall()
#endif

      ! Free memory allocated in nonzero_rows
      CALL nonzero_rows%reset

      ! Sort list of blocks to be sent and received
      DO i = 0, this%numnodes - 1
         IF (DEBUG) THEN
            PRINT *, "Rank", this%myrank, "has blocks", blocks_for_rank(i)%data, "for rank", i
            PRINT *, "Rank", this%myrank, "wants blocks", blocks_from_rank(i)%getall(), "from rank", i
         ENDIF
      ENDDO

      ! Workaround to make get calls on this%result_blocks_for_rank(...) threadsafe in other routines
      DO m = 0, this%numnodes - 1
         IF ((.NOT. this%result_blocks_for_rank(m)%sorted_up_to_date) .AND. (this%result_blocks_for_rank(m)%elements .GT. 0)) THEN
            CALL this%result_blocks_for_rank(m)%update_sorted
         ENDIF
      ENDDO

      ! Create and fill send buffers
      DO i = 0, this%numnodes - 1
         bufsize = 0
         DO j = 1, blocks_for_rank(i)%size
            k = blocks_for_rank(i)%data(j)
            bufsize = bufsize + this%col_blk_size(this%coo_cols(k))*this%col_blk_size(this%coo_rows(k))
         ENDDO
         CALL sendbufs(i)%alloc(bufsize)

         bufsize = 0
         CALL this%result_blocks_for_rank_buf_offsets(i)%alloc(this%result_blocks_for_rank(i)%elements)
         DO j = 1, this%result_blocks_for_rank(i)%elements
            k = this%result_blocks_for_rank(i)%get(j)
            this%result_blocks_for_rank_buf_offsets(i)%data(j) = bufsize
            bufsize = bufsize + this%col_blk_size(this%coo_cols(k))*this%col_blk_size(this%coo_rows(k))
         ENDDO
         CALL this%result_sendbufs(i)%alloc(bufsize)

         bufsize = 0
         DO j = 1, blocks_for_rank(i)%size
            k = blocks_for_rank(i)%data(j)
            CALL dbcsr_get_block_p(this%dbcsr_mat, row=this%coo_rows(k), col=this%coo_cols(k), block=blockp, found=valid)
            IF (.NOT. valid) THEN
               CPABORT("Block included in our COO and placed on our rank could not be fetched!")
            ENDIF
            bufsize_next = bufsize + SIZE(blockp)
            sendbufs(i)%data(bufsize + 1:bufsize_next) = blockp
            bufsize = bufsize_next
         ENDDO
      ENDDO

      ! Create receive buffers and mapping from blocks to memory locations
      DO i = 0, this%numnodes - 1
         bufsize = 0
         DO j = 1, blocks_from_rank(i)%elements
            k = blocks_from_rank(i)%get(j)
            bufsize = bufsize + this%col_blk_size(this%coo_cols(k))*this%col_blk_size(this%coo_rows(k))
         ENDDO
         CALL this%recvbufs(i)%alloc(bufsize)
         bufsize = 0
         DO j = 1, blocks_from_rank(i)%elements
            k = blocks_from_rank(i)%get(j)
            bufsize_next = bufsize + this%col_blk_size(this%coo_cols(k))*this%col_blk_size(this%coo_rows(k))
            this%coo_dptr(k)%target => this%recvbufs(i)%data(bufsize + 1:bufsize_next)
            bufsize = bufsize_next
         ENDDO
      ENDDO

      DO i = 0, this%numnodes - 1
         CALL blocks_for_rank(i)%dealloc
         CALL blocks_from_rank(i)%reset
      ENDDO
      DEALLOCATE (blocks_for_rank, blocks_from_rank)

#ifdef __parallel
      ! Communicate. We attempt to balance communication load in the network here by starting our sends with our right neighbor
      ! and first trying to receive from our left neighbor.
      DO i = 1, this%numnodes
         k = MODULO(this%myrank - i, this%numnodes) ! rank to receive from
         CALL mp_irecv(msgout=this%recvbufs(k)%data, source=k, comm=this%group, request=this%recvbufs(k)%mpi_request)
         k = MODULO(this%myrank + i, this%numnodes) ! rank to send to
         CALL mp_isend(msgin=sendbufs(k)%data, dest=k, comm=this%group, request=sendbufs(k)%mpi_request)
      ENDDO
      DO i = 0, this%numnodes - 1
         CALL mp_wait(sendbufs(i)%mpi_request)
         CALL mp_wait(this%recvbufs(i)%mpi_request)
         IF (DEBUG) PRINT *, this%myrank, "received from", i, ":", endl, this%recvbufs(i)%data
      ENDDO
#else
      this%recvbufs(0)%data = sendbufs(0)%data
#endif

      DO i = 0, this%numnodes - 1
         CALL sendbufs(i)%dealloc
      ENDDO
      DEALLOCATE (sendbufs)

      this%initialized = .TRUE.
   END SUBROUTINE submatrix_dissection_init

! **************************************************************************************************
!> \brief free all associated memory, afterwards submatrix_dissection_init needs to be called again
!> \param this - object of class submatrix_dissection_type
! **************************************************************************************************
   PURE SUBROUTINE submatrix_dissection_final(this)
      CLASS(submatrix_dissection_type), INTENT(INOUT) :: this
      INTEGER                                         :: i

      this%initialized = .FALSE.

      IF (ALLOCATED(this%submatrix_sizes)) DEALLOCATE (this%submatrix_sizes)
      IF (ALLOCATED(this%coo_cols_local)) DEALLOCATE (this%coo_cols_local)
      IF (ALLOCATED(this%coo_rows_local)) DEALLOCATE (this%coo_rows_local)
      IF (ALLOCATED(this%coo_col_offsets_local)) DEALLOCATE (this%coo_col_offsets_local)
      IF (ALLOCATED(this%result_blocks_for_rank_buf_offsets)) THEN
         DO i = 0, this%numnodes - 1
            CALL this%result_blocks_for_rank_buf_offsets(i)%dealloc
         ENDDO
         DEALLOCATE (this%result_blocks_for_rank_buf_offsets)
      ENDIF
      IF (ALLOCATED(this%recvbufs)) THEN
         DO i = 0, this%numnodes - 1
            CALL this%recvbufs(i)%dealloc
         ENDDO
         DEALLOCATE (this%recvbufs)
      ENDIF
      IF (ALLOCATED(this%result_sendbufs)) THEN
         DO i = 0, this%numnodes - 1
            CALL this%result_sendbufs(i)%dealloc
         ENDDO
         DEALLOCATE (this%result_sendbufs)
      ENDIF
      IF (ALLOCATED(this%result_blocks_for_rank)) THEN
         DO i = 0, this%numnodes - 1
            CALL this%result_blocks_for_rank(i)%reset
         ENDDO
         DEALLOCATE (this%result_blocks_for_rank)
      ENDIF
      IF (ALLOCATED(this%result_blocks_from_rank)) THEN
         DO i = 0, this%numnodes - 1
            CALL this%result_blocks_from_rank(i)%reset
         ENDDO
         DEALLOCATE (this%result_blocks_from_rank)
      ENDIF
      IF (ALLOCATED(this%coo_cols)) DEALLOCATE (this%coo_cols)
      IF (ALLOCATED(this%coo_rows)) DEALLOCATE (this%coo_rows)
      IF (ALLOCATED(this%coo_col_offsets)) DEALLOCATE (this%coo_col_offsets)
      IF (ALLOCATED(this%coo_dptr)) DEALLOCATE (this%coo_dptr)
      IF (ALLOCATED(this%submatrix_owners)) DEALLOCATE (this%submatrix_owners)

   END SUBROUTINE submatrix_dissection_final

! **************************************************************************************************
!> \brief generate a specific submatrix
!> \param this - object of class submatrix_dissection_type
!> \param sm_id - id of the submatrix to generate
!> \param sm - generated submatrix
! **************************************************************************************************
   SUBROUTINE submatrix_generate_sm(this, sm_id, sm)
      CLASS(submatrix_dissection_type), INTENT(IN)             :: this
      INTEGER, INTENT(IN)                                      :: sm_id
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT) :: sm

      INTEGER                                                  :: sm_dim, i, j, k, offset_x1, offset_x2, offset_y1, &
                                                                  offset_y2, k_limit_left, k_limit_right
      TYPE(set_type)                                           :: nonzero_rows

      IF (.NOT. this%initialized) THEN
         CPABORT("Submatrix dissection not initialized")
      ENDIF

      IF (this%myrank .NE. this%submatrix_owners(sm_id)) THEN
         CPABORT("This rank is not supposed to construct this submatrix")
      ENDIF

      ! TODO: Should we buffer the list of non-zero rows for each submatrix instead of recalculating it each time?
      CALL nonzero_rows%reset
      DO i = (sm_id - 1)*this%cols_per_sm + 1, sm_id*this%cols_per_sm     ! all colums that determine submatrix sm_id
         DO j = this%coo_col_offsets(i), this%coo_col_offsets(i + 1) - 1 ! all blocks that are within this column
            CALL nonzero_rows%insert(this%coo_rows(j))
         ENDDO
      ENDDO
      sm_dim = 0
      DO i = 1, nonzero_rows%elements
         sm_dim = sm_dim + this%col_blk_size(nonzero_rows%get(i))
      ENDDO

      ALLOCATE (sm(sm_dim, sm_dim))
      sm = 0

      offset_x1 = 0
      DO j = 1, nonzero_rows%elements
         offset_x2 = offset_x1 + this%col_blk_size(nonzero_rows%get(j))
         offset_y1 = 0
         k_limit_left = this%coo_col_offsets(nonzero_rows%get(j))
         k_limit_right = this%coo_col_offsets(nonzero_rows%get(j) + 1) - 1
         DO i = 1, nonzero_rows%elements
            offset_y2 = offset_y1 + this%col_blk_size(nonzero_rows%get(i))
            ! Copy block nonzero_rows(i),nonzero_rows(j) to sm(i,j) (if the block actually exists)
            k = k_limit_left
            DO WHILE (k .LE. k_limit_right)
               IF (this%coo_rows(k) .GE. nonzero_rows%get(i)) EXIT
               k = k + 1
            ENDDO
            k_limit_left = k
            IF (this%coo_rows(k) .EQ. nonzero_rows%get(i)) THEN ! it does exist and k is our block id
               IF (DEBUG) PRINT *, "uses block", k, ": (", nonzero_rows%get(i), ",", nonzero_rows%get(j), ")"
               IF (DEBUG) PRINT *, "block dimensions:", &
                  this%col_blk_size(nonzero_rows%get(i)), "x", this%col_blk_size(nonzero_rows%get(j))
               IF (DEBUG) PRINT *, this%coo_dptr(k)%target
               IF (DEBUG) PRINT *, "copying to sm(", offset_y1 + 1, ":", offset_y2, ",", offset_x1 + 1, ":", offset_x2, ")"
               sm(offset_y1 + 1:offset_y2, offset_x1 + 1:offset_x2) = RESHAPE(this%coo_dptr(k)%target, &
                                                                              (/offset_y2 - offset_y1, offset_x2 - offset_x1/))
            ENDIF
            offset_y1 = offset_y2
         ENDDO
         offset_x1 = offset_x2
      ENDDO

      ! Free memory allocated in nonzero_rows
      CALL nonzero_rows%reset

   END SUBROUTINE submatrix_generate_sm

! **************************************************************************************************
!> \brief determine submatrix ids that are handled by a specific rank
!> \param this - object of class submatrix_dissection_type
!> \param rank - rank id of interest
!> \param sm_ids - list of submatrix ids handled by that rank
! **************************************************************************************************
   SUBROUTINE submatrix_get_sm_ids_for_rank(this, rank, sm_ids)
      CLASS(submatrix_dissection_type), INTENT(IN)    :: this
      INTEGER, INTENT(IN)                             :: rank
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: sm_ids
      INTEGER                                         :: count, i

      IF (.NOT. this%initialized) THEN
         CPABORT("Submatrix dissection not initialized")
      ENDIF

      count = 0
      DO i = 1, this%number_of_submatrices
         IF (this%submatrix_owners(i) .EQ. rank) count = count + 1
      ENDDO

      ALLOCATE (sm_ids(count))

      count = 0
      DO i = 1, this%number_of_submatrices
         IF (this%submatrix_owners(i) .EQ. rank) THEN
            count = count + 1
            sm_ids(count) = i
         ENDIF
      ENDDO

   END SUBROUTINE submatrix_get_sm_ids_for_rank

! **************************************************************************************************
!> \brief copy result columns from a submatrix into result buffer
!> \param this - object of class submatrix_dissection_type
!> \param sm_id - id of the submatrix
!> \param sm - result-submatrix
! **************************************************************************************************
   SUBROUTINE submatrix_cpy_resultcol(this, sm_id, sm)
      CLASS(submatrix_dissection_type), INTENT(INOUT)         :: this
      INTEGER, INTENT(in)                                     :: sm_id
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE, INTENT(IN) :: sm

      TYPE(set_type)                                          :: nonzero_rows
      INTEGER                                                 :: i, j, k, m, sm_col_offset, offset_x1, offset_x2, offset_y1, &
                                                                 offset_y2, bufsize, bufsize_next, k_limit_left, k_limit_right
      INTEGER, DIMENSION(:), ALLOCATABLE                      :: buf_offset_idxs

      IF (.NOT. this%initialized) THEN
         CPABORT("Submatrix dissection is not initizalized")
      ENDIF

      IF (this%myrank .NE. this%submatrix_owners(sm_id)) THEN
         CPABORT("This rank is not supposed to operate on this submatrix")
      ENDIF

      ALLOCATE (buf_offset_idxs(0:this%numnodes - 1))
      buf_offset_idxs = 1

      ! TODO: Use get_first_relevant_col()

      ! TODO: Should we buffer the list of non-zero rows for each submatrix instead of recalculating it each time?
      sm_col_offset = 0
      DO i = (sm_id - 1)*this%cols_per_sm + 1, sm_id*this%cols_per_sm     ! all colums that determine submatrix sm_id
         DO j = this%coo_col_offsets(i), this%coo_col_offsets(i + 1) - 1 ! all blocks that are within this column
            CALL nonzero_rows%insert(this%coo_rows(j))
         ENDDO
      ENDDO

      sm_col_offset = 0
      DO i = 1, nonzero_rows%elements
         IF (nonzero_rows%get(i) .EQ. (sm_id - 1)*this%cols_per_sm + 1) THEN
            ! We just found the nonzero row that corresponds to the first inducing column of our submatrix
            sm_col_offset = i
            IF (DEBUG) PRINT *, "For sm", sm_id, "the relevant result columns start with column", sm_col_offset
            EXIT
         ENDIF
      ENDDO
      IF (sm_col_offset .EQ. 0) THEN
         CPABORT("Could not determine relevant result columns of submatrix")
      ENDIF

      offset_x1 = 0
      DO j = 1, nonzero_rows%elements
         offset_x2 = offset_x1 + this%col_blk_size(nonzero_rows%get(j))
         ! We only want to copy the blocks from the result columns
         IF ((j .GE. sm_col_offset) .AND. (j .LT. sm_col_offset + this%cols_per_sm)) THEN
            offset_y1 = 0
            k_limit_left = this%coo_col_offsets(nonzero_rows%get(j))
            k_limit_right = this%coo_col_offsets(nonzero_rows%get(j) + 1) - 1
            DO i = 1, nonzero_rows%elements
               offset_y2 = offset_y1 + this%col_blk_size(nonzero_rows%get(i))
               ! Check if sm(i,j), i.e., (nonzero_rows(i),nonzero_rows(j)) exists in the original matrix and if so, copy it into the
               ! correct send buffer.
               k = k_limit_left
               DO WHILE (k .LE. k_limit_right)
                  IF (this%coo_rows(k) .GE. nonzero_rows%get(i)) EXIT
                  k = k + 1
               ENDDO
               k_limit_left = k
               IF (this%coo_rows(k) .EQ. nonzero_rows%get(i)) THEN ! it does exist and k is our block id
                  CALL dbcsr_get_stored_coordinates(matrix=this%dbcsr_mat, row=this%coo_rows(k), column=this%coo_cols(k), &
                                                    processor=m)
                  DO WHILE (buf_offset_idxs(m) .LE. this%result_blocks_for_rank(m)%elements)
                     IF (this%result_blocks_for_rank(m)%get(buf_offset_idxs(m)) .GE. k) EXIT
                     buf_offset_idxs(m) = buf_offset_idxs(m) + 1
                  ENDDO
                  IF (this%result_blocks_for_rank(m)%get(buf_offset_idxs(m)) .NE. k) THEN
                     CPABORT("Could not determine buffer offset for block")
                  ENDIF
                  bufsize = this%result_blocks_for_rank_buf_offsets(m)%data(buf_offset_idxs(m))
                  bufsize_next = bufsize + this%col_blk_size(this%coo_cols(k))*this%col_blk_size(this%coo_rows(k))
                  this%result_sendbufs(m)%data(bufsize + 1:bufsize_next) = RESHAPE( &
                                                                           sm(offset_y1 + 1:offset_y2, offset_x1 + 1:offset_x2), &
                                                                           (/bufsize_next - bufsize/))
                  IF (DEBUG) PRINT *, "Copying block", k, "(", nonzero_rows%get(i), ", ", nonzero_rows%get(j), ")", &
                     "= sm(", offset_y1 + 1, ":", offset_y2, ",", offset_x1 + 1, ":", offset_x2, ") of sm", sm_id, &
                     "into result buffer for rank", m, &
                     "at pos", "(", bufsize + 1, ":", bufsize_next, ")"
               ENDIF
               offset_y1 = offset_y2
            ENDDO
         ENDIF
         offset_x1 = offset_x2
      ENDDO

      DEALLOCATE (buf_offset_idxs)

      ! Free memory in set
      CALL nonzero_rows%reset

   END SUBROUTINE submatrix_cpy_resultcol

! **************************************************************************************************
!> \brief finalize results back into a dbcsr matrix
!> \param this - object of class submatrix_dissection_type
!> \param resultmat - result dbcsr matrix
! **************************************************************************************************
   SUBROUTINE submatrix_communicate_results(this, resultmat)
      CLASS(submatrix_dissection_type), INTENT(INOUT)                 :: this
      TYPE(dbcsr_type), INTENT(INOUT)                                 :: resultmat

      INTEGER                                                         :: i, j, k, cur_row, cur_col, cur_sm, cur_buf, last_buf, &
                                                                         bufsize, bufsize_next
      TYPE(buffer_type), DIMENSION(:), ALLOCATABLE                    :: result_recvbufs

      ALLOCATE (result_recvbufs(0:this%numnodes - 1))
      DO i = 0, this%numnodes - 1
         bufsize = 0
         DO j = 1, this%result_blocks_from_rank(i)%elements
            k = this%result_blocks_from_rank(i)%get(j)
            bufsize = bufsize + this%col_blk_size(this%coo_cols(k))*this%col_blk_size(this%coo_rows(k))
         ENDDO
         CALL result_recvbufs(i)%alloc(bufsize)
      ENDDO

      ! Communicate. We attempt to balance communication load in the network here by starting our sends with our right neighbor
      ! and first trying to receive from our left neighbor.
#ifdef __parallel
      DO i = 1, this%numnodes
         k = MODULO(this%myrank - i, this%numnodes) ! rank to receive from
         CALL mp_irecv(msgout=result_recvbufs(k)%data, source=k, comm=this%group, request=result_recvbufs(k)%mpi_request)
         k = MODULO(this%myrank + i, this%numnodes) ! rank to send to
         CALL mp_isend(msgin=this%result_sendbufs(k)%data, dest=k, comm=this%group, request=this%result_sendbufs(k)%mpi_request)
      ENDDO
      DO i = 0, this%numnodes - 1
         CALL mp_wait(this%result_sendbufs(i)%mpi_request)
         CALL mp_wait(result_recvbufs(i)%mpi_request)
         IF (DEBUG) PRINT *, this%myrank, "received from", i, ":", endl, result_recvbufs(i)%data
      ENDDO
#else
      result_recvbufs(0)%data = this%result_sendbufs(0)%data
#endif

      last_buf = -1
      bufsize = 0
      DO i = 1, this%local_blocks
         cur_col = this%coo_cols_local(i)
         cur_row = this%coo_rows_local(i)
         cur_sm = (cur_col - 1)/this%cols_per_sm + 1
         cur_buf = this%submatrix_owners(cur_sm)
         IF (cur_buf .GT. last_buf) bufsize = 0
         bufsize_next = bufsize + this%col_blk_size(cur_row)*this%col_blk_size(cur_col)
         CALL dbcsr_put_block(matrix=resultmat, row=cur_row, col=cur_col, &
                              block=result_recvbufs(cur_buf)%data(bufsize + 1:bufsize_next))
         bufsize = bufsize_next
         last_buf = cur_buf
      ENDDO

      DO i = 0, this%numnodes - 1
         CALL result_recvbufs(i)%dealloc
      ENDDO
      DEALLOCATE (result_recvbufs)

      CALL dbcsr_finalize(resultmat)

   END SUBROUTINE submatrix_communicate_results

END MODULE submatrix_dissection
