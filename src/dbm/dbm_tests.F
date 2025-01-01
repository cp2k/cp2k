!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: BSD-3-Clause                                                          !
!--------------------------------------------------------------------------------------------------!

MODULE dbm_tests
   USE OMP_LIB,                         ONLY: omp_get_wtime
   USE dbm_api,                         ONLY: &
        dbm_checksum, dbm_copy, dbm_create, dbm_create_from_template, dbm_distribution_new, &
        dbm_distribution_obj, dbm_distribution_release, dbm_get_col_block_sizes, &
        dbm_get_row_block_sizes, dbm_get_stored_coordinates, dbm_multiply, dbm_put_block, &
        dbm_release, dbm_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE machine,                         ONLY: m_flush
   USE message_passing,                 ONLY: mp_cart_type,&
                                              mp_comm_type,&
                                              mp_dims_create
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: dbm_run_tests, generate_larnv_seed

   INTEGER, PRIVATE, SAVE                                      :: randmat_counter = 0

CONTAINS

! **************************************************************************************************
!> \brief Tests the DBM library.
!> \param mp_group         MPI communicator
!> \param io_unit          Unit to write to, if not negative
!> \param matrix_sizes     Size of matrices to test
!> \param trs              Transposes of the two matrices
!> \param bs_m             Block sizes of the 1. dimension
!> \param bs_n             Block sizes of the 2. dimension
!> \param bs_k             Block sizes of the 3. dimension
!> \param sparsities       Sparsities of matrices to create
!> \param alpha            Alpha value to use in multiply
!> \param beta             Beta value to use in multiply
!> \param n_loops          Number of repetition for each multiplication
!> \param eps              Epsilon value for filtering
!> \param retain_sparsity  Retain the result matrix's sparsity
!> \param always_checksum  Checksum after each multiplication
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_run_tests(mp_group, io_unit, matrix_sizes, trs, &
                            bs_m, bs_n, bs_k, sparsities, alpha, beta, &
                            n_loops, eps, retain_sparsity, always_checksum)

      CLASS(mp_comm_type), INTENT(IN)                     :: mp_group
      INTEGER, INTENT(IN)                                :: io_unit
      INTEGER, DIMENSION(:), INTENT(in)                  :: matrix_sizes
      LOGICAL, DIMENSION(2), INTENT(in)                  :: trs
      INTEGER, DIMENSION(:), POINTER                     :: bs_m, bs_n, bs_k
      REAL(kind=dp), DIMENSION(3), INTENT(in)            :: sparsities
      REAL(kind=dp), INTENT(in)                          :: alpha, beta
      INTEGER, INTENT(IN)                                :: n_loops
      REAL(kind=dp), INTENT(in)                          :: eps
      LOGICAL, INTENT(in)                                :: retain_sparsity, always_checksum

      CHARACTER(len=*), PARAMETER                        :: routineN = 'dbm_run_tests'

      INTEGER                                            :: bmax, bmin, handle
      INTEGER, CONTIGUOUS, DIMENSION(:), POINTER         :: col_dist_a, col_dist_b, col_dist_c, &
                                                            row_dist_a, row_dist_b, row_dist_c, &
                                                            sizes_k, sizes_m, sizes_n
      INTEGER, DIMENSION(2)                              :: npdims
      TYPE(dbm_distribution_obj)                         :: dist_a, dist_b, dist_c
      TYPE(dbm_type), TARGET                             :: matrix_a, matrix_b, matrix_c
      TYPE(mp_cart_type)                                 :: cart_group

      CALL timeset(routineN, handle)

      ! Create MPI processor grid.
      npdims(:) = 0
      CALL mp_dims_create(mp_group%num_pe, npdims)
      CALL cart_group%create(mp_group, 2, npdims)
      npdims = cart_group%num_pe_cart

      ! Initialize random number generator.
      randmat_counter = 12341313

      ! Create the row/column block sizes.
      NULLIFY (sizes_k, sizes_m, sizes_n)
      IF (ASSOCIATED(bs_m)) THEN
         bmin = MINVAL(bs_m(2::2))
         bmax = MAXVAL(bs_m(2::2))
         CALL generate_mixed_block_sizes(sizes_m, matrix_sizes(1), bs_m)
      ELSE
         CALL generate_mixed_block_sizes(sizes_m, matrix_sizes(1), (/1, 13, 2, 5/))
         bmin = 5; bmax = 13
      END IF
      IF (ASSOCIATED(bs_n)) THEN
         bmin = MIN(bmin, MINVAL(bs_n(2::2)))
         bmax = MAX(bmax, MAXVAL(bs_n(2::2)))
         CALL generate_mixed_block_sizes(sizes_n, matrix_sizes(2), bs_n)
      ELSE
         CALL generate_mixed_block_sizes(sizes_n, matrix_sizes(2), (/1, 13, 2, 5/))
         bmin = MIN(bmin, 5); bmax = MAX(bmax, 13)
      END IF
      IF (ASSOCIATED(bs_k)) THEN
         bmin = MIN(bmin, MINVAL(bs_k(2::2)))
         bmax = MAX(bmax, MAXVAL(bs_k(2::2)))
         CALL generate_mixed_block_sizes(sizes_k, matrix_sizes(3), bs_k)
      ELSE
         CALL generate_mixed_block_sizes(sizes_k, matrix_sizes(3), (/1, 13, 2, 5/))
         bmin = MIN(bmin, 5); bmax = MAX(bmax, 13)
      END IF

      ! Create Matrix C
      CALL generate_1d_dist(row_dist_c, SIZE(sizes_m), npdims(1), sizes_m)
      CALL generate_1d_dist(col_dist_c, SIZE(sizes_n), npdims(2), sizes_n)
      CALL dbm_distribution_new(dist_c, cart_group, row_dist_c, col_dist_c)
      CALL dbm_create(matrix_c, "Matrix C", dist_c, sizes_m, sizes_n)
      CALL fill_matrix(matrix_c, sparsity=sparsities(3), group=cart_group)
      CALL dbm_distribution_release(dist_c)

      ! Create Matrix A
      IF (trs(1)) THEN
         CALL generate_1d_dist(row_dist_a, SIZE(sizes_k), npdims(1), sizes_k)
         CALL generate_1d_dist(col_dist_a, SIZE(sizes_m), npdims(2), sizes_m)
         CALL dbm_distribution_new(dist_a, cart_group, row_dist_a, col_dist_a)
         CALL dbm_create(matrix_a, "Matrix A", dist_a, sizes_k, sizes_m)
         CALL fill_matrix(matrix_a, sparsity=sparsities(1), group=cart_group)
         DEALLOCATE (row_dist_a, col_dist_a)
      ELSE
         CALL generate_1d_dist(col_dist_a, SIZE(sizes_k), npdims(2), sizes_k)
         CALL dbm_distribution_new(dist_a, cart_group, row_dist_c, col_dist_a)
         CALL dbm_create(matrix_a, "Matrix A", dist_a, sizes_m, sizes_k)
         CALL fill_matrix(matrix_a, sparsity=sparsities(1), group=cart_group)
         DEALLOCATE (col_dist_a)
      END IF
      CALL dbm_distribution_release(dist_a)

      ! Create Matrix B
      IF (trs(2)) THEN
         CALL generate_1d_dist(row_dist_b, SIZE(sizes_n), npdims(1), sizes_n)
         CALL generate_1d_dist(col_dist_b, SIZE(sizes_k), npdims(2), sizes_k)
         CALL dbm_distribution_new(dist_b, cart_group, row_dist_b, col_dist_b)
         CALL dbm_create(matrix_b, "Matrix B", dist_b, sizes_n, sizes_k)
         CALL fill_matrix(matrix_b, sparsity=sparsities(2), group=cart_group)
         DEALLOCATE (row_dist_b, col_dist_b)
      ELSE
         CALL generate_1d_dist(row_dist_b, SIZE(sizes_k), npdims(1), sizes_k)
         CALL dbm_distribution_new(dist_b, cart_group, row_dist_b, col_dist_c)
         CALL dbm_create(matrix_b, "Matrix B", dist_b, sizes_k, sizes_n)
         CALL fill_matrix(matrix_b, sparsity=sparsities(2), group=cart_group)
         DEALLOCATE (row_dist_b)
      END IF
      CALL dbm_distribution_release(dist_b)
      DEALLOCATE (row_dist_c, col_dist_c, sizes_m, sizes_n, sizes_k)

      ! Prepare test parameters
      IF (io_unit > 0) THEN
         WRITE (io_unit, '(A,3(1X,I6),1X,A,2(1X,I5),1X,A,2(1X,L1))') &
            "Testing with sizes", matrix_sizes(1:3), &
            "min/max block sizes", bmin, bmax, "transposed?", trs(1:2)
      END IF

      CALL run_multiply_test(matrix_a, matrix_b, matrix_c, &
                             transa=trs(1), transb=trs(2), &
                             alpha=alpha, beta=beta, &
                             n_loops=n_loops, &
                             eps=eps, &
                             group=cart_group, &
                             io_unit=io_unit, &
                             always_checksum=always_checksum, &
                             retain_sparsity=retain_sparsity)

      CALL dbm_release(matrix_a)
      CALL dbm_release(matrix_b)
      CALL dbm_release(matrix_c)
      CALL cart_group%free()

      CALL timestop(handle)
   END SUBROUTINE dbm_run_tests

! **************************************************************************************************
!> \brief Runs the multiplication test.
!> \param matrix_a ...
!> \param matrix_b ...
!> \param matrix_c ...
!> \param transa ...
!> \param transb ...
!> \param alpha ...
!> \param beta ...
!> \param retain_sparsity ...
!> \param n_loops ...
!> \param eps ...
!> \param group ...
!> \param io_unit ...
!> \param always_checksum ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE run_multiply_test(matrix_a, matrix_b, matrix_c, transa, transb, alpha, beta, &
                                retain_sparsity, n_loops, eps, group, io_unit, always_checksum)
      TYPE(dbm_type), INTENT(in)                         :: matrix_a, matrix_b
      TYPE(dbm_type), INTENT(inout)                      :: matrix_c
      LOGICAL, INTENT(in)                                :: transa, transb
      REAL(kind=dp), INTENT(in)                          :: alpha, beta
      LOGICAL, INTENT(in)                                :: retain_sparsity
      INTEGER, INTENT(IN)                                :: n_loops
      REAL(kind=dp), INTENT(in)                          :: eps

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      INTEGER, INTENT(IN)                                :: io_unit
      LOGICAL, INTENT(in)                                :: always_checksum

      CHARACTER(len=*), PARAMETER                        :: routineN = 'run_multiply_test'

      INTEGER                                            :: handle, loop_iter
      INTEGER(kind=int_8)                                :: flop
      REAL(kind=dp)                                      :: cs, duration, flops_all, time_start
      TYPE(dbm_type)                                     :: matrix_c_orig

      CALL timeset(routineN, handle)

      CALL dbm_create_from_template(matrix_c_orig, "Original Matrix C", matrix_c)
      CALL dbm_copy(matrix_c_orig, matrix_c)

      ASSOCIATE (numnodes => group%num_pe)
         DO loop_iter = 1, n_loops
            CALL group%sync()
            time_start = omp_get_wtime()
            IF (eps < -0.0_dp) THEN
               CALL dbm_multiply(transa, transb, alpha, matrix_a, matrix_b, beta, matrix_c, &
                                 retain_sparsity=retain_sparsity, flop=flop)
            ELSE
               CALL dbm_multiply(transa, transb, alpha, matrix_a, matrix_b, beta, matrix_c, &
                                 retain_sparsity=retain_sparsity, flop=flop, filter_eps=eps)
            END IF
            duration = omp_get_wtime() - time_start

            CALL group%max(duration)
            CALL group%sum(flop)
            duration = MAX(duration, EPSILON(duration))  ! avoid division by zero
            flops_all = REAL(flop, KIND=dp)/duration/numnodes/(1024*1024)
            IF (io_unit .GT. 0) THEN
               WRITE (io_unit, '(A,I5,A,I5,A,F12.3,A,I9,A)') &
                  " loop ", loop_iter, " with ", numnodes, " MPI ranks: using ", &
                  duration, "s ", INT(flops_all), " Mflops/rank"
               CALL m_flush(io_unit)
            END IF

            IF (loop_iter .EQ. n_loops .OR. always_checksum) THEN
               cs = dbm_checksum(matrix_c)
               IF (io_unit > 0) THEN
                  WRITE (io_unit, *) "Final checksums", cs
               END IF
            END IF

            CALL dbm_copy(matrix_c, matrix_c_orig)
         END DO
      END ASSOCIATE

      CALL dbm_release(matrix_c_orig)
      CALL timestop(handle)
   END SUBROUTINE run_multiply_test

! **************************************************************************************************
!> \brief Fills give matrix with random blocks.
!> \param matrix ...
!> \param sparsity ...
!> \param group ...
! **************************************************************************************************
   SUBROUTINE fill_matrix(matrix, sparsity, group)
      TYPE(dbm_type), INTENT(inout)                      :: matrix
      REAL(kind=dp), INTENT(in)                          :: sparsity

      CLASS(mp_comm_type), INTENT(IN)                     :: group

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fill_matrix'

      INTEGER                                            :: block_node, col, handle, ncol, &
                                                            nrow, row
      INTEGER(KIND=int_8)                                :: counter, ele, increment, nmax
      INTEGER, DIMENSION(4)                              :: iseed, jseed
      INTEGER, DIMENSION(:), POINTER                     :: col_block_sizes, row_block_sizes
      REAL(kind=dp)                                      :: my_sparsity
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: block
      REAL(kind=dp), DIMENSION(1)                        :: value

      CALL timeset(routineN, handle)

      ! Check that the counter was initialised (or has not overflowed)
      CPASSERT(randmat_counter .NE. 0)
      ! the counter goes into the seed. Every new call gives a new random matrix
      randmat_counter = randmat_counter + 1

      IF (sparsity .GT. 1) THEN
         my_sparsity = sparsity/100.0
      ELSE
         my_sparsity = sparsity
      END IF

      row_block_sizes => dbm_get_row_block_sizes(matrix)
      col_block_sizes => dbm_get_col_block_sizes(matrix)
      nrow = SIZE(row_block_sizes)
      ncol = SIZE(col_block_sizes)
      nmax = INT(nrow, KIND=int_8)*INT(ncol, KIND=int_8)
      ele = -1
      counter = 0
      jseed = generate_larnv_seed(7, 42, 3, 42, randmat_counter)

      ASSOCIATE (mynode => group%mepos)
         DO
            ! find the next block to add, this is given by a geometrically distributed variable
            ! we number the blocks of the matrix and jump to the next one
            CALL dlarnv(1, jseed, 1, value)
            IF (my_sparsity > 0) THEN
               increment = 1 + FLOOR(LOG(value(1))/LOG(my_sparsity), KIND=int_8)
            ELSE
               increment = 1
            END IF
            ele = ele + increment
            IF (ele >= nmax) EXIT
            counter = counter + 1
            row = INT(ele/ncol) + 1
            col = INT(MOD(ele, INT(ncol, KIND=KIND(ele)))) + 1

            ! Only deal with the local blocks.
            CALL dbm_get_stored_coordinates(matrix, row, col, block_node)
            IF (block_node == mynode) THEN
               ! fill based on a block based seed, makes this the same values in parallel
               iseed = generate_larnv_seed(row, nrow, col, ncol, randmat_counter)
               ALLOCATE (block(row_block_sizes(row), col_block_sizes(col)))
               CALL dlarnv(1, iseed, SIZE(block), block)
               CALL dbm_put_block(matrix, row, col, block)
               DEALLOCATE (block)
            END IF
         END DO
      END ASSOCIATE

      CALL timestop(handle)
   END SUBROUTINE fill_matrix

! **************************************************************************************************
!> \brief Assigns given elements to bins. Uses given element_sizes for load balancing.
!> \param bin_dist Distribution of elements to bins
!> \param nelements Number of elements
!> \param nbins Number of bins
!> \param element_sizes sizes of elements
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE generate_1d_dist(bin_dist, nelements, nbins, element_sizes)
      INTEGER, DIMENSION(:), INTENT(OUT), POINTER        :: bin_dist
      INTEGER, INTENT(IN)                                :: nelements, nbins
      INTEGER, DIMENSION(:), INTENT(IN)                  :: element_sizes

      INTEGER                                            :: bin, i
      INTEGER, DIMENSION(nbins)                          :: bin_counts

      CPASSERT(SIZE(element_sizes) == nelements)
      ALLOCATE (bin_dist(nelements))

      bin_counts(:) = [(0, bin=0, nbins - 1)]
      DO i = 1, nelements
         bin = MINLOC(bin_counts, dim=1) ! greedy algorithm
         bin_dist(i) = bin - 1
         bin_counts(bin) = bin_counts(bin) + element_sizes(i)
      END DO
   END SUBROUTINE generate_1d_dist

! **************************************************************************************************
!> \brief Generates a block_sizes by "randomly" selecting from size_mix.
!> \param block_sizes ...
!> \param size_sum ...
!> \param size_mix ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE generate_mixed_block_sizes(block_sizes, size_sum, size_mix)
      INTEGER, DIMENSION(:), INTENT(inout), POINTER      :: block_sizes
      INTEGER, INTENT(in)                                :: size_sum
      INTEGER, DIMENSION(:), INTENT(in)                  :: size_mix

      INTEGER                                            :: block_size, current_sum, ipass, nblocks, &
                                                            nsize_mix, selector
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: mixer

      CPASSERT(.NOT. ASSOCIATED(block_sizes))
      nsize_mix = SIZE(size_mix)/2
      ALLOCATE (mixer(3, nsize_mix))

      ! 1st pass to compute nblocks and allocate block_sizes, 2nd pass to fill block_sizes.
      DO ipass = 1, 2
         mixer(1, :) = size_mix(1:nsize_mix*2 - 1:2)
         mixer(2, :) = size_mix(2:nsize_mix*2:2)
         mixer(3, :) = 1
         selector = 1
         nblocks = 0
         current_sum = 0
         DO WHILE (current_sum < size_sum)
            nblocks = nblocks + 1
            block_size = MIN(mixer(2, selector), size_sum - current_sum)
            IF (ipass == 2) THEN
               block_sizes(nblocks) = block_size
            END IF
            current_sum = current_sum + block_size
            mixer(3, selector) = mixer(3, selector) + 1
            IF (mixer(3, selector) > mixer(1, selector)) THEN
               mixer(3, selector) = 1
               selector = MOD(selector, nsize_mix) + 1
            END IF
         END DO
         IF (ipass == 1) THEN
            ALLOCATE (block_sizes(nblocks))
         END IF
      END DO

      current_sum = SUM(block_sizes)
      CPASSERT(current_sum == size_sum)
   END SUBROUTINE generate_mixed_block_sizes

! **************************************************************************************************
!> \brief Generate a seed respecting the lapack constraints,
!>        - values between 0..4095 (2**12-1)
!>        - iseed(4) odd
!>        also try to avoid iseed that are zero.
!>        Have but with a unique 2D mapping (irow,icol), and a 'mini-seed' ival
!>
!> \param irow 1..nrow
!> \param nrow ...
!> \param icol 1..ncol
!> \param ncol ...
!> \param ival mini-seed
!> \return a lapack compatible seed
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION generate_larnv_seed(irow, nrow, icol, ncol, ival) RESULT(iseed)

      INTEGER, INTENT(IN)                                :: irow, nrow, icol, ncol, ival
      INTEGER                                            :: iseed(4)

      INTEGER(KIND=int_8)                                :: map

      map = ((irow - 1 + icol*INT(nrow, int_8))*(1 + MODULO(ival, 2**16)))*2 + 1 + 0*ncol ! ncol used
      iseed(4) = INT(MODULO(map, 2_int_8**12)); map = map/2_int_8**12; ! keep odd
      iseed(3) = INT(MODULO(IEOR(map, 3541_int_8), 2_int_8**12)); map = map/2_int_8**12
      iseed(2) = INT(MODULO(IEOR(map, 1153_int_8), 2_int_8**12)); map = map/2_int_8**12
      iseed(1) = INT(MODULO(IEOR(map, 2029_int_8), 2_int_8**12)); map = map/2_int_8**12
   END FUNCTION generate_larnv_seed

END MODULE dbm_tests
