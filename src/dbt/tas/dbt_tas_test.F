!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief testing infrastructure for tall-and-skinny matrices
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_tas_test
   USE dbm_api,                         ONLY: &
        dbm_add, dbm_checksum, dbm_create, dbm_distribution_new, dbm_distribution_obj, &
        dbm_distribution_release, dbm_finalize, dbm_get_col_block_sizes, dbm_get_name, &
        dbm_get_row_block_sizes, dbm_maxabs, dbm_multiply, dbm_redistribute, dbm_release, &
        dbm_scale, dbm_type
   USE dbt_tas_base,                    ONLY: &
        dbt_tas_convert_to_dbm, dbt_tas_create, dbt_tas_distribution_new, dbt_tas_finalize, &
        dbt_tas_get_stored_coordinates, dbt_tas_info, dbt_tas_nblkcols_total, &
        dbt_tas_nblkrows_total, dbt_tas_put_block
   USE dbt_tas_global,                  ONLY: dbt_tas_blk_size_arb,&
                                              dbt_tas_default_distvec,&
                                              dbt_tas_dist_cyclic
   USE dbt_tas_mm,                      ONLY: dbt_tas_multiply
   USE dbt_tas_split,                   ONLY: dbt_tas_get_split_info,&
                                              dbt_tas_mp_comm
   USE dbt_tas_types,                   ONLY: dbt_tas_distribution_type,&
                                              dbt_tas_type
   USE dbt_tas_util,                    ONLY: generate_larnv_seed
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_cart_create,&
                                              mp_comm_free,&
                                              mp_environ
#include "../../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: &
      dbt_tas_benchmark_mm, &
      dbt_tas_checksum, &
      dbt_tas_random_bsizes, &
      dbt_tas_setup_test_matrix, &
      dbt_tas_test_mm, &
      dbt_tas_reset_randmat_seed

   INTEGER, SAVE :: randmat_counter = 0
   INTEGER, PARAMETER, PRIVATE :: rand_seed_init = 12341313

CONTAINS

! **************************************************************************************************
!> \brief Setup tall-and-skinny matrix for testing
!> \param matrix ...
!> \param mp_comm_out ...
!> \param mp_comm ...
!> \param nrows ...
!> \param ncols ...
!> \param rbsizes ...
!> \param cbsizes ...
!> \param dist_splitsize ...
!> \param name ...
!> \param sparsity ...
!> \param reuse_comm ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_setup_test_matrix(matrix, mp_comm_out, mp_comm, nrows, ncols, rbsizes, &
                                        cbsizes, dist_splitsize, name, sparsity, reuse_comm)
      TYPE(dbt_tas_type), INTENT(OUT)                    :: matrix
      INTEGER, INTENT(OUT)                               :: mp_comm_out
      INTEGER, INTENT(IN)                                :: mp_comm
      INTEGER(KIND=int_8), INTENT(IN)                    :: nrows, ncols
      INTEGER, DIMENSION(nrows), INTENT(IN)              :: rbsizes
      INTEGER, DIMENSION(ncols), INTENT(IN)              :: cbsizes
      INTEGER, DIMENSION(2), INTENT(IN)                  :: dist_splitsize
      CHARACTER(len=*), INTENT(IN)                       :: name
      REAL(KIND=dp), INTENT(IN)                          :: sparsity
      LOGICAL, INTENT(IN), OPTIONAL                      :: reuse_comm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_tas_setup_test_matrix'

      INTEGER                                            :: col_size, handle, max_col_size, max_nze, &
                                                            max_row_size, mynode, node_holds_blk, &
                                                            numnodes, nze, row_size
      INTEGER(KIND=int_8)                                :: col, col_s, ncol, nrow, row, row_s
      INTEGER, DIMENSION(2)                              :: pcoord, pdims
      INTEGER, DIMENSION(4)                              :: iseed, jseed
      LOGICAL                                            :: reuse_comm_prv, tr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: values
      REAL(KIND=dp), DIMENSION(1)                        :: rn
      TYPE(dbt_tas_blk_size_arb)                         :: cbsize_obj, rbsize_obj
      TYPE(dbt_tas_dist_cyclic)                          :: col_dist_obj, row_dist_obj
      TYPE(dbt_tas_distribution_type)                    :: dist

      ! we don't reserve blocks prior to putting them, so this time is meaningless and should not
      ! be considered in benchmark!
      CALL timeset(routineN, handle)

      ! Check that the counter was initialised (or has not overflowed)
      CPASSERT(randmat_counter .NE. 0)
      ! the counter goes into the seed. Every new call gives a new random matrix
      randmat_counter = randmat_counter + 1

      IF (PRESENT(reuse_comm)) THEN
         reuse_comm_prv = reuse_comm
      ELSE
         reuse_comm_prv = .FALSE.
      END IF

      IF (reuse_comm_prv) THEN
         mp_comm_out = mp_comm
      ELSE
         mp_comm_out = dbt_tas_mp_comm(mp_comm, nrows, ncols)
      END IF

      CALL mp_environ(numnodes, mynode, mp_comm_out)
      CALL mp_environ(numnodes, pdims, pcoord, mp_comm_out)

      row_dist_obj = dbt_tas_dist_cyclic(dist_splitsize(1), pdims(1), nrows)
      col_dist_obj = dbt_tas_dist_cyclic(dist_splitsize(2), pdims(2), ncols)

      rbsize_obj = dbt_tas_blk_size_arb(rbsizes)
      cbsize_obj = dbt_tas_blk_size_arb(cbsizes)

      CALL dbt_tas_distribution_new(dist, mp_comm_out, row_dist_obj, col_dist_obj)
      CALL dbt_tas_create(matrix, name, dist=dist, &
                          row_blk_size=rbsize_obj, col_blk_size=cbsize_obj, own_dist=.TRUE.)

      max_row_size = MAXVAL(rbsizes)
      max_col_size = MAXVAL(cbsizes)
      max_nze = max_row_size*max_col_size

      nrow = dbt_tas_nblkrows_total(matrix)
      ncol = dbt_tas_nblkcols_total(matrix)

      ALLOCATE (values(max_row_size, max_col_size))

      jseed = generate_larnv_seed(7, 42, 3, 42, randmat_counter)

      DO row = 1, dbt_tas_nblkrows_total(matrix)
         DO col = 1, dbt_tas_nblkcols_total(matrix)
            CALL dlarnv(1, jseed, 1, rn)
            IF (rn(1) .LT. sparsity) THEN
               tr = .FALSE.
               row_s = row; col_s = col
               CALL dbt_tas_get_stored_coordinates(matrix, row_s, col_s, node_holds_blk)

               IF (node_holds_blk .EQ. mynode) THEN
                  row_size = rbsize_obj%data(row_s)
                  col_size = cbsize_obj%data(col_s)
                  nze = row_size*col_size
                  iseed = generate_larnv_seed(INT(row_s), INT(nrow), INT(col_s), INT(ncol), randmat_counter)
                  CALL dlarnv(1, iseed, max_nze, values)
                  CALL dbt_tas_put_block(matrix, row_s, col_s, values(1:row_size, 1:col_size))
               END IF
            END IF
         END DO
      END DO

      CALL dbt_tas_finalize(matrix)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Benchmark routine. Due to random sparsity (as opposed to structured sparsity pattern),
!>        this may not be representative for actual applications.
!> \param transa ...
!> \param transb ...
!> \param transc ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param matrix_c ...
!> \param compare_dbm ...
!> \param filter_eps ...
!> \param io_unit ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_benchmark_mm(transa, transb, transc, matrix_a, matrix_b, matrix_c, compare_dbm, filter_eps, io_unit)

      LOGICAL, INTENT(IN)                                :: transa, transb, transc
      TYPE(dbt_tas_type), INTENT(INOUT)                  :: matrix_a, matrix_b, matrix_c
      LOGICAL, INTENT(IN)                                :: compare_dbm
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: filter_eps
      INTEGER, INTENT(IN), OPTIONAL                      :: io_unit

      INTEGER                                            :: comm_dbm, handle1, handle2, mp_comm
      INTEGER, CONTIGUOUS, DIMENSION(:), POINTER         :: cd_a, cd_b, cd_c, col_blk_size, rd_a, &
                                                            rd_b, rd_c, row_blk_size
      INTEGER, DIMENSION(2)                              :: myploc, npdims
      TYPE(dbm_distribution_obj)                         :: dist_a, dist_b, dist_c
      TYPE(dbm_type)                                     :: dbm_a, dbm_a_mm, dbm_b, dbm_b_mm, dbm_c, &
                                                            dbm_c_mm

      IF (PRESENT(io_unit)) THEN
      IF (io_unit > 0) THEN
         WRITE (io_unit, "(A)") "starting tall-and-skinny benchmark"
      END IF
      END IF
      CALL timeset("benchmark_tas_mm", handle1)
      CALL dbt_tas_multiply(transa, transb, transc, 1.0_dp, matrix_a, matrix_b, &
                            0.0_dp, matrix_c, &
                            filter_eps=filter_eps, unit_nr=io_unit)
      CALL timestop(handle1)
      IF (PRESENT(io_unit)) THEN
      IF (io_unit > 0) THEN
         WRITE (io_unit, "(A)") "tall-and-skinny benchmark completed"
      END IF
      END IF

      IF (compare_dbm) THEN
         CALL dbt_tas_convert_to_dbm(matrix_a, dbm_a)
         CALL dbt_tas_convert_to_dbm(matrix_b, dbm_b)
         CALL dbt_tas_convert_to_dbm(matrix_c, dbm_c)

         CALL dbt_tas_get_split_info(dbt_tas_info(matrix_a), mp_comm=mp_comm)
         npdims(:) = 0
         CALL mp_cart_create(mp_comm, 2, npdims, myploc, comm_dbm)

         ALLOCATE (rd_a(SIZE(dbm_get_row_block_sizes(dbm_a))))
         ALLOCATE (rd_b(SIZE(dbm_get_row_block_sizes(dbm_b))))
         ALLOCATE (rd_c(SIZE(dbm_get_row_block_sizes(dbm_c))))
         ALLOCATE (cd_a(SIZE(dbm_get_col_block_sizes(dbm_a))))
         ALLOCATE (cd_b(SIZE(dbm_get_col_block_sizes(dbm_b))))
         ALLOCATE (cd_c(SIZE(dbm_get_col_block_sizes(dbm_c))))

         CALL dbt_tas_default_distvec(INT(SIZE(dbm_get_row_block_sizes(dbm_a))), &
                                      npdims(1), dbm_get_row_block_sizes(dbm_a), rd_a)
         CALL dbt_tas_default_distvec(INT(SIZE(dbm_get_col_block_sizes(dbm_a))), &
                                      npdims(2), dbm_get_col_block_sizes(dbm_a), cd_a)
         CALL dbt_tas_default_distvec(INT(SIZE(dbm_get_row_block_sizes(dbm_b))), &
                                      npdims(1), dbm_get_row_block_sizes(dbm_b), rd_b)
         CALL dbt_tas_default_distvec(INT(SIZE(dbm_get_col_block_sizes(dbm_b))), &
                                      npdims(2), dbm_get_col_block_sizes(dbm_b), cd_b)
         CALL dbt_tas_default_distvec(INT(SIZE(dbm_get_row_block_sizes(dbm_c))), &
                                      npdims(1), dbm_get_row_block_sizes(dbm_c), rd_c)
         CALL dbt_tas_default_distvec(INT(SIZE(dbm_get_col_block_sizes(dbm_c))), &
                                      npdims(2), dbm_get_col_block_sizes(dbm_c), cd_c)

         CALL dbm_distribution_new(dist_a, comm_dbm, rd_a, cd_a)
         CALL dbm_distribution_new(dist_b, comm_dbm, rd_b, cd_b)
         CALL dbm_distribution_new(dist_c, comm_dbm, rd_c, cd_c)
         DEALLOCATE (rd_a, rd_b, rd_c, cd_a, cd_b, cd_c)

         row_blk_size => dbm_get_row_block_sizes(dbm_a)
         col_blk_size => dbm_get_col_block_sizes(dbm_a)
         CALL dbm_create(matrix=dbm_a_mm, &
                         name=dbm_get_name(dbm_a), &
                         dist=dist_a, &
                         row_block_sizes=row_blk_size, col_block_sizes=col_blk_size)
         row_blk_size => dbm_get_row_block_sizes(dbm_b)
         col_blk_size => dbm_get_col_block_sizes(dbm_b)
         CALL dbm_create(matrix=dbm_b_mm, &
                         name=dbm_get_name(dbm_b), &
                         dist=dist_b, &
                         row_block_sizes=row_blk_size, col_block_sizes=col_blk_size)
         row_blk_size => dbm_get_row_block_sizes(dbm_c)
         col_blk_size => dbm_get_col_block_sizes(dbm_c)
         CALL dbm_create(matrix=dbm_c_mm, &
                         name=dbm_get_name(dbm_c), &
                         dist=dist_c, &
                         row_block_sizes=row_blk_size, col_block_sizes=col_blk_size)

         CALL dbm_finalize(dbm_a_mm)
         CALL dbm_finalize(dbm_b_mm)
         CALL dbm_finalize(dbm_c_mm)

         CALL dbm_redistribute(dbm_a, dbm_a_mm)
         CALL dbm_redistribute(dbm_b, dbm_b_mm)
         IF (PRESENT(io_unit)) THEN
         IF (io_unit > 0) THEN
            WRITE (io_unit, "(A)") "starting dbm benchmark"
         END IF
         END IF
         CALL timeset("benchmark_block_mm", handle2)
         CALL dbm_multiply(transa, transb, 1.0_dp, dbm_a_mm, dbm_b_mm, &
                           0.0_dp, dbm_c_mm, filter_eps=filter_eps)
         CALL timestop(handle2)
         IF (PRESENT(io_unit)) THEN
         IF (io_unit > 0) THEN
            WRITE (io_unit, "(A)") "dbm benchmark completed"
         END IF
         END IF

         CALL dbm_release(dbm_a)
         CALL dbm_release(dbm_b)
         CALL dbm_release(dbm_c)
         CALL dbm_release(dbm_a_mm)
         CALL dbm_release(dbm_b_mm)
         CALL dbm_release(dbm_c_mm)
         CALL dbm_distribution_release(dist_a)
         CALL dbm_distribution_release(dist_b)
         CALL dbm_distribution_release(dist_c)

         CALL mp_comm_free(comm_dbm)
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief Test tall-and-skinny matrix multiplication for accuracy
!> \param transa ...
!> \param transb ...
!> \param transc ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param matrix_c ...
!> \param filter_eps ...
!> \param unit_nr ...
!> \param log_verbose ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_test_mm(transa, transb, transc, matrix_a, matrix_b, matrix_c, filter_eps, unit_nr, log_verbose)
      LOGICAL, INTENT(IN)                                :: transa, transb, transc
      TYPE(dbt_tas_type), INTENT(INOUT)                  :: matrix_a, matrix_b, matrix_c
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: filter_eps
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN), OPTIONAL                      :: log_verbose

      REAL(KIND=dp), PARAMETER                           :: test_tol = 1.0E-10_dp

      INTEGER                                            :: comm_dbm, io_unit, mp_comm, mynode, &
                                                            numnodes
      INTEGER, CONTIGUOUS, DIMENSION(:), POINTER         :: cd_a, cd_b, cd_c, rd_a, rd_b, rd_c
      INTEGER, DIMENSION(2)                              :: myploc, npdims
      LOGICAL                                            :: transa_prv, transb_prv
      REAL(KIND=dp)                                      :: norm, rc_cs, sq_cs
      TYPE(dbm_distribution_obj)                         :: dist_a, dist_b, dist_c
      TYPE(dbm_type)                                     :: dbm_a, dbm_a_mm, dbm_b, dbm_b_mm, dbm_c, &
                                                            dbm_c_mm, dbm_c_mm_check

!
! TODO: Dedup with code in dbt_tas_benchmark_mm.
!

      CALL dbt_tas_get_split_info(dbt_tas_info(matrix_a), mp_comm=mp_comm)
      CALL mp_environ(numnodes, mynode, mp_comm)
      io_unit = -1
      IF (mynode .EQ. 0) io_unit = unit_nr

      CALL dbt_tas_multiply(transa, transb, transc, 1.0_dp, matrix_a, matrix_b, &
                            0.0_dp, matrix_c, &
                            filter_eps=filter_eps, unit_nr=io_unit, log_verbose=log_verbose, optimize_dist=.TRUE.)

      CALL dbt_tas_convert_to_dbm(matrix_a, dbm_a)
      CALL dbt_tas_convert_to_dbm(matrix_b, dbm_b)
      CALL dbt_tas_convert_to_dbm(matrix_c, dbm_c)

      npdims(:) = 0
      CALL mp_cart_create(mp_comm, 2, npdims, myploc, comm_dbm)

      ALLOCATE (rd_a(SIZE(dbm_get_row_block_sizes(dbm_a))))
      ALLOCATE (rd_b(SIZE(dbm_get_row_block_sizes(dbm_b))))
      ALLOCATE (rd_c(SIZE(dbm_get_row_block_sizes(dbm_c))))
      ALLOCATE (cd_a(SIZE(dbm_get_col_block_sizes(dbm_a))))
      ALLOCATE (cd_b(SIZE(dbm_get_col_block_sizes(dbm_b))))
      ALLOCATE (cd_c(SIZE(dbm_get_col_block_sizes(dbm_c))))

      CALL dbt_tas_default_distvec(INT(SIZE(dbm_get_row_block_sizes(dbm_a))), &
                                   npdims(1), dbm_get_row_block_sizes(dbm_a), rd_a)
      CALL dbt_tas_default_distvec(INT(SIZE(dbm_get_col_block_sizes(dbm_a))), &
                                   npdims(2), dbm_get_col_block_sizes(dbm_a), cd_a)
      CALL dbt_tas_default_distvec(INT(SIZE(dbm_get_row_block_sizes(dbm_b))), &
                                   npdims(1), dbm_get_row_block_sizes(dbm_b), rd_b)
      CALL dbt_tas_default_distvec(INT(SIZE(dbm_get_col_block_sizes(dbm_b))), &
                                   npdims(2), dbm_get_col_block_sizes(dbm_b), cd_b)
      CALL dbt_tas_default_distvec(INT(SIZE(dbm_get_row_block_sizes(dbm_c))), &
                                   npdims(1), dbm_get_row_block_sizes(dbm_c), rd_c)
      CALL dbt_tas_default_distvec(INT(SIZE(dbm_get_col_block_sizes(dbm_c))), &
                                   npdims(2), dbm_get_col_block_sizes(dbm_c), cd_c)

      CALL dbm_distribution_new(dist_a, comm_dbm, rd_a, cd_a)
      CALL dbm_distribution_new(dist_b, comm_dbm, rd_b, cd_b)
      CALL dbm_distribution_new(dist_c, comm_dbm, rd_c, cd_c)
      DEALLOCATE (rd_a, rd_b, rd_c, cd_a, cd_b, cd_c)

      CALL dbm_create(matrix=dbm_a_mm, name="matrix a", dist=dist_a, &
                      row_block_sizes=dbm_get_row_block_sizes(dbm_a), &
                      col_block_sizes=dbm_get_col_block_sizes(dbm_a))

      CALL dbm_create(matrix=dbm_b_mm, name="matrix b", dist=dist_b, &
                      row_block_sizes=dbm_get_row_block_sizes(dbm_b), &
                      col_block_sizes=dbm_get_col_block_sizes(dbm_b))

      CALL dbm_create(matrix=dbm_c_mm, name="matrix c", dist=dist_c, &
                      row_block_sizes=dbm_get_row_block_sizes(dbm_c), &
                      col_block_sizes=dbm_get_col_block_sizes(dbm_c))

      CALL dbm_create(matrix=dbm_c_mm_check, name="matrix c check", &
                      dist=dist_c, &
                      row_block_sizes=dbm_get_row_block_sizes(dbm_c), &
                      col_block_sizes=dbm_get_col_block_sizes(dbm_c))

      CALL dbm_finalize(dbm_a_mm)
      CALL dbm_finalize(dbm_b_mm)
      CALL dbm_finalize(dbm_c_mm)
      CALL dbm_finalize(dbm_c_mm_check)

      CALL dbm_redistribute(dbm_a, dbm_a_mm)
      CALL dbm_redistribute(dbm_b, dbm_b_mm)
      CALL dbm_redistribute(dbm_c, dbm_c_mm_check)

      transa_prv = transa; transb_prv = transb

      IF (.NOT. transc) THEN
         CALL dbm_multiply(transa_prv, transb_prv, 1.0_dp, &
                           dbm_a_mm, dbm_b_mm, &
                           0.0_dp, dbm_c_mm, filter_eps=filter_eps)
      ELSE
         transa_prv = .NOT. transa_prv
         transb_prv = .NOT. transb_prv
         CALL dbm_multiply(transb_prv, transa_prv, 1.0_dp, &
                           dbm_b_mm, dbm_a_mm, &
                           0.0_dp, dbm_c_mm, filter_eps=filter_eps)
      END IF

      sq_cs = dbm_checksum(dbm_c_mm)
      rc_cs = dbm_checksum(dbm_c_mm_check)
      CALL dbm_scale(dbm_c_mm_check, -1.0_dp)
      CALL dbm_add(dbm_c_mm_check, dbm_c_mm)
      norm = dbm_maxabs(dbm_c_mm_check)

      IF (io_unit > 0) THEN
      IF (ABS(norm) .GT. test_tol) THEN
         WRITE (io_unit, '(A, A, A, A, A, 1X, A)') TRIM(dbm_get_name(matrix_a%matrix)), &
            transa, ' X ', &
            TRIM(dbm_get_name(matrix_b%matrix)), &
            transb, 'failed!'
         WRITE (io_unit, "(A,1X,E9.2,1X,E9.2)") "checksums", sq_cs, rc_cs
         WRITE (io_unit, "(A,1X,E9.2)") "difference norm", norm
         CPABORT("")
      ELSE
         WRITE (io_unit, '(A, A, A, A, A, 1X, A)') TRIM(dbm_get_name(matrix_a%matrix)), &
            transa, ' X ', &
            TRIM(dbm_get_name(matrix_b%matrix)), &
            transb, 'passed!'
         WRITE (io_unit, "(A,1X,E9.2,1X,E9.2)") "checksums", sq_cs, rc_cs
         WRITE (io_unit, "(A,1X,E9.2)") "difference norm", norm
      END IF
      END IF

      CALL dbm_release(dbm_a)
      CALL dbm_release(dbm_a_mm)
      CALL dbm_release(dbm_b)
      CALL dbm_release(dbm_b_mm)
      CALL dbm_release(dbm_c)
      CALL dbm_release(dbm_c_mm)
      CALL dbm_release(dbm_c_mm_check)

      CALL dbm_distribution_release(dist_a)
      CALL dbm_distribution_release(dist_b)
      CALL dbm_distribution_release(dist_c)

      CALL mp_comm_free(comm_dbm)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Calculate checksum of tall-and-skinny matrix consistent with dbm_checksum
!> \param matrix ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_tas_checksum(matrix)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      REAL(KIND=dp)                                      :: dbt_tas_checksum

      TYPE(dbm_type)                                     :: dbm_m

      CALL dbt_tas_convert_to_dbm(matrix, dbm_m)
      dbt_tas_checksum = dbm_checksum(dbm_m)
      CALL dbm_release(dbm_m)
   END FUNCTION

! **************************************************************************************************
!> \brief Create random block sizes
!> \param sizes ...
!> \param repeat ...
!> \param dbt_sizes ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_random_bsizes(sizes, repeat, dbt_sizes)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes
      INTEGER, INTENT(IN)                                :: repeat
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: dbt_sizes

      INTEGER                                            :: d, size_i

      DO d = 1, SIZE(dbt_sizes)
         size_i = MOD((d - 1)/repeat, SIZE(sizes)) + 1
         dbt_sizes(d) = sizes(size_i)
      END DO
   END SUBROUTINE

! **************************************************************************************************
!> \brief Reset the seed used for generating random matrices to default value
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_reset_randmat_seed()
      randmat_counter = rand_seed_init
   END SUBROUTINE

END MODULE
