!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief CSR matrix type, DBCSR-CSR conversions, and CSR file I/O built on cp_dbcsr_api
! **************************************************************************************************
MODULE cp_dbcsr_csr
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_complete_redistribute, dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, &
        dbcsr_distribution_get, dbcsr_distribution_new, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_get_info, dbcsr_get_readonly_block_p, dbcsr_has_symmetry, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_readonly_start, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_release, dbcsr_set, &
        dbcsr_type, dbcsr_type_no_symmetry, dbcsr_valid_index
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_comm_null,&
                                              mp_comm_type
   USE util,                            ONLY: sort
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_csr'

   INTEGER, PARAMETER, PUBLIC  :: dbcsr_csr_type_real_8 = 1
   INTEGER, PARAMETER, PRIVATE :: dbcsr_csr_type_complex_8 = 2

   INTEGER, PARAMETER, PUBLIC  :: dbcsr_csr_dbcsr_blkrow_dist = 1
   INTEGER, PARAMETER, PRIVATE :: dbcsr_csr_eqrow_ceil_dist = 2
   INTEGER, PARAMETER, PUBLIC  :: dbcsr_csr_eqrow_floor_dist = 3

   PUBLIC :: dbcsr_csr_type
   PUBLIC :: dbcsr_csr_p_type
   PUBLIC :: dbcsr_csr_create
   PUBLIC :: dbcsr_csr_destroy
   PUBLIC :: dbcsr_csr_create_from_dbcsr
   PUBLIC :: dbcsr_convert_dbcsr_to_csr
   PUBLIC :: dbcsr_convert_csr_to_dbcsr
   PUBLIC :: dbcsr_csr_create_and_convert_complex
   PUBLIC :: dbcsr_csr_print_sparsity
   PUBLIC :: dbcsr_csr_write

   TYPE csr_mapping_data
      PRIVATE
      TYPE(dbcsr_type) :: brd_mat = dbcsr_type()
      INTEGER          :: node_first_row = 1
   END TYPE csr_mapping_data

   TYPE csr_data_area_type
      REAL(KIND=dp), DIMENSION(:), POINTER    :: r_dp => NULL()
      COMPLEX(KIND=dp), DIMENSION(:), POINTER :: c_dp => NULL()
      INTEGER                                 :: data_type = -1
   END TYPE csr_data_area_type

   TYPE dbcsr_csr_type
      INTEGER                                :: nrows_total = 0
      INTEGER                                :: ncols_total = 0
      INTEGER                                :: nze_local = 0
      INTEGER                                :: nrows_local = 0
      TYPE(mp_comm_type)                     :: mp_group = mp_comm_null
      INTEGER(KIND=int_8)                    :: nze_total = 0_int_8
      INTEGER, DIMENSION(:), POINTER         :: rowptr_local => NULL()
      INTEGER, DIMENSION(:), POINTER         :: colind_local => NULL()
      INTEGER, DIMENSION(:), POINTER         :: nzerow_local => NULL()
      TYPE(csr_data_area_type)               :: nzval_local = csr_data_area_type()
      TYPE(csr_mapping_data)                 :: dbcsr_mapping = csr_mapping_data()
      LOGICAL                                :: has_mapping = .FALSE.
      LOGICAL                                :: valid = .FALSE.
      LOGICAL                                :: has_indices = .FALSE.
   END TYPE dbcsr_csr_type

   TYPE dbcsr_csr_p_type
      TYPE(dbcsr_csr_type), POINTER :: csr_mat => NULL()
   END TYPE dbcsr_csr_p_type

   INTERFACE dbcsr_csr_create
      MODULE PROCEDURE csr_create_new, csr_create_template
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Create a new CSR matrix and allocate internal storage (excluding DBCSR mapping)
!> \param csr_mat ...
!> \param nrows_total ...
!> \param ncols_total ...
!> \param nze_total ...
!> \param nze_local ...
!> \param nrows_local ...
!> \param mp_group ...
!> \param data_type ...
! **************************************************************************************************
   SUBROUTINE csr_create_new(csr_mat, nrows_total, ncols_total, nze_total, &
                             nze_local, nrows_local, mp_group, data_type)
      TYPE(dbcsr_csr_type), INTENT(OUT)                  :: csr_mat
      INTEGER, INTENT(IN)                                :: nrows_total, ncols_total
      INTEGER(KIND=int_8), INTENT(IN)                    :: nze_total
      INTEGER, INTENT(IN)                                :: nze_local, nrows_local, mp_group
      INTEGER, INTENT(IN), OPTIONAL                      :: data_type

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'csr_create_new'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (nrows_total < nrows_local) THEN
         CPABORT("local number of rows must not exceed total number of rows")
      END IF

      IF (nze_total < INT(nze_local, KIND=int_8)) THEN
         CPABORT("local number of non-zero elements must not exceed total number of non-zero elements")
      END IF

      IF (INT(nrows_total, KIND=int_8)*INT(ncols_total, KIND=int_8) < nze_total) THEN
         CPABORT("Total number of non-zero elements must not exceed total matrix size")
      END IF

      IF (INT(nrows_local, KIND=int_8)*INT(ncols_total, KIND=int_8) < INT(nze_local, KIND=int_8)) THEN
         CPABORT("Local number of non-zero elements must not exceed local matrix size")
      END IF

      csr_mat%ncols_total = ncols_total
      csr_mat%nrows_total = nrows_total
      csr_mat%nze_total = nze_total
      csr_mat%nze_local = nze_local
      csr_mat%nrows_local = nrows_local

      ALLOCATE (csr_mat%colind_local(nze_local))
      csr_mat%colind_local(:) = 0
      ALLOCATE (csr_mat%rowptr_local(nrows_local + 1))
      csr_mat%rowptr_local(:) = 1
      ALLOCATE (csr_mat%nzerow_local(nrows_local))
      csr_mat%nzerow_local(:) = 0

      IF (PRESENT(data_type)) THEN
         csr_mat%nzval_local%data_type = data_type
      ELSE
         csr_mat%nzval_local%data_type = dbcsr_csr_type_real_8
      END IF

      SELECT CASE (csr_mat%nzval_local%data_type)
      CASE (dbcsr_csr_type_real_8)
         ALLOCATE (csr_mat%nzval_local%r_dp(nze_local))
         csr_mat%nzval_local%r_dp(:) = 0.0_dp
      CASE (dbcsr_csr_type_complex_8)
         ALLOCATE (csr_mat%nzval_local%c_dp(nze_local))
         csr_mat%nzval_local%c_dp(:) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      CASE DEFAULT
         CPABORT("Invalid CSR matrix data type")
      END SELECT

      CALL csr_mat%mp_group%set_handle(mp_group)

      csr_mat%valid = .TRUE.
      csr_mat%has_mapping = .FALSE.
      csr_mat%has_indices = .FALSE.

      CALL timestop(handle)
   END SUBROUTINE csr_create_new

! **************************************************************************************************
!> \brief Create a new CSR matrix from an existing template CSR matrix
!> \param matrix_b ...
!> \param matrix_a ...
! **************************************************************************************************
   SUBROUTINE csr_create_template(matrix_b, matrix_a)
      TYPE(dbcsr_csr_type), INTENT(OUT)                  :: matrix_b
      TYPE(dbcsr_csr_type), INTENT(IN)                   :: matrix_a

      CHARACTER(LEN=*), PARAMETER :: routineN = 'csr_create_template'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (.NOT. matrix_a%valid) THEN
         CPABORT("Source CSR matrix must be created.")
      END IF

      CALL csr_create_new(matrix_b, matrix_a%nrows_total, matrix_a%ncols_total, &
                          matrix_a%nze_total, matrix_a%nze_local, matrix_a%nrows_local, &
                          matrix_a%mp_group%get_handle(), matrix_a%nzval_local%data_type)

      matrix_b%mp_group = matrix_a%mp_group
      matrix_b%has_mapping = matrix_a%has_mapping
      matrix_b%has_indices = matrix_a%has_indices

      IF (matrix_a%has_indices) THEN
         matrix_b%rowptr_local(:) = matrix_a%rowptr_local(:)
         matrix_b%nzerow_local(:) = matrix_a%nzerow_local(:)
         matrix_b%colind_local(:) = matrix_a%colind_local(:)
      END IF

      IF (matrix_a%has_mapping) THEN
         matrix_b%dbcsr_mapping%node_first_row = matrix_a%dbcsr_mapping%node_first_row
         CALL dbcsr_copy(matrix_b%dbcsr_mapping%brd_mat, matrix_a%dbcsr_mapping%brd_mat)
      END IF

      CALL timestop(handle)
   END SUBROUTINE csr_create_template

! **************************************************************************************************
!> \brief Destroy a CSR matrix and release all allocated memory
!> \param csr_mat ...
! **************************************************************************************************
   SUBROUTINE dbcsr_csr_destroy(csr_mat)
      TYPE(dbcsr_csr_type), INTENT(INOUT)                :: csr_mat

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbcsr_csr_destroy'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (.NOT. csr_mat%valid) THEN
         CPABORT("CSR matrix must be created before destroying it.")
      END IF

      IF (ASSOCIATED(csr_mat%rowptr_local)) DEALLOCATE (csr_mat%rowptr_local)
      IF (ASSOCIATED(csr_mat%nzerow_local)) DEALLOCATE (csr_mat%nzerow_local)
      IF (ASSOCIATED(csr_mat%colind_local)) DEALLOCATE (csr_mat%colind_local)

      IF (csr_mat%has_mapping) THEN
         CALL dbcsr_release(csr_mat%dbcsr_mapping%brd_mat)
      END IF

      IF (ASSOCIATED(csr_mat%nzval_local%r_dp)) DEALLOCATE (csr_mat%nzval_local%r_dp)
      IF (ASSOCIATED(csr_mat%nzval_local%c_dp)) DEALLOCATE (csr_mat%nzval_local%c_dp)

      csr_mat%has_mapping = .FALSE.
      csr_mat%valid = .FALSE.
      csr_mat%has_indices = .FALSE.
      csr_mat%nzval_local%data_type = -1
      csr_mat%dbcsr_mapping%node_first_row = 1

      CALL timestop(handle)
   END SUBROUTINE dbcsr_csr_destroy

! **************************************************************************************************
!> \brief Create CSR matrix including DBCSR mapping from an arbitrary DBCSR matrix
!> \param dbcsr_mat ...
!> \param csr_mat ...
!> \param dist_format ...
!> \param csr_sparsity ...
!> \param numnodes ...
! **************************************************************************************************
   SUBROUTINE dbcsr_csr_create_from_dbcsr(dbcsr_mat, csr_mat, dist_format, csr_sparsity, numnodes)
      TYPE(dbcsr_type), INTENT(IN)                       :: dbcsr_mat
      TYPE(dbcsr_csr_type), INTENT(OUT)                  :: csr_mat
      INTEGER, INTENT(IN)                                :: dist_format
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: csr_sparsity
      INTEGER, INTENT(IN), OPTIONAL                      :: numnodes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_csr_create_from_dbcsr'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      CALL csr_create_from_dbcsr_impl(dbcsr_mat, csr_mat, dist_format, dbcsr_csr_type_real_8, &
                                      csr_sparsity=csr_sparsity, numnodes=numnodes)
      CALL timestop(handle)
   END SUBROUTINE dbcsr_csr_create_from_dbcsr

! **************************************************************************************************
!> \brief Internal helper creating a real or complex CSR structure from a DBCSR matrix
!> \param dbcsr_mat ...
!> \param csr_mat ...
!> \param dist_format ...
!> \param data_type ...
!> \param csr_sparsity ...
!> \param numnodes ...
! **************************************************************************************************
   SUBROUTINE csr_create_from_dbcsr_impl(dbcsr_mat, csr_mat, dist_format, data_type, &
                                         csr_sparsity, numnodes)
      TYPE(dbcsr_type), INTENT(IN)                       :: dbcsr_mat
      TYPE(dbcsr_csr_type), INTENT(OUT)                  :: csr_mat
      INTEGER, INTENT(IN)                                :: dist_format, data_type
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: csr_sparsity
      INTEGER, INTENT(IN), OPTIONAL                      :: numnodes

      INTEGER :: blkcol, blkrow, col_offset, col_size, dbcsr_numnodes, local_row, m, mp_handle, n, &
         nfullcols_total, nfullrows_total, node_first_row, nrows_local, num_p, nze_local, nze_r, &
         pos, r, r_end, r_start, row_offset, row_size, sp_first_row, sp_nrows_local
      INTEGER(KIND=int_8)                                :: nze_total
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nzerow_priv, nzerow_tmp, perm, row_pos
      LOGICAL                                            :: equal_dist, floor_dist, found, use_sp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_ptr, sp_ptr
      TYPE(dbcsr_distribution_type)                      :: dist_current
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: brd_mat, csr_sparsity_brd, &
                                                            csr_sparsity_nosym, dbcsr_mat_nosym
      TYPE(mp_comm_type)                                 :: mp_group

      IF (.NOT. dbcsr_valid_index(dbcsr_mat)) THEN
         CPABORT("Invalid DBCSR matrix")
      END IF

      SELECT CASE (dist_format)
      CASE (dbcsr_csr_dbcsr_blkrow_dist)
         equal_dist = .FALSE.
         floor_dist = .FALSE.
      CASE (dbcsr_csr_eqrow_ceil_dist)
         equal_dist = .TRUE.
         floor_dist = .FALSE.
      CASE (dbcsr_csr_eqrow_floor_dist)
         equal_dist = .TRUE.
         floor_dist = .TRUE.
      CASE DEFAULT
         CPABORT("Invalid CSR distribution format")
      END SELECT

      IF (dbcsr_has_symmetry(dbcsr_mat)) THEN
         CALL dbcsr_desymmetrize(dbcsr_mat, dbcsr_mat_nosym)
      ELSE
         CALL dbcsr_copy(dbcsr_mat_nosym, dbcsr_mat)
      END IF

      CALL dbcsr_get_info(dbcsr_mat_nosym, distribution=dist_current, &
                          nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
      CALL dbcsr_distribution_get(dist_current, group=mp_handle, numnodes=dbcsr_numnodes)
      CALL mp_group%set_handle(mp_handle)

      IF (PRESENT(numnodes)) THEN
         IF (numnodes > dbcsr_numnodes) THEN
            CPABORT("Number of nodes used for CSR matrix must not exceed total number of nodes")
         END IF
         num_p = numnodes
      ELSE
         num_p = dbcsr_numnodes
      END IF

      CALL dbcsr_create_brd(dbcsr_mat_nosym, brd_mat, equal_dist, floor_dist, num_p, &
                            node_first_row, nrows_local)

      use_sp = PRESENT(csr_sparsity)
      IF (use_sp) THEN
         IF (dbcsr_has_symmetry(csr_sparsity)) THEN
            CALL dbcsr_desymmetrize(csr_sparsity, csr_sparsity_nosym)
         ELSE
            CALL dbcsr_copy(csr_sparsity_nosym, csr_sparsity)
         END IF
         CALL dbcsr_create_brd(csr_sparsity_nosym, csr_sparsity_brd, equal_dist, floor_dist, &
                               num_p, sp_first_row, sp_nrows_local)
         CPASSERT(sp_first_row == node_first_row .AND. sp_nrows_local == nrows_local)
         CALL dbcsr_release(csr_sparsity_nosym)
      END IF

      ALLOCATE (nzerow_tmp(nrows_local))
      nzerow_tmp(:) = 0

      ! Pass 1: count non-zero elements per local row
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(brd_mat, csr_sparsity_brd, use_sp, node_first_row, nrows_local, nzerow_tmp) &
!$OMP PRIVATE(iter, blkrow, blkcol, blk_ptr, sp_ptr, row_size, col_size, row_offset, &
!$OMP         col_offset, found, m, n, local_row, nzerow_priv)
      NULLIFY (blk_ptr, sp_ptr)
      ALLOCATE (nzerow_priv(nrows_local))
      nzerow_priv(:) = 0
      CALL dbcsr_iterator_readonly_start(iter, brd_mat)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, blkrow, blkcol, blk_ptr, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)
         IF (use_sp) THEN
            CALL dbcsr_get_readonly_block_p(csr_sparsity_brd, blkrow, blkcol, sp_ptr, found)
            IF (.NOT. found) CYCLE
            DO n = 1, col_size
               DO m = 1, row_size
                  IF (INT(sp_ptr(m, n)) /= 0) THEN
                     local_row = row_offset - node_first_row + m
                     nzerow_priv(local_row) = nzerow_priv(local_row) + 1
                  END IF
               END DO
            END DO
         ELSE
            DO m = 1, row_size
               local_row = row_offset - node_first_row + m
               nzerow_priv(local_row) = nzerow_priv(local_row) + col_size
            END DO
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP CRITICAL(csr_nzerow_crit)
      nzerow_tmp(:) = nzerow_tmp(:) + nzerow_priv(:)
!$OMP END CRITICAL(csr_nzerow_crit)
      DEALLOCATE (nzerow_priv)
!$OMP END PARALLEL

      nze_local = SUM(nzerow_tmp)
      nze_total = INT(nze_local, KIND=int_8)
      CALL mp_group%sum(nze_total)

      CALL csr_create_new(csr_mat, nfullrows_total, nfullcols_total, nze_total, &
                          nze_local, nrows_local, mp_handle, data_type)

      csr_mat%nzerow_local(:) = nzerow_tmp(:)
      csr_mat%rowptr_local(1) = 1
      DO r = 1, nrows_local
         csr_mat%rowptr_local(r + 1) = csr_mat%rowptr_local(r) + nzerow_tmp(r)
      END DO

      ! Pass 2: collect column indices per row and sort them in ascending order
      ALLOCATE (row_pos(nrows_local))
      IF (nrows_local > 0) row_pos(:) = csr_mat%rowptr_local(1:nrows_local)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(brd_mat, csr_sparsity_brd, use_sp, node_first_row, row_pos, csr_mat) &
!$OMP PRIVATE(iter, blkrow, blkcol, blk_ptr, sp_ptr, row_size, col_size, row_offset, &
!$OMP         col_offset, found, m, n, local_row, pos)
      NULLIFY (blk_ptr, sp_ptr)
      CALL dbcsr_iterator_readonly_start(iter, brd_mat)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, blkrow, blkcol, blk_ptr, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)
         IF (use_sp) THEN
            CALL dbcsr_get_readonly_block_p(csr_sparsity_brd, blkrow, blkcol, sp_ptr, found)
            IF (.NOT. found) CYCLE
            DO n = 1, col_size
               DO m = 1, row_size
                  IF (INT(sp_ptr(m, n)) /= 0) THEN
                     local_row = row_offset - node_first_row + m
!$OMP ATOMIC CAPTURE
                     pos = row_pos(local_row)
                     row_pos(local_row) = row_pos(local_row) + 1
!$OMP END ATOMIC
                     csr_mat%colind_local(pos) = col_offset + n - 1
                  END IF
               END DO
            END DO
         ELSE
            DO n = 1, col_size
               DO m = 1, row_size
                  local_row = row_offset - node_first_row + m
!$OMP ATOMIC CAPTURE
                  pos = row_pos(local_row)
                  row_pos(local_row) = row_pos(local_row) + 1
!$OMP END ATOMIC
                  csr_mat%colind_local(pos) = col_offset + n - 1
               END DO
            END DO
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      IF (nrows_local > 0) THEN
         ALLOCATE (perm(MAX(1, MAXVAL(nzerow_tmp))))
         DO r = 1, nrows_local
            nze_r = nzerow_tmp(r)
            IF (nze_r > 1) THEN
               r_start = csr_mat%rowptr_local(r)
               r_end = r_start + nze_r - 1
               CALL sort(csr_mat%colind_local(r_start:r_end), nze_r, perm)
            END IF
         END DO
         DEALLOCATE (perm)
      END IF

      DEALLOCATE (nzerow_tmp, row_pos)
      IF (use_sp) CALL dbcsr_release(csr_sparsity_brd)
      CALL dbcsr_release(dbcsr_mat_nosym)

      csr_mat%dbcsr_mapping%brd_mat = brd_mat
      csr_mat%dbcsr_mapping%node_first_row = node_first_row
      csr_mat%has_mapping = .TRUE.
      csr_mat%has_indices = .TRUE.
   END SUBROUTINE csr_create_from_dbcsr_impl

! **************************************************************************************************
!> \brief Convert a DBCSR matrix to an existing CSR matrix
!> \param dbcsr_mat ...
!> \param csr_mat ...
! **************************************************************************************************
   SUBROUTINE dbcsr_convert_dbcsr_to_csr(dbcsr_mat, csr_mat)
      TYPE(dbcsr_type), INTENT(IN)                       :: dbcsr_mat
      TYPE(dbcsr_csr_type), INTENT(INOUT)                :: csr_mat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_convert_dbcsr_to_csr'

      INTEGER                                            :: blkcol, blkrow, col_offset, col_size, &
                                                            handle, idx, local_row, m, n, &
                                                            node_first_row, r_end, r_start, &
                                                            row_offset, row_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_ptr
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: dbcsr_mat_nosym, tmp_brd

      CALL timeset(routineN, handle)

      IF (.NOT. dbcsr_valid_index(dbcsr_mat)) THEN
         CPABORT("Invalid DBCSR matrix")
      END IF
      IF (.NOT. csr_mat%has_mapping) THEN
         CPABORT("CSR_mat must contain mapping to DBCSR format")
      END IF
      IF (csr_mat%nzval_local%data_type /= dbcsr_csr_type_real_8) THEN
         CPABORT("DBCSR and CSR matrix must have same type")
      END IF

      CALL dbcsr_create(tmp_brd, template=csr_mat%dbcsr_mapping%brd_mat)
      IF (dbcsr_has_symmetry(dbcsr_mat)) THEN
         CALL dbcsr_desymmetrize(dbcsr_mat, dbcsr_mat_nosym)
         CALL dbcsr_complete_redistribute(dbcsr_mat_nosym, tmp_brd)
         CALL dbcsr_release(dbcsr_mat_nosym)
      ELSE
         CALL dbcsr_complete_redistribute(dbcsr_mat, tmp_brd)
      END IF

      csr_mat%nzval_local%r_dp(:) = 0.0_dp
      node_first_row = csr_mat%dbcsr_mapping%node_first_row

!$OMP PARALLEL DEFAULT(NONE) SHARED(tmp_brd, csr_mat, node_first_row) &
!$OMP PRIVATE(iter, blkrow, blkcol, blk_ptr, row_size, col_size, row_offset, col_offset, &
!$OMP         m, n, local_row, r_start, r_end, idx)
      NULLIFY (blk_ptr)
      CALL dbcsr_iterator_readonly_start(iter, tmp_brd)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, blkrow, blkcol, blk_ptr, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)
         DO m = 1, row_size
            local_row = row_offset - node_first_row + m
            r_start = csr_mat%rowptr_local(local_row)
            r_end = csr_mat%rowptr_local(local_row + 1) - 1
            IF (r_start > r_end) CYCLE
            idx = find_first_ge(csr_mat%colind_local, r_start, r_end, col_offset)
            DO WHILE (idx <= r_end)
               IF (csr_mat%colind_local(idx) >= col_offset + col_size) EXIT
               n = csr_mat%colind_local(idx) - col_offset + 1
               csr_mat%nzval_local%r_dp(idx) = blk_ptr(m, n)
               idx = idx + 1
            END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL dbcsr_release(tmp_brd)
      CALL timestop(handle)
   END SUBROUTINE dbcsr_convert_dbcsr_to_csr

! **************************************************************************************************
!> \brief Convert a CSR matrix back to an allocated DBCSR matrix
!> \param dbcsr_mat ...
!> \param csr_mat ...
! **************************************************************************************************
   SUBROUTINE dbcsr_convert_csr_to_dbcsr(dbcsr_mat, csr_mat)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: dbcsr_mat
      TYPE(dbcsr_csr_type), INTENT(INOUT)                :: csr_mat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_convert_csr_to_dbcsr'

      INTEGER                                            :: blkcol, blkrow, col_offset, col_size, &
                                                            handle, idx, local_row, m, n, &
                                                            node_first_row, r_end, r_start, &
                                                            row_offset, row_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_ptr
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      IF (.NOT. dbcsr_valid_index(dbcsr_mat)) THEN
         CPABORT("Invalid DBCSR matrix")
      END IF
      IF (.NOT. csr_mat%has_mapping) THEN
         CPABORT("CSR_mat must contain mapping to DBCSR format")
      END IF
      IF (csr_mat%nzval_local%data_type /= dbcsr_csr_type_real_8) THEN
         CPABORT("DBCSR and CSR matrix must have same type")
      END IF

      CALL dbcsr_set(csr_mat%dbcsr_mapping%brd_mat, 0.0_dp)
      node_first_row = csr_mat%dbcsr_mapping%node_first_row

!$OMP PARALLEL DEFAULT(NONE) SHARED(csr_mat, node_first_row) &
!$OMP PRIVATE(iter, blkrow, blkcol, blk_ptr, row_size, col_size, row_offset, col_offset, &
!$OMP         m, n, local_row, r_start, r_end, idx)
      NULLIFY (blk_ptr)
      CALL dbcsr_iterator_start(iter, csr_mat%dbcsr_mapping%brd_mat)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, blkrow, blkcol, blk_ptr, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)
         DO m = 1, row_size
            local_row = row_offset - node_first_row + m
            r_start = csr_mat%rowptr_local(local_row)
            r_end = csr_mat%rowptr_local(local_row + 1) - 1
            IF (r_start > r_end) CYCLE
            idx = find_first_ge(csr_mat%colind_local, r_start, r_end, col_offset)
            DO WHILE (idx <= r_end)
               IF (csr_mat%colind_local(idx) >= col_offset + col_size) EXIT
               n = csr_mat%colind_local(idx) - col_offset + 1
               blk_ptr(m, n) = csr_mat%nzval_local%r_dp(idx)
               idx = idx + 1
            END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL dbcsr_complete_redistribute(csr_mat%dbcsr_mapping%brd_mat, dbcsr_mat)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_convert_csr_to_dbcsr

! **************************************************************************************************
!> \brief Create and populate a complex CSR matrix from real and imaginary DBCSR matrices
!> \param rmatrix ...
!> \param imatrix ...
!> \param csr_mat ...
!> \param dist_format ...
! **************************************************************************************************
   SUBROUTINE dbcsr_csr_create_and_convert_complex(rmatrix, imatrix, csr_mat, dist_format)
      TYPE(dbcsr_type), INTENT(IN)                       :: rmatrix, imatrix
      TYPE(dbcsr_csr_type), INTENT(INOUT)                :: csr_mat
      INTEGER, INTENT(IN)                                :: dist_format

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_csr_create_and_convert_complex'

      INTEGER :: blkcol, blkrow, col_offset, col_size, handle, idx, ipart, local_row, m, n, &
         node_first_row, r_end, r_start, row_offset, row_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_ptr
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: brd_part, imat_nosym, rmat_nosym, &
                                                            union_mat

      CALL timeset(routineN, handle)

      IF (dbcsr_has_symmetry(rmatrix)) THEN
         CALL dbcsr_desymmetrize(rmatrix, rmat_nosym)
      ELSE
         CALL dbcsr_copy(rmat_nosym, rmatrix)
      END IF

      IF (dbcsr_has_symmetry(imatrix)) THEN
         CALL dbcsr_desymmetrize(imatrix, imat_nosym)
      ELSE
         CALL dbcsr_copy(imat_nosym, imatrix)
      END IF

      CALL dbcsr_copy(union_mat, rmat_nosym)
      CALL dbcsr_add(union_mat, imat_nosym, 1.0_dp, 1.0_dp)

      CALL csr_create_from_dbcsr_impl(union_mat, csr_mat, dist_format, dbcsr_csr_type_complex_8)
      CALL dbcsr_release(union_mat)

      csr_mat%nzval_local%c_dp(:) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      node_first_row = csr_mat%dbcsr_mapping%node_first_row

      CALL dbcsr_create(brd_part, template=csr_mat%dbcsr_mapping%brd_mat)

      DO ipart = 1, 2
         IF (ipart == 1) THEN
            CALL dbcsr_complete_redistribute(rmat_nosym, brd_part)
         ELSE
            CALL dbcsr_complete_redistribute(imat_nosym, brd_part)
         END IF

!$OMP PARALLEL DEFAULT(NONE) SHARED(brd_part, csr_mat, node_first_row, ipart) &
!$OMP PRIVATE(iter, blkrow, blkcol, blk_ptr, row_size, col_size, row_offset, col_offset, &
!$OMP         m, n, local_row, r_start, r_end, idx)
         NULLIFY (blk_ptr)
         CALL dbcsr_iterator_readonly_start(iter, brd_part)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, blkrow, blkcol, blk_ptr, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)
            DO m = 1, row_size
               local_row = row_offset - node_first_row + m
               r_start = csr_mat%rowptr_local(local_row)
               r_end = csr_mat%rowptr_local(local_row + 1) - 1
               IF (r_start > r_end) CYCLE
               idx = find_first_ge(csr_mat%colind_local, r_start, r_end, col_offset)
               DO WHILE (idx <= r_end)
                  IF (csr_mat%colind_local(idx) >= col_offset + col_size) EXIT
                  n = csr_mat%colind_local(idx) - col_offset + 1
                  IF (ipart == 1) THEN
                     csr_mat%nzval_local%c_dp(idx) = csr_mat%nzval_local%c_dp(idx) + &
                                                     CMPLX(blk_ptr(m, n), 0.0_dp, KIND=dp)
                  ELSE
                     csr_mat%nzval_local%c_dp(idx) = csr_mat%nzval_local%c_dp(idx) + &
                                                     CMPLX(0.0_dp, blk_ptr(m, n), KIND=dp)
                  END IF
                  idx = idx + 1
               END DO
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL
      END DO

      CALL dbcsr_release(brd_part)
      CALL dbcsr_release(rmat_nosym)
      CALL dbcsr_release(imat_nosym)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_csr_create_and_convert_complex

! **************************************************************************************************
!> \brief Print CSR sparsity statistics
!> \param csr_mat ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE dbcsr_csr_print_sparsity(csr_mat, unit_nr)
      TYPE(dbcsr_csr_type), INTENT(IN)                   :: csr_mat
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_csr_print_sparsity'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: nze_percentage

      CALL timeset(routineN, handle)

      IF (.NOT. csr_mat%valid) THEN
         CPABORT("CSR matrix must be created first")
      END IF

      IF (csr_mat%nrows_total > 0 .AND. csr_mat%ncols_total > 0) THEN
         nze_percentage = 100.0_dp*(REAL(csr_mat%nze_total, KIND=dp) &
                                    /REAL(csr_mat%nrows_total, KIND=dp)) &
                          /REAL(csr_mat%ncols_total, KIND=dp)
      ELSE
         nze_percentage = 0.0_dp
      END IF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, "(T15,A,T68,I13)") "Number of  CSR non-zero elements:", csr_mat%nze_total
         WRITE (unit_nr, "(T15,A,T75,F6.2)") "Percentage CSR non-zero elements:", nze_percentage
      END IF

      CALL timestop(handle)
   END SUBROUTINE dbcsr_csr_print_sparsity

! **************************************************************************************************
!> \brief Write a CSR matrix to file (formatted or unformatted)
!> \param csr_mat ...
!> \param unit_nr ...
!> \param upper_triangle ...
!> \param threshold ...
!> \param binary ...
! **************************************************************************************************
   SUBROUTINE dbcsr_csr_write(csr_mat, unit_nr, upper_triangle, threshold, binary)
      TYPE(dbcsr_csr_type), INTENT(IN)                   :: csr_mat
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN), OPTIONAL                      :: upper_triangle
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: threshold
      LOGICAL, INTENT(IN), OPTIONAL                      :: binary

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbcsr_csr_write'

      CHARACTER(LEN=default_string_length)               :: data_format
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: nzval_to_master_c_dp
      INTEGER                                            :: handle, i, ii, k, l, m, mynode, &
                                                            nrows_ii, numnodes, nze_ii, rowind, &
                                                            src, tag1, tag2, tag3
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: colind_to_master, nzerow_to_master, &
                                                            sizes_numrowlocal, sizes_nzelocal
      LOGICAL                                            :: bin, ut
      REAL(KIND=dp)                                      :: thld
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: nzval_to_master_r_dp

      CALL timeset(routineN, handle)

      ut = .FALSE.
      IF (PRESENT(upper_triangle)) ut = upper_triangle

      thld = 0.0_dp
      IF (PRESENT(threshold)) thld = threshold

      bin = .FALSE.
      IF (PRESENT(binary)) bin = binary

      IF (.NOT. csr_mat%valid) THEN
         CPABORT("can not write invalid CSR matrix")
      END IF

      numnodes = csr_mat%mp_group%num_pe
      mynode = csr_mat%mp_group%mepos

      ALLOCATE (sizes_nzelocal(numnodes))
      ALLOCATE (sizes_numrowlocal(numnodes))

      CALL csr_mat%mp_group%allgather(csr_mat%nze_local, sizes_nzelocal)
      CALL csr_mat%mp_group%allgather(csr_mat%nrows_local, sizes_numrowlocal)

      DO i = 0, numnodes - 1
         ii = i
         nze_ii = sizes_nzelocal(ii + 1)
         nrows_ii = sizes_numrowlocal(ii + 1)

         IF (mynode == 0) THEN
            ALLOCATE (colind_to_master(nze_ii))
            ALLOCATE (nzerow_to_master(nrows_ii))

            SELECT CASE (csr_mat%nzval_local%data_type)
            CASE (dbcsr_csr_type_real_8)
               data_format = "(2(I8),E23.14E3)"
               ALLOCATE (nzval_to_master_r_dp(nze_ii))
            CASE (dbcsr_csr_type_complex_8)
               data_format = "(2(I8),2(E23.14E3))"
               ALLOCATE (nzval_to_master_c_dp(nze_ii))
            END SELECT

            IF (ii == 0) THEN
               colind_to_master(:) = csr_mat%colind_local(:)
               nzerow_to_master(:) = csr_mat%nzerow_local(:)
               SELECT CASE (csr_mat%nzval_local%data_type)
               CASE (dbcsr_csr_type_real_8)
                  nzval_to_master_r_dp(:) = csr_mat%nzval_local%r_dp(:)
               CASE (dbcsr_csr_type_complex_8)
                  nzval_to_master_c_dp(:) = csr_mat%nzval_local%c_dp(:)
               END SELECT
            ELSE
               IF (nze_ii > 0) THEN
                  src = ii
                  tag1 = 0
                  CALL csr_mat%mp_group%recv(colind_to_master, src, tag1)
               END IF
               IF (nrows_ii > 0) THEN
                  src = ii
                  tag2 = 1
                  CALL csr_mat%mp_group%recv(nzerow_to_master, src, tag2)
               END IF
               IF (nze_ii > 0) THEN
                  src = ii
                  tag3 = 2
                  SELECT CASE (csr_mat%nzval_local%data_type)
                  CASE (dbcsr_csr_type_real_8)
                     CALL csr_mat%mp_group%recv(nzval_to_master_r_dp, src, tag3)
                  CASE (dbcsr_csr_type_complex_8)
                     CALL csr_mat%mp_group%recv(nzval_to_master_c_dp, src, tag3)
                  END SELECT
               END IF
            END IF
         END IF

         IF (mynode == ii .AND. ii /= 0) THEN
            IF (nze_ii > 0) THEN
               tag1 = 0
               CALL csr_mat%mp_group%send(csr_mat%colind_local, 0, tag1)
            END IF
            IF (nrows_ii > 0) THEN
               tag2 = 1
               CALL csr_mat%mp_group%send(csr_mat%nzerow_local, 0, tag2)
            END IF
            IF (nze_ii > 0) THEN
               tag3 = 2
               SELECT CASE (csr_mat%nzval_local%data_type)
               CASE (dbcsr_csr_type_real_8)
                  CALL csr_mat%mp_group%send(csr_mat%nzval_local%r_dp, 0, tag3)
               CASE (dbcsr_csr_type_complex_8)
                  CALL csr_mat%mp_group%send(csr_mat%nzval_local%c_dp, 0, tag3)
               END SELECT
            END IF
         END IF

         IF (mynode == 0) THEN
            IF (unit_nr > 0) THEN
               m = 0
               DO k = 1, nrows_ii
                  rowind = k + SUM(sizes_numrowlocal(1:ii))
                  DO l = 1, nzerow_to_master(k)
                     m = m + 1
                     IF ((.NOT. ut) .OR. (rowind <= colind_to_master(m))) THEN
                        SELECT CASE (csr_mat%nzval_local%data_type)
                        CASE (dbcsr_csr_type_real_8)
                           IF (ABS(nzval_to_master_r_dp(m)) >= thld) THEN
                              IF (bin) THEN
                                 WRITE (unit_nr) rowind, colind_to_master(m), nzval_to_master_r_dp(m)
                              ELSE
                                 WRITE (unit_nr, data_format) rowind, colind_to_master(m), &
                                    nzval_to_master_r_dp(m)
                              END IF
                           END IF
                        CASE (dbcsr_csr_type_complex_8)
                           IF (ABS(nzval_to_master_c_dp(m)) >= thld) THEN
                              IF (bin) THEN
                                 WRITE (unit_nr) rowind, colind_to_master(m), nzval_to_master_c_dp(m)
                              ELSE
                                 WRITE (unit_nr, data_format) rowind, colind_to_master(m), &
                                    nzval_to_master_c_dp(m)
                              END IF
                           END IF
                        END SELECT
                     END IF
                  END DO
               END DO
            END IF

            DEALLOCATE (colind_to_master, nzerow_to_master)
            SELECT CASE (csr_mat%nzval_local%data_type)
            CASE (dbcsr_csr_type_real_8)
               DEALLOCATE (nzval_to_master_r_dp)
            CASE (dbcsr_csr_type_complex_8)
               DEALLOCATE (nzval_to_master_c_dp)
            END SELECT
         END IF
      END DO

      DEALLOCATE (sizes_nzelocal, sizes_numrowlocal)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_csr_write

! **************************************************************************************************
!> \brief Convert a non-symmetric DBCSR matrix into a block-row-distributed (BRD) DBCSR matrix
!> \param dbcsr_mat ...
!> \param brd_mat ...
!> \param equal_dist ...
!> \param floor_dist ...
!> \param numnodes ...
!> \param node_first_row ...
!> \param nrows_local ...
! **************************************************************************************************
   SUBROUTINE dbcsr_create_brd(dbcsr_mat, brd_mat, equal_dist, floor_dist, numnodes, &
                               node_first_row, nrows_local)
      TYPE(dbcsr_type), INTENT(IN)                       :: dbcsr_mat
      TYPE(dbcsr_type), INTENT(OUT)                      :: brd_mat
      LOGICAL, INTENT(IN)                                :: equal_dist, floor_dist
      INTEGER, INTENT(IN)                                :: numnodes
      INTEGER, INTENT(OUT)                               :: node_first_row, nrows_local

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbcsr_create_brd'

      CHARACTER(LEN=default_string_length)               :: matrix_name
      INTEGER :: end_ind, handle, i, k, l, m, mp_group_handle, mynode, nblkcols_total, &
         nblkrows_total, nfullrows_local, nfullrows_total, node_size, numnodes_total, split_row, &
         start_ind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: rdist_tmp, row_blk_size_new_tmp
      INTEGER, CONTIGUOUS, DIMENSION(:), POINTER         :: cdist, col_blk_size, rdist, &
                                                            row_blk_size, row_blk_size_new
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      REAL(KIND=dp)                                      :: chunk_size
      TYPE(dbcsr_distribution_type)                      :: dist_current, dist_new

      CALL timeset(routineN, handle)

      NULLIFY (row_blk_size, col_blk_size, rdist, cdist, row_blk_size_new, pgrid)
      CALL dbcsr_get_info(dbcsr_mat, &
                          nblkrows_total=nblkrows_total, &
                          nblkcols_total=nblkcols_total, &
                          nfullrows_total=nfullrows_total, &
                          row_blk_size=row_blk_size, &
                          col_blk_size=col_blk_size, &
                          name=matrix_name, &
                          distribution=dist_current)

      CALL dbcsr_distribution_get(dist_current, group=mp_group_handle, &
                                  mynode=mynode, numnodes=numnodes_total)

      ALLOCATE (cdist(nblkcols_total))
      cdist(:) = 0

      IF (equal_dist) THEN
         IF (floor_dist) THEN
            nfullrows_local = FLOOR(REAL(nfullrows_total, KIND=dp)/REAL(numnodes, KIND=dp))
         ELSE
            nfullrows_local = CEILING(REAL(nfullrows_total, KIND=dp)/REAL(numnodes, KIND=dp))
         END IF

         ALLOCATE (rdist_tmp(nblkrows_total + numnodes))
         ALLOCATE (row_blk_size_new_tmp(nblkrows_total + numnodes))

         k = 0
         m = 0
         node_size = nfullrows_local
         IF (node_size > 0) THEN
            DO l = 1, nblkrows_total
               split_row = row_blk_size(l)
               DO WHILE (split_row >= node_size)
                  k = k + 1
                  m = m + 1
                  row_blk_size_new_tmp(k) = node_size
                  rdist_tmp(k) = m - 1
                  split_row = split_row - node_size
                  node_size = nfullrows_local
                  IF (floor_dist .AND. (m == numnodes - 1)) THEN
                     node_size = nfullrows_total - (numnodes - 1)*node_size
                  END IF
               END DO
               IF (split_row > 0) THEN
                  k = k + 1
                  row_blk_size_new_tmp(k) = split_row
                  rdist_tmp(k) = m
                  node_size = node_size - split_row
               END IF
            END DO
         ELSE
            rdist_tmp(1:nblkrows_total) = numnodes - 1
            row_blk_size_new_tmp(1:nblkrows_total) = row_blk_size(1:nblkrows_total)
            k = nblkrows_total
         END IF

         ALLOCATE (row_blk_size_new(k))
         row_blk_size_new(:) = row_blk_size_new_tmp(1:k)
         ALLOCATE (rdist(k))
         rdist(:) = rdist_tmp(1:k)
         DEALLOCATE (rdist_tmp, row_blk_size_new_tmp)
      ELSE
         ALLOCATE (rdist(nblkrows_total))
         rdist(:) = 0
         IF (numnodes > nblkrows_total) THEN
            DO i = 1, nblkrows_total
               rdist(i) = i - 1
            END DO
         ELSE
            chunk_size = REAL(nblkrows_total, KIND=dp)/REAL(numnodes, KIND=dp)
            DO i = 0, numnodes - 1
               start_ind = NINT(REAL(i, KIND=dp)*chunk_size) + 1
               end_ind = NINT(REAL(i + 1, KIND=dp)*chunk_size)
               IF (i == numnodes - 1) end_ind = nblkrows_total
               IF (start_ind <= end_ind) rdist(start_ind:end_ind) = i
            END DO
         END IF
         row_blk_size_new => row_blk_size
      END IF

      node_first_row = 1 + SUM(row_blk_size_new, MASK=(rdist < mynode))
      nrows_local = SUM(row_blk_size_new, MASK=(rdist == mynode))

      ALLOCATE (pgrid(0:numnodes_total - 1, 0:0))
      DO i = 0, numnodes_total - 1
         pgrid(i, 0) = i
      END DO

      CALL dbcsr_distribution_new(dist_new, group=mp_group_handle, pgrid=pgrid, &
                                  row_dist=rdist, col_dist=cdist)

      CALL dbcsr_create(brd_mat, TRIM(matrix_name)//" row-block distributed", &
                        dist_new, dbcsr_type_no_symmetry, row_blk_size_new, col_blk_size)
      CALL dbcsr_complete_redistribute(dbcsr_mat, brd_mat)

      CALL dbcsr_distribution_release(dist_new)
      DEALLOCATE (pgrid, rdist, cdist)
      IF (equal_dist) DEALLOCATE (row_blk_size_new)

      CALL timestop(handle)
   END SUBROUTINE dbcsr_create_brd

! **************************************************************************************************
!> \brief Binary search returning the smallest index i in [lo, hi] with arr(i) >= val, or hi + 1
!> \param arr ...
!> \param lo ...
!> \param hi ...
!> \param val ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION find_first_ge(arr, lo, hi, val) RESULT(res)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: arr
      INTEGER, INTENT(IN)                                :: lo, hi, val
      INTEGER                                            :: res

      INTEGER                                            :: l, mid, r

      l = lo
      r = hi + 1
      DO WHILE (l < r)
         mid = (l + r)/2
         IF (arr(mid) < val) THEN
            l = mid + 1
         ELSE
            r = mid
         END IF
      END DO
      res = l
   END FUNCTION find_first_ge

END MODULE cp_dbcsr_csr
