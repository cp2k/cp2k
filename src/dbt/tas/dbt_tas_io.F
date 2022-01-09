!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief tall-and-skinny matrices: Input / Output
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_tas_io
   USE dbm_api,                         ONLY: dbm_distribution_col_dist,&
                                              dbm_distribution_obj,&
                                              dbm_distribution_row_dist,&
                                              dbm_get_distribution
   USE dbt_tas_base,                    ONLY: dbt_tas_get_info,&
                                              dbt_tas_get_num_blocks,&
                                              dbt_tas_get_num_blocks_total,&
                                              dbt_tas_get_nze,&
                                              dbt_tas_get_nze_total,&
                                              dbt_tas_nblkcols_total,&
                                              dbt_tas_nblkrows_total
   USE dbt_tas_global,                  ONLY: dbt_tas_distribution,&
                                              dbt_tas_rowcol_data
   USE dbt_tas_split,                   ONLY: colsplit,&
                                              dbt_tas_get_split_info,&
                                              rowsplit
   USE dbt_tas_types,                   ONLY: dbt_tas_split_info,&
                                              dbt_tas_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_environ,&
                                              mp_max,&
                                              mp_sum
#include "../../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_tas_io'

   PUBLIC :: &
      dbt_tas_write_dist, &
      dbt_tas_write_matrix_info, &
      dbt_tas_write_split_info, &
      prep_output_unit

CONTAINS

! **************************************************************************************************
!> \brief Write basic infos of tall-and-skinny matrix:
!>        block dimensions, full dimensions, process grid dimensions
!> \param matrix ...
!> \param unit_nr ...
!> \param full_info Whether to print distribution and block size vectors
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_write_matrix_info(matrix, unit_nr, full_info)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN), OPTIONAL                      :: full_info

      INTEGER                                            :: unit_nr_prv
      INTEGER(KIND=int_8)                                :: nblkcols_total, nblkrows_total

      CLASS(dbt_tas_distribution), ALLOCATABLE :: proc_row_dist, proc_col_dist
      CLASS(dbt_tas_rowcol_data), ALLOCATABLE  :: row_blk_size, col_blk_size
      INTEGER(KIND=int_8)                      :: iblk
      CHARACTER(default_string_length)         :: name

      unit_nr_prv = prep_output_unit(unit_nr)
      IF (unit_nr_prv == 0) RETURN

      CALL dbt_tas_get_info(matrix, nblkrows_total=nblkrows_total, nblkcols_total=nblkcols_total, &
                            proc_row_dist=proc_row_dist, proc_col_dist=proc_col_dist, &
                            row_blk_size=row_blk_size, col_blk_size=col_blk_size, name=name)

      IF (unit_nr_prv > 0) THEN
         WRITE (unit_nr_prv, "(T2,A)") "GLOBAL INFO OF "//TRIM(name)
         WRITE (unit_nr_prv, "(T4,A,1X)", advance="no") "block dimensions:"
         WRITE (unit_nr_prv, "(I12,I12)", advance="no") nblkrows_total, nblkcols_total
         WRITE (unit_nr_prv, "(/T4,A,1X)", advance="no") "full dimensions:"
         WRITE (unit_nr_prv, "(I14,I14)", advance="no") row_blk_size%nfullrowcol, col_blk_size%nfullrowcol
         WRITE (unit_nr_prv, "(/T4,A,1X)", advance="no") "process grid dimensions:"
         WRITE (unit_nr_prv, "(I10,I10)", advance="no") proc_row_dist%nprowcol, proc_col_dist%nprowcol
         IF (PRESENT(full_info)) THEN
            IF (full_info) THEN
               WRITE (unit_nr_prv, '(/T4,A)', advance='no') "Block sizes:"
               WRITE (unit_nr_prv, '(/T8,A)', advance='no') 'Row:'
               DO iblk = 1, row_blk_size%nmrowcol
                  WRITE (unit_nr_prv, '(I4,1X)', advance='no') row_blk_size%data(iblk)
               END DO
               WRITE (unit_nr_prv, '(/T8,A)', advance='no') 'Column:'
               DO iblk = 1, col_blk_size%nmrowcol
                  WRITE (unit_nr_prv, '(I4,1X)', advance='no') col_blk_size%data(iblk)
               END DO
               WRITE (unit_nr_prv, '(/T4,A)', advance='no') "Block distribution:"
               WRITE (unit_nr_prv, '(/T8,A)', advance='no') 'Row:'
               DO iblk = 1, proc_row_dist%nmrowcol
                  WRITE (unit_nr_prv, '(I4,1X)', advance='no') proc_row_dist%dist(iblk)
               END DO
               WRITE (unit_nr_prv, '(/T8,A)', advance='no') 'Column:'
               DO iblk = 1, proc_col_dist%nmrowcol
                  WRITE (unit_nr_prv, '(I4,1X)', advance='no') proc_col_dist%dist(iblk)
               END DO

            END IF
         END IF
         WRITE (unit_nr_prv, *)
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief Write info on tall-and-skinny matrix distribution & load balance
!> \param matrix ...
!> \param unit_nr ...
!> \param full_info Whether to print subgroup DBM distribution
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_write_dist(matrix, unit_nr, full_info)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN), OPTIONAL                      :: full_info

      CHARACTER(default_string_length)                   :: name
      INTEGER :: icol, igroup, iproc, irow, mp_comm, mp_comm_group, nblock, ndbt_p_max, nelement, &
         nelement_p_max, ngroup, nproc, split_rowcol, unit_nr_prv
      INTEGER(KIND=int_8)                                :: ndbt_p_sum, ndbt_s, ndbt_s_max, &
                                                            ndbt_tot, nelement_p_sum, nelement_s, &
                                                            nelement_s_max
      INTEGER(KIND=int_8), DIMENSION(2)                  :: tmp_i8
      INTEGER, DIMENSION(2)                              :: tmp
      INTEGER, DIMENSION(:), POINTER                     :: coldist, rowdist
      REAL(KIND=dp)                                      :: occupation
      TYPE(dbm_distribution_obj)                         :: dist

      unit_nr_prv = prep_output_unit(unit_nr)
      IF (unit_nr_prv == 0) RETURN

      CALL dbt_tas_get_split_info(matrix%dist%info, mp_comm, ngroup, igroup, mp_comm_group, split_rowcol)
      CALL dbt_tas_get_info(matrix, name=name)
      CALL mp_environ(nproc, iproc, mp_comm)

      nblock = dbt_tas_get_num_blocks(matrix)
      nelement = dbt_tas_get_nze(matrix)

      ndbt_p_sum = dbt_tas_get_num_blocks_total(matrix)
      nelement_p_sum = dbt_tas_get_nze_total(matrix)

      tmp = (/nblock, nelement/)
      CALL mp_max(tmp, mp_comm)
      ndbt_p_max = tmp(1); nelement_p_max = tmp(2)

      ndbt_s = nblock
      nelement_s = nelement

      CALL mp_sum(ndbt_s, mp_comm_group)
      CALL mp_sum(nelement_s, mp_comm_group)

      tmp_i8 = (/ndbt_s, nelement_s/)
      CALL mp_max(tmp_i8, mp_comm)
      ndbt_s_max = tmp_i8(1); nelement_s_max = tmp_i8(2)

      ndbt_tot = dbt_tas_nblkrows_total(matrix)*dbt_tas_nblkcols_total(matrix)
      occupation = -1.0_dp
      IF (ndbt_tot .NE. 0) occupation = 100.0_dp*REAL(ndbt_p_sum, dp)/REAL(ndbt_tot, dp)

      dist = dbm_get_distribution(matrix%matrix)
      rowdist => dbm_distribution_row_dist(dist)
      coldist => dbm_distribution_col_dist(dist)

      IF (unit_nr_prv > 0) THEN
         WRITE (unit_nr_prv, "(T2,A)") &
            "DISTRIBUTION OF "//TRIM(name)
         WRITE (unit_nr_prv, "(T15,A,T68,I13)") "Number of non-zero blocks:", ndbt_p_sum
         WRITE (unit_nr_prv, "(T15,A,T75,F6.2)") "Percentage of non-zero blocks:", occupation
         WRITE (unit_nr_prv, "(T15,A,T68,I13)") "Average number of blocks per group:", (ndbt_p_sum + ngroup - 1)/ngroup
         WRITE (unit_nr_prv, "(T15,A,T68,I13)") "Maximum number of blocks per group:", ndbt_s_max
         WRITE (unit_nr_prv, "(T15,A,T68,I13)") "Average number of matrix elements per group:", (nelement_p_sum + ngroup - 1)/ngroup
         WRITE (unit_nr_prv, "(T15,A,T68,I13)") "Maximum number of matrix elements per group:", nelement_s_max
         WRITE (unit_nr_prv, "(T15,A,T68,I13)") "Average number of blocks per CPU:", (ndbt_p_sum + nproc - 1)/nproc
         WRITE (unit_nr_prv, "(T15,A,T68,I13)") "Maximum number of blocks per CPU:", ndbt_p_max
         WRITE (unit_nr_prv, "(T15,A,T68,I13)") "Average number of matrix elements per CPU:", (nelement_p_sum + nproc - 1)/nproc
         WRITE (unit_nr_prv, "(T15,A,T68,I13)") "Maximum number of matrix elements per CPU:", nelement_p_max
         IF (PRESENT(full_info)) THEN
            IF (full_info) THEN
               WRITE (unit_nr_prv, "(T15,A)") "Row distribution on subgroup:"
               WRITE (unit_nr_prv, '(T15)', advance='no')
               DO irow = 1, SIZE(rowdist)
                  WRITE (unit_nr_prv, '(I3, 1X)', advance='no') rowdist(irow)
               END DO
               WRITE (unit_nr_prv, "(/T15,A)") "Column distribution on subgroup:"
               WRITE (unit_nr_prv, '(T15)', advance='no')
               DO icol = 1, SIZE(coldist)
                  WRITE (unit_nr_prv, '(I3, 1X)', advance='no') coldist(icol)
               END DO
               WRITE (unit_nr_prv, *)
            END IF
         END IF
      END IF
   END SUBROUTINE

! **************************************************************************************************
!> \brief Print info on how matrix is split
!> \param info ...
!> \param unit_nr ...
!> \param name ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_write_split_info(info, unit_nr, name)
      TYPE(dbt_tas_split_info), INTENT(IN)               :: info
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: name

      CHARACTER(len=:), ALLOCATABLE                      :: name_prv
      INTEGER                                            :: groupsize, igroup, mp_comm, &
                                                            mp_comm_group, mynode, nsplit, &
                                                            numnodes, split_rowcol, unit_nr_prv
      INTEGER, DIMENSION(2)                              :: coord, dims, groupcoord, groupdims, &
                                                            pgrid_offset

      unit_nr_prv = prep_output_unit(unit_nr)
      IF (unit_nr_prv == 0) RETURN

      IF (PRESENT(name)) THEN
         ALLOCATE (name_prv, SOURCE=TRIM(name))
      ELSE
         ALLOCATE (name_prv, SOURCE="")
      END IF

      CALL dbt_tas_get_split_info(info, mp_comm, nsplit, igroup, mp_comm_group, split_rowcol, pgrid_offset)

      CALL mp_environ(numnodes, mynode, mp_comm)
      CALL mp_environ(numnodes, dims, coord, mp_comm)
      CALL mp_environ(groupsize, groupdims, groupcoord, mp_comm_group)

      IF (unit_nr_prv > 0) THEN
         SELECT CASE (split_rowcol)
         CASE (rowsplit)
            WRITE (unit_nr_prv, "(T4,A,I4,1X,A,I4)") name_prv//"splitting rows by factor", nsplit
         CASE (colsplit)
            WRITE (unit_nr_prv, "(T4,A,I4,1X,A,I4)") name_prv//"splitting columns by factor", nsplit
         END SELECT
         WRITE (unit_nr_prv, "(T4,A,I4,A1,I4)") name_prv//"global grid sizes:", dims(1), "x", dims(2)
      END IF

      IF (unit_nr_prv > 0) THEN
         WRITE (unit_nr_prv, "(T4,A,I4,A1,I4)") &
            name_prv//"grid sizes on subgroups:", &
            groupdims(1), "x", groupdims(2)
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param unit_nr ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION prep_output_unit(unit_nr) RESULT(unit_nr_out)
      INTEGER, INTENT(IN), OPTIONAL                      :: unit_nr
      INTEGER                                            :: unit_nr_out

      IF (PRESENT(unit_nr)) THEN
         unit_nr_out = unit_nr
      ELSE
         unit_nr_out = 0
      END IF

   END FUNCTION

END MODULE

