!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   [DESCRIPTION OF THE MODULE]
!> \author  [AUTHOR]
!> \date    [DATE]
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - [ADD MODIFICATION HISTORY or CREATION]
! *****************************************************************************
MODULE cp_bcsr_operations
  USE cp_bcsr_types,                   ONLY: &
       cp_bcsr_index_lookup, cp_bcsr_slot_col_blk_begin, &
       cp_bcsr_slot_col_blk_end, cp_bcsr_slot_row_blk_begin, &
       cp_bcsr_slot_row_blk_end, cp_bcsr_type, cp_bcsr_work_type, fetch_env, &
       fetch_info_empty, fetch_info_normal, fetch_info_local,&
       fetch_slot_col, fetch_slot_icols, fetch_slot_full_cols, &
       fetch_slot_data_begin, fetch_slot_data_end, fetch_slot_fetch_info, &
       fetch_slot_from, fetch_slot_index, fetch_slot_row,&
       fetch_slot_full_rows, fetch_slot_blocked_rows,&
       fetch_slot_tr, fetch_slot_offset
  USE cp_bcsr_util,                    ONLY: &
       addto_array, bcsr_gettime, cp_bcsr_addto_index_array, &
       cp_bcsr_checksum, cp_bcsr_create, cp_bcsr_create_dist2d_right, &
       cp_bcsr_destroy, cp_bcsr_destroy_fetch_env, cp_bcsr_fetch_pair, &
       cp_bcsr_fetch_pair_overlap, cp_bcsr_finalize, cp_bcsr_from_fm, &
       cp_bcsr_from_sm, cp_bcsr_printmat, cp_bcsr_release_fetch_slot,&
       cp_bcsr_release_fetch_slots,&
       cp_bcsr_set_debug, cp_bcsr_setup_fetch_env, cp_bcsr_verify_matrix, &
       cp_bcsr_work_create, ensure_array_size, cp_bcsr_data_window_complete,&
       cp_bcsr_switch_environments,&
       cp_bcsr_fetch_psubmatrix, cp_bcsr_fetch_sm_rows, checker_tr,&
       cp_bcsr_get_lookup, cp_bcsr_get_block_info, fetch_index_get_block_info,&
       aggregate_rows, swap, cp_bcsr_print, send_recv_fetch_data,&
       cp_bcsr_desymmetrize
  USE message_passing,                 ONLY: mp_window_fence, mp_sync, mp_sum,&
                                             mp_waitall
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type,&
                                             get_blacs_info,&
                                             cp_blacs_env_release
  USE cp_fm_basic_linalg,              ONLY: cp_fm_frobenius_norm
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_type, cp_fm_create, cp_fm_release
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_get,&
                                             distribution_2d_type,&
                                             distribution_2d_release
  USE kinds,                           ONLY: dp
  USE sparse_matrix_types,             ONLY: checksum_matrix,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE termination,                     ONLY: stop_memory
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_bcsr_operations'

  PUBLIC :: cp_bcsr_sm_sm_multiply,cp_bcsr_check

  PRIVATE
  !INCLUDE "mpif.h"
  REAL(KIND=dp) :: flop

CONTAINS

!> \par Information stored in the work matrix data structures.
!> For sm2_mult_general and _overlap, the %m%blk_p array (nominally) stores
!> offsets in the remote data pool. A negative value means the data block
!> has been transferred locally, its absolute value (i) indicating the offset
!> in the wm arrays. wm%row_i(i) stores the (negative) value of the data block
!> offset in the local buffer, wm%col_i(i) stores the original %m%blk_p value
!> (to revert when the local buffer is cleared), and wm%blk_p(i) stores the
!> block number of the original value (that is, where the original value was
!> stored in the %m%blk_p array).
!> \par For local matrices, the index is unchanged: the (all positive)
!> %m%blk_p point to the data block in the local data array.


! *****************************************************************************
!> \brief Performs a multiplication of two cp_bcsr_type matrices.
!> \par Matrices m_a and m_b are multiplied into the m_c product matrix. If the
!>      dist2d parameter is not specified, then a new distribution_2d is
!>      determined for it.
!> \param[in] matrix_a        left BCSR matrix
!> \param[in] matrix_b        right BCSR matrix
!> \param[out] matrix_c       resulting BCSR product matrix.
!> \param[in] dist2d          optional distribution_2d for the product matrix.
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE cp_bcsr_sm_sm_multiply(matrix_a, matrix_b, matrix_c, dist2d,&
       error)
    TYPE(cp_bcsr_type), POINTER              :: matrix_a, matrix_b, matrix_c
    TYPE(cp_bcsr_type), POINTER              :: matrix_afull
    TYPE(distribution_2d_type), POINTER      :: dist2d
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bcsr_sm_sm_multiply', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nblkcols_local, nblkrows_local
    INTEGER                                  :: timing_handle
    LOGICAL                                  :: dbg
    REAL(KIND=dp)                            :: checksum
    TYPE(cp_bcsr_work_type), POINTER         :: wm
!   ---------------------------------------------------------------------------
    dbg = .FALSE.
    CALL cp_bcsr_set_debug(dbg)
    CALL timeset(routineN, timing_handle)
    flop = 0.0_dp
    !> pre The matrices must be valid and they must have compatible shapes.
    ! Assign new dist2d
    ! Create & initialize the final & working matrix.
    IF (dbg) WRITE(*,*)routineP//'Creating 2-d distribution...'
    CALL mp_sync(matrix_a%dist2d%blacs_env%para_env%group)
    CALL distribution_2d_get(matrix_b%dist2d, error=error,&
         n_flat_local_rows=nblkrows_local, n_flat_local_cols=nblkcols_local)
    CALL cp_bcsr_create(matrix_c, "C", matrix_b%dist2d, .FALSE., .FALSE.,&
         (/0, 0,&
         matrix_b%nblkrows_total, matrix_b%nblkcols_total,&
         matrix_b%nfullrows_total, matrix_b%nfullrows_total,&
         matrix_b%nblkrows_local, matrix_b%nblkcols_local,&
         matrix_b%nfullrows_local, matrix_b%nfullcols_local/),&
         0, 0, error)
    IF (dbg) WRITE(*,*)routineP//'Updating new matrix...'
    CALL mp_sync(matrix_a%dist2d%blacs_env%para_env%group)
    CALL cp_bcsr_addto_index_array(matrix_c, cp_bcsr_slot_row_blk_begin,&
         matrix_a%row_blk_begin,&
         extra=matrix_a%nblkrows_total+2*matrix_b%nblkcols_total, error=error)
    CALL cp_bcsr_addto_index_array(matrix_c, cp_bcsr_slot_row_blk_end,&
         matrix_a%row_blk_end, error=error)
    CALL cp_bcsr_addto_index_array(matrix_c, cp_bcsr_slot_col_blk_begin,&
         matrix_b%col_blk_begin, error=error)
    CALL cp_bcsr_addto_index_array(matrix_c, cp_bcsr_slot_col_blk_end,&
         matrix_b%col_blk_end, error=error)
    IF (dbg) WRITE(*,*)routineP//' C initial index is', matrix_c%index
    IF (dbg) WRITE(*,*)routineP//'Creating work matrix...'
    CALL mp_sync(matrix_a%dist2d%blacs_env%para_env%group)
    CALL cp_bcsr_work_create(wm, matrix_c, matrix_c%nblkrows_local, sizedata_guess=SIZE(matrix_a%data), error=error)
    IF (dbg) WRITE(*,*)routineP//'Multiplying...'
    CALL mp_sync(matrix_a%dist2d%blacs_env%para_env%group)
    !CALL cp_bcsr_sm2_mult_overlap(matrix_a, matrix_b, wm, error)
    !CALL cp_bcsr_sm2_mult_general(matrix_a, matrix_b, wm, error)
    CALL cp_bcsr_sm2_mult_ls_rn(matrix_a, matrix_b, wm, error)
    IF (dbg) WRITE(*,*)routineP//' C prefinl index is', matrix_c%index
    CALL cp_bcsr_finalize(wm, .TRUE., error)
    IF (dbg) WRITE(*,*)routineP//' C final   index is', matrix_c%index
    IF (dbg) THEN
       CALL cp_bcsr_verify_matrix(matrix_c, error)
       CALL cp_bcsr_checksum(matrix_c, checksum, error)
       WRITE(*,*)routineP//' BCSR Sx Matrix checksum is', checksum
    ENDIF
    CALL mp_sum(flop, matrix_a%dist2d%blacs_env%para_env%group)
    WRITE(*,*)'FLOP=',flop
    CALL cp_bcsr_print(matrix_c, error)
    ! Perform multiplication
    CALL timestop(timing_handle)
  END SUBROUTINE cp_bcsr_sm_sm_multiply


! *****************************************************************************
!> \brief Perform actual multiplication of any two matrices.
!> \par The most general and slowest routine, does no packing of transfers.
!> \par Goes through the product matrix and for each element obtains the blocks
!>      from the two appropriate processers, multiplies them and places them
!>      into the matrix.
!> 
!> \param[in] m_a    left BCSR matrix
!> \param[in] m_b    right BCSR matrix
!> \param[out] m_c   resulting working BCSR product matrix.
!> \param error      cp2k error
! *****************************************************************************
  SUBROUTINE cp_bcsr_sm2_mult_general(m_a, m_b, m_c, error)
    TYPE(cp_bcsr_type), POINTER              :: m_a, m_b
    TYPE(cp_bcsr_work_type), POINTER         :: m_c
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bcsr_sm2_mult_general', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, fs_a, fs_b, k, mypcol_c, &
                                                myprow_c, row
    INTEGER, DIMENSION(:), POINTER           :: col_dist_c, row_dist_c
    LOGICAL                                  :: dbg, info
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_c
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(fetch_env), POINTER                 :: fetch_a1,&
                                                fetch_b1
!   ---------------------------------------------------------------------------
    dbg = .FALSE.
    info = .FALSE.
    CALL cp_bcsr_set_debug(dbg, info=info)
    CALL distribution_2d_get(m_c%m%dist2d, error=error, blacs_env=blacs_env_c,&
         row_distribution=row_dist_c, col_distribution=col_dist_c)
    CALL get_blacs_info(blacs_env_c, para_env=para_env,&
         my_process_row=myprow_c, my_process_column=mypcol_c)
    ! Initialize windows
    ! Initialize the fetch structures and communication environments.
    CALL cp_bcsr_setup_fetch_env(fetch_a1, para_env, m_a, error)
    CALL cp_bcsr_setup_fetch_env(fetch_b1, para_env, m_b, error)
    !> Go through my C-partition rows, within them columns, and for each
    !> (blocked) element, fetch the appropriate data and muliply it
    IF (dbg) WRITE(*,*)routineP//'End of setup'
    CALL mp_sync(para_env%group)
    IF (dbg) WRITE(*,*)routineP//'Starting transfers'
    DO row = 1, m_c%m%nblkrows_total
       IF(row_dist_c(row) .EQ. myprow_c) THEN
          DO col = 1, m_c%m%nblkcols_total
             IF(col_dist_c(col) .EQ. mypcol_c) THEN
                DO k = 1, m_a%nblkcols_total
                   CALL cp_bcsr_fetch_pair(fetch_a1, fetch_b1, row, k, col,&
                        fs_a, fs_b, error)
                   CALL do_simple_dgemm(m_c, fetch_a1, fetch_b1, fs_a, fs_b,&
                        error)
                   CALL cp_bcsr_release_fetch_slot(fetch_a1, fs_a,&
                        slotonly=.FALSE., error=error)
                   CALL cp_bcsr_release_fetch_slot(fetch_b1, fs_b,&
                        slotonly=.FALSE., error=error)
                ENDDO !k
             ENDIF
          ENDDO !col
       ENDIF
    ENDDO !row
    ! Destroy the fetch structure and communication environment.
    CALL cp_bcsr_destroy_fetch_env(fetch_a1, error)
    CALL cp_bcsr_destroy_fetch_env(fetch_b1, error)
  END SUBROUTINE cp_bcsr_sm2_mult_general


! *****************************************************************************
!> \brief Perform actual multiplication of any two matrices using
!>        overlapping computation and communication
!> \par The most general and slowest routine, does no packing of transfers.
!> \par Goes through the product matrix and for each element obtains the blocks
!>      from the two appropriate processers, multiplies them and places them
!>      into the matrix.
!> 
!> \param[in] m_a    left BCSR matrix
!> \param[in] m_b    right BCSR matrix
!> \param[out] m_c   resulting working BCSR product matrix.
!> \param error      cp2k error
! *****************************************************************************
  SUBROUTINE cp_bcsr_sm2_mult_overlap(m_a, m_b, m_c, error)
    TYPE(cp_bcsr_type), POINTER              :: m_a, m_b
    TYPE(cp_bcsr_work_type), POINTER         :: m_c
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bcsr_sm2_mult_overlap', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, k, &
                                                mypcol_c, myprow_c, &
                                                row, stat
    INTEGER, DIMENSION(:), POINTER           :: col_dist_c, row_dist_c
    LOGICAL                                  :: dbg = .FALSE., info = .FALSE.,&
                                                tocalc
    TYPE(cp_bcsr_index_lookup)               :: index_hint_a, index_hint_b
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_c
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(fetch_env), POINTER                 :: fetch_a1, fetch_a2,&
                                                fetch_a_comm, fetch_a_calc,&
                                                fetch_b1, fetch_b2,&
                                                fetch_b_comm, fetch_b_calc
!   ---------------------------------------------------------------------------
    CALL cp_bcsr_set_debug(dbg, info=info)
    CALL distribution_2d_get(m_c%m%dist2d, error=error, blacs_env=blacs_env_c,&
         row_distribution=row_dist_c, col_distribution=col_dist_c)
    CALL get_blacs_info(blacs_env_c, para_env=para_env,&
         my_process_row=myprow_c, my_process_column=mypcol_c)
    ! Initialize windows
    ! Initialize the fetch structure and communication environment.
    IF (dbg) WRITE(*,*)routineP//' Begin mult_overlap'
    CALL mp_sync(para_env%group)
    IF (dbg) WRITE(*,*)routineP//' Setup Environment A1'
    CALL mp_sync(para_env%group)
    CALL cp_bcsr_setup_fetch_env(fetch_a1, para_env, m_a, error)
    IF (dbg) WRITE(*,*)routineP//' Setup Environment A2'
    CALL mp_sync(para_env%group)
    CALL cp_bcsr_setup_fetch_env(fetch_a2, para_env, m_a, error, fetch_a1)
    IF (dbg) WRITE(*,*)routineP//' Setup Environment B1'
    CALL mp_sync(para_env%group)
    CALL cp_bcsr_setup_fetch_env(fetch_b1, para_env, m_b, error)
    IF (dbg) WRITE(*,*)routineP//' Setup Environment B2'
    CALL mp_sync(para_env%group)
    CALL cp_bcsr_setup_fetch_env(fetch_b2, para_env, m_b, error, fetch_b1)
    IF (dbg) WRITE(*,*)routineP//'End of environment setup'
    CALL mp_sync(para_env%group)
    IF (dbg) WRITE(*,*)routineP//'Starting transfers'
    NULLIFY(fetch_a_calc, fetch_b_calc)
    fetch_a_comm => fetch_a1
    fetch_b_comm => fetch_b1
    !
    index_hint_a%row = 0
    index_hint_b%row = 0
    ALLOCATE(fetch_a1%fs(m_a%nblkcols_total),&
         &   fetch_a2%fs(m_a%nblkcols_total), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'fetch_a1%fs or fetch_a2%fs',&
         m_a%nblkcols_total)
    fetch_a1%fs(:) = 0
    fetch_a2%fs(:) = 0
    ALLOCATE(fetch_b1%fs(m_a%nblkcols_total),&
         &   fetch_b2%fs(m_a%nblkcols_total), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'fetch_b1%fs or fetch_b2%fs',&
         m_a%nblkcols_total)
    fetch_b1%fs(:) = 0
    fetch_b2%fs(:) = 0
    tocalc = .FALSE. ! Nothing to calculate yet.
    ! Go through my C-partition rows, within them columns, and for each
    ! (blocked) element, fetch the appropriate data and muliply it
    DO row = 1, m_c%m%nblkrows_total
       IF(row_dist_c(row) .EQ. myprow_c) THEN
          DO col = 1, m_c%m%nblkcols_total
             IF(col_dist_c(col) .EQ. mypcol_c) THEN
                DO k = 1, m_a%nblkcols_total
                   CALL cp_bcsr_data_window_complete(fetch_a_comm,&
                        error, routineP)
                   CALL cp_bcsr_data_window_complete(fetch_b_comm,&
                        error, routineP)
                   IF (dbg)&
                        WRITE(*,*)routineP//' Inner loop COMM row, k, col',&
                        row, k, col
                   CALL cp_bcsr_fetch_pair_overlap(fetch_a_comm,fetch_b_comm,&
                        row, k, col,&
                        fetch_a_comm%fs(k), fetch_b_comm%fs(k), index_hint_a,&
                        index_hint_b, error)
                ENDDO
                IF (tocalc) THEN
                   DO k = 1, m_a%nblkcols_total
                      IF (dbg)&
                           WRITE(*,*)routineP//' Inner loop CALC row, k, col',&
                           row, k, col
                      CALL do_simple_dgemm(m_c, fetch_a_calc, fetch_b_calc,&
                           fetch_a_calc%fs(k), fetch_b_calc%fs(k), error)
                   ENDDO
                   !> Return the space used by the buffers.
                   CALL cp_bcsr_release_fetch_slots(fetch_a_calc, error)
                   CALL cp_bcsr_release_fetch_slots(fetch_b_calc, error)
                ENDIF
                !> Switch the two fetch windows
                CALL cp_bcsr_switch_environments(fetch_a_comm,&
                     fetch_a_calc, fetch_a1, fetch_a2, error)
                CALL cp_bcsr_switch_environments(fetch_b_comm,&
                     fetch_b_calc, fetch_b1, fetch_b2, error)
                tocalc = .TRUE.
             ENDIF
          ENDDO !col
       ENDIF
    ENDDO !row
    IF (dbg) WRITE(*,*)routineP//' Finalizing...'
    IF (tocalc) THEN
       !> Unlock windows, if necessary
       IF (dbg) WRITE(*,*)routineP//' Fencing...'
       CALL mp_window_fence(fetch_a_calc%window_data)
       CALL mp_window_fence(fetch_b_calc%window_data)
       IF (dbg) WRITE(*,*)routineP//' Fencing done.'
       CALL cp_bcsr_data_window_complete(fetch_a_calc, error, routineP)
       CALL cp_bcsr_data_window_complete(fetch_b_calc, error, routineP)
       IF (dbg) WRITE(*,*)routineP//' Fencing...'
       CALL mp_window_fence(fetch_a_calc%window_data)
       CALL mp_window_fence(fetch_b_calc%window_data)
       IF (dbg) WRITE(*,*)routineP//' Fencing done.'
       DO k = 1, m_a%nblkcols_total
          IF (dbg)&
               WRITE(*,*)routineP//' Inner loop CALC . k .', k
          CALL do_simple_dgemm(m_c, fetch_a_calc, fetch_b_calc,&
               fetch_a_calc%fs(k), fetch_b_calc%fs(k), error)
       ENDDO
    ENDIF
    ! Destroy the fetch structure and communication environment.
    DEALLOCATE(fetch_a1%fs, fetch_a2%fs, fetch_b1%fs, fetch_b2%fs)
    CALL cp_bcsr_destroy_fetch_env(fetch_a1, error)
    CALL cp_bcsr_destroy_fetch_env(fetch_a2, error, .TRUE.)
    CALL cp_bcsr_destroy_fetch_env(fetch_b1, error)
    CALL cp_bcsr_destroy_fetch_env(fetch_b2, error, .TRUE.)
  END SUBROUTINE cp_bcsr_sm2_mult_overlap


! *****************************************************************************
!> \brief Perform actual multiplication of a symmetric matrix and a non-
!>        symmetric, non-transposed one.
!> \par For the submatrix assigned to this processor, fetch the right matrix
!>      submatrix (having the same 2-D distribution as the product matrix),
!>      and the rows & columns of the left matrix needed for the product.
!> \param[in] m_a    left BCSR matrix
!> \param[in] m_b    right BCSR matrix
!> \param[out] m_c   resulting working BCSR product matrix.
!> \param error      cp2k error
! *****************************************************************************
  SUBROUTINE cp_bcsr_sm2_mult_ls_rn(m_a, m_b, m_c, error)
    TYPE(cp_bcsr_type), POINTER           :: m_a, m_b
    TYPE(cp_bcsr_work_type), POINTER      :: m_c
    TYPE(cp_error_type), INTENT(INOUT)    :: error
    !---------------------------------------
    TYPE(fetch_env), POINTER                 :: fetch_aN1, fetch_aT1,&
                                                fetch_aN2, fetch_aT2,&
                                                fetch_aN_comm, fetch_aT_comm,&
                                                fetch_aN_calc, fetch_aT_calc,&
                                                fetch_b1, fetch_b2,&
                                                fetch_b_comm, fetch_b_calc
    INTEGER           :: i
    INTEGER           :: p, row, col, k, nze, send_p, recv_p
    INTEGER           :: prow, pcol, pcol_i, nprows, npcols
    INTEGER           :: prow_a, pcol_a, nprows_a, npcols_a
    INTEGER           :: fs_b_comm, fs_b_calc, nfs_b_comm, nfs_b_calc,&
                         nfs_aN_comm, nfs_aT_comm, nfs_aN_calc, nfs_aT_calc
    INTEGER, DIMENSION(:), ALLOCATABLE, TARGET    :: fs_aN1, fs_aN2,&
                                                     fs_aT1, fs_aT2
    INTEGER, DIMENSION(:), POINTER        :: fs_aN_comm, fs_aT_comm,&
                                             fs_aN_calc, fs_aT_calc
    INTEGER, DIMENSION(:), POINTER        :: row_dist_c, col_dist_c,&
                                             row_dist_a, col_dist_a
    LOGICAL, DIMENSION(:), ALLOCATABLE    :: col_selection, row_selection
!    INTEGER, DIMENSION(:)                 :: col_list(m_a%nblkcols_total),&
!                                             row_list(m_a%nblkrows_total)
    INTEGER, DIMENSION(:)                 :: col_list(m_b%nblkrows_total),&
                                             row_list(m_b%nblkrows_total)
    INTEGER, DIMENSION(2)              :: send_req_aN_comm, send_req_aN_calc,&
                                          send_req_aT_comm, send_req_aT_calc,&
                                          recv_req_aN_comm, recv_req_aN_calc,&
                                          recv_req_aT_comm, recv_req_aT_calc
    INTEGER, DIMENSION(:), POINTER        :: local_rows_c, local_cols_c
    INTEGER                               :: n_local_rows_c, n_local_cols_c,&
                                             rowi, coli
    INTEGER           :: myprow_c, mypcol_c, myprow_a, mypcol_a
    INTEGER, DIMENSION(:,:), POINTER      :: blacs2mpi, blacs2mpi_a
    INTEGER           :: stat, timing_handle
    TYPE(cp_blacs_env_type), POINTER      :: blacs_env_c, blacs_env_a
    TYPE(cp_para_env_type), POINTER       :: para_env
    LOGICAL           :: failure = .TRUE., dosum, dbg, info, tocalc
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bcsr_sm2_mult_ls_rn',&
                                   routineP = moduleN//':'//routineN
!   ---------------------------------------------------------------------------
    dbg = .TRUE.
    info = .TRUE.
    CALL cp_bcsr_set_debug(dbg, info=info)
    CALL timeset(routineN, timing_handle)
    CALL distribution_2d_get(m_c%m%dist2d, error=error, blacs_env=blacs_env_c,&
         row_distribution=row_dist_c, col_distribution=col_dist_c,&
         n_flat_local_rows=n_local_rows_c, flat_local_rows=local_rows_c,&
         n_flat_local_cols=n_local_cols_c, flat_local_cols=local_cols_c)
    CALL distribution_2d_get(m_a%dist2d, error=error, blacs_env=blacs_env_a,&
         row_distribution=row_dist_a, col_distribution=col_dist_a)
    CALL get_blacs_info(blacs_env_c, para_env=para_env,&
         my_process_row=myprow_c, my_process_column=mypcol_c,&
         number_of_process_rows=nprows, number_of_process_columns=npcols,&
         blacs2mpi=blacs2mpi)
    CALL get_blacs_info(blacs_env_a,&
         my_process_row=myprow_a, my_process_column=mypcol_a,&
         number_of_process_rows=nprows_a,&
         number_of_process_columns=npcols_a, blacs2mpi=blacs2mpi_a)
    CALL ensure_array_size(m_c%row_i, ub=20*m_a%nblks, error=error)
    CALL ensure_array_size(m_c%col_i, ub=20*m_a%nblks, error=error)
    CALL ensure_array_size(m_c%blk_p, ub=20*m_a%nblks, error=error)
    ! Initialize the fetch structure and communication environment.
    CALL cp_bcsr_setup_fetch_env(fetch_aN1, para_env, m_a, error)
    CALL cp_bcsr_setup_fetch_env(fetch_aT1, para_env, m_a, error, fetch_aN1)
    CALL cp_bcsr_setup_fetch_env(fetch_aN2, para_env, m_a, error)
    CALL cp_bcsr_setup_fetch_env(fetch_aT2, para_env, m_a, error, fetch_aN2)
    CALL cp_bcsr_setup_fetch_env(fetch_b1, para_env, m_b, error)
    CALL cp_bcsr_setup_fetch_env(fetch_b2, para_env, m_b, error, fetch_b1)
    fetch_aN1%transposed = .FALSE. ; fetch_aN2%transposed = .FALSE.
    fetch_aT1%transposed = .TRUE.  ; fetch_aT2%transposed = .TRUE.
    fetch_aN_comm=>fetch_aN1
    fetch_aT_comm=>fetch_aT1
    fetch_b_comm=>fetch_b1
    NULLIFY(fetch_aN_calc) ; NULLIFY(fetch_aT_calc)
    NULLIFY(fetch_b_calc)
    !> For the B-partition, go through the my processor column and setup
    !> the fetches.
    ALLOCATE(row_selection(m_c%m%nblkrows_total), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'row_selection',&
         m_c%m%nblkrows_total)
    ALLOCATE(col_selection(m_c%m%nblkrows_total), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'col_selection',&
         m_c%m%nblkrows_total)
    ALLOCATE(fs_aN1(m_c%m%nblkrows_total*2), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'fs_aN1',&
         m_c%m%nblkrows_total*2)
    ALLOCATE(fs_aT1(m_c%m%nblkrows_total*2), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'fs_aT1',&
         m_c%m%nblkrows_total*2)
    ALLOCATE(fs_aN2(m_c%m%nblkrows_total*2), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'fs_aN2',&
         m_c%m%nblkrows_total*2)
    ALLOCATE(fs_aT2(m_c%m%nblkrows_total*2), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'fs_aT2',&
         m_c%m%nblkrows_total*2)
    fs_aN_comm => fs_aN1
    fs_aT_comm => fs_aT1
    NULLIFY (fs_aN_calc, fs_aT_calc)
    fs_aN_calc => fs_aN1
    fs_aT_calc => fs_aT1
    nfs_aN_comm = 0
    nfs_aN_calc = 0
    nfs_b_comm = 0
    nfs_b_calc = 0
    ! It must be true that the C & B distributions are the same!
    !col_selection(:) = col_dist_c(:) .EQ. mypcol_c
    ! The fetching cycle goes through B's (and my) process rows.
    ! Determine the (constant) non-transposed A matrix rows to fetch.
    !row_selection(:) = row_dist_a(:) .EQ. myprow_c
    IF (dbg) WRITE(*,*)routineP,' row dist', row_dist_c
    !IF (dbg) WRITE(*,*)routineP,' col dist', row_dist_c
    IF (dbg) WRITE(*,*)routineP,' c-nprows, mypcol_c', nprows, mypcol_c
    row_selection(:) = row_dist_c(:) .EQ. myprow_c
    !row_selection(:) = .TRUE.
    row_list(1:COUNT(row_selection)) =&
         PACK((/ (k, k = 1, SIZE(row_selection)) /), row_selection)
    tocalc = .FALSE.
    DO pcol_i = 0, nprows-1 ! Because npcols_a = nprows(_c)
       k = MOD(nprows+pcol_i+myprow_c+mypcol_c, nprows)
       p = blacs2mpi(k, mypcol_c)
       IF (dbg) WRITE(*,*)routineP//' pcol_i, k, p', pcol_i, k, p
       CALL mp_window_fence(fetch_b_comm%window_data,&
            no_pre_access=.TRUE., no_put_follows=.TRUE.)
       CALL mp_window_fence(fetch_aN_comm%window_data,&
            no_pre_access=.TRUE., no_put_follows=.TRUE.)
       CALL cp_bcsr_fetch_psubmatrix(fetch_b_comm, p, fs_b_comm, error)
       ! Fetch the non-transposed rows
       CALL aggregate_rows(fetch_aN_comm,&
            row_list(1:COUNT(row_selection)), error, k)
       col_selection(:) = row_dist_c(:) .EQ. k
       !col_selection(:) = .TRUE.
       IF (dbg) WRITE(*,*)routineP//' col selection', col_selection,'at k',k
       col_list(1:COUNT(col_selection)) =&
            PACK((/ (i, i = 1, SIZE(col_selection)) /), col_selection)
       CALL aggregate_rows(fetch_aT_comm,&
            col_list(1:COUNT(col_selection)), error, mypcol_a)
       IF (tocalc) THEN
          CALL mp_window_fence(fetch_aN_calc%window_data,&
               no_post_access=.TRUE., nothing_stored=.TRUE.)
          CALL mp_window_fence(fetch_b_calc%window_data,&
               no_post_access=.TRUE., nothing_stored=.TRUE.)
          CALL cp_bcsr_mult_aggregate_block(m_c, fetch_aN_calc, fetch_b_calc,&
               myprow_c, mypcol_c, k,&
               nfs_aN_calc, fs_aN_calc, 1, (/ fs_b_calc /), error)
          CALL cp_bcsr_mult_aggregate_block(m_c, fetch_aT_calc, fetch_b_calc,&
               myprow_c, mypcol_c, k,&
               nfs_aT_calc, fs_aT_calc, 1, (/ fs_b_calc /), error)
          CALL cp_bcsr_release_fetch_slot(fetch_b_calc, fs_b_calc,&
               reindex=.FALSE., slotonly=.FALSE., error=error)
       ENDIF
       CALL cp_bcsr_switch_environments(fetch_aN_comm, fetch_aN_calc,&
            fetch_aN1, fetch_aN2, error)
       CALL cp_bcsr_switch_environments(fetch_aT_comm, fetch_aT_calc,&
            fetch_aT1, fetch_aT2, error)
       CALL cp_bcsr_switch_environments(fetch_b_comm, fetch_b_calc,&
            fetch_b1, fetch_b2, error)
       IF (dbg) WRITE(*,*)routineP//' at pcol_i',pcol_i,'switched envs.'
       CALL swap(fs_b_comm, fs_b_calc)
       DO i = 1, 2
          CALL swap(send_req_aT_comm(i), send_req_aT_calc(i))
          CALL swap(send_req_aN_comm(i), send_req_aN_calc(i))
          CALL swap(recv_req_aT_comm(i), recv_req_aT_calc(i))          
          CALL swap(recv_req_aN_comm(i), recv_req_aN_calc(i))          
       ENDDO
       tocalc = .TRUE.
    ENDDO
    IF (dbg) WRITE(*,*)routineP//' Final calculation.'
    IF (tocalc) THEN
       CALL mp_window_fence(fetch_aN_calc%window_data,&
            no_post_access=.TRUE., nothing_stored=.TRUE.)
       CALL mp_window_fence(fetch_b_calc%window_data,&
            no_post_access=.TRUE., nothing_stored=.TRUE.)
       CALL cp_bcsr_mult_aggregate_block(m_c, fetch_aN_calc, fetch_b_calc,&
            myprow_c, mypcol_c, k,&
            nfs_aN_calc, fs_aN_calc, 1, (/ fs_b_calc /), error)
       CALL cp_bcsr_mult_aggregate_block(m_c, fetch_aT_calc, fetch_b_calc,&
            myprow_c, mypcol_c, k,&
            nfs_aT_calc, fs_aT_calc, 1, (/ fs_b_calc /), error)
       CALL cp_bcsr_release_fetch_slot(fetch_b_calc, fs_b_calc,&
            reindex=.FALSE., slotonly=.FALSE., error=error)
    ENDIF
    !CALL cp_bcsr_data_window_complete(fetch_a, error, routineP//'completion')
    !CALL cp_bcsr_data_window_complete(fetch_b, error, routineP//'completion')
    ! Destroy the fetch structure and communication environment.
    IF (dbg) WRITE(*,*)routineP//' Stopping.'
    CALL cp_bcsr_destroy_fetch_env(fetch_aN1, error)
    CALL cp_bcsr_destroy_fetch_env(fetch_aT1, error, .TRUE.)
    CALL cp_bcsr_destroy_fetch_env(fetch_aN2, error)
    CALL cp_bcsr_destroy_fetch_env(fetch_aT2, error, .TRUE.)
    CALL cp_bcsr_destroy_fetch_env(fetch_b1, error)
    CALL cp_bcsr_destroy_fetch_env(fetch_b2, error, .TRUE.)
    DEALLOCATE(row_selection, col_selection)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_bcsr_sm2_mult_ls_rn


! *****************************************************************************
!> \brief Multiplies an aggregated part of the left matrix and a whole part
!>        of the right matrix.
!> \param m_c       result matrix
!> \param[in] fetch_a         information on the left matrix row
!> \param[in] fetch_b         information on the right matrix column
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE cp_bcsr_mult_aggregate_block(m_c, fetch_a, fetch_b,&
       myprow_c, mypcol_c, pgrid_k, nfs_a, fs_a, nfs_b, fs_b, error)
    TYPE(cp_bcsr_work_type), POINTER    :: m_c
    TYPE(fetch_env), INTENT(INOUT)      :: fetch_a, fetch_b
    INTEGER, INTENT(IN)                 :: nfs_a, nfs_b, myprow_c, mypcol_c,&
                                           pgrid_k
    INTEGER, DIMENSION(:), INTENT(IN)   :: fs_a, fs_b
    TYPE(cp_error_type), INTENT(INOUT)  :: error
    !
    TYPE(cp_bcsr_type), POINTER    :: m_a, m_b
    INTEGER           :: i, k, iswap_tmp
    INTEGER           :: row_i, row_a, col_a, row_b, col_b, row_c, col_c
    INTEGER           :: col_i, virtual_row_a, virtual_col_a
    INTEGER           :: ifs_a, ifs_b, blk_a, blk_a_first, blk_a_last, blk_b, p_a, p_b, offset_a, offset_b
    INTEGER           :: fullrows, fullcols, lda, ldb, ldc
    INTEGER           :: datap, size_m, size_n, size_k, nze_a, nze_b, nze_c
    INTEGER, DIMENSION(:), POINTER      :: row_dist_c, col_dist_c
    INTEGER, DIMENSION(:, :), POINTER   :: fss, blacs2mpi_a, blacs2mpi_b
    TYPE(cp_blacs_env_type), POINTER    :: blacs_env_a, blacs_env_b
    REAL(KIND=dp)     :: beta
    REAL(KIND=dp), PARAMETER    :: alpha = 1.0
    REAL(KIND=dp), DIMENSION(:), POINTER  :: data_a, data_b
    CHARACTER         :: ctr_a, ctr_b
    INTEGER           :: timing_handle
    LOGICAL           :: failure = .TRUE., dosum, dbg, info, tr_a, tr_b,&
                         local_a, local_b, empty_a, empty_b, newblk
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bcsr_mult_aggregate_block',&
         &                         routineP = moduleN//':'//routineN
!   ---------------------------------------------------------------------------
    dbg = .FALSE.
    info = .TRUE.
    CALL cp_bcsr_set_debug(dbg, info=info)
    CALL timeset(routineN, timing_handle)
    IF (info) WRITE(*,*)routineP//' Multiplying for processor grid k',&
         pgrid_k,'transposed series',fetch_a%transposed
    CALL distribution_2d_get(m_c%m%dist2d, row_distribution=row_dist_c,&
         col_distribution=col_dist_c, error=error)
    CALL distribution_2d_get(fetch_a%m%dist2d, blacs_env=blacs_env_a,&
         error=error)
    CALL distribution_2d_get(fetch_b%m%dist2d, blacs_env=blacs_env_b,&
         error=error)
    CALL get_blacs_info(blacs_env_a, blacs2mpi=blacs2mpi_a)
    CALL get_blacs_info(blacs_env_b, blacs2mpi=blacs2mpi_b)
    p_b = fetch_b%fetch_slots(fetch_slot_from, fs_b(nfs_b))
    empty_b = fetch_b%fetch_slots(fetch_slot_fetch_info, fs_b(nfs_b)).EQ.&
         fetch_info_empty
    IF (dbg) WRITE(*,*)routineP//' B matrix'
    IF (dbg.AND..FALSE.) THEN
       WRITE(*,*)routineP//' B matrix row_p',fetch_b%m%row_p
       WRITE(*,*)routineP//' B matrix col_i',fetch_b%m%col_i
       WRITE(*,*)routineP//' B matrix blk_p',fetch_b%m%blk_p
       WRITE(*,*)routineP//' B matrix row_b',fetch_b%m%row_blk_begin
       WRITE(*,*)routineP//' B matrix row_e',fetch_b%m%row_blk_end
       WRITE(*,*)routineP//' B matrix col_b',fetch_b%m%col_blk_begin
       WRITE(*,*)routineP//' B matrix col_e',fetch_b%m%col_blk_end
       WRITE(*,*)routineP//' A matrix fetch index',fetch_a%index
    ENDIF
    !
    DO row_c = 1, m_c%m%nblkrows_total
       IF (row_dist_c(row_c) .NE. myprow_c) CYCLE
       DO col_c = 1, m_c%m%nblkcols_total
          IF (col_dist_c(col_c) .NE. mypcol_c) CYCLE
          !k_empty = .TRUE. ! To know when the first DGEMM is.
          newblk = .FALSE.
          dosum = .FALSE.
          DO k = 1, fetch_a%m%nblkrows_total
             IF (dbg) WRITE(*,*)routineP//'->Element',row_c,'x',col_c,'@',k,&
                  ' dosm:',dosum,'new block:',newblk
             ! First try to find the pointers to the data in the
             ! operand matrices.
             row_a = row_c
             col_a = k
             !CALL cp_bcsr_get_lookup(fetch_a%m, row_a, col_a, tr_a,&
             !     processor=p_a)
             tr_a = (col_a .GT. row_a) .EQV. BTEST(row_a+col_a, 0)
             IF (fetch_a%m%transposed) tr_a = .NOT. tr_a
             IF (tr_a) THEN
                iswap_tmp = row_a ; row_a = col_a ; col_a = iswap_tmp
             ENDIF
             p_a = fetch_a%m%dist2d%blacs_env%blacs2mpi(&
                  fetch_a%m%dist2d%row_distribution(row_a),&
                  fetch_a%m%dist2d%row_distribution(col_a))
             !
             IF (dbg) WRITE(*,*)routineP//' A lookup of',row_a,'x',col_a,&
                  ': tr_a, fetch env',tr_a,fetch_a%transposed
             local_a = p_a.EQ.fetch_a%para_env%mepos
             local_a = .FALSE.
             !tr_a = fetch_a%transposed
             !IF (tr_a.AND..NOT.fetch_a%transposed .OR.&
             !     .NOT.tr_a.AND.fetch_a%transposed) CYCLE
             IF (.NOT. tr_a.EQV.fetch_a%transposed) CYCLE
             CALL fetch_index_get_block_info(fetch_a, row_a, col_a, blk_a,&
                  offset_a)
             IF (dbg) THEN
                WRITE(*,*)routineP//' row col k:A row col p tr blk offset'
                WRITE(*,*)routineP, row_c, col_c, k,':',&
                     row_a, col_a, p_a, tr_a, blk_a, offset_a
             ENDIF
             IF (blk_a .EQ. 0) CYCLE
             row_b = k
             col_b = col_c
             !CALL cp_bcsr_get_lookup(fetch_b%m, row_b, col_b, tr_b,&
             !     processor=p_b)
             tr_b = (col_b .GT. row_b) .EQV. BTEST(row_b+col_b, 0)
             IF (fetch_b%m%transposed) tr_b = .NOT. tr_b
             IF (tr_b) THEN
                iswap_tmp = row_b ; row_b = col_b ; col_b = iswap_tmp
             ENDIF
             p_b = fetch_b%m%dist2d%blacs_env%blacs2mpi(&
                  fetch_b%m%dist2d%row_distribution(row_b),&
                  fetch_b%m%dist2d%row_distribution(col_b))
             !
             local_b = p_b.EQ.fetch_b%para_env%mepos
             CALL cp_bcsr_get_block_info(fetch_b%remote(p_b)%w%m,&
                  row_b, col_b, blk_b, data_offset=offset_b)
             IF (dbg) THEN
                WRITE(*,*)routineP//' row col k:B row col p tr blk offset'
                WRITE(*,*)routineP, row_c, col_c, k,':',&
                     row_b, col_b, p_b, tr_b, blk_b, offset_b
             ENDIF
             IF (blk_b .EQ. 0) CYCLE
             !IF (.NOT.local_b .AND. offset_b.GT.0) CYCLE
             offset_b = ABS(offset_b)
             ! Determine the sizes
             IF (dbg) WRITE(*,*)routineP//' Both are OK'
             lda = fetch_a%m%row_blk_end(row_a) -&
                  &fetch_a%m%row_blk_begin(row_a) + 1
             ldb = fetch_b%m%row_blk_end(row_b) -&
                  &fetch_b%m%row_blk_begin(row_b) + 1
             ldc = m_c%m%row_blk_end(row_c) -&
                  &m_c%m%row_blk_begin(row_c) + 1
             size_m = ldc
             size_k = fetch_b%m%row_blk_end(k) -&
                  &   fetch_b%m%row_blk_begin(k) + 1
             size_n = fetch_b%m%col_blk_end(col_c) -&
                  &   fetch_b%m%col_blk_begin(col_c) + 1
             IF (dbg) THEN
                !WRITE(*,*)'row',fetch_a%remote(p_a)%w%m%row_blk_begin
                !WRITE(*,*)'row',fetch_a%remote(p_a)%w%m%row_blk_end
                !WRITE(*,*)'col',fetch_a%remote(p_a)%w%m%col_blk_begin
                !WRITE(*,*)'col',fetch_a%remote(p_a)%w%m%col_blk_end
                WRITE(*,*)'A row',(fetch_a%m%row_blk_end-&
                     fetch_a%m%row_blk_begin)+1
                WRITE(*,*)'A col',(fetch_a%m%col_blk_end-&
                     fetch_a%m%col_blk_begin)+1
                WRITE(*,*)'B row',(fetch_b%m%row_blk_end-&
                     fetch_b%m%row_blk_begin)+1
                WRITE(*,*)'B col',(fetch_b%m%col_blk_end-&
                     fetch_b%m%col_blk_begin)+1
                WRITE(*,*)'C row',(m_c%m%row_blk_end-m_c%m%row_blk_begin)+1
                WRITE(*,*)'C col',(m_c%m%col_blk_end-m_c%m%col_blk_begin)+1
             ENDIF
             nze_a = lda * (fetch_a%m%col_blk_end(col_a) -&
                  &         fetch_a%m%col_blk_begin(col_a)+1)
             nze_b = ldb * (fetch_b%m%col_blk_end(col_b) -&
                  &         fetch_b%m%col_blk_begin(col_b)+1)
             nze_c = size_m * size_n
             ctr_a = 'N'
             IF (tr_a) ctr_a='T'
             ctr_b = 'N'
             IF (tr_b) ctr_b='T'
             datap = m_c%datasize+1
             IF (dbg) THEN
                WRITE(*,*)routineP//' lda, ldb, ldc=',lda,ldb,ldc
                WRITE(*,*)routineP//' sizes m, n, k', size_m, size_n, size_k
                WRITE(*,*)routineP//' nze a,b,c=', nze_a, nze_b, nze_c
                WRITE(*,*)routineP//' local; tr: A, B', local_a, local_b,';',&
                     ctr_a, ctr_b, datap
             ENDIF
             IF (local_a) THEN
                data_a => fetch_a%m%data(offset_a:offset_a+nze_a-1)
             ELSE
                data_a => fetch_a%buffer(offset_a:offset_a+nze_a-1)
             ENDIF
             IF (local_b) THEN
                data_b => fetch_b%m%data(offset_b:offset_b+nze_b-1)
             ELSE
                data_b => fetch_b%buffer(offset_b:offset_b+nze_b-1)
             ENDIF
             !> Prepare target.
             beta = 0.0
             ! Check history
             IF (m_c%lastblk .GT. 0) THEN
                DO i = m_c%lastblk, 1, -1
                   IF (m_c%row_i(i).EQ.row_c .AND. m_c%col_i(i).EQ.col_c) THEN
                      dosum = .TRUE.
                      datap = m_c%blk_p(i)
                      IF (dbg) WRITE(*,*)routineP//&
                           ' Found existing block at offset',datap,&
                           'pos',i,'of',m_c%lastblk
                      EXIT
                   ENDIF
                END DO
             ENDIF
             IF (dosum) THEN
                dosum = .TRUE.
                beta = 1.0
                flop = flop + size_m*size_n
             ELSE
                CALL ensure_array_size(m_c%data, ub=m_c%datasize+nze_c,&
                     nocopy=.FALSE., error=error)
                newblk = .TRUE.
             ENDIF
             flop = flop + size_m*size_n*size_k
             IF (dbg) THEN
                CALL cp_bcsr_printmat(data_a, size_m, size_k, tr=tr_a)
                CALL cp_bcsr_printmat(data_b, size_k, size_n, tr=tr_b)
                !IF (bcsr_verbose) WRITE(*,*)routineP//' B:',data_b
                WRITE(*,*)routineP//' A dp', DOT_PRODUCT(data_a, data_a)
                WRITE(*,*)routineP//' B dp', DOT_PRODUCT(data_b, data_b)
                IF (dosum) THEN
                   WRITE(*,*)routineP//' C dp', datap, &
                        DOT_PRODUCT(m_c%data(datap:datap+nze_c-1),&
                        &           m_c%data(datap:datap+nze_c-1))
                   !SUM(m_c%data(datap:datap+nze_c-1))
                ELSE
                   WRITE(*,*)routineP//' C will be cleared at', datap
                ENDIF
             ENDIF
             !
             CALL DGEMM(ctr_a, ctr_b, size_m, size_n, size_k,&
                  alpha, data_a, lda,&
                  data_b, ldb, beta,&
                  m_c%data(datap:datap+nze_c-1), ldc)
             IF (info) WRITE(*,*)routineP//' Done.',row_c,'x',col_c,'@',k,&
                  ' Summed:',dosum,'New block:',newblk, 'dp=',&
                  DOT_PRODUCT(m_c%data(datap:datap+nze_c-1),&
                  &           m_c%data(datap:datap+nze_c-1))
             IF (newblk) THEN
                m_c%datasize = m_c%datasize+nze_c
                m_c%lastblk = m_c%lastblk + 1
                CALL addto_array(m_c%row_i, (/row_c/), m_c%lastblk,&
                     error=error)
                CALL addto_array(m_c%col_i, (/col_c/), m_c%lastblk,&
                     error=error)
                CALL addto_array(m_c%blk_p, (/datap/), m_c%lastblk,&
                     error=error)
                !m_c%row_i(m_c%lastblk) = row_c
                !m_c%col_i(m_c%lastblk) = col_c
                !m_c%blk_p(m_c%lastblk) = datap
                IF (dbg) WRITE(*,*)routineP//' Added new block.', row_c,&
                     col_c, datap
                newblk = .FALSE.
             ELSE
                IF (dbg) WRITE(*,*)routineP//' No new block added.'
             ENDIF
          ENDDO ! k
          IF (dbg) WRITE(*,*)routineP//' Finished with',row_c,'x',col_c
       ENDDO ! row_c
    ENDDO ! col_c
    IF (dbg) WRITE(*,*)routineP//' Done.'
    CALL timestop(timing_handle)
  END SUBROUTINE cp_bcsr_mult_aggregate_block


! *****************************************************************************
!> \brief Performs simple DGEMM
! *****************************************************************************
  SUBROUTINE do_simple_dgemm(wm, fetch_a, fetch_b, fetchslot_a, fetchslot_b,&
       error)
    TYPE(cp_bcsr_work_type), POINTER         :: wm
    TYPE(fetch_env), POINTER                 :: fetch_a, fetch_b
    INTEGER, INTENT(IN)                      :: fetchslot_a, fetchslot_b
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'do_simple_dgemm', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: alpha = 1.0

    CHARACTER                                :: ctr_a, ctr_b
    INTEGER :: col, datap, stp, extent_i_a, extent_i_b, from_a, from_b, k, &
         lda, ldb, ldc, nze, row, size_k, size_m, size_n, start
    INTEGER                                  :: timing_handle
    LOGICAL                                  :: dbg, dosum, empty_a, &
         empty_b, failure = .TRUE., &
         info = .TRUE., local_a, local_b
    REAL(KIND=dp)                            :: beta
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_a, data_b
!   ---------------------------------------------------------------------------
    dbg = .FALSE.
    CALL cp_bcsr_set_debug(dbg, info=info)
    CALL timeset(routineN, timing_handle)
    IF (dbg) THEN
       WRITE(*,*)routineP//' slot A',fetchslot_a,'=',&
            fetch_a%fetch_slots(:,fetchslot_a)
       WRITE(*,*)routineP//' slot B',fetchslot_b,'=',&
            fetch_b%fetch_slots(:,fetchslot_b)
    ENDIF
    empty_a = fetch_a%fetch_slots(fetch_slot_fetch_info, fetchslot_a) .EQ.&
         fetch_info_empty
    empty_b = fetch_b%fetch_slots(fetch_slot_fetch_info, fetchslot_b) .EQ.&
         fetch_info_empty
    local_a = (fetch_a%fetch_slots(fetch_slot_fetch_info, fetchslot_a) .EQ.&
         &     fetch_info_local)
    local_b = (fetch_b%fetch_slots(fetch_slot_fetch_info, fetchslot_b) .EQ.&
         &     fetch_info_local)
    IF (dbg) WRITE(*,*)routineP//' Localness:', local_a, local_b
    IF (empty_a .OR. empty_b) THEN
       IF (dbg) WRITE(*,*)routineP//' Emptiness:', empty_a, empty_b
    ELSE
       !> Setup indices
       row = fetch_a%fetch_slots(fetch_slot_row, fetchslot_a)
       k = fetch_a%fetch_slots(fetch_slot_col, fetchslot_a)
       col = fetch_b%fetch_slots(fetch_slot_col, fetchslot_b)
       IF (dbg) THEN
          WRITE(*,*)routineP//' positions A', fetch_a%fetch_slots(1:2, fetchslot_a)
          WRITE(*,*)routineP//' positions B', fetch_b%fetch_slots(1:2, fetchslot_b)
          WRITE(*,*)routineP//' row, k, col=', row, k, col
       ENDIF
       CALL cp_assert(fetch_b%fetch_slots(fetch_slot_row, fetchslot_b) .EQ. k,&
            cp_warning_level,&
            cp_internal_error, routineP, 'Mismatch in common row and column',&
            error=error, failure=failure)
       lda = fetch_a%fetch_slots(fetch_slot_full_rows, fetchslot_a)
       !lda = fetch_a%m%row_blk_end(row) - fetch_a%m%row_blk_begin(row)+1
       IF (fetch_a%fetch_slots(fetch_slot_tr, fetchslot_a) .EQ. 0) THEN
          size_m = fetch_a%fetch_slots(fetch_slot_full_rows, fetchslot_a)
          size_k = fetch_a%fetch_slots(fetch_slot_full_cols, fetchslot_a)
          ctr_a = 'N'
       ELSE
          size_m = fetch_a%fetch_slots(fetch_slot_full_cols, fetchslot_a)
          size_k = fetch_a%fetch_slots(fetch_slot_full_rows, fetchslot_a)
          ctr_a = 'T'
       ENDIF
       ldb = fetch_b%fetch_slots(fetch_slot_full_rows, fetchslot_b)
       !ldb = fetch_b%m%row_blk_end(row) - fetch_b%m%row_blk_begin(row)+1
       IF (fetch_b%fetch_slots(fetch_slot_tr, fetchslot_b) .EQ. 0) THEN
          size_n = fetch_b%fetch_slots(fetch_slot_full_cols, fetchslot_b)
          ctr_b = 'N'
          IF(dbg) CALL cp_assert(size_k .EQ.&
               fetch_b%fetch_slots(fetch_slot_full_rows, fetchslot_b),&
               cp_warning_level, cp_internal_error, routineP,&
               'Row and column sizes do not match', error=error, failure=failure)
       ELSE
          size_n = fetch_b%fetch_slots(fetch_slot_full_rows, fetchslot_b)
          ctr_b = 'T'
          IF(dbg) CALL cp_assert(size_k .EQ.&
               fetch_b%fetch_slots(fetch_slot_full_cols, fetchslot_b),&
               cp_warning_level, cp_internal_error, routineP,&
               'Row and column sizes do not match', error=error, failure=failure)
       ENDIF
       nze = size_m * size_n
       !> Setup data pointers.
       from_a = fetch_a%fetch_slots(fetch_slot_from, fetchslot_a)
       start = fetch_a%fetch_slots(fetch_slot_data_begin, fetchslot_a)
       stp= fetch_a%fetch_slots(fetch_slot_data_end, fetchslot_a)
       extent_i_a = fetch_a%fetch_slots(fetch_slot_index, fetchslot_a)
       IF (dbg) WRITE(*,*)routineP//' A span',start,':',stp,'using extent',&
            extent_i_a!,&
       !'=',fetch%remote_a(from_a)%w%extents(1, i),'size',&
       !fetch%remote_a(from_a)%w%extents(2, i)
       !CALL cp_assert(start .EQ. fetch%remote_a(from_a)%w%extents(1, i),&
       !     cp_warning_level, cp_internal_error, routineP,&
       !     'Extent starts not equal', error=error, failure=failure)
       !CALL cp_assert(end-start+1 .EQ.&
       !     ABS(fetch%remote_a(from_a)%w%extents(2, i)),&
       !     cp_warning_level, cp_internal_error, routineP,&
       !     'Extent sizes not equal', error=error, failure=failure)
       IF (local_a) THEN
          data_a => fetch_a%remote(from_a)%w%data(start:stp)
       ELSE
          data_a => fetch_a%buffer(start:stp)
       ENDIF
       
       !
       from_b = fetch_b%fetch_slots(fetch_slot_from, fetchslot_b)
       start = fetch_b%fetch_slots(fetch_slot_data_begin, fetchslot_b)
       stp = fetch_b%fetch_slots(fetch_slot_data_end, fetchslot_b)
       extent_i_b = fetch_b%fetch_slots(fetch_slot_index, fetchslot_b)
       IF (dbg) WRITE(*,*)routineP//' B span',start,':',stp,'using extent',&
            extent_i_b!,&
       !'=',fetch%remote_b(from_b)%w%extents(1, i),'size',&
       !fetch%remote_b(from_a)%w%extents(2, i)
       !CALL cp_assert(start .EQ.&
       !     fetch%remote_b(from_b)%w%extents(1, i),&
       !     cp_warning_level, cp_internal_error, routineP,&
       !     'Extent starts not equal', error=error, failure=failure)
       !CALL cp_assert(end-start+1 .EQ.&
       !     ABS(fetch%remote_b(from_b)%w%extents(2, i)),&
       !     cp_warning_level, cp_internal_error, routineP,&
       !     'Extent sizes not equal', error=error, failure=failure)
       IF (local_b) THEN
          data_b => fetch_b%remote(from_b)%w%data(start:stp)
       ELSE
          data_b => fetch_b%buffer(start:stp)
       ENDIF
       !
       IF (dbg) WRITE(*,*)routineP//' Multipling row k col',row,k,col,&
            'sizes',size_m,size_k,size_n,'nze=',nze
       !IF (dbg) WRITE(*,*)routineP//'        '//' A (',fetch%nr_a,'x',&
       !     fetch%nc_a,fetch%tr_a,') B (',fetch%nr_b,'x',fetch%nc_b,fetch%tr_b,&
       !     ')'
       !
       dosum = .FALSE.
       datap = wm%datasize+1
       beta = 0.0
       IF(wm%lastblk .GT. 0) THEN
          IF(wm%row_i(wm%lastblk).EQ.row .AND. wm%col_i(wm%lastblk).EQ.col) THEN
             dosum = .TRUE.
             beta = 1.0
             datap = wm%blk_p(wm%lastblk)
          ENDIF
       ENDIF
       IF (.NOT.dosum) THEN
          wm%datasize = wm%datasize+nze
          wm%lastblk = wm%lastblk + 1
          CALL addto_array(wm%row_i, (/row/), wm%lastblk, error=error)
          CALL addto_array(wm%col_i, (/col/), wm%lastblk, error=error)
          CALL addto_array(wm%blk_p, (/datap/), wm%lastblk, error=error)
          CALL ensure_array_size(wm%data, ub=datap+nze-1, nocopy=.FALSE.,&
               error=error)
          !CALL ensure_array_size_i(wm%data, 1, datap+nze-1, 1.4, .FALSE., error)
       ENDIF
       ldc = size_m
       IF (dbg) WRITE(*,*)routineP//'        '//' In-place:',dosum,&
            'blks',wm%lastblk,'at pos',datap,'to',datap+nze-1,&
            'of size',wm%datasize,'&',SIZE(wm%data)
       IF (dbg) WRITE(*,*)routineP//'        '//' DGEMM parameters',&
            ctr_a, ctr_b, size_m, size_n, size_k,&
            alpha, 'fetch%buffer_a', lda,&
            'fetch%buffer_b', ldb, beta,&
            'wm%data(datap:datap+nze-1)', ldc
       !assert ldc = row-blk-size of C
       IF (dbg) THEN
          WRITE(*,*)routineP//'        '//' LDC=',ldc,'row',row,'size=',&
               wm%m%row_blk_end(row)-wm%m%row_blk_begin(row)+1
          CALL cp_assert(ldc.EQ.wm%m%row_blk_end(row)-wm%m%row_blk_begin(row)+1,&
               cp_warning_level, cp_internal_error, routineP,&
               'Target row size invalid', error=error, failure=failure)
       ENDIF
       IF (dbg) THEN
          CALL cp_bcsr_printmat(data_a,&
               fetch_a%fetch_slots(fetch_slot_full_rows, fetchslot_a),&
               fetch_a%fetch_slots(fetch_slot_full_cols, fetchslot_a),&
               'matrix A',&
               fetch_a%fetch_slots(fetch_slot_tr, fetchslot_a).NE.0)
          CALL cp_bcsr_printmat(data_b,&
               fetch_b%fetch_slots(fetch_slot_full_rows, fetchslot_b),&
               fetch_b%fetch_slots(fetch_slot_full_cols, fetchslot_b),&
               'matrix B',&
               fetch_a%fetch_slots(fetch_slot_tr, fetchslot_a).NE.0)
          IF (beta .NE. 0) CALL cp_bcsr_printmat(wm%data(datap:datap+nze-1),&
               size_m, size_n, 'matrix C-pre', .FALSE.)
          WRITE(*,*)routineP//' checksum A', DOT_PRODUCT(data_a, data_a)
          WRITE(*,*)routineP//' checksum B', DOT_PRODUCT(data_b, data_b)
          IF (beta .NE. 0) &
               WRITE(*,*)routineP//' checksum C pre', beta,&
               DOT_PRODUCT(wm%data(datap:datap+nze-1),&
               &           wm%data(datap:datap+nze-1))
       ENDIF
       CALL DGEMM(ctr_a, ctr_b, size_m, size_n, size_k,&
            alpha, data_a, lda,&
            data_b, ldb, beta,&
            wm%data(datap:datap+nze-1), ldc)
       k = nze
       nze = size_m * size_n * size_k
       IF (dosum) THEN
          flop = flop + (nze + size_m*size_k)
       ELSE
          flop = flop + nze
       ENDIF
       IF (dbg) THEN
          nze = k
          WRITE(*,*)routineP//' checksum C post', beta, &
               DOT_PRODUCT(wm%data(datap:datap+nze-1),&
               &           wm%data(datap:datap+nze-1))
       ENDIF
       !CALL cp_bcsr_printmat(wm%data(datap:datap+nze-1),&
       !     size_m, size_n, 'matrix C-post', .FALSE.)
    ENDIF ! empty
!    IF (dbg.AND..NOT.empty_a.AND..NOT.empty_b) THEN
!       !WRITE(*,*)routineP//' ',wm%data(datap:datap+nze-1)
!       WRITE(*,*)routineP//'        '//&
!            ' block checksum pos(',datap,'-',datap+nze-1,')=',&
!            DOT_PRODUCT(wm%data(datap:datap+nze-1), wm%data(datap:datap+nze-1))
!       !CALL cp_bcsr_printmat(wm%data(datap:datap+nze-1), size_m, size_n)
!    ENDIF
    CALL timestop(timing_handle)
  END SUBROUTINE do_simple_dgemm



! *****************************************************************************
!> \brief   Checks the new bcsr routines, to find bugs early
!> \author  Urban
!> \date    2008-10-07
! *****************************************************************************
  SUBROUTINE cp_bcsr_check(matrix_s, matrix_x, error)

    TYPE(cp_fm_type), POINTER                :: matrix_x, fmout
    TYPE(real_matrix_type), POINTER          :: matrix_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bcsr_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, k
    INTEGER, DIMENSION(:), POINTER           :: col_blk_begin, col_blk_end, &
                                                row_blk_begin, row_blk_end
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: checksum, oldchecksum, newtime,&
                                                walltime1, &
                                                walltime2, oldtime
    TYPE(cp_bcsr_type), POINTER              :: matrix_a, matrix_b, matrix_c
    TYPE(distribution_2d_type), POINTER      :: dist2d, dist2d_left
    TYPE(cp_blacs_env_type), POINTER         :: tmp_blacs_env
    LOGICAL                                  :: dbg = .FALSE.

    CALL timeset(routineN,handle)

    failure=.FALSE.

    !Do the new multiplication
    dist2d_left => matrix_s%distribution_2d
    IF (dbg) WRITE(*,*)'left-1',dist2d_left%ref_count,dist2d_left%blacs_env%ref_count,&
         dist2d_left%blacs_env%para_env%ref_count
    CALL cp_bcsr_from_sm(matrix_a, matrix_s, error)
    CALL cp_bcsr_verify_matrix(matrix_a, error)
    WRITE(*,*)'S checksum=',&
         checksum_matrix(matrix_s, matrix_s%distribution_2d%blacs_env%para_env)
    checksum = 0.0
    CALL cp_bcsr_checksum(matrix_a, checksum, error)
    WRITE(*,*)'BCSR S checksum',checksum
    !row_blk_begin => matrix_a%row_blk_begin
    !row_blk_end => matrix_a%row_blk_end
    CALL cp_fm_frobenius_norm(matrix_x, checksum, error)
    CALL cp_fm_get_info(matrix_x, ncol_global=k, error=error)
    CALL cp_bcsr_create_dist2d_right(dist2d, dist2d_left, k,&
         col_blk_begin, col_blk_end, error=error)
    IF (dbg) WRITE(*,*)'right0=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
         dist2d%blacs_env%para_env%ref_count
    CALL cp_bcsr_from_fm(matrix_b, matrix_x, 0.0d0, dist2d,&
         row_blk_begin=matrix_a%row_blk_begin,&
         row_blk_end=matrix_a%row_blk_end,&
         col_blk_begin=col_blk_begin, col_blk_end=col_blk_end,&
         error=error)
    DEALLOCATE(col_blk_begin, col_blk_end)
    IF (dbg) WRITE(*,*)'right1=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
         dist2d%blacs_env%para_env%ref_count
    CALL cp_bcsr_verify_matrix(matrix_b, error)
    WRITE(*,*)'x checksum',checksum*checksum
    CALL cp_bcsr_checksum(matrix_b, checksum, error)
    WRITE(*,*)'BCSR x checksum',checksum
    !
    walltime1 = bcsr_gettime()
    CALL cp_bcsr_desymmetrize(matrix_a, matrix_afull, error)
!!    CALL cp_bcsr_sm_sm_multiply(matrix_a, matrix_b, matrix_c, dist2d, error)
    IF (dbg) WRITE(*,*)'right2=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
         dist2d%blacs_env%para_env%ref_count
    walltime2 = bcsr_gettime()
    newtime = walltime2-walltime1
!    !
!    CALL cp_bcsr_verify_matrix(matrix_c, error)
!    !
!    CALL cp_sm_fm_multiply(matrix_s,matrix_x, &
!         matrix_sx, k,error=error)
!    walltime2 = bcsr_gettime()
!    oldtime = walltime2-walltime1
!    !
!    CALL cp_fm_frobenius_norm(matrix_sx, checksum, error)
!    write(*,*)'Sx checksum',checksum*checksum,'time=',oldtime
!!    CALL cp_bcsr_checksum(matrix_c, checksum, error)
!!    WRITE(*,*)'BCSR Sx checksum',checksum,'time=',newtime
    CALL cp_bcsr_destroy(matrix_a,error)
    CALL cp_bcsr_destroy(matrix_b,error)
!!    CALL cp_bcsr_destroy(matrix_c,error)
    IF (dbg) WRITE(*,*)'right3=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
         dist2d%blacs_env%para_env%ref_count
    tmp_blacs_env => dist2d%blacs_env
    CALL distribution_2d_release(dist2d, error)
    IF (dbg) WRITE(*,*)'right4=','"',0,'"',tmp_blacs_env%ref_count,&
         tmp_blacs_env%para_env%ref_count
    CALL cp_blacs_env_release(tmp_blacs_env, error)
    IF (dbg) WRITE(*,*)'left99',dist2d_left%ref_count,dist2d_left%blacs_env%ref_count,&
         dist2d_left%blacs_env%para_env%ref_count
    !WRITE(*,*)'right4=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
    !     dist2d%blacs_env%para_env%ref_count

    CALL timestop(handle)

    CALL timeset(routineN//'std', handle)

    CALL cp_fm_get_info(matrix_x, ncol_global=k, error=error)
    CALL cp_fm_create(fmout, matrix_x%matrix_struct, "FMout", error)
    walltime1 = bcsr_gettime()
    CALL cp_sm_fm_multiply(matrix_s, matrix_x, fmout, k, error=error)
    walltime2 = bcsr_gettime()
    oldtime = walltime2-walltime1
    CALL cp_fm_frobenius_norm(fmout, checksum, error)
    write(*,*)'Sx checksum',checksum*checksum,'time=',oldtime
    CALL cp_fm_release(fmout, error)

    CALL timestop(handle)


  END SUBROUTINE cp_bcsr_check


END MODULE cp_bcsr_operations
